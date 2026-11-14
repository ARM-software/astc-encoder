// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2026 Arm Limited
//
// Licensed under the Apache License, Version 2.0 (the "License"); you may not
// use this file except in compliance with the License. You may obtain a copy
// of the License at:
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
// WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
// License for the specific language governing permissions and limitations
// under the License.
// ----------------------------------------------------------------------------

/**
 * @file
 * @brief Fused RDO decoding and distortion for LDR blocks.
 */

#if !defined(ASTCENC_DECOMPRESS_ONLY)

#include <cstring>

#include "astcenc_internal.h"
#include "astcenc_rdo.h"

/** @brief Unpacked trial endpoints shared by RGB and RGBA distortion kernels. */
struct rdo_trial_endpoints
{
	/** @brief Low endpoints for the active partitions. */
	vint4 endpoint0[4];

	/** @brief Differences between each partition's high and low endpoints. */
	vint4 deltas[4];
};

/* See header for documentation. */
float rdo_get_mse_limit(
	float best_cost,
	float rate_cost,
	float mse_recip
) {
	if (rate_cost >= best_cost)
	{
		return 0.0f;
	}

	return (best_cost - rate_cost) * mse_recip;
}

/**
 * @brief Decode UNORM16 colors through the decoder's truncating FP16 intermediate.
 *
 * UNORM16 maps to p / 65536 with ten fraction bits, except 65535 maps to one.
 * Clearing thirteen FP32 fraction bits gives the same result, including FP16 subnormals.
 *
 * @param color   The unsigned 16-bit color values.
 *
 * @return The decoded LDR values in [0, 1].
 */
static ASTCENC_SIMD_INLINE vfloat rdo_decode_ldr_color(
	vint color
) {
	vfloat normalized = int_to_float(color) * (1.0f / 65536.0f);
	vfloat truncated = int_as_float(float_as_int(normalized) & vint(~0x1FFF));
	return select(truncated, vfloat(1.0f), color == vint(65535));
}

/**
 * @brief Decode four UNORM16 colors through the decoder's truncating FP16 intermediate.
 *
 * @param color   The unsigned 16-bit color values.
 *
 * @return The decoded LDR values in [0, 1].
 */
static ASTCENC_SIMD_INLINE vfloat4 rdo_decode_ldr_color_4(
	vint4 color
) {
	vfloat4 normalized = int_to_float(color) * (1.0f / 65536.0f);
	vfloat4 truncated = int_as_float(float_as_int(normalized) & vint4(~0x1FFF));
	return select(truncated, vfloat4(1.0f), color == vint4(65535));
}

/**
 * @brief Interpolate only the batch of weights needed by the next distortion calculation.
 *
 * @param[in] di                 The fixed decimation information.
 * @param[in] weights            The unquantized stored weights.
 * @param[in] table              The prepared lookup table, used only for decimated blocks.
 * @param     texel_offset       The first texel in the SIMD batch.
 * @param     plane_offset       The stored weight-plane offset, either zero or WEIGHTS_PLANE2_OFFSET.
 * @param     max_weight_count   The maximum contributing weight count for this batch.
 *
 * @return The integer interpolation weights for the batch.
 */
static ASTCENC_SIMD_INLINE vint rdo_interpolate_weights(
	const decimation_info& di,
	const uint8_t* weights,
	const vtable_64x8& table,
	unsigned int texel_offset,
	unsigned int plane_offset,
	unsigned int max_weight_count
) {
	if (di.weight_count == di.texel_count)
	{
		return vint(weights + plane_offset + texel_offset);
	}

	vint summed_value(8);
	promise(max_weight_count > 0);
	for (unsigned int j = 0; j < max_weight_count; j++)
	{
		vint texel_weights = vint(di.texel_weights_tr[j] + texel_offset) + vint(plane_offset);
		vint texel_weights_int(di.texel_weight_contribs_int_tr[j] + texel_offset);
		summed_value += vtable_lookup_32bit(table, texel_weights) * texel_weights_int;
	}

	return lsr<4>(summed_value);
}

/**
 * @brief Reuse original cached weights or interpolate a changed trial payload.
 *
 * @param[in] di             The fixed decimation information.
 * @param[in] weights        The trial's unquantized stored weights.
 * @param[in] table          The prepared lookup table, used only for decimated blocks.
 * @param     texel_offset   The first texel in the SIMD batch.
 * @param     plane_offset   The stored weight-plane offset, either zero or WEIGHTS_PLANE2_OFFSET.
 * @param[in] cache          The initialized original reconstruction cache.
 *
 * @return The integer interpolation weights for the batch.
 */
static ASTCENC_SIMD_INLINE vint rdo_interpolate_weights(
	const decimation_info& di,
	const uint8_t* weights,
	const vtable_64x8& table,
	unsigned int texel_offset,
	unsigned int plane_offset,
	const rdo_block_difference_cache& cache
) {
	if (cache.is_weights_unchanged)
	{
		unsigned int plane = 0;
		if (plane_offset != 0)
		{
			plane = 1;
		}

		return vint(cache.weights[plane] + texel_offset);
	}

	return rdo_interpolate_weights(di, weights, table, texel_offset, plane_offset,
	                               cache.weight_counts[texel_offset / ASTCENC_SIMD_WIDTH]);
}

/* See header for documentation. */
void rdo_prepare_difference_cache(
	const astcenc_config& config,
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& block,
	rdo_block_difference_cache& cache
) {
	assert(block.block_type == SYM_BTYPE_NONCONST);
	static_assert(BLOCK_MAX_TEXELS <= 224);
	static_assert(ASTCENC_SIMD_WIDTH >= 4);

	// Unchanged payloads can reuse reconstruction, independent of the bytes copied by a trial
	for (unsigned int p = 0; p < block.partition_count; p++)
	{
		vint4 endpoint0;
		vint4 endpoint1;
		bool is_rgb_lns;
		bool is_alpha_lns;
		[[maybe_unused]] bool is_valid = unpack_color_endpoints(
		    config.profile, block.color_formats[p], block.color_values[p],
		    is_rgb_lns, is_alpha_lns, endpoint0, endpoint1);
		assert(is_valid);

		store(endpoint0, cache.endpoints[p][0]);
		store(endpoint1, cache.endpoints[p][1]);
	}

	const auto& mode = bsd.get_block_mode(block.block_mode);
	const auto& di = bsd.get_decimation_info(mode.decimation_mode);
	unsigned int plane_count = 1;
	if (mode.is_dual_plane)
	{
		plane_count = 2;
	}

	// Contributor counts depend only on the fixed layout, so reduce them once for all trials
	for (unsigned int i = 0; i < bsd.texel_count; i += ASTCENC_SIMD_WIDTH)
	{
		cache.weight_counts[i / ASTCENC_SIMD_WIDTH] =
		    static_cast<uint8_t>(hmax_s(vint(di.texel_weight_count + i)));
	}

	vtable_64x8 table;
	if (di.weight_count != bsd.texel_count)
	{
		vtable_prepare(table, block.weights);
	}

	for (unsigned int plane = 0; plane < plane_count; plane++)
	{
		unsigned int plane_offset = plane * WEIGHTS_PLANE2_OFFSET;
		if (di.weight_count == bsd.texel_count)
		{
			std::memcpy(cache.weights[plane], block.weights + plane_offset, bsd.texel_count);
			size_t padded_count = round_up_to_simd_multiple_vla(bsd.texel_count);
			std::memset(cache.weights[plane] + bsd.texel_count, 0, padded_count - bsd.texel_count);
			continue;
		}

		for (unsigned int i = 0; i < bsd.texel_count; i += ASTCENC_SIMD_WIDTH)
		{
			vint weights = rdo_interpolate_weights(di, block.weights, table, i, plane_offset,
			                                       cache.weight_counts[i / ASTCENC_SIMD_WIDTH]);
			pack_and_store_low_bytes(weights, cache.weights[plane] + i);
		}
	}

	size_t padded_count = round_up_to_simd_multiple_vla(bsd.texel_count);
	std::fill_n(cache.texel_errors + bsd.texel_count, padded_count - bsd.texel_count, 0.0f);

	cache.is_colors_unchanged = true;
	cache.is_weights_unchanged = true;
}

/* See header for documentation. */
void rdo_prepare_difference_source(
	image_block& original,
	unsigned int texel_count
) {
	for (unsigned int i = 0; i < texel_count; i++)
	{
		vfloat4 texel = original.texel(i) * (255.0f / 65535.0f);
		original.data_rgba[i] = texel;
		original.data_r[i] = texel.lane<0>();
		original.data_g[i] = texel.lane<1>();
		original.data_b[i] = texel.lane<2>();
		original.data_a[i] = texel.lane<3>();
	}

	// Initialize inactive source lanes once per block for a padded final SIMD batch
	size_t padded_count = round_up_to_simd_multiple_vla(texel_count);
	for (size_t i = texel_count; i < padded_count; i++)
	{
		original.data_r[i] = 0.0f;
		original.data_g[i] = 0.0f;
		original.data_b[i] = 0.0f;
		original.data_a[i] = 0.0f;
	}

	original.data_mean = original.data_mean * (255.0f / 65535.0f);
	original.data_min = original.data_min * (255.0f / 65535.0f);
	original.data_max = original.data_max * (255.0f / 65535.0f);
}

/* See header for documentation. */
rdo_block_statistics rdo_compute_original_difference(
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& block,
	const image_block& original,
	rdo_block_difference_cache& cache
) {
	const auto& mode = bsd.get_block_mode(block.block_mode);
	const auto& partition = bsd.get_partition_info(block.partition_count, block.partition_index);
	vmask4 plane2_mask = vint4::lane_id() == vint4(block.plane2_component);
	vint4 endpoints[4];
	vint4 deltas[4];

	for (unsigned int p = 0; p < block.partition_count; p++)
	{
		endpoints[p] = vint4(cache.endpoints[p][0]);
		deltas[p] = vint4(cache.endpoints[p][1]) - endpoints[p];
	}

	// Reconstruct directly into the error cache without an intermediate decoded image block
	for (unsigned int i = 0; i < bsd.texel_count; i++)
	{
		unsigned int p = partition.partition_of_texel[i];
		vint4 weight(cache.weights[0][i]);
		if (mode.is_dual_plane)
		{
			weight = select(weight, vint4(cache.weights[1][i]), plane2_mask);
		}

		vint4 color = endpoints[p] + asr<6>(deltas[p] * weight + vint4(32));
		vfloat4 error = original.texel(i) - rdo_decode_ldr_color_4(color) * 255.0f;
		cache.texel_errors[i] = dot_s(error * error, original.channel_weight);
	}

	vfloat weight_r(original.channel_weight.lane<0>());
	vfloat weight_g(original.channel_weight.lane<1>());
	vfloat weight_b(original.channel_weight.lane<2>());
	vfloat weight_a(original.channel_weight.lane<3>());

	vfloat mean_r(original.data_mean.lane<0>());
	vfloat mean_g(original.data_mean.lane<1>());
	vfloat mean_b(original.data_mean.lane<2>());
	vfloat mean_a(original.data_mean.lane<3>());

	vfloatacc squared_error = vfloatacc::zero();
	vfloatacc variance = vfloatacc::zero();

	// Image edges use the source loader's replicated texels across the complete block footprint
	for (unsigned int i = 0; i < bsd.texel_count; i += ASTCENC_SIMD_WIDTH)
	{
		vmask active = (vint(i) + vint::lane_id()) < vint(bsd.texel_count);
		haccumulate(squared_error, vfloat(cache.texel_errors + i), active);

		vfloat deviation_r = vfloat(original.data_r + i) - mean_r;
		vfloat deviation_g = vfloat(original.data_g + i) - mean_g;
		vfloat deviation_b = vfloat(original.data_b + i) - mean_b;
		vfloat deviation_a = vfloat(original.data_a + i) - mean_a;

		vfloat metric_r = deviation_r * deviation_r * weight_r;
		vfloat metric_g = deviation_g * deviation_g * weight_g;
		vfloat metric_b = deviation_b * deviation_b * weight_b;
		vfloat metric_a = deviation_a * deviation_a * weight_a;

		haccumulate(variance, max(max(metric_r, metric_b), max(metric_g, metric_a)), active);
	}

	float texel_count = static_cast<float>(bsd.texel_count);

	return {
		.mse = hadd_s(squared_error) / texel_count,
		.std_dev = std::sqrt(hadd_s(variance) / texel_count)
	};
}

/**
 * @brief Unpack trial endpoints and identify an exactly zero alpha error.
 *
 * @param[in]  config      The LDR compressor configuration.
 * @param[in]  block       The nonconstant symbolic encoding.
 * @param[in]  original    The prepared source texels and normalized channel weights.
 * @param[in]  cache       The initialized original reconstruction cache.
 * @param[out] endpoints   The prepared trial endpoints and spans.
 *
 * @return True if every endpoint and source texel has exactly opaque alpha.
 */
static bool rdo_prepare_trial_endpoints(
	const astcenc_config& config,
	const symbolic_compressed_block& block,
	const image_block& original,
	const rdo_block_difference_cache& cache,
	rdo_trial_endpoints& endpoints
) {
	bool is_rgb = (original.data_min.lane<3>() == 255.0f)
	           && (original.data_max.lane<3>() == 255.0f);
	for (unsigned int p = 0; p < block.partition_count; p++)
	{
		vint4 endpoint0;
		vint4 endpoint1;
		if (cache.is_colors_unchanged)
		{
			endpoint0 = vint4(cache.endpoints[p][0]);
			endpoint1 = vint4(cache.endpoints[p][1]);
		}
		else
		{
			bool is_rgb_lns;
			bool is_alpha_lns;
			[[maybe_unused]] bool is_valid = unpack_color_endpoints(
			    config.profile, block.color_formats[p], block.color_values[p],
			    is_rgb_lns, is_alpha_lns, endpoint0, endpoint1);
			assert(is_valid);
		}

		endpoints.endpoint0[p] = endpoint0;
		endpoints.deltas[p] = endpoint1 - endpoint0;
		is_rgb = is_rgb && (endpoint0.lane<3>() == 65535) && (endpoint1.lane<3>() == 65535);
	}

	return is_rgb;
}

/**
 * @brief Score a dual-plane RGB block using SIMD across texels.
 *
 * @param[in] bsd         The block size descriptor.
 * @param[in] block       The nonconstant symbolic encoding.
 * @param[in] original    The prepared source texels and normalized channel weights.
 * @param     mse_limit   The MSE bound, or infinity for a complete score.
 * @param[in] cache       The initialized original reconstruction cache.
 * @param[in] endpoints   The prepared trial endpoints and spans.
 *
 * @return Weighted MSE, or maximum MSE for losing trials.
 */
static ASTCENC_SIMD_INLINE float rdo_compute_difference_2plane_rgb(
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& block,
	const image_block& original,
	float mse_limit,
	const rdo_block_difference_cache& cache,
	const rdo_trial_endpoints& endpoints
) {
	const auto& mode = bsd.get_block_mode(block.block_mode);
	const auto& decimation = bsd.get_decimation_info(mode.decimation_mode);
	vint4 endpoint0 = endpoints.endpoint0[0];
	vint4 endpoint_delta = endpoints.deltas[0];

	vtable_64x8 table;
	if ((decimation.weight_count != bsd.texel_count) && !cache.is_weights_unchanged)
	{
		vtable_prepare(table, block.weights);
	}

	vmask plane2_r(block.plane2_component == 0);
	vmask plane2_g(block.plane2_component == 1);
	vmask plane2_b(block.plane2_component == 2);

	vint endpoint_r(endpoint0.lane<0>());
	vint endpoint_g(endpoint0.lane<1>());
	vint endpoint_b(endpoint0.lane<2>());

	vint delta_r(endpoint_delta.lane<0>());
	vint delta_g(endpoint_delta.lane<1>());
	vint delta_b(endpoint_delta.lane<2>());

	vfloat channel_weight_r(original.channel_weight.lane<0>());
	vfloat channel_weight_g(original.channel_weight.lane<1>());
	vfloat channel_weight_b(original.channel_weight.lane<2>());

	unsigned int texel_count = bsd.texel_count;
	size_t vector_end = round_up_to_simd_multiple_vla(texel_count);

	vfloatacc squared_error = vfloatacc::zero();
	float error_limit = mse_limit * static_cast<float>(texel_count);
	unsigned int i = 0;

	// Reconstruct only the next SIMD batch so losing trials can stop before a full block decode
	for (; i < vector_end; i += ASTCENC_SIMD_WIDTH)
	{
		vint weight = rdo_interpolate_weights(decimation, block.weights, table, i, 0, cache);
		vint weight2 = rdo_interpolate_weights(decimation, block.weights, table, i,
		                                       WEIGHTS_PLANE2_OFFSET, cache);

		// Reuse cached errors when a weight-only change leaves this entire batch unchanged
		bool is_batch_unchanged = cache.is_colors_unchanged;
		if (is_batch_unchanged)
		{
			is_batch_unchanged = all(weight == vint(cache.weights[0] + i));
			is_batch_unchanged = is_batch_unchanged && all(weight2 == vint(cache.weights[1] + i));
		}

		if (is_batch_unchanged)
		{
			haccumulate(squared_error, vfloat(cache.texel_errors + i));

			if (hadd_s(squared_error) >= error_limit)
			{
				return std::numeric_limits<float>::max();
			}

			continue;
		}

		vint weight_r = weight;
		vint weight_g = weight;
		vint weight_b = weight;

		weight_r = select(weight, weight2, plane2_r);
		weight_g = select(weight, weight2, plane2_g);
		weight_b = select(weight, weight2, plane2_b);

		vint data_r = endpoint_r;
		vint data_g = endpoint_g;
		vint data_b = endpoint_b;

		vint texel_delta_r = delta_r;
		vint texel_delta_g = delta_g;
		vint texel_delta_b = delta_b;

		vint color_r = data_r + asr<6>(texel_delta_r * weight_r + vint(32));
		vint color_g = data_g + asr<6>(texel_delta_g * weight_g + vint(32));
		vint color_b = data_b + asr<6>(texel_delta_b * weight_b + vint(32));

		vfloat error_r = loada(original.data_r + i)
		               - rdo_decode_ldr_color(color_r) * 255.0f;
		vfloat error_g = loada(original.data_g + i)
		               - rdo_decode_ldr_color(color_g) * 255.0f;
		vfloat error_b = loada(original.data_b + i)
		               - rdo_decode_ldr_color(color_b) * 255.0f;

		// Preserve channel reduction order across all vector widths
		vfloat metric = error_r * error_r * channel_weight_r
		              + error_b * error_b * channel_weight_b
		              + error_g * error_g * channel_weight_g;

		// Exclude inactive texels from the invariant SIMD reduction
		if ((i + ASTCENC_SIMD_WIDTH) > texel_count)
		{
			vmask active = (vint(i) + vint::lane_id()) < vint(texel_count);
			metric = select(vfloat(0.0f), metric, active);
		}

		haccumulate(squared_error, metric);

		if (hadd_s(squared_error) >= error_limit)
		{
			return std::numeric_limits<float>::max();
		}
	}

	return hadd_s(squared_error) / static_cast<float>(texel_count);
}

/**
 * @brief Score a dual-plane RGBA block using SIMD across texels.
 *
 * @param[in] bsd         The block size descriptor.
 * @param[in] block       The nonconstant symbolic encoding.
 * @param[in] original    The prepared source texels and normalized channel weights.
 * @param     mse_limit   The MSE bound, or infinity for a complete score.
 * @param[in] cache       The initialized original reconstruction cache.
 * @param[in] endpoints   The prepared trial endpoints and spans.
 *
 * @return Weighted MSE, or maximum MSE for losing trials.
 */
static ASTCENC_SIMD_INLINE float rdo_compute_difference_2plane_rgba(
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& block,
	const image_block& original,
	float mse_limit,
	const rdo_block_difference_cache& cache,
	const rdo_trial_endpoints& endpoints
) {
	const auto& mode = bsd.get_block_mode(block.block_mode);
	const auto& decimation = bsd.get_decimation_info(mode.decimation_mode);
	vint4 endpoint0 = endpoints.endpoint0[0];
	vint4 endpoint_delta = endpoints.deltas[0];

	vtable_64x8 table;
	if ((decimation.weight_count != bsd.texel_count) && !cache.is_weights_unchanged)
	{
		vtable_prepare(table, block.weights);
	}

	vmask plane2_r(block.plane2_component == 0);
	vmask plane2_g(block.plane2_component == 1);
	vmask plane2_b(block.plane2_component == 2);
	vmask plane2_a(block.plane2_component == 3);

	vint endpoint_r(endpoint0.lane<0>());
	vint endpoint_g(endpoint0.lane<1>());
	vint endpoint_b(endpoint0.lane<2>());
	vint endpoint_a(endpoint0.lane<3>());

	vint delta_r(endpoint_delta.lane<0>());
	vint delta_g(endpoint_delta.lane<1>());
	vint delta_b(endpoint_delta.lane<2>());
	vint delta_a(endpoint_delta.lane<3>());

	vfloat channel_weight_r(original.channel_weight.lane<0>());
	vfloat channel_weight_g(original.channel_weight.lane<1>());
	vfloat channel_weight_b(original.channel_weight.lane<2>());
	vfloat channel_weight_a(original.channel_weight.lane<3>());

	unsigned int texel_count = bsd.texel_count;
	size_t vector_end = round_up_to_simd_multiple_vla(texel_count);

	vfloatacc squared_error = vfloatacc::zero();
	float error_limit = mse_limit * static_cast<float>(texel_count);
	unsigned int i = 0;

	// Reconstruct only the next SIMD batch so losing trials can stop before a full block decode
	for (; i < vector_end; i += ASTCENC_SIMD_WIDTH)
	{
		vint weight = rdo_interpolate_weights(decimation, block.weights, table, i, 0, cache);
		vint weight2 = rdo_interpolate_weights(decimation, block.weights, table, i,
		                                       WEIGHTS_PLANE2_OFFSET, cache);

		// Reuse cached errors when a weight-only change leaves this entire batch unchanged
		bool is_batch_unchanged = cache.is_colors_unchanged;
		if (is_batch_unchanged)
		{
			is_batch_unchanged = all(weight == vint(cache.weights[0] + i));
			is_batch_unchanged = is_batch_unchanged && all(weight2 == vint(cache.weights[1] + i));
		}

		if (is_batch_unchanged)
		{
			haccumulate(squared_error, vfloat(cache.texel_errors + i));

			if (hadd_s(squared_error) >= error_limit)
			{
				return std::numeric_limits<float>::max();
			}

			continue;
		}

		vint weight_r = weight;
		vint weight_g = weight;
		vint weight_b = weight;
		vint weight_a = weight;

		weight_r = select(weight, weight2, plane2_r);
		weight_g = select(weight, weight2, plane2_g);
		weight_b = select(weight, weight2, plane2_b);
		weight_a = select(weight, weight2, plane2_a);

		vint data_r = endpoint_r;
		vint data_g = endpoint_g;
		vint data_b = endpoint_b;
		vint data_a = endpoint_a;

		vint texel_delta_r = delta_r;
		vint texel_delta_g = delta_g;
		vint texel_delta_b = delta_b;
		vint texel_delta_a = delta_a;

		vint color_r = data_r + asr<6>(texel_delta_r * weight_r + vint(32));
		vint color_g = data_g + asr<6>(texel_delta_g * weight_g + vint(32));
		vint color_b = data_b + asr<6>(texel_delta_b * weight_b + vint(32));

		vfloat error_r = loada(original.data_r + i)
		               - rdo_decode_ldr_color(color_r) * 255.0f;
		vfloat error_g = loada(original.data_g + i)
		               - rdo_decode_ldr_color(color_g) * 255.0f;
		vfloat error_b = loada(original.data_b + i)
		               - rdo_decode_ldr_color(color_b) * 255.0f;
		vint color_a = data_a + asr<6>(texel_delta_a * weight_a + vint(32));

		vfloat error_a = loada(original.data_a + i) - rdo_decode_ldr_color(color_a) * 255.0f;

		// Preserve channel reduction order across all vector widths
		vfloat metric = (error_r * error_r * channel_weight_r + error_b * error_b * channel_weight_b)
		              + (error_g * error_g * channel_weight_g + error_a * error_a * channel_weight_a);

		// Exclude inactive texels from the invariant SIMD reduction
		if ((i + ASTCENC_SIMD_WIDTH) > texel_count)
		{
			vmask active = (vint(i) + vint::lane_id()) < vint(texel_count);
			metric = select(vfloat(0.0f), metric, active);
		}

		haccumulate(squared_error, metric);

		if (hadd_s(squared_error) >= error_limit)
		{
			return std::numeric_limits<float>::max();
		}
	}

	return hadd_s(squared_error) / static_cast<float>(texel_count);
}

/**
 * @brief Score a multipartition RGB block using SIMD across texels.
 *
 * @param[in] bsd         The block size descriptor.
 * @param[in] block       The nonconstant symbolic encoding.
 * @param[in] original    The prepared source texels and normalized channel weights.
 * @param     mse_limit   The MSE bound, or infinity for a complete score.
 * @param[in] cache       The initialized original reconstruction cache.
 * @param[in] endpoints   The prepared trial endpoints and spans.
 *
 * @return Weighted MSE, or maximum MSE for losing trials.
 */
static ASTCENC_SIMD_INLINE float rdo_compute_difference_1plane_rgb(
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& block,
	const image_block& original,
	float mse_limit,
	const rdo_block_difference_cache& cache,
	const rdo_trial_endpoints& endpoints
) {
	const auto& mode = bsd.get_block_mode(block.block_mode);
	const auto& decimation = bsd.get_decimation_info(mode.decimation_mode);
	vint4 endpoint0 = endpoints.endpoint0[0];
	vint4 endpoint_delta = endpoints.deltas[0];

	vtable_64x8 table;
	if ((decimation.weight_count != bsd.texel_count) && !cache.is_weights_unchanged)
	{
		vtable_prepare(table, block.weights);
	}

	vint endpoint_r(endpoint0.lane<0>());
	vint endpoint_g(endpoint0.lane<1>());
	vint endpoint_b(endpoint0.lane<2>());

	vint delta_r(endpoint_delta.lane<0>());
	vint delta_g(endpoint_delta.lane<1>());
	vint delta_b(endpoint_delta.lane<2>());

	vfloat channel_weight_r(original.channel_weight.lane<0>());
	vfloat channel_weight_g(original.channel_weight.lane<1>());
	vfloat channel_weight_b(original.channel_weight.lane<2>());

	unsigned int texel_count = bsd.texel_count;
	size_t vector_end = round_up_to_simd_multiple_vla(texel_count);

	vfloatacc squared_error = vfloatacc::zero();
	float error_limit = mse_limit * static_cast<float>(texel_count);
	unsigned int i = 0;

	// Reconstruct only the next SIMD batch so losing trials can stop before a full block decode
	for (; i < vector_end; i += ASTCENC_SIMD_WIDTH)
	{
		vint weight = rdo_interpolate_weights(decimation, block.weights, table, i, 0, cache);

		// Reuse cached errors when a weight-only change leaves this entire batch unchanged
		bool is_batch_unchanged = cache.is_colors_unchanged;
		if (is_batch_unchanged)
		{
			is_batch_unchanged = all(weight == vint(cache.weights[0] + i));
		}

		if (is_batch_unchanged)
		{
			haccumulate(squared_error, vfloat(cache.texel_errors + i));

			if (hadd_s(squared_error) >= error_limit)
			{
				return std::numeric_limits<float>::max();
			}

			continue;
		}

		vint weight_r = weight;
		vint weight_g = weight;
		vint weight_b = weight;

		vint data_r = endpoint_r;
		vint data_g = endpoint_g;
		vint data_b = endpoint_b;

		vint texel_delta_r = delta_r;
		vint texel_delta_g = delta_g;
		vint texel_delta_b = delta_b;

		const auto& partition = bsd.get_partition_info(block.partition_count, block.partition_index);
		vint partitions(partition.partition_of_texel + i);
		for (unsigned int p = 1; p < block.partition_count; p++)
		{
			vmask selected = partitions == vint(p);
			data_r = select(data_r, vint(endpoints.endpoint0[p].lane<0>()), selected);
			data_g = select(data_g, vint(endpoints.endpoint0[p].lane<1>()), selected);
			data_b = select(data_b, vint(endpoints.endpoint0[p].lane<2>()), selected);
			texel_delta_r = select(texel_delta_r, vint(endpoints.deltas[p].lane<0>()), selected);
			texel_delta_g = select(texel_delta_g, vint(endpoints.deltas[p].lane<1>()), selected);
			texel_delta_b = select(texel_delta_b, vint(endpoints.deltas[p].lane<2>()), selected);
		}

		vint color_r = data_r + asr<6>(texel_delta_r * weight_r + vint(32));
		vint color_g = data_g + asr<6>(texel_delta_g * weight_g + vint(32));
		vint color_b = data_b + asr<6>(texel_delta_b * weight_b + vint(32));

		vfloat error_r = loada(original.data_r + i)
		               - rdo_decode_ldr_color(color_r) * 255.0f;
		vfloat error_g = loada(original.data_g + i)
		               - rdo_decode_ldr_color(color_g) * 255.0f;
		vfloat error_b = loada(original.data_b + i)
		               - rdo_decode_ldr_color(color_b) * 255.0f;

		// Preserve channel reduction order across all vector widths
		vfloat metric = (error_r * error_r * channel_weight_r + error_b * error_b * channel_weight_b)
		              + error_g * error_g * channel_weight_g;
		// Exclude inactive texels from the invariant SIMD reduction
		if ((i + ASTCENC_SIMD_WIDTH) > texel_count)
		{
			vmask active = (vint(i) + vint::lane_id()) < vint(texel_count);
			metric = select(vfloat(0.0f), metric, active);
		}

		haccumulate(squared_error, metric);

		if (hadd_s(squared_error) >= error_limit)
		{
			return std::numeric_limits<float>::max();
		}
	}

	return hadd_s(squared_error) / static_cast<float>(texel_count);
}

/**
 * @brief Score a multipartition RGBA block using SIMD across texels.
 *
 * @param[in] bsd         The block size descriptor.
 * @param[in] block       The nonconstant symbolic encoding.
 * @param[in] original    The prepared source texels and normalized channel weights.
 * @param     mse_limit   The MSE bound, or infinity for a complete score.
 * @param[in] cache       The initialized original reconstruction cache.
 * @param[in] endpoints   The prepared trial endpoints and spans.
 *
 * @return Weighted MSE, or maximum MSE for losing trials.
 */
static ASTCENC_SIMD_INLINE float rdo_compute_difference_1plane_rgba(
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& block,
	const image_block& original,
	float mse_limit,
	const rdo_block_difference_cache& cache,
	const rdo_trial_endpoints& endpoints
) {
	const auto& mode = bsd.get_block_mode(block.block_mode);
	const auto& decimation = bsd.get_decimation_info(mode.decimation_mode);
	vint4 endpoint0 = endpoints.endpoint0[0];
	vint4 endpoint_delta = endpoints.deltas[0];

	vtable_64x8 table;
	if ((decimation.weight_count != bsd.texel_count) && !cache.is_weights_unchanged)
	{
		vtable_prepare(table, block.weights);
	}

	vint endpoint_r(endpoint0.lane<0>());
	vint endpoint_g(endpoint0.lane<1>());
	vint endpoint_b(endpoint0.lane<2>());
	vint endpoint_a(endpoint0.lane<3>());

	vint delta_r(endpoint_delta.lane<0>());
	vint delta_g(endpoint_delta.lane<1>());
	vint delta_b(endpoint_delta.lane<2>());
	vint delta_a(endpoint_delta.lane<3>());

	vfloat channel_weight_r(original.channel_weight.lane<0>());
	vfloat channel_weight_g(original.channel_weight.lane<1>());
	vfloat channel_weight_b(original.channel_weight.lane<2>());
	vfloat channel_weight_a(original.channel_weight.lane<3>());

	unsigned int texel_count = bsd.texel_count;
	size_t vector_end = round_up_to_simd_multiple_vla(texel_count);

	vfloatacc squared_error = vfloatacc::zero();
	float error_limit = mse_limit * static_cast<float>(texel_count);
	unsigned int i = 0;

	// Reconstruct only the next SIMD batch so losing trials can stop before a full block decode
	for (; i < vector_end; i += ASTCENC_SIMD_WIDTH)
	{
		vint weight = rdo_interpolate_weights(decimation, block.weights, table, i, 0, cache);

		// Reuse cached errors when a weight-only change leaves this entire batch unchanged
		bool is_batch_unchanged = cache.is_colors_unchanged;
		if (is_batch_unchanged)
		{
			is_batch_unchanged = all(weight == vint(cache.weights[0] + i));
		}

		if (is_batch_unchanged)
		{
			haccumulate(squared_error, vfloat(cache.texel_errors + i));

			if (hadd_s(squared_error) >= error_limit)
			{
				return std::numeric_limits<float>::max();
			}

			continue;
		}

		vint weight_r = weight;
		vint weight_g = weight;
		vint weight_b = weight;
		vint weight_a = weight;

		vint data_r = endpoint_r;
		vint data_g = endpoint_g;
		vint data_b = endpoint_b;
		vint data_a = endpoint_a;

		vint texel_delta_r = delta_r;
		vint texel_delta_g = delta_g;
		vint texel_delta_b = delta_b;
		vint texel_delta_a = delta_a;

		const auto& partition = bsd.get_partition_info(block.partition_count, block.partition_index);
		vint partitions(partition.partition_of_texel + i);
		for (unsigned int p = 1; p < block.partition_count; p++)
		{
			vmask selected = partitions == vint(p);
			data_r = select(data_r, vint(endpoints.endpoint0[p].lane<0>()), selected);
			data_g = select(data_g, vint(endpoints.endpoint0[p].lane<1>()), selected);
			data_b = select(data_b, vint(endpoints.endpoint0[p].lane<2>()), selected);
			texel_delta_r = select(texel_delta_r, vint(endpoints.deltas[p].lane<0>()), selected);
			texel_delta_g = select(texel_delta_g, vint(endpoints.deltas[p].lane<1>()), selected);
			texel_delta_b = select(texel_delta_b, vint(endpoints.deltas[p].lane<2>()), selected);
			data_a = select(data_a, vint(endpoints.endpoint0[p].lane<3>()), selected);
			texel_delta_a = select(texel_delta_a, vint(endpoints.deltas[p].lane<3>()), selected);
		}

		vint color_r = data_r + asr<6>(texel_delta_r * weight_r + vint(32));
		vint color_g = data_g + asr<6>(texel_delta_g * weight_g + vint(32));
		vint color_b = data_b + asr<6>(texel_delta_b * weight_b + vint(32));

		vfloat error_r = loada(original.data_r + i)
		               - rdo_decode_ldr_color(color_r) * 255.0f;
		vfloat error_g = loada(original.data_g + i)
		               - rdo_decode_ldr_color(color_g) * 255.0f;
		vfloat error_b = loada(original.data_b + i)
		               - rdo_decode_ldr_color(color_b) * 255.0f;
		vint color_a = data_a + asr<6>(texel_delta_a * weight_a + vint(32));

		vfloat error_a = loada(original.data_a + i) - rdo_decode_ldr_color(color_a) * 255.0f;

		// Preserve channel reduction order across all vector widths
		vfloat metric = (error_r * error_r * channel_weight_r + error_b * error_b * channel_weight_b)
		              + (error_g * error_g * channel_weight_g + error_a * error_a * channel_weight_a);
		// Exclude inactive texels from the invariant SIMD reduction
		if ((i + ASTCENC_SIMD_WIDTH) > texel_count)
		{
			vmask active = (vint(i) + vint::lane_id()) < vint(texel_count);
			metric = select(vfloat(0.0f), metric, active);
		}

		haccumulate(squared_error, metric);

		if (hadd_s(squared_error) >= error_limit)
		{
			return std::numeric_limits<float>::max();
		}
	}

	return hadd_s(squared_error) / static_cast<float>(texel_count);
}

/**
 * @brief Score a single-partition RGB block using SIMD across texels.
 *
 * @param[in] bsd         The block size descriptor.
 * @param[in] block       The nonconstant symbolic encoding.
 * @param[in] original    The prepared source texels and normalized channel weights.
 * @param     mse_limit   The MSE bound, or infinity for a complete score.
 * @param[in] cache       The initialized original reconstruction cache.
 * @param[in] endpoints   The prepared trial endpoints and spans.
 *
 * @return Weighted MSE, or maximum MSE for losing trials.
 */
static ASTCENC_SIMD_INLINE float rdo_compute_difference_1plane_1partition_rgb(
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& block,
	const image_block& original,
	float mse_limit,
	const rdo_block_difference_cache& cache,
	const rdo_trial_endpoints& endpoints
) {
	const auto& mode = bsd.get_block_mode(block.block_mode);
	const auto& decimation = bsd.get_decimation_info(mode.decimation_mode);
	vint4 endpoint0 = endpoints.endpoint0[0];
	vint4 endpoint_delta = endpoints.deltas[0];

	vtable_64x8 table;
	if ((decimation.weight_count != bsd.texel_count) && !cache.is_weights_unchanged)
	{
		vtable_prepare(table, block.weights);
	}

	vint endpoint_r(endpoint0.lane<0>());
	vint endpoint_g(endpoint0.lane<1>());
	vint endpoint_b(endpoint0.lane<2>());

	vint delta_r(endpoint_delta.lane<0>());
	vint delta_g(endpoint_delta.lane<1>());
	vint delta_b(endpoint_delta.lane<2>());

	vfloat channel_weight_r(original.channel_weight.lane<0>());
	vfloat channel_weight_g(original.channel_weight.lane<1>());
	vfloat channel_weight_b(original.channel_weight.lane<2>());

	unsigned int texel_count = bsd.texel_count;
	size_t vector_end = round_up_to_simd_multiple_vla(texel_count);

	vfloatacc squared_error = vfloatacc::zero();
	float error_limit = mse_limit * static_cast<float>(texel_count);
	unsigned int i = 0;

	// Reconstruct only the next SIMD batch so losing trials can stop before a full block decode
	for (; i < vector_end; i += ASTCENC_SIMD_WIDTH)
	{
		vint weight = rdo_interpolate_weights(decimation, block.weights, table, i, 0, cache);

		// Reuse cached errors when a weight-only change leaves this entire batch unchanged
		bool is_batch_unchanged = cache.is_colors_unchanged;
		if (is_batch_unchanged)
		{
			is_batch_unchanged = all(weight == vint(cache.weights[0] + i));
		}

		if (is_batch_unchanged)
		{
			haccumulate(squared_error, vfloat(cache.texel_errors + i));

			if (hadd_s(squared_error) >= error_limit)
			{
				return std::numeric_limits<float>::max();
			}

			continue;
		}

		vint weight_r = weight;
		vint weight_g = weight;
		vint weight_b = weight;

		vint data_r = endpoint_r;
		vint data_g = endpoint_g;
		vint data_b = endpoint_b;

		vint texel_delta_r = delta_r;
		vint texel_delta_g = delta_g;
		vint texel_delta_b = delta_b;

		vint color_r = data_r + asr<6>(texel_delta_r * weight_r + vint(32));
		vint color_g = data_g + asr<6>(texel_delta_g * weight_g + vint(32));
		vint color_b = data_b + asr<6>(texel_delta_b * weight_b + vint(32));

		vfloat error_r = loada(original.data_r + i)
		               - rdo_decode_ldr_color(color_r) * 255.0f;
		vfloat error_g = loada(original.data_g + i)
		               - rdo_decode_ldr_color(color_g) * 255.0f;
		vfloat error_b = loada(original.data_b + i)
		               - rdo_decode_ldr_color(color_b) * 255.0f;

		// Preserve channel reduction order across all vector widths
		vfloat metric = (error_r * error_r * channel_weight_r + error_b * error_b * channel_weight_b)
		              + error_g * error_g * channel_weight_g;
		// Exclude inactive texels from the invariant SIMD reduction
		if ((i + ASTCENC_SIMD_WIDTH) > texel_count)
		{
			vmask active = (vint(i) + vint::lane_id()) < vint(texel_count);
			metric = select(vfloat(0.0f), metric, active);
		}

		haccumulate(squared_error, metric);

		if (hadd_s(squared_error) >= error_limit)
		{
			return std::numeric_limits<float>::max();
		}
	}

	return hadd_s(squared_error) / static_cast<float>(texel_count);
}

/**
 * @brief Score a single-partition RGBA block using SIMD across texels.
 *
 * @param[in] bsd         The block size descriptor.
 * @param[in] block       The nonconstant symbolic encoding.
 * @param[in] original    The prepared source texels and normalized channel weights.
 * @param     mse_limit   The MSE bound, or infinity for a complete score.
 * @param[in] cache       The initialized original reconstruction cache.
 * @param[in] endpoints   The prepared trial endpoints and spans.
 *
 * @return Weighted MSE, or maximum MSE for losing trials.
 */
static ASTCENC_SIMD_INLINE float rdo_compute_difference_1plane_1partition_rgba(
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& block,
	const image_block& original,
	float mse_limit,
	const rdo_block_difference_cache& cache,
	const rdo_trial_endpoints& endpoints
) {
	const auto& mode = bsd.get_block_mode(block.block_mode);
	const auto& decimation = bsd.get_decimation_info(mode.decimation_mode);
	vint4 endpoint0 = endpoints.endpoint0[0];
	vint4 endpoint_delta = endpoints.deltas[0];

	vtable_64x8 table;
	if ((decimation.weight_count != bsd.texel_count) && !cache.is_weights_unchanged)
	{
		vtable_prepare(table, block.weights);
	}

	vint endpoint_r(endpoint0.lane<0>());
	vint endpoint_g(endpoint0.lane<1>());
	vint endpoint_b(endpoint0.lane<2>());
	vint endpoint_a(endpoint0.lane<3>());

	vint delta_r(endpoint_delta.lane<0>());
	vint delta_g(endpoint_delta.lane<1>());
	vint delta_b(endpoint_delta.lane<2>());
	vint delta_a(endpoint_delta.lane<3>());

	vfloat channel_weight_r(original.channel_weight.lane<0>());
	vfloat channel_weight_g(original.channel_weight.lane<1>());
	vfloat channel_weight_b(original.channel_weight.lane<2>());
	vfloat channel_weight_a(original.channel_weight.lane<3>());

	unsigned int texel_count = bsd.texel_count;
	size_t vector_end = round_up_to_simd_multiple_vla(texel_count);

	vfloatacc squared_error = vfloatacc::zero();
	float error_limit = mse_limit * static_cast<float>(texel_count);
	unsigned int i = 0;

	// Reconstruct only the next SIMD batch so losing trials can stop before a full block decode
	for (; i < vector_end; i += ASTCENC_SIMD_WIDTH)
	{
		vint weight = rdo_interpolate_weights(decimation, block.weights, table, i, 0, cache);

		// Reuse cached errors when a weight-only change leaves this entire batch unchanged
		bool is_batch_unchanged = cache.is_colors_unchanged;
		if (is_batch_unchanged)
		{
			is_batch_unchanged = all(weight == vint(cache.weights[0] + i));
		}

		if (is_batch_unchanged)
		{
			haccumulate(squared_error, vfloat(cache.texel_errors + i));

			if (hadd_s(squared_error) >= error_limit)
			{
				return std::numeric_limits<float>::max();
			}

			continue;
		}

		vint weight_r = weight;
		vint weight_g = weight;
		vint weight_b = weight;
		vint weight_a = weight;

		vint data_r = endpoint_r;
		vint data_g = endpoint_g;
		vint data_b = endpoint_b;
		vint data_a = endpoint_a;

		vint texel_delta_r = delta_r;
		vint texel_delta_g = delta_g;
		vint texel_delta_b = delta_b;
		vint texel_delta_a = delta_a;

		vint color_r = data_r + asr<6>(texel_delta_r * weight_r + vint(32));
		vint color_g = data_g + asr<6>(texel_delta_g * weight_g + vint(32));
		vint color_b = data_b + asr<6>(texel_delta_b * weight_b + vint(32));

		vfloat error_r = loada(original.data_r + i)
		               - rdo_decode_ldr_color(color_r) * 255.0f;
		vfloat error_g = loada(original.data_g + i)
		               - rdo_decode_ldr_color(color_g) * 255.0f;
		vfloat error_b = loada(original.data_b + i)
		               - rdo_decode_ldr_color(color_b) * 255.0f;
		vint color_a = data_a + asr<6>(texel_delta_a * weight_a + vint(32));

		vfloat error_a = loada(original.data_a + i) - rdo_decode_ldr_color(color_a) * 255.0f;

		// Preserve channel reduction order across all vector widths
		vfloat metric = (error_r * error_r * channel_weight_r + error_b * error_b * channel_weight_b)
		              + (error_g * error_g * channel_weight_g + error_a * error_a * channel_weight_a);
		// Exclude inactive texels from the invariant SIMD reduction
		if ((i + ASTCENC_SIMD_WIDTH) > texel_count)
		{
			vmask active = (vint(i) + vint::lane_id()) < vint(texel_count);
			metric = select(vfloat(0.0f), metric, active);
		}

		haccumulate(squared_error, metric);

		if (hadd_s(squared_error) >= error_limit)
		{
			return std::numeric_limits<float>::max();
		}
	}

	return hadd_s(squared_error) / static_cast<float>(texel_count);
}

/* See header for documentation. */
float rdo_compute_block_difference(
	const astcenc_config& config,
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& block,
	const image_block& original,
	float mse_limit,
	const rdo_block_difference_cache& cache
) {
	assert(block.block_type == SYM_BTYPE_NONCONST);
	rdo_trial_endpoints endpoints;
	bool is_rgb = rdo_prepare_trial_endpoints(config, block, original, cache, endpoints);

	// The core compressor emits dual-plane blocks only with a single partition
	if (block.partition_count > 1)
	{
		assert(!bsd.get_block_mode(block.block_mode).is_dual_plane);
		if (is_rgb)
		{
			return rdo_compute_difference_1plane_rgb(
			    bsd, block, original, mse_limit, cache, endpoints);
		}

		return rdo_compute_difference_1plane_rgba(
		    bsd, block, original, mse_limit, cache, endpoints);
	}

	const auto& mode = bsd.get_block_mode(block.block_mode);
	if (mode.is_dual_plane)
	{
		if (is_rgb)
		{
			return rdo_compute_difference_2plane_rgb(
			    bsd, block, original, mse_limit, cache, endpoints);
		}

		return rdo_compute_difference_2plane_rgba(
		    bsd, block, original, mse_limit, cache, endpoints);
	}

	if (is_rgb)
	{
		return rdo_compute_difference_1plane_1partition_rgb(
		    bsd, block, original, mse_limit, cache, endpoints);
	}

	return rdo_compute_difference_1plane_1partition_rgba(
	    bsd, block, original, mse_limit, cache, endpoints);
}

#endif
