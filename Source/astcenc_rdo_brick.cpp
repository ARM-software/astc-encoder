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
 * @brief Brick-local RDO integration for LDR image compression.
 */

#if !defined(ASTCENC_DECOMPRESS_ONLY)

#include <cmath>
#include <cstring>

#include "astcenc_internal.h"
#include "astcenc_rdo.h"

/** @brief Source data and scratch storage for one brick's distortion scoring. */
struct rdo_brick_context
{
	/** @brief The compressor context. */
	const astcenc_contexti& ctx;

	/** @brief The original image. */
	const astcenc_image& image;

	/** @brief The compression swizzle. */
	const astcenc_swizzle& swizzle;

	/** @brief Raster byte offsets indexed by compact brick position. */
	const size_t* offsets;

	/** @brief Original symbolic blocks indexed by compact brick position. */
	const symbolic_compressed_block* const* blocks;

	/** @brief The number of blocks in an image row. */
	size_t blocks_x;

	/** @brief Cached source texels in the 0..255 distortion domain, loaded once per block. */
	image_block original;

	/** @brief Reusable trial block retaining the original layout. */
	symbolic_compressed_block symbolic;

	/** @brief ISE layout prepared once for the current nonconstant block. */
	rdo_block_decode_state decode_state;

	/** @brief Original endpoint and interpolation results reused by unchanged payloads. */
	rdo_block_difference_cache difference_cache;
};

/**
 * @brief Test whether RDO supports and enables this configuration.
 *
 * @param[in] config   The compressor configuration.
 *
 * @return True if RDO is enabled for this configuration.
 */
static bool rdo_is_enabled(
	const astcenc_config& config
) {
	bool is_ldr = (config.profile == ASTCENC_PRF_LDR) ||
	              (config.profile == ASTCENC_PRF_LDR_SRGB);

	bool is_compression_enabled = !(config.flags & ASTCENC_FLG_DECOMPRESS_ONLY);

	return (config.rdo_lambda > 0.0f)
	    && is_ldr && (config.block_z == 1)
	    && is_compression_enabled;
}

/**
 * @brief Scale speculative search to the block's effective linear dimensions.
 *
 * This empirical fit targets about 0.3% combined DEFLATE growth on KodakMnt RGB at lambda 0.2.
 * Rate estimates still retain the full sliding history.
 *
 * @param block_x   The ASTC block width in texels.
 * @param block_y   The ASTC block height in texels.
 *
 * @return An even predecessor count within the configured search bounds.
 */
static unsigned int rdo_get_search_lookback(
	unsigned int block_x,
	unsigned int block_y
) {
	unsigned int texel_count = block_x * block_y;

	unsigned int perimeter_budget = rdo_search_dimension_scale * (block_x + block_y);

	unsigned int search_size = rdo_search_base_lookback
	                         + 2 * ((perimeter_budget + texel_count - 1) / texel_count);

	return astc::clamp(search_size, rdo_min_lookback, rdo_max_lookback);
}

/**
 * @brief Load a block's source data and normalized distortion weights.
 *
 * @param[in,out] ctx           The brick distortion context.
 * @param         block_index   The compact block index in brick traversal order.
 */
static void rdo_load_source(
	rdo_brick_context& ctx,
	unsigned int block_index
) {
	const auto& config = ctx.ctx.config;
	const auto& bsd = *ctx.ctx.bsd;
	size_t raster_index = ctx.offsets[block_index] / 16;

	const position origin {
		.x = (raster_index % ctx.blocks_x) * bsd.dim_x,
		.y = (raster_index / ctx.blocks_x) * bsd.dim_y,
		.z = 0
	};

	// Normalize source data once; every trial for this block uses the same distortion domain
	image_block& original = ctx.original;
	original.texel_count = bsd.texel_count;

	load_image_block(config.profile, ctx.image, original, bsd,
	                 origin.x, origin.y, origin.z, ctx.swizzle);

	vfloat4 weights(config.cw_r_weight, config.cw_g_weight,
	                config.cw_b_weight, config.cw_a_weight);

	if (config.flags & ASTCENC_FLG_USE_ALPHA_WEIGHT)
	{
		float alpha_scale = original.data_max.lane<3>() * (1.0f / 65535.0f);
		weights = weights * vfloat4(alpha_scale, alpha_scale, alpha_scale, 1.0f);
	}

	original.channel_weight = weights / hadd_s(weights);
	rdo_prepare_difference_source(original, bsd.texel_count);

	// Protected layout bits let trials reuse the encoder's mode and ISE sequence boundaries
	const auto& source = *ctx.blocks[block_index];
	assert(source.block_type == SYM_BTYPE_NONCONST);
	ctx.symbolic = source;
	ctx.decode_state = rdo_prepare_block_decode(bsd, source);

	// Initialize only the padding needed by full SIMD reads of the original weight planes
	size_t weight_count = ctx.decode_state.weight_count;
	size_t padding = round_up_to_simd_multiple_vla(weight_count) - weight_count;
	if (padding != 0)
	{
		std::memset(ctx.symbolic.weights + weight_count, 0, padding);
		if (ctx.decode_state.is_dual_plane)
		{
			std::memset(ctx.symbolic.weights + WEIGHTS_PLANE2_OFFSET + weight_count, 0, padding);
		}
	}
}

/* See header for documentation. */
rdo_block_statistics rdo_initialize_block(
	rdo_brick_context& ctx,
	unsigned int block_index
) {
	rdo_load_source(ctx, block_index);
	const auto& bsd = *ctx.ctx.bsd;
	rdo_prepare_difference_cache(
	    ctx.ctx.config, bsd, ctx.symbolic, ctx.difference_cache);

	return rdo_compute_original_difference(bsd, ctx.symbolic, ctx.original,
	                                      ctx.difference_cache);
}

/* See header for documentation. */
float rdo_block_difference(
	rdo_brick_context& ctx,
	const uint8_t* physical,
	unsigned int byte_mask,
	float mse_limit
) {
	const auto& bsd = *ctx.ctx.bsd;
	auto& symbolic = ctx.symbolic;
	rdo_decode_changed_payload(ctx.decode_state, physical, byte_mask, symbolic);

	ctx.difference_cache.is_weights_unchanged = !(byte_mask & ctx.decode_state.weight_byte_mask);
	ctx.difference_cache.is_colors_unchanged = !(byte_mask & ctx.decode_state.color_byte_mask);

	return rdo_compute_block_difference(ctx.ctx.config, bsd, symbolic, ctx.original,
	                                    mse_limit, ctx.difference_cache);
}

/* See header for documentation. */
void rdo_optimize_physical_brick(
	const astcenc_contexti& ctx,
	const astcenc_image& image,
	const astcenc_swizzle& swizzle,
	const brick_layout::validity_mask& valid_mask,
	const size_t* offsets,
	const symbolic_compressed_block* blocks,
	uint8_t* buffer
) {
	if (!rdo_is_enabled(ctx.config) || (image.dim_z != 1))
	{
		return;
	}

	// Compact valid blocks in the encoder brick's traversal order, without padding
	uint8_t physical[brick_layout::block_count * 16];
	uint8_t masks[brick_layout::block_count * 16];
	size_t compact_offsets[brick_layout::block_count];
	const symbolic_compressed_block* compact_blocks[brick_layout::block_count];
	unsigned int block_count = 0;

	for (unsigned int index = 0; index < brick_layout::block_count; index++)
	{
		if (!valid_mask[index])
		{
			continue;
		}

		compact_offsets[block_count] = offsets[index];
		compact_blocks[block_count] = blocks + index;
		std::memcpy(physical + block_count * 16, buffer + offsets[index], 16);
		rdo_get_structure_mask(*ctx.bsd, blocks[index], masks + block_count * 16);
		block_count++;
	}

	if (block_count < 2)
	{
		return;
	}

	rdo_brick_context context {
		.ctx = ctx,
		.image = image,
		.swizzle = swizzle,
		.offsets = compact_offsets,
		.blocks = compact_blocks,
		.blocks_x = astc::get_block_count_safe(image.dim_x, ctx.bsd->dim_x),
		.original = {},
		.symbolic = {},
		.decode_state = {},
		.difference_cache = {}
	};

	const rdo_parameters params {
		.lambda = ctx.config.rdo_lambda * rdo_lambda_scale,
		.search_lookback = rdo_get_search_lookback(ctx.bsd->dim_x, ctx.bsd->dim_y)
	};

	rdo_reduce_entropy(context, physical, block_count, masks, params);

	// Repacking through symbolic form could discard byte matches in ignored physical bits
	for (unsigned int i = 0; i < block_count; i++)
	{
		std::memcpy(buffer + compact_offsets[i], physical + i * 16, 16);
	}
}

#endif
