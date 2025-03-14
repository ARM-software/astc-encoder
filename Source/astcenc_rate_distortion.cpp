// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2011-2024 Arm Limited
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

#if !defined(ASTCENC_DECOMPRESS_ONLY)

#include "astcenc_internal_entry.h"
#include "ert.h"

#include <vector>

struct astcenc_rdo_context
{
	ert::reduce_entropy_params m_ert_params;
	std::vector<image_block> m_blocks;

	uint32_t m_image_x = 0;
	uint32_t m_image_y = 0;
	uint32_t m_image_z = 0;
};

#define ASTCENC_RDO_SPECIALIZE_DIFF 1

static constexpr uint32_t ASTCENC_BYTES_PER_BLOCK = 16;

template<typename T> T sqr(T v) { return v * v; }

extern "C" void rdo_progress_emitter(
	float value
) {
	static float previous_value = 100.0f;
	if (previous_value == 100.0f)
	{
		printf("\n\n");
		printf("Rate-distortion optimization\n");
		printf("============================\n\n");
	}
	previous_value = value;

	const unsigned int bar_size = 25;
	unsigned int parts = static_cast<int>(value / 4.0f);

	char buffer[bar_size + 3];
	buffer[0] = '[';

	for (unsigned int i = 0; i < parts; i++)
	{
		buffer[i + 1] = '=';
	}

	for (unsigned int i = parts; i < bar_size; i++)
	{
		buffer[i + 1] = ' ';
	}

	buffer[bar_size + 1] = ']';
	buffer[bar_size + 2] = '\0';

	printf("    Progress: %s %03.1f%%\r", buffer, static_cast<double>(value));
	fflush(stdout);
}

static uint32_t init_rdo_context(
		astcenc_contexti& ctx,
		const astcenc_image& image,
		const astcenc_swizzle& swz
) {
	ctx.rdo_context = new astcenc_rdo_context;
	astcenc_rdo_context& rdo_ctx = *ctx.rdo_context;

	uint32_t block_dim_x = ctx.bsd->xdim;
	uint32_t block_dim_y = ctx.bsd->ydim;
	uint32_t block_dim_z = ctx.bsd->zdim;
	uint32_t xblocks = (image.dim_x + block_dim_x - 1u) / block_dim_x;
	uint32_t yblocks = (image.dim_y + block_dim_y - 1u) / block_dim_y;
	uint32_t zblocks = (image.dim_z + block_dim_z - 1u) / block_dim_z;
	uint32_t total_blocks = xblocks * yblocks * zblocks;

	// Generate quality parameters
	auto& ert_params = rdo_ctx.m_ert_params;
	ert_params.m_lambda = ctx.config.rdo_quality;
	ert_params.m_lookback_window_size = ctx.config.rdo_lookback * ASTCENC_BYTES_PER_BLOCK;
	ert_params.m_smooth_block_max_mse_scale = ctx.config.rdo_max_smooth_block_error_scale;
	ert_params.m_max_smooth_block_std_dev = ctx.config.rdo_max_smooth_block_std_dev;

	ert_params.m_try_two_matches = true;

	rdo_ctx.m_blocks.resize(total_blocks);
	rdo_ctx.m_image_x = image.dim_x;
	rdo_ctx.m_image_y = image.dim_y;
	rdo_ctx.m_image_z = image.dim_z;

	vfloat4 channel_weight = vfloat4(ctx.config.cw_r_weight,
									 ctx.config.cw_g_weight,
									 ctx.config.cw_b_weight,
									 ctx.config.cw_a_weight);
	channel_weight = channel_weight / hadd_s(channel_weight);

	for (uint32_t block_z = 0, block_idx = 0; block_z < zblocks; ++block_z)
	{
		for (uint32_t block_y = 0; block_y < yblocks; ++block_y)
		{
			for (uint32_t block_x = 0; block_x < xblocks; ++block_x, ++block_idx)
			{
				image_block& blk = rdo_ctx.m_blocks[block_idx];
				blk.decode_unorm8 = ctx.config.flags & ASTCENC_FLG_USE_DECODE_UNORM8;
				blk.texel_count = ctx.bsd->texel_count;

				load_image_block(ctx.config.profile, image, blk, *ctx.bsd,
								 block_dim_x * block_x, block_dim_y * block_y, block_dim_z * block_z, swz);

				if (ctx.config.flags & ASTCENC_FLG_USE_ALPHA_WEIGHT)
				{
					float alpha_scale = blk.data_max.lane<3>() * (1.0f / 65535.0f);
					blk.channel_weight = vfloat4(ctx.config.cw_r_weight * alpha_scale,
												 ctx.config.cw_g_weight * alpha_scale,
												 ctx.config.cw_b_weight * alpha_scale,
												 ctx.config.cw_a_weight);
					blk.channel_weight = blk.channel_weight / hadd_s(blk.channel_weight);
				}
				else
				{
					blk.channel_weight = channel_weight;
				}
			}
		}
	}

	return total_blocks;
}

static float compute_block_std_dev(
	const image_block& blk,
	float scale
) {
	if (all(blk.data_min == blk.data_max)) return 0.0f;

	vfloatacc summav = vfloatacc::zero();
	vint lane_id = vint::lane_id();
	uint32_t texel_count = blk.texel_count;
	vfloat color_mean_r(blk.data_mean.lane<0>() * scale);
	vfloat color_mean_g(blk.data_mean.lane<1>() * scale);
	vfloat color_mean_b(blk.data_mean.lane<2>() * scale);
	vfloat color_mean_a(blk.data_mean.lane<3>() * scale);

	for (uint32_t i = 0; i < texel_count; i += ASTCENC_SIMD_WIDTH)
	{
		vfloat color_orig_r = loada(blk.data_r + i) * scale;
		vfloat color_orig_g = loada(blk.data_g + i) * scale;
		vfloat color_orig_b = loada(blk.data_b + i) * scale;
		vfloat color_orig_a = loada(blk.data_a + i) * scale;

		vfloat color_error_r = min(abs(color_orig_r - color_mean_r), vfloat(1e15f));
		vfloat color_error_g = min(abs(color_orig_g - color_mean_g), vfloat(1e15f));
		vfloat color_error_b = min(abs(color_orig_b - color_mean_b), vfloat(1e15f));
		vfloat color_error_a = min(abs(color_orig_a - color_mean_a), vfloat(1e15f));

		// Compute squared error metric
		color_error_r = color_error_r * color_error_r;
		color_error_g = color_error_g * color_error_g;
		color_error_b = color_error_b * color_error_b;
		color_error_a = color_error_a * color_error_a;

		vfloat metric = astc::max(color_error_r * blk.channel_weight.lane<0>(),
								  color_error_g * blk.channel_weight.lane<1>(),
								  color_error_b * blk.channel_weight.lane<2>(),
								  color_error_a * blk.channel_weight.lane<3>());

		// Mask off bad lanes
		vmask mask = lane_id < vint(texel_count);
		lane_id += vint(ASTCENC_SIMD_WIDTH);
		haccumulate(summav, metric, mask);
	}

	return sqrtf(hadd_s(summav) / texel_count);
}

#if ASTCENC_RDO_SPECIALIZE_DIFF
static float compute_symbolic_block_difference_constant(
		const astcenc_config& config,
		const block_size_descriptor& bsd,
		symbolic_compressed_block scb,
		const image_block& blk
) {
	vfloat4 color(0.0f);

	// UNORM16 constant color block
	if (scb.block_type == SYM_BTYPE_CONST_U16)
	{
		vint4 colori(scb.constant_color);

		// Determine the UNORM8 rounding on the decode
		vmask4 u8_mask = get_u8_component_mask(config.profile, blk);

		// The real decoder would just use the top 8 bits, but we rescale
		// in to a 16-bit value that rounds correctly.
		vint4 colori_u8 = asr<8>(colori) * 257;
		colori = select(colori, colori_u8, u8_mask);

		vint4 colorf16 = unorm16_to_sf16(colori);
		color = float16_to_float(colorf16);
	}
	// FLOAT16 constant color block
	else
	{
		switch (config.profile)
		{
			case ASTCENC_PRF_LDR_SRGB:
			case ASTCENC_PRF_LDR:
				return -ERROR_CALC_DEFAULT;
			case ASTCENC_PRF_HDR_RGB_LDR_A:
			case ASTCENC_PRF_HDR:
				// Constant-color block; unpack from FP16 to FP32.
				color = float16_to_float(vint4(scb.constant_color));
				break;
		}
	}

	if (all(blk.data_min == blk.data_max)) // Original block is also constant
	{
		vfloat4 color_error = min(abs(blk.origin_texel - color) * 65535.0f, vfloat4(1e15f));
		return dot_s(color_error * color_error, blk.channel_weight) * bsd.texel_count;
	}

	vfloatacc summav = vfloatacc::zero();
	vint lane_id = vint::lane_id();
	uint32_t texel_count = bsd.texel_count;

	vfloat color_r(color.lane<0>() * 65535.0f);
	vfloat color_g(color.lane<1>() * 65535.0f);
	vfloat color_b(color.lane<2>() * 65535.0f);
	vfloat color_a(color.lane<3>() * 65535.0f);

	for (uint32_t i = 0; i < texel_count; i += ASTCENC_SIMD_WIDTH)
	{
		vfloat color_orig_r = loada(blk.data_r + i);
		vfloat color_orig_g = loada(blk.data_g + i);
		vfloat color_orig_b = loada(blk.data_b + i);
		vfloat color_orig_a = loada(blk.data_a + i);

		vfloat color_error_r = min(abs(color_orig_r - color_r), vfloat(1e15f));
		vfloat color_error_g = min(abs(color_orig_g - color_g), vfloat(1e15f));
		vfloat color_error_b = min(abs(color_orig_b - color_b), vfloat(1e15f));
		vfloat color_error_a = min(abs(color_orig_a - color_a), vfloat(1e15f));

		// Compute squared error metric
		color_error_r = color_error_r * color_error_r;
		color_error_g = color_error_g * color_error_g;
		color_error_b = color_error_b * color_error_b;
		color_error_a = color_error_a * color_error_a;

		vfloat metric = color_error_r * blk.channel_weight.lane<0>()
					  + color_error_g * blk.channel_weight.lane<1>()
					  + color_error_b * blk.channel_weight.lane<2>()
					  + color_error_a * blk.channel_weight.lane<3>();

		// Mask off bad lanes
		vmask mask = lane_id < vint(texel_count);
		lane_id += vint(ASTCENC_SIMD_WIDTH);
		haccumulate(summav, metric, mask);
	}

	return hadd_s(summav);
}
#else
static float compute_block_mse(
	const image_block& orig,
	const image_block& cmp,
	const block_size_descriptor& bsd,
	uint32_t image_x,
	uint32_t image_y,
	uint32_t image_z,
	float orig_scale,
	float cmp_scale
) {
	vfloatacc summav = vfloatacc::zero();
	vint lane_id = vint::lane_id();
	uint32_t texel_count = orig.texel_count;

	uint32_t block_x = astc::min(image_x - orig.xpos, (uint32_t)bsd.xdim);
	uint32_t block_y = astc::min(image_y - orig.ypos, (uint32_t)bsd.ydim);
	uint32_t block_z = astc::min(image_z - orig.zpos, (uint32_t)bsd.zdim);

	for (uint32_t i = 0; i < texel_count; i += ASTCENC_SIMD_WIDTH, lane_id += vint(ASTCENC_SIMD_WIDTH))
	{
		vfloat color_orig_r = loada(orig.data_r + i) * orig_scale;
		vfloat color_orig_g = loada(orig.data_g + i) * orig_scale;
		vfloat color_orig_b = loada(orig.data_b + i) * orig_scale;
		vfloat color_orig_a = loada(orig.data_a + i) * orig_scale;

		vfloat color_cmp_r = loada(cmp.data_r + i) * cmp_scale;
		vfloat color_cmp_g = loada(cmp.data_g + i) * cmp_scale;
		vfloat color_cmp_b = loada(cmp.data_b + i) * cmp_scale;
		vfloat color_cmp_a = loada(cmp.data_a + i) * cmp_scale;

		vfloat color_error_r = min(abs(color_orig_r - color_cmp_r), vfloat(1e15f));
		vfloat color_error_g = min(abs(color_orig_g - color_cmp_g), vfloat(1e15f));
		vfloat color_error_b = min(abs(color_orig_b - color_cmp_b), vfloat(1e15f));
		vfloat color_error_a = min(abs(color_orig_a - color_cmp_a), vfloat(1e15f));

		// Compute squared error metric
		color_error_r = color_error_r * color_error_r;
		color_error_g = color_error_g * color_error_g;
		color_error_b = color_error_b * color_error_b;
		color_error_a = color_error_a * color_error_a;

		vfloat metric = color_error_r * orig.channel_weight.lane<0>()
					  + color_error_g * orig.channel_weight.lane<1>()
					  + color_error_b * orig.channel_weight.lane<2>()
					  + color_error_a * orig.channel_weight.lane<3>();

		// Mask off bad lanes
		vint lane_id_z(float_to_int(int_to_float(lane_id) / float(bsd.xdim * bsd.ydim)));
		vint rem_idx = lane_id - lane_id_z * vint(bsd.xdim * bsd.ydim);
		vint lane_id_y = float_to_int(int_to_float(rem_idx) / bsd.xdim);
		vint lane_id_x = rem_idx - lane_id_y * vint(bsd.xdim);
		vmask mask = (lane_id_x < vint(block_x)) & (lane_id_y < vint(block_y)) & (lane_id_z < vint(block_z));
		haccumulate(summav, metric, mask);
	}

	return hadd_s(summav) / (block_x * block_y * block_z);
}
#endif

struct local_rdo_context
{
	const astcenc_contexti* ctx;
	uint32_t base_offset;
};

static bool is_transparent(int v) { return (v & 0xFF) != 0xFF; }

static bool has_any_transparency(
	astcenc_profile decode_mode,
	const symbolic_compressed_block& scb
) {
	if (scb.block_type != SYM_BTYPE_NONCONST) return is_transparent(scb.constant_color[3]);

	vint4 ep0;
	vint4 ep1;
	bool rgb_lns;
	bool a_lns;

	for (int i = 0; i < scb.partition_count; i++)
	{
		unpack_color_endpoints(decode_mode, scb.color_formats[i], scb.color_values[i], rgb_lns, a_lns, ep0, ep1);
		if (is_transparent(ep0.lane<3>()) || is_transparent(ep1.lane<3>())) return true;
	}
	return false;
}

static float compute_block_difference(
	void* user_data,
	const uint8_t* pcb,
	uint32_t local_block_idx,
	float* out_max_std_dev
) {
	const local_rdo_context& local_ctx = *(local_rdo_context*)user_data;
	const astcenc_contexti& ctx = *local_ctx.ctx;

	symbolic_compressed_block scb;
	physical_to_symbolic(*ctx.bsd, pcb, scb);

	// Trial blocks may not be valid at all
	if (scb.block_type == SYM_BTYPE_ERROR) return -ERROR_CALC_DEFAULT;
	bool is_dual_plane = scb.block_type == SYM_BTYPE_NONCONST && ctx.bsd->get_block_mode(scb.block_mode).is_dual_plane;
	if (is_dual_plane && scb.partition_count != 1) return -ERROR_CALC_DEFAULT;
	if (ctx.config.cw_a_weight < 0.01f && has_any_transparency(ctx.config.profile, scb)) return -ERROR_CALC_DEFAULT;

	const astcenc_rdo_context& rdo_ctx = *ctx.rdo_context;
	uint32_t block_idx = local_block_idx + local_ctx.base_offset;
	const image_block& blk = rdo_ctx.m_blocks[block_idx];

	if (out_max_std_dev)
	{
		// ERT expects texel values to be in [0, 255]
		*out_max_std_dev = compute_block_std_dev(blk, 255.0f / 65535.0f);
	}

#if ASTCENC_RDO_SPECIALIZE_DIFF
	float squared_error = 0.0f;

	if (scb.block_type != SYM_BTYPE_NONCONST)
		squared_error = compute_symbolic_block_difference_constant(ctx.config, *ctx.bsd, scb, blk);
	else if (is_dual_plane)
		squared_error = compute_symbolic_block_difference_2plane(ctx.config, *ctx.bsd, scb, blk);
	else if (scb.partition_count == 1)
		squared_error = compute_symbolic_block_difference_1plane_1partition(ctx.config, *ctx.bsd, scb, blk);
	else
		squared_error = compute_symbolic_block_difference_1plane(ctx.config, *ctx.bsd, scb, blk);

	// ERT expects texel values to be in [0, 255]
	return squared_error / blk.texel_count * sqr(255.0f / 65535.0f);
#else
	image_block decoded_blk;
	decoded_blk.decode_unorm8 = blk.decode_unorm8;
	decoded_blk.texel_count = blk.texel_count;
	decoded_blk.channel_weight = blk.channel_weight;

	decompress_symbolic_block(ctx.config.profile, *ctx.bsd, blk.xpos, blk.ypos, blk.zpos, scb, decoded_blk);

	// ERT expects texel values to be in [0, 255]
	return compute_block_mse(blk, decoded_blk, *ctx.bsd, rdo_ctx.m_image_x, rdo_ctx.m_image_y, rdo_ctx.m_image_z, 255.0f / 65535.0f, 255.0f);
#endif
}

void rate_distortion_optimize(
	astcenc_context& ctxo,
	const astcenc_image& image,
	const astcenc_swizzle& swizzle,
	uint8_t* buffer
) {
	if (!ctxo.context.config.rdo_enabled)
	{
		return;
	}

	// Only the first thread actually runs the initializer
	ctxo.manage_rdo.init([&ctxo, &image, &swizzle]
		{
			return init_rdo_context(ctxo.context, image, swizzle);
		},
		ctxo.context.config.progress_callback ? rdo_progress_emitter : nullptr);

	const astcenc_contexti& ctx = ctxo.context;
	uint32_t xblocks = (image.dim_x + ctx.bsd->xdim - 1u) / ctx.bsd->xdim;
	uint32_t yblocks = (image.dim_y + ctx.bsd->ydim - 1u) / ctx.bsd->ydim;
	uint32_t zblocks = (image.dim_z + ctx.bsd->zdim - 1u) / ctx.bsd->zdim;
	uint32_t total_blocks = xblocks * yblocks * zblocks;

	uint32_t blocks_per_task = astc::min(ctx.config.rdo_lookback, total_blocks);
	// There is no way to losslessly partition the job (sequentially dependent on previous output)
	// So we reserve up to one task for each thread to minimize the quality impact.
	uint32_t partitions = ctx.config.rdo_partitions ? ctx.config.rdo_partitions : ctx.thread_count;
	blocks_per_task = astc::max(blocks_per_task, (total_blocks - 1) / partitions + 1);

	uint32_t total_modified = 0;
	while (true)
	{
		uint32_t count;
		uint32_t base = ctxo.manage_rdo.get_task_assignment(blocks_per_task, count);
		if (!count)
		{
			break;
		}

		local_rdo_context local_ctx{ &ctx, base };

		ert::reduce_entropy(buffer + base * ASTCENC_BYTES_PER_BLOCK, count,
							ASTCENC_BYTES_PER_BLOCK, ASTCENC_BYTES_PER_BLOCK,
							ctx.rdo_context->m_ert_params, total_modified,
							&compute_block_difference, &local_ctx);

		ctxo.manage_rdo.complete_task_assignment(count);
	}

	// Wait for rdo to complete before freeing memory
	ctxo.manage_rdo.wait();

	// Only the first thread to arrive actually runs the term
	ctxo.manage_rdo.term([&ctxo]
		{
			delete ctxo.context.rdo_context;
			ctxo.context.rdo_context = nullptr;
		});
}

#endif
