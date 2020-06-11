// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2011-2020 Arm Limited
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
 * @brief Functions for the library entrypoint.
 */

#include <cstring>
#include <new>

#include "astcenc.h"
#include "astcenc_internal.h"

__attribute__((visibility("default")))
astcenc_error astcenc_init_config(
	astcenc_profile profile,
	unsigned int block_x,
	unsigned int block_y,
	unsigned int block_z,
	astcenc_preset preset,
	unsigned int flags,
	astcenc_config& config
) {
	// Zero init all config fields; although most of will be over written
	std::memset(&config, 0, sizeof(config));

	// Process the block size
	if (block_z <= 1) {
		if (!is_legal_2d_block_size(block_x, block_y))
		{
			return ASTCENC_ERR_BAD_BLOCK_SIZE;
		}
		block_z = 1;
	}
	else if (!is_legal_3d_block_size(block_x, block_y, block_z))
	{
		return ASTCENC_ERR_BAD_BLOCK_SIZE;
	}

	config.block_x = block_x;
	config.block_y = block_y;
	config.block_z = block_z;

	float texels = block_x * block_y * block_z;
	float ltexels = logf(texels) / logf(10.0f);

	// Process the performance preset; note that this must be done before we
	// process any additional settings, such as color profile and flags, which
	// may replace some of these settings with more use case tuned values
	switch(preset) {
		case ASTCENC_PRE_FAST:
			config.tune_partition_limit = 4;
			config.tune_block_mode_limit = 50;
			config.tune_refinement_limit = 1;
			config.tune_db_limit = MAX(85 - 35 * ltexels, 63 - 19 * ltexels);
			config.tune_partition_early_out_limit = 1.0f;
			config.tune_two_plane_early_out_limit = 0.5f;
			break;
		case ASTCENC_PRE_MEDIUM:
			config.tune_partition_limit = 25;
			config.tune_block_mode_limit = 75;
			config.tune_refinement_limit = 2;
			config.tune_db_limit = MAX(95 - 35 * ltexels, 70 - 19 * ltexels);
			config.tune_partition_early_out_limit = 1.2f;
			config.tune_two_plane_early_out_limit = 0.75f;
			break;
		case ASTCENC_PRE_THOROUGH:
			config.tune_partition_limit = 100;
			config.tune_block_mode_limit = 95;
			config.tune_refinement_limit = 4;
			config.tune_db_limit = MAX(105 - 35 * ltexels, 77 - 19 * ltexels);
			config.tune_partition_early_out_limit = 2.5f;
			config.tune_two_plane_early_out_limit = 0.95f;
			break;
		case ASTCENC_PRE_EXHAUSTIVE:
			config.tune_partition_limit = 1024;
			config.tune_block_mode_limit = 100;
			config.tune_refinement_limit = 4;
			config.tune_db_limit = 999.0f;
			config.tune_partition_early_out_limit = 1000.0f;
			config.tune_two_plane_early_out_limit = 0.99f;
			break;
		default:
			return ASTCENC_ERR_BAD_PRESET;
	}

	// Set heuristics to the defaults for each color profile
	config.v_rgba_radius = 0;
	config.v_rgba_mean_stdev_mix = 0.0f;

	config.cw_r_weight = 1.0f;
	config.cw_g_weight = 1.0f;
	config.cw_b_weight = 1.0f;
	config.cw_a_weight = 1.0f;

	config.a_scale_radius = 0;

	config.b_deblock_weight = 0.0f;

	config.profile = profile;
	switch(profile) {
		case ASTCENC_PRF_LDR:
		case ASTCENC_PRF_LDR_SRGB:
			config.v_rgb_power = 1.0f;
			config.v_rgb_base = 1.0f;
			config.v_rgb_mean = 0.0f;
			config.v_rgb_stdev = 0.0f;

			config.v_a_power = 1.0f;
			config.v_a_base = 1.0f;
			config.v_a_mean = 0.0f;
			config.v_a_stdev = 0.0f;
			break;
		case ASTCENC_PRF_HDR_RGB_LDR_A:
			config.v_rgb_power = 0.75f;
			config.v_rgb_base = 0.0f;
			config.v_rgb_mean = 1.0f;
			config.v_rgb_stdev = 0.0f;

			config.v_a_power = 1.0f;
			config.v_a_base = 0.05f;
			config.v_a_mean = 0.0f;
			config.v_a_stdev = 0.0f;

			config.tune_db_limit = 999.0f;
			break;
		case ASTCENC_PRF_HDR:
			config.v_rgb_power = 0.75f;
			config.v_rgb_base = 0.0f;
			config.v_rgb_mean = 1.0f;
			config.v_rgb_stdev = 0.0f;

			config.v_a_power = 0.75f;
			config.v_a_base = 0.0f;
			config.v_a_mean = 1.0f;
			config.v_a_stdev = 0.0f;

			config.tune_db_limit = 999.0f;
			break;
		default:
			return ASTCENC_ERR_BAD_PROFILE;
	}

	// Flags field must not contain any unknown flag bits
	unsigned int exMask;
	exMask = ~ASTCENC_ALL_FLAGS;
	if (astc::popcount(flags & exMask) != 0)
	{
		return ASTCENC_ERR_BAD_FLAGS;
	}

	// Flags field must only contain at most a single map type
	exMask = ASTCENC_FLG_MAP_MASK | ASTCENC_FLG_MAP_NORMAL;
	if (astc::popcount(flags & exMask) > 1)
	{
		return ASTCENC_ERR_BAD_FLAGS;
	}

	if (flags & ASTCENC_FLG_MAP_NORMAL)
	{
		config.cw_r_weight = 1.0f;
		config.cw_g_weight = 0.0f;
		config.cw_b_weight = 0.0f;
		config.cw_a_weight = 1.0f;
		config.tune_partition_early_out_limit = 1000.0f;
		config.tune_two_plane_early_out_limit = 0.99f;

 		if (flags & ASTCENC_FLG_USE_PERCEPTUAL)
		{
			config.b_deblock_weight = 1.8f;
			config.v_rgba_radius = 3;
			config.v_rgba_mean_stdev_mix = 0.0f;
			config.v_rgb_mean = 0.0f;
			config.v_rgb_stdev = 50.0f;
			config.v_a_mean = 0.0f;
			config.v_a_stdev = 50.0f;
		}
	}

	if (flags & ASTCENC_FLG_MAP_MASK)
	{
		config.v_rgba_radius = 3;
		config.v_rgba_mean_stdev_mix = 0.03f;
		config.v_rgb_mean = 0.0f;
		config.v_rgb_stdev = 25.0f;
		config.v_a_mean = 0.0f;
		config.v_a_stdev = 25.0f;
	}

	// TODO: We don't support user threads yet
	if (flags & ASTCENC_FLG_USE_USER_THREADS)
	{
		return ASTCENC_ERR_NOT_IMPLEMENTED;
	}

	config.flags = flags;

	return ASTCENC_SUCCESS;
}

__attribute__((visibility("default")))
astcenc_error astcenc_context_alloc(
	astcenc_config const& config,
	int thread_count,
	astcenc_context** context
) {
	astcenc_context* ctx { nullptr };
	block_size_descriptor* bsd { nullptr };

	try
	{
		ctx = new astcenc_context;
		ctx->thread_count = thread_count;
		ctx->config = config;

		bsd = new block_size_descriptor;
		init_block_size_descriptor(config.block_x, config.block_y, config.block_z, bsd);
		ctx->bsd = bsd;
	}
	catch(const std::bad_alloc&)
	{
		term_block_size_descriptor(bsd);
		delete bsd;
		delete ctx;
		*context = nullptr;
		return ASTCENC_ERR_OUT_OF_MEM;
	}

	*context = ctx;

	// TODO: Currently static memory; should move to context memory
	prepare_angular_tables();
	build_quantization_mode_table();

	return ASTCENC_SUCCESS;
}

__attribute__((visibility("default")))
void astcenc_context_free(
	astcenc_context* context
) {
	if (context)
	{
		term_block_size_descriptor(context->bsd);
		delete context->bsd;
		delete context;
	}
}

struct compress_astc_image_info
{
	const block_size_descriptor* bsd;
	const error_weighting_params* ewp;
	uint8_t* buffer;
	int threadcount;
	astcenc_profile decode_mode;
	astcenc_swizzle swz_encode;
	const astc_codec_image* input_image;
};

static void encode_astc_image_threadfunc(
	int thread_count,
	int thread_id,
	void* vblk
) {
	const compress_astc_image_info *blk = (const compress_astc_image_info *)vblk;
	const block_size_descriptor *bsd = blk->bsd;
	int xdim = bsd->xdim;
	int ydim = bsd->ydim;
	int zdim = bsd->zdim;
	uint8_t *buffer = blk->buffer;
	const error_weighting_params *ewp = blk->ewp;
	astcenc_profile decode_mode = blk->decode_mode;
	astcenc_swizzle swz_encode = blk->swz_encode;
	const astc_codec_image *input_image = blk->input_image;

	imageblock pb;
	int ctr = thread_id;

	int x, y, z;
	int xsize = input_image->xsize;
	int ysize = input_image->ysize;
	int zsize = input_image->zsize;
	int xblocks = (xsize + xdim - 1) / xdim;
	int yblocks = (ysize + ydim - 1) / ydim;
	int zblocks = (zsize + zdim - 1) / zdim;

	//allocate memory for temporary buffers
	compress_symbolic_block_buffers temp_buffers;
	temp_buffers.ewb = new error_weight_block;
	temp_buffers.ewbo = new error_weight_block_orig;
	temp_buffers.tempblocks = new symbolic_compressed_block[4];
	temp_buffers.temp = new imageblock;
	temp_buffers.planes2 = new compress_fixed_partition_buffers;
	temp_buffers.planes2->ei1 = new endpoints_and_weights;
	temp_buffers.planes2->ei2 = new endpoints_and_weights;
	temp_buffers.planes2->eix1 = new endpoints_and_weights[MAX_DECIMATION_MODES];
	temp_buffers.planes2->eix2 = new endpoints_and_weights[MAX_DECIMATION_MODES];
	temp_buffers.planes2->decimated_quantized_weights = new float[2 * MAX_DECIMATION_MODES * MAX_WEIGHTS_PER_BLOCK];
	temp_buffers.planes2->decimated_weights = new float[2 * MAX_DECIMATION_MODES * MAX_WEIGHTS_PER_BLOCK];
	temp_buffers.planes2->flt_quantized_decimated_quantized_weights = new float[2 * MAX_WEIGHT_MODES * MAX_WEIGHTS_PER_BLOCK];
	temp_buffers.planes2->u8_quantized_decimated_quantized_weights = new uint8_t[2 * MAX_WEIGHT_MODES * MAX_WEIGHTS_PER_BLOCK];
	temp_buffers.plane1 = temp_buffers.planes2;

	for (z = 0; z < zblocks; z++)
	{
		for (y = 0; y < yblocks; y++)
		{
			for (x = 0; x < xblocks; x++)
			{
				if (ctr == 0)
				{
					int offset = ((z * yblocks + y) * xblocks + x) * 16;
					uint8_t *bp = buffer + offset;
					fetch_imageblock(decode_mode, input_image, &pb, bsd, x * xdim, y * ydim, z * zdim, swz_encode);
					symbolic_compressed_block scb;
					compress_symbolic_block(input_image, decode_mode, bsd, ewp, &pb, &scb, &temp_buffers);
					*(physical_compressed_block*) bp = symbolic_to_physical(bsd, &scb);
					ctr = thread_count - 1;
				}
				else
					ctr--;
			}
		}
	}

	delete[] temp_buffers.planes2->decimated_quantized_weights;
	delete[] temp_buffers.planes2->decimated_weights;
	delete[] temp_buffers.planes2->flt_quantized_decimated_quantized_weights;
	delete[] temp_buffers.planes2->u8_quantized_decimated_quantized_weights;
	delete[] temp_buffers.planes2->eix1;
	delete[] temp_buffers.planes2->eix2;
	delete   temp_buffers.planes2->ei1;
	delete   temp_buffers.planes2->ei2;
	delete   temp_buffers.planes2;
	delete[] temp_buffers.tempblocks;
	delete   temp_buffers.temp;
	delete   temp_buffers.ewbo;
	delete   temp_buffers.ewb;
}

__attribute__((visibility("default")))
astcenc_error astcenc_compress_image(
	astcenc_context* context,
	astcenc_image& image,
	astcenc_swizzle swizzle,
	uint8_t* data_out,
	size_t data_len,
	int thread_index
) {
	// TODO: Replace error_weighting_params in the core codec with the config / context structs
	error_weighting_params ewp;
	ewp.rgb_power = context->config.v_rgb_power;
	ewp.rgb_base_weight = context->config.v_rgb_base;
	ewp.rgb_mean_weight = context->config.v_rgb_mean;
	ewp.rgb_stdev_weight = context->config.v_rgb_stdev;
	ewp.alpha_power = context->config.v_a_power;
	ewp.alpha_base_weight = context->config.v_a_base;
	ewp.alpha_mean_weight = context->config.v_a_mean;
	ewp.alpha_stdev_weight = context->config.v_a_stdev;
	ewp.rgb_mean_and_stdev_mixing = context->config.v_rgba_mean_stdev_mix;
	ewp.mean_stdev_radius = context->config.v_rgba_radius;
	ewp.enable_rgb_scale_with_alpha = (context->config.flags & ASTCENC_FLG_USE_ALPHA_WEIGHT) == 0 ? 0 : 1;
	ewp.alpha_radius = context->config.a_scale_radius;
	ewp.ra_normal_angular_scale = (context->config.flags & ASTCENC_FLG_MAP_NORMAL) == 0 ? 0 : 1;
	ewp.block_artifact_suppression = context->config.b_deblock_weight;
	ewp.rgba_weights[0] = MAX(context->config.cw_r_weight, 0.001f);
	ewp.rgba_weights[1] = MAX(context->config.cw_g_weight, 0.001f);
	ewp.rgba_weights[2] = MAX(context->config.cw_b_weight, 0.001f);
	ewp.rgba_weights[3] = MAX(context->config.cw_a_weight, 0.001f);
	ewp.partition_search_limit = context->config.tune_partition_limit;
	ewp.block_mode_cutoff = context->config.tune_block_mode_limit / 100.0f;
	if ((context->config.profile == ASTCENC_PRF_LDR) || (context->config.profile == ASTCENC_PRF_LDR_SRGB))
	{
		ewp.texel_avg_error_limit  = powf(0.1f, context->config.tune_db_limit * 0.1f) * 65535.0f * 65535.0f;
	}
	else
	{
		ewp.texel_avg_error_limit = 0.0f;
	}

	ewp.partition_1_to_2_limit = context->config.tune_partition_early_out_limit;
	ewp.lowest_correlation_cutoff = context->config.tune_two_plane_early_out_limit;
	ewp.max_refinement_iters = context->config.tune_refinement_limit;

	// TODO: Replace astc_codec_image in the core codec with the astcenc_image struct
	astc_codec_image input_image;
	input_image.data8 = image.data8;
	input_image.data16 = image.data16;
	input_image.xsize = image.dim_x;
	input_image.ysize = image.dim_y;
	input_image.zsize = image.dim_z;
	input_image.padding = image.padding_texels;

	// Need to agree what we do with linearize sRGB
	input_image.linearize_srgb = (context->config.flags & ASTCENC_FLG_USE_LINEARIZED_SRGB) == 0 ? 0 : 1;

	input_image.input_averages = nullptr;
	input_image.input_variances = nullptr;
	input_image.input_alpha_averages = nullptr;

	if (image.padding_texels > 0 ||
	    ewp.rgb_mean_weight != 0.0f || ewp.rgb_stdev_weight != 0.0f ||
	    ewp.alpha_mean_weight != 0.0f || ewp.alpha_stdev_weight != 0.0f)
	{
		compute_averages_and_variances(&input_image, ewp.rgb_power, ewp.alpha_power,
		                               ewp.mean_stdev_radius, ewp.alpha_radius,
		                               input_image.linearize_srgb, swizzle, context->thread_count);
	}

	// TODO: This could be done once when the context is created
	expand_block_artifact_suppression(
		context->config.block_x, context->config.block_y, context->config.block_z, &ewp);

	compress_astc_image_info ai;
	ai.bsd = context->bsd;
	ai.buffer = data_out;
	ai.ewp = &ewp;
	ai.decode_mode = context->config.profile;
	ai.swz_encode = swizzle;
	ai.input_image = &input_image;

	// TODO: Add bounds checking
	(void)data_len;

	// TODO: Implement user-thread pools
	(void)thread_index;

	launch_threads(context->thread_count, encode_astc_image_threadfunc, &ai);

	return ASTCENC_SUCCESS;
}


__attribute__((visibility("default")))
astcenc_error astcenc_decompress_image(
	astcenc_context* context,
	uint8_t* data,
	size_t data_len,
	astcenc_image& image_out,
	astcenc_swizzle swizzle,
	int thread_index
) {
	unsigned int block_x = context->config.block_x;
	unsigned int block_y = context->config.block_y;
	unsigned int block_z = context->config.block_z;

	unsigned int xblocks = (image_out.dim_x + block_x - 1) / block_x;
	unsigned int yblocks = (image_out.dim_y + block_y - 1) / block_y;
	unsigned int zblocks = (image_out.dim_z + block_z - 1) / block_z;

	// TODO: Check output bounds
	(void)data_len;

	// TODO: Handle custom threading
	(void)thread_index;

	// TODO: Replace astc_codec_image in the core codec with the astcenc_image struct
	astc_codec_image image;
	image.data8 = image_out.data8;
	image.data16 = image_out.data16;
	image.xsize = image_out.dim_x;
	image.ysize = image_out.dim_y;
	image.zsize = image_out.dim_z;
	image.padding = image_out.padding_texels;

	// Need to agree what we do with linearize sRGB
	image.linearize_srgb = (context->config.flags & ASTCENC_FLG_USE_LINEARIZED_SRGB) == 0 ? 0 : 1;

	image.input_averages = nullptr;
	image.input_variances = nullptr;
	image.input_alpha_averages = nullptr;

	imageblock pb;
	for (unsigned int z = 0; z < zblocks; z++)
	{
		for (unsigned int y = 0; y < yblocks; y++)
		{
			for (unsigned int x = 0; x < xblocks; x++)
			{
				unsigned int offset = (((z * yblocks + y) * xblocks) + x) * 16;
				uint8_t* bp = data + offset;
				physical_compressed_block pcb = *(physical_compressed_block *) bp;
				symbolic_compressed_block scb;
				physical_to_symbolic(context->bsd, pcb, &scb);
				decompress_symbolic_block(context->config.profile, context->bsd,
				                          x * block_x, y * block_y, z * block_z,
				                          &scb, &pb);
				write_imageblock(&image, &pb, context->bsd,
				                 x * block_x, y * block_y, z * block_z, swizzle);
			}
		}
	}

	return ASTCENC_SUCCESS;
}


__attribute__((visibility("default")))
const char* astcenc_get_error_string(
	astcenc_error status
) {
	switch(status) {
		case ASTCENC_ERR_OUT_OF_MEM:
			return "ASTCENC_ERR_OUT_OF_MEM";
		case ASTCENC_ERR_BAD_PARAM:
			return "ASTCENC_ERR_BAD_PARAM";
		case ASTCENC_ERR_BAD_BLOCK_SIZE:
			return "ASTCENC_ERR_BAD_BLOCK_SIZE";
		case ASTCENC_ERR_BAD_PROFILE:
			return "ASTCENC_ERR_BAD_PROFILE";
		case ASTCENC_ERR_BAD_PRESET:
			return "ASTCENC_ERR_BAD_PRESET";
		case ASTCENC_ERR_BAD_FLAGS:
			return "ASTCENC_ERR_BAD_FLAGS";
		case ASTCENC_ERR_NOT_IMPLEMENTED:
			return "ASTCENC_ERR_NOT_IMPLEMENTED";
		default:
			return nullptr;
	}
}
