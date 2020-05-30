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

	return ASTCENC_SUCCESS;
}

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

struct encode_astc_image_info
{
	const block_size_descriptor* bsd;
	const error_weighting_params* ewp;
	uint8_t* buffer;
	int pack_and_unpack;
	int threadcount;
	astcenc_profile decode_mode;
	swizzlepattern swz_encode;
	swizzlepattern swz_decode;
	const astc_codec_image* input_image;
	astc_codec_image* output_image;
};

/*static*/ void encode_astc_image_threadfunc(
	int thread_count,
	int thread_id,
	void* vblk
) {
	const encode_astc_image_info *blk = (const encode_astc_image_info *)vblk;
	const block_size_descriptor *bsd = blk->bsd;
	int xdim = bsd->xdim;
	int ydim = bsd->ydim;
	int zdim = bsd->zdim;
	uint8_t *buffer = blk->buffer;
	const error_weighting_params *ewp = blk->ewp;
	int pack_and_unpack = blk->pack_and_unpack;
	astcenc_profile decode_mode = blk->decode_mode;
	swizzlepattern swz_encode = blk->swz_encode;
	swizzlepattern swz_decode = blk->swz_decode;
	const astc_codec_image *input_image = blk->input_image;
	astc_codec_image *output_image = blk->output_image;

	imageblock pb;
	int ctr = thread_id;
	int pctr = 0;

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
					fetch_imageblock(input_image, &pb, bsd, x * xdim, y * ydim, z * zdim, swz_encode);
					symbolic_compressed_block scb;
					compress_symbolic_block(input_image, decode_mode, bsd, ewp, &pb, &scb, &temp_buffers);
					if (pack_and_unpack)
					{
						decompress_symbolic_block(input_image, decode_mode, bsd, x * xdim, y * ydim, z * zdim, &scb, &pb);
						write_imageblock(output_image, &pb, bsd, x * xdim, y * ydim, z * zdim, swz_decode);
					}
					else
					{
						physical_compressed_block pcb;
						pcb = symbolic_to_physical(bsd, &scb);
						*(physical_compressed_block *) bp = pcb;
					}

					ctr = thread_count - 1;
					pctr++;
				}
				else
					ctr--;
			}
		}
	}

	// TODO - move this to the context creation
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

astcenc_error astcenc_compress_image(
	astcenc_context* context,
	astcenc_image const* image,
	size_t image_len,
	uint8_t* data_out,
	size_t data_len,
	int thread_index
) {
	(void)context;
	(void)image;
	(void)image_len;
	(void)data_out;
	(void)data_len;
	(void)thread_index;
	return ASTCENC_SUCCESS;
#if 0
	encode_astc_image_info ai;
	ai.bsd = contex->bsd;
	ai.buffer = data_out;
	ai.ewp = ewp;
	ai.pack_and_unpack = 0;
	ai.decode_mode = context.;
	ai.swz_encode = swz_encode;
	ai.swz_decode = swz_decode;
	ai.input_image = input_image;
	ai.output_image = output_image;

	launch_threads(threadcount, encode_astc_image_threadfunc, &ai);
#endif
}
