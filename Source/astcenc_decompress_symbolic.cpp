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
 * @brief Functions to decompress a symbolic block.
 */

#include "astcenc_internal.h"

#include <stdio.h>

static int compute_value_of_texel_int(
	int texel_to_get,
	const decimation_table* it,
	const int* weights
) {
	int summed_value = 8;
	int weights_to_evaluate = it->texel_num_weights[texel_to_get];
	for (int i = 0; i < weights_to_evaluate; i++)
	{
		summed_value += weights[it->texel_weights[texel_to_get][i]] * it->texel_weights_int[texel_to_get][i];
	}
	return summed_value >> 4;
}

static uint4 lerp_color_int(
	astcenc_profile decode_mode,
	uint4 color0,
	uint4 color1,
	int weight,
	int plane2_weight,
	int plane2_color_component	// -1 in 1-plane mode
) {
	uint4 weight1 = uint4(
		plane2_color_component == 0 ? plane2_weight : weight,
		plane2_color_component == 1 ? plane2_weight : weight,
		plane2_color_component == 2 ? plane2_weight : weight,
		plane2_color_component == 3 ? plane2_weight : weight);

	uint4 weight0 = uint4(64, 64, 64, 64) - weight1;

	if (decode_mode == ASTCENC_PRF_LDR_SRGB)
	{
		color0 = uint4(color0.x >> 8, color0.y >> 8, color0.z >> 8, color0.w >> 8);
		color1 = uint4(color1.x >> 8, color1.y >> 8, color1.z >> 8, color1.w >> 8);
	}

	uint4 color = (color0 * weight0) + (color1 * weight1) + uint4(32, 32, 32, 32);
	color = uint4(color.x >> 6, color.y >> 6, color.z >> 6, color.w >> 6);

	if (decode_mode == ASTCENC_PRF_LDR_SRGB)
	{
		color = color * 257;
	}

	return color;
}

void decompress_symbolic_block(
	astcenc_profile decode_mode,
	const block_size_descriptor* bsd,
	int xpos,
	int ypos,
	int zpos,
	const symbolic_compressed_block* scb,
	imageblock* blk
) {
	blk->xpos = xpos;
	blk->ypos = ypos;
	blk->zpos = zpos;

	// if we detected an error-block, blow up immediately.
	if (scb->error_block)
	{
		if (decode_mode == ASTCENC_PRF_LDR_SRGB)
		{
			for (int i = 0; i < bsd->texel_count; i++)
			{
				blk->orig_data[4 * i] = 1.0f;
				blk->orig_data[4 * i + 1] = 0.0f;
				blk->orig_data[4 * i + 2] = 1.0f;
				blk->orig_data[4 * i + 3] = 1.0f;
				blk->rgb_lns[i] = 0;
				blk->alpha_lns[i] = 0;
				blk->nan_texel[i] = 0;
			}
		}
		else
		{
			for (int i = 0; i < bsd->texel_count; i++)
			{
				blk->orig_data[4 * i] = 0.0f;
				blk->orig_data[4 * i + 1] = 0.0f;
				blk->orig_data[4 * i + 2] = 0.0f;
				blk->orig_data[4 * i + 3] = 0.0f;
				blk->rgb_lns[i] = 0;
				blk->alpha_lns[i] = 0;
				blk->nan_texel[i] = 1;
			}
		}

		imageblock_initialize_work_from_orig(blk, bsd->texel_count);
		update_imageblock_flags(blk, bsd->xdim, bsd->ydim, bsd->zdim);
		return;
	}

	if (scb->block_mode < 0)
	{
		float red = 0, green = 0, blue = 0, alpha = 0;
		int use_lns = 0;
		int use_nan = 0;

		if (scb->block_mode == -2)
		{
			int ired = scb->constant_color[0];
			int igreen = scb->constant_color[1];
			int iblue = scb->constant_color[2];
			int ialpha = scb->constant_color[3];

			// For sRGB decoding a real decoder would just use the top 8 bits
			// for color conversion. We don't color convert, so linearly scale
			// the top 8 bits into the full 16 bit dynamic range
			if (decode_mode == ASTCENC_PRF_LDR_SRGB)
			{
				ired = (ired >> 8) * 257;
				igreen = (igreen >> 8) * 257;
				iblue = (iblue >> 8) * 257;
				ialpha = (ialpha >> 8) * 257;
			}

			red = sf16_to_float(unorm16_to_sf16(ired));
			green = sf16_to_float(unorm16_to_sf16(igreen));
			blue = sf16_to_float(unorm16_to_sf16(iblue));
			alpha = sf16_to_float(unorm16_to_sf16(ialpha));
			use_lns = 0;
			use_nan = 0;
		}
		else
		{
			switch (decode_mode)
			{
			case ASTCENC_PRF_LDR_SRGB:
				red = 1.0f;
				green = 0.0f;
				blue = 1.0f;
				alpha = 1.0f;
				use_lns = 0;
				use_nan = 0;
				break;
			case ASTCENC_PRF_LDR:
				red = 0.0f;
				green = 0.0f;
				blue = 0.0f;
				alpha = 0.0f;
				use_lns = 0;
				use_nan = 1;
				break;
			case ASTCENC_PRF_HDR_RGB_LDR_A:
			case ASTCENC_PRF_HDR:
				// constant-color block; unpack from FP16 to FP32.
				red = sf16_to_float(scb->constant_color[0]);
				green = sf16_to_float(scb->constant_color[1]);
				blue = sf16_to_float(scb->constant_color[2]);
				alpha = sf16_to_float(scb->constant_color[3]);
				use_lns = 1;
				use_nan = 0;
				break;
			}
		}

		for (int i = 0; i < bsd->texel_count; i++)
		{
			blk->orig_data[4 * i] = red;
			blk->orig_data[4 * i + 1] = green;
			blk->orig_data[4 * i + 2] = blue;
			blk->orig_data[4 * i + 3] = alpha;
			blk->rgb_lns[i] = use_lns;
			blk->alpha_lns[i] = use_lns;
			blk->nan_texel[i] = use_nan;
		}

		imageblock_initialize_work_from_orig(blk, bsd->texel_count);
		update_imageblock_flags(blk, bsd->xdim, bsd->ydim, bsd->zdim);
		return;
	}

	// get the appropriate partition-table entry
	int partition_count = scb->partition_count;
	const partition_info *pt = get_partition_table(bsd, partition_count);
	pt += scb->partition_index;

	// get the appropriate block descriptor
	const decimation_table *const *ixtab2 = bsd->decimation_tables;

	const decimation_table *it = ixtab2[bsd->block_modes[scb->block_mode].decimation_mode];

	int is_dual_plane = bsd->block_modes[scb->block_mode].is_dual_plane;

	int weight_quantization_level = bsd->block_modes[scb->block_mode].quantization_mode;

	// decode the color endpoints
	uint4 color_endpoint0[4];
	uint4 color_endpoint1[4];
	int rgb_hdr_endpoint[4];
	int alpha_hdr_endpoint[4];
	int nan_endpoint[4];

	for (int i = 0; i < partition_count; i++)
	{
		unpack_color_endpoints(decode_mode,
		                       scb->color_formats[i],
		                       scb->color_quantization_level,
		                       scb->color_values[i],
		                       &(rgb_hdr_endpoint[i]),
		                       &(alpha_hdr_endpoint[i]),
		                       &(nan_endpoint[i]),
		                       &(color_endpoint0[i]),
		                       &(color_endpoint1[i]));
	}

	// first unquantize the weights
	int uq_plane1_weights[MAX_WEIGHTS_PER_BLOCK];
	int uq_plane2_weights[MAX_WEIGHTS_PER_BLOCK];
	int weight_count = it->num_weights;

	const quantization_and_transfer_table *qat = &(quant_and_xfer_tables[weight_quantization_level]);

	for (int i = 0; i < weight_count; i++)
	{
		uq_plane1_weights[i] = qat->unquantized_value[scb->plane1_weights[i]];
	}

	if (is_dual_plane)
	{
		for (int i = 0; i < weight_count; i++)
		{
			uq_plane2_weights[i] = qat->unquantized_value[scb->plane2_weights[i]];
		}
	}

	// then undecimate them.
	int weights[MAX_TEXELS_PER_BLOCK];
	int plane2_weights[MAX_TEXELS_PER_BLOCK];

	for (int i = 0; i < bsd->texel_count; i++)
	{
		weights[i] = compute_value_of_texel_int(i, it, uq_plane1_weights);
	}

	if (is_dual_plane)
	{
		for (int i = 0; i < bsd->texel_count; i++)
		{
			plane2_weights[i] = compute_value_of_texel_int(i, it, uq_plane2_weights);
		}
	}

	int plane2_color_component = scb->plane2_color_component;

	// now that we have endpoint colors and weights, we can unpack actual colors for
	// each texel.
	for (int i = 0; i < bsd->texel_count; i++)
	{
		int partition = pt->partition_of_texel[i];

		uint4 color = lerp_color_int(decode_mode,
		                             color_endpoint0[partition],
		                             color_endpoint1[partition],
		                             weights[i],
		                             plane2_weights[i],
		                             is_dual_plane ? plane2_color_component : -1);

		blk->rgb_lns[i] = rgb_hdr_endpoint[partition];
		blk->alpha_lns[i] = alpha_hdr_endpoint[partition];
		blk->nan_texel[i] = nan_endpoint[partition];

		blk->data_r[i] = (float)color.x;
		blk->data_g[i] = (float)color.y;
		blk->data_b[i] = (float)color.z;
		blk->data_a[i] = (float)color.w;
	}

	imageblock_initialize_orig_from_work(blk, bsd->texel_count);
	update_imageblock_flags(blk, bsd->xdim, bsd->ydim, bsd->zdim);
}

float compute_imageblock_difference(
	const block_size_descriptor* bsd,
	const imageblock* p1,
	const imageblock* p2,
	const error_weight_block* ewb
) {
	int texels_per_block = bsd->texel_count;
	float summa = 0.0f;

	for (int i = 0; i < texels_per_block; i++)
	{
		float rdiff = fabsf(p1->data_r[i] - p2->data_r[i]);
		float gdiff = fabsf(p1->data_g[i] - p2->data_g[i]);
		float bdiff = fabsf(p1->data_b[i] - p2->data_b[i]);
		float adiff = fabsf(p1->data_a[i] - p2->data_a[i]);
		rdiff = MIN(rdiff, 1e15f);
		gdiff = MIN(gdiff, 1e15f);
		bdiff = MIN(bdiff, 1e15f);
		adiff = MIN(adiff, 1e15f);

		summa += rdiff * rdiff * ewb->error_weights[i].x +
		         gdiff * gdiff * ewb->error_weights[i].y +
		         bdiff * bdiff * ewb->error_weights[i].z +
		         adiff * adiff * ewb->error_weights[i].w;
	}

	return summa;
}
