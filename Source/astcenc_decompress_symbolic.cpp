// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2011-2021 Arm Limited
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
#include <assert.h>

static int compute_value_of_texel_int(
	int texel_to_get,
	const decimation_table* dt,
	const int* weights
) {
	int summed_value = 8;
	int weights_to_evaluate = dt->texel_weight_count[texel_to_get];
	for (int i = 0; i < weights_to_evaluate; i++)
	{
		summed_value += weights[dt->texel_weights_t4[texel_to_get][i]]
		              * dt->texel_weights_int_t4[texel_to_get][i];
	}
	return summed_value >> 4;
}

static vint4 lerp_color_int(
	astcenc_profile decode_mode,
	vint4 color0,
	vint4 color1,
	int weight,
	int plane2_weight,
	vmask4 plane2_mask
) {
	vint4 weight1 = select(vint4(weight), vint4(plane2_weight), plane2_mask);
	vint4 weight0 = vint4(64) - weight1;

	if (decode_mode == ASTCENC_PRF_LDR_SRGB)
	{
		color0 = asr<8>(color0);
		color1 = asr<8>(color1);
	}

	vint4 color = (color0 * weight0) + (color1 * weight1) + vint4(32);
	color = asr<6>(color);

	if (decode_mode == ASTCENC_PRF_LDR_SRGB)
	{
		color = color * vint4(257);
	}

	return color;
}

// Turn packed unorm16 or LNS data into generic float data
static inline vfloat4 decode_texel(
	vint4 data,
	vmask4 lns_mask
) {
	vint4 color_lns = vint4::zero();
	vint4 color_unorm = vint4::zero();

	if (any(lns_mask))
	{
		color_lns = lns_to_sf16(data);
	}

	if (!all(lns_mask))
	{
		color_unorm = unorm16_to_sf16(data);
	}

	// Pick channels and then covert to FP16
	vint4 datai = select(color_unorm, color_lns, lns_mask);
	return float16_to_float(datai);
}

void unpack_weights(
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& scb,
	const decimation_table& dt,
	bool is_dual_plane,
	int weight_quant_level,
	int weights_plane1[MAX_TEXELS_PER_BLOCK],
	int weights_plane2[MAX_TEXELS_PER_BLOCK]
) {
	// First, unquantize the weights ...
	int uq_plane1_weights[MAX_WEIGHTS_PER_BLOCK];
	int uq_plane2_weights[MAX_WEIGHTS_PER_BLOCK];
	int weight_count = dt.weight_count;

	const quantization_and_transfer_table *qat = &(quant_and_xfer_tables[weight_quant_level]);

	for (int i = 0; i < weight_count; i++)
	{
		uq_plane1_weights[i] = qat->unquantized_value[scb.weights[i]];
	}

	if (is_dual_plane)
	{
		for (int i = 0; i < weight_count; i++)
		{
			uq_plane2_weights[i] = qat->unquantized_value[scb.weights[i + PLANE2_WEIGHTS_OFFSET]];
		}
	}

	// Second, undecimate the weights ...
	for (int i = 0; i < bsd.texel_count; i++)
	{
		weights_plane1[i] = compute_value_of_texel_int(i, &dt, uq_plane1_weights);
	}

	if (is_dual_plane)
	{
		for (int i = 0; i < bsd.texel_count; i++)
		{
			weights_plane2[i] = compute_value_of_texel_int(i, &dt, uq_plane2_weights);
		}
	}
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

	blk->data_min = vfloat4::zero();
	blk->data_max = vfloat4::zero();
	blk->grayscale = false;

	// if we detected an error-block, blow up immediately.
	if (scb->error_block)
	{
		// TODO: Check this - isn't linear LDR magenta too? Same below ...
		if (decode_mode == ASTCENC_PRF_LDR_SRGB)
		{
			for (int i = 0; i < bsd->texel_count; i++)
			{
				blk->data_r[i] = 1.0f;
				blk->data_g[i] = 0.0f;
				blk->data_b[i] = 1.0f;
				blk->data_a[i] = 1.0f;
				blk->rgb_lns[i] = 0;
				blk->alpha_lns[i] = 0;
				blk->nan_texel[i] = 0;
			}
		}
		else
		{
			for (int i = 0; i < bsd->texel_count; i++)
			{
				blk->data_r[i] = 0.0f;
				blk->data_g[i] = 0.0f;
				blk->data_b[i] = 0.0f;
				blk->data_a[i] = 0.0f;
				blk->rgb_lns[i] = 0;
				blk->alpha_lns[i] = 0;
				blk->nan_texel[i] = 1;
			}
		}

		return;
	}

	if (scb->block_mode < 0)
	{
		vfloat4 color;
		int use_lns = 0;
		int use_nan = 0;

		if (scb->block_mode == -2)
		{
			vint4 colori(scb->constant_color);

			// For sRGB decoding a real decoder would just use the top 8 bits
			// for color conversion. We don't color convert, so linearly scale
			// the top 8 bits into the full 16 bit dynamic range
			if (decode_mode == ASTCENC_PRF_LDR_SRGB)
			{
				colori = asr<8>(colori) * 257;
			}

			vint4 colorf16 = unorm16_to_sf16(colori);
			color = float16_to_float(colorf16);
		}
		else
		{
			switch (decode_mode)
			{
			case ASTCENC_PRF_LDR_SRGB:
				color = vfloat4(1.0f, 0.0f, 1.0f, 1.0f);
				break;
			case ASTCENC_PRF_LDR:
				color = vfloat4(0.0f);
				use_nan = 1;
				break;
			case ASTCENC_PRF_HDR_RGB_LDR_A:
			case ASTCENC_PRF_HDR:
				// constant-color block; unpack from FP16 to FP32.
				color = float16_to_float(vint4(scb->constant_color));
				use_lns = 1;
				break;
			}
		}

		// TODO: Skip this and add constant color transfer to img block?
		for (int i = 0; i < bsd->texel_count; i++)
		{
			blk->data_r[i] = color.lane<0>();
			blk->data_g[i] = color.lane<1>();
			blk->data_b[i] = color.lane<2>();
			blk->data_a[i] = color.lane<3>();
			blk->rgb_lns[i] = use_lns;
			blk->alpha_lns[i] = use_lns;
			blk->nan_texel[i] = use_nan;
		}

		return;
	}

	// get the appropriate partition-table entry
	int partition_count = scb->partition_count;
	const partition_info *pt = get_partition_table(bsd, partition_count);
	pt += scb->partition_index;

	// get the appropriate block descriptor
	const decimation_table *const *dts = bsd->decimation_tables;

	const int packed_index = bsd->block_mode_packed_index[scb->block_mode];
	assert(packed_index >= 0 && packed_index < bsd->block_mode_count);
	const block_mode& bm = bsd->block_modes[packed_index];
	const decimation_table *dt = dts[bm.decimation_mode];

	int is_dual_plane = bm.is_dual_plane;

	int weight_quant_level = bm.quant_mode;

	// decode the color endpoints
	vint4 color_endpoint0[4];
	vint4 color_endpoint1[4];
	int rgb_hdr_endpoint[4];
	int alpha_hdr_endpoint[4];
	int nan_endpoint[4];

	for (int i = 0; i < partition_count; i++)
	{
		unpack_color_endpoints(decode_mode,
		                       scb->color_formats[i],
		                       scb->color_quant_level,
		                       scb->color_values[i],
		                       &(rgb_hdr_endpoint[i]),
		                       &(alpha_hdr_endpoint[i]),
		                       &(nan_endpoint[i]),
		                       &(color_endpoint0[i]),
		                       &(color_endpoint1[i]));
	}

	// Unquantize and undecimate the weights
	int weights[MAX_TEXELS_PER_BLOCK];
	int plane2_weights[MAX_TEXELS_PER_BLOCK];
	unpack_weights(*bsd, *scb, *dt, is_dual_plane, weight_quant_level, weights, plane2_weights);

	// Now that we have endpoint colors and weights, we can unpack texel colors
	int plane2_color_component = is_dual_plane ? scb->plane2_color_component : -1;
	vmask4 plane2_mask = vint4::lane_id() == vint4(plane2_color_component);

	for (int i = 0; i < partition_count; i++)
	{
		vint4 ep0 = color_endpoint0[i];
		vint4 ep1 = color_endpoint1[i];
		bool rgb_lns = rgb_hdr_endpoint[i];
		bool nan = nan_endpoint[i];
		bool a_lns = alpha_hdr_endpoint[i];
		vmask4 lns_mask(rgb_lns, rgb_lns, rgb_lns, a_lns);

		int texel_count = pt->partition_texel_count[i];
		for (int j = 0; j < texel_count; j++)
		{
			int tix = pt->texels_of_partition[i][j];
			vint4 color = lerp_color_int(decode_mode,
			                             ep0,
			                             ep1,
			                             weights[tix],
			                             plane2_weights[tix],
			                             plane2_mask);

			vfloat4 colorf = decode_texel(color, lns_mask);

			blk->nan_texel[tix] = nan;
			blk->data_r[tix] = colorf.lane<0>();
			blk->data_g[tix] = colorf.lane<1>();
			blk->data_b[tix] = colorf.lane<2>();
			blk->data_a[tix] = colorf.lane<3>();
		}
	}
}

// Returns a negative error for encodings we want to reject as a part of a
// heuristic check, e.g. for RGBM textures which have zero M values.
float compute_symbolic_block_difference(
	const astcenc_config& config,
	const block_size_descriptor* bsd,
	const symbolic_compressed_block* scb,
	const imageblock* blk,
	const error_weight_block *ewb
) {
	// if we detected an error-block, blow up immediately.
	if (scb->error_block)
	{
		return 1e29f;
	}

	assert(scb->block_mode >= 0);

	// get the appropriate partition-table entry
	int partition_count = scb->partition_count;
	const partition_info *pt = get_partition_table(bsd, partition_count);
	pt += scb->partition_index;

	// get the appropriate block descriptor
	const decimation_table *const *dts = bsd->decimation_tables;

	const int packed_index = bsd->block_mode_packed_index[scb->block_mode];
	assert(packed_index >= 0 && packed_index < bsd->block_mode_count);
	const block_mode& bm = bsd->block_modes[packed_index];
	const decimation_table *dt = dts[bm.decimation_mode];

	int is_dual_plane = bm.is_dual_plane;
	int weight_quant_level = bm.quant_mode;

	int weight_count = dt->weight_count;
	int texel_count = bsd->texel_count;

	promise(partition_count > 0);
	promise(weight_count > 0);
	promise(texel_count > 0);

	// decode the color endpoints
	vint4 color_endpoint0[4];
	vint4 color_endpoint1[4];
	int rgb_hdr_endpoint[4];
	int alpha_hdr_endpoint[4];
	int nan_endpoint[4];

	for (int i = 0; i < partition_count; i++)
	{
		unpack_color_endpoints(config.profile,
		                       scb->color_formats[i],
		                       scb->color_quant_level,
		                       scb->color_values[i],
		                       &(rgb_hdr_endpoint[i]),
		                       &(alpha_hdr_endpoint[i]),
		                       &(nan_endpoint[i]),
		                       &(color_endpoint0[i]),
		                       &(color_endpoint1[i]));
	}

	// Unquantize and undecimate the weights
	int weights[MAX_TEXELS_PER_BLOCK];
	int plane2_weights[MAX_TEXELS_PER_BLOCK];
	unpack_weights(*bsd, *scb, *dt, is_dual_plane, weight_quant_level, weights, plane2_weights);

	// Now that we have endpoint colors and weights, we can unpack texel colors
	int plane2_color_component = is_dual_plane ? scb->plane2_color_component : -1;
	vmask4 plane2_mask = vint4::lane_id() == vint4(plane2_color_component);

	float summa = 0.0f;
	for (int i = 0; i < texel_count; i++)
	{
		int partition = pt->partition_of_texel[i];

		vint4 ep0 = color_endpoint0[partition];
		vint4 ep1 = color_endpoint1[partition];

		vint4 colori = lerp_color_int(config.profile,
		                              ep0,
		                              ep1,
		                              weights[i],
		                              plane2_weights[i],
		                              plane2_mask);

		vfloat4 color = int_to_float(colori);
		vfloat4 oldColor = blk->texel(i);

		if (config.flags & ASTCENC_FLG_MAP_RGBM)
		{
			// Fail encodings that result in zero weight M pixels. Note that
			// this can cause "interesting" artifacts if we reject all useful
			// encodings - we typically get max brightness encodings instead
			// which look just as bad. We recommend users apply a bias to their
			// stored M value, limiting the lower value to 16 or 32 to avoid
			// getting small M values post-quantization, but we can't prove it
			// would never happen, especially at low bit rates ...
			if (color.lane<3>() == 0.0f)
			{
				return -1e30f;
			}

			// Compute error based on decoded RGBM color
			color = vfloat4(
				color.lane<0>() * color.lane<3>() * config.rgbm_m_scale,
				color.lane<1>() * color.lane<3>() * config.rgbm_m_scale,
				color.lane<2>() * color.lane<3>() * config.rgbm_m_scale,
				1.0f
			);

			oldColor = vfloat4(
				oldColor.lane<0>() * oldColor.lane<3>() * config.rgbm_m_scale,
				oldColor.lane<1>() * oldColor.lane<3>() * config.rgbm_m_scale,
				oldColor.lane<2>() * oldColor.lane<3>() * config.rgbm_m_scale,
				1.0f
			);
		}

		vfloat4 error = oldColor - color;
		error = min(abs(error), 1e15f);
		error = error * error;

		float metric = dot_s(error, ewb->error_weights[i]);
		summa += astc::min(metric, 1e30f);
	}

	return summa;
}
