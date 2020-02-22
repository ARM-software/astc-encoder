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
 * @brief Functions to compress a symbolic block.
 */

#include "astc_codec_internals.h"

#include <string.h>
#include <stdio.h>

#ifdef DEBUG_CAPTURE_NAN
	#ifndef _GNU_SOURCE
		#define _GNU_SOURCE
	#endif

	#include <fenv.h>
#endif

/**
 * @brief Attempt to improve weights given a chosen configuration.
 *
 * Given a fixed weight grid decimation and weight value quantization, iterate
 * over all weights (per partition and per plane) and attempt to improve image
 * quality by moving each weight up by one or down by one quantization step.
 */
static int realign_weights(
	astc_decode_mode decode_mode,
	const block_size_descriptor* bsd,
	const imageblock* blk,
	const error_weight_block* ewb,
	symbolic_compressed_block* scb,
	uint8_t* plane1_weight_set8,
	uint8_t* plane2_weight_set8
) {
	// Get the partition descriptor
	int partition_count = scb->partition_count;
	const partition_info *pt = get_partition_table(bsd, partition_count);
	pt += scb->partition_index;

	// Get the quantization table
	int weight_quantization_level = bsd->block_modes[scb->block_mode].quantization_mode;
	const quantization_and_transfer_table *qat = &(quant_and_xfer_tables[weight_quantization_level]);

	// Get the decimation table
	const decimation_table *const *ixtab2 = bsd->decimation_tables;
	const decimation_table *it = ixtab2[bsd->block_modes[scb->block_mode].decimation_mode];
	int weight_count = it->num_weights;

	int max_plane = bsd->block_modes[scb->block_mode].is_dual_plane;
	int plane2_component = max_plane ? scb->plane2_color_component : 0;
	int plane_mask = max_plane ? 1 << plane2_component : 0;

	// Decode the color endpoints
	int rgb_hdr;
	int alpha_hdr;
	int nan_endpoint;
	int4 endpnt0[4];
	int4 endpnt1[4];
	float4 endpnt0f[4];
	float4 offset[4];

	for (int pa_idx = 0; pa_idx < partition_count; pa_idx++)
	{
		unpack_color_endpoints(decode_mode,
		                       scb->color_formats[pa_idx],
		                       scb->color_quantization_level,
		                       scb->color_values[pa_idx],
		                       &rgb_hdr, &alpha_hdr, &nan_endpoint,
		                       // TODO: Fix these casts ...
		                       reinterpret_cast<uint4*>(&endpnt0[pa_idx]),
		                       reinterpret_cast<uint4*>(&endpnt1[pa_idx]));
	}

	uint8_t uq_pl_weights[MAX_WEIGHTS_PER_BLOCK];
	uint8_t* weight_set8 = plane1_weight_set8;
	int adjustments = 0;

	// For each plane and partition ...
	for (int pl_idx = 0; pl_idx <= max_plane; pl_idx++)
	{
		for (int pa_idx = 0; pa_idx < partition_count; pa_idx++)
		{
			// Compute the endpoint delta for all channels in current plane
			int4 epd = endpnt1[pa_idx] - endpnt0[pa_idx];

			if (plane_mask & 1) epd.x = 0;
			if (plane_mask & 2) epd.y = 0;
			if (plane_mask & 4) epd.z = 0;
			if (plane_mask & 8) epd.w = 0;

			endpnt0f[pa_idx] = float4((float)endpnt0[pa_idx].x, (float)endpnt0[pa_idx].y,
			                          (float)endpnt0[pa_idx].z, (float)endpnt0[pa_idx].w);
			offset[pa_idx] = float4((float)epd.x, (float)epd.y, (float)epd.z, (float)epd.w);
			offset[pa_idx] = offset[pa_idx] * (1.0f / 64.0f);
		}

		// Create an unquantized weight grid for this decimation level
		for (int we_idx = 0; we_idx < weight_count; we_idx++)
		{
			uq_pl_weights[we_idx] = qat->unquantized_value[weight_set8[we_idx]];
		}

		// For each weight compute previous, current, and next errors
		for (int we_idx = 0; we_idx < weight_count; we_idx++)
		{
			int uqw = uq_pl_weights[we_idx];

			uint32_t prev_and_next = qat->prev_next_values[uqw];
			int prev_wt_uq = prev_and_next & 0xFF;
			int next_wt_uq = (prev_and_next >> 8) & 0xFF;

			int uqw_next_dif = next_wt_uq - uqw;
			int uqw_prev_dif = prev_wt_uq - uqw;

			float current_error = 0.0f;
			float up_error = 0.0f;
			float down_error = 0.0f;

			// Interpolate the colors to create the diffs
			int texels_to_evaluate = it->weight_num_texels[we_idx];
			for (int te_idx = 0; te_idx < texels_to_evaluate; te_idx++)
			{
				int texel = it->weight_texel[we_idx][te_idx];
				const uint8_t *texel_weights = it->texel_weights_texel[we_idx][te_idx];
				const float *texel_weights_float = it->texel_weights_float_texel[we_idx][te_idx];
				float twf0 = texel_weights_float[0];
				float weight_base =
				    ((uqw * twf0
				    + uq_pl_weights[texel_weights[1]]  * texel_weights_float[1])
				    + (uq_pl_weights[texel_weights[2]] * texel_weights_float[2]
				    + uq_pl_weights[texel_weights[3]]  * texel_weights_float[3]));

				int partition = pt->partition_of_texel[texel];

				weight_base = weight_base + 0.5f;
				float plane_weight = astc::flt_rd(weight_base);
				float plane_up_weight = astc::flt_rd(weight_base + uqw_next_dif * twf0) - plane_weight;
				float plane_down_weight = astc::flt_rd(weight_base + uqw_prev_dif * twf0) - plane_weight;

				float4 color_offset = offset[partition];
				float4 color_base   = endpnt0f[partition];

				float4 color = color_base + color_offset * plane_weight;

				float4 origcolor    = float4(blk->data_r[texel], blk->data_g[texel],
				                             blk->data_b[texel], blk->data_a[texel]);
				float4 error_weight = float4(ewb->texel_weight_r[texel], ewb->texel_weight_g[texel],
				                             ewb->texel_weight_b[texel], ewb->texel_weight_a[texel]);

				float4 colordiff       = color - origcolor;
				float4 color_up_diff   = colordiff + color_offset * plane_up_weight;
				float4 color_down_diff = colordiff + color_offset * plane_down_weight;
				current_error += dot(colordiff       * colordiff,       error_weight);
				up_error      += dot(color_up_diff   * color_up_diff,   error_weight);
				down_error    += dot(color_down_diff * color_down_diff, error_weight);
			}

			// Check if the prev or next error is better, and if so use it
			if ((up_error < current_error) && (up_error < down_error))
			{
				uq_pl_weights[we_idx] = next_wt_uq;
				weight_set8[we_idx] = (uint8_t)((prev_and_next >> 24) & 0xFF);
				adjustments++;
			}
			else if (down_error < current_error)
			{
				uq_pl_weights[we_idx] = prev_wt_uq;
				weight_set8[we_idx] = (uint8_t)((prev_and_next >> 16) & 0xFF);
				adjustments++;
			}

			// IQ loss: The v1 compressor iterated here multiple times, trying
			// multiple increments or decrements until the error stopped
			// improving. This was very expensive (~15% of the v1 compressor
			// coding time) for very small improvements in quality (typically
			// less than 0.005 dB PSNR), so we now only check one step in
			// either direction
		}

		// Prepare iteration for plane 2
		weight_set8 = plane2_weight_set8;
		plane_mask ^= 0xF;
	}

	return adjustments;
}

/*
	function for compressing a block symbolically, given that we have already decided on a partition
*/
static void compress_symbolic_block_fixed_partition_1_plane(
	astc_decode_mode decode_mode,
	float mode_cutoff,
	int max_refinement_iters,
	const block_size_descriptor* bsd,
	int partition_count, int partition_index,
	const imageblock* blk,
	const error_weight_block* ewb,
	symbolic_compressed_block* scb,
	compress_fixed_partition_buffers* tmpbuf
) {
	int i, j, k;

	static const int free_bits_for_partition_count[5] = { 0, 115 - 4, 111 - 4 - PARTITION_BITS, 108 - 4 - PARTITION_BITS, 105 - 4 - PARTITION_BITS };

	const partition_info *pi = get_partition_table(bsd, partition_count);
	pi += partition_index;

	// first, compute ideal weights and endpoint colors, under the assumption that
	// there is no quantization or decimation going on.
	endpoints_and_weights *ei = tmpbuf->ei1;
	endpoints_and_weights *eix = tmpbuf->eix1;
	compute_endpoints_and_ideal_weights_1_plane(bsd, pi, blk, ewb, ei);

	// next, compute ideal weights and endpoint colors for every decimation.
	const decimation_table *const *ixtab2 = bsd->decimation_tables;

	float *decimated_quantized_weights = tmpbuf->decimated_quantized_weights;
	float *decimated_weights = tmpbuf->decimated_weights;
	float *flt_quantized_decimated_quantized_weights = tmpbuf->flt_quantized_decimated_quantized_weights;
	uint8_t *u8_quantized_decimated_quantized_weights = tmpbuf->u8_quantized_decimated_quantized_weights;

	// for each decimation mode, compute an ideal set of weights
	// (that is, weights computed with the assumption that they are not quantized)
	for (i = 0; i < MAX_DECIMATION_MODES; i++)
	{
		if (bsd->permit_encode[i] == 0 || bsd->decimation_mode_maxprec_1plane[i] < 0 || bsd->decimation_mode_percentile[i] > mode_cutoff)
			continue;
		eix[i] = *ei;
		compute_ideal_weights_for_decimation_table(&(eix[i]), ixtab2[i], decimated_quantized_weights + i * MAX_WEIGHTS_PER_BLOCK, decimated_weights + i * MAX_WEIGHTS_PER_BLOCK);

	}

	// compute maximum colors for the endpoints and ideal weights.
	// for each endpoint-and-ideal-weight pair, compute the smallest weight value
	// that will result in a color value greater than 1.
	float4 min_ep = float4(10, 10, 10, 10);
	for (i = 0; i < partition_count; i++)
	{
		#ifdef DEBUG_CAPTURE_NAN
			fedisableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif

		float4 ep = float4(
			(1.0f - ei->ep.endpt0[i].x) / (ei->ep.endpt1[i].x - ei->ep.endpt0[i].x),
			(1.0f - ei->ep.endpt0[i].y) / (ei->ep.endpt1[i].y - ei->ep.endpt0[i].y),
			(1.0f - ei->ep.endpt0[i].z) / (ei->ep.endpt1[i].z - ei->ep.endpt0[i].z),
			(1.0f - ei->ep.endpt0[i].w) / (ei->ep.endpt1[i].w - ei->ep.endpt0[i].w));

		if (ep.x > 0.5f && ep.x < min_ep.x)
			min_ep.x = ep.x;
		if (ep.y > 0.5f && ep.y < min_ep.y)
			min_ep.y = ep.y;
		if (ep.z > 0.5f && ep.z < min_ep.z)
			min_ep.z = ep.z;
		if (ep.w > 0.5f && ep.w < min_ep.w)
			min_ep.w = ep.w;

		#ifdef DEBUG_CAPTURE_NAN
			feenableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif
	}

	float min_wt_cutoff = MIN(MIN(min_ep.x, min_ep.y), MIN(min_ep.z, min_ep.w));

	// for each mode, use the angular method to compute a shift.
	float weight_low_value[MAX_WEIGHT_MODES];
	float weight_high_value[MAX_WEIGHT_MODES];

	compute_angular_endpoints_1plane(mode_cutoff, bsd, decimated_quantized_weights, decimated_weights, weight_low_value, weight_high_value);

	// for each mode (which specifies a decimation and a quantization):
	// * compute number of bits needed for the quantized weights.
	// * generate an optimized set of quantized weights.
	// * compute quantization errors for the mode.
	int qwt_bitcounts[MAX_WEIGHT_MODES];
	float qwt_errors[MAX_WEIGHT_MODES];

	for (i = 0; i < MAX_WEIGHT_MODES; i++)
	{
		if (bsd->block_modes[i].permit_encode == 0 || bsd->block_modes[i].is_dual_plane != 0 || bsd->block_modes[i].percentile > mode_cutoff)
		{
			qwt_errors[i] = 1e38f;
			continue;
		}
		if (weight_high_value[i] > 1.02f * min_wt_cutoff)
			weight_high_value[i] = 1.0f;

		int decimation_mode = bsd->block_modes[i].decimation_mode;
		if (bsd->decimation_mode_percentile[decimation_mode] > mode_cutoff)
			ASTC_CODEC_INTERNAL_ERROR();

		// compute weight bitcount for the mode
		int bits_used_by_weights = compute_ise_bitcount(ixtab2[decimation_mode]->num_weights,
														(quantization_method) bsd->block_modes[i].quantization_mode);
		int bitcount = free_bits_for_partition_count[partition_count] - bits_used_by_weights;
		if (bitcount <= 0 || bits_used_by_weights < 24 || bits_used_by_weights > 96)
		{
			qwt_errors[i] = 1e38f;
			continue;
		}
		qwt_bitcounts[i] = bitcount;

		// then, generate the optimized set of weights for the weight mode.
		compute_ideal_quantized_weights_for_decimation_table(ixtab2[decimation_mode],
															 weight_low_value[i], weight_high_value[i],
															 decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * decimation_mode,
															 flt_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * i,
															 u8_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * i,
															 bsd->block_modes[i].quantization_mode);

		// then, compute weight-errors for the weight mode.
		qwt_errors[i] = compute_error_of_weight_set(&(eix[decimation_mode]), ixtab2[decimation_mode], flt_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * i);

		#ifdef DEBUG_PRINT_DIAGNOSTICS
			if (print_diagnostics)
				printf("Block mode %d -> weight error = %f\n", i, qwt_errors[i]);
		#endif
	}

	// for each weighting mode, determine the optimal combination of color endpoint encodings
	// and weight encodings; return results for the 4 best-looking modes.

	int partition_format_specifiers[4][4];
	int quantized_weight[4];
	int color_quantization_level[4];
	int color_quantization_level_mod[4];
	determine_optimal_set_of_endpoint_formats_to_use(bsd, pi, blk, ewb, &(ei->ep), -1,	// used to flag that we are in single-weight mode
													 qwt_bitcounts, qwt_errors, partition_format_specifiers, quantized_weight, color_quantization_level, color_quantization_level_mod);

	// then iterate over the 4 believed-to-be-best modes to find out which one is
	// actually best.
	for (i = 0; i < 4; i++)
	{
		uint8_t *u8_weight_src;
		int weights_to_copy;

		if (quantized_weight[i] < 0)
		{
			scb->error_block = 1;
			scb++;
			continue;
		}

		int decimation_mode = bsd->block_modes[quantized_weight[i]].decimation_mode;
		int weight_quantization_mode = bsd->block_modes[quantized_weight[i]].quantization_mode;
		const decimation_table *it = ixtab2[decimation_mode];

		#ifdef DEBUG_PRINT_DIAGNOSTICS
			if (print_diagnostics)
			{
				printf("Selected mode = %d\n", quantized_weight[i]);
				printf("Selected decimation mode = %d\n", decimation_mode);
				printf("Selected weight-quantization mode = %d\n", weight_quantization_mode);
			}
		#endif

		u8_weight_src = u8_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * quantized_weight[i];

		weights_to_copy = it->num_weights;

		// recompute the ideal color endpoints before storing them.
		float4 rgbs_colors[4];
		float4 rgbo_colors[4];

		int l;
		for (l = 0; l < max_refinement_iters; l++)
		{
			recompute_ideal_colors(weight_quantization_mode, &(eix[decimation_mode].ep), rgbs_colors, rgbo_colors, u8_weight_src, NULL, -1, pi, it, blk, ewb);

			// quantize the chosen color

			// store the colors for the block
			for (j = 0; j < partition_count; j++)
			{
				scb->color_formats[j] = pack_color_endpoints(eix[decimation_mode].ep.endpt0[j],
															 eix[decimation_mode].ep.endpt1[j],
															 rgbs_colors[j], rgbo_colors[j], partition_format_specifiers[i][j], scb->color_values[j], color_quantization_level[i]);
			}

			// if all the color endpoint modes are the same, we get a few more
			// bits to store colors; let's see if we can take advantage of this:
			// requantize all the colors and see if the endpoint modes remain the same;
			// if they do, then exploit it.
			scb->color_formats_matched = 0;

			if ((partition_count >= 2 && scb->color_formats[0] == scb->color_formats[1]
				&& color_quantization_level[i] != color_quantization_level_mod[i])
				&& (partition_count == 2 || (scb->color_formats[0] == scb->color_formats[2] && (partition_count == 3 || (scb->color_formats[0] == scb->color_formats[3])))))
			{
				int colorvals[4][12];
				int color_formats_mod[4] = { 0 };
				for (j = 0; j < partition_count; j++)
				{
					color_formats_mod[j] = pack_color_endpoints(eix[decimation_mode].ep.endpt0[j],
																eix[decimation_mode].ep.endpt1[j],
																rgbs_colors[j], rgbo_colors[j], partition_format_specifiers[i][j], colorvals[j], color_quantization_level_mod[i]);
				}
				if (color_formats_mod[0] == color_formats_mod[1]
					&& (partition_count == 2 || (color_formats_mod[0] == color_formats_mod[2] && (partition_count == 3 || (color_formats_mod[0] == color_formats_mod[3])))))
				{
					scb->color_formats_matched = 1;
					for (j = 0; j < 4; j++)
						for (k = 0; k < 12; k++)
							scb->color_values[j][k] = colorvals[j][k];
					for (j = 0; j < 4; j++)
						scb->color_formats[j] = color_formats_mod[j];
				}
			}

			// store header fields
			scb->partition_count = partition_count;
			scb->partition_index = partition_index;
			scb->color_quantization_level = scb->color_formats_matched ? color_quantization_level_mod[i] : color_quantization_level[i];
			scb->block_mode = quantized_weight[i];
			scb->error_block = 0;

			if (scb->color_quantization_level < 4)
			{
				scb->error_block = 1;	// should never happen, but cannot prove it impossible.
			}

			// perform a final pass over the weights to try to improve them.
			int adjustments = realign_weights(decode_mode,
											  bsd,
											  blk, ewb, scb,
											  u8_weight_src,
											  NULL);

			if (adjustments == 0)
				break;
		}

		for (j = 0; j < weights_to_copy; j++)
			scb->plane1_weights[j] = u8_weight_src[j];

		scb++;
	}
}

static void compress_symbolic_block_fixed_partition_2_planes(
	astc_decode_mode decode_mode,
	float mode_cutoff,
	int max_refinement_iters,
	const block_size_descriptor* bsd,
	int partition_count,
	int partition_index,
	int separate_component,
	const imageblock* blk,
	const error_weight_block* ewb,
	symbolic_compressed_block* scb,
	compress_fixed_partition_buffers* tmpbuf
) {
	int i, j, k;

	static const int free_bits_for_partition_count[5] =
		{ 0, 113 - 4, 109 - 4 - PARTITION_BITS, 106 - 4 - PARTITION_BITS, 103 - 4 - PARTITION_BITS };

	const partition_info *pi = get_partition_table(bsd, partition_count);
	pi += partition_index;

	// first, compute ideal weights and endpoint colors
	endpoints_and_weights *ei1 = tmpbuf->ei1;
	endpoints_and_weights *ei2 = tmpbuf->ei2;
	endpoints_and_weights *eix1 = tmpbuf->eix1;
	endpoints_and_weights *eix2 = tmpbuf->eix2;
	compute_endpoints_and_ideal_weights_2_planes(bsd, pi, blk, ewb, separate_component, ei1, ei2);

	// next, compute ideal weights and endpoint colors for every decimation.
	const decimation_table *const *ixtab2 = bsd->decimation_tables;

	float *decimated_quantized_weights = tmpbuf->decimated_quantized_weights;
	float *decimated_weights = tmpbuf->decimated_weights;
	float *flt_quantized_decimated_quantized_weights = tmpbuf->flt_quantized_decimated_quantized_weights;
	uint8_t *u8_quantized_decimated_quantized_weights = tmpbuf->u8_quantized_decimated_quantized_weights;

	// for each decimation mode, compute an ideal set of weights
	for (i = 0; i < MAX_DECIMATION_MODES; i++)
	{
		if (bsd->permit_encode[i] == 0 || bsd->decimation_mode_maxprec_2planes[i] < 0 || bsd->decimation_mode_percentile[i] > mode_cutoff)
			continue;

		eix1[i] = *ei1;
		eix2[i] = *ei2;
		compute_ideal_weights_for_decimation_table(&(eix1[i]), ixtab2[i], decimated_quantized_weights + (2 * i) * MAX_WEIGHTS_PER_BLOCK, decimated_weights + (2 * i) * MAX_WEIGHTS_PER_BLOCK);
		compute_ideal_weights_for_decimation_table(&(eix2[i]), ixtab2[i], decimated_quantized_weights + (2 * i + 1) * MAX_WEIGHTS_PER_BLOCK, decimated_weights + (2 * i + 1) * MAX_WEIGHTS_PER_BLOCK);
	}

	// compute maximum colors for the endpoints and ideal weights.
	// for each endpoint-and-ideal-weight pair, compute the smallest weight value
	// that will result in a color value greater than 1.

	float4 min_ep1 = float4(10, 10, 10, 10);
	float4 min_ep2 = float4(10, 10, 10, 10);
	for (i = 0; i < partition_count; i++)
	{
		#ifdef DEBUG_CAPTURE_NAN
			fedisableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif

		float4 ep1 = float4(
			(1.0f - ei1->ep.endpt0[i].x) / (ei1->ep.endpt1[i].x - ei1->ep.endpt0[i].x),
			(1.0f - ei1->ep.endpt0[i].y) / (ei1->ep.endpt1[i].y - ei1->ep.endpt0[i].y),
			(1.0f - ei1->ep.endpt0[i].z) / (ei1->ep.endpt1[i].z - ei1->ep.endpt0[i].z),
			(1.0f - ei1->ep.endpt0[i].w) / (ei1->ep.endpt1[i].w - ei1->ep.endpt0[i].w));

		if (ep1.x > 0.5f && ep1.x < min_ep1.x)
			min_ep1.x = ep1.x;
		if (ep1.y > 0.5f && ep1.y < min_ep1.y)
			min_ep1.y = ep1.y;
		if (ep1.z > 0.5f && ep1.z < min_ep1.z)
			min_ep1.z = ep1.z;
		if (ep1.w > 0.5f && ep1.w < min_ep1.w)
			min_ep1.w = ep1.w;

		float4 ep2 = float4(
			(1.0f - ei2->ep.endpt0[i].x) / (ei2->ep.endpt1[i].x - ei2->ep.endpt0[i].x),
			(1.0f - ei2->ep.endpt0[i].y) / (ei2->ep.endpt1[i].y - ei2->ep.endpt0[i].y),
			(1.0f - ei2->ep.endpt0[i].z) / (ei2->ep.endpt1[i].z - ei2->ep.endpt0[i].z),
			(1.0f - ei2->ep.endpt0[i].w) / (ei2->ep.endpt1[i].w - ei2->ep.endpt0[i].w));

		if (ep2.x > 0.5f && ep2.x < min_ep2.x)
			min_ep2.x = ep2.x;
		if (ep2.y > 0.5f && ep2.y < min_ep2.y)
			min_ep2.y = ep2.y;
		if (ep2.z > 0.5f && ep2.z < min_ep2.z)
			min_ep2.z = ep2.z;
		if (ep2.w > 0.5f && ep2.w < min_ep2.w)
			min_ep2.w = ep2.w;

		#ifdef DEBUG_CAPTURE_NAN
			feenableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif
	}

	float min_wt_cutoff1, min_wt_cutoff2;
	switch (separate_component)
	{
	case 0:
		min_wt_cutoff2 = min_ep2.x;
		min_ep1.x = 1e30f;
		break;
	case 1:
		min_wt_cutoff2 = min_ep2.y;
		min_ep1.y = 1e30f;
		break;
	case 2:
		min_wt_cutoff2 = min_ep2.z;
		min_ep1.z = 1e30f;
		break;
	case 3:
		min_wt_cutoff2 = min_ep2.w;
		min_ep1.w = 1e30f;
		break;
	default:
		min_wt_cutoff2 = 1e30f;
	}

	min_wt_cutoff1 = MIN(MIN(min_ep1.x, min_ep1.y), MIN(min_ep1.z, min_ep1.w));

	float weight_low_value1[MAX_WEIGHT_MODES];
	float weight_high_value1[MAX_WEIGHT_MODES];
	float weight_low_value2[MAX_WEIGHT_MODES];
	float weight_high_value2[MAX_WEIGHT_MODES];

	compute_angular_endpoints_2planes(mode_cutoff, bsd, decimated_quantized_weights, decimated_weights, weight_low_value1, weight_high_value1, weight_low_value2, weight_high_value2);

	// for each mode (which specifies a decimation and a quantization):
	// * generate an optimized set of quantized weights.
	// * compute quantization errors for each mode
	// * compute number of bits needed for the quantized weights.

	int qwt_bitcounts[MAX_WEIGHT_MODES];
	float qwt_errors[MAX_WEIGHT_MODES];
	for (i = 0; i < MAX_WEIGHT_MODES; i++)
	{
		if (bsd->block_modes[i].permit_encode == 0 || bsd->block_modes[i].is_dual_plane != 1 || bsd->block_modes[i].percentile > mode_cutoff)
		{
			qwt_errors[i] = 1e38f;
			continue;
		}
		int decimation_mode = bsd->block_modes[i].decimation_mode;

		if (weight_high_value1[i] > 1.02f * min_wt_cutoff1)
			weight_high_value1[i] = 1.0f;
		if (weight_high_value2[i] > 1.02f * min_wt_cutoff2)
			weight_high_value2[i] = 1.0f;

		// compute weight bitcount for the mode
		int bits_used_by_weights = compute_ise_bitcount(2 * ixtab2[decimation_mode]->num_weights,
														(quantization_method) bsd->block_modes[i].quantization_mode);
		int bitcount = free_bits_for_partition_count[partition_count] - bits_used_by_weights;
		if (bitcount <= 0 || bits_used_by_weights < 24 || bits_used_by_weights > 96)
		{
			qwt_errors[i] = 1e38f;
			continue;
		}
		qwt_bitcounts[i] = bitcount;

		// then, generate the optimized set of weights for the mode.
		compute_ideal_quantized_weights_for_decimation_table(ixtab2[decimation_mode],
															 weight_low_value1[i],
															 weight_high_value1[i],
															 decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * decimation_mode),
															 flt_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * i),
															 u8_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * i), bsd->block_modes[i].quantization_mode);

		compute_ideal_quantized_weights_for_decimation_table(ixtab2[decimation_mode],
															 weight_low_value2[i],
															 weight_high_value2[i],
															 decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * decimation_mode + 1),
															 flt_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * i + 1),
															 u8_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * i + 1), bsd->block_modes[i].quantization_mode);


		// then, compute quantization errors for the block mode.
		qwt_errors[i] =
			compute_error_of_weight_set(&(eix1[decimation_mode]),
									   ixtab2[decimation_mode],
									   flt_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * i))
			+ compute_error_of_weight_set(&(eix2[decimation_mode]), ixtab2[decimation_mode], flt_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * i + 1));
	}

	// decide the optimal combination of color endpoint encodings and weight encodings.
	int partition_format_specifiers[4][4];
	int quantized_weight[4];
	int color_quantization_level[4];
	int color_quantization_level_mod[4];

	endpoints epm;
	merge_endpoints(&(ei1->ep), &(ei2->ep), separate_component, &epm);

	determine_optimal_set_of_endpoint_formats_to_use(bsd, pi, blk, ewb,
													 &epm, separate_component, qwt_bitcounts, qwt_errors, partition_format_specifiers, quantized_weight, color_quantization_level, color_quantization_level_mod);

	for (i = 0; i < 4; i++)
	{
		if (quantized_weight[i] < 0)
		{
			scb->error_block = 1;
			scb++;
			continue;
		}

		uint8_t *u8_weight1_src;
		uint8_t *u8_weight2_src;
		int weights_to_copy;

		int decimation_mode = bsd->block_modes[quantized_weight[i]].decimation_mode;
		int weight_quantization_mode = bsd->block_modes[quantized_weight[i]].quantization_mode;
		const decimation_table *it = ixtab2[decimation_mode];

		u8_weight1_src = u8_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * quantized_weight[i]);
		u8_weight2_src = u8_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * quantized_weight[i] + 1);

		weights_to_copy = it->num_weights;

		// recompute the ideal color endpoints before storing them.
		merge_endpoints(&(eix1[decimation_mode].ep), &(eix2[decimation_mode].ep), separate_component, &epm);

		float4 rgbs_colors[4];
		float4 rgbo_colors[4];

		int l;
		for (l = 0; l < max_refinement_iters; l++)
		{
			recompute_ideal_colors(weight_quantization_mode, &epm, rgbs_colors, rgbo_colors, u8_weight1_src, u8_weight2_src, separate_component, pi, it, blk, ewb);

			// store the colors for the block
			for (j = 0; j < partition_count; j++)
			{
				scb->color_formats[j] = pack_color_endpoints(epm.endpt0[j],
															 epm.endpt1[j],
															 rgbs_colors[j], rgbo_colors[j], partition_format_specifiers[i][j], scb->color_values[j], color_quantization_level[i]);
			}
			scb->color_formats_matched = 0;

			if ((partition_count >= 2 && scb->color_formats[0] == scb->color_formats[1]
				&& color_quantization_level[i] != color_quantization_level_mod[i])
				&& (partition_count == 2 || (scb->color_formats[0] == scb->color_formats[2] && (partition_count == 3 || (scb->color_formats[0] == scb->color_formats[3])))))
			{
				int colorvals[4][12];
				int color_formats_mod[4] = { 0 };
				for (j = 0; j < partition_count; j++)
				{
					color_formats_mod[j] = pack_color_endpoints(epm.endpt0[j],
																epm.endpt1[j],
																rgbs_colors[j], rgbo_colors[j], partition_format_specifiers[i][j], colorvals[j], color_quantization_level_mod[i]);
				}
				if (color_formats_mod[0] == color_formats_mod[1]
					&& (partition_count == 2 || (color_formats_mod[0] == color_formats_mod[2] && (partition_count == 3 || (color_formats_mod[0] == color_formats_mod[3])))))
				{
					scb->color_formats_matched = 1;
					for (j = 0; j < 4; j++)
						for (k = 0; k < 12; k++)
							scb->color_values[j][k] = colorvals[j][k];
					for (j = 0; j < 4; j++)
						scb->color_formats[j] = color_formats_mod[j];
				}
			}

			// store header fields
			scb->partition_count = partition_count;
			scb->partition_index = partition_index;
			scb->color_quantization_level = scb->color_formats_matched ? color_quantization_level_mod[i] : color_quantization_level[i];
			scb->block_mode = quantized_weight[i];
			scb->plane2_color_component = separate_component;
			scb->error_block = 0;

			if (scb->color_quantization_level < 4)
			{
				scb->error_block = 1;	// should never happen, but cannot prove it impossible
			}

			int adjustments = realign_weights(decode_mode,
											  bsd,
											  blk, ewb, scb,
											  u8_weight1_src,
											  u8_weight2_src);

			if (adjustments == 0)
				break;
		}

		for (j = 0; j < weights_to_copy; j++)
		{
			scb->plane1_weights[j] = u8_weight1_src[j];
			scb->plane2_weights[j] = u8_weight2_src[j];
		}

		scb++;
	}
}

void expand_block_artifact_suppression(
	int xdim,
	int ydim,
	int zdim,
	error_weighting_params* ewp
) {
	int x, y, z;
	float centerpos_x = (xdim - 1) * 0.5f;
	float centerpos_y = (ydim - 1) * 0.5f;
	float centerpos_z = (zdim - 1) * 0.5f;
	float *bef = ewp->block_artifact_suppression_expanded;

	for (z = 0; z < zdim; z++)
	{
		for (y = 0; y < ydim; y++)
		{
			for (x = 0; x < xdim; x++)
			{
				float xdif = (x - centerpos_x) / xdim;
				float ydif = (y - centerpos_y) / ydim;
				float zdif = (z - centerpos_z) / zdim;

				float wdif = 0.36f;
				float dist = astc::sqrt(xdif * xdif + ydif * ydif + zdif * zdif + wdif * wdif);
				*bef = powf(dist, ewp->block_artifact_suppression);
				bef++;
			}
		}
	}
}

// Function to set error weights for each color component for each texel in a block.
// Returns the sum of all the error values set.
static float prepare_error_weight_block(
	const astc_codec_image* input_image,
	const block_size_descriptor* bsd,
	const error_weighting_params* ewp,
	const imageblock* blk,
	error_weight_block* ewb,
	error_weight_block_orig* ewbo
) {
	int idx = 0;

	int any_mean_stdev_weight =
		ewp->rgb_base_weight != 1.0f || ewp->alpha_base_weight != 1.0f || \
		ewp->rgb_mean_weight != 0.0f || ewp->rgb_stdev_weight != 0.0f || \
		ewp->alpha_mean_weight != 0.0f || ewp->alpha_stdev_weight != 0.0f;

	float4 derv[MAX_TEXELS_PER_BLOCK];
	imageblock_initialize_deriv(blk, bsd->texel_count, derv);
	float4 color_weights = float4(ewp->rgba_weights[0],
	                              ewp->rgba_weights[1],
	                              ewp->rgba_weights[2],
	                              ewp->rgba_weights[3]);

	ewb->contains_zeroweight_texels = 0;

	for (int z = 0; z < bsd->zdim; z++)
	{
		for (int y = 0; y < bsd->ydim; y++)
		{
			for (int x = 0; x < bsd->xdim; x++)
			{
				int xpos = x + blk->xpos;
				int ypos = y + blk->ypos;
				int zpos = z + blk->zpos;

				if (xpos >= input_image->xsize || ypos >= input_image->ysize || zpos >= input_image->zsize)
				{
					float4 weights = float4(1e-11f, 1e-11f, 1e-11f, 1e-11f);
					ewb->error_weights[idx] = weights;
					ewb->contains_zeroweight_texels = 1;
				}
				else
				{
					float4 error_weight = float4(ewp->rgb_base_weight,
					                             ewp->rgb_base_weight,
					                             ewp->rgb_base_weight,
					                             ewp->alpha_base_weight);

					int ydt = input_image->xsize;
					int zdt = input_image->xsize * input_image->ysize;

					if (any_mean_stdev_weight)
					{
						float4 avg = input_image->input_averages[zpos * zdt + ypos * ydt + xpos];
						if (avg.x < 6e-5f)
							avg.x = 6e-5f;
						if (avg.y < 6e-5f)
							avg.y = 6e-5f;
						if (avg.z < 6e-5f)
							avg.z = 6e-5f;
						if (avg.w < 6e-5f)
							avg.w = 6e-5f;

						avg = avg * avg;

						float4 variance = input_image->input_variances[zpos * zdt + ypos * ydt + xpos];
						variance = variance * variance;

						float favg = (avg.x + avg.y + avg.z) * (1.0f / 3.0f);
						float fvar = (variance.x + variance.y + variance.z) * (1.0f / 3.0f);

						float mixing = ewp->rgb_mean_and_stdev_mixing;
						avg.x = favg * mixing + avg.x * (1.0f - mixing);
						avg.y = favg * mixing + avg.y * (1.0f - mixing);
						avg.z = favg * mixing + avg.z * (1.0f - mixing);

						variance.x = fvar * mixing + variance.x * (1.0f - mixing);
						variance.y = fvar * mixing + variance.y * (1.0f - mixing);
						variance.z = fvar * mixing + variance.z * (1.0f - mixing);

						float4 stdev = float4(astc::sqrt(MAX(variance.x, 0.0f)),
											  astc::sqrt(MAX(variance.y, 0.0f)),
											  astc::sqrt(MAX(variance.z, 0.0f)),
											  astc::sqrt(MAX(variance.w, 0.0f)));

						avg.x *= ewp->rgb_mean_weight;
						avg.y *= ewp->rgb_mean_weight;
						avg.z *= ewp->rgb_mean_weight;
						avg.w *= ewp->alpha_mean_weight;

						stdev.x *= ewp->rgb_stdev_weight;
						stdev.y *= ewp->rgb_stdev_weight;
						stdev.z *= ewp->rgb_stdev_weight;
						stdev.w = stdev.w * ewp->alpha_stdev_weight;

						error_weight = error_weight + avg + stdev;

						error_weight = float4(1.0f / error_weight.x,
											  1.0f / error_weight.y,
											  1.0f / error_weight.z,
											  1.0f / error_weight.w);
					}

					if (ewp->ra_normal_angular_scale)
					{
						// Convert from 0 to 1 to -1 to +1 range.
						float xN = (blk->orig_data[4 * idx] - 0.5f) * 2.0f;
						float yN = (blk->orig_data[4 * idx + 3] - 0.5f) * 2.0f;

						float denom = 1.0f - xN * xN - yN * yN;
						if (denom < 0.1f)
							denom = 0.1f;
						denom = 1.0f / denom;
						error_weight.x *= 1.0f + xN * xN * denom;
						error_weight.w *= 1.0f + yN * yN * denom;
					}

					if (ewp->enable_rgb_scale_with_alpha)
					{
						float alpha_scale;
						if (ewp->alpha_radius != 0)
							alpha_scale = input_image->input_alpha_averages[zpos * zdt + ypos * ydt + xpos];
						else
							alpha_scale = blk->orig_data[4 * idx + 3];
						if (alpha_scale < 0.0001f)
							alpha_scale = 0.0001f;
						alpha_scale *= alpha_scale;
						error_weight.x *= alpha_scale;
						error_weight.y *= alpha_scale;
						error_weight.z *= alpha_scale;
					}
					error_weight = error_weight * color_weights;
					error_weight = error_weight * ewp->block_artifact_suppression_expanded[idx];

					// if we perform a conversion from linear to sRGB, then we multiply
					// the weight with the derivative of the linear->sRGB transform function.
					if (perform_srgb_transform)
					{
						float r = blk->orig_data[4 * idx];
						float g = blk->orig_data[4 * idx + 1];
						float b = blk->orig_data[4 * idx + 2];
						if (r < 0.0031308f)
							r = 12.92f;
						else
							r = 0.4396f * powf(r, -0.58333f);
						if (g < 0.0031308f)
							g = 12.92f;
						else
							g = 0.4396f * powf(g, -0.58333f);
						if (b < 0.0031308f)
							b = 12.92f;
						else
							b = 0.4396f * powf(b, -0.58333f);
						error_weight.x *= r;
						error_weight.y *= g;
						error_weight.z *= b;
					}

					// when we loaded the block to begin with, we applied a transfer function
					// and computed the derivative of the transfer function. However, the
					// error-weight computation so far is based on the original color values,
					// not the transfer-function values. As such, we must multiply the
					// error weights by the derivative of the inverse of the transfer function,
					// which is equivalent to dividing by the derivative of the transfer
					// function.

					ewbo->error_weights[idx] = error_weight;

					error_weight.x /= (derv[idx].x * derv[idx].x * 1e-10f);
					error_weight.y /= (derv[idx].y * derv[idx].y * 1e-10f);
					error_weight.z /= (derv[idx].z * derv[idx].z * 1e-10f);
					error_weight.w /= (derv[idx].w * derv[idx].w * 1e-10f);

					ewb->error_weights[idx] = error_weight;
					if (dot(error_weight, float4(1.0f, 1.0f, 1.0f, 1.0f)) < 1e-10f)
						ewb->contains_zeroweight_texels = 1;
				}
				idx++;
			}
		}
	}

	float4 error_weight_sum = float4(0.0f, 0.0f, 0.0f, 0.0f);
	int texels_per_block = bsd->texel_count;
	for (int i = 0; i < texels_per_block; i++)
	{
		error_weight_sum = error_weight_sum + ewb->error_weights[i];

		ewb->texel_weight_r[i] = ewb->error_weights[i].x;
		ewb->texel_weight_g[i] = ewb->error_weights[i].y;
		ewb->texel_weight_b[i] = ewb->error_weights[i].z;
		ewb->texel_weight_a[i] = ewb->error_weights[i].w;

		ewb->texel_weight_rg[i] = (ewb->error_weights[i].x + ewb->error_weights[i].y) * 0.5f;
		ewb->texel_weight_rb[i] = (ewb->error_weights[i].x + ewb->error_weights[i].z) * 0.5f;
		ewb->texel_weight_gb[i] = (ewb->error_weights[i].y + ewb->error_weights[i].z) * 0.5f;
		ewb->texel_weight_ra[i] = (ewb->error_weights[i].x + ewb->error_weights[i].w) * 0.5f;

		ewb->texel_weight_gba[i] = (ewb->error_weights[i].y + ewb->error_weights[i].z + ewb->error_weights[i].w) * 0.333333f;
		ewb->texel_weight_rba[i] = (ewb->error_weights[i].x + ewb->error_weights[i].z + ewb->error_weights[i].w) * 0.333333f;
		ewb->texel_weight_rga[i] = (ewb->error_weights[i].x + ewb->error_weights[i].y + ewb->error_weights[i].w) * 0.333333f;
		ewb->texel_weight_rgb[i] = (ewb->error_weights[i].x + ewb->error_weights[i].y + ewb->error_weights[i].z) * 0.333333f;
		ewb->texel_weight[i] = (ewb->error_weights[i].x + ewb->error_weights[i].y + ewb->error_weights[i].z + ewb->error_weights[i].w) * 0.25f;
	}

	return dot(error_weight_sum, float4(1.0f, 1.0f, 1.0f, 1.0f));
}

static void prepare_block_statistics(
	int texels_per_block,
	const imageblock * blk,
	const error_weight_block* ewb,
	int* is_normal_map,
	float* lowest_correl
) {
	int i;

	// compute covariance matrix, as a collection of 10 scalars
	// (that form the upper-triangular row of the matrix; the matrix is
	// symmetric, so this is all we need)
	float rs = 0.0f;
	float gs = 0.0f;
	float bs = 0.0f;
	float as = 0.0f;
	float rr_var = 0.0f;
	float gg_var = 0.0f;
	float bb_var = 0.0f;
	float aa_var = 0.0f;
	float rg_cov = 0.0f;
	float rb_cov = 0.0f;
	float ra_cov = 0.0f;
	float gb_cov = 0.0f;
	float ga_cov = 0.0f;
	float ba_cov = 0.0f;

	float weight_sum = 0.0f;

	for (i = 0; i < texels_per_block; i++)
	{
		float weight = ewb->texel_weight[i];
		if (weight < 0.0f)
			ASTC_CODEC_INTERNAL_ERROR();
		weight_sum += weight;

		float r = blk->data_r[i];
		float g = blk->data_g[i];
		float b = blk->data_b[i];
		float a = blk->data_a[i];

		float rw = r * weight;
		rs += rw;
		rr_var += r * rw;
		rg_cov += g * rw;
		rb_cov += b * rw;
		ra_cov += a * rw;

		float gw = g * weight;
		gs += gw;
		gg_var += g * gw;
		gb_cov += b * gw;
		ga_cov += a * gw;

		float bw = b * weight;
		bs += bw;
		bb_var += b * bw;
		ba_cov += a * bw;

		float aw = a * weight;
		as += aw;
		aa_var += a * aw;
	}

	float rpt = 1.0f / MAX(weight_sum, 1e-7f);

	rr_var -= rs * (rs * rpt);
	rg_cov -= gs * (rs * rpt);
	rb_cov -= bs * (rs * rpt);
	ra_cov -= as * (rs * rpt);

	gg_var -= gs * (gs * rpt);
	gb_cov -= bs * (gs * rpt);
	ga_cov -= as * (gs * rpt);

	bb_var -= bs * (bs * rpt);
	ba_cov -= as * (bs * rpt);

	aa_var -= as * (as * rpt);

	rg_cov *= astc::rsqrt(MAX(rr_var * gg_var, 1e-30f));
	rb_cov *= astc::rsqrt(MAX(rr_var * bb_var, 1e-30f));
	ra_cov *= astc::rsqrt(MAX(rr_var * aa_var, 1e-30f));
	gb_cov *= astc::rsqrt(MAX(gg_var * bb_var, 1e-30f));
	ga_cov *= astc::rsqrt(MAX(gg_var * aa_var, 1e-30f));
	ba_cov *= astc::rsqrt(MAX(bb_var * aa_var, 1e-30f));

	if (astc::isnan(rg_cov)) rg_cov = 1.0f;
	if (astc::isnan(rb_cov)) rb_cov = 1.0f;
	if (astc::isnan(ra_cov)) ra_cov = 1.0f;
	if (astc::isnan(gb_cov)) gb_cov = 1.0f;
	if (astc::isnan(ga_cov)) ga_cov = 1.0f;
	if (astc::isnan(ba_cov)) ba_cov = 1.0f;

	float lowest_correlation = MIN(fabsf(rg_cov), fabsf(rb_cov));
	lowest_correlation = MIN(lowest_correlation, fabsf(ra_cov));
	lowest_correlation = MIN(lowest_correlation, fabsf(gb_cov));
	lowest_correlation = MIN(lowest_correlation, fabsf(ga_cov));
	lowest_correlation = MIN(lowest_correlation, fabsf(ba_cov));
	*lowest_correl = lowest_correlation;

	// compute a "normal-map" factor
	// this factor should be exactly 0.0 for a normal map, while it may be all over the
	// place for anything that is NOT a normal map. We can probably assume that a factor
	// of less than 0.2f represents a normal map.

	float nf_sum = 0.0f;
	for (i = 0; i < texels_per_block; i++)
	{
		float3 val = float3(blk->orig_data[4 * i],
							blk->orig_data[4 * i + 1],
							blk->orig_data[4 * i + 2]);
		val = (val - float3(0.5f, 0.5f, 0.5f)) * 2.0f;
		float length_squared = dot(val, val);
		float nf = fabsf(length_squared - 1.0f);
		nf_sum += nf;
	}
	*is_normal_map = nf_sum < (0.2f * (float)texels_per_block);
}

#ifdef DEBUG_PRINT_DIAGNOSTICS
	int block_mode_histogram[2048];
#endif

float compress_symbolic_block(
	const astc_codec_image* input_image,
	astc_decode_mode decode_mode,
	const block_size_descriptor* bsd,
	const error_weighting_params* ewp,
	const imageblock* blk,
	symbolic_compressed_block* scb,
	compress_symbolic_block_buffers* tmpbuf)
{
	int i, j;
	int xpos = blk->xpos;
	int ypos = blk->ypos;
	int zpos = blk->zpos;

	#ifdef DEBUG_PRINT_DIAGNOSTICS
		if (print_diagnostics)
		{
			printf("Diagnostics of block of dimension %d x %d x %d\n\n", xdim, ydim, zdim);

			printf("XPos: %d  YPos: %d  ZPos: %d\n", xpos, ypos, zpos);

			printf("Red-min: %f   Red-max: %f\n", blk->red_min, blk->red_max);
			printf("Green-min: %f   Green-max: %f\n", blk->green_min, blk->green_max);
			printf("Blue-min: %f   Blue-max: %f\n", blk->blue_min, blk->blue_max);
			printf("Alpha-min: %f   Alpha-max: %f\n", blk->alpha_min, blk->alpha_max);
			printf("Grayscale: %d\n", blk->grayscale);

			for (int z = 0; z < zdim; z++)
				for (int y = 0; y < ydim; y++)
					for (int x = 0; x < xdim; x++)
					{
						int idx = ((z * ydim + y) * xdim + x) * 4;
						printf("Texel (%d %d %d) : orig=< %g, %g, %g, %g >, work=< %g, %g, %g, %g >\n",
							x, y, z,
							blk->orig_data[idx], blk->orig_data[idx + 1], blk->orig_data[idx + 2], blk->orig_data[idx + 3],
							blk->data_r[idx], blk->data_g[idx], blk->data_b[idx], blk->data_a[idx]);
					}
			printf("\n");
		}
	#endif

	if (blk->red_min == blk->red_max && blk->green_min == blk->green_max && blk->blue_min == blk->blue_max && blk->alpha_min == blk->alpha_max)
	{
		// detected a constant-color block. Encode as FP16 if using HDR
		scb->error_block = 0;

		if (rgb_force_use_of_hdr)
		{
			scb->block_mode = -1;
			scb->partition_count = 0;
			scb->constant_color[0] = float_to_sf16(blk->orig_data[0], SF_NEARESTEVEN);
			scb->constant_color[1] = float_to_sf16(blk->orig_data[1], SF_NEARESTEVEN);
			scb->constant_color[2] = float_to_sf16(blk->orig_data[2], SF_NEARESTEVEN);
			scb->constant_color[3] = float_to_sf16(blk->orig_data[3], SF_NEARESTEVEN);
		}
		else
		{
			// Encode as UNORM16 if NOT using HDR.
			scb->block_mode = -2;
			scb->partition_count = 0;
			float red = blk->orig_data[0];
			float green = blk->orig_data[1];
			float blue = blk->orig_data[2];
			float alpha = blk->orig_data[3];

			if (red < 0)
				red = 0;
			else if (red > 1)
				red = 1;

			if (green < 0)
				green = 0;
			else if (green > 1)
				green = 1;

			if (blue < 0)
				blue = 0;
			else if (blue > 1)
				blue = 1;

			if (alpha < 0)
				alpha = 0;
			else if (alpha > 1)
				alpha = 1;

			scb->constant_color[0] = astc::flt2int_rtn(red * 65535.0f);
			scb->constant_color[1] = astc::flt2int_rtn(green * 65535.0f);
			scb->constant_color[2] = astc::flt2int_rtn(blue * 65535.0f);
			scb->constant_color[3] = astc::flt2int_rtn(alpha * 65535.0f);
		}

		#ifdef DEBUG_PRINT_DIAGNOSTICS
			if (print_diagnostics)
			{
				printf("Block is single-color <%4.4X %4.4X %4.4X %4.4X>\n", scb->constant_color[0], scb->constant_color[1], scb->constant_color[2], scb->constant_color[3]);
			}

			if (print_tile_errors)
				printf("0\n");
		#endif

		physical_compressed_block psb = symbolic_to_physical(bsd, scb);
		physical_to_symbolic(bsd, psb, scb);

		return 0.0f;
	}

	error_weight_block *ewb = tmpbuf->ewb;
	error_weight_block_orig *ewbo = tmpbuf->ewbo;

	float error_weight_sum = prepare_error_weight_block(input_image, bsd, ewp, blk, ewb, ewbo);

	#ifdef DEBUG_PRINT_DIAGNOSTICS
		if (print_diagnostics)
		{
			printf("\n");
			for (int z = 0; z < zdim; z++)
				for (int y = 0; y < ydim; y++)
					for (int x = 0; x < xdim; x++)
					{
						int idx = (z * ydim + y) * xdim + x;
						printf("ErrorWeight (%d %d %d) : < %g, %g, %g, %g >\n", x, y, z, ewb->error_weights[idx].x, ewb->error_weights[idx].y, ewb->error_weights[idx].z, ewb->error_weights[idx].w);
					}
			printf("\n");
		}
	#endif

	symbolic_compressed_block *tempblocks = tmpbuf->tempblocks;

	float error_of_best_block = 1e20f;

	imageblock *temp = tmpbuf->temp;

	float best_errorvals_in_modes[17];
	for (i = 0; i < 17; i++)
		best_errorvals_in_modes[i] = 1e30f;

	int uses_alpha = imageblock_uses_alpha(blk);

	float mode_cutoff = ewp->block_mode_cutoff;

	// next, test mode #0. This mode uses 1 plane of weights and 1 partition.
	// we test it twice, first with a modecutoff of 0, then with the specified mode-cutoff.
	// This causes an early-out that speeds up encoding of "easy" content.

	float modecutoffs[2];
	float errorval_mult[2] = { 2.5, 1 };
	modecutoffs[0] = 0;
	modecutoffs[1] = mode_cutoff;

	float best_errorval_in_mode;
	for (i = 0; i < 2; i++)
	{
		compress_symbolic_block_fixed_partition_1_plane(decode_mode, modecutoffs[i], ewp->max_refinement_iters, bsd, 1,	// partition count
														0,	// partition index
														blk, ewb, tempblocks, tmpbuf->plane1);

		best_errorval_in_mode = 1e30f;
		for (j = 0; j < 4; j++)
		{
			if (tempblocks[j].error_block)
				continue;
			decompress_symbolic_block(decode_mode, bsd, xpos, ypos, zpos, tempblocks + j, temp);
			float errorval = compute_imageblock_difference(bsd, blk, temp, ewb) * errorval_mult[i];

			#ifdef DEBUG_PRINT_DIAGNOSTICS
				if (print_diagnostics)
				{
					printf("\n-----------------------------------\n");
					printf("Single-weight partition test 0 (1 partition) completed\n");
					printf("Resulting error value: %g\n", errorval);
				}
			#endif

			if (errorval < best_errorval_in_mode)
				best_errorval_in_mode = errorval;

			if (errorval < error_of_best_block)
			{
				#ifdef DEBUG_PRINT_DIAGNOSTICS
					if (print_diagnostics)
						printf("Accepted as better than previous-best-error, which was %g\n", error_of_best_block);
				#endif

				error_of_best_block = errorval;
				*scb = tempblocks[j];
			}

			#ifdef DEBUG_PRINT_DIAGNOSTICS
				if (print_diagnostics)
				{
					printf("-----------------------------------\n");
				}
			#endif
		}

		best_errorvals_in_modes[0] = best_errorval_in_mode;
		if ((error_of_best_block / error_weight_sum) < ewp->texel_avg_error_limit)
			goto END_OF_TESTS;
	}

	int is_normal_map;
	float lowest_correl;
	prepare_block_statistics(bsd->texel_count, blk, ewb, &is_normal_map, &lowest_correl);

	if (is_normal_map && lowest_correl < 0.99f)
		lowest_correl = 0.99f;

	// next, test the four possible 1-partition, 2-planes modes
	for (i = 0; i < 4; i++)
	{

		if (lowest_correl > ewp->lowest_correlation_cutoff)
			continue;

		if (blk->grayscale && i != 3)
			continue;

		if (!uses_alpha && i == 3)
			continue;

		compress_symbolic_block_fixed_partition_2_planes(decode_mode, mode_cutoff, ewp->max_refinement_iters,
														 bsd, 1,	// partition count
														 0,	// partition index
														 i,	// the color component to test a separate plane of weights for.
														 blk, ewb, tempblocks, tmpbuf->planes2);

		best_errorval_in_mode = 1e30f;
		for (j = 0; j < 4; j++)
		{
			if (tempblocks[j].error_block)
				continue;
			decompress_symbolic_block(decode_mode, bsd, xpos, ypos, zpos, tempblocks + j, temp);
			float errorval = compute_imageblock_difference(bsd, blk, temp, ewb);

			#ifdef DEBUG_PRINT_DIAGNOSTICS
				if (print_diagnostics)
				{
					printf("\n-----------------------------------\n");
					printf("Dual-weight partition test %d (1 partition) completed\n", i);
					printf("Resulting error value: %g\n", errorval);
				}
			#endif

			if (errorval < best_errorval_in_mode)
				best_errorval_in_mode = errorval;

			if (errorval < error_of_best_block)
			{
				#ifdef DEBUG_PRINT_DIAGNOSTICS
					if (print_diagnostics)
						printf("Accepted as better than previous-best-error, which was %g\n", error_of_best_block);
				#endif

				error_of_best_block = errorval;
				*scb = tempblocks[j];
			}

			#ifdef DEBUG_PRINT_DIAGNOSTICS
				if (print_diagnostics)
				{
					printf("-----------------------------------\n");
				}
			#endif

			best_errorvals_in_modes[i + 1] = best_errorval_in_mode;
		}

		if ((error_of_best_block / error_weight_sum) < ewp->texel_avg_error_limit)
			goto END_OF_TESTS;
	}

	// find best blocks for 2, 3 and 4 partitions
	int partition_count;
	for (partition_count = 2; partition_count <= 4; partition_count++)
	{
		int partition_indices_1plane[2];
		int partition_indices_2planes[2];

		find_best_partitionings(ewp->partition_search_limit,
								bsd, partition_count, blk, ewb, 1,
								&(partition_indices_1plane[0]), &(partition_indices_1plane[1]), &(partition_indices_2planes[0]));

		for (i = 0; i < 2; i++)
		{
			compress_symbolic_block_fixed_partition_1_plane(decode_mode, mode_cutoff, ewp->max_refinement_iters,
															bsd, partition_count, partition_indices_1plane[i], blk, ewb, tempblocks, tmpbuf->plane1);

			best_errorval_in_mode = 1e30f;
			for (j = 0; j < 4; j++)
			{
				if (tempblocks[j].error_block)
					continue;
				decompress_symbolic_block(decode_mode, bsd, xpos, ypos, zpos, tempblocks + j, temp);
				float errorval = compute_imageblock_difference(bsd, blk, temp, ewb);

				#ifdef DEBUG_PRINT_DIAGNOSTICS
					if (print_diagnostics)
					{
						printf("\n-----------------------------------\n");
						printf("Single-weight partition test %d (%d partitions) completed\n", i, partition_count);
						printf("Resulting error value: %g\n", errorval);
					}
				#endif

				if (errorval < best_errorval_in_mode)
					best_errorval_in_mode = errorval;

				if (errorval < error_of_best_block)
				{
					#ifdef DEBUG_PRINT_DIAGNOSTICS
						if (print_diagnostics)
							printf("Accepted as better than previous-best-error, which was %g\n", error_of_best_block);
					#endif

					error_of_best_block = errorval;
					*scb = tempblocks[j];
				}
			}

			best_errorvals_in_modes[4 * (partition_count - 2) + 5 + i] = best_errorval_in_mode;

			#ifdef DEBUG_PRINT_DIAGNOSTICS
				if (print_diagnostics)
				{
					printf("-----------------------------------\n");
				}
			#endif

			if ((error_of_best_block / error_weight_sum) < ewp->texel_avg_error_limit)
				goto END_OF_TESTS;
		}


		if (partition_count == 2 && !is_normal_map && MIN(best_errorvals_in_modes[5], best_errorvals_in_modes[6]) > (best_errorvals_in_modes[0] * ewp->partition_1_to_2_limit))
			goto END_OF_TESTS;

		// don't bother to check 4 partitions for dual plane of weights, ever.
		if (partition_count == 4)
			break;

		for (i = 0; i < 2; i++)
		{
			if (lowest_correl > ewp->lowest_correlation_cutoff)
				continue;
			compress_symbolic_block_fixed_partition_2_planes(decode_mode,
															 mode_cutoff,
															 ewp->max_refinement_iters,
															 bsd,
															 partition_count,
															 partition_indices_2planes[i] & (PARTITION_COUNT - 1), partition_indices_2planes[i] >> PARTITION_BITS,
															 blk, ewb, tempblocks, tmpbuf->planes2);

			best_errorval_in_mode = 1e30f;
			for (j = 0; j < 4; j++)
			{
				if (tempblocks[j].error_block)
					continue;
				decompress_symbolic_block(decode_mode, bsd, xpos, ypos, zpos, tempblocks + j, temp);

				float errorval = compute_imageblock_difference(bsd, blk, temp, ewb);

				#ifdef DEBUG_PRINT_DIAGNOSTICS
					if (print_diagnostics)
					{
						printf("\n-----------------------------------\n");
						printf("Dual-weight partition test %d (%d partitions) completed\n", i, partition_count);
						printf("Resulting error value: %g\n", errorval);
					}
				#endif

				if (errorval < best_errorval_in_mode)
					best_errorval_in_mode = errorval;

				if (errorval < error_of_best_block)
				{
					#ifdef DEBUG_PRINT_DIAGNOSTICS
						if (print_diagnostics)
							printf("Accepted as better than previous-best-error, which was %g\n", error_of_best_block);
					#endif

					error_of_best_block = errorval;
					*scb = tempblocks[j];
				}
			}

			best_errorvals_in_modes[4 * (partition_count - 2) + 5 + 2 + i] = best_errorval_in_mode;

			#ifdef DEBUG_PRINT_DIAGNOSTICS
				if (print_diagnostics)
				{
					printf("-----------------------------------\n");
				}
			#endif

			if ((error_of_best_block / error_weight_sum) < ewp->texel_avg_error_limit)
				goto END_OF_TESTS;
		}
	}

END_OF_TESTS:

	#ifdef DEBUG_PRINT_DIAGNOSTICS
		if (scb->block_mode >= 0)
			block_mode_histogram[scb->block_mode & 0x7ff]++;

		if (print_tile_errors)
			printf("%g\n", (double)error_of_best_block);
	#endif

	// compress/decompress to a physical block
	physical_compressed_block psb = symbolic_to_physical(bsd, scb);
	physical_to_symbolic(bsd, psb, scb);

	// mean squared error per color component.
	return error_of_best_block / ((float)(bsd->xdim * bsd->ydim * bsd->zdim));
}
