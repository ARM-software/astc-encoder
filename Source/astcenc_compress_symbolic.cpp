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

#if !defined(ASTCENC_DECOMPRESS_ONLY)

/**
 * @brief Functions to compress a symbolic block.
 */

#include "astcenc_internal.h"

#include <cassert>
#include <cstring>
#include <cstdio>

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
	astcenc_profile decode_mode,
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
	const int packed_index = bsd->block_mode_to_packed[scb->block_mode];
	assert(packed_index >= 0 && packed_index < bsd->block_mode_packed_count);
	const block_mode& bm = bsd->block_modes_packed[packed_index];
	int weight_quantization_level = bm.quantization_mode;
	const quantization_and_transfer_table *qat = &(quant_and_xfer_tables[weight_quantization_level]);

	// Get the decimation table
	const decimation_table *const *ixtab2 = bsd->decimation_tables;
	const decimation_table *it = ixtab2[bm.decimation_mode];
	int weight_count = it->num_weights;

	int max_plane = bm.is_dual_plane;
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

			if (plane_mask & 1) epd.r = 0;
			if (plane_mask & 2) epd.g = 0;
			if (plane_mask & 4) epd.b = 0;
			if (plane_mask & 8) epd.a = 0;

			endpnt0f[pa_idx] = float4((float)endpnt0[pa_idx].r, (float)endpnt0[pa_idx].g,
			                          (float)endpnt0[pa_idx].b, (float)endpnt0[pa_idx].a);
			offset[pa_idx] = float4((float)epd.r, (float)epd.g, (float)epd.b, (float)epd.a);
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
	astcenc_profile decode_mode,
	float mode_cutoff,
	int tune_candidate_limit,
	int max_refinement_iters,
	const block_size_descriptor* bsd,
	int partition_count, int partition_index,
	const imageblock* blk,
	const error_weight_block* ewb,
	symbolic_compressed_block* scb,
	compress_fixed_partition_buffers* tmpbuf
) {
	static const int free_bits_for_partition_count[5] = { 0, 115 - 4, 111 - 4 - PARTITION_BITS, 108 - 4 - PARTITION_BITS, 105 - 4 - PARTITION_BITS };

	const partition_info *pi = get_partition_table(bsd, partition_count);
	pi += partition_index;

	// first, compute ideal weights and endpoint colors, under the assumption that
	// there is no quantization or decimation going on.
	endpoints_and_weights *ei = &tmpbuf->ei1;
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
	for (int i = 0; i < MAX_DECIMATION_MODES; i++)
	{
		if (bsd->permit_encode[i] == 0 || bsd->decimation_mode_maxprec_1plane[i] < 0 || bsd->decimation_mode_percentile[i] > mode_cutoff)
		{
			continue;
		}
		eix[i] = *ei;
		compute_ideal_weights_for_decimation_table(&(eix[i]), ixtab2[i], decimated_quantized_weights + i * MAX_WEIGHTS_PER_BLOCK, decimated_weights + i * MAX_WEIGHTS_PER_BLOCK);

	}

	// compute maximum colors for the endpoints and ideal weights.
	// for each endpoint-and-ideal-weight pair, compute the smallest weight value
	// that will result in a color value greater than 1.
	float4 min_ep = float4(10.0f);
	for (int i = 0; i < partition_count; i++)
	{
		#ifdef DEBUG_CAPTURE_NAN
			fedisableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif

		float4 ep = float4(
			(1.0f - ei->ep.endpt0[i].r) / (ei->ep.endpt1[i].r - ei->ep.endpt0[i].r),
			(1.0f - ei->ep.endpt0[i].g) / (ei->ep.endpt1[i].g - ei->ep.endpt0[i].g),
			(1.0f - ei->ep.endpt0[i].b) / (ei->ep.endpt1[i].b - ei->ep.endpt0[i].b),
			(1.0f - ei->ep.endpt0[i].a) / (ei->ep.endpt1[i].a - ei->ep.endpt0[i].a));

		if (ep.r > 0.5f && ep.r < min_ep.r)
		{
			min_ep.r = ep.r;
		}

		if (ep.g > 0.5f && ep.g < min_ep.g)
		{
			min_ep.g = ep.g;
		}

		if (ep.b > 0.5f && ep.b < min_ep.b)
		{
			min_ep.b = ep.b;
		}

		if (ep.a > 0.5f && ep.a < min_ep.a)
		{
			min_ep.a = ep.a;
		}

		#ifdef DEBUG_CAPTURE_NAN
			feenableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif
	}

	float min_wt_cutoff = MIN(MIN(min_ep.r, min_ep.g), MIN(min_ep.b, min_ep.a));

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

	for (int i = 0, ni = bsd->block_mode_packed_count; i < ni; ++i)
	{
		const block_mode& bm = bsd->block_modes_packed[i];
		if (bm.is_dual_plane != 0 || bm.percentile > mode_cutoff)
		{
			qwt_errors[i] = 1e38f;
			continue;
		}

		if (weight_high_value[i] > 1.02f * min_wt_cutoff)
		{
			weight_high_value[i] = 1.0f;
		}

		int decimation_mode = bm.decimation_mode;

		// compute weight bitcount for the mode
		int bits_used_by_weights = compute_ise_bitcount(ixtab2[decimation_mode]->num_weights,
														(quantization_method) bm.quantization_mode);
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
															 bm.quantization_mode);

		// then, compute weight-errors for the weight mode.
		qwt_errors[i] = compute_error_of_weight_set(&(eix[decimation_mode]), ixtab2[decimation_mode], flt_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * i);
	}

	// for each weighting mode, determine the optimal combination of color endpoint encodings
	// and weight encodings; return results for the 4 best-looking modes.

	int partition_format_specifiers[TUNE_MAX_TRIAL_CANDIDATES][4];
	int quantized_weight[TUNE_MAX_TRIAL_CANDIDATES];
	int color_quantization_level[TUNE_MAX_TRIAL_CANDIDATES];
	int color_quantization_level_mod[TUNE_MAX_TRIAL_CANDIDATES];
	determine_optimal_set_of_endpoint_formats_to_use(
	    bsd, pi, blk, ewb, &(ei->ep), -1, qwt_bitcounts, qwt_errors,
	    tune_candidate_limit, partition_format_specifiers, quantized_weight,
	    color_quantization_level, color_quantization_level_mod);

	// then iterate over the tune_candidate_limit believed-to-be-best modes to
	// find out which one is actually best.
	for (int i = 0; i < tune_candidate_limit; i++)
	{
		uint8_t *u8_weight_src;
		int weights_to_copy;

		const int qw_packed_index = quantized_weight[i];
		if (qw_packed_index < 0)
		{
			scb->error_block = 1;
			scb++;
			continue;
		}

		assert(qw_packed_index >= 0 && qw_packed_index < bsd->block_mode_packed_count);
		const block_mode& qw_bm = bsd->block_modes_packed[qw_packed_index];

		int decimation_mode = qw_bm.decimation_mode;
		int weight_quantization_mode = qw_bm.quantization_mode;
		const decimation_table *it = ixtab2[decimation_mode];
		u8_weight_src = u8_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * qw_packed_index;

		weights_to_copy = it->num_weights;

		// recompute the ideal color endpoints before storing them.
		float4 rgbs_colors[4];
		float4 rgbo_colors[4];

		for (int l = 0; l < max_refinement_iters; l++)
		{
			recompute_ideal_colors(weight_quantization_mode, &(eix[decimation_mode].ep), rgbs_colors, rgbo_colors, u8_weight_src, nullptr, -1, pi, it, blk, ewb);

			// quantize the chosen color

			// store the colors for the block
			for (int j = 0; j < partition_count; j++)
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
				for (int j = 0; j < partition_count; j++)
				{
					color_formats_mod[j] = pack_color_endpoints(eix[decimation_mode].ep.endpt0[j],
																eix[decimation_mode].ep.endpt1[j],
																rgbs_colors[j], rgbo_colors[j], partition_format_specifiers[i][j], colorvals[j], color_quantization_level_mod[i]);
				}
				if (color_formats_mod[0] == color_formats_mod[1]
					&& (partition_count == 2 || (color_formats_mod[0] == color_formats_mod[2] && (partition_count == 3 || (color_formats_mod[0] == color_formats_mod[3])))))
				{
					scb->color_formats_matched = 1;
					for (int j = 0; j < 4; j++)
					{
						for (int k = 0; k < 12; k++)
						{
							scb->color_values[j][k] = colorvals[j][k];
						}
					}

					for (int j = 0; j < 4; j++)
					{
						scb->color_formats[j] = color_formats_mod[j];
					}
				}
			}

			// store header fields
			scb->partition_count = partition_count;
			scb->partition_index = partition_index;
			scb->color_quantization_level = scb->color_formats_matched ? color_quantization_level_mod[i] : color_quantization_level[i];
			scb->block_mode = qw_bm.mode_index;
			scb->error_block = 0;

			if (scb->color_quantization_level < 4)
			{
				scb->error_block = 1;	// should never happen, but cannot prove it impossible.
			}

			// perform a final pass over the weights to try to improve them.
			int adjustments = realign_weights(
				decode_mode, bsd, blk, ewb, scb, u8_weight_src, nullptr);

			if (adjustments == 0)
			{
				break;
			}
		}

		for (int j = 0; j < weights_to_copy; j++)
		{
			scb->plane1_weights[j] = u8_weight_src[j];
		}

		scb++;
	}
}

static void compress_symbolic_block_fixed_partition_2_planes(
	astcenc_profile decode_mode,
	float mode_cutoff,
	int tune_candidate_limit,
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
	static const int free_bits_for_partition_count[5] =
		{ 0, 113 - 4, 109 - 4 - PARTITION_BITS, 106 - 4 - PARTITION_BITS, 103 - 4 - PARTITION_BITS };

	const partition_info *pi = get_partition_table(bsd, partition_count);
	pi += partition_index;

	// first, compute ideal weights and endpoint colors
	endpoints_and_weights *ei1 = &tmpbuf->ei1;
	endpoints_and_weights *ei2 = &tmpbuf->ei2;
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
	for (int i = 0; i < MAX_DECIMATION_MODES; i++)
	{
		if (bsd->permit_encode[i] == 0 || bsd->decimation_mode_maxprec_2planes[i] < 0 || bsd->decimation_mode_percentile[i] > mode_cutoff)
		{
			continue;
		}

		eix1[i] = *ei1;
		eix2[i] = *ei2;
		compute_ideal_weights_for_decimation_table(&(eix1[i]), ixtab2[i], decimated_quantized_weights + (2 * i) * MAX_WEIGHTS_PER_BLOCK, decimated_weights + (2 * i) * MAX_WEIGHTS_PER_BLOCK);
		compute_ideal_weights_for_decimation_table(&(eix2[i]), ixtab2[i], decimated_quantized_weights + (2 * i + 1) * MAX_WEIGHTS_PER_BLOCK, decimated_weights + (2 * i + 1) * MAX_WEIGHTS_PER_BLOCK);
	}

	// compute maximum colors for the endpoints and ideal weights.
	// for each endpoint-and-ideal-weight pair, compute the smallest weight value
	// that will result in a color value greater than 1.

	float4 min_ep1 = float4(10.0f);
	float4 min_ep2 = float4(10.0f);
	for (int i = 0; i < partition_count; i++)
	{
		#ifdef DEBUG_CAPTURE_NAN
			fedisableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif

		float4 ep1 = float4(
			(1.0f - ei1->ep.endpt0[i].r) / (ei1->ep.endpt1[i].r - ei1->ep.endpt0[i].r),
			(1.0f - ei1->ep.endpt0[i].g) / (ei1->ep.endpt1[i].g - ei1->ep.endpt0[i].g),
			(1.0f - ei1->ep.endpt0[i].b) / (ei1->ep.endpt1[i].b - ei1->ep.endpt0[i].b),
			(1.0f - ei1->ep.endpt0[i].a) / (ei1->ep.endpt1[i].a - ei1->ep.endpt0[i].a));

		if (ep1.r > 0.5f && ep1.r < min_ep1.r)
		{
			min_ep1.r = ep1.r;
		}

		if (ep1.g > 0.5f && ep1.g < min_ep1.g)
		{
			min_ep1.g = ep1.g;
		}

		if (ep1.b > 0.5f && ep1.b < min_ep1.b)
		{
			min_ep1.b = ep1.b;
		}

		if (ep1.a > 0.5f && ep1.a < min_ep1.a)
		{
			min_ep1.a = ep1.a;
		}

		float4 ep2 = float4(
			(1.0f - ei2->ep.endpt0[i].r) / (ei2->ep.endpt1[i].r - ei2->ep.endpt0[i].r),
			(1.0f - ei2->ep.endpt0[i].g) / (ei2->ep.endpt1[i].g - ei2->ep.endpt0[i].g),
			(1.0f - ei2->ep.endpt0[i].b) / (ei2->ep.endpt1[i].b - ei2->ep.endpt0[i].b),
			(1.0f - ei2->ep.endpt0[i].a) / (ei2->ep.endpt1[i].a - ei2->ep.endpt0[i].a));

		if (ep2.r > 0.5f && ep2.r < min_ep2.r)
		{
			min_ep2.r = ep2.r;
		}

		if (ep2.g > 0.5f && ep2.g < min_ep2.g)
		{
			min_ep2.g = ep2.g;
		}

		if (ep2.b > 0.5f && ep2.b < min_ep2.b)
		{
			min_ep2.b = ep2.b;
		}

		if (ep2.a > 0.5f && ep2.a < min_ep2.a)
		{
			min_ep2.a = ep2.a;
		}

		#ifdef DEBUG_CAPTURE_NAN
			feenableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif
	}

	float min_wt_cutoff1, min_wt_cutoff2;
	switch (separate_component)
	{
	case 0:
		min_wt_cutoff2 = min_ep2.r;
		min_ep1.r = 1e30f;
		break;
	case 1:
		min_wt_cutoff2 = min_ep2.g;
		min_ep1.g = 1e30f;
		break;
	case 2:
		min_wt_cutoff2 = min_ep2.b;
		min_ep1.b = 1e30f;
		break;
	case 3:
		min_wt_cutoff2 = min_ep2.a;
		min_ep1.a = 1e30f;
		break;
	default:
		min_wt_cutoff2 = 1e30f;
	}

	min_wt_cutoff1 = MIN(MIN(min_ep1.r, min_ep1.g), MIN(min_ep1.b, min_ep1.a));

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
	for (int i = 0, ni = bsd->block_mode_packed_count; i < ni; ++i)
	{
		const block_mode& bm = bsd->block_modes_packed[i];
		if (bm.is_dual_plane != 1 || bm.percentile > mode_cutoff)
		{
			qwt_errors[i] = 1e38f;
			continue;
		}

		int decimation_mode = bm.decimation_mode;

		if (weight_high_value1[i] > 1.02f * min_wt_cutoff1)
		{
			weight_high_value1[i] = 1.0f;
		}

		if (weight_high_value2[i] > 1.02f * min_wt_cutoff2)
		{
			weight_high_value2[i] = 1.0f;
		}

		// compute weight bitcount for the mode
		int bits_used_by_weights = compute_ise_bitcount(2 * ixtab2[decimation_mode]->num_weights,
														(quantization_method) bm.quantization_mode);
		int bitcount = free_bits_for_partition_count[partition_count] - bits_used_by_weights;
		if (bitcount <= 0 || bits_used_by_weights < 24 || bits_used_by_weights > 96)
		{
			qwt_errors[i] = 1e38f;
			continue;
		}
		qwt_bitcounts[i] = bitcount;

		// then, generate the optimized set of weights for the mode.
		compute_ideal_quantized_weights_for_decimation_table(
		    ixtab2[decimation_mode],
		    weight_low_value1[i],
		    weight_high_value1[i],
		    decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * decimation_mode),
		    flt_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * i),
		    u8_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * i), bm.quantization_mode);

		compute_ideal_quantized_weights_for_decimation_table(
		    ixtab2[decimation_mode],
		    weight_low_value2[i],
		    weight_high_value2[i],
		    decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * decimation_mode + 1),
		    flt_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * i + 1),
		    u8_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * i + 1), bm.quantization_mode);


		// then, compute quantization errors for the block mode.
		qwt_errors[i] =	compute_error_of_weight_set(
		                    &(eix1[decimation_mode]),
		                    ixtab2[decimation_mode],
		                    flt_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * i))
		              + compute_error_of_weight_set(
		                    &(eix2[decimation_mode]),
		                    ixtab2[decimation_mode],
		                    flt_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * i + 1));
	}

	// decide the optimal combination of color endpoint encodings and weight encodings.
	int partition_format_specifiers[TUNE_MAX_TRIAL_CANDIDATES][4];
	int quantized_weight[TUNE_MAX_TRIAL_CANDIDATES];
	int color_quantization_level[TUNE_MAX_TRIAL_CANDIDATES];
	int color_quantization_level_mod[TUNE_MAX_TRIAL_CANDIDATES];

	endpoints epm;
	merge_endpoints(&(ei1->ep), &(ei2->ep), separate_component, &epm);

	determine_optimal_set_of_endpoint_formats_to_use(
	    bsd, pi, blk, ewb, &epm, separate_component, qwt_bitcounts, qwt_errors,
	    tune_candidate_limit, partition_format_specifiers, quantized_weight,
	    color_quantization_level, color_quantization_level_mod);

	for (int i = 0; i < tune_candidate_limit; i++)
	{
		const int qw_packed_index = quantized_weight[i];
		if (qw_packed_index < 0)
		{
			scb->error_block = 1;
			scb++;
			continue;
		}

		uint8_t *u8_weight1_src;
		uint8_t *u8_weight2_src;
		int weights_to_copy;

		assert(qw_packed_index >= 0 && qw_packed_index < bsd->block_mode_packed_count);
		const block_mode& qw_bm = bsd->block_modes_packed[qw_packed_index];

		int decimation_mode = qw_bm.decimation_mode;
		int weight_quantization_mode = qw_bm.quantization_mode;
		const decimation_table *it = ixtab2[decimation_mode];

		u8_weight1_src = u8_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * qw_packed_index);
		u8_weight2_src = u8_quantized_decimated_quantized_weights + MAX_WEIGHTS_PER_BLOCK * (2 * qw_packed_index + 1);

		weights_to_copy = it->num_weights;

		// recompute the ideal color endpoints before storing them.
		merge_endpoints(&(eix1[decimation_mode].ep), &(eix2[decimation_mode].ep), separate_component, &epm);

		float4 rgbs_colors[4];
		float4 rgbo_colors[4];

		for (int l = 0; l < max_refinement_iters; l++)
		{
			recompute_ideal_colors(
			    weight_quantization_mode, &epm, rgbs_colors, rgbo_colors,
			    u8_weight1_src, u8_weight2_src, separate_component, pi, it, blk, ewb);

			// store the colors for the block
			for (int j = 0; j < partition_count; j++)
			{
				scb->color_formats[j] = pack_color_endpoints(
				                            epm.endpt0[j], epm.endpt1[j],
				                            rgbs_colors[j], rgbo_colors[j],
				                            partition_format_specifiers[i][j],
				                            scb->color_values[j],
				                            color_quantization_level[i]);
			}
			scb->color_formats_matched = 0;

			if ((partition_count >= 2 && scb->color_formats[0] == scb->color_formats[1]
				&& color_quantization_level[i] != color_quantization_level_mod[i])
				&& (partition_count == 2 || (scb->color_formats[0] == scb->color_formats[2] && (partition_count == 3 || (scb->color_formats[0] == scb->color_formats[3])))))
			{
				int colorvals[4][12];
				int color_formats_mod[4] = { 0 };
				for (int j = 0; j < partition_count; j++)
				{
					color_formats_mod[j] = pack_color_endpoints(
					                           epm.endpt0[j], epm.endpt1[j],
					                           rgbs_colors[j], rgbo_colors[j],
					                           partition_format_specifiers[i][j],
					                           colorvals[j],
					                           color_quantization_level_mod[i]);
				}

				if (color_formats_mod[0] == color_formats_mod[1]
					&& (partition_count == 2 || (color_formats_mod[0] == color_formats_mod[2] && (partition_count == 3 || (color_formats_mod[0] == color_formats_mod[3])))))
				{
					scb->color_formats_matched = 1;
					for (int j = 0; j < 4; j++)
					{
						for (int k = 0; k < 12; k++)
						{
							scb->color_values[j][k] = colorvals[j][k];
						}
					}

					for (int j = 0; j < 4; j++)
					{
						scb->color_formats[j] = color_formats_mod[j];
					}
				}
			}

			// store header fields
			scb->partition_count = partition_count;
			scb->partition_index = partition_index;
			scb->color_quantization_level = scb->color_formats_matched ? color_quantization_level_mod[i] : color_quantization_level[i];
			scb->block_mode = qw_bm.mode_index;
			scb->plane2_color_component = separate_component;
			scb->error_block = 0;

			if (scb->color_quantization_level < 4)
			{
				scb->error_block = 1;	// should never happen, but cannot prove it impossible
			}

			int adjustments = realign_weights(
				decode_mode, bsd, blk, ewb, scb, u8_weight1_src, u8_weight2_src);

			if (adjustments == 0)
			{
				break;
			}
		}

		for (int j = 0; j < weights_to_copy; j++)
		{
			scb->plane1_weights[j] = u8_weight1_src[j];
			scb->plane2_weights[j] = u8_weight2_src[j];
		}

		scb++;
	}
}

void expand_deblock_weights(
	astcenc_context& ctx
) {
	unsigned int xdim = ctx.config.block_x;
	unsigned int ydim = ctx.config.block_y;
	unsigned int zdim = ctx.config.block_z;

	float centerpos_x = (xdim - 1) * 0.5f;
	float centerpos_y = (ydim - 1) * 0.5f;
	float centerpos_z = (zdim - 1) * 0.5f;
	float *bef = ctx.deblock_weights;

	for (unsigned int z = 0; z < zdim; z++)
	{
		for (unsigned int y = 0; y < ydim; y++)
		{
			for (unsigned int x = 0; x < xdim; x++)
			{
				float xdif = (x - centerpos_x) / xdim;
				float ydif = (y - centerpos_y) / ydim;
				float zdif = (z - centerpos_z) / zdim;

				float wdif = 0.36f;
				float dist = astc::sqrt(xdif * xdif + ydif * ydif + zdif * zdif + wdif * wdif);
				*bef = powf(dist, ctx.config.b_deblock_weight);
				bef++;
			}
		}
	}
}

// Function to set error weights for each color component for each texel in a block.
// Returns the sum of all the error values set.
static float prepare_error_weight_block(
	const astcenc_context& ctx,
	const astcenc_image& input_image,
	const block_size_descriptor* bsd,
	const imageblock* blk,
	error_weight_block* ewb
) {
	int idx = 0;
	int any_mean_stdev_weight =
		ctx.config.v_rgb_mean != 0.0f || ctx.config.v_rgb_stdev != 0.0f || \
		ctx.config.v_a_mean != 0.0f || ctx.config.v_a_stdev != 0.0f;

	float4 derv[MAX_TEXELS_PER_BLOCK];
	imageblock_initialize_deriv(blk, bsd->texel_count, derv);
	float4 color_weights = float4(ctx.config.cw_r_weight,
	                              ctx.config.cw_g_weight,
	                              ctx.config.cw_b_weight,
	                              ctx.config.cw_a_weight);

	ewb->contains_zeroweight_texels = 0;

	for (int z = 0; z < bsd->zdim; z++)
	{
		for (int y = 0; y < bsd->ydim; y++)
		{
			for (int x = 0; x < bsd->xdim; x++)
			{
				unsigned int xpos = x + blk->xpos;
				unsigned int ypos = y + blk->ypos;
				unsigned int zpos = z + blk->zpos;

				if (xpos >= input_image.dim_x || ypos >= input_image.dim_y || zpos >= input_image.dim_z)
				{
					float4 weights = float4(1e-11f);
					ewb->error_weights[idx] = weights;
					ewb->contains_zeroweight_texels = 1;
				}
				else
				{
					float4 error_weight = float4(ctx.config.v_rgb_base,
					                             ctx.config.v_rgb_base,
					                             ctx.config.v_rgb_base,
					                             ctx.config.v_a_base);

					int ydt = input_image.dim_x;
					int zdt = input_image.dim_x * input_image.dim_y;

					if (any_mean_stdev_weight)
					{
						float4 avg = ctx.input_averages[zpos * zdt + ypos * ydt + xpos];
						if (avg.r < 6e-5f)
							avg.r = 6e-5f;
						if (avg.g < 6e-5f)
							avg.g = 6e-5f;
						if (avg.b < 6e-5f)
							avg.b = 6e-5f;
						if (avg.a < 6e-5f)
							avg.a = 6e-5f;

						avg = avg * avg;

						float4 variance = ctx.input_variances[zpos * zdt + ypos * ydt + xpos];
						variance = variance * variance;

						float favg = (avg.r + avg.g + avg.b) * (1.0f / 3.0f);
						float fvar = (variance.r + variance.g + variance.b) * (1.0f / 3.0f);

						float mixing = ctx.config.v_rgba_mean_stdev_mix;
						avg.r = favg * mixing + avg.r * (1.0f - mixing);
						avg.g = favg * mixing + avg.g * (1.0f - mixing);
						avg.b = favg * mixing + avg.b * (1.0f - mixing);

						variance.r = fvar * mixing + variance.r * (1.0f - mixing);
						variance.g = fvar * mixing + variance.g * (1.0f - mixing);
						variance.b = fvar * mixing + variance.b * (1.0f - mixing);

						float4 stdev = float4(astc::sqrt(MAX(variance.r, 0.0f)),
						                      astc::sqrt(MAX(variance.g, 0.0f)),
						                      astc::sqrt(MAX(variance.b, 0.0f)),
						                      astc::sqrt(MAX(variance.a, 0.0f)));

						avg.r *= ctx.config.v_rgb_mean;
						avg.g *= ctx.config.v_rgb_mean;
						avg.b *= ctx.config.v_rgb_mean;
						avg.a *= ctx.config.v_a_mean;

						stdev.r *= ctx.config.v_rgb_stdev;
						stdev.g *= ctx.config.v_rgb_stdev;
						stdev.b *= ctx.config.v_rgb_stdev;
						stdev.a *= ctx.config.v_a_stdev;

						error_weight = error_weight + avg + stdev;

						error_weight = float4(1.0f / error_weight.r,
						                      1.0f / error_weight.g,
						                      1.0f / error_weight.b,
						                      1.0f / error_weight.a);
					}

					if (ctx.config.flags & ASTCENC_FLG_MAP_NORMAL)
					{
						// Convert from 0 to 1 to -1 to +1 range.
						float xN = ((blk->data_r[idx] * (1.0f / 65535.0f)) - 0.5f) * 2.0f;
						float yN = ((blk->data_a[idx] * (1.0f / 65535.0f)) - 0.5f) * 2.0f;

						float denom = 1.0f - xN * xN - yN * yN;
						if (denom < 0.1f)
							denom = 0.1f;
						denom = 1.0f / denom;
						error_weight.r *= 1.0f + xN * xN * denom;
						error_weight.a *= 1.0f + yN * yN * denom;
					}

					if (ctx.config.flags & ASTCENC_FLG_USE_ALPHA_WEIGHT)
					{
						float alpha_scale;
						if (ctx.config.a_scale_radius != 0)
						{
							alpha_scale = ctx.input_alpha_averages[zpos * zdt + ypos * ydt + xpos];
						}
						else
						{
							alpha_scale = blk->data_a[idx] * (1.0f / 65535.0f);
						}

						if (alpha_scale < 0.0001f)
						{
							alpha_scale = 0.0001f;
						}

						alpha_scale *= alpha_scale;
						error_weight.r *= alpha_scale;
						error_weight.g *= alpha_scale;
						error_weight.b *= alpha_scale;
					}

					error_weight = error_weight * color_weights;
					error_weight = error_weight * ctx.deblock_weights[idx];

					// when we loaded the block to begin with, we applied a transfer function
					// and computed the derivative of the transfer function. However, the
					// error-weight computation so far is based on the original color values,
					// not the transfer-function values. As such, we must multiply the
					// error weights by the derivative of the inverse of the transfer function,
					// which is equivalent to dividing by the derivative of the transfer
					// function.

					error_weight.r /= (derv[idx].r * derv[idx].r * 1e-10f);
					error_weight.g /= (derv[idx].g * derv[idx].g * 1e-10f);
					error_weight.b /= (derv[idx].b * derv[idx].b * 1e-10f);
					error_weight.a /= (derv[idx].a * derv[idx].a * 1e-10f);

					ewb->error_weights[idx] = error_weight;
					if (dot(error_weight, float4(1.0f, 1.0f, 1.0f, 1.0f)) < 1e-10f)
					{
						ewb->contains_zeroweight_texels = 1;
					}
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

		ewb->texel_weight_r[i] = ewb->error_weights[i].r;
		ewb->texel_weight_g[i] = ewb->error_weights[i].g;
		ewb->texel_weight_b[i] = ewb->error_weights[i].b;
		ewb->texel_weight_a[i] = ewb->error_weights[i].a;

		ewb->texel_weight_rg[i] = (ewb->error_weights[i].r + ewb->error_weights[i].g) * 0.5f;
		ewb->texel_weight_rb[i] = (ewb->error_weights[i].r + ewb->error_weights[i].b) * 0.5f;
		ewb->texel_weight_gb[i] = (ewb->error_weights[i].g + ewb->error_weights[i].b) * 0.5f;
		ewb->texel_weight_ra[i] = (ewb->error_weights[i].r + ewb->error_weights[i].a) * 0.5f;

		ewb->texel_weight_gba[i] = (ewb->error_weights[i].g + ewb->error_weights[i].b + ewb->error_weights[i].a) * 0.333333f;
		ewb->texel_weight_rba[i] = (ewb->error_weights[i].r + ewb->error_weights[i].b + ewb->error_weights[i].a) * 0.333333f;
		ewb->texel_weight_rga[i] = (ewb->error_weights[i].r + ewb->error_weights[i].g + ewb->error_weights[i].a) * 0.333333f;
		ewb->texel_weight_rgb[i] = (ewb->error_weights[i].r + ewb->error_weights[i].g + ewb->error_weights[i].b) * 0.333333f;

		ewb->texel_weight[i] = (ewb->error_weights[i].r + ewb->error_weights[i].g + ewb->error_weights[i].b + ewb->error_weights[i].a) * 0.25f;
	}

	return dot(error_weight_sum, float4(1.0f, 1.0f, 1.0f, 1.0f));
}

static float prepare_block_statistics(
	int texels_per_block,
	const imageblock * blk,
	const error_weight_block* ewb
) {
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

	for (int i = 0; i < texels_per_block; i++)
	{
		float weight = ewb->texel_weight[i];
		assert(weight >= 0.0f);
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
	return lowest_correlation;
}

void compress_block(
	const astcenc_context& ctx,
	const astcenc_image& input_image,
	const imageblock* blk,
	symbolic_compressed_block& scb,
	physical_compressed_block& pcb,
	compress_symbolic_block_buffers* tmpbuf)
{
	astcenc_profile decode_mode = ctx.config.profile;
	const block_size_descriptor* bsd = ctx.bsd;

	if (blk->red_min == blk->red_max && blk->green_min == blk->green_max && blk->blue_min == blk->blue_max && blk->alpha_min == blk->alpha_max)
	{
		// detected a constant-color block. Encode as FP16 if using HDR
		scb.error_block = 0;

		if ((decode_mode == ASTCENC_PRF_HDR) ||
		    (decode_mode == ASTCENC_PRF_HDR_RGB_LDR_A))
		{
			scb.block_mode = -1;
			scb.partition_count = 0;
			float4 orig_color = blk->origin_texel;
			scb.constant_color[0] = float_to_sf16(orig_color.r, SF_NEARESTEVEN);
			scb.constant_color[1] = float_to_sf16(orig_color.g, SF_NEARESTEVEN);
			scb.constant_color[2] = float_to_sf16(orig_color.b, SF_NEARESTEVEN);
			scb.constant_color[3] = float_to_sf16(orig_color.a, SF_NEARESTEVEN);
		}
		else
		{
			// Encode as UNORM16 if NOT using HDR.
			scb.block_mode = -2;
			scb.partition_count = 0;
			float4 orig_color = blk->origin_texel;
			float red   = orig_color.r;
			float green = orig_color.g;
			float blue  = orig_color.b;
			float alpha = orig_color.a;

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

			scb.constant_color[0] = astc::flt2int_rtn(red * 65535.0f);
			scb.constant_color[1] = astc::flt2int_rtn(green * 65535.0f);
			scb.constant_color[2] = astc::flt2int_rtn(blue * 65535.0f);
			scb.constant_color[3] = astc::flt2int_rtn(alpha * 65535.0f);
		}

		symbolic_to_physical(*bsd, scb, pcb);
		return;
	}

	error_weight_block *ewb = &tmpbuf->ewb;
	float error_weight_sum = prepare_error_weight_block(ctx, input_image, bsd, blk, ewb);

	symbolic_compressed_block *tempblocks = tmpbuf->tempblocks;

	float error_of_best_block = 1e20f;

	float best_errorvals_in_modes[13];
	for (int i = 0; i < 13; i++)
	{
		best_errorvals_in_modes[i] = 1e30f;
	}

	int uses_alpha = imageblock_uses_alpha(blk);

	float mode_cutoff = ctx.config.tune_block_mode_limit / 100.0f;

	// Trial using 1 plane of weights and 1 partition.

	// Most of the time we test it twice, first with a mode cutoff of 0 and
	// then with the specified mode cutoff. This causes an early-out that
	// speeds up encoding of easy blocks. However, this optimization is
	// disabled for 4x4 and 5x4 blocks where it nearly always slows down the
	// compression and slightly reduces image quality.

	float modecutoffs[2];
	float errorval_mult[2] = { 2.5, 1 };
	modecutoffs[0] = 0;
	modecutoffs[1] = mode_cutoff;

	float lowest_correl;
	float best_errorval_in_mode;

	int start_trial = bsd->texel_count < TUNE_MAX_TEXELS_MODE0_FASTPATH ? 1 : 0;
	for (int i = start_trial; i < 2; i++)
	{
		compress_symbolic_block_fixed_partition_1_plane(
		    decode_mode, modecutoffs[i],
		    ctx.config.tune_candidate_limit,
		    ctx.config.tune_refinement_limit,
		    bsd, 1, 0, blk, ewb, tempblocks, &tmpbuf->planes);

		best_errorval_in_mode = 1e30f;
		for (unsigned int j = 0; j < ctx.config.tune_candidate_limit; j++)
		{
			if (tempblocks[j].error_block)
			{
				continue;
			}

			float errorval = compute_symbolic_block_difference(decode_mode, bsd, tempblocks + j, blk, ewb);
			errorval *= errorval_mult[i];
			if (errorval < best_errorval_in_mode)
			{
				best_errorval_in_mode = errorval;
			}

			if (errorval < error_of_best_block)
			{
				error_of_best_block = errorval;
				scb = tempblocks[j];
			}
		}

		// Mode 0
		best_errorvals_in_modes[0] = best_errorval_in_mode;
		if ((error_of_best_block / error_weight_sum) < ctx.config.tune_db_limit)
		{
			goto END_OF_TESTS;
		}
	}

	lowest_correl = prepare_block_statistics(bsd->texel_count, blk, ewb);

	// next, test the four possible 1-partition, 2-planes modes
	for (int i = 0; i < 4; i++)
	{

		if (lowest_correl > ctx.config.tune_two_plane_early_out_limit)
		{
			continue;
		}

		if (blk->grayscale && i != 3)
		{
			continue;
		}

		if (!uses_alpha && i == 3)
		{
			continue;
		}

		compress_symbolic_block_fixed_partition_2_planes(
		    decode_mode, mode_cutoff,
		    ctx.config.tune_candidate_limit,
		    ctx.config.tune_refinement_limit,
		    bsd, 1,	// partition count
		    0,	// partition index
		    i,	// the color component to test a separate plane of weights for.
		    blk, ewb, tempblocks, &tmpbuf->planes);

		best_errorval_in_mode = 1e30f;
		for (unsigned int j = 0; j < ctx.config.tune_candidate_limit; j++)
		{
			if (tempblocks[j].error_block)
			{
				continue;
			}

			float errorval = compute_symbolic_block_difference(decode_mode, bsd, tempblocks + j, blk, ewb);
			if (errorval < best_errorval_in_mode)
			{
				best_errorval_in_mode = errorval;
			}

			if (errorval < error_of_best_block)
			{
				error_of_best_block = errorval;
				scb = tempblocks[j];
			}

			// Modes 1-4
			best_errorvals_in_modes[i + 1] = best_errorval_in_mode;
		}

		if ((error_of_best_block / error_weight_sum) < ctx.config.tune_db_limit)
		{
			goto END_OF_TESTS;
		}
	}

	// find best blocks for 2, 3 and 4 partitions
	for (int partition_count = 2; partition_count <= 4; partition_count++)
	{
		int partition_indices_1plane[2];
		int partition_index_2planes;

		find_best_partitionings(bsd, blk, ewb, partition_count,
		                        ctx.config.tune_partition_limit,
		                        &(partition_indices_1plane[0]),
		                        &(partition_indices_1plane[1]),
		                        &partition_index_2planes);

		for (int i = 0; i < 2; i++)
		{
			compress_symbolic_block_fixed_partition_1_plane(
			    decode_mode, mode_cutoff,
			    ctx.config.tune_candidate_limit,
			    ctx.config.tune_refinement_limit,
			    bsd, partition_count, partition_indices_1plane[i],
			    blk, ewb, tempblocks, &tmpbuf->planes);

			best_errorval_in_mode = 1e30f;
			for (unsigned int j = 0; j < ctx.config.tune_candidate_limit; j++)
			{
				if (tempblocks[j].error_block)
				{
					continue;
				}

				float errorval = compute_symbolic_block_difference(decode_mode, bsd, tempblocks + j, blk, ewb);
				if (errorval < best_errorval_in_mode)
				{
					best_errorval_in_mode = errorval;
				}

				if (errorval < error_of_best_block)
				{
					error_of_best_block = errorval;
					scb = tempblocks[j];
				}
			}

			// Modes 5, 6, 8, 9, 11, 12
			best_errorvals_in_modes[3 * (partition_count - 2) + 5 + i] = best_errorval_in_mode;

			if ((error_of_best_block / error_weight_sum) < ctx.config.tune_db_limit)
			{
				goto END_OF_TESTS;
			}
		}

		if (partition_count == 2 && MIN(best_errorvals_in_modes[5], best_errorvals_in_modes[6]) > (best_errorvals_in_modes[0] * ctx.config.tune_partition_early_out_limit))
		{
			goto END_OF_TESTS;
		}

		// Skip testing dual weight planes for:
		// * 4 partitions (can't be encoded by the format)
		// * Luminance only blocks (never need for a second plane)
		// * Blocks with higher component correlation than the tuning cutoff
		if ((partition_count == 4) ||
		    (blk->grayscale && !uses_alpha) ||
		    (lowest_correl > ctx.config.tune_two_plane_early_out_limit))
		{
			continue;
		}


		if (lowest_correl <= ctx.config.tune_two_plane_early_out_limit)
		{
			compress_symbolic_block_fixed_partition_2_planes(
			    decode_mode,
			    mode_cutoff,
			    ctx.config.tune_candidate_limit,
			    ctx.config.tune_refinement_limit,
			    bsd,
			    partition_count,
			    partition_index_2planes & (PARTITION_COUNT - 1),
			    partition_index_2planes >> PARTITION_BITS,
			    blk, ewb, tempblocks, &tmpbuf->planes);

			best_errorval_in_mode = 1e30f;
			for (unsigned int j = 0; j < ctx.config.tune_candidate_limit; j++)
			{
				if (tempblocks[j].error_block)
				{
					continue;
				}

				float errorval = compute_symbolic_block_difference(decode_mode, bsd, tempblocks + j, blk, ewb);
				if (errorval < best_errorval_in_mode)
				{
					best_errorval_in_mode = errorval;
				}

				if (errorval < error_of_best_block)
				{
					error_of_best_block = errorval;
					scb = tempblocks[j];
				}
			}

			// Modes 7, 10 (13 is unreachable)
			best_errorvals_in_modes[3 * (partition_count - 2) + 5 + 2] = best_errorval_in_mode;

			if ((error_of_best_block / error_weight_sum) < ctx.config.tune_db_limit)
			{
				goto END_OF_TESTS;
			}
		}
	}

END_OF_TESTS:
	// compress/decompress to a physical block
	symbolic_to_physical(*bsd, scb, pcb);
}

#endif
