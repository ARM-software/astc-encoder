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
 * @brief Functions for finding dominant direction of a set of colors.
 */
#if !defined(ASTCENC_DECOMPRESS_ONLY)

#include "astcenc_internal.h"

#include <cassert>

#ifdef DEBUG_CAPTURE_NAN
	#ifndef _GNU_SOURCE
		#define _GNU_SOURCE
	#endif

	#include <fenv.h>
#endif

// For a full block, functions to compute averages and dominant directions. The
// averages and directions are computed separately for each partition.
// We have separate versions for blocks with and without alpha, since the
// processing for blocks with alpha is significantly more expensive. The
// direction vectors it produces are NOT normalized.
void compute_averages_and_directions_rgba(
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	const float4* color_scalefactors,
	float4* averages,
	float4* directions_rgba
) {
	int partition_count = pt->partition_count;
	for (int partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *weights = pt->texels_of_partition[partition];
		int texelcount = pt->texels_per_partition[partition];

		float4 base_sum = float4(0.0f, 0.0f, 0.0f, 0.0f);
		float partition_weight = 0.0f;

		for (int i = 0; i < texelcount; i++)
		{
			int iwt = weights[i];
			float weight = ewb->texel_weight[iwt];
			float4 texel_datum = float4(blk->data_r[iwt],
			                            blk->data_g[iwt],
			                            blk->data_b[iwt],
			                            blk->data_a[iwt]) * weight;
			partition_weight += weight;

			base_sum = base_sum + texel_datum;
		}

		float4 average = base_sum * (1.0f / MAX(partition_weight, 1e-7f));
		averages[partition] = average * color_scalefactors[partition];

		float4 sum_xp = float4(0.0f, 0.0f, 0.0f, 0.0f);
		float4 sum_yp = float4(0.0f, 0.0f, 0.0f, 0.0f);
		float4 sum_zp = float4(0.0f, 0.0f, 0.0f, 0.0f);
		float4 sum_wp = float4(0.0f, 0.0f, 0.0f, 0.0f);

		for (int i = 0; i < texelcount; i++)
		{
			int iwt = weights[i];
			float weight = ewb->texel_weight[iwt];
			float4 texel_datum = float4(blk->data_r[iwt],
			                            blk->data_g[iwt],
			                            blk->data_b[iwt],
			                            blk->data_a[iwt]);
			texel_datum = (texel_datum - average) * weight;

			if (texel_datum.x > 0.0f)
			{
				sum_xp = sum_xp + texel_datum;
			}

			if (texel_datum.y > 0.0f)
			{
				sum_yp = sum_yp + texel_datum;
			}

			if (texel_datum.z > 0.0f)
			{
				sum_zp = sum_zp + texel_datum;
			}

			if (texel_datum.w > 0.0f)
			{
				sum_wp = sum_wp + texel_datum;
			}
		}

		float prod_xp = dot(sum_xp, sum_xp);
		float prod_yp = dot(sum_yp, sum_yp);
		float prod_zp = dot(sum_zp, sum_zp);
		float prod_wp = dot(sum_wp, sum_wp);

		float4 best_vector = sum_xp;
		float best_sum = prod_xp;

		if (prod_yp > best_sum)
		{
			best_vector = sum_yp;
			best_sum = prod_yp;
		}

		if (prod_zp > best_sum)
		{
			best_vector = sum_zp;
			best_sum = prod_zp;
		}

		if (prod_wp > best_sum)
		{
			best_vector = sum_wp;
		}

		directions_rgba[partition] = best_vector;
	}
}

void compute_averages_and_directions_rgb(
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	const float4* color_scalefactors,
	float3* averages,
	float3* directions_rgb
) {
	int partition_count = pt->partition_count;
	const float *texel_weights = ewb->texel_weight_rgb;

	for (int partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *weights = pt->texels_of_partition[partition];
		int texelcount = pt->texels_per_partition[partition];

		float3 base_sum = float3(0.0f, 0.0f, 0.0f);
		float partition_weight = 0.0f;

		for (int i = 0; i < texelcount; i++)
		{
			int iwt = weights[i];
			float weight = texel_weights[iwt];
			float3 texel_datum = float3(blk->data_r[iwt],
			                            blk->data_g[iwt],
			                            blk->data_b[iwt]) * weight;
			partition_weight += weight;

			base_sum = base_sum + texel_datum;
		}

		float4 csf = color_scalefactors[partition];
		float3 average = base_sum * (1.0f / MAX(partition_weight, 1e-7f));
		averages[partition] = average * float3(csf.x, csf.y, csf.z);

		float3 sum_xp = float3(0.0f, 0.0f, 0.0f);
		float3 sum_yp = float3(0.0f, 0.0f, 0.0f);
		float3 sum_zp = float3(0.0f, 0.0f, 0.0f);

		for (int i = 0; i < texelcount; i++)
		{
			int iwt = weights[i];
			float weight = texel_weights[iwt];
			float3 texel_datum = float3(blk->data_r[iwt],
			                            blk->data_g[iwt],
			                            blk->data_b[iwt]);
			texel_datum = (texel_datum - average) * weight;

			if (texel_datum.x > 0.0f)
			{
				sum_xp = sum_xp + texel_datum;
			}

			if (texel_datum.y > 0.0f)
			{
				sum_yp = sum_yp + texel_datum;
			}

			if (texel_datum.z > 0.0f)
			{
				sum_zp = sum_zp + texel_datum;
			}
		}

		float prod_xp = dot(sum_xp, sum_xp);
		float prod_yp = dot(sum_yp, sum_yp);
		float prod_zp = dot(sum_zp, sum_zp);

		float3 best_vector = sum_xp;
		float best_sum = prod_xp;

		if (prod_yp > best_sum)
		{
			best_vector = sum_yp;
			best_sum = prod_yp;
		}

		if (prod_zp > best_sum)
		{
			best_vector = sum_zp;
		}

		directions_rgb[partition] = best_vector;
	}
}

void compute_averages_and_directions_3_components(
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	const float3* color_scalefactors,
	int omitted_component,
	float3* averages,
	float3* directions
) {
	const float *texel_weights;
	const float* data_vr;
	const float* data_vg;
	const float* data_vb;

	if (omitted_component == 0)
	{
		texel_weights = ewb->texel_weight_gba;
		data_vr = blk->data_g;
		data_vg = blk->data_b;
		data_vb = blk->data_a;
	}
	else if (omitted_component == 1)
	{
		texel_weights = ewb->texel_weight_rba;
		data_vr = blk->data_r;
		data_vg = blk->data_b;
		data_vb = blk->data_a;
	}
	else if (omitted_component == 2)
	{
		texel_weights = ewb->texel_weight_rga;
		data_vr = blk->data_r;
		data_vg = blk->data_g;
		data_vb = blk->data_a;
	}
	else
	{
		assert(omitted_component == 3);
		texel_weights = ewb->texel_weight_rgb;
		data_vr = blk->data_r;
		data_vg = blk->data_g;
		data_vb = blk->data_b;
	}

	int partition_count = pt->partition_count;
	for (int partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *weights = pt->texels_of_partition[partition];
		int texelcount = pt->texels_per_partition[partition];

		float3 base_sum = float3(0.0f, 0.0f, 0.0f);
		float partition_weight = 0.0f;

		for (int i = 0; i < texelcount; i++)
		{
			int iwt = weights[i];
			float weight = texel_weights[iwt];
			float3 texel_datum = float3(data_vr[iwt],
			                            data_vg[iwt],
			                            data_vb[iwt]) * weight;
			partition_weight += weight;

			base_sum = base_sum + texel_datum;
		}

		float3 csf = color_scalefactors[partition];

		float3 average = base_sum * (1.0f / MAX(partition_weight, 1e-7f));
		averages[partition] = average * float3(csf.x, csf.y, csf.z);

		float3 sum_xp = float3(0.0f, 0.0f, 0.0f);
		float3 sum_yp = float3(0.0f, 0.0f, 0.0f);
		float3 sum_zp = float3(0.0f, 0.0f, 0.0f);

		for (int i = 0; i < texelcount; i++)
		{
			int iwt = weights[i];
			float weight = texel_weights[iwt];
			float3 texel_datum = float3(data_vr[iwt],
			                            data_vg[iwt],
			                            data_vb[iwt]);
			texel_datum = (texel_datum - average) * weight;

			if (texel_datum.x > 0.0f)
			{
				sum_xp = sum_xp + texel_datum;
			}

			if (texel_datum.y > 0.0f)
			{
				sum_yp = sum_yp + texel_datum;
			}

			if (texel_datum.z > 0.0f)
			{
				sum_zp = sum_zp + texel_datum;
			}
		}

		float prod_xp = dot(sum_xp, sum_xp);
		float prod_yp = dot(sum_yp, sum_yp);
		float prod_zp = dot(sum_zp, sum_zp);

		float3 best_vector = sum_xp;
		float best_sum = prod_xp;

		if (prod_yp > best_sum)
		{
			best_vector = sum_yp;
			best_sum = prod_yp;
		}

		if (prod_zp > best_sum)
		{
			best_vector = sum_zp;
		}

		if (dot(best_vector, best_vector) < 1e-18f)
		{
			best_vector = float3(1.0f, 1.0f, 1.0f);
		}

		directions[partition] = best_vector;
	}

}

void compute_averages_and_directions_2_components(
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	const float2* color_scalefactors,
	int component1,
	int component2,
	float2* averages,
	float2* directions
) {
	const float *texel_weights;
	const float* data_vr = nullptr;
	const float* data_vg = nullptr;

	if (component1 == 0 && component2 == 1)
	{
		texel_weights = ewb->texel_weight_rg;
		data_vr = blk->data_r;
		data_vg = blk->data_g;
	}
	else if (component1 == 0 && component2 == 2)
	{
		texel_weights = ewb->texel_weight_rb;
		data_vr = blk->data_r;
		data_vg = blk->data_b;
	}
	else // (component1 == 1 && component2 == 2)
	{
		assert(component1 == 1 && component2 == 2);
		texel_weights = ewb->texel_weight_gb;
		data_vr = blk->data_g;
		data_vg = blk->data_b;
	}

	int partition_count = pt->partition_count;
	for (int partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *weights = pt->texels_of_partition[partition];
		int texelcount = pt->texels_per_partition[partition];

		float2 base_sum = float2(0.0f, 0.0f);
		float partition_weight = 0.0f;

		for (int i = 0; i < texelcount; i++)
		{
			int iwt = weights[i];
			float weight = texel_weights[iwt];
			float2 texel_datum = float2(data_vr[iwt], data_vg[iwt]) * weight;
			partition_weight += weight;

			base_sum = base_sum + texel_datum;
		}

		float2 csf = color_scalefactors[partition];

		float2 average = base_sum * (1.0f / MAX(partition_weight, 1e-7f));
		averages[partition] = average * float2(csf.x, csf.y);

		float2 sum_xp = float2(0.0f, 0.0f);
		float2 sum_yp = float2(0.0f, 0.0f);

		for (int i = 0; i < texelcount; i++)
		{
			int iwt = weights[i];
			float weight = texel_weights[iwt];
			float2 texel_datum = float2(data_vr[iwt], data_vg[iwt]);
			texel_datum = (texel_datum - average) * weight;

			if (texel_datum.x > 0.0f)
			{
				sum_xp = sum_xp + texel_datum;
			}

			if (texel_datum.y > 0.0f)
			{
				sum_yp = sum_yp + texel_datum;
			}
		}

		float prod_xp = dot(sum_xp, sum_xp);
		float prod_yp = dot(sum_yp, sum_yp);

		float2 best_vector = sum_xp;
		float best_sum = prod_xp;

		if (prod_yp > best_sum)
		{
			best_vector = sum_yp;
		}

		directions[partition] = best_vector;
	}
}

void compute_error_squared_rgba(
	const partition_info* pt,    // the partition that we use when computing the squared-error.
	const imageblock* blk,
	const error_weight_block* ewb,
	const processed_line4* plines_uncorr,
	const processed_line4* plines_samechroma,
	const processed_line3* plines_separate_red,
	const processed_line3* plines_separate_green,
	const processed_line3* plines_separate_blue,
	const processed_line3* plines_separate_alpha,
	float* lengths_uncorr,
	float* lengths_samechroma,
	float4* lengths_separate,
	float* uncorr_errors,
	float* samechroma_errors,
	float4* separate_color_errors
) {
	float uncorr_errorsum = 0.0f;
	float samechroma_errorsum = 0.0f;
	float red_errorsum = 0.0f;
	float green_errorsum = 0.0f;
	float blue_errorsum = 0.0f;
	float alpha_errorsum = 0.0f;

	for (int partition = 0; partition < pt->partition_count; partition++)
	{
		// TODO: sort partitions by number of texels. For warp-architectures,
		// this can reduce the running time by about 25-50%.
		const uint8_t *weights = pt->texels_of_partition[partition];
		int texelcount = pt->texels_per_partition[partition];

		float uncorr_lowparam = 1e10f;
		float uncorr_highparam = -1e10f;

		float samechroma_lowparam = 1e10f;
		float samechroma_highparam = -1e10f;

		float4 separate_lowparam = float4(1e10f, 1e10f, 1e10f, 1e10f);
		float4 separate_highparam = float4(-1e10f, -1e10f, -1e10f, -1e10f);

		processed_line4 l_uncorr = plines_uncorr[partition];
		processed_line4 l_samechroma = plines_samechroma[partition];
		processed_line3 l_red = plines_separate_red[partition];
		processed_line3 l_green = plines_separate_green[partition];
		processed_line3 l_blue = plines_separate_blue[partition];
		processed_line3 l_alpha = plines_separate_alpha[partition];

		// TODO: split up this loop due to too many temporaries; in particular,
		// the six line functions will consume 18 vector registers
		for (int i = 0; i < texelcount; i++)
		{
			int iwt = weights[i];

			float texel_weight_rgba = ewb->texel_weight[iwt];
			if (texel_weight_rgba > 1e-20f)
			{
				float4 dat = float4(blk->data_r[iwt],
				                    blk->data_g[iwt],
				                    blk->data_b[iwt],
				                    blk->data_a[iwt]);

				float4 ews = ewb->error_weights[iwt];

				float uncorr_param = dot(dat, l_uncorr.bs);
				uncorr_lowparam = MIN(uncorr_param, uncorr_lowparam);
				uncorr_highparam = MAX(uncorr_param, uncorr_highparam);

				float samechroma_param = dot(dat, l_samechroma.bs);
				samechroma_lowparam = MIN(samechroma_param, samechroma_lowparam);
				samechroma_highparam = MAX(samechroma_param, samechroma_highparam);

				float4 separate_param = float4(dot(float3(dat.y, dat.z, dat.w), l_red.bs),
				                               dot(float3(dat.x, dat.z, dat.w), l_green.bs),
				                               dot(float3(dat.x, dat.y, dat.w), l_blue.bs),
				                               dot(float3(dat.x, dat.y, dat.z), l_alpha.bs));

				separate_lowparam = float4(MIN(separate_param.x, separate_lowparam.x),
				                           MIN(separate_param.y, separate_lowparam.y),
				                           MIN(separate_param.z, separate_lowparam.z),
				                           MIN(separate_param.w, separate_lowparam.w));

				separate_highparam = float4(MAX(separate_param.x, separate_highparam.x),
				                            MAX(separate_param.y, separate_highparam.y),
				                            MAX(separate_param.z, separate_highparam.z),
				                            MAX(separate_param.w, separate_highparam.w));

				float4 uncorr_dist  = (l_uncorr.amod - dat) + (uncorr_param * l_uncorr.bis);
				uncorr_errorsum += dot(ews, uncorr_dist * uncorr_dist);

				float4 samechroma_dist = (l_samechroma.amod - dat) +
				                         (samechroma_param * l_samechroma.bis);
				samechroma_errorsum += dot(ews, samechroma_dist * samechroma_dist);

				float3 red_dist = (l_red.amod - float3(dat.y, dat.z, dat.w)) +
				                  (separate_param.x * l_red.bis);
				red_errorsum += dot(float3(ews.y, ews.z, ews.w), red_dist * red_dist);

				float3 green_dist  = (l_green.amod - float3(dat.x, dat.z, dat.w)) +
				                     (separate_param.y * l_green.bis);
				green_errorsum += dot(float3(ews.x, ews.z, ews.w), green_dist * green_dist);

				float3 blue_dist  = (l_blue.amod - float3(dat.x, dat.y, dat.w)) +
				                    (separate_param.z * l_blue.bis);
				blue_errorsum += dot(float3(ews.x, ews.y, ews.w), blue_dist * blue_dist);

				float3 alpha_dist  = (l_alpha.amod - float3(dat.x, dat.y, dat.z)) +
				                     (separate_param.w * l_alpha.bis);
				alpha_errorsum += dot(float3(ews.x, ews.y, ews.z), alpha_dist * alpha_dist);
			}
		}

		float uncorr_linelen = uncorr_highparam - uncorr_lowparam;
		float samechroma_linelen = samechroma_highparam - samechroma_lowparam;
		float4 separate_linelen = separate_highparam - separate_lowparam;

		// Turn very small numbers and NaNs into a small number
		if (!(uncorr_linelen > 1e-7f))
		{
			uncorr_linelen = 1e-7f;
		}

		if (!(samechroma_linelen > 1e-7f))
		{
			samechroma_linelen = 1e-7f;
		}

		if (!(separate_linelen.x > 1e-7f))
		{
			separate_linelen.x = 1e-7f;
		}

		if (!(separate_linelen.y > 1e-7f))
		{
			separate_linelen.y = 1e-7f;
		}

		if (!(separate_linelen.z > 1e-7f))
		{
			separate_linelen.z = 1e-7f;
		}

		if (!(separate_linelen.w > 1e-7f))
		{
			separate_linelen.w = 1e-7f;
		}

		lengths_uncorr[partition] = uncorr_linelen;
		lengths_samechroma[partition] = samechroma_linelen;
		lengths_separate[partition] = separate_linelen;

		*uncorr_errors = uncorr_errorsum;
		*samechroma_errors = samechroma_errorsum;
		*separate_color_errors = float4(red_errorsum, green_errorsum, blue_errorsum, alpha_errorsum);
	}
}

void compute_error_squared_rgb(
	const partition_info *pt,    // the partition that we use when computing the squared-error.
	const imageblock *blk,
	const error_weight_block *ewb,
	const processed_line3 *plines_uncorr,
	const processed_line3 *plines_samechroma,
	const processed_line2 *plines_separate_red,
	const processed_line2 *plines_separate_green,
	const processed_line2 *plines_separate_blue,
	float *lengths_uncorr,
	float *lengths_samechroma,
	float3 *lengths_separate,
	float *uncorr_errors,
	float *samechroma_errors,
	float3 *separate_color_errors
) {
	float uncorr_errorsum = 0.0f;
	float samechroma_errorsum = 0.0f;
	float red_errorsum = 0.0f;
	float green_errorsum = 0.0f;
	float blue_errorsum = 0.0f;

	for (int partition = 0; partition < pt->partition_count; partition++)
	{
		// TODO: sort partitions by number of texels. For warp-architectures,
		// this can reduce the running time by about 25-50%.
		const uint8_t *weights = pt->texels_of_partition[partition];
		int texelcount = pt->texels_per_partition[partition];

		float uncorr_lowparam = 1e10f;
		float uncorr_highparam = -1e10f;

		float samechroma_lowparam = 1e10f;
		float samechroma_highparam = -1e10f;

		float3 separate_lowparam = float3(1e10f, 1e10f, 1e10f);
		float3 separate_highparam = float3(-1e10f, -1e10f, -1e10f);

		processed_line3 l_uncorr = plines_uncorr[partition];
		processed_line3 l_samechroma = plines_samechroma[partition];
		processed_line2 l_red = plines_separate_red[partition];
		processed_line2 l_green = plines_separate_green[partition];
		processed_line2 l_blue = plines_separate_blue[partition];

		// TODO: split up this loop due to too many temporaries; in
		// particular, the six line functions will consume 18 vector registers

		for (int i = 0; i < texelcount; i++)
		{
			int iwt = weights[i];

			float texel_weight_rgb = ewb->texel_weight_rgb[iwt];
			if (texel_weight_rgb > 1e-20f)
			{
				float3 dat = float3(blk->data_r[iwt],
				                    blk->data_g[iwt],
				                    blk->data_b[iwt]);

				float3 ews = float3(ewb->error_weights[iwt].x,
				                    ewb->error_weights[iwt].y,
				                    ewb->error_weights[iwt].z);

				float uncorr_param = dot(dat, l_uncorr.bs);
				uncorr_lowparam  = MIN(uncorr_param, uncorr_lowparam);
				uncorr_highparam = MAX(uncorr_param, uncorr_highparam);

				float samechroma_param = dot(dat, l_samechroma.bs);
				samechroma_lowparam  = MIN(samechroma_param, samechroma_lowparam);
				samechroma_highparam = MAX(samechroma_param, samechroma_highparam);

				float3 separate_param = float3(dot(float2(dat.y, dat.z), l_red.bs),
				                               dot(float2(dat.x, dat.z), l_green.bs),
				                               dot(float2(dat.x, dat.y), l_blue.bs));

				separate_lowparam  = float3(MIN(separate_param.x, separate_lowparam.x),
				                            MIN(separate_param.y, separate_lowparam.y),
				                            MIN(separate_param.z, separate_lowparam.z));

				separate_highparam  = float3(MAX(separate_param.x, separate_highparam.x),
				                             MAX(separate_param.y, separate_highparam.y),
				                             MAX(separate_param.z, separate_highparam.z));

				float3 uncorr_dist  = (l_uncorr.amod - dat) +
				                      (uncorr_param * l_uncorr.bis);
				uncorr_errorsum += dot(ews, uncorr_dist * uncorr_dist);

				float3 samechroma_dist = (l_samechroma.amod - dat) +
				                         (samechroma_param * l_samechroma.bis);
				samechroma_errorsum += dot(ews, samechroma_dist * samechroma_dist);

				float2 red_dist = (l_red.amod - float2(dat.y, dat.z)) +
				                  (separate_param.x * l_red.bis);
				red_errorsum += dot(float2(ews.y, ews.z), red_dist * red_dist);

				float2 green_dist = (l_green.amod - float2(dat.x, dat.z)) +
				                    (separate_param.y * l_green.bis);
				green_errorsum += dot(float2(ews.x, ews.z), green_dist * green_dist);

				float2 blue_dist = (l_blue.amod - float2(dat.x, dat.y)) +
				                   (separate_param.z * l_blue.bis);
				blue_errorsum += dot(float2(ews.x, ews.y), blue_dist * blue_dist);
			}
		}

		float uncorr_linelen = uncorr_highparam - uncorr_lowparam;
		float samechroma_linelen = samechroma_highparam - samechroma_lowparam;
		float3 separate_linelen = separate_highparam - separate_lowparam;

		// Turn very small numbers and NaNs into a small number
		if (!(uncorr_linelen > 1e-7f))
		{
			uncorr_linelen = 1e-7f;
		}

		if (!(samechroma_linelen > 1e-7f))
		{
			samechroma_linelen = 1e-7f;
		}

		if (!(separate_linelen.x > 1e-7f))
		{
			separate_linelen.x = 1e-7f;
		}

		if (!(separate_linelen.y > 1e-7f))
		{
			separate_linelen.y = 1e-7f;
		}

		if (!(separate_linelen.z > 1e-7f))
		{
			separate_linelen.z = 1e-7f;
		}

		lengths_uncorr[partition] = uncorr_linelen;
		lengths_samechroma[partition] = samechroma_linelen;
		lengths_separate[partition] = separate_linelen;

		*uncorr_errors = uncorr_errorsum;
		*samechroma_errors = samechroma_errorsum;
		*separate_color_errors = float3(red_errorsum, green_errorsum, blue_errorsum);
	}
}

// function to compute the error across a tile when using a particular line for
// a particular partition.
float compute_error_squared_rgb_single_partition(
	int partition_to_test,
	const block_size_descriptor* bsd,
	const partition_info* pt,	// the partition that we use when computing the squared-error.
	const imageblock* blk,
	const error_weight_block* ewb,
	const processed_line3* lin	// the line for the partition.
) {
	int texels_per_block = bsd->texel_count;
	float errorsum = 0.0f;

	for (int i = 0; i < texels_per_block; i++)
	{
		int partition = pt->partition_of_texel[i];
		float texel_weight = ewb->texel_weight_rgb[i];

		if (partition != partition_to_test || texel_weight < 1e-20f)
		{
			continue;
		}

		float3 point = float3(blk->data_r[i],
		                      blk->data_g[i],
		                      blk->data_b[i]);
		float param = dot(point, lin->bs);
		float3 rp1 = lin->amod + param * lin->bis;
		float3 dist = rp1 - point;
		float4 ews = ewb->error_weights[i];
		float3 ews3 = float3(ews.x, ews.y, ews.z);
		errorsum += dot(ews3, dist * dist);
	}

	return errorsum;
}

#endif
