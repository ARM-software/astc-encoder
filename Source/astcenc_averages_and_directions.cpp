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
	promise(partition_count > 0);

	for (int partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *weights = pt->texels_of_partition[partition];

		float4 base_sum = float4(0.0f);
		float partition_weight = 0.0f;

		int texel_count = pt->texels_per_partition[partition];
		promise(texel_count > 0);

		for (int i = 0; i < texel_count; i++)
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

		float4 average = base_sum * (1.0f / astc::max(partition_weight, 1e-7f));
		averages[partition] = average * color_scalefactors[partition];

		float4 sum_xp = float4(0.0f);
		float4 sum_yp = float4(0.0f);
		float4 sum_zp = float4(0.0f);
		float4 sum_wp = float4(0.0f);

		for (int i = 0; i < texel_count; i++)
		{
			int iwt = weights[i];
			float weight = ewb->texel_weight[iwt];
			float4 texel_datum = float4(blk->data_r[iwt],
			                            blk->data_g[iwt],
			                            blk->data_b[iwt],
			                            blk->data_a[iwt]);
			texel_datum = (texel_datum - average) * weight;

			if (texel_datum.r > 0.0f)
			{
				sum_xp = sum_xp + texel_datum;
			}

			if (texel_datum.g > 0.0f)
			{
				sum_yp = sum_yp + texel_datum;
			}

			if (texel_datum.b > 0.0f)
			{
				sum_zp = sum_zp + texel_datum;
			}

			if (texel_datum.a > 0.0f)
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
	const float *texel_weights = ewb->texel_weight_rgb;

	int partition_count = pt->partition_count;
	promise(partition_count > 0);

	for (int partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *weights = pt->texels_of_partition[partition];

		float3 base_sum = float3(0.0f, 0.0f, 0.0f);
		float partition_weight = 0.0f;

		int texel_count = pt->texels_per_partition[partition];
		promise(texel_count > 0);

		for (int i = 0; i < texel_count; i++)
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
		float3 average = base_sum * (1.0f / astc::max(partition_weight, 1e-7f));
		averages[partition] = average * float3(csf.r, csf.g, csf.b);

		float3 sum_xp = float3(0.0f);
		float3 sum_yp = float3(0.0f);
		float3 sum_zp = float3(0.0f);

		for (int i = 0; i < texel_count; i++)
		{
			int iwt = weights[i];
			float weight = texel_weights[iwt];
			float3 texel_datum = float3(blk->data_r[iwt],
			                            blk->data_g[iwt],
			                            blk->data_b[iwt]);
			texel_datum = (texel_datum - average) * weight;

			if (texel_datum.r > 0.0f)
			{
				sum_xp = sum_xp + texel_datum;
			}

			if (texel_datum.g > 0.0f)
			{
				sum_yp = sum_yp + texel_datum;
			}

			if (texel_datum.b > 0.0f)
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
	promise(partition_count > 0);

	for (int partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *weights = pt->texels_of_partition[partition];

		float3 base_sum = float3(0.0f);
		float partition_weight = 0.0f;

		int texel_count = pt->texels_per_partition[partition];
		promise(texel_count > 0);

		for (int i = 0; i < texel_count; i++)
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

		float3 average = base_sum * (1.0f / astc::max(partition_weight, 1e-7f));
		averages[partition] = average * float3(csf.r, csf.g, csf.b);

		float3 sum_xp = float3(0.0f);
		float3 sum_yp = float3(0.0f);
		float3 sum_zp = float3(0.0f);

		for (int i = 0; i < texel_count; i++)
		{
			int iwt = weights[i];
			float weight = texel_weights[iwt];
			float3 texel_datum = float3(data_vr[iwt],
			                            data_vg[iwt],
			                            data_vb[iwt]);
			texel_datum = (texel_datum - average) * weight;

			if (texel_datum.r > 0.0f)
			{
				sum_xp = sum_xp + texel_datum;
			}

			if (texel_datum.g > 0.0f)
			{
				sum_yp = sum_yp + texel_datum;
			}

			if (texel_datum.b > 0.0f)
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
	promise(partition_count > 0);

	for (int partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *weights = pt->texels_of_partition[partition];

		float2 base_sum = float2(0.0f);
		float partition_weight = 0.0f;

		int texel_count = pt->texels_per_partition[partition];
		promise(texel_count > 0);

		for (int i = 0; i < texel_count; i++)
		{
			int iwt = weights[i];
			float weight = texel_weights[iwt];
			float2 texel_datum = float2(data_vr[iwt], data_vg[iwt]) * weight;
			partition_weight += weight;

			base_sum = base_sum + texel_datum;
		}

		float2 csf = color_scalefactors[partition];

		float2 average = base_sum * (1.0f / astc::max(partition_weight, 1e-7f));
		averages[partition] = average * float2(csf.r, csf.g);

		float2 sum_xp = float2(0.0f);
		float2 sum_yp = float2(0.0f);

		for (int i = 0; i < texel_count; i++)
		{
			int iwt = weights[i];
			float weight = texel_weights[iwt];
			float2 texel_datum = float2(data_vr[iwt], data_vg[iwt]);
			texel_datum = (texel_datum - average) * weight;

			if (texel_datum.r > 0.0f)
			{
				sum_xp = sum_xp + texel_datum;
			}

			if (texel_datum.g > 0.0f)
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
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	const processed_line4* uncor_plines,
	const processed_line4* samec_plines,
	float* uncor_lengths,
	float* samec_lengths,
	float* uncor_errors,
	float* samec_errors
) {
	float uncor_errorsum = 0.0f;
	float samec_errorsum = 0.0f;

	int partition_count = pt->partition_count;
	promise(partition_count > 0);

	for (int partition = 0; partition < partition_count; partition++)
	{
		// TODO: sort partitions by number of texels. For warp-architectures,
		// this can reduce the running time by about 25-50%.
		const uint8_t *weights = pt->texels_of_partition[partition];

		float uncor_loparam = 1e10f;
		float uncor_hiparam = -1e10f;

		float samec_loparam = 1e10f;
		float samec_hiparam = -1e10f;

		processed_line4 l_uncor = uncor_plines[partition];
		processed_line4 l_samec = samec_plines[partition];

		// TODO: split up this loop due to too many temporaries; in particular,
		// the six line functions will consume 18 vector registers
		int texel_count = pt->texels_per_partition[partition];
		promise(texel_count > 0);

		int i = 0;

#if ASTCENC_SIMD_WIDTH > 1
		// This implementation is an example vectorization of this function.
		// It works for - the codec is a 2-4% faster than not vectorizing - but
		// the benefit is limited by the use of gathers and register pressure

		// Vectorize some useful scalar inputs
		vfloat l_uncor_bs0(l_uncor.bs.r);
		vfloat l_uncor_bs1(l_uncor.bs.g);
		vfloat l_uncor_bs2(l_uncor.bs.b);
		vfloat l_uncor_bs3(l_uncor.bs.a);

		vfloat l_uncor_amod0(l_uncor.amod.r);
		vfloat l_uncor_amod1(l_uncor.amod.g);
		vfloat l_uncor_amod2(l_uncor.amod.b);
		vfloat l_uncor_amod3(l_uncor.amod.a);

		vfloat l_uncor_bis0(l_uncor.bis.r);
		vfloat l_uncor_bis1(l_uncor.bis.g);
		vfloat l_uncor_bis2(l_uncor.bis.b);
		vfloat l_uncor_bis3(l_uncor.bis.a);

		vfloat l_samec_bs0(l_samec.bs.r);
		vfloat l_samec_bs1(l_samec.bs.g);
		vfloat l_samec_bs2(l_samec.bs.b);
		vfloat l_samec_bs3(l_samec.bs.a);

		vfloat l_samec_amod0(l_samec.amod.r);
		vfloat l_samec_amod1(l_samec.amod.g);
		vfloat l_samec_amod2(l_samec.amod.b);
		vfloat l_samec_amod3(l_samec.amod.a);

		vfloat l_samec_bis0(l_samec.bis.r);
		vfloat l_samec_bis1(l_samec.bis.g);
		vfloat l_samec_bis2(l_samec.bis.b);
		vfloat l_samec_bis3(l_samec.bis.a);

		vfloat uncor_loparamv(1e10f);
		vfloat uncor_hiparamv(-1e10f);
		vfloat uncor_errorsumv = vfloat::zero();

		vfloat samec_loparamv(1e10f);
		vfloat samec_hiparamv(-1e10f);
		vfloat samec_errorsumv = vfloat::zero();

		int clipped_texel_count = round_down_to_simd_multiple_vla(texel_count);
		for (/* */; i < clipped_texel_count; i += ASTCENC_SIMD_WIDTH)
		{
			vint texel_idxs(&(weights[i]));

			vfloat data_r = gatherf(blk->data_r, texel_idxs);
			vfloat data_g = gatherf(blk->data_g, texel_idxs);
			vfloat data_b = gatherf(blk->data_b, texel_idxs);
			vfloat data_a = gatherf(blk->data_a, texel_idxs);

			vfloat ew_r = gatherf(ewb->texel_weight_r, texel_idxs);
			vfloat ew_g = gatherf(ewb->texel_weight_g, texel_idxs);
			vfloat ew_b = gatherf(ewb->texel_weight_b, texel_idxs);
			vfloat ew_a = gatherf(ewb->texel_weight_a, texel_idxs);

			vfloat uncor_param  = (data_r * l_uncor_bs0)
			                    + (data_g * l_uncor_bs1)
			                    + (data_b * l_uncor_bs2)
			                    + (data_a * l_uncor_bs3);

			uncor_loparamv = min(uncor_param, uncor_loparamv);
			uncor_hiparamv = max(uncor_param, uncor_hiparamv);

			vfloat uncor_dist0 = (l_uncor_amod0 - data_r)
			                   + (uncor_param * l_uncor_bis0);
			vfloat uncor_dist1 = (l_uncor_amod1 - data_g)
			                   + (uncor_param * l_uncor_bis1);
			vfloat uncor_dist2 = (l_uncor_amod2 - data_b)
			                   + (uncor_param * l_uncor_bis2);
			vfloat uncor_dist3 = (l_uncor_amod3 - data_a)
			                   + (uncor_param * l_uncor_bis3);

			vfloat uncor_error = (ew_r * uncor_dist0 * uncor_dist0)
			                   + (ew_g * uncor_dist1 * uncor_dist1)
			                   + (ew_b * uncor_dist2 * uncor_dist2)
			                   + (ew_a * uncor_dist3 * uncor_dist3);

			uncor_errorsumv = uncor_errorsumv + uncor_error;

			// Process samechroma data
			vfloat samec_param = (data_r * l_samec_bs0)
			                   + (data_g * l_samec_bs1)
			                   + (data_b * l_samec_bs2)
			                   + (data_a * l_samec_bs3);

			samec_loparamv = min(samec_param, samec_loparamv);
			samec_hiparamv = max(samec_param, samec_hiparamv);


			vfloat samec_dist0 = (l_samec_amod0 - data_r)
			                   + (samec_param * l_samec_bis0);
			vfloat samec_dist1 = (l_samec_amod1 - data_g)
			                   + (samec_param * l_samec_bis1);
			vfloat samec_dist2 = (l_samec_amod2 - data_b)
			                   + (samec_param * l_samec_bis2);
			vfloat samec_dist3 = (l_samec_amod3 - data_a)
			                   + (samec_param * l_samec_bis3);

			vfloat samec_error = (ew_r * samec_dist0 * samec_dist0)
			                   + (ew_g * samec_dist1 * samec_dist1)
			                   + (ew_b * samec_dist2 * samec_dist2)
			                   + (ew_a * samec_dist3 * samec_dist3);

			samec_errorsumv = samec_errorsumv + samec_error;
		}

		uncor_loparam = hmin_s(uncor_loparamv);
		uncor_hiparam = hmax_s(uncor_hiparamv);
		uncor_errorsum += hadd_s(uncor_errorsumv);

		samec_loparam = hmin_s(samec_loparamv);
		samec_hiparam = hmax_s(samec_hiparamv);
		samec_errorsum += hadd_s(samec_errorsumv);
#endif

		// Loop tail
		for (/* */; i < texel_count; i++)
		{
			int iwt = weights[i];

			float4 dat = float4(blk->data_r[iwt],
			                    blk->data_g[iwt],
			                    blk->data_b[iwt],
			                    blk->data_a[iwt]);

			float4 ews = ewb->error_weights[iwt];

			float uncor_param = dot(dat, l_uncor.bs);
			uncor_loparam = astc::min(uncor_param, uncor_loparam);
			uncor_hiparam = astc::max(uncor_param, uncor_hiparam);

			float samec_param = dot(dat, l_samec.bs);
			samec_loparam = astc::min(samec_param, samec_loparam);
			samec_hiparam = astc::max(samec_param, samec_hiparam);

			float4 uncor_dist  = (l_uncor.amod - dat)
			                   + (uncor_param * l_uncor.bis);
			uncor_errorsum += dot(ews, uncor_dist * uncor_dist);

			float4 samec_dist = (l_samec.amod - dat)
			                  + (samec_param * l_samec.bis);
			samec_errorsum += dot(ews, samec_dist * samec_dist);
		}

		float uncor_linelen = uncor_hiparam - uncor_loparam;
		float samec_linelen = samec_hiparam - samec_loparam;

		// Turn very small numbers and NaNs into a small number
		uncor_linelen = astc::max(uncor_linelen, 1e-7f);
		samec_linelen = astc::max(samec_linelen, 1e-7f);

		uncor_lengths[partition] = uncor_linelen;
		samec_lengths[partition] = samec_linelen;
	}

	*uncor_errors = uncor_errorsum;
	*samec_errors = samec_errorsum;
}

void compute_error_squared_rgb(
	const partition_info *pt,
	const imageblock *blk,
	const error_weight_block *ewb,
	const processed_line3 *uncor_plines, // Uncorrelated channels
	const processed_line3 *samec_plines, // Same chroma channels
	float *uncor_lengths,
	float *samec_lengths,
	float *uncor_errors,
	float *samec_errors
) {
	float uncor_errorsum = 0.0f;
	float samec_errorsum = 0.0f;

	int partition_count = pt->partition_count;
	promise(partition_count > 0);

	for (int partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *weights = pt->texels_of_partition[partition];

		float uncor_loparam = 1e10f;
		float uncor_hiparam = -1e10f;

		float samec_loparam = 1e10f;
		float samec_hiparam = -1e10f;

		processed_line3 l_uncor = uncor_plines[partition];
		processed_line3 l_samec = samec_plines[partition];

		// TODO: split up this loop due to too many temporaries; in
		// particular, the six line functions will consume 18 vector registers
		int texel_count = pt->texels_per_partition[partition];
		promise(texel_count > 0);

		int i = 0;

#if ASTCENC_SIMD_WIDTH > 1
		// This implementation is an example vectorization of this function.
		// It works for - the codec is a 2-4% faster than not vectorizing - but
		// the benefit is limited by the use of gathers and register pressure

		// Vectorize some useful scalar inputs
		vfloat l_uncor_bs0(l_uncor.bs.r);
		vfloat l_uncor_bs1(l_uncor.bs.g);
		vfloat l_uncor_bs2(l_uncor.bs.b);

		vfloat l_uncor_amod0(l_uncor.amod.r);
		vfloat l_uncor_amod1(l_uncor.amod.g);
		vfloat l_uncor_amod2(l_uncor.amod.b);

		vfloat l_uncor_bis0(l_uncor.bis.r);
		vfloat l_uncor_bis1(l_uncor.bis.g);
		vfloat l_uncor_bis2(l_uncor.bis.b);

		vfloat l_samec_bs0(l_samec.bs.r);
		vfloat l_samec_bs1(l_samec.bs.g);
		vfloat l_samec_bs2(l_samec.bs.b);

		vfloat l_samec_amod0(l_samec.amod.r);
		vfloat l_samec_amod1(l_samec.amod.g);
		vfloat l_samec_amod2(l_samec.amod.b);

		vfloat l_samec_bis0(l_samec.bis.r);
		vfloat l_samec_bis1(l_samec.bis.g);
		vfloat l_samec_bis2(l_samec.bis.b);

		vfloat uncor_loparamv(1e10f);
		vfloat uncor_hiparamv(-1e10f);
		vfloat uncor_errorsumv = vfloat::zero();

		vfloat samec_loparamv(1e10f);
		vfloat samec_hiparamv(-1e10f);
		vfloat samec_errorsumv = vfloat::zero();

		int clipped_texel_count = round_down_to_simd_multiple_vla(texel_count);
		for (/* */; i < clipped_texel_count; i += ASTCENC_SIMD_WIDTH)
		{
			vint texel_idxs(&(weights[i]));

			vfloat data_r = gatherf(blk->data_r, texel_idxs);
			vfloat data_g = gatherf(blk->data_g, texel_idxs);
			vfloat data_b = gatherf(blk->data_b, texel_idxs);

			vfloat ew_r = gatherf(ewb->texel_weight_r, texel_idxs);
			vfloat ew_g = gatherf(ewb->texel_weight_g, texel_idxs);
			vfloat ew_b = gatherf(ewb->texel_weight_b, texel_idxs);

			vfloat uncor_param  = (data_r * l_uncor_bs0)
			                     + (data_g * l_uncor_bs1)
			                     + (data_b * l_uncor_bs2);

			uncor_loparamv = min(uncor_param, uncor_loparamv);
			uncor_hiparamv = max(uncor_param, uncor_hiparamv);

			vfloat uncor_dist0 = (l_uncor_amod0 - data_r)
			                   + (uncor_param * l_uncor_bis0);
			vfloat uncor_dist1 = (l_uncor_amod1 - data_g)
			                   + (uncor_param * l_uncor_bis1);
			vfloat uncor_dist2 = (l_uncor_amod2 - data_b)
			                   + (uncor_param * l_uncor_bis2);

			vfloat uncor_error = (ew_r * uncor_dist0 * uncor_dist0)
			                   + (ew_g * uncor_dist1 * uncor_dist1)
			                   + (ew_b * uncor_dist2 * uncor_dist2);

			uncor_errorsumv = uncor_errorsumv + uncor_error;

			// Process samechroma data
			vfloat samec_param = (data_r * l_samec_bs0)
			                   + (data_g * l_samec_bs1)
			                   + (data_b * l_samec_bs2);

			samec_loparamv = min(samec_param, samec_loparamv);
			samec_hiparamv = max(samec_param, samec_hiparamv);


			vfloat samec_dist0 = (l_samec_amod0 - data_r)
			                   + (samec_param * l_samec_bis0);
			vfloat samec_dist1 = (l_samec_amod1 - data_g)
			                   + (samec_param * l_samec_bis1);
			vfloat samec_dist2 = (l_samec_amod2 - data_b)
			                   + (samec_param * l_samec_bis2);

			vfloat samec_error = (ew_r * samec_dist0 * samec_dist0)
			                   + (ew_g * samec_dist1 * samec_dist1)
			                   + (ew_b * samec_dist2 * samec_dist2);

			samec_errorsumv = samec_errorsumv + samec_error;
		}

		uncor_loparam = hmin_s(uncor_loparamv);
		uncor_hiparam = hmax_s(uncor_hiparamv);
		uncor_errorsum += hadd_s(uncor_errorsumv);

		samec_loparam = hmin_s(samec_loparamv);
		samec_hiparam = hmax_s(samec_hiparamv);
		samec_errorsum += hadd_s(samec_errorsumv);
#endif

		// Loop tail
		for (/* */; i < texel_count; i++)
		{
			int iwt = weights[i];

			float3 dat = float3(blk->data_r[iwt],
			                    blk->data_g[iwt],
			                    blk->data_b[iwt]);

			float3 ews = float3(ewb->error_weights[iwt].r,
			                    ewb->error_weights[iwt].g,
			                    ewb->error_weights[iwt].b);

			float uncor_param = dot(dat, l_uncor.bs);
			uncor_loparam  = astc::min(uncor_param, uncor_loparam);
			uncor_hiparam = astc::max(uncor_param, uncor_hiparam);

			float samec_param = dot(dat, l_samec.bs);
			samec_loparam  = astc::min(samec_param, samec_loparam);
			samec_hiparam = astc::max(samec_param, samec_hiparam);

			float3 uncor_dist  = (l_uncor.amod - dat)
			                   + (uncor_param * l_uncor.bis);
			uncor_errorsum += dot(ews, uncor_dist * uncor_dist);

			float3 samec_dist = (l_samec.amod - dat)
			                  + (samec_param * l_samec.bis);
			samec_errorsum += dot(ews, samec_dist * samec_dist);
		}

		float uncor_linelen = uncor_hiparam - uncor_loparam;
		float samec_linelen = samec_hiparam - samec_loparam;

		// Turn very small numbers and NaNs into a small number
		uncor_linelen = astc::max(uncor_linelen, 1e-7f);
		samec_linelen = astc::max(samec_linelen, 1e-7f);

		uncor_lengths[partition] = uncor_linelen;
		samec_lengths[partition] = samec_linelen;
	}

	*uncor_errors = uncor_errorsum;
	*samec_errors = samec_errorsum;
}

#endif
