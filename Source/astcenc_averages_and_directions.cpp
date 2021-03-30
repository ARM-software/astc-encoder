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
void compute_avgs_and_dirs_4_comp(
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	partition_metrics pms[4]
) {
	int partition_count = pt->partition_count;
	promise(partition_count > 0);

	for (int partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *weights = pt->texels_of_partition[partition];

		vfloat4 base_sum = vfloat4::zero();
		float partition_weight = 0.0f;

		int texel_count = pt->partition_texel_count[partition];
		promise(texel_count > 0);

		for (int i = 0; i < texel_count; i++)
		{
			int iwt = weights[i];
			float weight = ewb->texel_weight[iwt];
			vfloat4 texel_datum = blk->texel(iwt);

			partition_weight += weight;
			base_sum += texel_datum * weight;
		}

		vfloat4 average = base_sum * (1.0f / astc::max(partition_weight, 1e-7f));
		pms[partition].avg = average * pms[partition].color_scale;

		vfloat4 sum_xp = vfloat4::zero();
		vfloat4 sum_yp = vfloat4::zero();
		vfloat4 sum_zp = vfloat4::zero();
		vfloat4 sum_wp = vfloat4::zero();

		for (int i = 0; i < texel_count; i++)
		{
			int iwt = weights[i];
			float weight = ewb->texel_weight[iwt];
			vfloat4 texel_datum = blk->texel(iwt);
			texel_datum = (texel_datum - average) * weight;

			vfloat4 zero = vfloat4::zero();

			vmask4 tdm0 = vfloat4(texel_datum.lane<0>()) > zero;
			sum_xp += select(zero, texel_datum, tdm0);

			vmask4 tdm1 = vfloat4(texel_datum.lane<1>()) > zero;
			sum_yp += select(zero, texel_datum, tdm1);

			vmask4 tdm2 = vfloat4(texel_datum.lane<2>()) > zero;
			sum_zp += select(zero, texel_datum, tdm2);

			vmask4 tdm3 = vfloat4(texel_datum.lane<3>()) > zero;
			sum_wp += select(zero, texel_datum, tdm3);
		}

		float prod_xp = dot_s(sum_xp, sum_xp);
		float prod_yp = dot_s(sum_yp, sum_yp);
		float prod_zp = dot_s(sum_zp, sum_zp);
		float prod_wp = dot_s(sum_wp, sum_wp);

		vfloat4 best_vector = sum_xp;
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

		pms[partition].dir = best_vector;
	}
}

void compute_avgs_and_dirs_3_comp(
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	int omitted_component,
	partition_metrics pm[4]
) {
	const float *texel_weights;
	const float* data_vr = blk->data_r;
	const float* data_vg = blk->data_g;
	const float* data_vb = blk->data_b;

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
		data_vg = blk->data_b;
		data_vb = blk->data_a;
	}
	else if (omitted_component == 2)
	{
		texel_weights = ewb->texel_weight_rga;
		data_vb = blk->data_a;
	}
	else
	{
		assert(omitted_component == 3);
		texel_weights = ewb->texel_weight_rgb;
	}

	int partition_count = pt->partition_count;
	promise(partition_count > 0);

	for (int partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *weights = pt->texels_of_partition[partition];

		vfloat4 base_sum = vfloat4::zero();
		float partition_weight = 0.0f;

		int texel_count = pt->partition_texel_count[partition];
		promise(texel_count > 0);

		for (int i = 0; i < texel_count; i++)
		{
			int iwt = weights[i];
			float weight = texel_weights[iwt];
			vfloat4 texel_datum(data_vr[iwt],
			                    data_vg[iwt],
			                    data_vb[iwt],
			                    0.0f);

			partition_weight += weight;
			base_sum += texel_datum * weight;
		}

		vfloat4 csf = pm[partition].color_scale;

		vfloat4 average = base_sum * (1.0f / astc::max(partition_weight, 1e-7f));
		pm[partition].avg = average * csf;

		vfloat4 sum_xp = vfloat4::zero();
		vfloat4 sum_yp = vfloat4::zero();
		vfloat4 sum_zp = vfloat4::zero();

		for (int i = 0; i < texel_count; i++)
		{
			int iwt = weights[i];
			float weight = texel_weights[iwt];
			vfloat4 texel_datum = vfloat4(data_vr[iwt],
			                              data_vg[iwt],
			                              data_vb[iwt],
			                              0.0f);
			texel_datum = (texel_datum - average) * weight;

			vfloat4 zero = vfloat4::zero();

			vmask4 tdm0 = vfloat4(texel_datum.lane<0>()) > zero;
			sum_xp += select(zero, texel_datum, tdm0);

			vmask4 tdm1 = vfloat4(texel_datum.lane<1>()) > zero;
			sum_yp += select(zero, texel_datum, tdm1);

			vmask4 tdm2 = vfloat4(texel_datum.lane<2>()) > zero;
			sum_zp += select(zero, texel_datum, tdm2);
		}

		float prod_xp = dot3_s(sum_xp, sum_xp);
		float prod_yp = dot3_s(sum_yp, sum_yp);
		float prod_zp = dot3_s(sum_zp, sum_zp);

		vfloat4 best_vector = sum_xp;
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

		if (dot3_s(best_vector, best_vector) < 1e-18f)
		{
			best_vector = vfloat4(1.0f, 1.0f, 1.0f, 0.0f);
		}

		pm[partition].dir = best_vector;
	}
}

void compute_avgs_and_dirs_2_comp(
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

		int texel_count = pt->partition_texel_count[partition];
		promise(texel_count > 0);

		for (int i = 0; i < texel_count; i++)
		{
			int iwt = weights[i];
			float weight = texel_weights[iwt];
			float2 texel_datum = float2(data_vr[iwt], data_vg[iwt]) * weight;
			partition_weight += weight;

			base_sum += texel_datum;
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
				sum_xp += texel_datum;
			}

			if (texel_datum.g > 0.0f)
			{
				sum_yp += texel_datum;
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
		const uint8_t *weights = pt->texels_of_partition[partition];

		float uncor_loparam = 1e10f;
		float uncor_hiparam = -1e10f;

		float samec_loparam = 1e10f;
		float samec_hiparam = -1e10f;

		processed_line4 l_uncor = uncor_plines[partition];
		processed_line4 l_samec = samec_plines[partition];

		int texel_count = pt->partition_texel_count[partition];
		promise(texel_count > 0);

		int i = 0;

		// Vectorize some useful scalar inputs
		vfloat l_uncor_bs0(l_uncor.bs.lane<0>());
		vfloat l_uncor_bs1(l_uncor.bs.lane<1>());
		vfloat l_uncor_bs2(l_uncor.bs.lane<2>());
		vfloat l_uncor_bs3(l_uncor.bs.lane<3>());

		vfloat l_uncor_amod0(l_uncor.amod.lane<0>());
		vfloat l_uncor_amod1(l_uncor.amod.lane<1>());
		vfloat l_uncor_amod2(l_uncor.amod.lane<2>());
		vfloat l_uncor_amod3(l_uncor.amod.lane<3>());

		vfloat l_uncor_bis0(l_uncor.bis.lane<0>());
		vfloat l_uncor_bis1(l_uncor.bis.lane<1>());
		vfloat l_uncor_bis2(l_uncor.bis.lane<2>());
		vfloat l_uncor_bis3(l_uncor.bis.lane<3>());

		vfloat l_samec_bs0(l_samec.bs.lane<0>());
		vfloat l_samec_bs1(l_samec.bs.lane<1>());
		vfloat l_samec_bs2(l_samec.bs.lane<2>());
		vfloat l_samec_bs3(l_samec.bs.lane<3>());

		assert(all(l_samec.amod == vfloat4(0.0f)));

		vfloat l_samec_bis0(l_samec.bis.lane<0>());
		vfloat l_samec_bis1(l_samec.bis.lane<1>());
		vfloat l_samec_bis2(l_samec.bis.lane<2>());
		vfloat l_samec_bis3(l_samec.bis.lane<3>());

		vfloat uncor_loparamv(1e10f);
		vfloat uncor_hiparamv(-1e10f);
		vfloat4 uncor_errorsumv = vfloat4::zero();

		vfloat samec_loparamv(1e10f);
		vfloat samec_hiparamv(-1e10f);
		vfloat4 samec_errorsumv = vfloat4::zero();

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

			haccumulate(uncor_errorsumv, uncor_error);

			// Process samechroma data
			vfloat samec_param = (data_r * l_samec_bs0)
			                   + (data_g * l_samec_bs1)
			                   + (data_b * l_samec_bs2)
			                   + (data_a * l_samec_bs3);

			samec_loparamv = min(samec_param, samec_loparamv);
			samec_hiparamv = max(samec_param, samec_hiparamv);

			vfloat samec_dist0 = samec_param * l_samec_bis0 - data_r;
			vfloat samec_dist1 = samec_param * l_samec_bis1 - data_g;
			vfloat samec_dist2 = samec_param * l_samec_bis2 - data_b;
			vfloat samec_dist3 = samec_param * l_samec_bis3 - data_a;

			vfloat samec_error = (ew_r * samec_dist0 * samec_dist0)
			                   + (ew_g * samec_dist1 * samec_dist1)
			                   + (ew_b * samec_dist2 * samec_dist2)
			                   + (ew_a * samec_dist3 * samec_dist3);

			haccumulate(samec_errorsumv, samec_error);
		}

		uncor_loparam = hmin_s(uncor_loparamv);
		uncor_hiparam = hmax_s(uncor_hiparamv);

		samec_loparam = hmin_s(samec_loparamv);
		samec_hiparam = hmax_s(samec_hiparamv);

		// Loop tail
		// Error is buffered and accumulated in blocks of 4 to ensure that
		// the partial sums added to the accumulator are invariant with the
		// vector implementation, irrespective of vector size ...
		alignas(16) float uncor_errorsum_tmp[4] { 0 };
		alignas(16) float samec_errorsum_tmp[4] { 0 };
		for (/* */; i < texel_count; i++)
		{
			int iwt = weights[i];

			vfloat4 dat = blk->texel(iwt);
			vfloat4 ews = ewb->error_weights[iwt];

			float uncor_param = dot_s(dat, l_uncor.bs);
			uncor_loparam = astc::min(uncor_param, uncor_loparam);
			uncor_hiparam = astc::max(uncor_param, uncor_hiparam);

			float samec_param = dot_s(dat, l_samec.bs);
			samec_loparam = astc::min(samec_param, samec_loparam);
			samec_hiparam = astc::max(samec_param, samec_hiparam);

			vfloat4 uncor_dist  = (l_uncor.amod - dat)
			                    + (uncor_param * l_uncor.bis);
			float uncor_error_tmp = dot_s(ews, uncor_dist * uncor_dist);

			vfloat4 samec_dist = samec_param * l_samec.bis - dat;
			float samec_error_tmp = dot_s(ews, samec_dist * samec_dist);

			// Accumulate error sum in the temporary array
			int error_index = i & 0x3;
			uncor_errorsum_tmp[error_index] = uncor_error_tmp;
			samec_errorsum_tmp[error_index] = samec_error_tmp;

#if ASTCENC_SIMD_WIDTH == 8
			// Zero the temporary staging buffer every 4 items unless last iter
			if ((i & 0x7) == 0x03)
			{
				haccumulate(uncor_errorsumv, vfloat4::loada(uncor_errorsum_tmp));
				storea(vfloat4::zero(), uncor_errorsum_tmp);

				haccumulate(samec_errorsumv, vfloat4::loada(samec_errorsum_tmp));
				storea(vfloat4::zero(), samec_errorsum_tmp);
			}
#endif
		}

		// Accumulate the loop tail using the vfloat4 swizzle
		haccumulate(uncor_errorsumv, vfloat4::loada(uncor_errorsum_tmp));
		haccumulate(samec_errorsumv, vfloat4::loada(samec_errorsum_tmp));

		// Resolve the final scalar accumulator sum
		haccumulate(uncor_errorsum, uncor_errorsumv);
		haccumulate(samec_errorsum, samec_errorsumv);

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
	partition_lines3 plines[4],
	float& uncor_error,
	float& samec_error
) {
	float uncor_errorsum = 0.0f;
	float samec_errorsum = 0.0f;

	int partition_count = pt->partition_count;
	promise(partition_count > 0);

	for (int partition = 0; partition < partition_count; partition++)
	{
		partition_lines3& pl = plines[partition];
		const uint8_t *weights = pt->texels_of_partition[partition];
		int texel_count = pt->partition_texel_count[partition];
		promise(texel_count > 0);

		float uncor_loparam = 1e10f;
		float uncor_hiparam = -1e10f;

		float samec_loparam = 1e10f;
		float samec_hiparam = -1e10f;

		processed_line3 l_uncor = pl.uncor_pline;
		processed_line3 l_samec = pl.samec_pline;

		int i = 0;

		// This implementation is an example vectorization of this function.
		// It works for - the codec is a 2-4% faster than not vectorizing - but
		// the benefit is limited by the use of gathers and register pressure

		// Vectorize some useful scalar inputs
		vfloat l_uncor_bs0(l_uncor.bs.lane<0>());
		vfloat l_uncor_bs1(l_uncor.bs.lane<1>());
		vfloat l_uncor_bs2(l_uncor.bs.lane<2>());

		vfloat l_uncor_amod0(l_uncor.amod.lane<0>());
		vfloat l_uncor_amod1(l_uncor.amod.lane<1>());
		vfloat l_uncor_amod2(l_uncor.amod.lane<2>());

		vfloat l_uncor_bis0(l_uncor.bis.lane<0>());
		vfloat l_uncor_bis1(l_uncor.bis.lane<1>());
		vfloat l_uncor_bis2(l_uncor.bis.lane<2>());

		vfloat l_samec_bs0(l_samec.bs.lane<0>());
		vfloat l_samec_bs1(l_samec.bs.lane<1>());
		vfloat l_samec_bs2(l_samec.bs.lane<2>());

		assert(all(l_samec.amod == vfloat4(0.0f)));

		vfloat l_samec_bis0(l_samec.bis.lane<0>());
		vfloat l_samec_bis1(l_samec.bis.lane<1>());
		vfloat l_samec_bis2(l_samec.bis.lane<2>());

		vfloat uncor_loparamv(1e10f);
		vfloat uncor_hiparamv(-1e10f);
		vfloat4 uncor_errorsumv = vfloat4::zero();

		vfloat samec_loparamv(1e10f);
		vfloat samec_hiparamv(-1e10f);
		vfloat4 samec_errorsumv = vfloat4::zero();

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

			vfloat uncor_err = (ew_r * uncor_dist0 * uncor_dist0)
			                 + (ew_g * uncor_dist1 * uncor_dist1)
			                 + (ew_b * uncor_dist2 * uncor_dist2);

			haccumulate(uncor_errorsumv, uncor_err);

			// Process samechroma data
			vfloat samec_param = (data_r * l_samec_bs0)
			                   + (data_g * l_samec_bs1)
			                   + (data_b * l_samec_bs2);

			samec_loparamv = min(samec_param, samec_loparamv);
			samec_hiparamv = max(samec_param, samec_hiparamv);


			vfloat samec_dist0 = samec_param * l_samec_bis0 - data_r;
			vfloat samec_dist1 = samec_param * l_samec_bis1 - data_g;
			vfloat samec_dist2 = samec_param * l_samec_bis2 - data_b;

			vfloat samec_err = (ew_r * samec_dist0 * samec_dist0)
			                 + (ew_g * samec_dist1 * samec_dist1)
			                 + (ew_b * samec_dist2 * samec_dist2);

			haccumulate(samec_errorsumv, samec_err);
		}

		uncor_loparam = hmin_s(uncor_loparamv);
		uncor_hiparam = hmax_s(uncor_hiparamv);

		samec_loparam = hmin_s(samec_loparamv);
		samec_hiparam = hmax_s(samec_hiparamv);

		// Loop tail
		// Error is buffered and accumulated in blocks of 4 to ensure that
		// the partial sums added to the accumulator are invariant with the
		// vector implementation, irrespective of vector size ...
		alignas(16) float uncor_errorsum_tmp[4] { 0 };
		alignas(16) float samec_errorsum_tmp[4] { 0 };
		for (/* */; i < texel_count; i++)
		{
			int iwt = weights[i];

			vfloat4 dat = blk->texel3(iwt);
			vfloat4 ews = ewb->error_weights[iwt];

			float uncor_param = dot3_s(dat, l_uncor.bs);
			uncor_loparam  = astc::min(uncor_param, uncor_loparam);
			uncor_hiparam = astc::max(uncor_param, uncor_hiparam);

			float samec_param = dot3_s(dat, l_samec.bs);
			samec_loparam  = astc::min(samec_param, samec_loparam);
			samec_hiparam = astc::max(samec_param, samec_hiparam);

			vfloat4 uncor_dist  = (l_uncor.amod - dat)
			                    + (uncor_param * l_uncor.bis);
			float uncor_error_tmp = dot3_s(ews, uncor_dist * uncor_dist);

			vfloat4 samec_dist = samec_param * l_samec.bis - dat;
			float samec_error_tmp = dot3_s(ews, samec_dist * samec_dist);

			// Accumulate error sum in the temporary array
			int error_index = i & 0x3;
			uncor_errorsum_tmp[error_index] = uncor_error_tmp;
			samec_errorsum_tmp[error_index] = samec_error_tmp;

#if ASTCENC_SIMD_WIDTH == 8
			// Emit the staging buffer every 4 items unless last iteration
			if ((i & 0x7) == 0x03)
			{
				haccumulate(uncor_errorsumv, vfloat4::loada(uncor_errorsum_tmp));
				storea(vfloat4::zero(), uncor_errorsum_tmp);

				haccumulate(samec_errorsumv, vfloat4::loada(samec_errorsum_tmp));
				storea(vfloat4::zero(), samec_errorsum_tmp);
			}
#endif
		}

		// Accumulate the loop tail using the vfloat4 swizzle
		haccumulate(uncor_errorsumv, vfloat4::loada(uncor_errorsum_tmp));
		haccumulate(samec_errorsumv, vfloat4::loada(samec_errorsum_tmp));

		// Resolve the final scalar accumulator sum
		haccumulate(uncor_errorsum, uncor_errorsumv);
		haccumulate(samec_errorsum, samec_errorsumv);

		float uncor_linelen = uncor_hiparam - uncor_loparam;
		float samec_linelen = samec_hiparam - samec_loparam;

		// Turn very small numbers and NaNs into a small number
		uncor_linelen = astc::max(uncor_linelen, 1e-7f);
		samec_linelen = astc::max(samec_linelen, 1e-7f);

		pl.uncor_line_len = uncor_linelen;
		pl.samec_line_len = samec_linelen;
	}

	uncor_error = uncor_errorsum;
	samec_error = samec_errorsum;
}

#endif
