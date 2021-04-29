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

// For a full block, functions to compute averages and dominant directions. The
// averages and directions are computed separately for each partition.
// We have separate versions for blocks with and without alpha, since the
// processing for blocks with alpha is significantly more expensive. The
// direction vectors it produces are NOT normalized.
void compute_avgs_and_dirs_4_comp(
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	partition_metrics pm[4]
) {
	int partition_count = pt->partition_count;
	promise(partition_count > 0);

	for (int partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *weights = pt->texels_of_partition[partition];

		vfloat4 error_sum = vfloat4::zero();
		vfloat4 base_sum = vfloat4::zero();
		vfloat4 rgba_min(1e38f);
		vfloat4 rgba_max(-1e38f);
		float partition_weight = 0.0f;

		int texel_count = pt->partition_texel_count[partition];
		promise(texel_count > 0);

		for (int i = 0; i < texel_count; i++)
		{
			int iwt = weights[i];
			float weight = ewb->texel_weight[iwt];
			vfloat4 texel_datum = blk->texel(iwt);
			vfloat4 error_weight = ewb->error_weights[iwt];

			if (weight > 1e-10f)
			{
				rgba_min = min(texel_datum, rgba_min);
				rgba_max = max(texel_datum, rgba_max);
			}

			partition_weight += weight;
			base_sum += texel_datum * weight;
			error_sum += error_weight;
		}

		error_sum = error_sum / texel_count;
		vfloat4 csf = normalize(sqrt(error_sum)) * 2.0f;

		vfloat4 average = base_sum * (1.0f / astc::max(partition_weight, 1e-7f));

		pm[partition].error_weight = error_sum;
		pm[partition].avg = average * csf;
		pm[partition].color_scale = csf;
		pm[partition].icolor_scale = 1.0f / max(csf, 1e-7f);
		vfloat4 range = max(rgba_max - rgba_min, 1e-10f);
		pm[partition].range_sq = range * range;

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

		pm[partition].dir = best_vector;
	}
}

void compute_avgs_and_dirs_3_comp(
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	int omitted_component,
	partition_metrics pm[4]
) {
	const float *texel_weights = ewb->texel_weight_rgb;

	const float* data_vr = blk->data_r;
	const float* data_vg = blk->data_g;
	const float* data_vb = blk->data_b;

	const float* error_vr = ewb->texel_weight_r;
	const float* error_vg = ewb->texel_weight_g;
	const float* error_vb = ewb->texel_weight_b;

	if (omitted_component == 0)
	{
		texel_weights = ewb->texel_weight_gba;

		data_vr = blk->data_g;
		data_vg = blk->data_b;
		data_vb = blk->data_a;

		error_vr = ewb->texel_weight_g;
		error_vg = ewb->texel_weight_b;
		error_vb = ewb->texel_weight_a;
	}
	else if (omitted_component == 1)
	{
		texel_weights = ewb->texel_weight_rba;

		data_vg = blk->data_b;
		data_vb = blk->data_a;

		error_vg = ewb->texel_weight_b;
		error_vb = ewb->texel_weight_a;
	}
	else if (omitted_component == 2)
	{
		texel_weights = ewb->texel_weight_rga;

		data_vb = blk->data_a;

		error_vb = ewb->texel_weight_a;
	}

	int partition_count = pt->partition_count;
	promise(partition_count > 0);

	for (int partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *weights = pt->texels_of_partition[partition];

		vfloat4 error_sum = vfloat4::zero();
		vfloat4 base_sum = vfloat4::zero();
		vfloat4 rgb_min(1e38f);
		vfloat4 rgb_max(-1e38f);
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

			vfloat4 error_weight(error_vr[iwt],
			                     error_vg[iwt],
			                     error_vb[iwt],
			                     0.0f);

			if (weight > 1e-10f)
			{
				rgb_min = min(texel_datum, rgb_min);
				rgb_max = max(texel_datum, rgb_max);
			}

			partition_weight += weight;
			base_sum += texel_datum * weight;
			error_sum += error_weight;
		}

		error_sum = error_sum / texel_count;
		vfloat4 csf = normalize(sqrt(error_sum)) * 1.73205080f;

		vfloat4 average = base_sum * (1.0f / astc::max(partition_weight, 1e-7f));

		pm[partition].error_weight = error_sum;
		pm[partition].avg = average * csf;
		pm[partition].color_scale = csf;
		pm[partition].icolor_scale = 1.0f / max(csf, 1e-7f);
		vfloat4 range = max(rgb_max - rgb_min, 1e-10f);
		pm[partition].range_sq = range * range;

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
	int component1,
	int component2,
	partition_metrics pm[4]
) {
	const float *texel_weights;

	const float* data_vr = nullptr;
	const float* data_vg = nullptr;

	const float* error_vr = nullptr;
	const float* error_vg = nullptr;

	if (component1 == 0 && component2 == 1)
	{
		texel_weights = ewb->texel_weight_rg;

		data_vr = blk->data_r;
		data_vg = blk->data_g;

		error_vr = ewb->texel_weight_r;
		error_vg = ewb->texel_weight_g;
	}
	else if (component1 == 0 && component2 == 2)
	{
		texel_weights = ewb->texel_weight_rb;

		data_vr = blk->data_r;
		data_vg = blk->data_b;

		error_vr = ewb->texel_weight_r;
		error_vg = ewb->texel_weight_b;
	}
	else // (component1 == 1 && component2 == 2)
	{
		assert(component1 == 1 && component2 == 2);
		texel_weights = ewb->texel_weight_gb;

		data_vr = blk->data_g;
		data_vg = blk->data_b;


		error_vr = ewb->texel_weight_g;
		error_vg = ewb->texel_weight_b;
	}

	int partition_count = pt->partition_count;
	promise(partition_count > 0);

	for (int partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *weights = pt->texels_of_partition[partition];

		float2 error_sum = float2(0.0f);
		float2 base_sum = float2(0.0f);
		float partition_weight = 0.0f;

		int texel_count = pt->partition_texel_count[partition];
		promise(texel_count > 0);

		for (int i = 0; i < texel_count; i++)
		{
			int iwt = weights[i];
			float weight = texel_weights[iwt];
			float2 texel_datum = float2(data_vr[iwt], data_vg[iwt]) * weight;

			float2 error_weight = float2(error_vr[iwt], error_vg[iwt]);

			partition_weight += weight;
			base_sum += texel_datum;
			error_sum += error_weight;
		}

		vfloat4 error_sum4 = vfloat4(error_sum.r, error_sum.g, 0.0f, 0.0f) * (1.0f / texel_count);
		vfloat4 csf = normalize(sqrt(error_sum4)) * 1.41421356f;
		vfloat4 average4 = vfloat4(base_sum.r, base_sum.g, 0.0f, 0.0f) * (1.0f / astc::max(partition_weight, 1e-7f));


		pm[partition].error_weight = error_sum4;
		pm[partition].avg = average4 * csf;
		float2 average = float2(average4.lane<0>(), average4.lane<1>());
		pm[partition].color_scale = csf;
		pm[partition].icolor_scale = 1.0f / max(csf, 1e-7f);

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

		pm[partition].dir = vfloat4(best_vector.r, best_vector.g, 0.0f, 0.0f);
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

		// This implementation over-shoots, but this is safe as we initialize the weights array
		// to extend the last value. This means min/max are not impacted, but we need to mask
		// out the dummy values when we compute the line weighting.
		vint lane_ids = vint::lane_id();
		for (/* */; i < texel_count; i += ASTCENC_SIMD_WIDTH)
		{
			vmask mask = lane_ids < vint(texel_count);
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

			vfloat uncor_err = (ew_r * uncor_dist0 * uncor_dist0)
			                 + (ew_g * uncor_dist1 * uncor_dist1)
			                 + (ew_b * uncor_dist2 * uncor_dist2)
			                 + (ew_a * uncor_dist3 * uncor_dist3);

			uncor_err = select(vfloat::zero(), uncor_err, mask);
			haccumulate(uncor_errorsumv, uncor_err);

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

			vfloat samec_err = (ew_r * samec_dist0 * samec_dist0)
			                 + (ew_g * samec_dist1 * samec_dist1)
			                 + (ew_b * samec_dist2 * samec_dist2)
			                 + (ew_a * samec_dist3 * samec_dist3);

			samec_err = select(vfloat::zero(), samec_err, mask);
			haccumulate(samec_errorsumv, samec_err);

			lane_ids = lane_ids + vint(ASTCENC_SIMD_WIDTH);
		}

		uncor_loparam = hmin_s(uncor_loparamv);
		uncor_hiparam = hmax_s(uncor_hiparamv);

		samec_loparam = hmin_s(samec_loparamv);
		samec_hiparam = hmax_s(samec_hiparamv);

		// Resolve the final scalar accumulator sum
		haccumulate(uncor_errorsum, uncor_errorsumv);
		haccumulate(samec_errorsum, samec_errorsumv);

		float uncor_linelen = uncor_hiparam - uncor_loparam;
		float samec_linelen = samec_hiparam - samec_loparam;

		// Turn very small numbers and NaNs into a small number
		uncor_lengths[partition] = astc::max(uncor_linelen, 1e-7f);
		samec_lengths[partition] = astc::max(samec_linelen, 1e-7f);
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

		// This implementation over-shoots, but this is safe as we initialize the weights array
		// to extend the last value. This means min/max are not impacted, but we need to mask
		// out the dummy values when we compute the line weighting.
		vint lane_ids = vint::lane_id();
		for (int i = 0; i < texel_count; i += ASTCENC_SIMD_WIDTH)
		{
			vmask mask = lane_ids < vint(texel_count);
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

			uncor_err = select(vfloat::zero(), uncor_err, mask);
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

			samec_err = select(vfloat::zero(), samec_err, mask);
			haccumulate(samec_errorsumv, samec_err);

			lane_ids = lane_ids + vint(ASTCENC_SIMD_WIDTH);
		}

		uncor_loparam = hmin_s(uncor_loparamv);
		uncor_hiparam = hmax_s(uncor_hiparamv);

		samec_loparam = hmin_s(samec_loparamv);
		samec_hiparam = hmax_s(samec_hiparamv);

		// Resolve the final scalar accumulator sum
		haccumulate(uncor_errorsum, uncor_errorsumv);
		haccumulate(samec_errorsum, samec_errorsumv);

		float uncor_linelen = uncor_hiparam - uncor_loparam;
		float samec_linelen = samec_hiparam - samec_loparam;

		// Turn very small numbers and NaNs into a small number
		pl.uncor_line_len = astc::max(uncor_linelen, 1e-7f);
		pl.samec_line_len = astc::max(samec_linelen, 1e-7f);
	}

	uncor_error = uncor_errorsum;
	samec_error = samec_errorsum;
}

#endif
