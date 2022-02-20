// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2011-2022 Arm Limited
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

/* See header for documentation. */
void compute_avgs_and_dirs_4_comp(
	const partition_info& pi,
	const image_block& blk,
	partition_metrics pm[BLOCK_MAX_PARTITIONS]
) {
	int partition_count = pi.partition_count;
	promise(partition_count > 0);

	float texel_weight = hadd_s(blk.channel_weight) / 4.0f;
	vfloat4 channel_weight = blk.channel_weight;

	for (int partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *texel_indexes = pi.texels_of_partition[partition];

		vfloat4 base_sum = vfloat4::zero();
		float partition_weight = 0.0f;

		unsigned int texel_count = pi.partition_texel_count[partition];
		promise(texel_count > 0);

		for (unsigned int i = 0; i < texel_count; i++)
		{
			int iwt = texel_indexes[i];
			vfloat4 texel_datum = blk.texel(iwt);
			partition_weight += texel_weight;
			base_sum += texel_datum * texel_weight;
		}

		vfloat4 error_sum = channel_weight;
		vfloat4 csf = normalize(sqrt(error_sum)) * 2.0f;
		vfloat4 average = base_sum / astc::max(partition_weight, 1e-7f);

		pm[partition].error_weight = error_sum;
		pm[partition].avg = average * csf;
		pm[partition].color_scale = csf;
		pm[partition].icolor_scale = 1.0f / max(csf, 1e-7f);

		vfloat4 sum_xp = vfloat4::zero();
		vfloat4 sum_yp = vfloat4::zero();
		vfloat4 sum_zp = vfloat4::zero();
		vfloat4 sum_wp = vfloat4::zero();

		for (unsigned int i = 0; i < texel_count; i++)
		{
			unsigned int iwt = texel_indexes[i];
			vfloat4 texel_datum = blk.texel(iwt);
			texel_datum = (texel_datum - average) * texel_weight;

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

/* See header for documentation. */
void compute_avgs_and_dirs_3_comp(
	const partition_info& pi,
	const image_block& blk,
	unsigned int omitted_component,
	partition_metrics pm[BLOCK_MAX_PARTITIONS]
) {
	float texel_weight = hadd_s(blk.channel_weight.swz<0, 1, 2>()) / 3.0f;
	vfloat4 channel_weight = blk.channel_weight.swz<0, 1, 2>();

	const float* data_vr = blk.data_r;
	const float* data_vg = blk.data_g;
	const float* data_vb = blk.data_b;

	if (omitted_component == 0)
	{
		texel_weight = hadd_s(blk.channel_weight.swz<1, 2, 3>()) / 3.0f;
		channel_weight = blk.channel_weight.swz<1, 2, 3>();

		data_vr = blk.data_g;
		data_vg = blk.data_b;
		data_vb = blk.data_a;
	}
	else if (omitted_component == 1)
	{
		texel_weight = hadd_s(blk.channel_weight.swz<0, 2, 3>()) / 3.0f;
		channel_weight = blk.channel_weight.swz<0, 2, 3>();

		data_vg = blk.data_b;
		data_vb = blk.data_a;
	}
	else if (omitted_component == 2)
	{
		texel_weight = hadd_s(blk.channel_weight.swz<0, 1, 3>()) / 3.0f;
		channel_weight = blk.channel_weight.swz<0, 1, 3>();

		data_vb = blk.data_a;
	}

	unsigned int partition_count = pi.partition_count;
	promise(partition_count > 0);

	for (unsigned int partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *texel_indexes = pi.texels_of_partition[partition];

		vfloat4 base_sum = vfloat4::zero();
		float partition_weight = 0.0f;

		unsigned int texel_count = pi.partition_texel_count[partition];
		promise(texel_count > 0);

		for (unsigned int i = 0; i < texel_count; i++)
		{
			unsigned int iwt = texel_indexes[i];

			vfloat4 texel_datum = vfloat3(data_vr[iwt],
			                              data_vg[iwt],
			                              data_vb[iwt]);

			partition_weight += texel_weight;
			base_sum += texel_datum * texel_weight;
		}

		vfloat4 error_sum = channel_weight;
		vfloat4 csf = normalize(sqrt(error_sum)) * 1.73205080f;

		vfloat4 average = base_sum / astc::max(partition_weight, 1e-7f);

		pm[partition].error_weight = error_sum;
		pm[partition].avg = average * csf;
		pm[partition].color_scale = csf;
		pm[partition].icolor_scale = 1.0f / max(csf, 1e-7f);

		vfloat4 sum_xp = vfloat4::zero();
		vfloat4 sum_yp = vfloat4::zero();
		vfloat4 sum_zp = vfloat4::zero();

		for (unsigned int i = 0; i < texel_count; i++)
		{
			unsigned int iwt = texel_indexes[i];

			vfloat4 texel_datum = vfloat3(data_vr[iwt],
			                              data_vg[iwt],
			                              data_vb[iwt]);

			texel_datum = (texel_datum - average) * texel_weight;

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

		pm[partition].dir = best_vector;
	}
}

/* See header for documentation. */
void compute_avgs_and_dirs_3_comp_rgb(
	const partition_info& pi,
	const image_block& blk,
	partition_metrics pm[BLOCK_MAX_PARTITIONS]
) {
	float texel_weight = hadd_s(blk.channel_weight.swz<0, 1, 2>()) / 3;
	vfloat4 channel_weight = blk.channel_weight.swz<0, 1, 2>();

	unsigned int partition_count = pi.partition_count;
	promise(partition_count > 0);

	for (unsigned int partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *texel_indexes = pi.texels_of_partition[partition];

		vfloat4 base_sum = vfloat4::zero();
		float partition_weight = 0.0f;

		unsigned int texel_count = pi.partition_texel_count[partition];
		promise(texel_count > 0);

		for (unsigned int i = 0; i < texel_count; i++)
		{
			unsigned int iwt = texel_indexes[i];

			vfloat4 texel_datum = blk.texel3(iwt);

			partition_weight += texel_weight;
			base_sum += texel_datum * texel_weight;
		}

		vfloat4 error_sum = channel_weight;
		vfloat4 csf = normalize(sqrt(error_sum)) * 1.73205080f;

		vfloat4 average = base_sum / astc::max(partition_weight, 1e-7f);

		pm[partition].error_weight = error_sum;
		pm[partition].avg = average * csf;
		pm[partition].color_scale = csf;
		pm[partition].icolor_scale = 1.0f / max(csf, 1e-7f);

		vfloat4 sum_xp = vfloat4::zero();
		vfloat4 sum_yp = vfloat4::zero();
		vfloat4 sum_zp = vfloat4::zero();

		for (unsigned int i = 0; i < texel_count; i++)
		{
			unsigned int iwt = texel_indexes[i];

			vfloat4 texel_datum = blk.texel3(iwt);

			texel_datum = (texel_datum - average) * texel_weight;

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

		pm[partition].dir = best_vector;
	}
}

/* See header for documentation. */
void compute_avgs_and_dirs_2_comp(
	const partition_info& pt,
	const image_block& blk,
	unsigned int component1,
	unsigned int component2,
	partition_metrics pm[BLOCK_MAX_PARTITIONS]
) {
	float texel_weight;
	vfloat4 channel_weight;

	const float* data_vr = nullptr;
	const float* data_vg = nullptr;

	if (component1 == 0 && component2 == 1)
	{
		texel_weight = hadd_s(blk.channel_weight.swz<0, 1>()) / 2.0f;
		channel_weight = blk.channel_weight.swz<0, 1>();

		data_vr = blk.data_r;
		data_vg = blk.data_g;
	}
	else if (component1 == 0 && component2 == 2)
	{
		texel_weight = hadd_s(blk.channel_weight.swz<0, 2>()) / 2.0f;
		channel_weight = blk.channel_weight.swz<0, 2>();

		data_vr = blk.data_r;
		data_vg = blk.data_b;
	}
	else // (component1 == 1 && component2 == 2)
	{
		assert(component1 == 1 && component2 == 2);

		texel_weight = hadd_s(blk.channel_weight.swz<1, 2>()) / 2.0f;
		channel_weight = blk.channel_weight.swz<1, 2>();

		data_vr = blk.data_g;
		data_vg = blk.data_b;
	}

	unsigned int partition_count = pt.partition_count;
	promise(partition_count > 0);

	for (unsigned int partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *texel_indexes = pt.texels_of_partition[partition];

		vfloat4 base_sum = vfloat4::zero();
		float partition_weight = 0.0f;

		unsigned int texel_count = pt.partition_texel_count[partition];
		promise(texel_count > 0);

		for (unsigned int i = 0; i < texel_count; i++)
		{
			unsigned int iwt = texel_indexes[i];
			vfloat4 texel_datum = vfloat2(data_vr[iwt], data_vg[iwt]) * texel_weight;

			partition_weight += texel_weight;
			base_sum += texel_datum;
		}

		vfloat4 error_sum = channel_weight;
		vfloat4 csf = normalize(sqrt(error_sum)) * 1.41421356f;
		vfloat4 average = base_sum / astc::max(partition_weight, 1e-7f);


		pm[partition].error_weight = error_sum;
		pm[partition].avg = average * csf;
		pm[partition].color_scale = csf;
		pm[partition].icolor_scale = 1.0f / max(csf, 1e-7f);

		vfloat4 sum_xp = vfloat4::zero();
		vfloat4 sum_yp = vfloat4::zero();

		for (unsigned int i = 0; i < texel_count; i++)
		{
			unsigned int iwt = texel_indexes[i];
			vfloat4 texel_datum = vfloat2(data_vr[iwt], data_vg[iwt]);
			texel_datum = (texel_datum - average) * texel_weight;

			vfloat4 zero = vfloat4::zero();

			vmask4 tdm0 = vfloat4(texel_datum.lane<0>()) > zero;
			sum_xp += select(zero, texel_datum, tdm0);

			vmask4 tdm1 = vfloat4(texel_datum.lane<1>()) > zero;
			sum_yp += select(zero, texel_datum, tdm1);
		}

		float prod_xp = dot_s(sum_xp, sum_xp);
		float prod_yp = dot_s(sum_yp, sum_yp);

		vfloat4 best_vector = sum_xp;
		float best_sum = prod_xp;

		if (prod_yp > best_sum)
		{
			best_vector = sum_yp;
		}

		pm[partition].dir = best_vector;
	}
}

/* See header for documentation. */
void compute_error_squared_rgba(
	const partition_info& pi,
	const image_block& blk,
	const processed_line4 uncor_plines[BLOCK_MAX_PARTITIONS],
	const processed_line4 samec_plines[BLOCK_MAX_PARTITIONS],
	float uncor_lengths[BLOCK_MAX_PARTITIONS],
	float samec_lengths[BLOCK_MAX_PARTITIONS],
	float& uncor_error,
	float& samec_error
) {
	unsigned int partition_count = pi.partition_count;
	promise(partition_count > 0);

	uncor_error = 0.0f;
	samec_error = 0.0f;

	for (unsigned int partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *texel_indexes = pi.texels_of_partition[partition];

		float uncor_loparam = 1e10f;
		float uncor_hiparam = -1e10f;

		float samec_loparam = 1e10f;
		float samec_hiparam = -1e10f;

		processed_line4 l_uncor = uncor_plines[partition];
		processed_line4 l_samec = samec_plines[partition];

		unsigned int texel_count = pi.partition_texel_count[partition];
		promise(texel_count > 0);

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

		vfloat ew_r(blk.channel_weight.lane<0>());
		vfloat ew_g(blk.channel_weight.lane<1>());
		vfloat ew_b(blk.channel_weight.lane<2>());
		vfloat ew_a(blk.channel_weight.lane<3>());

		// This implementation over-shoots, but this is safe as we initialize the texel_indexes
		// array to extend the last value. This means min/max are not impacted, but we need to mask
		// out the dummy values when we compute the line weighting.
		vint lane_ids = vint::lane_id();
		for (unsigned int i = 0; i < texel_count; i += ASTCENC_SIMD_WIDTH)
		{
			vmask mask = lane_ids < vint(texel_count);
			vint texel_idxs(&(texel_indexes[i]));

			vfloat data_r = gatherf(blk.data_r, texel_idxs);
			vfloat data_g = gatherf(blk.data_g, texel_idxs);
			vfloat data_b = gatherf(blk.data_b, texel_idxs);
			vfloat data_a = gatherf(blk.data_a, texel_idxs);

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

			lane_ids += vint(ASTCENC_SIMD_WIDTH);
		}

		uncor_loparam = hmin_s(uncor_loparamv);
		uncor_hiparam = hmax_s(uncor_hiparamv);

		samec_loparam = hmin_s(samec_loparamv);
		samec_hiparam = hmax_s(samec_hiparamv);

		// Resolve the final scalar accumulator sum
		haccumulate(uncor_error, uncor_errorsumv);
		haccumulate(samec_error, samec_errorsumv);

		float uncor_linelen = uncor_hiparam - uncor_loparam;
		float samec_linelen = samec_hiparam - samec_loparam;

		// Turn very small numbers and NaNs into a small number
		uncor_lengths[partition] = astc::max(uncor_linelen, 1e-7f);
		samec_lengths[partition] = astc::max(samec_linelen, 1e-7f);
	}
}

/* See header for documentation. */
void compute_error_squared_rgb(
	const partition_info& pi,
	const image_block& blk,
	partition_lines3 plines[BLOCK_MAX_PARTITIONS],
	float& uncor_error,
	float& samec_error
) {
	unsigned int partition_count = pi.partition_count;
	promise(partition_count > 0);

	uncor_error = 0.0f;
	samec_error = 0.0f;

	for (unsigned int partition = 0; partition < partition_count; partition++)
	{
		partition_lines3& pl = plines[partition];
		const uint8_t *texel_indexes = pi.texels_of_partition[partition];
		unsigned int texel_count = pi.partition_texel_count[partition];
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

		vfloat ew_r(blk.channel_weight.lane<0>());
		vfloat ew_g(blk.channel_weight.lane<1>());
		vfloat ew_b(blk.channel_weight.lane<2>());

		// This implementation over-shoots, but this is safe as we initialize the weights array
		// to extend the last value. This means min/max are not impacted, but we need to mask
		// out the dummy values when we compute the line weighting.
		vint lane_ids = vint::lane_id();
		for (unsigned int i = 0; i < texel_count; i += ASTCENC_SIMD_WIDTH)
		{
			vmask mask = lane_ids < vint(texel_count);
			vint texel_idxs(&(texel_indexes[i]));

			vfloat data_r = gatherf(blk.data_r, texel_idxs);
			vfloat data_g = gatherf(blk.data_g, texel_idxs);
			vfloat data_b = gatherf(blk.data_b, texel_idxs);

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

			lane_ids += vint(ASTCENC_SIMD_WIDTH);
		}

		uncor_loparam = hmin_s(uncor_loparamv);
		uncor_hiparam = hmax_s(uncor_hiparamv);

		samec_loparam = hmin_s(samec_loparamv);
		samec_hiparam = hmax_s(samec_hiparamv);

		// Resolve the final scalar accumulator sum
		haccumulate(uncor_error, uncor_errorsumv);
		haccumulate(samec_error, samec_errorsumv);

		float uncor_linelen = uncor_hiparam - uncor_loparam;
		float samec_linelen = samec_hiparam - samec_loparam;

		// Turn very small numbers and NaNs into a small number
		pl.uncor_line_len = astc::max(uncor_linelen, 1e-7f);
		pl.samec_line_len = astc::max(samec_linelen, 1e-7f);
	}
}

#endif
