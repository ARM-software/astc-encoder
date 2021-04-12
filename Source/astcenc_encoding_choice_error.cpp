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

#if !defined(ASTCENC_DECOMPRESS_ONLY)

/**
 * @brief Functions for finding color error post-compression.
 *
 * We assume there are two independent sources of error in any given partition.
 * - encoding choice errors
 * - quantization errors
 *
 * Encoding choice errors are caused by encoder decisions, such as:
 * - using luminance rather than RGB.
 * - using RGB+scale instead of two full RGB endpoints.
 * - dropping the alpha component.
 *
 * Quantization errors occur due to the limited precision we use for storage.
 * These errors generally scale with quantization level, but are not actually
 * independent of color encoding. In particular:
 * - if we can use offset encoding then quantization error is halved.
 * - if we can use blue-contraction, quantization error for RG is halved.
 * - quantization error is higher for the HDR endpoint modes.
 * Other than these errors, quantization error is assumed to be proportional to
 * the quantization step.
 */

#include "astcenc_internal.h"

// helper function to merge two endpoint-colors
void merge_endpoints(
	const endpoints* ep1,	// contains three of the color components
	const endpoints* ep2,	// contains the remaining color component
	int plane2_component,
	endpoints* res
) {
	int partition_count = ep1->partition_count;
	vmask4 sep_mask = vint4::lane_id() == vint4(plane2_component);

	res->partition_count = partition_count;
	promise(partition_count > 0);
	for (int i = 0; i < partition_count; i++)
	{
		res->endpt0[i] = select(ep1->endpt0[i], ep2->endpt0[i], sep_mask);
		res->endpt1[i] = select(ep1->endpt1[i], ep2->endpt1[i], sep_mask);
	}
}

// function to compute the error across a tile when using a particular line for
// a particular partition.
static void compute_error_squared_rgb_single_partition(
	int partition_to_test,
	const partition_info* pt,	// the partition that we use when computing the squared-error.
	const imageblock* blk,
	const error_weight_block* ewb,
	const processed_line3* uncor_pline,
	float* uncor_err,
	const processed_line3* samec_pline,
	float* samec_err,
	const processed_line3* rgbl_pline,
	float* rgbl_err,
	const processed_line3* l_pline,
	float* l_err,
	float* a_drop_err
) {
	float uncor_errorsum = 0.0f;
	float samec_errorsum = 0.0f;
	float rgbl_errorsum = 0.0f;
	float l_errorsum = 0.0f;
	float a_drop_errorsum = 0.0f;

	int texels_in_partition = pt->partition_texel_count[partition_to_test];
	promise(texels_in_partition > 0);

	for (int i = 0; i < texels_in_partition; i++)
	{
		int tix = pt->texels_of_partition[partition_to_test][i];
		float texel_weight = ewb->texel_weight_rgb[tix];
		if (texel_weight < 1e-20f)
		{
			continue;
		}

		vfloat4 point = blk->texel(tix);
		vfloat4 ews = ewb->error_weights[tix];

		// Compute the error that arises from just ditching alpha
		float default_alpha = imageblock_default_alpha(blk);
		float omalpha = point.lane<3>() - default_alpha;
		a_drop_errorsum += omalpha * omalpha * ews.lane<3>();

		{
			float param = dot3_s(point, uncor_pline->bs);
			vfloat4 rp1 = uncor_pline->amod + param * uncor_pline->bis;
			vfloat4 dist = rp1 - point;
			uncor_errorsum += dot3_s(ews, dist * dist);
		}

		{
			float param = dot3_s(point, samec_pline->bs);
			// No samec amod - we know it's always zero
			vfloat4 rp1 = /* samec_pline->amod + */ param * samec_pline->bis;
			vfloat4 dist = rp1 - point;
			samec_errorsum += dot3_s(ews, dist * dist);
		}

		{
			float param = dot3_s(point,  rgbl_pline->bs);
			vfloat4 rp1 = rgbl_pline->amod + param * rgbl_pline->bis;
			vfloat4 dist = rp1 - point;
			rgbl_errorsum += dot3_s(ews, dist * dist);
		}

		{
			float param = dot3_s(point, l_pline->bs);
			// No luma amod - we know it's always zero
			vfloat4 rp1 = /* l_pline->amod + */ param * l_pline->bis;
			vfloat4 dist = rp1 - point;
			l_errorsum += dot3_s(ews, dist * dist);
		}
	}

	*uncor_err = uncor_errorsum;
	*samec_err = samec_errorsum;
	*rgbl_err = rgbl_errorsum;
	*l_err = l_errorsum;
	*a_drop_err = a_drop_errorsum;
}

/*
   for a given set of input colors and a given partitioning, determine: color error that results
   from RGB-scale encoding (relevant for LDR only) color error that results from RGB-lumashift encoding
   (relevant for HDR only) color error that results from luminance-encoding color error that results
   form dropping alpha. whether we are eligible for offset encoding whether we are eligible for
   blue-contraction

   The input data are: color data partitioning error-weight data
 */
void compute_encoding_choice_errors(
	const block_size_descriptor* bsd,
	const imageblock* blk,
	const partition_info* pt,
	const error_weight_block* ewb,
	int plane2_component,	// component that is separated out in 2-plane mode, -1 in 1-plane mode
	encoding_choice_errors* eci)
{
	int partition_count = pt->partition_count;
	int texels_per_block = bsd->texel_count;

	promise(partition_count > 0);
	promise(texels_per_block > 0);

	partition_metrics pms[4];

	compute_partition_error_color_weightings(*ewb, *pt, pms);

	compute_avgs_and_dirs_3_comp(pt, blk, ewb, 3, pms);

	endpoints ep;
	if (plane2_component == -1)
	{
		endpoints_and_weights ei;
		compute_endpoints_and_ideal_weights_1_plane(bsd, pt, blk, ewb, &ei);
		ep = ei.ep;
	}
	else
	{
		endpoints_and_weights ei1, ei2;
		compute_endpoints_and_ideal_weights_2_planes(bsd, pt, blk, ewb, plane2_component, &ei1, &ei2);
		merge_endpoints(&(ei1.ep), &(ei2.ep), plane2_component, &ep);
	}

	for (int i = 0; i < partition_count; i++)
	{
		partition_metrics& pm = pms[i];

		// TODO: Can we skip rgb_luma_lines for LDR images?
		line3 uncor_rgb_lines;
		line3 samec_rgb_lines;	// for LDR-RGB-scale
		line3 rgb_luma_lines;	// for HDR-RGB-scale

		processed_line3 uncor_rgb_plines;
		processed_line3 samec_rgb_plines;	// for LDR-RGB-scale
		processed_line3 rgb_luma_plines;	// for HDR-RGB-scale
		processed_line3 luminance_plines;

		float uncorr_rgb_error;
		float samechroma_rgb_error;
		float rgb_luma_error;
		float luminance_rgb_error;
		float alpha_drop_error;

		vfloat4 csf = pm.color_scale;
		csf.set_lane<3>(0.0f);
		vfloat4 csfn = normalize(csf);

		vfloat4 icsf = pm.icolor_scale;
		icsf.set_lane<3>(0.0f);

		uncor_rgb_lines.a = pm.avg;
		uncor_rgb_lines.b = normalize_safe(pm.dir.swz<0, 1, 2>(), csfn);

		samec_rgb_lines.a = vfloat4::zero();
		samec_rgb_lines.b = normalize_safe(pm.avg.swz<0, 1, 2>(), csfn);

		rgb_luma_lines.a = pm.avg;
		rgb_luma_lines.b = csfn;

		uncor_rgb_plines.amod = (uncor_rgb_lines.a - uncor_rgb_lines.b * dot3(uncor_rgb_lines.a, uncor_rgb_lines.b)) * icsf;
		uncor_rgb_plines.bs   = uncor_rgb_lines.b * csf;
		uncor_rgb_plines.bis  = uncor_rgb_lines.b * icsf;

		// Same chroma always goes though zero, so this is simpler than the others
		samec_rgb_plines.amod = vfloat4::zero();
		samec_rgb_plines.bs   = samec_rgb_lines.b * csf;
		samec_rgb_plines.bis  = samec_rgb_lines.b * icsf;

		rgb_luma_plines.amod = (rgb_luma_lines.a - rgb_luma_lines.b * dot3(rgb_luma_lines.a, rgb_luma_lines.b)) * icsf;
		rgb_luma_plines.bs   = rgb_luma_lines.b * csf;
		rgb_luma_plines.bis  = rgb_luma_lines.b * icsf;

		// Luminance always goes though zero, so this is simpler than the others
		luminance_plines.amod = vfloat4::zero();
		luminance_plines.bs   = csfn * csf;
		luminance_plines.bis  = csfn * icsf;

		compute_error_squared_rgb_single_partition(
		    i, pt, blk, ewb,
		    &uncor_rgb_plines, &uncorr_rgb_error,
		    &samec_rgb_plines, &samechroma_rgb_error,
		    &rgb_luma_plines,  &rgb_luma_error,
		    &luminance_plines, &luminance_rgb_error,
		                       &alpha_drop_error);

		// Determine if we can offset encode RGB lanes
		vfloat4 endpt0 = ep.endpt0[i];
		vfloat4 endpt1 = ep.endpt1[i];
		vfloat4 endpt_diff = abs(endpt1 - endpt0);
		vmask4 endpt_can_offset = endpt_diff < vfloat4(0.12f * 65535.0f);
		bool can_offset_encode = (mask(endpt_can_offset) & 0x7) == 0x7;

		// Determine if we can blue contract encode RGB lanes
		vfloat4 endpt_diff_bc(
			endpt0.lane<0>() + (endpt0.lane<0>() - endpt0.lane<2>()),
			endpt1.lane<0>() + (endpt1.lane<0>() - endpt1.lane<2>()),
			endpt0.lane<1>() + (endpt0.lane<1>() - endpt0.lane<2>()),
			endpt1.lane<1>() + (endpt1.lane<1>() - endpt1.lane<2>())
		);

		vmask4 endpt_can_bc_lo = endpt_diff_bc > vfloat4(0.01f * 65535.0f);
		vmask4 endpt_can_bc_hi = endpt_diff_bc < vfloat4(0.99f * 65535.0f);
		bool can_blue_contract = (mask(endpt_can_bc_lo & endpt_can_bc_hi) & 0x7) == 0x7;

		// Store out the settings
		eci[i].rgb_scale_error = (samechroma_rgb_error - uncorr_rgb_error) * 0.7f;	// empirical
		eci[i].rgb_luma_error  = (rgb_luma_error - uncorr_rgb_error) * 1.5f;	// wild guess
		eci[i].luminance_error = (luminance_rgb_error - uncorr_rgb_error) * 3.0f;	// empirical
		eci[i].alpha_drop_error = alpha_drop_error * 3.0f;
		eci[i].can_offset_encode = can_offset_encode;
		eci[i].can_blue_contract = can_blue_contract;
	}
}

#endif
