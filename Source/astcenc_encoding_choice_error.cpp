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
 * @brief Functions for finding color error post-compression.
 *
 * We assume there are two independent sources of error in any given partition.
 * - encoding choice errors
 * - quantization errors
 *
 * Encoding choice errors are caused by encoder decisions, such as:
 * - using luminance rather than RGB.
 * - using RGB+scale instead of two full RGB endpoints.
 * - dropping the alpha channel.
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
	const endpoints * ep1,	// contains three of the color components
	const endpoints * ep2,	// contains the remaining color component
	int separate_component,
	endpoints * res
) {
	int partition_count = ep1->partition_count;
	res->partition_count = partition_count;
	for (int i = 0; i < partition_count; i++)
	{
		res->endpt0[i] = ep1->endpt0[i];
		res->endpt1[i] = ep1->endpt1[i];
	}

	switch (separate_component)
	{
	case 0:
		for (int i = 0; i < partition_count; i++)
		{
			res->endpt0[i].r = ep2->endpt0[i].r;
			res->endpt1[i].r = ep2->endpt1[i].r;
		}
		break;
	case 1:
		for (int i = 0; i < partition_count; i++)
		{
			res->endpt0[i].g = ep2->endpt0[i].g;
			res->endpt1[i].g = ep2->endpt1[i].g;
		}
		break;
	case 2:
		for (int i = 0; i < partition_count; i++)
		{
			res->endpt0[i].b = ep2->endpt0[i].b;
			res->endpt1[i].b = ep2->endpt1[i].b;
		}
		break;
	case 3:
		for (int i = 0; i < partition_count; i++)
		{
			res->endpt0[i].a = ep2->endpt0[i].a;
			res->endpt1[i].a = ep2->endpt1[i].a;
		}
		break;
	}
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
	const imageblock* pb,
	const partition_info* pi,
	const error_weight_block* ewb,
	int separate_component,	// component that is separated out in 2-plane mode, -1 in 1-plane mode
	encoding_choice_errors* eci)
{
	int partition_count = pi->partition_count;
	int texels_per_block = bsd->texel_count;

	float3 averages[4];
	float3 directions_rgb[4];
	float4 error_weightings[4];
	float4 color_scalefactors[4];
	float4 inverse_color_scalefactors[4];

	compute_partition_error_color_weightings(bsd, ewb, pi, error_weightings, color_scalefactors);
	compute_averages_and_directions_rgb(pi, pb, ewb, color_scalefactors, averages, directions_rgb);

	line3 uncorr_rgb_lines[4];
	line3 samechroma_rgb_lines[4];	// for LDR-RGB-scale
	line3 rgb_luma_lines[4];	// for HDR-RGB-scale
	line3 luminance_lines[4];

	processed_line3 proc_uncorr_rgb_lines[4];
	processed_line3 proc_samechroma_rgb_lines[4];	// for LDR-RGB-scale
	processed_line3 proc_rgb_luma_lines[4];	// for HDR-RGB-scale
	processed_line3 proc_luminance_lines[4];

	for (int i = 0; i < partition_count; i++)
	{
		inverse_color_scalefactors[i].r = 1.0f / MAX(color_scalefactors[i].r, 1e-7f);
		inverse_color_scalefactors[i].g = 1.0f / MAX(color_scalefactors[i].g, 1e-7f);
		inverse_color_scalefactors[i].b = 1.0f / MAX(color_scalefactors[i].b, 1e-7f);
		inverse_color_scalefactors[i].a = 1.0f / MAX(color_scalefactors[i].a, 1e-7f);

		float3 csf = float3(color_scalefactors[i].r, color_scalefactors[i].g, color_scalefactors[i].b);
		float3 icsf = float3(inverse_color_scalefactors[i].r, inverse_color_scalefactors[i].g, inverse_color_scalefactors[i].b);

		uncorr_rgb_lines[i].a = averages[i];
		if (dot(directions_rgb[i], directions_rgb[i]) == 0.0f)
		{
			uncorr_rgb_lines[i].b = normalize(csf);
		}
		else
		{
			uncorr_rgb_lines[i].b = normalize(directions_rgb[i]);
		}

		samechroma_rgb_lines[i].a = float3(0.0f);
		if (dot(averages[i], averages[i]) < 1e-20f)
		{
			samechroma_rgb_lines[i].b = normalize(csf);
		}
		else
		{
			samechroma_rgb_lines[i].b = normalize(averages[i]);
		}

		rgb_luma_lines[i].a = averages[i];
		rgb_luma_lines[i].b = normalize(csf);

		luminance_lines[i].a = float3(0.0f);
		luminance_lines[i].b = normalize(csf);

		proc_uncorr_rgb_lines[i].amod = (uncorr_rgb_lines[i].a - uncorr_rgb_lines[i].b * dot(uncorr_rgb_lines[i].a, uncorr_rgb_lines[i].b)) * icsf;
		proc_uncorr_rgb_lines[i].bs = uncorr_rgb_lines[i].b * csf;
		proc_uncorr_rgb_lines[i].bis = uncorr_rgb_lines[i].b * icsf;

		proc_samechroma_rgb_lines[i].amod = (samechroma_rgb_lines[i].a - samechroma_rgb_lines[i].b * dot(samechroma_rgb_lines[i].a, samechroma_rgb_lines[i].b)) * icsf;
		proc_samechroma_rgb_lines[i].bs = samechroma_rgb_lines[i].b * csf;
		proc_samechroma_rgb_lines[i].bis = samechroma_rgb_lines[i].b * icsf;

		proc_rgb_luma_lines[i].amod = (rgb_luma_lines[i].a - rgb_luma_lines[i].b * dot(rgb_luma_lines[i].a, rgb_luma_lines[i].b)) * icsf;
		proc_rgb_luma_lines[i].bs = rgb_luma_lines[i].b * csf;
		proc_rgb_luma_lines[i].bis = rgb_luma_lines[i].b * icsf;

		proc_luminance_lines[i].amod = (luminance_lines[i].a - luminance_lines[i].b * dot(luminance_lines[i].a, luminance_lines[i].b)) * icsf;
		proc_luminance_lines[i].bs = luminance_lines[i].b * csf;
		proc_luminance_lines[i].bis = luminance_lines[i].b * icsf;
	}

	float uncorr_rgb_error[4];
	float samechroma_rgb_error[4];
	float rgb_luma_error[4];
	float luminance_rgb_error[4];

	for (int i = 0; i < partition_count; i++)
	{
		uncorr_rgb_error[i] = compute_error_squared_rgb_single_partition(i, bsd, pi, pb, ewb, &(proc_uncorr_rgb_lines[i]));

		samechroma_rgb_error[i] = compute_error_squared_rgb_single_partition(i, bsd, pi, pb, ewb, &(proc_samechroma_rgb_lines[i]));

		rgb_luma_error[i] = compute_error_squared_rgb_single_partition(i, bsd, pi, pb, ewb, &(proc_rgb_luma_lines[i]));

		luminance_rgb_error[i] = compute_error_squared_rgb_single_partition(i, bsd, pi, pb, ewb, &(proc_luminance_lines[i]));
	}

	// compute the error that arises from just ditching alpha and RGB
	float alpha_drop_error[4];
	float rgb_drop_error[4];

	for (int i = 0; i < partition_count; i++)
	{
		alpha_drop_error[i] = 0;
		rgb_drop_error[i] = 0;
	}

	for (int i = 0; i < texels_per_block; i++)
	{
		int partition = pi->partition_of_texel[i];
		float alpha = pb->data_a[i];
		float default_alpha = pb->alpha_lns[i] ? (float)0x7800 : (float)0xFFFF;

		float omalpha = alpha - default_alpha;
		alpha_drop_error[partition] += omalpha * omalpha * ewb->error_weights[i].a;

		float red = pb->data_r[i];
		float green = pb->data_g[i];
		float blue = pb->data_b[i];
		rgb_drop_error[partition] += red * red * ewb->error_weights[i].r +
		                             green * green * ewb->error_weights[i].g +
		                             blue * blue * ewb->error_weights[i].b;
	}

	// check if we are eligible for blue-contraction and offset-encoding
	endpoints ep;
	if (separate_component == -1)
	{
		endpoints_and_weights ei;
		compute_endpoints_and_ideal_weights_1_plane(bsd, pi, pb, ewb, &ei);
		ep = ei.ep;
	}
	else
	{
		endpoints_and_weights ei1, ei2;
		compute_endpoints_and_ideal_weights_2_planes(bsd, pi, pb, ewb, separate_component, &ei1, &ei2);

		merge_endpoints(&(ei1.ep), &(ei2.ep), separate_component, &ep);
	}

	int eligible_for_offset_encode[4];
	int eligible_for_blue_contraction[4];
	for (int i = 0; i < partition_count; i++)
	{
		float4 endpt0 = ep.endpt0[i];
		float4 endpt1 = ep.endpt1[i];

		float4 endpt_dif = endpt1 - endpt0;
		if (fabsf(endpt_dif.r) < (0.12f * 65535.0f) &&
		    fabsf(endpt_dif.g) < (0.12f * 65535.0f) &&
		    fabsf(endpt_dif.b) < (0.12f * 65535.0f))
		{
			eligible_for_offset_encode[i] = 1;
		}
		else
		{
			eligible_for_offset_encode[i] = 0;
		}

		endpt0.r += (endpt0.r - endpt0.b);
		endpt0.g += (endpt0.g - endpt0.b);
		endpt1.r += (endpt1.r - endpt1.b);
		endpt1.g += (endpt1.g - endpt1.b);
		if (endpt0.r > (0.01f * 65535.0f) && endpt0.r < (0.99f * 65535.0f) &&
		    endpt1.r > (0.01f * 65535.0f) && endpt1.r < (0.99f * 65535.0f) &&
		    endpt0.g > (0.01f * 65535.0f) && endpt0.g < (0.99f * 65535.0f) &&
		    endpt1.g > (0.01f * 65535.0f) && endpt1.g < (0.99f * 65535.0f))
		{
			eligible_for_blue_contraction[i] = 1;
		}
		else
		{
			eligible_for_blue_contraction[i] = 0;
		}
	}

	// finally, gather up our results
	for (int i = 0; i < partition_count; i++)
	{
		eci[i].rgb_scale_error = (samechroma_rgb_error[i] - uncorr_rgb_error[i]) * 0.7f;	// empirical
		eci[i].rgb_luma_error = (rgb_luma_error[i] - uncorr_rgb_error[i]) * 1.5f;	// wild guess
		eci[i].luminance_error = (luminance_rgb_error[i] - uncorr_rgb_error[i]) * 3.0f;	// empirical
		eci[i].alpha_drop_error = alpha_drop_error[i] * 3.0f;
		eci[i].rgb_drop_error = rgb_drop_error[i] * 3.0f;
		eci[i].can_offset_encode = eligible_for_offset_encode[i];
		eci[i].can_blue_contract = eligible_for_blue_contraction[i];
	}
}

#endif
