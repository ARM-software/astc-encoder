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
 * @brief Functions for computing color endpoints and texel weights.
 */

#include "astc_codec_internals.h"

#ifdef DEBUG_PRINT_DIAGNOSTICS
	#include <stdio.h>
#endif

#ifdef DEBUG_CAPTURE_NAN
	#ifndef _GNU_SOURCE
		#define _GNU_SOURCE
	#endif

	#include <fenv.h>
#endif

static void compute_endpoints_and_ideal_weights_1_component(
	const block_size_descriptor* bsd,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	endpoints_and_weights* ei,
	int component
) {
	int i;

	int partition_count = pt->partition_count;
	ei->ep.partition_count = partition_count;

	float lowvalues[4], highvalues[4];
	float partition_error_scale[4];
	float linelengths_rcp[4];

	int texels_per_block = bsd->texel_count;

	const float *error_weights;
	const float* data_vr = nullptr;
	switch (component)
	{
	case 0:
		error_weights = ewb->texel_weight_r;
		data_vr = blk->data_r;
		break;
	case 1:
		error_weights = ewb->texel_weight_g;
		data_vr = blk->data_g;
		break;
	case 2:
		error_weights = ewb->texel_weight_b;
		data_vr = blk->data_b;
		break;
	case 3:
		error_weights = ewb->texel_weight_a;
		data_vr = blk->data_a;
		break;
	default:
		ASTC_CODEC_INTERNAL_ERROR();
	}

	for (i = 0; i < partition_count; i++)
	{
		lowvalues[i] = 1e10f;
		highvalues[i] = -1e10f;
	}

	for (i = 0; i < texels_per_block; i++)
	{
		if (error_weights[i] > 1e-10f)
		{
			float value = data_vr[i];
			int partition = pt->partition_of_texel[i];
			if (value < lowvalues[partition])
				lowvalues[partition] = value;
			if (value > highvalues[partition])
				highvalues[partition] = value;
		}
	}

	for (i = 0; i < partition_count; i++)
	{
		float diff = highvalues[i] - lowvalues[i];
		if (diff < 0)
		{
			lowvalues[i] = 0.0f;
			highvalues[i] = 0.0f;
		}
		if (diff < 1e-7f)
			diff = 1e-7f;
		partition_error_scale[i] = diff * diff;
		linelengths_rcp[i] = 1.0f / diff;
	}

	for (i = 0; i < texels_per_block; i++)
	{
		float value = data_vr[i];
		int partition = pt->partition_of_texel[i];
		value -= lowvalues[partition];
		value *= linelengths_rcp[partition];
		if (value > 1.0f)
			value = 1.0f;
		else if (!(value > 0.0f))
			value = 0.0f;

		ei->weights[i] = value;
		ei->weight_error_scale[i] = partition_error_scale[partition] * error_weights[i];
		if (astc::isnan(ei->weight_error_scale[i]))
		{
			ASTC_CODEC_INTERNAL_ERROR();
		}
	}

	for (i = 0; i < partition_count; i++)
	{
		ei->ep.endpt0[i] = float4(blk->red_min, blk->green_min, blk->blue_min, blk->alpha_min);
		ei->ep.endpt1[i] = float4(blk->red_max, blk->green_max, blk->blue_max, blk->alpha_max);
		switch (component)
		{
		case 0:				// red/x
			ei->ep.endpt0[i].x = lowvalues[i];
			ei->ep.endpt1[i].x = highvalues[i];
			break;
		case 1:				// green/y
			ei->ep.endpt0[i].y = lowvalues[i];
			ei->ep.endpt1[i].y = highvalues[i];
			break;
		case 2:				// blue/z
			ei->ep.endpt0[i].z = lowvalues[i];
			ei->ep.endpt1[i].z = highvalues[i];
			break;
		case 3:				// alpha/w
			ei->ep.endpt0[i].w = lowvalues[i];
			ei->ep.endpt1[i].w = highvalues[i];
			break;
		}
	}

	// print all the data that this function computes.
	#ifdef DEBUG_PRINT_DIAGNOSTICS
		if (print_diagnostics)
		{
			printf("%s: %dx%dx%d texels, %d partitions, component=%d\n", __func__, xdim, ydim, zdim, partition_count, component);
			printf("Endpoints:\n");
			for (i = 0; i < partition_count; i++)
			{
				printf("%d Low: <%g> => <%g %g %g %g>\n", i, lowvalues[i], ei->ep.endpt0[i].x, ei->ep.endpt0[i].y, ei->ep.endpt0[i].z, ei->ep.endpt0[i].w);
				printf("%d High: <%g> => <%g %g %g %g>\n", i, highvalues[i], ei->ep.endpt1[i].x, ei->ep.endpt1[i].y, ei->ep.endpt1[i].z, ei->ep.endpt1[i].w);
			}
			printf("Ideal-weights:\n");

			for (i = 0; i < texels_per_block; i++)
			{
				printf("%3d <%2d %2d %2d>=> %g (weight=%g)\n", i, i % xdim, (i / xdim) % ydim, i / (xdim * ydim), ei->weights[i], ei->weight_error_scale[i]);
			}
			printf("\n");
		}
	#endif
}

static void compute_endpoints_and_ideal_weights_2_components(
	const block_size_descriptor* bsd,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block * ewb,
	endpoints_and_weights* ei,
	int component1,
	int component2
) {
	int i;

	int partition_count = pt->partition_count;
	ei->ep.partition_count = partition_count;

	float4 error_weightings[4];
	float4 color_scalefactors[4];

	float2 scalefactors[4];

	const float *error_weights;
	const float* data_vr = nullptr;
	const float* data_vg = nullptr;
	if (component1 == 0 && component2 == 1)
	{
		error_weights = ewb->texel_weight_rg;
		data_vr = blk->data_r;
		data_vg = blk->data_g;
	}
	else if (component1 == 0 && component2 == 2)
	{
		error_weights = ewb->texel_weight_rb;
		data_vr = blk->data_r;
		data_vg = blk->data_b;
	}
	else if (component1 == 1 && component2 == 2)
	{
		error_weights = ewb->texel_weight_gb;
		data_vr = blk->data_g;
		data_vg = blk->data_b;
	}
	else
	{
		ASTC_CODEC_INTERNAL_ERROR();
	}

	int texels_per_block = bsd->texel_count;

	compute_partition_error_color_weightings(bsd, ewb, pt, error_weightings, color_scalefactors);

	for (i = 0; i < partition_count; i++)
	{
		float s1 = 0, s2 = 0;
		switch (component1)
		{
		case 0:
			s1 = color_scalefactors[i].x;
			break;
		case 1:
			s1 = color_scalefactors[i].y;
			break;
		case 2:
			s1 = color_scalefactors[i].z;
			break;
		case 3:
			s1 = color_scalefactors[i].w;
			break;
		}

		switch (component2)
		{
		case 0:
			s2 = color_scalefactors[i].x;
			break;
		case 1:
			s2 = color_scalefactors[i].y;
			break;
		case 2:
			s2 = color_scalefactors[i].z;
			break;
		case 3:
			s2 = color_scalefactors[i].w;
			break;
		}
		scalefactors[i] = normalize(float2(s1, s2)) * 1.41421356f;
	}

	float lowparam[4], highparam[4];

	float2 averages[4];
	float2 directions[4];

	line2 lines[4];
	float scale[4];
	float length_squared[4];

	for (i = 0; i < partition_count; i++)
	{
		lowparam[i] = 1e10;
		highparam[i] = -1e10;
	}

	compute_averages_and_directions_2_components(pt, blk, ewb, scalefactors, component1, component2, averages, directions);

	for (i = 0; i < partition_count; i++)
	{
		float2 egv = directions[i];
		if (egv.x + egv.y < 0.0f)
			directions[i] = float2(0, 0) - egv;
	}

	for (i = 0; i < partition_count; i++)
	{
		lines[i].a = averages[i];
		if (dot(directions[i], directions[i]) == 0.0f)
			lines[i].b = normalize(float2(1, 1));
		else
			lines[i].b = normalize(directions[i]);
	}

	for (i = 0; i < texels_per_block; i++)
	{
		if (error_weights[i] > 1e-10f)
		{
			int partition = pt->partition_of_texel[i];
			float2 point = float2(data_vr[i], data_vg[i]) * scalefactors[partition];
			line2 l = lines[partition];
			float param = dot(point - l.a, l.b);
			ei->weights[i] = param;
			if (param < lowparam[partition])
				lowparam[partition] = param;
			if (param > highparam[partition])
				highparam[partition] = param;
		}
		else
		{
			ei->weights[i] = -1e38f;
		}
	}

	float2 lowvalues[4];
	float2 highvalues[4];

	for (i = 0; i < partition_count; i++)
	{
		float length = highparam[i] - lowparam[i];
		if (length < 0.0f)			// case for when none of the texels had any weight
		{
			lowparam[i] = 0.0f;
			highparam[i] = 1e-7f;
		}

		// it is possible for a uniform-color partition to produce length=0; this
		// causes NaN-production and NaN-propagation later on. Set length to
		// a small value to avoid this problem.
		if (length < 1e-7f)
			length = 1e-7f;

		length_squared[i] = length * length;
		scale[i] = 1.0f / length;

		float2 ep0 = lines[i].a + lines[i].b * lowparam[i];
		float2 ep1 = lines[i].a + lines[i].b * highparam[i];

		ep0.x /= scalefactors[i].x;
		ep0.y /= scalefactors[i].y;

		ep1.x /= scalefactors[i].x;
		ep1.y /= scalefactors[i].y;

		lowvalues[i] = ep0;
		highvalues[i] = ep1;
	}

	for (i = 0; i < partition_count; i++)
	{
		ei->ep.endpt0[i] = float4(blk->red_min, blk->green_min, blk->blue_min, blk->alpha_min);
		ei->ep.endpt1[i] = float4(blk->red_max, blk->green_max, blk->blue_max, blk->alpha_max);

		float2 ep0 = lowvalues[i];
		float2 ep1 = highvalues[i];

		switch (component1)
		{
		case 0:
			ei->ep.endpt0[i].x = ep0.x;
			ei->ep.endpt1[i].x = ep1.x;
			break;
		case 1:
			ei->ep.endpt0[i].y = ep0.x;
			ei->ep.endpt1[i].y = ep1.x;
			break;
		case 2:
			ei->ep.endpt0[i].z = ep0.x;
			ei->ep.endpt1[i].z = ep1.x;
			break;
		case 3:
			ei->ep.endpt0[i].w = ep0.x;
			ei->ep.endpt1[i].w = ep1.x;
			break;
		}

		switch (component2)
		{
		case 0:
			ei->ep.endpt0[i].x = ep0.y;
			ei->ep.endpt1[i].x = ep1.y;
			break;
		case 1:
			ei->ep.endpt0[i].y = ep0.y;
			ei->ep.endpt1[i].y = ep1.y;
			break;
		case 2:
			ei->ep.endpt0[i].z = ep0.y;
			ei->ep.endpt1[i].z = ep1.y;
			break;
		case 3:
			ei->ep.endpt0[i].w = ep0.y;
			ei->ep.endpt1[i].w = ep1.y;
			break;
		}
	}

	for (i = 0; i < texels_per_block; i++)
	{
		int partition = pt->partition_of_texel[i];
		float idx = (ei->weights[i] - lowparam[partition]) * scale[partition];
		if (idx > 1.0f)
			idx = 1.0f;
		else if (!(idx > 0.0f))
			idx = 0.0f;

		ei->weights[i] = idx;
		ei->weight_error_scale[i] = length_squared[partition] * error_weights[i];
		if (astc::isnan(ei->weight_error_scale[i]))
		{
			ASTC_CODEC_INTERNAL_ERROR();
		}
	}

	// print all the data that this function computes.
	#ifdef DEBUG_PRINT_DIAGNOSTICS
		if (print_diagnostics)
		{
			printf("%s: %dx%dx%d texels, %d partitions, component1=%d, component2=%d\n", __func__, xdim, ydim, zdim, partition_count, component1, component2);
			printf("Endpoints:\n");
			for (i = 0; i < partition_count; i++)
			{
				printf("%d Low: <%g %g> => <%g %g %g %g>\n", i, lowvalues[i].x, lowvalues[i].y, ei->ep.endpt0[i].x, ei->ep.endpt0[i].y, ei->ep.endpt0[i].z, ei->ep.endpt0[i].w);
				printf("%d High: <%g %g> => <%g %g %g %g>\n", i, highvalues[i].x, highvalues[i].y, ei->ep.endpt1[i].x, ei->ep.endpt1[i].y, ei->ep.endpt1[i].z, ei->ep.endpt1[i].w);
			}
			printf("Ideal-weights:\n");

			for (i = 0; i < texels_per_block; i++)
			{
				printf("%3d <%2d %2d %2d>=> %g (weight=%g)\n", i, i % xdim, (i / xdim) % ydim, i / (xdim * ydim), ei->weights[i], ei->weight_error_scale[i]);
			}
			printf("\n");
		}
	#endif
}

static void compute_endpoints_and_ideal_weights_3_components(
	const block_size_descriptor* bsd,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	endpoints_and_weights* ei,
	int omittedComponent
) {
	int i;

	int partition_count = pt->partition_count;
	ei->ep.partition_count = partition_count;

	float4 error_weightings[4];
	float4 color_scalefactors[4];

	float3 scalefactors[4];

	int texels_per_block = bsd->texel_count;

	const float *error_weights;
	const float* data_vr = nullptr;
	const float* data_vg = nullptr;
	const float* data_vb = nullptr;
	if (omittedComponent == 0)
	{
		error_weights = ewb->texel_weight_gba;
		data_vr = blk->data_g;
		data_vg = blk->data_b;
		data_vb = blk->data_a;
	}
	else if (omittedComponent == 1)
	{
		error_weights = ewb->texel_weight_rba;
		data_vr = blk->data_r;
		data_vg = blk->data_b;
		data_vb = blk->data_a;
	}
	else if (omittedComponent == 2)
	{
		error_weights = ewb->texel_weight_rga;
		data_vr = blk->data_r;
		data_vg = blk->data_g;
		data_vb = blk->data_a;
	}
	else
	{
		error_weights = ewb->texel_weight_rgb;
		data_vr = blk->data_r;
		data_vg = blk->data_g;
		data_vb = blk->data_b;
	}

	compute_partition_error_color_weightings(bsd, ewb, pt, error_weightings, color_scalefactors);

	for (i = 0; i < partition_count; i++)
	{
		float s1 = 0, s2 = 0, s3 = 0;
		switch (omittedComponent)
		{
		case 0:
			s1 = color_scalefactors[i].y;
			s2 = color_scalefactors[i].z;
			s3 = color_scalefactors[i].w;
			break;
		case 1:
			s1 = color_scalefactors[i].x;
			s2 = color_scalefactors[i].z;
			s3 = color_scalefactors[i].w;
			break;
		case 2:
			s1 = color_scalefactors[i].x;
			s2 = color_scalefactors[i].y;
			s3 = color_scalefactors[i].w;
			break;
		case 3:
			s1 = color_scalefactors[i].x;
			s2 = color_scalefactors[i].y;
			s3 = color_scalefactors[i].z;
			break;
		}

		scalefactors[i] = normalize(float3(s1, s2, s3)) * 1.73205080f;
	}

	float lowparam[4], highparam[4];

	float3 averages[4];
	float3 directions[4];

	line3 lines[4];
	float scale[4];
	float length_squared[4];

	for (i = 0; i < partition_count; i++)
	{
		lowparam[i] = 1e10f;
		highparam[i] = -1e10f;
	}

	compute_averages_and_directions_3_components(pt, blk, ewb, scalefactors, omittedComponent, averages, directions);

	for (i = 0; i < partition_count; i++)
	{
		float3 direc = directions[i];
		if (direc.x + direc.y + direc.z < 0.0f)
			directions[i] = float3(0.0f, 0.0f, 0.0f) - direc;
	}

	for (i = 0; i < partition_count; i++)
	{
		lines[i].a = averages[i];
		if (dot(directions[i], directions[i]) == 0.0f)
			lines[i].b = normalize(float3(1, 1, 1));
		else
			lines[i].b = normalize(directions[i]);
	}

	for (i = 0; i < texels_per_block; i++)
	{
		if (error_weights[i] > 1e-10f)
		{
			int partition = pt->partition_of_texel[i];
			float3 point = float3(data_vr[i], data_vg[i], data_vb[i]) * scalefactors[partition];
			line3 l = lines[partition];
			float param = dot(point - l.a, l.b);
			ei->weights[i] = param;
			if (param < lowparam[partition])
				lowparam[partition] = param;
			if (param > highparam[partition])
				highparam[partition] = param;
		}
		else
		{
			ei->weights[i] = -1e38f;
		}
	}

	float3 lowvalues[4];
	float3 highvalues[4];

	for (i = 0; i < partition_count; i++)
	{
		float length = highparam[i] - lowparam[i];
		if (length < 0)			// case for when none of the texels had any weight
		{
			lowparam[i] = 0.0f;
			highparam[i] = 1e-7f;
		}

		// it is possible for a uniform-color partition to produce length=0; this
		// causes NaN-production and NaN-propagation later on. Set length to
		// a small value to avoid this problem.
		if (length < 1e-7f)
			length = 1e-7f;

		length_squared[i] = length * length;
		scale[i] = 1.0f / length;

		float3 ep0 = lines[i].a + lines[i].b * lowparam[i];
		float3 ep1 = lines[i].a + lines[i].b * highparam[i];

		ep0.x /= scalefactors[i].x;
		ep0.y /= scalefactors[i].y;
		ep0.z /= scalefactors[i].z;
		ep1.x /= scalefactors[i].x;
		ep1.y /= scalefactors[i].y;
		ep1.z /= scalefactors[i].z;

		lowvalues[i] = ep0;
		highvalues[i] = ep1;
	}

	for (i = 0; i < partition_count; i++)
	{
		ei->ep.endpt0[i] = float4(blk->red_min, blk->green_min, blk->blue_min, blk->alpha_min);
		ei->ep.endpt1[i] = float4(blk->red_max, blk->green_max, blk->blue_max, blk->alpha_max);


		float3 ep0 = lowvalues[i];
		float3 ep1 = highvalues[i];

		switch (omittedComponent)
		{
			case 0:
				ei->ep.endpt0[i].y = ep0.x;
				ei->ep.endpt0[i].z = ep0.y;
				ei->ep.endpt0[i].w = ep0.z;
				ei->ep.endpt1[i].y = ep1.x;
				ei->ep.endpt1[i].z = ep1.y;
				ei->ep.endpt1[i].w = ep1.z;
				break;
			case 1:
				ei->ep.endpt0[i].x = ep0.x;
				ei->ep.endpt0[i].z = ep0.y;
				ei->ep.endpt0[i].w = ep0.z;
				ei->ep.endpt1[i].x = ep1.x;
				ei->ep.endpt1[i].z = ep1.y;
				ei->ep.endpt1[i].w = ep1.z;
				break;
			case 2:
				ei->ep.endpt0[i].x = ep0.x;
				ei->ep.endpt0[i].y = ep0.y;
				ei->ep.endpt0[i].w = ep0.z;
				ei->ep.endpt1[i].x = ep1.x;
				ei->ep.endpt1[i].y = ep1.y;
				ei->ep.endpt1[i].w = ep1.z;
				break;
			case 3:
				ei->ep.endpt0[i].x = ep0.x;
				ei->ep.endpt0[i].y = ep0.y;
				ei->ep.endpt0[i].z = ep0.z;
				ei->ep.endpt1[i].x = ep1.x;
				ei->ep.endpt1[i].y = ep1.y;
				ei->ep.endpt1[i].z = ep1.z;
				break;
		}
	}

	for (i = 0; i < texels_per_block; i++)
	{
		int partition = pt->partition_of_texel[i];
		float idx = (ei->weights[i] - lowparam[partition]) * scale[partition];
		if (idx > 1.0f)
			idx = 1.0f;
		else if (!(idx > 0.0f))
			idx = 0.0f;

		ei->weights[i] = idx;
		ei->weight_error_scale[i] = length_squared[partition] * error_weights[i];
		if (astc::isnan(ei->weight_error_scale[i]))
		{
			ASTC_CODEC_INTERNAL_ERROR();
		}
	}

	// print all the data that this function computes.
	#ifdef DEBUG_PRINT_DIAGNOSTICS
		if (print_diagnostics)
		{
			printf("%s: %dx%dx%d texels, %d partitions, component1=%d, component2=%d, component3=%d\n", __func__, xdim, ydim, zdim, partition_count, component1, component2, component3);
			printf("Endpoints:\n");
			for (i = 0; i < partition_count; i++)
			{
				printf("%d Low: <%g %g %f> => <%g %g %g %g>\n", i, lowvalues[i].x, lowvalues[i].y, lowvalues[i].z, ei->ep.endpt0[i].x, ei->ep.endpt0[i].y, ei->ep.endpt0[i].z, ei->ep.endpt0[i].w);
				printf("%d High: <%g %g %g> => <%g %g %g %g>\n", i, highvalues[i].x, highvalues[i].y, highvalues[i].z, ei->ep.endpt1[i].x, ei->ep.endpt1[i].y, ei->ep.endpt1[i].z, ei->ep.endpt1[i].w);
			}
			printf("Ideal-weights:\n");

			for (i = 0; i < texels_per_block; i++)
			{
				printf("%3d <%2d %2d %2d>=> %g (weight=%g)\n", i, (i % xdim), (i / xdim) % ydim, i / (xdim * ydim), ei->weights[i], ei->weight_error_scale[i]);
			}
			printf("\n");
		}
	#endif
}

static void compute_endpoints_and_ideal_weights_rgba(
	const block_size_descriptor* bsd,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	endpoints_and_weights* ei
) {
	int i;

	const float *error_weights = ewb->texel_weight;

	int partition_count = pt->partition_count;
	float lowparam[4], highparam[4];
	for (i = 0; i < partition_count; i++)
	{
		lowparam[i] = 1e10;
		highparam[i] = -1e10;
	}

	float4 averages[4];
	float4 directions_rgba[4];

	line4 lines[4];

	float scale[4];
	float length_squared[4];

	float4 error_weightings[4];
	float4 color_scalefactors[4];
	float4 scalefactors[4];

	int texels_per_block = bsd->texel_count;

	compute_partition_error_color_weightings(bsd, ewb, pt, error_weightings, color_scalefactors);

	for (i = 0; i < partition_count; i++)
		scalefactors[i] = normalize(color_scalefactors[i]) * 2.0f;

	compute_averages_and_directions_rgba(pt, blk, ewb, scalefactors, averages, directions_rgba);

	// if the direction-vector ends up pointing from light to dark, FLIP IT!
	// this will make the first endpoint the darkest one.
	for (i = 0; i < partition_count; i++)
	{
		float4 direc = directions_rgba[i];
		if (direc.x + direc.y + direc.z < 0.0f)
			directions_rgba[i] = float4(0, 0, 0, 0) - direc;
	}

	for (i = 0; i < partition_count; i++)
	{
		lines[i].a = averages[i];
		if (dot(directions_rgba[i], directions_rgba[i]) == 0.0f)
			lines[i].b = normalize(float4(1, 1, 1, 1));
		else
			lines[i].b = normalize(directions_rgba[i]);
	}

	#ifdef DEBUG_PRINT_DIAGNOSTICS
		if (print_diagnostics)
		{
			for (i = 0; i < partition_count; i++)
			{
				printf("Direction-vector %d: <%f %f %f %f>\n", i, directions_rgba[i].x, directions_rgba[i].y, directions_rgba[i].z, directions_rgba[i].w);
				printf("Line %d A: <%f %f %f %f>\n", i, lines[i].a.x, lines[i].a.y, lines[i].a.z, lines[i].a.w);
				printf("Line %d B: <%f %f %f %f>\n", i, lines[i].b.x, lines[i].b.y, lines[i].b.z, lines[i].b.w);
				printf("Scalefactors %d: <%f %f %f %f>\n", i, scalefactors[i].x, scalefactors[i].y, scalefactors[i].z, scalefactors[i].w);
			}
		}
	#endif

	for (i = 0; i < texels_per_block; i++)
	{
		if (error_weights[i] > 1e-10f)
		{
			int partition = pt->partition_of_texel[i];

			float4 point = float4(blk->data_r[i], blk->data_g[i], blk->data_b[i], blk->data_a[i]) * scalefactors[partition];
			line4 l = lines[partition];

			float param = dot(point - l.a, l.b);
			ei->weights[i] = param;
			if (param < lowparam[partition])
				lowparam[partition] = param;
			if (param > highparam[partition])
				highparam[partition] = param;
		}
		else
		{
			ei->weights[i] = -1e38f;
		}
	}

	#ifdef DEBUG_PRINT_DIAGNOSTICS
		if (print_diagnostics)
		{
			for (i = 0; i < partition_count; i++)
				printf("Partition %d: Lowparam=%f Highparam=%f\n", i, lowparam[i], highparam[i]);
		}
	#endif

	for (i = 0; i < partition_count; i++)
	{
		float length = highparam[i] - lowparam[i];
		if (length < 0)
		{
			lowparam[i] = 0.0f;
			highparam[i] = 1e-7f;
		}

		// it is possible for a uniform-color partition to produce length=0; this
		// causes NaN-production and NaN-propagation later on. Set length to
		// a small value to avoid this problem.
		if (length < 1e-7f)
			length = 1e-7f;

		length_squared[i] = length * length;
		scale[i] = 1.0f / length;

		float4 ep0 = lines[i].a + lines[i].b * lowparam[i];
		float4 ep1 = lines[i].a + lines[i].b * highparam[i];

		ep0.x /= scalefactors[i].x;
		ep0.y /= scalefactors[i].y;
		ep0.z /= scalefactors[i].z;
		ep0.w /= scalefactors[i].w;

		ep1.x /= scalefactors[i].x;
		ep1.y /= scalefactors[i].y;
		ep1.z /= scalefactors[i].z;
		ep1.w /= scalefactors[i].w;

		ei->ep.endpt0[i] = ep0;
		ei->ep.endpt1[i] = ep1;
	}

	for (i = 0; i < texels_per_block; i++)
	{
		int partition = pt->partition_of_texel[i];
		float idx = (ei->weights[i] - lowparam[partition]) * scale[partition];
		if (idx > 1.0f)
			idx = 1.0f;
		else if (!(idx > 0.0f))
			idx = 0.0f;
		ei->weights[i] = idx;
		ei->weight_error_scale[i] = error_weights[i] * length_squared[partition];
		if (astc::isnan(ei->weight_error_scale[i]))
		{
			ASTC_CODEC_INTERNAL_ERROR();
		}
	}

	// print all the data that this function computes.
	#ifdef DEBUG_PRINT_DIAGNOSTICS
		if (print_diagnostics)
		{
			printf("%s: %dx%dx%d texels, %d partitions\n", __func__, xdim, ydim, zdim, partition_count);
			printf("Endpoints:\n");
			for (i = 0; i < partition_count; i++)
			{
				printf("%d Low: <%g %g %g %g>\n", i, ei->ep.endpt0[i].x, ei->ep.endpt0[i].y, ei->ep.endpt0[i].z, ei->ep.endpt0[i].w);
				printf("%d High: <%g %g %g %g>\n", i, ei->ep.endpt1[i].x, ei->ep.endpt1[i].y, ei->ep.endpt1[i].z, ei->ep.endpt1[i].w);
			}
			printf("\nIdeal-weights:\n");

			for (i = 0; i < texels_per_block; i++)
			{
				printf("%3d <%2d %2d %2d>=> %g (weight=%g)\n", i, i % xdim, (i / xdim) % ydim, i / (xdim * ydim), ei->weights[i], ei->weight_error_scale[i]);
			}
			printf("\n\n");
		}
	#endif

}

/*
	For a given partitioning, compute: for each partition, the ideal endpoint colors;
	these define a color line for the partition. for each pixel, the ideal position of the pixel on the partition's
	color line. for each pixel, the length of the color line.

	These data allow us to assess the error introduced by removing and quantizing the per-pixel weights.
 */
void compute_endpoints_and_ideal_weights_1_plane(
	const block_size_descriptor* bsd,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	endpoints_and_weights* ei
) {
	#ifdef DEBUG_PRINT_DIAGNOSTICS
		if (print_diagnostics)
			printf("%s: texels_per_block=%dx%dx%d\n\n", __func__, xdim, ydim, zdim);
	#endif

	int uses_alpha = imageblock_uses_alpha(blk);
	if (uses_alpha)
	{
		compute_endpoints_and_ideal_weights_rgba(bsd, pt, blk, ewb, ei);
	}
	else
	{
		compute_endpoints_and_ideal_weights_3_components(bsd, pt, blk, ewb, ei, 3);
	}
}

void compute_endpoints_and_ideal_weights_2_planes(
	const block_size_descriptor* bsd,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	int separate_component,
	endpoints_and_weights* ei1,
	endpoints_and_weights* ei2
) {
	#ifdef DEBUG_PRINT_DIAGNOSTICS
		if (print_diagnostics)
			printf("%s: texels_per_block=%dx%dx%d, separate_component=%d\n\n", __func__, xdim, ydim, zdim, separate_component);
	#endif

	int uses_alpha = imageblock_uses_alpha(blk);
	switch (separate_component)
	{
	case 0:					// separate weights for red
		if (uses_alpha == 1)
			compute_endpoints_and_ideal_weights_3_components(bsd, pt, blk, ewb, ei1, 0);
		else
			compute_endpoints_and_ideal_weights_2_components(bsd, pt, blk, ewb, ei1, 1, 2);
		compute_endpoints_and_ideal_weights_1_component(bsd, pt, blk, ewb, ei2, 0);
		break;

	case 1:					// separate weights for green
		if (uses_alpha == 1)
			compute_endpoints_and_ideal_weights_3_components(bsd, pt, blk, ewb, ei1, 1);
		else
			compute_endpoints_and_ideal_weights_2_components(bsd, pt, blk, ewb, ei1, 0, 2);
		compute_endpoints_and_ideal_weights_1_component(bsd, pt, blk, ewb, ei2, 1);
		break;

	case 2:					// separate weights for blue
		if (uses_alpha == 1)
			compute_endpoints_and_ideal_weights_3_components(bsd, pt, blk, ewb, ei1, 2);
		else
			compute_endpoints_and_ideal_weights_2_components(bsd, pt, blk, ewb, ei1, 0, 1);
		compute_endpoints_and_ideal_weights_1_component(bsd, pt, blk, ewb, ei2, 2);
		break;

	case 3:					// separate weights for alpha
		if (uses_alpha == 0)
		{
			ASTC_CODEC_INTERNAL_ERROR();
		}
		compute_endpoints_and_ideal_weights_3_components(bsd, pt, blk, ewb, ei1, 3);

		compute_endpoints_and_ideal_weights_1_component(bsd, pt, blk, ewb, ei2, 3);
		break;
	}
}

/*
   After having computed ideal weights for the case where a weight exists for
   every texel, we want to compute the ideal weights for the case where weights
   exist only for some texels.

   We do this with a steepest-descent grid solver; this works as follows:

   * First, for each actual weight, perform a weighted averaging based on the
     texels affected by the weight.
   * Then, set step size to <some initial value>
   * Then, repeat:
		1: First, compute for each weight how much the error will change
		   if we change the weight by an infinitesimal amount.
		2: This produces a vector that points the direction we should step in.
		   Normalize this vector.
		3: Perform a step
		4: Check if the step actually improved the error. If it did, perform
		   another step in the same direction; repeat until error no longer
		   improves. If the *first* step did not improve error, then we halve
		   the step size.
		5: If the step size dropped down below <some threshold value>,
		   then we quit, else we go back to #1.

   Subroutines: one routine to apply a step and compute the step's effect on
   the error one routine to compute the error change of an infinitesimal
   weight change

   Data structures needed:
   For every decimation pattern, we need:
   * For each weight, a list of <texel, weight> tuples that tell which texels
     the weight influences.
   * For each texel, a list of <texel, weight> tuples that tell which weights
     go into a given texel.
*/

static float compute_value_of_texel_flt(
	int texel_to_get,
	const decimation_table* it,
	const float* weights
) {
	const uint8_t *texel_weights = it->texel_weights[texel_to_get];
	const float *texel_weights_float = it->texel_weights_float[texel_to_get];

	return (weights[texel_weights[0]] * texel_weights_float[0] +
	        weights[texel_weights[1]] * texel_weights_float[1]) +
	       (weights[texel_weights[2]] * texel_weights_float[2] +
	        weights[texel_weights[3]] * texel_weights_float[3]);
}

static inline float compute_error_of_texel(
	const endpoints_and_weights * eai,
	int texel_to_get,
	const decimation_table* it,
	const float *weights
) {
	float current_value = compute_value_of_texel_flt(texel_to_get, it, weights);
	float valuedif = current_value - eai->weights[texel_to_get];
	return valuedif * valuedif * eai->weight_error_scale[texel_to_get];
}

float compute_error_of_weight_set(
	const endpoints_and_weights* eai,
	const decimation_table* it,
	const float* weights
) {
	int texel_count = it->num_texels;
	float error_summa = 0.0;
	for (int i = 0; i < texel_count; i++)
		error_summa += compute_error_of_texel(eai, i, it, weights);
	return error_summa;
}

/*
	Given a complete weight set and a decimation table, try to
	compute the optimal weight set (assuming infinite precision)
	given the selected decimation table.
*/
void compute_ideal_weights_for_decimation_table(
	const endpoints_and_weights* eai,
	const decimation_table* it,
	float* weight_set,
	float* weights
) {
	int texels_per_block = it->num_texels;
	int weight_count = it->num_weights;

	// perform a shortcut in the case of a complete decimation table
	if (texels_per_block == weight_count)
	{
		for (int i = 0; i < it->num_texels; i++)
		{
			int texel = it->weight_texel[i][0];
			weight_set[i] = eai->weights[texel];
			weights[i] = eai->weight_error_scale[texel];
		}
		return;
	}

	// if the shortcut is not available, we will instead compute a simple estimate
	// and perform a single iteration of refinement on that estimate.
	float infilled_weights[MAX_TEXELS_PER_BLOCK];

	// compute an initial average for each weight.
	for (int i = 0; i < weight_count; i++)
	{
		int texel_count = it->weight_num_texels[i];

		float weight_weight = 1e-10f;	// to avoid 0/0 later on
		float initial_weight = 0.0f;
		for (int j = 0; j < texel_count; j++)
		{
			int texel = it->weight_texel[i][j];
			float weight = it->weights_flt[i][j];
			float contrib_weight = weight * eai->weight_error_scale[texel];
			weight_weight += contrib_weight;
			initial_weight += eai->weights[texel] * contrib_weight;
		}

		weights[i] = weight_weight;
		weight_set[i] = initial_weight / weight_weight;	// this is the 0/0 that is to be avoided.
	}

	for (int i = 0; i < texels_per_block; i++)
	{
		const uint8_t *texel_weights = it->texel_weights[i];
		const float *texel_weights_float = it->texel_weights_float[i];
		infilled_weights[i] = (weight_set[texel_weights[0]] * texel_weights_float[0]
		                     + weight_set[texel_weights[1]] * texel_weights_float[1])
		                    + (weight_set[texel_weights[2]] * texel_weights_float[2]
		                     + weight_set[texel_weights[3]] * texel_weights_float[3]);
	}

	constexpr float stepsize = 0.25f;
	constexpr float ch0_scale = 4.0f * (stepsize * stepsize * (1.0f / (TEXEL_WEIGHT_SUM * TEXEL_WEIGHT_SUM)));
	constexpr float ch1_scale = -2.0f * (stepsize * (2.0f / TEXEL_WEIGHT_SUM));
	constexpr float chd_scale = (ch1_scale / ch0_scale) * stepsize;

	for (int i = 0; i < weight_count; i++)
	{
		float weight_val = weight_set[i];

		const uint8_t *weight_texel_ptr = it->weight_texel[i];
		const float *weights_ptr = it->weights_flt[i];

		// compute the two error changes that can occur from perturbing the current index.
		int num_weights = it->weight_num_texels[i];

		float error_change0 = 1e-10f; // done in order to ensure that this value isn't 0, in order to avoid a possible divide by zero later.
		float error_change1 = 0.0f;

		for (int k = 0; k < num_weights; k++)
		{
			uint8_t weight_texel = weight_texel_ptr[k];
			float weights2 = weights_ptr[k];

			float scale = eai->weight_error_scale[weight_texel] * weights2;
			float old_weight = infilled_weights[weight_texel];
			float ideal_weight = eai->weights[weight_texel];

			error_change0 += weights2 * scale;
			error_change1 += (old_weight - ideal_weight) * scale;
		}

		float step = (error_change1 * chd_scale) / error_change0;
		// clamp the step-value.
		if (step < -stepsize )
			step = -stepsize;
		else if(step > stepsize)
			step = stepsize;

		// update the weight
		weight_set[i] = weight_val + step;
	}
}

/*
	For a decimation table, try to compute an optimal weight set, assuming
	that the weights are quantized and subject to a transfer function.

	We do this as follows:
	First, we take the initial weights and quantize them. This is our initial estimate.
	Then, go through the weights one by one; try to perturb then up and down one weight at a
	time; apply any perturbations that improve overall error
	Repeat until we have made a complete processing pass over all weights without
	triggering any perturbations *OR* we have run 4 full passes.
*/
void compute_ideal_quantized_weights_for_decimation_table(
	const decimation_table* it,
	float low_bound,
	float high_bound,
	const float* weight_set_in,
	float* weight_set_out,
	uint8_t* quantized_weight_set,
	int quantization_level
) {
	int weight_count = it->num_weights;
	const quantization_and_transfer_table *qat = &(quant_and_xfer_tables[quantization_level]);

	#ifdef DEBUG_PRINT_DIAGNOSTICS
		if (print_diagnostics)
		{
			printf("%s : texels-per-block=%d,  weights=%d,  quantization-level=%d\n\n", __func__, texels_per_block, weight_count, quantization_level);

			printf("Weight values before quantization:\n");
			for (i = 0; i < weight_count; i++)
				printf("%3d : %g\n", i, weight_set_in[i]);

			printf("Low-bound: %f  High-bound: %f\n", low_bound, high_bound);
		}
	#endif

	static const int quant_levels[12] = { 2,3,4,5,6,8,10,12,16,20,24,32 };
	float quant_level_m1 = (float)(quant_levels[quantization_level] - 1);

	// Quantize the weight set using both the specified low/high bounds
	// and the standard 0..1 weight bounds.

	/*
	   TODO: WTF issue that we need to examine some time
	*/
	if (!((high_bound - low_bound) > 0.5f))
	{
		low_bound = 0.0f;
		high_bound = 1.0f;
	}

	float rscale = high_bound - low_bound;
	float scale = 1.0f / rscale;

	float scaled_low_bound = low_bound * scale;
	rscale *= 1.0f / 64.0f;

	int i = 0;

#if ASTC_AVX >= 2
	//  TODO: This is currently 4-wide. Could try 8?
	int clipped_weight_count = weight_count & ~3;
	__m128i shuf = _mm_set_epi8(-1, -1, -1, -1, -1, -1, -1, -1,
	                            -1, -1, -1, -1, 12,  8,  4,  0);
	__m128 scalev = _mm_set1_ps(scale);
	__m128 scaled_low_boundv = _mm_set1_ps(scaled_low_bound);
	for (/* Loop vector */; i < clipped_weight_count; i += 4)
	{
		__m128 ix = _mm_load_ps(&weight_set_in[i]);
		ix = _mm_mul_ps(ix, scalev);
		ix = _mm_sub_ps(ix, scaled_low_boundv);

		ix = _mm_max_ps(ix, _mm_setzero_ps());
		ix = _mm_min_ps(ix, _mm_set1_ps(1.0f));

		__m128 ix1 = _mm_mul_ps(ix, _mm_set1_ps(quant_level_m1));
		__m128i weight = _mm_cvtps_epi32(ix1);
		__m128 ixl = _mm_i32gather_ps(qat->unquantized_value_unsc, weight, 4);

		__m128i weight1 = _mm_add_epi32(weight, _mm_set1_epi32(1));
		__m128 ixh = _mm_i32gather_ps(qat->unquantized_value_unsc, weight1, 4);

		__m128 lhs = _mm_add_ps(ixl, ixh);
		__m128 rhs = _mm_mul_ps(ix, _mm_set1_ps(128.0f));
		__m128i mask = _mm_castps_si128(_mm_cmplt_ps(lhs, rhs));
		weight = _mm_blendv_epi8(weight, weight1, mask);
		ixl = _mm_blendv_ps(ixl, ixh, _mm_castsi128_ps(mask));

		// Invert the weight-scaling that was done initially
		__m128 wso = _mm_mul_ps(ixl, _mm_set1_ps(rscale));
		wso = _mm_add_ps(wso, _mm_set1_ps(low_bound));
		_mm_storeu_ps(&weight_set_out[i], wso);

		__m128i scm = _mm_i32gather_epi32(qat->scramble_map, weight, 4);
		__m128i scn = _mm_shuffle_epi8(scm, shuf);

		// This is a hack because _mm_storeu_si32 is still not implemented ...
		_mm_store_ss((float*)&quantized_weight_set[i], _mm_castsi128_ps(scn));
	}
#endif

	for (/*Loop tail */; i < weight_count; i++)
	{
		float ix = (weight_set_in[i] * scale) - scaled_low_bound;
		if (ix < 0.0f)
			ix = 0.0f;
		if (ix > 1.0f) // upper bound must be smaller than 1 to avoid an array overflow below.
			ix = 1.0f;

		// look up the two closest indexes and return the one that was closest.
		float ix1 = ix * quant_level_m1;
		int weight = (int)ix1;
		float ixl = qat->unquantized_value_unsc[weight];
		float ixh = qat->unquantized_value_unsc[weight + 1];

		if (ixl + ixh < 128.0f * ix)
		{
			weight++;
			ixl = ixh;
		}

		// Invert the weight-scaling that was done initially
		weight_set_out[i] = (ixl * rscale) + low_bound;
		quantized_weight_set[i] = (uint8_t)qat->scramble_map[weight];
	}
}

static inline float mat_square_sum(mat2 p)
{
	float a = p.v[0].x;
	float b = p.v[0].y;
	float c = p.v[1].x;
	float d = p.v[1].y;
	return a * a + b * b + c * c + d * d;
}

/* for a given weight set, we wish to recompute the colors so that they are optimal for a particular weight set. */
void recompute_ideal_colors(
	const block_size_descriptor* bsd,
	int weight_quantization_mode,
	endpoints* ep,	// contains the endpoints we wish to update
	float4* rgbs_vectors,	// used to return RGBS-vectors for endpoint mode #6
	float4* rgbo_vectors,	// used to return RGBS-vectors for endpoint mode #7
	const uint8_t* weight_set8,	// the current set of weight values
	const uint8_t* plane2_weight_set8,	// NULL if plane 2 is not actually used.
	int plane2_color_component,	// color component for 2nd plane of weights; -1 if the 2nd plane of weights is not present
	const partition_info* pi,
	const decimation_table* it,
	const imageblock* pb,	// picture-block containing the actual data.
	const error_weight_block* ewb
) {
	int texels_per_block = bsd->texel_count;

	const quantization_and_transfer_table *qat = &(quant_and_xfer_tables[weight_quantization_mode]);

	float weight_set[MAX_WEIGHTS_PER_BLOCK];
	float plane2_weight_set[MAX_WEIGHTS_PER_BLOCK];

	for (int i = 0; i < it->num_weights; i++)
	{
		weight_set[i] = qat->unquantized_value[weight_set8[i]] * (1.0f / 64.0f);;
	}

	if (plane2_weight_set8)
	{
		for (int i = 0; i < it->num_weights; i++)
			plane2_weight_set[i] = qat->unquantized_value[plane2_weight_set8[i]] * (1.0f / 64.0f);
	}

	int partition_count = pi->partition_count;

	#ifdef DEBUG_PRINT_DIAGNOSTICS
		if (print_diagnostics)
		{
			printf("%s : %dx%dx%d texels_per_block, %d partitions, plane2-color-component=%d\n\n", __func__, xdim, ydim, zdim, partition_count, plane2_color_component);

			printf("Pre-adjustment endpoint-colors: \n");
			for (i = 0; i < partition_count; i++)
			{
				printf("%d Low  <%g %g %g %g>\n", i, ep->endpt0[i].x, ep->endpt0[i].y, ep->endpt0[i].z, ep->endpt0[i].w);
				printf("%d High <%g %g %g %g>\n", i, ep->endpt1[i].x, ep->endpt1[i].y, ep->endpt1[i].z, ep->endpt1[i].w);
			}
		}
	#endif

	mat2 pmat1_red[4], pmat1_green[4], pmat1_blue[4], pmat1_alpha[4], pmat1_scale[4];	// matrices for plane of weights 1
	mat2 pmat2_red[4], pmat2_green[4], pmat2_blue[4], pmat2_alpha[4];	// matrices for plane of weights 2
	float2 red_vec[4];
	float2 green_vec[4];
	float2 blue_vec[4];
	float2 alpha_vec[4];
	float2 scale_vec[4];

	for (int i = 0; i < partition_count; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			pmat1_red[i].v[j] = float2(0, 0);
			pmat2_red[i].v[j] = float2(0, 0);
			pmat1_green[i].v[j] = float2(0, 0);
			pmat2_green[i].v[j] = float2(0, 0);
			pmat1_blue[i].v[j] = float2(0, 0);
			pmat2_blue[i].v[j] = float2(0, 0);
			pmat1_alpha[i].v[j] = float2(0, 0);
			pmat2_alpha[i].v[j] = float2(0, 0);
			pmat1_scale[i].v[j] = float2(0, 0);
		}

		red_vec[i] = float2(0, 0);
		green_vec[i] = float2(0, 0);
		blue_vec[i] = float2(0, 0);
		alpha_vec[i] = float2(0, 0);
		scale_vec[i] = float2(0, 0);
	}

	float wmin1[4], wmax1[4];
	float wmin2[4], wmax2[4];
	float red_weight_sum[4];
	float green_weight_sum[4];
	float blue_weight_sum[4];
	float alpha_weight_sum[4];
	float scale_weight_sum[4];

	float red_weight_weight_sum[4];
	float green_weight_weight_sum[4];
	float blue_weight_weight_sum[4];

	float psum[4];				// sum of (weight * qweight^2) across (red,green,blue)
	float qsum[4];				// sum of (weight * qweight * texelval) across (red,green,blue)

	for (int i = 0; i < partition_count; i++)
	{
		wmin1[i] = 1.0f;
		wmax1[i] = 0.0f;
		wmin2[i] = 1.0f;
		wmax2[i] = 0.0f;
		red_weight_sum[i] = 1e-17f;
		green_weight_sum[i] = 1e-17f;
		blue_weight_sum[i] = 1e-17f;
		alpha_weight_sum[i] = 1e-17f;

		scale_weight_sum[i] = 1e-17f;

		red_weight_weight_sum[i] = 1e-17f;
		green_weight_weight_sum[i] = 1e-17f;
		blue_weight_weight_sum[i] = 1e-17f;

		psum[i] = 1e-17f;
		qsum[i] = 1e-17f;
	}

	// for each partition, compute the direction that an RGB-scale color endpoint pair would have.
	float3 rgb_sum[4];
	float3 rgb_weight_sum[4];
	float3 scale_directions[4];
	float scale_min[4];
	float scale_max[4];

	for (int i = 0; i < partition_count; i++)
	{
		rgb_sum[i] = float3(1e-17f, 1e-17f, 1e-17f);
		rgb_weight_sum[i] = float3(1e-17f, 1e-17f, 1e-17f);
	}

	for (int i = 0; i < texels_per_block; i++)
	{
		float3 rgb = float3(pb->data_r[i], pb->data_g[i], pb->data_b[i]);
		float3 rgb_weight = float3(ewb->texel_weight_r[i],
		                           ewb->texel_weight_g[i],
		                           ewb->texel_weight_b[i]);

		int part = pi->partition_of_texel[i];
		rgb_sum[part] = rgb_sum[part] + (rgb * rgb_weight);
		rgb_weight_sum[part] = rgb_weight_sum[part] + rgb_weight;
	}

	for (int i = 0; i < partition_count; i++)
	{
		float3 tmp = float3(rgb_sum[i].x / rgb_weight_sum[i].x,
		                    rgb_sum[i].y / rgb_weight_sum[i].y,
		                    rgb_sum[i].z / rgb_weight_sum[i].z);
		scale_directions[i] = normalize(tmp);
		scale_max[i] = 0.0f;
		scale_min[i] = 1e10f;
	}

	for (int i = 0; i < texels_per_block; i++)
	{
		float r = pb->data_r[i];
		float g = pb->data_g[i];
		float b = pb->data_b[i];
		float a = pb->data_a[i];

		int part = pi->partition_of_texel[i];
		float idx0 = it ? compute_value_of_texel_flt(i, it, weight_set) : weight_set[i];
		float om_idx0 = 1.0f - idx0;

		if (idx0 > wmax1[part])
			wmax1[part] = idx0;
		if (idx0 < wmin1[part])
			wmin1[part] = idx0;

		float red_weight = ewb->texel_weight_r[i];
		float green_weight = ewb->texel_weight_g[i];
		float blue_weight = ewb->texel_weight_b[i];
		float alpha_weight = ewb->texel_weight_a[i];

		float scale_weight = (red_weight + green_weight + blue_weight);

		float3 scale_direction = scale_directions[part];
		float scale = dot(scale_direction, float3(r, g, b));
		if (scale < scale_min[part])
			scale_min[part] = scale;
		if (scale > scale_max[part])
			scale_max[part] = scale;

		red_weight_sum[part] += red_weight;
		green_weight_sum[part] += green_weight;
		blue_weight_sum[part] += blue_weight;
		alpha_weight_sum[part] += alpha_weight;
		scale_weight_sum[part] += scale_weight;

		pmat1_red[part].v[0].x += om_idx0 * om_idx0 * red_weight;
		pmat1_red[part].v[0].y += idx0 * om_idx0 * red_weight;
		pmat1_red[part].v[1].x += idx0 * om_idx0 * red_weight;
		pmat1_red[part].v[1].y += idx0 * idx0 * red_weight;

		pmat1_green[part].v[0].x += om_idx0 * om_idx0 * green_weight;
		pmat1_green[part].v[0].y += idx0 * om_idx0 * green_weight;
		pmat1_green[part].v[1].x += idx0 * om_idx0 * green_weight;
		pmat1_green[part].v[1].y += idx0 * idx0 * green_weight;

		pmat1_blue[part].v[0].x += om_idx0 * om_idx0 * blue_weight;
		pmat1_blue[part].v[0].y += idx0 * om_idx0 * blue_weight;
		pmat1_blue[part].v[1].x += idx0 * om_idx0 * blue_weight;
		pmat1_blue[part].v[1].y += idx0 * idx0 * blue_weight;

		pmat1_alpha[part].v[0].x += om_idx0 * om_idx0 * alpha_weight;
		pmat1_alpha[part].v[0].y += idx0 * om_idx0 * alpha_weight;
		pmat1_alpha[part].v[1].x += idx0 * om_idx0 * alpha_weight;
		pmat1_alpha[part].v[1].y += idx0 * idx0 * alpha_weight;

		pmat1_scale[part].v[0].x += om_idx0 * om_idx0 * scale_weight;
		pmat1_scale[part].v[0].y += idx0 * om_idx0 * scale_weight;
		pmat1_scale[part].v[1].x += idx0 * om_idx0 * scale_weight;
		pmat1_scale[part].v[1].y += idx0 * idx0 * scale_weight;

		float idx1 = 0.0f, om_idx1 = 0.0f;
		if (plane2_weight_set8)
		{
			idx1 = it ? compute_value_of_texel_flt(i, it, plane2_weight_set) : plane2_weight_set[i];
			om_idx1 = 1.0f - idx1;
			if (idx1 > wmax2[part])
				wmax2[part] = idx1;
			if (idx1 < wmin2[part])
				wmin2[part] = idx1;

			pmat2_red[part].v[0].x += om_idx1 * om_idx1 * red_weight;
			pmat2_red[part].v[0].y += idx1 * om_idx1 * red_weight;
			pmat2_red[part].v[1].x += idx1 * om_idx1 * red_weight;
			pmat2_red[part].v[1].y += idx1 * idx1 * red_weight;

			pmat2_green[part].v[0].x += om_idx1 * om_idx1 * green_weight;
			pmat2_green[part].v[0].y += idx1 * om_idx1 * green_weight;
			pmat2_green[part].v[1].x += idx1 * om_idx1 * green_weight;
			pmat2_green[part].v[1].y += idx1 * idx1 * green_weight;

			pmat2_blue[part].v[0].x += om_idx1 * om_idx1 * blue_weight;
			pmat2_blue[part].v[0].y += idx1 * om_idx1 * blue_weight;
			pmat2_blue[part].v[1].x += idx1 * om_idx1 * blue_weight;
			pmat2_blue[part].v[1].y += idx1 * idx1 * blue_weight;

			pmat2_alpha[part].v[0].x += om_idx1 * om_idx1 * alpha_weight;
			pmat2_alpha[part].v[0].y += idx1 * om_idx1 * alpha_weight;
			pmat2_alpha[part].v[1].x += idx1 * om_idx1 * alpha_weight;
			pmat2_alpha[part].v[1].y += idx1 * idx1 * alpha_weight;
		}

		float red_idx = (plane2_color_component == 0) ? idx1 : idx0;
		float green_idx = (plane2_color_component == 1) ? idx1 : idx0;
		float blue_idx = (plane2_color_component == 2) ? idx1 : idx0;
		float alpha_idx = (plane2_color_component == 3) ? idx1 : idx0;

		red_vec[part].x += (red_weight * r) * (1.0f - red_idx);
		green_vec[part].x += (green_weight * g) * (1.0f - green_idx);
		blue_vec[part].x += (blue_weight * b) * (1.0f - blue_idx);
		alpha_vec[part].x += (alpha_weight * a) * (1.0f - alpha_idx);
		scale_vec[part].x += (scale_weight * scale) * om_idx0;

		red_vec[part].y += (red_weight * r) * red_idx;
		green_vec[part].y += (green_weight * g) * green_idx;
		blue_vec[part].y += (blue_weight * b) * blue_idx;
		alpha_vec[part].y += (alpha_weight * a) * alpha_idx;
		scale_vec[part].y += (scale_weight * scale) * idx0;

		red_weight_weight_sum[part] += red_weight * red_idx;
		green_weight_weight_sum[part] += green_weight * green_idx;
		blue_weight_weight_sum[part] += blue_weight * blue_idx;

		psum[part] += red_weight * red_idx * red_idx + green_weight * green_idx * green_idx + blue_weight * blue_idx * blue_idx;
	}

	// calculations specific to mode #7, the HDR RGB-scale mode.
	float red_sum[4];
	float green_sum[4];
	float blue_sum[4];

	for (int i = 0; i < partition_count; i++)
	{
		red_sum[i] = red_vec[i].x + red_vec[i].y;
		green_sum[i] = green_vec[i].x + green_vec[i].y;
		blue_sum[i] = blue_vec[i].x + blue_vec[i].y;
		qsum[i] = red_vec[i].y + green_vec[i].y + blue_vec[i].y;
	}

	// RGB+offset for HDR endpoint mode #7
	int rgbo_fail[4];

	for (int i = 0; i < partition_count; i++)
	{
		mat4 mod7_mat;
		mod7_mat.v[0] = float4(red_weight_sum[i], 0.0f, 0.0f, red_weight_weight_sum[i]);
		mod7_mat.v[1] = float4(0.0f, green_weight_sum[i], 0.0f, green_weight_weight_sum[i]);
		mod7_mat.v[2] = float4(0.0f, 0.0f, blue_weight_sum[i], blue_weight_weight_sum[i]);
		mod7_mat.v[3] = float4(red_weight_weight_sum[i], green_weight_weight_sum[i], blue_weight_weight_sum[i], psum[i]);

		float4 vect = float4(red_sum[i], green_sum[i], blue_sum[i], qsum[i]);

		#ifdef DEBUG_CAPTURE_NAN
			fedisableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif

		mat4 rmod7_mat = invert(mod7_mat);
		float4 rgbovec = transform(rmod7_mat, vect);
		rgbo_vectors[i] = rgbovec;

		// we will occasionally get a failure due to a singular matrix. Record whether such a
		// failure has taken place; if it did, compute rgbo_vectors[] with a different method
		// later on.
		float chkval = dot(rgbovec, rgbovec);
		rgbo_fail[i] = chkval != chkval;

		#ifdef DEBUG_CAPTURE_NAN
			feenableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif
	}

	// initialize the luminance and scale vectors with a reasonable default,
	// just in case the subsequent calculation blows up.
	for (int i = 0; i < partition_count; i++)
	{
		#ifdef DEBUG_CAPTURE_NAN
			fedisableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif

		float scalediv = scale_min[i] / scale_max[i];
		if (!(scalediv > 0.0f))
			scalediv = 0.0f;	// set to zero if scalediv is zero, negative, or NaN.

		#ifdef DEBUG_CAPTURE_NAN
			feenableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif

		if (scalediv > 1.0f)
			scalediv = 1.0f;

		rgbs_vectors[i] = float4(scale_directions[i].x * scale_max[i],
		                         scale_directions[i].y * scale_max[i],
		                         scale_directions[i].z * scale_max[i],
		                         scalediv);
	}

	for (int i = 0; i < partition_count; i++)
	{
		if (wmin1[i] >= wmax1[i] * 0.999f)
		{
			// if all weights in the partition were equal, then just take average
			// of all colors in the partition and use that as both endpoint colors.
			float4 avg = float4((red_vec[i].x + red_vec[i].y) / red_weight_sum[i],
								(green_vec[i].x + green_vec[i].y) / green_weight_sum[i],
								(blue_vec[i].x + blue_vec[i].y) / blue_weight_sum[i],
								(alpha_vec[i].x + alpha_vec[i].y) / alpha_weight_sum[i]);

			if (plane2_color_component != 0 && avg.x == avg.x)
				ep->endpt0[i].x = ep->endpt1[i].x = avg.x;
			if (plane2_color_component != 1 && avg.y == avg.y)
				ep->endpt0[i].y = ep->endpt1[i].y = avg.y;
			if (plane2_color_component != 2 && avg.z == avg.z)
				ep->endpt0[i].z = ep->endpt1[i].z = avg.z;
			if (plane2_color_component != 3 && avg.w == avg.w)
				ep->endpt0[i].w = ep->endpt1[i].w = avg.w;

			rgbs_vectors[i] = float4(scale_directions[i].x * scale_max[i],
			                         scale_directions[i].y * scale_max[i],
			                         scale_directions[i].z * scale_max[i],
			                         1.0f);
		}
		else
		{
			// otherwise, complete the analytic calculation of ideal-endpoint-values
			// for the given set of texel weights and pixel colors.

			#ifdef DEBUG_CAPTURE_NAN
				fedisableexcept(FE_DIVBYZERO | FE_INVALID);
			#endif

			float red_det1 = determinant(pmat1_red[i]);
			float green_det1 = determinant(pmat1_green[i]);
			float blue_det1 = determinant(pmat1_blue[i]);
			float alpha_det1 = determinant(pmat1_alpha[i]);
			float scale_det1 = determinant(pmat1_scale[i]);

			float red_mss1 = mat_square_sum(pmat1_red[i]);
			float green_mss1 = mat_square_sum(pmat1_green[i]);
			float blue_mss1 = mat_square_sum(pmat1_blue[i]);
			float alpha_mss1 = mat_square_sum(pmat1_alpha[i]);
			float scale_mss1 = mat_square_sum(pmat1_scale[i]);

			#ifdef DEBUG_PRINT_DIAGNOSTICS
				if (print_diagnostics)
					printf("Plane-1 partition %d determinants: R=%g G=%g B=%g A=%g S=%g\n", i, red_det1, green_det1, blue_det1, alpha_det1, scale_det1);
			#endif

			pmat1_red[i] = invert(pmat1_red[i]);
			pmat1_green[i] = invert(pmat1_green[i]);
			pmat1_blue[i] = invert(pmat1_blue[i]);
			pmat1_alpha[i] = invert(pmat1_alpha[i]);
			pmat1_scale[i] = invert(pmat1_scale[i]);

			float4 ep0 = float4(dot(pmat1_red[i].v[0], red_vec[i]),
								dot(pmat1_green[i].v[0], green_vec[i]),
								dot(pmat1_blue[i].v[0], blue_vec[i]),
								dot(pmat1_alpha[i].v[0], alpha_vec[i]));
			float4 ep1 = float4(dot(pmat1_red[i].v[1], red_vec[i]),
								dot(pmat1_green[i].v[1], green_vec[i]),
								dot(pmat1_blue[i].v[1], blue_vec[i]),
								dot(pmat1_alpha[i].v[1], alpha_vec[i]));

			float scale_ep0 = dot(pmat1_scale[i].v[0], scale_vec[i]);
			float scale_ep1 = dot(pmat1_scale[i].v[1], scale_vec[i]);

			if (plane2_color_component != 0 && fabsf(red_det1) > (red_mss1 * 1e-4f) && ep0.x == ep0.x && ep1.x == ep1.x)
			{
				ep->endpt0[i].x = ep0.x;
				ep->endpt1[i].x = ep1.x;
			}
			if (plane2_color_component != 1 && fabsf(green_det1) > (green_mss1 * 1e-4f) && ep0.y == ep0.y && ep1.y == ep1.y)
			{
				ep->endpt0[i].y = ep0.y;
				ep->endpt1[i].y = ep1.y;
			}
			if (plane2_color_component != 2 && fabsf(blue_det1) > (blue_mss1 * 1e-4f) && ep0.z == ep0.z && ep1.z == ep1.z)
			{
				ep->endpt0[i].z = ep0.z;
				ep->endpt1[i].z = ep1.z;
			}
			if (plane2_color_component != 3 && fabsf(alpha_det1) > (alpha_mss1 * 1e-4f) && ep0.w == ep0.w && ep1.w == ep1.w)
			{
				ep->endpt0[i].w = ep0.w;
				ep->endpt1[i].w = ep1.w;
			}

			if (fabsf(scale_det1) > (scale_mss1 * 1e-4f) && scale_ep0 == scale_ep0 && scale_ep1 == scale_ep1 && scale_ep0 < scale_ep1)
			{
				float scalediv = scale_ep0 / scale_ep1;
				rgbs_vectors[i] =  float4(scale_directions[i].x * scale_ep1,
				                          scale_directions[i].y * scale_ep1,
				                          scale_directions[i].z * scale_ep1,
				                          scalediv);
			}

			#ifdef DEBUG_CAPTURE_NAN
				feenableexcept(FE_DIVBYZERO | FE_INVALID);
			#endif
		}

		if (plane2_weight_set8)
		{
			if (wmin2[i] >= wmax2[i] * 0.999f)
			{
				// if all weights in the partition were equal, then just take average
				// of all colors in the partition and use that as both endpoint colors.
				float4 avg = float4((red_vec[i].x + red_vec[i].y) / red_weight_sum[i],
									(green_vec[i].x + green_vec[i].y) / green_weight_sum[i],
									(blue_vec[i].x + blue_vec[i].y) / blue_weight_sum[i],
									(alpha_vec[i].x + alpha_vec[i].y) / alpha_weight_sum[i]);

				if (plane2_color_component == 0 && avg.x == avg.x)
					ep->endpt0[i].x = ep->endpt1[i].x = avg.x;
				if (plane2_color_component == 1 && avg.y == avg.y)
					ep->endpt0[i].y = ep->endpt1[i].y = avg.y;
				if (plane2_color_component == 2 && avg.z == avg.z)
					ep->endpt0[i].z = ep->endpt1[i].z = avg.z;
				if (plane2_color_component == 3 && avg.w == avg.w)
					ep->endpt0[i].w = ep->endpt1[i].w = avg.w;
			}
			else
			{
				#ifdef DEBUG_CAPTURE_NAN
					fedisableexcept(FE_DIVBYZERO | FE_INVALID);
				#endif

				// otherwise, complete the analytic calculation of ideal-endpoint-values
				// for the given set of texel weights and pixel colors.
				float red_det2 = determinant(pmat2_red[i]);
				float green_det2 = determinant(pmat2_green[i]);
				float blue_det2 = determinant(pmat2_blue[i]);
				float alpha_det2 = determinant(pmat2_alpha[i]);

				float red_mss2 = mat_square_sum(pmat2_red[i]);
				float green_mss2 = mat_square_sum(pmat2_green[i]);
				float blue_mss2 = mat_square_sum(pmat2_blue[i]);
				float alpha_mss2 = mat_square_sum(pmat2_alpha[i]);

				#ifdef DEBUG_PRINT_DIAGNOSTICS
					if (print_diagnostics)
						printf("Plane-2 partition %d determinants: R=%g G=%g B=%g A=%g\n", i, red_det2, green_det2, blue_det2, alpha_det2);
				#endif

				pmat2_red[i] = invert(pmat2_red[i]);
				pmat2_green[i] = invert(pmat2_green[i]);
				pmat2_blue[i] = invert(pmat2_blue[i]);
				pmat2_alpha[i] = invert(pmat2_alpha[i]);
				float4 ep0 = float4(dot(pmat2_red[i].v[0], red_vec[i]),
									dot(pmat2_green[i].v[0], green_vec[i]),
									dot(pmat2_blue[i].v[0], blue_vec[i]),
									dot(pmat2_alpha[i].v[0], alpha_vec[i]));
				float4 ep1 = float4(dot(pmat2_red[i].v[1], red_vec[i]),
									dot(pmat2_green[i].v[1], green_vec[i]),
									dot(pmat2_blue[i].v[1], blue_vec[i]),
									dot(pmat2_alpha[i].v[1], alpha_vec[i]));

				if (plane2_color_component == 0 && fabsf(red_det2) > (red_mss2 * 1e-4f) && ep0.x == ep0.x && ep1.x == ep1.x)
				{
					ep->endpt0[i].x = ep0.x;
					ep->endpt1[i].x = ep1.x;
				}

				if (plane2_color_component == 1 && fabsf(green_det2) > (green_mss2 * 1e-4f) && ep0.y == ep0.y && ep1.y == ep1.y)
				{
					ep->endpt0[i].y = ep0.y;
					ep->endpt1[i].y = ep1.y;
				}

				if (plane2_color_component == 2 && fabsf(blue_det2) > (blue_mss2 * 1e-4f) && ep0.z == ep0.z && ep1.z == ep1.z)
				{
					ep->endpt0[i].z = ep0.z;
					ep->endpt1[i].z = ep1.z;
				}

				if (plane2_color_component == 3 && fabsf(alpha_det2) > (alpha_mss2 * 1e-4f) && ep0.w == ep0.w && ep1.w == ep1.w)
				{
					ep->endpt0[i].w = ep0.w;
					ep->endpt1[i].w = ep1.w;
				}

				#ifdef DEBUG_CAPTURE_NAN
					feenableexcept(FE_DIVBYZERO | FE_INVALID);
				#endif
			}
		}
	}

	// if the calculation of an RGB-offset vector failed, try to compute
	// a somewhat-sensible value anyway
	for (int i = 0; i < partition_count; i++)
	{
		if (rgbo_fail[i])
		{
			float4 v0 = ep->endpt0[i];
			float4 v1 = ep->endpt1[i];
			float avgdif = dot(float3(v1.x - v0.x, v1.y - v0.y, v1.z - v0.z), float3(1, 1, 1)) * (1.0f / 3.0f);
			if (avgdif <= 0.0f)
				avgdif = 0.0f;
			float4 avg = (v0 + v1) * 0.5f;
			float4 ep0 = avg - float4(avgdif, avgdif, avgdif, avgdif) * 0.5f;

			rgbo_vectors[i] = float4(ep0.x, ep0.y, ep0.z, avgdif);
		}
	}

	#ifdef DEBUG_PRINT_DIAGNOSTICS
		if (print_diagnostics)
		{
			printf("Post-adjustment endpoint-colors: \n");
			for (i = 0; i < partition_count; i++)
			{
				printf("%d Low  <%g %g %g %g>\n", i,
				       (double)ep->endpt0[i].x, (double)ep->endpt0[i].y,
				       (double)ep->endpt0[i].z, (double)ep->endpt0[i].w);
				printf("%d High <%g %g %g %g>\n", i,
				       (double)ep->endpt1[i].x, (double)ep->endpt1[i].y,
				       (double)ep->endpt1[i].z, (double)ep->endpt1[i].w);
				printf("%d RGBS: <%g %g %g %g>\n", i,
				       (double)rgbs_vectors[i].x, (double)rgbs_vectors[i].y,
				       (double)rgbs_vectors[i].z, (double)rgbs_vectors[i].w);
				printf("%d RGBO <%g %g %g %g>\n", i,
				       (double)rgbo_vectors[i].x, (double)rgbo_vectors[i].y,
				       (double)rgbo_vectors[i].z, (double)rgbo_vectors[i].w);
			}
		}
	#endif
}
