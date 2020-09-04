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
 * @brief Functions for computing color endpoints and texel weights.
 */

#include <cassert>

#include "astcenc_internal.h"

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
	unsigned int component
) {
	int partition_count = pt->partition_count;
	ei->ep.partition_count = partition_count;

	float lowvalues[4], highvalues[4];
	float partition_error_scale[4];
	float linelengths_rcp[4];

	int texels_per_block = bsd->texel_count;

	const float *error_weights;
	const float* data_vr = nullptr;
	assert(component <= 3);
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
	}

	for (int i = 0; i < partition_count; i++)
	{
		lowvalues[i] = 1e10f;
		highvalues[i] = -1e10f;
	}

	for (int i = 0; i < texels_per_block; i++)
	{
		if (error_weights[i] > 1e-10f)
		{
			float value = data_vr[i];
			int partition = pt->partition_of_texel[i];

			if (value < lowvalues[partition])
			{
				lowvalues[partition] = value;
			}

			if (value > highvalues[partition])
			{
				highvalues[partition] = value;
			}
		}
	}

	for (int i = 0; i < partition_count; i++)
	{
		float diff = highvalues[i] - lowvalues[i];

		if (diff < 0)
		{
			lowvalues[i] = 0.0f;
			highvalues[i] = 0.0f;
		}

		if (diff < 1e-7f)
		{
			diff = 1e-7f;
		}

		partition_error_scale[i] = diff * diff;
		linelengths_rcp[i] = 1.0f / diff;
	}

	for (int i = 0; i < texels_per_block; i++)
	{
		float value = data_vr[i];
		int partition = pt->partition_of_texel[i];
		value -= lowvalues[partition];
		value *= linelengths_rcp[partition];

		if (value > 1.0f)
		{
			value = 1.0f;
		}
		else if (!(value > 0.0f))
		{
			value = 0.0f;
		}

		ei->weights[i] = value;
		ei->weight_error_scale[i] = partition_error_scale[partition] * error_weights[i];
		assert(!astc::isnan(ei->weight_error_scale[i]));
	}

	for (int i = 0; i < partition_count; i++)
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
	else // (component1 == 1 && component2 == 2)
	{
		error_weights = ewb->texel_weight_gb;
		data_vr = blk->data_g;
		data_vg = blk->data_b;
	}

	int texels_per_block = bsd->texel_count;

	compute_partition_error_color_weightings(bsd, ewb, pt, error_weightings, color_scalefactors);

	for (int i = 0; i < partition_count; i++)
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

	for (int i = 0; i < partition_count; i++)
	{
		lowparam[i] = 1e10;
		highparam[i] = -1e10;
	}

	compute_averages_and_directions_2_components(pt, blk, ewb, scalefactors, component1, component2, averages, directions);

	for (int i = 0; i < partition_count; i++)
	{
		float2 egv = directions[i];
		if (egv.x + egv.y < 0.0f)
			directions[i] = float2(0.0f, 0.0f) - egv;
	}

	for (int i = 0; i < partition_count; i++)
	{
		lines[i].a = averages[i];
		if (dot(directions[i], directions[i]) == 0.0f)
		{
			lines[i].b = normalize(float2(1.0f, 1.0f));
		}
		else
		{
			lines[i].b = normalize(directions[i]);
		}
	}

	for (int i = 0; i < texels_per_block; i++)
	{
		if (error_weights[i] > 1e-10f)
		{
			int partition = pt->partition_of_texel[i];
			float2 point = float2(data_vr[i], data_vg[i]) * scalefactors[partition];
			line2 l = lines[partition];
			float param = dot(point - l.a, l.b);
			ei->weights[i] = param;

			if (param < lowparam[partition])
			{
				lowparam[partition] = param;
			}

			if (param > highparam[partition])
			{
				highparam[partition] = param;
			}
		}
		else
		{
			ei->weights[i] = -1e38f;
		}
	}

	float2 lowvalues[4];
	float2 highvalues[4];

	for (int i = 0; i < partition_count; i++)
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
		{
			length = 1e-7f;
		}

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

	for (int i = 0; i < partition_count; i++)
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

	for (int i = 0; i < texels_per_block; i++)
	{
		int partition = pt->partition_of_texel[i];
		float idx = (ei->weights[i] - lowparam[partition]) * scale[partition];
		if (idx > 1.0f)
		{
			idx = 1.0f;
		}
		else if (!(idx > 0.0f))
		{
			idx = 0.0f;
		}

		ei->weights[i] = idx;
		ei->weight_error_scale[i] = length_squared[partition] * error_weights[i];
		assert(!astc::isnan(ei->weight_error_scale[i]));
	}
}

static void compute_endpoints_and_ideal_weights_3_components(
	const block_size_descriptor* bsd,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	endpoints_and_weights* ei,
	int omittedComponent
) {
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

	for (int i = 0; i < partition_count; i++)
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

	for (int i = 0; i < partition_count; i++)
	{
		lowparam[i] = 1e10f;
		highparam[i] = -1e10f;
	}

	compute_averages_and_directions_3_components(pt, blk, ewb, scalefactors, omittedComponent, averages, directions);

	for (int i = 0; i < partition_count; i++)
	{
		float3 direc = directions[i];
		if (direc.x + direc.y + direc.z < 0.0f)
		{
			directions[i] = float3(0.0f, 0.0f, 0.0f) - direc;
		}
	}

	for (int i = 0; i < partition_count; i++)
	{
		lines[i].a = averages[i];
		if (dot(directions[i], directions[i]) == 0.0f)
		{
			lines[i].b = normalize(float3(1.0f, 1.0f, 1.0f));
		}
		else
		{
			lines[i].b = normalize(directions[i]);
		}
	}

	for (int i = 0; i < texels_per_block; i++)
	{
		if (error_weights[i] > 1e-10f)
		{
			int partition = pt->partition_of_texel[i];
			float3 point = float3(data_vr[i], data_vg[i], data_vb[i]) * scalefactors[partition];
			line3 l = lines[partition];
			float param = dot(point - l.a, l.b);
			ei->weights[i] = param;

			if (param < lowparam[partition])
			{
				lowparam[partition] = param;
			}

			if (param > highparam[partition])
			{
				highparam[partition] = param;
			}
		}
		else
		{
			ei->weights[i] = -1e38f;
		}
	}

	float3 lowvalues[4];
	float3 highvalues[4];

	for (int i = 0; i < partition_count; i++)
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
		{
			length = 1e-7f;
		}

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

	for (int i = 0; i < partition_count; i++)
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

	for (int i = 0; i < texels_per_block; i++)
	{
		int partition = pt->partition_of_texel[i];
		float idx = (ei->weights[i] - lowparam[partition]) * scale[partition];
		if (idx > 1.0f)
		{
			idx = 1.0f;
		}
		else if (!(idx > 0.0f))
		{
			idx = 0.0f;
		}

		ei->weights[i] = idx;
		ei->weight_error_scale[i] = length_squared[partition] * error_weights[i];
		assert(!astc::isnan(ei->weight_error_scale[i]));
	}
}

static void compute_endpoints_and_ideal_weights_rgba(
	const block_size_descriptor* bsd,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	endpoints_and_weights* ei
) {
	const float *error_weights = ewb->texel_weight;

	int partition_count = pt->partition_count;
	float lowparam[4], highparam[4];
	for (int i = 0; i < partition_count; i++)
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

	for (int i = 0; i < partition_count; i++)
	{
		scalefactors[i] = normalize(color_scalefactors[i]) * 2.0f;
	}

	compute_averages_and_directions_rgba(pt, blk, ewb, scalefactors, averages, directions_rgba);

	// if the direction-vector ends up pointing from light to dark, FLIP IT!
	// this will make the first endpoint the darkest one.
	for (int i = 0; i < partition_count; i++)
	{
		float4 direc = directions_rgba[i];
		if (direc.x + direc.y + direc.z < 0.0f)
		{
			directions_rgba[i] = float4(0.0f, 0.0f, 0.0f, 0.0f) - direc;
		}
	}

	for (int i = 0; i < partition_count; i++)
	{
		lines[i].a = averages[i];
		if (dot(directions_rgba[i], directions_rgba[i]) == 0.0f)
		{
			lines[i].b = normalize(float4(1.0f, 1.0f, 1.0f, 1.0f));
		}
		else
		{
			lines[i].b = normalize(directions_rgba[i]);
		}
	}

	for (int i = 0; i < texels_per_block; i++)
	{
		if (error_weights[i] > 1e-10f)
		{
			int partition = pt->partition_of_texel[i];

			float4 point = float4(blk->data_r[i], blk->data_g[i], blk->data_b[i], blk->data_a[i]) * scalefactors[partition];
			line4 l = lines[partition];

			float param = dot(point - l.a, l.b);
			ei->weights[i] = param;

			if (param < lowparam[partition])
			{
				lowparam[partition] = param;
			}

			if (param > highparam[partition])
			{
				highparam[partition] = param;
			}
		}
		else
		{
			ei->weights[i] = -1e38f;
		}
	}

	for (int i = 0; i < partition_count; i++)
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
		{
			length = 1e-7f;
		}

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

	for (int i = 0; i < texels_per_block; i++)
	{
		int partition = pt->partition_of_texel[i];
		float idx = (ei->weights[i] - lowparam[partition]) * scale[partition];
		if (idx > 1.0f)
		{
			idx = 1.0f;
		}
		else if (!(idx > 0.0f))
		{
			idx = 0.0f;
		}
		ei->weights[i] = idx;
		ei->weight_error_scale[i] = error_weights[i] * length_squared[partition];
		assert(!astc::isnan(ei->weight_error_scale[i]));
	}
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
	int uses_alpha = imageblock_uses_alpha(blk);
	switch (separate_component)
	{
	case 0:					// separate weights for red
		if (uses_alpha == 1)
		{
			compute_endpoints_and_ideal_weights_3_components(bsd, pt, blk, ewb, ei1, 0);
		}
		else
		{
			compute_endpoints_and_ideal_weights_2_components(bsd, pt, blk, ewb, ei1, 1, 2);
		}
		compute_endpoints_and_ideal_weights_1_component(bsd, pt, blk, ewb, ei2, 0);
		break;

	case 1:					// separate weights for green
		if (uses_alpha == 1)
		{
			compute_endpoints_and_ideal_weights_3_components(bsd, pt, blk, ewb, ei1, 1);
		}
		else
		{
			compute_endpoints_and_ideal_weights_2_components(bsd, pt, blk, ewb, ei1, 0, 2);
		}
		compute_endpoints_and_ideal_weights_1_component(bsd, pt, blk, ewb, ei2, 1);
		break;

	case 2:					// separate weights for blue
		if (uses_alpha == 1)
		{
			compute_endpoints_and_ideal_weights_3_components(bsd, pt, blk, ewb, ei1, 2);
		}
		else
		{
			compute_endpoints_and_ideal_weights_2_components(bsd, pt, blk, ewb, ei1, 0, 1);
		}
		compute_endpoints_and_ideal_weights_1_component(bsd, pt, blk, ewb, ei2, 2);
		break;

	case 3:					// separate weights for alpha
		assert(uses_alpha != 0);
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
	{
		error_summa += compute_error_of_texel(eai, i, it, weights);
	}
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
		if (step < -stepsize)
		{
			step = -stepsize;
		}
		else if (step > stepsize)
		{
			step = stepsize;
		}

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

#if ASTCENC_AVX >= 2
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
		{
			ix = 0.0f;
		}
		if (ix > 1.0f) // upper bound must be smaller than 1 to avoid an array overflow below.
		{
			ix = 1.0f;
		}

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

static inline float4 compute_rgbovec(
	float4 rgba_weight_sum,
	float3 weight_weight_sum,
	float red_sum,
	float green_sum,
	float blue_sum,
	float psum,
	float qsum
) {
	// Compute the rgb+offset for HDR endpoint mode #7. Since the matrix needed
	// has a regular structure, we can simplify the inverse calculation. This
	// gives us ~24 multiplications, down from 96 for a generic inverse

	// mat[0] = float4(rgba_ws.x,      0.0f,      0.0f, wght_ws.x);
	// mat[1] = float4(     0.0f, rgba_ws.y,      0.0f, wght_ws.y);
	// mat[2] = float4(     0.0f,      0.0f, rgba_ws.z, wght_ws.z);
	// mat[3] = float4(wght_ws.x, wght_ws.y, wght_ws.z,      psum);
	// mat = invert(mat);

	float X = rgba_weight_sum.x;
	float Y = rgba_weight_sum.y;
	float Z = rgba_weight_sum.z;
	float P = weight_weight_sum.x;
	float Q = weight_weight_sum.y;
	float R = weight_weight_sum.z;
	float S = psum;

	float PP = P * P;
	float QQ = Q * Q;
	float RR = R * R;

	float SZmRR = S * Z - RR;
	float DT = SZmRR * Y - Z * QQ;
	float YP = Y * P;
	float QX = Q * X;
	float YX = Y * X;
	float mZYP = -Z * YP;
	float mZQX = -Z * QX;
	float mRYX = -R * YX;
	float ZQP = Z * Q * P;
	float RYP = R * YP;
	float RQX = R * QX;

	// Compute the reciprocal of matrix determinant.
	float rdet = 1.0f / (DT * X + mZYP * P);

	// Actually compute the adjugate matrix, not the inverse, and apply the
	// multiplication by 1/det to the vector separately.
	float4 mat0 = float4(DT, ZQP, RYP, mZYP);
	float4 mat1 = float4(ZQP, SZmRR * X - Z * PP, RQX, mZQX);
	float4 mat2 = float4(RYP, RQX, (S * Y - QQ) * X - Y * PP, mRYX);
	float4 mat3 = float4(mZYP, mZQX, mRYX, Z * YX);
	float4 vect = float4(red_sum, green_sum, blue_sum, qsum) * rdet;

	#ifdef DEBUG_CAPTURE_NAN
	    fedisableexcept(FE_DIVBYZERO | FE_INVALID);
	#endif

	return float4(dot(mat0, vect),
	              dot(mat1, vect),
	              dot(mat2, vect),
	              dot(mat3, vect));
}

/* for a given weight set, we wish to recompute the colors so that they are optimal for a particular weight set. */
void recompute_ideal_colors(
	int weight_quantization_mode,
	endpoints* ep,	// contains the endpoints we wish to update
	float4* rgbs_vectors,	// used to return RGBS-vectors for endpoint mode #6
	float4* rgbo_vectors,	// used to return RGBO-vectors for endpoint mode #7
	const uint8_t* weight_set8,	// the current set of weight values
	const uint8_t* plane2_weight_set8,	// nullptr if plane 2 is not actually used.
	int plane2_color_component,	// color component for 2nd plane of weights; -1 if the 2nd plane of weights is not present
	const partition_info* pi,
	const decimation_table* it,
	const imageblock* pb,	// picture-block containing the actual data.
	const error_weight_block* ewb
) {
	const quantization_and_transfer_table *qat = &(quant_and_xfer_tables[weight_quantization_mode]);

	float weight_set[MAX_WEIGHTS_PER_BLOCK];
	float plane2_weight_set[MAX_WEIGHTS_PER_BLOCK];

	for (int i = 0; i < it->num_weights; i++)
	{
		weight_set[i] = qat->unquantized_value[weight_set8[i]] * (1.0f / 64.0f);
	}

	if (plane2_weight_set8)
	{
		for (int i = 0; i < it->num_weights; i++)
		{
			plane2_weight_set[i] = qat->unquantized_value[plane2_weight_set8[i]] * (1.0f / 64.0f);
		}
	}

	int partition_count = pi->partition_count;

	for (int i = 0; i < partition_count; i++)
	{
		float4 rgba_sum        = float4(1e-17f, 1e-17f, 1e-17f, 1e-17f);
		float4 rgba_weight_sum = float4(1e-17f, 1e-17f, 1e-17f, 1e-17f);

		int texelcount = pi->texels_per_partition[i];
		const uint8_t *texel_indexes = pi->texels_of_partition[i];
		for (int j = 0; j < texelcount; j++)
		{
			int tix = texel_indexes[j];

			float4 rgba = float4(pb->data_r[tix], pb->data_g[tix], pb->data_b[tix], pb->data_a[tix]);
			float4 error_weight = float4(ewb->texel_weight_r[tix], ewb->texel_weight_g[tix], ewb->texel_weight_b[tix], ewb->texel_weight_a[tix]);

			rgba_sum = rgba_sum + (rgba * error_weight);
			rgba_weight_sum = rgba_weight_sum + error_weight;
		}

		float3 scale_direction = normalize(float3(
		        rgba_sum.x * (1.0f / rgba_weight_sum.x),
		        rgba_sum.y * (1.0f / rgba_weight_sum.y),
		        rgba_sum.z * (1.0f / rgba_weight_sum.z)));

		float scale_max = 0.0f;
		float scale_min = 1e10f;

		float wmin1 = 1.0f;
		float wmax1 = 0.0f;
		float wmin2 = 1.0f;
		float wmax2 = 0.0f;

		float4 left_sum    = float4(0.0f, 0.0f, 0.0f, 0.0f);
		float4 middle_sum  = float4(0.0f, 0.0f, 0.0f, 0.0f);
		float4 right_sum   = float4(0.0f, 0.0f, 0.0f, 0.0f);

		float4 left2_sum   = float4(0.0f, 0.0f, 0.0f, 0.0f);
		float4 middle2_sum = float4(0.0f, 0.0f, 0.0f, 0.0f);
		float4 right2_sum  = float4(0.0f, 0.0f, 0.0f, 0.0f);

		float3 lmrs_sum = float3(0.0f, 0.0f, 0.0f);

		float4 color_vec_x = float4(0.0f, 0.0f, 0.0f, 0.0f);
		float4 color_vec_y = float4(0.0f, 0.0f, 0.0f, 0.0f);

		float2 scale_vec = float2(0.0f, 0.0f);

		float3 weight_weight_sum = float3(1e-17f, 1e-17f, 1e-17f);
		float psum = 1e-17f;

		// FIXME: the loop below has too many responsibilities, making it inefficient.
		for (int j = 0; j < texelcount; j++)
		{
			int tix = texel_indexes[j];

			float4 rgba = float4(pb->data_r[tix], pb->data_g[tix], pb->data_b[tix], pb->data_a[tix]);
			float4 color_weight = float4(ewb->texel_weight_r[tix], ewb->texel_weight_g[tix], ewb->texel_weight_b[tix], ewb->texel_weight_a[tix]);

			float3 color_weight3 = float3(color_weight.x, color_weight.y, color_weight.z);
			float3 rgb = float3(rgba.x, rgba.y, rgba.z);

			// FIXME: move this calculation out to the color block.
			float ls_weight = (color_weight.x + color_weight.y + color_weight.z);

			const uint8_t *texel_weights = it->texel_weights[tix];
			const float *texel_weights_float = it->texel_weights_float[tix];
			float idx0 = (weight_set[texel_weights[0]] * texel_weights_float[0]
			            + weight_set[texel_weights[1]] * texel_weights_float[1])
			           + (weight_set[texel_weights[2]] * texel_weights_float[2]
			            + weight_set[texel_weights[3]] * texel_weights_float[3]);

			float om_idx0 = 1.0f - idx0;
			if (idx0 > wmax1)
			{
				wmax1 = idx0;
			}

			if (idx0 < wmin1)
			{
				wmin1 = idx0;
			}

			float scale = dot(scale_direction, rgb);
			if (scale < scale_min)
			{
				scale_min = scale;
			}

			if (scale > scale_max)
			{
				scale_max = scale;
			}

			float4 left   = color_weight * (om_idx0 * om_idx0);
			float4 middle = color_weight * (om_idx0 * idx0);
			float4 right  = color_weight * (idx0 * idx0);

			float3 lmrs = float3(om_idx0 * om_idx0,
			                     om_idx0 * idx0,
			                     idx0 * idx0) * ls_weight;

			left_sum   = left_sum + left;
			middle_sum = middle_sum + middle;
			right_sum  = right_sum + right;

			lmrs_sum = lmrs_sum + lmrs;

			float idx1 = 0.0f;
			float om_idx1 = 0.0f;

			if (plane2_weight_set8)
			{
				idx1 = (plane2_weight_set[texel_weights[0]] * texel_weights_float[0]
				      + plane2_weight_set[texel_weights[1]] * texel_weights_float[1])
				     + (plane2_weight_set[texel_weights[2]] * texel_weights_float[2]
				      + plane2_weight_set[texel_weights[3]] * texel_weights_float[3]);

				om_idx1 = 1.0f - idx1;
				if (idx1 > wmax2)
				{
					wmax2 = idx1;
				}

				if (idx1 < wmin2)
				{
					wmin2 = idx1;
				}

				float4 left2   = color_weight * (om_idx1 * om_idx1);
				float4 middle2 = color_weight * (om_idx1 * idx1);
				float4 right2  = color_weight * (idx1 * idx1);

				left2_sum   = left2_sum   + left2;
				middle2_sum = middle2_sum + middle2;
				right2_sum  = right2_sum  + right2;
			}

			float4 color_idx = float4((plane2_color_component == 0) ? idx1 : idx0,
			                          (plane2_color_component == 1) ? idx1 : idx0,
			                          (plane2_color_component == 2) ? idx1 : idx0,
			                          (plane2_color_component == 3) ? idx1 : idx0);

			float3 color_idx3 = float3(color_idx.x, color_idx.y, color_idx.z);

			float4 cwprod = color_weight * rgba;
			float4 cwiprod = cwprod * color_idx;

			color_vec_y = color_vec_y + cwiprod;
			color_vec_x = color_vec_x + (cwprod - cwiprod);

			scale_vec = scale_vec + float2(om_idx0, idx0) * (ls_weight * scale);

			weight_weight_sum = weight_weight_sum + (color_weight3 * color_idx3);

			psum += dot(color_weight3 * color_idx3, color_idx3);
		}

		// calculations specific to mode #7, the HDR RGB-scale mode.
		// FIXME: Can we skip this for LDR textures?
		float red_sum   = color_vec_x.x + color_vec_y.x;
		float green_sum = color_vec_x.y + color_vec_y.y;
		float blue_sum  = color_vec_x.z + color_vec_y.z;
		float qsum = color_vec_y.x + color_vec_y.y + color_vec_y.z;

		#ifdef DEBUG_CAPTURE_NAN
		    fedisableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif

		float4 rgbovec = compute_rgbovec(rgba_weight_sum, weight_weight_sum,
		                                 red_sum, green_sum, blue_sum, psum, qsum);
		rgbo_vectors[i] = rgbovec;

		// We will occasionally get a failure due to the use of a singular
		// (non-invertible) matrix. Record whether such a failure has taken
		// place; if it did, compute rgbo_vectors[] with a different method
		// later on.
		float chkval = dot(rgbovec, rgbovec);
		int rgbo_fail = chkval != chkval;

		// Initialize the luminance and scale vectors with a reasonable
		//  default, just in case the subsequent calculation blows up.
		#ifdef DEBUG_CAPTURE_NAN
			fedisableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif

		float scalediv = scale_min * (1.0f / MAX(scale_max, 1e-10f));
		if (!(scalediv > 0.0f))
		{
			scalediv = 0.0f;    // set to zero if scalediv is negative, or NaN.
		}

		if (scalediv > 1.0f)
		{
			scalediv = 1.0f;
		}

		#ifdef DEBUG_CAPTURE_NAN
			feenableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif

		float3 sds = scale_direction * scale_max;

		rgbs_vectors[i] = float4(sds.x, sds.y, sds.z, scalediv);

		if (wmin1 >= wmax1 * 0.999f)
		{
			// if all weights in the partition were equal, then just take average
			// of all colors in the partition and use that as both endpoint colors.
			float4 avg = (color_vec_x + color_vec_y) *
			             float4(1.0f / rgba_weight_sum.x,
			                    1.0f / rgba_weight_sum.y,
			                    1.0f / rgba_weight_sum.z,
			                    1.0f / rgba_weight_sum.w);

			if (plane2_color_component != 0 && avg.x == avg.x)
			{
				ep->endpt0[i].x = ep->endpt1[i].x = avg.x;
			}

			if (plane2_color_component != 1 && avg.y == avg.y)
			{
				ep->endpt0[i].y = ep->endpt1[i].y = avg.y;
			}

			if (plane2_color_component != 2 && avg.z == avg.z)
			{
				ep->endpt0[i].z = ep->endpt1[i].z = avg.z;
			}

			if (plane2_color_component != 3 && avg.w == avg.w)
			{
				ep->endpt0[i].w = ep->endpt1[i].w = avg.w;
			}

			rgbs_vectors[i] = float4( sds.x, sds.y, sds.z, 1.0f);
		}
		else
		{
			// otherwise, complete the analytic calculation of ideal-endpoint-values
			// for the given set of texel weights and pixel colors.

			#ifdef DEBUG_CAPTURE_NAN
			    fedisableexcept(FE_DIVBYZERO | FE_INVALID);
			#endif

			float4 color_det1 = (left_sum * right_sum) - (middle_sum * middle_sum);
			float4 color_rdet1 = float4(1.0f / color_det1.x,
			                            1.0f / color_det1.y,
			                            1.0f / color_det1.z,
			                            1.0f / color_det1.w );

			float ls_det1  = (lmrs_sum.x * lmrs_sum.z) - (lmrs_sum.y * lmrs_sum.y);
			float ls_rdet1 = 1.0f / ls_det1;

			float4 color_mss1 = (left_sum * left_sum)
			                  + (2.0f * middle_sum * middle_sum)
			                  + (right_sum * right_sum);

			float ls_mss1 = (lmrs_sum.x * lmrs_sum.x)
			              + (2.0f * lmrs_sum.y * lmrs_sum.y)
			              + (lmrs_sum.z * lmrs_sum.z);

			float4 ep0 = (right_sum * color_vec_x - middle_sum * color_vec_y) * color_rdet1;
			float4 ep1 = (left_sum * color_vec_y - middle_sum * color_vec_x) * color_rdet1;

			float scale_ep0 = (lmrs_sum.z * scale_vec.x - lmrs_sum.y * scale_vec.y) * ls_rdet1;
			float scale_ep1 = (lmrs_sum.x * scale_vec.y - lmrs_sum.y * scale_vec.x) * ls_rdet1;

			if (plane2_color_component != 0 && fabsf(color_det1.x) > (color_mss1.x * 1e-4f) && ep0.x == ep0.x && ep1.x == ep1.x)
			{
				ep->endpt0[i].x = ep0.x;
				ep->endpt1[i].x = ep1.x;
			}

			if (plane2_color_component != 1 && fabsf(color_det1.y) > (color_mss1.y * 1e-4f) && ep0.y == ep0.y && ep1.y == ep1.y)
			{
				ep->endpt0[i].y = ep0.y;
				ep->endpt1[i].y = ep1.y;
			}

			if (plane2_color_component != 2 && fabsf(color_det1.z) > (color_mss1.z * 1e-4f) && ep0.z == ep0.z && ep1.z == ep1.z)
			{
				ep->endpt0[i].z = ep0.z;
				ep->endpt1[i].z = ep1.z;
			}

			if (plane2_color_component != 3 && fabsf(color_det1.w) > (color_mss1.w * 1e-4f) && ep0.w == ep0.w && ep1.w == ep1.w)
			{
				ep->endpt0[i].w = ep0.w;
				ep->endpt1[i].w = ep1.w;
			}

			if (fabsf(ls_det1) > (ls_mss1 * 1e-4f) && scale_ep0 == scale_ep0 && scale_ep1 == scale_ep1 && scale_ep0 < scale_ep1)
			{
				float scalediv2 = scale_ep0 * (1.0f / scale_ep1);
				float3 sdsm = scale_direction * scale_ep1;
				rgbs_vectors[i] = float4(sdsm.x, sdsm.y, sdsm.z, scalediv2);
			}

			#ifdef DEBUG_CAPTURE_NAN
				feenableexcept(FE_DIVBYZERO | FE_INVALID);
			#endif
		}

		if (plane2_weight_set8)
		{
			if (wmin2 >= wmax2 * 0.999f)
			{
				// if all weights in the partition were equal, then just take average
				// of all colors in the partition and use that as both endpoint colors.
				float4 avg = (color_vec_x + color_vec_y)
				           * float4(1.0f / rgba_weight_sum.x,
				                    1.0f / rgba_weight_sum.y,
				                    1.0f / rgba_weight_sum.z,
				                    1.0f / rgba_weight_sum.w);

				if (plane2_color_component == 0 && avg.x == avg.x)
				{
					ep->endpt0[i].x = ep->endpt1[i].x = avg.x;
				}

				if (plane2_color_component == 1 && avg.y == avg.y)
				{
					ep->endpt0[i].y = ep->endpt1[i].y = avg.y;
				}

				if (plane2_color_component == 2 && avg.z == avg.z)
				{
					ep->endpt0[i].z = ep->endpt1[i].z = avg.z;
				}

				if (plane2_color_component == 3 && avg.w == avg.w)
				{
					ep->endpt0[i].w = ep->endpt1[i].w = avg.w;
				}
			}
			else
			{
				#ifdef DEBUG_CAPTURE_NAN
					fedisableexcept(FE_DIVBYZERO | FE_INVALID);
				#endif

				// otherwise, complete the analytic calculation of ideal-endpoint-values
				// for the given set of texel weigths and pixel colors.
				float4 color_det2 = (left2_sum * right2_sum) - (middle2_sum * middle2_sum);
				float4 color_rdet2 = float4(1.0f / color_det2.x,
				                            1.0f / color_det2.y,
				                            1.0f / color_det2.z,
				                            1.0f / color_det2.w);

				float4 color_mss2 = (left2_sum * left2_sum)
				                  + (2.0f * middle2_sum * middle2_sum)
				                  + (right2_sum * right2_sum);

				float4 ep0 = (right2_sum * color_vec_x - middle2_sum * color_vec_y) * color_rdet2;
				float4 ep1 = (left2_sum * color_vec_y - middle2_sum * color_vec_x) * color_rdet2;

				if (plane2_color_component == 0 && fabsf(color_det2.x) > (color_mss2.x * 1e-4f) && ep0.x == ep0.x && ep1.x == ep1.x)
				{
					ep->endpt0[i].x = ep0.x;
					ep->endpt1[i].x = ep1.x;
				}

				if (plane2_color_component == 1 && fabsf(color_det2.y) > (color_mss2.y * 1e-4f) && ep0.y == ep0.y && ep1.y == ep1.y)
				{
					ep->endpt0[i].y = ep0.y;
					ep->endpt1[i].y = ep1.y;
				}

				if (plane2_color_component == 2 && fabsf(color_det2.z) > (color_mss2.z * 1e-4f) && ep0.z == ep0.z && ep1.z == ep1.z)
				{
					ep->endpt0[i].z = ep0.z;
					ep->endpt1[i].z = ep1.z;
				}

				if (plane2_color_component == 3 && fabsf(color_det2.w) > (color_mss2.w * 1e-4f) && ep0.w == ep0.w && ep1.w == ep1.w)
				{
					ep->endpt0[i].w = ep0.w;
					ep->endpt1[i].w = ep1.w;
				}

				#ifdef DEBUG_CAPTURE_NAN
					feenableexcept(FE_DIVBYZERO | FE_INVALID);
				#endif
			}
		}

		// if the calculation of an RGB-offset vector failed, try to compute
		// a somewhat-sensible value anyway
		if (rgbo_fail)
		{
			float4 v0 = ep->endpt0[i];
			float4 v1 = ep->endpt1[i];
			float avgdif = ((v1.x - v0.x) + (v1.y - v0.y) + (v1.z - v0.z)) * (1.0f / 3.0f);

			if (avgdif <= 0.0f)
			{
				avgdif = 0.0f;
			}

			float4 avg = (v0 + v1) * 0.5f;
			float4 ep0 = avg - float4(avgdif, avgdif, avgdif, avgdif) * 0.5f;

			rgbo_vectors[i] = float4(ep0.x, ep0.y, ep0.z, avgdif);
		}
	}
}

#endif
