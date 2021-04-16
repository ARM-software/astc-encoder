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
 * @brief Functions for computing color endpoints and texel weights.
 */

#include <cassert>

#include "astcenc_internal.h"
#include "astcenc_vecmathlib.h"

#ifdef DEBUG_CAPTURE_NAN
	#ifndef _GNU_SOURCE
		#define _GNU_SOURCE
	#endif

	#include <fenv.h>
#endif

static void compute_endpoints_and_ideal_weights_1_comp(
	const block_size_descriptor* bsd,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	endpoints_and_weights* ei,
	unsigned int component
) {
	int partition_count = pt->partition_count;
	ei->ep.partition_count = partition_count;
	promise(partition_count > 0);

	int texel_count = bsd->texel_count;
	promise(texel_count > 0);

	float lowvalues[4] { 1e10f, 1e10f, 1e10f, 1e10f };
	float highvalues[4] { -1e10f, -1e10f, -1e10f, -1e10f };

	float partition_error_scale[4];
	float linelengths_rcp[4];

	const float *error_weights = nullptr;
	const float* data_vr = nullptr;

	assert(component < 4);
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
	default:
		error_weights = ewb->texel_weight_a;
		data_vr = blk->data_a;
		break;
	}

	for (int i = 0; i < texel_count; i++)
	{
		if (error_weights[i] > 1e-10f)
		{
			float value = data_vr[i];
			int partition = pt->partition_of_texel[i];

			lowvalues[partition] = astc::min(value, lowvalues[partition]);
			highvalues[partition] = astc::max(value, highvalues[partition]);
		}
	}

	vmask4 sep_mask = vint4::lane_id() == vint4(component);
	for (int i = 0; i < partition_count; i++)
	{
		float diff = highvalues[i] - lowvalues[i];

		if (diff < 0)
		{
			lowvalues[i] = 0.0f;
			highvalues[i] = 0.0f;
		}

		diff = astc::max(diff, 1e-7f);

		partition_error_scale[i] = diff * diff;
		linelengths_rcp[i] = 1.0f / diff;

		ei->ep.endpt0[i] = select(blk->data_min, vfloat4(lowvalues[i]), sep_mask);
		ei->ep.endpt1[i] = select(blk->data_max, vfloat4(highvalues[i]), sep_mask);
	}

	bool is_constant_wes = true;
	float constant_wes = partition_error_scale[pt->partition_of_texel[0]] * error_weights[0];

	for (int i = 0; i < texel_count; i++)
	{
		float value = data_vr[i];
		int partition = pt->partition_of_texel[i];
		value -= lowvalues[partition];
		value *= linelengths_rcp[partition];
		value = astc::clamp1f(value);

		ei->weights[i] = value;
		ei->weight_error_scale[i] = partition_error_scale[partition] * error_weights[i];
		assert(!astc::isnan(ei->weight_error_scale[i]));

		is_constant_wes = is_constant_wes && ei->weight_error_scale[i] == constant_wes;
	}

	ei->is_constant_weight_error_scale = is_constant_wes;
}

static void compute_endpoints_and_ideal_weights_2_comp(
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
	promise(partition_count > 0);

	int texel_count = bsd->texel_count;
	promise(texel_count > 0);

	partition_metrics pms[4];

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

	compute_partition_error_color_weightings(*ewb, *pt, pms);

	for (int i = 0; i < partition_count; i++)
	{
		float s1 = 0, s2 = 0;
		assert(component1 < 4);
		switch (component1)
		{
		case 0:
			s1 = pms[i].color_scale.lane<0>();
			break;
		case 1:
			s1 = pms[i].color_scale.lane<1>();
			break;
		case 2:
			s1 = pms[i].color_scale.lane<2>();
			break;
		default:
			s1 = pms[i].color_scale.lane<3>();
			break;
		}

		assert(component2 < 4);
		switch (component2)
		{
		case 0:
			s2 = pms[i].color_scale.lane<0>();
			break;
		case 1:
			s2 = pms[i].color_scale.lane<1>();
			break;
		case 2:
			s2 = pms[i].color_scale.lane<2>();
			break;
		default:
			s2 = pms[i].color_scale.lane<3>();
			break;
		}
		scalefactors[i] = normalize(float2(s1, s2)) * 1.41421356f;
	}

	float lowparam[4] { 1e10f, 1e10f, 1e10f, 1e10f };
	float highparam[4] { -1e10f, -1e10f, -1e10f, -1e10f };

	float2 averages[4];
	float2 directions[4];

	line2 lines[4];
	float scale[4];
	float length_squared[4];

	compute_avgs_and_dirs_2_comp(pt, blk, ewb, scalefactors, component1, component2, averages, directions);

	for (int i = 0; i < partition_count; i++)
	{
		float2 dir = directions[i];
		if (dir.r + dir.g < 0.0f)
		{
			dir = float2(0.0f) - dir;
		}

		lines[i].a = averages[i];
		if (dot(dir, dir) == 0.0f)
		{
			lines[i].b = normalize(float2(1.0f));
		}
		else
		{
			lines[i].b = normalize(dir);
		}
	}

	for (int i = 0; i < texel_count; i++)
	{
		if (error_weights[i] > 1e-10f)
		{
			int partition = pt->partition_of_texel[i];
			float2 point = float2(data_vr[i], data_vg[i]) * scalefactors[partition];
			line2 l = lines[partition];
			float param = dot(point - l.a, l.b);
			ei->weights[i] = param;

			lowparam[partition] = astc::min(param, lowparam[partition]);
			highparam[partition] = astc::max(param, highparam[partition]);
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
		length = astc::max(length, 1e-7f);
		length_squared[i] = length * length;
		scale[i] = 1.0f / length;

		float2 ep0 = lines[i].a + lines[i].b * lowparam[i];
		float2 ep1 = lines[i].a + lines[i].b * highparam[i];

		ep0.r /= scalefactors[i].r;
		ep0.g /= scalefactors[i].g;

		ep1.r /= scalefactors[i].r;
		ep1.g /= scalefactors[i].g;

		lowvalues[i] = ep0;
		highvalues[i] = ep1;
	}

	vmask4 comp1_mask = vint4::lane_id() == vint4(component1);
	vmask4 comp2_mask = vint4::lane_id() == vint4(component2);
	for (int i = 0; i < partition_count; i++)
	{
		vfloat4 ep0 = select(blk->data_min, vfloat4(lowvalues[i].r), comp1_mask);
		vfloat4 ep1 = select(blk->data_max, vfloat4(highvalues[i].r), comp1_mask);

		ei->ep.endpt0[i] = select(ep0, vfloat4(lowvalues[i].g), comp2_mask);
		ei->ep.endpt1[i] = select(ep1, vfloat4(highvalues[i].g), comp2_mask);
	}

	bool is_constant_wes = true;
	float constant_wes = length_squared[pt->partition_of_texel[0]] * error_weights[0];

	for (int i = 0; i < texel_count; i++)
	{
		int partition = pt->partition_of_texel[i];
		float idx = (ei->weights[i] - lowparam[partition]) * scale[partition];
		idx = astc::clamp1f(idx);

		ei->weights[i] = idx;
		ei->weight_error_scale[i] = length_squared[partition] * error_weights[i];
		assert(!astc::isnan(ei->weight_error_scale[i]));

		is_constant_wes = is_constant_wes && ei->weight_error_scale[i] == constant_wes;
	}

	ei->is_constant_weight_error_scale = is_constant_wes;
}

static void compute_endpoints_and_ideal_weights_3_comp(
	const block_size_descriptor* bsd,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	endpoints_and_weights* ei,
	int omitted_component
) {
	int partition_count = pt->partition_count;
	ei->ep.partition_count = partition_count;
	promise(partition_count > 0);

	int texel_count= bsd->texel_count;
	promise(texel_count > 0);

	partition_metrics pms[4];

	const float *error_weights;
	const float* data_vr = nullptr;
	const float* data_vg = nullptr;
	const float* data_vb = nullptr;
	if (omitted_component == 0)
	{
		error_weights = ewb->texel_weight_gba;
		data_vr = blk->data_g;
		data_vg = blk->data_b;
		data_vb = blk->data_a;
	}
	else if (omitted_component == 1)
	{
		error_weights = ewb->texel_weight_rba;
		data_vr = blk->data_r;
		data_vg = blk->data_b;
		data_vb = blk->data_a;
	}
	else if (omitted_component == 2)
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

	compute_partition_error_color_weightings(*ewb, *pt, pms);

	for (int i = 0; i < partition_count; i++)
	{
		vfloat4 color_scale;
		assert(omitted_component < 4);
		switch (omitted_component)
		{
		case 0:
			color_scale = pms[i].color_scale.swz<1, 2, 3>();
			break;
		case 1:
			color_scale = pms[i].color_scale.swz<0, 2, 3>();
			break;
		case 2:
			color_scale = pms[i].color_scale.swz<0, 1, 3>();
			break;
		default:
			color_scale = pms[i].color_scale.swz<0, 1, 2>();
			break;
		}

		pms[i].color_scale = normalize(color_scale) * 1.73205080f;
	}

	float lowparam[4] { 1e10f, 1e10f, 1e10f, 1e10f };
	float highparam[4] { -1e10f, -1e10f, -1e10f, -1e10f };

	line3 lines[4];
	float scale[4];
	float length_squared[4];

	compute_avgs_and_dirs_3_comp(pt, blk, ewb, omitted_component, pms);

	for (int i = 0; i < partition_count; i++)
	{
		vfloat4 dir = pms[i].dir;
		if (hadd_rgb_s(dir) < 0.0f)
		{
			dir = vfloat4(0.0f) - dir;
		}

		lines[i].a = pms[i].avg;
		if (dot3_s(dir, dir) == 0.0f)
		{
			lines[i].b = normalize(vfloat4(1.0f, 1.0f, 1.0f, 0.0f));
		}
		else
		{
			lines[i].b = normalize(dir);
		}
	}

	for (int i = 0; i < texel_count; i++)
	{
		if (error_weights[i] > 1e-10f)
		{
			int partition = pt->partition_of_texel[i];
			vfloat4 point = vfloat4(data_vr[i], data_vg[i], data_vb[i], 0.0f) * pms[partition].color_scale;
			line3 l = lines[partition];
			float param = dot3_s(point - l.a, l.b);
			ei->weights[i] = param;

			lowparam[partition] = astc::min(param, lowparam[partition]);
			highparam[partition] = astc::max(param, highparam[partition]);
		}
		else
		{
			ei->weights[i] = -1e38f;
		}
	}

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
		length = astc::max(length, 1e-7f);

		length_squared[i] = length * length;
		scale[i] = 1.0f / length;

		vfloat4 ep0 = lines[i].a + lines[i].b * lowparam[i];
		vfloat4 ep1 = lines[i].a + lines[i].b * highparam[i];

		ep0 = ep0 / pms[i].color_scale;
		ep1 = ep1 / pms[i].color_scale;

		vfloat4 bmin = blk->data_min;
		vfloat4 bmax = blk->data_max;

		// TODO: Probably a programmatic vector permute we can do here ...
		assert(omitted_component < 4);
		switch (omitted_component)
		{
			case 0:
				ei->ep.endpt0[i] = vfloat4(bmin.lane<0>(), ep0.lane<0>(), ep0.lane<1>(), ep0.lane<2>());
				ei->ep.endpt1[i] = vfloat4(bmax.lane<0>(), ep1.lane<0>(), ep1.lane<1>(), ep1.lane<2>());
				break;
			case 1:
				ei->ep.endpt0[i] = vfloat4(ep0.lane<0>(), bmin.lane<1>(), ep0.lane<1>(), ep0.lane<2>());
				ei->ep.endpt1[i] = vfloat4(ep1.lane<0>(), bmax.lane<1>(), ep1.lane<1>(), ep1.lane<2>());
				break;
			case 2:
				ei->ep.endpt0[i] = vfloat4(ep0.lane<0>(), ep0.lane<1>(), bmin.lane<2>(), ep0.lane<2>());
				ei->ep.endpt1[i] = vfloat4(ep1.lane<0>(), ep1.lane<1>(), bmax.lane<2>(), ep1.lane<2>());
				break;
			default:
				ei->ep.endpt0[i] = vfloat4(ep0.lane<0>(), ep0.lane<1>(), ep0.lane<2>(), bmin.lane<3>());
				ei->ep.endpt1[i] = vfloat4(ep1.lane<0>(), ep1.lane<1>(), ep1.lane<2>(), bmax.lane<3>());
				break;
		}
	}


	bool is_constant_wes = true;
	float constant_wes = length_squared[pt->partition_of_texel[0]] * error_weights[0];

	for (int i = 0; i < texel_count; i++)
	{
		int partition = pt->partition_of_texel[i];
		float idx = (ei->weights[i] - lowparam[partition]) * scale[partition];
		idx = astc::clamp1f(idx);

		ei->weights[i] = idx;
		ei->weight_error_scale[i] = length_squared[partition] * error_weights[i];
		assert(!astc::isnan(ei->weight_error_scale[i]));

		is_constant_wes = is_constant_wes && ei->weight_error_scale[i] == constant_wes;
	}

	ei->is_constant_weight_error_scale = is_constant_wes;
}

static void compute_endpoints_and_ideal_weights_4_comp(
	const block_size_descriptor* bsd,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	endpoints_and_weights* ei
) {
	const float *error_weights = ewb->texel_weight;

	int partition_count = pt->partition_count;

	int texel_count= bsd->texel_count;
	promise(texel_count > 0);
	promise(partition_count > 0);

	float lowparam[4] { 1e10, 1e10, 1e10, 1e10 };
	float highparam[4] {  -1e10,  -1e10,  -1e10, -1e10 };

	line4 lines[4];

	float scale[4];
	float length_squared[4];

	partition_metrics pms[4];

	compute_partition_error_color_weightings(*ewb, *pt, pms);

	for (int i = 0; i < partition_count; i++)
	{
		pms[i].color_scale = normalize(pms[i].color_scale) * 2.0f;
	}

	compute_avgs_and_dirs_4_comp(pt, blk, ewb, pms);

	// if the direction-vector ends up pointing from light to dark, FLIP IT!
	// this will make the first endpoint the darkest one.
	for (int i = 0; i < partition_count; i++)
	{
		vfloat4 dir = pms[i].dir;
		if (hadd_rgb_s(dir) < 0.0f)
		{
			dir = vfloat4::zero() - dir;
		}

		lines[i].a = pms[i].avg;
		if (dot_s(dir, dir) == 0.0f)
		{
			lines[i].b = normalize(vfloat4(1.0f));
		}
		else
		{
			lines[i].b = normalize(dir);
		}
	}

	for (int i = 0; i < texel_count; i++)
	{
		if (error_weights[i] > 1e-10f)
		{
			int partition = pt->partition_of_texel[i];

			vfloat4 point = blk->texel(i) * pms[partition].color_scale;
			line4 l = lines[partition];

			float param = dot_s(point - l.a, l.b);
			ei->weights[i] = param;

			lowparam[partition] = astc::min(param, lowparam[partition]);
			highparam[partition] = astc::max(param, highparam[partition]);
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
		length = astc::max(length, 1e-7f);

		length_squared[i] = length * length;
		scale[i] = 1.0f / length;

		vfloat4 ep0 = lines[i].a + lines[i].b * lowparam[i];
		vfloat4 ep1 = lines[i].a + lines[i].b * highparam[i];

		ei->ep.endpt0[i] = ep0 / pms[i].color_scale;
		ei->ep.endpt1[i] = ep1 / pms[i].color_scale;
	}

	bool is_constant_wes = true;
	float constant_wes = length_squared[pt->partition_of_texel[0]] * error_weights[0];

	for (int i = 0; i < texel_count; i++)
	{
		int partition = pt->partition_of_texel[i];
		float idx = (ei->weights[i] - lowparam[partition]) * scale[partition];
		idx = astc::clamp1f(idx);

		ei->weights[i] = idx;
		ei->weight_error_scale[i] = error_weights[i] * length_squared[partition];
		assert(!astc::isnan(ei->weight_error_scale[i]));

		is_constant_wes = is_constant_wes && ei->weight_error_scale[i] == constant_wes;
	}

	ei->is_constant_weight_error_scale = is_constant_wes;
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
		compute_endpoints_and_ideal_weights_4_comp(bsd, pt, blk, ewb, ei);
	}
	else
	{
		compute_endpoints_and_ideal_weights_3_comp(bsd, pt, blk, ewb, ei, 3);
	}
}

void compute_endpoints_and_ideal_weights_2_planes(
	const block_size_descriptor* bsd,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	int plane2_component,
	endpoints_and_weights* ei1,
	endpoints_and_weights* ei2
) {
	int uses_alpha = imageblock_uses_alpha(blk);

	assert(plane2_component < 4);
	switch (plane2_component)
	{
	case 0: // separate weights for red
		if (uses_alpha)
		{
			compute_endpoints_and_ideal_weights_3_comp(bsd, pt, blk, ewb, ei1, 0);
		}
		else
		{
			compute_endpoints_and_ideal_weights_2_comp(bsd, pt, blk, ewb, ei1, 1, 2);
		}
		compute_endpoints_and_ideal_weights_1_comp(bsd, pt, blk, ewb, ei2, 0);
		break;

	case 1: // separate weights for green
		if (uses_alpha)
		{
			compute_endpoints_and_ideal_weights_3_comp(bsd, pt, blk, ewb, ei1, 1);
		}
		else
		{
			compute_endpoints_and_ideal_weights_2_comp(bsd, pt, blk, ewb, ei1, 0, 2);
		}
		compute_endpoints_and_ideal_weights_1_comp(bsd, pt, blk, ewb, ei2, 1);
		break;

	case 2: // separate weights for blue
		if (uses_alpha)
		{
			compute_endpoints_and_ideal_weights_3_comp(bsd, pt, blk, ewb, ei1, 2);
		}
		else
		{
			compute_endpoints_and_ideal_weights_2_comp(bsd, pt, blk, ewb, ei1, 0, 1);
		}
		compute_endpoints_and_ideal_weights_1_comp(bsd, pt, blk, ewb, ei2, 2);
		break;

	default: // separate weights for alpha
		assert(uses_alpha);
		compute_endpoints_and_ideal_weights_3_comp(bsd, pt, blk, ewb, ei1, 3);
		compute_endpoints_and_ideal_weights_1_comp(bsd, pt, blk, ewb, ei2, 3);
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

float compute_error_of_weight_set(
	const endpoints_and_weights* eai,
	const decimation_table* dt,
	const float* weights
) {
	vfloat4 error_summav = vfloat4::zero();
	float error_summa = 0.0f;
	int texel_count = dt->texel_count;

	int i = 0;

	// Process SIMD-width texel coordinates at at time while we can
	int clipped_texel_count = round_down_to_simd_multiple_vla(texel_count);
	for (/* */; i < clipped_texel_count; i += ASTCENC_SIMD_WIDTH)
	{
		// Compute the bilinear interpolation of the decimated weight grid
		vfloat current_values = bilinear_infill_vla(*dt, weights, i);

		// Compute the error between the computed value and the ideal weight
		vfloat actual_values = loada(&(eai->weights[i]));
		vfloat diff = current_values - actual_values;
		vfloat significance = loada(&(eai->weight_error_scale[i]));
		vfloat error = diff * diff * significance;

		haccumulate(error_summav, error);
	}

	// Loop tail
	// Error is buffered and accumulated in blocks of 4 to ensure that
	// the partial sums added to the accumulator are invariant with the
	// vector implementation, irrespective of vector size ...
	alignas(16) float errorsum_tmp[4] { 0 };
	for (/* */; i < texel_count; i++)
	{
		float current_value = bilinear_infill(*dt, weights, i);
		float valuedif = current_value - eai->weights[i];
		float error = valuedif * valuedif * eai->weight_error_scale[i];

		// Accumulate error sum in the temporary array
		int error_index = i & 0x3;
		errorsum_tmp[error_index] = error;

#if ASTCENC_SIMD_WIDTH == 8
		// Zero the temporary staging buffer every 4 items unless last. Note
		// that this block can only trigger for 6x5 blocks, all other partials
		// tails are shorter than 4 ...
		if ((i & 0x7) == 0x03)
		{
			haccumulate(error_summav, vfloat4::loada(errorsum_tmp));
 			storea(vfloat4::zero(), errorsum_tmp);
		}
#endif
	}

	// Accumulate the loop tail using the vfloat4 swizzle
	haccumulate(error_summav, vfloat4::loada(errorsum_tmp));

	// Resolve the final scalar accumulator sum
	haccumulate(error_summa, error_summav);

	return error_summa;
}

/* See header for documentation. */
// Note: This function is vectorized, but needs to use gathers to access the
// decimation table structures so vectorization is currently only enabled for
// AVX2. The implementation loops over decimated weights, and then texels for
// each weight. We know the backing memory is "large enough" we can can
// overshoot the weight count to always use full vectors without a loop tail.
// The inner loop operates on 8 weights, each of which may have a different
// number of texels referenced by it. We iterate over the max reference count,
// and then use lane masks to disable lanes that are no longer in scope.
void compute_ideal_weights_for_decimation_table(
	const endpoints_and_weights& eai_in,
	endpoints_and_weights& eai_out,
	const decimation_table& dt,
	float* RESTRICT weight_set,
	float* RESTRICT weights
) {
	int i;
	int texel_count = dt.texel_count;
	int weight_count = dt.weight_count;

	promise(texel_count > 0);
	promise(weight_count > 0);

	// This function includes a copy of the epw from eai_in to eai_out. We do it
	// here because we want to load the data anyway, so we can avoid loading it
	// from memory twice.
	eai_out.ep = eai_in.ep;
	eai_out.is_constant_weight_error_scale = eai_in.is_constant_weight_error_scale;

	// If we have a 1:1 mapping just shortcut the computation - clone the
	// weights into both the weight set and the output epw copy.
	if (texel_count == weight_count)
	{
		for (i = 0; i < texel_count; i++)
		{
			assert(i == dt.weight_texel[0][i]);
			weight_set[i] = eai_in.weights[i];
			weights[i] = eai_in.weight_error_scale[i];

			eai_out.weights[i] = eai_in.weights[i];
			eai_out.weight_error_scale[i] = eai_in.weight_error_scale[i];
		}
		return;
	}
	// If we don't have a 1:1 mapping just clone the weights into the output
	// epw copy and then do the full algorithm to decimate weights.
	else
	{
		for (i = 0; i < texel_count; i++)
		{
			eai_out.weights[i] = eai_in.weights[i];
			eai_out.weight_error_scale[i] = eai_in.weight_error_scale[i];
		}
	}

	// Otherwise compute an estimate and perform single refinement iteration
	alignas(ASTCENC_VECALIGN) float infilled_weights[MAX_TEXELS_PER_BLOCK];

	// Compute an initial average for each decimated weight
	i = 0;

#if ASTCENC_SIMD_WIDTH >= 8
	int clipped_weight_count = round_down_to_simd_multiple_vla(weight_count);

	bool constant_wes = eai_in.is_constant_weight_error_scale;
	vfloat weight_error_scale(eai_in.weight_error_scale[0]);

	for (/* */; i < clipped_weight_count; i += ASTCENC_SIMD_WIDTH)
	{
		// Start with a small value to avoid div-by-zero later
		vfloat weight_weight(1e-10f);
		vfloat initial_weight = vfloat::zero();

		// Accumulate error weighting of all the texels using this weight
		vint weight_texel_count(dt.weight_texel_count + i);
		int max_texel_count = hmax(weight_texel_count).lane<0>();
		promise(max_texel_count > 0);

		for (int j = 0; j < max_texel_count; j++)
		{
			// Not all lanes may actually use j texels, so mask out if idle
			vmask active = weight_texel_count > vint(j);

			vint texel(dt.weight_texel[j] + i);
			texel = select(vint::zero(), texel, active);

			vfloat weight = loada(dt.weights_flt[j] + i);
			weight = select(vfloat::zero(), weight, active);

			if (!constant_wes)
			{
				weight_error_scale = gatherf(eai_in.weight_error_scale, texel);
			}

			vfloat contrib_weight = weight * weight_error_scale;

			weight_weight = weight_weight + contrib_weight;
			initial_weight = initial_weight + gatherf(eai_in.weights, texel) * contrib_weight;
		}

		storea(weight_weight, weights + i);
		storea(initial_weight / weight_weight, weight_set + i);
	}
#endif

	// Loop tail
	for (/* */; i < weight_count; i++)
	{
		// Start with a small value to avoid div-by-zero later
		float weight_weight = 1e-10f;
		float initial_weight = 0.0f;

		// Accumulate error weighting of all the texels using this weight
		int weight_texel_count = dt.weight_texel_count[i];
		promise(weight_texel_count > 0);

		for (int j = 0; j < weight_texel_count; j++)
		{
			int texel = dt.weight_texel[j][i];
			float weight = dt.weights_flt[j][i];
			float contrib_weight = weight * eai_in.weight_error_scale[texel];
			weight_weight += contrib_weight;
			initial_weight += eai_in.weights[texel] * contrib_weight;
		}

		weights[i] = weight_weight;
		weight_set[i] = initial_weight / weight_weight;
	}

	// Populate the interpolated weight grid based on the initital average
	i = 0;

#if ASTCENC_SIMD_WIDTH >= 8
	// Process SIMD-width texel coordinates at at time while we can
	int clipped_texel_count = round_down_to_simd_multiple_vla(texel_count);
	for (/* */; i < clipped_texel_count; i += ASTCENC_SIMD_WIDTH)
	{
		vfloat weight = bilinear_infill_vla(dt, weight_set, i);
		storea(weight, infilled_weights + i);
	}
#endif

	// Loop tail
	for (/* */; i < texel_count; i++)
	{
		infilled_weights[i] = bilinear_infill(dt, weight_set, i);
	}

	// Perform a single iteration of refinement
	constexpr float stepsize = 0.25f;
	constexpr float chd_scale = -TEXEL_WEIGHT_SUM;
	i = 0;

#if ASTCENC_SIMD_WIDTH >= 8
	for (/* */; i < clipped_weight_count; i += ASTCENC_SIMD_WIDTH)
	{
		// Start with a small value to avoid div-by-zero later
		vfloat weight_val = loada(weight_set + i);

		// Accumulate error weighting of all the texels using this weight
		vfloat error_change0(1e-10f);
		vfloat error_change1(0.0f);

		// Accumulate error weighting of all the texels using this weight
		vint weight_texel_count(dt.weight_texel_count + i);
		int max_texel_count = hmax(weight_texel_count).lane<0>();
		promise(max_texel_count > 0);

		for (int j = 0; j < max_texel_count; j++)
		{
			// Not all lanes may actually use j texels, so mask out if idle
			vmask active = weight_texel_count > vint(j);

			vint texel(dt.weight_texel[j] + i);
			texel = select(vint::zero(), texel, active);

			vfloat contrib_weight = loada(dt.weights_flt[j] + i);
			contrib_weight = select(vfloat::zero(), contrib_weight, active);

			if (!constant_wes)
			{
 				weight_error_scale = gatherf(eai_in.weight_error_scale, texel);
			}

			vfloat scale = weight_error_scale * contrib_weight;
			vfloat old_weight = gatherf(infilled_weights, texel);
			vfloat ideal_weight = gatherf(eai_in.weights, texel);

			error_change0 = error_change0 + contrib_weight * scale;
			error_change1 = error_change1 + (old_weight - ideal_weight) * scale;
		}

		vfloat step = (error_change1 * chd_scale) / error_change0;
		step = clamp(-stepsize, stepsize, step);

		// update the weight
		storea(weight_val + step, weight_set + i);
	}
#endif

	// Loop tail
	for (/* */; i < weight_count; i++)
	{
		float weight_val = weight_set[i];

		// Start with a small value to avoid div-by-zero later
		float error_change0 = 1e-10f;
		float error_change1 = 0.0f;

		// Compute the two error changes that occur from perturbing the current index
		int weight_texel_count = dt.weight_texel_count[i];
		promise(weight_texel_count > 0);
		for (int k = 0; k < weight_texel_count; k++)
		{
			uint8_t texel = dt.weight_texel[k][i];
			float contrib_weight = dt.weights_flt[k][i];

			float scale = eai_in.weight_error_scale[texel] * contrib_weight;
			float old_weight = infilled_weights[texel];
			float ideal_weight = eai_in.weights[texel];

			error_change0 +=  contrib_weight * scale;
			error_change1 += (old_weight - ideal_weight) * scale;
		}

		float step = (error_change1 * chd_scale) / error_change0;
		step = astc::clamp(step, -stepsize, stepsize);

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
void compute_quantized_weights_for_decimation_table(
	const decimation_table* dt,
	float low_bound,
	float high_bound,
	const float* weight_set_in,
	float* weight_set_out,
	uint8_t* quantized_weight_set,
	int quant_level
) {
	int weight_count = dt->weight_count;
	const quantization_and_transfer_table *qat = &(quant_and_xfer_tables[quant_level]);

	static const int quant_levels[12] { 2,3,4,5,6,8,10,12,16,20,24,32 };
	float quant_level_m1 = (float)(quant_levels[quant_level] - 1);

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

#if ASTCENC_SIMD_WIDTH > 1
	// SIMD loop; process weights in SIMD width batches while we can.
	vfloat scalev(scale);
	vfloat scaled_low_boundv(scaled_low_bound);
	vfloat quant_level_m1v(quant_level_m1);
	vfloat rscalev(rscale);
	vfloat low_boundv(low_bound);

	int clipped_weight_count = round_down_to_simd_multiple_vla(weight_count);
	for (/* */; i < clipped_weight_count; i += ASTCENC_SIMD_WIDTH)
	{
		vfloat ix = loada(&weight_set_in[i]) * scalev - scaled_low_boundv;
		ix = clampzo(ix);

		//Llook up the two closest indexes and return the one that was closest.
		vfloat ix1 = ix * quant_level_m1v;
		vint weight = float_to_int(ix1);
		vint weight1 = weight + vint(1);
		vfloat ixl = gatherf(qat->unquantized_value_unsc, weight);
		vfloat ixh = gatherf(qat->unquantized_value_unsc, weight1);

		vmask mask = (ixl + ixh) < (vfloat(128.0f) * ix);
		weight = select(weight, weight1, mask);
		ixl = select(ixl, ixh, mask);

		// Invert the weight-scaling that was done initially
		storea(ixl * rscalev + low_boundv, &weight_set_out[i]);
		vint scm = gatheri(qat->scramble_map, weight);
		vint scn = pack_low_bytes(scm);
		store_nbytes(scn, &quantized_weight_set[i]);
	}
#endif // #if ASTCENC_SIMD_WIDTH > 1

	// Loop tail
	for (/* */; i < weight_count; i++)
	{
		float ix = (weight_set_in[i] * scale) - scaled_low_bound;
		ix = astc::clamp1f(ix);

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

static inline vfloat4 compute_rgbovec(
	vfloat4 rgba_weight_sum,
	vfloat4 weight_weight_sum,
	vfloat4 rgbq_sum,
	float psum
) {
	// Compute the rgb+offset for HDR endpoint mode #7. Since the matrix needed
	// has a regular structure, we can simplify the inverse calculation. This
	// gives us ~24 multiplications, down from 96 for a generic inverse

	// mat[0] = vfloat4(rgba_ws.x,      0.0f,      0.0f, wght_ws.x);
	// mat[1] = vfloat4(     0.0f, rgba_ws.y,      0.0f, wght_ws.y);
	// mat[2] = vfloat4(     0.0f,      0.0f, rgba_ws.z, wght_ws.z);
	// mat[3] = vfloat4(wght_ws.x, wght_ws.y, wght_ws.z,      psum);
	// mat = invert(mat);

	float X = rgba_weight_sum.lane<0>();
	float Y = rgba_weight_sum.lane<1>();
	float Z = rgba_weight_sum.lane<2>();
	float P = weight_weight_sum.lane<0>();
	float Q = weight_weight_sum.lane<1>();
	float R = weight_weight_sum.lane<2>();
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
	vfloat4 mat0(DT, ZQP, RYP, mZYP);
	vfloat4 mat1(ZQP, SZmRR * X - Z * PP, RQX, mZQX);
	vfloat4 mat2(RYP, RQX, (S * Y - QQ) * X - Y * PP, mRYX);
	vfloat4 mat3(mZYP, mZQX, mRYX, Z * YX);
	vfloat4 vect = rgbq_sum * rdet;

	#ifdef DEBUG_CAPTURE_NAN
	    fedisableexcept(FE_DIVBYZERO | FE_INVALID);
	#endif

	return vfloat4(dot_s(mat0, vect),
	               dot_s(mat1, vect),
	               dot_s(mat2, vect),
	               dot_s(mat3, vect));
}

/* for a given weight set, we wish to recompute the colors so that they are optimal for a particular weight set. */
void recompute_ideal_colors_2planes(
	int weight_quant_mode,
	endpoints* ep,	// contains the endpoints we wish to update
	vfloat4* rgbs_vectors,	// used to return RGBS-vectors for endpoint mode #6
	vfloat4* rgbo_vectors,	// used to return RGBO-vectors for endpoint mode #7
	const uint8_t* weight_set8,	// the current set of weight values
	const uint8_t* plane2_weight_set8,	// nullptr if plane 2 is not actually used.
	int plane2_component,	// color component for 2nd plane of weights; -1 if the 2nd plane of weights is not present
	const partition_info* pt,
	const decimation_table* dt,
	const imageblock* blk,	// picture-block containing the actual data.
	const error_weight_block* ewb
) {
	const quantization_and_transfer_table *qat = &(quant_and_xfer_tables[weight_quant_mode]);

	float weight_set[MAX_WEIGHTS_PER_BLOCK];
	float plane2_weight_set[MAX_WEIGHTS_PER_BLOCK];

	for (int i = 0; i < dt->weight_count; i++)
	{
		weight_set[i] = qat->unquantized_value[weight_set8[i]] * (1.0f / 64.0f);
	}

	if (plane2_weight_set8)
	{
		for (int i = 0; i < dt->weight_count; i++)
		{
			plane2_weight_set[i] = qat->unquantized_value[plane2_weight_set8[i]] * (1.0f / 64.0f);
		}
	}

	int partition_count = pt->partition_count;

	for (int i = 0; i < partition_count; i++)
	{
		vfloat4 rgba_sum(1e-17f);
		vfloat4 rgba_weight_sum(1e-17f);

		int texelcount = pt->partition_texel_count[i];
		const uint8_t *texel_indexes = pt->texels_of_partition[i];
		for (int j = 0; j < texelcount; j++)
		{
			int tix = texel_indexes[j];

			vfloat4 rgba = blk->texel(tix);
			vfloat4 error_weight(ewb->texel_weight_r[tix], ewb->texel_weight_g[tix], ewb->texel_weight_b[tix], ewb->texel_weight_a[tix]);

			rgba_sum += rgba * error_weight;
			rgba_weight_sum += error_weight;
		}

		vfloat4 scale_direction = normalize((rgba_sum * (1.0f / rgba_weight_sum)).swz<0, 1, 2>());

		float scale_max = 0.0f;
		float scale_min = 1e10f;

		float wmin1 = 1.0f;
		float wmax1 = 0.0f;
		float wmin2 = 1.0f;
		float wmax2 = 0.0f;

		vfloat4 left_sum    = vfloat4::zero();
		vfloat4 middle_sum  = vfloat4::zero();
		vfloat4 right_sum   = vfloat4::zero();

		vfloat4 left2_sum   = vfloat4::zero();
		vfloat4 middle2_sum = vfloat4::zero();
		vfloat4 right2_sum  = vfloat4::zero();
		vfloat4 lmrs_sum    = vfloat4::zero();

		vfloat4 color_vec_x = vfloat4::zero();
		vfloat4 color_vec_y = vfloat4::zero();

		float2 scale_vec = float2(0.0f);

		vfloat4 weight_weight_sum = vfloat4(1e-17f);
		float psum = 1e-17f;

		// FIXME: the loop below has too many responsibilities, making it inefficient.
		for (int j = 0; j < texelcount; j++)
		{
			int tix = texel_indexes[j];

			vfloat4 rgba = blk->texel(tix);
			vfloat4 color_weight(ewb->texel_weight_r[tix], ewb->texel_weight_g[tix], ewb->texel_weight_b[tix], ewb->texel_weight_a[tix]);

			// FIXME: move this calculation out to the color block.
			float ls_weight = hadd_rgb_s(color_weight);

			float idx0 = bilinear_infill(*dt, weight_set, tix);

			float om_idx0 = 1.0f - idx0;
			wmin1 = astc::min(idx0, wmin1);
			wmax1 = astc::max(idx0, wmax1);

			float scale = dot3_s(scale_direction, rgba);
			scale_min = astc::min(scale, scale_min);
			scale_max = astc::max(scale, scale_max);

			vfloat4 left   = color_weight * (om_idx0 * om_idx0);
			vfloat4 middle = color_weight * (om_idx0 * idx0);
			vfloat4 right  = color_weight * (idx0 * idx0);

			vfloat4 lmrs = vfloat4(om_idx0 * om_idx0,
			                       om_idx0 * idx0,
			                       idx0 * idx0,
			                       0.0f) * ls_weight;

			left_sum   += left;
			middle_sum += middle;
			right_sum  += right;
			lmrs_sum   += lmrs;

			float idx1 = 0.0f;
			float om_idx1 = 0.0f;

			idx1 = bilinear_infill(*dt, plane2_weight_set, tix);

			om_idx1 = 1.0f - idx1;
			wmin2 = astc::min(idx1, wmin2);
			wmax2 = astc::max(idx1, wmax2);

			vfloat4 left2   = color_weight * (om_idx1 * om_idx1);
			vfloat4 middle2 = color_weight * (om_idx1 * idx1);
			vfloat4 right2  = color_weight * (idx1 * idx1);

			left2_sum   += left2;
			middle2_sum += middle2;
			right2_sum  += right2;

			vmask4 p2_mask = vint4::lane_id() == vint4(plane2_component);
			vfloat4 color_idx = select(vfloat4(idx0), vfloat4(idx1), p2_mask);

			vfloat4 cwprod = color_weight * rgba;
			vfloat4 cwiprod = cwprod * color_idx;

			color_vec_y += cwiprod;
			color_vec_x += cwprod - cwiprod;

			scale_vec += float2(om_idx0, idx0) * (ls_weight * scale);
			weight_weight_sum += (color_weight * color_idx);
			psum += dot3_s(color_weight * color_idx, color_idx);
		}

		// calculations specific to mode #7, the HDR RGB-scale mode.
		// FIXME: Can we skip this for LDR textures?
		vfloat4 rgbq_sum = color_vec_x + color_vec_y;
		rgbq_sum.set_lane<3>(hadd_rgb_s(color_vec_y));

		#ifdef DEBUG_CAPTURE_NAN
		    fedisableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif

		vfloat4 rgbovec = compute_rgbovec(rgba_weight_sum, weight_weight_sum,
		                                  rgbq_sum, psum);
		rgbo_vectors[i] = rgbovec;

		// We will occasionally get a failure due to the use of a singular
		// (non-invertible) matrix. Record whether such a failure has taken
		// place; if it did, compute rgbo_vectors[] with a different method
		// later on.
		float chkval = dot_s(rgbovec, rgbovec);
		int rgbo_fail = chkval != chkval;

		// Initialize the luminance and scale vectors with a reasonable
		//  default, just in case the subsequent calculation blows up.
		#ifdef DEBUG_CAPTURE_NAN
			fedisableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif

		float scalediv = scale_min * (1.0f / astc::max(scale_max, 1e-10f));
		scalediv = astc::clamp1f(scalediv);

		#ifdef DEBUG_CAPTURE_NAN
			feenableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif

		vfloat4 sds = scale_direction * scale_max;

		rgbs_vectors[i] = vfloat4(sds.lane<0>(), sds.lane<1>(), sds.lane<2>(), scalediv);

		if (wmin1 >= wmax1 * 0.999f)
		{
			// if all weights in the partition were equal, then just take average
			// of all colors in the partition and use that as both endpoint colors.
			vfloat4 avg = (color_vec_x + color_vec_y) * (1.0f / rgba_weight_sum);

			vmask4 p1_mask = vint4::lane_id() != vint4(plane2_component);
			vmask4 notnan_mask = avg == avg;
			vmask4 full_mask = p1_mask & notnan_mask;

			ep->endpt0[i] = select(ep->endpt0[i], avg, full_mask);
			ep->endpt1[i] = select(ep->endpt1[i], avg, full_mask);

			rgbs_vectors[i] = vfloat4(sds.lane<0>(), sds.lane<1>(), sds.lane<2>(), 1.0f);
		}
		else
		{
			// otherwise, complete the analytic calculation of ideal-endpoint-values
			// for the given set of texel weights and pixel colors.

			#ifdef DEBUG_CAPTURE_NAN
			    fedisableexcept(FE_DIVBYZERO | FE_INVALID);
			#endif

			vfloat4 color_det1 = (left_sum * right_sum) - (middle_sum * middle_sum);
			vfloat4 color_rdet1 = 1.0f / color_det1;

			float ls_det1  = (lmrs_sum.lane<0>() * lmrs_sum.lane<2>()) - (lmrs_sum.lane<1>() * lmrs_sum.lane<1>());
			float ls_rdet1 = 1.0f / ls_det1;

			vfloat4 color_mss1 = (left_sum * left_sum)
			                   + (2.0f * middle_sum * middle_sum)
			                   + (right_sum * right_sum);

			float ls_mss1 = (lmrs_sum.lane<0>() * lmrs_sum.lane<0>())
			              + (2.0f * lmrs_sum.lane<1>() * lmrs_sum.lane<1>())
			              + (lmrs_sum.lane<2>() * lmrs_sum.lane<2>());

			vfloat4 ep0 = (right_sum * color_vec_x - middle_sum * color_vec_y) * color_rdet1;
			vfloat4 ep1 = (left_sum * color_vec_y - middle_sum * color_vec_x) * color_rdet1;

			float scale_ep0 = (lmrs_sum.lane<2>() * scale_vec.r - lmrs_sum.lane<1>() * scale_vec.g) * ls_rdet1;
			float scale_ep1 = (lmrs_sum.lane<0>() * scale_vec.g - lmrs_sum.lane<1>() * scale_vec.r) * ls_rdet1;

			vmask4 p1_mask = vint4::lane_id() != vint4(plane2_component);
			vmask4 det_mask = abs(color_det1) > (color_mss1 * 1e-4f);
			vmask4 notnan_mask = (ep0 == ep0) & (ep1 == ep1);
			vmask4 full_mask = p1_mask & det_mask & notnan_mask;

			ep->endpt0[i] = select(ep->endpt0[i], ep0, full_mask);
			ep->endpt1[i] = select(ep->endpt1[i], ep1, full_mask);

			if (fabsf(ls_det1) > (ls_mss1 * 1e-4f) && scale_ep0 == scale_ep0 && scale_ep1 == scale_ep1 && scale_ep0 < scale_ep1)
			{
				float scalediv2 = scale_ep0 * (1.0f / scale_ep1);
				vfloat4 sdsm = scale_direction * scale_ep1;
				rgbs_vectors[i] = vfloat4(sdsm.lane<0>(), sdsm.lane<1>(), sdsm.lane<2>(), scalediv2);
			}

			#ifdef DEBUG_CAPTURE_NAN
				feenableexcept(FE_DIVBYZERO | FE_INVALID);
			#endif
		}

		if (wmin2 >= wmax2 * 0.999f)
		{
			// if all weights in the partition were equal, then just take average
			// of all colors in the partition and use that as both endpoint colors.
			vfloat4 avg = (color_vec_x + color_vec_y) * (1.0f / rgba_weight_sum);

			vmask4 p2_mask = vint4::lane_id() == vint4(plane2_component);
			vmask4 notnan_mask = avg == avg;
			vmask4 full_mask = p2_mask & notnan_mask;

			ep->endpt0[i] = select(ep->endpt0[i], avg, full_mask);
			ep->endpt1[i] = select(ep->endpt1[i], avg, full_mask);
		}
		else
		{
			#ifdef DEBUG_CAPTURE_NAN
				fedisableexcept(FE_DIVBYZERO | FE_INVALID);
			#endif

			// otherwise, complete the analytic calculation of ideal-endpoint-values
			// for the given set of texel weights and pixel colors.
			vfloat4 color_det2 = (left2_sum * right2_sum) - (middle2_sum * middle2_sum);
			vfloat4 color_rdet2 = 1.0f / color_det2;

			vfloat4 color_mss2 = (left2_sum * left2_sum)
			                   + (2.0f * middle2_sum * middle2_sum)
			                   + (right2_sum * right2_sum);

			vfloat4 ep0 = (right2_sum * color_vec_x - middle2_sum * color_vec_y) * color_rdet2;
			vfloat4 ep1 = (left2_sum * color_vec_y - middle2_sum * color_vec_x) * color_rdet2;

			vmask4 p2_mask = vint4::lane_id() == vint4(plane2_component);
			vmask4 det_mask = abs(color_det2) > (color_mss2 * 1e-4f);
			vmask4 notnan_mask = (ep0 == ep0) & (ep1 == ep1);
			vmask4 full_mask = p2_mask & det_mask & notnan_mask;

			ep->endpt0[i] = select(ep->endpt0[i], ep0, full_mask);
			ep->endpt1[i] = select(ep->endpt1[i], ep1, full_mask);

			#ifdef DEBUG_CAPTURE_NAN
				feenableexcept(FE_DIVBYZERO | FE_INVALID);
			#endif
		}

		// if the calculation of an RGB-offset vector failed, try to compute
		// a somewhat-sensible value anyway
		if (rgbo_fail)
		{
			vfloat4 v0 = ep->endpt0[i];
			vfloat4 v1 = ep->endpt1[i];

			float avgdif = hadd_rgb_s(v1 - v0) * (1.0f / 3.0f);
			avgdif = astc::max(avgdif, 0.0f);

			vfloat4 avg = (v0 + v1) * 0.5f;
			vfloat4 ep0 = avg - vfloat4(avgdif) * 0.5f;

			rgbo_vectors[i] = vfloat4(ep0.lane<0>(), ep0.lane<1>(), ep0.lane<2>(), avgdif);
		}
	}
}

/* for a given weight set, we wish to recompute the colors so that they are optimal for a particular weight set. */
void recompute_ideal_colors_1plane(
	int weight_quant_mode,
	endpoints* ep,	// contains the endpoints we wish to update
	vfloat4* rgbs_vectors,	// used to return RGBS-vectors for endpoint mode #6
	vfloat4* rgbo_vectors,	// used to return RGBO-vectors for endpoint mode #7
	const uint8_t* weight_set8,	// the current set of weight values
	const partition_info* pt,
	const decimation_table* dt,
	const imageblock* blk,	// picture-block containing the actual data.
	const error_weight_block* ewb
) {
	int weight_count = dt->weight_count;
	int partition_count = pt->partition_count;

	promise(weight_count > 0);
	promise(partition_count > 0);

	const quantization_and_transfer_table *qat = &(quant_and_xfer_tables[weight_quant_mode]);

	float weight_set[MAX_WEIGHTS_PER_BLOCK];
	for (int i = 0; i < weight_count; i++)
	{
		weight_set[i] = qat->unquantized_value[weight_set8[i]] * (1.0f / 64.0f);
	}

	for (int i = 0; i < partition_count; i++)
	{
		vfloat4 rgba_sum(1e-17f);
		vfloat4 rgba_weight_sum(1e-17f);

		int texelcount = pt->partition_texel_count[i];
		const uint8_t *texel_indexes = pt->texels_of_partition[i];

		promise(texelcount > 0);
		for (int j = 0; j < texelcount; j++)
		{
			int tix = texel_indexes[j];

			vfloat4 rgba = blk->texel(tix);
			vfloat4 error_weight(ewb->texel_weight_r[tix], ewb->texel_weight_g[tix], ewb->texel_weight_b[tix], ewb->texel_weight_a[tix]);

			rgba_sum += rgba * error_weight;
			rgba_weight_sum += error_weight;
		}

		vfloat4 scale_direction = normalize((rgba_sum * (1.0f / rgba_weight_sum)).swz<0, 1, 2>());

		float scale_max = 0.0f;
		float scale_min = 1e10f;

		float wmin1 = 1.0f;
		float wmax1 = 0.0f;

		vfloat4 left_sum    = vfloat4::zero();
		vfloat4 middle_sum  = vfloat4::zero();
		vfloat4 right_sum   = vfloat4::zero();
		vfloat4 lmrs_sum    = vfloat4::zero();

		vfloat4 color_vec_x = vfloat4::zero();
		vfloat4 color_vec_y = vfloat4::zero();

		float2 scale_vec = float2(0.0f);

		vfloat4 weight_weight_sum = vfloat4(1e-17f);
		float psum = 1e-17f;

		// FIXME: the loop below has too many responsibilities, making it inefficient.
		for (int j = 0; j < texelcount; j++)
		{
			int tix = texel_indexes[j];

			vfloat4 rgba = blk->texel(tix);
			vfloat4 color_weight(ewb->texel_weight_r[tix], ewb->texel_weight_g[tix], ewb->texel_weight_b[tix], ewb->texel_weight_a[tix]);

			// FIXME: move this calculation out to the color block.
			float ls_weight = hadd_rgb_s(color_weight);
			float idx0 = bilinear_infill(*dt, weight_set, tix);

			float om_idx0 = 1.0f - idx0;
			wmin1 = astc::min(idx0, wmin1);
			wmax1 = astc::max(idx0, wmax1);

			float scale = dot3_s(scale_direction, rgba);
			scale_min = astc::min(scale, scale_min);
			scale_max = astc::max(scale, scale_max);

			vfloat4 left   = color_weight * (om_idx0 * om_idx0);
			vfloat4 middle = color_weight * (om_idx0 * idx0);
			vfloat4 right  = color_weight * (idx0 * idx0);

			vfloat4 lmrs = vfloat4(om_idx0 * om_idx0,
			                       om_idx0 * idx0,
			                       idx0 * idx0,
			                       0.0f) * ls_weight;

			left_sum   += left;
			middle_sum += middle;
			right_sum  += right;
			lmrs_sum   += lmrs;

			vfloat4 color_idx(idx0);
			vfloat4 cwprod = color_weight * rgba;
			vfloat4 cwiprod = cwprod * color_idx;

			color_vec_y += cwiprod;
			color_vec_x += cwprod - cwiprod;

			scale_vec += (float2(om_idx0, idx0) * (ls_weight * scale));
			weight_weight_sum += color_weight * color_idx;
			psum += dot3_s(color_weight * color_idx, color_idx);
		}

		// calculations specific to mode #7, the HDR RGB-scale mode.
		// FIXME: Can we skip this for LDR textures?
		vfloat4 rgbq_sum = color_vec_x + color_vec_y;
		rgbq_sum.set_lane<3>(hadd_rgb_s(color_vec_y));

		#ifdef DEBUG_CAPTURE_NAN
		    fedisableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif

		vfloat4 rgbovec = compute_rgbovec(rgba_weight_sum, weight_weight_sum,
		                                  rgbq_sum, psum);
		rgbo_vectors[i] = rgbovec;

		// We will occasionally get a failure due to the use of a singular
		// (non-invertible) matrix. Record whether such a failure has taken
		// place; if it did, compute rgbo_vectors[] with a different method
		// later on.
		float chkval = dot_s(rgbovec, rgbovec);
		int rgbo_fail = chkval != chkval;

		// Initialize the luminance and scale vectors with a reasonable
		//  default, just in case the subsequent calculation blows up.
		#ifdef DEBUG_CAPTURE_NAN
			fedisableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif

		float scalediv = scale_min * (1.0f / astc::max(scale_max, 1e-10f));
		scalediv = astc::clamp1f(scalediv);

		#ifdef DEBUG_CAPTURE_NAN
			feenableexcept(FE_DIVBYZERO | FE_INVALID);
		#endif

		vfloat4 sds = scale_direction * scale_max;

		rgbs_vectors[i] = vfloat4(sds.lane<0>(), sds.lane<1>(), sds.lane<2>(), scalediv);

		if (wmin1 >= wmax1 * 0.999f)
		{
			// if all weights in the partition were equal, then just take average
			// of all colors in the partition and use that as both endpoint colors.
			vfloat4 avg = (color_vec_x + color_vec_y) * (1.0f / rgba_weight_sum);

			vmask4 notnan_mask = avg == avg;
			ep->endpt0[i] = select(ep->endpt0[i], avg, notnan_mask);
			ep->endpt1[i] = select(ep->endpt1[i], avg, notnan_mask);

			rgbs_vectors[i] = vfloat4(sds.lane<0>(), sds.lane<1>(), sds.lane<2>(), 1.0f);
		}
		else
		{
			// otherwise, complete the analytic calculation of ideal-endpoint-values
			// for the given set of texel weights and pixel colors.

			#ifdef DEBUG_CAPTURE_NAN
			    fedisableexcept(FE_DIVBYZERO | FE_INVALID);
			#endif

			vfloat4 color_det1 = (left_sum * right_sum) - (middle_sum * middle_sum);
			vfloat4 color_rdet1 = 1.0f / color_det1;

			float ls_det1  = (lmrs_sum.lane<0>() * lmrs_sum.lane<2>()) - (lmrs_sum.lane<1>() * lmrs_sum.lane<1>());
			float ls_rdet1 = 1.0f / ls_det1;

			vfloat4 color_mss1 = (left_sum * left_sum)
			                   + (2.0f * middle_sum * middle_sum)
			                   + (right_sum * right_sum);

			float ls_mss1 = (lmrs_sum.lane<0>() * lmrs_sum.lane<0>())
			              + (2.0f * lmrs_sum.lane<1>() * lmrs_sum.lane<1>())
			              + (lmrs_sum.lane<2>() * lmrs_sum.lane<2>());

			vfloat4 ep0 = (right_sum * color_vec_x - middle_sum * color_vec_y) * color_rdet1;
			vfloat4 ep1 = (left_sum * color_vec_y - middle_sum * color_vec_x) * color_rdet1;

			vmask4 det_mask = abs(color_det1) > (color_mss1 * 1e-4f);
			vmask4 notnan_mask = (ep0 == ep0) & (ep1 == ep1);
			vmask4 full_mask = det_mask & notnan_mask;

			ep->endpt0[i] = select(ep->endpt0[i], ep0, full_mask);
			ep->endpt1[i] = select(ep->endpt1[i], ep1, full_mask);

			float scale_ep0 = (lmrs_sum.lane<2>() * scale_vec.r - lmrs_sum.lane<1>() * scale_vec.g) * ls_rdet1;
			float scale_ep1 = (lmrs_sum.lane<0>() * scale_vec.g - lmrs_sum.lane<1>() * scale_vec.r) * ls_rdet1;

			if (fabsf(ls_det1) > (ls_mss1 * 1e-4f) && scale_ep0 == scale_ep0 && scale_ep1 == scale_ep1 && scale_ep0 < scale_ep1)
			{
				float scalediv2 = scale_ep0 * (1.0f / scale_ep1);
				vfloat4 sdsm = scale_direction * scale_ep1;
				rgbs_vectors[i] = vfloat4(sdsm.lane<0>(), sdsm.lane<1>(), sdsm.lane<2>(), scalediv2);
			}

			#ifdef DEBUG_CAPTURE_NAN
				feenableexcept(FE_DIVBYZERO | FE_INVALID);
			#endif
		}

		// if the calculation of an RGB-offset vector failed, try to compute
		// a somewhat-sensible value anyway
		if (rgbo_fail)
		{
			vfloat4 v0 = ep->endpt0[i];
			vfloat4 v1 = ep->endpt1[i];

			float avgdif = hadd_rgb_s(v1 - v0) * (1.0f / 3.0f);
			avgdif = astc::max(avgdif, 0.0f);

			vfloat4 avg = (v0 + v1) * 0.5f;
			vfloat4 ep0 = avg - vfloat4(avgdif) * 0.5f;

			rgbo_vectors[i] = vfloat4(ep0.lane<0>(), ep0.lane<1>(), ep0.lane<2>(), avgdif);
		}
	}
}

#endif
