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
 * @brief Functions for angular-sum algorithm for weight alignment.
 *
 * This algorithm works as follows:
 * - we compute a complex number P as (cos s*i, sin s*i) for each weight,
 *   where i is the input value and s is a scaling factor based on the spacing
 *   between the weights.
 * - we then add together complex numbers for all the weights.
 * - we then compute the length and angle of the resulting sum.
 *
 * This should produce the following results:
 * - perfect alignment results in a vector whose length is equal to the sum of
 *   lengths of all inputs
 * - even distribution results in a vector of length 0.
 * - all samples identical results in perfect alignment for every scaling.
 *
 * For each scaling factor within a given set, we compute an alignment factor
 * from 0 to 1. This should then result in some scalings standing out as having
 * particularly good alignment factors; we can use this to produce a set of
 * candidate scale/shift values for various quantization levels; we should then
 * actually try them and see what happens.
 *
 * Assuming N quantization steps, the scaling factor becomes s=2*PI*(N-1); we
 * should probably have about 1 scaling factor for every 1/4 quantization step
 * (perhaps 1/8 for low levels of quantization).
 */

#include "astcenc_internal.h"
#include "astcenc_vecmathlib.h"

#include <stdio.h>
#include <cassert>
#include <cstring>

#if ASTCENC_SIMD_WIDTH <= 4
    #define ANGULAR_STEPS 44
#elif ASTCENC_SIMD_WIDTH == 8
    // AVX code path loops over these tables 8 elements at a time,
    // so make sure to have their size a multiple of 8.
    #define ANGULAR_STEPS 48
#else
    #error Unknown SIMD width
#endif
static_assert((ANGULAR_STEPS % ASTCENC_SIMD_WIDTH) == 0, "ANGULAR_STEPS should be multiple of ASTCENC_SIMD_WIDTH");

alignas(ASTCENC_VECALIGN) static const float angular_steppings[ANGULAR_STEPS] = {
	 1.0f, 1.25f, 1.5f, 1.75f,

	 2.0f,  2.5f, 3.0f, 3.5f,
	 4.0f,  4.5f, 5.0f, 5.5f,
	 6.0f,  6.5f, 7.0f, 7.5f,

	 8.0f,  9.0f, 10.0f, 11.0f,
	12.0f, 13.0f, 14.0f, 15.0f,
	16.0f, 17.0f, 18.0f, 19.0f,
	20.0f, 21.0f, 22.0f, 23.0f,
	24.0f, 25.0f, 26.0f, 27.0f,
	28.0f, 29.0f, 30.0f, 31.0f,
	32.0f, 33.0f, 34.0f, 35.0f,
#if ANGULAR_STEPS >= 48
    // This is "redundant" and only used in more-than-4-wide
    // SIMD code paths, to make the steps table size
    // be a multiple of SIMD width. Values are replicated
    // from last entry so that AVX2 and SSE code paths
    // return the same results.
    35.0f, 35.0f, 35.0f, 35.0f,
#endif
};

alignas(ASTCENC_VECALIGN) static float stepsizes[ANGULAR_STEPS];
alignas(ASTCENC_VECALIGN) static float stepsizes_sqr[ANGULAR_STEPS];

static int max_angular_steps_needed_for_quant_level[13];

// Store a reduced sin/cos table for 64 possible weight values; this causes
// slight quality loss compared to using sin() and cos() directly. Must be 2^N.
#define SINCOS_STEPS 64

alignas(ASTCENC_VECALIGN) static float sin_table[SINCOS_STEPS][ANGULAR_STEPS];
alignas(ASTCENC_VECALIGN) static float cos_table[SINCOS_STEPS][ANGULAR_STEPS];

void prepare_angular_tables()
{
	int max_angular_steps_needed_for_quant_steps[40];
	for (int i = 0; i < ANGULAR_STEPS; i++)
	{
		stepsizes[i] = 1.0f / angular_steppings[i];
		stepsizes_sqr[i] = stepsizes[i] * stepsizes[i];

		for (int j = 0; j < SINCOS_STEPS; j++)
		{
			sin_table[j][i] = static_cast<float>(sinf((2.0f * astc::PI / (SINCOS_STEPS - 1.0f)) * angular_steppings[i] * j));
			cos_table[j][i] = static_cast<float>(cosf((2.0f * astc::PI / (SINCOS_STEPS - 1.0f)) * angular_steppings[i] * j));
		}

		int p = astc::flt2int_rd(angular_steppings[i]) + 1;
		max_angular_steps_needed_for_quant_steps[p] = MIN(i + 1, ANGULAR_STEPS - 1);
	}

	// yes, the next-to-last entry is supposed to have the value 33. This because under
	// ASTC, the 32-weight mode leaves a double-sized hole in the middle of the
	// weight space, so we are better off matching 33 weights than 32.
	static const int steps_of_level[] = { 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 33, 36 };

	for (int i = 0; i < 13; i++)
	{
		max_angular_steps_needed_for_quant_level[i] = max_angular_steps_needed_for_quant_steps[steps_of_level[i]];
	}
}

// function to compute angular sums; then, from the
// angular sums, compute alignment factor and offset.

static void compute_angular_offsets(
	int samplecount,
	const float* samples,
	const float* sample_weights,
	int max_angular_steps,
	float* offsets
) {
	alignas(ASTCENC_VECALIGN) float anglesum_x[ANGULAR_STEPS];
	alignas(ASTCENC_VECALIGN) float anglesum_y[ANGULAR_STEPS];
	std::memset(anglesum_x, 0, max_angular_steps*sizeof(anglesum_x[0]));
	std::memset(anglesum_y, 0, max_angular_steps*sizeof(anglesum_y[0]));

	// compute the angle-sums.
	for (int i = 0; i < samplecount; i++)
	{
		float sample = samples[i];
		float sample_weight = sample_weights[i];
		if32 p;
		p.f = (sample * (SINCOS_STEPS - 1.0f)) + 12582912.0f;
		unsigned int isample = p.u & (SINCOS_STEPS - 1);

		const float *sinptr = sin_table[isample];
		const float *cosptr = cos_table[isample];

		vfloat sample_weightv(sample_weight);
		for (int j = 0; j < max_angular_steps; j += ASTCENC_SIMD_WIDTH) // arrays are multiple of SIMD width (ANGULAR_STEPS), safe to overshoot max
		{
			vfloat cp = loada(&cosptr[j]);
			vfloat sp = loada(&sinptr[j]);
			vfloat ax = loada(&anglesum_x[j]) + cp * sample_weightv;
			vfloat ay = loada(&anglesum_y[j]) + sp * sample_weightv;
			storea(ax, &anglesum_x[j]);
			storea(ay, &anglesum_y[j]);
		}
	}

	// post-process the angle-sums
	vfloat mult = vfloat(1.0f / (2.0f * astc::PI));
	for (int i = 0; i < max_angular_steps; i += ASTCENC_SIMD_WIDTH) // arrays are multiple of SIMD width (ANGULAR_STEPS), safe to overshoot max
	{
		vfloat angle = atan2(loada(&anglesum_y[i]), loada(&anglesum_x[i]));
		vfloat ofs = angle * (loada(&stepsizes[i]) * mult);
		storea(ofs, &offsets[i]);
	}
}

// for a given step-size and a given offset, compute the
// lowest and highest weight that results from quantizing using the stepsize & offset.
// also, compute the resulting error.

static void compute_lowest_and_highest_weight(
	int samplecount,
	const float *samples,
	const float *sample_weights,
	int max_angular_steps,
	int max_quantization_steps,
	const float *offsets,
	int32_t * lowest_weight,
	int32_t * weight_span,
	float *error,
	float *cut_low_weight_error,
	float *cut_high_weight_error
) {
	// Arrays are always multiple of SIMD width (ANGULAR_STEPS), so this is safe even if overshoot max
	for (int sp = 0; sp < max_angular_steps; sp += ASTCENC_SIMD_WIDTH)
	{
		vint minidx(128);
		vint maxidx(-128);
		vfloat errval = vfloat::zero();
		vfloat cut_low_weight_err = vfloat::zero();
		vfloat cut_high_weight_err = vfloat::zero();
		vfloat rcp_stepsize = loada(&angular_steppings[sp]);
		vfloat offset = loada(&offsets[sp]);
		vfloat scaled_offset = rcp_stepsize * offset;
		for (int j = 0; j < samplecount; ++j)
		{
			vfloat wt = load1(&sample_weights[j]);
			vfloat sval = load1(&samples[j]) * rcp_stepsize - scaled_offset;
			vfloat svalrte = round(sval);
			vint idxv = float_to_int(svalrte);
			vfloat dif = sval - svalrte;
			vfloat dwt = dif * wt;
			errval = errval + dwt * dif;

			// Reset tracker on min hit.
			vmask mask = idxv < minidx;
			minidx = select(minidx, idxv, mask);
			cut_low_weight_err = select(cut_low_weight_err, vfloat::zero(), mask);

			// Accumulate on min hit.
			mask = idxv == minidx;
			minidx = select(minidx, idxv, mask);
			vfloat accum = cut_low_weight_err + wt - vfloat(2.0f) * dwt;
			cut_low_weight_err = select(cut_low_weight_err, accum, mask);

			// Reset tracker on max hit.
			mask = idxv > maxidx;
			maxidx = select(maxidx, idxv, mask);
			cut_high_weight_err = select(cut_high_weight_err, vfloat::zero(), mask);

			// Accumulate on max hit.
			mask = idxv == maxidx;
			accum = cut_high_weight_err + wt + vfloat(2.0f) * dwt;
			cut_high_weight_err = select(cut_high_weight_err, accum, mask);
		}

		// Write out min weight and weight span; clamp span to a usable range
		vint span = maxidx - minidx + vint(1);
		span = min(span, vint(max_quantization_steps + 3));
		span = max(span, vint(2));
		storea(minidx, &lowest_weight[sp]);
		storea(span, &weight_span[sp]);

		// The cut_(lowest/highest)_weight_error indicate the error that
		// results from  forcing samples that should have had the weight value
		// one step (up/down).
		vfloat errscale = loada(&stepsizes_sqr[sp]);
		storea(errval * errscale, &error[sp]);
		storea(cut_low_weight_err * errscale, &cut_low_weight_error[sp]);
		storea(cut_high_weight_err * errscale, &cut_high_weight_error[sp]);
	}
}

// main function for running the angular algorithm.
static void compute_angular_endpoints_for_quantization_levels(
	int samplecount,
	const float* samples,
	const float* sample_weights,
	int max_quantization_level,
	float low_value[12],
	float high_value[12]
) {
	static const int quantization_steps_for_level[13] = { 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 33, 36 };

	int max_quantization_steps = quantization_steps_for_level[max_quantization_level + 1];

	alignas(ASTCENC_VECALIGN) float angular_offsets[ANGULAR_STEPS];
	int max_angular_steps = max_angular_steps_needed_for_quant_level[max_quantization_level];
	compute_angular_offsets(samplecount, samples, sample_weights, max_angular_steps, angular_offsets);

	alignas(ASTCENC_VECALIGN) int32_t lowest_weight[ANGULAR_STEPS];
	alignas(ASTCENC_VECALIGN) int32_t weight_span[ANGULAR_STEPS];
	alignas(ASTCENC_VECALIGN) float error[ANGULAR_STEPS];
	alignas(ASTCENC_VECALIGN) float cut_low_weight_error[ANGULAR_STEPS];
	alignas(ASTCENC_VECALIGN) float cut_high_weight_error[ANGULAR_STEPS];

	compute_lowest_and_highest_weight(samplecount, samples, sample_weights,
	                                  max_angular_steps, max_quantization_steps,
	                                  angular_offsets, lowest_weight, weight_span, error,
	                                  cut_low_weight_error, cut_high_weight_error);

	// for each quantization level, find the best error terms.
	float best_errors[40];
	int best_scale[40];
	uint8_t cut_low_weight[40];
	for (int i = 0; i < (max_quantization_steps + 4); i++)
	{
		best_scale[i] = -1;	// Indicates no solution found
		best_errors[i] = 1e30f;
		cut_low_weight[i] = 0;
	}

	for (int i = 0; i < max_angular_steps; i++)
	{
		int idx_span = weight_span[i];

		if (best_errors[idx_span] > error[i])
		{
			best_errors[idx_span] = error[i];
			best_scale[idx_span] = i;
			cut_low_weight[idx_span] = 0;
		}

		float error_cut_low = error[i] + cut_low_weight_error[i];
		float error_cut_high = error[i] + cut_high_weight_error[i];
		float error_cut_low_high = error[i] + cut_low_weight_error[i] + cut_high_weight_error[i];

		if (best_errors[idx_span - 1] > error_cut_low)
		{
			best_errors[idx_span - 1] = error_cut_low;
			best_scale[idx_span - 1] = i;
			cut_low_weight[idx_span - 1] = 1;
		}

		if (best_errors[idx_span - 1] > error_cut_high)
		{
			best_errors[idx_span - 1] = error_cut_high;
			best_scale[idx_span - 1] = i;
			cut_low_weight[idx_span - 1] = 0;
		}

		if (best_errors[idx_span - 2] > error_cut_low_high)
		{
			best_errors[idx_span - 2] = error_cut_low_high;
			best_scale[idx_span - 2] = i;
			cut_low_weight[idx_span - 2] = 1;
		}

	}

	// if we got a better error-value for a low sample count than for a high one,
	// use the low sample count error value for the higher sample count as well.
	for (int i = 3; i <= max_quantization_steps; i++)
	{
		if (best_errors[i] > best_errors[i - 1])
		{
			best_errors[i] = best_errors[i - 1];
			best_scale[i] = best_scale[i - 1];
			cut_low_weight[i] = cut_low_weight[i - 1];
		}
	}

	for (int i = 0; i <= max_quantization_level; i++)
	{
		int q = quantization_steps_for_level[i];
		int bsi = best_scale[q];

		// Did we find anything?
		// TODO: Can we do better than bsi = 0 here. We should at least
		// propagate an error (and move the printf into the CLI).
		if (bsi < 0)
		{
			printf("WARNING: Unable to find encoding within specified error limit\n");
			bsi = 0;
		}

		float stepsize = stepsizes[bsi];
		int lwi = lowest_weight[bsi] + cut_low_weight[q];
		int hwi = lwi + q - 1;
		float offset = angular_offsets[bsi];

		low_value[i] = offset + lwi * stepsize;
		high_value[i] = offset + hwi * stepsize;
	}
}

// helper functions that will compute ideal angular-endpoints
// for a given set of weights and a given block size descriptors
void compute_angular_endpoints_1plane(
	float mode_cutoff,
	const block_size_descriptor* bsd,
	const float* decimated_quantized_weights,
	const float* decimated_weights,
	float low_value[MAX_WEIGHT_MODES],
	float high_value[MAX_WEIGHT_MODES]
) {
	float low_values[MAX_DECIMATION_MODES][12];
	float high_values[MAX_DECIMATION_MODES][12];

	for (int i = 0; i < MAX_DECIMATION_MODES; i++)
	{
		// TODO: Do this at build time and cache the result
		int samplecount = bsd->decimation_mode_samples[i];
		int quant_mode = bsd->decimation_mode_maxprec_1plane[i];
		float percentile = bsd->decimation_mode_percentile[i];
		int permit_encode = bsd->permit_encode[i];
		if (permit_encode == 0 || samplecount < 1 || quant_mode < 0 || percentile > mode_cutoff)
		{
			continue;
		}

		compute_angular_endpoints_for_quantization_levels(samplecount,
		                                                  decimated_quantized_weights + i * MAX_WEIGHTS_PER_BLOCK,
		                                                  decimated_weights + i * MAX_WEIGHTS_PER_BLOCK, quant_mode, low_values[i], high_values[i]);
	}

	for (int i = 0, ni = bsd->block_mode_packed_count; i < ni; ++i)
	{
		const block_mode& bm = bsd->block_modes_packed[i];
		if (bm.is_dual_plane != 0 || bm.percentile > mode_cutoff)
		{
			continue;
		}

		int quant_mode = bm.quantization_mode;
		int decim_mode = bm.decimation_mode;

		low_value[i] = low_values[decim_mode][quant_mode];
		high_value[i] = high_values[decim_mode][quant_mode];
	}
}

void compute_angular_endpoints_2planes(
	float mode_cutoff,
	const block_size_descriptor* bsd,
	const float* decimated_quantized_weights,
	const float* decimated_weights,
	float low_value1[MAX_WEIGHT_MODES],
	float high_value1[MAX_WEIGHT_MODES],
	float low_value2[MAX_WEIGHT_MODES],
	float high_value2[MAX_WEIGHT_MODES]
) {
	float low_values1[MAX_DECIMATION_MODES][12];
	float high_values1[MAX_DECIMATION_MODES][12];
	float low_values2[MAX_DECIMATION_MODES][12];
	float high_values2[MAX_DECIMATION_MODES][12];

	for (int i = 0; i < MAX_DECIMATION_MODES; i++)
	{
		// TODO: Do this at build time and cache the result
		int samplecount = bsd->decimation_mode_samples[i];
		int quant_mode = bsd->decimation_mode_maxprec_2planes[i];
		float percentile = bsd->decimation_mode_percentile[i];
		int permit_encode = bsd->permit_encode[i];

		if (permit_encode == 0 || samplecount < 1 || quant_mode < 0 || percentile > mode_cutoff)
		{
			continue;
		}

		compute_angular_endpoints_for_quantization_levels(samplecount,
		                                                  decimated_quantized_weights + 2 * i * MAX_WEIGHTS_PER_BLOCK,
		                                                  decimated_weights + 2 * i * MAX_WEIGHTS_PER_BLOCK, quant_mode, low_values1[i], high_values1[i]);

		compute_angular_endpoints_for_quantization_levels(samplecount,
		                                                  decimated_quantized_weights + (2 * i + 1) * MAX_WEIGHTS_PER_BLOCK,
		                                                  decimated_weights + (2 * i + 1) * MAX_WEIGHTS_PER_BLOCK, quant_mode, low_values2[i], high_values2[i]);
	}

	for (int i = 0, ni = bsd->block_mode_packed_count; i < ni; ++i)
	{
		const block_mode& bm = bsd->block_modes_packed[i];
		if (bm.is_dual_plane != 1 || bm.percentile > mode_cutoff)
		{
			continue;
		}

		int quant_mode = bm.quantization_mode;
		int decim_mode = bm.decimation_mode;

		low_value1[i] = low_values1[decim_mode][quant_mode];
		high_value1[i] = high_values1[decim_mode][quant_mode];
		low_value2[i] = low_values2[decim_mode][quant_mode];
		high_value2[i] = high_values2[decim_mode][quant_mode];
	}
}

#endif
