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

#include <stdio.h>
#include <cassert>

static const float angular_steppings[44] = {
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
	32.0f, 33.0f, 34.0f, 35.0f
};

#define ANGULAR_STEPS ((int)(sizeof(angular_steppings)/sizeof(angular_steppings[0])))

static float stepsizes[ANGULAR_STEPS];
static float stepsizes_sqr[ANGULAR_STEPS];

static int max_angular_steps_needed_for_quant_level[13];

// Store a reduced sin/cos table for 64 possible weight values; this causes
// slight quality loss compared to using sin() and cos() directly. Must be 2^N.
#define SINCOS_STEPS 64

static float sin_table[SINCOS_STEPS][ANGULAR_STEPS];
static float cos_table[SINCOS_STEPS][ANGULAR_STEPS];

void prepare_angular_tables()
{
	int max_angular_steps_needed_for_quant_steps[40];
	for (int i = 0; i < ANGULAR_STEPS; i++)
	{
		stepsizes[i] = 1.0f / angular_steppings[i];
		stepsizes_sqr[i] = stepsizes[i] * stepsizes[i];

		for (int j = 0; j < SINCOS_STEPS; j++)
		{
			sin_table[j][i] = static_cast<float>(sinf((2.0f * (float)M_PI / (SINCOS_STEPS - 1.0f)) * angular_steppings[i] * j));
			cos_table[j][i] = static_cast<float>(cosf((2.0f * (float)M_PI / (SINCOS_STEPS - 1.0f)) * angular_steppings[i] * j));
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
	float anglesum_x[ANGULAR_STEPS];
	float anglesum_y[ANGULAR_STEPS];

	for (int i = 0; i < max_angular_steps; i++)
	{
		anglesum_x[i] = 0;
		anglesum_y[i] = 0;
	}

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

		for (int j = 0; j < max_angular_steps; j++)
		{
			float cp = cosptr[j];
			float sp = sinptr[j];

			anglesum_x[j] += cp * sample_weight;
			anglesum_y[j] += sp * sample_weight;
		}
	}

	// post-process the angle-sums
	for (int i = 0; i < max_angular_steps; i++)
	{
		float angle = astc::atan2(anglesum_y[i], anglesum_x[i]);
		offsets[i] = angle * (stepsizes[i] * (1.0f / (2.0f * (float)M_PI)));
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
// TODO: Add AVX2 version of this; SSE4.2 vectorizes almost perfectly in terms
// of user-visible speedup. Might need to change the length of max_angular
// steps to be a multiple of 8 though ...
#if ASTCENC_SSE >= 42
	// Arrays are always multiple of 4, so this is safe even if overshoot max
	for (int sp = 0; sp < max_angular_steps; sp += 4)
	{
		__m128i minidx = _mm_set1_epi32(128);
		__m128i maxidx = _mm_set1_epi32(-128);
		__m128 errval = _mm_setzero_ps();
		__m128 cut_low_weight_err = _mm_setzero_ps();
		__m128 cut_high_weight_err = _mm_setzero_ps();

		__m128 rcp_stepsize = _mm_load_ps(&angular_steppings[sp]);
		__m128 offset = _mm_load_ps(&offsets[sp]);
		__m128 scaled_offset = _mm_mul_ps(rcp_stepsize, offset);

		for (int j = 0; j < samplecount; j++)
		{
			__m128 wt = _mm_load_ps1(&sample_weights[j]);
			__m128 sval = _mm_mul_ps(_mm_load_ps1(&samples[j]), rcp_stepsize);
			sval = _mm_sub_ps(sval, scaled_offset);
			const int flag = _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC;
			__m128 svalrte = _mm_round_ps(sval, flag);
			__m128i idxv = _mm_cvtps_epi32(svalrte);

			__m128 dif = _mm_sub_ps(sval, svalrte);
			__m128 dwt = _mm_mul_ps(dif, wt);
			errval = _mm_add_ps(errval, _mm_mul_ps(dwt, dif));

			__m128i mask;
			__m128 maskf;
			__m128 accum;

			/* Reset tracker on min hit. */
			mask = _mm_cmplt_epi32(idxv, minidx);
			minidx = _mm_blendv_epi8(minidx, idxv, mask);
			maskf = _mm_castsi128_ps(mask);
			cut_low_weight_err = _mm_blendv_ps(cut_low_weight_err, _mm_setzero_ps(), maskf);

			/* Accumulate on min hit. */
			mask = _mm_cmpeq_epi32(idxv, minidx);
			maskf = _mm_castsi128_ps(mask);
			accum = _mm_add_ps(cut_low_weight_err, _mm_sub_ps(wt, _mm_mul_ps(_mm_set_ps1(2.0f), dwt)));
			cut_low_weight_err = _mm_blendv_ps(cut_low_weight_err, accum, maskf);

			/* Reset tracker on max hit. */
			mask = _mm_cmpgt_epi32(idxv, maxidx);
			maxidx = _mm_blendv_epi8(maxidx, idxv, mask);
			maskf = _mm_castsi128_ps(mask);
			cut_high_weight_err = _mm_blendv_ps(cut_high_weight_err, _mm_setzero_ps(), maskf);

			/* Accumulate on max hit. */
			mask = _mm_cmpeq_epi32(idxv, maxidx);
			maskf = _mm_castsi128_ps(mask);
			accum = _mm_add_ps(cut_high_weight_err, _mm_add_ps(wt, _mm_mul_ps(_mm_set_ps1(2.0f), dwt)));
			cut_high_weight_err = _mm_blendv_ps(cut_high_weight_err, accum, maskf);
		}

		__m128i span = _mm_add_epi32(_mm_sub_epi32(maxidx, minidx), _mm_set1_epi32(1));
		__m128i spanmin = _mm_set1_epi32(2);
		span = _mm_max_epi32(span, spanmin);
		__m128i spanmax = _mm_set1_epi32(max_quantization_steps + 3);
		span = _mm_min_epi32(span, spanmax);

		// Write out min weight and weight span; clamp span to a usable range
		_mm_store_si128((__m128i*)&lowest_weight[sp], minidx);
		_mm_store_si128((__m128i*)&weight_span[sp], span);

		// The cut_(lowest/highest)_weight_error indicate the error that
		// results from  forcing samples that should have had the weight value
		// one step (up/down).
		__m128 errscale = _mm_load_ps(&stepsizes_sqr[sp]);
		_mm_store_ps(&error[sp], _mm_mul_ps(errval, errscale));
		_mm_store_ps(&cut_low_weight_error[sp], _mm_mul_ps(cut_low_weight_err, errscale));
		_mm_store_ps(&cut_high_weight_error[sp], _mm_mul_ps(cut_high_weight_err, errscale));
	}
#else
	for (int sp = 0; sp < max_angular_steps; sp++)
	{
		int minidx = 128;
		int maxidx = -128;
		float errval = 0.0f;
		float cut_low_weight_err = 0.0f;
		float cut_high_weight_err = 0.0f;

		float rcp_stepsize = angular_steppings[sp];
		float offset = offsets[sp];

		float scaled_offset = rcp_stepsize * offset;

		for (int j = 0; j < samplecount; j++)
		{
			float wt = sample_weights[j];
			float sval = (samples[j] * rcp_stepsize) - scaled_offset;
			int idxv = astc::flt2int_rtn(sval);
			float dif = sval - astc::flt2int_rtn(sval);
			float dwt = dif * wt;
			errval += dwt * dif;

			if (idxv < minidx)
			{
				minidx = idxv;
				cut_low_weight_err = wt - 2.0f * dwt;
			}
			else if (idxv == minidx)
			{
				cut_low_weight_err += wt - 2.0f * dwt;
			}

			if (idxv > maxidx)
			{
				maxidx = idxv;
				cut_high_weight_err = wt + 2.0f * dwt;
			}
			else if (idxv == maxidx)
			{
				cut_high_weight_err += wt + 2.0f * dwt;
			}
		}

		// Write out min weight and weight span; clamp span to a usable range
		int span = maxidx - minidx + 1;
		span = MIN(span, max_quantization_steps + 3);
		span = MAX(span, 2);
		lowest_weight[sp] = minidx;
		weight_span[sp] = span;

		// The cut_(lowest/highest)_weight_error indicate the error that
		// results from  forcing samples that should have had the weight value
		// one step (up/down).
		float errscale = stepsizes_sqr[sp];
		error[sp] = errval * errscale;
		cut_low_weight_error[sp] = cut_low_weight_err * errscale;
		cut_high_weight_error[sp] = cut_high_weight_err * errscale;
	}
#endif
}

// main function for running the angular algorithm.
void compute_angular_endpoints_for_quantization_levels(
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

	for (int i = 0; i < MAX_WEIGHT_MODES; i++)
	{
		if (bsd->block_modes[i].is_dual_plane != 0 || bsd->block_modes[i].percentile > mode_cutoff)
		{
			continue;
		}

		int quant_mode = bsd->block_modes[i].quantization_mode;
		int decim_mode = bsd->block_modes[i].decimation_mode;

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

	for (int i = 0; i < MAX_WEIGHT_MODES; i++)
	{
		if (bsd->block_modes[i].is_dual_plane != 1 || bsd->block_modes[i].percentile > mode_cutoff)
		{
			continue;
		}

		int quant_mode = bsd->block_modes[i].quantization_mode;
		int decim_mode = bsd->block_modes[i].decimation_mode;

		low_value1[i] = low_values1[decim_mode][quant_mode];
		high_value1[i] = high_values1[decim_mode][quant_mode];
		low_value2[i] = low_values2[decim_mode][quant_mode];
		high_value2[i] = high_values2[decim_mode][quant_mode];
	}
}

#endif
