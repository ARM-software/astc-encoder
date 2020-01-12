// ----------------------------------------------------------------------------
//  This confidential and proprietary software may be used only as authorised
//  by a licensing agreement from Arm Limited.
//      (C) COPYRIGHT 2011-2020 Arm Limited, ALL RIGHTS RESERVED
//  The entire notice above must be reproduced on all authorised copies and
//  copies may only be made to the extent permitted by a licensing agreement
//  from Arm Limited.
// ----------------------------------------------------------------------------

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

#include "astc_codec_internals.h"

#include <stdio.h>

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

// we store sine/cosine values for 64 possible weight values; this causes
// slight quality loss compared to using sin() and cos() directly.

#define SINCOS_STEPS 64

static float sin_table[SINCOS_STEPS][ANGULAR_STEPS];
static float cos_table[SINCOS_STEPS][ANGULAR_STEPS];

void prepare_angular_tables()
{
	int i, j;
	int max_angular_steps_needed_for_quant_steps[40];
	for (i = 0; i < ANGULAR_STEPS; i++)
	{
		stepsizes[i] = 1.0f / angular_steppings[i];
		stepsizes_sqr[i] = stepsizes[i] * stepsizes[i];

		for (j = 0; j < SINCOS_STEPS; j++)
		{
			sin_table[j][i] = static_cast<float>(sinf((2.0f * (float)M_PI / (SINCOS_STEPS - 1.0f)) * angular_steppings[i] * j));
			cos_table[j][i] = static_cast<float>(cosf((2.0f * (float)M_PI / (SINCOS_STEPS - 1.0f)) * angular_steppings[i] * j));
		}

		int p = static_cast < int >(floor(angular_steppings[i])) + 1;
		max_angular_steps_needed_for_quant_steps[p] = MIN(i + 1, ANGULAR_STEPS - 1);
	}

	// yes, the next-to-last entry is supposed to have the value 33. This because under
	// ASTC, the 32-weight mode leaves a double-sized hole in the middle of the
	// weight space, so we are better off matching 33 weights than 32.
	static const int steps_of_level[] = { 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 33, 36 };

	for (i = 0; i < 13; i++)
		max_angular_steps_needed_for_quant_level[i] = max_angular_steps_needed_for_quant_steps[steps_of_level[i]];
}

// function to compute angular sums; then, from the
// angular sums, compute alignment factor and offset.

void compute_angular_offsets(
	int samplecount,
	const float* samples,
	const float* sample_weights,
	int max_angular_steps,
	float* offsets
) {
	int i, j;

	float anglesum_x[ANGULAR_STEPS];
	float anglesum_y[ANGULAR_STEPS];

	for (i = 0; i < max_angular_steps; i++)
	{
		anglesum_x[i] = 0;
		anglesum_y[i] = 0;
	}

	// compute the angle-sums.
	for (i = 0; i < samplecount; i++)
	{
		float sample = samples[i];
		float sample_weight = sample_weights[i];
		if32 p;
		p.f = (sample * (SINCOS_STEPS - 1.0f)) + 12582912.0f;
		unsigned int isample = p.u & 0x3F;

		const float *sinptr = sin_table[isample];
		const float *cosptr = cos_table[isample];

		for (j = 0; j < max_angular_steps; j++)
		{
			float cp = cosptr[j];
			float sp = sinptr[j];

			anglesum_x[j] += cp * sample_weight;
			anglesum_y[j] += sp * sample_weight;
		}
	}

	// post-process the angle-sums
	for (i = 0; i < max_angular_steps; i++)
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
	int8_t * lowest_weight,
	int8_t * weight_span,
	float *error,
	float *cut_low_weight_error,
	float *cut_high_weight_error
) {
	// weight + 12
	static const unsigned int idxtab[128] = {
		12, 13, 14, 15, 16, 17, 18, 19,
		20, 21, 22, 23, 24, 25, 26, 27,
		28, 29, 30, 31, 32, 33, 34, 35,
		36, 37, 38, 39, 40, 41, 42, 43,
		44, 45, 46, 47, 48, 49, 50, 51,
		52, 53, 54, 55, 55, 55, 55, 55,
		55, 55, 55, 55, 55, 55, 55, 55,
		55, 55, 55, 55, 55, 55, 55, 55,
		 0,  0,  0,  0,  0,  0,  0,  0,
		 0,  0,  0,  0,  0,  0,  0,  0,
		 0,  0,  0,  0,  0,  0,  0,  0,
		 0,  0,  0,  0,  0,  0,  0,  0,
		 0,  0,  0,  0,  0,  0,  0,  0,
		 0,  0,  0,  0,  0,  0,  0,  0,
		 0,  0,  0,  0,  0,  1,  2,  3,
		 4,  5,  6,  7,  8,  9, 10, 11
	};

	for (int sp = 0; sp < max_angular_steps; sp++)
	{
		unsigned int minidx_bias12 = 55;
		unsigned int maxidx_bias12 = 0;

		float errval = 0.0f;
		float cut_low_weight_err = 0.0f;
		float cut_high_weight_err = 0.0f;

		float rcp_stepsize = angular_steppings[sp];
		float offset = offsets[sp];

		float scaled_offset = rcp_stepsize * offset;

		for (int j = 0; j < samplecount; j++)
		{
			float wt = sample_weights[j];
			if32 p;
			float sval = (samples[j] * rcp_stepsize) - scaled_offset;
			p.f = sval + 12582912.0f;	// FP representation abuse to avoid floor() and float->int conversion
			float isval = p.f - 12582912.0f;

			float dif = sval - isval;
			float dwt = dif * wt;
			errval += dwt * dif;

			unsigned int idx_bias12 = idxtab[p.u & 0x7F];

			if (idx_bias12 < minidx_bias12)
			{
				minidx_bias12 = idx_bias12;
				cut_low_weight_err = wt - 2.0f * dwt;
			}
			else if (idx_bias12 == minidx_bias12)
			{
				cut_low_weight_err += wt - 2.0f * dwt;
			}

			if (idx_bias12 > maxidx_bias12)
			{
				maxidx_bias12 = idx_bias12;
				cut_high_weight_err = wt + 2.0f * dwt;
			}
			else if (idx_bias12 == maxidx_bias12)
			{
				cut_high_weight_err += wt + 2.0f * dwt;
			}
		}

		int minIndex = minidx_bias12 - 12;
		int maxIndex = maxidx_bias12 - 12;
		int span = maxIndex - minIndex + 1;

		// Clamp the span to a usable range
		span = MIN(span, max_quantization_steps + 3);
		span = MAX(span, 2);
		lowest_weight[sp] = (int)minidx_bias12 - 12;
		weight_span[sp] = span;
		error[sp] = errval;

		// the cut_(lowest/highest)_weight_error indicate the error that results from
		// forcing samples that should have had the (lowest/highest) weight value
		// one step (up/down).
		cut_low_weight_error[sp] = cut_low_weight_err;
		cut_high_weight_error[sp] = cut_high_weight_err;
	}

	for (int sp = 0; sp < max_angular_steps; sp++)
	{
		float errscale = stepsizes_sqr[sp];
		error[sp] *= errscale;
		cut_low_weight_error[sp] *= errscale;
		cut_high_weight_error[sp] *= errscale;
	}
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
	int i;

	static const int quantization_steps_for_level[13] = { 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 33, 36 };

	int max_quantization_steps = quantization_steps_for_level[max_quantization_level + 1];

	float angular_offsets[ANGULAR_STEPS];
	int max_angular_steps = max_angular_steps_needed_for_quant_level[max_quantization_level];
	compute_angular_offsets(samplecount, samples, sample_weights, max_angular_steps, angular_offsets);

	// the +4 offsets are to allow for vectorization within compute_lowest_and_highest_weight().
	int8_t lowest_weight[ANGULAR_STEPS + 4];
	int8_t weight_span[ANGULAR_STEPS + 4];
	float error[ANGULAR_STEPS + 4];

	float cut_low_weight_error[ANGULAR_STEPS + 4];
	float cut_high_weight_error[ANGULAR_STEPS + 4];

	compute_lowest_and_highest_weight(samplecount, samples, sample_weights,
	                                  max_angular_steps, max_quantization_steps,
	                                  angular_offsets, lowest_weight, weight_span, error,
	                                  cut_low_weight_error, cut_high_weight_error);

	// for each quantization level, find the best error terms.
	float best_errors[40];
	int best_scale[40];
	uint8_t cut_low_weight[40];
	for (i = 0; i < (max_quantization_steps + 4); i++)
	{
		best_errors[i] = 1e30f;
		best_scale[i] = -1;	// Indicates no solution found
		cut_low_weight[i] = 0;
	}

	for (i = 0; i < max_angular_steps; i++)
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
	for (i = 3; i <= max_quantization_steps; i++)
	{
		if (best_errors[i] > best_errors[i - 1])
		{
			best_errors[i] = best_errors[i - 1];
			best_scale[i] = best_scale[i - 1];
			cut_low_weight[i] = cut_low_weight[i - 1];
		}
	}

	for (i = 0; i <= max_quantization_level; i++)
	{
		int q = quantization_steps_for_level[i];
		int bsi = best_scale[q];

		// Did we find anything?
		if (bsi < 0)
		{
			printf("ERROR: Unable to find an encoding within the specified error limits.\n");
			exit(1);
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
	int i;
	float low_values[MAX_DECIMATION_MODES][12];
	float high_values[MAX_DECIMATION_MODES][12];

	for (i = 0; i < MAX_DECIMATION_MODES; i++)
	{
		// TODO: Do this at build time and cache the result
		int samplecount = bsd->decimation_mode_samples[i];
		int quant_mode = bsd->decimation_mode_maxprec_1plane[i];
		float percentile = bsd->decimation_mode_percentile[i];
		int permit_encode = bsd->permit_encode[i];
		if (permit_encode == 0 || samplecount < 1 || quant_mode < 0 || percentile > mode_cutoff)
			continue;

		compute_angular_endpoints_for_quantization_levels(samplecount,
														  decimated_quantized_weights + i * MAX_WEIGHTS_PER_BLOCK,
														  decimated_weights + i * MAX_WEIGHTS_PER_BLOCK, quant_mode, low_values[i], high_values[i]);
	}

	for (i = 0; i < MAX_WEIGHT_MODES; i++)
	{
		if (bsd->block_modes[i].is_dual_plane != 0 || bsd->block_modes[i].percentile > mode_cutoff)
			continue;
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
	int i;
	float low_values1[MAX_DECIMATION_MODES][12];
	float high_values1[MAX_DECIMATION_MODES][12];
	float low_values2[MAX_DECIMATION_MODES][12];
	float high_values2[MAX_DECIMATION_MODES][12];

	for (i = 0; i < MAX_DECIMATION_MODES; i++)
	{
		// TODO: Do this at build time and cache the result
		int samplecount = bsd->decimation_mode_samples[i];
		int quant_mode = bsd->decimation_mode_maxprec_2planes[i];
		float percentile = bsd->decimation_mode_percentile[i];
		int permit_encode = bsd->permit_encode[i];
		if (permit_encode == 0 || samplecount < 1 || quant_mode < 0 || percentile > mode_cutoff)
			continue;

		compute_angular_endpoints_for_quantization_levels(samplecount,
														  decimated_quantized_weights + 2 * i * MAX_WEIGHTS_PER_BLOCK,
														  decimated_weights + 2 * i * MAX_WEIGHTS_PER_BLOCK, quant_mode, low_values1[i], high_values1[i]);

		compute_angular_endpoints_for_quantization_levels(samplecount,
														  decimated_quantized_weights + (2 * i + 1) * MAX_WEIGHTS_PER_BLOCK,
														  decimated_weights + (2 * i + 1) * MAX_WEIGHTS_PER_BLOCK, quant_mode, low_values2[i], high_values2[i]);
	}

	for (i = 0; i < MAX_WEIGHT_MODES; i++)
	{
		if (bsd->block_modes[i].is_dual_plane != 1 || bsd->block_modes[i].percentile > mode_cutoff)
			continue;
		int quant_mode = bsd->block_modes[i].quantization_mode;
		int decim_mode = bsd->block_modes[i].decimation_mode;

		low_value1[i] = low_values1[decim_mode][quant_mode];
		high_value1[i] = high_values1[decim_mode][quant_mode];
		low_value2[i] = low_values2[decim_mode][quant_mode];
		high_value2[i] = high_values2[decim_mode][quant_mode];
	}
}
