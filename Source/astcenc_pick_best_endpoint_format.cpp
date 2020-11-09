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
 * @brief Functions to pick best ASTC endpoint for a block.
 */

#include "astcenc_internal.h"
#include "astcenc_vecmathlib.h"

#include <assert.h>

/*
   functions to determine, for a given partitioning, which color endpoint formats are the best to use.
 */

// for a given partition, compute for every (integer-component-count, quantization-level)
// the color error.
static void compute_color_error_for_every_integer_count_and_quantization_level(
	int encode_hdr_rgb,	// 1 = perform HDR encoding, 0 = perform LDR encoding.
	int encode_hdr_alpha,
	int partition_index,
	const partition_info* pi,
	const encoding_choice_errors * eci,	// pointer to the structure for the CURRENT partition.
	 const endpoints * ep,
	 float4 error_weightings[4],
	// arrays to return results back through.
	float best_error[21][4],
	int format_of_choice[21][4]
) {
	int partition_size = pi->texels_per_partition[partition_index];

	static const float baseline_quant_error[21] = {
		(65536.0f * 65536.0f / 18.0f),				// 2 values, 1 step
		(65536.0f * 65536.0f / 18.0f) / (2 * 2),	// 3 values, 2 steps
		(65536.0f * 65536.0f / 18.0f) / (3 * 3),	// 4 values, 3 steps
		(65536.0f * 65536.0f / 18.0f) / (4 * 4),	// 5 values
		(65536.0f * 65536.0f / 18.0f) / (5 * 5),
		(65536.0f * 65536.0f / 18.0f) / (7 * 7),
		(65536.0f * 65536.0f / 18.0f) / (9 * 9),
		(65536.0f * 65536.0f / 18.0f) / (11 * 11),
		(65536.0f * 65536.0f / 18.0f) / (15 * 15),
		(65536.0f * 65536.0f / 18.0f) / (19 * 19),
		(65536.0f * 65536.0f / 18.0f) / (23 * 23),
		(65536.0f * 65536.0f / 18.0f) / (31 * 31),
		(65536.0f * 65536.0f / 18.0f) / (39 * 39),
		(65536.0f * 65536.0f / 18.0f) / (47 * 47),
		(65536.0f * 65536.0f / 18.0f) / (63 * 63),
		(65536.0f * 65536.0f / 18.0f) / (79 * 79),
		(65536.0f * 65536.0f / 18.0f) / (95 * 95),
		(65536.0f * 65536.0f / 18.0f) / (127 * 127),
		(65536.0f * 65536.0f / 18.0f) / (159 * 159),
		(65536.0f * 65536.0f / 18.0f) / (191 * 191),
		(65536.0f * 65536.0f / 18.0f) / (255 * 255)
	};

	float4 ep0 = ep->endpt0[partition_index];
	float4 ep1 = ep->endpt1[partition_index];

	float ep1_min = MIN(MIN(ep1.r, ep1.g), ep1.b);
	ep1_min = MAX(ep1_min, 0.0f);

	float4 error_weight = error_weightings[partition_index];

	float error_weight_rgbsum = error_weight.r + error_weight.g + error_weight.b;

	float range_upper_limit_rgb = encode_hdr_rgb ? 61440.0f : 65535.0f;
	float range_upper_limit_alpha = encode_hdr_alpha ? 61440.0f : 65535.0f;

	// it is possible to get endpoint colors significantly outside [0,upper-limit]
	// even if the input data are safely contained in [0,upper-limit];
	// we need to add an error term for this situation,
	float4 ep0_range_error_high;
	float4 ep1_range_error_high;
	float4 ep0_range_error_low;
	float4 ep1_range_error_low;

	ep0_range_error_high.r = MAX(0.0f, ep0.r - range_upper_limit_rgb);
	ep0_range_error_high.g = MAX(0.0f, ep0.g - range_upper_limit_rgb);
	ep0_range_error_high.b = MAX(0.0f, ep0.b - range_upper_limit_rgb);
	ep0_range_error_high.a = MAX(0.0f, ep0.a - range_upper_limit_alpha);

	ep1_range_error_high.r = MAX(0.0f, ep1.r - range_upper_limit_rgb);
	ep1_range_error_high.g = MAX(0.0f, ep1.g - range_upper_limit_rgb);
	ep1_range_error_high.b = MAX(0.0f, ep1.b - range_upper_limit_rgb);
	ep1_range_error_high.a = MAX(0.0f, ep1.a - range_upper_limit_alpha);

	ep0_range_error_low.r = MIN(0.0f, ep0.r);
	ep0_range_error_low.g = MIN(0.0f, ep0.g);
	ep0_range_error_low.b = MIN(0.0f, ep0.b);
	ep0_range_error_low.a = MIN(0.0f, ep0.a);

	ep1_range_error_low.r = MIN(0.0f, ep1.r);
	ep1_range_error_low.g = MIN(0.0f, ep1.g);
	ep1_range_error_low.b = MIN(0.0f, ep1.b);
	ep1_range_error_low.a = MIN(0.0f, ep1.a);

	float4 sum_range_error =
		(ep0_range_error_low * ep0_range_error_low) +
		(ep1_range_error_low * ep1_range_error_low) +
		(ep0_range_error_high * ep0_range_error_high) +
		(ep1_range_error_high * ep1_range_error_high);
	float rgb_range_error = dot(float3(sum_range_error.r, sum_range_error.g, sum_range_error.b),
	                            float3(error_weight.r, error_weight.g, error_weight.b)) * 0.5f * partition_size;
	float alpha_range_error = sum_range_error.a * error_weight.a * 0.5f * partition_size;

	if (encode_hdr_rgb)
	{

		// collect some statistics
		float af, cf;
		if (ep1.r > ep1.g && ep1.r > ep1.b)
		{
			af = ep1.r;
			cf = ep1.r - ep0.r;
		}
		else if (ep1.g > ep1.b)
		{
			af = ep1.g;
			cf = ep1.g - ep0.g;
		}
		else
		{
			af = ep1.b;
			cf = ep1.b - ep0.b;
		}

		float bf = af - ep1_min;	// estimate of color-component spread in high endpoint color
		float3 prd = float3(ep1.r, ep1.g, ep1.b) - float3(cf, cf, cf);
		float3 pdif = prd - float3(ep0.r, ep0.g, ep0.b);
		// estimate of color-component spread in low endpoint color
		float df = MAX(MAX(fabsf(pdif.r), fabsf(pdif.g)), fabsf(pdif.b));

		int b = (int)bf;
		int c = (int)cf;
		int d = (int)df;

		// determine which one of the 6 submodes is likely to be used in
		// case of an RGBO-mode
		int rgbo_mode = 5;		// 7 bits per component
		// mode 4: 8 7 6
		if (b < 32768 && c < 16384)
		{
			rgbo_mode = 4;
		}

		// mode 3: 9 6 7
		if (b < 8192 && c < 16384)
		{
			rgbo_mode = 3;
		}

		// mode 2: 10 5 8
		if (b < 2048 && c < 16384)
		{
			rgbo_mode = 2;
		}

		// mode 1: 11 6 5
		if (b < 2048 && c < 1024)
		{
			rgbo_mode = 1;
		}

		// mode 0: 11 5 7
		if (b < 1024 && c < 4096)
		{
			rgbo_mode = 0;
		}

		// determine which one of the 9 submodes is likely to be used in
		// case of an RGB-mode.
		int rgb_mode = 8;		// 8 bits per component, except 7 bits for blue

		// mode 0: 9 7 6 7
		if (b < 16384 && c < 8192 && d < 8192)
		{
			rgb_mode = 0;
		}

		// mode 1: 9 8 6 6
		if (b < 32768 && c < 8192 && d < 4096)
		{
			rgb_mode = 1;
		}

		// mode 2: 10 6 7 7
		if (b < 4096 && c < 8192 && d < 4096)
		{
			rgb_mode = 2;
		}

		// mode 3: 10 7 7 6
		if (b < 8192 && c < 8192 && d < 2048)
		{
			rgb_mode = 3;
		}

		// mode 4: 11 8 6 5
		if (b < 8192 && c < 2048 && d < 512)
		{
			rgb_mode = 4;
		}

		// mode 5: 11 6 8 6
		if (b < 2048 && c < 8192 && d < 1024)
		{
			rgb_mode = 5;
		}

		// mode 6: 12 7 7 5
		if (b < 2048 && c < 2048 && d < 256)
		{
			rgb_mode = 6;
		}

		// mode 7: 12 6 7 6
		if (b < 1024 && c < 2048 && d < 512)
		{
			rgb_mode = 7;
		}

		static const float rgbo_error_scales[6] = { 4.0f, 4.0f, 16.0f, 64.0f, 256.0f, 1024.0f };
		static const float rgb_error_scales[9] = { 64.0f, 64.0f, 16.0f, 16.0f, 4.0f, 4.0f, 1.0f, 1.0f, 384.0f };

		float mode7mult = rgbo_error_scales[rgbo_mode] * 0.0015f;	// empirically determined ....
		float mode11mult = rgb_error_scales[rgb_mode] * 0.010f;	// empirically determined ....


		float lum_high = (ep1.r + ep1.g + ep1.b) * (1.0f / 3.0f);
		float lum_low = (ep0.r + ep0.g + ep0.b) * (1.0f / 3.0f);
		float lumdif = lum_high - lum_low;
		float mode23mult = lumdif < 960 ? 4.0f : lumdif < 3968 ? 16.0f : 128.0f;

		mode23mult *= 0.0005f;	// empirically determined ....

		// pick among the available HDR endpoint modes
		for (int i = 0; i < 8; i++)
		{
			best_error[i][3] = 1e30f;
			format_of_choice[i][3] = encode_hdr_alpha ? FMT_HDR_RGBA : FMT_HDR_RGB_LDR_ALPHA;
			best_error[i][2] = 1e30f;
			format_of_choice[i][2] = FMT_HDR_RGB;
			best_error[i][1] = 1e30f;
			format_of_choice[i][1] = FMT_HDR_RGB_SCALE;
			best_error[i][0] = 1e30f;
			format_of_choice[i][0] = FMT_HDR_LUMINANCE_LARGE_RANGE;
		}

		for (int i = 8; i < 21; i++)
		{
			// base_quant_error should depend on the scale-factor that would be used
			// during actual encode of the color value.

			float base_quant_error = baseline_quant_error[i] * partition_size * 1.0f;
			float rgb_quantization_error = error_weight_rgbsum * base_quant_error * 2.0f;
			float alpha_quantization_error = error_weight.a * base_quant_error * 2.0f;
			float rgba_quantization_error = rgb_quantization_error + alpha_quantization_error;

			// for 8 integers, we have two encodings: one with HDR alpha and another one
			// with LDR alpha.

			float full_hdr_rgba_error = rgba_quantization_error + rgb_range_error + alpha_range_error;
			best_error[i][3] = full_hdr_rgba_error;
			format_of_choice[i][3] = encode_hdr_alpha ? FMT_HDR_RGBA : FMT_HDR_RGB_LDR_ALPHA;

			// for 6 integers, we have one HDR-RGB encoding
			float full_hdr_rgb_error = (rgb_quantization_error * mode11mult) + rgb_range_error + eci->alpha_drop_error;
			best_error[i][2] = full_hdr_rgb_error;
			format_of_choice[i][2] = FMT_HDR_RGB;

			// for 4 integers, we have one HDR-RGB-Scale encoding
			float hdr_rgb_scale_error = (rgb_quantization_error * mode7mult) + rgb_range_error + eci->alpha_drop_error + eci->rgb_luma_error;

			best_error[i][1] = hdr_rgb_scale_error;
			format_of_choice[i][1] = FMT_HDR_RGB_SCALE;

			// for 2 integers, we assume luminance-with-large-range
			float hdr_luminance_error = (rgb_quantization_error * mode23mult) + rgb_range_error + eci->alpha_drop_error + eci->luminance_error;
			best_error[i][0] = hdr_luminance_error;
			format_of_choice[i][0] = FMT_HDR_LUMINANCE_LARGE_RANGE;
		}
	}
	else
	{
		for (int i = 0; i < 4; i++)
		{
			best_error[i][3] = 1e30f;
			best_error[i][2] = 1e30f;
			best_error[i][1] = 1e30f;
			best_error[i][0] = 1e30f;

			format_of_choice[i][3] = FMT_RGBA;
			format_of_choice[i][2] = FMT_RGB;
			format_of_choice[i][1] = FMT_RGB_SCALE;
			format_of_choice[i][0] = FMT_LUMINANCE;
		}

		// pick among the available LDR endpoint modes
		for (int i = 4; i < 21; i++)
		{
			float base_quant_error = baseline_quant_error[i] * partition_size * 1.0f;
			float rgb_quantization_error = error_weight_rgbsum * base_quant_error;
			float alpha_quantization_error = error_weight.a * base_quant_error;
			float rgba_quantization_error = rgb_quantization_error + alpha_quantization_error;

			// for 8 integers, the available encodings are:
			// full LDR RGB-Alpha
			float full_ldr_rgba_error = rgba_quantization_error;

			if (eci->can_blue_contract)
			{
				full_ldr_rgba_error *= 0.625f;
			}

			if (eci->can_offset_encode && i <= 18)
			{
				full_ldr_rgba_error *= 0.5f;
			}

			full_ldr_rgba_error += rgb_range_error + alpha_range_error;

			best_error[i][3] = full_ldr_rgba_error;
			format_of_choice[i][3] = FMT_RGBA;

			// for 6 integers, we have:
			// - an LDR-RGB encoding
			// - an RGBS + Alpha encoding (LDR)

			float full_ldr_rgb_error = rgb_quantization_error;

			if (eci->can_blue_contract)
			{
				full_ldr_rgb_error *= 0.5f;
			}

			if (eci->can_offset_encode && i <= 18)
			{
				full_ldr_rgb_error *= 0.25f;
			}

			full_ldr_rgb_error += eci->alpha_drop_error + rgb_range_error;

			float rgbs_alpha_error = rgba_quantization_error + eci->rgb_scale_error + rgb_range_error + alpha_range_error;

			if (rgbs_alpha_error < full_ldr_rgb_error)
			{
				best_error[i][2] = rgbs_alpha_error;
				format_of_choice[i][2] = FMT_RGB_SCALE_ALPHA;
			}
			else
			{
				best_error[i][2] = full_ldr_rgb_error;
				format_of_choice[i][2] = FMT_RGB;
			}

			// for 4 integers, we have a Luminance-Alpha encoding and the RGBS encoding
			float ldr_rgbs_error = rgb_quantization_error + eci->alpha_drop_error + eci->rgb_scale_error + rgb_range_error;

			float lum_alpha_error = rgba_quantization_error + eci->luminance_error + rgb_range_error + alpha_range_error;

			if (ldr_rgbs_error < lum_alpha_error)
			{
				best_error[i][1] = ldr_rgbs_error;
				format_of_choice[i][1] = FMT_RGB_SCALE;
			}
			else
			{
				best_error[i][1] = lum_alpha_error;
				format_of_choice[i][1] = FMT_LUMINANCE_ALPHA;
			}

			// for 2 integers, we have a Luminance-encoding and an Alpha-encoding.
			float luminance_error = rgb_quantization_error + eci->alpha_drop_error + eci->luminance_error + rgb_range_error;

			best_error[i][0] = luminance_error;
			format_of_choice[i][0] = FMT_LUMINANCE;
		}
	}
}

// for 1 partition, find the best combination (one format + a quantization level) for a given bitcount
static void one_partition_find_best_combination_for_bitcount(
	float combined_best_error[21][4],
	int formats_of_choice[21][4],
	int bits_available,
	int* best_quantization_level,
	int* best_formats,
	float* error_of_best_combination
) {
	int best_integer_count = -1;
	float best_integer_count_error = 1e20f;
	for (int i = 0; i < 4; i++)
	{
		// compute the quantization level for a given number of integers and a given number of bits.
		int quantization_level = quantization_mode_table[i + 1][bits_available];

		if (quantization_level == -1)
		{
			continue;			// used to indicate the case where we don't have enough bits to represent a given endpoint format at all.
		}

		if (combined_best_error[quantization_level][i] < best_integer_count_error)
		{
			best_integer_count_error = combined_best_error[quantization_level][i];
			best_integer_count = i;
		}
	}

	int ql = quantization_mode_table[best_integer_count + 1][bits_available];

	*best_quantization_level = ql;
	*error_of_best_combination = best_integer_count_error;
	if (ql >= 0)
	{
		*best_formats = formats_of_choice[ql][best_integer_count];
	}
	else
	{
		*best_formats = FMT_LUMINANCE;
	}
}

// for 2 partitions, find the best format combinations for every (quantization-mode, integer-count) combination
static void two_partitions_find_best_combination_for_every_quantization_and_integer_count(
	float best_error[2][21][4],	// indexed by (partition, quant-level, integer-pair-count-minus-1)
	int format_of_choice[2][21][4],
	float combined_best_error[21][7],	// indexed by (quant-level, integer-pair-count-minus-2)
	int formats_of_choice[21][7][2]
) {
	for (int i = 0; i < 21; i++)
	{
		for (int j = 0; j < 7; j++)
		{
			combined_best_error[i][j] = 1e30f;
		}
	}

	for (int quant = 5; quant < 21; quant++)
	{
		for (int i = 0; i < 4; i++)	// integer-count for first endpoint-pair
		{
			for (int j = 0; j < 4; j++)	// integer-count for second endpoint-pair
			{
				int low2 = MIN(i, j);
				int high2 = MAX(i, j);
				if ((high2 - low2) > 1)
				{
					continue;
				}

				int intcnt = i + j;
				float errorterm = MIN(best_error[0][quant][i] + best_error[1][quant][j], 1e10f);
				if (errorterm <= combined_best_error[quant][intcnt])
				{
					combined_best_error[quant][intcnt] = errorterm;
					formats_of_choice[quant][intcnt][0] = format_of_choice[0][quant][i];
					formats_of_choice[quant][intcnt][1] = format_of_choice[1][quant][j];
				}
			}
		}
	}
}

// for 2 partitions, find the best combination (two formats + a quantization level) for a given bitcount
static void two_partitions_find_best_combination_for_bitcount(
	float combined_best_error[21][7],
	int formats_of_choice[21][7][2],
	int bits_available,
	int* best_quantization_level,
	int* best_quantization_level_mod,
	int* best_formats,
	float* error_of_best_combination
) {
	int best_integer_count = 0;
	float best_integer_count_error = 1e20f;

	for (int integer_count = 2; integer_count <= 8; integer_count++)
	{
		// compute the quantization level for a given number of integers and a given number of bits.
		int quantization_level = quantization_mode_table[integer_count][bits_available];

		if (quantization_level == -1)
		{
			break;				// used to indicate the case where we don't have enough bits to represent a given endpoint format at all.
		}

		float integer_count_error = combined_best_error[quantization_level][integer_count - 2];

		if (integer_count_error < best_integer_count_error)
		{
			best_integer_count_error = integer_count_error;
			best_integer_count = integer_count;
		}
	}

	int ql = quantization_mode_table[best_integer_count][bits_available];
	int ql_mod = quantization_mode_table[best_integer_count][bits_available + 2];

	*best_quantization_level = ql;
	*best_quantization_level_mod = ql_mod;
	*error_of_best_combination = best_integer_count_error;
	if (ql >= 0)
	{
		for (int i = 0; i < 2; i++)
		{
			best_formats[i] = formats_of_choice[ql][best_integer_count - 2][i];
		}
	}
	else
	{
		for (int i = 0; i < 2; i++)
		{
			best_formats[i] = FMT_LUMINANCE;
		}
	}
}

// for 3 partitions, find the best format combinations for every (quantization-mode, integer-count) combination
static void three_partitions_find_best_combination_for_every_quantization_and_integer_count(
	float best_error[3][21][4],	// indexed by (partition, quant-level, integer-count)
	int format_of_choice[3][21][4],
	float combined_best_error[21][10],
	int formats_of_choice[21][10][3]
) {
	for (int i = 0; i < 21; i++)
	{
		for (int j = 0; j < 10; j++)
		{
			combined_best_error[i][j] = 1e30f;
		}
	}

	for (int quant = 5; quant < 21; quant++)
	{
		for (int i = 0; i < 4; i++)	// integer-count for first endpoint-pair
		{
			for (int j = 0; j < 4; j++)	// integer-count for second endpoint-pair
			{
				int low2 = MIN(i, j);
				int high2 = MAX(i, j);
				if ((high2 - low2) > 1)
				{
					continue;
				}

				for (int k = 0; k < 4; k++)	// integer-count for third endpoint-pair
				{
					int low3 = MIN(k, low2);
					int high3 = MAX(k, high2);
					if ((high3 - low3) > 1)
					{
						continue;
					}

					int intcnt = i + j + k;
					float errorterm = MIN(best_error[0][quant][i] + best_error[1][quant][j] + best_error[2][quant][k], 1e10f);
					if (errorterm <= combined_best_error[quant][intcnt])
					{
						combined_best_error[quant][intcnt] = errorterm;
						formats_of_choice[quant][intcnt][0] = format_of_choice[0][quant][i];
						formats_of_choice[quant][intcnt][1] = format_of_choice[1][quant][j];
						formats_of_choice[quant][intcnt][2] = format_of_choice[2][quant][k];
					}
				}
			}
		}
	}
}

// for 3 partitions, find the best combination (three formats + a quantization level) for a given bitcount
static void three_partitions_find_best_combination_for_bitcount(
	float combined_best_error[21][10],
	int formats_of_choice[21][10][3],
	int bits_available,
	int* best_quantization_level,
	int* best_quantization_level_mod,
	int* best_formats,
	float* error_of_best_combination
) {
	int best_integer_count = 0;
	float best_integer_count_error = 1e20f;

	for (int integer_count = 3; integer_count <= 9; integer_count++)
	{
		// compute the quantization level for a given number of integers and a given number of bits.
		int quantization_level = quantization_mode_table[integer_count][bits_available];

		if (quantization_level == -1)
		{
			break;				// used to indicate the case where we don't have enough bits to represent a given endpoint format at all.
		}

		float integer_count_error = combined_best_error[quantization_level][integer_count - 3];

		if (integer_count_error < best_integer_count_error)
		{
			best_integer_count_error = integer_count_error;
			best_integer_count = integer_count;
		}
	}

	int ql = quantization_mode_table[best_integer_count][bits_available];
	int ql_mod = quantization_mode_table[best_integer_count][bits_available + 5];

	*best_quantization_level = ql;
	*best_quantization_level_mod = ql_mod;
	*error_of_best_combination = best_integer_count_error;
	if (ql >= 0)
	{
		for (int i = 0; i < 3; i++)
		{
			best_formats[i] = formats_of_choice[ql][best_integer_count - 3][i];
		}
	}
	else
	{
		for (int i = 0; i < 3; i++)
		{
			best_formats[i] = FMT_LUMINANCE;
		}
	}
}

// for 4 partitions, find the best format combinations for every (quantization-mode, integer-count) combination
static void four_partitions_find_best_combination_for_every_quantization_and_integer_count(
	float best_error[4][21][4],	// indexed by (partition, quant-level, integer-count)
	int format_of_choice[4][21][4],
	float combined_best_error[21][13],
	int formats_of_choice[21][13][4]
) {
	for (int i = 0; i < 21; i++)
	{
		for (int j = 0; j < 13; j++)
		{
			combined_best_error[i][j] = 1e30f;
		}
	}

	for (int quant = 5; quant < 21; quant++)
	{
		for (int i = 0; i < 4; i++)	// integer-count for first endpoint-pair
		{
			for (int j = 0; j < 4; j++)	// integer-count for second endpoint-pair
			{
				int low2 = MIN(i, j);
				int high2 = MAX(i, j);
				if ((high2 - low2) > 1)
				{
					continue;
				}

				for (int k = 0; k < 4; k++)	// integer-count for third endpoint-pair
				{
					int low3 = MIN(k, low2);
					int high3 = MAX(k, high2);
					if ((high3 - low3) > 1)
					{
						continue;
					}

					for (int l = 0; l < 4; l++)	// integer-count for fourth endpoint-pair
					{
						int low4 = MIN(l, low3);
						int high4 = MAX(l, high3);
						if ((high4 - low4) > 1)
						{
							continue;
						}

						int intcnt = i + j + k + l;
						float errorterm = MIN(best_error[0][quant][i] + best_error[1][quant][j] + best_error[2][quant][k] + best_error[3][quant][l], 1e10f);
						if (errorterm <= combined_best_error[quant][intcnt])
						{
							combined_best_error[quant][intcnt] = errorterm;
							formats_of_choice[quant][intcnt][0] = format_of_choice[0][quant][i];
							formats_of_choice[quant][intcnt][1] = format_of_choice[1][quant][j];
							formats_of_choice[quant][intcnt][2] = format_of_choice[2][quant][k];
							formats_of_choice[quant][intcnt][3] = format_of_choice[3][quant][l];
						}
					}
				}
			}
		}
	}
}

// for 4 partitions, find the best combination (four formats + a quantization level) for a given bitcount
static void four_partitions_find_best_combination_for_bitcount(
	float combined_best_error[21][13],
	int formats_of_choice[21][13][4],
	int bits_available,
	int* best_quantization_level,
	int* best_quantization_level_mod,
	int* best_formats,
	float* error_of_best_combination
) {
	int best_integer_count = 0;
	float best_integer_count_error = 1e20f;

	for (int integer_count = 4; integer_count <= 9; integer_count++)
	{
		// compute the quantization level for a given number of integers and a given number of bits.
		int quantization_level = quantization_mode_table[integer_count][bits_available];

		if (quantization_level == -1)
		{
			break;				// used to indicate the case where we don't have enough bits to represent a given endpoint format at all.
		}

		float integer_count_error = combined_best_error[quantization_level][integer_count - 4];

		if (integer_count_error < best_integer_count_error)
		{
			best_integer_count_error = integer_count_error;
			best_integer_count = integer_count;
		}
	}

	int ql = quantization_mode_table[best_integer_count][bits_available];
	int ql_mod = quantization_mode_table[best_integer_count][bits_available + 8];

	*best_quantization_level = ql;
	*best_quantization_level_mod = ql_mod;
	*error_of_best_combination = best_integer_count_error;
	if (ql >= 0)
	{
		for (int i = 0; i < 4; i++)
		{
			best_formats[i] = formats_of_choice[ql][best_integer_count - 4][i];
		}
	}
	else
	{
		for (int i = 0; i < 4; i++)
		{
			best_formats[i] = FMT_LUMINANCE;
		}
	}
}

/*
	The determine_optimal_set_of_endpoint_formats_to_use() function.

	It identifies, for each mode, which set of color endpoint encodings
	produces the best overall result. It then reports back which
	tune_candidate_limit,  modes look best, along with the ideal color encoding
	combination for each.

	It takes as input:
		a partitioning an imageblock,
		a set of color endpoints.
		for each mode, the number of bits available for color encoding and the error incurred by quantization.
		in case of 2 plane of weights, a specifier for which color component to use for the second plane of weights.

	It delivers as output for each of the tune_candidate_limit selected modes:
		format specifier
		for each partition
			quantization level to use
			modified quantization level to use
		(when all format specifiers are equal)
*/
void determine_optimal_set_of_endpoint_formats_to_use(
	const block_size_descriptor* bsd,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	const endpoints* ep,
	int separate_component,	// separate color component for 2-plane mode; -1 for single-plane mode
	 // bitcounts and errors computed for the various quantization methods
	const int* qwt_bitcounts,
	const float* qwt_errors,
	int tune_candidate_limit,
	// output data
	int partition_format_specifiers[TUNE_MAX_TRIAL_CANDIDATES][4],
	int quantized_weight[TUNE_MAX_TRIAL_CANDIDATES],
	int quantization_level[TUNE_MAX_TRIAL_CANDIDATES],
	int quantization_level_mod[TUNE_MAX_TRIAL_CANDIDATES]
) {
	int partition_count = pt->partition_count;

	int encode_hdr_rgb = blk->rgb_lns[0];
	int encode_hdr_alpha = blk->alpha_lns[0];

	// call a helper function to compute the errors that result from various
	// encoding choices (such as using luminance instead of RGB, discarding Alpha,
	// using RGB-scale in place of two separate RGB endpoints and so on)
	encoding_choice_errors eci[4];
	compute_encoding_choice_errors(bsd, blk, pt, ewb, separate_component, eci);

	// for each partition, compute the error weights to apply for that partition.
	float4 error_weightings[4];
	float4 dummied_color_scalefactors[4];	// only used to receive data
	compute_partition_error_color_weightings(bsd, ewb, pt, error_weightings, dummied_color_scalefactors);

	float best_error[4][21][4];
	int format_of_choice[4][21][4];
	for (int i = 0; i < partition_count; i++)
	{
		compute_color_error_for_every_integer_count_and_quantization_level(
		    encode_hdr_rgb, encode_hdr_alpha, i,
		    pt, &(eci[i]), ep, error_weightings, best_error[i],
		    format_of_choice[i]);
	}

	alignas(ASTCENC_VECALIGN) float errors_of_best_combination[MAX_WEIGHT_MODES];
	alignas(ASTCENC_VECALIGN) int best_quantization_levels[MAX_WEIGHT_MODES];
	int best_quantization_levels_mod[MAX_WEIGHT_MODES];
	int best_ep_formats[MAX_WEIGHT_MODES][4];

#if ASTCENC_SIMD_WIDTH > 1
	// have to ensure that the "overstep" of the last iteration in the vectorized
	// loop will contain data that will never be picked as best candidate
	const int packed_mode_count = bsd->block_mode_packed_count;
	const int packed_mode_count_simd_up = (packed_mode_count + ASTCENC_SIMD_WIDTH - 1) / ASTCENC_SIMD_WIDTH * ASTCENC_SIMD_WIDTH;
	for (int i = packed_mode_count; i < packed_mode_count_simd_up; ++i)
	{
		errors_of_best_combination[i] = 1e30f;
		best_quantization_levels[i] = 0;
		best_quantization_levels_mod[i] = 0;
	}
#endif // #if ASTCENC_SIMD_WIDTH > 1

	// code for the case where the block contains 1 partition
	if (partition_count == 1)
	{
		int best_quantization_level;
		int best_format;
		float error_of_best_combination;
		for (int i = 0, ni = bsd->block_mode_packed_count; i < ni; ++i)
		{
			if (qwt_errors[i] >= 1e29f)
			{
				errors_of_best_combination[i] = 1e30f;
				continue;
			}

			one_partition_find_best_combination_for_bitcount(
			    best_error[0], format_of_choice[0], qwt_bitcounts[i],
			    &best_quantization_level, &best_format, &error_of_best_combination);
			error_of_best_combination += qwt_errors[i];

			errors_of_best_combination[i] = error_of_best_combination;
			best_quantization_levels[i] = best_quantization_level;
			best_quantization_levels_mod[i] = best_quantization_level;
			best_ep_formats[i][0] = best_format;
		}
	}
	// code for the case where the block contains 2 partitions
	else if (partition_count == 2)
	{
		int best_quantization_level;
		int best_quantization_level_mod;
		int best_formats[2];
		float error_of_best_combination;

		float combined_best_error[21][7];
		int formats_of_choice[21][7][2];

		two_partitions_find_best_combination_for_every_quantization_and_integer_count(
		    best_error, format_of_choice, combined_best_error, formats_of_choice);


		for (int i = 0, ni = bsd->block_mode_packed_count; i < ni; ++i)
		{
			if (qwt_errors[i] >= 1e29f)
			{
				errors_of_best_combination[i] = 1e30f;
				continue;
			}

			two_partitions_find_best_combination_for_bitcount(
			    combined_best_error, formats_of_choice, qwt_bitcounts[i],
			    &best_quantization_level, &best_quantization_level_mod,
			    best_formats, &error_of_best_combination);

			error_of_best_combination += qwt_errors[i];

			errors_of_best_combination[i] = error_of_best_combination;
			best_quantization_levels[i] = best_quantization_level;
			best_quantization_levels_mod[i] = best_quantization_level_mod;
			best_ep_formats[i][0] = best_formats[0];
			best_ep_formats[i][1] = best_formats[1];
		}
	}
	// code for the case where the block contains 3 partitions
	else if (partition_count == 3)
	{
		int best_quantization_level;
		int best_quantization_level_mod;
		int best_formats[3];
		float error_of_best_combination;

		float combined_best_error[21][10];
		int formats_of_choice[21][10][3];

		three_partitions_find_best_combination_for_every_quantization_and_integer_count(
		    best_error, format_of_choice, combined_best_error, formats_of_choice);

		for (int i = 0, ni = bsd->block_mode_packed_count; i < ni; ++i)
		{
			if (qwt_errors[i] >= 1e29f)
			{
				errors_of_best_combination[i] = 1e30f;
				continue;
			}

			three_partitions_find_best_combination_for_bitcount(
			    combined_best_error, formats_of_choice, qwt_bitcounts[i],
			    &best_quantization_level, &best_quantization_level_mod,
			    best_formats, &error_of_best_combination);

			error_of_best_combination += qwt_errors[i];

			errors_of_best_combination[i] = error_of_best_combination;
			best_quantization_levels[i] = best_quantization_level;
			best_quantization_levels_mod[i] = best_quantization_level_mod;
			best_ep_formats[i][0] = best_formats[0];
			best_ep_formats[i][1] = best_formats[1];
			best_ep_formats[i][2] = best_formats[2];
		}
	}
	// code for the case where the block contains 4 partitions
	else if (partition_count == 4)
	{
		int best_quantization_level;
		int best_quantization_level_mod;
		int best_formats[4];
		float error_of_best_combination;

		float combined_best_error[21][13];
		int formats_of_choice[21][13][4];

		four_partitions_find_best_combination_for_every_quantization_and_integer_count(
		    best_error, format_of_choice, combined_best_error, formats_of_choice);

		for (int i = 0, ni = bsd->block_mode_packed_count; i < ni; ++i)
		{
			if (qwt_errors[i] >= 1e29f)
			{
				errors_of_best_combination[i] = 1e30f;
				continue;
			}

			four_partitions_find_best_combination_for_bitcount(
			    combined_best_error, formats_of_choice, qwt_bitcounts[i],
			    &best_quantization_level, &best_quantization_level_mod,
			    best_formats, &error_of_best_combination);

			error_of_best_combination += qwt_errors[i];

			errors_of_best_combination[i] = error_of_best_combination;
			best_quantization_levels[i] = best_quantization_level;
			best_quantization_levels_mod[i] = best_quantization_level_mod;
			best_ep_formats[i][0] = best_formats[0];
			best_ep_formats[i][1] = best_formats[1];
			best_ep_formats[i][2] = best_formats[2];
			best_ep_formats[i][3] = best_formats[3];
		}
	}

	// finally, go through the results and pick the best-looking modes.
	int best_error_weights[TUNE_MAX_TRIAL_CANDIDATES];
	for (int i = 0; i < tune_candidate_limit; i++)
	{
#if 0
		// reference; scalar code
		float best_ep_error = 1e30f;
		int best_error_index = -1;
		for (int j = 0, npack = bsd->block_mode_packed_count; j < npack; ++j)
		{
			if (errors_of_best_combination[j] < best_ep_error && best_quantization_levels[j] >= 5)
			{
				best_ep_error = errors_of_best_combination[j];
				best_error_index = j;
			}
		}
#else
		// find best mode, SIMD N-wide way
		static_assert((MAX_WEIGHT_MODES % ASTCENC_SIMD_WIDTH) == 0, "MAX_WEIGHT_MODES should be multiple of ASTCENC_SIMD_WIDTH");
		vint vbest_error_index(-1);
		vfloat vbest_ep_error(1e30f);
		vint lane_ids = vint::lane_id();
		for (int j = 0, npack = bsd->block_mode_packed_count; j < npack; j += ASTCENC_SIMD_WIDTH)
		{
			vfloat err = vfloat(&errors_of_best_combination[j]);
			vmask mask1 = err < vbest_ep_error;
			vmask mask2 = vint(&best_quantization_levels[j]) > vint(4);
			vmask mask = mask1 & mask2;
			vbest_ep_error = select(vbest_ep_error, err, mask);
			vbest_error_index = select(vbest_error_index, lane_ids, mask);
			lane_ids = lane_ids + vint(ASTCENC_SIMD_WIDTH);
		}

		// pick final best mode from the SIMD result.
		// note that if multiple SIMD lanes have "best" score,
		// we want to pick one with the lowest index, i.e. what
		// would happen if code was purely scalar.
		vmask lanes_with_min_error = vbest_ep_error == hmin(vbest_ep_error);
		// take smallest index from the SIMD lanes that had the best score
		vbest_error_index = select(vint(0x7fffffff), vbest_error_index, lanes_with_min_error);
		vbest_error_index = hmin(vbest_error_index);
		int best_error_index = vbest_error_index.lane(0);
#endif

		best_error_weights[i] = best_error_index;

		if (best_error_index >= 0)
		{
			errors_of_best_combination[best_error_index] = 1e30f;
		}
	}

	for (int i = 0; i < tune_candidate_limit; i++)
	{
		quantized_weight[i] = best_error_weights[i];
		if (quantized_weight[i] >= 0)
		{
			quantization_level[i] = best_quantization_levels[best_error_weights[i]];
			assert(quantization_level[i] >= 0 && quantization_level[i] < 21);
			quantization_level_mod[i] = best_quantization_levels_mod[best_error_weights[i]];
			for (int j = 0; j < partition_count; j++)
			{
				partition_format_specifiers[i][j] = best_ep_formats[best_error_weights[i]][j];
			}
		}
	}
}

#endif
