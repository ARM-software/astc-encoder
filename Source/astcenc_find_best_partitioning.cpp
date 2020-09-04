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
 * @brief Functions for finding best partition for a block.
 *
 * Major step 1:
 * - find best partitioning assuming uncorrelated colors
 * - find best partitioning assuming RGBS color representation
 *
 * Finding best partitioning for a block:
 *
 * foreach available partitioning:
 * - compute mean-color-value and dominant direction.
 * - this defines two lines, both of which go through the mean-color-value.
 * - one line has a direction defined by the dominant direction; this is used
 *   to assess the error from using an uncorrelated color representation.
 * - the other line goes through (0,0,0,1) and is used to assess the error from
 *   using an RGBS color representation.
 * - we then compute, as a sum across the block, the squared-errors that result
 *   from using the dominant-direction-lines and the squared-errors that result
 *   from using the 0001-lines.
 *
 *	Partition table representation:
 *	We have 3 tables, each with 1024 partitions
 *	(these correspond to the 3x128 hardware partitions crossed with all the
 *	partition-transform modes in the hardware.)
 *
 *	For each partitioning, we have:
 *	* a 4-entry table indicating how many texels there are in each of the 4
 *	  partitions. this may be from 2 to about 60 or so.
 *	* a 64-entry table indicating the partition index of each of the 64 texels
 *	  in the block. each index may be 0, 1, 2 or 3.
 *
 * each element in the table is an uint8_t indicating partition index (0, 1, 2 or 3)
 */

#include "astcenc_internal.h"

static void compute_alpha_range(
	int texels_per_block,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	float alpha_range[4]
) {
	float alpha_min[4];
	float alpha_max[4];

	int partition_count = pt->partition_count;
	for (int i = 0; i < partition_count; i++)
	{
		alpha_min[i] = 1e38f;
		alpha_max[i] = -1e38f;
	}

	for (int i = 0; i < texels_per_block; i++)
	{
		if (ewb->texel_weight[i] > 1e-10f)
		{
			int partition = pt->partition_of_texel[i];
			float alphaval = blk->data_a[i];

			if (alphaval > alpha_max[partition])
			{
				alpha_max[partition] = alphaval;
			}

			if (alphaval < alpha_min[partition])
			{
				alpha_min[partition] = alphaval;
			}
		}
	}

	for (int i = 0; i < partition_count; i++)
	{
		alpha_range[i] = alpha_max[i] - alpha_min[i];
		if (alpha_range[i] <= 0.0f)
		{
			alpha_range[i] = 1e-10f;
		}
	}
}

static void compute_rgb_range(
	int texels_per_block,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	float3 rgb_range[4]
) {
	float3 rgb_min[4];
	float3 rgb_max[4];

	int partition_count = pt->partition_count;
	for (int i = 0; i < partition_count; i++)
	{
		rgb_min[i] = float3(1e38f, 1e38f, 1e38f);
		rgb_max[i] = float3(-1e38f, -1e38f, -1e38f);
	}

	for (int i = 0; i < texels_per_block; i++)
	{
		if (ewb->texel_weight[i] > 1e-10f)
		{
			int partition = pt->partition_of_texel[i];

			float redval = blk->data_r[i];
			if (redval > rgb_max[partition].x)
			{
				rgb_max[partition].x = redval;
			}

			if (redval < rgb_min[partition].x)
			{
				rgb_min[partition].x = redval;
			}

			float greenval = blk->data_g[i];
			if (greenval > rgb_max[partition].y)
			{
				rgb_max[partition].y = greenval;
			}

			if (greenval < rgb_min[partition].y)
			{
				rgb_min[partition].y = greenval;
			}

			float blueval = blk->data_b[i];
			if (blueval > rgb_max[partition].z)
			{
				rgb_max[partition].z = blueval;
			}

			if (blueval < rgb_min[partition].z)
			{
				rgb_min[partition].z = blueval;
			}
		}
	}

	// Covert min/max into ranges forcing a min range of 1e-10
	// to avoid divide by zeros later ...
	for (int i = 0; i < partition_count; i++)
	{
		rgb_range[i].x = rgb_max[i].x - rgb_min[i].x;
		if (rgb_range[i].x <= 0.0f)
		{
			rgb_range[i].x = 1e-10f;
		}

		rgb_range[i].y = rgb_max[i].y - rgb_min[i].y;
		if (rgb_range[i].y <= 0.0f)
		{
			rgb_range[i].y = 1e-10f;
		}

		rgb_range[i].z = rgb_max[i].z - rgb_min[i].z;
		if (rgb_range[i].z <= 0.0f)
		{
			rgb_range[i].z = 1e-10f;
		}
	}
}

void compute_partition_error_color_weightings(
	const block_size_descriptor* bsd,
	const error_weight_block* ewb,
	const partition_info* pi,
	float4 error_weightings[4],
	float4 color_scalefactors[4]
) {
	int texels_per_block = bsd->texel_count;
	int pcnt = pi->partition_count;

	for (int i = 0; i < pcnt; i++)
	{
		error_weightings[i] = float4(1e-12f, 1e-12f, 1e-12f, 1e-12f);
	}

	for (int i = 0; i < texels_per_block; i++)
	{
		int part = pi->partition_of_texel[i];
		error_weightings[part] = error_weightings[part] + ewb->error_weights[i];
	}

	for (int i = 0; i < pcnt; i++)
	{
		error_weightings[i] = error_weightings[i] * (1.0f / pi->texels_per_partition[i]);
	}

	for (int i = 0; i < pcnt; i++)
	{
		color_scalefactors[i].x = astc::sqrt(error_weightings[i].x);
		color_scalefactors[i].y = astc::sqrt(error_weightings[i].y);
		color_scalefactors[i].z = astc::sqrt(error_weightings[i].z);
		color_scalefactors[i].w = astc::sqrt(error_weightings[i].w);
	}
}

/* main function to identify the best partitioning for a given number of texels */
void find_best_partitionings(
	int partition_search_limit,
	const block_size_descriptor* bsd,
	int partition_count,
	const imageblock* pb,
	const error_weight_block* ewb,
	int candidates_to_return,
	int* best_partitions_uncorrelated,
	int* best_partitions_samechroma,
	int* best_partitions_dual_weight_planes
) {
	// constant used to estimate quantization error for a given partitioning;
	// the optimal value for this constant depends on bitrate.
	// These constants have been determined empirically.
	int texels_per_block = bsd->texel_count;
	float weight_imprecision_estim = 100.0f;
	if (texels_per_block <= 20)
	{
		weight_imprecision_estim = 0.03f;
	}
	else if (texels_per_block <= 31)
	{
		weight_imprecision_estim = 0.04f;
	}
	else if (texels_per_block <= 41)
	{
		weight_imprecision_estim = 0.05f;
	}
	else
	{
		weight_imprecision_estim = 0.055f;
	}

	int partition_sequence[PARTITION_COUNT];

	kmeans_compute_partition_ordering(bsd, partition_count, pb, partition_sequence);

	float weight_imprecision_estim_squared = weight_imprecision_estim * weight_imprecision_estim;

	int uses_alpha = imageblock_uses_alpha(pb);

	const partition_info* ptab = get_partition_table(bsd, partition_count);

	// partitioning errors assuming uncorrelated-chrominance endpoints
	float uncorr_errors[PARTITION_COUNT];
	// partitioning errors assuming same-chrominance endpoints
	float samechroma_errors[PARTITION_COUNT];

	// partitioning errors assuming that one of the color channels
	// is uncorrelated from all the other ones
	float4 separate_errors[PARTITION_COUNT];

	int defacto_search_limit = PARTITION_COUNT - 1;

	if (uses_alpha)
	{

		for (int i = 0; i < PARTITION_COUNT; i++)
		{
			int partition = partition_sequence[i];
			int bk_partition_count = ptab[partition].partition_count;

			if (bk_partition_count < partition_count)
			{
				uncorr_errors[i] = 1e35f;
				samechroma_errors[i] = 1e35f;
				separate_errors[i] = float4(1e35f, 1e35f, 1e35f, 1e35f);
				continue;
			}

			// the sentinel value for partitions above the search limit must be smaller
			// than the sentinel value for invalid partitions
			if (i >= partition_search_limit)
			{
				defacto_search_limit = i;
				uncorr_errors[i] = 1e34f;
				samechroma_errors[i] = 1e34f;
				separate_errors[i] = float4(1e34f, 1e34f, 1e34f, 1e34f);
				break;
			}

			// compute the weighting to give to each color channel
			// in each partition.
			float4 error_weightings[4];
			float4 color_scalefactors[4];
			float4 inverse_color_scalefactors[4];
			compute_partition_error_color_weightings(bsd, ewb, ptab + partition, error_weightings, color_scalefactors);

			for (int j = 0; j < partition_count; j++)
			{
				inverse_color_scalefactors[j].x = 1.0f / MAX(color_scalefactors[j].x, 1e-7f);
				inverse_color_scalefactors[j].y = 1.0f / MAX(color_scalefactors[j].y, 1e-7f);
				inverse_color_scalefactors[j].z = 1.0f / MAX(color_scalefactors[j].z, 1e-7f);
				inverse_color_scalefactors[j].w = 1.0f / MAX(color_scalefactors[j].w, 1e-7f);
			}

			float4 averages[4];
			float4 directions_rgba[4];

			compute_averages_and_directions_rgba(ptab + partition, pb, ewb,
			                                     color_scalefactors, averages,
			                                     directions_rgba);

			line4 uncorr_lines[4];
			line4 samechroma_lines[4];
			line3 separate_red_lines[4];
			line3 separate_green_lines[4];
			line3 separate_blue_lines[4];
			line3 separate_alpha_lines[4];

			processed_line4 proc_uncorr_lines[4];
			processed_line4 proc_samechroma_lines[4];
			processed_line3 proc_separate_red_lines[4];
			processed_line3 proc_separate_green_lines[4];
			processed_line3 proc_separate_blue_lines[4];
			processed_line3 proc_separate_alpha_lines[4];

			float uncorr_linelengths[4];
			float samechroma_linelengths[4];
			float4 separate_linelengths[4];

			for (int j = 0; j < partition_count; j++)
			{
				uncorr_lines[j].a = averages[j];
				if (dot(directions_rgba[j], directions_rgba[j]) == 0.0f)
				{
					uncorr_lines[j].b = normalize(float4(1.0f, 1.0f, 1.0f, 1.0f));
				}
				else
				{
					uncorr_lines[j].b = normalize(directions_rgba[j]);
				}

				proc_uncorr_lines[j].amod = (uncorr_lines[j].a - uncorr_lines[j].b * dot(uncorr_lines[j].a, uncorr_lines[j].b)) * inverse_color_scalefactors[j];
				proc_uncorr_lines[j].bs = (uncorr_lines[j].b * color_scalefactors[j]);
				proc_uncorr_lines[j].bis = (uncorr_lines[j].b * inverse_color_scalefactors[j]);

				samechroma_lines[j].a = float4(0.0f, 0.0f, 0.0f, 0.0f);
				if (dot(averages[j], averages[j]) == 0.0f)
				{
					samechroma_lines[j].b = normalize(float4(1.0f, 1.0f, 1.0f, 1.0f));
				}
				else
				{
					samechroma_lines[j].b = normalize(averages[j]);
				}

				proc_samechroma_lines[j].amod = (samechroma_lines[j].a - samechroma_lines[j].b * dot(samechroma_lines[j].a, samechroma_lines[j].b)) * inverse_color_scalefactors[j];
				proc_samechroma_lines[j].bs = (samechroma_lines[j].b * color_scalefactors[j]);
				proc_samechroma_lines[j].bis = (samechroma_lines[j].b * inverse_color_scalefactors[j]);

				separate_red_lines[j].a = float3(averages[j].y, averages[j].z, averages[j].w);
				float3 dirs_gba = float3(directions_rgba[j].y, directions_rgba[j].z, directions_rgba[j].w);
				if (dot(dirs_gba, dirs_gba) == 0.0f)
				{
					separate_red_lines[j].b = normalize(float3(1.0f, 1.0f, 1.0f));
				}
				else
				{
					separate_red_lines[j].b = normalize(dirs_gba);
				}

				separate_green_lines[j].a = float3(averages[j].x, averages[j].z, averages[j].w);
				float3 dirs_rba = float3(directions_rgba[j].x, directions_rgba[j].z, directions_rgba[j].w);
				if (dot(dirs_rba, dirs_rba) == 0.0f)
				{
					separate_green_lines[j].b = normalize(float3(1.0f, 1.0f, 1.0f));
				}
				else
				{
					separate_green_lines[j].b = normalize(dirs_rba);
				}

				separate_blue_lines[j].a = float3(averages[j].x, averages[j].y, averages[j].w);
				float3 dirs_rga = float3(directions_rgba[j].x, directions_rgba[j].y, directions_rgba[j].w);
				if (dot(dirs_rga, dirs_rga) == 0.0f)
				{
					separate_blue_lines[j].b = normalize(float3(1.0f, 1.0f, 1.0f));
				}
				else
				{
					separate_blue_lines[j].b = normalize(dirs_rga);
				}

				separate_alpha_lines[j].a = float3(averages[j].x, averages[j].y, averages[j].z);
				float3 dirs_rgb = float3(directions_rgba[j].x, directions_rgba[j].y, directions_rgba[j].z);
				if (dot(dirs_rgb, dirs_rgb) == 0.0f)
				{
					separate_alpha_lines[j].b = normalize(float3(1.0f, 1.0f, 1.0f));
				}
				else
				{
					separate_alpha_lines[j].b = normalize(dirs_rgb);
				}

				proc_separate_red_lines[j].amod = (separate_red_lines[j].a - separate_red_lines[j].b * dot(separate_red_lines[j].a, separate_red_lines[j].b)) * float3(inverse_color_scalefactors[j].y, inverse_color_scalefactors[j].z, inverse_color_scalefactors[j].w);
				proc_separate_red_lines[j].bs = (separate_red_lines[j].b * float3(color_scalefactors[j].y, color_scalefactors[j].z, color_scalefactors[j].w));
				proc_separate_red_lines[j].bis = (separate_red_lines[j].b * float3(inverse_color_scalefactors[j].y, inverse_color_scalefactors[j].z, inverse_color_scalefactors[j].w));

				proc_separate_green_lines[j].amod =
					(separate_green_lines[j].a - separate_green_lines[j].b * dot(separate_green_lines[j].a, separate_green_lines[j].b)) * float3(inverse_color_scalefactors[j].x, inverse_color_scalefactors[j].z, inverse_color_scalefactors[j].w);
				proc_separate_green_lines[j].bs = (separate_green_lines[j].b * float3(color_scalefactors[j].x, color_scalefactors[j].z, color_scalefactors[j].w));
				proc_separate_green_lines[j].bis = (separate_green_lines[j].b * float3(inverse_color_scalefactors[j].x, inverse_color_scalefactors[j].z, inverse_color_scalefactors[j].w));

				proc_separate_blue_lines[j].amod = (separate_blue_lines[j].a - separate_blue_lines[j].b * dot(separate_blue_lines[j].a, separate_blue_lines[j].b)) * float3(inverse_color_scalefactors[j].x, inverse_color_scalefactors[j].y, inverse_color_scalefactors[j].w);
				proc_separate_blue_lines[j].bs = (separate_blue_lines[j].b * float3(color_scalefactors[j].x, color_scalefactors[j].y, color_scalefactors[j].w));
				proc_separate_blue_lines[j].bis = (separate_blue_lines[j].b * float3(inverse_color_scalefactors[j].x, inverse_color_scalefactors[j].y, inverse_color_scalefactors[j].w));

				proc_separate_alpha_lines[j].amod =
					(separate_alpha_lines[j].a - separate_alpha_lines[j].b * dot(separate_alpha_lines[j].a, separate_alpha_lines[j].b)) * float3(inverse_color_scalefactors[j].x, inverse_color_scalefactors[j].y, inverse_color_scalefactors[j].z);
				proc_separate_alpha_lines[j].bs = (separate_alpha_lines[j].b * float3(color_scalefactors[j].x, color_scalefactors[j].y, color_scalefactors[j].z));
				proc_separate_alpha_lines[j].bis = (separate_alpha_lines[j].b * float3(inverse_color_scalefactors[j].x, inverse_color_scalefactors[j].y, inverse_color_scalefactors[j].z));
			}

			float uncorr_error = 0.0f;
			float samechroma_error = 0.0f;
			float4 separate_error = float4(0.0f, 0.0f, 0.0f, 0.0f);
			compute_error_squared_rgba(ptab + partition,
			                           pb,
			                           ewb,
			                           proc_uncorr_lines,
			                           proc_samechroma_lines,
			                           proc_separate_red_lines,
			                           proc_separate_green_lines,
			                           proc_separate_blue_lines,
			                           proc_separate_alpha_lines,
			                           uncorr_linelengths,
			                           samechroma_linelengths,
			                           separate_linelengths,
			                           &uncorr_error,
			                           &samechroma_error,
			                           &separate_error);

			// compute minimum & maximum alpha values in each partition
			float3 rgb_range[4];
			float alpha_range[4];

			compute_alpha_range(bsd->texel_count, ptab + partition, pb, ewb, alpha_range);
			compute_rgb_range(bsd->texel_count, ptab + partition, pb, ewb, rgb_range);

			/*
			   Compute an estimate of error introduced by weight quantization imprecision.
			   This error is computed as follows, for each partition
			   1: compute the principal-axis vector (full length) in error-space
			   2: convert the principal-axis vector to regular RGB-space
			   3: scale the vector by a constant that estimates average quantization error
			   4: for each texel, square the vector, then do a dot-product with the texel's error weight;
			      sum up the results across all texels.
			   4(optimized): square the vector once, then do a dot-product with the average texel error,
			      then multiply by the number of texels.
			 */

			for (int j = 0; j < partition_count; j++)
			{
				float tpp = (float)(ptab[partition].texels_per_partition[j]);

				float4 ics = inverse_color_scalefactors[j];
				float4 error_weights = error_weightings[j] * (tpp * weight_imprecision_estim_squared);

				float4 uncorr_vector = (uncorr_lines[j].b * uncorr_linelengths[j]) * ics;
				float4 samechroma_vector = (samechroma_lines[j].b * samechroma_linelengths[j]) * ics;
				float3 separate_red_vector = (separate_red_lines[j].b * separate_linelengths[j].x) * float3(ics.y, ics.z, ics.w);
				float3 separate_green_vector = (separate_green_lines[j].b * separate_linelengths[j].y) * float3(ics.x, ics.z, ics.w);
				float3 separate_blue_vector = (separate_blue_lines[j].b * separate_linelengths[j].z) * float3(ics.x, ics.y, ics.w);
				float3 separate_alpha_vector = (separate_alpha_lines[j].b * separate_linelengths[j].w) * float3(ics.x, ics.y, ics.z);

				uncorr_vector = uncorr_vector * uncorr_vector;
				samechroma_vector = samechroma_vector * samechroma_vector;
				separate_red_vector = separate_red_vector * separate_red_vector;
				separate_green_vector = separate_green_vector * separate_green_vector;
				separate_blue_vector = separate_blue_vector * separate_blue_vector;
				separate_alpha_vector = separate_alpha_vector * separate_alpha_vector;

				uncorr_error += dot(uncorr_vector, error_weights);
				samechroma_error += dot(samechroma_vector, error_weights);
				separate_error.x += dot(separate_red_vector, float3(error_weights.y, error_weights.z, error_weights.w));
				separate_error.y += dot(separate_green_vector, float3(error_weights.x, error_weights.z, error_weights.w));
				separate_error.z += dot(separate_blue_vector, float3(error_weights.x, error_weights.y, error_weights.w));
				separate_error.w += dot(separate_alpha_vector, float3(error_weights.x, error_weights.y, error_weights.z));

				separate_error.x += rgb_range[j].x * rgb_range[j].x * error_weights.x;
				separate_error.y += rgb_range[j].y * rgb_range[j].y * error_weights.y;
				separate_error.z += rgb_range[j].z * rgb_range[j].z * error_weights.z;
				separate_error.w += alpha_range[j] * alpha_range[j] * error_weights.w;
			}

			uncorr_errors[i] = uncorr_error;
			samechroma_errors[i] = samechroma_error;
			separate_errors[i] = separate_error;
		}
	}
	else
	{
		for (int i = 0; i < PARTITION_COUNT; i++)
		{
			int partition = partition_sequence[i];

			int bk_partition_count = ptab[partition].partition_count;
			if (bk_partition_count < partition_count)
			{
				uncorr_errors[i] = 1e35f;
				samechroma_errors[i] = 1e35f;
				separate_errors[i] = float4(1e35f, 1e35f, 1e35f, 1e35f);
				continue;
			}

			// the sentinel value for valid partitions above the search limit must be smaller
			// than the sentinel value for invalid partitions
			if (i >= partition_search_limit)
			{
				defacto_search_limit = i;
				uncorr_errors[i] = 1e34f;
				samechroma_errors[i] = 1e34f;
				separate_errors[i] = float4(1e34f, 1e34f, 1e34f, 1e34f);
				break;
			}

			// compute the weighting to give to each color channel
			// in each partition.
			float4 error_weightings[4];
			float4 color_scalefactors[4];
			float4 inverse_color_scalefactors[4];

			compute_partition_error_color_weightings(bsd, ewb, ptab + partition, error_weightings, color_scalefactors);

			for (int j = 0; j < partition_count; j++)
			{
				inverse_color_scalefactors[j].x = 1.0f / MAX(color_scalefactors[j].x, 1e-7f);
				inverse_color_scalefactors[j].y = 1.0f / MAX(color_scalefactors[j].y, 1e-7f);
				inverse_color_scalefactors[j].z = 1.0f / MAX(color_scalefactors[j].z, 1e-7f);
				inverse_color_scalefactors[j].w = 1.0f / MAX(color_scalefactors[j].w, 1e-7f);
			}

			float3 averages[4];
			float3 directions_rgb[4];

			compute_averages_and_directions_rgb(ptab + partition, pb, ewb, color_scalefactors, averages, directions_rgb);

			line3 uncorr_lines[4];
			line3 samechroma_lines[4];
			line2 separate_red_lines[4];
			line2 separate_green_lines[4];
			line2 separate_blue_lines[4];

			processed_line3 proc_uncorr_lines[4];
			processed_line3 proc_samechroma_lines[4];

			processed_line2 proc_separate_red_lines[4];
			processed_line2 proc_separate_green_lines[4];
			processed_line2 proc_separate_blue_lines[4];

			float uncorr_linelengths[4];
			float samechroma_linelengths[4];
			float3 separate_linelengths[4];

			for (int j = 0; j < partition_count; j++)
			{
				uncorr_lines[j].a = averages[j];
				if (dot(directions_rgb[j], directions_rgb[j]) == 0.0f)
				{
					uncorr_lines[j].b = normalize(float3(1.0f, 1.0f, 1.0f));
				}
				else
				{
					uncorr_lines[j].b = normalize(directions_rgb[j]);
				}

				samechroma_lines[j].a = float3(0.0f, 0.0f, 0.0f);
				if (dot(averages[j], averages[j]) == 0.0f)
				{
					samechroma_lines[j].b = normalize(float3(1.0f, 1.0f, 1.0f));
				}
				else
				{
					samechroma_lines[j].b = normalize(averages[j]);
				}

				proc_uncorr_lines[j].amod = (uncorr_lines[j].a - uncorr_lines[j].b * dot(uncorr_lines[j].a, uncorr_lines[j].b)) * float3(inverse_color_scalefactors[j].x, inverse_color_scalefactors[j].y, inverse_color_scalefactors[j].z);
				proc_uncorr_lines[j].bs = (uncorr_lines[j].b * float3(color_scalefactors[j].x, color_scalefactors[j].y, color_scalefactors[j].z));
				proc_uncorr_lines[j].bis = (uncorr_lines[j].b * float3(inverse_color_scalefactors[j].x, inverse_color_scalefactors[j].y, inverse_color_scalefactors[j].z));

				proc_samechroma_lines[j].amod = (samechroma_lines[j].a - samechroma_lines[j].b * dot(samechroma_lines[j].a, samechroma_lines[j].b)) * float3(inverse_color_scalefactors[j].x, inverse_color_scalefactors[j].y, inverse_color_scalefactors[j].z);
				proc_samechroma_lines[j].bs = (samechroma_lines[j].b * float3(color_scalefactors[j].x, color_scalefactors[j].y, color_scalefactors[j].z));
				proc_samechroma_lines[j].bis = (samechroma_lines[j].b * float3(inverse_color_scalefactors[j].x, inverse_color_scalefactors[j].y, inverse_color_scalefactors[j].z));

				separate_red_lines[j].a = float2(averages[j].y, averages[j].z);
				float2 dirs_gb = float2(directions_rgb[j].y, directions_rgb[j].z);
				if (dot(dirs_gb, dirs_gb) == 0.0f)
				{
					separate_red_lines[j].b = normalize(float2(1.0f, 1.0f));
				}
				else
				{
					separate_red_lines[j].b = normalize(dirs_gb);
				}

				separate_green_lines[j].a = float2(averages[j].x, averages[j].z);
				float2 dirs_rb = float2(directions_rgb[j].x, directions_rgb[j].z);
				if (dot(dirs_rb, dirs_rb) == 0.0f)
				{
					separate_green_lines[j].b = normalize(float2(1.0f, 1.0f));
				}
				else
				{
					separate_green_lines[j].b = normalize(dirs_rb);
				}

				separate_blue_lines[j].a = float2(averages[j].x, averages[j].y);
				float2 dirs_rg = float2(directions_rgb[j].x, directions_rgb[j].y);
				if (dot(dirs_rg, dirs_rg) == 0.0f)
				{
					separate_blue_lines[j].b = normalize(float2(1.0f, 1.0f));
				}
				else
				{
					separate_blue_lines[j].b = normalize(dirs_rg);
				}

				proc_separate_red_lines[j].amod = (separate_red_lines[j].a - separate_red_lines[j].b * dot(separate_red_lines[j].a, separate_red_lines[j].b)) * float2(inverse_color_scalefactors[j].y, inverse_color_scalefactors[j].z);
				proc_separate_red_lines[j].bs = (separate_red_lines[j].b * float2(color_scalefactors[j].y, color_scalefactors[j].z));
				proc_separate_red_lines[j].bis = (separate_red_lines[j].b * float2(inverse_color_scalefactors[j].y, inverse_color_scalefactors[j].z));

				proc_separate_green_lines[j].amod = (separate_green_lines[j].a - separate_green_lines[j].b * dot(separate_green_lines[j].a, separate_green_lines[j].b)) * float2(inverse_color_scalefactors[j].x, inverse_color_scalefactors[j].z);
				proc_separate_green_lines[j].bs = (separate_green_lines[j].b * float2(color_scalefactors[j].x, color_scalefactors[j].z));
				proc_separate_green_lines[j].bis = (separate_green_lines[j].b * float2(inverse_color_scalefactors[j].x, inverse_color_scalefactors[j].z));

				proc_separate_blue_lines[j].amod = (separate_blue_lines[j].a - separate_blue_lines[j].b * dot(separate_blue_lines[j].a, separate_blue_lines[j].b)) * float2(inverse_color_scalefactors[j].x, inverse_color_scalefactors[j].y);
				proc_separate_blue_lines[j].bs = (separate_blue_lines[j].b * float2(color_scalefactors[j].x, color_scalefactors[j].y));
				proc_separate_blue_lines[j].bis = (separate_blue_lines[j].b * float2(inverse_color_scalefactors[j].x, inverse_color_scalefactors[j].y));
			}

			float uncorr_error = 0.0f;
			float samechroma_error = 0.0f;
			float3 separate_error = float3(0.0f, 0.0f, 0.0f);

			compute_error_squared_rgb(ptab + partition,
			                          pb,
			                          ewb,
			                          proc_uncorr_lines,
			                          proc_samechroma_lines,
			                          proc_separate_red_lines,
			                          proc_separate_green_lines,
			                          proc_separate_blue_lines,
			                          uncorr_linelengths,
			                          samechroma_linelengths,
			                          separate_linelengths,
			                          &uncorr_error,
			                          &samechroma_error,
			                          &separate_error);

			float3 rgb_range[4];
			compute_rgb_range(bsd->texel_count, ptab + partition, pb, ewb, rgb_range);

			/*
			   compute an estimate of error introduced by weight imprecision.
			   This error is computed as follows, for each partition
			   1: compute the principal-axis vector (full length) in error-space
			   2: convert the principal-axis vector to regular RGB-space
			   3: scale the vector by a constant that estimates average quantization error.
			   4: for each texel, square the vector, then do a dot-product with the texel's error weight;
			      sum up the results across all texels.
			   4(optimized): square the vector once, then do a dot-product with the average texel error,
			      then multiply by the number of texels.
			 */

			for (int j = 0; j < partition_count; j++)
			{
				float tpp = (float)(ptab[partition].texels_per_partition[j]);

				float3 ics = float3(inverse_color_scalefactors[j].x, inverse_color_scalefactors[j].y, inverse_color_scalefactors[j].z);
				float3 error_weights = float3(error_weightings[j].x, error_weightings[j].y, error_weightings[j].z) * (tpp * weight_imprecision_estim_squared);

				float3 uncorr_vector = (uncorr_lines[j].b * uncorr_linelengths[j]) * ics;
				float3 samechroma_vector = (samechroma_lines[j].b * samechroma_linelengths[j]) * ics;

				float2 separate_red_vector = (separate_red_lines[j].b * separate_linelengths[j].x) * float2(ics.y, ics.z);
				float2 separate_green_vector = (separate_green_lines[j].b * separate_linelengths[j].y) * float2(ics.x, ics.z);
				float2 separate_blue_vector = (separate_blue_lines[j].b * separate_linelengths[j].z) * float2(ics.x, ics.y);

				uncorr_vector = uncorr_vector * uncorr_vector;
				samechroma_vector = samechroma_vector * samechroma_vector;
				separate_red_vector = separate_red_vector * separate_red_vector;
				separate_green_vector = separate_green_vector * separate_green_vector;
				separate_blue_vector = separate_blue_vector * separate_blue_vector;

				uncorr_error += dot(uncorr_vector, error_weights);
				samechroma_error += dot(samechroma_vector, error_weights);
				separate_error.x += dot(separate_red_vector, float2(error_weights.y, error_weights.z));
				separate_error.y += dot(separate_green_vector, float2(error_weights.x, error_weights.z));
				separate_error.z += dot(separate_blue_vector, float2(error_weights.x, error_weights.y));

				separate_error.x += rgb_range[j].x * rgb_range[j].x * error_weights.x;
				separate_error.y += rgb_range[j].y * rgb_range[j].y * error_weights.y;
				separate_error.z += rgb_range[j].z * rgb_range[j].z * error_weights.z;
			}

			uncorr_errors[i] = uncorr_error;
			samechroma_errors[i] = samechroma_error;
			separate_errors[i] = float4(separate_error.x, separate_error.y, separate_error.z, 0.0f);
		}
	}

	for (int i = 0; i < candidates_to_return; i++)
	{
		int best_uncorr_partition = 0;
		int best_samechroma_partition = 0;
		float best_uncorr_error = 1e30f;
		float best_samechroma_error = 1e30f;

		for (int j = 0; j <= defacto_search_limit; j++)
		{
			if (uncorr_errors[j] < best_uncorr_error)
			{
				best_uncorr_partition = j;
				best_uncorr_error = uncorr_errors[j];
			}
		}

		best_partitions_uncorrelated[i] = partition_sequence[best_uncorr_partition];
		uncorr_errors[best_uncorr_partition] = 1e30f;
		samechroma_errors[best_uncorr_partition] = 1e30f;

		for (int j = 0; j <= defacto_search_limit; j++)
		{
			if (samechroma_errors[j] < best_samechroma_error)
			{
				best_samechroma_partition = j;
				best_samechroma_error = samechroma_errors[j];
			}
		}

		best_partitions_samechroma[i] = partition_sequence[best_samechroma_partition];
		samechroma_errors[best_samechroma_partition] = 1e30f;
		uncorr_errors[best_samechroma_partition] = 1e30f;
	}

	for (int i = 0; i < 2 * candidates_to_return; i++)
	{
		int best_partition = 0;
		int best_component = 0;
		float best_partition_error = 1e30f;

		for (int j = 0; j <= defacto_search_limit; j++)
		{
			if (separate_errors[j].x < best_partition_error)
			{
				best_component = 0;
				best_partition = j;
				best_partition_error = separate_errors[j].x;
			}

			if (separate_errors[j].y < best_partition_error)
			{
				best_component = 1;
				best_partition = j;
				best_partition_error = separate_errors[j].y;
			}

			if (separate_errors[j].z < best_partition_error)
			{
				best_component = 2;
				best_partition = j;
				best_partition_error = separate_errors[j].z;
			}

			if (uses_alpha)
			{
				if (separate_errors[j].w < best_partition_error)
				{
					best_component = 3;
					best_partition = j;
					best_partition_error = separate_errors[j].w;
				}
			}
		}

		switch(best_component)
		{
		case 0:
			separate_errors[best_partition].x = 1e30f;
			break;
		case 1:
			separate_errors[best_partition].y = 1e30f;
			break;
		case 2:
			separate_errors[best_partition].z = 1e30f;
			break;
		case 3:
			separate_errors[best_partition].w = 1e30f;
			break;
		}

		best_partition = (best_component << PARTITION_BITS) | partition_sequence[best_partition & (PARTITION_COUNT - 1)];
		best_partitions_dual_weight_planes[i] = best_partition;
	}
}

#endif
