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

// TODO: This slows down when inlining, but might change if we shrink caller
static void compute_rgba_range(
	int texels_per_block,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	vfloat4 rgba_range[4]
) {
	vfloat4 rgba_min[4];
	vfloat4 rgba_max[4];

	int partition_count = pt->partition_count;

	promise(partition_count > 0);
	for (int i = 0; i < partition_count; i++)
	{
		rgba_min[i] = vfloat4(1e38f);
		rgba_max[i] = vfloat4(-1e38f);
	}

	promise(texels_per_block > 0);
	for (int i = 0; i < texels_per_block; i++)
	{
		if (ewb->texel_weight[i] > 1e-10f)
		{
			int partition = pt->partition_of_texel[i];
			vfloat4 data = blk->texel(i);

			rgba_min[partition] = min(data, rgba_min[partition]);
			rgba_max[partition] = max(data, rgba_max[partition]);
		}
	}

	// Covert min/max into ranges forcing a min range of 1e-10
	// to avoid divide by zeros later ...
	for (int i = 0; i < partition_count; i++)
	{
		rgba_range[i] = max(rgba_max[i] - rgba_min[i], 1e-10f);
	}
}

void compute_partition_error_color_weightings(
	const block_size_descriptor* bsd,
	const error_weight_block* ewb,
	const partition_info* pi,
	vfloat4 error_weightings[4],
	vfloat4 color_scalefactors[4]
) {
	int texels_per_block = bsd->texel_count;
	int pcnt = pi->partition_count;

	for (int i = 0; i < pcnt; i++)
	{
		error_weightings[i] = vfloat4(1e-12f);
	}

	for (int i = 0; i < texels_per_block; i++)
	{
		int part = pi->partition_of_texel[i];
		error_weightings[part] = error_weightings[part] + ewb->error_weights[i];
	}

	for (int i = 0; i < pcnt; i++)
	{
		error_weightings[i] = error_weightings[i] * (1.0f / pi->texels_per_partition[i]);
		color_scalefactors[i] = sqrt(error_weightings[i]);
	}
}

/* main function to identify the best partitioning for a given number of texels */
void find_best_partitionings(
	const block_size_descriptor* bsd,
	const imageblock* blk,
	const error_weight_block* ewb,
	int partition_count,
	int partition_search_limit,
	int* best_partition_uncorrelated,
	int* best_partition_samechroma,
	int* best_partition_dualplane
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

	kmeans_compute_partition_ordering(bsd, partition_count, blk, partition_sequence);

	float weight_imprecision_estim_squared = weight_imprecision_estim * weight_imprecision_estim;

	int uses_alpha = imageblock_uses_alpha(blk);

	const partition_info* ptab = get_partition_table(bsd, partition_count);

	// Partitioning errors assuming uncorrelated-chrominance endpoints
	float uncorr_best_error { ERROR_CALC_DEFAULT };
	int uncorr_best_partition { 0 };

	// Partitioning errors assuming same-chrominance endpoints
	// Store two so we can always return one different to uncorr
	float samechroma_best_errors[2] { ERROR_CALC_DEFAULT, ERROR_CALC_DEFAULT };
	int samechroma_best_partitions[2] { 0, 0 };

	// Partitioning errors assuming that one color component is uncorrelated
	float separate_best_error { ERROR_CALC_DEFAULT };
	int separate_best_partition { 0 };
	int separate_best_component { 0 };

	if (uses_alpha)
	{

		for (int i = 0; i < partition_search_limit; i++)
		{
			int partition = partition_sequence[i];

			int bk_partition_count = ptab[partition].partition_count;
			if (bk_partition_count < partition_count)
			{
				continue;
			}

			// compute the weighting to give to each color channel
			// in each partition.
			vfloat4 error_weightings[4];
			vfloat4 color_scalefactors[4];
			vfloat4 inverse_color_scalefactors[4];
			compute_partition_error_color_weightings(bsd, ewb, ptab + partition, error_weightings, color_scalefactors);

			for (int j = 0; j < partition_count; j++)
			{
				inverse_color_scalefactors[j] = 1.0f / max(color_scalefactors[j], 1e-7f);
			}

			vfloat4 averages[4];
			vfloat4 directions_rgba[4];

			compute_averages_and_directions_rgba(ptab + partition, blk, ewb,
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

			for (int j = 0; j < partition_count; j++)
			{
				uncorr_lines[j].a = averages[j];
				if (dot_s(directions_rgba[j], directions_rgba[j]) == 0.0f)
				{
					uncorr_lines[j].b = normalize(vfloat4(1.0f));
				}
				else
				{
					uncorr_lines[j].b = normalize(directions_rgba[j]);
				}

				proc_uncorr_lines[j].amod = (uncorr_lines[j].a - uncorr_lines[j].b * dot(uncorr_lines[j].a, uncorr_lines[j].b)) * inverse_color_scalefactors[j];
				proc_uncorr_lines[j].bs = uncorr_lines[j].b * color_scalefactors[j];
				proc_uncorr_lines[j].bis = uncorr_lines[j].b * inverse_color_scalefactors[j];

				samechroma_lines[j].a = vfloat4::zero();
				if (dot_s(averages[j], averages[j]) == 0.0f)
				{
					samechroma_lines[j].b = normalize(vfloat4(1.0f));
				}
				else
				{
					samechroma_lines[j].b = normalize(averages[j]);
				}

				proc_samechroma_lines[j].amod = (samechroma_lines[j].a - samechroma_lines[j].b * dot(samechroma_lines[j].a, samechroma_lines[j].b)) * inverse_color_scalefactors[j];
				proc_samechroma_lines[j].bs = samechroma_lines[j].b * color_scalefactors[j];
				proc_samechroma_lines[j].bis = samechroma_lines[j].b * inverse_color_scalefactors[j];

				// TODO: Future: implement swizzle 3 as a vfloat 4
				separate_red_lines[j].a = averages[j].swz<1, 2, 3, 0>();
				separate_red_lines[j].a.set_lane<3>(0.0f);
				vfloat4 dirs_gba = directions_rgba[j].swz<1, 2, 3, 0>();
				if (dot3_s(dirs_gba, dirs_gba) == 0.0f)
				{
					separate_red_lines[j].b = normalize(vfloat4(1.0f, 1.0f, 1.0f, 0.0f));
				}
				else
				{
					separate_red_lines[j].b = normalize(dirs_gba);
				}

				separate_green_lines[j].a = averages[j].swz<0, 2, 3, 1>();
				separate_red_lines[j].a.set_lane<3>(0.0f);
				vfloat4 dirs_rba = directions_rgba[j].swz<0, 2, 3, 1>();
				if (dot3_s(dirs_rba, dirs_rba) == 0.0f)
				{
					separate_green_lines[j].b = normalize(vfloat4(1.0f, 1.0f, 1.0f, 0.0f));
				}
				else
				{
					separate_green_lines[j].b = normalize(dirs_rba);
				}

				separate_blue_lines[j].a = averages[j].swz<0, 1, 3, 2>();
				separate_blue_lines[j].a.set_lane<3>(0.0f);
				vfloat4 dirs_rga = directions_rgba[j].swz<0, 1, 3, 2>();
				if (dot3_s(dirs_rga, dirs_rga) == 0.0f)
				{
					separate_blue_lines[j].b = normalize(vfloat4(1.0f, 1.0f, 1.0f, 0.0f));
				}
				else
				{
					separate_blue_lines[j].b = normalize(dirs_rga);
				}

				separate_alpha_lines[j].a = averages[j].swz<0, 1, 2, 3>();
				separate_alpha_lines[j].a.set_lane<3>(0.0f);

				vfloat4 dirs_rgb = directions_rgba[j].swz<0, 1, 2, 3>();
				if (dot3_s(dirs_rgb, dirs_rgb) == 0.0f)
				{
					separate_alpha_lines[j].b = normalize(vfloat4(1.0f, 1.0f, 1.0f, 0.0f));
				}
				else
				{
					separate_alpha_lines[j].b = normalize(dirs_rgb);
				}

				proc_separate_red_lines[j].amod = (separate_red_lines[j].a - separate_red_lines[j].b * dot3_s(separate_red_lines[j].a, separate_red_lines[j].b)) * inverse_color_scalefactors[j].swz<1, 2, 3, 0>();
				proc_separate_red_lines[j].bs = (separate_red_lines[j].b * color_scalefactors[j].swz<1, 2, 3, 0>());
				proc_separate_red_lines[j].bis = (separate_red_lines[j].b * inverse_color_scalefactors[j].swz<1, 2, 3, 0>());

				proc_separate_green_lines[j].amod =
					(separate_green_lines[j].a - separate_green_lines[j].b * dot3_s(separate_green_lines[j].a, separate_green_lines[j].b)) * inverse_color_scalefactors[j].swz<0, 2, 3, 1>();
				proc_separate_green_lines[j].bs = (separate_green_lines[j].b * color_scalefactors[j].swz<0, 2, 3, 1>());
				proc_separate_green_lines[j].bis = (separate_green_lines[j].b * inverse_color_scalefactors[j].swz<0, 2, 3, 1>());

				proc_separate_blue_lines[j].amod = (separate_blue_lines[j].a - separate_blue_lines[j].b * dot3_s(separate_blue_lines[j].a, separate_blue_lines[j].b)) * inverse_color_scalefactors[j].swz<0, 1, 3, 2>();
				proc_separate_blue_lines[j].bs = (separate_blue_lines[j].b * color_scalefactors[j].swz<0, 1, 3, 2>());
				proc_separate_blue_lines[j].bis = (separate_blue_lines[j].b * inverse_color_scalefactors[j].swz<0, 1, 3, 2>());

				proc_separate_alpha_lines[j].amod =
					(separate_alpha_lines[j].a - separate_alpha_lines[j].b * dot3_s(separate_alpha_lines[j].a, separate_alpha_lines[j].b)) * inverse_color_scalefactors[j].swz<0, 1, 2, 3>();
				proc_separate_alpha_lines[j].bs = (separate_alpha_lines[j].b * color_scalefactors[j].swz<0, 1, 2, 3>());
				proc_separate_alpha_lines[j].bis = (separate_alpha_lines[j].b * inverse_color_scalefactors[j].swz<0, 1, 2, 3>());
			}

			float uncorr_error = 0.0f;
			float samechroma_error = 0.0f;
			vfloat4 separate_error = vfloat4::zero();
			compute_error_squared_rgba(ptab + partition,
			                           blk,
			                           ewb,
			                           proc_uncorr_lines,
			                           proc_samechroma_lines,
			                           uncorr_linelengths,
			                           samechroma_linelengths,
			                           &uncorr_error,
			                           &samechroma_error);

			// compute minimum & maximum alpha values in each partition
			vfloat4 rgba_range[4];
			compute_rgba_range(bsd->texel_count, ptab + partition, blk, ewb, rgba_range);

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

				vfloat4 ics = inverse_color_scalefactors[j];
				vfloat4 error_weights = error_weightings[j] * (tpp * weight_imprecision_estim_squared);

				vfloat4 uncorr_vector = uncorr_lines[j].b * uncorr_linelengths[j] * ics;
				vfloat4 samechroma_vector = samechroma_lines[j].b * samechroma_linelengths[j] * ics;
				vfloat4 separate_red_vector = separate_red_lines[j].b * ics.swz<1, 2, 3, 0>();
				vfloat4 separate_green_vector = separate_green_lines[j].b * ics.swz<0, 2, 3, 1>();
				vfloat4 separate_blue_vector = separate_blue_lines[j].b * ics.swz<0, 1, 3, 2>();
				vfloat4 separate_alpha_vector = separate_alpha_lines[j].b * ics.swz<0, 1, 2, 3>();

				uncorr_vector = uncorr_vector * uncorr_vector;
				samechroma_vector = samechroma_vector * samechroma_vector;
				separate_red_vector = separate_red_vector * separate_red_vector;
				separate_green_vector = separate_green_vector * separate_green_vector;
				separate_blue_vector = separate_blue_vector * separate_blue_vector;
				separate_alpha_vector = separate_alpha_vector * separate_alpha_vector;

				uncorr_error += dot_s(uncorr_vector, error_weights);
				samechroma_error += dot_s(samechroma_vector, error_weights);

				vfloat4 sep_err_inc(dot3_s(separate_red_vector, error_weights.swz<1, 2, 3, 0>()),
				                    dot3_s(separate_green_vector, error_weights.swz<0, 2, 3, 1>()),
				                    dot3_s(separate_blue_vector, error_weights.swz<0, 1, 3, 2>()),
				                    dot3_s(separate_alpha_vector, error_weights.swz<0, 1, 2, 3>()));

				separate_error = separate_error + sep_err_inc + (rgba_range[j] * rgba_range[j] * error_weights);
			}

			if (uncorr_error < uncorr_best_error)
			{
				uncorr_best_error = uncorr_error;
				uncorr_best_partition = partition;
			}

			if (samechroma_error < samechroma_best_errors[0])
			{
				samechroma_best_errors[1] = samechroma_best_errors[0];
				samechroma_best_partitions[1] = samechroma_best_partitions[0];

				samechroma_best_errors[0] = samechroma_error;
				samechroma_best_partitions[0] = partition;
			}
			else if (samechroma_error < samechroma_best_errors[1])
			{
				samechroma_best_errors[1] = samechroma_error;
				samechroma_best_partitions[1] = partition;
			}

			if (separate_error.lane<0>() < separate_best_error)
			{
				separate_best_error = separate_error.lane<0>();
				separate_best_partition = partition;
				separate_best_component = 0;
			}

			if (separate_error.lane<1>() < separate_best_error)
			{
				separate_best_error = separate_error.lane<1>();
				separate_best_partition = partition;
				separate_best_component = 1;
			}

			if (separate_error.lane<2>() < separate_best_error)
			{
				separate_best_error = separate_error.lane<2>();
				separate_best_partition = partition;
				separate_best_component = 2;
			}

			if (separate_error.lane<3>() < separate_best_error)
			{
				separate_best_error = separate_error.lane<3>();
				separate_best_partition = partition;
				separate_best_component = 3;
			}
		}
	}
	else
	{
		for (int i = 0; i < partition_search_limit; i++)
		{
			int partition = partition_sequence[i];

			int bk_partition_count = ptab[partition].partition_count;
			if (bk_partition_count < partition_count)
			{
				continue;
			}

			// compute the weighting to give to each color channel
			// in each partition.
			vfloat4 error_weightings[4];
			vfloat4 color_scalefactors[4];
			vfloat4 inverse_color_scalefactors[4];

			compute_partition_error_color_weightings(bsd, ewb, ptab + partition, error_weightings, color_scalefactors);

			for (int j = 0; j < partition_count; j++)
			{
				inverse_color_scalefactors[j] = 1.0f / max(color_scalefactors[j], 1e-7f);
			}

			vfloat4 averages[4];
			vfloat4 directions_rgb[4];

			compute_averages_and_directions_3_components(ptab + partition, blk, ewb, color_scalefactors, 3, averages, directions_rgb);

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

			for (int j = 0; j < partition_count; j++)
			{
				uncorr_lines[j].a = averages[j];
				if (dot3_s(directions_rgb[j], directions_rgb[j]) == 0.0f)
				{
					uncorr_lines[j].b = normalize(vfloat4(1.0f, 1.0f, 1.0f, 0.0f));
				}
				else
				{
					uncorr_lines[j].b = normalize(directions_rgb[j]);
				}

				samechroma_lines[j].a = vfloat4::zero();
				if (dot3_s(averages[j], averages[j]) == 0.0f)
				{
					samechroma_lines[j].b = normalize(vfloat4(1.0f, 1.0f, 1.0f, 0.0f));
				}
				else
				{
					samechroma_lines[j].b = normalize(averages[j]);
				}

				// TODO: dot3, not 3_s?
				proc_uncorr_lines[j].amod = (uncorr_lines[j].a - uncorr_lines[j].b * dot3_s(uncorr_lines[j].a, uncorr_lines[j].b)) * inverse_color_scalefactors[j].swz<0, 1, 2, 3>();
				proc_uncorr_lines[j].bs = (uncorr_lines[j].b * color_scalefactors[j].swz<0, 1, 2, 3>());
				proc_uncorr_lines[j].bis = (uncorr_lines[j].b * inverse_color_scalefactors[j].swz<0, 1, 2, 3>());

				proc_samechroma_lines[j].amod = (samechroma_lines[j].a - samechroma_lines[j].b * dot(samechroma_lines[j].a, samechroma_lines[j].b)) * inverse_color_scalefactors[j].swz<0, 1, 2, 3>();
				proc_samechroma_lines[j].bs = (samechroma_lines[j].b * color_scalefactors[j].swz<0, 1, 2, 3>());
				proc_samechroma_lines[j].bis = (samechroma_lines[j].b * inverse_color_scalefactors[j].swz<0, 1, 2, 3>());

				separate_red_lines[j].a = averages[j].swz<1, 2>();
				float2 dirs_gb = directions_rgb[j].swz<1, 2>();
				if (dot(dirs_gb, dirs_gb) == 0.0f)
				{
					separate_red_lines[j].b = normalize(float2(1.0f));
				}
				else
				{
					separate_red_lines[j].b = normalize(dirs_gb);
				}

				separate_green_lines[j].a = averages[j].swz<0, 2>();
				float2 dirs_rb = directions_rgb[j].swz<0, 2>();
				if (dot(dirs_rb, dirs_rb) == 0.0f)
				{
					separate_green_lines[j].b = normalize(float2(1.0f));
				}
				else
				{
					separate_green_lines[j].b = normalize(dirs_rb);
				}

				separate_blue_lines[j].a = averages[j].swz<0, 1>();
				float2 dirs_rg = directions_rgb[j].swz<0, 1>();
				if (dot(dirs_rg, dirs_rg) == 0.0f)
				{
					separate_blue_lines[j].b = normalize(float2(1.0f));
				}
				else
				{
					separate_blue_lines[j].b = normalize(dirs_rg);
				}

				proc_separate_red_lines[j].amod = (separate_red_lines[j].a - separate_red_lines[j].b * dot(separate_red_lines[j].a, separate_red_lines[j].b)) * inverse_color_scalefactors[j].swz<1, 2>();
				proc_separate_red_lines[j].bs = (separate_red_lines[j].b * color_scalefactors[j].swz<1, 2>());
				proc_separate_red_lines[j].bis = (separate_red_lines[j].b *  inverse_color_scalefactors[j].swz<1, 2>());

				proc_separate_green_lines[j].amod = (separate_green_lines[j].a - separate_green_lines[j].b * dot(separate_green_lines[j].a, separate_green_lines[j].b)) *  inverse_color_scalefactors[j].swz<0, 2>();
				proc_separate_green_lines[j].bs = (separate_green_lines[j].b * color_scalefactors[j].swz<0, 2>());
				proc_separate_green_lines[j].bis = (separate_green_lines[j].b *  inverse_color_scalefactors[j].swz<0, 2>());

				proc_separate_blue_lines[j].amod = (separate_blue_lines[j].a - separate_blue_lines[j].b * dot(separate_blue_lines[j].a, separate_blue_lines[j].b)) *  inverse_color_scalefactors[j].swz<0, 1>();
				proc_separate_blue_lines[j].bs = (separate_blue_lines[j].b * color_scalefactors[j].swz<0, 1>());
				proc_separate_blue_lines[j].bis = (separate_blue_lines[j].b *  inverse_color_scalefactors[j].swz<0, 1>());
			}

			float uncorr_error = 0.0f;
			float samechroma_error = 0.0f;
			vfloat4 separate_error = vfloat4(0.0f);

			compute_error_squared_rgb(ptab + partition,
			                          blk,
			                          ewb,
			                          proc_uncorr_lines,
			                          proc_samechroma_lines,
			                          uncorr_linelengths,
			                          samechroma_linelengths,
			                          &uncorr_error,
			                          &samechroma_error);

			vfloat4 rgba_range[4];
			compute_rgba_range(bsd->texel_count, ptab + partition, blk, ewb, rgba_range);

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

				vfloat4 ics = inverse_color_scalefactors[j];
				ics.set_lane<3>(0.0f);
				vfloat4 error_weights = error_weightings[j] * (tpp * weight_imprecision_estim_squared);
				error_weights.set_lane<3>(0.0f);

				vfloat4 uncorr_vector = (uncorr_lines[j].b * uncorr_linelengths[j]) * ics;
				vfloat4 samechroma_vector = (samechroma_lines[j].b * samechroma_linelengths[j]) * ics;

				float2 separate_red_vector = separate_red_lines[j].b * ics.swz<1, 2>();
				float2 separate_green_vector = separate_green_lines[j].b * ics.swz<0, 2>();
				float2 separate_blue_vector = separate_blue_lines[j].b * ics.swz<0, 1>();

				uncorr_vector = uncorr_vector * uncorr_vector;
				samechroma_vector = samechroma_vector * samechroma_vector;
				separate_red_vector = separate_red_vector * separate_red_vector;
				separate_green_vector = separate_green_vector * separate_green_vector;
				separate_blue_vector = separate_blue_vector * separate_blue_vector;

				uncorr_error += dot3_s(uncorr_vector, error_weights);
				samechroma_error += dot3_s(samechroma_vector, error_weights);
				separate_error.set_lane<0>(separate_error.lane<0>() + dot(separate_red_vector, error_weights.swz<1, 2>()));
				separate_error.set_lane<1>(separate_error.lane<1>() + dot(separate_green_vector, error_weights.swz<0, 2>()));
				separate_error.set_lane<2>(separate_error.lane<2>() + dot(separate_blue_vector, error_weights.swz<0, 1>()));

				separate_error.set_lane<0>(separate_error.lane<0>() + rgba_range[j].lane<0>() * rgba_range[j].lane<0>() * error_weights.lane<0>());
				separate_error.set_lane<1>(separate_error.lane<1>() + rgba_range[j].lane<1>() * rgba_range[j].lane<1>() * error_weights.lane<1>());
				separate_error.set_lane<2>(separate_error.lane<2>() + rgba_range[j].lane<2>() * rgba_range[j].lane<2>() * error_weights.lane<2>());
			}

			if (uncorr_error < uncorr_best_error)
			{
				uncorr_best_error = uncorr_error;
				uncorr_best_partition = partition;
			}

			if (samechroma_error < samechroma_best_errors[0])
			{
				samechroma_best_errors[1] = samechroma_best_errors[0];
				samechroma_best_partitions[1] = samechroma_best_partitions[0];

				samechroma_best_errors[0] = samechroma_error;
				samechroma_best_partitions[0] = partition;
			}
			else if (samechroma_error < samechroma_best_errors[1])
			{
				samechroma_best_errors[1] = samechroma_error;
				samechroma_best_partitions[1] = partition;
			}

			if (separate_error.lane<0>() < separate_best_error)
			{
				separate_best_error = separate_error.lane<0>();
				separate_best_partition = partition;
				separate_best_component = 0;
			}

			if (separate_error.lane<1>() < separate_best_error)
			{
				separate_best_error = separate_error.lane<1>();
				separate_best_partition = partition;
				separate_best_component = 1;
			}

			if (separate_error.lane<2>() < separate_best_error)
			{
				separate_best_error = separate_error.lane<2>();
				separate_best_partition = partition;
				separate_best_component = 2;
			}
		}
	}

	*best_partition_uncorrelated = uncorr_best_partition;

	int index { samechroma_best_partitions[0] != uncorr_best_partition ? 0 : 1 };
	*best_partition_samechroma = samechroma_best_partitions[index];

	*best_partition_dualplane = (separate_best_component << PARTITION_BITS) |
	                            (separate_best_partition);
}

#endif
