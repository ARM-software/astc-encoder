// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2011-2026 Arm Limited
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
 * @brief Functions for finding dominant direction of a set of colors.
 */
#if !defined(ASTCENC_DECOMPRESS_ONLY)

#include "astcenc_internal.h"

#include <cassert>

/**
 * @brief Compute the average RGB color of each partition.
 *
 * The algorithm here uses a vectorized sequential scan and per-partition
 * color accumulators, using select() to mask texel lanes in other partitions.
 *
 * We only accumulate sums for N-1 partitions during the scan; the value for
 * the last partition can be computed given that we know the block-wide average
 * already.
 *
 * Because of this we could reduce the loop iteration count so it "just" spans
 * the max texel index needed for the N-1 partitions, which could need fewer
 * iterations than the full block texel count. However, this makes the loop
 * count erratic and causes more branch mispredictions so is a net loss.
 *
 * @param      pi         The partitioning to use.
 * @param      blk        The block data to process.
 * @param[out] averages   The output averages. Unused partition indices will
 *                        not be initialized, and lane<3> will be zero.
 */
static void compute_partition_averages_rgb(
	const partition_info& pi,
	const image_block& blk,
	vfloat4 averages[BLOCK_MAX_PARTITIONS]
) {
	unsigned int partition_count = pi.partition_count;
	size_t texel_count = blk.texel_count;
	promise(texel_count > 0);

	// Use the precomputed mean for a single partition
	if (partition_count == 1)
	{
		averages[0] = blk.data_mean.swz<0, 1, 2>();
	}
	// For 2 partitions scan results for partition 0, compute partition 1
	else if (partition_count == 2)
	{
		vfloatacc pp_avg_rgb[3] {};

		vint lane_id = vint::lane_id();
		for (size_t i = 0; i < texel_count; i += ASTCENC_SIMD_WIDTH)
		{
			vint texel_partition(pi.partition_of_texel + i);

			vmask lane_mask = lane_id < vint_from_size(texel_count);
			lane_id += vint(ASTCENC_SIMD_WIDTH);

			vmask p0_mask = lane_mask & (texel_partition == vint(0));

			vfloat data_r = loada(blk.data_r + i);
			haccumulate(pp_avg_rgb[0], data_r, p0_mask);

			vfloat data_g = loada(blk.data_g + i);
			haccumulate(pp_avg_rgb[1], data_g, p0_mask);

			vfloat data_b = loada(blk.data_b + i);
			haccumulate(pp_avg_rgb[2], data_b, p0_mask);
		}

		vfloat4 block_total = blk.data_mean.swz<0, 1, 2>() * static_cast<float>(blk.texel_count);

		vfloat4 p0_total = vfloat3(hadd_s(pp_avg_rgb[0]),
		                           hadd_s(pp_avg_rgb[1]),
		                           hadd_s(pp_avg_rgb[2]));

		vfloat4 p1_total = block_total - p0_total;

		averages[0] = p0_total / static_cast<float>(pi.partition_texel_count[0]);
		averages[1] = p1_total / static_cast<float>(pi.partition_texel_count[1]);
	}
	// For 3 partitions scan results for partition 0/1, compute partition 2
	else if (partition_count == 3)
	{
		vfloatacc pp_avg_rgb[2][3] {};

		vint lane_id = vint::lane_id();
		for (size_t i = 0; i < texel_count; i += ASTCENC_SIMD_WIDTH)
		{
			vint texel_partition(pi.partition_of_texel + i);

			vmask lane_mask = lane_id < vint_from_size(texel_count);
			lane_id += vint(ASTCENC_SIMD_WIDTH);

			vmask p0_mask = lane_mask & (texel_partition == vint(0));
			vmask p1_mask = lane_mask & (texel_partition == vint(1));

			vfloat data_r = loada(blk.data_r + i);
			haccumulate(pp_avg_rgb[0][0], data_r, p0_mask);
			haccumulate(pp_avg_rgb[1][0], data_r, p1_mask);

			vfloat data_g = loada(blk.data_g + i);
			haccumulate(pp_avg_rgb[0][1], data_g, p0_mask);
			haccumulate(pp_avg_rgb[1][1], data_g, p1_mask);

			vfloat data_b = loada(blk.data_b + i);
			haccumulate(pp_avg_rgb[0][2], data_b, p0_mask);
			haccumulate(pp_avg_rgb[1][2], data_b, p1_mask);
		}

		vfloat4 block_total = blk.data_mean.swz<0, 1, 2>() * static_cast<float>(blk.texel_count);

		vfloat4 p0_total = vfloat3(hadd_s(pp_avg_rgb[0][0]),
		                           hadd_s(pp_avg_rgb[0][1]),
		                           hadd_s(pp_avg_rgb[0][2]));

		vfloat4 p1_total = vfloat3(hadd_s(pp_avg_rgb[1][0]),
		                           hadd_s(pp_avg_rgb[1][1]),
		                           hadd_s(pp_avg_rgb[1][2]));

		vfloat4 p2_total = block_total - p0_total - p1_total;

		averages[0] = p0_total / static_cast<float>(pi.partition_texel_count[0]);
		averages[1] = p1_total / static_cast<float>(pi.partition_texel_count[1]);
		averages[2] = p2_total / static_cast<float>(pi.partition_texel_count[2]);
	}
	else
	{
		// For 4 partitions scan results for partition 0/1/2, compute partition 3
		vfloatacc pp_avg_rgb[3][3] {};

		vint lane_id = vint::lane_id();
		for (size_t i = 0; i < texel_count; i += ASTCENC_SIMD_WIDTH)
		{
			vint texel_partition(pi.partition_of_texel + i);

			vmask lane_mask = lane_id < vint_from_size(texel_count);
			lane_id += vint(ASTCENC_SIMD_WIDTH);

			vmask p0_mask = lane_mask & (texel_partition == vint(0));
			vmask p1_mask = lane_mask & (texel_partition == vint(1));
			vmask p2_mask = lane_mask & (texel_partition == vint(2));

			vfloat data_r = loada(blk.data_r + i);
			haccumulate(pp_avg_rgb[0][0], data_r, p0_mask);
			haccumulate(pp_avg_rgb[1][0], data_r, p1_mask);
			haccumulate(pp_avg_rgb[2][0], data_r, p2_mask);

			vfloat data_g = loada(blk.data_g + i);
			haccumulate(pp_avg_rgb[0][1], data_g, p0_mask);
			haccumulate(pp_avg_rgb[1][1], data_g, p1_mask);
			haccumulate(pp_avg_rgb[2][1], data_g, p2_mask);

			vfloat data_b = loada(blk.data_b + i);
			haccumulate(pp_avg_rgb[0][2], data_b, p0_mask);
			haccumulate(pp_avg_rgb[1][2], data_b, p1_mask);
			haccumulate(pp_avg_rgb[2][2], data_b, p2_mask);
		}

		vfloat4 block_total = blk.data_mean.swz<0, 1, 2>() * static_cast<float>(blk.texel_count);

		vfloat4 p0_total = vfloat3(hadd_s(pp_avg_rgb[0][0]),
		                           hadd_s(pp_avg_rgb[0][1]),
		                           hadd_s(pp_avg_rgb[0][2]));

		vfloat4 p1_total = vfloat3(hadd_s(pp_avg_rgb[1][0]),
		                           hadd_s(pp_avg_rgb[1][1]),
		                           hadd_s(pp_avg_rgb[1][2]));

		vfloat4 p2_total = vfloat3(hadd_s(pp_avg_rgb[2][0]),
		                           hadd_s(pp_avg_rgb[2][1]),
		                           hadd_s(pp_avg_rgb[2][2]));

		vfloat4 p3_total = block_total - p0_total - p1_total- p2_total;

		averages[0] = p0_total / static_cast<float>(pi.partition_texel_count[0]);
		averages[1] = p1_total / static_cast<float>(pi.partition_texel_count[1]);
		averages[2] = p2_total / static_cast<float>(pi.partition_texel_count[2]);
		averages[3] = p3_total / static_cast<float>(pi.partition_texel_count[3]);
	}
}

/**
 * @brief Compute the average RGBA color of each partition.
 *
 * The algorithm here uses a vectorized sequential scan and per-partition
 * color accumulators, using select() to mask texel lanes in other partitions.
 *
 * We only accumulate sums for N-1 partitions during the scan; the value for
 * the last partition can be computed given that we know the block-wide average
 * already.
 *
 * Because of this we could reduce the loop iteration count so it "just" spans
 * the max texel index needed for the N-1 partitions, which could need fewer
 * iterations than the full block texel count. However, this makes the loop
 * count erratic and causes more branch mispredictions so is a net loss.
 *
 * @param      pi         The partitioning to use.
 * @param      blk        The block data to process.
 * @param[out] averages   The output averages. Unused partition indices will
 *                        not be initialized.
 */
static void compute_partition_averages_rgba(
	const partition_info& pi,
	const image_block& blk,
	vfloat4 averages[BLOCK_MAX_PARTITIONS]
) {
	unsigned int partition_count = pi.partition_count;
	size_t texel_count = blk.texel_count;
	promise(texel_count > 0);

	// Use the precomputed mean for a single partition
	if (partition_count == 1)
	{
		averages[0] = blk.data_mean;
	}
	// For 2 partitions scan results for partition 0, compute partition 1
	else if (partition_count == 2)
	{
		vfloat4 pp_avg_rgba[4] {};

		vint lane_id = vint::lane_id();
		for (size_t i = 0; i < texel_count; i += ASTCENC_SIMD_WIDTH)
		{
			vint texel_partition(pi.partition_of_texel + i);

			vmask lane_mask = lane_id < vint_from_size(texel_count);
			lane_id += vint(ASTCENC_SIMD_WIDTH);

			vmask p0_mask = lane_mask & (texel_partition == vint(0));

			vfloat data_r = loada(blk.data_r + i);
			haccumulate(pp_avg_rgba[0], data_r, p0_mask);

			vfloat data_g = loada(blk.data_g + i);
			haccumulate(pp_avg_rgba[1], data_g, p0_mask);

			vfloat data_b = loada(blk.data_b + i);
			haccumulate(pp_avg_rgba[2], data_b, p0_mask);

			vfloat data_a = loada(blk.data_a + i);
			haccumulate(pp_avg_rgba[3], data_a, p0_mask);
		}

		vfloat4 block_total = blk.data_mean * static_cast<float>(blk.texel_count);

		vfloat4 p0_total = vfloat4(hadd_s(pp_avg_rgba[0]),
		                           hadd_s(pp_avg_rgba[1]),
		                           hadd_s(pp_avg_rgba[2]),
		                           hadd_s(pp_avg_rgba[3]));

		vfloat4 p1_total = block_total - p0_total;

		averages[0] = p0_total / static_cast<float>(pi.partition_texel_count[0]);
		averages[1] = p1_total / static_cast<float>(pi.partition_texel_count[1]);
	}
	// For 3 partitions scan results for partition 0/1, compute partition 2
	else if (partition_count == 3)
	{
		vfloat4 pp_avg_rgba[2][4] {};

		vint lane_id = vint::lane_id();
		for (size_t i = 0; i < texel_count; i += ASTCENC_SIMD_WIDTH)
		{
			vint texel_partition(pi.partition_of_texel + i);

			vmask lane_mask = lane_id < vint_from_size(texel_count);
			lane_id += vint(ASTCENC_SIMD_WIDTH);

			vmask p0_mask = lane_mask & (texel_partition == vint(0));
			vmask p1_mask = lane_mask & (texel_partition == vint(1));

			vfloat data_r = loada(blk.data_r + i);
			haccumulate(pp_avg_rgba[0][0], data_r, p0_mask);
			haccumulate(pp_avg_rgba[1][0], data_r, p1_mask);

			vfloat data_g = loada(blk.data_g + i);
			haccumulate(pp_avg_rgba[0][1], data_g, p0_mask);
			haccumulate(pp_avg_rgba[1][1], data_g, p1_mask);

			vfloat data_b = loada(blk.data_b + i);
			haccumulate(pp_avg_rgba[0][2], data_b, p0_mask);
			haccumulate(pp_avg_rgba[1][2], data_b, p1_mask);

			vfloat data_a = loada(blk.data_a + i);
			haccumulate(pp_avg_rgba[0][3], data_a, p0_mask);
			haccumulate(pp_avg_rgba[1][3], data_a, p1_mask);
		}

		vfloat4 block_total = blk.data_mean * static_cast<float>(blk.texel_count);

		vfloat4 p0_total = vfloat4(hadd_s(pp_avg_rgba[0][0]),
		                           hadd_s(pp_avg_rgba[0][1]),
		                           hadd_s(pp_avg_rgba[0][2]),
		                           hadd_s(pp_avg_rgba[0][3]));

		vfloat4 p1_total = vfloat4(hadd_s(pp_avg_rgba[1][0]),
		                           hadd_s(pp_avg_rgba[1][1]),
		                           hadd_s(pp_avg_rgba[1][2]),
		                           hadd_s(pp_avg_rgba[1][3]));

		vfloat4 p2_total = block_total - p0_total - p1_total;

		averages[0] = p0_total / static_cast<float>(pi.partition_texel_count[0]);
		averages[1] = p1_total / static_cast<float>(pi.partition_texel_count[1]);
		averages[2] = p2_total / static_cast<float>(pi.partition_texel_count[2]);
	}
	else
	{
		// For 4 partitions scan results for partition 0/1/2, compute partition 3
		vfloat4 pp_avg_rgba[3][4] {};

		vint lane_id = vint::lane_id();
		for (size_t i = 0; i < texel_count; i += ASTCENC_SIMD_WIDTH)
		{
			vint texel_partition(pi.partition_of_texel + i);

			vmask lane_mask = lane_id < vint_from_size(texel_count);
			lane_id += vint(ASTCENC_SIMD_WIDTH);

			vmask p0_mask = lane_mask & (texel_partition == vint(0));
			vmask p1_mask = lane_mask & (texel_partition == vint(1));
			vmask p2_mask = lane_mask & (texel_partition == vint(2));

			vfloat data_r = loada(blk.data_r + i);
			haccumulate(pp_avg_rgba[0][0], data_r, p0_mask);
			haccumulate(pp_avg_rgba[1][0], data_r, p1_mask);
			haccumulate(pp_avg_rgba[2][0], data_r, p2_mask);

			vfloat data_g = loada(blk.data_g + i);
			haccumulate(pp_avg_rgba[0][1], data_g, p0_mask);
			haccumulate(pp_avg_rgba[1][1], data_g, p1_mask);
			haccumulate(pp_avg_rgba[2][1], data_g, p2_mask);

			vfloat data_b = loada(blk.data_b + i);
			haccumulate(pp_avg_rgba[0][2], data_b, p0_mask);
			haccumulate(pp_avg_rgba[1][2], data_b, p1_mask);
			haccumulate(pp_avg_rgba[2][2], data_b, p2_mask);

			vfloat data_a = loada(blk.data_a + i);
			haccumulate(pp_avg_rgba[0][3], data_a, p0_mask);
			haccumulate(pp_avg_rgba[1][3], data_a, p1_mask);
			haccumulate(pp_avg_rgba[2][3], data_a, p2_mask);
		}

		vfloat4 block_total = blk.data_mean * static_cast<float>(blk.texel_count);

		vfloat4 p0_total = vfloat4(hadd_s(pp_avg_rgba[0][0]),
		                           hadd_s(pp_avg_rgba[0][1]),
		                           hadd_s(pp_avg_rgba[0][2]),
		                           hadd_s(pp_avg_rgba[0][3]));

		vfloat4 p1_total = vfloat4(hadd_s(pp_avg_rgba[1][0]),
		                           hadd_s(pp_avg_rgba[1][1]),
		                           hadd_s(pp_avg_rgba[1][2]),
		                           hadd_s(pp_avg_rgba[1][3]));

		vfloat4 p2_total = vfloat4(hadd_s(pp_avg_rgba[2][0]),
		                           hadd_s(pp_avg_rgba[2][1]),
		                           hadd_s(pp_avg_rgba[2][2]),
		                           hadd_s(pp_avg_rgba[2][3]));

		vfloat4 p3_total = block_total - p0_total - p1_total- p2_total;

		averages[0] = p0_total / static_cast<float>(pi.partition_texel_count[0]);
		averages[1] = p1_total / static_cast<float>(pi.partition_texel_count[1]);
		averages[2] = p2_total / static_cast<float>(pi.partition_texel_count[2]);
		averages[3] = p3_total / static_cast<float>(pi.partition_texel_count[3]);
	}
}

/**
 * @brief Compute the average RGBA color of each partition for grayscale RGB.
 *
 * Accumulate only R and A, reusing the R total for G and B. Preserve the
 * general RGBA path's accumulation order and final-partition subtraction.
 *
 * @param      pi         The partitioning to use.
 * @param      blk        The block data to process, with identical RGB components.
 * @param[out] averages   The partition averages. Unused entries will not be initialized.
 */
static void compute_partition_averages_rgba_grayscale(
	const partition_info& pi,
	const image_block& blk,
	vfloat4 averages[BLOCK_MAX_PARTITIONS]
) {
	unsigned int partition_count = pi.partition_count;
	size_t texel_count = blk.texel_count;
	promise(texel_count > 0);

	// Use the precomputed mean for a single partition
	if (partition_count == 1)
	{
		averages[0] = blk.data_mean;
	}
	// For 2 partitions scan results for partition 0, compute partition 1
	else if (partition_count == 2)
	{
		vfloat4 pp_avg_ra[2] {};

		vint lane_id = vint::lane_id();
		for (size_t i = 0; i < texel_count; i += ASTCENC_SIMD_WIDTH)
		{
			vint texel_partition(pi.partition_of_texel + i);

			vmask lane_mask = lane_id < vint_from_size(texel_count);
			lane_id += vint(ASTCENC_SIMD_WIDTH);

			vmask p0_mask = lane_mask & (texel_partition == vint(0));

			vfloat data_r = loada(blk.data_r + i);
			haccumulate(pp_avg_ra[0], data_r, p0_mask);

			vfloat data_a = loada(blk.data_a + i);
			haccumulate(pp_avg_ra[1], data_a, p0_mask);
		}

		vfloat4 block_total = blk.data_mean * static_cast<float>(blk.texel_count);

		float p0_r = hadd_s(pp_avg_ra[0]);
		vfloat4 p0_total(p0_r, p0_r, p0_r, hadd_s(pp_avg_ra[1]));

		vfloat4 p1_total = block_total - p0_total;

		averages[0] = p0_total / static_cast<float>(pi.partition_texel_count[0]);
		averages[1] = p1_total / static_cast<float>(pi.partition_texel_count[1]);
	}
	// For 3 partitions scan results for partition 0/1, compute partition 2
	else if (partition_count == 3)
	{
		vfloat4 pp_avg_ra[2][2] {};

		vint lane_id = vint::lane_id();
		for (size_t i = 0; i < texel_count; i += ASTCENC_SIMD_WIDTH)
		{
			vint texel_partition(pi.partition_of_texel + i);

			vmask lane_mask = lane_id < vint_from_size(texel_count);
			lane_id += vint(ASTCENC_SIMD_WIDTH);

			vmask p0_mask = lane_mask & (texel_partition == vint(0));
			vmask p1_mask = lane_mask & (texel_partition == vint(1));

			vfloat data_r = loada(blk.data_r + i);
			haccumulate(pp_avg_ra[0][0], data_r, p0_mask);
			haccumulate(pp_avg_ra[1][0], data_r, p1_mask);

			vfloat data_a = loada(blk.data_a + i);
			haccumulate(pp_avg_ra[0][1], data_a, p0_mask);
			haccumulate(pp_avg_ra[1][1], data_a, p1_mask);
		}

		vfloat4 block_total = blk.data_mean * static_cast<float>(blk.texel_count);

		float p0_r = hadd_s(pp_avg_ra[0][0]);
		vfloat4 p0_total(p0_r, p0_r, p0_r, hadd_s(pp_avg_ra[0][1]));

		float p1_r = hadd_s(pp_avg_ra[1][0]);
		vfloat4 p1_total(p1_r, p1_r, p1_r, hadd_s(pp_avg_ra[1][1]));

		vfloat4 p2_total = block_total - p0_total - p1_total;

		averages[0] = p0_total / static_cast<float>(pi.partition_texel_count[0]);
		averages[1] = p1_total / static_cast<float>(pi.partition_texel_count[1]);
		averages[2] = p2_total / static_cast<float>(pi.partition_texel_count[2]);
	}
	else
	{
		// For 4 partitions scan results for partition 0/1/2, compute partition 3
		vfloat4 pp_avg_ra[3][2] {};

		vint lane_id = vint::lane_id();
		for (size_t i = 0; i < texel_count; i += ASTCENC_SIMD_WIDTH)
		{
			vint texel_partition(pi.partition_of_texel + i);

			vmask lane_mask = lane_id < vint_from_size(texel_count);
			lane_id += vint(ASTCENC_SIMD_WIDTH);

			vmask p0_mask = lane_mask & (texel_partition == vint(0));
			vmask p1_mask = lane_mask & (texel_partition == vint(1));
			vmask p2_mask = lane_mask & (texel_partition == vint(2));

			vfloat data_r = loada(blk.data_r + i);
			haccumulate(pp_avg_ra[0][0], data_r, p0_mask);
			haccumulate(pp_avg_ra[1][0], data_r, p1_mask);
			haccumulate(pp_avg_ra[2][0], data_r, p2_mask);

			vfloat data_a = loada(blk.data_a + i);
			haccumulate(pp_avg_ra[0][1], data_a, p0_mask);
			haccumulate(pp_avg_ra[1][1], data_a, p1_mask);
			haccumulate(pp_avg_ra[2][1], data_a, p2_mask);
		}

		vfloat4 block_total = blk.data_mean * static_cast<float>(blk.texel_count);

		float p0_r = hadd_s(pp_avg_ra[0][0]);
		vfloat4 p0_total(p0_r, p0_r, p0_r, hadd_s(pp_avg_ra[0][1]));

		float p1_r = hadd_s(pp_avg_ra[1][0]);
		vfloat4 p1_total(p1_r, p1_r, p1_r, hadd_s(pp_avg_ra[1][1]));

		float p2_r = hadd_s(pp_avg_ra[2][0]);
		vfloat4 p2_total(p2_r, p2_r, p2_r, hadd_s(pp_avg_ra[2][1]));

		vfloat4 p3_total = block_total - p0_total - p1_total- p2_total;

		averages[0] = p0_total / static_cast<float>(pi.partition_texel_count[0]);
		averages[1] = p1_total / static_cast<float>(pi.partition_texel_count[1]);
		averages[2] = p2_total / static_cast<float>(pi.partition_texel_count[2]);
		averages[3] = p3_total / static_cast<float>(pi.partition_texel_count[3]);
	}
}

/**
 * @brief Compute RGBA partition averages and directions for grayscale RGB.
 *
 * Skips computation of G and B, given we know it will be the same as R.
 *
 * @param      pi            The partitioning to use.
 * @param      blk           The image block color data to process.
 * @param[out] pm            The partition averages and directions.
 * @param[out] rgb_metrics   The RGB partition metrics, or @c nullptr to omit.
 */
static void compute_avgs_and_dirs_4_comp_grayscale(
	const partition_info& pi,
	const image_block& blk,
	partition_metrics pm[BLOCK_MAX_PARTITIONS],
	partition_metrics* rgb_metrics
) {
	size_t partition_count = pi.partition_count;
	promise(partition_count > 0);

	// Pre-compute partition_averages
	vfloat4 partition_averages[BLOCK_MAX_PARTITIONS];
	compute_partition_averages_rgba_grayscale(pi, blk, partition_averages);

	for (size_t partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *texel_indexes = pi.texels_of_partition[partition];
		size_t texel_count = pi.partition_texel_count[partition];
		promise(texel_count > 0);

		vfloat4 average = partition_averages[partition];
		pm[partition].avg = average;

		vfloat4 sum_xp = vfloat4::zero();
		vfloat4 sum_wp = vfloat4::zero();

		for (size_t i = 0; i < texel_count; i++)
		{
			unsigned int iwt = texel_indexes[i];
			float data_r = blk.data_r[iwt];
			float data_a = blk.data_a[iwt];
			vfloat4 texel_datum(data_r, data_r, data_r, data_a);
			texel_datum = texel_datum - average;

			vfloat4 zero = vfloat4::zero();

			vmask4 tdm0 = texel_datum.swz<0,0,0,0>() > zero;
			sum_xp += select(zero, texel_datum, tdm0);

			vmask4 tdm3 = texel_datum.swz<3,3,3,3>() > zero;
			sum_wp += select(zero, texel_datum, tdm3);
		}

		// Preserve the RGB-only direction before alpha can change the RGBA winner
		if (rgb_metrics)
		{
			rgb_metrics[partition].avg = average.swz<0, 1, 2>();
			rgb_metrics[partition].dir = sum_xp.swz<0, 1, 2>();
		}

		vfloat4 prod_xp = dot(sum_xp, sum_xp);
		vfloat4 prod_wp = dot(sum_wp, sum_wp);

		vfloat4 best_vector = sum_xp;

		vmask4 mask = prod_wp > prod_xp;
		best_vector = select(best_vector, sum_wp, mask);

		pm[partition].dir = best_vector;
	}
}

/**
 * @brief Compute RGBA partition averages and directions for general RGBA data.
 *
 * @param      pi    The partitioning to use.
 * @param      blk   The image block color data to process.
 * @param[out] pm    The partition averages and directions.
 */
static void compute_avgs_and_dirs_4_comp_general(
	const partition_info& pi,
	const image_block& blk,
	partition_metrics pm[BLOCK_MAX_PARTITIONS]
) {
	size_t partition_count = pi.partition_count;
	promise(partition_count > 0);

	// Pre-compute partition_averages
	vfloat4 partition_averages[BLOCK_MAX_PARTITIONS];
	compute_partition_averages_rgba(pi, blk, partition_averages);

	for (size_t partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *texel_indexes = pi.texels_of_partition[partition];
		size_t texel_count = pi.partition_texel_count[partition];
		promise(texel_count > 0);

		vfloat4 average = partition_averages[partition];
		pm[partition].avg = average;

		vfloat4 sum_xp = vfloat4::zero();
		vfloat4 sum_yp = vfloat4::zero();
		vfloat4 sum_zp = vfloat4::zero();
		vfloat4 sum_wp = vfloat4::zero();

		for (size_t i = 0; i < texel_count; i++)
		{
			unsigned int iwt = texel_indexes[i];
			vfloat4 texel_datum = blk.texel(iwt);
			texel_datum = texel_datum - average;

			vfloat4 zero = vfloat4::zero();

			vmask4 tdm0 = texel_datum.swz<0,0,0,0>() > zero;
			sum_xp += select(zero, texel_datum, tdm0);

			vmask4 tdm1 = texel_datum.swz<1,1,1,1>() > zero;
			sum_yp += select(zero, texel_datum, tdm1);

			vmask4 tdm2 = texel_datum.swz<2,2,2,2>() > zero;
			sum_zp += select(zero, texel_datum, tdm2);

			vmask4 tdm3 = texel_datum.swz<3,3,3,3>() > zero;
			sum_wp += select(zero, texel_datum, tdm3);
		}

		vfloat4 prod_xp = dot(sum_xp, sum_xp);
		vfloat4 prod_yp = dot(sum_yp, sum_yp);
		vfloat4 prod_zp = dot(sum_zp, sum_zp);
		vfloat4 prod_wp = dot(sum_wp, sum_wp);

		vfloat4 best_vector = sum_xp;
		vfloat4 best_sum = prod_xp;

		vmask4 mask = prod_yp > best_sum;
		best_vector = select(best_vector, sum_yp, mask);
		best_sum = select(best_sum, prod_yp, mask);

		mask = prod_zp > best_sum;
		best_vector = select(best_vector, sum_zp, mask);
		best_sum = select(best_sum, prod_zp, mask);

		mask = prod_wp > best_sum;
		best_vector = select(best_vector, sum_wp, mask);

		pm[partition].dir = best_vector;
	}
}

/* See header for documentation. */
void compute_avgs_and_dirs_4_comp(
	const partition_info& pi,
	const image_block& blk,
	partition_metrics pm[BLOCK_MAX_PARTITIONS],
	partition_metrics* rgb_metrics
) {
	if (blk.grayscale)
	{
		compute_avgs_and_dirs_4_comp_grayscale(pi, blk, pm, rgb_metrics);
	}
	else
	{
		compute_avgs_and_dirs_4_comp_general(pi, blk, pm);
	}
}

/* See header for documentation. */
void compute_avgs_and_dirs_3_comp(
	const partition_info& pi,
	const image_block& blk,
	unsigned int omitted_component,
	partition_metrics pm[BLOCK_MAX_PARTITIONS]
) {
	// Pre-compute partition_averages
	vfloat4 partition_averages[BLOCK_MAX_PARTITIONS];
	compute_partition_averages_rgba(pi, blk, partition_averages);

	const float* data_vr = blk.data_r;
	const float* data_vg = blk.data_g;
	const float* data_vb = blk.data_b;

	if (omitted_component == 0)
	{
		partition_averages[0] = partition_averages[0].swz<1, 2, 3>();
		partition_averages[1] = partition_averages[1].swz<1, 2, 3>();
		partition_averages[2] = partition_averages[2].swz<1, 2, 3>();
		partition_averages[3] = partition_averages[3].swz<1, 2, 3>();

		data_vr = blk.data_g;
		data_vg = blk.data_b;
		data_vb = blk.data_a;
	}
	else if (omitted_component == 1)
	{
		partition_averages[0] = partition_averages[0].swz<0, 2, 3>();
		partition_averages[1] = partition_averages[1].swz<0, 2, 3>();
		partition_averages[2] = partition_averages[2].swz<0, 2, 3>();
		partition_averages[3] = partition_averages[3].swz<0, 2, 3>();

		data_vg = blk.data_b;
		data_vb = blk.data_a;
	}
	else if (omitted_component == 2)
	{
		partition_averages[0] = partition_averages[0].swz<0, 1, 3>();
		partition_averages[1] = partition_averages[1].swz<0, 1, 3>();
		partition_averages[2] = partition_averages[2].swz<0, 1, 3>();
		partition_averages[3] = partition_averages[3].swz<0, 1, 3>();

		data_vb = blk.data_a;
	}
	else
	{
		partition_averages[0] = partition_averages[0].swz<0, 1, 2>();
		partition_averages[1] = partition_averages[1].swz<0, 1, 2>();
		partition_averages[2] = partition_averages[2].swz<0, 1, 2>();
		partition_averages[3] = partition_averages[3].swz<0, 1, 2>();
	}

	size_t partition_count = pi.partition_count;
	promise(partition_count > 0);

	for (size_t partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *texel_indexes = pi.texels_of_partition[partition];
		size_t texel_count = pi.partition_texel_count[partition];
		promise(texel_count > 0);

		vfloat4 average = partition_averages[partition];
		pm[partition].avg = average;

		vfloat4 sum_xp = vfloat4::zero();
		vfloat4 sum_yp = vfloat4::zero();
		vfloat4 sum_zp = vfloat4::zero();

		for (size_t i = 0; i < texel_count; i++)
		{
			unsigned int iwt = texel_indexes[i];

			vfloat4 texel_datum = vfloat3(data_vr[iwt],
			                              data_vg[iwt],
			                              data_vb[iwt]);
			texel_datum = texel_datum - average;

			vfloat4 zero = vfloat4::zero();

			vmask4 tdm0 = texel_datum.swz<0,0,0,0>() > zero;
			sum_xp += select(zero, texel_datum, tdm0);

			vmask4 tdm1 = texel_datum.swz<1,1,1,1>() > zero;
			sum_yp += select(zero, texel_datum, tdm1);

			vmask4 tdm2 = texel_datum.swz<2,2,2,2>() > zero;
			sum_zp += select(zero, texel_datum, tdm2);
		}

		vfloat4 prod_xp = dot(sum_xp, sum_xp);
		vfloat4 prod_yp = dot(sum_yp, sum_yp);
		vfloat4 prod_zp = dot(sum_zp, sum_zp);

		vfloat4 best_vector = sum_xp;
		vfloat4 best_sum = prod_xp;

		vmask4 mask = prod_yp > best_sum;
		best_vector = select(best_vector, sum_yp, mask);
		best_sum = select(best_sum, prod_yp, mask);

		mask = prod_zp > best_sum;
		best_vector = select(best_vector, sum_zp, mask);

		pm[partition].dir = best_vector;
	}
}

/* See header for documentation. */
void compute_avgs_and_dirs_3_comp_rgb(
	const partition_info& pi,
	const image_block& blk,
	partition_metrics pm[BLOCK_MAX_PARTITIONS]
) {
	size_t partition_count = pi.partition_count;
	promise(partition_count > 0);

	// Pre-compute partition_averages
	vfloat4 partition_averages[BLOCK_MAX_PARTITIONS];
	compute_partition_averages_rgb(pi, blk, partition_averages);

	for (size_t partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *texel_indexes = pi.texels_of_partition[partition];
		size_t texel_count = pi.partition_texel_count[partition];
		promise(texel_count > 0);

		vfloat4 average = partition_averages[partition];
		pm[partition].avg = average;

		vfloat4 sum_xp = vfloat4::zero();
		vfloat4 sum_yp = vfloat4::zero();
		vfloat4 sum_zp = vfloat4::zero();

		for (size_t i = 0; i < texel_count; i++)
		{
			unsigned int iwt = texel_indexes[i];

			vfloat4 texel_datum = blk.texel3(iwt);
			texel_datum = texel_datum - average;

			vfloat4 zero = vfloat4::zero();

			vmask4 tdm0 = texel_datum.swz<0,0,0,0>() > zero;
			sum_xp += select(zero, texel_datum, tdm0);

			vmask4 tdm1 = texel_datum.swz<1,1,1,1>() > zero;
			sum_yp += select(zero, texel_datum, tdm1);

			vmask4 tdm2 = texel_datum.swz<2,2,2,2>() > zero;
			sum_zp += select(zero, texel_datum, tdm2);
		}

		vfloat4 prod_xp = dot(sum_xp, sum_xp);
		vfloat4 prod_yp = dot(sum_yp, sum_yp);
		vfloat4 prod_zp = dot(sum_zp, sum_zp);

		vfloat4 best_vector = sum_xp;
		vfloat4 best_sum = prod_xp;

		vmask4 mask = prod_yp > best_sum;
		best_vector = select(best_vector, sum_yp, mask);
		best_sum = select(best_sum, prod_yp, mask);

		mask = prod_zp > best_sum;
		best_vector = select(best_vector, sum_zp, mask);

		pm[partition].dir = best_vector;
	}
}

/* See header for documentation. */
void compute_avgs_and_dirs_2_comp(
	const partition_info& pt,
	const image_block& blk,
	unsigned int component1,
	unsigned int component2,
	partition_metrics pm[BLOCK_MAX_PARTITIONS]
) {
	vfloat4 average;

	const float* data_vr = nullptr;
	const float* data_vg = nullptr;

	if (component1 == 0 && component2 == 1)
	{
		average = blk.data_mean.swz<0, 1>();

		data_vr = blk.data_r;
		data_vg = blk.data_g;
	}
	else if (component1 == 0 && component2 == 2)
	{
		average = blk.data_mean.swz<0, 2>();

		data_vr = blk.data_r;
		data_vg = blk.data_b;
	}
	else  // (component1 == 1 && component2 == 2)
	{
		assert(component1 == 1 && component2 == 2);

		average = blk.data_mean.swz<1, 2>();

		data_vr = blk.data_g;
		data_vg = blk.data_b;
	}

	size_t partition_count = pt.partition_count;
	promise(partition_count > 0);

	for (size_t partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *texel_indexes = pt.texels_of_partition[partition];
		size_t texel_count = pt.partition_texel_count[partition];
		promise(texel_count > 0);

		// Only compute a partition mean if more than one partition
		if (partition_count > 1)
		{
			average = vfloat4::zero();
			for (size_t i = 0; i < texel_count; i++)
			{
				unsigned int iwt = texel_indexes[i];
				average += vfloat2(data_vr[iwt], data_vg[iwt]);
			}

			average = average / static_cast<float>(texel_count);
		}

		pm[partition].avg = average;

		vfloat4 sum_xp = vfloat4::zero();
		vfloat4 sum_yp = vfloat4::zero();

		for (size_t i = 0; i < texel_count; i++)
		{
			unsigned int iwt = texel_indexes[i];
			vfloat4 texel_datum = vfloat2(data_vr[iwt], data_vg[iwt]);
			texel_datum = texel_datum - average;

			vfloat4 zero = vfloat4::zero();

			vmask4 tdm0 = texel_datum.swz<0,0,0,0>() > zero;
			sum_xp += select(zero, texel_datum, tdm0);

			vmask4 tdm1 = texel_datum.swz<1,1,1,1>() > zero;
			sum_yp += select(zero, texel_datum, tdm1);
		}

		vfloat4 prod_xp = dot(sum_xp, sum_xp);
		vfloat4 prod_yp = dot(sum_yp, sum_yp);

		vfloat4 best_vector = sum_xp;
		vfloat4 best_sum = prod_xp;

		vmask4 mask = prod_yp > best_sum;
		best_vector = select(best_vector, sum_yp, mask);

		pm[partition].dir = best_vector;
	}
}

/* See header for documentation. */
void compute_error_squared_rgba(
	const partition_info& pi,
	const image_block& blk,
	const processed_line4 uncor_plines[BLOCK_MAX_PARTITIONS],
	const processed_line4 samec_plines[BLOCK_MAX_PARTITIONS],
	float line_lengths[BLOCK_MAX_PARTITIONS],
	float& uncor_error,
	float& samec_error
) {
	size_t partition_count = pi.partition_count;
	promise(partition_count > 0);

	vfloatacc uncor_errorsumv = vfloatacc::zero();
	vfloatacc samec_errorsumv = vfloatacc::zero();

	vgatherf_table tr = vgatherf_load(blk.data_r, blk.texel_count);
	vgatherf_table tg = vgatherf_load(blk.data_g, blk.texel_count);
	vgatherf_table tb = vgatherf_load(blk.data_b, blk.texel_count);
	vgatherf_table ta = vgatherf_load(blk.data_a, blk.texel_count);

	for (size_t partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *texel_indexes = pi.texels_of_partition[partition];

		processed_line4 l_uncor = uncor_plines[partition];
		processed_line4 l_samec = samec_plines[partition];

		size_t texel_count = pi.partition_texel_count[partition];
		promise(texel_count > 0);

		// Vectorize some useful scalar inputs
		vfloat l_uncor_bs0(l_uncor.bs.lane<0>());
		vfloat l_uncor_bs1(l_uncor.bs.lane<1>());
		vfloat l_uncor_bs2(l_uncor.bs.lane<2>());
		vfloat l_uncor_bs3(l_uncor.bs.lane<3>());

		vfloat l_uncor_amod0(l_uncor.amod.lane<0>());
		vfloat l_uncor_amod1(l_uncor.amod.lane<1>());
		vfloat l_uncor_amod2(l_uncor.amod.lane<2>());
		vfloat l_uncor_amod3(l_uncor.amod.lane<3>());

		vfloat l_samec_bs0(l_samec.bs.lane<0>());
		vfloat l_samec_bs1(l_samec.bs.lane<1>());
		vfloat l_samec_bs2(l_samec.bs.lane<2>());
		vfloat l_samec_bs3(l_samec.bs.lane<3>());

		assert(all(l_samec.amod == vfloat4(0.0f)));

		vfloat uncor_loparamv(1e10f);
		vfloat uncor_hiparamv(-1e10f);

		vfloat ew_r(blk.channel_weight.lane<0>());
		vfloat ew_g(blk.channel_weight.lane<1>());
		vfloat ew_b(blk.channel_weight.lane<2>());
		vfloat ew_a(blk.channel_weight.lane<3>());

		// Padded texel indexes repeat the last value, preserving min/max. Mask padded lanes when
		// computing line weights.
		vint lane_ids = vint::lane_id();
		for (size_t i = 0; i < texel_count; i += ASTCENC_SIMD_WIDTH)
		{
			vmask mask = lane_ids < vint_from_size(texel_count);
			const uint8_t* texel_idxs = texel_indexes + i;

			vfloat data_r = gatherf(tr, texel_idxs);
			vfloat data_g = gatherf(tg, texel_idxs);
			vfloat data_b = gatherf(tb, texel_idxs);
			vfloat data_a = gatherf(ta, texel_idxs);

			vfloat uncor_param = (data_r * l_uncor_bs0)
			                   + (data_g * l_uncor_bs1)
			                   + (data_b * l_uncor_bs2)
			                   + (data_a * l_uncor_bs3);

			uncor_loparamv = min(uncor_param, uncor_loparamv);
			uncor_hiparamv = max(uncor_param, uncor_hiparamv);

			vfloat uncor_dist0 = (l_uncor_amod0 - data_r)
			                   + (uncor_param * l_uncor_bs0);
			vfloat uncor_dist1 = (l_uncor_amod1 - data_g)
			                   + (uncor_param * l_uncor_bs1);
			vfloat uncor_dist2 = (l_uncor_amod2 - data_b)
			                   + (uncor_param * l_uncor_bs2);
			vfloat uncor_dist3 = (l_uncor_amod3 - data_a)
			                   + (uncor_param * l_uncor_bs3);

			vfloat uncor_err = (ew_r * uncor_dist0 * uncor_dist0)
			                 + (ew_g * uncor_dist1 * uncor_dist1)
			                 + (ew_b * uncor_dist2 * uncor_dist2)
			                 + (ew_a * uncor_dist3 * uncor_dist3);

			haccumulate(uncor_errorsumv, uncor_err, mask);

			// Process samechroma data
			vfloat samec_param = (data_r * l_samec_bs0)
			                   + (data_g * l_samec_bs1)
			                   + (data_b * l_samec_bs2)
			                   + (data_a * l_samec_bs3);

			vfloat samec_dist0 = samec_param * l_samec_bs0 - data_r;
			vfloat samec_dist1 = samec_param * l_samec_bs1 - data_g;
			vfloat samec_dist2 = samec_param * l_samec_bs2 - data_b;
			vfloat samec_dist3 = samec_param * l_samec_bs3 - data_a;

			vfloat samec_err = (ew_r * samec_dist0 * samec_dist0)
			                 + (ew_g * samec_dist1 * samec_dist1)
			                 + (ew_b * samec_dist2 * samec_dist2)
			                 + (ew_a * samec_dist3 * samec_dist3);

			haccumulate(samec_errorsumv, samec_err, mask);

			lane_ids += vint(ASTCENC_SIMD_WIDTH);
		}

		// Turn very small numbers and NaNs into a small number
		float uncor_linelen = hmax_s(uncor_hiparamv) - hmin_s(uncor_loparamv);
		line_lengths[partition] = astc::max(uncor_linelen, 1e-7f);
	}

	uncor_error = hadd_s(uncor_errorsumv);
	samec_error = hadd_s(samec_errorsumv);
}

/* See header for documentation. */
void compute_error_squared_rgba_grayscale(
	const partition_info& pi,
	const image_block& blk,
	const processed_line4 uncor_plines[BLOCK_MAX_PARTITIONS],
	const processed_line4 samec_plines[BLOCK_MAX_PARTITIONS],
	float line_lengths[BLOCK_MAX_PARTITIONS],
	float& uncor_error,
	float& samec_error
) {
	assert(blk.grayscale);
	size_t partition_count = pi.partition_count;
	promise(partition_count > 0);

	vfloatacc uncor_errorsumv = vfloatacc::zero();
	vfloatacc samec_errorsumv = vfloatacc::zero();

	vgatherf_table tr = vgatherf_load(blk.data_r, blk.texel_count);
	vgatherf_table ta = vgatherf_load(blk.data_a, blk.texel_count);

	for (size_t partition = 0; partition < partition_count; partition++)
	{
		const uint8_t *texel_indexes = pi.texels_of_partition[partition];

		processed_line4 l_uncor = uncor_plines[partition];
		processed_line4 l_samec = samec_plines[partition];

		assert(all(l_uncor.bs == l_uncor.bs.swz<0, 0, 0, 3>()));
		assert(all(l_uncor.amod == l_uncor.amod.swz<0, 0, 0, 3>()));
		assert(all(l_samec.bs == l_samec.bs.swz<0, 0, 0, 3>()));
		assert(all(l_samec.amod == vfloat4(0.0f)));

		size_t texel_count = pi.partition_texel_count[partition];
		promise(texel_count > 0);

		// The line coefficients have identical RGB components
		vfloat l_uncor_bs0(l_uncor.bs.lane<0>());
		vfloat l_uncor_bs3(l_uncor.bs.lane<3>());
		vfloat l_uncor_amod0(l_uncor.amod.lane<0>());
		vfloat l_uncor_amod3(l_uncor.amod.lane<3>());
		vfloat l_samec_bs0(l_samec.bs.lane<0>());
		vfloat l_samec_bs3(l_samec.bs.lane<3>());

		vfloat uncor_loparamv(1e10f);
		vfloat uncor_hiparamv(-1e10f);

		vfloat ew_r(blk.channel_weight.lane<0>());
		vfloat ew_g(blk.channel_weight.lane<1>());
		vfloat ew_b(blk.channel_weight.lane<2>());
		vfloat ew_a(blk.channel_weight.lane<3>());

		// Padded texel indexes repeat the last value, preserving min/max. Mask these lanes when
		// accumulating errors.
		vint lane_ids = vint::lane_id();
		for (size_t i = 0; i < texel_count; i += ASTCENC_SIMD_WIDTH)
		{
			vmask mask = lane_ids < vint_from_size(texel_count);
			const uint8_t* texel_idxs = texel_indexes + i;

			vfloat data_r = gatherf(tr, texel_idxs);
			vfloat data_a = gatherf(ta, texel_idxs);

			// Reuse the RGB product, retaining all three additions in order
			vfloat uncor_rgb_param = data_r * l_uncor_bs0;
			vfloat uncor_param = uncor_rgb_param
			                   + uncor_rgb_param
			                   + uncor_rgb_param
			                   + (data_a * l_uncor_bs3);

			uncor_loparamv = min(uncor_param, uncor_loparamv);
			uncor_hiparamv = max(uncor_param, uncor_hiparamv);

			vfloat uncor_dist0 = (l_uncor_amod0 - data_r)
			                   + (uncor_param * l_uncor_bs0);
			vfloat uncor_dist3 = (l_uncor_amod3 - data_a)
			                   + (uncor_param * l_uncor_bs3);

			// Channel weights may differ, even though the RGB residuals match
			vfloat uncor_err = (ew_r * uncor_dist0 * uncor_dist0)
			                 + (ew_g * uncor_dist0 * uncor_dist0)
			                 + (ew_b * uncor_dist0 * uncor_dist0)
			                 + (ew_a * uncor_dist3 * uncor_dist3);

			haccumulate(uncor_errorsumv, uncor_err, mask);

			vfloat samec_rgb_param = data_r * l_samec_bs0;
			vfloat samec_param = samec_rgb_param
			                   + samec_rgb_param
			                   + samec_rgb_param
			                   + (data_a * l_samec_bs3);

			vfloat samec_dist0 = samec_param * l_samec_bs0 - data_r;
			vfloat samec_dist3 = samec_param * l_samec_bs3 - data_a;

			vfloat samec_err = (ew_r * samec_dist0 * samec_dist0)
			                 + (ew_g * samec_dist0 * samec_dist0)
			                 + (ew_b * samec_dist0 * samec_dist0)
			                 + (ew_a * samec_dist3 * samec_dist3);

			haccumulate(samec_errorsumv, samec_err, mask);
			lane_ids += vint(ASTCENC_SIMD_WIDTH);
		}

		float uncor_linelen = hmax_s(uncor_hiparamv) - hmin_s(uncor_loparamv);
		line_lengths[partition] = astc::max(uncor_linelen, 1e-7f);
	}

	uncor_error = hadd_s(uncor_errorsumv);
	samec_error = hadd_s(samec_errorsumv);
}

/* See header for documentation. */
void compute_error_squared_rgb(
	const partition_info& pi,
	const image_block& blk,
	partition_lines3 plines[BLOCK_MAX_PARTITIONS],
	float& uncor_error,
	float& samec_error
) {
	size_t partition_count = pi.partition_count;
	promise(partition_count > 0);

	vfloatacc uncor_errorsumv = vfloatacc::zero();
	vfloatacc samec_errorsumv = vfloatacc::zero();

	vgatherf_table tr = vgatherf_load(blk.data_r, blk.texel_count);
	vgatherf_table tg = vgatherf_load(blk.data_g, blk.texel_count);
	vgatherf_table tb = vgatherf_load(blk.data_b, blk.texel_count);

	for (size_t partition = 0; partition < partition_count; partition++)
	{
		partition_lines3& pl = plines[partition];
		const uint8_t *texel_indexes = pi.texels_of_partition[partition];
		size_t texel_count = pi.partition_texel_count[partition];
		promise(texel_count > 0);

		processed_line3 l_uncor = pl.uncor_pline;
		processed_line3 l_samec = pl.samec_pline;

		// Vectorize some useful scalar inputs
		vfloat l_uncor_bs0(l_uncor.bs.lane<0>());
		vfloat l_uncor_bs1(l_uncor.bs.lane<1>());
		vfloat l_uncor_bs2(l_uncor.bs.lane<2>());

		vfloat l_uncor_amod0(l_uncor.amod.lane<0>());
		vfloat l_uncor_amod1(l_uncor.amod.lane<1>());
		vfloat l_uncor_amod2(l_uncor.amod.lane<2>());

		vfloat l_samec_bs0(l_samec.bs.lane<0>());
		vfloat l_samec_bs1(l_samec.bs.lane<1>());
		vfloat l_samec_bs2(l_samec.bs.lane<2>());

		assert(all(l_samec.amod == vfloat4(0.0f)));

		vfloat uncor_loparamv(1e10f);
		vfloat uncor_hiparamv(-1e10f);

		vfloat ew_r(blk.channel_weight.lane<0>());
		vfloat ew_g(blk.channel_weight.lane<1>());
		vfloat ew_b(blk.channel_weight.lane<2>());

		// Padded weights repeat the last value, preserving min/max. Mask padded lanes when
		// computing line weights.
		vint lane_ids = vint::lane_id();
		for (size_t i = 0; i < texel_count; i += ASTCENC_SIMD_WIDTH)
		{
			vmask mask = lane_ids < vint_from_size(texel_count);
			const uint8_t* texel_idxs = texel_indexes + i;

			vfloat data_r = gatherf(tr, texel_idxs);
			vfloat data_g = gatherf(tg, texel_idxs);
			vfloat data_b = gatherf(tb, texel_idxs);

			vfloat uncor_param = (data_r * l_uncor_bs0)
			                   + (data_g * l_uncor_bs1)
			                   + (data_b * l_uncor_bs2);

			uncor_loparamv = min(uncor_param, uncor_loparamv);
			uncor_hiparamv = max(uncor_param, uncor_hiparamv);

			vfloat uncor_dist0 = (l_uncor_amod0 - data_r)
			                   + (uncor_param * l_uncor_bs0);
			vfloat uncor_dist1 = (l_uncor_amod1 - data_g)
			                   + (uncor_param * l_uncor_bs1);
			vfloat uncor_dist2 = (l_uncor_amod2 - data_b)
			                   + (uncor_param * l_uncor_bs2);

			vfloat uncor_err = (ew_r * uncor_dist0 * uncor_dist0)
			                 + (ew_g * uncor_dist1 * uncor_dist1)
			                 + (ew_b * uncor_dist2 * uncor_dist2);

			haccumulate(uncor_errorsumv, uncor_err, mask);

			// Process samechroma data
			vfloat samec_param = (data_r * l_samec_bs0)
			                   + (data_g * l_samec_bs1)
			                   + (data_b * l_samec_bs2);

			vfloat samec_dist0 = samec_param * l_samec_bs0 - data_r;
			vfloat samec_dist1 = samec_param * l_samec_bs1 - data_g;
			vfloat samec_dist2 = samec_param * l_samec_bs2 - data_b;

			vfloat samec_err = (ew_r * samec_dist0 * samec_dist0)
			                 + (ew_g * samec_dist1 * samec_dist1)
			                 + (ew_b * samec_dist2 * samec_dist2);

			haccumulate(samec_errorsumv, samec_err, mask);

			lane_ids += vint(ASTCENC_SIMD_WIDTH);
		}

		// Turn very small numbers and NaNs into a small number
		float uncor_linelen = hmax_s(uncor_hiparamv) - hmin_s(uncor_loparamv);
		pl.line_length = astc::max(uncor_linelen, 1e-7f);
	}

	uncor_error = hadd_s(uncor_errorsumv);
	samec_error = hadd_s(samec_errorsumv);
}

#if defined(ASTCENC_HAS_NATIVE_QUAD_OPS)

/* See header for documentation. */
ASTCENC_NO_INLINE void compute_partition_errors_rgb_batch(
	const partition_info* const* candidates,
	const image_block& blk,
	float weight_imprecision,
	float* uncor_errors,
	float* samec_errors
) {
	// This implementation is a batch-based implementation that merges the functionality of
	// compute_avgs_and_dirs_3_comp_rgb, compute_error_squared_rgb, and additional scoring logic.
	// It is designed to exploit wider vector hardware while keeping floating-point invariance with
	// the original implementation.
	//
	// The original code for direction finding processes a single candidate at a time, operating on
	// vectors of RGB0 (alpha unused) texel data. This path only uses 3 lanes out of each SIMD
	// vector, so gets no benefit from wider vectors.
	//
	// The original code for error calculation processes a single candidate at a time, computing
	// based on SIMD-width vectors of each color channel within a single partition. It uses linear
	// quad-by-quad accumulation at the end of the loop to ensure invariance across vector lengths.
	// This approach basically tries to use the vector length to brute-force process a single
	// candidate faster, but results in very under-filled vectors. Individual color partitions
	// contain relatively few texels, especially for smaller block sizes.
	//
	// This implementation attacks the problem differently. It fills wider vectors by processing
	// multiple partitioning candidates in parallel. Each "quad" of four lanes contains data
	// for a single partitioning candidate. As we now only group partitions into chunks of 4,
	// they are less likely to be idle, and the longer SIMD vector lengths get much higher average
	// occupancy. To vectorize this without spilling a lot of state between phases, we fuse the
	// direction finding and the error calculation within a single outer loop.
	//
	// Preserve floating-point invariance by keeping each candidate's arithmetic order unchanged.
	// Partition averages still use the original helper, and direction sums retain texel order.
	// During error scoring, each candidate has four running sums, matching the original
	// four-lane accumulator. Successive texel quads update these sums lane-wise in partition
	// order. The sums persist across partition boundaries and are reduced only after all
	// partitions, using (s0 + s2) + (s1 + s3), matching the original hadd_s() reduction. Weight
	// quantization errors are then added in the original partition order.

	constexpr unsigned int batch_size = ASTCENC_SIMD_WIDTH / 4;
	unsigned int partition_count = candidates[0]->partition_count;

	// Reuse the original per-channel texel accumulation and reduction for partition averages
	vfloat4 averages[batch_size][BLOCK_MAX_PARTITIONS];
	vfloat4 units[batch_size];
	for (unsigned int c = 0; c < batch_size; c++)
	{
		compute_partition_averages_rgb(*candidates[c], blk, averages[c]);
		units[c] = unit3();
	}

	vfloat safe = pack_quads(units);
	vfloat zero = vfloat::zero();

	// Preserve four running texel-error sums per candidate across partition boundaries
	vfloat uncor_sum = zero;
	vfloat samec_sum = zero;

	vfloat uncor_quant[BLOCK_MAX_PARTITIONS];
	vfloat samec_quant[BLOCK_MAX_PARTITIONS];

	vfloat ew_r(blk.channel_weight.lane<0>());
	vfloat ew_g(blk.channel_weight.lane<1>());
	vfloat ew_b(blk.channel_weight.lane<2>());

	vgatherf_table tr = vgatherf_load(blk.data_r, blk.texel_count);
	vgatherf_table tg = vgatherf_load(blk.data_g, blk.texel_count);
	vgatherf_table tb = vgatherf_load(blk.data_b, blk.texel_count);

	// Process partition p from every candidate together, preserving each candidate's
	// partition order
	for (unsigned int p = 0; p < partition_count; p++)
	{
		vfloat4 colors[batch_size];
		vfloat4 counts[batch_size];
		unsigned int max_count = 0;
		for (unsigned int c = 0; c < batch_size; c++)
		{
			unsigned int count = candidates[c]->partition_texel_count[p];
			max_count = astc::max(max_count, count);
			colors[c] = averages[c][p];
			counts[c] = vfloat4(static_cast<float>(count));
		}

		// Direction estimation uses RGB0 lanes per candidate, e.g. candidates 0 and 1 on AVX2:
		// [c0.R, c0.G, c0.B, 0 | c1.R, c1.G, c1.B, 0]
		vfloat mean = pack_quads(colors);
		vfloat countv = pack_quads(counts);
		vfloat sum_xp = zero;
		vfloat sum_yp = zero;
		vfloat sum_zp = zero;
		vfloat index = zero;

		for (unsigned int i = 0; i < max_count; i++)
		{
			for (unsigned int c = 0; c < batch_size; c++)
			{
				unsigned int last = candidates[c]->partition_texel_count[p] - 1u;
				unsigned int texel = candidates[c]->texels_of_partition[p][astc::min(i, last)];
				colors[c] = blk.texel3(texel);
			}

			// Preserve texel order and mask exhausted candidates without rebuilding RGB vectors
			vfloat datum = select(zero, pack_quads(colors) - mean, index < countv);
			index += vfloat(1.0f);

			sum_xp += select(zero, datum, swizzle_quads<0, 0, 0, 0>(datum) > zero);
			sum_yp += select(zero, datum, swizzle_quads<1, 1, 1, 1>(datum) > zero);
			sum_zp += select(zero, datum, swizzle_quads<2, 2, 2, 2>(datum) > zero);
		}

		// Reduce squared RGB0 components within each candidate
		vfloat prod_xp = dot_quads(sum_xp, sum_xp);
		vfloat prod_yp = dot_quads(sum_yp, sum_yp);
		vfloat prod_zp = dot_quads(sum_zp, sum_zp);

		vfloat best = sum_xp;
		vfloat best_sum = prod_xp;

		vmask mask = prod_yp > best_sum;
		best = select(best, sum_yp, mask);
		best_sum = select(best_sum, prod_yp, mask);

		mask = prod_zp > best_sum;
		best = select(best, sum_zp, mask);

		// Keep direction selection, normalization, and line setup in the same candidate groups
		vfloat bs = normalize_safe_quads(best, safe);
		vfloat samec_bs = normalize_safe_quads(mean, safe);
		vfloat amod = mean - bs * sum_rgb_quads(mean * bs);

		// Scoring uses four texels per candidate quad, with separate R, G, and B vectors
		// Broadcast each line component to all four texels of its candidate
		vfloat bs0 = swizzle_quads<0, 0, 0, 0>(bs);
		vfloat bs1 = swizzle_quads<1, 1, 1, 1>(bs);
		vfloat bs2 = swizzle_quads<2, 2, 2, 2>(bs);

		vfloat amod0 = swizzle_quads<0, 0, 0, 0>(amod);
		vfloat amod1 = swizzle_quads<1, 1, 1, 1>(amod);
		vfloat amod2 = swizzle_quads<2, 2, 2, 2>(amod);

		vfloat samec0 = swizzle_quads<0, 0, 0, 0>(samec_bs);
		vfloat samec1 = swizzle_quads<1, 1, 1, 1>(samec_bs);
		vfloat samec2 = swizzle_quads<2, 2, 2, 2>(samec_bs);
		vfloat lane_ids = int_to_float(vint::lane_id() & vint(3));

		vfloat low(1e10f);
		vfloat high(-1e10f);

		for (unsigned int i = 0; i < max_count; i += 4)
		{
			vint4 indices[batch_size];

			for (unsigned int c = 0; c < batch_size; c++)
			{
				// Repeat the final padded quad for exhausted candidates to preserve min/max
				unsigned int last = (candidates[c]->partition_texel_count[p] - 1u) & ~3u;
				indices[c] = vint4(candidates[c]->texels_of_partition[p] + astc::min(i, last));
			}

			vint texels = pack_quads(indices);
			vmask active = lane_ids < countv;
			lane_ids += vfloat(4.0f);

			vfloat r = gatherf(tr, texels);
			vfloat g = gatherf(tg, texels);
			vfloat b = gatherf(tb, texels);

			vfloat param = (r * bs0) + (g * bs1) + (b * bs2);

			low = min(param, low);
			high = max(param, high);

			vfloat d0 = (amod0 - r) + (param * bs0);
			vfloat d1 = (amod1 - g) + (param * bs1);
			vfloat d2 = (amod2 - b) + (param * bs2);

			vfloat error = (ew_r * d0 * d0) + (ew_g * d1 * d1) + (ew_b * d2 * d2);

			// Match haccumulate() by adding successive texel quads lane-wise, without reduction
			uncor_sum += select(zero, error, active);

			param = (r * samec0) + (g * samec1) + (b * samec2);

			d0 = param * samec0 - r;
			d1 = param * samec1 - g;
			d2 = param * samec2 - b;

			error = (ew_r * d0 * d0) + (ew_g * d1 * d1) + (ew_b * d2 * d2);

			samec_sum += select(zero, error, active);
		}

		// Defer adding quantization errors until after the texel sums are reduced
		vfloat length = max(max_quads(high) - min_quads(low), 1e-7f);
		vfloat error_weights = countv * vfloat(weight_imprecision);

		vfloat uncor_vector = bs * length;
		uncor_quant[p] = sum_rgb_quads((uncor_vector * uncor_vector) * error_weights);

		vfloat samec_vector = samec_bs * length;
		samec_quant[p] = sum_rgb_quads((samec_vector * samec_vector) * error_weights);
	}

	// Reduce each candidate's texel sums as (s0 + s2) + (s1 + s3), matching hadd_s(vfloat4)
	// Then add quantization errors in the original partition order
	uncor_sum = sum_quads(uncor_sum);
	samec_sum = sum_quads(samec_sum);
	for (unsigned int p = 0; p < partition_count; p++)
	{
		uncor_sum += uncor_quant[p];
		samec_sum += samec_quant[p];
	}

	store_lane0_quads(uncor_sum, uncor_errors);
	store_lane0_quads(samec_sum, samec_errors);
}
#endif

#endif
