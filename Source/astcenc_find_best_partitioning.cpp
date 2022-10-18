// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2011-2022 Arm Limited
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
 * The partition search operates in two stages. The first pass uses kmeans clustering to group
 * texels into an ideal partitioning for the requested partition count, and then compares that
 * against the 1024 partitionings generated by the ASTC partition hash function. The generated
 * partitions are then ranked by the number of texels in the wrong partition, compared to the ideal
 * clustering. All 1024 partitions are tested for similarity and ranked, apart from duplicates and
 * partitionings that actually generate fewer than the requested partition count, but only the top
 * N candidates are actually put through a more detailed search. N is determined by the compressor
 * quality preset.
 *
 * For the detailed search, each candidate is checked against two possible encoding methods:
 *
 *   - The best partitioning assuming different chroma colors (RGB + RGB or RGB + delta endpoints).
 *   - The best partitioning assuming same chroma colors (RGB + scale endpoints).
 *
 * This is implemented by computing the compute mean color and dominant direction for each
 * partition. This defines two lines, both of which go through the mean color value.
 *
 * - One line has a direction defined by the dominant direction; this is used to assess the error
 *   from using an uncorrelated color representation.
 * - The other line goes through (0,0,0,1) and is used to assess the error from using a same chroma
 *   (RGB + scale) color representation.
 *
 * The best candidate is selected by computing the squared-errors that result from using these
 * lines for endpoint selection.
 */

#include <limits>
#include "astcenc_internal.h"

/**
 * @brief Pick some initial kmeans cluster centers.
 *
 * @param      blk               The image block color data to compress.
 * @param      texel_count       The number of texels in the block.
 * @param      partition_count   The number of partitions in the block.
 * @param[out] cluster_centers   The initial partition cluster center colors.
 */
static void kmeans_init(
	const image_block& blk,
	unsigned int texel_count,
	unsigned int partition_count,
	vfloat4 cluster_centers[BLOCK_MAX_PARTITIONS]
) {
	promise(texel_count > 0);
	promise(partition_count > 0);

	unsigned int clusters_selected = 0;
	float distances[BLOCK_MAX_TEXELS];

	// Pick a random sample as first cluster center; 145897 from random.org
	unsigned int sample = 145897 % texel_count;
	vfloat4 center_color = blk.texel(sample);
	cluster_centers[clusters_selected] = center_color;
	clusters_selected++;

	// Compute the distance to the first cluster center
	float distance_sum = 0.0f;
	for (unsigned int i = 0; i < texel_count; i++)
	{
		vfloat4 color = blk.texel(i);
		vfloat4 diff = color - center_color;
		float distance = dot_s(diff * diff, blk.channel_weight);
		distance_sum += distance;
		distances[i] = distance;
	}

	// More numbers from random.org for weighted-random center selection
	const float cluster_cutoffs[9] {
		0.626220f, 0.932770f, 0.275454f,
		0.318558f, 0.240113f, 0.009190f,
		0.347661f, 0.731960f, 0.156391f
	};

	unsigned int cutoff = (clusters_selected - 1) + 3 * (partition_count - 2);

	// Pick the remaining samples as needed
	while (true)
	{
		// Pick the next center in a weighted-random fashion.
		float summa = 0.0f;
		float distance_cutoff = distance_sum * cluster_cutoffs[cutoff++];
		for (sample = 0; sample < texel_count; sample++)
		{
			summa += distances[sample];
			if (summa >= distance_cutoff)
			{
				break;
			}
		}

		// Clamp to a valid range and store the selected cluster center
		sample = astc::min(sample, texel_count - 1);

		center_color = blk.texel(sample);
		cluster_centers[clusters_selected++] = center_color;
		if (clusters_selected >= partition_count)
		{
			break;
		}

		// Compute the distance to the new cluster center, keep the min dist
		distance_sum = 0.0f;
		for (unsigned int i = 0; i < texel_count; i++)
		{
			vfloat4 color = blk.texel(i);
			vfloat4 diff = color - center_color;
			float distance = dot_s(diff * diff, blk.channel_weight);
			distance = astc::min(distance, distances[i]);
			distance_sum += distance;
			distances[i] = distance;
		}
	}
}

/**
 * @brief Assign texels to clusters, based on a set of chosen center points.
 *
 * @param      blk                  The image block color data to compress.
 * @param      texel_count          The number of texels in the block.
 * @param      partition_count      The number of partitions in the block.
 * @param      cluster_centers      The partition cluster center colors.
 * @param[out] partition_of_texel   The partition assigned for each texel.
 */
static void kmeans_assign(
	const image_block& blk,
	unsigned int texel_count,
	unsigned int partition_count,
	const vfloat4 cluster_centers[BLOCK_MAX_PARTITIONS],
	uint8_t partition_of_texel[BLOCK_MAX_TEXELS]
) {
	promise(texel_count > 0);
	promise(partition_count > 0);

	uint8_t partition_texel_count[BLOCK_MAX_PARTITIONS] { 0 };

	// Find the best partition for every texel
	for (unsigned int i = 0; i < texel_count; i++)
	{
		float best_distance = std::numeric_limits<float>::max();
		unsigned int best_partition = 0;

		vfloat4 color = blk.texel(i);
		for (unsigned int j = 0; j < partition_count; j++)
		{
			vfloat4 diff = color - cluster_centers[j];
			float distance = dot_s(diff * diff, blk.channel_weight);
			if (distance < best_distance)
			{
				best_distance = distance;
				best_partition = j;
			}
		}

		partition_of_texel[i] = static_cast<uint8_t>(best_partition);
		partition_texel_count[best_partition]++;
	}

	// It is possible to get a situation where a partition ends up without any texels. In this case,
	// assign texel N to partition N. This is silly, but ensures that every partition retains at
	// least one texel. Reassigning a texel in this manner may cause another partition to go empty,
	// so if we actually did a reassignment, run the whole loop over again.
	bool problem_case;
	do
	{
		problem_case = false;
		for (unsigned int i = 0; i < partition_count; i++)
		{
			if (partition_texel_count[i] == 0)
			{
				partition_texel_count[partition_of_texel[i]]--;
				partition_texel_count[i]++;
				partition_of_texel[i] = static_cast<uint8_t>(i);
				problem_case = true;
			}
		}
	} while (problem_case);
}

/**
 * @brief Compute new cluster centers based on their center of gravity.
 *
 * @param       blk                  The image block color data to compress.
 * @param       texel_count          The number of texels in the block.
 * @param       partition_count      The number of partitions in the block.
 * @param[out]  cluster_centers      The new cluster center colors.
 * @param       partition_of_texel   The partition assigned for each texel.
 */
static void kmeans_update(
	const image_block& blk,
	unsigned int texel_count,
	unsigned int partition_count,
	vfloat4 cluster_centers[BLOCK_MAX_PARTITIONS],
	const uint8_t partition_of_texel[BLOCK_MAX_TEXELS]
) {
	promise(texel_count > 0);
	promise(partition_count > 0);

	vfloat4 color_sum[BLOCK_MAX_PARTITIONS] {
		vfloat4::zero(),
		vfloat4::zero(),
		vfloat4::zero(),
		vfloat4::zero()
	};

	uint8_t partition_texel_count[BLOCK_MAX_PARTITIONS] { 0 };

	// Find the center-of-gravity in each cluster
	for (unsigned int i = 0; i < texel_count; i++)
	{
		uint8_t partition = partition_of_texel[i];
		color_sum[partition] += blk.texel(i);
		partition_texel_count[partition]++;
	}

	// Set the center of gravity to be the new cluster center
	for (unsigned int i = 0; i < partition_count; i++)
	{
		float scale = 1.0f / static_cast<float>(partition_texel_count[i]);
		cluster_centers[i] = color_sum[i] * scale;
	}
}

/**
 * @brief Compute bit-mismatch for partitioning in 2-partition mode.
 *
 * @param a   The texel assignment bitvector for the block.
 * @param b   The texel assignment bitvector for the partition table.
 *
 * @return    The number of bit mismatches.
 */
static inline unsigned int partition_mismatch2(
	const uint64_t a[2],
	const uint64_t b[2]
) {
	int v1 = popcount(a[0] ^ b[0]) + popcount(a[1] ^ b[1]);
	int v2 = popcount(a[0] ^ b[1]) + popcount(a[1] ^ b[0]);
	return astc::min(v1, v2);
}

/**
 * @brief Compute bit-mismatch for partitioning in 3-partition mode.
 *
 * @param a   The texel assignment bitvector for the block.
 * @param b   The texel assignment bitvector for the partition table.
 *
 * @return    The number of bit mismatches.
 */
static inline unsigned int partition_mismatch3(
	const uint64_t a[3],
	const uint64_t b[3]
) {
	int p00 = popcount(a[0] ^ b[0]);
	int p01 = popcount(a[0] ^ b[1]);
	int p02 = popcount(a[0] ^ b[2]);

	int p10 = popcount(a[1] ^ b[0]);
	int p11 = popcount(a[1] ^ b[1]);
	int p12 = popcount(a[1] ^ b[2]);

	int p20 = popcount(a[2] ^ b[0]);
	int p21 = popcount(a[2] ^ b[1]);
	int p22 = popcount(a[2] ^ b[2]);

	int s0 = p11 + p22;
	int s1 = p12 + p21;
	int v0 = astc::min(s0, s1) + p00;

	int s2 = p10 + p22;
	int s3 = p12 + p20;
	int v1 = astc::min(s2, s3) + p01;

	int s4 = p10 + p21;
	int s5 = p11 + p20;
	int v2 = astc::min(s4, s5) + p02;

	return astc::min(v0, v1, v2);
}

/**
 * @brief Compute bit-mismatch for partitioning in 4-partition mode.
 *
 * @param a   The texel assignment bitvector for the block.
 * @param b   The texel assignment bitvector for the partition table.
 *
 * @return    The number of bit mismatches.
 */
static inline unsigned int partition_mismatch4(
	const uint64_t a[4],
	const uint64_t b[4]
) {
	int p00 = popcount(a[0] ^ b[0]);
	int p01 = popcount(a[0] ^ b[1]);
	int p02 = popcount(a[0] ^ b[2]);
	int p03 = popcount(a[0] ^ b[3]);

	int p10 = popcount(a[1] ^ b[0]);
	int p11 = popcount(a[1] ^ b[1]);
	int p12 = popcount(a[1] ^ b[2]);
	int p13 = popcount(a[1] ^ b[3]);

	int p20 = popcount(a[2] ^ b[0]);
	int p21 = popcount(a[2] ^ b[1]);
	int p22 = popcount(a[2] ^ b[2]);
	int p23 = popcount(a[2] ^ b[3]);

	int p30 = popcount(a[3] ^ b[0]);
	int p31 = popcount(a[3] ^ b[1]);
	int p32 = popcount(a[3] ^ b[2]);
	int p33 = popcount(a[3] ^ b[3]);

	int mx23 = astc::min(p22 + p33, p23 + p32);
	int mx13 = astc::min(p21 + p33, p23 + p31);
	int mx12 = astc::min(p21 + p32, p22 + p31);
	int mx03 = astc::min(p20 + p33, p23 + p30);
	int mx02 = astc::min(p20 + p32, p22 + p30);
	int mx01 = astc::min(p21 + p30, p20 + p31);

	int v0 = p00 + astc::min(p11 + mx23, p12 + mx13, p13 + mx12);
	int v1 = p01 + astc::min(p10 + mx23, p12 + mx03, p13 + mx02);
	int v2 = p02 + astc::min(p11 + mx03, p10 + mx13, p13 + mx01);
	int v3 = p03 + astc::min(p11 + mx02, p12 + mx01, p10 + mx12);

	return astc::min(v0, v1, v2, v3);
}

using mismatch_dispatch = unsigned int (*)(const uint64_t*, const uint64_t*);

/**
 * @brief Count the partition table mismatches vs the data clustering.
 *
 * @param      bsd               The block size information.
 * @param      partition_count   The number of partitions in the block.
 * @param      bitmaps           The block texel partition assignment patterns.
 * @param[out] mismatch_counts   The array storing per partitioning mismatch counts.
 */
static void count_partition_mismatch_bits(
	const block_size_descriptor& bsd,
	unsigned int partition_count,
	const uint64_t bitmaps[BLOCK_MAX_PARTITIONS],
	unsigned int mismatch_counts[BLOCK_MAX_PARTITIONINGS]
) {
	unsigned int active_count = bsd.partitioning_count_selected[partition_count - 1];

	if (partition_count == 2)
	{
		for (unsigned int i = 0; i < active_count; i++)
		{
			mismatch_counts[i] = partition_mismatch2(bitmaps, bsd.coverage_bitmaps_2[i]);
		}
	}
	else if (partition_count == 3)
	{
		for (unsigned int i = 0; i < active_count; i++)
		{
			mismatch_counts[i] = partition_mismatch3(bitmaps, bsd.coverage_bitmaps_3[i]);
		}
	}
	else
	{
		for (unsigned int i = 0; i < active_count; i++)
		{
			mismatch_counts[i] = partition_mismatch4(bitmaps, bsd.coverage_bitmaps_4[i]);
		}
	}
}

/**
 * @brief Use counting sort on the mismatch array to sort partition candidates.
 *
 * @param      partitioning_count   The number of packed partitionings.
 * @param      mismatch_count       Partitioning mismatch counts, in index order.
 * @param[out] partition_ordering   Partition index values, in mismatch order.
 *
 * @return The number of active partitions in this selection.
 */
static unsigned int get_partition_ordering_by_mismatch_bits(
	unsigned int partitioning_count,
	const unsigned int mismatch_count[BLOCK_MAX_PARTITIONINGS],
	unsigned int partition_ordering[BLOCK_MAX_PARTITIONINGS]
) {
	unsigned int mscount[256] { 0 };

	// Create the histogram of mismatch counts
	for (unsigned int i = 0; i < partitioning_count; i++)
	{
		mscount[mismatch_count[i]]++;
	}

	unsigned int active_count = partitioning_count - mscount[255];

	// Create a running sum from the histogram array
	// Cells store previous values only; i.e. exclude self after sum
	unsigned int summa = 0;
	for (unsigned int i = 0; i < 256; i++)
	{
		unsigned int cnt = mscount[i];
		mscount[i] = summa;
		summa += cnt;
	}

	// Use the running sum as the index, incrementing after read to allow
	// sequential entries with the same count
	for (unsigned int i = 0; i < partitioning_count; i++)
	{
		unsigned int idx = mscount[mismatch_count[i]]++;
		partition_ordering[idx] = i;
	}

	return active_count;
}

/**
 * @brief Use k-means clustering to compute a partition ordering for a block..
 *
 * @param      bsd                  The block size information.
 * @param      blk                  The image block color data to compress.
 * @param      partition_count      The desired number of partitions in the block.
 * @param[out] partition_ordering   The list of recommended partition indices, in priority order.
 *
 * @return The number of active partitionings in this selection.
 */
static unsigned int compute_kmeans_partition_ordering(
	const block_size_descriptor& bsd,
	const image_block& blk,
	unsigned int partition_count,
	unsigned int partition_ordering[BLOCK_MAX_PARTITIONINGS]
) {
	vfloat4 cluster_centers[BLOCK_MAX_PARTITIONS];
	uint8_t texel_partitions[BLOCK_MAX_TEXELS];

	// Use three passes of k-means clustering to partition the block data
	for (unsigned int i = 0; i < 3; i++)
	{
		if (i == 0)
		{
			kmeans_init(blk, bsd.texel_count, partition_count, cluster_centers);
		}
		else
		{
			kmeans_update(blk, bsd.texel_count, partition_count, cluster_centers, texel_partitions);
		}

		kmeans_assign(blk, bsd.texel_count, partition_count, cluster_centers, texel_partitions);
	}

	// Construct the block bitmaps of texel assignments to each partition
	uint64_t bitmaps[BLOCK_MAX_PARTITIONS] { 0 };
	unsigned int texels_to_process = astc::min(bsd.texel_count, BLOCK_MAX_KMEANS_TEXELS);
	promise(texels_to_process > 0);
	for (unsigned int i = 0; i < texels_to_process; i++)
	{
		unsigned int idx = bsd.kmeans_texels[i];
		bitmaps[texel_partitions[idx]] |= 1ULL << i;
	}

	// Count the mismatch between the block and the format's partition tables
	unsigned int mismatch_counts[BLOCK_MAX_PARTITIONINGS];
	count_partition_mismatch_bits(bsd, partition_count, bitmaps, mismatch_counts);

	// Sort the partitions based on the number of mismatched bits
	return get_partition_ordering_by_mismatch_bits(
	    bsd.partitioning_count_selected[partition_count - 1],
	    mismatch_counts, partition_ordering);
}

/* See header for documentation. */
void find_best_partition_candidates(
	const block_size_descriptor& bsd,
	const image_block& blk,
	unsigned int partition_count,
	unsigned int partition_search_limit,
	unsigned int best_partitions[2]
) {
	// Constant used to estimate quantization error for a given partitioning; the optimal value for
	// this depends on bitrate. These values have been determined empirically.
	unsigned int texels_per_block = bsd.texel_count;
	float weight_imprecision_estim = 0.055f;
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

	promise(partition_count > 0);
	promise(partition_search_limit > 0);

	weight_imprecision_estim = weight_imprecision_estim * weight_imprecision_estim;

	unsigned int partition_sequence[BLOCK_MAX_PARTITIONINGS];
	unsigned int sequence_len = compute_kmeans_partition_ordering(bsd, blk, partition_count, partition_sequence);
	partition_search_limit = astc::min(partition_search_limit, sequence_len);

	bool uses_alpha = !blk.is_constant_channel(3);

	// Partitioning errors assuming uncorrelated-chrominance endpoints
	float uncor_best_error { ERROR_CALC_DEFAULT };
	unsigned int uncor_best_partition { 0 };

	// Partitioning errors assuming same-chrominance endpoints
	// Store two so we can always return one different to uncorr
	float samec_best_errors[2] { ERROR_CALC_DEFAULT, ERROR_CALC_DEFAULT };
	unsigned int samec_best_partitions[2] { 0, 0 };

	if (uses_alpha)
	{
		for (unsigned int i = 0; i < partition_search_limit; i++)
		{
			unsigned int partition = partition_sequence[i];
			const auto& pi = bsd.get_raw_partition_info(partition_count, partition);

			// Compute weighting to give to each component in each partition
			partition_metrics pms[BLOCK_MAX_PARTITIONS];

			compute_avgs_and_dirs_4_comp(pi, blk, pms);

			line4 uncor_lines[BLOCK_MAX_PARTITIONS];
			line4 samec_lines[BLOCK_MAX_PARTITIONS];

			processed_line4 uncor_plines[BLOCK_MAX_PARTITIONS];
			processed_line4 samec_plines[BLOCK_MAX_PARTITIONS];

			float uncor_line_lens[BLOCK_MAX_PARTITIONS];
			float samec_line_lens[BLOCK_MAX_PARTITIONS];

			for (unsigned int j = 0; j < partition_count; j++)
			{
				partition_metrics& pm = pms[j];

				uncor_lines[j].a = pm.avg;
				uncor_lines[j].b = normalize_safe(pm.dir, unit4());

				uncor_plines[j].amod = uncor_lines[j].a - uncor_lines[j].b * dot(uncor_lines[j].a, uncor_lines[j].b);
				uncor_plines[j].bs = uncor_lines[j].b;

				samec_lines[j].a = vfloat4::zero();
				samec_lines[j].b = normalize_safe(pm.avg, unit4());

				samec_plines[j].amod = vfloat4::zero();
				samec_plines[j].bs = samec_lines[j].b;
			}

			float uncor_error = 0.0f;
			float samec_error = 0.0f;

			compute_error_squared_rgba(pi,
			                           blk,
			                           uncor_plines,
			                           samec_plines,
			                           uncor_line_lens,
			                           samec_line_lens,
			                           uncor_error,
			                           samec_error);

			// Compute an estimate of error introduced by weight quantization imprecision.
			// This error is computed as follows, for each partition
			//     1: compute the principal-axis vector (full length) in error-space
			//     2: convert the principal-axis vector to regular RGB-space
			//     3: scale the vector by a constant that estimates average quantization error
			//     4: for each texel, square the vector, then do a dot-product with the texel's
			//        error weight; sum up the results across all texels.
			//     4(optimized): square the vector once, then do a dot-product with the average
			//        texel error, then multiply by the number of texels.

			for (unsigned int j = 0; j < partition_count; j++)
			{
				float tpp = static_cast<float>(pi.partition_texel_count[j]);
				vfloat4 error_weights(tpp * weight_imprecision_estim);

				vfloat4 uncor_vector = uncor_lines[j].b * uncor_line_lens[j];
				vfloat4 samec_vector = samec_lines[j].b * samec_line_lens[j];

				uncor_error += dot_s(uncor_vector * uncor_vector, error_weights);
				samec_error += dot_s(samec_vector * samec_vector, error_weights);
			}

			if (uncor_error < uncor_best_error)
			{
				uncor_best_error = uncor_error;
				uncor_best_partition = partition;
			}

			if (samec_error < samec_best_errors[0])
			{
				samec_best_errors[1] = samec_best_errors[0];
				samec_best_partitions[1] = samec_best_partitions[0];

				samec_best_errors[0] = samec_error;
				samec_best_partitions[0] = partition;
			}
			else if (samec_error < samec_best_errors[1])
			{
				samec_best_errors[1] = samec_error;
				samec_best_partitions[1] = partition;
			}
		}
	}
	else
	{
		for (unsigned int i = 0; i < partition_search_limit; i++)
		{
			unsigned int partition = partition_sequence[i];
			const auto& pi = bsd.get_raw_partition_info(partition_count, partition);

			// Compute weighting to give to each component in each partition
			partition_metrics pms[BLOCK_MAX_PARTITIONS];
			compute_avgs_and_dirs_3_comp_rgb(pi, blk, pms);

			partition_lines3 plines[BLOCK_MAX_PARTITIONS];

			for (unsigned int j = 0; j < partition_count; j++)
			{
				partition_metrics& pm = pms[j];
				partition_lines3& pl = plines[j];

				pl.uncor_line.a = pm.avg;
				pl.uncor_line.b = normalize_safe(pm.dir.swz<0, 1, 2>(), unit3());

				pl.samec_line.a = vfloat4::zero();
				pl.samec_line.b = normalize_safe(pm.avg.swz<0, 1, 2>(), unit3());

				pl.uncor_pline.amod = pl.uncor_line.a - pl.uncor_line.b * dot3(pl.uncor_line.a, pl.uncor_line.b);
				pl.uncor_pline.bs   = pl.uncor_line.b;

				pl.samec_pline.amod = vfloat4::zero();
				pl.samec_pline.bs   = pl.samec_line.b;
			}

			float uncor_error = 0.0f;
			float samec_error = 0.0f;

			compute_error_squared_rgb(pi,
			                          blk,
			                          plines,
			                          uncor_error,
			                          samec_error);

			// Compute an estimate of error introduced by weight quantization imprecision.
			// This error is computed as follows, for each partition
			//     1: compute the principal-axis vector (full length) in error-space
			//     2: convert the principal-axis vector to regular RGB-space
			//     3: scale the vector by a constant that estimates average quantization error
			//     4: for each texel, square the vector, then do a dot-product with the texel's
			//        error weight; sum up the results across all texels.
			//     4(optimized): square the vector once, then do a dot-product with the average
			//        texel error, then multiply by the number of texels.

			for (unsigned int j = 0; j < partition_count; j++)
			{
				partition_lines3& pl = plines[j];

				float tpp = static_cast<float>(pi.partition_texel_count[j]);
				vfloat4 error_weights(tpp * weight_imprecision_estim);

				vfloat4 uncor_vector = pl.uncor_line.b * pl.uncor_line_len;
				vfloat4 samec_vector = pl.samec_line.b * pl.samec_line_len;

				uncor_error += dot3_s(uncor_vector * uncor_vector, error_weights);
				samec_error += dot3_s(samec_vector * samec_vector, error_weights);
			}

			if (uncor_error < uncor_best_error)
			{
				uncor_best_error = uncor_error;
				uncor_best_partition = partition;
			}

			if (samec_error < samec_best_errors[0])
			{
				samec_best_errors[1] = samec_best_errors[0];
				samec_best_partitions[1] = samec_best_partitions[0];

				samec_best_errors[0] = samec_error;
				samec_best_partitions[0] = partition;
			}
			else if (samec_error < samec_best_errors[1])
			{
				samec_best_errors[1] = samec_error;
				samec_best_partitions[1] = partition;
			}
		}
	}

	// Same partition is best for both, so use this first unconditionally
	if (uncor_best_partition == samec_best_partitions[0])
	{
		best_partitions[0] = samec_best_partitions[0];
		best_partitions[1] = samec_best_partitions[1];
	}
	// Uncor is best
	else if (uncor_best_error <= samec_best_errors[0])
	{
		best_partitions[0] = uncor_best_partition;
		best_partitions[1] = samec_best_partitions[0];
	}
	// Samec is best
	else
	{
		best_partitions[0] = samec_best_partitions[0];
		best_partitions[1] = uncor_best_partition;
	}

	// Convert these back into canonical partition IDs for the rest of the codec
	best_partitions[0] = bsd.get_raw_partition_info(partition_count, best_partitions[0]).partition_index;
	best_partitions[1] = bsd.get_raw_partition_info(partition_count, best_partitions[1]).partition_index;
}

#endif
