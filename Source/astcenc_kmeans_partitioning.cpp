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
 * @brief Functions for approximate partitioning by kmeans clustering.
 *
 * Do this in 2 stages:
 * 1: basic clustering, a couple of passes just to get a few clusters
 * 2: clustering based on line, a few passes until it seems to stabilize.
 *
 * After clustering is done, we use the clustering result to construct one
 * bitmap for each partition. We then scan though the partition table, counting
 * how well the bitmaps matched.
  */

#include "astcenc_internal.h"

// for k++ means, we need pseudo-random numbers, however using random numbers
// directly results in unreproducible encoding results. As such, we will
// instead just supply a handful of numbers from random.org, and apply an
// algorithm similar to XKCD #221. (http://xkcd.com/221/)

// cluster the texels using the k++ means clustering initialization algorithm.
static void kpp_initialize(
	int xdim,
	int ydim,
	int zdim,
	int partition_count,
	const imageblock* blk,
	float4* cluster_centers
) {
	int texels_per_block = xdim * ydim * zdim;

	int cluster_center_samples[4];
	// pick a random sample as first center-point.
	cluster_center_samples[0] = 145897 /* number from random.org */  % texels_per_block;
	int samples_selected = 1;

	float distances[MAX_TEXELS_PER_BLOCK];

	// compute the distance to the first point.
	int sample = cluster_center_samples[0];
	float4 center_color = float4(blk->data_r[sample],
								 blk->data_g[sample],
								 blk->data_b[sample],
								 blk->data_a[sample]);

	float distance_sum = 0.0f;
	for (int i = 0; i < texels_per_block; i++)
	{
		float4 color = float4(blk->data_r[i],
							  blk->data_g[i],
							  blk->data_b[i],
							  blk->data_a[i]);
		float4 diff = color - center_color;
		float distance = dot(diff, diff);
		distance_sum += distance;
		distances[i] = distance;
	}

	// more numbers from random.org
	float cluster_cutoffs[25] = {
		0.952312f, 0.206893f, 0.835984f, 0.507813f, 0.466170f,
		0.872331f, 0.488028f, 0.866394f, 0.363093f, 0.467905f,
		0.812967f, 0.626220f, 0.932770f, 0.275454f, 0.832020f,
		0.362217f, 0.318558f, 0.240113f, 0.009190f, 0.983995f,
		0.566812f, 0.347661f, 0.731960f, 0.156391f, 0.297786f
	};

	while (1)
	{
		// pick a point in a weighted-random fashion.
		float summa = 0.0f;
		float distance_cutoff = distance_sum * cluster_cutoffs[samples_selected + 5 * partition_count];
		for (sample = 0; sample < texels_per_block; sample++)
		{
			summa += distances[sample];
			if (summa >= distance_cutoff)
			{
				break;
			}
		}

		if (sample >= texels_per_block)
		{
			sample = texels_per_block - 1;
		}

		cluster_center_samples[samples_selected] = sample;
		samples_selected++;
		if (samples_selected >= partition_count)
		{
			break;
		}

		// update the distances with the new point.
		center_color = float4(blk->data_r[sample],
		                      blk->data_g[sample],
		                      blk->data_b[sample],
		                      blk->data_a[sample]);

		distance_sum = 0.0f;
		for (int i = 0; i < texels_per_block; i++)
		{
			float4 color = float4(blk->data_r[i],
			                      blk->data_g[i],
			                      blk->data_b[i],
			                      blk->data_a[i]);
			float4 diff = color - center_color;
			float distance = dot(diff, diff);
			distance = MIN(distance, distances[i]);
			distance_sum += distance;
			distances[i] = distance;
		}
	}

	// finally, gather up the results.
	for (int i = 0; i < partition_count; i++)
	{
		int center_sample = cluster_center_samples[i];
		float4 color = float4(blk->data_r[center_sample],
		                      blk->data_g[center_sample],
		                      blk->data_b[center_sample],
		                      blk->data_a[center_sample]);
		cluster_centers[i] = color;
	}
}

// basic K-means clustering: given a set of cluster centers,
// assign each texel to a partition
static void basic_kmeans_assign_pass(
	int xdim,
	int ydim,
	int zdim,
	int partition_count,
	const imageblock* blk,
	const float4* cluster_centers,
	int* partition_of_texel
) {
	int texels_per_block = xdim * ydim * zdim;

	float distances[MAX_TEXELS_PER_BLOCK];

	int texels_per_partition[4];

	texels_per_partition[0] = texels_per_block;
	for (int i = 1; i < partition_count; i++)
	{
		texels_per_partition[i] = 0;
	}

	for (int i = 0; i < texels_per_block; i++)
	{
		float4 color = float4(blk->data_r[i],
							  blk->data_g[i],
							  blk->data_b[i],
							  blk->data_a[i]);
		float4 diff = color - cluster_centers[0];
		float distance = dot(diff, diff);
		distances[i] = distance;
		partition_of_texel[i] = 0;
	}

	for (int j = 1; j < partition_count; j++)
	{
		float4 center_color = cluster_centers[j];

		for (int i = 0; i < texels_per_block; i++)
		{
			float4 color = float4(blk->data_r[i],
								  blk->data_g[i],
								  blk->data_b[i],
								  blk->data_a[i]);
			float4 diff = color - center_color;
			float distance = dot(diff, diff);
			if (distance < distances[i])
			{
				distances[i] = distance;
				texels_per_partition[partition_of_texel[i]]--;
				texels_per_partition[j]++;
				partition_of_texel[i] = j;
			}
		}
	}

	// it is possible to get a situation where one of the partitions ends up
	// without any texels. In this case, we assign texel N to partition N;
	// this is silly, but ensures that every partition retains at least one texel.
	// Reassigning a texel in this manner may cause another partition to go empty,
	// so if we actually did a reassignment, we run the whole loop over again.
	int problem_case;
	do
	{
		problem_case = 0;
		for (int i = 0; i < partition_count; i++)
		{
			if (texels_per_partition[i] == 0)
			{
				texels_per_partition[partition_of_texel[i]]--;
				texels_per_partition[i]++;
				partition_of_texel[i] = i;
				problem_case = 1;
			}
		}
	}
	while (problem_case != 0);
}

// basic k-means clustering: given a set of cluster assignments
// for the texels, find the center position of each cluster.
static void basic_kmeans_update(
	int xdim,
	int ydim,
	int zdim,
	int partition_count,
	const imageblock* blk,
	const int* partition_of_texel,
	float4* cluster_centers
) {
	int texels_per_block = xdim * ydim * zdim;

	float4 color_sum[4];
	int weight_sum[4];

	for (int i = 0; i < partition_count; i++)
	{
		color_sum[i] = float4(0.0f, 0.0f, 0.0f, 0.0f);
		weight_sum[i] = 0;
	}

	// first, find the center-of-gravity in each cluster
	for (int i = 0; i < texels_per_block; i++)
	{
		float4 color = float4(blk->data_r[i],
		                      blk->data_g[i],
		                      blk->data_b[i],
		                      blk->data_a[i]);
		int part = partition_of_texel[i];
		color_sum[part] = color_sum[part] + color;
		weight_sum[part]++;
	}

	for (int i = 0; i < partition_count; i++)
	{
		cluster_centers[i] = color_sum[i] * (1.0f / weight_sum[i]);
	}
}

// compute the bit-mismatch for a partitioning in 2-partition mode
static inline int partition_mismatch2(
	uint64_t a0,
	uint64_t a1,
	uint64_t b0,
	uint64_t b1
) {
	int v1 = astc::popcount(a0 ^ b0) + astc::popcount(a1 ^ b1);
	int v2 = astc::popcount(a0 ^ b1) + astc::popcount(a1 ^ b0);
	return MIN(v1, v2);
}

// compute the bit-mismatch for a partitioning in 3-partition mode
static inline int partition_mismatch3(
	uint64_t a0,
	uint64_t a1,
	uint64_t a2,
	uint64_t b0,
	uint64_t b1,
	uint64_t b2
) {
	int p00 = astc::popcount(a0 ^ b0);
	int p01 = astc::popcount(a0 ^ b1);
	int p02 = astc::popcount(a0 ^ b2);

	int p10 = astc::popcount(a1 ^ b0);
	int p11 = astc::popcount(a1 ^ b1);
	int p12 = astc::popcount(a1 ^ b2);

	int p20 = astc::popcount(a2 ^ b0);
	int p21 = astc::popcount(a2 ^ b1);
	int p22 = astc::popcount(a2 ^ b2);

	int s0 = p11 + p22;
	int s1 = p12 + p21;
	int v0 = MIN(s0, s1) + p00;

	int s2 = p10 + p22;
	int s3 = p12 + p20;
	int v1 = MIN(s2, s3) + p01;

	int s4 = p10 + p21;
	int s5 = p11 + p20;
	int v2 = MIN(s4, s5) + p02;

	if (v1 < v0)
		v0 = v1;
	if (v2 < v0)
		v0 = v2;

	return v0;
}

static inline int MIN3(
	int a,
	int b,
	int c
) {
	int d = MIN(a, b);
	return MIN(c, d);
}

// compute the bit-mismatch for a partitioning in 4-partition mode
static inline int partition_mismatch4(
	uint64_t a0,
	uint64_t a1,
	uint64_t a2,
	uint64_t a3,
	uint64_t b0,
	uint64_t b1,
	uint64_t b2,
	uint64_t b3
) {
	int p00 = astc::popcount(a0 ^ b0);
	int p01 = astc::popcount(a0 ^ b1);
	int p02 = astc::popcount(a0 ^ b2);
	int p03 = astc::popcount(a0 ^ b3);

	int p10 = astc::popcount(a1 ^ b0);
	int p11 = astc::popcount(a1 ^ b1);
	int p12 = astc::popcount(a1 ^ b2);
	int p13 = astc::popcount(a1 ^ b3);

	int p20 = astc::popcount(a2 ^ b0);
	int p21 = astc::popcount(a2 ^ b1);
	int p22 = astc::popcount(a2 ^ b2);
	int p23 = astc::popcount(a2 ^ b3);

	int p30 = astc::popcount(a3 ^ b0);
	int p31 = astc::popcount(a3 ^ b1);
	int p32 = astc::popcount(a3 ^ b2);
	int p33 = astc::popcount(a3 ^ b3);

	int mx23 = MIN(p22 + p33, p23 + p32);
	int mx13 = MIN(p21 + p33, p23 + p31);
	int mx12 = MIN(p21 + p32, p22 + p31);
	int mx03 = MIN(p20 + p33, p23 + p30);
	int mx02 = MIN(p20 + p32, p22 + p30);
	int mx01 = MIN(p21 + p30, p20 + p31);

	int v0 = p00 + MIN3(p11 + mx23, p12 + mx13, p13 + mx12);
	int v1 = p01 + MIN3(p10 + mx23, p12 + mx03, p13 + mx02);
	int v2 = p02 + MIN3(p11 + mx03, p10 + mx13, p13 + mx01);
	int v3 = p03 + MIN3(p11 + mx02, p12 + mx01, p10 + mx12);

	int x0 = MIN(v0, v1);
	int x1 = MIN(v2, v3);
	return MIN(x0, x1);
}

static void count_partition_mismatch_bits(
	const block_size_descriptor* bsd,
	int partition_count,
	const uint64_t bitmaps[4],
	int bitcounts[PARTITION_COUNT]
) {
	const partition_info *pi = get_partition_table(bsd, partition_count);

	if (partition_count == 2)
	{
		uint64_t bm0 = bitmaps[0];
		uint64_t bm1 = bitmaps[1];
		for (int i = 0; i < PARTITION_COUNT; i++)
		{
			if (pi->partition_count == 2)
			{
				bitcounts[i] = partition_mismatch2(bm0, bm1, pi->coverage_bitmaps[0], pi->coverage_bitmaps[1]);
			}
			else
			{
				bitcounts[i] = 255;
			}
			pi++;
		}
	}
	else if (partition_count == 3)
	{
		uint64_t bm0 = bitmaps[0];
		uint64_t bm1 = bitmaps[1];
		uint64_t bm2 = bitmaps[2];
		for (int i = 0; i < PARTITION_COUNT; i++)
		{
			if (pi->partition_count == 3)
			{
				bitcounts[i] = partition_mismatch3(bm0, bm1, bm2, pi->coverage_bitmaps[0], pi->coverage_bitmaps[1], pi->coverage_bitmaps[2]);
			}
			else
			{
				bitcounts[i] = 255;
			}
			pi++;
		}
	}
	else if (partition_count == 4)
	{
		uint64_t bm0 = bitmaps[0];
		uint64_t bm1 = bitmaps[1];
		uint64_t bm2 = bitmaps[2];
		uint64_t bm3 = bitmaps[3];
		for (int i = 0; i < PARTITION_COUNT; i++)
		{
			if (pi->partition_count == 4)
			{
				bitcounts[i] = partition_mismatch4(bm0, bm1, bm2, bm3, pi->coverage_bitmaps[0], pi->coverage_bitmaps[1], pi->coverage_bitmaps[2], pi->coverage_bitmaps[3]);
			}
			else
			{
				bitcounts[i] = 255;
			}
			pi++;
		}
	}

}

// counting-sort on the mismatch-bits, thereby
// sorting the partitions into an ordering.
static void get_partition_ordering_by_mismatch_bits(
	const int mismatch_bits[PARTITION_COUNT],
	int partition_ordering[PARTITION_COUNT]
) {
	int mscount[256];
	for (int i = 0; i < 256; i++)
	{
		mscount[i] = 0;
	}

	for (int i = 0; i < PARTITION_COUNT; i++)
	{
		mscount[mismatch_bits[i]]++;
	}

	int summa = 0;
	for (int i = 0; i < 256; i++)
	{
		int cnt = mscount[i];
		mscount[i] = summa;
		summa += cnt;
	}

	for (int i = 0; i < PARTITION_COUNT; i++)
	{
		int idx = mscount[mismatch_bits[i]]++;
		partition_ordering[idx] = i;
	}
}

void kmeans_compute_partition_ordering(
	const block_size_descriptor* bsd,
	int partition_count,
	const imageblock* blk,
	int* ordering
) {
	float4 cluster_centers[4];
	int partition_of_texel[MAX_TEXELS_PER_BLOCK];

	// 3 passes of plain k-means partitioning
	for (int i = 0; i < 3; i++)
	{
		if (i == 0)
		{
			kpp_initialize(bsd->xdim, bsd->ydim, bsd->zdim, partition_count, blk, cluster_centers);
		}
		else
		{
			basic_kmeans_update(bsd->xdim, bsd->ydim, bsd->zdim, partition_count, blk, partition_of_texel, cluster_centers);
		}

		basic_kmeans_assign_pass(bsd->xdim, bsd->ydim, bsd->zdim, partition_count, blk, cluster_centers, partition_of_texel);
	}

	// at this point, we have a near-ideal partitioning.

	// construct bitmaps
	uint64_t bitmaps[4];
	for (int i = 0; i < 4; i++)
	{
		bitmaps[i] = 0ULL;
	}

	int texels_to_process = bsd->texelcount_for_bitmap_partitioning;
	for (int i = 0; i < texels_to_process; i++)
	{
		int idx = bsd->texels_for_bitmap_partitioning[i];
		bitmaps[partition_of_texel[idx]] |= 1ULL << i;
	}

	int bitcounts[PARTITION_COUNT];
	// for each entry in the partition table, count bits of partition-mismatch.
	count_partition_mismatch_bits(bsd, partition_count, bitmaps, bitcounts);

	// finally, sort the partitions by bits-of-partition-mismatch
	get_partition_ordering_by_mismatch_bits(bitcounts, ordering);
}

#endif
