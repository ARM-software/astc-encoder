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
 * @brief Compression quality preset configurations.
 */

#if !defined(ASTCENC_PRESET_CONFIGS_INCLUDED)
#define ASTCENC_PRESET_CONFIGS_INCLUDED

#include <array>

#include "astcenc.h"

/**
 * @brief Record of the quality tuning parameter values.
 *
 * See the @c astcenc_config structure for detailed parameter documentation.
 *
 * Note that the mse_overshoot entries are scaling factors relative to the base MSE to hit db_limit.
 * A 20% overshoot is harder to hit for a higher base db_limit, so we may actually use lower ratios
 * for the more thorough search presets because the underlying db_limit is so much higher.
 */
struct astcenc_preset_config
{
	/** @brief The quality level represented by this preset. */
	float quality;

	/** @brief The maximum number of partitions searched. */
	unsigned int tune_partition_count_limit;

	/** @brief The maximum number of two-partition indices searched. */
	unsigned int tune_2partition_index_limit;

	/** @brief The maximum number of three-partition indices searched. */
	unsigned int tune_3partition_index_limit;

	/** @brief The maximum number of four-partition indices searched. */
	unsigned int tune_4partition_index_limit;

	/** @brief The maximum centile for block modes searched. */
	unsigned int tune_block_mode_limit;

	/** @brief The maximum number of iterative refinements applied. */
	unsigned int tune_refinement_limit;

	/** @brief The number of trial candidates per mode search. */
	unsigned int tune_candidate_limit;

	/** @brief The number of trial two-partitionings per search. */
	unsigned int tune_2partitioning_candidate_limit;

	/** @brief The number of trial three-partitionings per search. */
	unsigned int tune_3partitioning_candidate_limit;

	/** @brief The number of trial four-partitionings per search. */
	unsigned int tune_4partitioning_candidate_limit;

	/** @brief The dB limit base for the 35 * log10(texel_count) adjustment. */
	float tune_db_limit_a_base;

	/** @brief The dB limit base for the 19 * log10(texel_count) adjustment. */
	float tune_db_limit_b_base;

	/** @brief The MSE overshoot needed to early-out trials. */
	float tune_mse_overshoot;

	/** @brief The threshold for skipping three- and four-partition trials. */
	float tune_2partition_early_out_limit_factor;

	/** @brief The threshold for skipping four-partition trials. */
	float tune_3partition_early_out_limit_factor;

	/** @brief The threshold for skipping two weight planes. */
	float tune_2plane_early_out_limit_correlation;

	/** @brief The enable threshold for the mode0 fast-path search. */
	float tune_search_mode0_enable;
};

/**
 * @brief The static presets for high bandwidth encodings (x < 25 texels per block).
 */
static const std::array<astcenc_preset_config, 6> preset_configs_high {{
	// Settings for high bitrate, fastest
	{
		.quality                                 = ASTCENC_PRE_FASTEST,
		.tune_partition_count_limit              = 2,
		.tune_2partition_index_limit             = 10,
		.tune_3partition_index_limit             = 6,
		.tune_4partition_index_limit             = 4,
		.tune_block_mode_limit                   = 43,
		.tune_refinement_limit                   = 2,
		.tune_candidate_limit                    = 2,
		.tune_2partitioning_candidate_limit      = 2,
		.tune_3partitioning_candidate_limit      = 2,
		.tune_4partitioning_candidate_limit      = 2,
		.tune_db_limit_a_base                    = 85.2f,
		.tune_db_limit_b_base                    = 63.2f,
		.tune_mse_overshoot                      = 3.5f,
		.tune_2partition_early_out_limit_factor  = 1.00f,
		.tune_3partition_early_out_limit_factor  = 1.00f,
		.tune_2plane_early_out_limit_correlation = 0.85f,
		.tune_search_mode0_enable                = 0.0f,
	},
	// Settings for high bitrate, fast
	{
		.quality                                 = ASTCENC_PRE_FAST,
		.tune_partition_count_limit              = 3,
		.tune_2partition_index_limit             = 18,
		.tune_3partition_index_limit             = 10,
		.tune_4partition_index_limit             = 8,
		.tune_block_mode_limit                   = 55,
		.tune_refinement_limit                   = 3,
		.tune_candidate_limit                    = 3,
		.tune_2partitioning_candidate_limit      = 2,
		.tune_3partitioning_candidate_limit      = 2,
		.tune_4partitioning_candidate_limit      = 2,
		.tune_db_limit_a_base                    = 85.2f,
		.tune_db_limit_b_base                    = 63.2f,
		.tune_mse_overshoot                      = 3.5f,
		.tune_2partition_early_out_limit_factor  = 1.00f,
		.tune_3partition_early_out_limit_factor  = 1.00f,
		.tune_2plane_early_out_limit_correlation = 0.90f,
		.tune_search_mode0_enable                = 0.0f,
	},
	// Settings for high bitrate, medium
	{
		.quality                                 = ASTCENC_PRE_MEDIUM,
		.tune_partition_count_limit              = 4,
		.tune_2partition_index_limit             = 34,
		.tune_3partition_index_limit             = 28,
		.tune_4partition_index_limit             = 16,
		.tune_block_mode_limit                   = 77,
		.tune_refinement_limit                   = 3,
		.tune_candidate_limit                    = 3,
		.tune_2partitioning_candidate_limit      = 2,
		.tune_3partitioning_candidate_limit      = 2,
		.tune_4partitioning_candidate_limit      = 2,
		.tune_db_limit_a_base                    = 95.0f,
		.tune_db_limit_b_base                    = 70.0f,
		.tune_mse_overshoot                      = 2.5f,
		.tune_2partition_early_out_limit_factor  = 1.10f,
		.tune_3partition_early_out_limit_factor  = 1.05f,
		.tune_2plane_early_out_limit_correlation = 0.95f,
		.tune_search_mode0_enable                = 0.0f,
	},
	// Settings for high bitrate, thorough
	{
		.quality                                 = ASTCENC_PRE_THOROUGH,
		.tune_partition_count_limit              = 4,
		.tune_2partition_index_limit             = 82,
		.tune_3partition_index_limit             = 60,
		.tune_4partition_index_limit             = 30,
		.tune_block_mode_limit                   = 94,
		.tune_refinement_limit                   = 4,
		.tune_candidate_limit                    = 4,
		.tune_2partitioning_candidate_limit      = 3,
		.tune_3partitioning_candidate_limit      = 2,
		.tune_4partitioning_candidate_limit      = 2,
		.tune_db_limit_a_base                    = 105.0f,
		.tune_db_limit_b_base                    = 77.0f,
		.tune_mse_overshoot                      = 10.0f,
		.tune_2partition_early_out_limit_factor  = 1.35f,
		.tune_3partition_early_out_limit_factor  = 1.15f,
		.tune_2plane_early_out_limit_correlation = 0.97f,
		.tune_search_mode0_enable                = 0.0f,
	},
	// Settings for high bitrate, very thorough
	{
		.quality                                 = ASTCENC_PRE_VERYTHOROUGH,
		.tune_partition_count_limit              = 4,
		.tune_2partition_index_limit             = 256,
		.tune_3partition_index_limit             = 128,
		.tune_4partition_index_limit             = 64,
		.tune_block_mode_limit                   = 98,
		.tune_refinement_limit                   = 4,
		.tune_candidate_limit                    = 6,
		.tune_2partitioning_candidate_limit      = 8,
		.tune_3partitioning_candidate_limit      = 6,
		.tune_4partitioning_candidate_limit      = 4,
		.tune_db_limit_a_base                    = 200.0f,
		.tune_db_limit_b_base                    = 200.0f,
		.tune_mse_overshoot                      = 10.0f,
		.tune_2partition_early_out_limit_factor  = 1.60f,
		.tune_3partition_early_out_limit_factor  = 1.40f,
		.tune_2plane_early_out_limit_correlation = 0.98f,
		.tune_search_mode0_enable                = 0.0f,
	},
	// Settings for high bitrate, exhaustive
	{
		.quality                                 = ASTCENC_PRE_EXHAUSTIVE,
		.tune_partition_count_limit              = 4,
		.tune_2partition_index_limit             = 512,
		.tune_3partition_index_limit             = 512,
		.tune_4partition_index_limit             = 512,
		.tune_block_mode_limit                   = 100,
		.tune_refinement_limit                   = 4,
		.tune_candidate_limit                    = 8,
		.tune_2partitioning_candidate_limit      = 8,
		.tune_3partitioning_candidate_limit      = 8,
		.tune_4partitioning_candidate_limit      = 8,
		.tune_db_limit_a_base                    = 200.0f,
		.tune_db_limit_b_base                    = 200.0f,
		.tune_mse_overshoot                      = 10.0f,
		.tune_2partition_early_out_limit_factor  = 2.00f,
		.tune_3partition_early_out_limit_factor  = 2.00f,
		.tune_2plane_early_out_limit_correlation = 0.99f,
		.tune_search_mode0_enable                = 0.0f,
	}
}};

/**
 * @brief The static presets for medium bandwidth encodings (25 <= x < 64 texels per block).
 */
static const std::array<astcenc_preset_config, 6> preset_configs_mid {{
	// Settings for medium bitrate, fastest
	{
		.quality                                 = ASTCENC_PRE_FASTEST,
		.tune_partition_count_limit              = 2,
		.tune_2partition_index_limit             = 10,
		.tune_3partition_index_limit             = 6,
		.tune_4partition_index_limit             = 4,
		.tune_block_mode_limit                   = 43,
		.tune_refinement_limit                   = 2,
		.tune_candidate_limit                    = 2,
		.tune_2partitioning_candidate_limit      = 2,
		.tune_3partitioning_candidate_limit      = 2,
		.tune_4partitioning_candidate_limit      = 2,
		.tune_db_limit_a_base                    = 85.2f,
		.tune_db_limit_b_base                    = 63.2f,
		.tune_mse_overshoot                      = 3.5f,
		.tune_2partition_early_out_limit_factor  = 1.00f,
		.tune_3partition_early_out_limit_factor  = 1.00f,
		.tune_2plane_early_out_limit_correlation = 0.80f,
		.tune_search_mode0_enable                = 1.0f,
	},
	// Settings for medium bitrate, fast
	{
		.quality                                 = ASTCENC_PRE_FAST,
		.tune_partition_count_limit              = 3,
		.tune_2partition_index_limit             = 18,
		.tune_3partition_index_limit             = 12,
		.tune_4partition_index_limit             = 10,
		.tune_block_mode_limit                   = 55,
		.tune_refinement_limit                   = 3,
		.tune_candidate_limit                    = 3,
		.tune_2partitioning_candidate_limit      = 2,
		.tune_3partitioning_candidate_limit      = 2,
		.tune_4partitioning_candidate_limit      = 2,
		.tune_db_limit_a_base                    = 85.2f,
		.tune_db_limit_b_base                    = 63.2f,
		.tune_mse_overshoot                      = 3.5f,
		.tune_2partition_early_out_limit_factor  = 1.00f,
		.tune_3partition_early_out_limit_factor  = 1.00f,
		.tune_2plane_early_out_limit_correlation = 0.85f,
		.tune_search_mode0_enable                = 1.0f,
	},
	// Settings for medium bitrate, medium
	{
		.quality                                 = ASTCENC_PRE_MEDIUM,
		.tune_partition_count_limit              = 3,
		.tune_2partition_index_limit             = 34,
		.tune_3partition_index_limit             = 28,
		.tune_4partition_index_limit             = 16,
		.tune_block_mode_limit                   = 77,
		.tune_refinement_limit                   = 3,
		.tune_candidate_limit                    = 3,
		.tune_2partitioning_candidate_limit      = 2,
		.tune_3partitioning_candidate_limit      = 2,
		.tune_4partitioning_candidate_limit      = 2,
		.tune_db_limit_a_base                    = 95.0f,
		.tune_db_limit_b_base                    = 70.0f,
		.tune_mse_overshoot                      = 3.0f,
		.tune_2partition_early_out_limit_factor  = 1.10f,
		.tune_3partition_early_out_limit_factor  = 1.05f,
		.tune_2plane_early_out_limit_correlation = 0.90f,
		.tune_search_mode0_enable                = 1.0f,
	},
	// Settings for medium bitrate, thorough
	{
		.quality                                 = ASTCENC_PRE_THOROUGH,
		.tune_partition_count_limit              = 4,
		.tune_2partition_index_limit             = 82,
		.tune_3partition_index_limit             = 60,
		.tune_4partition_index_limit             = 30,
		.tune_block_mode_limit                   = 94,
		.tune_refinement_limit                   = 4,
		.tune_candidate_limit                    = 4,
		.tune_2partitioning_candidate_limit      = 3,
		.tune_3partitioning_candidate_limit      = 2,
		.tune_4partitioning_candidate_limit      = 2,
		.tune_db_limit_a_base                    = 105.0f,
		.tune_db_limit_b_base                    = 77.0f,
		.tune_mse_overshoot                      = 10.0f,
		.tune_2partition_early_out_limit_factor  = 1.40f,
		.tune_3partition_early_out_limit_factor  = 1.20f,
		.tune_2plane_early_out_limit_correlation = 0.95f,
		.tune_search_mode0_enable                = 0.0f,
	},
	// Settings for medium bitrate, very thorough
	{
		.quality                                 = ASTCENC_PRE_VERYTHOROUGH,
		.tune_partition_count_limit              = 4,
		.tune_2partition_index_limit             = 256,
		.tune_3partition_index_limit             = 128,
		.tune_4partition_index_limit             = 64,
		.tune_block_mode_limit                   = 98,
		.tune_refinement_limit                   = 4,
		.tune_candidate_limit                    = 6,
		.tune_2partitioning_candidate_limit      = 8,
		.tune_3partitioning_candidate_limit      = 6,
		.tune_4partitioning_candidate_limit      = 3,
		.tune_db_limit_a_base                    = 200.0f,
		.tune_db_limit_b_base                    = 200.0f,
		.tune_mse_overshoot                      = 10.0f,
		.tune_2partition_early_out_limit_factor  = 1.60f,
		.tune_3partition_early_out_limit_factor  = 1.40f,
		.tune_2plane_early_out_limit_correlation = 0.98f,
		.tune_search_mode0_enable                = 0.0f,
	},
	// Settings for medium bitrate, exhaustive
	{
		.quality                                 = ASTCENC_PRE_EXHAUSTIVE,
		.tune_partition_count_limit              = 4,
		.tune_2partition_index_limit             = 256,
		.tune_3partition_index_limit             = 256,
		.tune_4partition_index_limit             = 256,
		.tune_block_mode_limit                   = 100,
		.tune_refinement_limit                   = 4,
		.tune_candidate_limit                    = 8,
		.tune_2partitioning_candidate_limit      = 8,
		.tune_3partitioning_candidate_limit      = 8,
		.tune_4partitioning_candidate_limit      = 8,
		.tune_db_limit_a_base                    = 200.0f,
		.tune_db_limit_b_base                    = 200.0f,
		.tune_mse_overshoot                      = 10.0f,
		.tune_2partition_early_out_limit_factor  = 2.00f,
		.tune_3partition_early_out_limit_factor  = 2.00f,
		.tune_2plane_early_out_limit_correlation = 0.99f,
		.tune_search_mode0_enable                = 0.0f,
	}
}};

/**
 * @brief The static presets for low bandwidth encodings (64 <= x texels per block).
 */
static const std::array<astcenc_preset_config, 6> preset_configs_low {{
	// Settings for low bitrate, fastest
	{
		.quality                                 = ASTCENC_PRE_FASTEST,
		.tune_partition_count_limit              = 2,
		.tune_2partition_index_limit             = 10,
		.tune_3partition_index_limit             = 6,
		.tune_4partition_index_limit             = 4,
		.tune_block_mode_limit                   = 40,
		.tune_refinement_limit                   = 2,
		.tune_candidate_limit                    = 2,
		.tune_2partitioning_candidate_limit      = 2,
		.tune_3partitioning_candidate_limit      = 2,
		.tune_4partitioning_candidate_limit      = 2,
		.tune_db_limit_a_base                    = 85.0f,
		.tune_db_limit_b_base                    = 63.0f,
		.tune_mse_overshoot                      = 3.5f,
		.tune_2partition_early_out_limit_factor  = 1.00f,
		.tune_3partition_early_out_limit_factor  = 1.00f,
		.tune_2plane_early_out_limit_correlation = 0.80f,
		.tune_search_mode0_enable                = 1.0f,
	},
	// Settings for low bitrate, fast
	{
		.quality                                 = ASTCENC_PRE_FAST,
		.tune_partition_count_limit              = 2,
		.tune_2partition_index_limit             = 18,
		.tune_3partition_index_limit             = 12,
		.tune_4partition_index_limit             = 10,
		.tune_block_mode_limit                   = 55,
		.tune_refinement_limit                   = 3,
		.tune_candidate_limit                    = 3,
		.tune_2partitioning_candidate_limit      = 2,
		.tune_3partitioning_candidate_limit      = 2,
		.tune_4partitioning_candidate_limit      = 2,
		.tune_db_limit_a_base                    = 85.0f,
		.tune_db_limit_b_base                    = 63.0f,
		.tune_mse_overshoot                      = 3.5f,
		.tune_2partition_early_out_limit_factor  = 1.00f,
		.tune_3partition_early_out_limit_factor  = 1.00f,
		.tune_2plane_early_out_limit_correlation = 0.85f,
		.tune_search_mode0_enable                = 1.0f,
	},
	// Settings for low bitrate, medium
	{
		.quality                                 = ASTCENC_PRE_MEDIUM,
		.tune_partition_count_limit              = 3,
		.tune_2partition_index_limit             = 34,
		.tune_3partition_index_limit             = 28,
		.tune_4partition_index_limit             = 16,
		.tune_block_mode_limit                   = 77,
		.tune_refinement_limit                   = 3,
		.tune_candidate_limit                    = 3,
		.tune_2partitioning_candidate_limit      = 2,
		.tune_3partitioning_candidate_limit      = 2,
		.tune_4partitioning_candidate_limit      = 2,
		.tune_db_limit_a_base                    = 95.0f,
		.tune_db_limit_b_base                    = 70.0f,
		.tune_mse_overshoot                      = 3.5f,
		.tune_2partition_early_out_limit_factor  = 1.10f,
		.tune_3partition_early_out_limit_factor  = 1.05f,
		.tune_2plane_early_out_limit_correlation = 0.90f,
		.tune_search_mode0_enable                = 1.0f,
	},
	// Settings for low bitrate, thorough
	{
		.quality                                 = ASTCENC_PRE_THOROUGH,
		.tune_partition_count_limit              = 4,
		.tune_2partition_index_limit             = 82,
		.tune_3partition_index_limit             = 60,
		.tune_4partition_index_limit             = 30,
		.tune_block_mode_limit                   = 93,
		.tune_refinement_limit                   = 4,
		.tune_candidate_limit                    = 4,
		.tune_2partitioning_candidate_limit      = 3,
		.tune_3partitioning_candidate_limit      = 2,
		.tune_4partitioning_candidate_limit      = 2,
		.tune_db_limit_a_base                    = 105.0f,
		.tune_db_limit_b_base                    = 77.0f,
		.tune_mse_overshoot                      = 10.0f,
		.tune_2partition_early_out_limit_factor  = 1.30f,
		.tune_3partition_early_out_limit_factor  = 1.20f,
		.tune_2plane_early_out_limit_correlation = 0.97f,
		.tune_search_mode0_enable                = 1.0f,
	},
	// Settings for low bitrate, very thorough
	{
		.quality                                 = ASTCENC_PRE_VERYTHOROUGH,
		.tune_partition_count_limit              = 4,
		.tune_2partition_index_limit             = 256,
		.tune_3partition_index_limit             = 128,
		.tune_4partition_index_limit             = 64,
		.tune_block_mode_limit                   = 98,
		.tune_refinement_limit                   = 4,
		.tune_candidate_limit                    = 6,
		.tune_2partitioning_candidate_limit      = 8,
		.tune_3partitioning_candidate_limit      = 5,
		.tune_4partitioning_candidate_limit      = 2,
		.tune_db_limit_a_base                    = 200.0f,
		.tune_db_limit_b_base                    = 200.0f,
		.tune_mse_overshoot                      = 10.0f,
		.tune_2partition_early_out_limit_factor  = 1.60f,
		.tune_3partition_early_out_limit_factor  = 1.40f,
		.tune_2plane_early_out_limit_correlation = 0.98f,
		.tune_search_mode0_enable                = 1.0f,
	},
	// Settings for low bitrate, exhaustive
	{
		.quality                                 = ASTCENC_PRE_EXHAUSTIVE,
		.tune_partition_count_limit              = 4,
		.tune_2partition_index_limit             = 256,
		.tune_3partition_index_limit             = 256,
		.tune_4partition_index_limit             = 256,
		.tune_block_mode_limit                   = 100,
		.tune_refinement_limit                   = 4,
		.tune_candidate_limit                    = 8,
		.tune_2partitioning_candidate_limit      = 8,
		.tune_3partitioning_candidate_limit      = 8,
		.tune_4partitioning_candidate_limit      = 8,
		.tune_db_limit_a_base                    = 200.0f,
		.tune_db_limit_b_base                    = 200.0f,
		.tune_mse_overshoot                      = 10.0f,
		.tune_2partition_early_out_limit_factor  = 2.00f,
		.tune_3partition_early_out_limit_factor  = 2.00f,
		.tune_2plane_early_out_limit_correlation = 0.99f,
		.tune_search_mode0_enable                = 1.0f,
	}
}};

#endif  // #ifndef ASTCENC_PRESET_CONFIGS_INCLUDED
