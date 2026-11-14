// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2026 Arm Limited
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
 * @file
 * @brief RDO fixed-offset byte search using DEFLATE-aware rate estimates.
 *
 * Based on Rich Geldreich's public-domain ERT implementation in bc7enc_rdo.
 * See https://github.com/richgel999/bc7enc_rdo/blob/master/ert.cpp.
 */

#if !defined(ASTCENC_DECOMPRESS_ONLY)

#include "astcenc_internal.h"
#include "astcenc_rdo.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstring>
#include <limits>

/** @brief The physical block size in bytes. */
static constexpr unsigned int block_size = 16;

/** @brief The minimum LZ match length in bytes. */
static constexpr unsigned int min_match_length = 3;

static_assert((brick_layout::block_count > 0) && (brick_layout::block_count <= 256),
              "RDO block tags must fit in eight bits");

static_assert(brick_layout::block_count <= std::numeric_limits<int>::max() / block_size,
              "RDO byte offsets must fit in a signed int");

/** @brief Shared scoring inputs and the best result for one block's byte searches. */
struct rdo_search_state
{
	/** @brief The original physical bytes, unchanged until both searches finish. */
	const uint8_t* original;

	/** @brief The distortion and history settings. */
	const rdo_parameters& params;

	/** @brief Rate estimates prepared for the current block. */
	const rdo_rate_model& rates;

	/** @brief Source data and reconstruction caches used for distortion scoring. */
	rdo_brick_context& ctx;

	/** @brief The current block's index in the compact brick stream. */
	unsigned int block_index;

	/** @brief The original block's estimated storage cost in bits. */
	float original_rate;

	/** @brief The distortion multiplier protecting smooth blocks. */
	float mse_scale;

	/** @brief The reciprocal of the distortion multiplier. */
	float mse_recip;

	/** @brief The maximum trial MSE allowed by the RMS guard. */
	float max_mse;

	/** @brief The lowest rate-distortion cost found so far. */
	float best_cost;

	/** @brief The physical bytes producing the best cost. */
	uint8_t best_block[block_size];
};

/**
 * @brief Build byte prefixes for selecting a copied interval.
 *
 * @return The masks for prefixes of zero to sixteen bytes.
 */
static constexpr std::array<std::array<uint8_t, 16>, 17> rdo_make_prefix_masks()
{
	std::array<std::array<uint8_t, 16>, 17> masks {};
	for (unsigned int end = 0; end <= 16; end++)
	{
		for (unsigned int i = 0; i < end; i++)
		{
			masks[end][i] = 0xFF;
		}
	}

	return masks;
}

/** @brief Byte prefixes used to select a copied interval without a variable-size memcpy. */
static constexpr auto prefix_masks = rdo_make_prefix_masks();

/** @brief Estimated cost of continuing the previous block's match, in bits. */
static constexpr float continuation_bits = 1.0f;

/** @brief Polynomial base for reusable substring hashes. */
static constexpr uint32_t hash_base = 16777619u;

/**
 * @brief Build powers for removing a substring's preceding bytes.
 *
 * @return The polynomial base raised to powers zero through sixteen.
 */
static constexpr std::array<uint32_t, 17> rdo_make_hash_powers()
{
	std::array<uint32_t, 17> powers {};
	powers[0] = 1;
	for (unsigned int i = 1; i <= 16; i++)
	{
		powers[i] = powers[i - 1] * hash_base;
	}

	return powers;
}

/** @brief Polynomial factors for extracting each substring length. */
static constexpr auto hash_powers = rdo_make_hash_powers();

/* See header for documentation. */
void rdo_prepare_substring_hashes(
	const uint8_t data[16],
	uint32_t prefixes[17]
) {
	prefixes[0] = 0;
	for (unsigned int i = 0; i < block_size; i++)
	{
		prefixes[i + 1] = prefixes[i] * hash_base + data[i] + 1;
	}
}

/* See header for documentation. */
void rdo_hash_substrings(
	const uint32_t* prefixes,
	unsigned int length,
	uint32_t hashes[16]
) {
	assert((length >= min_match_length) && (length <= block_size));
	unsigned int count = block_size - length + 1;
	for (unsigned int offset = 0; offset < count; offset += ASTCENC_SIMD_WIDTH)
	{
		vint starts(reinterpret_cast<const int*>(prefixes + offset));
		vint ends(reinterpret_cast<const int*>(prefixes + offset + length));
		vint hash = ends - starts * vint(static_cast<int>(hash_powers[length]));
		vint offsets = vint(offset) + vint::lane_id();
		hash += offsets * vint(0x1E35A7BD) + vint(length);
		hash = (hash ^ lsr<16>(hash)) * vint(0x7FEB352D);
		hash = (hash ^ lsr<15>(hash)) * vint(0x2C1B3C6D);
		hash = hash ^ lsr<16>(hash);
		store_lanes_masked(reinterpret_cast<uint8_t*>(hashes + offset), hash, offsets < vint(count));
	}
}

/**
 * @brief Find bytes whose replacement would change protected layout bits.
 *
 * @param[in] original         The current target block.
 * @param[in] source           The predecessor block supplying copied bytes.
 * @param[in] protected_mask   The target block's protected bits.
 *
 * @return A bit mask with one bit for each incompatible source byte.
 */
static uint32_t rdo_get_incompatible_bytes(
	const uint8_t original[16],
	const uint8_t source[16],
	const uint8_t protected_mask[16]
) {
	static_assert((block_size % ASTCENC_SIMD_WIDTH) == 0);
	uint32_t incompatible = 0;
	for (unsigned int i = 0; i < block_size; i += ASTCENC_SIMD_WIDTH)
	{
		vint changed = (vint(original + i) ^ vint(source + i)) & vint(protected_mask + i);
		incompatible |= mask(changed != vint(0)) << i;
	}

	return incompatible;
}

/**
 * @brief Test whether a trial exceeds the physical bit-change search budget.
 *
 * @param[in] original   The original physical block.
 * @param[in] trial      The trial physical block.
 *
 * @return True if the configured bit-change budget is exceeded.
 */
static bool rdo_is_distant_trial(
	const uint8_t original[16],
	const uint8_t trial[16]
) {
	uint64_t original_words[2];
	uint64_t trial_words[2];
	std::memcpy(original_words, original, sizeof(original_words));
	std::memcpy(trial_words, trial, sizeof(trial_words));
	unsigned int changed_bits = std::popcount(original_words[0] ^ trial_words[0])
	                          + std::popcount(original_words[1] ^ trial_words[1]);
	return changed_bits > rdo_max_changed_bits;
}

/** @brief Padded row offsets for each cached substring length. */
static constexpr unsigned int substring_offsets[17] {
	0, 0, 0, 0, 16, 32, 44, 56, 68, 80, 88, 96, 104, 112, 116, 120, 124
};

/**
 * @brief Score a trial and retain it if it improves the block's rate-distortion cost.
 *
 * @param[in,out] search      The current block's scoring state and best result.
 * @param[in]     trial       The candidate physical bytes.
 * @param         rate        The candidate's estimated storage cost in bits.
 * @param         byte_mask   The physical bytes that may differ from the original.
 *
 * @return True if the trial improves the best cost.
 */
static bool rdo_evaluate_trial(
	rdo_search_state& search,
	const uint8_t trial[block_size],
	float rate,
	unsigned int byte_mask
) {
	// Spend distortion work only on trials predicted to reduce storage
	if (rate >= search.original_rate)
	{
		return false;
	}

	// Convert remaining cost headroom into a distortion cutoff
	float rate_cost = rate * search.params.lambda;
	float mse_limit = std::min(search.max_mse,
	                           rdo_get_mse_limit(search.best_cost, rate_cost, search.mse_recip));
	if (mse_limit <= 0.0f)
	{
		return false;
	}

	// Large physical changes rarely preserve useful endpoints and weights
	if (rdo_is_distant_trial(search.original, trial))
	{
		return false;
	}

	float trial_mse = rdo_block_difference(search.ctx, trial, byte_mask, mse_limit);
	if (trial_mse >= search.max_mse)
	{
		return false;
	}

	float cost = trial_mse * search.mse_scale + rate_cost;
	if (!(cost < search.best_cost))
	{
		return false;
	}

	search.best_cost = cost;
	std::memcpy(search.best_block, trial, block_size);
	return true;
}

/**
 * @brief Search predecessor blocks for the best single copied byte range.
 *
 * @param[in,out] search                The current block's scoring state and best result.
 * @param[in]     data                  The compact brick stream with finalized predecessors.
 * @param[in]     mask                  The current block's protected layout bits.
 * @param         search_start          The first predecessor to search.
 * @param         continuation_source   The next source byte of a preceding match, or -1.
 * @param[in]     source_hashes         Cached substring hashes for finalized predecessors.
 * @param[in,out] hashes                Payload deduplication entries tagged by block index.
 * @param[out]    incompatible          Incompatible byte masks for searched predecessors.
 *
 * @return The best first match, or an empty match if no trial improves the cost.
 */
static rdo_byte_match rdo_search_first_match(
	rdo_search_state& search,
	const uint8_t* data,
	const uint8_t mask[block_size],
	unsigned int search_start,
	int continuation_source,
	const uint32_t source_hashes[][128],
	uint32_t hashes[rdo_trial_hash_entries],
	uint32_t incompatible[brick_layout::block_count]
) {
	rdo_byte_match first_match;
	vint4 original_bytes = vint4::load(search.original);
	unsigned int block_tag = search.block_index;

	// Copy aligned byte ranges from finalized predecessors, trying long matches first
	for (unsigned int previous = search.block_index; previous > search_start;)
	{
		previous--;
		const uint8_t* source = data + previous * block_size;
		vint4 changed_bytes = vint4::load(source) ^ original_bytes;
		uint32_t invalid_bytes = rdo_get_incompatible_bytes(search.original, source, mask);
		incompatible[previous] = invalid_bytes;
		unsigned int distance = (search.block_index - previous) * block_size;
		// Intersect byte compatibility to enumerate only valid starts for each match length
		uint32_t compatible = (~invalid_bytes) & 0xFFFFu;
		uint32_t valid_offsets[17];
		uint32_t valid = compatible;
		for (unsigned int length = 1; length <= block_size; length++)
		{
			valid_offsets[length] = valid;
			valid &= compatible >> length;
		}

		for (unsigned int length = block_size; length >= min_match_length; length--)
		{
			uint32_t offsets = valid_offsets[length];
			if (offsets == 0)
			{
				continue;
			}

			float bits = search.rates.get_match_cost(distance, length);
			uint32_t match_mask = (1u << length) - 1;
			const uint32_t* substring_hashes = source_hashes[previous] + substring_offsets[length];

			for (; offsets != 0; offsets &= offsets - 1)
			{
				unsigned int offset = std::countr_zero(offsets);
				unsigned int source_offset = previous * block_size + offset;
				float trial_bits = bits;
				uint32_t hash = substring_hashes[offset];
				uint32_t hash_entry = (hash & 0xFFFFFF00) | block_tag;
				if ((static_cast<int>(source_offset) == continuation_source) && (offset == 0))
				{
					trial_bits = continuation_bits;
				}
				else
				{
					// Deduplicate payloads except continuations with a cheaper match price
					uint32_t entry = hashes[hash & (rdo_trial_hash_entries - 1)];
					if (entry == hash_entry)
					{
						continue;
					}
				}

				hashes[hash & (rdo_trial_hash_entries - 1)] = hash_entry;

				uint8_t trial[block_size];
				vint4 copied_bytes = vint4::load(prefix_masks[offset + length].data())
				                   ^ vint4::load(prefix_masks[offset].data());
				store(original_bytes ^ (changed_bytes & copied_bytes), trial);
				float rate = search.rates.estimate({ length, source_offset, offset, trial_bits }, {});
				if (rdo_evaluate_trial(search, trial, rate, match_mask << offset))
				{
					first_match = { length, source_offset, offset, trial_bits };
				}
			}
		}
	}

	return first_match;
}

/**
 * @brief Hold the best first match fixed and search for a second disjoint byte range.
 *
 * @param[in,out] search         The current block's scoring state and best result.
 * @param[in]     data           The compact brick stream with finalized predecessors.
 * @param         search_start   The first predecessor to search.
 * @param[in]     first_match    The best single copied range.
 * @param[in]     incompatible   Incompatible byte masks from the first search.
 *
 * @return The best second match, or an empty match if no trial improves the cost.
 */
static rdo_byte_match rdo_search_second_match(
	rdo_search_state& search,
	const uint8_t* data,
	unsigned int search_start,
	const rdo_byte_match& first_match,
	const uint32_t incompatible[brick_layout::block_count]
) {
	rdo_byte_match second_match;
	if ((!rdo_is_two_matches_enabled) || (first_match.length == 0)
	    || (first_match.length > block_size - min_match_length))
	{
		return second_match;
	}

	uint8_t first_block[block_size];
	std::memcpy(first_block, search.best_block, block_size);
	vint4 first_bytes = vint4::load(first_block);
	float first_rate = search.rates.estimate(first_match, {});
	unsigned int first_mask = ((1u << first_match.length) - 1) << first_match.offset;
	for (unsigned int previous = search.block_index; previous > search_start;)
	{
		previous--;
		const uint8_t* source = data + previous * block_size;
		vint4 changed_bytes = vint4::load(source) ^ first_bytes;
		uint32_t invalid_bytes = incompatible[previous];
		unsigned int distance = (search.block_index - previous) * block_size;
		for (unsigned int length = min_match_length;
		     length <= block_size - first_match.length; length++)
		{
			float bits = search.rates.get_match_cost(distance, length);
			uint32_t match_mask = (1u << length) - 1;

			for (unsigned int offset = 0; offset <= block_size - length; offset++)
			{
				if ((offset < first_match.offset + first_match.length)
				    && (offset + length > first_match.offset))
				{
					continue;
				}

				if ((invalid_bytes & (match_mask << offset)) != 0)
				{
					continue;
				}

				uint8_t trial[block_size];
				vint4 copied_bytes = vint4::load(prefix_masks[offset + length].data())
				                   ^ vint4::load(prefix_masks[offset].data());
				store(first_bytes ^ (changed_bytes & copied_bytes), trial);
				rdo_byte_match match { length, previous * block_size + offset, offset, bits };
				float trial_rate = search.rates.estimate(first_match, match);
				// Require meaningful estimated savings before decoding second matches
				if ((first_rate - trial_rate) < rdo_min_second_match_saving)
				{
					continue;
				}

				unsigned int byte_mask = first_mask | (match_mask << offset);
				if (rdo_evaluate_trial(search, trial, trial_rate, byte_mask))
				{
					second_match = { length, previous * block_size + offset, offset, 0.0f };
				}
			}
		}
	}

	return second_match;
}

/* See header for documentation. */
void rdo_reduce_entropy(
	rdo_brick_context& ctx,
	uint8_t* data,
	unsigned int block_count,
	const uint8_t* masks,
	const rdo_parameters& params
) {
	assert(block_count <= brick_layout::block_count);

	// Learn brick-wide token prices before optimizing the accepted prefix in stream order
	rdo_rate_model rates;
	rates.init(data, block_count, rdo_history_bytes / 16);

	constexpr unsigned int hash_size = rdo_trial_hash_entries;
	uint32_t hashes[hash_size] {};
	uint32_t source_hashes[brick_layout::block_count][128];
	uint32_t incompatible[brick_layout::block_count];
	int continuation_source = -1;

	for (unsigned int block_index = 0; block_index < block_count; block_index++)
	{
		if (block_index > 0)
		{
			// A predecessor's bytes remain final, so its substring hashes can serve every later block
			// Zero padding permits full-width prefix loads for the masked final hash batch
			uint32_t prefixes[16 + ASTCENC_SIMD_WIDTH] {};
			rdo_prepare_substring_hashes(data + (block_index - 1) * block_size, prefixes);
			for (unsigned int length = min_match_length; length <= block_size; length++)
			{
				rdo_hash_substrings(prefixes, length, source_hashes[block_index - 1] + substring_offsets[length]);
			}
		}

		uint8_t* original = data + block_index * block_size;

		// Refresh token prices once optimized history covers half of the brick
		if (rdo_is_midpoint_refresh_enabled && (block_index > 0) && (block_index == block_count / 2))
		{
			rates.init(data, block_count, rdo_history_bytes / 16);
		}

		// Keep constant blocks in byte history without spending work on their distortion or trials
		if ((original[0] == 0xFC) && ((original[1] & 1) != 0))
		{
			continuation_source = -1;
			continue;
		}

		const uint8_t* mask = masks + block_index * block_size;
		const auto statistics = rdo_initialize_block(ctx, block_index);

		// Penalize distortion more strongly in smooth regions where small errors are more visible
		float smooth = std::clamp(statistics.std_dev / rdo_max_smooth_std_dev, 0.0f, 1.0f);
		smooth *= smooth;
		float mse_scale = rdo_smooth_mse_scale + (1.0f - rdo_smooth_mse_scale) * smooth;
		float mse_recip = 1.0f / mse_scale;

		// Rate scoring retains the full window even when speculative copies search fewer blocks
		unsigned int history_count = std::min(block_index, rdo_history_bytes / 16);
		unsigned int history_start = block_index - history_count;
		unsigned int search_start = block_index - std::min(history_count, params.search_lookback);

		rates.prepare(data, block_index, history_start);
		float original_rate = rates.estimate({}, {});

		float original_cost = statistics.mse * mse_scale + original_rate * params.lambda;

		// Preserve the original RMS guard, including its error floor for exact blocks
		float max_mse = rdo_max_rms_increase * rdo_max_rms_increase
		              * std::max(statistics.mse, 1.0f);

		rdo_search_state search {
			.original = original,
			.params = params,
			.rates = rates,
			.ctx = ctx,
			.block_index = block_index,
			.original_rate = original_rate,
			.mse_scale = mse_scale,
			.mse_recip = mse_recip,
			.max_mse = max_mse,
			.best_cost = original_cost,
			.best_block = {}
		};

		std::memcpy(search.best_block, original, block_size);

		rdo_byte_match first_match = rdo_search_first_match(search, data, mask, search_start,
		                                                    continuation_source, source_hashes,
		                                                    hashes, incompatible);

		rdo_byte_match second_match = rdo_search_second_match(search, data, search_start,
		                                                      first_match, incompatible);

		// Preserve exact accepted bytes and track matches that can continue into the next block
		continuation_source = -1;
		if (search.best_cost < original_cost)
		{
			std::memcpy(original, search.best_block, block_size);
			rdo_byte_match last_match = first_match;
			if ((second_match.length != 0) && (second_match.offset > first_match.offset))
			{
				last_match = second_match;
			}

			if (last_match.offset + last_match.length == block_size)
			{
				continuation_source = static_cast<int>(last_match.source + last_match.length);
			}
		}
	}
}

#endif
