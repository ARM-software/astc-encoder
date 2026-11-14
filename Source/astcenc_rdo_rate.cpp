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
 * @brief DEFLATE-aware token costs for the RDO entropy reduction transform.
 *
 * Token alphabets and extra-bit counts follow RFC 1951, sections 3.2.5 and 3.2.7.
 */

#if !defined(ASTCENC_DECOMPRESS_ONLY)

#include "astcenc_rdo.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <limits>
#include <queue>
#include <vector>

/** @brief The base length represented by each DEFLATE length symbol. */
static constexpr unsigned int length_base[29] {
	3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27,
	31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258
};

/** @brief Extra bits following each DEFLATE length symbol. */
static constexpr unsigned int length_extra[29] {
	0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2,
	2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0
};

/** @brief The base distance represented by each DEFLATE distance symbol. */
static constexpr unsigned int distance_base[30] {
	1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129,
	193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097,
	6145, 8193, 12289, 16385, 24577
};

/** @brief Extra bits following each DEFLATE distance symbol. */
static constexpr unsigned int distance_extra[30] {
	0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6,
	6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13
};

/**
 * @brief Find the DEFLATE length symbol index.
 *
 * @param length   The match length, in the range 3..258.
 *
 * @return The length symbol index, excluding the literal alphabet and end marker.
 */
static unsigned int rdo_get_length_code(
	unsigned int length
) {
	unsigned int code = 0;
	while (((code + 1) < 29) && (length >= length_base[code + 1]))
	{
		code++;
	}

	return code;
}

/**
 * @brief Find the DEFLATE distance symbol index.
 *
 * @param distance   The backward distance, in the range 1..32768.
 *
 * @return The distance symbol index.
 */
static unsigned int rdo_get_distance_code(
	unsigned int distance
) {
	unsigned int code = 0;
	while (((code + 1) < 30) && (distance >= distance_base[code + 1]))
	{
		code++;
	}

	return code;
}

/**
 * @brief Hash three bytes into the local match dictionary.
 *
 * @param[in] data   The three source bytes.
 *
 * @return The hash bucket index.
 */
static unsigned int rdo_hash_bytes(
	const uint8_t* data
) {
	unsigned int value = static_cast<unsigned int>(data[0]) << 16;
	value |= static_cast<unsigned int>(data[1]) << 8;
	value |= data[2];
	return (value * 0x1E35A7BDu) >> (32 - rdo_rate_hash_bits);
}

/**
 * @brief Build deterministic Huffman code lengths with a small unseen-symbol prior.
 *
 * Large bricks can exceed DEFLATE's 15-bit limit; callers then use fixed code lengths.
 *
 * @param[in]  frequencies   Token frequencies, with at most 286 entries.
 * @param      count         The alphabet size.
 * @param[out] lengths       The resulting code lengths.
 *
 * @return True if every code length fits DEFLATE's limit.
 */
static bool rdo_build_code_lengths(
	const unsigned int* frequencies,
	unsigned int count,
	uint8_t* lengths
) {
	unsigned int weights[571];
	unsigned int parents[571];
	std::fill_n(parents, 2 * count - 1, ~0u);
	for (unsigned int i = 0; i < count; i++)
	{
		weights[i] = frequencies[i] * 4 + 1;
	}

	// Break equal-frequency ties by symbol index to keep prices invariant across platforms
	auto compare = [&](unsigned int a, unsigned int b)
	{
		if (weights[a] != weights[b])
		{
			return weights[a] > weights[b];
		}

		return a > b;
	};

	std::priority_queue<unsigned int, std::vector<unsigned int>, decltype(compare)> heap(compare);
	for (unsigned int i = 0; i < count; i++)
	{
		heap.push(i);
	}

	unsigned int next = count;
	while (heap.size() > 1)
	{
		unsigned int a = heap.top();
		heap.pop();
		unsigned int b = heap.top();
		heap.pop();
		weights[next] = weights[a] + weights[b];
		parents[a] = next;
		parents[b] = next;
		heap.push(next++);
	}

	for (unsigned int i = 0; i < count; i++)
	{
		unsigned int depth = 0;
		for (unsigned int node = i; parents[node] != ~0u; node = parents[node])
		{
			depth++;
		}

		if (depth > 15)
		{
			return false;
		}

		lengths[i] = static_cast<uint8_t>(depth);
	}

	return true;
}

/* See header for documentation. */
void rdo_rate_model::init(
	const uint8_t* data,
	unsigned int block_count,
	unsigned int lookback
) {
	assert(block_count <= brick_layout::block_count);
	unsigned int byte_count = block_count * 16;
	history_bytes = std::min(lookback, block_count) * 16;
	unsigned int literal_length_frequencies[286] {};
	unsigned int distance_frequencies[30] {};

	std::fill_n(heads, hash_size, -1);

	// Train on emitted tokens rather than raw byte frequencies
	for (unsigned int position = 0; position < byte_count;)
	{
		unsigned int longest = 0;
		unsigned int best_distance = 0;
		unsigned int history_start = 0;

		if (position > history_bytes)
		{
			history_start = position - history_bytes;
		}

		unsigned int limit = std::min(258u, byte_count - position);
		if (limit >= 3)
		{
			int reference = heads[rdo_hash_bytes(data + position)];

			while (reference >= static_cast<int>(history_start))
			{
				unsigned int source = static_cast<unsigned int>(reference);
				unsigned int length = 0;

				while ((length < limit) && (data[source + length] == data[position + length]))
				{
					length++;
				}

				if (length > longest)
				{
					longest = length;
					best_distance = position - source;
				}

				reference = previous[source];
			}
		}

		unsigned int consumed = 1;
		if (longest >= 3)
		{
			literal_length_frequencies[257 + rdo_get_length_code(longest)]++;
			distance_frequencies[rdo_get_distance_code(best_distance)]++;
			consumed = longest;
		}
		else
		{
			literal_length_frequencies[data[position]]++;
		}

		unsigned int end = position + consumed;
		for (unsigned int source = position; (source < end) && ((source + 2) < byte_count); source++)
		{
			unsigned int hash = rdo_hash_bytes(data + source);
			previous[source] = heads[hash];
			heads[hash] = static_cast<int>(source);
		}

		position = end;
	}

	// Include the end marker and price unseen symbols so new trial tokens remain representable
	literal_length_frequencies[256] = 1;
	bool is_literal_tree_valid = rdo_build_code_lengths(literal_length_frequencies, 286,
	                                                    literal_length_bits);
	if (!is_literal_tree_valid)
	{
		for (unsigned int symbol = 0; symbol < 286; symbol++)
		{
			unsigned int bits = 8;
			if ((symbol >= 144) && (symbol <= 255))
			{
				bits = 9;
			}
			else if ((symbol >= 256) && (symbol <= 279))
			{
				bits = 7;
			}

			literal_length_bits[symbol] = static_cast<uint8_t>(bits);
		}
	}

	uint8_t distance_bits[30];
	bool is_distance_tree_valid = rdo_build_code_lengths(distance_frequencies, 30, distance_bits);
	if (!is_distance_tree_valid)
	{
		std::fill_n(distance_bits, 30, 5);
	}

	// Fold Huffman lengths and DEFLATE extra bits into tables used by the inner search
	for (unsigned int length = 3; length <= 16; length++)
	{
		unsigned int code = rdo_get_length_code(length);
		length_costs[length] = static_cast<uint8_t>(literal_length_bits[257 + code] + length_extra[code]);
	}

	for (unsigned int distance = 1; distance <= max_bytes; distance++)
	{
		unsigned int code = rdo_get_distance_code(distance);
		distance_costs[distance] = static_cast<uint8_t>(distance_bits[code] + distance_extra[code]);
	}

	// Rebuild the dictionary incrementally from accepted bytes rather than the training stream
	std::fill_n(heads, hash_size, -1);
	indexed_bytes = 0;
}

/* See header for documentation. */
void rdo_rate_model::prepare(
	const uint8_t* data,
	unsigned int block_index,
	unsigned int history_start
) {
	unsigned int block_start = block_index * 16;
	unsigned int history_byte_start = history_start * 16;
	std::memcpy(original, data + block_start, 16);

	// Enter only completed source substrings; accepted prefix bytes no longer change
	while ((indexed_bytes + 2) < block_start)
	{
		unsigned int hash = rdo_hash_bytes(data + indexed_bytes);
		previous[indexed_bytes] = heads[hash];
		heads[hash] = static_cast<int>(indexed_bytes);
		indexed_bytes++;
	}

	// Find existing unaligned matches so unchanged bytes are not overcharged as literals
	float matches[16][17];
	for (unsigned int start = 0; start < 16; start++)
	{
		std::fill_n(matches[start], 17, std::numeric_limits<float>::infinity());
		unsigned int limit = 16 - start;
		if (limit < 3)
		{
			continue;
		}

		// Advance expiry per byte, including unaligned matches in the oldest history block
		unsigned int first_source = history_byte_start;
		if ((block_start + start) > history_bytes)
		{
			first_source = std::max(first_source, block_start + start - history_bytes);
		}

		int reference = heads[rdo_hash_bytes(original + start)];
		while (reference >= static_cast<int>(first_source))
		{
			unsigned int source = static_cast<unsigned int>(reference);
			unsigned int max_length = std::min(limit, block_start - source);
			unsigned int length = 0;

			while ((length < max_length) && (data[source + length] == original[start + length]))
			{
				length++;
			}

			unsigned int distance = block_start + start - source;
			for (unsigned int count = 3; count <= length; count++)
			{
				matches[start][count] = std::min(matches[start][count], get_match_cost(distance, count));
			}

			reference = previous[source];
		}
	}

	// Cache the cheapest parse of every interval that a trial can leave untouched
	for (unsigned int start = 0; start <= 16; start++)
	{
		intervals[start][start] = 0.0f;
		for (unsigned int end = start + 1; end <= 16; end++)
		{
			float best = intervals[start][end - 1] + get_literal_cost(original[end - 1]);
			for (unsigned int length = 3; length <= (end - start); length++)
			{
				best = std::min(best, intervals[start][end - length] + matches[end - length][length]);
			}

			intervals[start][end] = best;
		}
	}
}

/* See header for documentation. */
float rdo_rate_model::get_match_cost(
	unsigned int distance,
	unsigned int length
) const {
	assert((distance > 0) && (distance <= max_bytes));
	assert((length >= 3) && (length <= 16));
	return static_cast<float>(length_costs[length] + distance_costs[distance]);
}

/* See header for documentation. */
float rdo_rate_model::get_literal_cost(
	uint8_t value
) const {
	return static_cast<float>(literal_length_bits[value]);
}

/* See header for documentation. */
float rdo_rate_model::estimate(
	const rdo_byte_match& first,
	const rdo_byte_match& second
) const {
	if (first.length == 0)
	{
		return intervals[0][16];
	}

	// Splice forced matches between cached optimal parses of the untouched intervals
	const rdo_byte_match* left = &first;
	const rdo_byte_match* right = &second;
	if ((second.length > 0) && (second.offset < first.offset))
	{
		std::swap(left, right);
	}

	unsigned int end = left->offset + left->length;
	float bits = intervals[0][left->offset] + left->bits;
	if (second.length > 0)
	{
		bits += intervals[end][right->offset] + right->bits;
		end = right->offset + right->length;
	}

	bits += intervals[end][16];
	return bits;
}

#endif
