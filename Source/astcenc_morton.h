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
 * @brief Brick layouts for encoder scheduling and RDO.
 *
 * These dimensions are independent of the fixed 8x8 bricks used by the .astcm file format.
 */

#ifndef ASTCENC_MORTON_INCLUDED
#define ASTCENC_MORTON_INCLUDED

#include <bit>
#include <bitset>
#include <cstddef>
#include <cstdint>
#include <limits>

/** @brief A 3D position, in units specified by the caller. */
struct position
{
	/** @brief The X coordinate. */
	size_t x;

	/** @brief The Y coordinate. */
	size_t y;

	/** @brief The Z coordinate. */
	size_t z;
};

/**
 * @brief A square or 2:1 Morton brick with power-of-two dimensions.
 *
 * @tparam dim_x   The brick width in ASTC blocks.
 * @tparam dim_y   The brick height in ASTC blocks.
 */
template<unsigned int dim_x, unsigned int dim_y>
struct morton_brick_layout
{
	static_assert(std::has_single_bit(dim_x) && std::has_single_bit(dim_y),
	              "Morton brick dimensions must be non-zero powers of two");

	static_assert((dim_x == dim_y) || ((dim_x / 2 == dim_y) && (dim_x > dim_y)),
	              "Morton bricks must be 1:1 or 2:1 width-to-height ratio");

	static_assert(dim_x <= 16,
	              "Morton brick must be less than or equal to 16 blocks");

	static_assert(dim_y <= 16,
	              "Morton brick must be less than or equal to 16 blocks");

	/** @brief The brick width in ASTC blocks. */
	static constexpr unsigned int width = dim_x;

	/** @brief The brick height in ASTC blocks. */
	static constexpr unsigned int height = dim_y;

	/** @brief The number of Morton positions in the brick. */
	static constexpr unsigned int block_count = dim_x * dim_y;

	/** @brief A validity bit for each Morton position, including sparse edge bricks. */
	using validity_mask = std::bitset<block_count>;

	/**
	 * @brief Decode a dense Morton index into local block coordinates.
	 *
	 * A 2:1 brick appends a high X bit after the interleaved square coordinates.
	 *
	 * @param index   The Morton index in the range [0, block_count).
	 *
	 * @return The local block coordinates, with Z set to zero.
	 */
	static constexpr position get_position(
		unsigned int index
	) {
		unsigned int x = 0;
		unsigned int y = 0;
		for (unsigned int bit = 1; bit < height; bit <<= 1)
		{
			x |= (index & 1) * bit;
			y |= ((index >> 1) & 1) * bit;
			index >>= 2;
		}

		if constexpr (width > height)
		{
			x |= (index & 1) * height;
		}

		return { .x = x, .y = y, .z = 0 };
	}
};

/** @brief The selected encoder and RDO brick layout, independent of .astcm file tiles. */
struct brick_layout
{
	/** @brief The number of block positions in either layout. */
	static constexpr unsigned int block_count = 128;

	/** @brief A validity bit for each position, including sparse edge bricks. */
	using validity_mask = std::bitset<block_count>;

	/** @brief The brick width in ASTC blocks. */
	unsigned int width { 128 };

	/** @brief The brick height in ASTC blocks. */
	unsigned int height { 1 };

	/**
	 * @brief Select the brick dimensions and traversal order.
	 *
	 * @param is_morton   Whether to use Morton order instead of scanline order.
	 */
	explicit brick_layout(bool is_morton)
	{
		if (is_morton)
		{
			width = 16;
			height = 8;
		}
	}

	/**
	 * @brief Decode a brick index into local block coordinates in the selected order.
	 *
	 * @param index   The block index in the range [0, block_count).
	 *
	 * @return The local block coordinates, with Z set to zero.
	 */
	position get_position(unsigned int index) const
	{
		if (height == 1)
		{
			return { .x = index, .y = 0, .z = 0 };
		}

		return morton_brick_layout<16, 8>::get_position(index);
	}
};

#endif
