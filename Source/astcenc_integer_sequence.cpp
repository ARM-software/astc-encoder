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
 * @brief Functions for encoding/decoding Bounded Integer Sequence Encoding.
 */

#include "astcenc_internal.h"

#include <array>

/** @brief Unpacked quint triplets <low,middle,high> for each packed value. */
static constexpr uint8_t quints_of_integer[128][3] {
	{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0},
	{4, 0, 0}, {0, 4, 0}, {4, 4, 0}, {4, 4, 4},
	{0, 1, 0}, {1, 1, 0}, {2, 1, 0}, {3, 1, 0},
	{4, 1, 0}, {1, 4, 0}, {4, 4, 1}, {4, 4, 4},
	{0, 2, 0}, {1, 2, 0}, {2, 2, 0}, {3, 2, 0},
	{4, 2, 0}, {2, 4, 0}, {4, 4, 2}, {4, 4, 4},
	{0, 3, 0}, {1, 3, 0}, {2, 3, 0}, {3, 3, 0},
	{4, 3, 0}, {3, 4, 0}, {4, 4, 3}, {4, 4, 4},
	{0, 0, 1}, {1, 0, 1}, {2, 0, 1}, {3, 0, 1},
	{4, 0, 1}, {0, 4, 1}, {4, 0, 4}, {0, 4, 4},
	{0, 1, 1}, {1, 1, 1}, {2, 1, 1}, {3, 1, 1},
	{4, 1, 1}, {1, 4, 1}, {4, 1, 4}, {1, 4, 4},
	{0, 2, 1}, {1, 2, 1}, {2, 2, 1}, {3, 2, 1},
	{4, 2, 1}, {2, 4, 1}, {4, 2, 4}, {2, 4, 4},
	{0, 3, 1}, {1, 3, 1}, {2, 3, 1}, {3, 3, 1},
	{4, 3, 1}, {3, 4, 1}, {4, 3, 4}, {3, 4, 4},
	{0, 0, 2}, {1, 0, 2}, {2, 0, 2}, {3, 0, 2},
	{4, 0, 2}, {0, 4, 2}, {2, 0, 4}, {3, 0, 4},
	{0, 1, 2}, {1, 1, 2}, {2, 1, 2}, {3, 1, 2},
	{4, 1, 2}, {1, 4, 2}, {2, 1, 4}, {3, 1, 4},
	{0, 2, 2}, {1, 2, 2}, {2, 2, 2}, {3, 2, 2},
	{4, 2, 2}, {2, 4, 2}, {2, 2, 4}, {3, 2, 4},
	{0, 3, 2}, {1, 3, 2}, {2, 3, 2}, {3, 3, 2},
	{4, 3, 2}, {3, 4, 2}, {2, 3, 4}, {3, 3, 4},
	{0, 0, 3}, {1, 0, 3}, {2, 0, 3}, {3, 0, 3},
	{4, 0, 3}, {0, 4, 3}, {0, 0, 4}, {1, 0, 4},
	{0, 1, 3}, {1, 1, 3}, {2, 1, 3}, {3, 1, 3},
	{4, 1, 3}, {1, 4, 3}, {0, 1, 4}, {1, 1, 4},
	{0, 2, 3}, {1, 2, 3}, {2, 2, 3}, {3, 2, 3},
	{4, 2, 3}, {2, 4, 3}, {0, 2, 4}, {1, 2, 4},
	{0, 3, 3}, {1, 3, 3}, {2, 3, 3}, {3, 3, 3},
	{4, 3, 3}, {3, 4, 3}, {0, 3, 4}, {1, 3, 4}
};

/** @brief Packed quint values for each unpacked value, indexed [hi][mid][lo]. */
static const uint8_t integer_of_quints[5][5][5] {
	{
		{0, 1, 2, 3, 4},
		{8, 9, 10, 11, 12},
		{16, 17, 18, 19, 20},
		{24, 25, 26, 27, 28},
		{5, 13, 21, 29, 6}
	},
	{
		{32, 33, 34, 35, 36},
		{40, 41, 42, 43, 44},
		{48, 49, 50, 51, 52},
		{56, 57, 58, 59, 60},
		{37, 45, 53, 61, 14}
	},
	{
		{64, 65, 66, 67, 68},
		{72, 73, 74, 75, 76},
		{80, 81, 82, 83, 84},
		{88, 89, 90, 91, 92},
		{69, 77, 85, 93, 22}
	},
	{
		{96, 97, 98, 99, 100},
		{104, 105, 106, 107, 108},
		{112, 113, 114, 115, 116},
		{120, 121, 122, 123, 124},
		{101, 109, 117, 125, 30}
	},
	{
		{102, 103, 70, 71, 38},
		{110, 111, 78, 79, 46},
		{118, 119, 86, 87, 54},
		{126, 127, 94, 95, 62},
		{39, 47, 55, 63, 31}
	}
};

/** @brief Unpacked trit quintuplets <low,...,high> for each packed value. */
static constexpr uint8_t trits_of_integer[256][5] {
	{0, 0, 0, 0, 0}, {1, 0, 0, 0, 0}, {2, 0, 0, 0, 0}, {0, 0, 2, 0, 0},
	{0, 1, 0, 0, 0}, {1, 1, 0, 0, 0}, {2, 1, 0, 0, 0}, {1, 0, 2, 0, 0},
	{0, 2, 0, 0, 0}, {1, 2, 0, 0, 0}, {2, 2, 0, 0, 0}, {2, 0, 2, 0, 0},
	{0, 2, 2, 0, 0}, {1, 2, 2, 0, 0}, {2, 2, 2, 0, 0}, {2, 0, 2, 0, 0},
	{0, 0, 1, 0, 0}, {1, 0, 1, 0, 0}, {2, 0, 1, 0, 0}, {0, 1, 2, 0, 0},
	{0, 1, 1, 0, 0}, {1, 1, 1, 0, 0}, {2, 1, 1, 0, 0}, {1, 1, 2, 0, 0},
	{0, 2, 1, 0, 0}, {1, 2, 1, 0, 0}, {2, 2, 1, 0, 0}, {2, 1, 2, 0, 0},
	{0, 0, 0, 2, 2}, {1, 0, 0, 2, 2}, {2, 0, 0, 2, 2}, {0, 0, 2, 2, 2},
	{0, 0, 0, 1, 0}, {1, 0, 0, 1, 0}, {2, 0, 0, 1, 0}, {0, 0, 2, 1, 0},
	{0, 1, 0, 1, 0}, {1, 1, 0, 1, 0}, {2, 1, 0, 1, 0}, {1, 0, 2, 1, 0},
	{0, 2, 0, 1, 0}, {1, 2, 0, 1, 0}, {2, 2, 0, 1, 0}, {2, 0, 2, 1, 0},
	{0, 2, 2, 1, 0}, {1, 2, 2, 1, 0}, {2, 2, 2, 1, 0}, {2, 0, 2, 1, 0},
	{0, 0, 1, 1, 0}, {1, 0, 1, 1, 0}, {2, 0, 1, 1, 0}, {0, 1, 2, 1, 0},
	{0, 1, 1, 1, 0}, {1, 1, 1, 1, 0}, {2, 1, 1, 1, 0}, {1, 1, 2, 1, 0},
	{0, 2, 1, 1, 0}, {1, 2, 1, 1, 0}, {2, 2, 1, 1, 0}, {2, 1, 2, 1, 0},
	{0, 1, 0, 2, 2}, {1, 1, 0, 2, 2}, {2, 1, 0, 2, 2}, {1, 0, 2, 2, 2},
	{0, 0, 0, 2, 0}, {1, 0, 0, 2, 0}, {2, 0, 0, 2, 0}, {0, 0, 2, 2, 0},
	{0, 1, 0, 2, 0}, {1, 1, 0, 2, 0}, {2, 1, 0, 2, 0}, {1, 0, 2, 2, 0},
	{0, 2, 0, 2, 0}, {1, 2, 0, 2, 0}, {2, 2, 0, 2, 0}, {2, 0, 2, 2, 0},
	{0, 2, 2, 2, 0}, {1, 2, 2, 2, 0}, {2, 2, 2, 2, 0}, {2, 0, 2, 2, 0},
	{0, 0, 1, 2, 0}, {1, 0, 1, 2, 0}, {2, 0, 1, 2, 0}, {0, 1, 2, 2, 0},
	{0, 1, 1, 2, 0}, {1, 1, 1, 2, 0}, {2, 1, 1, 2, 0}, {1, 1, 2, 2, 0},
	{0, 2, 1, 2, 0}, {1, 2, 1, 2, 0}, {2, 2, 1, 2, 0}, {2, 1, 2, 2, 0},
	{0, 2, 0, 2, 2}, {1, 2, 0, 2, 2}, {2, 2, 0, 2, 2}, {2, 0, 2, 2, 2},
	{0, 0, 0, 0, 2}, {1, 0, 0, 0, 2}, {2, 0, 0, 0, 2}, {0, 0, 2, 0, 2},
	{0, 1, 0, 0, 2}, {1, 1, 0, 0, 2}, {2, 1, 0, 0, 2}, {1, 0, 2, 0, 2},
	{0, 2, 0, 0, 2}, {1, 2, 0, 0, 2}, {2, 2, 0, 0, 2}, {2, 0, 2, 0, 2},
	{0, 2, 2, 0, 2}, {1, 2, 2, 0, 2}, {2, 2, 2, 0, 2}, {2, 0, 2, 0, 2},
	{0, 0, 1, 0, 2}, {1, 0, 1, 0, 2}, {2, 0, 1, 0, 2}, {0, 1, 2, 0, 2},
	{0, 1, 1, 0, 2}, {1, 1, 1, 0, 2}, {2, 1, 1, 0, 2}, {1, 1, 2, 0, 2},
	{0, 2, 1, 0, 2}, {1, 2, 1, 0, 2}, {2, 2, 1, 0, 2}, {2, 1, 2, 0, 2},
	{0, 2, 2, 2, 2}, {1, 2, 2, 2, 2}, {2, 2, 2, 2, 2}, {2, 0, 2, 2, 2},
	{0, 0, 0, 0, 1}, {1, 0, 0, 0, 1}, {2, 0, 0, 0, 1}, {0, 0, 2, 0, 1},
	{0, 1, 0, 0, 1}, {1, 1, 0, 0, 1}, {2, 1, 0, 0, 1}, {1, 0, 2, 0, 1},
	{0, 2, 0, 0, 1}, {1, 2, 0, 0, 1}, {2, 2, 0, 0, 1}, {2, 0, 2, 0, 1},
	{0, 2, 2, 0, 1}, {1, 2, 2, 0, 1}, {2, 2, 2, 0, 1}, {2, 0, 2, 0, 1},
	{0, 0, 1, 0, 1}, {1, 0, 1, 0, 1}, {2, 0, 1, 0, 1}, {0, 1, 2, 0, 1},
	{0, 1, 1, 0, 1}, {1, 1, 1, 0, 1}, {2, 1, 1, 0, 1}, {1, 1, 2, 0, 1},
	{0, 2, 1, 0, 1}, {1, 2, 1, 0, 1}, {2, 2, 1, 0, 1}, {2, 1, 2, 0, 1},
	{0, 0, 1, 2, 2}, {1, 0, 1, 2, 2}, {2, 0, 1, 2, 2}, {0, 1, 2, 2, 2},
	{0, 0, 0, 1, 1}, {1, 0, 0, 1, 1}, {2, 0, 0, 1, 1}, {0, 0, 2, 1, 1},
	{0, 1, 0, 1, 1}, {1, 1, 0, 1, 1}, {2, 1, 0, 1, 1}, {1, 0, 2, 1, 1},
	{0, 2, 0, 1, 1}, {1, 2, 0, 1, 1}, {2, 2, 0, 1, 1}, {2, 0, 2, 1, 1},
	{0, 2, 2, 1, 1}, {1, 2, 2, 1, 1}, {2, 2, 2, 1, 1}, {2, 0, 2, 1, 1},
	{0, 0, 1, 1, 1}, {1, 0, 1, 1, 1}, {2, 0, 1, 1, 1}, {0, 1, 2, 1, 1},
	{0, 1, 1, 1, 1}, {1, 1, 1, 1, 1}, {2, 1, 1, 1, 1}, {1, 1, 2, 1, 1},
	{0, 2, 1, 1, 1}, {1, 2, 1, 1, 1}, {2, 2, 1, 1, 1}, {2, 1, 2, 1, 1},
	{0, 1, 1, 2, 2}, {1, 1, 1, 2, 2}, {2, 1, 1, 2, 2}, {1, 1, 2, 2, 2},
	{0, 0, 0, 2, 1}, {1, 0, 0, 2, 1}, {2, 0, 0, 2, 1}, {0, 0, 2, 2, 1},
	{0, 1, 0, 2, 1}, {1, 1, 0, 2, 1}, {2, 1, 0, 2, 1}, {1, 0, 2, 2, 1},
	{0, 2, 0, 2, 1}, {1, 2, 0, 2, 1}, {2, 2, 0, 2, 1}, {2, 0, 2, 2, 1},
	{0, 2, 2, 2, 1}, {1, 2, 2, 2, 1}, {2, 2, 2, 2, 1}, {2, 0, 2, 2, 1},
	{0, 0, 1, 2, 1}, {1, 0, 1, 2, 1}, {2, 0, 1, 2, 1}, {0, 1, 2, 2, 1},
	{0, 1, 1, 2, 1}, {1, 1, 1, 2, 1}, {2, 1, 1, 2, 1}, {1, 1, 2, 2, 1},
	{0, 2, 1, 2, 1}, {1, 2, 1, 2, 1}, {2, 2, 1, 2, 1}, {2, 1, 2, 2, 1},
	{0, 2, 1, 2, 2}, {1, 2, 1, 2, 2}, {2, 2, 1, 2, 2}, {2, 1, 2, 2, 2},
	{0, 0, 0, 1, 2}, {1, 0, 0, 1, 2}, {2, 0, 0, 1, 2}, {0, 0, 2, 1, 2},
	{0, 1, 0, 1, 2}, {1, 1, 0, 1, 2}, {2, 1, 0, 1, 2}, {1, 0, 2, 1, 2},
	{0, 2, 0, 1, 2}, {1, 2, 0, 1, 2}, {2, 2, 0, 1, 2}, {2, 0, 2, 1, 2},
	{0, 2, 2, 1, 2}, {1, 2, 2, 1, 2}, {2, 2, 2, 1, 2}, {2, 0, 2, 1, 2},
	{0, 0, 1, 1, 2}, {1, 0, 1, 1, 2}, {2, 0, 1, 1, 2}, {0, 1, 2, 1, 2},
	{0, 1, 1, 1, 2}, {1, 1, 1, 1, 2}, {2, 1, 1, 1, 2}, {1, 1, 2, 1, 2},
	{0, 2, 1, 1, 2}, {1, 2, 1, 1, 2}, {2, 2, 1, 1, 2}, {2, 1, 2, 1, 2},
	{0, 2, 2, 2, 2}, {1, 2, 2, 2, 2}, {2, 2, 2, 2, 2}, {2, 1, 2, 2, 2}
};

/** @brief Packed trit values for each unpacked value, indexed [hi][][][][lo]. */
static const uint8_t integer_of_trits[3][3][3][3][3] {
	{
		{
			{
				{0, 1, 2},
				{4, 5, 6},
				{8, 9, 10}
			},
			{
				{16, 17, 18},
				{20, 21, 22},
				{24, 25, 26}
			},
			{
				{3, 7, 15},
				{19, 23, 27},
				{12, 13, 14}
			}
		},
		{
			{
				{32, 33, 34},
				{36, 37, 38},
				{40, 41, 42}
			},
			{
				{48, 49, 50},
				{52, 53, 54},
				{56, 57, 58}
			},
			{
				{35, 39, 47},
				{51, 55, 59},
				{44, 45, 46}
			}
		},
		{
			{
				{64, 65, 66},
				{68, 69, 70},
				{72, 73, 74}
			},
			{
				{80, 81, 82},
				{84, 85, 86},
				{88, 89, 90}
			},
			{
				{67, 71, 79},
				{83, 87, 91},
				{76, 77, 78}
			}
		}
	},
	{
		{
			{
				{128, 129, 130},
				{132, 133, 134},
				{136, 137, 138}
			},
			{
				{144, 145, 146},
				{148, 149, 150},
				{152, 153, 154}
			},
			{
				{131, 135, 143},
				{147, 151, 155},
				{140, 141, 142}
			}
		},
		{
			{
				{160, 161, 162},
				{164, 165, 166},
				{168, 169, 170}
			},
			{
				{176, 177, 178},
				{180, 181, 182},
				{184, 185, 186}
			},
			{
				{163, 167, 175},
				{179, 183, 187},
				{172, 173, 174}
			}
		},
		{
			{
				{192, 193, 194},
				{196, 197, 198},
				{200, 201, 202}
			},
			{
				{208, 209, 210},
				{212, 213, 214},
				{216, 217, 218}
			},
			{
				{195, 199, 207},
				{211, 215, 219},
				{204, 205, 206}
			}
		}
	},
	{
		{
			{
				{96, 97, 98},
				{100, 101, 102},
				{104, 105, 106}
			},
			{
				{112, 113, 114},
				{116, 117, 118},
				{120, 121, 122}
			},
			{
				{99, 103, 111},
				{115, 119, 123},
				{108, 109, 110}
			}
		},
		{
			{
				{224, 225, 226},
				{228, 229, 230},
				{232, 233, 234}
			},
			{
				{240, 241, 242},
				{244, 245, 246},
				{248, 249, 250}
			},
			{
				{227, 231, 239},
				{243, 247, 251},
				{236, 237, 238}
			}
		},
		{
			{
				{28, 29, 30},
				{60, 61, 62},
				{92, 93, 94}
			},
			{
				{156, 157, 158},
				{188, 189, 190},
				{220, 221, 222}
			},
			{
				{31, 63, 127},
				{159, 191, 255},
				{252, 253, 254}
			}
		}
	}
};

/**
 * @brief The number of bits, trits, and quints needed for a quant level.
 */
struct btq_count
{
	/** @brief The number of bits. */
	uint8_t bits: 6;

	/** @brief The number of trits. */
	uint8_t trits: 1;

	/** @brief The number of quints. */
	uint8_t quints: 1;
};

/**
 * @brief The table of bits, trits, and quints needed for a quant encode.
 */
static constexpr std::array<btq_count, 21> btq_counts {{
	{ 1, 0, 0 },  // QUANT_2
	{ 0, 1, 0 },  // QUANT_3
	{ 2, 0, 0 },  // QUANT_4
	{ 0, 0, 1 },  // QUANT_5
	{ 1, 1, 0 },  // QUANT_6
	{ 3, 0, 0 },  // QUANT_8
	{ 1, 0, 1 },  // QUANT_10
	{ 2, 1, 0 },  // QUANT_12
	{ 4, 0, 0 },  // QUANT_16
	{ 2, 0, 1 },  // QUANT_20
	{ 3, 1, 0 },  // QUANT_24
	{ 5, 0, 0 },  // QUANT_32
	{ 3, 0, 1 },  // QUANT_40
	{ 4, 1, 0 },  // QUANT_48
	{ 6, 0, 0 },  // QUANT_64
	{ 4, 0, 1 },  // QUANT_80
	{ 5, 1, 0 },  // QUANT_96
	{ 7, 0, 0 },  // QUANT_128
	{ 5, 0, 1 },  // QUANT_160
	{ 6, 1, 0 },  // QUANT_192
	{ 8, 0, 0 }   // QUANT_256
}};

/**
 * @brief The sequence scale, round, and divisors needed to compute sizing.
 *
 * The length of a quantized sequence in bits is:
 *     (scale * <sequence_len> + round) / divisor
 */
struct ise_size
{
	/** @brief The scaling parameter. */
	uint8_t scale:6;

	/** @brief The divisor parameter. */
	uint8_t divisor:2;
};

/**
 * @brief The table of scale, round, and divisors needed for quant sizing.
 */
static const std::array<ise_size, 21> ise_sizes {{
	{  1, 0 },  // QUANT_2
	{  8, 2 },  // QUANT_3
	{  2, 0 },  // QUANT_4
	{  7, 1 },  // QUANT_5
	{ 13, 2 },  // QUANT_6
	{  3, 0 },  // QUANT_8
	{ 10, 1 },  // QUANT_10
	{ 18, 2 },  // QUANT_12
	{  4, 0 },  // QUANT_16
	{ 13, 1 },  // QUANT_20
	{ 23, 2 },  // QUANT_24
	{  5, 0 },  // QUANT_32
	{ 16, 1 },  // QUANT_40
	{ 28, 2 },  // QUANT_48
	{  6, 0 },  // QUANT_64
	{ 19, 1 },  // QUANT_80
	{ 33, 2 },  // QUANT_96
	{  7, 0 },  // QUANT_128
	{ 22, 1 },  // QUANT_160
	{ 38, 2 },  // QUANT_192
	{  8, 0 }   // QUANT_256
}};

/* See header for documentation. */
unsigned int get_ise_sequence_bitcount(
	unsigned int character_count,
	quant_method quant_level
) {
	// Cope with out-of bounds values - input might be invalid
	if (static_cast<size_t>(quant_level) >= ise_sizes.size())
	{
		// Arbitrary large number that's more than an ASTC block can hold
		return 1024;
	}

	auto& entry = ise_sizes[quant_level];
	unsigned int divisor = (entry.divisor << 1) + 1;
	return (entry.scale * character_count + divisor - 1) / divisor;
}

/**
 * @brief Write up to 8 bits at an arbitrary bit offset.
 *
 * The stored value is at most 8 bits, but can be stored at an offset of between 0 and 7 bits so may
 * span two separate bytes in memory.
 *
 * @param         value       The value to write.
 * @param         bitcount    The number of bits to write, starting from LSB.
 * @param         bitoffset   The bit offset to store at, between 0 and 7.
 * @param[in,out] ptr         The data pointer to write to.
 */
static inline void write_bits(
	unsigned int value,
	unsigned int bitcount,
	unsigned int bitoffset,
	uint8_t ptr[2]
) {
	unsigned int mask = (1 << bitcount) - 1;
	value &= mask;
	ptr += bitoffset >> 3;
	bitoffset &= 7;
	value <<= bitoffset;
	mask <<= bitoffset;
	mask = ~mask;

	ptr[0] &= mask;
	ptr[0] |= value;
	ptr[1] &= mask >> 8;
	ptr[1] |= value >> 8;
}

/**
 * @brief Read up to 16 bits from two bytes.
 *
 * This function reads a packed N-bit field from two bytes in memory. The stored value must exist
 * within the two bytes, but can start at an arbitrary bit offset and span the two bytes in memory.
 *
 * @param bitcount    The number of bits to read.
 * @param bitoffset   The bit offset to read from, between 0 and 7.
 * @param ptr         The data pointer to read from.
 *
 * @return The read value.
 */
[[maybe_unused]] static inline unsigned int read_bits(
	unsigned int bitcount,
	unsigned int bitoffset,
	const uint8_t* ptr
) {
	unsigned int mask = (1 << bitcount) - 1;
	ptr += bitoffset >> 3;
	bitoffset &= 7;
	unsigned int value = ptr[0] | (ptr[1] << 8);
	value >>= bitoffset;
	value &= mask;
	return value;
}

/* See header for documentation. */
void encode_ise(
	quant_method quant_level,
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	promise(character_count > 0);

	unsigned int bits = btq_counts[quant_level].bits;
	unsigned int trits = btq_counts[quant_level].trits;
	unsigned int quints = btq_counts[quant_level].quints;
	unsigned int mask = (1 << bits) - 1;

	// Write out trits and bits
	if (trits)
	{
		unsigned int i = 0;
		unsigned int full_trit_blocks = character_count / 5;

		for (unsigned int j = 0; j < full_trit_blocks; j++)
		{
			unsigned int i4 = input_data[i + 4] >> bits;
			unsigned int i3 = input_data[i + 3] >> bits;
			unsigned int i2 = input_data[i + 2] >> bits;
			unsigned int i1 = input_data[i + 1] >> bits;
			unsigned int i0 = input_data[i + 0] >> bits;

			uint8_t T = integer_of_trits[i4][i3][i2][i1][i0];

			// Trits are max of 6-bits, so can pack a single MX value with following 1-2 T bits
			uint8_t pack;

			// Element 0 + T0 + T1
			pack = (input_data[i++] & mask) | (((T >> 0) & 0x3) << bits);
			write_bits(pack, bits + 2, bit_offset, output_data);
			bit_offset += bits + 2;

			// Element 1 + T2 + T3
			pack = (input_data[i++] & mask) | (((T >> 2) & 0x3) << bits);
			write_bits(pack, bits + 2, bit_offset, output_data);
			bit_offset += bits + 2;

			// Element 2 + T4
			pack = (input_data[i++] & mask) | (((T >> 4) & 0x1) << bits);
			write_bits(pack, bits + 1, bit_offset, output_data);
			bit_offset += bits + 1;

			// Element 3 + T5 + T6
			pack = (input_data[i++] & mask) | (((T >> 5) & 0x3) << bits);
			write_bits(pack, bits + 2, bit_offset, output_data);
			bit_offset += bits + 2;

			// Element 4 + T7
			pack = (input_data[i++] & mask) | (((T >> 7) & 0x1) << bits);
			write_bits(pack, bits + 1, bit_offset, output_data);
			bit_offset += bits + 1;
		}

		// Loop tail for a partial block
		if (i != character_count)
		{
			// A partial block has no i4
			// A nonempty block always has i0
			unsigned int i4 =                            0;
			unsigned int i3 = i + 3 >= character_count ? 0 : input_data[i + 3] >> bits;
			unsigned int i2 = i + 2 >= character_count ? 0 : input_data[i + 2] >> bits;
			unsigned int i1 = i + 1 >= character_count ? 0 : input_data[i + 1] >> bits;
			unsigned int i0 =                                input_data[i + 0] >> bits;

			uint8_t T = integer_of_trits[i4][i3][i2][i1][i0];

			for (unsigned int j = 0; i < character_count; i++, j++)
			{
				// Truncated table as this iteration is always partial
				static const uint8_t tbits[4]  { 2, 2, 1, 2 };
				static const uint8_t tshift[4] { 0, 2, 4, 5 };

				uint8_t pack = (input_data[i] & mask) |
				               (((T >> tshift[j]) & ((1 << tbits[j]) - 1)) << bits);

				write_bits(pack, bits + tbits[j], bit_offset, output_data);
				bit_offset += bits + tbits[j];
			}
		}
	}
	// Write out quints and bits
	else if (quints)
	{
		unsigned int i = 0;
		unsigned int full_quint_blocks = character_count / 3;

		for (unsigned int j = 0; j < full_quint_blocks; j++)
		{
			unsigned int i2 = input_data[i + 2] >> bits;
			unsigned int i1 = input_data[i + 1] >> bits;
			unsigned int i0 = input_data[i + 0] >> bits;

			uint8_t T = integer_of_quints[i2][i1][i0];


			// Quints are max of 5-bits, so can pack a single M value with following 2-3 T bits
			uint8_t pack;

			// Element 0
			pack = (input_data[i++] & mask) | (((T >> 0) & 0x7) << bits);
			write_bits(pack, bits + 3, bit_offset, output_data);
			bit_offset += bits + 3;

			// Element 1
			pack = (input_data[i++] & mask) | (((T >> 3) & 0x3) << bits);
			write_bits(pack, bits + 2, bit_offset, output_data);
			bit_offset += bits + 2;

			// Element 2
			pack = (input_data[i++] & mask) | (((T >> 5) & 0x3) << bits);
			write_bits(pack, bits + 2, bit_offset, output_data);
			bit_offset += bits + 2;
		}

		// Loop tail for a partial block
		if (i != character_count)
		{
			// A partial block has no i2
			// A nonempty block always has i0
			unsigned int i2 =                            0;
			unsigned int i1 = i + 1 >= character_count ? 0 : input_data[i + 1] >> bits;
			unsigned int i0 =                                input_data[i + 0] >> bits;

			uint8_t T = integer_of_quints[i2][i1][i0];

			for (unsigned int j = 0; i < character_count; i++, j++)
			{
				// Truncated table as this iteration is always partial
				static const uint8_t tbits[2]  { 3, 2 };
				static const uint8_t tshift[2] { 0, 3 };

				uint8_t pack = (input_data[i] & mask) |
				               (((T >> tshift[j]) & ((1 << tbits[j]) - 1)) << bits);

				write_bits(pack, bits + tbits[j], bit_offset, output_data);
				bit_offset += bits + tbits[j];
			}
		}
	}
	// Write the remaining bits
	else
	{
		for (unsigned int i = 0; i < character_count; i++)
		{
			write_bits(input_data[i], bits, bit_offset, output_data);
			bit_offset += bits;
		}
	}
}

// Packed word loads and stores require little-endian byte order
#if !defined(ASTCENC_BIG_ENDIAN)

/**
 * @brief Pack radix digits into byte lanes for a single table load per ISE group.
 *
 * @tparam EntryCount   The number of radix codes.
 * @tparam ValueCount   The number of digits per code.
 *
 * @param[in] values   The decoded radix digits.
 *
 * @return The digit groups packed into little-endian words.
 */
template<size_t EntryCount, size_t ValueCount>
static constexpr std::array<uint64_t, EntryCount> pack_radix_table(
	const uint8_t (&values)[EntryCount][ValueCount]
) {
	std::array<uint64_t, EntryCount> packed {};
	for (size_t i = 0; i < EntryCount; i++)
	{
		for (size_t j = 0; j < ValueCount; j++)
		{
			packed[i] |= static_cast<uint64_t>(values[i][j]) << (8 * j);
		}
	}

	return packed;
}

/** @brief Three quint digits stored in byte lanes. */
static constexpr auto packed_quints = pack_radix_table(quints_of_integer);

/** @brief Five trit digits stored in byte lanes. */
static constexpr auto packed_trits = pack_radix_table(trits_of_integer);

/**
 * @brief Read a little-endian ISE group without accessing bytes outside the group.
 *
 * The group must span at most seven bytes, including the initial bit offset within a byte.
 *
 * @param[in] input_data   The packed sequence.
 * @param[in] bit_offset   The initial bit offset in the packed sequence.
 * @param[in] bit_count    The number of group bits, at most 56.
 *
 * @return The extracted group bits.
 */
static ASTCENC_SIMD_INLINE uint64_t read_ise_group_packed(
	const uint8_t* input_data,
	unsigned int bit_offset,
	unsigned int bit_count
) {
	assert((bit_offset & 7) + bit_count <= 56);

	unsigned int byte_count = ((bit_offset & 7) + bit_count + 7) / 8;
	const uint8_t* source = input_data + (bit_offset >> 3);
	uint64_t packed = 0;
	unsigned int shift = 0;

	// Input is not padded, so limit reads to the bytes containing this group
	// Fixed-size copies encourage inline loads; a variable-size memcpy can emit a library call
	if (byte_count >= 4)
	{
		uint32_t value;
		std::memcpy(&value, source, 4);
		packed = value;
		source += 4;
		byte_count -= 4;
		shift = 32;
	}

	if (byte_count >= 2)
	{
		uint16_t value;
		std::memcpy(&value, source, 2);
		packed |= static_cast<uint64_t>(value) << shift;
		source += 2;
		byte_count -= 2;
		shift += 16;
	}

	if (byte_count != 0)
	{
		packed |= static_cast<uint64_t>(*source) << shift;
	}

	return extract_bit_range(packed, bit_offset & 7, bit_count);
}

/**
 * @brief Decode binary-only ISE values using packed byte-lane expansion.
 *
 * @tparam Bits   The number of bits per value, in the range 1..8.
 *
 * @param[in]  count         The number of values to decode.
 * @param[in]  input_data    The packed sequence.
 * @param[out] output_data   The decoded values.
 * @param[in]  bit_offset    The initial bit offset in the packed sequence.
 */
template<unsigned int Bits>
static ASTCENC_SIMD_INLINE void decode_ise_bits_packed(
	unsigned int count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	// Seven eight-bit values still fit a word when the sequence starts between bytes
	unsigned int group_size = 8;
	if constexpr (Bits == 8)
	{
		group_size = 7;
	}

	// Select Bits low positions per byte so deposited values occupy separate byte lanes
	uint64_t lane_mask = 0x0101010101010101ull * ((1ull << Bits) - 1);

	// Include the initial offset when finding the end of the readable input
	unsigned int input_bytes = (bit_offset + count * Bits + 7) / 8;

	for (unsigned int i = 0; i < count; i += group_size)
	{
		// Limit the last group to the encoded values that remain
		unsigned int value_count = astc::min(group_size, count - i);
		promise(value_count <= group_size);

		unsigned int bit_count = value_count * Bits;
		uint64_t packed;

		if (input_bytes >= 8)
		{
			// Move the final word load backward when needed to stay inside the input
			unsigned int byte_offset = astc::min(bit_offset >> 3, input_bytes - 8);
			std::memcpy(&packed, input_data + byte_offset, 8);
			packed = extract_bit_range(packed, bit_offset - 8 * byte_offset, bit_count);
		}
		else
		{
			// Read only the available bytes when the input is smaller than a word
			packed = read_ise_group_packed(input_data, bit_offset, bit_count);
		}

		// Spread consecutive input fields into the low bits of separate output bytes
		uint64_t decoded = deposit_bits(packed, lane_mask);
		if (value_count == 8)
		{
			std::memcpy(output_data + i, &decoded, 8);
		}
		else
		{
			// Preserve caller bytes beyond the decoded values when storing a partial word
			vint4 tail(static_cast<int>(decoded), static_cast<int>(decoded >> 32), 0, 0);
			store_nbytes(tail, output_data + i, value_count);
		}

		bit_offset += bit_count;
	}
}

/**
 * @brief Decode fixed-size radix groups using packed bit extraction.
 *
 * @tparam Bits        The number of low binary bits per value.
 * @tparam GroupSize   The number of values in a complete radix group.
 * @tparam RadixBits   The number of radix-code bits in a complete group.
 *
 * @param[in]  count         The number of values to decode.
 * @param[in]  input_data    The packed sequence.
 * @param[out] output_data   The decoded values.
 * @param[in]  bit_offset    The initial bit offset in the packed sequence.
 * @param[in]  low_mask      The group bit positions containing low binary bits.
 * @param[in]  radix_table   The decoded radix digits packed into byte lanes.
 */
template<unsigned int Bits, unsigned int GroupSize, unsigned int RadixBits>
static ASTCENC_SIMD_INLINE void decode_ise_radix_packed(
	unsigned int count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset,
	uint64_t low_mask,
	const uint64_t* radix_table
) {
	// Select the shared radix-code bits between the independent low binary fields
	constexpr unsigned int group_bits = GroupSize * Bits + RadixBits;
	constexpr uint64_t group_mask = (1ull << group_bits) - 1;
	uint64_t code_mask = group_mask & ~low_mask;

	// Select the low binary bits in each unpacked output byte lane
	constexpr uint64_t value_mask = (1ull << Bits) - 1;
	constexpr uint64_t byte_mask = 0xFFFFFFFFFFFFFFFFull >> (64 - GroupSize * 8);
	constexpr uint64_t lane_mask = byte_mask / 255 * value_mask;

	for (unsigned int i = 0; i < count; i += GroupSize)
	{
		unsigned int value_count = astc::min(GroupSize, count - i);
		promise(value_count <= GroupSize);

		// Bound partial-group reads so omitted high radix-code bits decode as zero
		unsigned int radix_bit_count = value_count * RadixBits;
		radix_bit_count = (radix_bit_count + GroupSize - 1) / GroupSize;
		unsigned int bit_count = value_count * Bits + radix_bit_count;
		uint64_t packed = read_ise_group_packed(input_data, bit_offset, bit_count);

		// Expand the radix digits into high bits and merge the low binary fields
		uint64_t radix_code = compress_bits(packed, code_mask);
		unsigned int code = static_cast<unsigned int>(radix_code);
		uint64_t decoded = radix_table[code] << Bits;
		uint64_t low_bits = compress_bits(packed, low_mask);
		decoded |= deposit_bits(low_bits, lane_mask);

		// Preserve caller bytes beyond the end of the decoded sequence
		int low_word = static_cast<int>(decoded);
		int high_word = static_cast<int>(decoded >> 32);
		vint4 values(low_word, high_word, 0, 0);
		store_nbytes(values, output_data + i, value_count);

		bit_offset += bit_count;
	}
}

/**
 * @brief Decode an ISE sequence with a fixed quantization level.
 *
 * @tparam QuantLevel   The sequence quantization level.
 *
 * @param[in]  character_count   The number of encoded values.
 * @param[in]  input_data        The packed input sequence.
 * @param[out] output_data       The decoded values.
 * @param[in]  bit_offset        The initial bit offset in the packed sequence.
 */
template<quant_method QuantLevel>
static void decode_ise_quant(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	promise(character_count > 0);

	constexpr unsigned int bits = btq_counts[QuantLevel].bits;
	constexpr bool is_trit_sequence = btq_counts[QuantLevel].trits != 0;
	constexpr bool is_quint_sequence = btq_counts[QuantLevel].quints != 0;

	if constexpr (is_trit_sequence)
	{
		// Select binary fields before radix fragments of 2, 2, 1, 2, and 1 bits
		static constexpr uint64_t low_masks[7] {
			0x0000000000ull, 0x0000000949ull, 0x0000019B33ull, 0x000039DCE7ull,
			0x00079EF3CFull, 0x00F9F7CF9Full, 0x1F9FBF3F3Full
		};

		decode_ise_radix_packed<bits, 5, 8>(character_count, input_data, output_data,
		                                    bit_offset, low_masks[bits], packed_trits.data());
	}
	else if constexpr (is_quint_sequence)
	{
		// Select binary fields before radix fragments of 3, 2, and 2 bits
		static constexpr uint64_t low_masks[6] {
			0x0000000000ull, 0x0000000091ull, 0x0000000663ull,
			0x00000039C7ull, 0x000001E78Full, 0x00000F9F1Full
		};

		decode_ise_radix_packed<bits, 3, 7>(character_count, input_data, output_data,
		                                    bit_offset, low_masks[bits], packed_quints.data());
	}
	else
	{
		decode_ise_bits_packed<bits>(character_count, input_data, output_data, bit_offset);
	}
}

/* See header for documentation. */
void decode_ise(
	quant_method quant_level,
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	static constexpr std::array<decltype(&decode_ise_quant<QUANT_2>), 21> decoders {{
		decode_ise_quant<QUANT_2>,
		decode_ise_quant<QUANT_3>,
		decode_ise_quant<QUANT_4>,
		decode_ise_quant<QUANT_5>,
		decode_ise_quant<QUANT_6>,
		decode_ise_quant<QUANT_8>,
		decode_ise_quant<QUANT_10>,
		decode_ise_quant<QUANT_12>,
		decode_ise_quant<QUANT_16>,
		decode_ise_quant<QUANT_20>,
		decode_ise_quant<QUANT_24>,
		decode_ise_quant<QUANT_32>,
		decode_ise_quant<QUANT_40>,
		decode_ise_quant<QUANT_48>,
		decode_ise_quant<QUANT_64>,
		decode_ise_quant<QUANT_80>,
		decode_ise_quant<QUANT_96>,
		decode_ise_quant<QUANT_128>,
		decode_ise_quant<QUANT_160>,
		decode_ise_quant<QUANT_192>,
		decode_ise_quant<QUANT_256>,
	}};

	decoders[quant_level](character_count, input_data, output_data, bit_offset);
}

#else

/* See header for documentation. */
void decode_ise(
	quant_method quant_level,
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	promise(character_count > 0);

	// Trit/quint unpacking can write beyond character_count.
	// Reserve four padding entries after the maximum 64 results.
	uint8_t results[68];
	uint8_t tq_blocks[22] { 0 };  // Trit-blocks or quint-blocks, must be zeroed

	unsigned int bits = btq_counts[quant_level].bits;
	unsigned int trits = btq_counts[quant_level].trits;
	unsigned int quints = btq_counts[quant_level].quints;

	unsigned int lcounter = 0;
	unsigned int hcounter = 0;

	// Collect bits for each element, as well as bits for any trit-blocks and quint-blocks
	for (unsigned int i = 0; i < character_count; i++)
	{
		results[i] = static_cast<uint8_t>(read_bits(bits, bit_offset, input_data));
		bit_offset += bits;

		if (trits)
		{
			static const uint8_t bits_to_read[5]  { 2, 2, 1, 2, 1 };
			static const uint8_t block_shift[5]   { 0, 2, 4, 5, 7 };
			static const uint8_t next_lcounter[5] { 1, 2, 3, 4, 0 };
			static const uint8_t hcounter_incr[5] { 0, 0, 0, 0, 1 };
			unsigned int tdata = read_bits(bits_to_read[lcounter], bit_offset, input_data);
			bit_offset += bits_to_read[lcounter];
			tq_blocks[hcounter] |= tdata << block_shift[lcounter];
			hcounter += hcounter_incr[lcounter];
			lcounter = next_lcounter[lcounter];
		}

		if (quints)
		{
			static const uint8_t bits_to_read[3]  { 3, 2, 2 };
			static const uint8_t block_shift[3]   { 0, 3, 5 };
			static const uint8_t next_lcounter[3] { 1, 2, 0 };
			static const uint8_t hcounter_incr[3] { 0, 0, 1 };
			unsigned int tdata = read_bits(bits_to_read[lcounter], bit_offset, input_data);
			bit_offset += bits_to_read[lcounter];
			tq_blocks[hcounter] |= tdata << block_shift[lcounter];
			hcounter += hcounter_incr[lcounter];
			lcounter = next_lcounter[lcounter];
		}
	}

	// Unpack trit-blocks or quint-blocks as needed
	if (trits)
	{
		unsigned int trit_blocks = (character_count + 4) / 5;
		promise(trit_blocks > 0);
		for (unsigned int i = 0; i < trit_blocks; i++)
		{
			const uint8_t *tritptr = trits_of_integer[tq_blocks[i]];
			results[5 * i    ] |= tritptr[0] << bits;
			results[5 * i + 1] |= tritptr[1] << bits;
			results[5 * i + 2] |= tritptr[2] << bits;
			results[5 * i + 3] |= tritptr[3] << bits;
			results[5 * i + 4] |= tritptr[4] << bits;
		}
	}

	if (quints)
	{
		unsigned int quint_blocks = (character_count + 2) / 3;
		promise(quint_blocks > 0);
		for (unsigned int i = 0; i < quint_blocks; i++)
		{
			const uint8_t *quintptr = quints_of_integer[tq_blocks[i]];
			results[3 * i    ] |= quintptr[0] << bits;
			results[3 * i + 1] |= quintptr[1] << bits;
			results[3 * i + 2] |= quintptr[2] << bits;
		}
	}

	for (unsigned int i = 0; i < character_count; i++)
	{
		output_data[i] = results[i];
	}
}

#endif
