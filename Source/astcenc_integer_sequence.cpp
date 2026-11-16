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

#if defined(ASTCENC_RDO)
#include "astcenc_rdo.h"
#endif

#include <array>
#include <cstring>

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
	uint8_t bits:6;

	/** @brief The number of trits. */
	uint8_t trits:1;

	/** @brief The number of quints. */
	uint8_t quints:1;
};

/**
 * @brief The table of bits, trits, and quints needed for a quant encode.
 */
static const std::array<btq_count, 21> btq_counts {{
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
static inline unsigned int read_bits(
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

/**
 * @brief Pack radix digits into byte lanes for a single table load per ISE group.
 *
 * @tparam EntryCount   The number of radix codes.
 * @tparam ValueCount   The number of digits per code.
 * @param  values       The decoded radix digits.
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

#if defined(ASTCENC_RDO)
/* See astcenc_rdo.h for documentation. */
void rdo_decode_ise_radix_weights(
	unsigned int quant_level,
	unsigned int count,
	const uint8_t* input_data,
	uint8_t* output_data
) {
	if (quant_level == QUANT_3)
	{
		unsigned int i = 0;
		for (; i + 5 <= count; i += 5)
		{
			uint64_t weights = packed_trits[input_data[i / 5]] << 5;
			store_le64_prefix(weights, output_data + i, 5);
		}

		unsigned int tail = count - i;
		if (tail != 0)
		{
			unsigned int bits = (tail * 8 + 4) / 5;
			unsigned int code = input_data[i / 5] & ((1u << bits) - 1);
			uint64_t weights = packed_trits[code] << 5;
			for (unsigned int j = 0; j < tail; j++)
			{
				output_data[i + j] = static_cast<uint8_t>(weights >> (8 * j));
			}
		}
	}
	else
	{
		assert(quant_level == QUANT_5);
		unsigned int i = 0;
		unsigned int offset = 0;
		for (; i + 3 <= count; i += 3)
		{
			unsigned int code = read_bits(7, offset, input_data);
			uint64_t weights = packed_quints[code] << 4;
			store_le64_prefix(weights, output_data + i, 3);
			offset += 7;
		}

		unsigned int tail = count - i;
		if (tail != 0)
		{
			unsigned int bits = (tail * 7 + 2) / 3;
			unsigned int code = read_bits(bits, offset, input_data);
			uint64_t weights = packed_quints[code] << 4;
			for (unsigned int j = 0; j < tail; j++)
			{
				output_data[i + j] = static_cast<uint8_t>(weights >> (8 * j));
			}
		}
	}
}
#endif

/**
 * @brief Read a little-endian ISE group without accessing bytes outside the group.
 *
 * @param input_data   The packed sequence.
 * @param bit_offset   The initial bit offset in the packed sequence.
 * @param bit_count    The number of group bits, at most 56.
 *
 * @return The extracted group bits.
 */
static ASTCENC_SIMD_INLINE uint64_t read_ise_group_packed(
	const uint8_t* input_data,
	unsigned int bit_offset,
	unsigned int bit_count
) {
	unsigned int byte_count = ((bit_offset & 7) + bit_count + 7) / 8;
	const uint8_t* source = input_data + (bit_offset >> 3);
	uint64_t packed = 0;
	unsigned int shift = 0;
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
 * @param      bits          The number of bits per value, in the range 1..8.
 * @param      count         The number of values to decode.
 * @param      input_data    The packed sequence.
 * @param[out] output_data   The decoded values.
 * @param      bit_offset    The initial bit offset in the packed sequence.
 */
[[maybe_unused]] static ASTCENC_SIMD_INLINE void decode_ise_bits_packed(
	unsigned int bits,
	unsigned int count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	// Seven eight-bit values still fit a word when the sequence starts between bytes
	unsigned int group_size = 8;
	if (bits == 8)
	{
		group_size = 7;
	}

	uint64_t lane_mask = 0x0101010101010101ull * mask_low_bits(~0ull, bits);
	unsigned int input_bytes = (bit_offset + count * bits + 7) / 8;
	for (unsigned int i = 0; i < count; i += group_size)
	{
		unsigned int value_count = astc::min(group_size, count - i);
		unsigned int bit_count = value_count * bits;
		uint64_t packed;
		if (input_bytes >= 8)
		{
			unsigned int byte_offset = astc::min(bit_offset >> 3, input_bytes - 8);
			std::memcpy(&packed, input_data + byte_offset, 8);
			packed = extract_bit_range(packed, bit_offset - 8 * byte_offset, bit_count);
		}
		else
		{
			packed = read_ise_group_packed(input_data, bit_offset, bit_count);
		}

		uint64_t decoded = deposit_bits(packed, lane_mask);
		if (value_count == 8)
		{
			std::memcpy(output_data + i, &decoded, 8);
		}
		else
		{
			for (unsigned int j = 0; j < value_count; j++)
			{
				output_data[i + j] = static_cast<uint8_t>(decoded >> (8 * j));
			}
		}

		bit_offset += bit_count;
	}
}

/**
 * @brief Decode interleaved five-value trit groups using packed bit extraction.
 *
 * @tparam HasBits   True for mixed binary and radix values.
 *
 * @param      bits          The number of low binary bits per value.
 * @param      count         The number of values to decode.
 * @param      input_data    The packed sequence.
 * @param[out] output_data   The decoded values.
 * @param      bit_offset    The initial bit offset in the packed sequence.
 */
template<bool HasBits>
static ASTCENC_SIMD_INLINE void decode_ise_trits_packed(
	unsigned int bits,
	unsigned int count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	if constexpr (!HasBits)
	{
		bits = 0;
	}

	// Low binary bits and shared radix-code bits occupy disjoint positions in each group
	static constexpr uint64_t low_masks[7] {
		0x0000000000ull, 0x0000000949ull, 0x0000019B33ull, 0x000039DCE7ull,
		0x00079EF3CFull, 0x00F9F7CF9Full, 0x1F9FBF3F3Full
	};
	uint64_t low_mask = low_masks[bits];
	uint64_t lane_mask = 0x0101010101ull * mask_low_bits(~0ull, bits);
	uint64_t code_mask = mask_low_bits(~low_mask, 5 * bits + 8);
	unsigned int full_count = (count / 5) * 5;
	unsigned int bit_count = 5 * bits + 8;
	unsigned int i = 0;
	for (; i < full_count; i += 5)
	{
		uint64_t packed = read_ise_group_packed(input_data, bit_offset, bit_count);
		unsigned int code = static_cast<unsigned int>(compress_bits(packed, code_mask));
		uint64_t decoded = packed_trits[code];
		decoded = (decoded << bits) | deposit_bits(compress_bits(packed, low_mask), lane_mask);
		std::memcpy(output_data + i, &decoded, 5);
		bit_offset += bit_count;
	}

	unsigned int tail_count = count - i;
	if (tail_count != 0)
	{
		bit_count = tail_count * bits + (tail_count * 8 + 4) / 5;
		uint64_t packed = read_ise_group_packed(input_data, bit_offset, bit_count);
		unsigned int code = static_cast<unsigned int>(compress_bits(packed, code_mask));
		uint64_t decoded = packed_trits[code];
		decoded = (decoded << bits) | deposit_bits(compress_bits(packed, low_mask), lane_mask);
		for (unsigned int j = 0; j < tail_count; j++)
		{
			output_data[i + j] = static_cast<uint8_t>(decoded >> (8 * j));
		}
	}
}

/**
 * @brief Decode interleaved three-value quint groups using packed bit extraction.
 *
 * @tparam HasBits   True for mixed binary and radix values.
 *
 * @param      bits          The number of low binary bits per value.
 * @param      count         The number of values to decode.
 * @param      input_data    The packed sequence.
 * @param[out] output_data   The decoded values.
 * @param      bit_offset    The initial bit offset in the packed sequence.
 */
template<bool HasBits>
static ASTCENC_SIMD_INLINE void decode_ise_quints_packed(
	unsigned int bits,
	unsigned int count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	if constexpr (!HasBits)
	{
		bits = 0;
	}

	// Low binary bits and shared radix-code bits occupy disjoint positions in each group
	static constexpr uint64_t low_masks[6] {
		0x0000000000ull, 0x0000000091ull, 0x0000000663ull,
		0x00000039C7ull, 0x000001E78Full, 0x00000F9F1Full
	};
	uint64_t low_mask = low_masks[bits];
	uint64_t lane_mask = 0x010101ull * mask_low_bits(~0ull, bits);
	uint64_t code_mask = mask_low_bits(~low_mask, 3 * bits + 7);
	unsigned int full_count = (count / 3) * 3;
	unsigned int bit_count = 3 * bits + 7;
	unsigned int i = 0;
	for (; i < full_count; i += 3)
	{
		uint64_t packed = read_ise_group_packed(input_data, bit_offset, bit_count);
		unsigned int code = static_cast<unsigned int>(compress_bits(packed, code_mask));
		uint64_t decoded = packed_quints[code];
		decoded = (decoded << bits) | deposit_bits(compress_bits(packed, low_mask), lane_mask);
		std::memcpy(output_data + i, &decoded, 3);
		bit_offset += bit_count;
	}

	unsigned int tail_count = count - i;
	if (tail_count != 0)
	{
		bit_count = tail_count * bits + (tail_count * 7 + 2) / 3;
		uint64_t packed = read_ise_group_packed(input_data, bit_offset, bit_count);
		unsigned int code = static_cast<unsigned int>(compress_bits(packed, code_mask));
		uint64_t decoded = packed_quints[code];
		decoded = (decoded << bits) | deposit_bits(compress_bits(packed, low_mask), lane_mask);
		for (unsigned int j = 0; j < tail_count; j++)
		{
			output_data[i + j] = static_cast<uint8_t>(decoded >> (8 * j));
		}
	}
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

/**
 * @brief Decode an ISE sequence with explicit binary and radix parameters.
 *
 * @param      bits                The binary bits per value.
 * @param      is_trit_sequence    True if the sequence has trit groups.
 * @param      is_quint_sequence   True if the sequence has quint groups.
 * @param      character_count     The number of encoded values.
 * @param      input_data          The packed input sequence.
 * @param[out] output_data         The decoded values.
 * @param      bit_offset          The initial input bit offset.
 */
static ASTCENC_SIMD_INLINE void decode_ise_fixed(
	unsigned int bits,
	bool is_trit_sequence,
	bool is_quint_sequence,
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	promise(character_count > 0);

	if constexpr (is_bit_permute_fast)
	{
		if (is_trit_sequence)
		{
			if (bits == 0)
			{
				decode_ise_trits_packed<false>(bits, character_count, input_data, output_data, bit_offset);
			}
			else
			{
				decode_ise_trits_packed<true>(bits, character_count, input_data, output_data, bit_offset);
			}

			return;
		}

		if (is_quint_sequence)
		{
			if (bits == 0)
			{
				decode_ise_quints_packed<false>(bits, character_count, input_data, output_data, bit_offset);
			}
			else
			{
				decode_ise_quints_packed<true>(bits, character_count, input_data, output_data, bit_offset);
			}

			return;
		}

		decode_ise_bits_packed(bits, character_count, input_data, output_data, bit_offset);
	}
	else
	{
		// Trit/quint unpacking can write beyond character_count. Reserve four padding entries after
		// the maximum 64 results.
		uint8_t results[68];
		uint8_t tq_blocks[22] { 0 };  // Trit-blocks or quint-blocks, must be zeroed

		unsigned int lcounter = 0;
		unsigned int hcounter = 0;

		// Collect bits for each element, as well as bits for any trit-blocks and quint-blocks
		for (unsigned int i = 0; i < character_count; i++)
		{
			results[i] = static_cast<uint8_t>(read_bits(bits, bit_offset, input_data));
			bit_offset += bits;

			if (is_trit_sequence)
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

			if (is_quint_sequence)
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
		if (is_trit_sequence)
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

		if (is_quint_sequence)
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
}

/* See header for documentation. */
void decode_ise(
	quant_method quant_level,
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	const auto& counts = btq_counts[quant_level];
	decode_ise_fixed(counts.bits, counts.trits != 0, counts.quints != 0,
	                 character_count, input_data, output_data, bit_offset);
}

#if defined(ASTCENC_RDO)

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_2(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(1, false, false, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_3(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(0, true, false, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_4(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(2, false, false, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_5(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(0, false, true, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_6(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(1, true, false, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_8(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(3, false, false, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_10(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(1, false, true, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_12(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(2, true, false, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_16(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(4, false, false, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_20(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(2, false, true, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_24(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(3, true, false, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_32(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(5, false, false, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_40(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(3, false, true, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_48(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(4, true, false, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_64(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(6, false, false, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_80(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(4, false, true, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_96(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(5, true, false, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_128(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(7, false, false, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_160(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(5, false, true, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_192(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(6, true, false, character_count,
	                 input_data, output_data, bit_offset);
}

/** @copydoc rdo_ise_decode_func */
static void rdo_decode_ise_quant_256(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset
) {
	decode_ise_fixed(8, false, false, character_count,
	                 input_data, output_data, bit_offset);
}

/* See header for documentation. */
rdo_ise_decode_func* rdo_get_ise_decoder(unsigned int quant_level)
{
	static rdo_ise_decode_func* const decoders[] {
		rdo_decode_ise_quant_2,
		rdo_decode_ise_quant_3,
		rdo_decode_ise_quant_4,
		rdo_decode_ise_quant_5,
		rdo_decode_ise_quant_6,
		rdo_decode_ise_quant_8,
		rdo_decode_ise_quant_10,
		rdo_decode_ise_quant_12,
		rdo_decode_ise_quant_16,
		rdo_decode_ise_quant_20,
		rdo_decode_ise_quant_24,
		rdo_decode_ise_quant_32,
		rdo_decode_ise_quant_40,
		rdo_decode_ise_quant_48,
		rdo_decode_ise_quant_64,
		rdo_decode_ise_quant_80,
		rdo_decode_ise_quant_96,
		rdo_decode_ise_quant_128,
		rdo_decode_ise_quant_160,
		rdo_decode_ise_quant_192,
		rdo_decode_ise_quant_256,
	};
	assert(quant_level < 21);
	return decoders[quant_level];
}
#endif
