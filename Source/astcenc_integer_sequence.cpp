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

/**
 * @brief Functions for encoding/decoding Bounded Integer Sequence Encoding.
 */

#include "astcenc_internal.h"

#include <array>

// unpacked quint triplets <low,middle,high> for each packed-quint value
static const uint8_t quints_of_integer[128][3] = {
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

// packed quint-value for every unpacked quint-triplet
// indexed by [high][middle][low]
static const uint8_t integer_of_quints[5][5][5] = {
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

// unpacked trit quintuplets <low,_,_,_,high> for each packed-quint value
static const uint8_t trits_of_integer[256][5] = {
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

// packed trit-value for every unpacked trit-quintuplet
// indexed by [high][][][][low]
static const uint8_t integer_of_trits[3][3][3][3][3] = {
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
struct btq_count {
	/**< The quantization level. */
	uint8_t quant;

	/**< The number of bits. */
	uint8_t bits;

	/**< The number of trits. */
	uint8_t trits;

	/**< The number of quints. */
	uint8_t quints;
};

/**
 * @brief The table of bits, trits, and quints needed for a quant encode.
 */
static const std::array<btq_count, 21> btq_counts = {{
	{   QUANT_2, 1, 0, 0 },
	{   QUANT_3, 0, 1, 0 },
	{   QUANT_4, 2, 0, 0 },
	{   QUANT_5, 0, 0, 1 },
	{   QUANT_6, 1, 1, 0 },
	{   QUANT_8, 3, 0, 0 },
	{  QUANT_10, 1, 0, 1 },
	{  QUANT_12, 2, 1, 0 },
	{  QUANT_16, 4, 0, 0 },
	{  QUANT_20, 2, 0, 1 },
	{  QUANT_24, 3, 1, 0 },
	{  QUANT_32, 5, 0, 0 },
	{  QUANT_40, 3, 0, 1 },
	{  QUANT_48, 4, 1, 0 },
	{  QUANT_64, 6, 0, 0 },
	{  QUANT_80, 4, 0, 1 },
	{  QUANT_96, 5, 1, 0 },
	{ QUANT_128, 7, 0, 0 },
	{ QUANT_160, 5, 0, 1 },
	{ QUANT_192, 6, 1, 0 },
	{ QUANT_256, 8, 0, 0 }
}};

/**
 * @brief The sequence scale, round, and divisors needed to compute sizing.
 *
 * The length of a quantized sequence in bits is:
 *     (scale * <sequence_len> + round) / divisor
 */
struct ise_size {
	/**< The quantization level. */
	uint8_t quant;

	/**< The scaling parameter. */
	uint8_t scale;

	/**< The rounding parameter. */
	uint8_t round;

	/**< The divisor parameter. */
	uint8_t divisor;
};

/**
 * @brief The table of scale, round, and divisors needed for quant sizing.
 */
static const std::array<ise_size, 21> ise_sizes = {{
	{   QUANT_2,  1, 0, 1 },
	{   QUANT_3,  8, 4, 5 },
	{   QUANT_4,  2, 0, 1 },
	{   QUANT_5,  7, 2, 3 },
	{   QUANT_6, 13, 4, 5 },
	{   QUANT_8,  3, 0, 1 },
	{  QUANT_10, 10, 2, 3 },
	{  QUANT_12, 18, 4, 5 },
	{  QUANT_16,  4, 0, 1 },
	{  QUANT_20, 13, 2, 3 },
	{  QUANT_24, 23, 4, 5 },
	{  QUANT_32,  5, 0, 1 },
	{  QUANT_40, 16, 2, 3 },
	{  QUANT_48, 28, 4, 5 },
	{  QUANT_64,  6, 0, 1 },
	{  QUANT_80, 19, 2, 3 },
	{  QUANT_96, 33, 4, 5 },
	{ QUANT_128,  7, 0, 1 },
	{ QUANT_160, 22, 2, 3 },
	{ QUANT_192, 38, 4, 5 },
	{ QUANT_256,  8, 0, 1 }
}};

/* See header for documentation. */
int get_ise_sequence_bitcount(
	int items,
	quant_method quant
) {
	// Cope with out-of bounds values - input might be invalid
	if (static_cast<size_t>(quant) >= ise_sizes.size())
	{
		// Arbitrary large number that's more than an ASTC block can hold
		return 1024;
	}

	auto& entry = ise_sizes[quant];
	return (entry.scale * items + entry.round) / entry.divisor;
}

// routine to write up to 8 bits
static inline void write_bits(
	int value,
	int bitcount,
	int bitoffset,
	uint8_t* ptr
) {
	int mask = (1 << bitcount) - 1;
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

// routine to read up to 8 bits
static inline int read_bits(
	int bitcount,
	int bitoffset,
	const uint8_t* ptr
) {
	int mask = (1 << bitcount) - 1;
	ptr += bitoffset >> 3;
	bitoffset &= 7;
	int value = ptr[0] | (ptr[1] << 8);
	value >>= bitoffset;
	value &= mask;
	return value;
}

void encode_ise(
	int quant_level,
	int elements,
	const uint8_t* input_data,
	uint8_t* output_data,
	int bit_offset
) {
	uint8_t lowparts[64];
	uint8_t highparts[69];		// 64 elements + 5 elements for padding
	uint8_t tq_blocks[22];		// trit-blocks or quint-blocks

	int bits = btq_counts[quant_level].bits;
	int trits = btq_counts[quant_level].trits;
	int quints = btq_counts[quant_level].quints;

	promise(elements > 0);
	for (int i = 0; i < elements; i++)
	{
		lowparts[i] = input_data[i] & ((1 << bits) - 1);
		highparts[i] = input_data[i] >> bits;
	}

	for (int i = elements; i < elements + 5; i++)
	{
		highparts[i] = 0;		// padding before we start constructing trit-blocks or quint-blocks
	}

	// construct trit-blocks or quint-blocks as necessary
	if (trits)
	{
		int trit_blocks = (elements + 4) / 5;
		for (int i = 0; i < trit_blocks; i++)
		{
			tq_blocks[i] = integer_of_trits[highparts[5 * i + 4]][highparts[5 * i + 3]][highparts[5 * i + 2]][highparts[5 * i + 1]][highparts[5 * i]];
		}
	}

	if (quints)
	{
		int quint_blocks = (elements + 2) / 3;
		for (int i = 0; i < quint_blocks; i++)
		{
			tq_blocks[i] = integer_of_quints[highparts[3 * i + 2]][highparts[3 * i + 1]][highparts[3 * i]];
		}
	}

	// then, write out the actual bits.
	int lcounter = 0;
	int hcounter = 0;
	for (int i = 0; i < elements; i++)
	{
		write_bits(lowparts[i], bits, bit_offset, output_data);
		bit_offset += bits;

		if (trits)
		{
			static const int bits_to_write[5] = { 2, 2, 1, 2, 1 };
			static const int block_shift[5]   = { 0, 2, 4, 5, 7 };
			static const int next_lcounter[5] = { 1, 2, 3, 4, 0 };
			static const int hcounter_incr[5] = { 0, 0, 0, 0, 1 };
			write_bits(tq_blocks[hcounter] >> block_shift[lcounter], bits_to_write[lcounter], bit_offset, output_data);
			bit_offset += bits_to_write[lcounter];
			hcounter += hcounter_incr[lcounter];
			lcounter = next_lcounter[lcounter];
		}

		if (quints)
		{
			static const int bits_to_write[3] = { 3, 2, 2 };
			static const int block_shift[3]   = { 0, 3, 5 };
			static const int next_lcounter[3] = { 1, 2, 0 };
			static const int hcounter_incr[3] = { 0, 0, 1 };
			write_bits(tq_blocks[hcounter] >> block_shift[lcounter], bits_to_write[lcounter], bit_offset, output_data);
			bit_offset += bits_to_write[lcounter];
			hcounter += hcounter_incr[lcounter];
			lcounter = next_lcounter[lcounter];
		}
	}
}

void decode_ise(
	int quant_level,
	int elements,
	const uint8_t* input_data,
	uint8_t* output_data,
	int bit_offset
) {
	// note: due to how the trit/quint-block unpacking is done in this function,
	// we may write more temporary results than the number of outputs
	// The maximum actual number of results is 64 bit, but we keep 4 additional elements
	// of padding.
	uint8_t results[68];
	uint8_t tq_blocks[22];		// trit-blocks or quint-blocks

	int bits = btq_counts[quant_level].bits;
	int trits = btq_counts[quant_level].trits;
	int quints = btq_counts[quant_level].quints;

	int lcounter = 0;
	int hcounter = 0;

	// trit-blocks or quint-blocks must be zeroed out before we collect them in the loop below.
	for (int i = 0; i < 22; i++)
	{
		tq_blocks[i] = 0;
	}

	// collect bits for each element, as well as bits for any trit-blocks and quint-blocks.
	for (int i = 0; i < elements; i++)
	{
		results[i] = read_bits(bits, bit_offset, input_data);
		bit_offset += bits;

		if (trits)
		{
			static const int bits_to_read[5]  = { 2, 2, 1, 2, 1 };
			static const int block_shift[5]   = { 0, 2, 4, 5, 7 };
			static const int next_lcounter[5] = { 1, 2, 3, 4, 0 };
			static const int hcounter_incr[5] = { 0, 0, 0, 0, 1 };
			int tdata = read_bits(bits_to_read[lcounter], bit_offset, input_data);
			bit_offset += bits_to_read[lcounter];
			tq_blocks[hcounter] |= tdata << block_shift[lcounter];
			hcounter += hcounter_incr[lcounter];
			lcounter = next_lcounter[lcounter];
		}

		if (quints)
		{
			static const int bits_to_read[3]  = { 3, 2, 2 };
			static const int block_shift[3]   = { 0, 3, 5 };
			static const int next_lcounter[3] = { 1, 2, 0 };
			static const int hcounter_incr[3] = { 0, 0, 1 };
			int tdata = read_bits(bits_to_read[lcounter], bit_offset, input_data);
			bit_offset += bits_to_read[lcounter];
			tq_blocks[hcounter] |= tdata << block_shift[lcounter];
			hcounter += hcounter_incr[lcounter];
			lcounter = next_lcounter[lcounter];
		}
	}

	// unpack trit-blocks or quint-blocks as needed
	if (trits)
	{
		int trit_blocks = (elements + 4) / 5;
		for (int i = 0; i < trit_blocks; i++)
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
		int quint_blocks = (elements + 2) / 3;
		for (int i = 0; i < quint_blocks; i++)
		{
			const uint8_t *quintptr = quints_of_integer[tq_blocks[i]];
			results[3 * i    ] |= quintptr[0] << bits;
			results[3 * i + 1] |= quintptr[1] << bits;
			results[3 * i + 2] |= quintptr[2] << bits;
		}
	}

	for (int i = 0; i < elements; i++)
	{
		output_data[i] = results[i];
	}
}
