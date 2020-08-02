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

#include "astcenc_mathlib.h"

float3 cross(float3 p, float3 q)
{
	return float3(p.y * q.z - p.z * q.y,
	              p.z * q.x - p.x * q.z,
	              p.x * q.y - p.y * q.x);
}

/* Public function, see header file for detailed documentation */
float astc::log2(float val)
{
	if32 p;
	p.f = val;
	if (p.s < 0x800000)
		p.s = 0x800000; // negative, 0, denormal get clamped to non-denormal.

	// normalize mantissa to range [0.66, 1.33] and extract an exponent
	// in such a way that 1.0 returns 0.
	p.s -= 0x3f2aaaab;
	int expo = p.s >> 23;
	p.s &= 0x7fffff;
	p.s += 0x3f2aaaab;

	float x = p.f - 1.0f;

	// taylor polynomial that, with horner's-rule style evaluation,
	// gives sufficient precision for our use
	// (relative error of about 1 in 10^6)

	float res = (float)expo
	          + x * ( 1.442695040888963f
	          + x * (-0.721347520444482f
	          + x * ( 0.480898346962988f
	          + x * (-0.360673760222241f
	          + x * ( 0.288539008177793f
	          + x * (-0.240449173481494f
	          + x * ( 0.206099291555566f
	          + x * (-0.180336880111120f
	          + x * ( 0.160299448987663f
	          )))))))));
	return res;
}

/* Public function, see header file for detailed documentation */
float astc::atan2(
	float y,
	float x
) {
	const float PI = (float)M_PI;
	const float PI_2 = PI / 2.0f;

	// Handle the discontinuity at x == 0
	if (x == 0.0f)
	{
		if (y > 0.0f)
		{
			return PI_2;
		}
		else if (y == 0.0f)
		{
			return 0.0f;
		}
		return -PI_2;
	}

	float z = y / x;
	float z2 = z * z;
	if (std::fabs(z) < 1.0f)
	{
		float atan = z / (1.0f + (0.28f * z2));
		if (x < 0.0f)
		{
			if (y < 0.0f)
			{
				return atan - PI;
			}
			else
			{
				return atan + PI;
			}
		}
		return atan;
	}
	else
	{
		float atan = PI_2 - (z / (z2 + 0.28f));
		if (y < 0.0f)
		{
			return atan - PI;
		}
		else
		{
			return atan;
		}
	}
}

/**
 * @brief 64-bit rotate left.
 *
 * @param val   The value to rotate.
 * @param count The rotation, in bits.
 */
static inline uint64_t rotl(uint64_t val, int count)
{
	return (val << count) | (val >> (64 - count));
}

/* Public function, see header file for detailed documentation */
void astc::rand_init(uint64_t state[2])
{
	state[0] = 0xfaf9e171cea1ec6bULL;
	state[1] = 0xf1b318cc06af5d71ULL;
}

/* Public function, see header file for detailed documentation */
uint64_t astc::rand(uint64_t state[2])
{
	uint64_t s0 = state[0];
	uint64_t s1 = state[1];
	uint64_t res = s0 + s1;
	s1 ^= s0;
	state[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16);
	state[1] = rotl(s1, 37);
	return res;
}
