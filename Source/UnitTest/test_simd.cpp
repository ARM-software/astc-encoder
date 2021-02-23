// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2020-2021 Arm Limited
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
 * @brief Unit tests for the vectorized SIMD functionality.
 *
 * This test suite is a partial implementation, focussing on 4-wide vectors.
 * We're adding things as we touch related parts of the code, but there is some
 * technical debt to catch up on to get full coverage.
 */

#include <limits>

#include "gtest/gtest.h"

#include "../astcenc_internal.h"
#include "../astcenc_vecmathlib.h"

namespace astcenc
{

// Misc utility tests - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static int round_down(int x)
{
	int remainder = x % ASTCENC_SIMD_WIDTH;
	return x - remainder;
}

static int round_up(int x)
{
	int remainder = x % ASTCENC_SIMD_WIDTH;
	if (!remainder)
	{
		return x;
	}

	return x - remainder + ASTCENC_SIMD_WIDTH;
}

/** @brief Test VLA loop limit round down. */
TEST(misc, RoundDownVLA)
{
	// Static ones which are valid for all VLA widths
	EXPECT_EQ(round_down_to_simd_multiple_vla(0),  0);
	EXPECT_EQ(round_down_to_simd_multiple_vla(8),  8);
	EXPECT_EQ(round_down_to_simd_multiple_vla(16), 16);

	// Variable ones which depend on VLA width
	EXPECT_EQ(round_down_to_simd_multiple_vla(3),   round_down(3));
	EXPECT_EQ(round_down_to_simd_multiple_vla(5),   round_down(5));
	EXPECT_EQ(round_down_to_simd_multiple_vla(7),   round_down(7));
	EXPECT_EQ(round_down_to_simd_multiple_vla(231), round_down(231));
}

/** @brief Test VLA loop limit round up. */
TEST(misc, RoundUpVLA)
{
	// Static ones which are valid for all VLA widths
	EXPECT_EQ(round_up_to_simd_multiple_vla(0),  0);
	EXPECT_EQ(round_up_to_simd_multiple_vla(8),  8);
	EXPECT_EQ(round_up_to_simd_multiple_vla(16), 16);

	// Variable ones which depend on VLA width
	EXPECT_EQ(round_up_to_simd_multiple_vla(3),   round_up(3));
	EXPECT_EQ(round_up_to_simd_multiple_vla(5),   round_up(5));
	EXPECT_EQ(round_up_to_simd_multiple_vla(7),   round_up(7));
	EXPECT_EQ(round_up_to_simd_multiple_vla(231), round_up(231));
}

#if ASTCENC_SIMD_WIDTH == 1

// VLA (1-wide) tests - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

/** @brief Test VLA change_sign. */
TEST(vfloat, ChangeSign)
{
	vfloat a0(-1.0f);
	vfloat b0(-1.0f);
	vfloat r0 = change_sign(a0, b0);
	EXPECT_EQ(r0.lane<0>(), 1.0f);

	vfloat a1( 1.0f);
	vfloat b1(-1.0f);
	vfloat r1 = change_sign(a1, b1);
	EXPECT_EQ(r1.lane<0>(), -1.0f);

	vfloat a2(-3.12f);
	vfloat b2( 3.12f);
	vfloat r2 = change_sign(a2, b2);
	EXPECT_EQ(r2.lane<0>(), -3.12f);

	vfloat a3( 3.12f);
	vfloat b3( 3.12f);
	vfloat r3 = change_sign(a3, b3);
	EXPECT_EQ(r3.lane<0>(), 3.12f);
}

/** @brief Test VLA atan. */
TEST(vfloat, Atan)
{
	vfloat a0(-0.15f);
	vfloat r0 = atan(a0);
	EXPECT_NEAR(r0.lane<0>(), -0.149061f, 0.005f);

	vfloat a1(0.0f);
	vfloat r1 = atan(a1);
	EXPECT_NEAR(r1.lane<0>(),  0.000000f, 0.005f);

	vfloat a2(0.9f);
	vfloat r2 = atan(a2);
	EXPECT_NEAR(r2.lane<0>(),  0.733616f, 0.005f);

	vfloat a3(2.1f);
	vfloat r3 = atan(a3);
	EXPECT_NEAR(r3.lane<0>(),  1.123040f, 0.005f);
}

/** @brief Test VLA atan2. */
TEST(vfloat, Atan2)
{
	vfloat a0(-0.15f);
	vfloat b0( 1.15f);
	vfloat r0 = atan2(a0, b0);
	EXPECT_NEAR(r0.lane<0>(), -0.129816f, 0.005f);

	vfloat a1( 0.0f);
	vfloat b1(-3.0f);
	vfloat r1 = atan2(a1, b1);
	EXPECT_NEAR(r1.lane<0>(),  3.141592f, 0.005f);

	vfloat a2( 0.9f);
	vfloat b2(-0.9f);
	vfloat r2 = atan2(a2, b2);
	EXPECT_NEAR(r2.lane<0>(),  2.360342f, 0.005f);

	vfloat a3( 2.1f);
	vfloat b3( 1.1f);
	vfloat r3 = atan2(a3, b3);
	EXPECT_NEAR(r3.lane<0>(),  1.084357f, 0.005f);
}

#elif ASTCENC_SIMD_WIDTH == 4

// VLA (4-wide) tests - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

/** @brief Test VLA change_sign. */
TEST(vfloat, ChangeSign)
{
	vfloat a(-1.0f,  1.0f, -3.12f, 3.12f);
	vfloat b(-1.0f, -1.0f,  3.12f, 3.12f);
	vfloat r = change_sign(a, b);
	EXPECT_EQ(r.lane<0>(),  1.0f);
	EXPECT_EQ(r.lane<1>(), -1.0f);
	EXPECT_EQ(r.lane<2>(), -3.12f);
	EXPECT_EQ(r.lane<3>(),  3.12f);
}

/** @brief Test VLA atan. */
TEST(vfloat, Atan)
{
	vfloat a(-0.15f, 0.0f, 0.9f, 2.1f);
	vfloat r = atan(a);
	EXPECT_NEAR(r.lane<0>(), -0.149061f, 0.005f);
	EXPECT_NEAR(r.lane<1>(),  0.000000f, 0.005f);
	EXPECT_NEAR(r.lane<2>(),  0.733616f, 0.005f);
	EXPECT_NEAR(r.lane<3>(),  1.123040f, 0.005f);
}

/** @brief Test VLA atan2. */
TEST(vfloat, Atan2)
{
	vfloat a(-0.15f, 0.0f, 0.9f, 2.1f);
	vfloat b(1.15f, -3.0f, -0.9f, 1.1f);
	vfloat r = atan2(a, b);
	EXPECT_NEAR(r.lane<0>(), -0.129816f, 0.005f);
	EXPECT_NEAR(r.lane<1>(),  3.141592f, 0.005f);
	EXPECT_NEAR(r.lane<2>(),  2.360342f, 0.005f);
	EXPECT_NEAR(r.lane<3>(),  1.084357f, 0.005f);
}

#elif ASTCENC_SIMD_WIDTH == 8

// VLA (8-wide) tests - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

/** @brief Test VLA change_sign. */
TEST(vfloat, ChangeSign)
{
	vfloat a(-1.0f,  1.0f, -3.12f, 3.12f, -1.0f,  1.0f, -3.12f, 3.12f);
	vfloat b(-1.0f, -1.0f,  3.12f, 3.12f, -1.0f, -1.0f,  3.12f, 3.12f);
	vfloat r = change_sign(a, b);
	EXPECT_EQ(r.lane<0>(),  1.0f);
	EXPECT_EQ(r.lane<1>(), -1.0f);
	EXPECT_EQ(r.lane<2>(), -3.12f);
	EXPECT_EQ(r.lane<3>(),  3.12f);
	EXPECT_EQ(r.lane<4>(),  1.0f);
	EXPECT_EQ(r.lane<5>(), -1.0f);
	EXPECT_EQ(r.lane<6>(), -3.12f);
	EXPECT_EQ(r.lane<7>(),  3.12f);
}

/** @brief Test VLA atan. */
TEST(vfloat, Atan)
{
	vfloat a(-0.15f, 0.0f, 0.9f, 2.1f, -0.15f, 0.0f, 0.9f, 2.1f);
	vfloat r = atan(a);
	EXPECT_NEAR(r.lane<0>(), -0.149061f, 0.005f);
	EXPECT_NEAR(r.lane<1>(),  0.000000f, 0.005f);
	EXPECT_NEAR(r.lane<2>(),  0.733616f, 0.005f);
	EXPECT_NEAR(r.lane<3>(),  1.123040f, 0.005f);
	EXPECT_NEAR(r.lane<4>(), -0.149061f, 0.005f);
	EXPECT_NEAR(r.lane<5>(),  0.000000f, 0.005f);
	EXPECT_NEAR(r.lane<6>(),  0.733616f, 0.005f);
	EXPECT_NEAR(r.lane<7>(),  1.123040f, 0.005f);
}

/** @brief Test VLA atan2. */
TEST(vfloat, Atan2)
{
	vfloat a(-0.15f, 0.0f, 0.9f, 2.1f, -0.15f, 0.0f, 0.9f, 2.1f);
	vfloat b(1.15f, -3.0f, -0.9f, 1.1f, 1.15f, -3.0f, -0.9f, 1.1f);
	vfloat r = atan2(a, b);
	EXPECT_NEAR(r.lane<0>(), -0.129816f, 0.005f);
	EXPECT_NEAR(r.lane<1>(),  3.141592f, 0.005f);
	EXPECT_NEAR(r.lane<2>(),  2.360342f, 0.005f);
	EXPECT_NEAR(r.lane<3>(),  1.084357f, 0.005f);
	EXPECT_NEAR(r.lane<4>(), -0.129816f, 0.005f);
	EXPECT_NEAR(r.lane<5>(),  3.141592f, 0.005f);
	EXPECT_NEAR(r.lane<6>(),  2.360342f, 0.005f);
	EXPECT_NEAR(r.lane<7>(),  1.084357f, 0.005f);
}

#endif

static const float qnan = std::numeric_limits<float>::quiet_NaN();

alignas(32) static const float f32_data[9] {
	0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f
};

alignas(32) static const int s32_data[9] {
	0, 1, 2, 3, 4, 5 , 6, 7, 8
};

alignas(32) static const uint8_t u8_data[9] {
	0, 1, 2, 3, 4, 5 , 6, 7, 8
};

// VFLOAT4 tests - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

/** @brief Test unaligned vfloat4 data load. */
TEST(vfloat4, UnalignedLoad)
{
	vfloat4 a(&(f32_data[1]));
	EXPECT_EQ(a.lane<0>(), 1.0f);
	EXPECT_EQ(a.lane<1>(), 2.0f);
	EXPECT_EQ(a.lane<2>(), 3.0f);
	EXPECT_EQ(a.lane<3>(), 4.0f);
}

/** @brief Test scalar duplicated vfloat4 load. */
TEST(vfloat4, ScalarDupLoad)
{
	vfloat4 a(1.1f);
	EXPECT_EQ(a.lane<0>(), 1.1f);
	EXPECT_EQ(a.lane<1>(), 1.1f);
	EXPECT_EQ(a.lane<2>(), 1.1f);
	EXPECT_EQ(a.lane<3>(), 1.1f);
}

/** @brief Test scalar vfloat4 load. */
TEST(vfloat4, ScalarLoad)
{
	vfloat4 a(1.1f, 2.2f, 3.3f, 4.4f);
	EXPECT_EQ(a.lane<0>(), 1.1f);
	EXPECT_EQ(a.lane<1>(), 2.2f);
	EXPECT_EQ(a.lane<2>(), 3.3f);
	EXPECT_EQ(a.lane<3>(), 4.4f);
}

/** @brief Test copy vfloat4 load. */
TEST(vfloat4, CopyLoad)
{
	vfloat4 s(1.1f, 2.2f, 3.3f, 4.4f);
	vfloat4 a(s.m);
	EXPECT_EQ(a.lane<0>(), 1.1f);
	EXPECT_EQ(a.lane<1>(), 2.2f);
	EXPECT_EQ(a.lane<2>(), 3.3f);
	EXPECT_EQ(a.lane<3>(), 4.4f);
}

/** @brief Test vfloat4 scalar lane set. */
TEST(vfloat4, SetLane)
{
	vfloat4 a(0.0f);

	a.set_lane<0>(1.0f);
	EXPECT_EQ(a.lane<0>(), 1.0f);
	EXPECT_EQ(a.lane<1>(), 0.0f);
	EXPECT_EQ(a.lane<2>(), 0.0f);
	EXPECT_EQ(a.lane<3>(), 0.0f);

	a.set_lane<1>(2.0f);
	EXPECT_EQ(a.lane<0>(), 1.0f);
	EXPECT_EQ(a.lane<1>(), 2.0f);
	EXPECT_EQ(a.lane<2>(), 0.0f);
	EXPECT_EQ(a.lane<3>(), 0.0f);

	a.set_lane<2>(3.0f);
	EXPECT_EQ(a.lane<0>(), 1.0f);
	EXPECT_EQ(a.lane<1>(), 2.0f);
	EXPECT_EQ(a.lane<2>(), 3.0f);
	EXPECT_EQ(a.lane<3>(), 0.0f);

	a.set_lane<3>(4.0f);
	EXPECT_EQ(a.lane<0>(), 1.0f);
	EXPECT_EQ(a.lane<1>(), 2.0f);
	EXPECT_EQ(a.lane<2>(), 3.0f);
	EXPECT_EQ(a.lane<3>(), 4.0f);
}

/** @brief Test vfloat4 zero. */
TEST(vfloat4, Zero)
{
	vfloat4 a = vfloat4::zero();
	EXPECT_EQ(a.lane<0>(), 0.0f);
	EXPECT_EQ(a.lane<1>(), 0.0f);
	EXPECT_EQ(a.lane<2>(), 0.0f);
	EXPECT_EQ(a.lane<3>(), 0.0f);
}

/** @brief Test vfloat4 load1. */
TEST(vfloat4, Load1)
{
	float s = 3.14f;
	vfloat4 a = vfloat4::load1(&s);
	EXPECT_EQ(a.lane<0>(), 3.14f);
	EXPECT_EQ(a.lane<1>(), 3.14f);
	EXPECT_EQ(a.lane<2>(), 3.14f);
	EXPECT_EQ(a.lane<3>(), 3.14f);
}

/** @brief Test vfloat4 loada. */
TEST(vfloat4, Loada)
{
	vfloat4 a(&(f32_data[0]));
	EXPECT_EQ(a.lane<0>(), 0.0f);
	EXPECT_EQ(a.lane<1>(), 1.0f);
	EXPECT_EQ(a.lane<2>(), 2.0f);
	EXPECT_EQ(a.lane<3>(), 3.0f);
}

/** @brief Test vfloat4 lane_id. */
TEST(vfloat4, LaneID)
{
	vfloat4 a = vfloat4::lane_id();
	EXPECT_EQ(a.lane<0>(), 0.0f);
	EXPECT_EQ(a.lane<1>(), 1.0f);
	EXPECT_EQ(a.lane<2>(), 2.0f);
	EXPECT_EQ(a.lane<3>(), 3.0f);
}

/** @brief Test vfloat4 swz to float4. */
TEST(vfloat4, swz4)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 r = a.swz<0, 3, 2, 1>();
	EXPECT_EQ(r.lane<0>(), 1.0f);
	EXPECT_EQ(r.lane<1>(), 4.0f);
	EXPECT_EQ(r.lane<2>(), 3.0f);
	EXPECT_EQ(r.lane<3>(), 2.0f);

	r = a.swz<3, 1, 1, 0>();
	EXPECT_EQ(r.lane<0>(), 4.0f);
	EXPECT_EQ(r.lane<1>(), 2.0f);
	EXPECT_EQ(r.lane<2>(), 2.0f);
	EXPECT_EQ(r.lane<3>(), 1.0f);
}

/** @brief Test vfloat4 swz to float3. */
TEST(vfloat4, swz3)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 r = a.swz<0, 3, 2>();
	EXPECT_EQ(r.lane<0>(), 1.0f);
	EXPECT_EQ(r.lane<1>(), 4.0f);
	EXPECT_EQ(r.lane<2>(), 3.0f);
	EXPECT_EQ(r.lane<3>(), 0.0f);

	r = a.swz<3, 1, 1>();
	EXPECT_EQ(r.lane<0>(), 4.0f);
	EXPECT_EQ(r.lane<1>(), 2.0f);
	EXPECT_EQ(r.lane<2>(), 2.0f);
	EXPECT_EQ(r.lane<3>(), 0.0f);
}

/** @brief Test vfloat4 swz to float2. */
TEST(vfloat4, swz2)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	float2 r = a.swz<0, 3>();
	EXPECT_EQ(r.r, 1.0f);
	EXPECT_EQ(r.g, 4.0f);

	r = a.swz<2, 1>();
	EXPECT_EQ(r.r, 3.0f);
	EXPECT_EQ(r.g, 2.0f);
}

/** @brief Test vfloat4 add. */
TEST(vfloat4, vadd)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b(0.1f, 0.2f, 0.3f, 0.4f);
	a = a + b;
	EXPECT_EQ(a.lane<0>(), 1.0f + 0.1f);
	EXPECT_EQ(a.lane<1>(), 2.0f + 0.2f);
	EXPECT_EQ(a.lane<2>(), 3.0f + 0.3f);
	EXPECT_EQ(a.lane<3>(), 4.0f + 0.4f);
}

/** @brief Test vfloat4 sub. */
TEST(vfloat4, vsub)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b(0.1f, 0.2f, 0.3f, 0.4f);
	a = a - b;
	EXPECT_EQ(a.lane<0>(), 1.0f - 0.1f);
	EXPECT_EQ(a.lane<1>(), 2.0f - 0.2f);
	EXPECT_EQ(a.lane<2>(), 3.0f - 0.3f);
	EXPECT_EQ(a.lane<3>(), 4.0f - 0.4f);
}

/** @brief Test vfloat4 mul. */
TEST(vfloat4, vmul)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b(0.1f, 0.2f, 0.3f, 0.4f);
	a = a * b;
	EXPECT_EQ(a.lane<0>(), 1.0f * 0.1f);
	EXPECT_EQ(a.lane<1>(), 2.0f * 0.2f);
	EXPECT_EQ(a.lane<2>(), 3.0f * 0.3f);
	EXPECT_EQ(a.lane<3>(), 4.0f * 0.4f);
}

/** @brief Test vfloat4 mul. */
TEST(vfloat4, vsmul)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	float b = 3.14f;
	a = a * b;
	EXPECT_EQ(a.lane<0>(), 1.0f * 3.14f);
	EXPECT_EQ(a.lane<1>(), 2.0f * 3.14f);
	EXPECT_EQ(a.lane<2>(), 3.0f * 3.14f);
	EXPECT_EQ(a.lane<3>(), 4.0f * 3.14f);
}

/** @brief Test vfloat4 mul. */
TEST(vfloat4, svmul)
{
	float a = 3.14f;
	vfloat4 b(1.0f, 2.0f, 3.0f, 4.0f);
	b = a * b;
	EXPECT_EQ(b.lane<0>(), 3.14f * 1.0f);
	EXPECT_EQ(b.lane<1>(), 3.14f * 2.0f);
	EXPECT_EQ(b.lane<2>(), 3.14f * 3.0f);
	EXPECT_EQ(b.lane<3>(), 3.14f * 4.0f);
}

/** @brief Test vfloat4 div. */
TEST(vfloat4, vdiv)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b(0.1f, 0.2f, 0.3f, 0.4f);
	a = a / b;
	EXPECT_EQ(a.lane<0>(), 1.0f / 0.1f);
	EXPECT_EQ(a.lane<1>(), 2.0f / 0.2f);
	EXPECT_EQ(a.lane<2>(), 3.0f / 0.3f);
	EXPECT_EQ(a.lane<3>(), 4.0f / 0.4f);
}

/** @brief Test vfloat4 div. */
TEST(vfloat4, vsdiv)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	float b = 0.3f;
	a = a / b;
	EXPECT_EQ(a.lane<0>(), 1.0f / 0.3f);
	EXPECT_EQ(a.lane<1>(), 2.0f / 0.3f);
	EXPECT_EQ(a.lane<2>(), 3.0f / 0.3f);
	EXPECT_EQ(a.lane<3>(), 4.0f / 0.3f);
}

/** @brief Test vfloat4 div. */
TEST(vfloat4, svdiv)
{
	float a = 3.0f;
	vfloat4 b(0.1f, 0.2f, 0.3f, 0.4f);
	b = a / b;
	EXPECT_EQ(b.lane<0>(), 3.0f / 0.1f);
	EXPECT_EQ(b.lane<1>(), 3.0f / 0.2f);
	EXPECT_EQ(b.lane<2>(), 3.0f / 0.3f);
	EXPECT_EQ(b.lane<3>(), 3.0f / 0.4f);
}

/** @brief Test vfloat4 ceq. */
TEST(vfloat4, ceq)
{
	vfloat4 a1(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b1(0.1f, 0.2f, 0.3f, 0.4f);
	vmask4 r1 = a1 == b1;
	EXPECT_EQ(0, mask(r1));
	EXPECT_EQ(false, any(r1));
	EXPECT_EQ(false, all(r1));

	vfloat4 a2(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b2(1.0f, 0.2f, 0.3f, 0.4f);
	vmask4 r2 = a2 == b2;
	EXPECT_EQ(0x1, mask(r2));
	EXPECT_EQ(true, any(r2));
	EXPECT_EQ(false, all(r2));

	vfloat4 a3(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b3(1.0f, 0.2f, 3.0f, 0.4f);
	vmask4 r3 = a3 == b3;
	EXPECT_EQ(0x5, mask(r3));
	EXPECT_EQ(true, any(r3));
	EXPECT_EQ(false, all(r3));

	vfloat4 a4(1.0f, 2.0f, 3.0f, 4.0f);
	vmask4 r4 = a4 == a4;
	EXPECT_EQ(0xF, mask(r4));
	EXPECT_EQ(true, any(r4));
	EXPECT_EQ(true, all(r4));
}

/** @brief Test vfloat4 cne. */
TEST(vfloat4, cne)
{
	vfloat4 a1(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b1(0.1f, 0.2f, 0.3f, 0.4f);
	vmask4 r1 = a1 != b1;
	EXPECT_EQ(0xF, mask(r1));
	EXPECT_EQ(true, any(r1));
	EXPECT_EQ(true, all(r1));

	vfloat4 a2(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b2(1.0f, 0.2f, 0.3f, 0.4f);
	vmask4 r2 = a2 != b2;
	EXPECT_EQ(0xE, mask(r2));
	EXPECT_EQ(true, any(r2));
	EXPECT_EQ(false, all(r2));

	vfloat4 a3(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b3(1.0f, 0.2f, 3.0f, 0.4f);
	vmask4 r3 = a3 != b3;
	EXPECT_EQ(0xA, mask(r3));
	EXPECT_EQ(true, any(r3));
	EXPECT_EQ(false, all(r3));

	vfloat4 a4(1.0f, 2.0f, 3.0f, 4.0f);
	vmask4 r4 = a4 != a4;
	EXPECT_EQ(0, mask(r4));
	EXPECT_EQ(false, any(r4));
	EXPECT_EQ(false, all(r4));
}

/** @brief Test vfloat4 clt. */
TEST(vfloat4, clt)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b(0.9f, 2.1f, 3.0f, 4.1f);
	vmask4 r = a < b;
	EXPECT_EQ(0xA, mask(r));
}

/** @brief Test vfloat4 cle. */
TEST(vfloat4, cle)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b(0.9f, 2.1f, 3.0f, 4.1f);
	vmask4 r = a <= b;
	EXPECT_EQ(0xE, mask(r));
}

/** @brief Test vfloat4 cgt. */
TEST(vfloat4, cgt)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b(0.9f, 2.1f, 3.0f, 4.1f);
	vmask4 r = a > b;
	EXPECT_EQ(0x1, mask(r));
}

/** @brief Test vfloat4 cge. */
TEST(vfloat4, cge)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b(0.9f, 2.1f, 3.0f, 4.1f);
	vmask4 r = a >= b;
	EXPECT_EQ(0x5, mask(r));
}

/** @brief Test vfloat4 min. */
TEST(vfloat4, min)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b(0.9f, 2.1f, 3.0f, 4.1f);
	vfloat4 r = min(a, b);
	EXPECT_EQ(r.lane<0>(), 0.9f);
	EXPECT_EQ(r.lane<1>(), 2.0f);
	EXPECT_EQ(r.lane<2>(), 3.0f);
	EXPECT_EQ(r.lane<3>(), 4.0f);

	float c = 0.3f;
	r = min(a, c);
	EXPECT_EQ(r.lane<0>(), 0.3f);
	EXPECT_EQ(r.lane<1>(), 0.3f);
	EXPECT_EQ(r.lane<2>(), 0.3f);
	EXPECT_EQ(r.lane<3>(), 0.3f);

	float d = 1.5f;
	r = min(a, d);
	EXPECT_EQ(r.lane<0>(), 1.0f);
	EXPECT_EQ(r.lane<1>(), 1.5f);
	EXPECT_EQ(r.lane<2>(), 1.5f);
	EXPECT_EQ(r.lane<3>(), 1.5f);
}

/** @brief Test vfloat4 max. */
TEST(vfloat4, max)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b(0.9f, 2.1f, 3.0f, 4.1f);
	vfloat4 r = max(a, b);
	EXPECT_EQ(r.lane<0>(), 1.0f);
	EXPECT_EQ(r.lane<1>(), 2.1f);
	EXPECT_EQ(r.lane<2>(), 3.0f);
	EXPECT_EQ(r.lane<3>(), 4.1f);

	float c = 4.3f;
	r = max(a, c);
	EXPECT_EQ(r.lane<0>(), 4.3f);
	EXPECT_EQ(r.lane<1>(), 4.3f);
	EXPECT_EQ(r.lane<2>(), 4.3f);
	EXPECT_EQ(r.lane<3>(), 4.3f);

	float d = 1.5f;
	r = max(a, d);
	EXPECT_EQ(r.lane<0>(), 1.5f);
	EXPECT_EQ(r.lane<1>(), 2.0f);
	EXPECT_EQ(r.lane<2>(), 3.0f);
	EXPECT_EQ(r.lane<3>(), 4.0f);
}

/** @brief Test vfloat4 clamp. */
TEST(vfloat4, clamp)
{
	vfloat4 a1(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 r1 = clamp(2.1f, 3.0f, a1);
	EXPECT_EQ(r1.lane<0>(), 2.1f);
	EXPECT_EQ(r1.lane<1>(), 2.1f);
	EXPECT_EQ(r1.lane<2>(), 3.0f);
	EXPECT_EQ(r1.lane<3>(), 3.0f);

	vfloat4 a2(1.0f, 2.0f, qnan, 4.0f);
	vfloat4 r2 = clamp(2.1f, 3.0f, a2);
	EXPECT_EQ(r2.lane<0>(), 2.1f);
	EXPECT_EQ(r2.lane<1>(), 2.1f);
	EXPECT_EQ(r2.lane<2>(), 2.1f);
	EXPECT_EQ(r2.lane<3>(), 3.0f);
}

/** @brief Test vfloat4 clampz. */
TEST(vfloat4, clampz)
{
	vfloat4 a1(-1.0f, 0.0f, 0.1f, 4.0f);
	vfloat4 r1 = clampz(3.0f, a1);
	EXPECT_EQ(r1.lane<0>(), 0.0f);
	EXPECT_EQ(r1.lane<1>(), 0.0f);
	EXPECT_EQ(r1.lane<2>(), 0.1f);
	EXPECT_EQ(r1.lane<3>(), 3.0f);

	vfloat4 a2(-1.0f, 0.0f, qnan, 4.0f);
	vfloat4 r2 = clampz(3.0f, a2);
	EXPECT_EQ(r2.lane<0>(), 0.0f);
	EXPECT_EQ(r2.lane<1>(), 0.0f);
	EXPECT_EQ(r2.lane<2>(), 0.0f);
	EXPECT_EQ(r2.lane<3>(), 3.0f);
}

/** @brief Test vfloat4 clampz. */
TEST(vfloat4, clampzo)
{
	vfloat4 a1(-1.0f, 0.0f, 0.1f, 4.0f);
	vfloat4 r1 = clampzo(a1);
	EXPECT_EQ(r1.lane<0>(), 0.0f);
	EXPECT_EQ(r1.lane<1>(), 0.0f);
	EXPECT_EQ(r1.lane<2>(), 0.1f);
	EXPECT_EQ(r1.lane<3>(), 1.0f);

	vfloat4 a2(-1.0f, 0.0f, qnan, 4.0f);
	vfloat4 r2 = clampzo(a2);
	EXPECT_EQ(r2.lane<0>(), 0.0f);
	EXPECT_EQ(r2.lane<1>(), 0.0f);
	EXPECT_EQ(r2.lane<2>(), 0.0f);
	EXPECT_EQ(r2.lane<3>(), 1.0f);
}

/** @brief Test vfloat4 abs. */
TEST(vfloat4, abs)
{
	vfloat4 a(-1.0f, 0.0f, 0.1f, 4.0f);
	vfloat4 r = abs(a);
	EXPECT_EQ(r.lane<0>(), 1.0f);
	EXPECT_EQ(r.lane<1>(), 0.0f);
	EXPECT_EQ(r.lane<2>(), 0.1f);
	EXPECT_EQ(r.lane<3>(), 4.0f);
}

/** @brief Test vfloat4 round. */
TEST(vfloat4, round)
{
	vfloat4 a(1.1f, 1.5f, 1.6f, 4.0f);
	vfloat4 r = round(a);
	EXPECT_EQ(r.lane<0>(), 1.0f);
	EXPECT_EQ(r.lane<1>(), 2.0f);
	EXPECT_EQ(r.lane<2>(), 2.0f);
	EXPECT_EQ(r.lane<3>(), 4.0f);
}

/** @brief Test vfloat4 hmin. */
TEST(vfloat4, hmin)
{
	vfloat4 a1(1.1f, 1.5f, 1.6f, 4.0f);
	vfloat4 r1 = hmin(a1);
	EXPECT_EQ(r1.lane<0>(), 1.1f);
	EXPECT_EQ(r1.lane<1>(), 1.1f);
	EXPECT_EQ(r1.lane<2>(), 1.1f);
	EXPECT_EQ(r1.lane<3>(), 1.1f);

	vfloat4 a2(1.1f, 1.5f, 1.6f, 0.2f);
	vfloat4 r2 = hmin(a2);
	EXPECT_EQ(r2.lane<0>(), 0.2f);
	EXPECT_EQ(r2.lane<1>(), 0.2f);
	EXPECT_EQ(r2.lane<2>(), 0.2f);
	EXPECT_EQ(r2.lane<3>(), 0.2f);
}

/** @brief Test vfloat4 hmin_s. */
TEST(vfloat4, hmin_s)
{
	vfloat4 a1(1.1f, 1.5f, 1.6f, 4.0f);
	float r1 = hmin_s(a1);
	EXPECT_EQ(r1, 1.1f);

	vfloat4 a2(1.1f, 1.5f, 1.6f, 0.2f);
	float r2 = hmin_s(a2);
	EXPECT_EQ(r2, 0.2f);
}

/** @brief Test vfloat4 hmin_rgb_s. */
TEST(vfloat4, hmin_rgb_s)
{
	vfloat4 a1(1.1f, 1.5f, 1.6f, 0.2f);
	float r1 = hmin_rgb_s(a1);
	EXPECT_EQ(r1, 1.1f);

	vfloat4 a2(1.5f, 0.9f, 1.6f, 1.2f);
	float r2 = hmin_rgb_s(a2);
	EXPECT_EQ(r2, 0.9f);
}

/** @brief Test vfloat4 hmax. */
TEST(vfloat4, hmax)
{
	vfloat4 a1(1.1f, 1.5f, 1.6f, 4.0f);
	vfloat4 r1 = hmax(a1);
	EXPECT_EQ(r1.lane<0>(), 4.0f);
	EXPECT_EQ(r1.lane<1>(), 4.0f);
	EXPECT_EQ(r1.lane<2>(), 4.0f);
	EXPECT_EQ(r1.lane<3>(), 4.0f);

	vfloat4 a2(1.1f, 1.5f, 1.6f, 0.2f);
	vfloat4 r2 = hmax(a2);
	EXPECT_EQ(r2.lane<0>(), 1.6f);
	EXPECT_EQ(r2.lane<1>(), 1.6f);
	EXPECT_EQ(r2.lane<2>(), 1.6f);
	EXPECT_EQ(r2.lane<3>(), 1.6f);
}

/** @brief Test vfloat4 hmax_s. */
TEST(vfloat4, hmax_s)
{
	vfloat4 a1(1.1f, 1.5f, 1.6f, 4.0f);
	float r1 = hmax_s(a1);
	EXPECT_EQ(r1, 4.0f);

	vfloat4 a2(1.1f, 1.5f, 1.6f, 0.2f);
	float r2 = hmax_s(a2);
	EXPECT_EQ(r2, 1.6f);
}

/** @brief Test vfloat4 hadd_s. */
TEST(vfloat4, hadd_s)
{
	vfloat4 a1(1.1f, 1.5f, 1.6f, 4.0f);
	float sum = 1.1f + 1.5f + 1.6f + 4.0f;
	float r = hadd_s(a1);
	EXPECT_NEAR(r, sum, 0.005f);
}

/** @brief Test vfloat4 hadd_rgb_s. */
TEST(vfloat4, hadd_rgb_s)
{
	vfloat4 a1(1.1f, 1.5f, 1.6f, 4.0f);
	float sum = 1.1f + 1.5f + 1.6f;
	float r = hadd_rgb_s(a1);
	EXPECT_NEAR(r, sum, 0.005f);
}

/** @brief Test vfloat4 sqrt. */
TEST(vfloat4, sqrt)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 r = sqrt(a);
	EXPECT_EQ(r.lane<0>(), std::sqrt(1.0f));
	EXPECT_EQ(r.lane<1>(), std::sqrt(2.0f));
	EXPECT_EQ(r.lane<2>(), std::sqrt(3.0f));
	EXPECT_EQ(r.lane<3>(), std::sqrt(4.0f));
}

/** @brief Test vfloat4 select. */
TEST(vfloat4, select)
{
	vfloat4 m1(1.0f, 1.0f, 1.0f, 1.0f);
	vfloat4 m2(1.0f, 2.0f, 1.0f, 2.0f);
	vmask4 cond = m1 == m2;

	vfloat4 a(1.0f, 3.0f, 3.0f, 1.0f);
	vfloat4 b(4.0f, 2.0f, 2.0f, 4.0f);

	// Select in one direction
	vfloat4 r1 = select(a, b, cond);
	EXPECT_EQ(r1.lane<0>(), 4.0f);
	EXPECT_EQ(r1.lane<1>(), 3.0f);
	EXPECT_EQ(r1.lane<2>(), 2.0f);
	EXPECT_EQ(r1.lane<3>(), 1.0f);

	// Select in the other
	vfloat4 r2 = select(b, a, cond);
	EXPECT_EQ(r2.lane<0>(), 1.0f);
	EXPECT_EQ(r2.lane<1>(), 2.0f);
	EXPECT_EQ(r2.lane<2>(), 3.0f);
	EXPECT_EQ(r2.lane<3>(), 4.0f);
}

/** @brief Test vfloat4 select MSB only. */
TEST(vfloat4, select_msb)
{
	vint4 msb(0x80000000, 0, 0x80000000, 0);
	vmask4 cond(msb.m);

	vfloat4 a(1.0f, 3.0f, 3.0f, 1.0f);
	vfloat4 b(4.0f, 2.0f, 2.0f, 4.0f);

	// Select in one direction
	vfloat4 r1 = select(a, b, cond);
	EXPECT_EQ(r1.lane<0>(), 4.0f);
	EXPECT_EQ(r1.lane<1>(), 3.0f);
	EXPECT_EQ(r1.lane<2>(), 2.0f);
	EXPECT_EQ(r1.lane<3>(), 1.0f);

	// Select in the other
	vfloat4 r2 = select(b, a, cond);
	EXPECT_EQ(r2.lane<0>(), 1.0f);
	EXPECT_EQ(r2.lane<1>(), 2.0f);
	EXPECT_EQ(r2.lane<2>(), 3.0f);
	EXPECT_EQ(r2.lane<3>(), 4.0f);
}

/** @brief Test vfloat4 gatherf. */
TEST(vfloat4, gatherf)
{
	vint4 indices(0, 4, 3, 2);
	vfloat4 r = gatherf(f32_data, indices);
	EXPECT_EQ(r.lane<0>(), 0.0f);
	EXPECT_EQ(r.lane<1>(), 4.0f);
	EXPECT_EQ(r.lane<2>(), 3.0f);
	EXPECT_EQ(r.lane<3>(), 2.0f);
}

/** @brief Test vfloat4 storea. */
TEST(vfloat4, storea)
{
	alignas(16) float out[4];
	vfloat4 a(f32_data);
	storea(a, out);
	EXPECT_EQ(out[0], 0.0f);
	EXPECT_EQ(out[1], 1.0f);
	EXPECT_EQ(out[2], 2.0f);
	EXPECT_EQ(out[3], 3.0f);
}

/** @brief Test vfloat4 store. */
TEST(vfloat4, store)
{
	alignas(16) float out[5];
	vfloat4 a(f32_data);
	store(a, &(out[1]));
	EXPECT_EQ(out[1], 0.0f);
	EXPECT_EQ(out[2], 1.0f);
	EXPECT_EQ(out[3], 2.0f);
	EXPECT_EQ(out[4], 3.0f);
}

/** @brief Test vfloat4 dot. */
TEST(vfloat4, dot)
{
	vfloat4 a(1.0f, 2.0f, 4.0f, 8.0f);
	vfloat4 b(1.0f, 0.5f, 0.25f, 0.125f);
	vfloat4 r = dot(a, b);
	EXPECT_EQ(r.lane<0>(), 4.0f);
	EXPECT_EQ(r.lane<1>(), 4.0f);
	EXPECT_EQ(r.lane<2>(), 4.0f);
	EXPECT_EQ(r.lane<3>(), 4.0f);
}

/** @brief Test vfloat4 dot_s. */
TEST(vfloat4, dot_s)
{
	vfloat4 a(1.0f, 2.0f, 4.0f, 8.0f);
	vfloat4 b(1.0f, 0.5f, 0.25f, 0.125f);
	float r = dot_s(a, b);
	EXPECT_EQ(r, 4.0f);
}

/** @brief Test vfloat4 dot3. */
TEST(vfloat4, dot3)
{
	vfloat4 a(1.0f, 2.0f, 4.0f, 8.0f);
	vfloat4 b(1.0f, 0.5f, 0.25f, 0.125f);
	vfloat4 r = dot3(a, b);
	EXPECT_EQ(r.lane<0>(), 3.0f);
	EXPECT_EQ(r.lane<1>(), 3.0f);
	EXPECT_EQ(r.lane<2>(), 3.0f);
	EXPECT_EQ(r.lane<3>(), 0.0f);
}

/** @brief Test vfloat4 dot3_s. */
TEST(vfloat4, dot3_s)
{
	vfloat4 a(1.0f, 2.0f, 4.0f, 8.0f);
	vfloat4 b(1.0f, 0.5f, 0.25f, 0.125f);
	float r = dot3_s(a, b);
	EXPECT_EQ(r, 3.0f);
}

/** @brief Test vfloat4 normalize. */
TEST(vfloat4, normalize)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 r = normalize(a);
	EXPECT_NEAR(r.lane<0>(), 1.0f / astc::sqrt(30.0f), 0.0005f);
	EXPECT_NEAR(r.lane<1>(), 2.0f / astc::sqrt(30.0f), 0.0005f);
	EXPECT_NEAR(r.lane<2>(), 3.0f / astc::sqrt(30.0f), 0.0005f);
	EXPECT_NEAR(r.lane<3>(), 4.0f / astc::sqrt(30.0f), 0.0005f);
}

/** @brief Test vfloat4 normalize_safe. */
TEST(vfloat4, normalize_safe)
{
	vfloat4 s(-1.0f, -1.0f, -1.0f, -1.0f);

	vfloat4 a1(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 r1 = normalize_safe(a1, s);
	EXPECT_NEAR(r1.lane<0>(), 1.0f / astc::sqrt(30.0f), 0.0005f);
	EXPECT_NEAR(r1.lane<1>(), 2.0f / astc::sqrt(30.0f), 0.0005f);
	EXPECT_NEAR(r1.lane<2>(), 3.0f / astc::sqrt(30.0f), 0.0005f);
	EXPECT_NEAR(r1.lane<3>(), 4.0f / astc::sqrt(30.0f), 0.0005f);

	vfloat4 a2(0.0f, 0.0f, 0.0f, 0.0f);
	vfloat4 r2 = normalize_safe(a2, s);
	EXPECT_EQ(r2.lane<0>(), -1.0f);
	EXPECT_EQ(r2.lane<1>(), -1.0f);
	EXPECT_EQ(r2.lane<2>(), -1.0f);
	EXPECT_EQ(r2.lane<3>(), -1.0f);
}

/** @brief Test vfloat4 float_to_int. */
TEST(vfloat4, float_to_int)
{
	vfloat4 a(1.1f, 1.5f, 1.6f, 4.0f);
	vint4 r = float_to_int(a);
	EXPECT_EQ(r.lane<0>(), 1);
	EXPECT_EQ(r.lane<1>(), 1);
	EXPECT_EQ(r.lane<2>(), 1);
	EXPECT_EQ(r.lane<3>(), 4);
}

/** @brief Test vfloat4 round. */
TEST(vfloat4, float_to_int_rtn)
{
	vfloat4 a(1.1f, 1.5f, 1.6f, 4.0f);
	vint4 r = float_to_int_rtn(a);
	EXPECT_EQ(r.lane<0>(), 1);
	EXPECT_EQ(r.lane<1>(), 2);
	EXPECT_EQ(r.lane<2>(), 2);
	EXPECT_EQ(r.lane<3>(), 4);
}

/** @brief Test vfloat4 round. */
TEST(vfloat4, int_to_float)
{
	vint4 a(1, 2, 3, 4);
	vfloat4 r = int_to_float(a);
	EXPECT_EQ(r.lane<0>(), 1.0f);
	EXPECT_EQ(r.lane<1>(), 2.0f);
	EXPECT_EQ(r.lane<2>(), 3.0f);
	EXPECT_EQ(r.lane<3>(), 4.0f);
}

/** @brief Test vfloat4 float to fp16 conversion. */
TEST(vfloat4, float_to_float16)
{
	vfloat4 a(1.5, 234.5, 345345.0, qnan);
	vint4 r = float_to_float16(a);

	// Normal numbers
	EXPECT_EQ(r.lane<0>(), 0x3E00);
	EXPECT_EQ(r.lane<1>(), 0x5B54);

	// Large numbers convert to infinity
	EXPECT_EQ(r.lane<2>(), 0x7C00);

	// NaN must convert to any valid NaN encoding
	EXPECT_EQ((r.lane<3>() >> 10) & 0x1F, 0x1F); // Exponent must be all 1s
	EXPECT_NE(r.lane<3>() & (0x3FF), 0);         // Mantissa must be non-zero
}

/** @brief Test float to fp16 conversion. */
TEST(sfloat, float_to_float16)
{
	int r = float_to_float16(234.5);
	EXPECT_EQ(r, 0x5B54);
}

/** @brief Test vfloat4 fp16 to float conversion. */
TEST(vfloat4, float16_to_float)
{	vint4 a(0x3E00, 0x5B54, 0x7C00, 0xFFFF);
	vfloat4 r = float16_to_float(a);

	// Normal numbers
	EXPECT_EQ(r.lane<0>(), 1.5);
	EXPECT_EQ(r.lane<1>(), 234.5);

	// Infinities must be preserved
	EXPECT_NE(std::isinf(r.lane<2>()), 0);

	// NaNs must be preserved
	EXPECT_NE(std::isnan(r.lane<3>()), 0);
}

/** @brief Test fp16 to float conversion. */
TEST(sfloat, float16_to_float)
{
	float r = float16_to_float(0x5B54);
	EXPECT_EQ(r, 234.5);
}

// VINT4 tests - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

/** @brief Test unaligned vint4 data load. */
TEST(vint4, UnalignedLoad)
{
	vint4 a(&(s32_data[1]));
	EXPECT_EQ(a.lane<0>(), 1);
	EXPECT_EQ(a.lane<1>(), 2);
	EXPECT_EQ(a.lane<2>(), 3);
	EXPECT_EQ(a.lane<3>(), 4);
}

/** @brief Test unaligned vint4 data load. */
TEST(vint4, UnalignedLoad8)
{
	vint4 a(&(u8_data[1]));
	EXPECT_EQ(a.lane<0>(), 1);
	EXPECT_EQ(a.lane<1>(), 2);
	EXPECT_EQ(a.lane<2>(), 3);
	EXPECT_EQ(a.lane<3>(), 4);
}

/** @brief Test scalar duplicated vint4 load. */
TEST(vint4, ScalarDupLoad)
{
	vint4 a(42);
	EXPECT_EQ(a.lane<0>(), 42);
	EXPECT_EQ(a.lane<1>(), 42);
	EXPECT_EQ(a.lane<2>(), 42);
	EXPECT_EQ(a.lane<3>(), 42);
}

/** @brief Test scalar vint4 load. */
TEST(vint4, ScalarLoad)
{
	vint4 a(11, 22, 33, 44);
	EXPECT_EQ(a.lane<0>(), 11);
	EXPECT_EQ(a.lane<1>(), 22);
	EXPECT_EQ(a.lane<2>(), 33);
	EXPECT_EQ(a.lane<3>(), 44);
}

/** @brief Test copy vint4 load. */
TEST(vint4, CopyLoad)
{
	vint4 s(11, 22, 33, 44);
	vint4 a(s.m);
	EXPECT_EQ(a.lane<0>(), 11);
	EXPECT_EQ(a.lane<1>(), 22);
	EXPECT_EQ(a.lane<2>(), 33);
	EXPECT_EQ(a.lane<3>(), 44);
}

/** @brief Test vint4 scalar lane set. */
TEST(int4, SetLane)
{
	vint4 a(0);

	a.set_lane<0>(1);
	EXPECT_EQ(a.lane<0>(), 1);
	EXPECT_EQ(a.lane<1>(), 0);
	EXPECT_EQ(a.lane<2>(), 0);
	EXPECT_EQ(a.lane<3>(), 0);

	a.set_lane<1>(2);
	EXPECT_EQ(a.lane<0>(), 1);
	EXPECT_EQ(a.lane<1>(), 2);
	EXPECT_EQ(a.lane<2>(), 0);
	EXPECT_EQ(a.lane<3>(), 0);

	a.set_lane<2>(3);
	EXPECT_EQ(a.lane<0>(), 1);
	EXPECT_EQ(a.lane<1>(), 2);
	EXPECT_EQ(a.lane<2>(), 3);
	EXPECT_EQ(a.lane<3>(), 0);

	a.set_lane<3>(4);
	EXPECT_EQ(a.lane<0>(), 1);
	EXPECT_EQ(a.lane<1>(), 2);
	EXPECT_EQ(a.lane<2>(), 3);
	EXPECT_EQ(a.lane<3>(), 4);
}

/** @brief Test vint4 zero. */
TEST(vint4, Zero)
{
	vint4 a = vint4::zero();
	EXPECT_EQ(a.lane<0>(), 0);
	EXPECT_EQ(a.lane<1>(), 0);
	EXPECT_EQ(a.lane<2>(), 0);
	EXPECT_EQ(a.lane<3>(), 0);
}

/** @brief Test vint4 load1. */
TEST(vint4, Load1)
{
	int s = 42;
	vint4 a = vint4::load1(&s);
	EXPECT_EQ(a.lane<0>(), 42);
	EXPECT_EQ(a.lane<1>(), 42);
	EXPECT_EQ(a.lane<2>(), 42);
	EXPECT_EQ(a.lane<3>(), 42);
}

/** @brief Test vint4 loada. */
TEST(vint4, Loada)
{
	vint4 a(&(s32_data[0]));
	EXPECT_EQ(a.lane<0>(), 0);
	EXPECT_EQ(a.lane<1>(), 1);
	EXPECT_EQ(a.lane<2>(), 2);
	EXPECT_EQ(a.lane<3>(), 3);
}

/** @brief Test vint4 lane_id. */
TEST(vint4, LaneID)
{
	vint4 a = vint4::lane_id();
	EXPECT_EQ(a.lane<0>(), 0);
	EXPECT_EQ(a.lane<1>(), 1);
	EXPECT_EQ(a.lane<2>(), 2);
	EXPECT_EQ(a.lane<3>(), 3);
}

/** @brief Test vint4 add. */
TEST(vint4, vadd)
{
	vint4 a(1, 2, 3, 4);
	vint4 b(2, 3, 4, 5);
	a = a + b;
	EXPECT_EQ(a.lane<0>(), 1 + 2);
	EXPECT_EQ(a.lane<1>(), 2 + 3);
	EXPECT_EQ(a.lane<2>(), 3 + 4);
	EXPECT_EQ(a.lane<3>(), 4 + 5);
}

/** @brief Test vint4 add. */
TEST(vint4, vsadd)
{
	vint4 a(1, 2, 3, 4);
	int b = 5;
	a = a + b;
	EXPECT_EQ(a.lane<0>(), 1 + 5);
	EXPECT_EQ(a.lane<1>(), 2 + 5);
	EXPECT_EQ(a.lane<2>(), 3 + 5);
	EXPECT_EQ(a.lane<3>(), 4 + 5);
}

/** @brief Test vint4 sub. */
TEST(vint4, vsub)
{
	vint4 a(1, 2, 4, 4);
	vint4 b(2, 3, 3, 5);
	a = a - b;
	EXPECT_EQ(a.lane<0>(), 1 - 2);
	EXPECT_EQ(a.lane<1>(), 2 - 3);
	EXPECT_EQ(a.lane<2>(), 4 - 3);
	EXPECT_EQ(a.lane<3>(), 4 - 5);
}

/** @brief Test vint4 sub. */
TEST(vint4, vssub)
{
	vint4 a(1, 2, 4, 4);
	int b = 5;
	a = a - b;
	EXPECT_EQ(a.lane<0>(), 1 - 5);
	EXPECT_EQ(a.lane<1>(), 2 - 5);
	EXPECT_EQ(a.lane<2>(), 4 - 5);
	EXPECT_EQ(a.lane<3>(), 4 - 5);
}

/** @brief Test vint4 mul. */
TEST(vint4, vmul)
{
	vint4 a(1, 2, 4, 4);
	vint4 b(2, 3, 3, 5);
	a = a * b;
	EXPECT_EQ(a.lane<0>(), 1 * 2);
	EXPECT_EQ(a.lane<1>(), 2 * 3);
	EXPECT_EQ(a.lane<2>(), 4 * 3);
	EXPECT_EQ(a.lane<3>(), 4 * 5);
}

/** @brief Test vint4 mul. */
TEST(vint4, vsmul)
{
	vint4 a(1, 2, 4, 4);
	a = a * 3;
	EXPECT_EQ(a.lane<0>(), 1 * 3);
	EXPECT_EQ(a.lane<1>(), 2 * 3);
	EXPECT_EQ(a.lane<2>(), 4 * 3);
	EXPECT_EQ(a.lane<3>(), 4 * 3);

	vint4 b(1, 2, -4, 4);
	b = b * -3;
	EXPECT_EQ(b.lane<0>(), 1 * -3);
	EXPECT_EQ(b.lane<1>(), 2 * -3);
	EXPECT_EQ(b.lane<2>(), -4 * -3);
	EXPECT_EQ(b.lane<3>(), 4 * -3);
}

/** @brief Test vint4 bitwise invert. */
TEST(vint4, bit_invert)
{
	vint4 a(-1, 0, 1, 2);
	a = ~a;
	EXPECT_EQ(a.lane<0>(), ~-1);
	EXPECT_EQ(a.lane<1>(), ~0);
	EXPECT_EQ(a.lane<2>(), ~1);
	EXPECT_EQ(a.lane<3>(), ~2);
}

/** @brief Test vint4 bitwise or. */
TEST(vint4, bit_vor)
{
	vint4 a(1, 2, 3, 4);
	vint4 b(2, 3, 4, 5);
	a = a | b;
	EXPECT_EQ(a.lane<0>(), 3);
	EXPECT_EQ(a.lane<1>(), 3);
	EXPECT_EQ(a.lane<2>(), 7);
	EXPECT_EQ(a.lane<3>(), 5);
}

TEST(vint4, bit_vsor)
{
	vint4 a(1, 2, 3, 4);
	int b = 2;
	a = a | b;
	EXPECT_EQ(a.lane<0>(), 3);
	EXPECT_EQ(a.lane<1>(), 2);
	EXPECT_EQ(a.lane<2>(), 3);
	EXPECT_EQ(a.lane<3>(), 6);
}

/** @brief Test vint4 bitwise and. */
TEST(vint4, bit_vand)
{
	vint4 a(1, 2, 3, 4);
	vint4 b(2, 3, 4, 5);
	a = a & b;
	EXPECT_EQ(a.lane<0>(), 0);
	EXPECT_EQ(a.lane<1>(), 2);
	EXPECT_EQ(a.lane<2>(), 0);
	EXPECT_EQ(a.lane<3>(), 4);
}

/** @brief Test vint4 bitwise and. */
TEST(vint4, bit_vsand)
{
	vint4 a(1, 2, 3, 4);
	int b = 2;
	a = a & b;
	EXPECT_EQ(a.lane<0>(), 0);
	EXPECT_EQ(a.lane<1>(), 2);
	EXPECT_EQ(a.lane<2>(), 2);
	EXPECT_EQ(a.lane<3>(), 0);
}

/** @brief Test vint4 bitwise xor. */
TEST(vint4, bit_vxor)
{
	vint4 a(1, 2, 3, 4);
	vint4 b(2, 3, 4, 5);
	a = a ^ b;
	EXPECT_EQ(a.lane<0>(), 3);
	EXPECT_EQ(a.lane<1>(), 1);
	EXPECT_EQ(a.lane<2>(), 7);
	EXPECT_EQ(a.lane<3>(), 1);
}

/** @brief Test vint4 bitwise xor. */
TEST(vint4, bit_vsxor)
{
	vint4 a(1, 2, 3, 4);
	int b = 2;
	a = a ^ b;
	EXPECT_EQ(a.lane<0>(), 3);
	EXPECT_EQ(a.lane<1>(), 0);
	EXPECT_EQ(a.lane<2>(), 1);
	EXPECT_EQ(a.lane<3>(), 6);
}

/** @brief Test vint4 ceq. */
TEST(vint4, ceq)
{
	vint4 a1(1, 2, 3, 4);
	vint4 b1(0, 1, 2, 3);
	vmask4 r1 = a1 == b1;
	EXPECT_EQ(0, mask(r1));
	EXPECT_EQ(false, any(r1));
	EXPECT_EQ(false, all(r1));

	vint4 a2(1, 2, 3, 4);
	vint4 b2(1, 0, 0, 0);
	vmask4 r2 = a2 == b2;
	EXPECT_EQ(0x1, mask(r2));
	EXPECT_EQ(true, any(r2));
	EXPECT_EQ(false, all(r2));

	vint4 a3(1, 2, 3, 4);
	vint4 b3(1, 0, 3, 0);
	vmask4 r3 = a3 == b3;
	EXPECT_EQ(0x5, mask(r3));
	EXPECT_EQ(true, any(r3));
	EXPECT_EQ(false, all(r3));

	vint4 a4(1, 2, 3, 4);
	vmask4 r4 = a4 == a4;
	EXPECT_EQ(0xF, mask(r4));
	EXPECT_EQ(true, any(r4));
	EXPECT_EQ(true, all(r4));
}

/** @brief Test vint4 cne. */
TEST(vint4, cne)
{
	vint4 a1(1, 2, 3, 4);
	vint4 b1(0, 1, 2, 3);
	vmask4 r1 = a1 != b1;
	EXPECT_EQ(0xF, mask(r1));
	EXPECT_EQ(true, any(r1));
	EXPECT_EQ(true, all(r1));

	vint4 a2(1, 2, 3, 4);
	vint4 b2(1, 0, 0, 0);
	vmask4 r2 = a2 != b2;
	EXPECT_EQ(0xE, mask(r2));
	EXPECT_EQ(true, any(r2));
	EXPECT_EQ(false, all(r2));

	vint4 a3(1, 2, 3, 4);
	vint4 b3(1, 0, 3, 0);
	vmask4 r3 = a3 != b3;
	EXPECT_EQ(0xA, mask(r3));
	EXPECT_EQ(true, any(r3));
	EXPECT_EQ(false, all(r3));

	vint4 a4(1, 2, 3, 4);
	vmask4 r4 = a4 != a4;
	EXPECT_EQ(0, mask(r4));
	EXPECT_EQ(false, any(r4));
	EXPECT_EQ(false, all(r4));
}

/** @brief Test vint4 clt. */
TEST(vint4, clt)
{
	vint4 a(1, 2, 3, 4);
	vint4 b(0, 3, 3, 5);
	vmask4 r = a < b;
	EXPECT_EQ(0xA, mask(r));
}

/** @brief Test vint4 cgt. */
TEST(vint4, cle)
{
	vint4 a(1, 2, 3, 4);
	vint4 b(0, 3, 3, 5);
	vmask4 r = a > b;
	EXPECT_EQ(0x1, mask(r));
}

/** @brief Test vint4 lsl. */
TEST(vint4, lsl)
{
	vint4 a(1, 2, 4, 4);
	a = lsl<0>(a);
	EXPECT_EQ(a.lane<0>(), 1);
	EXPECT_EQ(a.lane<1>(), 2);
	EXPECT_EQ(a.lane<2>(), 4);
	EXPECT_EQ(a.lane<3>(), 4);

	a = lsl<1>(a);
	EXPECT_EQ(a.lane<0>(), 2);
	EXPECT_EQ(a.lane<1>(), 4);
	EXPECT_EQ(a.lane<2>(), 8);
	EXPECT_EQ(a.lane<3>(), 8);

	a = lsl<2>(a);
	EXPECT_EQ(a.lane<0>(), 8);
	EXPECT_EQ(a.lane<1>(), 16);
	EXPECT_EQ(a.lane<2>(), 32);
	EXPECT_EQ(a.lane<3>(), 32);
}

/** @brief Test vint4 asr. */
TEST(vint4, asr)
{
	vint4 a(1, 2, 4, -4);
	a = asr<0>(a);
	EXPECT_EQ(a.lane<0>(),  1);
	EXPECT_EQ(a.lane<1>(),  2);
	EXPECT_EQ(a.lane<2>(),  4);
	EXPECT_EQ(a.lane<3>(), -4);

	a = asr<1>(a);
	EXPECT_EQ(a.lane<0>(),  0);
	EXPECT_EQ(a.lane<1>(),  1);
	EXPECT_EQ(a.lane<2>(),  2);
	EXPECT_EQ(a.lane<3>(), -2);

	// Note - quirk of asr is that you will get "stuck" at -1
	a = asr<2>(a);
	EXPECT_EQ(a.lane<0>(),  0);
	EXPECT_EQ(a.lane<1>(),  0);
	EXPECT_EQ(a.lane<2>(),  0);
	EXPECT_EQ(a.lane<3>(), -1);
}

/** @brief Test vint4 min. */
TEST(vint4, min)
{
	vint4 a(1, 2, 3, 4);
	vint4 b(0, 3, 3, 5);
	vint4 r = min(a, b);
	EXPECT_EQ(r.lane<0>(), 0);
	EXPECT_EQ(r.lane<1>(), 2);
	EXPECT_EQ(r.lane<2>(), 3);
	EXPECT_EQ(r.lane<3>(), 4);
}

/** @brief Test vint4 max. */
TEST(vint4, max)
{
	vint4 a(1, 2, 3, 4);
	vint4 b(0, 3, 3, 5);
	vint4 r = max(a, b);
	EXPECT_EQ(r.lane<0>(), 1);
	EXPECT_EQ(r.lane<1>(), 3);
	EXPECT_EQ(r.lane<2>(), 3);
	EXPECT_EQ(r.lane<3>(), 5);
}

/** @brief Test vint4 clamp. */
TEST(vint4, clamp)
{
	vint4 a(1, 2, 3, 4);
	vint4 r = clamp(2, 3, a);
	EXPECT_EQ(r.lane<0>(), 2);
	EXPECT_EQ(r.lane<1>(), 2);
	EXPECT_EQ(r.lane<2>(), 3);
	EXPECT_EQ(r.lane<3>(), 3);
}

/** @brief Test vint4 hmin. */
TEST(vint4, hmin)
{
	vint4 a1(1, 2, 1, 2);
	vint4 r1 = hmin(a1);
	EXPECT_EQ(r1.lane<0>(), 1);
	EXPECT_EQ(r1.lane<1>(), 1);
	EXPECT_EQ(r1.lane<2>(), 1);
	EXPECT_EQ(r1.lane<3>(), 1);

	vint4 a2(1, 2, -1, 5);
	vint4 r2 = hmin(a2);
	EXPECT_EQ(r2.lane<0>(), -1);
	EXPECT_EQ(r2.lane<1>(), -1);
	EXPECT_EQ(r2.lane<2>(), -1);
	EXPECT_EQ(r2.lane<3>(), -1);
}

/** @brief Test vint4 hmax. */
TEST(vint4, hmax)
{
	vint4 a1(1, 3, 1, 2);
	vint4 r1 = hmax(a1);
	EXPECT_EQ(r1.lane<0>(), 3);
	EXPECT_EQ(r1.lane<1>(), 3);
	EXPECT_EQ(r1.lane<2>(), 3);
	EXPECT_EQ(r1.lane<3>(), 3);

	vint4 a2(1, 2, -1, 5);
	vint4 r2 = hmax(a2);
	EXPECT_EQ(r2.lane<0>(), 5);
	EXPECT_EQ(r2.lane<1>(), 5);
	EXPECT_EQ(r2.lane<2>(), 5);
	EXPECT_EQ(r2.lane<3>(), 5);
}

/** @brief Test vint4 hadd_s. */
TEST(vint4, hadd_s)
{
	vint4 a1(1, 3, 5, 7);
	int r1 = hadd_s(a1);
	EXPECT_EQ(r1, 16);

	vint4 a2(1, 2, -1, 5);
	int r2 = hadd_s(a2);
	EXPECT_EQ(r2, 7);
}

/** @brief Test vint4 hadd_s. */
TEST(vint4, hadd_rgb_s)
{
	vint4 a1(1, 3, 5, 7);
	int r1 = hadd_rgb_s(a1);
	EXPECT_EQ(r1, 9);

	vint4 a2(1, 2, -1, 5);
	int r2 = hadd_rgb_s(a2);
	EXPECT_EQ(r2, 2);
}

/** @brief Test vint4 storea. */
TEST(vint4, storea)
{
	alignas(16) int out[4];
	vint4 a(s32_data);
	storea(a, out);
	EXPECT_EQ(out[0], 0);
	EXPECT_EQ(out[1], 1);
	EXPECT_EQ(out[2], 2);
	EXPECT_EQ(out[3], 3);
}

/** @brief Test vint4 store. */
TEST(vint4, store)
{
	alignas(16) int out[5];
	vint4 a(s32_data);
	store(a, &(out[1]));
	EXPECT_EQ(out[1], 0);
	EXPECT_EQ(out[2], 1);
	EXPECT_EQ(out[3], 2);
	EXPECT_EQ(out[4], 3);
}

/** @brief Test vint4 store_nbytes. */
TEST(vint4, store_nbytes)
{
	alignas(16) int out;
	vint4 a(42, 314, 75, 90);
	store_nbytes(a, (uint8_t*)&out);
	EXPECT_EQ(out, 42);
}

/** @brief Test vint4 gatheri. */
TEST(vint4, gatheri)
{
	vint4 indices(0, 4, 3, 2);
	vint4 r = gatheri(s32_data, indices);
	EXPECT_EQ(r.lane<0>(), 0);
	EXPECT_EQ(r.lane<1>(), 4);
	EXPECT_EQ(r.lane<2>(), 3);
	EXPECT_EQ(r.lane<3>(), 2);
}

/** @brief Test vint4 pack_low_bytes. */
TEST(vint4, pack_low_bytes)
{
	vint4 a(1, 2, 3, 4);
	vint4 r = pack_low_bytes(a);
	EXPECT_EQ(r.lane<0>(), (4 << 24) | (3 << 16) | (2  << 8) | (1 << 0));
}

/** @brief Test vint4 select. */
TEST(vint4, select)
{
	vint4 m1(1, 1, 1, 1);
	vint4 m2(1, 2, 1, 2);
	vmask4 cond = m1 == m2;

	vint4 a(1, 3, 3, 1);
	vint4 b(4, 2, 2, 4);

	vint4 r1 = select(a, b, cond);
	EXPECT_EQ(r1.lane<0>(), 4);
	EXPECT_EQ(r1.lane<1>(), 3);
	EXPECT_EQ(r1.lane<2>(), 2);
	EXPECT_EQ(r1.lane<3>(), 1);

	vint4 r2 = select(b, a, cond);
	EXPECT_EQ(r2.lane<0>(), 1);
	EXPECT_EQ(r2.lane<1>(), 2);
	EXPECT_EQ(r2.lane<2>(), 3);
	EXPECT_EQ(r2.lane<3>(), 4);
}

/** @brief Test vint4 select MSB. */
TEST(vint4, select_msb)
{
	vint4 msb(0x80000000, 0, 0x80000000, 0);
	vmask4 cond(msb.m);

	vint4 a(1, 3, 3, 1);
	vint4 b(4, 2, 2, 4);

	vint4 r1 = select(a, b, cond);
	EXPECT_EQ(r1.lane<0>(), 4);
	EXPECT_EQ(r1.lane<1>(), 3);
	EXPECT_EQ(r1.lane<2>(), 2);
	EXPECT_EQ(r1.lane<3>(), 1);

	vint4 r2 = select(b, a, cond);
	EXPECT_EQ(r2.lane<0>(), 1);
	EXPECT_EQ(r2.lane<1>(), 2);
	EXPECT_EQ(r2.lane<2>(), 3);
	EXPECT_EQ(r2.lane<3>(), 4);
}

// VMASK4 tests - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/** @brief Test vmask4 scalar constructor. */
TEST(vmask4, scalar_construct)
{
	vfloat4 m1a(0, 0, 0, 0);
	vfloat4 m1b(1, 1, 1, 1);
	vmask4 m1(true, false, true, false);

	vfloat4 r = select(m1a, m1b, m1);

	EXPECT_EQ(r.lane<0>(), 1);
	EXPECT_EQ(r.lane<1>(), 0);
	EXPECT_EQ(r.lane<2>(), 1);
	EXPECT_EQ(r.lane<3>(), 0);
}

/** @brief Test vmask4 or. */
TEST(vmask4, or)
{
	vfloat4 m1a(0, 1, 0, 1);
	vfloat4 m1b(1, 1, 1, 1);
	vmask4 m1 = m1a == m1b;

	vfloat4 m2a(1, 1, 0, 0);
	vfloat4 m2b(1, 1, 1, 1);
	vmask4 m2 = m2a == m2b;

	vmask4 r = m1 | m2;
	EXPECT_EQ(mask(r), 0xB);
}

/** @brief Test vmask4 and. */
TEST(vmask4, and)
{
	vfloat4 m1a(0, 1, 0, 1);
	vfloat4 m1b(1, 1, 1, 1);
	vmask4 m1 = m1a == m1b;

	vfloat4 m2a(1, 1, 0, 0);
	vfloat4 m2b(1, 1, 1, 1);
	vmask4 m2 = m2a == m2b;

	vmask4 r = m1 & m2;
	EXPECT_EQ(mask(r), 0x2);
}

/** @brief Test vmask4 xor. */
TEST(vmask4, xor)
{
	vfloat4 m1a(0, 1, 0, 1);
	vfloat4 m1b(1, 1, 1, 1);
	vmask4 m1 = m1a == m1b;

	vfloat4 m2a(1, 1, 0, 0);
	vfloat4 m2b(1, 1, 1, 1);
	vmask4 m2 = m2a == m2b;

	vmask4 r = m1 ^ m2;
	EXPECT_EQ(mask(r), 0x9);
}

/** @brief Test vmask4 not. */
TEST(vmask4, not)
{
	vfloat4 m1a(0, 1, 0, 1);
	vfloat4 m1b(1, 1, 1, 1);
	vmask4 m1 = m1a == m1b;
	vmask4 r = ~m1;
	EXPECT_EQ(mask(r), 0x5);
}

# if ASTCENC_SIMD_WIDTH == 8

// VFLOAT8 tests - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

/** @brief Test unaligned vfloat8 data load. */
TEST(vfloat8, UnalignedLoad)
{
	vfloat8 a(&(f32_data[1]));
	EXPECT_EQ(a.lane<0>(), 1.0f);
	EXPECT_EQ(a.lane<1>(), 2.0f);
	EXPECT_EQ(a.lane<2>(), 3.0f);
	EXPECT_EQ(a.lane<3>(), 4.0f);
	EXPECT_EQ(a.lane<4>(), 5.0f);
	EXPECT_EQ(a.lane<5>(), 6.0f);
	EXPECT_EQ(a.lane<6>(), 7.0f);
	EXPECT_EQ(a.lane<7>(), 8.0f);
}

/** @brief Test scalar duplicated vfloat8 load. */
TEST(vfloat8, ScalarDupLoad)
{
	vfloat8 a(1.1f);
	EXPECT_EQ(a.lane<0>(), 1.1f);
	EXPECT_EQ(a.lane<1>(), 1.1f);
	EXPECT_EQ(a.lane<2>(), 1.1f);
	EXPECT_EQ(a.lane<3>(), 1.1f);
	EXPECT_EQ(a.lane<4>(), 1.1f);
	EXPECT_EQ(a.lane<5>(), 1.1f);
	EXPECT_EQ(a.lane<6>(), 1.1f);
	EXPECT_EQ(a.lane<7>(), 1.1f);
}

/** @brief Test scalar vfloat8 load. */
TEST(vfloat8, ScalarLoad)
{
	vfloat8 a(1.1f, 2.2f, 3.3f, 4.4f, 5.5f, 6.6f, 7.7f, 8.8f);
	EXPECT_EQ(a.lane<0>(), 1.1f);
	EXPECT_EQ(a.lane<1>(), 2.2f);
	EXPECT_EQ(a.lane<2>(), 3.3f);
	EXPECT_EQ(a.lane<3>(), 4.4f);
	EXPECT_EQ(a.lane<4>(), 5.5f);
	EXPECT_EQ(a.lane<5>(), 6.6f);
	EXPECT_EQ(a.lane<6>(), 7.7f);
	EXPECT_EQ(a.lane<7>(), 8.8f);
}

/** @brief Test copy vfloat8 load. */
TEST(vfloat8, CopyLoad)
{
	vfloat8 s(1.1f, 2.2f, 3.3f, 4.4f, 5.5f, 6.6f, 7.7f, 8.8f);
	vfloat8 a(s.m);
	EXPECT_EQ(a.lane<0>(), 1.1f);
	EXPECT_EQ(a.lane<1>(), 2.2f);
	EXPECT_EQ(a.lane<2>(), 3.3f);
	EXPECT_EQ(a.lane<3>(), 4.4f);
	EXPECT_EQ(a.lane<4>(), 5.5f);
	EXPECT_EQ(a.lane<5>(), 6.6f);
	EXPECT_EQ(a.lane<6>(), 7.7f);
	EXPECT_EQ(a.lane<7>(), 8.8f);
}

/** @brief Test vfloat8 zero. */
TEST(vfloat8, Zero)
{
	vfloat8 a = vfloat8::zero();
	EXPECT_EQ(a.lane<0>(), 0.0f);
	EXPECT_EQ(a.lane<1>(), 0.0f);
	EXPECT_EQ(a.lane<2>(), 0.0f);
	EXPECT_EQ(a.lane<3>(), 0.0f);
	EXPECT_EQ(a.lane<4>(), 0.0f);
	EXPECT_EQ(a.lane<5>(), 0.0f);
	EXPECT_EQ(a.lane<6>(), 0.0f);
	EXPECT_EQ(a.lane<7>(), 0.0f);
}

/** @brief Test vfloat8 load1. */
TEST(vfloat8, Load1)
{
	float s = 3.14f;
	vfloat8 a = vfloat8::load1(&s);
	EXPECT_EQ(a.lane<0>(), 3.14f);
	EXPECT_EQ(a.lane<1>(), 3.14f);
	EXPECT_EQ(a.lane<2>(), 3.14f);
	EXPECT_EQ(a.lane<3>(), 3.14f);
	EXPECT_EQ(a.lane<4>(), 3.14f);
	EXPECT_EQ(a.lane<5>(), 3.14f);
	EXPECT_EQ(a.lane<6>(), 3.14f);
	EXPECT_EQ(a.lane<7>(), 3.14f);
}

/** @brief Test vfloat8 loada. */
TEST(vfloat8, Loada)
{
	vfloat8 a(&(f32_data[0]));
	EXPECT_EQ(a.lane<0>(), 0.0f);
	EXPECT_EQ(a.lane<1>(), 1.0f);
	EXPECT_EQ(a.lane<2>(), 2.0f);
	EXPECT_EQ(a.lane<3>(), 3.0f);
	EXPECT_EQ(a.lane<4>(), 4.0f);
	EXPECT_EQ(a.lane<5>(), 5.0f);
	EXPECT_EQ(a.lane<6>(), 6.0f);
	EXPECT_EQ(a.lane<7>(), 7.0f);
}

/** @brief Test vfloat8 lane_id. */
TEST(vfloat8, LaneID)
{
	vfloat8 a = vfloat8::lane_id();
	EXPECT_EQ(a.lane<0>(), 0.0f);
	EXPECT_EQ(a.lane<1>(), 1.0f);
	EXPECT_EQ(a.lane<2>(), 2.0f);
	EXPECT_EQ(a.lane<3>(), 3.0f);
	EXPECT_EQ(a.lane<4>(), 4.0f);
	EXPECT_EQ(a.lane<5>(), 5.0f);
	EXPECT_EQ(a.lane<6>(), 6.0f);
	EXPECT_EQ(a.lane<7>(), 7.0f);
}

/** @brief Test vfloat8 add. */
TEST(vfloat8, vadd)
{
	vfloat8 a(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
	vfloat8 b(0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f);
	a = a + b;
	EXPECT_EQ(a.lane<0>(), 1.0f + 0.1f);
	EXPECT_EQ(a.lane<1>(), 2.0f + 0.2f);
	EXPECT_EQ(a.lane<2>(), 3.0f + 0.3f);
	EXPECT_EQ(a.lane<3>(), 4.0f + 0.4f);
	EXPECT_EQ(a.lane<4>(), 5.0f + 0.5f);
	EXPECT_EQ(a.lane<5>(), 6.0f + 0.6f);
	EXPECT_EQ(a.lane<6>(), 7.0f + 0.7f);
	EXPECT_EQ(a.lane<7>(), 8.0f + 0.8f);
}

/** @brief Test vfloat8 sub. */
TEST(vfloat8, vsub)
{
	vfloat8 a(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
	vfloat8 b(0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f);
	a = a - b;
	EXPECT_EQ(a.lane<0>(), 1.0f - 0.1f);
	EXPECT_EQ(a.lane<1>(), 2.0f - 0.2f);
	EXPECT_EQ(a.lane<2>(), 3.0f - 0.3f);
	EXPECT_EQ(a.lane<3>(), 4.0f - 0.4f);
	EXPECT_EQ(a.lane<4>(), 5.0f - 0.5f);
	EXPECT_EQ(a.lane<5>(), 6.0f - 0.6f);
	EXPECT_EQ(a.lane<6>(), 7.0f - 0.7f);
	EXPECT_EQ(a.lane<7>(), 8.0f - 0.8f);
}

/** @brief Test vfloat8 mul. */
TEST(vfloat8, vmul)
{
	vfloat8 a(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
	vfloat8 b(0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f);
	a = a * b;
	EXPECT_EQ(a.lane<0>(), 1.0f * 0.1f);
	EXPECT_EQ(a.lane<1>(), 2.0f * 0.2f);
	EXPECT_EQ(a.lane<2>(), 3.0f * 0.3f);
	EXPECT_EQ(a.lane<3>(), 4.0f * 0.4f);
	EXPECT_EQ(a.lane<4>(), 5.0f * 0.5f);
	EXPECT_EQ(a.lane<5>(), 6.0f * 0.6f);
	EXPECT_EQ(a.lane<6>(), 7.0f * 0.7f);
	EXPECT_EQ(a.lane<7>(), 8.0f * 0.8f);
}

/** @brief Test vfloat8 mul. */
TEST(vfloat8, vsmul)
{
	vfloat8 a(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
	float b = 3.14f;
	a = a * b;
	EXPECT_EQ(a.lane<0>(), 1.0f * 3.14f);
	EXPECT_EQ(a.lane<1>(), 2.0f * 3.14f);
	EXPECT_EQ(a.lane<2>(), 3.0f * 3.14f);
	EXPECT_EQ(a.lane<3>(), 4.0f * 3.14f);
	EXPECT_EQ(a.lane<4>(), 5.0f * 3.14f);
	EXPECT_EQ(a.lane<5>(), 6.0f * 3.14f);
	EXPECT_EQ(a.lane<6>(), 7.0f * 3.14f);
	EXPECT_EQ(a.lane<7>(), 8.0f * 3.14f);
}

/** @brief Test vfloat8 mul. */
TEST(vfloat8, svmul)
{
	float a = 3.14f;
	vfloat8 b(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
	b = a * b;
	EXPECT_EQ(b.lane<0>(), 3.14f * 1.0f);
	EXPECT_EQ(b.lane<1>(), 3.14f * 2.0f);
	EXPECT_EQ(b.lane<2>(), 3.14f * 3.0f);
	EXPECT_EQ(b.lane<3>(), 3.14f * 4.0f);
	EXPECT_EQ(b.lane<4>(), 3.14f * 5.0f);
	EXPECT_EQ(b.lane<5>(), 3.14f * 6.0f);
	EXPECT_EQ(b.lane<6>(), 3.14f * 7.0f);
	EXPECT_EQ(b.lane<7>(), 3.14f * 8.0f);
}

/** @brief Test vfloat8 div. */
TEST(vfloat8, vdiv)
{
	vfloat8 a(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
	vfloat8 b(0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f);
	a = a / b;
	EXPECT_EQ(a.lane<0>(), 1.0f / 0.1f);
	EXPECT_EQ(a.lane<1>(), 2.0f / 0.2f);
	EXPECT_EQ(a.lane<2>(), 3.0f / 0.3f);
	EXPECT_EQ(a.lane<3>(), 4.0f / 0.4f);
	EXPECT_EQ(a.lane<4>(), 5.0f / 0.5f);
	EXPECT_EQ(a.lane<5>(), 6.0f / 0.6f);
	EXPECT_EQ(a.lane<6>(), 7.0f / 0.7f);
	EXPECT_EQ(a.lane<7>(), 8.0f / 0.8f);
}

/** @brief Test vfloat8 div. */
TEST(vfloat8, vsdiv)
{
	vfloat8 a(0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f);
	float b = 3.14f;
	vfloat8 r = a / b;

	EXPECT_EQ(r.lane<0>(), 0.1f / 3.14f);
	EXPECT_EQ(r.lane<1>(), 0.2f / 3.14f);
	EXPECT_EQ(r.lane<2>(), 0.3f / 3.14f);
	EXPECT_EQ(r.lane<3>(), 0.4f / 3.14f);
	EXPECT_EQ(r.lane<4>(), 0.5f / 3.14f);
	EXPECT_EQ(r.lane<5>(), 0.6f / 3.14f);
	EXPECT_EQ(r.lane<6>(), 0.7f / 3.14f);
	EXPECT_EQ(r.lane<7>(), 0.8f / 3.14f);
}

/** @brief Test vfloat8 div. */
TEST(vfloat8, svdiv)
{
	float a = 3.14f;
	vfloat8 b(0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f);
	vfloat8 r = a / b;

	EXPECT_EQ(r.lane<0>(), 3.14f / 0.1f);
	EXPECT_EQ(r.lane<1>(), 3.14f / 0.2f);
	EXPECT_EQ(r.lane<2>(), 3.14f / 0.3f);
	EXPECT_EQ(r.lane<3>(), 3.14f / 0.4f);
	EXPECT_EQ(r.lane<4>(), 3.14f / 0.5f);
	EXPECT_EQ(r.lane<5>(), 3.14f / 0.6f);
	EXPECT_EQ(r.lane<6>(), 3.14f / 0.7f);
	EXPECT_EQ(r.lane<7>(), 3.14f / 0.8f);
}

/** @brief Test vfloat8 ceq. */
TEST(vfloat8, ceq)
{
	vfloat8 a1(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
	vfloat8 b1(0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f);
	vmask8 r1 = a1 == b1;
	EXPECT_EQ(0, mask(r1));
	EXPECT_EQ(false, any(r1));
	EXPECT_EQ(false, all(r1));

	vfloat8 a2(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
	vfloat8 b2(1.0f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f);
	vmask8 r2 = a2 == b2;
	EXPECT_EQ(0x1, mask(r2));
	EXPECT_EQ(true, any(r2));
	EXPECT_EQ(false, all(r2));

	vfloat8 a3(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
	vfloat8 b3(1.0f, 0.2f, 3.0f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f);
	vmask8 r3 = a3 == b3;
	EXPECT_EQ(0x5, mask(r3));
	EXPECT_EQ(true, any(r3));
	EXPECT_EQ(false, all(r3));

	vfloat8 a4(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
	vmask8 r4 = a4 == a4;
	EXPECT_EQ(0xFF, mask(r4));
	EXPECT_EQ(true, any(r4));
	EXPECT_EQ(true, all(r4));
}

/** @brief Test vfloat8 cne. */
TEST(vfloat8, cne)
{
	vfloat8 a1(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
	vfloat8 b1(0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f);
	vmask8 r1 = a1 != b1;
	EXPECT_EQ(0xFF, mask(r1));
	EXPECT_EQ(true, any(r1));
	EXPECT_EQ(true, all(r1));

	vfloat8 a2(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
	vfloat8 b2(1.0f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f);
	vmask8 r2 = a2 != b2;
	EXPECT_EQ(0xFE, mask(r2));
	EXPECT_EQ(true, any(r2));
	EXPECT_EQ(false, all(r2));

	vfloat8 a3(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
	vfloat8 b3(1.0f, 0.2f, 3.0f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f);
	vmask8 r3 = a3 != b3;
	EXPECT_EQ(0xFA, mask(r3));
	EXPECT_EQ(true, any(r3));
	EXPECT_EQ(false, all(r3));

	vfloat8 a4(1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f);
	vmask8 r4 = a4 != a4;
	EXPECT_EQ(0, mask(r4));
	EXPECT_EQ(false, any(r4));
	EXPECT_EQ(false, all(r4));
}

/** @brief Test vfloat8 clt. */
TEST(vfloat8, clt)
{
	vfloat8 a(1.0f, 2.0f, 3.0f, 4.0f, 1.0f, 2.0f, 3.0f, 4.0f);
	vfloat8 b(0.9f, 2.1f, 3.0f, 4.1f, 0.9f, 2.1f, 3.0f, 4.1f);
	vmask8 r = a < b;
	EXPECT_EQ(0xAA, mask(r));
}

/** @brief Test vfloat8 cle. */
TEST(vfloat8, cle)
{
	vfloat8 a(1.0f, 2.0f, 3.0f, 4.0f, 1.0f, 2.0f, 3.0f, 4.0f);
	vfloat8 b(0.9f, 2.1f, 3.0f, 4.1f, 0.9f, 2.1f, 3.0f, 4.1f);
	vmask8 r = a <= b;
	EXPECT_EQ(0xEE, mask(r));
}

/** @brief Test vfloat8 cgt. */
TEST(vfloat8, cgt)
{
	vfloat8 a(1.0f, 2.0f, 3.0f, 4.0f, 1.0f, 2.0f, 3.0f, 4.0f);
	vfloat8 b(0.9f, 2.1f, 3.0f, 4.1f, 0.9f, 2.1f, 3.0f, 4.1f);
	vmask8 r = a > b;
	EXPECT_EQ(0x11, mask(r));
}

/** @brief Test vfloat8 cge. */
TEST(vfloat8, cge)
{
	vfloat8 a(1.0f, 2.0f, 3.0f, 4.0f, 1.0f, 2.0f, 3.0f, 4.0f);
	vfloat8 b(0.9f, 2.1f, 3.0f, 4.1f, 0.9f, 2.1f, 3.0f, 4.1f);
	vmask8 r = a >= b;
	EXPECT_EQ(0x55, mask(r));
}

/** @brief Test vfloat8 min. */
TEST(vfloat8, min)
{
	vfloat8 a(1.0f, 2.0f, 3.0f, 4.0f, 1.0f, 2.0f, 3.0f, 4.0f);
	vfloat8 b(0.9f, 2.1f, 3.0f, 4.1f, 0.9f, 2.1f, 3.0f, 4.1f);
	vfloat8 r = min(a, b);
	EXPECT_EQ(r.lane<0>(), 0.9f);
	EXPECT_EQ(r.lane<1>(), 2.0f);
	EXPECT_EQ(r.lane<2>(), 3.0f);
	EXPECT_EQ(r.lane<3>(), 4.0f);
	EXPECT_EQ(r.lane<4>(), 0.9f);
	EXPECT_EQ(r.lane<5>(), 2.0f);
	EXPECT_EQ(r.lane<6>(), 3.0f);
	EXPECT_EQ(r.lane<7>(), 4.0f);
}

/** @brief Test vfloat8 max. */
TEST(vfloat8, max)
{
	vfloat8 a(1.0f, 2.0f, 3.0f, 4.0f, 1.0f, 2.0f, 3.0f, 4.0f);
	vfloat8 b(0.9f, 2.1f, 3.0f, 4.1f, 0.9f, 2.1f, 3.0f, 4.1f);
	vfloat8 r = max(a, b);
	EXPECT_EQ(r.lane<0>(), 1.0f);
	EXPECT_EQ(r.lane<1>(), 2.1f);
	EXPECT_EQ(r.lane<2>(), 3.0f);
	EXPECT_EQ(r.lane<3>(), 4.1f);
	EXPECT_EQ(r.lane<4>(), 1.0f);
	EXPECT_EQ(r.lane<5>(), 2.1f);
	EXPECT_EQ(r.lane<6>(), 3.0f);
	EXPECT_EQ(r.lane<7>(), 4.1f);
}

/** @brief Test vfloat8 clamp. */
TEST(vfloat8, clamp)
{
	vfloat8 a1(1.0f, 2.0f, 3.0f, 4.0f, 1.0f, 2.0f, 3.0f, 4.0f);
	vfloat8 r1 = clamp(2.1f, 3.0f, a1);
	EXPECT_EQ(r1.lane<0>(), 2.1f);
	EXPECT_EQ(r1.lane<1>(), 2.1f);
	EXPECT_EQ(r1.lane<2>(), 3.0f);
	EXPECT_EQ(r1.lane<3>(), 3.0f);
	EXPECT_EQ(r1.lane<4>(), 2.1f);
	EXPECT_EQ(r1.lane<5>(), 2.1f);
	EXPECT_EQ(r1.lane<6>(), 3.0f);
	EXPECT_EQ(r1.lane<7>(), 3.0f);

	vfloat8 a2(1.0f, 2.0f, qnan, 4.0f, 1.0f, 2.0f, qnan, 4.0f);
	vfloat8 r2 = clamp(2.1f, 3.0f, a2);
	EXPECT_EQ(r2.lane<0>(), 2.1f);
	EXPECT_EQ(r2.lane<1>(), 2.1f);
	EXPECT_EQ(r2.lane<2>(), 2.1f);
	EXPECT_EQ(r2.lane<3>(), 3.0f);
	EXPECT_EQ(r2.lane<4>(), 2.1f);
	EXPECT_EQ(r2.lane<5>(), 2.1f);
	EXPECT_EQ(r2.lane<6>(), 2.1f);
	EXPECT_EQ(r2.lane<7>(), 3.0f);
}

/** @brief Test vfloat8 clampz. */
TEST(vfloat8, clampz)
{
	vfloat8 a1(-1.0f, 0.0f, 0.1f, 4.0f, -1.0f, 0.0f, 0.1f, 4.0f);
	vfloat8 r1 = clampz(3.0f, a1);
	EXPECT_EQ(r1.lane<0>(), 0.0f);
	EXPECT_EQ(r1.lane<1>(), 0.0f);
	EXPECT_EQ(r1.lane<2>(), 0.1f);
	EXPECT_EQ(r1.lane<3>(), 3.0f);
	EXPECT_EQ(r1.lane<4>(), 0.0f);
	EXPECT_EQ(r1.lane<5>(), 0.0f);
	EXPECT_EQ(r1.lane<6>(), 0.1f);
	EXPECT_EQ(r1.lane<7>(), 3.0f);

	vfloat8 a2(-1.0f, 0.0f, qnan, 4.0f, -1.0f, 0.0f, qnan, 4.0f);
	vfloat8 r2 = clampz(3.0f, a2);
	EXPECT_EQ(r2.lane<0>(), 0.0f);
	EXPECT_EQ(r2.lane<1>(), 0.0f);
	EXPECT_EQ(r2.lane<2>(), 0.0f);
	EXPECT_EQ(r2.lane<3>(), 3.0f);
	EXPECT_EQ(r2.lane<4>(), 0.0f);
	EXPECT_EQ(r2.lane<5>(), 0.0f);
	EXPECT_EQ(r2.lane<6>(), 0.0f);
	EXPECT_EQ(r2.lane<7>(), 3.0f);
}

/** @brief Test vfloat8 clampz. */
TEST(vfloat8, clampzo)
{
	vfloat8 a1(-1.0f, 0.0f, 0.1f, 4.0f, -1.0f, 0.0f, 0.1f, 4.0f);
	vfloat8 r1 = clampzo(a1);
	EXPECT_EQ(r1.lane<0>(), 0.0f);
	EXPECT_EQ(r1.lane<1>(), 0.0f);
	EXPECT_EQ(r1.lane<2>(), 0.1f);
	EXPECT_EQ(r1.lane<3>(), 1.0f);
	EXPECT_EQ(r1.lane<4>(), 0.0f);
	EXPECT_EQ(r1.lane<5>(), 0.0f);
	EXPECT_EQ(r1.lane<6>(), 0.1f);
	EXPECT_EQ(r1.lane<7>(), 1.0f);

	vfloat8 a2(-1.0f, 0.0f, qnan, 4.0f, -1.0f, 0.0f, qnan, 4.0f);
	vfloat8 r2 = clampzo(a2);
	EXPECT_EQ(r2.lane<0>(), 0.0f);
	EXPECT_EQ(r2.lane<1>(), 0.0f);
	EXPECT_EQ(r2.lane<2>(), 0.0f);
	EXPECT_EQ(r2.lane<3>(), 1.0f);
	EXPECT_EQ(r2.lane<4>(), 0.0f);
	EXPECT_EQ(r2.lane<5>(), 0.0f);
	EXPECT_EQ(r2.lane<6>(), 0.0f);
	EXPECT_EQ(r2.lane<7>(), 1.0f);
}

/** @brief Test vfloat8 abs. */
TEST(vfloat8, abs)
{
	vfloat8 a(-1.0f, 0.0f, 0.1f, 4.0f, -1.0f, 0.0f, 0.1f, 4.0f);
	vfloat8 r = abs(a);
	EXPECT_EQ(r.lane<0>(), 1.0f);
	EXPECT_EQ(r.lane<1>(), 0.0f);
	EXPECT_EQ(r.lane<2>(), 0.1f);
	EXPECT_EQ(r.lane<3>(), 4.0f);
	EXPECT_EQ(r.lane<4>(), 1.0f);
	EXPECT_EQ(r.lane<5>(), 0.0f);
	EXPECT_EQ(r.lane<6>(), 0.1f);
	EXPECT_EQ(r.lane<7>(), 4.0f);
}

/** @brief Test vfloat8 round. */
TEST(vfloat8, round)
{
	vfloat8 a(1.1f, 1.5f, 1.6f, 4.0f, 1.1f, 1.5f, 1.6f, 4.0f);
	vfloat8 r = round(a);
	EXPECT_EQ(r.lane<0>(), 1.0f);
	EXPECT_EQ(r.lane<1>(), 2.0f);
	EXPECT_EQ(r.lane<2>(), 2.0f);
	EXPECT_EQ(r.lane<3>(), 4.0f);
	EXPECT_EQ(r.lane<4>(), 1.0f);
	EXPECT_EQ(r.lane<5>(), 2.0f);
	EXPECT_EQ(r.lane<6>(), 2.0f);
	EXPECT_EQ(r.lane<7>(), 4.0f);
}

/** @brief Test vfloat8 hmin. */
TEST(vfloat8, hmin)
{
	vfloat8 a1(1.1f, 1.5f, 1.6f, 4.0f, 1.1f, 1.5f, 1.6f, 4.0f);
	vfloat8 r1 = hmin(a1);
	EXPECT_EQ(r1.lane<0>(), 1.1f);
	EXPECT_EQ(r1.lane<1>(), 1.1f);
	EXPECT_EQ(r1.lane<2>(), 1.1f);
	EXPECT_EQ(r1.lane<3>(), 1.1f);
	EXPECT_EQ(r1.lane<4>(), 1.1f);
	EXPECT_EQ(r1.lane<5>(), 1.1f);
	EXPECT_EQ(r1.lane<6>(), 1.1f);
	EXPECT_EQ(r1.lane<7>(), 1.1f);

	vfloat8 a2(1.1f, 1.5f, 1.6f, 0.2f, 1.1f, 1.5f, 1.6f, 0.2f);
	vfloat8 r2 = hmin(a2);
	EXPECT_EQ(r2.lane<0>(), 0.2f);
	EXPECT_EQ(r2.lane<1>(), 0.2f);
	EXPECT_EQ(r2.lane<2>(), 0.2f);
	EXPECT_EQ(r2.lane<3>(), 0.2f);
	EXPECT_EQ(r2.lane<4>(), 0.2f);
	EXPECT_EQ(r2.lane<5>(), 0.2f);
	EXPECT_EQ(r2.lane<6>(), 0.2f);
	EXPECT_EQ(r2.lane<7>(), 0.2f);
}

/** @brief Test vfloat8 hmin_s. */
TEST(vfloat8, hmin_s)
{
	vfloat8 a1(1.1f, 1.5f, 1.6f, 4.0f, 1.1f, 1.5f, 1.6f, 4.0f);
	float r1 = hmin_s(a1);
	EXPECT_EQ(r1, 1.1f);

	vfloat8 a2(1.1f, 1.5f, 1.6f, 0.2f, 1.1f, 1.5f, 1.6f, 0.2f);
	float r2 = hmin_s(a2);
	EXPECT_EQ(r2, 0.2f);
}

/** @brief Test vfloat8 hmax. */
TEST(vfloat8, hmax)
{
	vfloat8 a1(1.1f, 1.5f, 1.6f, 4.0f, 1.1f, 1.5f, 1.6f, 4.0f);
	vfloat8 r1 = hmax(a1);
	EXPECT_EQ(r1.lane<0>(), 4.0f);
	EXPECT_EQ(r1.lane<1>(), 4.0f);
	EXPECT_EQ(r1.lane<2>(), 4.0f);
	EXPECT_EQ(r1.lane<3>(), 4.0f);
	EXPECT_EQ(r1.lane<4>(), 4.0f);
	EXPECT_EQ(r1.lane<5>(), 4.0f);
	EXPECT_EQ(r1.lane<6>(), 4.0f);
	EXPECT_EQ(r1.lane<7>(), 4.0f);

	vfloat8 a2(1.1f, 1.5f, 1.6f, 0.2f, 1.1f, 1.5f, 1.6f, 0.2f);
	vfloat8 r2 = hmax(a2);
	EXPECT_EQ(r2.lane<0>(), 1.6f);
	EXPECT_EQ(r2.lane<1>(), 1.6f);
	EXPECT_EQ(r2.lane<2>(), 1.6f);
	EXPECT_EQ(r2.lane<3>(), 1.6f);
	EXPECT_EQ(r2.lane<4>(), 1.6f);
	EXPECT_EQ(r2.lane<5>(), 1.6f);
	EXPECT_EQ(r2.lane<6>(), 1.6f);
	EXPECT_EQ(r2.lane<7>(), 1.6f);
}

/** @brief Test vfloat8 hmax_s. */
TEST(vfloat8, hmax_s)
{
	vfloat8 a1(1.1f, 1.5f, 1.6f, 4.0f, 1.1f, 1.5f, 1.6f, 4.0f);
	float r1 = hmax_s(a1);
	EXPECT_EQ(r1, 4.0f);

	vfloat8 a2(1.1f, 1.5f, 1.6f, 0.2f, 1.1f, 1.5f, 1.6f, 0.2f);
	float r2 = hmax_s(a2);
	EXPECT_EQ(r2, 1.6f);
}

/** @brief Test vfloat8 hadd_s. */
TEST(vfloat8, hadd_s)
{
	vfloat8 a1(1.1f, 1.5f, 1.6f, 4.0f, 1.1f, 1.5f, 1.6f, 4.0f);
	float sum = 1.1f + 1.5f + 1.6f + 4.0f + 1.1f + 1.5f + 1.6f + 4.0f;
	float r = hadd_s(a1);
	EXPECT_NEAR(r, sum, 0.005f);
}

/** @brief Test vfloat8 sqrt. */
TEST(vfloat8, sqrt)
{
	vfloat8 a(1.0f, 2.0f, 3.0f, 4.0f, 1.0f, 2.0f, 3.0f, 4.0f);
	vfloat8 r = sqrt(a);
	EXPECT_EQ(r.lane<0>(), std::sqrt(1.0f));
	EXPECT_EQ(r.lane<1>(), std::sqrt(2.0f));
	EXPECT_EQ(r.lane<2>(), std::sqrt(3.0f));
	EXPECT_EQ(r.lane<3>(), std::sqrt(4.0f));
	EXPECT_EQ(r.lane<4>(), std::sqrt(1.0f));
	EXPECT_EQ(r.lane<5>(), std::sqrt(2.0f));
	EXPECT_EQ(r.lane<6>(), std::sqrt(3.0f));
	EXPECT_EQ(r.lane<7>(), std::sqrt(4.0f));
}

/** @brief Test vfloat8 select. */
TEST(vfloat8, select)
{
	vfloat8 m1(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f);
	vfloat8 m2(1.0f, 2.0f, 1.0f, 2.0f, 1.0f, 2.0f, 1.0f, 2.0f);
	vmask8 cond = m1 == m2;

	vfloat8 a(1.0f, 3.0f, 3.0f, 1.0f, 1.0f, 3.0f, 3.0f, 1.0);
	vfloat8 b(4.0f, 2.0f, 2.0f, 4.0f, 4.0f, 2.0f, 2.0f, 4.0);

	// Select in one direction
	vfloat8 r1 = select(a, b, cond);
	EXPECT_EQ(r1.lane<0>(), 4.0f);
	EXPECT_EQ(r1.lane<1>(), 3.0f);
	EXPECT_EQ(r1.lane<2>(), 2.0f);
	EXPECT_EQ(r1.lane<3>(), 1.0f);
	EXPECT_EQ(r1.lane<4>(), 4.0f);
	EXPECT_EQ(r1.lane<5>(), 3.0f);
	EXPECT_EQ(r1.lane<6>(), 2.0f);
	EXPECT_EQ(r1.lane<7>(), 1.0f);

	// Select in the other
	vfloat8 r2 = select(b, a, cond);
	EXPECT_EQ(r2.lane<0>(), 1.0f);
	EXPECT_EQ(r2.lane<1>(), 2.0f);
	EXPECT_EQ(r2.lane<2>(), 3.0f);
	EXPECT_EQ(r2.lane<3>(), 4.0f);
	EXPECT_EQ(r2.lane<4>(), 1.0f);
	EXPECT_EQ(r2.lane<5>(), 2.0f);
	EXPECT_EQ(r2.lane<6>(), 3.0f);
	EXPECT_EQ(r2.lane<7>(), 4.0f);
}

/** @brief Test vfloat8 select MSB only. */
TEST(vfloat8, select_msb)
{
	vint8 msb(0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000, 0);
	vmask8 cond(msb.m);

	vfloat8 a(1.0f, 3.0f, 3.0f, 1.0f, 1.0f, 3.0f, 3.0f, 1.0f);
	vfloat8 b(4.0f, 2.0f, 2.0f, 4.0f, 4.0f, 2.0f, 2.0f, 4.0f);

	// Select in one direction
	vfloat8 r1 = select(a, b, cond);
	EXPECT_EQ(r1.lane<0>(), 4.0f);
	EXPECT_EQ(r1.lane<1>(), 3.0f);
	EXPECT_EQ(r1.lane<2>(), 2.0f);
	EXPECT_EQ(r1.lane<3>(), 1.0f);
	EXPECT_EQ(r1.lane<4>(), 4.0f);
	EXPECT_EQ(r1.lane<5>(), 3.0f);
	EXPECT_EQ(r1.lane<6>(), 2.0f);
	EXPECT_EQ(r1.lane<7>(), 1.0f);


	// Select in the other
	vfloat8 r2 = select(b, a, cond);
	EXPECT_EQ(r2.lane<0>(), 1.0f);
	EXPECT_EQ(r2.lane<1>(), 2.0f);
	EXPECT_EQ(r2.lane<2>(), 3.0f);
	EXPECT_EQ(r2.lane<3>(), 4.0f);
	EXPECT_EQ(r2.lane<4>(), 1.0f);
	EXPECT_EQ(r2.lane<5>(), 2.0f);
	EXPECT_EQ(r2.lane<6>(), 3.0f);
	EXPECT_EQ(r2.lane<7>(), 4.0f);
}

/** @brief Test vfloat8 gatherf. */
TEST(vfloat8, gatherf)
{
	vint8 indices(0, 4, 3, 2, 7, 4, 3, 2);
	vfloat8 r = gatherf(f32_data, indices);
	EXPECT_EQ(r.lane<0>(), 0.0f);
	EXPECT_EQ(r.lane<1>(), 4.0f);
	EXPECT_EQ(r.lane<2>(), 3.0f);
	EXPECT_EQ(r.lane<3>(), 2.0f);
	EXPECT_EQ(r.lane<4>(), 7.0f);
	EXPECT_EQ(r.lane<5>(), 4.0f);
	EXPECT_EQ(r.lane<6>(), 3.0f);
	EXPECT_EQ(r.lane<7>(), 2.0f);
}

/** @brief Test vfloat8 store. */
TEST(vfloat8, store)
{
	alignas(32) float out[9];
	vfloat8 a(f32_data);
	store(a, &(out[1]));
	EXPECT_EQ(out[1], 0.0f);
	EXPECT_EQ(out[2], 1.0f);
	EXPECT_EQ(out[3], 2.0f);
	EXPECT_EQ(out[4], 3.0f);
	EXPECT_EQ(out[5], 4.0f);
	EXPECT_EQ(out[6], 5.0f);
	EXPECT_EQ(out[7], 6.0f);
	EXPECT_EQ(out[8], 7.0f);
}

/** @brief Test vfloat8 storea. */
TEST(vfloat8, storea)
{
	alignas(32) float out[9];
	vfloat8 a(f32_data);
	store(a, out);
	EXPECT_EQ(out[0], 0.0f);
	EXPECT_EQ(out[1], 1.0f);
	EXPECT_EQ(out[2], 2.0f);
	EXPECT_EQ(out[3], 3.0f);
	EXPECT_EQ(out[4], 4.0f);
	EXPECT_EQ(out[5], 5.0f);
	EXPECT_EQ(out[6], 6.0f);
	EXPECT_EQ(out[7], 7.0f);
}

/** @brief Test vfloat8 float_to_int. */
TEST(vfloat8, float_to_int)
{
	vfloat8 a(1.1f, 1.5f, 1.6f, 4.0f, 1.1f, 1.5f, 1.6f, 4.0f);
	vint8 r = float_to_int(a);
	EXPECT_EQ(r.lane<0>(), 1);
	EXPECT_EQ(r.lane<1>(), 1);
	EXPECT_EQ(r.lane<2>(), 1);
	EXPECT_EQ(r.lane<3>(), 4);
	EXPECT_EQ(r.lane<4>(), 1);
	EXPECT_EQ(r.lane<5>(), 1);
	EXPECT_EQ(r.lane<6>(), 1);
	EXPECT_EQ(r.lane<7>(), 4);
}

// vint8 tests - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

/** @brief Test unaligned vint8 data load. */
TEST(vint8, UnalignedLoad)
{
	vint8 a(&(s32_data[1]));
	EXPECT_EQ(a.lane<0>(), 1);
	EXPECT_EQ(a.lane<1>(), 2);
	EXPECT_EQ(a.lane<2>(), 3);
	EXPECT_EQ(a.lane<3>(), 4);
	EXPECT_EQ(a.lane<4>(), 5);
	EXPECT_EQ(a.lane<5>(), 6);
	EXPECT_EQ(a.lane<6>(), 7);
	EXPECT_EQ(a.lane<7>(), 8);
}

/** @brief Test unaligned vint8 data load. */
TEST(vint8, UnalignedLoad8)
{
	vint8 a(&(u8_data[1]));
	EXPECT_EQ(a.lane<0>(), 1);
	EXPECT_EQ(a.lane<1>(), 2);
	EXPECT_EQ(a.lane<2>(), 3);
	EXPECT_EQ(a.lane<3>(), 4);
	EXPECT_EQ(a.lane<4>(), 5);
	EXPECT_EQ(a.lane<5>(), 6);
	EXPECT_EQ(a.lane<6>(), 7);
	EXPECT_EQ(a.lane<7>(), 8);
}


/** @brief Test scalar duplicated vint8 load. */
TEST(vint8, ScalarDupLoad)
{
	vint8 a(42);
	EXPECT_EQ(a.lane<0>(), 42);
	EXPECT_EQ(a.lane<1>(), 42);
	EXPECT_EQ(a.lane<2>(), 42);
	EXPECT_EQ(a.lane<3>(), 42);
	EXPECT_EQ(a.lane<4>(), 42);
	EXPECT_EQ(a.lane<5>(), 42);
	EXPECT_EQ(a.lane<6>(), 42);
	EXPECT_EQ(a.lane<7>(), 42);
}

/** @brief Test scalar vint8 load. */
TEST(vint8, ScalarLoad)
{
	vint8 a(11, 22, 33, 44, 55, 66, 77, 88);
	EXPECT_EQ(a.lane<0>(), 11);
	EXPECT_EQ(a.lane<1>(), 22);
	EXPECT_EQ(a.lane<2>(), 33);
	EXPECT_EQ(a.lane<3>(), 44);
	EXPECT_EQ(a.lane<4>(), 55);
	EXPECT_EQ(a.lane<5>(), 66);
	EXPECT_EQ(a.lane<6>(), 77);
	EXPECT_EQ(a.lane<7>(), 88);
}

/** @brief Test copy vint8 load. */
TEST(vint8, CopyLoad)
{
	vint8 s(11, 22, 33, 44, 55, 66, 77, 88);
	vint8 a(s.m);
	EXPECT_EQ(a.lane<0>(), 11);
	EXPECT_EQ(a.lane<1>(), 22);
	EXPECT_EQ(a.lane<2>(), 33);
	EXPECT_EQ(a.lane<3>(), 44);
	EXPECT_EQ(a.lane<4>(), 55);
	EXPECT_EQ(a.lane<5>(), 66);
	EXPECT_EQ(a.lane<6>(), 77);
	EXPECT_EQ(a.lane<7>(), 88);
}

/** @brief Test vint8 zero. */
TEST(vint8, Zero)
{
	vint8 a = vint8::zero();
	EXPECT_EQ(a.lane<0>(), 0);
	EXPECT_EQ(a.lane<1>(), 0);
	EXPECT_EQ(a.lane<2>(), 0);
	EXPECT_EQ(a.lane<3>(), 0);
	EXPECT_EQ(a.lane<4>(), 0);
	EXPECT_EQ(a.lane<5>(), 0);
	EXPECT_EQ(a.lane<6>(), 0);
	EXPECT_EQ(a.lane<7>(), 0);
}

/** @brief Test vint8 load1. */
TEST(vint8, Load1)
{
	int s = 42;
	vint8 a = vint8::load1(&s);
	EXPECT_EQ(a.lane<0>(), 42);
	EXPECT_EQ(a.lane<1>(), 42);
	EXPECT_EQ(a.lane<2>(), 42);
	EXPECT_EQ(a.lane<3>(), 42);
	EXPECT_EQ(a.lane<4>(), 42);
	EXPECT_EQ(a.lane<5>(), 42);
	EXPECT_EQ(a.lane<6>(), 42);
	EXPECT_EQ(a.lane<7>(), 42);
}

/** @brief Test vint8 loada. */
TEST(vint8, Loada)
{
	vint8 a(&(s32_data[0]));
	EXPECT_EQ(a.lane<0>(), 0);
	EXPECT_EQ(a.lane<1>(), 1);
	EXPECT_EQ(a.lane<2>(), 2);
	EXPECT_EQ(a.lane<3>(), 3);
	EXPECT_EQ(a.lane<4>(), 4);
	EXPECT_EQ(a.lane<5>(), 5);
	EXPECT_EQ(a.lane<6>(), 6);
	EXPECT_EQ(a.lane<7>(), 7);
}

/** @brief Test vint8 lane_id. */
TEST(vint8, LaneID)
{
	vint8 a = vint8::lane_id();
	EXPECT_EQ(a.lane<0>(), 0);
	EXPECT_EQ(a.lane<1>(), 1);
	EXPECT_EQ(a.lane<2>(), 2);
	EXPECT_EQ(a.lane<3>(), 3);
	EXPECT_EQ(a.lane<4>(), 4);
	EXPECT_EQ(a.lane<5>(), 5);
	EXPECT_EQ(a.lane<6>(), 6);
	EXPECT_EQ(a.lane<7>(), 7);
}

/** @brief Test vint8 add. */
TEST(vint8, vadd)
{
	vint8 a(1, 2, 3, 4, 1, 2, 3, 4);
	vint8 b(2, 3, 4, 5, 2, 3, 4, 5);
	a = a + b;
	EXPECT_EQ(a.lane<0>(), 1 + 2);
	EXPECT_EQ(a.lane<1>(), 2 + 3);
	EXPECT_EQ(a.lane<2>(), 3 + 4);
	EXPECT_EQ(a.lane<3>(), 4 + 5);
	EXPECT_EQ(a.lane<4>(), 1 + 2);
	EXPECT_EQ(a.lane<5>(), 2 + 3);
	EXPECT_EQ(a.lane<6>(), 3 + 4);
	EXPECT_EQ(a.lane<7>(), 4 + 5);
}

/** @brief Test vint8 sub. */
TEST(vint8, vsub)
{
	vint8 a(1, 2, 4, 4, 1, 2, 4, 4);
	vint8 b(2, 3, 3, 5, 2, 3, 3, 5);
	a = a - b;
	EXPECT_EQ(a.lane<0>(), 1 - 2);
	EXPECT_EQ(a.lane<1>(), 2 - 3);
	EXPECT_EQ(a.lane<2>(), 4 - 3);
	EXPECT_EQ(a.lane<3>(), 4 - 5);
	EXPECT_EQ(a.lane<4>(), 1 - 2);
	EXPECT_EQ(a.lane<5>(), 2 - 3);
	EXPECT_EQ(a.lane<6>(), 4 - 3);
	EXPECT_EQ(a.lane<7>(), 4 - 5);
}

/** @brief Test vint8 bitwise invert. */
TEST(vint8, bit_invert)
{
	vint8 a(-1, 0, 1, 2, -1, 0, 1, 2);
	a = ~a;
	EXPECT_EQ(a.lane<0>(), ~-1);
	EXPECT_EQ(a.lane<1>(), ~0);
	EXPECT_EQ(a.lane<2>(), ~1);
	EXPECT_EQ(a.lane<3>(), ~2);
	EXPECT_EQ(a.lane<4>(), ~-1);
	EXPECT_EQ(a.lane<5>(), ~0);
	EXPECT_EQ(a.lane<6>(), ~1);
	EXPECT_EQ(a.lane<7>(), ~2);
}

/** @brief Test vint8 bitwise or. */
TEST(vint8, bit_vor)
{
	vint8 a(1, 2, 3, 4, 1, 2, 3, 4);
	vint8 b(2, 3, 4, 5, 2, 3, 4, 5);
	a = a | b;
	EXPECT_EQ(a.lane<0>(), 3);
	EXPECT_EQ(a.lane<1>(), 3);
	EXPECT_EQ(a.lane<2>(), 7);
	EXPECT_EQ(a.lane<3>(), 5);
	EXPECT_EQ(a.lane<4>(), 3);
	EXPECT_EQ(a.lane<5>(), 3);
	EXPECT_EQ(a.lane<6>(), 7);
	EXPECT_EQ(a.lane<7>(), 5);
}

/** @brief Test vint8 bitwise and. */
TEST(vint8, bit_vand)
{
	vint8 a(1, 2, 3, 4, 1, 2, 3, 4);
	vint8 b(2, 3, 4, 5, 2, 3, 4, 5);
	a = a & b;
	EXPECT_EQ(a.lane<0>(), 0);
	EXPECT_EQ(a.lane<1>(), 2);
	EXPECT_EQ(a.lane<2>(), 0);
	EXPECT_EQ(a.lane<3>(), 4);
	EXPECT_EQ(a.lane<4>(), 0);
	EXPECT_EQ(a.lane<5>(), 2);
	EXPECT_EQ(a.lane<6>(), 0);
	EXPECT_EQ(a.lane<7>(), 4);
}

/** @brief Test vint8 bitwise xor. */
TEST(vint8, bit_vxor)
{
	vint8 a(1, 2, 3, 4, 1, 2, 3, 4);
	vint8 b(2, 3, 4, 5, 2, 3, 4, 5);
	a = a ^ b;
	EXPECT_EQ(a.lane<0>(), 3);
	EXPECT_EQ(a.lane<1>(), 1);
	EXPECT_EQ(a.lane<2>(), 7);
	EXPECT_EQ(a.lane<3>(), 1);
	EXPECT_EQ(a.lane<4>(), 3);
	EXPECT_EQ(a.lane<5>(), 1);
	EXPECT_EQ(a.lane<6>(), 7);
	EXPECT_EQ(a.lane<7>(), 1);
}

/** @brief Test vint8 ceq. */
TEST(vint8, ceq)
{
	vint8 a1(1, 2, 3, 4, 1, 2, 3, 4);
	vint8 b1(0, 1, 2, 3, 0, 1, 2, 3);
	vmask8 r1 = a1 == b1;
	EXPECT_EQ(0, mask(r1));
	EXPECT_EQ(false, any(r1));
	EXPECT_EQ(false, all(r1));

	vint8 a2(1, 2, 3, 4, 1, 2, 3, 4);
	vint8 b2(1, 0, 0, 0, 1, 0, 0, 0);
	vmask8 r2 = a2 == b2;
	EXPECT_EQ(0x11, mask(r2));
	EXPECT_EQ(true, any(r2));
	EXPECT_EQ(false, all(r2));

	vint8 a3(1, 2, 3, 4, 1, 2, 3, 4);
	vint8 b3(1, 0, 3, 0, 1, 0, 3, 0);
	vmask8 r3 = a3 == b3;
	EXPECT_EQ(0x55, mask(r3));
	EXPECT_EQ(true, any(r3));
	EXPECT_EQ(false, all(r3));

	vint8 a4(1, 2, 3, 4, 1, 2, 3, 4);
	vmask8 r4 = a4 == a4;
	EXPECT_EQ(0xFF, mask(r4));
	EXPECT_EQ(true, any(r4));
	EXPECT_EQ(true, all(r4));
}

/** @brief Test vint8 cne. */
TEST(vint8, cne)
{
	vint8 a1(1, 2, 3, 4, 1, 2, 3, 4);
	vint8 b1(0, 1, 2, 3, 0, 1, 2, 3);
	vmask8 r1 = a1 != b1;
	EXPECT_EQ(0xFF, mask(r1));
	EXPECT_EQ(true, any(r1));
	EXPECT_EQ(true, all(r1));

	vint8 a2(1, 2, 3, 4, 1, 2, 3, 4);
	vint8 b2(1, 0, 0, 0, 1, 0, 0, 0);
	vmask8 r2 = a2 != b2;
	EXPECT_EQ(0xEE, mask(r2));
	EXPECT_EQ(true, any(r2));
	EXPECT_EQ(false, all(r2));

	vint8 a3(1, 2, 3, 4, 1, 2, 3, 4);
	vint8 b3(1, 0, 3, 0, 1, 0, 3, 0);
	vmask8 r3 = a3 != b3;
	EXPECT_EQ(0xAA, mask(r3));
	EXPECT_EQ(true, any(r3));
	EXPECT_EQ(false, all(r3));

	vint8 a4(1, 2, 3, 4, 1, 2, 3, 4);
	vmask8 r4 = a4 != a4;
	EXPECT_EQ(0, mask(r4));
	EXPECT_EQ(false, any(r4));
	EXPECT_EQ(false, all(r4));
}

/** @brief Test vint8 clt. */
TEST(vint8, clt)
{
	vint8 a(1, 2, 3, 4, 1, 2, 3, 4);
	vint8 b(0, 3, 3, 5, 0, 3, 3, 5);
	vmask8 r = a < b;
	EXPECT_EQ(0xAA, mask(r));
}

/** @brief Test vint8 cgt. */
TEST(vint8, cle)
{
	vint8 a(1, 2, 3, 4, 1, 2, 3, 4);
	vint8 b(0, 3, 3, 5, 0, 3, 3, 5);
	vmask8 r = a > b;
	EXPECT_EQ(0x11, mask(r));
}

/** @brief Test vint8 min. */
TEST(vint8, min)
{
	vint8 a(1, 2, 3, 4, 1, 2, 3, 4);
	vint8 b(0, 3, 3, 5, 0, 3, 3, 5);
	vint8 r = min(a, b);
	EXPECT_EQ(r.lane<0>(), 0);
	EXPECT_EQ(r.lane<1>(), 2);
	EXPECT_EQ(r.lane<2>(), 3);
	EXPECT_EQ(r.lane<3>(), 4);
	EXPECT_EQ(r.lane<4>(), 0);
	EXPECT_EQ(r.lane<5>(), 2);
	EXPECT_EQ(r.lane<6>(), 3);
	EXPECT_EQ(r.lane<7>(), 4);
}

/** @brief Test vint8 max. */
TEST(vint8, max)
{
	vint8 a(1, 2, 3, 4, 1, 2, 3, 4);
	vint8 b(0, 3, 3, 5, 0, 3, 3, 5);
	vint8 r = max(a, b);
	EXPECT_EQ(r.lane<0>(), 1);
	EXPECT_EQ(r.lane<1>(), 3);
	EXPECT_EQ(r.lane<2>(), 3);
	EXPECT_EQ(r.lane<3>(), 5);
	EXPECT_EQ(r.lane<4>(), 1);
	EXPECT_EQ(r.lane<5>(), 3);
	EXPECT_EQ(r.lane<6>(), 3);
	EXPECT_EQ(r.lane<7>(), 5);
}

/** @brief Test vint8 hmin. */
TEST(vint8, hmin)
{
	vint8 a1(1, 2, 1, 2, 1, 2, 1, 2);
	vint8 r1 = hmin(a1);
	EXPECT_EQ(r1.lane<0>(), 1);
	EXPECT_EQ(r1.lane<1>(), 1);
	EXPECT_EQ(r1.lane<2>(), 1);
	EXPECT_EQ(r1.lane<3>(), 1);
	EXPECT_EQ(r1.lane<4>(), 1);
	EXPECT_EQ(r1.lane<5>(), 1);
	EXPECT_EQ(r1.lane<6>(), 1);
	EXPECT_EQ(r1.lane<7>(), 1);

	vint8 a2(1, 2, -1, 5, 1, 2, -1, 5);
	vint8 r2 = hmin(a2);
	EXPECT_EQ(r2.lane<0>(), -1);
	EXPECT_EQ(r2.lane<1>(), -1);
	EXPECT_EQ(r2.lane<2>(), -1);
	EXPECT_EQ(r2.lane<3>(), -1);
	EXPECT_EQ(r2.lane<4>(), -1);
	EXPECT_EQ(r2.lane<5>(), -1);
	EXPECT_EQ(r2.lane<6>(), -1);
	EXPECT_EQ(r2.lane<7>(), -1);
}

/** @brief Test vint8 hmax. */
TEST(vint8, hmax)
{
	vint8 a1(1, 2, 1, 2, 1, 3, 1, 2);
	vint8 r1 = hmax(a1);
	EXPECT_EQ(r1.lane<0>(), 3);
	EXPECT_EQ(r1.lane<1>(), 3);
	EXPECT_EQ(r1.lane<2>(), 3);
	EXPECT_EQ(r1.lane<3>(), 3);
	EXPECT_EQ(r1.lane<4>(), 3);
	EXPECT_EQ(r1.lane<5>(), 3);
	EXPECT_EQ(r1.lane<6>(), 3);
	EXPECT_EQ(r1.lane<7>(), 3);

	vint8 a2(1, 2, -1, 5, 1, 2, -1, 5);
	vint8 r2 = hmax(a2);
	EXPECT_EQ(r2.lane<0>(), 5);
	EXPECT_EQ(r2.lane<1>(), 5);
	EXPECT_EQ(r2.lane<2>(), 5);
	EXPECT_EQ(r2.lane<3>(), 5);
	EXPECT_EQ(r2.lane<4>(), 5);
	EXPECT_EQ(r2.lane<5>(), 5);
	EXPECT_EQ(r2.lane<6>(), 5);
	EXPECT_EQ(r2.lane<7>(), 5);
}

/** @brief Test vint8 storea. */
TEST(vint8, storea)
{
	alignas(32) int out[8];
	vint8 a(s32_data);
	storea(a, out);
	EXPECT_EQ(out[0], 0);
	EXPECT_EQ(out[1], 1);
	EXPECT_EQ(out[2], 2);
	EXPECT_EQ(out[3], 3);
	EXPECT_EQ(out[4], 4);
	EXPECT_EQ(out[5], 5);
	EXPECT_EQ(out[6], 6);
	EXPECT_EQ(out[7], 7);
}

/** @brief Test vint8 store_nbytes. */
TEST(vint8, store_nbytes)
{
	alignas(32) int out[2];
	vint8 a(42, 314, 75, 90, 42, 314, 75, 90);
	store_nbytes(a, (uint8_t*)&out);
	EXPECT_EQ(out[0], 42);
	EXPECT_EQ(out[1], 314);
}

/** @brief Test vint8 gatheri. */
TEST(vint8, gatheri)
{
	vint8 indices(0, 4, 3, 2, 7, 4, 3, 2);
	vint8 r = gatheri(s32_data, indices);
	EXPECT_EQ(r.lane<0>(), 0);
	EXPECT_EQ(r.lane<1>(), 4);
	EXPECT_EQ(r.lane<2>(), 3);
	EXPECT_EQ(r.lane<3>(), 2);
	EXPECT_EQ(r.lane<4>(), 7);
	EXPECT_EQ(r.lane<5>(), 4);
	EXPECT_EQ(r.lane<6>(), 3);
	EXPECT_EQ(r.lane<7>(), 2);
}

/** @brief Test vint8 pack_low_bytes. */
TEST(vint8, pack_low_bytes)
{
	vint8 a(1, 2, 3, 4, 2, 3, 4, 5);
	vint8 r = pack_low_bytes(a);
	EXPECT_EQ(r.lane<0>(), (4 << 24) | (3 << 16) | (2  << 8) | (1 << 0));
	EXPECT_EQ(r.lane<1>(), (5 << 24) | (4 << 16) | (3  << 8) | (2 << 0));
}

/** @brief Test vint8 select. */
TEST(vint8, select)
{
	vint8 m1(1, 1, 1, 1, 1, 1, 1, 1);
	vint8 m2(1, 2, 1, 2, 1, 2, 1, 2);
	vmask8 cond = m1 == m2;

	vint8 a(1, 3, 3, 1, 1, 3, 3, 1);
	vint8 b(4, 2, 2, 4, 4, 2, 2, 4);

	vint8 r1 = select(a, b, cond);
	EXPECT_EQ(r1.lane<0>(), 4);
	EXPECT_EQ(r1.lane<1>(), 3);
	EXPECT_EQ(r1.lane<2>(), 2);
	EXPECT_EQ(r1.lane<3>(), 1);
	EXPECT_EQ(r1.lane<4>(), 4);
	EXPECT_EQ(r1.lane<5>(), 3);
	EXPECT_EQ(r1.lane<6>(), 2);
	EXPECT_EQ(r1.lane<7>(), 1);

	vint8 r2 = select(b, a, cond);
	EXPECT_EQ(r2.lane<0>(), 1);
	EXPECT_EQ(r2.lane<1>(), 2);
	EXPECT_EQ(r2.lane<2>(), 3);
	EXPECT_EQ(r2.lane<3>(), 4);
	EXPECT_EQ(r2.lane<4>(), 1);
	EXPECT_EQ(r2.lane<5>(), 2);
	EXPECT_EQ(r2.lane<6>(), 3);
	EXPECT_EQ(r2.lane<7>(), 4);
}

/** @brief Test vint8 select MSB. */
TEST(vint8, select_msb)
{
	vint8 msb(0x80000000, 0, 0x80000000, 0, 0x80000000, 0, 0x80000000, 0);
	vmask8 cond(msb.m);

	vint8 a(1, 3, 3, 1, 1, 3, 3, 1);
	vint8 b(4, 2, 2, 4, 4, 2, 2, 4);

	vint8 r1 = select(a, b, cond);
	EXPECT_EQ(r1.lane<0>(), 4);
	EXPECT_EQ(r1.lane<1>(), 3);
	EXPECT_EQ(r1.lane<2>(), 2);
	EXPECT_EQ(r1.lane<3>(), 1);
	EXPECT_EQ(r1.lane<4>(), 4);
	EXPECT_EQ(r1.lane<5>(), 3);
	EXPECT_EQ(r1.lane<6>(), 2);
	EXPECT_EQ(r1.lane<7>(), 1);

	vint8 r2 = select(b, a, cond);
	EXPECT_EQ(r2.lane<0>(), 1);
	EXPECT_EQ(r2.lane<1>(), 2);
	EXPECT_EQ(r2.lane<2>(), 3);
	EXPECT_EQ(r2.lane<3>(), 4);
	EXPECT_EQ(r2.lane<4>(), 1);
	EXPECT_EQ(r2.lane<5>(), 2);
	EXPECT_EQ(r2.lane<6>(), 3);
	EXPECT_EQ(r2.lane<7>(), 4);
}

// vmask8 tests - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
/** @brief Test vmask8 or. */
TEST(vmask8, or)
{
	vfloat8 m1a(0, 1, 0, 1, 0, 1, 0, 1);
	vfloat8 m1b(1, 1, 1, 1, 1, 1, 1, 1);
	vmask8 m1 = m1a == m1b;

	vfloat8 m2a(1, 1, 0, 0, 1, 1, 0, 0);
	vfloat8 m2b(1, 1, 1, 1, 1, 1, 1, 1);
	vmask8 m2 = m2a == m2b;

	vmask8 r = m1 | m2;
	EXPECT_EQ(mask(r), 0xBB);
}

/** @brief Test vmask8 and. */
TEST(vmask8, and)
{
	vfloat8 m1a(0, 1, 0, 1, 0, 1, 0, 1);
	vfloat8 m1b(1, 1, 1, 1, 1, 1, 1, 1);
	vmask8 m1 = m1a == m1b;

	vfloat8 m2a(1, 1, 0, 0, 1, 1, 0, 0);
	vfloat8 m2b(1, 1, 1, 1, 1, 1, 1, 1);
	vmask8 m2 = m2a == m2b;

	vmask8 r = m1 & m2;
	EXPECT_EQ(mask(r), 0x22);
}

/** @brief Test vmask8 xor. */
TEST(vmask8, xor)
{
	vfloat8 m1a(0, 1, 0, 1, 0, 1, 0, 1);
	vfloat8 m1b(1, 1, 1, 1, 1, 1, 1, 1);
	vmask8 m1 = m1a == m1b;

	vfloat8 m2a(1, 1, 0, 0, 1, 1, 0, 0);
	vfloat8 m2b(1, 1, 1, 1, 1, 1, 1, 1);
	vmask8 m2 = m2a == m2b;

	vmask8 r = m1 ^ m2;
	EXPECT_EQ(mask(r), 0x99);
}

/** @brief Test vmask8 not. */
TEST(vmask8, not)
{
	vfloat8 m1a(0, 1, 0, 1, 0, 1, 0, 1);
	vfloat8 m1b(1, 1, 1, 1, 1, 1, 1, 1);
	vmask8 m1 = m1a == m1b;
	vmask8 r = ~m1;
	EXPECT_EQ(mask(r), 0x55);
}

#endif

}
