// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2020 Arm Limited
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

#include <limits>

#include "gtest/gtest.h"

#include "../astcenc_internal.h"
#include "../astcenc_vecmathlib.h"

static const float qnan = std::numeric_limits<float>::quiet_NaN();
alignas(16) static const float f32x4_data[5] { 0.0f, 1.0f, 2.0f, 3.0f, 4.0f };
alignas(16) static const int s32x4_data[5] {      0,    1,    2,    3,    4 };

namespace astcenc
{

#if ASTCENC_SIMD_WIDTH >= 4

// VFLOAT4 tests - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

/** \brief Test unaligned vfloat4 data load. */
TEST(vfloat4, UnalignedLoad)
{
	vfloat4 a(&(f32x4_data[1]));
	EXPECT_EQ(a.lane<0>(), 1.0f);
	EXPECT_EQ(a.lane<1>(), 2.0f);
	EXPECT_EQ(a.lane<2>(), 3.0f);
	EXPECT_EQ(a.lane<3>(), 4.0f);
}

/** \brief Test scalar duplicated vfloat4 load. */
TEST(vfloat4, ScalarDupLoad)
{
	vfloat4 a(1.1f);
	EXPECT_EQ(a.lane<0>(), 1.1f);
	EXPECT_EQ(a.lane<1>(), 1.1f);
	EXPECT_EQ(a.lane<2>(), 1.1f);
	EXPECT_EQ(a.lane<3>(), 1.1f);
}

/** \brief Test scalar vfloat4 load. */
TEST(vfloat4, ScalarLoad)
{
	vfloat4 a(1.1f, 2.2f, 3.3f, 4.4f);
	EXPECT_EQ(a.lane<0>(), 1.1f);
	EXPECT_EQ(a.lane<1>(), 2.2f);
	EXPECT_EQ(a.lane<2>(), 3.3f);
	EXPECT_EQ(a.lane<3>(), 4.4f);
}

/** \brief Test copy vfloat4 load. */
TEST(vfloat4, CopyLoad)
{
	vfloat4 s(1.1f, 2.2f, 3.3f, 4.4f);
	vfloat4 a(s.m);
	EXPECT_EQ(a.lane<0>(), 1.1f);
	EXPECT_EQ(a.lane<1>(), 2.2f);
	EXPECT_EQ(a.lane<2>(), 3.3f);
	EXPECT_EQ(a.lane<3>(), 4.4f);
}

/** \brief Test vfloat4 scalar lane set. */
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

/** \brief Test vfloat4 zero. */
TEST(vfloat4, Zero)
{
	vfloat4 a = vfloat4::zero();
	EXPECT_EQ(a.lane<0>(), 0.0f);
	EXPECT_EQ(a.lane<1>(), 0.0f);
	EXPECT_EQ(a.lane<2>(), 0.0f);
	EXPECT_EQ(a.lane<3>(), 0.0f);
}

/** \brief Test vfloat4 load1. */
TEST(vfloat4, Load1)
{
	float s = 3.14f;
	vfloat4 a = vfloat4::load1(&s);
	EXPECT_EQ(a.lane<0>(), 3.14f);
	EXPECT_EQ(a.lane<1>(), 3.14f);
	EXPECT_EQ(a.lane<2>(), 3.14f);
	EXPECT_EQ(a.lane<3>(), 3.14f);
}

/** \brief Test vfloat4 loada. */
TEST(vfloat4, Loada)
{
	vfloat4 a(&(f32x4_data[0]));
	EXPECT_EQ(a.lane<0>(), 0.0f);
	EXPECT_EQ(a.lane<1>(), 1.0f);
	EXPECT_EQ(a.lane<2>(), 2.0f);
	EXPECT_EQ(a.lane<3>(), 3.0f);
}

/** \brief Test vfloat4 lane_id. */
TEST(vfloat4, LaneID)
{
	vfloat4 a = vfloat4::lane_id();
	EXPECT_EQ(a.lane<0>(), 0.0f);
	EXPECT_EQ(a.lane<1>(), 1.0f);
	EXPECT_EQ(a.lane<2>(), 2.0f);
	EXPECT_EQ(a.lane<3>(), 3.0f);
}

/** \brief Test vfloat4 add. */
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

/** \brief Test vfloat4 sub. */
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

/** \brief Test vfloat4 mul. */
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

/** \brief Test vfloat4 mul. */
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

/** \brief Test vfloat4 mul. */
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

/** \brief Test vfloat4 div. */
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

/** \brief Test vfloat4 ceq. */
TEST(vfloat4, ceq)
{
	vfloat4 a1(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b1(0.1f, 0.2f, 0.3f, 0.4f);
	vmask r1 = a1 == b1;
	EXPECT_EQ(0, mask(r1));
	EXPECT_EQ(false, any(r1));
	EXPECT_EQ(false, all(r1));

	vfloat4 a2(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b2(1.0f, 0.2f, 0.3f, 0.4f);
	vmask r2 = a2 == b2;
	EXPECT_EQ(0x1, mask(r2));
	EXPECT_EQ(true, any(r2));
	EXPECT_EQ(false, all(r2));

	vfloat4 a3(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b3(1.0f, 0.2f, 3.0f, 0.4f);
	vmask r3 = a3 == b3;
	EXPECT_EQ(0x5, mask(r3));
	EXPECT_EQ(true, any(r3));
	EXPECT_EQ(false, all(r3));

	vfloat4 a4(1.0f, 2.0f, 3.0f, 4.0f);
	vmask r4 = a4 == a4;
	EXPECT_EQ(0xF, mask(r4));
	EXPECT_EQ(true, any(r4));
	EXPECT_EQ(true, all(r4));
}

/** \brief Test vfloat4 cne. */
TEST(vfloat4, cne)
{
	vfloat4 a1(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b1(0.1f, 0.2f, 0.3f, 0.4f);
	vmask r1 = a1 != b1;
	EXPECT_EQ(0xF, mask(r1));
	EXPECT_EQ(true, any(r1));
	EXPECT_EQ(true, all(r1));

	vfloat4 a2(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b2(1.0f, 0.2f, 0.3f, 0.4f);
	vmask r2 = a2 != b2;
	EXPECT_EQ(0xE, mask(r2));
	EXPECT_EQ(true, any(r2));
	EXPECT_EQ(false, all(r2));

	vfloat4 a3(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b3(1.0f, 0.2f, 3.0f, 0.4f);
	vmask r3 = a3 != b3;
	EXPECT_EQ(0xA, mask(r3));
	EXPECT_EQ(true, any(r3));
	EXPECT_EQ(false, all(r3));

	vfloat4 a4(1.0f, 2.0f, 3.0f, 4.0f);
	vmask r4 = a4 != a4;
	EXPECT_EQ(0, mask(r4));
	EXPECT_EQ(false, any(r4));
	EXPECT_EQ(false, all(r4));
}

/** \brief Test vfloat4 clt. */
TEST(vfloat4, clt)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b(0.9f, 2.1f, 3.0f, 4.1f);
	vmask r = a < b;
	EXPECT_EQ(0xA, mask(r));
}

/** \brief Test vfloat4 cle. */
TEST(vfloat4, cle)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b(0.9f, 2.1f, 3.0f, 4.1f);
	vmask r = a <= b;
	EXPECT_EQ(0xE, mask(r));
}

/** \brief Test vfloat4 cgt. */
TEST(vfloat4, cgt)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b(0.9f, 2.1f, 3.0f, 4.1f);
	vmask r = a > b;
	EXPECT_EQ(0x1, mask(r));
}

/** \brief Test vfloat4 cge. */
TEST(vfloat4, cge)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b(0.9f, 2.1f, 3.0f, 4.1f);
	vmask r = a >= b;
	EXPECT_EQ(0x5, mask(r));
}

/** \brief Test vfloat4 min. */
TEST(vfloat4, min)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b(0.9f, 2.1f, 3.0f, 4.1f);
	vfloat4 r = min(a, b);
	EXPECT_EQ(r.lane<0>(), 0.9f);
	EXPECT_EQ(r.lane<1>(), 2.0f);
	EXPECT_EQ(r.lane<2>(), 3.0f);
	EXPECT_EQ(r.lane<3>(), 4.0f);
}

/** \brief Test vfloat4 max. */
TEST(vfloat4, max)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 b(0.9f, 2.1f, 3.0f, 4.1f);
	vfloat4 r = max(a, b);
	EXPECT_EQ(r.lane<0>(), 1.0f);
	EXPECT_EQ(r.lane<1>(), 2.1f);
	EXPECT_EQ(r.lane<2>(), 3.0f);
	EXPECT_EQ(r.lane<3>(), 4.1f);
}

/** \brief Test vfloat4 clamp. */
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

/** \brief Test vfloat4 clampz. */
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

/** \brief Test vfloat4 clampz. */
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

/** \brief Test vfloat4 abs. */
TEST(vfloat4, abs)
{
	vfloat4 a(-1.0f, 0.0f, 0.1f, 4.0f);
	vfloat4 r = abs(a);
	EXPECT_EQ(r.lane<0>(), 1.0f);
	EXPECT_EQ(r.lane<1>(), 0.0f);
	EXPECT_EQ(r.lane<2>(), 0.1f);
	EXPECT_EQ(r.lane<3>(), 4.0f);
}

/** \brief Test vfloat4 round. */
TEST(vfloat4, round)
{
	vfloat4 a(1.1f, 1.5f, 1.6f, 4.0f);
	vfloat4 r = round(a);
	EXPECT_EQ(r.lane<0>(), 1.0f);
	EXPECT_EQ(r.lane<1>(), 2.0f);
	EXPECT_EQ(r.lane<2>(), 2.0f);
	EXPECT_EQ(r.lane<3>(), 4.0f);
}

/** \brief Test vfloat4 hmin. */
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

/** \brief Test vfloat4 sqrt. */
TEST(vfloat4, sqrt)
{
	vfloat4 a(1.0f, 2.0f, 3.0f, 4.0f);
	vfloat4 r = sqrt(a);
	EXPECT_EQ(r.lane<0>(), std::sqrt(1.0f));
	EXPECT_EQ(r.lane<1>(), std::sqrt(2.0f));
	EXPECT_EQ(r.lane<2>(), std::sqrt(3.0f));
	EXPECT_EQ(r.lane<3>(), std::sqrt(4.0f));
}

/** \brief Test vfloat4 select. */
TEST(vfloat4, select)
{
	vfloat4 m1(1.0f, 1.0f, 1.0f, 1.0f);
	vfloat4 m2(1.0f, 2.0f, 1.0f, 2.0f);
	vmask4 cond = m1 == m2;

	vfloat4 a(1.0f, 3.0f, 3.0f, 1.0f);
	vfloat4 b(4.0f, 2.0f, 2.0f, 4.0f);

	vfloat4 r1 = select(a, b, cond);
	EXPECT_EQ(r1.lane<0>(), 4.0f);
	EXPECT_EQ(r1.lane<1>(), 3.0f);
	EXPECT_EQ(r1.lane<2>(), 2.0f);
	EXPECT_EQ(r1.lane<3>(), 1.0f);

	vfloat4 r2 = select(b, a, cond);
	EXPECT_EQ(r2.lane<0>(), 1.0f);
	EXPECT_EQ(r2.lane<1>(), 2.0f);
	EXPECT_EQ(r2.lane<2>(), 3.0f);
	EXPECT_EQ(r2.lane<3>(), 4.0f);
}

/** \brief Test vfloat4 gatherf. */
TEST(vfloat4, gatherf)
{
	vint4 indices(0, 4, 3, 2);
	vfloat4 r = gatherf(f32x4_data, indices);
	EXPECT_EQ(r.lane<0>(), 0.0f);
	EXPECT_EQ(r.lane<1>(), 4.0f);
	EXPECT_EQ(r.lane<2>(), 3.0f);
	EXPECT_EQ(r.lane<3>(), 2.0f);
}

/** \brief Test vfloat4 store. */
TEST(vfloat4, store)
{
	alignas(16) float out[5];
	vfloat4 a(f32x4_data);
	store(a, &(out[1]));

	EXPECT_EQ(out[1], 0.0f);
	EXPECT_EQ(out[2], 1.0f);
	EXPECT_EQ(out[3], 2.0f);
	EXPECT_EQ(out[4], 3.0f);
}

/** \brief Test vfloat4 storea. */
TEST(vfloat4, storea)
{
	alignas(16) float out[4];
	vfloat4 a(f32x4_data);
	store(a, out);

	EXPECT_EQ(out[0], 0.0f);
	EXPECT_EQ(out[1], 1.0f);
	EXPECT_EQ(out[2], 2.0f);
	EXPECT_EQ(out[3], 3.0f);
}

/** \brief Test vfloat4 dot. */
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

/** \brief Test vfloat4 float_to_int. */
TEST(vfloat4, float_to_int)
{
	vfloat4 a(1.1f, 1.5f, 1.6f, 4.0f);
	vint4 r = float_to_int(a);
	EXPECT_EQ(r.lane<0>(), 1);
	EXPECT_EQ(r.lane<1>(), 1);
	EXPECT_EQ(r.lane<2>(), 1);
	EXPECT_EQ(r.lane<3>(), 4);
}

/** \brief Test vfloat4 round. */
TEST(vfloat4, float_to_int_rtn)
{
	vfloat4 a(1.1f, 1.5f, 1.6f, 4.0f);
	vint4 r = float_to_int_rtn(a);
	EXPECT_EQ(r.lane<0>(), 1);
	EXPECT_EQ(r.lane<1>(), 2);
	EXPECT_EQ(r.lane<2>(), 2);
	EXPECT_EQ(r.lane<3>(), 4);
}


// VINT4 tests - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

/** \brief Test unaligned vint4 data load. */
TEST(vint4, UnalignedLoad)
{
	vint4 a(&(s32x4_data[1]));
	EXPECT_EQ(a.lane<0>(), 1);
	EXPECT_EQ(a.lane<1>(), 2);
	EXPECT_EQ(a.lane<2>(), 3);
	EXPECT_EQ(a.lane<3>(), 4);
}

/** \brief Test scalar duplicated vint4 load. */
TEST(vint4, ScalarDupLoad)
{
	vint4 a(42);
	EXPECT_EQ(a.lane<0>(), 42);
	EXPECT_EQ(a.lane<1>(), 42);
	EXPECT_EQ(a.lane<2>(), 42);
	EXPECT_EQ(a.lane<3>(), 42);
}

/** \brief Test scalar vint4 load. */
TEST(vint4, ScalarLoad)
{
	vint4 a(11, 22, 33, 44);
	EXPECT_EQ(a.lane<0>(), 11);
	EXPECT_EQ(a.lane<1>(), 22);
	EXPECT_EQ(a.lane<2>(), 33);
	EXPECT_EQ(a.lane<3>(), 44);
}

/** \brief Test copy vint4 load. */
TEST(vint4, CopyLoad)
{
	vint4 s(11, 22, 33, 44);
	vint4 a(s.m);
	EXPECT_EQ(a.lane<0>(), 11);
	EXPECT_EQ(a.lane<1>(), 22);
	EXPECT_EQ(a.lane<2>(), 33);
	EXPECT_EQ(a.lane<3>(), 44);
}

/** \brief Test vint4 lane_id. */
TEST(vint4, LaneID)
{
	vint4 a = vint4::lane_id();
	EXPECT_EQ(a.lane<0>(), 0);
	EXPECT_EQ(a.lane<1>(), 1);
	EXPECT_EQ(a.lane<2>(), 2);
	EXPECT_EQ(a.lane<3>(), 3);
}

#endif

}
