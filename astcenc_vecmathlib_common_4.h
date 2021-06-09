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
 * @brief Generic 4x32-bit vector functions.
 *
 * This module implements generic 4-wide vector functions that are valid for
 * all instruction sets, typically implemented using lower level 4-wide
 * operations that are ISA-specific.
 */

#ifndef ASTC_VECMATHLIB_COMMON_4_H_INCLUDED
#define ASTC_VECMATHLIB_COMMON_4_H_INCLUDED

#ifndef ASTCENC_SIMD_INLINE
	#error "Include astcenc_vecmathlib.h, do not include directly"
#endif

#include <cstdio>

// ============================================================================
// vmask4 operators and functions
// ============================================================================

/**
 * @brief True if any lanes are enabled, false otherwise.
 */
ASTCENC_SIMD_INLINE bool any(vmask4 a)
{
	return mask(a) != 0;
}

/**
 * @brief True if all lanes are enabled, false otherwise.
 */
ASTCENC_SIMD_INLINE bool all(vmask4 a)
{
	return mask(a) == 0xF;
}

// ============================================================================
// vint4 operators and functions
// ============================================================================

/**
 * @brief Overload: vector by scalar addition.
 */
ASTCENC_SIMD_INLINE vint4 operator+(vint4 a, int b)
{
	return a + vint4(b);
}

/**
 * @brief Overload: vector by vector incremental addition.
 */
ASTCENC_SIMD_INLINE vint4& operator+=(vint4& a, const vint4& b)
{
	a = a + b;
	return a;
}

/**
 * @brief Overload: vector by scalar subtraction.
 */
ASTCENC_SIMD_INLINE vint4 operator-(vint4 a, int b)
{
	return a - vint4(b);
}

/**
 * @brief Overload: vector by scalar multiplication.
 */
ASTCENC_SIMD_INLINE vint4 operator*(vint4 a, int b)
{
	return a * vint4(b);
}

/**
 * @brief Overload: vector by scalar bitwise or.
 */
ASTCENC_SIMD_INLINE vint4 operator|(vint4 a, int b)
{
	return a | vint4(b);
}

/**
 * @brief Overload: vector by scalar bitwise and.
 */
ASTCENC_SIMD_INLINE vint4 operator&(vint4 a, int b)
{
	return a & vint4(b);
}

/**
 * @brief Overload: vector by scalar bitwise xor.
 */
ASTCENC_SIMD_INLINE vint4 operator^(vint4 a, int b)
{
	return a ^ vint4(b);
}

/**
 * @brief Return the clamped value between min and max.
 */
ASTCENC_SIMD_INLINE vint4 clamp(int minv, int maxv, vint4 a)
{
	return min(max(a, vint4(minv)), vint4(maxv));
}

/**
 * @brief Return the horizontal sum of RGB vector lanes as a scalar.
 */
ASTCENC_SIMD_INLINE int hadd_rgb_s(vint4 a)
{
	return a.lane<0>() + a.lane<1>() + a.lane<2>();
}

/**
 * @brief Debug function to print a vector of ints.
 */
ASTCENC_SIMD_INLINE void print(vint4 a)
{
	alignas(16) int v[4];
	storea(a, v);
	printf("v4_i32:\n  %8d %8d %8d %8d\n",
	       v[0], v[1], v[2], v[3]);
}

// ============================================================================
// vfloat4 operators and functions
// ============================================================================

/**
 * @brief Overload: vector by vector incremental addition.
 */
ASTCENC_SIMD_INLINE vfloat4& operator+=(vfloat4& a, const vfloat4& b)
{
	a = a + b;
	return a;
}

/**
 * @brief Overload: vector by scalar addition.
 */
ASTCENC_SIMD_INLINE vfloat4 operator+(vfloat4 a, float b)
{
	return a + vfloat4(b);
}

/**
 * @brief Overload: vector by scalar subtraction.
 */
ASTCENC_SIMD_INLINE vfloat4 operator-(vfloat4 a, float b)
{
	return a - vfloat4(b);
}

/**
 * @brief Overload: vector by scalar multiplication.
 */
ASTCENC_SIMD_INLINE vfloat4 operator*(vfloat4 a, float b)
{
	return a * vfloat4(b);
}

/**
 * @brief Overload: scalar by vector multiplication.
 */
ASTCENC_SIMD_INLINE vfloat4 operator*(float a, vfloat4 b)
{
	return vfloat4(a) * b;
}

/**
 * @brief Overload: vector by scalar division.
 */
ASTCENC_SIMD_INLINE vfloat4 operator/(vfloat4 a, float b)
{
	return a / vfloat4(b);
}

/**
 * @brief Overload: scalar by vector division.
 */
ASTCENC_SIMD_INLINE vfloat4 operator/(float a, vfloat4 b)
{
	return vfloat4(a) / b;
}

/**
 * @brief Return the min vector of a vector and a scalar.
 *
 * If either lane value is NaN, @c b will be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat4 min(vfloat4 a, float b)
{
	return min(a, vfloat4(b));
}

/**
 * @brief Return the max vector of a vector and a scalar.
 *
 * If either lane value is NaN, @c b will be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat4 max(vfloat4 a, float b)
{
	return max(a, vfloat4(b));
}

/**
 * @brief Return the clamped value between min and max.
 *
 * It is assumed that neither @c min nor @c max are NaN values. If @c a is NaN
 * then @c min will be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat4 clamp(float minv, float maxv, vfloat4 a)
{
	// Do not reorder - second operand will return if either is NaN
	return min(max(a, minv), maxv);
}

/**
 * @brief Return the clamped value between 0.0f and max.
 *
 * It is assumed that  @c max is not a NaN value. If @c a is NaN then zero will
 * be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat4 clampz(float maxv, vfloat4 a)
{
	// Do not reorder - second operand will return if either is NaN
	return min(max(a, vfloat4::zero()), maxv);
}

/**
 * @brief Return the clamped value between 0.0f and 1.0f.
 *
 * If @c a is NaN then zero will be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat4 clampzo(vfloat4 a)
{
	// Do not reorder - second operand will return if either is NaN
	return min(max(a, vfloat4::zero()), 1.0f);
}

/**
 * @brief Return the horizontal minimum of a vector.
 */
ASTCENC_SIMD_INLINE float hmin_s(vfloat4 a)
{
	return hmin(a).lane<0>();
}

/**
 * @brief Return the horizontal min of RGB vector lanes as a scalar.
 */
ASTCENC_SIMD_INLINE float hmin_rgb_s(vfloat4 a)
{
	a.set_lane<3>(a.lane<0>());
	return hmin_s(a);
}

/**
 * @brief Return the horizontal maximum of a vector.
 */
ASTCENC_SIMD_INLINE float hmax_s(vfloat4 a)
{
	return hmax(a).lane<0>();
}

/**
 * @brief Accumulate the full horizontal sum of a vector.
 */
ASTCENC_SIMD_INLINE void haccumulate(float& accum, vfloat4 a)
{
	accum += hadd_s(a);
}

/**
 * @brief Accumulate lane-wise sums for a vector.
 */
ASTCENC_SIMD_INLINE void haccumulate(vfloat4& accum, vfloat4 a)
{
	accum = accum + a;
}

/**
 * @brief Return the horizontal sum of RGB vector lanes as a scalar.
 */
ASTCENC_SIMD_INLINE float hadd_rgb_s(vfloat4 a)
{
	return a.lane<0>() + a.lane<1>() + a.lane<2>();
}

/**
 * @brief Return the dot product for the full 4 lanes, returning scalar.
 */
ASTCENC_SIMD_INLINE float dot_s(vfloat4 a, vfloat4 b)
{
	vfloat4 m = a * b;
	return hadd_s(m);
}

/**
 * @brief Return the dot product for the full 4 lanes, returning vector.
 */
ASTCENC_SIMD_INLINE vfloat4 dot(vfloat4 a, vfloat4 b)
{
	vfloat4 m = a * b;
	return vfloat4(hadd_s(m));
}

/**
 * @brief Return the dot product for the bottom 3 lanes, returning scalar.
 */
ASTCENC_SIMD_INLINE float dot3_s(vfloat4 a, vfloat4 b)
{
	vfloat4 m = a * b;
	return hadd_rgb_s(m);
}

/**
 * @brief Return the dot product for the full 4 lanes, returning vector.
 */
ASTCENC_SIMD_INLINE vfloat4 dot3(vfloat4 a, vfloat4 b)
{
	vfloat4 m = a * b;
	float d3 = hadd_rgb_s(m);
	return vfloat4(d3, d3, d3, 0.0f);
}

/**
 * @brief Generate a reciprocal of a vector.
 */
ASTCENC_SIMD_INLINE vfloat4 recip(vfloat4 b)
{
	return 1.0f / b;
}

/**
 * @brief Debug function to print a vector of floats.
 */
ASTCENC_SIMD_INLINE void print(vfloat4 a)
{
	alignas(16) float v[4];
	storea(a, v);
	printf("v4_f32:\n  %0.4f %0.4f %0.4f %0.4f\n",
	       (double)v[0], (double)v[1], (double)v[2], (double)v[3]);
}

#endif // #ifndef ASTC_VECMATHLIB_COMMON_4_H_INCLUDED
