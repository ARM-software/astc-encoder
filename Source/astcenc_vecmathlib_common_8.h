// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2019-2026 Arm Limited
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
 * @brief Generic 8x32-bit vector functions.
 *
 * This module implements generic 8-wide vector functions that are valid for
 * all instruction sets, typically implemented using lower level 8-wide
 * operations that are ISA-specific.
 */

#if !defined(ASTC_VECMATHLIB_COMMON_8_H_INCLUDED)
#define ASTC_VECMATHLIB_COMMON_8_H_INCLUDED

#if !defined(ASTCENC_SIMD_INLINE)
	#error "Include astcenc_vecmathlib.h, do not include directly"
#endif

#include <cstdio>
#include <cstring>
#include <limits>

// ============================================================================
// vint8 operators and functions
// ============================================================================

/**
 * @brief Overload: vector by vector incremental addition.
 *
 * @param[in,out] a   The accumulator to update.
 * @param         b   The values to add.
 */
ASTCENC_SIMD_INLINE vint8& operator+=(vint8& a, const vint8& b)
{
	a = a + b;
	return a;
}

/**
 * @brief Generate a vint8 from a size_t.
 *
 * @param a   The scalar size to replicate across all lanes.
 */
ASTCENC_SIMD_INLINE vint8 vint8_from_size(size_t a)
{
	assert(a <= std::numeric_limits<int>::max());
	return vint8(static_cast<int>(a));
}

// ============================================================================
// vfloat8 operators and functions
// ============================================================================

/**
 * @brief Overload: vector by vector incremental addition.
 *
 * @param[in,out] a   The accumulator to update.
 * @param         b   The values to add.
 */
ASTCENC_SIMD_INLINE vfloat8& operator+=(vfloat8& a, const vfloat8& b)
{
	a = a + b;
	return a;
}

#if !defined(ASTCENC_USE_NATIVE_VFLOAT8_SCALAR_ARITH)

/**
 * @brief Overload: vector by scalar multiplication.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat8 operator*(vfloat8 a, float b)
{
	return a * vfloat8(b);
}

/**
 * @brief Overload: scalar by vector multiplication.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat8 operator*(float a, vfloat8 b)
{
	return vfloat8(a) * b;
}

/**
 * @brief Overload: vector by scalar division.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat8 operator/(vfloat8 a, float b)
{
	return a / vfloat8(b);
}

/**
 * @brief Overload: scalar by vector division.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat8 operator/(float a, vfloat8 b)
{
	return vfloat8(a) / b;
}

#endif

/**
 * @brief Return the min vector of a vector and a scalar.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vfloat8 min(vfloat8 a, float b)
{
	return min(a, vfloat8(b));
}

/**
 * @brief Return the max vector of a vector and a scalar.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vfloat8 max(vfloat8 a, float b)
{
	return max(a, vfloat8(b));
}

/**
 * @brief Return the clamped value between min and max.
 *
 * It is assumed that neither @c min nor @c max are NaN values. If @c a is NaN
 * then @c min will be returned for that lane.
 *
 * @param minv   The inclusive lower bound.
 * @param maxv   The inclusive upper bound.
 * @param a      The values to clamp.
 */
ASTCENC_SIMD_INLINE vfloat8 clamp(float minv, float maxv, vfloat8 a)
{
	// Preserve operation and operand order to give desired NaN handling
	return min(max(a, minv), maxv);
}

/**
 * @brief Return the clamped value between 0.0f and 1.0f.
 *
 * If @c a is NaN then zero will be returned for that lane.
 *
 * @param a   The values to clamp.
 */
ASTCENC_SIMD_INLINE vfloat8 clampzo(vfloat8 a)
{
	return clamp(0.0f, 1.0f, a);
}

/**
 * @brief Accumulate lane-wise sums for a vector.
 *
 * This is NOT invariant with 4-wide implementations.
 *
 * @param[in,out] accum   The accumulator to update.
 * @param         a       The values to accumulate.
 */
ASTCENC_SIMD_INLINE void haccumulate(vfloat8& accum, vfloat8 a)
{
	accum += a;
}

/**
 * @brief Accumulate masked lane-wise sums for a vector, folded 4-wide.
 *
 * This is invariant with 4-wide implementations.
 *
 * @param[in,out] accum   The accumulator to update.
 * @param         a       The values to accumulate.
 * @param         m       The lanes to include in the accumulation.
 */
ASTCENC_SIMD_INLINE void haccumulate(vfloat4& accum, vfloat8 a, vmask8 m)
{
	a = select(vfloat8::zero(), a, m);
	haccumulate(accum, a);
}

#if !defined(ASTCENC_USE_NATIVE_FLOAT_TO_INT_RTN_8)

/**
 * @brief Return integers rounded to nearest for nonnegative float values.
 *
 * @param a   The nonnegative floating-point values to convert.
 */
ASTCENC_SIMD_INLINE vint8 float_to_int_rtn(vfloat8 a)
{
	return float_to_int(a + vfloat8(0.5f));
}

#endif

/**
 * @brief Return a vector of interleaved RGBA data.
 *
 * Input vectors have the value stored in the bottom 8 bits of each lane,
 * with high bits set to zero.
 *
 * Output vector stores a single RGBA texel packed in each lane.
 *
 * @param r   The red channel values.
 * @param g   The green channel values.
 * @param b   The blue channel values.
 * @param a   The alpha channel values.
 */
ASTCENC_SIMD_INLINE vint8 interleave_rgba8(vint8 r, vint8 g, vint8 b, vint8 a)
{
	return r + lsl<8>(g) + lsl<16>(b) + lsl<24>(a);
}

/**
 * @brief Debug function to print a vector of ints.
 *
 * @param a   The vector to print.
 */
ASTCENC_SIMD_INLINE void print(vint8 a)
{
	ASTCENC_ALIGNAS int v[8];
	storea(a, v);
	printf("v8_i32:\n  %8d %8d %8d %8d %8d %8d %8d %8d\n",
	       v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]);
}

/**
 * @brief Debug function to print a vector of ints.
 *
 * @param a   The vector to print.
 */
ASTCENC_SIMD_INLINE void printx(vint8 a)
{
	ASTCENC_ALIGNAS int v[8];
	storea(a, v);

	unsigned int uv[8];
	std::memcpy(uv, v, sizeof(int) * 8);

	printf("v8_i32:\n  %08x %08x %08x %08x %08x %08x %08x %08x\n",
	       uv[0], uv[1], uv[2], uv[3], uv[4], uv[5], uv[6], uv[7]);
}

/**
 * @brief Debug function to print a vector of floats.
 *
 * @param a   The vector to print.
 */
ASTCENC_SIMD_INLINE void print(vfloat8 a)
{
	ASTCENC_ALIGNAS float v[8];
	storea(a, v);
	printf("v8_f32:\n  %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f\n",
	       static_cast<double>(v[0]), static_cast<double>(v[1]),
	       static_cast<double>(v[2]), static_cast<double>(v[3]),
	       static_cast<double>(v[4]), static_cast<double>(v[5]),
	       static_cast<double>(v[6]), static_cast<double>(v[7]));
}

/**
 * @brief Debug function to print a vector of masks.
 *
 * @param a   The vector to print.
 */
ASTCENC_SIMD_INLINE void print(vmask8 a)
{
	print(select(vint8(0), vint8(1), a));
}

#endif  // #ifndef ASTC_VECMATHLIB_COMMON_8_H_INCLUDED
