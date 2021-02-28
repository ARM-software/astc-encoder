// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2019-2021 Arm Limited
// Copyright 2008 Jose Fonseca
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

/*
 * This module implements vector support for floats, ints, and vector lane
 * control masks. It provides access to both explicit vector width types, and
 * flexible N-wide types where N can be determined at compile time.
 *
 * The design of this module encourages use of vector length agnostic code, via
 * the vint, vfloat, and vmask types. These will take on the widest SIMD vector
 * with that is available at compile time. The current vector width is
 * accessible for e.g. loop strides via the ASTCENC_SIMD_WIDTH constant.
 *
 * Explicit scalar types are acessible via the vint1, vfloat1, vmask1 types.
 * These are provided primarily for prototyping and algorithm debug of VLA
 * implementations.
 *
 * Explicit 4-wide types are accessible via the vint4, vfloat4, and vmask4
 * types. These are provided for use by VLA code, but are also expected to be
 * used as a fixed-width type and will supported a reference C++ fallback for
 * use on platforms without SIMD intrinsics.
 *
 * Explicit 8-wide types are accessible via the vint8, vfloat8, and vmask8
 * types. These are provide for use by VLA code, and are not expected to be
 * used as a fixed-width type in normal code. No reference C implementation is
 * provided on platforms without underlying SIMD intrinsics.
 *
 * With the current implementation ISA support is provided for:
 *
 *     * 1-wide for scalar reference.
 *     * 4-wide for Armv8-A NEON.
 *     * 4-wide for x86-64 SSE2.
 *     * 4-wide for x86-64 SSE4.1.
 *     * 8-wide for x86-64 AVX2.
 *
 */

#ifndef ASTC_VECMATHLIB_H_INCLUDED
#define ASTC_VECMATHLIB_H_INCLUDED

#if ASTCENC_SSE != 0 || ASTCENC_AVX != 0
	#include <immintrin.h>
#elif ASTCENC_NEON != 0
	#include <arm_neon.h>
#endif

#if !defined(__clang__) && defined(_MSC_VER)
	#define ASTCENC_SIMD_INLINE __forceinline
#elif defined(__GNUC__) && !defined(__clang__)
	#define ASTCENC_SIMD_INLINE __attribute__((always_inline)) inline
#else
	#define ASTCENC_SIMD_INLINE __attribute__((always_inline, nodebug)) inline
#endif

#if ASTCENC_AVX >= 2
	/* If we have AVX2 expose 8-wide VLA. */
	#include "astcenc_vecmathlib_sse_4.h"
	#include "astcenc_vecmathlib_avx2_8.h"

	#define ASTCENC_SIMD_WIDTH 8

	using vfloat = vfloat8;
	using vint = vint8;
	using vmask = vmask8;

	constexpr auto loada = vfloat8::loada;
	constexpr auto load1 = vfloat8::load1;

#elif ASTCENC_SSE >= 20
	/* If we have SSE expose 4-wide VLA, and 4-wide fixed width. */
	#include "astcenc_vecmathlib_sse_4.h"

	#define ASTCENC_SIMD_WIDTH 4

	using vfloat = vfloat4;
	using vint = vint4;
	using vmask = vmask4;

	constexpr auto loada = vfloat4::loada;
	constexpr auto load1 = vfloat4::load1;

#elif ASTCENC_NEON > 0
	/* If we have NEON expose 4-wide VLA. */
	#include "astcenc_vecmathlib_neon_4.h"

	#define ASTCENC_SIMD_WIDTH 4

	using vfloat = vfloat4;
	using vint = vint4;
	using vmask = vmask4;

	constexpr auto loada = vfloat4::loada;
	constexpr auto load1 = vfloat4::load1;

#else
	// If we have nothing expose 4-wide VLA, and 4-wide fixed width.

	// Note: We no longer expose the 1-wide scalar fallback because it is not
	// invariant with the 4-wide path due to algorithms that use horizontal
	// operations that accumulate a local vector sum before accumulating into
	// a running sum.
	//
	// For 4 items adding into an accumulator using 1-wide vectors the sum is:
	//
	//     result = ((((sum + l0) + l1) + l2) + l3)
	//
    // ... whereas the accumulator for a 4-wide vector sum is:
	//
	//     result = sum + ((l0 + l2) + (l1 + l3))
	//
	// In "normal maths" this is the same, but the floating point reassociation
	// differences mean that these will not produce the same result.

	#include "astcenc_vecmathlib_none_4.h"

	#define ASTCENC_SIMD_WIDTH 4

	using vfloat = vfloat4;
	using vint = vint4;
	using vmask = vmask4;

	constexpr auto loada = vfloat4::loada;
	constexpr auto load1 = vfloat4::load1;
#endif

/**
 * @brief Round a count down to the largest multiple of 8.
 *
 * @param count   The unrounded value.
 *
 * @return The rounded value.
 */
ASTCENC_SIMD_INLINE int round_down_to_simd_multiple_8(int count)
{
	return count & ~(8 - 1);
}

/**
 * @brief Round a count down to the largest multiple of 4.
 *
 * @param count   The unrounded value.
 *
 * @return The rounded value.
 */
ASTCENC_SIMD_INLINE int round_down_to_simd_multiple_4(int count)
{
	return count & ~(4 - 1);
}

/**
 * @brief Round a count down to the largest multiple of the SIMD width.
 *
 * Assumption that the vector width is a power of two ...
 *
 * @param count   The unrounded value.
 *
 * @return The rounded value.
 */
ASTCENC_SIMD_INLINE int round_down_to_simd_multiple_vla(int count)
{
	return count & ~(ASTCENC_SIMD_WIDTH - 1);
}

/**
 * @brief Round a count up to the largest multiple of the SIMD width.
 *
 * Assumption that the vector width is a power of two ...
 *
 * @param count   The unrounded value.
 *
 * @return The rounded value.
 */
ASTCENC_SIMD_INLINE int round_up_to_simd_multiple_vla(int count)
{
	int multiples = (count + ASTCENC_SIMD_WIDTH - 1) / ASTCENC_SIMD_WIDTH;
	return multiples * ASTCENC_SIMD_WIDTH;
}

/**
 * @brief Return @c a with lanes negated if the @c b lane is negative.
 */
ASTCENC_SIMD_INLINE vfloat change_sign(vfloat a, vfloat b)
{
	vint ia = float_as_int(a);
	vint ib = float_as_int(b);
	vint sign_mask((int)0x80000000);
	vint r = ia ^ (ib & sign_mask);
	return int_as_float(r);
}

/**
 * @brief Return fast, but approximate, vector atan(x).
 *
 * Max error of this implementaiton is 0.004883.
 */
ASTCENC_SIMD_INLINE vfloat atan(vfloat x)
{
	vmask c = abs(x) > vfloat(1.0f);
	vfloat z = change_sign(vfloat(astc::PI_OVER_TWO), x);
	vfloat y = select(x, vfloat(1.0f) / x, c);
	y = y / (y * y * vfloat(0.28f) + vfloat(1.0f));
	return select(y, z - y, c);
}

/**
 * @brief Return fast, but approximate, vector atan2(x, y).
 */
ASTCENC_SIMD_INLINE vfloat atan2(vfloat y, vfloat x)
{
	vfloat z = atan(abs(y / x));
	vmask xmask = vmask(float_as_int(x).m);
	return change_sign(select(z, vfloat(astc::PI) - z, xmask), y);
}

/*
 * @brief Factory that returns a unit length 4 component vfloat4.
 */
static ASTCENC_SIMD_INLINE vfloat4 unit4()
{
	return vfloat4(0.5f);
}

/**
 * @brief Factory that returns a unit length 3 component vfloat4.
 */
static ASTCENC_SIMD_INLINE vfloat4 unit3()
{
	return vfloat4(0.57735f, 0.57735f, 0.57735f, 0.0f);
}

/**
 * @brief Normalize a non-zero length vector to unit length.
 */
static ASTCENC_SIMD_INLINE vfloat4 normalize(vfloat4 a)
{
	vfloat4 length = dot(a, a);
	return a / sqrt(length);
}

/**
 * @brief Normalize a vector, returning @c safe if len is zero.
 */
static ASTCENC_SIMD_INLINE vfloat4 normalize_safe(vfloat4 a, vfloat4 safe)
{
	vfloat4 length = dot(a, a);
	if (length.lane<0>() != 0.0f)
	{
		return a / sqrt(length);
	}

	return safe;
}

// log2 and exp2 based on 5th degree minimax polynomials, ported from this blog
// https://jrfonseca.blogspot.com/2008/09/fast-sse2-pow-tables-or-polynomials.html

#define POLY0(x, c0)                     (                                     c0)
#define POLY1(x, c0, c1)                 ((POLY0(x, c1) * x)                 + c0)
#define POLY2(x, c0, c1, c2)             ((POLY1(x, c1, c2) * x)             + c0)
#define POLY3(x, c0, c1, c2, c3)         ((POLY2(x, c1, c2, c3) * x)         + c0)
#define POLY4(x, c0, c1, c2, c3, c4)     ((POLY3(x, c1, c2, c3, c4) * x)     + c0)
#define POLY5(x, c0, c1, c2, c3, c4, c5) ((POLY4(x, c1, c2, c3, c4, c5) * x) + c0)

static ASTCENC_SIMD_INLINE vfloat4 exp2(vfloat4 x)
{
	x = clamp(-126.99999f, 129.0f, x);

	vint4 ipart = float_to_int(x - 0.5f);
	vfloat4 fpart = x - int_to_float(ipart);

	// Integer contrib, using 1 << ipart
	vfloat4 iexp = int_as_float(lsl<23>(ipart + 127));

	// Fractional contrib, using polynomial fit of 2^x in range [-0.5, 0.5)
	vfloat4 fexp = POLY5(fpart,
	                     9.9999994e-1f,
	                     6.9315308e-1f,
	                     2.4015361e-1f,
	                     5.5826318e-2f,
	                     8.9893397e-3f,
	                     1.8775767e-3f);

	return iexp * fexp;
}

static ASTCENC_SIMD_INLINE vfloat4 log2(vfloat4 x)
{
	vint4 exp(0x7F800000);
	vint4 mant(0x007FFFFF);
	vint4 one(0x3F800000);

	vint4 i = float_as_int(x);

	vfloat4 e = int_to_float(lsr<23>(i & exp) - 127);

	vfloat4 m = int_as_float((i & mant) | one);

	// Polynomial fit of log2(x)/(x - 1), for x in range [1, 2)
	vfloat4 p = POLY4(m,
	                  2.8882704548164776201f,
	                 -2.52074962577807006663f,
	                  1.48116647521213171641f,
	                 -0.465725644288844778798f,
	                  0.0596515482674574969533f);

	// Increases the polynomial degree, but ensures that log2(1) == 0
	p = p * (m - 1.0f);

	return p + e;
}

static ASTCENC_SIMD_INLINE vfloat4 pow(vfloat4 x, vfloat4 y)
{
	vmask4 zero_mask = y == vfloat4(0.0f);
	vfloat4 estimate = exp2(log2(x) * y);

	// Guarantee that y == 0 returns exactly 1.0f
	return select(estimate, vfloat4(1.0f), zero_mask);
}

namespace astc
{

static ASTCENC_SIMD_INLINE float pow(float x, float y)
{
	return pow(vfloat4(x), vfloat4(y)).lane<0>();
}

}

#endif // #ifndef ASTC_VECMATHLIB_H_INCLUDED
