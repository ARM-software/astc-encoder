// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2019-2026 Arm Limited
// Copyright 2008 Jose Fonseca
// Copyright 2026 Olaf Bernstein
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
 * Explicit scalar types are accessible via the vint1, vfloat1, vmask1 types.
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
 *     * 1-wide for scalar reference
 *     * 4-wide for Armv8-A NEON
 *     * 4-wide for Armv8-A NEON augented with 128-bit SVE
 *     * 4-wide for x86-64 SSE2
 *     * 4-wide for x86-64 SSE4.1
 *     * 8-wide for Armv8-A 256-bit SVE
 *     * 8-wide for x86-64 AVX2
 *     * 16-wide for x86-64 AVX-512
 */

#if !defined(ASTC_VECMATHLIB_H_INCLUDED)
#define ASTC_VECMATHLIB_H_INCLUDED

#if ASTCENC_SSE != 0 || ASTCENC_AVX != 0
	#include <immintrin.h>
#endif

#if ASTCENC_SVE != 0
	#include <arm_sve.h>
	#include <arm_neon_sve_bridge.h>
#endif

#if ASTCENC_NEON != 0
	#include <arm_neon.h>
#endif

#if !defined(__clang__) && defined(_MSC_VER)
	#define ASTCENC_SIMD_INLINE __forceinline
	#define ASTCENC_NO_INLINE
#elif defined(__GNUC__) && !defined(__clang__)
	#define ASTCENC_SIMD_INLINE __attribute__((always_inline)) inline
	#define ASTCENC_NO_INLINE __attribute__ ((noinline))
#else
	#define ASTCENC_SIMD_INLINE __attribute__((always_inline, nodebug)) inline
	#define ASTCENC_NO_INLINE __attribute__ ((noinline))
#endif

/**
 * @brief Gather float values using byte-sized indices.
 *
 * @tparam T         The vector type of the gathered values.
 *
 * @param  base      The source float table.
 * @param  indices   The byte-sized table indices, one per lane.
 */
template<typename T> T gatherf_byte_inds(const float* base, const uint8_t* indices);

#if ASTCENC_AVX >= 4
	// If we have 512-bit AVX-512, expose 16-wide VLA and 4-wide fixed width
	#include "astcenc_vecmathlib_sse_4.h"
	#include "astcenc_vecmathlib_common_4.h"
	#include "astcenc_vecmathlib_avx512_16.h"

	#define ASTCENC_SIMD_WIDTH 16

	using vfloat = vfloat16;

	#if defined(ASTCENC_NO_INVARIANCE)
		using vfloatacc = vfloat16;
	#else
		using vfloatacc = vfloat4;
	#endif

	using vint = vint16;
	using vmask = vmask16;

	using vtable_16x8 = vtable16_16x8;
	using vtable_32x8 = vtable16_32x8;
	using vtable_64x8 = vtable16_64x8;

	constexpr auto loada = vfloat16::loada;
	constexpr auto load1 = vfloat16::load1;
	constexpr auto vint_from_size = vint16_from_size;

#elif ASTCENC_AVX >= 2
	// If we have AVX2 or 256-bit AVX-512, expose 8-wide VLA and 4-wide fixed width
	#include "astcenc_vecmathlib_sse_4.h"
	#include "astcenc_vecmathlib_common_4.h"
	#include "astcenc_vecmathlib_avx2_8.h"
	#include "astcenc_vecmathlib_common_8.h"

	#define ASTCENC_SIMD_WIDTH 8

	using vfloat = vfloat8;

	#if defined(ASTCENC_NO_INVARIANCE)
		using vfloatacc = vfloat8;
	#else
		using vfloatacc = vfloat4;
	#endif

	using vint = vint8;
	using vmask = vmask8;

	using vtable_16x8 = vtable8_16x8;
	using vtable_32x8 = vtable8_32x8;
	using vtable_64x8 = vtable8_64x8;

	constexpr auto loada = vfloat8::loada;
	constexpr auto load1 = vfloat8::load1;
	constexpr auto vint_from_size = vint8_from_size;

#elif ASTCENC_SSE >= 20
	// If we have SSE, expose 4-wide VLA and 4-wide fixed width
	#include "astcenc_vecmathlib_sse_4.h"
	#include "astcenc_vecmathlib_common_4.h"

	#define ASTCENC_SIMD_WIDTH 4

	using vfloat = vfloat4;
	using vfloatacc = vfloat4;
	using vint = vint4;
	using vmask = vmask4;

	using vtable_16x8 = vtable4_16x8;
	using vtable_32x8 = vtable4_32x8;
	using vtable_64x8 = vtable4_64x8;

	constexpr auto loada = vfloat4::loada;
	constexpr auto load1 = vfloat4::load1;
	constexpr auto vint_from_size = vint4_from_size;

#elif ASTCENC_SVE == 8
	// Check the compiler is configured with fixed-length 256-bit SVE
	#if !defined(__ARM_FEATURE_SVE_BITS) || (__ARM_FEATURE_SVE_BITS != 256)
		#error "__ARM_FEATURE_SVE_BITS is not set to 256 bits"
	#endif

	// If we have 256-bit SVE, expose 8-wide VLA and 4-wide fixed width
	#include "astcenc_vecmathlib_neon_4.h"
	#include "astcenc_vecmathlib_common_4.h"
	#include "astcenc_vecmathlib_sve_8.h"
	#include "astcenc_vecmathlib_common_8.h"

	#define ASTCENC_SIMD_WIDTH 8

	using vfloat = vfloat8;

	#if defined(ASTCENC_NO_INVARIANCE)
		using vfloatacc = vfloat8;
	#else
		using vfloatacc = vfloat4;
	#endif

	using vint = vint8;
	using vmask = vmask8;

	using vtable_16x8 = vtable8_16x8;
	using vtable_32x8 = vtable8_32x8;
	using vtable_64x8 = vtable8_64x8;

	constexpr auto loada = vfloat8::loada;
	constexpr auto load1 = vfloat8::load1;
	constexpr auto vint_from_size = vint8_from_size;

#elif ASTCENC_NEON > 0
	// If we have NEON, expose 4-wide VLA and 4-wide fixed width
	#include "astcenc_vecmathlib_neon_4.h"
	#include "astcenc_vecmathlib_common_4.h"

	#define ASTCENC_SIMD_WIDTH 4

	using vfloat = vfloat4;
	using vfloatacc = vfloat4;
	using vint = vint4;
	using vmask = vmask4;

	using vtable_16x8 = vtable4_16x8;
	using vtable_32x8 = vtable4_32x8;
	using vtable_64x8 = vtable4_64x8;

	constexpr auto loada = vfloat4::loada;
	constexpr auto load1 = vfloat4::load1;
	constexpr auto vint_from_size = vint4_from_size;

#elif defined(__riscv_v_fixed_vlen) && (__riscv_v_fixed_vlen == 256)
	// If we have 256-bit RVV, expose 8-wide VLA and 4-wide fixed width
	#include <riscv_vector.h>
	#include "astcenc_vecmathlib_none_4.h"
	#include "astcenc_vecmathlib_common_4.h"
	#include "astcenc_vecmathlib_rvv_n.h"

	#if defined(ASTCENC_NO_INVARIANCE)
		using vfloatacc = vfloat;
	#else
		using vfloatacc = vfloat4;
	#endif

	constexpr auto loada = vfloat::loada;
	constexpr auto load1 = vfloat::load1;

	// For unit-tests
	using vfloat8 = vfloat;
	using vint8 = vint;
	using vmask8 = vmask;
	using vtable8_16x8 = vtable_16x8;
	using vtable8_32x8 = vtable_32x8;
	using vtable8_64x8 = vtable_64x8;

#else
	// If we have no SIMD, expose 4-wide VLA and 4-wide fixed width

	// The scalar fallback uses four lanes to preserve the SIMD accumulation order:
	//
	// One-lane accumulation:
	//     result = ((((sum + l0) + l1) + l2) + l3)
	//
	// Four-lane accumulation:
	//     result = sum + ((l0 + l2) + (l1 + l3))
	//
	// These are algebraically equivalent, but floating-point reassociation can change the result

	#include "astcenc_vecmathlib_none_4.h"
	#include "astcenc_vecmathlib_common_4.h"

	#define ASTCENC_SIMD_WIDTH 4

	using vfloat = vfloat4;
	using vfloatacc = vfloat4;
	using vint = vint4;
	using vmask = vmask4;

	using vtable_16x8 = vtable4_16x8;
	using vtable_32x8 = vtable4_32x8;
	using vtable_64x8 = vtable4_64x8;

	constexpr auto loada = vfloat4::loada;
	constexpr auto load1 = vfloat4::load1;
	constexpr auto vint_from_size = vint4_from_size;
#endif

#if !defined(ASTCENC_HAS_NATIVE_QUAD_OPS)

/**
 * @brief Store lane zero from each quad as consecutive floats.
 *
 * @param      a   The input vector.
 * @param[out] p   The unaligned output array of ASTCENC_SIMD_WIDTH / 4 floats.
 */
ASTCENC_SIMD_INLINE void store_lane0_quads(vfloat a, float* p)
{
#if ASTCENC_SIMD_WIDTH == 4
	p[0] = a.lane<0>();
#else
	ASTCENC_ALIGNAS float values[ASTCENC_SIMD_WIDTH];
	storea(a, values);
	for (unsigned int c = 0; c < ASTCENC_SIMD_WIDTH / 4; c++)
	{
		p[c] = values[c * 4];
	}
#endif
}

/**
 * @brief Pack consecutive four-lane float vectors into a VLA vector.
 *
 * A quad is four consecutive lanes. Quad c occupies lanes 4*c through 4*c+3.
 * Four-wide backends return quad zero.
 *
 * @param groups   The ASTCENC_SIMD_WIDTH / 4 source groups.
 *
 * @return The packed vector with source groups in consecutive quads.
 */
ASTCENC_SIMD_INLINE vfloat pack_quads(const vfloat4* groups)
{
#if ASTCENC_SIMD_WIDTH == 4
	return groups[0];
#else
	ASTCENC_ALIGNAS float values[ASTCENC_SIMD_WIDTH];
	for (unsigned int c = 0; c < ASTCENC_SIMD_WIDTH / 4; c++)
	{
		storea(groups[c], values + c * 4);
	}
	return vfloat::loada(values);
#endif
}

/**
 * @brief Pack consecutive four-lane integer vectors into a VLA vector.
 *
 * A quad is four consecutive lanes. Quad c occupies lanes 4*c through 4*c+3.
 * Four-wide backends return quad zero.
 *
 * @param groups   The ASTCENC_SIMD_WIDTH / 4 source groups.
 *
 * @return The packed vector with source groups in consecutive quads.
 */
ASTCENC_SIMD_INLINE vint pack_quads(const vint4* groups)
{
#if ASTCENC_SIMD_WIDTH == 4
	return groups[0];
#else
	ASTCENC_ALIGNAS int values[ASTCENC_SIMD_WIDTH];
	for (unsigned int c = 0; c < ASTCENC_SIMD_WIDTH / 4; c++)
	{
		storea(groups[c], values + c * 4);
	}
	return vint::loada(values);
#endif
}

/**
 * @brief Apply the same swizzle independently within each quad of four lanes.
 *
 * Indices are relative to each group; values never cross group boundaries.
 *
 * @tparam l0   The source lane for output lane 0 in each group.
 * @tparam l1   The source lane for output lane 1 in each group.
 * @tparam l2   The source lane for output lane 2 in each group.
 * @tparam l3   The source lane for output lane 3 in each group.
 *
 * @param a   The input vector.
 *
 * @return The vector with lanes swizzled within each quad.
 */
template <int l0, int l1, int l2, int l3>
ASTCENC_SIMD_INLINE vfloat swizzle_quads(vfloat a)
{
	static_assert(l0 >= 0 && l0 < 4 && l1 >= 0 && l1 < 4 &&
	              l2 >= 0 && l2 < 4 && l3 >= 0 && l3 < 4, "Invalid four-lane swizzle");
#if ASTCENC_SIMD_WIDTH == 4
	return a.swz<l0, l1, l2, l3>();
#else
	ASTCENC_ALIGNAS float values[ASTCENC_SIMD_WIDTH];
	ASTCENC_ALIGNAS float result[ASTCENC_SIMD_WIDTH];
	storea(a, values);
	for (unsigned int c = 0; c < ASTCENC_SIMD_WIDTH; c += 4)
	{
		result[c] = values[c + l0];
		result[c + 1] = values[c + l1];
		result[c + 2] = values[c + l2];
		result[c + 3] = values[c + l3];
	}
	return vfloat::loada(result);
#endif
}
#endif

/**
 * @brief Sum each quad independently, replicating the result within its four lanes.
 *
 * The reduction order matches hadd_s(vfloat4).
 *
 * @param a   The input vector.
 *
 * @return The sums replicated within their respective quads.
 */
ASTCENC_SIMD_INLINE vfloat sum_quads(vfloat a)
{
	vfloat pairs = a + swizzle_quads<2, 3, 0, 1>(a);
	return pairs + swizzle_quads<1, 0, 3, 2>(pairs);
}

/**
 * @brief Compute a four-lane dot product independently for each quad.
 *
 * Each result is replicated within its quad using the invariant four-lane reduction order.
 *
 * @param a   The first vector.
 * @param b   The second vector.
 *
 * @return The dot products replicated within their respective quads.
 */
ASTCENC_SIMD_INLINE vfloat dot_quads(vfloat a, vfloat b)
{
	return sum_quads(a * b);
}

/**
 * @brief Normalize each quad, returning @c safe for zero-length quads.
 *
 * @param a      The vector to normalize.
 * @param safe   The fallback vector for zero-length quads.
 *
 * @return The normalized vectors or fallback vectors in each quad.
 */
ASTCENC_SIMD_INLINE vfloat normalize_safe_quads(vfloat a, vfloat safe)
{
	vfloat length = dot_quads(a, a);
	vmask is_zero = length == vfloat::zero();
	vfloat normalized = a / sqrt(length);
	return select(normalized, safe, is_zero);
}

/**
 * @brief Sum RGB lanes in each quad, replicating the result within its four lanes.
 *
 * The reduction order matches hadd_rgb_s(vfloat4).
 *
 * @param a   The input vector.
 *
 * @return The RGB sums replicated within their respective quads.
 */
ASTCENC_SIMD_INLINE vfloat sum_rgb_quads(vfloat a)
{
	return (swizzle_quads<0, 0, 0, 0>(a) + swizzle_quads<1, 1, 1, 1>(a)) + swizzle_quads<2, 2, 2, 2>(a);
}

/**
 * @brief Find the minimum in each quad, replicating the result within its four lanes.
 *
 * @param a   The input vector.
 *
 * @return The minima replicated within their respective quads.
 */
ASTCENC_SIMD_INLINE vfloat min_quads(vfloat a)
{
	vfloat pairs = min(a, swizzle_quads<2, 3, 0, 1>(a));
	return min(pairs, swizzle_quads<1, 0, 3, 2>(pairs));
}

/**
 * @brief Find the maximum in each quad, replicating the result within its four lanes.
 *
 * @param a   The input vector.
 *
 * @return The maxima replicated within their respective quads.
 */
ASTCENC_SIMD_INLINE vfloat max_quads(vfloat a)
{
	vfloat pairs = max(a, swizzle_quads<2, 3, 0, 1>(a));
	return max(pairs, swizzle_quads<1, 0, 3, 2>(pairs));
}

#if !defined(ASTCENC_HAS_VGATHERF_TABLE)

/**
 * @brief Portable float-table gather: keep the pointer, gather from memory.
 *
 * AVX-512 replaces this with a register-resident float table in
 * astcenc_vecmathlib_avx512_16.h.
 */
struct vgatherf_table {
	/** @brief The source table, used for gathers from memory. */
	const float* base;
};

/**
 * @brief Prepare a float table for portable gathers from memory.
 *
 * The source table must remain valid while the returned table is in use.
 *
 * @param base    The source float table.
 * @param count   The number of entries; unused by this implementation.
 *
 * @return The table representation to use with gatherf().
 */
ASTCENC_SIMD_INLINE vgatherf_table vgatherf_load(const float* base, unsigned int count)
{
	(void)count;
	return vgatherf_table{base};
}

/**
 * @brief Gather from a prepared float table using vector indices.
 *
 * @param t         The prepared float table.
 * @param indices   The table index for each lane.
 */
ASTCENC_SIMD_INLINE vfloat gatherf(const vgatherf_table& t, vint indices)
{
	return gatherf(t.base, indices);
}

/**
 * @brief Gather from a prepared float table using byte indices from memory.
 *
 * @param t         The prepared float table.
 * @param indices   The byte-sized table indices, one per lane.
 */
ASTCENC_SIMD_INLINE vfloat gatherf(const vgatherf_table& t, const uint8_t* indices)
{
	return gatherf_byte_inds<vfloat>(t.base, indices);
}
#endif

/**
 * @brief Round a count down to the largest multiple of the SIMD width.
 *
 * Assumption that the vector width is a power of two ...
 *
 * @param count   The unrounded value.
 *
 * @return The rounded value.
 */
ASTCENC_SIMD_INLINE size_t round_down_to_simd_multiple_vla(size_t count)
{
	return count & static_cast<size_t>(~(ASTCENC_SIMD_WIDTH - 1));
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
ASTCENC_SIMD_INLINE size_t round_up_to_simd_multiple_vla(size_t count)
{
	size_t multiples = (count + ASTCENC_SIMD_WIDTH - 1) / ASTCENC_SIMD_WIDTH;
	return multiples * ASTCENC_SIMD_WIDTH;
}

/**
 * @brief Return @c a with lanes negated if the @c b lane is negative.
 *
 * @param a   The values to conditionally negate.
 * @param b   The values whose sign bits control negation.
 */
ASTCENC_SIMD_INLINE vfloat change_sign(vfloat a, vfloat b)
{
	vint ia = float_as_int(a);
	vint ib = float_as_int(b);
	vint sign_mask(static_cast<int>(0x80000000));
	vint r = ia ^ (ib & sign_mask);
	return int_as_float(r);
}

/**
 * @brief Return fast, but approximate, vector atan(x).
 *
 * Max error of this implementation is 0.004883.
 *
 * @param x   The values whose arctangents are required.
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
 *
 * @param y   The y coordinates.
 * @param x   The x coordinates.
 */
ASTCENC_SIMD_INLINE vfloat atan2(vfloat y, vfloat x)
{
	vfloat z = atan(abs(y / x));
	vmask xmask = x < vfloat::zero();
	return change_sign(select(z, vfloat(astc::PI) - z, xmask), y);
}

/**
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
	float val = 0.577350258827209473f;
	return vfloat4(val, val, val, 0.0f);
}

/**
 * @brief Factory that returns a unit length 2 component vfloat4.
 */
static ASTCENC_SIMD_INLINE vfloat4 unit2()
{
	float val = 0.707106769084930420f;
	return vfloat4(val, val, 0.0f, 0.0f);
}

/**
 * @brief Factory that returns a 3 component vfloat4.
 *
 * @param a   The value for lane 0.
 * @param b   The value for lane 1.
 * @param c   The value for lane 2.
 */
static ASTCENC_SIMD_INLINE vfloat4 vfloat3(float a, float b, float c)
{
	return vfloat4(a, b, c, 0.0f);
}

/**
 * @brief Factory that returns a 2 component vfloat4.
 *
 * @param a   The value for lane 0.
 * @param b   The value for lane 1.
 */
static ASTCENC_SIMD_INLINE vfloat4 vfloat2(float a, float b)
{
	return vfloat4(a, b, 0.0f, 0.0f);
}

/**
 * @brief Normalize a non-zero length vector to unit length.
 *
 * @param a   The non-zero vector to normalize.
 */
static ASTCENC_SIMD_INLINE vfloat4 normalize(vfloat4 a)
{
	vfloat4 length = dot(a, a);
	return a / sqrt(length);
}

/**
 * @brief Normalize a vector, returning @c safe if len is zero.
 *
 * @param a      The vector to normalize.
 * @param safe   The fallback vector if @c a has zero length.
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



/**
 * @brief Evaluate a degree-0 polynomial using Horner's method.
 *
 * @param x    The polynomial argument; unused for a constant.
 * @param c0   The constant coefficient.
 */
#define POLY0(x, c0)                     (                                     c0)

/**
 * @brief Evaluate a degree-1 polynomial using Horner's method.
 *
 * @param x    The polynomial argument.
 * @param c0   The constant coefficient.
 * @param c1   The coefficient of x^1.
 */
#define POLY1(x, c0, c1)                 ((POLY0(x, c1) * x)                 + c0)

/**
 * @brief Evaluate a degree-2 polynomial using Horner's method.
 *
 * @param x    The polynomial argument.
 * @param c0   The constant coefficient.
 * @param c1   The coefficient of x^1.
 * @param c2   The coefficient of x^2.
 */
#define POLY2(x, c0, c1, c2)             ((POLY1(x, c1, c2) * x)             + c0)

/**
 * @brief Evaluate a degree-3 polynomial using Horner's method.
 *
 * @param x    The polynomial argument.
 * @param c0   The constant coefficient.
 * @param c1   The coefficient of x^1.
 * @param c2   The coefficient of x^2.
 * @param c3   The coefficient of x^3.
 */
#define POLY3(x, c0, c1, c2, c3)         ((POLY2(x, c1, c2, c3) * x)         + c0)

/**
 * @brief Evaluate a degree-4 polynomial using Horner's method.
 *
 * @param x    The polynomial argument.
 * @param c0   The constant coefficient.
 * @param c1   The coefficient of x^1.
 * @param c2   The coefficient of x^2.
 * @param c3   The coefficient of x^3.
 * @param c4   The coefficient of x^4.
 */
#define POLY4(x, c0, c1, c2, c3, c4)     ((POLY3(x, c1, c2, c3, c4) * x)     + c0)

/**
 * @brief Evaluate a degree-5 polynomial using Horner's method.
 *
 * @param x    The polynomial argument.
 * @param c0   The constant coefficient.
 * @param c1   The coefficient of x^1.
 * @param c2   The coefficient of x^2.
 * @param c3   The coefficient of x^3.
 * @param c4   The coefficient of x^4.
 * @param c5   The coefficient of x^5.
 */
#define POLY5(x, c0, c1, c2, c3, c4, c5) ((POLY4(x, c1, c2, c3, c4, c5) * x) + c0)

/**
 * @brief Compute an approximate exp2(x) for each lane in the vector.
 *
 * Based on 5th degree minimax polynomials, ported from this blog
 * https://jrfonseca.blogspot.com/2008/09/fast-sse2-pow-tables-or-polynomials.html
 *
 * @param x   The exponents.
 */
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

/**
 * @brief Compute an approximate log2(x) for each lane in the vector.
 *
 * Based on 5th degree minimax polynomials, ported from this blog
 * https://jrfonseca.blogspot.com/2008/09/fast-sse2-pow-tables-or-polynomials.html
 *
 * @param x   The values whose logarithms are required.
 */
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

/**
 * @brief Compute an approximate pow(x, y) for each lane in the vector.
 *
 * Power function based on the exp2(log2(x) * y) transform.
 *
 * @param x   The bases.
 * @param y   The exponents.
 */
static ASTCENC_SIMD_INLINE vfloat4 pow(vfloat4 x, vfloat4 y)
{
	vmask4 zero_mask = y == vfloat4(0.0f);
	vfloat4 estimate = exp2(log2(x) * y);

	// Guarantee that y == 0 returns exactly 1.0f
	return select(estimate, vfloat4(1.0f), zero_mask);
}

/**
 * @brief Count the leading zeros for each lane in @c a.
 *
 * Valid for all data values of @c a; will return a per-lane value [0, 32].
 *
 * @param a   The values whose leading zeros are counted.
 */
static ASTCENC_SIMD_INLINE vint4 clz(vint4 a)
{
	// Derive the highest set bit from the integer's floating-point exponent

	// Keep the top eight bits and clear the ninth to prevent rounding up during float conversion
	a = (~lsr<8>(a)) & a;
	a = float_as_int(int_to_float(a));

	// Extract and unbias exponent
	a = vint4(127 + 31) - lsr<23>(a);

	// Clamp result to a valid 32-bit range
	return clamp(0, 32, a);
}

/**
 * @brief Return lanewise 2^a for each lane in @c a.
 *
 * Use of signed int means that this is only valid for values in range [0, 31].
 *
 * @param a   The exponents, in the range [0, 30].
 */
static ASTCENC_SIMD_INLINE vint4 two_to_the_n(vint4 a)
{
	// 2^30 is the largest signed number than can be represented
	assert(all(a < vint4(31)));

	// This function is a horrible abuse of floating point to use the exponent and float conversion
	// to generate a 2^N multiple

	// Bias the exponent
	vint4 exp = a + 127;
	exp = lsl<23>(exp);

	// Reinterpret the bits as a float, and then convert to an int
	vfloat4 f = int_as_float(exp);
	return float_to_int(f);
}

/**
 * @brief Convert unorm16 [0, 65535] to float16 in range [0, 1].
 *
 * @param p   The unsigned normalized 16-bit values.
 */
static ASTCENC_SIMD_INLINE vint4 unorm16_to_sf16(vint4 p)
{
	vint4 fp16_one = vint4(0x3C00);
	vint4 fp16_small = lsl<8>(p);

	vmask4 is_one = p == vint4(0xFFFF);
	vmask4 is_small = p < vint4(4);

	// Manually inline clz() on Visual Studio to avoid release build codegen bug
	// see https://github.com/ARM-software/astc-encoder/issues/259
#if !defined(__clang__) && defined(_MSC_VER)
	vint4 a = (~lsr<8>(p)) & p;
	a = float_as_int(int_to_float(a));
	a = vint4(127 + 31) - lsr<23>(a);
	vint4 lz = clamp(0, 32, a) - 16;
#else
	vint4 lz = clz(p) - 16;
#endif

	p = p * two_to_the_n(lz + 1);
	p = p & vint4(0xFFFF);

	p = lsr<6>(p);

	p = p | lsl<10>(vint4(14) - lz);

	vint4 r = select(p, fp16_one, is_one);
	r = select(r, fp16_small, is_small);
	return r;
}

/**
 * @brief Convert 16-bit LNS to float16.
 *
 * @param p   The 16-bit LNS values.
 */
static ASTCENC_SIMD_INLINE vint4 lns_to_sf16(vint4 p)
{
	vint4 mc = p & 0x7FF;
	vint4 ec = lsr<11>(p);

	vint4 mc_512 = mc * 3;
	vmask4 mask_512 = mc < vint4(512);

	vint4 mc_1536 = mc * 4 - 512;
	vmask4 mask_1536 = mc < vint4(1536);

	vint4 mc_else = mc * 5 - 2048;

	vint4 mt = mc_else;
	mt = select(mt, mc_1536, mask_1536);
	mt = select(mt, mc_512, mask_512);

	vint4 res = lsl<10>(ec) | lsr<3>(mt);
	return min(res, vint4(0x7BFF));
}

/**
 * @brief Extract mantissa and exponent of a float value.
 *
 * @param      a     The input value.
 * @param[out] exp   The output exponent.
 *
 * @return The mantissa.
 */
static ASTCENC_SIMD_INLINE vfloat4 frexp(vfloat4 a, vint4& exp)
{
	// Interpret the bits as an integer
	vint4 ai = float_as_int(a);

	// Extract and unbias the exponent
	exp = (lsr<23>(ai) & 0xFF) - 126;

	// Extract and unbias the mantissa
	vint4 manti = (ai &  static_cast<int>(0x807FFFFF)) | 0x3F000000;
	return int_as_float(manti);
}

/**
 * @brief Convert float to 16-bit LNS.
 *
 * @param a   The floating-point values to convert.
 */
static ASTCENC_SIMD_INLINE vfloat4 float_to_lns(vfloat4 a)
{
	vint4 exp;
	vfloat4 mant = frexp(a, exp);

	// Classify underflow, NaNs, and infinity before transforming the input
	vmask4 mask_underflow_nan = ~(a > vfloat4(1.0f / 67108864.0f));
	vmask4 mask_infinity = a >= vfloat4(65536.0f);

	// If input is smaller than 2^-14, multiply by 2^25 and don't bias
	vmask4 exp_lt_m13 = exp < vint4(-13);

	vfloat4 a1a = a * 33554432.0f;
	vint4 expa = vint4::zero();

	vfloat4 a1b = (mant - 0.5f) * 4096;
	vint4 expb = exp + 14;

	a = select(a1b, a1a, exp_lt_m13);
	exp = select(expb, expa, exp_lt_m13);

	vmask4 a_lt_384 = a < vfloat4(384.0f);
	vmask4 a_lt_1408 = a <= vfloat4(1408.0f);

	vfloat4 a2a = a * (4.0f / 3.0f);
	vfloat4 a2b = a + 128.0f;
	vfloat4 a2c = (a + 512.0f) * (4.0f / 5.0f);

	a = a2c;
	a = select(a, a2b, a_lt_1408);
	a = select(a, a2a, a_lt_384);

	a = a + (int_to_float(exp) * 2048.0f) + 1.0f;

	a = select(a, vfloat4(65535.0f), mask_infinity);
	a = select(a, vfloat4::zero(), mask_underflow_nan);

	return a;
}

namespace astc
{

/**
 * @brief Compute an approximate scalar pow(x, y).
 *
 * @param x   The base.
 * @param y   The exponent.
 */
static ASTCENC_SIMD_INLINE float pow(float x, float y)
{
	return pow(vfloat4(x), vfloat4(y)).lane<0>();
}

}

#endif  // #ifndef ASTC_VECMATHLIB_H_INCLUDED
