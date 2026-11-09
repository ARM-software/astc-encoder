// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2020-2026 Arm Limited
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

#if !defined(ASTC_VECMATHLIB_COMMON_4_H_INCLUDED)
#define ASTC_VECMATHLIB_COMMON_4_H_INCLUDED

#if !defined(ASTCENC_SIMD_INLINE)
	#error "Include astcenc_vecmathlib.h, do not include directly"
#endif

#include <cstdio>
#include <cstring>
#include <limits>

// ============================================================================
// vmask4 operators and functions
// ============================================================================

#if !defined(ASTCENC_USE_NATIVE_VMASK_ANY_ALL)

/**
 * @brief True if any lanes are enabled, false otherwise.
 *
 * @param a   The input mask.
 */
ASTCENC_SIMD_INLINE bool any(vmask4 a)
{
	return mask(a) != 0;
}

/**
 * @brief True if all lanes are enabled, false otherwise.
 *
 * @param a   The input mask.
 */
ASTCENC_SIMD_INLINE bool all(vmask4 a)
{
	return mask(a) == 0xF;
}

#endif

// ============================================================================
// vint4 operators and functions
// ============================================================================

/**
 * @brief Overload: vector by scalar addition.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint4 operator+(vint4 a, int b)
{
	return a + vint4(b);
}

/**
 * @brief Overload: vector by vector incremental addition.
 *
 * @param[in,out] a   The accumulator to update.
 * @param         b   The values to add.
 */
ASTCENC_SIMD_INLINE vint4& operator+=(vint4& a, const vint4& b)
{
	a = a + b;
	return a;
}

/**
 * @brief Overload: vector by scalar subtraction.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint4 operator-(vint4 a, int b)
{
	return a - vint4(b);
}

/**
 * @brief Overload: vector by scalar multiplication.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint4 operator*(vint4 a, int b)
{
	return a * vint4(b);
}

/**
 * @brief Overload: vector by scalar bitwise or.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint4 operator|(vint4 a, int b)
{
	return a | vint4(b);
}

/**
 * @brief Overload: vector by scalar bitwise and.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint4 operator&(vint4 a, int b)
{
	return a & vint4(b);
}

/**
 * @brief Overload: vector by scalar bitwise xor.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint4 operator^(vint4 a, int b)
{
	return a ^ vint4(b);
}

/**
 * @brief Return the clamped value between min and max.
 *
 * @param minv   The inclusive lower bound.
 * @param maxv   The inclusive upper bound.
 * @param a      The values to clamp.
 */
ASTCENC_SIMD_INLINE vint4 clamp(int minv, int maxv, vint4 a)
{
	return min(max(a, vint4(minv)), vint4(maxv));
}

/**
 * @brief Return the horizontal sum of RGB vector lanes as a scalar.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE int hadd_rgb_s(vint4 a)
{
	return a.lane<0>() + a.lane<1>() + a.lane<2>();
}

/**
 * @brief Return the horizontal minimum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE int hmin_s(vint4 a)
{
	return hmin(a).lane<0>();
}

/**
 * @brief Generate a vint4 from a size_t.
 *
 * @param a   The scalar size to replicate across all lanes.
 */
ASTCENC_SIMD_INLINE vint4 vint4_from_size(size_t a)
{
	assert(a <= std::numeric_limits<int>::max());
	return vint4(static_cast<int>(a));
}

/**
 * @brief Return the horizontal maximum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE int hmax_s(vint4 a)
{
	return hmax(a).lane<0>();
}

// ============================================================================
// vfloat4 operators and functions
// ============================================================================

/**
 * @brief Overload: vector by vector incremental addition.
 *
 * @param[in,out] a   The accumulator to update.
 * @param         b   The values to add.
 */
ASTCENC_SIMD_INLINE vfloat4& operator+=(vfloat4& a, const vfloat4& b)
{
	a = a + b;
	return a;
}

/**
 * @brief Overload: vector by scalar addition.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat4 operator+(vfloat4 a, float b)
{
	return a + vfloat4(b);
}

/**
 * @brief Overload: vector by scalar subtraction.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat4 operator-(vfloat4 a, float b)
{
	return a - vfloat4(b);
}

/**
 * @brief Overload: vector by scalar multiplication.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat4 operator*(vfloat4 a, float b)
{
	return a * vfloat4(b);
}

/**
 * @brief Overload: scalar by vector multiplication.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat4 operator*(float a, vfloat4 b)
{
	return vfloat4(a) * b;
}

/**
 * @brief Overload: vector by scalar division.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat4 operator/(vfloat4 a, float b)
{
	return a / vfloat4(b);
}

/**
 * @brief Overload: scalar by vector division.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat4 operator/(float a, vfloat4 b)
{
	return vfloat4(a) / b;
}

/**
 * @brief Return the min vector of a vector and a scalar.
 *
 * If either lane value is NaN, @c b will be returned for that lane.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vfloat4 min(vfloat4 a, float b)
{
	return min(a, vfloat4(b));
}

/**
 * @brief Return the max vector of a vector and a scalar.
 *
 * If either lane value is NaN, @c b will be returned for that lane.
 *
 * @param a   The first operand.
 * @param b   The second operand.
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
 *
 * @param minv   The inclusive lower bound.
 * @param maxv   The inclusive upper bound.
 * @param a      The values to clamp.
 */
ASTCENC_SIMD_INLINE vfloat4 clamp(float minv, float maxv, vfloat4 a)
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
ASTCENC_SIMD_INLINE vfloat4 clampzo(vfloat4 a)
{
	return clamp(0.0f, 1.0f, a);
}

/**
 * @brief Return the horizontal minimum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE float hmin_s(vfloat4 a)
{
	return hmin(a).lane<0>();
}

/**
 * @brief Return the horizontal min of RGB vector lanes as a scalar.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE float hmin_rgb_s(vfloat4 a)
{
	a.set_lane<3>(a.lane<0>());
	return hmin_s(a);
}

/**
 * @brief Return the horizontal maximum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE float hmax_s(vfloat4 a)
{
	return hmax(a).lane<0>();
}

/**
 * @brief Accumulate lane-wise sums for a vector.
 *
 * @param[in,out] accum   The accumulator to update.
 * @param         a       The values to accumulate.
 */
ASTCENC_SIMD_INLINE void haccumulate(vfloat4& accum, vfloat4 a)
{
	accum += a;
}

/**
 * @brief Accumulate masked lane-wise sums for a vector.
 *
 * @param[in,out] accum   The accumulator to update.
 * @param         a       The values to accumulate.
 * @param         m       The lanes to include in the accumulation.
 */
ASTCENC_SIMD_INLINE void haccumulate(vfloat4& accum, vfloat4 a, vmask4 m)
{
	a = select(vfloat4::zero(), a, m);
	haccumulate(accum, a);
}

#if !defined(ASTCENC_USE_NATIVE_GATHERF_BYTE_INDS)

/**
 * @brief Load a vector of gathered results from an array using byte indices from memory.
 *
 * @param base      The source float table.
 * @param indices   The byte-sized table indices, one per lane.
 */
template<>
ASTCENC_SIMD_INLINE vfloat4 gatherf_byte_inds<vfloat4>(const float* base, const uint8_t* indices)
{
	// Experimentally, in this particular use case (byte indices in memory), using 4 separate scalar
	// loads is appreciably faster than using gathers even if they're available, on every x86 uArch
	// tried, so always do the separate loads even when ASTCENC_X86_GATHERS is enabled
	//
	// Tested on:
	//   - Intel Skylake-X, Coffee Lake, Crestmont, Redwood Cove
	//   - AMD Zen 2, Zen 4
	return vfloat4(base[indices[0]], base[indices[1]], base[indices[2]], base[indices[3]]);
}

#endif

#if !defined(ASTCENC_USE_NATIVE_FLOAT_TO_INT_RTN)

/**
 * @brief Return integers rounded to nearest for nonnegative float values.
 *
 * @param a   The nonnegative floating-point values to convert.
 */
ASTCENC_SIMD_INLINE vint4 float_to_int_rtn(vfloat4 a)
{
	return float_to_int(a + vfloat4(0.5f));
}

#endif

/**
 * @brief Return the horizontal sum of RGB vector lanes as a scalar.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE float hadd_rgb_s(vfloat4 a)
{
	return a.lane<0>() + a.lane<1>() + a.lane<2>();
}

#if !defined(ASTCENC_USE_NATIVE_DOT_PRODUCT)

/**
 * @brief Return the dot product for the full 4 lanes, returning scalar.
 *
 * @param a   The first vector.
 * @param b   The second vector.
 */
ASTCENC_SIMD_INLINE float dot_s(vfloat4 a, vfloat4 b)
{
	vfloat4 m = a * b;
	return hadd_s(m);
}

/**
 * @brief Return the dot product for the full 4 lanes, returning vector.
 *
 * @param a   The first vector.
 * @param b   The second vector.
 */
ASTCENC_SIMD_INLINE vfloat4 dot(vfloat4 a, vfloat4 b)
{
	vfloat4 m = a * b;
	return vfloat4(hadd_s(m));
}

/**
 * @brief Return the dot product for the bottom 3 lanes, returning scalar.
 *
 * @param a   The first vector.
 * @param b   The second vector.
 */
ASTCENC_SIMD_INLINE float dot3_s(vfloat4 a, vfloat4 b)
{
	vfloat4 m = a * b;
	return hadd_rgb_s(m);
}

/**
 * @brief Return the dot product for the bottom 3 lanes, returning vector.
 *
 * @param a   The first vector.
 * @param b   The second vector.
 */
ASTCENC_SIMD_INLINE vfloat4 dot3(vfloat4 a, vfloat4 b)
{
	vfloat4 m = a * b;
	float d3 = hadd_rgb_s(m);
	return vfloat4(d3, d3, d3, 0.0f);
}

#endif

/**
 * @brief Apply signed bit transfer.
 *
 * @param[in,out] input0   The first encoded endpoint.
 * @param[in,out] input1   The second encoded endpoint.
 */
static ASTCENC_SIMD_INLINE void bit_transfer_signed(
	vint4& input0,
	vint4& input1
) {
	input1 = lsr<1>(input1) | (input0 & 0x80);
	input0 = lsr<1>(input0) & 0x3F;

	vmask4 mask = (input0 & 0x20) != vint4::zero();
	input0 = select(input0, input0 - 0x40, mask);
}

#if !defined(ASTCENC_USE_NATIVE_STORE_LANE)

/**
 * @brief Store a single vector lane to an unaligned address.
 *
 * @param[out] base   The destination address.
 * @param      data   The lane value to store.
 */
ASTCENC_SIMD_INLINE void store_lane(uint8_t* base, int data)
{
	std::memcpy(base, &data, sizeof(int));
}

#endif

#if !defined(ASTCENC_USE_NATIVE_STORE_LANES_MASKED)

/**
 * @brief Store a vector, skipping masked lanes.
 *
 * All masked lanes must be at the end of vector, after all non-masked lanes.
 *
 * @param[out] base   The destination address.
 * @param      data   The vector to store.
 * @param      mask   The lanes to store.
 */
ASTCENC_SIMD_INLINE void store_lanes_masked(uint8_t* base, vint4 data, vmask4 mask)
{
	// Do not use masked vector stores here. Some implementations do not guarantee fault suppression
	// for masked lanes, which can cause page faults at the end of an image.
	if (mask.lane<3>() != 0.0f)
	{
		store(data, base);
	}
	else if (mask.lane<2>() != 0.0f)
	{
		store_lane(base + 0, data.lane<0>());
		store_lane(base + 4, data.lane<1>());
		store_lane(base + 8, data.lane<2>());
	}
	else if (mask.lane<1>() != 0.0f)
	{
		store_lane(base + 0, data.lane<0>());
		store_lane(base + 4, data.lane<1>());
	}
	else if (mask.lane<0>() != 0.0f)
	{
		store_lane(base + 0, data.lane<0>());
	}
}

#endif

/**
 * @brief Debug function to print a vector of ints.
 *
 * @param a   The vector to print.
 */
ASTCENC_SIMD_INLINE void print(vint4 a)
{
	ASTCENC_ALIGNAS int v[4];
	storea(a, v);
	printf("v4_i32:\n  %8d %8d %8d %8d\n",
	       v[0], v[1], v[2], v[3]);
}

/**
 * @brief Debug function to print a vector of ints.
 *
 * @param a   The vector to print.
 */
ASTCENC_SIMD_INLINE void printx(vint4 a)
{
	ASTCENC_ALIGNAS int v[4];
	storea(a, v);

	unsigned int uv[4];
	std::memcpy(uv, v, sizeof(int) * 4);

	printf("v4_i32:\n  %08x %08x %08x %08x\n",
		uv[0], uv[1], uv[2], uv[3]);
}

/**
 * @brief Debug function to print a vector of floats.
 *
 * @param a   The vector to print.
 */
ASTCENC_SIMD_INLINE void print(vfloat4 a)
{
	ASTCENC_ALIGNAS float v[4];
	storea(a, v);
	printf("v4_f32:\n  %0.4f %0.4f %0.4f %0.4f\n",
	       static_cast<double>(v[0]), static_cast<double>(v[1]),
	       static_cast<double>(v[2]), static_cast<double>(v[3]));
}

/**
 * @brief Debug function to print a vector of masks.
 *
 * @param a   The vector to print.
 */
ASTCENC_SIMD_INLINE void print(vmask4 a)
{
	print(select(vint4(0), vint4(1), a));
}

#endif  // #ifndef ASTC_VECMATHLIB_COMMON_4_H_INCLUDED
