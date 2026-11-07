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

#ifndef ASTC_VECMATHLIB_COMMON_4_H_INCLUDED
#define ASTC_VECMATHLIB_COMMON_4_H_INCLUDED

#ifndef ASTCENC_SIMD_INLINE
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

#endif

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
 * @brief Return the horizontal minimum of a vector.
 */
ASTCENC_SIMD_INLINE int hmin_s(vint4 a)
{
	return hmin(a).lane<0>();
}

/**
 * @brief Generate a vint4 from a size_t.
 */
ASTCENC_SIMD_INLINE vint4 vint4_from_size(size_t a)
{
	assert(a <= std::numeric_limits<int>::max());
	return vint4(static_cast<int>(a));
}

/**
 * @brief Return the horizontal maximum of a vector.
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
	// Preserve operation and operand order to give desired NaN handling
	return min(max(a, minv), maxv);
}

/**
 * @brief Return the clamped value between 0.0f and 1.0f.
 *
 * If @c a is NaN then zero will be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat4 clampzo(vfloat4 a)
{
	return clamp(0.0f, 1.0f, a);
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
 * @brief Accumulate lane-wise sums for a vector.
 */
ASTCENC_SIMD_INLINE void haccumulate(vfloat4& accum, vfloat4 a)
{
	accum += a;
}

/**
 * @brief Accumulate lane-wise sums for a masked vector.
 */
ASTCENC_SIMD_INLINE void haccumulate(vfloat4& accum, vfloat4 a, vmask4 m)
{
	a = select(vfloat4::zero(), a, m);
	haccumulate(accum, a);
}

#if !defined(ASTCENC_USE_NATIVE_GATHERF_BYTE_INDS)

/**
 * @brief Load a vector of gathered results from an array using byte indices from memory.
 */
template<>
ASTCENC_SIMD_INLINE vfloat4 gatherf_byte_inds<vfloat4>(const float* base, const uint8_t* indices)
{
	// Experimentally, in this particular use case (byte indices in memory),
	// using 4 separate scalar loads is appreciably faster than using gathers
	// even if they're available, on every x86 uArch tried, so always do the
	// separate loads even when ASTCENC_X86_GATHERS is enabled.
	//
	// Tested on:
	//   - Intel Skylake-X, Coffee Lake, Crestmont, Redwood Cove
	//   - AMD Zen 2, Zen 4
	return vfloat4(base[indices[0]], base[indices[1]], base[indices[2]], base[indices[3]]);
}

#endif

#if !defined(ASTCENC_USE_NATIVE_FLOAT_TO_INT_RTN)

/**
 * @brief Return an integer value for a float vector, using round-to-nearest.
 */
ASTCENC_SIMD_INLINE vint4 float_to_int_rtn(vfloat4 a)
{
	return float_to_int(a + vfloat4(0.5f));
}

#endif

/**
 * @brief Return the horizontal sum of RGB vector lanes as a scalar.
 */
ASTCENC_SIMD_INLINE float hadd_rgb_s(vfloat4 a)
{
	return a.lane<0>() + a.lane<1>() + a.lane<2>();
}

#if !defined(ASTCENC_USE_NATIVE_DOT_PRODUCT)

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
 * @brief Return the dot product for the bottom 3 lanes, returning vector.
 */
ASTCENC_SIMD_INLINE vfloat4 dot3(vfloat4 a, vfloat4 b)
{
	vfloat4 m = a * b;
	float d3 = hadd_rgb_s(m);
	return vfloat4(d3, d3, d3, 0.0f);
}

#endif

#if !defined(ASTCENC_USE_NATIVE_POPCOUNT)

/**
 * @brief Population bit count.
 *
 * @param v   The value to population count.
 *
 * @return The number of 1 bits.
 */
static inline int popcount(uint64_t v)
{
	uint64_t mask1 = 0x5555555555555555ULL;
	uint64_t mask2 = 0x3333333333333333ULL;
	uint64_t mask3 = 0x0F0F0F0F0F0F0F0FULL;
	v -= (v >> 1) & mask1;
	v = (v & mask2) + ((v >> 2) & mask2);
	v += v >> 4;
	v &= mask3;
	v *= 0x0101010101010101ULL;
	v >>= 56;
	return static_cast<int>(v);
}

#endif

#if !defined(ASTCENC_USE_NATIVE_POPCOUNT_XOR)

/**
 * @brief Count differing bits against four 64-bit values.
 *
 * @param a   The value to compare.
 * @param b   The four values to compare against.
 *
 * @return Per-lane counts of differing bits, in the range 0 to 64.
 */
ASTCENC_SIMD_INLINE vint4 popcount_xor(uint64_t a, const uint64_t b[4])
{
	return vint4(popcount(a ^ b[0]), popcount(a ^ b[1]),
	             popcount(a ^ b[2]), popcount(a ^ b[3]));
}

#endif

/**
 * @brief Apply signed bit transfer.
 *
 * @param input0   The first encoded endpoint.
 * @param input1   The second encoded endpoint.
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
 */
ASTCENC_SIMD_INLINE void store_lanes_masked(uint8_t* base, vint4 data, vmask4 mask)
{
	// Do not use masked vector stores here. Some implementations do not
	// guarantee fault suppression for masked lanes, which can cause page faults
	// at the end of an image.
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
 */
ASTCENC_SIMD_INLINE void print(vmask4 a)
{
	print(select(vint4(0), vint4(1), a));
}

#endif // #ifndef ASTC_VECMATHLIB_COMMON_4_H_INCLUDED
