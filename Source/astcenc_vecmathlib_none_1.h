// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2019-2020 Arm Limited
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
 * @brief 1x32-bit vectors, implemented using plain C++.
 *
 * This module implements 1-wide 32-bit float, int, and mask vectors. This
 * module provides a scalar fallback for VLA code, primarily useful for
 * debugging VLA algorithms without the complexity of handling SIMD. Only the
 * baseline level of functionality needed to support VLA is provided.
 *
 * Note that the vector conditional operators implemented by this module are
 * designed to behave like SIMD conditional operators that generate lane masks.
 * Rather than returning 0/1 booleans like normal C++ code they will return
 * 0/-1 to give a full lane-width bitmask.
 *
 * Note that the documentation for this module still talks about "vectors" to
 * help developers think about the implied VLA behavior when writing optimized
 * paths.
 */

#ifndef ASTC_VECMATHLIB_NONE_1_H_INCLUDED
#define ASTC_VECMATHLIB_NONE_1_H_INCLUDED

#ifndef ASTCENC_SIMD_INLINE
	#error "Include astcenc_vecmathlib.h, do not include directly"
#endif

#include <cstring>

// ============================================================================
// vfloat1 data type
// ============================================================================

/**
 * @brief Data type for 1-wide floats.
 */
struct vfloat1
{
	/**
	 * @brief Construct from zero-initialized value.
	 */
	ASTCENC_SIMD_INLINE vfloat1() {}

	/**
	 * @brief Construct from 1 value loaded from an unaligned address.
	 *
	 * Consider using loada() which is better with wider VLA vectors if data is
	 * aligned to vector length.
	 */
	ASTCENC_SIMD_INLINE explicit vfloat1(const float* p)
	{
		m = *p;
	}

	/**
	 * @brief Construct from 1 scalar value replicated across all lanes.
	 *
	 * Consider using zero() for constexpr zeros.
	 */
	ASTCENC_SIMD_INLINE explicit vfloat1(float a)
	{
		m = a;
	}

	/**
	 * @brief Get the scalar value of a single lane.
	 */
	template <int l> ASTCENC_SIMD_INLINE float lane() const
	{
		return m;
	}

	/**
	 * @brief Factory that returns a vector of zeros.
	 */
	static ASTCENC_SIMD_INLINE vfloat1 zero()
	{
		return vfloat1(0.0f);
	}

	/**
	 * @brief Factory that returns a replicated scalar loaded from memory.
	 */
	static ASTCENC_SIMD_INLINE vfloat1 load1(const float* p)
	{
		return vfloat1(*p);
	}

	/**
	 * @brief Factory that returns a vector loaded from aligned memory.
	 */
	static ASTCENC_SIMD_INLINE vfloat1 loada(const float* p)
	{
		return vfloat1(*p);
	}

	/**
	 * @brief Factory that returns a vector containing the lane IDs.
	 */
	static ASTCENC_SIMD_INLINE vfloat1 lane_id()
	{
		return vfloat1(0.0f);
	}

	/**
	 * @brief The vector ...
	 */
	float m;
};

// ============================================================================
// vint1 data type
// ============================================================================

/**
 * @brief Data type for 1-wide ints.
 */
struct vint1
{
	/**
	 * @brief Construct from zero-initialized value.
	 */
	ASTCENC_SIMD_INLINE vint1() {}

	/**
	 * @brief Construct from 1 value loaded from an unaligned address.
	 *
	 * Consider using vint1::loada() which is better with wider VLA vectors
	 * if data is aligned.
	 */
	ASTCENC_SIMD_INLINE explicit vint1(const int* p)
	{
		m = *p;
	}

	/**
	 * @brief Construct from 1 scalar value replicated across all lanes.
	 *
	 * Consider using vint1::zero() for constexpr zeros.
	 */
	ASTCENC_SIMD_INLINE explicit vint1(int a)
	{
		m = a;
	}

	/**
	 * @brief Get the scalar value of a single lane.
	 */
	template <int l> ASTCENC_SIMD_INLINE int lane() const
	{
		return m;
	}

	/**
	 * @brief Factory that returns a vector containing the lane IDs.
	 */
	static ASTCENC_SIMD_INLINE vint1 lane_id()
	{
		return vint1(0);
	}

	/**
	 * @brief The vector ...
	 */
	int m;
};

// ============================================================================
// vmask1 data type
// ============================================================================

/**
 * @brief Data type for 1-wide control plane masks.
 */
struct vmask1
{
	/**
	 * @brief Construct from an existing mask value.
	 */
	ASTCENC_SIMD_INLINE explicit vmask1(int v)
	{
		m = v;
	}

	/**
	 * @brief The vector ...
	 */
	int m;
};

// ============================================================================
// vmask1 operators and functions
// ============================================================================

/**
 * @brief Overload: mask union (or).
 */
ASTCENC_SIMD_INLINE vmask1 operator|(vmask1 a, vmask1 b)
{
	return vmask1(a.m | b.m);
}

/**
 * @brief Overload: mask intersect (and).
 */
ASTCENC_SIMD_INLINE vmask1 operator&(vmask1 a, vmask1 b)
{
	return vmask1(a.m & b.m);
}

/**
 * @brief Overload: mask difference (xor).
 */
ASTCENC_SIMD_INLINE vmask1 operator^(vmask1 a, vmask1 b)
{
	return vmask1(a.m ^ b.m);
}

/**
 * @brief Overload: mask invert (not).
 */
ASTCENC_SIMD_INLINE vmask1 operator~(vmask1 a)
{
	return vmask1(~a.m);
}

/**
 * @brief Return a 1-bit mask code indicating mask status.
 *
 * bit0 = lane 0
 */
ASTCENC_SIMD_INLINE unsigned int mask(vmask1 a)
{
	return (a.m >> 31) & 0x1;
}

/**
 * @brief True if any lanes are enabled, false otherwise.
 */
ASTCENC_SIMD_INLINE bool any(vmask1 a)
{
	return mask(a) != 0;
}

/**
 * @brief True if all lanes are enabled, false otherwise.
 */
ASTCENC_SIMD_INLINE bool all(vmask1 a)
{
	return mask(a) != 0;
}

// ============================================================================
// vint1 operators and functions
// ============================================================================

/**
 * @brief Overload: vector by vector addition.
 */
ASTCENC_SIMD_INLINE vint1 operator+(vint1 a, vint1 b)
{
	return vint1(a.m + b.m);
}

/**
 * @brief Overload: vector by vector subtraction.
 */
ASTCENC_SIMD_INLINE vint1 operator-(vint1 a, vint1 b)
{
	return vint1(a.m - b.m);
}

/**
 * @brief Overload: vector bit invert.
 */
ASTCENC_SIMD_INLINE vint1 operator~(vint1 a)
{
	return vint1(~a.m);
}

/**
 * @brief Overload: vector by vector bitwise or.
 */
ASTCENC_SIMD_INLINE vint1 operator|(vint1 a, vint1 b)
{
	return vint1(a.m | b.m);
}

/**
 * @brief Overload: vector by vector bitwise and.
 */
ASTCENC_SIMD_INLINE vint1 operator&(vint1 a, vint1 b)
{
	return vint1(a.m & b.m);
}

/**
 * @brief Overload: vector by vector bitwise xor.
 */
ASTCENC_SIMD_INLINE vint1 operator^(vint1 a, vint1 b)
{
	return vint1(a.m ^ b.m);
}

/**
 * @brief Overload: vector by vector equality.
 */
ASTCENC_SIMD_INLINE vmask1 operator==(vint1 a, vint1 b)
{
	return vmask1(a.m == b.m ? 0xFFFFFFFF : 0);
}

/**
 * @brief Overload: vector by vector inequality.
 */
ASTCENC_SIMD_INLINE vmask1 operator!=(vint1 a, vint1 b)
{
	return vmask1(a.m != b.m ? 0xFFFFFFFF : 0);
}

/**
 * @brief Overload: vector by vector less than.
 */
ASTCENC_SIMD_INLINE vmask1 operator<(vint1 a, vint1 b)
{
	return vmask1(a.m < b.m ? 0xFFFFFFFF : 0);
}

/**
 * @brief Overload: vector by vector greater than.
 */
ASTCENC_SIMD_INLINE vmask1 operator>(vint1 a, vint1 b)
{
	return vmask1(a.m > b.m ? 0xFFFFFFFF : 0);
}

/**
 * @brief Return the min vector of two vectors.
 */
ASTCENC_SIMD_INLINE vint1 min(vint1 a, vint1 b)
{
	return vint1(a.m < b.m ? a.m : b.m);
}

/**
 * @brief Return the min vector of two vectors.
 */
ASTCENC_SIMD_INLINE vint1 max(vint1 a, vint1 b)
{
	return vint1(a.m > b.m ? a.m : b.m);
}

/**
 * @brief Return the horizontal minimum of a single vector.
 */
ASTCENC_SIMD_INLINE vint1 hmin(vint1 v)
{
	return v;
}

/**
 * @brief Store a vector to an aligned memory address.
 */
ASTCENC_SIMD_INLINE void storea(vint1 a, int* p)
{
	*p = a.m;
}

/**
 * @brief Store lowest N (vector width) bytes into an unaligned address.
 */
ASTCENC_SIMD_INLINE void store_nbytes(vint1 a, uint8_t* p)
{
	*p = (uint8_t)a.m;
}

/**
 * @brief Gather N (vector width) indices from the array.
 */
ASTCENC_SIMD_INLINE vint1 gatheri(const int* base, vint1 indices)
{
	return vint1(base[indices.m]);
}

/**
 * @brief Pack low 8 bits of N (vector width) lanes into bottom of vector.
 */
ASTCENC_SIMD_INLINE vint1 pack_low_bytes(vint1 a)
{
	return a;
}

/**
 * @brief Return lanes from @c b if MSB of @c cond is set, else @c a.
 */
ASTCENC_SIMD_INLINE vint1 select(vint1 a, vint1 b, vmask1 cond)
{
	return vint1((cond.m & 0x80000000) ? b : a);
}

// ============================================================================
// vfloat1 operators and functions
// ============================================================================

/**
 * @brief Overload: vector by vector addition.
 */
ASTCENC_SIMD_INLINE vfloat1 operator+(vfloat1 a, vfloat1 b)
{
	return vfloat1(a.m + b.m);
}

/**
 * @brief Overload: vector by vector subtraction.
 */
ASTCENC_SIMD_INLINE vfloat1 operator-(vfloat1 a, vfloat1 b)
{
	return vfloat1(a.m - b.m);
}

/**
 * @brief Overload: vector by vector multiplication.
 */
ASTCENC_SIMD_INLINE vfloat1 operator*(vfloat1 a, vfloat1 b)
{
	return vfloat1(a.m * b.m);
}

/**
 * @brief Overload: vector by scalar multiplication.
 */
ASTCENC_SIMD_INLINE vfloat1 operator*(vfloat1 a, float b)
{
	return vfloat1(a.m * b);
}

/**
 * @brief Overload: scalar by vector multiplication.
 */
ASTCENC_SIMD_INLINE vfloat1 operator*(float a, vfloat1 b)
{
	return vfloat1(a * b.m);
}

/**
 * @brief Overload: vector by vector division.
 */
ASTCENC_SIMD_INLINE vfloat1 operator/(vfloat1 a, vfloat1 b)
{
	return vfloat1(a.m / b.m);
}

/**
 * @brief Overload: vector by vector equality.
 */
ASTCENC_SIMD_INLINE vmask1 operator==(vfloat1 a, vfloat1 b)
{
	return vmask1(a.m == b.m ? 0xFFFFFFFF : 0);
}

/**
 * @brief Overload: vector by vector inequality.
 */
ASTCENC_SIMD_INLINE vmask1 operator!=(vfloat1 a, vfloat1 b)
{
	return vmask1(a.m != b.m ? 0xFFFFFFFF : 0);
}

/**
 * @brief Overload: vector by vector less than.
 */
ASTCENC_SIMD_INLINE vmask1 operator<(vfloat1 a, vfloat1 b)
{
	return vmask1(a.m < b.m ? 0xFFFFFFFF : 0);
}

/**
 * @brief Overload: vector by vector greater than.
 */
ASTCENC_SIMD_INLINE vmask1 operator>(vfloat1 a, vfloat1 b)
{
	return vmask1(a.m > b.m ? 0xFFFFFFFF : 0);
}

/**
 * @brief Overload: vector by vector less than or equal.
 */
ASTCENC_SIMD_INLINE vmask1 operator<=(vfloat1 a, vfloat1 b)
{
	return vmask1(a.m <= b.m ? 0xFFFFFFFF : 0);
}

/**
 * @brief Overload: vector by vector greater than or equal.
 */
ASTCENC_SIMD_INLINE vmask1 operator>=(vfloat1 a, vfloat1 b)
{
	return vmask1(a.m >= b.m ? 0xFFFFFFFF : 0);
}

/**
 * @brief Return the min vector of two vectors.
 *
 * If either lane value is NaN, @c b will be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat1 min(vfloat1 a, vfloat1 b)
{
	return vfloat1(a.m < b.m ? a.m : b.m);
}

/**
 * @brief Return the max vector of two vectors.
 *
 * If either lane value is NaN, @c b will be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat1 max(vfloat1 a, vfloat1 b)
{
	return vfloat1(a.m > b.m ? a.m : b.m);
}

/**
 * @brief Return the clamped value between min and max.
 *
 * It is assumed that neither @c minv nor @c maxv are NaN values. If @c a is
 * NaN then @c minv will be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat1 clamp(float minv, float maxv, vfloat1 a)
{
	return  min(max(a, vfloat1(minv)), vfloat1(maxv));;
}

/**
 * @brief Return a clamped value between 0.0f and max.
 *
 * It is assumed that @c maxv is not a NaN value. If @c a is NaN then zero will
 * be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat1 clampz(float maxv, vfloat1 a)
{
	return min(max(a, vfloat1(0.0f)), vfloat1(maxv));
}

/**
 * @brief Return a clamped value between 0.0f and 1.0f.
 *
 * If @c a is NaN then zero will be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat1 clampzo(vfloat1 a)
{
	return min(max(a, vfloat1(0.0f)), vfloat1(1.0f));
}

/**
 * @brief Return the absolute value of the float vector.
 */
ASTCENC_SIMD_INLINE vfloat1 abs(vfloat1 a)
{
	return vfloat1(std::abs(a.m));
}

/**
 * @brief Return a float rounded to the nearest integer value.
 */
ASTCENC_SIMD_INLINE vfloat1 round(vfloat1 a)
{
	return vfloat1(std::floor(a.m + 0.5f));
}

/**
 * @brief Return the horizontal minimum of a vector.
 */
ASTCENC_SIMD_INLINE vfloat1 hmin(vfloat1 a)
{
	return a;
}

/**
 * @brief Return lanes from @c b if MSB of @c cond is set, else @c a.
 */
ASTCENC_SIMD_INLINE vfloat1 select(vfloat1 a, vfloat1 b, vmask1 cond)
{
	return vfloat1((cond.m & 0x80000000) ? b : a);
}

/**
 * @brief Load a vector of gathered results from an array;
 */
ASTCENC_SIMD_INLINE vfloat1 gatherf(const float* base, vint1 indices)
{
	return vfloat1(base[indices.m]);
}

/**
 * @brief Store a vector to an aligned memory address.
 */
ASTCENC_SIMD_INLINE void storea(vfloat1 v, float* ptr)
{
	*ptr = v.m;
}

/**
 * @brief Return a integer value for a float vector, using truncation.
 */
ASTCENC_SIMD_INLINE vint1 float_to_int(vfloat1 v)
{
	return vint1(v.m);
}

/**
 * @brief Return a float value as an integer bit pattern (i.e. no conversion).
 *
 * It is a common trick to convert floats into integer bit patterns, perform
 * some bit hackery based on knowledge they are IEEE 754 layout, and then
 * convert them back again. This is the first half of that flip.
 */
ASTCENC_SIMD_INLINE vint1 float_as_int(vfloat1 v)
{
	vint1 r;
	memcpy(&r.m, &v.m, 4);
	return r;
}

/**
 * @brief Return a integer value as a float bit pattern (i.e. no conversion).
 *
 * It is a common trick to convert floats into integer bit patterns, perform
 * some bit hackery based on knowledge they are IEEE 754 layout, and then
 * convert them back again. This is the second half of that flip.
 */
ASTCENC_SIMD_INLINE vfloat1 int_as_float(vint1 v)
{
	vfloat1 r;
	memcpy(&r.m, &v.m, 4);
	return r;
}

#endif // #ifndef ASTC_VECMATHLIB_NONE_1_H_INCLUDED
