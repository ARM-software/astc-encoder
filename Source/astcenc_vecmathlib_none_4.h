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
 * @brief 4x32-bit vectors, implemented using plain C++.
 *
 * This module implements 4-wide 32-bit float, int, and mask vectors. This
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

#ifndef ASTC_VECMATHLIB_NONE_4_H_INCLUDED
#define ASTC_VECMATHLIB_NONE_4_H_INCLUDED

#ifndef ASTCENC_SIMD_INLINE
	#error "Include astcenc_vecmathlib.h, do not include directly"
#endif

#include <cstring>

// ============================================================================
// vfloat4 data type
// ============================================================================

/**
 * @brief Data type for 4-wide floats.
 */
struct vfloat4
{
	/**
	 * @brief Construct from zero-initialized value.
	 */
	ASTCENC_SIMD_INLINE vfloat4() {}

	/**
	 * @brief Construct from 4 values loaded from an unaligned address.
	 *
	 * Consider using loada() which is better with wider VLA vectors if data is
	 * aligned to vector length.
	 */
	ASTCENC_SIMD_INLINE explicit vfloat4(const float* p)
	{
		m[0] = p[0];
		m[1] = p[1];
		m[2] = p[2];
		m[3] = p[3];
	}

	/**
	 * @brief Construct from 4 scalar values replicated across all lanes.
	 *
	 * Consider using zero() for constexpr zeros.
	 */
	ASTCENC_SIMD_INLINE explicit vfloat4(float a)
	{
		m[0] = a;
		m[1] = a;
		m[2] = a;
		m[3] = a;
	}

	/**
	 * @brief Construct from 4 scalar values.
	 *
	 * The value of @c a is stored to lane 0 (LSB) in the SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vfloat4(float a, float b, float c, float d)
	{
		m[0] = a;
		m[1] = b;
		m[2] = c;
		m[3] = d;
	}

	/**
	 * @brief Get the scalar value of a single lane.
	 */
	template <int l> ASTCENC_SIMD_INLINE float lane() const
	{
		return m[l];
	}

	/**
	 * @brief Set the scalar value of a single lane.
	 */
	template <int l> ASTCENC_SIMD_INLINE void set_lane(float a)
	{
		m[l] = a;
	}

	/**
	 * @brief Factory that returns a vector of zeros.
	 */
	static ASTCENC_SIMD_INLINE vfloat4 zero()
	{
		return vfloat4(0.0f);
	}

	/**
	 * @brief Factory that returns a replicated scalar loaded from memory.
	 */
	static ASTCENC_SIMD_INLINE vfloat4 load1(const float* p)
	{
		return vfloat4(*p);
	}

	/**
	 * @brief Factory that returns a vector loaded from aligned memory.
	 */
	static ASTCENC_SIMD_INLINE vfloat4 loada(const float* p)
	{
		return vfloat4(p);
	}

	/**
	 * @brief Factory that returns a vector containing the lane IDs.
	 */
	static ASTCENC_SIMD_INLINE vfloat4 lane_id()
	{
		return vfloat4(0.0f, 1.0f, 2.0f, 3.0f);
	}

	/**
	 * @brief The vector ...
	 */
	float m[4];
};

// ============================================================================
// vint4 data type
// ============================================================================

/**
 * @brief Data type for 4-wide ints.
 */
struct vint4
{
	/**
	 * @brief Construct from zero-initialized value.
	 */
	ASTCENC_SIMD_INLINE vint4() {}

	/**
	 * @brief Construct from 4 values loaded from an unaligned address.
	 *
	 * Consider using vint4::loada() which is better with wider VLA vectors
	 * if data is aligned.
	 */
	ASTCENC_SIMD_INLINE explicit vint4(const int* p)
	{
		m[0] = p[0];
		m[1] = p[1];
		m[2] = p[2];
		m[3] = p[3];
	}

	/**
	 * @brief Construct from 4 scalar values.
	 *
	 * The value of @c a is stored to lane 0 (LSB) in the SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vint4(int a, int b, int c, int d)
	{
		m[0] = a;
		m[1] = b;
		m[2] = c;
		m[3] = d;
	}


	/**
	 * @brief Construct from 4 scalar values replicated across all lanes.
	 *
	 * Consider using vint4::zero() for constexpr zeros.
	 */
	ASTCENC_SIMD_INLINE explicit vint4(int a)
	{
		m[0] = a;
		m[1] = a;
		m[2] = a;
		m[3] = a;
	}

	/**
	 * @brief Get the scalar value of a single lane.
	 */
	template <int l> ASTCENC_SIMD_INLINE int lane() const
	{
		return m[l];
	}

	/**
	 * @brief Factory that returns a vector containing the lane IDs.
	 */
	static ASTCENC_SIMD_INLINE vint4 lane_id()
	{
		return vint4(0, 1, 2, 3);
	}

	/**
	 * @brief The vector ...
	 */
	int m[4];
};

// ============================================================================
// vmask4 data type
// ============================================================================

/**
 * @brief Data type for 4-wide control plane masks.
 */
struct vmask4
{
	/**
	 * @brief Construct from an existing mask value.
	 */
	ASTCENC_SIMD_INLINE explicit vmask4(int* v)
	{
		m[0] = v[0];
		m[1] = v[1];
		m[2] = v[2];
		m[3] = v[3];
	}

	/**
	 * @brief Construct from 4 scalar values.
	 *
	 * The value of @c a is stored to lane 0 (LSB) in the SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vmask4(int a, int b, int c, int d)
	{
		m[0] = a;
		m[1] = b;
		m[2] = c;
		m[3] = d;
	}

	/**
	 * @brief The vector ...
	 */
	int m[4];
};

// ============================================================================
// vmask4 operators and functions
// ============================================================================

/**
 * @brief Overload: mask union (or).
 */
ASTCENC_SIMD_INLINE vmask4 operator|(vmask4 a, vmask4 b)
{
	return vmask4(a.m[0] | b.m[0],
	              a.m[1] | b.m[1],
	              a.m[2] | b.m[2],
	              a.m[3] | b.m[3]);
}

/**
 * @brief Overload: mask intersect (and).
 */
ASTCENC_SIMD_INLINE vmask4 operator&(vmask4 a, vmask4 b)
{
	return vmask4(a.m[0] & b.m[0],
	              a.m[1] & b.m[1],
	              a.m[2] & b.m[2],
	              a.m[3] & b.m[3]);
}

/**
 * @brief Overload: mask difference (xor).
 */
ASTCENC_SIMD_INLINE vmask4 operator^(vmask4 a, vmask4 b)
{
	return vmask4(a.m[0] ^ b.m[0],
	              a.m[1] ^ b.m[1],
	              a.m[2] ^ b.m[2],
	              a.m[3] ^ b.m[3]);
}

/**
 * @brief Overload: mask invert (not).
 */
ASTCENC_SIMD_INLINE vmask4 operator~(vmask4 a)
{
	return vmask4(~a.m[0],
	              ~a.m[1],
	              ~a.m[2],
	              ~a.m[3]);
}

/**
 * @brief Return a 1-bit mask code indicating mask status.
 *
 * bit0 = lane 0
 */
ASTCENC_SIMD_INLINE unsigned int mask(vmask4 a)
{
	return ((a.m[0] >> 31) & 0x1) |
	       ((a.m[1] >> 30) & 0x2) |
	       ((a.m[2] >> 29) & 0x4) |
	       ((a.m[3] >> 28) & 0x8);
}

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
 * @brief Overload: vector by vector addition.
 */
ASTCENC_SIMD_INLINE vint4 operator+(vint4 a, vint4 b)
{
	return vint4(a.m[0] + b.m[0],
	             a.m[1] + b.m[1],
	             a.m[2] + b.m[2],
	             a.m[3] + b.m[3]);
}

/**
 * @brief Overload: vector by vector subtraction.
 */
ASTCENC_SIMD_INLINE vint4 operator-(vint4 a, vint4 b)
{
	return vint4(a.m[0] - b.m[0],
	             a.m[1] - b.m[1],
	             a.m[2] - b.m[2],
	             a.m[3] - b.m[3]);
}

/**
 * @brief Overload: vector bit invert.
 */
ASTCENC_SIMD_INLINE vint4 operator~(vint4 a)
{
	return vint4(~a.m[0],
	             ~a.m[1],
	             ~a.m[2],
	             ~a.m[3]);
}

/**
 * @brief Overload: vector by vector bitwise or.
 */
ASTCENC_SIMD_INLINE vint4 operator|(vint4 a, vint4 b)
{
	return vint4(a.m[0] | b.m[0],
	             a.m[1] | b.m[1],
	             a.m[2] | b.m[2],
	             a.m[3] | b.m[3]);
}

/**
 * @brief Overload: vector by vector bitwise and.
 */
ASTCENC_SIMD_INLINE vint4 operator&(vint4 a, vint4 b)
{
	return vint4(a.m[0] & b.m[0],
	             a.m[1] & b.m[1],
	             a.m[2] & b.m[2],
	             a.m[3] & b.m[3]);
}

/**
 * @brief Overload: vector by vector bitwise xor.
 */
ASTCENC_SIMD_INLINE vint4 operator^(vint4 a, vint4 b)
{
	return vint4(a.m[0] ^ b.m[0],
	             a.m[1] ^ b.m[1],
	             a.m[2] ^ b.m[2],
	             a.m[3] ^ b.m[3]);
}

/**
 * @brief Overload: vector by vector equality.
 */
ASTCENC_SIMD_INLINE vmask4 operator==(vint4 a, vint4 b)
{
	return vmask4(a.m[0] == b.m[0] ? 0xFFFFFFFF : 0,
	              a.m[1] == b.m[1] ? 0xFFFFFFFF : 0,
	              a.m[2] == b.m[2] ? 0xFFFFFFFF : 0,
	              a.m[3] == b.m[3] ? 0xFFFFFFFF : 0);
}

/**
 * @brief Overload: vector by vector inequality.
 */
ASTCENC_SIMD_INLINE vmask4 operator!=(vint4 a, vint4 b)
{
	return vmask4(a.m[0] != b.m[0] ? 0xFFFFFFFF : 0,
	              a.m[1] != b.m[1] ? 0xFFFFFFFF : 0,
	              a.m[2] != b.m[2] ? 0xFFFFFFFF : 0,
	              a.m[3] != b.m[3] ? 0xFFFFFFFF : 0);
}

/**
 * @brief Overload: vector by vector less than.
 */
ASTCENC_SIMD_INLINE vmask4 operator<(vint4 a, vint4 b)
{
	return vmask4(a.m[0] < b.m[0] ? 0xFFFFFFFF : 0,
	              a.m[1] < b.m[1] ? 0xFFFFFFFF : 0,
	              a.m[2] < b.m[2] ? 0xFFFFFFFF : 0,
	              a.m[3] < b.m[3] ? 0xFFFFFFFF : 0);
}

/**
 * @brief Overload: vector by vector greater than.
 */
ASTCENC_SIMD_INLINE vmask4 operator>(vint4 a, vint4 b)
{
	return vmask4(a.m[0] > b.m[0] ? 0xFFFFFFFF : 0,
	              a.m[1] > b.m[1] ? 0xFFFFFFFF : 0,
	              a.m[2] > b.m[2] ? 0xFFFFFFFF : 0,
	              a.m[3] > b.m[3] ? 0xFFFFFFFF : 0);
}

/**
 * @brief Return the min vector of two vectors.
 */
ASTCENC_SIMD_INLINE vint4 min(vint4 a, vint4 b)
{
	return vint4(a.m[0] < b.m[0] ? a.m[0] : b.m[0],
	             a.m[1] < b.m[1] ? a.m[1] : b.m[1],
	             a.m[2] < b.m[2] ? a.m[2] : b.m[2],
	             a.m[3] < b.m[3] ? a.m[3] : b.m[3]);
}

/**
 * @brief Return the min vector of two vectors.
 */
ASTCENC_SIMD_INLINE vint4 max(vint4 a, vint4 b)
{
	return vint4(a.m[0] > b.m[0] ? a.m[0] : b.m[0],
	             a.m[1] > b.m[1] ? a.m[1] : b.m[1],
	             a.m[2] > b.m[2] ? a.m[2] : b.m[2],
	             a.m[3] > b.m[3] ? a.m[3] : b.m[3]);
}

/**
 * @brief Return the horizontal minimum of a single vector.
 */
ASTCENC_SIMD_INLINE vint4 hmin(vint4 v)
{
	int a = std::min(v.m[0], v.m[1]);
	int b = std::min(v.m[2], v.m[3]);
	return vint4(std::min(a, b));
}

/**
 * @brief Store a vector to an aligned memory address.
 */
ASTCENC_SIMD_INLINE void storea(vint4 a, int* p)
{
	p[0] = a.m[0];
	p[1] = a.m[1];
	p[2] = a.m[2];
	p[3] = a.m[3];
}

/**
 * @brief Store lowest N (vector width) bytes into an unaligned address.
 */
ASTCENC_SIMD_INLINE void store_nbytes(vint4 a, uint8_t* p)
{
	int* pi = (int*)p;
	*pi = a.m[0];
}

/**
 * @brief Gather N (vector width) indices from the array.
 */
ASTCENC_SIMD_INLINE vint4 gatheri(const int* base, vint4 indices)
{
	return vint4(base[indices.m[0]],
	             base[indices.m[1]],
	             base[indices.m[2]],
	             base[indices.m[3]]);
}

/**
 * @brief Pack low 8 bits of N (vector width) lanes into bottom of vector.
 */
ASTCENC_SIMD_INLINE vint4 pack_low_bytes(vint4 a)
{
	int b0 = a.m[0] & 0xFF;
	int b1 = a.m[1] & 0xFF;
	int b2 = a.m[2] & 0xFF;
	int b3 = a.m[3] & 0xFF;

	int b = b0 | (b1 << 8) | (b2 << 16) | (b3 << 24);
	return vint4(b, 0, 0, 0);
}

/**
 * @brief Return lanes from @c b if MSB of @c cond is set, else @c a.
 */
ASTCENC_SIMD_INLINE vint4 select(vint4 a, vint4 b, vmask4 cond)
{
	return vint4((cond.m[0] & 0x80000000) ? b.m[0] : a.m[0],
	             (cond.m[1] & 0x80000000) ? b.m[1] : a.m[1],
	             (cond.m[2] & 0x80000000) ? b.m[2] : a.m[2],
	             (cond.m[3] & 0x80000000) ? b.m[3] : a.m[3]);
}

// ============================================================================
// vfloat4 operators and functions
// ============================================================================

/**
 * @brief Overload: vector by vector addition.
 */
ASTCENC_SIMD_INLINE vfloat4 operator+(vfloat4 a, vfloat4 b)
{
	return vfloat4(a.m[0] + b.m[0],
	               a.m[1] + b.m[1],
	               a.m[2] + b.m[2],
	               a.m[3] + b.m[3]);
}

/**
 * @brief Overload: vector by vector subtraction.
 */
ASTCENC_SIMD_INLINE vfloat4 operator-(vfloat4 a, vfloat4 b)
{
	return vfloat4(a.m[0] - b.m[0],
	               a.m[1] - b.m[1],
	               a.m[2] - b.m[2],
	               a.m[3] - b.m[3]);
}

/**
 * @brief Overload: vector by vector multiplication.
 */
ASTCENC_SIMD_INLINE vfloat4 operator*(vfloat4 a, vfloat4 b)
{
	return vfloat4(a.m[0] * b.m[0],
	               a.m[1] * b.m[1],
	               a.m[2] * b.m[2],
	               a.m[3] * b.m[3]);
}

/**
 * @brief Overload: vector by scalar multiplication.
 */
ASTCENC_SIMD_INLINE vfloat4 operator*(vfloat4 a, float b)
{
	return vfloat4(a.m[0] * b,
	               a.m[1] * b,
	               a.m[2] * b,
	               a.m[3] * b);
}

/**
 * @brief Overload: scalar by vector multiplication.
 */
ASTCENC_SIMD_INLINE vfloat4 operator*(float a, vfloat4 b)
{
	return vfloat4(a * b.m[0],
	               a * b.m[1],
	               a * b.m[2],
	               a * b.m[3]);
}

/**
 * @brief Overload: vector by vector division.
 */
ASTCENC_SIMD_INLINE vfloat4 operator/(vfloat4 a, vfloat4 b)
{
	return vfloat4(a.m[0] / b.m[0],
	               a.m[1] / b.m[1],
	               a.m[2] / b.m[2],
	               a.m[3] / b.m[3]);
}

/**
 * @brief Overload: vector by vector equality.
 */
ASTCENC_SIMD_INLINE vmask4 operator==(vfloat4 a, vfloat4 b)
{
	return vmask4(a.m[0] == b.m[0] ? 0xFFFFFFFF : 0,
	              a.m[1] == b.m[1] ? 0xFFFFFFFF : 0,
	              a.m[2] == b.m[2] ? 0xFFFFFFFF : 0,
	              a.m[3] == b.m[3] ? 0xFFFFFFFF : 0);
}

/**
 * @brief Overload: vector by vector inequality.
 */
ASTCENC_SIMD_INLINE vmask4 operator!=(vfloat4 a, vfloat4 b)
{
	return vmask4(a.m[0] != b.m[0] ? 0xFFFFFFFF : 0,
	              a.m[1] != b.m[1] ? 0xFFFFFFFF : 0,
	              a.m[2] != b.m[2] ? 0xFFFFFFFF : 0,
	              a.m[3] != b.m[3] ? 0xFFFFFFFF : 0);
}

/**
 * @brief Overload: vector by vector less than.
 */
ASTCENC_SIMD_INLINE vmask4 operator<(vfloat4 a, vfloat4 b)
{
	return vmask4(a.m[0] < b.m[0] ? 0xFFFFFFFF : 0,
	              a.m[1] < b.m[1] ? 0xFFFFFFFF : 0,
	              a.m[2] < b.m[2] ? 0xFFFFFFFF : 0,
	              a.m[3] < b.m[3] ? 0xFFFFFFFF : 0);
}

/**
 * @brief Overload: vector by vector greater than.
 */
ASTCENC_SIMD_INLINE vmask4 operator>(vfloat4 a, vfloat4 b)
{
	return vmask4(a.m[0] > b.m[0] ? 0xFFFFFFFF : 0,
	              a.m[1] > b.m[1] ? 0xFFFFFFFF : 0,
	              a.m[2] > b.m[2] ? 0xFFFFFFFF : 0,
	              a.m[3] > b.m[3] ? 0xFFFFFFFF : 0);
}

/**
 * @brief Overload: vector by vector less than or equal.
 */
ASTCENC_SIMD_INLINE vmask4 operator<=(vfloat4 a, vfloat4 b)
{
	return vmask4(a.m[0] <= b.m[0] ? 0xFFFFFFFF : 0,
	              a.m[1] <= b.m[1] ? 0xFFFFFFFF : 0,
	              a.m[2] <= b.m[2] ? 0xFFFFFFFF : 0,
	              a.m[3] <= b.m[3] ? 0xFFFFFFFF : 0);
}

/**
 * @brief Overload: vector by vector greater than or equal.
 */
ASTCENC_SIMD_INLINE vmask4 operator>=(vfloat4 a, vfloat4 b)
{
	return vmask4(a.m[0] >= b.m[0] ? 0xFFFFFFFF : 0,
	              a.m[1] >= b.m[1] ? 0xFFFFFFFF : 0,
	              a.m[2] >= b.m[2] ? 0xFFFFFFFF : 0,
	              a.m[3] >= b.m[3] ? 0xFFFFFFFF : 0);
}

/**
 * @brief Return the min vector of two vectors.
 *
 * If either lane value is NaN, @c b will be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat4 min(vfloat4 a, vfloat4 b)
{
	return vfloat4(a.m[0] < b.m[0] ? a.m[0] : b.m[0],
	               a.m[1] < b.m[1] ? a.m[1] : b.m[1],
	               a.m[2] < b.m[2] ? a.m[2] : b.m[2],
	               a.m[3] < b.m[3] ? a.m[3] : b.m[3]);
}

/**
 * @brief Return the max vector of two vectors.
 *
 * If either lane value is NaN, @c b will be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat4 max(vfloat4 a, vfloat4 b)
{
	return vfloat4(a.m[0] > b.m[0] ? a.m[0] : b.m[0],
	               a.m[1] > b.m[1] ? a.m[1] : b.m[1],
	               a.m[2] > b.m[2] ? a.m[2] : b.m[2],
	               a.m[3] > b.m[3] ? a.m[3] : b.m[3]);
}

/**
 * @brief Return the clamped value between min and max.
 *
 * It is assumed that neither @c minv nor @c maxv are NaN values. If @c a is
 * NaN then @c minv will be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat4 clamp(float minv, float maxv, vfloat4 a)
{
	return  min(max(a, vfloat4(minv)), vfloat4(maxv));
}

/**
 * @brief Return a clamped value between 0.0f and max.
 *
 * It is assumed that @c maxv is not a NaN value. If @c a is NaN then zero will
 * be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat4 clampz(float maxv, vfloat4 a)
{
	return min(max(a, vfloat4(0.0f)), vfloat4(maxv));
}

/**
 * @brief Return a clamped value between 0.0f and 1.0f.
 *
 * If @c a is NaN then zero will be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat4 clampzo(vfloat4 a)
{
	return min(max(a, vfloat4(0.0f)), vfloat4(1.0f));
}

/**
 * @brief Return the absolute value of the float vector.
 */
ASTCENC_SIMD_INLINE vfloat4 abs(vfloat4 a)
{
	return vfloat4(std::abs(a.m[0]),
	               std::abs(a.m[1]),
	               std::abs(a.m[2]),
	               std::abs(a.m[3]));
}

/**
 * @brief Return a float rounded to the nearest integer value.
 */
ASTCENC_SIMD_INLINE vfloat4 round(vfloat4 a)
{
	return vfloat4(std::floor(a.m[0] + 0.5f),
	               std::floor(a.m[1] + 0.5f),
	               std::floor(a.m[2] + 0.5f),
	               std::floor(a.m[3] + 0.5f));
}

/**
 * @brief Return the horizontal minimum of a vector.
 */
ASTCENC_SIMD_INLINE vfloat4 hmin(vfloat4 a)
{
	float tmp1 = std::min(a.m[0], a.m[1]);
	float tmp2 = std::min(a.m[2], a.m[3]);
	return vfloat4(std::min(tmp1, tmp2));
}

/**
 * @brief Return the sqrt of the lanes in the vector.
 */
ASTCENC_SIMD_INLINE vfloat4 sqrt(vfloat4 a)
{
	return vfloat4(std::sqrt(a.m[0]),
	               std::sqrt(a.m[1]),
	               std::sqrt(a.m[2]),
	               std::sqrt(a.m[3]));
}

/**
 * @brief Return lanes from @c b if MSB of @c cond is set, else @c a.
 */
ASTCENC_SIMD_INLINE vfloat4 select(vfloat4 a, vfloat4 b, vmask4 cond)
{
	return vfloat4((cond.m[0] & 0x80000000) ? b.m[0] : a.m[0],
	               (cond.m[1] & 0x80000000) ? b.m[1] : a.m[1],
	               (cond.m[2] & 0x80000000) ? b.m[2] : a.m[2],
	               (cond.m[3] & 0x80000000) ? b.m[3] : a.m[3]);
}

/**
 * @brief Load a vector of gathered results from an array;
 */
ASTCENC_SIMD_INLINE vfloat4 gatherf(const float* base, vint4 indices)
{
	return vfloat4(base[indices.m[0]],
	               base[indices.m[1]],
	               base[indices.m[2]],
	               base[indices.m[3]]);
}

/**
 * @brief Store a vector to an unaligned memory address.
 */
ASTCENC_SIMD_INLINE void store(vfloat4 v, float* ptr)
{
	ptr[0] = v.m[0];
	ptr[1] = v.m[1];
	ptr[2] = v.m[2];
	ptr[3] = v.m[3];
}

/**
 * @brief Store a vector to an aligned memory address.
 */
ASTCENC_SIMD_INLINE void storea(vfloat4 v, float* ptr)
{
	ptr[0] = v.m[0];
	ptr[1] = v.m[1];
	ptr[2] = v.m[2];
	ptr[3] = v.m[3];
}

/**
 * @brief Return the dot product for the full 4 lanes, returning vector.
 */
ASTCENC_SIMD_INLINE vfloat4 dot(vfloat4 a, vfloat4 b)
{
	float s = a.m[0] * b.m[0] +
	          a.m[1] * b.m[1] +
	          a.m[2] * b.m[2] +
	          a.m[3] * b.m[3];

	return vfloat4(s);
}

/**
 * @brief Return a integer value for a float vector, using truncation.
 */
ASTCENC_SIMD_INLINE vint4 float_to_int(vfloat4 v)
{
	return vint4(v.m[0],
	             v.m[1],
	             v.m[2],
	             v.m[3]);
}

/**
 * @brief Return a integer value for a float vector, using truncation.
 */
ASTCENC_SIMD_INLINE vint4 float_to_int_rtn(vfloat4 v)
{
	return vint4((int)(v.m[0] + 0.5f),
	             (int)(v.m[1] + 0.5f),
	             (int)(v.m[2] + 0.5f),
	             (int)(v.m[3] + 0.5f));
}

/**
 * @brief Return a float value as an integer bit pattern (i.e. no conversion).
 *
 * It is a common trick to convert floats into integer bit patterns, perform
 * some bit hackery based on knowledge they are IEEE 754 layout, and then
 * convert them back again. This is the first half of that flip.
 */
ASTCENC_SIMD_INLINE vint4 float_as_int(vfloat4 v)
{
	vint4 r;
	memcpy(r.m, v.m, 4 * 4);
	return r;
}

/**
 * @brief Return a integer value as a float bit pattern (i.e. no conversion).
 *
 * It is a common trick to convert floats into integer bit patterns, perform
 * some bit hackery based on knowledge they are IEEE 754 layout, and then
 * convert them back again. This is the second half of that flip.
 */
ASTCENC_SIMD_INLINE vfloat4 int_as_float(vint4 v)
{
	vfloat4 r;
	memcpy(r.m, v.m, 4 * 4);
	return r;
}

#endif // #ifndef ASTC_VECMATHLIB_NONE_1_H_INCLUDED
