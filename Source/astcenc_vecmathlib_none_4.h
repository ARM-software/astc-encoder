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

#if !defined(ASTC_VECMATHLIB_NONE_4_H_INCLUDED)
#define ASTC_VECMATHLIB_NONE_4_H_INCLUDED

#if !defined(ASTCENC_SIMD_INLINE)
	#error "Include astcenc_vecmathlib.h, do not include directly"
#endif

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cfenv>

// ============================================================================
// vfloat4 data type
// ============================================================================

/**
 * @brief Data type for 4-wide floats.
 */
struct vfloat4
{
	/**
	 * @brief Construct with default initialization.
	 */
	ASTCENC_SIMD_INLINE vfloat4() = default;

	/**
	 * @brief Construct from 4 values loaded from an unaligned address.
	 *
	 * Consider using loada() which is better with wider VLA vectors if data is
	 * aligned to vector length.
	 *
	 * @param p   The source values.
	 */
	ASTCENC_SIMD_INLINE explicit vfloat4(const float* p)
	{
		m[0] = p[0];
		m[1] = p[1];
		m[2] = p[2];
		m[3] = p[3];
	}

	/**
	 * @brief Construct from 1 scalar value replicated across all lanes.
	 *
	 * Consider using zero() for constexpr zeros.
	 *
	 * @param a   The value to replicate across all lanes.
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
	 *
	 * @param a   The value for lane 0.
	 * @param b   The value for lane 1.
	 * @param c   The value for lane 2.
	 * @param d   The value for lane 3.
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
	 *
	 * @tparam l   The lane index to read.
	 */
	template <int l> ASTCENC_SIMD_INLINE float lane() const
	{
		return m[l];
	}

	/**
	 * @brief Set the scalar value of a single lane.
	 *
	 * @tparam l   The lane index to update.
	 *
	 * @param  a   The value to store in the selected lane.
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
	 *
	 * @param p   The address of the scalar to replicate.
	 */
	static ASTCENC_SIMD_INLINE vfloat4 load1(const float* p)
	{
		return vfloat4(*p);
	}

	/**
	 * @brief Factory that returns a vector loaded from 16B aligned memory.
	 *
	 * @param p   The aligned source address.
	 */
	static ASTCENC_SIMD_INLINE vfloat4 loada(const float* p)
	{
		return vfloat4(p);
	}

	/**
	 * @brief Return a swizzled float 2.
	 *
	 * @tparam l0   The source lane for output lane 0.
	 * @tparam l1   The source lane for output lane 1.
	 */
	template <int l0, int l1> ASTCENC_SIMD_INLINE vfloat4 swz() const
	{
		return  vfloat4(lane<l0>(), lane<l1>(), 0.0f, 0.0f);
	}

	/**
	 * @brief Return a swizzled float 3.
	 *
	 * @tparam l0   The source lane for output lane 0.
	 * @tparam l1   The source lane for output lane 1.
	 * @tparam l2   The source lane for output lane 2.
	 */
	template <int l0, int l1, int l2> ASTCENC_SIMD_INLINE vfloat4 swz() const
	{
		return vfloat4(lane<l0>(), lane<l1>(), lane<l2>(), 0.0f);
	}

	/**
	 * @brief Return a swizzled float 4.
	 *
	 * @tparam l0   The source lane for output lane 0.
	 * @tparam l1   The source lane for output lane 1.
	 * @tparam l2   The source lane for output lane 2.
	 * @tparam l3   The source lane for output lane 3.
	 */
	template <int l0, int l1, int l2, int l3> ASTCENC_SIMD_INLINE vfloat4 swz() const
	{
		return vfloat4(lane<l0>(), lane<l1>(), lane<l2>(), lane<l3>());
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
	 * @brief Construct with default initialization.
	 */
	ASTCENC_SIMD_INLINE vint4() = default;

	/**
	 * @brief Construct from 4 values loaded from an unaligned address.
	 *
	 * Consider using vint4::loada() which is better with wider VLA vectors
	 * if data is aligned.
	 *
	 * @param p   The source values.
	 */
	ASTCENC_SIMD_INLINE explicit vint4(const int* p)
	{
		m[0] = p[0];
		m[1] = p[1];
		m[2] = p[2];
		m[3] = p[3];
	}

	/**
	 * @brief Construct from 4 uint8_t loaded from an unaligned address.
	 *
	 * @param p   The source values.
	 */
	ASTCENC_SIMD_INLINE explicit vint4(const uint8_t *p)
	{
		m[0] = p[0];
		m[1] = p[1];
		m[2] = p[2];
		m[3] = p[3];
	}

	/**
	 * @brief Construct from 1 scalar value replicated across all lanes.
	 *
	 * Consider using zero() for constexpr zeros.
	 *
	 * @param a   The value to replicate across all lanes.
	 */
	ASTCENC_SIMD_INLINE explicit vint4(int a)
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
	 *
	 * @param a   The value for lane 0.
	 * @param b   The value for lane 1.
	 * @param c   The value for lane 2.
	 * @param d   The value for lane 3.
	 */
	ASTCENC_SIMD_INLINE explicit vint4(int a, int b, int c, int d)
	{
		m[0] = a;
		m[1] = b;
		m[2] = c;
		m[3] = d;
	}

	/**
	 * @brief Get the scalar value of a single lane.
	 *
	 * @tparam l   The lane index to read.
	 */
	template <int l> ASTCENC_SIMD_INLINE int lane() const
	{
		return m[l];
	}

	/**
	 * @brief Set the scalar value of a single lane.
	 *
	 * @tparam l   The lane index to update.
	 *
	 * @param  a   The value to store in the selected lane.
	 */
	template <int l> ASTCENC_SIMD_INLINE void set_lane(int a)
	{
		m[l] = a;
	}

	/**
	 * @brief Factory that returns a vector of zeros.
	 */
	static ASTCENC_SIMD_INLINE vint4 zero()
	{
		return vint4(0);
	}

	/**
	 * @brief Factory that returns a replicated scalar loaded from memory.
	 *
	 * @param p   The address of the scalar to replicate.
	 */
	static ASTCENC_SIMD_INLINE vint4 load1(const int* p)
	{
		return vint4(*p);
	}

	/**
	 * @brief Factory that returns a vector loaded from unaligned memory.
	 *
	 * @param p   The source bytes.
	 */
	static ASTCENC_SIMD_INLINE vint4 load(const uint8_t* p)
	{
		vint4 data;
		std::memcpy(&data.m, p, 4 * sizeof(int));
		return data;
	}

	/**
	 * @brief Factory that returns a vector loaded from 16B aligned memory.
	 *
	 * @param p   The aligned source address.
	 */
	static ASTCENC_SIMD_INLINE vint4 loada(const int* p)
	{
		return vint4(p);
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
	 *
	 * @param p   The source values.
	 */
	ASTCENC_SIMD_INLINE explicit vmask4(int* p)
	{
		m[0] = p[0];
		m[1] = p[1];
		m[2] = p[2];
		m[3] = p[3];
	}

	/**
	 * @brief Construct from 1 scalar value.
	 *
	 * @param a   The value to replicate across all lanes.
	 */
	ASTCENC_SIMD_INLINE explicit vmask4(bool a)
	{
		m[0] = a == false ? 0 : -1;
		m[1] = a == false ? 0 : -1;
		m[2] = a == false ? 0 : -1;
		m[3] = a == false ? 0 : -1;
	}

	/**
	 * @brief Construct from 4 scalar values.
	 *
	 * The value of @c a is stored to lane 0 (LSB) in the SIMD register.
	 *
	 * @param a   The value for lane 0.
	 * @param b   The value for lane 1.
	 * @param c   The value for lane 2.
	 * @param d   The value for lane 3.
	 */
	ASTCENC_SIMD_INLINE explicit vmask4(bool a, bool b, bool c, bool d)
	{
		m[0] = a == false ? 0 : -1;
		m[1] = b == false ? 0 : -1;
		m[2] = c == false ? 0 : -1;
		m[3] = d == false ? 0 : -1;
	}

	/**
	 * @brief Get the scalar value of a single lane.
	 *
	 * @tparam l   The lane index to read.
	 */
	template <int l> ASTCENC_SIMD_INLINE bool lane() const
	{
		return m[l] != 0;
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
 *
 * @param a   The left operand.
 * @param b   The right operand.
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
 *
 * @param a   The left operand.
 * @param b   The right operand.
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
 *
 * @param a   The left operand.
 * @param b   The right operand.
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
 *
 * @param a   The input mask.
 */
ASTCENC_SIMD_INLINE vmask4 operator~(vmask4 a)
{
	return vmask4(~a.m[0],
	              ~a.m[1],
	              ~a.m[2],
	              ~a.m[3]);
}

/**
 * @brief Return a 4-bit mask code indicating mask status.
 *
 * bit0 = lane 0
 *
 * @param a   The input mask.
 */
ASTCENC_SIMD_INLINE unsigned int mask(vmask4 a)
{
	return (a.m[0] & 0x1) |
	       (a.m[1] & 0x2) |
	       (a.m[2] & 0x4) |
	       (a.m[3] & 0x8);
}

#define ASTCENC_USE_NATIVE_VMASK_ANY_ALL 1

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

// ============================================================================
// vint4 operators and functions
// ============================================================================

/**
 * @brief Overload: vector by vector addition.
 *
 * @param a   The left operand.
 * @param b   The right operand.
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
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint4 operator-(vint4 a, vint4 b)
{
	return vint4(a.m[0] - b.m[0],
	             a.m[1] - b.m[1],
	             a.m[2] - b.m[2],
	             a.m[3] - b.m[3]);
}

/**
 * @brief Overload: vector by vector multiplication.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint4 operator*(vint4 a, vint4 b)
{
	return vint4(a.m[0] * b.m[0],
	             a.m[1] * b.m[1],
	             a.m[2] * b.m[2],
	             a.m[3] * b.m[3]);
}

/**
 * @brief Overload: vector bit invert.
 *
 * @param a   The input vector.
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
 *
 * @param a   The left operand.
 * @param b   The right operand.
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
 *
 * @param a   The left operand.
 * @param b   The right operand.
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
 *
 * @param a   The left operand.
 * @param b   The right operand.
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
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask4 operator==(vint4 a, vint4 b)
{
	return vmask4(a.m[0] == b.m[0],
	              a.m[1] == b.m[1],
	              a.m[2] == b.m[2],
	              a.m[3] == b.m[3]);
}

/**
 * @brief Overload: vector by vector inequality.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask4 operator!=(vint4 a, vint4 b)
{
	return vmask4(a.m[0] != b.m[0],
	              a.m[1] != b.m[1],
	              a.m[2] != b.m[2],
	              a.m[3] != b.m[3]);
}

/**
 * @brief Overload: vector by vector less than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask4 operator<(vint4 a, vint4 b)
{
	return vmask4(a.m[0] < b.m[0],
	              a.m[1] < b.m[1],
	              a.m[2] < b.m[2],
	              a.m[3] < b.m[3]);
}

/**
 * @brief Overload: vector by vector greater than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask4 operator>(vint4 a, vint4 b)
{
	return vmask4(a.m[0] > b.m[0],
	              a.m[1] > b.m[1],
	              a.m[2] > b.m[2],
	              a.m[3] > b.m[3]);
}

/**
 * @brief Logical shift left.
 *
 * @tparam s   The shift distance, in bits.
 *
 * @param  a   The values to shift.
 */
template <int s> ASTCENC_SIMD_INLINE vint4 lsl(vint4 a)
{
	// Cast to unsigned to avoid shift in/out of sign bit undefined behavior
	unsigned int as0 = static_cast<unsigned int>(a.m[0]) << s;
	unsigned int as1 = static_cast<unsigned int>(a.m[1]) << s;
	unsigned int as2 = static_cast<unsigned int>(a.m[2]) << s;
	unsigned int as3 = static_cast<unsigned int>(a.m[3]) << s;

	return vint4(static_cast<int>(as0),
	             static_cast<int>(as1),
	             static_cast<int>(as2),
	             static_cast<int>(as3));
}

/**
 * @brief Logical shift right.
 *
 * @tparam s   The shift distance, in bits.
 *
 * @param  a   The values to shift.
 */
template <int s> ASTCENC_SIMD_INLINE vint4 lsr(vint4 a)
{
	// Cast to unsigned to avoid shift in/out of sign bit undefined behavior
	unsigned int as0 = static_cast<unsigned int>(a.m[0]) >> s;
	unsigned int as1 = static_cast<unsigned int>(a.m[1]) >> s;
	unsigned int as2 = static_cast<unsigned int>(a.m[2]) >> s;
	unsigned int as3 = static_cast<unsigned int>(a.m[3]) >> s;

	return vint4(static_cast<int>(as0),
	             static_cast<int>(as1),
	             static_cast<int>(as2),
	             static_cast<int>(as3));
}

/**
 * @brief Arithmetic shift right.
 *
 * @tparam s   The shift distance, in bits.
 *
 * @param  a   The values to shift.
 */
template <int s> ASTCENC_SIMD_INLINE vint4 asr(vint4 a)
{
	return vint4(a.m[0] >> s,
	             a.m[1] >> s,
	             a.m[2] >> s,
	             a.m[3] >> s);
}

/**
 * @brief Return the min vector of two vectors.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vint4 min(vint4 a, vint4 b)
{
	return vint4(a.m[0] < b.m[0] ? a.m[0] : b.m[0],
	             a.m[1] < b.m[1] ? a.m[1] : b.m[1],
	             a.m[2] < b.m[2] ? a.m[2] : b.m[2],
	             a.m[3] < b.m[3] ? a.m[3] : b.m[3]);
}

/**
 * @brief Return the max vector of two vectors.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vint4 max(vint4 a, vint4 b)
{
	return vint4(a.m[0] > b.m[0] ? a.m[0] : b.m[0],
	             a.m[1] > b.m[1] ? a.m[1] : b.m[1],
	             a.m[2] > b.m[2] ? a.m[2] : b.m[2],
	             a.m[3] > b.m[3] ? a.m[3] : b.m[3]);
}

/**
 * @brief Return the horizontal minimum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE vint4 hmin(vint4 a)
{
	int b = std::min(a.m[0], a.m[1]);
	int c = std::min(a.m[2], a.m[3]);
	return vint4(std::min(b, c));
}

/**
 * @brief Return the horizontal maximum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE vint4 hmax(vint4 a)
{
	int b = std::max(a.m[0], a.m[1]);
	int c = std::max(a.m[2], a.m[3]);
	return vint4(std::max(b, c));
}

/**
 * @brief Store a vector to a 16B aligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The aligned destination address.
 */
ASTCENC_SIMD_INLINE void storea(vint4 a, int* p)
{
	p[0] = a.m[0];
	p[1] = a.m[1];
	p[2] = a.m[2];
	p[3] = a.m[3];
}

/**
 * @brief Store a vector to an unaligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void store(vint4 a, int* p)
{
	p[0] = a.m[0];
	p[1] = a.m[1];
	p[2] = a.m[2];
	p[3] = a.m[3];
}

/**
 * @brief Store a vector to an unaligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void store(vint4 a, uint8_t* p)
{
	std::memcpy(p, a.m, sizeof(int) * 4);
}

/**
 * @brief Store lowest N (vector width) bytes into an unaligned address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void store_nbytes(vint4 a, uint8_t* p)
{
	std::memcpy(p, a.m, sizeof(uint8_t) * 4);
}

/**
 * @brief Pack the low 8 bits of each lane and store to an unaligned address.
 *
 * @param      v   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void pack_and_store_low_bytes(vint4 v, uint8_t* p)
{
	int b0 = v.m[0] & 0xFF;
	int b1 = v.m[1] & 0xFF;
	int b2 = v.m[2] & 0xFF;
	int b3 = v.m[3] & 0xFF;

#if !defined(ASTCENC_BIG_ENDIAN)
	int b = b0 | (b1 << 8) | (b2 << 16) | (b3 << 24);
#else
	int b = b3 | (b2 << 8) | (b1 << 16) | (b0 << 24);
#endif
	v = vint4(b, 0, 0, 0);
	store_nbytes(v, p);
}

/**
 * @brief Return lanes from @c b if MSB of @c cond is set, else @c a.
 *
 * @param a      The vector selected where @c cond is false.
 * @param b      The vector selected where @c cond is true.
 * @param cond   The selection mask.
 */
ASTCENC_SIMD_INLINE vint4 select(vint4 a, vint4 b, vmask4 cond)
{
	return vint4((cond.m[0] & static_cast<int>(0x80000000)) ? b.m[0] : a.m[0],
	             (cond.m[1] & static_cast<int>(0x80000000)) ? b.m[1] : a.m[1],
	             (cond.m[2] & static_cast<int>(0x80000000)) ? b.m[2] : a.m[2],
	             (cond.m[3] & static_cast<int>(0x80000000)) ? b.m[3] : a.m[3]);
}

// ============================================================================
// vfloat4 operators and functions
// ============================================================================

/**
 * @brief Overload: vector by vector addition.
 *
 * @param a   The left operand.
 * @param b   The right operand.
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
 *
 * @param a   The left operand.
 * @param b   The right operand.
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
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat4 operator*(vfloat4 a, vfloat4 b)
{
	return vfloat4(a.m[0] * b.m[0],
	               a.m[1] * b.m[1],
	               a.m[2] * b.m[2],
	               a.m[3] * b.m[3]);
}

/**
 * @brief Overload: vector by vector division.
 *
 * @param a   The left operand.
 * @param b   The right operand.
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
 *
 * Returns vector of false mask values if a or b is NaN.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask4 operator==(vfloat4 a, vfloat4 b)
{
	return vmask4(a.m[0] == b.m[0],
	              a.m[1] == b.m[1],
	              a.m[2] == b.m[2],
	              a.m[3] == b.m[3]);
}

/**
 * @brief Overload: vector by vector inequality.
 *
 * Returns vector of true mask values if a or b is NaN.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask4 operator!=(vfloat4 a, vfloat4 b)
{
	return vmask4(a.m[0] != b.m[0],
	              a.m[1] != b.m[1],
	              a.m[2] != b.m[2],
	              a.m[3] != b.m[3]);
}

/**
 * @brief Overload: vector by vector less than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask4 operator<(vfloat4 a, vfloat4 b)
{
	return vmask4(a.m[0] < b.m[0],
	              a.m[1] < b.m[1],
	              a.m[2] < b.m[2],
	              a.m[3] < b.m[3]);
}

/**
 * @brief Overload: vector by vector greater than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask4 operator>(vfloat4 a, vfloat4 b)
{
	return vmask4(a.m[0] > b.m[0],
	              a.m[1] > b.m[1],
	              a.m[2] > b.m[2],
	              a.m[3] > b.m[3]);
}

/**
 * @brief Overload: vector by vector less than or equal.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask4 operator<=(vfloat4 a, vfloat4 b)
{
	return vmask4(a.m[0] <= b.m[0],
	              a.m[1] <= b.m[1],
	              a.m[2] <= b.m[2],
	              a.m[3] <= b.m[3]);
}

/**
 * @brief Overload: vector by vector greater than or equal.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask4 operator>=(vfloat4 a, vfloat4 b)
{
	return vmask4(a.m[0] >= b.m[0],
	              a.m[1] >= b.m[1],
	              a.m[2] >= b.m[2],
	              a.m[3] >= b.m[3]);
}

/**
 * @brief Return the min vector of two vectors.
 *
 * If either lane value is NaN, @c b will be returned for that lane.
 *
 * @param a   The first operand.
 * @param b   The second operand.
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
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vfloat4 max(vfloat4 a, vfloat4 b)
{
	return vfloat4(a.m[0] > b.m[0] ? a.m[0] : b.m[0],
	               a.m[1] > b.m[1] ? a.m[1] : b.m[1],
	               a.m[2] > b.m[2] ? a.m[2] : b.m[2],
	               a.m[3] > b.m[3] ? a.m[3] : b.m[3]);
}

/**
 * @brief Return the absolute value of the float vector.
 *
 * @param a   The values whose absolute values are required.
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
 *
 * @param a   The values to round.
 */
ASTCENC_SIMD_INLINE vfloat4 round(vfloat4 a)
{
	assert(std::fegetround() == FE_TONEAREST);
	return vfloat4(std::nearbyint(a.m[0]),
	               std::nearbyint(a.m[1]),
	               std::nearbyint(a.m[2]),
	               std::nearbyint(a.m[3]));
}

/**
 * @brief Return the horizontal minimum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE vfloat4 hmin(vfloat4 a)
{
	float tmp1 = std::min(a.m[0], a.m[1]);
	float tmp2 = std::min(a.m[2], a.m[3]);
	return vfloat4(std::min(tmp1, tmp2));
}

/**
 * @brief Return the horizontal maximum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE vfloat4 hmax(vfloat4 a)
{
	float tmp1 = std::max(a.m[0], a.m[1]);
	float tmp2 = std::max(a.m[2], a.m[3]);
	return vfloat4(std::max(tmp1, tmp2));
}

/**
 * @brief Return the horizontal sum of a vector as a scalar.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE float hadd_s(vfloat4 a)
{
	// Use halving add, gives invariance with SIMD versions
	return (a.m[0] + a.m[2]) + (a.m[1] + a.m[3]);
}

/**
 * @brief Return the sqrt of the lanes in the vector.
 *
 * @param a   The values whose square roots are required.
 */
ASTCENC_SIMD_INLINE vfloat4 sqrt(vfloat4 a)
{
	return vfloat4(std::sqrt(a.m[0]),
	               std::sqrt(a.m[1]),
	               std::sqrt(a.m[2]),
	               std::sqrt(a.m[3]));
}

/**
 * @brief Return lanes from @c b if @c cond is set, else @c a.
 *
 * @param a      The vector selected where @c cond is false.
 * @param b      The vector selected where @c cond is true.
 * @param cond   The selection mask.
 */
ASTCENC_SIMD_INLINE vfloat4 select(vfloat4 a, vfloat4 b, vmask4 cond)
{
	return vfloat4((cond.m[0] & static_cast<int>(0x80000000)) ? b.m[0] : a.m[0],
	               (cond.m[1] & static_cast<int>(0x80000000)) ? b.m[1] : a.m[1],
	               (cond.m[2] & static_cast<int>(0x80000000)) ? b.m[2] : a.m[2],
	               (cond.m[3] & static_cast<int>(0x80000000)) ? b.m[3] : a.m[3]);
}

/**
 * @brief Gather integers from a table using vector indices.
 *
 * @param base      The source integer table.
 * @param indices   The table index for each lane.
 *
 * @return The selected table entries.
 */
ASTCENC_SIMD_INLINE vint4 gatheri(const int32_t* base, vint4 indices)
{
	return vint4(base[indices.m[0]],
	             base[indices.m[1]],
	             base[indices.m[2]],
	             base[indices.m[3]]);
}

/**
 * @brief Load a vector of gathered results from an array.
 *
 * @param base      The source float table.
 * @param indices   The table index for each lane.
 */
ASTCENC_SIMD_INLINE vfloat4 gatherf(const float* base, vint4 indices)
{
	return vfloat4(base[indices.m[0]],
	               base[indices.m[1]],
	               base[indices.m[2]],
	               base[indices.m[3]]);
}

#define ASTCENC_USE_NATIVE_GATHERF_BYTE_INDS 1

/**
 * @brief Load a vector of gathered results from an array using byte indices from memory.
 *
 * @param base      The source float table.
 * @param indices   The byte-sized table indices, one per lane.
 */
template<>
ASTCENC_SIMD_INLINE vfloat4 gatherf_byte_inds<vfloat4>(const float* base, const uint8_t* indices)
{
	return vfloat4(base[indices[0]],
	               base[indices[1]],
	               base[indices[2]],
	               base[indices[3]]);
}

/**
 * @brief Store a vector to an unaligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void store(vfloat4 a, float* p)
{
	p[0] = a.m[0];
	p[1] = a.m[1];
	p[2] = a.m[2];
	p[3] = a.m[3];
}

/**
 * @brief Store a vector to a 16B aligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The aligned destination address.
 */
ASTCENC_SIMD_INLINE void storea(vfloat4 a, float* p)
{
	p[0] = a.m[0];
	p[1] = a.m[1];
	p[2] = a.m[2];
	p[3] = a.m[3];
}

/**
 * @brief Return an integer value for a float vector, using truncation.
 *
 * @param a   The floating-point values to convert.
 */
ASTCENC_SIMD_INLINE vint4 float_to_int(vfloat4 a)
{
	return vint4(static_cast<int>(a.m[0]),
	             static_cast<int>(a.m[1]),
	             static_cast<int>(a.m[2]),
	             static_cast<int>(a.m[3]));
}

#define ASTCENC_USE_NATIVE_FLOAT_TO_INT_RTN 1

/**
 * @brief Return integers rounded to nearest for nonnegative float values.
 *
 * @param a   The nonnegative floating-point values to convert.
 */
ASTCENC_SIMD_INLINE vint4 float_to_int_rtn(vfloat4 a)
{
	a = a + vfloat4(0.5f);
	return vint4(static_cast<int>(a.m[0]),
	             static_cast<int>(a.m[1]),
	             static_cast<int>(a.m[2]),
	             static_cast<int>(a.m[3]));
}

/**
 * @brief Return a float value for an integer vector.
 *
 * @param a   The integer values to convert.
 */
ASTCENC_SIMD_INLINE vfloat4 int_to_float(vint4 a)
{
	return vfloat4(static_cast<float>(a.m[0]),
	               static_cast<float>(a.m[1]),
	               static_cast<float>(a.m[2]),
	               static_cast<float>(a.m[3]));
}

/**
 * @brief Return a float16 value for a float vector, using round-to-nearest.
 *
 * @param a   The floating-point values to convert.
 */
ASTCENC_SIMD_INLINE vint4 float_to_float16(vfloat4 a)
{
	return vint4(
		float_to_sf16(a.lane<0>()),
		float_to_sf16(a.lane<1>()),
		float_to_sf16(a.lane<2>()),
		float_to_sf16(a.lane<3>()));
}

/**
 * @brief Return a float16 value for a float scalar, using round-to-nearest.
 *
 * @param a   The floating-point value to convert.
 */
static inline uint16_t float_to_float16(float a)
{
	return float_to_sf16(a);
}

/**
 * @brief Return a float value for a float16 vector.
 *
 * @param a   The half-precision bit patterns to convert.
 */
ASTCENC_SIMD_INLINE vfloat4 float16_to_float(vint4 a)
{
	return vfloat4(
		sf16_to_float(static_cast<uint16_t>(a.lane<0>())),
		sf16_to_float(static_cast<uint16_t>(a.lane<1>())),
		sf16_to_float(static_cast<uint16_t>(a.lane<2>())),
		sf16_to_float(static_cast<uint16_t>(a.lane<3>())));
}

/**
 * @brief Return a float value for a float16 scalar.
 *
 * @param a   The half-precision bit pattern to convert.
 */
ASTCENC_SIMD_INLINE float float16_to_float(uint16_t a)
{
	return sf16_to_float(a);
}

/**
 * @brief Return a float value as an integer bit pattern (i.e. no conversion).
 *
 * It is a common trick to convert floats into integer bit patterns, perform
 * some bit hackery based on knowledge they are IEEE 754 layout, and then
 * convert them back again. This is the first half of that flip.
 *
 * @param a   The floating-point vector to reinterpret.
 */
ASTCENC_SIMD_INLINE vint4 float_as_int(vfloat4 a)
{
	return std::bit_cast<vint4>(a);
}

/**
 * @brief Return an integer value as a float bit pattern (i.e. no conversion).
 *
 * It is a common trick to convert floats into integer bit patterns, perform
 * some bit hackery based on knowledge they are IEEE 754 layout, and then
 * convert them back again. This is the second half of that flip.
 *
 * @param a   The integer vector to reinterpret.
 */
ASTCENC_SIMD_INLINE vfloat4 int_as_float(vint4 a)
{
	return std::bit_cast<vfloat4>(a);
}

/*
 * Table structure for a 16x 8-bit entry table.
 */
struct vtable4_16x8 {
	const uint8_t* data;
};

/*
 * Table structure for a 32x 8-bit entry table.
 */
struct vtable4_32x8 {
	const uint8_t* data;
};

/*
 * Table structure for a 64x 8-bit entry table.
 */
struct vtable4_64x8 {
	const uint8_t* data;
};

/**
 * @brief Prepare a vtable lookup table for 16x 8-bit entries.
 *
 * @param[out] table   The lookup table to prepare.
 * @param      data    The source table bytes.
 */
ASTCENC_SIMD_INLINE void vtable_prepare(
	vtable4_16x8& table,
	const uint8_t* data
) {
	table.data = data;
}

/**
 * @brief Prepare a vtable lookup table for 32x 8-bit entries.
 *
 * @param[out] table   The lookup table to prepare.
 * @param      data    The source table bytes.
 */
ASTCENC_SIMD_INLINE void vtable_prepare(
	vtable4_32x8& table,
	const uint8_t* data
) {
	table.data = data;
}

/**
 * @brief Prepare a vtable lookup table for 64x 8-bit entries.
 *
 * @param[out] table   The lookup table to prepare.
 * @param      data    The source table bytes.
 */
ASTCENC_SIMD_INLINE void vtable_prepare(
	vtable4_64x8& table,
	const uint8_t* data
) {
	table.data = data;
}

/**
 * @brief Perform a vtable lookup in a 16x 8-bit table with 32-bit indices.
 *
 * @param tbl   The prepared lookup table.
 * @param idx   The table index for each lane.
 */
ASTCENC_SIMD_INLINE vint4 vtable_lookup_32bit(
	const vtable4_16x8& tbl,
	vint4 idx
) {
	return vint4(tbl.data[idx.lane<0>()],
	             tbl.data[idx.lane<1>()],
	             tbl.data[idx.lane<2>()],
	             tbl.data[idx.lane<3>()]);
}

/**
 * @brief Perform a vtable lookup in a 32x 8-bit table with 32-bit indices.
 *
 * @param tbl   The prepared lookup table.
 * @param idx   The table index for each lane.
 */
ASTCENC_SIMD_INLINE vint4 vtable_lookup_32bit(
	const vtable4_32x8& tbl,
	vint4 idx
) {
	return vint4(tbl.data[idx.lane<0>()],
	             tbl.data[idx.lane<1>()],
	             tbl.data[idx.lane<2>()],
	             tbl.data[idx.lane<3>()]);
}

/**
 * @brief Perform a vtable lookup in a 64x 8-bit table with 32-bit indices.
 *
 * @param tbl   The prepared lookup table.
 * @param idx   The table index for each lane.
 */
ASTCENC_SIMD_INLINE vint4 vtable_lookup_32bit(
	const vtable4_64x8& tbl,
	vint4 idx
) {
	return vint4(tbl.data[idx.lane<0>()],
	             tbl.data[idx.lane<1>()],
	             tbl.data[idx.lane<2>()],
	             tbl.data[idx.lane<3>()]);
}

/**
 * @brief Return a vector of interleaved RGBA data.
 *
 * Input vectors have the value stored in the bottom 8 bits of each lane,
 * with high  bits set to zero.
 *
 * Output vector stores a single RGBA texel packed in each lane.
 *
 * @param r   The red channel values.
 * @param g   The green channel values.
 * @param b   The blue channel values.
 * @param a   The alpha channel values.
 */
ASTCENC_SIMD_INLINE vint4 interleave_rgba8(vint4 r, vint4 g, vint4 b, vint4 a)
{
#if !defined(ASTCENC_BIG_ENDIAN)
	return r + lsl<8>(g) + lsl<16>(b) + lsl<24>(a);
#else
	return a + lsl<8>(b) + lsl<16>(g) + lsl<24>(r);
#endif
}

#define ASTCENC_USE_NATIVE_STORE_LANE 1
#define ASTCENC_USE_NATIVE_STORE_LANES_MASKED 1

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

/**
 * @brief Store a vector, skipping masked lanes.
 *
 * All masked lanes must be at the end of vector, after all non-masked lanes.
 * Input is a byte array of at least 4 bytes per unmasked entry.
 *
 * @param[out] base   The destination address.
 * @param      data   The vector to store.
 * @param      mask   The lanes to store.
 */
ASTCENC_SIMD_INLINE void store_lanes_masked(uint8_t* base, vint4 data, vmask4 mask)
{
	if (mask.m[3])
	{
		store(data, base);
	}
	else if (mask.m[2])
	{
		store_lane(base + 0, data.lane<0>());
		store_lane(base + 4, data.lane<1>());
		store_lane(base + 8, data.lane<2>());
	}
	else if (mask.m[1])
	{
		store_lane(base + 0, data.lane<0>());
		store_lane(base + 4, data.lane<1>());
	}
	else if (mask.m[0])
	{
		store_lane(base + 0, data.lane<0>());
	}
}

#endif  // #ifndef ASTC_VECMATHLIB_NONE_4_H_INCLUDED
