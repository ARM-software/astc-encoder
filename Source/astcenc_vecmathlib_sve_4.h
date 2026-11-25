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
 * @brief 4x32-bit vectors using NEON data registers and SVE predicate registers.
 *
 * This module implements 4-wide 32-bit float, int, and mask vectors for
 * Arm SVE and SVE2, retaining NEON registers for four-lane arithmetic.
 *
 * There is a baseline level of functionality provided by all vector widths and
 * implementations. This is implemented using identical function signatures,
 * modulo data type, so we can use them as substitutable implementations in VLA
 * code.
 *
 * The 4-wide vectors are also used as a fixed-width type, and significantly
 * extend the functionality above that available to VLA code.
 */

#if !defined(ASTC_VECMATHLIB_SVE_4_H_INCLUDED)
#define ASTC_VECMATHLIB_SVE_4_H_INCLUDED

#if !defined(ASTCENC_SIMD_INLINE)
	#error "Include astcenc_vecmathlib.h, do not include directly"
#endif

#include <cstdio>
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
	 * @brief Construct with default initialization.
	 */
	ASTCENC_SIMD_INLINE vfloat4() = default;

	/**
	 * @brief Construct from 4 values loaded from an unaligned address.
	 *
	 * Consider using loada() which is better with vectors if data is aligned
	 * to vector length.
	 *
	 * @param p   The source values.
	 */
	ASTCENC_SIMD_INLINE explicit vfloat4(const float *p)
	{
		m = vld1q_f32(p);
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
		m = vdupq_n_f32(a);
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
		float v[4] { a, b, c, d };
		m = vld1q_f32(v);
	}

	/**
	 * @brief Construct from an existing SIMD register.
	 *
	 * @param a   The SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vfloat4(float32x4_t a)
	{
		m = a;
	}

	/**
	 * @brief Get the scalar value of a single lane.
	 *
	 * @tparam l   The lane index to read.
	 *
	 * @return The floating-point result.
	 */
	template <int l> ASTCENC_SIMD_INLINE float lane() const
	{
		return vgetq_lane_f32(m, l);
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
		m = vsetq_lane_f32(a, m, l);
	}

	/**
	 * @brief Factory that returns a vector of zeros.
	 *
	 * @return The floating-point vector result.
	 */
	static ASTCENC_SIMD_INLINE vfloat4 zero()
	{
		return vfloat4(0.0f);
	}

	/**
	 * @brief Factory that returns a replicated scalar loaded from memory.
	 *
	 * @param p   The address of the scalar to replicate.
	 *
	 * @return The floating-point vector result.
	 */
	static ASTCENC_SIMD_INLINE vfloat4 load1(const float* p)
	{
		return vfloat4(vld1q_dup_f32(p));
	}

	/**
	 * @brief Factory that returns a vector loaded from 16B aligned memory.
	 *
	 * @param p   The aligned source address.
	 *
	 * @return The floating-point vector result.
	 */
	static ASTCENC_SIMD_INLINE vfloat4 loada(const float* p)
	{
		return vfloat4(vld1q_f32(p));
	}

	/**
	 * @brief Return a swizzled float 2.
	 *
	 * @tparam l0   The source lane for output lane 0.
	 * @tparam l1   The source lane for output lane 1.
	 *
	 * @return The floating-point vector result.
	 */
	template <int l0, int l1> ASTCENC_SIMD_INLINE vfloat4 swz() const
	{
		return vfloat4(lane<l0>(), lane<l1>(), 0.0f, 0.0f);
	}

	/**
	 * @brief Return a swizzled float 3.
	 *
	 * @tparam l0   The source lane for output lane 0.
	 * @tparam l1   The source lane for output lane 1.
	 * @tparam l2   The source lane for output lane 2.
	 *
	 * @return The floating-point vector result.
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
	 *
	 * @return The floating-point vector result.
	 */
	template <int l0, int l1, int l2, int l3> ASTCENC_SIMD_INLINE vfloat4 swz() const
	{
		return vfloat4(lane<l0>(), lane<l1>(), lane<l2>(), lane<l3>());
	}

	/** @brief The NEON register holding the four floating-point lanes. */
	float32x4_t m;
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
	 * Consider using loada() which is better with vectors if data is aligned
	 * to vector length.
	 *
	 * @param p   The source values.
	 */
	ASTCENC_SIMD_INLINE explicit vint4(const int *p)
	{
		m = vld1q_s32(p);
	}

	/**
	 * @brief Construct from 4 uint8_t loaded from an unaligned address.
	 *
	 * @param p   The source values.
	 */
	ASTCENC_SIMD_INLINE explicit vint4(const uint8_t *p)
	{
		svint32_t data = svld1ub_s32(svptrue_pat_b32(SV_VL4), p);
		m = svget_neonq(data);
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
		m = vdupq_n_s32(a);
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
		int v[4] { a, b, c, d };
		m = vld1q_s32(v);
	}

	/**
	 * @brief Construct from an existing SIMD register.
	 *
	 * @param a   The SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vint4(int32x4_t a)
	{
		m = a;
	}

	/**
	 * @brief Get the scalar value of a single lane.
	 *
	 * @tparam l   The lane index to read.
	 *
	 * @return The selected integer lane.
	 */
	template <int l> ASTCENC_SIMD_INLINE int lane() const
	{
		return vgetq_lane_s32(m, l);
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
		m = vsetq_lane_s32(a, m, l);
	}

	/**
	 * @brief Factory that returns a vector of zeros.
	 *
	 * @return The integer vector result.
	 */
	static ASTCENC_SIMD_INLINE vint4 zero()
	{
		return vint4(0);
	}

	/**
	 * @brief Factory that returns a replicated scalar loaded from memory.
	 *
	 * @param p   The address of the scalar to replicate.
	 *
	 * @return The integer vector result.
	 */
	static ASTCENC_SIMD_INLINE vint4 load1(const int* p)
	{
		return vint4(*p);
	}

	/**
	 * @brief Factory that returns a vector loaded from unaligned memory.
	 *
	 * @param p   The source bytes.
	 *
	 * @return The integer vector result.
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
	 *
	 * @return The integer vector result.
	 */
	static ASTCENC_SIMD_INLINE vint4 loada(const int* p)
	{
		return vint4(p);
	}

	/**
	 * @brief Factory that returns a vector containing the lane IDs.
	 *
	 * @return The integer vector result.
	 */
	static ASTCENC_SIMD_INLINE vint4 lane_id()
	{
		alignas(16) static const int data[4] { 0, 1, 2, 3 };
		return vint4(vld1q_s32(data));
	}

	/** @brief The NEON register holding the four integer lanes. */
	int32x4_t m;
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
	 * @brief Copy the predicate without accessing unused storage bytes.
	 *
	 * @param a   The source mask.
	 */
	ASTCENC_SIMD_INLINE vmask4(const vmask4& a)
	{
		*reinterpret_cast<svbool_t*>(m) = a.predicate();
	}

	/**
	 * @brief Assign the predicate without accessing unused storage bytes.
	 *
	 * @param a   The source mask.
	 *
	 * @return This mask after assignment.
	 */
	ASTCENC_SIMD_INLINE vmask4& operator=(const vmask4& a)
	{
		*reinterpret_cast<svbool_t*>(m) = a.predicate();
		return *this;
	}

	/**
	 * @brief Construct from an SVE predicate, using only the first four 32-bit lanes.
	 *
	 * @param a   The source predicate.
	 */
	ASTCENC_SIMD_INLINE explicit vmask4(svbool_t a)
	{
		*reinterpret_cast<svbool_t*>(m) = a;
	}

	/**
	 * @brief Construct from one scalar boolean replicated across four lanes.
	 *
	 * @param a   The value to replicate.
	 */
	ASTCENC_SIMD_INLINE explicit vmask4(bool a)
	{
		svbool_t value = svpfalse_b();
		if (a)
		{
			value = svptrue_pat_b32(SV_VL4);
		}

		*reinterpret_cast<svbool_t*>(m) = value;
	}

	/**
	 * @brief Construct from four scalar booleans.
	 *
	 * @param a   The value for lane 0.
	 * @param b   The value for lane 1.
	 * @param c   The value for lane 2.
	 * @param d   The value for lane 3.
	 */
	ASTCENC_SIMD_INLINE explicit vmask4(bool a, bool b, bool c, bool d)
	{
		int values[4] { static_cast<int>(a), static_cast<int>(b),
		                static_cast<int>(c), static_cast<int>(d) };
		svbool_t pg = svptrue_pat_b32(SV_VL4);
		*reinterpret_cast<svbool_t*>(m) = svcmpne_n_s32(pg, svld1_s32(pg, values), 0);
	}

	/**
	 * @brief Return the scalable predicate represented by this mask.
	 *
	 * @return The mask predicate; only its first four 32-bit lanes are meaningful.
	 */
	ASTCENC_SIMD_INLINE svbool_t predicate() const
	{
		return *reinterpret_cast<const svbool_t*>(m);
	}

	/**
	 * @brief Get the boolean value of a single lane.
	 *
	 * @tparam l   The lane index to read.
	 *
	 * @return The boolean value of the selected lane.
	 */
	template <int l> ASTCENC_SIMD_INLINE bool lane() const
	{
		static_assert((l >= 0) && (l < 4));
		svbool_t pg = svptrue_pat_b32(SV_VL4);
		svbool_t lane_mask = svcmpeq_n_u32(pg, svindex_u32(0, 1), l);
		return svptest_any(lane_mask, predicate());
	}

	/**
	 * @brief Storage for a predicate at the architectural maximum vector length of 2048 bits.
	 *
	 * Sizeless predicates cannot be class members. Inlined access lets the compiler retain this
	 * value in a predicate register without assuming a fixed hardware vector length.
	 */
	alignas(4) uint8_t m[32];
};

// ============================================================================
// vmask4 operators and functions
// ============================================================================

/**
 * @brief Overload: mask union (or).
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The comparison or logical mask.
 */
ASTCENC_SIMD_INLINE vmask4 operator|(const vmask4& a, const vmask4& b)
{
	return vmask4(svorr_z(svptrue_pat_b32(SV_VL4), a.predicate(), b.predicate()));
}

/**
 * @brief Overload: mask intersect (and).
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The comparison or logical mask.
 */
ASTCENC_SIMD_INLINE vmask4 operator&(const vmask4& a, const vmask4& b)
{
	return vmask4(svand_z(svptrue_pat_b32(SV_VL4), a.predicate(), b.predicate()));
}

/**
 * @brief Overload: mask difference (xor).
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The comparison or logical mask.
 */
ASTCENC_SIMD_INLINE vmask4 operator^(const vmask4& a, const vmask4& b)
{
	return vmask4(sveor_z(svptrue_pat_b32(SV_VL4), a.predicate(), b.predicate()));
}

/**
 * @brief Overload: mask invert (not).
 *
 * @param a   The input mask.
 *
 * @return The comparison or logical mask.
 */
ASTCENC_SIMD_INLINE vmask4 operator~(const vmask4& a)
{
	return vmask4(svnot_z(svptrue_pat_b32(SV_VL4), a.predicate()));
}

/**
 * @brief Return a 4-bit mask code indicating mask status.
 *
 * bit0 = lane 0
 *
 * @param a   The input mask.
 *
 * @return The four-bit lane mask.
 */
ASTCENC_SIMD_INLINE unsigned int mask(const vmask4& a)
{
	svbool_t pg = svptrue_pat_b32(SV_VL4);
	svuint32_t weights = svlsl_u32_x(pg, svdup_u32(1), svindex_u32(0, 1));
	svuint32_t active = svsel_u32(a.predicate(), weights, svdup_u32(0));
	return svaddv_u32(pg, active);
}

#define ASTCENC_USE_NATIVE_VMASK_ANY_ALL 1

/**
 * @brief True if any lanes are enabled, false otherwise.
 *
 * @param a   The input mask.
 *
 * @return The result of the mask test.
 */
ASTCENC_SIMD_INLINE bool any(const vmask4& a)
{
	return svptest_any(svptrue_pat_b32(SV_VL4), a.predicate());
}

/**
 * @brief True if all lanes are enabled, false otherwise.
 *
 * @param a   The input mask.
 *
 * @return The result of the mask test.
 */
ASTCENC_SIMD_INLINE bool all(const vmask4& a)
{
	svbool_t pg = svptrue_pat_b32(SV_VL4);
	return !svptest_any(pg, svnot_z(pg, a.predicate()));
}

// ============================================================================
// vint4 operators and functions
// ============================================================================

#define ASTCENC_HAS_NATIVE_MATCH8

/**
 * @brief Compare eight consecutive bytes against a value.
 *
 * @param p       The unaligned source array of eight bytes.
 * @param value   The byte to match.
 *
 * @return An eight-bit mask, with bit i set when p[i] equals value.
 */
ASTCENC_SIMD_INLINE unsigned int match8(const uint8_t* p, uint8_t value)
{
	static const uint8_t lane_bits[8] { 1, 2, 4, 8, 16, 32, 64, 128 };
	uint8x8_t matches = vceq_u8(vld1_u8(p), vdup_n_u8(value));
	return vaddv_u8(vand_u8(matches, vld1_u8(lane_bits)));
}

/**
 * @brief Overload: vector by vector addition.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 operator+(vint4 a, vint4 b)
{
	return vint4(vaddq_s32(a.m, b.m));
}

/**
 * @brief Overload: vector by vector subtraction.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 operator-(vint4 a, vint4 b)
{
	return vint4(vsubq_s32(a.m, b.m));
}

/**
 * @brief Overload: vector by vector multiplication.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 operator*(vint4 a, vint4 b)
{
	return vint4(vmulq_s32(a.m, b.m));
}

/**
 * @brief Overload: vector bit invert.
 *
 * @param a   The input vector.
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 operator~(vint4 a)
{
	return vint4(vmvnq_s32(a.m));
}

/**
 * @brief Overload: vector by vector bitwise or.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 operator|(vint4 a, vint4 b)
{
	return vint4(vorrq_s32(a.m, b.m));
}

/**
 * @brief Overload: vector by vector bitwise and.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 operator&(vint4 a, vint4 b)
{
	return vint4(vandq_s32(a.m, b.m));
}

/**
 * @brief Overload: vector by vector bitwise xor.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 operator^(vint4 a, vint4 b)
{
	return vint4(veorq_s32(a.m, b.m));
}

/**
 * @brief Overload: vector by vector equality.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The comparison or logical mask.
 */
ASTCENC_SIMD_INLINE vmask4 operator==(vint4 a, vint4 b)
{
	svbool_t pg = svptrue_pat_b32(SV_VL4);
	svint32_t lhs = svset_neonq_s32(svundef_s32(), a.m);
	svint32_t rhs = svset_neonq_s32(svundef_s32(), b.m);
	return vmask4(svcmpeq_s32(pg, lhs, rhs));
}

/**
 * @brief Overload: vector by vector inequality.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The comparison or logical mask.
 */
ASTCENC_SIMD_INLINE vmask4 operator!=(vint4 a, vint4 b)
{
	svbool_t pg = svptrue_pat_b32(SV_VL4);
	svint32_t lhs = svset_neonq_s32(svundef_s32(), a.m);
	svint32_t rhs = svset_neonq_s32(svundef_s32(), b.m);
	return vmask4(svcmpne_s32(pg, lhs, rhs));
}

/**
 * @brief Overload: vector by vector less than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The comparison or logical mask.
 */
ASTCENC_SIMD_INLINE vmask4 operator<(vint4 a, vint4 b)
{
	svbool_t pg = svptrue_pat_b32(SV_VL4);
	svint32_t lhs = svset_neonq_s32(svundef_s32(), a.m);
	svint32_t rhs = svset_neonq_s32(svundef_s32(), b.m);
	return vmask4(svcmplt_s32(pg, lhs, rhs));
}

/**
 * @brief Overload: vector by vector greater than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The comparison or logical mask.
 */
ASTCENC_SIMD_INLINE vmask4 operator>(vint4 a, vint4 b)
{
	svbool_t pg = svptrue_pat_b32(SV_VL4);
	svint32_t lhs = svset_neonq_s32(svundef_s32(), a.m);
	svint32_t rhs = svset_neonq_s32(svundef_s32(), b.m);
	return vmask4(svcmpgt_s32(pg, lhs, rhs));
}

/**
 * @brief Logical shift left.
 *
 * @tparam s   The shift distance, in bits.
 *
 * @param  a   The values to shift.
 *
 * @return The integer vector result.
 */
template <int s> ASTCENC_SIMD_INLINE vint4 lsl(vint4 a)
{
	return vint4(vshlq_s32(a.m, vdupq_n_s32(s)));
}

/**
 * @brief Logical shift right.
 *
 * @tparam s   The shift distance, in bits.
 *
 * @param  a   The values to shift.
 *
 * @return The integer vector result.
 */
template <int s> ASTCENC_SIMD_INLINE vint4 lsr(vint4 a)
{
	uint32x4_t ua = vreinterpretq_u32_s32(a.m);
	ua = vshlq_u32(ua, vdupq_n_s32(-s));
	return vint4(vreinterpretq_s32_u32(ua));
}

/**
 * @brief Arithmetic shift right.
 *
 * @tparam s   The shift distance, in bits.
 *
 * @param  a   The values to shift.
 *
 * @return The integer vector result.
 */
template <int s> ASTCENC_SIMD_INLINE vint4 asr(vint4 a)
{
	return vint4(vshlq_s32(a.m, vdupq_n_s32(-s)));
}

/**
 * @brief Return the min vector of two vectors.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 min(vint4 a, vint4 b)
{
	return vint4(vminq_s32(a.m, b.m));
}

/**
 * @brief Return the max vector of two vectors.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 max(vint4 a, vint4 b)
{
	return vint4(vmaxq_s32(a.m, b.m));
}

/**
 * @brief Return the horizontal minimum of a vector.
 *
 * @param a   The vector to reduce.
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 hmin(vint4 a)
{
	return vint4(vminvq_s32(a.m));
}

/**
 * @brief Return the horizontal maximum of a vector.
 *
 * @param a   The vector to reduce.
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 hmax(vint4 a)
{
	return vint4(vmaxvq_s32(a.m));
}

/**
 * @brief Store a vector to a 16B aligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The aligned destination address.
 */
ASTCENC_SIMD_INLINE void storea(vint4 a, int* p)
{
	vst1q_s32(p, a.m);
}

/**
 * @brief Store a vector to an unaligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void store(vint4 a, int* p)
{
	vst1q_s32(p, a.m);
}

/**
 * @brief Store a vector to an unaligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void store(vint4 a, uint8_t* p)
{
	std::memcpy(p, &a.m, sizeof(int) * 4);
}

/**
 * @brief Store lowest N (vector width) bytes into an unaligned address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void store_nbytes(vint4 a, uint8_t* p)
{
	vst1q_lane_s32(reinterpret_cast<int32_t*>(p), a.m, 0);
}

#define ASTCENC_HAS_NATIVE_STORE_NBYTES_COUNT

/**
 * @brief Store the lowest packed bytes to an unaligned address.
 *
 * @param[in]  a       The vector to store.
 * @param[out] p       The destination address.
 * @param[in]  count   The number of bytes to store, in [0, 16].
 */
ASTCENC_SIMD_INLINE void store_nbytes(vint4 a, uint8_t* p, unsigned int count)
{
	// Predicate only the defined NEON bytes, even when SVE vectors are wider than 128 bits
	uint8x16_t bytes = vreinterpretq_u8_s32(a.m);
	svuint8_t data = svset_neonq_u8(svundef_u8(), bytes);
	svbool_t mask = svwhilelt_b8_u32(0, count);
	svst1_u8(mask, p, data);
}

/**
 * @brief Pack the low 8 bits of each lane and store to an unaligned address.
 *
 * @param      v   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void pack_and_store_low_bytes(vint4 v, uint8_t* p)
{
	// Limit the store to the four NEON lanes, even on wider SVE targets
	svuint32_t av = svset_neonq_u32(svundef_u32(), vreinterpretq_u32_s32(v.m));
	svst1b_u32(svptrue_pat_b32(SV_VL4), p, av);
}

/**
 * @brief Return lanes from @c b if @c cond is set, else @c a.
 *
 * @param a      The vector selected where @c cond is false.
 * @param b      The vector selected where @c cond is true.
 * @param cond   The selection mask.
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 select(vint4 a, vint4 b, const vmask4& cond)
{
	svint32_t lhs = svset_neonq_s32(svundef_s32(), a.m);
	svint32_t rhs = svset_neonq_s32(svundef_s32(), b.m);
	return vint4(svget_neonq_s32(svsel_s32(cond.predicate(), rhs, lhs)));
}

/**
 * @brief Apply optional UNORM8 rounding to UNORM16 color lanes.
 *
 * @param a         The color values, in [0, 65535].
 * @param control   The control from prepare_unorm8_rounding().
 *
 * @return The colors rounded to UNORM8 and rescaled to UNORM16, or unchanged.
 */
ASTCENC_SIMD_INLINE vint4 apply_unorm8_rounding(vint4 a, vint4 control)
{
	uint8x16_t bytes = vreinterpretq_u8_s32(a.m);
	uint8x16_t indices = vreinterpretq_u8_s32(control.m);
	return vint4(vreinterpretq_s32_u8(vqtbl1q_u8(bytes, indices)));
}

// ============================================================================
// vfloat4 operators and functions
// ============================================================================

/**
 * @brief Overload: vector by vector addition.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The floating-point vector result.
 */
ASTCENC_SIMD_INLINE vfloat4 operator+(vfloat4 a, vfloat4 b)
{
	return vfloat4(vaddq_f32(a.m, b.m));
}

/**
 * @brief Overload: vector by vector subtraction.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The floating-point vector result.
 */
ASTCENC_SIMD_INLINE vfloat4 operator-(vfloat4 a, vfloat4 b)
{
	return vfloat4(vsubq_f32(a.m, b.m));
}

/**
 * @brief Overload: vector by vector multiplication.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The floating-point vector result.
 */
ASTCENC_SIMD_INLINE vfloat4 operator*(vfloat4 a, vfloat4 b)
{
	return vfloat4(vmulq_f32(a.m, b.m));
}

/**
 * @brief Overload: vector by vector division.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The floating-point vector result.
 */
ASTCENC_SIMD_INLINE vfloat4 operator/(vfloat4 a, vfloat4 b)
{
	return vfloat4(vdivq_f32(a.m, b.m));
}

/**
 * @brief Overload: vector by vector equality.
 *
 * Returns vector of false mask values if a or b is NaN.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The comparison or logical mask.
 */
ASTCENC_SIMD_INLINE vmask4 operator==(vfloat4 a, vfloat4 b)
{
	svbool_t pg = svptrue_pat_b32(SV_VL4);
	svfloat32_t lhs = svset_neonq_f32(svundef_f32(), a.m);
	svfloat32_t rhs = svset_neonq_f32(svundef_f32(), b.m);
	return vmask4(svcmpeq_f32(pg, lhs, rhs));
}

/**
 * @brief Overload: vector by vector inequality.
 *
 * Returns vector of true mask values if a or b is NaN.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The comparison or logical mask.
 */
ASTCENC_SIMD_INLINE vmask4 operator!=(vfloat4 a, vfloat4 b)
{
	svbool_t pg = svptrue_pat_b32(SV_VL4);
	svfloat32_t lhs = svset_neonq_f32(svundef_f32(), a.m);
	svfloat32_t rhs = svset_neonq_f32(svundef_f32(), b.m);
	return vmask4(svcmpne_f32(pg, lhs, rhs));
}

/**
 * @brief Overload: vector by vector less than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The comparison or logical mask.
 */
ASTCENC_SIMD_INLINE vmask4 operator<(vfloat4 a, vfloat4 b)
{
	svbool_t pg = svptrue_pat_b32(SV_VL4);
	svfloat32_t lhs = svset_neonq_f32(svundef_f32(), a.m);
	svfloat32_t rhs = svset_neonq_f32(svundef_f32(), b.m);
	return vmask4(svcmplt_f32(pg, lhs, rhs));
}

/**
 * @brief Overload: vector by vector greater than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The comparison or logical mask.
 */
ASTCENC_SIMD_INLINE vmask4 operator>(vfloat4 a, vfloat4 b)
{
	svbool_t pg = svptrue_pat_b32(SV_VL4);
	svfloat32_t lhs = svset_neonq_f32(svundef_f32(), a.m);
	svfloat32_t rhs = svset_neonq_f32(svundef_f32(), b.m);
	return vmask4(svcmpgt_f32(pg, lhs, rhs));
}

/**
 * @brief Overload: vector by vector less than or equal.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The comparison or logical mask.
 */
ASTCENC_SIMD_INLINE vmask4 operator<=(vfloat4 a, vfloat4 b)
{
	svbool_t pg = svptrue_pat_b32(SV_VL4);
	svfloat32_t lhs = svset_neonq_f32(svundef_f32(), a.m);
	svfloat32_t rhs = svset_neonq_f32(svundef_f32(), b.m);
	return vmask4(svcmple_f32(pg, lhs, rhs));
}

/**
 * @brief Overload: vector by vector greater than or equal.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 *
 * @return The comparison or logical mask.
 */
ASTCENC_SIMD_INLINE vmask4 operator>=(vfloat4 a, vfloat4 b)
{
	svbool_t pg = svptrue_pat_b32(SV_VL4);
	svfloat32_t lhs = svset_neonq_f32(svundef_f32(), a.m);
	svfloat32_t rhs = svset_neonq_f32(svundef_f32(), b.m);
	return vmask4(svcmpge_f32(pg, lhs, rhs));
}

/**
 * @brief Return the min vector of two vectors.
 *
 * If either lane value is NaN, @c b will be returned for that lane.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 *
 * @return The floating-point vector result.
 */
ASTCENC_SIMD_INLINE vfloat4 min(vfloat4 a, vfloat4 b)
{
	// Do not reorder - second operand will return if either is NaN
	return vfloat4(vminnmq_f32(a.m, b.m));
}

/**
 * @brief Return the max vector of two vectors.
 *
 * If either lane value is NaN, @c b will be returned for that lane.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 *
 * @return The floating-point vector result.
 */
ASTCENC_SIMD_INLINE vfloat4 max(vfloat4 a, vfloat4 b)
{
	// Do not reorder - second operand will return if either is NaN
	return vfloat4(vmaxnmq_f32(a.m, b.m));
}

/**
 * @brief Return the absolute value of the float vector.
 *
 * @param a   The values whose absolute values are required.
 *
 * @return The floating-point vector result.
 */
ASTCENC_SIMD_INLINE vfloat4 abs(vfloat4 a)
{
	float32x4_t zero = vdupq_n_f32(0.0f);
	float32x4_t inv = vsubq_f32(zero, a.m);
	return vfloat4(vmaxq_f32(a.m, inv));
}

/**
 * @brief Return a float rounded to the nearest integer value.
 *
 * @param a   The values to round.
 *
 * @return The floating-point vector result.
 */
ASTCENC_SIMD_INLINE vfloat4 round(vfloat4 a)
{
	return vfloat4(vrndnq_f32(a.m));
}

/**
 * @brief Return the horizontal minimum of a vector.
 *
 * @param a   The vector to reduce.
 *
 * @return The floating-point vector result.
 */
ASTCENC_SIMD_INLINE vfloat4 hmin(vfloat4 a)
{
	return vfloat4(vminvq_f32(a.m));
}

/**
 * @brief Return the horizontal maximum of a vector.
 *
 * @param a   The vector to reduce.
 *
 * @return The floating-point vector result.
 */
ASTCENC_SIMD_INLINE vfloat4 hmax(vfloat4 a)
{
	return vfloat4(vmaxvq_f32(a.m));
}

/**
 * @brief Return the horizontal sum of a vector as a scalar.
 *
 * @param a   The vector to reduce.
 *
 * @return The floating-point result.
 */
ASTCENC_SIMD_INLINE float hadd_s(vfloat4 a)
{
	// Match the x86 reduction order (0 + 2) + (1 + 3); vaddqv uses a different association
	float32x2_t t = vadd_f32(vget_high_f32(a.m), vget_low_f32(a.m));
	return vget_lane_f32(vpadd_f32(t, t), 0);
}

/**
 * @brief Return the sqrt of the lanes in the vector.
 *
 * @param a   The values whose square roots are required.
 *
 * @return The floating-point vector result.
 */
ASTCENC_SIMD_INLINE vfloat4 sqrt(vfloat4 a)
{
	return vfloat4(vsqrtq_f32(a.m));
}

/**
 * @brief Return lanes from @c b if @c cond is set, else @c a.
 *
 * @param a      The vector selected where @c cond is false.
 * @param b      The vector selected where @c cond is true.
 * @param cond   The selection mask.
 *
 * @return The floating-point vector result.
 */
ASTCENC_SIMD_INLINE vfloat4 select(vfloat4 a, vfloat4 b, const vmask4& cond)
{
	svfloat32_t lhs = svset_neonq_f32(svundef_f32(), a.m);
	svfloat32_t rhs = svset_neonq_f32(svundef_f32(), b.m);
	return vfloat4(svget_neonq_f32(svsel_f32(cond.predicate(), rhs, lhs)));
}

/**
 * @brief Load a vector of gathered results from an array.
 *
 * @param base      The source float table.
 * @param indices   The table index for each lane.
 *
 * @return The floating-point vector result.
 */
ASTCENC_SIMD_INLINE vfloat4 gatherf(const float* base, vint4 indices)
{
	svint32_t offsets = svset_neonq_s32(svundef_s32(), indices.m);
	svfloat32_t data = svld1_gather_s32index_f32(svptrue_pat_b32(SV_VL4), base, offsets);
	return vfloat4(svget_neonq_f32(data));
}

/**
 * @brief Load a vector of gathered integers from an array.
 *
 * @param base      The source integer table.
 * @param indices   The table index for each lane.
 *
 * @return The selected table entries.
 */
ASTCENC_SIMD_INLINE vint4 gatheri(const int32_t* base, vint4 indices)
{
	svint32_t offsets = svset_neonq_s32(svundef_s32(), indices.m);
	svint32_t data = svld1_gather_s32index_s32(svptrue_pat_b32(SV_VL4), base, offsets);
	return vint4(svget_neonq_s32(data));
}

#define ASTCENC_USE_NATIVE_GATHERF_BYTE_INDS 1

/**
 * @brief Load a vector of gathered results from an array using byte indices from memory.
 *
 * @param base      The source float table.
 * @param indices   The byte-sized table indices, one per lane.
 *
 * @return The gathered floating-point vector.
 */
template<>
ASTCENC_SIMD_INLINE vfloat4 gatherf_byte_inds<vfloat4>(const float* base, const uint8_t* indices)
{
	svint32_t offsets = svld1ub_s32(svptrue_pat_b32(SV_VL4), indices);
	svfloat32_t data = svld1_gather_s32index_f32(svptrue_pat_b32(SV_VL4), base, offsets);
	return vfloat4(svget_neonq_f32(data));
}

/**
 * @brief Store a vector to an unaligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void store(vfloat4 a, float* p)
{
	vst1q_f32(p, a.m);
}

/**
 * @brief Store a vector to a 16B aligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The aligned destination address.
 */
ASTCENC_SIMD_INLINE void storea(vfloat4 a, float* p)
{
	vst1q_f32(p, a.m);
}

/**
 * @brief Return an integer value for a float vector, using truncation.
 *
 * @param a   The floating-point values to convert.
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 float_to_int(vfloat4 a)
{
	return vint4(vcvtq_s32_f32(a.m));
}

#define ASTCENC_USE_NATIVE_FLOAT_TO_INT_RTN 1

/**
 * @brief Return integers rounded to nearest for nonnegative float values.
 *
 * @param a   The nonnegative floating-point values to convert.
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 float_to_int_rtn(vfloat4 a)
{
	a = a + vfloat4(0.5f);
	return vint4(vcvtq_s32_f32(a.m));
}

/**
 * @brief Return a float value for an integer vector.
 *
 * @param a   The integer values to convert.
 *
 * @return The floating-point vector result.
 */
ASTCENC_SIMD_INLINE vfloat4 int_to_float(vint4 a)
{
	return vfloat4(vcvtq_f32_s32(a.m));
}

/**
 * @brief Return a float16 value for a float vector, using round-to-nearest.
 *
 * @param a   The floating-point values to convert.
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 float_to_float16(vfloat4 a)
{
	// Generate float16 value
	float16x4_t f16 = vcvt_f16_f32(a.m);

	// Convert each 16-bit float pattern to a 32-bit pattern
	uint16x4_t u16 = vreinterpret_u16_f16(f16);
	uint32x4_t u32 = vmovl_u16(u16);
	return vint4(vreinterpretq_s32_u32(u32));
}

/**
 * @brief Return a float16 value for a float scalar, using round-to-nearest.
 *
 * @param a   The floating-point value to convert.
 *
 * @return The half-precision bit pattern.
 */
static inline uint16_t float_to_float16(float a)
{
	vfloat4 av(a);
	return static_cast<uint16_t>(float_to_float16(av).lane<0>());
}

/**
 * @brief Return a float value for a float16 vector.
 *
 * @param a   The half-precision bit patterns to convert.
 *
 * @return The floating-point vector result.
 */
ASTCENC_SIMD_INLINE vfloat4 float16_to_float(vint4 a)
{
	// Convert each 32-bit float pattern to a 16-bit pattern
	uint32x4_t u32 = vreinterpretq_u32_s32(a.m);
	uint16x4_t u16 = vmovn_u32(u32);
	float16x4_t f16 = vreinterpret_f16_u16(u16);

	// Generate float16 value
	return vfloat4(vcvt_f32_f16(f16));
}

/**
 * @brief Return a float value for a float16 scalar.
 *
 * @param a   The half-precision bit pattern to convert.
 *
 * @return The floating-point result.
 */
ASTCENC_SIMD_INLINE float float16_to_float(uint16_t a)
{
	vint4 av(a);
	return float16_to_float(av).lane<0>();
}

/**
 * @brief Return a float value as an integer bit pattern (i.e. no conversion).
 *
 * It is a common trick to convert floats into integer bit patterns, perform
 * some bit hackery based on knowledge they are IEEE 754 layout, and then
 * convert them back again. This is the first half of that flip.
 *
 * @param a   The floating-point vector to reinterpret.
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 float_as_int(vfloat4 a)
{
	return vint4(vreinterpretq_s32_f32(a.m));
}

/**
 * @brief Return an integer value as a float bit pattern (i.e. no conversion).
 *
 * It is a common trick to convert floats into integer bit patterns, perform
 * some bit hackery based on knowledge they are IEEE 754 layout, and then
 * convert them back again. This is the second half of that flip.
 *
 * @param a   The integer vector to reinterpret.
 *
 * @return The floating-point vector result.
 */
ASTCENC_SIMD_INLINE vfloat4 int_as_float(vint4 a)
{
	return vfloat4(vreinterpretq_f32_s32(a.m));
}

/** @brief Register table containing 16 byte entries. */
struct vtable4_16x8 {
	/** @brief The NEON registers holding the table bytes. */
	uint8x16_t t0;
};

/** @brief Register table containing 32 byte entries. */
struct vtable4_32x8 {
	/** @brief The NEON registers holding the table bytes. */
	uint8x16x2_t t01;
};

/** @brief Register table containing 64 byte entries. */
struct vtable4_64x8 {
	/** @brief The NEON registers holding the table bytes. */
	uint8x16x4_t t0123;
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
	table.t0 = vld1q_u8(data);
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
	table.t01 = uint8x16x2_t {
		vld1q_u8(data),
		vld1q_u8(data + 16)
	};
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
	table.t0123 = uint8x16x4_t {
		vld1q_u8(data),
		vld1q_u8(data + 16),
		vld1q_u8(data + 32),
		vld1q_u8(data + 48)
	};
}

/**
 * @brief Perform a vtable lookup in a 16x 8-bit table with 32-bit indices.
 *
 * @param tbl   The prepared lookup table.
 * @param idx   The table index for each lane.
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 vtable_lookup_32bit(
	const vtable4_16x8& tbl,
	vint4 idx
) {
	// Set index byte above max index for unused bytes so table lookup returns zero
	int32x4_t idx_masked = vorrq_s32(idx.m, vdupq_n_s32(0xFFFFFF00));
	uint8x16_t idx_bytes = vreinterpretq_u8_s32(idx_masked);

	return vint4(vreinterpretq_s32_u8(vqtbl1q_u8(tbl.t0, idx_bytes)));
}

/**
 * @brief Perform a vtable lookup in a 32x 8-bit table with 32-bit indices.
 *
 * @param tbl   The prepared lookup table.
 * @param idx   The table index for each lane.
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 vtable_lookup_32bit(
	const vtable4_32x8& tbl,
	vint4 idx
) {
	// Set index byte above max index for unused bytes so table lookup returns zero
	int32x4_t idx_masked = vorrq_s32(idx.m, vdupq_n_s32(0xFFFFFF00));
	uint8x16_t idx_bytes = vreinterpretq_u8_s32(idx_masked);

	return vint4(vreinterpretq_s32_u8(vqtbl2q_u8(tbl.t01, idx_bytes)));
}

/**
 * @brief Perform a vtable lookup in a 64x 8-bit table with 32-bit indices.
 *
 * @param tbl   The prepared lookup table.
 * @param idx   The table index for each lane.
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 vtable_lookup_32bit(
	const vtable4_64x8& tbl,
	vint4 idx
) {
	// Set index byte above max index for unused bytes so table lookup returns zero
	int32x4_t idx_masked = vorrq_s32(idx.m, vdupq_n_s32(0xFFFFFF00));
	uint8x16_t idx_bytes = vreinterpretq_u8_s32(idx_masked);

	return vint4(vreinterpretq_s32_u8(vqtbl4q_u8(tbl.t0123, idx_bytes)));
}

#define ASTCENC_HAS_NATIVE_BIT_EXTRACT

/**
 * @brief Extract a contiguous bit range from a word.
 *
 * @param[in] value    The source word.
 * @param[in] offset   The first bit, in [0, 63].
 * @param[in] count    The number of bits, in [0, 64 - offset].
 *
 * @return The extracted bits, starting at bit zero.
 */
ASTCENC_SIMD_INLINE uint64_t extract_bit_range(
	uint64_t value,
	unsigned int offset,
	unsigned int count
) {
	// NEON shifts by 64 produce zero, defining the mask for both empty and full-width fields
	uint64_t mask = ~vshld_u64(~0ull, static_cast<int64_t>(count));
	return (value >> offset) & mask;
}

#if (ASTCENC_SVE >= 2)

#define ASTCENC_HAS_NATIVE_BIT_PERMUTE

/**
 * @brief Compress selected bits into consecutive low bits, preserving their order.
 *
 * @param[in] value   The source word.
 * @param[in] mask    The source bit positions to select.
 *
 * @return The selected bits, starting at bit zero.
 */
ASTCENC_SIMD_INLINE uint64_t compress_bits(uint64_t value, uint64_t mask)
{
	svuint64_t result = svbext_u64(svdup_u64(value), svdup_u64(mask));
	return svlastb_u64(svptrue_pat_b64(SV_VL1), result);
}

/**
 * @brief Deposit consecutive low bits into selected positions, preserving their order.
 *
 * @param[in] value   The source word.
 * @param[in] mask    The destination bit positions to fill.
 *
 * @return The deposited bits, with all other positions zero.
 */
ASTCENC_SIMD_INLINE uint64_t deposit_bits(uint64_t value, uint64_t mask)
{
	svuint64_t result = svbdep_u64(svdup_u64(value), svdup_u64(mask));
	return svlastb_u64(svptrue_pat_b64(SV_VL1), result);
}

#endif

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
 *
 * @return The integer vector result.
 */
ASTCENC_SIMD_INLINE vint4 interleave_rgba8(vint4 r, vint4 g, vint4 b, vint4 a)
{
	// Insert each channel above the already packed low bytes
	int32x4_t rg = vsliq_n_s32(r.m, g.m, 8);
	int32x4_t rgb = vsliq_n_s32(rg, b.m, 16);
	return vint4(vsliq_n_s32(rgb, a.m, 24));
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
 *
 * @param[out] base   The destination address.
 * @param      data   The vector to store.
 * @param      mask   The lanes to store.
 */
ASTCENC_SIMD_INLINE void store_lanes_masked(uint8_t* base, vint4 data, const vmask4& mask)
{
	svbool_t pg = svptrue_pat_b32(SV_VL4);
	svbool_t active = svand_z(pg, mask.predicate(), pg);
	svint32_t values = svset_neonq_s32(svundef_s32(), data.m);
	svst1_s32(active, reinterpret_cast<int32_t*>(base), values);
}

#endif  // #ifndef ASTC_VECMATHLIB_SVE_4_H_INCLUDED
