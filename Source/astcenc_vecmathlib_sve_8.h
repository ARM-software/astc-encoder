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
 * @brief 8x32-bit vectors, implemented using SVE.
 *
 * This module implements 8-wide 32-bit float, int, and mask vectors for Arm
 * SVE.
 *
 * There is a baseline level of functionality provided by all vector widths and
 * implementations. This is implemented using identical function signatures,
 * modulo data type, so we can use them as substitutable implementations in VLA
 * code.
 */

#if !defined(ASTC_VECMATHLIB_SVE_8_H_INCLUDED)
#define ASTC_VECMATHLIB_SVE_8_H_INCLUDED

#if !defined(ASTCENC_SIMD_INLINE)
	#error "Include astcenc_vecmathlib.h, do not include directly"
#endif

#include <cstdio>

typedef svbool_t svbool_8_t __attribute__((arm_sve_vector_bits(256)));
typedef svuint8_t svuint8_8_t __attribute__((arm_sve_vector_bits(256)));
typedef svuint16_t svuint16_8_t __attribute__((arm_sve_vector_bits(256)));
typedef svuint32_t svuint32_8_t __attribute__((arm_sve_vector_bits(256)));
typedef svint32_t svint32_8_t __attribute__((arm_sve_vector_bits(256)));
typedef svfloat32_t svfloat32_8_t __attribute__((arm_sve_vector_bits(256)));

// ============================================================================
// vfloat8 data type
// ============================================================================

/**
 * @brief Data type for 8-wide floats.
 */
struct vfloat8
{
	/**
	 * @brief Construct with default initialization.
	 */
	ASTCENC_SIMD_INLINE vfloat8() = default;

	/**
	 * @brief Construct from 8 values loaded from an unaligned address.
	 *
	 * Consider using loada() which is better with vectors if data is aligned
	 * to vector length.
	 *
	 * @param p   The source values.
	 */
	ASTCENC_SIMD_INLINE explicit vfloat8(const float *p)
	{
		m = svld1_f32(svptrue_b32(), p);
	}

	/**
	 * @brief Construct from 1 scalar value replicated across all lanes.
	 *
	 * Consider using zero() for constexpr zeros.
	 *
	 * @param a   The value to replicate across all lanes.
	 */
	ASTCENC_SIMD_INLINE explicit vfloat8(float a)
	{
		m = svdup_f32(a);
	}

	/**
	 * @brief Construct from an existing SIMD register.
	 *
	 * @param a   The SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vfloat8(svfloat32_8_t a)
	{
		m = a;
	}

	/**
	 * @brief Factory that returns a vector of zeros.
	 */
	static ASTCENC_SIMD_INLINE vfloat8 zero()
	{
		return vfloat8(0.0f);
	}

	/**
	 * @brief Factory that returns a replicated scalar loaded from memory.
	 *
	 * @param p   The address of the scalar to replicate.
	 */
	static ASTCENC_SIMD_INLINE vfloat8 load1(const float* p)
	{
		return vfloat8(*p);
	}

	/**
	 * @brief Factory that returns a vector loaded from 32B aligned memory.
	 *
	 * @param p   The aligned source address.
	 */
	static ASTCENC_SIMD_INLINE vfloat8 loada(const float* p)
	{
		return vfloat8(p);
	}

	/**
	 * @brief The vector ...
	 */
	svfloat32_8_t m;
};

// ============================================================================
// vint8 data type
// ============================================================================

/**
 * @brief Data type for 8-wide ints.
 */
struct vint8
{
	/**
	 * @brief Construct with default initialization.
	 */
	ASTCENC_SIMD_INLINE vint8() = default;

	/**
	 * @brief Construct from 8 values loaded from an unaligned address.
	 *
	 * Consider using loada() which is better with vectors if data is aligned
	 * to vector length.
	 *
	 * @param p   The source values.
	 */
	ASTCENC_SIMD_INLINE explicit vint8(const int *p)
	{
		m = svld1_s32(svptrue_b32(), p);
	}

	/**
	 * @brief Construct from 8 uint8_t loaded from an unaligned address.
	 *
	 * @param p   The source values.
	 */
	ASTCENC_SIMD_INLINE explicit vint8(const uint8_t *p)
	{
		// Load 8-bit values and expand to 32-bits
		m = svld1ub_s32(svptrue_b32(), p);
	}

	/**
	 * @brief Construct from 1 scalar value replicated across all lanes.
	 *
	 * Consider using zero() for constexpr zeros.
	 *
	 * @param a   The value to replicate across all lanes.
	 */
	ASTCENC_SIMD_INLINE explicit vint8(int a)
	{
		m = svdup_s32(a);
	}

	/**
	 * @brief Construct from an existing SIMD register.
	 *
	 * @param a   The SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vint8(svint32_8_t a)
	{
		m = a;
	}

	/**
	 * @brief Factory that returns a vector of zeros.
	 */
	static ASTCENC_SIMD_INLINE vint8 zero()
	{
		return vint8(0);
	}

	/**
	 * @brief Factory that returns a replicated scalar loaded from memory.
	 *
	 * @param p   The address of the scalar to replicate.
	 */
	static ASTCENC_SIMD_INLINE vint8 load1(const int* p)
	{
		return vint8(*p);
	}

	/**
	 * @brief Factory that returns a vector loaded from unaligned memory.
	 *
	 * @param p   The source bytes.
	 */
	static ASTCENC_SIMD_INLINE vint8 load(const uint8_t* p)
	{
		svuint8_8_t data = svld1_u8(svptrue_b8(), p);
		return vint8(svreinterpret_s32_u8(data));
	}

	/**
	 * @brief Factory that returns a vector loaded from 32B aligned memory.
	 *
	 * @param p   The aligned source address.
	 */
	static ASTCENC_SIMD_INLINE vint8 loada(const int* p)
	{
		return vint8(p);
	}

	/**
	 * @brief Factory that returns a vector containing the lane IDs.
	 */
	static ASTCENC_SIMD_INLINE vint8 lane_id()
	{
		return vint8(svindex_s32(0, 1));
	}

	/**
	 * @brief The vector ...
	 */
	svint32_8_t m;
};

// ============================================================================
// vmask8 data type
// ============================================================================

/**
 * @brief Data type for 8-wide control plane masks.
 */
struct vmask8
{
	/**
	 * @brief Construct from an existing SIMD register.
	 *
	 * @param a   The SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vmask8(svbool_8_t a)
	{
		m = a;
	}

	/**
	 * @brief Construct from 1 scalar value.
	 *
	 * @param a   The value to replicate across all lanes.
	 */
	ASTCENC_SIMD_INLINE explicit vmask8(bool a)
	{
		m = svdup_b32(a);
	}

	/**
	 * @brief The vector ...
	 */
	svbool_8_t m;
};

// ============================================================================
// Quad operations (four consecutive lanes)
// ============================================================================

#define ASTCENC_HAS_NATIVE_QUAD_OPS

/**
 * @brief Store lane zero from each quad as two consecutive floats.
 *
 * @param      a   The input vector.
 * @param[out] p   The unaligned output array of two floats.
 */
ASTCENC_SIMD_INLINE void store_lane0_quads(vfloat8 a, float* p)
{
	svfloat32_t compacted = svtbl_f32(a.m, svindex_u32(0, 4));
	svst1_f32(svptrue_pat_b32(SV_VL2), p, compacted);
}

/**
 * @brief Pack two four-lane float vectors into one eight-lane vector.
 *
 * @param groups   The two source groups, in lane order.
 *
 * @return The packed vector with source groups in consecutive quads.
 */
ASTCENC_SIMD_INLINE vfloat8 pack_quads(const vfloat4* groups)
{
	svfloat32_t lo = svset_neonq_f32(svundef_f32(), groups[0].m);
	svfloat32_t hi = svset_neonq_f32(svundef_f32(), groups[1].m);
	return vfloat8(svsplice_f32(svptrue_pat_b32(SV_VL4), lo, hi));
}

/**
 * @brief Pack two four-lane integer vectors into one eight-lane vector.
 *
 * @param groups   The two source groups, in lane order.
 *
 * @return The packed vector with source groups in consecutive quads.
 */
ASTCENC_SIMD_INLINE vint8 pack_quads(const vint4* groups)
{
	svint32_t lo = svset_neonq_s32(svundef_s32(), groups[0].m);
	svint32_t hi = svset_neonq_s32(svundef_s32(), groups[1].m);
	return vint8(svsplice_s32(svptrue_pat_b32(SV_VL4), lo, hi));
}

/**
 * @brief Apply the same swizzle independently within each quad of four lanes.
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
ASTCENC_SIMD_INLINE vfloat8 swizzle_quads(vfloat8 a)
{
	static_assert(l0 >= 0 && l0 < 4 && l1 >= 0 && l1 < 4 &&
	              l2 >= 0 && l2 < 4 && l3 >= 0 && l3 < 4, "Invalid four-lane swizzle");
	svuint32_t groups = svand_n_u32_x(svptrue_b32(), svindex_u32(0, 1), ~3u);
	svuint32_t lanes = svdupq_n_u32(l0, l1, l2, l3);
	return vfloat8(svtbl_f32(a.m, svorr_u32_x(svptrue_b32(), groups, lanes)));
}

// ============================================================================
// vmask8 operators and functions
// ============================================================================

/**
 * @brief Overload: mask union (or).
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator|(vmask8 a, vmask8 b)
{
	return vmask8(svorr_z(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Overload: mask intersect (and).
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator&(vmask8 a, vmask8 b)
{
	return vmask8(svand_z(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Overload: mask difference (xor).
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator^(vmask8 a, vmask8 b)
{
	return vmask8(sveor_z(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Overload: mask invert (not).
 *
 * @param a   The input mask.
 */
ASTCENC_SIMD_INLINE vmask8 operator~(vmask8 a)
{
	return vmask8(svnot_z(svptrue_b32(), a.m));
}

/**
 * @brief Return an 8-bit mask code indicating mask status.
 *
 * bit0 = lane 0
 *
 * @param a   The input mask.
 */
ASTCENC_SIMD_INLINE unsigned int mask(vmask8 a)
{
	alignas(32) const int shifta[8] { 0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80 };
	svint32_8_t template_vals = svld1_s32(svptrue_b32(), shifta);
	svint32_8_t active_vals = svsel_s32(a.m, template_vals, svdup_s32(0));
	return static_cast<unsigned int>(svaddv_s32(svptrue_b32(), active_vals));
}

/**
 * @brief True if any lanes are enabled, false otherwise.
 *
 * @param a   The input mask.
 */
ASTCENC_SIMD_INLINE bool any(vmask8 a)
{
	return svptest_any(svptrue_b32(), a.m);
}

/**
 * @brief True if all lanes are enabled, false otherwise.
 *
 * @param a   The input mask.
 */
ASTCENC_SIMD_INLINE bool all(vmask8 a)
{
	return !svptest_any(svptrue_b32(), (~a).m);
}

// ============================================================================
// vint8 operators and functions
// ============================================================================
/**
 * @brief Overload: vector by vector addition.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint8 operator+(vint8 a, vint8 b)
{
	return vint8(svadd_s32_x(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Overload: vector by vector subtraction.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint8 operator-(vint8 a, vint8 b)
{
	return vint8(svsub_s32_x(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Overload: vector by vector multiplication.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint8 operator*(vint8 a, vint8 b)
{
	return vint8(svmul_s32_x(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Overload: vector bit invert.
 *
 * @param a   The input vector.
 */
ASTCENC_SIMD_INLINE vint8 operator~(vint8 a)
{
	return vint8(svnot_s32_x(svptrue_b32(), a.m));
}

/**
 * @brief Overload: vector by vector bitwise or.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint8 operator|(vint8 a, vint8 b)
{
	return vint8(svorr_s32_x(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Overload: vector by vector bitwise and.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint8 operator&(vint8 a, vint8 b)
{
	return vint8(svand_s32_x(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Overload: vector by vector bitwise xor.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint8 operator^(vint8 a, vint8 b)
{
	return vint8(sveor_s32_x(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Overload: vector by vector equality.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator==(vint8 a, vint8 b)
{
	return vmask8(svcmpeq_s32(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Overload: vector by vector inequality.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator!=(vint8 a, vint8 b)
{
	return vmask8(svcmpne_s32(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Overload: vector by vector less than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator<(vint8 a, vint8 b)
{
	return vmask8(svcmplt_s32(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Overload: vector by vector greater than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator>(vint8 a, vint8 b)
{
	return vmask8(svcmpgt_s32(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Logical shift left.
 *
 * @tparam s   The shift distance, in bits.
 *
 * @param  a   The values to shift.
 */
template <int s> ASTCENC_SIMD_INLINE vint8 lsl(vint8 a)
{
	return vint8(svlsl_n_s32_x(svptrue_b32(), a.m, s));
}

/**
 * @brief Logical shift right.
 *
 * @tparam s   The shift distance, in bits.
 *
 * @param  a   The values to shift.
 */
template <int s> ASTCENC_SIMD_INLINE vint8 lsr(vint8 a)
{
	svuint32_8_t r = svreinterpret_u32_s32(a.m);
	r = svlsr_n_u32_x(svptrue_b32(), r, s);
	return vint8(svreinterpret_s32_u32(r));
}

/**
 * @brief Arithmetic shift right.
 *
 * @tparam s   The shift distance, in bits.
 *
 * @param  a   The values to shift.
 */
template <int s> ASTCENC_SIMD_INLINE vint8 asr(vint8 a)
{
	return vint8(svasr_n_s32_x(svptrue_b32(), a.m, s));
}

/**
 * @brief Return the min vector of two vectors.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vint8 min(vint8 a, vint8 b)
{
	return vint8(svmin_s32_x(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Return the max vector of two vectors.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vint8 max(vint8 a, vint8 b)
{
	return vint8(svmax_s32_x(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Return the horizontal minimum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE vint8 hmin(vint8 a)
{
	return vint8(svminv_s32(svptrue_b32(), a.m));
}

/**
 * @brief Return the horizontal minimum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE int hmin_s(vint8 a)
{
	return svminv_s32(svptrue_b32(), a.m);
}

/**
 * @brief Return the horizontal maximum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE vint8 hmax(vint8 a)
{
	return vint8(svmaxv_s32(svptrue_b32(), a.m));
}

/**
 * @brief Return the horizontal maximum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE int hmax_s(vint8 a)
{
	return svmaxv_s32(svptrue_b32(), a.m);
}

/**
 * @brief Store a vector to a 32B aligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The aligned destination address.
 */
ASTCENC_SIMD_INLINE void storea(vint8 a, int* p)
{
	svst1_s32(svptrue_b32(), p, a.m);
}

/**
 * @brief Store a vector to an unaligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void store(vint8 a, int* p)
{
	svst1_s32(svptrue_b32(), p, a.m);
}

/**
 * @brief Store lowest N (vector width) bytes into an unaligned address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void store_nbytes(vint8 a, uint8_t* p)
{
	svuint8_8_t r = svreinterpret_u8_s32(a.m);
	svst1_u8(svptrue_pat_b8(SV_VL8), p, r);
}

/**
 * @brief Pack the low 8 bits of each lane and store to an unaligned address.
 *
 * @param      v   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void pack_and_store_low_bytes(vint8 v, uint8_t* p)
{
	svuint32_8_t data = svreinterpret_u32_s32(v.m);
	svst1b_u32(svptrue_b32(), p, data);
}

/**
 * @brief Return lanes from @c b if @c cond is set, else @c a.
 *
 * @param a      The vector selected where @c cond is false.
 * @param b      The vector selected where @c cond is true.
 * @param cond   The selection mask.
 */
ASTCENC_SIMD_INLINE vint8 select(vint8 a, vint8 b, vmask8 cond)
{
	return vint8(svsel_s32(cond.m, b.m, a.m));
}

// ============================================================================
// vfloat8 operators and functions
// ============================================================================

/**
 * @brief Overload: vector by vector addition.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat8 operator+(vfloat8 a, vfloat8 b)
{
	return vfloat8(svadd_f32_x(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Overload: vector by vector subtraction.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat8 operator-(vfloat8 a, vfloat8 b)
{
	return vfloat8(svsub_f32_x(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Overload: vector by vector multiplication.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat8 operator*(vfloat8 a, vfloat8 b)
{
	return vfloat8(svmul_f32_x(svptrue_b32(), a.m, b.m));
}

// Retain native forms to avoid fixed-length SVE code generation regressions
#define ASTCENC_USE_NATIVE_VFLOAT8_SCALAR_ARITH 1

/**
 * @brief Overload: vector by scalar multiplication.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat8 operator*(vfloat8 a, float b)
{
	return vfloat8(svmul_f32_x(svptrue_b32(), a.m, svdup_f32(b)));
}

/**
 * @brief Overload: scalar by vector multiplication.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat8 operator*(float a, vfloat8 b)
{
	return vfloat8(svmul_f32_x(svptrue_b32(), svdup_f32(a), b.m));
}

/**
 * @brief Overload: vector by vector division.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat8 operator/(vfloat8 a, vfloat8 b)
{
	return vfloat8(svdiv_f32_x(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Overload: vector by scalar division.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat8 operator/(vfloat8 a, float b)
{
	return vfloat8(svdiv_f32_x(svptrue_b32(), a.m, svdup_f32(b)));
}

/**
 * @brief Overload: scalar by vector division.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat8 operator/(float a, vfloat8 b)
{
	return vfloat8(svdiv_f32_x(svptrue_b32(), svdup_f32(a), b.m));
}

/**
 * @brief Overload: vector by vector equality.
 *
 * Returns vector of false mask values if a or b is NaN.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator==(vfloat8 a, vfloat8 b)
{
	return vmask8(svcmpeq_f32(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Overload: vector by vector inequality.
 *
 * Returns vector of true mask values if a or b is NaN.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator!=(vfloat8 a, vfloat8 b)
{
	return vmask8(svcmpne_f32(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Overload: vector by vector less than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator<(vfloat8 a, vfloat8 b)
{
	return vmask8(svcmplt_f32(svptrue_b32(), a.m, b.m));;
}

/**
 * @brief Overload: vector by vector greater than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator>(vfloat8 a, vfloat8 b)
{
	return vmask8(svcmpgt_f32(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Overload: vector by vector less than or equal.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator<=(vfloat8 a, vfloat8 b)
{
	return vmask8(svcmple_f32(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Overload: vector by vector greater than or equal.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator>=(vfloat8 a, vfloat8 b)
{
	return vmask8(svcmpge_f32(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Return the min vector of two vectors.
 *
 * If either lane value is NaN, the other lane will be returned.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vfloat8 min(vfloat8 a, vfloat8 b)
{
	return vfloat8(svminnm_f32_x(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Return the max vector of two vectors.
 *
 * If either lane value is NaN, the other lane will be returned.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vfloat8 max(vfloat8 a, vfloat8 b)
{
	return vfloat8(svmaxnm_f32_x(svptrue_b32(), a.m, b.m));
}

/**
 * @brief Return the absolute value of the float vector.
 *
 * @param a   The values whose absolute values are required.
 */
ASTCENC_SIMD_INLINE vfloat8 abs(vfloat8 a)
{
	return vfloat8(svabs_f32_x(svptrue_b32(), a.m));
}

/**
 * @brief Return a float rounded to the nearest integer value.
 *
 * @param a   The values to round.
 */
ASTCENC_SIMD_INLINE vfloat8 round(vfloat8 a)
{
	return vfloat8(svrintn_f32_x(svptrue_b32(), a.m));
}

/**
 * @brief Return the horizontal minimum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE vfloat8 hmin(vfloat8 a)
{
	return vfloat8(svminnmv_f32(svptrue_b32(), a.m));
}

/**
 * @brief Return the horizontal minimum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE float hmin_s(vfloat8 a)
{
	return svminnmv_f32(svptrue_b32(), a.m);
}

/**
 * @brief Return the horizontal maximum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE vfloat8 hmax(vfloat8 a)
{
	return vfloat8(svmaxnmv_f32(svptrue_b32(), a.m));
}

/**
 * @brief Return the horizontal maximum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE float hmax_s(vfloat8 a)
{
	return svmaxnmv_f32(svptrue_b32(), a.m);
}

/**
 * @brief Return the horizontal sum of a vector as a scalar.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE float hadd_s(vfloat8 a)
{
	// Can't use svaddv - it's not invariant
	vfloat4 lo(svget_neonq_f32(a.m));
	vfloat4 hi(svget_neonq_f32(svext_f32(a.m, a.m, 4)));
	return hadd_s(lo) + hadd_s(hi);
}

/**
 * @brief Return the sqrt of the lanes in the vector.
 *
 * @param a   The values whose square roots are required.
 */
ASTCENC_SIMD_INLINE vfloat8 sqrt(vfloat8 a)
{
	return vfloat8(svsqrt_f32_x(svptrue_b32(), a.m));
}

/**
 * @brief Return lanes from @c b if @c cond is set, else @c a.
 *
 * @param a      The vector selected where @c cond is false.
 * @param b      The vector selected where @c cond is true.
 * @param cond   The selection mask.
 */
ASTCENC_SIMD_INLINE vfloat8 select(vfloat8 a, vfloat8 b, vmask8 cond)
{
	return vfloat8(svsel_f32(cond.m, b.m, a.m));
}

/**
 * @brief Accumulate lane-wise sums for a vector, folded 4-wide.
 *
 * This is invariant with 4-wide implementations.
 *
 * @param[in,out] accum   The accumulator to update.
 * @param         a       The values to accumulate.
 */
ASTCENC_SIMD_INLINE void haccumulate(vfloat4& accum, vfloat8 a)
{
	vfloat4 lo(svget_neonq_f32(a.m));
	haccumulate(accum, lo);

	vfloat4 hi(svget_neonq_f32(svext_f32(a.m, a.m, 4)));
	haccumulate(accum, hi);
}

/**
 * @brief Accumulate masked lane-wise sums for a vector, folded 4-wide.
 *
 * This is NOT invariant with 4-wide implementations.
 *
 * @param[in,out] accum   The accumulator to update.
 * @param         a       The values to accumulate.
 * @param         m       The lanes to include in the accumulation.
 */
ASTCENC_SIMD_INLINE void haccumulate(vfloat8& accum, vfloat8 a, vmask8 m)
{
	accum.m = svadd_f32_m(m.m, accum.m, a.m);
}

/**
 * @brief Gather integers from a table using vector indices.
 *
 * @param base      The source integer table.
 * @param indices   The table index for each lane.
 *
 * @return The selected table entries.
 */
ASTCENC_SIMD_INLINE vint8 gatheri(const int32_t* base, vint8 indices)
{
	return vint8(svld1_gather_s32index_s32(svptrue_b32(), base, indices.m));
}

/**
 * @brief Load a vector of gathered results from an array.
 *
 * @param base      The source float table.
 * @param indices   The table index for each lane.
 */
ASTCENC_SIMD_INLINE vfloat8 gatherf(const float* base, vint8 indices)
{
	return vfloat8(svld1_gather_s32index_f32(svptrue_b32(), base, indices.m));
}

/**
 * @brief Load a vector of gathered results from an array using byte indices from memory.
 *
 * @param base      The source float table.
 * @param indices   The byte-sized table indices, one per lane.
 */
template<>
ASTCENC_SIMD_INLINE vfloat8 gatherf_byte_inds<vfloat8>(const float* base, const uint8_t* indices)
{
	svint32_t offsets = svld1ub_s32(svptrue_b32(), indices);
	return vfloat8(svld1_gather_s32index_f32(svptrue_b32(), base, offsets));
}

/**
 * @brief Store a vector to an unaligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void store(vfloat8 a, float* p)
{
	svst1_f32(svptrue_b32(), p, a.m);
}

/**
 * @brief Store a vector to a 32B aligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The aligned destination address.
 */
ASTCENC_SIMD_INLINE void storea(vfloat8 a, float* p)
{
	svst1_f32(svptrue_b32(), p, a.m);
}

/**
 * @brief Return an integer value for a float vector, using truncation.
 *
 * @param a   The floating-point values to convert.
 */
ASTCENC_SIMD_INLINE vint8 float_to_int(vfloat8 a)
{
	return vint8(svcvt_s32_f32_x(svptrue_b32(), a.m));
}

// Retain the native form to avoid fixed-length SVE code generation regressions
#define ASTCENC_USE_NATIVE_FLOAT_TO_INT_RTN_8 1

/**
 * @brief Return integers rounded to nearest for nonnegative float values.
 *
 * @param a   The nonnegative floating-point values to convert.
 */
ASTCENC_SIMD_INLINE vint8 float_to_int_rtn(vfloat8 a)
{
	a = a + vfloat8(0.5f);
	return vint8(svcvt_s32_f32_x(svptrue_b32(), a.m));
}

/**
 * @brief Return a float value for an integer vector.
 *
 * @param a   The integer values to convert.
 */
ASTCENC_SIMD_INLINE vfloat8 int_to_float(vint8 a)
{
	return vfloat8(svcvt_f32_s32_x(svptrue_b32(), a.m));
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
ASTCENC_SIMD_INLINE vint8 float_as_int(vfloat8 a)
{
	return vint8(svreinterpret_s32_f32(a.m));
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
ASTCENC_SIMD_INLINE vfloat8 int_as_float(vint8 a)
{
	return vfloat8(svreinterpret_f32_s32(a.m));
}

/*
 * Table structure for a 16x 8-bit entry table.
 */
struct vtable8_16x8 {
	svuint8_8_t t0;
};

/*
 * Table structure for a 32x 8-bit entry table.
 */
struct vtable8_32x8 {
	svuint8_8_t t0;
};

/*
 * Table structure for a 64x 8-bit entry table.
 */
struct vtable8_64x8 {
	svuint8_8_t t0;
	svuint8_8_t t1;
};

/**
 * @brief Prepare a vtable lookup table for 16x 8-bit entries.
 *
 * @param[out] table   The lookup table to prepare.
 * @param      data    The source table bytes.
 */
ASTCENC_SIMD_INLINE void vtable_prepare(
	vtable8_16x8& table,
	const uint8_t* data
) {
	// Top half of register will be zeros
	table.t0 = svld1_u8(svptrue_pat_b8(SV_VL16), data);
}

/**
 * @brief Prepare a vtable lookup table for 32x 8-bit entries.
 *
 * @param[out] table   The lookup table to prepare.
 * @param      data    The source table bytes.
 */
ASTCENC_SIMD_INLINE void vtable_prepare(
	vtable8_32x8& table,
	const uint8_t* data
) {
	table.t0 = svld1_u8(svptrue_b8(), data);
}

/**
 * @brief Prepare a vtable lookup table for 64x 8-bit entries.
 *
 * @param[out] table   The lookup table to prepare.
 * @param      data    The source table bytes.
 */
ASTCENC_SIMD_INLINE void vtable_prepare(
	vtable8_64x8& table,
	const uint8_t* data
) {
	table.t0 = svld1_u8(svptrue_b8(), data);
	table.t1 = svld1_u8(svptrue_b8(), data + 32);
}

/**
 * @brief Perform a vtable lookup in a 16x 8-bit table with 32-bit indices.
 *
 * @param tbl   The prepared lookup table.
 * @param idx   The table index for each lane.
 */
ASTCENC_SIMD_INLINE vint8 vtable_lookup_32bit(
	const vtable8_16x8& tbl,
	vint8 idx
) {
	// Set index byte above max index for unused bytes so table lookup returns zero
	svint32_8_t idx_masked = svorr_s32_x(svptrue_b32(), idx.m, svdup_s32(0xFFFFFF00));
	svuint8_8_t idx_bytes = svreinterpret_u8_s32(idx_masked);

	svuint8_8_t result = svtbl_u8(tbl.t0, idx_bytes);
	return vint8(svreinterpret_s32_u8(result));
}

/**
 * @brief Perform a vtable lookup in a 32x 8-bit table with 32-bit indices.
 *
 * @param tbl   The prepared lookup table.
 * @param idx   The table index for each lane.
 */
ASTCENC_SIMD_INLINE vint8 vtable_lookup_32bit(
	const vtable8_32x8& tbl,
	vint8 idx
) {
	// Set index byte above max index for unused bytes so table lookup returns zero
	svint32_8_t idx_masked = svorr_s32_x(svptrue_b32(), idx.m, svdup_s32(0xFFFFFF00));
	svuint8_8_t idx_bytes = svreinterpret_u8_s32(idx_masked);

	svuint8_8_t result = svtbl_u8(tbl.t0, idx_bytes);
	return vint8(svreinterpret_s32_u8(result));
}

/**
 * @brief Perform a vtable lookup in a 64x 8-bit table with 32-bit indices.
 *
 * Future: SVE2 can directly do svtbl2_u8() for a two register table.
 *
 * @param tbl   The prepared lookup table.
 * @param idx   The table index for each lane.
 */
ASTCENC_SIMD_INLINE vint8 vtable_lookup_32bit(
	const vtable8_64x8& tbl,
	vint8 idx
) {
	// Set index byte above max index for unused bytes so table lookup returns zero
	svint32_8_t idxm = svorr_s32_x(svptrue_b32(), idx.m, svdup_s32(0xFFFFFF00));

	svuint8_8_t idxm8 = svreinterpret_u8_s32(idxm);
	svuint8_8_t t0_lookup = svtbl_u8(tbl.t0, idxm8);

	idxm8 = svsub_u8_x(svptrue_b8(), idxm8, svdup_u8(32));
	svuint8_8_t t1_lookup = svtbl_u8(tbl.t1, idxm8);

	svuint8_8_t result = svorr_u8_x(svptrue_b32(), t0_lookup, t1_lookup);
	return vint8(svreinterpret_s32_u8(result));
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
ASTCENC_SIMD_INLINE void store_lanes_masked(uint8_t* base, vint8 data, vmask8 mask)
{
	svst1_s32(mask.m, reinterpret_cast<int32_t*>(base), data.m);
}

#endif  // #ifndef ASTC_VECMATHLIB_SVE_8_H_INCLUDED
