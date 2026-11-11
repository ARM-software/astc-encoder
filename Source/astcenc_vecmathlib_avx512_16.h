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
 * @brief 16x32-bit vectors, implemented using AVX-512 and VBMI.
 *
 * This module implements 16-wide 32-bit float, int, and mask vectors for x86
 * AVX-512.
 *
 * There is a baseline level of functionality provided by all vector widths and
 * implementations. This is implemented using identical function signatures,
 * modulo data type, so we can use them as substitutable implementations in VLA
 * code.
 *
 * gatherf uses native 32-bit permutes for tables that fit in 1-4 ZMM registers,
 * and scalar loads for larger index ranges. It does not use vgatherdps.
 * vgatherf_load hoists those masked ZMM loads so callers can permute many
 * times from one table. vtable_lookup_32bit uses VBMI vpermb over packed byte
 * tables instead of AVX2-style lane-local pshufb.
 */

#if !defined(ASTC_VECMATHLIB_AVX512_16_H_INCLUDED)
#define ASTC_VECMATHLIB_AVX512_16_H_INCLUDED

#if !defined(ASTCENC_SIMD_INLINE)
	#error "Include astcenc_vecmathlib.h, do not include directly"
#endif

#include <cstdio>
#include <cstring>
#include <limits>

// ============================================================================
// vfloat16 data type
// ============================================================================

/**
 * @brief Data type for 16-wide floats.
 */
struct vfloat16
{
	/**
	 * @brief Construct with default initialization.
	 */
	ASTCENC_SIMD_INLINE vfloat16() = default;

	/**
	 * @brief Construct from 16 values loaded from an unaligned address.
	 *
	 * Consider using loada() which is better with vectors if data is aligned
	 * to vector length.
	 *
	 * @param p   The source values.
	 */
	ASTCENC_SIMD_INLINE explicit vfloat16(const float *p)
	{
		m = _mm512_loadu_ps(p);
	}

	/**
	 * @brief Construct from 1 scalar value replicated across all lanes.
	 *
	 * Consider using zero() for constexpr zeros.
	 *
	 * @param a   The value to replicate across all lanes.
	 */
	ASTCENC_SIMD_INLINE explicit vfloat16(float a)
	{
		m = _mm512_set1_ps(a);
	}

	/**
	 * @brief Construct from an existing SIMD register.
	 *
	 * @param a   The SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vfloat16(__m512 a)
	{
		m = a;
	}

	/**
	 * @brief Factory that returns a vector of zeros.
	 */
	static ASTCENC_SIMD_INLINE vfloat16 zero()
	{
		return vfloat16(_mm512_setzero_ps());
	}

	/**
	 * @brief Factory that returns a replicated scalar loaded from memory.
	 *
	 * @param p   The address of the scalar to replicate.
	 */
	static ASTCENC_SIMD_INLINE vfloat16 load1(const float* p)
	{
		return vfloat16(_mm512_set1_ps(*p));
	}

	/**
	 * @brief Factory that returns a vector loaded from 64B aligned memory.
	 *
	 * @param p   The aligned source address.
	 */
	static ASTCENC_SIMD_INLINE vfloat16 loada(const float* p)
	{
		return vfloat16(_mm512_load_ps(p));
	}

	/** @brief The vector register. */
	__m512 m;
};

// ============================================================================
// vint16 data type
// ============================================================================

/**
 * @brief Data type for 16-wide ints.
 */
struct vint16
{
	/**
	 * @brief Construct with default initialization.
	 */
	ASTCENC_SIMD_INLINE vint16() = default;

	/**
	 * @brief Construct from 16 values loaded from an unaligned address.
	 *
	 * Consider using loada() which is better with vectors if data is aligned
	 * to vector length.
	 *
	 * @param p   The source values.
	 */
	ASTCENC_SIMD_INLINE explicit vint16(const int *p)
	{
		m = _mm512_loadu_si512(p);
	}

	/**
	 * @brief Construct from 16 uint8_t loaded from an unaligned address.
	 *
	 * @param p   The source values.
	 */
	ASTCENC_SIMD_INLINE explicit vint16(const uint8_t *p)
	{
		__m128i bytes;
		std::memcpy(&bytes, p, sizeof(bytes));
		m = _mm512_cvtepu8_epi32(bytes);
	}

	/**
	 * @brief Construct from 1 scalar value replicated across all lanes.
	 *
	 * Consider using zero() for constexpr zeros.
	 *
	 * @param a   The value to replicate across all lanes.
	 */
	ASTCENC_SIMD_INLINE explicit vint16(int a)
	{
		m = _mm512_set1_epi32(a);
	}

	/**
	 * @brief Construct from an existing SIMD register.
	 *
	 * @param a   The SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vint16(__m512i a)
	{
		m = a;
	}

	/**
	 * @brief Factory that returns a vector of zeros.
	 */
	static ASTCENC_SIMD_INLINE vint16 zero()
	{
		return vint16(_mm512_setzero_si512());
	}

	/**
	 * @brief Factory that returns a replicated scalar loaded from memory.
	 *
	 * @param p   The address of the scalar to replicate.
	 */
	static ASTCENC_SIMD_INLINE vint16 load1(const int* p)
	{
		return vint16(_mm512_set1_epi32(*p));
	}

	/**
	 * @brief Factory that returns a vector loaded from 64 unaligned bytes.
	 *
	 * @param p   The source bytes.
	 */
	static ASTCENC_SIMD_INLINE vint16 load(const uint8_t* p)
	{
		return vint16(_mm512_loadu_si512(p));
	}

	/**
	 * @brief Factory that returns a vector loaded from 64B aligned memory.
	 *
	 * @param p   The aligned source address.
	 */
	static ASTCENC_SIMD_INLINE vint16 loada(const int* p)
	{
		return vint16(_mm512_load_si512(p));
	}

	/**
	 * @brief Factory that returns a vector containing the lane IDs.
	 */
	static ASTCENC_SIMD_INLINE vint16 lane_id()
	{
		return vint16(_mm512_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7,
		                                8, 9, 10, 11, 12, 13, 14, 15));
	}

	/** @brief The vector register. */
	__m512i m;
};

// ============================================================================
// vmask16 data type (AVX-512 k-mask)
// ============================================================================

/**
 * @brief Data type for 16-wide control plane masks.
 */
struct vmask16
{
	/**
	 * @brief Construct from an existing SIMD mask register.
	 *
	 * @param a   The SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vmask16(__mmask16 a)
	{
		m = a;
	}

	/**
	 * @brief Construct from 1 scalar value.
	 *
	 * @param a   The value to replicate across all lanes.
	 */
	ASTCENC_SIMD_INLINE explicit vmask16(bool a)
	{
		if (a)
		{
			m = static_cast<__mmask16>(0xFFFF);
		}
		else
		{
			m = static_cast<__mmask16>(0);
		}
	}

	/** @brief The lane mask register. */
	__mmask16 m;
};

// ============================================================================
// Quad operations (four consecutive lanes)
// ============================================================================

#define ASTCENC_HAS_NATIVE_QUAD_OPS

/**
 * @brief Store lane zero from each quad as four consecutive floats.
 *
 * @param      a   The input vector.
 * @param[out] p   The unaligned output array of four floats.
 */
ASTCENC_SIMD_INLINE void store_lane0_quads(vfloat16 a, float* p)
{
	__m512i indices = _mm512_setr_epi32(0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12, 0, 4, 8, 12);
	__m512 compacted = _mm512_permutexvar_ps(indices, a.m);
	_mm_storeu_ps(p, _mm512_castps512_ps128(compacted));
}

/**
 * @brief Pack four four-lane float vectors into one sixteen-lane vector.
 *
 * @param groups   The four source groups, in lane order.
 *
 * @return The packed vector with source groups in consecutive quads.
 */
ASTCENC_SIMD_INLINE vfloat16 pack_quads(const vfloat4* groups)
{
	__m512 values = _mm512_castps128_ps512(groups[0].m);
	values = _mm512_insertf32x4(values, groups[1].m, 1);
	values = _mm512_insertf32x4(values, groups[2].m, 2);
	values = _mm512_insertf32x4(values, groups[3].m, 3);
	return vfloat16(values);
}

/**
 * @brief Pack four four-lane integer vectors into one sixteen-lane vector.
 *
 * @param groups   The four source groups, in lane order.
 *
 * @return The packed vector with source groups in consecutive quads.
 */
ASTCENC_SIMD_INLINE vint16 pack_quads(const vint4* groups)
{
	__m512i values = _mm512_castsi128_si512(groups[0].m);
	values = _mm512_inserti32x4(values, groups[1].m, 1);
	values = _mm512_inserti32x4(values, groups[2].m, 2);
	values = _mm512_inserti32x4(values, groups[3].m, 3);
	return vint16(values);
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
ASTCENC_SIMD_INLINE vfloat16 swizzle_quads(vfloat16 a)
{
	static_assert(l0 >= 0 && l0 < 4 && l1 >= 0 && l1 < 4 &&
	              l2 >= 0 && l2 < 4 && l3 >= 0 && l3 < 4, "Invalid four-lane swizzle");
	return vfloat16(_mm512_permute_ps(a.m, _MM_SHUFFLE(l3, l2, l1, l0)));
}

// ============================================================================
// vmask16 operators and functions
// ============================================================================

/**
 * @brief Overload: mask union (or).
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask16 operator|(vmask16 a, vmask16 b)
{
	return vmask16(static_cast<__mmask16>(a.m | b.m));
}

/**
 * @brief Overload: mask intersect (and).
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask16 operator&(vmask16 a, vmask16 b)
{
	return vmask16(static_cast<__mmask16>(a.m & b.m));
}

/**
 * @brief Overload: mask difference (xor).
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask16 operator^(vmask16 a, vmask16 b)
{
	return vmask16(static_cast<__mmask16>(a.m ^ b.m));
}

/**
 * @brief Overload: mask invert (not).
 *
 * @param a   The input mask.
 */
ASTCENC_SIMD_INLINE vmask16 operator~(vmask16 a)
{
	return vmask16(static_cast<__mmask16>(a.m ^ 0xFFFF));
}

/**
 * @brief Return a 16-bit mask code indicating mask status.
 *
 * bit0 = lane 0
 *
 * @param a   The input mask.
 */
ASTCENC_SIMD_INLINE unsigned int mask(vmask16 a)
{
	return static_cast<unsigned int>(a.m);
}

/**
 * @brief True if any lanes are enabled, false otherwise.
 *
 * @param a   The input mask.
 */
ASTCENC_SIMD_INLINE bool any(vmask16 a)
{
	return a.m != 0;
}

/**
 * @brief True if all lanes are enabled, false otherwise.
 *
 * @param a   The input mask.
 */
ASTCENC_SIMD_INLINE bool all(vmask16 a)
{
	return a.m == 0xFFFF;
}

// ============================================================================
// vint16 operators and functions
// ============================================================================

/**
 * @brief Overload: vector by vector addition.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint16 operator+(vint16 a, vint16 b)
{
	return vint16(_mm512_add_epi32(a.m, b.m));
}

/**
 * @brief Overload: vector by vector incremental addition.
 *
 * @param[in,out] a   The accumulator to update.
 * @param         b   The values to add.
 */
ASTCENC_SIMD_INLINE vint16& operator+=(vint16& a, const vint16& b)
{
	a = a + b;
	return a;
}

/**
 * @brief Overload: vector by vector subtraction.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint16 operator-(vint16 a, vint16 b)
{
	return vint16(_mm512_sub_epi32(a.m, b.m));
}

/**
 * @brief Overload: vector by vector multiplication.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint16 operator*(vint16 a, vint16 b)
{
	return vint16(_mm512_mullo_epi32(a.m, b.m));
}

/**
 * @brief Overload: vector bit invert.
 *
 * @param a   The input vector.
 */
ASTCENC_SIMD_INLINE vint16 operator~(vint16 a)
{
	return vint16(_mm512_xor_si512(a.m, _mm512_set1_epi32(-1)));
}

/**
 * @brief Overload: vector by vector bitwise or.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint16 operator|(vint16 a, vint16 b)
{
	return vint16(_mm512_or_si512(a.m, b.m));
}

/**
 * @brief Overload: vector by vector bitwise and.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint16 operator&(vint16 a, vint16 b)
{
	return vint16(_mm512_and_si512(a.m, b.m));
}

/**
 * @brief Overload: vector by vector bitwise xor.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint16 operator^(vint16 a, vint16 b)
{
	return vint16(_mm512_xor_si512(a.m, b.m));
}

/**
 * @brief Overload: vector by vector equality.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask16 operator==(vint16 a, vint16 b)
{
	return vmask16(_mm512_cmpeq_epi32_mask(a.m, b.m));
}

/**
 * @brief Overload: vector by vector inequality.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask16 operator!=(vint16 a, vint16 b)
{
	return vmask16(_mm512_cmpneq_epi32_mask(a.m, b.m));
}

/**
 * @brief Overload: vector by vector less than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask16 operator<(vint16 a, vint16 b)
{
	return vmask16(_mm512_cmplt_epi32_mask(a.m, b.m));
}

/**
 * @brief Overload: vector by vector greater than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask16 operator>(vint16 a, vint16 b)
{
	return vmask16(_mm512_cmpgt_epi32_mask(a.m, b.m));
}

/**
 * @brief Logical shift left.
 *
 * @tparam s   The shift distance, in bits.
 *
 * @param  a   The values to shift.
 */
template <int s> ASTCENC_SIMD_INLINE vint16 lsl(vint16 a)
{
	return vint16(_mm512_slli_epi32(a.m, s));
}

/**
 * @brief Logical shift right.
 *
 * @tparam s   The shift distance, in bits.
 *
 * @param  a   The values to shift.
 */
template <int s> ASTCENC_SIMD_INLINE vint16 lsr(vint16 a)
{
	return vint16(_mm512_srli_epi32(a.m, s));
}

/**
 * @brief Arithmetic shift right.
 *
 * @tparam s   The shift distance, in bits.
 *
 * @param  a   The values to shift.
 */
template <int s> ASTCENC_SIMD_INLINE vint16 asr(vint16 a)
{
	return vint16(_mm512_srai_epi32(a.m, s));
}

/**
 * @brief Return the min vector of two vectors.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vint16 min(vint16 a, vint16 b)
{
	return vint16(_mm512_min_epi32(a.m, b.m));
}

/**
 * @brief Return the max vector of two vectors.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vint16 max(vint16 a, vint16 b)
{
	return vint16(_mm512_max_epi32(a.m, b.m));
}

/**
 * @brief Return the horizontal minimum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE vint16 hmin(vint16 a)
{
	return vint16(_mm512_reduce_min_epi32(a.m));
}

/**
 * @brief Return the horizontal minimum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE int hmin_s(vint16 a)
{
	return _mm512_reduce_min_epi32(a.m);
}

/**
 * @brief Return the horizontal maximum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE vint16 hmax(vint16 a)
{
	return vint16(_mm512_reduce_max_epi32(a.m));
}

/**
 * @brief Return the horizontal maximum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE int hmax_s(vint16 a)
{
	return _mm512_reduce_max_epi32(a.m);
}

/**
 * @brief Generate a vint16 from a size_t.
 *
 * @param a   The scalar size to replicate across all lanes.
 */
ASTCENC_SIMD_INLINE vint16 vint16_from_size(size_t a)
{
	assert(a <= std::numeric_limits<int>::max());
	return vint16(static_cast<int>(a));
}

/**
 * @brief Store a vector to a 64B aligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The aligned destination address.
 */
ASTCENC_SIMD_INLINE void storea(vint16 a, int* p)
{
	_mm512_store_si512(p, a.m);
}

/**
 * @brief Store a vector to an unaligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void store(vint16 a, int* p)
{
	_mm512_storeu_si512(p, a.m);
}

/**
 * @brief Store lowest N (vector width) bytes into an unaligned address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void store_nbytes(vint16 a, uint8_t* p)
{
	_mm_storeu_si128(reinterpret_cast<__m128i*>(p),
	                 _mm512_extracti32x4_epi32(a.m, 0));
}

/**
 * @brief Pack the low 8 bits of each lane and store to an unaligned address.
 *
 * @param      v   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void pack_and_store_low_bytes(vint16 v, uint8_t* p)
{
	_mm_storeu_si128(reinterpret_cast<__m128i*>(p), _mm512_cvtepi32_epi8(v.m));
}

/**
 * @brief Return lanes from @c b if @c cond is set, else @c a.
 *
 * @param a      The vector selected where @c cond is false.
 * @param b      The vector selected where @c cond is true.
 * @param cond   The selection mask.
 */
ASTCENC_SIMD_INLINE vint16 select(vint16 a, vint16 b, vmask16 cond)
{
	return vint16(_mm512_mask_blend_epi32(cond.m, a.m, b.m));
}

// ============================================================================
// vfloat16 operators and functions
// ============================================================================

/**
 * @brief Overload: vector by vector addition.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat16 operator+(vfloat16 a, vfloat16 b)
{
	return vfloat16(_mm512_add_ps(a.m, b.m));
}

/**
 * @brief Overload: vector by vector incremental addition.
 *
 * @param[in,out] a   The accumulator to update.
 * @param         b   The values to add.
 */
ASTCENC_SIMD_INLINE vfloat16& operator+=(vfloat16& a, const vfloat16& b)
{
	a = a + b;
	return a;
}

/**
 * @brief Overload: vector by vector subtraction.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat16 operator-(vfloat16 a, vfloat16 b)
{
	return vfloat16(_mm512_sub_ps(a.m, b.m));
}

/**
 * @brief Overload: vector by vector multiplication.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat16 operator*(vfloat16 a, vfloat16 b)
{
	return vfloat16(_mm512_mul_ps(a.m, b.m));
}

/**
 * @brief Overload: vector by scalar multiplication.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat16 operator*(vfloat16 a, float b)
{
	return vfloat16(_mm512_mul_ps(a.m, _mm512_set1_ps(b)));
}

/**
 * @brief Overload: scalar by vector multiplication.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat16 operator*(float a, vfloat16 b)
{
	return vfloat16(_mm512_mul_ps(_mm512_set1_ps(a), b.m));
}

/**
 * @brief Overload: vector by vector division.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat16 operator/(vfloat16 a, vfloat16 b)
{
	return vfloat16(_mm512_div_ps(a.m, b.m));
}

/**
 * @brief Overload: vector by scalar division.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat16 operator/(vfloat16 a, float b)
{
	return vfloat16(_mm512_div_ps(a.m, _mm512_set1_ps(b)));
}

/**
 * @brief Overload: scalar by vector division.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat16 operator/(float a, vfloat16 b)
{
	return vfloat16(_mm512_div_ps(_mm512_set1_ps(a), b.m));
}

/**
 * @brief Overload: vector by vector equality.
 *
 * Returns vector of false mask values if a or b is NaN.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask16 operator==(vfloat16 a, vfloat16 b)
{
	return vmask16(_mm512_cmp_ps_mask(a.m, b.m, _CMP_EQ_OQ));
}

/**
 * @brief Overload: vector by vector inequality.
 *
 * Returns vector of true mask values if a or b is NaN.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask16 operator!=(vfloat16 a, vfloat16 b)
{
	return vmask16(_mm512_cmp_ps_mask(a.m, b.m, _CMP_NEQ_UQ));
}

/**
 * @brief Overload: vector by vector less than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask16 operator<(vfloat16 a, vfloat16 b)
{
	return vmask16(_mm512_cmp_ps_mask(a.m, b.m, _CMP_LT_OQ));
}

/**
 * @brief Overload: vector by vector greater than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask16 operator>(vfloat16 a, vfloat16 b)
{
	return vmask16(_mm512_cmp_ps_mask(a.m, b.m, _CMP_GT_OQ));
}

/**
 * @brief Overload: vector by vector less than or equal.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask16 operator<=(vfloat16 a, vfloat16 b)
{
	return vmask16(_mm512_cmp_ps_mask(a.m, b.m, _CMP_LE_OQ));
}

/**
 * @brief Overload: vector by vector greater than or equal.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask16 operator>=(vfloat16 a, vfloat16 b)
{
	return vmask16(_mm512_cmp_ps_mask(a.m, b.m, _CMP_GE_OQ));
}

/**
 * @brief Return the min vector of two vectors.
 *
 * If either lane value is NaN, @c b will be returned for that lane.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vfloat16 min(vfloat16 a, vfloat16 b)
{
	return vfloat16(_mm512_min_ps(a.m, b.m));
}

/**
 * @brief Return the min vector of a vector and a scalar.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vfloat16 min(vfloat16 a, float b)
{
	return min(a, vfloat16(b));
}

/**
 * @brief Return the max vector of two vectors.
 *
 * If either lane value is NaN, @c b will be returned for that lane.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vfloat16 max(vfloat16 a, vfloat16 b)
{
	return vfloat16(_mm512_max_ps(a.m, b.m));
}

/**
 * @brief Return the max vector of a vector and a scalar.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vfloat16 max(vfloat16 a, float b)
{
	return max(a, vfloat16(b));
}

/**
 * @brief Return the clamped value between min and max.
 *
 * It is assumed that neither @c minv nor @c maxv are NaN values. If @c a is
 * NaN then @c minv will be returned for that lane.
 *
 * @param minv   The inclusive lower bound.
 * @param maxv   The inclusive upper bound.
 * @param a      The values to clamp.
 */
ASTCENC_SIMD_INLINE vfloat16 clamp(float minv, float maxv, vfloat16 a)
{
	a.m = _mm512_max_ps(a.m, _mm512_set1_ps(minv));
	a.m = _mm512_min_ps(a.m, _mm512_set1_ps(maxv));
	return a;
}

/**
 * @brief Return the clamped value between 0.0f and 1.0f.
 *
 * If @c a is NaN then zero will be returned for that lane.
 *
 * @param a   The values to clamp.
 */
ASTCENC_SIMD_INLINE vfloat16 clampzo(vfloat16 a)
{
	a.m = _mm512_max_ps(a.m, _mm512_setzero_ps());
	a.m = _mm512_min_ps(a.m, _mm512_set1_ps(1.0f));
	return a;
}

/**
 * @brief Return the absolute value of the float vector.
 *
 * @param a   The values whose absolute values are required.
 */
ASTCENC_SIMD_INLINE vfloat16 abs(vfloat16 a)
{
	return vfloat16(_mm512_and_ps(a.m, _mm512_castsi512_ps(_mm512_set1_epi32(0x7fffffff))));
}

/**
 * @brief Return a float rounded to the nearest integer value.
 *
 * @param a   The values to round.
 */
ASTCENC_SIMD_INLINE vfloat16 round(vfloat16 a)
{
	return vfloat16(_mm512_roundscale_ps(a.m, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC));
}

/**
 * @brief Return the horizontal minimum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE vfloat16 hmin(vfloat16 a)
{
	return vfloat16(_mm512_reduce_min_ps(a.m));
}

/**
 * @brief Return the horizontal minimum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE float hmin_s(vfloat16 a)
{
	return _mm512_reduce_min_ps(a.m);
}

/**
 * @brief Return the horizontal maximum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE vfloat16 hmax(vfloat16 a)
{
	return vfloat16(_mm512_reduce_max_ps(a.m));
}

/**
 * @brief Return the horizontal maximum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE float hmax_s(vfloat16 a)
{
	return _mm512_reduce_max_ps(a.m);
}

/**
 * @brief Return the horizontal sum of a vector as a scalar.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE float hadd_s(vfloat16 a)
{
	// Reduce each 4-wide chunk, then combine pairs of chunk sums
	vfloat4 q0(_mm512_extractf32x4_ps(a.m, 0));
	vfloat4 q1(_mm512_extractf32x4_ps(a.m, 1));
	vfloat4 q2(_mm512_extractf32x4_ps(a.m, 2));
	vfloat4 q3(_mm512_extractf32x4_ps(a.m, 3));
	return (hadd_s(q0) + hadd_s(q1)) + (hadd_s(q2) + hadd_s(q3));
}

/**
 * @brief Return the sqrt of the lanes in the vector.
 *
 * @param a   The values whose square roots are required.
 */
ASTCENC_SIMD_INLINE vfloat16 sqrt(vfloat16 a)
{
	return vfloat16(_mm512_sqrt_ps(a.m));
}

/**
 * @brief Return lanes from @c b if @c cond is set, else @c a.
 *
 * @param a      The vector selected where @c cond is false.
 * @param b      The vector selected where @c cond is true.
 * @param cond   The selection mask.
 */
ASTCENC_SIMD_INLINE vfloat16 select(vfloat16 a, vfloat16 b, vmask16 cond)
{
	return vfloat16(_mm512_mask_blend_ps(cond.m, a.m, b.m));
}

/**
 * @brief Accumulate lane-wise sums for a vector, folded 4-wide.
 *
 * This is invariant with 4-wide implementations.
 *
 * @param[in,out] accum   The accumulator to update.
 * @param         a       The values to accumulate.
 */
ASTCENC_SIMD_INLINE void haccumulate(vfloat4& accum, vfloat16 a)
{
	haccumulate(accum, vfloat4(_mm512_extractf32x4_ps(a.m, 0)));
	haccumulate(accum, vfloat4(_mm512_extractf32x4_ps(a.m, 1)));
	haccumulate(accum, vfloat4(_mm512_extractf32x4_ps(a.m, 2)));
	haccumulate(accum, vfloat4(_mm512_extractf32x4_ps(a.m, 3)));
}

/**
 * @brief Accumulate lane-wise sums for a vector.
 *
 * This is NOT invariant with 4-wide implementations.
 *
 * @param[in,out] accum   The accumulator to update.
 * @param         a       The values to accumulate.
 */
ASTCENC_SIMD_INLINE void haccumulate(vfloat16& accum, vfloat16 a)
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
ASTCENC_SIMD_INLINE void haccumulate(vfloat4& accum, vfloat16 a, vmask16 m)
{
	a = select(vfloat16::zero(), a, m);
	haccumulate(accum, a);
}

/**
 * @brief Accumulate masked lane-wise sums for a vector.
 *
 * This is NOT invariant with 4-wide implementations.
 *
 * @param[in,out] accum   The accumulator to update.
 * @param         a       The values to accumulate.
 * @param         m       The lanes to include in the accumulation.
 */
ASTCENC_SIMD_INLINE void haccumulate(vfloat16& accum, vfloat16 a, vmask16 m)
{
	a = select(vfloat16::zero(), a, m);
	haccumulate(accum, a);
}

/**
 * @brief Return a mask for the live lanes of a 16-entry float table load.
 *
 * @param count    The number of entries in the table.
 * @param offset   The first table entry covered by the load.
 *
 * @return The lane mask, with no lanes enabled if offset is past the table.
 */
ASTCENC_SIMD_INLINE __mmask16 vbmi_load_mask(unsigned int count, unsigned int offset)
{
	if (count <= offset)
	{
		return 0;
	}

	unsigned int n = count - offset;
	if (n >= 16)
	{
		return static_cast<__mmask16>(0xFFFF);
	}

	return static_cast<__mmask16>((1u << n) - 1u);
}

/**
 * @brief Gather integers from a table using vector indices.
 *
 * @param base      The source integer table.
 * @param indices   The table index for each lane.
 *
 * @return The selected table entries.
 */
ASTCENC_SIMD_INLINE vint16 gatheri(const int32_t* base, vint16 indices)
{
	return vint16(_mm512_i32gather_epi32(indices.m, base, 4));
}

/**
 * @brief Load a vector of gathered results from an array using scalar loads.
 *
 * @param base      The source float table.
 * @param indices   The table index for each lane.
 */
ASTCENC_SIMD_INLINE vfloat16 gatherf_scalar(const float* base, vint16 indices)
{
	ASTCENC_ALIGNAS int idx[16];
	_mm512_store_si512(idx, indices.m);
	return vfloat16(_mm512_setr_ps(
	    base[idx[0]],  base[idx[1]],  base[idx[2]],  base[idx[3]],
	    base[idx[4]],  base[idx[5]],  base[idx[6]],  base[idx[7]],
	    base[idx[8]],  base[idx[9]],  base[idx[10]], base[idx[11]],
	    base[idx[12]], base[idx[13]], base[idx[14]], base[idx[15]]));
}

/**
 * @brief Register-resident float table for gatherf (up to 64 entries).
 *
 * Load once with vgatherf_load; permute many times with gatherf. Tables larger
 * than 64 keep the pointer and use the scalar-load gather.
 */
struct vgatherf_table {
	/** @brief Table entries 0-15, with unused lanes zeroed. */
	__m512 t0;

	/** @brief Table entries 16-31, with unused lanes zeroed. */
	__m512 t1;

	/** @brief Table entries 32-47, with unused lanes zeroed. */
	__m512 t2;

	/** @brief Table entries 48-63, with unused lanes zeroed. */
	__m512 t3;

	/** @brief The source table, used for scalar gathers. */
	const float* base;

	/** @brief The number of entries in the source table. */
	unsigned int count;
};

#define ASTCENC_HAS_VGATHERF_TABLE 1

// GCC 13 reports array bounds errors in unreachable loads after inlining when the table span is
// inferred from a vector of constant gather indices
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ == 13
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#endif

/**
 * @brief Load a float table into registers when it contains up to 64 entries.
 *
 * Larger tables retain the source pointer for scalar gathers. The source
 * table must remain valid while the returned table is in use.
 *
 * @param base    The source float table.
 * @param count   The number of entries in the table.
 *
 * @return The table representation to use with gatherf().
 */
ASTCENC_SIMD_INLINE vgatherf_table vgatherf_load(const float* base, unsigned int count)
{
	vgatherf_table t;
	t.t0 = _mm512_setzero_ps();
	t.t1 = _mm512_setzero_ps();
	t.t2 = _mm512_setzero_ps();
	t.t3 = _mm512_setzero_ps();
	t.base = base;
	t.count = count;

	if (count > 0 && count <= 64)
	{
		t.t0 = _mm512_maskz_loadu_ps(vbmi_load_mask(count, 0), base);

		if (count > 16)
		{
			t.t1 = _mm512_maskz_loadu_ps(vbmi_load_mask(count, 16), base + 16);
		}

		if (count > 32)
		{
			t.t2 = _mm512_maskz_loadu_ps(vbmi_load_mask(count, 32), base + 32);
		}

		if (count > 48)
		{
			t.t3 = _mm512_maskz_loadu_ps(vbmi_load_mask(count, 48), base + 48);
		}
	}

	return t;
}

#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ == 13
#pragma GCC diagnostic pop
#endif

/**
 * @brief Gather 16 floats from a preloaded float table, else scalar loads.
 *
 * Use one 32-bit permute for up to 32 entries. Larger tables use two permutes
 * with the same low 5 index bits, selecting the half with index bit 5.
 *
 * @param t         The prepared float table.
 * @param indices   The table index for each lane.
 */
ASTCENC_SIMD_INLINE vfloat16 gatherf(const vgatherf_table& t, vint16 indices)
{
	if (t.count == 0 || t.count > 64)
	{
		return gatherf_scalar(t.base, indices);
	}

	if (t.count <= 16)
	{
		return vfloat16(_mm512_permutexvar_ps(indices.m, t.t0));
	}

	if (t.count <= 32)
	{
		return vfloat16(_mm512_permutex2var_ps(t.t0, indices.m, t.t1));
	}

	__m512 lo = _mm512_permutex2var_ps(t.t0, indices.m, t.t1);
	__m512 hi = _mm512_permutex2var_ps(t.t2, indices.m, t.t3);
	__mmask16 ge32 = _mm512_cmpge_epi32_mask(indices.m, _mm512_set1_epi32(32));
	return vfloat16(_mm512_mask_blend_ps(ge32, lo, hi));
}

/**
 * @brief Gather from a preloaded table using 16 byte indices from memory.
 *
 * @param t         The prepared float table.
 * @param indices   The byte-sized table indices, one per lane.
 */
ASTCENC_SIMD_INLINE vfloat16 gatherf(const vgatherf_table& t, const uint8_t* indices)
{
	return gatherf(t, vint16(indices));
}

/**
 * @brief Gather 16 floats with AVX-512 32-bit permutes, else scalar loads.
 *
 * Convenience wrapper: infers the live table span from the index vector
 * (hmax + 1) then permutes. Prefer vgatherf_load when the same table is
 * gathered many times. Masked table loads use only [0, max_index] so
 * mid-array callers stay in bounds.
 *
 * @param base      The source float table.
 * @param indices   The table index for each lane.
 */
ASTCENC_SIMD_INLINE vfloat16 gatherf(const float* base, vint16 indices)
{
	int mx = hmax_s(indices);
	if (mx < 0)
	{
		return gatherf_scalar(base, indices);
	}

	return gatherf(vgatherf_load(base, static_cast<unsigned int>(mx) + 1u), indices);
}

/**
 * @brief Load a vector of gathered results from an array using byte indices from memory.
 *
 * @param base      The source float table.
 * @param indices   The byte-sized table indices, one per lane.
 */
template<>
ASTCENC_SIMD_INLINE vfloat16 gatherf_byte_inds<vfloat16>(const float* base, const uint8_t* indices)
{
	return gatherf(base, vint16(indices));
}

/**
 * @brief Store a vector to an unaligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void store(vfloat16 a, float* p)
{
	_mm512_storeu_ps(p, a.m);
}

/**
 * @brief Store a vector to a 64B aligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The aligned destination address.
 */
ASTCENC_SIMD_INLINE void storea(vfloat16 a, float* p)
{
	_mm512_store_ps(p, a.m);
}

/**
 * @brief Return an integer value for a float vector, using truncation.
 *
 * @param a   The floating-point values to convert.
 */
ASTCENC_SIMD_INLINE vint16 float_to_int(vfloat16 a)
{
	return vint16(_mm512_cvttps_epi32(a.m));
}

/**
 * @brief Return integers rounded to nearest for nonnegative float values.
 *
 * Halfway values round upwards.
 *
 * @param a   The nonnegative floating-point values to convert.
 */
ASTCENC_SIMD_INLINE vint16 float_to_int_rtn(vfloat16 a)
{
	a = a + vfloat16(0.5f);
	return vint16(_mm512_cvttps_epi32(a.m));
}

/**
 * @brief Return a float value for an integer vector.
 *
 * @param a   The integer values to convert.
 */
ASTCENC_SIMD_INLINE vfloat16 int_to_float(vint16 a)
{
	return vfloat16(_mm512_cvtepi32_ps(a.m));
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
ASTCENC_SIMD_INLINE vint16 float_as_int(vfloat16 a)
{
	return vint16(_mm512_castps_si512(a.m));
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
ASTCENC_SIMD_INLINE vfloat16 int_as_float(vint16 a)
{
	return vfloat16(_mm512_castsi512_ps(a.m));
}

/**
 * @brief Table structure for 16x 8-bit entries.
 *
 * VBMI can index any of 64 bytes in one register, so each table size uses a
 * single register without the lane replication required by AVX2.
 */
struct vtable16_16x8 {
	/** @brief The packed byte table. */
	vint16 t0;
};

/**
 * @brief Table structure for 32x 8-bit entries.
 */
struct vtable16_32x8 {
	/** @brief The packed byte table. */
	vint16 t0;
};

/**
 * @brief Table structure for 64x 8-bit entries.
 */
struct vtable16_64x8 {
	/** @brief The packed byte table. */
	vint16 t0;
};

/**
 * @brief Look up packed bytes using VBMI, returning one byte per 32-bit lane.
 *
 * @param table   The packed byte table.
 * @param idx     The table index for each lane.
 */
ASTCENC_SIMD_INLINE vint16 vtable_lookup_32bit_vbmi(vint16 table, vint16 idx)
{
	__m512i r = _mm512_permutexvar_epi8(idx.m, table.m);
	return vint16(_mm512_and_si512(r, _mm512_set1_epi32(0xFF)));
}

/**
 * @brief Prepare a vtable lookup table for 16x 8-bit entries.
 *
 * @param[out] table   The lookup table to prepare.
 * @param      data    The source table bytes.
 */
ASTCENC_SIMD_INLINE void vtable_prepare(
    vtable16_16x8& table,
    const uint8_t* data
) {
	table.t0 = vint16(_mm512_zextsi128_si512(_mm_loadu_si128(
	    reinterpret_cast<const __m128i*>(data))));
}

/**
 * @brief Prepare a vtable lookup table for 32x 8-bit entries.
 *
 * @param[out] table   The lookup table to prepare.
 * @param      data    The source table bytes.
 */
ASTCENC_SIMD_INLINE void vtable_prepare(
    vtable16_32x8& table,
    const uint8_t* data
) {
	table.t0 = vint16(_mm512_zextsi256_si512(_mm256_loadu_si256(
	    reinterpret_cast<const __m256i*>(data))));
}

/**
 * @brief Prepare a vtable lookup table for 64x 8-bit entries.
 *
 * @param[out] table   The lookup table to prepare.
 * @param      data    The source table bytes.
 */
ASTCENC_SIMD_INLINE void vtable_prepare(
    vtable16_64x8& table,
    const uint8_t* data
) {
	table.t0 = vint16(_mm512_loadu_si512(data));
}

/**
 * @brief Perform a vtable lookup in a 16x 8-bit table with 32-bit indices.
 *
 * @param tbl   The prepared lookup table.
 * @param idx   The table index for each lane.
 */
ASTCENC_SIMD_INLINE vint16 vtable_lookup_32bit(
    const vtable16_16x8& tbl,
    vint16 idx
) {
	return vtable_lookup_32bit_vbmi(tbl.t0, idx);
}

/**
 * @brief Perform a vtable lookup in a 32x 8-bit table with 32-bit indices.
 *
 * @param tbl   The prepared lookup table.
 * @param idx   The table index for each lane.
 */
ASTCENC_SIMD_INLINE vint16 vtable_lookup_32bit(
    const vtable16_32x8& tbl,
    vint16 idx
) {
	return vtable_lookup_32bit_vbmi(tbl.t0, idx);
}

/**
 * @brief Perform a vtable lookup in a 64x 8-bit table with 32-bit indices.
 *
 * @param tbl   The prepared lookup table.
 * @param idx   The table index for each lane.
 */
ASTCENC_SIMD_INLINE vint16 vtable_lookup_32bit(
    const vtable16_64x8& tbl,
    vint16 idx
) {
	return vtable_lookup_32bit_vbmi(tbl.t0, idx);
}

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
ASTCENC_SIMD_INLINE vint16 interleave_rgba8(vint16 r, vint16 g, vint16 b, vint16 a)
{
	return r + lsl<8>(g) + lsl<16>(b) + lsl<24>(a);
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
ASTCENC_SIMD_INLINE void store_lanes_masked(uint8_t* base, vint16 data, vmask16 mask)
{
	_mm512_mask_storeu_epi32(base, mask.m, data.m);
}

/**
 * @brief Debug function to print a vector of ints.
 *
 * @param a   The vector to print.
 */
ASTCENC_SIMD_INLINE void print(vint16 a)
{
	ASTCENC_ALIGNAS int v[16];
	storea(a, v);
	printf("v16_i32:\n  %8d %8d %8d %8d %8d %8d %8d %8d\n  %8d %8d %8d %8d %8d %8d %8d %8d\n",
	       v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7],
	       v[8], v[9], v[10], v[11], v[12], v[13], v[14], v[15]);
}

/**
 * @brief Debug function to print a vector of ints.
 *
 * @param a   The vector to print.
 */
ASTCENC_SIMD_INLINE void printx(vint16 a)
{
	ASTCENC_ALIGNAS int v[16];
	storea(a, v);

	unsigned int uv[16];
	std::memcpy(uv, v, sizeof(uv));

	printf("v16_i32:\n  %08x %08x %08x %08x %08x %08x %08x %08x\n  %08x %08x %08x %08x %08x %08x %08x %08x\n",
	       uv[0], uv[1], uv[2], uv[3], uv[4], uv[5], uv[6], uv[7],
	       uv[8], uv[9], uv[10], uv[11], uv[12], uv[13], uv[14], uv[15]);
}

/**
 * @brief Debug function to print a vector of floats.
 *
 * @param a   The vector to print.
 */
ASTCENC_SIMD_INLINE void print(vfloat16 a)
{
	ASTCENC_ALIGNAS float v[16];
	storea(a, v);
	printf("v16_f32:\n");
	for (int i = 0; i < 16; i++)
	{
		printf("  %0.4f", static_cast<double>(v[i]));
		if ((i & 7) == 7)
		{
			printf("\n");
		}
	}
}

/**
 * @brief Debug function to print a vector of masks.
 *
 * @param a   The vector to print.
 */
ASTCENC_SIMD_INLINE void print(vmask16 a)
{
	print(select(vint16(0), vint16(1), a));
}

#endif  // #ifndef ASTC_VECMATHLIB_AVX512_16_H_INCLUDED
