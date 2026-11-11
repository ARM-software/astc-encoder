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
 * @brief 8x32-bit vectors, implemented using AVX2.
 *
 * This module implements 8-wide 32-bit float, int, and mask vectors for x86
 * AVX2, optionally enhanced with 256-bit AVX-512 instructions for AVX256 builds.
 *
 * There is a baseline level of functionality provided by all vector widths and
 * implementations. This is implemented using identical function signatures,
 * modulo data type, so we can use them as substitutable implementations in VLA
 * code.
 */

#if !defined(ASTC_VECMATHLIB_AVX2_8_H_INCLUDED)
#define ASTC_VECMATHLIB_AVX2_8_H_INCLUDED

#if !defined(ASTCENC_SIMD_INLINE)
	#error "Include astcenc_vecmathlib.h, do not include directly"
#endif

#include <cstdio>

// Define convenience intrinsics that are missing on older compilers
/**
 * @brief Combine two 128-bit integer vectors into one 256-bit vector.
 *
 * @param m   The upper 128-bit vector.
 * @param n   The lower 128-bit vector.
 */
#define astcenc_mm256_set_m128i(m, n) _mm256_insertf128_si256(_mm256_castsi128_si256((n)), (m), 1)

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
		m = _mm256_loadu_ps(p);
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
		m = _mm256_set1_ps(a);
	}

	/**
	 * @brief Construct from an existing SIMD register.
	 *
	 * @param a   The SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vfloat8(__m256 a)
	{
		m = a;
	}

	/**
	 * @brief Factory that returns a vector of zeros.
	 */
	static ASTCENC_SIMD_INLINE vfloat8 zero()
	{
		return vfloat8(_mm256_setzero_ps());
	}

	/**
	 * @brief Factory that returns a replicated scalar loaded from memory.
	 *
	 * @param p   The address of the scalar to replicate.
	 */
	static ASTCENC_SIMD_INLINE vfloat8 load1(const float* p)
	{
		return vfloat8(_mm256_broadcast_ss(p));
	}

	/**
	 * @brief Factory that returns a vector loaded from 32B aligned memory.
	 *
	 * @param p   The aligned source address.
	 */
	static ASTCENC_SIMD_INLINE vfloat8 loada(const float* p)
	{
		return vfloat8(_mm256_load_ps(p));
	}

	/**
	 * @brief The vector ...
	 */
	__m256 m;
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
		m = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));
	}

	/**
	 * @brief Construct from 8 uint8_t loaded from an unaligned address.
	 *
	 * @param p   The source values.
	 */
	ASTCENC_SIMD_INLINE explicit vint8(const uint8_t *p)
	{
		// Copy through a uint32_t to avoid load through a reinterpreted pointer, which is both
		// unaligned and an aliasing violation
		uint64_t tmp;
		std::memcpy(&tmp, p, sizeof(tmp));

		// _mm_loadu_si64 would be nicer syntax, but missing on older GCC and generates broken code
		// on GCC 11.x before 11.3
		m = _mm256_cvtepu8_epi32(_mm_cvtsi64_si128(tmp));
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
		m = _mm256_set1_epi32(a);
	}

	/**
	 * @brief Construct from an existing SIMD register.
	 *
	 * @param a   The SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vint8(__m256i a)
	{
		m = a;
	}

	/**
	 * @brief Factory that returns a vector of zeros.
	 */
	static ASTCENC_SIMD_INLINE vint8 zero()
	{
		return vint8(_mm256_setzero_si256());
	}

	/**
	 * @brief Factory that returns a replicated scalar loaded from memory.
	 *
	 * @param p   The address of the scalar to replicate.
	 */
	static ASTCENC_SIMD_INLINE vint8 load1(const int* p)
	{
		__m128i a = _mm_set1_epi32(*p);
		return vint8(_mm256_broadcastd_epi32(a));
	}

	/**
	 * @brief Factory that returns a vector loaded from unaligned memory.
	 *
	 * @param p   The source bytes.
	 */
	static ASTCENC_SIMD_INLINE vint8 load(const uint8_t* p)
	{
		return vint8(_mm256_lddqu_si256(reinterpret_cast<const __m256i*>(p)));
	}

	/**
	 * @brief Factory that returns a vector loaded from 32B aligned memory.
	 *
	 * @param p   The aligned source address.
	 */
	static ASTCENC_SIMD_INLINE vint8 loada(const int* p)
	{
		return vint8(_mm256_load_si256(reinterpret_cast<const __m256i*>(p)));
	}

	/**
	 * @brief Factory that returns a vector containing the lane IDs.
	 */
	static ASTCENC_SIMD_INLINE vint8 lane_id()
	{
		return vint8(_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0));
	}

	/**
	 * @brief The vector ...
	 */
	__m256i m;
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
	ASTCENC_SIMD_INLINE explicit vmask8(__m256 a)
	{
		m = a;
	}

	/**
	 * @brief Construct from an existing SIMD register.
	 *
	 * @param a   The SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vmask8(__m256i a)
	{
		m = _mm256_castsi256_ps(a);
	}

	/**
	 * @brief Construct from 1 scalar value.
	 *
	 * @param a   The value to replicate across all lanes.
	 */
	ASTCENC_SIMD_INLINE explicit vmask8(bool a)
	{
		vint8 mask(a == false ? 0 : -1);
		m = _mm256_castsi256_ps(mask.m);
	}

	/**
	 * @brief The vector ...
	 */
	__m256 m;
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
	__m256 compacted = _mm256_permutevar8x32_ps(a.m, _mm256_setr_epi32(0, 4, 0, 4, 0, 4, 0, 4));
	__m128i lo = _mm_castps_si128(_mm256_castps256_ps128(compacted));
	_mm_storel_epi64(reinterpret_cast<__m128i*>(p), lo);
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
	return vfloat8(_mm256_insertf128_ps(_mm256_castps128_ps256(groups[0].m), groups[1].m, 1));
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
	return vint8(_mm256_insertf128_si256(_mm256_castsi128_si256(groups[0].m), groups[1].m, 1));
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
	return vfloat8(_mm256_permute_ps(a.m, _MM_SHUFFLE(l3, l2, l1, l0)));
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
	return vmask8(_mm256_or_ps(a.m, b.m));
}

/**
 * @brief Overload: mask intersect (and).
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator&(vmask8 a, vmask8 b)
{
	return vmask8(_mm256_and_ps(a.m, b.m));
}

/**
 * @brief Overload: mask difference (xor).
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator^(vmask8 a, vmask8 b)
{
	return vmask8(_mm256_xor_ps(a.m, b.m));
}

/**
 * @brief Overload: mask invert (not).
 *
 * @param a   The input mask.
 */
ASTCENC_SIMD_INLINE vmask8 operator~(vmask8 a)
{
	return vmask8(_mm256_xor_si256(_mm256_castps_si256(a.m), _mm256_set1_epi32(-1)));
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
	return static_cast<unsigned int>(_mm256_movemask_ps(a.m));
}

/**
 * @brief True if any lanes are enabled, false otherwise.
 *
 * @param a   The input mask.
 */
ASTCENC_SIMD_INLINE bool any(vmask8 a)
{
	return mask(a) != 0;
}

/**
 * @brief True if all lanes are enabled, false otherwise.
 *
 * @param a   The input mask.
 */
ASTCENC_SIMD_INLINE bool all(vmask8 a)
{
	return mask(a) == 0xFF;
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
	return vint8(_mm256_add_epi32(a.m, b.m));
}

/**
 * @brief Overload: vector by vector subtraction.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint8 operator-(vint8 a, vint8 b)
{
	return vint8(_mm256_sub_epi32(a.m, b.m));
}

/**
 * @brief Overload: vector by vector multiplication.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint8 operator*(vint8 a, vint8 b)
{
	return vint8(_mm256_mullo_epi32(a.m, b.m));
}

/**
 * @brief Overload: vector bit invert.
 *
 * @param a   The input vector.
 */
ASTCENC_SIMD_INLINE vint8 operator~(vint8 a)
{
	return vint8(_mm256_xor_si256(a.m, _mm256_set1_epi32(-1)));
}

/**
 * @brief Overload: vector by vector bitwise or.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint8 operator|(vint8 a, vint8 b)
{
	return vint8(_mm256_or_si256(a.m, b.m));
}

/**
 * @brief Overload: vector by vector bitwise and.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint8 operator&(vint8 a, vint8 b)
{
	return vint8(_mm256_and_si256(a.m, b.m));
}

/**
 * @brief Overload: vector by vector bitwise xor.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vint8 operator^(vint8 a, vint8 b)
{
	return vint8(_mm256_xor_si256(a.m, b.m));
}

/**
 * @brief Overload: vector by vector equality.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator==(vint8 a, vint8 b)
{
	return vmask8(_mm256_cmpeq_epi32(a.m, b.m));
}

/**
 * @brief Overload: vector by vector inequality.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator!=(vint8 a, vint8 b)
{
	return ~vmask8(_mm256_cmpeq_epi32(a.m, b.m));
}

/**
 * @brief Overload: vector by vector less than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator<(vint8 a, vint8 b)
{
	return vmask8(_mm256_cmpgt_epi32(b.m, a.m));
}

/**
 * @brief Overload: vector by vector greater than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator>(vint8 a, vint8 b)
{
	return vmask8(_mm256_cmpgt_epi32(a.m, b.m));
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
	return vint8(_mm256_slli_epi32(a.m, s));
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
	return vint8(_mm256_srli_epi32(a.m, s));
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
	return vint8(_mm256_srai_epi32(a.m, s));
}

/**
 * @brief Return the min vector of two vectors.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vint8 min(vint8 a, vint8 b)
{
	return vint8(_mm256_min_epi32(a.m, b.m));
}

/**
 * @brief Return the max vector of two vectors.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vint8 max(vint8 a, vint8 b)
{
	return vint8(_mm256_max_epi32(a.m, b.m));
}

/**
 * @brief Return the horizontal minimum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE vint8 hmin(vint8 a)
{
	// Build min within groups of 2, then 4, then 8
	__m256i m = _mm256_min_epi32(a.m, _mm256_shuffle_epi32(a.m, _MM_SHUFFLE(2, 3, 0, 1)));
	m = _mm256_min_epi32(m, _mm256_shuffle_epi32(m, _MM_SHUFFLE(1, 0, 3, 2)));
	m = _mm256_min_epi32(m, _mm256_permute2x128_si256(m, m, 0x01));

	vint8 vmin(m);
	return vmin;
}

/**
 * @brief Return the horizontal minimum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE int hmin_s(vint8 a)
{
	return _mm256_cvtsi256_si32(hmin(a).m);
}

/**
 * @brief Return the horizontal maximum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE vint8 hmax(vint8 a)
{
	// Build max within groups of 2, then 4, then 8
	__m256i m = _mm256_max_epi32(a.m, _mm256_shuffle_epi32(a.m, _MM_SHUFFLE(2, 3, 0, 1)));
	m = _mm256_max_epi32(m, _mm256_shuffle_epi32(m, _MM_SHUFFLE(1, 0, 3, 2)));
	m = _mm256_max_epi32(m, _mm256_permute2x128_si256(m, m, 0x01));

	vint8 vmax(m);
	return vmax;
}

/**
 * @brief Return the horizontal maximum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE int hmax_s(vint8 a)
{
	return _mm256_cvtsi256_si32(hmax(a).m);
}

/**
 * @brief Store a vector to a 32B aligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The aligned destination address.
 */
ASTCENC_SIMD_INLINE void storea(vint8 a, int* p)
{
	_mm256_store_si256(reinterpret_cast<__m256i*>(p), a.m);
}

/**
 * @brief Store a vector to an unaligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void store(vint8 a, int* p)
{
	_mm256_storeu_si256(reinterpret_cast<__m256i*>(p), a.m);
}

/**
 * @brief Store lowest N (vector width) bytes into an unaligned address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void store_nbytes(vint8 a, uint8_t* p)
{
	// Use _mm_storel_epi64 for compatibility; _mm_storeu_si64 requires GCC 9 or Clang 9
	_mm_storel_epi64(reinterpret_cast<__m128i*>(p), _mm256_extracti128_si256(a.m, 0));
}

/**
 * @brief Pack the low 8 bits of each lane and store to an unaligned address.
 *
 * @param      v   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void pack_and_store_low_bytes(vint8 v, uint8_t* p)
{
#if ASTCENC_AVX >= 3
	_mm256_mask_cvtepi32_storeu_epi8(p, 0xFF, v.m);
#else
	__m256i shuf = _mm256_set_epi8(0, 0, 0, 0,  0,  0,  0,  0,
	                               0, 0, 0, 0, 28, 24, 20, 16,
	                               0, 0, 0, 0,  0,  0,  0,  0,
	                               0, 0, 0, 0, 12,  8,  4,  0);
	__m256i a = _mm256_shuffle_epi8(v.m, shuf);
	__m128i a0 = _mm256_extracti128_si256(a, 0);
	__m128i a1 = _mm256_extracti128_si256(a, 1);
	__m128i b = _mm_unpacklo_epi32(a0, a1);

	__m256i r = astcenc_mm256_set_m128i(b, b);

	store_nbytes(vint8(r), p);
#endif
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
	__m256i condi = _mm256_castps_si256(cond.m);
	return vint8(_mm256_blendv_epi8(a.m, b.m, condi));
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
	return vfloat8(_mm256_add_ps(a.m, b.m));
}

/**
 * @brief Overload: vector by vector subtraction.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat8 operator-(vfloat8 a, vfloat8 b)
{
	return vfloat8(_mm256_sub_ps(a.m, b.m));
}

/**
 * @brief Overload: vector by vector multiplication.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat8 operator*(vfloat8 a, vfloat8 b)
{
	return vfloat8(_mm256_mul_ps(a.m, b.m));
}

/**
 * @brief Overload: vector by vector division.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vfloat8 operator/(vfloat8 a, vfloat8 b)
{
	return vfloat8(_mm256_div_ps(a.m, b.m));
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
	return vmask8(_mm256_cmp_ps(a.m, b.m, _CMP_EQ_OQ));
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
	return vmask8(_mm256_cmp_ps(a.m, b.m, _CMP_NEQ_UQ));
}

/**
 * @brief Overload: vector by vector less than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator<(vfloat8 a, vfloat8 b)
{
	return vmask8(_mm256_cmp_ps(a.m, b.m, _CMP_LT_OQ));
}

/**
 * @brief Overload: vector by vector greater than.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator>(vfloat8 a, vfloat8 b)
{
	return vmask8(_mm256_cmp_ps(a.m, b.m, _CMP_GT_OQ));
}

/**
 * @brief Overload: vector by vector less than or equal.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator<=(vfloat8 a, vfloat8 b)
{
	return vmask8(_mm256_cmp_ps(a.m, b.m, _CMP_LE_OQ));
}

/**
 * @brief Overload: vector by vector greater than or equal.
 *
 * @param a   The left operand.
 * @param b   The right operand.
 */
ASTCENC_SIMD_INLINE vmask8 operator>=(vfloat8 a, vfloat8 b)
{
	return vmask8(_mm256_cmp_ps(a.m, b.m, _CMP_GE_OQ));
}

/**
 * @brief Return the min vector of two vectors.
 *
 * If either lane value is NaN, @c b will be returned for that lane.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vfloat8 min(vfloat8 a, vfloat8 b)
{
	return vfloat8(_mm256_min_ps(a.m, b.m));
}

/**
 * @brief Return the max vector of two vectors.
 *
 * If either lane value is NaN, @c b will be returned for that lane.
 *
 * @param a   The first operand.
 * @param b   The second operand.
 */
ASTCENC_SIMD_INLINE vfloat8 max(vfloat8 a, vfloat8 b)
{
	return vfloat8(_mm256_max_ps(a.m, b.m));
}

/**
 * @brief Return the absolute value of the float vector.
 *
 * @param a   The values whose absolute values are required.
 */
ASTCENC_SIMD_INLINE vfloat8 abs(vfloat8 a)
{
	__m256 msk = _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff));
	return vfloat8(_mm256_and_ps(a.m, msk));
}

/**
 * @brief Return a float rounded to the nearest integer value.
 *
 * @param a   The values to round.
 */
ASTCENC_SIMD_INLINE vfloat8 round(vfloat8 a)
{
	constexpr int flags = _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC;
	return vfloat8(_mm256_round_ps(a.m, flags));
}

/**
 * @brief Return the horizontal minimum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE vfloat8 hmin(vfloat8 a)
{
	__m128 vlow = _mm256_castps256_ps128(a.m);
	__m128 vhigh = _mm256_extractf128_ps(a.m, 1);
	vlow = _mm_min_ps(vlow, vhigh);

	// First do an horizontal reduction
	__m128 shuf = _mm_shuffle_ps(vlow, vlow, _MM_SHUFFLE(2, 3, 0, 1));
	__m128 mins = _mm_min_ps(vlow, shuf);
	shuf = _mm_movehl_ps(shuf, mins);
	mins = _mm_min_ss(mins, shuf);

	// Use insertf128 for compatibility; _mm256_set_m128 requires GCC 9 or Clang 9
	__m256 r = _mm256_insertf128_ps(_mm256_castps128_ps256(mins), mins, 1);

	return vfloat8(_mm256_permute_ps(r, 0));
}

/**
 * @brief Return the horizontal minimum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE float hmin_s(vfloat8 a)
{
	return _mm256_cvtss_f32(hmin(a).m);
}

/**
 * @brief Return the horizontal maximum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE vfloat8 hmax(vfloat8 a)
{
	__m128 vlow = _mm256_castps256_ps128(a.m);
	__m128 vhigh = _mm256_extractf128_ps(a.m, 1);
	vhigh = _mm_max_ps(vlow, vhigh);

	// First do an horizontal reduction
	__m128 shuf = _mm_shuffle_ps(vhigh, vhigh, _MM_SHUFFLE(2, 3, 0, 1));
	__m128 maxs = _mm_max_ps(vhigh, shuf);
	shuf = _mm_movehl_ps(shuf,maxs);
	maxs = _mm_max_ss(maxs, shuf);

	// Use insertf128 for compatibility; _mm256_set_m128 requires GCC 9 or Clang 9
	__m256 r = _mm256_insertf128_ps(_mm256_castps128_ps256(maxs), maxs, 1);
	return vfloat8(_mm256_permute_ps(r, 0));
}

/**
 * @brief Return the horizontal maximum of a vector.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE float hmax_s(vfloat8 a)
{
	return _mm256_cvtss_f32(hmax(a).m);
}

/**
 * @brief Return the horizontal sum of a vector as a scalar.
 *
 * @param a   The vector to reduce.
 */
ASTCENC_SIMD_INLINE float hadd_s(vfloat8 a)
{
	// Two sequential 4-wide adds gives invariance with 4-wide code
	vfloat4 lo(_mm256_extractf128_ps(a.m, 0));
	vfloat4 hi(_mm256_extractf128_ps(a.m, 1));
	return hadd_s(lo) + hadd_s(hi);
}

/**
 * @brief Return the sqrt of the lanes in the vector.
 *
 * @param a   The values whose square roots are required.
 */
ASTCENC_SIMD_INLINE vfloat8 sqrt(vfloat8 a)
{
	return vfloat8(_mm256_sqrt_ps(a.m));
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
	return vfloat8(_mm256_blendv_ps(a.m, b.m, cond.m));
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
	vfloat4 lo(_mm256_extractf128_ps(a.m, 0));
	haccumulate(accum, lo);

	vfloat4 hi(_mm256_extractf128_ps(a.m, 1));
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
	a = select(vfloat8::zero(), a, m);
	accum = accum + a;
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
	return vint8(_mm256_i32gather_epi32(base, indices.m, 4));
}

/**
 * @brief Load a vector of gathered results from an array.
 *
 * @param base      The source float table.
 * @param indices   The table index for each lane.
 */
ASTCENC_SIMD_INLINE vfloat8 gatherf(const float* base, vint8 indices)
{
	return vfloat8(_mm256_i32gather_ps(base, indices.m, 4));
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
#if ASTCENC_X86_GATHERS == 0
	// Gather with scalar loads in two dependency chains, merging late. MSVC preserves this layout;
	// Clang instead uses memory-operand inserts into 128-bit halves, which benchmark slower.
	__m256 m0 = _mm256_broadcast_ss(base + indices[0]);
	__m256 m1 = _mm256_broadcast_ss(base + indices[1]);
	m0 = _mm256_blend_ps(m0, _mm256_broadcast_ss(base + indices[2]), 1 << 2);
	m1 = _mm256_blend_ps(m1, _mm256_broadcast_ss(base + indices[3]), 1 << 3);
	m0 = _mm256_blend_ps(m0, _mm256_broadcast_ss(base + indices[4]), 1 << 4);
	m1 = _mm256_blend_ps(m1, _mm256_broadcast_ss(base + indices[5]), 1 << 5);
	m0 = _mm256_blend_ps(m0, _mm256_broadcast_ss(base + indices[6]), 1 << 6);
	m1 = _mm256_blend_ps(m1, _mm256_broadcast_ss(base + indices[7]), 1 << 7);

	return vfloat8(_mm256_blend_ps(m0, m1, 0xaa));
#else
	vint8 inds(indices);
	return gatherf(base, inds);
#endif
}

/**
 * @brief Store a vector to an unaligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The destination address.
 */
ASTCENC_SIMD_INLINE void store(vfloat8 a, float* p)
{
	_mm256_storeu_ps(p, a.m);
}

/**
 * @brief Store a vector to a 32B aligned memory address.
 *
 * @param      a   The vector to store.
 * @param[out] p   The aligned destination address.
 */
ASTCENC_SIMD_INLINE void storea(vfloat8 a, float* p)
{
	_mm256_store_ps(p, a.m);
}

/**
 * @brief Return an integer value for a float vector, using truncation.
 *
 * @param a   The floating-point values to convert.
 */
ASTCENC_SIMD_INLINE vint8 float_to_int(vfloat8 a)
{
	return vint8(_mm256_cvttps_epi32(a.m));
}

/**
 * @brief Return a float value for an integer vector.
 *
 * @param a   The integer values to convert.
 */
ASTCENC_SIMD_INLINE vfloat8 int_to_float(vint8 a)
{
	return vfloat8(_mm256_cvtepi32_ps(a.m));
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
	return vint8(_mm256_castps_si256(a.m));
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
	return vfloat8(_mm256_castsi256_ps(a.m));
}

/** @brief Prepared lookup table for 16x 8-bit entries. */
struct vtable8_16x8 {
	/** @brief The table entries, replicated in each 128-bit half. */
	vint8 t0;
};

/** @brief Prepared lookup table for 32x 8-bit entries. */
struct vtable8_32x8 {
	/** @brief The table entries, or the replicated first 16 entries for AVX2. */
	vint8 t0;
#if ASTCENC_AVX < 3
	/** @brief The second 16 entries, XOR-chained with the first. */
	vint8 t1;
#endif
};

/** @brief Prepared lookup table for 64x 8-bit entries. */
struct vtable8_64x8 {
	/** @brief The first 32 entries, or replicated first 16 entries for AVX2. */
	vint8 t0;
	/** @brief The second 32 entries, or the second XOR-chained AVX2 row. */
	vint8 t1;
#if ASTCENC_AVX < 3
	/** @brief The third XOR-chained AVX2 row. */
	vint8 t2;
	/** @brief The fourth XOR-chained AVX2 row. */
	vint8 t3;
#endif
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
	// AVX2 tables duplicate table entries in each 128-bit half-register
	vint4 d0 = vint4::load(data);

	table.t0 = vint8(astcenc_mm256_set_m128i(d0.m, d0.m));
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
#if ASTCENC_AVX >= 3
	// VBMI permutes can index across the full 256-bit register
	table.t0 = vint8::load(data);
#else
	// AVX2 tables duplicate table entries in each 128-bit half-register
	vint4 d0 = vint4::load(data);
	vint4 d1 = vint4::load(data + 16);

	table.t0 = vint8(astcenc_mm256_set_m128i(d0.m, d0.m));
	table.t1 = vint8(astcenc_mm256_set_m128i(d1.m, d1.m));

	// XOR chain the high rows to allow table emulation
	table.t1 = table.t1 ^ table.t0;
#endif
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
#if ASTCENC_AVX >= 3
	// Two full-width registers hold all 64 table entries
	table.t0 = vint8::load(data);
	table.t1 = vint8::load(data + 32);
#else
	// AVX2 tables duplicate table entries in each 128-bit half-register
	vint4 d0 = vint4::load(data);
	vint4 d1 = vint4::load(data + 16);
	vint4 d2 = vint4::load(data + 32);
	vint4 d3 = vint4::load(data + 48);

	table.t0 = vint8(astcenc_mm256_set_m128i(d0.m, d0.m));
	table.t1 = vint8(astcenc_mm256_set_m128i(d1.m, d1.m));
	table.t2 = vint8(astcenc_mm256_set_m128i(d2.m, d2.m));
	table.t3 = vint8(astcenc_mm256_set_m128i(d3.m, d3.m));

	// XOR chain the high rows to allow table emulation
	table.t3 = table.t3 ^ table.t2;
	table.t2 = table.t2 ^ table.t1;
	table.t1 = table.t1 ^ table.t0;
#endif
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
	// Set index byte MSB to 1 for unused bytes so shuffle returns zero
	__m256i idxx = _mm256_or_si256(idx.m, _mm256_set1_epi32(static_cast<int>(0xFFFFFF00)));

	__m256i result = _mm256_shuffle_epi8(tbl.t0.m, idxx);
	return vint8(result);
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
#if ASTCENC_AVX >= 3
	// Select one byte per dword and zero the three unused result bytes
	return vint8(_mm256_maskz_permutexvar_epi8(0x11111111, idx.m, tbl.t0.m));
#else
	// Set index byte MSB to 1 for unused bytes so shuffle returns zero
	__m256i idxx = _mm256_or_si256(idx.m, _mm256_set1_epi32(static_cast<int>(0xFFFFFF00)));

	__m256i result = _mm256_shuffle_epi8(tbl.t0.m, idxx);
	idxx = _mm256_sub_epi8(idxx, _mm256_set1_epi8(16));

	__m256i result2 = _mm256_shuffle_epi8(tbl.t1.m, idxx);
	result = _mm256_xor_si256(result, result2);
	return vint8(result);
#endif
}

/**
 * @brief Perform a vtable lookup in a 64x 8-bit table with 32-bit indices.
 *
 * @param tbl   The prepared lookup table.
 * @param idx   The table index for each lane.
 */
ASTCENC_SIMD_INLINE vint8 vtable_lookup_32bit(
	const vtable8_64x8& tbl,
	vint8 idx
) {
#if ASTCENC_AVX >= 3
	// Index bit 5 selects between the two 32-byte table registers
	return vint8(_mm256_maskz_permutex2var_epi8(0x11111111, tbl.t0.m, idx.m, tbl.t1.m));
#else
	// Set index byte MSB to 1 for unused bytes so shuffle returns zero
	__m256i idxx = _mm256_or_si256(idx.m, _mm256_set1_epi32(static_cast<int>(0xFFFFFF00)));

	__m256i result = _mm256_shuffle_epi8(tbl.t0.m, idxx);
	idxx = _mm256_sub_epi8(idxx, _mm256_set1_epi8(16));

	__m256i result2 = _mm256_shuffle_epi8(tbl.t1.m, idxx);
	result = _mm256_xor_si256(result, result2);
	idxx = _mm256_sub_epi8(idxx, _mm256_set1_epi8(16));

	result2 = _mm256_shuffle_epi8(tbl.t2.m, idxx);
	result = _mm256_xor_si256(result, result2);
	idxx = _mm256_sub_epi8(idxx, _mm256_set1_epi8(16));

	result2 = _mm256_shuffle_epi8(tbl.t3.m, idxx);
	result = _mm256_xor_si256(result, result2);

	return vint8(result);
#endif
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
	_mm256_maskstore_epi32(reinterpret_cast<int*>(base), _mm256_castps_si256(mask.m), data.m);
}

#endif  // #ifndef ASTC_VECMATHLIB_AVX2_8_H_INCLUDED
