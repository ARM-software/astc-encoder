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
 * @brief 4x32-bit vectors, implemented using SSE.
 *
 * This module implements 4-wide 32-bit float, int, and mask vectors for x86
 * SSE. The implementation requires at least SSE2, but higher levels of SSE can
 * be selected at compile time to improve performance.
 *
 * There is a baseline level of functionality provided by all vector widths and
 * implementations. This is implemented using identical function signatures,
 * modulo data type, so we can use them as substitutable implementations in VLA
 * code.
 *
 * The 4-wide vectors are also used as a fixed-width type, and significantly
 * extend the functionality above that available to VLA code.
 */

#ifndef ASTC_VECMATHLIB_SSE_4_H_INCLUDED
#define ASTC_VECMATHLIB_SSE_4_H_INCLUDED

#ifndef ASTCENC_SIMD_INLINE
	#error "Include astcenc_vecmathlib.h, do not include directly"
#endif

#include <cstdio>

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
	 * Consider using loada() which is better with vectors if data is aligned
	 * to vector length.
	 */
	ASTCENC_SIMD_INLINE explicit vfloat4(const float *p)
	{
		m = _mm_loadu_ps(p);
	}

	/**
	 * @brief Construct from 1 scalar value replicated across all lanes.
	 *
	 * Consider using zero() for constexpr zeros.
	 */
	ASTCENC_SIMD_INLINE explicit vfloat4(float a)
	{
		m = _mm_set1_ps(a);
	}

	/**
	 * @brief Construct from 4 scalar values.
	 *
	 * The value of @c a is stored to lane 0 (LSB) in the SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vfloat4(float a, float b, float c, float d)
	{
		m = _mm_set_ps(d, c, b, a);
	}

	/**
	 * @brief Construct from an existing SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vfloat4(__m128 a)
	{
		m = a;
	}

	/**
	 * @brief Get the scalar value of a single lane.
	 *
	 * TODO: Can we do better for lane0, which is the common case for VLA?
	 */
	template <int l> ASTCENC_SIMD_INLINE float lane() const
	{
		return _mm_cvtss_f32(_mm_shuffle_ps(m, m, l));
	}

	/**
	 * @brief Set the scalar value of a single lane.
	 */
	template <int l> ASTCENC_SIMD_INLINE void set_lane(float a)
	{
#if ASTCENC_SSE >= 41
		__m128 v = _mm_set1_ps(a);
		m = _mm_insert_ps(m, v, l << 6 | l << 4);
#else
		alignas(16) float idx[4];
		_mm_store_ps(idx, m);
		idx[l] = a;
		m = _mm_load_ps(idx);
#endif
	}

	/**
	 * @brief Factory that returns a vector of zeros.
	 */
	static ASTCENC_SIMD_INLINE vfloat4 zero()
	{
		return vfloat4(_mm_setzero_ps());
	}

	/**
	 * @brief Factory that returns a replicated scalar loaded from memory.
	 */
	static ASTCENC_SIMD_INLINE vfloat4 load1(const float* p)
	{
		return vfloat4(_mm_load_ps1(p));
	}

	/**
	 * @brief Factory that returns a vector loaded from 16B aligned memory.
	 */
	static ASTCENC_SIMD_INLINE vfloat4 loada(const float* p)
	{
		return vfloat4(_mm_load_ps(p));
	}

	/**
	 * @brief Factory that returns a vector containing the lane IDs.
	 */
	static ASTCENC_SIMD_INLINE vfloat4 lane_id()
	{
		return vfloat4(_mm_set_ps(3, 2, 1, 0));
	}

	/**
	 * @brief The vector ...
	 */
	__m128 m;
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
	 * Consider using loada() which is better with vectors if data is aligned
	 * to vector length.
	 */
	ASTCENC_SIMD_INLINE explicit vint4(const int *p)
	{
		m = _mm_loadu_si128((const __m128i*)p);
	}

	/**
	 * @brief Construct from 1 scalar value replicated across all lanes.
	 *
	 * Consider using vfloat4::zero() for constexpr zeros.
	 */
	ASTCENC_SIMD_INLINE explicit vint4(int a)
	{
		m = _mm_set1_epi32(a);
	}

	/**
	 * @brief Construct from 4 scalar values.
	 *
	 * The value of @c a is stored to lane 0 (LSB) in the SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vint4(int a, int b, int c, int d)
	{
		m = _mm_set_epi32(d, c, b, a);
	}

	/**
	 * @brief Construct from an existing SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vint4(__m128i a)
	{
		m = a;
	}

	/**
	 * @brief Get the scalar from a single lane.
	 */
	template <int l> ASTCENC_SIMD_INLINE int lane() const
	{
		return _mm_cvtsi128_si32(_mm_shuffle_epi32(m, l));
	}

	/**
	 * @brief Factory that returns a vector containing the lane IDs.
	 */
	static ASTCENC_SIMD_INLINE vint4 lane_id()
	{
		return vint4(_mm_set_epi32(3, 2, 1, 0));
	}

	/**
	 * @brief The vector ...
	 */
	__m128i m;
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
	 * @brief Construct from an existing SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vmask4(__m128 a)
	{
		m = a;
	}

	/**
	 * @brief Construct from an existing SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vmask4(__m128i a)
	{
		m = _mm_castsi128_ps(a);
	}

	/**
	 * @brief The vector ...
	 */
	__m128 m;
};

// ============================================================================
// vmask4 operators and functions
// ============================================================================

/**
 * @brief Overload: mask union (or).
 */
ASTCENC_SIMD_INLINE vmask4 operator|(vmask4 a, vmask4 b)
{
	return vmask4(_mm_or_ps(a.m, b.m));
}

/**
 * @brief Overload: mask intersect (and).
 */
ASTCENC_SIMD_INLINE vmask4 operator&(vmask4 a, vmask4 b)
{
	return vmask4(_mm_and_ps(a.m, b.m));
}

/**
 * @brief Overload: mask difference (xor).
 */
ASTCENC_SIMD_INLINE vmask4 operator^(vmask4 a, vmask4 b)
{
	return vmask4(_mm_xor_ps(a.m, b.m));
}

/**
 * @brief Overload: mask invert (not).
 */
ASTCENC_SIMD_INLINE vmask4 operator~(vmask4 a)
{
	return vmask4(_mm_xor_si128(_mm_castps_si128(a.m), _mm_set1_epi32(-1)));
}

/**
 * @brief Return a 4-bit mask code indicating mask status.
 *
 * bit0 = lane 0
 */
ASTCENC_SIMD_INLINE unsigned int mask(vmask4 a)
{
	return _mm_movemask_ps(a.m);
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
	return vint4(_mm_add_epi32(a.m, b.m));
}

/**
 * @brief Overload: vector by vector subtraction.
 */
ASTCENC_SIMD_INLINE vint4 operator-(vint4 a, vint4 b)
{
	return vint4(_mm_sub_epi32(a.m, b.m));
}

/**
 * @brief Overload: vector bit invert.
 */
ASTCENC_SIMD_INLINE vint4 operator~(vint4 a)
{
	return vint4(_mm_xor_si128(a.m, _mm_set1_epi32(-1)));
}

/**
 * @brief Overload: vector by vector bitwise or.
 */
ASTCENC_SIMD_INLINE vint4 operator|(vint4 a, vint4 b)
{
	return vint4(_mm_or_si128(a.m, b.m));
}

/**
 * @brief Overload: vector by vector bitwise and.
 */
ASTCENC_SIMD_INLINE vint4 operator&(vint4 a, vint4 b)
{
	return vint4(_mm_and_si128(a.m, b.m));
}

/**
 * @brief Overload: vector by vector bitwise xor.
 */
ASTCENC_SIMD_INLINE vint4 operator^(vint4 a, vint4 b)
{
	return vint4(_mm_xor_si128(a.m, b.m));
}

/**
 * @brief Overload: vector by vector equality.
 */
ASTCENC_SIMD_INLINE vmask4 operator==(vint4 a, vint4 b)
{
	return vmask4(_mm_cmpeq_epi32(a.m, b.m));
}

/**
 * @brief Overload: vector by vector inequality.
 */
ASTCENC_SIMD_INLINE vmask4 operator!=(vint4 a, vint4 b)
{
	return ~vmask4(_mm_cmpeq_epi32(a.m, b.m));
}

/**
 * @brief Overload: vector by vector less than.
 */
ASTCENC_SIMD_INLINE vmask4 operator<(vint4 a, vint4 b)
{
	return vmask4(_mm_cmplt_epi32(a.m, b.m));
}

/**
 * @brief Overload: vector by vector greater than.
 */
ASTCENC_SIMD_INLINE vmask4 operator>(vint4 a, vint4 b)
{
	return vmask4(_mm_cmpgt_epi32(a.m, b.m));
}

/**
 * @brief Return the min vector of two vectors.
 */
ASTCENC_SIMD_INLINE vint4 min(vint4 a, vint4 b)
{
#if ASTCENC_SSE >= 41
	return vint4(_mm_min_epi32(a.m, b.m));
#else
	vmask4 d = a < b;
	__m128i ap = _mm_and_si128(_mm_castps_si128(d.m), a.m);
	__m128i bp = _mm_andnot_si128(_mm_castps_si128(d.m), b.m);
	return vint4(_mm_or_si128(ap,bp));
#endif
}

/**
 * @brief Return the max vector of two vectors.
 */
ASTCENC_SIMD_INLINE vint4 max(vint4 a, vint4 b)
{
#if ASTCENC_SSE >= 41
	return vint4(_mm_max_epi32(a.m, b.m));
#else
	vmask4 d = a > b;
	__m128i ap = _mm_and_si128(_mm_castps_si128(d.m), a.m);
	__m128i bp = _mm_andnot_si128(_mm_castps_si128(d.m), b.m);
	return vint4(_mm_or_si128(ap,bp));
#endif
}

/**
 * @brief Return the horizontal minimum of a vector.
 */
ASTCENC_SIMD_INLINE vint4 hmin(vint4 a)
{
	a = min(a, vint4(_mm_shuffle_epi32(a.m, _MM_SHUFFLE(0, 0, 3, 2))));
	a = min(a, vint4(_mm_shuffle_epi32(a.m, _MM_SHUFFLE(0, 0, 0, 1))));
	return vint4(_mm_shuffle_epi32(a.m, _MM_SHUFFLE(0, 0, 0, 0)));
}

/**
 * @brief Store a vector to a 16B aligned memory address.
 */
ASTCENC_SIMD_INLINE void storea(vint4 a, int* p)
{
	_mm_store_si128((__m128i*)p, a.m);
}

/**
 * @brief Store lowest N (vector width) bytes into an unaligned address.
 */
ASTCENC_SIMD_INLINE void store_nbytes(vint4 a, uint8_t* p)
{
	// This is the most logical implementation, but the convenience intrinsic
	// is missing on older compilers (supported in g++ 9 and clang++ 9).
	// _mm_storeu_si32(ptr, v.m);
	_mm_store_ss((float*)p, _mm_castsi128_ps(a.m));
}

/**
 * @brief Gather N (vector width) indices from the array.
 */
ASTCENC_SIMD_INLINE vint4 gatheri(const int* base, vint4 indices)
{
	alignas(16) int idx[4];
	storea(indices, idx);
	return vint4(base[idx[0]], base[idx[1]], base[idx[2]], base[idx[3]]);
}

/**
 * @brief Pack low 8 bits of N (vector width) lanes into bottom of vector.
 */
ASTCENC_SIMD_INLINE vint4 pack_low_bytes(vint4 a)
{
#if ASTCENC_SSE >= 41
	__m128i shuf = _mm_set_epi8(0,0,0,0, 0,0,0,0, 0,0,0,0, 12,8,4,0);
	return vint4(_mm_shuffle_epi8(a.m, shuf));
#else
	__m128i va = _mm_unpacklo_epi8(a.m, _mm_shuffle_epi32(a.m, _MM_SHUFFLE(1,1,1,1)));
	__m128i vb = _mm_unpackhi_epi8(a.m, _mm_shuffle_epi32(a.m, _MM_SHUFFLE(3,3,3,3)));
	return vint4(_mm_unpacklo_epi16(va, vb));
#endif
}

/**
 * @brief Return lanes from @c b if MSB of @c cond is set, else @c a.
 */
ASTCENC_SIMD_INLINE vint4 select(vint4 a, vint4 b, vmask4 cond)
{
#if ASTCENC_SSE >= 41
	// Don't use _mm_blendv_epi8 directly, as it doesn't give the select on
	// float sign-bit in the mask behavior which is useful. Performance is the
	// same, these casts are free.
	__m128 av = _mm_castsi128_ps(a.m);
	__m128 bv = _mm_castsi128_ps(b.m);
	return vint4(_mm_castps_si128(_mm_blendv_ps(av, bv, cond.m)));
#else
	__m128i d = _mm_srai_epi32(_mm_castps_si128(cond.m), 31);
	return vint4(_mm_or_si128(_mm_and_si128(d, b.m), _mm_andnot_si128(d, a.m)));
#endif
}

/**
 * @brief Debug function to print a vector of ints.
 */
ASTCENC_SIMD_INLINE void print(vint4 a)
{
	alignas(16) int v[4];
	storea(a, v);
	printf("v4_i32:\n  %8u %8u %8u %8u\n",
	       v[0], v[1], v[2], v[3]);
}

// ============================================================================
// vfloat4 operators and functions
// ============================================================================

/**
 * @brief Overload: vector by vector addition.
 */
ASTCENC_SIMD_INLINE vfloat4 operator+(vfloat4 a, vfloat4 b)
{
	return vfloat4(_mm_add_ps(a.m, b.m));
}

/**
 * @brief Overload: vector by vector subtraction.
 */
ASTCENC_SIMD_INLINE vfloat4 operator-(vfloat4 a, vfloat4 b)
{
	return vfloat4(_mm_sub_ps(a.m, b.m));
}

/**
 * @brief Overload: vector by vector multiplication.
 */
ASTCENC_SIMD_INLINE vfloat4 operator*(vfloat4 a, vfloat4 b)
{
	return vfloat4(_mm_mul_ps(a.m, b.m));
}

/**
 * @brief Overload: vector by scalar multiplication.
 */
ASTCENC_SIMD_INLINE vfloat4 operator*(vfloat4 a, float b)
{
	return vfloat4(_mm_mul_ps(a.m, _mm_set1_ps(b)));
}

/**
 * @brief Overload: scalar by vector multiplication.
 */
ASTCENC_SIMD_INLINE vfloat4 operator*(float a, vfloat4 b)
{
	return vfloat4(_mm_mul_ps(_mm_set1_ps(a), b.m));
}

/**
 * @brief Overload: vector by vector division.
 */
ASTCENC_SIMD_INLINE vfloat4 operator/(vfloat4 a, vfloat4 b)
{
	return vfloat4(_mm_div_ps(a.m, b.m));
}

/**
 * @brief Overload: vector by vector equality.
 */
ASTCENC_SIMD_INLINE vmask4 operator==(vfloat4 a, vfloat4 b)
{
	return vmask4(_mm_cmpeq_ps(a.m, b.m));
}

/**
 * @brief Overload: vector by vector inequality.
 */
ASTCENC_SIMD_INLINE vmask4 operator!=(vfloat4 a, vfloat4 b)
{
	return vmask4(_mm_cmpneq_ps(a.m, b.m));
}

/**
 * @brief Overload: vector by vector less than.
 */
ASTCENC_SIMD_INLINE vmask4 operator<(vfloat4 a, vfloat4 b)
{
	return vmask4(_mm_cmplt_ps(a.m, b.m));
}

/**
 * @brief Overload: vector by vector greater than.
 */
ASTCENC_SIMD_INLINE vmask4 operator>(vfloat4 a, vfloat4 b)
{
	return vmask4(_mm_cmpgt_ps(a.m, b.m));
}

/**
 * @brief Overload: vector by vector less than or equal.
 */
ASTCENC_SIMD_INLINE vmask4 operator<=(vfloat4 a, vfloat4 b)
{
	return vmask4(_mm_cmple_ps(a.m, b.m));
}

/**
 * @brief Overload: vector by vector greater than or equal.
 */
ASTCENC_SIMD_INLINE vmask4 operator>=(vfloat4 a, vfloat4 b)
{
	return vmask4(_mm_cmpge_ps(a.m, b.m));
}

/**
 * @brief Return the min vector of two vectors.
 *
 * If either lane value is NaN, @c b will be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat4 min(vfloat4 a, vfloat4 b)
{
	// Do not reorder - second operand will return if either is NaN
	return vfloat4(_mm_min_ps(a.m, b.m));
}

/**
 * @brief Return the max vector of two vectors.
 *
 * If either lane value is NaN, @c b will be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat4 max(vfloat4 a, vfloat4 b)
{
	// Do not reorder - second operand will return if either is NaN
	return vfloat4(_mm_max_ps(a.m, b.m));
}

/**
 * @brief Return the clamped value between min and max.
 *
 * It is assumed that neither @c min nor @c max are NaN values. If @c a is NaN
 * then @c min will be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat4 clamp(float min, float max, vfloat4 a)
{
	// Do not reorder - second operand will return if either is NaN
	a.m = _mm_max_ps(a.m, _mm_set1_ps(min));
	a.m = _mm_min_ps(a.m, _mm_set1_ps(max));
	return a;
}

/**
 * @brief Return a clamped value between 0.0f and max.
 *
 * It is assumed that @c max is not a NaN value. If @c a is NaN then zero will
 * be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat4 clampz(float max, vfloat4 a)
{
	// Do not reorder - second operand will return if either is NaN
	a.m = _mm_max_ps(a.m, _mm_setzero_ps());
	a.m = _mm_min_ps(a.m, _mm_set1_ps(max));
	return a;
}

/**
 * @brief Return a clamped value between 0.0f and 1.0f.
 *
 * If @c a is NaN then zero will be returned for that lane.
 */
ASTCENC_SIMD_INLINE vfloat4 clampzo(vfloat4 a)
{
	a.m = _mm_max_ps(a.m, _mm_setzero_ps());
	a.m = _mm_min_ps(a.m, _mm_set1_ps(1.0f));
	return a;
}

/**
 * @brief Return the absolute value of the float vector.
 */
ASTCENC_SIMD_INLINE vfloat4 abs(vfloat4 a)
{
	return vfloat4(_mm_max_ps(_mm_sub_ps(_mm_setzero_ps(), a.m), a.m));
}

/**
 * @brief Return a float rounded to the nearest integer value.
 *
 * TODO: Can we do a better fallback here, if we exploit the fact that we
 *       can assume that values are positive?
 */
ASTCENC_SIMD_INLINE vfloat4 round(vfloat4 a)
{
#if ASTCENC_SSE >= 41
	constexpr int flags = _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC;
	return vfloat4(_mm_round_ps(a.m, flags));
#else
	__m128 v = a.m;
	__m128 neg_zero = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
	__m128 no_fraction = _mm_set1_ps(8388608.0f);
	__m128 abs_mask = _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF));
	__m128 sign = _mm_and_ps(v, neg_zero);
	__m128 s_magic = _mm_or_ps(no_fraction, sign);
	__m128 r1 = _mm_add_ps(v, s_magic);
	r1 = _mm_sub_ps(r1, s_magic);
	__m128 r2 = _mm_and_ps(v, abs_mask);
	__m128 mask = _mm_cmple_ps(r2, no_fraction);
	r2 = _mm_andnot_ps(mask, v);
	r1 = _mm_and_ps(r1, mask);
	return vfloat4(_mm_xor_ps(r1, r2));
#endif
}

/**
 * @brief Return the horizontal minimum of a vector.
 */
ASTCENC_SIMD_INLINE vfloat4 hmin(vfloat4 a)
{
	a = min(a, vfloat4(_mm_shuffle_ps(a.m, a.m, _MM_SHUFFLE(0, 0, 3, 2))));
	a = min(a, vfloat4(_mm_shuffle_ps(a.m, a.m, _MM_SHUFFLE(0, 0, 0, 1))));
	return vfloat4(_mm_shuffle_ps(a.m, a.m, _MM_SHUFFLE(0, 0, 0, 0)));
}

/**
 * @brief Return the sqrt of the lanes in the vector.
 */
ASTCENC_SIMD_INLINE vfloat4 sqrt(vfloat4 a)
{
	return vfloat4(_mm_sqrt_ps(a.m));
}

/**
 * @brief Return lanes from @c b if MSB of @c cond is set, else @c a.
 */
ASTCENC_SIMD_INLINE vfloat4 select(vfloat4 a, vfloat4 b, vmask4 cond)
{
#if ASTCENC_SSE >= 41
	return vfloat4(_mm_blendv_ps(a.m, b.m, cond.m));
#else
	__m128 d = _mm_castsi128_ps(_mm_srai_epi32(_mm_castps_si128(cond.m), 31));
	return vfloat4(_mm_or_ps(_mm_and_ps(d, b.m), _mm_andnot_ps(d, a.m)));
#endif
}

/**
 * @brief Load a vector of gathered results from an array;
 */
ASTCENC_SIMD_INLINE vfloat4 gatherf(const float* base, vint4 indices)
{
	alignas(16) int idx[4];
	storea(indices, idx);
	return vfloat4(base[idx[0]], base[idx[1]], base[idx[2]], base[idx[3]]);
}

/**
 * @brief Store a vector to an unaligned memory address.
 */
ASTCENC_SIMD_INLINE void store(vfloat4 a, float* p)
{
	_mm_storeu_ps(p, a.m);
}

/**
 * @brief Store a vector to a 16B aligned memory address.
 */
ASTCENC_SIMD_INLINE void storea(vfloat4 a, float* p)
{
	_mm_store_ps(p, a.m);
}

/**
 * @brief Return the dot product for the full 4 lanes, returning vector.
 */
ASTCENC_SIMD_INLINE vfloat4 dot(vfloat4 a, vfloat4 b)
{
#if (ASTCENC_SSE >= 41) && (ASTCENC_ISA_INVARIANCE == 0)
	return vfloat4(_mm_dp_ps(a.m, b.m, 0xFF));
#else
	alignas(16) float av[4];
	alignas(16) float bv[4];
	storea(a, av);
	storea(b, bv);
	float s = av[0] * bv[0] + av[1] * bv[1] + av[2] * bv[2] + av[3] * bv[3];
	return vfloat4(s);
#endif
}

/**
 * @brief Return a integer value for a float vector, using truncation.
 */
ASTCENC_SIMD_INLINE vint4 float_to_int(vfloat4 a)
{
	return vint4(_mm_cvttps_epi32(a.m));
}

/**
 * @brief Return a integer value for a float vector, using round-to-nearest.
 */
ASTCENC_SIMD_INLINE vint4 float_to_int_rtn(vfloat4 a)
{
	a = round(a);
	return vint4(_mm_cvttps_epi32(a.m));
}

/**
 * @brief Return a float value as an integer bit pattern (i.e. no conversion).
 *
 * It is a common trick to convert floats into integer bit patterns, perform
 * some bit hackery based on knowledge they are IEEE 754 layout, and then
 * convert them back again. This is the first half of that flip.
 */
ASTCENC_SIMD_INLINE vint4 float_as_int(vfloat4 a)
{
	return vint4(_mm_castps_si128(a.m));
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
	return vfloat4(_mm_castsi128_ps(v.m));
}

/**
 * @brief Debug function to print a vector of floats.
 */
ASTCENC_SIMD_INLINE void print(vfloat4 a)
{
	alignas(16) float v[4];
	storea(a, v);
	printf("v4_f32:\n  %0.4f %0.4f %0.4f %0.4f\n",
	       (double)v[0], (double)v[1], (double)v[2], (double)v[3]);
}

#endif // #ifndef ASTC_VECMATHLIB_SSE_4_H_INCLUDED
