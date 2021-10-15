// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2019-2021 Arm Limited
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
	ASTCENC_SIMD_INLINE vfloat4() = default;

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
	 * @brief Return a swizzled float 2.
	 */
	template <int l0, int l1> ASTCENC_SIMD_INLINE vfloat4 swz() const
	{
		vfloat4 result(_mm_shuffle_ps(m, m, l0 | l1 << 2));
		result.set_lane<2>(0.0f);
		result.set_lane<3>(0.0f);
		return result;
	}

	/**
	 * @brief Return a swizzled float 3.
	 */
	template <int l0, int l1, int l2> ASTCENC_SIMD_INLINE vfloat4 swz() const
	{
		vfloat4 result(_mm_shuffle_ps(m, m, l0 | l1 << 2 | l2 << 4));
		result.set_lane<3>(0.0f);
		return result;
	}

	/**
	 * @brief Return a swizzled float 4.
	 */
	template <int l0, int l1, int l2, int l3> ASTCENC_SIMD_INLINE vfloat4 swz() const
	{
		return vfloat4(_mm_shuffle_ps(m, m, l0 | l1 << 2 | l2 << 4 | l3 << 6));
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
	ASTCENC_SIMD_INLINE vint4() = default;

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
	 * @brief Construct from 4 uint8_t loaded from an unaligned address.
	 */
	ASTCENC_SIMD_INLINE explicit vint4(const uint8_t *p)
	{
		// _mm_loadu_si32 would be nicer syntax, but missing on older GCC
		__m128i t = _mm_cvtsi32_si128(*(const int*)p);

#if ASTCENC_SSE >= 41
		m = _mm_cvtepu8_epi32(t);
#else
		t = _mm_unpacklo_epi8(t, _mm_setzero_si128());
		m = _mm_unpacklo_epi16(t, _mm_setzero_si128());
#endif
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
	 * @brief Set the scalar value of a single lane.
	 */
	template <int l> ASTCENC_SIMD_INLINE void set_lane(int a)
	{
#if ASTCENC_SSE >= 41
		m = _mm_insert_epi32(m, a, l);
#else
		alignas(16) int idx[4];
		_mm_store_si128((__m128i*)idx, m);
		idx[l] = a;
		m = _mm_load_si128((const __m128i*)idx);
#endif
	}

	/**
	 * @brief Factory that returns a vector of zeros.
	 */
	static ASTCENC_SIMD_INLINE vint4 zero()
	{
		return vint4(_mm_setzero_si128());
	}

	/**
	 * @brief Factory that returns a replicated scalar loaded from memory.
	 */
	static ASTCENC_SIMD_INLINE vint4 load1(const int* p)
	{
		return vint4(*p);
	}

	/**
	 * @brief Factory that returns a vector loaded from 16B aligned memory.
	 */
	static ASTCENC_SIMD_INLINE vint4 loada(const int* p)
	{
		return vint4(_mm_load_si128((const __m128i*)p));
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
	 * @brief Construct from 1 scalar value.
	 */
	ASTCENC_SIMD_INLINE explicit vmask4(bool a)
	{
		vint4 mask(a == false ? 0 : -1);
		m = _mm_castsi128_ps(mask.m);
	}

	/**
	 * @brief Construct from 4 scalar values.
	 *
	 * The value of @c a is stored to lane 0 (LSB) in the SIMD register.
	 */
	ASTCENC_SIMD_INLINE explicit vmask4(bool a, bool b, bool c, bool d)
	{
		vint4 mask(a == false ? 0 : -1,
		           b == false ? 0 : -1,
		           c == false ? 0 : -1,
		           d == false ? 0 : -1);

		m = _mm_castsi128_ps(mask.m);
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
 * @brief Overload: vector by vector multiplication.
 */
ASTCENC_SIMD_INLINE vint4 operator*(vint4 a, vint4 b)
{
#if ASTCENC_SSE >= 41
	return vint4(_mm_mullo_epi32 (a.m, b.m));
#else
	__m128i t1 = _mm_mul_epu32(a.m, b.m);
	__m128i t2 = _mm_mul_epu32(
	                 _mm_srli_si128(a.m, 4),
	                 _mm_srli_si128(b.m, 4));
	__m128i r =  _mm_unpacklo_epi32(
	                 _mm_shuffle_epi32(t1, _MM_SHUFFLE (0, 0, 2, 0)),
	                 _mm_shuffle_epi32(t2, _MM_SHUFFLE (0, 0, 2, 0)));
	return vint4(r);
#endif
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
 * @brief Logical shift left.
 */
template <int s> ASTCENC_SIMD_INLINE vint4 lsl(vint4 a)
{
	return vint4(_mm_slli_epi32(a.m, s));
}

/**
 * @brief Logical shift right.
 */
template <int s> ASTCENC_SIMD_INLINE vint4 lsr(vint4 a)
{
	return vint4(_mm_srli_epi32(a.m, s));
}

/**
 * @brief Arithmetic shift right.
 */
template <int s> ASTCENC_SIMD_INLINE vint4 asr(vint4 a)
{
	return vint4(_mm_srai_epi32(a.m, s));
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

/*
 * @brief Return the horizontal maximum of a vector.
 */
ASTCENC_SIMD_INLINE vint4 hmax(vint4 a)
{
	a = max(a, vint4(_mm_shuffle_epi32(a.m, _MM_SHUFFLE(0, 0, 3, 2))));
	a = max(a, vint4(_mm_shuffle_epi32(a.m, _MM_SHUFFLE(0, 0, 0, 1))));
	return vint4(_mm_shuffle_epi32(a.m, _MM_SHUFFLE(0, 0, 0, 0)));
}

/**
 * @brief Return the horizontal sum of a vector as a scalar.
 */
ASTCENC_SIMD_INLINE int hadd_s(vint4 a)
{
	// Add top and bottom halves, lane 1/0
	__m128i fold = _mm_castps_si128(_mm_movehl_ps(_mm_castsi128_ps(a.m),
	                                              _mm_castsi128_ps(a.m)));
	__m128i t = _mm_add_epi32(a.m, fold);

	// Add top and bottom halves, lane 0 (_mm_hadd_ps exists but slow)
	t = _mm_add_epi32(t, _mm_shuffle_epi32(t, 0x55));

	return _mm_cvtsi128_si32(t);
}

/**
 * @brief Store a vector to a 16B aligned memory address.
 */
ASTCENC_SIMD_INLINE void storea(vint4 a, int* p)
{
	_mm_store_si128((__m128i*)p, a.m);
}

/**
 * @brief Store a vector to an unaligned memory address.
 */
ASTCENC_SIMD_INLINE void store(vint4 a, int* p)
{
	// Cast due to missing intrinsics
	_mm_storeu_ps((float*)p, _mm_castsi128_ps(a.m));
}

/**
 * @brief Store lowest N (vector width) bytes into an unaligned address.
 */
ASTCENC_SIMD_INLINE void store_nbytes(vint4 a, uint8_t* p)
{
	// Cast due to missing intrinsics
	_mm_store_ss((float*)p, _mm_castsi128_ps(a.m));
}

/**
 * @brief Gather N (vector width) indices from the array.
 */
ASTCENC_SIMD_INLINE vint4 gatheri(const int* base, vint4 indices)
{
#if ASTCENC_AVX >= 2
	return vint4(_mm_i32gather_epi32(base, indices.m, 4));
#else
	alignas(16) int idx[4];
	storea(indices, idx);
	return vint4(base[idx[0]], base[idx[1]], base[idx[2]], base[idx[3]]);
#endif
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
 * @brief Return the absolute value of the float vector.
 */
ASTCENC_SIMD_INLINE vfloat4 abs(vfloat4 a)
{
	return vfloat4(_mm_max_ps(_mm_sub_ps(_mm_setzero_ps(), a.m), a.m));
}

/**
 * @brief Return a float rounded to the nearest integer value.
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
 * @brief Return the horizontal maximum of a vector.
 */
ASTCENC_SIMD_INLINE vfloat4 hmax(vfloat4 a)
{
	a = max(a, vfloat4(_mm_shuffle_ps(a.m, a.m, _MM_SHUFFLE(0, 0, 3, 2))));
	a = max(a, vfloat4(_mm_shuffle_ps(a.m, a.m, _MM_SHUFFLE(0, 0, 0, 1))));
	return vfloat4(_mm_shuffle_ps(a.m, a.m, _MM_SHUFFLE(0, 0, 0, 0)));
}

/**
 * @brief Return the horizontal sum of a vector as a scalar.
 */
ASTCENC_SIMD_INLINE float hadd_s(vfloat4 a)
{
	// Add top and bottom halves, lane 1/0
	__m128 t = _mm_add_ps(a.m, _mm_movehl_ps(a.m, a.m));

	// Add top and bottom halves, lane 0 (_mm_hadd_ps exists but slow)
	t = _mm_add_ss(t, _mm_shuffle_ps(t, t, 0x55));

	return _mm_cvtss_f32(t);
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
#if ASTCENC_AVX >= 2
	return vfloat4(_mm_i32gather_ps(base, indices.m, 4));
#else
	alignas(16) int idx[4];
	storea(indices, idx);
	return vfloat4(base[idx[0]], base[idx[1]], base[idx[2]], base[idx[3]]);
#endif
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
 * @brief Return a float value for an integer vector.
 */
ASTCENC_SIMD_INLINE vfloat4 int_to_float(vint4 a)
{
	return vfloat4(_mm_cvtepi32_ps(a.m));
}

/**
 * @brief Return a float16 value for a float vector, using round-to-nearest.
 */
ASTCENC_SIMD_INLINE vint4 float_to_float16(vfloat4 a)
{
#if ASTCENC_F16C >= 1
	__m128i packedf16 = _mm_cvtps_ph(a.m, 0);
	__m128i f16 = _mm_cvtepu16_epi32(packedf16);
	return vint4(f16);
#else
	return vint4(
		float_to_sf16(a.lane<0>()),
		float_to_sf16(a.lane<1>()),
		float_to_sf16(a.lane<2>()),
		float_to_sf16(a.lane<3>()));
#endif
}

/**
 * @brief Return a float16 value for a float scalar, using round-to-nearest.
 */
static inline uint16_t float_to_float16(float a)
{
#if ASTCENC_F16C >= 1
	__m128i f16 = _mm_cvtps_ph(_mm_set1_ps(a), 0);
	return  (uint16_t)_mm_cvtsi128_si32(f16);
#else
	return float_to_sf16(a);
#endif
}

/**
 * @brief Return a float value for a float16 vector.
 */
ASTCENC_SIMD_INLINE vfloat4 float16_to_float(vint4 a)
{
#if ASTCENC_F16C >= 1
	__m128i packed = _mm_packs_epi32(a.m, a.m);
	__m128 f32 = _mm_cvtph_ps(packed);
	return vfloat4(f32);
#else
	return vfloat4(
		sf16_to_float(a.lane<0>()),
		sf16_to_float(a.lane<1>()),
		sf16_to_float(a.lane<2>()),
		sf16_to_float(a.lane<3>()));
#endif
}

/**
 * @brief Return a float value for a float16 scalar.
 */
ASTCENC_SIMD_INLINE float float16_to_float(uint16_t a)
{
#if ASTCENC_F16C >= 1
	__m128i packed = _mm_set1_epi16(a);
	__m128 f32 = _mm_cvtph_ps(packed);
	return _mm_cvtss_f32(f32);
#else
	return sf16_to_float(a);
#endif
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

#endif // #ifndef ASTC_VECMATHLIB_SSE_4_H_INCLUDED
