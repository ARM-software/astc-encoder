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
 * @brief 16x32-bit vectors, implemented using AVX-512F.
 *
 * Experimental compile-time AVX-512 ISA: ASTCENC_SIMD_WIDTH is 16. This is a
 * VLA backend sibling of astcenc_vecmathlib_avx2_8.h, not a runtime overlay.
 * gatherf is implemented with AVX-512 VBMI (vpermb / vpermt2b) for tables that
 * fit in 1–4 zmm, and a scalar-load fallback for larger index ranges. It does
 * not use AVX-512F vgatherdps. vgatherf_load hoists those masked zmm loads so
 * callers can permute many times from one table. vtable_lookup_32bit uses
 * vpermb over packed byte tables instead of AVX2-style lane-local pshufb.
 */

#ifndef ASTC_VECMATHLIB_AVX512_16_H_INCLUDED
#define ASTC_VECMATHLIB_AVX512_16_H_INCLUDED

#ifndef ASTCENC_SIMD_INLINE
	#error "Include astcenc_vecmathlib.h, do not include directly"
#endif

#include <cstdio>
#include <cstring>
#include <limits>

// ============================================================================
// vfloat16 data type
// ============================================================================

struct vfloat16
{
	ASTCENC_SIMD_INLINE vfloat16() = default;

	ASTCENC_SIMD_INLINE explicit vfloat16(const float *p)
	{
		m = _mm512_loadu_ps(p);
	}

	ASTCENC_SIMD_INLINE explicit vfloat16(float a)
	{
		m = _mm512_set1_ps(a);
	}

	ASTCENC_SIMD_INLINE explicit vfloat16(__m512 a)
	{
		m = a;
	}

	static ASTCENC_SIMD_INLINE vfloat16 zero()
	{
		return vfloat16(_mm512_setzero_ps());
	}

	static ASTCENC_SIMD_INLINE vfloat16 load1(const float* p)
	{
		return vfloat16(_mm512_set1_ps(*p));
	}

	static ASTCENC_SIMD_INLINE vfloat16 loada(const float* p)
	{
		return vfloat16(_mm512_load_ps(p));
	}

	__m512 m;
};

// ============================================================================
// vint16 data type
// ============================================================================

struct vint16
{
	ASTCENC_SIMD_INLINE vint16() = default;

	ASTCENC_SIMD_INLINE explicit vint16(const int *p)
	{
		m = _mm512_loadu_si512(p);
	}

	ASTCENC_SIMD_INLINE explicit vint16(const uint8_t *p)
	{
		__m128i bytes;
		std::memcpy(&bytes, p, sizeof(bytes));
		m = _mm512_cvtepu8_epi32(bytes);
	}

	ASTCENC_SIMD_INLINE explicit vint16(int a)
	{
		m = _mm512_set1_epi32(a);
	}

	ASTCENC_SIMD_INLINE explicit vint16(__m512i a)
	{
		m = a;
	}

	static ASTCENC_SIMD_INLINE vint16 zero()
	{
		return vint16(_mm512_setzero_si512());
	}

	static ASTCENC_SIMD_INLINE vint16 load1(const int* p)
	{
		return vint16(_mm512_set1_epi32(*p));
	}

	static ASTCENC_SIMD_INLINE vint16 load(const uint8_t* p)
	{
		return vint16(_mm512_loadu_si512(p));
	}

	static ASTCENC_SIMD_INLINE vint16 loada(const int* p)
	{
		return vint16(_mm512_load_si512(p));
	}

	static ASTCENC_SIMD_INLINE vint16 lane_id()
	{
		return vint16(_mm512_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7,
		                                8, 9, 10, 11, 12, 13, 14, 15));
	}

	__m512i m;
};

// ============================================================================
// vmask16 data type (AVX-512 k-mask)
// ============================================================================

struct vmask16
{
	ASTCENC_SIMD_INLINE explicit vmask16(__mmask16 a)
	{
		m = a;
	}

	ASTCENC_SIMD_INLINE explicit vmask16(bool a)
	{
		m = a ? static_cast<__mmask16>(0xFFFF) : static_cast<__mmask16>(0);
	}

	__mmask16 m;
};

// ============================================================================
// vmask16 operators and functions
// ============================================================================

ASTCENC_SIMD_INLINE vmask16 operator|(vmask16 a, vmask16 b)
{
	return vmask16(static_cast<__mmask16>(a.m | b.m));
}

ASTCENC_SIMD_INLINE vmask16 operator&(vmask16 a, vmask16 b)
{
	return vmask16(static_cast<__mmask16>(a.m & b.m));
}

ASTCENC_SIMD_INLINE vmask16 operator^(vmask16 a, vmask16 b)
{
	return vmask16(static_cast<__mmask16>(a.m ^ b.m));
}

ASTCENC_SIMD_INLINE vmask16 operator~(vmask16 a)
{
	return vmask16(static_cast<__mmask16>(a.m ^ 0xFFFF));
}

ASTCENC_SIMD_INLINE unsigned int mask(vmask16 a)
{
	return static_cast<unsigned int>(a.m);
}

ASTCENC_SIMD_INLINE bool any(vmask16 a)
{
	return a.m != 0;
}

ASTCENC_SIMD_INLINE bool all(vmask16 a)
{
	return a.m == 0xFFFF;
}

// ============================================================================
// vint16 operators and functions
// ============================================================================

ASTCENC_SIMD_INLINE vint16 operator+(vint16 a, vint16 b)
{
	return vint16(_mm512_add_epi32(a.m, b.m));
}

ASTCENC_SIMD_INLINE vint16& operator+=(vint16& a, const vint16& b)
{
	a = a + b;
	return a;
}

ASTCENC_SIMD_INLINE vint16 operator-(vint16 a, vint16 b)
{
	return vint16(_mm512_sub_epi32(a.m, b.m));
}

ASTCENC_SIMD_INLINE vint16 operator*(vint16 a, vint16 b)
{
	return vint16(_mm512_mullo_epi32(a.m, b.m));
}

ASTCENC_SIMD_INLINE vint16 operator~(vint16 a)
{
	return vint16(_mm512_xor_si512(a.m, _mm512_set1_epi32(-1)));
}

ASTCENC_SIMD_INLINE vint16 operator|(vint16 a, vint16 b)
{
	return vint16(_mm512_or_si512(a.m, b.m));
}

ASTCENC_SIMD_INLINE vint16 operator&(vint16 a, vint16 b)
{
	return vint16(_mm512_and_si512(a.m, b.m));
}

ASTCENC_SIMD_INLINE vint16 operator^(vint16 a, vint16 b)
{
	return vint16(_mm512_xor_si512(a.m, b.m));
}

ASTCENC_SIMD_INLINE vmask16 operator==(vint16 a, vint16 b)
{
	return vmask16(_mm512_cmpeq_epi32_mask(a.m, b.m));
}

ASTCENC_SIMD_INLINE vmask16 operator!=(vint16 a, vint16 b)
{
	return vmask16(_mm512_cmpneq_epi32_mask(a.m, b.m));
}

ASTCENC_SIMD_INLINE vmask16 operator<(vint16 a, vint16 b)
{
	return vmask16(_mm512_cmplt_epi32_mask(a.m, b.m));
}

ASTCENC_SIMD_INLINE vmask16 operator>(vint16 a, vint16 b)
{
	return vmask16(_mm512_cmpgt_epi32_mask(a.m, b.m));
}

template <int s> ASTCENC_SIMD_INLINE vint16 lsl(vint16 a)
{
	return vint16(_mm512_slli_epi32(a.m, s));
}

template <int s> ASTCENC_SIMD_INLINE vint16 asr(vint16 a)
{
	return vint16(_mm512_srai_epi32(a.m, s));
}

template <int s> ASTCENC_SIMD_INLINE vint16 lsr(vint16 a)
{
	return vint16(_mm512_srli_epi32(a.m, s));
}

ASTCENC_SIMD_INLINE vint16 min(vint16 a, vint16 b)
{
	return vint16(_mm512_min_epi32(a.m, b.m));
}

ASTCENC_SIMD_INLINE vint16 max(vint16 a, vint16 b)
{
	return vint16(_mm512_max_epi32(a.m, b.m));
}

ASTCENC_SIMD_INLINE int hmin_s(vint16 a)
{
	return _mm512_reduce_min_epi32(a.m);
}

ASTCENC_SIMD_INLINE vint16 hmin(vint16 a)
{
	return vint16(hmin_s(a));
}

ASTCENC_SIMD_INLINE int hmax_s(vint16 a)
{
	return _mm512_reduce_max_epi32(a.m);
}

ASTCENC_SIMD_INLINE vint16 hmax(vint16 a)
{
	return vint16(hmax_s(a));
}

ASTCENC_SIMD_INLINE vint16 vint16_from_size(size_t a)
{
	assert(a <= std::numeric_limits<int>::max());
	return vint16(static_cast<int>(a));
}

ASTCENC_SIMD_INLINE void storea(vint16 a, int* p)
{
	_mm512_store_si512(p, a.m);
}

ASTCENC_SIMD_INLINE void store(vint16 a, int* p)
{
	_mm512_storeu_si512(p, a.m);
}

ASTCENC_SIMD_INLINE void store_nbytes(vint16 a, uint8_t* p)
{
	_mm_storeu_si128(reinterpret_cast<__m128i*>(p),
	                 _mm512_extracti32x4_epi32(a.m, 0));
}

ASTCENC_SIMD_INLINE void pack_and_store_low_bytes(vint16 v, uint8_t* p)
{
	_mm_storeu_si128(reinterpret_cast<__m128i*>(p), _mm512_cvtepi32_epi8(v.m));
}

ASTCENC_SIMD_INLINE vint16 select(vint16 a, vint16 b, vmask16 cond)
{
	return vint16(_mm512_mask_blend_epi32(cond.m, a.m, b.m));
}

// ============================================================================
// vfloat16 operators and functions
// ============================================================================

ASTCENC_SIMD_INLINE vfloat16 operator+(vfloat16 a, vfloat16 b)
{
	return vfloat16(_mm512_add_ps(a.m, b.m));
}

ASTCENC_SIMD_INLINE vfloat16& operator+=(vfloat16& a, const vfloat16& b)
{
	a = a + b;
	return a;
}

ASTCENC_SIMD_INLINE vfloat16 operator-(vfloat16 a, vfloat16 b)
{
	return vfloat16(_mm512_sub_ps(a.m, b.m));
}

ASTCENC_SIMD_INLINE vfloat16 operator*(vfloat16 a, vfloat16 b)
{
	return vfloat16(_mm512_mul_ps(a.m, b.m));
}

ASTCENC_SIMD_INLINE vfloat16 operator*(vfloat16 a, float b)
{
	return vfloat16(_mm512_mul_ps(a.m, _mm512_set1_ps(b)));
}

ASTCENC_SIMD_INLINE vfloat16 operator*(float a, vfloat16 b)
{
	return vfloat16(_mm512_mul_ps(_mm512_set1_ps(a), b.m));
}

ASTCENC_SIMD_INLINE vfloat16 operator/(vfloat16 a, vfloat16 b)
{
	return vfloat16(_mm512_div_ps(a.m, b.m));
}

ASTCENC_SIMD_INLINE vfloat16 operator/(vfloat16 a, float b)
{
	return vfloat16(_mm512_div_ps(a.m, _mm512_set1_ps(b)));
}

ASTCENC_SIMD_INLINE vfloat16 operator/(float a, vfloat16 b)
{
	return vfloat16(_mm512_div_ps(_mm512_set1_ps(a), b.m));
}

ASTCENC_SIMD_INLINE vmask16 operator==(vfloat16 a, vfloat16 b)
{
	return vmask16(_mm512_cmp_ps_mask(a.m, b.m, _CMP_EQ_OQ));
}

ASTCENC_SIMD_INLINE vmask16 operator!=(vfloat16 a, vfloat16 b)
{
	return vmask16(_mm512_cmp_ps_mask(a.m, b.m, _CMP_NEQ_UQ));
}

ASTCENC_SIMD_INLINE vmask16 operator<(vfloat16 a, vfloat16 b)
{
	return vmask16(_mm512_cmp_ps_mask(a.m, b.m, _CMP_LT_OQ));
}

ASTCENC_SIMD_INLINE vmask16 operator>(vfloat16 a, vfloat16 b)
{
	return vmask16(_mm512_cmp_ps_mask(a.m, b.m, _CMP_GT_OQ));
}

ASTCENC_SIMD_INLINE vmask16 operator<=(vfloat16 a, vfloat16 b)
{
	return vmask16(_mm512_cmp_ps_mask(a.m, b.m, _CMP_LE_OQ));
}

ASTCENC_SIMD_INLINE vmask16 operator>=(vfloat16 a, vfloat16 b)
{
	return vmask16(_mm512_cmp_ps_mask(a.m, b.m, _CMP_GE_OQ));
}

ASTCENC_SIMD_INLINE vfloat16 min(vfloat16 a, vfloat16 b)
{
	return vfloat16(_mm512_min_ps(a.m, b.m));
}

ASTCENC_SIMD_INLINE vfloat16 min(vfloat16 a, float b)
{
	return min(a, vfloat16(b));
}

ASTCENC_SIMD_INLINE vfloat16 max(vfloat16 a, vfloat16 b)
{
	return vfloat16(_mm512_max_ps(a.m, b.m));
}

ASTCENC_SIMD_INLINE vfloat16 max(vfloat16 a, float b)
{
	return max(a, vfloat16(b));
}

ASTCENC_SIMD_INLINE vfloat16 clamp(float minv, float maxv, vfloat16 a)
{
	a.m = _mm512_max_ps(a.m, _mm512_set1_ps(minv));
	a.m = _mm512_min_ps(a.m, _mm512_set1_ps(maxv));
	return a;
}

ASTCENC_SIMD_INLINE vfloat16 clampzo(vfloat16 a)
{
	a.m = _mm512_max_ps(a.m, _mm512_setzero_ps());
	a.m = _mm512_min_ps(a.m, _mm512_set1_ps(1.0f));
	return a;
}

ASTCENC_SIMD_INLINE vfloat16 abs(vfloat16 a)
{
	return vfloat16(_mm512_and_ps(a.m, _mm512_castsi512_ps(_mm512_set1_epi32(0x7fffffff))));
}

ASTCENC_SIMD_INLINE vfloat16 round(vfloat16 a)
{
	return vfloat16(_mm512_roundscale_ps(a.m, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC));
}

ASTCENC_SIMD_INLINE float hmin_s(vfloat16 a)
{
	return _mm512_reduce_min_ps(a.m);
}

ASTCENC_SIMD_INLINE vfloat16 hmin(vfloat16 a)
{
	return vfloat16(hmin_s(a));
}

ASTCENC_SIMD_INLINE float hmax_s(vfloat16 a)
{
	return _mm512_reduce_max_ps(a.m);
}

ASTCENC_SIMD_INLINE vfloat16 hmax(vfloat16 a)
{
	return vfloat16(hmax_s(a));
}

ASTCENC_SIMD_INLINE float hadd_s(vfloat16 a)
{
	// Four sequential 4-wide adds: invariant with 4-wide accumulation order
	// of pairing 128-bit lanes, matching AVX2's two 4-wide adds extended.
	vfloat4 q0(_mm512_extractf32x4_ps(a.m, 0));
	vfloat4 q1(_mm512_extractf32x4_ps(a.m, 1));
	vfloat4 q2(_mm512_extractf32x4_ps(a.m, 2));
	vfloat4 q3(_mm512_extractf32x4_ps(a.m, 3));
	return (hadd_s(q0) + hadd_s(q1)) + (hadd_s(q2) + hadd_s(q3));
}

ASTCENC_SIMD_INLINE vfloat16 select(vfloat16 a, vfloat16 b, vmask16 cond)
{
	return vfloat16(_mm512_mask_blend_ps(cond.m, a.m, b.m));
}

ASTCENC_SIMD_INLINE void haccumulate(vfloat4& accum, vfloat16 a)
{
	haccumulate(accum, vfloat4(_mm512_extractf32x4_ps(a.m, 0)));
	haccumulate(accum, vfloat4(_mm512_extractf32x4_ps(a.m, 1)));
	haccumulate(accum, vfloat4(_mm512_extractf32x4_ps(a.m, 2)));
	haccumulate(accum, vfloat4(_mm512_extractf32x4_ps(a.m, 3)));
}

ASTCENC_SIMD_INLINE void haccumulate(vfloat16& accum, vfloat16 a)
{
	accum += a;
}

ASTCENC_SIMD_INLINE void haccumulate(vfloat4& accum, vfloat16 a, vmask16 m)
{
	a = select(vfloat16::zero(), a, m);
	haccumulate(accum, a);
}

ASTCENC_SIMD_INLINE void haccumulate(vfloat16& accum, vfloat16 a, vmask16 m)
{
	a = select(vfloat16::zero(), a, m);
	haccumulate(accum, a);
}

ASTCENC_SIMD_INLINE vfloat16 sqrt(vfloat16 a)
{
	return vfloat16(_mm512_sqrt_ps(a.m));
}

/**
 * @brief Expand 16 dword indices into VBMI byte indices for a float table.
 *
 * Lane i with float index k becomes bytes {4k, 4k+1, 4k+2, 4k+3}.
 */
ASTCENC_SIMD_INLINE __m512i vbmi_float_byte_idx(__m512i idx32)
{
	return _mm512_add_epi32(
		_mm512_mullo_epi32(idx32, _mm512_set1_epi32(0x04040404)),
		_mm512_set1_epi32(0x03020100));
}

ASTCENC_SIMD_INLINE __mmask16 vbmi_load_mask(unsigned int count, unsigned int offset)
{
	if (count <= offset)
	{
		return 0;
	}
	unsigned int n = count - offset;
	return n >= 16 ? static_cast<__mmask16>(0xFFFF)
	               : static_cast<__mmask16>((1u << n) - 1u);
}

ASTCENC_SIMD_INLINE vfloat16 gatherf_scalar(const float* base, vint16 indices)
{
	alignas(64) int idx[16];
	_mm512_store_si512(idx, indices.m);
	return vfloat16(_mm512_setr_ps(
		base[idx[0]],  base[idx[1]],  base[idx[2]],  base[idx[3]],
		base[idx[4]],  base[idx[5]],  base[idx[6]],  base[idx[7]],
		base[idx[8]],  base[idx[9]],  base[idx[10]], base[idx[11]],
		base[idx[12]], base[idx[13]], base[idx[14]], base[idx[15]]));
}

/**
 * @brief Register-resident float table for VBMI gatherf (up to 64 entries).
 *
 * Load once with vgatherf_load; permute many times with gatherf. Tables larger
 * than 64 keep the pointer and use the scalar-load gather.
 */
struct vgatherf_table {
	__m512 t0;
	__m512 t1;
	__m512 t2;
	__m512 t3;
	const float* base;
	unsigned int count;
};

#define ASTCENC_HAS_VGATHERF_TABLE 1

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

/**
 * @brief Gather 16 floats from a preloaded VBMI table, else scalar loads.
 *
 * Replaces AVX-512F vgatherdps: convert dword indices to byte indices and
 * vpermb (16-entry table) or vpermt2b (32-entry, or two vpermt2b for 64).
 */
ASTCENC_SIMD_INLINE vfloat16 gatherf(const vgatherf_table& t, vint16 indices)
{
	if (t.count == 0 || t.count > 64)
	{
		return gatherf_scalar(t.base, indices);
	}

	if (t.count <= 16)
	{
		__m512i r = _mm512_permutexvar_epi8(
			vbmi_float_byte_idx(indices.m), _mm512_castps_si512(t.t0));
		return vfloat16(_mm512_castsi512_ps(r));
	}

	if (t.count <= 32)
	{
		__m512i r = _mm512_permutex2var_epi8(
			_mm512_castps_si512(t.t0),
			vbmi_float_byte_idx(indices.m),
			_mm512_castps_si512(t.t1));
		return vfloat16(_mm512_castsi512_ps(r));
	}

	__m512i bidx_lo = vbmi_float_byte_idx(indices.m);
	__m512i bidx_hi = vbmi_float_byte_idx(
		_mm512_sub_epi32(indices.m, _mm512_set1_epi32(32)));
	__m512i lo = _mm512_permutex2var_epi8(
		_mm512_castps_si512(t.t0), bidx_lo, _mm512_castps_si512(t.t1));
	__m512i hi = _mm512_permutex2var_epi8(
		_mm512_castps_si512(t.t2), bidx_hi, _mm512_castps_si512(t.t3));
	__mmask16 ge32 = _mm512_cmpge_epi32_mask(indices.m, _mm512_set1_epi32(32));
	return vfloat16(_mm512_mask_blend_ps(ge32, _mm512_castsi512_ps(lo),
	                                      _mm512_castsi512_ps(hi)));
}

ASTCENC_SIMD_INLINE vfloat16 gatherf(const vgatherf_table& t, const uint8_t* indices)
{
	return gatherf(t, vint16(indices));
}

/**
 * @brief Gather 16 floats with AVX-512 VBMI permutes, else scalar loads.
 *
 * Convenience wrapper: infers the live table span from the index vector
 * (hmax + 1) then permutes. Prefer vgatherf_load when the same table is
 * gathered many times. Masked table loads use only [0, max_index] so
 * mid-array callers stay in bounds.
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

template<>
ASTCENC_SIMD_INLINE vfloat16 gatherf_byte_inds<vfloat16>(const float* base, const uint8_t* indices)
{
	return gatherf(base, vint16(indices));
}

ASTCENC_SIMD_INLINE void store(vfloat16 a, float* p)
{
	_mm512_storeu_ps(p, a.m);
}

ASTCENC_SIMD_INLINE void storea(vfloat16 a, float* p)
{
	_mm512_store_ps(p, a.m);
}

ASTCENC_SIMD_INLINE vint16 float_to_int(vfloat16 a)
{
	return vint16(_mm512_cvttps_epi32(a.m));
}

ASTCENC_SIMD_INLINE vint16 float_to_int_rtn(vfloat16 a)
{
	a = a + vfloat16(0.5f);
	return vint16(_mm512_cvttps_epi32(a.m));
}

ASTCENC_SIMD_INLINE vfloat16 int_to_float(vint16 a)
{
	return vfloat16(_mm512_cvtepi32_ps(a.m));
}

ASTCENC_SIMD_INLINE vint16 float_as_int(vfloat16 a)
{
	return vint16(_mm512_castps_si512(a.m));
}

ASTCENC_SIMD_INLINE vfloat16 int_as_float(vint16 a)
{
	return vfloat16(_mm512_castsi512_ps(a.m));
}

/*
 * Packed 8-bit tables in a zmm. VBMI vpermb indexes any of 64 bytes in one
 * instruction, so 16/32/64-entry tables are a single register (no AVX2-style
 * 16-byte broadcast or XOR-chain pshufb).
 */
struct vtable16_16x8 {
	vint16 t0;
};

struct vtable16_32x8 {
	vint16 t0;
};

struct vtable16_64x8 {
	vint16 t0;
};

ASTCENC_SIMD_INLINE vint16 vtable_lookup_32bit_vbmi(vint16 table, vint16 idx)
{
	__m512i r = _mm512_permutexvar_epi8(idx.m, table.m);
	return vint16(_mm512_and_si512(r, _mm512_set1_epi32(0xFF)));
}

ASTCENC_SIMD_INLINE void vtable_prepare(
	vtable16_16x8& table,
	const uint8_t* data
) {
	table.t0 = vint16(_mm512_zextsi128_si512(_mm_loadu_si128(
		reinterpret_cast<const __m128i*>(data))));
}

ASTCENC_SIMD_INLINE void vtable_prepare(
	vtable16_32x8& table,
	const uint8_t* data
) {
	table.t0 = vint16(_mm512_zextsi256_si512(_mm256_loadu_si256(
		reinterpret_cast<const __m256i*>(data))));
}

ASTCENC_SIMD_INLINE void vtable_prepare(
	vtable16_64x8& table,
	const uint8_t* data
) {
	table.t0 = vint16(_mm512_loadu_si512(data));
}

ASTCENC_SIMD_INLINE vint16 vtable_lookup_32bit(
	const vtable16_16x8& tbl,
	vint16 idx
) {
	return vtable_lookup_32bit_vbmi(tbl.t0, idx);
}

ASTCENC_SIMD_INLINE vint16 vtable_lookup_32bit(
	const vtable16_32x8& tbl,
	vint16 idx
) {
	return vtable_lookup_32bit_vbmi(tbl.t0, idx);
}

ASTCENC_SIMD_INLINE vint16 vtable_lookup_32bit(
	const vtable16_64x8& tbl,
	vint16 idx
) {
	return vtable_lookup_32bit_vbmi(tbl.t0, idx);
}

ASTCENC_SIMD_INLINE vint16 interleave_rgba8(vint16 r, vint16 g, vint16 b, vint16 a)
{
	return r + lsl<8>(g) + lsl<16>(b) + lsl<24>(a);
}

ASTCENC_SIMD_INLINE void store_lanes_masked(uint8_t* base, vint16 data, vmask16 mask)
{
	_mm512_mask_storeu_epi32(base, mask.m, data.m);
}

ASTCENC_SIMD_INLINE void print(vint16 a)
{
	alignas(64) int v[16];
	storea(a, v);
	printf("v16_i32:\n  %8d %8d %8d %8d %8d %8d %8d %8d\n  %8d %8d %8d %8d %8d %8d %8d %8d\n",
	       v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7],
	       v[8], v[9], v[10], v[11], v[12], v[13], v[14], v[15]);
}

ASTCENC_SIMD_INLINE void printx(vint16 a)
{
	alignas(64) int v[16];
	storea(a, v);
	unsigned int uv[16];
	std::memcpy(uv, v, sizeof(uv));
	printf("v16_i32:\n  %08x %08x %08x %08x %08x %08x %08x %08x\n  %08x %08x %08x %08x %08x %08x %08x %08x\n",
	       uv[0], uv[1], uv[2], uv[3], uv[4], uv[5], uv[6], uv[7],
	       uv[8], uv[9], uv[10], uv[11], uv[12], uv[13], uv[14], uv[15]);
}

ASTCENC_SIMD_INLINE void print(vfloat16 a)
{
	alignas(64) float v[16];
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

ASTCENC_SIMD_INLINE void print(vmask16 a)
{
	print(select(vint16(0), vint16(1), a));
}

#endif // ASTC_VECMATHLIB_AVX512_16_H_INCLUDED
