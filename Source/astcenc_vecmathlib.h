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

/*
 * This module implements N-wide float and integer vectors and common
 * operations on them; with N depending on underlying ISA (e.g. 4 for SSE, 1 for
 * pure scalar).
 */

#ifndef ASTC_VECMATHLIB_H_INCLUDED
#define ASTC_VECMATHLIB_H_INCLUDED

#if ASTCENC_SSE != 0 || ASTCENC_AVX != 0
	#include <immintrin.h>
#endif

#if defined(_MSC_VER)
	#define SIMD_INLINE __forceinline
#else
	#define SIMD_INLINE __attribute__((unused, always_inline, nodebug)) inline
#endif

#if ASTCENC_AVX >= 2
	#define ASTCENC_SIMD_ISA_AVX2
#elif ASTCENC_SSE >= 20
	#define ASTCENC_SIMD_ISA_SSE
#else
	#define ASTCENC_SIMD_ISA_SCALAR
#endif


// ----------------------------------------------------------------------------
// AVX2 8-wide implementation

#ifdef ASTCENC_SIMD_ISA_AVX2

#define ASTCENC_SIMD_WIDTH 8

// N-wide float
struct vfloat
{
	SIMD_INLINE vfloat() {}
	// Initialize with N floats from an unaligned memory address.
	// Using loada() when address is aligned might be more optimal.
	SIMD_INLINE explicit vfloat(const float *p) { m = _mm256_loadu_ps(p); }
	// Initialize with the same given float value in all lanes.
	SIMD_INLINE explicit vfloat(float v) { m = _mm256_set1_ps(v); }

	SIMD_INLINE explicit vfloat(__m256 v) { m = v; }

	// Get SIMD lane #i value.
	SIMD_INLINE float lane(int i) const
	{
		#ifdef _MSC_VER
		return m.m256_f32[i];
		#else
		union { __m256 m; float f[ASTCENC_SIMD_WIDTH]; } cvt;
		cvt.m = m;
		return cvt.f[i];
		#endif
	}

	// Float vector with all zero values
	static SIMD_INLINE vfloat zero() { return vfloat(_mm256_setzero_ps()); }

	// Float vector with each lane having the lane index (0, 1, 2, ...)
	static SIMD_INLINE vfloat lane_id() { return vfloat(_mm256_set_ps(7, 6, 5, 4, 3, 2, 1, 0)); }

	__m256 m;
};

// N-wide integer (32 bit in each lane)
struct vint
{
	SIMD_INLINE vint() {}
	// Initialize with N ints from an unaligned memory address.
	SIMD_INLINE explicit vint(const int *p) { m = _mm256_loadu_si256((const __m256i*)p); }
	// Initialize with the same given integer value in all lanes.
	SIMD_INLINE explicit vint(int v) { m = _mm256_set1_epi32(v); }

	SIMD_INLINE explicit vint(__m256i v) { m = v; }

	// Get SIMD lane #i value
	SIMD_INLINE int lane(int i) const
	{
		#ifdef _MSC_VER
		return m.m256i_i32[i];
		#else
		union { __m256i m; int f[ASTCENC_SIMD_WIDTH]; } cvt;
		cvt.m = m;
		return cvt.f[i];
		#endif
	}

	// Integer vector with each lane having the lane index (0, 1, 2, ...)
	static SIMD_INLINE vint lane_id() { return vint(_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0)); }

	__m256i m;
};

// N-wide comparison mask. vmask is a result of comparison operators,
// and an argument for select() function below.
struct vmask
{
	SIMD_INLINE explicit vmask(__m256 v) { m = v; }
	SIMD_INLINE explicit vmask(__m256i v) { m = _mm256_castsi256_ps(v); }
	__m256 m;
};

// Initialize with one float in all SIMD lanes, from an aligned memory address.
SIMD_INLINE vfloat load1a(const float* p) { return vfloat(_mm256_broadcast_ss(p)); }
// Initialize with N floats from an aligned memory address.
SIMD_INLINE vfloat loada(const float* p) { return vfloat(_mm256_load_ps(p)); }

// Per-lane float arithmetic operations
SIMD_INLINE vfloat operator+ (vfloat a, vfloat b) { a.m = _mm256_add_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat operator- (vfloat a, vfloat b) { a.m = _mm256_sub_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat operator* (vfloat a, vfloat b) { a.m = _mm256_mul_ps(a.m, b.m); return a; }

// Per-lane float comparison operations
SIMD_INLINE vmask operator==(vfloat a, vfloat b) { return vmask(_mm256_cmp_ps(a.m, b.m, _CMP_EQ_OQ)); }
SIMD_INLINE vmask operator!=(vfloat a, vfloat b) { return vmask(_mm256_cmp_ps(a.m, b.m, _CMP_NEQ_OQ)); }
SIMD_INLINE vmask operator< (vfloat a, vfloat b) { return vmask(_mm256_cmp_ps(a.m, b.m, _CMP_LT_OQ)); }
SIMD_INLINE vmask operator> (vfloat a, vfloat b) { return vmask(_mm256_cmp_ps(a.m, b.m, _CMP_GT_OQ)); }
SIMD_INLINE vmask operator<=(vfloat a, vfloat b) { return vmask(_mm256_cmp_ps(a.m, b.m, _CMP_LE_OQ)); }
SIMD_INLINE vmask operator>=(vfloat a, vfloat b) { return vmask(_mm256_cmp_ps(a.m, b.m, _CMP_GE_OQ)); }

// Logical operations on comparison mask values
SIMD_INLINE vmask operator| (vmask a, vmask b) { return vmask(_mm256_or_ps(a.m, b.m)); }
SIMD_INLINE vmask operator& (vmask a, vmask b) { return vmask(_mm256_and_ps(a.m, b.m)); }
SIMD_INLINE vmask operator^ (vmask a, vmask b) { return vmask(_mm256_xor_ps(a.m, b.m)); }

// Returns a 8-bit code where bit0..bit7 map to lanes
SIMD_INLINE unsigned mask(vmask v) { return _mm256_movemask_ps(v.m); }
// Whether any lane in the comparison mask is set
SIMD_INLINE bool any(vmask v) { return mask(v) != 0; }
// Whether all lanes in the comparison mask are set
SIMD_INLINE bool all(vmask v) { return mask(v) == 0xFF; }

// Per-lane float min & max
SIMD_INLINE vfloat min(vfloat a, vfloat b) { a.m = _mm256_min_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat max(vfloat a, vfloat b) { a.m = _mm256_max_ps(a.m, b.m); return a; }

// Per-lane clamp to 0..1 range
SIMD_INLINE vfloat saturate(vfloat a)
{
	__m256 zero = _mm256_setzero_ps();
	__m256 one = _mm256_set1_ps(1.0f);
	return vfloat(_mm256_min_ps(_mm256_max_ps(a.m, zero), one));
}

// Round to nearest integer (nearest even for .5 cases)
SIMD_INLINE vfloat round(vfloat v)
{
	return vfloat(_mm256_round_ps(v.m, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC));
}

// Per-lane convert to integer (truncate)
SIMD_INLINE vint floatToInt(vfloat v) { return vint(_mm256_cvttps_epi32(v.m)); }

// Reinterpret-bitcast integer vector as a float vector (this is basically a no-op on the CPU)
SIMD_INLINE vfloat intAsFloat(vint v) { return vfloat(_mm256_castsi256_ps(v.m)); }

SIMD_INLINE vint operator~ (vint a) { return vint(_mm256_xor_si256(a.m, _mm256_set1_epi32(-1))); }
SIMD_INLINE vmask operator~ (vmask a) { return vmask(_mm256_xor_si256(_mm256_castps_si256(a.m), _mm256_set1_epi32(-1))); }

// Per-lane arithmetic integer operations
SIMD_INLINE vint operator+ (vint a, vint b) { a.m = _mm256_add_epi32(a.m, b.m); return a; }
SIMD_INLINE vint operator- (vint a, vint b) { a.m = _mm256_sub_epi32(a.m, b.m); return a; }

// Per-lane integer comparison operations
SIMD_INLINE vmask operator< (vint a, vint b) { return vmask(_mm256_cmpgt_epi32(b.m, a.m)); }
SIMD_INLINE vmask operator> (vint a, vint b) { return vmask(_mm256_cmpgt_epi32(a.m, b.m)); }
SIMD_INLINE vmask operator==(vint a, vint b) { return vmask(_mm256_cmpeq_epi32(a.m, b.m)); }
SIMD_INLINE vmask operator!=(vint a, vint b) { return ~vmask(_mm256_cmpeq_epi32(a.m, b.m)); }

// Per-lane integer min & max
SIMD_INLINE vint min(vint a, vint b) { a.m = _mm256_min_epi32(a.m, b.m); return a; }
SIMD_INLINE vint max(vint a, vint b) { a.m = _mm256_max_epi32(a.m, b.m); return a; }

// Horizontal minimum - returns vector with all lanes
// set to the minimum value of the input vector.
SIMD_INLINE vfloat hmin(vfloat v)
{
	__m128 vlow  = _mm256_castps256_ps128(v.m);
	__m128 vhigh = _mm256_extractf128_ps(v.m, 1);
		   vlow  = _mm_min_ps(vlow, vhigh);

	// First do an horizontal reduction.                                // v = [ D C | B A ]
	__m128 shuf = _mm_shuffle_ps(vlow, vlow, _MM_SHUFFLE(2, 3, 0, 1));  //     [ C D | A B ]
	__m128 mins = _mm_min_ps(vlow, shuf);                            // mins = [ D+C C+D | B+A A+B ]
	shuf        = _mm_movehl_ps(shuf, mins);                         //        [   C   D | D+C C+D ]
	mins        = _mm_min_ss(mins, shuf);

	vfloat vmin(_mm256_permute_ps(_mm256_set_m128(mins, mins), 0)); // _MM256_PERMUTE(0, 0, 0, 0, 0, 0, 0, 0)
	return vmin;
}

SIMD_INLINE vint hmin(vint v)
{
	__m128i m = _mm_min_epi32(_mm256_extracti128_si256(v.m, 0), _mm256_extracti128_si256(v.m, 1));
	m = _mm_min_epi32(m, _mm_shuffle_epi32(m, _MM_SHUFFLE(0,0,3,2)));
	m = _mm_min_epi32(m, _mm_shuffle_epi32(m, _MM_SHUFFLE(0,0,0,1)));
	m = _mm_shuffle_epi32(m, _MM_SHUFFLE(0,0,0,0));
	vint vmin(_mm256_set_m128i(m, m));
	return vmin;
}

// Store float vector into an aligned address.
SIMD_INLINE void store(vfloat v, float* ptr) { _mm256_store_ps(ptr, v.m); }
// Store integer vector into an aligned address.
SIMD_INLINE void store(vint v, int* ptr) { _mm256_store_si256((__m256i*)ptr, v.m); }

// Store lowest N (simd width) bytes of integer vector into an unaligned address.
SIMD_INLINE void store_nbytes(vint v, uint8_t* ptr) { _mm_storeu_si64(ptr, _mm256_extracti128_si256(v.m, 0)); }

// SIMD "gather" - load each lane with base[indices[i]]
SIMD_INLINE vfloat gatherf(const float* base, vint indices)
{
	return vfloat(_mm256_i32gather_ps(base, indices.m, 4));
}
SIMD_INLINE vint gatheri(const int* base, vint indices)
{
	return vint(_mm256_i32gather_epi32(base, indices.m, 4));
}

// Pack low 8 bits of each lane into low 64 bits of result.
SIMD_INLINE vint pack_low_bytes(vint v)
{
	__m256i shuf = _mm256_set_epi8(0,0,0,0,0,0,0,0, 0,0,0,0,28,24,20,16, 0,0,0,0,0,0,0,0, 0,0,0,0,12,8,4,0);
	__m256i a = _mm256_shuffle_epi8(v.m, shuf);
	__m128i a0 = _mm256_extracti128_si256(a, 0);
	__m128i a1 = _mm256_extracti128_si256(a, 1);
	__m128i b = _mm_unpacklo_epi32(a0, a1);
	return vint(_mm256_set_m128i(b, b));
}

// "select", i.e. highbit(cond) ? b : a
SIMD_INLINE vfloat select(vfloat a, vfloat b, vmask cond)
{
	return vfloat(_mm256_blendv_ps(a.m, b.m, cond.m));
}
SIMD_INLINE vint select(vint a, vint b, vmask cond)
{
	return vint(_mm256_blendv_epi8(a.m, b.m, _mm256_castps_si256(cond.m)));
}

SIMD_INLINE void print(vfloat a)
{
	alignas(ASTCENC_VECALIGN) float v[8];
	store(a, v);
	printf("v8_f32:\n  %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f\n",
	       (double)v[0], (double)v[1], (double)v[2], (double)v[3],
	       (double)v[4], (double)v[5], (double)v[6], (double)v[7]);
}

SIMD_INLINE void print(vint a)
{
	alignas(ASTCENC_VECALIGN) int v[8];
	store(a, v);
	printf("v8_i32:\n  %8u %8u %8u %8u %8u %8u %8u %8u\n",
	       v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]);
}

#endif // #ifdef ASTCENC_SIMD_ISA_AVX2


// ----------------------------------------------------------------------------
// SSE 4-wide implementation
// Uses SSE2 as baseline, optionally SSE4.x instructions based on ASTCENC_SSE value

#ifdef ASTCENC_SIMD_ISA_SSE

#define ASTCENC_SIMD_WIDTH 4

struct vfloat
{
	SIMD_INLINE vfloat() {}
	SIMD_INLINE explicit vfloat(const float *p) { m = _mm_loadu_ps(p); }
	SIMD_INLINE explicit vfloat(float v) { m = _mm_set_ps1(v); }
	SIMD_INLINE explicit vfloat(__m128 v) { m = v; }
	SIMD_INLINE float lane(int i) const
	{
		#ifdef _MSC_VER
		return m.m128_f32[i];
		#else
		union { __m128 m; float f[ASTCENC_SIMD_WIDTH]; } cvt;
		cvt.m = m;
		return cvt.f[i];
		#endif
	}
	static SIMD_INLINE vfloat zero() { return vfloat(_mm_setzero_ps()); }
	static SIMD_INLINE vfloat lane_id() { return vfloat(_mm_set_ps(3, 2, 1, 0)); }
	__m128 m;
};

struct vint
{
	SIMD_INLINE vint() {}
	SIMD_INLINE explicit vint(const int *p) { m = _mm_load_si128((const __m128i*)p); }
	SIMD_INLINE explicit vint(int v) { m = _mm_set1_epi32(v); }
	SIMD_INLINE explicit vint(__m128i v) { m = v; }
	SIMD_INLINE int lane(int i) const
	{
		#ifdef _MSC_VER
		return m.m128i_i32[i];
		#else
		union { __m128i m; int f[ASTCENC_SIMD_WIDTH]; } cvt;
		cvt.m = m;
		return cvt.f[i];
		#endif
	}
	static SIMD_INLINE vint lane_id() { return vint(_mm_set_epi32(3, 2, 1, 0)); }
	__m128i m;
};

struct vmask
{
	SIMD_INLINE explicit vmask(__m128 v) { m = v; }
	SIMD_INLINE explicit vmask(__m128i v) { m = _mm_castsi128_ps(v); }
	__m128 m;
};


SIMD_INLINE vfloat load1a(const float* p) { return vfloat(_mm_load_ps1(p)); }
SIMD_INLINE vfloat loada(const float* p) { return vfloat(_mm_load_ps(p)); }

SIMD_INLINE vfloat operator+ (vfloat a, vfloat b) { a.m = _mm_add_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat operator- (vfloat a, vfloat b) { a.m = _mm_sub_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat operator* (vfloat a, vfloat b) { a.m = _mm_mul_ps(a.m, b.m); return a; }
SIMD_INLINE vmask operator==(vfloat a, vfloat b) { return vmask(_mm_cmpeq_ps(a.m, b.m)); }
SIMD_INLINE vmask operator!=(vfloat a, vfloat b) { return vmask(_mm_cmpneq_ps(a.m, b.m)); }
SIMD_INLINE vmask operator< (vfloat a, vfloat b) { return vmask(_mm_cmplt_ps(a.m, b.m)); }
SIMD_INLINE vmask operator> (vfloat a, vfloat b) { return vmask(_mm_cmpgt_ps(a.m, b.m)); }
SIMD_INLINE vmask operator<=(vfloat a, vfloat b) { return vmask(_mm_cmple_ps(a.m, b.m)); }
SIMD_INLINE vmask operator>=(vfloat a, vfloat b) { return vmask(_mm_cmpge_ps(a.m, b.m)); }
SIMD_INLINE vmask operator| (vmask a, vmask b) { return vmask(_mm_or_ps(a.m, b.m)); }
SIMD_INLINE vmask operator& (vmask a, vmask b) { return vmask(_mm_and_ps(a.m, b.m)); }
SIMD_INLINE vmask operator^ (vmask a, vmask b) { return vmask(_mm_xor_ps(a.m, b.m)); }
// Returns a 4-bit code where bit0..bit3 is X..W
SIMD_INLINE unsigned mask(vmask v) { return _mm_movemask_ps(v.m); }
SIMD_INLINE bool any(vmask v) { return mask(v) != 0; }
SIMD_INLINE bool all(vmask v) { return mask(v) == 0xF; }

SIMD_INLINE vfloat min(vfloat a, vfloat b) { a.m = _mm_min_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat max(vfloat a, vfloat b) { a.m = _mm_max_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat saturate(vfloat a)
{
	__m128 zero = _mm_setzero_ps();
	__m128 one = _mm_set1_ps(1.0f);
	return vfloat(_mm_min_ps(_mm_max_ps(a.m, zero), one));
}

SIMD_INLINE vfloat round(vfloat v)
{
#if ASTCENC_SSE >= 41
	return vfloat(_mm_round_ps(v.m, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC));
#else
	__m128 V = v.m;
	__m128 negZero = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
	__m128 noFraction = _mm_set_ps1(8388608.0f);
	__m128 absMask = _mm_castsi128_ps(_mm_set1_epi32(0x7FFFFFFF));
	__m128 sign = _mm_and_ps(V, negZero);
	__m128 sMagic = _mm_or_ps(noFraction, sign);
	__m128 R1 = _mm_add_ps(V, sMagic);
	R1 = _mm_sub_ps(R1, sMagic);
	__m128 R2 = _mm_and_ps(V, absMask);
	__m128 mask = _mm_cmple_ps(R2, noFraction);
	R2 = _mm_andnot_ps(mask, V);
	R1 = _mm_and_ps(R1, mask);
	return vfloat(_mm_xor_ps(R1, R2));
#endif
}

SIMD_INLINE vint floatToInt(vfloat v) { return vint(_mm_cvttps_epi32(v.m)); }

SIMD_INLINE vfloat intAsFloat(vint v) { return vfloat(_mm_castsi128_ps(v.m)); }

SIMD_INLINE vint operator~ (vint a) { return vint(_mm_xor_si128(a.m, _mm_set1_epi32(-1))); }
SIMD_INLINE vmask operator~ (vmask a) { return vmask(_mm_xor_si128(_mm_castps_si128(a.m), _mm_set1_epi32(-1))); }

SIMD_INLINE vint operator+ (vint a, vint b) { a.m = _mm_add_epi32(a.m, b.m); return a; }
SIMD_INLINE vint operator- (vint a, vint b) { a.m = _mm_sub_epi32(a.m, b.m); return a; }
SIMD_INLINE vmask operator< (vint a, vint b) { return vmask(_mm_cmplt_epi32(a.m, b.m)); }
SIMD_INLINE vmask operator> (vint a, vint b) { return vmask(_mm_cmpgt_epi32(a.m, b.m)); }
SIMD_INLINE vmask operator==(vint a, vint b) { return vmask(_mm_cmpeq_epi32(a.m, b.m)); }
SIMD_INLINE vmask operator!=(vint a, vint b) { return ~vmask(_mm_cmpeq_epi32(a.m, b.m)); }
SIMD_INLINE vint min(vint a, vint b) {
#if ASTCENC_SSE >= 41
	a.m = _mm_min_epi32(a.m, b.m);
#else
	vmask d = a < b;
	a.m = _mm_or_si128(_mm_and_si128(_mm_castps_si128(d.m), a.m), _mm_andnot_si128(_mm_castps_si128(d.m), b.m));
#endif
	return a;
}

SIMD_INLINE vint max(vint a, vint b) {
#if ASTCENC_SSE >= 41
	a.m = _mm_max_epi32(a.m, b.m);
#else
	vmask d = a > b;
	a.m = _mm_or_si128(_mm_and_si128(_mm_castps_si128(d.m), a.m), _mm_andnot_si128(_mm_castps_si128(d.m), b.m));
#endif
	return a;
}

#define ASTCENC_SHUFFLE4F(V, X,Y,Z,W) vfloat(_mm_shuffle_ps((V).m, (V).m, _MM_SHUFFLE(W,Z,Y,X)))
#define ASTCENC_SHUFFLE4I(V, X,Y,Z,W) vint(_mm_shuffle_epi32((V).m, _MM_SHUFFLE(W,Z,Y,X)))

SIMD_INLINE vfloat hmin(vfloat v)
{
	v = min(v, ASTCENC_SHUFFLE4F(v, 2, 3, 0, 0));
	v = min(v, ASTCENC_SHUFFLE4F(v, 1, 0, 0, 0));
	return ASTCENC_SHUFFLE4F(v, 0,0,0,0);
}
SIMD_INLINE vint hmin(vint v)
{
	v = min(v, ASTCENC_SHUFFLE4I(v, 2, 3, 0, 0));
	v = min(v, ASTCENC_SHUFFLE4I(v, 1, 0, 0, 0));
	return ASTCENC_SHUFFLE4I(v, 0,0,0,0);
}

SIMD_INLINE void store(vfloat v, float* ptr) { _mm_store_ps(ptr, v.m); }
SIMD_INLINE void store(vint v, int* ptr) { _mm_store_si128((__m128i*)ptr, v.m); }

SIMD_INLINE void store_nbytes(vint v, uint8_t* ptr) { _mm_storeu_si32(ptr, v.m); }

SIMD_INLINE vfloat gatherf(const float* base, vint indices)
{
	int idx[4];
	store(indices, idx);
	return vfloat(_mm_set_ps(base[idx[3]], base[idx[2]], base[idx[1]], base[idx[0]]));
}

SIMD_INLINE vint gatheri(const int* base, vint indices)
{
	int idx[4];
	store(indices, idx);
	return vint(_mm_set_epi32(base[idx[3]], base[idx[2]], base[idx[1]], base[idx[0]]));
}

// packs low 8 bits of each lane into low 32 bits of result
SIMD_INLINE vint pack_low_bytes(vint v)
{
	#if ASTCENC_SSE >= 41
	__m128i shuf = _mm_set_epi8(0,0,0,0, 0,0,0,0, 0,0,0,0, 12,8,4,0);
	return vint(_mm_shuffle_epi8(v.m, shuf));
	#else
	__m128i va = _mm_unpacklo_epi8(v.m, _mm_shuffle_epi32(v.m, _MM_SHUFFLE(1,1,1,1)));
	__m128i vb = _mm_unpackhi_epi8(v.m, _mm_shuffle_epi32(v.m, _MM_SHUFFLE(3,3,3,3)));
	return vint(_mm_unpacklo_epi16(va, vb));
	#endif
}

// "select", i.e. highbit(cond) ? b : a
// on SSE4.1 and up this can be done easily via "blend" instruction;
// on older SSEs we have to do some hoops, see
// https://fgiesen.wordpress.com/2016/04/03/sse-mind-the-gap/
SIMD_INLINE vfloat select(vfloat a, vfloat b, vmask cond)
{
#if ASTCENC_SSE >= 41
	a.m = _mm_blendv_ps(a.m, b.m, cond.m);
#else
	__m128 d = _mm_castsi128_ps(_mm_srai_epi32(_mm_castps_si128(cond.m), 31));
	a.m = _mm_or_ps(_mm_and_ps(d, b.m), _mm_andnot_ps(d, a.m));
#endif
	return a;
}

SIMD_INLINE vint select(vint a, vint b, vmask cond)
{
#if ASTCENC_SSE >= 41
	return vint(_mm_blendv_epi8(a.m, b.m, _mm_castps_si128(cond.m)));
#else
	__m128i d = _mm_srai_epi32(_mm_castps_si128(cond.m), 31);
	return vint(_mm_or_si128(_mm_and_si128(d, b.m), _mm_andnot_si128(d, a.m)));
#endif
}

SIMD_INLINE void print(vfloat a)
{
	alignas(ASTCENC_VECALIGN) float v[4];
	store(a, v);
	printf("v4_f32:\n  %0.4f %0.4f %0.4f %0.4f\n",
	       (double)v[0], (double)v[1], (double)v[2], (double)v[3]);
}

SIMD_INLINE void print(vint a)
{
	alignas(ASTCENC_VECALIGN) int v[4];
	store(a, v);
	printf("v4_i32:\n  %8u %8u %8u %8u\n",
	       v[0], v[1], v[2], v[3]);
}


#endif // #ifdef ASTCENC_SIMD_ISA_SSE


// ----------------------------------------------------------------------------
// Pure scalar, 1-wide implementation

#ifdef ASTCENC_SIMD_ISA_SCALAR

#define ASTCENC_SIMD_WIDTH 1

struct vfloat
{
	SIMD_INLINE vfloat() {}
	SIMD_INLINE explicit vfloat(const float *p) { m = *p; }
	SIMD_INLINE explicit vfloat(float v) { m = v; }
	SIMD_INLINE float lane(int i) const { return m; }
	static SIMD_INLINE vfloat zero() { return vfloat(0.0f); }
	static SIMD_INLINE vfloat lane_id() { return vfloat(0.0f); }
	float m;
};

struct vint
{
	SIMD_INLINE vint() {}
	SIMD_INLINE explicit vint(const int *p) { m = *p; }
	SIMD_INLINE explicit vint(int v) { m = v; }
	SIMD_INLINE int lane(int i) const { return m; }
	static SIMD_INLINE vint lane_id() { return vint(0); }
	int m;
};

struct vmask
{
	SIMD_INLINE explicit vmask(bool v) { m = v; }
	bool m;
};


SIMD_INLINE vfloat load1a(const float* p) { return vfloat(*p); }
SIMD_INLINE vfloat loada(const float* p) { return vfloat(*p); }

SIMD_INLINE vfloat operator+ (vfloat a, vfloat b) { a.m = a.m + b.m; return a; }
SIMD_INLINE vfloat operator- (vfloat a, vfloat b) { a.m = a.m - b.m; return a; }
SIMD_INLINE vfloat operator* (vfloat a, vfloat b) { a.m = a.m * b.m; return a; }
SIMD_INLINE vmask operator==(vfloat a, vfloat b) { return vmask(a.m = a.m == b.m); }
SIMD_INLINE vmask operator!=(vfloat a, vfloat b) { return vmask(a.m = a.m != b.m); }
SIMD_INLINE vmask operator< (vfloat a, vfloat b) { return vmask(a.m = a.m < b.m); }
SIMD_INLINE vmask operator> (vfloat a, vfloat b) { return vmask(a.m = a.m > b.m); }
SIMD_INLINE vmask operator<=(vfloat a, vfloat b) { return vmask(a.m = a.m <= b.m); }
SIMD_INLINE vmask operator>=(vfloat a, vfloat b) { return vmask(a.m = a.m >= b.m); }
SIMD_INLINE vmask operator| (vmask a, vmask b) { return vmask(a.m || b.m); }
SIMD_INLINE vmask operator& (vmask a, vmask b) { return vmask(a.m && b.m); }
SIMD_INLINE vmask operator^ (vmask a, vmask b) { return vmask(a.m ^ b.m); }
SIMD_INLINE unsigned mask(vmask v) { return v.m; }
SIMD_INLINE bool any(vmask v) { return mask(v) != 0; }
SIMD_INLINE bool all(vmask v) { return mask(v) != 0; }

SIMD_INLINE vfloat min(vfloat a, vfloat b) { a.m = a.m < b.m ? a.m : b.m; return a; }
SIMD_INLINE vfloat max(vfloat a, vfloat b) { a.m = a.m > b.m ? a.m : b.m; return a; }
SIMD_INLINE vfloat saturate(vfloat a) { return vfloat(std::min(std::max(a.m,0.0f), 1.0f)); }

SIMD_INLINE vfloat round(vfloat v)
{
	return vfloat(std::floor(v.m + 0.5f));
}

SIMD_INLINE vint floatToInt(vfloat v) { return vint(v.m); }

SIMD_INLINE vfloat intAsFloat(vint v) { vfloat r; memcpy(&r.m, &v.m, 4); return r; }

SIMD_INLINE vint operator~ (vint a) { a.m = ~a.m; return a; }
SIMD_INLINE vint operator+ (vint a, vint b) { a.m = a.m + b.m; return a; }
SIMD_INLINE vint operator- (vint a, vint b) { a.m = a.m - b.m; return a; }
SIMD_INLINE vmask operator< (vint a, vint b) { return vmask(a.m = a.m < b.m); }
SIMD_INLINE vmask operator> (vint a, vint b) { return vmask(a.m = a.m > b.m); }
SIMD_INLINE vmask operator==(vint a, vint b) { return vmask(a.m = a.m == b.m); }
SIMD_INLINE vmask operator!=(vint a, vint b) { return vmask(a.m = a.m != b.m); }
SIMD_INLINE vint min(vint a, vint b) { a.m = a.m < b.m ? a.m : b.m; return a; }
SIMD_INLINE vint max(vint a, vint b) { a.m = a.m > b.m ? a.m : b.m; return a; }

SIMD_INLINE vfloat hmin(vfloat v) { return v; }
SIMD_INLINE vint hmin(vint v) { return v; }

SIMD_INLINE void store(vfloat v, float* ptr) { *ptr = v.m; }
SIMD_INLINE void store(vint v, int* ptr) { *ptr = v.m; }

SIMD_INLINE void store_nbytes(vint v, uint8_t* ptr) { *ptr = (uint8_t)v.m; }

SIMD_INLINE vfloat gatherf(const float* base, vint indices)
{
	return vfloat(base[indices.m]);
}
SIMD_INLINE vint gatheri(const int* base, vint indices)
{
	return vint(base[indices.m]);
}

// packs low 8 bits of each lane into low 8 bits of result (a no-op in scalar code path)
SIMD_INLINE vint pack_low_bytes(vint v)
{
	return v;
}


// "select", i.e. highbit(cond) ? b : a
SIMD_INLINE vfloat select(vfloat a, vfloat b, vmask cond)
{
	return cond.m ? b : a;
}
SIMD_INLINE vint select(vint a, vint b, vmask cond)
{
	return cond.m ? b : a;
}


#endif // #ifdef ASTCENC_SIMD_ISA_SCALAR

#endif // #ifndef ASTC_VECMATHLIB_H_INCLUDED
