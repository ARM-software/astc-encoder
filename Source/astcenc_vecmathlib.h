// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2011-2020 Arm Limited
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

struct vfloat
{
    SIMD_INLINE vfloat() {}
    SIMD_INLINE explicit vfloat(const float *p) { m = _mm256_loadu_ps(p); }
    SIMD_INLINE explicit vfloat(float v) { m = _mm256_set1_ps(v); }
    SIMD_INLINE explicit vfloat(__m256 v) { m = v; }
    static vfloat zero() { return vfloat(_mm256_setzero_ps()); }
    __m256 m;
};

struct vint
{
    SIMD_INLINE vint() {}
    SIMD_INLINE explicit vint(const int *p) { m = _mm256_loadu_si256((const __m256i*)p); }
    SIMD_INLINE explicit vint(int v) { m = _mm256_set1_epi32(v); }
    SIMD_INLINE explicit vint(__m256i v) { m = v; }
    __m256i m;
};

SIMD_INLINE vfloat load1a(const float* p) { return vfloat(_mm256_broadcast_ss(p)); }
SIMD_INLINE vfloat loada(const float* p) { return vfloat(_mm256_load_ps(p)); }

SIMD_INLINE vfloat operator+ (vfloat a, vfloat b) { a.m = _mm256_add_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat operator- (vfloat a, vfloat b) { a.m = _mm256_sub_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat operator* (vfloat a, vfloat b) { a.m = _mm256_mul_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat operator==(vfloat a, vfloat b) { a.m = _mm256_cmp_ps(a.m, b.m, _CMP_EQ_OQ); return a; }
SIMD_INLINE vfloat operator!=(vfloat a, vfloat b) { a.m = _mm256_cmp_ps(a.m, b.m, _CMP_NEQ_OQ); return a; }
SIMD_INLINE vfloat operator< (vfloat a, vfloat b) { a.m = _mm256_cmp_ps(a.m, b.m, _CMP_LT_OQ); return a; }
SIMD_INLINE vfloat operator> (vfloat a, vfloat b) { a.m = _mm256_cmp_ps(a.m, b.m, _CMP_GT_OQ); return a; }
SIMD_INLINE vfloat operator<=(vfloat a, vfloat b) { a.m = _mm256_cmp_ps(a.m, b.m, _CMP_LE_OQ); return a; }
SIMD_INLINE vfloat operator>=(vfloat a, vfloat b) { a.m = _mm256_cmp_ps(a.m, b.m, _CMP_GE_OQ); return a; }

SIMD_INLINE vfloat min(vfloat a, vfloat b) { a.m = _mm256_min_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat max(vfloat a, vfloat b) { a.m = _mm256_max_ps(a.m, b.m); return a; }

SIMD_INLINE vfloat round(vfloat v)
{
    return vfloat(_mm256_round_ps(v.m, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC));
}

SIMD_INLINE vint floatToInt(vfloat v) { return vint(_mm256_cvtps_epi32(v.m)); }

SIMD_INLINE vfloat intAsFloat(vint v) { return vfloat(_mm256_castsi256_ps(v.m)); }

SIMD_INLINE vint operator~ (vint a) { return vint(_mm256_xor_si256(a.m, _mm256_set1_epi32(-1))); }
SIMD_INLINE vint operator+ (vint a, vint b) { a.m = _mm256_add_epi32(a.m, b.m); return a; }
SIMD_INLINE vint operator- (vint a, vint b) { a.m = _mm256_sub_epi32(a.m, b.m); return a; }
SIMD_INLINE vint operator<(vint a, vint b) { a.m = _mm256_cmpgt_epi32(b.m, a.m); return a; }
SIMD_INLINE vint operator>(vint a, vint b) { a.m = _mm256_cmpgt_epi32(a.m, b.m); return a; }
SIMD_INLINE vint operator==(vint a, vint b) { a.m = _mm256_cmpeq_epi32(a.m, b.m); return a; }
SIMD_INLINE vint operator!=(vint a, vint b) { a.m = _mm256_cmpeq_epi32(a.m, b.m); return ~a; }
SIMD_INLINE vint min(vint a, vint b) { a.m = _mm256_min_epi32(a.m, b.m); return a; }
SIMD_INLINE vint max(vint a, vint b) { a.m = _mm256_max_epi32(a.m, b.m); return a; }

SIMD_INLINE void store(vfloat v, float* ptr) { _mm256_store_ps(ptr, v.m); }
SIMD_INLINE void store(vint v, int* ptr) { _mm256_store_si256((__m256i*)ptr, v.m); }


// "select", i.e. highbit(cond) ? b : a
SIMD_INLINE vfloat select(vfloat a, vfloat b, vfloat cond)
{
    return vfloat(_mm256_blendv_ps(a.m, b.m, cond.m));
}
SIMD_INLINE vint select(vint a, vint b, vint cond)
{
    return vint(_mm256_blendv_epi8(a.m, b.m, cond.m));
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
    static vfloat zero() { return vfloat(_mm_setzero_ps()); }
    __m128 m;
};

struct vint
{
    SIMD_INLINE vint() {}
    SIMD_INLINE explicit vint(const int *p) { m = _mm_load_si128((const __m128i*)p); }
    SIMD_INLINE explicit vint(int v) { m = _mm_set1_epi32(v); }
    SIMD_INLINE explicit vint(__m128i v) { m = v; }
    __m128i m;
};

SIMD_INLINE vfloat load1a(const float* p) { return vfloat(_mm_load_ps1(p)); }
SIMD_INLINE vfloat loada(const float* p) { return vfloat(_mm_load_ps(p)); }

SIMD_INLINE vfloat operator+ (vfloat a, vfloat b) { a.m = _mm_add_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat operator- (vfloat a, vfloat b) { a.m = _mm_sub_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat operator* (vfloat a, vfloat b) { a.m = _mm_mul_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat operator==(vfloat a, vfloat b) { a.m = _mm_cmpeq_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat operator!=(vfloat a, vfloat b) { a.m = _mm_cmpneq_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat operator< (vfloat a, vfloat b) { a.m = _mm_cmplt_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat operator> (vfloat a, vfloat b) { a.m = _mm_cmpgt_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat operator<=(vfloat a, vfloat b) { a.m = _mm_cmple_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat operator>=(vfloat a, vfloat b) { a.m = _mm_cmpge_ps(a.m, b.m); return a; }

SIMD_INLINE vfloat min(vfloat a, vfloat b) { a.m = _mm_min_ps(a.m, b.m); return a; }
SIMD_INLINE vfloat max(vfloat a, vfloat b) { a.m = _mm_max_ps(a.m, b.m); return a; }

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

SIMD_INLINE vint floatToInt(vfloat v) { return vint(_mm_cvtps_epi32(v.m)); }

SIMD_INLINE vfloat intAsFloat(vint v) { return vfloat(_mm_castsi128_ps(v.m)); }

SIMD_INLINE vint operator~ (vint a) { return vint(_mm_xor_si128(a.m, _mm_set1_epi32(-1))); }
SIMD_INLINE vint operator+ (vint a, vint b) { a.m = _mm_add_epi32(a.m, b.m); return a; }
SIMD_INLINE vint operator- (vint a, vint b) { a.m = _mm_sub_epi32(a.m, b.m); return a; }
SIMD_INLINE vint operator<(vint a, vint b) { a.m = _mm_cmplt_epi32(a.m, b.m); return a; }
SIMD_INLINE vint operator>(vint a, vint b) { a.m = _mm_cmpgt_epi32(a.m, b.m); return a; }
SIMD_INLINE vint operator==(vint a, vint b) { a.m = _mm_cmpeq_epi32(a.m, b.m); return a; }
SIMD_INLINE vint operator!=(vint a, vint b) { a.m = _mm_cmpeq_epi32(a.m, b.m); return ~a; }
SIMD_INLINE vint min(vint a, vint b) { a.m = _mm_min_epi32(a.m, b.m); return a; }
SIMD_INLINE vint max(vint a, vint b) { a.m = _mm_max_epi32(a.m, b.m); return a; }

SIMD_INLINE void store(vfloat v, float* ptr) { _mm_store_ps(ptr, v.m); }
SIMD_INLINE void store(vint v, int* ptr) { _mm_store_si128((__m128i*)ptr, v.m); }


// "select", i.e. highbit(cond) ? b : a
// on SSE4.1 and up this can be done easily via "blend" instruction;
// on older SSEs we have to do some hoops, see
// https://fgiesen.wordpress.com/2016/04/03/sse-mind-the-gap/
SIMD_INLINE vfloat select(vfloat a, vfloat b, vfloat cond)
{
#if ASTCENC_SSE >= 41
    a.m = _mm_blendv_ps(a.m, b.m, cond.m);
#else
    __m128 d = _mm_castsi128_ps(_mm_srai_epi32(_mm_castps_si128(cond.m), 31));
    a.m = _mm_or_ps(_mm_and_ps(d, b.m), _mm_andnot_ps(d, a.m));
#endif
    return a;
}
SIMD_INLINE vint select(vint a, vint b, vint cond)
{
#if ASTCENC_SSE >= 41
    return vint(_mm_blendv_epi8(a.m, b.m, cond.m));
#else
    __m128i d = _mm_srai_epi32(cond.m, 31);
    return vint(_mm_or_si128(_mm_and_si128(d, b.m), _mm_andnot_si128(d, a.m)));
#endif
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
    static vfloat zero() { return vfloat(0.0f); }
    float m;
};

struct vint
{
    SIMD_INLINE vint() {}
    SIMD_INLINE explicit vint(const int *p) { m = *p; }
    SIMD_INLINE explicit vint(int v) { m = v; }
    int m;
};

SIMD_INLINE vfloat load1a(const float* p) { return vfloat(*p); }
SIMD_INLINE vfloat loada(const float* p) { return vfloat(*p); }

SIMD_INLINE vfloat operator+ (vfloat a, vfloat b) { a.m = a.m + b.m; return a; }
SIMD_INLINE vfloat operator- (vfloat a, vfloat b) { a.m = a.m - b.m; return a; }
SIMD_INLINE vfloat operator* (vfloat a, vfloat b) { a.m = a.m * b.m; return a; }
SIMD_INLINE vfloat operator==(vfloat a, vfloat b) { a.m = a.m == b.m; return a; }
SIMD_INLINE vfloat operator!=(vfloat a, vfloat b) { a.m = a.m != b.m; return a; }
SIMD_INLINE vfloat operator< (vfloat a, vfloat b) { a.m = a.m < b.m; return a; }
SIMD_INLINE vfloat operator> (vfloat a, vfloat b) { a.m = a.m > b.m; return a; }
SIMD_INLINE vfloat operator<=(vfloat a, vfloat b) { a.m = a.m <= b.m; return a; }
SIMD_INLINE vfloat operator>=(vfloat a, vfloat b) { a.m = a.m >= b.m; return a; }

SIMD_INLINE vfloat min(vfloat a, vfloat b) { a.m = a.m < b.m ? a.m : b.m; return a; }
SIMD_INLINE vfloat max(vfloat a, vfloat b) { a.m = a.m > b.m ? a.m : b.m; return a; }

SIMD_INLINE vfloat round(vfloat v)
{
    return vfloat(std::floor(v.m + 0.5f));
}

SIMD_INLINE vint floatToInt(vfloat v) { return vint(v.m); }

SIMD_INLINE vfloat intAsFloat(vint v) { vfloat r; memcpy(&r.m, &v.m, 4); return r; }

SIMD_INLINE vint operator~ (vint a) { a.m = ~a.m; return a; }
SIMD_INLINE vint operator+ (vint a, vint b) { a.m = a.m + b.m; return a; }
SIMD_INLINE vint operator- (vint a, vint b) { a.m = a.m - b.m; return a; }
SIMD_INLINE vint operator<(vint a, vint b) { a.m = a.m < b.m; return a; }
SIMD_INLINE vint operator>(vint a, vint b) { a.m = a.m > b.m; return a; }
SIMD_INLINE vint operator==(vint a, vint b) { a.m = a.m == b.m; return a; }
SIMD_INLINE vint operator!=(vint a, vint b) { a.m = a.m != b.m; return a; }
SIMD_INLINE vint min(vint a, vint b) { a.m = a.m < b.m ? a.m : b.m; return a; }
SIMD_INLINE vint max(vint a, vint b) { a.m = a.m > b.m ? a.m : b.m; return a; }

SIMD_INLINE void store(vfloat v, float* ptr) { *ptr = v.m; }
SIMD_INLINE void store(vint v, int* ptr) { *ptr = v.m; }


// "select", i.e. highbit(cond) ? b : a
SIMD_INLINE vfloat select(vfloat a, vfloat b, vfloat cond)
{
    return cond.m ? b : a;
}
SIMD_INLINE vint select(vint a, vint b, vint cond)
{
    return cond.m ? b : a;
}


#endif // #ifdef ASTCENC_SIMD_ISA_SCALAR

#endif // #ifndef ASTC_VECMATHLIB_H_INCLUDED
