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
 * This module implements a variety of mathematical data types and library
 * functions used by the codec.
 */

#ifndef ASTC_MATHLIB_H_INCLUDED
#define ASTC_MATHLIB_H_INCLUDED

#include <cstdint>
#include <cmath>

#if ASTCENC_SSE != 0 || ASTCENC_AVX != 0
	#include <immintrin.h>
#endif


#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
  Fast math library; note that many of the higher-order functions in this set
  use approximations which are less accurate, but faster, than <cmath> standard
  library equivalents.

  Note: Many of these are not necessarily faster than simple C versions when
  used on a single scalar value, but are included for testing purposes as most
  have an option based on SSE intrinsics and therefore provide an obvious route
  to future vectorization.
============================================================================ */

// We support scalar versions of many maths functions which use SSE intrinsics
// as an "optimized" path, using just one lane from the SIMD hardware. In
// reality these are often slower than standard C due to setup and scheduling
// overheads, and the fact that we're not offsetting that cost with any actual
// vectorization.
//
// These variants are only included as a means to test that the accuracy of an
// SSE implementation would be acceptable before refactoring code paths to use
// an actual vectorized implementation which gets some advantage from SSE. It
// is therefore expected that the code will go *slower* with this macro
// set to 1 ...
#define USE_SCALAR_SSE 0

// These are namespaced to avoid colliding with C standard library functions.
namespace astc
{

/**
 * @brief Fast approximation of log2(x)
 *
 * This does not produce correct results for special cases such as
 * zero/inf/nan/denormal/negative inputs:
 *
 *   * Any negative, zero, or denormal will get clamped to smallest-normal,
 *     resulting in a logarithm of -126.
 *   * +Inf and +NaN get treated as an extension of largest-finite values,
 *     which should result in a logarithm value between 128 and 129.
 */
float log2(float val);

/**
 * @brief Fast approximation of atan2.
 *
 * TODO: This implementation is reasonably accurate and a lot faster than the
 * standard library, but quite branch heavy which makes it difficult to
 * vectorize effectively. If we need to vectorize in future, consider using a
 * different approximation algorithm.
 *
 * @param y The proportion of the Y coordinate.
 * @param x The proportion of the X coordinate.
 *
 * @return The approximation of atan2().
 */
float atan2(float y, float x);

/**
 * @brief SP float absolute value.
 *
 * @param val The value to make absolute.
 *
 * @return The absolute value.
 */
static inline float fabs(float val)
{
#if (ASTCENC_SSE >= 20) && USE_SCALAR_SSE
	static const union {
		uint32_t u[4];
		__m128 v;
	} mask = { { 0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff } };
	return _mm_cvtss_f32(_mm_and_ps(_mm_set_ss(val), mask.v));
#else
	return std::fabs(val);
#endif
}

/**
 * @brief SP float min.
 *
 * Note: GCC versions prior to 7.x assume input arguments to the intrinsics are
 * commutative, which isn't true for NaN handling, so it's best not to rely on
 * argument order unless you are very sure about the compiler ...
 *
 * @param valA The first value to compare.
 * @param valB The second value to compare.
 *
 * @return The smallest value.
 */
static inline float fmin(float p, float q)
{
#if (ASTCENC_SSE >= 20) && USE_SCALAR_SSE
	return _mm_cvtss_f32(_mm_min_ss(_mm_set_ss(p),_mm_set_ss(q)));
#else
	return p < q ? p : q;
#endif
}

/**
 * @brief SP float max.
 *
 * Note: GCC versions prior to 7.x assume input arguments to the intrinsics are
 * commutative, which isn't true for NaN handling, so it's best not to rely on
 * argument order unless you are very sure about the compiler ...
 *
 * @param valA The first value to compare.
 * @param valB The second value to compare.
 *
 * @return The largest value.
 */
static inline float fmax(float p, float q)
{
#if (ASTCENC_SSE >= 20) && USE_SCALAR_SSE
    return _mm_cvtss_f32(_mm_max_ss(_mm_set_ss(p),_mm_set_ss(q)));
#else
    return q < p ? p : q;
#endif
}

/**
 * @brief Test if a float value is a nan.
 *
 * @param val The value test.
 *
 * @return Zero is not a NaN, non-zero otherwise.
 */
static inline int isnan(float val)
{
	return val != val;
}

/**
 * @brief Clamp a float value between 0.0f and 1.0f.
 *
 * NaNs are turned into 0.0f.
 *
 * @param val The value clamp.
 *
 * @return The clamped value.
 */
static inline float clamp1f(float val)
{
	// Do not reorder these, correct NaN handling relies on the fact that
	// any comparison with NaN returns false so will fall-though to the 0.0f.
	if (val > 1.0f) return 1.0f;
	if (val > 0.0f) return val;
	return 0.0f;
}

/**
 * @brief Clamp a float value between 0.0f and 255.0f.
 *
 * NaNs are turned into 0.0f.
 *
 * @param val The value clamp.
 *
 * @return The clamped value.
 */
static inline float clamp255f(float val)
{
	// Do not reorder these, correct NaN handling relies on the fact that
	// any comparison with NaN returns false so will fall-though to the 0.0f.
	if (val > 255.0f) return 255.0f;
	if (val > 0.0f) return val;
	return 0.0f;
}

/**
 * @brief Clamp a float value between 0.0f and 65504.0f.
 *
 * NaNs are turned into 0.0f.
 *
 * @param val The value to clamp
 *
 * @return The clamped value
 */
static inline float clamp64Kf(float val)
{
	// Do not reorder these, correct NaN handling relies on the fact that
	// any comparison with NaN returns false so will fall-though to the 0.0f.
	if (val > 65504.0f) return 65504.0f;
	if (val > 0.0f) return val;
	return 0.0f;
}

/**
 * @brief Clamp an integer between two specified limits.
 *
 * @param val The value clamp.
 *
 * @return The clamped value.
 */
static inline int clampi(int val, int low, int high)
{
	if (val < low) return low;
	if (val > high) return high;
	return val;
}

/**
 * @brief SP float round-to-nearest.
 *
 * @param val The value to round.
 *
 * @return The rounded value.
 */
static inline float flt_rte(float val)
{
#if (ASTCENC_SSE >= 42)
	const int flag = _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC;
	__m128 tmp = _mm_set_ss(val);
	tmp = _mm_round_ss(tmp, tmp, flag);
	return _mm_cvtss_f32(tmp);
#else
	return std::floor(val + 0.5f);
#endif
}

/**
 * @brief SP float round-down.
 *
 * @param val The value to round.
 *
 * @return The rounded value.
 */
static inline float flt_rd(float val)
{
#if (ASTCENC_SSE >= 42)
	const int flag = _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC;
	__m128 tmp = _mm_set_ss(val);
	tmp = _mm_round_ss(tmp, tmp, flag);
	return _mm_cvtss_f32(tmp);
#else
	return std::floor(val);
#endif
}

/**
 * @brief SP float round-to-nearest and convert to integer.
 *
 * @param val The value to round.
 *
 * @return The rounded value.
 */
static inline int flt2int_rtn(float val)
{
#if (ASTCENC_SSE >= 42) && USE_SCALAR_SSE
	return _mm_cvt_ss2si(_mm_set_ss(val));
#else
	return (int)(val + 0.5f);
#endif
}

/**
 * @brief SP float round down and convert to integer.
 *
 * @param val The value to round.
 *
 * @return The rounded value.
 */
static inline int flt2int_rd(float val)
{
#if (ASTCENC_SSE >= 42) && USE_SCALAR_SSE
	return _mm_cvt_ss2si(_mm_set_ss(val));
#else
	return (int)(val);
#endif
}

/**
 * @brief Population bit count.
 *
 * @param val The value to count.
 *
 * @return The number of 1 bits.
 */
static inline int popcount(uint64_t p)
{
#if ASTCENC_POPCNT >= 1
	return (int)_mm_popcnt_u64(p);
#else
	uint64_t mask1 = 0x5555555555555555ULL;
	uint64_t mask2 = 0x3333333333333333ULL;
	uint64_t mask3 = 0x0F0F0F0F0F0F0F0FULL;
	p -= (p >> 1) & mask1;
	p = (p & mask2) + ((p >> 2) & mask2);
	p += p >> 4;
	p &= mask3;
	p *= 0x0101010101010101ULL;
	p >>= 56;
	return (int)p;
#endif
}

/**
 * @brief Fast approximation of 1.0 / sqrt(val).
 *
 * @param val The input value.
 *
 * @return The approximated result.
 */
static inline float rsqrt(float val)
{
#if (ASTCENC_SSE >= 20) && USE_SCALAR_SSE
	// FIXME: setting val = 99 causes a crash, which it really shouldn't.
	return _mm_cvtss_f32(_mm_rsqrt_ss(_mm_set_ss(val)));
#else
	return 1.0f / std::sqrt(val);
#endif
}

/**
 * @brief Fast approximation of sqrt(val).
 *
 * @param val The input value.
 *
 * @return The approximated result.
 */
static inline float sqrt(float val)
{
#if (ASTCENC_SSE >= 20) && USE_SCALAR_SSE
	return 1.0f * astc::rsqrt(val);
#else
	return std::sqrt(val);
#endif
}

/**
 * @brief Fast approximation of 1.0 / val.
 *
 * @param val The input value.
 *
 * @return The approximated result.
 */
static inline float recip(float val)
{
#if (ASTCENC_SSE >= 20) && USE_SCALAR_SSE
	return _mm_cvtss_f32(_mm_rcp_ss(_mm_set_ss(val)));
#else
	return 1.0f / val;
#endif
}

/**
 * @brief Log base 2, linearized from 2^-14.
 *
 * @param val The value to log2.
 *
 * @return The approximated result.
 */
static inline float xlog2(float val)
{
	if (val >= 0.00006103515625f)
	{
		return astc::log2(val);
	}

	// Linearized region
	return -15.44269504088896340735f + val * 23637.11554992477646609062f;
}

/**
 * @brief Initialize the seed structure for a random number generator.
 *
 * Important note: For the purposes of ASTC we want sets of random numbers to
 * use the codec, but we want the same seed value across instances and threads
 * to ensure that image output is stable across compressor runs and across
 * platforms. Every PRNG created by this call will therefore return the same
 * sequence of values ...
 *
 * @param state The state structure to initialize.
 */
void rand_init(uint64_t state[2]);

/**
 * @brief Return the next random number from the generator.
 *
 * This RNG is an implementation of the "xoroshoro-128+ 1.0" PRNG, based on the
 * public-domain implementation given by David Blackman & Sebastiano Vigna at
 * http://vigna.di.unimi.it/xorshift/xoroshiro128plus.c
 *
 * @param state The state structure to use/update.
 */
uint64_t rand(uint64_t state[2]);

}

/* ============================================================================
  Utility vector template classes with basic operations
============================================================================ */

template <typename T> class vtype2
{
public:
	T x, y;
	vtype2() {}
	vtype2(T p, T q)         : x(p),   y(q)   {}
	vtype2(const vtype2 & p) : x(p.x), y(p.y) {}
	vtype2 &operator =(const vtype2 &s) {
		this->x = s.x;
		this->y = s.y;
		return *this;
	}
};

template <typename T> class vtype3
{
public:
	T x, y, z;
	vtype3() {}
	vtype3(T p, T q, T r)    : x(p),   y(q),   z(r)   {}
	vtype3(const vtype3 & p) : x(p.x), y(p.y), z(p.z) {}
	vtype3 &operator =(const vtype3 &s) {
		this->x = s.x;
		this->y = s.y;
		this->z = s.z;
		return *this;
	}
};

template <typename T> class vtype4
{
public:
	T x, y, z, w;
	vtype4() {}
	vtype4(T p, T q, T r, T s) : x(p),   y(q),   z(r),   w(s)   {}
	vtype4(const vtype4 & p)   : x(p.x), y(p.y), z(p.z), w(p.w) {}
	vtype4 &operator =(const vtype4 &s) {
		this->x = s.x;
		this->y = s.y;
		this->z = s.z;
		this->w = s.w;
		return *this;
	}
};

typedef vtype2<float>        float2;
typedef vtype3<float>        float3;
typedef vtype4<float>        float4;
typedef vtype3<int>          int3;
typedef vtype4<int>          int4;
typedef vtype4<unsigned int> uint4;

static inline float2  operator+(float2 p,  float2 q)   { return float2(  p.x + q.x, p.y + q.y ); }
static inline float3  operator+(float3 p,  float3 q)   { return float3(  p.x + q.x, p.y + q.y, p.z + q.z ); }
static inline float4  operator+(float4 p,  float4 q)   { return float4(  p.x + q.x, p.y + q.y, p.z + q.z, p.w + q.w ); }
static inline int4    operator+(int4 p,    int4 q)     { return int4(    p.x + q.x, p.y + q.y, p.z + q.z, p.w + q.w ); }
static inline uint4   operator+(uint4 p,   uint4 q)    { return uint4(   p.x + q.x, p.y + q.y, p.z + q.z, p.w + q.w ); }

static inline float2  operator-(float2 p,  float2 q)   { return float2(  p.x - q.x, p.y - q.y ); }
static inline float3  operator-(float3 p,  float3 q)   { return float3(  p.x - q.x, p.y - q.y, p.z - q.z ); }
static inline float4  operator-(float4 p,  float4 q)   { return float4(  p.x - q.x, p.y - q.y, p.z - q.z, p.w - q.w ); }
static inline int4    operator-(int4 p,    int4 q)     { return int4(    p.x - q.x, p.y - q.y, p.z - q.z, p.w - q.w ); }
static inline uint4   operator-(uint4 p,   uint4 q)    { return uint4(   p.x - q.x, p.y - q.y, p.z - q.z, p.w - q.w ); }

static inline float2  operator*(float2 p,  float2 q)   { return float2(  p.x * q.x, p.y * q.y ); }
static inline float3  operator*(float3 p,  float3 q)   { return float3(  p.x * q.x, p.y * q.y, p.z * q.z ); }
static inline float4  operator*(float4 p,  float4 q)   { return float4(  p.x * q.x, p.y * q.y, p.z * q.z, p.w * q.w ); }
static inline int4    operator*(int4 p,    int4 q)     { return int4(    p.x * q.x, p.y * q.y, p.z * q.z, p.w * q.w ); }
static inline uint4   operator*(uint4 p,   uint4 q)    { return uint4(   p.x * q.x, p.y * q.y, p.z * q.z, p.w * q.w ); }

static inline float2  operator*(float2 p,  float q)    { return float2(  p.x * q, p.y * q ); }
static inline float3  operator*(float3 p,  float q)    { return float3(  p.x * q, p.y * q, p.z * q ); }
static inline float4  operator*(float4 p,  float q)    { return float4(  p.x * q, p.y * q, p.z * q, p.w * q ); }
static inline int4    operator*(int4 p,    int q)      { return int4(    p.x * q, p.y * q, p.z * q, p.w * q ); }
static inline uint4   operator*(uint4 p,   uint32_t q) { return uint4(   p.x * q, p.y * q, p.z * q, p.w * q ); }

static inline float2  operator*(float p,    float2 q)  { return q * p; }
static inline float3  operator*(float p,    float3 q)  { return q * p; }
static inline float4  operator*(float p,    float4 q)  { return q * p; }
static inline int4    operator*(int p,      int4 q)    { return q * p; }
static inline uint4   operator*(uint32_t p, uint4 q)   { return q * p; }

static inline float dot(float2 p, float2 q)  { return p.x * q.x + p.y * q.y; }
static inline float dot(float3 p, float3 q)  { return p.x * q.x + p.y * q.y + p.z * q.z; }
static inline float dot(float4 p, float4 q)  { return p.x * q.x + p.y * q.y + p.z * q.z + p.w * q.w; }

static inline float2 normalize(float2 p) { return p * astc::rsqrt(dot(p, p)); }
static inline float3 normalize(float3 p) { return p * astc::rsqrt(dot(p, p)); }
static inline float4 normalize(float4 p) { return p * astc::rsqrt(dot(p, p)); }

#ifndef MIN
	#define MIN(x,y) ((x)<(y)?(x):(y))
#endif

#ifndef MAX
	#define MAX(x,y) ((x)>(y)?(x):(y))
#endif

/* ============================================================================
  Softfloat library with fp32 and fp16 conversion functionality.
============================================================================ */
typedef union if32_
{
	uint32_t u;
	int32_t s;
	float f;
} if32;

uint32_t clz32(uint32_t p);

/*	sized soft-float types. These are mapped to the sized integer
    types of C99, instead of C's floating-point types; this is because
    the library needs to maintain exact, bit-level control on all
    operations on these data types. */
typedef uint16_t sf16;
typedef uint32_t sf32;

/* the five rounding modes that IEEE-754r defines */
typedef enum
{
	SF_UP = 0,				/* round towards positive infinity */
	SF_DOWN = 1,			/* round towards negative infinity */
	SF_TOZERO = 2,			/* round towards zero */
	SF_NEARESTEVEN = 3,		/* round toward nearest value; if mid-between, round to even value */
	SF_NEARESTAWAY = 4		/* round toward nearest value; if mid-between, round away from zero */
} roundmode;

/* narrowing float->float conversions */
sf16 sf32_to_sf16(sf32, roundmode);

/* widening float->float conversions */
sf32 sf16_to_sf32(sf16);

sf16 float_to_sf16(float, roundmode);

float sf16_to_float(sf16);


/*********************************
  Declaration of line types
*********************************/
// parametric line, 2D: The line is given by line = a + b*t.

struct line2
{
	float2 a;
	float2 b;
};

// parametric line, 3D
struct line3
{
	float3 a;
	float3 b;
};

struct line4
{
	float4 a;
	float4 b;
};


struct processed_line2
{
	float2 amod;
	float2 bs;
	float2 bis;
};

struct processed_line3
{
	float3 amod;
	float3 bs;
	float3 bis;
};

struct processed_line4
{
	float4 amod;
	float4 bs;
	float4 bis;
};

#endif
