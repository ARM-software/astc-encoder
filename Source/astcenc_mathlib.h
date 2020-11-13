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

#include <cassert>
#include <cstdint>
#include <cmath>

#if ASTCENC_SSE != 0 || ASTCENC_AVX != 0
	#include <immintrin.h>
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

// These are namespaced to avoid colliding with C standard library functions.
namespace astc
{

static const float PI          = 3.14159265358979323846f;
static const float PI_OVER_TWO = 1.57079632679489661923f;

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
 * @brief SP float absolute value.
 *
 * @param val The value to make absolute.
 *
 * @return The absolute value.
 */
static inline float fabs(float val)
{
	return std::fabs(val);
}

/**
 * @brief SP float min.
 *
 * @param valA The first value to compare.
 * @param valB The second value to compare.
 *
 * @return The smallest value.
 */
static inline float fmin(float p, float q)
{
	return p < q ? p : q;
}

/**
 * @brief SP float max.
 *
 * @param valA The first value to compare.
 * @param valB The second value to compare.
 *
 * @return The largest value.
 */
static inline float fmax(float p, float q)
{
	return q < p ? p : q;
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
 * @brief Clamp a value value between mn and mx
 *
 * For floats, NaNs are turned into mn.
 *
 * @param val The value clamp.
 * @param mn  The min value (inclusive).
 * @param mx  The max value (inclusive).
 *
 * @return The clamped value.
 */
template<typename T>
inline T clamp(T val, T mn, T mx)
{
	// Do not reorder; correct NaN handling relies on the fact that comparison
	// with NaN returns false and will fall-though to the "min" value.
	if (val > mx) return mx;
	if (val > mn) return val;
	return mn;
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
	return std::floor(val + 0.5f);
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
	return std::floor(val);
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

	return (int)(val + 0.5f);
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
	return (int)(val);
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
	return 1.0f / std::sqrt(val);
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
	return std::sqrt(val);
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
	// Data storage
	T r, g;

	// Default constructor
	vtype2() {}

	// Initialize from 1 scalar
	vtype2(T p) : r(p), g(p) {}

	// Initialize from N scalars
	vtype2(T p, T q) : r(p), g(q) {}

	// Initialize from another vector
	vtype2(const vtype2 & p) : r(p.r), g(p.g) {}

	// Assignment operator
	vtype2& operator=(const vtype2 &s) {
		this->r = s.r;
		this->g = s.g;
		return *this;
	}
};

// Vector by vector addition
template <typename T>
vtype2<T> operator+(vtype2<T> p, vtype2<T> q) {
	return vtype2<T> { p.r + q.r, p.g + q.g };
}

// Vector by vector subtraction
template <typename T>
vtype2<T> operator-(vtype2<T> p, vtype2<T> q) {
	return vtype2<T> { p.r - q.r, p.g - q.g };
}

// Vector by vector multiplication operator
template <typename T>
vtype2<T> operator*(vtype2<T> p, vtype2<T> q) {
	return vtype2<T> { p.r * q.r, p.g * q.g };
}

// Vector by scalar multiplication operator
template <typename T>
vtype2<T> operator*(vtype2<T> p, T q) {
	return vtype2<T> { p.r * q, p.g * q };
}

// Scalar by vector multiplication operator
template <typename T>
vtype2<T> operator*(T p, vtype2<T> q){
	return vtype2<T> { p * q.r, p * q.g };
}

template <typename T> class vtype3
{
public:
	// Data storage
	T r, g, b;

	// Default constructor
	vtype3() {}

	// Initialize from 1 scalar
	vtype3(T p) : r(p), g(p), b(p) {}

	// Initialize from N scalars
	vtype3(T p, T q, T s) : r(p), g(q), b(s) {}

	// Initialize from another vector
	vtype3(const vtype3 & p) : r(p.r), g(p.g), b(p.b) {}

	// Assignment operator
	vtype3& operator=(const vtype3 &s) {
		this->r = s.r;
		this->g = s.g;
		this->b = s.b;
		return *this;
	}
};

// Vector by vector addition
template <typename T>
vtype3<T> operator+(vtype3<T> p, vtype3<T> q) {
	return vtype3<T> { p.r + q.r, p.g + q.g, p.b + q.b };
}

// Vector by vector subtraction
template <typename T>
vtype3<T> operator-(vtype3<T> p, vtype3<T> q) {
	return vtype3<T> { p.r - q.r, p.g - q.g, p.b - q.b };
}

// Vector by vector multiplication operator
template <typename T>
vtype3<T> operator*(vtype3<T> p, vtype3<T> q) {
	return vtype3<T> { p.r * q.r, p.g * q.g, p.b * q.b };
}

// Vector by scalar multiplication operator
template <typename T>
vtype3<T> operator*(vtype3<T> p, T q) {
	return vtype3<T> { p.r * q, p.g * q, p.b * q };
}

// Scalar by vector multiplication operator
template <typename T>
vtype3<T> operator*(T p, vtype3<T> q){
	return vtype3<T> { p * q.r, p * q.g, p * q.b };
}

template <typename T> class alignas(16) vtype4
{
public:
	// Data storage
	T r, g, b, a;

	// Default constructor
	vtype4() {}

	// Initialize from 1 scalar
	vtype4(T p) : r(p), g(p), b(p), a(p) {}

	// Initialize from N scalars
	vtype4(T p, T q, T s, T t) : r(p), g(q), b(s), a(t) {}

	// Initialize from another vector
	vtype4(const vtype4 & p) : r(p.r), g(p.g), b(p.b), a(p.a) {}

	// Assignment operator
	vtype4& operator=(const vtype4 &s) {
		this->r = s.r;
		this->g = s.g;
		this->b = s.b;
		this->a = s.a;
		return *this;
	}
};

// Vector by vector addition
template <typename T>
vtype4<T> operator+(vtype4<T> p, vtype4<T> q) {
	return vtype4<T> { p.r + q.r, p.g + q.g, p.b + q.b, p.a + q.a };
}

// Vector by vector subtraction
template <typename T>
vtype4<T> operator-(vtype4<T> p, vtype4<T> q) {
	return vtype4<T> { p.r - q.r, p.g - q.g, p.b - q.b, p.a - q.a };
}

// Vector by vector multiplication operator
template <typename T>
vtype4<T> operator*(vtype4<T> p, vtype4<T> q) {
	return vtype4<T> { p.r * q.r, p.g * q.g, p.b * q.b, p.a * q.a };
}

// Vector by scalar multiplication operator
template <typename T>
vtype4<T> operator*(vtype4<T> p, T q) {
	return vtype4<T> { p.r * q, p.g * q, p.b * q, p.a * q };
}

// Scalar by vector multiplication operator
template <typename T>
vtype4<T> operator*(T p, vtype4<T> q){
	return vtype4<T> { p * q.r, p * q.g, p * q.b, p * q.a };
}

typedef vtype2<float>        float2;
typedef vtype3<float>        float3;
typedef vtype4<float>        float4;
typedef vtype3<int>          int3;
typedef vtype4<int>          int4;
typedef vtype4<unsigned int> uint4;

static inline float dot(float2 p, float2 q)  { return p.r * q.r + p.g * q.g; }
static inline float dot(float3 p, float3 q)  { return p.r * q.r + p.g * q.g + p.b * q.b; }
static inline float dot(float4 p, float4 q)  {
#if (ASTCENC_SSE >= 42) && (ASTCENC_ISA_INVARIANCE == 0)
	__m128 pv = _mm_load_ps((float*)&p);
	__m128 qv = _mm_load_ps((float*)&q);
	__m128 t  = _mm_dp_ps(pv, qv, 0xFF);
	return _mm_cvtss_f32(t);
#else
	return p.r * q.r + p.g * q.g + p.b * q.b  + p.a * q.a;
#endif
}

static inline float2 normalize(float2 p) { return p * astc::rsqrt(dot(p, p)); }
static inline float3 normalize(float3 p) { return p * astc::rsqrt(dot(p, p)); }
static inline float4 normalize(float4 p) { return p * astc::rsqrt(dot(p, p)); }

static inline float4 sqrt(float4 p) {
	float4 r;
#if ASTCENC_SSE >= 20
	__m128 pv = _mm_load_ps((float*)&p);
	__m128 t  = _mm_sqrt_ps(pv);
	_mm_store_ps((float*)&r, t);
#else
	r.r = std::sqrt(p.r);
	r.g = std::sqrt(p.g);
	r.b = std::sqrt(p.b);
	r.a = std::sqrt(p.a);
#endif
	return r;
}

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
// parametric line, 2D: The line is given by line = a + b * t.

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
