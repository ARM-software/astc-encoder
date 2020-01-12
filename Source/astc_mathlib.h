// ----------------------------------------------------------------------------
//  This confidential and proprietary software may be used only as authorised
//  by a licensing agreement from Arm Limited.
//      (C) COPYRIGHT 2011-2020 Arm Limited, ALL RIGHTS RESERVED
//  The entire notice above must be reproduced on all authorised copies and
//  copies may only be made to the extent permitted by a licensing agreement
//  from Arm Limited.
// ----------------------------------------------------------------------------

/*
 * This module implements a variety of mathematical data types and library
 * functions used by the codec.
 */

#ifndef ASTC_MATHLIB_H_INCLUDED
#define ASTC_MATHLIB_H_INCLUDED

#include <cstdint>
#include <cmath>

#include <immintrin.h>

#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif

// **********************************************************************
// vector data types and basic add/subtract/multiply functions on them
// **********************************************************************

template <typename T> class vtype2
{
public:
	T x, y;
	vtype2() {};
	vtype2(T p, T q)         : x(p),   y(q)   {};
	vtype2(const vtype2 & p) : x(p.x), y(p.y) {};
};

template <typename T> class vtype3
{
public:
	T x, y, z;
	vtype3() {};
	vtype3(T p, T q, T r)    : x(p),   y(q),   z(r)   {};
	vtype3(const vtype3 & p) : x(p.x), y(p.y), z(p.z) {};
};

template <typename T> class vtype4
{
public:
	T x, y, z, w;
	vtype4() {};
	vtype4(T p, T q, T r, T s) : x(p),   y(q),   z(r),   w(s)   {};
	vtype4(const vtype4 & p)   : x(p.x), y(p.y), z(p.z), w(p.w) {};
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

static inline float2 normalize(float2 p) { return p * (1.0f / sqrtf(dot(p,p))); }
static inline float3 normalize(float3 p) { return p * (1.0f / sqrtf(dot(p,p))); }
static inline float4 normalize(float4 p) { return p * (1.0f / sqrtf(dot(p,p))); }

#ifndef MIN
	#define MIN(x,y) ((x)<(y)?(x):(y))
#endif

#ifndef MAX
	#define MAX(x,y) ((x)>(y)?(x):(y))
#endif

#define astc_isnan(p) ((p)!=(p))

// Clamp an input value to [0,1]; NaN is turned into 0.
static inline float clamp01(float val)
{
	if (val > 1.0f) return 1.0f;
	if (val > 0.0f) return val;
	return 0.0f;
}

static inline int iclamp(int val, int low, int high)
{
	if (val < low) return low;
	if (val > high) return high;
	return val;
}

static inline float srgb_transform(float val)
{
	if (val <= 0.04045f) return val * (1.0f / 12.92f);
	if (val <= 1) return powf((val + 0.055f) * (1.0f / 1.055f), 2.4f);
	return val;
}


/*******************************************************
  softfloat library: conversion between FP16 and FP32
*******************************************************/

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

struct mat2
{
	float2 v[2];
};

struct mat4
{
	float4 v[4];
};


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

float determinant(mat2 p);

float2 transform(mat2 p, float2 q);
float4 transform(mat4 p, float4 q);

mat2 invert(mat2 p);
mat4 invert(mat4 p);

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

/**
 * @brief SP float absolute value.
 *
 * @param val The value to make absolute.
 *
 * @return The absolute value.
 */
static inline float fabs(float val)
{
#if ASTC_SSE >= 20
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
 * @brief SP float round-to-nearest.
 *
 * @param val The value to round.
 *
 * @return The rounded value.
 */
static inline float flt_rte(float val)
{
// round to integer, round-to-nearest
#if ASTC_SSE >= 42
	const int flag = _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC;
	__m128 tmp = _mm_set_ss(val);
	tmp = _mm_round_ss(tmp, tmp, flag);
	return _mm_cvtss_f32(tmp);
#else
	return floorf(val + 0.5f);
#endif
}

/**
 * @brief SP float round-to-nearest and convert to integer.
 *
 * @param val The value to round.
 *
 * @return The rounded value.
 */
static inline int flt2int_rte(float val)
{ // convert to integer, round-to-nearest
#if ASTC_SSE >= 20
	return _mm_cvt_ss2si(_mm_set_ss(val));
#else
	return (int)(floorf(val + 0.5f));
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
#if ASTC_SSE >= 42
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
static inline float atan2(float y, float x)
{
	const float PI = (float)M_PI;
	const float PI_2 = PI / 2.f;

	// Handle the discontinuity at x == 0
	if (x == 0.0f)
	{
		if (y > 0.0f)
		{
			return PI_2;
		}
		else if (y == 0.0f)
		{
			return 0.0f;
		}
		return -PI_2;
	}

	float z = y / x;
	float z2 = z * z;
	if (std::fabs(z) < 1.0f)
	{
		float atan = z / (1.0f + (0.28f * z2));
		if (x < 0.0f)
		{
			if (y < 0.0f)
			{
				return atan - PI;
			}
			else
			{
				return atan + PI;
			}
		}
		return atan;
	}
	else
	{
		float atan = PI_2 - (z / (z2 + 0.28f));
		if (y < 0.0f)
		{
			return atan - PI;
		}
		else
		{
			return atan;
		}
	}
}

}

#endif
