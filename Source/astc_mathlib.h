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

#endif
