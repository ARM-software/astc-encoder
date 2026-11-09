// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2011-2026 Arm Limited
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

#if !defined(ASTC_MATHLIB_H_INCLUDED)
#define ASTC_MATHLIB_H_INCLUDED

#include <cassert>
#include <bit>
#include <cstdint>
#include <cmath>
#include <cstring>

#if !defined(ASTCENC_POPCNT)
  #if defined(__POPCNT__)
    #define ASTCENC_POPCNT 1
  #else
    #define ASTCENC_POPCNT 0
  #endif
#endif

#if !defined(ASTCENC_F16C)
  #if defined(__F16C__)
    #define ASTCENC_F16C 1
  #else
    #define ASTCENC_F16C 0
  #endif
#endif

#if !defined(ASTCENC_SSE)
  #if defined(__SSE4_2__)
    #define ASTCENC_SSE 42
  #elif defined(__SSE4_1__)
    #define ASTCENC_SSE 41
  #elif defined(__SSE2__) || (defined(_M_AMD64) && !defined(_M_ARM64EC))
    #define ASTCENC_SSE 20
  #else
    #define ASTCENC_SSE 0
  #endif
#endif

// AVX levels: 1 = AVX, 2 = AVX2, 3 = AVX512 8-wide, 4 = AVX512 16-wide
#if !defined(ASTCENC_AVX)
  #if defined(__AVX512F__) && defined(__AVX512VBMI__)
    #define ASTCENC_AVX 4
    #define ASTCENC_X86_GATHERS 1
  #elif defined(__AVX2__)
    #define ASTCENC_AVX 2
    #define ASTCENC_X86_GATHERS 1
  #elif defined(__AVX__)
    #define ASTCENC_AVX 1
    #define ASTCENC_X86_GATHERS 1
  #else
    #define ASTCENC_AVX 0
  #endif
#endif

#if !defined(ASTCENC_NEON)
  #if defined(__aarch64__) || defined(_M_ARM64) || defined(_M_ARM64EC)
    #define ASTCENC_NEON 1
  #else
    #define ASTCENC_NEON 0
  #endif
#endif

#if !defined(ASTCENC_SVE)
  #if defined(__ARM_FEATURE_SVE)
    #if defined(__ARM_FEATURE_SVE_BITS) && __ARM_FEATURE_SVE_BITS == 256
      #define ASTCENC_SVE 8
    // Auto-detected SVE can only assume vector width of 4 is available, but must also allow for
    // hardware being longer and so all use of intrinsics must explicitly use predicate masks to
    // limit to 4-wide
    #else
      #define ASTCENC_SVE 4
    #endif
    #else
    #define ASTCENC_SVE 0
  #endif
#endif

// Force vector-sized SIMD alignment
#if ASTCENC_AVX >= 4
  #define ASTCENC_VECALIGN 64
#elif ASTCENC_AVX || ASTCENC_SVE == 8
  #define ASTCENC_VECALIGN 32
#elif ASTCENC_SSE || ASTCENC_NEON || ASTCENC_SVE == 4
  #define ASTCENC_VECALIGN 16
// Use default alignment for non-SIMD builds
#else
  #define ASTCENC_VECALIGN 0
#endif

// C++11 states that alignas(0) should be ignored but GCC doesn't do this on some versions, so
// workaround and avoid emitting alignas(0)
#if ASTCENC_VECALIGN > 0
	#define ASTCENC_ALIGNAS alignas(ASTCENC_VECALIGN)
#else
	#define ASTCENC_ALIGNAS
#endif

#if ASTCENC_SSE != 0 || ASTCENC_AVX != 0 || ASTCENC_POPCNT != 0
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

// These are namespaced to avoid colliding with C standard library functions
namespace astc
{

static const float PI          = 3.14159265358979323846f;
static const float PI_OVER_TWO = 1.57079632679489661923f;

/**
 * @brief SP float absolute value.
 *
 * @param v   The value to make absolute.
 *
 * @return The absolute value.
 */
static inline float fabs(float v)
{
	return std::fabs(v);
}

/**
 * @brief Test if a float value is a NaN.
 *
 * @param v   The value to test.
 *
 * @return @c true if the value is a NaN, @c false otherwise.
 */
static inline bool isnan(float v)
{
	return v != v;
}

/**
 * @brief Return the minimum of two values.
 *
 * For floats, an unordered comparison selects @c q.
 *
 * @tparam T   The value type.
 *
 * @param  p   The first value to compare.
 * @param  q   The second value to compare.
 *
 * @return The smallest value.
 */
template<typename T>
static inline T min(T p, T q)
{
	return p < q ? p : q;
}

/**
 * @brief Return the minimum of three values.
 *
 * For floats, each pairwise min selects its second operand when the
 * comparison is unordered.
 *
 * @tparam T   The value type.
 *
 * @param  p   The first value to compare.
 * @param  q   The second value to compare.
 * @param  r   The third value to compare.
 *
 * @return The smallest value.
 */
template<typename T>
static inline T min(T p, T q, T r)
{
	return min(min(p, q), r);
}

/**
 * @brief Return the minimum of four values.
 *
 * For floats, each pairwise min selects its second operand when the
 * comparison is unordered.
 *
 * @tparam T   The value type.
 *
 * @param  p   The first value to compare.
 * @param  q   The second value to compare.
 * @param  r   The third value to compare.
 * @param  s   The fourth value to compare.
 *
 * @return The smallest value.
 */
template<typename T>
static inline T min(T p, T q, T r, T s)
{
	return min(min(p, q), min(r, s));
}

/**
 * @brief Return the maximum of two values.
 *
 * For floats, an unordered comparison selects @c q.
 *
 * @tparam T   The value type.
 *
 * @param  p   The first value to compare.
 * @param  q   The second value to compare.
 *
 * @return The largest value.
 */
template<typename T>
static inline T max(T p, T q)
{
	return p > q ? p : q;
}

/**
 * @brief Return the maximum of three values.
 *
 * For floats, each pairwise max selects its second operand when the
 * comparison is unordered.
 *
 * @tparam T   The value type.
 *
 * @param  p   The first value to compare.
 * @param  q   The second value to compare.
 * @param  r   The third value to compare.
 *
 * @return The largest value.
 */
template<typename T>
static inline T max(T p, T q, T r)
{
	return max(max(p, q), r);
}

/**
 * @brief Return the maximum of four values.
 *
 * For floats, each pairwise max selects its second operand when the
 * comparison is unordered.
 *
 * @tparam T   The value type.
 *
 * @param  p   The first value to compare.
 * @param  q   The second value to compare.
 * @param  r   The third value to compare.
 * @param  s   The fourth value to compare.
 *
 * @return The largest value.
 */
template<typename T>
static inline T max(T p, T q, T r, T s)
{
	return max(max(p, q), max(r, s));
}

/**
 * @brief Clamp a value between @c mn and @c mx.
 *
 * For floats, NaNs are turned into @c mn.
 *
 * @tparam T    The value type.
 *
 * @param  v    The value to clamp.
 * @param  mn   The min value (inclusive).
 * @param  mx   The max value (inclusive).
 *
 * @return The clamped value.
 */
template<typename T>
inline T clamp(T v, T mn, T mx)
{
	// Preserve comparison order so NaN falls through to the minimum
	if (v > mx) return mx;
	if (v > mn) return v;
	return mn;
}

/**
 * @brief Clamp a float value between 0.0f and 1.0f.
 *
 * NaNs are turned into 0.0f.
 *
 * @param v   The value to clamp.
 *
 * @return The clamped value.
 */
static inline float clamp1f(float v)
{
	return astc::clamp(v, 0.0f, 1.0f);
}

/**
 * @brief Clamp a float value between 0.0f and 255.0f.
 *
 * NaNs are turned into 0.0f.
 *
 * @param v   The value to clamp.
 *
 * @return The clamped value.
 */
static inline float clamp255f(float v)
{
	return astc::clamp(v, 0.0f, 255.0f);
}

/**
 * @brief SP float round-down.
 *
 * @param v   The value to round.
 *
 * @return The rounded value.
 */
static inline float flt_rd(float v)
{
	return std::floor(v);
}

/**
 * @brief SP float round-to-nearest and convert to integer.
 *
 * @param v   The nonnegative value to round; halfway cases round upward.
 *
 * @return The rounded value.
 */
static inline int flt2int_rtn(float v)
{

	return static_cast<int>(v + 0.5f);
}

/**
 * @brief SP float round down and convert to integer.
 *
 * @param v   The value to round.
 *
 * @return The rounded value.
 */
static inline int flt2int_rd(float v)
{
	return static_cast<int>(v);
}

/**
 * @brief Compute the reciprocal square root.
 *
 * @param v   The input value.
 *
 * @return The reciprocal square root.
 */
static inline float rsqrt(float v)
{
	return 1.0f / std::sqrt(v);
}

/**
 * @brief Compute the square root.
 *
 * @param v   The input value.
 *
 * @return The square root.
 */
static inline float sqrt(float v)
{
	return std::sqrt(v);
}

/**
 * @brief Extract mantissa and exponent of a float value.
 *
 * @param      v      The finite, normalized input value.
 * @param[out] expo   The output exponent.
 *
 * @return The mantissa.
 */
static inline float frexp(float v, int* expo)
{
	unsigned int iv = std::bit_cast<unsigned int>(v);
	*expo = ((iv >> 23) & 0xFF) - 126;
	iv = (iv & 0x807fffff) | 0x3f000000;
	return std::bit_cast<float>(iv);
}

/**
 * @brief Compute the product of two sizes.
 *
 * This function is implemented to indicate if overflow has occurred, which may
 * occur when input values are not trusted. Implementation is obviously slower
 * than one that does not do this, so don't use for values we know cannot
 * overflow.
 *
 * Overflow signaling is sticky, so calling code can check at the end of a
 * sequence of multiplies.
 *
 * @param         val_a      The first value to multiply.
 * @param         val_b      The second value to multiply.
 * @param[in,out] overflow   Did previous or this calculation overflow?
 *
 * @return The multiplication result, which may have overflowed.
 */
static inline size_t mul_safe(
	size_t val_a,
	size_t val_b,
	bool& overflow
) {
	size_t result = val_a * val_b;
	overflow = overflow || ((val_b != 0) && ((result / val_b) != val_a));
	return result;
}

/**
 * @brief Get the number of blocks along a single axis.
 *
 * This function is implemented so that intermediate values will not overflow,
 * which may occur when input values are not trusted. Implementation is
 * obviously slower than one that does not do this, so don't use for values
 * we know cannot overflow.
 *
 * @param dim_axis    The axis dimension, in pixels.
 * @param dim_block   The block dimension, in pixels.
 *
 * @return The number of blocks needed in this dimension.
 */
static inline size_t get_block_count_safe(
	size_t dim_axis,
	size_t dim_block
) {
	// Compute number of whole blocks
	size_t blocks = dim_axis / dim_block;

	// Add in any residual partial block
	if (dim_axis != (dim_block * blocks))
	{
		blocks++;
	}

	return blocks;
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
 * @param[out] state   The state structure to initialize.
 */
void rand_init(uint64_t state[2]);

/**
 * @brief Return the next random number from the generator.
 *
 * This RNG is an implementation of the "xoroshiro-128+ 1.0" PRNG, based on the
 * public-domain implementation given by David Blackman & Sebastiano Vigna at
 * http://vigna.di.unimi.it/xorshift/xoroshiro128plus.c
 *
 * @param[in,out] state   The state structure to use/update.
 */
uint64_t rand(uint64_t state[2]);

}

/* ============================================================================
  Softfloat library with fp32 and fp16 conversion functionality.
============================================================================ */
#if (ASTCENC_F16C == 0) && (ASTCENC_NEON == 0)
	/* narrowing float->float conversions */
	uint16_t float_to_sf16(float val);
	float sf16_to_float(uint16_t val);
#endif

/*********************************
  Vector library
*********************************/
#include "astcenc_vecmathlib.h"

/*********************************
  Declaration of line types
*********************************/
/** @brief A 2D parametric line, @c a + @c b * @c t. */
struct line2
{
	/** @brief A point on the line. */
	vfloat4 a;

	/** @brief The line direction vector. */
	vfloat4 b;
};

/** @brief A 3D parametric line, @c a + @c b * @c t. */
struct line3
{
	/** @brief A point on the line. */
	vfloat4 a;

	/** @brief The line direction vector. */
	vfloat4 b;
};

/** @brief A 4D parametric line, @c a + @c b * @c t. */
struct line4
{
	/** @brief A point on the line. */
	vfloat4 a;

	/** @brief The line direction vector. */
	vfloat4 b;
};

/** @brief A 2D line preprocessed for projecting points onto it. */
struct processed_line2
{
	/** @brief The point on the line closest to the origin. */
	vfloat4 amod;

	/** @brief The unit-length line direction vector. */
	vfloat4 bs;
};

/** @brief A 3D line preprocessed for projecting points onto it. */
struct processed_line3
{
	/** @brief The point on the line closest to the origin. */
	vfloat4 amod;

	/** @brief The unit-length line direction vector. */
	vfloat4 bs;
};

/** @brief A 4D line preprocessed for projecting points onto it. */
struct processed_line4
{
	/** @brief The point on the line closest to the origin. */
	vfloat4 amod;

	/** @brief The unit-length line direction vector. */
	vfloat4 bs;
};

#endif
