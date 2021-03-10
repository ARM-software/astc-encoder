// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2021 Arm Limited
//
// Licensed under the Apache License, Version 2.0 (the "License"); you may not
// use this file except in compliance with the License. You may obtain a copy
// of the License at:
//
//	 http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
// WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
// License for the specific language governing permissions and limitations
// under the License.
// ----------------------------------------------------------------------------

/**
 * @brief Intrinsics for Armv7 NEON.
 *
 * This module implements a few Armv7-compatible intrinsics indentical to Armv8
 * ones. Thus, astcenc can be compiled using Armv7 architecture.
 */

#ifndef ASTC_VECMATHLIB_NEON_ARMV7_4_H_INCLUDED
#define ASTC_VECMATHLIB_NEON_ARMV7_4_H_INCLUDED

#ifndef ASTCENC_SIMD_INLINE
	#error "Include astcenc_vecmathlib.h, do not include directly"
#endif

#include <algorithm>

/**
 * @brief Return the max vector of two vectors.
 *
 * If one vector element is numeric and the other is a quiet NaN,
 * the result placed in the vector is the numerical value.
 */
ASTCENC_SIMD_INLINE float32x4_t vmaxnmq_f32(float32x4_t a, float32x4_t b)
{
	uint32x4_t amask = vceqq_f32(a, a);
	uint32x4_t bmask = vceqq_f32(b, b);
	a = vbslq_f32(amask, a, b);
	b = vbslq_f32(bmask, b, a);
	return vmaxq_f32(a, b);
}

/**
 * @brief Return the min vector of two vectors.
 *
 * If one vector element is numeric and the other is a quiet NaN,
 * the result placed in the vector is the numerical value.
 */
ASTCENC_SIMD_INLINE float32x4_t vminnmq_f32(float32x4_t a, float32x4_t b)
{
	uint32x4_t amask = vceqq_f32(a, a);
	uint32x4_t bmask = vceqq_f32(b, b);
	a = vbslq_f32(amask, a, b);
	b = vbslq_f32(bmask, b, a);
	return vminq_f32(a, b);
}

/**
 * @brief Return the horizontal maximum of a vector.
 */
ASTCENC_SIMD_INLINE float32_t vmaxvq_f32(float32x4_t a)
{
	float32_t a0 = vgetq_lane_f32(a, 0);
	float32_t a1 = vgetq_lane_f32(a, 1);
	float32_t a2 = vgetq_lane_f32(a, 2);
	float32_t a3 = vgetq_lane_f32(a, 3);
	return std::max(std::max(a0, a1), std::max(a2, a3));
}

/**
 * @brief Return the horizontal maximum of a vector.
 */
ASTCENC_SIMD_INLINE float32_t vminvq_f32(float32x4_t a)
{
	float32_t a0 = vgetq_lane_f32(a, 0);
	float32_t a1 = vgetq_lane_f32(a, 1);
	float32_t a2 = vgetq_lane_f32(a, 2);
	float32_t a3 = vgetq_lane_f32(a, 3);
	return std::min(std::min(a0, a1), std::min(a2, a3));
}

/**
 * @brief Return the horizontal maximum of a vector.
 */
ASTCENC_SIMD_INLINE int32_t vmaxvq_s32(int32x4_t a)
{
	int32_t a0 = vgetq_lane_s32(a, 0);
	int32_t a1 = vgetq_lane_s32(a, 1);
	int32_t a2 = vgetq_lane_s32(a, 2);
	int32_t a3 = vgetq_lane_s32(a, 3);
	return std::max(std::max(a0, a1), std::max(a2, a3));
}

/**
 * @brief Return the horizontal maximum of a vector.
 */
ASTCENC_SIMD_INLINE int32_t vminvq_s32(int32x4_t a)
{
	int32_t a0 = vgetq_lane_s32(a, 0);
	int32_t a1 = vgetq_lane_s32(a, 1);
	int32_t a2 = vgetq_lane_s32(a, 2);
	int32_t a3 = vgetq_lane_s32(a, 3);
	return std::min(std::min(a0, a1), std::min(a2, a3));
}

/**
 * @brief Return the sqrt of the lanes in the vector.
 */
ASTCENC_SIMD_INLINE float32x4_t vsqrtq_f32(float32x4_t a)
{
	float32_t a0 = std::sqrt(vgetq_lane_f32(a, 0));
	float32_t a1 = std::sqrt(vgetq_lane_f32(a, 1));
	float32_t a2 = std::sqrt(vgetq_lane_f32(a, 2));
	float32_t a3 = std::sqrt(vgetq_lane_f32(a, 3));
	return { a0, a1, a2, a3 };
}

/**
 * @brief Return a float rounded to the nearest integer value.
 */
ASTCENC_SIMD_INLINE float32x4_t vrndnq_f32(float32x4_t a)
{
	assert(::fegetround() == FE_TONEAREST);
	float32_t a0 = std::nearbyintf(vgetq_lane_f32(a, 0));
	float32_t a1 = std::nearbyintf(vgetq_lane_f32(a, 1));
	float32_t a2 = std::nearbyintf(vgetq_lane_f32(a, 2));
	float32_t a3 = std::nearbyintf(vgetq_lane_f32(a, 3));
	return { a0, a1, a2, a3 };
}

/**
 * @brief vector by vector division.
 */
ASTCENC_SIMD_INLINE float32x4_t vdivq_f32(float32x4_t a, float32x4_t b)
{
	float32x4_t x = vrecpeq_f32(b);
	x = vmulq_f32(vrecpsq_f32(b, x), x);
	x = vmulq_f32(vrecpsq_f32(b, x), x);
	return vmulq_f32(a, x);
}

/**
 * @brief Table vector Lookup.
 */
ASTCENC_SIMD_INLINE int8x16_t vqtbl1q_s8(int8x16_t t, uint8x16_t idx)
{
	int8x8x2_t tab;
	tab.val[0] = vget_low_s8(t);
	tab.val[1] = vget_high_s8(t);
	int8x16_t id = vreinterpretq_s8_u8(idx);
	return vcombine_s8(
		vtbl2_s8(tab, vget_low_s8(id)),
		vtbl2_s8(tab, vget_high_s8(id)));
}

#endif