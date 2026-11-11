// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2025-2026 Arm Limited
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
 * @brief Unit tests for the vectorized SIMD functionality.
 */

#include <limits>

#include "gtest/gtest.h"

#include "../astcenc_internal.h"

namespace astcenc
{

/** @brief Preserve HDR mode bits when endpoint quantization rounds down across a bit boundary. */
TEST(compress, hdr_quantization_mode_bits)
{
	const vfloat4 colors[] {
		vfloat4(31806.0f, 27290.0f, 27570.0f, 65535.0f),
		vfloat4(30545.0f, 26180.0f, 28316.0f, 65535.0f)
	};
	const float offsets[] { 1209.0f, 1246.0f };
	const int formats[] { FMT_HDR_RGB_SCALE, FMT_HDR_RGB };
	const float error_limits[] { 5000000.0f, 300000.0f };

	for (size_t i = 0; i < 2; i++)
	{
		vfloat4 ep0 = colors[i];
		vfloat4 ep1 = ep0 + vfloat4(offsets[i], offsets[i], offsets[i], 0.0f);
		vfloat4 rgbo = ep0;
		rgbo.set_lane<3>(offsets[i]);

		uint8_t values[8];
		int format = pack_color_endpoints(ep0, ep1, vfloat4::zero(), rgbo,
		                                  formats[i], values, QUANT_16);
		EXPECT_EQ(format, formats[i]);

		bool is_rgb_hdr;
		bool is_alpha_hdr;
		vint4 decoded0;
		vint4 decoded1;
		unpack_color_endpoints(ASTCENC_PRF_HDR_RGB_LDR_A, format, values,
		                       is_rgb_hdr, is_alpha_hdr, decoded0, decoded1);
		EXPECT_TRUE(is_rgb_hdr);
		EXPECT_FALSE(is_alpha_hdr);

		vfloat4 diff0 = ep0 - int_to_float(decoded0);
		vfloat4 diff1 = ep1 - int_to_float(decoded1);
		float error = hadd_rgb_s(diff0 * diff0 + diff1 * diff1);
		EXPECT_LT(error, error_limits[i]);
	}
}

/** @brief Input overflow tests for xy * z. */
TEST(compress, overflow_in_z)
{
	astcenc_error status;
	astcenc_config config;
	astcenc_context* context;

	static const astcenc_swizzle swizzle {
		ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
	};

	astcenc_config_init(ASTCENC_PRF_LDR, 4, 4, 1, ASTCENC_PRE_MEDIUM, 0, &config);
	status = astcenc_context_alloc(&config, 1, &context, nullptr);
	EXPECT_EQ(status, ASTCENC_SUCCESS);

	// Arrays are too short, but should never be touched
	uint8_t data[1];
	uint8_t input[1];

	uint8_t* slices = input;

	// x * y always fits in 64-bit size_t, but xy * z will overflow
	astcenc_image image {
		.dim_x = 0xFFFFFFFFu,
		.dim_y = 0xFFFFFFFFu,
		.dim_z = 0xFFFFFFFFu,
		.data_type = ASTCENC_TYPE_U8,
		.data = reinterpret_cast<void**>(&slices),
		.row_stride = 0
	};

	status = astcenc_compress_image(context, &image, &swizzle, data, -1, 0);
	EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);

	astcenc_context_free(context);
}

/** @brief Input validation for config.block_z = 0. */
TEST(compress, zero_block_z)
{
	astcenc_error status;
	astcenc_config config;
	astcenc_context* context;

	static const astcenc_swizzle swizzle {
		ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
	};

	astcenc_config_init(ASTCENC_PRF_LDR, 4, 4, 1, ASTCENC_PRE_MEDIUM, 0, &config);
	config.block_z = 0;
	status = astcenc_context_alloc(&config, 1, &context, nullptr);
	EXPECT_EQ(status, ASTCENC_SUCCESS);

	uint8_t data[16];
	uint8_t input[4 * 4 * 4];

	uint8_t* slices = input;
	astcenc_image image {
		.dim_x = 4,
		.dim_y = 4,
		.dim_z = 1,
		.data_type = ASTCENC_TYPE_U8,
		.data = reinterpret_cast<void**>(&slices),
		.row_stride = 0
	};

	// Compression must succeed despite zero config.block_z value
	status = astcenc_compress_image(context, &image, &swizzle, data, sizeof(data), 0);
	EXPECT_EQ(status, ASTCENC_SUCCESS);

	astcenc_context_free(context);
}

/** @brief Compressed byte size overflow with a valid uncompressed byte size. */
TEST(compress, overflow_in_16)
{
	astcenc_error status;
	astcenc_config config;
	astcenc_context* context;

	static const astcenc_swizzle swizzle {
		ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
	};

	astcenc_config_init(ASTCENC_PRF_LDR, 4, 4, 1, ASTCENC_PRE_MEDIUM, 0, &config);
	status = astcenc_context_alloc(&config, 1, &context, nullptr);
	EXPECT_EQ(status, ASTCENC_SUCCESS);

	// Arrays are too short, but should never be touched
	uint8_t data[1];
	uint8_t input[1];

	uint8_t* slices = input;

	// Rounding up to a full block makes the compressed byte size overflow
	astcenc_image image {
		.dim_x = std::numeric_limits<size_t>::max() / 4,
		.dim_y = 1,
		.dim_z = 1,
		.data_type = ASTCENC_TYPE_U8,
		.data = reinterpret_cast<void**>(&slices),
		.row_stride = 0
	};

	status = astcenc_compress_image(context, &image, &swizzle, data, -1, 0);
	EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);

	astcenc_context_free(context);
}

/** @brief Input test for an out-of-range alpha-scale radius. */
TEST(compress, ascale_radius_clamped)
{
	astcenc_error status;
	astcenc_config config;
	astcenc_context* context;

	static const astcenc_swizzle swizzle {
		ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
	};

	astcenc_config_init(ASTCENC_PRF_LDR, 4, 4, 1, ASTCENC_PRE_MEDIUM, 0, &config);

	// Set to largest possible value to check input validation
	config.a_scale_radius = 0xFFFFFFFFu;
	status = astcenc_context_alloc(&config, 1, &context, nullptr);
	EXPECT_EQ(status, ASTCENC_SUCCESS);

	uint8_t data[16];
	uint8_t input[4 * 4 * 4];

	uint8_t* slices = input;
	astcenc_image image {
		.dim_x = 4,
		.dim_y = 4,
		.dim_z = 1,
		.data_type = ASTCENC_TYPE_U8,
		.data = reinterpret_cast<void**>(&slices),
		.row_stride = 0
	};

	// Compression should succeed without bad_alloc or ASAN failure
	status = astcenc_compress_image(context, &image, &swizzle, data, sizeof(data), 0);
	EXPECT_EQ(status, ASTCENC_SUCCESS);

	astcenc_context_free(context);
}

/** @brief Image storage overflow must be rejected before accessing the data. */
TEST(compress, image_size_overflow)
{
	astcenc_error status;
	astcenc_config config;
	astcenc_context* context;

	// Use 3D blocks so compressed storage does not mask RGBA byte overflow
	astcenc_config_init(ASTCENC_PRF_LDR, 6, 6, 6, ASTCENC_PRE_MEDIUM, 0, &config);
	status = astcenc_context_alloc(&config, 1, &context, nullptr);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	const astcenc_swizzle swizzle {
		ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
	};

	for (size_t component_size : { 1u, 2u, 4u })
	{
		astcenc_type type = component_size == 1 ? ASTCENC_TYPE_U8 :
		                    component_size == 2 ? ASTCENC_TYPE_F16 : ASTCENC_TYPE_F32;

		for (size_t axis = 0; axis < 3; axis++)
		{
			size_t dims[3] { 1, 1, 1 };
			dims[axis] = std::numeric_limits<size_t>::max() / (4 * component_size) + 1;
			astcenc_image image {
				.dim_x = dims[0],
				.dim_y = dims[1],
				.dim_z = dims[2],
				.data_type = type,
				.data = nullptr,
				.row_stride = 0
			};
			status = astcenc_compress_image(context, &image, &swizzle, nullptr, 0, 0);
			EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);

			// Include a multi-row or multi-slice product overflow
			dims[axis] = std::numeric_limits<size_t>::max() / (8 * component_size) + 1;
			dims[(axis + 1) % 3] = 2;
			image = {
				.dim_x = dims[0],
				.dim_y = dims[1],
				.dim_z = dims[2],
				.data_type = type,
				.data = nullptr,
				.row_stride = 0
			};

			status = astcenc_compress_image(context, &image, &swizzle, nullptr, 0, 0);
			EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);
		}
	}

	// The slice pointer array can overflow even when U8 texel storage fits
	astcenc_image image {
		.dim_x = 1,
		.dim_y = 1,
		.dim_z = std::numeric_limits<size_t>::max() / sizeof(void*) + 1,
		.data_type = ASTCENC_TYPE_U8,
		.data = nullptr,
		.row_stride = 0
	};
	status = astcenc_compress_image(context, &image, &swizzle, nullptr, 0, 0);
	EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);

	astcenc_context_free(context);
}

/** @brief Dimensions wider than unsigned int must reach the buffer capacity check. */
TEST(compress, large_dimensions)
{
	if (sizeof(size_t) <= sizeof(unsigned int))
	{
		GTEST_SKIP();
	}

	astcenc_error status;
	astcenc_config config;
	astcenc_context* context;
	astcenc_config_init(ASTCENC_PRF_LDR, 4, 4, 1, ASTCENC_PRE_MEDIUM, 0, &config);
	status = astcenc_context_alloc(&config, 1, &context, nullptr);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	const astcenc_swizzle swizzle {
		ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
	};

	for (size_t axis = 0; axis < 3; axis++)
	{
		size_t dims[3] { 1, 1, 1 };
		dims[axis] = static_cast<size_t>(std::numeric_limits<unsigned int>::max()) + 1;
		astcenc_image image {
			.dim_x = dims[0],
			.dim_y = dims[1],
			.dim_z = dims[2],
			.data_type = ASTCENC_TYPE_U8,
			.data = nullptr,
			.row_stride = 0
		};

		status = astcenc_compress_image(context, &image, &swizzle, nullptr, 0, 0);
		EXPECT_EQ(status, ASTCENC_ERR_OUT_OF_MEM);
	}

	astcenc_image image {
		.dim_x = 1,
		.dim_y = 1,
		.dim_z = 1,
		.data_type = ASTCENC_TYPE_U8,
		.data = nullptr,
		.row_stride = 0
	};

	status = astcenc_compress_image(
		context, &image, &swizzle, nullptr, 0, std::numeric_limits<size_t>::max());
	EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);

	astcenc_context_free(context);
}

/** @brief Per-thread scratch storage must not wrap when allocating a context. */
TEST(compress, thread_count_overflow)
{
	astcenc_error status;
	astcenc_config config;
	astcenc_context* context = nullptr;
	astcenc_config_init(ASTCENC_PRF_LDR, 4, 4, 1, ASTCENC_PRE_MEDIUM, 0, &config);
	status = astcenc_context_alloc(&config, std::numeric_limits<size_t>::max(), &context, nullptr);
	EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);
	EXPECT_EQ(context, nullptr);
}

/** @brief Input data buffer overrun tests. */
TEST(compress, data_buffer_exceeded)
{
	astcenc_error status;
	astcenc_config config;
	astcenc_context* context;

	static const astcenc_swizzle swizzle {
		ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
	};

	astcenc_config_init(ASTCENC_PRF_LDR, 4, 4, 1, ASTCENC_PRE_MEDIUM, 0, &config);
	status = astcenc_context_alloc(&config, 1, &context, nullptr);
	EXPECT_EQ(status, ASTCENC_SUCCESS);

	// Arrays are too short, but should never be touched
	uint8_t data[1];
	uint8_t input[1];

	uint8_t* slices = input;
	astcenc_image image {
		.dim_x = 0x4u,
		.dim_y = 0x4u,
		.dim_z = 0x1u,
		.data_type = ASTCENC_TYPE_U8,
		.data = reinterpret_cast<void**>(&slices),
		.row_stride = 0
	};

	// Data size is 1 byte too short so this should error
	status = astcenc_compress_image(context, &image, &swizzle, data, 15, 0);
	EXPECT_EQ(status, ASTCENC_ERR_OUT_OF_MEM);

	astcenc_context_free(context);
}

/**
 * @brief Helper to run a compression on a single 4x4 block.
 *
 * @param profile   The color profile to test.
 * @param input     The 4x4 block of interleaved RGBA float values.
 */
static void run_positive_test(
	astcenc_profile profile,
	float* input
) {
	astcenc_error status;
	astcenc_config config;
	astcenc_context* context;

	static const astcenc_swizzle swizzle {
		ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
	};

	astcenc_config_init(profile, 4, 4, 1, ASTCENC_PRE_MEDIUM, 0, &config);
	status = astcenc_context_alloc(&config, 1, &context, nullptr);
	EXPECT_EQ(status, ASTCENC_SUCCESS);

	// Single 16 byte block output data
	uint8_t data[16];

	float* slices = input;
	astcenc_image image {
		.dim_x = 4u,
		.dim_y = 4u,
		.dim_z = 1u,
		.data_type = ASTCENC_TYPE_F32,
		.data = reinterpret_cast<void**>(&slices),
		.row_stride = 0
	};

	// Should not error or trip any UBSAN errors
	status = astcenc_compress_image(context, &image, &swizzle, data, 16, 0);
	EXPECT_EQ(status, ASTCENC_SUCCESS);

	astcenc_context_free(context);
}

/** @brief Input data contains negative infinities. */
TEST(compress, data_neg_inf_ldr)
{
	// Cycle through a single bad value in each color channel
	for(int offset = 0; offset < 4; offset++)
	{
		// 4x4 block input data
		std::vector<float> input(4 * 4 * 4, 0.5f);
		input[0 + offset] = -std::numeric_limits<float>::infinity();

		// Make data more complex to hit more compressor paths
		input[4 + offset] = 0.75f;
		input[8 + offset] = 0.35f;
		input[12 + offset] = 0.0f;

		run_positive_test(ASTCENC_PRF_LDR, input.data());
	}
}

/** @brief Input data contains negative infinities. */
TEST(compress, data_neg_inf_hdr)
{
	// Cycle through a single bad value in each color channel
	for(int offset = 0; offset < 4; offset++)
	{
		// 4x4 block input data
		std::vector<float> input(4 * 4 * 4, 0.5f);
		input[0 + offset] = -std::numeric_limits<float>::infinity();

		run_positive_test(ASTCENC_PRF_HDR, input.data());
	}
}

/** @brief Input data contains negative infinities. */
TEST(compress, data_neg_inf_hdr_ldra)
{
	// Cycle through a single bad value in each color channel
	for(int offset = 0; offset < 4; offset++)
	{
		// 4x4 block input data
		std::vector<float> input(4 * 4 * 4, 0.5f);
		input[0 + offset] = -std::numeric_limits<float>::infinity();

		run_positive_test(ASTCENC_PRF_HDR_RGB_LDR_A, input.data());
	}
}

/** @brief Input data contains positive infinities. */
TEST(compress, data_pos_inf_ldr)
{
	// Cycle through a single bad value in each color channel
	for(int offset = 0; offset < 4; offset++)
	{
		// 4x4 block input data
		std::vector<float> input(4 * 4 * 4, 0.5f);
		input[0 + offset] = std::numeric_limits<float>::infinity();

		run_positive_test(ASTCENC_PRF_LDR, input.data());
	}
}

/** @brief Input data contains positive infinities. */
TEST(compress, data_pos_inf_hdr)
{
	// Cycle through a single bad value in each color channel
	for(int offset = 0; offset < 4; offset++)
	{
		// 4x4 block input data
		std::vector<float> input(4 * 4 * 4, 0.5f);
		input[0 + offset] = std::numeric_limits<float>::infinity();

		run_positive_test(ASTCENC_PRF_HDR, input.data());
	}
}

/** @brief Input data contains positive infinities. */
TEST(compress, data_pos_inf_hdr_ldra)
{
	// Cycle through a single bad value in each color channel
	for(int offset = 0; offset < 4; offset++)
	{
		// 4x4 block input data
		std::vector<float> input(4 * 4 * 4, 0.5f);
		input[0 + offset] = std::numeric_limits<float>::infinity();

		run_positive_test(ASTCENC_PRF_HDR_RGB_LDR_A, input.data());
	}
}


/** @brief Input data contains NaN. */
TEST(compress, data_nan_ldr)
{
	// Cycle through a single bad value in each color channel
	for(int offset = 0; offset < 4; offset++)
	{
		// 4x4 block input data
		std::vector<float> input(4 * 4 * 4, 0.5f);
		input[0 + offset] = std::numeric_limits<float>::quiet_NaN();

		run_positive_test(ASTCENC_PRF_LDR, input.data());
	}
}

/** @brief Input data contains NaN. */
TEST(compress, data_nan_hdr)
{
	// Cycle through a single bad value in each color channel
	for(int offset = 0; offset < 4; offset++)
	{
		// 4x4 block input data
		std::vector<float> input(4 * 4 * 4, 0.5f);
		input[0 + offset] = std::numeric_limits<float>::quiet_NaN();

		run_positive_test(ASTCENC_PRF_HDR, input.data());
	}
}

/** @brief Input data contains NaN. */
TEST(compress, data_nan_hdr_ldra)
{
	// Cycle through a single bad value in each color channel
	for(int offset = 0; offset < 4; offset++)
	{
		// 4x4 block input data
		std::vector<float> input(4 * 4 * 4, 0.5f);
		input[0 + offset] = std::numeric_limits<float>::quiet_NaN();

		run_positive_test(ASTCENC_PRF_HDR_RGB_LDR_A, input.data());
	}
}

}
