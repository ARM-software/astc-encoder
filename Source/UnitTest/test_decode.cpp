// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2023-2026 Arm Limited
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

#include "../astcenc.h"

namespace astcenc
{

/** @brief Input overflow tests for xy * z. */
TEST(decompress, overflow_in_z)
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
	uint8_t output[1];

	uint8_t* slices = output;

	// x * y always fits in 64-bit size_t, but xy * z will overflow
	astcenc_image image {
		.dim_x = 0xFFFFFFFFu,
		.dim_y = 0xFFFFFFFFu,
		.dim_z = 0xFFFFFFFFu,
		.data_type = ASTCENC_TYPE_U8,
		.data = reinterpret_cast<void**>(&slices)
	};

	status = astcenc_decompress_image(context, data, -1, &image, &swizzle, 0);
	EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);

	astcenc_context_free(context);
}


/** @brief Compressed byte size overflow with a valid uncompressed byte size. */
TEST(decompress, overflow_in_16)
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
	uint8_t output[1];

	uint8_t* slices = output;

	// Rounding up to a full block makes the compressed byte size overflow
	astcenc_image image {
		.dim_x = std::numeric_limits<size_t>::max() / 4,
		.dim_y = 1,
		.dim_z = 1,
		.data_type = ASTCENC_TYPE_U8,
		.data = reinterpret_cast<void**>(&slices)
	};

	status = astcenc_decompress_image(context, data, -1, &image, &swizzle, 0);
	EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);

	astcenc_context_free(context);
}

/** @brief Image storage overflow must be rejected before accessing the data. */
TEST(decompress, image_size_overflow)
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
				.data = nullptr
			};

			status = astcenc_decompress_image(context, nullptr, 0, &image, &swizzle, 0);
			EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);

			// Include a multi-row or multi-slice product overflow
			dims[axis] = std::numeric_limits<size_t>::max() / (8 * component_size) + 1;
			dims[(axis + 1) % 3] = 2;
			image = {
				.dim_x = dims[0],
				.dim_y = dims[1],
				.dim_z = dims[2],
				.data_type = type,
				.data = nullptr
			};

			status = astcenc_decompress_image(context, nullptr, 0, &image, &swizzle, 0);
			EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);
		}
	}

	// The slice pointer array can overflow even when U8 texel storage fits
	astcenc_image image {
		.dim_x = 1,
		.dim_y = 1,
		.dim_z = std::numeric_limits<size_t>::max() / sizeof(void*) + 1,
		.data_type = ASTCENC_TYPE_U8,
		.data = nullptr
	};

	status = astcenc_decompress_image(context, nullptr, 0, &image, &swizzle, 0);
	EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);

	astcenc_context_free(context);
}

/** @brief Dimensions wider than unsigned int must reach the buffer capacity check. */
TEST(decompress, large_dimensions)
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
			.data = nullptr
		};

		status = astcenc_decompress_image(context, nullptr, 0, &image, &swizzle, 0);
		EXPECT_EQ(status, ASTCENC_ERR_OUT_OF_MEM);
	}

	astcenc_image image {
		.dim_x = 1,
		.dim_y = 1,
		.dim_z = 1,
		.data_type = ASTCENC_TYPE_U8,
		.data = nullptr
	};
	status = astcenc_decompress_image(
		context, nullptr, 0, &image, &swizzle, std::numeric_limits<size_t>::max());
	EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);

	astcenc_context_free(context);
}

/** @brief Input data buffer overrun tests. */
TEST(decompress, data_buffer_exceeded)
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
	uint8_t output[1];

	uint8_t* slices = output;
	astcenc_image image {
		.dim_x = 0x4u,
		.dim_y = 0x4u,
		.dim_z = 0x1u,
		.data_type = ASTCENC_TYPE_U8,
		.data = reinterpret_cast<void**>(&slices)
	};

	// Data size is 1 byte too short so this should error
	status = astcenc_decompress_image(context, data, 15, &image, &swizzle, 0);
	EXPECT_EQ(status, ASTCENC_ERR_OUT_OF_MEM);

	astcenc_context_free(context);
}

/** @brief Test harness for exploring issue #447. */
TEST(decompress, decompress12x12)
{
	astcenc_error status;
	astcenc_config config;
	astcenc_context* context;

	static const astcenc_swizzle swizzle {
		ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
	};

	uint8_t data[16] {
		0x29,0x00,0x1A,0x97,0x01,0x00,0x00,0x00,
		0x00,0x00,0x00,0x00,0x00,0xCF,0x97,0x86
	};

	uint8_t output[12*12*4];
	astcenc_config_init(ASTCENC_PRF_LDR, 12, 12, 1, ASTCENC_PRE_MEDIUM, 0, &config);

	status = astcenc_context_alloc(&config, 1, &context, nullptr);
	EXPECT_EQ(status, ASTCENC_SUCCESS);

	uint8_t* slices = output;
	astcenc_image image {
		.dim_x = 12,
		.dim_y = 12,
		.dim_z = 1,
		.data_type = ASTCENC_TYPE_U8,
		.data = reinterpret_cast<void**>(&slices)
	};

	status = astcenc_decompress_image(context, data, 16, &image, &swizzle, 0);
	EXPECT_EQ(status, ASTCENC_SUCCESS);

	astcenc_context_free(context);
}

/** @brief Test harness for context inheritance. */
TEST(decompress, context_inherit)
{
	astcenc_error status;
	astcenc_config config;
	astcenc_context* parent_context;
	astcenc_context* context;
	astcenc_context* error_context;

	static const astcenc_swizzle swizzle {
		ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
	};

	uint8_t data[16] {
		0x29,0x00,0x1A,0x97,0x01,0x00,0x00,0x00,
		0x00,0x00,0x00,0x00,0x00,0xCF,0x97,0x86
	};

	uint8_t output[12*12*4];
	astcenc_config_init(ASTCENC_PRF_LDR, 12, 12, 1, ASTCENC_PRE_MEDIUM, 0, &config);

	status = astcenc_context_alloc(&config, 1, &parent_context, nullptr);
	EXPECT_EQ(status, ASTCENC_SUCCESS);

	status = astcenc_context_alloc(nullptr, 1, &context, parent_context);
	EXPECT_EQ(status, ASTCENC_SUCCESS);

	status = astcenc_context_alloc(&config, 1, &error_context, parent_context);
	EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);

	uint8_t* slices = output;
	astcenc_image image {
		.dim_x = 12,
		.dim_y = 12,
		.dim_z = 1,
		.data_type = ASTCENC_TYPE_U8,
		.data = reinterpret_cast<void**>(&slices)
	};

	status = astcenc_decompress_image(context, data, 16, &image, &swizzle, 0);
	EXPECT_EQ(status, ASTCENC_SUCCESS);

	astcenc_context_free(context);
	astcenc_context_free(parent_context);
}

/** @brief Decode an unsupported partitioning with a self-decompress-only context. */
TEST(decompress, self_decompress_unsupported_partition)
{
	astcenc_error status;
	astcenc_config config;
	astcenc_context* context;

	static const astcenc_swizzle swizzle {
		ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
	};

	// A 3-partition block. An 8x8 fastest self-decompress-only context keeps no 3-partition tables,
	// so this selects a partitioning that is not present.
	uint8_t data[16] {
		0xCF, 0x50, 0x80, 0x64, 0x84, 0xCD, 0xD8, 0xB7,
		0xB0, 0xF1, 0x1E, 0x6A, 0x4F, 0x96, 0xC6, 0xBC
	};

	astcenc_config_init(ASTCENC_PRF_LDR, 8, 8, 1, ASTCENC_PRE_FASTEST,
	                    ASTCENC_FLG_DECOMPRESS_ONLY | ASTCENC_FLG_SELF_DECOMPRESS_ONLY,
	                    &config);

	status = astcenc_context_alloc(&config, 1, &context, nullptr);
	EXPECT_EQ(status, ASTCENC_SUCCESS);

	uint8_t output[8 * 8 * 4];
	uint8_t* slices = output;
	astcenc_image image {
		.dim_x = 8,
		.dim_y = 8,
		.dim_z = 1,
		.data_type = ASTCENC_TYPE_U8,
		.data = reinterpret_cast<void**>(&slices)
	};

	status = astcenc_decompress_image(context, data, 16, &image, &swizzle, 0);
	EXPECT_EQ(status, ASTCENC_SUCCESS);

	// An absent partitioning must decode as error color block (magenta for LDR)
	for (int i = 0; i < 8 * 8; i++)
	{
		EXPECT_EQ(output[4 * i + 0], 0xFF);
		EXPECT_EQ(output[4 * i + 1], 0x00);
		EXPECT_EQ(output[4 * i + 2], 0xFF);
		EXPECT_EQ(output[4 * i + 3], 0xFF);
	}

	astcenc_context_free(context);
}

}
