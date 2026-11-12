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

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <thread>

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
		.data = reinterpret_cast<void**>(&slices),
		.row_stride = 0
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
		.data = reinterpret_cast<void**>(&slices),
		.row_stride = 0
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
				.data = nullptr,
				.row_stride = 0
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
				.data = nullptr,
				.row_stride = 0
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
		.data = nullptr,
		.row_stride = 0
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
			.data = nullptr,
			.row_stride = 0
		};

		status = astcenc_decompress_image(context, nullptr, 0, &image, &swizzle, 0);
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
		.data = reinterpret_cast<void**>(&slices),
		.row_stride = 0
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
		.data = reinterpret_cast<void**>(&slices),
		.row_stride = 0
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
		.data = reinterpret_cast<void**>(&slices),
		.row_stride = 0
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
		.data = reinterpret_cast<void**>(&slices),
		.row_stride = 0
	};

	status = astcenc_decompress_image(context, data, 16, &image, &swizzle, 0);
	EXPECT_EQ(status, ASTCENC_ERR_BAD_BLOCK);

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

/** @brief Report malformed blocks without confusing valid magenta content with error colors. */
TEST(decompress, invalid_blocks)
{
	const uint8_t invalid_data[][16] {
		{ 0 },
		{ 0xFC, 0x01, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF },
		{ 0xFC, 0x0D, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 }
	};

	const uint8_t magenta_data[16] {
		0xFC, 0xFD, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0xFF, 0xFF, 0x00, 0x00, 0xFF, 0xFF, 0xFF, 0xFF
	};

	const astcenc_swizzle swizzle {
		ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
	};

	const astcenc_profile profiles[] {
		ASTCENC_PRF_LDR, ASTCENC_PRF_LDR_SRGB, ASTCENC_PRF_HDR_RGB_LDR_A, ASTCENC_PRF_HDR
	};

	const astcenc_type types[] { ASTCENC_TYPE_U8, ASTCENC_TYPE_F16, ASTCENC_TYPE_F32 };

	for (auto profile : profiles)
	{
		astcenc_config config;
		astcenc_error status = astcenc_config_init(profile, 4, 4, 1, ASTCENC_PRE_FASTEST,
		                                           ASTCENC_FLG_DECOMPRESS_ONLY, &config);
		ASSERT_EQ(status, ASTCENC_SUCCESS);

		astcenc_context* context;
		status = astcenc_context_alloc(&config, 1, &context, nullptr);
		ASSERT_EQ(status, ASTCENC_SUCCESS);

		uint8_t output_u8[4 * 4 * 4];
		uint16_t output_f16[4 * 4 * 4];
		float output_f32[4 * 4 * 4];
		void* outputs[] { output_u8, output_f16, output_f32 };
		for (size_t type = 0; type < 3; type++)
		{
			astcenc_image image {
				.dim_x = 4,
				.dim_y = 4,
				.dim_z = 1,
				.data_type = types[type],
				.data = &outputs[type],
				.row_stride = 0
			};

			for (const auto& data : invalid_data)
			{
				status = astcenc_decompress_image(context, data, sizeof(data), &image, &swizzle, 0);
				EXPECT_EQ(status, ASTCENC_ERR_BAD_BLOCK);
				for (size_t i = 0; i < 4 * 4 * 4; i++)
				{
					if (types[type] == ASTCENC_TYPE_U8)
					{
						uint8_t expected = 0xFF;
						if ((i % 4) == 1)
						{
							expected = 0;
						}

						EXPECT_EQ(output_u8[i], expected);
					}
					else if (types[type] == ASTCENC_TYPE_F16)
					{
						EXPECT_EQ(output_f16[i], 0xFFFF);
					}
					else
					{
						bool is_nan = std::isnan(output_f32[i]);
						EXPECT_TRUE(is_nan);
					}
				}

				// Single-threaded contexts reset implicitly after a failed image
				status = astcenc_decompress_image(context, magenta_data, sizeof(magenta_data),
				                                  &image, &swizzle, 0);
				EXPECT_EQ(status, ASTCENC_SUCCESS);
			}
		}

		astcenc_context_free(context);
	}

	const char* error_string = astcenc_get_error_string(ASTCENC_ERR_BAD_BLOCK);
	EXPECT_STREQ(error_string, "ASTCENC_ERR_BAD_BLOCK");
}

/** @brief Report HDR constant and endpoint encodings when decoded with an LDR profile. */
TEST(decompress, unsupported_hdr_encodings)
{
	const uint8_t data[][16] {
		{ 0xFC, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		  0x80, 0x3E, 0x60, 0x3F, 0x40, 0x40, 0x00, 0x3C },
		{ 0x53, 0xE2, 0xAC, 0x61, 0x81, 0x9B, 0x69, 0xF5,
		  0x34, 0xDA, 0x1D, 0x23, 0xBF, 0xF3, 0xC1, 0xB8 }
	};

	const astcenc_swizzle swizzle {
		ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
	};

	const astcenc_profile profiles[] {
		ASTCENC_PRF_LDR, ASTCENC_PRF_LDR_SRGB, ASTCENC_PRF_HDR_RGB_LDR_A, ASTCENC_PRF_HDR
	};

	for (auto profile : profiles)
	{
		astcenc_config config;
		astcenc_error status = astcenc_config_init(profile, 4, 4, 1, ASTCENC_PRE_FASTEST,
		                                           ASTCENC_FLG_DECOMPRESS_ONLY, &config);
		ASSERT_EQ(status, ASTCENC_SUCCESS);

		astcenc_context* context;
		status = astcenc_context_alloc(&config, 1, &context, nullptr);
		ASSERT_EQ(status, ASTCENC_SUCCESS);

		astcenc_error expected = ASTCENC_SUCCESS;
		if ((profile == ASTCENC_PRF_LDR) || (profile == ASTCENC_PRF_LDR_SRGB))
		{
			expected = ASTCENC_ERR_BAD_BLOCK;
		}

		float output[4 * 4 * 4];
		void* slice = output;
		astcenc_image image {
			.dim_x = 4,
			.dim_y = 4,
			.dim_z = 1,
			.data_type = ASTCENC_TYPE_F32,
			.data = &slice,
			.row_stride = 0
		};

		for (const auto& block : data)
		{
			status = astcenc_decompress_image(context, block, sizeof(block), &image, &swizzle, 0);
			EXPECT_EQ(status, expected);
		}

		astcenc_context_free(context);
	}
}

/** @brief Finish decoding mixed valid and invalid images across workers, and reset for reuse. */
TEST(decompress, invalid_blocks_threaded)
{
	const uint8_t black_data[16] {
		0xFC, 0xFD, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xFF, 0xFF
	};

	std::array<uint8_t, 257 * 16> data;
	for (size_t i = 0; i < 257; i++)
	{
		std::copy_n(black_data, 16, data.data() + i * 16);
	}

	const astcenc_swizzle swizzle {
		ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
	};

	astcenc_config config;
	astcenc_error status = astcenc_config_init(ASTCENC_PRF_LDR, 4, 4, 1, ASTCENC_PRE_FASTEST,
	                                           ASTCENC_FLG_DECOMPRESS_ONLY, &config);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	astcenc_context* context;
	status = astcenc_context_alloc(&config, 2, &context, nullptr);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	std::array<uint8_t, 257 * 4 * 4 * 4> output;
	void* slice = output.data();
	astcenc_image image {
		.dim_x = 257 * 4,
		.dim_y = 4,
		.dim_z = 1,
		.data_type = ASTCENC_TYPE_U8,
		.data = &slice,
		.row_stride = 0
	};

	for (bool is_invalid : { true, false })
	{
		if (is_invalid)
		{
			std::fill_n(data.data(), 16, 0);
		}
		else
		{
			std::copy_n(black_data, 16, data.data());
		}

		output.fill(0x55);
		astcenc_error results[2];
		std::thread worker([&] {
			results[1] = astcenc_decompress_image(context, data.data(), data.size(), &image, &swizzle, 1);
		});

		results[0] = astcenc_decompress_image(context, data.data(), data.size(), &image, &swizzle, 0);
		worker.join();

		size_t failures = 0;
		for (auto result : results)
		{
			if (result == ASTCENC_ERR_BAD_BLOCK)
			{
				failures++;
			}
			else
			{
				EXPECT_EQ(result, ASTCENC_SUCCESS);
			}
		}

		if (is_invalid)
		{
			EXPECT_GT(failures, 0u);
		}
		else
		{
			EXPECT_EQ(failures, 0u);
		}

		for (size_t i = 0; i < output.size(); i++)
		{
			size_t x = (i / 4) % image.dim_x;
			size_t channel = i % 4;
			uint8_t expected = 0;
			if ((channel == 3) || (is_invalid && (x < 4) && (channel != 1)))
			{
				expected = 0xFF;
			}

			EXPECT_EQ(output[i], expected);
		}

		status = astcenc_decompress_reset(context);
		EXPECT_EQ(status, ASTCENC_SUCCESS);
	}

	astcenc_context_free(context);
}

}
