// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2026 Arm Limited
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
 * @brief Unit tests for image storage and block traversal.
 */

#include <algorithm>
#include <limits>
#include <memory>
#include <thread>
#include <type_traits>
#include <vector>

#include "gtest/gtest.h"

#include "../astcenc.h"
#include "../astcenc_morton.h"

namespace astcenc
{

/** @brief An owning pointer to a codec context. */
using astcenc_context_unq_ptr = std::unique_ptr<astcenc_context, decltype(&astcenc_context_free)>;

/**
 * @brief Compare packed and strided compression and decompression.
 *
 * @tparam T                 The component storage type.
 *
 * @param  type              The component format.
 * @param  depth             The number of image planes.
 * @param  is_alpha_scaled   Whether to exercise swizzling and alpha preprocessing.
 */
template<typename T>
static void test_row_stride(
	astcenc_type type,
	size_t depth,
	bool is_alpha_scaled
) {
	constexpr size_t width = 5;
	constexpr size_t height = 7;
	constexpr size_t row_components = 4 * width + 3;
	constexpr size_t offset = row_components + 4;
	constexpr size_t slice_components = offset + height * row_components + 4;
	constexpr T sentinel = static_cast<T>(42);

	std::vector<T> packed(width * height * depth * 4);
	std::vector<T> strided(slice_components * depth, sentinel);
	std::vector<void*> packed_slices(depth);
	std::vector<void*> strided_slices(depth);
	for (size_t z = 0; z < depth; z++)
	{
		packed_slices[z] = packed.data() + z * width * height * 4;
		strided_slices[z] = strided.data() + z * slice_components + offset;
		for (size_t y = 0; y < height; y++)
		{
			for (size_t x = 0; x < 4 * width; x++)
			{
				size_t value = (x * 3 + y * 2 + z * 7) % 5;
				T component;
				if constexpr (std::is_same_v<T, uint16_t>)
				{
					constexpr uint16_t values[] { 0, 0x3400, 0x3800, 0x3A00, 0x3C00 };
					component = values[value];
				}
				else if constexpr (std::is_same_v<T, uint8_t>)
				{
					component = static_cast<T>(value * 63);
				}
				else
				{
					component = static_cast<T>(value) * 0.25f;
				}

				packed[(z * height + y) * width * 4 + x] = component;
				strided[z * slice_components + offset + y * row_components + x] = component;
			}
		}
	}

	astcenc_image packed_image {
		.dim_x = width,
		.dim_y = height,
		.dim_z = depth,
		.data_type = type,
		.data = packed_slices.data(),
		.row_stride = 0
	};

	astcenc_image strided_image {
		.dim_x = width,
		.dim_y = height,
		.dim_z = depth,
		.data_type = type,
		.data = strided_slices.data(),
		.row_stride = row_components * sizeof(T)
	};

	astcenc_config config;
	unsigned int block_xy = 3;
	unsigned int block_z = 3;
	if (depth == 1)
	{
		block_xy = 4;
		block_z = 1;
	}

	astcenc_error status = astcenc_config_init(
	    ASTCENC_PRF_LDR, block_xy, block_xy, block_z, ASTCENC_PRE_FASTEST, 0, &config);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	if (is_alpha_scaled)
	{
		config.a_scale_radius = 1;
	}

	astcenc_context* context_raw;
	status = astcenc_context_alloc(&config, 1, &context_raw, nullptr);
	ASSERT_EQ(status, ASTCENC_SUCCESS);
	astcenc_context_unq_ptr context(context_raw, astcenc_context_free);

	astcenc_swizzle swizzle {
		.r = ASTCENC_SWZ_R,
		.g = ASTCENC_SWZ_G,
		.b = ASTCENC_SWZ_B,
		.a = ASTCENC_SWZ_A
	};

	if (is_alpha_scaled)
	{
		swizzle.r = ASTCENC_SWZ_B;
		swizzle.b = ASTCENC_SWZ_R;
	}

	size_t compressed_size = ((width + block_xy - 1) / block_xy)
	                       * ((height + block_xy - 1) / block_xy) * 16;
	std::vector<uint8_t> packed_compressed(compressed_size);
	std::vector<uint8_t> strided_compressed(compressed_size);

	status = astcenc_compress_image(context.get(), &packed_image, &swizzle,
	                                packed_compressed.data(), compressed_size, 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);
	EXPECT_EQ(packed_image.row_stride, 0u);

	status = astcenc_compress_reset(context.get());
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	status = astcenc_compress_image(context.get(), &strided_image, &swizzle,
	                                strided_compressed.data(), compressed_size, 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);
	EXPECT_EQ(strided_compressed, packed_compressed);

	// Explicit packed pitch must match the zero-stride default
	packed_image.row_stride = width * 4 * sizeof(T);
	status = astcenc_compress_reset(context.get());
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	status = astcenc_compress_image(context.get(), &packed_image, &swizzle,
	                                packed_compressed.data(), compressed_size, 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);
	EXPECT_EQ(strided_compressed, packed_compressed);

	std::fill(strided.begin(), strided.end(), sentinel);
	packed_image.row_stride = 0;

	status = astcenc_decompress_image(context.get(), packed_compressed.data(), compressed_size,
	                                  &packed_image, &swizzle, 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);
	EXPECT_EQ(packed_image.row_stride, 0u);

	status = astcenc_decompress_reset(context.get());
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	status = astcenc_decompress_image(context.get(), packed_compressed.data(), compressed_size,
	                                  &strided_image, &swizzle, 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	std::vector<T> expected(strided.size(), sentinel);
	for (size_t z = 0; z < depth; z++)
	{
		for (size_t y = 0; y < height; y++)
		{
			std::copy_n(packed.data() + (z * height + y) * width * 4, width * 4,
			            expected.data() + z * slice_components + offset + y * row_components);
		}
	}
	EXPECT_EQ(strided, expected);
}

/** @brief Row stride supports every component type, cropped origins, and 3D images. */
TEST(image, row_stride)
{
	for (size_t depth : { 1u, 3u })
	{
		for (bool is_alpha_scaled : { false, true })
		{
			SCOPED_TRACE(::testing::Message() << "Depth " << depth << ", alpha " << is_alpha_scaled);
			test_row_stride<uint8_t>(ASTCENC_TYPE_U8, depth, is_alpha_scaled);
			test_row_stride<uint16_t>(ASTCENC_TYPE_F16, depth, is_alpha_scaled);
			test_row_stride<float>(ASTCENC_TYPE_F32, depth, is_alpha_scaled);
		}
	}
}

/**
 * @brief Check partial brick traversal preserves raster output without padding blocks.
 *
 * @param is_morton   Whether to use Morton order instead of scanline order.
 */
static void test_bricks(bool is_morton)
{
	const unsigned int block_sizes[][3] { { 4, 4, 1 }, { 10, 6, 1 }, { 3, 3, 3 } };
	const astcenc_swizzle swizzle {
		ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
	};
	const uint8_t constant_header[8] { 0xFC, 0xFD, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF };

	for (const auto& block : block_sizes)
	{
		size_t brick_width = 128;
		if (is_morton)
		{
			brick_width = 16;
		}

		size_t width = block[0] * brick_width + 1;
		size_t height = block[1] * 8 + 1;
		size_t depth = 1;
		if (block[2] > 1)
		{
			depth = block[2] * 2 + 1;
		}

		size_t blocks_x = (width + block[0] - 1) / block[0];
		size_t blocks_y = (height + block[1] - 1) / block[1];
		size_t blocks_z = (depth + block[2] - 1) / block[2];
		size_t block_count = blocks_x * blocks_y * blocks_z;
		size_t data_size = block_count * 16;

		// Each raster block has a distinct constant color, with guards outside the output payload
		std::vector<uint8_t> expected(data_size + 32, 0x55);
		for (size_t i = 0; i < block_count; i++)
		{
			uint8_t* output = expected.data() + 16 + i * 16;
			std::copy_n(constant_header, 8, output);
			for (size_t channel = 0; channel < 4; channel++)
			{
				uint8_t value = static_cast<uint8_t>(((i >> (channel * 3)) + channel * 47) & 0xFF);
				output[8 + channel * 2] = value;
				output[9 + channel * 2] = value;
			}
		}

		std::vector<uint8_t> pixels(width * height * depth * 4);
		std::vector<void*> slices(depth);
		for (size_t z = 0; z < depth; z++)
		{
			slices[z] = pixels.data() + z * width * height * 4;
			for (size_t y = 0; y < height; y++)
			{
				for (size_t x = 0; x < width; x++)
				{
					size_t index = ((z / block[2]) * blocks_y + y / block[1]) * blocks_x
					             + x / block[0];
					for (size_t channel = 0; channel < 4; channel++)
					{
						pixels[((z * height + y) * width + x) * 4 + channel] =
						    expected[16 + index * 16 + 8 + channel * 2];
					}
				}
			}
		}

		astcenc_image image {
			.dim_x = width,
			.dim_y = height,
			.dim_z = depth,
			.data_type = ASTCENC_TYPE_U8,
			.data = slices.data(),
			.row_stride = 0
		};

		unsigned int flags = 0;
		if (is_morton)
		{
			flags |= ASTCENC_FLG_USE_MORTON;
		}

		astcenc_config config;
		astcenc_error status = astcenc_config_init(
		    ASTCENC_PRF_LDR, block[0], block[1], block[2], ASTCENC_PRE_FASTEST, flags, &config);
		ASSERT_EQ(status, ASTCENC_SUCCESS);

		for (unsigned int thread_count : { 1u, 4u })
		{
			SCOPED_TRACE(::testing::Message() << block[0] << 'x' << block[1] << 'x' << block[2]
			                                << ", threads " << thread_count);
			astcenc_context* context_raw;
			status = astcenc_context_alloc(&config, thread_count, &context_raw, nullptr);
			ASSERT_EQ(status, ASTCENC_SUCCESS);
			astcenc_context_unq_ptr context(context_raw, astcenc_context_free);

			std::vector<uint8_t> output(data_size + 32, 0x55);
			std::vector<astcenc_error> results(thread_count);
			std::vector<std::thread> workers;
			for (unsigned int thread_index = 0; thread_index < thread_count; thread_index++)
			{
				workers.emplace_back([&, thread_index] {
					results[thread_index] = astcenc_compress_image(
					    context.get(), &image, &swizzle,
					    output.data() + 16, data_size, thread_index);
				});
			}

			for (auto& worker : workers)
			{
				worker.join();
			}

			for (astcenc_error result : results)
			{
				ASSERT_EQ(result, ASTCENC_SUCCESS);
			}

			EXPECT_EQ(output, expected);
		}
	}
}

/** @brief Morton bricks cover partial edges and preserve raster output. */
TEST(image, morton_bricks)
{
	test_bricks(true);
}

/** @brief Scanline bricks cover partial edges and preserve raster output. */
TEST(image, scanline_bricks)
{
	test_bricks(false);
}

/** @brief Both traversal orders cover their brick once, in the expected block order. */
TEST(image, brick_layout)
{
	for (bool is_morton : { false, true })
	{
		brick_layout layout(is_morton);
		unsigned int expected_width = 128;
		unsigned int expected_height = 1;
		if (is_morton)
		{
			expected_width = 16;
			expected_height = 8;
		}

		EXPECT_EQ(layout.width, expected_width);
		EXPECT_EQ(layout.height, expected_height);
		brick_layout::validity_mask visited;
		for (unsigned int index = 0; index < brick_layout::block_count; index++)
		{
			auto position = layout.get_position(index);
			ASSERT_LT(position.x, expected_width);
			ASSERT_LT(position.y, expected_height);
			EXPECT_EQ(position.z, 0u);
			size_t raster_index = position.y * expected_width + position.x;
			bool is_visited = visited[raster_index];
			EXPECT_FALSE(is_visited);
			visited.set(raster_index);

			size_t expected_index = position.x;
			if (is_morton)
			{
				size_t x = position.x;
				size_t y = position.y;
				expected_index = (x & 1) | ((y & 1) << 1) | ((x & 2) << 1) |
				                 ((y & 2) << 2) | ((x & 4) << 2) | ((y & 4) << 3) |
				                 ((x & 8) << 3);
			}

			EXPECT_EQ(index, expected_index);
		}

		bool is_complete = visited.all();
		EXPECT_TRUE(is_complete);
	}
}

/** @brief Invalid strides must be rejected before accessing image storage. */
TEST(image, invalid_row_stride)
{
	astcenc_config config;
	astcenc_error status = astcenc_config_init(
	    ASTCENC_PRF_LDR, 4, 4, 1, ASTCENC_PRE_FASTEST, 0, &config);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	astcenc_context* context_raw;
	status = astcenc_context_alloc(&config, 1, &context_raw, nullptr);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	astcenc_context_unq_ptr context(context_raw, astcenc_context_free);

	const astcenc_swizzle swizzle {
		ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
	};

	for (size_t component_size : { 1u, 2u, 4u })
	{
		astcenc_type type = ASTCENC_TYPE_F32;
		if (component_size == 1)
		{
			type = ASTCENC_TYPE_U8;
		}
		else if (component_size == 2)
		{
			type = ASTCENC_TYPE_F16;
		}

		size_t row_bytes = 4 * 4 * component_size;
		size_t max_size = std::numeric_limits<size_t>::max();
		const size_t strides[] {
			row_bytes - 1,
			max_size - component_size + 1,
			(max_size / 2 / component_size) * component_size
		};

		for (size_t stride : strides)
		{
			astcenc_image image {
				.dim_x = 4,
				.dim_y = 4,
				.dim_z = 1,
				.data_type = type,
				.data = nullptr,
				.row_stride = stride
			};

			status = astcenc_compress_image(context.get(), &image, &swizzle, nullptr, 0, 0);
			EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);

			status = astcenc_decompress_image(context.get(), nullptr, 0, &image, &swizzle, 0);
			EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);
		}

		// Two rows can overflow only when adding the final row width
		astcenc_image bad_image {
			.dim_x = 4,
			.dim_y = 2,
			.dim_z = 1,
			.data_type = type,
			.data = nullptr,
			.row_stride = (max_size / component_size) * component_size
		};

		status = astcenc_compress_image(context.get(), &bad_image, &swizzle, nullptr, 0, 0);
		EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);

		status = astcenc_decompress_image(context.get(), nullptr, 0, &bad_image, &swizzle, 0);
		EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);

		if (component_size > 1)
		{
			astcenc_image image {
				.dim_x = 4,
				.dim_y = 1,
				.dim_z = 1,
				.data_type = type,
				.data = nullptr,
				.row_stride = row_bytes + 1
			};

			status = astcenc_compress_image(context.get(), &image, &swizzle, nullptr, 0, 0);
			EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);

			status = astcenc_decompress_image(context.get(), nullptr, 0, &image, &swizzle, 0);
			EXPECT_EQ(status, ASTCENC_ERR_BAD_PARAM);
		}
	}
}

}
