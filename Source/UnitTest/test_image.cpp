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
 * @brief Unit tests for strided image storage.
 */

#include <algorithm>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>

#include "gtest/gtest.h"

#include "../astcenc.h"

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
