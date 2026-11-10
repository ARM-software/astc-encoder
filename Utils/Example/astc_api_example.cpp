// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2021-2026 Arm Limited
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

// Minimal astcenc library example
//
// This sample shows how to include the astcenc library in your CMake project as an external
// dependency, and how to compress and decompress images using the C library API.
//
// The sample uses a minimal command line and fixed compression settings. Comments indicate
// extension points. Error paths log failures and exit without freeing resources, keeping the
// example simple.

#include <limits>
#include <stdio.h>

#include "astcenc.h"

#define STB_IMAGE_IMPLEMENTATION
#include "ThirdParty/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "ThirdParty/stb_image_write.h"

int main(int argc, char **argv)
{
	// Parse command line
	if (argc != 3)
	{
		printf("Usage:\n"
		       "   %s <source> <dest>\n\n"
		       "   <source> : Uncompressed LDR source image.\n"
		       "   <dest>   : Uncompressed LDR destination image (png).\n"
		       , argv[0]);
		return 1;
	}

	// ------------------------------------------------------------------------
	// Use fixed compressor settings for this example
	static const size_t thread_count = 1;
	static const unsigned int block_x = 6;
	static const unsigned int block_y = 6;
	static const unsigned int block_z = 1;
	static const astcenc_profile profile = ASTCENC_PRF_LDR;
	static const float quality = ASTCENC_PRE_MEDIUM;
	static const astcenc_swizzle swizzle {
		ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
	};

	// ------------------------------------------------------------------------
	// Load input image, forcing 4 components
	int image_x, image_y, image_c;
	uint8_t *image_data = (uint8_t*)stbi_load(argv[1], &image_x, &image_y, &image_c, 4);
	if (!image_data)
	{
		printf("Failed to load image \"%s\"\n", argv[1]);
		return 1;
	}

	// Compute the number of ASTC blocks in each dimension
	size_t block_count_x = (static_cast<size_t>(image_x) + block_x - 1) / block_x;
	size_t block_count_y = (static_cast<size_t>(image_y) + block_y - 1) / block_y;

	// ------------------------------------------------------------------------
	// Initialize the default configuration for the block size and quality
	astcenc_config config;
	astcenc_error status;
	status = astcenc_config_init(profile, block_x, block_y, block_z, quality, 0, &config);
	if (status != ASTCENC_SUCCESS)
	{
		printf("ERROR: Codec config init failed: %s\n", astcenc_get_error_string(status));
		return 1;
	}

	// Power users can customize any config settings after calling config_init() and before calling
	// context alloc()

	// ------------------------------------------------------------------------
	// Create a context based on the configuration
	astcenc_context* context;
	status = astcenc_context_alloc(&config, thread_count, &context, nullptr);
	if (status != ASTCENC_SUCCESS)
	{
		printf("ERROR: Codec context alloc failed: %s\n", astcenc_get_error_string(status));
		return 1;
	}

	// ------------------------------------------------------------------------
	// Compress the image
	uint8_t* slices = image_data;
	astcenc_image image {
		.dim_x = static_cast<size_t>(image_x),
		.dim_y = static_cast<size_t>(image_y),
		.dim_z = 1,
		.data_type = ASTCENC_TYPE_U8,
		.data = reinterpret_cast<void**>(&slices),
		.row_stride = 0
	};

	// Space needed for 16 bytes of output per compressed block
	if (block_count_x > std::numeric_limits<size_t>::max() / 16 / block_count_y)
	{
		printf("ERROR: Compressed image size overflow\n");
		return 1;
	}

	size_t comp_len = block_count_x * block_count_y * 16;
	uint8_t* comp_data = new uint8_t[comp_len];

	status = astcenc_compress_image(context, &image, &swizzle, comp_data, comp_len, 0);
	if (status != ASTCENC_SUCCESS)
	{
		printf("ERROR: Codec compress failed: %s\n", astcenc_get_error_string(status));
		return 1;
	}

	// The comp_data array contains the raw compressed data you would pass to the graphics API, or
	// pack into a wrapper format such as a KTX file

	// If using multithreaded compression to sequentially compress multiple images you should reuse
	// the same context, calling the function astcenc_compress_reset() between each image in
	// the series

	// ------------------------------------------------------------------------
	// Decompress the image
	// Reuse the image structure for decompressed output
	status = astcenc_decompress_image(context, comp_data, comp_len, &image, &swizzle, 0);
	if (status != ASTCENC_SUCCESS)
	{
		printf("ERROR: Codec decompress failed: %s\n", astcenc_get_error_string(status));
		return 1;
	}

	// If using multithreaded decompression to sequentially decompress multiple images you should
	// reuse the same context, calling the function astcenc_decompress_reset() between each image in
	// the series

	// ------------------------------------------------------------------------
	// Store the result back to disk
	stbi_write_png(argv[2], image_x, image_y, 4, image_data, 4 * image_x);

	// ------------------------------------------------------------------------
	// Cleanup library resources
	stbi_image_free(image_data);
	astcenc_context_free(context);
	delete[] comp_data;

	return 0;
}
