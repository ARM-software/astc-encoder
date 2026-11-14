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
 * @brief Load and store Morton-ordered ASTC blocks without changing the codec's raster layout.
 */

#include <cstring>
#include <fstream>
#include <new>
#include <utility>

#include "astcenccli_internal.h"

/** @brief The ASTCM header, using the ASTC header layout with a distinct magic. */
struct astcm_header
{
	/** @brief The ASCII format identifier "ASCM". */
	uint8_t magic[4];

	/** @brief The block width in texels. */
	uint8_t block_x;

	/** @brief The block height in texels. */
	uint8_t block_y;

	/** @brief The block depth, which must be one. */
	uint8_t block_z;

	/** @brief The image width as a little-endian 24-bit integer. */
	uint8_t dim_x[3];

	/** @brief The image height as a little-endian 24-bit integer. */
	uint8_t dim_y[3];

	/** @brief The image depth, which must be one. */
	uint8_t dim_z[3];
};

static_assert(sizeof(astcm_header) == 16);

/** @brief The ASTCM file identifier. */
static constexpr uint8_t astcm_magic[4] { 'A', 'S', 'C', 'M' };

/**
 * @brief Decode a little-endian 24-bit integer.
 *
 * @param bytes   The three encoded bytes.
 *
 * @return The decoded value.
 */
static size_t unpack_u24(
	const uint8_t bytes[3]
) {
	return static_cast<size_t>(bytes[0]) |
	       (static_cast<size_t>(bytes[1]) << 8) |
	       (static_cast<size_t>(bytes[2]) << 16);
}

/**
 * @brief Validate the image layout and compute its block counts and payload size.
 *
 * @param      img         The image metadata to validate.
 * @param[out] blocks_x    The number of blocks across the image.
 * @param[out] blocks_y    The number of blocks down the image.
 * @param[out] data_size   The payload size in bytes.
 *
 * @return @c true on error, @c false on success.
 */
static bool get_astcm_layout(
	const astc_compressed_image& img,
	size_t& blocks_x,
	size_t& blocks_y,
	size_t& data_size
) {
	if ((img.block_z != 1) || (img.dim_z != 1))
	{
		print_error("ERROR: ASTCM only supports 2D images and blocks\n");
		return true;
	}

	if ((img.block_x == 0) || (img.block_x > 255) ||
	    (img.block_y == 0) || (img.block_y > 255) ||
	    (img.dim_x == 0) || (img.dim_x > 0xFFFFFF) ||
	    (img.dim_y == 0) || (img.dim_y > 0xFFFFFF))
	{
		print_error("ERROR: Invalid ASTCM image or block dimensions\n");
		return true;
	}

	blocks_x = astc::get_block_count_safe(img.dim_x, img.block_x);
	blocks_y = astc::get_block_count_safe(img.dim_y, img.block_y);
	bool is_overflow = false;
	data_size = astc::mul_safe(blocks_x, blocks_y, is_overflow);
	data_size = astc::mul_safe(data_size, 16, is_overflow);
	if (is_overflow || (data_size > img.data.max_size()))
	{
		print_error("ERROR: ASTCM image size overflow\n");
		return true;
	}

	return false;
}

/**
 * @brief Find the raster byte offsets of a brick's blocks in Morton order.
 *
 * @param      blocks_x   The number of blocks across the image.
 * @param      blocks_y   The number of blocks down the image.
 * @param      brick_x    The brick origin in blocks.
 * @param      brick_y    The brick origin in blocks.
 * @param[out] offsets    The byte offsets of in-range blocks in the raster payload.
 *
 * @return The number of in-range blocks in this brick.
 */
static size_t get_brick_offsets(
	size_t blocks_x,
	size_t blocks_y,
	size_t brick_x,
	size_t brick_y,
	size_t offsets[64]
) {
	size_t count = 0;
	for (size_t morton = 0; morton < 64; morton++)
	{
		// Even Morton bits encode X; odd bits encode Y
		size_t x = brick_x + ((morton & 1) | ((morton >> 1) & 2) | ((morton >> 2) & 4));
		size_t y = brick_y + (((morton >> 1) & 1) | ((morton >> 2) & 2) | ((morton >> 3) & 4));
		if ((x < blocks_x) && (y < blocks_y))
		{
			offsets[count++] = (y * blocks_x + x) * 16;
		}
	}

	return count;
}

/* See header for documentation. */
bool load_astcm_compressed_image(
	const char* filename,
	astc_compressed_image& img
) {
	std::ifstream file(filename, std::ios::binary);
	if (!file)
	{
		print_error("ERROR: File open failed '%s'\n", filename);
		return true;
	}

	astcm_header header;
	file.read(reinterpret_cast<char*>(&header), sizeof(header));
	if (!file || (std::memcmp(header.magic, astcm_magic, sizeof(astcm_magic)) != 0))
	{
		print_error("ERROR: Invalid ASTCM header '%s'\n", filename);
		return true;
	}

	astc_compressed_image loaded {
		.block_x = header.block_x,
		.block_y = header.block_y,
		.block_z = header.block_z,
		.dim_x = unpack_u24(header.dim_x),
		.dim_y = unpack_u24(header.dim_y),
		.dim_z = unpack_u24(header.dim_z),
		.data = {}
	};

	size_t blocks_x;
	size_t blocks_y;
	size_t data_size;
	if (get_astcm_layout(loaded, blocks_x, blocks_y, data_size))
	{
		return true;
	}

	// Check the file length before allocating memory based on untrusted dimensions
	file.seekg(0, std::ios::end);
	std::streamoff file_size = file.tellg();
	std::streamoff header_size = sizeof(header);
	if ((file_size < header_size) || (static_cast<uintmax_t>(file_size - header_size) != data_size))
	{
		print_error("ERROR: ASTCM payload size mismatch '%s'\n", filename);
		return true;
	}

	file.seekg(sizeof(header), std::ios::beg);
	try
	{
		loaded.data.resize(data_size);
	}
	catch (const std::bad_alloc&)
	{
		print_error("ERROR: ASTCM image allocation failed\n");
		return true;
	}

	// Transfer one brick at a time to avoid a second image-sized buffer
	uint8_t brick[64 * 16];
	size_t offsets[64];
	for (size_t y = 0; y < blocks_y; y += 8)
	{
		for (size_t x = 0; x < blocks_x; x += 8)
		{
			size_t count = get_brick_offsets(blocks_x, blocks_y, x, y, offsets);
			file.read(reinterpret_cast<char*>(brick), static_cast<std::streamsize>(count * 16));
			if (!file)
			{
				print_error("ERROR: ASTCM payload read failed '%s'\n", filename);
				return true;
			}

			for (size_t i = 0; i < count; i++)
			{
				std::memcpy(loaded.data.data() + offsets[i], brick + i * 16, 16);
			}
		}
	}

	img = std::move(loaded);
	return false;
}

/* See header for documentation. */
bool write_astcm_compressed_image(
	const astc_compressed_image& img,
	std::ostream& stream
) {
	size_t blocks_x;
	size_t blocks_y;
	size_t data_size;
	if (get_astcm_layout(img, blocks_x, blocks_y, data_size))
	{
		return true;
	}

	if (img.data.size() != data_size)
	{
		print_error("ERROR: ASTCM payload size does not match image dimensions\n");
		return true;
	}

	astcm_header header {};
	std::memcpy(header.magic, astcm_magic, sizeof(astcm_magic));
	header.block_x = static_cast<uint8_t>(img.block_x);
	header.block_y = static_cast<uint8_t>(img.block_y);
	header.block_z = 1;
	header.dim_z[0] = 1;
	for (unsigned int i = 0; i < 3; i++)
	{
		header.dim_x[i] = static_cast<uint8_t>((img.dim_x >> (i * 8)) & 0xFF);
		header.dim_y[i] = static_cast<uint8_t>((img.dim_y >> (i * 8)) & 0xFF);
	}

	stream.write(reinterpret_cast<const char*>(&header), sizeof(header));
	uint8_t brick[64 * 16];
	size_t offsets[64];
	for (size_t y = 0; y < blocks_y; y += 8)
	{
		for (size_t x = 0; x < blocks_x; x += 8)
		{
			size_t count = get_brick_offsets(blocks_x, blocks_y, x, y, offsets);
			for (size_t i = 0; i < count; i++)
			{
				std::memcpy(brick + i * 16, img.data.data() + offsets[i], 16);
			}

			stream.write(reinterpret_cast<const char*>(brick),
			             static_cast<std::streamsize>(count * 16));
			if (!stream)
			{
				return true;
			}
		}
	}

	return stream.fail();
}

/* See header for documentation. */
bool store_astcm_compressed_image(
	const astc_compressed_image& img,
	const char* filename
) {
	// Validate before opening the output to avoid truncating files for unsupported image layouts
	size_t blocks_x;
	size_t blocks_y;
	size_t data_size;
	if (get_astcm_layout(img, blocks_x, blocks_y, data_size))
	{
		return true;
	}

	if (img.data.size() != data_size)
	{
		print_error("ERROR: ASTCM payload size does not match image dimensions\n");
		return true;
	}

	std::ofstream file(filename, std::ios::binary);
	if (!file)
	{
		print_error("ERROR: File open failed '%s'\n", filename);
		return true;
	}

	bool is_write_error = write_astcm_compressed_image(img, file);
	file.close();
	return is_write_error || file.fail();
}
