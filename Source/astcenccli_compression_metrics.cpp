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
 * @brief Measure lossless compression of serialized ASTC and Morton-ordered ASTCM images.
 */

#if defined(ASTCENC_RDO_SIZE)

#include <algorithm>
#include <limits>
#include <sstream>
#include <string_view>

#include <zlib.h>

#include "astcenccli_internal.h"

/**
 * @brief Measure a level 9 raw DEFLATE stream without retaining the compressed output.
 *
 * @param      data              The serialized file to compress.
 * @param[out] comp_size   The compressed size in bytes, excluding ZIP container metadata.
 *
 * @return @c true on error, @c false on success.
 */
static bool get_deflate_size(
	std::string_view data,
	size_t& comp_size
) {
	z_stream stream {};
	int status = deflateInit2(&stream, Z_BEST_COMPRESSION, Z_DEFLATED, -MAX_WBITS,
	                          8, Z_DEFAULT_STRATEGY);
	if (status != Z_OK)
	{
		print_error("ERROR: zlib initialization failed: %d\n", status);
		return true;
	}

	uint8_t output[64 * 1024];
	size_t offset = 0;
	comp_size = 0;
	do
	{
		// Chunk input to fit zlib's unsigned int buffer lengths
		if (stream.avail_in == 0)
		{
			size_t count = std::min(data.size() - offset,
			                        static_cast<size_t>(std::numeric_limits<uInt>::max()));
			stream.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(data.data() + offset));
			stream.avail_in = static_cast<uInt>(count);
			offset += count;
		}

		int flush = Z_NO_FLUSH;
		if (offset == data.size())
		{
			flush = Z_FINISH;
		}

		stream.next_out = output;
		stream.avail_out = sizeof(output);
		status = deflate(&stream, flush);
		comp_size += sizeof(output) - stream.avail_out;
	} while (status == Z_OK);

	int end_status = deflateEnd(&stream);
	if ((status != Z_STREAM_END) || (end_status != Z_OK))
	{
		print_error("ERROR: zlib compression failed: %d\n", status);
		return true;
	}

	return false;
}

/* See header for documentation. */
bool print_compression_metrics(
	const astc_compressed_image& img,
	bool is_morton
) {
	// Use file serializers so headers, block order, and partial tiles match the selected format
	std::ostringstream stream(std::ios::out | std::ios::binary);
	const char* format = "ASTC";
	bool is_write_error;
	if (is_morton)
	{
		format = "ASTCM";
		is_write_error = write_astcm_compressed_image(img, stream);
	}
	else
	{
		is_write_error = write_cimage(img, stream);
	}

	if (is_write_error)
	{
		print_error("ERROR: Failed to serialize %s for compression metrics\n", format);
		return true;
	}

	std::string_view data = stream.view();
	size_t comp_size;
	bool is_compression_error = get_deflate_size(data, comp_size);
	if (is_compression_error)
	{
		return true;
	}

	printf("\nStorage metrics\n");
	printf("===============\n\n");

	const char* file_label = "ASTC file size:";
	const char* deflate_label = "ASTC DEFLATE size:";
	if (is_morton)
	{
		file_label = "ASTCM file size:";
		deflate_label = "ASTCM DEFLATE size:";
	}

	double data_sized = static_cast<double>(data.size()) / 1024.0;
	double comp_sized = static_cast<double>(comp_size) / 1024.0;
	double comp_percent = 100.0 * comp_sized / data_sized;

	printf("    %-27s%8.1f KB\n",          file_label,    data_sized);
	printf("    %-27s%8.1f KB (%.2f%%)\n", deflate_label, comp_sized, comp_percent);
	printf("    DEFLATE level:             %8d\n", Z_BEST_COMPRESSION);
	return false;
}

#endif
