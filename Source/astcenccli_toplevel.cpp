// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2011-2020 Arm Limited
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
 * @brief Functions for codec library front-end.
 */

#include "astcenc.h"
#include "astcenc_internal.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <vector>

#ifdef DEBUG_CAPTURE_NAN
	#ifndef _GNU_SOURCE
		#define _GNU_SOURCE
	#endif

	#include <fenv.h>
#endif

static double start_time;
static double end_time;
static double start_coding_time;
static double end_coding_time;

static const uint32_t MAGIC_FILE_CONSTANT = 0x5CA1AB13;

struct astc_header
{
	uint8_t magic[4];
	uint8_t blockdim_x;
	uint8_t blockdim_y;
	uint8_t blockdim_z;
	uint8_t xsize[3];			// x-size = xsize[0] + xsize[1] + xsize[2]
	uint8_t ysize[3];			// x-size, y-size and z-size are given in texels;
	uint8_t zsize[3];			// block count is inferred
};

struct astc_compressed_image
{
	int block_x;
	int block_y;
	int block_z;
	int dim_x;
	int dim_y;
	int dim_z;
	uint8_t* data;
};

static uint32_t unpack_bytes(
	uint8_t a,
	uint8_t b,
	uint8_t c,
	uint8_t d
) {
	return ((uint32_t)(a))       +
	       ((uint32_t)(b) << 8)  +
	       ((uint32_t)(c) << 16) +
	       ((uint32_t)(d) << 24);
}


int load_astc_file(
	const char *filename,
	astc_compressed_image& out_image
) {
 	std::ifstream file(filename, std::ios::in | std::ios::binary);
	if (!file)
	{
		printf("ERROR: File open failed '%s'\n", filename);
		return 1;
	}

	astc_header hdr;
	file.read((char*)&hdr, sizeof(astc_header));
	if (!file)
	{
		printf("ERROR: File read failed '%s'\n", filename);
		return 1;
	}

	uint32_t magicval = unpack_bytes(hdr.magic[0], hdr.magic[1], hdr.magic[2], hdr.magic[3]);
	if (magicval != MAGIC_FILE_CONSTANT)
	{
		printf("ERROR: File not recognized '%s'\n", filename);
		return 1;
	}

	int block_x = hdr.blockdim_x;
	int block_y = hdr.blockdim_y;
	int block_z = hdr.blockdim_z;

	if (((block_z == 1) && !is_legal_2d_block_size(block_x, block_y)) ||
	    ((block_z > 1) && !is_legal_3d_block_size(block_x, block_y, block_z)))
	{
		printf("ERROR: File corrupt '%s'\n", filename);
		return 1;
	}

	int xsize = unpack_bytes(hdr.xsize[0], hdr.xsize[1], hdr.xsize[2], 0);
	int ysize = unpack_bytes(hdr.ysize[0], hdr.ysize[1], hdr.ysize[2], 0);
	int zsize = unpack_bytes(hdr.zsize[0], hdr.zsize[1], hdr.zsize[2], 0);

	if (xsize == 0 || ysize == 0 || zsize == 0)
	{
		printf("ERROR: File corrupt '%s'\n", filename);
		return 1;
	}

	int xblocks = (xsize + block_x - 1) / block_x;
	int yblocks = (ysize + block_y - 1) / block_y;
	int zblocks = (zsize + block_z - 1) / block_z;

	size_t data_size = xblocks * yblocks * zblocks * 16;
	uint8_t *buffer = new uint8_t[data_size];

	file.read((char*)buffer, data_size);
	if (!file)
	{
		printf("ERROR: File read failed '%s'\n", filename);
		return 1;
	}

	out_image.data = buffer;
	out_image.block_x = block_x;
	out_image.block_y = block_y;
	out_image.block_z = block_z;
	out_image.dim_x = xsize;
	out_image.dim_y = ysize;
	out_image.dim_z = zsize;
	return 0;
}


int store_astc_file(
	const astc_compressed_image& comp_img,
	const char* filename
) {
	int xblocks = (comp_img.dim_x + comp_img.block_x - 1) / comp_img.block_x;
	int yblocks = (comp_img.dim_y + comp_img.block_y - 1) / comp_img.block_y;
	int zblocks = (comp_img.dim_z + comp_img.block_z - 1) / comp_img.block_z;
	size_t data_bytes = xblocks * yblocks * zblocks * 16;

	astc_header hdr;
	hdr.magic[0] = MAGIC_FILE_CONSTANT & 0xFF;
	hdr.magic[1] = (MAGIC_FILE_CONSTANT >> 8) & 0xFF;
	hdr.magic[2] = (MAGIC_FILE_CONSTANT >> 16) & 0xFF;
	hdr.magic[3] = (MAGIC_FILE_CONSTANT >> 24) & 0xFF;
	hdr.blockdim_x = comp_img.block_x;
	hdr.blockdim_y = comp_img.block_y;
	hdr.blockdim_z = comp_img.block_z;
	hdr.xsize[0] =  comp_img.dim_x & 0xFF;
	hdr.xsize[1] = (comp_img.dim_x >> 8) & 0xFF;
	hdr.xsize[2] = (comp_img.dim_x >> 16) & 0xFF;
	hdr.ysize[0] =  comp_img.dim_y & 0xFF;
	hdr.ysize[1] = (comp_img.dim_y >> 8) & 0xFF;
	hdr.ysize[2] = (comp_img.dim_y >> 16) & 0xFF;
	hdr.zsize[0] =  comp_img.dim_z & 0xFF;
	hdr.zsize[1] = (comp_img.dim_z >> 8) & 0xFF;
	hdr.zsize[2] = (comp_img.dim_z >> 16) & 0xFF;

 	std::ofstream file(filename, std::ios::out | std::ios::binary);
	if (!file)
	{
		printf("ERROR: File open failed '%s'\n", filename);
		return 1;
	}

	file.write((char*)&hdr, sizeof(astc_header));
	file.write((char*)comp_img.data, data_bytes);
	return 0;
}

// The ASTC codec is written with the assumption that a float threaded through
// the "if32" union will in fact be stored and reloaded as a 32-bit IEEE-754 single-precision
// float, stored with round-to-nearest rounding. This is always the case in an
// IEEE-754 compliant system, however not every system is actually IEEE-754 compliant
// in the first place. As such, we run a quick test to check that this is actually the case
// (e.g. gcc on 32-bit x86 will typically fail unless -msse2 -mfpmath=sse2 is specified).

// TODO: Move to mathlib
volatile float xprec_testval = 2.51f;
static void test_inappropriate_extended_precision()
{
	if32 p;
	p.f = xprec_testval + 12582912.0f;
	float q = p.f - 12582912.0f;
	if (q != 3.0f)
	{
		printf("CPU support error: float math is not IEEE-754 compliant.\n");
		printf("    Please recompile with IEEE-754 support.\n");
		exit(1);
	}
}

static void test_inappropriate_cpu_extensions()
{
	#if ASTC_SSE >= 42
		if (!cpu_supports_sse42()) {
			printf("CPU support error: host lacks SSE 4.2 support.\n");
			printf("    Please recompile with VEC=sse2.\n");
			exit(1);
		}
	#endif

	#if ASTC_POPCNT >= 1
		if (!cpu_supports_popcnt()) {
			printf("CPU support error: host lacks POPCNT support.\n");
			printf("    Please recompile with VEC=sse2.\n");
			exit(1);
		}
	#endif

	#if ASTC_AVX >= 2
		if (!cpu_supports_avx2()) {
			printf("CPU support error: host lacks AVX2 support.\n");
			printf("    Please recompile with VEC=sse4.2 or VEC=sse2.\n");
			exit(1);
		}
	#endif
}

/**
 * @brief Utility to generate a slice file name from a pattern.
 *
 * Convert "foo/bar.png" in to "foo/bar_<slice>.png"
 *
 * @param basename The base pattern; must contain a file extension.
 * @param index    The slice index.
 * @param error    Set to true on success, false on error (no extension found).
 *
 * @return The slice file name.
 */
static std::string get_slice_filename(
	const std::string& basename,
	unsigned int index,
	bool& error
) {
	error = false;

	size_t sep = basename.find_last_of(".");
	if (sep == std::string::npos)
	{
		error = true;
		return "";
	}

	std::string base = basename.substr(0, sep);
	std::string ext = basename.substr(sep);

	std::string name = base + "_" + std::to_string(index) + ext;
	return name;
}

/**
 * @brief Load a non-astc image file from memory.
 *
 * @param filename            The file to load, or a pattern for array loads.
 * @param dim_z               The number of slices to load.
 * @param padding             The number of texels of padding.
 * @param y_flip              Should this image be Y flipped?
 * @param linearize_srgb      Should this image be converted to linear from sRGB?
 * @param[out] is_hdr         Is the loaded image HDR?
 * @param[out] num_components The number of components in the loaded image.
 *
 * @return The astc image file, or nullptr on error.
 */
static astc_codec_image* load_uncomp_file(
	const char* filename,
	unsigned int dim_z,
	int padding,
	bool y_flip,
	bool linearize_srgb,
	bool& is_hdr,
	int& num_components
) {
	astc_codec_image *image = nullptr;

	// For a 2D image just load the image directly
	if (dim_z == 1)
	{
		image = astc_codec_load_image(filename, padding, y_flip, linearize_srgb,
		                              is_hdr, num_components);
	}
	else
	{
		bool slice_is_hdr;
		int slice_num_components;
		astc_codec_image* slice = nullptr;
		std::vector<astc_codec_image*> slices;

		// For a 3D image load an array of slices
		for (unsigned int image_index = 0; image_index < dim_z; image_index++)
		{
			bool error;
			std::string slice_name = get_slice_filename(filename, image_index, error);
			if (error)
			{
				printf("ERROR: Image pattern does not contain an extension: %s\n", filename);
				break;
			}

			slice = astc_codec_load_image(slice_name.c_str(), padding, y_flip, linearize_srgb,
			                              slice_is_hdr, slice_num_components);
			if (!slice)
			{
				break;
			}

			slices.push_back(slice);

			// Check it is not a 3D image
			if (slice->zsize != 1)
			{
				printf("ERROR: Image arrays do not support 3D sources: %s\n", slice_name.c_str());
				break;
			}

			// Check slices are consistent with each other
			if (image_index != 0)
			{
				if ((is_hdr != slice_is_hdr) || (num_components != slice_num_components))
				{
					printf("ERROR: Image array[0] and [%d] are different formats\n", image_index);
					break;
				}

				if ((slices[0]->xsize != slice->xsize) ||
				    (slices[0]->ysize != slice->ysize) ||
				    (slices[0]->zsize != slice->zsize))
				{
					printf("ERROR: Image array[0] and [%d] are different dimensions\n", image_index);
					break;
				}
			}
			else
			{
				is_hdr = slice_is_hdr;
				num_components = slice_num_components;
			}
		}

		// If all slices loaded correctly then repack them into a single image
		if (slices.size() == dim_z)
		{
			unsigned int dim_x = slices[0]->xsize;
			unsigned int dim_y = slices[0]->ysize;
			int bitness = is_hdr ? 16 : 8;
			int slice_size = (dim_x + (2 * padding)) * (dim_y + (2 * padding));

			image = alloc_image(bitness, dim_x, dim_y, dim_y, padding);

			// Combine 2D source images into one 3D image; skipping padding slices
			for (unsigned int z = padding; z < dim_z + padding; z++)
			{
				if (bitness == 8)
				{
					size_t copy_size = slice_size * 4 * sizeof(uint8_t);
					memcpy(*image->data8[z], *slices[z - padding]->data8[0], copy_size);
				}
				else
				{
					size_t copy_size = slice_size * 4 * sizeof(uint16_t);
					memcpy(*image->data16[z], *slices[z - padding]->data16[0], copy_size);
				}
			}

			// Fill in the padding slices with clamped data
			fill_image_padding_area(image);
		}

		for (auto &i : slices)
		{
			free_image(i);
		}
	}

	return image;
}


int astc_main(
	int argc,
	char **argv
) {
	// TODO: These should move into the backend code and fail on context_alloc.
	test_inappropriate_extended_precision();
	test_inappropriate_cpu_extensions();

	start_time = get_time();

	#ifdef DEBUG_CAPTURE_NAN
		feenableexcept(FE_DIVBYZERO | FE_INVALID);
	#endif

	if (argc < 2)
	{
		astcenc_print_shorthelp();
		return 0;
	}

	enum astc_op_mode {
		ASTC_ENCODE,
		ASTC_DECODE,
		ASTC_ENCODE_AND_DECODE,
		ASTC_PRINT_LONGHELP,
		ASTC_PRINT_VERSION,
		ASTC_UNRECOGNIZED
	};

	astcenc_profile decode_mode = ASTCENC_PRF_LDR;
	astc_op_mode op_mode = ASTC_UNRECOGNIZED;

	struct {
		const char *opt;
		astc_op_mode op_mode;
		astcenc_profile decode_mode;
	} modes[] = {
		{"-cl",      ASTC_ENCODE,            ASTCENC_PRF_LDR},
		{"-dl",      ASTC_DECODE,            ASTCENC_PRF_LDR},
		{"-tl",      ASTC_ENCODE_AND_DECODE, ASTCENC_PRF_LDR},
		{"-cs",      ASTC_ENCODE,            ASTCENC_PRF_LDR_SRGB},
		{"-ds",      ASTC_DECODE,            ASTCENC_PRF_LDR_SRGB},
		{"-ts",      ASTC_ENCODE_AND_DECODE, ASTCENC_PRF_LDR_SRGB},
		{"-ch",      ASTC_ENCODE,            ASTCENC_PRF_HDR},
		{"-dh",      ASTC_DECODE,            ASTCENC_PRF_HDR},
		{"-th",      ASTC_ENCODE_AND_DECODE, ASTCENC_PRF_HDR},
		{"-cH",      ASTC_ENCODE,            ASTCENC_PRF_HDR_RGB_LDR_A},
		{"-dH",      ASTC_DECODE,            ASTCENC_PRF_HDR_RGB_LDR_A},
		{"-tH",      ASTC_ENCODE_AND_DECODE, ASTCENC_PRF_HDR_RGB_LDR_A},
		{"-h",       ASTC_PRINT_LONGHELP,    ASTCENC_PRF_HDR},
		{"-help",    ASTC_PRINT_LONGHELP,    ASTCENC_PRF_HDR},
		{"-v",       ASTC_PRINT_VERSION,     ASTCENC_PRF_HDR},
		{"-version", ASTC_PRINT_VERSION,     ASTCENC_PRF_HDR}
	};

	int modes_count = sizeof(modes) / sizeof(modes[0]);
	for (int i = 0; i < modes_count; i++)
	{
		if (!strcmp(argv[1], modes[i].opt))
		{
			op_mode = modes[i].op_mode;
			decode_mode = modes[i].decode_mode;
			break;
		}
	}

	switch (op_mode)
	{
		case ASTC_UNRECOGNIZED:
			printf("ERROR: Unrecognized operation '%s'\n", argv[1]);
			return 1;
		case ASTC_PRINT_LONGHELP:
			astcenc_print_longhelp();
			return 0;
		case ASTC_PRINT_VERSION:
			astcenc_print_header();
			return 0;
		default:
			break;
	}

	int array_size = 1;

	const char *input_filename = argc >= 3 ? argv[2] : nullptr;
	const char *output_filename = argc >= 4 ? argv[3] : nullptr;

	int silentmode = 0;
	int y_flip = 0;
	int linearize_srgb = 0;

	error_weighting_params ewp;

	ewp.rgb_power = 1.0f;
	ewp.alpha_power = 1.0f;
	ewp.rgb_base_weight = 1.0f;
	ewp.alpha_base_weight = 1.0f;
	ewp.rgb_mean_weight = 0.0f;
	ewp.rgb_stdev_weight = 0.0f;
	ewp.alpha_mean_weight = 0.0f;
	ewp.alpha_stdev_weight = 0.0f;

	ewp.rgb_mean_and_stdev_mixing = 0.0f;
	ewp.mean_stdev_radius = 0;
	ewp.enable_rgb_scale_with_alpha = 0;
	ewp.alpha_radius = 0;

	ewp.block_artifact_suppression = 0.0f;
	ewp.rgba_weights[0] = 1.0f;
	ewp.rgba_weights[1] = 1.0f;
	ewp.rgba_weights[2] = 1.0f;
	ewp.rgba_weights[3] = 1.0f;
	ewp.ra_normal_angular_scale = 0;

	astcenc_swizzle swz_encode {
		ASTCENC_SWZ_R,
		ASTCENC_SWZ_G,
		ASTCENC_SWZ_B,
		ASTCENC_SWZ_A
	};

	astcenc_swizzle swz_decode {
		ASTCENC_SWZ_R,
		ASTCENC_SWZ_G,
		ASTCENC_SWZ_B,
		ASTCENC_SWZ_A
	};

	int thread_count = 0;		// default value

	int plimit_set = -1;
	float dblimit_set = 0.0f;
	float oplimit_set = 0.0f;
	float mincorrel_set = 0.0f;
	float bmc_set = 0.0f;
	int maxiters_set = 0;

	int block_x = 0;
	int block_y = 0;
	int block_z = 1;
	float log10_texels = 0.0f;

	int low_fstop = -10;
	int high_fstop = 10;

	// parse the command line's encoding options.
	int argidx;
	if (op_mode == ASTC_ENCODE || op_mode == ASTC_ENCODE_AND_DECODE)
	{
		if (argc < 5)
		{
			printf("ERROR: Block size not specified\n");
			return 1;
		}

		if (argc < 6)
		{
			printf("ERROR: Search quality preset not specified\n");
			return 1;
		}

		int cnt2D, cnt3D;
		int dimensions = sscanf(argv[4], "%dx%d%nx%d%n", &block_x, &block_y, &cnt2D, &block_z, &cnt3D);
		switch (dimensions)
		{
		case 2:
			// Character after the last match should be a NUL
			if (argv[4][cnt2D] || !is_legal_2d_block_size(block_x, block_y))
			{
				printf("ERROR: Invalid block size %s specified\n", argv[4]);
				return 1;
			}

			assert(block_z == 1);
			break;
		case 3:
			// Character after the last match should be a NUL
			if (argv[4][cnt3D] || !is_legal_3d_block_size(block_x, block_y, block_z))
			{
				printf("ERROR: Invalid block size %s specified\n", argv[4]);
				return 1;
			}
			break;
		default:
			printf("ERROR: Invalid block size %s specified\n", argv[4]);
			return 1;
		}

		log10_texels = logf((float)(block_x * block_y * block_z)) / logf(10.0f);

		if (!strcmp(argv[5], "-fast"))
		{
			plimit_set = 4;
			oplimit_set = 1.0f;
			mincorrel_set = 0.5f;
			dblimit_set = MAX(85.0f - 35.0f * log10_texels, 63.0f - 19.0f * log10_texels);
			bmc_set = 50.0f;
			maxiters_set = 1;
		}
		else if (!strcmp(argv[5], "-medium"))
		{
			plimit_set = 25;
			oplimit_set = 1.2f;
			mincorrel_set = 0.75f;
			dblimit_set = MAX(95.0f - 35.0f * log10_texels, 70.0f - 19.0f * log10_texels);
			bmc_set = 75.0f;
			maxiters_set = 2;
		}
		else if (!strcmp(argv[5], "-thorough"))
		{
			plimit_set = 100;
			oplimit_set = 2.5f;
			mincorrel_set = 0.95f;
			dblimit_set = MAX(105.0f - 35.0f * log10_texels, 77.0f - 19.0f * log10_texels);
			bmc_set = 95.0f;
			maxiters_set = 4;
		}
		else if (!strcmp(argv[5], "-exhaustive"))
		{
			plimit_set = PARTITION_COUNT;
			oplimit_set = 1000.0f;
			mincorrel_set = 0.99f;
			dblimit_set = 999.0f;
			bmc_set = 100.0f;
			maxiters_set = 4;
		}
		else
		{
			printf("ERROR: Unknown search preset\n");
			return 1;
		}

		if (decode_mode == ASTCENC_PRF_HDR)
		{
			ewp.rgb_power = 0.75f;
			ewp.rgb_base_weight = 0.0f;
			ewp.rgb_mean_weight = 1.0f;

			ewp.alpha_power = 0.75f;
			ewp.alpha_base_weight = 0.0f;
			ewp.alpha_mean_weight = 1.0f;

			dblimit_set = 999.0f;
		}
		else if (decode_mode == ASTCENC_PRF_HDR_RGB_LDR_A)
		{
			ewp.rgb_power = 0.75f;
			ewp.rgb_base_weight = 0.0f;
			ewp.rgb_mean_weight = 1.0f;

			ewp.alpha_base_weight = 0.05f;

			dblimit_set = 999.0f;
		}

		argidx = 6;
	}
	else
	{
		// For decode the block size and preset are not needed.
		argidx = 4;
	}

	while (argidx < argc)
	{
		if (!strcmp(argv[argidx], "-silent"))
		{
			argidx++;
			silentmode = 1;
		}
		else if (!strcmp(argv[argidx], "-v"))
		{
			argidx += 7;
			if (argidx > argc)
			{
				printf("-v switch with less than 6 arguments, quitting\n");
				return 1;
			}
			ewp.mean_stdev_radius = atoi(argv[argidx - 6]);
			ewp.rgb_power = static_cast<float>(atof(argv[argidx - 5]));
			ewp.rgb_base_weight = static_cast<float>(atof(argv[argidx - 4]));
			ewp.rgb_mean_weight = static_cast<float>(atof(argv[argidx - 3]));
			ewp.rgb_stdev_weight = static_cast<float>(atof(argv[argidx - 2]));
			ewp.rgb_mean_and_stdev_mixing = static_cast<float>(atof(argv[argidx - 1]));
		}
		else if (!strcmp(argv[argidx], "-va"))
		{
			argidx += 5;
			if (argidx > argc)
			{
				printf("-va switch with less than 4 arguments, quitting\n");
				return 1;
			}
			ewp.alpha_power = static_cast<float>(atof(argv[argidx - 4]));
			ewp.alpha_base_weight = static_cast<float>(atof(argv[argidx - 3]));
			ewp.alpha_mean_weight = static_cast<float>(atof(argv[argidx - 2]));
			ewp.alpha_stdev_weight = static_cast<float>(atof(argv[argidx - 1]));
		}
		else if (!strcmp(argv[argidx], "-cw"))
		{
			argidx += 5;
			if (argidx > argc)
			{
				printf("-cw switch with less than 4 arguments\n");
				return 1;
			}
			ewp.rgba_weights[0] = static_cast<float>(atof(argv[argidx - 4]));
			ewp.rgba_weights[1] = static_cast<float>(atof(argv[argidx - 3]));
			ewp.rgba_weights[2] = static_cast<float>(atof(argv[argidx - 2]));
			ewp.rgba_weights[3] = static_cast<float>(atof(argv[argidx - 1]));
		}
		else if (!strcmp(argv[argidx], "-a"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-a switch with no argument\n");
				return 1;
			}
			ewp.enable_rgb_scale_with_alpha = 1;
			ewp.alpha_radius = atoi(argv[argidx - 1]);
		}
		else if (!strcmp(argv[argidx], "-b"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-b switch with no argument\n");
				return 1;
			}
			ewp.block_artifact_suppression = static_cast<float>(atof(argv[argidx - 1]));
		}
		else if (!strcmp(argv[argidx], "-esw"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-esw switch with no argument\n");
				return 1;
			}

			if (strlen(argv[argidx - 1]) != 4)
			{
				printf("Swizzle pattern for the -esw switch must have exactly 4 characters\n");
				return 1;
			}

			astcenc_swz swizzle_components[4];
			for (int i = 0; i < 4; i++)
			{
				switch (argv[argidx - 1][i])
				{
				case 'r':
					swizzle_components[i] = ASTCENC_SWZ_R;
					break;
				case 'g':
					swizzle_components[i] = ASTCENC_SWZ_G;
					break;
				case 'b':
					swizzle_components[i] = ASTCENC_SWZ_B;
					break;
				case 'a':
					swizzle_components[i] = ASTCENC_SWZ_A;
					break;
				case '0':
					swizzle_components[i] = ASTCENC_SWZ_0;
					break;
				case '1':
					swizzle_components[i] = ASTCENC_SWZ_1;
					break;
				default:
					printf("Character '%c' is not a valid swizzle-character\n", argv[argidx - 1][i]);
					return 1;
				}
			}
			swz_encode.r = swizzle_components[0];
			swz_encode.g = swizzle_components[1];
			swz_encode.b = swizzle_components[2];
			swz_encode.a = swizzle_components[3];
		}
		else if (!strcmp(argv[argidx], "-dsw"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-dsw switch with no argument\n");
				return 1;
			}

			if (strlen(argv[argidx - 1]) != 4)
			{
				printf("Swizzle pattern for the -dsw switch must have exactly 4 characters\n");
				return 1;
			}

			astcenc_swz swizzle_components[4];
			for (int i = 0; i < 4; i++)
			{
				switch (argv[argidx - 1][i])
				{
				case 'r':
					swizzle_components[i] = ASTCENC_SWZ_R;
					break;
				case 'g':
					swizzle_components[i] = ASTCENC_SWZ_G;
					break;
				case 'b':
					swizzle_components[i] = ASTCENC_SWZ_B;
					break;
				case 'a':
					swizzle_components[i] = ASTCENC_SWZ_A;
					break;
				case '0':
					swizzle_components[i] = ASTCENC_SWZ_0;
					break;
				case '1':
					swizzle_components[i] = ASTCENC_SWZ_1;
					break;
				case 'z':
					swizzle_components[i] =  ASTCENC_SWZ_Z;
					break;
				default:
					printf("Character '%c' is not a valid swizzle-character\n", argv[argidx - 1][i]);
					return 1;
				}
			}
			swz_decode.r = swizzle_components[0];
			swz_decode.g = swizzle_components[1];
			swz_decode.b = swizzle_components[2];
			swz_decode.a = swizzle_components[3];
		}
		// presets begin here
		else if (!strcmp(argv[argidx], "-normal_psnr"))
		{
			argidx++;
			ewp.rgba_weights[0] = 1.0f;
			ewp.rgba_weights[1] = 0.0f;
			ewp.rgba_weights[2] = 0.0f;
			ewp.rgba_weights[3] = 1.0f;
			ewp.ra_normal_angular_scale = 1;

			swz_encode.r = ASTCENC_SWZ_R;
			swz_encode.g = ASTCENC_SWZ_R;
			swz_encode.b = ASTCENC_SWZ_R;
			swz_encode.a = ASTCENC_SWZ_G;

			swz_decode.r = ASTCENC_SWZ_R;
			swz_decode.g = ASTCENC_SWZ_A;
			swz_decode.b = ASTCENC_SWZ_Z;
			swz_decode.a = ASTCENC_SWZ_1;

			oplimit_set = 1000.0f;
			mincorrel_set = 0.99f;
		}
		else if (!strcmp(argv[argidx], "-normal_percep"))
		{
			argidx++;
			ewp.rgba_weights[0] = 1.0f;
			ewp.rgba_weights[1] = 0.0f;
			ewp.rgba_weights[2] = 0.0f;
			ewp.rgba_weights[3] = 1.0f;
			ewp.ra_normal_angular_scale = 1;

			swz_encode.r = ASTCENC_SWZ_R;
			swz_encode.g = ASTCENC_SWZ_R;
			swz_encode.b = ASTCENC_SWZ_R;
			swz_encode.a = ASTCENC_SWZ_G;

			swz_decode.r = ASTCENC_SWZ_R;
			swz_decode.g = ASTCENC_SWZ_A;
			swz_decode.b = ASTCENC_SWZ_Z;
			swz_decode.a = ASTCENC_SWZ_1;

			oplimit_set = 1000.0f;
			mincorrel_set = 0.99f;
			dblimit_set = 999.0f;

			ewp.block_artifact_suppression = 1.8f;
			ewp.mean_stdev_radius = 3;
			ewp.rgb_mean_weight = 0.0f;
			ewp.rgb_stdev_weight = 50.0f;
			ewp.rgb_mean_and_stdev_mixing = 0.0f;
			ewp.alpha_mean_weight = 0.0f;
			ewp.alpha_stdev_weight = 50.0f;
		}
		else if (!strcmp(argv[argidx], "-mask"))
		{
			argidx++;
			ewp.mean_stdev_radius = 3;
			ewp.rgb_mean_weight = 0.0f;
			ewp.rgb_stdev_weight = 25.0f;
			ewp.rgb_mean_and_stdev_mixing = 0.03f;
			ewp.alpha_mean_weight = 0.0f;
			ewp.alpha_stdev_weight = 25.0f;
		}
		else if (!strcmp(argv[argidx], "-blockmodelimit"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-blockmodelimit switch with no argument\n");
				return 1;
			}
			float cutoff = (float)atof(argv[argidx - 1]);
			if (cutoff > 100.0f || !(cutoff >= 0.0f))
				cutoff = 100.0f;
			bmc_set = cutoff;
		}
		else if (!strcmp(argv[argidx], "-partitionlimit"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-partitionlimit switch with no argument\n");
				return 1;
			}
			plimit_set = atoi(argv[argidx - 1]);
		}
		else if (!strcmp(argv[argidx], "-dblimit"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-dblimit switch with no argument\n");
				return 1;
			}
			dblimit_set = static_cast<float>(atof(argv[argidx - 1]));
		}
		else if (!strcmp(argv[argidx], "-partitionearlylimit"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-partitionearlylimit switch with no argument\n");
				return 1;
			}
			oplimit_set = static_cast<float>(atof(argv[argidx - 1]));
		}
		else if (!strcmp(argv[argidx], "-planecorlimit"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-planecorlimit switch with no argument\n");
				return 1;
			}
			mincorrel_set = static_cast<float>(atof(argv[argidx - 1]));
		}
		else if (!strcmp(argv[argidx], "-refinementlimit"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-refinementlimit switch with no argument\n");
				return 1;
			}
			maxiters_set = atoi(argv[argidx - 1]);
		}
		else if (!strcmp(argv[argidx], "-j"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-j switch with no argument\n");
				return 1;
			}
			thread_count = atoi(argv[argidx - 1]);
		}
		else if (!strcmp(argv[argidx], "-linsrgb"))
		{
			argidx++;
			linearize_srgb = 1;
		}
		else if (!strcmp(argv[argidx], "-yflip"))
		{
			argidx++;
			y_flip = 1;
		}
		else if (!strcmp(argv[argidx], "-mpsnr"))
		{
			argidx += 3;
			if (argidx > argc)
			{
				printf("-mpsnr switch with less than 2 arguments\n");
				return 1;
			}
			low_fstop = atoi(argv[argidx - 2]);
			high_fstop = atoi(argv[argidx - 1]);
			if (high_fstop < low_fstop)
			{
				printf("For -mpsnr switch, the <low> argument cannot be greater than the\n" "high argument.\n");
				return 1;
			}
		}
		// Option: Encode a 3D image from an array of 2D images.
		else if (!strcmp(argv[argidx], "-array"))
		{
			// Only supports compressing (not decompressing or comparison).
			if (op_mode != ASTC_ENCODE)
			{
				printf("-array switch given when not compressing files - decompression and comparison of arrays not supported.\n");
				return 1;
			}

			// Image depth must be specified.
			if (argidx + 2 > argc)
			{
				printf("-array switch given, but no array size (image depth) given.\n");
				return 1;
			}
			argidx++;

			// Read array size (image depth).
			if (!sscanf(argv[argidx], "%d", &array_size) || array_size == 0)
			{
				printf("Invalid array size (image depth) given with -array option: \"%s\".\n", argv[argidx]);
				return 1;
			}
			argidx++;
		}
		else
		{
			printf("Commandline argument \"%s\" not recognized\n", argv[argidx]);
			return 1;
		}
	}

	if (!input_filename)
	{
		printf("ERROR: Input file not specified\n");
		return 1;
	}

	if (!output_filename)
	{
		printf("ERROR: Output file not specified\n");
		return 1;
	}

	if (op_mode == ASTC_ENCODE || op_mode == ASTC_ENCODE_AND_DECODE)
	{
		ewp.max_refinement_iters = maxiters_set;
		ewp.block_mode_cutoff = bmc_set / 100.0f;

		ewp.texel_avg_error_limit = 0.0f;
		if ((decode_mode == ASTCENC_PRF_LDR) || (decode_mode == ASTCENC_PRF_LDR_SRGB))
		{
			ewp.texel_avg_error_limit = powf(0.1f, dblimit_set * 0.1f) * 65535.0f * 65535.0f;
		}

		ewp.partition_1_to_2_limit = oplimit_set;
		ewp.lowest_correlation_cutoff = mincorrel_set;
		ewp.partition_search_limit = astc::clampi(plimit_set, 1, PARTITION_COUNT);

		if (thread_count < 1)
		{
			thread_count = get_cpu_count();
		}

		// Specifying the error weight of a color component as 0 is not allowed.
		// If weights are 0, then they are instead set to a small positive value.
		float max_color_component_weight = MAX(MAX(ewp.rgba_weights[0], ewp.rgba_weights[1]),
		                                       MAX(ewp.rgba_weights[2], ewp.rgba_weights[3]));
		ewp.rgba_weights[0] = MAX(ewp.rgba_weights[0], max_color_component_weight / 1000.0f);
		ewp.rgba_weights[1] = MAX(ewp.rgba_weights[1], max_color_component_weight / 1000.0f);
		ewp.rgba_weights[2] = MAX(ewp.rgba_weights[2], max_color_component_weight / 1000.0f);
		ewp.rgba_weights[3] = MAX(ewp.rgba_weights[3], max_color_component_weight / 1000.0f);

		// print all encoding settings unless specifically told otherwise.
		if (!silentmode)
		{
			printf("Compressor settings\n");
			printf("===================\n\n");

			switch(decode_mode)
			{
				case ASTCENC_PRF_LDR:
					printf("    Color profile:              LDR linear\n");
					break;
				case ASTCENC_PRF_LDR_SRGB:
					printf("    Color profile:              LDR sRGB\n");
					break;
				case ASTCENC_PRF_HDR_RGB_LDR_A:
					printf("    Color profile:              HDR RGB + LDR A\n");
					break;
				case ASTCENC_PRF_HDR:
					printf("    Color profile:              HDR RGBA\n");
					break;
			}

			if (block_z == 1)
			{
				printf("    Block size:                 %dx%d\n", block_x, block_y);
			}
			else
			{
				printf("    Block size:                 %dx%dx%d\n", block_x, block_y, block_z);
			}

			printf("    Bitrate:                    %3.2f bpp\n", 128.0 / (block_x * block_y * block_z));

			printf("    Radius mean/stdev:          %d texels\n", ewp.mean_stdev_radius);
			printf("    RGB power:                  %g\n", (double)ewp.rgb_power);
			printf("    RGB base weight:            %g\n", (double)ewp.rgb_base_weight);
			printf("    RGB mean weight:            %g\n", (double)ewp.rgb_mean_weight);
			printf("    RGB stdev weight:           %g\n", (double)ewp.rgb_stdev_weight);
			printf("    RGB mean/stdev mixing:      %g\n", (double)ewp.rgb_mean_and_stdev_mixing);
			printf("    Alpha power:                %g\n", (double)ewp.alpha_power);
			printf("    Alpha base weight:          %g\n", (double)ewp.alpha_base_weight);
			printf("    Alpha mean weight:          %g\n", (double)ewp.alpha_mean_weight);
			printf("    Alpha stdev weight:         %g\n", (double)ewp.alpha_stdev_weight);
			printf("    RGB alpha scale weight:     %d\n", ewp.enable_rgb_scale_with_alpha);
			if (ewp.enable_rgb_scale_with_alpha)
			{
				printf("    Radius RGB alpha scale:     %d texels\n", ewp.alpha_radius);
			}

			printf("    R channel weight:           %g\n",(double)ewp.rgba_weights[0]);
			printf("    G channel weight:           %g\n",(double)ewp.rgba_weights[1]);
			printf("    B channel weight:           %g\n",(double)ewp.rgba_weights[2]);
			printf("    A channel weight:           %g\n",(double)ewp.rgba_weights[3]);
			printf("    Deblock artifact setting:   %g\n", (double)ewp.block_artifact_suppression);
			printf("    Block partition cutoff:     %d partitions\n", ewp.partition_search_limit);
			printf("    PSNR cutoff:                %g dB\n", (double)dblimit_set);
			printf("    1->2 partition cutoff:      %g\n", (double)ewp.partition_1_to_2_limit);
			printf("    2 plane correlation cutoff: %g\n", (double)ewp.lowest_correlation_cutoff);
			printf("    Block mode centile cutoff:  %g%%\n", (double)(ewp.block_mode_cutoff * 100.0f));
			printf("    Max refinement cutoff:      %d iterations\n", ewp.max_refinement_iters);
			printf("    Compressor thread count:    %d\n", thread_count);
			printf("\n");
		}
	}

	int padding = MAX(ewp.mean_stdev_radius, ewp.alpha_radius);


	// Flatten out the list of operations we need to perform
	bool stage_compress = (op_mode == ASTC_ENCODE) || (op_mode == ASTC_ENCODE_AND_DECODE);
	bool stage_decompress = (op_mode == ASTC_DECODE) || (op_mode == ASTC_ENCODE_AND_DECODE);
	bool stage_load_uncomp = stage_compress;
	bool stage_load_comp = op_mode == ASTC_DECODE;
	bool stage_compare = op_mode == ASTC_ENCODE_AND_DECODE;
	bool stage_store_comp = op_mode == ASTC_ENCODE;
	bool stage_store_decomp = stage_decompress;

	astc_codec_image* image_uncomp_in = nullptr ;
	int image_uncomp_in_num_chan = 0;
	bool image_uncomp_in_is_hdr = false;

	astc_compressed_image image_comp;

	astc_codec_image* image_decomp_out = nullptr;

	// TODO: Handle RAII resources so they get freed when out of scope
	// Load the compressed input file if needed

	// This has to come first, as the block size is in the file header
	if (stage_load_comp)
	{
		int error = load_astc_file(input_filename, image_comp);
		if (error)
		{
			return 1;
		}
	}

	astcenc_error    codec_status;
	astcenc_config   codec_config;
	astcenc_context* codec_context;

	// TODO: Temporary code until we get command line parser populating config natively
	codec_config.profile = decode_mode;

	codec_config.flags = 0;

	if (ewp.ra_normal_angular_scale) {
		codec_config.flags |= ASTCENC_FLG_MAP_NORMAL;
	}

	if (ewp.enable_rgb_scale_with_alpha) {
		codec_config.flags |= ASTCENC_FLG_USE_ALPHA_WEIGHT;
	}

	if (linearize_srgb) {
		codec_config.flags |= ASTCENC_FLG_USE_LINEARIZED_SRGB;
	}

	if (stage_load_uncomp)
	{
		codec_config.block_x = block_x;
		codec_config.block_y = block_y;
		codec_config.block_z = block_z;
	}
	else
	{
		codec_config.block_x = image_comp.block_x;
		codec_config.block_y = image_comp.block_y;
		codec_config.block_z = image_comp.block_z;
	}

	codec_config.v_rgba_radius = ewp.mean_stdev_radius;
	codec_config.v_rgba_mean_stdev_mix = ewp.rgb_mean_and_stdev_mixing;
	codec_config.v_rgb_power = ewp.rgb_power;
	codec_config.v_rgb_base = ewp.rgb_base_weight;
	codec_config.v_rgb_mean = ewp.rgb_mean_weight;
	codec_config.v_rgb_stdev = ewp.rgb_stdev_weight;
	codec_config.v_a_power = ewp.alpha_power;
	codec_config.v_a_base = ewp.alpha_base_weight;
	codec_config.v_a_mean = ewp.alpha_mean_weight;
	codec_config.v_a_stdev = ewp.alpha_stdev_weight;
	codec_config.cw_r_weight = ewp.rgba_weights[0];
	codec_config.cw_g_weight = ewp.rgba_weights[1];
	codec_config.cw_b_weight = ewp.rgba_weights[2];
	codec_config.cw_a_weight = ewp.rgba_weights[3];
	codec_config.a_scale_radius = ewp.alpha_radius;
	codec_config.b_deblock_weight = ewp.block_artifact_suppression;
	codec_config.tune_partition_limit = ewp.partition_search_limit;
	codec_config.tune_block_mode_limit = ewp.block_mode_cutoff * 100.0f;
	codec_config.tune_refinement_limit = ewp.max_refinement_iters;
	// TODO: This isn't exactly right - we should move computing avg_error_limit in to the codec
	codec_config.tune_db_limit = ewp.texel_avg_error_limit;
	codec_config.tune_partition_early_out_limit = ewp.partition_1_to_2_limit;
	codec_config.tune_two_plane_early_out_limit = ewp.lowest_correlation_cutoff;

	codec_status = astcenc_context_alloc(codec_config, thread_count, &codec_context);
	if (codec_status != ASTCENC_SUCCESS) {
		printf("ERROR: Codec context alloc failed: %s\n", astcenc_get_error_string(codec_status));
		return 1;
	}

	// Load the uncompressed input file if needed
	if (stage_load_uncomp)
	{
		// TODO: This can be moved to command line parsing
		if ((array_size > 1) && (block_z == 1))
		{
			printf("ERROR: 3D input data for a 2D ASTC block format\n");
			return 1;
		}

		image_uncomp_in = load_uncomp_file(input_filename, array_size, padding, y_flip, linearize_srgb,
		                                   image_uncomp_in_is_hdr, image_uncomp_in_num_chan);
		if (!image_uncomp_in)
		{
			return 1;
		}

		if (!silentmode)
		{
			printf("Source image\n");
			printf("============\n\n");
			printf("    Source:                     %s\n", input_filename);
			printf("    Color profile:              %s\n", image_uncomp_in_is_hdr ? "HDR" : "LDR");
			if (image_uncomp_in->zsize > 1)
			{
				printf("    Dimensions:                 3D, %d x %d x %d\n",
				       image_uncomp_in->xsize, image_uncomp_in->ysize, image_uncomp_in->zsize);
			}
			else
			{
				printf("    Dimensions:                 2D, %d x %d\n",
				       image_uncomp_in->xsize, image_uncomp_in->ysize);
			}
			printf("    Channels:                   %d\n\n", image_uncomp_in_num_chan);
		}
	}

	start_coding_time = get_time();

	// Compress an image
	if (stage_compress)
	{
		// TODO: Loaded functions should return astcenc_image natively
		astcenc_image image;
		image.data8 = image_uncomp_in->data8;
		image.data16 = image_uncomp_in->data16;
		image.dim_x = image_uncomp_in->xsize;
		image.dim_y = image_uncomp_in->ysize;
		image.dim_z = image_uncomp_in->zsize;
		image.padding_texels = image_uncomp_in->padding;

		int xdims = image.dim_x;
		int ydims = image.dim_y;
		int zdims = image.dim_z;

		int xblocks = (xdims + block_x - 1) / block_x;
		int yblocks = (ydims + block_y - 1) / block_y;
		int zblocks = (zdims + block_z - 1) / block_z;

		size_t buffer_size = xblocks * yblocks * zblocks * 16;
		uint8_t* buffer = new uint8_t[buffer_size];

		codec_status = astcenc_compress_image(codec_context, image, swz_encode, buffer, buffer_size, 0);
		if (codec_status != ASTCENC_SUCCESS) {
			printf("ERROR: Codec compress failed: %s\n", astcenc_get_error_string(codec_status));
			return 1;
		}

		image_comp.block_x = block_x;
		image_comp.block_y = block_y;
		image_comp.block_z = block_z;
		image_comp.dim_x = xdims;
		image_comp.dim_y = ydims;
		image_comp.dim_z = zdims;
		image_comp.data = buffer;
	}

	// Decompress an image
	if (stage_decompress)
	{
		int out_bitness = get_output_filename_enforced_bitness(output_filename);
		if (out_bitness == -1)
		{
			bool is_hdr = (decode_mode == ASTCENC_PRF_HDR) || (decode_mode == ASTCENC_PRF_HDR_RGB_LDR_A);
			out_bitness = is_hdr ? 16 : 8;
		}

		image_decomp_out = alloc_image(
		    out_bitness, image_comp.dim_x, image_comp.dim_y, image_comp.dim_z, 0);

		// TODO: Standardized on astcenc_image everywhere in the CLI
		astcenc_image image;
		image.dim_x = image_comp.dim_x;
		image.dim_y = image_comp.dim_y;
		image.dim_z = image_comp.dim_z;
		image.data8 = image_decomp_out->data8;
		image.data16 = image_decomp_out->data16;
		image.padding_texels = image_decomp_out->padding;

		// TODO: Pass through data len to avoid out-of-bounds reads
		codec_status = astcenc_decompress_image(codec_context, image_comp.data, 0, image, swz_decode, 0);
		if (codec_status != ASTCENC_SUCCESS) {
			printf("ERROR: Codec decompress failed: %s\n", astcenc_get_error_string(codec_status));
			return 1;
		}
	}

	end_coding_time = get_time();

	// Print metrics in comparison mode
	if (stage_compare)
	{
		compute_error_metrics(image_uncomp_in_is_hdr, image_uncomp_in_num_chan, image_uncomp_in,
		                      image_decomp_out, low_fstop, high_fstop);
	}

	// Store compressed image
	if (stage_store_comp)
	{
		int error = store_astc_file(image_comp, output_filename);
		if (error)
		{
			return 1;
		}
	}

	// Store decompressed image
	if (stage_store_decomp)
	{
		int store_result = -1;
		const char *format_string = "";

		store_result = astc_codec_store_image(image_decomp_out, output_filename, &format_string, y_flip);
		if (store_result < 0)
		{
			printf("ERROR: Failed to write output image %s\n", output_filename);
			return 1;
		}
	}

	free_image(image_uncomp_in);
	free_image(image_decomp_out);
	astcenc_context_free(codec_context);
	delete[] image_comp.data;

	end_time = get_time();

	if (stage_compare || !silentmode)
	{
		printf("Coding time\n");
		printf("===========\n\n");
		printf("    Total time:                %6.2f s\n", end_time - start_time);
		printf("    Coding time:               %6.2f s\n", end_coding_time - start_coding_time);
	}

	return 0;
}
