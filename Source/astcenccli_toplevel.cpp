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
#include "astcenccli_internal.h"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <vector>

// TODO: Move these into astc_main().
static double start_time;
static double end_time;
static double start_coding_time;
static double end_coding_time;

static const uint32_t MAGIC_FILE_CONSTANT = 0x5CA1AB13;

struct astc_header
{
	uint8_t magic[4];
	uint8_t block_x;
	uint8_t block_y;
	uint8_t block_z;
	uint8_t dim_x[3];			// dims = dim[0] + (dim[1] << 8) + (dim[2] << 16)
	uint8_t dim_y[3];			// Sizes are given in texels;
	uint8_t dim_z[3];			// block count is inferred
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

enum astc_op_mode {
	ASTC_ENCODE,
	ASTC_DECODE,
	ASTC_ENCODE_AND_DECODE,
	ASTC_PRINT_LONGHELP,
	ASTC_PRINT_VERSION,
	ASTC_UNRECOGNIZED
};

int load_astc_file(
	const char* filename,
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

	int block_x = hdr.block_x;
	int block_y = hdr.block_y;
	int block_z = MAX(hdr.block_z, 1);

	if (((block_z == 1) && !is_legal_2d_block_size(block_x, block_y)) ||
	    ((block_z > 1) && !is_legal_3d_block_size(block_x, block_y, block_z)))
	{
		printf("ERROR: File corrupt '%s'\n", filename);
		return 1;
	}

	int dim_x = unpack_bytes(hdr.dim_x[0], hdr.dim_x[1], hdr.dim_x[2], 0);
	int dim_y = unpack_bytes(hdr.dim_y[0], hdr.dim_y[1], hdr.dim_y[2], 0);
	int dim_z = unpack_bytes(hdr.dim_z[0], hdr.dim_z[1], hdr.dim_z[2], 0);

	if (dim_x == 0 || dim_z == 0 || dim_z == 0)
	{
		printf("ERROR: File corrupt '%s'\n", filename);
		return 1;
	}

	int xblocks = (dim_x + block_x - 1) / block_x;
	int yblocks = (dim_y + block_y - 1) / block_y;
	int zblocks = (dim_z + block_z - 1) / block_z;

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
	out_image.dim_x = dim_x;
	out_image.dim_y = dim_y;
	out_image.dim_z = dim_z;
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
	hdr.block_x = comp_img.block_x;
	hdr.block_y = comp_img.block_y;
	hdr.block_z = comp_img.block_z;
	hdr.dim_x[0] =  comp_img.dim_x & 0xFF;
	hdr.dim_x[1] = (comp_img.dim_x >> 8) & 0xFF;
	hdr.dim_x[2] = (comp_img.dim_x >> 16) & 0xFF;
	hdr.dim_y[0] =  comp_img.dim_y & 0xFF;
	hdr.dim_y[1] = (comp_img.dim_y >> 8) & 0xFF;
	hdr.dim_y[2] = (comp_img.dim_y >> 16) & 0xFF;
	hdr.dim_z[0] =  comp_img.dim_z & 0xFF;
	hdr.dim_z[1] = (comp_img.dim_z >> 8) & 0xFF;
	hdr.dim_z[2] = (comp_img.dim_z >> 16) & 0xFF;

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
		if (!cpu_supports_sse42())
		{
			printf("CPU support error: host lacks SSE 4.2 support.\n");
			printf("    Please recompile with VEC=sse2.\n");
			exit(1);
		}
	#endif

	#if ASTC_POPCNT >= 1
		if (!cpu_supports_popcnt())
		{
			printf("CPU support error: host lacks POPCNT support.\n");
			printf("    Please recompile with VEC=sse2.\n");
			exit(1);
		}
	#endif

	#if ASTC_AVX >= 2
		if (!cpu_supports_avx2())
		{
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
static astcenc_image* load_uncomp_file(
	const char* filename,
	unsigned int dim_z,
	unsigned int dim_pad,
	bool y_flip,
	bool& is_hdr,
	unsigned int& num_components
) {
	astcenc_image *image = nullptr;

	// For a 2D image just load the image directly
	if (dim_z == 1)
	{
		image = astc_codec_load_image(filename, dim_pad, y_flip, is_hdr, num_components);
	}
	else
	{
		bool slice_is_hdr;
		unsigned int slice_num_components;
		astcenc_image* slice = nullptr;
		std::vector<astcenc_image*> slices;

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

			slice = astc_codec_load_image(slice_name.c_str(), dim_pad, y_flip,
			                              slice_is_hdr, slice_num_components);
			if (!slice)
			{
				break;
			}

			slices.push_back(slice);

			// Check it is not a 3D image
			if (slice->dim_z != 1)
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

				if ((slices[0]->dim_x != slice->dim_x) ||
				    (slices[0]->dim_y != slice->dim_y) ||
				    (slices[0]->dim_z != slice->dim_z))
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
			unsigned int dim_x = slices[0]->dim_x;
			unsigned int dim_y = slices[0]->dim_y;
			int bitness = is_hdr ? 16 : 8;
			int slice_size = (dim_x + (2 * dim_pad)) * (dim_y + (2 * dim_pad));

			image = alloc_image(bitness, dim_x, dim_y, dim_z, dim_pad);

			// Combine 2D source images into one 3D image; skipping padding slices
			for (unsigned int z = dim_pad; z < dim_z + dim_pad; z++)
			{
				if (bitness == 8)
				{
					size_t copy_size = slice_size * 4 * sizeof(uint8_t);
					memcpy(*image->data8[z], *slices[z - dim_pad]->data8[0], copy_size);
				}
				else
				{
					size_t copy_size = slice_size * 4 * sizeof(uint16_t);
					memcpy(*image->data16[z], *slices[z - dim_pad]->data16[0], copy_size);
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

/**
 * @brief Parse the command line and read op_mode, profile
 *        input and output file names
 *
 * @param argc
 * @param argv
 * @param op_mode              ASTC operation mode
 * @param profile              ASTC profile
 * @param input_filename       The input file name
 * @param output_filename      The output file name
 *
 * @return 0 if everything is Okay, 1 if there is some error
 */
int parse_commandline_options(
	int argc,
	char **argv,
	astc_op_mode& op_mode,
	astcenc_profile& profile,
	std::string& input_filename,
	std::string& output_filename
) {

	profile = ASTCENC_PRF_LDR;
	op_mode = ASTC_UNRECOGNIZED;

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
			profile = modes[i].decode_mode;
			break;
		}
	}

	if (op_mode == ASTC_UNRECOGNIZED)
	{
		printf("ERROR: Unrecognized operation '%s'\n", argv[1]);
		return 1;
	}

	input_filename = argc >= 3 ? argv[2] : "";
	output_filename = argc >= 4 ? argv[3] : "";

	if (input_filename.empty())
	{
		printf("ERROR: Input file not specified\n");
		return 1;
	}

	if (output_filename.empty())
	{
		printf("ERROR: Output file not specified\n");
		return 1;
	}

	return 0;
}

/**
 * @brief Initialize the astcenc_config
 *
 * @param argc
 * @param argv
 * @param op_mode          ASTC operation mode
 * @param comp_image
 * @param config           The astcenc configuration
 *
 * @return 0 if everything is Okay, 1 if there is some error
 */
int init_astcenc_config(
	int argc,
	char **argv,
	astcenc_profile profile,
	astc_op_mode op_mode,
	astc_compressed_image& comp_image,
	astcenc_config& config
) {
	unsigned int block_x = 0;
	unsigned int block_y = 0;
	unsigned int block_z = 1;

	// For decode the block size is set by the incoming image.
	if (op_mode == ASTC_DECODE)
	{
		block_x = comp_image.block_x;
		block_y = comp_image.block_y;
		block_z = comp_image.block_y;
	}

	astcenc_preset preset = ASTCENC_PRE_FAST;

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
		int dimensions = sscanf(argv[4], "%ux%u%nx%u%n", &block_x, &block_y, &cnt2D, &block_z, &cnt3D);
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

		if (!strcmp(argv[5], "-fast"))
		{
			preset = ASTCENC_PRE_FAST;
		}
		else if (!strcmp(argv[5], "-medium"))
		{
			preset = ASTCENC_PRE_MEDIUM;
		}
		else if (!strcmp(argv[5], "-thorough"))
		{
			preset = ASTCENC_PRE_THOROUGH;
		}
		else if (!strcmp(argv[5], "-exhaustive"))
		{
			preset = ASTCENC_PRE_EXHAUSTIVE;
		}
		else
		{
			printf("ERROR: Unknown search preset\n");
			return 1;
		}

		argidx = 6;
	}
	else
	{
		// For decode the block size and preset are not needed.
		argidx = 4;
	}

	unsigned int flags = 0;

	while (argidx < argc)
	{
		if (!strcmp(argv[argidx], "-a"))
		{
			argidx += 2;
			flags |= ASTCENC_FLG_USE_ALPHA_WEIGHT;
		}
		// presets begin here
		else if (!strcmp(argv[argidx], "-normal_psnr"))
		{
			argidx++;
			flags |= ASTCENC_FLG_MAP_NORMAL;
		}
		else if (!strcmp(argv[argidx], "-normal_percep"))
		{
			argidx++;
			flags |= ASTCENC_FLG_MAP_NORMAL;
			flags |= ASTCENC_FLG_USE_PERCEPTUAL;
		}
		else if (!strcmp(argv[argidx], "-mask"))
		{
			argidx++;
			flags |= ASTCENC_FLG_MAP_MASK;
		}
		else if (!strcmp(argv[argidx], "-linsrgb"))
		{
			argidx++;
			flags |= ASTCENC_FLG_USE_LINEARIZED_SRGB;
		}
		else
		{
			argidx ++;
		}
	}

	astcenc_error status = astcenc_init_config(profile, block_x, block_y, block_z, preset, flags, config);
	if (status != ASTCENC_SUCCESS)
	{
		printf ("ERROR: Failed to initialize configuration '%s'\n", astcenc_get_error_string(status));
		return 1;
	}

	return 0;
}

/**
 * @brief Edit the astcenc_config
 *
 * @param argc
 * @param argv
 * @param op_mode               ASTC operation mode
 * @param cli_config            Command line config options
 * @param config                The astcenc configuration
 *
 * @return 0 if everything is Okay, 1 if there is some error
 */
int edit_astcenc_config(
	int argc,
	char **argv,
	const astc_op_mode op_mode,
	cli_config_options& cli_config,
	astcenc_config& config
) {

	int argidx = (op_mode == ASTC_ENCODE || op_mode == ASTC_ENCODE_AND_DECODE) ? 6 : 4;

	while (argidx < argc)
	{
		if (!strcmp(argv[argidx], "-silent"))
		{
			argidx++;
			cli_config.silentmode = 1;
		}
		else
			if (!strcmp(argv[argidx], "-v"))
		{
			argidx += 7;
			if (argidx > argc)
			{
				printf("-v switch with less than 6 arguments, quitting\n");
				return 1;
			}

			config.v_rgba_radius = atoi(argv[argidx - 6]);
			config.v_rgb_power = static_cast<float>(atof(argv[argidx - 5]));
			config.v_rgb_base = static_cast<float>(atof(argv[argidx - 4]));
			config.v_rgb_mean = static_cast<float>(atof(argv[argidx - 3]));
			config.v_rgb_stdev = static_cast<float>(atof(argv[argidx - 2]));
			config.v_rgba_mean_stdev_mix = static_cast<float>(atof(argv[argidx - 1]));
		}
		else if (!strcmp(argv[argidx], "-va"))
		{
			argidx += 5;
			if (argidx > argc)
			{
				printf("-va switch with less than 4 arguments, quitting\n");
				return 1;
			}

			config.v_a_power= static_cast<float>(atof(argv[argidx - 4]));
			config.v_a_base = static_cast<float>(atof(argv[argidx - 3]));
			config.v_a_mean = static_cast<float>(atof(argv[argidx - 2]));
			config.v_a_stdev = static_cast<float>(atof(argv[argidx - 1]));
		}
		else if (!strcmp(argv[argidx], "-cw"))
		{
			argidx += 5;
			if (argidx > argc)
			{
				printf("-cw switch with less than 4 arguments\n");
				return 1;
			}

			config.cw_r_weight = static_cast<float>(atof(argv[argidx - 4]));
			config.cw_g_weight = static_cast<float>(atof(argv[argidx - 3]));
			config.cw_b_weight = static_cast<float>(atof(argv[argidx - 2]));
			config.cw_a_weight = static_cast<float>(atof(argv[argidx - 1]));
		}
		else if (!strcmp(argv[argidx], "-a"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-a switch with no argument\n");
				return 1;
			}

			config.a_scale_radius = atoi(argv[argidx - 1]);
		}
		else if (!strcmp(argv[argidx], "-b"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-b switch with no argument\n");
				return 1;
			}

			config.b_deblock_weight = static_cast<float>(atof(argv[argidx - 1]));
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

			cli_config.swz_encode.r = swizzle_components[0];
			cli_config.swz_encode.g = swizzle_components[1];
			cli_config.swz_encode.b = swizzle_components[2];
			cli_config.swz_encode.a = swizzle_components[3];
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

			cli_config.swz_decode.r = swizzle_components[0];
			cli_config.swz_decode.g = swizzle_components[1];
			cli_config.swz_decode.b = swizzle_components[2];
			cli_config.swz_decode.a = swizzle_components[3];
		}
		// presets begin here
		else if (!strcmp(argv[argidx], "-normal_psnr"))
		{
			argidx++;

			cli_config.swz_encode.r = ASTCENC_SWZ_R;
			cli_config.swz_encode.g = ASTCENC_SWZ_R;
			cli_config.swz_encode.b = ASTCENC_SWZ_R;
			cli_config.swz_encode.a = ASTCENC_SWZ_G;

			cli_config.swz_decode.r = ASTCENC_SWZ_R;
			cli_config.swz_decode.g = ASTCENC_SWZ_A;
			cli_config.swz_decode.b = ASTCENC_SWZ_Z;
			cli_config.swz_decode.a = ASTCENC_SWZ_1;
		}
		else if (!strcmp(argv[argidx], "-normal_percep"))
		{
			argidx++;

			cli_config.swz_encode.r = ASTCENC_SWZ_R;
			cli_config.swz_encode.g = ASTCENC_SWZ_R;
			cli_config.swz_encode.b = ASTCENC_SWZ_R;
			cli_config.swz_encode.a = ASTCENC_SWZ_G;

			cli_config.swz_decode.r = ASTCENC_SWZ_R;
			cli_config.swz_decode.g = ASTCENC_SWZ_A;
			cli_config.swz_decode.b = ASTCENC_SWZ_Z;
			cli_config.swz_decode.a = ASTCENC_SWZ_1;
		}
		else if (!strcmp(argv[argidx], "-mask"))
		{
			argidx++;
		}
		else if (!strcmp(argv[argidx], "-blockmodelimit"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-blockmodelimit switch with no argument\n");
				return 1;
			}

			config.tune_block_mode_limit = atoi(argv[argidx - 1]);
		}
		else if (!strcmp(argv[argidx], "-partitionlimit"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-partitionlimit switch with no argument\n");
				return 1;
			}

			config.tune_partition_limit = atoi(argv[argidx - 1]);
		}
		else if (!strcmp(argv[argidx], "-dblimit"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-dblimit switch with no argument\n");
				return 1;
			}

			if ((config.profile == ASTCENC_PRF_LDR) || (config.profile == ASTCENC_PRF_LDR_SRGB))
			{
				config.tune_db_limit = static_cast<float>(atof(argv[argidx - 1]));
			}
		}
		else if (!strcmp(argv[argidx], "-partitionearlylimit"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-partitionearlylimit switch with no argument\n");
				return 1;
			}

			config.tune_partition_early_out_limit = static_cast<float>(atof(argv[argidx - 1]));
		}
		else if (!strcmp(argv[argidx], "-planecorlimit"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-planecorlimit switch with no argument\n");
				return 1;
			}

			config.tune_two_plane_early_out_limit = static_cast<float>(atof(argv[argidx - 1]));
		}
		else if (!strcmp(argv[argidx], "-refinementlimit"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-refinementlimit switch with no argument\n");
				return 1;
			}

			config.tune_refinement_limit = atoi(argv[argidx - 1]);
		}
		else if (!strcmp(argv[argidx], "-j"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-j switch with no argument\n");
				return 1;
			}

			cli_config.thread_count = atoi(argv[argidx - 1]);
		}
		else if (!strcmp(argv[argidx], "-linsrgb"))
		{
			argidx++;
			cli_config.linearize_srgb = 1;
		}
		else if (!strcmp(argv[argidx], "-yflip"))
		{
			argidx++;
			cli_config.y_flip = 1;
		}
		else if (!strcmp(argv[argidx], "-mpsnr"))
		{
			argidx += 3;
			if (argidx > argc)
			{
				printf("-mpsnr switch with less than 2 arguments\n");
				return 1;
			}

			cli_config.low_fstop = atoi(argv[argidx - 2]);
			cli_config.high_fstop = atoi(argv[argidx - 1]);
			if (cli_config.high_fstop < cli_config.low_fstop)
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
			if (!sscanf(argv[argidx], "%d", &cli_config.array_size) || cli_config.array_size == 0)
			{
				printf("Invalid array size (image depth) given with -array option: \"%s\".\n", argv[argidx]);
				return 1;
			}

			if ((op_mode == ASTC_ENCODE || op_mode == ASTC_ENCODE_AND_DECODE)
					&& (cli_config.array_size > 1) && (config.block_z == 1))
			{
				printf("ERROR: 3D input data for a 2D ASTC block format\n");
				return 1;
			}
			argidx++;
		}
		else // check others as well
		{
			printf("Commandline argument \"%s\" not recognized\n", argv[argidx]);
			return 1;
		}
	}

	if (op_mode == ASTC_ENCODE || op_mode == ASTC_ENCODE_AND_DECODE)
	{
		if (cli_config.thread_count < 1)
		{
			cli_config.thread_count = get_cpu_count();
		}

		// print all encoding settings unless specifically told otherwise.
		if (!cli_config.silentmode)
		{
			printf("Compressor settings\n");
			printf("===================\n\n");

			switch(config.profile)
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

			if (config.block_z == 1)
			{
				printf("    Block size:                 %ux%u\n", config.block_x, config.block_y);
			}
			else
			{
				printf("    Block size:                 %ux%ux%u\n", config.block_x, config.block_y, config.block_z);
			}

			printf("    Bitrate:                    %3.2f bpp\n", 128.0 / (config.block_x * config.block_y * config.block_z));

			printf("    Radius mean/stdev:          %d texels\n", config.v_rgba_radius);
			printf("    RGB power:                  %g\n", (double)config.v_rgb_power );
			printf("    RGB base weight:            %g\n", (double)config.v_rgb_base);
			printf("    RGB mean weight:            %g\n", (double)config.v_rgb_mean);
			printf("    RGB stdev weight:           %g\n", (double)config.v_rgba_mean_stdev_mix);
			printf("    RGB mean/stdev mixing:      %g\n", (double)config.v_rgba_mean_stdev_mix);
			printf("    Alpha power:                %g\n", (double)config.v_a_power);
			printf("    Alpha base weight:          %g\n", (double)config.v_a_base);
			printf("    Alpha mean weight:          %g\n", (double)config.v_a_mean);
			printf("    Alpha stdev weight:         %g\n", (double)config.v_a_stdev);
			printf("    RGB alpha scale weight:     %d\n", (config.flags & ASTCENC_FLG_MAP_NORMAL));
			if ((config.flags & ASTCENC_FLG_MAP_NORMAL))
			{
				printf("    Radius RGB alpha scale:     %d texels\n", config.a_scale_radius);
			}

			printf("    R channel weight:           %g\n",(double)config.cw_r_weight);
			printf("    G channel weight:           %g\n",(double)config.cw_g_weight);
			printf("    B channel weight:           %g\n",(double)config.cw_b_weight);
			printf("    A channel weight:           %g\n",(double)config.cw_a_weight);
			printf("    Deblock artifact setting:   %g\n", (double)config.b_deblock_weight);
			printf("    Block partition cutoff:     %d partitions\n", config.tune_partition_limit);
			printf("    PSNR cutoff:                %g dB\n", (double)config.tune_db_limit);
			printf("    1->2 partition cutoff:      %g\n", (double)config.tune_partition_early_out_limit);
			printf("    2 plane correlation cutoff: %g\n", (double)config.tune_two_plane_early_out_limit);
			printf("    Block mode centile cutoff:  %g%%\n", (double)(config.tune_block_mode_limit));
			printf("    Max refinement cutoff:      %d iterations\n", config.tune_refinement_limit);
			printf("    Compressor thread count:    %d\n", cli_config.thread_count);
			printf("\n");
		}
	}

	return 0;
}

int main(
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

	astc_op_mode op_mode;
	astcenc_profile profile;

	std::string input_filename;
	std::string output_filename;

	int error = parse_commandline_options(argc, argv, op_mode, profile, input_filename, output_filename);
	if (error)
	{
		return 1;
	}

	switch (op_mode)
	{
		case ASTC_PRINT_LONGHELP:
			astcenc_print_longhelp();
			return 0;
		case ASTC_PRINT_VERSION:
			astcenc_print_header();
			return 0;
		default:
			break;
	}

	// TODO: Handle RAII resources so they get freed when out of scope
	// Load the compressed input file if needed

	// This has to come first, as the block size is in the file header
	astc_compressed_image image_comp;
	if (op_mode == ASTC_DECODE)
	{
		error = load_astc_file(input_filename.c_str(), image_comp);
		if (error)
		{
			return 1;
		}
	}

	astcenc_config config;

	error = init_astcenc_config(argc, argv, profile, op_mode, image_comp, config);
	if(error)
	{
		printf("ERROR: Astcenc configuration intialization failed.\n");
		return 1;
	}

	// Initialize cli_config_options with default values
	cli_config_options cli_config { 1, 0, 0, 0, 0, -10, 10,
		{ ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A },
		{ ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A } };

	error = edit_astcenc_config(argc, argv, op_mode, cli_config, config);
	if(error)
	{
		printf("ERROR: Astcenc configuration intialization failed.\n");
		return 1;
	}

	int padding = MAX(config.v_rgba_radius, config.a_scale_radius);

	// Flatten out the list of operations we need to perform
	bool stage_compress = (op_mode == ASTC_ENCODE) || (op_mode == ASTC_ENCODE_AND_DECODE);
	bool stage_decompress = (op_mode == ASTC_DECODE) || (op_mode == ASTC_ENCODE_AND_DECODE);
	bool stage_load_uncomp = stage_compress;
	bool stage_load_comp = op_mode == ASTC_DECODE;
	bool stage_compare = op_mode == ASTC_ENCODE_AND_DECODE;
	bool stage_store_comp = op_mode == ASTC_ENCODE;
	bool stage_store_decomp = stage_decompress;

	astcenc_image* image_uncomp_in = nullptr ;
	unsigned int image_uncomp_in_num_chan = 0;
	bool image_uncomp_in_is_hdr = false;

	astcenc_image* image_decomp_out = nullptr;

	// TODO: Handle RAII resources so they get freed when out of scope
	// Load the compressed input file if needed

	// This has to come first, as the block size is in the file header
	if (stage_load_comp)
	{
		error = load_astc_file(input_filename.c_str(), image_comp);
		if (error)
		{
			return 1;
		}
	}

	astcenc_error    codec_status;
	astcenc_context* codec_context;

	if (!stage_load_uncomp)
	{
		config.block_x = image_comp.block_x;
		config.block_y = image_comp.block_y;
		config.block_z = image_comp.block_z;
	}

	codec_status = astcenc_context_alloc(config, cli_config.thread_count, &codec_context);

	if (codec_status != ASTCENC_SUCCESS)
	{
		printf("ERROR: Codec context alloc failed: %s\n", astcenc_get_error_string(codec_status));
		return 1;
	}

	// Load the uncompressed input file if needed
	if (stage_load_uncomp)
	{
		image_uncomp_in = load_uncomp_file(input_filename.c_str(), cli_config.array_size, padding,
		                                   cli_config.y_flip, image_uncomp_in_is_hdr, image_uncomp_in_num_chan);
		if (!image_uncomp_in)
		{
			printf ("ERROR: Failed to load uncompressed image file\n");
			return 1;
		}

		if (!cli_config.silentmode)
		{
			printf("Source image\n");
			printf("============\n\n");
			printf("    Source:                     %s\n", input_filename.c_str());
			printf("    Color profile:              %s\n", image_uncomp_in_is_hdr ? "HDR" : "LDR");
			if (image_uncomp_in->dim_z > 1)
			{
				printf("    Dimensions:                 3D, %ux%ux%u\n",
				       image_uncomp_in->dim_x, image_uncomp_in->dim_y, image_uncomp_in->dim_z);
			}
			else
			{
				printf("    Dimensions:                 2D, %ux%u\n",
				       image_uncomp_in->dim_x, image_uncomp_in->dim_y);
			}
			printf("    Channels:                   %d\n\n", image_uncomp_in_num_chan);
		}
	}

	start_coding_time = get_time();

	// Compress an image
	if (stage_compress)
	{
		unsigned int blocks_x = (image_uncomp_in->dim_x + config.block_x - 1) / config.block_x;
		unsigned int blocks_y = (image_uncomp_in->dim_y + config.block_y - 1) / config.block_y;
		unsigned int blocks_z = (image_uncomp_in->dim_z + config.block_z - 1) / config.block_z;
		size_t buffer_size = blocks_x * blocks_y * blocks_z * 16;
		uint8_t* buffer = new uint8_t[buffer_size];

		codec_status = astcenc_compress_image(codec_context, *image_uncomp_in, cli_config.swz_encode,
		                                      buffer, buffer_size, 0);
		if (codec_status != ASTCENC_SUCCESS)
		{
			printf("ERROR: Codec compress failed: %s\n", astcenc_get_error_string(codec_status));
			return 1;
		}

		image_comp.block_x = config.block_x;
		image_comp.block_y = config.block_y;
		image_comp.block_z = config.block_z;
		image_comp.dim_x = image_uncomp_in->dim_x;
		image_comp.dim_y = image_uncomp_in->dim_y;
		image_comp.dim_z = image_uncomp_in->dim_z;
		image_comp.data = buffer;
	}

	// Decompress an image
	if (stage_decompress)
	{
		int out_bitness = get_output_filename_enforced_bitness(output_filename.c_str());
		if (out_bitness == -1)
		{
			bool is_hdr = (config.profile == ASTCENC_PRF_HDR) || (config.profile == ASTCENC_PRF_HDR_RGB_LDR_A);
			out_bitness = is_hdr ? 16 : 8;
		}

		image_decomp_out = alloc_image(
		    out_bitness, image_comp.dim_x, image_comp.dim_y, image_comp.dim_z, 0);

		// TODO: Pass through data len to avoid out-of-bounds reads
		codec_status = astcenc_decompress_image(codec_context, image_comp.data, 0,
		                                        *image_decomp_out, cli_config.swz_decode, 0);
		if (codec_status != ASTCENC_SUCCESS)
		{
			printf("ERROR: Codec decompress failed: %s\n", astcenc_get_error_string(codec_status));
			return 1;
		}
	}

	end_coding_time = get_time();

	// Print metrics in comparison mode
	if (stage_compare)
	{
		compute_error_metrics(image_uncomp_in_is_hdr, image_uncomp_in_num_chan, image_uncomp_in,
		                      image_decomp_out, cli_config.low_fstop, cli_config.high_fstop);
	}

	// Store compressed image
	if (stage_store_comp)
	{
		error = store_astc_file(image_comp, output_filename.c_str());
		if (error)
		{
			printf ("ERROR: Failed to store compressed image\n");
			return 1;
		}
	}

	// Store decompressed image
	if (stage_store_decomp)
	{
		int store_result = -1;
		const char *format_string = "";

		store_result = astc_codec_store_image(image_decomp_out, output_filename.c_str(),
		                                      &format_string, cli_config.y_flip);
		if (store_result < 0)
		{
			printf("ERROR: Failed to write output image %s\n", output_filename.c_str());
			return 1;
		}
	}

	free_image(image_uncomp_in);
	free_image(image_decomp_out);
	astcenc_context_free(codec_context);

	delete[] image_comp.data;

	end_time = get_time();

	if (stage_compare || !cli_config.silentmode)
	{
		printf("Coding time\n");
		printf("===========\n\n");
		printf("    Total time:                %6.2f s\n", end_time - start_time);
		printf("    Coding time:               %6.2f s\n", end_coding_time - start_coding_time);
	}

	return 0;
}
