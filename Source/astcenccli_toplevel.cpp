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
#include <cstring>
#include <string>
#include <vector>

/* ============================================================================
	Data structure definitions
============================================================================ */

typedef unsigned int astcenc_operation;

struct mode_entry {
	std::string opt;
	astcenc_operation operation;
	astcenc_profile decode_mode;
};

/* ============================================================================
	Constants and literals
============================================================================ */

static const unsigned int ASTCENC_STAGE_LD_COMP    = 1 << 0;
static const unsigned int ASTCENC_STAGE_ST_COMP    = 1 << 1;

static const unsigned int ASTCENC_STAGE_LD_NCOMP   = 1 << 2;
static const unsigned int ASTCENC_STAGE_ST_NCOMP   = 1 << 3;

static const unsigned int ASTCENC_STAGE_COMPRESS   = 1 << 4;
static const unsigned int ASTCENC_STAGE_DECOMPRESS = 1 << 5;
static const unsigned int ASTCENC_STAGE_COMPARE    = 1 << 6;

static const astcenc_operation ASTCENC_OP_UNKNOWN  = 0;
static const astcenc_operation ASTCENC_OP_HELP     = 1 << 7;
static const astcenc_operation ASTCENC_OP_VERSION  = 1 << 8;

static const astcenc_operation ASTCENC_OP_COMPRESS =
                               ASTCENC_STAGE_LD_NCOMP |
                               ASTCENC_STAGE_COMPRESS |
                               ASTCENC_STAGE_ST_COMP;

static const astcenc_operation ASTCENC_OP_DECOMPRESS =
                               ASTCENC_STAGE_LD_COMP |
                               ASTCENC_STAGE_DECOMPRESS |
                               ASTCENC_STAGE_ST_NCOMP;

static const astcenc_operation ASTCENC_OP_TEST =
                               ASTCENC_STAGE_LD_NCOMP |
                               ASTCENC_STAGE_COMPRESS |
                               ASTCENC_STAGE_DECOMPRESS |
                               ASTCENC_STAGE_COMPARE |
                               ASTCENC_STAGE_ST_NCOMP;

static const mode_entry modes[] = {
	{"-cl",      ASTCENC_OP_COMPRESS,   ASTCENC_PRF_LDR},
	{"-dl",      ASTCENC_OP_DECOMPRESS, ASTCENC_PRF_LDR},
	{"-tl",      ASTCENC_OP_TEST,       ASTCENC_PRF_LDR},
	{"-cs",      ASTCENC_OP_COMPRESS,   ASTCENC_PRF_LDR_SRGB},
	{"-ds",      ASTCENC_OP_DECOMPRESS, ASTCENC_PRF_LDR_SRGB},
	{"-ts",      ASTCENC_OP_TEST,       ASTCENC_PRF_LDR_SRGB},
	{"-ch",      ASTCENC_OP_COMPRESS,   ASTCENC_PRF_HDR},
	{"-dh",      ASTCENC_OP_DECOMPRESS, ASTCENC_PRF_HDR},
	{"-th",      ASTCENC_OP_TEST,       ASTCENC_PRF_HDR},
	{"-cH",      ASTCENC_OP_COMPRESS,   ASTCENC_PRF_HDR_RGB_LDR_A},
	{"-dH",      ASTCENC_OP_DECOMPRESS, ASTCENC_PRF_HDR_RGB_LDR_A},
	{"-tH",      ASTCENC_OP_TEST,       ASTCENC_PRF_HDR_RGB_LDR_A},
	{"-h",       ASTCENC_OP_HELP,       ASTCENC_PRF_HDR},
	{"-help",    ASTCENC_OP_HELP,       ASTCENC_PRF_HDR},
	{"-v",       ASTCENC_OP_VERSION,    ASTCENC_PRF_HDR},
	{"-version", ASTCENC_OP_VERSION,    ASTCENC_PRF_HDR}
};

struct compression_workload {
	astcenc_context* context;
	astcenc_image* image;
	astcenc_swizzle swizzle;
	uint8_t* data_out;
	size_t data_len;
	astcenc_error error;
};

static bool ends_with(const std::string& str, const std::string& suffix)
{
	return (str.size() >= suffix.size()) &&
	       (0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix));
}

static void compression_workload_runner(
	int thread_count,
	int thread_id,
	void* payload
) {
	(void)thread_count;

	compression_workload* work = static_cast<compression_workload*>(payload);
	astcenc_error error = astcenc_compress_image(
	                       work->context, *work->image, work->swizzle,
	                       work->data_out, work->data_len, thread_id);

	// This is a racy update, so which error gets returned is a random, but it
	// will reliably report an error if an error occurs
	if (error != ASTCENC_SUCCESS) {
		work->error = error;
	}
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
	size_t sep = basename.find_last_of(".");
	if (sep == std::string::npos)
	{
		error = true;
		return "";
	}

	std::string base = basename.substr(0, sep);
	std::string ext = basename.substr(sep);
	std::string name = base + "_" + std::to_string(index) + ext;
	error = false;
	return name;
}

/**
 * @brief Load a non-astc image file from memory.
 *
 * @param filename            The file to load, or a pattern for array loads.
 * @param dim_z               The number of slices to load.
 * @param padding             The number of texels of padding.
 * @param y_flip              Should this image be Y flipped?
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
		image = load_ncimage(filename, dim_pad, y_flip, is_hdr, num_components);
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
				printf("ERROR: Image pattern does not contain file extension: %s\n", filename);
				break;
			}

			slice = load_ncimage(slice_name.c_str(), dim_pad, y_flip,
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
				if (image->data_type == ASTCENC_TYPE_U8)
				{
					uint8_t*** data8 = static_cast<uint8_t***>(image->data);
					uint8_t*** data8src = static_cast<uint8_t***>(slices[z - dim_pad]->data);
					size_t copy_size = slice_size * 4 * sizeof(uint8_t);
					memcpy(*data8[z], *data8src[0], copy_size);
				}
				else // if (image->data_type == ASTCENC_TYPE_F16)
				{
					assert(image->data_type == ASTCENC_TYPE_F16);
					uint16_t*** data16 = static_cast<uint16_t***>(image->data);
					uint16_t*** data16src = static_cast<uint16_t***>(slices[z - dim_pad]->data);
					size_t copy_size = slice_size * 4 * sizeof(uint16_t);
					memcpy(*data16[z], *data16src[0], copy_size);
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
 * @brief Parse the command line and read operation, profile
 *        input and output file names
 *
 * @param argc
 * @param argv
 * @param operation            ASTC operation mode
 * @param profile              ASTC profile
 *
 * @return 0 if everything is Okay, 1 if there is some error
 */
int parse_commandline_options(
	int argc,
	char **argv,
	astcenc_operation& operation,
	astcenc_profile& profile
) {
	assert(argc >= 2); (void)argc;

	profile = ASTCENC_PRF_LDR;
	operation = ASTCENC_OP_UNKNOWN;

	int modes_count = sizeof(modes) / sizeof(modes[0]);
	for (int i = 0; i < modes_count; i++)
	{
		if (modes[i].opt == argv[1])
		{
			operation = modes[i].operation;
			profile = modes[i].decode_mode;
			break;
		}
	}

	if (operation == ASTCENC_OP_UNKNOWN)
	{
		printf("ERROR: Unrecognized operation '%s'\n", argv[1]);
		return 1;
	}

	return 0;
}

/**
 * @brief Initialize the astcenc_config
 *
 * @param argc
 * @param argv
 * @param operation          ASTC operation mode
 * @param comp_image
 * @param config           The astcenc configuration
 *
 * @return 0 if everything is Okay, 1 if there is some error
 */
int init_astcenc_config(
	int argc,
	char **argv,
	astcenc_profile profile,
	astcenc_operation operation,
	astc_compressed_image& comp_image,
	astcenc_config& config
) {
	unsigned int block_x = 0;
	unsigned int block_y = 0;
	unsigned int block_z = 1;

	// For decode the block size is set by the incoming image.
	if (operation == ASTCENC_OP_DECOMPRESS)
	{
		block_x = comp_image.block_x;
		block_y = comp_image.block_y;
		block_z = comp_image.block_z;
	}

	astcenc_preset preset = ASTCENC_PRE_FAST;

	// parse the command line's encoding options.
	int argidx = 4;
	if (operation & ASTCENC_STAGE_COMPRESS)
	{
		// Read and decode block size
		if (argc < 5)
		{
			printf("ERROR: Block size must be specified\n");
			return 1;
		}

		int cnt2D, cnt3D;
		int dimensions = sscanf(argv[4], "%ux%u%nx%u%n",
		                        &block_x, &block_y, &cnt2D, &block_z, &cnt3D);
		// Character after the last match should be a NUL
		if (!(((dimensions == 2) && !argv[4][cnt2D]) || ((dimensions == 3) && !argv[4][cnt3D])))
		{
			printf("ERROR: Block size '%s' is invalid\n", argv[4]);
			return 1;
		}

		// Read and decode search preset
		if (argc < 6)
		{
			printf("ERROR: Search preset must be specified\n");
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
			printf("ERROR: Search preset '%s' is invalid\n", argv[5]);
			return 1;
		}

		argidx = 6;
	}

	unsigned int flags = 0;

	// Gather the flags that we need
	while (argidx < argc)
	{
		if (!strcmp(argv[argidx], "-a"))
		{
			// Skip over the data value for now
			argidx++;
			flags |= ASTCENC_FLG_USE_ALPHA_WEIGHT;
		}
		else if (!strcmp(argv[argidx], "-normal"))
		{
			flags |= ASTCENC_FLG_MAP_NORMAL;
		}
		else if (!strcmp(argv[argidx], "-perceptual"))
		{
			flags |= ASTCENC_FLG_USE_PERCEPTUAL;
		}
		else if (!strcmp(argv[argidx], "-mask"))
		{

			flags |= ASTCENC_FLG_MAP_MASK;
		}
		argidx ++;
	}

#if defined(ASTCENC_DECOMPRESS_ONLY)
	flags |= ASTCENC_FLG_DECOMPRESS_ONLY;
#endif

	astcenc_error status = astcenc_config_init(profile, block_x, block_y, block_z, preset, flags, config);
	if (status == ASTCENC_ERR_BAD_BLOCK_SIZE)
	{
		printf("ERROR: Block size '%s' is invalid\n", argv[4]);
		return 1;
	}
	else if (status == ASTCENC_ERR_BAD_CPU_ISA)
	{
		printf("ERROR: Required SIMD ISA support missing on this CPU\n");
		return 1;
	}
	else if (status == ASTCENC_ERR_BAD_CPU_FLOAT)
	{
		printf("ERROR: astcenc must not be compiled with -ffast-math\n");
		return 1;
	}
	else if (status != ASTCENC_SUCCESS)
	{
		printf("ERROR: Init config failed with %s\n", astcenc_get_error_string(status));
		return 1;
	}

	return 0;
}

/**
 * @brief Edit the astcenc_config
 *
 * @param argc
 * @param argv
 * @param operation               ASTC operation mode
 * @param cli_config            Command line config options
 * @param config                The astcenc configuration
 *
 * @return 0 if everything is Okay, 1 if there is some error
 */
int edit_astcenc_config(
	int argc,
	char **argv,
	const astcenc_operation operation,
	cli_config_options& cli_config,
	astcenc_config& config
) {

	int argidx = (operation & ASTCENC_STAGE_COMPRESS) ? 6 : 4;

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
				printf("ERROR: -v switch with less than 6 arguments\n");
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
				printf("ERROR: -va switch with less than 4 arguments\n");
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
				printf("ERROR: -cw switch with less than 4 arguments\n");
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
				printf("ERROR: -a switch with no argument\n");
				return 1;
			}

			config.a_scale_radius = atoi(argv[argidx - 1]);
		}
		else if (!strcmp(argv[argidx], "-b"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("ERROR: -b switch with no argument\n");
				return 1;
			}

			config.b_deblock_weight = static_cast<float>(atof(argv[argidx - 1]));
		}
		else if (!strcmp(argv[argidx], "-esw"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("ERROR: -esw switch with no argument\n");
				return 1;
			}

			if (strlen(argv[argidx - 1]) != 4)
			{
				printf("ERROR: -esw pattern does not contain 4 characters\n");
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
					printf("ERROR: -esw channel '%c' is not valid\n", argv[argidx - 1][i]);
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
				printf("ERROR: -dsw switch with no argument\n");
				return 1;
			}

			if (strlen(argv[argidx - 1]) != 4)
			{
				printf("ERROR: -dsw switch does not contain 4 characters\n");
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
					printf("ERROR: ERROR: -dsw channel '%c' is not valid\n", argv[argidx - 1][i]);
					return 1;
				}
			}

			cli_config.swz_decode.r = swizzle_components[0];
			cli_config.swz_decode.g = swizzle_components[1];
			cli_config.swz_decode.b = swizzle_components[2];
			cli_config.swz_decode.a = swizzle_components[3];
		}
		// presets begin here
		else if (!strcmp(argv[argidx], "-normal"))
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
		else if (!strcmp(argv[argidx], "-perceptual"))
		{
			argidx++;
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
				printf("ERROR: -blockmodelimit switch with no argument\n");
				return 1;
			}

			config.tune_block_mode_limit = atoi(argv[argidx - 1]);
		}
		else if (!strcmp(argv[argidx], "-partitionlimit"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("ERROR: -partitionlimit switch with no argument\n");
				return 1;
			}

			config.tune_partition_limit = atoi(argv[argidx - 1]);
		}
		else if (!strcmp(argv[argidx], "-dblimit"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("ERROR: -dblimit switch with no argument\n");
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
				printf("ERROR: -partitionearlylimit switch with no argument\n");
				return 1;
			}

			config.tune_partition_early_out_limit = static_cast<float>(atof(argv[argidx - 1]));
		}
		else if (!strcmp(argv[argidx], "-planecorlimit"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("ERROR: -planecorlimit switch with no argument\n");
				return 1;
			}

			config.tune_two_plane_early_out_limit = static_cast<float>(atof(argv[argidx - 1]));
		}
		else if (!strcmp(argv[argidx], "-refinementlimit"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("ERROR: -refinementlimit switch with no argument\n");
				return 1;
			}

			config.tune_refinement_limit = atoi(argv[argidx - 1]);
		}
		else if (!strcmp(argv[argidx], "-j"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("ERROR: -j switch with no argument\n");
				return 1;
			}

			cli_config.thread_count = atoi(argv[argidx - 1]);
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
				printf("ERROR: -mpsnr switch with less than 2 arguments\n");
				return 1;
			}

			cli_config.low_fstop = atoi(argv[argidx - 2]);
			cli_config.high_fstop = atoi(argv[argidx - 1]);
			if (cli_config.high_fstop < cli_config.low_fstop)
			{
				printf("ERROR: -mpsnr switch <low> is greater than the <high>\n");
				return 1;
			}
		}
		// Option: Encode a 3D image from an array of 2D images.
		else if (!strcmp(argv[argidx], "-array"))
		{
			// Only supports compressing
			if (!(operation & ASTCENC_STAGE_COMPRESS))
			{
				printf("ERROR: -array switch is only valid for compression\n");
				return 1;
			}

			// Image depth must be specified.
			if (argidx + 2 > argc)
			{
				printf("ERROR: -array switch with no argument\n");
				return 1;
			}
			argidx++;

			// Read array size (image depth).
			if (!sscanf(argv[argidx], "%u", &cli_config.array_size) || cli_config.array_size == 0)
			{
				printf("ERROR: -array size '%s' is invalid\n", argv[argidx]);
				return 1;
			}

			if ((cli_config.array_size > 1) && (config.block_z == 1))
			{
				printf("ERROR: -array with 3D input data for a 2D output format\n");
				return 1;
			}
			argidx++;
		}
		else // check others as well
		{
			printf("ERROR: Argument '%s' not recognized\n", argv[argidx]);
			return 1;
		}
	}

	if (cli_config.thread_count <= 0)
	{
		cli_config.thread_count = get_cpu_count();
	}

#if defined(ASTCENC_DECOMPRESS_ONLY)
	cli_config.thread_count = 1;
#endif

	if (operation & ASTCENC_STAGE_COMPRESS)
	{
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

			printf("    Radius mean/stdev:          %u texels\n", config.v_rgba_radius);
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
				printf("    Radius RGB alpha scale:     %u texels\n", config.a_scale_radius);
			}

			printf("    R channel weight:           %g\n",(double)config.cw_r_weight);
			printf("    G channel weight:           %g\n",(double)config.cw_g_weight);
			printf("    B channel weight:           %g\n",(double)config.cw_b_weight);
			printf("    A channel weight:           %g\n",(double)config.cw_a_weight);
			printf("    Deblock artifact setting:   %g\n", (double)config.b_deblock_weight);
			printf("    Block partition cutoff:     %u partitions\n", config.tune_partition_limit);
			printf("    PSNR cutoff:                %g dB\n", (double)config.tune_db_limit);
			printf("    1->2 partition cutoff:      %g\n", (double)config.tune_partition_early_out_limit);
			printf("    2 plane correlation cutoff: %g\n", (double)config.tune_two_plane_early_out_limit);
			printf("    Block mode centile cutoff:  %g%%\n", (double)(config.tune_block_mode_limit));
			printf("    Max refinement cutoff:      %u iterations\n", config.tune_refinement_limit);
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
	double start_time = get_time();

	#ifdef DEBUG_CAPTURE_NAN
		feenableexcept(FE_DIVBYZERO | FE_INVALID);
	#endif

	if (argc < 2)
	{
		astcenc_print_shorthelp();
		return 0;
	}

	astcenc_operation operation;
	astcenc_profile profile;
	int error = parse_commandline_options(argc, argv, operation, profile);
	if (error)
	{
		return 1;
	}

	switch (operation)
	{
	case ASTCENC_OP_HELP:
		astcenc_print_longhelp();
		return 0;
	case ASTCENC_OP_VERSION:
		astcenc_print_header();
		return 0;
	default:
		break;
	}


	std::string input_filename = argc >= 3 ? argv[2] : "";
	std::string output_filename = argc >= 4 ? argv[3] : "";

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

	// TODO: Handle RAII resources so they get freed when out of scope
	// Load the compressed input file if needed

	// This has to come first, as the block size is in the file header
	astc_compressed_image image_comp {};
	if (operation & ASTCENC_STAGE_LD_COMP)
	{
		if (ends_with(input_filename, ".astc"))
		{
			// TODO: Just pass on a std::string
			error = load_cimage(input_filename.c_str(), image_comp);
			if (error)
			{
				return 1;
			}
		}
		else if (ends_with(input_filename, ".ktx"))
		{
			// TODO: Just pass on a std::string
			bool is_srgb;
			error = load_ktx_compressed_image(input_filename.c_str(), is_srgb, image_comp);
			if (error)
			{
				return 1;
			}

			if (is_srgb && (profile != ASTCENC_PRF_LDR_SRGB))
			{
				printf("WARNING: Input file is sRGB, but decompressing as linear\n");
			}

			if (!is_srgb && (profile == ASTCENC_PRF_LDR_SRGB))
			{
				printf("WARNING: Input file is linear, but decompressing as sRGB\n");
			}
		}
		else
		{
			printf("ERROR: Unknown compressed input file type\n");
		}
	}

	astcenc_config config {};
	error = init_astcenc_config(argc, argv, profile, operation, image_comp, config);
	if (error)
	{
		return 1;
	}

	// Initialize cli_config_options with default values
	cli_config_options cli_config { 0, 1, false, false, -10, 10,
		{ ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A },
		{ ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A } };

	error = edit_astcenc_config(argc, argv, operation, cli_config, config);
	if (error)
	{
		return 1;
	}

	int padding = MAX(config.v_rgba_radius, config.a_scale_radius);
	astcenc_image* image_uncomp_in = nullptr ;
	unsigned int image_uncomp_in_num_chan = 0;
	bool image_uncomp_in_is_hdr = false;
	astcenc_image* image_decomp_out = nullptr;

	// TODO: Handle RAII resources so they get freed when out of scope
	astcenc_error    codec_status;
	astcenc_context* codec_context;

	codec_status = astcenc_context_alloc(config, cli_config.thread_count, &codec_context);
	if (codec_status != ASTCENC_SUCCESS)
	{
		printf("ERROR: Codec context alloc failed: %s\n", astcenc_get_error_string(codec_status));
		return 1;
	}

	// Load the uncompressed input file if needed
	if (operation & ASTCENC_STAGE_LD_NCOMP)
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

	double start_coding_time = get_time();

	// Compress an image
	if (operation & ASTCENC_STAGE_COMPRESS)
	{
		unsigned int blocks_x = (image_uncomp_in->dim_x + config.block_x - 1) / config.block_x;
		unsigned int blocks_y = (image_uncomp_in->dim_y + config.block_y - 1) / config.block_y;
		unsigned int blocks_z = (image_uncomp_in->dim_z + config.block_z - 1) / config.block_z;
		size_t buffer_size = blocks_x * blocks_y * blocks_z * 16;
		uint8_t* buffer = new uint8_t[buffer_size];

		compression_workload work;
		work.context = codec_context;
		work.image = image_uncomp_in;
		work.swizzle = cli_config.swz_encode;
		work.data_out = buffer;
		work.data_len = buffer_size;
		work.error = ASTCENC_SUCCESS;

		// Only launch worker threads for multi-threaded use - it makes basic
		// single-threaded profiling and debugging a little less convoluted
		if (cli_config.thread_count > 1)
		{
			launch_threads(cli_config.thread_count, compression_workload_runner, &work);
		}
		else
		{
			work.error = astcenc_compress_image(
			    work.context, *work.image, work.swizzle,
			    work.data_out, work.data_len, 0);
		}

		if (work.error != ASTCENC_SUCCESS)
		{
			printf("ERROR: Codec compress failed: %s\n", astcenc_get_error_string(work.error));
			return 1;
		}

		image_comp.block_x = config.block_x;
		image_comp.block_y = config.block_y;
		image_comp.block_z = config.block_z;
		image_comp.dim_x = image_uncomp_in->dim_x;
		image_comp.dim_y = image_uncomp_in->dim_y;
		image_comp.dim_z = image_uncomp_in->dim_z;
		image_comp.data = buffer;
		image_comp.data_len = buffer_size;
	}

	// Decompress an image
	if (operation & ASTCENC_STAGE_DECOMPRESS)
	{
		int out_bitness = get_output_filename_enforced_bitness(output_filename.c_str());
		if (out_bitness == -1)
		{
			bool is_hdr = (config.profile == ASTCENC_PRF_HDR) || (config.profile == ASTCENC_PRF_HDR_RGB_LDR_A);
			out_bitness = is_hdr ? 16 : 8;
		}

		image_decomp_out = alloc_image(
		    out_bitness, image_comp.dim_x, image_comp.dim_y, image_comp.dim_z, 0);

		codec_status = astcenc_decompress_image(codec_context, image_comp.data, image_comp.data_len,
		                                        *image_decomp_out, cli_config.swz_decode);
		if (codec_status != ASTCENC_SUCCESS)
		{
			printf("ERROR: Codec decompress failed: %s\n", astcenc_get_error_string(codec_status));
			return 1;
		}
	}

	double end_coding_time = get_time();

	// Print metrics in comparison mode
	if (operation & ASTCENC_STAGE_COMPARE)
	{
		compute_error_metrics(image_uncomp_in_is_hdr, image_uncomp_in_num_chan, image_uncomp_in,
		                      image_decomp_out, cli_config.low_fstop, cli_config.high_fstop);
	}

	// Store compressed image
	if (operation & ASTCENC_STAGE_ST_COMP)
	{
		if (ends_with(output_filename, ".astc"))
		{
			error = store_cimage(image_comp, output_filename.c_str());
			if (error)
			{
				printf ("ERROR: Failed to store compressed image\n");
				return 1;
			}
		}
		else if (ends_with(output_filename, ".ktx"))
		{
			bool srgb = profile == ASTCENC_PRF_LDR_SRGB;
			error = store_ktx_compressed_image(image_comp, output_filename.c_str(), srgb);
			if (error)
			{
				printf ("ERROR: Failed to store compressed image\n");
				return 1;
			}
		}
		else
		{
			printf("ERROR: Unknown compressed output file type\n");
			return 1;
		}
	}

	// Store decompressed image
	if (operation & ASTCENC_STAGE_ST_NCOMP)
	{
		int store_result = -1;
		const char *format_string = "";

		store_result = store_ncimage(image_decomp_out, output_filename.c_str(),
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

	double end_time = get_time();

	if ((operation & ASTCENC_STAGE_COMPARE) || (!cli_config.silentmode))
	{
		printf("Coding time\n");
		printf("===========\n\n");
		printf("    Total time:                %6.2f s\n", end_time - start_time);
		printf("    Coding time:               %6.2f s\n", end_coding_time - start_coding_time);
	}

	return 0;
}
