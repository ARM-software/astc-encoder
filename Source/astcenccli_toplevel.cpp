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

#ifdef DEBUG_CAPTURE_NAN
	#ifndef _GNU_SOURCE
		#define _GNU_SOURCE
	#endif

	#include <fenv.h>
#endif

#ifdef DEBUG_PRINT_DIAGNOSTICS
	int print_diagnostics = 0;
	int diagnostics_tile = -1;
	int print_tile_errors = 0;
	int print_statistics = 0;
#endif

static double start_time;
static double end_time;
static double start_coding_time;
static double end_coding_time;


NORETURN __attribute__((visibility("default")))
void astc_codec_internal_error(
	const char *filename,
	int line
) {
	printf("Internal error: File=%s Line=%d\n", filename, line);
	exit(1);
}

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

astc_codec_image* decompress_astc_image(
	const astc_compressed_image& data,
	int bitness,
	astc_decode_mode decode_mode,
	swizzlepattern swz_decode,
	int linearize_srgb,
	int rgb_force_use_of_hdr,
	int alpha_force_use_of_hdr
) {
	astc_codec_image *img = alloc_image(bitness, data.dim_x, data.dim_y, data.dim_z, 0);
	// TODO: Can be null?
	img->linearize_srgb = linearize_srgb;
	img->rgb_force_use_of_hdr = rgb_force_use_of_hdr;
	img->alpha_force_use_of_hdr = alpha_force_use_of_hdr;

	block_size_descriptor bsd;
	init_block_size_descriptor(data.block_x, data.block_y, data.block_z, &bsd);

	int xblocks = (data.dim_x + data.block_x - 1) / data.block_x;
	int yblocks = (data.dim_y + data.block_y - 1) / data.block_y;
	int zblocks = (data.dim_z + data.block_z - 1) / data.block_z;

	imageblock pb;
	for (int z = 0; z < zblocks; z++)
	{
		for (int y = 0; y < yblocks; y++)
		{
			for (int x = 0; x < xblocks; x++)
			{
				int offset = (((z * yblocks + y) * xblocks) + x) * 16;
				uint8_t *bp = data.data + offset;
				physical_compressed_block pcb = *(physical_compressed_block *) bp;
				symbolic_compressed_block scb;
				physical_to_symbolic(&bsd, pcb, &scb);
				decompress_symbolic_block(img, decode_mode, &bsd, x * data.block_x, y * data.block_y, z * data.block_z, &scb, &pb);
				write_imageblock(img, &pb, &bsd, x * data.block_x, y * data.block_y, z * data.block_z, swz_decode);
			}
		}
	}

	term_block_size_descriptor(&bsd);
	return img;
}

struct compress_astc_image_info
{
	const block_size_descriptor* bsd;
	const error_weighting_params* ewp;
	uint8_t* buffer;
	int threadcount;
	astc_decode_mode decode_mode;
	swizzlepattern swz_encode;
	const astc_codec_image* input_image;
};

static void compress_astc_image_threadfunc(
	int thread_count,
	int thread_id,
	void* vblk
) {
	const compress_astc_image_info *blk = (const compress_astc_image_info *)vblk;
	const block_size_descriptor *bsd = blk->bsd;
	int xdim = bsd->xdim;
	int ydim = bsd->ydim;
	int zdim = bsd->zdim;
	uint8_t *buffer = blk->buffer;
	const error_weighting_params *ewp = blk->ewp;
	astc_decode_mode decode_mode = blk->decode_mode;
	swizzlepattern swz_encode = blk->swz_encode;
	const astc_codec_image *input_image = blk->input_image;

	imageblock pb;
	int ctr = thread_id;
	int pctr = 0;

	int x, y, z;
	int xsize = input_image->xsize;
	int ysize = input_image->ysize;
	int zsize = input_image->zsize;
	int xblocks = (xsize + xdim - 1) / xdim;
	int yblocks = (ysize + ydim - 1) / ydim;
	int zblocks = (zsize + zdim - 1) / zdim;

	//allocate memory for temporary buffers
	compress_symbolic_block_buffers temp_buffers;
	temp_buffers.ewb = new error_weight_block;
	temp_buffers.ewbo = new error_weight_block_orig;
	temp_buffers.tempblocks = new symbolic_compressed_block[4];
	temp_buffers.temp = new imageblock;
	temp_buffers.planes2 = new compress_fixed_partition_buffers;
	temp_buffers.planes2->ei1 = new endpoints_and_weights;
	temp_buffers.planes2->ei2 = new endpoints_and_weights;
	temp_buffers.planes2->eix1 = new endpoints_and_weights[MAX_DECIMATION_MODES];
	temp_buffers.planes2->eix2 = new endpoints_and_weights[MAX_DECIMATION_MODES];
	temp_buffers.planes2->decimated_quantized_weights = new float[2 * MAX_DECIMATION_MODES * MAX_WEIGHTS_PER_BLOCK];
	temp_buffers.planes2->decimated_weights = new float[2 * MAX_DECIMATION_MODES * MAX_WEIGHTS_PER_BLOCK];
	temp_buffers.planes2->flt_quantized_decimated_quantized_weights = new float[2 * MAX_WEIGHT_MODES * MAX_WEIGHTS_PER_BLOCK];
	temp_buffers.planes2->u8_quantized_decimated_quantized_weights = new uint8_t[2 * MAX_WEIGHT_MODES * MAX_WEIGHTS_PER_BLOCK];
	temp_buffers.plane1 = temp_buffers.planes2;

	for (z = 0; z < zblocks; z++)
	{
		for (y = 0; y < yblocks; y++)
		{
			for (x = 0; x < xblocks; x++)
			{
				if (ctr == 0)
				{
					int offset = ((z * yblocks + y) * xblocks + x) * 16;
					uint8_t *bp = buffer + offset;
				#ifdef DEBUG_PRINT_DIAGNOSTICS
					if (diagnostics_tile < 0 || diagnostics_tile == pctr)
					{
						print_diagnostics = (diagnostics_tile == pctr) ? 1 : 0;
				#endif
						fetch_imageblock(input_image, &pb, bsd, x * xdim, y * ydim, z * zdim, swz_encode);
						symbolic_compressed_block scb;
						compress_symbolic_block(input_image, decode_mode, bsd, ewp, &pb, &scb, &temp_buffers);
						*(physical_compressed_block*) bp = symbolic_to_physical(bsd, &scb);
				#ifdef DEBUG_PRINT_DIAGNOSTICS
					}
				#endif

					ctr = thread_count - 1;
					pctr++;
				}
				else
					ctr--;
			}
		}
	}

	delete[] temp_buffers.planes2->decimated_quantized_weights;
	delete[] temp_buffers.planes2->decimated_weights;
	delete[] temp_buffers.planes2->flt_quantized_decimated_quantized_weights;
	delete[] temp_buffers.planes2->u8_quantized_decimated_quantized_weights;
	delete[] temp_buffers.planes2->eix1;
	delete[] temp_buffers.planes2->eix2;
	delete   temp_buffers.planes2->ei1;
	delete   temp_buffers.planes2->ei2;
	delete   temp_buffers.planes2;
	delete[] temp_buffers.tempblocks;
	delete   temp_buffers.temp;
	delete   temp_buffers.ewbo;
	delete   temp_buffers.ewb;
}

static void compress_astc_image(
	const astc_codec_image* input_image,
	int xdim,
	int ydim,
	int zdim,
	const error_weighting_params* ewp,
	astc_decode_mode decode_mode,
	swizzlepattern swz_encode,
	uint8_t* buffer,
	int threadcount
) {
	// before entering into the multi-threaded routine, ensure that the block size descriptors
	// and the partition table descriptors needed actually exist.
	block_size_descriptor* bsd = new block_size_descriptor;
	init_block_size_descriptor(xdim, ydim, zdim, bsd);
	get_partition_table(bsd, 0);

	compress_astc_image_info ai;
	ai.bsd = bsd;
	ai.buffer = buffer;
	ai.ewp = ewp;
	ai.decode_mode = decode_mode;
	ai.swz_encode = swz_encode;
	ai.input_image = input_image;

	launch_threads(threadcount, compress_astc_image_threadfunc, &ai);

	term_block_size_descriptor(bsd);
	delete bsd;
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

int astc_main(
	int argc,
	char **argv
) {
	test_inappropriate_extended_precision();

	test_inappropriate_cpu_extensions();

	// initialization routines
	prepare_angular_tables();
	build_quantization_mode_table();

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

	astc_decode_mode decode_mode = DECODE_LDR;
	astc_op_mode op_mode = ASTC_UNRECOGNIZED;

	struct {
		const char *opt;
		astc_op_mode op_mode;
		astc_decode_mode decode_mode;
	} modes[] = {
		{"-cl",      ASTC_ENCODE,            DECODE_LDR},
		{"-dl",      ASTC_DECODE,            DECODE_LDR},
		{"-tl",      ASTC_ENCODE_AND_DECODE, DECODE_LDR},
		{"-cs",      ASTC_ENCODE,            DECODE_LDR_SRGB},
		{"-ds",      ASTC_DECODE,            DECODE_LDR_SRGB},
		{"-ts",      ASTC_ENCODE_AND_DECODE, DECODE_LDR_SRGB},
		{"-ch",      ASTC_ENCODE,            DECODE_HDR},
		{"-dh",      ASTC_DECODE,            DECODE_HDR},
		{"-th",      ASTC_ENCODE_AND_DECODE, DECODE_HDR},
		{"-cH",      ASTC_ENCODE,            DECODE_HDRA},
		{"-dH",      ASTC_DECODE,            DECODE_HDRA},
		{"-tH",      ASTC_ENCODE_AND_DECODE, DECODE_HDRA},
		{"-h",       ASTC_PRINT_LONGHELP,    DECODE_HDR},
		{"-help",    ASTC_PRINT_LONGHELP,    DECODE_HDR},
		{"-v",       ASTC_PRINT_VERSION,     DECODE_HDR},
		{"-version", ASTC_PRINT_VERSION,     DECODE_HDR}
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
			printf("Unrecognized operation mode \"%s\"\n", argv[1]);
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
	int rgb_force_use_of_hdr = 0;
	int alpha_force_use_of_hdr = 0;

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

	swizzlepattern swz_encode = { 0, 1, 2, 3 };
	swizzlepattern swz_decode = { 0, 1, 2, 3 };

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

		if (decode_mode == DECODE_HDRA)
		{
			ewp.rgb_power = 0.75f;
			ewp.rgb_base_weight = 0.0f;
			ewp.rgb_mean_weight = 1.0f;

			ewp.alpha_power = 0.75f;
			ewp.alpha_base_weight = 0.0f;
			ewp.alpha_mean_weight = 1.0f;

			rgb_force_use_of_hdr = 1;
			alpha_force_use_of_hdr = 1;

			dblimit_set = 999.0f;
		}
		else if (decode_mode == DECODE_HDR)
		{
			ewp.rgb_power = 0.75f;
			ewp.rgb_base_weight = 0.0f;
			ewp.rgb_mean_weight = 1.0f;

			ewp.alpha_base_weight = 0.05f;

			rgb_force_use_of_hdr = 1;
			alpha_force_use_of_hdr = 0;

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

			int swizzle_components[4];
			for (int i = 0; i < 4; i++)
				switch (argv[argidx - 1][i])
				{
				case 'r':
					swizzle_components[i] = 0;
					break;
				case 'g':
					swizzle_components[i] = 1;
					break;
				case 'b':
					swizzle_components[i] = 2;
					break;
				case 'a':
					swizzle_components[i] = 3;
					break;
				case '0':
					swizzle_components[i] = 4;
					break;
				case '1':
					swizzle_components[i] = 5;
					break;
				default:
					printf("Character '%c' is not a valid swizzle-character\n", argv[argidx - 1][i]);
					return 1;
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

			int swizzle_components[4];
			for (int i = 0; i < 4; i++)
			{
				switch (argv[argidx - 1][i])
				{
				case 'r':
					swizzle_components[i] = 0;
					break;
				case 'g':
					swizzle_components[i] = 1;
					break;
				case 'b':
					swizzle_components[i] = 2;
					break;
				case 'a':
					swizzle_components[i] = 3;
					break;
				case '0':
					swizzle_components[i] = 4;
					break;
				case '1':
					swizzle_components[i] = 5;
					break;
				case 'z':
					swizzle_components[i] = 6;
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
			swz_encode.r = 0;	// r <- red
			swz_encode.g = 0;	// g <- red
			swz_encode.b = 0;	// b <- red
			swz_encode.a = 1;	// a <- green
			swz_decode.r = 0;	// r <- red
			swz_decode.g = 3;	// g <- alpha
			swz_decode.b = 6;	// b <- reconstruct
			swz_decode.a = 5;	// 1.0

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
			swz_encode.r = 0;	// r <- red
			swz_encode.g = 0;	// g <- red
			swz_encode.b = 0;	// b <- red
			swz_encode.a = 1;	// a <- green
			swz_decode.r = 0;	// r <- red
			swz_decode.g = 3;	// g <- alpha
			swz_decode.b = 6;	// b <- reconstruct
			swz_decode.a = 5;	// 1.0

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
	#ifdef DEBUG_PRINT_DIAGNOSTICS
		else if (!strcmp(argv[argidx], "-diag"))
		{
			argidx += 2;
			if (argidx > argc)
			{
				printf("-diag switch with no argument\n");
				return 1;
			}

			diagnostics_tile = atoi(argv[argidx - 1]);
		}
		else if (!strcmp(argv[argidx], "-pte"))
		{
			argidx++;
			print_tile_errors = 1;
		}
		else if (!strcmp(argv[argidx], "-stats"))
		{
			argidx++;
			print_statistics = 1;
		}
	#endif
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

	float texel_avg_error_limit = 0.0f;
	if (op_mode == ASTC_ENCODE || op_mode == ASTC_ENCODE_AND_DECODE)
	{
		int partitions_to_test = plimit_set;
		float oplimit = oplimit_set;
		float mincorrel = mincorrel_set;

		int maxiters = maxiters_set;
		ewp.max_refinement_iters = maxiters;
		ewp.block_mode_cutoff = bmc_set / 100.0f;

		if (rgb_force_use_of_hdr == 0)
		{
			texel_avg_error_limit = powf(0.1f, dblimit_set * 0.1f) * 65535.0f * 65535.0f;
		}
		else
		{
			texel_avg_error_limit = 0.0f;
		}

		ewp.partition_1_to_2_limit = oplimit;
		ewp.lowest_correlation_cutoff = mincorrel;
		ewp.partition_search_limit = astc::clampi(partitions_to_test, 1, PARTITION_COUNT);

		// if diagnostics are run, force the thread count to 1.
		#ifdef DEBUG_PRINT_DIAGNOSTICS
			if (diagnostics_tile >= 0 ||
				print_tile_errors > 0 ||
				print_statistics > 0)
			{
				thread_count = 1;
			}
		#endif

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
				case DECODE_LDR:
					printf("    Color profile: LDR linear\n");
					break;
				case DECODE_LDR_SRGB:
					printf("    Color profile: LDR sRGB\n");
					break;
				case DECODE_HDR:
					printf("    Color profile: HDR RGB + LDR A\n");
					break;
				case DECODE_HDRA:
					printf("    Color profile: HDR RGBA\n");
					break;
			}

			if (block_z == 1)
			{
				printf("    Block size: %dx%d (%.2f bpp)\n", block_x, block_y, 128.0 / (block_x * block_y));
			}
			else
			{
				printf("    Block size: %dx%dx%d (%.2f bpp)\n", block_x, block_y, block_z, 128.0 / (block_x * block_y * block_z));
			}

			printf("    Radius mean/stdev: %d texels\n", ewp.mean_stdev_radius);
			printf("    RGB power: %g\n", (double)ewp.rgb_power);
			printf("    RGB base weight: %g\n", (double)ewp.rgb_base_weight);
			printf("    RGB mean weight: %g\n", (double)ewp.rgb_mean_weight);
			printf("    RGB stdev weight: %g\n", (double)ewp.rgb_stdev_weight);
			printf("    RGB mean/stdev channel mixing: %g\n", (double)ewp.rgb_mean_and_stdev_mixing);
			printf("    Alpha power: %g\n", (double)ewp.alpha_power);
			printf("    Alpha base weight: %g\n", (double)ewp.alpha_base_weight);
			printf("    Alpha mean weight: %g\n", (double)ewp.alpha_mean_weight);
			printf("    Alpha stdev weight: %g\n", (double)ewp.alpha_stdev_weight);
			printf("    RGB alpha scale weight: %d\n", ewp.enable_rgb_scale_with_alpha);
			if (ewp.enable_rgb_scale_with_alpha)
			{
				printf("    RGB alpha scale radius: %d)\n", ewp.alpha_radius);
			}

			printf("    Color channel weight: R=%g G=%g B=%g A=%g\n", (double)ewp.rgba_weights[0], (double)ewp.rgba_weights[1], (double)ewp.rgba_weights[2], (double)ewp.rgba_weights[3]);
			printf("    Deblock artifact weight: %g\n", (double)ewp.block_artifact_suppression);
			printf("    Block partition cutoff: %d\n", ewp.partition_search_limit);
			printf("    PSNR decibel cutoff: %g\n", (double)dblimit_set);
			printf("    1->2 partition cutoff: %g\n", (double)oplimit);
			printf("    2 plane correlation cutoff: %g\n", (double)mincorrel_set);
			printf("    Block mode centile cutoff: %g\n", (double)(ewp.block_mode_cutoff * 100.0f));
			printf("    Max refinement cutoff: %d\n", ewp.max_refinement_iters);
			printf("    Compressor thread count: %d\n", thread_count);
			printf("\n");
		}
	}

	int padding = MAX(ewp.mean_stdev_radius, ewp.alpha_radius);

	// determine encoding bitness as follows:
	// if enforced by the output format, follow the output format's result
	// else use decode_mode to pick bitness.
	int out_bitness = (op_mode == ASTC_DECODE || op_mode == ASTC_ENCODE_AND_DECODE) ? get_output_filename_enforced_bitness(output_filename) : -1;
	if (out_bitness == -1)
	{
		out_bitness = (decode_mode == DECODE_HDR) ? 16 : 8;
	}

	// Temporary image array (for merging multiple 2D images into one 3D image).
	int load_result = 0;
	astc_codec_image *input_decomp_img = nullptr;
	//astc_compressed_image *input_comp_img = nullptr;
	astc_codec_image *output_decomp_img = nullptr;
	//astc_compressed_image *output_comp_img = nullptr;

	int input_components = 0;
	int input_image_is_hdr = 0;

	// Load uncompressed inputs
	if (op_mode == ASTC_ENCODE || op_mode == ASTC_ENCODE_AND_DECODE)
	{
		// Allocate temporary arrays for image data and load results
		int* load_results = new int [array_size];
		astc_codec_image** input_images = new astc_codec_image* [array_size];

		// Iterate over all input images
		for (int image_index = 0; image_index < array_size; image_index++)
		{
			// 2D input data
			if (array_size == 1)
			{
				input_images[image_index] = astc_codec_load_image(
					input_filename, padding, y_flip, linearize_srgb,
					rgb_force_use_of_hdr, alpha_force_use_of_hdr,
					&load_results[image_index]);
			}
			// 3D input data - multiple 2D images
			else
			{
				// TODO: These can overflow if the file name is longer than 256
				char new_input_filename[256];

				// Check for extension: <name>.<extension>
				if (!strrchr(input_filename, '.'))
				{
					printf("ERROR: Unable to determine file extension: %s\n", input_filename);
					return 1;
				}

				// Construct new file name and load: <name>_N.<extension>
				strcpy(new_input_filename, input_filename);
				sprintf(strrchr(new_input_filename, '.'), "_%d%s", image_index, strrchr(input_filename, '.'));
				input_images[image_index] = astc_codec_load_image
					(new_input_filename, padding, y_flip, linearize_srgb,
					rgb_force_use_of_hdr, alpha_force_use_of_hdr,
					&load_results[image_index]);

				// If image loaded correctly, check image is not 3D.
				if ((load_results[image_index] >= 0) &&
				    (input_images[image_index]->zsize != 1))
				{
					printf("3D source images not supported with -array option: %s\n", new_input_filename);
					return 1;
				}
			}

			// Check load result.
			if (load_results[image_index] < 0)
			{
				return 1;
			}

			// Check format matches other slices.
			if (load_results[image_index] != load_results[0])
			{
				printf("Mismatching image format - image 0 and %d are a different format\n", image_index);
				return 1;
			}
		}

		load_result = load_results[0];

		// Assign input image.
		if (array_size == 1)
		{
			input_decomp_img = input_images[0];
		}
		// Merge input image data.
		else
		{
			int z, xsize, ysize, zsize, bitness, slice_size;

			xsize = input_images[0]->xsize;
			ysize = input_images[0]->ysize;
			zsize = array_size;
			bitness = (load_result & 0x80) ? 16 : 8;
			slice_size = (xsize + (2 * padding)) * (ysize + (2 * padding));

			// Allocate image memory.
			input_decomp_img = alloc_image(bitness, xsize, ysize, zsize, padding);

			// Combine 2D source images into one 3D image (skip padding slices as these don't exist in 2D textures).
			for (z = padding; z < zsize + padding; z++)
			{
				if (bitness == 8)
				{
					memcpy(*input_decomp_img->data8[z], *input_images[z - padding]->data8[0], slice_size * 4 * sizeof(uint8_t));
				}
				else
				{
					memcpy(*input_decomp_img->data16[z], *input_images[z - padding]->data16[0], slice_size * 4 * sizeof(uint16_t));
				}
			}

			// Clean up temporary images.
			for (int i = 0; i < array_size; i++)
			{
				free_image(input_images[i]);
			}

			// Clamp texels outside the actual image area.
			fill_image_padding_area(input_decomp_img);
		}

		delete[] input_images;
		input_images = nullptr;

		delete[] load_results;
		load_results = nullptr;

		input_components = load_result & 7;
		input_image_is_hdr = (load_result & 0x80) ? 1 : 0;

		if ((input_decomp_img->zsize > 1) && (block_z == 1))
		{
			printf("ERROR: 3D input data for a 2D ASTC block format\n");
			exit(1);
		}

		ewp.texel_avg_error_limit = texel_avg_error_limit;
		expand_block_artifact_suppression(block_x, block_y, block_z, &ewp);

		if (!silentmode)
		{
			printf("Source image\n");
			printf("============\n\n");
			printf("    Source: %s\n", input_filename);
			printf("    Color profile: %s\n", input_image_is_hdr ? "HDR" : "LDR");
			if (input_decomp_img->zsize > 1)
			{
				printf("    Dimensions: 3D, %d x %d x %d\n",
				       input_decomp_img->xsize, input_decomp_img->ysize, input_decomp_img->zsize);
			}
			else
			{
				printf("    Dimensions: 2D, %d x %d\n",
				       input_decomp_img->xsize, input_decomp_img->ysize);
			}
			printf("    Channels: %d\n\n", load_result & 7);
		}

		if (padding > 0 ||
		    ewp.rgb_mean_weight != 0.0f || ewp.rgb_stdev_weight != 0.0f ||
		    ewp.alpha_mean_weight != 0.0f || ewp.alpha_stdev_weight != 0.0f)
		{
			if (!silentmode)
			{
				printf("Computing texel-neighborhood means and variances ... ");
				fflush(stdout);
			}

			compute_averages_and_variances(
				input_decomp_img,
				ewp.rgb_power,
				ewp.alpha_power,
				ewp.mean_stdev_radius,
				ewp.alpha_radius,
				linearize_srgb,
				swz_encode,
				thread_count);

			if (!silentmode)
			{
				printf("done\n");
			}
		}
	}

	start_coding_time = get_time();

	if (op_mode == ASTC_DECODE)
	{
		astc_compressed_image comp_img;
		int error = load_astc_file(input_filename, comp_img);
		if (error)
		{
			return 1;
		}

		output_decomp_img = decompress_astc_image(
		    comp_img, out_bitness, decode_mode, swz_decode,
		    linearize_srgb, rgb_force_use_of_hdr, alpha_force_use_of_hdr);

		delete[] comp_img.data;
	}

	// process image, if relevant
	if (op_mode == ASTC_ENCODE_AND_DECODE)
	{
		int xsize = input_decomp_img->xsize;
		int ysize = input_decomp_img->ysize;
		int zsize = input_decomp_img->zsize;
		int xblocks = (xsize + block_x - 1) / block_x;
		int yblocks = (ysize + block_y - 1) / block_y;
		int zblocks = (zsize + block_z - 1) / block_z;
		uint8_t* buffer = new uint8_t [xblocks * yblocks * zblocks * 16];
		compress_astc_image(input_decomp_img, block_x, block_y, block_z, &ewp, decode_mode, swz_encode, buffer, thread_count);

		astc_compressed_image comp_img {
			block_x, block_y, block_z,
			xsize, ysize, zsize,
			buffer
		};

		output_decomp_img = decompress_astc_image(
		    comp_img, out_bitness, decode_mode, swz_decode,
		    0, rgb_force_use_of_hdr, alpha_force_use_of_hdr);

		delete[] buffer;
	}

	end_coding_time = get_time();

	// print PSNR if encoding
	if (op_mode == ASTC_ENCODE_AND_DECODE)
	{
		compute_error_metrics(input_image_is_hdr, input_components, input_decomp_img, output_decomp_img, low_fstop, high_fstop);
	}

	// store image
	if (op_mode == ASTC_DECODE || op_mode == ASTC_ENCODE_AND_DECODE)
	{
		int store_result = -1;
		const char *format_string = "";

		store_result = astc_codec_store_image(output_decomp_img, output_filename, &format_string, y_flip);
		if (store_result < 0)
		{
			printf("ERROR: Failed to write output image %s\n", output_filename);
			return 1;
		}
	}

	if (op_mode == ASTC_ENCODE)
	{
		int xsize = input_decomp_img->xsize;
		int ysize = input_decomp_img->ysize;
		int zsize = input_decomp_img->zsize;
		int xblocks = (xsize + block_x - 1) / block_x;
		int yblocks = (ysize + block_y - 1) / block_y;
		int zblocks = (zsize + block_z - 1) / block_z;
		uint8_t* buffer = new uint8_t [xblocks * yblocks * zblocks * 16];
		compress_astc_image(input_decomp_img, block_x, block_y, block_z, &ewp, decode_mode, swz_encode, buffer, thread_count);

		end_coding_time = get_time();

		astc_compressed_image comp_img {
			block_x, block_y, block_z,
			xsize, ysize, zsize,
			buffer
		};

		int error = store_astc_file(comp_img, output_filename);
		if (error)
		{
			return 1;
		}
	}

	free_image(input_decomp_img);
	end_time = get_time();

	if (op_mode == ASTC_ENCODE_AND_DECODE || !silentmode)
	{
		printf("\nCoding time\n");
		printf("-----------\n");
		printf("Total time:             %6.2f s\n", end_time - start_time);
		printf("Coding time:            %6.2f s\n", end_coding_time - start_coding_time);
	}

	return 0;
}
