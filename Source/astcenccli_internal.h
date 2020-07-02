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
 * @brief Functions and data declarations.
 */

#ifndef ASTCENCCLI_INTERNAL_INCLUDED
#define ASTCENCCLI_INTERNAL_INCLUDED

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include "astcenc.h"
#include "astcenc_mathlib.h"

struct astc_compressed_image
{
	unsigned int block_x;
	unsigned int block_y;
	unsigned int block_z;
	unsigned int dim_x;
	unsigned int dim_y;
	unsigned int dim_z;
	uint8_t* data;
};

// Config options to be read from command line
struct cli_config_options
{
	unsigned int thread_count;
	unsigned int array_size;
	bool silentmode;
	bool y_flip;
	bool linearize_srgb;
	int low_fstop;
	int high_fstop;
	astcenc_swizzle swz_encode;
	astcenc_swizzle swz_decode;
};

/**
 * Functions to load image from file.
 *
 * @param filename            The file path on disk.
 * @param padding             The texel padding needed around the image.
 * @param y_flip              Should this image be Y flipped?
 * @param[out] is_hdr         Is the loaded image HDR?
 * @param[out] num_components The number of components in the loaded image.
 *
 * @return The astc image file, or nullptr on error.
 */
astcenc_image* load_ncimage(
	const char* filename,
	unsigned int dim_pad,
	bool y_flip,
	bool& is_hdr,
	unsigned int& num_components);

int store_ncimage(
	const astcenc_image* output_image,
	const char* output_filename,
	const char** file_format_name,
	int y_flip);

int get_output_filename_enforced_bitness(
	const char* filename);


astcenc_image* alloc_image(
	unsigned int bitness,
	unsigned int dim_x,
	unsigned int dim_y,
	unsigned int dim_z,
	unsigned int dim_pad);

void free_image(
	astcenc_image* img);

void fill_image_padding_area(
	astcenc_image* img);

int determine_image_channels(
	const astcenc_image* img);

int load_cimage(
	const char* filename,
	astc_compressed_image& out_image);

int store_cimage(
	const astc_compressed_image& comp_img,
	const char* filename);

// helper functions to prepare an ASTC image object from a flat array
// Used by the image loaders in "astc_file_load_store.cpp"
astcenc_image* astc_img_from_floatx4_array(
	const float* data,
	unsigned int dim_x,
	unsigned int dim_y,
	unsigned int dim_pad,
	bool y_flip);

astcenc_image*astc_img_from_unorm8x4_array(
	const uint8_t* data,
	unsigned int dim_x,
	unsigned int dim_y,
	unsigned int dim_pad,
	bool y_flip);

// helper functions to prepare a flat array from an ASTC image object.
// the array is allocated with new[], and must be freed with delete[].
float* floatx4_array_from_astc_img(
	const astcenc_image* img,
	bool y_flip);

uint8_t* unorm8x4_array_from_astc_img(
	const astcenc_image* img,
	bool y_flip);

/* ============================================================================
  Functions for printing build info and help messages
============================================================================ */
void astcenc_print_header();

void astcenc_print_shorthelp();

void astcenc_print_longhelp();

/**
 * @brief Compute error metrics comparing two images.
 *
 * @param compute_hdr_metrics Non-zero if HDR metrics should be computed.
 * @param input_components    The number of input color components.
 * @param img1                The original iamge.
 * @param img2                The compressed image.
 * @param fstop_lo            The low exposure fstop (HDR only).
 * @param fstop_hi            The high exposure fstop (HDR only).
 * @param show_psnr           Non-zero if metrics should be logged to stdout.
 */
void compute_error_metrics(
	int compute_hdr_metrics,
	int input_components,
	const astcenc_image* img1,
	const astcenc_image* img2,
	int fstop_lo,
	int fstop_hi);

/**
 * @brief Get the current time.
 *
 * @returns The current time in seconds since arbitrary epoch.
 */
double get_time();

#endif
