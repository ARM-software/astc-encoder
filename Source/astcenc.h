// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2020 Arm Limited
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
 * @brief The core astcenc codec library interface.
 *
 * This interface is the entry point to the core astcenc codec. It aims to be
 * easy to use for non-experts, but also to allow experts to have fine control
 * over the compressor heuristics. The core codec only handles compression and
 * decompression, transferring all inputs and outputs via memory buffers.
 *
 * While the aim is that we keep this interface mostly stable, it should be
 * viewed as a mutable interface tied to a specific source version. We are not
 * trying to maintain backwards compatibility across codec versions.
 *
 * The API state management is based around an explicit context object. The
 * context is the container for all allocated resources -- memory and threads
 * -- needed to compress and decompress for a single image configuration. A
 * context keeps the state needed to process one image at a time, but can be
 * used to sequentially compress multiple images using the same configuration.
 * This allows setup overheads to be amortized over multiple images, which is
 * particularly important when images are small. An application can process
 * multiple images concurrently by allocating multiple contexts.
 *
 * The API allows compression of a single image by using multiple threads
 * within the scope of a single context. The API supports two mechanisms for
 * this: automatic and manual. Automatic threading uses threads created and
 * managed by the codec itself; no user setup needed - just configure the
 * thread count and call compress from a single application thread.
 * Manual threading uses user-created threads, which is a common pattern in
 * some developer tools using persistent pools of worker threads. For this
 * configure the codec with the maximum thread count you'll use and then invoke
 * compress from each worker thread.
 *
 * To catch obvious input/output buffer sizing issues, which can cause security
 * and stability problems, all transfer buffers are explicitly sized.
 *
 * Automatic threading
 * ===================
 *
 * In pseudocode, the usage for automatic threading looks like this:
 *
 *     // Configure the compressor run
 *     astcenc_config my_config;
 *     astcenc_init_config(..., &my_config);
 *
 *     // Allocate working state given config and thread_count
 *     astcenc_context* my_context;
 *     astcenc_context_alloc(&my_config, thread_count, &my_context);
 *
 *     // Compress the image with automatic multithreading
 *     foreach image:
 *         astc_image input;
 *         uint8_t* output;
 *         astcenc_compress_image(&my_context, &input, output, -1);
 *
 *     // Clean up
 *     astcenc_context_free(my_context);
 *
 * Manual threading
 * ================
 *
 * In pseudocode, the usage for manual user threading looks like this:
 *
 *     // Configure the compressor run
 *     astcenc_config my_config;
 *     astcenc_flags my_flags = ASTCENC_USE_USER_THREADS | ...;
 *     astcenc_init_config(..., my_flags, &my_config);
 *
 *     // Allocate working state given config and thread_count
 *     astcenc_context* my_context;
 *     astcenc_context_alloc(&my_config, thread_count, &my_context);
 *
 *     // Compress the image with manually managed multithreading
 *     foreach image:
 *         astc_image input;
 *         uint8_t* output;
 *         for i in range(0, thread_count):
 *             astcenc_compress_image(&my_context, &input, output, i);
 *
 *         // Wait for all threads to complete one image
 *         barrier();
 *
 *     // Clean up
 *     astcenc_context_free(my_context);
 */

#ifndef ASTCENC_INCLUDED
#define ASTCENC_INCLUDED

#include <cstddef>
#include <cstdint>

/* ============================================================================
    Data declarations
============================================================================ */
struct astcenc_context;

// Return codes
enum astcenc_error {
	ASTCENC_SUCCESS = 0,
	ASTCENC_ERR_OUT_OF_MEM,
	ASTCENC_ERR_BAD_PARAM,
	ASTCENC_ERR_BAD_BLOCK_SIZE,
	ASTCENC_ERR_BAD_PROFILE,
	ASTCENC_ERR_BAD_PRESET,
	ASTCENC_ERR_BAD_FLAGS
};

// Compression color feature profile
enum astcenc_profile {
	ASTCENC_PRF_LDR_RGBA = 0,
	ASTCENC_PRF_LDR_SRGB_RGBA,
	ASTCENC_PRF_HDR_RGB_LDR_A,
	ASTCENC_PRF_HDR_RGBA
};

// Compression quality preset
enum astcenc_preset {
	ASTCENC_PRE_FAST = 0,
	ASTCENC_PRE_MEDIUM,
	ASTCENC_PRE_THOROUGH,
	ASTCENC_PRE_EXHAUSTIVE
};

// Incoming data type
enum astcenc_type {
	ASTCENC_TYP_U8 = 0,
	ASTCENC_TYP_F16 = 1,
	ASTCENC_TYP_F32 = 2
};

// Incoming data type
enum astcenc_swz {
	ASTCENC_SWZ_R = 0,
	ASTCENC_SWZ_G = 1,
	ASTCENC_SWZ_B = 2,
	ASTCENC_SWZ_A = 3,
	ASTCENC_SWZ_0 = 4,
	ASTCENC_SWZ_1 = 5,
	ASTCENC_SWZ_Z = 6
};

// Incoming data type
struct astcenc_swizzle {
	astcenc_swz r;
	astcenc_swz g;
	astcenc_swz b;
	astcenc_swz a;
};

// Config mode flags
static const unsigned int ASTCENC_FLG_MAP_NORMAL       = 1 << 0;
static const unsigned int ASTCENC_FLG_MAP_MASK         = 1 << 1;
static const unsigned int ASTCENC_FLG_WEIGHT_ALPHA     = 1 << 2;
static const unsigned int ASTCENC_FLG_USE_PERCEPTUAL   = 1 << 3;
static const unsigned int ASTCENC_FLG_USE_USER_THREADS = 1 << 4;
static const unsigned int ASTCENC_FLG_ENABLE_LOGGING   = 1 << 5;

static const unsigned int ASTCENC_ALL_FLAGS =
                              ASTCENC_FLG_MAP_NORMAL |
                              ASTCENC_FLG_MAP_MASK |
                              ASTCENC_FLG_WEIGHT_ALPHA |
                              ASTCENC_FLG_USE_PERCEPTUAL |
                              ASTCENC_FLG_USE_USER_THREADS |
                              ASTCENC_FLG_ENABLE_LOGGING;

// Config structure
struct astcenc_config {
	astcenc_profile profile;
	unsigned int flags;
	unsigned int block_x;
	unsigned int block_y;
	unsigned int block_z;
	unsigned int v_rgba_radius;
	float v_rgba_mean_stdev_mix;
	float v_rgb_power;
	float v_rgb_base;
	float v_rgb_mean;
	float v_rgb_stdev;
	float v_a_power;
	float v_a_base;
	float v_a_mean;
	float v_a_stdev;
	float cw_r_weight;
	float cw_g_weight;
	float cw_b_weight;
	float cw_a_weight;
	unsigned int a_scale_radius;
	float b_deblock_weight;
	unsigned int tune_partition_limit;
	unsigned int tune_block_mode_limit;
	unsigned int tune_refinement_limit;
	float tune_db_limit;
	float tune_partition_early_out_limit;
	float tune_two_plane_early_out_limit;
};

/**
 * TODO: What memory layout do we pass this in as? Current code uses 3D arrays
 * with unique allocations per row, which is relatively heavy, but flexible.
 */
struct astcenc_image {
	unsigned int dim_x;
	unsigned int dim_y;
	astcenc_type type;
	astcenc_swizzle swizzle;
	void* data;        // Data stored as RGBA texels of "type" primitive.
	size_t row_stride; // Row stride in bytes; 0 for tightly packed.
	size_t data_size;  // Data length in bytes.
};

/**
 * Populate a compressor config based on default settings.
 *
 * Note: Advanced users can edit the returned config struct to apply manual
 * fine tuning before creating the context. Non-power users will not need to.
 *
 * Will error if any of the inputs are invalid, either individually or in
 * combination.
 */
astcenc_error astcenc_init_config(
	astcenc_profile profile,
	unsigned int block_x,
	unsigned int block_y,
	unsigned int block_z,
	astcenc_preset preset,
	unsigned int flags,
	astcenc_config* config);

/**
 * Allocate a new compressor context based on the settings in the config.
 *
 * The persistent parts of the config are copied by this function; it can be
 * freed after this function is called.
 *
 * A single context can be used to sequentially compress or decompress multiple
 * images using the same configuration settings.
 *
 * This will configure the context for N threads, allocated all of the data
 * tables and per-thread working buffers.
 */
astcenc_error astcenc_context_alloc(
	astcenc_config const* config,
	int thread_count,
	astcenc_context** context);

/**
 * Compress the image.
 *
 * If the config did not specify USE_USER_THREADS this function will
 * automatically create and start the N threads needed to run the compressor.
 */
astcenc_error astcenc_compress_image(
	astcenc_context* context,
	astcenc_image const* image,
	size_t image_len, // Length of image array (volumetric image Z dimension)
	uint8_t* data_out,
	size_t data_len, // Length of output array for bounds checking.
	int thread_index);

/**
 * Decompress the image.
 *
 * If the config did not specify USE_USER_THREADS this function will
 * automatically use the context thread count to decompress the image.
 * Otherwise the user must manually invoke this function with thread_index
 * set between [0, thread_count - 1], once for each index.
 */
astcenc_error astcenc_decompress_image(
	astcenc_context * context,
	void* data,
	size_t data_len,
	astcenc_image * image_out,
	size_t image_out_len, // Length of image array (volumetric image Z dimension)
	int thread_index);

/**
 * Free the compressor context.
 *
 * This must only be called when no other threads are inside a codec function
 * using this context.
 */
void astcenc_context_free(
	astcenc_context* context);

#endif
