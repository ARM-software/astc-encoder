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
 * over the compressor heuristics if needed. The core codec only handles
 * compression and decompression, transferring all inputs and outputs via
 * memory buffers. To catch obvious input/output buffer sizing issues, which
 * can cause security and stability problems, all transfer buffers are
 * explicitly sized.
 *
 * While the aim is that we keep this interface mostly stable, it should be
 * viewed as a mutable interface tied to a specific source version. We are not
 * trying to maintain backwards compatibility across codec versions.
 *
 * The API state management is based around an explicit context object, which
 * is the context for all allocated memory resources needed to compress and
 * decompress a single image. A context can be used to sequentially compress
 * multiple images using the same configuration, allowing setup overheads to be
 * amortized over multiple images, which is particularly important when images
 * are small.
 *
 * Multi-threading can be used two ways.
 *
 *     * An application wishing to process multiple images in parallel can
 *       allocate multiple contexts and assign each context to a thread.
 *     * An application wishing to process a single image in using multiple
 *       threads can configure the context for multi-threaded use, and invoke
 *       astcenc_compress() once per thread for faster compression. The caller
 *       is responsible for creating the worker threads. Note that
 *       decompression is always single-threaded.
 *
 * Threading
 * =========
 *
 * In pseudocode, the usage for manual user threading looks like this:
 *
 *     // Configure the compressor run
 *     astcenc_config my_config;
 *     astcenc_config_init(..., &my_config);
 *
 *     // Power users can tune the tweak <my_config> settings here ...
 *
 *     // Allocate working state given config and thread_count
 *     astcenc_context* my_context;
 *     astcenc_context_alloc(&my_config, thread_count, &my_context);
 *
 *     // Compress each image using these config settings
 *     foreach image:
 *         // For each thread in the thread pool
 *         for i in range(0, thread_count):
 *             astcenc_compress_image(my_context, &my_input, my_output, i);
 *
 *         astcenc_compress_reset(my_context);
 *
 *     // Clean up
 *     astcenc_context_free(my_context);
 *
 * Images
 * ======
 *
 * Images are passed in as a astcenc_image structure. Inputs can be either
 * 8-bit unorm inputs (passed in via the data8 pointer), or 16-bit floating
 * point inputs (passed in via the data16 pointer). The unused pointer should
 * be set to nullptr.
 *
 * Images can be any dimension; there is no requirement for them to be a
 * multiple of the ASTC block size.
 *
 * Data is always passed in as 4 color channels, and accessed as 3D array
 * indexed using e.g.
 *
 *     data8[z_coord][y_coord][x_coord * 4    ]   // Red
 *     data8[z_coord][y_coord][x_coord * 4 + 1]   // Green
 *     data8[z_coord][y_coord][x_coord * 4 + 2]   // Blue
 *     data8[z_coord][y_coord][x_coord * 4 + 3]   // Alpha
 *
 * If a region-based error heuristic is used for compression, the input image
 * must be padded on all sides by pad_dim texels in x and y dimensions (and z
 * dimensions for 3D images). The padding region must be filled by
 * extrapolating the nearest edge color. The required padding size is given by
 * the following config settings:
 *
 *     max(config.v_rgba_radius, config.a_scale_radius)
 *
 * This can be programatically determined by reading the config containing the
 * values passed into astcenc_context_alloc().
 *
 * Common compressor usage
 * =======================
 *
 * One of the most important things for coding image quality is to align the
 * input data channel count with the ASTC color endpoint mode. This avoids
 * wasting bits encoding channels you don't need in the endpoint colors.
 *
 *         | Input data | Encoding swizzle | Sampling swizzle |
 *         | ---------- | ---------------- | ---------------- |
 *         | 1 channel  | RRR1             | .g               |
 *         | 2 channels | RRRG             | .ga              |
 *         | 3 channels | RGB1             | .rgb             |
 *         | 4 channels | RGBA             | .rgba            |
 *
 * The 1 and 2 channel modes recommend sampling from "g" to recover the
 * luminance value as this provide best compatibility with ETC1. For ETC the
 * luminance data will be stored as RGB565 where the green channel has the
 * best quality. For ASTC any of the rgb channels can be used - the same data
 * will be returned for all three.
 *
 * When using the normal map compression mode ASTC will store normals as a two
 * channel X+Y map. Input images must contain unit-length normalized and should
 * be passed in using the RRRG swizzle. The Z component can be programatically
 * recovered in shader code, using knowledge that the vector is unit length.
 *
 * Decompress-only usage
 * =====================
 *
 * For some use cases it is useful to have a cut-down context and/or library
 * which supports decompression but not compression.
 *
 * A context can be made decompress-only using the ASTCENC_FLG_DECOMPRESS_ONLY
 * flag when the context is allocated. These contexts have lower dynamic memory
 * footprint than a full context.
 *
 * The entire library can be made decompress-only by building the files with
 * the define ASTCENC_DECOMPRESS_ONLY set. In this build the context will be
 * smaller, and the library will exclude the functionality which is only needed
 * for compression. This reduces the binary size by ~180KB. For these builds
 * contexts must be created with the ASTCENC_FLG_DECOMPRESS_ONLY flag.
 *
 * Note that context structures returned by a library built as decompress-only
 * are incompatible with a library built with compression included, and visa
 * versa, as they have different sizes and memory layout.
 */

#ifndef ASTCENC_INCLUDED
#define ASTCENC_INCLUDED

#include <cstddef>
#include <cstdint>

/* ============================================================================
    Data declarations
============================================================================ */

/**
 * @brief An opaque structure; see astcenc_internal.h for definition.
 */
struct astcenc_context;

/**
 * @brief A codec API error code.
 */
enum astcenc_error {
	/** @brief The call was successful. */
	ASTCENC_SUCCESS = 0,
	/** @brief The call failed due to low memory, or undersized I/O buffers. */
	ASTCENC_ERR_OUT_OF_MEM,
	/** @brief The call failed due to the build using fast math. */
	ASTCENC_ERR_BAD_CPU_FLOAT,
	/** @brief The call failed due to the build using an unsupported ISA. */
	ASTCENC_ERR_BAD_CPU_ISA,
	/** @brief The call failed due to an out-of-spec parameter. */
	ASTCENC_ERR_BAD_PARAM,
	/** @brief The call failed due to an out-of-spec block size. */
	ASTCENC_ERR_BAD_BLOCK_SIZE,
	/** @brief The call failed due to an out-of-spec color profile. */
	ASTCENC_ERR_BAD_PROFILE,
	/** @brief The call failed due to an out-of-spec quality preset. */
	ASTCENC_ERR_BAD_PRESET,
	/** @brief The call failed due to an out-of-spec channel swizzle. */
	ASTCENC_ERR_BAD_SWIZZLE,
	/** @brief The call failed due to an out-of-spec flag set. */
	ASTCENC_ERR_BAD_FLAGS,
	/** @brief The call failed due to the context not supporting the operation. */
	ASTCENC_ERR_BAD_CONTEXT,
	/** @brief The call failed due to unimplemented functionality. */
	ASTCENC_ERR_NOT_IMPLEMENTED
};

/**
 * @brief A codec color profile.
 */
enum astcenc_profile {
	/** @brief The LDR sRGB color profile. */
	ASTCENC_PRF_LDR_SRGB = 0,
	/** @brief The LDR linear color profile. */
	ASTCENC_PRF_LDR,
	/** @brief The HDR RGB with LDR alpha color profile. */
	ASTCENC_PRF_HDR_RGB_LDR_A,
	/** @brief The HDR RGBA color profile. */
	ASTCENC_PRF_HDR
};

/**
 * @brief A codec quality preset.
 */
enum astcenc_preset {
	/** @brief The fast, lowest quality, search preset. */
	ASTCENC_PRE_FAST = 0,
	/** @brief The medium quality search preset. */
	ASTCENC_PRE_MEDIUM,
	/** @brief The throrough quality search preset. */
	ASTCENC_PRE_THOROUGH,
	/** @brief The exhaustive quality search preset. */
	ASTCENC_PRE_EXHAUSTIVE
};

/**
 * @brief A codec channel swizzle selector.
 */
enum astcenc_swz {
	/** @brief Select the red channel. */
	ASTCENC_SWZ_R = 0,
	/** @brief Select the green channel. */
	ASTCENC_SWZ_G = 1,
	/** @brief Select the blue channel. */
	ASTCENC_SWZ_B = 2,
	/** @brief Select the alpha channel. */
	ASTCENC_SWZ_A = 3,
	/** @brief Use a constant zero channel. */
	ASTCENC_SWZ_0 = 4,
	/** @brief Use a constant one channel. */
	ASTCENC_SWZ_1 = 5,
	/** @brief Use a reconstructed normal vector Z channel. */
	ASTCENC_SWZ_Z = 6
};

/**
 * @brief A texel channel swizzle.
 */
struct astcenc_swizzle {
	/** @brief The red channel selector. */
	astcenc_swz r;
	/** @brief The green channel selector. */
	astcenc_swz g;
	/** @brief The blue channel selector. */
	astcenc_swz b;
	/** @brief The alpha channel selector. */
	astcenc_swz a;
};

/**
 * @brief A texel channel data format.
 */
enum astcenc_type {
	/** @brief Unorm 8-bit data per channel. */
	ASTCENC_TYPE_U8 = 0,
	/** @brief 16-bit float per channel. */
	ASTCENC_TYPE_F16 = 1,
	/** @brief 32-bit float per channel. */
	ASTCENC_TYPE_F32 = 2
};

/**
 * @brief Enable normal map compression.
 *
 * Input data will be treated a two channel normal map, storing X and Y, and
 * the codec will optimize for angular error rather than simple linear PSNR.
 * In this mode the input swizzle should be e.g. rrrg (the default ordering for
 * ASTC normals on the command line) or gggr (the ordering used by BC5).
 */
static const unsigned int ASTCENC_FLG_MAP_NORMAL          = 1 << 0;

/**
 * @brief Enable mask map compression.
 *
 * Input data will be treated a multi-layer mask map, where is is desirable for
 * the color channels to be treated independently for the purposes of error
 * analysis.
 */
static const unsigned int ASTCENC_FLG_MAP_MASK            = 1 << 1;

/**
 * @brief Enable alpha weighting.
 *
 * The input alpha value is used for transparency, so errors in the RGB
 * channels are weighted by the transparency level. This allows the codec to
 * more accurately encode the alpha value in areas where the color value
 * is less significant.
 */
static const unsigned int ASTCENC_FLG_USE_ALPHA_WEIGHT    = 1 << 2;

/**
 * @brief Enable perceptual error metrics.
 *
 * This mode enables perceptual compression mode, which will optimize for
 * perceptual error rather than best PSNR. Only some input modes support
 * perceptual error metrics.
 */
static const unsigned int ASTCENC_FLG_USE_PERCEPTUAL      = 1 << 3;

/**
 * @brief Create a decompression-only context.
 *
 * This mode enables context allocation to skip some transient buffer
 * allocation, resulting in a lower-memory footprint.
 */
static const unsigned int ASTCENC_FLG_DECOMPRESS_ONLY     = 1 << 4;

/**
 * @brief The bit mask of all valid flags.
 */
static const unsigned int ASTCENC_ALL_FLAGS =
                              ASTCENC_FLG_MAP_NORMAL |
                              ASTCENC_FLG_MAP_MASK |
                              ASTCENC_FLG_USE_ALPHA_WEIGHT |
                              ASTCENC_FLG_USE_PERCEPTUAL |
                              ASTCENC_FLG_DECOMPRESS_ONLY;

/**
 * @brief The config structure.
 *
 * This structure will initially be populated by a call to astcenc_config_init,
 * but power users may modify it before calling astcenc_context_alloc. See
 * astcenccli_toplevel_help.cpp for full user documentation of the power-user
 * settings.
 *
 * Note for any settings which are associated with a specific color channel,
 * the value in the config applies to the channel that exists after any
 * compression data swizzle is applied.
 */
struct astcenc_config {
	/** @brief The color profile. */
	astcenc_profile profile;

	/** @brief The set of set flags. */
	unsigned int flags;

	/** @brief The ASTC block size X dimension. */
	unsigned int block_x;

	/** @brief The ASTC block size Y dimension. */
	unsigned int block_y;

	/** @brief The ASTC block size Z dimension. */
	unsigned int block_z;

	/** @brief The size of the texel kernel for error weighting (-v). */
	unsigned int v_rgba_radius;

	/** @brief The mean and stdev channel mix for error weighting (-v). */
	float v_rgba_mean_stdev_mix;

	/** @brief The texel RGB power for error weighting (-v). */
	float v_rgb_power;

	/** @brief The texel RGB base weight for error weighting (-v). */
	float v_rgb_base;

	/** @brief The texel RGB mean weight for error weighting (-v). */
	float v_rgb_mean;

	/** @brief The texel RGB stdev for error weighting (-v). */
	float v_rgb_stdev;

	/** @brief The texel A power for error weighting (-va). */
	float v_a_power;

	/** @brief The texel A base weight for error weighting (-va). */
	float v_a_base;

	/** @brief The texel A mean weight for error weighting (-va). */
	float v_a_mean;

	/** @brief The texel A stdev for error weighting (-va). */
	float v_a_stdev;

	/** @brief The red channel weight scale for error weighting (-cw). */
	float cw_r_weight;

	/** @brief The green channel weight scale for error weighting (-cw). */
	float cw_g_weight;

	/** @brief The blue channel weight scale for error weighting (-cw). */
	float cw_b_weight;

	/** @brief The alpha channel weight scale for error weighting (-cw). */
	float cw_a_weight;

	/**
	 * @brief The radius for any alpha-weight scaling (-a).
	 *
	 * It is recommended that this is set to 1 when using FLG_USE_ALPHA_WEIGHT
	 * on a texture that will be sampled using linear texture filtering to
	 * minimize color bleed out of transparent texels that are adjcent to
	 * non-transparent texels.
	 */
	unsigned int a_scale_radius;

	/**
	 * @brief The additional weight for block edge texels (-b).
	 *
	 * This is generic tool for reducing artefacts visible on block changes.
	 */
	float b_deblock_weight;

	/**
	 * @brief The maximum number of partitions searched (-partitionlimit).
	 *
	 * Valid values are between 1 and 1024.
	 */
	unsigned int tune_partition_limit;

	/**
	 * @brief The maximum centile for block modes searched (-blockmodelimit).
	 *
	 * Valid values are between 1 and 100.
	 */
	unsigned int tune_block_mode_limit;

	/**
	 * @brief The maximum iterative refinements applied (-refinementlimit).
	 *
	 * Valid values are between 1 and N; there is no technical upper limit
	 * but little benefit is expected after N=4.
	 */
	unsigned int tune_refinement_limit;

	/**
	 * @brief The dB threshold for stopping block search (-dblimit).
	 *
	 * This option is ineffective for HDR textures.
	 */
	float tune_db_limit;

	/**
	 * @brief The threshold for skipping 3+ partitions (-partitionearlylimit).
	 *
	 * This option is ineffective for normal maps.
	 */
	float tune_partition_early_out_limit;

	/**
	 * @brief The threshold for skipping 2 weight planess (-planecorlimit).
	 *
	 * This option is ineffective for normal maps.
	 */
	float tune_two_plane_early_out_limit;
};

/**
 * @brief An uncompressed 2D or 3D image.
 *
 * Inputs can be either 8-bit unorm inputs (passed in via the data8 pointer),
 * or 16-bit floating point inputs (passed in via the data16 pointer). The
 * unused pointer must be set to nullptr. Data is always passed in as 4 color
 * channels, and accessed as 3D array indexed using [Z][Y][(X * 4) + (0..3)].
 */
struct astcenc_image {
	/** @brief The X dimension of the image, in texels. */
	unsigned int dim_x;
	/** @brief The Y dimension of the image, in texels. */
	unsigned int dim_y;
	/** @brief The X dimension of the image, in texels. */
	unsigned int dim_z;
	/** @brief The border padding dimensions, in texels. */
	unsigned int dim_pad;
	/** @brief The data type per channel. */
	astcenc_type data_type;
	/** @brief The data; actually of type <t>***. */
	void *data;
};

/**
 * Populate a codec config based on default settings.
 *
 * Power users can edit the returned config struct to apply manual fine tuning
 * before allocating the context.
 *
 * @param      profile   Color profile.
 * @param      block_x   ASTC block size X dimension.
 * @param      block_y   ASTC block size Y dimension.
 * @param      block_z   ASTC block size Z dimension.
 * @param      preset    Search quality preset.
 * @param      flags     A valid set of ASTCENC_FLG_* flag bits.
 * @param[out] config    Output config struct to populate.
 *
 * @return ASTCENC_SUCCESS on success, or an error if the inputs are invalid
 * either individually, or in combination.
 */
astcenc_error astcenc_config_init(
	astcenc_profile profile,
	unsigned int block_x,
	unsigned int block_y,
	unsigned int block_z,
	astcenc_preset preset,
	unsigned int flags,
	astcenc_config& config);

/**
 * @brief Allocate a new codec context based on a config.
 *
 * This function allocates all of the memory resources and threads needed by
 * the codec. This can be slow, so it is recommended that contexts are reused
 * to serially compress or decompress multiple images to amortize setup cost.
 *
 * Contexts can be allocated to support only decompression by setting the
 * ASTCENC_FLG_DECOMPRESS_ONLY flag when creating the configuration. These
 * contexts must be allocated with a thread count of 1 (decompression is always
 * single threaded), and the compression functions will fail if invoked. For
 * a decompress-only library build the ASTCENC_FLG_DECOMPRESS_ONLY flag must
 * be set when creating ay context.
 *
 * @param[in]  config         Codec config.
 * @param      thread_count   Thread count to configure for. Decompress-only
 *                            contexts must have a thread_count of 1.
 * @param[out] context        Location to store an opaque context pointer.
 *
 * @return ASTCENC_SUCCESS on success, or an error if context creation failed.
 */
astcenc_error astcenc_context_alloc(
	const astcenc_config& config,
	unsigned int thread_count,
	astcenc_context** context);

/**
 * @brief Compress an image.
 *
 * A single context can only compress or decompress a single image at a time.
 *
 * For a context configured for multi-threading, any set of the N threads can
 * call this function. Work will be dynamically scheduled across the threads
 * available. Each thread must have a unique thread_index.
 *
 * @param         context        Codec context.
 * @param[in,out] image          Input image.
 * @param         swizzle        Compression data swizzle.
 * @param[out]    data_out       Pointer to output data array.
 * @param         data_len       Length of the output data array.
 * @param         thread_index   Thread index [0..N-1] of calling thread.
 *
 * @return ASTCENC_SUCCESS on success, or an error if compression failed.
 */
astcenc_error astcenc_compress_image(
	astcenc_context* context,
	astcenc_image& image,
	astcenc_swizzle swizzle,
	uint8_t* data_out,
	size_t data_len,
	unsigned int thread_index);

/**
 * @brief Reset the compressor state for a new compression.
 *
 * The caller is responsible for synchronizing threads in the worker thread
 * pool. This function must only be called when all threads have exited the
 * astcenc_compress_image() function for image N, but before any thread enters
 * it for image N + 1.
 *
 * @param context   Codec context.
 *
 * @return ASTCENC_SUCCESS on success, or an error if reset failed.
 */
astcenc_error astcenc_compress_reset(
	astcenc_context* context);

/**
 * @brief Decompress an image.
 *
 * @param         context      Codec context.
 * @param[in]     data         Pointer to compressed data.
 * @param         data_len     Length of the compressed data, in bytes.
 * @param[in,out] image_out    Output image.
 * @param         swizzle      Decompression data swizzle.
 *
 * @return ASTCENC_SUCCESS on success, or an error if decompression failed.
 */
astcenc_error astcenc_decompress_image(
	astcenc_context* context,
	const uint8_t* data,
	size_t data_len,
	astcenc_image& image_out,
	astcenc_swizzle swizzle);

/**
 * Free the compressor context.
 *
 * @param context   The codec context.
 */
void astcenc_context_free(
	astcenc_context* context);

/**
 * @brief Get a printable string for specific status code.
 *
 * @param status   The status value.
 *
 * @return A human readable nul-terminated string.
 */
const char* astcenc_get_error_string(
	astcenc_error status);

#endif
