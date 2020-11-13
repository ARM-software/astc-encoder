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

#ifndef ASTCENC_INTERNAL_INCLUDED
#define ASTCENC_INTERNAL_INCLUDED

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <condition_variable>
#include <functional>
#include <mutex>

#ifndef ASTCENC_SSE
#error ERROR: ASTCENC_SSE not defined
#endif

#ifndef ASTCENC_POPCNT
#error ERROR: ASTCENC_POPCNT not defined
#endif

#ifndef ASTCENC_AVX
#error ERROR: ASTCENC_AVX not defined
#endif

#ifndef ASTCENC_ISA_INVARIANCE
#error ERROR: ASTCENC_ISA_INVARIANCE not defined
#endif

#include "astcenc.h"
#include "astcenc_mathlib.h"

/* ============================================================================
  Constants
============================================================================ */
#define MAX_TEXELS_PER_BLOCK 216
#define MAX_WEIGHTS_PER_BLOCK 64
#define MIN_WEIGHT_BITS_PER_BLOCK 24
#define MAX_WEIGHT_BITS_PER_BLOCK 96
#define PARTITION_BITS 10
#define PARTITION_COUNT (1 << PARTITION_BITS)

// the sum of weights for one texel.
#define TEXEL_WEIGHT_SUM 16
#define MAX_DECIMATION_MODES 87
#define MAX_WEIGHT_MODES 2048

// A high default error value
static const float ERROR_CALC_DEFAULT { 1e30f };

/* ============================================================================
  Compile-time tuning parameters
============================================================================ */
// The max texel count in a block which can try the one partition fast path.
// Default: enabled for 4x4 and 5x4 blocks.
static const int TUNE_MAX_TEXELS_MODE0_FASTPATH { 24 };

// The maximum number of candidate encodings returned for each encoding mode.
// Default: depends on quality preset
static const int TUNE_MAX_TRIAL_CANDIDATES { 4 };

/* ============================================================================
  Other configuration parameters
============================================================================ */

// Uncomment to enable checking for inappropriate NaNs; Linux only, and slow!
// #define DEBUG_CAPTURE_NAN

/* ============================================================================
  Parallel execution control
============================================================================ */

/**
 * @brief A simple counter-based manager for parallel task execution.
 *
 * The task processing execution consists of:
 *
 *     * A single-threaded init stage.
 *     * A multi-threaded processing stage.
 *     * A condition variable so threads can wait for processing completion.
 *
 * The init stage will be executed by the first thread to arrive in the
 * critical section, there is no master thread in the thread pool.
 *
 * The processing stage uses dynamic dispatch to assign task tickets to threads
 * on an on-demand basis. Threads may each therefore executed different numbers
 * of tasks, depending on their processing complexity. The task queue and the
 * task tickets are just counters; the caller must map these integers to an
 * actual processing partition in a specific problem domain.
 *
 * The exit wait condition is needed to ensure processing has finished before
 * a worker thread can progress to the next stage of the pipeline. Specifically
 * a worker may exit the processing stage because there are no new tasks to
 * assign to it while other worker threads are still processing. Calling wait()
 * will ensure that all other worker have finished before the thread can
 * proceed.
 *
 * The basic usage model:
 *
 *     // --------- From single-threaded code ---------
 *
 *     // Reset the tracker state
 *     manager->reset()
 *
 *     // --------- From multi-threaded code ---------
 *
 *     // Run the stage init; only first thread actually runs the lambda
 *     manager->init(<lambda>)
 *
 *     do
 *     {
 *         // Request a task assignment
 *         uint task_count;
 *         uint base_index = manager->get_tasks(<granule>, task_count);
 *
 *         // Process any tasks we were given (task_count <= granule size)
 *         if (task_count)
 *         {
 *             // Run the user task processing code for N tasks here
 *             ...
 *
 *             // Flag these tasks as complete
 *             manager->complete_tasks(task_count);
 *         }
 *     } while (task_count);
 *
 *     // Wait for all threads to complete tasks before progressing
 *     manager->wait()
 *
  *     // Run the stage term; only first thread actually runs the lambda
 *     manager->term(<lambda>)
 */
class ParallelManager
{
private:
	/** \brief Lock used for critical section and condition synchronization. */
	std::mutex m_lock;

	/** \brief True if the stage init() step has been executed. */
	bool m_init_done;

	/** \brief True if the stage term() step has been executed. */
	bool m_term_done;

	/** \brief Contition variable for tracking stage processing completion. */
	std::condition_variable m_complete;

	/** \brief Number of tasks started, but not necessarily finished. */
	unsigned int m_start_count;

	/** \brief Number of tasks finished. */
	unsigned int m_done_count;

	/** \brief Number of tasks that need to be processed. */
	unsigned int m_task_count;

public:
	/** \brief Create a new ParallelManager. */
	ParallelManager()
	{
		reset();
	}

	/**
	 * \brief Reset the tracker for a new processing batch.
	 *
	 * This must be called from single-threaded code before starting the
	 * multi-threaded procesing operations.
	 */
	void reset()
	{
		m_init_done = false;
		m_term_done = false;
		m_start_count = 0;
		m_done_count = 0;
		m_task_count = 0;
	}

	/**
	 * \brief Trigger the pipeline stage init step.
	 *
	 * This can be called from multi-threaded code. The first thread to
	 * hit this will process the initialization. Other threads will block
	 * and wait for it to complete.
	 *
	 * \param init_func    Callable which executes the stage initialization.
	 *                     Must return the number of tasks in the stage.
	 */
	void init(std::function<unsigned int(void)> init_func)
	{
		std::lock_guard<std::mutex> lck(m_lock);
		if (!m_init_done)
		{
			m_task_count = init_func();
			m_init_done = true;
		}
	}

	/**
	 * \brief Trigger the pipeline stage init step.
	 *
	 * This can be called from multi-threaded code. The first thread to
	 * hit this will process the initialization. Other threads will block
	 * and wait for it to complete.
	 *
	 * \param task_count   Total number of tasks needing processing.
	 */
	void init(unsigned int task_count)
	{
		std::lock_guard<std::mutex> lck(m_lock);
		if (!m_init_done)
		{
			m_task_count = task_count;
			m_init_done = true;
		}
	}

	/**
	 * \brief Request a task assignment.
	 *
	 * Assign up to \c granule tasks to the caller for processing.
	 *
	 * \param      granule   Maximum number of tasks that can be assigned.
	 * \param[out] count     Actual number of tasks assigned, or zero if
	 *                       no tasks were assigned.
	 *
	 * \return Task index of the first assigned task; assigned tasks
	 *         increment from this.
	 */
	unsigned int get_task_assignment(unsigned int granule, unsigned int& count)
	{
		std::lock_guard<std::mutex> lck(m_lock);
		unsigned int base = m_start_count;
		count = std::min(granule, m_task_count - m_start_count);
		m_start_count += count;
		return base;
	}

	/**
	 * \brief Complete a task assignment.
	 *
	 * Mark \c count tasks as complete. This will notify all threads blocked
	 * on \c wait() if this completes the processing of the stage.
	 *
	 * \param count   The number of completed tasks.
	 */
	void complete_task_assignment(unsigned int count)
	{
		std::unique_lock<std::mutex> lck(m_lock);
		this->m_done_count += count;
		if (m_done_count == m_task_count)
		{
			lck.unlock();
			m_complete.notify_all();
		}
	}

	/**
	 * \brief Wait for stage processing to complete.
	 */
	void wait()
	{
		std::unique_lock<std::mutex> lck(m_lock);
		m_complete.wait(lck, [this]{ return m_done_count == m_task_count; });
	}

	/**
	 * \brief Trigger the pipeline stage term step.
	 *
	 * This can be called from multi-threaded code. The first thread to
	 * hit this will process the thread termintion. Caller must have called
	 * wait() prior to calling this function to ensure processing is complete.
	 *
	 * \param term_func   Callable which executes the stage termination.
	 */
	void term(std::function<void(void)> term_func)
	{
		std::lock_guard<std::mutex> lck(m_lock);
		if (!m_term_done)
		{
			term_func();
			m_term_done = true;
		}
	}
};


/*
	Partition table representation:
	For each block size, we have 3 tables, each with 1024 partitionings;
	these three tables correspond to 2, 3 and 4 partitions respectively.
	For each partitioning, we have:
	* a 4-entry table indicating how many texels there are in each of the 4 partitions.
	  This may be from 0 to a very large value.
	* a table indicating the partition index of each of the texels in the block.
	  Each index may be 0, 1, 2 or 3.
	* Each element in the table is an uint8_t indicating partition index (0, 1, 2 or 3)
*/

struct partition_info
{
	int partition_count;
	uint8_t texels_per_partition[4];
	uint8_t partition_of_texel[MAX_TEXELS_PER_BLOCK];
	uint8_t texels_of_partition[4][MAX_TEXELS_PER_BLOCK];
	uint64_t coverage_bitmaps[4];
};

/*
   In ASTC, we don't necessarily provide a weight for every texel.
   As such, for each block size, there are a number of patterns where some texels
   have their weights computed as a weighted average of more than 1 weight.
   As such, the codec uses a data structure that tells us: for each texel, which
   weights it is a combination of for each weight, which texels it contributes to.
   The decimation_table is this data structure.
*/
struct decimation_table
{
	int num_texels;
	int num_weights;
	uint8_t texel_num_weights[MAX_TEXELS_PER_BLOCK];	// number of indices that go into the calculation for a texel
	uint8_t texel_weights_int[MAX_TEXELS_PER_BLOCK][4];	// the weight to assign to each weight
	float texel_weights_float[MAX_TEXELS_PER_BLOCK][4];	// the weight to assign to each weight
	uint8_t texel_weights[MAX_TEXELS_PER_BLOCK][4];	// the weights that go into a texel calculation
	uint8_t weight_num_texels[MAX_WEIGHTS_PER_BLOCK];	// the number of texels that a given weight contributes to
	uint8_t weight_texel[MAX_WEIGHTS_PER_BLOCK][MAX_TEXELS_PER_BLOCK];	// the texels that the weight contributes to
	uint8_t weights_int[MAX_WEIGHTS_PER_BLOCK][MAX_TEXELS_PER_BLOCK];	// the weights that the weight contributes to a texel.
	float weights_flt[MAX_WEIGHTS_PER_BLOCK][MAX_TEXELS_PER_BLOCK];	// the weights that the weight contributes to a texel.

	// folded data structures:
	//  * texel_weights_texel[i][j] = texel_weights[weight_texel[i][j]];
	//  * texel_weights_float_texel[i][j] = texel_weights_float[weight_texel[i][j]
	uint8_t texel_weights_texel[MAX_WEIGHTS_PER_BLOCK][MAX_TEXELS_PER_BLOCK][4];
	float texel_weights_float_texel[MAX_WEIGHTS_PER_BLOCK][MAX_TEXELS_PER_BLOCK][4];
};

/*
   data structure describing information that pertains to a block size and its associated block modes.
*/
struct block_mode
{
	int8_t decimation_mode;
	int8_t quantization_mode;
	int8_t is_dual_plane;
	int16_t mode_index;
	float percentile;
};

struct block_size_descriptor
{
	int xdim;
	int ydim;
	int zdim;
	int texel_count;

	int decimation_mode_count;
	int decimation_mode_samples[MAX_DECIMATION_MODES];
	int decimation_mode_maxprec_1plane[MAX_DECIMATION_MODES];
	int decimation_mode_maxprec_2planes[MAX_DECIMATION_MODES];
	float decimation_mode_percentile[MAX_DECIMATION_MODES];
	int permit_encode[MAX_DECIMATION_MODES];
	const decimation_table *decimation_tables[MAX_DECIMATION_MODES];

	// out of all possible 2048 weight modes, only a subset is
	// actually valid for the current configuration (e.g. 6x6
	// 2D LDR has 370 valid modes); the valid ones are packed into
	// block_modes_packed array.
	block_mode block_modes_packed[MAX_WEIGHT_MODES];
	int block_mode_packed_count;
	// get index of block mode inside the block_modes_packed array,
	// or -1 if mode is not valid for the current configuration.
	int16_t block_mode_to_packed[MAX_WEIGHT_MODES];

	// for the k-means bed bitmap partitioning algorithm, we don't
	// want to consider more than 64 texels; this array specifies
	// which 64 texels (if that many) to consider.
	int texelcount_for_bitmap_partitioning;
	int texels_for_bitmap_partitioning[64];

	// All the partitioning information for this block size
	partition_info partitions[(3*PARTITION_COUNT)+1];
};

// data structure representing one block of an image.
// it is expanded to float prior to processing to save some computation time
// on conversions to/from uint8_t (this also allows us to handle HDR textures easily)
struct imageblock
{
	float data_r[MAX_TEXELS_PER_BLOCK];  // the data that we will compress, either linear or LNS (0..65535 in both cases)
	float data_g[MAX_TEXELS_PER_BLOCK];
	float data_b[MAX_TEXELS_PER_BLOCK];
	float data_a[MAX_TEXELS_PER_BLOCK];
	float4 origin_texel;

	uint8_t rgb_lns[MAX_TEXELS_PER_BLOCK];      // 1 if RGB data are being treated as LNS
	uint8_t alpha_lns[MAX_TEXELS_PER_BLOCK];    // 1 if Alpha data are being treated as LNS
	uint8_t nan_texel[MAX_TEXELS_PER_BLOCK];    // 1 if the texel is a NaN-texel.

	float red_min, red_max;
	float green_min, green_max;
	float blue_min, blue_max;
	float alpha_min, alpha_max;
	int grayscale;				// 1 if R=G=B for every pixel, 0 otherwise

	int xpos, ypos, zpos;
};

static inline int imageblock_uses_alpha(const imageblock * pb)
{
	return pb->alpha_max != pb->alpha_min;
}

void update_imageblock_flags(
	imageblock* pb,
	int xdim,
	int ydim,
	int zdim);

void imageblock_initialize_orig_from_work(
	imageblock * pb,
	int pixelcount);

void imageblock_initialize_work_from_orig(
	imageblock * pb,
	int pixelcount);

/*
	Data structure representing error weighting for one block of an image. this is used as
	a multiplier for the error weight to apply to each color component when computing PSNR.

	This weighting has several uses: it's usable for RA, GA, BA, A weighting, which is useful
	for alpha-textures it's usable for HDR textures, where weighting should be approximately inverse to
	luminance it's usable for perceptual weighting, where we assign higher weight to low-variability
	regions than to high-variability regions. it's usable for suppressing off-edge block content in
	case the texture doesn't actually extend to the edge of the block.

	For the default case (everything is evenly weighted), every weight is 1. For the RA,GA,BA,A case,
	we multiply the R,G,B weights with that of the alpha.

	Putting the same weight in every component should result in the default case.
	The following relations should hold:

	texel_weight_rg[i] = (texel_weight_r[i] + texel_weight_g[i]) / 2
	texel_weight_lum[i] = (texel_weight_r[i] + texel_weight_g[i] + texel_weight_b[i]) / 3
	texel_weight[i] = (texel_weight_r[i] + texel_weight_g[i] + texel_weight_b[i] + texel_weight_a[i] / 4
 */

struct error_weight_block
{
	float4 error_weights[MAX_TEXELS_PER_BLOCK];
	float texel_weight[MAX_TEXELS_PER_BLOCK];
	float texel_weight_gba[MAX_TEXELS_PER_BLOCK];
	float texel_weight_rba[MAX_TEXELS_PER_BLOCK];
	float texel_weight_rga[MAX_TEXELS_PER_BLOCK];
	float texel_weight_rgb[MAX_TEXELS_PER_BLOCK];

	float texel_weight_rg[MAX_TEXELS_PER_BLOCK];
	float texel_weight_rb[MAX_TEXELS_PER_BLOCK];
	float texel_weight_gb[MAX_TEXELS_PER_BLOCK];
	float texel_weight_ra[MAX_TEXELS_PER_BLOCK];

	float texel_weight_r[MAX_TEXELS_PER_BLOCK];
	float texel_weight_g[MAX_TEXELS_PER_BLOCK];
	float texel_weight_b[MAX_TEXELS_PER_BLOCK];
	float texel_weight_a[MAX_TEXELS_PER_BLOCK];

	int contains_zeroweight_texels;
};

// enumeration of all the quantization methods we support under this format.
enum quantization_method
{
	QUANT_2 = 0,
	QUANT_3 = 1,
	QUANT_4 = 2,
	QUANT_5 = 3,
	QUANT_6 = 4,
	QUANT_8 = 5,
	QUANT_10 = 6,
	QUANT_12 = 7,
	QUANT_16 = 8,
	QUANT_20 = 9,
	QUANT_24 = 10,
	QUANT_32 = 11,
	QUANT_40 = 12,
	QUANT_48 = 13,
	QUANT_64 = 14,
	QUANT_80 = 15,
	QUANT_96 = 16,
	QUANT_128 = 17,
	QUANT_160 = 18,
	QUANT_192 = 19,
	QUANT_256 = 20
};

/**
 * @brief Weight quantization transfer table.
 *
 * ASTC can store texel weights at many quantization levels, so for performance
 * we store essential information about each level as a precomputed data
 * structure.
 *
 * Unquantized weights are integers in the range [0, 64], or floats [0, 1].
 *
 * This structure provides the following information:
 * A table, used to estimate the closest quantized
	weight for a given floating-point weight. For each quantized weight, the corresponding unquantized
	and floating-point values. For each quantized weight, a previous-value and a next-value.
*/
struct quantization_and_transfer_table
{
	/** The quantization level used */
	quantization_method method;
	/** The unscrambled unquantized value. */
	// TODO: Converted to floats to support AVX gathers
	float unquantized_value_unsc[33];
	/** The scrambling order: value[map[i]] == value_unsc[i] */
	// TODO: Converted to u32 to support AVX gathers
	int32_t scramble_map[32];
	/** The scrambled unquantized values. */
	uint8_t unquantized_value[32];
	/**
	 * An encoded table of previous-and-next weight values, indexed by the
	 * current unquantized value.
	 *  * bits 7:0 = previous-index, unquantized
	 *  * bits 15:8 = next-index, unquantized
	 *  * bits 23:16 = previous-index, quantized
	 *  * bits 31:24 = next-index, quantized
	 */
	uint32_t prev_next_values[65];
};

extern const quantization_and_transfer_table quant_and_xfer_tables[12];

enum endpoint_formats
{
	FMT_LUMINANCE = 0,
	FMT_LUMINANCE_DELTA = 1,
	FMT_HDR_LUMINANCE_LARGE_RANGE = 2,
	FMT_HDR_LUMINANCE_SMALL_RANGE = 3,
	FMT_LUMINANCE_ALPHA = 4,
	FMT_LUMINANCE_ALPHA_DELTA = 5,
	FMT_RGB_SCALE = 6,
	FMT_HDR_RGB_SCALE = 7,
	FMT_RGB = 8,
	FMT_RGB_DELTA = 9,
	FMT_RGB_SCALE_ALPHA = 10,
	FMT_HDR_RGB = 11,
	FMT_RGBA = 12,
	FMT_RGBA_DELTA = 13,
	FMT_HDR_RGB_LDR_ALPHA = 14,
	FMT_HDR_RGBA = 15
};

struct symbolic_compressed_block
{
	int error_block;			// 1 marks error block, 0 marks non-error-block.
	int block_mode;				// 0 to 2047. Negative value marks constant-color block (-1: FP16, -2:UINT16)
	int partition_count;		// 1 to 4; Zero marks a constant-color block.
	int partition_index;		// 0 to 1023
	int color_formats[4];		// color format for each endpoint color pair.
	int color_formats_matched;	// color format for all endpoint pairs are matched.
	int color_values[4][12];	// quantized endpoint color pairs.
	int color_quantization_level;
	uint8_t plane1_weights[MAX_WEIGHTS_PER_BLOCK];	// quantized and decimated weights
	uint8_t plane2_weights[MAX_WEIGHTS_PER_BLOCK];
	int plane2_color_component;	// color component for the secondary plane of weights
	int constant_color[4];		// constant-color, as FP16 or UINT16. Used for constant-color blocks only.
};

struct physical_compressed_block
{
	uint8_t data[16];
};

/* ============================================================================
  Functions and data pertaining to quantization and encoding
============================================================================ */

/**
 * @brief Populate the blocksize descriptor for the target block size.
 *
 * This will also initialize the partition table metadata, which is stored
 * as part of the BSD structure.
 *
 * @param xdim The x axis size of the block.
 * @param ydim The y axis size of the block.
 * @param zdim The z axis size of the block.
 * @param bsd  The structure to populate.
 */
void init_block_size_descriptor(
	int xdim,
	int ydim,
	int zdim,
	block_size_descriptor* bsd);

void term_block_size_descriptor(
	block_size_descriptor* bsd);

/**
 * @brief Populate the partition tables for the target block size.
 *
 * Note the block_size_size descriptor must be initialized before calling this
 * function.
 *
 * @param bsd  The structure to populate.
 */
void init_partition_tables(
	block_size_descriptor* bsd);

static inline const partition_info *get_partition_table(
	const block_size_descriptor* bsd,
	int partition_count
) {
	if (partition_count == 1) {
		partition_count = 5;
	}
	int index = (partition_count - 2) * PARTITION_COUNT;
	return bsd->partitions + index;
}

/**
 * @brief Get the percentile table for 2D block modes.
 *
 * This is an empirically determined prioritization of which block modes to
 * use in the search in terms of their centile (lower centiles = more useful).
 *
 * Returns a dynamically allocated array; caller must free with delete[].
 *
 * @param xdim The block x size.
 * @param ydim The block y size.
 *
 * @return The unpacked table.
 */
const float *get_2d_percentile_table(
	int xdim,
	int ydim);

/**
 * @brief Query if a 2D block size is legal.
 *
 * @return A non-zero value if legal, zero otherwise.
 */
int is_legal_2d_block_size(
	int xdim,
	int ydim);

/**
 * @brief Query if a 3D block size is legal.
 *
 * @return A non-zero value if legal, zero otherwise.
 */
int is_legal_3d_block_size(
	int xdim,
	int ydim,
	int zdim);

// ***********************************************************
// functions and data pertaining to quantization and encoding
// **********************************************************

extern const uint8_t color_quantization_tables[21][256];
extern const uint8_t color_unquantization_tables[21][256];
extern int quantization_mode_table[17][128];

void encode_ise(
	int quantization_level,
	int elements,
	const uint8_t* input_data,
	uint8_t* output_data,
	int bit_offset);

void decode_ise(
	int quantization_level,
	int elements,
	const uint8_t* input_data,
	uint8_t* output_data,
	int bit_offset);

int compute_ise_bitcount(
	int items,
	quantization_method quant);

void build_quantization_mode_table(void);

// **********************************************
// functions and data pertaining to partitioning
// **********************************************

// functions to compute color averages and dominant directions
// for each partition in a block

void compute_averages_and_directions_rgb(
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	const float4* color_scalefactors,
	float3* averages,
	float3* directions_rgb);

void compute_averages_and_directions_rgba(
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	const float4* color_scalefactors,
	float4* averages,
	float4* directions_rgba);

void compute_averages_and_directions_3_components(
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	const float3 * color_scalefactors,
	int omittedComponent,
	float3* averages,
	float3* directions);

void compute_averages_and_directions_2_components(
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	const float2* color_scalefactors,
	int component1,
	int component2,
	float2* averages,
	float2* directions);

void compute_error_squared_rgba(
	const partition_info* pt,	// the partition that we use when computing the squared-error.
	const imageblock* blk,
	const error_weight_block* ewb,
	const processed_line4* plines_uncorr,
	const processed_line4* plines_samechroma,
	const processed_line3* plines_separate_red,
	const processed_line3* plines_separate_green,
	const processed_line3* plines_separate_blue,
	const processed_line3* plines_separate_alpha,
	float* length_uncorr,
	float* length_samechroma,
	float4* length_separate,
	float* uncorr_error,
	float* samechroma_error,
	float4* separate_color_error);

void compute_error_squared_rgb(
	const partition_info* pt,	// the partition that we use when computing the squared-error.
	const imageblock* blk,
	const error_weight_block* ewb,
	const processed_line3* plines_uncorr,
	const processed_line3* plines_samechroma,
	const processed_line2* plines_separate_red,
	const processed_line2* plines_separate_green,
	const processed_line2* plines_separate_blue,
	float* length_uncorr,
	float* length_samechroma,
	float3* length_separate,
	float* uncorr_error,
	float* samechroma_error,
	float3* separate_color_error);

// functions to compute error value across a tile for a particular line function
// for a single partition.
float compute_error_squared_rgb_single_partition(
	int partition_to_test,
	const block_size_descriptor* bsd,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	const processed_line3* lin	// the line for the partition.
);

// for each partition, compute its color weightings.
void compute_partition_error_color_weightings(
	const block_size_descriptor* bsd,
	const error_weight_block * ewb,
	const partition_info* pi,
	float4 error_weightings[4],
	float4 color_scalefactors[4]);

/**
 * \brief Find the best set of partitions to trial for a given block.
 *
 * On return \c best_partition_uncorrelated contains the best partition
 * assuming the data has noncorrelated chroma, \c best_partition_samechroma
 * contains the best partition assuming the data has corelated chroma, and
 * \c best_partition_dualplane contains the best partition assuming the data
 * has one uncorrelated color component.
 *
 * \c best_partition_dualplane is stored packed; bits [9:0] contain the
 * best partition, bits [11:10] contain the best color component.
 */
void find_best_partitionings(
	const block_size_descriptor* bsd,
	const imageblock* blk,
	const error_weight_block* ewb,
	int partition_count,
	int partition_search_limit,
	int* best_partition_uncorrelated,
	int* best_partition_samechroma,
	int* best_partition_dualplane);

// use k-means clustering to compute a partition ordering for a block.
void kmeans_compute_partition_ordering(
	const block_size_descriptor* bsd,
	int partition_count,
	const imageblock* blk,
	int *ordering);

// *********************************************************
// functions and data pertaining to images and imageblocks
// *********************************************************

/**
 * @brief Parameter structure for compute_pixel_region_variance().
 *
 * This function takes a structure to avoid spilling arguments to the stack
 * on every function invocation, as there are a lot of parameters.
 */
struct pixel_region_variance_args
{
	/** The image to analyze. */
	const astcenc_image* img;
	/** The RGB channel power adjustment. */
	float rgb_power;
	/** The alpha channel power adjustment. */
	float alpha_power;
	/** The channel swizzle pattern. */
	astcenc_swizzle swz;
	/** Should the algorithm bother with Z axis processing? */
	int have_z;
	/** The kernel radius for average and variance. */
	int avg_var_kernel_radius;
	/** The kernel radius for alpha processing. */
	int alpha_kernel_radius;
	/** The size of the working data to process. */
	int3 size;
	/** The position of first src and dst data in the data set. */
	int3 offset;
	/** The working memory buffer. */
	float4 *work_memory;
};

/**
 * @brief Parameter structure for compute_averages_and_variances_proc().
 */
struct avg_var_args
{
	/** The arguments for the nested variance computation. */
	pixel_region_variance_args arg;
	/** The image dimensions. */
	int3 img_size;
	/** The maximum working block dimensions. */
	int3 blk_size;
	/** The working block memory size. */
	int work_memory_size;
};

/**
 * @brief Compute regional averages and variances in an image.
 *
 * Results are written back into img->input_averages, img->input_variances,
 * and img->input_alpha_averages.
 *
 * @param img                   The input image data, also holds output data.
 * @param rgb_power             The RGB channel power.
 * @param alpha_power           The A channel power.
 * @param avg_var_kernel_radius The kernel radius (in pixels) for avg and var.
 * @param alpha_kernel_radius   The kernel radius (in pixels) for alpha mods.
 * @param swz                   Input data channel swizzle.
 * @param thread_count          The number of threads to use.
 *
 * @return The number of tasks in the processing stage.
 */
unsigned int init_compute_averages_and_variances(
	astcenc_image& img,
	float rgb_power,
	float alpha_power,
	int avg_var_kernel_radius,
	int alpha_kernel_radius,
	astcenc_swizzle swz,
	pixel_region_variance_args& arg,
	avg_var_args& ag);

void compute_averages_and_variances(
	astcenc_context& ctx,
	const avg_var_args& ag);

// fetch an image-block from the input file
void fetch_imageblock(
	astcenc_profile decode_mode,
	const astcenc_image& img,
	imageblock* pb,	// picture-block to initialize with image data
	const block_size_descriptor* bsd,
	// position in picture to fetch block from
	int xpos,
	int ypos,
	int zpos,
	astcenc_swizzle swz);

// write an image block to the output file buffer.
// the data written are taken from orig_data.
void write_imageblock(
	astcenc_image& img,
	const imageblock* pb,	// picture-block to initialize with image data
	const block_size_descriptor* bsd,
	// position in picture to write block to.
	int xpos,
	int ypos,
	int zpos,
	astcenc_swizzle swz);

// helper function to check whether a given picture-block has alpha that is not
// just uniformly 1.
int imageblock_uses_alpha(
	const imageblock * pb);

float compute_symbolic_block_difference(
	astcenc_profile decode_mode,
	const block_size_descriptor* bsd,
	const symbolic_compressed_block* scb,
	const imageblock* pb,
	const error_weight_block *ewb) ;

// ***********************************************************
// functions pertaining to computing texel weights for a block
// ***********************************************************
struct endpoints
{
	int partition_count;
	float4 endpt0[4];
	float4 endpt1[4];
};

struct endpoints_and_weights
{
	endpoints ep;
	float weights[MAX_TEXELS_PER_BLOCK];
	float weight_error_scale[MAX_TEXELS_PER_BLOCK];
};

void compute_endpoints_and_ideal_weights_1_plane(
	const block_size_descriptor* bsd,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	endpoints_and_weights* ei);

void compute_endpoints_and_ideal_weights_2_planes(
	const block_size_descriptor* bsd,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	int separate_component,
	endpoints_and_weights* ei1, // primary plane weights
	endpoints_and_weights* ei2); // secondary plane weights

void compute_ideal_weights_for_decimation_table(
	const endpoints_and_weights* eai,
	const decimation_table* it,
	float* weight_set,
	float* weights);

void compute_ideal_quantized_weights_for_decimation_table(
	const decimation_table* it,
	float low_bound,
	float high_bound,
	const float* weight_set_in,
	float* weight_set_out,
	uint8_t* quantized_weight_set,
	int quantization_level);

float compute_error_of_weight_set(
	const endpoints_and_weights* eai,
	const decimation_table* it,
	const float *weights);

void merge_endpoints(
	const endpoints* ep1,	// contains three of the color components
	const endpoints* ep2,	// contains the remaining color component
	int separate_component, endpoints* res);

// functions dealing with color endpoints

// function to pack a pair of color endpoints into a series of integers.
// the format used may or may not match the format specified;
// the return value is the format actually used.
int pack_color_endpoints(
	float4 color0,
	float4 color1,
	float4 rgbs_color,
	float4 rgbo_color,
	int format,
	int* output,
	int quantization_level);

// unpack a pair of color endpoints from a series of integers.
void unpack_color_endpoints(
	astcenc_profile decode_mode,
	int format,
	int quantization_level,
	const int* input,
	int* rgb_hdr,
	int* alpha_hdr,
	int* nan_endpoint,
	uint4* output0,
	uint4* output1);

struct encoding_choice_errors
{
	float rgb_scale_error;		// error of using LDR RGB-scale instead of complete endpoints.
	float rgb_luma_error;		// error of using HDR RGB-scale instead of complete endpoints.
	float luminance_error;		// error of using luminance instead of RGB
	float alpha_drop_error;		// error of discarding alpha
	float rgb_drop_error;		// error of discarding RGB
	int can_offset_encode;
	int can_blue_contract;
};

// buffers used to store intermediate data in compress_symbolic_block_fixed_partition_*()
struct alignas(ASTCENC_VECALIGN) compress_fixed_partition_buffers
{
	endpoints_and_weights ei1;
	endpoints_and_weights ei2;
	endpoints_and_weights eix1[MAX_DECIMATION_MODES];
	endpoints_and_weights eix2[MAX_DECIMATION_MODES];
	alignas(ASTCENC_VECALIGN) float decimated_quantized_weights[2 * MAX_DECIMATION_MODES * MAX_WEIGHTS_PER_BLOCK];
	alignas(ASTCENC_VECALIGN) float decimated_weights[2 * MAX_DECIMATION_MODES * MAX_WEIGHTS_PER_BLOCK];
	alignas(ASTCENC_VECALIGN) float flt_quantized_decimated_quantized_weights[2 * MAX_WEIGHT_MODES * MAX_WEIGHTS_PER_BLOCK];
	alignas(ASTCENC_VECALIGN) uint8_t u8_quantized_decimated_quantized_weights[2 * MAX_WEIGHT_MODES * MAX_WEIGHTS_PER_BLOCK];
};

struct compress_symbolic_block_buffers
{
	error_weight_block ewb;
	symbolic_compressed_block tempblocks[TUNE_MAX_TRIAL_CANDIDATES];
	compress_fixed_partition_buffers planes;
};

void compute_encoding_choice_errors(
	const block_size_descriptor* bsd,
	const imageblock* pb,
	const partition_info* pi,
	const error_weight_block* ewb,
	int separate_component,	// component that is separated out in 2-plane mode, -1 in 1-plane mode
	encoding_choice_errors* eci);

void determine_optimal_set_of_endpoint_formats_to_use(
	const block_size_descriptor* bsd,
	const partition_info* pt,
	const imageblock* blk,
	const error_weight_block* ewb,
	const endpoints* ep,
	int separate_component,	// separate color component for 2-plane mode; -1 for single-plane mode
	 // bitcounts and errors computed for the various quantization methods
	const int* qwt_bitcounts,
	const float* qwt_errors,
	int tune_candidate_limit,
	// output data
	int partition_format_specifiers[4][4],
	int quantized_weight[4],
	int quantization_level[4],
	int quantization_level_mod[4]);

void recompute_ideal_colors(
	int weight_quantization_mode,
	endpoints* ep,	// contains the endpoints we wish to update
	float4* rgbs_vectors,	// used to return RGBS-vectors for endpoint mode #6
	float4* rgbo_vectors,	// used to return RGBS-vectors for endpoint mode #7
	const uint8_t* weight_set8,	// the current set of weight values
	const uint8_t* plane2_weight_set8,	// nullptr if plane 2 is not actually used.
	int plane2_color_component,	// color component for 2nd plane of weights; -1 if the 2nd plane of weights is not present
	const partition_info* pi,
	const decimation_table* it,
	const imageblock* pb,	// picture-block containing the actual data.
	const error_weight_block* ewb);

void expand_deblock_weights(
	astcenc_context& ctx);

// functions pertaining to weight alignment
void prepare_angular_tables();

void imageblock_initialize_deriv(
	const imageblock* pb,
	int pixelcount,
	float4* dptr);

void compute_angular_endpoints_1plane(
	float mode_cutoff,
	const block_size_descriptor* bsd,
	const float* decimated_quantized_weights,
	const float* decimated_weights,
	float low_value[MAX_WEIGHT_MODES],
	float high_value[MAX_WEIGHT_MODES]);

void compute_angular_endpoints_2planes(
	float mode_cutoff,
	const block_size_descriptor * bsd,
	const float* decimated_quantized_weights,
	const float* decimated_weights,
	float low_value1[MAX_WEIGHT_MODES],
	float high_value1[MAX_WEIGHT_MODES],
	float low_value2[MAX_WEIGHT_MODES],
	float high_value2[MAX_WEIGHT_MODES]);

/* *********************************** high-level encode and decode functions ************************************ */

void compress_block(
	const astcenc_context& ctx,
	const astcenc_image& image,
	const imageblock* blk,
	symbolic_compressed_block& scb,
	physical_compressed_block& pcb,
	compress_symbolic_block_buffers* tmpbuf);

void decompress_symbolic_block(
	astcenc_profile decode_mode,
	const block_size_descriptor* bsd,
	int xpos,
	int ypos,
	int zpos,
	const symbolic_compressed_block* scb,
	imageblock* blk);

void symbolic_to_physical(
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& scb,
	physical_compressed_block& pcb);

void physical_to_symbolic(
	const block_size_descriptor& bsd,
	const physical_compressed_block& pcb,
	symbolic_compressed_block& scb);

uint16_t unorm16_to_sf16(
	uint16_t p);

struct astcenc_context
{
	astcenc_config config;
	unsigned int thread_count;
	block_size_descriptor* bsd;

	// Fields below here are not needed in a decompress-only build, but some
	// remain as they are small and it avoids littering the code with #ifdefs.
	// The most significant contributors to large structure size are omitted.

	// Regional average-and-variance information, initialized by
	// compute_averages_and_variances() only if the astc encoder
	// is requested to do error weighting based on averages and variances.
	float4 *input_averages;
	float4 *input_variances;
	float *input_alpha_averages;

	compress_symbolic_block_buffers* working_buffers;

#if !defined(ASTCENC_DECOMPRESS_ONLY)
	pixel_region_variance_args arg;
	avg_var_args ag;

	float deblock_weights[MAX_TEXELS_PER_BLOCK];

	ParallelManager manage_avg_var;
	ParallelManager manage_compress;
#endif
};

/* ============================================================================
  Platform-specific functions
============================================================================ */
/**
 * @brief Run-time detection if the host CPU supports SSE 4.2.
 * @returns Zero if not supported, positive value if it is.
 */
int cpu_supports_sse42();

/**
 * @brief Run-time detection if the host CPU supports popcnt.
 * @returns Zero if not supported, positive value if it is.
 */
int cpu_supports_popcnt();

/**
 * @brief Run-time detection if the host CPU supports avx2.
 * @returns Zero if not supported, positive value if it is.
 */
int cpu_supports_avx2();


/**
 * @brief Allocate an aligned memory buffer.
 *
 * Allocated memory must be freed by aligned_free;
 *
 * @param size    The desired buffer size.
 * @param align   The desired buffer alignment; must be 2^N.
 *
 * @returns The memory buffer pointer or nullptr on allocation failure.
 */
template<typename T>
T* aligned_malloc(size_t size, size_t align)
{
	void* ptr;
	int error = 0;

#if defined(_WIN32)
	ptr = _aligned_malloc(size, align);
#else
	error = posix_memalign(&ptr, align, size);
#endif

	if (error || (!ptr))
	{
		return nullptr;
	}

	return static_cast<T*>(ptr);
}

/**
 * @brief Free an aligned memory buffer.
 *
 * @param ptr   The buffer to free.
 */
template<typename T>
void aligned_free(T* ptr)
{
#if defined(_WIN32)
	_aligned_free(ptr);
#else
	free(ptr);
#endif
}

#endif
