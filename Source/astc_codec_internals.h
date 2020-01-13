// ----------------------------------------------------------------------------
//  This confidential and proprietary software may be used only as authorised
//  by a licensing agreement from Arm Limited.
//      (C) COPYRIGHT 2011-2020 Arm Limited, ALL RIGHTS RESERVED
//  The entire notice above must be reproduced on all authorised copies and
//  copies may only be made to the extent permitted by a licensing agreement
//  from Arm Limited.
// ----------------------------------------------------------------------------

/**
 * @brief Functions and data declarations.
 */

#ifndef ASTC_CODEC_INTERNALS_INCLUDED
#define ASTC_CODEC_INTERNALS_INCLUDED

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include "astc_mathlib.h"

// Temporary workaround to build machine still running VS2013
// Version 1910 is Visual Studio 2017
#if !defined(_MSC_VER) || (_MSC_VER >= 1910)
    #define NORETURN [[noreturn]]
#else
    #define NORETURN
    #define __func__ __FUNCTION__
#endif

// ASTC parameters
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

// error reporting for codec internal errors.
#define ASTC_CODEC_INTERNAL_ERROR() astc_codec_internal_error(__FILE__, __LINE__)
NORETURN void astc_codec_internal_error(const char *filename, int linenumber);

// uncomment this macro to enable checking for inappropriate NaNs;
// works on Linux only, and slows down encoding significantly.
// #define DEBUG_CAPTURE_NAN

// uncomment this macro to enable the ability to log diagnostics using the
// -diag command line switch, which is useful for bug hunting in the encoder
// #define DEBUG_PRINT_DIAGNOSTICS

#ifdef DEBUG_PRINT_DIAGNOSTICS
	extern int print_diagnostics;
#endif

extern int print_tile_errors;
extern int print_statistics;

extern int perform_srgb_transform;
extern int rgb_force_use_of_hdr;
extern int alpha_force_use_of_hdr;

enum astc_decode_mode
{
	DECODE_LDR_SRGB,
	DECODE_LDR,
	DECODE_HDR
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
	int8_t permit_encode;
	int8_t permit_decode;
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
	const decimation_table *decimation_tables[MAX_DECIMATION_MODES + 1];
	block_mode block_modes[MAX_WEIGHT_MODES];

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
	float orig_data[MAX_TEXELS_PER_BLOCK * 4];  // original input data
	float work_data[MAX_TEXELS_PER_BLOCK * 4];  // the data that we will compress, either linear or LNS (0..65535 in both cases)

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

struct error_weighting_params
{
	float rgb_power;
	float rgb_base_weight;
	float rgb_mean_weight;
	float rgb_stdev_weight;
	float alpha_power;
	float alpha_base_weight;
	float alpha_mean_weight;
	float alpha_stdev_weight;
	float rgb_mean_and_stdev_mixing;
	int mean_stdev_radius;
	int enable_rgb_scale_with_alpha;
	int alpha_radius;
	int ra_normal_angular_scale;
	float block_artifact_suppression;
	float rgba_weights[4];

	float block_artifact_suppression_expanded[MAX_TEXELS_PER_BLOCK];

	// parameters that deal with heuristic codec speedups
	int partition_search_limit;
	float block_mode_cutoff;
	float texel_avg_error_limit;
	float partition_1_to_2_limit;
	float lowest_correlation_cutoff;
	int max_refinement_iters;
};

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

struct error_weight_block_orig
{
	float4 error_weights[MAX_TEXELS_PER_BLOCK];
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
	uint8_t unquantized_value_unsc[33];
	/** The scrambling order: value[map[i]] == value_unsc[i] */
	uint8_t scramble_map[32];
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
	FMT_HDR_RGBA = 15,
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
	int component1,
	int component2,
	int component3,
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

// function to find the best partitioning for a given block.

void find_best_partitionings(
	int partition_search_limit,
	const block_size_descriptor* bsd,
	int partition_count,
	const imageblock* pb,
	const error_weight_block* ewb,
	int candidates_to_return,
	int* best_partitions_uncorrelated,
	int* best_partitions_samechroma,
	int* best_partitions_dual_weight_planes);

// use k-means clustering to compute a partition ordering for a block.
void kmeans_compute_partition_ordering(
	const block_size_descriptor* bsd,
	int partition_count,
	const imageblock* blk,
	int *ordering);

// *********************************************************
// functions and data pertaining to images and imageblocks
// *********************************************************

struct astc_codec_image
{
	uint8_t ***data8;
	uint16_t ***data16;
	int xsize;
	int ysize;
	int zsize;
	int padding;

	// Regional average-and-variance information, initialized by
	// compute_averages_and_variances() only if the astc encoder
	// is requested to do error weighting based on averages and variances.
	float4 *input_averages;
	float4 *input_variances;
	float *input_alpha_averages;
};

astc_codec_image *allocate_image(
	int bitness,
	int xsize,
	int ysize,
	int zsize,
	int padding);

void initialize_image(
	astc_codec_image * img);

void destroy_image(
	astc_codec_image * img);

void fill_image_padding_area(
	astc_codec_image * img);

int determine_image_channels(
	const astc_codec_image * img);

// the entries here : 0=red, 1=green, 2=blue, 3=alpha, 4=0.0, 5=1.0
struct swizzlepattern
{
	uint8_t r;
	uint8_t g;
	uint8_t b;
	uint8_t a;
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
 * @param need_srgb_transform   Do we need srgb transform or not?
 * @param swz                   Input data channel swizzle.
 * @param thread_count          The number of threads to use.
 */
void compute_averages_and_variances(
	astc_codec_image * img,
	float rgb_power,
	float alpha_power,
	int avg_var_kernel_radius,
	int alpha_kernel_radius,
	int need_srgb_transform,
	swizzlepattern swz,
	int thread_count);

/*
	Functions to load image from file.
	If successful, return an astc_codec_image object.
	If unsuccessful, returns NULL.

	*result is used to return a result. In case of a successfully loaded image, bits[2:0]
	of *result indicate how many components are present, and bit[7] indicate whether
	the input image was LDR or HDR (0=LDR, 1=HDR).

	In case of failure, *result is given a negative value.
*/

astc_codec_image *load_ktx_uncompressed_image(const char *filename, int padding, int *result);
astc_codec_image *load_dds_uncompressed_image(const char *filename, int padding, int *result);
astc_codec_image *load_tga_image(const char *tga_filename, int padding, int *result);
astc_codec_image *load_image_with_stb(const char *filename, int padding, int *result);

astc_codec_image *astc_codec_load_image(const char *filename, int padding, int *result);

// function to store image to file
// If successful, returns the number of channels in input image
// If unsuccessful, returns a negative number.
int store_ktx_uncompressed_image(const astc_codec_image * img, const char *filename, int bitness);
int store_dds_uncompressed_image(const astc_codec_image * img, const char *filename, int bitness);
int store_tga_image(const astc_codec_image * img, const char *tga_filename, int bitness);

int astc_codec_store_image(const astc_codec_image * img, const char *filename, int bitness, const char **format_string);

int get_output_filename_enforced_bitness(const char *filename);

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
	const astc_codec_image* img1,
	const astc_codec_image* img2,
	int fstop_lo,
	int fstop_hi,
	int show_psnr);

// fetch an image-block from the input file
void fetch_imageblock(
	const astc_codec_image* img,
	imageblock* pb,	// picture-block to initialize with image data
	const block_size_descriptor* bsd,
	// position in picture to fetch block from
	int xpos,
	int ypos,
	int zpos,
	swizzlepattern swz);

// write an image block to the output file buffer.
// the data written are taken from orig_data.
void write_imageblock(
	astc_codec_image* img,
	const imageblock* pb,	// picture-block to initialize with image data
	const block_size_descriptor* bsd,
	// position in picture to write block to.
	int xpos,
	int ypos,
	int zpos,
	swizzlepattern swz);

// helper function to check whether a given picture-block has alpha that is not
// just uniformly 1.
int imageblock_uses_alpha(
	const imageblock * pb);

float compute_imageblock_difference(
	const block_size_descriptor* bsd,
	const imageblock * p1,
	const imageblock * p2,
	const error_weight_block * ewb);

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

float compute_value_of_texel_flt(
	int texel_to_get,
	const decimation_table * it,
	const float *weights);

int compute_value_of_texel_int(
	int texel_to_get,
	const decimation_table* it,
	const int* weights);

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
	astc_decode_mode decode_mode,
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
struct compress_fixed_partition_buffers
{
	endpoints_and_weights* ei1;
	endpoints_and_weights* ei2;
	endpoints_and_weights* eix1;
	endpoints_and_weights* eix2;
	float* decimated_quantized_weights;
	float* decimated_weights;
	float* flt_quantized_decimated_quantized_weights;
	uint8_t* u8_quantized_decimated_quantized_weights;
};

struct compress_symbolic_block_buffers
{
	error_weight_block* ewb;
	error_weight_block_orig* ewbo;
	symbolic_compressed_block* tempblocks;
	imageblock* temp;
	compress_fixed_partition_buffers* plane1;
	compress_fixed_partition_buffers* planes2;
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
	// output data
	int partition_format_specifiers[4][4],
	int quantized_weight[4],
	int quantization_level[4],
	int quantization_level_mod[4]);

void recompute_ideal_colors(
	const block_size_descriptor* bsd,
	int weight_quantization_mode,
	endpoints* ep,	// contains the endpoints we wish to update
	float4* rgbs_vectors,	// used to return RGBS-vectors for endpoint mode #6
	float4* rgbo_vectors,	// used to return RGBS-vectors for endpoint mode #7
	const uint8_t* weight_set8,	// the current set of weight values
	const uint8_t* plane2_weight_set8,	// NULL if plane 2 is not actually used.
	int plane2_color_component,	// color component for 2nd plane of weights; -1 if the 2nd plane of weights is not present
	const partition_info* pi,
	const decimation_table* it,
	const imageblock* pb,	// picture-block containing the actual data.
	const error_weight_block* ewb);

void expand_block_artifact_suppression(
	int xdim,
	int ydim,
	int zdim,
	error_weighting_params* ewp);

// functions pertaining to weight alignment
void prepare_angular_tables(void);

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

float compress_symbolic_block(
	const astc_codec_image* input_image,
	astc_decode_mode decode_mode,
	const block_size_descriptor* bsd,
	const error_weighting_params* ewp,
	const imageblock* blk,
	symbolic_compressed_block* scb,
	compress_symbolic_block_buffers* tmpbuf);

float4 lerp_color_flt(
	const float4 color0,
	const float4 color1,
	float weight,	// 0..1
	float plane2_weight,	// 0..1
	int plane2_color_component);	// 0..3; -1 if only one plane of weights is present.

uint4 lerp_color_int(
	astc_decode_mode decode_mode,
	uint4 color0,
	uint4 color1,
	int weight,	// 0..64
	int plane2_weight,	// 0..64
	int plane2_color_component); // 0..3; -1 if only one plane of weights is present.

void decompress_symbolic_block(
	astc_decode_mode decode_mode,
	const block_size_descriptor* bsd,
	int xpos,
	int ypos,
	int zpos,
	const symbolic_compressed_block* scb,
	imageblock* blk);

physical_compressed_block symbolic_to_physical(
	const block_size_descriptor* bsd,
	const symbolic_compressed_block* sc);

void physical_to_symbolic(
	const block_size_descriptor* bsd,
	physical_compressed_block pb,
	symbolic_compressed_block* res);

uint16_t unorm16_to_sf16(
	uint16_t p);

uint16_t lns_to_sf16(
	uint16_t p);

/* ============================================================================
  Platform-specific functions
============================================================================ */

/**
 * @brief Get the current time.
 *
 * @returns The current time in seconds since arbitrary epoch.
 */
double get_time();

/**
 * @brief Get the number of CPU cores.
 *
 * @returns The number of online or onlineable CPU cores in the system.
 */
int get_cpu_count();

/**
 * @brief Delete (unlink) a file on the filesystem.
 *
|* @param filename The file to delete.
 *
 * @returns Zero on success, non-zero otherwise.
 */
int unlink_file(const char *filename);

/**
 * @brief Launch N worker threads and wait for them to complete.
 *
 * All threads run the same thread function, and have the same thread payload,
 * but are given a unique thread ID (0 .. N-1) as a parameter to the run
 * function to allow thread-specific behavior.
 *
|* @param thread_count The number of threads to spawn.
 * @param func         The function to execute. Must have the signature:
 *                     void (int thread_count, int thread_id, void* payload)
 * @param payload      Pointer to an opaque thread payload object.
 */
void launch_threads(
	int thread_count,
	void (*func)(int, int, void*),
	void *payload);

/**
 * @brief The main entry point.
 *
 * @param argc The number of arguments.
 * @param argb The array of command line arguments.
 *
 * @returns Zero on success, non-zero otherwise.
 */
int astc_main(
	int argc,
	char** argv);

#endif
