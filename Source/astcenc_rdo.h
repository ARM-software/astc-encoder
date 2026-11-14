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
 * @file
 * @brief RDO brick integration, byte search, rate estimation, and distortion interfaces.
 *
 * Rate-distortion optimization (RDO) trades image quality for smaller distribution files. It
 * improves the compression ratio that can be achieved when applying lossless compression, such as
 * ZIP or LZ4, to an ASTC compressed image.
 *
 * RDO is a functionally transparent optimization; the output of the RDO process is still a normal
 * compressed ASTC block that can be used directly by the GPU. The data has just been conditioned
 * so that bit patterns are more similar to surrounding blocks, making it more compressible.
 *
 * Entropy reduction transform
 * ===========================
 *
 * The algorithm used by ASTCEnc is based on Rich Geldreich's public-domain Entropy Reduction
 * Transform (ERT) implementation in bc7enc_rdo.
 *
 *     See https://github.com/richgel999/bc7enc_rdo/blob/master/ert.cpp.
 *
 * The original ERT implementation is a generic implementation that is oblivious to the data format
 * that is it transforming. This makes it trivial to integrate with any data source, as it only
 * needs some cost function to use, but it makes applying the ERT very slow. Compressing the Kodak
 * suite using ASTC -thorough compression takes around 1.6 seconds. Applying the out-of-the-box ERT
 * transform with a 1KB history buffer takes 60 seconds at lambda=0.2. This is really too slow to
 * apply at scale on a large project.
 *
 * This implementation is a heavily optimized version of the original ERT transform. It includes
 * both generic optimizations, and format-specific optimizations that exploit the fact that we know
 * it is transforming compressed ASTC blocks. On our benchmarking machine the specialized ASTC ERT
 * (AERT) implementation takes just 3.5 seconds for similar PSNR and ZIP-compressed file size.
 *
 * Algorithm
 * =========
 *
 * When using RDO the main ASTC compression will generate a "brick" containing 2KB of compressed
 * data (128 compressed ASTC blocks). Bricks are either arranged linearly, or as a 16x8 Morton
 * order block. Each encoder brick is AERT optimized independently.
 *
 * The AERT process is:
 *
 * 1. Gather valid physical blocks into a compact stream in brick traversal order.
 * 2. Identify protected bits that cannot be mutated by the ERT. We allow only color endpoints and
 *    weights to change, protecting everything else.
 * 2. Learn DEFLATE literal, length, and distance Huffman prices from a greedy parse of the brick.
 *    A sliding history follows the accepted prefix; token prices are refreshed at the midpoint.
 * 3. Visit blocks in stream order. Try copying byte ranges from the same byte offsets in nearby
 *    predecessors, then optionally add a second non-overlapping range to the best first match.
 *    Structural masks, candidate hashes, and search budgets prune trials before decoding.
 * 4. Minimize weighted MSE * smoothness_scale + internal_lambda * estimated_bits. MSE uses
 *    normalized channel weights in the 0..255 domain; smooth blocks receive extra protection.
 *    Decode only affected payloads and stop scoring once the trial cannot beat the best cost.
 *    Internal lambda is the user lambda multiplied by rdo_lambda_scale; zero disables RDO.
 * 5. Keep each winning block's exact physical bytes as history for subsequent blocks, then
 *    scatter the optimized brick back to the output stream. .
 *
 * The rate model estimates storage cost without running a lossless compressor for every trial.
 * It caches optimal parses of unchanged byte intervals and adds the proposed match costs.
 * Matches can continue across ASTC block boundaries, but history resets at brick boundaries.
 * Estimated improvements need not reduce the final DEFLATE file size by the same amount.
 *
 * Limitaitons
 * ===========
 *
 * In the current implementation, only LDR textures are supported.
 */

#ifndef ASTCENC_RDO_INCLUDED
#define ASTCENC_RDO_INCLUDED

#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>

#include "astcenc_morton.h"

struct astcenc_config;
struct astcenc_contexti;
struct astcenc_image;
struct astcenc_swizzle;
struct block_size_descriptor;
struct image_block;
struct rdo_brick_context;
struct symbolic_compressed_block;

/* ============================================================================
  Compile-time tuning

  Lambda and brick traversal order are runtime controls in astcenc_config.
  Brick layouts are defined by brick_layout in astcenc_morton.h.
============================================================================ */

/** @brief Scale user lambda to roughly match the original ERT quality and size tradeoff. */
static constexpr float rdo_lambda_scale { 3.0f };

/** @brief Sliding rate-model history in bytes; a positive multiple of 16. */
static constexpr unsigned int rdo_history_bytes { 1024 };

/** @brief Minimum predecessor scan budget for large ASTC footprints. */
static constexpr unsigned int rdo_min_lookback { 32 };

/** @brief Maximum predecessor scan budget for small ASTC footprints. */
static constexpr unsigned int rdo_max_lookback { 64 };

/** @brief Base scan budget before the inverse block-dimension adjustment; an even count. */
static constexpr unsigned int rdo_search_base_lookback { 22 };

/** @brief Inverse-dimension scale in predecessor pairs; larger values search more history. */
static constexpr unsigned int rdo_search_dimension_scale { 30 };

/** @brief Minimum estimated saving over the first match required for a second match, in bits. */
static constexpr float rdo_min_second_match_saving { 16.0f };

/** @brief Maximum physical bit changes allowed in a nonconstant trial, in [0, 128]. */
static constexpr unsigned int rdo_max_changed_bits { 32 };

/** @brief Trial RMS ceiling relative to the original, with an original MSE floor of one. */
static constexpr float rdo_max_rms_increase { 10.0f };

/** @brief Source standard deviation below which smooth-block protection applies. */
static constexpr float rdo_max_smooth_std_dev { 18.0f };

/** @brief Maximum distortion multiplier for smooth blocks. */
static constexpr float rdo_smooth_mse_scale { 10.0f };

/** @brief True to search for a second non-overlapping copied byte range. */
static constexpr bool rdo_is_two_matches_enabled { true };

/** @brief True to refresh Huffman prices after optimizing half the brick. */
static constexpr bool rdo_is_midpoint_refresh_enabled { true };

/** @brief Candidate deduplication entries; a power of two trading cache size against collisions. */
static constexpr unsigned int rdo_trial_hash_entries { 8192 };

/** @brief Rate-dictionary index bits, in [1, 16]; more bits shorten collision chains. */
static constexpr unsigned int rdo_rate_hash_bits { 12 };

static_assert(rdo_max_smooth_std_dev > 0.0f);
static_assert(rdo_smooth_mse_scale > 0.0f);
static_assert((rdo_history_bytes > 0) && ((rdo_history_bytes % 16) == 0));
static_assert((rdo_min_lookback > 0) && (rdo_min_lookback <= rdo_max_lookback));
static_assert(((rdo_min_lookback | rdo_max_lookback) & 1u) == 0);
static_assert((rdo_search_base_lookback & 1u) == 0);
static_assert(rdo_max_changed_bits <= 128);
static_assert(std::has_single_bit(rdo_trial_hash_entries));
static_assert((rdo_rate_hash_bits >= 1) && (rdo_rate_hash_bits <= 16));

/* ============================================================================
  Type aliases
============================================================================ */

/**
 * @brief Decode a bounded integer sequence using a fixed quantization alphabet.
 *
 * @param      character_count   The number of encoded values.
 * @param[in]  input_data        The packed input sequence.
 * @param[out] output_data       The decoded values.
 * @param      bit_offset        The initial bit offset in the input sequence.
 */
using rdo_ise_decode_func = void(
	unsigned int character_count,
	const uint8_t* input_data,
	uint8_t* output_data,
	unsigned int bit_offset);

/* ============================================================================
  Structures
============================================================================ */

/** @brief Parameters for the entropy reduction transform. */
struct rdo_parameters
{
	/** @brief Rate weight in the weighted MSE plus estimated bits objective. */
	float lambda { 0.0f };

	/** @brief Maximum number of predecessors searched for speculative copied ranges. */
	unsigned int search_lookback { rdo_max_lookback };
};

/** @brief Original block distortion and source statistics used to initialize the search. */
struct rdo_block_statistics
{
	/** @brief The original block's weighted MSE in the 0..255 domain. */
	float mse;

	/** @brief The source standard deviation used for smooth-block protection. */
	float std_dev;
};

/** @brief A byte match selected for the current block. */
struct rdo_byte_match
{
	/** @brief The match length in bytes; zero indicates no match. */
	unsigned int length { 0 };

	/** @brief The absolute source byte offset in the input stream. */
	unsigned int source { 0 };

	/** @brief The destination byte offset within the current block. */
	unsigned int offset { 0 };

	/** @brief The estimated cost of the match in bits. */
	float bits { 0.0f };
};

/** @brief Fixed ISE layout for trials retaining the original nonconstant block structure. */
struct rdo_block_decode_state
{
	/** @brief Fixed-quantization ISE decoder for weight indices. */
	rdo_ise_decode_func* weight_decode;

	/** @brief Fixed-quantization ISE decoder for endpoint indices. */
	rdo_ise_decode_func* color_decode;

	/** @brief The number of weights in each plane. */
	unsigned int weight_count;

	/** @brief The total number of encoded weights. */
	unsigned int real_weight_count;

	/** @brief The weight quantization alphabet index. */
	unsigned int weight_quant_mode;

	/** @brief The total number of endpoint integers. */
	unsigned int color_count;

	/** @brief The endpoint sequence's bit offset. */
	unsigned int color_offset;

	/** @brief The number of endpoint integers in each partition. */
	unsigned int color_counts[4];

	/** @brief Bytes containing weight ISE bits. */
	unsigned int weight_byte_mask;

	/** @brief Bytes containing endpoint ISE bits. */
	unsigned int color_byte_mask;

	/** @brief True if the block uses two weight planes. */
	bool is_dual_plane;
};

/**
 * @brief Original reconstructed values reusable when a trial leaves a payload unchanged.
 *
 * The caller updates the unchanged flags for each trial using conservative payload byte masks.
 */
struct rdo_block_difference_cache
{
	/** @brief Unpacked endpoint pairs for each partition. */
	int endpoints[4][2][4];

	/** @brief Interpolated weight planes, padded for the largest supported SIMD width. */
	uint8_t weights[2][224];

	/** @brief Original weighted squared errors in raster texel order. */
	float texel_errors[224];

	/** @brief Maximum contributing weight count for each fixed-layout SIMD batch. */
	uint8_t weight_counts[224 / 4];

	/** @brief True if this trial preserves the original endpoint payload. */
	bool is_colors_unchanged;

	/** @brief True if this trial preserves the original weight payload. */
	bool is_weights_unchanged;
};

/* ============================================================================
  Classes
============================================================================ */

/** @brief Brick-local Huffman prices and cached rates for unchanged byte intervals. */
class rdo_rate_model
{
public:
	/**
	 * @brief Learn token prices from a greedy parse of the original brick.
	 *
	 * Prices remain fixed until the next init() call. The search retrains them halfway through
	 * the brick; the match dictionary follows accepted changes between retraining points.
	 *
	 * @param[in] data          The original compact physical block stream.
	 * @param     block_count   The number of valid blocks.
	 * @param     lookback      The maximum history length in blocks.
	 */
	void init(
		const uint8_t* data,
		unsigned int block_count,
		unsigned int lookback);

	/**
	 * @brief Prepare rates using the current optimized prefix.
	 *
	 * Blocks must be prepared in increasing index order between init() calls. The stream prefix
	 * before the current block must contain all accepted changes and remain unchanged afterward.
	 *
	 * @param[in] data            The compact physical block stream.
	 * @param     block_index     The block to optimize.
	 * @param     history_start   The first retained history block.
	 */
	void prepare(
		const uint8_t* data,
		unsigned int block_index,
		unsigned int history_start);

	/**
	 * @brief Estimate a match's token cost.
	 *
	 * @param distance   The backward distance in bytes.
	 * @param length     The match length in bytes, in the range 3..16.
	 *
	 * @return The match cost in bits.
	 */
	float get_match_cost(
		unsigned int distance,
		unsigned int length) const;

	/**
	 * @brief Estimate a block using zero, one, or two copied ranges.
	 *
	 * Bytes outside the copied ranges must match the original block passed to prepare().
	 * The two ranges must not overlap. A zero-length first match requests the original cost.
	 * Estimates at or above that original cost cannot improve the search and must be rejected.
	 *
	 * @param[in] first    The first match, or an empty match for the original block.
	 * @param[in] second   The second match, or an empty match for a single-match trial.
	 *
	 * @return The estimated block cost in bits.
	 */
	float estimate(
		const rdo_byte_match& first,
		const rdo_byte_match& second) const;

	/**
	 * @brief Get the frozen Huffman price of a literal byte.
	 *
	 * @param value   The literal byte.
	 *
	 * @return The literal cost in bits.
	 */
	float get_literal_cost(
		uint8_t value) const;

private:
	static_assert(brick_layout::block_count <= 32768 / 16,
	              "RDO brick distances must fit in the DEFLATE window");

	/** @brief The number of input bytes in a full encoder brick. */
	static constexpr unsigned int max_bytes = brick_layout::block_count * 16;

	/** @brief The number of buckets in the three-byte match dictionary. */
	static constexpr unsigned int hash_size = 1u << rdo_rate_hash_bits;

	/** @brief Frozen code lengths for literal, length, and end symbols. */
	uint8_t literal_length_bits[286];

	/** @brief Length-symbol and extra-bit costs, indexed by match length. */
	uint8_t length_costs[17];

	/** @brief Distance-symbol and extra-bit costs, indexed by byte distance. */
	uint8_t distance_costs[max_bytes + 1];

	/** @brief Newest completed source substring in each hash bucket. */
	int heads[hash_size];

	/** @brief Earlier source substring in the same hash bucket. */
	int previous[max_bytes];

	/** @brief The first source byte not yet entered into the match dictionary. */
	unsigned int indexed_bytes;

	/** @brief Maximum match distance in the sliding history, in bytes. */
	unsigned int history_bytes;

	/** @brief The current block's unchanged physical bytes. */
	uint8_t original[16];

	/** @brief Minimum token cost for each unchanged half-open byte interval. */
	float intervals[17][17];
};

/* ============================================================================
  Brick integration
============================================================================ */

/**
 * @brief Optimize a brick in its selected traversal order, then scatter accepted physical bytes.
 *
 * Disabled RDO and unsupported configurations leave the stream unchanged.
 *
 * @param[in]     ctx          The compressor context.
 * @param[in]     image        The source image.
 * @param[in]     swizzle      The compression swizzle.
 * @param[in]     valid_mask   The bitmask of valid brick positions in the brick.
 * @param[in]     offsets      The raster byte offset at each valid brick position.
 * @param[in]     blocks       The original symbolic blocks at each valid brick position.
 * @param[in,out] buffer       The raster-ordered physical output stream.
 */
void rdo_optimize_physical_brick(
	const astcenc_contexti& ctx,
	const astcenc_image& image,
	const astcenc_swizzle& swizzle,
	const brick_layout::validity_mask& valid_mask,
	const size_t* offsets,
	const symbolic_compressed_block* blocks,
	uint8_t* buffer);

/**
 * @brief Load an original block and initialize its distortion caches and source statistics.
 *
 * @param[in,out] ctx           The brick distortion context.
 * @param         block_index   The compact block index in brick traversal order.
 *
 * @return The original weighted MSE and source standard deviation.
 */
rdo_block_statistics rdo_initialize_block(
	rdo_brick_context& ctx,
	unsigned int block_index);

/**
 * @brief Evaluate a physical trial using decoded weighted MSE.
 *
 * Initialize the current block with rdo_initialize_block() before scoring trials. byte_mask may
 * include unchanged bytes.
 *
 * @param[in,out] ctx         The initialized brick distortion context.
 * @param[in]     physical    The trial physical block.
 * @param         byte_mask   A conservative mask of changed physical bytes.
 * @param         mse_limit   The MSE bound, or infinity for a complete score.
 *
 * @return Weighted MSE in the 0..255 domain, or @c std::numeric_limits<float>::max()
 *         for losing trials.
 */
float rdo_block_difference(
	rdo_brick_context& ctx,
	const uint8_t* physical,
	unsigned int byte_mask,
	float mse_limit);

/* ============================================================================
  Search and rate estimation
============================================================================ */

/**
 * @brief Prepare polynomial prefixes for all substrings of one physical block.
 *
 * @param[in]  data       The source block bytes.
 * @param[out] prefixes   The seventeen polynomial prefixes, including the empty prefix.
 */
void rdo_prepare_substring_hashes(
	const uint8_t data[16],
	uint32_t prefixes[17]);

/**
 * @brief Compute salted hashes for every fixed-offset substring of one length.
 *
 * @param[in]  prefixes   The prefixes, padded to 16 + ASTCENC_SIMD_WIDTH initialized entries.
 * @param      length     The substring length, in the range 3..16.
 * @param[out] hashes     The resulting 17 - length hashes; remaining entries are untouched.
 */
void rdo_hash_substrings(
	const uint32_t* prefixes,
	unsigned int length,
	uint32_t hashes[16]);

/**
 * @brief Improve estimated compressibility by copying byte ranges from preceding blocks.
 *
 * Input blocks must be valid encoder output. History and match continuation state are local to
 * this call. Constant blocks remain unchanged. Only trials preserving the protected layout bits
 * reach distortion scoring.
 *
 * @param[in,out] ctx           The brick distortion context.
 * @param[in,out] data          The contiguous stream of 16-byte compressed blocks.
 * @param         block_count   The number of blocks, at most brick_layout::block_count.
 * @param[in]     masks         The protected layout bits for each 16-byte block.
 * @param[in]     params        The distortion and history settings.
 */
void rdo_reduce_entropy(
	rdo_brick_context& ctx,
	uint8_t* data,
	unsigned int block_count,
	const uint8_t* masks,
	const rdo_parameters& params);

/* ============================================================================
  Physical layout and payload decoding
============================================================================ */

/**
 * @brief Select a dedicated decoder for a fixed quantization alphabet.
 *
 * @param quant_level   The quantization alphabet index.
 *
 * @return The fixed-quantization ISE decoder.
 */
rdo_ise_decode_func* rdo_get_ise_decoder(unsigned int quant_level);

/**
 * @brief Decode pure trit or quint weight sequences without an intermediate index array.
 *
 * @param      quant_level   QUANT_3 or QUANT_5.
 * @param      count         The number of stored weights.
 * @param[in]  input_data    The reversed physical block, with the sequence at bit zero.
 * @param[out] output_data   The unquantized weights, with exactly count bytes written.
 */
void rdo_decode_ise_radix_weights(
	unsigned int quant_level,
	unsigned int count,
	const uint8_t* input_data,
	uint8_t* output_data);

/**
 * @brief Prepare payload decoding from the encoder's original symbolic layout.
 *
 * @param[in] bsd     The block size descriptor.
 * @param[in] block   The original valid nonconstant symbolic block.
 *
 * @return The fixed payload layout.
 */
rdo_block_decode_state rdo_prepare_block_decode(
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& block);

/**
 * @brief Decode changed payload bytes, leaving untouched payload storage unspecified.
 *
 * Layout bits must match the original block. The byte mask may overestimate changes, but must
 * include every changed byte to avoid reusing a stale payload.
 *
 * @param[in]     state       The fixed nonconstant layout.
 * @param[in]     physical    The candidate physical bytes.
 * @param         byte_mask   A conservative mask of changed physical bytes.
 * @param[in,out] block       The trial symbolic block with its fixed layout initialized.
 */
void rdo_decode_changed_payload(
	const rdo_block_decode_state& state,
	const uint8_t physical[16],
	unsigned int byte_mask,
	symbolic_compressed_block& block);

/**
 * @brief Identify bits selecting the block layout rather than endpoint or weight values.
 *
 * @param[in]  bsd     The block size descriptor.
 * @param[in]  block   The original valid symbolic encoding.
 * @param[out] mask    The protected layout bits, or all bits for constant blocks.
 */
void rdo_get_structure_mask(
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& block,
	uint8_t mask[16]);

/* ============================================================================
  Distortion scoring
============================================================================ */

/**
 * @brief Bound the MSE of a candidate that can improve the current rate-distortion cost.
 *
 * @param best_cost   The current best cost.
 * @param rate_cost   The candidate's weighted rate cost.
 * @param mse_recip   The reciprocal of the positive distortion multiplier.
 *
 * @return The MSE limit, or zero if the rate alone cannot improve the cost.
 */
float rdo_get_mse_limit(
	float best_cost,
	float rate_cost,
	float mse_recip);

/**
 * @brief Scale cached source texels into the distortion metric's 0..255 domain.
 *
 * @param[in,out] original      The loaded source block, initially in the 0..65535 domain.
 * @param         texel_count   The number of source texels.
 */
void rdo_prepare_difference_source(
	image_block& original,
	unsigned int texel_count);

/**
 * @brief Cache reconstructed endpoints and interpolated weights for an original block.
 *
 * @param[in]  config   The LDR compressor configuration.
 * @param[in]  bsd      The block size descriptor.
 * @param[in]  block    The valid nonconstant symbolic encoding from the core compressor.
 * @param[out] cache    The reconstructed values to reuse.
 */
void rdo_prepare_difference_cache(
	const astcenc_config& config,
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& block,
	rdo_block_difference_cache& cache);

/**
 * @brief Score the original reconstruction and cache its errors without decoding it again.
 *
 * @param[in]     bsd        The block size descriptor.
 * @param[in]     block      The original valid nonconstant symbolic block.
 * @param[in]     original   The prepared source texels and normalized channel weights.
 * @param[in,out] cache      The prepared endpoints and weights, receiving original errors.
 *
 * @return The original weighted MSE and source standard deviation.
 */
rdo_block_statistics rdo_compute_original_difference(
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& block,
	const image_block& original,
	rdo_block_difference_cache& cache);

/**
 * @brief Decode and score a full LDR block, stopping when a trial cannot win.
 *
 * The block must be nonconstant. The source uses normalized channel weights.
 *
 * @param[in] config      The LDR compressor configuration.
 * @param[in] bsd         The block size descriptor.
 * @param[in] block       The symbolic encoding.
 * @param[in] original    The prepared 0..255 source texels and normalized channel weights.
 * @param     mse_limit   The MSE bound, or infinity for a complete score.
 * @param[in] cache       The initialized original reconstruction cache, including per-texel errors.
 *
 * @return Weighted MSE in the 0..255 domain, or @c std::numeric_limits<float>::max()
 *         for losing trials.
 */
float rdo_compute_block_difference(
	const astcenc_config& config,
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& block,
	const image_block& original,
	float mse_limit,
	const rdo_block_difference_cache& cache);

#endif
