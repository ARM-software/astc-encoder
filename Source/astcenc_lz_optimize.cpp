/**
 * @brief Implementation of LZ backend rate-distortion optimization.
 *
 * This file implements a rate-distortion optimization algorithm for ASTC compressed textures.
 * The algorithm aims to improve compression by exploiting redundancy between neighboring blocks
 * while maintaining visual quality.
 *
 * We keep a Move-To-Front (MTF) list of recently used endpoints (including header information)
 * and weight bits and try to build new blocks using them. Whenever we do, this allows a backend
 * LZ coder (like Deflate, Zstd, LZ4, LZMA, ...) to produce a match, decreasing the compressed size
 * significantly. The trade-off is that the blocks we get this way are worse in terms of error than
 * the blocks we get if we optimize the parameters for the block contents.
 *
 * We use a combined rate-distortion score to weight this trade-off, controlled by a single
 * parameter (lambda). Low values of lambda heavily prioritize minimizing error (maximum quality).
 * Higher values of lambda care more about the bit rate estimate and will generally lead to higher
 * errors.
 *
 * It is important for visual quality to not just use raw squared error or equivalently PSNR for
 * this. Since we are intentionally increasing error, we should only do so where errors are less
 * likely to be noticed.
 *
 * To this end, we determine per-texel weights for the image to determine where errors are more or
 * less likely to be noticed. In general, the errors that result from the types of block
 * modifications we're evaluating are most apparent in smooth regions, less so in noisy areas or
 * near edges.
 *
 * Therefore, we identify the amount of high-frequency energy near any given texel; texels in smooth
 * regions try to keep the error low, while texels in high-frequency areas use looser tolerances.
 */

#if !defined(ASTCENC_DECOMPRESS_ONLY)

#include "astcenc_internal_entry.h"
#include <cstring>

/** @brief Size of an ASTC block in bytes. */
static constexpr unsigned int ASTC_BLOCK_BYTES{16};
/** @brief Number of best candidates to keep for each block during optimization. */
static constexpr unsigned int BEST_CANDIDATES_COUNT{8};
/** @brief Maximum number of blocks to process per worker thread item, should be a power of 2. */
static constexpr unsigned int MAX_BLOCKS_PER_ITEM{4096};
/** @brief Cache size for block decompression results, must be a power of 2. */
static constexpr unsigned int CACHE_SIZE{MAX_BLOCKS_PER_ITEM};

// Check for known little-endian platforms. Windows is always LE.
#if defined(_WIN32) || defined(__LITTLE_ENDIAN__) || (defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)

/**
 * @brief Load a 64-bit little-endian value.
 */
static inline uint64_t load64le(const uint8_t* p)
{
	uint64_t result;
	std::memcpy(&result, p, sizeof(result));
	return result;
}

/**
 * @brief Store a 64-bit little-endian value.
 */
static inline void store64le(uint64_t value, uint8_t* p)
{
	std::memcpy(p, &value, sizeof(value));
}

#else

/**
 * @brief Load a 64-bit little-endian value.
 */
static inline uint64_t load64le(const uint8_t* p)
{
	uint64_t result = 0;
	for (int i = 0; i < 8; i++)
	{
		result |= static_cast<uint64_t>(p[i]) << (i * 8);
	}

	return result;
}

/**
 * @brief Store a 64-bit little-endian value.
 */
static inline void store64le(uint64_t value, uint8_t* p)
{
	for (int i = 0; i < 8; i++)
	{
		p[i] = static_cast<uint8_t>(value >> (i * 8));
	}
}

#endif

/**
 * @brief Representation of ASTC blocks. Includes facilities for some basic bitwise operations.
 */
struct physical_block
{
	union
	{
		/** @brief The block bytes. */
		uint8_t bytes[16];
		/** @brief The bytes viewed as uint64s. Endian-dependent, internal use only. */
		uint64_t uint64[2];
	};

	/**
	 * @brief Construct an all-zero physical block.
	 */
	physical_block() : uint64 { 0, 0 }
	{
	}

	/**
	 * @brief Construct from raw byte data.
	 */
	explicit physical_block(const uint8_t* data)
	{
		std::memcpy(bytes, data, 16);
	}

	/**
	 * @brief Creates a bit mask for the specified "count" top-most bits.
	 */
	static physical_block top_bits_mask(unsigned int count)
	{
		if (count == 0)
		{
			return physical_block { 0, 0 };
		}
		else if (count <= 64)
		{
			// count is in [1,64], so 64 - count is in [0,63]
			physical_block block;
			store64le(~uint64_t(0) << (64 - count), block.bytes + 8);
			return block;
		}
		else
		{
			if (count > 128)
				count = 128;

			// count is in [65,128], so 128 - count is in [0,63]
			physical_block block { 0, ~uint64_t(0) };
			store64le(~uint64_t(0) << (128 - count), block.bytes);
			return block;
		}
	}

	// NOTE: bitwise operations as well as equality/inequality comparison are naturally
	// endian-agnostic.

	/**
	 * @brief Bitwise AND
	 */
	physical_block operator &(const physical_block& other) const
	{
		return physical_block { uint64[0] & other.uint64[0], uint64[1] & other.uint64[1] };
	}

	/**
	 * @brief Bitwise OR
	 */
	physical_block operator |(const physical_block& other) const
	{
		return physical_block { uint64[0] | other.uint64[0], uint64[1] | other.uint64[1] };
	}

	/**
	 * @brief Bitwise NOT
	 */
	physical_block operator ~() const
	{
		return physical_block { ~uint64[0], ~uint64[1] };
	}

	/**
	 * @brief Equality comparison
	 */
	bool operator ==(const physical_block& other) const
	{
		return uint64[0] == other.uint64[0] && uint64[1] == other.uint64[1];
	}

	/**
	 * @brief Inequality comparison
	 */
	bool operator !=(const physical_block& other) const
	{
		return !(*this == other);
	}

	/**
	 * @brief Compute hash for block
	 */
	uint32_t hash() const
	{
		// FNV-1a inspired constants for 64-bit operations
		const uint64_t PRIME64_1 = 0x9E3779B185EBCA87ULL;
		const uint64_t PRIME64_2 = 0xC2B2AE3D27D4EB4FULL;

		// Get the two 64-bit words directly.
		// We want the little-endian interpretation.
		uint64_t h1 = load64le(bytes + 0);
		uint64_t h2 = load64le(bytes + 8);

		// Mix the first word
		h1 *= PRIME64_1;
		h1 = (h1 << 31) | (h1 >> 33);
		h1 *= PRIME64_2;

		// Mix the second word
		h2 *= PRIME64_2;
		h2 = (h2 << 29) | (h2 >> 35);
		h2 *= PRIME64_1;

		// Combine the results
		uint32_t result = (uint32_t)(h1 ^ h2);

		// Final avalanche
		result ^= result >> 15;
		result *= 0x85ebca6b;
		result ^= result >> 13;
		result *= 0xc2b2ae35;
		result ^= result >> 16;

		return result;
	}

private:
	/**
	 * @brief Constructor from 64-bit halves.
	 *
	 * This is endian-dependent, so for internal use only.
	 */
	physical_block(uint64_t lo, uint64_t hi)
		: uint64 { lo, hi }
	{
	}
};

/** @brief Division that rounds up. */
static size_t ceil_div(size_t a, size_t b)
{
	return (a + b - 1) / b;
}

/** @brief Decode a raster-scan order index into 3D coordinates given X/Y dimension of an array. */
template<typename T>
static void decode_index_to_coords(
	T index,
	T xdim,
	T ydim,
	T& x,
	T& y,
	T& z
) {
	T quot = index / xdim;
	x = index % xdim;
	z = quot / ydim;
	y = quot % ydim;
}

/** @brief Approximate log2, same approximation as astcenc_vecmathlib.h log2 (but scalar). */
static inline float log2_approx(float v)
{
	union { float f; int i; } u;

	u.f = v;
	float e = static_cast<float>((u.i >> 23) - 127);
	u.i = (u.i & 0x007FFFFF) | 0x3F800000;
	float m = u.f;

	// Polynomial fit of log2(x)/(x - 1), for x in range [1, 2)
	float p = POLY4(m,
	                2.8882704548164776201f,
	               -2.52074962577807006663f,
	                1.48116647521213171641f,
	               -0.465725644288844778798f,
	                0.0596515482674574969533f);

	// Increases the polynomial degree, but ensures that log2(1) == 0
	p = p * (m - 1.0f);

	return p + e;
}

/**
 * @brief Histograms for the first few bytes of blocks that select the mode.
 */
struct mode_byte_histogram
{
	/** @brief Number of initial bytes to histogram. */
	static constexpr unsigned int NHISTO = 2;

	/** @brief Everyt his many updates, decay existing frequncies for adaptation. */
	static constexpr unsigned int DECAY_PERIOD = 256;

	/** @brief Implicit count in every bucket so nothing has 0 frequency. */
	static constexpr float BIAS = 0.25f;

	/** @brief Histograms for the byte values. Due to periodic decay, uint16_t is sufficient. */
	uint16_t counts[NHISTO][256];

	/** @brief Sum of all values in the given histogram slice. */
	uint16_t totals[NHISTO];

	/** @brief Counter to keep track of when to decay. */
	unsigned int counter;

	/** @brief Reset the histogram. */
	void reset()
	{
		std::memset(counts, 0, sizeof(counts));
		std::memset(totals, 0, sizeof(totals));
		counter = 0;
	}

	/** @brief Update the histogram with the stats for a given block. */
	void update(const physical_block& value)
	{
		// Periodically decay the histogram
		counter++;
		if (counter == DECAY_PERIOD)
		{
			counter = 0;

			for (unsigned int j = 0; j < NHISTO; j++)
			{
				uint16_t new_total = 0;
				for (unsigned int i = 0; i < 256; i++)
				{
					counts[j][i] >>= 1;
					new_total += counts[j][i];
				}

				totals[j] = new_total;
			}
		}

		for (unsigned int j = 0; j < NHISTO; j++)
		{
			counts[j][value.bytes[j]]++;
			totals[j]++;
		}

	}

	/** @brief Estimates the cost of the histogrammed bytes for a given block. */
	float cost(const physical_block& value) const
	{
		float cost_numer = 1.0f;
		float cost_denom = 1.0f;
		for (unsigned int i = 0; i < NHISTO; i++)
		{
			cost_numer *= totals[i] + 256 * BIAS; // implicit extra count of BIAS in every bucket
			cost_denom *= counts[i][value.bytes[i]] + BIAS;
		}

		return log2_approx(cost_numer / cost_denom);
	}
};

/**
 * @brief Move-to-front (MTF) list structure keeping track of recently used values.
 *
 * We keep track of recent endpoint and weight values to try and reuse them for other blocks.
 */
struct mtf_list
{
	/** @brief The list of blocks. Most recent is in front. */
	physical_block* list;

	/** @brief Current number of items in the list. Limited by capacity. */
	unsigned int size;

	/** @brief Maximum number of items in the list. */
	unsigned int capacity;

	explicit mtf_list(unsigned int list_capacity)
		: list(nullptr), size(0), capacity(list_capacity)
	{
		list = new physical_block[capacity];
	}

	~mtf_list()
	{
		delete[] list;
	}

	/** @brief Reset the list to empty state. */
	void reset()
	{
		size = 0;
	}

	/** @brief Find a given value in the list (masked by "mask"), return index or -1 if not found. */
	int search(const physical_block& value, const physical_block& mask) const
	{
		physical_block masked_value = value & mask;
		for (unsigned int i = 0; i < size; i++)
		{
			if ((list[i] & mask) == masked_value)
			{
				return i;
			}
		}

		return -1; // Not found
	}

	/** @brief Update MTF list after using a given value. pos is the value returned by search. */
	int update(const physical_block& value, int pos)
	{
		if (pos == -1)
		{
			// If not found, insert at the end, growing up to capacity
			if (size < capacity)
			{
				size++;
			}
			pos = size - 1;
		}

		// Move the found value to the front of the list
		for (int i = pos; i > 0; i--)
		{
			list[i] = list[i - 1];
		}
		list[0] = value;

		return pos;
	}

	/** @brief Performs a combined "search" and "update". */
	int encode(const physical_block& value, const physical_block& mask)
	{
		return update(value, search(value, mask));
	}
};

/**
 * @brief Computes approximate cost in bits for encoding a block.
 *
 * @param mtf_pos_endpoints  MTF list position for endpoints/config portion of block.
 * @param mtf_pos_weights    MTF list position for weights portion of block.
 * @param block              The block itself.
 * @param weight_bits        Size of weights portion of the block in bits.
 * @param mode_histo         The mode byte histogram used to estimate cost of the first few bytes.
 */
static float calculate_bit_cost_simple(
	int mtf_pos_endpoints,
	int mtf_pos_weights,
	const physical_block& block,
	int weight_bits,
	const mode_byte_histogram& mode_histo
) {
	float total_cost;
	int endpoint_bits = 128 - weight_bits;

	if (mtf_pos_endpoints == -1)
	{
		// Assume bytes 2 and onwards are random, but score first 2 using histo
		total_cost = endpoint_bits - 16.0f + mode_histo.cost(block);
	}
	else
	{
		// Heuristic match cost, don't over-tune this
		total_cost = 10.0f + log2_approx(mtf_pos_endpoints + 32.0f);

		// If both matches are from the same block, it's just one big match
		if (mtf_pos_weights == mtf_pos_endpoints)
		{
			return total_cost;
		}

	}

	// If endpoints ended in the middle of a byte and we don't have a full match
	// treat the remaining bits in that byte as literals
	if (endpoint_bits & 7)
	{
		total_cost += 8.0f - (endpoint_bits & 7);
	}

	if (mtf_pos_weights == -1)
	{
		// Bill the remaining weight bits (the first byte, if partial, is already accounted for)
		total_cost += static_cast<float>(weight_bits & ~7);
	}
	else
	{
		// Heuristic match cost, don't over-tune this
		total_cost += 10.0f + log2_approx(mtf_pos_weights + 32.0f);
	}

	return total_cost;
}

/**
 * @brief Calculate Sum of Squared Differenced (SSD) with per-texel weights.
 *
 * @param img1             Interleaved RGBA data for the first block.
 * @param img2             Interleaved RGBA data for the second block.
 * @param texel_count      Number of texels in the block.
 * @param weights          The per-texel weights (from activity calculation).
 * @param channel_weights  Per-channel importance weights.
 */
static inline float calculate_ssd_weighted(
	const float* img1,
	const float* img2,
	unsigned int texel_count,
	const float* weights,
	const vfloat4& channel_weights
) {
	vfloat4 sum = vfloat4::zero();
	for (unsigned int i = 0; i < texel_count; i++)
	{
		vfloat4 diff = vfloat4(img1 + i*4) - vfloat4(img2 + i*4);
		haccumulate(sum, diff * diff * vfloat4::load1(weights + i));
	}
	return dot_s(sum, channel_weights);
}

/** @brief Extracts the block mode ID from an encoded block. */
static inline unsigned int get_block_mode(const uint8_t* block)
{
	return (block[0] | (block[1] << 8)) & 0x7ff;
}

/**
 * @brief Determine the weight bits count for an encoded block.
 *
 * @param block   The encoded block.
 * @param bsd     The block size descriptor.
 */
static inline int get_weight_bits(
	const uint8_t* block,
	const block_size_descriptor* bsd
) {
	return bsd->weight_bits_for_mode[get_block_mode(block)];
}

/** @brief Simple Xorshift32 RNG used for block seeding. */
static inline uint32_t xorshift32(uint32_t& state)
{
	state ^= state << 13;
	state ^= state >> 17;
	state ^= state << 5;
	return state;
}

/**
 * @brief Applies UNORM8-style rounding to 16-bit integer values.
 *
 * @param colori The input color value (ints in [0,65535])
 * @param u8_mask The mask of channels that should have 8-bit rounding applied.
 * @return The rounded color value in [0,65535]
 */
static inline vint4 apply_u8_rounding(
	vint4 colori,
	vmask4 u8_mask
) {
	// The real decoder would just use the top 8 bits, but we rescale
	// into a 16-bit value that rounds correctly.
	vint4 colori_u8 = asr<8>(colori) * 257;
	return select(colori, colori_u8, u8_mask);
}

/**
 * @brief Decode symbolic ASTC block to internal value range.
 *
 * We want to do our error measurement the way the encoder does, which is to say, in
 * the internal working 16-bit integer value space that the decoder produces before
 * the final output conversion.
 *
 * Regular "decompress_symbolic_block" does the final output conversion and produces
 * a deinterleaved RGBA result; for our purposes, interleaved RGBA without the final
 * conversion is actually better.
 *
 * @param config         The encoder config to use.
 * @param bsd            The block size descriptor.
 * @param scb            The symbolic compressed block we want to decode.
 * @param[out] out_rgba  Where to write the output texels, 4 floats (R, G, B, A) per texel.
 */
static void decompress_symbolic_block_raw(
	const astcenc_config& config,
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& scb,
	float* out_rgba
) {
	// The unorm8 rounding for the decode.
	vmask4 u8_mask(false);
	if (config.flags & ASTCENC_FLG_USE_DECODE_UNORM8)
	{
		u8_mask = vmask4(true);
	}
	else if (config.profile == ASTCENC_PRF_LDR_SRGB)
	{
		u8_mask = vmask4(true, true, true, false);
	}

	if (scb.block_type != SYM_BTYPE_NONCONST)
	{
		// Error and constant blocks output the same value for every texel
		// for error and SYM_BTYPE_CONST_F16 blocks (which are not legal in the profiles this
		// module supports), produce values far outside the nominal [0,65535] range.
		vfloat4 color(-65536.0f);

		if (scb.block_type == SYM_BTYPE_CONST_U16)
		{
			vint4 colori(scb.constant_color);

			colori = apply_u8_rounding(colori, u8_mask);
			color = int_to_float(colori);
		}

		// Output the same color value for every pixel in the block
		for (unsigned int i = 0; i < bsd.texel_count; i++)
		{
			store(color, out_rgba + i * 4);
		}

		return;
	}

	// Get the appropriate partition-table entry
	int partition_count = scb.partition_count;
	const auto& pi = bsd.get_partition_info(partition_count, scb.partition_index);

	// Get the appropriate block descriptors
	const auto& bm = bsd.get_block_mode(scb.block_mode);
	const auto& di = bsd.get_decimation_info(bm.decimation_mode);

	bool is_dual_plane = static_cast<bool>(bm.is_dual_plane);

	// Unquantize and undecimate the weights
	int plane1_weights[BLOCK_MAX_TEXELS];
	int plane2_weights[BLOCK_MAX_TEXELS];
	unpack_weights(bsd, scb, di, is_dual_plane, plane1_weights, plane2_weights);

	// Now that we have endpoint colors and weights, we can unpack texel colors
	int plane2_component = scb.plane2_component;
	vmask4 plane2_mask = vint4::lane_id() == vint4(plane2_component);

	for (int i = 0; i < partition_count; i++)
	{
		// Decode the color endpoints for this partition
		vint4 ep0;
		vint4 ep1;
		bool rgb_lns;
		bool a_lns;

		unpack_color_endpoints(config.profile,
			scb.color_formats[i],
			scb.color_values[i],
			rgb_lns, a_lns,
			ep0, ep1);

		vint4 diff = ep1 - ep0;
		int texel_count = pi.partition_texel_count[i];

		for (int j = 0; j < texel_count; j++)
		{
			unsigned int tix = pi.texels_of_partition[i][j];
			vint4 weight = select(vint4(plane1_weights[tix]), vint4(plane2_weights[tix]), plane2_mask);

			// Interpolate the color
			vint4 colori = ep0 + asr<6>(diff * weight + vint4(32));
			colori = apply_u8_rounding(colori, u8_mask);

			store(int_to_float(colori), out_rgba + tix * 4);
		}
	}
}

/**
  * @brief Block error calculation from physical block bits.
  *
  * We keep a cache (hash table) of pre-decoded RGBA blocks for input block bit patterns
  * since the block-splicing we do frequently ends up making the same block multiple ways,
  * and decoding a block from physical bits is fairly expensive.
  *
  * The hash table is organized in "rows". The hash code of the block bytes determines the row,
  * each row has multiple entries with FIFO eviction inside a row.
  *
  * Rows also keep tag bytes. Each entry has a corresponding tag byte storing some extra hash
  * bits (to catch most hash collisions before looking at "entries"), and we use one extra byte per
  * hash table row to keep track of the FIFO evict position.
  */
struct block_error_calculator
{
	/** @brief Determines row size. Rows have (1 << ROW_SHIFT) elements. */
	static const unsigned int ROW_SHIFT = 3;

	/** @brief Items per row, derived from ROW_SHIFT. */
	static const unsigned int ITEMS_PER_ROW = 1 << ROW_SHIFT;

	/** @brief Tag bytes per row, including the FIFO counter. */
	static const unsigned int TAGS_PER_ROW = ITEMS_PER_ROW + 1;

	/** @brief The block size descriptor used. */
	const block_size_descriptor* bsd;

	/** @brief The encoder configuration. */
	const astcenc_config& config;

	/** @brief The channel weights to use during the error calculation. */
	vfloat4 channel_weights;

	/** @brief The number of rows in the hash table. Must be a power of 2. */
	unsigned int row_count;

	/** @brief Number of floats worth of pixel data stored per hash. */
	unsigned int floats_per_entry;

	/** @brief The tag bytes for the hash table. */
	uint8_t* tags;

	/** @brief The block bits for the hash table entries. */
	physical_block* entries;

	/** @brief The decoded pixel values. floats_per_entry values per hash table entry. */
	float* pixels;

	/** @brief Pointer to the block texel weights. */
	const float* block_texel_weights;

	/** @brief Interleaved RGBA values for the original pixels in the current block. */
	float original_rgba[BLOCK_MAX_TEXELS * 4];

	/**
	 * @brief Initializes the block error calculator.
	 *
	 * @param in_config        The encoder configuration being used.
	 * @param num_entries      The number of entries in the hash table (must be pow2).
	 * @param in_bsd           The block size descriptor.
	 */
	block_error_calculator(
		const astcenc_config& in_config,
		unsigned int num_entries,
		const block_size_descriptor* in_bsd
	) : bsd(in_bsd), config(in_config)
	{
		assert((num_entries & (num_entries - 1)) == 0);
		assert(num_entries >= ITEMS_PER_ROW);
		row_count = num_entries >> ROW_SHIFT;
		floats_per_entry = bsd->texel_count * 4;

		unsigned int entry_count = row_count * ITEMS_PER_ROW;
		unsigned int tag_count = row_count * TAGS_PER_ROW;

		tags = new uint8_t[tag_count];
		entries = new physical_block[entry_count];
		pixels = new float[entry_count * floats_per_entry];
		block_texel_weights = nullptr;

		// Clear all the tags (marking entries as invalid)
		std::memset(tags, 0, tag_count * sizeof(*tags));

	}

	~block_error_calculator()
	{
		delete[] tags;
		delete[] entries;
		delete[] pixels;
	}

	/**
	 * @brief Changes the current active block being encoded.
	 *
	 * @param blk      Source pixel data in image_block format.
	 * @param weights  The per-texel weights.
	 */
	void set_current_block(
		const image_block& blk,
		const float* weights
	) {
		// Convert to interleaved texels
		for (unsigned int i = 0; i < bsd->texel_count; i++)
		{
			original_rgba[i * 4 + 0] = blk.data_r[i];
			original_rgba[i * 4 + 1] = blk.data_g[i];
			original_rgba[i * 4 + 2] = blk.data_b[i];
			original_rgba[i * 4 + 3] = blk.data_a[i];
		}

		channel_weights = blk.channel_weight;
		block_texel_weights = weights;
	}

	/**
	 * @brief Evaluate the error for the given input block against the current target pixels.
	 *
	 * @param candidate_bits       The block being evaluated.
	 */
	float eval(const physical_block& candidate_bits)
	{
		unsigned int hash = candidate_bits.hash();

		// Use low bits of hash to determine row index
		// and high bits to determine tags
		uint8_t our_tag = (uint8_t) ((hash >> 25) | 0x80);
		unsigned int row = hash & (row_count - 1);
		float* matched_pixels = nullptr;

		uint8_t* row_tags = tags + row * TAGS_PER_ROW;
		physical_block* row_entries = entries + row * ITEMS_PER_ROW;

		// Look for matches in current hash row
		for (unsigned int i = 0; i < ITEMS_PER_ROW; i++)
		{
			if (row_tags[i] == our_tag && row_entries[i] == candidate_bits)
			{
				matched_pixels = pixels + (row * ITEMS_PER_ROW + i) * floats_per_entry;
				break;
			}
		}

		// If we got no match, insert a new entry
		if (!matched_pixels)
		{
			uint8_t lru_slot = row_tags[ITEMS_PER_ROW] & (ITEMS_PER_ROW - 1);
			row_tags[ITEMS_PER_ROW]++;

			row_tags[lru_slot] = our_tag;
			row_entries[lru_slot] = candidate_bits;
			matched_pixels = pixels + (row * ITEMS_PER_ROW + lru_slot) * floats_per_entry;

			// Turn the physical encoding back to symbolic
			symbolic_compressed_block scb;
			physical_to_symbolic(*bsd, candidate_bits.bytes, scb);

			// Decompress the block in raw form (without applying final color transform)
			decompress_symbolic_block_raw(config, *bsd, scb, matched_pixels);
		}

		return calculate_ssd_weighted(original_rgba, matched_pixels, bsd->texel_count, block_texel_weights, channel_weights);
	}
};

/** @brief A list of the best few candidate blocks in the running. */
struct candidate_list
{
	/** @brief One of our current best blocks. */
	struct item
	{
		/** @brief The physical block. */
		physical_block bits;
		/** @brief Rate-distortion cost. */
		float rd_cost;
		/** @brief MTF position for this candidate. */
		int mtf_position;
		/** @brief Number of weight bits for this candidate. */
		int weight_bits;  // Weight bits
	};

	/** @brief The list of best candidates, in ascending order of rd_cost. */
	item list[BEST_CANDIDATES_COUNT];

	/** @brief The current number of candidates in the list. */
	unsigned int count = 0;

	/**
	 * @brief Quick rejection check for candidates.
	 *
	 * If the candidate list is full and the provided cost lower bound is
	 * bigger than our current worst, can skip further evaluation.
	 *
	 * @param cost_lower_bound  The cost lower bound to check for a candidate.
	 */
	bool quick_reject(float cost_lower_bound) const
	{
		// If all candidate slots are currently populated and this is worse than the
		// current worst, don't bother
		return count == BEST_CANDIDATES_COUNT && cost_lower_bound >= list[BEST_CANDIDATES_COUNT - 1].rd_cost;
	}

	/**
	 * @brief Try adding a new candidate to the list.
	 *
	 * Candidates are kept in order of increasing cost.
	 *
	 * @param bits          The candidate encoding.
	 * @param rd_cost       The encoding's rate-distortion cost.
	 * @param mtf_position  The MTF list position for the endpoint/weight set used to make this candidate.
	 * @param bsd           The block size descriptor.
	 */
	void add(
		const physical_block& bits,
		float rd_cost,
		int mtf_position,
		const block_size_descriptor* bsd
	) {
		if (quick_reject(rd_cost))
		{
			return;
		}

		unsigned int insert_pos = count < BEST_CANDIDATES_COUNT ? count : BEST_CANDIDATES_COUNT - 1;

		// Find the position to insert, shifting elements out of the way
		while (insert_pos > 0 && rd_cost < list[insert_pos - 1].rd_cost)
		{
			list[insert_pos] = list[insert_pos - 1];
			insert_pos--;
		}

		// Determine number of weight bits
		int weight_bits = get_weight_bits(bits.bytes, bsd);

		// Insert the candidate into the list
		list[insert_pos] = {bits, rd_cost, mtf_position, weight_bits};

		// Increment count if not yet at capacity
		if (count < BEST_CANDIDATES_COUNT)
		{
			count++;
		}
	}

	/** @brief Iterator begin() for range-based loops. */
	const item* begin() const
	{
		return list;
	}

	/** @brief Iterator end() for range-based loops. */
	const item* end() const
	{
		return list + count;
	}
};

/**
 * @brief Compact representation of a simple symbolic block for endpoint substitution,
 * without its weights.
 *
 * "Simple" blocks are non-constant, 1-plane, 1-partition blocks
 */
struct compact_simple_endpoints
{
	/** @brief The block mode. */
	uint16_t block_mode;

	/** @brief The endpoint color quant mode (actually quant_method enum). */
	uint8_t quant_mode;

	/** @brief The endpoint color format for the single partition. */
	uint8_t color_format;

	/** @brief The quantized endpoint color pairs. */
	uint8_t color_values[8];

	/** @brief Initialize from an appropriate symbolic_compressed_block. */
	void from_symbolic(const symbolic_compressed_block& blk)
	{
		assert(blk.block_type == SYM_BTYPE_NONCONST);
		assert(blk.partition_count == 1);
		// Even though color formats in 1-partition blocks are "matched" (by default),
		// physical_to_symbolic flags them as 0, so that's what we do.
		assert(blk.color_formats_matched == 0);
		assert(blk.plane2_component == -1);

		block_mode = blk.block_mode;
		color_format = blk.color_formats[0];
		quant_mode = static_cast<uint8_t>(blk.quant_mode);

		// Make the color values we encode canonical so we can do equality compares on them
		// (by setting unused values to 0)
		unsigned int num_color_values = (color_format >> 2) * 2 + 2;

		for (unsigned int i = 0; i < 8; i++)
		{
			uint8_t value = blk.color_values[0][i];
			color_values[i] = (i < num_color_values) ? value : 0;
		}
	}

	/** @brief Convert back to symbolic_compressed_block form. */
	void to_symbolic(symbolic_compressed_block& blk) const
	{
		blk.block_type = SYM_BTYPE_NONCONST;
		blk.partition_count = 1;
		blk.color_formats_matched = 0;
		blk.plane2_component = -1;
		blk.block_mode = block_mode;
		blk.color_formats[0] = color_format;
		blk.quant_mode = static_cast<quant_method>(quant_mode);
		std::memcpy(blk.color_values[0], color_values, sizeof(blk.color_values[0]));
	}

	/** @brief Compare two compact simple endpoints for equality. */
	bool operator ==(const compact_simple_endpoints& x) const
	{
		if (block_mode != x.block_mode ||
			quant_mode != x.quant_mode ||
			color_format != x.color_format)
		{
			return false;
		}

		uint64_t colors, colors_x;
		static_assert(sizeof(color_values) == sizeof(uint64_t), "Need sizes to match");
		memcpy(&colors, color_values, sizeof(color_values));
		memcpy(&colors_x, x.color_values, sizeof(x.color_values));
		return colors == colors_x;
	}

	/**
	 * @brief Determine normalized endpoint axis direction.
	 *
	 * Used for an early-out to quickly rule out endpoint pairs.
	 */
	vfloat4 normalized_axis() const
	{
		vint4 color0(color_values + 0);
		vint4 color1(color_values + 4);
		return normalize_safe(int_to_float(color0 - color1), vfloat4(0.0f));
	}
};

/**
 * @brief Endpoint substitution move-to-front list.
 */
struct endpoint_subst_mtf
{
	/** @brief An entry in the substitution list. */
	struct entry
	{
		/** @brief The endpoints and block mode etc., in compact form. */
		compact_simple_endpoints endpoints;
		/** @brief Cached normalized_axis() for this candidate. */
		vfloat4 normalized_axis;
	};

	/** @brief The list of entries in the list. Entries don't move after inserting. */
	entry* entries;
	/** @brief List of indices defining an index permutation. We shuffle indices instead of elements. */
	unsigned int* order;
	/** @brief Number of elements currently in the move-to-front list. */
	unsigned int size;
	/** @brief Maximum capacity of move-to-front list. */
	unsigned int capacity;

	/** @brief Initializes the list with a given capacity. */
	explicit endpoint_subst_mtf(unsigned int in_capacity)
		: size(0), capacity(in_capacity)
	{
		entries = new entry[capacity];
		order = new unsigned int[capacity];

		// initialize order table to identity map
		for (unsigned int i = 0; i < capacity; i++)
		{
			order[i] = i;
		}
	}

	~endpoint_subst_mtf()
	{
		delete[] entries;
		delete[] order;
	}

	/** @brief Reset the list back to empty. */
	void reset()
	{
		size = 0;
	}

	/** @brief Indexes entries in their logical order. */
	const entry& operator[](size_t index) const
	{
		assert(index < size);
		return entries[order[index]];
	}

	/** @brief Searches for a compact endpoint encoding in the list. */
	int search(const compact_simple_endpoints& endpoints)
	{
		for (unsigned int i = 0; i < size; i++)
		{
			if (entries[order[i]].endpoints == endpoints)
			{
				return (int)i;
			}
		}

		return -1;
	}

	/**
	 * @brief Update the MTF list after using an element.
	 *
	 * @param mtf_index    The MTF index (search() result) for the new element, or -1 if it wasn't a match.
	 * @param endpoints    The new endpoints. Ignored unless mtf_index == -1.
	 */
	void update(
		int mtf_index,
		const compact_simple_endpoints& endpoints
	) {
		unsigned int id;

		if (mtf_index == -1)
		{
			// New entry. If we're not yet at capacity, grow the list, else replace oldest.
			if (size < capacity)
			{
				size++;
			}

			mtf_index = size - 1;
			id = order[mtf_index];
			entries[id].endpoints = endpoints;
			entries[id].normalized_axis = endpoints.normalized_axis();
		}
		else
		{
			id = order[mtf_index];
		}

		// Move accessed block to front
		for (int i = mtf_index; i > 0; i--)
		{
			order[i] = order[i - 1];
		}
		order[0] = id;
	}

	/** @brief Combined search and update. */
	void encode(const compact_simple_endpoints& endpoints)
	{
		update(search(endpoints), endpoints);
	}
};

/** @brief Determines whether a physical block is "simple" (1-plane, 1-partition). */
static bool is_simple_block(
	const uint8_t* bytes,
	const block_size_descriptor* bsd
) {
	int which_mode = get_block_mode(bytes);

	unsigned int npart = (bytes[1] >> 3) & 3;
	if (npart != 0)
	{
		return false;
	}

	const block_mode& bm = bsd->get_block_mode(which_mode);
	return !bm.is_dual_plane;
}

/**
 * @brief Determines a set of weights for a block given known mode and endpoints.
 *
 * Only supports "simple" blocks (1-plane, 1-partition).
 *
 * @param profile    The encoder profile being used.
 * @param blk        The block to be encoded.
 * @param scb        The symbolic compressed block holding the mode and endpoints.
 * @param[out] ei    The output weights (as well as endpoints in the expected form) get written here.
 */
static void compute_known_endpoint_weights(
	astcenc_profile profile,
	const image_block& blk,
	const symbolic_compressed_block& scb,
	endpoints_and_weights& ei
) {
	bool rgb_hdr, alpha_hdr;
	vint4 ep0, ep1;

	unpack_color_endpoints(
		profile,
		scb.color_formats[0],
		scb.color_values[0],
		rgb_hdr, alpha_hdr,
		ep0, ep1);

	ei.ep.partition_count = 1;
	ei.ep.endpt0[0] = int_to_float(ep0);
	ei.ep.endpt1[0] = int_to_float(ep1);

	vfloat4 basef = ei.ep.endpt0[0];
	vfloat4 dirf = ei.ep.endpt1[0] - basef;

	float length_squared = dot_s(dirf, dirf);

	if (length_squared > 1e-7f)
	{
		// Scale direction to cancel out squared length
		dirf = dirf * vfloat4(1.0f / length_squared);
	}

	// NOTE: compare compute_ideal_colors_and_weights_4_comp;
	// this path normalizes channel weights so hadd_s(blk.channel_weight) == 4.0f always
	const float error_weight = 1.0f;
	const float wes = length_squared * error_weight;

	for (size_t i = 0; i < blk.texel_count; i++)
	{
		vfloat4 point = blk.texel(i);
		float param = dot_s(point - basef, dirf);

		ei.weights[i] = astc::clamp1f(param);
		ei.weight_error_scale[i] = wes;
	}

	// Zero initialize any SIMD over-fetch
	size_t texel_count_simd = round_up_to_simd_multiple_vla(blk.texel_count);
	for (size_t i = blk.texel_count; i < texel_count_simd; i++)
	{
		ei.weights[i] = 0.0f;
		ei.weight_error_scale[i] = 0.0f;
	}

	ei.is_constant_weight_error_scale = true;
}

/**
 * @brief Runs the LZ-RDO optimization passes.
 *
 * This is the meat of the LZ-RDO encode. Over the course of potentially multiple passes,
 * try encoding blocks of the image by reusing previous endpoints or weights, using the
 * result if they score favorably in our combined rate-distortion metric.
 *
 * @param ctxo                 The compression context.
 * @param thread_index         Index for the thread this runs on.
 * @param data                 Pointer to the encoded blocks, initially contains results of regular ASTC encode.
 * @param lambda               Lambda parameter for rate-distortion optimization.
 * @param init_pct             Initial progress callback percentage.
 * @param img                  The image being encoded.
 * @param swz                  The component swizzle for the image being encoded.
 */
static void lz_rdo_optimization_passes(
	astcenc_context& ctxo,
	size_t thread_index,
	uint8_t* data,
	float lambda,
	float init_pct,
	const astcenc_image& img,
	const astcenc_swizzle& swz
) {
	astcenc_contexti& ctx = ctxo.context;
	const astcenc_config& config = ctx.config;
	const block_size_descriptor* bsd = ctx.bsd;

	float* per_texel_weights = ctx.lz_rdo_per_texel_weights;
	astcenc_profile profile = config.profile;

	// HDR not currently supported. Just assert() here, user-level validation
	// for API callers happens in astc_entry.cpp validate_flags().
	assert(profile == ASTCENC_PRF_LDR_SRGB || profile == ASTCENC_PRF_LDR);

	size_t blocks_x = ceil_div(img.dim_x, config.block_x);
	size_t blocks_y = ceil_div(img.dim_y, config.block_y);
	size_t blocks_z = ceil_div(img.dim_z, config.block_z);

	// Stuff for encoding
	auto& tmpbuf = ctx.working_buffers[thread_index];
	// can look into HDR opt later
	auto load_func = load_image_block;
	image_block blk;

	// Set up block re-encoding parameters
	blk.texel_count = bsd->texel_count;
	blk.decode_unorm8 = config.flags & ASTCENC_FLG_USE_DECODE_UNORM8;

	// Allocate thread-local resources
	uint8_t* original_blocks = new uint8_t[MAX_BLOCKS_PER_ITEM * ASTC_BLOCK_BYTES];
	mtf_list mtf_weights { config.tune_lz_rdo_weight_history_size }; // MTF list for weights
	mtf_list mtf_endpoints { config.tune_lz_rdo_endpoint_history_size * 2 }; // MTF list for endpoints and constant block bits
	endpoint_subst_mtf mtf_simple_endpoints { config.tune_lz_rdo_endpoint_history_size }; // MTF list for simple endpoints
	mode_byte_histogram hist; // Histogram for mode bytes.

	// Determine channel weights
	blk.channel_weight = vfloat4 { config.cw_r_weight, config.cw_g_weight, config.cw_b_weight, config.cw_a_weight };

	// Unlike normal encoding, we care about the absolute scale of errors, not just their relative values, so
	// normalize channel_weight to always have a sum of 4 (which is what we get with the default of all-1s)
	blk.channel_weight = blk.channel_weight * vfloat4(4.0f / hadd_s(blk.channel_weight));

	// Set up the error calculator
	block_error_calculator err_calc(config, CACHE_SIZE, bsd);

	ctxo.manage_lz_rdo_optimize.init(blocks_x * blocks_y * blocks_z, config.progress_callback, init_pct, 100.0f);

	// Main thread loop to process work items from the queue
	while (true)
	{
		size_t block_count;
		size_t block_start = ctxo.manage_lz_rdo_optimize.get_task_assignment(MAX_BLOCKS_PER_ITEM, block_count);
		if (!block_count)
			break;

		// Save the original block encodings we had so far
		std::memcpy(original_blocks, data + block_start * ASTC_BLOCK_BYTES, block_count * ASTC_BLOCK_BYTES);

		// Set up pointers to other block encodings for this slice of blocks
		uint8_t* other_encodings[2];
		other_encodings[0] = ctx.lz_rdo_restricted_blocks + block_start * ASTC_BLOCK_BYTES;
	 	other_encodings[1] = original_blocks;

		// Initialize current encodings for this slice to start with the restricted encodings
		std::memcpy(data + block_start * ASTC_BLOCK_BYTES, other_encodings[0], block_count * ASTC_BLOCK_BYTES);

		// Outer pass loop
		for (unsigned int pass = 0; pass < config.tune_lz_rdo_num_passes; pass++)
		{
			// Is current pass forward or not? (Odd passes go over blocks in reverse order.)
			bool is_forward = (pass & 1) == 0;

			// Reset MTF lists and histogram
			mtf_weights.reset();
			mtf_endpoints.reset();
			mtf_simple_endpoints.reset();
			hist.reset();

			// Seed the structures with random blocks
			// Use block start index and pass index as part of the seed for variety between chunks
			// multiply by large prime for better mixing
			uint32_t rng_state = (static_cast<uint32_t>(block_start) * 3677199193u + pass * 7u + 23857527u) | 1u;
			constexpr size_t MAX_SAMPLES = 64;
			const size_t num_samples = astc::min(block_count, MAX_SAMPLES);

			for (size_t i = 0; i < num_samples; i++)
			{
				unsigned int rng = xorshift32(rng_state);

				// Draw a random block index within our group
				size_t block_idx = (num_samples == block_count) ? i : (rng % block_count);

				// Get block data
				physical_block block_bits;

				if (pass == 0)
				{
					// Initial pass draws from seed blocks
					unsigned int which = (rng >> 31) & 1;
					block_bits = physical_block(other_encodings[which] + block_idx * ASTC_BLOCK_BYTES);
				}
				else
				{
					// Other passes draw from existing blocks
					block_bits = physical_block(data + (block_idx + block_start) * ASTC_BLOCK_BYTES);
				}

				int block_weight_bits = get_weight_bits(block_bits.bytes, bsd);
				physical_block weights_mask = physical_block::top_bits_mask(block_weight_bits);

				// Update structures
				mtf_weights.encode(block_bits, weights_mask);
				mtf_endpoints.encode(block_bits, ~weights_mask);
				hist.update(block_bits);

				if (is_simple_block(block_bits.bytes, bsd))
				{
					symbolic_compressed_block scb;
					physical_to_symbolic(*bsd, block_bits.bytes, scb);

					compact_simple_endpoints compact;
					compact.from_symbolic(scb);

					mtf_simple_endpoints.encode(compact);
				}
			}

			// Process the blocks in either forward or reverse order
			for (size_t block_iter = 0; block_iter < block_count; block_iter++)
			{
				size_t block_index = block_start + (is_forward ? block_iter : block_count - 1 - block_iter);

				// Calculate block coordinates
				size_t x, y, z;
				decode_index_to_coords(block_index, blocks_x, blocks_y, x, y, z);

				// Get current block data and compute its weight bits
				uint8_t* current_block = data + block_index * ASTC_BLOCK_BYTES;
				physical_block current_bits(current_block);
				int current_weight_bits = get_weight_bits(current_bits.bytes, bsd);
				physical_block best = current_bits;

				// Don't process blocks with no weight bits, accept void-extent as is
				if (current_weight_bits == 0)
				{
					mtf_weights.encode(current_bits, physical_block {});
					mtf_endpoints.encode(current_bits, ~physical_block {});
					hist.update(current_bits);
					continue;
				}

				// Get the source pixels
				load_func(profile, img, blk, *bsd, x * bsd->dim_x, y * bsd->dim_y, z * bsd->dim_z, swz);

				// Set up error calculation for the current block
				const float* block_texel_weights = per_texel_weights + block_index * bsd->texel_count;
				err_calc.set_current_block(blk, block_texel_weights);

				// Decode the original block to compute initial error, as well as the error for the restricted approximation
				float original_err = err_calc.eval(current_bits);
				float restricted_err = ERROR_CALC_DEFAULT;

				// Calculate masks for weights and endpoints
				physical_block current_weights_mask = physical_block::top_bits_mask(current_weight_bits);
				int mtf_weights_pos = mtf_weights.search(current_bits, current_weights_mask);
				int mtf_endpoints_pos = mtf_endpoints.search(current_bits, ~current_weights_mask);

				// Figure out rate to determine initial RD cost
				float original_rate = calculate_bit_cost_simple(mtf_endpoints_pos, mtf_weights_pos, current_bits, current_weight_bits, hist);
				float best_rd_cost = original_err + lambda * original_rate;

				// Candidate structure for storing best candidates
				candidate_list best_endpoints;

				// Add the current block to the candidates
				best_endpoints.add(current_bits, best_rd_cost, mtf_endpoints_pos, bsd);

				// Try the other candidate encodings
				for (const uint8_t* other_encoding : other_encodings)
				{
					physical_block other_bits(other_encoding + (block_index - block_start) * ASTC_BLOCK_BYTES);
					int other_weight_bits = get_weight_bits(other_bits.bytes, bsd);
					physical_block other_weight_mask = physical_block::top_bits_mask(other_weight_bits);

					int mtf_weights_pos_other = mtf_weights.search(other_bits, other_weight_mask);
					int mtf_endpoints_pos_other = mtf_endpoints.search(other_bits, ~other_weight_mask);

					float other_err = err_calc.eval(other_bits);
					if (restricted_err == ERROR_CALC_DEFAULT)
					{
						restricted_err = other_err;
					}

					float other_rate = calculate_bit_cost_simple(mtf_endpoints_pos_other, mtf_weights_pos_other, other_bits, other_weight_bits, hist);
					float other_rd_cost = other_err + lambda * other_rate;

					if (other_rd_cost < best_rd_cost)
					{
						best = other_bits;
						best_rd_cost = other_rd_cost;
					}

					best_endpoints.add(other_bits, other_rd_cost, mtf_endpoints_pos_other, bsd);
				}

				// Find best simple endpoint candidates
				const uint8_t* simple_ref_block = other_encodings[0] + (block_index - block_start) * ASTC_BLOCK_BYTES;
				if (is_simple_block(simple_ref_block, bsd))
				{
					symbolic_compressed_block scb;
					compact_simple_endpoints ref_compact;

					// Decode reference block to determine target normalized_axis
					physical_to_symbolic(*bsd, simple_ref_block, scb);
					ref_compact.from_symbolic(scb);

					vfloat4 target_axis = ref_compact.normalized_axis();
					float best_subst_err = ERROR_CALC_DEFAULT;

					for (size_t k = 0; k < mtf_simple_endpoints.size; k++)
					{
						const endpoint_subst_mtf::entry& entry = mtf_simple_endpoints[k];

						// Don't bother with trying endpoint substitution if the candidate's endpoint
						// axis isn't at least somewhat aligned with our reference block
						//
						// This uses restricted_err (error of the unconstrained 1-plane 1-subset encoding) as a
						// proxy for how well these types of blocks can get. If the axes align perfectly (absolute
						// value of dot product is 1), we can expect to get close to restricted_err. If they don't,
						// we expect to be worse.
						//
						// Somewhat arbitrarily, bill bad dot product matches (dot product 0), as 2x the restricted_err,
						// and perfect matches as 1x the restricted_err; such candidates are interesting to evaluate if
						// they're not much worse than best_subst_err, the best substitution we've found so far. They
						// are allowed to be a bit worse, because we track multiple candidates.
						float quick_error_estimate = (2.0f - astc::fabs(dot_s(entry.normalized_axis, target_axis))) * restricted_err;
						if (quick_error_estimate >= best_subst_err * 1.25f)
						{
							continue;
						}

						entry.endpoints.to_symbolic(scb);
						const block_mode& bm = bsd->get_block_mode(scb.block_mode);

						// Compute weights for our block using those endpoints
						endpoints_and_weights& ei = tmpbuf.ei1;
						compute_known_endpoint_weights(profile, blk, scb, ei);

						// Compute ideal weights for our decimation mode
						float* dec_weights_ideal = tmpbuf.dec_weights_ideal;
						const auto& di = bsd->get_decimation_info(bm.decimation_mode);
						compute_ideal_weights_for_decimation(ei, di, dec_weights_ideal);

						// Quantize them
						uint8_t* dec_weights_uquant = tmpbuf.dec_weights_uquant;
						compute_quantized_weights_for_decimation(
							di,
							0.0f, 1.0f,
							dec_weights_ideal,
							tmpbuf.dec_weights_ideal + BLOCK_MAX_WEIGHTS, // ignored, just pass something
							dec_weights_uquant,
							bm.get_weight_quant_mode());

						// Copy quantized weights into scb
						for (unsigned int wi = 0; wi < di.weight_count; wi++)
						{
							scb.weights[wi] = dec_weights_uquant[wi];
						}

						// Determine the corresponding error
						float mse = compute_symbolic_block_difference_1plane_1partition_weighted(
							ctx.config,
							*bsd,
							scb,
							blk,
							block_texel_weights
						);

						// No need to perform RD scoring if we blew past our target on the distortion term alone
						if (best_endpoints.quick_reject(mse))
						{
							continue;
						}

						// Keep track of best substitution we've found so far
						best_subst_err = astc::min(best_subst_err, mse);

						// Emit the corresponding block so we have the bits for the rate calc
						physical_block candidate;
						symbolic_to_physical(*bsd, scb, candidate.bytes);

						// Calculate the rate
						int endpoints_weight_bits = get_weight_bits(candidate.bytes, bsd);
						physical_block weights_mask = physical_block::top_bits_mask(endpoints_weight_bits);

						// Find the corresponding position in the MTF lists.
						// If it's not actually in the endpoint MTF list anymore, skip.
						// (This is possible if we e.g. have a long run of distinct constant-color
						// blocks that crowd out the matched endpoints, but very unlikely.)
						int endpoint_pos = mtf_endpoints.search(candidate, ~weights_mask);
						if (endpoint_pos == -1)
						{
							continue;
						}

						int weight_pos = mtf_weights.search(candidate, weights_mask);

						float bit_cost = calculate_bit_cost_simple(endpoint_pos, weight_pos, candidate, endpoints_weight_bits, hist);
						float rd_cost = mse + lambda * bit_cost;
						if (rd_cost < best_rd_cost)
						{
							best = candidate;
							best_rd_cost = rd_cost;
						}

						// Insert into best_endpoints if it's one of the best candidates
						best_endpoints.add(candidate, rd_cost, endpoint_pos, bsd);
					}
				}

				// Find best weight candidates
				for (size_t k = 0; k < mtf_weights.size; k++)
				{
					physical_block candidate_weights = mtf_weights.list[k];
					int weights_weight_bits = get_weight_bits(candidate_weights.bytes, bsd);
					if (weights_weight_bits == 0)
					{
						continue;
					}

					physical_block weights_mask = physical_block::top_bits_mask(weights_weight_bits);
					physical_block just_weight_bits = candidate_weights & weights_mask;
					physical_block endpoints_mask = ~weights_mask;

					// Try every endpoint candidate that matches in weight bits
					for (const candidate_list::item& endpoints : best_endpoints)
					{
						if (weights_weight_bits == endpoints.weight_bits)
						{
							physical_block combined_bits = just_weight_bits | (endpoints.bits & endpoints_mask);

							float err = err_calc.eval(combined_bits);
							float bit_cost = calculate_bit_cost_simple(endpoints.mtf_position, static_cast<int>(k), combined_bits, weights_weight_bits, hist);
							float rd_cost = err + lambda * bit_cost;
							if (rd_cost < best_rd_cost)
							{
								best = combined_bits;
								best_rd_cost = rd_cost;
							}
						}
					}
				}

				// Write back best candidate
				std::memcpy(current_block, best.bytes, 16);

				// Recalculate masks for the best match
				int best_weight_bits = get_weight_bits(current_block, bsd);
				physical_block best_weights_mask = physical_block::top_bits_mask(best_weight_bits);

				// Update histogram with literal mask
				int best_mtf_weights_pos = mtf_weights.search(best, best_weights_mask);
				int best_mtf_endpoints_pos = mtf_endpoints.search(best, ~best_weights_mask);

				// Update statistics
				mtf_weights.update(best, best_mtf_weights_pos);
				mtf_endpoints.update(best, best_mtf_endpoints_pos);
				hist.update(best);

				if (is_simple_block(best.bytes, bsd))
				{
					// Simple blocks can be compactly encoded
					symbolic_compressed_block scb;
					physical_to_symbolic(*bsd, best.bytes, scb);

					compact_simple_endpoints compact;
					compact.from_symbolic(scb);

					mtf_simple_endpoints.encode(compact);
				}
			}
		}

		ctxo.manage_lz_rdo_optimize.complete_task_assignment(block_count);
	}

	delete[] original_blocks;
	ctxo.manage_lz_rdo_optimize.wait();
}

/** @brief Parameters for a convolution filter kernel. */
struct filter_kernel
{
	/** @brief Radius of the kernel. We have 2*radius + 1 taps. */
	int radius;
	/** @brief The filter coefficients. */
	float coeffs[15];
};

/** @brief Gaussian kernel with sigma=2.20. */
static const filter_kernel kernel_initial_lpf {
	7,
	{1.14907966e-03f, 4.40146100e-03f, 1.37123950e-02f, 3.47455753e-02f, 7.16069925e-02f, 1.20027593e-01f, 1.63635110e-01f, 1.81443588e-01f,
	 1.63635110e-01f, 1.20027593e-01f, 7.16069925e-02f, 3.47455753e-02f, 1.37123950e-02f, 4.40146100e-03f, 1.14907966e-03f}
};

/** @brief Gaussian kernel with sigma=1.25. */
static const filter_kernel kernel_spread {
	4,
	{1.90769133e-03f, 1.79195767e-02f, 8.87562444e-02f, 2.31804370e-01f, 3.19224234e-01f, 2.31804370e-01f, 8.87562444e-02f, 1.79195767e-02f,
	 1.90769133e-03f}
};

/**
 * @brief Apply a convolution filter along 1D slices of a 3D array.
 *
 * All arrays are xdim * ydim * zdim values and are expected to have
 * at least ASTCENC_SIMD_WIDTH floats of padding at the end.
 *
 * @param input        The input values.
 * @param[out] output  Where output values get written.
 * @param xdim         Size of the array in the X dimension.
 * @param ydim         Size of the array in the Y dimension.
 * @param zdim         Size of the array in the Z dimension.
 * @param filter       The filter kernel to use.
 * @param direction    The direction to work in (0=X, 1=Y, 2=Z).
 */
static void apply_1d_convolution_3d(
	float* input,
	float* output,
	size_t xdim,
	size_t ydim,
	size_t zdim,
	const filter_kernel& filter,
	int direction
) {
	int radius = filter.radius;
	size_t radius_u = filter.radius;
	const float* coeffs = filter.coeffs + radius;

	if (direction == 0)
	{
		// Direction is X.

		// The first and last radius_u pixels in each line include
		// border pixels and need careful handling.
		//
		// Additionally, shrink safe region on the right by ASTCENC_SIMD_WIDTH - 1
		// to account for SIMD width (we need the rightmost pixel of the span to
		// still be inside the safe region).
		size_t excluded_region_width = 2 * radius_u + ASTCENC_SIMD_WIDTH - 1;

		// From this, compute safe region width (clamps at 0).
		size_t safe_width = xdim - astc::min(excluded_region_width, xdim);

		// Can process Y and Z in one loop, we only care about horizontal scan lines.
		for (size_t yz = 0; yz < ydim * zdim; yz++)
		{
			// NOTE: no x++ (intentional!) - see inside of loop
			for (size_t x = 0; x < xdim; )
			{
				size_t pixel_index = yz * xdim + x;

				// Unsigned subtract makes one-sided test work:
				if ((x - radius_u) < safe_width)
				{
					// We're in the safe region, can process ASTCENC_SIMD_WIDTH pixels
					// at once and don't need to worry about re-normalizing
					vfloat sum = vfloat::zero();
					for (int k = -radius; k <= radius; k++)
					{
						vfloat weight(coeffs[k]);
						sum += weight * vfloat(input + (pixel_index + k));
					}

					store(sum, output + pixel_index);
					x += ASTCENC_SIMD_WIDTH;
				}
				else
				{
					// Near boundary, work one pixel at a time.
					float sum = 0.0f;
					float weight_sum = 0.0f;

					for (int k = -radius; k <= radius; k++)
					{
						// Compute displaced position. Negatives wrap around,
						// which shows up as >=width below.
						size_t sx = x + k;

						if (sx < xdim)
						{
							float weight = coeffs[k];
							weight_sum += weight;
							sum += weight * input[pixel_index + k];
						}
					}

					// Normalize by the actual sum of kernel weights used
					sum /= weight_sum;
					output[pixel_index] = sum;

					x++;
				}
			}
		}
	}
	else
	{
		// Direction is Y or Z
		intptr_t stride = (direction == 1) ? xdim : xdim * ydim;
		size_t max_in_dir = (direction == 1) ? (ydim - 1) : (zdim - 1);

		// Clear input padding
		store(vfloat::zero(), input + xdim * ydim * zdim);

		for (size_t z = 0; z < zdim; z++)
		{
			for (size_t y = 0; y < ydim; y++)
			{
				// We can do the kernel clamping and sum calc once for the whole scan line
				size_t pos_in_dir = (direction == 1) ? y : z;

				// Clamp kernel to active bounds
				int k0 = -static_cast<int>(astc::min(pos_in_dir, radius_u));
				int k1 = static_cast<int>(astc::min(max_in_dir - pos_in_dir, radius_u));

				// If we cut off pixels, determine normalization factor
				vfloat overall_scale(1.0f);
				if (k0 != -radius || k1 != radius)
				{
					float kernel_sum = 0.0f;
					for (int k = k0; k <= k1; k++)
					{
						kernel_sum += coeffs[k];
					}
					overall_scale = vfloat(1.0f / kernel_sum);
				}

				// Process the pixels
				// this runs off the ends of scan lines but input != output and both
				// have padding (that we just initialized to zero), so this is OK.
				size_t row_index = (z * ydim + y) * xdim;
				for (size_t x = 0; x < xdim; x += ASTCENC_SIMD_WIDTH)
				{
					vfloat sum = vfloat::zero();
					size_t pixel_index = row_index + x;

					for (int k = k0; k <= k1; k++)
					{
						vfloat weight(coeffs[k]);
						sum = sum + weight * vfloat(input + pixel_index + k * stride);
					}

					store(overall_scale * sum, output + pixel_index);
				}
			}
		}
	}
}

/**
 * @brief Separable convolution filter for 3D arrays.
 *
 * All arrays are xdim * ydim * zdim values and are expected to have
 * at least ASTCENC_SIMD_WIDTH floats of padding at the end.
 *
 * @param input			   The input values.
 * @param[out] output	   Where the output values are written.
 * @param[out] workspace   Workspace memory.
 * @param xdim             Size of the array in the X dimension.
 * @param ydim             Size of the array in the Y dimension.
 * @param zdim             Size of the array in the Z dimension.
 * @param filter           The filter kernel to use.
 */
static void separable_convolve_3d(
	float* input,
	float* output,
	float* workspace,
	size_t xdim,
	size_t ydim,
	size_t zdim,
	const filter_kernel& filter
) {
	// Actual 3D blur?
	if (zdim > 1u)
	{
		// 3D: X, Y, then Z
		apply_1d_convolution_3d(input, output, xdim, ydim, zdim, filter, 0);
		apply_1d_convolution_3d(output, workspace, xdim, ydim, zdim, filter, 1);
		apply_1d_convolution_3d(workspace, output, xdim, ydim, zdim, filter, 2);
	}
	else
	{
		// 2D only: X then Y
		apply_1d_convolution_3d(input, workspace, xdim, ydim, zdim, filter, 0);
		apply_1d_convolution_3d(workspace, output, xdim, ydim, zdim, filter, 1);
	}
}

/** @brief Maximum number of texels in a tile. */
static constexpr unsigned int WEIGHTS_TILE_MAX_TEXELS { 256 * 256 };

/** @brief Actual size of a tile buffer, intentionally larger to allow for padding and avoid pow2 aliasing conflicts. */
static constexpr unsigned int WEIGHTS_TILE_BUFFER_TEXELS { WEIGHTS_TILE_MAX_TEXELS + 128 + ASTCENC_SIMD_WIDTH };

/** @brief Tile buffers for per-texel weight calculation. */
struct weights_tile_buffers
{
	/** @brief Input pixels organized by color channel */
	ASTCENC_ALIGNAS float chan[4][WEIGHTS_TILE_BUFFER_TEXELS];
	/** @brief Filter work space */
	ASTCENC_ALIGNAS float work[WEIGHTS_TILE_BUFFER_TEXELS];
	/** @brief Filter result space */
	ASTCENC_ALIGNAS float filtered[WEIGHTS_TILE_BUFFER_TEXELS];
};

/**
 * @brief LZ-RDO setup pass does initial computations to determine per-texel weights.
 *
 * This pass allocates memory for global working buffers, then chops the image into
 * independent fixed-size tiles (cuboids for 3D textures) and computes a per-texel
 * weight that weights down regions with a lot of high-frequency energy.
 *
 * Quantization errors are much more visible and objectionable in smooth regions
 * (flat areas or slow gradients). Where the image is noisy or has sharp edges,
 * small differences in pixel values are less apparent.
 *
 * @param ctxo    The ASTC encoder context.
 * @param img     The image being encoded.
 * @param swz     The component swizzle for the image being encoded.
 */
void lz_rdo_setup(
	astcenc_context& ctxo,
	const astcenc_image& img,
	const astcenc_swizzle& swz
) {
	const astcenc_config& config = ctxo.context.config;
	const block_size_descriptor* bsd = ctxo.context.bsd;

	// Choose tile size to be an integer multiple of block size
	size_t tile_texels_x, tile_texels_y, tile_texels_z;
	if (bsd->dim_z == 1)
	{
		// 2D blocks - use approximately 256x256 tiles
		tile_texels_x = (256 / bsd->dim_x) * bsd->dim_x;
		tile_texels_y = (256 / bsd->dim_y) * bsd->dim_y;
		tile_texels_z = 1;
	}
	else
	{
		// 3D blocks - use approximately 32x32x32 tiles
		tile_texels_x = (32 / bsd->dim_x) * bsd->dim_x;
		tile_texels_y = (32 / bsd->dim_y) * bsd->dim_y;
		tile_texels_z = (32 / bsd->dim_z) * bsd->dim_z;
	}

	assert(tile_texels_x * tile_texels_y * tile_texels_z <= WEIGHTS_TILE_MAX_TEXELS);

	// Calculate number of tiles needed in each dimension
	size_t tiles_x = ceil_div(img.dim_x, tile_texels_x);
	size_t tiles_y = ceil_div(img.dim_y, tile_texels_y);
	size_t tiles_z = ceil_div(img.dim_z, tile_texels_z);
	size_t total_tile_count = tiles_x * tiles_y * tiles_z;

	size_t num_blocks =
		ceil_div(img.dim_x, bsd->dim_x) *
		ceil_div(img.dim_y, bsd->dim_y) *
		ceil_div(img.dim_z, bsd->dim_z);

	// Only the first thread actually runs the initializer
	auto init_func = [&ctxo, total_tile_count, num_blocks]()
	{
		ctxo.context.lz_rdo_restricted_blocks = new uint8_t[num_blocks * ASTC_BLOCK_BYTES];
		ctxo.context.lz_rdo_per_texel_weights = new float[num_blocks * ctxo.context.bsd->texel_count + ASTCENC_SIMD_WIDTH]; // Padded
		return total_tile_count;
	};

	ctxo.manage_lz_rdo_setup.init(init_func);

	float* block_texel_weights = ctxo.context.lz_rdo_per_texel_weights;

	// Image dimensions in blocks (only need X and Y)
	size_t img_blocks_x = ceil_div(img.dim_x, bsd->dim_x);
	size_t img_blocks_y = ceil_div(img.dim_y, bsd->dim_y);

	// Thread-local buffers for tile data. Allocated on demand.
	weights_tile_buffers* buffers = nullptr;

	// Set up image block
	image_block blk;
	blk.texel_count = bsd->texel_count;
	blk.decode_unorm8 = config.flags & ASTCENC_FLG_USE_DECODE_UNORM8;

	auto load_func = load_image_block;

	while (true)
	{
		// Grab next tile
		size_t tile_count;
		size_t current_tile_index = ctxo.manage_lz_rdo_setup.get_task_assignment(1, tile_count);
		if (!tile_count)
		{
			break;
		}

		// First time this thread gets a tile, allocate the workspace buffers.
		// Do this lazily in case we have more workers than tiles to work on
		// (common for small textures)
		if (!buffers)
		{
			buffers = aligned_malloc<weights_tile_buffers>(sizeof(weights_tile_buffers), ASTCENC_VECALIGN);
		}

		// Determine tile coordinates
		size_t tx, ty, tz;
		decode_index_to_coords(current_tile_index, tiles_x, tiles_y, tx, ty, tz);

		// Calculate tile bounds
		size_t tile_xpos = tx * tile_texels_x;
		size_t tile_ypos = ty * tile_texels_y;
		size_t tile_zpos = tz * tile_texels_z;

		// Number of active texels in a tile ...
		size_t tile_active_texels_x = astc::min(tile_texels_x, img.dim_x - tile_xpos);
		size_t tile_active_texels_y = astc::min(tile_texels_y, img.dim_y - tile_ypos);
		size_t tile_active_texels_z = astc::min(tile_texels_z, img.dim_z - tile_zpos);

		// ... rounded up to an integer block multiple
		size_t tile_xdim = ceil_div(tile_active_texels_x, bsd->dim_x) * bsd->dim_x;
		size_t tile_ydim = ceil_div(tile_active_texels_y, bsd->dim_y) * bsd->dim_y;
		size_t tile_zdim = ceil_div(tile_active_texels_z, bsd->dim_z) * bsd->dim_z;

		// Load input pixels
		for (size_t bz = 0; bz < tile_zdim; bz += bsd->dim_z)
		{
			for (size_t by = 0; by < tile_ydim; by += bsd->dim_y)
			{
				for (size_t bx = 0; bx < tile_xdim; bx += bsd->dim_x)
				{
					// Load the block
					load_func(
						config.profile,
						img, blk, *bsd,
						tile_xpos + bx,
						tile_ypos + by,
						tile_zpos + bz,
						swz
					);

					// Copy complete blocks; the loader extrudes texels beyond the image edge
					size_t src_index = 0;

					for (size_t z = 0; z < bsd->dim_z; z++)
					{
						for (size_t y = 0; y < bsd->dim_y; y++)
						{
							size_t dst_index = ((bz + z) * tile_ydim + (by + y)) * tile_xdim + bx;
							std::memcpy(buffers->chan[0] + dst_index, blk.data_r + src_index, bsd->dim_x * sizeof(float));
							std::memcpy(buffers->chan[1] + dst_index, blk.data_g + src_index, bsd->dim_x * sizeof(float));
							std::memcpy(buffers->chan[2] + dst_index, blk.data_b + src_index, bsd->dim_x * sizeof(float));
							std::memcpy(buffers->chan[3] + dst_index, blk.data_a + src_index, bsd->dim_x * sizeof(float));

							src_index += bsd->dim_x;
						}
					}
				}
			}
		}

		size_t tile_texels = tile_xdim * tile_ydim * tile_zdim;

		// Compute the low-pass filtered input and, from there, the high-pass energy
		for (unsigned int c = 0; c < 4; c++)
		{
			float* channel = buffers->chan[c];
			float* filtered = buffers->filtered;

			separable_convolve_3d(
				channel,
				filtered,
				buffers->work,
				tile_xdim, tile_ydim, tile_zdim,
				kernel_initial_lpf
			);

			// Clear padding
			store(vfloat::zero(), channel + tile_texels);
			store(vfloat::zero(), filtered + tile_texels);

			// Compute highpass energy (squared difference between original and low-pass filtered pixels)
			// store it in the space for channel 0 which is freed up as we're going along.
			float* hp_energy = buffers->chan[0];

			for (size_t i = 0; i < tile_texels; i += ASTCENC_SIMD_WIDTH)
			{
				vfloat cur_energy = vfloat::zero();
				if (c != 0)
				{
					cur_energy = loada(hp_energy + i);
				}

				vfloat diff = loada(channel + i) - loada(filtered + i);
				cur_energy = cur_energy + diff * diff;
				storea(cur_energy, hp_energy + i);
			}
		}

		// Apply second blur to spread the high-pass energy around a bit
		float* tile_output = buffers->filtered;

		separable_convolve_3d(
			buffers->chan[0],
			tile_output,
			buffers->work,
			tile_xdim, tile_ydim, tile_zdim,
			kernel_spread
		);

		// Clear padding
		store(vfloat::zero(), tile_output + tile_texels);

		// Transform into final activity weights
		vfloat C1 { 1.0f / 257.0f };
		vfloat C2 { 257.0f / 256.0f };
		vfloat activity_scale { 4.0f / 255.0f };

		for (size_t i = 0; i < tile_texels; i += ASTCENC_SIMD_WIDTH)
		{
			vfloat values(tile_output + i);
			vfloat weights = C1 / (C2 + activity_scale * sqrt(values));
			store(weights, tile_output + i);
		}

		// Copy to output in block layout
		size_t tile_bx = tile_xpos / bsd->dim_x;
		size_t tile_by = tile_ypos / bsd->dim_y;
		size_t tile_bz = tile_zpos / bsd->dim_z;
		size_t tile_blocks_x = tile_xdim / bsd->dim_x;
		size_t tile_blocks_y = tile_ydim / bsd->dim_y;
		size_t tile_blocks_z = tile_zdim / bsd->dim_z;

		for (size_t bz = 0; bz < tile_blocks_z; bz++)
		{
			for (size_t by = 0; by < tile_blocks_y; by++)
			{
				for (size_t bx = 0; bx < tile_blocks_x; bx++)
				{
					size_t block_index = ((bz + tile_bz) * img_blocks_y + (by + tile_by)) * img_blocks_x + (bx + tile_bx);
					float* block_weights = block_texel_weights + block_index * bsd->texel_count;

					// Iterate over texels in the block
					size_t block_weight_index = 0;

					for (size_t z = 0; z < bsd->dim_z; z++)
					{
						size_t zt = bz * bsd->dim_z + z;
						for (size_t y = 0; y < bsd->dim_y; y++)
						{
							size_t yt = by * bsd->dim_y + y;
							size_t src_row_idx = (zt * tile_ydim + yt) * tile_xdim;

							for (size_t x = 0; x < bsd->dim_x; x++)
							{
								size_t xt = bx * bsd->dim_x + x;
								if (((tile_xpos + xt) < img.dim_x) &&
								    ((tile_ypos + yt) < img.dim_y) &&
								    ((tile_zpos + zt) < img.dim_z))
								{
									// Copy over weight
									block_weights[block_weight_index] = tile_output[src_row_idx + xt];
								}
								else
								{
									// Texels outside image bounds get zero weight
									block_weights[block_weight_index] = 0.0f;
								}

								block_weight_index++;
							}
						}
					}
				}
			}
		}

		ctxo.manage_lz_rdo_setup.complete_task_assignment(tile_count);
	}

	// Clean up
	if (buffers)
	{
		aligned_free<weights_tile_buffers>(buffers);
	}
}

/* See header for documentation. */
void lz_rdo_optimize(
	astcenc_context& ctxo,
	size_t thread_index,
	astcenc_image& img,
	const astcenc_swizzle& swz,
	float init_pct,
	uint8_t* data
) {
	astcenc_contexti& ctx = ctxo.context;

	float lambda = ctx.config.lz_rdo_lambda;
	assert(lambda > 0.0f); // Caller must ensure this.

	// Global lambda scaling fudge factor
	lambda *= 0.75f;

	// Squared error in blocks scales with texel count; this makes lambda mostly
	// independent of block size.
	lambda *= static_cast<float>(ctx.bsd->texel_count) / 16.0f;

	// Run the optimization passes!
	// calls ctxo.manage_lz_rdo_optimize.wait()
	// so all work completed on return
	lz_rdo_optimization_passes(
		ctxo, thread_index, data, lambda, init_pct,
		img, swz
	);

	auto term_optimize = [&ctx]()
	{
		delete[] ctx.lz_rdo_restricted_blocks;
		ctx.lz_rdo_restricted_blocks = nullptr;

		delete[] ctx.lz_rdo_per_texel_weights;
		ctx.lz_rdo_per_texel_weights = nullptr;
	};

	ctxo.manage_lz_rdo_optimize.term(term_optimize);
}

#endif

