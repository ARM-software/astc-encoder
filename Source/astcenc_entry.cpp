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
 * @brief Functions for the library entrypoint.
 */

#include <cstring>
#include <new>

#include "astcenc.h"
#include "astcenc_internal.h"

// The ASTC codec is written with the assumption that a float threaded through
// the "if32" union will in fact be stored and reloaded as a 32-bit IEEE-754 single-precision
// float, stored with round-to-nearest rounding. This is always the case in an
// IEEE-754 compliant system, however not every system is actually IEEE-754 compliant
// in the first place. As such, we run a quick test to check that this is actually the case
// (e.g. gcc on 32-bit x86 will typically fail unless -msse2 -mfpmath=sse2 is specified).
static astcenc_error validate_cpu_float()
{
	if32 p;
	volatile float xprec_testval = 2.51f;
	p.f = xprec_testval + 12582912.0f;
	float q = p.f - 12582912.0f;

	if (q != 3.0f)
	{
		return ASTCENC_ERR_BAD_CPU_FLOAT;
	}

	return ASTCENC_SUCCESS;
}

static astcenc_error validate_cpu_isa()
{
	#if ASTCENC_SSE >= 42
		if (!cpu_supports_sse42())
		{
			return ASTCENC_ERR_BAD_CPU_ISA;
		}
	#endif

	#if ASTCENC_POPCNT >= 1
		if (!cpu_supports_popcnt())
		{
			return ASTCENC_ERR_BAD_CPU_ISA;
		}
	#endif

	#if ASTCENC_AVX >= 2
		if (!cpu_supports_avx2())
		{
			return ASTCENC_ERR_BAD_CPU_ISA;
		}
	#endif

	return ASTCENC_SUCCESS;
}

static astcenc_error validate_profile(
	astcenc_profile profile
) {
	switch(profile)
	{
	case ASTCENC_PRF_LDR_SRGB:
	case ASTCENC_PRF_LDR:
	case ASTCENC_PRF_HDR_RGB_LDR_A:
	case ASTCENC_PRF_HDR:
		return ASTCENC_SUCCESS;
	default:
		return ASTCENC_ERR_BAD_PROFILE;
	}
}

static astcenc_error validate_block_size(
	unsigned int block_x,
	unsigned int block_y,
	unsigned int block_z
) {
	if (((block_z <= 1) && is_legal_2d_block_size(block_x, block_y)) ||
	    ((block_z >= 2) && is_legal_3d_block_size(block_x, block_y, block_z)))
	{
		return ASTCENC_SUCCESS;
	}

	return ASTCENC_ERR_BAD_BLOCK_SIZE;
}

static astcenc_error validate_flags(
	unsigned int flags
) {
	// Flags field must not contain any unknown flag bits
	unsigned int exMask = ~ASTCENC_ALL_FLAGS;
	if (astc::popcount(flags & exMask) != 0)
	{
		return ASTCENC_ERR_BAD_FLAGS;
	}

	// Flags field must only contain at most a single map type
	exMask = ASTCENC_FLG_MAP_MASK | ASTCENC_FLG_MAP_NORMAL;
	if (astc::popcount(flags & exMask) > 1)
	{
		return ASTCENC_ERR_BAD_FLAGS;
	}

	return ASTCENC_SUCCESS;
}

#if !defined(ASTCENC_DECOMPRESS_ONLY)
static astcenc_error validate_compression_swz(
	astcenc_swz swizzle
) {
	switch(swizzle)
	{
	case ASTCENC_SWZ_R:
	case ASTCENC_SWZ_G:
	case ASTCENC_SWZ_B:
	case ASTCENC_SWZ_A:
	case ASTCENC_SWZ_0:
	case ASTCENC_SWZ_1:
		return ASTCENC_SUCCESS;
	default:
		return ASTCENC_ERR_BAD_SWIZZLE;
	}
}

static astcenc_error validate_compression_swizzle(
	astcenc_swizzle swizzle
) {
	if (validate_compression_swz(swizzle.r) ||
	    validate_compression_swz(swizzle.g) ||
	    validate_compression_swz(swizzle.b) ||
	    validate_compression_swz(swizzle.a))
	{
		return ASTCENC_ERR_BAD_SWIZZLE;
	}

	return ASTCENC_SUCCESS;
}
#endif

static astcenc_error validate_decompression_swz(
	astcenc_swz swizzle
) {
	switch(swizzle)
	{
	case ASTCENC_SWZ_R:
	case ASTCENC_SWZ_G:
	case ASTCENC_SWZ_B:
	case ASTCENC_SWZ_A:
	case ASTCENC_SWZ_0:
	case ASTCENC_SWZ_1:
	case ASTCENC_SWZ_Z:
		return ASTCENC_SUCCESS;
	default:
		return ASTCENC_ERR_BAD_SWIZZLE;
	}
}

static astcenc_error validate_decompression_swizzle(
	astcenc_swizzle swizzle
) {
	if (validate_decompression_swz(swizzle.r) ||
	    validate_decompression_swz(swizzle.g) ||
	    validate_decompression_swz(swizzle.b) ||
	    validate_decompression_swz(swizzle.a))
	{
		return ASTCENC_ERR_BAD_SWIZZLE;
	}

	return ASTCENC_SUCCESS;
}

/**
 * Validate that an incoming configuration is in-spec.
 *
 * This function can respond in two ways:
 *
 *   * Numerical inputs that have valid ranges are clamped to those valid
 *     ranges. No error is thrown for out-of-range inputs in this case.
 *   * Numerical inputs and logic inputs are are logically invalid and which
 *     make no sense algorithmically will return an error.
 */
static astcenc_error validate_config(
	astcenc_config &config,
	unsigned int thread_count
) {
	astcenc_error status;

	status = validate_profile(config.profile);
	if (status != ASTCENC_SUCCESS)
	{
		return status;
	}

	status = validate_flags(config.flags);
	if (status != ASTCENC_SUCCESS)
	{
		return status;
	}

	status = validate_block_size(config.block_x, config.block_y, config.block_z);
	if (status != ASTCENC_SUCCESS)
	{
		return status;
	}

	// Decompress-only contexts must be single threaded
	if ((config.flags & ASTCENC_FLG_DECOMPRESS_ONLY) && (thread_count > 1))
	{
		return ASTCENC_ERR_BAD_PARAM;
	}

#if defined(ASTCENC_DECOMPRESS_ONLY)
	// Decompress-only builds only support decompress-only contexts
	if (!(config.flags & ASTCENC_FLG_DECOMPRESS_ONLY))
	{
		return ASTCENC_ERR_BAD_PARAM;
	}
#endif

	config.v_rgba_mean_stdev_mix = MAX(config.v_rgba_mean_stdev_mix, 0.0f);
	config.v_rgb_power = MAX(config.v_rgb_power, 0.0f);
	config.v_rgb_base = MAX(config.v_rgb_base, 0.0f);
	config.v_rgb_mean = MAX(config.v_rgb_mean, 0.0f);
	config.v_rgb_stdev = MAX(config.v_rgb_stdev, 0.0f);
	config.v_a_power = MAX(config.v_a_power, 0.0f);
	config.v_a_base = MAX(config.v_a_base, 0.0f);
	config.v_a_mean = MAX(config.v_a_mean, 0.0f);
	config.v_a_stdev = MAX(config.v_a_stdev, 0.0f);

	config.b_deblock_weight = MAX(config.b_deblock_weight, 0.0f);

	config.tune_partition_limit = astc::clampi(config.tune_partition_limit, 1, PARTITION_COUNT);
	config.tune_block_mode_limit = astc::clampi(config.tune_block_mode_limit, 1, 100);
	config.tune_refinement_limit = MAX(config.tune_refinement_limit, 1);
	config.tune_db_limit = MAX(config.tune_db_limit, 0.0f);
	config.tune_partition_early_out_limit = MAX(config.tune_partition_early_out_limit, 0.0f);
	config.tune_two_plane_early_out_limit = MAX(config.tune_two_plane_early_out_limit, 0.0f);

	// Specifying a zero weight color component is not allowed; force to small value
	float max_weight = MAX(MAX(config.cw_r_weight, config.cw_g_weight),
	                       MAX(config.cw_b_weight, config.cw_a_weight));
	if (max_weight > 0.0f)
	{
		max_weight /= 1000.0f;
		config.cw_r_weight = MAX(config.cw_r_weight, max_weight);
		config.cw_g_weight = MAX(config.cw_g_weight, max_weight);
		config.cw_b_weight = MAX(config.cw_b_weight, max_weight);
		config.cw_a_weight = MAX(config.cw_a_weight, max_weight);
	}
	// If all color components error weights are zero then return an error
	else
	{
		return ASTCENC_ERR_BAD_PARAM;
	}

	return ASTCENC_SUCCESS;
}

astcenc_error astcenc_config_init(
	astcenc_profile profile,
	unsigned int block_x,
	unsigned int block_y,
	unsigned int block_z,
	astcenc_preset preset,
	unsigned int flags,
	astcenc_config& config
) {
	astcenc_error status;

	// Zero init all config fields; although most of will be over written
	std::memset(&config, 0, sizeof(config));

	// Process the block size
	block_z = MAX(block_z, 1); // For 2D blocks Z==0 is accepted, but convert to 1
	status = validate_block_size(block_x, block_y, block_z);
	if (status != ASTCENC_SUCCESS)
	{
		return status;
	}

	config.block_x = block_x;
	config.block_y = block_y;
	config.block_z = block_z;

	float texels = static_cast<float>(block_x * block_y * block_z);
	float ltexels = logf(texels) / logf(10.0f);

	// Process the performance preset; note that this must be done before we
	// process any additional settings, such as color profile and flags, which
	// may replace some of these settings with more use case tuned values
	switch(preset)
	{
	case ASTCENC_PRE_FAST:
		config.tune_partition_limit = 4;
		config.tune_block_mode_limit = 50;
		config.tune_refinement_limit = 1;
		config.tune_db_limit = MAX(85 - 35 * ltexels, 63 - 19 * ltexels);
		config.tune_partition_early_out_limit = 1.0f;
		config.tune_two_plane_early_out_limit = 0.5f;
		break;
	case ASTCENC_PRE_MEDIUM:
		config.tune_partition_limit = 25;
		config.tune_block_mode_limit = 75;
		config.tune_refinement_limit = 2;
		config.tune_db_limit = MAX(95 - 35 * ltexels, 70 - 19 * ltexels);
		config.tune_partition_early_out_limit = 1.2f;
		config.tune_two_plane_early_out_limit = 0.75f;
		break;
	case ASTCENC_PRE_THOROUGH:
		config.tune_partition_limit = 100;
		config.tune_block_mode_limit = 95;
		config.tune_refinement_limit = 4;
		config.tune_db_limit = MAX(105 - 35 * ltexels, 77 - 19 * ltexels);
		config.tune_partition_early_out_limit = 2.5f;
		config.tune_two_plane_early_out_limit = 0.95f;
		break;
	case ASTCENC_PRE_EXHAUSTIVE:
		config.tune_partition_limit = 1024;
		config.tune_block_mode_limit = 100;
		config.tune_refinement_limit = 4;
		config.tune_db_limit = 999.0f;
		config.tune_partition_early_out_limit = 1000.0f;
		config.tune_two_plane_early_out_limit = 0.99f;
		break;
	default:
		return ASTCENC_ERR_BAD_PRESET;
	}

	// Set heuristics to the defaults for each color profile
	config.v_rgba_radius = 0;
	config.v_rgba_mean_stdev_mix = 0.0f;

	config.cw_r_weight = 1.0f;
	config.cw_g_weight = 1.0f;
	config.cw_b_weight = 1.0f;
	config.cw_a_weight = 1.0f;

	config.a_scale_radius = 0;

	config.b_deblock_weight = 0.0f;

	config.profile = profile;
	switch(profile)
	{
	case ASTCENC_PRF_LDR:
	case ASTCENC_PRF_LDR_SRGB:
		config.v_rgb_power = 1.0f;
		config.v_rgb_base = 1.0f;
		config.v_rgb_mean = 0.0f;
		config.v_rgb_stdev = 0.0f;

		config.v_a_power = 1.0f;
		config.v_a_base = 1.0f;
		config.v_a_mean = 0.0f;
		config.v_a_stdev = 0.0f;
		break;
	case ASTCENC_PRF_HDR_RGB_LDR_A:
		config.v_rgb_power = 0.75f;
		config.v_rgb_base = 0.0f;
		config.v_rgb_mean = 1.0f;
		config.v_rgb_stdev = 0.0f;

		config.v_a_power = 1.0f;
		config.v_a_base = 0.05f;
		config.v_a_mean = 0.0f;
		config.v_a_stdev = 0.0f;

		config.tune_db_limit = 999.0f;
		break;
	case ASTCENC_PRF_HDR:
		config.v_rgb_power = 0.75f;
		config.v_rgb_base = 0.0f;
		config.v_rgb_mean = 1.0f;
		config.v_rgb_stdev = 0.0f;

		config.v_a_power = 0.75f;
		config.v_a_base = 0.0f;
		config.v_a_mean = 1.0f;
		config.v_a_stdev = 0.0f;

		config.tune_db_limit = 999.0f;
		break;
	default:
		return ASTCENC_ERR_BAD_PROFILE;
	}

	// Flags field must not contain any unknown flag bits
	status = validate_flags(flags);
	if (status != ASTCENC_SUCCESS)
	{
		return status;
	}

	if (flags & ASTCENC_FLG_MAP_NORMAL)
	{
		config.cw_r_weight = 1.0f;
		config.cw_g_weight = 0.0f;
		config.cw_b_weight = 0.0f;
		config.cw_a_weight = 1.0f;
		config.tune_partition_early_out_limit = 1000.0f;
		config.tune_two_plane_early_out_limit = 0.99f;

 		if (flags & ASTCENC_FLG_USE_PERCEPTUAL)
		{
			config.b_deblock_weight = 1.8f;
			config.v_rgba_radius = 3;
			config.v_rgba_mean_stdev_mix = 0.0f;
			config.v_rgb_mean = 0.0f;
			config.v_rgb_stdev = 50.0f;
			config.v_a_mean = 0.0f;
			config.v_a_stdev = 50.0f;
		}
	}

	if (flags & ASTCENC_FLG_MAP_MASK)
	{
		config.v_rgba_radius = 3;
		config.v_rgba_mean_stdev_mix = 0.03f;
		config.v_rgb_mean = 0.0f;
		config.v_rgb_stdev = 25.0f;
		config.v_a_mean = 0.0f;
		config.v_a_stdev = 25.0f;
	}

	config.flags = flags;

	return ASTCENC_SUCCESS;
}

astcenc_error astcenc_context_alloc(
	astcenc_config const& config,
	unsigned int thread_count,
	astcenc_context** context
) {
	astcenc_error status;
	astcenc_context* ctx = nullptr;
	block_size_descriptor* bsd = nullptr;

	status = validate_cpu_float();
	if (status != ASTCENC_SUCCESS)
	{
		return status;
	}

	status = validate_cpu_isa();
	if (status != ASTCENC_SUCCESS)
	{
		return status;
	}

	if (thread_count == 0)
	{
		return ASTCENC_ERR_BAD_PARAM;
	}

	ctx = new astcenc_context;
	ctx->thread_count = thread_count;
	ctx->config = config;
	ctx->working_buffers = nullptr;

	// These are allocated per-compress, as they depend on image size
	ctx->input_averages = nullptr;
	ctx->input_variances = nullptr;
	ctx->input_alpha_averages = nullptr;

	// Copy the config first and validate the copy (we may modify it)
	status = validate_config(ctx->config, thread_count);
	if (status != ASTCENC_SUCCESS)
	{
		delete ctx;
		return status;
	}

	bsd = new block_size_descriptor;
	init_block_size_descriptor(config.block_x, config.block_y, config.block_z, bsd);
	ctx->bsd = bsd;

#if !defined(ASTCENC_DECOMPRESS_ONLY)
	// Do setup only needed by compression
	if (!(status & ASTCENC_FLG_DECOMPRESS_ONLY))
	{
		// Expand deblock supression into a weight scale per texel in the block
		expand_deblock_weights(*ctx);

		// Turn a dB limit into a per-texel error for faster use later
		if ((ctx->config.profile == ASTCENC_PRF_LDR) || (ctx->config.profile == ASTCENC_PRF_LDR_SRGB))
		{
			ctx->config.tune_db_limit = powf(0.1f, ctx->config.tune_db_limit * 0.1f) * 65535.0f * 65535.0f;
		}
		else
		{
			ctx->config.tune_db_limit = 0.0f;
		}

		size_t worksize = sizeof(compress_symbolic_block_buffers) * thread_count;
		ctx->working_buffers = aligned_malloc<compress_symbolic_block_buffers>(worksize , 32);
		if (!ctx->working_buffers)
		{
			goto error_oom;
		}
	}
#endif

	*context = ctx;

	// TODO: Currently static memory; should move to context memory
#if !defined(ASTCENC_DECOMPRESS_ONLY)
	prepare_angular_tables();
#endif
	build_quantization_mode_table();

	return ASTCENC_SUCCESS;

error_oom:
	term_block_size_descriptor(bsd);
	delete bsd;
	delete ctx;
	*context = nullptr;
	return ASTCENC_ERR_OUT_OF_MEM;
}

void astcenc_context_free(
	astcenc_context* ctx
) {
	if (ctx)
	{
		aligned_free<compress_symbolic_block_buffers>(ctx->working_buffers);
		term_block_size_descriptor(ctx->bsd);
		delete ctx->bsd;
		delete ctx;
	}
}

#if !defined(ASTCENC_DECOMPRESS_ONLY)
static void compress_image(
	astcenc_context& ctx,
	unsigned int thread_index,
	const astcenc_image& image,
	astcenc_swizzle swizzle,
	uint8_t* buffer
) {
	const block_size_descriptor *bsd = ctx.bsd;
	int block_x = bsd->xdim;
	int block_y = bsd->ydim;
	int block_z = bsd->zdim;
	astcenc_profile decode_mode = ctx.config.profile;

	imageblock pb;
	int dim_x = image.dim_x;
	int dim_y = image.dim_y;
	int dim_z = image.dim_z;
	int xblocks = (dim_x + block_x - 1) / block_x;
	int yblocks = (dim_y + block_y - 1) / block_y;
	int zblocks = (dim_z + block_z - 1) / block_z;

	int row_blocks = xblocks;
	int plane_blocks = xblocks * yblocks;

	// Use preallocated scratch buffer
	auto temp_buffers = &(ctx.working_buffers[thread_index]);

	// Only the first thread actually runs the initializer
	ctx.manage_compress.init(zblocks * yblocks * xblocks);

	// All threads run this processing loop until there is no work remaining
	while (true)
	{
		unsigned int count;
		unsigned int base = ctx.manage_compress.get_task_assignment(4, count);
		if (!count)
		{
			break;
		}

		for (unsigned int i = base; i < base + count; i++)
		{
			// Decode i into x, y, z block indices
			int z = i / plane_blocks;
			unsigned int rem = i - (z * plane_blocks);
			int y = rem / row_blocks;
			int x = rem - (y * row_blocks);

			// Decompress
			int offset = ((z * yblocks + y) * xblocks + x) * 16;
			const uint8_t *bp = buffer + offset;
			fetch_imageblock(decode_mode, image, &pb, bsd, x * block_x, y * block_y, z * block_z, swizzle);
			symbolic_compressed_block scb;
			compress_symbolic_block(ctx, image, decode_mode, bsd, &pb, &scb, temp_buffers);
			*(physical_compressed_block*) bp = symbolic_to_physical(bsd, &scb);
		}

		ctx.manage_compress.complete_task_assignment(count);
	};
}
#endif

astcenc_error astcenc_compress_image(
	astcenc_context* ctx,
	astcenc_image& image,
	astcenc_swizzle swizzle,
	uint8_t* data_out,
	size_t data_len,
	unsigned int thread_index
) {
#if defined(ASTCENC_DECOMPRESS_ONLY)
	(void)ctx;
	(void)image;
	(void)swizzle;
	(void)data_out;
	(void)data_len;
	(void)thread_index;
	return ASTCENC_ERR_BAD_CONTEXT;
#else
	astcenc_error status;

	if (ctx->config.flags & ASTCENC_FLG_DECOMPRESS_ONLY)
	{
		return ASTCENC_ERR_BAD_CONTEXT;
	}

	status = validate_compression_swizzle(swizzle);
	if (status != ASTCENC_SUCCESS)
	{
		return status;
	}

	if (thread_index >= ctx->thread_count)
	{
		return ASTCENC_ERR_BAD_PARAM;
	}

	unsigned int block_x = ctx->config.block_x;
	unsigned int block_y = ctx->config.block_y;
	unsigned int block_z = ctx->config.block_z;

	unsigned int xblocks = (image.dim_x + block_x - 1) / block_x;
	unsigned int yblocks = (image.dim_y + block_y - 1) / block_y;
	unsigned int zblocks = (image.dim_z + block_z - 1) / block_z;

	// Check we have enough output space (16 bytes per block)
	size_t size_needed = xblocks * yblocks * zblocks * 16;
	if (data_len < size_needed)
	{
		return ASTCENC_ERR_OUT_OF_MEM;
	}

	if (image.dim_pad > 0 ||
	    ctx->config.v_rgb_mean != 0.0f || ctx->config.v_rgb_stdev != 0.0f ||
	    ctx->config.v_a_mean != 0.0f || ctx->config.v_a_stdev != 0.0f)
	{
		// First thread to enter will do setup, other threads will subsequently
		// enter the critical section but simply skip over the initialization
		auto init_avg_var = [ctx, &image, swizzle]() {
			// Perform memory allocations for the destination buffers
			size_t texel_count = image.dim_x * image.dim_y * image.dim_z;
			ctx->input_averages = new float4[texel_count];
			ctx->input_variances = new float4[texel_count];
			ctx->input_alpha_averages = new float[texel_count];

			return init_compute_averages_and_variances(
				image, ctx->config.v_rgb_power, ctx->config.v_a_power,
				ctx->config.v_rgba_radius, ctx->config.a_scale_radius, swizzle,
				ctx->arg, ctx->ag);
		};

		// Only the first thread actually runs the initializer
		ctx->manage_avg_var.init(init_avg_var);

		// All threads will enter this function and dynamically grab work
		compute_averages_and_variances(*ctx, ctx->ag);
	}

	// Wait for compute_averages_and_variances to complete before compressing
	ctx->manage_avg_var.wait();

	compress_image(*ctx, thread_index, image, swizzle, data_out);

	// Wait for compress to complete before freeing memory
	ctx->manage_compress.wait();

	auto term_compress = [ctx]() {
		delete[] ctx->input_averages;
		ctx->input_averages = nullptr;

		delete[] ctx->input_variances;
		ctx->input_variances = nullptr;

		delete[] ctx->input_alpha_averages;
		ctx->input_alpha_averages = nullptr;
	};

	// Only the first thread to arrive actually runs the term
	ctx->manage_compress.term(term_compress);

	return ASTCENC_SUCCESS;
#endif
}

astcenc_error astcenc_compress_reset(
	astcenc_context* ctx
) {
#if defined(ASTCENC_DECOMPRESS_ONLY)
	(void)ctx;
	return ASTCENC_ERR_BAD_CONTEXT;
#else
	if (ctx->config.flags & ASTCENC_FLG_DECOMPRESS_ONLY)
	{
		return ASTCENC_ERR_BAD_CONTEXT;
	}

	ctx->manage_avg_var.reset();
	ctx->manage_compress.reset();
	return ASTCENC_SUCCESS;
#endif
}

astcenc_error astcenc_decompress_image(
	astcenc_context* context,
	const uint8_t* data,
	size_t data_len,
	astcenc_image& image_out,
	astcenc_swizzle swizzle
) {
	astcenc_error status;

	status = validate_decompression_swizzle(swizzle);
	if (status != ASTCENC_SUCCESS)
	{
		return status;
	}

	unsigned int block_x = context->config.block_x;
	unsigned int block_y = context->config.block_y;
	unsigned int block_z = context->config.block_z;

	unsigned int xblocks = (image_out.dim_x + block_x - 1) / block_x;
	unsigned int yblocks = (image_out.dim_y + block_y - 1) / block_y;
	unsigned int zblocks = (image_out.dim_z + block_z - 1) / block_z;

	// Check we have enough output space (16 bytes per block)
	size_t size_needed = xblocks * yblocks * zblocks * 16;
	if (data_len < size_needed)
	{
		return ASTCENC_ERR_OUT_OF_MEM;
	}

	imageblock pb;

	for (unsigned int z = 0; z < zblocks; z++)
	{
		for (unsigned int y = 0; y < yblocks; y++)
		{
			for (unsigned int x = 0; x < xblocks; x++)
			{
				unsigned int offset = (((z * yblocks + y) * xblocks) + x) * 16;
				const uint8_t* bp = data + offset;
				physical_compressed_block pcb = *(physical_compressed_block *) bp;
				symbolic_compressed_block scb;
				physical_to_symbolic(context->bsd, pcb, &scb);
				decompress_symbolic_block(context->config.profile, context->bsd,
				                          x * block_x, y * block_y, z * block_z,
				                          &scb, &pb);
				write_imageblock(image_out, &pb, context->bsd,
				                 x * block_x, y * block_y, z * block_z, swizzle);
			}
		}
	}

	return ASTCENC_SUCCESS;
}

const char* astcenc_get_error_string(
	astcenc_error status
) {
	switch(status)
	{
	case ASTCENC_ERR_OUT_OF_MEM:
		return "ASTCENC_ERR_OUT_OF_MEM";
	case ASTCENC_ERR_BAD_CPU_FLOAT:
		return "ASTCENC_ERR_BAD_CPU_FLOAT";
	case ASTCENC_ERR_BAD_CPU_ISA:
		return "ASTCENC_ERR_BAD_CPU_ISA";
	case ASTCENC_ERR_BAD_PARAM:
		return "ASTCENC_ERR_BAD_PARAM";
	case ASTCENC_ERR_BAD_BLOCK_SIZE:
		return "ASTCENC_ERR_BAD_BLOCK_SIZE";
	case ASTCENC_ERR_BAD_PROFILE:
		return "ASTCENC_ERR_BAD_PROFILE";
	case ASTCENC_ERR_BAD_PRESET:
		return "ASTCENC_ERR_BAD_PRESET";
	case ASTCENC_ERR_BAD_FLAGS:
		return "ASTCENC_ERR_BAD_FLAGS";
	case ASTCENC_ERR_BAD_SWIZZLE:
		return "ASTCENC_ERR_BAD_SWIZZLE";
	case ASTCENC_ERR_BAD_CONTEXT:
		return "ASTCENC_ERR_BAD_CONTEXT";
	case ASTCENC_ERR_NOT_IMPLEMENTED:
		return "ASTCENC_ERR_NOT_IMPLEMENTED";
	default:
		return nullptr;
	}
}
