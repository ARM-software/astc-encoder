// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2025 Arm Limited
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
 * @brief Unit tests for guided compression.
 *
 * These tests verify that compress_block_guided() can reproduce encoding
 * quality comparable to the normal compress_block() path when given the
 * same block mode, partition, and format class parameters extracted from
 * a normal encode's output.
 */

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "gtest/gtest.h"

#include "../astcenc.h"
#include "../astcenc_internal.h"
#include "../astcenc_internal_entry.h"

namespace astcenc
{

// ============================================================================
// Test fixture
// ============================================================================

class GuidedCompressionTest : public ::testing::Test
{
protected:
	static constexpr unsigned int BLOCK_X = 8;
	static constexpr unsigned int BLOCK_Y = 8;
	static constexpr unsigned int TEXEL_COUNT = BLOCK_X * BLOCK_Y;

	static const astcenc_swizzle swizzle;

	astcenc_context* ctx_thorough = nullptr;
	astcenc_context* ctx_fastest = nullptr;

	void SetUp() override
	{
		astcenc_error status;
		astcenc_config config;

		// Create thorough context
		status = astcenc_config_init(
			ASTCENC_PRF_LDR, BLOCK_X, BLOCK_Y, 1,
			ASTCENC_PRE_THOROUGH, 0, &config);
		ASSERT_EQ(status, ASTCENC_SUCCESS);

		status = astcenc_context_alloc(&config, 1, &ctx_thorough);
		ASSERT_EQ(status, ASTCENC_SUCCESS);

		// Create fastest context
		status = astcenc_config_init(
			ASTCENC_PRF_LDR, BLOCK_X, BLOCK_Y, 1,
			ASTCENC_PRE_FASTEST, 0, &config);
		ASSERT_EQ(status, ASTCENC_SUCCESS);

		status = astcenc_context_alloc(&config, 1, &ctx_fastest);
		ASSERT_EQ(status, ASTCENC_SUCCESS);
	}

	void TearDown() override
	{
		if (ctx_thorough)
		{
			astcenc_context_free(ctx_thorough);
		}
		if (ctx_fastest)
		{
			astcenc_context_free(ctx_fastest);
		}
	}

	/**
	 * @brief Fill an image_block from a pixel array using load_image_block.
	 *
	 * @param ctx      The context (provides profile and BSD).
	 * @param pixels   RGBA U8 pixel data, BLOCK_X * BLOCK_Y * 4 bytes.
	 * @param blk      The image block to populate.
	 */
	void fill_image_block(
		astcenc_context* ctx,
		const uint8_t* pixels,
		image_block& blk
	) {
		astcenc_contexti& ctxi = ctx->context;

		// Set up a fake astcenc_image pointing at our pixel data
		astcenc_image img;
		img.dim_x = BLOCK_X;
		img.dim_y = BLOCK_Y;
		img.dim_z = 1;
		img.data_type = ASTCENC_TYPE_U8;
		// The data pointer is an array of z-slice pointers
		void* slice = const_cast<uint8_t*>(pixels);
		img.data = &slice;

		blk.texel_count = static_cast<uint8_t>(BLOCK_X * BLOCK_Y);

		load_image_block(
			ctxi.config.profile,
			img, blk, *ctxi.bsd,
			0, 0, 0, swizzle);

		blk.channel_weight = vfloat4(
			ctxi.config.cw_r_weight,
			ctxi.config.cw_g_weight,
			ctxi.config.cw_b_weight,
			ctxi.config.cw_a_weight);

		blk.decode_unorm8 = (ctxi.config.flags & ASTCENC_FLG_USE_DECODE_UNORM8) != 0;
	}

	/**
	 * @brief Compress a block normally and return the SCB.
	 */
	void compress_and_decode(
		astcenc_context* ctx,
		image_block& blk,
		uint8_t pcb[16],
		symbolic_compressed_block& scb
	) {
		astcenc_contexti& ctxi = ctx->context;
		compression_working_buffers& tmpbuf = ctxi.working_buffers[0];

		compress_block(ctxi, blk, pcb, tmpbuf);
		physical_to_symbolic(*ctxi.bsd, pcb, scb);
	}

	/**
	 * @brief Extract guide parameters from an SCB.
	 */
	static void extract_guide_params(
		const symbolic_compressed_block& scb,
		unsigned int& block_mode,
		unsigned int& partition_count,
		unsigned int& partition_index,
		uint8_t format_classes[BLOCK_MAX_PARTITIONS]
	) {
		block_mode = scb.block_mode;
		partition_count = scb.partition_count;
		partition_index = scb.partition_index;

		for (unsigned int j = 0; j < BLOCK_MAX_PARTITIONS; j++)
		{
			format_classes[j] = scb.color_formats[j] >> 2;
		}
	}

	/**
	 * @brief Compress with guided path and decode the result.
	 */
	void guided_compress_and_decode(
		astcenc_context* ctx,
		image_block& blk,
		unsigned int guide_block_mode,
		unsigned int guide_partition_count,
		unsigned int guide_partition_index,
		const uint8_t guide_format_classes[BLOCK_MAX_PARTITIONS],
		uint8_t pcb[16],
		symbolic_compressed_block& scb
	) {
		astcenc_contexti& ctxi = ctx->context;
		compression_working_buffers& tmpbuf = ctxi.working_buffers[0];

		compress_block_guided(
			ctxi, blk,
			guide_block_mode, guide_partition_count,
			guide_partition_index, guide_format_classes,
			pcb, tmpbuf);

		physical_to_symbolic(*ctxi.bsd, pcb, scb);
	}

	/**
	 * @brief Compute sum of squared errors between original pixels and decoded block.
	 *
	 * Decompresses the PCB and computes per-texel RGBA error against the original U8 pixels.
	 */
	float compute_sse_vs_original(
		astcenc_context* ctx,
		const uint8_t pcb[16],
		const uint8_t* original_pixels
	) {
		astcenc_contexti& ctxi = ctx->context;

		// Decompress the block
		symbolic_compressed_block scb;
		physical_to_symbolic(*ctxi.bsd, pcb, scb);

		image_block decoded_blk;
		decompress_symbolic_block(
			ctxi.config.profile, *ctxi.bsd,
			0, 0, 0, scb, decoded_blk);

		// Compare against original
		float sse = 0.0f;
		for (unsigned int i = 0; i < TEXEL_COUNT; i++)
		{
			// Decoded block data is in 0-65535 range for LDR
			float dr = decoded_blk.data_r[i] / 65535.0f * 255.0f - static_cast<float>(original_pixels[i * 4 + 0]);
			float dg = decoded_blk.data_g[i] / 65535.0f * 255.0f - static_cast<float>(original_pixels[i * 4 + 1]);
			float db = decoded_blk.data_b[i] / 65535.0f * 255.0f - static_cast<float>(original_pixels[i * 4 + 2]);
			float da = decoded_blk.data_a[i] / 65535.0f * 255.0f - static_cast<float>(original_pixels[i * 4 + 3]);
			sse += dr * dr + dg * dg + db * db + da * da;
		}

		return sse;
	}

	/**
	 * @brief Print diagnostic info about an SCB.
	 */
	static void dump_scb(const char* label, const symbolic_compressed_block& scb)
	{
		fprintf(stderr, "\n=== %s ===\n", label);
		fprintf(stderr, "  block_type=%u  block_mode=%u  partition_count=%u  partition_index=%u\n",
			scb.block_type, scb.block_mode, scb.partition_count, scb.partition_index);
		fprintf(stderr, "  plane2_component=%d  quant_mode=%d  color_formats_matched=%u\n",
			(int)scb.plane2_component, (int)scb.quant_mode, scb.color_formats_matched);
		fprintf(stderr, "  errorval=%.2f\n", scb.errorval);

		for (unsigned int j = 0; j < scb.partition_count; j++)
		{
			fprintf(stderr, "  partition[%u]: color_format=%u (class=%u)  values=[",
				j, scb.color_formats[j], scb.color_formats[j] >> 2);
			for (int k = 0; k < 8; k++)
			{
				fprintf(stderr, "%s%u", k ? "," : "", scb.color_values[j][k]);
			}
			fprintf(stderr, "]\n");
		}

		fprintf(stderr, "  weights[0..15]=[");
		for (int k = 0; k < 16; k++)
		{
			fprintf(stderr, "%s%u", k ? "," : "", scb.weights[k]);
		}
		fprintf(stderr, "]\n");

		if (scb.plane2_component >= 0)
		{
			fprintf(stderr, "  weights_p2[0..15]=[");
			for (int k = 0; k < 16; k++)
			{
				fprintf(stderr, "%s%u", k ? "," : "", scb.weights[k + WEIGHTS_PLANE2_OFFSET]);
			}
			fprintf(stderr, "]\n");
		}
	}

	// ---------- Test pixel pattern generators ----------

	void generate_rgb_gradient(uint8_t* pixels)
	{
		for (unsigned int y = 0; y < BLOCK_Y; y++)
		{
			for (unsigned int x = 0; x < BLOCK_X; x++)
			{
				unsigned int idx = (y * BLOCK_X + x) * 4;
				pixels[idx + 0] = static_cast<uint8_t>(x * 32);    // R
				pixels[idx + 1] = static_cast<uint8_t>(y * 32);    // G
				pixels[idx + 2] = static_cast<uint8_t>((x + y) * 16); // B
				pixels[idx + 3] = 255;                              // A
			}
		}
	}

	void generate_decorrelated_alpha(uint8_t* pixels)
	{
		// RGB varies in one direction, alpha varies orthogonally → forces dual-plane
		for (unsigned int y = 0; y < BLOCK_Y; y++)
		{
			for (unsigned int x = 0; x < BLOCK_X; x++)
			{
				unsigned int idx = (y * BLOCK_X + x) * 4;
				pixels[idx + 0] = static_cast<uint8_t>(x * 32);    // R gradient horizontal
				pixels[idx + 1] = static_cast<uint8_t>(x * 32);    // G gradient horizontal
				pixels[idx + 2] = static_cast<uint8_t>(x * 32);    // B gradient horizontal
				pixels[idx + 3] = static_cast<uint8_t>(y * 32);    // A gradient vertical (decorrelated)
			}
		}
	}

	void generate_grayscale_ramp(uint8_t* pixels)
	{
		for (unsigned int y = 0; y < BLOCK_Y; y++)
		{
			for (unsigned int x = 0; x < BLOCK_X; x++)
			{
				unsigned int idx = (y * BLOCK_X + x) * 4;
				uint8_t val = static_cast<uint8_t>((x + y) * 16);
				pixels[idx + 0] = val;  // R
				pixels[idx + 1] = val;  // G
				pixels[idx + 2] = val;  // B
				pixels[idx + 3] = 255;  // A
			}
		}
	}

	void generate_random_rgba(uint8_t* pixels, unsigned int seed = 42)
	{
		srand(seed);
		for (unsigned int i = 0; i < TEXEL_COUNT * 4; i++)
		{
			pixels[i] = static_cast<uint8_t>(rand() % 256);
		}
	}
};

const astcenc_swizzle GuidedCompressionTest::swizzle {
	ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
};

// ============================================================================
// Test A: Fastest self-guide (sanity baseline)
// ============================================================================

TEST_F(GuidedCompressionTest, RoundTrip_Fastest_SelfGuide)
{
	uint8_t pixels[TEXEL_COUNT * 4];
	generate_rgb_gradient(pixels);

	image_block blk;
	fill_image_block(ctx_fastest, pixels, blk);

	// Normal compress
	uint8_t pcb1[16];
	symbolic_compressed_block scb1;
	compress_and_decode(ctx_fastest, blk, pcb1, scb1);

	// Skip if constant-color block
	if (scb1.block_type != SYM_BTYPE_NONCONST)
	{
		GTEST_SKIP() << "Normal encode produced constant-color block";
	}

	// Extract guide params
	unsigned int guide_bm, guide_pc, guide_pi;
	uint8_t guide_fc[BLOCK_MAX_PARTITIONS];
	extract_guide_params(scb1, guide_bm, guide_pc, guide_pi, guide_fc);

	// Guided compress using same context
	uint8_t pcb2[16];
	symbolic_compressed_block scb2;
	guided_compress_and_decode(
		ctx_fastest, blk,
		guide_bm, guide_pc, guide_pi, guide_fc,
		pcb2, scb2);

	ASSERT_EQ(scb2.block_type, SYM_BTYPE_NONCONST)
		<< "Guided encode produced non-NONCONST block";

	float sse1 = compute_sse_vs_original(ctx_fastest, pcb1, pixels);
	float sse2 = compute_sse_vs_original(ctx_fastest, pcb2, pixels);

	fprintf(stderr, "\n[Fastest Self-Guide] normal_sse=%.2f  guided_sse=%.2f  ratio=%.3f\n",
		sse1, sse2, sse2 / (sse1 + 1e-10f));

	if (sse2 > sse1 * 1.5f)
	{
		dump_scb("Normal (fastest)", scb1);
		dump_scb("Guided (fastest)", scb2);
	}

	// Guided should be within 50% of normal for self-guide
	EXPECT_LE(sse2, sse1 * 1.5f + 1.0f)
		<< "Guided encode is significantly worse than normal encode on same preset";
}

// ============================================================================
// Test B: Thorough self-guide (the critical test)
// ============================================================================

TEST_F(GuidedCompressionTest, RoundTrip_Thorough_SelfGuide)
{
	uint8_t pixels[TEXEL_COUNT * 4];
	generate_rgb_gradient(pixels);

	image_block blk;
	fill_image_block(ctx_thorough, pixels, blk);

	// Normal compress with thorough
	uint8_t pcb1[16];
	symbolic_compressed_block scb1;
	compress_and_decode(ctx_thorough, blk, pcb1, scb1);

	if (scb1.block_type != SYM_BTYPE_NONCONST)
	{
		GTEST_SKIP() << "Normal encode produced constant-color block";
	}

	unsigned int guide_bm, guide_pc, guide_pi;
	uint8_t guide_fc[BLOCK_MAX_PARTITIONS];
	extract_guide_params(scb1, guide_bm, guide_pc, guide_pi, guide_fc);

	// Guided compress with same thorough context
	uint8_t pcb2[16];
	symbolic_compressed_block scb2;
	guided_compress_and_decode(
		ctx_thorough, blk,
		guide_bm, guide_pc, guide_pi, guide_fc,
		pcb2, scb2);

	ASSERT_EQ(scb2.block_type, SYM_BTYPE_NONCONST)
		<< "Guided encode produced non-NONCONST block";

	float sse1 = compute_sse_vs_original(ctx_thorough, pcb1, pixels);
	float sse2 = compute_sse_vs_original(ctx_thorough, pcb2, pixels);

	fprintf(stderr, "\n[Thorough Self-Guide] normal_sse=%.2f  guided_sse=%.2f  ratio=%.3f\n",
		sse1, sse2, sse2 / (sse1 + 1e-10f));

	// Always dump diagnostics for the critical test
	dump_scb("Normal (thorough)", scb1);
	dump_scb("Guided (thorough)", scb2);

	// Guided should be within 2x of normal when using same preset
	EXPECT_LE(sse2, sse1 * 2.0f + 1.0f)
		<< "Guided encode quality is significantly degraded with thorough self-guide";
}

// ============================================================================
// Test C: Thorough 1-plane isolation
// ============================================================================

TEST_F(GuidedCompressionTest, RoundTrip_Thorough_1Plane)
{
	// Grayscale ramp forces luminance format, highly correlated channels → 1 plane
	uint8_t pixels[TEXEL_COUNT * 4];
	generate_grayscale_ramp(pixels);

	image_block blk;
	fill_image_block(ctx_thorough, pixels, blk);

	uint8_t pcb1[16];
	symbolic_compressed_block scb1;
	compress_and_decode(ctx_thorough, blk, pcb1, scb1);

	if (scb1.block_type != SYM_BTYPE_NONCONST)
	{
		GTEST_SKIP() << "Normal encode produced constant-color block";
	}

	fprintf(stderr, "\n[1-Plane Test] normal: plane2_component=%d\n", (int)scb1.plane2_component);

	ASSERT_EQ(scb1.plane2_component, -1)
		<< "Expected grayscale ramp to encode as 1-plane";

	unsigned int guide_bm, guide_pc, guide_pi;
	uint8_t guide_fc[BLOCK_MAX_PARTITIONS];
	extract_guide_params(scb1, guide_bm, guide_pc, guide_pi, guide_fc);

	uint8_t pcb2[16];
	symbolic_compressed_block scb2;
	guided_compress_and_decode(
		ctx_thorough, blk,
		guide_bm, guide_pc, guide_pi, guide_fc,
		pcb2, scb2);

	ASSERT_EQ(scb2.block_type, SYM_BTYPE_NONCONST)
		<< "Guided encode produced non-NONCONST block";

	float sse1 = compute_sse_vs_original(ctx_thorough, pcb1, pixels);
	float sse2 = compute_sse_vs_original(ctx_thorough, pcb2, pixels);

	fprintf(stderr, "[1-Plane Test] normal_sse=%.2f  guided_sse=%.2f  ratio=%.3f\n",
		sse1, sse2, sse2 / (sse1 + 1e-10f));

	if (sse2 > sse1 * 2.0f)
	{
		dump_scb("Normal (1-plane)", scb1);
		dump_scb("Guided (1-plane)", scb2);
	}

	EXPECT_LE(sse2, sse1 * 2.0f + 1.0f)
		<< "1-plane guided encode quality is significantly degraded";
}

// ============================================================================
// Test D: Thorough 2-plane isolation
// ============================================================================

TEST_F(GuidedCompressionTest, RoundTrip_Thorough_2Plane)
{
	// Decorrelated alpha: constant RGB with alpha gradient → forces dual-plane
	uint8_t pixels[TEXEL_COUNT * 4];
	generate_decorrelated_alpha(pixels);

	image_block blk;
	fill_image_block(ctx_thorough, pixels, blk);

	uint8_t pcb1[16];
	symbolic_compressed_block scb1;
	compress_and_decode(ctx_thorough, blk, pcb1, scb1);

	if (scb1.block_type != SYM_BTYPE_NONCONST)
	{
		GTEST_SKIP() << "Normal encode produced constant-color block";
	}

	fprintf(stderr, "\n[2-Plane Test] normal: plane2_component=%d\n", (int)scb1.plane2_component);

	// This should use dual-plane; if it doesn't, the test still works but skip the 2-plane assertion
	if (scb1.plane2_component < 0)
	{
		fprintf(stderr, "WARNING: Thorough encoder chose 1-plane for decorrelated alpha. "
		        "Test still runs but doesn't isolate dual-plane path.\n");
	}

	unsigned int guide_bm, guide_pc, guide_pi;
	uint8_t guide_fc[BLOCK_MAX_PARTITIONS];
	extract_guide_params(scb1, guide_bm, guide_pc, guide_pi, guide_fc);

	uint8_t pcb2[16];
	symbolic_compressed_block scb2;
	guided_compress_and_decode(
		ctx_thorough, blk,
		guide_bm, guide_pc, guide_pi, guide_fc,
		pcb2, scb2);

	ASSERT_EQ(scb2.block_type, SYM_BTYPE_NONCONST)
		<< "Guided encode produced non-NONCONST block";

	float sse1 = compute_sse_vs_original(ctx_thorough, pcb1, pixels);
	float sse2 = compute_sse_vs_original(ctx_thorough, pcb2, pixels);

	fprintf(stderr, "[2-Plane Test] normal_sse=%.2f  guided_sse=%.2f  ratio=%.3f\n",
		sse1, sse2, sse2 / (sse1 + 1e-10f));

	if (scb1.plane2_component >= 0)
	{
		fprintf(stderr, "  normal plane2_component=%d  guided plane2_component=%d\n",
			(int)scb1.plane2_component, (int)scb2.plane2_component);
	}

	if (sse2 > sse1 * 2.0f)
	{
		dump_scb("Normal (2-plane)", scb1);
		dump_scb("Guided (2-plane)", scb2);
	}

	EXPECT_LE(sse2, sse1 * 2.0f + 1.0f)
		<< "2-plane guided encode quality is significantly degraded";
}

// ============================================================================
// Test E: Detailed diagnostic comparison (the "smoking gun" test)
// ============================================================================

TEST_F(GuidedCompressionTest, ABComparison_DetailedDiag)
{
	uint8_t pixels[TEXEL_COUNT * 4];
	generate_rgb_gradient(pixels);

	image_block blk;
	fill_image_block(ctx_thorough, pixels, blk);

	// Normal thorough encode
	uint8_t pcb1[16];
	symbolic_compressed_block scb1;
	compress_and_decode(ctx_thorough, blk, pcb1, scb1);

	if (scb1.block_type != SYM_BTYPE_NONCONST)
	{
		GTEST_SKIP() << "Normal encode produced constant-color block";
	}

	unsigned int guide_bm, guide_pc, guide_pi;
	uint8_t guide_fc[BLOCK_MAX_PARTITIONS];
	extract_guide_params(scb1, guide_bm, guide_pc, guide_pi, guide_fc);

	// Guided encode with same thorough context
	uint8_t pcb2[16];
	symbolic_compressed_block scb2;
	guided_compress_and_decode(
		ctx_thorough, blk,
		guide_bm, guide_pc, guide_pi, guide_fc,
		pcb2, scb2);

	float sse1 = compute_sse_vs_original(ctx_thorough, pcb1, pixels);
	float sse2 = compute_sse_vs_original(ctx_thorough, pcb2, pixels);

	fprintf(stderr, "\n[AB Comparison] normal_sse=%.2f  guided_sse=%.2f  ratio=%.3f\n",
		sse1, sse2, sse2 / (sse1 + 1e-10f));

	// If quality diverges significantly, dump detailed diagnostics
	if (sse2 > sse1 * 2.0f + 1.0f)
	{
		fprintf(stderr, "\n*** QUALITY DIVERGENCE DETECTED ***\n");
		dump_scb("Normal", scb1);
		dump_scb("Guided", scb2);

		// Print block mode details from the BSD
		const block_size_descriptor& bsd = *ctx_thorough->context.bsd;
		unsigned int packed_idx = bsd.block_mode_packed_index[scb1.block_mode];
		if (packed_idx != BLOCK_BAD_BLOCK_MODE)
		{
			const block_mode& bm = bsd.get_block_mode(scb1.block_mode);
			const auto& di = bsd.get_decimation_info(bm.decimation_mode);
			fprintf(stderr, "\nBlock mode details (mode_index=%u):\n", scb1.block_mode);
			fprintf(stderr, "  is_dual_plane=%d  weight_bits=%u  quant_mode=%u\n",
				(int)bm.is_dual_plane, bm.weight_bits, bm.quant_mode);
			fprintf(stderr, "  weight_count=%u (decimation of %u texels)\n",
				di.weight_count, bsd.texel_count);
		}

		// Compare physical bytes
		fprintf(stderr, "\nPhysical block bytes:\n  normal: ");
		for (int i = 0; i < 16; i++)
		{
			fprintf(stderr, "%02x ", pcb1[i]);
		}
		fprintf(stderr, "\n  guided: ");
		for (int i = 0; i < 16; i++)
		{
			fprintf(stderr, "%02x ", pcb2[i]);
		}
		fprintf(stderr, "\n");

		// Key comparison: same quant_mode?
		fprintf(stderr, "\nQuant mode: normal=%d  guided=%d  %s\n",
			(int)scb1.quant_mode, (int)scb2.quant_mode,
			scb1.quant_mode == scb2.quant_mode ? "MATCH" : "MISMATCH");

		// Key comparison: same format class?
		for (unsigned int j = 0; j < scb1.partition_count; j++)
		{
			fprintf(stderr, "Format[%u]: normal=%u (class %u)  guided=%u (class %u)  %s\n",
				j, scb1.color_formats[j], scb1.color_formats[j] >> 2,
				scb2.color_formats[j], scb2.color_formats[j] >> 2,
				(scb1.color_formats[j] >> 2) == (scb2.color_formats[j] >> 2) ? "MATCH" : "MISMATCH");
		}
	}

	// This is a diagnostic test; the assertion is lenient
	EXPECT_LE(sse2, sse1 * 5.0f + 10.0f)
		<< "Guided encode quality is dramatically worse than normal";
}

// ============================================================================
// Test F: Bit budget / quant mode verification
// ============================================================================

TEST_F(GuidedCompressionTest, BitBudget_QuantMode)
{
	// Test multiple pixel patterns to exercise different block modes
	uint8_t pixels_grad[TEXEL_COUNT * 4];
	uint8_t pixels_gray[TEXEL_COUNT * 4];
	uint8_t pixels_rand[TEXEL_COUNT * 4];
	generate_rgb_gradient(pixels_grad);
	generate_grayscale_ramp(pixels_gray);
	generate_random_rgba(pixels_rand);

	const uint8_t* pixel_sets[] = { pixels_grad, pixels_gray, pixels_rand };
	const char* set_names[] = { "gradient", "grayscale", "random" };

	for (int s = 0; s < 3; s++)
	{
		image_block blk;
		fill_image_block(ctx_thorough, pixel_sets[s], blk);

		uint8_t pcb1[16];
		symbolic_compressed_block scb1;
		compress_and_decode(ctx_thorough, blk, pcb1, scb1);

		if (scb1.block_type != SYM_BTYPE_NONCONST)
		{
			continue;
		}

		unsigned int guide_bm, guide_pc, guide_pi;
		uint8_t guide_fc[BLOCK_MAX_PARTITIONS];
		extract_guide_params(scb1, guide_bm, guide_pc, guide_pi, guide_fc);

		// Guided encode
		uint8_t pcb2[16];
		symbolic_compressed_block scb2;
		guided_compress_and_decode(
			ctx_thorough, blk,
			guide_bm, guide_pc, guide_pi, guide_fc,
			pcb2, scb2);

		if (scb2.block_type != SYM_BTYPE_NONCONST)
		{
			fprintf(stderr, "[BitBudget %s] Guided produced non-NONCONST block\n", set_names[s]);
			EXPECT_EQ(scb2.block_type, SYM_BTYPE_NONCONST);
			continue;
		}

		// The block_mode in the guided output should match the guide
		EXPECT_EQ(scb2.block_mode, scb1.block_mode)
			<< "Block mode mismatch for " << set_names[s];

		// The partition count and index should match
		EXPECT_EQ(scb2.partition_count, scb1.partition_count)
			<< "Partition count mismatch for " << set_names[s];

		if (scb1.partition_count >= 2)
		{
			EXPECT_EQ(scb2.partition_index, scb1.partition_index)
				<< "Partition index mismatch for " << set_names[s];
		}

		// The format class should match (exact format may differ within a class)
		for (unsigned int j = 0; j < scb1.partition_count; j++)
		{
			uint8_t class1 = scb1.color_formats[j] >> 2;
			uint8_t class2 = scb2.color_formats[j] >> 2;
			EXPECT_EQ(class2, class1)
				<< "Format class mismatch for " << set_names[s] << " partition " << j;
		}

		float sse1 = compute_sse_vs_original(ctx_thorough, pcb1, pixel_sets[s]);
		float sse2 = compute_sse_vs_original(ctx_thorough, pcb2, pixel_sets[s]);

		fprintf(stderr, "[BitBudget %s] bm=%u pc=%u dp=%d quant: normal=%d guided=%d  "
		        "sse: normal=%.2f guided=%.2f ratio=%.3f\n",
			set_names[s], scb1.block_mode, scb1.partition_count,
			(int)scb1.plane2_component,
			(int)scb1.quant_mode, (int)scb2.quant_mode,
			sse1, sse2, sse2 / (sse1 + 1e-10f));
	}
}

// ============================================================================
// Test G: Cross-preset guide (thorough → fastest) — the real use case
// ============================================================================

TEST_F(GuidedCompressionTest, CrossPreset_ThoroughToFastest)
{
	uint8_t pixels[TEXEL_COUNT * 4];
	generate_rgb_gradient(pixels);

	// Compress with thorough to get guide params
	image_block blk_thorough;
	fill_image_block(ctx_thorough, pixels, blk_thorough);

	uint8_t pcb_thorough[16];
	symbolic_compressed_block scb_thorough;
	compress_and_decode(ctx_thorough, blk_thorough, pcb_thorough, scb_thorough);

	if (scb_thorough.block_type != SYM_BTYPE_NONCONST)
	{
		GTEST_SKIP() << "Thorough encode produced constant-color block";
	}

	unsigned int guide_bm, guide_pc, guide_pi;
	uint8_t guide_fc[BLOCK_MAX_PARTITIONS];
	extract_guide_params(scb_thorough, guide_bm, guide_pc, guide_pi, guide_fc);

	// Guided compress with fastest context using thorough's guide params
	image_block blk_fastest;
	fill_image_block(ctx_fastest, pixels, blk_fastest);

	uint8_t pcb_guided[16];
	symbolic_compressed_block scb_guided;
	guided_compress_and_decode(
		ctx_fastest, blk_fastest,
		guide_bm, guide_pc, guide_pi, guide_fc,
		pcb_guided, scb_guided);

	float sse_thorough = compute_sse_vs_original(ctx_thorough, pcb_thorough, pixels);
	float sse_guided = compute_sse_vs_original(ctx_fastest, pcb_guided, pixels);

	// Also get baseline: fastest without guide
	uint8_t pcb_fastest[16];
	symbolic_compressed_block scb_fastest;
	compress_and_decode(ctx_fastest, blk_fastest, pcb_fastest, scb_fastest);
	float sse_fastest = compute_sse_vs_original(ctx_fastest, pcb_fastest, pixels);

	fprintf(stderr, "\n[Cross-Preset] thorough_sse=%.2f  fastest_sse=%.2f  guided_sse=%.2f\n",
		sse_thorough, sse_fastest, sse_guided);
	fprintf(stderr, "  guided/thorough ratio=%.3f  fastest/thorough ratio=%.3f\n",
		sse_guided / (sse_thorough + 1e-10f),
		sse_fastest / (sse_thorough + 1e-10f));

	if (scb_guided.block_type == SYM_BTYPE_NONCONST)
	{
		dump_scb("Thorough (guide source)", scb_thorough);
		dump_scb("Guided (fastest ctx)", scb_guided);
	}
	else
	{
		fprintf(stderr, "WARNING: Guided encode produced block_type=%u (not NONCONST)\n",
			scb_guided.block_type);
		dump_scb("Thorough (guide source)", scb_thorough);
	}

	// Guided with thorough guide should be better than unguided fastest
	// (or at least not catastrophically worse than thorough)
	EXPECT_LE(sse_guided, sse_thorough * 5.0f + 10.0f)
		<< "Cross-preset guided encode is catastrophically worse than thorough baseline";
}

// ============================================================================
// Pipeline-level tests: guide generation → consumption round-trip
// These test the full astcenc_generate_guide / astcenc_compress_image_guided
// path to find bugs in packing, multi-block dispatch, or state pollution.
// ============================================================================

/**
 * @brief Test fixture for full-pipeline guided compression.
 *
 * Uses larger images (multi-block) to exercise the guide file
 * generation and consumption APIs end-to-end.
 */
class GuidedPipelineTest : public ::testing::Test
{
protected:
	static constexpr unsigned int BLOCK_X = 8;
	static constexpr unsigned int BLOCK_Y = 8;
	// 2x2 grid of blocks = 16x16 image
	static constexpr unsigned int IMG_W = 16;
	static constexpr unsigned int IMG_H = 16;
	static constexpr unsigned int BLOCK_COUNT = (IMG_W / BLOCK_X) * (IMG_H / BLOCK_Y);

	static const astcenc_swizzle swizzle;

	astcenc_context* ctx_thorough = nullptr;
	astcenc_context* ctx_fastest = nullptr;

	// Pixel data for the test image
	uint8_t pixels[IMG_W * IMG_H * 4];

	// Slice pointer for astcenc_image (must outlive all uses of the image)
	void* slice_ptr = nullptr;

	void SetUp() override
	{
		astcenc_error status;
		astcenc_config config;

		status = astcenc_config_init(
			ASTCENC_PRF_LDR, BLOCK_X, BLOCK_Y, 1,
			ASTCENC_PRE_THOROUGH, 0, &config);
		ASSERT_EQ(status, ASTCENC_SUCCESS);

		status = astcenc_context_alloc(&config, 1, &ctx_thorough);
		ASSERT_EQ(status, ASTCENC_SUCCESS);

		status = astcenc_config_init(
			ASTCENC_PRF_LDR, BLOCK_X, BLOCK_Y, 1,
			ASTCENC_PRE_FASTEST, 0, &config);
		ASSERT_EQ(status, ASTCENC_SUCCESS);

		status = astcenc_context_alloc(&config, 1, &ctx_fastest);
		ASSERT_EQ(status, ASTCENC_SUCCESS);
	}

	void TearDown() override
	{
		if (ctx_thorough)
		{
			astcenc_context_free(ctx_thorough);
		}
		if (ctx_fastest)
		{
			astcenc_context_free(ctx_fastest);
		}
	}

	/**
	 * @brief Set up an astcenc_image pointing at the pixel buffer.
	 */
	void make_image(astcenc_image& img, uint8_t* data)
	{
		img.dim_x = IMG_W;
		img.dim_y = IMG_H;
		img.dim_z = 1;
		img.data_type = ASTCENC_TYPE_U8;
		// Store slice pointer in the fixture so it outlives the call
		slice_ptr = data;
		img.data = &slice_ptr;
	}

	/**
	 * @brief Fill a 16x16 image with varied content per block.
	 *
	 * Each 8x8 block gets a different pattern to exercise different
	 * block modes and ensure cross-block pollution would be visible.
	 */
	void generate_varied_image(uint8_t* pix)
	{
		for (unsigned int y = 0; y < IMG_H; y++)
		{
			for (unsigned int x = 0; x < IMG_W; x++)
			{
				unsigned int bx = x / BLOCK_X;
				unsigned int by = y / BLOCK_Y;
				unsigned int block_id = by * (IMG_W / BLOCK_X) + bx;
				unsigned int lx = x % BLOCK_X;
				unsigned int ly = y % BLOCK_Y;
				unsigned int idx = (y * IMG_W + x) * 4;

				switch (block_id)
				{
				case 0: // RGB gradient
					pix[idx + 0] = static_cast<uint8_t>(lx * 32);
					pix[idx + 1] = static_cast<uint8_t>(ly * 32);
					pix[idx + 2] = static_cast<uint8_t>((lx + ly) * 16);
					pix[idx + 3] = 255;
					break;
				case 1: // Decorrelated alpha (force dual-plane)
					pix[idx + 0] = static_cast<uint8_t>(lx * 32);
					pix[idx + 1] = static_cast<uint8_t>(lx * 32);
					pix[idx + 2] = static_cast<uint8_t>(lx * 32);
					pix[idx + 3] = static_cast<uint8_t>(ly * 32);
					break;
				case 2: // Grayscale ramp
				{
					uint8_t val = static_cast<uint8_t>((lx + ly) * 16);
					pix[idx + 0] = val;
					pix[idx + 1] = val;
					pix[idx + 2] = val;
					pix[idx + 3] = 255;
					break;
				}
				case 3: // Random-ish (deterministic)
				{
					unsigned int seed = lx * 37 + ly * 113 + 7;
					pix[idx + 0] = static_cast<uint8_t>((seed * 1103515245 + 12345) >> 16);
					pix[idx + 1] = static_cast<uint8_t>((seed * 1103515245 + 54321) >> 16);
					pix[idx + 2] = static_cast<uint8_t>((seed * 1103515245 + 98765) >> 16);
					pix[idx + 3] = static_cast<uint8_t>((seed * 1103515245 + 11111) >> 16);
					break;
				}
				default:
					pix[idx + 0] = 128;
					pix[idx + 1] = 128;
					pix[idx + 2] = 128;
					pix[idx + 3] = 255;
					break;
				}
			}
		}
	}

	/**
	 * @brief Compute per-pixel MSE between decompressed ASTC and original.
	 */
	float compute_image_mse(
		astcenc_context* ctx,
		const uint8_t* astc_data,
		const uint8_t* original
	) {
		// Decompress
		uint8_t decoded[IMG_W * IMG_H * 4];
		astcenc_image dec_img;
		dec_img.dim_x = IMG_W;
		dec_img.dim_y = IMG_H;
		dec_img.dim_z = 1;
		dec_img.data_type = ASTCENC_TYPE_U8;
		uint8_t* dec_slice = decoded;
		dec_img.data = reinterpret_cast<void**>(&dec_slice);

		astcenc_error status = astcenc_decompress_image(
			ctx, astc_data, BLOCK_COUNT * 16,
			&dec_img, &swizzle, 0);

		if (status != ASTCENC_SUCCESS)
		{
			return 1e30f;
		}

		float sse = 0.0f;
		for (unsigned int i = 0; i < IMG_W * IMG_H * 4; i++)
		{
			float d = static_cast<float>(decoded[i]) - static_cast<float>(original[i]);
			sse += d * d;
		}

		return sse / static_cast<float>(IMG_W * IMG_H);
	}
};

const astcenc_swizzle GuidedPipelineTest::swizzle {
	ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
};

// ============================================================================
// Test H: Guide record pack/unpack round-trip
//
// Compress with thorough, generate guide, then verify each guide record
// matches the SCB from the original encode. Catches bit-packing bugs.
// ============================================================================

TEST_F(GuidedPipelineTest, GuideRecordRoundTrip)
{
	generate_varied_image(pixels);

	astcenc_image img;
	make_image(img, pixels);

	// Compress with thorough
	uint8_t astc_data[BLOCK_COUNT * 16];
	astcenc_error status = astcenc_compress_image(
		ctx_thorough, &img, &swizzle,
		astc_data, sizeof(astc_data), 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// Generate guide
	size_t guide_len = 0;
	status = astcenc_generate_guide(
		ctx_thorough, &img, astc_data, sizeof(astc_data),
		nullptr, &guide_len);
	ASSERT_EQ(status, ASTCENC_SUCCESS);
	ASSERT_GT(guide_len, static_cast<size_t>(0));

	std::vector<uint8_t> guide_data(guide_len);
	status = astcenc_generate_guide(
		ctx_thorough, &img, astc_data, sizeof(astc_data),
		guide_data.data(), &guide_len);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// Verify header
	uint32_t magic = static_cast<uint32_t>(guide_data[0])
	               | (static_cast<uint32_t>(guide_data[1]) << 8)
	               | (static_cast<uint32_t>(guide_data[2]) << 16)
	               | (static_cast<uint32_t>(guide_data[3]) << 24);
	EXPECT_EQ(magic, 0x44475341u) << "Guide magic mismatch";
	EXPECT_EQ(guide_data[4], 1u) << "Guide version mismatch";
	EXPECT_EQ(guide_data[5], BLOCK_X) << "Guide block_x mismatch";
	EXPECT_EQ(guide_data[6], BLOCK_Y) << "Guide block_y mismatch";

	// For each block, decode the ASTC data to SCB and compare with guide record
	const block_size_descriptor& bsd = *ctx_thorough->context.bsd;
	const uint8_t* records = guide_data.data() + 24; // GUIDE_HEADER_SIZE

	fprintf(stderr, "\n[GuideRecordRoundTrip] Checking %u blocks\n", BLOCK_COUNT);

	for (unsigned int i = 0; i < BLOCK_COUNT; i++)
	{
		// Decode original ASTC block
		symbolic_compressed_block scb;
		physical_to_symbolic(bsd, astc_data + i * 16, scb);

		// Read guide record
		const uint8_t* rec = records + i * 4;
		uint32_t packed = static_cast<uint32_t>(rec[0])
		                | (static_cast<uint32_t>(rec[1]) << 8)
		                | (static_cast<uint32_t>(rec[2]) << 16)
		                | (static_cast<uint32_t>(rec[3]) << 24);

		bool guide_is_constant = (packed >> 31) & 1;
		bool scb_is_constant = (scb.block_type == SYM_BTYPE_CONST_F16) ||
		                       (scb.block_type == SYM_BTYPE_CONST_U16) ||
		                       (scb.block_type == SYM_BTYPE_ERROR);

		EXPECT_EQ(guide_is_constant, scb_is_constant)
			<< "Block " << i << ": constant flag mismatch";

		if (!guide_is_constant && !scb_is_constant)
		{
			unsigned int guide_bm = (packed >> 20) & 0x7FF;
			unsigned int guide_pc = ((packed >> 18) & 0x3) + 1;
			unsigned int guide_pi = (packed >> 8) & 0x3FF;
			uint8_t guide_fc[4];
			guide_fc[0] = packed & 0x3;
			guide_fc[1] = (packed >> 2) & 0x3;
			guide_fc[2] = (packed >> 4) & 0x3;
			guide_fc[3] = (packed >> 6) & 0x3;

			EXPECT_EQ(guide_bm, scb.block_mode)
				<< "Block " << i << ": block_mode mismatch";

			EXPECT_EQ(guide_pc, scb.partition_count)
				<< "Block " << i << ": partition_count mismatch";

			if (scb.partition_count > 1)
			{
				EXPECT_EQ(guide_pi, scb.partition_index)
					<< "Block " << i << ": partition_index mismatch";
			}

			for (unsigned int j = 0; j < scb.partition_count; j++)
			{
				uint8_t expected_fc = scb.color_formats[j] >> 2;
				EXPECT_EQ(guide_fc[j], expected_fc)
					<< "Block " << i << " partition " << j << ": format_class mismatch"
					<< " (guide=" << (int)guide_fc[j] << " expected=" << (int)expected_fc << ")";
			}

			fprintf(stderr, "  block[%u]: bm=%u pc=%u pi=%u dp=%d fc=[%u,%u,%u,%u] OK\n",
				i, guide_bm, guide_pc, guide_pi, (int)scb.plane2_component,
				guide_fc[0], guide_fc[1], guide_fc[2], guide_fc[3]);
		}
		else
		{
			fprintf(stderr, "  block[%u]: constant\n", i);
		}
	}
}

// ============================================================================
// Test I: Full pipeline — thorough self-guide (multi-block)
//
// Compress thorough → generate guide → guided compress thorough → compare.
// This is the exact flow that shows 24.74 dB. If it fails here, the bug
// is in the pipeline glue, not in compress_block_guided.
// ============================================================================

TEST_F(GuidedPipelineTest, FullPipeline_ThoroughSelfGuide)
{
	generate_varied_image(pixels);

	astcenc_image img;
	make_image(img, pixels);

	// Step 1: Normal thorough compress
	uint8_t astc_normal[BLOCK_COUNT * 16];
	astcenc_error status = astcenc_compress_image(
		ctx_thorough, &img, &swizzle,
		astc_normal, sizeof(astc_normal), 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	float mse_normal = compute_image_mse(ctx_thorough, astc_normal, pixels);

	// Step 2: Generate guide from thorough encode
	size_t guide_len = 0;
	status = astcenc_generate_guide(
		ctx_thorough, &img, astc_normal, sizeof(astc_normal),
		nullptr, &guide_len);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	std::vector<uint8_t> guide_data(guide_len);
	status = astcenc_generate_guide(
		ctx_thorough, &img, astc_normal, sizeof(astc_normal),
		guide_data.data(), &guide_len);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// Step 3: Guided compress with thorough context
	astcenc_compress_reset(ctx_thorough);

	uint8_t astc_guided[BLOCK_COUNT * 16];
	status = astcenc_compress_image_guided(
		ctx_thorough, &img, &swizzle,
		guide_data.data(), guide_len,
		astc_guided, sizeof(astc_guided), 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	float mse_guided = compute_image_mse(ctx_thorough, astc_guided, pixels);

	fprintf(stderr, "\n[FullPipeline ThoroughSelfGuide] normal_mse=%.2f  guided_mse=%.2f  ratio=%.3f\n",
		mse_normal, mse_guided, mse_guided / (mse_normal + 1e-10f));

	// Compare per-block using the actual U8 decompressed images
	const block_size_descriptor& bsd = *ctx_thorough->context.bsd;

	// Decompress both to U8 for per-block analysis
	uint8_t dec_normal_full[IMG_W * IMG_H * 4];
	uint8_t dec_guided_full[IMG_W * IMG_H * 4];
	{
		astcenc_image dn_img, dg_img;
		dn_img.dim_x = dg_img.dim_x = IMG_W;
		dn_img.dim_y = dg_img.dim_y = IMG_H;
		dn_img.dim_z = dg_img.dim_z = 1;
		dn_img.data_type = dg_img.data_type = ASTCENC_TYPE_U8;

		void* dn_slice = dec_normal_full;
		void* dg_slice = dec_guided_full;
		dn_img.data = &dn_slice;
		dg_img.data = &dg_slice;

		ASSERT_EQ(astcenc_decompress_image(ctx_thorough, astc_normal,
			sizeof(astc_normal), &dn_img, &swizzle, 0), ASTCENC_SUCCESS);
		ASSERT_EQ(astcenc_decompress_image(ctx_thorough, astc_guided,
			sizeof(astc_guided), &dg_img, &swizzle, 0), ASTCENC_SUCCESS);
	}

	unsigned int bad_blocks = 0;

	for (unsigned int i = 0; i < BLOCK_COUNT; i++)
	{
		symbolic_compressed_block scb_normal, scb_guided;
		physical_to_symbolic(bsd, astc_normal + i * 16, scb_normal);
		physical_to_symbolic(bsd, astc_guided + i * 16, scb_guided);

		// Compute per-block SSE from U8 decompressed data
		float sse_n = 0.0f, sse_g = 0.0f;
		unsigned int bx = (i % (IMG_W / BLOCK_X)) * BLOCK_X;
		unsigned int by = (i / (IMG_W / BLOCK_X)) * BLOCK_Y;

		for (unsigned int ly = 0; ly < BLOCK_Y; ly++)
		{
			for (unsigned int lx = 0; lx < BLOCK_X; lx++)
			{
				unsigned int px = ((by + ly) * IMG_W + (bx + lx)) * 4;
				for (int c = 0; c < 4; c++)
				{
					float orig = static_cast<float>(pixels[px + c]);
					float d_n = static_cast<float>(dec_normal_full[px + c]) - orig;
					float d_g = static_cast<float>(dec_guided_full[px + c]) - orig;
					sse_n += d_n * d_n;
					sse_g += d_g * d_g;
				}
			}
		}

		float ratio = sse_g / (sse_n + 1e-10f);

		bool is_bad = ratio > 1.5f && sse_g > 100.0f;
		if (is_bad)
		{
			bad_blocks++;
		}

		fprintf(stderr, "  block[%u]: normal_sse=%.1f guided_sse=%.1f ratio=%.3f"
		        " bm_n=%u/%u dp_n=%d/%d type_n=%u/%u%s\n",
			i, sse_n, sse_g, ratio,
			scb_normal.block_mode, scb_guided.block_mode,
			(int)scb_normal.plane2_component, (int)scb_guided.plane2_component,
			scb_normal.block_type, scb_guided.block_type,
			is_bad ? " *** BAD ***" : "");

		// Dump SCB details for bad blocks
		if (is_bad)
		{
			fprintf(stderr, "    NORMAL: pc=%u pi=%u qm=%d cfm=%u\n",
				scb_normal.partition_count, scb_normal.partition_index,
				(int)scb_normal.quant_mode, scb_normal.color_formats_matched);
			for (unsigned int j = 0; j < scb_normal.partition_count; j++)
			{
				fprintf(stderr, "      part[%u]: fmt=%u vals=[%u,%u,%u,%u,%u,%u,%u,%u]\n",
					j, scb_normal.color_formats[j],
					scb_normal.color_values[j][0], scb_normal.color_values[j][1],
					scb_normal.color_values[j][2], scb_normal.color_values[j][3],
					scb_normal.color_values[j][4], scb_normal.color_values[j][5],
					scb_normal.color_values[j][6], scb_normal.color_values[j][7]);
			}
			fprintf(stderr, "      weights[0..15]=[");
			for (int k = 0; k < 16; k++)
				fprintf(stderr, "%s%u", k ? "," : "", scb_normal.weights[k]);
			fprintf(stderr, "]\n");

			fprintf(stderr, "    GUIDED: pc=%u pi=%u qm=%d cfm=%u\n",
				scb_guided.partition_count, scb_guided.partition_index,
				(int)scb_guided.quant_mode, scb_guided.color_formats_matched);
			for (unsigned int j = 0; j < scb_guided.partition_count; j++)
			{
				fprintf(stderr, "      part[%u]: fmt=%u vals=[%u,%u,%u,%u,%u,%u,%u,%u]\n",
					j, scb_guided.color_formats[j],
					scb_guided.color_values[j][0], scb_guided.color_values[j][1],
					scb_guided.color_values[j][2], scb_guided.color_values[j][3],
					scb_guided.color_values[j][4], scb_guided.color_values[j][5],
					scb_guided.color_values[j][6], scb_guided.color_values[j][7]);
			}
			fprintf(stderr, "      weights[0..15]=[");
			for (int k = 0; k < 16; k++)
				fprintf(stderr, "%s%u", k ? "," : "", scb_guided.weights[k]);
			fprintf(stderr, "]\n");
		}
	}

	fprintf(stderr, "  bad_blocks: %u / %u\n", bad_blocks, BLOCK_COUNT);

	// Check raw ASTC bytes
	bool astc_identical = true;
	for (unsigned int i = 0; i < BLOCK_COUNT * 16; i++)
	{
		if (astc_normal[i] != astc_guided[i])
		{
			astc_identical = false;
			break;
		}
	}
	fprintf(stderr, "  ASTC bytes identical: %s\n", astc_identical ? "YES" : "NO");

	if (!astc_identical)
	{
		for (unsigned int b = 0; b < BLOCK_COUNT; b++)
		{
			bool block_same = true;
			for (int j = 0; j < 16; j++)
			{
				if (astc_normal[b * 16 + j] != astc_guided[b * 16 + j])
				{
					block_same = false;
				}
			}
			if (!block_same)
			{
				fprintf(stderr, "  block[%u] ASTC differs:\n    normal: ", b);
				for (int j = 0; j < 16; j++)
				{
					fprintf(stderr, "%02x ", astc_normal[b * 16 + j]);
				}
				fprintf(stderr, "\n    guided: ");
				for (int j = 0; j < 16; j++)
				{
					fprintf(stderr, "%02x ", astc_guided[b * 16 + j]);
				}
				fprintf(stderr, "\n");
			}
		}
	}

	// Compare decompressed U8 pixels
	float sse_n_img = 0.0f, sse_g_img = 0.0f;
	unsigned int diff_pixels = 0;
	for (unsigned int i = 0; i < IMG_W * IMG_H * 4; i++)
	{
		float d_n = static_cast<float>(dec_normal_full[i]) - static_cast<float>(pixels[i]);
		float d_g = static_cast<float>(dec_guided_full[i]) - static_cast<float>(pixels[i]);
		sse_n_img += d_n * d_n;
		sse_g_img += d_g * d_g;
		if (dec_normal_full[i] != dec_guided_full[i])
		{
			diff_pixels++;
		}
	}
	fprintf(stderr, "  U8 decompress: normal_sse=%.1f guided_sse=%.1f diff_samples=%u/%u\n",
		sse_n_img, sse_g_img, diff_pixels, IMG_W * IMG_H * 4);
	fprintf(stderr, "  U8 decompress: normal_mse=%.2f guided_mse=%.2f\n",
		sse_n_img / (IMG_W * IMG_H), sse_g_img / (IMG_W * IMG_H));

	// Dump first few differing pixels to understand the nature of the difference
	unsigned int dumped = 0;
	for (unsigned int y = 0; y < IMG_H && dumped < 8; y++)
	{
		for (unsigned int x = 0; x < IMG_W && dumped < 8; x++)
		{
			unsigned int idx = (y * IMG_W + x) * 4;
			bool differs = false;
			for (int c = 0; c < 4; c++)
			{
				if (dec_normal_full[idx + c] != dec_guided_full[idx + c])
				{
					differs = true;
				}
			}
			if (differs)
			{
				unsigned int blk_x = x / BLOCK_X;
				unsigned int blk_y = y / BLOCK_Y;
				unsigned int block_id = blk_y * (IMG_W / BLOCK_X) + blk_x;
				fprintf(stderr, "  pixel[%u,%u] (block %u): orig=(%u,%u,%u,%u) "
				        "normal=(%u,%u,%u,%u) guided=(%u,%u,%u,%u)\n",
					x, y, block_id,
					pixels[idx], pixels[idx+1], pixels[idx+2], pixels[idx+3],
					dec_normal_full[idx], dec_normal_full[idx+1], dec_normal_full[idx+2], dec_normal_full[idx+3],
					dec_guided_full[idx], dec_guided_full[idx+1], dec_guided_full[idx+2], dec_guided_full[idx+3]);
				dumped++;
			}
		}
	}

	// The guided encode should be within 3x MSE of normal (generous)
	float actual_guided_mse = sse_g_img / static_cast<float>(IMG_W * IMG_H);
	float actual_normal_mse = sse_n_img / static_cast<float>(IMG_W * IMG_H);
	EXPECT_LE(actual_guided_mse, actual_normal_mse * 3.0f + 10.0f)
		<< "Full pipeline guided encode is significantly worse than normal";
}

// ============================================================================
// Test J: Full pipeline — thorough → fastest (the real use case)
// ============================================================================

TEST_F(GuidedPipelineTest, FullPipeline_ThoroughToFastest)
{
	generate_varied_image(pixels);

	astcenc_image img;
	make_image(img, pixels);

	// Step 1: Thorough compress (baseline quality)
	uint8_t astc_thorough[BLOCK_COUNT * 16];
	astcenc_error status = astcenc_compress_image(
		ctx_thorough, &img, &swizzle,
		astc_thorough, sizeof(astc_thorough), 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	float mse_thorough = compute_image_mse(ctx_thorough, astc_thorough, pixels);

	// Step 2: Fastest compress (baseline speed)
	uint8_t astc_fastest[BLOCK_COUNT * 16];
	status = astcenc_compress_image(
		ctx_fastest, &img, &swizzle,
		astc_fastest, sizeof(astc_fastest), 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	float mse_fastest = compute_image_mse(ctx_fastest, astc_fastest, pixels);

	// Step 3: Generate guide from thorough encode
	size_t guide_len = 0;
	status = astcenc_generate_guide(
		ctx_thorough, &img, astc_thorough, sizeof(astc_thorough),
		nullptr, &guide_len);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	std::vector<uint8_t> guide_data(guide_len);
	status = astcenc_generate_guide(
		ctx_thorough, &img, astc_thorough, sizeof(astc_thorough),
		guide_data.data(), &guide_len);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// Step 4: Guided compress with fastest context using thorough guide
	astcenc_compress_reset(ctx_fastest);

	uint8_t astc_guided[BLOCK_COUNT * 16];
	status = astcenc_compress_image_guided(
		ctx_fastest, &img, &swizzle,
		guide_data.data(), guide_len,
		astc_guided, sizeof(astc_guided), 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	float mse_guided = compute_image_mse(ctx_fastest, astc_guided, pixels);

	fprintf(stderr, "\n[FullPipeline ThoroughToFastest]\n");
	fprintf(stderr, "  thorough_mse=%.2f  fastest_mse=%.2f  guided_mse=%.2f\n",
		mse_thorough, mse_fastest, mse_guided);
	fprintf(stderr, "  guided/thorough=%.3f  fastest/thorough=%.3f\n",
		mse_guided / (mse_thorough + 1e-10f),
		mse_fastest / (mse_thorough + 1e-10f));

	// Per-block comparison
	const block_size_descriptor& bsd_t = *ctx_thorough->context.bsd;
	const block_size_descriptor& bsd_f = *ctx_fastest->context.bsd;

	unsigned int bad_mode_blocks = 0;
	unsigned int constant_fallback_blocks = 0;

	for (unsigned int i = 0; i < BLOCK_COUNT; i++)
	{
		symbolic_compressed_block scb_t, scb_g;
		physical_to_symbolic(bsd_t, astc_thorough + i * 16, scb_t);
		physical_to_symbolic(bsd_f, astc_guided + i * 16, scb_g);

		// Check if the thorough block mode is available in the fastest BSD
		bool mode_available = true;
		if (scb_t.block_type == SYM_BTYPE_NONCONST)
		{
			unsigned int packed_idx = bsd_f.block_mode_packed_index[scb_t.block_mode];
			mode_available = (packed_idx != BLOCK_BAD_BLOCK_MODE);
			if (!mode_available)
			{
				bad_mode_blocks++;
			}
		}

		bool guided_is_constant = (scb_g.block_type == SYM_BTYPE_CONST_F16) ||
		                          (scb_g.block_type == SYM_BTYPE_CONST_U16);
		if (guided_is_constant && scb_t.block_type == SYM_BTYPE_NONCONST)
		{
			constant_fallback_blocks++;
		}

		fprintf(stderr, "  block[%u]: thorough bm=%u dp=%d type=%u | "
		        "guided bm=%u dp=%d type=%u | mode_avail=%d\n",
			i, scb_t.block_mode, (int)scb_t.plane2_component, scb_t.block_type,
			scb_g.block_mode, (int)scb_g.plane2_component, scb_g.block_type,
			(int)mode_available);
	}

	fprintf(stderr, "  bad_mode_blocks=%u  constant_fallback=%u\n",
		bad_mode_blocks, constant_fallback_blocks);

	// Guided should not be catastrophically worse than thorough
	EXPECT_LE(mse_guided, mse_thorough * 5.0f + 10.0f)
		<< "Cross-preset full pipeline guided encode is catastrophically worse";

	// No blocks should fall back to constant-color due to BAD_BLOCK_MODE
	EXPECT_EQ(constant_fallback_blocks, 0u)
		<< "Some blocks fell back to constant-color because thorough block modes "
		   "are not available in the fastest BSD";
}

// ============================================================================
// Test K: Guide record packing edge cases
//
// Test that specific known block_mode values round-trip through pack/unpack
// correctly, including high bit patterns that could collide with is_constant.
// ============================================================================

TEST_F(GuidedPipelineTest, GuidePackUnpack_EdgeCases)
{
	// Test pack/unpack for various parameter combinations
	struct TestCase
	{
		unsigned int block_mode;      // 11 bits max (0-2047)
		unsigned int partition_count;  // 1-4
		unsigned int partition_index;  // 10 bits max (0-1023)
		uint8_t format_classes[4];     // 2 bits each (0-3)
	};

	TestCase cases[] = {
		// Normal cases
		{ 53,   3, 743, { 2, 2, 2, 0 } },  // 3-partition RGB, like our gradient
		{ 1090, 1, 0,   { 1, 0, 0, 0 } },  // Dual-plane compatible mode
		{ 563,  1, 0,   { 0, 0, 0, 0 } },  // Luminance
		{ 1348, 1, 0,   { 3, 0, 0, 0 } },  // RGBA
		{ 324,  2, 512, { 3, 3, 0, 0 } },  // 2-partition RGBA

		// Edge cases: max values
		{ 2047, 4, 1023, { 3, 3, 3, 3 } },
		{ 0,    1, 0,    { 0, 0, 0, 0 } },
		// High block_mode that uses upper bits
		{ 1965, 1, 0, { 2, 0, 0, 0 } },
		{ 2046, 1, 0, { 3, 0, 0, 0 } },

		// Block modes with bit 10 set (tests upper bit doesn't leak into is_constant)
		{ 1024, 1, 0, { 0, 0, 0, 0 } },
		{ 2047, 1, 0, { 0, 0, 0, 0 } },
	};

	for (const auto& tc : cases)
	{
		// Pack (same logic as astcenc_generate_guide)
		uint32_t format_classes_packed = 0;
		for (unsigned int j = 0; j < tc.partition_count; j++)
		{
			format_classes_packed |= static_cast<uint32_t>(tc.format_classes[j] & 0x3) << (j * 2);
		}

		uint32_t packed = (static_cast<uint32_t>(tc.block_mode) << 20)
		                | (static_cast<uint32_t>(tc.partition_count - 1) << 18)
		                | (static_cast<uint32_t>(tc.partition_index) << 8)
		                | format_classes_packed;

		// Verify is_constant is not set
		EXPECT_EQ(packed >> 31, 0u)
			<< "Non-constant block has is_constant flag set for bm=" << tc.block_mode;

		// Unpack (same logic as compress_image_guided)
		unsigned int unpacked_bm = (packed >> 20) & 0x7FF;
		unsigned int unpacked_pc = ((packed >> 18) & 0x3) + 1;
		unsigned int unpacked_pi = (packed >> 8) & 0x3FF;
		uint8_t unpacked_fc[4];
		unpacked_fc[0] = packed & 0x3;
		unpacked_fc[1] = (packed >> 2) & 0x3;
		unpacked_fc[2] = (packed >> 4) & 0x3;
		unpacked_fc[3] = (packed >> 6) & 0x3;

		EXPECT_EQ(unpacked_bm, tc.block_mode)
			<< "block_mode round-trip failed for bm=" << tc.block_mode;
		EXPECT_EQ(unpacked_pc, tc.partition_count)
			<< "partition_count round-trip failed for bm=" << tc.block_mode;

		if (tc.partition_count > 1)
		{
			EXPECT_EQ(unpacked_pi, tc.partition_index)
				<< "partition_index round-trip failed for bm=" << tc.block_mode;
		}

		for (unsigned int j = 0; j < tc.partition_count; j++)
		{
			EXPECT_EQ(unpacked_fc[j], tc.format_classes[j])
				<< "format_class[" << j << "] round-trip failed for bm=" << tc.block_mode;
		}
	}

	// Test constant-color record
	uint32_t const_packed = 1u << 31;
	EXPECT_EQ((const_packed >> 31) & 1u, 1u) << "Constant flag not set";
}

// ============================================================================
// Test L: State pollution between blocks
//
// Compress blocks in different orders and verify the guided output doesn't
// change. If temp buffer state leaks between blocks, order would matter.
// ============================================================================

TEST_F(GuidedPipelineTest, StatePollution_BlockOrder)
{
	generate_varied_image(pixels);

	astcenc_image img;
	make_image(img, pixels);

	// Compress with thorough
	uint8_t astc_data[BLOCK_COUNT * 16];
	astcenc_error status = astcenc_compress_image(
		ctx_thorough, &img, &swizzle,
		astc_data, sizeof(astc_data), 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// Generate guide
	size_t guide_len = 0;
	status = astcenc_generate_guide(
		ctx_thorough, &img, astc_data, sizeof(astc_data),
		nullptr, &guide_len);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	std::vector<uint8_t> guide_data(guide_len);
	status = astcenc_generate_guide(
		ctx_thorough, &img, astc_data, sizeof(astc_data),
		guide_data.data(), &guide_len);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// Guided compress — run it twice to check determinism
	astcenc_compress_reset(ctx_thorough);

	uint8_t astc_guided1[BLOCK_COUNT * 16];
	status = astcenc_compress_image_guided(
		ctx_thorough, &img, &swizzle,
		guide_data.data(), guide_len,
		astc_guided1, sizeof(astc_guided1), 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	astcenc_compress_reset(ctx_thorough);

	uint8_t astc_guided2[BLOCK_COUNT * 16];
	status = astcenc_compress_image_guided(
		ctx_thorough, &img, &swizzle,
		guide_data.data(), guide_len,
		astc_guided2, sizeof(astc_guided2), 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// The two guided compresses should produce bit-identical output
	bool identical = true;
	for (unsigned int i = 0; i < BLOCK_COUNT; i++)
	{
		bool block_match = true;
		for (int b = 0; b < 16; b++)
		{
			if (astc_guided1[i * 16 + b] != astc_guided2[i * 16 + b])
			{
				block_match = false;
				identical = false;
			}
		}

		if (!block_match)
		{
			fprintf(stderr, "Block %u differs between two guided runs!\n", i);
			fprintf(stderr, "  run1: ");
			for (int b = 0; b < 16; b++)
			{
				fprintf(stderr, "%02x ", astc_guided1[i * 16 + b]);
			}
			fprintf(stderr, "\n  run2: ");
			for (int b = 0; b < 16; b++)
			{
				fprintf(stderr, "%02x ", astc_guided2[i * 16 + b]);
			}
			fprintf(stderr, "\n");
		}
	}

	EXPECT_TRUE(identical) << "Guided compression is non-deterministic across runs";

	// Now compare guided output to normal block-by-block
	// Each block should use the SAME block_mode as in the guide (unless constant)
	const block_size_descriptor& bsd = *ctx_thorough->context.bsd;
	const uint8_t* records = guide_data.data() + 24;

	for (unsigned int i = 0; i < BLOCK_COUNT; i++)
	{
		const uint8_t* rec = records + i * 4;
		uint32_t packed = static_cast<uint32_t>(rec[0])
		                | (static_cast<uint32_t>(rec[1]) << 8)
		                | (static_cast<uint32_t>(rec[2]) << 16)
		                | (static_cast<uint32_t>(rec[3]) << 24);

		bool guide_is_constant = (packed >> 31) & 1;
		if (guide_is_constant)
		{
			continue;
		}

		unsigned int expected_bm = (packed >> 20) & 0x7FF;

		symbolic_compressed_block scb;
		physical_to_symbolic(bsd, astc_guided1 + i * 16, scb);

		if (scb.block_type == SYM_BTYPE_NONCONST)
		{
			EXPECT_EQ(scb.block_mode, expected_bm)
				<< "Block " << i << ": guided output block_mode doesn't match guide record";
		}
	}
}

// ============================================================================
// Test M: Constant-color block pipeline round-trip
//
// Create an image where all blocks are constant-color (all pixels identical
// within each block). Verify the guided pipeline reproduces identical output.
// ============================================================================

TEST_F(GuidedPipelineTest, ConstantColor_AllBlocks)
{
	// Fill entire image with a single constant color
	for (unsigned int i = 0; i < IMG_W * IMG_H * 4; i += 4)
	{
		pixels[i + 0] = 100;
		pixels[i + 1] = 150;
		pixels[i + 2] = 200;
		pixels[i + 3] = 255;
	}

	astcenc_image img;
	make_image(img, pixels);

	// Normal compress
	uint8_t astc_normal[BLOCK_COUNT * 16];
	astcenc_error status = astcenc_compress_image(
		ctx_thorough, &img, &swizzle,
		astc_normal, sizeof(astc_normal), 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// Generate guide
	size_t guide_len = 0;
	status = astcenc_generate_guide(
		ctx_thorough, &img, astc_normal, sizeof(astc_normal),
		nullptr, &guide_len);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	std::vector<uint8_t> guide_data(guide_len);
	status = astcenc_generate_guide(
		ctx_thorough, &img, astc_normal, sizeof(astc_normal),
		guide_data.data(), &guide_len);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// Verify all guide records are constant
	const uint8_t* records = guide_data.data() + 24;
	for (unsigned int i = 0; i < BLOCK_COUNT; i++)
	{
		const uint8_t* rec = records + i * 4;
		uint32_t packed = static_cast<uint32_t>(rec[0])
		                | (static_cast<uint32_t>(rec[1]) << 8)
		                | (static_cast<uint32_t>(rec[2]) << 16)
		                | (static_cast<uint32_t>(rec[3]) << 24);
		EXPECT_TRUE((packed >> 31) & 1)
			<< "Block " << i << ": expected constant in guide";
	}

	// Guided compress
	astcenc_compress_reset(ctx_thorough);

	uint8_t astc_guided[BLOCK_COUNT * 16];
	status = astcenc_compress_image_guided(
		ctx_thorough, &img, &swizzle,
		guide_data.data(), guide_len,
		astc_guided, sizeof(astc_guided), 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// Compare ASTC output — should be bit-identical for constant blocks
	for (unsigned int b = 0; b < BLOCK_COUNT; b++)
	{
		bool block_match = true;
		for (int j = 0; j < 16; j++)
		{
			if (astc_normal[b * 16 + j] != astc_guided[b * 16 + j])
			{
				block_match = false;
			}
		}
		EXPECT_TRUE(block_match) << "Block " << b << ": constant-color block differs";
	}

	// Verify MSE = 0 for guided
	float mse = compute_image_mse(ctx_thorough, astc_guided, pixels);
	EXPECT_FLOAT_EQ(mse, 0.0f) << "Constant-color guided MSE should be 0";
}

// ============================================================================
// Test N: Transparent blocks (all-zero) mixed with gradient blocks
//
// Simulates a sprite atlas: some blocks are fully transparent, others
// have content. Verifies both constant and non-constant paths work.
// ============================================================================

TEST_F(GuidedPipelineTest, MixedConstantNonConstant)
{
	// Block 0: fully transparent (0,0,0,0)
	// Block 1: gradient with content
	// Block 2: fully transparent
	// Block 3: gradient with content
	for (unsigned int y = 0; y < IMG_H; y++)
	{
		for (unsigned int x = 0; x < IMG_W; x++)
		{
			unsigned int bx = x / BLOCK_X;
			unsigned int by = y / BLOCK_Y;
			unsigned int block_id = by * (IMG_W / BLOCK_X) + bx;
			unsigned int lx = x % BLOCK_X;
			unsigned int ly = y % BLOCK_Y;
			unsigned int idx = (y * IMG_W + x) * 4;

			if (block_id % 2 == 0)
			{
				// Fully transparent
				pixels[idx + 0] = 0;
				pixels[idx + 1] = 0;
				pixels[idx + 2] = 0;
				pixels[idx + 3] = 0;
			}
			else
			{
				// Colored gradient
				pixels[idx + 0] = static_cast<uint8_t>(lx * 32);
				pixels[idx + 1] = static_cast<uint8_t>(ly * 32);
				pixels[idx + 2] = static_cast<uint8_t>((lx + ly) * 16);
				pixels[idx + 3] = 255;
			}
		}
	}

	astcenc_image img;
	make_image(img, pixels);

	// Normal compress
	uint8_t astc_normal[BLOCK_COUNT * 16];
	astcenc_error status = astcenc_compress_image(
		ctx_thorough, &img, &swizzle,
		astc_normal, sizeof(astc_normal), 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// Generate guide
	size_t guide_len = 0;
	status = astcenc_generate_guide(
		ctx_thorough, &img, astc_normal, sizeof(astc_normal),
		nullptr, &guide_len);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	std::vector<uint8_t> guide_data(guide_len);
	status = astcenc_generate_guide(
		ctx_thorough, &img, astc_normal, sizeof(astc_normal),
		guide_data.data(), &guide_len);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// Guided compress
	astcenc_compress_reset(ctx_thorough);

	uint8_t astc_guided[BLOCK_COUNT * 16];
	status = astcenc_compress_image_guided(
		ctx_thorough, &img, &swizzle,
		guide_data.data(), guide_len,
		astc_guided, sizeof(astc_guided), 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// Per-block analysis
	const uint8_t* records = guide_data.data() + 24;

	// Decompress both
	uint8_t dec_normal[IMG_W * IMG_H * 4];
	uint8_t dec_guided[IMG_W * IMG_H * 4];
	{
		astcenc_image dn_img, dg_img;
		dn_img.dim_x = dg_img.dim_x = IMG_W;
		dn_img.dim_y = dg_img.dim_y = IMG_H;
		dn_img.dim_z = dg_img.dim_z = 1;
		dn_img.data_type = dg_img.data_type = ASTCENC_TYPE_U8;
		void* dn_slice = dec_normal;
		void* dg_slice = dec_guided;
		dn_img.data = &dn_slice;
		dg_img.data = &dg_slice;
		ASSERT_EQ(astcenc_decompress_image(ctx_thorough, astc_normal,
			sizeof(astc_normal), &dn_img, &swizzle, 0), ASTCENC_SUCCESS);
		ASSERT_EQ(astcenc_decompress_image(ctx_thorough, astc_guided,
			sizeof(astc_guided), &dg_img, &swizzle, 0), ASTCENC_SUCCESS);
	}

	for (unsigned int i = 0; i < BLOCK_COUNT; i++)
	{
		const uint8_t* rec = records + i * 4;
		uint32_t packed = static_cast<uint32_t>(rec[0])
		                | (static_cast<uint32_t>(rec[1]) << 8)
		                | (static_cast<uint32_t>(rec[2]) << 16)
		                | (static_cast<uint32_t>(rec[3]) << 24);
		bool guide_is_const = (packed >> 31) & 1;

		// Compute per-block SSE
		float sse_n = 0.0f, sse_g = 0.0f;
		unsigned int bx = (i % (IMG_W / BLOCK_X)) * BLOCK_X;
		unsigned int by = (i / (IMG_W / BLOCK_X)) * BLOCK_Y;

		for (unsigned int ly = 0; ly < BLOCK_Y; ly++)
		{
			for (unsigned int lx = 0; lx < BLOCK_X; lx++)
			{
				unsigned int px = ((by + ly) * IMG_W + (bx + lx)) * 4;
				for (int c = 0; c < 4; c++)
				{
					float orig = static_cast<float>(pixels[px + c]);
					float d_n = static_cast<float>(dec_normal[px + c]) - orig;
					float d_g = static_cast<float>(dec_guided[px + c]) - orig;
					sse_n += d_n * d_n;
					sse_g += d_g * d_g;
				}
			}
		}

		fprintf(stderr, "  block[%u]: const=%d normal_sse=%.1f guided_sse=%.1f ratio=%.3f\n",
			i, (int)guide_is_const, sse_n, sse_g, sse_g / (sse_n + 1e-10f));

		if (guide_is_const)
		{
			// Constant blocks should be identical
			EXPECT_LE(sse_g, 1.0f)
				<< "Block " << i << ": constant block has non-zero SSE";
		}
		else
		{
			// Non-constant blocks should be within tolerance
			EXPECT_LE(sse_g, sse_n * 3.0f + 10.0f)
				<< "Block " << i << ": guided SSE too high";
		}
	}
}

// ============================================================================
// Test O: Larger image stress test (64x64 = 64 blocks)
//
// Tests more block diversity to catch edge cases that don't appear
// in the small 4-block test image.
// ============================================================================

class GuidedLargeImageTest : public ::testing::Test
{
protected:
	static constexpr unsigned int BLOCK_X = 8;
	static constexpr unsigned int BLOCK_Y = 8;
	static constexpr unsigned int IMG_W = 64;
	static constexpr unsigned int IMG_H = 64;
	static constexpr unsigned int BLOCK_COUNT = (IMG_W / BLOCK_X) * (IMG_H / BLOCK_Y);  // 64

	static const astcenc_swizzle swizzle;

	astcenc_context* ctx_exhaustive = nullptr;
	astcenc_context* ctx_fastest = nullptr;

	uint8_t pixels[IMG_W * IMG_H * 4];
	void* slice_ptr = nullptr;

	void SetUp() override
	{
		astcenc_error status;
		astcenc_config config;

		status = astcenc_config_init(
			ASTCENC_PRF_LDR, BLOCK_X, BLOCK_Y, 1,
			ASTCENC_PRE_EXHAUSTIVE, 0, &config);
		ASSERT_EQ(status, ASTCENC_SUCCESS);

		status = astcenc_context_alloc(&config, 1, &ctx_exhaustive);
		ASSERT_EQ(status, ASTCENC_SUCCESS);

		status = astcenc_config_init(
			ASTCENC_PRF_LDR, BLOCK_X, BLOCK_Y, 1,
			ASTCENC_PRE_FASTEST, 0, &config);
		ASSERT_EQ(status, ASTCENC_SUCCESS);

		status = astcenc_context_alloc(&config, 1, &ctx_fastest);
		ASSERT_EQ(status, ASTCENC_SUCCESS);
	}

	void TearDown() override
	{
		if (ctx_exhaustive) astcenc_context_free(ctx_exhaustive);
		if (ctx_fastest) astcenc_context_free(ctx_fastest);
	}

	void make_image(astcenc_image& img, uint8_t* data)
	{
		img.dim_x = IMG_W;
		img.dim_y = IMG_H;
		img.dim_z = 1;
		img.data_type = ASTCENC_TYPE_U8;
		slice_ptr = data;
		img.data = &slice_ptr;
	}

	/**
	 * @brief Generate a sprite-atlas-like image.
	 *
	 * Mix of constant (transparent), smooth gradient, noisy, and
	 * partially transparent blocks to simulate real sprite atlas content.
	 */
	void generate_sprite_atlas(uint8_t* pix)
	{
		for (unsigned int y = 0; y < IMG_H; y++)
		{
			for (unsigned int x = 0; x < IMG_W; x++)
			{
				unsigned int bx = x / BLOCK_X;
				unsigned int by = y / BLOCK_Y;
				unsigned int block_id = by * (IMG_W / BLOCK_X) + bx;
				unsigned int lx = x % BLOCK_X;
				unsigned int ly = y % BLOCK_Y;
				unsigned int idx = (y * IMG_W + x) * 4;

				unsigned int pattern = block_id % 8;

				switch (pattern)
				{
				case 0: // Fully transparent
					pix[idx + 0] = 0;
					pix[idx + 1] = 0;
					pix[idx + 2] = 0;
					pix[idx + 3] = 0;
					break;
				case 1: // Solid color
					pix[idx + 0] = 180;
					pix[idx + 1] = 60;
					pix[idx + 2] = 30;
					pix[idx + 3] = 255;
					break;
				case 2: // RGB gradient
					pix[idx + 0] = static_cast<uint8_t>(lx * 32);
					pix[idx + 1] = static_cast<uint8_t>(ly * 32);
					pix[idx + 2] = static_cast<uint8_t>((lx + ly) * 16);
					pix[idx + 3] = 255;
					break;
				case 3: // Decorrelated alpha (horizontal RGB, vertical alpha)
					pix[idx + 0] = static_cast<uint8_t>(lx * 32);
					pix[idx + 1] = static_cast<uint8_t>(lx * 32);
					pix[idx + 2] = static_cast<uint8_t>(lx * 32);
					pix[idx + 3] = static_cast<uint8_t>(ly * 32);
					break;
				case 4: // Noisy content (deterministic)
				{
					unsigned int h = (block_id * 7 + lx * 37 + ly * 113 + 42);
					h = ((h * 1103515245u) + 12345u);
					pix[idx + 0] = static_cast<uint8_t>((h >> 16) & 0xFF);
					pix[idx + 1] = static_cast<uint8_t>((h >> 8) & 0xFF);
					pix[idx + 2] = static_cast<uint8_t>(h & 0xFF);
					pix[idx + 3] = 255;
					break;
				}
				case 5: // Grayscale ramp
				{
					uint8_t val = static_cast<uint8_t>((lx + ly) * 16);
					pix[idx + 0] = val;
					pix[idx + 1] = val;
					pix[idx + 2] = val;
					pix[idx + 3] = 255;
					break;
				}
				case 6: // Edge (left half dark, right half bright)
				{
					uint8_t val = (lx < 4) ? 30 : 220;
					pix[idx + 0] = val;
					pix[idx + 1] = val;
					pix[idx + 2] = val;
					pix[idx + 3] = 255;
					break;
				}
				case 7: // Partially transparent (semi-circle of opacity)
				{
					float cx = lx - 3.5f;
					float cy = ly - 3.5f;
					float r = cx * cx + cy * cy;
					uint8_t alpha = (r < 9.0f) ? 255 : 0;
					pix[idx + 0] = static_cast<uint8_t>(lx * 32);
					pix[idx + 1] = static_cast<uint8_t>(ly * 32);
					pix[idx + 2] = 128;
					pix[idx + 3] = alpha;
					break;
				}
				}
			}
		}
	}
};

const astcenc_swizzle GuidedLargeImageTest::swizzle {
	ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A
};

TEST_F(GuidedLargeImageTest, SpriteAtlas_ExhaustiveSelfGuide)
{
	generate_sprite_atlas(pixels);

	astcenc_image img;
	make_image(img, pixels);

	// Exhaustive compress
	uint8_t astc_normal[BLOCK_COUNT * 16];
	astcenc_error status = astcenc_compress_image(
		ctx_exhaustive, &img, &swizzle,
		astc_normal, sizeof(astc_normal), 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// Generate guide
	size_t guide_len = 0;
	status = astcenc_generate_guide(
		ctx_exhaustive, &img, astc_normal, sizeof(astc_normal),
		nullptr, &guide_len);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	std::vector<uint8_t> guide_data(guide_len);
	status = astcenc_generate_guide(
		ctx_exhaustive, &img, astc_normal, sizeof(astc_normal),
		guide_data.data(), &guide_len);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// Guided compress with SAME exhaustive context (self-guide)
	astcenc_compress_reset(ctx_exhaustive);

	uint8_t astc_guided[BLOCK_COUNT * 16];
	status = astcenc_compress_image_guided(
		ctx_exhaustive, &img, &swizzle,
		guide_data.data(), guide_len,
		astc_guided, sizeof(astc_guided), 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// Decompress both
	uint8_t dec_normal[IMG_W * IMG_H * 4];
	uint8_t dec_guided[IMG_W * IMG_H * 4];
	{
		astcenc_image dn_img, dg_img;
		dn_img.dim_x = dg_img.dim_x = IMG_W;
		dn_img.dim_y = dg_img.dim_y = IMG_H;
		dn_img.dim_z = dg_img.dim_z = 1;
		dn_img.data_type = dg_img.data_type = ASTCENC_TYPE_U8;
		void* dn_slice = dec_normal;
		void* dg_slice = dec_guided;
		dn_img.data = &dn_slice;
		dg_img.data = &dg_slice;
		ASSERT_EQ(astcenc_decompress_image(ctx_exhaustive, astc_normal,
			sizeof(astc_normal), &dn_img, &swizzle, 0), ASTCENC_SUCCESS);
		ASSERT_EQ(astcenc_decompress_image(ctx_exhaustive, astc_guided,
			sizeof(astc_guided), &dg_img, &swizzle, 0), ASTCENC_SUCCESS);
	}

	// Analyze per-block quality
	const block_size_descriptor& bsd = *ctx_exhaustive->context.bsd;
	const uint8_t* records = guide_data.data() + 24;

	unsigned int const_blocks = 0;
	unsigned int good_blocks = 0;
	unsigned int bad_blocks = 0;
	float total_sse_normal = 0.0f;
	float total_sse_guided = 0.0f;

	for (unsigned int i = 0; i < BLOCK_COUNT; i++)
	{
		const uint8_t* rec = records + i * 4;
		uint32_t packed = static_cast<uint32_t>(rec[0])
		                | (static_cast<uint32_t>(rec[1]) << 8)
		                | (static_cast<uint32_t>(rec[2]) << 16)
		                | (static_cast<uint32_t>(rec[3]) << 24);
		bool guide_is_const = (packed >> 31) & 1;

		// Compute per-block SSE
		float sse_n = 0.0f, sse_g = 0.0f;
		unsigned int bx = (i % (IMG_W / BLOCK_X)) * BLOCK_X;
		unsigned int by = (i / (IMG_W / BLOCK_X)) * BLOCK_Y;

		for (unsigned int ly = 0; ly < BLOCK_Y; ly++)
		{
			for (unsigned int lx = 0; lx < BLOCK_X; lx++)
			{
				unsigned int px = ((by + ly) * IMG_W + (bx + lx)) * 4;
				for (int c = 0; c < 4; c++)
				{
					float orig = static_cast<float>(pixels[px + c]);
					float d_n = static_cast<float>(dec_normal[px + c]) - orig;
					float d_g = static_cast<float>(dec_guided[px + c]) - orig;
					sse_n += d_n * d_n;
					sse_g += d_g * d_g;
				}
			}
		}

		total_sse_normal += sse_n;
		total_sse_guided += sse_g;

		if (guide_is_const)
		{
			const_blocks++;
		}

		bool is_bad = (sse_g > sse_n * 5.0f + 50.0f);
		if (is_bad) bad_blocks++;
		else good_blocks++;

		if (is_bad || guide_is_const)
		{
			symbolic_compressed_block scb_n, scb_g;
			physical_to_symbolic(bsd, astc_normal + i * 16, scb_n);
			physical_to_symbolic(bsd, astc_guided + i * 16, scb_g);

			unsigned int guide_bm = (packed >> 20) & 0x7FF;
			unsigned int guide_pc = ((packed >> 18) & 0x3) + 1;
			uint8_t guide_fc[4];
			guide_fc[0] = packed & 0x3;
			guide_fc[1] = (packed >> 2) & 0x3;
			guide_fc[2] = (packed >> 4) & 0x3;
			guide_fc[3] = (packed >> 6) & 0x3;

			fprintf(stderr, "  block[%u] %s: const=%d sse_n=%.1f sse_g=%.1f ratio=%.3f\n",
				i, is_bad ? "BAD" : "CONST", (int)guide_is_const, sse_n, sse_g,
				sse_g / (sse_n + 1e-10f));
			fprintf(stderr, "    guide: bm=%u pc=%u fc=[%u,%u,%u,%u]\n",
				guide_bm, guide_pc, guide_fc[0], guide_fc[1], guide_fc[2], guide_fc[3]);
			fprintf(stderr, "    normal: type=%u bm=%u pc=%u dp=%d qm=%d cfm=%u\n",
				scb_n.block_type, scb_n.block_mode, scb_n.partition_count,
				(int)scb_n.plane2_component, (int)scb_n.quant_mode,
				scb_n.color_formats_matched);
			fprintf(stderr, "    guided: type=%u bm=%u pc=%u dp=%d qm=%d cfm=%u\n",
				scb_g.block_type, scb_g.block_mode, scb_g.partition_count,
				(int)scb_g.plane2_component, (int)scb_g.quant_mode,
				scb_g.color_formats_matched);

			if (is_bad)
			{
				for (unsigned int j = 0; j < scb_n.partition_count; j++)
				{
					fprintf(stderr, "    normal part[%u]: fmt=%u vals=[%u,%u,%u,%u,%u,%u,%u,%u]\n",
						j, scb_n.color_formats[j],
						scb_n.color_values[j][0], scb_n.color_values[j][1],
						scb_n.color_values[j][2], scb_n.color_values[j][3],
						scb_n.color_values[j][4], scb_n.color_values[j][5],
						scb_n.color_values[j][6], scb_n.color_values[j][7]);
				}
				for (unsigned int j = 0; j < scb_g.partition_count; j++)
				{
					fprintf(stderr, "    guided part[%u]: fmt=%u vals=[%u,%u,%u,%u,%u,%u,%u,%u]\n",
						j, scb_g.color_formats[j],
						scb_g.color_values[j][0], scb_g.color_values[j][1],
						scb_g.color_values[j][2], scb_g.color_values[j][3],
						scb_g.color_values[j][4], scb_g.color_values[j][5],
						scb_g.color_values[j][6], scb_g.color_values[j][7]);
				}
			}
		}
	}

	float mse_normal = total_sse_normal / (IMG_W * IMG_H);
	float mse_guided = total_sse_guided / (IMG_W * IMG_H);

	fprintf(stderr, "\n[SpriteAtlas ExhaustiveSelfGuide] %u blocks: %u const, %u good, %u bad\n",
		BLOCK_COUNT, const_blocks, good_blocks, bad_blocks);
	fprintf(stderr, "  mse_normal=%.2f  mse_guided=%.2f  ratio=%.3f\n",
		mse_normal, mse_guided, mse_guided / (mse_normal + 1e-10f));

	EXPECT_EQ(bad_blocks, 0u) << "Some blocks have dramatically worse quality in guided mode";
	EXPECT_LE(mse_guided, mse_normal * 3.0f + 5.0f)
		<< "Overall guided MSE is too high";
}

// ============================================================================
// Test P: Cross-preset sprite atlas (exhaustive → fastest)
//
// The primary use case: use exhaustive guide with fastest re-encode.
// Tests that block modes available in exhaustive but not fastest are handled.
// ============================================================================

TEST_F(GuidedLargeImageTest, SpriteAtlas_ExhaustiveToFastest)
{
	generate_sprite_atlas(pixels);

	astcenc_image img;
	make_image(img, pixels);

	// Exhaustive compress (guide source)
	uint8_t astc_exhaustive[BLOCK_COUNT * 16];
	astcenc_error status = astcenc_compress_image(
		ctx_exhaustive, &img, &swizzle,
		astc_exhaustive, sizeof(astc_exhaustive), 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// Generate guide
	size_t guide_len = 0;
	status = astcenc_generate_guide(
		ctx_exhaustive, &img, astc_exhaustive, sizeof(astc_exhaustive),
		nullptr, &guide_len);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	std::vector<uint8_t> guide_data(guide_len);
	status = astcenc_generate_guide(
		ctx_exhaustive, &img, astc_exhaustive, sizeof(astc_exhaustive),
		guide_data.data(), &guide_len);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// Guided compress with FASTEST context
	uint8_t astc_guided[BLOCK_COUNT * 16];
	status = astcenc_compress_image_guided(
		ctx_fastest, &img, &swizzle,
		guide_data.data(), guide_len,
		astc_guided, sizeof(astc_guided), 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// Also compute unguided fastest for comparison
	uint8_t astc_fastest_normal[BLOCK_COUNT * 16];
	astcenc_compress_reset(ctx_fastest);
	status = astcenc_compress_image(
		ctx_fastest, &img, &swizzle,
		astc_fastest_normal, sizeof(astc_fastest_normal), 0);
	ASSERT_EQ(status, ASTCENC_SUCCESS);

	// Decompress all three
	uint8_t dec_exhaustive[IMG_W * IMG_H * 4];
	uint8_t dec_guided[IMG_W * IMG_H * 4];
	uint8_t dec_fastest[IMG_W * IMG_H * 4];
	{
		auto decompress = [&](astcenc_context* ctx, const uint8_t* astc, uint8_t* dec) {
			astcenc_image d_img;
			d_img.dim_x = IMG_W;
			d_img.dim_y = IMG_H;
			d_img.dim_z = 1;
			d_img.data_type = ASTCENC_TYPE_U8;
			void* d_slice = dec;
			d_img.data = &d_slice;
			return astcenc_decompress_image(ctx, astc, BLOCK_COUNT * 16, &d_img, &swizzle, 0);
		};
		ASSERT_EQ(decompress(ctx_exhaustive, astc_exhaustive, dec_exhaustive), ASTCENC_SUCCESS);
		ASSERT_EQ(decompress(ctx_fastest, astc_guided, dec_guided), ASTCENC_SUCCESS);
		ASSERT_EQ(decompress(ctx_fastest, astc_fastest_normal, dec_fastest), ASTCENC_SUCCESS);
	}

	// Compute image-level MSE
	float sse_exh = 0.0f, sse_guided = 0.0f, sse_fastest = 0.0f;
	for (unsigned int i = 0; i < IMG_W * IMG_H * 4; i++)
	{
		float d_e = static_cast<float>(dec_exhaustive[i]) - static_cast<float>(pixels[i]);
		float d_g = static_cast<float>(dec_guided[i]) - static_cast<float>(pixels[i]);
		float d_f = static_cast<float>(dec_fastest[i]) - static_cast<float>(pixels[i]);
		sse_exh += d_e * d_e;
		sse_guided += d_g * d_g;
		sse_fastest += d_f * d_f;
	}

	float mse_exh = sse_exh / (IMG_W * IMG_H);
	float mse_guided_val = sse_guided / (IMG_W * IMG_H);
	float mse_fastest_val = sse_fastest / (IMG_W * IMG_H);

	fprintf(stderr, "\n[SpriteAtlas ExhaustiveToFastest]\n");
	fprintf(stderr, "  exhaustive_mse=%.2f  fastest_mse=%.2f  guided_mse=%.2f\n",
		mse_exh, mse_fastest_val, mse_guided_val);
	fprintf(stderr, "  guided/exhaustive=%.3f  fastest/exhaustive=%.3f\n",
		mse_guided_val / (mse_exh + 1e-10f),
		mse_fastest_val / (mse_exh + 1e-10f));

	// Count bad blocks (>10x worse than exhaustive per-block)
	const block_size_descriptor& bsd_e = *ctx_exhaustive->context.bsd;
	const block_size_descriptor& bsd_f = *ctx_fastest->context.bsd;
	const uint8_t* records = guide_data.data() + 24;

	unsigned int bad_mode_blocks = 0;
	unsigned int bad_quality_blocks = 0;

	for (unsigned int i = 0; i < BLOCK_COUNT; i++)
	{
		const uint8_t* rec = records + i * 4;
		uint32_t packed = static_cast<uint32_t>(rec[0])
		                | (static_cast<uint32_t>(rec[1]) << 8)
		                | (static_cast<uint32_t>(rec[2]) << 16)
		                | (static_cast<uint32_t>(rec[3]) << 24);
		bool guide_is_const = (packed >> 31) & 1;

		if (guide_is_const)
		{
			continue;
		}

		unsigned int guide_bm = (packed >> 20) & 0x7FF;

		// Check if the exhaustive block mode is available in fastest BSD
		unsigned int fastest_packed = bsd_f.block_mode_packed_index[guide_bm];
		if (fastest_packed == BLOCK_BAD_BLOCK_MODE)
		{
			bad_mode_blocks++;
		}

		// Per-block quality check
		float sse_n = 0.0f, sse_g = 0.0f;
		unsigned int bx = (i % (IMG_W / BLOCK_X)) * BLOCK_X;
		unsigned int by = (i / (IMG_W / BLOCK_X)) * BLOCK_Y;

		for (unsigned int ly = 0; ly < BLOCK_Y; ly++)
		{
			for (unsigned int lx = 0; lx < BLOCK_X; lx++)
			{
				unsigned int px = ((by + ly) * IMG_W + (bx + lx)) * 4;
				for (int c = 0; c < 4; c++)
				{
					float orig = static_cast<float>(pixels[px + c]);
					float d_n = static_cast<float>(dec_exhaustive[px + c]) - orig;
					float d_g = static_cast<float>(dec_guided[px + c]) - orig;
					sse_n += d_n * d_n;
					sse_g += d_g * d_g;
				}
			}
		}

		if (sse_g > sse_n * 10.0f + 100.0f)
		{
			bad_quality_blocks++;
			symbolic_compressed_block scb_e, scb_g;
			physical_to_symbolic(bsd_e, astc_exhaustive + i * 16, scb_e);
			physical_to_symbolic(bsd_f, astc_guided + i * 16, scb_g);

			unsigned int guide_pc = ((packed >> 18) & 0x3) + 1;
			uint8_t guide_fc[4];
			guide_fc[0] = packed & 0x3;
			guide_fc[1] = (packed >> 2) & 0x3;
			guide_fc[2] = (packed >> 4) & 0x3;
			guide_fc[3] = (packed >> 6) & 0x3;

			fprintf(stderr, "  BAD block[%u]: sse_e=%.1f sse_g=%.1f bm=%u dp=%d pc=%u "
			        "fc=[%u,%u,%u,%u] mode_avail=%d\n",
				i, sse_n, sse_g, guide_bm, (int)scb_e.plane2_component,
				guide_pc, guide_fc[0], guide_fc[1], guide_fc[2], guide_fc[3],
				(int)(fastest_packed != BLOCK_BAD_BLOCK_MODE));
			fprintf(stderr, "    exhaustive: type=%u bm=%u qm=%d cfm=%u\n",
				scb_e.block_type, scb_e.block_mode, (int)scb_e.quant_mode,
				scb_e.color_formats_matched);
			fprintf(stderr, "    guided:     type=%u bm=%u qm=%d cfm=%u\n",
				scb_g.block_type, scb_g.block_mode, (int)scb_g.quant_mode,
				scb_g.color_formats_matched);
		}
	}

	fprintf(stderr, "  bad_mode_blocks=%u  bad_quality_blocks=%u\n",
		bad_mode_blocks, bad_quality_blocks);

	// Guided should be significantly better than unguided fastest
	// (or at worst not catastrophically worse than exhaustive)
	EXPECT_LE(mse_guided_val, mse_exh * 5.0f + 20.0f)
		<< "Cross-preset guided encode is catastrophically worse than exhaustive";

	EXPECT_EQ(bad_quality_blocks, 0u)
		<< "Some blocks have >10x worse quality in guided mode";
}

} // namespace astcenc
