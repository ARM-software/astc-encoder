//Copyright Tencent Timi-J1&F1 Studio, Inc. All Rights Reserved

#if !defined(ASTCENC_DECOMPRESS_ONLY)

#include "astcenc_internal_entry.h"
#include "ert.h"

struct astcenc_rdo_context
{
	ert::reduce_entropy_params m_ert_params;
	std::atomic<uint32_t> m_total_modified;

	uint32_t m_xblocks = 0;
	uint32_t m_yblocks = 0;
	uint32_t m_zblocks = 0;
	uint32_t m_total_blocks = 0;
	uint32_t m_image_dimx = 0;
	uint32_t m_image_dimy = 0;
	uint32_t m_image_dimz = 0;
	uint32_t m_bytes_per_texel = 0;
	uint32_t m_component_count = 0;
	astcenc_swizzle m_swizzle;

	std::vector<uint8_t> m_block_pixels;
};

static constexpr uint8_t ASTCENC_RDO_PADDING_PIXEL = 0xff;
static constexpr uint32_t ASTCENC_BYTES_PER_BLOCK = 16;

template<typename T, typename U> T lerp(T a, T b, U s) { return (T)(a + (b - a) * s); }

extern "C" void rdo_progress_emitter(
	float value
) {
	static float previous_value = 100.0f;
	const uint32_t bar_size = 25;
	uint32_t parts = static_cast<uint32_t>(value / 4.0f);

	char buffer[bar_size + 3];
	buffer[0] = '[';

	for (uint32_t i = 0; i < parts; i++)
	{
		buffer[i + 1] = '=';
	}

	for (uint32_t i = parts; i < bar_size; i++)
	{
		buffer[i + 1] = ' ';
	}

	buffer[bar_size + 1] = ']';
	buffer[bar_size + 2] = '\0';

	if (previous_value == 100.0f)
	{
		printf("\n\n");
		printf("Rate-distortion optimization\n");
		printf("============================\n\n");
	}
	previous_value = value;

	printf("    Progress: %s %03.1f%%\r", buffer, static_cast<double>(value));
	fflush(stdout);
}

static uint32_t init_rdo_context(
	astcenc_contexti& ctx,
	const astcenc_image& image,
	const astcenc_swizzle& swizzle
) {
	ctx.rdo_context = new astcenc_rdo_context;
	astcenc_rdo_context& rdo_ctx = *ctx.rdo_context;

	// Perform memory allocations for intermediary buffers
	uint32_t block_dim_x = ctx.bsd->xdim;
	uint32_t block_dim_y = ctx.bsd->ydim;
	uint32_t block_dim_z = ctx.bsd->zdim;
	uint32_t xblocks = (image.dim_x + block_dim_x - 1u) / block_dim_x;
	uint32_t yblocks = (image.dim_y + block_dim_y - 1u) / block_dim_y;
	uint32_t zblocks = (image.dim_z + block_dim_z - 1u) / block_dim_z;
	uint32_t total_blocks = xblocks * yblocks * zblocks;

	uint32_t bytes_per_texel = 4u;
	if (image.data_type == ASTCENC_TYPE_F16) bytes_per_texel *= 2;
	else if (image.data_type == ASTCENC_TYPE_F32) bytes_per_texel *= 4;
	uint32_t texels_per_block = block_dim_x * block_dim_y * block_dim_z;

	rdo_ctx.m_block_pixels.resize(total_blocks * ctx.bsd->texel_count * bytes_per_texel);
	for (uint32_t block_z = 0; block_z < zblocks; ++block_z)
	{
		const uint32_t valid_z = astc::min(image.dim_z - block_z * block_dim_z, block_dim_z);
		const uint32_t padding_z = block_dim_z - valid_z;
	
		for (uint32_t offset_z = 0; offset_z < valid_z; offset_z++)
		{
			const auto* src_data = static_cast<const uint8_t*>(image.data[block_z * block_dim_z + offset_z]);
	
			for (uint32_t block_y = 0, block_idx = 0; block_y < yblocks; ++block_y)
			{
				const uint8_t* src_block_row = src_data + block_y * block_dim_y * image.dim_x * bytes_per_texel;
				const uint32_t valid_y = astc::min(image.dim_y - block_y * block_dim_y, block_dim_y);
				const uint32_t padding_y = block_dim_y - valid_y;
	
				for (uint32_t block_x = 0; block_x < xblocks; ++block_x, ++block_idx)
				{
					const uint8_t* src_block = src_block_row + block_x * block_dim_x * bytes_per_texel;
					const uint32_t valid_x = astc::min(image.dim_x - block_x * block_dim_x, block_dim_x);
					const uint32_t padding_x = block_dim_x - valid_x;
	
					uint8_t* dst_block = rdo_ctx.m_block_pixels.data() + block_idx * texels_per_block * bytes_per_texel;
					for (uint32_t offset_y = 0; offset_y < valid_y; ++offset_y)
					{
						memcpy(dst_block, src_block, valid_x * bytes_per_texel);
						if (padding_x) memset(dst_block + valid_x * bytes_per_texel, ASTCENC_RDO_PADDING_PIXEL, padding_x * bytes_per_texel);
						src_block += image.dim_x * bytes_per_texel;
						dst_block += block_dim_x * bytes_per_texel;
					}
					if (padding_y)
					{
						memset(dst_block, ASTCENC_RDO_PADDING_PIXEL, padding_y * block_dim_x * bytes_per_texel);
						dst_block += padding_y * block_dim_x * bytes_per_texel;
					}
					if (padding_z) memset(dst_block, ASTCENC_RDO_PADDING_PIXEL, padding_z * block_dim_y * block_dim_x * bytes_per_texel);
				}
			}
		}
	}

	// Generate quality parameters
	ert::reduce_entropy_params& ert_params = rdo_ctx.m_ert_params;
	ert_params.m_lambda = astc::clamp(ctx.config.rdo_level, 0.0f, 1.0f);
	ert_params.m_color_weights[0] = static_cast<uint32_t>(ctx.config.cw_r_weight * 255.0f);
	ert_params.m_color_weights[1] = static_cast<uint32_t>(ctx.config.cw_g_weight * 255.0f);
	ert_params.m_color_weights[2] = static_cast<uint32_t>(ctx.config.cw_b_weight * 255.0f);
	ert_params.m_color_weights[3] = static_cast<uint32_t>(ctx.config.cw_a_weight * 255.0f);
	ert_params.m_lookback_window_size = astc::min(ctx.config.rdo_lookback, total_blocks) * ASTCENC_BYTES_PER_BLOCK;

	ert_params.m_try_two_matches = true;
	// ert_params.m_allow_relative_movement = true;
	// ert_params.m_debug_output = true;

	rdo_ctx.m_image_dimx = image.dim_x;
	rdo_ctx.m_image_dimy = image.dim_y;
	rdo_ctx.m_image_dimz = image.dim_z;
	rdo_ctx.m_xblocks = xblocks;
	rdo_ctx.m_yblocks = yblocks;
	rdo_ctx.m_zblocks = zblocks;
	rdo_ctx.m_total_blocks = total_blocks;
	rdo_ctx.m_bytes_per_texel = bytes_per_texel;
	rdo_ctx.m_component_count = 4;
	rdo_ctx.m_swizzle = swizzle;

	return rdo_ctx.m_total_blocks;
}

static bool store_image_block_u8(
	uint8_t* output,
	const image_block& blk,
	const block_size_descriptor& bsd,
	uint32_t x_start,
	uint32_t y_start,
	uint32_t z_start,
	uint32_t x_size,
	uint32_t y_size,
	uint32_t z_size,
	const astcenc_swizzle& swz,
	uint32_t& x_count,
	uint32_t& y_count
) {
	uint32_t x_end = astc::min(x_size, x_start + bsd.xdim);
	x_count = x_end - x_start;
	uint32_t x_nudge = bsd.xdim - x_count;

	uint32_t y_end = astc::min(y_size, y_start + bsd.ydim);
	y_count = y_end - y_start;
	uint32_t y_nudge = (bsd.ydim - y_count) * bsd.xdim;

	uint32_t z_end = astc::min(z_size, z_start + bsd.zdim);

	// True if any non-identity swizzle
	bool needs_swz = (swz.r != ASTCENC_SWZ_R) || (swz.g != ASTCENC_SWZ_G) ||
					 (swz.b != ASTCENC_SWZ_B) || (swz.a != ASTCENC_SWZ_A);

	// True if any swizzle uses Z reconstruct
	bool needs_z = (swz.r == ASTCENC_SWZ_Z) || (swz.g == ASTCENC_SWZ_Z) ||
				   (swz.b == ASTCENC_SWZ_Z) || (swz.a == ASTCENC_SWZ_Z);

	int idx = 0;
	for (uint32_t z = z_start; z < z_end; z++)
	{
		uint8_t* data_slice = output + z * bsd.xdim * bsd.ydim * 4;

		for (uint32_t y = y_start; y < y_end; ++y)
		{
			for (uint32_t x = 0; x < x_count; x += ASTCENC_SIMD_WIDTH)
			{
				uint32_t max_texels = ASTCENC_SIMD_WIDTH;
				uint32_t used_texels = astc::min(x_count - x, max_texels);

				// Unaligned load as rows are not always SIMD_WIDTH long
				vfloat data_r(blk.data_r + idx);
				vfloat data_g(blk.data_g + idx);
				vfloat data_b(blk.data_b + idx);
				vfloat data_a(blk.data_a + idx);

				// Errors are NaN encoded - abort this sample immediately
				// Branch is OK here - it is almost never true so predicts well
				vmask nan_mask = data_r != data_r;
				if (any(nan_mask))
				{
					return false;
				}

				vint data_ri = float_to_int_rtn(min(data_r, 1.0f) * 255.0f);
				vint data_gi = float_to_int_rtn(min(data_g, 1.0f) * 255.0f);
				vint data_bi = float_to_int_rtn(min(data_b, 1.0f) * 255.0f);
				vint data_ai = float_to_int_rtn(min(data_a, 1.0f) * 255.0f);

				if (needs_swz)
				{
					vint swizzle_table[7];
					swizzle_table[ASTCENC_SWZ_0] = vint(0);
					swizzle_table[ASTCENC_SWZ_1] = vint(255);
					swizzle_table[ASTCENC_SWZ_R] = data_ri;
					swizzle_table[ASTCENC_SWZ_G] = data_gi;
					swizzle_table[ASTCENC_SWZ_B] = data_bi;
					swizzle_table[ASTCENC_SWZ_A] = data_ai;

					if (needs_z)
					{
						vfloat data_x = (data_r * vfloat(2.0f)) - vfloat(1.0f);
						vfloat data_y = (data_a * vfloat(2.0f)) - vfloat(1.0f);
						vfloat data_z = vfloat(1.0f) - (data_x * data_x) - (data_y * data_y);
						data_z = max(data_z, 0.0f);
						data_z = (sqrt(data_z) * vfloat(0.5f)) + vfloat(0.5f);

						swizzle_table[ASTCENC_SWZ_Z] = float_to_int_rtn(min(data_z, 1.0f) * 255.0f);
					}

					data_ri = swizzle_table[swz.r];
					data_gi = swizzle_table[swz.g];
					data_bi = swizzle_table[swz.b];
					data_ai = swizzle_table[swz.a];
				}

				vint data_rgbai = interleave_rgba8(data_ri, data_gi, data_bi, data_ai);
				vmask store_mask = vint::lane_id() < vint(used_texels);
				store_lanes_masked(data_slice + idx * 4, data_rgbai, store_mask);

				idx += used_texels;
			}

			if (x_nudge)
			{
				memset(data_slice + idx * 4, ASTCENC_RDO_PADDING_PIXEL, x_nudge * 4);
				idx += x_nudge;
			}
		}

		if (y_nudge)
		{
			memset(data_slice + idx * 4, ASTCENC_RDO_PADDING_PIXEL, y_nudge * 4);
			idx += y_nudge;
		}
	}

	return true;
}

struct local_rdo_context
{
	const astcenc_contexti* ctx;
	uint32_t base_block_idx = 0;
};

static bool unpack_block(
	const void* block,
	ert::color_rgba* pixels,
	uint32_t local_block_idx,
	uint32_t& active_x,
	uint32_t& active_y,
	void* user_data
) {
	const local_rdo_context& local_ctx = *(local_rdo_context*)user_data;
	const astcenc_contexti& ctx = *local_ctx.ctx;
	const astcenc_rdo_context& rdo_ctx = *ctx.rdo_context;
	uint32_t block_idx = local_ctx.base_block_idx + local_block_idx;

	uint32_t block_dim_x = ctx.bsd->xdim;
	uint32_t block_dim_y = ctx.bsd->ydim;
	assert(block_idx < rdo_ctx.m_xblocks * rdo_ctx.m_yblocks);

	uint32_t block_y = block_idx / rdo_ctx.m_xblocks;
	uint32_t block_x = block_idx - block_y * rdo_ctx.m_xblocks;
	assert(block_y * rdo_ctx.m_xblocks + block_x == block_idx);

	uint8_t pcb[ASTCENC_BYTES_PER_BLOCK];
	memcpy(pcb, block, sizeof(pcb));
	symbolic_compressed_block scb;
	physical_to_symbolic(*ctx.bsd, pcb, scb);

	if (scb.block_type == SYM_BTYPE_ERROR)
	{
		return false;
	}

	image_block blk;
	decompress_symbolic_block(ctx.config.profile, *ctx.bsd,
		block_x * block_dim_x, block_y * block_dim_y, 0,
		scb, blk);

	return store_image_block_u8(reinterpret_cast<uint8_t*>(pixels), blk, *ctx.bsd,
		block_x * block_dim_x, block_y * block_dim_y, 0,
		rdo_ctx.m_image_dimx, rdo_ctx.m_image_dimy, 1, rdo_ctx.m_swizzle, active_x, active_y);
}

void rate_distortion_optimize(
	astcenc_context& ctxo,
	const astcenc_image& image,
	const astcenc_swizzle& swizzle,
	uint8_t* buffer
) {
	if (ctxo.context.config.rdo_level == 0.0f || image.data_type != ASTCENC_TYPE_U8 || image.dim_z != 1)
	{
		return;
	}

	// Only the first thread actually runs the initializer
	ctxo.manage_rdo.init([&ctxo, &image, &swizzle]
		{
			return init_rdo_context(ctxo.context, image, swizzle);
		},
		ctxo.context.config.progress_callback ? rdo_progress_emitter : nullptr);

	uint32_t block_dim_x = ctxo.context.bsd->xdim;
	uint32_t block_dim_y = ctxo.context.bsd->ydim;
	uint32_t block_dim_z = ctxo.context.bsd->zdim;
	uint32_t texels_per_block = block_dim_x * block_dim_y * block_dim_z;
	astcenc_rdo_context& rdo_ctx = *ctxo.context.rdo_context;
	uint32_t blocks_per_task = ctxo.context.config.rdo_lookback;
	if (!blocks_per_task) blocks_per_task = rdo_ctx.m_total_blocks; // Effectively single-threaded

	uint32_t total_modified = 0;
	while (true)
	{
		uint32_t count;
		uint32_t base = ctxo.manage_rdo.get_task_assignment(blocks_per_task, count);
		if (!count)
		{
			break;
		}

		local_rdo_context local_context{ &ctxo.context, base };

		ert::reduce_entropy(buffer + base * ASTCENC_BYTES_PER_BLOCK, count, ASTCENC_BYTES_PER_BLOCK, ASTCENC_BYTES_PER_BLOCK,
			block_dim_x, block_dim_y, rdo_ctx.m_component_count, 
			reinterpret_cast<const ert::color_rgba*>(rdo_ctx.m_block_pixels.data() + base * texels_per_block * rdo_ctx.m_bytes_per_texel),
			rdo_ctx.m_ert_params, total_modified, unpack_block, &local_context
		);

		ctxo.manage_rdo.complete_task_assignment(count);
	}

	rdo_ctx.m_total_modified.fetch_add(total_modified, std::memory_order_relaxed);

	// Wait for rdo to complete before freeing memory
	ctxo.manage_rdo.wait();

	// Only the first thread to arrive actually runs the term
	ctxo.manage_rdo.term([&ctxo]
		{
			delete ctxo.context.rdo_context;
			ctxo.context.rdo_context = nullptr;
		});
}

#endif
