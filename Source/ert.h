// SPDX-License-Identifier: Unlicense
// ----------------------------------------------------------------------------
// This is free and unencumbered software released into the public domain.
// Anyone is free to copy, modify, publish, use, compile, sell, or distribute this
// software, either in source code form or as a compiled binary, for any purpose,
// commercial or non - commercial, and by any means.
// In jurisdictions that recognize copyright laws, the author or authors of this
// software dedicate any and all copyright interest in the software to the public
// domain. We make this dedication for the benefit of the public at large and to
// the detriment of our heirs and successors.We intend this dedication to be an
// overt act of relinquishment in perpetuity of all present and future rights to
// this software under copyright law.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
// AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
// ----------------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdio>

// Based on https://github.com/richgel999/bc7enc_rdo/blob/master/ert.h
// With interface tweaks that can make format integration more general & performant.

namespace ert
{
	struct reduce_entropy_params
	{
		// m_lambda: The post-processor tries to reduce distortion*smooth_block_scale + rate*lambda (rate is approximate LZ bits and distortion is scaled MS error multiplied against the smooth block MSE weighting factor).
		// Larger values push the postprocessor towards optimizing more for lower rate, and smaller values more for distortion. 0=minimal distortion.
		float m_lambda;

		// m_lookback_window_size: The number of bytes the encoder can look back from each block to find matches. The larger this value, the slower the encoder but the higher the quality per LZ compressed bit.
		uint32_t m_lookback_window_size;

		// m_max_allowed_rms_increase_ratio: How much the RMS error of a block is allowed to increase before a trial is rejected. 1.0=no increase allowed, 1.05=5% increase allowed, etc.
		float m_max_allowed_rms_increase_ratio;

		float m_max_smooth_block_std_dev;
		float m_smooth_block_max_mse_scale;

		bool m_try_two_matches;
		bool m_allow_relative_movement;
		bool m_skip_zero_mse_blocks;
		bool m_debug_output;

		reduce_entropy_params() { clear(); }

		void clear()
		{
			m_lookback_window_size = 256;
			m_lambda = 1.0f;
			m_max_allowed_rms_increase_ratio = 10.0f;
			m_max_smooth_block_std_dev = 18.0f;
			m_smooth_block_max_mse_scale = 10.0f;
			m_try_two_matches = false;
			m_allow_relative_movement = false;
			m_skip_zero_mse_blocks = false;
			m_debug_output = false;
		}

		void print() const
		{
			printf("lambda: %f\n", (double)m_lambda);
			printf("Lookback window size: %u\n", m_lookback_window_size);
			printf("Max allowed RMS increase ratio: %f\n", (double)m_max_allowed_rms_increase_ratio);
			printf("Max smooth block std dev: %f\n", (double)m_max_smooth_block_std_dev);
			printf("Smooth block max MSE scale: %f\n", (double)m_smooth_block_max_mse_scale);
			printf("Try two matches: %u\n", m_try_two_matches);
			printf("Allow relative movement: %u\n", m_allow_relative_movement);
			printf("Skip zero MSE blocks: %u\n", m_skip_zero_mse_blocks);
		}
	};

	/**
	 * @brief Callback to compute trial block differences.
	 *
	 * All comparing texel values should be in range [0, 255].
	 *
	 * @param[out] out_max_std_dev The channel-wise maximum standard deviation for the original block,
	 * 	can be @c nullptr if not requested.
	 *
	 * @return Should return the mean squared error for current trial block, or any negative value to indicate errors.
	 */
	typedef float diff_block_func_type(void* pUser_data, const uint8_t* pBlock, uint32_t block_index, float* out_max_std_dev);

	// BC7 entropy reduction transform with Deflate/LZMA/LZHAM optimizations
	bool reduce_entropy(uint8_t* pBlock_bytes, uint32_t num_blocks,
		uint32_t total_block_stride_in_bytes, uint32_t block_size_to_optimize_in_bytes,
		const reduce_entropy_params& params, uint32_t& total_modified,
		diff_block_func_type* pDiff_block_func, void* pDiff_block_func_user_data,
		const float* pBlock_mse_scales = nullptr);

} // namespace ert
