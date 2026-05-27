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

#include "ert.h"

#include <cassert>
#include <cstring>

#define ERT_FAVOR_CONT_AND_REP0_MATCHES (1)
#define ERT_FAVOR_REP0_MATCHES (0)
#define ERT_ENABLE_DEBUG (0)

namespace ert
{
	const uint32_t MAX_BLOCK_SIZE_IN_BYTES = 256;
	const uint32_t MIN_MATCH_LEN = 3;
	const float LITERAL_BITS = 13.0f;
	const float MATCH_CONTINUE_BITS = 1.0f;
	const float MATCH_REP0_BITS = 4.0f;

	static inline float clampf(float value, float low, float high) { if (value < low) value = low; else if (value > high) value = high;	return value; }
	template<typename F> inline F lerp(F a, F b, F s) { return a + (b - a) * s; }

	static const uint8_t g_tdefl_small_dist_extra[512] =
	{
		0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
		7, 7, 7, 7, 7, 7, 7, 7
	};

	static const uint8_t g_tdefl_large_dist_extra[128] =
	{
		0, 0, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
		12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13,
		13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13
	};

	static inline uint32_t compute_match_cost_estimate(uint32_t dist, uint32_t match_len_in_bytes)
	{
		assert(match_len_in_bytes <= 258);

		uint32_t len_cost = 6;
		if (match_len_in_bytes >= 12)
			len_cost = 9;
		else if (match_len_in_bytes >= 8)
			len_cost = 8;
		else if (match_len_in_bytes >= 6)
			len_cost = 7;

		uint32_t dist_cost = 5;
		if (dist < 512)
			dist_cost += g_tdefl_small_dist_extra[dist & 511];
		else
		{
			dist_cost += g_tdefl_large_dist_extra[std::min<uint32_t>(dist, 32767) >> 8];
			while (dist >= 32768)
			{
				dist_cost++;
				dist >>= 1;
			}
		}
		return len_cost + dist_cost;
	}

	uint32_t hash_hsieh(const uint8_t* pBuf, size_t len, uint32_t salt)
	{
		if (!pBuf || !len)
			return 0;

		uint32_t h = static_cast<uint32_t>(len + (salt << 16));

		const uint32_t bytes_left = len & 3;
		len >>= 2;

		while (len--)
		{
			const uint16_t* pWords = reinterpret_cast<const uint16_t*>(pBuf);

			h += pWords[0];

			const uint32_t t = (pWords[1] << 11) ^ h;
			h = (h << 16) ^ t;

			pBuf += sizeof(uint32_t);

			h += h >> 11;
		}

		switch (bytes_left)
		{
		case 1:
			h += *reinterpret_cast<const signed char*>(pBuf);
			h ^= h << 10;
			h += h >> 1;
			break;
		case 2:
			h += *reinterpret_cast<const uint16_t*>(pBuf);
			h ^= h << 11;
			h += h >> 17;
			break;
		case 3:
			h += *reinterpret_cast<const uint16_t*>(pBuf);
			h ^= h << 16;
			h ^= (static_cast<signed char>(pBuf[sizeof(uint16_t)])) << 18;
			h += h >> 11;
			break;
		default:
			break;
		}

		h ^= h << 3;
		h += h >> 5;
		h ^= h << 4;
		h += h >> 17;
		h ^= h << 25;
		h += h >> 6;

		return h;
	}

	// BC7 entropy reduction transform with Deflate/LZMA/LZHAM optimizations
	bool reduce_entropy(uint8_t* pBlock_bytes, uint32_t num_blocks,
		uint32_t total_block_stride_in_bytes, uint32_t block_size_to_optimize_in_bytes,
		const reduce_entropy_params& params, uint32_t& total_modified,
		diff_block_func_type* pDiff_block_func, void* pDiff_block_func_user_data,
		const float* pBlock_mse_scales)
	{
		assert(total_block_stride_in_bytes && block_size_to_optimize_in_bytes);
		assert(total_block_stride_in_bytes >= block_size_to_optimize_in_bytes);

		assert((block_size_to_optimize_in_bytes >= MIN_MATCH_LEN) && (block_size_to_optimize_in_bytes <= MAX_BLOCK_SIZE_IN_BYTES));
		if ((block_size_to_optimize_in_bytes < MIN_MATCH_LEN) || (block_size_to_optimize_in_bytes > MAX_BLOCK_SIZE_IN_BYTES))
			return false;

		const int total_blocks_to_check = std::max<uint32_t>(1U, params.m_lookback_window_size / total_block_stride_in_bytes);

#if ERT_ENABLE_DEBUG
		uint32_t len_hist[MAX_BLOCK_SIZE_IN_BYTES + 1];
		uint32_t second_len_hist[MAX_BLOCK_SIZE_IN_BYTES + 1];
		uint32_t total_second_matches = 0;
		uint32_t total_smooth_blocks = 0;
#endif

		int prev_match_window_ofs_to_favor_cont = -1, prev_match_dist_to_favor = -1;

		const uint32_t HASH_SIZE = 8192;
		uint32_t hash[HASH_SIZE];
				
		for (uint32_t block_index = 0; block_index < num_blocks; block_index++)
		{
			if ((block_index & 0xFF) == 0)
				memset(hash, 0, sizeof(hash));

			uint8_t* pOrig_block = &pBlock_bytes[block_index * total_block_stride_in_bytes];

			float max_std_dev = 0.0f;
			float cur_mse = (*pDiff_block_func)(pDiff_block_func_user_data, pOrig_block, block_index, &max_std_dev);
			if (cur_mse < 0.0f)
				return false;

			if ((params.m_skip_zero_mse_blocks) && (cur_mse == 0.0f))
				continue;

			float yl = clampf(max_std_dev / params.m_max_smooth_block_std_dev, 0.0f, 1.0f);
			yl = yl * yl;
			float smooth_block_mse_scale = lerp(params.m_smooth_block_max_mse_scale, 1.0f, yl);

			if (pBlock_mse_scales)
			{
				if (pBlock_mse_scales[block_index] > 0.0f)
				{
					smooth_block_mse_scale = pBlock_mse_scales[block_index];
				}
			}
			
#if ERT_ENABLE_DEBUG
			if (smooth_block_mse_scale > 1.0f)
				total_smooth_blocks++;
#endif

			float cur_bits = (LITERAL_BITS * block_size_to_optimize_in_bytes);
			float cur_t = cur_mse * smooth_block_mse_scale + cur_bits * params.m_lambda;

			int first_block_to_check = std::max<int>(0, block_index - total_blocks_to_check);
			int last_block_to_check = block_index - 1;

			uint8_t best_block[MAX_BLOCK_SIZE_IN_BYTES];
			memcpy(best_block, pOrig_block, block_size_to_optimize_in_bytes);

			float best_t = cur_t;
			uint32_t best_match_len = 0, best_match_src_window_ofs = 0, best_match_dst_window_ofs = 0, best_match_dst_block_ofs = 0;
			float best_match_bits = 0;

			// Don't let thresh_ms_err be 0 to let zero error blocks have slightly increased distortion
			const float thresh_ms_err = params.m_max_allowed_rms_increase_ratio * params.m_max_allowed_rms_increase_ratio * std::max(cur_mse, 1.0f);
			
			for (int prev_block_index = last_block_to_check; prev_block_index >= first_block_to_check; --prev_block_index)
			{
				const uint8_t* pPrev_blk = &pBlock_bytes[prev_block_index * total_block_stride_in_bytes];

				for (uint32_t len = block_size_to_optimize_in_bytes; len >= MIN_MATCH_LEN; len--)
				{
					if (params.m_allow_relative_movement)
					{
						for (uint32_t src_ofs = 0; src_ofs <= (block_size_to_optimize_in_bytes - len); src_ofs++)
						{
							assert(len + src_ofs <= block_size_to_optimize_in_bytes);
							
							const uint32_t src_match_window_ofs = prev_block_index * total_block_stride_in_bytes + src_ofs;

							for (uint32_t dst_ofs = 0; dst_ofs <= (block_size_to_optimize_in_bytes - len); dst_ofs++)
							{
								assert(len + dst_ofs <= block_size_to_optimize_in_bytes);
								
								const uint32_t dst_match_window_ofs = block_index * total_block_stride_in_bytes + dst_ofs;

								const uint32_t match_dist = dst_match_window_ofs - src_match_window_ofs;
																
								float trial_match_bits, trial_total_bits;

								uint32_t hs = hash_hsieh(pPrev_blk + src_ofs, len, dst_ofs);

#if ERT_FAVOR_CONT_AND_REP0_MATCHES
								// Continue a previous match (which would cross block boundaries)
								if (((int)src_match_window_ofs == prev_match_window_ofs_to_favor_cont) && (dst_ofs == 0))
								{
									trial_match_bits = MATCH_CONTINUE_BITS;
									trial_total_bits = (block_size_to_optimize_in_bytes - len) * LITERAL_BITS + MATCH_CONTINUE_BITS;
								}
								// Exploit REP0 matches
								else if ((prev_match_dist_to_favor != -1) && (src_match_window_ofs == (dst_match_window_ofs - prev_match_dist_to_favor)))
								{
									trial_match_bits = MATCH_REP0_BITS;
									trial_total_bits = (block_size_to_optimize_in_bytes - len) * LITERAL_BITS + MATCH_REP0_BITS;
								}
								else
								{
									trial_match_bits = (float)compute_match_cost_estimate(match_dist, len);
									trial_total_bits = (block_size_to_optimize_in_bytes - len) * LITERAL_BITS + trial_match_bits;
										
									uint32_t hash_check = hash[hs & (HASH_SIZE - 1)];
									if ((hash_check & 0xFF) == (block_index & 0xFF))
									{
										if ((hash_check >> 8) == (hs >> 8))
											continue;
									}
								}
#else
								uint32_t hash_check = hash[hs & (HASH_SIZE - 1)];
								if ((hash_check & 0xFF) == (block_index & 0xFF))
								{
									if ((hash_check >> 8) == (hs >> 8))
										continue;
								}
#endif

								hash[hs & (HASH_SIZE - 1)] = (hs & 0xFFFFFF00) | (block_index & 0xFF);

								const float trial_total_bits_times_lambda = trial_total_bits * params.m_lambda;
								
								uint8_t trial_block[MAX_BLOCK_SIZE_IN_BYTES];
								memcpy(trial_block, pOrig_block, block_size_to_optimize_in_bytes);
								memcpy(trial_block + dst_ofs, pPrev_blk + src_ofs, len);

								float trial_mse = (*pDiff_block_func)(pDiff_block_func_user_data, trial_block, block_index, nullptr);
								if (trial_mse < 0.0f)
									continue;

								if (trial_mse < thresh_ms_err)
								{
									float t = trial_mse * smooth_block_mse_scale + trial_total_bits_times_lambda;

									if (t < best_t)
									{
										best_t = t;
										memcpy(best_block, trial_block, block_size_to_optimize_in_bytes);
										best_match_len = len;
										best_match_src_window_ofs = src_match_window_ofs;
										best_match_dst_window_ofs = dst_match_window_ofs;
										best_match_dst_block_ofs = dst_ofs;
										best_match_bits = trial_match_bits;
									}
								}

							} // dst_ofs
						} // src_ofs
					}
					else
					{
						const uint32_t match_dist = (block_index - prev_block_index) * total_block_stride_in_bytes;

						// Assume the block has 1 match and block_size_to_optimize_in_bytes-match_len literals.
						const float trial_match_bits = (float)compute_match_cost_estimate(match_dist, len);
						const float trial_total_bits = (block_size_to_optimize_in_bytes - len) * LITERAL_BITS + trial_match_bits;
						const float trial_total_bits_times_lambda = trial_total_bits * params.m_lambda;

						for (uint32_t ofs = 0; ofs <= (block_size_to_optimize_in_bytes - len); ofs++)
						{
							assert(len + ofs <= block_size_to_optimize_in_bytes);
							
							const uint32_t dst_match_window_ofs = block_index * total_block_stride_in_bytes + ofs;
							const uint32_t src_match_window_ofs = prev_block_index * total_block_stride_in_bytes + ofs;

							float trial_match_bits_to_use = trial_match_bits;
							float trial_total_bits_times_lambda_to_use = trial_total_bits_times_lambda;
														
							uint32_t hs = hash_hsieh(pPrev_blk + ofs, len, ofs);

#if ERT_FAVOR_CONT_AND_REP0_MATCHES
							// Continue a previous match (which would cross block boundaries)
							if (((int)src_match_window_ofs == prev_match_window_ofs_to_favor_cont) && (ofs == 0))
							{
								float continue_match_trial_bits = (block_size_to_optimize_in_bytes - len) * LITERAL_BITS + MATCH_CONTINUE_BITS;
								trial_match_bits_to_use = MATCH_CONTINUE_BITS;
								trial_total_bits_times_lambda_to_use = continue_match_trial_bits * params.m_lambda;
							}
							// Exploit REP0 matches
							else if ((prev_match_dist_to_favor != -1) && (src_match_window_ofs == (dst_match_window_ofs - prev_match_dist_to_favor)))
							{
								float continue_match_trial_bits = (block_size_to_optimize_in_bytes - len) * LITERAL_BITS + MATCH_REP0_BITS;
								trial_match_bits_to_use = MATCH_REP0_BITS;
								trial_total_bits_times_lambda_to_use = continue_match_trial_bits * params.m_lambda;
							}
							else
							{
								uint32_t hash_check = hash[hs & (HASH_SIZE - 1)];
								if ((hash_check & 0xFF) == (block_index & 0xFF))
								{
									if ((hash_check >> 8) == (hs >> 8))
										continue;
								}
							}
#else
							uint32_t hash_check = hash[hs & (HASH_SIZE - 1)];
							if ((hash_check & 0xFF) == (block_index & 0xFF))
							{
								if ((hash_check >> 8) == (hs >> 8))
									continue;
							}
#endif

							hash[hs & (HASH_SIZE - 1)] = (hs & 0xFFFFFF00) | (block_index & 0xFF);

							uint8_t trial_block[MAX_BLOCK_SIZE_IN_BYTES];
							memcpy(trial_block, pOrig_block, block_size_to_optimize_in_bytes);
							memcpy(trial_block + ofs, pPrev_blk + ofs, len);

							float trial_mse = (*pDiff_block_func)(pDiff_block_func_user_data, trial_block, block_index, nullptr);
							if (trial_mse < 0.0f)
								continue;

							if (trial_mse < thresh_ms_err)
							{
								float t = trial_mse * smooth_block_mse_scale + trial_total_bits_times_lambda_to_use;
								
								if (t < best_t)
								{
									best_t = t;
									memcpy(best_block, trial_block, block_size_to_optimize_in_bytes);
									best_match_len = len;
									best_match_src_window_ofs = src_match_window_ofs;
									best_match_dst_window_ofs = dst_match_window_ofs;
									best_match_dst_block_ofs = ofs;
									best_match_bits = trial_match_bits_to_use;
								}
							}
						} // ofs
					}

				} // len

			} // prev_block_index

			if (best_t < cur_t)
			{
				uint32_t best_second_match_len = 0, best_second_match_src_window_ofs = 0, best_second_match_dst_window_ofs = 0, best_second_match_dst_block_ofs = 0;
								
				// Try injecting a second match, being sure it does't overlap with the first.
				if ((params.m_try_two_matches) && (best_match_len <= (block_size_to_optimize_in_bytes - 3)))
				{
					uint8_t matched_flags[MAX_BLOCK_SIZE_IN_BYTES]{};
					memset(matched_flags + best_match_dst_block_ofs, 1, best_match_len);

					uint8_t orig_best_block[MAX_BLOCK_SIZE_IN_BYTES];
					memcpy(orig_best_block, best_block, block_size_to_optimize_in_bytes);
										
					for (int prev_block_index = last_block_to_check; prev_block_index >= first_block_to_check; --prev_block_index)
					{
						const uint8_t* pPrev_blk = &pBlock_bytes[prev_block_index * total_block_stride_in_bytes];

						const uint32_t match_dist = (block_index - prev_block_index) * total_block_stride_in_bytes;

						for (uint32_t len = 3; len <= (block_size_to_optimize_in_bytes - best_match_len); len++)
						{
							const float trial_total_bits = (block_size_to_optimize_in_bytes - len - best_match_len) * LITERAL_BITS + compute_match_cost_estimate(match_dist, len) + best_match_bits;

							const float trial_total_bits_times_lambda = trial_total_bits * params.m_lambda;

							for (uint32_t ofs = 0; ofs <= (block_size_to_optimize_in_bytes - len); ofs++)
							{
								int i;
								for (i = 0; i < (int)len; i++)
									if (matched_flags[ofs + i])
										break;
								if (i != (int)len)
									continue;

								assert(len + ofs <= block_size_to_optimize_in_bytes);

								const uint32_t dst_match_window_ofs = block_index * total_block_stride_in_bytes + ofs;
								const uint32_t src_match_window_ofs = prev_block_index * total_block_stride_in_bytes + ofs;

								uint8_t trial_block[MAX_BLOCK_SIZE_IN_BYTES];
								memcpy(trial_block, orig_best_block, block_size_to_optimize_in_bytes);
								memcpy(trial_block + ofs, pPrev_blk + ofs, len);

								float trial_mse = (*pDiff_block_func)(pDiff_block_func_user_data, trial_block, block_index, nullptr);
								if (trial_mse < 0.0f)
									continue;

								if (trial_mse < thresh_ms_err)
								{
									float t = trial_mse * smooth_block_mse_scale + trial_total_bits_times_lambda;

									if (t < best_t)
									{
										best_t = t;
										memcpy(best_block, trial_block, block_size_to_optimize_in_bytes);
										best_second_match_len = len;
										best_second_match_src_window_ofs = src_match_window_ofs;
										best_second_match_dst_window_ofs = dst_match_window_ofs;
										best_second_match_dst_block_ofs = ofs;
									}
								}
							}
						}
					}
				}

				memcpy(pOrig_block, best_block, block_size_to_optimize_in_bytes);
				total_modified++;

				if ((best_second_match_len == 0) || (best_match_dst_window_ofs > best_second_match_dst_window_ofs))
				{
					int best_match_dist = best_match_dst_window_ofs - best_match_src_window_ofs;
					assert(best_match_dist >= 1);
					(void)best_match_dist;

					if (block_size_to_optimize_in_bytes == total_block_stride_in_bytes)
					{
						// If the match goes all the way to the end of a block, we can try to continue it on the next encoded block.
						if ((best_match_dst_block_ofs + best_match_len) == total_block_stride_in_bytes)
							prev_match_window_ofs_to_favor_cont = best_match_src_window_ofs + best_match_len;
						else
							prev_match_window_ofs_to_favor_cont = -1;
					}

#if ERT_FAVOR_REP0_MATCHES
					// Compute the window offset where a cheaper REP0 match would be available
					prev_match_dist_to_favor = best_match_dist;
#endif
				}
				else
				{
					int best_match_dist = best_second_match_dst_window_ofs - best_second_match_src_window_ofs;
					assert(best_match_dist >= 1);
					(void)best_match_dist;

					if (block_size_to_optimize_in_bytes == total_block_stride_in_bytes)
					{
						// If the match goes all the way to the end of a block, we can try to continue it on the next encoded block.
						if ((best_second_match_dst_block_ofs + best_second_match_len) == total_block_stride_in_bytes)
							prev_match_window_ofs_to_favor_cont = best_second_match_src_window_ofs + best_second_match_len;
						else
							prev_match_window_ofs_to_favor_cont = -1;
					}

#if ERT_FAVOR_REP0_MATCHES
					// Compute the window offset where a cheaper REP0 match would be available
					prev_match_dist_to_favor = best_match_dist;
#endif
				}

#if ERT_ENABLE_DEBUG
				len_hist[best_match_len]++;

				if (best_second_match_len)
				{
					second_len_hist[best_second_match_len]++;
					total_second_matches++;
				}
#endif
			}
			else
			{
				prev_match_window_ofs_to_favor_cont = -1;
			}
						
		} // block_index
#if ERT_ENABLE_DEBUG
		if (params.m_debug_output)
		{
			printf("Total smooth blocks: %3.2f%%\n", total_smooth_blocks * 100.0f / num_blocks);

			printf("Match length histogram:\n");
			for (uint32_t i = MIN_MATCH_LEN; i <= block_size_to_optimize_in_bytes; i++)
				printf("%u%c", len_hist[i], (i < block_size_to_optimize_in_bytes) ? ',' : '\n');

			printf("Total second matches: %u %3.2f%%\n", total_second_matches, total_second_matches * 100.0f / num_blocks);
			printf("Secod match length histogram:\n");
			for (uint32_t i = MIN_MATCH_LEN; i <= block_size_to_optimize_in_bytes; i++)
				printf("%u%c", second_len_hist[i], (i < block_size_to_optimize_in_bytes) ? ',' : '\n');
		}
#endif
		return true;
	}

} // namespace ert
