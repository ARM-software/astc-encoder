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
 * @brief RDO payload decoding using the encoder's original symbolic layout.
 */

#if !defined(ASTCENC_DECOMPRESS_ONLY)

#include <algorithm>
#include <cstring>

#include "astcenc_internal.h"
#include "astcenc_rdo.h"

/**
 * @brief Mark a contiguous range of physical bits.
 *
 * @param[in,out] mask    The physical bit mask.
 * @param         first   The first bit to mark.
 * @param         count   The number of bits to mark.
 */
static void rdo_mark_bits(
	uint8_t mask[16],
	unsigned int first,
	unsigned int count
) {
	for (unsigned int bit = first; bit < first + count; bit++)
	{
		mask[bit / 8] |= static_cast<uint8_t>(1u << (bit % 8));
	}
}

/* See header for documentation. */
void rdo_get_structure_mask(
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& block,
	uint8_t mask[16]
) {
	std::fill_n(mask, 16, 0);
	if (block.block_type != SYM_BTYPE_NONCONST)
	{
		// Constant blocks are retained unchanged in the byte history
		std::fill_n(mask, 16, 0xFF);
		return;
	}

	// Layout fields surround both payloads; protecting only the low header would miss CEM bits
	unsigned int header_bits = 17;
	if (block.partition_count > 1)
	{
		header_bits = 19 + PARTITION_INDEX_BITS;
	}

	rdo_mark_bits(mask, 0, header_bits);
	const auto& mode = bsd.get_block_mode(block.block_mode);
	unsigned int below_weights = 128 - mode.weight_bits;
	if ((block.partition_count > 1) && !block.color_formats_matched)
	{
		unsigned int format_bits = 3 * block.partition_count - 4;
		below_weights -= format_bits;
		rdo_mark_bits(mask, below_weights, format_bits);
	}

	if (mode.is_dual_plane)
	{
		rdo_mark_bits(mask, below_weights - 2, 2);
	}
}

/* See header for documentation. */
rdo_block_decode_state rdo_prepare_block_decode(
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& block
) {
	assert(block.block_type == SYM_BTYPE_NONCONST);
	const auto& mode = bsd.get_block_mode(block.block_mode);
	const auto& decimation = bsd.get_decimation_info(mode.decimation_mode);

	rdo_block_decode_state state {
		.weight_decode = rdo_get_ise_decoder(mode.quant_mode),
		.color_decode = rdo_get_ise_decoder(block.quant_mode),
		.weight_count = decimation.weight_count,
		.real_weight_count = decimation.weight_count,
		.weight_quant_mode = mode.quant_mode,
		.color_count = 0,
		.color_offset = 17,
		.color_counts = {},
		.weight_byte_mask = 0,
		.color_byte_mask = 0,
		.is_dual_plane = mode.is_dual_plane != 0
	};

	if (state.is_dual_plane)
	{
		state.real_weight_count *= 2;
	}

	if (block.partition_count > 1)
	{
		state.color_offset = 19 + PARTITION_INDEX_BITS;
	}

	for (unsigned int i = 0; i < block.partition_count; i++)
	{
		state.color_counts[i] = 2 * (block.color_formats[i] >> 2) + 2;
		state.color_count += state.color_counts[i];
	}

	// Byte masks may overlap at shared payload boundaries, conservatively decoding both streams
	unsigned int first_weight_byte = (128 - mode.weight_bits) / 8;
	state.weight_byte_mask = 0xFFFFu & (~0u << first_weight_byte);

	unsigned int color_end = state.color_offset
	                       + get_ise_sequence_bitcount(state.color_count, block.quant_mode);

	state.color_byte_mask = ((1u << ((color_end + 7) / 8)) - 1)
	                      & (~0u << (state.color_offset / 8));
	return state;
}

/**
 * @brief Decode and unquantize the fixed-layout weight sequence.
 *
 * @param[in]     state      The fixed nonconstant layout.
 * @param[in]     physical   The physical block bytes.
 * @param[in,out] block      The symbolic block receiving the decoded values.
 */
static void rdo_decode_weight_payload(
	const rdo_block_decode_state& state,
	const uint8_t physical[16],
	symbolic_compressed_block& block
) {
	// Weight ISE runs backward from bit 127; layout bits beyond the sequence are ignored
	uint8_t reversed[16];
	for (unsigned int i = 0; i < 16; i++)
	{
		unsigned int value = physical[15 - i];
		value = ((value & 0x0F) << 4) | ((value >> 4) & 0x0F);
		value = ((value & 0x33) << 2) | ((value >> 2) & 0x33);
		value = ((value & 0x55) << 1) | ((value >> 1) & 0x55);
		reversed[i] = static_cast<uint8_t>(value);
	}

	if (!state.is_dual_plane
	    && ((state.weight_quant_mode == QUANT_3) || (state.weight_quant_mode == QUANT_5)))
	{
		// Pure trit and quint sequences can combine ISE decoding with weight unquantization
		rdo_decode_ise_radix_weights(state.weight_quant_mode, state.weight_count, reversed, block.weights);
		return;
	}

	// Binary quantizers can expand directly into mapped weight bytes
	static constexpr uint8_t binary_bits[12] { 1, 0, 2, 0, 0, 3, 0, 0, 4, 0, 0, 5 };
	unsigned int bits = binary_bits[state.weight_quant_mode];
	if (!state.is_dual_plane && (bits != 0))
	{
		const auto& transfer = quant_and_xfer_tables[state.weight_quant_mode];
		vtable4_16x8 low_table;
		vtable_prepare(low_table, transfer.unscramble_and_unquant_map);
		vtable4_32x8 table;
		vtable_prepare(table, transfer.unscramble_and_unquant_map);

		for (unsigned int i = 0; i < state.weight_count; i += 8)
		{
			unsigned int bit_offset = i * bits;
			unsigned int byte_offset = std::min(bit_offset >> 3, 8u);

			uint64_t packed = load_le64(reversed + byte_offset);
			packed >>= bit_offset - byte_offset * 8;

			vint4 indices = unpack_byte_fields(packed, bits);

			vint4 weights;
			if (bits == 5)
			{
				weights = vtable_lookup_8bit(table, indices);
			}
			else
			{
				weights = vtable_lookup_8bit(low_table, indices);
			}

			// The rounded batch remains inside the 64 stored weight slots
			store_low_8bytes(weights, block.weights + i);
		}

		return;
	}

	// Mixed-radix and dual-plane layouts share the general ISE decoder
	uint8_t indices[80];
	auto weight_quant = static_cast<quant_method>(state.weight_quant_mode);
	state.weight_decode(state.real_weight_count, reversed, indices, 0);
	const auto& transfer = quant_and_xfer_tables[weight_quant];
	if (state.is_dual_plane)
	{
		// Padding keeps the final eight-pair batch inside allocated, unused weight slots
		static_assert(WEIGHTS_PLANE2_OFFSET == 32);
		std::memset(indices + state.real_weight_count, 0, 16);
		vtable4_32x8 table;
		vtable_prepare(table, transfer.unscramble_and_unquant_map);
		for (unsigned int i = 0; i < state.weight_count; i += 8)
		{
			vint4 values = vint4::load(indices + 2 * i);
			vint4 weights = vtable_lookup_8bit(table, values);
			store_deinterleaved_bytes(weights, block.weights + i, block.weights + WEIGHTS_PLANE2_OFFSET + i);
		}
	}
	else
	{
		vtable_32x8 unquant_table;
		vtable_prepare(unquant_table, transfer.unscramble_and_unquant_map);
		unsigned int i = 0;
		for (; i + ASTCENC_SIMD_WIDTH <= state.weight_count; i += ASTCENC_SIMD_WIDTH)
		{
			vint weights = vtable_lookup_32bit(unquant_table, vint(indices + i));
			pack_and_store_low_bytes(weights, block.weights + i);
		}

		for (; i < state.weight_count; i++)
		{
			block.weights[i] = transfer.unscramble_and_unquant_map[indices[i]];
		}
	}
}

/**
 * @brief Decode and unquantize the fixed-layout endpoint sequence.
 *
 * @param[in]     state      The fixed nonconstant layout.
 * @param[in]     physical   The physical block bytes.
 * @param[in,out] block      The symbolic block receiving the decoded values.
 */
static void rdo_decode_color_payload(
	const rdo_block_decode_state& state,
	const uint8_t physical[16],
	symbolic_compressed_block& block
) {
	uint8_t values[32];
	state.color_decode(state.color_count, physical, values, state.color_offset);
	const auto* unpack = color_scrambled_pquant_to_uquant_tables[block.quant_mode - QUANT_6];
	unsigned int offset = 0;
	for (unsigned int i = 0; i < block.partition_count; i++)
	{
		for (unsigned int j = 0; j < state.color_counts[i]; j++)
		{
			block.color_values[i][j] = unpack[values[offset + j]];
		}

		offset += state.color_counts[i];
	}
}

/* See header for documentation. */
void rdo_decode_changed_payload(
	const rdo_block_decode_state& state,
	const uint8_t physical[16],
	unsigned int byte_mask,
	symbolic_compressed_block& block
) {
	if (byte_mask & state.weight_byte_mask)
	{
		rdo_decode_weight_payload(state, physical, block);
	}

	if (byte_mask & state.color_byte_mask)
	{
		rdo_decode_color_payload(state, physical, block);
	}
}

#endif
