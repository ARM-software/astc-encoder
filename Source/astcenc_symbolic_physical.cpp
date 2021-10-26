// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2011-2021 Arm Limited
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
 * @brief Functions for converting between symbolic and physical encodings.
 */

#include "astcenc_internal.h"

#include <cassert>

/**
 * @brief Write up to 8 bits at an arbitrary bit offset.
 *
 * The stored value is at most 8 bits, but can be stored at an offset of between 0 and 7 bits so
 * may span two separate bytes in memory.
 *
 * @param         value       The value to write.
 * @param         bitcount    The number of bits to write, starting from LSB.
 * @param         bitoffset   The bit offset to store at, between 0 and 7.
 * @param[in,out] ptr         The data pointer to write to.
 */
static inline void write_bits(
	int value,
	int bitcount,
	int bitoffset,
	uint8_t* ptr
) {
	int mask = (1 << bitcount) - 1;
	value &= mask;
	ptr += bitoffset >> 3;
	bitoffset &= 7;
	value <<= bitoffset;
	mask <<= bitoffset;
	mask = ~mask;

	ptr[0] &= mask;
	ptr[0] |= value;
	ptr[1] &= mask >> 8;
	ptr[1] |= value >> 8;
}

/**
 * @brief Read up to 8 bits at an arbitrary bit offset.
 *
 * The stored value is at most 8 bits, but can be stored at an offset of between 0 and 7 bits so may
 * span two separate bytes in memory.
 *
 * @param         bitcount    The number of bits to read.
 * @param         bitoffset   The bit offset to read from, between 0 and 7.
 * @param[in,out] ptr         The data pointer to read from.
 *
 * @return The read value.
 */
static inline int read_bits(
	int bitcount,
	int bitoffset,
	const uint8_t* ptr
) {
	int mask = (1 << bitcount) - 1;
	ptr += bitoffset >> 3;
	bitoffset &= 7;
	int value = ptr[0] | (ptr[1] << 8);
	value >>= bitoffset;
	value &= mask;
	return value;
}

/**
 * @brief Reverse bits in a byte.
 *
 * @param p   The value to reverse.
  *
 * @return The reversed result.
 */
static inline int bitrev8(int p)
{
	p = ((p & 0x0F) << 4) | ((p >> 4) & 0x0F);
	p = ((p & 0x33) << 2) | ((p >> 2) & 0x33);
	p = ((p & 0x55) << 1) | ((p >> 1) & 0x55);
	return p;
}

/* See header for documentation. */
void symbolic_to_physical(
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& scb,
	physical_compressed_block& pcb
) {
	assert(scb.block_type != SYM_BTYPE_ERROR);

	// Constant color block using UNORM16 colors
	if (scb.block_type == SYM_BTYPE_CONST_U16)
	{
		// There is currently no attempt to coalesce larger void-extents
		static const uint8_t cbytes[8] { 0xFC, 0xFD, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF };
		for (unsigned int i = 0; i < 8; i++)
		{
			pcb.data[i] = cbytes[i];
		}

		for (unsigned int i = 0; i < BLOCK_MAX_COMPONENTS; i++)
		{
			pcb.data[2 * i + 8] = scb.constant_color[i] & 0xFF;
			pcb.data[2 * i + 9] = (scb.constant_color[i] >> 8) & 0xFF;
		}

		return;
	}

	// Constant color block using FP16 colors
	if (scb.block_type == SYM_BTYPE_CONST_F16)
	{
		// There is currently no attempt to coalesce larger void-extents
		static const uint8_t cbytes[8]  { 0xFC, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF };
		for (unsigned int i = 0; i < 8; i++)
		{
			pcb.data[i] = cbytes[i];
		}

		for (unsigned int i = 0; i < BLOCK_MAX_COMPONENTS; i++)
		{
			pcb.data[2 * i + 8] = scb.constant_color[i] & 0xFF;
			pcb.data[2 * i + 9] = (scb.constant_color[i] >> 8) & 0xFF;
		}

		return;
	}

	unsigned int partition_count = scb.partition_count;

	// Compress the weights.
	// They are encoded as an ordinary integer-sequence, then bit-reversed
	uint8_t weightbuf[16] { 0 };

	const auto& bm = bsd.get_block_mode(scb.block_mode);
	const auto& di = bsd.get_decimation_info(bm.decimation_mode);
	int weight_count = di.weight_count;
	quant_method weight_quant_method = bm.get_weight_quant_mode();
	int is_dual_plane = bm.is_dual_plane;

	int real_weight_count = is_dual_plane ? 2 * weight_count : weight_count;

	int bits_for_weights = get_ise_sequence_bitcount(real_weight_count, weight_quant_method);

	if (is_dual_plane)
	{
		uint8_t weights[64];
		for (int i = 0; i < weight_count; i++)
		{
			weights[2 * i] = scb.weights[i];
			weights[2 * i + 1] = scb.weights[i + WEIGHTS_PLANE2_OFFSET];
		}
		encode_ise(weight_quant_method, real_weight_count, weights, weightbuf, 0);
	}
	else
	{
		encode_ise(weight_quant_method, weight_count, scb.weights, weightbuf, 0);
	}

	for (int i = 0; i < 16; i++)
	{
		pcb.data[i] = static_cast<uint8_t>(bitrev8(weightbuf[15 - i]));
	}

	write_bits(scb.block_mode, 11, 0, pcb.data);
	write_bits(partition_count - 1, 2, 11, pcb.data);

	int below_weights_pos = 128 - bits_for_weights;

	// Encode partition index and color endpoint types for blocks with 2+ partitions
	if (partition_count > 1)
	{
		write_bits(scb.partition_index, 6, 13, pcb.data);
		write_bits(scb.partition_index >> 6, PARTITION_INDEX_BITS - 6, 19, pcb.data);

		if (scb.color_formats_matched)
		{
			write_bits(scb.color_formats[0] << 2, 6, 13 + PARTITION_INDEX_BITS, pcb.data);
		}
		else
		{
			// Check endpoint types for each partition to determine the lowest class present
			int low_class = 4;

			for (unsigned int i = 0; i < partition_count; i++)
			{
				int class_of_format = scb.color_formats[i] >> 2;
				low_class = astc::min(class_of_format, low_class);
			}

			if (low_class == 3)
			{
				low_class = 2;
			}

			int encoded_type = low_class + 1;
			int bitpos = 2;

			for (unsigned int i = 0; i < partition_count; i++)
			{
				int classbit_of_format = (scb.color_formats[i] >> 2) - low_class;
				encoded_type |= classbit_of_format << bitpos;
				bitpos++;
			}

			for (unsigned int i = 0; i < partition_count; i++)
			{
				int lowbits_of_format = scb.color_formats[i] & 3;
				encoded_type |= lowbits_of_format << bitpos;
				bitpos += 2;
			}

			int encoded_type_lowpart = encoded_type & 0x3F;
			int encoded_type_highpart = encoded_type >> 6;
			int encoded_type_highpart_size = (3 * partition_count) - 4;
			int encoded_type_highpart_pos = 128 - bits_for_weights - encoded_type_highpart_size;
			write_bits(encoded_type_lowpart, 6, 13 + PARTITION_INDEX_BITS, pcb.data);
			write_bits(encoded_type_highpart, encoded_type_highpart_size, encoded_type_highpart_pos, pcb.data);
			below_weights_pos -= encoded_type_highpart_size;
		}
	}
	else
	{
		write_bits(scb.color_formats[0], 4, 13, pcb.data);
	}

	// In dual-plane mode, encode the color component of the second plane of weights
	if (is_dual_plane)
	{
		write_bits(scb.plane2_component, 2, below_weights_pos - 2, pcb.data);
	}

	// Encode the color components
	uint8_t values_to_encode[32];
	int valuecount_to_encode = 0;
	for (unsigned int i = 0; i < scb.partition_count; i++)
	{
		int vals = 2 * (scb.color_formats[i] >> 2) + 2;
		assert(vals <= 8);
		for (int j = 0; j < vals; j++)
		{
			values_to_encode[j + valuecount_to_encode] = scb.color_values[i][j];
		}
		valuecount_to_encode += vals;
	}

	encode_ise(scb.get_color_quant_mode(), valuecount_to_encode, values_to_encode, pcb.data,
	           scb.partition_count == 1 ? 17 : 19 + PARTITION_INDEX_BITS);
}

/* See header for documentation. */
void physical_to_symbolic(
	const block_size_descriptor& bsd,
	const physical_compressed_block& pcb,
	symbolic_compressed_block& scb
) {
	uint8_t bswapped[16];

	scb.block_type = SYM_BTYPE_NONCONST;

	// Extract header fields
	int block_mode = read_bits(11, 0, pcb.data);
	if ((block_mode & 0x1FF) == 0x1FC)
	{
		// Constant color block

		// Check what format the data has
		if (block_mode & 0x200)
		{
			scb.block_type = SYM_BTYPE_CONST_F16;
		}
		else
		{
			scb.block_type = SYM_BTYPE_CONST_U16;
		}

		scb.partition_count = 0;
		for (int i = 0; i < 4; i++)
		{
			scb.constant_color[i] = pcb.data[2 * i + 8] | (pcb.data[2 * i + 9] << 8);
		}

		// Additionally, check that the void-extent
		if (bsd.zdim == 1)
		{
			// 2D void-extent
			int rsvbits = read_bits(2, 10, pcb.data);
			if (rsvbits != 3)
			{
				scb.block_type = SYM_BTYPE_ERROR;
				return;
			}

			int vx_low_s = read_bits(8, 12, pcb.data) | (read_bits(5, 12 + 8, pcb.data) << 8);
			int vx_high_s = read_bits(8, 25, pcb.data) | (read_bits(5, 25 + 8, pcb.data) << 8);
			int vx_low_t = read_bits(8, 38, pcb.data) | (read_bits(5, 38 + 8, pcb.data) << 8);
			int vx_high_t = read_bits(8, 51, pcb.data) | (read_bits(5, 51 + 8, pcb.data) << 8);

			int all_ones = vx_low_s == 0x1FFF && vx_high_s == 0x1FFF && vx_low_t == 0x1FFF && vx_high_t == 0x1FFF;

			if ((vx_low_s >= vx_high_s || vx_low_t >= vx_high_t) && !all_ones)
			{
				scb.block_type = SYM_BTYPE_ERROR;
				return;
			}
		}
		else
		{
			// 3D void-extent
			int vx_low_s = read_bits(9, 10, pcb.data);
			int vx_high_s = read_bits(9, 19, pcb.data);
			int vx_low_t = read_bits(9, 28, pcb.data);
			int vx_high_t = read_bits(9, 37, pcb.data);
			int vx_low_p = read_bits(9, 46, pcb.data);
			int vx_high_p = read_bits(9, 55, pcb.data);

			int all_ones = vx_low_s == 0x1FF && vx_high_s == 0x1FF && vx_low_t == 0x1FF && vx_high_t == 0x1FF && vx_low_p == 0x1FF && vx_high_p == 0x1FF;

			if ((vx_low_s >= vx_high_s || vx_low_t >= vx_high_t || vx_low_p >= vx_high_p) && !all_ones)
			{
				scb.block_type = SYM_BTYPE_ERROR;
				return;
			}
		}

		return;
	}

	unsigned int packed_index = bsd.block_mode_packed_index[block_mode];
	if (packed_index == BLOCK_BAD_BLOCK_MODE)
	{
		scb.block_type = SYM_BTYPE_ERROR;
		return;
	}

	const auto& bm = bsd.get_block_mode(block_mode);
	const auto& di = bsd.get_decimation_info(bm.decimation_mode);

	int weight_count = di.weight_count;
	quant_method weight_quant_method = (quant_method)bm.quant_mode;
	int is_dual_plane = bm.is_dual_plane;

	int real_weight_count = is_dual_plane ? 2 * weight_count : weight_count;

	int partition_count = read_bits(2, 11, pcb.data) + 1;

	scb.block_mode = static_cast<uint16_t>(block_mode);
	scb.partition_count = static_cast<uint8_t>(partition_count);

	for (int i = 0; i < 16; i++)
	{
		bswapped[i] = static_cast<uint8_t>(bitrev8(pcb.data[15 - i]));
	}

	int bits_for_weights = get_ise_sequence_bitcount(real_weight_count, weight_quant_method);

	int below_weights_pos = 128 - bits_for_weights;

	if (is_dual_plane)
	{
		uint8_t indices[64];
		decode_ise(weight_quant_method, real_weight_count, bswapped, indices, 0);
		for (int i = 0; i < weight_count; i++)
		{
			scb.weights[i] = indices[2 * i];
			scb.weights[i + WEIGHTS_PLANE2_OFFSET] = indices[2 * i + 1];
		}
	}
	else
	{
		decode_ise(weight_quant_method, weight_count, bswapped, scb.weights, 0);
	}

	if (is_dual_plane && partition_count == 4)
	{
		scb.block_type = SYM_BTYPE_ERROR;
		return;
	}

	scb.color_formats_matched = 0;

	// Determine the format of each endpoint pair
	int color_formats[BLOCK_MAX_PARTITIONS];
	int encoded_type_highpart_size = 0;
	if (partition_count == 1)
	{
		color_formats[0] = read_bits(4, 13, pcb.data);
		scb.partition_index = 0;
	}
	else
	{
		encoded_type_highpart_size = (3 * partition_count) - 4;
		below_weights_pos -= encoded_type_highpart_size;
		int encoded_type = read_bits(6, 13 + PARTITION_INDEX_BITS, pcb.data) | (read_bits(encoded_type_highpart_size, below_weights_pos, pcb.data) << 6);
		int baseclass = encoded_type & 0x3;
		if (baseclass == 0)
		{
			for (int i = 0; i < partition_count; i++)
			{
				color_formats[i] = (encoded_type >> 2) & 0xF;
			}

			below_weights_pos += encoded_type_highpart_size;
			scb.color_formats_matched = 1;
			encoded_type_highpart_size = 0;
		}
		else
		{
			int bitpos = 2;
			baseclass--;

			for (int i = 0; i < partition_count; i++)
			{
				color_formats[i] = (((encoded_type >> bitpos) & 1) + baseclass) << 2;
				bitpos++;
			}

			for (int i = 0; i < partition_count; i++)
			{
				color_formats[i] |= (encoded_type >> bitpos) & 3;
				bitpos += 2;
			}
		}
		scb.partition_index = static_cast<uint16_t>(read_bits(6, 13, pcb.data) | (read_bits(PARTITION_INDEX_BITS - 6, 19, pcb.data) << 6));
	}

	for (int i = 0; i < partition_count; i++)
	{
		scb.color_formats[i] = static_cast<uint8_t>(color_formats[i]);
	}

	// Determine number of color endpoint integers
	int color_integer_count = 0;
	for (int i = 0; i < partition_count; i++)
	{
		int endpoint_class = color_formats[i] >> 2;
		color_integer_count += (endpoint_class + 1) * 2;
	}

	if (color_integer_count > 18)
	{
		scb.block_type = SYM_BTYPE_ERROR;
		return;
	}

	// Determine the color endpoint format to use
	static const int color_bits_arr[5] { -1, 115 - 4, 113 - 4 - PARTITION_INDEX_BITS, 113 - 4 - PARTITION_INDEX_BITS, 113 - 4 - PARTITION_INDEX_BITS };
	int color_bits = color_bits_arr[partition_count] - bits_for_weights - encoded_type_highpart_size;
	if (is_dual_plane)
	{
		color_bits -= 2;
	}

	if (color_bits < 0)
	{
		color_bits = 0;
	}

	int color_quant_level = quant_mode_table[color_integer_count >> 1][color_bits];
	if (color_quant_level < QUANT_6)
	{
		scb.block_type = SYM_BTYPE_ERROR;
		return;
	}

	// Unpack the integer color values and assign to endpoints
	scb.quant_mode = (quant_method)color_quant_level;
	uint8_t values_to_decode[32];
	decode_ise((quant_method)color_quant_level, color_integer_count, pcb.data, values_to_decode, (partition_count == 1 ? 17 : 19 + PARTITION_INDEX_BITS));

	int valuecount_to_decode = 0;
	for (int i = 0; i < partition_count; i++)
	{
		int vals = 2 * (color_formats[i] >> 2) + 2;
		for (int j = 0; j < vals; j++)
		{
			scb.color_values[i][j] = values_to_decode[j + valuecount_to_decode];
		}
		valuecount_to_decode += vals;
	}

	// Fetch component for second-plane in the case of dual plane of weights.
	if (is_dual_plane)
	{
		scb.plane2_component = static_cast<int8_t>(read_bits(2, below_weights_pos - 2, pcb.data));
	}
}
