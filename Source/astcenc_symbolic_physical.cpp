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
 * @brief Functions for converting between symbolic and physical encodings.
 */

#include "astcenc_internal.h"

#include <assert.h>

// routine to write up to 8 bits
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

// routine to read up to 8 bits
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

static inline int bitrev8(int p)
{
	p = ((p & 0xF) << 4) | ((p >> 4) & 0xF);
	p = ((p & 0x33) << 2) | ((p >> 2) & 0x33);
	p = ((p & 0x55) << 1) | ((p >> 1) & 0x55);
	return p;
}

void symbolic_to_physical(
	const block_size_descriptor& bsd,
	const symbolic_compressed_block& scb,
	physical_compressed_block& pcb
) {
	if (scb.block_mode == -2)
	{
		// UNORM16 constant-color block.
		// This encodes separate constant-color blocks. There is currently
		// no attempt to coalesce them into larger void-extents.

		static const uint8_t cbytes[8] = { 0xFC, 0xFD, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF };
		for (int i = 0; i < 8; i++)
		{
			pcb.data[i] = cbytes[i];
		}

		for (int i = 0; i < 4; i++)
		{
			pcb.data[2 * i + 8] = scb.constant_color[i] & 0xFF;
			pcb.data[2 * i + 9] = (scb.constant_color[i] >> 8) & 0xFF;
		}

		return;
	}

	if (scb.block_mode == -1)
	{
		// FP16 constant-color block.
		// This encodes separate constant-color blocks. There is currently
		// no attempt to coalesce them into larger void-extents.

		static const uint8_t cbytes[8] = { 0xFC, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF };
		for (int i = 0; i < 8; i++)
		{
			pcb.data[i] = cbytes[i];
		}

		for (int i = 0; i < 4; i++)
		{
			pcb.data[2 * i + 8] = scb.constant_color[i] & 0xFF;
			pcb.data[2 * i + 9] = (scb.constant_color[i] >> 8) & 0xFF;
		}

		return;
	}

	int partition_count = scb.partition_count;

	// first, compress the weights. They are encoded as an ordinary
	// integer-sequence, then bit-reversed
	uint8_t weightbuf[16];
	for (int i = 0; i < 16; i++)
	{
		weightbuf[i] = 0;
	}

	const decimation_table *const *ixtab2 = bsd.decimation_tables;
	
	const int packed_index = bsd.block_mode_to_packed[scb.block_mode];
	assert(packed_index >= 0 && packed_index < bsd.block_mode_packed_count);
	const block_mode& bm = bsd.block_modes_packed[packed_index];

	int weight_count = ixtab2[bm.decimation_mode]->num_weights;
	int weight_quantization_method = bm.quantization_mode;
	int is_dual_plane = bm.is_dual_plane;

	int real_weight_count = is_dual_plane ? 2 * weight_count : weight_count;

	int bits_for_weights = compute_ise_bitcount(real_weight_count,
	                                            (quantization_method) weight_quantization_method);

	if (is_dual_plane)
	{
		uint8_t weights[64];
		for (int i = 0; i < weight_count; i++)
		{
			weights[2 * i] = scb.plane1_weights[i];
			weights[2 * i + 1] = scb.plane2_weights[i];
		}
		encode_ise(weight_quantization_method, real_weight_count, weights, weightbuf, 0);
	}
	else
	{
		encode_ise(weight_quantization_method, weight_count, scb.plane1_weights, weightbuf, 0);
	}

	for (int i = 0; i < 16; i++)
	{
		pcb.data[i] = bitrev8(weightbuf[15 - i]);
	}

	write_bits(scb.block_mode, 11, 0, pcb.data);
	write_bits(partition_count - 1, 2, 11, pcb.data);

	int below_weights_pos = 128 - bits_for_weights;

	// encode partition index and color endpoint types for blocks with
	// 2 or more partitions.
	if (partition_count > 1)
	{
		write_bits(scb.partition_index, 6, 13, pcb.data);
		write_bits(scb.partition_index >> 6, PARTITION_BITS - 6, 19, pcb.data);

		if (scb.color_formats_matched)
		{
			write_bits(scb.color_formats[0] << 2, 6, 13 + PARTITION_BITS, pcb.data);
		}
		else
		{
			// go through the selected endpoint type classes for each partition
			// in order to determine the lowest class present.
			int low_class = 4;

			for (int i = 0; i < partition_count; i++)
			{
				int class_of_format = scb.color_formats[i] >> 2;
				if (class_of_format < low_class)
				{
					low_class = class_of_format;
				}
			}

			if (low_class == 3)
			{
				low_class = 2;
			}

			int encoded_type = low_class + 1;
			int bitpos = 2;

			for (int i = 0; i < partition_count; i++)
			{
				int classbit_of_format = (scb.color_formats[i] >> 2) - low_class;
				encoded_type |= classbit_of_format << bitpos;
				bitpos++;
			}

			for (int i = 0; i < partition_count; i++)
			{
				int lowbits_of_format = scb.color_formats[i] & 3;
				encoded_type |= lowbits_of_format << bitpos;
				bitpos += 2;
			}

			int encoded_type_lowpart = encoded_type & 0x3F;
			int encoded_type_highpart = encoded_type >> 6;
			int encoded_type_highpart_size = (3 * partition_count) - 4;
			int encoded_type_highpart_pos = 128 - bits_for_weights - encoded_type_highpart_size;
			write_bits(encoded_type_lowpart, 6, 13 + PARTITION_BITS, pcb.data);
			write_bits(encoded_type_highpart, encoded_type_highpart_size, encoded_type_highpart_pos, pcb.data);
			below_weights_pos -= encoded_type_highpart_size;
		}
	}
	else
	{
		write_bits(scb.color_formats[0], 4, 13, pcb.data);
	}

	// in dual-plane mode, encode the color component of the second plane of weights
	if (is_dual_plane)
	{
		write_bits(scb.plane2_color_component, 2, below_weights_pos - 2, pcb.data);
	}

	// finally, encode the color bits
	// first, get hold of all the color components to encode
	uint8_t values_to_encode[32];
	int valuecount_to_encode = 0;
	for (int i = 0; i < scb.partition_count; i++)
	{
		int vals = 2 * (scb.color_formats[i] >> 2) + 2;
		for (int j = 0; j < vals; j++)
		{
			values_to_encode[j + valuecount_to_encode] = scb.color_values[i][j];
		}
		valuecount_to_encode += vals;
	}

	// then, encode an ISE based on them.
	encode_ise(scb.color_quantization_level, valuecount_to_encode, values_to_encode, pcb.data, (scb.partition_count == 1 ? 17 : 19 + PARTITION_BITS));
}

void physical_to_symbolic(
	const block_size_descriptor& bsd,
	const physical_compressed_block& pcb,
	symbolic_compressed_block& scb
) {
	uint8_t bswapped[16];

	scb.error_block = 0;

	// get hold of the decimation tables.
	const decimation_table *const *ixtab2 = bsd.decimation_tables;

	// extract header fields
	int block_mode = read_bits(11, 0, pcb.data);
	if ((block_mode & 0x1FF) == 0x1FC)
	{
		// void-extent block!

		// check what format the data has
		if (block_mode & 0x200)
		{
			scb.block_mode = -1;	// floating-point
		}
		else
		{
			scb.block_mode = -2;	// unorm16.
		}

		scb.partition_count = 0;
		for (int i = 0; i < 4; i++)
		{
			scb.constant_color[i] = pcb.data[2 * i + 8] | (pcb.data[2 * i + 9] << 8);
		}

		// additionally, check that the void-extent
		if (bsd.zdim == 1)
		{
			// 2D void-extent
			int rsvbits = read_bits(2, 10, pcb.data);
			if (rsvbits != 3)
			{
				scb.error_block = 1;
			}

			int vx_low_s = read_bits(8, 12, pcb.data) | (read_bits(5, 12 + 8, pcb.data) << 8);
			int vx_high_s = read_bits(8, 25, pcb.data) | (read_bits(5, 25 + 8, pcb.data) << 8);
			int vx_low_t = read_bits(8, 38, pcb.data) | (read_bits(5, 38 + 8, pcb.data) << 8);
			int vx_high_t = read_bits(8, 51, pcb.data) | (read_bits(5, 51 + 8, pcb.data) << 8);

			int all_ones = vx_low_s == 0x1FFF && vx_high_s == 0x1FFF && vx_low_t == 0x1FFF && vx_high_t == 0x1FFF;

			if ((vx_low_s >= vx_high_s || vx_low_t >= vx_high_t) && !all_ones)
			{
				scb.error_block = 1;
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
				scb.error_block = 1;
			}
		}

		return;
	}

	const int packed_index = bsd.block_mode_to_packed[block_mode];
	if (packed_index < 0)
	{
		scb.error_block = 1;
		return;
	}
	assert(packed_index >= 0 && packed_index < bsd.block_mode_packed_count);
	const struct block_mode& bm = bsd.block_modes_packed[packed_index];

	int weight_count = ixtab2[bm.decimation_mode]->num_weights;
	int weight_quantization_method = bm.quantization_mode;
	int is_dual_plane = bm.is_dual_plane;

	int real_weight_count = is_dual_plane ? 2 * weight_count : weight_count;

	int partition_count = read_bits(2, 11, pcb.data) + 1;

	scb.block_mode = block_mode;
	scb.partition_count = partition_count;

	for (int i = 0; i < 16; i++)
	{
		bswapped[i] = bitrev8(pcb.data[15 - i]);
	}

	int bits_for_weights = compute_ise_bitcount(real_weight_count,
												(quantization_method) weight_quantization_method);

	int below_weights_pos = 128 - bits_for_weights;

	if (is_dual_plane)
	{
		uint8_t indices[64];
		decode_ise(weight_quantization_method, real_weight_count, bswapped, indices, 0);
		for (int i = 0; i < weight_count; i++)
		{
			scb.plane1_weights[i] = indices[2 * i];
			scb.plane2_weights[i] = indices[2 * i + 1];
		}
	}
	else
	{
		decode_ise(weight_quantization_method, weight_count, bswapped, scb.plane1_weights, 0);
	}

	if (is_dual_plane && partition_count == 4)
	{
		scb.error_block = 1;
	}

	scb.color_formats_matched = 0;

	// then, determine the format of each endpoint pair
	int color_formats[4];
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
		int encoded_type = read_bits(6, 13 + PARTITION_BITS, pcb.data) | (read_bits(encoded_type_highpart_size, below_weights_pos, pcb.data) << 6);
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
		scb.partition_index = read_bits(6, 13, pcb.data) | (read_bits(PARTITION_BITS - 6, 19, pcb.data) << 6);
	}

	for (int i = 0; i < partition_count; i++)
	{
		scb.color_formats[i] = color_formats[i];
	}

	// then, determine the number of integers we need to unpack for the endpoint pairs
	int color_integer_count = 0;
	for (int i = 0; i < partition_count; i++)
	{
		int endpoint_class = color_formats[i] >> 2;
		color_integer_count += (endpoint_class + 1) * 2;
	}

	if (color_integer_count > 18)
	{
		scb.error_block = 1;
	}

	// then, determine the color endpoint format to use for these integers
	static const int color_bits_arr[5] = { -1, 115 - 4, 113 - 4 - PARTITION_BITS, 113 - 4 - PARTITION_BITS, 113 - 4 - PARTITION_BITS };
	int color_bits = color_bits_arr[partition_count] - bits_for_weights - encoded_type_highpart_size;
	if (is_dual_plane)
	{
		color_bits -= 2;
	}

	if (color_bits < 0)
	{
		color_bits = 0;
	}

	int color_quantization_level = quantization_mode_table[color_integer_count >> 1][color_bits];
	scb.color_quantization_level = color_quantization_level;
	if (color_quantization_level < 4)
	{
		scb.error_block = 1;
	}

	// then unpack the integer-bits
	uint8_t values_to_decode[32];
	decode_ise(color_quantization_level, color_integer_count, pcb.data, values_to_decode, (partition_count == 1 ? 17 : 19 + PARTITION_BITS));

	// and distribute them over the endpoint types
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

	// get hold of color component for second-plane in the case of dual plane of weights.
	if (is_dual_plane)
	{
		scb.plane2_color_component = read_bits(2, below_weights_pos - 2, pcb.data);
	}
}
