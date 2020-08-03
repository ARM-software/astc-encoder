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
 * @brief Functions to generate block size descriptor and decimation tables.
 */

#include "astcenc_internal.h"

// return 0 on invalid mode, 1 on valid mode.
static int decode_block_mode_2d(
	int blockmode,
	int* Nval,
	int* Mval,
	int* dual_weight_plane,
	int* quant_mode
) {
	int base_quant_mode = (blockmode >> 4) & 1;
	int H = (blockmode >> 9) & 1;
	int D = (blockmode >> 10) & 1;

	int A = (blockmode >> 5) & 0x3;

	int N = 0, M = 0;

	if ((blockmode & 3) != 0)
	{
		base_quant_mode |= (blockmode & 3) << 1;
		int B = (blockmode >> 7) & 3;
		switch ((blockmode >> 2) & 3)
		{
		case 0:
			N = B + 4;
			M = A + 2;
			break;
		case 1:
			N = B + 8;
			M = A + 2;
			break;
		case 2:
			N = A + 2;
			M = B + 8;
			break;
		case 3:
			B &= 1;
			if (blockmode & 0x100)
			{
				N = B + 2;
				M = A + 2;
			}
			else
			{
				N = A + 2;
				M = B + 6;
			}
			break;
		}
	}
	else
	{
		base_quant_mode |= ((blockmode >> 2) & 3) << 1;
		if (((blockmode >> 2) & 3) == 0)
		{
			return 0;
		}

		int B = (blockmode >> 9) & 3;
		switch ((blockmode >> 7) & 3)
		{
		case 0:
			N = 12;
			M = A + 2;
			break;
		case 1:
			N = A + 2;
			M = 12;
			break;
		case 2:
			N = A + 6;
			M = B + 6;
			D = 0;
			H = 0;
			break;
		case 3:
			switch ((blockmode >> 5) & 3)
			{
			case 0:
				N = 6;
				M = 10;
				break;
			case 1:
				N = 10;
				M = 6;
				break;
			case 2:
			case 3:
				return 0;
			}
			break;
		}
	}

	int weight_count = N * M * (D + 1);
	int qmode = (base_quant_mode - 2) + 6 * H;

	int weightbits = compute_ise_bitcount(weight_count, (quantization_method) qmode);
	if (weight_count > MAX_WEIGHTS_PER_BLOCK || weightbits < MIN_WEIGHT_BITS_PER_BLOCK || weightbits > MAX_WEIGHT_BITS_PER_BLOCK)
	{
		return 0;
	}

	*Nval = N;
	*Mval = M;
	*dual_weight_plane = D;
	*quant_mode = qmode;
	return 1;
}

static int decode_block_mode_3d(
	int blockmode,
	int* Nval,
	int* Mval,
	int* Qval,
	int* dual_weight_plane,
	int* quant_mode
) {
	int base_quant_mode = (blockmode >> 4) & 1;
	int H = (blockmode >> 9) & 1;
	int D = (blockmode >> 10) & 1;

	int A = (blockmode >> 5) & 0x3;

	int N = 0, M = 0, Q = 0;

	if ((blockmode & 3) != 0)
	{
		base_quant_mode |= (blockmode & 3) << 1;
		int B = (blockmode >> 7) & 3;
		int C = (blockmode >> 2) & 0x3;
		N = A + 2;
		M = B + 2;
		Q = C + 2;
	}
	else
	{
		base_quant_mode |= ((blockmode >> 2) & 3) << 1;
		if (((blockmode >> 2) & 3) == 0)
		{
			return 0;
		}

		int B = (blockmode >> 9) & 3;
		if (((blockmode >> 7) & 3) != 3)
		{
			D = 0;
			H = 0;
		}
		switch ((blockmode >> 7) & 3)
		{
		case 0:
			N = 6;
			M = B + 2;
			Q = A + 2;
			break;
		case 1:
			N = A + 2;
			M = 6;
			Q = B + 2;
			break;
		case 2:
			N = A + 2;
			M = B + 2;
			Q = 6;
			break;
		case 3:
			N = 2;
			M = 2;
			Q = 2;
			switch ((blockmode >> 5) & 3)
			{
			case 0:
				N = 6;
				break;
			case 1:
				M = 6;
				break;
			case 2:
				Q = 6;
				break;
			case 3:
				return 0;
			}
			break;
		}
	}

	int weight_count = N * M * Q * (D + 1);
	int qmode = (base_quant_mode - 2) + 6 * H;

	int weightbits = compute_ise_bitcount(weight_count, (quantization_method) qmode);
	if (weight_count > MAX_WEIGHTS_PER_BLOCK ||
	    weightbits < MIN_WEIGHT_BITS_PER_BLOCK ||
	    weightbits > MAX_WEIGHT_BITS_PER_BLOCK)
	{
		return 0;
	}

	*Nval = N;
	*Mval = M;
	*Qval = Q;
	*dual_weight_plane = D;
	*quant_mode = qmode;
	return 1;
}

static void initialize_decimation_table_2d(
	int xdim,
	int ydim,
	int x_weights,
	int y_weights,
	decimation_table* dt
) {
	int texels_per_block = xdim * ydim;
	int weights_per_block = x_weights * y_weights;

	int weightcount_of_texel[MAX_TEXELS_PER_BLOCK];
	int grid_weights_of_texel[MAX_TEXELS_PER_BLOCK][4];
	int weights_of_texel[MAX_TEXELS_PER_BLOCK][4];

	int texelcount_of_weight[MAX_WEIGHTS_PER_BLOCK];
	int texels_of_weight[MAX_WEIGHTS_PER_BLOCK][MAX_TEXELS_PER_BLOCK];
	int texelweights_of_weight[MAX_WEIGHTS_PER_BLOCK][MAX_TEXELS_PER_BLOCK];

	for (int i = 0; i < weights_per_block; i++)
	{
		texelcount_of_weight[i] = 0;
	}

	for (int i = 0; i < texels_per_block; i++)
	{
		weightcount_of_texel[i] = 0;
	}

	for (int y = 0; y < ydim; y++)
	{
		for (int x = 0; x < xdim; x++)
		{
			int texel = y * xdim + x;

			int x_weight = (((1024 + xdim / 2) / (xdim - 1)) * x * (x_weights - 1) + 32) >> 6;
			int y_weight = (((1024 + ydim / 2) / (ydim - 1)) * y * (y_weights - 1) + 32) >> 6;

			int x_weight_frac = x_weight & 0xF;
			int y_weight_frac = y_weight & 0xF;
			int x_weight_int = x_weight >> 4;
			int y_weight_int = y_weight >> 4;
			int qweight[4];
			int weight[4];
			qweight[0] = x_weight_int + y_weight_int * x_weights;
			qweight[1] = qweight[0] + 1;
			qweight[2] = qweight[0] + x_weights;
			qweight[3] = qweight[2] + 1;

			// truncated-precision bilinear interpolation.
			int prod = x_weight_frac * y_weight_frac;

			weight[3] = (prod + 8) >> 4;
			weight[1] = x_weight_frac - weight[3];
			weight[2] = y_weight_frac - weight[3];
			weight[0] = 16 - x_weight_frac - y_weight_frac + weight[3];

			for (int i = 0; i < 4; i++)
			{
				if (weight[i] != 0)
				{
					grid_weights_of_texel[texel][weightcount_of_texel[texel]] = qweight[i];
					weights_of_texel[texel][weightcount_of_texel[texel]] = weight[i];
					weightcount_of_texel[texel]++;
					texels_of_weight[qweight[i]][texelcount_of_weight[qweight[i]]] = texel;
					texelweights_of_weight[qweight[i]][texelcount_of_weight[qweight[i]]] = weight[i];
					texelcount_of_weight[qweight[i]]++;
				}
			}
		}
	}

	for (int i = 0; i < texels_per_block; i++)
	{
		dt->texel_num_weights[i] = weightcount_of_texel[i];

		// ensure that all 4 entries are actually initialized.
		// This allows a branch-free implementation of compute_value_of_texel_flt()
		for (int j = 0; j < 4; j++)
		{
			dt->texel_weights_int[i][j] = 0;
			dt->texel_weights_float[i][j] = 0.0f;
			dt->texel_weights[i][j] = 0;
		}

		for (int j = 0; j < weightcount_of_texel[i]; j++)
		{
			dt->texel_weights_int[i][j] = (uint8_t)weights_of_texel[i][j];
			dt->texel_weights_float[i][j] = ((float)weights_of_texel[i][j]) * (1.0f / TEXEL_WEIGHT_SUM);
			dt->texel_weights[i][j] = (uint8_t)grid_weights_of_texel[i][j];
		}
	}

	for (int i = 0; i < weights_per_block; i++)
	{
		dt->weight_num_texels[i] = texelcount_of_weight[i];

		for (int j = 0; j < texelcount_of_weight[i]; j++)
		{
			int texel = texels_of_weight[i][j];
			dt->weight_texel[i][j] = (uint8_t)texel;
			dt->weights_int[i][j] = (uint8_t)texelweights_of_weight[i][j];
			dt->weights_flt[i][j] = (float)texelweights_of_weight[i][j];

			// perform a layer of array unrolling. An aspect of this unrolling is that
			// one of the texel-weight indexes is an identity-mapped index; we will use this
			// fact to reorder the indexes so that the first one is the identity index.
			int swap_idx = -1;
			for (int k = 0; k < 4; k++)
			{
				int dttw = dt->texel_weights[texel][k];
				float dttwf = dt->texel_weights_float[texel][k];
				if (dttw == i && dttwf != 0.0f)
				{
					swap_idx = k;
				}
				dt->texel_weights_texel[i][j][k] = (uint8_t)dttw;
				dt->texel_weights_float_texel[i][j][k] = dttwf;
			}

			if (swap_idx != 0)
			{
				int vi = dt->texel_weights_texel[i][j][0];
				float vf = dt->texel_weights_float_texel[i][j][0];
				dt->texel_weights_texel[i][j][0] = dt->texel_weights_texel[i][j][swap_idx];
				dt->texel_weights_float_texel[i][j][0] = dt->texel_weights_float_texel[i][j][swap_idx];
				dt->texel_weights_texel[i][j][swap_idx] = (uint8_t)vi;
				dt->texel_weights_float_texel[i][j][swap_idx] = vf;
			}
		}
	}

	dt->num_texels = texels_per_block;
	dt->num_weights = weights_per_block;
}

static void initialize_decimation_table_3d(
	int xdim,
	int ydim,
	int zdim,
	int x_weights,
	int y_weights,
	int z_weights,
	decimation_table* dt
) {
	int texels_per_block = xdim * ydim * zdim;
	int weights_per_block = x_weights * y_weights * z_weights;

	int weightcount_of_texel[MAX_TEXELS_PER_BLOCK];
	int grid_weights_of_texel[MAX_TEXELS_PER_BLOCK][4];
	int weights_of_texel[MAX_TEXELS_PER_BLOCK][4];

	int texelcount_of_weight[MAX_WEIGHTS_PER_BLOCK];
	int texels_of_weight[MAX_WEIGHTS_PER_BLOCK][MAX_TEXELS_PER_BLOCK];
	int texelweights_of_weight[MAX_WEIGHTS_PER_BLOCK][MAX_TEXELS_PER_BLOCK];

	for (int i = 0; i < weights_per_block; i++)
	{
		texelcount_of_weight[i] = 0;
	}

	for (int i = 0; i < texels_per_block; i++)
	{
		weightcount_of_texel[i] = 0;
	}

	for (int z = 0; z < zdim; z++)
	{
		for (int y = 0; y < ydim; y++)
		{
			for (int x = 0; x < xdim; x++)
			{
				int texel = (z * ydim + y) * xdim + x;

				int x_weight = (((1024 + xdim / 2) / (xdim - 1)) * x * (x_weights - 1) + 32) >> 6;
				int y_weight = (((1024 + ydim / 2) / (ydim - 1)) * y * (y_weights - 1) + 32) >> 6;
				int z_weight = (((1024 + zdim / 2) / (zdim - 1)) * z * (z_weights - 1) + 32) >> 6;

				int x_weight_frac = x_weight & 0xF;
				int y_weight_frac = y_weight & 0xF;
				int z_weight_frac = z_weight & 0xF;
				int x_weight_int = x_weight >> 4;
				int y_weight_int = y_weight >> 4;
				int z_weight_int = z_weight >> 4;
				int qweight[4];
				int weight[4];
				qweight[0] = (z_weight_int * y_weights + y_weight_int) * x_weights + x_weight_int;
				qweight[3] = ((z_weight_int + 1) * y_weights + (y_weight_int + 1)) * x_weights + (x_weight_int + 1);

				// simplex interpolation
				int fs = x_weight_frac;
				int ft = y_weight_frac;
				int fp = z_weight_frac;

				int cas = ((fs > ft) << 2) + ((ft > fp) << 1) + ((fs > fp));
				int N = x_weights;
				int NM = x_weights * y_weights;

				int s1, s2, w0, w1, w2, w3;
				switch (cas)
				{
				case 7:
					s1 = 1;
					s2 = N;
					w0 = 16 - fs;
					w1 = fs - ft;
					w2 = ft - fp;
					w3 = fp;
					break;
				case 3:
					s1 = N;
					s2 = 1;
					w0 = 16 - ft;
					w1 = ft - fs;
					w2 = fs - fp;
					w3 = fp;
					break;
				case 5:
					s1 = 1;
					s2 = NM;
					w0 = 16 - fs;
					w1 = fs - fp;
					w2 = fp - ft;
					w3 = ft;
					break;
				case 4:
					s1 = NM;
					s2 = 1;
					w0 = 16 - fp;
					w1 = fp - fs;
					w2 = fs - ft;
					w3 = ft;
					break;
				case 2:
					s1 = N;
					s2 = NM;
					w0 = 16 - ft;
					w1 = ft - fp;
					w2 = fp - fs;
					w3 = fs;
					break;
				case 0:
					s1 = NM;
					s2 = N;
					w0 = 16 - fp;
					w1 = fp - ft;
					w2 = ft - fs;
					w3 = fs;
					break;
				default:
					s1 = NM;
					s2 = N;
					w0 = 16 - fp;
					w1 = fp - ft;
					w2 = ft - fs;
					w3 = fs;
					break;
				}

				qweight[1] = qweight[0] + s1;
				qweight[2] = qweight[1] + s2;
				weight[0] = w0;
				weight[1] = w1;
				weight[2] = w2;
				weight[3] = w3;

				for (int i = 0; i < 4; i++)
				{
					if (weight[i] != 0)
					{
						grid_weights_of_texel[texel][weightcount_of_texel[texel]] = qweight[i];
						weights_of_texel[texel][weightcount_of_texel[texel]] = weight[i];
						weightcount_of_texel[texel]++;
						texels_of_weight[qweight[i]][texelcount_of_weight[qweight[i]]] = texel;
						texelweights_of_weight[qweight[i]][texelcount_of_weight[qweight[i]]] = weight[i];
						texelcount_of_weight[qweight[i]]++;
					}
				}
			}
		}
	}

	for (int i = 0; i < texels_per_block; i++)
	{
		dt->texel_num_weights[i] = weightcount_of_texel[i];

		// ensure that all 4 entries are actually initialized.
		// This allows a branch-free implementation of compute_value_of_texel_flt()
		for (int j = 0; j < 4; j++)
		{
			dt->texel_weights_int[i][j] = 0;
			dt->texel_weights_float[i][j] = 0.0f;
			dt->texel_weights[i][j] = 0;
		}

		for (int j = 0; j < weightcount_of_texel[i]; j++)
		{
			dt->texel_weights_int[i][j] = (uint8_t)weights_of_texel[i][j];
			dt->texel_weights_float[i][j] = ((float)weights_of_texel[i][j]) * (1.0f / TEXEL_WEIGHT_SUM);
			dt->texel_weights[i][j] = (uint8_t)grid_weights_of_texel[i][j];
		}
	}

	for (int i = 0; i < weights_per_block; i++)
	{
		dt->weight_num_texels[i] = texelcount_of_weight[i];
		for (int j = 0; j < texelcount_of_weight[i]; j++)
		{
			int texel = texels_of_weight[i][j];
			dt->weight_texel[i][j] = (uint8_t)texel;
			dt->weights_int[i][j] = (uint8_t)texelweights_of_weight[i][j];
			dt->weights_flt[i][j] = (float)texelweights_of_weight[i][j];

			// perform a layer of array unrolling. An aspect of this unrolling is that
			// one of the texel-weight indexes is an identity-mapped index; we will use this
			// fact to reorder the indexes so that the first one is the identity index.
			int swap_idx = -1;
			for (int k = 0; k < 4; k++)
			{
				int dttw = dt->texel_weights[texel][k];
				float dttwf = dt->texel_weights_float[texel][k];
				if (dttw == i && dttwf != 0.0f)
				{
					swap_idx = k;
				}
				dt->texel_weights_texel[i][j][k] = (uint8_t)dttw;
				dt->texel_weights_float_texel[i][j][k] = dttwf;
			}

			if (swap_idx != 0)
			{
				int vi = dt->texel_weights_texel[i][j][0];
				float vf = dt->texel_weights_float_texel[i][j][0];
				dt->texel_weights_texel[i][j][0] = dt->texel_weights_texel[i][j][swap_idx];
				dt->texel_weights_float_texel[i][j][0] = dt->texel_weights_float_texel[i][j][swap_idx];
				dt->texel_weights_texel[i][j][swap_idx] = (uint8_t)vi;
				dt->texel_weights_float_texel[i][j][swap_idx] = vf;
			}
		}
	}

	dt->num_texels = texels_per_block;
	dt->num_weights = weights_per_block;
}

static void construct_block_size_descriptor_2d(
	int xdim,
	int ydim,
	block_size_descriptor* bsd
) {
	int decimation_mode_index[256];	// for each of the 256 entries in the decim_table_array, its index
	int decimation_mode_count = 0;

	bsd->xdim = xdim;
	bsd->ydim = ydim;
	bsd->zdim = 1;
	bsd->texel_count = xdim * ydim;

	for (int i = 0; i < 256; i++)
	{
		decimation_mode_index[i] = -1;
	}

	// gather all the infill-modes that can be used with the current block size
	for (int x_weights = 2; x_weights <= 12; x_weights++)
	{
		for (int y_weights = 2; y_weights <= 12; y_weights++)
		{
			if (x_weights * y_weights > MAX_WEIGHTS_PER_BLOCK)
			{
				continue;
			}

			decimation_table *dt = new decimation_table;
			decimation_mode_index[y_weights * 16 + x_weights] = decimation_mode_count;
			initialize_decimation_table_2d(xdim, ydim, x_weights, y_weights, dt);

			int weight_count = x_weights * y_weights;

			int maxprec_1plane = -1;
			int maxprec_2planes = -1;
			for (int i = 0; i < 12; i++)
			{
				int bits_1plane = compute_ise_bitcount(weight_count, (quantization_method) i);
				int bits_2planes = compute_ise_bitcount(2 * weight_count, (quantization_method) i);

				if (bits_1plane >= MIN_WEIGHT_BITS_PER_BLOCK && bits_1plane <= MAX_WEIGHT_BITS_PER_BLOCK)
				{
					maxprec_1plane = i;
				}

				if (bits_2planes >= MIN_WEIGHT_BITS_PER_BLOCK && bits_2planes <= MAX_WEIGHT_BITS_PER_BLOCK)
				{
					maxprec_2planes = i;
				}
			}

			if (2 * x_weights * y_weights > MAX_WEIGHTS_PER_BLOCK)
			{
				maxprec_2planes = -1;
			}

			bsd->permit_encode[decimation_mode_count] = (x_weights <= xdim && y_weights <= ydim);

			bsd->decimation_mode_samples[decimation_mode_count] = weight_count;
			bsd->decimation_mode_maxprec_1plane[decimation_mode_count] = maxprec_1plane;
			bsd->decimation_mode_maxprec_2planes[decimation_mode_count] = maxprec_2planes;
			bsd->decimation_tables[decimation_mode_count] = dt;

			decimation_mode_count++;
		}
	}

	for (int i = 0; i < MAX_DECIMATION_MODES; i++)
	{
		bsd->decimation_mode_percentile[i] = 1.0f;
	}

	for (int i = decimation_mode_count; i < MAX_DECIMATION_MODES; i++)
	{
		bsd->permit_encode[i] = 0;
		bsd->decimation_mode_samples[i] = 0;
		bsd->decimation_mode_maxprec_1plane[i] = -1;
		bsd->decimation_mode_maxprec_2planes[i] = -1;
	}

	bsd->decimation_mode_count = decimation_mode_count;

#if !defined(ASTCENC_DECOMPRESS_ONLY)
	const float *percentiles = get_2d_percentile_table(xdim, ydim);
#endif

	// then construct the list of block formats
	for (int i = 0; i < 2048; i++)
	{
		int x_weights, y_weights;
		int is_dual_plane;
		int quantization_mode;
		int fail = 0;
		int permit_encode = 1;

		if (decode_block_mode_2d(i, &x_weights, &y_weights, &is_dual_plane, &quantization_mode))
		{
			if (x_weights > xdim || y_weights > ydim)
			{
				permit_encode = 0;
			}
		}
		else
		{
			fail = 1;
			permit_encode = 0;
		}

		if (fail)
		{
			bsd->block_modes[i].decimation_mode = -1;
			bsd->block_modes[i].quantization_mode = -1;
			bsd->block_modes[i].is_dual_plane = -1;
			bsd->block_modes[i].permit_encode = 0;
			bsd->block_modes[i].permit_decode = 0;
			bsd->block_modes[i].percentile = 1.0f;
		}
		else
		{
			int decimation_mode = decimation_mode_index[y_weights * 16 + x_weights];
			bsd->block_modes[i].decimation_mode = decimation_mode;
			bsd->block_modes[i].quantization_mode = quantization_mode;
			bsd->block_modes[i].is_dual_plane = is_dual_plane;
			bsd->block_modes[i].permit_encode = permit_encode;
			bsd->block_modes[i].permit_decode = permit_encode;	// disallow decode of grid size larger than block size.

#if !defined(ASTCENC_DECOMPRESS_ONLY)
			bsd->block_modes[i].percentile = percentiles[i];
			if (bsd->decimation_mode_percentile[decimation_mode] > percentiles[i])
			{
				bsd->decimation_mode_percentile[decimation_mode] = percentiles[i];
			}
#else
			bsd->block_modes[i].percentile = 0.0f;
#endif
		}
	}

#if !defined(ASTCENC_DECOMPRESS_ONLY)
	delete[] percentiles;
#endif

	if (xdim * ydim <= 64)
	{
		bsd->texelcount_for_bitmap_partitioning = xdim * ydim;
		for (int i = 0; i < xdim * ydim; i++)
		{
			bsd->texels_for_bitmap_partitioning[i] = i;
		}
	}
	else
	{
		uint64_t rng_state[2];
		astc::rand_init(rng_state);

		// pick 64 random texels for use with bitmap partitioning.
		int arr[MAX_TEXELS_PER_BLOCK];
		for (int i = 0; i < xdim * ydim; i++)
		{
			arr[i] = 0;
		}

		int arr_elements_set = 0;
		while (arr_elements_set < 64)
		{
			unsigned int idx = (unsigned int)astc::rand(rng_state);
			idx %= xdim * ydim;
			if (arr[idx] == 0)
			{
				arr_elements_set++;
				arr[idx] = 1;
			}
		}

		int texel_weights_written = 0;
		int idx = 0;
		while (texel_weights_written < 64)
		{
			if (arr[idx])
			{
				bsd->texels_for_bitmap_partitioning[texel_weights_written++] = idx;
			}
			idx++;
		}

		bsd->texelcount_for_bitmap_partitioning = 64;
	}
}

static void construct_block_size_descriptor_3d(
	int xdim,
	int ydim,
	int zdim,
	block_size_descriptor * bsd
) {
	int decimation_mode_index[512];	// for each of the 512 entries in the decim_table_array, its index
	int decimation_mode_count = 0;

	bsd->xdim = xdim;
	bsd->ydim = ydim;
	bsd->zdim = zdim;
	bsd->texel_count = xdim * ydim * zdim;

	for (int i = 0; i < 512; i++)
	{
		decimation_mode_index[i] = -1;
	}

	// gather all the infill-modes that can be used with the current block size
	for (int x_weights = 2; x_weights <= 6; x_weights++)
	{
		for (int y_weights = 2; y_weights <= 6; y_weights++)
		{
			for (int z_weights = 2; z_weights <= 6; z_weights++)
			{
				if ((x_weights * y_weights * z_weights) > MAX_WEIGHTS_PER_BLOCK)
				{
					continue;
				}

				decimation_table *dt = new decimation_table;
				decimation_mode_index[z_weights * 64 + y_weights * 8 + x_weights] = decimation_mode_count;
				initialize_decimation_table_3d(xdim, ydim, zdim, x_weights, y_weights, z_weights, dt);

				int weight_count = x_weights * y_weights * z_weights;

				int maxprec_1plane = -1;
				int maxprec_2planes = -1;
				for (int i = 0; i < 12; i++)
				{
					int bits_1plane = compute_ise_bitcount(weight_count, (quantization_method) i);
					int bits_2planes = compute_ise_bitcount(2 * weight_count, (quantization_method) i);

					if (bits_1plane >= MIN_WEIGHT_BITS_PER_BLOCK && bits_1plane <= MAX_WEIGHT_BITS_PER_BLOCK)
					{
						maxprec_1plane = i;
					}

					if (bits_2planes >= MIN_WEIGHT_BITS_PER_BLOCK && bits_2planes <= MAX_WEIGHT_BITS_PER_BLOCK)
					{
						maxprec_2planes = i;
					}
				}

				if ((2 * x_weights * y_weights * z_weights) > MAX_WEIGHTS_PER_BLOCK)
				{
					maxprec_2planes = -1;
				}

				bsd->permit_encode[decimation_mode_count] = (x_weights <= xdim && y_weights <= ydim && z_weights <= zdim);

				bsd->decimation_mode_samples[decimation_mode_count] = weight_count;
				bsd->decimation_mode_maxprec_1plane[decimation_mode_count] = maxprec_1plane;
				bsd->decimation_mode_maxprec_2planes[decimation_mode_count] = maxprec_2planes;
				bsd->decimation_tables[decimation_mode_count] = dt;

				decimation_mode_count++;
			}
		}
	}

	for (int i = 0; i < MAX_DECIMATION_MODES; i++)
	{
		bsd->decimation_mode_percentile[i] = 1.0f;
	}

	for (int i = decimation_mode_count; i < MAX_DECIMATION_MODES; i++)
	{
		bsd->permit_encode[i] = 0;
		bsd->decimation_mode_samples[i] = 0;
		bsd->decimation_mode_maxprec_1plane[i] = -1;
		bsd->decimation_mode_maxprec_2planes[i] = -1;
	}

	bsd->decimation_mode_count = decimation_mode_count;

	// then construct the list of block formats
	for (int i = 0; i < 2048; i++)
	{
		int x_weights, y_weights, z_weights;
		int is_dual_plane;
		int quantization_mode;
		int fail = 0;
		int permit_encode = 1;

		if (decode_block_mode_3d(i, &x_weights, &y_weights, &z_weights, &is_dual_plane, &quantization_mode))
		{
			if (x_weights > xdim || y_weights > ydim || z_weights > zdim)
			{
				permit_encode = 0;
			}
		}
		else
		{
			fail = 1;
			permit_encode = 0;
		}
		if (fail)
		{
			bsd->block_modes[i].decimation_mode = -1;
			bsd->block_modes[i].quantization_mode = -1;
			bsd->block_modes[i].is_dual_plane = -1;
			bsd->block_modes[i].permit_encode = 0;
			bsd->block_modes[i].permit_decode = 0;
			bsd->block_modes[i].percentile = 1.0f;
		}
		else
		{
			int decimation_mode = decimation_mode_index[z_weights * 64 + y_weights * 8 + x_weights];
			bsd->block_modes[i].decimation_mode = decimation_mode;
			bsd->block_modes[i].quantization_mode = quantization_mode;
			bsd->block_modes[i].is_dual_plane = is_dual_plane;
			bsd->block_modes[i].permit_encode = permit_encode;
			bsd->block_modes[i].permit_decode = permit_encode;

			bsd->block_modes[i].percentile = 0.0f; // No percentile table
			if (bsd->decimation_mode_percentile[decimation_mode] > 0.0f)
			{
				bsd->decimation_mode_percentile[decimation_mode] = 0.0f;
			}
		}
	}

	if (xdim * ydim * zdim <= 64)
	{
		bsd->texelcount_for_bitmap_partitioning = xdim * ydim * zdim;
		for (int i = 0; i < xdim * ydim * zdim; i++)
		{
			bsd->texels_for_bitmap_partitioning[i] = i;
		}
	}
	else
	{
		uint64_t rng_state[2];
		astc::rand_init(rng_state);

		// pick 64 random texels for use with bitmap partitioning.
		int arr[MAX_TEXELS_PER_BLOCK];
		for (int i = 0; i < xdim * ydim * zdim; i++)
		{
			arr[i] = 0;
		}

		int arr_elements_set = 0;
		while (arr_elements_set < 64)
		{
			unsigned int idx = (unsigned int)astc::rand(rng_state);
			idx %= xdim * ydim * zdim;
			if (arr[idx] == 0)
			{
				arr_elements_set++;
				arr[idx] = 1;
			}
		}

		int texel_weights_written = 0;
		int idx = 0;
		while (texel_weights_written < 64)
		{
			if (arr[idx])
			{
				bsd->texels_for_bitmap_partitioning[texel_weights_written++] = idx;
			}
			idx++;
		}
		bsd->texelcount_for_bitmap_partitioning = 64;
	}
}

/* Public function, see header file for detailed documentation */
void init_block_size_descriptor(
	int xdim,
	int ydim,
	int zdim,
	block_size_descriptor* bsd
) {
	if (zdim > 1)
	{
		construct_block_size_descriptor_3d(xdim, ydim, zdim, bsd);
	}
	else
	{
		construct_block_size_descriptor_2d(xdim, ydim, bsd);
	}

	init_partition_tables(bsd);
}

void term_block_size_descriptor(
	block_size_descriptor* bsd)
{
	for (int i = 0; i < bsd->decimation_mode_count; i++)
	{
		delete bsd->decimation_tables[i];
	}
}
