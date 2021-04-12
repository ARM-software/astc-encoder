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
 * @brief Functions for creating in-memory ASTC image structures.
 */

#include <cassert>
#include <cstring>

#include "astcenc_internal.h"

void imageblock_initialize_deriv(
	const imageblock* blk,
	int pixelcount,
	vfloat4* dptr
) {
	// TODO: For LDR on the current codec we can skip this if no LNS and just
	// early-out as we use the same LNS settings everywhere ...
	for (int i = 0; i < pixelcount; i++)
	{
		vfloat4 derv_unorm(65535.0f);
		vfloat4 derv_lns = vfloat4::zero();

		// TODO: Pack these into bits and avoid the disjoint fetch
		int rgb_lns = blk->rgb_lns[i];
		int a_lns = blk->alpha_lns[i];

		// Compute derivatives if we have any use of LNS
		if (rgb_lns || a_lns)
		{
			vfloat4 data = blk->texel(i);
			vint4 datai = lns_to_sf16(float_to_int(data));

			vfloat4 dataf = float16_to_float(datai);
			dataf = max(dataf, 6e-5f);

			vfloat4 data_lns1 = dataf * 1.05f;
			data_lns1 = float_to_lns(data_lns1);

			vfloat4 data_lns2 = dataf;
			data_lns2 = float_to_lns(data_lns2);

			vfloat4 divisor_lns = dataf * 0.05f;

			// Clamp derivatives between 1/32 and 2^25
			float lo = 1.0f / 32.0f;
			float hi = 33554432.0f;
			derv_lns = clamp(lo, hi, (data_lns1 - data_lns2) / divisor_lns);
		}

		vint4 use_lns(rgb_lns, rgb_lns, rgb_lns, a_lns);
		vmask4 lns_mask = use_lns != vint4::zero();
		*dptr = select(derv_unorm, derv_lns, lns_mask);
		dptr++;
	}
}

// helper function to initialize the work-data from the orig-data
static void imageblock_initialize_work_from_orig(
	imageblock* blk,
	int pixelcount
) {
	blk->origin_texel = blk->texel(0);

	vfloat4 data_min(1e38f);
	vfloat4 data_max(-1e38f);
	bool grayscale = true;

	for (int i = 0; i < pixelcount; i++)
	{
		vfloat4 data = blk->texel(i);
		vfloat4 color_lns = vfloat4::zero();
		vfloat4 color_unorm = data * 65535.0f;

		int rgb_lns = blk->rgb_lns[i];
		int a_lns = blk->alpha_lns[i];

		if (rgb_lns || a_lns)
		{
			color_lns = float_to_lns(data);
		}

		vint4 use_lns(rgb_lns, rgb_lns, rgb_lns, a_lns);
		vmask4 lns_mask = use_lns != vint4::zero();
		data = select(color_unorm, color_lns, lns_mask);

		// Compute block metadata
		data_min = min(data_min, data);
		data_max = max(data_max, data);

		if (grayscale && (data.lane<0>() != data.lane<1>() || data.lane<0>() != data.lane<2>()))
		{
			grayscale = false;
		}

		// Store block data
		blk->data_r[i] = data.lane<0>();
		blk->data_g[i] = data.lane<1>();
		blk->data_b[i] = data.lane<2>();
		blk->data_a[i] = data.lane<3>();
	}

	// Store block metadata
	blk->data_min = data_min;
	blk->data_max = data_max;
	blk->grayscale = grayscale;
}

// fetch an imageblock from the input file.
void fetch_imageblock(
	astcenc_profile decode_mode,
	const astcenc_image& img,
	imageblock* blk,	// picture-block to initialize with image data
	const block_size_descriptor* bsd,
	// position in texture.
	int xpos,
	int ypos,
	int zpos,
	astcenc_swizzle swz
) {
	int xsize = img.dim_x;
	int ysize = img.dim_y;
	int zsize = img.dim_z;

	blk->xpos = xpos;
	blk->ypos = ypos;
	blk->zpos = zpos;

	// True if any non-identity swizzle
	bool needs_swz = (swz.r != ASTCENC_SWZ_R) || (swz.g != ASTCENC_SWZ_G) ||
	                 (swz.b != ASTCENC_SWZ_B) || (swz.a != ASTCENC_SWZ_A);

	int idx = 0;
	if (img.data_type == ASTCENC_TYPE_U8)
	{
		uint8_t data[6];
		data[ASTCENC_SWZ_0] = 0x00;
		data[ASTCENC_SWZ_1] = 0xFF;

		for (int z = 0; z < bsd->zdim; z++)
		{
			int zi = astc::min(zpos + z, zsize - 1);
			uint8_t* data8 = static_cast<uint8_t*>(img.data[zi]);

			for (int y = 0; y < bsd->ydim; y++)
			{
				int yi = astc::min(ypos + y, ysize - 1);

				for (int x = 0; x < bsd->xdim; x++)
				{
					int xi = astc::min(xpos + x, xsize - 1);

					int r = data8[(4 * xsize * yi) + (4 * xi    )];
					int g = data8[(4 * xsize * yi) + (4 * xi + 1)];
					int b = data8[(4 * xsize * yi) + (4 * xi + 2)];
					int a = data8[(4 * xsize * yi) + (4 * xi + 3)];

					if (needs_swz)
					{
						data[ASTCENC_SWZ_R] = r;
						data[ASTCENC_SWZ_G] = g;
						data[ASTCENC_SWZ_B] = b;
						data[ASTCENC_SWZ_A] = a;

						r = data[swz.r];
						g = data[swz.g];
						b = data[swz.b];
						a = data[swz.a];
					}

					blk->data_r[idx] = static_cast<float>(r) / 255.0f;
					blk->data_g[idx] = static_cast<float>(g) / 255.0f;
					blk->data_b[idx] = static_cast<float>(b) / 255.0f;
					blk->data_a[idx] = static_cast<float>(a) / 255.0f;
					idx++;
				}
			}
		}
	}
	else if (img.data_type == ASTCENC_TYPE_F16)
	{
		uint16_t data[6];
		data[ASTCENC_SWZ_0] = 0x0000;
		data[ASTCENC_SWZ_1] = 0x3C00;

		for (int z = 0; z < bsd->zdim; z++)
		{
			int zi = astc::min(zpos + z, zsize - 1);
			uint16_t* data16 = static_cast<uint16_t*>(img.data[zi]);

			for (int y = 0; y < bsd->ydim; y++)
			{
				int yi = astc::min(ypos + y, ysize - 1);

				for (int x = 0; x < bsd->xdim; x++)
				{
					int xi = astc::min(xpos + x, xsize - 1);

					int r = data16[(4 * xsize * yi) + (4 * xi    )];
					int g = data16[(4 * xsize * yi) + (4 * xi + 1)];
					int b = data16[(4 * xsize * yi) + (4 * xi + 2)];
					int a = data16[(4 * xsize * yi) + (4 * xi + 3)];

					if (needs_swz)
					{
						data[ASTCENC_SWZ_R] = r;
						data[ASTCENC_SWZ_G] = g;
						data[ASTCENC_SWZ_B] = b;
						data[ASTCENC_SWZ_A] = a;

						r = data[swz.r];
						g = data[swz.g];
						b = data[swz.b];
						a = data[swz.a];
					}

					vfloat4 dataf = max(float16_to_float(vint4(r, g, b, a)), 1e-8f);
					blk->data_r[idx] = dataf.lane<0>();
					blk->data_g[idx] = dataf.lane<1>();
					blk->data_b[idx] = dataf.lane<2>();
					blk->data_a[idx] = dataf.lane<3>();
					idx++;
				}
			}
		}
	}
	else // if (img.data_type == ASTCENC_TYPE_F32)
	{
		assert(img.data_type == ASTCENC_TYPE_F32);

		float data[6];
		data[ASTCENC_SWZ_0] = 0.0f;
		data[ASTCENC_SWZ_1] = 1.0f;

		for (int z = 0; z < bsd->zdim; z++)
		{
			int zi = astc::min(zpos + z, zsize - 1);
			float* data32 = static_cast<float*>(img.data[zi]);

			for (int y = 0; y < bsd->ydim; y++)
			{
				int yi = astc::min(ypos + y, ysize - 1);

				for (int x = 0; x < bsd->xdim; x++)
				{
					int xi = astc::min(xpos + x, xsize - 1);

					float r = data32[(4 * xsize * yi) + (4 * xi    )];
					float g = data32[(4 * xsize * yi) + (4 * xi + 1)];
					float b = data32[(4 * xsize * yi) + (4 * xi + 2)];
					float a = data32[(4 * xsize * yi) + (4 * xi + 3)];

					if (needs_swz)
					{
						data[ASTCENC_SWZ_R] = r;
						data[ASTCENC_SWZ_G] = g;
						data[ASTCENC_SWZ_B] = b;
						data[ASTCENC_SWZ_A] = a;

						r = data[swz.r];
						g = data[swz.g];
						b = data[swz.b];
						a = data[swz.a];
					}

					blk->data_r[idx] = astc::max(r, 1e-8f);
					blk->data_g[idx] = astc::max(g, 1e-8f);
					blk->data_b[idx] = astc::max(b, 1e-8f);
					blk->data_a[idx] = astc::max(a, 1e-8f);
					idx++;
				}
			}
		}
	}

	int rgb_lns = (decode_mode == ASTCENC_PRF_HDR) || (decode_mode == ASTCENC_PRF_HDR_RGB_LDR_A);
	int alpha_lns = decode_mode == ASTCENC_PRF_HDR;

	// impose the choice on every pixel when encoding.
	for (int i = 0; i < bsd->texel_count; i++)
	{
		blk->rgb_lns[i] = rgb_lns;
		blk->alpha_lns[i] = alpha_lns;
	}

	imageblock_initialize_work_from_orig(blk, bsd->texel_count);
}

void write_imageblock(
	astcenc_image& img,
	const imageblock* blk,	// picture-block to initialize with image data. We assume that orig_data is valid.
	const block_size_descriptor* bsd,
	// position to write the block to
	int xpos,
	int ypos,
	int zpos,
	astcenc_swizzle swz
) {
	int xsize = img.dim_x;
	int ysize = img.dim_y;
	int zsize = img.dim_z;

	int x_start = xpos;
	int x_end = std::min(xsize, xpos + bsd->xdim);
	int x_nudge = bsd->xdim - (x_end - x_start);

	int y_start = ypos;
	int y_end = std::min(ysize, ypos + bsd->ydim);
	int y_nudge = (bsd->ydim - (y_end - y_start)) * bsd->xdim;

	int z_start = zpos;
	int z_end = std::min(zsize, zpos + bsd->zdim);

	float data[7];
	data[ASTCENC_SWZ_0] = 0.0f;
	data[ASTCENC_SWZ_1] = 1.0f;

	// True if any non-identity swizzle
	bool needs_swz = (swz.r != ASTCENC_SWZ_R) || (swz.g != ASTCENC_SWZ_G) ||
	                 (swz.b != ASTCENC_SWZ_B) || (swz.a != ASTCENC_SWZ_A);

	// True if any swizzle uses Z reconstruct
	bool needs_z = (swz.r == ASTCENC_SWZ_Z) || (swz.g == ASTCENC_SWZ_Z) ||
	               (swz.b == ASTCENC_SWZ_Z) || (swz.a == ASTCENC_SWZ_Z);

	int idx = 0;
	if (img.data_type == ASTCENC_TYPE_U8)
	{
		for (int z = z_start; z < z_end; z++)
		{
			// Fetch the image plane
			uint8_t* data8 = static_cast<uint8_t*>(img.data[z]);

			for (int y = y_start; y < y_end; y++)
			{
				for (int x = x_start; x < x_end; x++)
				{
					vint4 colori = vint4::zero();

					if (blk->data_r[idx] == std::numeric_limits<float>::quiet_NaN())
					{
						// Can't display NaN - show magenta error color
						colori = vint4(0xFF, 0x00, 0xFF, 0xFF);
					}
					else if (needs_swz)
					{
						data[ASTCENC_SWZ_R] = blk->data_r[idx];
						data[ASTCENC_SWZ_G] = blk->data_g[idx];
						data[ASTCENC_SWZ_B] = blk->data_b[idx];
						data[ASTCENC_SWZ_A] = blk->data_a[idx];

						if (needs_z)
						{
							float xcoord = (data[0] * 2.0f) - 1.0f;
							float ycoord = (data[3] * 2.0f) - 1.0f;
							float zcoord = 1.0f - xcoord * xcoord - ycoord * ycoord;
							if (zcoord < 0.0f)
							{
								zcoord = 0.0f;
							}
							data[ASTCENC_SWZ_Z] = (astc::sqrt(zcoord) * 0.5f) + 0.5f;
						}

						vfloat4 color = vfloat4(data[swz.r], data[swz.g], data[swz.b], data[swz.a]);
						colori = float_to_int_rtn(min(color, 1.0f) * 255.0f);
					}
					else
					{
						vfloat4 color = blk->texel(idx);
						colori = float_to_int_rtn(min(color, 1.0f) * 255.0f);
					}

					colori = pack_low_bytes(colori);
					store_nbytes(colori, data8 + (4 * xsize * y) + (4 * x    ));

					idx++;
				}
				idx += x_nudge;
			}
			idx += y_nudge;
		}
	}
	else if (img.data_type == ASTCENC_TYPE_F16)
	{
		for (int z = z_start; z < z_end; z++)
		{
			// Fetch the image plane
			uint16_t* data16 = static_cast<uint16_t*>(img.data[z]);

			for (int y = y_start; y < y_end; y++)
			{
				for (int x = x_start; x < x_end; x++)
				{
					vint4 color;

					if (blk->data_r[idx] == std::numeric_limits<float>::quiet_NaN())
					{
						color = vint4(0xFFFF);
					}
					else if (needs_swz)
					{
						data[ASTCENC_SWZ_R] = blk->data_r[idx];
						data[ASTCENC_SWZ_G] = blk->data_g[idx];
						data[ASTCENC_SWZ_B] = blk->data_b[idx];
						data[ASTCENC_SWZ_A] = blk->data_a[idx];

						if (needs_z)
						{
							float xN = (data[0] * 2.0f) - 1.0f;
							float yN = (data[3] * 2.0f) - 1.0f;
							float zN = 1.0f - xN * xN - yN * yN;
							if (zN < 0.0f)
							{
								zN = 0.0f;
							}
							data[ASTCENC_SWZ_Z] = (astc::sqrt(zN) * 0.5f) + 0.5f;
						}

						vfloat4 colorf(data[swz.r], data[swz.g], data[swz.b], data[swz.a]);
						color = float_to_float16(colorf);
					}
					else
					{
						vfloat4 colorf = blk->texel(idx);
						color = float_to_float16(colorf);
					}

					data16[(4 * xsize * y) + (4 * x    )] = (uint16_t)color.lane<0>();
					data16[(4 * xsize * y) + (4 * x + 1)] = (uint16_t)color.lane<1>();
					data16[(4 * xsize * y) + (4 * x + 2)] = (uint16_t)color.lane<2>();
					data16[(4 * xsize * y) + (4 * x + 3)] = (uint16_t)color.lane<3>();

					idx++;
				}
				idx += x_nudge;
			}
			idx += y_nudge;
		}
	}
	else // if (img.data_type == ASTCENC_TYPE_F32)
	{
		assert(img.data_type == ASTCENC_TYPE_F32);

		for (int z = z_start; z < z_end; z++)
		{
			// Fetch the image plane
			float* data32 = static_cast<float*>(img.data[z]);

			for (int y = y_start; y < y_end; y++)
			{
				for (int x = x_start; x < x_end; x++)
				{
					vfloat4 color = blk->texel(idx);

					if (color.lane<0>() == std::numeric_limits<float>::quiet_NaN())
					{
						color = vfloat4(std::numeric_limits<float>::quiet_NaN());
					}
					else if (needs_swz)
					{
						data[ASTCENC_SWZ_R] = color.lane<0>();
						data[ASTCENC_SWZ_G] = color.lane<1>();
						data[ASTCENC_SWZ_B] = color.lane<2>();
						data[ASTCENC_SWZ_A] = color.lane<3>();

						if (needs_z)
						{
							float xN = (data[0] * 2.0f) - 1.0f;
							float yN = (data[3] * 2.0f) - 1.0f;
							float zN = 1.0f - xN * xN - yN * yN;
							if (zN < 0.0f)
							{
								zN = 0.0f;
							}
							data[ASTCENC_SWZ_Z] = (astc::sqrt(zN) * 0.5f) + 0.5f;
						}

						color = vfloat4(data[swz.r], data[swz.g], data[swz.b], data[swz.a]);
					}

					store(color, data32 + (4 * xsize * y) + (4 * x    ));

					idx++;
				}
				idx += x_nudge;
			}
			idx += y_nudge;
		}
	}
}
