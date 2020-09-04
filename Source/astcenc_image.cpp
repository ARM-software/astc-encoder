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
 * @brief Functions for creating in-memory ASTC image structures.
 */

#include <cassert>
#include <cstring>

#include "astcenc_internal.h"

// conversion functions between the LNS representation and the FP16 representation.
float float_to_lns(float p)
{
	if (astc::isnan(p) || p <= 1.0f / 67108864.0f)
	{
		// underflow or NaN value, return 0.
		// We count underflow if the input value is smaller than 2^-26.
		return 0;
	}

	if (fabsf(p) >= 65536.0f)
	{
		// overflow, return a +INF value
		return 65535;
	}

	int expo;
	float normfrac = frexpf(p, &expo);
	float p1;
	if (expo < -13)
	{
		// input number is smaller than 2^-14. In this case, multiply by 2^25.
		p1 = p * 33554432.0f;
		expo = 0;
	}
	else
	{
		expo += 14;
		p1 = (normfrac - 0.5f) * 4096.0f;
	}

	if (p1 < 384.0f)
		p1 *= 4.0f / 3.0f;
	else if (p1 <= 1408.0f)
		p1 += 128.0f;
	else
		p1 = (p1 + 512.0f) * (4.0f / 5.0f);

	p1 += expo * 2048.0f;
	return p1 + 1.0f;
}

uint16_t lns_to_sf16(uint16_t p)
{
	uint16_t mc = p & 0x7FF;
	uint16_t ec = p >> 11;
	uint16_t mt;
	if (mc < 512)
		mt = 3 * mc;
	else if (mc < 1536)
		mt = 4 * mc - 512;
	else
		mt = 5 * mc - 2048;

	uint16_t res = (ec << 10) | (mt >> 3);
	if (res >= 0x7BFF)
		res = 0x7BFF;
	return res;
}

// conversion function from 16-bit LDR value to FP16.
// note: for LDR interpolation, it is impossible to get a denormal result;
// this simplifies the conversion.
// FALSE; we can receive a very small UNORM16 through the constant-block.
uint16_t unorm16_to_sf16(uint16_t p)
{
	if (p == 0xFFFF)
		return 0x3C00;			// value of 1.0 .
	if (p < 4)
		return p << 8;

	int lz = clz32(p) - 16;
	p <<= (lz + 1);
	p >>= 6;
	p |= (14 - lz) << 10;
	return p;
}

void imageblock_initialize_deriv(
	const imageblock* pb,
	int pixelcount,
	float4* dptr
) {
	const float *fptr = pb->orig_data;
	for (int i = 0; i < pixelcount; i++)
	{
		// compute derivatives for RGB first
		if (pb->rgb_lns[i])
		{
			float r = MAX(fptr[0], 6e-5f);
			float g = MAX(fptr[1], 6e-5f);
			float b = MAX(fptr[2], 6e-5f);

			float rderiv = (float_to_lns(r * 1.05f) - float_to_lns(r)) / (r * 0.05f);
			float gderiv = (float_to_lns(g * 1.05f) - float_to_lns(g)) / (g * 0.05f);
			float bderiv = (float_to_lns(b * 1.05f) - float_to_lns(b)) / (b * 0.05f);

			// the derivative may not actually take values smaller than 1/32 or larger than 2^25;
			// if it does, we clamp it.
			if (rderiv < (1.0f / 32.0f))
			{
				rderiv = (1.0f / 32.0f);
			}
			else if (rderiv > 33554432.0f)
			{
				rderiv = 33554432.0f;
			}

			if (gderiv < (1.0f / 32.0f))
			{
				gderiv = (1.0f / 32.0f);
			}
			else if (gderiv > 33554432.0f)
			{
				gderiv = 33554432.0f;
			}

			if (bderiv < (1.0f / 32.0f))
			{
				bderiv = (1.0f / 32.0f);
			}
			else if (bderiv > 33554432.0f)
			{
				bderiv = 33554432.0f;
			}

			dptr->x = rderiv;
			dptr->y = gderiv;
			dptr->z = bderiv;
		}
		else
		{
			dptr->x = 65535.0f;
			dptr->y = 65535.0f;
			dptr->z = 65535.0f;
		}

		// then compute derivatives for Alpha
		if (pb->alpha_lns[i])
		{
			float a = MAX(fptr[3], 6e-5f);
			float aderiv = (float_to_lns(a * 1.05f) - float_to_lns(a)) / (a * 0.05f);
			// the derivative may not actually take values smaller than 1/32 or larger than 2^25;
			// if it does, we clamp it.
			if (aderiv < (1.0f / 32.0f))
			{
				aderiv = (1.0f / 32.0f);
			}
			else if (aderiv > 33554432.0f)
			{
				aderiv = 33554432.0f;
			}

			dptr->w = aderiv;
		}
		else
		{
			dptr->w = 65535.0f;
		}

		fptr += 4;
		dptr += 1;
	}
}

// helper function to initialize the work-data from the orig-data
void imageblock_initialize_work_from_orig(
	imageblock* pb,
	int pixelcount
) {
	float *fptr = pb->orig_data;

	for (int i = 0; i < pixelcount; i++)
	{
		if (pb->rgb_lns[i])
		{
			pb->data_r[i] = float_to_lns(fptr[0]);
			pb->data_g[i] = float_to_lns(fptr[1]);
			pb->data_b[i] = float_to_lns(fptr[2]);
		}
		else
		{
			pb->data_r[i] = fptr[0] * 65535.0f;
			pb->data_g[i] = fptr[1] * 65535.0f;
			pb->data_b[i] = fptr[2] * 65535.0f;
		}

		if (pb->alpha_lns[i])
		{
			pb->data_a[i] = float_to_lns(fptr[3]);
		}
		else
		{
			pb->data_a[i] = fptr[3] * 65535.0f;
		}

		fptr += 4;
	}
}

// helper function to initialize the orig-data from the work-data
void imageblock_initialize_orig_from_work(
	imageblock* pb,
	int pixelcount
) {
	float *fptr = pb->orig_data;

	for (int i = 0; i < pixelcount; i++)
	{
		if (pb->rgb_lns[i])
		{
			fptr[0] = sf16_to_float(lns_to_sf16((uint16_t)pb->data_r[i]));
			fptr[1] = sf16_to_float(lns_to_sf16((uint16_t)pb->data_g[i]));
			fptr[2] = sf16_to_float(lns_to_sf16((uint16_t)pb->data_b[i]));
		}
		else
		{
			fptr[0] = sf16_to_float(unorm16_to_sf16((uint16_t)pb->data_r[i]));
			fptr[1] = sf16_to_float(unorm16_to_sf16((uint16_t)pb->data_g[i]));
			fptr[2] = sf16_to_float(unorm16_to_sf16((uint16_t)pb->data_b[i]));
		}

		if (pb->alpha_lns[i])
		{
			fptr[3] = sf16_to_float(lns_to_sf16((uint16_t)pb->data_a[i]));
		}
		else
		{
			fptr[3] = sf16_to_float(unorm16_to_sf16((uint16_t)pb->data_a[i]));
		}

		fptr += 4;
	}
}

// fetch an imageblock from the input file.
void fetch_imageblock(
	astcenc_profile decode_mode,
	const astcenc_image& img,
	imageblock* pb,	// picture-block to initialize with image data
	const block_size_descriptor* bsd,
	// position in texture.
	int xpos,
	int ypos,
	int zpos,
	astcenc_swizzle swz
) {
	float *fptr = pb->orig_data;
	int xsize = img.dim_x + 2 * img.dim_pad;
	int ysize = img.dim_y + 2 * img.dim_pad;
	int zsize = (img.dim_z == 1) ? 1 : img.dim_z + 2 * img.dim_pad;

	pb->xpos = xpos;
	pb->ypos = ypos;
	pb->zpos = zpos;

	xpos += img.dim_pad;
	ypos += img.dim_pad;
	if (img.dim_z > 1)
	{
		zpos += img.dim_pad;
	}

	float data[6];
	data[4] = 0;
	data[5] = 1;

	if (img.data_type == ASTCENC_TYPE_U8)
	{
		uint8_t*** data8 = static_cast<uint8_t***>(img.data);
		for (int z = 0; z < bsd->zdim; z++)
		{
			for (int y = 0; y < bsd->ydim; y++)
			{
				for (int x = 0; x < bsd->xdim; x++)
				{
					int xi = xpos + x;
					int yi = ypos + y;
					int zi = zpos + z;
					// clamp XY coordinates to the picture.
					if (xi < 0)
						xi = 0;
					if (yi < 0)
						yi = 0;
					if (zi < 0)
						zi = 0;
					if (xi >= xsize)
						xi = xsize - 1;
					if (yi >= ysize)
						yi = ysize - 1;
					if (zi >= zsize)
						zi = zsize - 1;

					int r = data8[zi][yi][4 * xi];
					int g = data8[zi][yi][4 * xi + 1];
					int b = data8[zi][yi][4 * xi + 2];
					int a = data8[zi][yi][4 * xi + 3];

					data[0] = r / 255.0f;
					data[1] = g / 255.0f;
					data[2] = b / 255.0f;
					data[3] = a / 255.0f;

					fptr[0] = data[swz.r];
					fptr[1] = data[swz.g];
					fptr[2] = data[swz.b];
					fptr[3] = data[swz.a];

					fptr += 4;
				}
			}
		}
	}
	else // if (img.data_type == ASTCENC_TYPE_F16)
	{
		assert(img.data_type == ASTCENC_TYPE_F16);
		uint16_t*** data16 = static_cast<uint16_t***>(img.data);
		for (int z = 0; z < bsd->zdim; z++)
		{
			for (int y = 0; y < bsd->ydim; y++)
			{
				for (int x = 0; x < bsd->xdim; x++)
				{
					int xi = xpos + x;
					int yi = ypos + y;
					int zi = zpos + z;
					// clamp XY coordinates to the picture.
					if (xi < 0)
						xi = 0;
					if (yi < 0)
						yi = 0;
					if (zi < 0)
						zi = 0;
					if (xi >= xsize)
						xi = xsize - 1;
					if (yi >= ysize)
						yi = ysize - 1;
					if (zi >= ysize)
						zi = zsize - 1;

					int r = data16[zi][yi][4 * xi];
					int g = data16[zi][yi][4 * xi + 1];
					int b = data16[zi][yi][4 * xi + 2];
					int a = data16[zi][yi][4 * xi + 3];

					float rf = sf16_to_float(r);
					float gf = sf16_to_float(g);
					float bf = sf16_to_float(b);
					float af = sf16_to_float(a);

					// equalize the color components somewhat, and get rid of negative values.
					rf = MAX(rf, 1e-8f);
					gf = MAX(gf, 1e-8f);
					bf = MAX(bf, 1e-8f);
					af = MAX(af, 1e-8f);

					data[0] = rf;
					data[1] = gf;
					data[2] = bf;
					data[3] = af;

					fptr[0] = data[swz.r];
					fptr[1] = data[swz.g];
					fptr[2] = data[swz.b];
					fptr[3] = data[swz.a];
					fptr += 4;
				}
			}
		}
	}

	int rgb_lns = (decode_mode == ASTCENC_PRF_HDR) || (decode_mode == ASTCENC_PRF_HDR_RGB_LDR_A);
	int alpha_lns = decode_mode == ASTCENC_PRF_HDR;

	// impose the choice on every pixel when encoding.
	for (int i = 0; i < bsd->texel_count; i++)
	{
		pb->rgb_lns[i] = rgb_lns;
		pb->alpha_lns[i] = alpha_lns;
		pb->nan_texel[i] = 0;
	}

	imageblock_initialize_work_from_orig(pb,bsd->texel_count);
	update_imageblock_flags(pb, bsd->xdim, bsd->ydim, bsd->zdim);
}

void write_imageblock(
	astcenc_image& img,
	const imageblock* pb,	// picture-block to initialize with image data. We assume that orig_data is valid.
	const block_size_descriptor* bsd,
	// position to write the block to
	int xpos,
	int ypos,
	int zpos,
	astcenc_swizzle swz
) {
	const float *fptr = pb->orig_data;
	const uint8_t *nptr = pb->nan_texel;
	int xsize = img.dim_x;
	int ysize = img.dim_y;
	int zsize = img.dim_z;

	float data[7];
	data[4] = 0.0f;
	data[5] = 1.0f;

	if (img.data_type == ASTCENC_TYPE_U8)
	{
		uint8_t*** data8 = static_cast<uint8_t***>(img.data);
		for (int z = 0; z < bsd->zdim; z++)
		{
			for (int y = 0; y < bsd->ydim; y++)
			{
				for (int x = 0; x < bsd->xdim; x++)
				{
					int xi = xpos + x;
					int yi = ypos + y;
					int zi = zpos + z;

					if (xi >= 0 && yi >= 0 && zi >= 0 && xi < xsize && yi < ysize && zi < zsize)
					{
						if (*nptr)
						{
							// NaN-pixel, but we can't display it. Display purple instead.
							data8[zi][yi][4 * xi] = 0xFF;
							data8[zi][yi][4 * xi + 1] = 0x00;
							data8[zi][yi][4 * xi + 2] = 0xFF;
							data8[zi][yi][4 * xi + 3] = 0xFF;
						}
						else
						{
							data[0] = fptr[0];
							data[1] = fptr[1];
							data[2] = fptr[2];
							data[3] = fptr[3];

							float xcoord = (data[0] * 2.0f) - 1.0f;
							float ycoord = (data[3] * 2.0f) - 1.0f;
							float zcoord = 1.0f - xcoord * xcoord - ycoord * ycoord;
							if (zcoord < 0.0f)
							{
								zcoord = 0.0f;
							}
							data[6] = (astc::sqrt(zcoord) * 0.5f) + 0.5f;

							// clamp to [0,1]
							if (data[0] > 1.0f)
							{
								data[0] = 1.0f;
							}

							if (data[1] > 1.0f)
							{
								data[1] = 1.0f;
							}

							if (data[2] > 1.0f)
							{
								data[2] = 1.0f;
							}

							if (data[3] > 1.0f)
							{
								data[3] = 1.0f;
							}

							// pack the data
							int ri = astc::flt2int_rtn(data[swz.r] * 255.0f);
							int gi = astc::flt2int_rtn(data[swz.g] * 255.0f);
							int bi = astc::flt2int_rtn(data[swz.b] * 255.0f);
							int ai = astc::flt2int_rtn(data[swz.a] * 255.0f);

							data8[zi][yi][4 * xi] = ri;
							data8[zi][yi][4 * xi + 1] = gi;
							data8[zi][yi][4 * xi + 2] = bi;
							data8[zi][yi][4 * xi + 3] = ai;
						}
					}
					fptr += 4;
					nptr++;
				}
			}
		}
	}
	else // if (img.data_type == ASTCENC_TYPE_F16)
	{
		assert(img.data_type == ASTCENC_TYPE_F16);
		uint16_t*** data16 = static_cast<uint16_t***>(img.data);
		for (int z = 0; z < bsd->zdim; z++)
		{
			for (int y = 0; y < bsd->ydim; y++)
			{
				for (int x = 0; x < bsd->xdim; x++)
				{
					int xi = xpos + x;
					int yi = ypos + y;
					int zi = zpos + z;

					if (xi >= 0 && yi >= 0 && zi >= 0 && xi < xsize && yi < ysize && zi < zsize)
					{
						if (*nptr)
						{
							data16[zi][yi][4 * xi] = 0xFFFF;
							data16[zi][yi][4 * xi + 1] = 0xFFFF;
							data16[zi][yi][4 * xi + 2] = 0xFFFF;
							data16[zi][yi][4 * xi + 3] = 0xFFFF;
						}

						else
						{
							data[0] = fptr[0];
							data[1] = fptr[1];
							data[2] = fptr[2];
							data[3] = fptr[3];

							float xN = (data[0] * 2.0f) - 1.0f;
							float yN = (data[3] * 2.0f) - 1.0f;
							float zN = 1.0f - xN * xN - yN * yN;
							if (zN < 0.0f)
							{
								zN = 0.0f;
							}
							data[6] = (astc::sqrt(zN) * 0.5f) + 0.5f;

							int r = float_to_sf16(data[swz.r], SF_NEARESTEVEN);
							int g = float_to_sf16(data[swz.g], SF_NEARESTEVEN);
							int b = float_to_sf16(data[swz.b], SF_NEARESTEVEN);
							int a = float_to_sf16(data[swz.a], SF_NEARESTEVEN);
							data16[zi][yi][4 * xi] = r;
							data16[zi][yi][4 * xi + 1] = g;
							data16[zi][yi][4 * xi + 2] = b;
							data16[zi][yi][4 * xi + 3] = a;
						}
					}
					fptr += 4;
					nptr++;
				}
			}
		}
	}
}

/*
   For an imageblock, update its flags.
   The updating is done based on data, not orig_data.
*/
void update_imageblock_flags(
	imageblock* pb,
	int xdim,
	int ydim,
	int zdim
) {
	float red_min = 1e38f, red_max = -1e38f;
	float green_min = 1e38f, green_max = -1e38f;
	float blue_min = 1e38f, blue_max = -1e38f;
	float alpha_min = 1e38f, alpha_max = -1e38f;

	int texels_per_block = xdim * ydim * zdim;

	int grayscale = 1;

	for (int i = 0; i < texels_per_block; i++)
	{
		float red = pb->data_r[i];
		float green = pb->data_g[i];
		float blue = pb->data_b[i];
		float alpha = pb->data_a[i];
		if (red < red_min)
			red_min = red;
		if (red > red_max)
			red_max = red;
		if (green < green_min)
			green_min = green;
		if (green > green_max)
			green_max = green;
		if (blue < blue_min)
			blue_min = blue;
		if (blue > blue_max)
			blue_max = blue;
		if (alpha < alpha_min)
			alpha_min = alpha;
		if (alpha > alpha_max)
			alpha_max = alpha;

		if (grayscale == 1 && (red != green || red != blue))
		{
			grayscale = 0;
		}
	}

	pb->red_min = red_min;
	pb->red_max = red_max;
	pb->green_min = green_min;
	pb->green_max = green_max;
	pb->blue_min = blue_min;
	pb->blue_max = blue_max;
	pb->alpha_min = alpha_min;
	pb->alpha_max = alpha_max;
	pb->grayscale = grayscale;
}
