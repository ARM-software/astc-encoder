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

astc_codec_image *alloc_image(
	int bitness,
	int xsize,
	int ysize,
	int zsize,
	int padding
) {
	int i, j;
	astc_codec_image *img = new astc_codec_image;
	img->xsize = xsize;
	img->ysize = ysize;
	img->zsize = zsize;
	img->padding = padding;

	img->input_averages = nullptr;
	img->input_variances = nullptr;
	img->input_alpha_averages = nullptr;

	img->linearize_srgb = 0;

	int exsize = xsize + 2 * padding;
	int eysize = ysize + 2 * padding;
	int ezsize = (zsize == 1) ? 1 : zsize + 2 * padding;

	assert(bitness == 8 || bitness == 16);
	if (bitness == 8)
	{
		img->data8 = new uint8_t **[ezsize];
		img->data8[0] = new uint8_t *[ezsize * eysize];
		img->data8[0][0] = new uint8_t[4 * ezsize * eysize * exsize];
		memset(img->data8[0][0], 0, 4 * ezsize * eysize * exsize);

		for (i = 1; i < ezsize; i++)
		{
			img->data8[i] = img->data8[0] + i * eysize;
			img->data8[i][0] = img->data8[0][0] + 4 * i * exsize * eysize;
		}

		for (i = 0; i < ezsize; i++)
		{
			for (j = 1; j < eysize; j++)
			{
				img->data8[i][j] = img->data8[i][0] + 4 * j * exsize;
			}
		}

		img->data16 = nullptr;
	}
	else if (bitness == 16)
	{
		img->data16 = new uint16_t **[ezsize];
		img->data16[0] = new uint16_t *[ezsize * eysize];
		img->data16[0][0] = new uint16_t[4 * ezsize * eysize * exsize];
		memset(img->data16[0][0], 0, 8 * ezsize * eysize * exsize);

		for (i = 1; i < ezsize; i++)
		{
			img->data16[i] = img->data16[0] + i * eysize;
			img->data16[i][0] = img->data16[0][0] + 4 * i * exsize * eysize;
		}

		for (i = 0; i < ezsize; i++)
		{
			for (j = 1; j < eysize; j++)
			{
				img->data16[i][j] = img->data16[i][0] + 4 * j * exsize;
			}
		}

		img->data8 = nullptr;
	}

	return img;
}

void free_image(astc_codec_image * img)
{
	if (img == nullptr)
	{
		return;
	}

	if (img->data8)
	{
		delete[] img->data8[0][0];
		delete[] img->data8[0];
		delete[] img->data8;
	}

	if (img->data16)
	{
		delete[] img->data16[0][0];
		delete[] img->data16[0];
		delete[] img->data16;
	}

	if (img->input_averages)
	{
		delete[] img->input_averages;
	}

	if (img->input_variances)
	{
		delete[] img->input_variances;
	}

	if (img->input_alpha_averages)
	{
		delete[] img->input_alpha_averages;
	}

	delete img;
}


// fill the padding area of the input-file buffer with clamp-to-edge data
// Done inefficiently, in that it will overwrite all the interior data at least once;
// this is not considered a problem, since this makes up a very small part of total
// running time.
void fill_image_padding_area(astc_codec_image * img)
{
	if (img->padding == 0)
	{
		return;
	}

	int x, y, z, i;
	int exsize = img->xsize + 2 * img->padding;
	int eysize = img->ysize + 2 * img->padding;
	int ezsize = (img->zsize == 1) ? 1 : (img->zsize + 2 * img->padding);

	int xmin = img->padding;
	int ymin = img->padding;
	int zmin = (img->zsize == 1) ? 0 : img->padding;
	int xmax = img->xsize + img->padding - 1;
	int ymax = img->ysize + img->padding - 1;
	int zmax = (img->zsize == 1) ? 0 : img->zsize + img->padding - 1;

	// This is a very simple implementation. Possible optimizations include:
	// * Testing if texel is outside the edge.
	// * Looping over texels that we know are outside the edge.
	if (img->data8)
	{
		for (z = 0; z < ezsize; z++)
		{
			int zc = MIN(MAX(z, zmin), zmax);
			for (y = 0; y < eysize; y++)
			{
				int yc = MIN(MAX(y, ymin), ymax);
				for (x = 0; x < exsize; x++)
				{
					int xc = MIN(MAX(x, xmin), xmax);
					for (i = 0; i < 4; i++)
					{
						img->data8[z][y][4 * x + i] = img->data8[zc][yc][4 * xc + i];
					}
				}
			}
		}
	}
	else if (img->data16)
	{
		for (z = 0; z < ezsize; z++)
		{
			int zc = MIN(MAX(z, zmin), zmax);
			for (y = 0; y < eysize; y++)
			{
				int yc = MIN(MAX(y, ymin), ymax);
				for (x = 0; x < exsize; x++)
				{
					int xc = MIN(MAX(x, xmin), xmax);
					for (i = 0; i < 4; i++)
					{
						img->data16[z][y][4 * x + i] = img->data16[zc][yc][4 * xc + i];
					}
				}
			}
		}
	}
}

int determine_image_channels(const astc_codec_image * img)
{
	int x, y, z;

	int xsize = img->xsize;
	int ysize = img->ysize;
	int zsize = img->zsize;
	// scan through the image data
	// to determine how many color channels the image has.

	int lum_mask;
	int alpha_mask;
	int alpha_mask_ref;
	if (img->data8)
	{
		alpha_mask_ref = 0xFF;
		alpha_mask = 0xFF;
		lum_mask = 0;
		for (z = 0; z < zsize; z++)
		{
			for (y = 0; y < ysize; y++)
			{
				for (x = 0; x < xsize; x++)
				{
					int r = img->data8[z][y][4 * x];
					int g = img->data8[z][y][4 * x + 1];
					int b = img->data8[z][y][4 * x + 2];
					int a = img->data8[z][y][4 * x + 3];
					lum_mask |= (r ^ g) | (r ^ b);
					alpha_mask &= a;
				}
			}
		}
	}
	else // (bitness == 16)
	{
		alpha_mask_ref = 0xFFFF;
		alpha_mask = 0xFFFF;
		lum_mask = 0;
		for (z = 0; z < zsize; z++)
		{
			for (y = 0; y < ysize; y++)
			{
				for (x = 0; x < xsize; x++)
				{
					int r = img->data16[z][y][4 * x];
					int g = img->data16[z][y][4 * x + 1];
					int b = img->data16[z][y][4 * x + 2];
					int a = img->data16[z][y][4 * x + 3];
					lum_mask |= (r ^ g) | (r ^ b);
					alpha_mask &= (a ^ 0xC3FF);	// a ^ 0xC3FF returns FFFF if and only if the input is 1.0
				}
			}
		}
	}

	int image_channels = 1 + (lum_mask == 0 ? 0 : 2) + (alpha_mask == alpha_mask_ref ? 0 : 1);

	return image_channels;
}

// initialize an astc_codec_image data structure from a 2D array of RGBA float*4
astc_codec_image* astc_img_from_floatx4_array(
	const float* image,
	int xsize,
	int ysize,
	int padding,
	int y_flip
) {
	astc_codec_image* astc_img = alloc_image(16, xsize, ysize, 1, padding);

	for (int y = 0; y < ysize; y++)
	{
		int y_dst = y + padding;
		int y_src = y_flip ? (ysize - y - 1) : y;
		const float* src = image + 4 * xsize * y_src;

		for (int x = 0; x < xsize; x++)
		{
			int x_dst = x + padding;
			astc_img->data16[0][y_dst][4 * x_dst]     = float_to_sf16(src[4 * x],     SF_NEARESTEVEN);
			astc_img->data16[0][y_dst][4 * x_dst + 1] = float_to_sf16(src[4 * x + 1], SF_NEARESTEVEN);
			astc_img->data16[0][y_dst][4 * x_dst + 2] = float_to_sf16(src[4 * x + 2], SF_NEARESTEVEN);
			astc_img->data16[0][y_dst][4 * x_dst + 3] = float_to_sf16(src[4 * x + 3], SF_NEARESTEVEN);
		}
	}

	fill_image_padding_area(astc_img);
	return astc_img;
}

// initialize an astc_codec_image data structure from a 2D array of UNORM8
astc_codec_image* astc_img_from_unorm8x4_array(
	const uint8_t* imageptr,
	int xsize,
	int ysize,
	int padding,
	int y_flip
) {
	astc_codec_image *astc_img = alloc_image(8, xsize, ysize, 1, padding);

	for (int y = 0; y < ysize; y++)
	{
		int y_dst = y + padding;
		int y_src = y_flip ? (ysize - y - 1) : y;
		const uint8_t* src = imageptr + 4 * xsize * y_src;

		for (int x = 0; x < xsize; x++)
		{
			int x_dst = x + padding;
			astc_img->data8[0][y_dst][4 * x_dst]     = src[4 * x];
			astc_img->data8[0][y_dst][4 * x_dst + 1] = src[4 * x + 1];
			astc_img->data8[0][y_dst][4 * x_dst + 2] = src[4 * x + 2];
			astc_img->data8[0][y_dst][4 * x_dst + 3] = src[4 * x + 3];
		}
	}

	fill_image_padding_area(astc_img);
	return astc_img;
}

// initialize a flattened array of float4 values from an ASTC codec image
// The returned array is allocated with new[] and must be deleted with delete[].
float* floatx4_array_from_astc_img(
	const astc_codec_image* img,
	int y_flip
) {
	int xsize = img->xsize;
	int ysize = img->ysize;

	float *buf = new float[4 * xsize * ysize];
	if (img->data8)
	{
		for (int y = 0; y < ysize; y++)
		{
			int ymod = y_flip ? ysize - y - 1 : y;
			const uint8_t* src = img->data8[0][ymod + img->padding] + (4 * img->padding);
			float* dst = buf + y * xsize * 4;
			for (int x = 0; x < xsize; x++)
			{
				dst[4 * x]     = src[4 * x]     * (1.0f / 255.0f);
				dst[4 * x + 1] = src[4 * x + 1] * (1.0f / 255.0f);
				dst[4 * x + 2] = src[4 * x + 2] * (1.0f / 255.0f);
				dst[4 * x + 3] = src[4 * x + 3] * (1.0f / 255.0f);
			}
		}
	}
	else
	{
		for (int y = 0; y < ysize; y++)
		{
			int ymod = y_flip ? ysize - y - 1 : y;
			const uint16_t *src = img->data16[0][ymod + img->padding] + (4 * img->padding);
			float *dst = buf + y * xsize * 4;
			for (int x = 0; x < xsize; x++)
			{
				dst[4 * x]     = sf16_to_float(src[4 * x]);
				dst[4 * x + 1] = sf16_to_float(src[4 * x + 1]);
				dst[4 * x + 2] = sf16_to_float(src[4 * x + 2]);
				dst[4 * x + 3] = sf16_to_float(src[4 * x + 3]);
			}
		}
	}

	return buf;
}

// initialize a flattened array of unorm8x4 values from an ASTC codec image
// The returned array is allocated with new[] and must be deleted with delete[].
uint8_t* unorm8x4_array_from_astc_img(
	const astc_codec_image* img,
	int y_flip
) {
	int xsize = img->xsize;
	int ysize = img->ysize;

	uint8_t* buf = new uint8_t[4 * xsize * ysize];

	if (img->data8)
	{
		for (int y = 0; y < ysize; y++)
		{
			int ymod = y_flip ? ysize-y-1 : y;
			const uint8_t* src = img->data8[0][ymod + img->padding] + (4 * img->padding);
			uint8_t* dst = buf + y * xsize * 4;
			for (int x = 0; x < xsize; x++)
			{
				dst[4 * x]     = src[4 * x];
				dst[4 * x + 1] = src[4 * x + 1];
				dst[4 * x + 2] = src[4 * x + 2];
				dst[4 * x + 3] = src[4 * x + 3];
			}
		}
	}
	else
	{
		for (int y = 0; y < ysize; y++)
		{
			int ymod = y_flip ? ysize-y-1 : y;
			const uint16_t* src = img->data16[0][ymod + img->padding] + (4 * img->padding);
			uint8_t* dst = buf + y * xsize * 4;
			for (int x = 0; x < xsize; x++)
			{
				dst[4 * x]     = (uint8_t)astc::flt2int_rtn(astc::clamp1f(sf16_to_float(src[4*x]))   * 255.0f);
				dst[4 * x + 1] = (uint8_t)astc::flt2int_rtn(astc::clamp1f(sf16_to_float(src[4*x+1])) * 255.0f);
				dst[4 * x + 2] = (uint8_t)astc::flt2int_rtn(astc::clamp1f(sf16_to_float(src[4*x+2])) * 255.0f);
				dst[4 * x + 3] = (uint8_t)astc::flt2int_rtn(astc::clamp1f(sf16_to_float(src[4*x+3])) * 255.0f);
			}
		}
	}

	return buf;
}
