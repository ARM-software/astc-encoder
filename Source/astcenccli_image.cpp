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

#include "astcenccli_internal.h"

astcenc_image *alloc_image(
	unsigned int bitness,
	unsigned int dim_x,
	unsigned int dim_y,
	unsigned int dim_z
) {
	astcenc_image *img = new astcenc_image;
	img->dim_x = dim_x;
	img->dim_y = dim_y;
	img->dim_z = dim_z;

	if (bitness == 8)
	{
		uint8_t*** data8 = new uint8_t **[dim_z];
		data8[0] = new uint8_t *[dim_z * dim_y];
		data8[0][0] = new uint8_t[4 * dim_z * dim_y * dim_x];
		memset(data8[0][0], 0, 4 * dim_z * dim_y * dim_x);

		for (unsigned int z = 1; z < dim_z; z++)
		{
			data8[z] = data8[0] + z * dim_y;
			data8[z][0] = data8[0][0] + 4 * z * dim_x * dim_y;
		}

		for (unsigned int z = 0; z < dim_z; z++)
		{
			for (unsigned int y = 1; y < dim_y; y++)
			{
				data8[z][y] = data8[z][0] + 4 * y * dim_x;
			}
		}

		img->data_type = ASTCENC_TYPE_U8;
		img->data = static_cast<void*>(data8);
	}
	else if (bitness == 16)
	{
		uint16_t*** data16 = new uint16_t **[dim_z];
		data16[0] = new uint16_t *[dim_z * dim_y];
		data16[0][0] = new uint16_t[4 * dim_z * dim_y * dim_x];
		memset(data16[0][0], 0, 8 * dim_z * dim_y * dim_x);

		for (unsigned int z = 1; z < dim_z; z++)
		{
			data16[z] = data16[0] + z * dim_y;
			data16[z][0] = data16[0][0] + 4 * z * dim_x * dim_y;
		}

		for (unsigned int z = 0; z < dim_z; z++)
		{
			for (unsigned int y = 1; y < dim_y; y++)
			{
				data16[z][y] = data16[z][0] + 4 * y * dim_x;
			}
		}

		img->data_type = ASTCENC_TYPE_F16;
		img->data = static_cast<void*>(data16);
	}
	else // if (bitness == 32)
	{
		assert(bitness == 32);
		float*** data32 = new float**[dim_z];
		data32[0] = new float*[dim_z * dim_y];
		data32[0][0] = new float[4 * dim_z * dim_y * dim_x];
		memset(data32[0][0], 0, 8 * dim_z * dim_y * dim_x);

		for (unsigned int z = 1; z < dim_z; z++)
		{
			data32[z] = data32[0] + z * dim_y;
			data32[z][0] = data32[0][0] + 4 * z * dim_x * dim_y;
		}

		for (unsigned int z = 0; z < dim_z; z++)
		{
			for (unsigned int y = 1; y < dim_y; y++)
			{
				data32[z][y] = data32[z][0] + 4 * y * dim_x;
			}
		}

		img->data_type = ASTCENC_TYPE_F32;
		img->data = static_cast<void*>(data32);
	}

	return img;
}

void free_image(astcenc_image * img)
{
	if (img == nullptr)
	{
		return;
	}

	if (img->data_type == ASTCENC_TYPE_U8)
	{
		uint8_t*** data8 = static_cast<uint8_t***>(img->data);
		delete[] data8[0][0];
		delete[] data8[0];
		delete[] data8;
	}
	else if (img->data_type == ASTCENC_TYPE_F16)
	{
		uint16_t*** data16 = static_cast<uint16_t***>(img->data);
		delete[] data16[0][0];
		delete[] data16[0];
		delete[] data16;
	}
	else // if (img->data_type == ASTCENC_TYPE_F32)
	{
		assert(img->data_type == ASTCENC_TYPE_F32);
		float*** data32 = static_cast<float***>(img->data);
		delete[] data32[0][0];
		delete[] data32[0];
		delete[] data32;
	}

	delete img;
}

int determine_image_channels(const astcenc_image * img)
{
	unsigned int dim_x = img->dim_x;
	unsigned int dim_y = img->dim_y;
	unsigned int dim_z = img->dim_z;

	// scan through the image data
	// to determine how many color channels the image has.

	bool is_luma = true;
	bool has_alpha = false;

	if (img->data_type == ASTCENC_TYPE_U8)
	{
		uint8_t*** data8 = static_cast<uint8_t***>(img->data);

		for (unsigned int z = 0; z < dim_z; z++)
		{
			for (unsigned int y = 0; y < dim_y; y++)
			{
				for (unsigned int x = 0; x < dim_x; x++)
				{
					int r = data8[z][y][4 * x];
					int g = data8[z][y][4 * x + 1];
					int b = data8[z][y][4 * x + 2];
					int a = data8[z][y][4 * x + 3];

					is_luma = is_luma && (r == g) && (r == b);
					has_alpha = has_alpha || (a != 0xFF);
				}
			}
		}
	}
	else if (img->data_type == ASTCENC_TYPE_F16)
	{
		uint16_t*** data16 = static_cast<uint16_t***>(img->data);
		for (unsigned int z = 0; z < dim_z; z++)
		{
			for (unsigned int y = 0; y < dim_y; y++)
			{
				for (unsigned int x = 0; x < dim_x; x++)
				{
					int r = data16[z][y][4 * x];
					int g = data16[z][y][4 * x + 1];
					int b = data16[z][y][4 * x + 2];
					int a = data16[z][y][4 * x + 3];

					is_luma = is_luma && (r == g) && (r == b);
					has_alpha = has_alpha || ((a ^ 0xC3FF) != 0xFFFF);
					// a ^ 0xC3FF returns FFFF if and only if the input is 1.0
				}
			}
		}
	}
	else // if (img->data_type == ASTCENC_TYPE_F32)
	{
		assert(img->data_type == ASTCENC_TYPE_F32);
		float*** data32 = static_cast<float***>(img->data);
		for (unsigned int z = 0; z < dim_z; z++)
		{
			for (unsigned int y = 0; y < dim_y; y++)
			{
				for (unsigned int x = 0; x < dim_x; x++)
				{
					float r = data32[z][y][4 * x];
					float g = data32[z][y][4 * x + 1];
					float b = data32[z][y][4 * x + 2];
					float a = data32[z][y][4 * x + 3];

					is_luma = is_luma && (r == g) && (r == b);
					has_alpha = has_alpha || (a != 1.0f);
				}
			}
		}
	}

	int image_channels = 1 + (is_luma == 0 ? 0 : 2) + (has_alpha ? 0 : 1);
	return image_channels;
}

// initialize an astcenc_image data structure from a 2D array of RGBA float*4
astcenc_image* astc_img_from_floatx4_array(
	const float* data,
	unsigned int dim_x,
	unsigned int dim_y,
	bool y_flip
) {
	// TODO: Make this 32 to use direct passthough as float
	astcenc_image* img = alloc_image(16, dim_x, dim_y, 1);

	for (unsigned int y = 0; y < dim_y; y++)
	{
#if 0
		float*** data32 = static_cast<float***>(img->data);
		unsigned int y_src = y_flip ? (dim_y - y - 1) : y;
		const float* src = data + 4 * dim_y * y_src;

		for (unsigned int x = 0; x < dim_x; x++)
		{
			data32[0][y][4 * x    ] = src[4 * x    ];
			data32[0][y][4 * x + 1] = src[4 * x + 1];
			data32[0][y][4 * x + 2] = src[4 * x + 2];
			data32[0][y][4 * x + 3] = src[4 * x + 3];
		}
#else
		uint16_t*** data16 = static_cast<uint16_t***>(img->data);
		unsigned int y_src = y_flip ? (dim_y - y - 1) : y;
		const float* src = data + 4 * dim_y * y_src;

		for (unsigned int x = 0; x < dim_x; x++)
		{
			data16[0][y][4 * x    ] = float_to_sf16(src[4 * x    ], SF_NEARESTEVEN);
			data16[0][y][4 * x + 1] = float_to_sf16(src[4 * x + 1], SF_NEARESTEVEN);
			data16[0][y][4 * x + 2] = float_to_sf16(src[4 * x + 2], SF_NEARESTEVEN);
			data16[0][y][4 * x + 3] = float_to_sf16(src[4 * x + 3], SF_NEARESTEVEN);
		}
#endif
	}

	return img;
}

// initialize an astcenc_image data structure from a 2D array of UNORM8
astcenc_image* astc_img_from_unorm8x4_array(
	const uint8_t* data,
	unsigned int dim_x,
	unsigned int dim_y,
	bool y_flip
) {
	astcenc_image* img = alloc_image(8, dim_x, dim_y, 1);

	for (unsigned int y = 0; y < dim_y; y++)
	{
		uint8_t*** data8 = static_cast<uint8_t***>(img->data);
		unsigned int y_src = y_flip ? (dim_y - y - 1) : y;
		const uint8_t* src = data + 4 * dim_x * y_src;

		for (unsigned int x = 0; x < dim_x; x++)
		{
			data8[0][y][4 * x    ] = src[4 * x    ];
			data8[0][y][4 * x + 1] = src[4 * x + 1];
			data8[0][y][4 * x + 2] = src[4 * x + 2];
			data8[0][y][4 * x + 3] = src[4 * x + 3];
		}
	}

	return img;
}

// initialize a flattened array of float4 values from an ASTC codec image
// The returned array is allocated with new[] and must be deleted with delete[].
float* floatx4_array_from_astc_img(
	const astcenc_image* img,
	bool y_flip
) {
	unsigned int dim_x = img->dim_x;
	unsigned int dim_y = img->dim_y;
	float *buf = new float[4 * dim_x * dim_y];

	if (img->data_type == ASTCENC_TYPE_U8)
	{
		uint8_t*** data8 = static_cast<uint8_t***>(img->data);
		for (unsigned int y = 0; y < dim_y; y++)
		{
			unsigned int ymod = y_flip ? dim_y - y - 1 : y;
			const uint8_t* src = data8[0][ymod];
			float* dst = buf + y * dim_x * 4;

			for (unsigned int x = 0; x < dim_x; x++)
			{
				dst[4 * x    ] = src[4 * x    ] * (1.0f / 255.0f);
				dst[4 * x + 1] = src[4 * x + 1] * (1.0f / 255.0f);
				dst[4 * x + 2] = src[4 * x + 2] * (1.0f / 255.0f);
				dst[4 * x + 3] = src[4 * x + 3] * (1.0f / 255.0f);
			}
		}
	}
	else if (img->data_type == ASTCENC_TYPE_F16)
	{
		uint16_t*** data16 = static_cast<uint16_t***>(img->data);
		for (unsigned int y = 0; y < dim_y; y++)
		{
			unsigned int ymod = y_flip ? dim_y - y - 1 : y;
			const uint16_t *src = data16[0][ymod];
			float *dst = buf + y * dim_x * 4;

			for (unsigned int x = 0; x < dim_x; x++)
			{
				dst[4 * x    ] = sf16_to_float(src[4 * x    ]);
				dst[4 * x + 1] = sf16_to_float(src[4 * x + 1]);
				dst[4 * x + 2] = sf16_to_float(src[4 * x + 2]);
				dst[4 * x + 3] = sf16_to_float(src[4 * x + 3]);
			}
		}
	}
	else // if (img->data_type == ASTCENC_TYPE_F32)
	{
		assert(img->data_type == ASTCENC_TYPE_F32);
		float*** data32 = static_cast<float***>(img->data);
		for (unsigned int y = 0; y < dim_y; y++)
		{
			unsigned int ymod = y_flip ? dim_y - y - 1 : y;
			const float *src = data32[0][ymod];
			float *dst = buf + y * dim_x * 4;

			for (unsigned int x = 0; x < dim_x; x++)
			{
				dst[4 * x    ] = src[4 * x    ];
				dst[4 * x + 1] = src[4 * x + 1];
				dst[4 * x + 2] = src[4 * x + 2];
				dst[4 * x + 3] = src[4 * x + 3];
			}
		}
	}

	return buf;
}

// initialize a flattened array of unorm8x4 values from an ASTC codec image
// The returned array is allocated with new[] and must be deleted with delete[].
uint8_t* unorm8x4_array_from_astc_img(
	const astcenc_image* img,
	bool y_flip
) {
	unsigned int dim_x = img->dim_x;
	unsigned int dim_y = img->dim_y;
	uint8_t* buf = new uint8_t[4 * dim_x * dim_y];

	if (img->data_type == ASTCENC_TYPE_U8)
	{
		uint8_t*** data8 = static_cast<uint8_t***>(img->data);
		for (unsigned int y = 0; y < dim_y; y++)
		{
			unsigned int ymod = y_flip ? dim_y - y - 1 : y;
			const uint8_t* src = data8[0][ymod];
			uint8_t* dst = buf + y * dim_x * 4;

			for (unsigned int x = 0; x < dim_x; x++)
			{
				dst[4 * x    ] = src[4 * x    ];
				dst[4 * x + 1] = src[4 * x + 1];
				dst[4 * x + 2] = src[4 * x + 2];
				dst[4 * x + 3] = src[4 * x + 3];
			}
		}
	}
	else if (img->data_type == ASTCENC_TYPE_F16)
	{
		uint16_t*** data16 = static_cast<uint16_t***>(img->data);
		for (unsigned int y = 0; y < dim_y; y++)
		{
			unsigned int ymod = y_flip ? dim_y - y - 1 : y;
			const uint16_t* src = data16[0][ymod];
			uint8_t* dst = buf + y * dim_x * 4;

			for (unsigned int x = 0; x < dim_x; x++)
			{
				dst[4 * x   ]  = (uint8_t)astc::flt2int_rtn(astc::clamp1f(sf16_to_float(src[4 * x    ])) * 255.0f);
				dst[4 * x + 1] = (uint8_t)astc::flt2int_rtn(astc::clamp1f(sf16_to_float(src[4 * x + 1])) * 255.0f);
				dst[4 * x + 2] = (uint8_t)astc::flt2int_rtn(astc::clamp1f(sf16_to_float(src[4 * x + 2])) * 255.0f);
				dst[4 * x + 3] = (uint8_t)astc::flt2int_rtn(astc::clamp1f(sf16_to_float(src[4 * x + 3])) * 255.0f);
			}
		}
	}
	else // if (img->data_type == ASTCENC_TYPE_F32)
	{
		assert(img->data_type == ASTCENC_TYPE_F32);
		float*** data32 = static_cast<float***>(img->data);
		for (unsigned int y = 0; y < dim_y; y++)
		{
			unsigned int ymod = y_flip ? dim_y - y - 1 : y;
			const float* src = data32[0][ymod];
			uint8_t* dst = buf + y * dim_x * 4;

			for (unsigned int x = 0; x < dim_x; x++)
			{
				dst[4 * x    ] = (uint8_t)astc::flt2int_rtn(astc::clamp1f(src[4 * x    ]) * 255.0f);
				dst[4 * x + 1] = (uint8_t)astc::flt2int_rtn(astc::clamp1f(src[4 * x + 1]) * 255.0f);
				dst[4 * x + 2] = (uint8_t)astc::flt2int_rtn(astc::clamp1f(src[4 * x + 2]) * 255.0f);
				dst[4 * x + 3] = (uint8_t)astc::flt2int_rtn(astc::clamp1f(src[4 * x + 3]) * 255.0f);
			}
		}
	}

	return buf;
}
