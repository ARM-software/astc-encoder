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
 * @brief Functions for loading/storing ASTC compressed images.
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "astcenccli_internal.h"

// on windows/msvc, compile stb and tinyexr together with this file;
// on other systems, use makefile to compile them separately.
#ifdef _MSC_VER
	#define STB_IMAGE_IMPLEMENTATION
	#define STB_IMAGE_WRITE_IMPLEMENTATION
	#define STBI_MSC_SECURE_CRT
	#define TINYEXR_IMPLEMENTATION
	#define STBI_NO_GIF
	#define STBI_NO_PIC
	#define STBI_NO_PNM
	#define STBI_NO_PSD
#else
	#define STBI_HEADER_FILE_ONLY
#endif

#include "stb_image.h"
#include "stb_image_write.h"
#include "tinyexr.h"

/*******************************************************************
Image load and store through the stb_iamge and tinyexr libraries
*******************************************************************/

static astcenc_image* load_image_with_tinyexr(
	const char* filename,
	unsigned int dim_pad,
	bool y_flip,
	bool& is_hdr,
	unsigned int& num_components
) {
	int dim_x, dim_y;
	float* image;
	const char* err;

	int load_res = LoadEXR(&image, &dim_x, &dim_y, filename, &err);
	if (load_res != TINYEXR_SUCCESS)
	{
		printf("ERROR: Failed to load image %s (%s)\n", filename, err);
		free((void*)err);
		return nullptr;
	}

	astcenc_image* res_img = astc_img_from_floatx4_array(image, dim_x, dim_y, dim_pad, y_flip);
	free(image);

	is_hdr = true;
	num_components = 4;
	return res_img;
}

static astcenc_image* load_image_with_stb(
	const char* filename,
	unsigned int dim_pad,
	bool y_flip,
	bool& is_hdr,
	unsigned int& num_components
) {
	int dim_x, dim_y;

	if (stbi_is_hdr(filename))
	{
		float* data = stbi_loadf(filename, &dim_x, &dim_y, nullptr, STBI_rgb_alpha);
		if (data)
		{
			astcenc_image* img = astc_img_from_floatx4_array(data, dim_x, dim_y, dim_pad, y_flip);
			stbi_image_free(data);
			is_hdr = true;
			num_components = 4;
			return img;
		}
	}
	else
	{
		uint8_t* data = stbi_load(filename, &dim_x, &dim_y, nullptr, STBI_rgb_alpha);
		if (data)
		{
			astcenc_image* img = astc_img_from_unorm8x4_array(data, dim_x, dim_y, dim_pad, y_flip);
			stbi_image_free(data);
			is_hdr = false;
			num_components = 4;
			return img;
		}
	}

	printf("ERROR: Failed to load image %s (%s)\n", filename, stbi_failure_reason());
	return nullptr;
}

static int store_exr_image_with_tinyexr(
	const astcenc_image* img,
	const char* filename,
	int y_flip
) {
	float *buf = floatx4_array_from_astc_img(img, y_flip);
	int res = SaveEXR(buf, img->dim_x, img->dim_y, 4, 1, filename, nullptr);
	delete[] buf;
	return (res == 0) ? 4 : res;
}

static int store_png_image_with_stb(
	const astcenc_image* img,
	const char* filename,
	int y_flip
) {
	uint8_t* buf = unorm8x4_array_from_astc_img(img, y_flip);
	int res = stbi_write_png(filename, img->dim_x, img->dim_y, 4, buf, img->dim_x * 4);
	delete[] buf;
	return (res == 0) ? -1 : 4;
}

static int store_tga_image_with_stb(
	const astcenc_image* img,
	const char* filename,
	int y_flip
) {
	uint8_t* buf = unorm8x4_array_from_astc_img(img, y_flip);
	int res = stbi_write_tga(filename, img->dim_x, img->dim_y, 4, buf);
	delete[] buf;
	return (res == 0) ? -1 : 4;
}

static int store_bmp_image_with_stb(
	const astcenc_image* img,
	const char* filename,
	int y_flip
) {
	uint8_t* buf = unorm8x4_array_from_astc_img(img, y_flip);
	int res = stbi_write_bmp(filename, img->dim_x, img->dim_y, 4, buf);
	delete[] buf;
	return (res == 0) ? -1 : 4;
}

/*********************************************************************
Native Load and store of KTX and DDS file formats.

Unlike "regular" 2D image formats, which are mostly supported
through stb_image and tinyexr, these formats are supported directly;
this involves a relatively large number of pixel formats.

The following restrictions apply to loading of these file formats:
 * Only uncompressed data supported
 * Only first mipmap in mipmap pyramid supported
 * KTX: Cube-map arrays are not supported
*********************************************************************/
enum scanline_copy_method
{
	R8_TO_RGBA8,
	RG8_TO_RGBA8,
	RGB8_TO_RGBA8,
	RGBA8_TO_RGBA8,
	BGR8_TO_RGBA8,
	BGRA8_TO_RGBA8,
	L8_TO_RGBA8,
	LA8_TO_RGBA8,

	RGBX8_TO_RGBA8,
	BGRX8_TO_RGBA8,

	R16_TO_RGBA16F,
	RG16_TO_RGBA16F,
	RGB16_TO_RGBA16F,
	RGBA16_TO_RGBA16F,
	BGR16_TO_RGBA16F,
	BGRA16_TO_RGBA16F,
	L16_TO_RGBA16F,
	LA16_TO_RGBA16F,

	R16F_TO_RGBA16F,
	RG16F_TO_RGBA16F,
	RGB16F_TO_RGBA16F,
	RGBA16F_TO_RGBA16F,
	BGR16F_TO_RGBA16F,
	BGRA16F_TO_RGBA16F,
	L16F_TO_RGBA16F,
	LA16F_TO_RGBA16F,

	R32F_TO_RGBA16F,
	RG32F_TO_RGBA16F,
	RGB32F_TO_RGBA16F,
	RGBA32F_TO_RGBA16F,
	BGR32F_TO_RGBA16F,
	BGRA32F_TO_RGBA16F,
	L32F_TO_RGBA16F,
	LA32F_TO_RGBA16F
};

// scanline copying function: this function expands data to RGBA, either U8 or FP16.
static void copy_scanline(
	void* dst,
	const void* src,
	int pixels,
	int method
) {

#define id(x) (x)
#define u16_sf16(x) float_to_sf16(x * (1.0f/65535.0f), SF_NEARESTEVEN)
#define f32_sf16(x) sf32_to_sf16(x, SF_NEARESTEVEN)

#define COPY_R(dsttype, srctype, convfunc, oneval) \
	do { \
		srctype *s = (srctype *)src; \
		dsttype *d = (dsttype *)dst; \
		for (i=0;i<pixels;i++)\
		{\
			d[4*i] = convfunc(s[i]); \
			d[4*i+1] = 0; \
			d[4*i+2] = 0; \
			d[4*i+3] = oneval; \
		} \
	} while (0); \
	break;

#define COPY_RG(dsttype, srctype, convfunc, oneval) \
	do { \
		srctype *s = (srctype *)src; \
		dsttype *d = (dsttype *)dst; \
		for (i=0;i<pixels;i++)\
		{\
			d[4*i] = convfunc(s[2*i]); \
			d[4*i+1] = convfunc(s[2*i+1]); \
			d[4*i+2] = 0; \
			d[4*i+3] = oneval; \
		} \
	} while (0); \
	break;

#define COPY_RGB(dsttype, srctype, convfunc, oneval) \
	do { \
		srctype *s = (srctype *)src; \
		dsttype *d = (dsttype *)dst; \
		for (i=0;i<pixels;i++)\
		{\
			d[4*i] = convfunc(s[3*i]); \
			d[4*i+1] = convfunc(s[3*i+1]); \
			d[4*i+2] = convfunc(s[3*i+2]); \
			d[4*i+3] = oneval; \
		} \
	} while (0); \
	break;

#define COPY_BGR(dsttype, srctype, convfunc, oneval) \
	do { \
		srctype *s = (srctype *)src; \
		dsttype *d = (dsttype *)dst; \
		for (i=0;i<pixels;i++)\
		{\
			d[4*i] = convfunc(s[3*i+2]); \
			d[4*i+1] = convfunc(s[3*i+1]); \
			d[4*i+2] = convfunc(s[3*i]); \
			d[4*i+3] = oneval; \
		} \
	} while (0); \
	break;

#define COPY_RGBX(dsttype, srctype, convfunc, oneval) \
	do { \
		srctype *s = (srctype *)src; \
		dsttype *d = (dsttype *)dst; \
		for (i=0;i<pixels;i++)\
		{\
			d[4*i] = convfunc(s[4*i]); \
			d[4*i+1] = convfunc(s[4*i+1]); \
			d[4*i+2] = convfunc(s[4*i+2]); \
			d[4*i+3] = oneval; \
		} \
	} while (0); \
	break;

#define COPY_BGRX(dsttype, srctype, convfunc, oneval) \
	do { \
		srctype *s = (srctype *)src; \
		dsttype *d = (dsttype *)dst; \
		for (i=0;i<pixels;i++)\
		{\
			d[4*i] = convfunc(s[4*i+2]); \
			d[4*i+1] = convfunc(s[4*i+1]); \
			d[4*i+2] = convfunc(s[4*i]); \
			d[4*i+3] = oneval; \
		} \
	} while (0); \
	break;

#define COPY_RGBA(dsttype, srctype, convfunc, oneval) \
	do { \
		srctype *s = (srctype *)src; \
		dsttype *d = (dsttype *)dst; \
		for (i=0;i<pixels;i++)\
		{\
			d[4*i] = convfunc(s[4*i]); \
			d[4*i+1] = convfunc(s[4*i+1]); \
			d[4*i+2] = convfunc(s[4*i+2]); \
			d[4*i+3] = convfunc(s[4*i+3]); \
		} \
	} while (0); \
	break;

#define COPY_BGRA(dsttype, srctype, convfunc, oneval) \
	do { \
		srctype *s = (srctype *)src; \
		dsttype *d = (dsttype *)dst; \
		for (i=0;i<pixels;i++)\
		{\
			d[4*i] = convfunc(s[4*i+2]); \
			d[4*i+1] = convfunc(s[4*i+1]); \
			d[4*i+2] = convfunc(s[4*i]); \
			d[4*i+3] = convfunc(s[4*i+3]); \
		} \
	} while (0); \
	break;

#define COPY_L(dsttype, srctype, convfunc, oneval) \
	do { \
		srctype *s = (srctype *)src; \
		dsttype *d = (dsttype *)dst; \
		for (i=0;i<pixels;i++)\
		{\
			d[4*i] = convfunc(s[i]); \
			d[4*i+1] = convfunc(s[i]); \
			d[4*i+2] = convfunc(s[i]); \
			d[4*i+3] = oneval; \
		} \
	} while (0); \
	break;

#define COPY_LA(dsttype, srctype, convfunc, oneval) \
	do { \
		srctype *s = (srctype *)src; \
		dsttype *d = (dsttype *)dst; \
		for (i=0;i<pixels;i++)\
		{\
			d[4*i] = convfunc(s[2*i]); \
			d[4*i+1] = convfunc(s[2*i]); \
			d[4*i+2] = convfunc(s[2*i]); \
			d[4*i+3] = convfunc(s[2*i+1]); \
		} \
	} while (0); \
	break;

	int i;
	switch (method)
	{
	case R8_TO_RGBA8:
		COPY_R(uint8_t, uint8_t, id, 0xFF);
	case RG8_TO_RGBA8:
		COPY_RG(uint8_t, uint8_t, id, 0xFF);
	case RGB8_TO_RGBA8:
		COPY_RGB(uint8_t, uint8_t, id, 0xFF);
	case RGBA8_TO_RGBA8:
		COPY_RGBA(uint8_t, uint8_t, id, 0xFF);
	case BGR8_TO_RGBA8:
		COPY_BGR(uint8_t, uint8_t, id, 0xFF);
	case BGRA8_TO_RGBA8:
		COPY_BGRA(uint8_t, uint8_t, id, 0xFF);
	case RGBX8_TO_RGBA8:
		COPY_RGBX(uint8_t, uint8_t, id, 0xFF);
	case BGRX8_TO_RGBA8:
		COPY_BGRX(uint8_t, uint8_t, id, 0xFF);
	case L8_TO_RGBA8:
		COPY_L(uint8_t, uint8_t, id, 0xFF);
	case LA8_TO_RGBA8:
		COPY_LA(uint8_t, uint8_t, id, 0xFF);

	case R16F_TO_RGBA16F:
		COPY_R(uint16_t, uint16_t, id, 0x3C00);
	case RG16F_TO_RGBA16F:
		COPY_RG(uint16_t, uint16_t, id, 0x3C00);
	case RGB16F_TO_RGBA16F:
		COPY_RGB(uint16_t, uint16_t, id, 0x3C00);
	case RGBA16F_TO_RGBA16F:
		COPY_RGBA(uint16_t, uint16_t, id, 0x3C00);
	case BGR16F_TO_RGBA16F:
		COPY_BGR(uint16_t, uint16_t, id, 0x3C00);
	case BGRA16F_TO_RGBA16F:
		COPY_BGRA(uint16_t, uint16_t, id, 0x3C00);
	case L16F_TO_RGBA16F:
		COPY_L(uint16_t, uint16_t, id, 0x3C00);
	case LA16F_TO_RGBA16F:
		COPY_LA(uint16_t, uint16_t, id, 0x3C00);

	case R16_TO_RGBA16F:
		COPY_R(uint16_t, uint16_t, u16_sf16, 0x3C00);
	case RG16_TO_RGBA16F:
		COPY_RG(uint16_t, uint16_t, u16_sf16, 0x3C00);
	case RGB16_TO_RGBA16F:
		COPY_RGB(uint16_t, uint16_t, u16_sf16, 0x3C00);
	case RGBA16_TO_RGBA16F:
		COPY_RGBA(uint16_t, uint16_t, u16_sf16, 0x3C00);
	case BGR16_TO_RGBA16F:
		COPY_BGR(uint16_t, uint16_t, u16_sf16, 0x3C00);
	case BGRA16_TO_RGBA16F:
		COPY_BGRA(uint16_t, uint16_t, u16_sf16, 0x3C00);
	case L16_TO_RGBA16F:
		COPY_L(uint16_t, uint16_t, u16_sf16, 0x3C00);
	case LA16_TO_RGBA16F:
		COPY_LA(uint16_t, uint16_t, u16_sf16, 0x3C00);

	case R32F_TO_RGBA16F:
		COPY_R(uint16_t, uint32_t, f32_sf16, 0x3C00);
	case RG32F_TO_RGBA16F:
		COPY_RG(uint16_t, uint32_t, f32_sf16, 0x3C00);
	case RGB32F_TO_RGBA16F:
		COPY_RGB(uint16_t, uint32_t, f32_sf16, 0x3C00);
	case RGBA32F_TO_RGBA16F:
		COPY_RGBA(uint16_t, uint32_t, f32_sf16, 0x3C00);
	case BGR32F_TO_RGBA16F:
		COPY_BGR(uint16_t, uint32_t, f32_sf16, 0x3C00);
	case BGRA32F_TO_RGBA16F:
		COPY_BGRA(uint16_t, uint32_t, f32_sf16, 0x3C00);
	case L32F_TO_RGBA16F:
		COPY_L(uint16_t, uint32_t, f32_sf16, 0x3C00);
	case LA32F_TO_RGBA16F:
		COPY_LA(uint16_t, uint32_t, f32_sf16, 0x3C00);
	};
}

// perform endianness switch on raw data
static void switch_endianness2(
	void* dataptr,
	int bytes
) {
	int i;
	uint8_t *data = (uint8_t *) dataptr;

	for (i = 0; i < bytes / 2; i++)
	{
		uint8_t d0 = data[0];
		uint8_t d1 = data[1];
		data[0] = d1;
		data[1] = d0;
		data += 2;
	}
}

static void switch_endianness4(
	void* dataptr,
	int bytes
) {
	int i;
	uint8_t *data = (uint8_t *) dataptr;

	for (i = 0; i < bytes / 4; i++)
	{
		uint8_t d0 = data[0];
		uint8_t d1 = data[1];
		uint8_t d2 = data[2];
		uint8_t d3 = data[3];
		data[0] = d3;
		data[1] = d2;
		data[2] = d1;
		data[3] = d0;
		data += 4;
	}
}

static uint32_t u32_byterev(uint32_t v)
{
	return (v >> 24) | ((v >> 8) & 0xFF00) | ((v << 8) & 0xFF0000) | (v << 24);
}

/*
	Notes about KTX:

	After the header and the key/value data area, the actual image data follows.
	Each image starts with a 4-byte "imageSize" value indicating the number of bytes of image data follow.
	(For cube-maps, this value appears only after first image; the remaining 5 images are all of equal size.)
	If the size of an image is not a multiple of 4, then it is padded to the next multiple of 4.
	Note that this padding is NOT included in the "imageSize" field.
	In a cubemap, the padding appears after each face note that in a 2D/3D texture, padding does
	NOT appear between the lines/planes of the texture!

	In a KTX file, there may be multiple images; they are organized as follows:

	For each mipmap_level in numberOfMipmapLevels
		UInt32 imageSize;
		For each array_element in numberOfArrayElements
		* for each face in numberOfFaces
			* for each z_slice in pixelDepth
				* for each row or row_of_blocks in pixelHeight
					* for each pixel or block_of_pixels in pixelWidth
						Byte data[format-specific-number-of-bytes]
					* end
				* end
			*end
			Byte cubePadding[0-3]
		*end
		Byte mipPadding[3 - ((imageSize+ 3) % 4)]
	*end

	In the ASTC codec, we will, for the time being only harvest the first image,
	and we will support only a limited set of formats:

	gl_type: UNSIGNED_BYTE UNSIGNED_SHORT HALF_FLOAT FLOAT UNSIGNED_INT_8_8_8_8 UNSIGNED_INT_8_8_8_8_REV
	gl_format: RED, RG. RGB, RGBA BGR, BGRA
	gl_internal_format: used for upload to OpenGL; we can ignore it on uncompressed-load, but
		need to provide a reasonable value on store: RGB8 RGBA8 RGB16F RGBA16F
	gl_base_internal_format: same as gl_format unless texture is compressed (well, BGR is turned into RGB)
		RED, RG, RGB, RGBA
 */

// enums copied from GL/GL.h
#define GL_RED    0x1903
#define GL_RG     0x8227
#define GL_RGB    0x1907
#define GL_RGBA   0x1908
#define GL_BGR    0x80E0
#define GL_BGRA   0x80E1
#define GL_LUMINANCE        0x1909
#define GL_LUMINANCE_ALPHA  0x190A

#define GL_UNSIGNED_BYTE   0x1401
#define GL_UNSIGNED_SHORT  0x1403
#define GL_HALF_FLOAT      0x140B
#define GL_FLOAT           0x1406

struct ktx_header
{
	uint8_t magic[12];
	uint32_t endianness;		// should be 0x04030201; if it is instead 0x01020304, then the endianness of everything must be switched.
	uint32_t gl_type;			// 0 for compressed textures, otherwise value from table 3.2 (page 162) of OpenGL 4.0 spec
	uint32_t gl_type_size;		// size of data elements to do endianness swap on (1=endian-neutral data)
	uint32_t gl_format;			// 0 for compressed textures, otherwise value from table 3.3 (page 163) of OpenGLl spec
	uint32_t gl_internal_format;	// sized-internal-format, corresponding to table 3.12 to 3.14 (pages 182-185) of OpenGL spec
	uint32_t gl_base_internal_format;	// unsized-internal-format: corresponding to table 3.11 (page 179) of OpenGL spec
	uint32_t pixel_width;		// texture dimensions; not rounded up to block size for compressed.
	uint32_t pixel_height;		// must be 0 for 1D textures.
	uint32_t pixel_depth;		// must be 0 for 1D, 2D and cubemap textures.
	uint32_t number_of_array_elements;	// 0 if not a texture array
	uint32_t number_of_faces;	// 6 for cubemaps, 1 for non-cubemaps
	uint32_t number_of_mipmap_levels;	// 0 or 1 for non-mipmapped textures; 0 indicates that auto-mipmap-gen should be done at load time.
	uint32_t bytes_of_key_value_data;	// size in bytes of the key-and-value area immediately following the header.
};

// magic 12-byte sequence that must appear at the beginning of every KTX file.
uint8_t ktx_magic[12] = {
	0xAB, 0x4B, 0x54, 0x58, 0x20, 0x31, 0x31, 0xBB, 0x0D, 0x0A, 0x1A, 0x0A
};

static void ktx_header_switch_endianness(ktx_header * kt)
{
	#define REV(x) kt->x = u32_byterev(kt->x)
	REV(endianness);
	REV(gl_type);
	REV(gl_type_size);
	REV(gl_format);
	REV(gl_internal_format);
	REV(gl_base_internal_format);
	REV(pixel_width);
	REV(pixel_height);
	REV(pixel_depth);
	REV(number_of_array_elements);
	REV(number_of_faces);
	REV(number_of_mipmap_levels);
	REV(bytes_of_key_value_data);
	#undef REV
}

static astcenc_image* load_ktx_uncompressed_image(
	const char* filename,
	unsigned int dim_pad,
	bool y_flip,
	bool& is_hdr,
	unsigned int& num_components
) {
	FILE *f = fopen(filename, "rb");
	if (!f)
	{
		printf("Failed to open file %s\n", filename);
		return nullptr;
	}

	ktx_header hdr;
	size_t header_bytes_read = fread(&hdr, 1, sizeof(hdr), f);

	if (header_bytes_read != sizeof(hdr))
	{
		printf("Failed to read header of KTX file %s\n", filename);
		fclose(f);
		return nullptr;
	}

	if (memcmp(hdr.magic, ktx_magic, 12) != 0 || (hdr.endianness != 0x04030201 && hdr.endianness != 0x01020304))
	{
		printf("File %s does not have a valid KTX header\n", filename);
		fclose(f);
		return nullptr;
	}

	int switch_endianness = 0;
	if (hdr.endianness == 0x01020304)
	{
		ktx_header_switch_endianness(&hdr);
		switch_endianness = 1;
	}

	if (hdr.gl_type == 0 || hdr.gl_format == 0)
	{
		printf("File %s appears to be compressed, not supported as input\n", filename);
		fclose(f);
		return nullptr;
	}

	// the formats we support are:

	// Cartesian product of gl_type=(UNSIGNED_BYTE, UNSIGNED_SHORT, HALF_FLOAT, FLOAT) x gl_format=(RED, RG, RGB, RGBA, BGR, BGRA)

	int components;
	switch (hdr.gl_format)
	{
	case GL_RED:
		components = 1;
		break;
	case GL_RG:
		components = 2;
		break;
	case GL_RGB:
		components = 3;
		break;
	case GL_RGBA:
		components = 4;
		break;
	case GL_BGR:
		components = 3;
		break;
	case GL_BGRA:
		components = 4;
		break;
	case GL_LUMINANCE:
		components = 1;
		break;
	case GL_LUMINANCE_ALPHA:
		components = 2;
		break;
	default:
		printf("KTX file %s has unsupported GL type\n", filename);
		fclose(f);
		return nullptr;
	};

	// Although these are set up later, we include a default initializer to remove warnings
	int bytes_per_component = 1;	// bytes per component in the KTX file.
	int bitness = 8;			// internal precision we will use in the codec.
	scanline_copy_method cm = R8_TO_RGBA8;

	switch (hdr.gl_type)
	{
	case GL_UNSIGNED_BYTE:
		{
			bitness = 8;
			bytes_per_component = 1;
			switch (hdr.gl_format)
			{
			case GL_RED:
				cm = R8_TO_RGBA8;
				break;
			case GL_RG:
				cm = RG8_TO_RGBA8;
				break;
			case GL_RGB:
				cm = RGB8_TO_RGBA8;
				break;
			case GL_RGBA:
				cm = RGBA8_TO_RGBA8;
				break;
			case GL_BGR:
				cm = BGR8_TO_RGBA8;
				break;
			case GL_BGRA:
				cm = BGRA8_TO_RGBA8;
				break;
			case GL_LUMINANCE:
				cm = L8_TO_RGBA8;
				break;
			case GL_LUMINANCE_ALPHA:
				cm = LA8_TO_RGBA8;
				break;
			}
			break;
		}
	case GL_UNSIGNED_SHORT:
		{
			bitness = 16;
			bytes_per_component = 2;
			switch (hdr.gl_format)
			{
			case GL_RED:
				cm = R16_TO_RGBA16F;
				break;
			case GL_RG:
				cm = RG16_TO_RGBA16F;
				break;
			case GL_RGB:
				cm = RGB16_TO_RGBA16F;
				break;
			case GL_RGBA:
				cm = RGBA16_TO_RGBA16F;
				break;
			case GL_BGR:
				cm = BGR16_TO_RGBA16F;
				break;
			case GL_BGRA:
				cm = BGRA16_TO_RGBA16F;
				break;
			case GL_LUMINANCE:
				cm = L16_TO_RGBA16F;
				break;
			case GL_LUMINANCE_ALPHA:
				cm = LA16_TO_RGBA16F;
				break;
			}
			break;
		}
	case GL_HALF_FLOAT:
		{
			bitness = 16;
			bytes_per_component = 2;
			switch (hdr.gl_format)
			{
			case GL_RED:
				cm = R16F_TO_RGBA16F;
				break;
			case GL_RG:
				cm = RG16F_TO_RGBA16F;
				break;
			case GL_RGB:
				cm = RGB16F_TO_RGBA16F;
				break;
			case GL_RGBA:
				cm = RGBA16F_TO_RGBA16F;
				break;
			case GL_BGR:
				cm = BGR16F_TO_RGBA16F;
				break;
			case GL_BGRA:
				cm = BGRA16F_TO_RGBA16F;
				break;
			case GL_LUMINANCE:
				cm = L16F_TO_RGBA16F;
				break;
			case GL_LUMINANCE_ALPHA:
				cm = LA16F_TO_RGBA16F;
				break;
			}
			break;
		}
	case GL_FLOAT:
		{
			bitness = 16;
			bytes_per_component = 4;
			switch (hdr.gl_format)
			{
			case GL_RED:
				cm = R32F_TO_RGBA16F;
				break;
			case GL_RG:
				cm = RG32F_TO_RGBA16F;
				break;
			case GL_RGB:
				cm = RGB32F_TO_RGBA16F;
				break;
			case GL_RGBA:
				cm = RGBA32F_TO_RGBA16F;
				break;
			case GL_BGR:
				cm = BGR32F_TO_RGBA16F;
				break;
			case GL_BGRA:
				cm = BGRA32F_TO_RGBA16F;
				break;
			case GL_LUMINANCE:
				cm = L32F_TO_RGBA16F;
				break;
			case GL_LUMINANCE_ALPHA:
				cm = LA32F_TO_RGBA16F;
				break;
			}
			break;
		}
	default:
		printf("KTX file %s has unsupported GL format\n", filename);
		fclose(f);
		return nullptr;
	}

	if (hdr.number_of_mipmap_levels > 1)
		printf("WARNING: KTX file %s has %d mipmap levels; only the first one will be encoded.\n", filename, hdr.number_of_mipmap_levels);

	if (hdr.number_of_array_elements > 1)
		printf("WARNING: KTX file %s contains a texture array with %d layers; only the first one will be encoded.\n", filename, hdr.number_of_array_elements);

	if (hdr.number_of_faces > 1)
		printf("WARNING: KTX file %s contains a cubemap with 6 faces; only the first one will be encoded.\n", filename);


	unsigned int dim_x = hdr.pixel_width;
	unsigned int dim_y = MAX(hdr.pixel_height, 1);
	unsigned int dim_z = MAX(hdr.pixel_depth, 1);

	// ignore the key/value data
	fseek(f, hdr.bytes_of_key_value_data, SEEK_CUR);

	uint32_t specified_bytes_of_surface = 0;
	size_t sb_read = fread(&specified_bytes_of_surface, 1, 4, f);
	if (sb_read != 4)
	{
		printf("Failed to read header of KTX file %s\n", filename);
		fclose(f);
		return nullptr;
	}

	if (switch_endianness)
		specified_bytes_of_surface = u32_byterev(specified_bytes_of_surface);

	// read the surface
	uint32_t xstride = bytes_per_component * components * dim_x;
	uint32_t ystride = xstride * dim_y;
	uint32_t computed_bytes_of_surface = dim_z * ystride;
	if (computed_bytes_of_surface != specified_bytes_of_surface)
	{
		fclose(f);
		printf("%s: KTX file inconsistency: computed surface size is %d bytes, but specified size is %d bytes\n", filename, computed_bytes_of_surface, specified_bytes_of_surface);
		return nullptr;
	}

	uint8_t *buf = new uint8_t[specified_bytes_of_surface];
	size_t bytes_read = fread(buf, 1, specified_bytes_of_surface, f);
	fclose(f);
	if (bytes_read != specified_bytes_of_surface)
	{
		delete[] buf;
		printf("Failed to read file %s\n", filename);
		return nullptr;
	}

	// perform an endianness swap on the surface if needed.
	if (switch_endianness)
	{
		if (hdr.gl_type_size == 2)
			switch_endianness2(buf, specified_bytes_of_surface);
		if (hdr.gl_type_size == 4)
			switch_endianness4(buf, specified_bytes_of_surface);
	}

	// then transfer data from the surface to our own image-data-structure.
	astcenc_image *astc_img = alloc_image(bitness, dim_x, dim_y, dim_z, dim_pad);

	for (unsigned int z = 0; z < dim_z; z++)
	{
		unsigned int zdst = (dim_z == 1) ? z : z + dim_pad;

		for (unsigned int y = 0; y < dim_y; y++)
		{
			unsigned int ymod = y_flip ? dim_y - y - 1 : y;
			unsigned int ydst = ymod + dim_pad;
			void *dst;

			if (bitness == 16)
				dst = (void *)(astc_img->data16[zdst][ydst] + 4 * dim_pad);
			else
				dst = (void *)(astc_img->data8[zdst][ydst] + 4 * dim_pad);

			uint8_t *src = buf + (z * ystride) + (y * xstride);
			copy_scanline(dst, src, dim_x, cm);
		}
	}

	delete[] buf;
	fill_image_padding_area(astc_img);
	is_hdr = bitness == 16;
	num_components = components;
	return astc_img;
}

static int store_ktx_uncompressed_image(
	const astcenc_image* img,
	const char* ktx_filename,
	int y_flip
) {
	unsigned int dim_x = img->dim_x;
	unsigned int dim_y = img->dim_y;
	unsigned int dim_z = img->dim_z;

	int bitness = img->data16 == nullptr ? 8 : 16;
	int image_channels = determine_image_channels(img);

	ktx_header hdr;

	int gl_format_of_channels[4] = { GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA };

	memcpy(hdr.magic, ktx_magic, 12);
	hdr.endianness = 0x04030201;
	hdr.gl_type = (bitness == 16) ? GL_HALF_FLOAT : GL_UNSIGNED_BYTE;
	hdr.gl_type_size = bitness / 8;
	hdr.gl_format = gl_format_of_channels[image_channels - 1];
	hdr.gl_internal_format = gl_format_of_channels[image_channels - 1];
	hdr.gl_base_internal_format = gl_format_of_channels[image_channels - 1];
	hdr.pixel_width = dim_x;
	hdr.pixel_height = dim_y;
	hdr.pixel_depth = (dim_z == 1) ? 0 : dim_z;
	hdr.number_of_array_elements = 0;
	hdr.number_of_faces = 1;
	hdr.number_of_mipmap_levels = 1;
	hdr.bytes_of_key_value_data = 0;

	// collect image data to write
	uint8_t ***row_pointers8 = nullptr;
	uint16_t ***row_pointers16 = nullptr;
	if (bitness == 8)
	{
		row_pointers8 = new uint8_t **[dim_z];
		row_pointers8[0] = new uint8_t *[dim_y * dim_z];
		row_pointers8[0][0] = new uint8_t[dim_x * dim_y * dim_z * image_channels + 3];

		for (unsigned int z = 1; z < dim_z; z++)
		{
			row_pointers8[z] = row_pointers8[0] + dim_y * z;
			row_pointers8[z][0] = row_pointers8[0][0] + dim_y * dim_x * image_channels * z;
		}

		for (unsigned int z = 0; z < dim_z; z++)
		{
			for (unsigned int y = 1; y < dim_y; y++)
			{
				row_pointers8[z][y] = row_pointers8[z][0] + dim_x * image_channels * y;
			}
		}

		for (unsigned int z = 0; z < dim_z; z++)
		{
			for (unsigned int y = 0; y < dim_y; y++)
			{
				int ym = y_flip ? dim_y - y - 1 : y;
				switch (image_channels)
				{
				case 1:		// single-component, treated as Luminance
					for (unsigned int x = 0; x < dim_x; x++)
					{
						row_pointers8[z][y][x] = img->data8[z][ym][4 * x];
					}
					break;
				case 2:		// two-component, treated as Luminance-Alpha
					for (unsigned int x = 0; x < dim_x; x++)
					{
						row_pointers8[z][y][2 * x]     = img->data8[z][ym][4 * x];
						row_pointers8[z][y][2 * x + 1] = img->data8[z][ym][4 * x + 3];
					}
					break;
				case 3:		// three-component, treated as RGB
					for (unsigned int x = 0; x < dim_x; x++)
					{
						row_pointers8[z][y][3 * x]     = img->data8[z][ym][4 * x];
						row_pointers8[z][y][3 * x + 1] = img->data8[z][ym][4 * x + 1];
						row_pointers8[z][y][3 * x + 2] = img->data8[z][ym][4 * x + 2];
					}
					break;
				case 4:		// four-component, treated as RGBA
					for (unsigned int x = 0; x < dim_x; x++)
					{
						row_pointers8[z][y][4 * x]     = img->data8[z][ym][4 * x];
						row_pointers8[z][y][4 * x + 1] = img->data8[z][ym][4 * x + 1];
						row_pointers8[z][y][4 * x + 2] = img->data8[z][ym][4 * x + 2];
						row_pointers8[z][y][4 * x + 3] = img->data8[z][ym][4 * x + 3];
					}
					break;
				}
			}
		}
	}
	else						// if bitness == 16
	{
		row_pointers16 = new uint16_t **[dim_z];
		row_pointers16[0] = new uint16_t *[dim_y * dim_z];
		row_pointers16[0][0] = new uint16_t[dim_x * dim_y * dim_z * image_channels + 1];

		for (unsigned int z = 1; z < dim_z; z++)
		{
			row_pointers16[z] = row_pointers16[0] + dim_y * z;
			row_pointers16[z][0] = row_pointers16[0][0] + dim_y * dim_x * image_channels * z;
		}

		for (unsigned int z = 0; z < dim_z; z++)
		{
			for (unsigned int y = 1; y < dim_y; y++)
			{
				row_pointers16[z][y] = row_pointers16[z][0] + dim_x * image_channels * y;
			}
		}

		for (unsigned int z = 0; z < dim_z; z++)
		{
			for (unsigned int y = 0; y < dim_y; y++)
			{
				int ym = y_flip ? dim_y - y - 1 : y;
				switch (image_channels)
				{
				case 1:		// single-component, treated as Luminance
					for (unsigned int x = 0; x < dim_x; x++)
					{
						row_pointers16[z][y][x] = img->data16[z][ym][4 * x];
					}
					break;
				case 2:		// two-component, treated as Luminance-Alpha
					for (unsigned int x = 0; x < dim_x; x++)
					{
						row_pointers16[z][y][2 * x]     = img->data16[z][ym][4 * x];
						row_pointers16[z][y][2 * x + 1] = img->data16[z][ym][4 * x + 3];
					}
					break;
				case 3:		// three-component, treated as RGB
					for (unsigned int x = 0; x < dim_x; x++)
					{
						row_pointers16[z][y][3 * x]     = img->data16[z][ym][4 * x];
						row_pointers16[z][y][3 * x + 1] = img->data16[z][ym][4 * x + 1];
						row_pointers16[z][y][3 * x + 2] = img->data16[z][ym][4 * x + 2];
					}
					break;
				case 4:		// four-component, treated as RGBA
					for (unsigned int x = 0; x < dim_x; x++)
					{
						row_pointers16[z][y][4 * x]     = img->data16[z][ym][4 * x];
						row_pointers16[z][y][4 * x + 1] = img->data16[z][ym][4 * x + 1];
						row_pointers16[z][y][4 * x + 2] = img->data16[z][ym][4 * x + 2];
						row_pointers16[z][y][4 * x + 3] = img->data16[z][ym][4 * x + 3];
					}
					break;
				}
			}
		}
	}

	int retval = image_channels + (bitness == 16 ? 0x80 : 0);
	uint32_t image_bytes = dim_x * dim_y * dim_z * image_channels * (bitness / 8);
	uint32_t image_write_bytes = (image_bytes + 3) & ~3;

	FILE *wf = fopen(ktx_filename, "wb");
	if (wf)
	{
		void *dataptr = (bitness == 16) ? (void *)(row_pointers16[0][0]) : (void *)(row_pointers8[0][0]);

		size_t expected_bytes_written = sizeof(ktx_header) + image_write_bytes + 4;
		size_t hdr_bytes_written = fwrite(&hdr, 1, sizeof(ktx_header), wf);
		size_t bytecount_bytes_written = fwrite(&image_bytes, 1, 4, wf);
		size_t data_bytes_written = fwrite(dataptr, 1, image_write_bytes, wf);
		fclose(wf);
		if (hdr_bytes_written + bytecount_bytes_written + data_bytes_written != expected_bytes_written)
			retval = -1;
	}
	else
	{
		retval = -1;
	}

	if (row_pointers8)
	{
		delete[] row_pointers8[0][0];
		delete[] row_pointers8[0];
		delete[] row_pointers8;
	}

	if (row_pointers16)
	{
		delete[] row_pointers16[0][0];
		delete[] row_pointers16[0];
		delete[] row_pointers16;
	}

	return retval;
}

/*
	Loader for DDS files.

	Note that after the header, data are densely packed with no padding;
	in the case of multiple surfaces, they appear one after another in
	the file, again with no padding.

	This code is NOT endian-neutral.
*/
struct dds_pixelformat
{
	uint32_t size;				// structure size, set to 32.
	/*
	   flags bits are a combination of the following: 0x1 : Texture contains alpha data 0x2 : ---- (older files: texture contains alpha data, for Alpha-only texture) 0x4 : The fourcc field is valid,
	   indicating a compressed or DX10 texture format 0x40 : texture contains uncompressed RGB data 0x200 : ---- (YUV in older files) 0x20000 : Texture contains Luminance data (can be combined with
	   0x1 for Lum-Alpha) */
	uint32_t flags;
	uint32_t fourcc;			// "DX10" to indicate a DX10 format, "DXTn" for the DXT formats
	uint32_t rgbbitcount;		// number of bits per texel; up to 32 for non-DX10 formats.
	uint32_t rbitmask;			// bitmap indicating position of red/luminance color component
	uint32_t gbitmask;			// bitmap indicating position of green color component
	uint32_t bbitmask;			// bitmap indicating position of blue color component
	uint32_t abitmask;			// bitmap indicating position of alpha color component
};

struct dds_header
{
	uint32_t size;				// header size; must be exactly 124.
	/*
	   flag field is an OR or the following bits, that indicate fields containing valid data:
		1: caps/caps2/caps3/caps4 (set in all DDS files, ignore on read)
		2: height (set in all DDS files, ignore on read)
		4: width (set in all DDS files, ignore on read)
		8: pitch (for uncompressed texture)
		0x1000: the pixel format field (set in all DDS files, ignore on read)
		0x20000: mipmap count (for mipmapped textures with >1 level)
		0x80000: pitch (for compressed texture)
		0x800000: depth (for 3d textures)
	*/
	uint32_t flags;
	uint32_t height;
	uint32_t width;
	uint32_t pitch_or_linear_size;	// scanline pitch for uncompressed; total size in bytes for compressed
	uint32_t depth;
	uint32_t mipmapcount;
	// unused, set to 0
	uint32_t reserved1[11];
	dds_pixelformat ddspf;
	/*
	   caps field is an OR of the following values:
		8 : should be set for a file that contains more than 1 surface (ignore on read)
		0x400000 : should be set for a mipmapped texture
		0x1000 : should be set if the surface is a texture at all (all DDS files, ignore on read)
	*/
	uint32_t caps;
	/*
	   caps2 field is an OR of the following values:
		0x200 : texture is cubemap
		0x400 : +X face of cubemap is present
		0x800 : -X face of cubemap is present
		0x1000 : +Y face of cubemap is present
		0x2000 : -Y face of cubemap is present
		0x4000 : +Z face of cubemap is present
		0x8000 : -Z face of cubemap is present
		0x200000 : texture is a 3d texture.
	*/
	uint32_t caps2;
	// unused, set to 0
	uint32_t caps3;
	// unused, set to 0
	uint32_t caps4;
	// unused, set to 0
	uint32_t reserved2;
};

struct dds_header_dx10
{
	uint32_t dxgi_format;
	uint32_t resource_dimension;	// 2=1d-texture, 3=2d-texture or cubemap, 4=3d-texture
	uint32_t misc_flag;			// 4 if cubemap, else 0
	uint32_t array_size;		// size of array in case of a texture array; set to 1 for a non-array
	uint32_t reserved;			// set to 0.
};

#define DDS_MAGIC 0x20534444
#define DX10_MAGIC 0x30315844

astcenc_image* load_dds_uncompressed_image(
	const char* filename,
	unsigned int dim_pad,
	bool y_flip,
	bool& is_hdr,
	unsigned int& num_components
) {
	FILE *f = fopen(filename, "rb");
	if (!f)
	{
		printf("Failed to open file %s\n", filename);
		return nullptr;
	}

	uint8_t magic[4];

	dds_header hdr;
	size_t magic_bytes_read = fread(magic, 1, 4, f);
	size_t header_bytes_read = fread(&hdr, 1, sizeof(hdr), f);
	if (magic_bytes_read != 4 || header_bytes_read != sizeof(hdr))
	{
		printf("Failed to read header of DDS file %s\n", filename);
		fclose(f);
		return nullptr;
	}

	uint32_t magicx = magic[0] | (magic[1] << 8) | (magic[2] << 16) | (magic[3] << 24);

	if (magicx != DDS_MAGIC || hdr.size != 124)
	{
		printf("File %s does not have a valid DDS header\n", filename);
		fclose(f);
		return nullptr;
	}

	int use_dx10_header = 0;
	if (hdr.ddspf.flags & 4)
	{
		if (hdr.ddspf.fourcc == DX10_MAGIC)
		{
			use_dx10_header = 1;
		}
		else
		{
			printf("DDS file %s is compressed, not supported\n", filename);
			fclose(f);
			return nullptr;
		}
	}

	dds_header_dx10 dx10_header;
	if (use_dx10_header)
	{
		size_t dx10_header_bytes_read = fread(&dx10_header, 1, sizeof(dx10_header), f);
		if (dx10_header_bytes_read != sizeof(dx10_header))
		{
			printf("Failed to read header of DDS file %s\n", filename);
			fclose(f);
			return nullptr;
		}
	}

	unsigned int dim_x = hdr.width;
	unsigned int dim_y = hdr.height;
	unsigned int dim_z = (hdr.flags & 0x800000) ? hdr.depth : 1;

	int bitness;				// the bitcount that we will use internally in the codec
	int bytes_per_component;	// the bytes per component in the DDS file itself
	int components;
	int copy_method;

	// figure out the format actually used in the DDS file.
	if (use_dx10_header)
	{
		// DX10 header present; use the DXGI format.
		#define DXGI_FORMAT_R32G32B32A32_FLOAT   2
		#define DXGI_FORMAT_R32G32B32_FLOAT      6
		#define DXGI_FORMAT_R16G16B16A16_FLOAT  10
		#define DXGI_FORMAT_R16G16B16A16_UNORM  11
		#define DXGI_FORMAT_R32G32_FLOAT        16
		#define DXGI_FORMAT_R8G8B8A8_UNORM      28
		#define DXGI_FORMAT_R16G16_FLOAT    34
		#define DXGI_FORMAT_R16G16_UNORM    35
		#define DXGI_FORMAT_R32_FLOAT       41
		#define DXGI_FORMAT_R8G8_UNORM      49
		#define DXGI_FORMAT_R16_FLOAT       54
		#define DXGI_FORMAT_R16_UNORM       56
		#define DXGI_FORMAT_R8_UNORM        61
		#define DXGI_FORMAT_B8G8R8A8_UNORM  86
		#define DXGI_FORMAT_B8G8R8X8_UNORM  87

		struct dxgi_params
		{
			int bitness;
			int bytes_per_component;
			int components;
			int copy_method;
			uint32_t dxgi_format_number;
		};

		static const dxgi_params format_params[] = {
			{16, 4, 4, RGBA32F_TO_RGBA16F, DXGI_FORMAT_R32G32B32A32_FLOAT},
			{16, 4, 3, RGB32F_TO_RGBA16F, DXGI_FORMAT_R32G32B32_FLOAT},
			{16, 2, 4, RGBA16F_TO_RGBA16F, DXGI_FORMAT_R16G16B16A16_FLOAT},
			{16, 2, 4, RGBA16_TO_RGBA16F, DXGI_FORMAT_R16G16B16A16_UNORM},
			{16, 4, 2, RG32F_TO_RGBA16F, DXGI_FORMAT_R32G32_FLOAT},
			{8, 1, 4, RGBA8_TO_RGBA8, DXGI_FORMAT_R8G8B8A8_UNORM},
			{16, 2, 2, RG16F_TO_RGBA16F, DXGI_FORMAT_R16G16_FLOAT},
			{16, 2, 2, RG16_TO_RGBA16F, DXGI_FORMAT_R16G16_UNORM},
			{16, 4, 1, R32F_TO_RGBA16F, DXGI_FORMAT_R32_FLOAT},
			{8, 1, 2, RG8_TO_RGBA8, DXGI_FORMAT_R8G8_UNORM},
			{16, 2, 1, R16F_TO_RGBA16F, DXGI_FORMAT_R16_FLOAT},
			{16, 2, 1, R16_TO_RGBA16F, DXGI_FORMAT_R16_UNORM},
			{8, 1, 1, R8_TO_RGBA8, DXGI_FORMAT_R8_UNORM},
			{8, 1, 4, BGRA8_TO_RGBA8, DXGI_FORMAT_B8G8R8A8_UNORM},
			{8, 1, 4, BGRX8_TO_RGBA8, DXGI_FORMAT_B8G8R8X8_UNORM},
		};

		int dxgi_modes_supported = sizeof(format_params) / sizeof(format_params[0]);
		int did_select_format = 0;
		for (int i = 0; i < dxgi_modes_supported; i++)
		{
			if (dx10_header.dxgi_format == format_params[i].dxgi_format_number)
			{
				bitness = format_params[i].bitness;
				bytes_per_component = format_params[i].bytes_per_component;
				components = format_params[i].components;
				copy_method = format_params[i].copy_method;
				did_select_format = 1;
				break;
			}
		}

		if (!did_select_format)
		{
			printf("DDS file %s: DXGI format not supported by codec\n", filename);
			fclose(f);
			return nullptr;
		}
	}
	else
	{
		// No DX10 header present. Then try to match the bitcount and bitmask against
		// a set of prepared patterns.
		uint32_t flags = hdr.ddspf.flags;
		uint32_t bitcount = hdr.ddspf.rgbbitcount;
		uint32_t rmask = hdr.ddspf.rbitmask;
		uint32_t gmask = hdr.ddspf.gbitmask;
		uint32_t bmask = hdr.ddspf.bbitmask;
		uint32_t amask = hdr.ddspf.abitmask;

		// RGBA-unorm8
		if ((flags & 0x41) == 0x41 && bitcount == 32 && rmask == 0xFF && gmask == 0xFF00 && bmask == 0xFF0000 && amask == 0xFF000000)
		{
			bytes_per_component = 1;
			components = 4;
			copy_method = RGBA8_TO_RGBA8;
		}
		// BGRA-unorm8
		else if ((flags & 0x41) == 0x41 && bitcount == 32 && rmask == 0xFF0000 && gmask == 0xFF00 && bmask == 0xFF && amask == 0xFF000000)
		{
			bytes_per_component = 1;
			components = 4;
			copy_method = BGRA8_TO_RGBA8;
		}
		// RGBX-unorm8
		else if ((flags & 0x40) && bitcount == 32 && rmask == 0xFF && gmask == 0xFF00 && bmask == 0xFF0000)
		{
			bytes_per_component = 1;
			components = 4;
			copy_method = RGBX8_TO_RGBA8;
		}
		// BGRX-unorm8
		else if ((flags & 0x40) && bitcount == 32 && rmask == 0xFF0000 && gmask == 0xFF00 && bmask == 0xFF)
		{
			bytes_per_component = 1;
			components = 4;
			copy_method = BGRX8_TO_RGBA8;
		}
		// RGB-unorm8
		else if ((flags & 0x40) && bitcount == 24 && rmask == 0xFF && gmask == 0xFF00 && bmask == 0xFF0000)
		{
			bytes_per_component = 1;
			components = 3;
			copy_method = RGB8_TO_RGBA8;
		}
		// BGR-unorm8
		else if ((flags & 0x40) && bitcount == 24 && rmask == 0xFF0000 && gmask == 0xFF00 && bmask == 0xFF)
		{
			bytes_per_component = 1;
			components = 3;
			copy_method = BGR8_TO_RGBA8;
		}
		// RG-unorm16
		else if ((flags & 0x40) && bitcount == 16 && rmask == 0xFFFF && gmask == 0xFFFF0000)
		{
			bytes_per_component = 2;
			components = 2;
			copy_method = RG16_TO_RGBA16F;
		}
		// A8L8
		else if ((flags & 0x20001) == 0x20001 && bitcount == 16 && rmask == 0xFF && amask == 0xFF00)
		{
			bytes_per_component = 1;
			components = 2;
			copy_method = LA8_TO_RGBA8;
		}
		// L8
		else if ((flags & 0x20000) && bitcount == 8 && rmask == 0xFF)
		{
			bytes_per_component = 1;
			components = 1;
			copy_method = L8_TO_RGBA8;
		}
		// L16
		else if ((flags & 0x20000) && bitcount == 16 && rmask == 0xFFFF)
		{
			bytes_per_component = 2;
			components = 1;
			copy_method = L16_TO_RGBA16F;
		}
		else
		{
			printf("DDS file %s: Non-DXGI format not supported by codec\n", filename);
			fclose(f);
			return nullptr;
		}

		bitness = bytes_per_component * 8;
	}

	// then, load the actual file.
	uint32_t xstride = bytes_per_component * components * dim_x;
	uint32_t ystride = xstride * dim_y;
	uint32_t bytes_of_surface = ystride * dim_z;

	uint8_t *buf = new uint8_t[bytes_of_surface];
	size_t bytes_read = fread(buf, 1, bytes_of_surface, f);
	fclose(f);
	if (bytes_read != bytes_of_surface)
	{
		delete[] buf;
		printf("Failed to read file %s\n", filename);
		return nullptr;
	}

	// then transfer data from the surface to our own image-data-structure.
	astcenc_image *astc_img = alloc_image(bitness, dim_x, dim_y, dim_z, dim_pad);

	for (unsigned int z = 0; z < dim_z; z++)
	{
		unsigned int zdst = dim_z == 1 ? z : z + dim_pad;

		for (unsigned int y = 0; y < dim_y; y++)
		{
			unsigned int ymod = y_flip ? dim_y - y - 1 : y;
			unsigned int ydst = ymod + dim_pad;
			void* dst;

			if (bitness == 16)
			{
				dst = (void*)(astc_img->data16[zdst][ydst] + 4 * dim_pad);
			}
			else
			{
				dst = (void*)(astc_img->data8[zdst][ydst] + 4 * dim_pad);
			}

			uint8_t *src = buf + (z * ystride) + (y * xstride);
			copy_scanline(dst, src, dim_x, copy_method);
		}
	}

	delete[] buf;
	fill_image_padding_area(astc_img);
	is_hdr = bitness == 16;
	num_components = components;
	return astc_img;
}

static int store_dds_uncompressed_image(
	const astcenc_image* img,
	const char* dds_filename,
	int y_flip
) {
	unsigned int dim_x = img->dim_x;
	unsigned int dim_y = img->dim_y;
	unsigned int dim_z = img->dim_z;

	int bitness = img->data16 == nullptr ? 8 : 16;
	int image_channels = (bitness == 16) ? 4 : determine_image_channels(img);


	// DDS-pixel-format structures to use when storing LDR image with 1,2,3 or 4 components.
	static const dds_pixelformat format_of_image_channels[4] =
	{
		{32, 0x20000, 0, 8, 0xFF, 0, 0, 0},	// luminance
		{32, 0x20001, 0, 16, 0xFF, 0, 0, 0xFF00},	// L8A8
		{32, 0x40, 0, 24, 0xFF, 0xFF00, 0xFF0000, 0},	// RGB8
		{32, 0x41, 0, 32, 0xFF, 0xFF00, 0xFF0000, 0xFF000000}	// RGBA8
	};

	// DDS-pixel-format structures to use when storing HDR image.
	static const dds_pixelformat dxt10_diverter =
	{
		32, 4, DX10_MAGIC, 0, 0, 0, 0, 0
	};

	// header handling. We will write:
	// * DDS magic value
	// * DDS header
	// * DDS DX10 header, if the file is floating-point
	// * pixel data.

	// main header data
	dds_header hdr;
	hdr.size = 124;
	hdr.flags = 0x100F | (dim_z > 1 ? 0x800000 : 0);
	hdr.height = dim_y;
	hdr.width = dim_x;
	hdr.pitch_or_linear_size = image_channels * (bitness / 8) * dim_x;
	hdr.depth = dim_z;
	hdr.mipmapcount = 1;
	for (unsigned int i = 0; i < 11; i++)
	{
		hdr.reserved1[i] = 0;
	}
	hdr.caps = 0x1000;
	hdr.caps2 = (dim_z > 1) ? 0x200000 : 0;
	hdr.caps3 = 0;
	hdr.caps4 = 0;

	// pixel-format data
	if (bitness == 8)
	{
		hdr.ddspf = format_of_image_channels[image_channels - 1];
	}
	else
	{
		hdr.ddspf = dxt10_diverter;
	}

	// DX10 data
	dds_header_dx10 dx10;
	dx10.dxgi_format = DXGI_FORMAT_R16G16B16A16_FLOAT;
	dx10.resource_dimension = (dim_z > 1) ? 4 : 3;
	dx10.misc_flag = 0;
	dx10.array_size = 1;
	dx10.reserved = 0;

	// collect image data to write
	uint8_t ***row_pointers8 = nullptr;
	uint16_t ***row_pointers16 = nullptr;

	if (bitness == 8)
	{
		row_pointers8 = new uint8_t **[dim_z];
		row_pointers8[0] = new uint8_t *[dim_y * dim_z];
		row_pointers8[0][0] = new uint8_t[dim_x * dim_y * dim_z * image_channels];

		for (unsigned int z = 1; z < dim_z; z++)
		{
			row_pointers8[z] = row_pointers8[0] + dim_y * z;
			row_pointers8[z][0] = row_pointers8[0][0] + dim_y * dim_z * image_channels * z;
		}

		for (unsigned int z = 0; z < dim_z; z++)
		{
			for (unsigned int y = 1; y < dim_y; y++)
			{
				row_pointers8[z][y] = row_pointers8[z][0] + dim_x * image_channels * y;
			}
		}

		for (unsigned int z = 0; z < dim_z; z++)
		{
			for (unsigned int y = 0; y < dim_y; y++)
			{
				int ym = y_flip ? dim_y - y - 1 : y;
				switch (image_channels)
				{
				case 1:		// single-component, treated as Luminance
					for (unsigned int x = 0; x < dim_x; x++)
					{
						row_pointers8[z][y][x] = img->data8[z][ym][4 * x];
					}
					break;
				case 2:		// two-component, treated as Luminance-Alpha
					for (unsigned int x = 0; x < dim_x; x++)
					{
						row_pointers8[z][y][2 * x]     = img->data8[z][ym][4 * x];
						row_pointers8[z][y][2 * x + 1] = img->data8[z][ym][4 * x + 3];
					}
					break;
				case 3:		// three-component, treated as RGB
					for (unsigned int x = 0; x < dim_x; x++)
					{
						row_pointers8[z][y][3 * x]     = img->data8[z][ym][4 * x];
						row_pointers8[z][y][3 * x + 1] = img->data8[z][ym][4 * x + 1];
						row_pointers8[z][y][3 * x + 2] = img->data8[z][ym][4 * x + 2];
					}
					break;
				case 4:		// four-component, treated as RGBA
					for (unsigned int x = 0; x < dim_x; x++)
					{
						row_pointers8[z][y][4 * x]     = img->data8[z][ym][4 * x];
						row_pointers8[z][y][4 * x + 1] = img->data8[z][ym][4 * x + 1];
						row_pointers8[z][y][4 * x + 2] = img->data8[z][ym][4 * x + 2];
						row_pointers8[z][y][4 * x + 3] = img->data8[z][ym][4 * x + 3];
					}
					break;
				}
			}
		}
	}
	else						// if bitness == 16
	{
		row_pointers16 = new uint16_t **[dim_z];
		row_pointers16[0] = new uint16_t *[dim_y * dim_z];
		row_pointers16[0][0] = new uint16_t[dim_x * dim_y * dim_z * image_channels];

		for (unsigned int z = 1; z < dim_z; z++)
		{
			row_pointers16[z] = row_pointers16[0] + dim_y * z;
			row_pointers16[z][0] = row_pointers16[0][0] + dim_y * dim_x * image_channels * z;
		}

		for (unsigned int z = 0; z < dim_z; z++)
		{
			for (unsigned int y = 1; y < dim_y; y++)
			{
				row_pointers16[z][y] = row_pointers16[z][0] + dim_x * image_channels * y;
			}
		}

		for (unsigned int z = 0; z < dim_z; z++)
		{
			for (unsigned int y = 0; y < dim_y; y++)
			{
				int ym = y_flip ? dim_y - y - 1: y;
				switch (image_channels)
				{
				case 1:		// single-component, treated as Luminance
					for (unsigned int x = 0; x < dim_x; x++)
					{
						row_pointers16[z][y][x] = img->data16[z][ym][4 * x];
					}
					break;
				case 2:		// two-component, treated as Luminance-Alpha
					for (unsigned int x = 0; x < dim_x; x++)
					{
						row_pointers16[z][y][2 * x]     = img->data16[z][ym][4 * x];
						row_pointers16[z][y][2 * x + 1] = img->data16[z][ym][4 * x + 3];
					}
					break;
				case 3:		// three-component, treated as RGB
					for (unsigned int x = 0; x < dim_x; x++)
					{
						row_pointers16[z][y][3 * x]     = img->data16[z][ym][4 * x];
						row_pointers16[z][y][3 * x + 1] = img->data16[z][ym][4 * x + 1];
						row_pointers16[z][y][3 * x + 2] = img->data16[z][ym][4 * x + 2];
					}
					break;
				case 4:		// four-component, treated as RGBA
					for (unsigned int x = 0; x < dim_x; x++)
					{
						row_pointers16[z][y][4 * x]     = img->data16[z][ym][4 * x];
						row_pointers16[z][y][4 * x + 1] = img->data16[z][ym][4 * x + 1];
						row_pointers16[z][y][4 * x + 2] = img->data16[z][ym][4 * x + 2];
						row_pointers16[z][y][4 * x + 3] = img->data16[z][ym][4 * x + 3];
					}
					break;
				}
			}
		}
	}

	int retval = image_channels;
	uint32_t image_bytes = dim_x * dim_y * dim_z * image_channels * (bitness / 8);

	uint32_t dds_magic = DDS_MAGIC;

	FILE *wf = fopen(dds_filename, "wb");
	if (wf)
	{
		void *dataptr = (bitness == 16) ? (void *)(row_pointers16[0][0]) : (void *)(row_pointers8[0][0]);

		size_t expected_bytes_written = 4 + sizeof(dds_header) + (bitness > 8 ? sizeof(dds_header_dx10) : 0) + image_bytes;

		size_t magic_bytes_written = fwrite(&dds_magic, 1, 4, wf);
		size_t hdr_bytes_written = fwrite(&hdr, 1, sizeof(dds_header), wf);

		size_t dx10_bytes_written;
		if (bitness > 8)
			dx10_bytes_written = fwrite(&dx10, 1, sizeof(dx10), wf);
		else
			dx10_bytes_written = 0;

		size_t data_bytes_written = fwrite(dataptr, 1, image_bytes, wf);

		fclose(wf);
		if (magic_bytes_written + hdr_bytes_written + dx10_bytes_written + data_bytes_written != expected_bytes_written)
			retval = -1;
	}
	else
	{
		retval = -1;
	}

	if (row_pointers8)
	{
		delete[] row_pointers8[0][0];
		delete[] row_pointers8[0];
		delete[] row_pointers8;
	}

	if (row_pointers16)
	{
		delete[] row_pointers16[0][0];
		delete[] row_pointers16[0];
		delete[] row_pointers16;
	}

	return retval;
}

/***************************************************************************

Main image load/store functions

We have specialized loaders for DDS and KTX; for other formats,
we use stb_image. This image loader will choose one based on filename.

We have specialized image storer for DDS and KTX; for OpenEXR,
we use tinyexr; for TGA, BMP and PNG, we use stb_image_write.

***************************************************************************/






// descriptors for each image/texture file format that we support loading of - endings and
// loader function. The last entry is a catch-all to use when nothing else matches;
// this will result in an attempt to use stb_image to load the image.
static const struct {
	const char* ending1;
	const char* ending2;
	astcenc_image* (*loader_func)(const char*, unsigned int, bool, bool&, unsigned int&);
} loader_descs[] = {
	// HDR formats
	{".exr",   ".EXR",  load_image_with_tinyexr },
	// Container formats
	{".ktx",   ".KTX",  load_ktx_uncompressed_image },
	{".dds",   ".DDS",  load_dds_uncompressed_image },
	// Generic catch all; this one must be last in the list
	{ nullptr, nullptr, load_image_with_stb }
};

static const int loader_descr_count = sizeof(loader_descs) / sizeof(loader_descs[0]);

// descriptors for each image/texture file format that we support storing to to - endings,
// enforced-bitness, and storer function.
static const struct
{
	const char *ending1;
	const char *ending2;
	const char *file_format_name;
	int enforced_bitness;
	int (*storer_func)(const astcenc_image *output_image, const char *output_filename, int y_flip);
} storer_descs[] = {
	// LDR formats
	{".bmp", ".BMP", "BMP",             8, store_bmp_image_with_stb},
	{".png", ".PNG", "PNG",             8, store_png_image_with_stb},
	{".tga", ".TGA", "Targa",           8, store_tga_image_with_stb},
	// HDR formats
	{".exr", ".EXR", "OpenEXR",        16, store_exr_image_with_tinyexr},
	// Container formats
	{".dds", ".DDS", "DirectDraw DDS", -1, store_dds_uncompressed_image},
	{".ktx", ".KTX", "Khronos KTX",    -1, store_ktx_uncompressed_image}
};

static const int storer_descr_count = sizeof(storer_descs) / sizeof(storer_descs[0]);

// check from filename ending what the enforced bitness of the format-to-store is.
// May return:
//  8:  enforced 8-bit UNOR8
//  16: enforced 16-bit FP16
//  -1: no format enforced
// If the format has an unrecognized ending, an error message is produced.
// Lack of an ending is likely to result from a write to /dev/null
// or some other non-file thing; for this, we use the KTX format.

int get_output_filename_enforced_bitness(
	const char*output_filename
) {
	const char *eptr = strrchr(output_filename, '.');
	if (!eptr)
	{
		return -1;
	}

	for (int i = 0; i < storer_descr_count; i++)
	{
		if (strcmp(eptr, storer_descs[i].ending1) == 0
		 || strcmp(eptr, storer_descs[i].ending2) == 0)
		{
			return storer_descs[i].enforced_bitness;
		}
	}

	printf("ERROR: Unknown file extension for output file: %s\n", eptr);
	exit(1);
}

astcenc_image* astc_codec_load_image(
	const char* filename,
	unsigned int dim_pad,
	bool y_flip,
	bool& is_hdr,
	unsigned int& num_components
) {
	// Get the file extension
	const char* eptr = strrchr(filename, '.');
	if (!eptr)
	{
		eptr = filename;
	}

	// Scan through descriptors until a matching loader is found
	for (unsigned int i = 0; i < loader_descr_count; i++)
	{
		if (loader_descs[i].ending1 == nullptr
			|| strcmp(eptr, loader_descs[i].ending1) == 0
			|| strcmp(eptr, loader_descs[i].ending2) == 0)
		{
			return loader_descs[i].loader_func(filename, dim_pad, y_flip, is_hdr, num_components);
		}
	}

	// Should never reach here - stb_image provides a generic handler
	return nullptr;
}

int astc_codec_store_image(
	const astcenc_image* output_image,
	const char* output_filename,
	const char** file_format_name,
	int y_flip
) {
	const char* eptr = strrchr(output_filename, '.');
	if (!eptr)
	{
		eptr = ".ktx"; // use KTX file format if we don't have an ending.
	}

	for (int i=0; i < storer_descr_count; i++)
	{
		if (strcmp(eptr, storer_descs[i].ending1) == 0
		 || strcmp(eptr, storer_descs[i].ending2) == 0)
		{
			*file_format_name = storer_descs[i].file_format_name;
			return storer_descs[i].storer_func(output_image, output_filename, y_flip);
		}
	}

	// Should never reach here - get_output_filename_enforced_bitness should
	// have acted as a preflight check
	return -1;
}
