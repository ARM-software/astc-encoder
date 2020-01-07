// ----------------------------------------------------------------------------
//  This confidential and proprietary software may be used only as authorised
//  by a licensing agreement from Arm Limited.
//      (C) COPYRIGHT 2011-2020 Arm Limited, ALL RIGHTS RESERVED
//  The entire notice above must be reproduced on all authorised copies and
//  copies may only be made to the extent permitted by a licensing agreement
//  from Arm Limited.
// ----------------------------------------------------------------------------

/**
 * @brief Functions for loading/storing ASTC compressed images.
 */


#include "astc_codec_internals.h"

#include <cstdio>
#include <cstring>

void destroy_image(astc_codec_image * img)
{
	if (img == NULL)
		return;

	if (img->data8)
	{
		delete[]img->data8[0][0];
		delete[]img->data8[0];
		delete[]img->data8;
	}

	if (img->data16)
	{
		delete[]img->data16[0][0];
		delete[]img->data16[0];
		delete[]img->data16;
	}

	if(img->input_averages)
		delete[] img->input_averages;
	if(img->input_variances )
		delete[] img->input_variances;
	if(img->input_alpha_averages )
		delete[] img->input_alpha_averages;

	delete img;
}

astc_codec_image *allocate_image(
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

	img->input_averages = NULL;
	img->input_variances = NULL;
	img->input_alpha_averages = NULL;

	int exsize = xsize + 2 * padding;
	int eysize = ysize + 2 * padding;
	int ezsize = (zsize == 1) ? 1 : zsize + 2 * padding;

	if (bitness == 8)
	{
		img->data8 = new uint8_t **[ezsize];
		img->data8[0] = new uint8_t *[ezsize * eysize];
		img->data8[0][0] = new uint8_t[4 * ezsize * eysize * exsize];
		for (i = 1; i < ezsize; i++)
		{
			img->data8[i] = img->data8[0] + i * eysize;
			img->data8[i][0] = img->data8[0][0] + 4 * i * exsize * eysize;
		}
		for (i = 0; i < ezsize; i++)
			for (j = 1; j < eysize; j++)
				img->data8[i][j] = img->data8[i][0] + 4 * j * exsize;

		img->data16 = NULL;
	}
	else if (bitness == 16)
	{
		img->data16 = new uint16_t **[ezsize];
		img->data16[0] = new uint16_t *[ezsize * eysize];
		img->data16[0][0] = new uint16_t[4 * ezsize * eysize * exsize];
		for (i = 1; i < ezsize; i++)
		{
			img->data16[i] = img->data16[0] + i * eysize;
			img->data16[i][0] = img->data16[0][0] + 4 * i * exsize * eysize;
		}
		for (i = 0; i < ezsize; i++)
			for (j = 1; j < eysize; j++)
				img->data16[i][j] = img->data16[i][0] + 4 * j * exsize;

		img->data8 = NULL;
	}
	else
	{
		ASTC_CODEC_INTERNAL_ERROR();
	}

	return img;
}

void initialize_image(astc_codec_image * img)
{
	int x, y, z;

	int exsize = img->xsize + 2 * img->padding;
	int eysize = img->ysize + 2 * img->padding;
	int ezsize = (img->zsize == 1) ? 1 : img->zsize + 2 * img->padding;

	if (img->data8)
	{
		for (z = 0; z < ezsize; z++)
			for (y = 0; y < eysize; y++)
				for (x = 0; x < exsize; x++)
				{
					img->data8[z][y][4 * x] = 0;
					img->data8[z][y][4 * x + 1] = 0;
					img->data8[z][y][4 * x + 2] = 0;
					img->data8[z][y][4 * x + 3] = 0xFF;
				}
	}
	else if (img->data16)
	{
		for (z = 0; z < ezsize; z++)
			for (y = 0; y < eysize; y++)
				for (x = 0; x < exsize; x++)
				{
					img->data16[z][y][4 * x] = 0;
					img->data16[z][y][4 * x + 1] = 0;
					img->data16[z][y][4 * x + 2] = 0;
					img->data16[z][y][4 * x + 3] = 0x3C00;
				}
	}
	else
	{
		ASTC_CODEC_INTERNAL_ERROR();
	}
}

// fill the padding area of the input-file buffer with clamp-to-edge data
// Done inefficiently, in that it will overwrite all the interior data at least once;
// this is not considered a problem, since this makes up a very small part of total
// running time.
void fill_image_padding_area(astc_codec_image * img)
{
	if (img->padding == 0)
		return;

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
	else						// if( bitness == 16 )
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

// conversion functions between the LNS representation and the FP16 representation.
float float_to_lns(float p)
{
	if (astc_isnan(p) || p <= 1.0f / 67108864.0f)
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
	float normfrac = frexp(p, &expo);
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

void imageblock_initialize_deriv_from_work_and_orig(
	imageblock* pb,
	int pixelcount
) {
	int i;

	const float *fptr = pb->orig_data;
	const float *wptr = pb->work_data;
	float *dptr = pb->deriv_data;

	for (i = 0; i < pixelcount; i++)
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
				rderiv = (1.0f / 32.0f);
			else if (rderiv > 33554432.0f)
				rderiv = 33554432.0f;

			if (gderiv < (1.0f / 32.0f))
				gderiv = (1.0f / 32.0f);
			else if (gderiv > 33554432.0f)
				gderiv = 33554432.0f;

			if (bderiv < (1.0f / 32.0f))
				bderiv = (1.0f / 32.0f);
			else if (bderiv > 33554432.0f)
				bderiv = 33554432.0f;

			dptr[0] = rderiv;
			dptr[1] = gderiv;
			dptr[2] = bderiv;
		}
		else
		{
			dptr[0] = 65535.0f;
			dptr[1] = 65535.0f;
			dptr[2] = 65535.0f;
		}

		// then compute derivatives for Alpha
		if (pb->alpha_lns[i])
		{
			float a = MAX(fptr[3], 6e-5f);
			float aderiv = (float_to_lns(a * 1.05f) - float_to_lns(a)) / (a * 0.05f);
			// the derivative may not actually take values smaller than 1/32 or larger than 2^25;
			// if it does, we clamp it.
			if (aderiv < (1.0f / 32.0f))
				aderiv = (1.0f / 32.0f);
			else if (aderiv > 33554432.0f)
				aderiv = 33554432.0f;

			dptr[3] = aderiv;
		}
		else
		{
			dptr[3] = 65535.0f;
		}

		fptr += 4;
		wptr += 4;
		dptr += 4;
	}
}

// helper function to initialize the work-data from the orig-data
void imageblock_initialize_work_from_orig(
	imageblock* pb,
	int pixelcount
) {
	int i;
	float *fptr = pb->orig_data;
	float *wptr = pb->work_data;

	for (i = 0; i < pixelcount; i++)
	{
		if (pb->rgb_lns[i])
		{
			wptr[0] = float_to_lns(fptr[0]);
			wptr[1] = float_to_lns(fptr[1]);
			wptr[2] = float_to_lns(fptr[2]);
		}
		else
		{
			wptr[0] = fptr[0] * 65535.0f;
			wptr[1] = fptr[1] * 65535.0f;
			wptr[2] = fptr[2] * 65535.0f;
		}

		if (pb->alpha_lns[i])
		{
			wptr[3] = float_to_lns(fptr[3]);
		}
		else
		{
			wptr[3] = fptr[3] * 65535.0f;
		}

		fptr += 4;
		wptr += 4;
	}

	imageblock_initialize_deriv_from_work_and_orig(pb, pixelcount);
}

// helper function to initialize the orig-data from the work-data
void imageblock_initialize_orig_from_work(
	imageblock* pb,
	int pixelcount
) {
	int i;
	float *fptr = pb->orig_data;
	float *wptr = pb->work_data;

	for (i = 0; i < pixelcount; i++)
	{
		if (pb->rgb_lns[i])
		{
			fptr[0] = sf16_to_float(lns_to_sf16((uint16_t) wptr[0]));
			fptr[1] = sf16_to_float(lns_to_sf16((uint16_t) wptr[1]));
			fptr[2] = sf16_to_float(lns_to_sf16((uint16_t) wptr[2]));
		}
		else
		{
			fptr[0] = sf16_to_float(unorm16_to_sf16((uint16_t) wptr[0]));
			fptr[1] = sf16_to_float(unorm16_to_sf16((uint16_t) wptr[1]));
			fptr[2] = sf16_to_float(unorm16_to_sf16((uint16_t) wptr[2]));
		}

		if (pb->alpha_lns[i])
		{
			fptr[3] = sf16_to_float(lns_to_sf16((uint16_t) wptr[3]));
		}
		else
		{
			fptr[3] = sf16_to_float(unorm16_to_sf16((uint16_t) wptr[3]));
		}

		fptr += 4;
		wptr += 4;
	}

	imageblock_initialize_deriv_from_work_and_orig(pb, pixelcount);
}

// fetch an imageblock from the input file.
void fetch_imageblock(
	const astc_codec_image* img,
	imageblock* pb,	// picture-block to initialize with image data
	const block_size_descriptor* bsd,
	// position in texture.
	int xpos,
	int ypos,
	int zpos,
	swizzlepattern swz
) {
	float *fptr = pb->orig_data;
	int xsize = img->xsize + 2 * img->padding;
	int ysize = img->ysize + 2 * img->padding;
	int zsize = (img->zsize == 1) ? 1 : img->zsize + 2 * img->padding;

	int x, y, z, i;

	pb->xpos = xpos;
	pb->ypos = ypos;
	pb->zpos = zpos;

	xpos += img->padding;
	ypos += img->padding;
	if (img->zsize > 1)
		zpos += img->padding;

	float data[6];
	data[4] = 0;
	data[5] = 1;

	if (img->data8)
	{
		for (z = 0; z < bsd->zdim; z++)
		{
			for (y = 0; y < bsd->ydim; y++)
			{
				for (x = 0; x < bsd->xdim; x++)
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

					int r = img->data8[zi][yi][4 * xi];
					int g = img->data8[zi][yi][4 * xi + 1];
					int b = img->data8[zi][yi][4 * xi + 2];
					int a = img->data8[zi][yi][4 * xi + 3];

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
	else if (img->data16)
	{
		for (z = 0; z < bsd->zdim; z++)
		{
			for (y = 0; y < bsd->ydim; y++)
			{
				for (x = 0; x < bsd->xdim; x++)
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

					int r = img->data16[zi][yi][4 * xi];
					int g = img->data16[zi][yi][4 * xi + 1];
					int b = img->data16[zi][yi][4 * xi + 2];
					int a = img->data16[zi][yi][4 * xi + 3];

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

	// perform sRGB-to-linear transform on input data, if requested.
	int pixelcount = bsd->texel_count;

	if (perform_srgb_transform)
	{
		fptr = pb->orig_data;
		for (i = 0; i < pixelcount; i++)
		{
			float r = fptr[0];
			float g = fptr[1];
			float b = fptr[2];

			if (r <= 0.04045f)
				r = r * (1.0f / 12.92f);
			else if (r <= 1)
				r = pow((r + 0.055f) * (1.0f / 1.055f), 2.4f);

			if (g <= 0.04045f)
				g = g * (1.0f / 12.92f);
			else if (g <= 1)
				g = pow((g + 0.055f) * (1.0f / 1.055f), 2.4f);

			if (b <= 0.04045f)
				b = b * (1.0f / 12.92f);
			else if (b <= 1)
				b = pow((b + 0.055f) * (1.0f / 1.055f), 2.4f);

			fptr[0] = r;
			fptr[1] = g;
			fptr[2] = b;

			fptr += 4;
		}
	}

	// collect color max-value, in order to determine whether to use LDR or HDR
	// interpolation.
	float max_red, max_green, max_blue, max_alpha;
	max_red = 0.0f;
	max_green = 0.0f;
	max_blue = 0.0f;
	max_alpha = 0.0f;

	fptr = pb->orig_data;
	for (i = 0; i < pixelcount; i++)
	{
		float r = fptr[0];
		float g = fptr[1];
		float b = fptr[2];
		float a = fptr[3];

		if (r > max_red)
			max_red = r;
		if (g > max_green)
			max_green = g;
		if (b > max_blue)
			max_blue = b;
		if (a > max_alpha)
			max_alpha = a;

		fptr += 4;
	}

	float max_rgb = MAX(max_red, MAX(max_green, max_blue));

	// use LNS if:
	// * RGB-maximum is less than 0.15
	// * RGB-maximum is greater than 1
	// * Alpha-maximum is greater than 1
	int rgb_lns = (max_rgb < 0.15f || max_rgb > 1.0f || max_alpha > 1.0f) ? 1 : 0;
	int alpha_lns = rgb_lns ? (max_alpha > 1.0f || max_alpha < 0.15f) : 0;

	// not yet though; for the time being, just obey the command line.
	rgb_lns = rgb_force_use_of_hdr;
	alpha_lns = alpha_force_use_of_hdr;

	// impose the choice on every pixel when encoding.
	for (i = 0; i < pixelcount; i++)
	{
		pb->rgb_lns[i] = rgb_lns;
		pb->alpha_lns[i] = alpha_lns;
		pb->nan_texel[i] = 0;
	}

	imageblock_initialize_work_from_orig(pb, pixelcount);
	update_imageblock_flags(pb, bsd->xdim, bsd->ydim, bsd->zdim);
}

void write_imageblock(
	astc_codec_image* img,
	const imageblock* pb,	// picture-block to initialize with image data. We assume that orig_data is valid.
	const block_size_descriptor* bsd,
	// position to write the block to
	int xpos,
	int ypos,
	int zpos,
	swizzlepattern swz
) {
	const float *fptr = pb->orig_data;
	const uint8_t *nptr = pb->nan_texel;
	int xsize = img->xsize;
	int ysize = img->ysize;
	int zsize = img->zsize;
	int x, y, z;

	float data[7];
	data[4] = 0.0f;
	data[5] = 1.0f;

	if (img->data8)
	{
		for (z = 0; z < bsd->zdim; z++)
		{
			for (y = 0; y < bsd->ydim; y++)
			{
				for (x = 0; x < bsd->xdim; x++)
				{
					int xi = xpos + x;
					int yi = ypos + y;
					int zi = zpos + z;

					if (xi >= 0 && yi >= 0 && zi >= 0 && xi < xsize && yi < ysize && zi < zsize)
					{
						if (*nptr)
						{
							// NaN-pixel, but we can't display it. Display purple instead.
							img->data8[zi][yi][4 * xi] = 0xFF;
							img->data8[zi][yi][4 * xi + 1] = 0x00;
							img->data8[zi][yi][4 * xi + 2] = 0xFF;
							img->data8[zi][yi][4 * xi + 3] = 0xFF;
						}
						else
						{
							// apply swizzle
							if (perform_srgb_transform)
							{
								float r = fptr[0];
								float g = fptr[1];
								float b = fptr[2];

								if (r <= 0.0031308f)
									r = r * 12.92f;
								else if (r <= 1)
									r = 1.055f * powf(r, (1.0f / 2.4f)) - 0.055f;

								if (g <= 0.0031308f)
									g = g * 12.92f;
								else if (g <= 1)
									g = 1.055f * powf(g, (1.0f / 2.4f)) - 0.055f;

								if (b <= 0.0031308f)
									b = b * 12.92f;
								else if (b <= 1)
									b = 1.055f * powf(b, (1.0f / 2.4f)) - 0.055f;

								data[0] = r;
								data[1] = g;
								data[2] = b;
							}
							else
							{
								float r = fptr[0];
								float g = fptr[1];
								float b = fptr[2];

								data[0] = r;
								data[1] = g;
								data[2] = b;
							}
							data[3] = fptr[3];

							float xcoord = (data[0] * 2.0f) - 1.0f;
							float ycoord = (data[3] * 2.0f) - 1.0f;
							float zcoord = 1.0f - xcoord * xcoord - ycoord * ycoord;
							if (zcoord < 0.0f)
								zcoord = 0.0f;
							data[6] = (sqrtf(zcoord) * 0.5f) + 0.5f;

							// clamp to [0,1]
							if (data[0] > 1.0f)
								data[0] = 1.0f;
							if (data[1] > 1.0f)
								data[1] = 1.0f;
							if (data[2] > 1.0f)
								data[2] = 1.0f;
							if (data[3] > 1.0f)
								data[3] = 1.0f;

							// pack the data
							int ri = static_cast < int >(floor(data[swz.r] * 255.0f + 0.5f));
							int gi = static_cast < int >(floor(data[swz.g] * 255.0f + 0.5f));
							int bi = static_cast < int >(floor(data[swz.b] * 255.0f + 0.5f));
							int ai = static_cast < int >(floor(data[swz.a] * 255.0f + 0.5f));

							img->data8[zi][yi][4 * xi] = ri;
							img->data8[zi][yi][4 * xi + 1] = gi;
							img->data8[zi][yi][4 * xi + 2] = bi;
							img->data8[zi][yi][4 * xi + 3] = ai;
						}
					}
					fptr += 4;
					nptr++;
				}
			}
		}
	}
	else if (img->data16)
	{
		for (z = 0; z < bsd->zdim; z++)
		{
			for (y = 0; y < bsd->ydim; y++)
			{
				for (x = 0; x < bsd->xdim; x++)
				{
					int xi = xpos + x;
					int yi = ypos + y;
					int zi = zpos + z;

					if (xi >= 0 && yi >= 0 && zi >= 0 && xi < xsize && yi < ysize && zi < zsize)
					{
						if (*nptr)
						{
							img->data16[zi][yi][4 * xi] = 0xFFFF;
							img->data16[zi][yi][4 * xi + 1] = 0xFFFF;
							img->data16[zi][yi][4 * xi + 2] = 0xFFFF;
							img->data16[zi][yi][4 * xi + 3] = 0xFFFF;
						}

						else
						{
							// apply swizzle
							if (perform_srgb_transform)
							{
								float r = fptr[0];
								float g = fptr[1];
								float b = fptr[2];

								if (r <= 0.0031308f)
									r = r * 12.92f;
								else if (r <= 1)
									r = 1.055f * powf(r, (1.0f / 2.4f)) - 0.055f;
								if (g <= 0.0031308f)
									g = g * 12.92f;
								else if (g <= 1)
									g = 1.055f * powf(g, (1.0f / 2.4f)) - 0.055f;
								if (b <= 0.0031308f)
									b = b * 12.92f;
								else if (b <= 1)
									b = 1.055f * powf(b, (1.0f / 2.4f)) - 0.055f;

								data[0] = r;
								data[1] = g;
								data[2] = b;
							}
							else
							{
								data[0] = fptr[0];
								data[1] = fptr[1];
								data[2] = fptr[2];
							}
							data[3] = fptr[3];

							float xN = (data[0] * 2.0f) - 1.0f;
							float yN = (data[3] * 2.0f) - 1.0f;
							float zN = 1.0f - xN * xN - yN * yN;
							if (zN < 0.0f)
								zN = 0.0f;
							data[6] = (sqrtf(zN) * 0.5f) + 0.5f;

							int r = float_to_sf16(data[swz.r], SF_NEARESTEVEN);
							int g = float_to_sf16(data[swz.g], SF_NEARESTEVEN);
							int b = float_to_sf16(data[swz.b], SF_NEARESTEVEN);
							int a = float_to_sf16(data[swz.a], SF_NEARESTEVEN);
							img->data16[zi][yi][4 * xi] = r;
							img->data16[zi][yi][4 * xi + 1] = g;
							img->data16[zi][yi][4 * xi + 2] = b;
							img->data16[zi][yi][4 * xi + 3] = a;
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
   The updating is done based on work_data, not orig_data.
*/
void update_imageblock_flags(
	imageblock* pb,
	int xdim,
	int ydim,
	int zdim
) {
	int i;
	float red_min = 1e38f, red_max = -1e38f;
	float green_min = 1e38f, green_max = -1e38f;
	float blue_min = 1e38f, blue_max = -1e38f;
	float alpha_min = 1e38f, alpha_max = -1e38f;

	int texels_per_block = xdim * ydim * zdim;

	int grayscale = 1;

	for (i = 0; i < texels_per_block; i++)
	{
		float red = pb->work_data[4 * i];
		float green = pb->work_data[4 * i + 1];
		float blue = pb->work_data[4 * i + 2];
		float alpha = pb->work_data[4 * i + 3];
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
			grayscale = 0;
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

/*
	Main image loader function.

	We have specialized loaders for DDS, KTX and HTGA; for other formats, we use stb_image.
	This image loader will choose one based on filename.
*/
astc_codec_image *astc_codec_load_image(
	const char* input_filename,
	int padding,
	int* load_result
) {
	#define LOAD_HTGA 0
	#define LOAD_KTX 1
	#define LOAD_DDS 2
	#define LOAD_STB_IMAGE 3

	// check the ending of the input filename
	int load_fileformat = LOAD_STB_IMAGE;
	size_t filename_len = strlen(input_filename);

	const char *eptr = input_filename + filename_len - 5;
	if (eptr > input_filename && (strcmp(eptr, ".htga") == 0 || strcmp(eptr, ".HTGA") == 0))
		load_fileformat = LOAD_HTGA;
	eptr = input_filename + filename_len - 4;
	if (eptr > input_filename && (strcmp(eptr, ".ktx") == 0 || strcmp(eptr, ".KTX") == 0))
		load_fileformat = LOAD_KTX;
	if (eptr > input_filename && (strcmp(eptr, ".dds") == 0 || strcmp(eptr, ".DDS") == 0))
		load_fileformat = LOAD_DDS;

	// OpenEXR support: call exr_to_htga to convert from EXR to HTGA.
	char htga_load_filename[300];
	int load_exr = 0;
	if (eptr > input_filename && (strcmp(eptr, ".exr") == 0 || strcmp(eptr, ".EXR") == 0))
	{
		// don't support filenames longer than 250 characters; this way, we
		// cannot get a buffer overflow from the sprintfs below.
		if (filename_len > 250)
		{
			*load_result = -1;
			return NULL;
		}

		char exr_to_htga_command[550];
		sprintf(htga_load_filename, "%s.htga", input_filename);
		sprintf(exr_to_htga_command, "exr_to_htga -q %s %s", input_filename, htga_load_filename);

		int retval = system(exr_to_htga_command);
		if (retval != 0)
		{
			printf("Failed to run exr_to_htga to convert input .exr file.\n");
			exit(1);
		}
		input_filename = htga_load_filename;
		load_fileformat = LOAD_HTGA;
		load_exr = 1;
	}

	astc_codec_image *input_image;

	switch (load_fileformat)
	{
	case LOAD_KTX:
		input_image = load_ktx_uncompressed_image(input_filename, padding, load_result);
		break;
	case LOAD_DDS:
		input_image = load_dds_uncompressed_image(input_filename, padding, load_result);
		break;
	case LOAD_HTGA:
		input_image = load_tga_image(input_filename, padding, load_result);
		break;
	case LOAD_STB_IMAGE:
		input_image = load_image_with_stb(input_filename, padding, load_result);
		break;
	default:
		ASTC_CODEC_INTERNAL_ERROR();
	}

	if (load_exr)
		astc_codec_unlink(htga_load_filename);

	return input_image;
}

int get_output_filename_enforced_bitness(const char *output_filename)
{
	if (output_filename == NULL)
		return -1;

	size_t filename_len = strlen(output_filename);
	const char *eptr = output_filename + filename_len - 5;

	if (eptr > output_filename && (strcmp(eptr, ".htga") == 0 || strcmp(eptr, ".HTGA") == 0))
	{
		return 16;
	}

	eptr = output_filename + filename_len - 4;
	if (eptr > output_filename && (strcmp(eptr, ".tga") == 0 || strcmp(eptr, ".TGA") == 0))
	{
		return 8;
	}
	if (eptr > output_filename && (strcmp(eptr, ".exr") == 0 || strcmp(eptr, ".EXR") == 0))
	{
		return 16;
	}

	// file formats that don't match any of the templates above are capable of accommodating
	// both 8-bit and 16-bit data (DDS, KTX)
	return -1;
}

int astc_codec_store_image(
	const astc_codec_image* output_image,
	const char* output_filename,
	int bitness,
	const char**format_string
) {
	#define STORE_TGA 0
	#define STORE_HTGA 1
	#define STORE_KTX 2
	#define STORE_DDS 3
	#define STORE_EXR 4

	size_t filename_len = strlen(output_filename);

	int store_fileformat = STORE_TGA;
	const char *eptr = output_filename + filename_len - 5;
	if (eptr > output_filename && (strcmp(eptr, ".htga") == 0 || strcmp(eptr, ".HTGA") == 0))
	{
		store_fileformat = STORE_HTGA;
	}
	eptr = output_filename + filename_len - 4;

	if (eptr > output_filename && (strcmp(eptr, ".ktx") == 0 || strcmp(eptr, ".KTX") == 0))
	{
		store_fileformat = STORE_KTX;
	}

	if (eptr > output_filename && (strcmp(eptr, ".dds") == 0 || strcmp(eptr, ".DDS") == 0))
	{
		store_fileformat = STORE_DDS;
	}

	if (eptr > output_filename && (strcmp(eptr, ".exr") == 0 || strcmp(eptr, ".EXR") == 0))
	{
		store_fileformat = STORE_EXR;
	}

	if (store_fileformat == STORE_TGA && bitness == 16)
		store_fileformat = STORE_HTGA;

	// guard against OpenEXR files with too-long names
	if (store_fileformat == STORE_EXR && filename_len > 250)
	{
		*format_string = "EXR";
		return -1;
	}

	char htga_output_filename[300];
	char htga_output_command[550];
	int system_retval;

	int store_result = -1;
	switch (store_fileformat)
	{
	case STORE_TGA:
	case STORE_HTGA:
		*format_string = bitness == 16 ? "HTGA" : "TGA";
		store_result = store_tga_image(output_image, output_filename, bitness);
		break;
	case STORE_KTX:
		*format_string = "KTX";
		store_result = store_ktx_uncompressed_image(output_image, output_filename, bitness);
		break;
	case STORE_DDS:
		*format_string = "DDS";
		store_result = store_dds_uncompressed_image(output_image, output_filename, bitness);
		break;
	case STORE_EXR:
		*format_string = "EXR";
		sprintf(htga_output_filename, "%s.htga", output_filename);
		store_result = store_tga_image(output_image, htga_output_filename, 16);
		sprintf(htga_output_command, "exr_to_htga -e %s %s", htga_output_filename, output_filename);
		system_retval = system(htga_output_command);
		astc_codec_unlink(htga_output_filename);
		if (system_retval != 0)
			store_result = -99;
		break;
	default:
		ASTC_CODEC_INTERNAL_ERROR();
	};

	return store_result;
}
