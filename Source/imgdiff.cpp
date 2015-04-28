/*----------------------------------------------------------------------------*/
/**
 *	This confidential and proprietary software may be used only as
 *	authorised by a licensing agreement from ARM Limited
 *	(C) COPYRIGHT 2011-2012 ARM Limited
 *	ALL RIGHTS RESERVED
 *
 *	The entire notice above must be reproduced on all authorised
 *	copies and copies may only be made to the extent permitted
 *	by a licensing agreement from ARM Limited.
 *
 *	@brief	Program to load two TGA images and compute the difference between
 *			them (psnr)
 */
/*----------------------------------------------------------------------------*/

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct image
{
	uint8_t **data;
	int xsize;
	int ysize;
	int channels;
};

/* 
   given a TGA filename, read in a TGA file and create test-vectors from it. */

struct tga_header
{
	uint8_t identsize;
	uint8_t colormaptype;
	uint8_t imagetype;
	uint8_t dummied[5];
	uint16_t xstart;
	uint16_t ystart;
	uint16_t xsize;
	uint16_t ysize;
	uint8_t bitsperpixel;
	uint8_t descriptor;
};



image *allocate_input_image_space(int xsize, int ysize)
{
	int i;
	image *img = new image;

	img->data = new uint8_t *[ysize];
	img->data[0] = new uint8_t[ysize * xsize * 4];
	for (i = 1; i < ysize; i++)
		img->data[i] = img->data[0] + i * xsize * 4;
	img->xsize = xsize;
	img->ysize = ysize;
	img->channels = 0;
	return img;
}


/* 
   return image object on load, NULL on failed load */

image *load_tga_image(const char *tga_filename)
{
	int x, y;
	int i;

	FILE *f = fopen(tga_filename, "rb");
	if (!f)
		return NULL;

	tga_header hdr;
	int bytes_read = fread(&hdr, 1, 18, f);
	if (bytes_read != 18)
	{
		fclose(f);
		return NULL;
	}
	if (hdr.colormaptype != 0)
	{
		fclose(f);
		return NULL;
	}

	// do a quick test for RLE-pictures so that we reject them
	if ((hdr.imagetype == 10 && hdr.bitsperpixel == 32)
		|| (hdr.imagetype == 10 && hdr.bitsperpixel == 24) || (hdr.imagetype == 11 && hdr.bitsperpixel == 16) || (hdr.imagetype == 11 && hdr.bitsperpixel == 8))
	{
		fclose(f);
		return NULL;
	}

	// support 4 formats (non-RLE only):
	// 8-bit grayscale
	// 8-bit grayscale + 8-bit alpha
	// RGB 8:8:8
	// RGBA 8:8:8:8

	if (!(hdr.imagetype == 2 && hdr.bitsperpixel == 32)
		&& !(hdr.imagetype == 2 && hdr.bitsperpixel == 24) && !(hdr.imagetype == 3 && hdr.bitsperpixel == 16) && !(hdr.imagetype == 3 && hdr.bitsperpixel == 8))
	{
		fclose(f);
		return NULL;
	}

	if (hdr.identsize != 0)		// skip ID field if it present.
		fseek(f, hdr.identsize, SEEK_CUR);

	int bytesperpixel = hdr.bitsperpixel / 8;


	// OK, it seems we have a legit TGA file of a format we understand.
	// Now, let's read it.

	uint8_t **row_pointers = new uint8_t *[hdr.ysize];
	row_pointers[0] = new uint8_t[hdr.xsize * hdr.ysize * bytesperpixel];
	for (i = 1; i < hdr.ysize; i++)
		row_pointers[i] = row_pointers[0] + hdr.xsize * bytesperpixel * i;

	int bytestoread = hdr.xsize * hdr.ysize * bytesperpixel;
	bytes_read = fread(row_pointers[0], 1, bytestoread, f);
	fclose(f);
	if (bytes_read != bytestoread)
	{
		delete[]row_pointers[0];
		delete[]row_pointers;
		return NULL;
	}


	// OK, at this point, we can expand the image data to RGBA.
	int ysize = hdr.ysize;
	int xsize = hdr.xsize;

	image *img = allocate_input_image_space(xsize, ysize);
	for (y = 0; y < ysize; y++)
	{
		switch (bytesperpixel)
		{
		case 1:				// single-component, treated as Luminance
			for (x = 0; x < xsize; x++)
			{
				img->data[y][4 * x] = row_pointers[y][x];
				img->data[y][4 * x + 1] = row_pointers[y][x];
				img->data[y][4 * x + 2] = row_pointers[y][x];
				img->data[y][4 * x + 3] = 0xFF;
			}
			break;
		case 2:				// two-component, treated as Luminance-Alpha
			for (x = 0; x < xsize; x++)
			{
				img->data[y][4 * x] = row_pointers[y][2 * x];
				img->data[y][4 * x + 1] = row_pointers[y][2 * x];
				img->data[y][4 * x + 2] = row_pointers[y][2 * x];
				img->data[y][4 * x + 3] = row_pointers[y][2 * x + 1];
			}
			break;
		case 3:				// three-component, treated as RGB
			for (x = 0; x < xsize; x++)
			{
				img->data[y][4 * x] = row_pointers[y][3 * x + 2];	// tga uses BGR, we use RGB
				img->data[y][4 * x + 1] = row_pointers[y][3 * x + 1];
				img->data[y][4 * x + 2] = row_pointers[y][3 * x];
				img->data[y][4 * x + 3] = 0xFF;
			}
			break;
		case 4:				// three-component, treated as RGB
			for (x = 0; x < xsize; x++)
			{
				img->data[y][4 * x] = row_pointers[y][4 * x + 2];	// tga uses BGR, we use RGB
				img->data[y][4 * x + 1] = row_pointers[y][4 * x + 1];
				img->data[y][4 * x + 2] = row_pointers[y][4 * x];
				img->data[y][4 * x + 3] = row_pointers[y][4 * x + 3];
			}
			break;
		}
	}

	delete[]row_pointers[0];
	delete[]row_pointers;

	img->channels = bytesperpixel;
	return img;
}




void compute_psnr(const image * m1, const image * m2)
{
	int x, y;
	int channelmask;
	if ((m1->channels == 1 || m1->channels == 3) && (m2->channels == 1 || m2->channels == 3))
	{
		channelmask = 0x7;		// luminance or RGB case
		printf("Comparing as RGB images\n");
	}
	else if (m1->channels == 2 && m2->channels == 2)
	{
		channelmask = 0xC;		// luminance-alpha case
		printf("Comparing as Lum-Alpha images\n");
	}
	else
	{
		channelmask = 0xF;
		printf("Comparing as RGBA images\n");
	}

	uint64_t rsum = 0;
	uint64_t gsum = 0;
	uint64_t bsum = 0;
	uint64_t asum = 0;
	for (y = 0; y < m1->ysize; y++)
		for (x = 0; x < m1->xsize; x++)
		{
			int rdiff = m1->data[y][4 * x] - m2->data[y][4 * x];
			int gdiff = m1->data[y][4 * x + 1] - m2->data[y][4 * x + 1];
			int bdiff = m1->data[y][4 * x + 2] - m2->data[y][4 * x + 2];
			int adiff = m1->data[y][4 * x + 3] - m2->data[y][4 * x + 3];
			rsum += rdiff * rdiff;
			gsum += gdiff * gdiff;
			bsum += bdiff * bdiff;
			asum += adiff * adiff;
		}
	uint64_t limit = m1->ysize * m1->xsize * (uint64_t) (255 * 255);

	uint64_t num = 0;
	uint64_t denom = 0;
	if (channelmask & 1)
	{
		num += rsum;
		denom += limit;
	}
	if (channelmask & 2)
	{
		num += gsum;
		denom += limit;
	}
	if (channelmask & 4)
	{
		num += bsum;
		denom += limit;
	}
	if (channelmask & 8)
	{
		num += asum;
		denom += limit;
	}

	if (num == 0)
		printf("Images are identical\n");
	else
	{
		double psnr = 10.0 * log10((double)denom / (double)num);
		printf("PSNR: %f\n", psnr);
	}
}




int main(int argc, char **argv)
{
	if (argc < 3)
	{
		printf("Image comparison tool\n" "Usage:\n" "    %s <image1> <image2>\n", argv[0]);
		exit(1);
	}
	image *m1 = load_tga_image(argv[1]);
	image *m2 = load_tga_image(argv[2]);

	if (!m1)
	{
		printf("Failed to load image %s\n", argv[1]);
		exit(1);
	}
	if (!m2)
	{
		printf("Failed to load image %s\n", argv[2]);
		exit(1);
	}
	if (m1->xsize != m2->xsize || m1->ysize != m2->ysize)
	{
		printf("Image dimension mismatch:\n" "%s: %dx%d   %s: %dx%d\n", argv[1], m1->xsize, m1->ysize, argv[2], m2->xsize, m2->ysize);
		exit(1);
	}

	compute_psnr(m1, m2);
}
