// ----------------------------------------------------------------------------
//  This confidential and proprietary software may be used only as authorised
//  by a licensing agreement from Arm Limited.
//      (C) COPYRIGHT 2011-2020 Arm Limited, ALL RIGHTS RESERVED
//  The entire notice above must be reproduced on all authorised copies and
//  copies may only be made to the extent permitted by a licensing agreement
//  from Arm Limited.
// ----------------------------------------------------------------------------

/**
 * @brief Functions for computing image error metrics.
 */

#include "astc_codec_internals.h"

#include <cstdio>

/**
 * @brief An accumulator using Kahan compensated floating-point summation.
 *
 * This method keeps higher precision than direct summation by keeping track of
 * the error compensation factor, c, which can be added into the next
 * calculation. This allows single precision floats to be used in places that
 * would otherwise need double precision, which is useful when vectorizing.
 */
class kahan_accum4
{
public:
	/** The running sum. */
	float4 sum;

	/** The current compensation factor. */
	float4 comp;

	/**
	 * @brief Create a new Kahan accumulator
	 */
	kahan_accum4() {
		sum = float4(0.0f, 0.0f, 0.0f, 0.0f);
		comp = float4(0.0f, 0.0f, 0.0f, 0.0f);
	}
};

/**
 * @brief The incremental addition operator for Kahan summation.
 *
 * @param val The Kahan accumulator to increment
 * @param inc The increment to apply
 *
 * @return The updated accumulator
 */
static kahan_accum4 &operator+=(kahan_accum4 &val, float4 inc)
{
	float4 y = inc - val.comp;
	float4 t = val.sum + y;
	val.comp = (t - val.sum) - y;
	val.sum = t;
	return val;
}

/**
 * @brief Clamp value in range [0.0f, 65504.0f], with NaNs converted to zero.
 *
 * @param val The value to clamp
 *
 * @return The clamped value
 */
static float clamp(float val)
{
	// Do not reorder these, correct NaN handling relies on the fact that
	// any comparison with NaN returns false so will fall-though to the 0.0f.
	if (val > 65504.0f) return 65504.0f;
	if (val > 0.0f) return val;
	return 0.0f;
}

/**
 * @brief Logarithm, linearized from 2^-14.
 *
 * @param val The value to log2
 *
 * @return log2(val)
 */
static float xlog2(float val)
{
	if (val >= 0.00006103515625f)
		return log(val) * 1.44269504088896340735f;	// log(x)/log(2)
	else
		return -15.44269504088896340735f + val * 23637.11554992477646609062f;
}

/**
 * @brief mPSNR tonemapping operator for HDR images.
 *
 * @param val The color value to tone map
 * @param fstop The exposure fstop; should be in range [-125, 125]
 *
 * @return The mapped color value in [0.0f, 255.0f] range
 */
static float mpsnr_operator(float val, int fstop)
{
	if32 p;
	p.u = 0x3f800000 + (fstop << 23);  // 0x3f800000 is 1.0f
	val *= p.f;
	val = pow(val, (1.0f / 2.2f));
	val *= 255.0f;

	// Do not reorder these, correct NaN handling relies on the fact that
	// any comparison with NaN returns false so will fall-though to the 0.0f.
	if (val > 255.0f) return 255.0f;
	if (val > 0.0f) return val;
	return 0.0f;
}

/**
 * @brief mPSNR difference between two values.
 *
 * Differences are given as "val1 - val2".
 *
 * @param val1     The first color value
 * @param val2     The second color value
 * @param fstop_lo The low exposure fstop; should be in range [-125, 125]
 * @param fstop_hi The high exposure fstop; should be in range [-125, 125]
 *
 * @return The summed mPSNR difference across all active fstop levels
 */
static float mpsnr_sumdiff(float val1, float val2, int fstop_lo, int fstop_hi)
{
	float summa = 0.0f;
	for (int i = fstop_lo; i <= fstop_hi; i++)
	{
		float mval1 = mpsnr_operator(val1, i);
		float mval2 = mpsnr_operator(val2, i);
		float mdiff = mval1 - mval2;
		summa += mdiff * mdiff;
	}
	return summa;
}

/* Public function, see header file for detailed documentation */
void compute_error_metrics(
	int compute_hdr_metrics,
	int input_components,
	const astc_codec_image * img1,
	const astc_codec_image * img2,
	int fstop_lo,
	int fstop_hi,
	int show_psnr
) {
	static int channelmasks[5] = { 0x00, 0x07, 0x0C, 0x07, 0x0F };
	int channelmask = channelmasks[input_components];

	kahan_accum4 errorsum;
	kahan_accum4 alpha_scaled_errorsum;
	kahan_accum4 log_errorsum;
	kahan_accum4 mpsnr_errorsum;

	int xsize = MIN(img1->xsize, img2->xsize);
	int ysize = MIN(img1->ysize, img2->ysize);
	int zsize = MIN(img1->zsize, img2->zsize);

	if (img1->xsize != img2->xsize ||
	    img1->ysize != img2->ysize ||
	    img1->zsize != img2->zsize)
	{
		printf("Warning: Only intersection of images will be compared:\n"
		       "  Image 1: %dx%dx%d\n"
		       "  Image 2: %dx%dx%d\n",
		       img1->xsize, img1->ysize, img1->zsize,
		       img2->xsize, img2->ysize, img2->zsize);
	}

	// HDR metric computation across all fstops is slow, so warn the user
	if (compute_hdr_metrics)
	{
		printf("Computing error metrics ... ");
		fflush(stdout);
	}

	int img1pad = img1->padding;
	int img2pad = img2->padding;
	float rgb_peak = 0.0f;

	for (int z = 0; z < zsize; z++)
	{
		for (int y = 0; y < ysize; y++)
		{
			int ze1 = (img1->zsize == 1) ? z : z + img1pad;
			int ze2 = (img2->zsize == 1) ? z : z + img2pad;

			int ye1 = y + img1pad;
			int ye2 = y + img2pad;

			for (int x = 0; x < xsize; x++)
			{
				float4 color1;
				float4 color2;

				int xe1 = 4 * x + 4 * img1pad;
				int xe2 = 4 * x + 4 * img2pad;

				if (img1->data8)
				{
					color1 = float4(
					    img1->data8[ze1][ye1][xe1]     * (1.0f / 255.0f),
					    img1->data8[ze1][ye1][xe1 + 1] * (1.0f / 255.0f),
					    img1->data8[ze1][ye1][xe1 + 2] * (1.0f / 255.0f),
					    img1->data8[ze1][ye1][xe1 + 3] * (1.0f / 255.0f));
				}
				else
				{
					color1 = float4(
					    clamp(sf16_to_float(img1->data16[ze1][ye1][xe1])),
					    clamp(sf16_to_float(img1->data16[ze1][ye1][xe1 + 1])),
					    clamp(sf16_to_float(img1->data16[ze1][ye1][xe1 + 2])),
					    clamp(sf16_to_float(img1->data16[ze1][ye1][xe1 + 3])));
				}

				if (img2->data8)
				{
					color2 = float4(
					    img2->data8[ze2][ye2][xe2]     * (1.0f / 255.0f),
					    img2->data8[ze2][ye2][xe2 + 1] * (1.0f / 255.0f),
					    img2->data8[ze2][ye2][xe2 + 2] * (1.0f / 255.0f),
					    img2->data8[ze2][ye2][xe2 + 3] * (1.0f / 255.0f));
				}
				else
				{
					color2 = float4(
					    clamp(sf16_to_float(img2->data16[ze2][ye2][xe2])),
					    clamp(sf16_to_float(img2->data16[ze2][ye2][xe2 + 1])),
					    clamp(sf16_to_float(img2->data16[ze2][ye2][xe2 + 2])),
					    clamp(sf16_to_float(img2->data16[ze2][ye2][xe2 + 3])));
				}

				rgb_peak = MAX(MAX(color1.x, color1.y),
				               MAX(color1.z, rgb_peak));

				float4 diffcolor = color1 - color2;
				errorsum += diffcolor * diffcolor;

				float4 alpha_scaled_diffcolor = float4(
				    diffcolor.x * color1.w,
				    diffcolor.y * color1.w,
				    diffcolor.z * color1.w,
				    diffcolor.w);

				alpha_scaled_errorsum += alpha_scaled_diffcolor * \
				                         alpha_scaled_diffcolor;

				if (compute_hdr_metrics)
				{
					float4 log_input_color1 = float4(
					    xlog2(color1.x),
					    xlog2(color1.y),
					    xlog2(color1.z),
					    xlog2(color1.w));

					float4 log_input_color2 = float4(
					    xlog2(color2.x),
					    xlog2(color2.y),
					    xlog2(color2.z),
					    xlog2(color2.w));

					float4 log_diffcolor = log_input_color1 - log_input_color2;

					log_errorsum += log_diffcolor * log_diffcolor;

					float4 mpsnr_error = float4(
					    mpsnr_sumdiff(color1.x, color2.x, fstop_lo, fstop_hi),
					    mpsnr_sumdiff(color1.y, color2.y, fstop_lo, fstop_hi),
					    mpsnr_sumdiff(color1.z, color2.z, fstop_lo, fstop_hi),
					    mpsnr_sumdiff(color1.w, color2.w, fstop_lo, fstop_hi));

					mpsnr_errorsum += mpsnr_error;
				}
			}
		}
	}

	if (compute_hdr_metrics)
	{
		printf("done\n");
	}

	float pixels = (float)(xsize * ysize * zsize);
	float num = 0.0f;
	float alpha_num = 0.0f;
	float log_num = 0.0f;
	float mpsnr_num = 0.0f;
	float samples = 0.0f;

	if (channelmask & 1)
	{
		num += errorsum.sum.x;
		alpha_num += alpha_scaled_errorsum.sum.x;
		log_num += log_errorsum.sum.x;
		mpsnr_num += mpsnr_errorsum.sum.x;
		samples += pixels;
	}

	if (channelmask & 2)
	{
		num += errorsum.sum.y;
		alpha_num += alpha_scaled_errorsum.sum.y;
		log_num += log_errorsum.sum.y;
		mpsnr_num += mpsnr_errorsum.sum.y;
		samples += pixels;
	}

	if (channelmask & 4)
	{
		num += errorsum.sum.z;
		alpha_num += alpha_scaled_errorsum.sum.z;
		log_num += log_errorsum.sum.z;
		mpsnr_num += mpsnr_errorsum.sum.z;
		samples += pixels;
	}

	if (channelmask & 8)
	{
		num += errorsum.sum.w;
		alpha_num += alpha_scaled_errorsum.sum.w;
		samples += pixels;
	}

	float denom = samples;
	float stopcount = (float)(fstop_hi - fstop_lo + 1);
	float mpsnr_denom = pixels * 3.0f * stopcount * 255.0f * 255.0f;

	float psnr;
	if (num == 0.0f)
		psnr = 999.0f;
	else
		psnr = 10.0f * log10(denom / num);

	float rgb_psnr = psnr;
	if (show_psnr)
	{
		if (channelmask & 8)
		{
			printf("PSNR (LDR-RGBA): %.6f dB\n", psnr);

			float alpha_psnr;
			if (alpha_num == 0.0f)
				alpha_psnr = 999.0f;
			else
				alpha_psnr = 10.0f * log10(denom / alpha_num);
			printf("Alpha-Weighted PSNR: %.6f dB\n", alpha_psnr);

			float rgb_num = errorsum.sum.x + errorsum.sum.y + errorsum.sum.z;
			if (rgb_num == 0.0f)
				rgb_psnr = 999.0f;
			else
				rgb_psnr = 10.0f * log10( pixels * 3.0f / rgb_num);
			printf("PSNR (LDR-RGB): %.6f dB\n", rgb_psnr);
		}
		else
			printf("PSNR (LDR-RGB): %.6f dB\n", psnr);

		if (compute_hdr_metrics)
		{
			printf("Color peak value: %f\n", rgb_peak);
			printf("PSNR (RGB normalized to peak): %f dB\n",
			       rgb_psnr + 20.0f * log10(rgb_peak));

			float mpsnr;
			if (mpsnr_num == 0.0f)
				mpsnr = 999.0f;
			else
				mpsnr = 10.0f * log10(mpsnr_denom / mpsnr_num);
			printf("mPSNR (RGB) [fstops: %+d to %+d] : %.6f dB\n",
			       fstop_lo, fstop_hi, mpsnr);

			float logrmse = sqrt(log_num / pixels);
			printf("LogRMSE (RGB): %.6f\n", logrmse);
		}
	}
}
