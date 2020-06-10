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
 * @brief Functions for computing image error metrics.
 */

#include <cstdio>

#include "astcenccli_internal.h"

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
static kahan_accum4 &operator+=(
	kahan_accum4 &val,
	float4 inc
) {
	float4 y = inc - val.comp;
	float4 t = val.sum + y;
	val.comp = (t - val.sum) - y;
	val.sum = t;
	return val;
}

/**
 * @brief mPSNR tonemapping operator for HDR images.
 *
 * @param val The color value to tone map
 * @param fstop The exposure fstop; should be in range [-125, 125]
 *
 * @return The mapped color value in [0.0f, 255.0f] range
 */
static float mpsnr_operator(
	float val,
	int fstop
) {
	if32 p;
	p.u = 0x3f800000 + (fstop << 23);  // 0x3f800000 is 1.0f
	val *= p.f;
	val = powf(val, (1.0f / 2.2f));
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
static float mpsnr_sumdiff(
	float val1,
	float val2,
	int fstop_lo,
	int fstop_hi
) {
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
	const astcenc_image* img1,
	const astcenc_image* img2,
	int fstop_lo,
	int fstop_hi
) {
	static int channelmasks[5] = { 0x00, 0x07, 0x0C, 0x07, 0x0F };
	int channelmask = channelmasks[input_components];

	kahan_accum4 errorsum;
	kahan_accum4 alpha_scaled_errorsum;
	kahan_accum4 log_errorsum;
	kahan_accum4 mpsnr_errorsum;

	unsigned int dim_x = MIN(img1->dim_x, img2->dim_x);
	unsigned int dim_y = MIN(img1->dim_y, img2->dim_y);
	unsigned int dim_z = MIN(img1->dim_z, img2->dim_z);

	if (img1->dim_x != img2->dim_x ||
	    img1->dim_y != img2->dim_y ||
	    img1->dim_z != img2->dim_z)
	{
		printf("WARNING: Only intersection of images will be compared:\n"
		       "  Image 1: %dx%dx%d\n"
		       "  Image 2: %dx%dx%d\n",
		       img1->dim_x, img1->dim_y, img1->dim_z,
		       img2->dim_x, img2->dim_y, img2->dim_z);
	}

	int img1pad = img1->dim_pad;
	int img2pad = img2->dim_pad;
	float rgb_peak = 0.0f;

	for (unsigned int z = 0; z < dim_z; z++)
	{
		for (unsigned int y = 0; y < dim_y; y++)
		{
			int ze1 = (img1->dim_z == 1) ? z : z + img1pad;
			int ze2 = (img2->dim_z == 1) ? z : z + img2pad;

			int ye1 = y + img1pad;
			int ye2 = y + img2pad;

			for (unsigned int x = 0; x < dim_x; x++)
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
					    astc::clamp64Kf(sf16_to_float(img1->data16[ze1][ye1][xe1])),
					    astc::clamp64Kf(sf16_to_float(img1->data16[ze1][ye1][xe1 + 1])),
					    astc::clamp64Kf(sf16_to_float(img1->data16[ze1][ye1][xe1 + 2])),
					    astc::clamp64Kf(sf16_to_float(img1->data16[ze1][ye1][xe1 + 3])));
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
					    astc::clamp64Kf(sf16_to_float(img2->data16[ze2][ye2][xe2])),
					    astc::clamp64Kf(sf16_to_float(img2->data16[ze2][ye2][xe2 + 1])),
					    astc::clamp64Kf(sf16_to_float(img2->data16[ze2][ye2][xe2 + 2])),
					    astc::clamp64Kf(sf16_to_float(img2->data16[ze2][ye2][xe2 + 3])));
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
					    astc::xlog2(color1.x),
					    astc::xlog2(color1.y),
					    astc::xlog2(color1.z),
					    astc::xlog2(color1.w));

					float4 log_input_color2 = float4(
					    astc::xlog2(color2.x),
					    astc::xlog2(color2.y),
					    astc::xlog2(color2.z),
					    astc::xlog2(color2.w));

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

	float pixels = (float)(dim_x * dim_y * dim_z);
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
		psnr = 10.0f * log10f(denom / num);

	float rgb_psnr = psnr;

	printf("Error metrics\n");
	printf("=============\n\n");

	if (channelmask & 8)
	{
		printf("    PSNR (LDR-RGBA):            %9.6f dB\n", (double)psnr);

		float alpha_psnr;
		if (alpha_num == 0.0f)
			alpha_psnr = 999.0f;
		else
			alpha_psnr = 10.0f * log10f(denom / alpha_num);
		printf("    Alpha-weighted PSNR:        %9.6f dB\n", (double)alpha_psnr);

		float rgb_num = errorsum.sum.x + errorsum.sum.y + errorsum.sum.z;
		if (rgb_num == 0.0f)
			rgb_psnr = 999.0f;
		else
			rgb_psnr = 10.0f * log10f(pixels * 3.0f / rgb_num);
		printf("    PSNR (LDR-RGB):             %9.6f dB\n", (double)rgb_psnr);
	}
	else
	{
		printf("    PSNR (LDR-RGB):             %9.6f dB\n", (double)psnr);
	}

	if (compute_hdr_metrics)
	{
		printf("    PSNR (RGB norm to peak):    %9.6f dB (peak %f)\n",
				(double)(rgb_psnr + 20.0f * log10f(rgb_peak)),
				(double)rgb_peak);

		float mpsnr;
		if (mpsnr_num == 0.0f)
			mpsnr = 999.0f;
		else
			mpsnr = 10.0f * log10f(mpsnr_denom / mpsnr_num);
		printf("    mPSNR (RGB):                %9.6f dB (fstops %+d to %+d)\n",
				(double)mpsnr, fstop_lo, fstop_hi);

		float logrmse = astc::sqrt(log_num / pixels);
		printf("    LogRMSE (RGB):              %9.6f\n", (double)logrmse);
	}

	printf("\n");
}
