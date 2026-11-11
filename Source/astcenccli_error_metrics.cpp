// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2011-2026 Arm Limited
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

#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>

#include "astcenccli_internal.h"

/**
 * @brief An accumulator for errors.
 */
class error_accum4
{
public:
	/** @brief The running sum. */
	double sum_r { 0.0 };
	double sum_g { 0.0 };
	double sum_b { 0.0 };
	double sum_a { 0.0 };
};

/**
 * @brief SSIM statistics for display in the CLI.
 */
struct ssim_stats
{
	/** @brief Mean RGB SSIM, averaged over exposures for HDR images. */
	double rgb { 0.0 };

	/** @brief Mean SSIM of the selected LDR components, including alpha. */
	double rgba { 0.0 };
};

/**
 * @brief Accumulate positive log luminance for an HDR exposure estimate.
 *
 * @param         color               The linear RGBA color.
 * @param[in,out] log_luminance_sum   The accumulated log2 luminance.
 * @param[in,out] luminance_count     The number of positive luminance samples.
 */
static void accumulate_hdr_luminance(
	vfloat4 color,
	double& log_luminance_sum,
	size_t& luminance_count
) {
	// Ignore alpha and black texels when accumulating the geometric mean
	vfloat4 luminance_weights(0.2126f, 0.7152f, 0.0722f, 0.0f);
	float luminance = dot3_s(color, luminance_weights);
	if (luminance > 0.0f)
	{
		log_luminance_sum += std::log2(static_cast<double>(luminance));
		luminance_count++;
	}
}

/**
 * @brief Select an HDR metric exposure using the original image's geometric mean luminance.
 *
 * @param img   The original image, including any crop view.
 *
 * @return The nearest whole fstop mapping the geometric mean to 18% grey, clamped to [-123, 123],
 *         or zero for black images.
 */
static int compute_hdr_exposure(
	const astcenc_image* img
) {
	// Average log luminance in a fixed order; exclude black texels and alpha
	double log_luminance_sum = 0.0;
	size_t luminance_count = 0;

	// Select the source format outside the loops and respect cropped row strides
	if (img->data_type == ASTCENC_TYPE_U8)
	{
		size_t stride = img->row_stride;
		if (stride == 0)
		{
			stride = img->dim_x * 4;
		}

		for (size_t z = 0; z < img->dim_z; z++)
		{
			const uint8_t* data = static_cast<const uint8_t*>(img->data[z]);
			for (size_t y = 0; y < img->dim_y; y++)
			{
				for (size_t x = 0; x < img->dim_x; x++)
				{
					const uint8_t* texel = data + y * stride + 4 * x;
					vfloat4 color = vfloat4(texel[0], texel[1], texel[2], texel[3]) / 255.0f;
					accumulate_hdr_luminance(color, log_luminance_sum, luminance_count);
				}
			}
		}
	}
	else if (img->data_type == ASTCENC_TYPE_F16)
	{
		size_t stride = img->row_stride / sizeof(uint16_t);
		if (stride == 0)
		{
			stride = img->dim_x * 4;
		}

		for (size_t z = 0; z < img->dim_z; z++)
		{
			const uint16_t* data = static_cast<const uint16_t*>(img->data[z]);
			for (size_t y = 0; y < img->dim_y; y++)
			{
				for (size_t x = 0; x < img->dim_x; x++)
				{
					const uint16_t* texel = data + y * stride + 4 * x;
					vfloat4 color = float16_to_float(vint4(texel[0], texel[1], texel[2], texel[3]));
					color = clamp(0.0f, 65504.0f, color);
					accumulate_hdr_luminance(color, log_luminance_sum, luminance_count);
				}
			}
		}
	}
	else
	{
		assert(img->data_type == ASTCENC_TYPE_F32);
		size_t stride = img->row_stride / sizeof(float);
		if (stride == 0)
		{
			stride = img->dim_x * 4;
		}

		for (size_t z = 0; z < img->dim_z; z++)
		{
			const float* data = static_cast<const float*>(img->data[z]);
			for (size_t y = 0; y < img->dim_y; y++)
			{
				for (size_t x = 0; x < img->dim_x; x++)
				{
					const float* texel = data + y * stride + 4 * x;
					vfloat4 color = clamp(0.0f, 65504.0f, vfloat4(texel));
					accumulate_hdr_luminance(color, log_luminance_sum, luminance_count);
				}
			}
		}
	}

	// Keep a neutral exposure for black images and room for the surrounding two stops
	if (luminance_count == 0)
	{
		return 0;
	}

	double fstop = std::round(std::log2(0.18) - log_luminance_sum / luminance_count);
	return static_cast<int>(astc::clamp(fstop, -123.0, 123.0));
}

/**
 * @brief Incremental addition operator for error accumulators.
 *
 * @param[in,out] val   The accumulator to increment.
 * @param         inc   The increment to apply.
 *
 * @return The updated accumulator.
 */
static error_accum4& operator+=(
	error_accum4 &val,
	vfloat4 inc
) {
	val.sum_r += static_cast<double>(inc.lane<0>());
	val.sum_g += static_cast<double>(inc.lane<1>());
	val.sum_b += static_cast<double>(inc.lane<2>());
	val.sum_a += static_cast<double>(inc.lane<3>());
	return val;
}

/**
 * @brief mPSNR tone-mapping operator for HDR images.
 *
 * @param val     The color value to tone map.
 * @param fstop   The exposure fstop; should be in range [-125, 125].
 *
 * @return The mapped color value in [0.0f, 255.0f] range.
 */
static float mpsnr_operator(
	float val,
	int fstop
) {
	// Fast implementation of pow(2.0, fstop), assuming IEEE float layout
	// Use unsigned arithmetic to construct the exponent bit pattern
	unsigned int uscale = 0x3f800000u + (std::bit_cast<unsigned int>(fstop) << 23);
	float scale = std::bit_cast<float>(uscale);

	val = powf(val * scale, (1.0f / 2.2f));
	return astc::clamp(val * 255.0f, 0.0f, 255.0f);
}

/**
 * @brief Sum squared tone-mapped differences between two values.
 *
 * Each difference is squared after applying the mPSNR tone mapping.
 *
 * @param val1       The first color value.
 * @param val2       The second color value.
 * @param fstop_lo   The low exposure fstop; should be in range [-125, 125].
 * @param fstop_hi   The high exposure fstop; should be in range [-125, 125].
 *
 * @return The summed squared difference across all active fstop levels.
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

/**
 * @brief Load and pad one row of a single color channel for the SSIM filter.
 *
 * Samples are shifted by -0.5 before filtering to reduce cancellation in variance calculations.
 *
 * @param      img         The source image.
 * @param      component   The color component index, in RGBA order.
 * @param      z           The depth slice index.
 * @param      y           The row index.
 * @param      width       The comparison width.
 * @param      radius      The filter radius.
 * @param      is_hdr      True to apply the mPSNR tone map.
 * @param      fstop       The exposure fstop for HDR tone mapping.
 * @param[out] row         The padded output row, centered around zero.
 */
static void load_ssim_row(
	const astcenc_image* img,
	unsigned int component,
	size_t z,
	size_t y,
	size_t width,
	size_t radius,
	bool is_hdr,
	int fstop,
	std::vector<float>& row
) {
	// Select the source type and exposure mode outside the pixel loops
	if (img->data_type == ASTCENC_TYPE_U8)
	{
		size_t stride = img->row_stride;
		if (stride == 0)
		{
			stride = img->dim_x * 4;
		}

		// Normalize UNORM samples before optional display tone mapping
		const uint8_t* data = static_cast<const uint8_t*>(img->data[z]) + y * stride + component;
		if (is_hdr)
		{
			for (size_t x = 0; x < width; x++)
			{
				float value = data[4 * x] / 255.0f;
				row[x + radius] = mpsnr_operator(value, fstop) / 255.0f - 0.5f;
			}
		}
		else
		{
			for (size_t x = 0; x < width; x++)
			{
				row[x + radius] = data[4 * x] / 255.0f - 0.5f;
			}
		}
	}
	else if (img->data_type == ASTCENC_TYPE_F16)
	{
		size_t stride = img->row_stride / sizeof(uint16_t);
		if (stride == 0)
		{
			stride = img->dim_x * 4;
		}

		// Convert FP16 samples and clamp unconstrained inputs before optional tone mapping
		const uint16_t* data = static_cast<const uint16_t*>(img->data[z]) + y * stride + component;
		if (is_hdr)
		{
			for (size_t x = 0; x < width; x++)
			{
				vint4 valuei(data[4 * x]);
				float value = float16_to_float(valuei).lane<0>();
				value = astc::clamp(value, 0.0f, 65504.0f);
				row[x + radius] = mpsnr_operator(value, fstop) / 255.0f - 0.5f;
			}
		}
		else
		{
			for (size_t x = 0; x < width; x++)
			{
				vint4 valuei(data[4 * x]);
				float value = float16_to_float(valuei).lane<0>();
				row[x + radius] = astc::clamp(value, 0.0f, 65504.0f) - 0.5f;
			}
		}
	}
	else
	{
		assert(img->data_type == ASTCENC_TYPE_F32);
		size_t stride = img->row_stride / sizeof(float);
		if (stride == 0)
		{
			stride = img->dim_x * 4;
		}

		// Clamp unconstrained FP32 samples before optional display tone mapping
		const float* data = static_cast<const float*>(img->data[z]) + y * stride + component;
		if (is_hdr)
		{
			for (size_t x = 0; x < width; x++)
			{
				float value = astc::clamp(data[4 * x], 0.0f, 65504.0f);
				row[x + radius] = mpsnr_operator(value, fstop) / 255.0f - 0.5f;
			}
		}
		else
		{
			for (size_t x = 0; x < width; x++)
			{
				row[x + radius] = astc::clamp(data[4 * x], 0.0f, 65504.0f) - 0.5f;
			}
		}
	}

	// Replicate edge values for the filter window and padded SIMD lanes
	for (size_t x = 0; x < radius; x++)
	{
		row[x] = row[radius];
	}

	for (size_t x = width + radius; x < row.size(); x++)
	{
		row[x] = row[width + radius - 1];
	}
}

/**
 * @brief Compute single-channel SSIM using an 11x11 Gaussian window with sigma 1.5.
 *
 * Uses population covariance, K1 = 0.01, K2 = 0.03, and a unit data range. Image edges are extended
 * by replication, including for images smaller than the window. All pixels in the image
 * intersection contribute equally, with each depth slice filtered independently.
 *
 * @param img1        The original image.
 * @param img2        The compressed image.
 * @param component   The color component index, in RGBA order.
 * @param is_hdr      True to apply the mPSNR tone map before computing SSIM.
 * @param fstop       The exposure fstop for HDR tone mapping, in range [-125, 125].
 *
 * @return The mean SSIM, where 1 indicates identical images.
 */
static float compute_ssim(
	const astcenc_image* img1,
	const astcenc_image* img2,
	unsigned int component,
	bool is_hdr,
	int fstop
) {
	// Compare the image intersection, filtering each depth slice independently
	assert(component < 4);
	size_t width = astc::min(img1->dim_x, img2->dim_x);
	size_t height = astc::min(img1->dim_y, img2->dim_y);
	size_t depth = astc::min(img1->dim_z, img2->dim_z);
	assert((width > 0) && (height > 0) && (depth > 0));

	// Pad rows so each vector lane filters one pixel, including the partial final vector
	constexpr size_t radius { 5 };
	constexpr size_t window { 2 * radius + 1 };
	size_t padded_width = round_up_to_simd_multiple_vla(width);

	// Build a normalized Gaussian kernel for both passes of the separable filter
	std::array<float, window> weights;

	float weight_sum = 0.0f;
	for (size_t k = 0; k < window; k++)
	{
		float offset = static_cast<float>(k) - static_cast<float>(radius);
		weights[k] = std::exp(-(offset * offset) / (2.0f * 1.5f * 1.5f));
		weight_sum += weights[k];
	}

	for (float& weight : weights)
	{
		weight /= weight_sum;
	}

	// Cache horizontally filtered means, sum of second moments, and squared differences
	// A ring of window rows keeps storage under control
	std::array<std::vector<float>, 4> rows;
	for (auto& row : rows)
	{
		row.resize(padded_width * window);
	}

	// Reuse padded source rows and pixel losses across output rows and depth slices
	std::vector<float> row1(padded_width + 2 * radius);
	std::vector<float> row2(row1.size());
	std::vector<float> row_loss(padded_width);
	vfloat4 loss_sum = vfloat4::zero();
	for (size_t z = 0; z < depth; z++)
	{
		size_t next_y = 0;
		size_t next_row = 0;
		for (size_t y = 0; y < height; y++)
		{
			// Cache horizontally filtered rows until the vertical window is available
			size_t last_y = astc::min(y + radius, height - 1);
			for (; next_y <= last_y; next_y++)
			{
				load_ssim_row(img1, component, z, next_y, width, radius, is_hdr, fstop, row1);
				load_ssim_row(img2, component, z, next_y, width, radius, is_hdr, fstop, row2);

				// Filter the four moments horizontally, with one texel per SIMD lane
				size_t row_offset = next_row * padded_width;
				for (size_t x = 0; x < padded_width; x += ASTCENC_SIMD_WIDTH)
				{
					vfloat mean1 = vfloat::zero();
					vfloat mean2 = vfloat::zero();
					vfloat squares = vfloat::zero();
					vfloat differences = vfloat::zero();
					for (size_t k = 0; k < window; k++)
					{
						vfloat a(row1.data() + x + k);
						vfloat b(row2.data() + x + k);
						vfloat weight(weights[k]);
						vfloat diff = a - b;
						mean1 += weight * a;
						mean2 += weight * b;
						squares += weight * (a * a + b * b);
						differences += weight * (diff * diff);
					}

					store(mean1, rows[0].data() + row_offset + x);
					store(mean2, rows[1].data() + row_offset + x);
					store(squares, rows[2].data() + row_offset + x);
					store(differences, rows[3].data() + row_offset + x);
				}

				next_row++;
				if (next_row == window)
				{
					next_row = 0;
				}
			}

			// Resolve the ring slots once per row, including replicated top and bottom edges
			std::array<size_t, window> row_offsets;
			for (size_t k = 0; k < window; k++)
			{
				size_t sy = astc::min(y + k, height - 1 + radius);
				sy = astc::max(sy, radius) - radius;
				size_t slot = next_row + window - 1 - (last_y - sy);
				if (slot >= window)
				{
					slot -= window;
				}

				row_offsets[k] = slot * padded_width;
			}

			// Apply the vertical Gaussian pass to recover local moments for each pixel
			for (size_t x = 0; x < padded_width; x += ASTCENC_SIMD_WIDTH)
			{
				vfloat mean1 = vfloat::zero();
				vfloat mean2 = vfloat::zero();
				vfloat squares = vfloat::zero();
				vfloat differences = vfloat::zero();

				for (size_t k = 0; k < window; k++)
				{
					size_t row_offset = row_offsets[k] + x;
					vfloat weight(weights[k]);
					mean1 += weight * vfloat(rows[0].data() + row_offset);
					mean2 += weight * vfloat(rows[1].data() + row_offset);
					squares += weight * vfloat(rows[2].data() + row_offset);
					differences += weight * vfloat(rows[3].data() + row_offset);
				}

				// Express SSIM as a loss to avoid cancellation for nearly identical images
				vfloat diff = mean1 - mean2;
				vfloat mean_diff_sq = diff * diff;
				vfloat variance = max(squares - (mean1 * mean1 + mean2 * mean2), vfloat::zero());
				vfloat diff_variance = max(differences - mean_diff_sq, vfloat::zero());

				// Restore uncentered means for luminance; variances are unchanged by centering
				mean1 += vfloat(0.5f);
				mean2 += vfloat(0.5f);
				vfloat luminance_loss = mean_diff_sq / (mean1 * mean1 + mean2 * mean2 + vfloat(0.0001f));
				vfloat structure_loss = diff_variance / (variance + vfloat(0.0009f));
				vfloat loss = luminance_loss + structure_loss - luminance_loss * structure_loss;
				store(loss, row_loss.data() + x);
			}

			// Reduce rows in fixed-width groups to limit drift and preserve invariance
			for (size_t x = width; x < padded_width; x++)
			{
				row_loss[x] = 0.0f;
			}

			vfloat4 row_sum = vfloat4::zero();
			for (size_t x = 0; x < width; x += 4)
			{
				row_sum += vfloat4(row_loss.data() + x);
			}

			loss_sum += row_sum;
		}
	}

	// Convert mean loss over all compared texels back to an SSIM score
	return 1.0f - hadd_s(loss_sum) / static_cast<float>(width * height * depth);
}

/**
 * @brief Compute SSIM statistics for reporting.
 *
 * @param img1             The original image.
 * @param img2             The compressed image.
 * @param component_mask   The mask of input color components.
 * @param is_hdr           True to average tone-mapped RGB scores over HDR exposures.
 * @param fstop_lo         The low exposure fstop (HDR only).
 * @param fstop_hi         The high exposure fstop (HDR only).
 *
 * @return The aggregated RGB and LDR RGBA SSIM statistics.
 */
static ssim_stats compute_ssim_stats(
	const astcenc_image* img1,
	const astcenc_image* img2,
	int component_mask,
	bool is_hdr,
	int fstop_lo,
	int fstop_hi
) {
	// HDR SSIM averages the same display exposures as mPSNR; alpha is only scored for LDR
	int ssim_lo = 0;
	int ssim_hi = 0;
	unsigned int component_count = 3;

	if (is_hdr)
	{
		ssim_lo = fstop_lo;
		ssim_hi = fstop_hi;
	}
	else if (component_mask & 8)
	{
		component_count = 4;
	}

	// Average each channel over exposures; LDR uses a single pass without tone mapping
	double channel_ssim[4] {};
	for (unsigned int c = 0; c < component_count; c++)
	{
		for (int fstop = ssim_lo; fstop <= ssim_hi; fstop++)
		{
			channel_ssim[c] += static_cast<double>(compute_ssim(img1, img2, c, is_hdr, fstop));
		}

		channel_ssim[c] /= static_cast<double>(ssim_hi - ssim_lo + 1);
	}

	// Match PSNR's component mask, including luminance-alpha inputs
	double rgba_ssim = 0.0;
	if ((!is_hdr) && (component_mask & 8))
	{
		unsigned int count = 0;
		for (unsigned int c = 0; c < 4; c++)
		{
			if (component_mask & (1 << c))
			{
				rgba_ssim += channel_ssim[c];
				count++;
			}
		}

		rgba_ssim /= count;
	}

	// RGB always averages three channels, including replicated luminance for greyscale inputs
	return ssim_stats {
		.rgb = (channel_ssim[0] + channel_ssim[1] + channel_ssim[2]) / 3.0,
		.rgba = rgba_ssim
	};
}

/* See header for documentation */
void compute_error_metrics(
	bool is_hdr,
	bool is_hdr_fstop_auto,
	bool is_normal_map,
	int input_components,
	const astcenc_image* img1,
	const astcenc_image* img2,
	int fstop_lo,
	int fstop_hi
) {
	// Use one reference-derived exposure range for both metrics and both images
	if (is_hdr && is_hdr_fstop_auto)
	{
		int fstop = compute_hdr_exposure(img1);
		fstop_lo = fstop - 2;
		fstop_hi = fstop + 2;
	}

	static const int componentmasks[5] { 0x00, 0x07, 0x0C, 0x07, 0x0F };
	int component_mask = componentmasks[input_components];

	error_accum4 errorsum;
	error_accum4 alpha_scaled_errorsum;
	error_accum4 log_errorsum;
	error_accum4 mpsnr_errorsum;
	double mean_angular_errorsum = 0.0;
	double worst_angular_errorsum = 0.0;

	size_t dim_x = astc::min(img1->dim_x, img2->dim_x);
	size_t dim_y = astc::min(img1->dim_y, img2->dim_y);
	size_t dim_z = astc::min(img1->dim_z, img2->dim_z);

	if (img1->dim_x != img2->dim_x ||
	    img1->dim_y != img2->dim_y ||
	    img1->dim_z != img2->dim_z)
	{
		printf("WARNING: Only intersection of images will be compared:\n"
		       "  Image 1: %zux%zux%zu\n"
		       "  Image 2: %zux%zux%zu\n",
		       img1->dim_x, img1->dim_y, img1->dim_z,
		       img2->dim_x, img2->dim_y, img2->dim_z);
	}

	double rgb_peak = 0.0;
	size_t row_stride1 = img1->row_stride;
	if (img1->data_type == ASTCENC_TYPE_F16)
	{
		row_stride1 /= sizeof(uint16_t);
	}
	else if (img1->data_type == ASTCENC_TYPE_F32)
	{
		row_stride1 /= sizeof(float);
	}

	size_t row_stride2 = img2->row_stride;
	if (img2->data_type == ASTCENC_TYPE_F16)
	{
		row_stride2 /= sizeof(uint16_t);
	}
	else if (img2->data_type == ASTCENC_TYPE_F32)
	{
		row_stride2 /= sizeof(float);
	}

	for (size_t z = 0; z < dim_z; z++)
	{
		for (size_t y = 0; y < dim_y; y++)
		{
			for (size_t x = 0; x < dim_x; x++)
			{
				vfloat4 color1;
				vfloat4 color2;

				if (img1->data_type == ASTCENC_TYPE_U8)
				{
					uint8_t* data8 = static_cast<uint8_t*>(img1->data[z]);

					color1 = vfloat4(
					    data8[(row_stride1 * y) + (4 * x    )],
					    data8[(row_stride1 * y) + (4 * x + 1)],
					    data8[(row_stride1 * y) + (4 * x + 2)],
					    data8[(row_stride1 * y) + (4 * x + 3)]);

					color1 = color1 / 255.0f;
				}
				else if (img1->data_type == ASTCENC_TYPE_F16)
				{
					uint16_t* data16 = static_cast<uint16_t*>(img1->data[z]);

					vint4 color1i = vint4(
					    data16[(row_stride1 * y) + (4 * x    )],
					    data16[(row_stride1 * y) + (4 * x + 1)],
					    data16[(row_stride1 * y) + (4 * x + 2)],
					    data16[(row_stride1 * y) + (4 * x + 3)]);

					color1 = float16_to_float(color1i);
					color1 = clamp(0, 65504.0f, color1);
				}
				else  // if (img1->data_type == ASTCENC_TYPE_F32)
				{
					assert(img1->data_type == ASTCENC_TYPE_F32);
					float* data32 = static_cast<float*>(img1->data[z]);

					color1 = vfloat4(
					    data32[(row_stride1 * y) + (4 * x    )],
					    data32[(row_stride1 * y) + (4 * x + 1)],
					    data32[(row_stride1 * y) + (4 * x + 2)],
					    data32[(row_stride1 * y) + (4 * x + 3)]);

					color1 = clamp(0, 65504.0f, color1);
				}

				if (img2->data_type == ASTCENC_TYPE_U8)
				{
					uint8_t* data8 = static_cast<uint8_t*>(img2->data[z]);

					color2 = vfloat4(
					    data8[(row_stride2 * y) + (4 * x    )],
					    data8[(row_stride2 * y) + (4 * x + 1)],
					    data8[(row_stride2 * y) + (4 * x + 2)],
					    data8[(row_stride2 * y) + (4 * x + 3)]);

					color2 = color2 / 255.0f;
				}
				else if (img2->data_type == ASTCENC_TYPE_F16)
				{
					uint16_t* data16 = static_cast<uint16_t*>(img2->data[z]);

					vint4 color2i = vint4(
					    data16[(row_stride2 * y) + (4 * x    )],
					    data16[(row_stride2 * y) + (4 * x + 1)],
					    data16[(row_stride2 * y) + (4 * x + 2)],
					    data16[(row_stride2 * y) + (4 * x + 3)]);

					color2 = float16_to_float(color2i);
					color2 = clamp(0, 65504.0f, color2);
				}
				else  // if (img2->data_type == ASTCENC_TYPE_F32)
				{
					assert(img2->data_type == ASTCENC_TYPE_F32);
					float* data32 = static_cast<float*>(img2->data[z]);

					color2 = vfloat4(
					    data32[(row_stride2 * y) + (4 * x    )],
					    data32[(row_stride2 * y) + (4 * x + 1)],
					    data32[(row_stride2 * y) + (4 * x + 2)],
					    data32[(row_stride2 * y) + (4 * x + 3)]);

					color2 = clamp(0, 65504.0f, color2);
				}

				rgb_peak = astc::max(static_cast<double>(color1.lane<0>()),
				                     static_cast<double>(color1.lane<1>()),
				                     static_cast<double>(color1.lane<2>()),
				                     rgb_peak);

				vfloat4 diffcolor = color1 - color2;
				vfloat4 diffcolor_sq = diffcolor * diffcolor;
				errorsum += diffcolor_sq;

				vfloat4 alpha_scaled_diffcolor = vfloat4(
				    diffcolor.lane<0>() * color1.lane<3>(),
				    diffcolor.lane<1>() * color1.lane<3>(),
				    diffcolor.lane<2>() * color1.lane<3>(),
				    diffcolor.lane<3>());

				vfloat4 alpha_scaled_diffcolor_sq = alpha_scaled_diffcolor * alpha_scaled_diffcolor;
				alpha_scaled_errorsum += alpha_scaled_diffcolor_sq;

				if (is_hdr)
				{
					vfloat4 log_input_color1 = log2(color1);
					vfloat4 log_input_color2 = log2(color2);

					vfloat4 log_diffcolor = log_input_color1 - log_input_color2;

					log_errorsum += log_diffcolor * log_diffcolor;

					vfloat4 mpsnr_error = vfloat4(
					    mpsnr_sumdiff(color1.lane<0>(), color2.lane<0>(), fstop_lo, fstop_hi),
					    mpsnr_sumdiff(color1.lane<1>(), color2.lane<1>(), fstop_lo, fstop_hi),
					    mpsnr_sumdiff(color1.lane<2>(), color2.lane<2>(), fstop_lo, fstop_hi),
					    mpsnr_sumdiff(color1.lane<3>(), color2.lane<3>(), fstop_lo, fstop_hi));

					mpsnr_errorsum += mpsnr_error;
				}

				if (is_normal_map)
				{
					// Decode the normal vector
					vfloat4 normal1 = (color1 - 0.5f) * 2.0f;
					normal1 = normalize_safe(normal1.swz<0, 1, 2>(), unit3());

					vfloat4 normal2 = (color2 - 0.5f) * 2.0f;
					normal2 = normalize_safe(normal2.swz<0, 1, 2>(), unit3());

					// Float error can push this outside of valid range for acos, so clamp to avoid NaN issues
					float normal_cos = clamp(-1.0f, 1.0f, dot3(normal1, normal2)).lane<0>();
					float rad_to_degrees = 180.0f / astc::PI;
					double error_degrees = std::acos(static_cast<double>(normal_cos)) * static_cast<double>(rad_to_degrees);

					mean_angular_errorsum += error_degrees / (dim_x * dim_y * dim_z);
					worst_angular_errorsum = astc::max(worst_angular_errorsum, error_degrees);
				}
			}
		}
	}

	double pixels = static_cast<double>(dim_x * dim_y * dim_z);
	double samples = 0.0;

	double num = 0.0;
	double alpha_num = 0.0;
	double log_num = 0.0;
	double mpsnr_num = 0.0;

	if (component_mask & 1)
	{
		num += errorsum.sum_r;
		alpha_num += alpha_scaled_errorsum.sum_r;
		log_num += log_errorsum.sum_r;
		mpsnr_num += mpsnr_errorsum.sum_r;
		samples += pixels;
	}

	if (component_mask & 2)
	{
		num += errorsum.sum_g;
		alpha_num += alpha_scaled_errorsum.sum_g;
		log_num += log_errorsum.sum_g;
		mpsnr_num += mpsnr_errorsum.sum_g;
		samples += pixels;
	}

	if (component_mask & 4)
	{
		num += errorsum.sum_b;
		alpha_num += alpha_scaled_errorsum.sum_b;
		log_num += log_errorsum.sum_b;
		mpsnr_num += mpsnr_errorsum.sum_b;
		samples += pixels;
	}

	if (component_mask & 8)
	{
		num += errorsum.sum_a;
		alpha_num += alpha_scaled_errorsum.sum_a;
		samples += pixels;
	}

	double denom = samples;
	double stopcount = static_cast<double>(fstop_hi - fstop_lo + 1);
	double mpsnr_denom = pixels * 3.0 * stopcount * 255.0 * 255.0;

	double psnr;
	if (num == 0.0)
	{
		psnr = 999.0;
	}
	else
	{
		psnr = 10.0 * log10(denom / num);
	}

	double rgb_psnr = psnr;
	ssim_stats ssim = compute_ssim_stats(img1, img2, component_mask, is_hdr, fstop_lo, fstop_hi);

	printf("Quality metrics\n");
	printf("===============\n\n");

	if (component_mask & 8)
	{
		printf("    PSNR (LDR-RGBA):          %9.4f dB\n", psnr);

		double alpha_psnr;
		if (alpha_num == 0.0)
		{
			alpha_psnr = 999.0;
		}
		else
		{
			alpha_psnr = 10.0 * log10(denom / alpha_num);
		}
		printf("    Alpha-weighted PSNR:      %9.4f dB\n", alpha_psnr);

		double rgb_num = errorsum.sum_r + errorsum.sum_g + errorsum.sum_b;
		if (rgb_num == 0.0)
		{
			rgb_psnr = 999.0;
		}
		else
		{
			rgb_psnr = 10.0 * log10(pixels * 3.0 / rgb_num);
		}
		printf("    PSNR (LDR-RGB):           %9.4f dB\n", rgb_psnr);
	}
	else
	{
		printf("    PSNR (LDR-RGB):           %9.4f dB\n", psnr);
	}

	if (is_hdr)
	{
		printf("    PSNR (RGB norm to peak):  %9.4f dB (peak %f)\n",
		       rgb_psnr + 20.0 * log10(rgb_peak), rgb_peak);

		double mpsnr;
		if (mpsnr_num == 0.0)
		{
			mpsnr = 999.0;
		}
		else
		{
			mpsnr = 10.0 * log10(mpsnr_denom / mpsnr_num);
		}

		printf("    mPSNR (RGB):              %9.4f dB (fstops %+d to %+d)\n",
		       mpsnr, fstop_lo, fstop_hi);

		double logrmse = sqrt(log_num / pixels);
		printf("    LogRMSE (RGB):            %9.4f\n", logrmse);
	}

	if (is_hdr)
	{
		printf("    mSSIM (RGB):              %9.4f (fstops %+d to %+d)\n",
		       ssim.rgb, fstop_lo, fstop_hi);
	}
	else
	{
		if (component_mask & 8)
		{
			printf("    SSIM (LDR-RGBA):          %9.4f\n", ssim.rgba);
		}

		printf("    SSIM (LDR-RGB):           %9.4f\n", ssim.rgb);
	}

	if (is_normal_map)
	{
		printf("    Mean Angular Error:       %9.4f degrees\n", mean_angular_errorsum);
		printf("    Worst Angular Error:      %9.4f degrees\n", worst_angular_errorsum);
	}

	printf("\n");
}
