// ----------------------------------------------------------------------------
//  This confidential and proprietary software may be used only as authorised
//  by a licensing agreement from Arm Limited.
//      (C) COPYRIGHT 2011-2020 Arm Limited, ALL RIGHTS RESERVED
//  The entire notice above must be reproduced on all authorised copies and
//  copies may only be made to the extent permitted by a licensing agreement
//  from Arm Limited.
// ----------------------------------------------------------------------------

/**
 * @brief Functions to calculate variance per channel in a NxN footprint.
 *
 * We need N to be parametric, so the routine below uses summed area tables in
 * order to execute in O(1) time independent of how big N is.
 *
 * The addition uses a Brent-Kung-based parallel prefix adder. This uses the
 * prefix tree to first perform a binary reduction, and then distributes the
 * results. This method means that there is no serial dependency between a
 * given element and the next one, and also significantly improves numerical
 * stability allowing us to use floats rather than doubles.
 */

#include "astc_codec_internals.h"
#include <cassert>

/**
 * @brief Parameter structure for compute_pixel_region_variance().
 *
 * This function takes a structure to avoid spilling arguments to the stack
 * on every function invocation, as there are a lot of parameters.
 */
struct pixel_region_variance_args
{
	/** The image to analyze. */
	const astc_codec_image *img;
	/** The RGB channel power adjustment. */
	float rgb_power;
	/** The alpha channel power adjustment. */
	float alpha_power;
	/** The RGB data should be treated as sRGB. */
	int need_srgb_transform;
	/** The channel swizzle pattern. */
	swizzlepattern swz;
	/** Should the algorithm bother with Z axis processing? */
	int have_z;
	/** The kernel radius for average and variance. */
	int avg_var_kernel_radius;
	/** The kernel radius for alpha processing. */
	int alpha_kernel_radius;
	/** The size of the working data to process. */
	int3 size;
	/** The position of first src data in the data set. */
	int3 src_offset;
	/** The position of first dst data in the data set. */
	int3 dst_offset;
	/** The working memory buffer. */
	float4 *work_memory;
};

/**
 * @brief Parameter structure for compute_averages_and_variances_proc().
 */
struct avg_var_args
{
	/** The arguments for the nested variance computation. */
	pixel_region_variance_args arg;
	/** The image dimensions. */
	int3 img_size;
	/** The maximum working block dimensions. */
	int3 blk_size;
	/** The working block memory size. */
	int work_memory_size;
};

/**
 * @brief Generate a prefix-sum array using Brent-Kung algorithm.
 *
 * This will take an input array of the form:
 *     v0, v1, v2, ...
 * ... and modify in-place to turn it into a prefix-sum array of the form:
 *     v0, v0+v1, v0+v1+v2, ...
 *
 * @param d      The array to prefix-sum.
 * @param items  The number of items in the array.
 * @param stride The item spacing in the array; i.e. dense arrays should use 1.
 */
static void brent_kung_prefix_sum(
	float4* d,
	size_t items,
	int stride
) {
	if (items < 2)
		return;

	size_t lc_stride = 2;
	size_t log2_stride = 1;

	// The reduction-tree loop
	do {
		size_t step = lc_stride >> 1;
		size_t start = lc_stride - 1;
		size_t iters = items >> log2_stride;

		float4 *da = d + (start * stride);
		ptrdiff_t ofs = -(step * stride);
		size_t ofs_stride = stride << log2_stride;

		while (iters)
		{
			*da = *da + da[ofs];
			da += ofs_stride;
			iters--;
		}

		log2_stride += 1;
		lc_stride <<= 1;
	} while (lc_stride <= items);

	// The expansion-tree loop
	do {
		log2_stride -= 1;
		lc_stride >>= 1;

		size_t step = lc_stride >> 1;
		size_t start = step + lc_stride - 1;
		size_t iters = (items - step) >> log2_stride;

		float4 *da = d + (start * stride);
		ptrdiff_t ofs = -(step * stride);
		size_t ofs_stride = stride << log2_stride;

		while (iters)
		{
			*da = *da + da[ofs];
			da += ofs_stride;
			iters--;
		}
	} while (lc_stride > 2);
}

/**
 * @brief Compute averages and variances for a pixel region.
 *
 * The routine computes both in a single pass, using a summed-area table to
 * decouple the running time from the averaging/variance kernel size.
 *
 * @param arg The input parameter structure.
 */
static void compute_pixel_region_variance(
	const pixel_region_variance_args* arg
) {
	// Unpack the memory structure into local variables
	const astc_codec_image *img = arg->img;
	float rgb_power = arg->rgb_power;
	float alpha_power = arg->alpha_power;
	int need_srgb_transform = arg->need_srgb_transform;
	swizzlepattern swz = arg->swz;
	int have_z = arg->have_z;

	int size_x = arg->size.x;
	int size_y = arg->size.y;
	int size_z = arg->size.z;

	int src_offset_x = arg->src_offset.x;
	int src_offset_y = arg->src_offset.y;
	int src_offset_z = arg->src_offset.z;

	int dst_offset_x = arg->dst_offset.x;
	int dst_offset_y = arg->dst_offset.y;
	int dst_offset_z = arg->dst_offset.z;

	int avg_var_kernel_radius = arg->avg_var_kernel_radius;
	int alpha_kernel_radius = arg->alpha_kernel_radius;

	float  *input_alpha_averages = img->input_alpha_averages;
	float4 *input_averages = img->input_averages;
	float4 *input_variances = img->input_variances;
	float4 *work_memory = arg->work_memory;

	// Compute memory sizes and dimensions that we need
	int kernel_radius = MAX(avg_var_kernel_radius, alpha_kernel_radius);
	int kerneldim = 2 * kernel_radius + 1;
	int kernel_radius_xy = kernel_radius;
	int kernel_radius_z = have_z ? kernel_radius : 0;

	int padsize_x = size_x + kerneldim;
	int padsize_y = size_y + kerneldim;
	int padsize_z = size_z + (have_z ? kerneldim : 0);
	int sizeprod = padsize_x * padsize_y * padsize_z;

	int zd_start = have_z ? 1 : 0;
	int are_powers_1 = (rgb_power == 1.0f) && (alpha_power == 1.0f);

	float4 *varbuf1 = work_memory;
	float4 *varbuf2 = work_memory + sizeprod;

	// Scaling factors to apply to Y and Z for accesses into the work buffers
	int yst = padsize_x;
	int zst = padsize_x * padsize_y;

	// Scaling factors to apply to Y and Z for accesses into result buffers
	int ydt = img->xsize;
	int zdt = img->xsize * img->ysize;

	// Macros to act as accessor functions for the work-memory
	#define VARBUF1(z, y, x) varbuf1[z * zst + y * yst + x]
	#define VARBUF2(z, y, x) varbuf2[z * zst + y * yst + x]

	// Load N and N^2 values into the work buffers
	if (img->data8)
	{
		// Swizzle data structure 4 = ZERO, 5 = ONE
		uint8_t data[6];
		data[4] = 0;
		data[5] = 255;

		for (int z = zd_start; z < padsize_z; z++)
		{
			int z_src = (z - zd_start) + src_offset_z - kernel_radius_z;
			for (int y = 1; y < padsize_y; y++)
			{
				int y_src = (y - 1) + src_offset_y - kernel_radius_xy;
				for (int x = 1; x < padsize_x; x++)
				{
					int x_src = (x - 1) + src_offset_x - kernel_radius_xy;
					data[0] = img->data8[z_src][y_src][4 * x_src + 0];
					data[1] = img->data8[z_src][y_src][4 * x_src + 1];
					data[2] = img->data8[z_src][y_src][4 * x_src + 2];
					data[3] = img->data8[z_src][y_src][4 * x_src + 3];

					uint8_t r = data[swz.r];
					uint8_t g = data[swz.g];
					uint8_t b = data[swz.b];
					uint8_t a = data[swz.a];

					float4 d = float4 (r * (1.0f / 255.0f),
					                   g * (1.0f / 255.0f),
					                   b * (1.0f / 255.0f),
					                   a * (1.0f / 255.0f));

					if (need_srgb_transform)
					{
						d.x = srgb_transform(d.x);
						d.y = srgb_transform(d.y);
						d.z = srgb_transform(d.z);
					}

					if (!are_powers_1)
					{
						d.x = powf(MAX(d.x, 1e-6f), rgb_power);
						d.y = powf(MAX(d.y, 1e-6f), rgb_power);
						d.z = powf(MAX(d.z, 1e-6f), rgb_power);
						d.w = powf(MAX(d.w, 1e-6f), alpha_power);
					}

					VARBUF1(z, y, x) = d;
					VARBUF2(z, y, x) = d * d;
				}
			}
		}
	}
	else
	{
		// Swizzle data structure 4 = ZERO, 5 = ONE (in FP16)
		uint16_t data[6];
		data[4] = 0;
		data[5] = 0x3C00;

		for (int z = zd_start; z < padsize_z; z++)
		{
			int z_src = (z - zd_start) + src_offset_z - kernel_radius_z;

			for (int y = 1; y < padsize_y; y++)
			{
				int y_src = (y - 1) + src_offset_y - kernel_radius_xy;
				for (int x = 1; x < padsize_x; x++)
				{
					int x_src = (x - 1) + src_offset_x - kernel_radius_xy;
					data[0] = img->data16[z_src][y_src][4 * x_src];
					data[1] = img->data16[z_src][y_src][4 * x_src + 1];
					data[2] = img->data16[z_src][y_src][4 * x_src + 2];
					data[3] = img->data16[z_src][y_src][4 * x_src + 3];

					uint16_t r = data[swz.r];
					uint16_t g = data[swz.g];
					uint16_t b = data[swz.b];
					uint16_t a = data[swz.a];

					float4 d = float4(sf16_to_float(r),
					                  sf16_to_float(g),
					                  sf16_to_float(b),
					                  sf16_to_float(a));

					if (need_srgb_transform)
					{
						d.x = srgb_transform(d.x);
						d.y = srgb_transform(d.y);
						d.z = srgb_transform(d.z);
					}

					if (!are_powers_1)
					{
						d.x = powf(MAX(d.x, 1e-6f), rgb_power);
						d.y = powf(MAX(d.y, 1e-6f), rgb_power);
						d.z = powf(MAX(d.z, 1e-6f), rgb_power);
						d.w = powf(MAX(d.w, 1e-6f), alpha_power);
					}

					VARBUF1(z, y, x) = d;
					VARBUF2(z, y, x) = d * d;
				}
			}
		}
	}

	// Pad with an extra layer of 0s; this forms the edge of the SAT tables
	float4 vbz = float4(0, 0, 0, 0);
	for (int z = 0; z < padsize_z; z++)
	{
		for (int y = 0; y < padsize_y; y++)
		{
			VARBUF1(z, y, 0) = vbz;
			VARBUF2(z, y, 0) = vbz;
		}

		for (int x = 0; x < padsize_x; x++)
		{
			VARBUF1(z, 0, x) = vbz;
			VARBUF2(z, 0, x) = vbz;
		}
	}

	if (have_z)
	{
		for (int y = 0; y < padsize_y; y++)
		{
			for (int x = 0; x < padsize_x; x++)
			{
				VARBUF1(0, y, x) = vbz;
				VARBUF2(0, y, x) = vbz;
			}
		}
	}

	// Generate summed-area tables for N and N^2; this is done in-place, using
	// a Brent-Kung parallel-prefix based algorithm to minimize precision loss
	for (int z = zd_start; z < padsize_z; z++)
	{
		for (int y = 1; y < padsize_y; y++)
		{
			brent_kung_prefix_sum(&(VARBUF1(z, y, 1)), padsize_x - 1, 1);
			brent_kung_prefix_sum(&(VARBUF2(z, y, 1)), padsize_x - 1, 1);
		}
	}

	for (int z = zd_start; z < padsize_z; z++)
	{
		for (int x = 1; x < padsize_x; x++)
		{
			brent_kung_prefix_sum(&(VARBUF1(z, 1, x)), padsize_y - 1, yst);
			brent_kung_prefix_sum(&(VARBUF2(z, 1, x)), padsize_y - 1, yst);
		}
	}

	if (have_z)
	{
		for (int y = 1; y < padsize_y; y++)
		{
			for (int x = 1; x < padsize_x; x++)
			{
				brent_kung_prefix_sum(&(VARBUF1(1, y, x)), padsize_z - 1, zst);
				brent_kung_prefix_sum(&(VARBUF2(1, y, x)), padsize_z - 1, zst);
			}
		}
	}

	int avg_var_kdim = 2 * avg_var_kernel_radius + 1;
	int alpha_kdim = 2 * alpha_kernel_radius + 1;

	// Compute a few constants used in the variance-calculation.
	float avg_var_samples;
	float alpha_rsamples;
	float mul1;

	if (have_z)
	{
		avg_var_samples = (float)(avg_var_kdim * avg_var_kdim * avg_var_kdim);
		alpha_rsamples = 1.0f / (float)(alpha_kdim * alpha_kdim * alpha_kdim);
	}
	else
	{
		avg_var_samples = (float)(avg_var_kdim * avg_var_kdim);
		alpha_rsamples = 1.0f / (float)(alpha_kdim * alpha_kdim);
	}

	float avg_var_rsamples = 1.0f / avg_var_samples;
	if (avg_var_samples == 1)
	{
		mul1 = 1.0f;
	}
	else
	{
		mul1 = 1.0f / (float)(avg_var_samples * (avg_var_samples - 1));
	}

	float mul2 = avg_var_samples * mul1;

	// Use the summed-area tables to compute variance for each neighborhood
	if (have_z)
	{
		for (int z = 0; z < size_z; z++)
		{
			int z_src = z + kernel_radius_z;
			int z_dst = z + dst_offset_z;
			for (int y = 0; y < size_y; y++)
			{
				int y_src = y + kernel_radius_xy;
				int y_dst = y + dst_offset_y;

				for (int x = 0; x < size_x; x++)
				{
					int x_src = x + kernel_radius_xy;
					int x_dst = x + dst_offset_x;

					int x_low  = x_src - alpha_kernel_radius;
					int x_high = x_src + alpha_kernel_radius + 1;
					int y_low  = y_src - alpha_kernel_radius;
					int y_high = y_src + alpha_kernel_radius + 1;
					int z_low  = z_src - alpha_kernel_radius;
					int z_high = z_src + alpha_kernel_radius + 1;

					// Summed-area table lookups for alpha average
					float vasum = (  VARBUF1(z_high, y_low,  x_low).w
					               - VARBUF1(z_high, y_low,  x_high).w
					               - VARBUF1(z_high, y_high, x_low).w
					               + VARBUF1(z_high, y_high, x_high).w) -
					              (  VARBUF1(z_low,  y_low,  x_low).w
					               - VARBUF1(z_low,  y_low,  x_high).w
					               - VARBUF1(z_low,  y_high, x_low).w
					               + VARBUF1(z_low,  y_high, x_high).w);

					int out_index = z_dst * zdt + y_dst * ydt + x_dst;
					input_alpha_averages[out_index] = (vasum * alpha_rsamples);

					x_low  = x_src - avg_var_kernel_radius;
					x_high = x_src + avg_var_kernel_radius + 1;
					y_low  = y_src - avg_var_kernel_radius;
					y_high = y_src + avg_var_kernel_radius + 1;
					z_low  = z_src - avg_var_kernel_radius;
					z_high = z_src + avg_var_kernel_radius + 1;

					// Summed-area table lookups for RGBA average and variance
					float4 v1sum = (  VARBUF1(z_high, y_low,  x_low)
					                - VARBUF1(z_high, y_low,  x_high)
					                - VARBUF1(z_high, y_high, x_low)
					                + VARBUF1(z_high, y_high, x_high)) -
					               (  VARBUF1(z_low,  y_low,  x_low)
					                - VARBUF1(z_low,  y_low,  x_high)
					                - VARBUF1(z_low,  y_high, x_low)
					                + VARBUF1(z_low,  y_high, x_high));

					float4 v2sum = (  VARBUF2(z_high, y_low,  x_low)
					                - VARBUF2(z_high, y_low,  x_high)
					                - VARBUF2(z_high, y_high, x_low)
					                + VARBUF2(z_high, y_high, x_high)) -
					               (  VARBUF2(z_low,  y_low,  x_low)
					                - VARBUF2(z_low,  y_low,  x_high)
					                - VARBUF2(z_low,  y_high, x_low)
					                + VARBUF2(z_low,  y_high, x_high));

					// Compute and emit the average
					float4 avg = v1sum * avg_var_rsamples;
					input_averages[out_index] = avg;

					// Compute and emit the actual variance
					float4 variance = mul2 * v2sum - mul1 * (v1sum * v1sum);
					input_variances[out_index] = variance;
				}
			}
		}
	}
	else
	{
		for (int y = 0; y < size_y; y++)
		{
			int y_src = y + kernel_radius_xy;
			int y_dst = y + dst_offset_y;

			for (int x = 0; x < size_x; x++)
			{
				int x_src = x + kernel_radius_xy;
				int x_dst = x + dst_offset_x;

				int x_low  = x_src - alpha_kernel_radius;
				int x_high = x_src + alpha_kernel_radius + 1;
				int y_low  = y_src - alpha_kernel_radius;
				int y_high = y_src + alpha_kernel_radius + 1;

				// Summed-area table lookups for alpha average
				float vasum = VARBUF1(0, y_low,  x_low ).w
				            - VARBUF1(0, y_low,  x_high).w
				            - VARBUF1(0, y_high, x_low ).w
				            + VARBUF1(0, y_high, x_high).w;

				int out_index = y_dst * ydt + x_dst;
				input_alpha_averages[out_index] = (vasum * alpha_rsamples);

				x_low  = x_src - avg_var_kernel_radius;
				x_high = x_src + avg_var_kernel_radius + 1;
				y_low  = y_src - avg_var_kernel_radius;
				y_high = y_src + avg_var_kernel_radius + 1;

				// summed-area table lookups for RGBA average and variance
				float4 v1sum = VARBUF1(0, y_low,  x_low)
				             - VARBUF1(0, y_low,  x_high)
				             - VARBUF1(0, y_high, x_low)
				             + VARBUF1(0, y_high, x_high);

				float4 v2sum = VARBUF2(0, y_low,  x_low)
				             - VARBUF2(0, y_low,  x_high)
				             - VARBUF2(0, y_high, x_low)
				             + VARBUF2(0, y_high, x_high);

				// Compute and emit the average
				float4 avg = v1sum * avg_var_rsamples;
				input_averages[out_index] = avg;

				// Compute and emit the actual variance
				float4 variance = mul2 * v2sum - mul1 * (v1sum * v1sum);
				input_variances[out_index] = variance;
			}
		}
	}
}

static void compute_averages_and_variances_proc(
	int thread_count,
	int thread_id,
	void* payload
) {
	const struct avg_var_args *ag = (const struct avg_var_args *)payload;
	struct pixel_region_variance_args arg = ag->arg;
	arg.work_memory = new float4[ag->work_memory_size];

	int block_counter = 0;

	int size_x = ag->img_size.x;
	int size_y = ag->img_size.y;
	int size_z = ag->img_size.z;

	int step_x = ag->blk_size.x;
	int step_y = ag->blk_size.y;
	int step_z = ag->blk_size.z;

	int padding_xy = arg.img->padding;
	int padding_z = arg.have_z ? padding_xy : 0;

	for (int z = 0; z < size_z; z += step_z)
	{
		arg.size.z = MIN(step_z, size_z - z);
		arg.dst_offset.z = z;
		arg.src_offset.z = z + padding_z;

		for (int y = 0; y < size_y; y += step_y)
		{
			if (block_counter == thread_id)
			{
				arg.size.y = MIN(step_y, size_y - y);
				arg.dst_offset.y = y;
				arg.src_offset.y = y + padding_xy;

				for (int x = 0; x < size_x; x += step_x)
				{
						arg.size.x = MIN(step_x, size_x - x);
						arg.dst_offset.x = x;
						arg.src_offset.x = x + padding_xy;
						compute_pixel_region_variance(&arg);
				}
			}

			block_counter++;
			if (block_counter >= thread_count)
				block_counter = 0;
		}
	}

	delete[] arg.work_memory;
}

/* Public function, see header file for detailed documentation */
void compute_averages_and_variances(
	astc_codec_image* img,
	float rgb_power,
	float alpha_power,
	int avg_var_kernel_radius,
	int alpha_kernel_radius,
	int need_srgb_transform,
	swizzlepattern swz,
	int thread_count
) {
	int size_x = img->xsize;
	int size_y = img->ysize;
	int size_z = img->zsize;
	int pixel_count = size_x * size_y * size_z;

	// Perform memory allocations for the destination buffers
	if (img->input_averages)       delete[] img->input_averages;
	if (img->input_variances)      delete[] img->input_variances;
	if (img->input_alpha_averages) delete[] img->input_alpha_averages;

	img->input_averages = new float4[pixel_count];
	img->input_variances = new float4[pixel_count];
	img->input_alpha_averages = new float[pixel_count];

	// Compute maximum block size and from that the working memory buffer size
	int kernel_radius = MAX(avg_var_kernel_radius, alpha_kernel_radius);
	int kerneldim = 2 * kernel_radius + 1;

	int have_z = (size_z > 1);
	int max_blk_size_xy = have_z ? 16 : 32;
	int max_blk_size_z = MIN(size_z, have_z ? 16 : 1);

	int max_padsize_xy = max_blk_size_xy + kerneldim;
	int max_padsize_z = max_blk_size_z + (have_z ? kerneldim : 0);

	// Perform block-wise averages-and-variances calculations across the image
	struct pixel_region_variance_args arg;
	arg.img = img;
	arg.rgb_power = rgb_power;
	arg.alpha_power = alpha_power;
	arg.need_srgb_transform = need_srgb_transform;
	arg.swz = swz;
	arg.have_z = have_z;
	arg.avg_var_kernel_radius = avg_var_kernel_radius;
	arg.alpha_kernel_radius = alpha_kernel_radius;

	struct avg_var_args ag;
	ag.arg = arg;
	ag.img_size = int3(size_x, size_y, size_z);
	ag.blk_size = int3(max_blk_size_xy, max_blk_size_xy, max_blk_size_z);
	ag.work_memory_size = 2 * max_padsize_xy * max_padsize_xy * max_padsize_z;

	launch_threads(thread_count, compute_averages_and_variances_proc, &ag);
}
