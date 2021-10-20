// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2021 Arm Limited
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

// This is a utility tool to measure angular errors in a normal map.

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "astcenc_mathlib.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define MODE_ENCODE 0
#define MODE_DECODE 1

int main(int argc, char **argv)
{
	// Parse command line
	if (argc != 3)
	{
		printf("Usage: astc_angular_error <reference> <test>\n");
		exit(1);
	}

	const char* ref_file = argv[1];
	int dim_ref_x;
	int dim_ref_y;
	const uint8_t* data_ref = stbi_load(ref_file, &dim_ref_x, &dim_ref_y, nullptr, 4);
	if (!data_ref)
	{
		printf("ERROR: Failed to load reference image.\n");
		exit(1);
	}

	const char* tst_file = argv[2];
	int dim_tst_x;
	int dim_tst_y;
	const uint8_t* data_tst = stbi_load(tst_file, &dim_tst_x, &dim_tst_y, nullptr, 4);
	if (!data_tst)
	{
		printf("ERROR: Failed to load test image.\n");
		exit(1);
	}

	if (dim_ref_x != dim_tst_x || dim_ref_y != dim_tst_y)
	{
		printf("ERROR: Reference and test images are of different dimensions.\n");
		exit(1);
	}

	double mean_angular_error = 0.0;
	double total_angular_error = 0.0;

	// For each pixel compute angular error
	for (int y = 0; y < dim_ref_y; y++)
	{
		const uint8_t* row_ref = data_ref + (4 * dim_ref_x * y);
		const uint8_t* row_tst = data_tst + (4 * dim_tst_x * y);

		for (int x = 0; x < dim_ref_x; x++)
		{
			const uint8_t* pixel_ref = row_ref + 4 * x;
			const uint8_t* pixel_tst = row_tst + 4 * x;

			// Compute normal components in -1 to +1 range
			float x_ref = ((static_cast<float>(pixel_ref[0]) / 255.0f) - 0.5f) * 2.0f;
			float y_ref = ((static_cast<float>(pixel_ref[1]) / 255.0f) - 0.5f) * 2.0f;
			float z_ref = ((static_cast<float>(pixel_ref[2]) / 255.0f) - 0.5f) * 2.0f;

			float x_tst = ((static_cast<float>(pixel_tst[0]) / 255.0f) - 0.5f) * 2.0f;
			float y_tst = ((static_cast<float>(pixel_tst[1]) / 255.0f) - 0.5f) * 2.0f;
			float z_tst = ((static_cast<float>(pixel_tst[2]) / 255.0f) - 0.5f) * 2.0f;

			vfloat4 ref(x_ref, y_ref, z_ref, 0.0f);
			ref = normalize_safe(ref, unit3());

			vfloat4 tst(x_tst, y_tst, z_tst, 0.0f);
			tst = normalize_safe(tst, unit3());

			float angle = std::acos(dot3_s(ref, tst)) * (180.0f / astc::PI);
			if (std::isnan(angle))
			{
				continue;
			}


			printf("%4.4f\n", angle);
			mean_angular_error += static_cast<double>(angle) / (dim_ref_x * dim_ref_y);
			total_angular_error += static_cast<double>(angle);
		}
	}

	printf("\nMean angular error:  %4.4f degrees\n", mean_angular_error);
	printf("\nTotal angular error: %4.4f degrees\n", total_angular_error);

	return 0;
}
