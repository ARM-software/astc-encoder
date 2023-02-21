// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2023 Arm Limited
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
 * @brief Functions for the library entrypoint veneer layer.
 *
 * This veneer is a minimal shim compiled for a base instruction set that is
 * guaranteed to be supported, e.g. SSE2 on x86_64. This prevents use of any
 * unsupported enhanced instruction set before the ISA compatability error
 * check is performed.
 */

#include "astcenc.h"
#include "astcenc_internal_entry.h"

/* Non-veneer function implementation. */
astcenc_error astcenc_config_init_actual(
	astcenc_profile profile,
	unsigned int block_x,
	unsigned int block_y,
	unsigned int block_z,
	float quality,
	unsigned int flags,
	astcenc_config* config);

/* Non-veneer function implementation. */
astcenc_error astcenc_context_alloc_actual(
	const astcenc_config* config,
	unsigned int thread_count,
	astcenc_context** context);

/**
 * @brief Validate CPU ISA support meets the requirements of this build of the library.
 *
 * Each library build is statically compiled for a particular set of CPU ISA features, such as the
 * SIMD support or other ISA extensions such as POPCNT. This function checks that the host CPU
 * actually supports everything this build needs.
 *
 * @return Return @c ASTCENC_SUCCESS if validated, otherwise an error on failure.
 */
static astcenc_error validate_cpu_isa()
{
	#if ASTCENC_SSE >= 41
		if (!cpu_supports_sse41())
		{
			return ASTCENC_ERR_BAD_CPU_ISA;
		}
	#endif

	#if ASTCENC_POPCNT >= 1
		if (!cpu_supports_popcnt())
		{
			return ASTCENC_ERR_BAD_CPU_ISA;
		}
	#endif

	#if ASTCENC_F16C >= 1
		if (!cpu_supports_f16c())
		{
			return ASTCENC_ERR_BAD_CPU_ISA;
		}
	#endif

	#if ASTCENC_AVX >= 2
		if (!cpu_supports_avx2())
		{
			return ASTCENC_ERR_BAD_CPU_ISA;
		}
	#endif

	return ASTCENC_SUCCESS;
}

/* See header for documentation. */
astcenc_error astcenc_config_init(
	astcenc_profile profile,
	unsigned int block_x,
	unsigned int block_y,
	unsigned int block_z,
	float quality,
	unsigned int flags,
	astcenc_config* configp
) {
	// Check ISA compatability in the veneer before handing off
	astcenc_error status = validate_cpu_isa();
	if (status != ASTCENC_SUCCESS)
	{
		return status;
	}

	cpu_isb();

	return astcenc_config_init_actual(
		profile, block_x, block_y, block_z, quality, flags, configp);
}

/* See header for documentation. */
astcenc_error astcenc_context_alloc(
	const astcenc_config* configp,
	unsigned int thread_count,
	astcenc_context** context
) {
	// Check ISA compatability in the veneer before handing off
	astcenc_error status = validate_cpu_isa();
	if (status != ASTCENC_SUCCESS)
	{
		return status;
	}

	cpu_isb();

	return astcenc_context_alloc_actual(
		configp, thread_count, context);
}
