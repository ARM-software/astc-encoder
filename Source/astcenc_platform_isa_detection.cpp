// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2020 Arm Limited
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

#if (ASTCENC_SSE > 0) || (ASTCENC_AVX > 0) || (ASTCENC_POPCNT > 0)

/**
 * @brief Platform-specific function implementations.
 *
 * This module contains functions for querying the host extended ISA support.
 */

#include "astcenc_internal.h"

static int g_cpu_has_sse42 = -1;
static int g_cpu_has_avx2 = -1;
static int g_cpu_has_popcnt = -1;

/* ============================================================================
   Platform code for Visual Studio
============================================================================ */
#if defined(_MSC_VER)
#include <intrin.h>

static void detect_cpu_isa()
{
	int data[4];

	__cpuid(data, 0);
	int num_id = data[0];

	g_cpu_has_sse42 = 0;
	g_cpu_has_popcnt = 0;
	if (num_id >= 1)
	{
		__cpuidex(data, 1, 0);
		// SSE42 = Bank 1, ECX, bit 20
		g_cpu_has_sse42 = data[2] & (1 << 20) ? 1 : 0;
		// POPCNT = Bank 1, ECX, bit 23
		g_cpu_has_popcnt = data[2] & (1 << 23) ? 1 : 0;
	}

	g_cpu_has_avx2 = 0;
	if (num_id >= 7)
	{
		__cpuidex(data, 7, 0);
		// AVX2 = Bank 7, EBX, bit 5
		g_cpu_has_avx2 = data[1] & (1 << 5) ? 1 : 0;
	}
}

/* ============================================================================
   Platform code for GCC and Clang
============================================================================ */
#else
#include <cpuid.h>

static void detect_cpu_isa()
{
	unsigned int data[4];

	g_cpu_has_sse42 = 0;
	g_cpu_has_popcnt = 0;
	if (__get_cpuid_count(1, 0, &data[0], &data[1], &data[2], &data[3]))
	{
		// SSE42 = Bank 1, ECX, bit 20
		g_cpu_has_sse42 = data[2] & (1 << 20) ? 1 : 0;
		// POPCNT = Bank 1, ECX, bit 23
		g_cpu_has_popcnt = data[2] & (1 << 23) ? 1 : 0;
	}

	g_cpu_has_avx2 = 0;
	if (__get_cpuid_count(7, 0, &data[0], &data[1], &data[2], &data[3]))
	{
		// AVX2 = Bank 7, EBX, bit 5
		g_cpu_has_avx2 = data[1] & (1 << 5) ? 1 : 0;
	}
}
#endif

/* Public function, see header file for detailed documentation */
int cpu_supports_sse42()
{
	if (g_cpu_has_sse42 == -1)
	{
		detect_cpu_isa();
	}

	return g_cpu_has_sse42;
}

/* Public function, see header file for detailed documentation */
int cpu_supports_popcnt()
{
	if (g_cpu_has_popcnt == -1)
	{
		detect_cpu_isa();
	}

	return g_cpu_has_popcnt;
}

/* Public function, see header file for detailed documentation */
int cpu_supports_avx2()
{
	if (g_cpu_has_avx2 == -1)
	{
		detect_cpu_isa();
	}

	return g_cpu_has_avx2;
}

#endif
