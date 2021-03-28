// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2020-2021 Arm Limited
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
 * @brief Platform-specific function implementations.
 *
 * This module contains functions for querying the host extended ISA support.
 */

// Include before the defines below to pick up any auto-setup based on compiler
// built-in config, if not being set explicitly by the build system
#include "astcenc_internal.h"

#if (ASTCENC_SSE > 0)    || (ASTCENC_AVX > 0) || \
    (ASTCENC_POPCNT > 0) || (ASTCENC_F16C > 0)

static int g_cpu_has_sse41 = -1;
static int g_cpu_has_avx2 = -1;
static int g_cpu_has_popcnt = -1;
static int g_cpu_has_f16c = -1;

/* ============================================================================
   Platform code for Visual Studio
============================================================================ */
#if !defined(__clang__) && defined(_MSC_VER)
#include <intrin.h>

static void detect_cpu_isa()
{
	int data[4];

	__cpuid(data, 0);
	int num_id = data[0];

	g_cpu_has_sse41 = 0;
	g_cpu_has_popcnt = 0;
	g_cpu_has_f16c = 0;
	if (num_id >= 1)
	{
		__cpuidex(data, 1, 0);
		// SSE41 = Bank 1, ECX, bit 19
		g_cpu_has_sse41 = data[2] & (1 << 19) ? 1 : 0;
		// POPCNT = Bank 1, ECX, bit 23
		g_cpu_has_popcnt = data[2] & (1 << 23) ? 1 : 0;
		// F16C = Bank 1, ECX, bit 29
		g_cpu_has_f16c = data[2] & (1 << 29) ? 1 : 0;
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

	g_cpu_has_sse41 = 0;
	g_cpu_has_popcnt = 0;
	g_cpu_has_f16c = 0;
	if (__get_cpuid_count(1, 0, &data[0], &data[1], &data[2], &data[3]))
	{
		// SSE41 = Bank 1, ECX, bit 19
		g_cpu_has_sse41 = data[2] & (1 << 19) ? 1 : 0;
		// POPCNT = Bank 1, ECX, bit 23
		g_cpu_has_popcnt = data[2] & (1 << 23) ? 1 : 0;
		// F16C = Bank 1, ECX, bit 29
		g_cpu_has_f16c = data[2] & (1 << 29) ? 1 : 0;
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
int cpu_supports_sse41()
{
	if (g_cpu_has_sse41 == -1)
	{
		detect_cpu_isa();
	}

	return g_cpu_has_sse41;
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
int cpu_supports_f16c()
{
	if (g_cpu_has_f16c == -1)
	{
		detect_cpu_isa();
	}

	return g_cpu_has_f16c;
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
