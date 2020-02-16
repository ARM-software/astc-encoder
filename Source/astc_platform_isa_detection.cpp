// ----------------------------------------------------------------------------
//  This confidential and proprietary software may be used only as authorised
//  by a licensing agreement from Arm Limited.
//      (C) COPYRIGHT 2011-2020 Arm Limited, ALL RIGHTS RESERVED
//  The entire notice above must be reproduced on all authorised copies and
//  copies may only be made to the extent permitted by a licensing agreement
//  from Arm Limited.
// ----------------------------------------------------------------------------

/**
 * @brief Platform-specific function implementations.
 *
 * This module contains functions with strongly OS-dependent implementations:
 *
 *  * CPU count queries
 *  * Threading
 *  * Time
 *
 * In addition to the basic thread abstraction (which is native pthreads on
 * all platforms, except Windows where it is an emulation of pthreads), a
 * utility function to create N threads and wait for them to complete a batch
 * task has also been provided.
 */

#include "astc_codec_internals.h"

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

	if (num_id >= 1)
	{
		__cpuidex(data, 1, 0);
		// SSE42 = Bank 1, ECX, bit 20
		g_cpu_has_sse42 = data[2] & (1 << 20) ? 1 : 0;
		// POPCNT = Bank 1, ECX, bit 23
		g_cpu_has_popcnt = data[2] & (1 << 23) ? 1 : 0;
	}

	if (num_id >= 7) {
		__cpuidex(data, 7, 0);
		// AVX2 = Bank 7, EBX, bit 5
		g_cpu_has_avx2 = data[1] & (1 << 5) ? 1 : 0;
	}
}

/* ============================================================================
   Platform code for GCC and Clang
============================================================================ */
#else
static void detect_cpu_isa()
{
	g_cpu_has_sse42 = __builtin_cpu_supports("sse4.2") ? 1 : 0;
	g_cpu_has_popcnt = __builtin_cpu_supports("popcnt") ? 1 : 0;
	g_cpu_has_avx2 = __builtin_cpu_supports("avx2") ? 1 : 0;
}
#endif

/* Public function, see header file for detailed documentation */
int cpu_supports_sse42()
{
	if (g_cpu_has_sse42 == -1)
		detect_cpu_isa();

	return g_cpu_has_sse42;
}

/* Public function, see header file for detailed documentation */
int cpu_supports_popcnt()
{
	if (g_cpu_has_popcnt == -1)
		detect_cpu_isa();

	return g_cpu_has_popcnt;
}

/* Public function, see header file for detailed documentation */
int cpu_supports_avx2()
{
	if (g_cpu_has_avx2 == -1)
		detect_cpu_isa();

	return g_cpu_has_avx2;
}