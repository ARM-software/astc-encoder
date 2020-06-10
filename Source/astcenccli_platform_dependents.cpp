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

#include "astcenccli_internal.h"

/* ============================================================================
   Platform code for Windows using the Win32 APIs.
============================================================================ */
#if defined(_WIN32) && !defined(__CYGWIN__)

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

/* Public function, see header file for detailed documentation */
double get_time()
{
	FILETIME tv;
	GetSystemTimeAsFileTime(&tv);
	unsigned long long ticks = tv.dwHighDateTime;
	ticks = (ticks << 32) | tv.dwLowDateTime;
	return ((double)ticks) / 1.0e7;
}

/* Public function, see header file for detailed documentation */
int get_cpu_count()
{
	SYSTEM_INFO sysinfo;
	GetSystemInfo(&sysinfo);
	return sysinfo.dwNumberOfProcessors;
}

/* Public function, see header file for detailed documentation */
int unlink_file(const char *filename)
{
	BOOL res = DeleteFileA(filename);
	return (res ? 0 : -1);
}

/* ============================================================================
   Platform code for an platform using POSIX APIs.
============================================================================ */
#else

#include <sys/time.h>
#include <unistd.h>

/* Public function, see header file for detailed documentation */
double get_time()
{
	timeval tv;
	gettimeofday(&tv, 0);
	return (double)tv.tv_sec + (double)tv.tv_usec * 1.0e-6;
}

/* Public function, see header file for detailed documentation */
int get_cpu_count()
{
	return sysconf(_SC_NPROCESSORS_ONLN);
}

/* Public function, see header file for detailed documentation */
int unlink_file(const char *filename)
{
	return unlink(filename);
}

#endif
