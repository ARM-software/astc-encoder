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

/**
 * @brief Functions for the library entrypoint.
 */

#if defined(ASTCENC_DIAGNOSTICS)

#include <cassert>
#include <cstdarg>
#include <cstdio>

#include "astcenc_diagnostic_trace.h"


/** @brief The global trace logger. */
static TraceLog* g_TraceLog;

TraceLog::TraceLog(const char* file_name): m_file(file_name)
{
	assert(!g_TraceLog);
	g_TraceLog = this;
}

TraceLog::~TraceLog()
{
	assert(g_TraceLog == this);
	g_TraceLog = nullptr;
}

TraceNode::TraceNode(const char* format, ...)
{
	constexpr size_t bufsz = 256;
	char buffer[bufsz];

	va_list args;
	va_start (args, format);
	vsnprintf (buffer, bufsz, format, args);
	va_end (args);

	// Guarantee there is a nul termintor
	buffer[bufsz - 1] = 0;

	printf("%s\n", buffer);
}

TraceNode::~TraceNode()
{

}

void trace_add_data(const char* key, const char* format, ...)
{
	constexpr size_t bufsz = 256;
	char buffer[bufsz];

	va_list args;
	va_start (args, format);
	vsnprintf (buffer, bufsz, format, args);
	va_end (args);

	// Guarantee there is a nul termintor
	buffer[bufsz - 1] = 0;

	printf("%s => %s\n", key, buffer);
}

void trace_add_data(const char* key, float value)
{
	printf("%s => %g\n", key, (double)value);
}

void trace_add_data(const char* key, int value)
{
	printf("%s => %d\n", key, value);
}

void trace_add_data(const char* key, unsigned int value)
{
	printf("%s => %u\n", key, value);
}

#endif
