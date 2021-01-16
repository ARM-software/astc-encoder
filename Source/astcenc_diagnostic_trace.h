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
 * @brief This module provides a set of diagnostic tracing utilities.
 */

#ifndef ASTCENC_DIAGNOSTIC_TRACE_INCLUDED
#define ASTCENC_DIAGNOSTIC_TRACE_INCLUDED

#if defined(ASTCENC_DIAGNOSTICS)

#include <iostream>
#include <fstream>
#include <vector>

class TraceNode
{
public:
	TraceNode(const char* format, ...);
	void add_attrib(std::string type, std::string key, std::string value);
	~TraceNode();
	unsigned int m_attrib_count { 0 };
};


class TraceLog
{
public:
	TraceLog(const char* file_name);
	~TraceLog();

	TraceNode* get_current_leaf();
	int get_depth();

 	std::ofstream m_file;
	std::vector<TraceNode*> m_stack;
	TraceNode* m_root;
};

#define TRACE_NODE(name, ...) TraceNode name(__VA_ARGS__);

void trace_add_data(const char* key, const char* format, ...);

void trace_add_data(const char* key, float value);

void trace_add_data(const char* key, int value);

void trace_add_data(const char* key, unsigned int value);

#else

#define TRACE_NODE(name, ...)

#define trace_add_data(...)

#endif

#endif
