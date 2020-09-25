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

/**
 * @brief Fuzz target for physical_to_symbolic.
 *
 * This function is the main entrypoint for decompressing a 16 byte block of
 * input ASTC data from disk. The 16 bytes can contain arbitrary data; they
 * are read from an external source, but the block size used must be a valid
 * ASTC block footprint.
 *
 * Fuzz design improvements - this function itself is quite small, but it
 * requires use to create a BSD, which actually takes quite a bit of time
 * compared to physical_to_symbolic function we want to fuzz. Possible to
 * avoid that?
 */

#include "astcenc_internal.h"

#include <fuzzer/FuzzedDataProvider.h>
#include <array>
#include <vector>

struct BlockSizes {
	int x;
	int y;
	int z;
};

std::array<BlockSizes, 3> testSz {{
	{ 4,  4, 1}, // Highest bitrate
	{12, 12, 1}, // Largest 2D block
	{6,  6,  6}  // Largest 3D block
}};

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
	// Must have 4 (int - block size select) and 16 (blob - payload) bytes
	if (size < 4 + 16)
	{
		return 0;
	}

	FuzzedDataProvider stream(data, size);

	// Select a block size to test
	// TODO - building a BSD is slow compared to the physical_to_symbolic step.
	//        Is it possible to cache these without triggering leaksanitizer?
	block_size_descriptor bsd;
	int i = stream.ConsumeIntegralInRange<int>(0, testSz.size() - 1);
	init_block_size_descriptor(testSz[i].x, testSz[i].y, testSz[i].z, &bsd);

	physical_compressed_block pcb;
	std::vector<uint8_t> buffer = stream.ConsumeBytes<uint8_t>(16);
	std::memcpy(&pcb, buffer.data(), 16);

	symbolic_compressed_block scb;
	physical_to_symbolic(&bsd, pcb, &scb);

	term_block_size_descriptor(&bsd);
	return 0;
}
