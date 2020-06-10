# SPDX-License-Identifier: Apache-2.0
# -----------------------------------------------------------------------------
# Copyright 2019-2020 Arm Limited
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License. You may obtain a copy
# of the License at:
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations
# under the License.
# -----------------------------------------------------------------------------

# Configure the vectorization intrinsics support; valid values are:
#
#   * nointrin - allow use of sse2 by the compiler, but no manual instrinsics
#   * sse2     - allow use of sse2
#   * sse4.2   - allow use of sse4.2 and popcnt
#   * avx2     - allow use of avx2, sse4.2, and popcnt
#
# Note that we always enable at least sse2 support for compiler generated code,
# as it is guaranteed to be present in x86-64, we only allow disabling of our
# manually written instrinsic functions.
#
# Also note that we currently assume that sse4.2 implies support for popcnt,
# which is the default GCC behavior and true on currently shipping hardware.
VEC?=avx2

DBG?=0

SOURCES_LIBRARY = \
    astcenc_averages_and_directions.cpp \
    astcenc_block_sizes2.cpp \
    astcenc_color_quantize.cpp \
    astcenc_color_unquantize.cpp \
    astcenc_compress_symbolic.cpp \
    astcenc_compute_variance.cpp \
    astcenc_decompress_symbolic.cpp \
    astcenc_encoding_choice_error.cpp \
    astcenc_entry.cpp \
    astcenc_find_best_partitioning.cpp \
    astcenc_ideal_endpoints_and_weights.cpp \
    astcenc_image.cpp \
    astcenc_integer_sequence.cpp \
    astcenc_kmeans_partitioning.cpp \
    astcenc_mathlib.cpp \
    astcenc_mathlib_softfloat.cpp \
    astcenc_partition_tables.cpp \
    astcenc_percentile_tables.cpp \
    astcenc_pick_best_endpoint_format.cpp \
    astcenc_platform_dependents.cpp \
    astcenc_quantization.cpp \
    astcenc_symbolic_physical.cpp \
    astcenc_weight_align.cpp \
    astcenc_weight_quant_xfer_tables.cpp

SOURCES_CLI = \
    astcenc_mathlib.cpp \
    astcenc_mathlib_softfloat.cpp \
    astcenccli_error_metrics.cpp \
    astcenccli_image.cpp \
    astcenccli_image_load_store.cpp \
    astcenccli_main.cpp \
    astcenccli_platform_isa_detection.cpp \
    astcenccli_platform_dependents.cpp \
    astcenccli_toplevel.cpp \
    astcenccli_toplevel_help.cpp

EXTERNAL_SOURCES = \
    stb_image.h \
    stb_image_write.h \
    tinyexr.h

HEADERS = \
    astcenc.h \
    astcenc_internal.h \
    astcenc_mathlib.h \
    stb_image.h \
    stb_image_write.h \
    tinyexr.h

OBJECTS_LIBRARY = $(SOURCES_LIBRARY:.cpp=-l-$(VEC).o)

OBJECTS_CLI = $(SOURCES_CLI:.cpp=-c-$(VEC).o)

EXTERNAL_OBJECTS = $(EXTERNAL_SOURCES:.h=-c-$(VEC).o)

BINARY = astcenc-$(VEC)

LIBRARY = libastcenc-$(VEC).so

CXXFLAGS = -std=c++14 -fvisibility=hidden -mfpmath=sse \
           -Wall -Wextra -Wpedantic -Werror -Werror=shadow -Wdouble-promotion

# Validate that the DBG parameter is a supported value, and patch CXXFLAGS
ifeq ($(DBG),0)
CXXFLAGS += -O3 -DNDEBUG
else
CXXFLAGS += -O0 -g
endif

# Validate that the VEC parameter is a supported value, and patch CXXFLAGS
ifeq ($(VEC),nointrin)
CXXFLAGS += -msse2 -DASTC_SSE=0 -DASTC_AVX=0 -DASTC_POPCNT=0
else
ifeq ($(VEC),sse2)
CXXFLAGS += -msse2 -DASTC_SSE=20 -DASTC_AVX=0 -DASTC_POPCNT=0
else
ifeq ($(VEC),sse4.2)
CXXFLAGS += -msse4.2 -mpopcnt -DASTC_SSE=42 -DASTC_AVX=0 -DASTC_POPCNT=1
else
ifeq ($(VEC),avx2)
CXXFLAGS += -mavx2 -mpopcnt -DASTC_SSE=42 -DASTC_AVX=2 -DASTC_POPCNT=1
else
$(error Unsupported VEC target, use VEC=nointrin/sse2/sse4.2/avx2)
endif
endif
endif
endif

# Disable necessary optimizations and warnings for third-party source files
CXXFLAGS_EXTERNAL = \
    $(CXXFLAGS) \
    -Wno-unused-parameter \
    -Wno-double-promotion \
    -fno-strict-aliasing

# ==================================================
# Build rules for the command line wrapper
$(BINARY): $(EXTERNAL_OBJECTS) $(OBJECTS_CLI) | $(LIBRARY)
	@$(CXX) -o $@ $^ $(CXXFLAGS) -L. -l$(BINARY) -lpthread -Wl,-rpath,'$$ORIGIN'
	@echo "[Link] $@ (using $(VEC), debug=$(DBG))"

# Note: ensure NDEBUG is undefined; all three libraries are weak at runtime
# handling of corrupt files, relying on asserts to handle bd file input.
# Compiling with NDEBUG is a potential security risk.
stb_image-c-$(VEC).o: stb_image.h
	@$(CXX) -c -x c++ -o $@ $< $(CXXFLAGS_EXTERNAL) -DSTB_IMAGE_IMPLEMENTATION \
		-DSTBI_NO_PSD -DSTBI_NO_GIF -DSTBI_NO_PIC -DSTBI_NO_PNM -UNDEBUG
	@echo "[C++] $<"

stb_image_write-c-$(VEC).o: stb_image_write.h
	@$(CXX) -c -x c++ -o $@ $< $(CXXFLAGS_EXTERNAL) -DSTB_IMAGE_WRITE_IMPLEMENTATION \
		-UNDEBUG
	@echo "[C++] $<"

tinyexr-c-$(VEC).o: tinyexr.h
	@$(CXX) -c -x c++ -o $@ $< $(CXXFLAGS) -DTINYEXR_IMPLEMENTATION -UNDEBUG
	@echo "[C++] $<"

$(OBJECTS_CLI): %-c-$(VEC).o: %.cpp $(HEADERS)
	@$(CXX) -c -o $@ $(abspath $<) $(CXXFLAGS)
	@echo "[C++] $<"

# ==================================================
# Build rules for the codec library shared object
$(LIBRARY): $(OBJECTS_LIBRARY)
	@$(CXX) -shared -o $@ $^ $(CXXFLAGS)
	@echo "[Link] $@ (using $(VEC), debug=$(DBG))"

$(OBJECTS_LIBRARY): %-l-$(VEC).o: %.cpp $(HEADERS)
	@$(CXX) -c -fpic -o $@ $(abspath $<) $(CXXFLAGS)
	@echo "[C++] $<"

# ==================================================
# Generic utility rules

clean:
	@rm -f *-$(VEC).o
	@rm -f *-$(VEC).so
	@rm -f *-$(VEC).obj
	@rm -f $(BINARY)
	@rm -f $(BINARY).exe
	@echo "[Clean] Removing $(VEC) build outputs"

batchbuild:
	@+$(MAKE) -s VEC=nointrin
	@+$(MAKE) -s VEC=sse2
	@+$(MAKE) -s VEC=sse4.2
	@+$(MAKE) -s VEC=avx2

batchclean:
	@+$(MAKE) -s clean VEC=nointrin
	@+$(MAKE) -s clean VEC=sse2
	@+$(MAKE) -s clean VEC=sse4.2
	@+$(MAKE) -s clean VEC=avx2
