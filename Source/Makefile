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

SOURCES = \
    astc_averages_and_directions.cpp \
    astc_block_sizes2.cpp \
    astc_color_quantize.cpp \
    astc_color_unquantize.cpp \
    astc_compress_symbolic.cpp \
    astc_compute_variance.cpp \
    astc_decompress_symbolic.cpp \
    astc_encoding_choice_error.cpp \
    astc_error_metrics.cpp \
    astc_find_best_partitioning.cpp \
    astc_ideal_endpoints_and_weights.cpp \
    astc_image_load_store.cpp \
    astc_integer_sequence.cpp \
    astc_kmeans_partitioning.cpp \
    astc_ktx_dds.cpp \
    astc_main.cpp \
    astc_mathlib.cpp \
    astc_mathlib_softfloat.cpp \
    astc_partition_tables.cpp \
    astc_percentile_tables.cpp \
    astc_pick_best_endpoint_format.cpp \
    astc_platform_isa_detection.cpp \
    astc_platform_dependents.cpp \
    astc_quantization.cpp \
    astc_stb_tga.cpp \
    astc_symbolic_physical.cpp \
    astc_toplevel.cpp \
    astc_toplevel_help.cpp \
    astc_weight_align.cpp \
    astc_weight_quant_xfer_tables.cpp

EXTERNAL_SOURCES = \
    stb_image.h \
    stb_image_write.h \
    tinyexr.h

HEADERS = \
    astc_codec_internals.h \
    astc_mathlib.h \
    stb_image.h \
    stb_image_write.h \
    tinyexr.h

OBJECTS = $(SOURCES:.cpp=-$(VEC).o)

EXTERNAL_OBJECTS = $(EXTERNAL_SOURCES:.h=-$(VEC).o)

BINARY = astcenc-$(VEC)

CXXFLAGS = -std=c++14 -mfpmath=sse \
           -Wall -Wextra -Wpedantic -Werror -Werror=shadow -Wdouble-promotion

# Validate that the DBG parameter is a supported value, and patch CXXFLAGS
ifeq ($(DBG),0)
CXXFLAGS += -O3
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

$(BINARY): $(EXTERNAL_OBJECTS) $(OBJECTS)
	@$(CXX) -o $@ $^ $(CXXFLAGS) -lpthread
	@echo "[Link] $@ (using $(VEC), debug=$(DBG))"

stb_image-$(VEC).o: stb_image.h
	@$(CXX) -c -x c++ -o $@ $< $(CXXFLAGS_EXTERNAL) -DSTB_IMAGE_IMPLEMENTATION \
		-DSTBI_NO_PSD -DSTBI_NO_GIF -DSTBI_NO_PIC -DSTBI_NO_PNM
	@echo "[C++] $<"

stb_image_write-$(VEC).o: stb_image_write.h
	@$(CXX) -c -x c++ -o $@ $< $(CXXFLAGS_EXTERNAL) -DSTB_IMAGE_WRITE_IMPLEMENTATION
	@echo "[C++] $<"

tinyexr-$(VEC).o: tinyexr.h
	@$(CXX) -c -x c++ -o $@ $< $(CXXFLAGS) -DTINYEXR_IMPLEMENTATION
	@echo "[C++] $<"

$(OBJECTS): %-$(VEC).o: %.cpp $(HEADERS)
	@$(CXX) -c -o $@ $< $(CXXFLAGS)
	@echo "[C++] $<"

clean:
	@rm -f *-$(VEC).o
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
