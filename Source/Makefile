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

# Configure the application to build
#
#   * astcenc - a full compressor/decompressor
#   * astcdec - a decompressor without compression support
#
# Note that "astcdec" is not a fully customized command line application. It
# will look like a normal astcenc build except for the binary name, and all
# operations requring a compression will fail. It is only intended for use as
# test vehicle for the decompress-only library testing.
APP ?= astcenc

# Configure the SIMD ISA and intrinsics support; valid values are:
#
#   * neon   - compile for Arm aarch64 + NEON
#   * sse2   - allow use of x86-64 + sse2
#   * sse4.2 - allow use of x86-64 + sse4.2 and popcnt
#   * avx2   - allow use of x86-64 + avx2, sse4.2, and popcnt
VEC ?= avx2

# Configure the build type; valid values are:
#
#   * release - build a fully optimized build without symbols
#   * profile - build a moderately optimized build, with symbols
#   * debug   - build an unoptimized build with symbols
BUILD ?= release

# ==================================================
# File listings

SOURCES = \
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
    astcenc_platform_isa_detection.cpp \
    astcenc_quantization.cpp \
    astcenc_symbolic_physical.cpp \
    astcenc_weight_align.cpp \
    astcenc_weight_quant_xfer_tables.cpp \
    astcenccli_error_metrics.cpp \
    astcenccli_image.cpp \
    astcenccli_image_load_store.cpp \
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
    astcenc_vecmathlib.h \
    astcenccli_internal.h \
    stb_image.h \
    stb_image_write.h \
    tinyexr.h

BINARY = $(APP)-$(VEC)

OBJECTS = $(SOURCES:.cpp=-$(BINARY).o)

EXTERNAL_OBJECTS = $(EXTERNAL_SOURCES:.h=-$(BINARY).o)

# ==================================================
# CXXFLAGS setup (and input validation). This appends on top of whatever
# CXXFLAGS are set on the command line, which is useful for integration
# with the oss-fuzz build infrastructure for fuzz testing

CXXFLAGS += -std=c++14 -fvisibility=hidden \
            -Wall -Wextra -Wpedantic -Werror -Werror=shadow -Wdouble-promotion

# Validate that the APP parameter is a supported value, and patch CXXFLAGS
ifeq ($(APP),astcenc)
# Nothing to set up
else
ifeq ($(APP),astcdec)
CXXFLAGS += -DASTCENC_DECOMPRESS_ONLY
else
$(error Unsupported app, use APP=astcenc/astcdec)
endif
endif

# Validate that the BUILD parameter is a supported value, and patch CXXFLAGS
ifeq ($(BUILD),release)
CXXFLAGS += -O3 -DNDEBUG -flto
else
ifeq ($(BUILD),profile)
CXXFLAGS += -O2 -g
else
ifeq ($(BUILD),debug)
CXXFLAGS += -O0 -g
else
$(error Unsupported build, use BUILD=release/profile/debug)
endif
endif
endif

# Validate that the VEC parameter is a supported value, and patch CXXFLAGS
ifeq ($(VEC),neon)
# NEON is on by default; no enabled needed
CXXFLAGS += -DASTCENC_SSE=0 -DASTCENC_AVX=0 -DASTCENC_POPCNT=0 -DASTCENC_VECALIGN=16
else
ifeq ($(VEC),sse2)
CXXFLAGS += -mfpmath=sse -msse2 -DASTCENC_SSE=20 -DASTCENC_AVX=0 -DASTCENC_POPCNT=0 -DASTCENC_VECALIGN=16
else
ifeq ($(VEC),sse4.2)
CXXFLAGS += -mfpmath=sse -msse4.2 -mpopcnt -DASTCENC_SSE=42 -DASTCENC_AVX=0 -DASTCENC_POPCNT=1 -DASTCENC_VECALIGN=16
else
ifeq ($(VEC),avx2)
CXXFLAGS += -mfpmath=sse -mavx2 -mpopcnt -DASTCENC_SSE=42 -DASTCENC_AVX=2 -DASTCENC_POPCNT=1 -DASTCENC_VECALIGN=32
else
$(error Unsupported VEC target, use VEC=neon/sse2/sse4.2/avx2)
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
# Exports for child make invocations

export CXX
export BUILD

# ==================================================
# Build rules for the command line wrapper

$(BINARY): $(EXTERNAL_OBJECTS) $(OBJECTS)
	@$(CXX) -o $@ $^ $(CXXFLAGS) -lpthread
	@echo "[Link] $@ (using $(VEC) $(BUILD))"

# Note: ensure NDEBUG is undefined; all three libraries are weak at runtime
# handling of corrupt files, relying on asserts to handle bad file input.
# Compiling with NDEBUG is a potential security risk.
stb_image-$(BINARY).o: stb_image.h
	@$(CXX) -c -x c++ -o $@ $< $(CXXFLAGS_EXTERNAL) -DSTB_IMAGE_IMPLEMENTATION \
		-DSTBI_NO_PSD -DSTBI_NO_GIF -DSTBI_NO_PIC -DSTBI_NO_PNM -UNDEBUG
	@echo "[C++] $<"

stb_image_write-$(BINARY).o: stb_image_write.h
	@$(CXX) -c -x c++ -o $@ $< $(CXXFLAGS_EXTERNAL) -DSTB_IMAGE_WRITE_IMPLEMENTATION \
		-UNDEBUG
	@echo "[C++] $<"

tinyexr-$(BINARY).o: tinyexr.h
	@$(CXX) -c -x c++ -o $@ $< $(CXXFLAGS) -DTINYEXR_IMPLEMENTATION -UNDEBUG
	@echo "[C++] $<"

$(OBJECTS): %-$(BINARY).o: %.cpp $(HEADERS)
	@$(CXX) -c -o $@ $(abspath $<) $(CXXFLAGS)
	@echo "[C++] $<"

# ==================================================
# Generic utility rules

clean:
	@rm -f *-$(BINARY).o
	@rm -f *-$(BINARY).obj
	@rm -f $(BINARY)
	@rm -f $(BINARY).exe
	@echo "[Clean] Removing $(VEC) build outputs"

batchbuild:
	# NEON is deliberately not here - needs a different CXX setting
	@+$(MAKE) -s VEC=sse2
	@+$(MAKE) -s VEC=sse4.2
	@+$(MAKE) -s VEC=avx2

batchclean:
	@+$(MAKE) -s clean VEC=neon
	@+$(MAKE) -s clean VEC=sse2
	@+$(MAKE) -s clean VEC=sse4.2
	@+$(MAKE) -s clean VEC=avx2
