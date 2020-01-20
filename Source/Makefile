#  ----------------------------------------------------------------------------
#   This confidential and proprietary software may be used only as authorised
#   by a licensing agreement from Arm Limited.
#       (C) COPYRIGHT 2011-2020 Arm Limited, ALL RIGHTS RESERVED
#   The entire notice above must be reproduced on all authorised copies and
#   copies may only be made to the extent permitted by a licensing agreement
#   from Arm Limited.
#  ----------------------------------------------------------------------------

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
    astc_platform_dependents.cpp \
    astc_quantization.cpp \
    astc_stb_tga.cpp \
    astc_symbolic_physical.cpp \
    astc_toplevel.cpp \
    astc_toplevel_help.cpp \
    astc_weight_align.cpp \
    astc_weight_quant_xfer_tables.cpp

HEADERS = \
    astc_codec_internals.h \
    astc_mathlib.h \
    stb_image.h

OBJECTS = $(SOURCES:.cpp=.o)

CPPFLAGS = -std=c++14 -O3 -mfpmath=sse \
           -Wall -Wextra -Wpedantic -Werror -Werror=shadow -Wdouble-promotion

# Validate that the VEC parameter is a supported value, and patch CPPFLAGS
ifeq ($(VEC),nointrin)
CPPFLAGS +=  -msse2 -DASTC_SSE=0 -DASTC_AVX=0 -DASTC_POPCNT=0
else
ifeq ($(VEC),sse2)
CPPFLAGS += -msse2 -DASTC_SSE=20 -DASTC_AVX=0 -DASTC_POPCNT=0
else
ifeq ($(VEC),sse4.2)
CPPFLAGS += -msse4.2 -mpopcnt -DASTC_SSE=42 -DASTC_AVX=0 -DASTC_POPCNT=0
else
ifeq ($(VEC),avx2)
CPPFLAGS += -mavx2 -mpopcnt -DASTC_SSE=42 -DASTC_AVX=2 -DASTC_POPCNT=1
else
$(error Unsupported VEC target, use VEC=nointrin/sse2/sse4.2/avx2)
endif
endif
endif
endif

astcenc: $(OBJECTS)
	@g++ -o $@ $^ $(CPPFLAGS) -lpthread
	@echo "[Link] $@ (using $(VEC))"

$(OBJECTS): %.o: %.cpp $(HEADERS)
	@g++ -c -o $@ $< $(CPPFLAGS)
	@echo "[C++] $<"

clean:
	@rm -f *.o
	@rm -f *.obj
	@rm -f astcenc
	@rm -f astcenc.exe
	@echo "[Clean] Removing build outputs"
