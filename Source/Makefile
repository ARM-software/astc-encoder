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

SOURCES = \
    astc_averages_and_directions.cpp \
    astc_block_sizes2.cpp \
    astc_color_quantize.cpp \
    astc_color_unquantize.cpp \
    astc_compress_symbolic.cpp \
    astc_compute_variance.cpp \
    astc_decompress_symbolic.cpp \
    astc_encoding_choice_error.cpp \
    astc_find_best_partitioning.cpp \
    astc_ideal_endpoints_and_weights.cpp \
    astc_image_load_store.cpp \
    astc_integer_sequence.cpp \
    astc_kmeans_partitioning.cpp \
    astc_ktx_dds.cpp \
    astc_main.cpp \
    astc_partition_tables.cpp \
    astc_percentile_tables.cpp \
    astc_pick_best_endpoint_format.cpp \
    astc_quantization.cpp \
    astc_stb_tga.cpp \
    astc_symbolic_physical.cpp \
    astc_toplevel.cpp \
    astc_weight_align.cpp \
    astc_weight_quant_xfer_tables.cpp \
    mathlib.cpp \
    softfloat.cpp

HEADERS = \
    astc_codec_internals.h \
    mathlib.h \
    softfloat.h \
    stb_image.h

OBJECTS = $(SOURCES:.cpp=.o)

CPPFLAGS = -std=c++11 -O3 -msse2 -mfpmath=sse \
           -Wall -Wextra -Wpedantic -Werror -Werror=shadow

astcenc: $(OBJECTS)
	@g++ -o $@ $^ $(CPPFLAGS) -lpthread
	@echo "[Link] $@"

$(OBJECTS): %.o: %.cpp $(HEADERS)
	@g++ -c -o $@ $< $(CPPFLAGS)
	@echo "[C++] $<"

clean:
	@rm -f *.o
	@rm -f *.obj
	@rm -f astcenc
	@rm -f astcenc.exe
	@echo "[Clean] Removing build outputs"
