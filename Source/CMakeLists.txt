#  SPDX-License-Identifier: Apache-2.0
#  ----------------------------------------------------------------------------
#  Copyright 2020 Arm Limited
#
#  Licensed under the Apache License, Version 2.0 (the "License"); you may not
#  use this file except in compliance with the License. You may obtain a copy
#  of the License at:
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
#  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
#  License for the specific language governing permissions and limitations
#  under the License.
#  ----------------------------------------------------------------------------

project(astcenc)

add_executable(astcenc-${ISA_SIMD})

target_sources(astcenc-${ISA_SIMD}
    PRIVATE
        astcenc_averages_and_directions.cpp
        astcenc_block_sizes2.cpp
        astcenc_color_quantize.cpp
        astcenc_color_unquantize.cpp
        astcenc_compress_symbolic.cpp
        astcenc_compute_variance.cpp
        astcenc_decompress_symbolic.cpp
        astcenc_encoding_choice_error.cpp
        astcenc_entry.cpp
        astcenc_find_best_partitioning.cpp
        astcenc_ideal_endpoints_and_weights.cpp
        astcenc_image.cpp
        astcenc_integer_sequence.cpp
        astcenc_kmeans_partitioning.cpp
        astcenc_mathlib.cpp
        astcenc_mathlib_softfloat.cpp
        astcenc_partition_tables.cpp
        astcenc_percentile_tables.cpp
        astcenc_pick_best_endpoint_format.cpp
        astcenc_platform_isa_detection.cpp
        astcenc_quantization.cpp
        astcenc_symbolic_physical.cpp
        astcenc_weight_align.cpp
        astcenc_weight_quant_xfer_tables.cpp
        astcenccli_error_metrics.cpp
        astcenccli_image.cpp
        astcenccli_image_external.cpp
        astcenccli_image_load_store.cpp
        astcenccli_platform_dependents.cpp
        astcenccli_toplevel.cpp
        astcenccli_toplevel_help.cpp)

target_compile_features(astcenc-${ISA_SIMD}
    PRIVATE
        cxx_std_14)

target_compile_definitions(astcenc-${ISA_SIMD}
    PRIVATE
        # Common defines
        ASTCENC_ISA_INVARIANCE=$<BOOL:${ISA_INVARIANCE}>
        # MSVC defines
        $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>)

target_compile_options(astcenc-${ISA_SIMD}
    PRIVATE
        # MSVC defines
        $<$<CXX_COMPILER_ID:MSVC>:/EHsc>
        # G++ and Clang++ defines
        $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall>
        $<$<CXX_COMPILER_ID:GNU,Clang>:-Wextra>
        $<$<CXX_COMPILER_ID:GNU,Clang>:-Wpedantic>
        $<$<CXX_COMPILER_ID:GNU,Clang>:-Werror>
        $<$<CXX_COMPILER_ID:GNU,Clang>:-Wshadow>
        $<$<CXX_COMPILER_ID:GNU,Clang>:-Wdouble-promotion>)

target_include_directories(astcenc-${ISA_SIMD}
    PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(astcenc-${ISA_SIMD}
    PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:pthread>)

# Enable LTO on release builds
set_property(TARGET astcenc-${ISA_SIMD} PROPERTY
             INTERPROCEDURAL_OPTIMIZATION_RELEASE True)

# Use a static runtime on MSVC builds (ignored on non-MSVC compilers)
set_property(TARGET astcenc-${ISA_SIMD} PROPERTY
             MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")

    string(CONCAT EXTERNAL_CXX_FLAGS "-Wno-unused-parameter"
                                     " -Wno-double-promotion"
                                     " -fno-strict-aliasing")

    set_source_files_properties(astcenccli_image_external.cpp
        PROPERTIES COMPILE_FLAGS ${EXTERNAL_CXX_FLAGS})
endif()

# Set up configuration for SIMD ISA builds
if(${ISA_SIMD} MATCHES "none")
    target_compile_definitions(astcenc-${ISA_SIMD}
        PRIVATE
            ASTCENC_SSE=0
            ASTCENC_AVX=0
            ASTCENC_POPCNT=0
            ASTCENC_VECALIGN=16)

    target_compile_options(astcenc-${ISA_SIMD}
        PRIVATE
            $<$<CXX_COMPILER_ID:GNU,Clang>:-mfpmath=sse -msse2>)

elseif(${ISA_SIMD} MATCHES "sse2")
    target_compile_definitions(astcenc-${ISA_SIMD}
        PRIVATE
            ASTCENC_SSE=20
            ASTCENC_AVX=0
            ASTCENC_POPCNT=0
            ASTCENC_VECALIGN=16)

    target_compile_options(astcenc-${ISA_SIMD}
        PRIVATE
            $<$<CXX_COMPILER_ID:GNU,Clang>:-mfpmath=sse -msse2>)

elseif(${ISA_SIMD} MATCHES "sse4.1")
    target_compile_definitions(astcenc-${ISA_SIMD}
        PRIVATE
            ASTCENC_SSE=41
            ASTCENC_AVX=0
            ASTCENC_POPCNT=1
            ASTCENC_VECALIGN=16)

    target_compile_options(astcenc-${ISA_SIMD}
        PRIVATE
            $<$<CXX_COMPILER_ID:GNU,Clang>:-mfpmath=sse -msse4.1 -mpopcnt>)

elseif(${ISA_SIMD} MATCHES "avx2")
    target_compile_definitions(astcenc-${ISA_SIMD}
        PRIVATE
            ASTCENC_SSE=41
            ASTCENC_AVX=2
            ASTCENC_POPCNT=1
            ASTCENC_VECALIGN=32)

    target_compile_options(astcenc-${ISA_SIMD}
        PRIVATE
            $<$<CXX_COMPILER_ID:GNU,Clang>:-mfpmath=sse -mavx2 -mpopcnt>
            $<$<CXX_COMPILER_ID:MSVC>:/arch:AVX2>)
endif()

install(TARGETS astcenc-${ISA_SIMD} DESTINATION ${PACKAGE_ROOT})
