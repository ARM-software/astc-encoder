#  SPDX-License-Identifier: Apache-2.0
#  ----------------------------------------------------------------------------
#  Copyright 2020-2021 Arm Limited
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

if(${UNIVERSAL_BUILD})
    set(ASTC_TARGET astc${CODEC})
else()
    set(ASTC_TARGET astc${CODEC}-${ISA_SIMD})
endif()

project(${ASTC_TARGET})

set(GNU_LIKE "GNU,Clang,AppleClang")
set(CLANG_LIKE "Clang,AppleClang")

add_library(${ASTC_TARGET}-static
    STATIC
        astcenc_averages_and_directions.cpp
        astcenc_block_sizes.cpp
        astcenc_color_quantize.cpp
        astcenc_color_unquantize.cpp
        astcenc_compress_symbolic.cpp
        astcenc_compute_variance.cpp
        astcenc_decompress_symbolic.cpp
        astcenc_diagnostic_trace.cpp
        astcenc_entry.cpp
        astcenc_find_best_partitioning.cpp
        astcenc_ideal_endpoints_and_weights.cpp
        astcenc_image.cpp
        astcenc_integer_sequence.cpp
        astcenc_mathlib.cpp
        astcenc_mathlib_softfloat.cpp
        astcenc_partition_tables.cpp
        astcenc_percentile_tables.cpp
        astcenc_pick_best_endpoint_format.cpp
        astcenc_platform_isa_detection.cpp
        astcenc_quantization.cpp
        astcenc_symbolic_physical.cpp
        astcenc_weight_align.cpp
        astcenc_weight_quant_xfer_tables.cpp)

target_include_directories(${ASTC_TARGET}-static
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:.>)

if(${CLI})
    add_executable(${ASTC_TARGET}
        astcenccli_error_metrics.cpp
        astcenccli_image.cpp
        astcenccli_image_external.cpp
        astcenccli_image_load_store.cpp
        astcenccli_platform_dependents.cpp
        astcenccli_toplevel.cpp
        astcenccli_toplevel_help.cpp)

    target_link_libraries(${ASTC_TARGET}
        PRIVATE
            ${ASTC_TARGET}-static)
endif()

macro(astcenc_set_properties NAME)

    target_compile_features(${NAME}
        PRIVATE
            cxx_std_14)

    target_compile_definitions(${NAME}
        PRIVATE
            # MSVC defines
            $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>)

    if(${DECOMPRESSOR})
        target_compile_definitions(${NAME}
            PRIVATE
                ASTCENC_DECOMPRESS_ONLY)
    endif()

    if(${DIAGNOSTICS})
        target_compile_definitions(${NAME}
            PUBLIC
                ASTCENC_DIAGNOSTICS)
    endif()

    target_compile_options(${NAME}
        PRIVATE
            # Use pthreads on Linux/macOS
            $<$<PLATFORM_ID:Linux,Darwin>:-pthread>

            # MSVC compiler defines
            $<$<CXX_COMPILER_ID:MSVC>:/EHsc>
            $<$<CXX_COMPILER_ID:MSVC>:/fp:strict>

            # G++ and Clang++ compiler defines
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wextra>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wpedantic>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Werror>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wshadow>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wdouble-promotion>

            # Hide noise thrown up by Clang 10 and clang-cl
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wno-unknown-warning-option>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wno-c++98-compat-pedantic>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wno-c++98-c++11-compat-pedantic>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wno-float-equal>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wno-deprecated-declarations>

            # Clang 10 also throws up warnings we need to investigate (ours)
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wno-old-style-cast>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wno-cast-align>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wno-sign-conversion>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wno-implicit-int-conversion>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wno-shift-sign-overflow>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wno-format-nonliteral>

            $<$<CXX_COMPILER_ID:Clang>:-Wdocumentation>)


    target_link_options(${NAME}
        PRIVATE
            # Use pthreads on Linux/macOS
            $<$<PLATFORM_ID:Linux,Darwin>:-pthread>)

    if(${CLI})
        # Enable LTO on release builds
        set_property(TARGET ${NAME}
            PROPERTY
                INTERPROCEDURAL_OPTIMIZATION_RELEASE True)

        # Use a static runtime on MSVC builds (ignored on non-MSVC compilers)
        set_property(TARGET ${NAME}
            PROPERTY
                MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif()

    # Set up configuration for SIMD ISA builds
    if(${ISA_SIMD} MATCHES "none")
        if(NOT ${UNIVERSAL_BUILD})
            target_compile_definitions(${NAME}
                PRIVATE
                    ASTCENC_NEON=0
                    ASTCENC_SSE=0
                    ASTCENC_AVX=0
                    ASTCENC_POPCNT=0
                    ASTCENC_F16C=0)
        endif()

    elseif(${ISA_SIMD} MATCHES "neon")
        if(NOT ${UNIVERSAL_BUILD})
            target_compile_definitions(${NAME}
                PRIVATE
                    ASTCENC_NEON=1
                    ASTCENC_SSE=0
                    ASTCENC_AVX=0
                    ASTCENC_POPCNT=0
                    ASTCENC_F16C=0)
        endif()

    elseif((${ISA_SIMD} MATCHES "sse2") OR (${UNIVERSAL_BUILD} AND ${ISA_SSE2}))
        if(NOT ${UNIVERSAL_BUILD})
            target_compile_definitions(${NAME}
                PRIVATE
                    ASTCENC_NEON=0
                    ASTCENC_SSE=20
                    ASTCENC_AVX=0
                    ASTCENC_POPCNT=0
                    ASTCENC_F16C=0)
        endif()

        # These settings are needed on AppleClang as SSE4.1 is on by default
        # Suppress unused argument for macOS universal build behavior
        target_compile_options(${NAME}
            PRIVATE
                $<$<CXX_COMPILER_ID:AppleClang>:-msse2>
                $<$<CXX_COMPILER_ID:AppleClang>:-mno-sse4.1>
                $<$<CXX_COMPILER_ID:AppleClang>:-Wno-unused-command-line-argument>)

    elseif((${ISA_SIMD} MATCHES "sse4.1") OR (${UNIVERSAL_BUILD} AND ${ISA_SSE41}))
        if(NOT ${UNIVERSAL_BUILD})
            target_compile_definitions(${NAME}
                PRIVATE
                    ASTCENC_NEON=0
                    ASTCENC_SSE=41
                    ASTCENC_AVX=0
                    ASTCENC_POPCNT=1
                    ASTCENC_F16C=0)
        endif()

        # Suppress unused argument for macOS universal build behavior
        target_compile_options(${NAME}
            PRIVATE
                $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-msse4.1 -mpopcnt>
                $<$<CXX_COMPILER_ID:AppleClang>:-Wno-unused-command-line-argument>)

    elseif((${ISA_SIMD} MATCHES "avx2") OR (${UNIVERSAL_BUILD} AND ${ISA_AVX2}))
        if(NOT ${UNIVERSAL_BUILD})
            target_compile_definitions(${NAME}
                PRIVATE
                    ASTCENC_NEON=0
                    ASTCENC_SSE=41
                    ASTCENC_AVX=2
                    ASTCENC_POPCNT=1
                    ASTCENC_F16C=1)
        endif()

        # Suppress unused argument for macOS universal build behavior
        target_compile_options(${NAME}
            PRIVATE
                $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-mavx2 -mpopcnt -mf16c>
                $<$<CXX_COMPILER_ID:MSVC>:/arch:AVX2>
                $<$<CXX_COMPILER_ID:AppleClang>:-Wno-unused-command-line-argument>)

    endif()

endmacro()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    string(CONCAT EXTERNAL_CXX_FLAGS
            " $<$<NOT:$<CXX_COMPILER_ID:MSVC>>: -fno-strict-aliasing>"
            " $<$<NOT:$<CXX_COMPILER_ID:MSVC>>: -Wno-unused-parameter>"
            " $<$<NOT:$<CXX_COMPILER_ID:MSVC>>: -Wno-double-promotion>"
            " $<$<NOT:$<CXX_COMPILER_ID:MSVC>>: -Wno-zero-as-null-pointer-constant>"
            " $<$<NOT:$<CXX_COMPILER_ID:MSVC>>: -Wno-disabled-macro-expansion>"
            " $<$<NOT:$<CXX_COMPILER_ID:MSVC>>: -Wno-reserved-id-macro>"
            " $<$<NOT:$<CXX_COMPILER_ID:MSVC>>: -Wno-extra-semi-stmt>"
            " $<$<NOT:$<CXX_COMPILER_ID:MSVC>>: -Wno-implicit-fallthrough>"
            " $<$<NOT:$<CXX_COMPILER_ID:MSVC>>: -Wno-tautological-type-limit-compare>"
            " $<$<NOT:$<CXX_COMPILER_ID:MSVC>>: -Wno-cast-qual>"
            " $<$<CXX_COMPILER_ID:Clang>: -Wno-missing-prototypes>")

    set_source_files_properties(astcenccli_image_external.cpp
        PROPERTIES
            COMPILE_FLAGS ${EXTERNAL_CXX_FLAGS})
endif()

astcenc_set_properties(${ASTC_TARGET}-static)

if(${CLI})
    astcenc_set_properties(${ASTC_TARGET})

    string(TIMESTAMP astcencoder_YEAR "%Y")

    configure_file(
        astcenccli_version.h.in
        astcenccli_version.h
        ESCAPE_QUOTES @ONLY)

    target_include_directories(${ASTC_TARGET}
        PRIVATE
            ${CMAKE_CURRENT_BINARY_DIR})

    install(TARGETS ${ASTC_TARGET} DESTINATION ${PACKAGE_ROOT})
endif()
