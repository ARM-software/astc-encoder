#  SPDX-License-Identifier: Apache-2.0
#  ----------------------------------------------------------------------------
#  Copyright 2020-2026 Arm Limited
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

set(ASTCENC_TARGET astc${ASTCENC_CODEC}-${ASTCENC_ISA_SIMD})

project(${ASTCENC_TARGET})

include(cmake_compiler.cmake)

add_library(${ASTCENC_TARGET}-static
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
        astcenc_quantization.cpp
        astcenc_symbolic_physical.cpp
        astcenc_weight_align.cpp
        astcenc_weight_quant_xfer_tables.cpp
        astcenc_lz_optimize.cpp)

target_include_directories(${ASTCENC_TARGET}-static
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:.>)

if(${ASTCENC_SHAREDLIB})
    add_library(${ASTCENC_TARGET}-shared
        SHARED
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
            astcenc_quantization.cpp
            astcenc_symbolic_physical.cpp
            astcenc_weight_align.cpp
            astcenc_weight_quant_xfer_tables.cpp
            astcenc_lz_optimize.cpp)

    target_include_directories(${ASTCENC_TARGET}-shared
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
            $<INSTALL_INTERFACE:.>)
endif()

if(${ASTCENC_CLI})
    # Veneer is compiled without any extended ISA so we can safely do
    # ISA compatability checks without triggering a SIGILL
    add_library(${ASTCENC_TARGET}-veneer1
        astcenccli_entry.cpp)

    # Veneer is compiled with extended ISA but without vector length overrides
    # so we can safely do SVE vector length compatibility checks
    add_library(${ASTCENC_TARGET}-veneer2
        astcenccli_entry2.cpp)

    add_executable(${ASTCENC_TARGET}
        astcenccli_compression_metrics.cpp
        astcenccli_error_metrics.cpp
        astcenccli_image.cpp
        astcenccli_image_external.cpp
        astcenccli_image_load_store.cpp
        astcenccli_platform_dependents.cpp
        astcenccli_toplevel.cpp
        astcenccli_toplevel_help.cpp)

    set_property(TARGET ${ASTCENC_TARGET}
        PROPERTY
            OUTPUT_NAME "${ASTCENC_TARGET}${ASTCENC_CLI_POSTFIX}")

    target_link_libraries(${ASTCENC_TARGET}
        PRIVATE
            ${ASTCENC_TARGET}-veneer1
            ${ASTCENC_TARGET}-veneer2
            ${ASTCENC_TARGET}-static)

    if(ASTCENC_RDO_SIZE)
        target_compile_definitions(${ASTCENC_TARGET}
            PRIVATE
                ASTCENC_RDO_SIZE)

        target_link_libraries(${ASTCENC_TARGET}
            PRIVATE
                ZLIB::ZLIB)
    endif()
endif()

macro(astcenc_set_properties ASTCENC_TARGET_NAME ASTCENC_VENEER_TYPE)

    if(${ASTCENC_DECOMPRESSOR})
        target_compile_definitions(${ASTCENC_TARGET_NAME}
            PRIVATE
                ASTCENC_DECOMPRESS_ONLY)
    endif()

    if(${ASTCENC_BLOCK_MAX_TEXELS})
        target_compile_definitions(${ASTCENC_TARGET_NAME}
            PRIVATE
                ASTCENC_BLOCK_MAX_TEXELS=${ASTCENC_BLOCK_MAX_TEXELS})
    endif()

    if(${ASTCENC_DIAGNOSTICS})
        target_compile_definitions(${ASTCENC_TARGET_NAME}
            PUBLIC
                ASTCENC_DIAGNOSTICS)
    endif()

    astcenc_set_compiler_options(${ASTCENC_TARGET_NAME} ${ASTCENC_ISA_SIMD} ${ASTCENC_VENEER_TYPE})

    if(${ASTCENC_CLI})
        # Enable LTO on release builds
        set_property(TARGET ${ASTCENC_TARGET_NAME}
            PROPERTY
                INTERPROCEDURAL_OPTIMIZATION_RELEASE True)

        # Use a static runtime on MSVC builds (ignored on non-MSVC compilers)
        set_property(TARGET ${ASTCENC_TARGET_NAME}
            PROPERTY
                MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif()

endmacro()

string(CONCAT EXTERNAL_CXX_FLAGS
       " $<${is_gnu_fe}: -fno-strict-aliasing>"
       " $<${is_gnu_fe}: -Wno-pedantic>"
       " $<${is_gnu_fe}: -Wno-unused-function>"
       " $<${is_clangcl}: -Xclang -Wno-unused-function>"
       " $<${is_gnu_fe}: -Wno-unused-parameter>"
       " $<${is_gnu_fe}: -Wno-old-style-cast>"
       " $<${is_gnu_fe}: -Wno-double-promotion>"
       " $<${is_gnu_fe}: -Wno-calloc-transposed-args>"
       " $<${is_gnu_fe}: -Wno-zero-as-null-pointer-constant>"
       " $<${is_gnu_fe}: -Wno-disabled-macro-expansion>"
       " $<${is_gnu_fe}: -Wno-reserved-id-macro>"
       " $<${is_gnu_fe}: -Wno-extra-semi-stmt>"
       " $<${is_gnu_fe}: -Wno-implicit-fallthrough>"
       " $<${is_gnu_fe}: -Wno-tautological-type-limit-compare>"
       " $<${is_gnu_fe}: -Wno-cast-qual>"
       " $<${is_gnu_fe}: -Wno-reserved-identifier>"
       " $<${is_clang}: -Wno-missing-prototypes>"
       " $<${is_gnu_fe}: -Wno-missing-field-initializers>"
       " $<${is_gnu_fe}: -Wno-suggest-override>"
       " $<${is_gnu_fe}: -Wno-used-but-marked-unused>"
       " $<${is_gnu_fe}: -Wno-noexcept-type>"
       " $<${is_gnu_fe}: -Wno-comma>"
       " $<${is_gnu_fe}: -Wno-c99-extensions>")

set_source_files_properties(astcenccli_image_external.cpp
    PROPERTIES
        COMPILE_FLAGS ${EXTERNAL_CXX_FLAGS})

astcenc_set_properties(${ASTCENC_TARGET}-static OFF)

# The monkey test requires compression support and excludes diagnostic tracing
if(ASTCENC_UNITTEST AND NOT ASTCENC_DECOMPRESSOR AND NOT ASTCENC_DIAGNOSTICS)
    set(ASTCENC_MONKEY_TARGET test-monkey-${ASTCENC_ISA_SIMD})
    add_executable(${ASTCENC_MONKEY_TARGET} UnitTest/astc_test_threading.cpp)
    target_link_libraries(${ASTCENC_MONKEY_TARGET} PRIVATE ${ASTCENC_TARGET}-static)
    astcenc_set_properties(${ASTCENC_MONKEY_TARGET} OFF)
endif()

target_compile_options(${ASTCENC_TARGET}-static
    PRIVATE
        $<${is_msvc_fe}:/W4>)

if(${ASTCENC_SHAREDLIB})
    astcenc_set_properties(${ASTCENC_TARGET}-shared OFF)

    target_compile_definitions(${ASTCENC_TARGET}-shared
        PRIVATE
            ASTCENC_DYNAMIC_LIBRARY=1)

    target_compile_options(${ASTCENC_TARGET}-shared
        PRIVATE
            $<${is_gnu_fe}:-fvisibility=hidden>
            $<${is_msvc_fe}:/W4>)

    if(NOT ${ASTCENC_UNIVERSAL_BUILD})
        install(TARGETS ${ASTCENC_TARGET}-shared)
    endif()
endif()

if(${ASTCENC_CLI})
    astcenc_set_properties(${ASTCENC_TARGET}-veneer1 1)
    astcenc_set_properties(${ASTCENC_TARGET}-veneer2 2)
    astcenc_set_properties(${ASTCENC_TARGET} 0)

    target_compile_options(${ASTCENC_TARGET}-veneer1
        PRIVATE
            $<${is_msvc_fe}:/W3>)

    target_compile_options(${ASTCENC_TARGET}-veneer2
        PRIVATE
            $<${is_msvc_fe}:/W3>)

    target_compile_options(${ASTCENC_TARGET}
        PRIVATE
            $<${is_msvc_fe}:/W3>)

    string(TIMESTAMP astcencoder_YEAR "%Y")

    configure_file(
        astcenccli_version.h.in
        astcenccli_version.h
        ESCAPE_QUOTES @ONLY)

    target_include_directories(${ASTCENC_TARGET}
        PRIVATE
            ${CMAKE_CURRENT_BINARY_DIR})

    if(NOT ${ASTCENC_UNIVERSAL_BUILD})
        install(TARGETS ${ASTCENC_TARGET})
    endif()
endif()
