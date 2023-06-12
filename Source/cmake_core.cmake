#  SPDX-License-Identifier: Apache-2.0
#  ----------------------------------------------------------------------------
#  Copyright 2020-2023 Arm Limited
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

# On CMake 3.25 or older CXX_COMPILER_FRONTEND_VARIANT is not always set
if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "")
    set(CMAKE_CXX_COMPILER_FRONTEND_VARIANT "${CMAKE_CXX_COMPILER_ID}")
endif()

# Compiler accepts MSVC-style command line options
set(is_msvc_fe "$<STREQUAL:${CMAKE_CXX_COMPILER_FRONTEND_VARIANT},MSVC>")
# Compiler accepts GNU-style command line options
set(is_gnu_fe1 "$<STREQUAL:${CMAKE_CXX_COMPILER_FRONTEND_VARIANT},GNU>")
# Compiler accepts AppleClang-style command line options, which is also GNU-style
set(is_gnu_fe2 "$<STREQUAL:${CMAKE_CXX_COMPILER_FRONTEND_VARIANT},AppleClang>")
# Compiler accepts GNU-style command line options
set(is_gnu_fe "$<OR:${is_gnu_fe1},${is_gnu_fe2}>")

# Compiler is Visual Studio cl.exe
set(is_msvccl "$<AND:${is_msvc_fe},$<CXX_COMPILER_ID:MSVC>>")
# Compiler is Visual Studio clangcl.exe
set(is_clangcl "$<AND:${is_msvc_fe},$<CXX_COMPILER_ID:Clang>>")
# Compiler is upstream clang with the standard frontend
set(is_clang "$<AND:${is_gnu_fe},$<CXX_COMPILER_ID:Clang,AppleClang>>")

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
        astcenc_weight_quant_xfer_tables.cpp)

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
            astcenc_weight_quant_xfer_tables.cpp)

    target_include_directories(${ASTCENC_TARGET}-shared
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
            $<INSTALL_INTERFACE:.>)
endif()

if(${ASTCENC_CLI})
    # Veneer is compiled without any extended ISA so we can safely do
    # ISA compatability checks without triggering a SIGILL
    add_library(${ASTCENC_TARGET}-veneer
        astcenccli_entry.cpp)

    add_executable(${ASTCENC_TARGET}
        astcenccli_error_metrics.cpp
        astcenccli_image.cpp
        astcenccli_image_external.cpp
        astcenccli_image_load_store.cpp
        astcenccli_platform_dependents.cpp
        astcenccli_toplevel.cpp
        astcenccli_toplevel_help.cpp)

    target_link_libraries(${ASTCENC_TARGET}
        PRIVATE
            ${ASTCENC_TARGET}-veneer
            ${ASTCENC_TARGET}-static)
endif()

macro(astcenc_set_properties ASTCENC_TARGET_NAME ASTCENC_IS_VENEER)

    target_compile_features(${ASTCENC_TARGET_NAME}
        PRIVATE
            cxx_std_14)

    target_compile_definitions(${ASTCENC_TARGET_NAME}
        PRIVATE
            $<${is_msvc_fe}:_CRT_SECURE_NO_WARNINGS>)

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

    target_compile_options(${ASTCENC_TARGET_NAME}
        PRIVATE
            # Use pthreads on Linux/macOS
            $<$<PLATFORM_ID:Linux,Darwin>:-pthread>

            # MSVC compiler defines
            $<${is_msvc_fe}:/EHsc>
            $<${is_msvccl}:/wd4324>

            # G++ and Clang++ compiler defines
            $<${is_gnu_fe}:-Wall>
            $<${is_gnu_fe}:-Wextra>
            $<${is_gnu_fe}:-Wpedantic>
            $<${is_gnu_fe}:-Werror>
            $<${is_gnu_fe}:-Wshadow>
            $<${is_gnu_fe}:-Wdouble-promotion>
            $<${is_clang}:-Wdocumentation>

            # Hide noise thrown up by Clang 10 and clang-cl
            $<${is_gnu_fe}:-Wno-unknown-warning-option>
            $<${is_gnu_fe}:-Wno-c++98-compat-pedantic>
            $<${is_gnu_fe}:-Wno-c++98-c++11-compat-pedantic>
            $<${is_gnu_fe}:-Wno-float-equal>
            $<${is_gnu_fe}:-Wno-deprecated-declarations>
            $<${is_gnu_fe}:-Wno-atomic-implicit-seq-cst>

            # Clang 10 also throws up warnings we need to investigate (ours)
            $<${is_gnu_fe}:-Wno-cast-align>
            $<${is_gnu_fe}:-Wno-sign-conversion>
            $<${is_gnu_fe}:-Wno-implicit-int-conversion>
            $<${is_gnu_fe}:-Wno-shift-sign-overflow>
            $<${is_gnu_fe}:-Wno-format-nonliteral>
            $<${is_gnu_fe}:-Wno-reserved-identifier>
            $<${is_gnu_fe}:-Wno-cast-function-type>

            # Force DWARF4 for Valgrind profiling
            $<$<AND:$<PLATFORM_ID:Linux,Darwin>,${is_clang}>:-gdwarf-4>

            # Disable non-portable Windows.h warning (fixing it fails builds on MinGW)
            $<$<AND:$<PLATFORM_ID:Windows>,${is_clang}>:-Wno-nonportable-system-include-path>)

    target_link_options(${ASTCENC_TARGET_NAME}
        PRIVATE
            # Use pthreads on Linux/macOS
            $<$<PLATFORM_ID:Linux,Darwin>:-pthread>)

    if(${ASTCENC_ASAN})
        target_compile_options(${ASTCENC_TARGET_NAME}
            PRIVATE
                $<${is_clang}:-fsanitize=address>)

        target_link_options(${ASTCENC_TARGET_NAME}
            PRIVATE
                $<${is_clang}:-fsanitize=address>)
    endif()

    if(NOT ${ASTCENC_INVARIANCE})
        target_compile_definitions(${ASTCENC_TARGET_NAME}
            PRIVATE
                ASTCENC_NO_INVARIANCE=1)

        # For Visual Studio prior to 2022 (compiler < 19.30) /fp:precise
        # For Visual Studio 2022 (compiler >= 19.30) /fp:precise and /fp:contract

        # For Visual Studio 2022 ClangCL seems to have accidentally enabled contraction by default,
        # so behaves differently to CL.exe. Use the -Xclang argument to workaround and allow access
        # GNU-style switch to control contraction on the assumption this gets fixed and disabled.
        # Note ClangCL does not accept /fp:contract as an argument as of v15.0.7.
        target_compile_options(${ASTCENC_TARGET_NAME}
            PRIVATE
                $<${is_msvccl}:/fp:precise>
                $<${is_clangcl}:/fp:precise>
                $<$<AND:${is_msvccl},$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,19.30>>:/fp:contract>
                $<$<AND:${is_clangcl},$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,14.0.0>>:-Xclang -ffp-contract=fast>
                $<$<AND:${is_clang},$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,10.0.0>>:-ffp-model=precise>
                $<${is_gnu_fe}:-ffp-contract=fast>)
    else()
        # For Visual Studio prior to 2022 (compiler < 19.30) /fp:strict
        # For Visual Studio 2022 (compiler >= 19.30) /fp:precise

        # For Visual Studio 2022 ClangCL seems to have accidentally enabled contraction by default,
        # so behaves differently to CL.exe. Use the -Xclang argument to workaround and allow access
        # GNU-style switch to control contraction and force disable.
        target_compile_options(${ASTCENC_TARGET_NAME}
            PRIVATE
                $<$<AND:${is_msvccl},$<VERSION_LESS:$<CXX_COMPILER_VERSION>,19.30>>:/fp:strict>
                $<$<AND:${is_msvccl},$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,19.30>>:/fp:precise>
                $<${is_clangcl}:/fp:precise>
                $<$<AND:${is_clangcl},$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,14.0.0>>:-Xclang -ffp-contract=off>
                $<$<AND:${is_clang},$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,10.0.0>>:-ffp-model=precise>
                $<${is_gnu_fe}:-ffp-contract=off>)
    endif()

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

    # Set up configuration for SIMD ISA builds
    if(${ASTCENC_ISA_SIMD} MATCHES "none")
        target_compile_definitions(${ASTCENC_TARGET_NAME}
            PRIVATE
                ASTCENC_NEON=0
                ASTCENC_SSE=0
                ASTCENC_AVX=0
                ASTCENC_POPCNT=0
                ASTCENC_F16C=0)

    elseif(${ASTCENC_ISA_SIMD} MATCHES "neon")
        target_compile_definitions(${ASTCENC_TARGET_NAME}
            PRIVATE
                ASTCENC_NEON=1
                ASTCENC_SSE=0
                ASTCENC_AVX=0
                ASTCENC_POPCNT=0
                ASTCENC_F16C=0)

        # Workaround MSVC codegen bug for NEON builds on VS 2022 17.2 or older
        # https://developercommunity.visualstudio.com/t/inlining-turns-constant-into-register-operand-for/1394798
        if((CMAKE_CXX_COMPILER_ID MATCHES "MSVC") AND (MSVC_VERSION LESS 1933))
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    $<${is_msvccl}:/d2ssa-cfg-sink->)
        endif()

    elseif(${ASTCENC_ISA_SIMD} MATCHES "sse2")
        target_compile_definitions(${ASTCENC_TARGET_NAME}
            PRIVATE
                ASTCENC_NEON=0
                ASTCENC_SSE=20
                ASTCENC_AVX=0
                ASTCENC_POPCNT=0
                ASTCENC_F16C=0)

        # Force SSE2 on AppleClang (normally SSE4.1 is the default)
        target_compile_options(${ASTCENC_TARGET_NAME}
            PRIVATE
                $<${is_clangcl}:-msse2>
                $<${is_gnu_fe}:-msse2>
                $<${is_gnu_fe}:-mno-sse4.1>
                $<${is_gnu_fe}:-Wno-unused-command-line-argument>)

    elseif(${ASTCENC_ISA_SIMD} MATCHES "sse4.1")
        target_compile_definitions(${ASTCENC_TARGET_NAME}
            PRIVATE
                ASTCENC_NEON=0
                ASTCENC_SSE=41
                ASTCENC_AVX=0
                ASTCENC_POPCNT=1
                ASTCENC_F16C=0)

        if (${ASTCENC_IS_VENEER})
            # Force SSE2 on AppleClang (normally SSE4.1 is the default)
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    $<${is_gnu_fe}:-msse2>
                    $<${is_gnu_fe}:-mno-sse4.1>
                    $<${is_gnu_fe}:-Wno-unused-command-line-argument>)
        else()
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    $<${is_clangcl}:-msse4.1 -mpopcnt>
                    $<${is_gnu_fe}:-msse4.1 -mpopcnt>
                    $<${is_gnu_fe}:-Wno-unused-command-line-argument>)
        endif()

    elseif(${ASTCENC_ISA_SIMD} MATCHES "avx2")
        target_compile_definitions(${ASTCENC_TARGET_NAME}
            PRIVATE
                ASTCENC_NEON=0
                ASTCENC_SSE=41
                ASTCENC_AVX=2
                ASTCENC_POPCNT=1
                ASTCENC_F16C=1)

        if (${ASTCENC_IS_VENEER})
            # Force SSE2 on AppleClang (normally SSE4.1 is the default)
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    $<${is_gnu_fe}:-msse2>
                    $<${is_gnu_fe}:-mno-sse4.1>
                    $<${is_gnu_fe}:-Wno-unused-command-line-argument>)
        else()
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    $<${is_msvc_fe}:/arch:AVX2>
                    $<${is_clangcl}:-mavx2 -mpopcnt -mf16c>
                    $<${is_gnu_fe}:-mavx2 -mpopcnt -mf16c>
                    $<${is_gnu_fe}:-Wno-unused-command-line-argument>)
        endif()

        # Non-invariant builds enable us to loosen the compiler constraints on
        # floating point, but this is only worth doing on CPUs with AVX2 because
        # this implies we can also enable the FMA instruction set extensions
        # which significantly improve performance. Note that this DOES reduce
        # image quality by up to 0.2 dB (normally much less), but buys an
        # average of 10-15% performance improvement ...
        if((NOT ${ASTCENC_INVARIANCE}) AND (NOT ${ASTCENC_IS_VENEER}))
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    $<${is_gnu_fe}:-mfma>)
        endif()

    endif()

endmacro()

string(CONCAT EXTERNAL_CXX_FLAGS
       " $<${is_gnu_fe}: -fno-strict-aliasing>"
       " $<${is_gnu_fe}: -Wno-unused-parameter>"
       " $<${is_gnu_fe}: -Wno-old-style-cast>"
       " $<${is_gnu_fe}: -Wno-double-promotion>"
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
    astcenc_set_properties(${ASTCENC_TARGET}-veneer ON)
    astcenc_set_properties(${ASTCENC_TARGET} OFF)

    target_compile_options(${ASTCENC_TARGET}
        PRIVATE
            $<${is_msvc_fe}:/W3>)

    target_compile_options(${ASTCENC_TARGET}-veneer
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
