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

# @brief Configure compiler options shared by codec and unit-test targets.
#
# @param ASTCENC_TARGET_NAME   The target to configure.
# @param ASTCENC_ISA_SIMD      The target SIMD instruction set.
# @param ASTCENC_VENEER_TYPE   0 for the codec, 1 for the ISA check, or 2 for the width check.
function(astcenc_set_compiler_options ASTCENC_TARGET_NAME ASTCENC_ISA_SIMD ASTCENC_VENEER_TYPE)

    target_compile_definitions(${ASTCENC_TARGET_NAME}
        PRIVATE
            $<${is_msvc_fe}:_CRT_SECURE_NO_WARNINGS>)

    target_compile_options(${ASTCENC_TARGET_NAME}
        PRIVATE
            # Use pthreads on Linux/macOS
            $<$<PLATFORM_ID:Linux,Darwin>:-pthread>

            # MSVC compiler defines
            $<${is_msvc_fe}:/EHsc>
            $<$<AND:$<BOOL:${ASTCENC_WERROR}>,${is_msvc_fe}>:/WX>
            $<${is_msvccl}:/wd4324>

            # G++ and Clang++ compiler defines
            $<${is_gnu_fe}:-Wall>
            $<${is_gnu_fe}:-Wextra>
            $<${is_gnu_fe}:-Wpedantic>
            $<$<AND:$<BOOL:${ASTCENC_WERROR}>,${is_gnu_fe}>:-Werror>
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
            $<${is_clang}:-Wno-overriding-option>

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
                $<${is_gnu_fe}:-fsanitize=address>
                $<${is_gnu_fe}:-fno-sanitize-recover=all>)

        target_link_options(${ASTCENC_TARGET_NAME}
            PRIVATE
                $<${is_gnu_fe}:-fsanitize=address>
                $<${is_clang}:-fuse-ld=lld>)
    endif()

    if(${ASTCENC_TSAN})
        target_compile_options(${ASTCENC_TARGET_NAME}
            PRIVATE
                $<${is_gnu_fe}:-fsanitize=thread>
                $<${is_gnu_fe}:-fno-sanitize-recover=all>)

        target_link_options(${ASTCENC_TARGET_NAME}
            PRIVATE
                $<${is_gnu_fe}:-fsanitize=thread>
                $<${is_clang}:-fuse-ld=lld>)
    endif()

    if(${ASTCENC_UBSAN})
        target_compile_options(${ASTCENC_TARGET_NAME}
            PRIVATE
                $<${is_gnu_fe}:-fsanitize=undefined>
                $<${is_gnu_fe}:-fno-sanitize-recover=all>)

        target_link_options(${ASTCENC_TARGET_NAME}
            PRIVATE
                $<${is_gnu_fe}:-fsanitize=undefined>
                $<${is_clang}:-fuse-ld=lld>)
    endif()

    # Clang's -fno-unsafe-math-optimizations enables FP trapping
    # Disable traps explicitly without relaxing floating-point value semantics
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
                $<$<AND:${is_msvccl},$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,19.30>>:/fp:contract>

                $<${is_clangcl}:/fp:precise>
                $<$<AND:${is_clangcl},$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,14.0.0>>:-Xclang -ffp-contract=fast>

                $<$<AND:${is_clang},$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,10.0.0>>:-fno-unsafe-math-optimizations>
                $<$<AND:${is_clang},$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,10.0.0>>:-fno-trapping-math>

                $<${is_gnu_fe}:-ffp-contract=fast>
                $<${is_gnu_fe}:-fno-math-errno>)
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

                $<$<AND:${is_clang},$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,10.0.0>>:-fno-unsafe-math-optimizations>
                $<$<AND:${is_clang},$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,10.0.0>>:-fno-trapping-math>

                $<${is_gnu_fe}:-ffp-contract=off>
                $<${is_gnu_fe}:-fno-math-errno>)
    endif()

    # Set up configuration for SIMD ISA builds
    if(${ASTCENC_ISA_SIMD} MATCHES "none")
        target_compile_definitions(${ASTCENC_TARGET_NAME}
            PRIVATE
                ASTCENC_NEON=0
                ASTCENC_SVE=0
                ASTCENC_SVE_WIDTH=0
                ASTCENC_SSE=0
                ASTCENC_AVX=0
                ASTCENC_POPCNT=0
                ASTCENC_F16C=0)

        if(${ASTCENC_BIG_ENDIAN})
            target_compile_definitions(${ASTCENC_TARGET_NAME}
                PRIVATE
                    ASTCENC_BIG_ENDIAN=1)
        endif()

    elseif(${ASTCENC_ISA_SIMD} MATCHES "neon")
        target_compile_definitions(${ASTCENC_TARGET_NAME}
            PRIVATE
                ASTCENC_NEON=1
                ASTCENC_SVE=0
                ASTCENC_SVE_WIDTH=0
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

    elseif(${ASTCENC_ISA_SIMD} MATCHES "^sve2_256$")
        target_compile_definitions(${ASTCENC_TARGET_NAME}
            PRIVATE
                ASTCENC_NEON=1
                ASTCENC_SVE=2
                ASTCENC_SVE_WIDTH=8
                ASTCENC_SSE=0
                ASTCENC_AVX=0
                ASTCENC_POPCNT=0
                ASTCENC_F16C=0)

        # Enable SVE2 and BitPerm in the core library
        if (NOT ${ASTCENC_VENEER_TYPE})
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    -march=armv8-a+sve2+sve2-bitperm -msve-vector-bits=256)

        # Enable SVE2 and BitPerm without fixed vector length in the veneer
        elseif (${ASTCENC_VENEER_TYPE} EQUAL 2)
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    -march=armv8-a+sve2+sve2-bitperm)
        endif()

    elseif(${ASTCENC_ISA_SIMD} MATCHES "^sve2_128$")
        target_compile_definitions(${ASTCENC_TARGET_NAME}
            PRIVATE
                ASTCENC_NEON=1
                ASTCENC_SVE=2
                ASTCENC_SVE_WIDTH=4
                ASTCENC_SSE=0
                ASTCENC_AVX=0
                ASTCENC_POPCNT=0
                ASTCENC_F16C=0)

        # Enable SVE2 and BitPerm in the core library
        # Note that for 128-bit SVE the generated code is actually
        # vector-length agnostic, but any manual intrinsics used in the
        # four-wide SVE library use 128-bit data width predicates
        if (NOT ${ASTCENC_VENEER_TYPE})
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    -march=armv8-a+sve2+sve2-bitperm)

        # Enable SVE2 and BitPerm without fixed vector length in the veneer
        elseif (${ASTCENC_VENEER_TYPE} EQUAL 2)
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    -march=armv8-a+sve2+sve2-bitperm)
        endif()

    elseif(${ASTCENC_ISA_SIMD} MATCHES "^sve_256$")
        target_compile_definitions(${ASTCENC_TARGET_NAME}
            PRIVATE
                ASTCENC_NEON=1
                ASTCENC_SVE=1
                ASTCENC_SVE_WIDTH=8
                ASTCENC_SSE=0
                ASTCENC_AVX=0
                ASTCENC_POPCNT=0
                ASTCENC_F16C=0)

        # Enable SVE in the core library
        if (NOT ${ASTCENC_VENEER_TYPE})
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    -march=armv8-a+sve -msve-vector-bits=256)

        # Enable SVE without fixed vector length in the veneer
        elseif (${ASTCENC_VENEER_TYPE} EQUAL 2)
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    -march=armv8-a+sve)
        endif()

    elseif(${ASTCENC_ISA_SIMD} MATCHES "^sve_128$")
        target_compile_definitions(${ASTCENC_TARGET_NAME}
            PRIVATE
                ASTCENC_NEON=1
                ASTCENC_SVE=1
                ASTCENC_SVE_WIDTH=4
                ASTCENC_SSE=0
                ASTCENC_AVX=0
                ASTCENC_POPCNT=0
                ASTCENC_F16C=0)

        # Enable SVE in the core library
        # Note that for 128-bit SVE the generated code is actually
        # vector-length agnostic, but any manual intrinsics used in the
        # four-wide SVE library use 128-bit data width predicates
        if (NOT ${ASTCENC_VENEER_TYPE})
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    -march=armv8-a+sve)

        # Enable SVE without fixed vector length in the veneer
        elseif (${ASTCENC_VENEER_TYPE} EQUAL 2)
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    -march=armv8-a+sve)
        endif()

    elseif(${ASTCENC_ISA_SIMD} MATCHES "sse2")
        target_compile_definitions(${ASTCENC_TARGET_NAME}
            PRIVATE
                ASTCENC_NEON=0
                ASTCENC_SVE=0
                ASTCENC_SVE_WIDTH=0
                ASTCENC_SSE=20
                ASTCENC_AVX=0
                ASTCENC_POPCNT=0
                ASTCENC_F16C=0
                ASTCENC_X86_GATHERS=0)

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
                ASTCENC_SVE=0
                ASTCENC_SVE_WIDTH=0
                ASTCENC_SSE=41
                ASTCENC_AVX=0
                ASTCENC_POPCNT=1
                ASTCENC_F16C=0
                ASTCENC_X86_GATHERS=0)

        if (${ASTCENC_VENEER_TYPE} GREATER 0)
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

    elseif(${ASTCENC_ISA_SIMD} MATCHES "avx512")
        target_compile_definitions(${ASTCENC_TARGET_NAME}
            PRIVATE
                ASTCENC_NEON=0
                ASTCENC_SVE=0
                ASTCENC_SVE_WIDTH=0
                ASTCENC_SSE=41
                ASTCENC_AVX=3
                ASTCENC_X86_GATHERS=$<BOOL:${ASTCENC_X86_GATHERS}>
                ASTCENC_POPCNT=1
                ASTCENC_F16C=1)

        if (${ASTCENC_VENEER_TYPE} GREATER 0)
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    $<${is_gnu_fe}:-msse2>
                    $<${is_gnu_fe}:-mno-sse4.1>
                    $<${is_gnu_fe}:-Wno-unused-command-line-argument>)
        else()
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    -mavx512f -mavx512bw -mavx512dq -mavx512vl -mavx512vbmi -mbmi -mbmi2 -mpopcnt -mf16c
                    $<${is_gnu_fe}:-Wno-unused-command-line-argument>)
        endif()

        # Non-invariant builds enable FMA (included in AVX-512F). This reduces
        # image quality by up to 0.2 dB (normally much less) in exchange for
        # higher encode throughput.
        if((NOT ${ASTCENC_INVARIANCE}) AND (NOT ${ASTCENC_VENEER_TYPE}))
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    $<${is_gnu_fe}:-mfma>)
        endif()

    elseif(${ASTCENC_ISA_SIMD} MATCHES "avx2")
        # Gathers are quite slow on many x86 microarchitectures, to the point where
        # it can be significantly faster to just avoid them use scalar loads.

        target_compile_definitions(${ASTCENC_TARGET_NAME}
            PRIVATE
                ASTCENC_NEON=0
                ASTCENC_SVE=0
                ASTCENC_SVE_WIDTH=0
                ASTCENC_SSE=41
                ASTCENC_AVX=2
                ASTCENC_X86_GATHERS=$<BOOL:${ASTCENC_X86_GATHERS}>
                ASTCENC_POPCNT=1
                ASTCENC_F16C=1)

        if (${ASTCENC_VENEER_TYPE} GREATER 0)
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
                    $<${is_clangcl}:-mavx2 -mbmi -mbmi2 -mpopcnt -mf16c>
                    $<${is_gnu_fe}:-mavx2 -mbmi -mbmi2 -mpopcnt -mf16c>
                    $<${is_gnu_fe}:-Wno-unused-command-line-argument>)
        endif()

        # Non-invariant builds enable us to loosen the compiler constraints on
        # floating point, but this is only worth doing on CPUs with AVX2 because
        # this implies we can also enable the FMA instruction set extensions
        # which significantly improve performance. Note that this DOES reduce
        # image quality by up to 0.2 dB (normally much less), but buys an
        # average of 10-15% performance improvement ...
        if((NOT ${ASTCENC_INVARIANCE}) AND (NOT ${ASTCENC_VENEER_TYPE}))
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    $<${is_gnu_fe}:-mfma>)
        endif()

    elseif(${ASTCENC_ISA_SIMD} MATCHES "native")
        target_compile_definitions(${ASTCENC_TARGET_NAME}
            PRIVATE)

        if (${ASTCENC_VENEER_TYPE} GREATER 0)
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    $<${is_gnu_fe}:-Wno-unused-command-line-argument>)
        else()
            # Clang can report conflicting AVX10 and AVX-512 features for
            # -march=native on Granite Rapids; keep this diagnostic nonfatal
            target_compile_options(${ASTCENC_TARGET_NAME}
                PRIVATE
                    $<${is_clangcl}:-march=native>
                    $<${is_gnu_fe}:-march=native>
                    $<$<AND:$<CXX_COMPILER_ID:Clang,AppleClang>,$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,18.0.0>>:-Wno-error=invalid-feature-combination>
                    $<${is_gnu_fe}:-Wno-unused-command-line-argument>)
        endif()
    endif()

endfunction()
