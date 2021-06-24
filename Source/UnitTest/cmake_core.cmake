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

if (${UNIVERSAL_BUILD})
    set(astc_test test-simd)
else()
    set(astc_test test-simd-${ISA_SIMD})
endif()

add_executable(${astc_test})

target_sources(${astc_test}
    PRIVATE
        test_simd.cpp
        ../astcenc_mathlib_softfloat.cpp)

target_include_directories(${astc_test}
    PRIVATE
        ${gtest_SOURCE_DIR}/include)

target_compile_options(${astc_test}
    PRIVATE
        # Use pthreads on Linux/macOS
        $<$<PLATFORM_ID:Linux,Darwin>:-pthread>

        # MSVC compiler defines
        $<$<CXX_COMPILER_ID:MSVC>:/EHsc>

        # G++ and Clang++ compiler defines
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wextra>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wpedantic>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Werror>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wshadow>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wdouble-promotion>)

# Set up configuration for SIMD ISA builds
if(${ISA_SIMD} MATCHES "none")
    if (NOT ${UNIVERSAL_BUILD})
        target_compile_definitions(${astc_test}
            PRIVATE
                ASTCENC_NEON=0
                ASTCENC_SSE=0
                ASTCENC_AVX=0
                ASTCENC_POPCNT=0
                ASTCENC_F16C=0)
    endif()

    if (${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
        target_compile_options(${astc_test}
            PRIVATE
                $<$<CXX_COMPILER_ID:${GNU_LIKE}>:-mfpmath=sse -msse2>)

    endif()

elseif(${ISA_SIMD} MATCHES "neon")
    if (NOT ${UNIVERSAL_BUILD})
        target_compile_definitions(${astc_test}
            PRIVATE
                ASTCENC_NEON=1
                ASTCENC_SSE=0
                ASTCENC_AVX=0
                ASTCENC_POPCNT=0
                ASTCENC_F16C=0)
    endif()

elseif(${ISA_SIMD} MATCHES "sse2")
    if (NOT ${UNIVERSAL_BUILD})
        target_compile_definitions(${astc_test}
            PRIVATE
                ASTCENC_NEON=0
                ASTCENC_SSE=20
                ASTCENC_AVX=0
                ASTCENC_POPCNT=0
                ASTCENC_F16C=0)
    endif()

    target_compile_options(${astc_test}
        PRIVATE
        $<$<CXX_COMPILER_ID:${GNU_LIKE}>:-mfpmath=sse -msse2>)

elseif(${ISA_SIMD} MATCHES "sse4.1")
    if (NOT ${UNIVERSAL_BUILD})
        target_compile_definitions(${astc_test}
            PRIVATE
                ASTCENC_NEON=0
                ASTCENC_SSE=41
                ASTCENC_AVX=0
                ASTCENC_POPCNT=1
                ASTCENC_F16C=0)
    endif()

    target_compile_options(${astc_test}
        PRIVATE
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-mfpmath=sse -msse4.1 -mpopcnt>)

elseif(${ISA_SIMD} MATCHES "avx2")
    if (NOT ${UNIVERSAL_BUILD})
        target_compile_definitions(${astc_test}
            PRIVATE
                ASTCENC_NEON=0
                ASTCENC_SSE=41
                ASTCENC_AVX=2
                ASTCENC_POPCNT=1
                ASTCENC_F16C=1)
    endif()

    target_compile_options(${astc_test}
        PRIVATE
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-mfpmath=sse -mavx2 -mpopcnt -mf16c>
            $<$<CXX_COMPILER_ID:MSVC>:/arch:AVX2>)

endif()

target_link_libraries(${astc_test}
    PRIVATE
        gtest_main)

add_test(NAME ${astc_test}
         COMMAND ${astc_test})

install(TARGETS ${astc_test} DESTINATION ${PACKAGE_ROOT})
