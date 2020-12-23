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

add_executable(test-simd-${ISA_SIMD})

target_sources(test-simd-${ISA_SIMD}
    PRIVATE
        test_simd.cpp)

target_include_directories(test-simd-${ISA_SIMD}
    PRIVATE
        ${gtest_SOURCE_DIR}/include)

target_compile_options(test-simd-${ISA_SIMD}
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
    target_compile_definitions(test-simd-${ISA_SIMD}
        PRIVATE
            ASTCENC_NEON=0
            ASTCENC_SSE=0
            ASTCENC_AVX=0
            ASTCENC_POPCNT=0)

    if (${ARCH} MATCHES x64)
        target_compile_options(astcenc-${ISA_SIMD}
            PRIVATE
                $<$<CXX_COMPILER_ID:${GNU_LIKE}>:-mfpmath=sse -msse2>)
    endif()

elseif(${ISA_SIMD} MATCHES "neon")
    target_compile_definitions(test-simd-${ISA_SIMD}
        PRIVATE
            ASTCENC_NEON=1
            ASTCENC_SSE=0
            ASTCENC_AVX=0
            ASTCENC_POPCNT=0)

elseif(${ISA_SIMD} MATCHES "sse2")
    target_compile_definitions(test-simd-${ISA_SIMD}
        PRIVATE
            ASTCENC_NEON=0
            ASTCENC_SSE=20
            ASTCENC_AVX=0
            ASTCENC_POPCNT=0)

    target_compile_options(test-simd-${ISA_SIMD}
        PRIVATE
        $<$<CXX_COMPILER_ID:${GNU_LIKE}>:-mfpmath=sse -msse2>)

elseif(${ISA_SIMD} MATCHES "sse4.1")
    target_compile_definitions(test-simd-${ISA_SIMD}
        PRIVATE
            ASTCENC_NEON=0
            ASTCENC_SSE=41
            ASTCENC_AVX=0
            ASTCENC_POPCNT=1)

    target_compile_options(test-simd-${ISA_SIMD}
        PRIVATE
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-mfpmath=sse -msse4.1 -mpopcnt>)

elseif(${ISA_SIMD} MATCHES "avx2")
    target_compile_definitions(test-simd-${ISA_SIMD}
        PRIVATE
            ASTCENC_NEON=0
            ASTCENC_SSE=41
            ASTCENC_AVX=2
            ASTCENC_POPCNT=1)

    target_compile_options(test-simd-${ISA_SIMD}
        PRIVATE
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-mfpmath=sse -mavx2 -mpopcnt>
            $<$<CXX_COMPILER_ID:MSVC>:/arch:AVX2>)
endif()

target_link_libraries(test-simd-${ISA_SIMD}
    PRIVATE
        gtest_main)

add_test(NAME test-simd-${ISA_SIMD}
         COMMAND test-simd-${ISA_SIMD})

install(TARGETS test-simd-${ISA_SIMD} DESTINATION ${PACKAGE_ROOT})
