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

include(../cmake_compiler.cmake)

set(ASTCENC_TEST test-unit-${ASTCENC_ISA_SIMD})

add_executable(${ASTCENC_TEST})

# Enable LTO under the conditions where the codec library will use LTO.
# The library link will fail if the settings don't match
if(${ASTCENC_CLI})
    set_property(TARGET ${ASTCENC_TEST}
        PROPERTY
            INTERPROCEDURAL_OPTIMIZATION_RELEASE True)
endif()

# Use a static runtime on MSVC builds (ignored on non-MSVC compilers)
set_property(TARGET ${ASTCENC_TEST}
    PROPERTY
        MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")


target_sources(${ASTCENC_TEST}
    PRIVATE
        test_simd.cpp
        test_softfloat.cpp
        test_decode.cpp
        test_encode.cpp
        test_image.cpp)

target_include_directories(${ASTCENC_TEST}
    PRIVATE
        ${gtest_SOURCE_DIR}/include)

target_link_libraries(${ASTCENC_TEST}
    PRIVATE
        astcenc-${ASTCENC_ISA_SIMD}-static)

astcenc_set_compiler_options(${ASTCENC_TEST} ${ASTCENC_ISA_SIMD} 0)

# Ignore warnings from the GoogleTest headers
target_compile_options(${ASTCENC_TEST}
    PRIVATE
        $<${is_gnu_fe}:-Wno-double-promotion>
        $<${is_gnu_fe}:-Wno-undef>
        $<${is_gnu_fe}:-Wno-global-constructors>
        $<${is_gnu_fe}:-Wno-unsafe-buffer-usage>
        $<${is_clang}:-Wno-switch-default>)

target_link_libraries(${ASTCENC_TEST}
    PRIVATE
        gtest_main)

add_test(NAME ${ASTCENC_TEST}
         COMMAND ${ASTCENC_TEST})
