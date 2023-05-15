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

# Overwrite the LTO flags to force fat LTO; worth 3-4% performance
# See https://gitlab.kitware.com/cmake/cmake/-/issues/16808
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" AND ${CLI})
    set(CMAKE_CXX_COMPILE_OPTIONS_IPO "-flto")
endif()

if(${ASTCENC_DECOMPRESSOR})
    set(ASTCENC_CODEC dec)
else()
    set(ASTCENC_CODEC enc)
endif()

if(${ASTCENC_UNIVERSAL_BUILD})
    if(${ASTCENC_ISA_AVX2})
        set(ASTCENC_ISA_SIMD "avx2")
    elseif(${ASTCENC_ISA_SSE41})
        set(ASTCENC_ISA_SIMD "sse4.1")
    elseif(${ASTCENC_ISA_SSE2})
        set(ASTCENC_ISA_SIMD "sse2")
    endif()
    include(cmake_core.cmake)
else()
    set(ARTIFACTS native none neon avx2 sse4.1 sse2)
    set(CONFIGS ${ASTCENC_ISA_NATIVE} ${ASTCENC_ISA_NONE} ${ASTCENC_ISA_NEON} ${ASTCENC_ISA_AVX2} ${ASTCENC_ISA_SSE41} ${ASTCENC_ISA_SSE2})
    list(LENGTH ARTIFACTS ARTIFACTS_LEN)
    math(EXPR ARTIFACTS_LEN "${ARTIFACTS_LEN} - 1")

    foreach(INDEX RANGE ${ARTIFACTS_LEN})
        list(GET ARTIFACTS ${INDEX} ARTIFACT)
        list(GET CONFIGS ${INDEX} CONFIG)
        if(${CONFIG})
            set(ISA_SIMD ${ARTIFACT})
            include(cmake_core.cmake)
        endif()
    endforeach()
endif()

# - - - - - - - - - - - - - - - - - -
# Unit testing
if(${ASTCENC_UNITTEST})
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    add_subdirectory(GoogleTest)
    enable_testing()
    add_subdirectory(UnitTest)
endif()
