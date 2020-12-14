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


# CMake configuration
cmake_minimum_required(VERSION 3.15)
cmake_policy(SET CMP0069 NEW)  # LTO support
cmake_policy(SET CMP0091 NEW)  # MSVC runtime support

# Project configuration
project(astcenc VERSION 2.2.0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

set(PACKAGE_ROOT astcenc)

# Validate command line parameters
set(ISA_SIMD "avx2" CACHE STRING "SIMD instruction set")
message("-- ISA SIMD is set to ${ISA_SIMD}")
set(VALID_ISA_SIMD none sse2 sse4.1 avx2)
list(FIND VALID_ISA_SIMD ${ISA_SIMD} index)
if(index EQUAL -1)
    message(FATAL_ERROR "ISA_SIMD must be one of ${VALID_ISA_SIMD}")
endif()

set(ISA_INVARIANCE OFF CACHE BOOL "ISA invariance")
message("-- ISA invariance is set to ${ISA_INVARIANCE}")
set(VALID_ISA_INVARIANCE ON OFF)
list(FIND VALID_ISA_INVARIANCE ${ISA_INVARIANCE} index)
if(index EQUAL -1)
    message(FATAL_ERROR "ISA_INVARIANCE must be one of ${VALID_ISA_INVARIANCE}")
endif()


# Include subcomponenets
add_subdirectory(Source)

# Configure package archive
if("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
    set(PACKAGE_PLATFORM "macOS")
else()
    set(PACKAGE_PLATFORM ${CMAKE_SYSTEM_NAME})
endif()

set(CPACK_PACKAGE_FILE_NAME "astcenc_${PACKAGE_PLATFORM}_x64")
set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY FALSE)
if("${CMAKE_SYSTEM_NAME}" STREQUAL "Windows")
    set(CPACK_GENERATOR ZIP)
else()
    set(CPACK_GENERATOR TGZ)
endif()

include(CPack) # Must be included after CPack configuration.
