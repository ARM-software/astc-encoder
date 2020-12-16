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

# Command line configuration
option(ISA_AVX2 "Enable builds for AVX2 SIMD")

message("-- Selecting backend build type(s)")
set(ANY_ISA 0)
if(${ISA_AVX2})
    set(ANY_ISA 1)
    message("  -- AVX2 backend: ON")
endif()

option(ISA_SSE41 "Enable builds for SSE4.1 SIMD")
if(${ISA_SSE41})
    set(ANY_ISA 1)
    message("  -- SSE4.1 backend: ON")
endif()

option(ISA_SSE2 "Enable builds for SSE2 SIMD")
if(${ISA_SSE2})
    set(ANY_ISA 1)
    message("  -- SSE2 backend: ON")
endif()

option(ISA_NONE "Enable builds for noSIMD")
if(${ISA_NONE})
    set(ANY_ISA 1)
    message("  -- No SIMD backend: ON")
endif()

option(ISA_INVARIANCE "Enable builds for ISA invariance")
if(${ISA_INVARIANCE})
    message("  -- ISA invariance: ON")
endif()

if(NOT ${ANY_ISA})
    message(FATAL_ERROR "At least one backend ISA must be enabled")
endif()

# Subcomponents
add_subdirectory(Source)

# Configure package archive
if("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
    string(TOLOWER "macOS" PKG_OS)
else()
    string(TOLOWER ${CMAKE_SYSTEM_NAME} PKG_OS)
endif()

set(PKG_VER ${CMAKE_PROJECT_VERSION_MAJOR}.${CMAKE_PROJECT_VERSION_MINOR})

set(CPACK_PACKAGE_FILE_NAME "astcenc-${PKG_VER}-${PKG_OS}-x64")
set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY FALSE)
set(CPACK_PACKAGE_CHECKSUM SHA256)
set(CPACK_GENERATOR ZIP)

include(CPack) # Must be included after CPack configuration.
