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

# CMake configuration
cmake_minimum_required(VERSION 3.15)
cmake_policy(SET CMP0069 NEW)  # LTO support
cmake_policy(SET CMP0091 NEW)  # MSVC runtime support

project(astcencoder VERSION 3.0.0)

option(ISA_AVX2 "Enable builds for AVX2 SIMD")
option(ISA_SSE41 "Enable builds for SSE4.1 SIMD")
option(ISA_SSE2 "Enable builds for SSE2 SIMD")
option(ISA_NEON "Enable builds for NEON SIMD")
option(ISA_NONE "Enable builds for no SIMD")
# Manually provide this as a tristate so we can distinguish the default
# not set state from explicitly ON and OFF ...
set(ISA_NATIVE AUTO CACHE STRING "Enable builds for the native local CPU")
set_property(CACHE ISA_NATIVE PROPERTY STRINGS AUTO ON OFF)

set(UNIVERSAL_BUILD OFF)
set(X64_ISA_COUNT 0)

set(CONFIGS ${ISA_NONE} ${ISA_NATIVE} ${ISA_AVX2} ${ISA_SSE41} ${ISA_SSE2})
foreach(CONFIG ${CONFIGS})
    if(${CONFIG})
        math(EXPR X64_ISA_COUNT "${X64_ISA_COUNT} + 1")
    endif()
endforeach()

math(EXPR BUILD_NON_UNIVERSAL "${X64_ISA_COUNT} >> 1")

# Case when more than one x64 ISA and NEON is specified
# This prevents something like ISA_AVX2 ISA_SSE41 ISA_NEON all specified
if (${BUILD_NON_UNIVERSAL} AND ${ISA_NEON})
        message(FATAL_ERROR "Conflicting ISAs specified. Only one of ISA_NONE|ISA_AVX2|ISA_SSE41|ISA_SSE2 is valid with ISA_NEON.")
endif()

# When universal build, either provide both ISA_<x64> and ISA_NEON or none
if (CMAKE_OSX_ARCHITECTURES)
    if ((${X64_ISA_COUNT} EQUAL 1) AND (NOT ${ISA_NEON}))
            message(FATAL_ERROR "Incomplete ISAs specified for universal build. Use ISA_NEON and one of ISA_NONE|ISA_AVX2|ISA_SSE41|ISA_SSE2 for universal build.")
    endif()

    if ((${X64_ISA_COUNT} EQUAL 0) AND ${ISA_NEON})
            message(FATAL_ERROR "Incomplete ISAs specified for universal build. Use ISA_NEON and one of ISA_NONE|ISA_AVX2|ISA_SSE41|ISA_SSE2 for universal build.")
    endif()

    if (${BUILD_NON_UNIVERSAL})
        message(FATAL_ERROR "Conflicting ISAs specified for universal build. Use ISA_NEON and one of ISA_NONE|ISA_AVX2|ISA_SSE41|ISA_SSE2 for universal build.")
    endif()
endif()

# Build non-universal if only one ISA_<x64> is specified
if ((${X64_ISA_COUNT} EQUAL 1) AND
    (NOT CMAKE_OSX_ARCHITECTURES) AND (NOT ${ISA_NEON}))
        set(BUILD_NON_UNIVERSAL 1)
endif()

if ((${X64_ISA_COUNT} EQUAL 0) AND
    (NOT CMAKE_OSX_ARCHITECTURES) AND ${ISA_NEON})
        set(BUILD_NON_UNIVERSAL 1)
endif()

if (NOT ${BUILD_NON_UNIVERSAL})
    # Try to work out if universal build is requested
    if (CMAKE_OSX_ARCHITECTURES)
        # Default to SSE2 if no other x64 ISA is specified
        if ((NOT ${ISA_NONE}) OR (NOT ${ISA_NATIVE}) OR (NOT ${ISA_AVX2}) OR (NOT ${ISA_SSE41}))
            set(ISA_SSE2 ON)
        endif()
        set(ISA_NEON ON)
        set(UNIVERSAL_BUILD ON)
    else()
        if ((${X64_ISA_COUNT} EQUAL 1) AND (${ISA_NEON}))
            set(UNIVERSAL_BUILD ON)
            set(CMAKE_OSX_ARCHITECTURES "arm64;x86_64")
        endif()
    endif()

    if (${UNIVERSAL_BUILD})
        message(STATUS "Building universal binary for ${CMAKE_OSX_ARCHITECTURES}")
    endif()

    foreach(value ${CMAKE_OSX_ARCHITECTURES})
        string(APPEND ARCH_STRING "_${value}")
    endforeach()
else()
    if (${ISA_NEON})
        string(APPEND ARCH_STRING "_aarch64")
        set(CMAKE_OSX_ARCHITECTURES "arm64")
    else()
        string(APPEND ARCH_STRING "_x64")
        set(CMAKE_OSX_ARCHITECTURES "x86_64")
    endif()
endif()

function(printopt optName optVal)
    if(${optVal})
        message(STATUS "  ${optName} backend - ON")
    else()
        message(STATUS "  ${optName} backend - OFF")
    endif()
endfunction()

printopt("AVX2"  ${ISA_AVX2})
printopt("SSE41" ${ISA_SSE41})
printopt("SSE2"  ${ISA_SSE2})
printopt("NEON"  ${ISA_NEON})
printopt("NONE"  ${ISA_NONE})

if((${ISA_NATIVE} STREQUAL AUTO) OR (${ISA_NATIVE} STREQUAL ON))
    set(ISA_NATIVE 1)
    message(STATUS "  Native SIMD backend - ON")
else()
    set(ISA_NATIVE 0)
    message(STATUS "  Native SIMD backend - OFF")
endif()

option(DECOMPRESSOR "Enable builds for decompression only")
if(${DECOMPRESSOR})
    message(STATUS "  Decompress-only backend - ON")
else()
    message(STATUS "  Decompress-only backend - OFF")
endif()

option(DIAGNOSTICS "Enable builds for diagnostic trace")
if(${DIAGNOSTICS})
    message(STATUS "  Diagnostic tracing - ON")
else()
    message(STATUS "  Diagnostic tracing - OFF")
endif()

option(UNITTEST "Enable builds for unit tests")
if(${UNITTEST})
    message(STATUS "  Unit tests - ON")
else()
    message(STATUS "  Unit tests - OFF")
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(PACKAGE_ROOT astcenc)
include(CTest)

# Subcomponents
add_subdirectory(Source)

# Configure package archive
if("${CMAKE_SYSTEM_NAME}" STREQUAL "Darwin")
    string(TOLOWER "macOS" PKG_OS)
else()
    string(TOLOWER ${CMAKE_SYSTEM_NAME} PKG_OS)
endif()

set(PKG_VER ${CMAKE_PROJECT_VERSION_MAJOR}.${CMAKE_PROJECT_VERSION_MINOR})

set(CPACK_PACKAGE_FILE_NAME "astcenc-${PKG_VER}-${PKG_OS}${ARCH_STRING}")
set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY FALSE)
set(CPACK_PACKAGE_CHECKSUM SHA256)
set(CPACK_GENERATOR ZIP)

include(CPack) # Must be included after CPack configuration.
