#  SPDX-License-Identifier: Apache-2.0
#  ----------------------------------------------------------------------------
#  Copyright 2021 Arm Limited
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
include(ExternalProject)

project(astcencoder_example VERSION 1.0.0)

# Command line configuration
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# Add the external project and pull out the project directories we need

# Add e.g. "CMAKE_CACHE_ARGS -DISA_SSE41:String=ON" to ExternalProject_Add to
# set configure settings usually found on the astcenc configure command line.
# Remember that changing ISA changes the name of the library binary produced.

ExternalProject_Add(astcencoder
    GIT_REPOSITORY https://github.com/ARM-software/astc-encoder
    GIT_TAG main
    CMAKE_CACHE_ARGS -DCLI:String=OFF
    INSTALL_COMMAND "")

ExternalProject_Get_property(astcencoder
    SOURCE_DIR)

ExternalProject_Get_property(astcencoder
    BINARY_DIR)

# Build the command line
add_executable(astcenc_example astc_api_example.cpp)

# ... with astcencoder as a dependency
add_dependencies(astcenc_example astcencoder)

# ... with astcencoder Source dir on the include path
target_include_directories(astcenc_example
    PRIVATE
        ${SOURCE_DIR}/Source)

# ... with astcencoder Binary dir on the library path and as a library dep
target_link_directories(astcenc_example
    PRIVATE
        ${BINARY_DIR}/Source)

target_link_libraries(astcenc_example
    PRIVATE
        astcenc-native-static)
