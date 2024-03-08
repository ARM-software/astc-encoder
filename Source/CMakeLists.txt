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
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" AND ${ASTCENC_CLI})
    set(CMAKE_CXX_COMPILE_OPTIONS_IPO "-flto")
endif()

if(${ASTCENC_DECOMPRESSOR})
    set(ASTCENC_CODEC dec)
else()
    set(ASTCENC_CODEC enc)
endif()

if (${ASTCENC_ISA_NEON})
    set(ASTCENC_ISA_SIMD "neon")
elseif(${ASTCENC_ISA_AVX2})
    set(ASTCENC_ISA_SIMD "avx2")
elseif(${ASTCENC_ISA_SSE41})
    set(ASTCENC_ISA_SIMD "sse4.1")
elseif(${ASTCENC_ISA_SSE2})
    set(ASTCENC_ISA_SIMD "sse2")
elseif(${ASTCENC_ISA_NONE})
    set(ASTCENC_ISA_SIMD "none")
else()
    set(ASTCENC_ISA_SIMD "native")
endif()

set(ASTCENC_SOURCES
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

set(ASTCENC_CLI_SOURCES
    astcenccli_entry.cpp
    astcenccli_error_metrics.cpp
    astcenccli_image.cpp
    astcenccli_image_external.cpp
    astcenccli_image_load_store.cpp
    astcenccli_platform_dependents.cpp
    astcenccli_toplevel.cpp
    astcenccli_toplevel_help.cpp)


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

string(TIMESTAMP astcencoder_YEAR "%Y")
configure_file(
    astcenccli_version.h.in
    astcenccli_version.h
    ESCAPE_QUOTES @ONLY)


# Install the public interface
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Include/
        DESTINATION include
        FILES_MATCHING PATTERN "*.h")

target_include_directories(common INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Include>
    $<INSTALL_INTERFACE:./include>)

# For universal builds, we build a static library, a shared library, and a cli (linked against the shared, so the cli is standalone and can be copied around)
if (${ASTCENC_UNIVERSAL_BUILD})
    message(STATUS "Foo")
    # For legacy naming compatibility
    set(ASTCENC_TARGET astc${ASTCENC_CODEC})

    set(ASTCENC_ARCHITECTURES arm64 x86_64h x86_64)
    # Create all the arch-specific stargets
    foreach(ASTCENC_ARCHITECTURE ${ASTCENC_ARCHITECTURES})
        set(CMAKE_OSX_ARCHITECTURES ${ASTCENC_ARCHITECTURE})
        # shared library
        set(TARGET_NAME astcenc-${ASTCENC_ARCHITECTURE})
        add_library(${TARGET_NAME} SHARED ${ASTCENC_SOURCES})
        target_simd_options_by_arch(${TARGET_NAME} ${ASTCENC_ARCHITECTURE})
        target_compile_definitions(${TARGET_NAME} PRIVATE ASTCENC_DYNAMIC_LIBRARY=1)

        # static library (for executable)
        set(TARGET_NAME astcenc-static-${ASTCENC_ARCHITECTURE})
        add_library(${TARGET_NAME} STATIC ${ASTCENC_SOURCES})
        target_simd_options_by_arch(${TARGET_NAME} ${ASTCENC_ARCHITECTURE})

        if(${ASTCENC_CLI})
            set(TARGET_NAME cli-${ASTCENC_ARCHITECTURE})
            add_executable(${TARGET_NAME} ${ASTCENC_CLI_SOURCES})
            target_include_directories(${TARGET_NAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
            target_compile_options(${TARGET_NAME} PRIVATE $<${is_msvc_fe}:/W3>)
            target_simd_options_by_arch(${TARGET_NAME} ${ASTCENC_ARCHITECTURE})
            target_link_libraries(${TARGET_NAME} PUBLIC astcenc-static-${ASTCENC_ARCHITECTURE})
        endif()
    endforeach()


    # Create the output binary using lipo
    set(ASTCENC_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/${ASTCENC_TARGET}-shared.dylib)
    add_custom_target(astcenc
        ALL
        COMMAND lipo -create -output ${ASTCENC_LIBRARY} -arch x86_64 $<TARGET_FILE:astcenc-x86_64> -arch x86_64h $<TARGET_FILE:astcenc-x86_64h> -arch arm64 $<TARGET_FILE:astcenc-arm64>
    )
    add_dependencies(astcenc
        astcenc-x86_64
        astcenc-x86_64h
        astcenc-arm64
    )
    install(PROGRAMS ${ASTCENC_LIBRARY} DESTINATION bin)

    # Handle the CLI
    if(${ASTCENC_CLI})
        set(CLI_EXEC ${CMAKE_CURRENT_BINARY_DIR}/${ASTCENC_TARGET})
        add_custom_target(cli ALL
            COMMAND lipo -create -output ${CLI_EXEC} -arch x86_64 $<TARGET_FILE:cli-x86_64> -arch x86_64h $<TARGET_FILE:cli-x86_64h> -arch arm64 $<TARGET_FILE:cli-arm64>
        )
        add_dependencies(cli
            astcenc
            cli-x86_64
            cli-x86_64h
            cli-arm64
        )
        install(PROGRAMS ${CLI_EXEC} DESTINATION bin)
    endif()
else()
    # For legacy naming compatibility
    set(ASTCENC_TARGET astc${ASTCENC_CODEC}-${ASTCENC_ISA_SIMD})
    add_library(astcenc ${ASTCENC_SOURCES})
    target_include_directories(astcenc PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Include>
        $<INSTALL_INTERFACE:./include>)

    if (${ASTCENC_LEGACY_NAMING})
        if (${BUILD_SHARED_LIBS})
            set_target_properties(astcenc PROPERTIES OUTPUT_NAME ${ASTCENC_TARGET}-shared)
        else()
            set_target_properties(astcenc PROPERTIES OUTPUT_NAME ${ASTCENC_TARGET})
        endif()
    endif()
    target_simd_options_by_name(astcenc ${ASTCENC_ISA_SIMD})

    if(${ASTCENC_CLI})
        add_executable(cli ${ASTCENC_CLI_SOURCES})
        if (${ASTCENC_LEGACY_NAMING})
            set_target_properties(cli PROPERTIES OUTPUT_NAME ${ASTCENC_TARGET})
        else()
            set_target_properties(cli PROPERTIES OUTPUT_NAME astcenc)
        endif()
        target_simd_options_by_name(cli ${ASTCENC_ISA_SIMD})
        target_include_directories(cli PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
        target_link_libraries(cli PUBLIC astcenc)
    endif()

    include(CMakePackageConfigHelpers)

    # Create an export type "astcencTargets" detailing the targetable artifacts created
    # by this project.
    # install(TARGETS
    #         EXPORT simdTargets)
    install(TARGETS astcenc common simd-avx2 simd-neon simd-sse2 simd-sse41 simd-none
            EXPORT astcencTargets)

    # Install the export targets as a CMake config file in the share/astcenc folder so
    # that they can referenced  by downstream projects as `astcenc::astcenc` after a
    # successful `find_package` call.
    install(EXPORT astcencTargets
            NAMESPACE astcenc::
            FILE astcencTargets.cmake
            DESTINATION share/astcenc/)

    # Generate the config file that includes the exports.
    configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/../cmake/Config.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/astcencConfig.cmake"
        INSTALL_DESTINATION "share/astcenc/"
        NO_SET_AND_CHECK_MACRO
        NO_CHECK_REQUIRED_COMPONENTS_MACRO)

    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/astcencConfig.cmake
            DESTINATION "share/astcenc/")

    if(${ASTCENC_CLI})
        install(TARGETS cli)
    endif()
endif()


# - - - - - - - - - - - - - - - - - -
# Unit testing
if(${ASTCENC_UNITTEST})
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    add_subdirectory(GoogleTest)
    enable_testing()
    add_subdirectory(UnitTest)
endif()
