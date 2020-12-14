# This confidential and proprietary software may be used only as
# authorised by a licensing agreement from Arm Limited.
#    Copyright 2020 Arm Ltd. All Rights Reserved.
# The entire notice above must be reproduced on all authorised
# copies and copies may only be made to the extent permitted
# by a licensing agreement from Arm Limited.

# Set a default build type if none were specified.
# With a single-config generator (Make, Ninja, etc) CMAKE_BUILD_TYPE needs to
# be set to an appropriate value.
# For multi-config generators (Visual Studio, XCode, etc) CMAKE_CONFIGURATION_TYPES
# will be set to a list of available configurations and CMAKE_BUILD_TYPE is not used.
if(NOT CMAKE_BUILD_TYPE AND NOT GENERATOR_IS_MULTI_CONFIG)
    # It is necessary to force the value here because CMake will have already set a
    # value "" in the cache and we need to override it.
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Set the build type" FORCE)
endif()

if(GENERATOR_IS_MULTI_CONFIG)
    message(STATUS "Multi-config generator with build types: ${CMAKE_CONFIGURATION_TYPES}")
else()
    message(STATUS "Build type is set to ${CMAKE_BUILD_TYPE}")
endif()

# Configure settings for running Python scripts
if(${CMAKE_SYSTEM_NAME} STREQUAL Linux)
    # Jenkins Linux builds run in a docker container with the dependency libraries installed
    # and we want the build to use dependencies installed on the system. This can be overridden
    # for developer local builds to allow CMake to manage the Python dependencies.
    set(USE_PYTHON_VENV_DEFAULT OFF)
else()
    set(USE_PYTHON_VENV_DEFAULT ON)
endif()

option(USE_PYTHON_VENV "Run Python scripts inside a virtual environment" ${USE_PYTHON_VENV_DEFAULT})

if(${USE_PYTHON_VENV})
    message(STATUS "Python scripts will be run inside a virtual environment")
else()
    message(STATUS "Python scripts will run in the system environment")
endif()

# Configure code signing
if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    option(ENABLE_CODE_SIGNING "Enable cryptographic signing of libraries and executables" OFF)
    if(${ENABLE_CODE_SIGNING})
        message(STATUS "Code signing enabled")
    else()
        message(STATUS "Code signing disabled")
    endif()
else()
    set(ENABLE_CODE_SIGNING OFF CACHE INTERNAL "")
endif()
