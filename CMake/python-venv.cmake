# This confidential and proprietary software may be used only as
# authorised by a licensing agreement from Arm Limited.
#    Copyright 2020 Arm Ltd. All Rights Reserved.
# The entire notice above must be reproduced on all authorised
# copies and copies may only be made to the extent permitted
# by a licensing agreement from Arm Limited.
#

# CMake support code to manage running Python scripts.

find_package(Python3 COMPONENTS Interpreter)

if(${USE_PYTHON_VENV})
    # Create a virtual environment to contain any Python dependencies.
    # This is run once at configure time.
    execute_process(COMMAND ${Python3_EXECUTABLE} -m venv python-venv
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

    set(PYTHON_VENV_DIR ${CMAKE_BINARY_DIR}/python-venv)
    if(${CMAKE_SYSTEM_NAME} STREQUAL Windows)
        set(PYTHON_EXECUTABLE ${PYTHON_VENV_DIR}/Scripts/python.exe)
    else()
        set(PYTHON_EXECUTABLE ${PYTHON_VENV_DIR}/bin/python)
    endif()

    set(PYTHON_VENV_INSTALLED_PACKAGES "" CACHE INTERNAL "")

    # Install a package into the venv at configure time.
    function(pip_install)
        foreach(PACKAGE ${ARGV})
            list(FIND PYTHON_VENV_INSTALLED_PACKAGES ${PACKAGE} INDEX)
            if(INDEX EQUAL -1)
                # We haven't installed this package already.
                set(PYTHON_VENV_INSTALLED_PACKAGES "${PYTHON_VENV_INSTALLED_PACKAGES};${PACKAGE}"
                    CACHE INTERNAL "")
                list(APPEND PYTHON_VENV_INSTALLED_PACKAGES ${PACKAGE})
                message("Installing ${PACKAGE} into python-venv")
                execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install ${PACKAGE}
                                WORKING_DIRECTORY ${PYTHON_VENV_DIR})
            endif()
        endforeach()
    endfunction()
else()

    set(PYTHON_FOUND_PACKAGES "" CACHE INTERNAL "")

    # Check's whether a package is available, but doesn't install it.
    function(pip_install)
        foreach(PACKAGE ${ARGV})
            list(FIND PYTHON_FOUND_PACKAGES ${PACKAGE} INDEX)
            if(INDEX EQUAL -1)
                # Check whether the package is installed
                execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip show --quiet ${PACKAGE}
                                RESULT_VARIABLE RET_CODE)
                if(NOT ${RET_CODE} EQUAL 0)
                    message(SEND_ERROR "Missing Python package ${PACKAGE}")
                else()
                    message(STATUS "Found Python package ${PACKAGE}")
                    set(PYTHON_FOUND_PACKAGES "${PYTHON_FOUND_PACKAGES};${PACKAGE}"
                        CACHE INTERNAL "")
                endif()
            endif()
        endforeach()

    endfunction()
    # VENV is disabled so just use the system Python executable.
    set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
endif()

# Function which wraps add_custom_command() to run a Python script within the venv.
function(add_python_command)
    cmake_parse_arguments(APC # Output variable prefix
                          "" # Option arguments
                          "" # Single-arg keywords
                          "SCRIPT;PACKAGES" # Multi-arg keywords
                          ${ARGV})
    if(NOT DEFINED APC_SCRIPT)
        message(SEND_ERROR "add_python_command: SCRIPT must be specified.")
    endif()
    if(DEFINED APC_PACKAGES)
        pip_install(${APC_PACKAGES})
    endif()
    add_custom_command(COMMAND ${PYTHON_EXECUTABLE} ${APC_SCRIPT}
                       ${APC_UNPARSED_ARGUMENTS})
endfunction()
