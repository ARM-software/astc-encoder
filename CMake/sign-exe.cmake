# This confidential and proprietary software may be used only as
# authorised by a licensing agreement from Arm Limited.
#    Copyright 2020 Arm Ltd. All Rights Reserved.
# The entire notice above must be reproduced on all authorised
# copies and copies may only be made to the extent permitted
# by a licensing agreement from Arm Limited.

# Sign an executable as appropriate for platform.
# Usage: sign_target(<target> [USERNAME <username>] [DESCRIPTION <description>]
#                   [DIGEST_ALGORITHM <digest-algorithm>])
# Currently only signs executables on Windows. Does nothing on other platforms.
function(sign_target EXE_TARGET)
    cmake_parse_arguments(PARSE_ARGV 1
                          ST # Output variable prefix
                          "" # Option arguments
                          "USERNAME;DESCRIPTION;DIGEST_ALGORITHM" # Single-arg keywords
                          "") # Multi-arg keywords
    if(DEFINED ST_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "sign_target: Invalid arguments: ${ST_UNPARSED_ARGUMENTS}")
    endif()
    if(DEFINED ST_KEYWORDS_MISSING_VALUES)
        foreach(KEYWORD ${ST_KEYWORDS_MISSING_VALUES})
            message(FATAL_ERROR "sign_target: ${KEYWORD} is missing an argument")
        endforeach()
    endif()
    if(NOT TARGET ${EXE_TARGET})
        message(FATAL_ERROR "sign_target: First argument must be a CMake target. Did you meant to use sign_file()")
    endif()

    if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows" AND ${ENABLE_CODE_SIGNING})
        set(ST_SCRIPT_ARGS "")
        if(DEFINED ST_USERNAME)
            set(ST_SCRIPT_ARGS ${ST_SCRIPT_ARGS} --username ${ST_USERNAME})
        endif()
        if(DEFINED ST_DESCRIPTION)
            set(ST_SCRIPT_ARGS ${ST_SCRIPT_ARGS} --description ${ST_DESCRIPTION})
        endif()
        if(DEFINED ST_DIGEST_ALGORITHM)
            set(ST_SCRIPT_ARGS ${ST_SCRIPT_ARGS} --digest-algorithm ${ST_DIGEST_ALGORITHM})
        endif()

        # add_custom_command(TARGET...) only works on targets created in the same directory
        # as the add_custom_command command (for unexplained reasons).
        if(TARGET ${EXE_TARGET})
            get_property(TARGET_SOURCE_DIR
                         TARGET ${EXE_TARGET}
                         PROPERTY SOURCE_DIR)
            if(${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${TARGET_SOURCE_DIR})
                set(SIGN_TARGET ${EXE_TARGET})
                set(WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
            else()
                # The target is not from the current directory (could be from a sub-project)
                # so we need to create a dummy target to attach the command to.
                set(SIGN_TARGET sign_${EXE_TARGET})
                add_custom_target(${SIGN_TARGET} ALL)
                add_dependencies(${SIGN_TARGET} ${EXE_TARGET})
                get_property(WORKING_DIRECTORY
                             TARGET ${EXE_TARGET}
                             PROPERTY BINARY_DIR)
            endif()
        endif()
        add_python_command(TARGET ${SIGN_TARGET} POST_BUILD
                           WORKING_DIRECTORY ${WORKING_DIRECTORY}
                           COMMENT "Signing ${EXE_TARGET}"
                           SCRIPT ${CMAKE_SOURCE_DIR}/scripts/sign-exe.py
                                  # Overwrite the original exe.
                                  $<TARGET_FILE:${EXE_TARGET}> $<TARGET_FILE:${EXE_TARGET}>
                                  ${ST_SCRIPT_ARGS}
                           PACKAGES requests)
    endif()
endfunction()

# Sign an executable as appropriate for platform.
# Usage: sign_file(<input_filename> <output_filename> [USERNAME <username>]
#                  [DESCRIPTION <description>] [DIGEST_ALGORITHM <digest-algorithm>])
# Currently only signs executables on Windows. On other platforms the input file will
# just be copied to the output file.
# In order to ensure that signing is carried out the caller must create a custom target
# which depends on the output file:
#   sign_file(foo.exe foo.exe.signed)
#   add_custom_target(sign_foo ALL DEPENDS foo.exe.signed)
function(sign_file INPUT_FILENAME OUTPUT_FILENAME)
    cmake_parse_arguments(PARSE_ARGV 2
                          SF # Output variable prefix
                          "" # Option arguments
                          "USERNAME;DESCRIPTION;DIGEST_ALGORITHM" # Single-arg keywords
                          "DEPENDS") # Multi-arg keywords
    if(DEFINED SF_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "sign_file: Invalid arguments: ${SF_UNPARSED_ARGUMENTS}")
    endif()
    if(DEFINED SF_KEYWORDS_MISSING_VALUES)
        foreach(KEYWORD ${SF_KEYWORDS_MISSING_VALUES})
            message(FATAL_ERROR "sign_file: ${KEYWORD} is missing an argument")
        endforeach()
    endif()
    file(TO_NATIVE_PATH ${INPUT_FILENAME} INPUT_FILENAME)
    file(TO_NATIVE_PATH ${OUTPUT_FILENAME} OUTPUT_FILENAME)

    if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows" AND ${ENABLE_CODE_SIGNING})
        set(SF_SCRIPT_ARGS "")
        if(DEFINED SF_USERNAME)
            set(SF_SCRIPT_ARGS ${SF_SCRIPT_ARGS} --username ${SF_USERNAME})
        endif()
        if(DEFINED SF_DESCRIPTION)
            set(SF_SCRIPT_ARGS ${SF_SCRIPT_ARGS} --description ${SF_DESCRIPTION})
        endif()
        if(DEFINED SF_DIGEST_ALGORITHM)
            set(SF_SCRIPT_ARGS ${SF_SCRIPT_ARGS} --digest-algorithm ${SF_DIGEST_ALGORITHM})
        endif()

        add_python_command(OUTPUT ${OUTPUT_FILENAME}
                           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                           DEPENDS ${SF_DEPENDS}
                           COMMENT "Signing ${INPUT_FILENAME}"
                           SCRIPT ${CMAKE_SOURCE_DIR}/scripts/sign-exe.py
                                  ${INPUT_FILENAME} ${OUTPUT_FILENAME}
                                  ${SF_SCRIPT_ARGS}
                           PACKAGES requests)
    else()
        add_custom_command(OUTPUT ${OUTPUT_FILENAME}
                           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                           DEPENDS ${SF_DEPENDS}
                           COMMENT "Not signing ${OUTPUT_FILENAME} (signing disabled or not available on platform)"
                           COMMAND ${CMAKE_COMMAND} -E copy
                                   ${INPUT_FILENAME} ${OUTPUT_FILENAME})
    endif()

endfunction()
