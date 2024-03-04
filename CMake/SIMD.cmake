
# On CMake 3.25 or older CXX_COMPILER_FRONTEND_VARIANT is not always set
if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "")
    set(CMAKE_CXX_COMPILER_FRONTEND_VARIANT "${CMAKE_CXX_COMPILER_ID}")
endif()

# Compiler accepts MSVC-style command line options
set(is_msvc_fe "$<STREQUAL:${CMAKE_CXX_COMPILER_FRONTEND_VARIANT},MSVC>")
# Compiler accepts GNU-style command line options
set(is_gnu_fe1 "$<STREQUAL:${CMAKE_CXX_COMPILER_FRONTEND_VARIANT},GNU>")
# Compiler accepts AppleClang-style command line options, which is also GNU-style
set(is_gnu_fe2 "$<STREQUAL:${CMAKE_CXX_COMPILER_FRONTEND_VARIANT},AppleClang>")
# Compiler accepts GNU-style command line options
set(is_gnu_fe "$<OR:${is_gnu_fe1},${is_gnu_fe2}>")

# Compiler is Visual Studio cl.exe
set(is_msvccl "$<AND:${is_msvc_fe},$<CXX_COMPILER_ID:MSVC>>")
# Compiler is Visual Studio clangcl.exe
set(is_clangcl "$<AND:${is_msvc_fe},$<CXX_COMPILER_ID:Clang>>")
# Compiler is upstream clang with the standard frontend
set(is_clang "$<AND:${is_gnu_fe},$<CXX_COMPILER_ID:Clang,AppleClang>>")

# All the common compiler options for all targets
add_library(common INTERFACE)
target_compile_definitions(common INTERFACE
    $<${is_msvc_fe}:_CRT_SECURE_NO_WARNINGS>)

if(${ASTCENC_DECOMPRESSOR})
        target_compile_definitions(common INTERFACE
            ASTCENC_DECOMPRESS_ONLY)
endif()

if(${ASTCENC_BLOCK_MAX_TEXELS})
    target_compile_definitions(common INTERFACE
        ASTCENC_BLOCK_MAX_TEXELS=${ASTCENC_BLOCK_MAX_TEXELS})
endif()

if(${ASTCENC_DIAGNOSTICS})
    target_compile_definitions(common INTERFACE
        ASTCENC_DIAGNOSTICS)
endif()

target_compile_options(common
    INTERFACE
    # Use pthreads on Linux/macOS
    $<$<PLATFORM_ID:Linux,Darwin>:-pthread>

    # MSVC compiler defines
    $<${is_msvc_fe}:/EHsc>
    $<${is_msvccl}:/wd4324>

    # G++ and Clang++ compiler defines
    $<${is_gnu_fe}:-Wall>
    $<${is_gnu_fe}:-Wextra>
    $<${is_gnu_fe}:-Wpedantic>
    $<${is_gnu_fe}:-Werror>
    $<${is_gnu_fe}:-Wshadow>
    $<${is_gnu_fe}:-Wdouble-promotion>
    $<${is_clang}:-Wdocumentation>

    # Hide noise thrown up by Clang 10 and clang-cl
    $<${is_gnu_fe}:-Wno-unknown-warning-option>
    $<${is_gnu_fe}:-Wno-c++98-compat-pedantic>
    $<${is_gnu_fe}:-Wno-c++98-c++11-compat-pedantic>
    $<${is_gnu_fe}:-Wno-float-equal>
    $<${is_gnu_fe}:-Wno-deprecated-declarations>
    $<${is_gnu_fe}:-Wno-atomic-implicit-seq-cst>

    # Clang 10 also throws up warnings we need to investigate (ours)
    $<${is_gnu_fe}:-Wno-cast-align>
    $<${is_gnu_fe}:-Wno-sign-conversion>
    $<${is_gnu_fe}:-Wno-implicit-int-conversion>
    $<${is_gnu_fe}:-Wno-shift-sign-overflow>
    $<${is_gnu_fe}:-Wno-format-nonliteral>
    $<${is_gnu_fe}:-Wno-reserved-identifier>
    $<${is_gnu_fe}:-Wno-cast-function-type>

    # Force DWARF4 for Valgrind profiling
    $<$<AND:$<PLATFORM_ID:Linux,Darwin>,${is_clang}>:-gdwarf-4>

    # Disable non-portable Windows.h warning (fixing it fails builds on MinGW)
    $<$<AND:$<PLATFORM_ID:Windows>,${is_clang}>:-Wno-nonportable-system-include-path>)

target_link_options(common INTERFACE
        # Use pthreads on Linux/macOS
        $<$<PLATFORM_ID:Linux,Darwin>:-pthread>)

target_compile_options(common INTERFACE $<${is_msvc_fe}:/W4>)

if(${ASTCENC_ASAN})
    target_compile_options(common INTERFACE
            $<${is_clang}:-fsanitize=address>)

    target_link_options(common INTERFACE
            $<${is_clang}:-fsanitize=address>)
endif()

if(NOT ${ASTCENC_INVARIANCE})
    target_compile_definitions(common INTERFACE
        ASTCENC_NO_INVARIANCE=1)

        # For Visual Studio prior to 2022 (compiler < 19.30) /fp:precise
    # For Visual Studio 2022 (compiler >= 19.30) /fp:precise and /fp:contract

    # For Visual Studio 2022 ClangCL seems to have accidentally enabled contraction by default,
    # so behaves differently to CL.exe. Use the -Xclang argument to workaround and allow access
    # GNU-style switch to control contraction on the assumption this gets fixed and disabled.
    # Note ClangCL does not accept /fp:contract as an argument as of v15.0.7.
    target_compile_options(common INTERFACE
            $<${is_msvccl}:/fp:precise>
            $<${is_clangcl}:/fp:precise>
            $<$<AND:${is_msvccl},$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,19.30>>:/fp:contract>
            $<$<AND:${is_clangcl},$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,14.0.0>>:-Xclang -ffp-contract=fast>
            $<$<AND:${is_clang},$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,10.0.0>>:-ffp-model=precise>
            $<${is_gnu_fe}:-ffp-contract=fast>)
else()
    # For Visual Studio prior to 2022 (compiler < 19.30) /fp:strict
    # For Visual Studio 2022 (compiler >= 19.30) /fp:precise

    # For Visual Studio 2022 ClangCL seems to have accidentally enabled contraction by default,
    # so behaves differently to CL.exe. Use the -Xclang argument to workaround and allow access
    # GNU-style switch to control contraction and force disable.
    target_compile_options(common INTERFACE
            $<$<AND:${is_msvccl},$<VERSION_LESS:$<CXX_COMPILER_VERSION>,19.30>>:/fp:strict>
            $<$<AND:${is_msvccl},$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,19.30>>:/fp:precise>
            $<${is_clangcl}:/fp:precise>
            $<$<AND:${is_clangcl},$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,14.0.0>>:-Xclang -ffp-contract=off>
            $<$<AND:${is_clang},$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,10.0.0>>:-ffp-model=precise>
            $<${is_gnu_fe}:-ffp-contract=off>)
endif()


if(${ASTCENC_CLI})
    # Enable LTO on release builds
    set_property(TARGET common
        PROPERTY
            INTERPROCEDURAL_OPTIMIZATION_RELEASE True)
endif()

add_library(simd-none INTERFACE)
target_link_libraries(simd-none INTERFACE common)
target_compile_definitions(simd-none INTERFACE
    ASTCENC_NEON=0
    ASTCENC_SSE=0
    ASTCENC_AVX=0
    ASTCENC_POPCNT=0
    ASTCENC_F16C=0)

# Force SSE2 on AppleClang (normally SSE4.1 is the default)
target_compile_options(simd-none INTERFACE
    $<${is_gnu_fe}:-mno-sse2>
    $<${is_gnu_fe}:-mno-sse4.1>
    $<${is_gnu_fe}:-Wno-unused-command-line-argument>
)


add_library(simd-neon INTERFACE)
target_link_libraries(simd-neon INTERFACE common)
target_compile_definitions(simd-neon INTERFACE
    ASTCENC_NEON=1
    ASTCENC_SSE=0
    ASTCENC_AVX=0
    ASTCENC_POPCNT=0
    ASTCENC_F16C=0)

# target_compile_options(simd-neon INTERFACE
#     $<${is_clang}:-mfpu=neon -mfloat-abi=hard>
#     $<${is_gnu_fe}:-Wno-unused-command-line-argument>)


# Workaround MSVC codegen bug for NEON builds on VS 2022 17.2 or older
# https://developercommunity.visualstudio.com/t/inlining-turns-constant-into-register-operand-for/1394798
if (MSVC_VERSION LESS 1933)
    target_compile_options(simd-neon INTERFACE
        $<${is_msvccl}:/d2ssa-cfg-sink->)
endif()

# SSE 2.0 compiler options
add_library(simd-sse2 INTERFACE)
target_link_libraries(simd-sse2 INTERFACE common)
target_compile_definitions(simd-sse2 INTERFACE
    ASTCENC_NEON=0
    ASTCENC_SSE=20
    ASTCENC_AVX=0
    ASTCENC_POPCNT=0
    ASTCENC_F16C=0)

# Force SSE2 on AppleClang (normally SSE4.1 is the default)
target_compile_options(simd-sse2 INTERFACE
    $<${is_clangcl}:-msse2>
    $<${is_gnu_fe}:-msse2>
    $<${is_gnu_fe}:-mno-sse4.1>
    $<${is_gnu_fe}:-Wno-unused-command-line-argument>)

# SSE 4.1 compiler options
add_library(simd-sse41 INTERFACE)
target_link_libraries(simd-sse41 INTERFACE common)
target_compile_definitions(simd-sse41 INTERFACE
    ASTCENC_NEON=0
    ASTCENC_SSE=41
    ASTCENC_AVX=0
    ASTCENC_POPCNT=1
    ASTCENC_F16C=0)

target_compile_options(simd-sse41 INTERFACE
    $<${is_clangcl}:-msse4.1 -mpopcnt>
    $<${is_gnu_fe}:-msse4.1 -mpopcnt>
    $<${is_gnu_fe}:-Wno-unused-command-line-argument>)

# AVX2 compiler options
add_library(simd-avx2 INTERFACE)
target_link_libraries(simd-avx2 INTERFACE common)
target_compile_definitions(simd-avx2 INTERFACE
    ASTCENC_NEON=0
    ASTCENC_SSE=41
    ASTCENC_AVX=2
    ASTCENC_POPCNT=1
    ASTCENC_F16C=1)

target_compile_options(simd-avx2 INTERFACE
        $<${is_msvc_fe}:/arch:AVX2>
        $<${is_clangcl}:-mavx2 -mpopcnt -mf16c>
        $<${is_gnu_fe}:-mavx2 -mpopcnt -mf16c>
        $<${is_gnu_fe}:-Wno-unused-command-line-argument>)


# Non-invariant builds enable us to loosen the compiler constraints on
# floating point, but this is only worth doing on CPUs with AVX2 because
# this implies we can also enable the FMA instruction set extensions
# which significantly improve performance. Note that this DOES reduce
# image quality by up to 0.2 dB (normally much less), but buys an
# average of 10-15% performance improvement ...
if((NOT ${ASTCENC_INVARIANCE}) AND (NOT ${ASTCENC_IS_VENEER}))
    target_compile_options(simd-avx2 INTERFACE
        $<${is_gnu_fe}:-mfma>)
endif()

macro(TARGET_SIMD_OPTIONS_BY_ARCH SIMD_TARGET_NAME ARCH_NAME)
    if (${ARCH_NAME} STREQUAL "arm64")
        target_link_libraries(${SIMD_TARGET_NAME} PRIVATE simd-neon)
    elseif (${ARCH_NAME} STREQUAL "x86_64h")
        target_link_libraries(${SIMD_TARGET_NAME} PRIVATE simd-avx2)
    elseif(${ARCH_NAME} STREQUAL "x86_64")
        target_link_libraries(${SIMD_TARGET_NAME} PRIVATE simd-sse41)
    endif()
endmacro()

macro(TARGET_SIMD_OPTIONS_BY_NAME SIMD_TARGET_NAME SIMD_NAME)
    if (${SIMD_NAME} STREQUAL "neon")
        target_link_libraries(${SIMD_TARGET_NAME} PRIVATE simd-neon)
    elseif (${SIMD_NAME} STREQUAL "avx2")
        target_link_libraries(${SIMD_TARGET_NAME} PRIVATE simd-avx2)
    elseif(${SIMD_NAME} STREQUAL "sse4.1")
        target_link_libraries(${SIMD_TARGET_NAME} PRIVATE simd-sse41)
    elseif(${SIMD_NAME} STREQUAL "sse2")
        target_link_libraries(${SIMD_TARGET_NAME} PRIVATE simd-sse2)
    elseif(${SIMD_NAME} STREQUAL "none")
        target_link_libraries(${SIMD_TARGET_NAME} PRIVATE common)
    endif()


    if(${ASTCENC_ISA_SIMD} MATCHES "neon")
        set_target_properties(${SIMD_TARGET_NAME} PROPERTIES OSX_ARCHITECTURES arm64)
    elseif(${ASTCENC_ISA_SIMD} MATCHES "avx2")
        set_target_properties(${SIMD_TARGET_NAME} PROPERTIES OSX_ARCHITECTURES x86_64h)
    elseif(NOT ${ASTCENC_ISA_SIMD} MATCHES "none")
        set_target_properties(${SIMD_TARGET_NAME} PROPERTIES OSX_ARCHITECTURES x86_64)
    endif()

endmacro()
