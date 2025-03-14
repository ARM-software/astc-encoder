# Building Big Endian ASTC Encoder

**NOTE:** The current code is _NOT_ BE compatible, but we plan to fix this!

We don't officially support a big-endian compressor build, but ensuring the
code is BE-compatible is good practice and some of the open source
distributions still support BE platforms.

Even though Arm64 can run in a BE mode, it's now very rare in practice. It's no
longer supported out of the box in the latest Arm upstream compiler releases,
and getting hold of a sysroot is increasingly difficult. To test BE builds, I
therefore cross-compile Linux builds for MIPS64 and use `qemu-user` to run
them. This doesn't use a real sysroot, and so everything must be compiled with
`-static` linkage.

## Host software

Install the following host software:

```bash
# Compiler
sudo apt-get install g++-mips64-linux-gnuabi64

# Multi-arch libraries
sudo apt-get install g++-multilib-mips64-linux-gnuabi64

# QEMU
sudo apt-get install qemu-user-static
sudo mkdir /etc/qemu-binfmt
sudo ln -s /usr/mips64-linux-gnuabi64 /etc/qemu-binfmt/mips64
```

## CMake toolchain file

Cross-compiling needs a correctly configured CMake, and the easiest way to
do this consistently is to use a toolchain file. Create a `CMake-BE.toolchain`
file in the root of the project, with the following content:

```
# Operating system
set(CMAKE_SYSTEM_NAME Linux)

# Cross-compilers for C and C++
set(CMAKE_C_COMPILER mips64-linux-gnuabi64-gcc)
set(CMAKE_CXX_COMPILER mips64-linux-gnuabi64-g++)

# Compiler environment
set(CMAKE_FIND_ROOT_PATH /usr/mips64-linux-gnuabi64)

# Default compiler and linker flags to use
set(CMAKE_C_FLAGS_INIT -static)
set(CMAKE_CXX_FLAGS_INIT -static)
set(CMAKE_EXE_LINKER_FLAGS_INIT -static)
set(CMAKE_SHARED_LINKER_FLAGS_INIT -static)
set(CMAKE_MODULE_LINKER_FLAGS_INIT -static)

# Never match host tools
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)

# Only match headers and libraries in the compiler environment
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
```

## Build astcenc

Building uses CMake as normal, with the additional specification of the
toolchain file to configure the build for cross-compilation. We don't have any
SIMD implementations for big-endian architectures so these builds must compile for the reference C SIMD implementation, `ASTCENC_ISA_NONE`.

```
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../ -DASTCENC_ISA_NONE=ON -DCMAKE_TOOLCHAIN_FILE=../CMake-BE.toolchain ..
```

## Run astcenc

The cross-compiled `astcenc` binary runs as normal, and can access host files, but must run through QEMU to do the instruction set translation.

```
qemu-mips64 ./bin/astcenc-none ...
```

- - -

_Copyright © 2025, Arm Limited and contributors._
