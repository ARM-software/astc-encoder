# Building ASTC Encoder

This page provides instructions for building `astcenc` from the sources in
this repository.

Builds use CMake 3.15 or higher as the build system generator. The examples on
this page only show how to use it to target NMake (Windows) and Make
(Linux and macOS), but CMake supports other build system backends.

## Windows

Builds for Windows are tested with CMake 3.17 and Visual Studio 2019.

### Configuring the build

To use CMake you must first configure the build. Create a build directory
in the root of the astenc checkout, and then run `cmake` inside that directory
to generate the build system.

```shell
# Create a build directory
mkdir build
cd build

# Configure your build of choice, for example:

# x86-64 using NMake
cmake -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=.\ ^
    -DISA_AVX2=ON -DISA_SSE41=ON -DISA_SSE2=ON ..

# x86-64 using Visual Studio solution
cmake -G "Visual Studio 16 2019" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=.\ ^
    -DISA_AVX2=ON -DISA_SSE41=ON -DISA_SSE2=ON ..

```

This example shows all SIMD variants being enabled. It is possible to build a
subset of the supported variants by enabling only the ones you require. At
least one variant must be enabled.

### Building

Once you have configured the build you can use NMake to compile the project
from your build dir, and install to your target install directory.

```shell
# Run a build and install build outputs in `${CMAKE_INSTALL_PREFIX}/astcenc/`
cd build
nmake install
```

## macOS and Linux

Builds for macOS and Linux are tested with CMake 3.17 and clang++ 9.0.

> Compiling using g++ is supported, but clang++ builds are faster by ~15%.

### Configuring the build

To use CMake you must first configure the build. Create a build directory
in the root of the astenc checkout, and then run `cmake` inside that directory
to generate the build system.

```shell
# Select your compiler (clang++ recommended, but g++ works)
export CXX=clang++

# Create a build directory
mkdir build
cd build

# Configure your build of choice, for example:

# Arm arch64
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./ \
    -DARCH=aarch64 -DISA_NEON=ON ..

# x86-64
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./ \
    -DISA_AVX2=ON -DISA_SSE41=ON -DISA_SSE2=ON ..

# Host builds optimized for the local CPU's microarchitecture
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./ \
    -DISA_NATIVE=ON ..
```

This example shows all SIMD variants being enabled. It is possible to build a
subset of the supported variants by enabling only the ones you require. If no
variant is explicitly specified the build will compile a native binary for the
build machine.

:warning: Compiling a native binary for the build machine usually gives the
fastest binaries for that specific CPU type but may produce binaries that do
not work (due to using ISA extensions) or run slowly (due to poor instruction
scheduling) on other processors.

### Building

Once you have configured the build you can use Make to compile the project from
your build dir, and install to your target install directory.

```shell
# Run a build and install build outputs in `${CMAKE_INSTALL_PREFIX}/astcenc/`
cd build
make install -j16
```

## Advanced build options

For codec developers there are a number of useful features in the build system.

### No intrinsics build

All normal builds will use SIMD accelerated code paths using intrinsics, as all
target architectures (x86-64 and aarch64) guarantee SIMD availability. For
development purposes it is possible to build an intrinsic-free build which uses
no explicit SIMD acceleration (the compiler may still auto-vectorize).

To enable this binary variant add `-DISA_NONE=ON` to the CMake command line
when configuring. It is NOT recommended to use this for production; it is
significantly slower than the vectorized SIMD builds.

### Build Types

We support and test the following `CMAKE_BUILD_TYPE` options.

| Value            | Description                                              |
| ---------------- | -------------------------------------------------------- |
| Release          | Optimized release build                                  |
| RelWithDebInfo   | Optimized release build with debug info                  |
| Debug            | Unoptimized debug build with debug info                  |

Note that optimized release builds are compiled with link-time optimization,
which can make profiling more challenging ...

### Testing

We support building unit tests.

These builds use the `googletest` framework, which is pulled in though a git
submodule. On first use, you must fetch the submodule dependency:

```shell
git submodule init
git submodule update
```

To build unit tests add `-DUNITTEST=ON` to the CMake command line when
configuring.

To run unit tests use the CMake `ctest` utility from your build directory after
you have built the tests.

```shell
cd build
ctest --verbose
```

### Packaging

We support building a release bundle of all enabled binary configurations in
the current CMake configuration using the `package` build target

```shell
# Run a build and package build outputs in `./astcenc-<ver>-<os>-<arch>.<fmt>`
cd build
make package -j16
```

Windows packages will use the `.zip` format, other packages will use the
`.tar.gz` format.
