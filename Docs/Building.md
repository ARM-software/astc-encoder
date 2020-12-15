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

# Create the build system
cmake -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./ \
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
nmake install -j16
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

# Create the build system
cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./ \
    -DISA_AVX2=ON -DISA_SSE41=ON -DISA_SSE2=ON ..
```

This example shows all SIMD variants being enabled. It is possible to build a
subset of the supported variants by enabling only the ones you require. At
least one variant must be enabled.

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

All normal builds will use SIMD accelerated code paths using instrinsics, as
x86-64 guarantees availability of at least SSE2. For development purposes it
is possible to build an intrinsic-free build which uses no explicit SIMD
acceleration (the compiler may still auto-vectorize).

To enable this binary variant add `-DISA_NONE=ON` to the CMake command line
when configuring. It is NOT

### ISA Invariance

Normal builds are not ISA invariant, meaning that builds for different
instruction sets on the same CPU hardware can produce subtly different outputs.
This is caused by precision differences between, e.g, FMA and DOT hardware
instructions and their equivalent C code.

To build an ISA invariant build add `-DISA_INVARIANCE=ON` to the CMake command
line when configuring. Note that this will reduce performance, as it will
disable use of hardware instructions that cannot be matched by the reference
functionality available in SSE2.

Note that even with ISA invariance enabled we do not guarantee invariant output
across CPU implementations or compilers. The C specification does not require
bit-exact implementations for many maths library functions (`sin()`, `cos()`,
etc) and there are known precision differences across vendors.

### Build Types

We support and test the following `CMAKE_BUILD_TYPE` options.

| Value            | Description                                              |
| ---------------- | -------------------------------------------------------- |
| Release          | Optimized release build                                  |
| RelWithDebInfo   | Optimized release build with debug info                  |
| Debug            | Unoptimized debug build with debug info                  |

Note that optimized release builds are compiled with link-time optimization,
which can make profiling more challenging ...

### Packaging

We support building a release bundle of all enabled binary configurations in
the current CMake configuration using the `package` build target

```bash
# Run a build and package build outputs in `./astcenc-<ver>-<os>-<arch>.<fmt>`
cd build
make package -j16
```

Windows packages will use the `.zip` format, other packages will use the
`.tar.gz` format.
