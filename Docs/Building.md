# Building ASTC Encoder

This page provides instructions for building `astcenc` from the sources in
this repository. The current `master` branch is configured to statically build
binaries which each use a specific level of SIMD support (SSE2, SSE4.1,
or AVX2) selected at compile time. Binaries are produced with a name postfix
indicating the SIMD type in use; e.g. `astcenc-avx2` for the AVX2 binary.

## Windows

Builds for Windows use Visual Studio 2019, using the solution and/or project
files located in the `Source/VS2019/` directory.

### Single variants

To compile a single SIMD variant from the command line, it is possible to use
the `msbuild` utility from the Visual Studio installation, targeting the
specific project for the variant desired.

```
msbuild astcenc-<variant>.vcxproj /p:Configuration=Release /p:Platform=x64
```

### Multiple variants

To compile all supported SIMD variants from the command line, it is possible
to use the `msbuild` utility from the Visual Studio installation, targeting the
solution.

```
msbuild astcenc.sln /p:Configuration=Release /p:Platform=x64
```

## macOS and Linux

Builds for macOS and Linux use GCC or Clang and Make. They are tested using
Clang 9.0, GCC 7.4, and Make 3.82. Using Clang 9.0 is recommended, as it
out-performs GCC by 15-20% in benchmarked test runs.

### Single variants

To compile a single SIMD variant compile with:

```
make -C Source CXX=clang++ VEC=[sse2|sse4.1|avx2] -j8
```

... and use:

```
make -C Source CXX=clang++ VEC=[sse2|sse4.1|avx2] clean
```

... to clean the build.

### All variants

To compile all supported SIMD variants compile with:

```
make -sC Source CXX=clang++ batchbuild -j8
```

... and use:

```
make -sC Source batchclean -j8
```

... to clean the build.

### Developer build options

The default build, `BUILD=release`, is heavily optimized using link-time
optimization (`-O3 -flto`) and does not include symbols.

For debugging, add `BUILD=debug` to the Make command line. This will build a
non-optimized build (`-O0`) with symbols included (`-g`).

For easier profiling, add `BUILD=profile` to the Make command line. This will
build a moderately optimized build (`-O2`) with symbols include (`-g`), but
without link-time optimization (no `-flto`).

Normal builds are not ISA invariant, builds for different instruction sets on
the same CPU hardware can produce subtly different outputs. To build an ISA
invariant build set `-DASTCENC_ISA_INVARIANCE=1`. For make builds this can be
achieved by setting `ISA_INV=1` on the command line. This will reduce
performance, as optimizations will be disabled to keep alignment with SSE2.
