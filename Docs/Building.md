# Building ASTC Encoder

This page provides instructions for building `astcenc` from the sources in
this repository.

**Note:** The current `master` branch is configured to statically build
binaries which each use a specific level of SIMD support (SSE2, SSE4.2,
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

Builds for macOS and Linux use GCC and Make, and are tested with GCC 7.4 and
GNU Make 3.82.

### Single variants

To compile a single SIMD variant compile with:

```
make -sC Source VEC=[nointrin|sse2|sse4.2|avx2] -j8
```

... and use:

```
make -sC Source VEC=[nointrin|sse2|sse4.2|avx2] clean
```

... to clean the build.

### All variants

To compile all supported SIMD variants compile with:

```
make -sC Source batchbuild -j8
```

... and use:

```
make -sC Source batchclean -j8
```

... to clean the build.

