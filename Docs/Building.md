# Building ASTC Encoder

This page provides instructions for building `astcenc` from the sources in
this repository.

## Windows

Builds for Windows use Visual Studio 2017, using the solution file located in
the `Source/win32-2017/astcenc/` directory. To compile a release build from
the command line, it is possible to use the `msbuild` utility:

```
msbuild .\Source\win32-2017\astcenc\astcenc.sln /p:Configuration=Release /p:Platform=x64
```

## macOS

Builds for macOS use GCC and Make, and are tested with GCC 4.6 and GNU Make
3.82.

```
cd Source
make
```

## Linux

Builds for Linux use GCC and Make, and are tested with GCC 4.6 and GNU Make
3.82.

```
cd Source
make
```

