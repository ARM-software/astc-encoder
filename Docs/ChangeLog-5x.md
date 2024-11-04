# 5.x series change log

This page summarizes the major functional and performance changes in each
release of the 5.x series.

All performance data on this page is measured on an Intel Core i5-9600K
clocked at 4.2 GHz, running `astcenc` using AVX2 and 6 threads.


<!-- ---------------------------------------------------------------------- -->
## 5.1.0

**Status:** In development.

The 5.1.0 release is a maintenance release.

* **General:**
  * **Feature:** Added a new CMake build option to control use of native
    gathers, as they can be slower than scalar loads on some common x86
    microarchitectures. Build with `-DASTCENC_X86_GATHERS=OFF` to disable use
    of native gathers in AVX2 builds.
  * **Optimization:** Added new `gather()` abstraction for gathers using byte
    indices, allowing implementations without gather hardware to skip the
    byte-to-int index conversion.
  * **Optimization:** Added improved intrinsics sequence for SSE and AVX2
    `hmin()` and `hmax()`.

<!-- ---------------------------------------------------------------------- -->
## 5.0.0

**Status:** November 2024

The 5.0.0 release is the first stable release in the 5.x series. The main new
feature is support for the Arm Scalable Vector Extensions (SVE) SIMD instruction
set.

* **General:**
  * **Bug fix:** Fixed incorrect return type in "None" vector library
    reference implementation.
  * **Bug fix:** Fixed sincos table index under/overflow.
  * **Feature:** Changed `ASTCENC_ISA_NATIVE` builds to use `-march=native` and
    `-mcpu=native`.
  * **Feature:** Added backend for Arm SVE fixed-width 256-bit builds. These
    can only run on hardware implementing 256-bit SVE.
  * **Feature:** Added backend for Arm SVE 128-bit builds. These are portable
    builds and can run on hardware implementing any SVE vector length, but the
    explicit SVE use is augmented NEON and will only use the bottom 128-bits of
    each SVE vector.
  * **Feature:** Optimized NEON mask `any()` and `all()` functions.
  * **Feature:** Migrated build and test to GitHub Actions pipelines.

- - -

_Copyright © 2022-2024, Arm Limited and contributors. All rights reserved._
