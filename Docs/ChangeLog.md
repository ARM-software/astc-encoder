# 2.x series change log

This page summarizes the major functional and performance changes in each
release of the 2.x series.

All performance data on this page is measured on an Intel Core i5-9600K
clocked at 4.2 GHz, running astcenc using 6 threads.

<!-- ---------------------------------------------------------------------- -->
## 2.2

**Status:** :warning: In development (ETA February 2021)

The 2.2 release is the third release in the 2.x series. It includes ...

Reminder for users of the library interface - the API is not designed to be
stable across versions, and this release is not compatible with 2.1. Please
recompile your client-side code using the updated `astcenc.h` header.

* **General:**
  * **Improvement:** SSE4.2 feature profile changed to SSE4.1, which more
    accurately reflects the feature set used.
  * **Improvement:** Build system changed to use CMake for all platforms.
* **Binary releases:**
  * **Improvement:** Linux binaries changed to use use Clang 9.0, which gives
    up to 15% performance improvement.
* **Command Line:**
  * **Feature:** New image preprocess `-pp-normalize` option added. This forces
    normal vectors to be unit length, which is useful when compressing source
    textures that use normal length to encode an NDF, which is incompatible
    with ASTC's two channel encoding.
  * **Feature:** New image preprocess `-pp-premultiply` option added. This
    scales RGB values by the alpha value. This can be useful to minimize
    cross-channel color bleed caused by GPU post-multiply filtering/blending.
  * **Improvements:** Command line tool cleanly traps and reports errors for
    corrupt input images rather than relying on hard standard library
    `assert()` calls.
* **Core API:**
  * **API Change:** Images using region-based metrics no longer need to include
    padding; all input images should be tightly packed and `dim_pad` is removed
    from the `astcenc_image` structure. This makes it easier to directly use
    images loaded from other libraries.
  * **API Change:** Image `data` is no longer a 3D array accessed using
    `data[z][y][x]` indexing, it's an array of 2D slices. This makes it easier
    to directly use images loaded from other libraries.
  * **API Change:** New `ASTCENC_FLG_SELF_DECOMPRESS_ONLY` flag added to the
    codec config. Using this flag enables additional optimizations that
    aggressively exploit implementation- and configuration-specific, behavior
	to gain performance. When using this flag the codec can only reliably
	decompress images that were compressed in the same context session. Images
	produced via other means may fail to decompress correctly, even if they are
	otherwise valid ASTC files.

<!-- ---------------------------------------------------------------------- -->
## 2.1

**Status:** Released, November 2020

The 2.1 release is the second release in the 2.x series. It includes a number
of performance optimizations and new features.

Reminder for users of the library interface - the API is not designed to be
stable across versions, and this release is not compatible with 2.0. Please
recompile your client-side code using the updated `astcenc.h` header.

### Features:

* **Command line:**
  * **Bug fix:** The meaning of the `-tH\cH\dH` and `-th\ch\dh` compression
    modes was inverted. They now match the documentation; use `-*H` for HDR
    RGBA, and `-*h` for HDR RGB with LDR alpha.
  * **Feature:** A new `-fastest` quality preset is now available. This is
    designed for fast "roughing out" of new content, and sacrifices significant
    image quality compared to `-fast`. We do not recommend its use for
    production builds.
  * **Feature:** A new `-candidatelimit` compression tuning option is now
    available. This is a power-user control to determine how many candidates
    are returned for each block mode encoding trial. This feature is used
	automatically by the search presets; see `-help` for details.
  * **Improvement:** The compression test modes (`-tl\ts\th\tH`) now emit a
    MTex/s performance metric, in addition to coding time.
* **Core API:**
  * **Feature:** A new quality preset `ASTCENC_PRE_FASTEST` is available. See
    `-fastest` above for details.
  * **Feature:** A new tuning option `tune_candidate_limit` is available in
    the config structure. See `-candidatelimit` above for details.
  * **Feature:** Image input/output can now use `ASTCENC_TYPE_F32` data types.
* **Stability:**
  * **Feature:** The SSE2, SSE4.2, and AVX2 variants now produce identical
    compressed output when run on the same CPU when compiled with the
    preprocessor define `ASTCENC_ISA_INVARIANCE=1`. For Make builds this can
    be set on the command line by setting `ISA_INV=1`. ISA invariance is off
    by default; it reduces performance by 1-3%.

### Performance

Key for performance charts:

* Color = block size (see legend).
* Letter = image format (N = normal map, G = greyscale, L = LDR, H = HDR).

**Absolute performance vs 2.0 release:**

![Absolute scores 2.1 vs 2.0](./ChangeLogImg/absolute-2.0-to-2.1.png)

**Relative performance vs 2.0 release:**

![Relative scores 2.1 vs 2.0](./ChangeLogImg/relative-2.0-to-2.1.png)

<!-- ---------------------------------------------------------------------- -->
## 2.0

**Status:** Released, August 2020

The 2.0 release is first release in the 2.x series. It includes a number of
major changes over the earlier 1.7 series, and is not command-line compatible.

### Features:

* The core codec can be built as a library, exposed via a new codec API.
* The core codec supports accelerated SIMD paths for SSE2, SSE4.2, and AVX2.
* The command line syntax has a clearer mapping to Khronos feature profiles.

### Performance:

Key for performance charts

* Color = block size (see legend).
* Letter = image format (N = normal map, G = greyscale, L = LDR, H = HDR).

**Absolute performance vs 1.7 release:**

![Absolute scores 2.0 vs 1.7](./ChangeLogImg/absolute-1.7-to-2.0.png)

**Relative performance vs 1.7 release:**

![Relative scores 2.0 vs 1.7](./ChangeLogImg/relative-1.7-to-2.0.png)
