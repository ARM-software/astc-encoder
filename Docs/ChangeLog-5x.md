# 5.x series change log

This page summarizes the major functional and performance changes in each
release of the 5.x series.

All performance data on this page is measured on an Intel Xeon w7-3465X clocked
at 3.9 GHz, running `astcenc` using AVX2 and 12 threads.

<!-- ---------------------------------------------------------------------- -->
## 5.8.0

**Status:** In development.

The 5.8.0 release is a minor release, focusing on performance optimization.

* **Codec library updates:**
  * **Optimization:** Hoist loop invariant weight diff calculation out of the
    `compute_ideal_weights_for_decimation()` refinement loop.
  * **Optimization:** Precompute an infilled weight set for the realign loop
    in `realign_weights_decimated()`, rebuilding set entries if a decimated
    weight changes, rather than always recomputing the infilled value every
    loop iteration.
  * **Optimization:** Use analytical calculation of weight change error diffs
    in `realign_weights_decimated()` and `realign_weights_undecimated()`.
    Changes compute simpler linear and quadratic error coefficients, and then
    compute the up and down candidate deltas.
  * **Optimization:** Hoist loop invariant channel weight out of loop in
    `realign_weights_decimated()`.
  * **Optimization:** Use analytical calculation of moments in
    `recompute_ideal_colors_1plane()` after the loop, reducing the amount of
    calculation in the loop.
  * **Optimization:** Reuse RGB partition metrics computed during ideal weight
    generation when selecting endpoint formats, avoiding a duplicate texel
    scan.
  * **Optimization:** Interpolate deltas in `lerp_color_int()` and unrolled
    specialization, saving one multiply.
  * **Optimization:** Hoist loop invariant error accumulator calculation out of
    deeper nested `..._find_best_combination...()` selector loops.
  * **Optimization:** Directly copy non-decimated weight grids when unpacking
    symbolic blocks, avoiding redundant interpolation.
  * **Optimization:** Specialize `compute_avgs_and_dirs_4_comp()` for blocks
    with identical `RGB` components, skipping duplicate `GB` calculations.
    This significantly improves normal map compression performance.
  * **Optimization:** Specialize `compute_partition_averages_rgba()` for blocks
    with identical `RGB` components, skipping duplicate `GB` calculations.
    This significantly improves normal map compression performance.
  * **Optimization:** Specialize RGBA partition-error scoring for blocks with
    identical `RGB` components, reusing channel gathers, projection products,
    and residual calculations while preserving channel weights and error sums.

<!-- ---------------------------------------------------------------------- -->
## 5.7.1

**Status:** Released November 2026.

The 5.7.1 release is a minor maintenance release that addresses a small number
of bugs. These issues are unlikely to affect normal usage, but the fixes
improve robustness and reliability.

* **Codec library updates:**
  * **Bug fix:** Capped alpha-scale (`-a`) option radius to 64 texels to
    avoid very large intermediate buffer memory allocation.
  * **Bug fix:** Avoid UB caused by SIMD code accessing uninitialized
    entries in the `image_block` input data arrays.
  * **Bug fix:** Avoid `float16_to_float()` decoding negative values as NaN in
    builds using x86-64 F16C instructions.
  * **Bug fix:** Avoid div-by-zero error by normalizing `block_z` inside
    `astcenc_context_alloc()` in cases where user manually zeroed the dimension
    after calling `astcenc_config_init()`.
* **Command line tool updates:**
  * **Bug fix:** Improved command line option parameter validation to more
    reliably detect malformed or out-of-range values and emit an error.

<!-- ---------------------------------------------------------------------- -->
## 5.7.0

**Status:** Released August 2026.

The 5.7.0 release is a minor maintenance release.

* **Command line tool updates:**
  * **Bug fix:** Avoid pixel count overflow when loading input PNG file.
  * **Bug fix:** Error if payload is undersized for stated image dimensions
    when loading compressed input KTX file.
  * **Improvement:** Top-level `-help` text as been rewritten to improve
    clarity and consistency.
  * **Change:** Removed `-repeats` benchmarking option, because we prefer to
    repeat by invoking the command line multiple times. This gives performance
    results closer to what a real user sees.
  * **Change:** Removed `-dimage` diagnostic images generation option from
    normal builds. Now only available for builds configured with
    `ASTCENC_DIAGNOSTICS=ON`.

<!-- ---------------------------------------------------------------------- -->
## 5.6.0

**Status:** Released July 2026

The 5.6.0 release is a minor maintenance release.

* **Codec library updates:**
  * **Bug fix:** Avoid undefined behavior caused by passing floating point
    values outside of the [0.0, 1.0] range as input data for a UNORM color
    channel.
  * **Bug fix:** Avoid undefined behavior caused by unaligned access in SSE,
    AVX2, and NEON SIMD library implementations.
* **Command line tool updates:**
  * **Bug fix:** Fixed incorrect plane stride when writing an uncompressed 3D
    LDR image to a DDS container.
  * **Bug fix:** Fixed potential integer overflow when storing very large
    uncompressed images to a `dds` or `.ktx` output image format.
  * **Bug fix:** Fixed missing initialization of the `reserved2` field in the
    DDS file header when writing uncompressed `.dds` output images.
  * **Bug fix:** Fixed typos and grammar issues in the command line `-help`
    text.

<!-- ---------------------------------------------------------------------- -->
## 5.5.0

**Status:** Released June 2026

The 5.5.0 release is a minor maintenance release, fixing minor functional
issues and improving robustness when processing invalid images.

* **Codec library updates:**
  * **Bug fix:** Correct `vfloat8` NaN inequality comparisons for AVX2 builds,
    so that decompressing an invalid encoding returns the correct max magenta
    error color for UNORM output formats instead of black.
  * **Bug fix:** Add overflow detection so that the API returns errors for
    cases when input image size would result in integer overflow in
    calculations of total texel count, total block count, or compressed image
    byte size.
  * **Bug fix:** Fix bad offset calculation in `compute_averages()` when
    using 3D block sizes.
  * **Bug fix:** Add missing low-clamp when storing decompressed HDR data into
    a UNORM8 image. Prior to this fix, output colors would be incorrect if a
    HDR void-extent block contained a negative FP16 constant color component.
  * **Improvement:** Return standard error color when decompression requires an
    unavailable partitioning when using an `ASTCENC_FLG_SELF_DECOMPRESS_ONLY`
    context. Note that this is still not allowed API usage and may trigger
    asserts in a debug build.
* **Command line tool updates:**
  * **Update:** Update stb_image to v2.30.
  * **Update:** Update Wuffs to v0.3.4.
  * **Update:** Update TinyEXR to v1.0.13.
  * **Bug fix:** Use C++ RAII to manage resource lifetime, removing resource
    leaks when an error condition is encountered.
  * **Bug fix:** Report errors for `.astc`, `.dds`, and `.ktx` input image
    sizes that would result in integer overflow in calculations of total block
    count or data size.
  * **Bug fix:** Avoid undefined behavior on NaN to int conversion triggered by
    `atan2()` calls when both inputs are zero.
  * **Bug fix:** Avoid undefined behavior on left shift of negative integer,
    that is triggered for some HDR encodings in `hdr_alpha_unpack()`.
  * **Bug fix:** Redirect error logging for failed `.dds` and `.ktx` input
    image loading to stderr, and improve message consistency.

<!-- ---------------------------------------------------------------------- -->
## 5.4.0

**Status:** Released May 2026

The 5.4.0 release is a minor feature release.

This release includes changes to the public interface in the `astcenc.h`
header. We always recommend rebuilding your client-side code using the
header from the same release to avoid compatibility issues.

* **General:**
  * **Improvement:** The interface header `astcenc.h` is now C compliant to
    make it usable from C programs.
  * **Improvement:** Contexts using the same configuration can now share
    read-only data tables. This can significantly reduce the amount of memory
    needed for applications that parallelize by processing multiple images
    in parallel instead of slicing a single image in parallel.
  * **Improvement:** Decompressor (`astcdec`) builds, which lack compression
    support, now use a smaller `block_size_descriptor` by omitting fields that
    are only needed for compression. This reduces the size of a decompressor
    context by more than 10MB!
  * **Optimization:** A SIMD backend for the RISC-V Vector extensions has been
    added, and is auto-selected when configured with `ASTCENC_ISA_NONE` and a
    core with a 256-bit vector width. See [Building.md](Building.md) for
    details.
  * **Bug fix:** Avoid using an undefined `quant_weight` value if all one
    partition trials return an error block.
  * **Bug fix:** Remove remaining instances of type aliasing through unions.
  * **Bug fix:** Avoid compiler double definition warning for `NOMINMAX` when
    compiling with MinGW.
  * **Bug fix:** Avoid compiler floating point model override warning when
    compiling with Clang 20.

<!-- ---------------------------------------------------------------------- -->
## 5.3.0

**Status:** Released March 2025

The 5.3.0 release is a minor maintenance release.

* **General:**
  * **Feature:** Reference C builds (`ASTCENC_ISA_NONE`) now support compiling
    for big-endian CPUs. Compile with `-DASTCENC_BIG_ENDIAN=ON` when compiling
    for a big-endian target; it is not auto-detected.
  * **Improvement:** Builds using GCC now specify `-flto=auto` to allow
    parallel link steps, and remove the log warnings about not setting a CPU
    count parameter value.
  * **Bug fix:** Builds using MSVC `cl.exe` that do not specify an explicit
    ISA using the preprocessor configuration defines will now correctly
    default to the SSE2 backend on x86-64 and the NEON backend on Arm64.
    Previously they would have defaulted to the reference C implementation,
    which is around 3.25 times slower.


<!-- ---------------------------------------------------------------------- -->
## 5.2.0

**Status:** Released February 2025

The 5.2.0 release is a minor maintenance release.

This release includes changes to the public interface in the `astcenc.h`
header. We always recommend rebuilding your client-side code using the
header from the same release to avoid compatibility issues.

* **General:**
  * **Change:** Changed sRGB alpha channel endpoint expansion to match the
    revised Khronos Data Format Specification (v1.4.0), which reverts an
    unintended specification change. Compared to previous releases, this change
    can cause LSB bit differences in the alpha channel of compressed images.
  * **Feature:** Arm64 builds for Linux added to the GitHub Actions builds, and
    Arm64 binaries for NEON, 128-bit SVE and 256-bit SVE added to release
    builds.
  * **Feature:** Added a new codec API, `astcenc_compress_cancel()`, which can
    be used to cancel an in-flight compression. This is designed to help make
    it easier to integrate the codec into an interactive user interface that
    can respond to user events with low latency.
  * **Bug fix:** Removed incorrect `static` variable qualifier, which could
    result in an incorrect `tune_mse_overshoot` heuristic threshold being used
    if a user ran multiple concurrent compressions with different settings.

<!-- ---------------------------------------------------------------------- -->
## 5.1.0

**Status:** Released November 2024

The 5.1.0 release is an optimization release, giving moderate performance
improvements on all platforms. There are no image quality differences.

* **General:**
  * **Feature:** Added a new CMake build option to control use of native
    gathers, as they can be slower than scalar loads on some common x86
    microarchitectures. Build with `-DASTCENC_X86_GATHERS=OFF` to disable use
    of native gathers in AVX2 builds.
  * **Optimization:** Added new `gather()` abstraction for gathers using byte
    indices, allowing implementations without gather hardware to skip the
    byte-to-int index conversion.
  * **Optimization:** Optimized `compute_lowest_and_highest_weight()` to
    pre-compute min/max outside of the main loop.
  * **Optimization:** Added improved intrinsics sequence for SSE and AVX2
    integer `hmin()` and `hmax()`.
  * **Optimization:** Added improved intrinsics sequence for `vint4(uint8_t*)`
    on systems implementing Arm SVE.

<!-- ---------------------------------------------------------------------- -->
## 5.0.0

**Status:** Released November 2024

The 5.0.0 release is the first stable release in the 5.x series. The main new
feature is support for the Arm Scalable Vector Extensions (SVE) SIMD
instruction set.

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
    explicit SVE use is augmented NEON and will only use the bottom 128 bits of
    each SVE vector.
  * **Feature:** Optimized NEON mask `any()` and `all()` functions.
  * **Feature:** Migrated build and test to GitHub Actions pipelines.

- - -

_Copyright © 2022-2026, Arm Limited and contributors._
