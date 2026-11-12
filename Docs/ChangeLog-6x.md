# 6.x series change log

This page summarizes the major functional and performance changes in each
release of the 6.x series.

All performance data on this page is measured on an Intel Xeon w7-3465X clocked
at 3.9 GHz, running `astcenc` using AVX2 and 12 threads.

<!-- ---------------------------------------------------------------------- -->
## 6.0.0

**Status:** In development.

The 6.0.0 release is a major release, adding support for two flavors of AVX-512
and making significant performance optimizations to all of the ISA variants.

This release includes changes to the public interface in the `astcenc.h`
header. We always recommend rebuilding your client-side code using the
header from the same release to avoid compatibility issues.

* **Major API and feature updates:**
  * **API Change:** Add a row stride to `astcenc_image` for compression input
    and decompression output. A zero stride selects tightly packed rows.
  * **API Change:** Use `size_t` for image dimensions, context thread counts,
    and thread indices in the library API.
  * **API Change:** Decompressing images containing invalid encodings, or HDR
    encodings in an LDR profile, now return an identifiable
    `ASTCENC_ERR_BAD_BLOCK` return code. The image will still be decompressed,
    using the specified error color for the format in use for the bad blocks.
    The command line tool will write the image, emit a warning, and exit with
    code 2 for these failures. Functional errors continue to use exit code 1.
  * **Feature:** Add AVX-512 SIMD backends to library, with corresponding
    command line tool variants. Implementation provides "AVX-256", which uses
    an enhanced version of the AVX2 backend operating on 256-bit vectors, and
    full AVX-512 operating on 512-bit vectors. Note that full AVX-512
    performance can be slower than AVX2 for smaller block sizes.
  * **Change:** Require C++20 to build the codec, replacing local library
    functions with new standard library operations where possible.
  * **Change:** Minimum supported compiler versions increased to Clang 18, GCC
    13, and Visual Studio 2022. Minimum supported CMake version increased to
    3.21.
* **Codec library updates:**
  * **Bug fix:** Preserve HDR endpoint mode bits when quantization rounds down
    across a mode bit boundary.
  * **Bug fix:** Wait for outstanding tasks to finish before releasing shared
    buffers when compression is cancelled.
  * **Bug fix:** Remove a data race when reporting compression progress.
  * **Bug fix:** Remove data races when initializing configurations and
    allocating contexts concurrently.
  * **Bug fix:** Preserve the error threshold inherited by child contexts,
    avoiding a second conversion from dB.
  * **Image quality:** Reduce visible HDR artifacts with perceptual mode
    enabled by weighting logarithmic channel errors according to geomean
    of channel brightness.
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
  * **Optimization:** Reuse RGB partition metrics from grayscale RGBA statistics
    during single-plane endpoint selection, improving normal map compression
    performance.
  * **Optimization:** Skip RGB+lumashift endpoint error calculations for LDR
    blocks, as this encoding is only available for HDR.
  * **Optimization:** Interpolate deltas in `lerp_color_int()` and unrolled
    specialization, saving one multiply.
  * **Optimization:** Hoist loop invariant error accumulator calculation out of
    deeper nested `..._find_best_combination...()` selector loops.
  * **Optimization:** Directly copy non-decimated weight grids when unpacking
    symbolic blocks, avoiding redundant interpolation.
  * **Optimization:** Specialize `compute_avgs_and_dirs_4_comp()` for blocks
    with identical `RGB` components, skipping duplicate `GB` calculations.
    This significantly improves normal map compression performance.
  * **Optimization:** Use SoA-style SIMD batches for k-means cluster assignment.
  * **Optimization:** Batch partition candidate ranking for RGB and RRRA,
    improving performance for wider SIMD vector hardware.
  * **Optimization:** Specialize `compute_partition_averages_rgba()` for blocks
    with identical `RGB` components, skipping duplicate `GB` calculations.
    This significantly improves normal map compression performance.
  * **Optimization:** Specialize RGBA partition-error scoring for blocks with
    identical `RGB` components, reusing channel gathers, projection products,
    and residual calculations while preserving channel weights and error sums.
  * **Optimization:** Specialize endpoint refitting for LDR blocks with
    identical `RGB` components and zero requested `GB` weights, fitting only `R`
    and `A`. This improves normal map compression performance.
  * **Optimization:** Avoid constructing three- and four-partition
    endpoint-format combinations that exceed the allowed endpoint integer limit.
  * **Optimization:** Split multi-partition endpoint-format searches by
    quantization level to skip unreachable high-precision combinations.
  * **Optimization:** Use NEON combined shift-and-insert for byte packing in
    `interleave_rgba8()`.
  * **Optimization:** Use SVE narrowing stores in the 4-wide implementation of
    `pack_and_store_low_bytes()`.
* **Command line tool updates:**
  * **Feature:** Report SSIM for LDR images and exposure-averaged mSSIM for
    HDR images in test mode. HDR exposures follow the `-mpsnr` range.
  * **Change:** Automatically center HDR image quality metric exposures on the
    original image's geometric mean luminance, covering two stops either side.
    An explicit `-mpsnr` range overrides automatic exposure selection.
  * **Feature:** Add crop region support for compression. Crop regions are
    extracted using image coordinate space before flipping or preprocessing.
* **Build system updates:**
  * **Bug fix:** Restore Clang compression performance by disabling unintended
    floating-point trapping, retaining precise arithmetic and the configured
    FMA contraction policy.

- - -

_Copyright © 2026, Arm Limited and contributors._
