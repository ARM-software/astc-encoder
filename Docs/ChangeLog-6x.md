# 6.x series change log

This page summarizes the major functional and performance changes in each
release of the 6.x series.

All performance data on this page is measured on an Intel Xeon w7-3465X clocked
at 3.9 GHz, running `astcenc` using AVX2 and 12 threads.

<!-- ---------------------------------------------------------------------- -->
## 6.0.0

**Status:** In development.

The 6.0.0 release is a major release, adding support for Rate-Distortion
Optimization (RDO) to both CLI and library, and support for AVX-512 as a
library backend. There is also significant performance uplifit across all
backends, and numerous small improvements to image quality and stability.

Thank you to Epic Games for the contribution of the LZ RDO implementation, and
to AMD for the contribution of the AVX-512 backend.

This release includes changes to the public interface in the `astcenc.h`
header. We always recommend rebuilding your client-side code using the
header from the same release to avoid compatibility issues.

* **Major API and feature updates:**
  * **API Change:** Add a row stride to `astcenc_image` for compression input
    and decompression output. A zero stride selects tightly packed rows.
  * **API Change:** Use `size_t` for image dimensions, context thread counts,
    and thread indices in the library API.
  * **API Change:** Add support for rate-distortion optimization to reduce
    size of bundled ASTC texture assets. Enable by setting the config flag
    `ASTCENC_FLG_USE_LZ_RDO` and setting the appropriate `lz_rdo_lambda`
    level (0 = no RDO, 20 = light RDO, 40 = moderate RDO ...).
  * **API Change:** Decompressing images containing invalid encodings, or HDR
    encodings in an LDR profile, now return an identifiable
    `ASTCENC_ERR_BAD_BLOCK` return code. The image will still be decompressed,
    using the specified block error color for the ASTC profile in use. The
    command line tool will write the image, emit a warning, and exit with code
    2 for these failures. Functional errors continue to use exit code 1.
  * **CLI Change:** Add support for rate-distortion optimization to reduce
    size of bundled ASTC texture assets. Enable by passing the `-rdo-lambda`
    command line option. See API Change above for details.
  * **Feature:** Add an AVX-512 SIMD backend to the library, with a
    corresponding command line tool variant operating on 512-bit vectors.
    AVX-512 performance can be slower than AVX2 for smaller block sizes.
  * **Build Change:** AVX2 builds now also require BMI and BMI2.
    functions with new standard library operations where possible.
  * **Build Change:** Require C++20 to build the codec, replacing local library
    functions with new standard library operations where possible.
  * **Build Change:** Minimum supported compiler versions increased to Clang
    18, GCC 13, and Visual Studio 2022. Minimum supported CMake version
    increased to 3.21.
* **Codec library updates:**
  * **Bug fix:** Disable dual-plane encoding at a zero correlation cutoff.
  * **Bug fix:** Preserve HDR endpoint mode bits when quantization rounds down
    across a mode bit boundary.
  * **Bug fix:** Wait for outstanding tasks to finish before releasing shared
    buffers when compression is cancelled.
  * **Bug fix:** Remove a data race when reporting compression progress.
  * **Bug fix:** Remove incorrectly shared data that could cause concurrently
    initializing contexts to use an incorrect configuration.
  * **Bug fix:** Preserve the error threshold inherited by child contexts,
    avoiding a second conversion from dB.
  * **Image quality:** Reduce visible HDR artifacts with perceptual mode
    enabled by weighting logarithmic channel errors according to geomean
    of channel brightness.
  * **Optimization:** Restructuture partition candidate selection by processing
    multiple candidates in parallel on wider vector hardware, one candidate for
    every 128-bits in the available SIMD vectors. This significantly improves
    performance on 256-bit and 512-bit vector hardware.
  * **Optimization:** Specialize a simpler partition candidate selection for
    luminance-only blocks, as we know that direction is always through the
    origin.
  * **Optimization:** Specialize many data paths for luminance+alpha data,
    which is commonly used for normal maps.
  * **Optimization:** Use analytical calculation of weight change error diffs
    in `realign_weights_decimated()` and `realign_weights_undecimated()`.
    Changes compute simpler linear and quadratic error coefficients, and then
    compute the up and down candidate deltas.
  * **Optimization:** Precompute an infilled weight set for the realign loop
    in `realign_weights_decimated()`, rebuilding set entries if a decimated
    weight changes, rather than always recomputing the infilled value every
    loop iteration.
  * **Optimization:** Reuse RGB partition metrics computed during ideal weight
    generation when selecting endpoint formats, avoiding a duplicate texel
    scan.
  * **Optimization:** Skip RGB+lumashift endpoint error calculations for LDR
    blocks, as this encoding is only available for HDR.
  * **Optimization:** Interpolate deltas in `lerp_color_int()` and unrolled
    specialization, saving one multiply.
  * **Optimization:** Directly copy non-decimated weight grids when unpacking
    symbolic blocks, avoiding redundant interpolation.
  * **Optimization:** Use SoA-style SIMD batches for k-means cluster assignment.
  * **Optimization:** Avoid constructing three- and four-partition
    endpoint-format combinations that exceed the allowed endpoint integer limit.
  * **Optimization:** Split multi-partition endpoint-format searches by
    quantization level to skip unreachable high-precision combinations.
  * **Optimization:** Specialize ISE decoding at compile-time into functions
    per quant level to help with RDO performance.
  * **Optimization:** Optimized implementaiton of unorm16_to_sf16 to help with
    RDO performance.
  * **Optimization:** Use NEON combined shift-and-insert for byte packing in
    `interleave_rgba8()`.
  * **Optimization:** Use SVE narrowing stores in the 4-wide implementation of
    `pack_and_store_low_bytes()`.
* **Command line tool updates:**
  * **CLI Change:** Writing to an LDR output file format, such as PNG, when
    test compressing an HDR image is now supported. The tool will output
    multiple files, one for each active integer f-stop exposure level.
  * **Feature:** Report SSIM for LDR images and exposure-averaged mSSIM for
    HDR images in test mode. HDR exposures follow the `-mpsnr` range.
  * **Feature:** Report DEFLATE compression when using RDO. Note - requires
    CLI to be built with zlib and is not enabled in our release binaries.
  * **Change:** Automatically center HDR image quality metric exposures on the
    original image's geometric mean luminance, covering two stops either side.
    An explicit `-mpsnr` range overrides automatic exposure selection.
  * **Feature:** Add crop region support for compression. Crop regions are
    extracted in the image coordinate space before flipping or preprocessing.
* **Build system updates:**
  * **Bug fix:** Restore Clang compression performance by disabling unintended
    floating-point trapping, retaining precise arithmetic and the configured
    FMA contraction policy.

- - -

_Copyright © 2026, Arm Limited and contributors._
