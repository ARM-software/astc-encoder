# 6.x series change log

This page summarizes the major functional and performance changes in each
release of the 6.x series.

All performance data on this page is measured on an Intel Xeon w7-3465X clocked
at 3.9 GHz, running `astcenc` using AVX2 and 12 threads.

<!-- ---------------------------------------------------------------------- -->
## 6.0.0

**Status:** In development.

The 6.0.0 release is a major release, adding support for two flavors of AVX-512
and making significant performance optimizations. Changes that break
compatibility with the 5.x library API are also planned.

* **General updates:**
  * **Feature:** Add AVX-512 SIMD backends to library, with corresponding
    command line tool variants. Implementation provides "AVX-256", which uses
    an enhanced version of the AVX2 backend operating on 256-bit vectors, and
    full AVX-512 operating on 512-bit vectors. Note that full AVX-512
    performance can be slower than AVX2 for smaller block sizes.
* **Codec library updates:**
  * **Bug fix:** Remove data races when initializing configurations and
    allocating contexts concurrently.
  * **Bug fix:** Preserve the error threshold inherited by child contexts,
    avoiding a second conversion from dB.
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
* **Build system updates:**
  * **Bug fix:** Restore Clang compression performance by disabling unintended
    floating-point trapping, retaining precise arithmetic and the configured
    FMA contraction policy.

- - -

_Copyright © 2026, Arm Limited and contributors._
