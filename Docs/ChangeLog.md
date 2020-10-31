# 2.x series change log

This page summarizes the major changes in each release of the 2.x series.

## 2.1

**Status:** :warning: In development (ETA November 2020)

The 2.1 release is the second release in the 2.x series. It includes a number
of performance optimizations and new features.
### Features:

* **Command line:**
  * **Bug fix:** The meaning of the `-tH\cH\dH` and `-th\ch\dh` compression
    modes was inverted. They now match the documentation; use `-*H` for HDR
    RGBA, and `-*h` for HDR RGB with LDR alpha.
  * **Feature:** A new `-fastest` quality preset is now available. This is
    designed for fast "roughing out" of new content, and sacrifices significant
    image quality compared to `-fast`. We do not recommend it's use for
    production builds.
  * **Feature:** A new `-candidatelimit` compression tuning option is now
    available. This is a power-user control to determine how many candidates
    are returned for each block mode encoding trial. This feature is used
	automatically by the search presets; see `-help` for details.
* **Core API:**
  * **Feature:** A new quality preset `ASTCENC_PRE_FASTEST` is available. See
    `-fastest` above for details.
  * **Feature:** A new tuning option `tune_candidate_limit` is available in
    the config structure. See `-candidatelimit` above for details.
  * **Feature:** Image input/output can now use `ASTCENC_TYPE_F32` data types.
* **Stability:**
  * **Improvement:** The SSE and AVX variants now produce identical output when
    run on the same CPU.

### Performance

* Average compression performance is 1.2x - 2x faster than version 2.0,
  depending on search preset and block size.
* Average image quality is similar to 2.0, with only minor differences of
  up to 0.05dB in either direction.

## 2.0

**Status:** Released, August 2020

The 2.0 release is first release in the 2.x series. It includes a number of
major changes over the earlier 1.7 series, and is not command-line compatible.

### Features:

* The core codec can be built as a library, exposed via a new codec API.
* The core codec supports accelerated SIMD paths for SSE2, SSE4.2, and AVX2.
* The command line syntax has a clearer mapping to Khronos feature profiles.

### Performance:

* Average compression performance is between 2 and 3x faster than version 1.7.
* Average image quality is lower by up to 0.1dB than version 1.7.
