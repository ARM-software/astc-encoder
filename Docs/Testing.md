# Testing astcenc

The repository contains a suite of tests which can be used to validate source
code changes to the compressor. The test suite is split into four main test
types:

* Python tests that drive the command line functional interface tests.
* Python tests that drive command line performance and image quality tests.
* C++ unit tests for some internal codec components.
* C++ threading monkey tests for concurrent use of the library API.

Most library API coverage comes from the command line functional interface
tests. A separate C++ threading monkey test exercises concurrent library use,
including context inheritance and compression cancellation.

In addition to the tests, we have automated Python script lint, typing, and
style checkers.

## Prerequisites

Running the Python tests requires Python 3.12 to be installed on the host
machine.

# Running Python-based CLI functional tests

Compile the build of `astcenc` you wish to test, and install it into the
project `./bin` directory. This can be achieved by configuring the CMake build
using the install prefix `-DCMAKE_INSTALL_PREFIX=../` and then running a build
with the `install` build target.

Run the functional tests against the appropriate encoder variant using:

```sh
python3 ./Test/astc_test_functional.py -v --encoder <encoder_variant>
```

To test a binary directly without installing it into `./bin`, pass its path:

```sh
python3 ./Test/astc_test_functional.py -v --binary <path_to_astcenc>
```

`--binary` overrides the path selected by `--encoder`. Without `--binary`,
the default remains `./bin/astcenc-avx2` (or the selected encoder variant).

# Running Python-based CLI performance and quality tests

Compile the build of `astcenc` you wish to test and profile, and install it
into the project `./bin` directory. This can be achieved by configuring the
CMake build using the install prefix `-DCMAKE_INSTALL_PREFIX=../` and then
running a build with the `install` build target.

To run the image test suite run the following command from the root directory
of the repository:

```sh
python3 ./Test/astc_test_image.py
```

This will run through a series of image compression tests, comparing the image
PSNR against a set of reference results from the last stable baseline. The test
will fail if any reduction in PSNR above a set threshold is detected.

The `--test-set` parameter can be passed to select which test images you use.
The `--color-format`, `--block-size`, and `--test-quality` options allow you
to subset which compressor options are tested. The `-j` option lets you control
how many threads the compressor is allowed to use.

For a smaller smoke test, pass `--pairwise`. This groups images by color profile
and format, covers every valid pair of image group, block size, and selected
quality, and tests each image with a compatible block size at least once. For
the Small set, all four default presets take 46 runs per encoder, compared with
596 for the full matrix.

By default output images are discarded during testing. To store images for
debug purposes, the decompressed output images can be kept by passing
`--keep-output`. The resulting images are stored in the `TestOutput`
subdirectory, using the same test set directory structure as the `Test/Images`
folder.

## Benchmarking

All runs of this test suite will report performance data. To improve result
stability, the `--repeats` option makes the test suite process each image
multiple times and return the best measured performance.

Performance regressions will not automatically cause test failure, and must be
manually reviewed.

**Note:** The upstream reference data contains performance results measured on
an Intel Xeon w7-3465X running at 3.9 GHz, running each test 5 times.

## Updating reference data

The reference quality and performance scores are stored in CSVs committed to
the repository. This data is created by running the tests using the last stable
release on a standard test machine we use for performance testing builds.

It can be useful for developers to rebuild the reference results for their
local machine, in particular for measuring performance improvements.

To build new reference CSVs for an official release, checkout and build a
release build of the baseline release you want to compare against and place it
in to the `./Binaries/<version>/` directory, where `<version>` must be the
`major.minor` identifier used by the test script, e.g. `5.7`. For AVX2, place
the binary at `./Binaries/5.7/astcenc-avx2` (`astcenc-avx2.exe` on Windows).

Once this is done, run the command (e.g.):

```sh
python3 ./Test/astc_test_image.py --encoder ref-5.7-avx2 --reference rebuild ...
```

... for whichever test set you want to build a reference for. This can take
some time to complete, especially if your command line enables many test sets,
high quality levels, or uses many repeats to improve result stability.

Add `--quality-update` to update existing CSV records only when their PSNR
changes. Records with unchanged PSNR are not updated to make it easier to
review impact of changes in PRs.

**WARNING:** This may take some hours to complete, and it is best done when the
test suite gets exclusive use of the machine to avoid other processing slowing
down the compression and disturbing the performance data. It is recommended to
shut down or disable any background applications that are running.

# Running C++ unit tests

We support a small (but growing) number of C++ unit tests, which are written
using the googletest framework and integrated in the CMake CTest test runner.

To build unit tests pull the `googletest` git submodule and add
`-DASTCENC_UNITTEST=ON` to the CMake command line when configuring.

To run unit tests use the CMake `ctest` utility from your build directory after
you have built the tests.

```shell
cd build
ctest --verbose
```

# Running the threading monkey test

The monkey test runs independent workers that repeatedly initialize
configurations, create, reuse, and free contexts, compress images, and
decompress images compressed by other workers. It uses standalone contexts
and parent-child chains, including creating children from active parents and
freeing intermediate parents while their descendants remain alive. Root
contexts remain alive until all of their descendants have been freed.

Enabling `ASTCENC_UNITTEST` also builds the monkey test as a separate executable
for each selected ISA in full codec builds without diagnostic tracing. It is
not registered with CTest. Build it for the ISA you want to test:

```sh
cmake -S . -B build-monkey -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DASTCENC_ISA_AVX2=ON -DASTCENC_CLI=OFF -DASTCENC_UNITTEST=ON
cmake --build build-monkey -j
./build-monkey/Source/test-monkey-avx2 --duration 60 --workers 4 --seed 1
```

`--duration` sets the randomized loop duration in seconds, from 1 to 86400,
defaulting to 30.
Reference generation and finishing operations already in progress add some
time. A watchdog terminates the process with a failure if the test exceeds
the requested duration plus 30 seconds, including reference generation.
`--workers` sets the number of independent workers, from 2 to 64, defaulting
to four. `--max-threads` limits the configured thread count per context, from
1 to 16, defaulting to four. Invalid numeric values report the option,
rejected value, and allowed range. Operations randomly use one or more of
each context's thread indices, with staggered starts. Each context processes
one image at a time, and all participating threads join before it is reset
or freed.

The test reuses four generated images across a fixed settings matrix
covering LDR, sRGB, HDR, normal maps, alpha weighting and preprocessing, 2D
and 3D block sizes, several quality presets, and decompress-only contexts.
It first records serial configurations and single-threaded compressed and
decompressed reference bytes. Newly initialized configurations must match
their serial references field-for-field. Every completed compression and
decompression must match its reference bytes exactly. Cancelled output is
discarded; the same context is then reset and used for a complete
compression, which must match the reference. Cancellation is requested
externally and from progress callbacks. Only compression cancellation is
supported by the API.

Aggregate context and operation counts are logged every five seconds and at
shutdown. These include validated configurations, compressions and
decompressions, transfers between workers, and cancellation requests. Final
output also shows compression and decompression counts for each
image/settings pair. Settings 0 through 7 are:

| Setting | Profile | Block size | Quality | Additional settings |
| --- | --- | --- | --- | --- |
| 0 | LDR | 4x4 | Fastest | None |
| 1 | LDR | 6x6 | Fast | Alpha weighting, alpha-scale radius 4 |
| 2 | sRGB | 8x8 | Medium | UNORM8 decode rounding |
| 3 | LDR | 4x4 | Medium | Normal map |
| 4 | HDR | 6x6 | Fast | None |
| 5 | HDR RGB, LDR alpha | 4x4 | Fastest | None |
| 6 | LDR | 4x4x4 | Fastest | None |
| 7 | LDR | 6x6 | Fast | Decompress-only context |

Images 0 through 3 have dimensions 33x29, 127x65, 257x513, and 25x19x5.
Each 2D setting uses all three 2D images; setting 6 uses the 3D image.
Setting 7 has no compression counts because its contexts only decompress.

Failures report the seed, worker, iteration, and reference index. Watchdog
timeouts also report each worker's current operation. Reusing a seed gives
each worker the same initial random state, but scheduling and cross-worker
exchanges affect later choices, so a seed does not reproduce an exact run.
The executable returns nonzero on failure and is kept separate from the CTest
unit suite.

ASAN, TSAN, and UBSAN can be enabled using the options below. For a ThreadSanitizer
build with Clang or GCC on Linux, configure a separate build directory:

```sh
cmake -S . -B build-monkey-tsan -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DASTCENC_ISA_AVX2=ON -DASTCENC_CLI=OFF -DASTCENC_UNITTEST=ON \
    -DASTCENC_TSAN=ON
cmake --build build-monkey-tsan -j
TSAN_OPTIONS=halt_on_error=1 \
    ./build-monkey-tsan/Source/test-monkey-avx2 --duration 60 --seed 1
```

# Using ASAN, TSAN, and UBSAN

Running builds compiled with ASAN (address sanitizer), TSAN (thread sanitizer),
and UBSAN (undefined behavior sanitizer) can detect memory errors, data races,
and undefined behavior. Support is integrated into the build system for GCC
and Clang on Linux targets.

Build a release build with debug information using
`-DCMAKE_BUILD_TYPE=RelWithDebInfo` and specify one of
`-DASTCENC_ASAN=ON`, `-DASTCENC_TSAN=ON`, or `-DASTCENC_UBSAN=ON`.
ASAN and TSAN require separate builds and cannot be enabled together.

Run test commands as normal with the generated binaries.

# Using Valgrind memcheck

Running builds using the Valgrind memcheck tool is another way to validate that
we have not introduced memory errors.

Build a release build with debug information using
`-DCMAKE_BUILD_TYPE=RelWithDebInfo`.

Run test commands through Valgrind:

    valgrind --tool=memcheck --track-origins=yes <command>

# OSS-Fuzz

ASTC-Encoder has been integrated into the Google OSS-Fuzz program, which
performs API fuzz testing on a Google-hosted CI infrastructure.

Our OSS-Fuzz test harnesses can be found in the
[/Source/Fuzzers](../Source/Fuzzers/) directory, and we have some short
[documentation](TestingOSSFuzz.md) about how to run the fuzzers locally.

The OSS-Fuzz Project is hosted on GitHub at [google/oss-fuzz][1].

The OSS-Fuzz CI is visible to maintainers at [oss-fuzz.com][2].

[1]: https://github.com/google/oss-fuzz
[2]: https://oss-fuzz.com

- - -

_Copyright © 2019-2026, Arm Limited and contributors._
