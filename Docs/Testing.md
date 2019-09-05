# Testing ASTC Encoder

The repository contains a small suite of application tests, which can be used
to sanity check source code changes to the compressor. It must be noted that
this test suite is relatively limited in scope and does not cover every feature
or bitrate of the standard.

# Required Software

Running the tests requires Python 3.7 to be installed on the host machine, with
the following packages installed:

    python3 -m pip install junit_xml
    python3 -m pip install pillow

# Running the Test Suite

To run the full test suite first build the 64-bit release configuration of
`astcenc` for the host OS, and then run the following command from the root of
directory of the repository:

    python3 ./Test/runner.py

This will run though a series of image compression tests, comparing the
resulting PSNR against the reference results from the last stable release. The
test will fail if any reduction in PSNR is detected.

All decompressed test output images are stored in the `TestOutput` directory,
as is a summary result report in JUnit XML format which also includes test
compression time. Note that while compression speed is reported by the test
suite, a performance regression will not cause a test to fail.

## Smoke tests

To quickly sanity check changes you can run a smaller test list, but that
stills runs a few tests from each image category, using the following commands:

    python3 ./Test/runner.py --test-level smoke

You can further filter tests by selecting runs for specific profiles, data
formats, or block sizes. See the following options in the `--help` for
more details:

* `--dynamic-range` : select a single set of either LDR or HDR tests.
* `--format` : select a single input data format.
* `--block-size` : select a single block size.

## Performance tests

To provide less noisy performance results the test suite supports running each
compression pass multiple times and returning the average performance. To
enable this mode use the following two options:

* `--warmup <N>` : Run N warmup compression passes which are not timed.
* `--repeats <M>` : Run M test compression passes which are timed.

**Note:**  The reference CSV contains performance results measured on an Intel
Core i5 9600K running at 4.3GHz, running each test 10 times after a single
warmup pass.

# Updating Reference Scores

The pass and fail conditions are stored in a reference result CSV, which is
based on the image quality and performance of the latest stable tag. The
runner can be used to regenerate the CSV file using the 64-bit release build
binary in the [Binary directory](/Binary/).

* `--rebuild-ref-csv` : regenerate the whole reference, rerunning all tests.
* `--update-ref-csv` : patch the reference to add new test images, but keep
  any existing results.

# Known limitations

The current test suite is viewed as a set of bare-essentials, but has a number
of significant omissions:

* Only square block sizes from 4x4 (8bpp) up to 8x8 (2pp) are tested.
* Only `-thorough` compression speed tested.
* Few optional compressor options are tested other than the use of
  `-normal_psnr` for the `xy` data test set.
* No LDR profile coverage of luminance-only input textures.
* Limited HDR profile coverage of HDR input textures.
* No Full profile coverage of 3D input textures.

It is intended that pair-wise test coverage should be used in future to allow
us to have have wider coverage without excessive test runtime.
