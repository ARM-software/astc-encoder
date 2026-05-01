#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# -----------------------------------------------------------------------------
# Copyright 2019-2026 Arm Limited
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License. You may obtain a copy
# of the License at:
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations
# under the License.
# -----------------------------------------------------------------------------
"""
The image test runner is used for image quality and performance testing.

It is designed to process directories of arbitrary test images, using the
directory structure and path naming conventions to self-describe how each image
is to be compressed. Some built-in test sets are provided in the ./Test/Images
directory, and others can be downloaded by running the astc_test_image_dl
script.

Attributes:
    RESULT_THRESHOLD_WARN: The result threshold (dB) for getting a WARN.
    RESULT_THRESHOLD_FAIL: The result threshold (dB) for getting a FAIL.
    TEST_BLOCK_SIZES: The block sizes we can test. This is a subset of the
        block sizes supported by ASTC, simply to keep test run times
        manageable.
"""

import argparse
import os
from pathlib import Path
import sys

import testlib.encoder as te
import testlib.testset as tts
import testlib.resultset as trs

# Require bit exact with reference scores
RESULT_THRESHOLD_WARN = -0.00
RESULT_THRESHOLD_FAIL = -0.00
RESULT_THRESHOLD_3D_FAIL = -0.00


TEST_BLOCK_SIZES = [
    "4x4", "5x5", "6x6", "8x8", "12x12", "3x3x3", "6x6x6"
]

TEST_QUALITIES = [
    "fastest", "fast", "medium", "thorough", "verythorough", "exhaustive"
]


def is_3d(block_size):
    """
    Is the given block size a 3D block type?

    Args:
        block_size (str): The block size.

    Return:
        bool: ``True`` if the block string is a 3D block size, ``False`` if 2D.
    """
    return block_size.count("x") == 2


def count_test_set(test_set, block_sizes):
    """
    Count the number of test executions needed for a test set.

    Args:
        test_set (TestSet): The test set to run.
        block_sizes (list(str)): The block sizes to run.

    Return:
        int: The number of test executions needed.
    """
    count = 0
    for block_size in block_sizes:
        for image in test_set.tests:
            # 3D block sizes require 3D images
            if is_3d(block_size) != image.is_3d:
                continue

            count += 1

    return count


def determine_result(image, reference, result):
    """
    Determine a test result against a reference and thresholds.

    Args:
        image (TestImage): The image being compressed.
        reference (Record): The reference result to compare against.
        result (Record): The test result.

    Return:
        ResultStatus: The result code.
    """
    dPSNR = result.psnr - reference.psnr

    if (dPSNR < RESULT_THRESHOLD_FAIL) and (not image.is_3d):
        return trs.ResultStatus.FAIL

    if (dPSNR < RESULT_THRESHOLD_3D_FAIL) and image.is_3d:
        return trs.ResultStatus.FAIL

    if dPSNR < RESULT_THRESHOLD_WARN:
        return trs.ResultStatus.WARN

    return trs.ResultStatus.PASS


def format_solo_result(image, result):
    """
    Format a metrics string for a single (no compare) result.

    Args:
        image (TestImage): The image being tested.
        result (Record): The test result.

    Return:
        str: The metrics string.
    """
    name = "%5s %s" % (result.block_size, result.name)
    tPSNR = "%2.3f dB" % result.psnr
    tTTime = "%.3f s" % result.total_time
    tCTime = "%.3f s" % result.coding_time
    tCMTS = "%.3f MT/s" % result.coding_rate

    return "%-32s | %8s | %9s | %9s | %11s" % \
        (name, tPSNR, tTTime, tCTime, tCMTS)


def format_result(image, reference, result):
    """
    Format a metrics string for a comparison result.

    Args:
        image (TestImage): The image being tested.
        reference (Record): The reference result to compare against.
        result (Record): The test result.

    Return:
        str: The metrics string.
    """
    dPSNR = result.psnr - reference.psnr

    try:
        sTTime = reference.total_time / result.total_time
    except ZeroDivisionError:
        sTTime = float('NaN')

    try:
        sCTime = reference.coding_time / result.coding_time
    except ZeroDivisionError:
        sCTime = float('NaN')

    name = "%5s %s" % (result.block_size, result.name)
    tPSNR = "%2.3f dB (% 1.3f dB)" % (result.psnr, dPSNR)
    tTTime = "%.3f s (%1.2fx)" % (result.total_time, sTTime)
    tCTime = "%.3f s (%1.2fx)" % (result.coding_time, sCTime)
    tCMTS = "%.3f MT/s" % (result.coding_rate)
    result = determine_result(image, reference, result)

    return "%-32s | %22s | %15s | %15s | %11s | %s" % \
           (name, tPSNR, tTTime, tCTime, tCMTS, result.name)


def run_test_set(encoder, testRef, test_set, quality, block_sizes, testRuns,
                 keep_output, threads):
    """
    Execute all tests in the test set.

    Args:
        encoder (EncoderBase): The encoder to use.
        testRef (ResultSet): The test reference results.
        test_set (TestSet): The test set.
        quality (str): The quality level to execute the test against.
        block_sizes (list(str)): The block sizes to execute each test against.
        testRuns (int): The number of test repeats to run for each image test.
        keep_output (bool): Should the test preserve output images? This is
            only a hint and discarding output may be ignored if the encoder
            version used can't do it natively.
        threads (int or None): The thread count to use.

    Return:
        ResultSet: The test results.
    """
    resultSet = trs.ResultSet(test_set.name)

    curCount = 0
    maxCount = count_test_set(test_set, block_sizes)

    dat = (test_set.name, encoder.name, quality)
    title = "Test Set: %s / Encoder: %s -%s" % dat
    print(title)
    print("=" * len(title))

    for block_size in block_sizes:
        for image in test_set.tests:
            # 3D block sizes require 3D images
            if is_3d(block_size) != image.is_3d:
                continue

            curCount += 1

            dat = (curCount, maxCount, block_size, image.file_name)
            print("Running %u/%u %s %s ... " % dat, end='', flush=True)
            res = encoder.run_test(
                image, block_size, "-%s" % quality, testRuns,
                keep_output, threads)

            res = trs.Record(block_size, image.file_name, *res)
            resultSet.add_record(res)

            if testRef:
                refResult = testRef.get_matching_record(res)
                res.set_status(determine_result(image, refResult, res))

                try:
                    res.total_time_rel = refResult.total_time / res.total_time
                except ZeroDivisionError:
                    res.total_time_rel = float('NaN')

                try:
                    res.coding_time_rel = \
                        refResult.coding_time / res.coding_time
                except ZeroDivisionError:
                    res.coding_time_rel = float('NaN')

                res.psnr_rel = res.psnr - refResult.psnr

                res = format_result(image, refResult, res)
            else:
                res = format_solo_result(image, res)

            print("\r[%3u] %s" % (curCount, res))

    return resultSet


def get_encoder_params(encoderName, referenceName, imageSet):
    """
    The the encoder and image set parameters for a test run.

    Args:
        encoderName (str): The encoder name.
        referenceName (str): The reference encoder name.
        imageSet (str): The test image set.

    Return:
        tuple(EncoderBase, str, str, str): The test parameters for the
        requested encoder and test set. An instance of the encoder wrapper
        class, the output data name, the output result directory, and the
        reference to use.
    """
    if encoderName.startswith("ref"):
        _, version, simd = encoderName.split("-")

        # 2.x, 3.x, and 4.x variants
        compatible2xPrefixes = ["2.", "3.", "4.", "5."]
        if any(True for x in compatible2xPrefixes if version.startswith(x)):
            encoder = te.Encoder2xRel(version, simd)
            name = f"reference-{version}-{simd}"
            outDir = "Test/Images/%s" % imageSet
            refName = None
            return (encoder, name, outDir, refName)

        # Latest main
        if version == "main":
            encoder = te.Encoder2x(simd)
            name = f"reference-{version}-{simd}"
            outDir = "Test/Images/%s" % imageSet
            refName = None
            return (encoder, name, outDir, refName)

        assert False, f"Encoder {encoderName} not recognized"

    encoder = te.Encoder2x(encoderName)
    name = "develop-%s" % encoderName
    outDir = "TestOutput/%s" % imageSet
    refName = referenceName.replace("ref", "reference")
    return (encoder, name, outDir, refName)


def parse_command_line():
    """
    Parse the command line.

    Return:
        Namespace: The parsed command line container.
    """
    parser = argparse.ArgumentParser()

    # All reference encoders
    refcoders = [
        "ref-5.0-neon",
        "ref-5.0-sse2", "ref-5.0-sse4.1", "ref-5.0-avx2",

        "ref-main-neon", "ref-main-sve_256", "ref-main-sve_128",
        "ref-main-sse2", "ref-main-sse4.1", "ref-main-avx2"
    ]

    # All test encoders
    testcoders = [
        "none", "native", "universal",
        "neon", "sve_256", "sve_128",
        "sse2", "sse4.1", "avx2"
    ]

    testcodersAArch64 = [
        "neon", "sve_256", "sve_128"
    ]

    testcodersX86 = [
        "sse2", "sse4.1", "avx2"
    ]

    coders = refcoders + testcoders + ["all-aarch64", "all-x86"]

    parser.add_argument("--encoder", dest="encoders",
                        default="avx2",
                        choices=coders, help="test encoder variant")

    parser.add_argument("--reference", default="ref-main-avx2",
                        choices=refcoders, help="reference encoder variant")

    astcProfile = ["ldr", "ldrs", "hdr", "all"]
    parser.add_argument("--color-profile", dest="profiles", default="all",
                        choices=astcProfile, help="test color profile")

    imgFormat = ["l", "xy", "rgb", "rgba", "all"]
    parser.add_argument("--color-format", dest="formats", default="all",
                        choices=imgFormat, help="test color format")

    choices = list(TEST_BLOCK_SIZES) + ["all"]
    parser.add_argument("--block-size", dest="block_sizes",
                        action="append", choices=choices,
                        help="test block size")

    testDir = os.path.dirname(__file__)
    testDir = os.path.join(testDir, "Images")
    testSets = []
    for path in os.listdir(testDir):
        fqPath = os.path.join(testDir, path)
        if os.path.isdir(fqPath):
            testSets.append(path)
    testSets.append("all")

    parser.add_argument("--test-set", dest="testSets", default="Small",
                        choices=testSets, help="test image test set")

    parser.add_argument("--test-image", dest="testImage", default=None,
                        help="select a specific test image from the test set")

    choices = list(TEST_QUALITIES) + ["all", "all+"]
    parser.add_argument("--test-quality", dest="testQual", default="thorough",
                        choices=choices, help="select a specific test quality")

    parser.add_argument("--repeats", dest="testRepeats", default=1,
                        type=int, help="test iteration count")

    parser.add_argument("--keep-output", dest="keep_output", default=False,
                        action="store_true", help="keep image output")

    parser.add_argument("-j", dest="threads", default=None,
                        type=int, help="thread count")

    args = parser.parse_args()

    # Turn things into canonical format lists
    if args.encoders == "all-aarch64":
        args.encoders = testcodersAArch64
    elif args.encoders == "all-x86":
        args.encoders = testcodersX86
    else:
        args.encoders = [args.encoders]

    if args.testQual == "all+":
        args.testQual = TEST_QUALITIES
    elif args.testQual == "all":
        args.testQual = TEST_QUALITIES
        args.testQual.remove("verythorough")
        args.testQual.remove("exhaustive")
    else:
        args.testQual = [args.testQual]

    if not args.block_sizes or ("all" in args.block_sizes):
        args.block_sizes = TEST_BLOCK_SIZES

    args.testSets = testSets[:-1] if args.testSets == "all" \
        else [args.testSets]

    args.profiles = astcProfile[:-1] if args.profiles == "all" \
        else [args.profiles]

    args.formats = imgFormat[:-1] if args.formats == "all" \
        else [args.formats]

    return args


def main() -> int:
    """
    The main function.

    Return:
        int: The process return code.
    """
    # Parse command lines
    args = parse_command_line()

    testSetCount = 0
    worstResult = trs.ResultStatus.NOT_RUN

    for quality in args.testQual:
        for imageSet in args.testSets:
            for encoderName in args.encoders:
                (encoder, name, outDir, refName) = \
                    get_encoder_params(encoderName, args.reference, imageSet)

                testDir = "Test/Images/%s" % imageSet
                testRes = "%s/astc_%s_%s_results.csv" % (outDir, name, quality)

                testRef = None
                if refName:
                    dat = (testDir, refName, quality)
                    testRefPath = "%s/astc_%s_%s_results.csv" % dat
                    testRef = trs.ResultSet(imageSet)
                    testRef.load_from_file(Path(testRefPath))

                testSetCount += 1
                test_set = tts.TestSet(
                    imageSet, Path(testDir),
                    args.profiles, args.formats, args.testImage)

                resultSet = run_test_set(
                    encoder, testRef, test_set, quality,
                    args.block_sizes, args.testRepeats,
                    args.keep_output, args.threads)

                resultSet.save_to_file(Path(testRes))

                if refName:
                    summary = resultSet.get_results_summary()
                    worstResult = max(summary.get_worst_result(), worstResult)
                    print(summary)

        if (testSetCount > 1) and (worstResult != trs.ResultStatus.NOT_RUN):
            print("OVERALL STATUS: %s" % worstResult.name)

    if worstResult == trs.ResultStatus.FAIL:
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
