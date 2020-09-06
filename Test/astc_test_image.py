#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# -----------------------------------------------------------------------------
# Copyright 2019-2020 Arm Limited
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
import sys

import testlib.encoder as te
import testlib.testset as tts
import testlib.resultset as trs


# Set the reference comparison against either 1.7 or latest 2.0 build
COMPARE_WITH_1_7 = False

# If we are comparing with 1.7 then allow some relaxed thresholds, note that
# some failures are expected here as some images degrade by more than others
if COMPARE_WITH_1_7:
    RESULT_THRESHOLD_WARN = -0.1
    RESULT_THRESHOLD_FAIL = -0.2
    RESULT_THRESHOLD_3D_FAIL = -0.6
    RESULT_REF_NAME = "reference-1.7"
else:
    RESULT_THRESHOLD_WARN = -0.01
    RESULT_THRESHOLD_FAIL = -0.05
    RESULT_THRESHOLD_3D_FAIL = -0.02
    RESULT_REF_NAME = "reference-2.0-avx2"


TEST_BLOCK_SIZES = ["4x4", "5x5", "6x6", "8x8", "12x12",
                    "3x3x3", "6x6x6"]


def is_3d(blockSize):
    """
    Is the given block size a 3D block type?

    Args:
        blockSize (str): The block size.

    Returns:
        bool: ``True`` if the block string is a 3D block size, ``False`` if 2D.
    """
    return blockSize.count("x") == 2


def count_test_set(testSet, blockSizes):
    """
    Count the number of test executions needed for a test set.

    Args:
        testSet (TestSet): The test set to run.
        blockSizes (list(str)): The block sizes to run.

    Returns:
        int: The number of test executions needed.
    """
    count = 0
    for blkSz in blockSizes:
        for image in testSet.tests:
            # 3D block sizes require 3D images
            if is_3d(blkSz) != image.is3D:
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

    Returns:
        Result: The result code.
    """
    dPSNR = result.psnr - reference.psnr

    if (dPSNR < RESULT_THRESHOLD_FAIL) and (not image.is3D):
        return trs.Result.FAIL

    if (dPSNR < RESULT_THRESHOLD_3D_FAIL) and image.is3D:
        return trs.Result.FAIL

    if dPSNR < RESULT_THRESHOLD_WARN:
        return trs.Result.WARN

    return trs.Result.PASS


def format_solo_result(image, result):
    """
    Format a metrics string for a single (no compare) result.

    Args:
        image (TestImage): The image being tested.
        result (Record): The test result.

    Returns:
        str: The metrics string.
    """
    imSize = image.get_size()
    if imSize:
        mpix = float(imSize[0] * imSize[1]) / 1000000.0
        tCMPS = "%3.3f MP/s" % (mpix / result.cTime)
    else:
        tCMPS = "?"

    name = "%5s %s" % (result.blkSz, result.name)
    tPSNR = "%2.3f dB" % result.psnr
    tTTime = "%.2f s" % result.tTime
    tCTime = "%.2f s" % result.cTime

    return "%-32s | %8s | %8s | %8s | %10s" % \
        (name, tPSNR, tTTime, tCTime, tCMPS)


def format_result(image, reference, result):
    """
    Format a metrics string for a comparison result.

    Args:
        image (TestImage): The image being tested.
        reference (Record): The reference result to compare against.
        result (Record): The test result.

    Returns:
        str: The metrics string.
    """
    imSize = image.get_size()
    if imSize:
        mpix = float(imSize[0] * imSize[1]) / 1000000.0
        tCMPS = "%3.3f MP/s" % (mpix / result.cTime)
    else:
        tCMPS = "?"

    dPSNR = result.psnr - reference.psnr
    sTTime = reference.tTime / result.tTime
    sCTime = reference.cTime / result.cTime

    name = "%5s %s" % (result.blkSz, result.name)
    tPSNR = "%2.3f dB (% 1.3f dB)" % (result.psnr, dPSNR)
    tTTime = "%.2f s (%1.2fx)" % (result.tTime, sTTime)
    tCTime = "%.2f s (%1.2fx)" % (result.cTime, sCTime)
    result = determine_result(image, reference, result)

    return "%-32s | %22s | %14s | %14s | %10s | %s" % \
           (name, tPSNR, tTTime, tCTime, tCMPS, result.name)


def run_test_set(encoder, testRef, testSet, blockSizes, testRuns):
    """
    Execute all tests in the test set.

    Args:
        encoder (EncoderBase): The encoder to use.
        testRef (ResultSet): The test reference results.
        testSet (TestSet): The test set.
        blockSizes (list(str)): The block sizes to execute each test against.
        testRuns (int): The number of test repeats to run for each image test.

    Returns:
        ResultSet: The test results.
    """
    resultSet = trs.ResultSet(testSet.name)

    curCount = 0
    maxCount = count_test_set(testSet, blockSizes)

    title = "Test Set: %s / Encoder: %s" % (testSet.name, encoder.name)
    print(title)
    print("=" * len(title))

    for blkSz in blockSizes:
        for image in testSet.tests:
            # 3D block sizes require 3D images
            if is_3d(blkSz) != image.is3D:
                continue

            curCount += 1

            dat = (curCount, maxCount, blkSz, image.testFile)
            print("Running %u/%u %s %s ... " % dat, end='', flush=True)
            res = encoder.run_test(image, blkSz, "-thorough", testRuns)
            res = trs.Record(blkSz, image.testFile, res[0], res[1], res[2])
            resultSet.add_record(res)

            if testRef:
                refResult = testRef.get_matching_record(res)
                res.set_status(determine_result(image, refResult, res))

                res.tTimeRel = refResult.tTime / res.tTime
                res.cTimeRel = refResult.cTime / res.cTime
                res.psnrRel = res.psnr - refResult.psnr

                res = format_result(image, refResult, res)
            else:
                res = format_solo_result(image, res)

            print("\r[%3u] %s" % (curCount, res))

    return resultSet


def get_encoder_params(encoderName, imageSet):
    """
    The the encoder and image set parameters for a test run.

    Args:
        encoderName (str): The encoder name.
        imageSet (str): The test image set.

    Returns:
        tuple(EncoderBase, str, str, str): The test parameters for the
        requested encoder and test set. An instance of the encoder wrapper
        class, the output data name, the output result directory, and the
        reference to use.
    """
    if encoderName == "ref-1.7":
        encoder = te.Encoder1x()
        name = "reference-1.7"
        outDir = "Test/Images/%s" % imageSet
        refName = None
    elif encoderName == "ref-2.0":
        # Note this option rebuilds a new reference test set using the
        # user's locally build encoder.
        encoder = te.Encoder2x("avx2")
        name = "reference-2.0-avx2"
        outDir = "Test/Images/%s" % imageSet
        refName = None
    elif encoderName == "ref-prototype":
        encoder = te.EncoderProto()
        name = "reference-prototype"
        outDir = "Test/Images/%s" % imageSet
        refName = None
    else:
        encoder = te.Encoder2x(encoderName)
        name = "develop-%s" % encoderName
        outDir = "TestOutput/%s" % imageSet
        refName = RESULT_REF_NAME

    return (encoder, name, outDir, refName)


def parse_command_line():
    """
    Parse the command line.

    Returns:
        Namespace: The parsed command line container.
    """
    parser = argparse.ArgumentParser()

    refcoders = ["ref-1.7", "ref-2.0", "ref-prototype"]
    testcoders = ["sse2", "sse4.2", "avx2"]
    coders = refcoders + testcoders + ["all"]
    parser.add_argument("--encoder", dest="encoders", default="avx2",
                        choices=coders, help="test encoder variant")

    astcProfile = ["ldr", "ldrs", "hdr", "all"]
    parser.add_argument("--color-profile", dest="profiles", default="all",
                        choices=astcProfile, help="test color profile")

    imgFormat = ["l", "xy", "rgb", "rgba", "all"]
    parser.add_argument("--color-format", dest="formats", default="all",
                        choices=imgFormat, help="test color format")

    choices = list(TEST_BLOCK_SIZES) + ["all"]
    parser.add_argument("--block-size", dest="blockSizes",
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

    parser.add_argument("--repeats", dest="testRepeats", default=1,
                        type=int, help="test iteration count")

    args = parser.parse_args()

    # Turn things into canonical format lists
    args.encoders = testcoders if args.encoders == "all" \
        else [args.encoders]

    if not args.blockSizes or ("all" in args.blockSizes):
        args.blockSizes = TEST_BLOCK_SIZES

    args.testSets = testSets[:-1] if args.testSets == "all" \
        else [args.testSets]

    args.profiles = astcProfile[:-1] if args.profiles == "all" \
        else [args.profiles]

    args.formats = imgFormat[:-1] if args.formats == "all" \
        else [args.formats]

    return args


def main():
    """
    The main function.

    Returns:
        int: The process return code.
    """
    # Parse command lines
    args = parse_command_line()

    testSetCount = 0
    worstResult = trs.Result.NOTRUN

    for imageSet in args.testSets:
        for encoderName in args.encoders:
            (encoder, name, outDir, refName) = \
                get_encoder_params(encoderName, imageSet)

            testDir = "Test/Images/%s" % imageSet
            testRes = "%s/astc_%s_results.csv" % (outDir, name)

            testRef = None
            if refName:
                testRefPath = "%s/astc_%s_results.csv" % (testDir, refName)
                testRef = trs.ResultSet(imageSet)
                testRef.load_from_file(testRefPath)

            testSetCount += 1
            testSet = tts.TestSet(imageSet, testDir,
                                  args.profiles, args.formats, args.testImage)

            resultSet = run_test_set(encoder, testRef, testSet,
                                     args.blockSizes, args.testRepeats)

            resultSet.save_to_file(testRes)

            if refName:
                summary = resultSet.get_results_summary()
                worstResult = max(summary.get_worst_result(), worstResult)
                print(summary)

    if (testSetCount > 1) and (worstResult != trs.Result.NOTRUN):
        print("OVERALL STATUS: %s" % worstResult.name)

    if worstResult == trs.Result.FAIL:
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
