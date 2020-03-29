#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# -----------------------------------------------------------------------------
# Copyright 2020 Arm Limited
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
The image test runner is designed to process directories of arbitrary test
images, using the directory structure and naming conventions to self-describe
how each image is to be compressed.
"""
import filecmp
import os
import signal
import subprocess as sp
import sys
import tempfile
import unittest

import numpy
from PIL import Image

# Enable these to always write out, irrespective of test result
ASTCENC_CLI_ALWAYS = False
ASTCENC_LOG_ALWAYS = False

# Enable these to write out on failure for positive tests
ASTCENC_CLI_ON_ERROR = True
ASTCENC_LOG_ON_ERROR = True

# Enable these to write out on failure for negative tests
ASTCENC_CLI_ON_ERROR_NEG = True
ASTCENC_LOG_ON_ERROR_NEG = True


class CLITestBase(unittest.TestCase):
    """
    Command line interface base class.

    These tests are designed to test the command line is handled correctly.
    They are not detailed tests of the codec itself; only basic sanity checks
    that some type of processing occurred are used.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # TODO: Replace this with a cross-platform solution
        self.binary = os.path.join(".", "Source", "astcenc-avx2")

    def setUp(self):
        """
        Set up a test case.
        """
        self.tempDir = tempfile.TemporaryDirectory()

    def tearDown(self):
        """
        Tear down a test case.
        """
        self.tempDir.cleanup()
        self.tempDir = None

    @staticmethod
    def get_ref_image_path(profile, mode, image):
        """
        Get the path of a reference image on disk.

        Args:
            profile (str): The color profile.
            mode (str): The type of image to load.
            image (str): The image variant to load.

        Returns:
            str: The path to the test image file on disk.
        """
        nameMux = {
            "LDR": {
                "input": "png",
                "comp": "astc"
            },
            "LDRS": {
                "input": "png",
                "comp": "astc"
            },
            "HDR": {
                "input": "exr",
                "comp": "astc"
            }
        }

        assert profile in nameMux.keys()
        assert mode in nameMux["LDR"].keys()

        scriptDir = os.path.dirname(__file__)
        fileName = "%s-%s-1x1.%s" % (profile, image, nameMux[profile][mode])
        return os.path.join(scriptDir, "Data", fileName)

    def get_tmp_image_path(self, profile, mode):
        """
        Get the path of a temporary output image on disk.

        Temporary files are automatically cleaned up when the test tearDown
        occurs.

        Args:
            profile (str): The color profile.
            mode (str): The type of image to load.

        Returns:
            str: The path to the test image file on disk.
        """
        nameMux = {
            "LDR": {
                "comp": ".astc",
                "decomp": ".png",
                "bad": ".foo"
            },
            "LDRS": {
                "comp": ".astc",
                "decomp": ".png",
                "bad": ".foo"
            },
            "HDR": {
                "comp": ".astc",
                "decomp": ".exr",
                "bad": ".foo"
            }
        }

        assert profile in nameMux.keys()
        assert mode in nameMux["LDR"].keys()

        suffix = nameMux[profile][mode]
        tmpFile, tmpPath = tempfile.mkstemp(suffix, dir=self.tempDir.name)
        os.close(tmpFile)
        os.remove(tmpPath)
        return tmpPath


class CLIPTest(CLITestBase):
    """
    Command line interface positive tests.

    These tests are designed to test the command line is handled correctly.
    They are not detailed tests of the codec itself; only basic sanity checks
    that some type of processing occurred are used.
    """

    def compare(self, image1, image2):
        """
        Utility function to compare two images.

        Note that this comparison tests only decoded color values; any file
        metadata is ignored, and encoding methods are not compared.

        Args:
            image1 (str): Path to the first image.
            image2 (str): Path to the second image.

        Returns:
            bool: ``True` if the images are the same, ``False`` otherwise.
        """
        img1 = Image.open(image1)
        img2 = Image.open(image2)

        # Images must have same size
        if img1.size != img2.size:
            print("Size")
            return False

        # Images must have same number of color channels
        if img1.getbands() != img2.getbands():
            # ... unless the only different is alpha
            self.assertEqual(img1.getbands(), ("R", "G", "B"))
            self.assertEqual(img2.getbands(), ("R", "G", "B", "A"))

            # ... and the alpha is always one
            bands = img2.split()
            alphaHist = bands[3].histogram()
            self.assertEqual(sum(alphaHist[:-1]), 0)

            # Generate a version of img2 without alpha
            img2 = Image.merge("RGB", (bands[0], bands[1], bands[2]))

        # Compute sum of absolute differences
        dat1 = numpy.array(img1)
        dat2 = numpy.array(img2)
        sad = numpy.sum(numpy.abs(dat1 - dat2))

        if sad != 0:
            print(img1.load()[0, 0])
            print(img2.load()[0, 0])

        return sad == 0

    def exec(self, command):
        """
        Execute a positive test.

        Test will automatically fail if the subprocess return code is any value
        other than zero.

        Args:
            command (list(str)): The command to execute.

        Returns:
            str: The stdout output of the child process.
        """
        try:
            result = sp.run(command, stdout=sp.PIPE, stderr=sp.PIPE,
                            universal_newlines=True, check=True)
            error = False
        except sp.CalledProcessError as ex:
            result = ex
            error = True

        # Emit debug logging if needed
        if ASTCENC_CLI_ALWAYS or (error and ASTCENC_CLI_ON_ERROR):
            print(" ".join(command))

        if ASTCENC_LOG_ALWAYS or (error and ASTCENC_LOG_ON_ERROR):
            print(result.stdout)

        rcode = result.returncode

        if rcode < 0:
            msg = "Exec died with signal %s" % signal.Signals(-rcode).name
            self.assertGreaterEqual(rcode, 0, msg)

        if rcode > 0:
            msg = "Exec died with application error %u" % rcode
            self.assertEqual(rcode, 0, msg)

        return result.stdout

    def test_ldr_compress(self):
        """
        Test basic LDR compression.
        """
        imIn = self.get_ref_image_path("LDR", "input", "A")
        imOut = self.get_tmp_image_path("LDR", "comp")
        imRef = self.get_ref_image_path("LDR", "comp", "A")

        command = [self.binary, "-cl", imIn, imOut, "6x6", "-exhaustive"]
        self.exec(command)
        self.assertTrue(filecmp.cmp(imRef, imOut, False))

    def test_srgb_compress(self):
        """
        Test basic LDR sRGB compression.
        """
        imIn = self.get_ref_image_path("LDRS", "input", "A")
        imOut = self.get_tmp_image_path("LDRS", "comp")
        imRef = self.get_ref_image_path("LDRS", "comp", "A")

        command = [self.binary, "-cs", imIn, imOut, "6x6", "-exhaustive"]
        self.exec(command)
        self.assertTrue(filecmp.cmp(imRef, imOut, False))

    @unittest.skip("Not implemented")
    def test_hdr_compress(self):
        """
        Test basic HDR compression.
        """

    def test_ldr_decompress(self):
        """
        Test basic LDR decompression.
        """
        imIn = self.get_ref_image_path("LDR", "comp", "A")
        imOut = self.get_tmp_image_path("LDR", "decomp")
        imRef = self.get_ref_image_path("LDR", "input", "A")

        command = [self.binary, "-dl", imIn, imOut]
        self.exec(command)
        self.assertTrue(self.compare(imRef, imOut))

    def test_srgb_decompress(self):
        """
        Test basic LDR sRGB decompression.
        """
        imIn = self.get_ref_image_path("LDRS", "comp", "A")
        imOut = self.get_tmp_image_path("LDRS", "decomp")
        imRef = self.get_ref_image_path("LDRS", "input", "A")

        command = [self.binary, "-ds", imIn, imOut]
        self.exec(command)
        self.assertTrue(self.compare(imRef, imOut))

    @unittest.skip("Not implemented")
    def test_hdr_decompress(self):
        """
        Test basic HDR decompression.
        """

    def test_ldr_roundtrip(self):
        """
        Test basic LDR round-trip
        """
        imIn = self.get_ref_image_path("LDR", "input", "A")
        imOut = self.get_tmp_image_path("LDR", "decomp")

        command = [self.binary, "-tl", imIn, imOut, "6x6", "-exhaustive"]
        self.exec(command)
        self.assertTrue(self.compare(imIn, imOut))

    def test_srgb_roundtrip(self):
        """
        Test basic LDR sRGB round-trip
        """
        imIn = self.get_ref_image_path("LDRS", "input", "A")
        imOut = self.get_tmp_image_path("LDRS", "decomp")

        command = [self.binary, "-ts", imIn, imOut, "6x6", "-exhaustive"]
        self.exec(command)
        self.assertTrue(self.compare(imIn, imOut))

    @unittest.skip("Not implemented")
    def test_hdr_roundtrip(self):
        """
        Test basic HDR round-trip.
        """


class CLINTest(CLITestBase):
    """
    Command line interface negative tests.

    These tests are designed to test that bad inputs to the command line are
    handled correctly, and that errors are correctly thrown.
    """

    def exec(self, command, expectPass=False):
        """
        Execute a negative test.

        Test will automatically fail if:

        * The subprocess return code is zero, unless ``expectPass==True``.
        * The subprocess correctly returned non-zero, but without any error
          message.
        * The subprocess dies with any kind of signal.

        Args:
            command (list(str)): The command to execute.
            expectPass (bool): ``True`` if this command is actually expected to
                pass, which is used to validate commands before mutating them.
        """
        try:
            result = sp.run(command, stdout=sp.PIPE, stderr=sp.PIPE,
                            universal_newlines=True, check=True)
            error = False
        except sp.CalledProcessError as ex:
            # Pop out of the CPE scope to handle the error, as this reduces
            # test log verbosity on failure by avoiding nested exceptions
            result = ex
            error = True

        # Emit debug logging if needed
        badResult = error == expectPass
        if ASTCENC_CLI_ALWAYS or (badResult and ASTCENC_CLI_ON_ERROR_NEG):
            print(" ".join(command))

        if ASTCENC_LOG_ALWAYS or (badResult and ASTCENC_LOG_ON_ERROR_NEG):
            print(result.stdout)

        rcode = result.returncode

        # If we expected a pass, then rcode == 0
        if expectPass:
            self.assertEqual(rcode, 0, "Exec did not pass as expected")
            self.assertNotIn("ERROR", result.stdout)
            return

        # If we got a negative that's always bad (signal of some kind)
        if rcode < 0:
            msg = "Exec died with signal %s" % signal.Signals(-rcode).name
            self.assertGreaterEqual(rcode, 0, msg)

        # Otherwise just assert that we got an error log, and some positive
        # return code value was returned

        # TODO: Disabled until we fix issue #100
        # self.assertIn("ERROR", result.stdout)

        self.assertGreater(rcode, 0, "Exec did not fail as expected")

    def test_cl_missing_args(self):
        """
        Test -cl with missing arguments.
        """
        # Build a valid command
        command = [
            self.binary, "-cl",
            self.get_ref_image_path("LDR", "input", "A"),
            self.get_tmp_image_path("LDR", "comp"),
            "4x4", "-fast"]

        # Run the command, incrementally omitting arguments
        commandLen = len(command)
        for subLen in range(2, commandLen + 1):
            omit = len(command) - subLen
            with self.subTest(omit=omit):
                testCommand = command[:subLen]

                # For the last run we omit no arguments; make sure this works
                # to ensure that the underlying test is actually valid and
                # we're not failing for a different reason
                expectPass = omit == 0
                self.exec(testCommand, expectPass)

    @unittest.skip("Bug #93")
    def test_cl_missing_input(self):
        """
        Test -cl with a missing input file.
        """
        # Build a valid command with a missing input file
        command = [
            self.binary, "-cl",
            "./Test/Data/missing.png",
            self.get_tmp_image_path("LDR", "comp"),
            "4x4", "-fast"]

        self.exec(command)

    @unittest.skip("Bug #93")
    def test_cl_unknown_input(self):
        """
        Test -cl with an unknown input file extension.
        """
        # Build an otherwise valid command with the test flaw
        command = [
            self.binary, "-cl",
            "./Test/Data/empty.unk",
            self.get_tmp_image_path("LDR", "comp"),
            "4x4", "-fast"]

        self.exec(command)

    @unittest.skip("Bug #93")
    def test_cl_missing_output(self):
        """
        Test -cl with a missing output directory.
        """
        # Build an otherwise valid command with the test flaw
        command = [
            self.binary, "-cl",
            self.get_ref_image_path("LDR", "input", "A"),
            "./DoesNotExist/test.astc",
            "4x4", "-fast"]

        self.exec(command)

    @unittest.skip("Bug #99")
    def test_cl_bad_block_size(self):
        """
        Test -cl with an invalid block size.
        """
        badSizes = [
            "4x5",      # Illegal 2D block size
            "3x3x4",    # Illegal 3D block size
            "4x4x4x4",  # Too many dimensions
            "4x",       # Incomplete 2D block size
            "4x4x",     # Incomplete 3D block size
            "4xe",      # Illegal non-numeric character
            "4x4e"      # Additional non-numeric character
        ]

        # Build an otherwise valid command with the test flaw
        command = [
            self.binary, "-cl",
            self.get_ref_image_path("LDR", "input", "A"),
            self.get_tmp_image_path("LDR", "comp"),
            "4x4", "-fast"]

        # Test that the underlying command is valid
        self.exec(command, True)

        blockIndex = command.index("4x4")
        for badSize in badSizes:
            with self.subTest(blockSize=badSize):
                command[blockIndex] = badSize
                self.exec(command)

    def test_cl_bad_preset(self):
        """
        Test -cl with an invalid encoding preset.
        """
        # Build an otherwise valid command with the test flaw
        command = [
            self.binary, "-cl",
            self.get_ref_image_path("LDR", "input", "A"),
            self.get_tmp_image_path("LDR", "comp"),
            "4x4", "-fastt"]

        self.exec(command)

    def test_cl_bad_argument(self):
        """
        Test -cl with an unknown additional argument.
        """
        # Build an otherwise valid command with the test flaw
        command = [
            self.binary, "-cl",
            self.get_ref_image_path("LDR", "input", "A"),
            self.get_tmp_image_path("LDR", "comp"),
            "4x4", "-fast", "-unknown"]

        self.exec(command)

    def test_tl_missing_args(self):
        """
        Test -tl with missing arguments.
        """
        # Build a valid command
        command = [
            self.binary, "-tl",
            self.get_ref_image_path("LDR", "input", "A"),
            self.get_tmp_image_path("LDR", "decomp"),
            "4x4", "-fast"]

        # Run the command, incrementally omitting arguments
        commandLen = len(command)
        for subLen in range(2, commandLen + 1):
            omit = len(command) - subLen
            with self.subTest(omit=omit):
                testCommand = command[:subLen]

                # For the last run we omit no arguments; make sure this works
                # to ensure that the underlying test is actually valid and
                # we're not failing for a different reason
                expectPass = omit == 0
                self.exec(testCommand, expectPass)

    @unittest.skip("Bug #93")
    def test_tl_missing_input(self):
        """
        Test -tl with a missing input file.
        """
        # Build a valid command with a missing input file
        command = [
            self.binary, "-tl",
            "./Test/Data/missing.png",
            self.get_tmp_image_path("LDR", "decomp"),
            "4x4", "-fast"]

        self.exec(command)

    @unittest.skip("Bug #93")
    def test_tl_unknown_input(self):
        """
        Test -tl with an unknown input file extension.
        """
        # Build an otherwise valid command with the test flaw
        command = [
            self.binary, "-tl",
            "./Test/Data/empty.unk",
            self.get_tmp_image_path("LDR", "decomp"),
            "4x4", "-fast"]

        self.exec(command)

    @unittest.skip("Bug #93")
    def test_tl_missing_output(self):
        """
        Test -tl with a missing output directory.
        """
        # Build an otherwise valid command with the test flaw
        command = [
            self.binary, "-tl",
            self.get_ref_image_path("LDR", "input", "A"),
            "./DoesNotExist/test.png",
            "4x4", "-fast"]

        self.exec(command)

    @unittest.skip("Bug #99")
    def test_tl_bad_block_size(self):
        """
        Test -tl with an invalid block size.
        """
        badSizes = [
            "4x5",      # Illegal 2D block size
            "3x3x4",    # Illegal 3D block size
            "4x4x4x4",  # Too many dimensions
            "4x",       # Incomplete 2D block size
            "4x4x",     # Incomplete 3D block size
            "4xe",      # Illegal non-numeric character
            "4x4e"      # Additional non-numeric character
        ]

        # Build an otherwise valid command with the test flaw
        command = [
            self.binary, "-tl",
            self.get_ref_image_path("LDR", "input", "A"),
            self.get_tmp_image_path("LDR", "decomp"),
            "4x4", "-fast"]

        # Test that the underlying command is valid
        self.exec(command, True)

        blockIndex = command.index("4x4")
        for badSize in badSizes:
            with self.subTest(blockSize=badSize):
                command[blockIndex] = badSize
                self.exec(command)

    def test_tl_bad_preset(self):
        """
        Test -tl with an invalid encoding preset.
        """
        # Build an otherwise valid command with the test flaw
        command = [
            self.binary, "-tl",
            self.get_ref_image_path("LDR", "input", "A"),
            self.get_tmp_image_path("LDR", "decomp"),
            "4x4", "-fastt"]

        self.exec(command)

    def test_tl_bad_argument(self):
        """
        Test -tl with an unknown additional argument.
        """
        # Build an otherwise valid command with the test flaw
        command = [
            self.binary, "-tl",
            self.get_ref_image_path("LDR", "input", "A"),
            self.get_tmp_image_path("LDR", "decomp"),
            "4x4", "-fast", "-unknown"]

        self.exec(command)

    @unittest.skip("Bug #93")
    def test_dl_missing_args(self):
        """
        Test -dl with missing arguments.
        """
        # Build a valid command
        command = [
            self.binary, "-dl",
            self.get_ref_image_path("LDR", "comp", "A"),
            self.get_tmp_image_path("LDR", "decomp")]

        # Run the command, incrementally omitting arguments
        commandLen = len(command)
        for subLen in range(2, commandLen + 1):
            omit = len(command) - subLen
            with self.subTest(omit=omit):
                testCommand = command[:subLen]

                # For the last run we omit no arguments; make sure this works
                # to ensure that the underlying test is actually valid and
                # we're not failing for a different reason
                expectPass = omit == 0
                self.exec(testCommand, expectPass)

    @unittest.skip("Bug #93")
    def test_compare_no_args(self):
        """
        Test -compare with no arguments.
        """
        command = [self.binary, "-compare"]
        self.exec(command)


def main():
    """
    The main function.

    Returns:
        int: The process return code.
    """
    results = unittest.main(exit=False)
    return 0 if results.result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main())
