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
These classes provide an abstraction around the astcenc command line tool,
allowing the rest of the image test suite to ignore changes in the command line
interface.
"""

import os
import re
import subprocess as sp
import sys


class EncoderBase():
    """
    This class is a Python wrapper for the `astcenc` binary, providing an
    abstract means to set command line options and parse key results.

    This is an abstract base class providing some generic helper functionality
    used by concrete instantiations of subclasses.

    Attributes:
        binary: The encoder binary path.
        variant: The encoder SIMD variant being tested.
        name: The encoder name to use in reports.
        VERSION: The encoder version or branch.
        SWITCHES: Dict of switch replacements for different color formats.
        OUTPUTS: Dict of output file extensions for different color formats.
    """

    VERSION = None
    SWITCHES = None
    OUTPUTS = None

    def __init__(self, name, variant, binary):
        """
        Create a new encoder instance.

        Args:
            name (str): The name of the encoder.
            variant (str): The SIMD variant of the encoder.
            binary (str): The path to the binary on the file system.
        """
        self.name = name
        self.variant = variant
        self.binary = binary

    def build_cli(self, image, blockSize="6x6", preset="-thorough"):
        """
        Build the command line needed for the given test.

        Args:
            image (TestImage): The test image to compress.
            blockSize (str): The block size to use.
            preset (str): The quality-performance preset to use.

        Returns:
            list(str): A list of command line arguments.
        """
        # pylint: disable=unused-argument,no-self-use,redundant-returns-doc
        assert False, "Missing subclass implementation"

    def execute(self, command):
        """
        Run a subprocess with the specified command.

        Args:
            command (list(str)): The list of command line arguments.

        Returns:
            list(str): The output log (stdout) split into lines.
        """
        # pylint: disable=no-self-use
        try:
            result = sp.run(command, stdout=sp.PIPE, stderr=sp.PIPE,
                            check=True, universal_newlines=True)
        except (OSError, sp.CalledProcessError):
            print("ERROR: Test run failed")
            print("  + %s" % " ".join(command))
            sys.exit(1)

        return result.stdout.splitlines()

    def parse_output(self, image, output):
        """
        Parse the log output for PSNR and performance metrics.

        Args:
            image (TestImage): The test image to compress.
            output (list(str)): The output log from the compression process.

        Returns:
            tuple(float, float, float): PSNR in dB, TotalTime in seconds, and
            CodingTime in seconds.
        """
        # Regex pattern for image quality
        patternPSNR = re.compile(self.get_psnr_pattern(image))
        patternTTime = re.compile(self.get_total_time_pattern())
        patternCTime = re.compile(self.get_coding_time_pattern())

        # Extract results from the log
        runPSNR = None
        runTTime = None
        runCTime = None

        for line in output:
            match = patternPSNR.match(line)
            if match:
                runPSNR = float(match.group(1))

            match = patternTTime.match(line)
            if match:
                runTTime = float(match.group(1))

            match = patternCTime.match(line)
            if match:
                runCTime = float(match.group(1))

        stdout = "\n".join(output)
        assert runPSNR is not None, "No coding PSNR found %s" % stdout
        assert runTTime is not None, "No total time found %s" % stdout
        assert runCTime is not None, "No coding time found %s" % stdout
        return (runPSNR, runTTime, runCTime)

    def get_psnr_pattern(self, image):
        """
        Get the regex pattern to match the image quality metric.

        Note, while this function is called PSNR for some images we may choose
        to match another metric (e.g. mPSNR for HDR images).

        Args:
            image (TestImage): The test image we are compressing.

        Returns:
            str: The string for a regex pattern.
        """
        # pylint: disable=unused-argument,no-self-use,redundant-returns-doc
        assert False, "Missing subclass implementation"

    def get_total_time_pattern(self):
        """
        Get the regex pattern to match the total compression time.

        Returns:
            str: The string for a regex pattern.
        """
        # pylint: disable=unused-argument,no-self-use,redundant-returns-doc
        assert False, "Missing subclass implementation"

    def get_coding_time_pattern(self):
        """
        Get the regex pattern to match the coding compression time.

        Returns:
            str: The string for a regex pattern.
        """
        # pylint: disable=unused-argument,no-self-use,redundant-returns-doc
        assert False, "Missing subclass implementation"

    def run_test(self, image, blockSize, preset, testRuns):
        """
        Run the test N times.

        Args:
            image (TestImage): The test image to compress.
            blockSize (str): The block size to use.
            preset (str): The quality-performance preset to use.
            testRuns (int): The number of test runs.

        Returns:
            tuple(float, float, float): Returns the best results from the N
            test runs, as PSNR (dB), total time (seconds), and coding time
            (seconds).
        """
        # pylint: disable=assignment-from-no-return
        command = self.build_cli(image, blockSize, preset)

        # Execute test runs
        bestPSNR = 0
        bestTTime = sys.float_info.max
        bestCTime = sys.float_info.max

        for _ in range(0, testRuns):
            output = self.execute(command)
            result = self.parse_output(image, output)

            # Keep the best results (highest PSNR, lowest times)
            bestPSNR = max(bestPSNR, result[0])
            bestTTime = min(bestTTime, result[1])
            bestCTime = min(bestCTime, result[2])

        return (bestPSNR, bestTTime, bestCTime)


class Encoder2x(EncoderBase):
    """
    This class wraps the latest `astcenc` 2.x series binaries from the master
    branch.
    """
    VERSION = "master"

    SWITCHES = {
        "ldr": "-tl",
        "ldrs": "-ts",
        "hdr": "-th",
        "hdra": "-tH"
    }

    OUTPUTS = {
        "ldr": ".png",
        "ldrs": ".png",
        "hdr": ".exr",
        "hdra": ".exr"
    }

    def __init__(self, variant):
        name = "astcenc-%s-%s" % (variant, self.VERSION)
        if os.name == 'nt':
            dat = (variant, variant)
            binary = "./Source/VS2019/astcenc-%s-Release/astcenc-%s.exe" % dat
        else:
            binary = "./Source/astcenc-%s" % variant
        super().__init__(name, variant, binary)

    def build_cli(self, image, blockSize="6x6", preset="-thorough"):
        opmode = self.SWITCHES[image.colorProfile]
        srcPath = image.filePath

        dstPath = image.outFilePath + self.OUTPUTS[image.colorProfile]
        dstDir = os.path.dirname(dstPath)
        dstFile = os.path.basename(dstPath)
        dstPath = os.path.join(dstDir, self.name, blockSize, dstFile)

        dstDir = os.path.dirname(dstPath)
        os.makedirs(dstDir, exist_ok=True)

        command = [
            self.binary, opmode, srcPath, dstPath,
            blockSize, preset, "-silent"
        ]

        if image.colorFormat == "xy":
            command.append("-normal_psnr")

        if image.isMask:
            command.append("-mask")

        if image.isAlphaScaled:
            command.append("-a")
            command.append("1")

        return command

    def get_psnr_pattern(self, image):
        if image.colorProfile != "hdr":
            if image.colorFormat != "rgba":
                patternPSNR = r"PSNR \(LDR-RGB\):\s*([0-9.]*) dB"
            else:
                patternPSNR = r"PSNR \(LDR-RGBA\):\s*([0-9.]*) dB"
        else:
            patternPSNR = r"mPSNR \(RGB\)(?: \[.*?\] )?:\s*([0-9.]*) dB.*"
        return patternPSNR

    def get_total_time_pattern(self):
        return r"Total time:\s*([0-9.]*) s"

    def get_coding_time_pattern(self):
        return r"Coding time:\s*([0-9.]*) s"


class Encoder1x(EncoderBase):
    """
    This class wraps the latest `astcenc` 1.x series binaries from the 1.x
    branch.
    """
    VERSION = "1.7"

    SWITCHES = {
        "ldr": "-tl",
        "ldrs": "-ts",
        "hdr": "-t"
    }

    OUTPUTS = {
        "ldr": ".tga",
        "ldrs": ".tga",
        "hdr": ".htga"
    }

    def __init__(self, binary=None):
        name = "astcenc-%s" % self.VERSION
        if not binary:
            if os.name == 'nt':
                binary = "./Binaries/1.7/astcenc.exe"
            else:
                binary = "./Binaries/1.7/astcenc"
        super().__init__(name, None, binary)

    def build_cli(self, image, blockSize="6x6", preset="-thorough"):
        opmode = self.SWITCHES[image.colorProfile]
        srcPath = image.filePath

        dstPath = image.outFilePath + self.OUTPUTS[image.colorProfile]
        dstDir = os.path.dirname(dstPath)
        dstFile = os.path.basename(dstPath)
        dstPath = os.path.join(dstDir, self.name, blockSize, dstFile)

        dstDir = os.path.dirname(dstPath)
        os.makedirs(dstDir, exist_ok=True)

        command = [
            self.binary, opmode, srcPath, dstPath,
            blockSize, preset, "-silentmode", "-time", "-showpsnr"
        ]

        if image.colorFormat == "xy":
            command.append("-normal_psnr")

        if image.colorProfile == "hdr":
            command.append("-hdr")

        if image.isMask:
            command.append("-mask")

        if image.isAlphaScaled:
            command.append("-alphablend")

        return command

    def get_psnr_pattern(self, image):
        if image.colorProfile != "hdr":
            if image.colorFormat != "rgba":
                patternPSNR = r"PSNR \(LDR-RGB\):\s*([0-9.]*) dB"
            else:
                patternPSNR = r"PSNR \(LDR-RGBA\):\s*([0-9.]*) dB"
        else:
            patternPSNR = r"mPSNR \(RGB\)(?: \[.*?\] )?:\s*([0-9.]*) dB.*"
        return patternPSNR

    def get_total_time_pattern(self):
        return r"Elapsed time:\s*([0-9.]*) seconds.*"

    def get_coding_time_pattern(self):
        return r".* coding time: \s*([0-9.]*) seconds"


class EncoderProto(Encoder1x):
    """
    This class wraps a prototype `astcenc` binary variant of the 1.x branch.
    """

    VERSION = "Proto"

    OUTPUTS = {
        "ldr": ".png",
        "ldrs": ".png",
        "hdr": ".exr",
        "hdra": ".exr"
    }

    def __init__(self):
        name = "astcenc-%s" % self.VERSION
        assert os.name != 'nt', "Windows builds not available"
        binary = "./Binaries/Prototype/astcenc"
        super().__init__(binary)
