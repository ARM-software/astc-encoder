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

import argparse
import os
from PIL import Image
import re
import shutil
import subprocess as sp
import sys


class TestImage():
    """
    A single test definition, and the test results if it is actually run.
    """

    def __init__(self, filePath):
        """
        Construct a new test definition.
        """
        # Remove any leading "./" or ".\" prefix
        if filePath.startswith("./") or filePath.startswith(".\\"):
            filePath = filePath[2:]

        name = os.path.basename(filePath)
        nameParts = name.split("-")

        self.filePath = filePath
        self.dynamicRange = nameParts[0]
        self.format = nameParts[1]

    def run_once(self, testBinary, blockSize, profile):
        """
        Run a single compression pass.
        """
        pathParts = splitall(self.filePath)

        # Create the test output directory if it doesn't exist
        outDir = os.path.join("TestOutput", "Images", pathParts[2], blockSize)
        os.makedirs(outDir,  exist_ok=True)

        if self.dynamicRange == "ldr":
            outFile = pathParts[-1].replace(".png", ".tga")
        else:
            outFile = pathParts[-1].replace(".hdr", ".htga")
        outFilePath = os.path.join(outDir, outFile)

        if self.dynamicRange == "ldr":
            outFile = pathParts[-1].replace(".png", "-out.png")
            outFilePath2 = os.path.join(outDir, outFile)

        # Run the compressor
        args = [testBinary, "-t", self.filePath, outFilePath,
                blockSize, "-thorough", "-time", "-showpsnr", "-silentmode"]

        # Switch normal maps into angular error metrics
        if self.format == "xy":
            args.append("-normal_psnr")

        # Inject the callgrind profiler prefix
        if profile:
            args.insert(0, "--callgrind-out-file=callgrind.out")
            args.insert(0, "--tool=callgrind")
            args.insert(0, "valgrind")
            if os.path.exists("callgrind.out"):
                os.remove("callgrind.out")
            if os.path.exists("callgrind.dot"):
                os.remove("callgrind.dot")
            if os.path.exists("callgrind.png"):
                os.remove("callgrind.png")

        try:
            result = sp.run(args, stdout=sp.PIPE, stderr=sp.PIPE,
                            check=True, universal_newlines=True)
        except (OSError, sp.CalledProcessError):
            print("ERROR: Test run failed")
            print("  + %s" % " ".join(args))
            sys.exit(1)

        # Convert the TGA to PNG and delete the TGA (LDR only)
        if self.dynamicRange == "ldr":
            im = Image.open(outFilePath)
            im.save(outFilePath2)
            os.remove(outFilePath)
            outFilePath = outFilePath2

        # Create log parsing patterns
        if self.dynamicRange == "ldr":
            if self.format in ("rgb", "xy"):
                patternPSNR = "PSNR \\(LDR-RGB\\): ([0-9.]*) dB"
            elif self.format == "rgba":
                patternPSNR = "PSNR \\(LDR-RGBA\\): ([0-9.]*) dB"
            else:
                assert False, "Unsupported LDR color format %s" % self.format
        else:
            patternPSNR = "PSNR \\(RGB normalized to peak\\): ([0-9.]*) dB"

        patternPSNR = re.compile(patternPSNR)
        patternTime = re.compile(".* coding time: ([0-9.]*) seconds")

        # Extract results from the log
        runPSNR = None
        runTime = None

        for line in result.stdout.splitlines():
            match = patternPSNR.match(line)
            if match:
                runPSNR = float(match.group(1))

            match = patternTime.match(line)
            if match:
                runTime = float(match.group(1))

        # Convert the callgrind log to an image
        if profile:
            os.system("gprof2dot ./callgrind.out -f callgrind -s -z \"astc_main(int, char*)\" > callgrind.dot")
            os.system("dot -Tpng callgrind.dot -o callgrind.png")
            os.remove("callgrind.dot")

        assert runPSNR is not None, "No coding PSNR found %s" % result.stdout
        assert runTime is not None, "No coding time found %s" % result.stdout
        return (runPSNR, runTime, outFilePath)


def splitall(path):
    """
    Completely tokenize a path into its component pieces.
    """
    allparts = []
    while True:
        parts = os.path.split(path)
        if parts[0] == path:
            allparts.insert(0, parts[0])
            break
        elif parts[1] == path:
            allparts.insert(0, parts[1])
            break
        else:
            path = parts[0]
            allparts.insert(0, parts[1])

    return allparts


def get_test_binary():
    """
    Return the test binary path for the current host machine.
    """
    if "linux" in sys.platform:
        return "./Source/astcenc"
    elif sys.platform == "darwin":
        return "./Source/astcenc"
    elif sys.platform == "win32":
        return "./Source/VS2017/Release/astcenc.exe"

    assert False, "Unknown operating system %s" % sys.platform


def get_reference_binary():
    """
    Return the reference binary path for the current host machine.
    """
    if "linux" in sys.platform:
        return "./Binary/linux-x64/astcenc"
    elif sys.platform == "darwin":
        return "./Binary/mac-x64/astcenc"
    elif sys.platform == "win32":
        return "./Binary/windows-x64/astcenc.exe"

    assert False, "Unknown operating system %s" % sys.platform


def parse_command_line():
    """
    Parse the command line.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("image", help="image to compress")

    parser.add_argument("block", nargs="?", default="6x6",
                        choices=["4x4", "5x5", "6x6", "8x8"],
                        help="block size")

    parser.add_argument("--profile", dest="profile", default=False,
                        action="store_true", help="use callgrind profiler")

    parser.add_argument("--reference", dest="useReference", default=False,
                        action="store_true", help="use reference binary")

    parser.add_argument("--repeats", dest="testRepeats", default=1,
                        type=int, help="test iteration count")

    parser.add_argument("--warmup", dest="testWarmups", default=0,
                        type=int, help="warmup iteration count")

    args = parser.parse_args()
    return args


def main():
    """
    The main function.
    """
    # Parse command lines
    args = parse_command_line()

    app = get_test_binary()
    if args.useReference:
        app = get_reference_binary()

    image = TestImage(args.image)

    # Run the test scenario
    times = []
    for _ in range(0, args.testRepeats + args.testWarmups):
        psnr, secs, output = image.run_once(app, args.block, args.profile)
        times.append(secs)

    # Strip off the warmup times and average the rest
    times = times[args.testWarmups:]
    secs = sum(times) / len(times)

    # Print summary results
    print("PSNR: %0.3f dB" % psnr)
    print("Time: %0.3f s (avg of %u runs)" % (secs, len(times)))
    print("Image: %s" % output)


if __name__ == "__main__":
    sys.exit(main())
