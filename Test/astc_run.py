#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# This confidential and proprietary software may be used only as authorised by
# a licensing agreement from Arm Limited.
#     (C) COPYRIGHT 2019-2020 Arm Limited, ALL RIGHTS RESERVED
# The entire notice above must be reproduced on all authorised copies and
# copies may only be made to the extent permitted by a licensing agreement from
# Arm Limited.
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

    def rewrite_args_for_old_cli(self, args):
        replacements = [
            ("-silent", "-silentmode")
        ]

        extensions = [
            ("-t", ("-showpsnr", "-time")),
            ("-tl", ("-showpsnr", "-time")),
            ("-ts", ("-showpsnr", "-time"))
        ]

        for new, old in replacements:
            args = [old if x == new else x for x in args]

        for new, exts in extensions:
            if new in args:
                args.extend(exts)

        return args

    def run_once(self, testBinary, blockSize, profile, verbose, newCLI):
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

        # Switch sRGB images into sRGB mode
        if self.format in ("srgb", "srgba"):
            opmode = "-ts"
        elif self.dynamicRange == "ldr":
            opmode = "-tl"
        else:
            opmode = "-t"

        # Run the compressor
        args = [testBinary, opmode, self.filePath, outFilePath,
                blockSize, "-thorough", "-silent"]

        # Switch normal maps into angular error metrics
        if self.format == "xy":
            args.append("-normal_psnr")

        # Switch HDR data formats into HDR compression mode; note that this
        # mode assumes that the alpha channel is non-correlated
        if self.dynamicRange == "hdr":
            args.append("-hdr")

        if verbose:
            args.append("-j")
            args.append("1")

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


        if not newCLI:
            args = self.rewrite_args_for_old_cli(args)

        try:
            result = sp.run(args, stdout=sp.PIPE, stderr=sp.PIPE,
                            check=True, universal_newlines=True)
        except sp.CalledProcessError as e:
            print("ERROR: Test run failed")
            print("  + %s" % " ".join(args))

            if verbose:
                print(e.stderr)
                print(e.stdout)

            sys.exit(1)

        # Convert the TGA to PNG and delete the TGA (LDR only)
        if self.dynamicRange == "ldr":
            im = Image.open(outFilePath)
            im.save(outFilePath2)
            os.remove(outFilePath)
            outFilePath = outFilePath2

        # Create log parsing patterns
        if self.dynamicRange == "ldr":
            if self.format in ("rgb", "xy", "l"):
                patternPSNR = r"PSNR \(LDR-RGB\):\s*([0-9.]*) dB"
            elif self.format in ("srgba", "rgba"):
                patternPSNR = r"PSNR \(LDR-RGBA\):\s*([0-9.]*) dB"
            else:
                assert False, "Unsupported LDR color format %s" % self.format
        else:
            patternPSNR = r"mPSNR \(RGB\)(?: \[.*?\] )?:\s*([0-9.]*) dB.*"

        patternPSNR = re.compile(patternPSNR)
        patternTime1 = re.compile(".*[Cc]oding time:\s*([0-9.]*) s.*")
        patternTime2 = re.compile(".*[Tt]otal time:\s*([0-9.]*) s.*")

        # Extract results from the log
        runPSNR = None
        runTime = None

        for line in result.stdout.splitlines():
            match = patternPSNR.match(line)
            if match:
                runPSNR = float(match.group(1))

            match = patternTime1.match(line)
            if match:
                runTime = float(match.group(1))

            match = patternTime2.match(line)
            if match:
                allTime = float(match.group(1))

        # Convert the callgrind log to an image
        if profile:
            os.system("gprof2dot ./callgrind.out -f callgrind -s -z \"astc_main(int, char*)\" > callgrind.dot")
            os.system("dot -Tpng callgrind.dot -o callgrind.png")
            os.remove("callgrind.dot")

        if verbose:
            print(result.stdout)

        assert runPSNR is not None, "No coding PSNR found %s" % result.stdout
        assert runTime is not None, "No coding time found %s" % result.stdout
        assert allTime is not None, "No total time found %s" % result.stdout
        return (runPSNR, runTime, allTime, outFilePath)


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


def get_binary(binary):
    """
    Return the reference binary path for the current host machine.
    """
    appmap = {
        "reference": "./Binary/linux-x64/astcenc",
        "new": "./Source/astcenc",
        "original": "./SourceOrig/astcenc",
        "prototype": "./SourceProto/astcenc",
    }

    return appmap[binary]


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

    choices = ("reference", "original", "prototype", "new")
    parser.add_argument("--binary", dest="useBinary", default="new",
                        choices=choices, help="select binary")

    parser.add_argument("--repeats", dest="testRepeats", default=1,
                        type=int, help="test iteration count")

    parser.add_argument("--warmup", dest="testWarmups", default=0,
                        type=int, help="warmup iteration count")

    parser.add_argument("-v", dest="verbose", default=False,
                        action="store_true", help="enable verbose logging")

    args = parser.parse_args()
    return args


def main():
    """
    The main function.
    """
    # Parse command lines
    args = parse_command_line()

    app = get_binary(args.useBinary)

    image = TestImage(args.image)

    # Run the test scenario
    allTimes = []
    codeTimes = []
    for _ in range(0, args.testRepeats + args.testWarmups):
        psnr, codeSecs, allSecs, output = \
            image.run_once(app, args.block, args.profile, args.verbose,
            args.useBinary == "new")
        codeTimes.append(codeSecs)
        allTimes.append(allSecs)


    # Strip off the warmup times and average the rest
    allTimes = allTimes[args.testWarmups:]
    allSecs = sum(allTimes) / len(allTimes)

    codeTimes = codeTimes[args.testWarmups:]
    codeSecs = sum(codeTimes) / len(codeTimes)

    # Print summary results
    print("Binary: %s" % args.useBinary)
    print("PSNR: %0.3f dB" % psnr)
    print("Coding time: %0.3f s (avg of %u runs) %0.3f s best" % (codeSecs, len(codeTimes), min(codeTimes)))
    print("All time:    %0.3f s (avg of %u runs) %0.3f s best" % (allSecs, len(allTimes), min(allTimes)))
    print("Image: %s" % output)


if __name__ == "__main__":
    sys.exit(main())
