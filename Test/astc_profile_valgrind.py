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
A simple wrapper utility to run a callgrind profile over a test image, and
post-process the output into an call graph image.

Only runs on Linux and requires the following tools available on the PATH:

  * valgrind
  * gprof2dot
  * dot
"""


import argparse
import os
import subprocess as sp
import sys

def run_pass(image, encoder, blocksize, quality):
    """
    Run Valgrind on a single binary.

    Args:
        image (str): The path of the image to compress.
        encoder (str): The name of the encoder variant to run.
        blocksize (str): The block size to use.
        quality (str): The encoding quality to use.

    Raises:
        CalledProcessException: Any subprocess failed.
    """
    binary =  "./Source/astcenc-%s" % encoder
    qualityFlag = "-%s" % quality
    args = ["valgrind", "--tool=callgrind", "--callgrind-out-file=callgrind.txt",
            binary, "-cl", image, "out.astc", blocksize, qualityFlag, "-j", "1"]

    print(" ".join(args))
    result = sp.run(args, check=True, universal_newlines=True)

    args = ["gprof2dot", "--format=callgrind", "--output=out.dot", "callgrind.txt",
            "-s", "-z", "compress_block(astcenc_context const&, astcenc_image const&, imageblock const*, symbolic_compressed_block&, physical_compressed_block&, compress_symbolic_block_buffers*)"]

    result = sp.run(args, check=True, universal_newlines=True)

    args = ["dot", "-Tpng", "out.dot", "-o", "%s.png" % quality]
    result = sp.run(args, check=True, universal_newlines=True)

    os.remove("out.astc")
    os.remove("out.dot")
    os.remove("callgrind.txt")

def parse_command_line():
    """
    Parse the command line.

    Returns:
        Namespace: The parsed command line container.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("img", type=argparse.FileType("r"),
                        help="The image file to test")

    testencoders = ["sse2", "sse4.2", "avx2"]
    encoders = testencoders + ["all"]
    parser.add_argument("--encoder", dest="encoders", default="avx2",
                        choices=encoders, help="select encoder variant")

    testqualities = ["fast", "medium", "thorough"]
    qualities = testqualities + ["all"]
    parser.add_argument("--test-quality", dest="qualities", default="medium",
                        choices=qualities, help="select compression quality")

    args = parser.parse_args()

    if args.encoders == "all":
        args.encoders = testencoders
    else:
        args.encoders = [args.encoders]

    if args.qualities == "all":
        args.qualities = testqualities
    else:
        args.qualities = [args.qualities]

    return args


def main():
    """
    The main function.

    Returns:
        int: The process return code.
    """
    args = parse_command_line()

    for quality in args.qualities:
        for encoder in args.encoders:
            run_pass(args.img.name, encoder, "6x6", quality)

    return 0


if __name__ == "__main__":
    sys.exit(main())
