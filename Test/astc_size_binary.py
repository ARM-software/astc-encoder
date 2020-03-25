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
Utility to show the section sizes of a binary, or compare the section sizes of
two binaries. This script uses the Linux "size" utility, which must be on the
PATH.
"""


import argparse
import subprocess as sp
import sys


def run_size(binary):
    """
    Run size on a single binary.

    Args:
        binary: The path.

    Returns:
        A tuple of (code size, read-only data size, and zero-init data size).
    """
    args = ["size", "--format=sysv", binary]
    result = sp.run(args, stdout=sp.PIPE, stderr=sp.PIPE,
                    check=True, universal_newlines=True)

    data = {}
    patterns = {"Code": ".text", "RO": ".rodata", "ZI": ".bss"}

    lines = result.stdout.splitlines()
    for line in lines:
        for key, value in patterns.items():
            if line.startswith(value):
                size = float(line.split()[1])
                data[key] = size

    return (data["Code"], data["RO"], data["ZI"])


def parse_command_line():
    """
    Parse the command line.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("bin", type=argparse.FileType("r"),
                        help="The new binary to size")
    parser.add_argument("ref", nargs="?", type=argparse.FileType("r"),
                        help="The reference binary to compare against")

    return parser.parse_args()


def main():
    """
    The main function.
    """
    args = parse_command_line()

    # Collect the data
    newSize = run_size(args.bin.name)
    if args.ref:
        refSize = run_size(args.ref.name)

    # Print the basic table of absolute values
    print("%8s  % 8s  % 8s  % 8s" % ("", "Code", "RO", "ZI"))
    if args.ref:
        print("%8s  % 8u  % 8u  % 8u" % ("Ref", *refSize))
    print("%8s  % 8u  % 8u  % 8u" % ("New", *newSize))

    # Print the difference if we have a reference
    if args.ref:
        diffAbs = []
        diffRel = []
        for refVal, newVal in zip(refSize, newSize):
            diff = newVal - refVal
            diffAbs.append(diff)
            diffRel.append((diff / refVal) * 100.0)

        dat = ("Abs D", diffAbs[0], diffAbs[1], diffAbs[2])
        print("%8s  % 8u  % 8u  % 8u" % dat)
        dat = ("Rel D", diffRel[0], diffRel[1], diffRel[2])
        print("%8s  % 7.2f%%  % 7.2f%%  % 7.2f%%" % dat)

    return 0


if __name__ == "__main__":
    sys.exit(main())
