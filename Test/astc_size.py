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
import shutil
import subprocess as sp
import sys


def get_test_binary():
    """
    Return the test binary path for the current host machine.
    """
    return "./Source/astcenc"


def get_reference_binary():
    """
    Return the reference binary path for the current host machine.
    """
    return "./Binary/linux-x64/astcenc"


def run_size(binary):
    args = ["size", binary]
    result = sp.run(args, stdout=sp.PIPE, stderr=sp.PIPE,
                    check=True, universal_newlines=True)

    lines = result.stdout.splitlines()
    assert len(lines) == 2
    values = lines[1].split()

    codeSection = float(values[0])
    roSection = float(values[1])
    ziSection = float(values[2])
    return (codeSection, roSection, ziSection)


def main():
    """
    The main function.
    """
    refApp = get_reference_binary()
    newApp = get_test_binary()

    refSize = run_size(refApp)
    newSize = run_size(newApp)

    print("%8s  % 8s  % 8s  % 8s" % ("", "Code", "RO", "ZI"))
    print("%8s  % 8u  % 8u  % 8u" % ("Ref", *refSize))
    print("%8s  % 8u  % 8u  % 8u" % ("New", *newSize))

    diffAbs = []
    diffRel = []
    for a, b in zip(refSize, newSize):
        diff = b - a
        diffAbs.append(diff)
        rel = diff / a
        diffRel.append((diff / a) * 100.0)

    print("%8s  % 8u  % 8u  % 8u" % ("Abs D", *diffAbs))
    print("%8s  % 7.2f%%  % 7.2f%%  % 7.2f%%" % ("Rel D", *diffRel))


if __name__ == "__main__":
    sys.exit(main())
