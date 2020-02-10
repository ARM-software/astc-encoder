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
