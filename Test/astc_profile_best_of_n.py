#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# -----------------------------------------------------------------------------
# Copyright 2026 Arm Limited
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
'''
Utility script to compare the best performance of two astcenc binaries.

The binaries are run in alternating order and the best reported processing time for each binary is
retained. The command line options following the script options are passed through to astcenc
unchanged.
'''

import argparse
import sys

from testlib.cli import parse_comparison_args
from testlib.profile import profile_binaries
from testlib.profile import PROFILE_ERRORS
from testlib.profile import report_error


def parse_command_line() -> argparse.Namespace:
    '''
    Parse the binary comparison command line.

    Return:
        The binary paths, codec arguments, and validated sampling settings.
    '''

    return parse_comparison_args(__doc__, is_ttest=False)


def main() -> int:
    '''
    Run the comparison and print the best result.

    Return:
        The process return code.
    '''
    args = parse_command_line()

    try:
        binary_a, binary_b, samples_a, samples_b = profile_binaries(
            args.binary_a, args.binary_b, args.arguments,
            args.warmups, args.runs)
    except PROFILE_ERRORS as ex:
        return report_error(ex)

    best_a = min(sample.time for sample in samples_a)
    best_b = min(sample.time for sample in samples_b)

    print()
    print(f'A best: {best_a:.6f} s ({binary_a})')
    print(f'B best: {best_b:.6f} s ({binary_b})')

    if best_a == best_b:
        print('Result: equal best processing time')

        return 0

    if best_a < best_b:
        winner = 'A'
        faster_time = best_a
        slower_time = best_b
    else:
        winner = 'B'
        faster_time = best_b
        slower_time = best_a

    reduction = (slower_time - faster_time) / slower_time * 100.0
    print(f'Result: {winner} is faster ({reduction:.2f}% lower best processing time)')

    return 0


if __name__ == '__main__':
    sys.exit(main())
