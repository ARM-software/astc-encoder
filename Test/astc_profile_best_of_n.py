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

The binaries are run in alternating order and the best reported processing
time for each binary is retained. The command line options following the
script options are passed through to astcenc unchanged.
'''

# Pylint reports duplicate-code at module scope, so it cannot be suppressed
# only on parse_command_line(). Keep the command-line interface local and
# aligned with astc_profile_ttest.py.
# pylint: disable=duplicate-code

import argparse
import sys

from testlib.profile import profile_binaries
from testlib.profile import PROFILE_ERRORS
from testlib.profile import report_error


def confidence_value(value: str) -> float:
    '''
    Parse and validate a confidence level argument.

    This option is accepted for command-line compatibility with
    astc_profile_ttest.py, but is not otherwise used.

    Args:
        value: The confidence level string.

    Return:
        The parsed confidence level.

    Raises:
        ArgumentTypeError: The confidence level is outside the valid range.
    '''
    result = float(value)
    if not 0.0 < result < 1.0:
        raise argparse.ArgumentTypeError('Must be between 0 and 1')
    return result


def parse_command_line() -> argparse.Namespace:
    '''
    Parse the command line.

    Return:
        The parsed command-line arguments.
    '''
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('binary_a', help='the first astcenc binary')

    parser.add_argument('binary_b', help='the second astcenc binary')

    parser.add_argument('arguments', nargs=argparse.REMAINDER,
                        help='arguments passed through to both binaries')

    parser.add_argument('--runs', type=int, default=10,
                        help='number of measured pairs (default: %(default)s)')

    parser.add_argument('--warmups', type=int, default=1,
                        help='number of warm-up pairs (default: %(default)s)')

    parser.add_argument('--confidence', type=confidence_value, default=0.95,
                        help='accepted for compatibility; otherwise ignored')

    args = parser.parse_args()

    if args.runs < 1:
        parser.error('--runs must be at least 1')

    if args.warmups < 0:
        parser.error('--warmups cannot be negative')

    if not args.arguments:
        parser.error('no astcenc arguments specified')

    return args


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
    print(f'Result: {winner} is faster '
          f'({reduction:.2f}% lower best processing time)')

    return 0


if __name__ == '__main__':
    sys.exit(main())
