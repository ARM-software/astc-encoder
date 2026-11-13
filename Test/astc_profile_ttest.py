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
Utility script to compare the performance and quality of two astcenc binaries.

The binaries are run in alternating order and their reported processing times are analyzed with a
paired, two-sided Student's t-test. The command line options following the script options are passed
through to astcenc unchanged. After at least two measured pairs, stop when the paired test reaches
the requested confidence level (95% by default), or the --repeats limit is reached. Repeated checks
use nominal confidence, not a sequentially adjusted guarantee.

Use an astcenc test mode (-tl, -ts, -th, or -tH) to report image quality. Test image names (e.g.
ldr-rgb-... or ldr-rgba-...) select RGB or RGBA PSNR; other names report both when available.
Alpha-weighted compression uses alpha-weighted PSNR, HDR images use mPSNR, and normal maps also
report mean and worst angular errors. SSIM is included when both binaries report it. Quality scores
are averaged over measured runs, excluding warmups; changes are reported as B minus A.
'''

import argparse
import statistics
import sys

from testlib.cli import parse_comparison_args
from testlib.timing import paired_t_test
from testlib.profile import profile_binaries
from testlib.profile import PROFILE_ERRORS
from testlib.reporting import report_quality
from testlib.profile import report_error


def parse_command_line() -> argparse.Namespace:
    '''
    Parse the paired t-test command line.

    Return:
        The binary paths, codec arguments, and validated sampling settings.
    '''

    return parse_comparison_args(__doc__, is_ttest=True)


def main() -> int:
    '''
    Run the comparison and print the statistical result.

    Return:
        The process return code.
    '''
    args = parse_command_line()

    try:
        binary_a, binary_b, samples_a, samples_b = profile_binaries(
            args.binary_a, args.binary_b, args.arguments,
            args.warmups, args.repeats, args.confidence)
        times_a = [sample.time for sample in samples_a]
        times_b = [sample.time for sample in samples_b]
        statistic, p_value = paired_t_test(times_a, times_b)
        report_quality(args.arguments, samples_a, samples_b)
    except PROFILE_ERRORS as ex:
        return report_error(ex)

    mean_a = statistics.mean(times_a)
    mean_b = statistics.mean(times_b)

    print()
    print(f'Measured pairs: {len(samples_a)} (maximum {args.repeats})')
    print(f'A mean: {mean_a:.6f} s ({binary_a})')
    print(f'B mean: {mean_b:.6f} s ({binary_b})')
    print(f'Paired t-test: t={statistic:.4f}, p={p_value:.6g}')

    if p_value >= 1.0 - args.confidence or mean_a == mean_b:
        print(f'Result: inconclusive at {args.confidence:.1%} nominal confidence')

        return 0

    if mean_a < mean_b:
        winner = 'A'
        faster_time = mean_a
        slower_time = mean_b
    else:
        winner = 'B'
        faster_time = mean_b
        slower_time = mean_a

    reduction = (slower_time - faster_time) / slower_time * 100.0
    print(f'Result: {winner} is faster at '
          f'{args.confidence:.1%} nominal confidence '
          f'({reduction:.2f}% lower mean processing time)')

    return 0


if __name__ == '__main__':
    sys.exit(main())
