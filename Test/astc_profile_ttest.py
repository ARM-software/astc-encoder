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

The binaries are run in alternating order and their reported processing times
are analyzed with a paired, two-sided Student's t-test. The command line
options following the script options are passed through to astcenc unchanged.

Use an astcenc test mode (-tl, -ts, -th, or -tH) to report image quality.
Test image names (e.g. ldr-rgb-... or ldr-rgba-...) select RGB or RGBA PSNR;
other names report both when available. Alpha-weighted compression uses
alpha-weighted PSNR, HDR images use mPSNR, and normal maps also report mean
and worst angular errors. SSIM is included when both binaries report it.
Quality scores are averaged over measured runs,
excluding warmups; changes are reported as B minus A.
'''

import argparse
from pathlib import Path
import statistics
import sys

from scipy import stats

from testlib.profile import CodecResult
from testlib.profile import profile_binaries
from testlib.profile import PROFILE_ERRORS
from testlib.profile import QUALITY_METRICS
from testlib.profile import report_error


def select_quality_metrics(arguments: list[str],
                           available: dict[str, float]) -> list[str]:
    '''
    Choose quality metrics for the input image and compression settings.

    Args:
        arguments: The codec command-line arguments.
        available: The quality scores reported by the codec.

    Return:
        The relevant metric names, or an empty list outside test mode.
    '''
    if arguments[0] not in ('-tl', '-ts', '-th', '-tH'):
        return []

    if '-normal' in arguments:
        return ['PSNR (LDR-RGB)', 'Mean Angular Error', 'Worst Angular Error']

    image = Path(arguments[1])
    parts = image.stem.lower().split('-')
    if ('mPSNR (RGB)' in available or
            image.suffix.lower() in ('.hdr', '.exr') or
            (len(parts) >= 3 and parts[0] == 'hdr')):
        return ['mPSNR (RGB)']

    if '-a' in arguments:
        return ['Alpha-weighted PSNR']

    if len(parts) >= 3 and parts[0] in ('ldr', 'ldrs', 'hdr'):
        if parts[1] in ('rgb', 'l', 'xy', 'rgba', 'la'):
            channels = 'RGBA' if parts[1] in ('rgba', 'la') else 'RGB'
            return [f'PSNR (LDR-{channels})']

    # Generic filenames do not encode which channels are relevant. Show both
    # instead of silently including or ignoring alpha in a single score.
    return [name for name in ('PSNR (LDR-RGB)', 'PSNR (LDR-RGBA)')
            if name in available]


def report_quality(arguments: list[str], samples_a: list[CodecResult],
                   samples_b: list[CodecResult]) -> None:
    '''
    Report mean image quality for each binary and the signed change.

    Args:
        arguments: The codec command-line arguments.
        samples_a: The measured results for binary A.
        samples_b: The measured results for binary B.

    Raises:
        ValueError: A test-mode run did not report a required quality score.
    '''
    if arguments[0] not in ('-tl', '-ts', '-th', '-tH'):
        print('Quality: unavailable (use -tl, -ts, -th, or -tH test mode)')
        return

    metrics = select_quality_metrics(arguments, samples_a[0].quality)
    if not metrics:
        raise ValueError('Codec did not report image quality scores')

    # Older binaries do not report SSIM; compare it when both sides support it
    ssim_metrics = {
        'PSNR (LDR-RGB)': 'SSIM (LDR-RGB)',
        'PSNR (LDR-RGBA)': 'SSIM (LDR-RGBA)',
        'Alpha-weighted PSNR': 'SSIM (LDR-RGBA)',
        'mPSNR (RGB)': 'mSSIM (RGB)'
    }
    for metric in list(metrics):
        ssim = ssim_metrics.get(metric)
        if ssim and all(ssim in samples[0].quality
                        for samples in (samples_a, samples_b)):
            metrics.append(ssim)

    print('Image quality (mean over measured runs; B - A):')
    for metric in metrics:
        if any(metric not in sample.quality
               for sample in samples_a + samples_b):
            raise ValueError(f'Codec did not report {metric} in every run')

        mean_a = statistics.mean(
            sample.quality[metric] for sample in samples_a)
        mean_b = statistics.mean(
            sample.quality[metric] for sample in samples_b)
        unit, direction = QUALITY_METRICS[metric]
        precision = 4
        if not unit:
            precision = 6

        print(f'  {metric}: A {mean_a:.{precision}f}, '
              f'B {mean_b:.{precision}f}, '
              f'change {mean_b - mean_a:+.{precision}f} {unit} '
              f'({direction} is better)')


def paired_t_test(sample_a: list[float],
                  sample_b: list[float]) -> tuple[float, float]:
    '''
    Run a paired, two-sided Student's t-test.

    Args:
        sample_a: The timing samples for binary A.
        sample_b: The timing samples for binary B.

    Return:
        The t-statistic and two-sided p-value.

    Raises:
        ValueError: The sample sizes differ or contain fewer than two values.
    '''
    if len(sample_a) != len(sample_b) or len(sample_a) < 2:
        raise ValueError('Paired t-test requires equal sample set sizes')

    result = stats.ttest_rel(sample_a, sample_b)
    return float(result.statistic), float(result.pvalue)


def confidence_value(value: str) -> float:
    '''
    Parse and validate a confidence level argument.

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
                        help='confidence level (default: %(default)s)')

    args = parser.parse_args()

    if args.runs < 2:
        parser.error('--runs must be at least 2')

    if args.warmups < 0:
        parser.error('--warmups cannot be negative')

    if not args.arguments:
        parser.error('no astcenc arguments specified')

    return args


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
            args.warmups, args.runs)
        times_a = [sample.time for sample in samples_a]
        times_b = [sample.time for sample in samples_b]
        statistic, p_value = paired_t_test(times_a, times_b)
        report_quality(args.arguments, samples_a, samples_b)
    except PROFILE_ERRORS as ex:
        return report_error(ex)

    mean_a = statistics.mean(times_a)
    mean_b = statistics.mean(times_b)

    print()
    print(f'A mean: {mean_a:.6f} s ({binary_a})')
    print(f'B mean: {mean_b:.6f} s ({binary_b})')
    print(f'Paired t-test: t={statistic:.4f}, p={p_value:.6g}')

    if p_value >= 1.0 - args.confidence or mean_a == mean_b:
        print(f'Result: inconclusive at {args.confidence:.1%} confidence')
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
          f'{args.confidence:.1%} confidence '
          f'({reduction:.2f}% lower mean processing time)')

    return 0


if __name__ == '__main__':
    sys.exit(main())
