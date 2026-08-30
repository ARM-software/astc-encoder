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
Utility script to compare the compression performance of two astcenc binaries.

The binaries are run in alternating order and their reported coding times are
analyzed with a paired, two-sided Student's t-test.
'''

import argparse
from pathlib import Path
import re
import shutil
import statistics
import subprocess as sp
import sys

from scipy import stats


CODING_TIME_PATTERN = re.compile(r'^\s*Coding time:\s*([0-9.]+)\s+s\s*$')
QUALITY_LEVELS = [
    'fastest', 'fast', 'medium', 'thorough', 'verythorough', 'exhaustive'
]


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


def resolve_binary(binary: str) -> str:
    '''
    Resolve a binary path, including commands found on PATH.

    Args:
        binary: The binary path or command name.

    Return:
        The resolved binary path.

    Raises:
        FileNotFoundError: The binary cannot be found.
    '''
    path = Path(binary)
    if path.is_file():
        return str(path)

    resolved = shutil.which(binary)
    if not resolved:
        raise FileNotFoundError(binary)

    return resolved


def run_encoder(binary: str, image_in: str, image_out: str,
                block_size: str, quality: str) -> float:
    '''
    Run one astcenc test-mode compression.

    Args:
        binary: The astcenc binary path.
        image_in: The input image path.
        image_out: The output image path.
        block_size: The ASTC block size.
        quality: The compression quality level.

    Return:
        The reported coding time in seconds.

    Raises:
        CalledProcessError: The encoder process failed.
        ValueError: The encoder did not report a coding time.
    '''
    command = [
        binary, '-tl', image_in, image_out, block_size, f'-{quality}'
    ]
    result = sp.run(command, stdout=sp.PIPE, stderr=sp.PIPE,
                    check=True, text=True)

    for line in result.stdout.splitlines():
        if match := CODING_TIME_PATTERN.match(line):
            return float(match.group(1))

    raise ValueError(f'{binary} did not report a coding time')


def collect_samples(binary_a: str, binary_b: str, image_in: str,
                    image_out: str, block_size: str, quality: str,
                    warmups: int,
                    runs: int) -> tuple[list[float], list[float]]:
    '''
    Warm up the encoders and collect interleaved paired timing samples.

    Args:
        binary_a: The binary A path.
        binary_b: The binary B path.
        image_in: The input image path.
        image_out: The output image path.
        block_size: The ASTC block size.
        quality: The compression quality level.
        warmups: The number of warm-up pairs.
        runs: The number of measured pairs.

    Return:
        The timing samples for binary A and binary B.
    '''
    for _ in range(warmups):
        run_encoder(binary_a, image_in, image_out, block_size, quality)
        run_encoder(binary_b, image_in, image_out, block_size, quality)

    samples_a = []
    samples_b = []
    for run in range(runs):
        if run % 2 == 0:
            time_a = run_encoder(binary_a, image_in, image_out,
                                 block_size, quality)
            time_b = run_encoder(binary_b, image_in, image_out,
                                 block_size, quality)
        else:
            time_b = run_encoder(binary_b, image_in, image_out,
                                 block_size, quality)
            time_a = run_encoder(binary_a, image_in, image_out,
                                 block_size, quality)

        samples_a.append(time_a)
        samples_b.append(time_b)
        print(f'Pair {run + 1:>2}/{runs}: '
              f'A {time_a:.4f} s, '
              f'B {time_b:.4f} s')

    return samples_a, samples_b


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
    '''Parse the command line.

    Return:
        The parsed command-line arguments.
    '''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('binary_a', help='the first astcenc binary')
    parser.add_argument('binary_b', help='the second astcenc binary')
    parser.add_argument('input', help='the input image')
    parser.add_argument('output', help='the output image')
    parser.add_argument('--block-size', default='6x6',
                        help='the ASTC block size (default: %(default)s)')
    parser.add_argument('--quality', choices=QUALITY_LEVELS, default='medium',
                        help='the compression quality (default: %(default)s)')
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
    if not Path(args.input).is_file():
        parser.error(f'input image does not exist: {args.input}')
    return args


def main() -> int:
    '''Run the comparison and print the statistical result.

    Return:
        The process return code.
    '''
    args = parse_command_line()

    try:
        binary_a = resolve_binary(args.binary_a)
        binary_b = resolve_binary(args.binary_b)
        samples_a, samples_b = collect_samples(
            binary_a, binary_b, args.input, args.output,
            args.block_size, args.quality, args.warmups, args.runs)
        statistic, p_value = paired_t_test(samples_a, samples_b)

    except FileNotFoundError as ex:
        print(f'ERROR: Binary not found: {ex}', file=sys.stderr)
        return 1

    except OSError as ex:
        print(f'ERROR: Could not execute encoder: {ex}', file=sys.stderr)
        return 1

    except sp.CalledProcessError as ex:
        print(f'ERROR: Encoder failed: {" ".join(ex.cmd)}', file=sys.stderr)
        output = ex.stderr.strip() or ex.stdout.strip()
        if output:
            print(output, file=sys.stderr)
        return 1

    except ValueError as ex:
        print(f'ERROR: {ex}', file=sys.stderr)
        return 1

    mean_a = statistics.mean(samples_a)
    mean_b = statistics.mean(samples_b)
    alpha = 1.0 - args.confidence

    print()
    print(f'A mean: {mean_a:.6f} s ({binary_a})')
    print(f'B mean: {mean_b:.6f} s ({binary_b})')
    print(f'Paired t-test: t={statistic:.4f}, p={p_value:.6g}')

    if p_value >= alpha or mean_a == mean_b:
        print(f'Result: inconclusive at {args.confidence:.1%} confidence')

    else:
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
              f'({reduction:.2f}% lower mean coding time)')

    return 0


if __name__ == '__main__':
    sys.exit(main())
