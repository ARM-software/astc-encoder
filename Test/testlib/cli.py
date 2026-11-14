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
Common command-line arguments for codec profiling tools.
'''

import argparse


def confidence_value(value: str) -> float:
    '''
    Parse a confidence level expressed as a fraction.

    Args:
        value: The command-line confidence value.

    Return:
        The confidence level in (0, 1).

    Raises:
        ArgumentTypeError: The confidence level is outside its valid range.
    '''
    result = float(value)

    if not 0.0 < result < 1.0:
        raise argparse.ArgumentTypeError('Must be between 0 and 1')

    return result


def add_codec_arguments(parser: argparse.ArgumentParser, is_pair: bool = False) -> None:
    '''
    Add binary paths followed by arguments passed through unchanged to the codec.

    Args:
        parser: The argument parser receiving these arguments.
        is_pair: Whether to accept two binaries for a comparison.
    '''
    if is_pair:
        parser.add_argument(
            'binary_a',
            help='the first astcenc binary')

        parser.add_argument(
            'binary_b',
            help='the second astcenc binary')
    else:
        parser.add_argument(
            'binary',
            help='the astcenc binary to profile')

    parser.add_argument(
        'arguments', nargs=argparse.REMAINDER,
        help='arguments passed through to the binary or binaries')


def parse_comparison_args(description: str | None, *, is_ttest: bool) -> argparse.Namespace:
    '''
    Parse paired timing arguments for a t-test or best-of-N comparison.

    Args:
        description: The script description shown in command-line help.
        is_ttest: Whether to select the paired t-test interface and minimum repeat count.

    Return:
        The binary paths, codec arguments, and validated sampling settings.
    '''
    parser = argparse.ArgumentParser(description=description)
    add_codec_arguments(parser, is_pair=True)

    # Preserve each tool's count flag and its statistical minimum
    if is_ttest:
        count_name = 'repeats'
        minimum = 2
        count_help = 'maximum measured pairs, at least 2 (default: %(default)s)'
        confidence_help = 'confidence level (default: %(default)s)'
    else:
        count_name = 'runs'
        minimum = 1
        count_help = 'number of measured pairs (default: %(default)s)'
        confidence_help = 'accepted for compatibility; otherwise ignored'

    parser.add_argument(
        f'--{count_name}', type=int, default=10,
        help=count_help)

    parser.add_argument(
        '--warmups', type=int, default=1,
        help='number of warm-up pairs (default: %(default)s)')

    parser.add_argument(
        '--confidence', type=confidence_value, default=0.95,
        help=confidence_help)

    args = parser.parse_args()

    if getattr(args, count_name) < minimum:
        parser.error(f'--{count_name} must be at least {minimum}')

    if args.warmups < 0:
        parser.error('--warmups cannot be negative')

    if not args.arguments:
        parser.error('no astcenc arguments specified')

    return args
