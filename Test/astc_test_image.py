#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# -----------------------------------------------------------------------------
# Copyright 2019-2026 Arm Limited
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
Run image quality and performance tests using reusable image selection and execution helpers.

Images are discovered under Test/Images using their directory and filename metadata. Results are
written as CSV files, optionally compared against references from an earlier encoder build.
'''

import argparse
from dataclasses import dataclass
import math
from pathlib import Path
import re
import sys

import testlib.encoder as te
from testlib.encoder import RunOptions
from testlib.resultset import ResultSet, ResultStatus
from testlib.runner import run_test_set
from testlib.selection import TEST_BLOCK_SIZES, TEST_QUALITIES, select_test_cases
from testlib.testset import TestSet


TEST_ENCODERS = ('none', 'native', 'universal', 'neon', 'sve_256', 'sve_128',
                 'sse2', 'sse4.1', 'avx2', 'avx256', 'avx512')
# AVX-512 variants stay opt-in so CI can probe host support before running them
ENCODER_GROUPS = {
    'all-aarch64': ('neon', 'sve_256', 'sve_128'),
    'all-x86': ('sse2', 'sse4.1', 'avx2'),
}
PROFILES = ('ldr', 'ldrs', 'hdr')
FORMATS = ('l', 'xy', 'rgb', 'rgba')


@dataclass(frozen=True)
class EncoderTarget:
    '''
    An encoder and its result file naming settings.

    Attributes:
        encoder: The encoder to run.
        name: The output CSV name prefix.
        out_dir: The directory receiving result CSVs.
        reference_name: The reference CSV name prefix, or None when rebuilding references.
    '''
    encoder: te.EncoderBase
    name: str
    out_dir: Path
    reference_name: str | None


def get_encoder_params(encoder_name: str, reference_name: str | None,
                       test_set_name: str) -> EncoderTarget:
    '''
    Get the encoder and image set parameters for a test run.

    Args:
        encoder_name: The encoder name.
        reference_name: The reference encoder name.
        test_set_name: The test image set.

    Return:
        The encoder and its output and reference naming settings.
    '''
    encoder: te.EncoderBase

    if encoder_name.startswith('ref'):
        _, version, simd = encoder_name.split('-')

        # Version-specific variants
        if re.match(r'[0-9]+\.', version):
            encoder = te.Encoder2xRel(version, simd)
            name = f'reference-{version}-{simd}'
            out_dir = f'Test/Images/{test_set_name}'

            return EncoderTarget(encoder=encoder, name=name, out_dir=Path(out_dir),
                                 reference_name=reference_name)

        # Latest main
        if version == 'main':
            encoder = te.Encoder2x(simd)
            name = f'reference-{version}-{simd}'
            out_dir = f'Test/Images/{test_set_name}'

            return EncoderTarget(encoder=encoder, name=name, out_dir=Path(out_dir),
                                 reference_name=reference_name)

        assert False, f'Encoder {encoder_name} not recognized'

    encoder = te.Encoder2x(encoder_name)
    name = f'develop-{encoder_name}'
    out_dir = f'TestOutput/{test_set_name}'

    if reference_name is not None:
        reference_name = reference_name.replace('ref', 'reference')

    return EncoderTarget(encoder=encoder, name=name, out_dir=Path(out_dir),
                         reference_name=reference_name)


def find_reference_encoders(test_dir: Path) -> list[str]:
    '''
    Find reference encoders with result CSVs in the image tree.

    Args:
        test_dir: The root of the test image tree.

    Return:
        The sorted, unique reference encoder names.
    '''
    quality_pattern = '|'.join(TEST_QUALITIES)
    result_pattern = re.compile(rf'astc_reference-(.+)_(?:{quality_pattern})_results\.csv')
    reference_encoders = set()

    for result_file in test_dir.rglob('astc_reference-*_results.csv'):
        match = result_pattern.fullmatch(result_file.name)

        if match and result_file.is_file():
            reference_encoders.add(f'ref-{match.group(1)}')

    return sorted(reference_encoders)


def add_selection_arguments(parser: argparse.ArgumentParser, test_dir: Path) -> None:
    '''
    Add encoder, image, and coverage selection arguments.

    Args:
        parser: The argument parser receiving these options.
        test_dir: The root directory containing image test sets.
    '''
    # Encoder aliases and stored references define the available comparison targets
    reference_encoders = ['rebuild', *find_reference_encoders(test_dir)]

    parser.add_argument('--encoder', dest='encoders', default='avx2',
                        help='test encoder variant or '
                             'ref-<version>-<variant>')

    parser.add_argument('--reference', default='ref-main-avx2',
                        choices=reference_encoders,
                        help='reference encoder variant')

    # Image filters apply before the test plan is constructed
    parser.add_argument('--color-profile', dest='profiles', default='all',
                        choices=[*PROFILES, 'all'], help='test color profile')

    parser.add_argument('--color-format', dest='formats', default='all',
                        choices=[*FORMATS, 'all'], help='test color format')

    parser.add_argument('--block-size', dest='block_sizes', action='append',
                        choices=[*TEST_BLOCK_SIZES, 'all'], help='test block size')

    test_set_choices = [x.name for x in test_dir.iterdir() if x.is_dir()]
    test_set_choices.sort()
    test_set_choices.append('all')
    parser.add_argument('--test-set', dest='test_sets', default='Small',
                        choices=test_set_choices, help='test image set')

    parser.add_argument('--test-image', dest='test_image', default=None,
                        help='select a specific test image from the test set')

    # Coverage selection spans all requested qualities
    quality_choices = list(TEST_QUALITIES)
    quality_choices.append('all')
    quality_choices.append('all+')
    parser.add_argument('--test-quality', dest='test_qual', action='append',
                        choices=quality_choices, help='test quality')

    parser.add_argument('--pairwise', action='store_true',
                        help='cover all profile-format/block/quality pairs '
                             'and images')

    parser.add_argument('--quality-update', action='store_true',
                        help='only update records with PSNR changes')


def add_run_arguments(parser: argparse.ArgumentParser) -> None:
    '''
    Add execution and timing precision arguments.

    Args:
        parser: The argument parser receiving these options.
    '''
    parser.add_argument('--repeats', dest='repeats', default=1,
                        type=int, help='maximum test repeats; stop early when '
                                       'the timing confidence interval meets '
                                       'the requested precision')

    parser.add_argument('--confidence', type=float, default=95.0,
                        help='timing confidence level in percent '
                             '(default: %(default)s)')

    parser.add_argument('--precision', type=float, default=1.0,
                        help='allowed timing confidence interval half-width '
                             'as a percent of the mean (default: %(default)s)')

    parser.add_argument('--keep-output', dest='keep_output', default=False,
                        action='store_true', help='keep image output')

    parser.add_argument('-j', dest='threads', default=None, type=int, help='encoder thread count')


def expand_all(value: str, choices: tuple[str, ...] | list[str]) -> list[str]:
    '''
    Expand an all selection into a fresh list of supported values.

    Args:
        value: The selected value, or all.
        choices: The values included by all.

    Return:
        The selected values without sharing mutable defaults.
    '''

    if value == 'all':
        return list(choices)

    return [value]


def normalize_selection(args: argparse.Namespace, test_dir: Path) -> None:
    '''
    Convert selection aliases into explicit lists without modifying shared constants.

    Args:
        args: The parsed arguments to normalize.
        test_dir: The root directory containing image test sets.
    '''
    # Copy shared selections so repeated parses cannot change later defaults
    args.encoders = list(ENCODER_GROUPS.get(args.encoders, (args.encoders,)))

    if (args.test_qual is not None) and ('all+' in args.test_qual):
        args.test_qual = list(TEST_QUALITIES)
    elif (args.test_qual is None) or ('all' in args.test_qual):
        args.test_qual = list(TEST_QUALITIES[:-2])

    if not args.block_sizes or ('all' in args.block_sizes):
        args.block_sizes = list(TEST_BLOCK_SIZES)

    test_sets = sorted(path.name for path in test_dir.iterdir() if path.is_dir())
    args.test_sets = expand_all(args.test_sets, test_sets)
    args.profiles = expand_all(args.profiles, PROFILES)
    args.formats = expand_all(args.formats, FORMATS)

    if args.reference == 'rebuild':
        args.reference = None


def parse_command_line(argv: list[str] | None = None) -> argparse.Namespace:
    '''
    Parse and normalize the image test command line.

    Args:
        argv: The command-line arguments, or None to read the process arguments.

    Return:
        The normalized selections and validated run options.
    '''
    parser = argparse.ArgumentParser()
    test_dir = Path(__file__).parent / 'Images'
    add_selection_arguments(parser, test_dir)
    add_run_arguments(parser)
    args = parser.parse_args(argv)

    # Keep validation errors in the units accepted by the command line
    if args.repeats < 1:
        parser.error('--repeats must be at least one')

    if not 0.0 < args.confidence < 100.0:
        parser.error('--confidence must be between 0 and 100, exclusive')

    if (not math.isfinite(args.precision)) or (args.precision <= 0.0):
        parser.error('--precision must be a finite positive percentage')

    # Allow new reference encoders before their result CSVs have been generated
    if args.encoders not in (*TEST_ENCODERS, *ENCODER_GROUPS):
        pattern = r'ref-(?:main|[0-9]+\.[^-]+)-(.+)'
        match = re.fullmatch(pattern, args.encoders)

        if (match is None) or (match.group(1) not in TEST_ENCODERS):
            parser.error(f'--encoder: invalid choice: {args.encoders!r}')

    # Shared execution code uses fractions rather than CLI percentages
    args.run_options = RunOptions(repeats=args.repeats, keep_output=args.keep_output,
                                  threads=args.threads, confidence=args.confidence / 100.0,
                                  precision=args.precision / 100.0)

    normalize_selection(args, test_dir)

    return args


def run_configuration(args: argparse.Namespace, quality: str, test_set_name: str,
                      encoder_name: str) -> ResultStatus:
    '''
    Run one encoder, test set, and quality preset and save its result CSV.

    Args:
        args: The normalized image selections and run settings.
        quality: The quality preset to run.
        test_set_name: The image test set to load.
        encoder_name: The encoder variant to run.

    Return:
        The worst comparison status, or NOT_RUN when no reference is supplied.
    '''
    target = get_encoder_params(encoder_name, args.reference, test_set_name)
    test_dir = Path('Test/Images') / test_set_name

    # Reference CSVs provide comparison data without influencing sampling
    reference = None

    if target.reference_name:
        ref_csv = test_dir / f'astc_{target.reference_name}_{quality}_results.csv'
        reference = ResultSet(test_set_name)
        reference.load_from_file(ref_csv)

    test_set = TestSet(test_set_name, test_dir, args.profiles, args.formats, args.test_image)

    # Select across all requested qualities before grouping to preserve pairwise coverage
    cases = select_test_cases(test_set, args.block_sizes, args.test_qual, args.pairwise)
    cases = [case for case in cases if case.quality == quality]

    results = run_test_set(target.encoder, test_set, quality, cases, args.run_options, reference)

    # Quality-only updates retain stored timing data for unchanged scores
    test_csv = target.out_dir / f'astc_{target.name}_{quality}_results.csv'

    if args.quality_update:
        results.save_quality_changes(test_csv)
    else:
        results.save_to_file(test_csv)

    if reference is not None:
        summary = results.get_results_summary()
        print(summary)

        return summary.get_worst_result()

    return ResultStatus.NOT_RUN


def main() -> int:
    '''
    Run the selected configurations and report the overall quality status.

    Return:
        One if a quality comparison fails, otherwise zero.
    '''
    args = parse_command_line()
    test_set_count = 0
    worst_result = ResultStatus.NOT_RUN

    # Preserve the quality, image-set, and encoder ordering used by existing reports
    for quality in args.test_qual:
        for test_set_name in args.test_sets:
            for encoder_name in args.encoders:
                status = run_configuration(args, quality, test_set_name, encoder_name)
                test_set_count += 1
                worst_result = max(status, worst_result)

        if (test_set_count > 1) and (worst_result != ResultStatus.NOT_RUN):
            print(f'OVERALL STATUS: {worst_result.name}')

    if worst_result == ResultStatus.FAIL:
        return 1

    return 0


if __name__ == '__main__':
    sys.exit(main())
