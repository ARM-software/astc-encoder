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
The image test runner is used for image quality and performance testing.

It is designed to process directories of arbitrary test images, using the
directory structure and path naming conventions to self-describe how each image
is to be compressed. Some built-in test sets are provided in the ./Test/Images
directory, and others can be downloaded by running the astc_test_image_dl
script.

Attributes:
    RESULT_THRESHOLD_WARN: The result threshold (dB) for getting a WARN.
    RESULT_THRESHOLD_FAIL: The result threshold (dB) for getting a FAIL.
    TEST_BLOCK_SIZES: The block sizes we can test. This is a subset of the
        block sizes supported by ASTC, simply to keep test run times
        manageable.
'''

import argparse
import os
from pathlib import Path
import sys
from typing import Optional

import testlib.encoder as te
import testlib.testset as tts
import testlib.resultset as trs

# Require bit exact with reference scores
RESULT_THRESHOLD_WARN = -0.00
RESULT_THRESHOLD_FAIL = -0.00
RESULT_THRESHOLD_3D_FAIL = -0.00

TEST_BLOCK_SIZES = [
    '4x4', '5x5', '6x6', '8x8', '12x12', '3x3x3', '6x6x6'
]

TEST_QUALITIES = [
    'fastest', 'fast', 'medium', 'thorough', 'verythorough', 'exhaustive'
]


def is_3d(block_size: str) -> bool:
    '''
    Is the given block size a 3D block type?

    Args:
        block_size: The block size.

    Return:
        True if the block size if 3D, False if 2D.
    '''
    return block_size.count('x') == 2


def count_test_set(test_set: tts.TestSet, block_sizes: list[str]) -> int:
    '''
    Count the number of test executions needed for a test set.

    Args:
        test_set: The test set to run.
        block_sizes: The block sizes to run.

    Return:
        The number of test executions needed.
    '''
    count = 0
    for block_size in block_sizes:
        for image in test_set.tests:
            # 3D block sizes require 3D images
            if is_3d(block_size) != image.is_3d:
                continue

            count += 1

    return count


def determine_result(image: tts.TestImage, reference: trs.Record,
                     result: trs.Record) -> trs.ResultStatus:
    '''
    Determine a test result against a reference and thresholds.

    Args:
        image: The image being compressed.
        reference: The reference result to compare against.
        result: The test result.

    Return:
        The result code.
    '''
    psnr_diff = result.psnr - reference.psnr

    if (psnr_diff < RESULT_THRESHOLD_FAIL) and (not image.is_3d):
        return trs.ResultStatus.FAIL

    if (psnr_diff < RESULT_THRESHOLD_3D_FAIL) and image.is_3d:
        return trs.ResultStatus.FAIL

    if psnr_diff < RESULT_THRESHOLD_WARN:
        return trs.ResultStatus.WARN

    return trs.ResultStatus.PASS


def format_solo_result(image: tts.TestImage, result: trs.Record) -> str:
    '''
    Format a metrics string for a single (no compare) result.

    Args:
        image: The image being tested.
        result: The test result.

    Return:
        The metrics report.
    '''
    del image  # Unused

    name = f'{result.block_size:>5} {result.name}'
    psnr = f'{result.psnr:2.3f} dB'
    total_time = f'{result.total_time:.3f} s'
    coding_time = f'{result.coding_time:.3f} s'
    coding_rate = f'{result.coding_rate:.3f} MT/s'

    fields = (
        f'{name:<32}',
        f'{psnr:>11}',
        f'{total_time:>9}',
        f'{coding_time:>9}',
        f'{coding_rate:>11}',
    )

    return ' | '.join(fields)


def format_compare_result(image: tts.TestImage, reference: trs.Record,
                          result: trs.Record) -> str:
    '''
    Format a metrics string for a comparison result.

    Args:
        image: The image being tested.
        reference: The reference result to compare against.
        result: The test result.

    Return:
        The metrics report.
    '''
    psnr_diff = result.psnr - reference.psnr

    try:
        total_time_rel = reference.total_time / result.total_time
    except ZeroDivisionError:
        total_time_rel = float('NaN')

    try:
        coding_time_rel = reference.coding_time / result.coding_time
    except ZeroDivisionError:
        coding_time_rel = float('NaN')

    name = f'{result.block_size:>5} {result.name}'
    psnr = f'{result.psnr:2.3f} dB ({psnr_diff:1.3f} dB)'
    total_time = f'{result.total_time:.3f} s ({total_time_rel:1.2f}x)'
    coding_time = f'{result.coding_time:.3f} s ({coding_time_rel:1.2f}x)'
    coding_rate = f'{result.coding_rate:.3f} MT/s'
    result_status = determine_result(image, reference, result)

    fields = (
        f'{name:<32}',
        f'{psnr:>22}',
        f'{total_time:>15}',
        f'{coding_time:>15}',
        f'{coding_rate:>11}',
        f'{coding_rate:>11}',
        result_status.name
    )

    return ' | '.join(fields)


def run_test_set(encoder, reference, test_set, quality, block_sizes, repeats,
                 keep_output, threads):
    '''
    Execute all tests in the test set.

    Args:
        encoder (EncoderBase): The encoder to use.
        reference (ResultSet): The test reference results.
        test_set (TestSet): The test set.
        quality (str): The quality level to execute the test against.
        block_sizes (list(str)): The block sizes to execute each test against.
        repeats (int): The number of test repeats to run for each image test.
        keep_output (bool): Should the test preserve output images? This is
            only a hint and discarding output may be ignored if the encoder
            version used can't do it natively.
        threads (int or None): The thread count to use.

    Return:
        ResultSet: The test results.
    '''
    result_set = trs.ResultSet(test_set.name)

    current_test_index = 0
    max_test_index = count_test_set(test_set, block_sizes)

    title = f'Test Set: {test_set.name} / Encoder: {encoder.name} -{quality}'
    print(title)
    print('=' * len(title))

    for block_size in block_sizes:
        for image in test_set.tests:
            # 3D block sizes require 3D images
            if is_3d(block_size) != image.is_3d:
                continue

            current_test_index += 1

            progress = f'{current_test_index}/{max_test_index}'
            msg = f'Running {progress} {block_size} {image.file_name} ... '
            print(msg, end='', flush=True)

            res = encoder.run_test(
                image, block_size, f'-{quality}', repeats,
                keep_output, threads)

            res = trs.Record(block_size, image.file_name, *res)
            result_set.add_record(res)

            if reference:
                result_ref = reference.get_matching_record(res)
                res.set_status(determine_result(image, result_ref, res))

                try:
                    res.total_time_rel = result_ref.total_time / res.total_time
                except ZeroDivisionError:
                    res.total_time_rel = float('NaN')

                try:
                    res.coding_time_rel = \
                        result_ref.coding_time / res.coding_time
                except ZeroDivisionError:
                    res.coding_time_rel = float('NaN')

                res.psnr_rel = res.psnr - result_ref.psnr

                res = format_compare_result(image, result_ref, res)
            else:
                res = format_solo_result(image, res)

            print(f'\r[{current_test_index:>3}] {res}')

    return result_set


def get_encoder_params(encoder_name: str, reference_name: str,
                       test_set_name: str) \
                       -> tuple[te.EncoderBase, str, str, str]:
    '''
    The the encoder and image set parameters for a test run.

    Args:
        encoder_name: The encoder name.
        reference_name: The reference encoder name.
        test_set_name: The test image set.

    Return:
        The test parameters for the requested encoder and test set. An instance
        of the encoder wrapper class, the output data name, the output result
        directory, and the reference to use.
    '''
    used_reference_name: Optional[str] = None
    encoder: te.EncoderBase

    if encoder_name.startswith('ref'):
        _, version, simd = encoder_name.split('-')

        # Version-specific variants
        compatible_prefixes = ['5.']
        if any(True for x in compatible_prefixes if version.startswith(x)):
            encoder = te.Encoder2xRel(version, simd)
            name = f'reference-{version}-{simd}'
            out_dir = f'Test/Images/{test_set_name}'

            return (encoder, name, out_dir, reference_name)

        # Latest main
        if version == 'main':
            encoder = te.Encoder2x(simd)
            name = f'reference-{version}-{simd}'
            out_dir = f'Test/Images/{test_set_name}'
            return (encoder, name, out_dir, reference_name)

        assert False, f'Encoder {encoder_name} not recognized'

    encoder = te.Encoder2x(encoder_name)
    name = f'develop-{encoder_name}'
    out_dir = f'TestOutput/{test_set_name}'
    used_reference_name = reference_name.replace('ref', 'reference')
    return (encoder, name, out_dir, used_reference_name)


def parse_command_line():
    '''
    Parse the command line.

    Return:
        Namespace: The parsed command line container.
    '''
    parser = argparse.ArgumentParser()

    # All reference encoders
    reference_encoders = [
        'ref-5.0-neon',
        'ref-5.0-sse2', 'ref-5.0-sse4.1', 'ref-5.0-avx2',

        'ref-main-neon', 'ref-main-sve_256', 'ref-main-sve_128',
        'ref-main-sse2', 'ref-main-sse4.1', 'ref-main-avx2'
    ]

    # All test encoders
    test_encoders = [
        'none', 'native', 'universal',
        'neon', 'sve_256', 'sve_128',
        'sse2', 'sse4.1', 'avx2'
    ]

    test_encoders_arm64 = [
        'neon', 'sve_256', 'sve_128'
    ]

    test_encoders_x86 = [
        'sse2', 'sse4.1', 'avx2'
    ]

    coders = reference_encoders + test_encoders + ['all-aarch64', 'all-x86']

    parser.add_argument('--encoder', dest='encoders',
                        default='avx2',
                        choices=coders, help='test encoder variant')

    parser.add_argument('--reference', default='ref-main-avx2',
                        choices=reference_encoders,
                        help='reference encoder variant')

    astc_profile = ['ldr', 'ldrs', 'hdr', 'all']
    parser.add_argument('--color-profile', dest='profiles', default='all',
                        choices=astc_profile, help='test color profile')

    color_format = ['l', 'xy', 'rgb', 'rgba', 'all']
    parser.add_argument('--color-format', dest='formats', default='all',
                        choices=color_format, help='test color format')

    choices = list(TEST_BLOCK_SIZES) + ['all']
    parser.add_argument('--block-size', dest='block_sizes',
                        action='append', choices=choices,
                        help='test block size')

    test_dir = os.path.dirname(__file__)
    test_dir = os.path.join(test_dir, 'Images')
    test_sets = []
    for path in os.listdir(test_dir):
        full_path = os.path.join(test_dir, path)
        if os.path.isdir(full_path):
            test_sets.append(path)
    test_sets.append('all')

    parser.add_argument('--test-set', dest='test_sets', default='Small',
                        choices=test_sets, help='test image test set')

    parser.add_argument('--test-image', dest='test_image', default=None,
                        help='select a specific test image from the test set')

    choices = list(TEST_QUALITIES) + ['all', 'all+']
    parser.add_argument('--test-quality', dest='test_qual', default='thorough',
                        choices=choices, help='select a specific test quality')

    parser.add_argument('--repeats', dest='repeats', default=1,
                        type=int, help='test iteration count')

    parser.add_argument('--keep-output', dest='keep_output', default=False,
                        action='store_true', help='keep image output')

    parser.add_argument('-j', dest='threads', default=None,
                        type=int, help='thread count')

    args = parser.parse_args()

    # Turn things into canonical format lists
    if args.encoders == 'all-aarch64':
        args.encoders = test_encoders_arm64
    elif args.encoders == 'all-x86':
        args.encoders = test_encoders_x86
    else:
        args.encoders = [args.encoders]

    if args.test_qual == 'all+':
        args.test_qual = TEST_QUALITIES
    elif args.test_qual == 'all':
        args.test_qual = TEST_QUALITIES
        args.test_qual.remove('verythorough')
        args.test_qual.remove('exhaustive')
    else:
        args.test_qual = [args.test_qual]

    if not args.block_sizes or ('all' in args.block_sizes):
        args.block_sizes = TEST_BLOCK_SIZES

    args.test_sets = test_sets[:-1] if args.test_sets == 'all' \
        else [args.test_sets]

    args.profiles = astc_profile[:-1] if args.profiles == 'all' \
        else [args.profiles]

    args.formats = color_format[:-1] if args.formats == 'all' \
        else [args.formats]

    return args


def main() -> int:
    '''
    The main function.

    Return:
        int: The process return code.
    '''
    # Parse command lines
    args = parse_command_line()

    test_set_count = 0
    worst_result = trs.ResultStatus.NOT_RUN

    for quality in args.test_qual:
        for test_set_name in args.test_sets:
            for encoder_name in args.encoders:
                (encoder, name, out_dir, reference_name) = \
                    get_encoder_params(
                        encoder_name, args.reference, test_set_name)

                test_dir = f'Test/Images/{test_set_name}'
                test_csv = f'{out_dir}/astc_{name}_{quality}_results.csv'

                reference = None
                if reference_name:
                    ref_csv = f'astc_{reference_name}_{quality}_results.csv'
                    ref_csv = f'{test_dir}/{ref_csv}'
                    reference = trs.ResultSet(test_set_name)
                    reference.load_from_file(Path(ref_csv))

                test_set_count += 1
                test_set = tts.TestSet(
                    test_set_name, Path(test_dir),
                    args.profiles, args.formats, args.test_image)

                result_set = run_test_set(
                    encoder, reference, test_set, quality,
                    args.block_sizes, args.repeats,
                    args.keep_output, args.threads)

                result_set.save_to_file(Path(test_csv))

                if reference_name:
                    summary = result_set.get_results_summary()
                    this_result = summary.get_worst_result()
                    worst_result = max(this_result, worst_result)
                    print(summary)

        if (test_set_count > 1) and (worst_result != trs.ResultStatus.NOT_RUN):
            print(f'OVERALL STATUS: {worst_result.name}')

    if worst_result == trs.ResultStatus.FAIL:
        return 1

    return 0


if __name__ == '__main__':
    sys.exit(main())
