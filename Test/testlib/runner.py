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
Execute image test plans and report progress independently of command-line parsing.
'''

from .encoder import EncoderBase, RunOptions
from .reporting import format_compare_result, format_solo_result
from .resultset import Record, ResultSet
from .selection import ImageCase
from .testset import TestSet


def run_test_image(encoder: EncoderBase, case: ImageCase, options: RunOptions,
                   reference: ResultSet | None = None) -> Record:
    '''
    Compress one image and evaluate its metrics against an optional reference.

    Args:
        encoder: The encoder to test.
        case: The image, block size, and quality preset to test.
        options: The run settings and timing thresholds.
        reference: The reference records, or None for an independent run.

    Return:
        The measured record, including comparison metrics when a reference is supplied.
    '''
    result = encoder.run_test(case.image, case.block_size, f'-{case.quality}', options)

    # Convert measured metrics into the persistent CSV record format
    record = Record(block_size=case.block_size, name=case.image.file_name,
                    psnr=result.psnr, ssim=result.ssim, total_time=result.total_time,
                    coding_time=result.coding_time, coding_rate=result.coding_rate)

    if reference is not None:
        ref = reference.get_matching_record(record)
        record.compare_to(ref, case.image.is_3d)

    return record


def run_test_set(encoder: EncoderBase, test_set: TestSet, quality: str, cases: list[ImageCase],
                 options: RunOptions, reference: ResultSet | None = None) -> ResultSet:
    '''
    Run a test plan for one quality preset and print progress and result rows.

    Args:
        encoder: The encoder to test.
        test_set: The source image set.
        quality: The quality preset shared by the selected cases.
        cases: The selected cases, all using the same quality preset.
        options: The run settings and timing thresholds.
        reference: The reference records, or None for an independent run.

    Return:
        The measured records in test plan order.
    '''
    results = ResultSet(test_set.name)

    title = f'Test Set: {test_set.name} / Encoder: {encoder.name} -{quality}'
    print(title)
    print('=' * len(title))

    for index, case in enumerate(cases, start=1):
        progress = f'{index}/{len(cases)}'
        print(f'Running {progress} {case.block_size} {case.image.file_name} ... ',
              end='', flush=True)

        record = run_test_image(encoder, case, options, reference)
        results.add_record(record)

        # Rendering consumes the comparison already stored on the record
        if reference is None:
            report = format_solo_result(record)
        else:
            report = format_compare_result(record)

        print(f'\r[{index:>3}] {report} | {encoder.last_run_count} runs')

    return results
