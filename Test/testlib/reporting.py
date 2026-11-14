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
Console reporting for image tests and profiling comparisons.
'''

import statistics

from .profile import CodecResult, QUALITY_METRICS, select_quality_metrics
from .resultset import Record


def format_solo_result(result: Record) -> str:
    '''
    Format a metrics string for a single (no compare) result.

    Args:
        result: The test result.

    Return:
        The metrics report.
    '''
    name = f'{result.block_size:>5} {result.name}'
    psnr = f'{result.psnr:2.3f} dB'
    ssim = f'{result.ssim:.3f}'
    total_time = f'{result.total_time:.3f} s'
    coding_time = f'{result.coding_time:.3f} s'
    coding_rate = f'{result.coding_rate:.3f} MT/s'

    fields = (
        f'{name:<33}',
        f'{psnr:>11}',
        f'{ssim:>5}',
        f'{total_time:>9}',
        f'{coding_time:>9}',
        f'{coding_rate:>12}',
    )

    return ' | '.join(fields)


def format_compare_result(result: Record) -> str:
    '''
    Format a metrics string for a comparison result.

    Args:
        result: The test result.

    Return:
        The metrics report.
    '''
    assert result.psnr_rel is not None
    assert result.ssim_rel is not None
    assert result.total_time_rel is not None
    assert result.coding_time_rel is not None

    name = f'{result.block_size:>5} {result.name}'
    psnr = f'{result.psnr:2.3f} dB ({result.psnr_rel:1.3f} dB)'
    ssim = f'{result.ssim:.3f} ({result.ssim_rel:+.3f})'
    total_time = f'{result.total_time:.3f} s ({result.total_time_rel:1.2f}x)'
    coding_time = f'{result.coding_time:.3f} s ({result.coding_time_rel:1.2f}x)'
    coding_rate = f'{result.coding_rate:.3f} MT/s'

    fields = (
        f'{name:<33}',
        f'{psnr:>22}',
        f'{ssim:>14}',
        f'{total_time:>15}',
        f'{coding_time:>15}',
        f'{coding_rate:>12}',
        result.status.name
    )

    return ' | '.join(fields)


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
        'mPSNR (RGB)': 'mSSIM (RGB)',
    }

    for metric in list(metrics):
        ssim = ssim_metrics.get(metric)

        if ssim and all(ssim in samples[0].quality for samples in (samples_a, samples_b)):
            metrics.append(ssim)

    print('Image quality (mean over measured runs; B - A):')

    for metric in metrics:
        if any(metric not in sample.quality for sample in samples_a + samples_b):
            raise ValueError(f'Codec did not report {metric} in every run')

        mean_a = statistics.mean(sample.quality[metric] for sample in samples_a)
        mean_b = statistics.mean(sample.quality[metric] for sample in samples_b)
        unit, direction = QUALITY_METRICS[metric]
        precision = 4

        if not unit:
            precision = 6

        print(f'  {metric}: A {mean_a:.{precision}f}, '
              f'B {mean_b:.{precision}f}, '
              f'change {mean_b - mean_a:+.{precision}f} {unit} '
              f'({direction} is better)')
