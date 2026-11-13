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
Benchmark three astcenc binaries and plot RGB, normal-map, and HDR throughput.

Requires Python 3.10 or newer and the packages in requirements.txt. Shared profiling helpers live
in testlib beside this script. The report layout does not depend on an external template.

Example, from the project root:

  python3 Test/astc_profile_report.py ./release-avx2 ./main-avx2 ./main-avx512

Default images are the KodakMnt RGB and normal montages and HDRIHaven Canarywharf, under Images
beside this script. Use --rgb-image, --normal-image, and --hdr-image when copying the script
elsewhere with testlib. The HDR input is a Radiance .hdr image. Binaries can be paths or commands
on PATH. Use matched build options and give the script exclusive use of the machine during
measurement.

Runs 4x4, 6x6, and 8x8 at -thorough, -medium, and -fast for all three images. Each case has one
warmup and at most ten measured rounds, with 12 codec threads. After at least two measured rounds,
stop a case when every binary's 95% timing confidence interval has a half-width within 1% of its
mean. Repeated checks use nominal confidence, not a sequentially adjusted guarantee. The three
binaries run serially in rotating order. All three comparisons reuse the same samples, keeping each
binary's mean consistent between tables. HDR compression uses -perceptual. HDR metrics use a fixed
-3 to +1 f-stop window across binary versions, matching the current compressor's automatic exposure
for Canarywharf.

Writes report.md, report.json, six {rgb,nml,hdr}-{perf_abs,perf_rel}.png files, and three
{rgb,nml,hdr}-qual.png files in the output directory. Each quality
plot shows signed B-A changes for baseline versus new AVX2 in two metric panels:
PSNR and SSIM for RGB, PSNR and mean angular error for normals, and mPSNR and mSSIM for HDR. Metrics
missing from older binaries are marked unavailable. The default output directory is
x_profiles/report-YYMMDD-HH.MM. Existing files with these names are overwritten. Tables report
coding seconds; plots report MTex/s or throughput relative to A = 100%. Throughput is image texels
divided by the unrounded mean coding time. Warmups are excluded from means. There is no checkpoint
or resume support.
'''

import argparse
from dataclasses import asdict
from datetime import datetime
import json
from pathlib import Path
import platform
import re
import statistics
import sys
import tempfile
from typing import Any

from PIL import Image

from testlib.profile import binary_metadata, CodecResult, collect_precision_samples, fingerprint
from testlib.profile import PROFILE_ERRORS, report_error
from testlib.profile import run_codec as execute_codec


BINARY_NAMES = {
    'release-avx2': '5.7.1 AVX2',
    'main-avx2': 'main AVX2',
    'main-avx512': 'main AVX512',
}
SETS = (
    ('release-avx2', 'main-avx2'),
    ('release-avx2', 'main-avx512'),
    ('main-avx2', 'main-avx512'),
)
SCENARIOS = [(block, preset) for block in ('4x4', '6x6', '8x8')
             for preset in ('-thorough', '-medium', '-fast')]
QUALITY = ('PSNR (LDR-RGB)', 'SSIM (LDR-RGB)', 'Mean Angular Error',
           'Worst Angular Error', 'mPSNR (RGB)', 'mSSIM (RGB)')
HDR_FSTOPS = (-3, 1)
CONFIDENCE = 0.95
PRECISION = 0.01
QUALITY_PLOTS = {
    'RGB': (('PSNR (LDR-RGB)', 'PSNR', 'dB', 'higher'),
            ('SSIM (LDR-RGB)', 'SSIM', '', 'higher')),
    'Normal': (('PSNR (LDR-RGB)', 'PSNR', 'dB', 'higher'),
               ('Mean Angular Error', 'Mean angular error', 'degrees', 'lower')),
    'HDR': (('mPSNR (RGB)', 'mPSNR', 'dB', 'higher'),
            ('mSSIM (RGB)', 'mSSIM', '', 'higher')),
}
# Use the older report's light A and dark B shades for each block size
BLOCK_COLORS = {
    '4x4': ('#a8c7e3', '#35689a'),  # Light and dark blue
    '6x6': ('#f1bb99', '#ce6a2c'),  # Light and dark orange
    '8x8': ('#9dd5c6', '#2b8a74'),  # Light and dark teal
}
BACKGROUND_COLOR = 'white'  # White figure background
ANNOTATION_COLOR = '#64748b'  # Slate gray for reference lines and missing metrics


MEASUREMENT_NOTES = '''## Measurement notes

* Started: {started}; finished: {finished}.
* {rounds}, capped at {runs}, with {warmups} warmup rounds and {threads} threads.
* Stop after at least two measured rounds when every binary's mean timing CI meets the precision
  target ({thresholds}); otherwise continue to the cap. Confidence is nominal under repeated checks.
* Invocations run serially in rotating order. All three sets reuse the same per-binary samples.
* Times are arithmetic means of reported coding times. B-A (%) = 100 x (B/A - 1); negative means B
  takes less time.
* Quality plots show baseline versus new AVX2 signed B-A changes in the metric's own units. Higher
  PSNR/SSIM is better; lower angular error is better. Metrics not reported by a comparison binary
  are unavailable.
* HDR compression uses -perceptual. Quality metrics use a fixed -3 to +1 f-stop window for every
  binary, matching Canarywharf's current auto-exposure.
* Host: {host}.
* Binary versions, SHA-256 hashes, image dimensions, quality scores, and raw samples:
  [report.json](report.json).
'''


def parse_command_line() -> argparse.Namespace:
    '''
    Parse executable paths, image paths, and measurement settings.

    Return:
        The validated command-line arguments.
    '''
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('release_avx2', help='5.7.1 AVX2 binary')
    parser.add_argument('main_avx2', help='main AVX2 binary')
    parser.add_argument('main_avx512', help='main AVX512 binary')
    images = Path(__file__).resolve().parent / 'Images'
    kodak = images / 'KodakMnt'
    parser.add_argument('--rgb-image', type=Path, default=kodak / 'LDR-RGB/ldr-rgb-montage.png')
    parser.add_argument('--normal-image', type=Path, default=kodak / 'LDR-XY/ldr-xy-montage.png')
    hdr = images / 'HDRIHaven/HDR-RGB/hdr-rgb-canarywharf.hdr'
    parser.add_argument('--hdr-image', type=Path, default=hdr,
                        help='Radiance HDR image (default: Canarywharf)')
    stamp = datetime.now().strftime('%y%m%d-%H.%M')
    out_dir = Path('x_profiles') / f'report-{stamp}'
    parser.add_argument('--output', type=Path, default=out_dir,
                        help='output directory (default: %(default)s)')
    parser.add_argument('--runs', type=int, default=10,
                        help='maximum rounds per case, 2 to 10 (default: %(default)s)')
    parser.add_argument('--warmups', type=int, default=1,
                        help='warmup rounds per case (default: %(default)s)')
    parser.add_argument('--threads', type=int, default=12,
                        help='codec threads (default: %(default)s)')
    args = parser.parse_args()

    if args.runs < 2 or args.runs > 10:
        parser.error('--runs must be between 2 and 10')

    if args.threads < 1 or args.warmups < 0:
        parser.error('threads must be positive; warmups cannot be negative')

    return args


def image_size(path: Path) -> tuple[int, int]:
    '''
    Read image dimensions, including Radiance HDR headers unsupported by Pillow.

    Args:
        path: The input image path.

    Return:
        The image width and height.

    Raises:
        OSError: The image cannot be read.
        ValueError: The Radiance header or dimensions are invalid.
    '''

    if path.suffix.lower() != '.hdr':
        with Image.open(path) as image:
            return image.size

    with path.open('rb') as stream:
        if stream.readline().strip() not in (b'#?RADIANCE', b'#?RGBE'):
            raise ValueError(f'Invalid Radiance header: {path}')

        while True:
            line = stream.readline()

            if not line:
                raise ValueError(f'Missing Radiance dimensions: {path}')

            if not line.strip():
                break

        line = stream.readline()
        match = re.fullmatch(rb'[+-]Y\s+(\d+)\s+[+-]X\s+(\d+)\s*', line)

        if not match:
            raise ValueError(f'Invalid Radiance dimensions: {path}')

        height, width = map(int, match.groups())

        if width <= 0 or height <= 0:
            raise ValueError(f'Invalid Radiance dimensions: {path}')

    return width, height


def run_codec(binary: str, arguments: list[str]) -> CodecResult:
    '''
    Run the shared codec reader and require the metrics used by this report.

    Args:
        binary: The executable path.
        arguments: The codec command-line arguments.

    Return:
        The coding time and available report quality metrics.

    Raises:
        ValueError: A required report metric is missing.
    '''
    result = execute_codec(binary, arguments)
    required = ['PSNR (LDR-RGB)']

    if arguments[0] in ('-th', '-tH'):
        required = ['mPSNR (RGB)']

    if '-normal' in arguments:
        required.append('Mean Angular Error')

    # Required metrics fail immediately; optional metrics remain unavailable on older binaries
    for name in required:
        if name not in result.quality:
            raise ValueError(f'{binary} did not report {name}')

    quality = {name: value for name, value in result.quality.items() if name in QUALITY}

    return CodecResult(time=result.time, quality=quality)


def measure_case(binaries: dict[str, Any], args: list[str], runs: int, warmups: int,
                 confidence: float, precision: float) -> dict[str, Any]:
    '''
    Collect one scenario using shared serial precision sampling.

    Args:
        binaries: Binary paths and identities keyed by role.
        args: The codec arguments shared by all binaries.
        runs: The maximum number of measured rounds.
        warmups: The number of excluded warmup rounds.
        confidence: The confidence level as a fraction in (0, 1).
        precision: The interval half-width target as a fraction of the mean.

    Return:
        The samples, intervals, and stopping reason in report JSON format.
    '''
    paths = {name: binaries[name]['path'] for name in BINARY_NAMES}
    samples = collect_precision_samples(paths, args, runs, warmups, confidence, precision,
                                        run_codec)

    return asdict(samples)


def codec_args(image: str, path: str, output: Path, block: str,
               preset: str, threads: int) -> list[str]:
    '''
    Build codec arguments for one image and compression setting.

    Args:
        image: The image type.
        path: The source image path.
        output: The temporary output path without its extension.
        block: The ASTC block size.
        preset: The compression quality preset.
        threads: The codec thread count.

    Return:
        The codec arguments, including fixed HDR settings when applicable.
    '''
    mode = '-tl'
    suffix = '.png'

    if image == 'HDR':
        mode = '-th'
        suffix = '.exr'

    args = [mode, path, str(output.with_suffix(suffix)), block, preset, '-j', str(threads)]

    if image == 'Normal':
        args.append('-normal')

    if image == 'HDR':
        # Freeze Canarywharf's current auto-exposure across binary versions
        low, high = HDR_FSTOPS
        args += ['-perceptual', '-mpsnr', str(low), str(high)]

    return args


def collect_samples(binaries: dict[str, Any], images: dict[str, Any],
                    runs: int, warmups: int, threads: int,
                    confidence: float = CONFIDENCE,
                    precision: float = PRECISION) -> dict[str, Any]:
    '''
    Measure every case with one codec process at a time.

    Args:
        binaries: Binary paths and identities, keyed by role.
        images: Image paths, dimensions, and identities, keyed by type.
        runs: The maximum measured rounds per case, capped at ten.
        warmups: The number of warmup rounds per case.
        threads: The codec thread count.
        confidence: The confidence level as a fraction between zero and one.
        precision: The allowed interval half-width as a fraction of the mean.

    Return:
        The arguments, warmups, and raw samples for every scenario.
    '''
    runs = min(runs, 10)
    cases: dict[str, Any] = {}
    total = len(images) * len(SCENARIOS)
    with tempfile.TemporaryDirectory(prefix='astc-profile-') as tmp:
        for image, meta in images.items():
            for block, preset in SCENARIOS:
                key = f'{image}-{block}-{preset.removeprefix("-")}'
                args = codec_args(image, meta['path'], Path(tmp) / 'output', block, preset, threads)
                print(f'Case {len(cases) + 1}/{total}: {key}', flush=True)
                cases[key] = measure_case(binaries, args, runs, warmups, confidence, precision)

    return cases


def mean_time(case: dict[str, Any], binary: str) -> float:
    '''
    Compute the mean coding time, excluding warmups.

    Args:
        case: The measured scenario.
        binary: The binary role to summarize.

    Return:
        The arithmetic mean coding time in seconds.
    '''

    return statistics.mean(r['results'][binary]['time'] for r in case['rounds'])


def mean_quality(case: dict[str, Any], binary: str, metric: str) -> float | None:
    '''
    Compute a mean quality score when every measured round reports the metric.

    Args:
        case: The measured scenario.
        binary: The binary role to summarize.
        metric: The quality score to summarize.

    Return:
        The mean score, or None if the metric is missing from any round.
    '''
    scores = [r['results'][binary]['quality'].get(metric) for r in case['rounds']]

    if any(score is None for score in scores):
        return None

    return statistics.mean(scores)


def round_count_text(data: dict[str, Any], image: str | None = None) -> str:
    '''
    Describe the actual measured-round counts across report scenarios.

    Args:
        data: The recorded scenarios and their measured rounds.
        image: Optional image type to restrict the count to one plot.

    Return:
        A measured-round count or range, excluding warmups.
    '''
    counts = []

    for key, case in data['cases'].items():
        if image is not None and not key.startswith(image + '-'):
            continue

        counts.append(len(case['rounds']))

    low, high = min(counts), max(counts)

    if low == high:
        return f'{low} measured rounds per case'

    return f'{low}–{high} measured rounds per case'


def sampling_threshold_text(data: dict[str, Any]) -> str:
    '''
    Describe the sampling thresholds recorded for a report.

    Args:
        data: The recorded confidence and relative precision thresholds.

    Return:
        The confidence and precision percentages, or N/A for missing metadata.
    '''
    ci_text = 'CI: N/A'
    precision_text = 'Precision: N/A'

    confidence = data.get('confidence')
    precision = data.get('precision')

    if confidence is not None:
        percent = confidence * 100.0
        ci_text = f'CI: {percent:g}%'

    if precision is not None:
        percent = precision * 100.0
        precision_text = f'Precision: ±{percent:g}%'

    return f'{ci_text} | {precision_text}'


def timing_rows(data: dict[str, Any], image: str, a: str, b: str) -> list[str]:
    '''
    Format the timing comparison rows for one image and binary pair.

    Args:
        data: The measured scenarios.
        image: The image type.
        a: The baseline binary role.
        b: The candidate binary role.

    Return:
        The Markdown rows in scenario order.
    '''
    lines = []

    for block, preset in SCENARIOS:
        key = f'{image}-{block}-{preset.removeprefix("-")}'
        case = data['cases'][key]
        time_a = mean_time(case, a)
        time_b = mean_time(case, b)
        change = 100 * (time_b / time_a - 1)
        count = len(case['rounds'])
        lines.append(f'| {block:<10} | {preset:<9} | {time_a:11.2f} | '
                     f'{time_b:11.2f} | {change:7.2f} | {count:4} |')

    return lines


def write_report(out_dir: Path, data: dict[str, Any]) -> None:
    '''
    Write nine comparison tables and their measurement provenance.

    Args:
        out_dir: The directory containing the report artifacts.
        data: The settings, identities, and measurements.
    '''
    lines = []

    for index, (a, b) in enumerate(SETS, start=1):
        lines += [f'# Set {index}', '',
                  f'* Revision A: {BINARY_NAMES[a]}',
                  f'* Revision B: {BINARY_NAMES[b]}',
                  f'* Threads: {data["threads"]}', '', '']

        for image, metadata in data['images'].items():
            lines += [f'## Test: {metadata["title"]}', '',
                      '| Block size | Quality   | A (seconds) | B (seconds) | B-A (%) | Runs |',
                      '| ---------- | --------- | ----------: | ----------: | ------: | ---: |']
            lines += timing_rows(data, image, a, b)
            lines.append('')

    notes = MEASUREMENT_NOTES.format(
        started=data['started'], finished=data['finished'], rounds=round_count_text(data),
        runs=data['runs'], warmups=data['warmups'], threads=data['threads'],
        thresholds=sampling_threshold_text(data), host=data['host'])
    lines.append(notes)

    output = out_dir / 'report.md'
    output.write_text('\n'.join(lines), encoding='utf-8')


def add_comparison_legend(ax: Any, a: str, b: str, is_change: bool = False) -> None:
    '''
    Add the shared block-size and binary-role legend to a comparison panel.

    Args:
        ax: The plot axes receiving the legend.
        a: The baseline binary role.
        b: The candidate binary role.
        is_change: Show only candidate block colors for signed-change bars.
    '''
    # Keep Matplotlib optional until plots are requested
    from matplotlib.patches import Patch  # pylint: disable=import-outside-toplevel

    roles = [(0, 'A'), (1, 'B')]

    if is_change:
        roles = [(1, 'B')]

    handles = [Patch(facecolor=pair[i], label=f'{block} {role}')
               for block, pair in BLOCK_COLORS.items()
               for i, role in roles]

    ax.legend(handles=handles, loc='upper right', frameon=False,
              ncol=3, bbox_to_anchor=(1.0, 1.24), fontsize=9,
              title=f'A: {BINARY_NAMES[a]} | B: {BINARY_NAMES[b]}')


def new_figure(size: tuple[int, int], panels: int) -> tuple[Any, Any]:
    '''
    Create a figure with vertically stacked panels and a headless rendering canvas.

    Args:
        size: The figure width and height in inches.
        panels: The number of panels.

    Return:
        The figure and its axes.
    '''
    # pylint: disable=import-outside-toplevel

    # Keep Matplotlib optional until plots are requested
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    figure = Figure(figsize=size, facecolor=BACKGROUND_COLOR)
    FigureCanvasAgg(figure)

    return figure, figure.subplots(panels, 1)


def perf_rates(data: dict[str, Any], image: str, pair: tuple[str, str],
               is_relative: bool) -> list[list[float]]:
    '''
    Compute the baseline and candidate throughput values for a panel.

    Args:
        data: The settings, dimensions, and measurements.
        image: The image type.
        pair: The baseline and candidate binary roles.
        is_relative: Report throughput relative to the baseline.

    Return:
        The baseline and candidate throughput series.
    '''
    a, b = pair
    pixels = data['images'][image]['width'] * data['images'][image]['height']
    rates: list[list[float]] = [[], []]

    for block, preset in SCENARIOS:
        key = f'{image}-{block}-{preset.removeprefix("-")}'
        case = data['cases'][key]
        time_a = mean_time(case, a)
        time_b = mean_time(case, b)

        if is_relative:
            rates[0].append(100.0)
            rates[1].append(100.0 * time_a / time_b)
        else:
            rates[0].append(pixels / time_a / 1e6)
            rates[1].append(pixels / time_b / 1e6)

    return rates


def plot_perf_panel(ax: Any, data: dict[str, Any], image: str, pair: tuple[str, str],
                    is_relative: bool) -> None:
    '''
    Draw one absolute or relative performance comparison.

    Args:
        ax: The plot axes receiving the panel.
        data: The settings, dimensions, and measurements.
        image: The image type.
        pair: The baseline and candidate binary roles.
        is_relative: Report throughput relative to the baseline.
    '''
    a, b = pair
    rates = perf_rates(data, image, pair, is_relative)

    for offset, values in enumerate(rates):
        positions = [i + (offset - 0.5) * 0.36 for i in range(9)]
        shades = [BLOCK_COLORS[block][offset] for block, _ in SCENARIOS]
        bars = ax.bar(positions, values, width=0.36, color=shades)

        if is_relative and offset == 0:
            continue

        fmt = '%.2f'

        if is_relative:
            fmt = '%.2f%%'

        ax.bar_label(bars, fmt=fmt, padding=2, fontsize=9)

    ax.set_title(f'Set {SETS.index(pair) + 1}: {BINARY_NAMES[a]} vs {BINARY_NAMES[b]}',
                 loc='left', pad=35)
    unit = 'Throughput (MTex/s)'

    if is_relative:
        unit = 'Throughput (% of A)'

    ax.set_ylabel(unit)
    labels = [f'{block}\n{preset}' for block, preset in SCENARIOS]
    ax.set_xticks(range(9), labels, fontsize=9)
    ax.set_ylim(0, max(rates[0] + rates[1]) * 1.19)
    ax.set_axisbelow(True)
    ax.grid(axis='y', alpha=0.25)
    ax.spines[['top', 'right']].set_visible(False)
    add_comparison_legend(ax, a, b)

    if is_relative:
        ax.axhline(100, color=ANNOTATION_COLOR, linewidth=0.8, linestyle='--')


def plot_charts(out_dir: Path, data: dict[str, Any]) -> None:
    '''
    Write RGB, normal-map, and HDR absolute and relative throughput plots.

    Args:
        out_dir: The directory containing the report artifacts.
        data: The settings, dimensions, and measurements.
    '''
    thresholds = sampling_threshold_text(data)

    for image, metadata in data['images'].items():
        tag = image.lower()

        if tag == 'normal':
            tag = 'nml'

        for is_relative in (False, True):
            suffix = 'perf_abs'
            scale = 'absolute'

            if is_relative:
                suffix = 'perf_rel'
                scale = 'relative'

            figure, axes = new_figure((14, 15), 3)
            figure.subplots_adjust(left=0.075, right=0.985, bottom=0.075, top=0.895, hspace=0.57)
            title = metadata['title']
            figure.suptitle(f'{title} | {scale} compression throughput',
                            fontsize=22, y=0.975, fontweight='bold')
            figure.text(0.5, 0.945, f'{data["threads"]} threads | {thresholds}', ha='center')

            for index, pair in enumerate(SETS):
                plot_perf_panel(axes[index], data, image, pair, is_relative)

            output = out_dir / f'{tag}-{suffix}.png'
            figure.savefig(output, dpi=180)
            print(f'Wrote {output}', flush=True)


def quality_values(data: dict[str, Any], image: str, metric: str) \
        -> tuple[list[int], list[float], list[str], set[str]]:
    '''
    Collect signed quality differences and the roles with missing scores.

    Args:
        data: The measured quality samples.
        image: The image type.
        metric: The quality metric to compare.

    Return:
        The bar positions, differences, colors, and roles with missing scores.
    '''
    a, b = SETS[0]
    missing = set()
    positions = []
    values = []
    shades = []

    for pos, (block, preset) in enumerate(SCENARIOS):
        key = f'{image}-{block}-{preset.removeprefix("-")}'
        case = data['cases'][key]
        score_a = mean_quality(case, a, metric)
        score_b = mean_quality(case, b, metric)

        if score_a is None:
            missing.add('A')

        if score_b is None:
            missing.add('B')

        if score_a is None or score_b is None:
            continue

        positions.append(pos)
        values.append(score_b - score_a)
        shades.append(BLOCK_COLORS[block][1])

    return positions, values, shades, missing


def set_change_limits(ax: Any, values: list[float], unit: str) -> None:
    '''
    Set quality-change axis limits with room for zero and the bar labels.

    Args:
        ax: The plot axes receiving the limits.
        values: The available signed quality differences.
        unit: The metric unit, or an empty string for dimensionless scores.
    '''

    if not values:
        return

    low = min(0, *values)
    high = max(0, *values)
    span = high - low

    if span == 0:
        span = 0.001

        if not unit:
            span = 0.000001

    ax.set_ylim(low - 0.35 * span, high + 0.35 * span)


def plot_quality_panel(ax: Any, data: dict[str, Any], image: str,
                       info: tuple[str, str, str, str]) -> None:
    '''
    Draw one signed quality-change comparison.

    Args:
        ax: The plot axes receiving the panel.
        data: The image metadata and measured quality samples.
        image: The image type.
        info: The metric name, display name, unit, and preferred direction.
    '''
    metric, name, unit, direction = info
    positions, values, shades, missing = quality_values(data, image, metric)
    bars = ax.bar(positions, values, width=0.72, color=shades)
    fmt = '%+.4f'

    if not unit:
        fmt = '%+.6f'

    ax.bar_label(bars, fmt=fmt, padding=3, fontsize=9)
    ax.set_title(f'{name} change ({direction} is better)', loc='left', pad=35)
    ylabel = f'Δ{name}'

    if unit:
        ylabel += f' ({unit})'

    ax.set_ylabel(ylabel)
    labels = [f'{block}\n{preset}' for block, preset in SCENARIOS]
    ax.set_xticks(range(len(SCENARIOS)), labels, fontsize=9)
    ax.set_xlim(-0.6, len(SCENARIOS) - 0.4)
    ax.axhline(0, color=ANNOTATION_COLOR, linewidth=0.8)
    set_change_limits(ax, values, unit)
    ax.set_axisbelow(True)
    ax.grid(axis='y', alpha=0.25)
    ax.spines[['top', 'right']].set_visible(False)
    add_comparison_legend(ax, *SETS[0], is_change=True)

    if missing:
        ax.text(0.5, 0.5, 'Unavailable: ' + ', '.join(sorted(missing)),
                transform=ax.transAxes, ha='center', fontsize=10, color=ANNOTATION_COLOR)


def plot_quality_changes(out_dir: Path, data: dict[str, Any]) -> None:
    '''
    Plot two signed quality changes for baseline versus new AVX2.

    Args:
        out_dir: The directory containing the report artifacts.
        data: The image metadata and measured quality samples.
    '''
    thresholds = sampling_threshold_text(data)

    for image, metadata in data['images'].items():
        tag = image.lower()

        if tag == 'normal':
            tag = 'nml'

        figure, axes = new_figure((14, 11), 2)
        figure.subplots_adjust(left=0.09, right=0.985, bottom=0.08, top=0.85, hspace=0.65)
        title = metadata['title']
        figure.suptitle(f'{title} | quality change vs baseline',
                        fontsize=22, y=0.985, fontweight='bold')
        figure.text(0.5, 0.93, 'Signed B − A differences; baseline A = 0 | ' + thresholds,
                    ha='center')

        for ax, info in zip(axes, QUALITY_PLOTS[image]):
            plot_quality_panel(ax, data, image, info)

        output = out_dir / f'{tag}-qual.png'
        figure.savefig(output, dpi=180)
        print(f'Wrote {output}', flush=True)


def describe_images(args: argparse.Namespace) -> dict[str, Any]:
    '''
    Capture the paths, dimensions, and identities of report input images.

    Args:
        args: The image paths selected on the command line.

    Return:
        The image metadata keyed by report image role.
    '''
    hdr_name = args.hdr_image.stem.removeprefix('hdr-rgb-').capitalize()

    images = [
        ('RGB', 'Kodak RGB', args.rgb_image),
        ('Normal', 'Kodak Normal', args.normal_image),
        ('HDR', f'{hdr_name} HDR', args.hdr_image)
    ]

    metadata = {}

    for name, title, image in images:
        path = image.expanduser().resolve()
        width, height = image_size(path)
        metadata[name] = {
            'path': str(path),
            'title': title,
            'sha256': fingerprint(path),
            'width': width,
            'height': height,
        }

    return metadata


def save_report(out_dir: Path, data: dict[str, Any]) -> None:
    '''
    Save raw measurements and generate the report tables and charts.

    Args:
        out_dir: The directory receiving all report artifacts.
        data: The completed measurements and provenance.
    '''
    # Save raw data before rendering so plotting failures leave the measurements available
    raw = json.dumps(data, indent=2, allow_nan=False) + '\n'
    raw_path = out_dir / 'report.json'
    raw_path.write_text(raw, encoding='utf-8')

    write_report(out_dir, data)
    plot_charts(out_dir, data)
    plot_quality_changes(out_dir, data)
    report_path = out_dir / 'report.md'
    print(f'Report: {report_path}\nRaw measurements: {raw_path}')


def main() -> int:
    '''
    Measure the matrix and generate the report artifacts.

    Return:
        Zero on success, one on failure, or 130 on interruption.
    '''
    args = parse_command_line()
    data = {
        'started': datetime.now().astimezone().isoformat(),
        'host': platform.platform(),
        'runs': args.runs,
        'confidence': CONFIDENCE,
        'precision': PRECISION,
        'stopping_rule': 'mean_ci_precision',
        'warmups': args.warmups,
        'threads': args.threads,
        'binaries': {},
        'images': {},
    }

    try:
        # Identify inputs before collecting samples, then verify them before saving
        data['binaries'] = {
            name: binary_metadata(getattr(args, name.replace('-', '_')))
            for name in BINARY_NAMES
        }

        data['images'] = describe_images(args)
        out_dir = args.output.expanduser().resolve()
        out_dir.mkdir(parents=True, exist_ok=True)

        data['cases'] = collect_samples(data['binaries'], data['images'],
                                        args.runs, args.warmups, args.threads,
                                        data['confidence'], data['precision'])

        for metadata in [*data['binaries'].values(), *data['images'].values()]:
            if fingerprint(Path(metadata['path'])) != metadata['sha256']:
                raise ValueError(f'Input changed during measurement: {metadata["path"]}')

        data['finished'] = datetime.now().astimezone().isoformat()
        save_report(out_dir, data)
    except KeyboardInterrupt:
        print('\nBenchmark interrupted', file=sys.stderr)

        return 130
    except PROFILE_ERRORS as error:
        return report_error(error)

    return 0


if __name__ == '__main__':
    sys.exit(main())
