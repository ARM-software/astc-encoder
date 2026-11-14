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
Compare ASTC file compressibility and image quality for two astcenc binaries.

Rate distortion optimization (RDO) aims to change the ASTC encoding to make the data file more
compressible when bundled into a distribution format, such as an Android APK or asset bundle. This
script tests the impact of RDO on both data compressibility and image quality.

A is the baseline and B is the candidate. Both currently receive identical encoding settings;
candidate-only RDO options will be added when the encoder supports them. Each case runs compression
to ASTC and compression/decompression to PNG for both binaries, then compresses the ASTC files
using ZIP and LZ4. All sizes include file headers.

Profiling timing is disabled by default. With --repeats N, run each ASTC encoding N times and
report its mean compression time (including any RDO pass) in seconds.
'''

import argparse
import csv
from dataclasses import dataclass
import io
import math
from pathlib import Path
import statistics
import subprocess as sp
import sys
import zipfile

import lz4.frame

from testlib.profile import select_quality_metrics


BLOCK_SIZES = ('4x4', '5x4', '5x5', '6x5', '6x6', '8x5', '8x6', '10x5', '10x6', '8x8')

PROFILES = {
    'ldr': 'l',
    'srgb': 's',
    'hdr': 'h',
    'hdra':
    'H'
}

QUALITIES = ('fastest', 'fast', 'medium', 'thorough', 'verythorough', 'exhaustive')

COLUMNS = (
    ('image', None, 'Image', 0),
    ('block_size', None, 'Block', 0),
    ('compressor', None, 'Codec', 0),
    ('compressed_b_kb', 'compressed_diff_kb', 'Cmp Sz (KB)', 2),
    ('compressed_b_pct', 'compressed_diff_pct', 'Cmp %', 2),
    ('time_b', 'time_diff', 'Time (s)', 4),
    ('psnr_b', 'psnr_diff', 'PSNR (dB)', 2),
    ('ssim_b', 'ssim_diff', 'SSIM', 2),
)

METADATA = (
    'source_image',
    'binary_a',
    'binary_b',
    'format_a',
    'format_b',
    'profile',
    'quality',
    'threads',
    'repeats',
    'zip_level',
    'lz4_level'
)


@dataclass(frozen=True)
class EncodingResult:
    '''
    An ASTC bitstream and the quality reported by its matching test-mode invocation.

    Attributes:
        data: The complete ASTC or ASTCM file, including its header.
        quality: The PSNR and SSIM scores keyed by their CLI labels.
        coding_time: The mean ASTC encoding time in seconds, or None when timing is disabled.
    '''
    data: bytes
    quality: dict[str, float]
    coding_time: float | None = None


def parse_quality(output: str) -> dict[str, float]:
    '''
    Read PSNR and SSIM scores without depending on reported processing times.

    Args:
        output: The CLI test-mode output.

    Return:
        The available PSNR and SSIM scores keyed by label.

    Raises:
        ValueError: A reported score is invalid.
    '''
    quality = {}
    labels = ('PSNR (LDR-RGB)', 'PSNR (LDR-RGBA)', 'mPSNR (RGB)',
              'SSIM (LDR-RGB)', 'SSIM (LDR-RGBA)', 'mSSIM (RGB)')

    for line in output.splitlines():
        name, _, value = line.strip().partition(':')

        if name in labels:
            if not value.strip():
                raise ValueError(f'Codec reported an empty {name}')

            score = float(value.split()[0])

            if not math.isfinite(score):
                raise ValueError(f'Codec reported a non-finite {name}')

            quality[name] = score

    return quality


def parse_coding_time(output: str) -> float:
    '''
    Read the CLI encoding time, allowing short runs rounded to zero seconds.

    Args:
        output: The CLI compression-mode output.

    Return:
        The reported encoding time in seconds.

    Raises:
        ValueError: The encoding time is missing, negative, or non-finite.
    '''
    for line in output.splitlines():
        name, _, value = line.strip().partition(':')

        if name == 'Coding time':
            fields = value.split()

            if len(fields) != 2 or fields[1] != 's':
                raise ValueError('Codec reported an invalid Coding time')

            coding_time = float(fields[0])

            if not math.isfinite(coding_time) or coding_time < 0.0:
                raise ValueError('Codec reported an invalid Coding time')

            return coding_time

    raise ValueError('Codec did not report Coding time')


def encode_image(args: argparse.Namespace, binary: str, image: Path,
                 block_size: str, output_path: Path) -> EncodingResult:
    '''
    Encode a texture file, optionally measure repeated encodings, and save a PNG and quality report.

    Args:
        args: The shared encoding settings.
        binary: The codec executable.
        image: The source image.
        block_size: The ASTC block size.
        output_path: The output bitstream path, with an .astc or .astcm extension.

    Return:
        The encoded file bytes, reported quality scores, and optional mean encoding time.

    Raises:
        CalledProcessError: A codec invocation failed; its output is saved in the matching log.
    '''
    output = ''
    profile = PROFILES[args.profile]
    coding_times = []

    for operation, target in (('c', output_path), ('t', output_path.with_suffix('.png'))):
        repeats = 1

        if operation == 'c' and args.repeats > 0:
            repeats = args.repeats

        for repeat in range(repeats):
            command = [binary, f'-{operation}{profile}', str(image),
                       str(target), block_size, f'-{args.quality}',
                       '-j', str(args.threads)]
            result = sp.run(command, stdout=sp.PIPE, stderr=sp.STDOUT, text=True, check=False)
            log_suffix = f'.{operation}.log'

            if operation == 'c' and args.repeats > 0:
                log_suffix = f'.c.{repeat + 1}.log'

            output_path.with_suffix(log_suffix).write_text(result.stdout, encoding='utf-8')
            result.check_returncode()
            output = result.stdout

            if operation == 'c' and args.repeats > 0:
                coding_times.append(parse_coding_time(output))

    coding_time = None

    if coding_times:
        coding_time = statistics.mean(coding_times)

    return EncodingResult(data=output_path.read_bytes(),
                          quality=parse_quality(output), coding_time=coding_time)


def compress_astc(data: bytes, zip_level: int, lz4_level: int) -> dict[str, bytes]:
    '''
    Compress a complete ASTC or ASTCM file into a deterministic ZIP archive and an LZ4 frame.

    Args:
        data: The encoded file bytes.
        zip_level: The ZIP Deflate compression level.
        lz4_level: The LZ4 compression level.

    Return:
        The complete compressed files keyed by their format names.
    '''
    buffer = io.BytesIO()

    # A fixed neutral name prevents filenames, extensions, or timestamps from biasing comparisons
    with zipfile.ZipFile(buffer, 'w') as archive:
        archive.writestr(zipfile.ZipInfo('texture.bin'), data,
                         compress_type=zipfile.ZIP_DEFLATED, compresslevel=zip_level)

    frame = lz4.frame.compress(data, compression_level=lz4_level,
                               block_size=lz4.frame.BLOCKSIZE_MAX64KB, block_linked=True,
                               content_checksum=False, block_checksum=False, store_size=True)

    return {'ZIP': buffer.getvalue(), 'LZ4': frame}


def compare_quality(image: Path, profile: str, result_a: EncodingResult,
                    result_b: EncodingResult) -> list[dict[str, str | float]]:
    '''
    Compare matching PSNR and SSIM scores using the existing image-test channel policy.

    Args:
        image: The source image.
        profile: The ASTC color profile.
        result_a: The baseline encoding result.
        result_b: The candidate encoding result.

    Return:
        The quality comparison for the selected channel group.

    Raises:
        ValueError: The binaries do not both report the required PSNR and SSIM scores.
    '''
    metrics = select_quality_metrics([f'-t{PROFILES[profile]}', str(image)], result_a.quality)

    if not metrics:
        raise ValueError('Codec did not report PSNR scores')

    rows: list[dict[str, str | float]] = []

    # Prefer RGBA when a generic filename reports both channel groups
    for psnr_label in metrics[-1:]:
        ssim_label = psnr_label.replace('PSNR', 'SSIM')

        for side, result in (('A', result_a), ('B', result_b)):
            for label in (psnr_label, ssim_label):
                if label not in result.quality:
                    raise ValueError(f'{side} did not report {label}; '
                                     'use a binary with PSNR and SSIM support')

        psnr_a = result_a.quality[psnr_label]
        psnr_b = result_b.quality[psnr_label]
        ssim_a = result_a.quality[ssim_label]
        ssim_b = result_b.quality[ssim_label]

        rows.append({
            'psnr_a': psnr_a, 'psnr_b': psnr_b, 'psnr_diff': psnr_b - psnr_a,
            'ssim_a': ssim_a, 'ssim_b': ssim_b, 'ssim_diff': ssim_b - ssim_a,
        })

    return rows


def compare_case(args: argparse.Namespace, image: Path, block_size: str,
                 directory: Path) -> list[dict[str, str | int | float]]:
    '''
    Run both binaries and compare their compressed sizes and matching quality metrics.

    Args:
        args: The benchmark settings.
        image: The source image.
        block_size: The ASTC block size.
        directory: The output directory for this case.

    Return:
        One report row per lossless compressor and selected channel group.
    '''
    directory.mkdir(parents=True, exist_ok=True)
    result_a = encode_image(args, args.binary_a, image, block_size,
                            directory / f'a.{args.format_a}')
    result_b = encode_image(args, args.binary_b, image, block_size,
                            directory / f'b.{args.format_b}')
    qualities = compare_quality(image, args.profile, result_a, result_b)
    archives_a = compress_astc(result_a.data, args.zip_level, args.lz4_level)
    archives_b = compress_astc(result_b.data, args.zip_level, args.lz4_level)
    timing = {}

    if args.repeats > 0:
        assert result_a.coding_time is not None
        assert result_b.coding_time is not None
        timing = {'time_a': result_a.coding_time, 'time_b': result_b.coding_time,
                  'time_diff': result_b.coding_time - result_a.coding_time}

    rows: list[dict[str, str | int | float]] = []

    for compressor, data_a in archives_a.items():
        data_b = archives_b[compressor]
        (directory / f'a.{args.format_a}.{compressor.lower()}').write_bytes(data_a)
        (directory / f'b.{args.format_b}.{compressor.lower()}').write_bytes(data_b)

        for quality in qualities:
            rows.append({
                **quality,
                **timing,
                'image': image.name,
                'block_size': block_size,
                'compressor': compressor,
                'compressed_a_kb': len(data_a) / 1024.0,
                'compressed_b_kb': len(data_b) / 1024.0,
                'compressed_diff_kb': (len(data_b) - len(data_a)) / 1024.0,
                'compressed_a_pct': 100.0 * len(data_a) / len(result_a.data),
                'compressed_b_pct': 100.0 * len(data_b) / len(result_b.data),
                'compressed_diff_pct': (100.0 * len(data_b) / len(result_b.data)
                                        - 100.0 * len(data_a) / len(result_a.data)),
                'source_image': str(image.resolve()),
                'binary_a': args.binary_a,
                'binary_b': args.binary_b,
                'format_a': args.format_a,
                'format_b': args.format_b,
                'profile': args.profile,
                'quality': args.quality,
                'threads': args.threads,
                'repeats': args.repeats,
                'zip_level': args.zip_level,
                'lz4_level': args.lz4_level,
            })

    return rows


def format_comparison(value: float, difference: float, precision: int) -> str:
    '''
    Format a B value with its signed difference from A in brackets.

    Args:
        value: The candidate B value.
        difference: The candidate value minus the baseline A value.
        precision: The number of decimal places for both values.

    Return:
        The candidate value followed by its signed difference in brackets.
    '''
    sign = '+'

    if difference < 0.0:
        sign = '-'

    return f'{value:.{precision}f} ({sign} {abs(difference):.{precision}f})'


def report_results(rows: list[dict[str, str | int | float]], csv_path: Path | None) -> None:
    '''
    Print the comparison table and optionally save the results and settings as CSV.

    Args:
        rows: The benchmark results.
        csv_path: The optional destination CSV file.
    '''
    is_timed = bool(rows and 'time_b' in rows[0])
    columns = [column for column in COLUMNS if is_timed or column[0] != 'time_b']
    table = [[title for _, _, title, _ in columns]]

    for row in rows:
        cells = []

        for key, difference_key, _, precision in columns:
            if difference_key is not None:
                cell = format_comparison(float(row[key]), float(row[difference_key]), precision)
            else:
                cell = str(row[key])

            cells.append(cell)

        table.append(cells)

    widths = [max(len(cell) for cell in column) for column in zip(*table)]
    print('Values are B (B - A); 1 KB = 1024 bytes.')
    print('Cmp % is compressed / encoded file size x 100; changes are percentage points.')

    if is_timed:
        print('Times are mean ASTC encoding seconds from the CLI Coding time report.')

    for cells in table:
        print('| ' + ' | '.join(cell.rjust(width) for cell, width in zip(cells, widths)) + ' |')

    if csv_path is not None:
        csv_path.parent.mkdir(parents=True, exist_ok=True)

        with csv_path.open('w', encoding='utf-8', newline='') as handle:
            fields = []

            for key, difference_key, _, _ in columns:
                fields.append(key)

                if difference_key is not None:
                    fields.append(difference_key)

            fields.extend(METADATA)
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()

            for row in rows:
                record = {key: row[key] for key in fields}

                for key, value in record.items():
                    if isinstance(value, float) and key.startswith(('psnr_', 'ssim_')):
                        record[key] = f'{value:.2f}'

                writer.writerow(record)


def parse_command_line(argv: list[str] | None = None) -> argparse.Namespace:
    '''
    Parse the input images and shared A/B encoding and compression settings.

    Args:
        argv: Arguments to parse, or None to use the process command line.

    Return:
        The validated command-line settings.
    '''
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        'binary_a', type=argparse.FileType('rb'),
        help='baseline astcenc executable')

    parser.add_argument(
        'binary_b', type=argparse.FileType('rb'),
        help='candidate astcenc executable; may be identical to A')

    parser.add_argument(
        'images', nargs='+', type=Path, help='input images')

    parser.add_argument(
        '--format-a', choices=('astc', 'astcm'), default='astc',
        help='baseline bitstream file format (default: %(default)s)')

    parser.add_argument(
        '--format-b', choices=('astc', 'astcm'), default='astcm',
        help='candidate bitstream file format (default: %(default)s)')

    parser.add_argument(
        '--block-size', dest='block_sizes', action='append', choices=BLOCK_SIZES,
        help='block size to test (repeatable, default: all < 64 texels)')

    parser.add_argument(
        '--quality', choices=QUALITIES, default='medium',
        help='encoding quality preset (default: %(default)s)')

    parser.add_argument(
        '--profile', choices=PROFILES, default='ldr',
        help='compressor color profile (default: %(default)s)')

    parser.add_argument(
        '-j', '--threads', type=int, default=12,
        help='compressor threads (default: %(default)s)')

    parser.add_argument(
        '--repeats', type=int, default=0,
        help='test repeats for timing, 0 disables timing (default: %(default)s)')

    parser.add_argument(
        '--zip-level', type=int, choices=range(10), default=7,
        help='ZIP Deflate compression level (default: %(default)s)')

    parser.add_argument(
        '--lz4-level', type=int, choices=range(17), default=9,
        help='LZ4 compression level, 0 is fast and 3+ is HC (default: %(default)s)')

    # TODO: Remove this and make a TempDir in future
    parser.add_argument(
        '--output-dir', type=Path, default=Path('TestOutput/RDO'),
        help='directory for bitstreams, PNGs, archives, and logs (default: %(default)s)')

    parser.add_argument(
        '--csv',
        type=Path, help='optional output CSV with results and settings')

    args = parser.parse_args(argv)

    for name in ('binary_a', 'binary_b'):
        # Retain the filename and close the FileType handle when the "with" block exits
        with getattr(args, name) as binary:
            setattr(args, name, binary.name)

    if args.threads < 1:
        parser.error('--threads must be positive')

    if args.repeats < 0:
        parser.error('--repeats cannot be negative')

    args.block_sizes = list(dict.fromkeys(args.block_sizes or BLOCK_SIZES))

    for image in args.images:
        if not image.is_file():
            parser.error(f'cannot find input image: {image}')

    return args


def main() -> int:
    '''
    Run each image and block-size comparison sequentially and report the results.

    Return:
        Zero on success, or one if encoding, compression, or reporting fails.
    '''
    args = parse_command_line()
    rows = []

    print(f'A: {args.binary_a} (.{args.format_a})\n'
          f'B: {args.binary_b} (.{args.format_b})', file=sys.stderr)
    print(f'Profile: {args.profile}, quality: {args.quality}, threads: {args.threads}; '
          f'ZIP level: {args.zip_level}, LZ4 level: {args.lz4_level}', file=sys.stderr)

    try:
        for index, image in enumerate(args.images):
            for block_size in args.block_sizes:
                print(f'Encoding {image} at {block_size}', file=sys.stderr, flush=True)
                directory = args.output_dir / f'{index:04d}-{image.stem}' / block_size
                rows.extend(compare_case(args, image, block_size, directory))

        report_results(rows, args.csv)
    except (OSError, ValueError, RuntimeError, sp.CalledProcessError) as ex:
        print(f'ERROR: {ex}', file=sys.stderr)

        if isinstance(ex, sp.CalledProcessError) and ex.output:
            print(ex.output, file=sys.stderr)

        return 1

    return 0


if __name__ == '__main__':
    sys.exit(main())
