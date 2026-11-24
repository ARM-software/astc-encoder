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
Plot KodakMnt RGB RDO curves for lambda 0 through 100 in steps of 10.

Defaults to the AVX2 build, -medium, 12 threads, and 4x4, 5x5, 6x6, and 8x8 blocks.
Plots RGB PSNR against raw DEFLATE level 9 bits per texel, excluding the ASTC header.
Use --metric rmse to plot RGB RMSE on the 0 to 255 scale instead.
The graph uses the report's block colors, with purple for 5x5, and is 800 pixels wide.
All codec invocations run serially. Intermediate files are temporary; only the PNG is saved.

Example:
    python3 Test/astc_profile_report_rdo.py --output x_profiles/rdo-rgb.png
'''

import argparse
import math
from pathlib import Path
import sys
import tempfile
import zlib

from PIL import Image, ImageChops

from astc_profile_report import ANNOTATION_COLOR, BLOCK_COLORS, PLOT_WIDTH, new_figure
from testlib.profile import PROFILE_ERRORS, report_error, resolve_binary, run_codec


ROOT = Path(__file__).resolve().parent.parent
BLOCK_SIZES = ('4x4', '5x5', '6x6', '8x8')
LAMBDAS = range(0, 101, 10)


def measure_case(binary: str, image: Path, original: Image.Image, block: str,
                 strength: int, args: argparse.Namespace,
                 temporary: Path) -> tuple[float, float]:
    '''
    Encode one case and measure the saved bitstream's compressed size and RGB error.

    Args:
        binary: The codec executable path.
        image: The source image path.
        original: The source image converted to RGB.
        block: The ASTC block footprint.
        strength: The RDO lambda value.
        args: The compression preset, thread count, and quality metric.
        temporary: The directory for temporary encoded and decoded files.

    Return:
        The DEFLATE bitrate in bits per texel and RGB PSNR in dB or RMSE on the 0 to 255 scale.
    '''
    encoded = temporary / 'image.astc'
    decoded = temporary / 'image.png'
    run_codec(binary, ['-cl', str(image), str(encoded), block, f'-{args.preset}',
                       '-rdo-lambda', str(strength), '-j', str(args.threads)])
    run_codec(binary, ['-dl', str(encoded), str(decoded), '-j', str(args.threads)])

    # Match the CLI storage metric by excluding the 16-byte ASTC file header
    compressor = zlib.compressobj(level=9, wbits=-zlib.MAX_WBITS)
    compressed = compressor.compress(encoded.read_bytes()[16:]) + compressor.flush()
    rate = len(compressed) * 8 / (original.width * original.height)

    # Measure the decoded file so quality and size describe the same encoding
    with Image.open(decoded) as restored:
        histogram = ImageChops.difference(original, restored.convert('RGB')).histogram()

    squared_error = sum(count * (index % 256) ** 2 for index, count in enumerate(histogram))
    mse = squared_error / (original.width * original.height * 3)
    if args.metric == 'rmse':
        return rate, math.sqrt(mse)

    if mse <= 0.0:
        raise ValueError('Cannot plot infinite PSNR for a lossless encoding')

    return rate, 10 * math.log10(255.0 ** 2 / mse)


def parse_command_line() -> argparse.Namespace:
    '''
    Parse the sweep settings and output path.

    Return:
        The validated command-line arguments.
    '''
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--binary', default=str(ROOT / 'bin/astcenc-avx2'))
    parser.add_argument('--image', type=Path,
                        default=ROOT / 'Test/Images/KodakMnt/LDR-RGB/ldr-rgb-montage.png')
    parser.add_argument('--preset', choices=('fastest', 'fast', 'medium', 'thorough',
                                             'verythorough', 'exhaustive'), default='medium')
    parser.add_argument('--threads', type=int, default=12)
    parser.add_argument('--metric', choices=('psnr', 'rmse'), default='psnr',
                        help='Quality metric: PSNR in dB (default) or RMSE on the 0 to 255 scale')
    parser.add_argument('--output', type=Path, default=ROOT / 'x_profiles/rdo-rgb.png')
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('--threads must be positive')

    return args


def main() -> int:
    '''
    Run the complete sweep and save an 800-pixel-wide graph.

    Return:
        Zero on success, or one if measurement or plotting fails.
    '''
    args = parse_command_line()
    try:
        binary = resolve_binary(args.binary)
        with Image.open(args.image) as image:
            original = image.convert('RGB')

        colors = {block: pair[1] for block, pair in BLOCK_COLORS.items()}
        colors['5x5'] = '#8767a5'  # Purple distinguishes the additional footprint
        figure, axes = new_figure((12, 7.5), 1)
        ax = axes[0]
        with tempfile.TemporaryDirectory(prefix='astc-rdo-') as temporary:
            for block in BLOCK_SIZES:
                points = []
                for strength in LAMBDAS:
                    print(f'{block}, lambda {strength}: encoding', flush=True)
                    points.append(measure_case(binary, args.image, original, block, strength,
                                               args, Path(temporary)))

                rates, scores = zip(*points)
                color = colors[block]
                ax.plot(rates, scores, '-o', color=color, linewidth=2, markersize=4, label=block)
                ax.plot(rates[0], scores[0], 'o', color=color,
                        markerfacecolor='white', markersize=7)
                ax.plot(rates[-1], scores[-1], 's', color=color, markersize=6)

        figure.suptitle('KodakMnt RGB · RDO sweep', fontsize=22, fontweight='bold', y=0.97)
        ax.set_title(f'AVX2 · -{args.preset} · {args.threads} threads',
                     fontsize=12, color=ANNOTATION_COLOR, pad=16)
        ax.set_xlabel('DEFLATE-compressed bitrate (bits / texel)', fontsize=12)
        if args.metric == 'rmse':
            ax.set_ylabel('RGB RMSE (0–255 scale)', fontsize=12)
        else:
            ax.set_ylabel('RGB PSNR (dB)', fontsize=12)

        ax.grid(alpha=0.25)
        ax.set_axisbelow(True)
        ax.spines[['top', 'right']].set_visible(False)
        ax.margins(x=0.06, y=0.08)
        ax.legend(title='ASTC block size', frameon=False)
        figure.text(0.5, 0.035,
                    'λ = 0 to 100 in steps of 10 · Open circles: λ = 0 · Squares: λ = 100\n'
                    'Raw DEFLATE level 9 over ASTC payload; headers excluded',
                    ha='center', fontsize=10, color=ANNOTATION_COLOR)
        figure.subplots_adjust(left=0.10, right=0.97, bottom=0.18, top=0.84)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(args.output, dpi=PLOT_WIDTH / figure.get_figwidth())
        print(f'Saved graph to {args.output}')
    except PROFILE_ERRORS as error:
        return report_error(error)

    return 0


if __name__ == '__main__':
    sys.exit(main())
