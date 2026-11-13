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
This script is a utility to download test images that are available online, avoiding the need to
duplicate them in the repository.
'''

from pathlib import Path
import sys
from typing import Optional
import urllib.request

import numpy
from PIL import Image


TEST_IMAGE_DIR = Path('Test/Images')

NORMAL_MAP_STRENGTH = 8.0

# Alpha image for each RGB image in the non-correlated alpha set. This is a stable derangement, so
# every source is used exactly once and never with its matching RGB image.
KODAK_NCA_MAP = {
    1: 2,
    2: 23,
    3: 17,
    4: 10,
    5: 19,
    6: 5,
    7: 15,
    8: 3,
    9: 18,
    10: 1,
    11: 12,
    12: 16,
    13: 9,
    14: 21,
    15: 14,
    16: 13,
    17: 6,
    18: 8,
    19: 7,
    20: 11,
    21: 22,
    22: 24,
    23: 20,
    24: 4,
}


def download(test_set: str, index: int, src_url: str, dst_path: Path) -> None:
    '''
    Download a single image.

    Args:
        test_set: The test set name.
        index: The download image index.
        src_url: The download URL.
        dst_path: The destination file path.
    '''
    # Skip downloads if the file already exists
    if dst_path.exists():
        print(f'{test_set} image {index}: Skipping')
        return

    dir_name = dst_path.parent
    dir_name.mkdir(parents=True, exist_ok=True)

    print(f'{test_set} image {index}: Downloading')
    urllib.request.urlretrieve(src_url, dst_path)


def make_landscape(file_path: Path) -> None:
    '''
    Modify an existing image on disk to make it landscape aspect.

    Args:
        file_path: The path of the image on disk.
    '''
    image = Image.open(file_path)
    if image.size[0] < image.size[1]:
        rotated_image = image.rotate(90, expand=True)
        rotated_image.save(file_path)


def make_mixed_image(srca_path: Path, srcb_path: Path, dst_path: Path) -> None:
    '''
    Make image consisting of RGB from A, and alpha from B.

    Args:
        srca_path: The path of input A on disk.
        srcb_path: The path of input B on disk.
        dst_path: The path of the destination.
    '''
    dir_name = dst_path.parent
    dir_name.mkdir(parents=True, exist_ok=True)

    image_rgb = Image.open(srca_path)
    image_alpha = Image.open(srcb_path).convert('L')

    image_rgb.putalpha(image_alpha)
    image_rgb.save(dst_path)


def make_normal_map(image: Image.Image) -> Image.Image:
    '''
    Convert an image to an XYZ unit-vector normal map.

    Args:
        image: The source image.

    Return:
        The RGB encoded normal map.
    '''
    height = numpy.asarray(image.convert('L'), dtype=numpy.float32) / 255.0
    gradient_y, gradient_x = numpy.gradient(height)

    normal_x = -gradient_x * NORMAL_MAP_STRENGTH
    normal_y = -gradient_y * NORMAL_MAP_STRENGTH
    normal_z = numpy.ones_like(height)

    length = numpy.sqrt(normal_x * normal_x + normal_y * normal_y + normal_z * normal_z)

    normal = numpy.stack((normal_x, normal_y, normal_z), axis=2)
    normal /= length[:, :, numpy.newaxis]

    encoded = numpy.rint((normal * 0.5 + 0.5) * 255.0)
    return Image.fromarray(encoded.astype(numpy.uint8))


def make_montage(src_dir_path: Path, dst_path: Path,
                 alpha_images: Optional[list[Path]] = None,
                 image_format: str = 'RGB') -> None:
    '''
    Make a single mosaic montage consisting of all of the Kodak images.

    Args:
        src_dir_path: The directory path of the Kodak images on disk.
        dst_path: The file path of the resulting montage.
        alpha_images: Optional image paths to use as alpha sources.
        image_format: The output image format, one of L, RGB, or XY.
    '''
    cols = 6
    rows = 4

    width = 768
    height = 512

    images = list(src_dir_path.glob('*.png'))
    images.sort()

    if image_format not in ('L', 'RGB', 'XY'):
        raise ValueError('Unsupported montage image format')

    if alpha_images is not None and image_format != 'RGB':
        raise ValueError('Alpha images require RGB image format')

    if alpha_images is not None and len(alpha_images) != len(images):
        raise ValueError('RGB and alpha image counts must match')

    montage = Image.new(
        'RGBA' if alpha_images is not None else
        ('RGB' if image_format == 'XY' else image_format),
        (width * cols, height * rows))

    for i, src_path in enumerate(images):
        image = Image.open(src_path).convert('RGB')
        if image_format == 'L':
            image = image.convert('L')
        elif image_format == 'XY':
            image = make_normal_map(image)

        if alpha_images is not None:
            image_alpha = Image.open(alpha_images[i]).convert('L')
            image.putalpha(image_alpha)

        col = i % cols
        row = i // cols
        montage.paste(image, (width * col, height * row))

    dst_path.parent.mkdir(parents=True, exist_ok=True)

    montage.save(dst_path)


def retrieve_kodak_set():
    '''
    Download the public domain Kodak image set.

    To make test set mosaics easier to build we rotate images to make everything landscape.
    '''
    test_set = 'Kodak'

    # Download the original RGB images
    for i in range(1, 25):
        src_url = f'http://r0k.us/graphics/kodak/kodak/kodim{i:02}.png'

        dst_file = f'ldr-rgb-kodak{i:02}.png'
        dst_path = TEST_IMAGE_DIR / 'Kodak' / 'LDR-RGB' / dst_file
        download(test_set, i, src_url, dst_path)

        # Canonicalize image aspect
        make_landscape(dst_path)

    # Make some correlated alpha RGBA images
    src_dir = TEST_IMAGE_DIR / 'Kodak' / 'LDR-RGB'
    dst_dir = TEST_IMAGE_DIR / 'KodakSim' / 'LDR-RGBA'

    for i in (22, 23):
        src_path = src_dir / f'ldr-rgb-kodak{i:02}.png'
        dst_path = dst_dir / f'ldr-rgba-kodak{i:02}+ca.png'
        make_mixed_image(src_path, src_path, dst_path)

    # Make some non-correlated alpha RGBA images
    for i in (22, 23):
        j = KODAK_NCA_MAP[i]
        srca_path = src_dir / f'ldr-rgb-kodak{i:02}.png'
        srcb_path = src_dir / f'ldr-rgb-kodak{j:02}.png'

        dst_path = dst_dir / f'ldr-rgba-kodak{i:02}+{j:02}+nca.png'
        make_mixed_image(srca_path, srcb_path, dst_path)

    # Make a large RGB montage
    dst_path = TEST_IMAGE_DIR / 'KodakMnt' / 'LDR-RGB' / 'ldr-rgb-montage.png'
    make_montage(src_dir, dst_path)

    # Make a large correlated alpha montage
    images = list(src_dir.glob('*.png'))
    images.sort()

    dst_dir = TEST_IMAGE_DIR / 'KodakMnt' / 'LDR-RGBA'
    dst_path = dst_dir / 'ldr-rgba-montage+ca.png'
    make_montage(src_dir, dst_path, images)

    # Make a large non-correlated alpha montage
    alpha_images = [src_dir / f'ldr-rgb-kodak{KODAK_NCA_MAP[i]:02}.png' for i in range(1, 25)]
    dst_path = dst_dir / 'ldr-rgba-montage+nca.png'
    make_montage(src_dir, dst_path, alpha_images)

    # Make a large luminance montage
    dst_path = TEST_IMAGE_DIR / 'KodakMnt' / 'LDR-L' / 'ldr-l-montage.png'
    make_montage(src_dir, dst_path, image_format='L')

    # Make a large XYZ normal map montage
    dst_path = TEST_IMAGE_DIR / 'KodakMnt' / 'LDR-XY' / 'ldr-xy-montage.png'
    make_montage(src_dir, dst_path, image_format='XY')


def main() -> int:
    '''
    The main function.

    Return:
        The process return code.
    '''
    retrieve_kodak_set()
    return 0


if __name__ == '__main__':
    sys.exit(main())
