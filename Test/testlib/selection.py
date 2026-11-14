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
Deterministic image test selection with full or pairwise coverage.
'''

from collections.abc import Sequence
from dataclasses import dataclass

from .testset import TestImage, TestSet


TEST_BLOCK_SIZES = ('4x4', '5x5', '6x6', '8x8', '12x12', '3x3x3', '6x6x6')
TEST_QUALITIES = ('fastest', 'fast', 'medium', 'thorough', 'verythorough', 'exhaustive')


@dataclass(frozen=True)
class ImageCase:
    '''
    One image and its compression settings in a test plan.

    Attributes:
        image: The source image and its metadata.
        block_size: The ASTC block size.
        quality: The quality preset without its leading dash.
    '''
    image: TestImage
    block_size: str
    quality: str


def is_3d(block_size: str) -> bool:
    '''
    Is the given block size a 3D block type?

    Args:
        block_size: The block size.

    Return:
        True if the block size is 3D, False if 2D.
    '''
    return block_size.count('x') == 2


def select_test_cases(test_set: TestSet, block_sizes: Sequence[str],
                      qualities: Sequence[str], pairwise: bool = False) -> list[ImageCase]:
    '''
    Select image, block size, and quality combinations to test.

    Pairwise selection groups images by color profile and format, then covers every feasible pair of
    image group, block size, and quality. Every image with a compatible block size is tested at
    least once. Ties are resolved in input order to keep selection deterministic.

    Args:
        test_set: The test set to run.
        block_sizes: The block sizes to run.
        qualities: The quality levels to run.
        pairwise: Select pairwise coverage instead of all combinations?

    Return:
        The selected image cases in deterministic execution order.
    '''
    cases = [
        ImageCase(image=image, block_size=block_size, quality=quality)
        for block_size in block_sizes
        for image in test_set.tests
        if is_3d(block_size) == image.is_3d
        for quality in qualities
    ]

    if not pairwise:
        return cases

    coverage = [
        {
            ('image', str(case.image.file_path), ''),
            ('format-block', (case.image.color_profile, case.image.color_format),
             case.block_size),
            ('format-quality', (case.image.color_profile, case.image.color_format),
             case.quality),
            ('block-quality', case.block_size, case.quality)
        }
        for case in cases
    ]

    uncovered = set().union(*coverage)
    selected = set()

    while uncovered:
        index = max(range(len(cases)), key=lambda i: len(coverage[i] & uncovered))
        selected.add(index)
        uncovered.difference_update(coverage[index])

    return [case for index, case in enumerate(cases) if index in selected]
