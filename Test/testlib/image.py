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
An Image provides a wrapper around an image file on the filesystem which allows queries to be made
against that image via a Python API.

ImageMagick is used to provide some of the heavy lifting, as it supports more file formats than PIL.
'''

import os
from pathlib import Path
import subprocess as sp

Color = tuple[float, float, float, float]

CONVERT_BINARY = ['convert']


# ImageMagick emits names for some exact colors instead of channel tuples
NAMED_COLORS: dict[str, Color] = {
    'black': (0.0, 0.0, 0.0, 1.0),
    'white': (1.0, 1.0, 1.0, 1.0),
    'red': (1.0, 0.0, 0.0, 1.0),
    'blue': (0.0, 0.0, 1.0, 1.0),
}


def decode_color(raw_color: str) -> Color:
    '''
    Decode an ImageMagick named color or RGB/RGBA channel tuple.

    Args:
        raw_color: The pixel value emitted by ImageMagick.

    Return:
        The normalized RGBA channels, with opaque alpha when none is reported.

    Raises:
        AssertionError: The color representation is not supported.
    '''
    if raw_color in NAMED_COLORS:
        return NAMED_COLORS[raw_color]

    kind, _, channels = raw_color.partition('(')
    assert kind in ('rgb', 'rgba', 'srgb', 'srgba'), f'Unknown color format {raw_color}'
    assert channels.endswith(')'), f'Invalid color tuple {raw_color}'
    values = channels[:-1].split(',')
    color = [0.0] * 4

    # RGB percentages and byte values use different scales; alpha is already normalized
    for index, value in enumerate(values):
        if index >= 3:
            color[index] = float(value)
        elif value.endswith('%'):
            color[index] = float(value[:-1]) / 100.0
        else:
            color[index] = float(value) / 255.0

    if kind in ('rgb', 'srgb'):
        color[3] = 1.0

    return (color[0], color[1], color[2], color[3])


class Image:
    '''
    Wrapper around an image on the file system.
    '''

    # No support for KTX yet because ImageMagick doesn't support it
    SUPPORTED_LDR = ['bmp', 'dds', 'jpg', 'png', 'tga']
    SUPPORTED_HDR = ['exr', 'hdr']

    @classmethod
    def is_format_supported(cls, file_format: str, profile: str):
        '''
        Test if a given file format is supported by the library.

        Args:
            file_format: The target file extension, excluding the '.'.
            profile: The color profile (ldr or hdr) of the image.

        Return:
            True if the image file format is supported, False otherwise.
        '''
        assert profile in ['ldr', 'hdr']

        if profile == 'ldr':
            return file_format in cls.SUPPORTED_LDR

        return file_format in cls.SUPPORTED_HDR

    def __init__(self, file_path: Path):
        '''
        Construct a new Image wrapper.

        Args:
            file_path: The path of the image on disk.
        '''
        self.file_path = file_path

    def get_colors(self, coord: tuple[int, int]) -> Color:
        '''
        Get the image colors at the given coordinate.

        Args:
            coord: An image coordinate to sample.

        Return:
            The sampled color. Colors are returned as float values between 0.0 and 1.0 for LDR, and
            float values which may exceed 1.0 for HDR.
        '''
        x = coord[0]
        y = coord[1]

        command = list(CONVERT_BINARY)
        command.append(str(self.file_path))

        # Ensure convert factors in format origin if needed
        command.append('-auto-orient')

        command += ['-format', '%%[pixel:p{%u,%u}]' % (x, y), 'info:']

        if os.name == 'nt':
            command.insert(0, 'magick')

        result = sp.run(command, stdout=sp.PIPE, stderr=sp.PIPE, check=True, text=True)

        color = list(decode_color(result.stdout.strip()))
        suffix = self.file_path.suffix.lower()

        # ImageMagick decodes DDS as BGRA not RGBA; manually correct
        if suffix == '.dds':
            color[0], color[2] = color[2], color[0]

        # EXR samples are premultiplied; fully transparent pixels cannot be unpremultiplied
        if (suffix == '.exr') and (color[3] != 0.0):
            color[0] /= color[3]
            color[1] /= color[3]
            color[2] /= color[3]

        return (color[0], color[1], color[2], color[3])
