# SPDX-License-Identifier: Apache-2.0
# -----------------------------------------------------------------------------
# Copyright 2019-2022 Arm Limited
#
# Licensed under the Apache License, Version 2.0 (the 'License'); you may not
# use this file except in compliance with the License. You may obtain a copy
# of the License at:
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an 'AS IS' BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations
# under the License.
# -----------------------------------------------------------------------------
'''
TODO
'''

from collections.abc import Iterable
import os
import re
import subprocess as sp

import testlib.misc as misc


CONVERT_BINARY = ['convert']


class Image():
    '''
    Wrapper around an image on the file system.
    '''

    # TODO: We don't support KTX yet, as ImageMagick doesn't.
    SUPPORTED_LDR = ['bmp', 'jpg', 'png', 'tga']
    SUPPORTED_HDR = ['exr', 'hdr']

    @classmethod
    def is_format_supported(cls, fileFormat, profile=None):
        '''
        Test if a given file format is supported by the library.

        Args:
            fileFormat (str): The file extension (excluding the '.').
            profile (str or None): The profile (ldr or hdr) of the image.

        Returns:
            bool: `True` if the image is supported, `False` otherwise.
        '''
        assert profile in [None, 'ldr', 'hdr']

        if profile == 'ldr':
            return fileFormat in cls.SUPPORTED_LDR

        if profile == 'hdr':
            return fileFormat in cls.SUPPORTED_HDR

        return fileFormat in cls.SUPPORTED_LDR or \
            fileFormat in cls.SUPPORTED_HDR

    def __init__(self, file_path):
        '''
        Construct a new Image.

        Args:
            file_path (str): The path to the image on disk.
        '''
        self.file_path = file_path
        self.proxyPath = None

    def get_colors(self, coords):
        '''
        Get the image colors at the given coordinate.

        Args:
            coords (tuple or list): A single coordinate, or a list of
                coordinates to sample.

        Returns:
            tuple: A single sample color (if `coords` was a coordinate).
            list: A list of sample colors (if `coords` was a list).

            Colors are returned as float values between 0.0 and 1.0 for LDR,
            and float values which may exceed 1.0 for HDR.
        '''
        colors = []

        # We accept both a list of positions and a single position;
        # canonicalize here so the main processing only handles lists
        isList = len(coords) != 0 and isinstance(coords[0], Iterable)

        if not isList:
            coords = [coords]

        for (x, y) in coords:
            command = list(CONVERT_BINARY)
            command += [self.file_path]

            # Ensure convert factors in format origin if needed
            command += ['-auto-orient']

            command += [
                '-format', '%%[pixel:p{%u,%u}]' % (x, y),
                'info:'
            ]

            if os.name == 'nt':
                command.insert(0, 'magick')

            result = sp.run(command, stdout=sp.PIPE, stderr=sp.PIPE,
                            check=True, universal_newlines=True)

            rawcolor = result.stdout.strip()

            # Decode ImageMagick's annoying named color outputs. Note that this
            # only handles 'known' cases triggered by our test images, we don't
            # support the entire ImageMagick named color table.
            if rawcolor == 'black':
                colors.append([0.0, 0.0, 0.0, 1.0])
            elif rawcolor == 'white':
                colors.append([1.0, 1.0, 1.0, 1.0])
            elif rawcolor == 'red':
                colors.append([1.0, 0.0, 0.0, 1.0])
            elif rawcolor == 'blue':
                colors.append([0.0, 0.0, 1.0, 1.0])

            # Decode ImageMagick's format tuples
            elif rawcolor.startswith('srgba'):
                rawcolor = rawcolor[6:]
                rawcolor = rawcolor[:-1]
                channels = rawcolor.split(',')
                for i, channel in enumerate(channels):
                    if (i < 3) and channel.endswith('%'):
                        channels[i] = float(channel[:-1]) / 100.0
                    elif (i < 3) and not channel.endswith('%'):
                        channels[i] = float(channel) / 255.0
                    else:
                        channels[i] = float(channel)
                colors.append(channels)
            elif rawcolor.startswith('srgb'):
                rawcolor = rawcolor[5:]
                rawcolor = rawcolor[:-1]
                channels = rawcolor.split(',')
                for i, channel in enumerate(channels):
                    if (i < 3) and channel.endswith('%'):
                        channels[i] = float(channel[:-1]) / 100.0
                    if (i < 3) and not channel.endswith('%'):
                        channels[i] = float(channel) / 255.0
                channels.append(1.0)
                colors.append(channels)
            elif rawcolor.startswith('rgba'):
                rawcolor = rawcolor[5:]
                rawcolor = rawcolor[:-1]
                channels = rawcolor.split(',')
                for i, channel in enumerate(channels):
                    if (i < 3) and channel.endswith('%'):
                        channels[i] = float(channel[:-1]) / 100.0
                    elif (i < 3) and not channel.endswith('%'):
                        channels[i] = float(channel) / 255.0
                    else:
                        channels[i] = float(channel)
                colors.append(channels)
            elif rawcolor.startswith('rgb'):
                rawcolor = rawcolor[4:]
                rawcolor = rawcolor[:-1]
                channels = rawcolor.split(',')
                for i, channel in enumerate(channels):
                    if (i < 3) and channel.endswith('%'):
                        channels[i] = float(channel[:-1]) / 100.0
                    if (i < 3) and not channel.endswith('%'):
                        channels[i] = float(channel) / 255.0
                channels.append(1.0)
                colors.append(channels)
            else:
                print(x, y)
                print(rawcolor)
                assert False

        # ImageMagick decodes DDS files as BGRA not RGBA; manually correct
        if self.file_path.endswith('dds'):
            for color in colors:
                tmp = color[0]
                color[0] = color[2]
                color[2] = tmp

        # ImageMagick decodes EXR files with premult alpha; manually correct
        if self.file_path.endswith('exr'):
            for color in colors:
                color[0] /= color[3]
                color[1] /= color[3]
                color[2] /= color[3]

        # Undo list canonicalization if we were given a single scalar coord
        if not isList:
            return colors[0]

        return colors
