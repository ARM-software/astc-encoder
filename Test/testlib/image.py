# SPDX-License-Identifier: Apache-2.0
# -----------------------------------------------------------------------------
# Copyright 2019-2020 Arm Limited
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
"""
This module contains code for loading image metadata from a file path on disk.

The directory path is structured:

    TestSetName/TestFormat/FileName

... and the file name is structured:

    colorProfile-colorFormat-name[-flags].extension
"""

import os

import testlib.misc as misc


class ImageException(Exception):
    """
    Exception thrown for bad image specification.
    """


class TestImage():
    """
    Objects of this type contain metadata for a single test image on disk.

    Attributes:
        filePath: The path of the file on disk.
        outFilePath: The path of the output file on disk.
        testSet: The name of the test set.
        testFormat: The test format group.
        testFile: The test file name.
        colorProfile: The image compression color profile.
        colorFormat: The image color format.
        name: The image human name.
        is3D: True if the image is 3D, else False.
        isMask: True if the image is a non-correlated mask texture, else False.
        isAlphaScaled: True if the image wants alpha scaling, else False.
        TEST_EXTS: Expected test image extensions.
        PROFILES: Tuple of valid color profile values.
        FORMATS: Tuple of valid color format values.
        FLAGS: Map of valid flags (key) and their meaning (value).
    """
    TEST_EXTS = (".png", ".dds", ".hdr")

    PROFILES = ("ldr", "ldrs", "hdr")

    FORMATS = ("l", "la", "xy", "rgb", "rgba")

    FLAGS = {
        # Flags for image compression control
        "3": "3D image",
        "m": "Mask image",
        "a": "Alpha scaled image"
    }

    def __init__(self, filePath):
        """
        Create a new image definition, based on a structured file path.

        Args:
            filePath: The path of the image on disk.
        """
        self.filePath = os.path.abspath(filePath)
        if not os.path.exists(self.filePath):
            raise ImageException("Image doesn't exist (%s)" % filePath)

        # Decode the path
        scriptDir = os.path.dirname(__file__)
        rootInDir = os.path.join(scriptDir, "..", "Images")
        partialPath = os.path.relpath(self.filePath, rootInDir)
        parts = misc.path_splitall(partialPath)
        if len(parts) != 3:
            raise ImageException("Image path not path triplet (%s)" % parts)
        self.testSet = parts[0]
        self.testFormat = parts[1]
        self.testFile = parts[2]

        # Decode the file name
        self.decode_file_name(self.testFile)

        # Output file path (store base without extension)
        rootOutDir = os.path.join(scriptDir, "..", "..", "TestOutput")
        outFilePath = os.path.join(rootOutDir, partialPath)
        outFilePath = os.path.abspath(outFilePath)
        outFilePath = os.path.splitext(outFilePath)[0]
        self.outFilePath = outFilePath

    def decode_file_name(self, fileName):
        """
        Utility function to decode metadata from an encoded file name.

        Args:
            fileName: The file name to tokenize.
        """
        # Strip off the extension
        rootName = os.path.splitext(fileName)[0]

        parts = rootName.split("-")

        # Decode the mandatory fields
        if len(parts) >= 3:
            self.colorProfile = parts[0]
            if self.colorProfile not in self.PROFILES:
                raise ImageException("Unknown color profile (%s)" % parts[0])

            self.colorFormat = parts[1]
            if self.colorFormat not in self.FORMATS:
                raise ImageException("Unknown color format (%s)" % parts[1])

            # Consistency check between directory and file names
            reencode = "%s-%s" % (self.colorProfile, self.colorFormat)
            compare = self.testFormat.lower()
            if reencode != compare:
                dat = (self.testFormat, reencode)
                raise ImageException("Mismatched test and image (%s:%s)" % dat)

            self.name = parts[2]

        # Set default values for the optional fields
        self.is3D = False
        self.isMask = False
        self.isAlphaScaled = False

        # Decode the flags field if present
        if len(parts) >= 4:
            flags = parts[3]
            seenFlags = set()
            for flag in flags:
                if flag in seenFlags:
                    raise ImageException("Duplicate flag (%s)" % flag)
                if flag not in self.FLAGS:
                    raise ImageException("Unknown flag (%s)" % flag)
                seenFlags.add(flag)

            self.is3D = "3" in seenFlags
            self.isMask = "m" in seenFlags
            self.isAlphaScaled = "a" in seenFlags
