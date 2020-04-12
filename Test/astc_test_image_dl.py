#!/usr/bin/env python3
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
The ``astc_test_image_dl`` utility provides a means to programatically download
test images that are available online, avoiding the need to duplicate them in
the git repository.
"""


import os
import sys
import urllib.request


TEST_IMAGE_DIR = os.path.join("Test", "Images")


def download(testSet, index, srcUrl, dstPath):
    """
    Download a single image.

    Args:
        testSet (str): The test set name.
        index (int): The download index.
        srcUrl (str): The download URL.
        dstPath (str): The destination path.
    """
    dirName = os.path.dirname(dstPath)
    if not os.path.exists(dirName):
        os.makedirs(dirName)

    # Skip downloads if the file already exists
    if not os.path.exists(dstPath):
        print("%s image %u: Downloading" % (testSet, index))
        urllib.request.urlretrieve(srcUrl, dstPath)
    else:
        print("%s image %u: Skipping" % (testSet, index))


def retrieve_kodak_set():
    """
    Download the public domain Kodak image set.
    """
    testSet = "Kodak"
    for i in range(1, 25):
        fle = "ldr-rgb-kodak%02u.png" % i
        dst = os.path.join(TEST_IMAGE_DIR, "Kodak_Images", "LDR-RGB", fle)
        src = "http://r0k.us/graphics/kodak/kodak/kodim%02u.png" % i
        download(testSet, i, src, dst)


def main():
    """
    The main function.

    Returns:
        int: The process return code.
    """
    retrieve_kodak_set()
    return 0


if __name__ == "__main__":
    sys.exit(main())
