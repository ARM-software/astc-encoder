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
Download test image sets, if images are missing already.
"""
import os
import sys
import urllib.request

def retrieve_kodak_set():
    for i in range(1, 25):
        fle = "ldr-rgb-kodim%02u.png" % i
        dst = os.path.join("Test", "Kodak_Images", "LDR-RGB", fle)
        src = "http://r0k.us/graphics/kodak/kodak/kodim%02u.png" % i

        if not os.path.exists(dst):
            print("Kodak image %u: Downloading" % i)
            urllib.request.urlretrieve(src, dst)
        else:
            print("Kodak image %u: Skipping" % i)

def main():
    """
    The main function.
    """
    retrieve_kodak_set()
    return 0

if __name__ == "__main__":
    sys.exit(main())
