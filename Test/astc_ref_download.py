#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# This confidential and proprietary software may be used only as authorised by
# a licensing agreement from Arm Limited.
#     (C) COPYRIGHT 2019-2020 Arm Limited, ALL RIGHTS RESERVED
# The entire notice above must be reproduced on all authorised copies and
# copies may only be made to the extent permitted by a licensing agreement from
# Arm Limited.
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
