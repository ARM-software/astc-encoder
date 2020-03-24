#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# -----------------------------------------------------------------------------
# Copyright 2020 Arm Limited
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
Utility to show info about an image.
"""

import argparse
from PIL import Image
import sys


def main_pipette(args):
    """
    Main function for the "pipette" mode.

    This mode prints the color at a pixel coordinate in the image, in a variety
    of color format (decimal, HTML string, float).

    Returns:
        The process return code.
    """
    retCode = 0

    for i, image in enumerate(args.images):
        if i != 0:
            print("")

        img = Image.open(image.name)
        x = args.location[0]
        y = args.location[1]

        print(image.name)
        print("=" * len(image.name))

        if (x >= img.size[0]) or (y >= img.size[1]):
            print("- ERROR: location out-of-bounds [%ux%u]" % img.size)
            retCode = 1
        else:
            color = img.getpixel((x, y))
            # Print byte values
            print("+ Byte: %s" % str(color))

            # Print hex values
            fmtString = "+ Hex: #" + ("%02X" * len(color))
            print(fmtString % color)

            # Print float values
            parts = ["%g"] * len(color)
            parts = ", ".join(parts)
            fmtString = "+ Float: (" + parts + ")"
            print(fmtString % tuple(float(x)/255.0 for x in color))

    return retCode


def main_info(args):
    """
    Main function for the "info" mode.

    This mode prints the basic metadata of an image:

        - the overall image size.
        - the number of color channels.
        - the min/max value in each color channel.

    Returns:
        The process return code.
    """
    for i, image in enumerate(args.images):
        if i != 0:
            print("")

        img = Image.open(image.name)
        minmax = img.getextrema()

        print(image.name)
        print("=" * len(image.name))
        print("+ Size: %ux%u" % (img.size[0], img.size[1]))
        print("+ Channels: %s" % ("".join(img.getbands())))
        for j, channel in enumerate(img.getbands()):
            print("  + %s: %u - %u" % (channel, *minmax[j]))

    return 0


def parse_loc(value):
    """
    Command line argument parser for position arguments.

    Args:
        value: The command line argument string to parse. Must be of the form
            "<int>x<int>", where both integers must be zero or positive.

    Returns:
        The parsed value [int, int].

    Raises:
        ArgumentTypeError: The value is not a valid location.
    """
    error = argparse.ArgumentTypeError("%s is an invalid location" % value)
    svalue = value.split("x")

    if len(svalue) != 2:
        raise error

    try:
        ivalue = [int(x) for x in svalue if int(x) >= 0]
    except ValueError:
        raise error

    if len(ivalue) != len(svalue):
        raise error

    return ivalue


def parse_command_line():
    """
    Parse the command line.

    Returns:
        The parsed command line container.
    """
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(
        title="Operations")

    # Create the parser for the "pipette" command
    parser_a = subparsers.add_parser(
        "pipette",
        help="Print color at given coordinate")

    parser_a.set_defaults(func=main_pipette)

    parser_a.add_argument(
        "location", metavar="loc", type=parse_loc,
        help="The location spec XxY")

    parser_a.add_argument(
        "images", metavar="image", nargs="+", type=argparse.FileType("r"),
        help="The images to query")

    # Create the parser for the "size" command
    parser_b = subparsers.add_parser(
        "info",
        help="Print image metadata info")

    parser_b.set_defaults(func=main_info)

    parser_b.add_argument(
        "images", metavar="image", nargs="+", type=argparse.FileType("r"),
        help="The images to query")

    # Cope with the user failing to specify any sub-command. Note on Python 3.8
    # we could use required=True on the add_subparsers call, but we cannot do
    # this on 3.6 which is our current min-spec.
    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        return None

    return args


def main():
    """
    The main function.
    """
    args = parse_command_line()
    if args:
        return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
