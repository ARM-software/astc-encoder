#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# -----------------------------------------------------------------------------
# Copyright 2020-2026 Arm Limited
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
This script is a wrapper around Valgrind that will capture a callgrind profile, and then use
gprof2dot and dot to produce an annotated call graph image. The command line options following the
script options are passed through to astcenc unchanged.

This script requires the following tools available on the PATH:

  * valgrind
  * dot (part of graphviz)

... and the following Python modules:

  * gprof2dot
'''

import argparse
from collections import defaultdict
import os
import re
import subprocess as sp
import sys

from testlib.cli import add_codec_arguments
from testlib.profile import PROFILE_ERRORS, report_error, resolve_binary


def rank_hot_functions(lines: list[str]) -> str:
    '''
    Generate a summary of hot functions in a callgrind_annotate log.

    Args:
        lines: The output of callgrind_annotate.

    Return:
        The report ranking the hot functions.
    '''
    # Parse pattern for useful callgrind lines
    parse_pattern = re.compile(r'^\s*([0-9,]+)\s+\([ 0-9.]+%\)\s+Source/(\S+):(\S+)\(.*\).*$')

    function_cost: dict[str, float] = defaultdict(float)

    for line in lines:
        line = line.strip()
        match = parse_pattern.match(line)

        if not match:
            continue

        cost = float(match.group(1).replace(',', ''))
        source_file = match.group(2)
        function = match.group(3)

        # Filter out third-party library code
        if source_file.startswith('ThirdParty'):
            continue

        # Accumulate the scores from functions in multiple call chains
        function_cost[function] += cost

    # Sum the total cost of the profile
    total_cost = sum(function_cost.values())

    # Sort functions into most expensive first
    function_table = list(function_cost.items())
    function_table.sort(key=lambda x: x[1], reverse=True)

    report: list[str] = []
    total_reported = 0.0

    for func, func_cost in function_table:
        func_percent = (func_cost / total_cost) * 100.0

        # Stop emitting entries when we get to less than 1% contribution
        if func_percent < 1.0:
            break

        total_reported += func_percent
        report.append(f'{func_percent:5.2f}%  {func}')

    report.append('======')
    report.append(f'{total_reported:5.2f}%')

    return '\n'.join(report)


def run_pass(binary: str, arguments: list[str], output: str, show_cli: bool) -> None:
    '''
    Run Valgrind on a single binary.

    Args:
        binary: The astcenc binary path.
        arguments: The astcenc command-line arguments.
        output: The output file base name.
        show_cli: Show command line wrapper tool cost in the profile.

    Raise:
        CalledProcessException: Any subprocess failed.
    '''
    # Run Valgrind
    args = [
        'valgrind',
        '--tool=callgrind',
        '--callgrind-out-file=callgrind.txt',
        binary,
        *arguments
    ]

    sp.run(args, check=True)

    # Run callgrind analysis
    args = ['callgrind_annotate', 'callgrind.txt']

    ret = sp.run(args, stdout=sp.PIPE, check=True, encoding='utf-8')
    lines = ret.stdout.splitlines()

    # Write raw annotated callgrind log
    report = '\n'.join(lines)
    with open(f'{output}_cga.txt', 'w', encoding='utf=8') as handle:
        handle.write(report)

    # Write our text summary of the top functions
    report = rank_hot_functions(lines)
    with open(f'{output}.txt', 'w', encoding='utf=8') as handle:
        handle.write(report)

    # Write the visual call graph
    entry_func = 'main'

    if not show_cli:
        entry_func = 'compress_block(' \
            'astcenc_contexti const&, ' \
            'image_block const&, ' \
            'unsigned char*, ' \
            'compression_working_buffers&)'

    args = [
        'gprof2dot',
        '--format=callgrind',
        '--output=out.dot',
        'callgrind.txt',
        '-s', '-z', entry_func
    ]

    sp.run(args, check=True, text=True)

    args = ['dot', '-Tpng', 'out.dot', '-o', f'{output}.png']

    sp.run(args, check=True, text=True)

    os.remove('out.dot')
    os.remove('callgrind.txt')


def parse_command_line() -> argparse.Namespace:
    '''
    Parse the command line.

    Return:
        The parsed command line.
    '''
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('--show-cli', dest='show_cli', default=False,
                        action='store_true',
                        help='enable CLI wrapper nodes in call graph')

    parser.add_argument('--output', default='perf',
                        help='output file base name (default: %(default)s)')

    add_codec_arguments(parser)

    args = parser.parse_args()

    if not args.arguments:
        parser.error('no astcenc arguments specified')

    return args


def main() -> int:
    '''
    The main function.

    Return:
        int: The process return code.
    '''
    args = parse_command_line()

    try:
        binary = resolve_binary(args.binary)
        run_pass(binary, args.arguments, args.output, args.show_cli)
    except PROFILE_ERRORS as error:
        return report_error(error)

    return 0


if __name__ == '__main__':
    sys.exit(main())
