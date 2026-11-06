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
Shared helpers for the astcenc performance profiling scripts.
'''

import math
from pathlib import Path
import re
import shutil
import subprocess as sp
import sys
from typing import NamedTuple


PROCESSING_TIME_PATTERN = re.compile(
    r'^\s*(?:Coding|Decoding) time:\s*([0-9.]+)\s+s\s*$')

QUALITY_METRICS = {
    'PSNR (LDR-RGB)': ('dB', 'higher'),
    'PSNR (LDR-RGBA)': ('dB', 'higher'),
    'Alpha-weighted PSNR': ('dB', 'higher'),
    'mPSNR (RGB)': ('dB', 'higher'),
    'Mean Angular Error': ('degrees', 'lower'),
    'Worst Angular Error': ('degrees', 'lower')
}

PROFILE_ERRORS = (
    FileNotFoundError,
    OSError,
    sp.CalledProcessError,
    ValueError
)


class CodecResult(NamedTuple):
    '''
    Metrics reported by one codec invocation.

    Attributes:
        time: The processing time in seconds.
        quality: The available image quality scores, keyed by metric name.
    '''
    time: float
    quality: dict[str, float]


def resolve_binary(binary: str) -> str:
    '''
    Resolve a binary path, including commands found on PATH.

    Args:
        binary: The binary path or command name.

    Return:
        The resolved binary path.

    Raises:
        FileNotFoundError: The binary cannot be found.
    '''
    path = Path(binary)
    if path.is_file():
        return str(path)

    resolved = shutil.which(binary)
    if not resolved:
        raise FileNotFoundError(binary)

    return resolved


def run_codec(binary: str, arguments: list[str]) -> CodecResult:
    '''
    Run one astcenc command.

    Args:
        binary: The astcenc binary path.
        arguments: The astcenc command-line arguments.

    Return:
        The reported processing time and available image quality scores.

    Raises:
        CalledProcessError: The codec process failed.
        ValueError: The time is missing or a quality score is invalid.
    '''
    command = [binary] + arguments

    result = sp.run(command, stdout=sp.PIPE, stderr=sp.PIPE,
                    check=True, text=True)

    processing_time = None
    quality = {}
    for line in result.stdout.splitlines():
        if processing_time is None:
            if match := PROCESSING_TIME_PATTERN.match(line):
                processing_time = float(match.group(1))

        name, _, value = line.strip().partition(':')
        # Older encoders can include an f-stop range in the mPSNR label.
        name = re.sub(r'\s+\[.*?\]\s*$', '', name)
        if name in QUALITY_METRICS:
            if not value.strip():
                raise ValueError(f'Codec reported an invalid {name}')
            score = float(value.split()[0])
            if not math.isfinite(score):
                raise ValueError(f'Codec reported a non-finite {name}')
            quality[name] = score

    if processing_time is not None:
        return CodecResult(processing_time, quality)

    raise ValueError(f'{binary} did not report a processing time')


def collect_samples(binary_a: str, binary_b: str, arguments: list[str],
                    warmups: int,
                    runs: int) -> tuple[list[CodecResult], list[CodecResult]]:
    '''
    Warm up the codecs and collect interleaved paired timing samples.

    Args:
        binary_a: The binary A path.
        binary_b: The binary B path.
        arguments: The astcenc command-line arguments.
        warmups: The number of warm-up pairs.
        runs: The number of measured pairs.

    Return:
        The timing and quality samples for binary A and binary B.
    '''
    for _ in range(warmups):
        run_codec(binary_a, arguments)
        run_codec(binary_b, arguments)

    samples_a = []
    samples_b = []
    for run in range(runs):
        if run % 2 == 0:
            result_a = run_codec(binary_a, arguments)
            result_b = run_codec(binary_b, arguments)
        else:
            result_b = run_codec(binary_b, arguments)
            result_a = run_codec(binary_a, arguments)

        samples_a.append(result_a)
        samples_b.append(result_b)
        print(f'Pair {run + 1:>2}/{runs}: '
              f'A {result_a.time:.4f} s, '
              f'B {result_b.time:.4f} s')

    return samples_a, samples_b


def profile_binaries(
    binary_a: str,
    binary_b: str,
    arguments: list[str],
    warmups: int,
    runs: int
) -> tuple[str, str, list[CodecResult], list[CodecResult]]:
    '''
    Resolve the binaries and collect timing samples.

    Args:
        binary_a: The binary A path or command name.
        binary_b: The binary B path or command name.
        arguments: The astcenc command-line arguments.
        warmups: The number of warm-up pairs.
        runs: The number of measured pairs.

    Return:
        The two binary paths followed by their timing and quality samples.
    '''
    resolved_a = resolve_binary(binary_a)
    resolved_b = resolve_binary(binary_b)

    samples_a, samples_b = collect_samples(
        resolved_a, resolved_b, arguments, warmups, runs)

    return resolved_a, resolved_b, samples_a, samples_b


def report_error(error: Exception) -> int:
    '''
    Report an expected profiling error.

    Args:
        error: The profiling error.

    Return:
        The process return code.
    '''
    output = ''
    if isinstance(error, FileNotFoundError):
        message = f'Binary not found: {error}'
    elif isinstance(error, sp.CalledProcessError):
        message = f'Codec failed: {" ".join(error.cmd)}'
        output = error.stderr.strip() or error.stdout.strip()
    elif isinstance(error, OSError):
        message = f'Could not execute codec: {error}'
    else:
        message = str(error)

    print(f'ERROR: {message}', file=sys.stderr)
    if output:
        print(output, file=sys.stderr)

    return 1
