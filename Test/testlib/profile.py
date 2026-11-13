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

from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
import hashlib
import math
from pathlib import Path
import re
import shutil
import subprocess as sp
import sys

from .timing import paired_t_test, relative_ci_half_width


PROCESSING_TIME_PATTERN = re.compile(r'^\s*(?:Coding|Decoding) time:\s*([0-9.]+)\s+s\s*$')

QUALITY_METRICS = {
    'PSNR (LDR-RGB)': ('dB', 'higher'),
    'PSNR (LDR-RGBA)': ('dB', 'higher'),
    'Alpha-weighted PSNR': ('dB', 'higher'),
    'mPSNR (RGB)': ('dB', 'higher'),
    'SSIM (LDR-RGB)': ('', 'higher'),
    'SSIM (LDR-RGBA)': ('', 'higher'),
    'mSSIM (RGB)': ('', 'higher'),
    'Mean Angular Error': ('degrees', 'lower'),
    'Worst Angular Error': ('degrees', 'lower'),
}

PROFILE_ERRORS = (FileNotFoundError, OSError, sp.CalledProcessError, ValueError)


@dataclass(frozen=True)
class CodecResult:
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
    path = Path(binary).expanduser()

    if path.is_file():
        return str(path.resolve())

    resolved = shutil.which(binary)

    if not resolved:
        raise FileNotFoundError(binary)

    return str(Path(resolved).resolve())


def parse_codec_output(output: str) -> CodecResult:
    '''
    Parse processing time and available quality metrics from codec output.

    Args:
        output: The codec standard output.

    Return:
        The processing time and named image quality scores.

    Raises:
        ValueError: A required timing value is missing or a reported metric is invalid.
    '''
    processing_time = None
    quality = {}

    # Retain named metrics so each caller can choose its channel and profile policy
    for line in output.splitlines():
        if processing_time is None:
            if match := PROCESSING_TIME_PATTERN.match(line):
                processing_time = float(match.group(1))

        name, _, value = line.strip().partition(':')
        # Older encoders can include an f-stop range in the mPSNR label
        name = re.sub(r'\s+\[.*?\]\s*$', '', name)

        if name in QUALITY_METRICS:
            if not value.strip():
                raise ValueError(f'Codec reported an invalid {name}')

            score = float(value.split()[0])

            if not math.isfinite(score):
                raise ValueError(f'Codec reported a non-finite {name}')

            quality[name] = score

    # Reject unusable timings before they reach statistics or throughput calculations
    if processing_time is None:
        raise ValueError('Codec did not report a processing time')

    if (not math.isfinite(processing_time)) or (processing_time <= 0.0):
        raise ValueError(f'Codec reported an invalid processing time: {processing_time}')

    return CodecResult(time=processing_time, quality=quality)


def run_codec(binary: str, arguments: list[str]) -> CodecResult:
    '''
    Run one codec invocation and parse its metrics.

    Args:
        binary: The astcenc binary path.
        arguments: The codec command-line arguments.

    Return:
        The processing time and available image quality scores.

    Raises:
        CalledProcessError: The codec process failed.
        ValueError: A reported metric is missing or invalid.
    '''
    result = sp.run([binary, *arguments], stdout=sp.PIPE, stderr=sp.PIPE, check=True, text=True)

    return parse_codec_output(result.stdout)


def fingerprint(path: Path) -> str:
    '''
    Hash a binary or image to detect changes during a measurement.

    Args:
        path: The file to hash.

    Return:
        The hexadecimal SHA-256 digest.
    '''

    return hashlib.sha256(path.read_bytes()).hexdigest()


def binary_metadata(binary: str) -> dict[str, str]:
    '''
    Resolve and identify a codec binary for a report.

    Args:
        binary: The executable path or command name.

    Return:
        The absolute path, SHA-256 digest, and codec version output.

    Raises:
        FileNotFoundError: The binary cannot be found.
        CalledProcessError: The version query failed.
    '''
    path = resolve_binary(binary)
    version = sp.run([path, '-version'], check=True, capture_output=True, text=True)

    return {
        'path': path,
        'sha256': fingerprint(Path(path)),
        'version': version.stdout.strip(),
    }


def select_quality_metrics(arguments: list[str], available: dict[str, float]) -> list[str]:
    '''
    Choose quality metrics for the input image and compression settings.

    Args:
        arguments: The codec command-line arguments.
        available: The quality scores reported by the codec.

    Return:
        The relevant metric names, or an empty list outside test mode.
    '''

    if arguments[0] not in ('-tl', '-ts', '-th', '-tH'):
        return []

    if '-normal' in arguments:
        return ['PSNR (LDR-RGB)', 'Mean Angular Error', 'Worst Angular Error']

    image = Path(arguments[1])
    parts = image.stem.lower().split('-')

    if ('mPSNR (RGB)' in available or
            image.suffix.lower() in ('.hdr', '.exr') or
            (len(parts) >= 3 and parts[0] == 'hdr')):
        return ['mPSNR (RGB)']

    if '-a' in arguments:
        return ['Alpha-weighted PSNR']

    if len(parts) >= 3 and parts[0] in ('ldr', 'ldrs', 'hdr'):
        if parts[1] in ('rgb', 'l', 'xy', 'rgba', 'la'):
            channels = 'RGB'

            if parts[1] in ('rgba', 'la'):
                channels = 'RGBA'

            return [f'PSNR (LDR-{channels})']

    # Generic filenames do not encode which channels are relevant. Show both instead of silently
    # including or ignoring alpha in a single score.
    return [name for name in ('PSNR (LDR-RGB)', 'PSNR (LDR-RGBA)') if name in available]


def collect_samples(binary_a: str, binary_b: str, arguments: list[str],
                    warmups: int,
                    runs: int,
                    confidence: float | None = None
                    ) -> tuple[list[CodecResult], list[CodecResult]]:
    '''
    Warm up the codecs and collect interleaved paired timing samples.

    Args:
        binary_a: The binary A path.
        binary_b: The binary B path.
        arguments: The astcenc command-line arguments.
        warmups: The number of warm-up pairs.
        runs: The maximum number of measured pairs.
        confidence: Optional confidence level for early stopping.

    Return:
        The timing and quality samples for binary A and binary B.
    '''

    # Warmups use the same codec settings but never enter the measured sample sets
    for _ in range(warmups):
        run_codec(binary_a, arguments)
        run_codec(binary_b, arguments)

    samples_a = []
    samples_b = []

    # Alternate the leading binary to reduce systematic ordering effects
    for run in range(runs):
        if run % 2 == 0:
            result_a = run_codec(binary_a, arguments)
            result_b = run_codec(binary_b, arguments)
        else:
            result_b = run_codec(binary_b, arguments)
            result_a = run_codec(binary_a, arguments)

        samples_a.append(result_a)
        samples_b.append(result_b)
        print(f'Pair {run + 1:>2}/{runs}: A {result_a.time:.4f} s, B {result_b.time:.4f} s')

        if confidence is not None and len(samples_a) >= 2:
            times_a = [sample.time for sample in samples_a]
            times_b = [sample.time for sample in samples_b]
            _, p = paired_t_test(times_a, times_b)

            if p < 1.0 - confidence:
                print(f'Stopping after {len(samples_a)} pairs: '
                      f'{confidence:.1%} nominal confidence reached (p={p:.6g})')
                break

    return samples_a, samples_b


def profile_binaries(
    binary_a: str,
    binary_b: str,
    arguments: list[str],
    warmups: int,
    runs: int,
    confidence: float | None = None
) -> tuple[str, str, list[CodecResult], list[CodecResult]]:
    '''
    Resolve the binaries and collect timing samples.

    Args:
        binary_a: The binary A path or command name.
        binary_b: The binary B path or command name.
        arguments: The astcenc command-line arguments.
        warmups: The number of warm-up pairs.
        runs: The maximum number of measured pairs.
        confidence: Optional confidence level for early stopping.

    Return:
        The two binary paths followed by their timing and quality samples.
    '''
    resolved_a = resolve_binary(binary_a)
    resolved_b = resolve_binary(binary_b)

    samples_a, samples_b = collect_samples(
        resolved_a, resolved_b, arguments, warmups, runs, confidence)

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
        output = (error.stderr or '').strip() or (error.stdout or '').strip()
    elif isinstance(error, OSError):
        message = f'Could not execute codec: {error}'
    else:
        message = str(error)

    print(f'ERROR: {message}', file=sys.stderr)

    if output:
        print(output, file=sys.stderr)

    return 1


@dataclass
class SampleRound:
    '''
    One serial invocation of every binary in a specified order.

    Attributes:
        order: The binary roles in execution order.
        results: The parsed results keyed by binary role.
    '''
    order: list[str]
    results: dict[str, CodecResult] = field(default_factory=dict)


@dataclass
class PrecisionSamples:
    '''
    Warmups and measured rounds for independent per-binary timing precision checks.

    Attributes:
        arguments: The codec arguments shared by every binary.
        warmups: The excluded warmup rounds.
        rounds: The measured rounds in execution order.
        relative_ci_half_widths: The final relative confidence interval half-widths.
        stop_reason: Whether the run limit or precision target ended collection.
    '''
    arguments: list[str]
    warmups: list[SampleRound] = field(default_factory=list)
    rounds: list[SampleRound] = field(default_factory=list)
    relative_ci_half_widths: dict[str, float] = field(default_factory=dict)
    stop_reason: str = 'run_limit'

    def update_intervals(self, names: Sequence[str], confidence: float, precision: float) -> bool:
        '''
        Update confidence intervals and check the precision target for every binary.

        Args:
            names: The binary roles in report order.
            confidence: The confidence level as a fraction in (0, 1).
            precision: The allowed interval half-width as a fraction of the mean.

        Return:
            True when at least two measured rounds meet the target for every binary.
        '''

        if len(self.rounds) < 2:
            return False

        # Each build must meet the target using only its own measured samples
        self.relative_ci_half_widths = {
            name: relative_ci_half_width([r.results[name].time for r in self.rounds], confidence)
            for name in names
        }

        return all(width <= precision for width in self.relative_ci_half_widths.values())


def collect_precision_samples(binaries: dict[str, str], arguments: list[str], runs: int,
                              warmups: int, confidence: float = 0.95, precision: float = 0.01,
                              execute: Callable[[str, list[str]], CodecResult] | None = None
                              ) -> PrecisionSamples:
    '''
    Collect serial, rotating rounds until every binary meets its timing precision target.

    Args:
        binaries: The codec paths keyed by binary role, in report order.
        arguments: The shared codec command-line arguments.
        runs: The maximum number of measured rounds.
        warmups: The number of excluded warmup rounds.
        confidence: The confidence level as a fraction in (0, 1).
        precision: The allowed interval half-width as a fraction of the mean.
        execute: The invocation function, or None to use the shared codec runner.

    Return:
        The raw samples, final intervals, and reason for stopping.

    Raises:
        ValueError: The binary set is empty or the run counts are invalid.
    '''

    if not binaries or (runs < 1) or (warmups < 0):
        raise ValueError('Precision sampling requires binaries and valid run counts')

    if execute is None:
        execute = run_codec

    # Rotate the leading binary, then reverse the rotations to balance execution order
    names = list(binaries)
    orders = [names[i:] + names[:i] for i in range(len(names))]
    orders += [list(reversed(order)) for order in orders]
    samples = PrecisionSamples(arguments=arguments)

    for index in range(warmups + runs):
        sample = SampleRound(order=orders[index % len(orders)])

        for name in sample.order:
            sample.results[name] = execute(binaries[name], arguments)

        if index < warmups:
            samples.warmups.append(sample)
            print(f'  Warmup {index + 1}/{warmups}', flush=True)
        else:
            samples.rounds.append(sample)

            timings = ' '.join(f'{name}={sample.results[name].time:.4f}s' for name in names)
            print(f'  Round {index - warmups + 1}/{runs}: {timings}', flush=True)

            if samples.update_intervals(names, confidence, precision):
                samples.stop_reason = 'precision'
                count = len(samples.rounds)
                print(f'  Stopping after {count} rounds: '
                      f'all {confidence:.0%} timing CIs within '
                      f'+/-{precision:.0%} of their means', flush=True)
                break

    return samples
