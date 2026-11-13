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
These classes provide an abstraction around the astcenc command line tool, allowing the rest of the
image test suite to ignore changes in the command line interface.

Currently this module supports the latest 2.x branch onwards; the 1.x branch is no longer supported.

The base class handles execution and metric aggregation. Concrete wrappers select binary paths
and metric patterns for the supported versions.
'''

from abc import ABC, abstractmethod
from dataclasses import dataclass
import math
import os
from pathlib import Path
import re
import subprocess as sp
import sys
from typing import Optional

from .testset import TestImage
from .timing import is_timing_stable


@dataclass(frozen=True)
class RunOptions:
    '''
    Settings shared by image compression runs.

    Attributes:
        repeats: The maximum number of measured runs, at least one.
        keep_output: Whether to preserve decompressed output images.
        threads: The codec thread count, or None for automatic selection.
        confidence: The timing confidence level as a fraction in (0, 1).
        precision: The allowed interval half-width as a fraction of the mean.
    '''
    repeats: int = 1
    keep_output: bool = False
    threads: int | None = None
    confidence: float = 0.95
    precision: float = 0.01

    def __post_init__(self) -> None:
        '''
        Validate the repeat limit and timing thresholds.

        Raises:
            ValueError: A run count or timing threshold is outside its valid range.
        '''

        if self.repeats < 1:
            raise ValueError('The repeat limit must be at least one')

        if not 0.0 < self.confidence < 1.0:
            raise ValueError('Confidence must be between zero and one, exclusive')

        if (not math.isfinite(self.precision)) or (self.precision <= 0.0):
            raise ValueError('Precision must be a finite positive fraction')


@dataclass(frozen=True)
class RunResult:
    '''
    Image quality and performance metrics from encoder runs.

    Attributes:
        psnr: The PSNR or mPSNR score in dB.
        ssim: The SSIM or mSSIM score, or zero if not reported.
        total_time: The total elapsed time in seconds.
        coding_time: The compression time in seconds.
        coding_rate: The compression throughput in MT/s.
    '''
    psnr: float
    ssim: float
    total_time: float
    coding_time: float
    coding_rate: float

    def best_of(self, other: 'RunResult') -> 'RunResult':
        '''
        Keep the best value of each metric across two results.

        Args:
            other: The result to merge with this result.

        Return:
            The highest quality scores and rate, and the lowest elapsed times.
        '''

        return RunResult(
            psnr=max(self.psnr, other.psnr),
            ssim=max(self.ssim, other.ssim),
            total_time=min(self.total_time, other.total_time),
            coding_time=min(self.coding_time, other.coding_time),
            coding_rate=max(self.coding_rate, other.coding_rate))


class EncoderBase(ABC):
    '''
    A wrapper around the astcenc binary, abstracting the command line so that the rest of the test
    suite does not need to worry about changes to the command line interface over versions.

    This is an abstract base class providing some generic helper functionality used by concrete
    instantiations of subclasses.

    Attributes:
        name: The encoder name to use in reports.
        binary_path: The encoder binary path.
        last_run_count: The number of samples collected by the latest test.

    Class attributes:
        version: Encoder version or branch.
        remap_color_switches: Dict of color format to profile switches.
        remap_output_types: Dict of color format to output file extensions.
    '''

    # Subclasses override these if they need them
    version = ''
    remap_color_switches: dict[str, str] = {}
    remap_output_types: dict[str, str] = {}

    def __init__(self, name: str, binary_path: Path):
        '''
        Create a new encoder instance.

        Args:
            name: The name of the encoder.
            binary_path: The path to the binary on the file system.
        '''
        self.name = name
        self.binary_path = binary_path
        self.last_run_count = 0

    @abstractmethod
    def build_cli(self, image: TestImage, block_size: str = '6x6',
                  preset: str = '-thorough', keep_output: bool = True,
                  threads: Optional[int] = None) -> list[str]:
        '''
        Build the command line needed for the given test.

        Args:
            image: The test image to compress.
            block_size: The block size to use.
            preset: The quality preset to use.
            keep_output: Should the test preserve output images? This is a hint and may be ignored
                if astcenc version used can't do it.
            threads: The thread count to use.

        Return:
            The command line arguments for this encoder version.
        '''

    def execute(self, command: list[str]) -> list[str]:
        '''
        Run a subprocess with the specified command.

        Args:
            command: The command line arguments to use.

        Return:
            The output log (stdout) split into lines.
        '''
        try:
            result = sp.run(command, stdout=sp.PIPE, stderr=sp.PIPE, check=True, text=True)

        except (OSError, sp.CalledProcessError):
            print(f'  + {" ".join(command)}')
            assert False, 'ERROR: Test run failed'

        return result.stdout.splitlines()

    def parse_output(self, image: TestImage, output: list[str]) -> RunResult:
        '''
        Parse the log output for quality and performance metrics.

        Args:
            image: The test image that was compressed.
            output: The astcenc compression output log.

        Return:
            The quality scores, elapsed times, and coding throughput.
        '''
        # Regex patterns. provided by this particular subclass
        pattern_psnr = self.get_psnr_pattern(image)
        pattern_ssim = self.get_ssim_pattern(image)
        pattern_total_time = self.get_total_time_pattern()
        pattern_coding_time = self.get_coding_time_pattern()
        pattern_coding_rate = self.get_coding_rate_pattern()

        # Extract results from the log
        psnr = None
        ssim = 0.0
        total_time = None
        coding_time = None
        coding_rate = None

        for line in output:
            if match := pattern_psnr.match(line):
                psnr = float(match.group(1))
                continue

            if match := pattern_ssim.match(line):
                ssim = float(match.group(1))
                continue

            if match := pattern_total_time.match(line):
                total_time = float(match.group(1))
                continue

            if match := pattern_coding_time.match(line):
                coding_time = float(match.group(1))
                continue

            if match := pattern_coding_rate.match(line):
                coding_rate = float(match.group(1))
                continue

        stdout = '\n'.join(output)
        assert psnr is not None, f'Missing PSNR {stdout}'
        assert total_time is not None, f'Missing total time {stdout}'
        assert coding_time is not None, f'Missing coding time {stdout}'
        assert coding_rate is not None, f'Missing coding rate {stdout}'

        return RunResult(psnr=psnr, ssim=ssim, total_time=total_time,
                         coding_time=coding_time, coding_rate=coding_rate)

    @abstractmethod
    def get_psnr_pattern(self, image: TestImage) -> re.Pattern:
        '''
        Get the regex pattern to match the image quality metric.

        Note that, although this function is called PSNR, for some images we may choose to match
        another metric (e.g. mPSNR for HDR images).

        Args:
            image: The test image that was compressed.

        Return:
            The regex pattern.
        '''

    @abstractmethod
    def get_ssim_pattern(self, image: TestImage) -> re.Pattern:
        '''
        Get the regex pattern to match SSIM or mSSIM for the image.

        Older encoder versions may omit this metric.

        Args:
            image: The test image that was compressed.

        Return:
            The regex pattern.
        '''

    @abstractmethod
    def get_total_time_pattern(self) -> re.Pattern:
        '''
        Get the regex pattern to match the total compression time.

        Return:
            The regex pattern.
        '''

    @abstractmethod
    def get_coding_time_pattern(self) -> re.Pattern:
        '''
        Get the regex pattern to match the coding compression time.

        Return:
            The regex pattern.
        '''

    @abstractmethod
    def get_coding_rate_pattern(self) -> re.Pattern:
        '''
        Get the regex pattern to match the coding rate.

        Return:
            The regex pattern.
        '''

    def run_test(self, image: TestImage, block_size: str, preset: str,
                 options: RunOptions) -> RunResult:
        '''
        Run until the timing precision target or repeat limit is reached.

        Args:
            image: The test image to compress.
            block_size: The block size to use.
            preset: The quality-performance preset to use.
            options: The repeat limit, codec settings, and timing thresholds.

        Return:
            The best value of each quality and performance metric across measured runs.
        '''
        # Keep codec arguments fixed across repeated measurements
        command = self.build_cli(image, block_size, preset, options.keep_output, options.threads)

        # Accumulate each metric independently; SSIM can be negative
        best = RunResult(psnr=0.0, ssim=-1.0, total_time=sys.float_info.max,
                         coding_time=sys.float_info.max, coding_rate=0.0)
        times = []
        self.last_run_count = 0

        for _ in range(options.repeats):
            output = self.execute(command)
            result = self.parse_output(image, output)
            best = best.best_of(result)

            # Check this build's timing precision independently of the CSV
            times.append(result.coding_time)
            self.last_run_count = len(times)

            if is_timing_stable(times, options.confidence, options.precision):
                break

        return best


class Encoder2x(EncoderBase):
    '''
    This class wraps the latest astcenc interface, supported from the 2.0 release or later, which
    broke argument compatibility with the earlier 1.x series.
    '''
    version = 'main'

    remap_color_switches = {
        'ldr': '-tl',
        'ldrs': '-ts',
        'hdr': '-th',
        'hdra': '-tH',
    }

    remap_output_types = {
        'ldr': '.png',
        'ldrs': '.png',
        'hdr': '.exr',
        'hdra': '.exr',
    }

    def __init__(self, variant: str, binary_path: Optional[Path] = None):
        name = f'astcenc-{variant}-{self.version}'

        if binary_path is None:
            binary = 'astcenc'
            post = f'-{variant}' if variant != 'universal' else ''
            ext = '.exe' if os.name == 'nt' else ''
            binary_path = Path('./') / 'bin' / f'{binary}{post}{ext}'

        super().__init__(name, binary_path)

    def build_cli(self, image, block_size='6x6', preset='-thorough',
                  keep_output=True, threads=None):
        operation = self.remap_color_switches[image.color_profile]

        src_path = image.file_path

        if keep_output:
            dst_path = str(self.get_output_path(image, block_size, preset))

        else:
            if sys.platform == 'win32':
                dst_path = 'nul'
            else:
                dst_path = '/dev/null'

        command = [
            str(self.binary_path), operation, str(src_path), str(dst_path),
            block_size, preset, '-silent'
        ]

        if image.color_format == 'xy':
            command.append('-normal')

        if image.is_alpha_scaled:
            command.append('-a')
            command.append('1')

        if threads is not None:
            command.append('-j')
            command.append(f'{threads}')

        return command

    def get_output_path(self, image: TestImage, block_size: str, preset: str) -> Path:
        '''
        Choose an output path and create its parent directory for retained test images.

        Args:
            image: The source image metadata.
            block_size: The ASTC block size.
            preset: The quality preset including its leading dash.

        Return:
            The output path grouped by image set, format, encoder, preset, and block size.
        '''
        root = Path(__file__).parent.parent / 'TestOutput'
        base = root / image.test_set / image.test_format / image.file_path.stem
        directory = base.parent / self.name / preset[1:] / block_size
        directory.mkdir(parents=True, exist_ok=True)
        extension = self.remap_output_types[image.color_profile]

        return directory / f'{base.stem}{extension}'

    def get_psnr_pattern(self, image: TestImage) -> re.Pattern:
        # HDR profile
        if image.color_profile == 'hdr':
            pattern_psnr = r'\s*mPSNR \(RGB\)(?: \[.*?\] )?:\s*([0-9.]*) dB.*'
        # LDR profile color formats
        elif image.color_format == 'rgba':
            pattern_psnr = r'\s*PSNR \(LDR-RGBA\):\s*([0-9.]*) dB'
        else:
            pattern_psnr = r'\s*PSNR \(LDR-RGB\):\s*([0-9.]*) dB'

        return re.compile(pattern_psnr)

    def get_ssim_pattern(self, image: TestImage) -> re.Pattern:
        '''
        Get the SSIM pattern for the image profile and color format.

        Args:
            image: The test image that was compressed.

        Return:
            The regex pattern.
        '''

        if image.color_profile in ('hdr', 'hdra'):
            label = 'mSSIM (RGB)'
        elif image.color_format == 'rgba':
            label = 'SSIM (LDR-RGBA)'
        else:
            label = 'SSIM (LDR-RGB)'

        return re.compile(rf'\s*{re.escape(label)}:\s*([-+0-9.eE]+)(?:\s|$)')

    def get_total_time_pattern(self) -> re.Pattern:
        return re.compile(r'\s*Total time:\s*([0-9.]*) s')

    def get_coding_time_pattern(self) -> re.Pattern:
        return re.compile(r'\s*Coding time:\s*([0-9.]*) s')

    def get_coding_rate_pattern(self) -> re.Pattern:
        return re.compile(r'\s*Coding rate:\s*([0-9.]*) MT/s')


class Encoder2xRel(Encoder2x):
    '''
    This class wraps a released 2.x series binary.
    '''
    def __init__(self, version: str, variant: str):

        self.version = version

        binary = 'astcenc'
        post = f'-{variant}' if variant != 'universal' else ''
        ext = '.exe' if os.name == 'nt' else ''
        binary_path = Path('./') / 'Binaries' / version / f'{binary}{post}{ext}'

        super().__init__(variant, binary_path)
