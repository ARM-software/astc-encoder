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
The python test runner is designed to run some basic tests against the Python
test code base.
'''

import io
import os
from pathlib import Path
import re
import sys
import tempfile
import unittest

import mypy
import mypy.api
import pycodestyle
import pylint
import pylint.lint
from pylint.reporters.text import TextReporter

import astc_test_image as ati
from testlib.resultset import Record, ResultSet
from testlib.testset import TestSet


class ImageSelectionTests(unittest.TestCase):
    '''
    Check pairwise image coverage against the full test matrix.

    Attributes:
        test_set: The Small image test set.
    '''

    def setUp(self):
        '''
        Load the Small test images without running the compressor.
        '''
        self.test_set = TestSet(
            'Small', Path(__file__).parent / 'Images' / 'Small',
            ['ldr', 'ldrs', 'hdr'], ['l', 'xy', 'rgb', 'rgba'])

    def test_pairwise_coverage(self):
        '''
        Cover all feasible pairs and every image without duplicate runs.
        '''
        qualities = ['fastest', 'fast', 'medium', 'thorough']
        full = ati.select_test_cases(
            self.test_set, ati.TEST_BLOCK_SIZES, qualities)
        selected = ati.select_test_cases(
            self.test_set, ati.TEST_BLOCK_SIZES, qualities, True)

        full_values = [
            ((image.color_profile, image.color_format), block, quality)
            for image, block, quality in full
        ]
        selected_values = [
            ((image.color_profile, image.color_format), block, quality)
            for image, block, quality in selected
        ]

        for first, second in [(0, 1), (0, 2), (1, 2)]:
            self.assertEqual(
                {(case[first], case[second]) for case in full_values},
                {(case[first], case[second]) for case in selected_values})

        self.assertEqual(
            {case[0].file_path for case in full},
            {case[0].file_path for case in selected})
        self.assertEqual(len(selected), len(set(selected)))
        self.assertTrue(set(selected).issubset(full))
        self.assertLess(len(selected), len(full))
        self.assertEqual(
            selected, ati.select_test_cases(
                self.test_set, ati.TEST_BLOCK_SIZES, qualities, True))

    def test_filtered_selection(self):
        '''
        Respect block size and quality filters, including empty selections.
        '''
        for blocks in [[], ['4x4'], ['3x3x3']]:
            with self.subTest(blocks=blocks):
                full = ati.select_test_cases(self.test_set, blocks, ['fast'])
                selected = ati.select_test_cases(
                    self.test_set, blocks, ['fast'], True)
                self.assertEqual(full, selected)

    def test_full_selection(self):
        '''
        Preserve the default execution order and dimensionality constraints.
        '''
        expected = [
            (image, block, 'fast')
            for block in ati.TEST_BLOCK_SIZES
            for image in self.test_set.tests
            if image.is_3d == (block.count('x') == 2)
        ]

        self.assertEqual(
            expected, ati.select_test_cases(
                self.test_set, ati.TEST_BLOCK_SIZES, ['fast']))


class ResultSetUpdateTests(unittest.TestCase):
    '''
    Check quality-only CSV updates preserve existing reference data.
    '''

    def test_unchanged_quality(self):
        '''
        Ignore timing and SSIM changes when PSNR rounds to the same value.
        '''
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'results.csv'
            existing = ResultSet('Small')
            existing.add_record(Record('4x4', 'image.png', 30.1235,
                                       0.9, 1.0, 0.8, 2.0))
            existing.save_to_file(path)
            original = path.read_bytes()
            os.utime(path, ns=(1000000000, 1000000000))

            updated = ResultSet('Small')
            updated.add_record(Record('4x4', 'image.png', 30.123456,
                                      0.99, 5.0, 4.0, 0.4))
            updated.save_quality_changes(path)

            self.assertEqual(path.read_bytes(), original)
            self.assertEqual(path.stat().st_mtime_ns, 1000000000)

    def test_partial_quality_update(self):
        '''
        Replace changed scores, append new cases, and keep untested rows.
        '''
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'results.csv'
            existing = ResultSet('Small')
            existing.add_record(Record('4x4', 'same.png', 30.0,
                                       0.9, 1.0, 0.8, 2.0))
            existing.add_record(Record('4x4', 'better.png', 30.0,
                                       0.9, 1.0, 0.8, 2.0))
            existing.add_record(Record('4x4', 'worse.png', 30.0,
                                       0.9, 1.0, 0.8, 2.0))
            existing.add_record(Record('6x6', 'same.png', 28.0,
                                       0.8, 1.0, 0.8, 2.0))
            existing.save_to_file(path)

            updated = ResultSet('Small')
            updated.add_record(Record('4x4', 'same.png', 30.0,
                                      0.99, 5.0, 4.0, 0.4))
            updated.add_record(Record('4x4', 'better.png', 31.0,
                                      0.95, 2.0, 1.6, 1.0))
            updated.add_record(Record('4x4', 'worse.png', 29.0,
                                      0.85, 2.0, 1.6, 1.0))
            updated.add_record(Record('8x8', 'same.png', 27.0,
                                      0.75, 2.0, 1.6, 1.0))
            updated.save_quality_changes(path)

            expected = ResultSet('Small')
            expected.add_record(existing.records[0])
            expected.add_record(updated.records[1])
            expected.add_record(updated.records[2])
            expected.add_record(existing.records[3])
            expected.add_record(updated.records[3])
            expected_path = Path(directory) / 'expected.csv'
            expected.save_to_file(expected_path)
            self.assertEqual(path.read_bytes(), expected_path.read_bytes())

    def test_new_quality_csv(self):
        '''
        Create missing CSVs and preserve the default full-replacement mode.
        '''
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'nested' / 'results.csv'
            results = ResultSet('Small')
            results.add_record(Record('4x4', 'image.png', 30.0,
                                      0.9, 1.0, 0.8, 2.0))
            results.save_quality_changes(path)
            original = path.read_bytes()
            results.records[0].total_time = 5.0
            results.save_to_file(path)

            self.assertNotEqual(path.read_bytes(), original)
            loaded = ResultSet('Small')
            loaded.load_from_file(path)
            self.assertEqual(len(loaded.records), 1)
            self.assertEqual(loaded.records[0].total_time, 5.0)


class PythonTests(unittest.TestCase):
    '''
    Some basic Python static analysis and style checks.
    '''

    def test_pylint(self):
        '''
        Run pylint over the codebase.
        '''
        # Run Pylint
        stream = io.StringIO()
        reporter = TextReporter(stream)
        pylint.lint.Run(['./Test'], reporter, False)
        output = stream.getvalue()

        # Write the Pylint log
        with open('pylint.log', 'w', encoding='utf-8') as handle:
            handle.write(output)

        # Analyze the results
        pattern = re.compile(r'Your code has been rated at (.*?)/10')
        match = pattern.search(output)
        self.assertIsNotNone(match)
        score = float(match.group(1))

        # This target is currently low but we will increase over time
        self.assertGreaterEqual(score, 9.95, 'Found Pylint score regression')

    def test_pycodestyle(self):
        '''
        Run pycodestyle over the codebase.
        '''
        style = pycodestyle.StyleGuide()

        # Write the Pycodestyle log
        with open('pycodestyle.log', 'w', encoding='utf-8') as handle:
            # Redirect stdout to file so we can capture output
            old_stdout = sys.stdout
            sys.stdout = handle

            result = style.check_files(['./Test'])

            # Restore stdout
            sys.stdout = old_stdout

        self.assertEqual(result.total_errors, 0,
                         'Found pycodestyle warnings or errors.')

    def test_mypy(self):
        '''
        Run mypy over the codebase.
        '''
        # pylint: disable=c-extension-no-member
        result = mypy.api.run(['./Test'])

        # Write the mypy log
        if result[0]:
            with open('mypy.log', 'w', encoding='utf-8') as handle:
                handle.write(result[0])

            # Extract actual result lines from the log
            lines = result[0].splitlines()
            error_lines = []
            for line in lines:
                # Skip lines that are just status
                if line.startswith('Success'):
                    continue
                if line.startswith('Found'):
                    continue

                error_lines.append(line)

            self.assertEqual(len(error_lines), 0,
                             'Found mypy warnings or errors.')

        self.assertFalse(bool(result[1]),
                         'mypy failed to run.')


def main():
    '''
    The main function.

    Return:
        int: The process return code.
    '''
    results = unittest.main(exit=False)
    return 0 if results.result.wasSuccessful() else 1


if __name__ == '__main__':
    sys.exit(main())
