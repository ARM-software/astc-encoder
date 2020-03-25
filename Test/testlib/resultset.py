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
A ResultSet stores a set of results about the performance of a TestSet. Each
set keeps result Records for each image and block size tested, that store the
PSNR and coding time.

ResultSets are often backed by a CSV file on disk, and a ResultSet can be
compared against a set of reference results created by an earlier test run.
"""


import csv
import enum


@enum.unique
class Result(enum.IntEnum):
    """
    An enumeration of test result status values.

    Attributes:
        NOTRUN: The test has not been run.
        PASS: The test passed.
        WARN: The test image quality was below the pass threshold but above
            the fail threshold.
        FAIL: The test image quality was below the fail threshold.
    """
    NOTRUN = 0
    PASS = 1
    WARN = 2
    FAIL = 3


class ResultSummary():
    """
    An result summary data container, storing number of results of each type.

    Attributes:
        notruns: The number of tests that did not run.
        passes: The number of tests that passed.
        warnings: The number of tests that produced a warning.
        fails: The number of tests that failed.
    """

    def __init__(self):
        """
        Create a new result summary.
        """
        self.notruns = 0
        self.passes = 0
        self.warnings = 0
        self.fails = 0

    def add_record(self, record):
        """
        Add a record to this summary.

        Args:
            record: The Record to add.
        """
        if record.status == Result.PASS:
            self.passes += 1
        elif record.status == Result.WARN:
            self.warnings += 1
        elif record.status == Result.FAIL:
            self.fails += 1
        else:
            self.notruns += 1

    def get_worst_result(self):
        """
        Get the worst result in this set.

        Returns:
            The worst Result value.
        """
        if self.fails:
            return Result.FAIL

        if self.warnings:
            return Result.WARN

        if self.passes:
            return Result.PASS

        return Result.NOTRUN

    def __str__(self):
        overall = self.get_worst_result().name
        dat = (overall, self.passes, self.warnings, self.fails)
        return "Set Status: %s (Pass: %u | Warn: %u | Fail: %u)" % dat


class Record():
    """
    A single result record, sotring results for a singel image and block size.

    Attributes:
        blkSz: The block size.
        name: The test image name.
        psnr: The image quality (PSNR dB)
        tTime: The total compression time.
        cTime: The coding compression time.
        status: The test Result.
    """

    def __init__(self, blkSz, name, psnr, tTime, cTime):
        """
        Create a result record, initially in the NOTRUN status.

        Args:
            blkSz: The block size.
            name: The test image name.
            psnr: The image quality (PSNR dB)
            tTime: The total compression time.
            cTime: The coding compression time.
        """
        self.blkSz = blkSz
        self.name = name
        self.psnr = psnr
        self.tTime = tTime
        self.cTime = cTime
        self.status = Result.NOTRUN

    def set_status(self, result):
        """
        Set the result status.

        Args:
            result: The test Result.
        """
        self.status = result

    def __str__(self):
        return "'%s' / '%s'" % (self.blkSz, self.name)


class ResultSet():
    """
    A set of results for a TestSet, across one or more block sizes.

    Attributes:
        testSet: The test set these results are linked to.
        records: The list of test results.
    """

    def __init__(self, testSet):
        """
        Create a new empty ResultSet.

        Args:
            testSet: The test set these results are linked to.
        """
        self.testSet = testSet
        self.records = []

    def add_record(self, record):
        """
        Add a new test record to this result set.

        Args:
            record: The test record to add.
        """
        self.records.append(record)

    def get_record(self, testSet, blkSz, name):
        """
        Get a record matching the arguments.

        Args:
            testSet: The test set to get results from.
            blkSz: The block size.
            name: The test name.

        Returns:
            The result Record, if present.

        Raises:
            KeyError: No match could be found.
        """
        if testSet != self.testSet:
            raise KeyError()

        for record in self.records:
            if record.blkSz == blkSz and record.name == name:
                return record

        raise KeyError()

    def get_matching_record(self, other):
        """
        Get a record matching the config of another record.

        Args:
            other: The pattern result Record.

        Returns:
            The result Record, if present.

        Raises:
            KeyError: No match could be found.
        """
        for record in self.records:
            if record.blkSz == other.blkSz and record.name == other.name:
                return record

        raise KeyError()

    def get_results_summary(self):
        """
        Get a results summary of all the records in this result set.

        Returns:
            The result summary.
        """
        summary = ResultSummary()
        for record in self.records:
            summary.add_record(record)
        return summary

    def save_to_file(self, filePath):
        """
        Save this result set to a CSV file.

        Args:
            filePath: The output file path.
        """
        with open(filePath, "w") as csvfile:
            writer = csv.writer(csvfile)
            self._save_header(writer)
            for record in self.records:
                self._save_record(writer, record)

    @staticmethod
    def _save_header(writer):
        """
        Write the header to the CSV file.

        Args:
            writer: The CSV writer.
        """
        row = ["Image Set", "Block Size", "Name",
               "PSNR", "Total Time", "Coding Time"]
        writer.writerow(row)

    def _save_record(self, writer, record):
        """
        Write a record to the CSV file.

        Args:
            writer: The CSV writer.
            record: The record to write.
        """
        row = [self.testSet,
               record.blkSz,
               record.name,
               "%0.5f" % record.psnr,
               "%0.3f" % record.tTime,
               "%0.3f" % record.cTime]
        writer.writerow(row)

    def load_from_file(self, filePath):
        """
        Load a reference result set from a CSV file on disk.

        Args:
            filePath: The input file path.
        """
        with open(filePath, "r") as csvfile:
            reader = csv.reader(csvfile)
            # Skip the header
            next(reader)
            for row in reader:
                assert row[0] == self.testSet
                record = Record(row[1], row[2],
                                float(row[3]), float(row[4]), float(row[5]))
                self.add_record(record)
