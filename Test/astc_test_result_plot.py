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
The ``astc_test_result_plot.py`` script consolidates all current sets of
reference results into a single graphical plot.
"""

import re
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

import testlib.resultset as trs
from collections import defaultdict as ddict


def find_reference_results():
    """
    Scrape the Test/Images directory for result CSV files and return an
    mapping of the result sets.

    Returns:
        Returns a three deep tree of dictionaries, with the final dict
        pointing at a `ResultSet` object. The hierarchy is:

            imageSet => quality => encoder => result
    """
    scriptDir = os.path.dirname(__file__)
    imageDir = os.path.join(scriptDir, "Images")

    # Pattern for extracting useful data from the CSV file name
    filePat = re.compile(r"astc_reference-(.+)_(.+)_results\.csv")

    # Build a three level dictionary we can write into
    results = ddict(lambda: ddict(lambda: ddict()))

    # Final all CSVs, load them and store them in the dict tree
    for root, dirs, files in os.walk(imageDir):
        for name in files:
            match = filePat.match(name)
            if match:
                fullPath = os.path.join(root, name)

                encoder = match.group(1)
                quality = match.group(2)
                imageSet = os.path.basename(root)

                testRef = trs.ResultSet(imageSet)
                testRef.load_from_file(fullPath)

                results[imageSet][quality]["ref-%s" % encoder] = testRef

    return results


def get_series(results, tgtEncoder, tgtQuality, resFilter=lambda x: True):
    psnrData = []
    mtsData = []
    marker = []
    records = []

    for imageSet, iResults in results.items():

        for quality, qResults in iResults.items():
            if quality != tgtQuality:
                continue

            for encoder, eResults in qResults.items():
                if encoder != tgtEncoder:
                    continue

                for record in eResults.records:
                    if resFilter(record):
                        records.append(record)
                        psnrData.append(record.psnr)
                        mtsData.append(record.cRate)

                        if "ldr-xy" in record.name:
                            marker.append('$N$')
                        elif "ldr-l" in record.name:
                            marker.append('$G$')
                        elif "ldr" in record.name:
                            marker.append('$L$')
                        elif "hdr" in record.name:
                            marker.append('$H$')
                        else:
                            marker.append('$?$')


    return mtsData, psnrData, marker, records


def get_series_rel(results, refEncoder, refQuality, tgtEncoder, tgtQuality, resFilter=lambda x: True):

    mts1, psnr1, marker1, rec1 = get_series(results, tgtEncoder, tgtQuality, resFilter)

    if refEncoder is None:
        refEncoder = tgtEncoder

    if refQuality is None:
        refQuality = tgtQuality

    mts2, psnr2, marker2, rec2 = get_series(results, refEncoder, refQuality, resFilter)

    mtsm  = [x/mts2[i] for i, x in enumerate(mts1)]
    psnrm = [x - psnr2[i] for i, x in enumerate(psnr1)]

    return mtsm, psnrm, marker1, rec1


def get_human_eq_name(encoder, quality):
    parts = encoder.split("-")
    if len(parts) == 2:
        return "astcenc %s -%s" % (parts[1], quality)
    else:
        return "astcenc-%s %s -%s" % (parts[2], parts[1], quality)


def get_human_e_name(encoder):
    parts = encoder.split("-")
    if len(parts) == 2:
        return "astcenc %s" % parts[1]
    else:
        return "astcenc-%s %s" % (parts[2], parts[1])


def get_human_q_name(quality):
    return "-%s" % quality


def plot(results, chartRows, chartCols, blockSizes,
         relative, pivotEncoder, pivotQuality, fileName, limits):

    fig, axs = plt.subplots(nrows=len(chartRows), ncols=len(chartCols),
                            sharex=True, sharey=True, figsize=(15, 8.43))

    for a in fig.axes:
        a.tick_params(
            axis="x", which="both",
            bottom=True, top=False, labelbottom=True)

        a.tick_params(
            axis="y", which="both",
            left=True, right=False, labelleft=True)

    for i, row in enumerate(chartRows):
        for j, col in enumerate(chartCols):
            if row == "fastest" and "ref-2.1" not in col:
                if len(chartCols) == 1:
                    fig.delaxes(axs[i])
                else:
                    fig.delaxes(axs[i][j])
                continue

            if len(chartRows) == 1 and len(chartCols) == 1:
                ax = axs
            elif len(chartCols) == 1:
                ax = axs[i]
            else:
                ax = axs[i, j]

            title = get_human_eq_name(col, row)

            if not relative:
                ax.set_title(title, y=0.97, backgroundcolor="white")
                ax.set_xlabel('Coding performance (MTex/s)')
                ax.set_ylabel('PSNR (dB)')
            else:
                if pivotEncoder and pivotQuality:
                    tag = get_human_eq_name(pivotEncoder, pivotQuality)
                elif pivotEncoder:
                    tag = get_human_e_name(pivotEncoder)
                else:
                    assert(pivotQuality)
                    tag = get_human_q_name(pivotQuality)

                ax.set_title("%s vs. %s" % (title, tag), y=0.97, backgroundcolor="white")
                ax.set_xlabel('Performance scaling')
                ax.set_ylabel('PSNR delta (dB)')

            for k, series in enumerate(blockSizes):
                fn = lambda x: x.blkSz == series

                if not relative:
                    x, y, m, r = get_series(results, col, row, fn)
                else:
                    x, y, m, r = get_series_rel(results, pivotEncoder, pivotQuality,
                                                col, row, fn)

                color = None
                label = "%s blocks" % series
                for xp, yp, mp in zip(x, y, m):
                    ax.scatter([xp],[yp], s=16, marker=mp,
                               color="C%u" % k, label=label)
                    label = None

            if i == 0 and j == 0:
                ax.legend(loc="lower right")

    for i, row in enumerate(chartRows):
        for j, col in enumerate(chartCols):

            if len(chartRows) == 1 and len(chartCols) == 1:
                ax = axs
            elif len(chartCols) == 1:
                ax = axs[i]
            else:
                ax = axs[i, j]

            ax.grid(ls=':')

            if not relative:
                ax.set_xlim(left=0, right=limits[0])
            else:
                ax.set_xlim(left=1, right=limits[0])

    fig.tight_layout()
    fig.savefig(fileName)


def main():
    """
    The main function.

    Returns:
        int: The process return code.
    """

    charts = [
        [
            # Plot headline 1.7 to 2.1 scores
            ["medium"],
            ["ref-2.1-avx2"],
            ["4x4", "6x6", "8x8"],
            True,
            "ref-1.7",
            None,
            "relative-1.7-to-2.1.png",
            (8, None)
        ], [
            # --------------------------------------------------------
            # Plot all absolute scores
            ["thorough", "medium", "fast", "fastest"],
            ["ref-2.0-avx2", "ref-2.1-avx2"],
            ["4x4", "5x5", "6x6", "8x8"],
            False,
            None,
            None,
            "absolute-all.png",
            (30, None)
        ], [
            # Plot 1.7 to 2.0 absolute scores
            ["thorough", "medium", "fast", "fastest"],
            ["ref-1.7", "ref-2.0-avx2"],
            ["4x4", "5x5", "6x6", "8x8"],
            False,
            None,
            None,
            "absolute-1.7-to-2.0.png",
            (30, None)
        ], [
            # Plot 2.0 to 2.1 absolute scores
            ["thorough", "medium", "fast", "fastest"],
            ["ref-2.0-avx2", "ref-2.1-avx2"],
            ["4x4", "5x5", "6x6", "8x8"],
            False,
            None,
            None,
            "absolute-2.0-to-2.1.png",
            (30, None)
        ], [
            # --------------------------------------------------------
            # Plot all relative scores
            ["thorough", "medium", "fast"],
            ["ref-2.0-avx2", "ref-2.1-avx2"],
            ["4x4", "5x5", "6x6", "8x8"],
            True,
            "ref-1.7",
            None,
            "relative-all.png",
            (None, None)
        ], [
            # Plot 1.7 to 2.0 relative scores
            ["thorough", "medium", "fast"],
            ["ref-2.0-avx2"],
            ["4x4", "5x5", "6x6", "8x8"],
            True,
            "ref-1.7",
            None,
            "relative-1.7-to-2.0.png",
            (None, None)
        ], [
            # Plot 2.0 to 2.1 relative scores
            ["thorough", "medium", "fast"],
            ["ref-2.1-avx2"],
            ["4x4", "5x5", "6x6", "8x8"],
            True,
            "ref-2.0-avx2",
            None,
            "relative-2.0-to-2.1.png",
            (None, None)
        ], [
            # Plot relative scores of ISAs of latest
            ["thorough", "medium", "fast", "fastest"],
            ["ref-2.1-sse4.2", "ref-2.1-avx2"],
            ["4x4", "5x5", "6x6", "8x8"],
            True,
            "ref-2.1-sse2",
            None,
            "relative-2.1-isas.png",
            (None, None)
        ], [
            # Plot relative scores of qualities of latest
            ["medium", "fast", "fastest"],
            ["ref-2.1-avx2"],
            ["4x4", "5x5", "6x6", "8x8"],
            True,
            None,
            "thorough",
            "relative-2.1-qualities.png",
            (None, None)
        ]
    ]

    results = find_reference_results()

    for chart in charts:
        plot(results, *chart)

    return 0


if __name__ == "__main__":
    sys.exit(main())
