#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# -----------------------------------------------------------------------------
# Copyright 2021 Arm Limited
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
The ``astc_trace_analysis`` utility provides a tool to analyze trace files.

WARNING: Trace files are an engineering tool, and not part of the standard
product, so traces and their associated tools are volatile and may change
significantly without notice.
"""

import argparse
import json
import sys

QUANT_TABLE = {
	 0:   2,
	 1:   3,
	 2:   4,
	 3:   5,
	 4:   6,
	 5:   8,
	 6:  10,
	 7:  12,
	 8:  16,
	 9:  20,
	10:  24,
	11:  32
}


def rev_enumerate(seq):
    return zip(reversed(range(len(seq))), reversed(seq))


class Trace:

    def __init__(self, block_x, block_y, block_z):
        self.block_x = block_x
        self.block_y = block_y
        self.block_z = block_z
        self.blocks = []

    def add_block(self, block):
        self.blocks.append(block)

    def __getitem__(self, i):
       return self.blocks[i]

    def __delitem__(self, i):
       del self.blocks[i]

    def __len__(self):
       return len(self.blocks)

class Block:

    def __init__(self, pos_x, pos_y, pos_z):
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.pos_z = pos_z

        self.passes = []
        self.qualityHit = None

    def add_pass(self, pas):
        self.passes.append(pas)

    def __getitem__(self, i):
       return self.passes[i]

    def __delitem__(self, i):
       del self.passes[i]

    def __len__(self):
       return len(self.passes)

class Pass:

    def __init__(self, partitions, planes, target_hit, mode=None, channel=None):
        self.partitions = partitions
        self.planes = planes
        self.plane2_channel = channel
        self.target_hit = target_hit
        self.search_mode = mode
        self.candidates = []

    def add_candidate(self, candidate):
        self.candidates.append(candidate)

    def __getitem__(self, i):
       return self.candidates[i]

    def __delitem__(self, i):
       del self.candidates[i]

    def __len__(self):
       return len(self.candidates)


class Candidate:

    def __init__(self, weight_x, weight_y, weight_z, weight_quant):
        self.weight_x = weight_x
        self.weight_y = weight_y
        self.weight_z = weight_z
        self.weight_quant = weight_quant
        self.refinement_errors = []

    def add_refinement(self, errorval):
        self.refinement_errors.append(errorval)


def get_attrib(data, name, multiple=False, hard_fail=True):
    results = []
    for attrib in data:
        if len(attrib) == 2 and attrib[0] == name:
            results.append(attrib[1])

    if not results:
        if hard_fail:
            print(json.dumps(data, indent=2))
            assert False, "Attribute %s not found" % name
        return None

    if not multiple:
        if len(results) > 1:
            print(json.dumps(data, indent=2))
            assert False, "Attribute %s found %u times" % (name, len(results))
        return results[0]

    return results


def get_node(data, name, multiple=False, hard_fail=True):
    results = []
    for attrib in data:
        if len(attrib) == 3 and attrib[0] == "node" and attrib[1] == name:
            results.append(attrib[2])

    if not results:
        if hard_fail:
            print(json.dumps(data, indent=2))
            assert False, "Node %s not found" % name
        return None

    if not multiple:
        if len(results) > 1:
            print(json.dumps(data, indent=2))
            assert False, "Node %s found %u times" % (name, len(results))
        return results[0]

    return results


def find_best_pass_and_candidate(block):
    explicit_pass = None

    best_error = 1e30
    best_pass = None
    best_candidate = None

    for pas in block:
        # Special case for constant color blocks - no trial candidates
        if pas.target_hit and pas.partitions == 0:
            return (pas, None)

        for candidate in pas:
            errorval = candidate.refinement_errors[-1]
            if errorval <= best_error:
                best_error = errorval
                best_pass = pas
                best_candidate = candidate

    # Every other return type must have both best pass and best candidate
    assert (best_pass and best_candidate)
    return (best_pass, best_candidate)


def generate_database(data):
    # Skip header
    assert(data[0] == "node")
    assert(data[1] == "root")
    data = data[2]

    bx = get_attrib(data, "block_x")
    by = get_attrib(data, "block_y")
    bz = get_attrib(data, "block_z")
    dbStruct = Trace(bx, by, bz)

    for block in get_node(data, "block", True):
        bx = get_attrib(block, "x_pos")
        by = get_attrib(block, "y_pos")
        bz = get_attrib(block, "z_pos")
        blockStruct = Block(bx, by, bz)
        dbStruct.add_block(blockStruct)

        for pas in get_node(block, "pass", True):
            # Don't copy across passes we skipped due to heuristics
            skipped = get_attrib(pas, "skip", False, False)
            if skipped:
                continue

            prts = get_attrib(pas, "partition_count")
            plns = get_attrib(pas, "plane_count")
            chan = get_attrib(pas, "plane_channel", False, plns > 2)
            mode = get_attrib(pas, "search_mode", False, False)
            ehit = get_attrib(pas, "exit", False, False) == "quality hit"

            passStruct = Pass(prts, plns, ehit, mode, chan)
            blockStruct.add_pass(passStruct)

            # Constant color blocks don't have any candidates
            if prts == 0:
                continue

            for candidate in get_node(pas, "candidate", True):
                # Don't copy across candidates we couldn't encode
                failed = get_attrib(candidate, "failed", False, False)
                if failed:
                    continue

                wx = get_attrib(candidate, "weight_x")
                wy = get_attrib(candidate, "weight_y")
                wz = get_attrib(candidate, "weight_z")
                wq = QUANT_TABLE[get_attrib(candidate, "weight_quant")]
                epre = get_attrib(candidate, "error_prerealign", True)
                epst = get_attrib(candidate, "error_postrealign")

                candStruct = Candidate(wx, wy, wz, wq)
                passStruct.add_candidate(candStruct)
                for value in epre:
                    candStruct.add_refinement(value)
                candStruct.add_refinement(epst)

    return dbStruct


def filter_database(data):

    for block in data:
        best_pass, best_candidate = find_best_pass_and_candidate(block)

        for i, pas in rev_enumerate(block):
            if pas != best_pass:
                del block[i]
                continue

            if best_candidate is None:
                continue

            for j, candidate in rev_enumerate(pas):
                if candidate != best_candidate:
                    del pas[j]


def parse_command_line():
    """
    Parse the command line.

    Returns:
        Namespace: The parsed command line container.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("trace", type=argparse.FileType("r"),
                        help="The trace file to analyze")

    return parser.parse_args()


def main():
    """
    The main function.

    Returns:
        int: The process return code.
    """
    args = parse_command_line()

    data = json.load(args.trace)
    db = generate_database(data)
    filter_database(db)



    #count_blocks(data)

    return 0


if __name__ == "__main__":
    sys.exit(main())
