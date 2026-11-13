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
Confidence intervals and paired tests for codec timing measurements.
'''

import math
import statistics

from scipy import stats


def ci_half_width(times: list[float], confidence: float) -> float:
    '''
    Compute the two-sided Student's t confidence interval half-width.

    Args:
        times: At least two measured coding times.
        confidence: The confidence level as a fraction in (0, 1).

    Return:
        The confidence interval half-width in seconds.
    '''
    critical = stats.t.ppf((1.0 + confidence) * 0.5, len(times) - 1)

    return float(critical * stats.sem(times))


def relative_ci_half_width(times: list[float], confidence: float) -> float:
    '''
    Compute the timing interval half-width relative to the sample mean.

    Args:
        times: At least two positive measured coding times.
        confidence: The confidence level as a fraction in (0, 1).

    Return:
        The interval half-width as a fraction of the mean coding time.
    '''

    return ci_half_width(times, confidence) / statistics.mean(times)


def is_timing_stable(times: list[float], confidence: float = 0.95, precision: float = 0.01) -> bool:
    '''
    Check whether a build's coding times have a sufficiently precise mean.

    Args:
        times: Coding time samples in seconds from the current build.
        confidence: The confidence level as a fraction between zero and one.
        precision: The allowed interval half-width as a fraction of the mean.

    Return:
        True if the two-sided Student's t confidence interval has a half-width within the requested
        precision of the positive sample mean.
    '''
    count = len(times)

    if count < 2:
        return False

    mean = sum(times) / count

    if mean <= 0.0:
        return False

    half_width = ci_half_width(times, confidence)

    return bool(half_width <= precision * mean)


def paired_t_test(sample_a: list[float], sample_b: list[float]) -> tuple[float, float]:
    '''
    Run a paired, two-sided Student's t-test, including constant differences.

    Args:
        sample_a: The timing samples for binary A.
        sample_b: The timing samples for binary B.

    Return:
        The t-statistic and two-sided p-value. Identical samples give p = 1.

    Raises:
        ValueError: The sample sizes differ or contain fewer than two values.
    '''

    if len(sample_a) != len(sample_b) or len(sample_a) < 2:
        raise ValueError('Paired t-test requires equal sample set sizes')

    diffs = [a - b for a, b in zip(sample_a, sample_b)]

    if all(diff == diffs[0] for diff in diffs):
        if diffs[0] == 0:
            return 0.0, 1.0

        return math.copysign(math.inf, diffs[0]), 0.0

    result = stats.ttest_rel(sample_a, sample_b)

    return float(result.statistic), float(result.pvalue)
