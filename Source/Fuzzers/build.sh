#!/bin/bash -eu

#  SPDX-License-Identifier: Apache-2.0
#  ----------------------------------------------------------------------------
#  Copyright 2020 Arm Limited
#  Copyright 2020 Google Inc.
#
#  Licensed under the Apache License, Version 2.0 (the "License"); you may not
#  use this file except in compliance with the License. You may obtain a copy
#  of the License at:
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
#  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
#  License for the specific language governing permissions and limitations
#  under the License.
#  ----------------------------------------------------------------------------

# Build the core project
make -j$(nproc) BUILD=debug VEC=sse2

# Package up a library for fuzzers to link against
ar -qc libastcenc.a  *.o

# Build project local fuzzers
for fuzzer in $SRC/astc-encoder/Source/Fuzzers/fuzz_*.cpp; do
  $CXX $CXXFLAGS \
      -DASTCENC_SSE=0 \
      -DASTCENC_AVX=0 \
      -DASTCENC_POPCNT=0 \
      -DASTCENC_VECALIGN=16 \
      -DASTCENC_ISA_INVARIANCE=0 \
      -I. -std=c++14 $fuzzer $LIB_FUZZING_ENGINE $SRC/astc-encoder/Source/libastcenc.a \
      -o $OUT/$(basename -s .cpp $fuzzer)
done
