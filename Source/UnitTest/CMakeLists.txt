#  SPDX-License-Identifier: Apache-2.0
#  ----------------------------------------------------------------------------
#  Copyright 2020 Arm Limited
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

if (${UNIVERSAL_BUILD})
    include(cmake_core.cmake)
else()
    set(artefacts none native neon avx2 sse4.1 sse2)
    set(configs ${ISA_NONE} ${ISA_NATIVE} ${ISA_NEON} ${ISA_AVX2} ${ISA_SSE41} ${ISA_SSE2})
    list(LENGTH artefacts artefacts_len)
    math(EXPR artefacts_len "${artefacts_len} - 1")

    foreach(index RANGE ${artefacts_len})
        list(GET artefacts ${index} artefact)
        list(GET configs ${index} config)
        if(${config})
            set(ISA_SIMD ${artefact})
            include(cmake_core.cmake)
        endif()
    endforeach()
endif()