#!/bin/bash

./astcenc/astcenc-avx2-good -tl in.png out-good.png 6x6 -medium

./astcenc/astcenc-avx2 -tl in.png out-bad.png 6x6 -medium

./astcenc/astc_test_autoextract 6x6 in.png out-good.png out-bad.png in-min.png
