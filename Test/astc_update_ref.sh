#!/bin/bash

if [ -z "$1" ]; then
    echo "ERROR: Missing reference version; use 'main' or '5.7'"
    exit 1
fi

echo "Generating new ref-$1 results"

if [ "$1" = "main" ]; then
    binary_path="./bin/"
else
    binary_path="./Binaries/${1}/"
fi

if [ ! -d "$binary_path" ]; then
    echo "ERROR: Binary directory not found: $binary_path" >&2
    exit 1
fi

echo "Using binaries from $binary_path"
echo ""


#python3 ./Test/astc_test_image.py --test-set all --block-size all --test-quality all --repeats 6 --reference rebuild --encoder ref-$1-sve2_256 -j12
#python3 ./Test/astc_test_image.py --test-set all --block-size all --test-quality all --repeats 6 --reference rebuild --encoder ref-$1-sve2_128 -j12
#python3 ./Test/astc_test_image.py --test-set all --block-size all --test-quality all --repeats 6 --reference rebuild --encoder ref-$1-sve_256 -j12
#python3 ./Test/astc_test_image.py --test-set all --block-size all --test-quality all --repeats 6 --reference rebuild --encoder ref-$1-sve_128 -j12
#python3 ./Test/astc_test_image.py --test-set all --block-size all --test-quality all --repeats 6 --reference rebuild --encoder ref-$1-neon -j12


#python3 ./Test/astc_test_image.py --test-set all --block-size all --test-quality all --repeats 8 --reference rebuild --encoder "ref-${1}-avx512" -j12
python3 ./Test/astc_test_image.py --test-set all --block-size all --test-quality all --repeats 8 --reference rebuild --encoder "ref-${1}-avx2" -j12
#python3 ./Test/astc_test_image.py --test-set all --block-size all --test-quality all --repeats 8 --reference rebuild --encoder ref-$1-sse4.1 -j12
#python3 ./Test/astc_test_image.py --test-set all --block-size all --test-quality all --repeats 8 --reference rebuild --encoder ref-$1-sse2 -j12

