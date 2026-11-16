# Rate-distortion optimization

Rate-distortion optimization (RDO) trades image quality for smaller files after
applying additional lossless compression to compressed ASTC files. It is used
to reduce the size of application distribution bundles, providing finer-grained
control over storage size than simply reducing ASTC block size.

> [!NOTE]
> The RDO implementation in `astcenc` currently supports 2D LDR textures, using
> either the linear or sRGB color profile.

## How RDO works

ASTC compression is a fixed-rate compression that encodes every block into 16
bytes of data. For a given image size and block footprint, the ASTC payload
data has a fixed size. RDO is a post-process that changes the compressed ASTC
bitstream to make byte patterns more similar to nearby bytes in the stream.
This makes the data more compressible when applying subsequent lossless
compression, such as ZIP or LZ4.

The RDO transform is lossy, sacrificing image quality to improve data
compressibility. The amount of allowed quality loss is controlled by a lambda
setting, where a higher lambda values allow more image distortion in exchange
for a more compressible pattern. Lambda values between 0.0 (no RDO) and 0.3
(moderate RDO) give a reasoanble tradeoff between quality and storage. Lambda
values above 0.3 quickly loose image quality for relatively little saving in
packaged file size, and it is usally a better option to choose a larger ASTC
block footprint.

The ASTC blocks after RDO has been applied are still standard ASTC data blocks.
RDO does not require a special GPU decoder, and it does not change the number of
blocks or the size of each block. To realize the storage saving, additional
lossless compression must be applied in your asset packaging pipeline. For
Android, this is compression is commonly applied when generating an Application
APK (application bundle) or AAB (asset bundle).

### Using RDO versus changing the ASTC block size

When using ASTC, you have two controls that can reduce packaged file size. You
can choose to use a larger ASTC block size, which simply reduces the amount of
data that needs to be packaged, or you can choose to use RDO to make data more
compressible.

Where possible we recommend using a larger ASTC block size. This gives you
reliable savings to run-time memory footprint and memory bandwidth, in
addition to the packging size reduction.

We recommend using RDO to fine-tune packaging size when dropping to the next
lowest ASTC bitrate is not an option. This could be because the next bitrate
option loses too much detail, or suffers from more visible block artifacts.

## Building with RDO

RDO support in the codec is disabled by default. Enable it using the
`ASTCENC_RDO` CMake option. For example, from the repository root on an
x86 machine supporting AVX2:

```shell
cmake -S . -B build-rdo -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER=clang++ -DASTCENC_ISA_AVX2=ON -DASTCENC_RDO=ON

cmake --build build-rdo --config Release --parallel
```

> [!NOTE]
> See [Building ASTC Encoder](Building.md) for compiler requirements and other
> build configurations.

On Linux, optionally add `-DASTCENC_RDO_SIZE=ON` when configuring to report the
serialized file size and its DEFLATE-compressed size. This requires the zlib
development package to be installed in your system:

```shell
sudo apt-get install zlib1g-dev
```

Reporting the saving using DEFLATE does not actually compress the output file;
your packaging pipeline must still perform lossless compression.

## Compressing with RDO

Use the `-rdo-lambda` command line option to enable RDO and set its strength.
Setting 0.0 disables RDO, higher values increase the amount of quality loss
permitted.

```shell
# Reference encode without RDO
build-rdo/Source/astcenc-avx2 -cl input.png reference.astc 5x5 -thorough \
    -rdo-lambda 0

# Same footprint and search preset, with RDO enabled
build-rdo/Source/astcenc-avx2 -cl input.png rdo.astc 5x5 -thorough \
    -rdo-lambda 0.2
```

The RDO lambda value is a heuristic to control tradeoffs, not a specific target
bitrate or PSNR threshold. Lower positive values favor fidelity; higher values
place more weight on storage savings. A useful initial experiment is to compare
0.0, 0.1, 0.2, and 0.5 at the same block footprint and quality preset. Measure
quality and final packaged size for your own content.

## Morton order file layout

Ordinary `.astc` files, and other formats such as KTX, store ASTC blocks in
scanline order. Nearby rows in the image, that are likely to have more similar
colors, can therefore be far apart in the byte stream.

The `.astcm` file format is a variant of the `.astc` format that stores the
astc blocks in the image as a sequence of 2D "tiles" of ASTC blocks arranged
in Morton order. This groups spatial neighbors together, which improves the
efficiency of the RDO data conditioning and the effectiveness of the
packaging-level compression.

Select Morton order by writing to an `.astcm` file:

```shell
build-rdo/Source/astcenc-avx2 -cl input.png rdo.astcm 5x5 -thorough \
    -rdo-lambda 0.2
```

The `.astcm` file layout stores 8x8 tiles of ASTC blocks, with tiles being
stored sequentially in scanline order. Blocks within each tile are stored in
Morton order. Partial edge tiles omit blocks that fall outside the image.

> [!NOTE]
> See [File formats](FileFormat.md#the-astcm-format) for the header and exact
> block traversal.

### Data loading

The `.astcm` data ordering is a storage layout that improves packing
efficiency, but it is not the layout that the graphcis APs expect.
Before you can upload the data to OpenGL ES or Vulkan you must convert
the data back into scanline block order:

### Morton order in the library API

When passing `ASTCENC_FLG_USE_MORTON` to the library API, it will issue tasks
to the compressor in blocks of one or more 8x8 Morton-order "tiles". This
allows the RDO transform to exploit 2D similarity in the image.

However, it is important to note that the compressed ASTC data returned by the
codec API will be returned in block-by-block scanline order. Tool file
serialization must also apply the 8x8 Morton transform to realize the packaging
benefits. If you are not using Morton-order storage files, it is recommended
that you do not pass `ASTCENC_FLG_USE_MORTON` to the library API,

## RDO Entropy Reduction Transform

The ASTC Entroy Reduction Transform (AERT) used in our RDO implementation is
based on the public-domain Entropy Reduction Transform (ERT) implementation
developed by Rich Geldrich for [BC7Enc RDO][1]. The underlying algorithm is
the same, but it has been heavily optimized for ASTC to improve RDO
processing performance.

- - -

_Copyright © 2026, Arm Limited and contributors._

[1]: https://github.com/richgel999/bc7enc_rdo/blob/master/ert.h
