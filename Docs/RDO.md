# Rate-distortion optimization

Rate-distortion optimization (RDO) trades image quality for smaller files after
applying additional lossless compression to compressed ASTC files. It is used
to reduce the size of application distribution bundles, providing finer-grained
control over storage size than simply reducing ASTC block size.

> [!NOTE]
> The RDO implementation in `astcenc` currently only supports LDR textures.

## How RDO works

ASTC compression is a fixed-rate compression that encodes every block into 16
bytes of data. For a given image size and block footprint, the ASTC payload
data has a fixed size. RDO is a post-process that changes the compressed ASTC
bitstream to make byte patterns more similar to nearby bytes in the stream.
This makes the data more compressible when applying subsequent lossless
compression, such as ZIP or LZ4.

The RDO transform is lossy, sacrificing image quality to improve data
compressibility. The amount of allowed quality loss is controlled by a lambda
setting, `-rdo-lambda` for the command line tool, where a higher lambda values
allow more image distortion in exchange for a more compressible pattern. A
lambda value of zero disables the RDO transform. We recommend using a value
of 20 as a starting point, which give a reasonable tradeoff between quality
and storage. Normal maps are more sensitive to introduced error, so a lower
lambda value might be required.

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
additional savings to run-time memory bandwidth, in addition to the packging
size reduction.

We recommend using RDO to fine-tune packaging size when dropping to the next
lowest ASTC bitrate is not an option. This could be because the next bitrate
option loses too much detail, or suffers from more visible block artifacts.

## Compressing with RDO

Use the `-rdo-lambda` command line option to enable RDO and set its strength.
Setting 0 disables RDO, higher values increase the amount of quality loss
permitted.

```shell
# Reference encode without RDO
build-rdo/Source/astcenc-avx2 -cl input.png reference.astc 5x5 -thorough \
    -rdo-lambda 0

# Same footprint and search preset, with a light RDO pass enabled
build-rdo/Source/astcenc-avx2 -cl input.png rdo.astc 5x5 -thorough \
    -rdo-lambda 20
```

The RDO lambda value is a heuristic to control tradeoffs, not a specific target
bitrate or PSNR threshold. Lower positive values favor fidelity; higher values
place more weight on storage savings. You much measure quality and check that
it is suitable for your content.

- - -

_Copyright © 2026, Arm Limited and contributors._
