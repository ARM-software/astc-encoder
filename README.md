# About

This is the official repository for the Arm® Adaptive Scalable Texture
Compression (ASTC) Encoder, `astcenc`, a command-line tool for compressing
and decompressing images using the ASTC texture compression standard.

## The ASTC format

The ASTC compressed data format, developed by ARM® and AMD, has been adopted as
an official extension to the Open GL®, OpenGL ES, and Vulkan graphics APIs. It
provides a major step forward both in terms of image quality at a given
bitrate, and in terms of the format and bitrate flexibility available to
content creators. This allows more assets to use compression, often at a
reduced bitrate compared to legacy formats, reducing memory bandwidth and
energy consumption.

The ASTC format specification is available here:

* [ASTC Specification v1.0 (PDF)](./Documentation/ASTC%20Specification%20v1.0.pdf)

# Encoder feature support

The encoder supports compression of PNG, TGA and KTX input images into ASTC
format output images. The encoder supports decompression of ASTC input images
into TGA or KTX format output images.

The encoder allows control over the compression time/quality tradeoff with
`exhaustive`, `thorough`, `medium`, `fast`, and `very fast` encoding speeds.

The encoder allows compression time and quality analysis by reporting the
compression time, and the Peak Signal-to-Noise Ratio (PSNR) between the input
image and the compressed output.

## ASTC format support

The ASTC specification allows three profiles of implementation:

* 2D Low Dynamic Range (LDR)
* 2D LDR and High Dynamic Range (HDR)
* 2D and 3D, LDR and HDR

The `astcenc` compressor supports generation of images for all three profiles.
In addition it also supports all of the ASTC block sizes and compression
modes, allowing content creators access the full spectrum of quality-to-bitrate
trade options ranging from 0.89 bits/pixel up to 8bits/pixel.

# License

By downloading any component from this repository you acknowledge that you
accept the End User Licence Agreement for the ASTC Encoder. See the
[license.txt](license.txt) file for details.

# Prebuilt binaries

Prebuilt release build binaries for Windows (x86 and x64), Linux (x86 and x64),
and macOS (x64) are available here:

* [Binary directory](/Binary/).

# Building from source

Builds for Linux and macOS use GCC and Make, and are tested with GCC 4.6 and
GNU Make 3.82.

```
cd Source
make
```

Builds for Windows platforms use Visual Studio 2017, using the solution file
located in the `Source/win32-2017/astcenc/` directory.

# Getting started

Open a terminal, change to the appropriate directory for your system, and run
the astcenc encoder program, like this on Linux or Mac OS:

    ./astcenc

Or like this on Windows:

    astcenc

Invoking the tool with no arguments gives an extensive help message, including
usage instructions, and details of all the available options.

## Compressing an image

Compress an image using the `-c` option:

    astcenc -c example.png example.astc 6x6 -medium

This compresses `example.png` using the 6x6 block footprint (3.55 bits/pixel)
and a `medium` compression speed, storing the compressed output to
`example.astc`.

## Decompressing an image

Decompress an image using the `-d` option:

    astcenc -d example.astc example.tga

This decompresses `example.astc` storing the decompressed output to
`example.tga`.

## Measuring image quality

Review the compression quality using the `-t` option:

    astcenc -t example.png example.tga

This is equivalent to compressing and then immediately decompressing the
image, allowing a visual inspection of the decompression quality. In addition
this mode also prints out the PSNR quality of the compressed image to the
console.

## Experimenting

Efficient real-time graphics benefits from low bandwidth texture access,
as it reduces bandwidth, saves energy, and can improve texture cache
efficiency. However, like any lossy compression format there will come a point
where the output compression quality is unacceptable because there are simply
not enough bits to represent the output at the quality needed. We recommend
experimenting with the block footprint to find the optimum balance between size
and quality, as the finely adjustable compression ratio is one of major
strengths of the ASTC format.

The compression speed can be controlled from `-veryfast`, through `-fast`,
`-medium` and `-thorough`, up to `-exhaustive`. In general, the more time the
encoder has to spend looking for good encodings the better the results, but it
does result in increasingly small improvements for the amount of time required.

There are many other command line options for tuning the encoder parameters
which can be used to fine tune the compression algorithm. See the command line
help message for more details.

# Support

If you have issues with the `astcenc` encoder, or questions about the ASTC
texture format itself, please raise them in the GitHub issue tracker.

If you have any questions about Arm Mali GPUs or other graphics technology
please submit your questions to the [Arm Mali Graphics
forums](https://community.arm.com/graphics/).

- - -

_Copyright (c) 2013-2019, Arm Limited and Contributors. All rights reserved._
