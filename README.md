# About

This is the official repository for the Arm® Adaptive Scalable Texture
Compression (ASTC) Encoder, `astcenc`, a command-line tool for compressing
and decompressing images using the ASTC texture compression standard.

## The ASTC format

The ASTC compressed data format, developed by Arm® and AMD, has been adopted as
an official extension to the Open GL®, OpenGL ES, and Vulkan® graphics APIs. It
provides a major step forward both in terms of image quality at a given
bitrate, and in terms of the format and bitrate flexibility available to
content creators. This allows more assets to use compression, often at a
reduced bitrate compared to legacy formats, reducing memory bandwidth and
energy consumption.

Read our [ASTC Format Overview](./Docs/FormatOverview.md) for a quick
introduction, or read the full data format specification:

* [Khronos Data Format Specification v1.2 # ASTC](https://www.khronos.org/registry/DataFormat/specs/1.2/dataformat.1.2.html#ASTC)

## License

This project is licensed under the Apache 2.0 license. By downloading any
component from this repository you acknowledge that you accept terms specified
in the [LICENSE.txt](LICENSE.txt) file.

# Branches

The `master` branch is an active development branch for the next major release
of the compressor; version 2.x. It aims to be a stable branch, but as it is
still under development expect changes to both the command line and the
quality-performance trade offs the compressor is making.

The `1.x` branch is a maintenance branch for the 1.x release series. It is
stable and we will now only land bug fixes for this branch; no new
functionality or performance improvements should be expected.

# Encoder feature support

The encoder supports compression of low dynamic range (BMP, JPEG, PNG, TGA) and
high dynamic range (EXR, HDR) images, as well as a subset of image data wrapped
in the DDS and KTX container formats, into ASTC or KTX format output images.

The decoder supports decompression of ASTC or KTX format input images into low
dynamic range (BMP, PNG, TGA), high dynamic range (EXR, HDR), or DDS and KTX
wrapped output images.

The encoder allows control over the compression time/quality tradeoff with
`exhaustive`, `thorough`, `medium`, and `fast` encoding quality presets.

The encoder allows compression time and quality analysis by reporting the
compression time, and the Peak Signal-to-Noise Ratio (PSNR) between the input
image and the compressed output.

## ASTC format support

The ASTC specification allows three profiles of implementation:

* 2D Low Dynamic Range (LDR profile)
* 2D LDR and High Dynamic Range (HDR profile)
* 2D and 3D, LDR and HDR (Full profile)

The `astcenc` compressor supports generation of images for all three profiles.
In addition it also supports all of the ASTC block sizes and compression
modes, allowing content creators access the full spectrum of quality-to-bitrate
options ranging from 0.89 bits/pixel up to 8 bits/pixel.

# Prebuilt binaries

Prebuilt release build binaries of `astcenc` for 64-bit Linux, macOS, and
Windows are available in the
[GitHub Releases page](https://github.com/ARM-software/astc-encoder/releases).
Note that currently no 2.x series pre-built binaries are available.

# Getting started

Open a terminal, change to the appropriate directory for your system, and run
the astcenc encoder program, like this on Linux or macOS:

    ./astcenc

... or like this on Windows:

    astcenc

Invoking `astcenc -help` gives an extensive help message, including usage
instructions and details of all available command line options. A summary of
the main encoder options are shown below.

## Compressing an image

Compress an image using the `-cl` \ `-cs` \ `-ch` \ `-cH` modes. For example:

    astcenc -cl example.png example.astc 6x6 -medium

This compresses `example.png` using the LDR color profile and a 6x6 block
footprint (3.55 bits/pixel). The `-medium` quality preset gives a reasonable
image quality for a relatively fast compression speed. The output is stored to
a linear color space compressed image, `example.astc`.

The modes available are:

* `-cl` : use the linear LDR color profile.
* `-cs` : use the sRGB LDR color profile.
* `-ch` : use the HDR color profile, tuned for HDR RGB and LDR A.
* `-cH` : use the HDR color profile, tuned for HDR RGBA.

## Decompressing an image

Decompress an image using the `-dl` \ `-ds` \ `-dh` \ `-dH` modes. For
example:

    astcenc -dh example.astc example.tga

This decompresses `example.astc` using the full HDR feature profile, storing
the decompressed output to `example.tga`.

The modes available are:

* `-dl` : use the linear LDR color profile.
* `-ds` : use the sRGB LDR color profile.
* `-dh` and `-dH` : use the HDR color profile.

Note that for decompression there is no difference between the two HDR modes,
they are both provided simply to maintain symmetry across operations.

## Measuring image quality

Review the compression quality using the `-tl` \ `-ts` \ -`th` \ -`tH` modes.
For example:

    astcenc -tl example.png example.tga 5x5 -thorough

This is equivalent to using using the LDR color profile and a 5x5 block size
to compress the image, using the `-thorough` quality preset, and then
immediately decompressing the image and saving the result. This can be used
to enable a visual inspection of the compressed image quality. In addition
this mode also prints out some image quality metrics to the console.

The modes available are:

* `-tl` : use the linear LDR color profile.
* `-ts` : use the sRGB LDR color profile.
* `-th` : use the HDR color profile, tuned for HDR RGB and LDR A.
* `-tH` : use the HDR color profile, tuned for HDR RGBA.

## Experimenting

Efficient real-time graphics benefits from minimizing compressed texture size,
as it reduces memory bandwidth, saves energy, and can improve texture cache
efficiency. However, like any lossy compression format there will come a point
where the compressed image quality is unacceptable because there are simply
not enough bits to represent the output with the precision needed. We
recommend experimenting with the block footprint to find the optimum balance
between size and quality, as the finely adjustable compression ratio is one of
major strengths of the ASTC format.

The compression speed can be controlled from `-fast`, through `-medium` and
`-thorough`, up to `-exhaustive`. In general, the more time the encoder has to
spend looking for good encodings the better the results, but it does result in
increasingly small improvements for the amount of time required.

There are many other command line options for tuning the encoder parameters
which can be used to fine tune the compression algorithm. See the command line
help message for more details.

# Documentation

The [ASTC Format Overview](./Docs/FormatOverview.md) page provides a high level
overview of the ASTC data format, how it encodes data, and why it is both
flexible and efficient.

The [Effective ASTC Encoding](./Docs/Encoding.md) page looks at some of the
guidelines that should be followed when compressing data using `astcenc`.
It covers:

* How to efficiently encode data with fewer than 4 channels.
* How to efficiently encode normal maps, sRGB data, and HDR data.
* Coding equivalents to other compression formats.

The [Building ASTC Encoder](./Docs/Building.md) page provides instructions on
how to build `astcenc` from the sources in this repository.

The [Testing ASTC Encoder](./Docs/Testing.md) page provides instructions on
how to test any modifications to the source code in this repository.

# Support

If you have issues with the `astcenc` encoder, or questions about the ASTC
texture format itself, please raise them in the GitHub issue tracker.

If you have any questions about Arm Mali GPUs, application development for Arm
Mali GPUs, or general graphics technology please submit them on the [Arm Mali
Graphics forums](https://community.arm.com/graphics/).

- - -

_Copyright (c) 2013-2020, Arm Limited and contributors. All rights reserved._
