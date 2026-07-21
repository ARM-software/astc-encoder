// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2011-2026 Arm Limited
//
// Licensed under the Apache License, Version 2.0 (the "License"); you may not
// use this file except in compliance with the License. You may obtain a copy
// of the License at:
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
// WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
// License for the specific language governing permissions and limitations
// under the License.
// ----------------------------------------------------------------------------

/**
 * @brief Functions for printing build info and help messages.
 */

#include "astcenccli_internal.h"
#include "astcenccli_version.h"

/** @brief The version header. */
static const char *astcenc_copyright_string =
R"(astcenc v%s, %u-bit %s%s%s
Copyright (c) 2011-%s Arm Limited. All rights reserved.
)";

/** @brief The short-form help text. */
static const char *astcenc_short_help =
R"(
Basic usage:

To compress an image use:
    astcenc {-cl|-cs|-ch|-cH} <in> <out> <blockdim> <quality> [options]

e.g. using LDR profile, 8x6 blocks, and the thorough quality preset:
    astcenc -cl kodim01.png kodim01.astc 8x6 -thorough

To decompress an image use:
    astcenc {-dl|-ds|-dh|-dH} <in> <out>

e.g. using LDR profile:
    astcenc -dl kodim01.astc kodim01.png

To perform a compression test, writing back the decompressed output, use:
    astcenc {-tl|-ts|-th|-tH} <in> <out> <blockdim> <quality> [options]

e.g. using LDR profile, 8x6 blocks, and the thorough quality preset:
    astcenc -tl kodim01.png kodim01-test.png 8x6 -thorough

The -*l options are used to configure the codec to support only the linear
LDR profile, preventing use of the HDR encoding features.

The -*s options are used to configure the codec to support only
the sRGB LDR profile, preventing use of the HDR encoding features. Input
texture data must be encoded in the sRGB color space for this option to
provide correct output results.

The -*h/-*H options are used to configure the codec to support the HDR ASTC
color profile. Textures compressed with this profile may fail to decompress
correctly on GPU hardware without HDR profile support. The -*h options
configure the compressor for HDR RGB components and an LDR alpha component.
The -*H options configure the compressor for HDR across all 4 components.

For full help documentation run 'astcenc -help'.
)";

/** @brief The long-form help text. */
static const char* astcenc_long_help = R"(
NAME

    astcenc - compress or decompress images using the ASTC format.

SYNOPSIS

    astcenc {-h|-help}
    astcenc {-v|-version}
    astcenc {-cl|-cs|-ch|-cH} <in> <out> <blocksize> <quality> [options]
    astcenc {-dl|-ds|-dh|-dH} <in> <out>
    astcenc {-tl|-ts|-th|-tH} <in> <out> <blocksize> <quality> [options]

DESCRIPTION

    astcenc compresses image files into the Adaptive Scalable Texture
    Compression (ASTC) data format, a lossy compression format designed for
    use in real-time graphics applications. astcenc supports all of the
    compression profiles and block sizes allowed by the ASTC format:

    * All color profiles (LDR linear, LDR sRGB, and HDR)
    * All 2D block sizes (4x4 through to 12x12)
    * All 3D block sizes (3x3x3 through to 6x6x6)

    astcenc compression converts an input image into an ASTC image. It
    provides a flexible tradeoff between compression time and image quality
    using a general quality setting as well as many advanced control
    options for fine tuning.

    astcenc decompression converts an ASTC image into an output image.

    astcenc test mode can perform image quality analysis for a given
    set of compression options.

COMPRESSION

    To compress an image using astcenc you must specify the color profile,
    the input file name, the output file name, the target block size, and
    the quality level.

    The color profile is specified using the -cl (LDR linear), -cs (LDR
    sRGB), -ch (HDR RGB, LDR A), or -cH (HDR RGBA) encoder options. GPUs
    that support ASTC must support both LDR profiles, but support for the
    HDR profile is optional.

    The input and output file paths and extensions must specify a valid
    file for the compression operation. See the COMPRESSION FILE FORMATS
    section for the list of supported formats.

    The block size must be a valid ASTC block size. Every block compresses
    into 128 bits of compressed output, so the block size determines the
    compressed data bitrate.

    Supported 2D block sizes are:

    * 4x4: 8.00 bpp
    * 5x4: 6.40 bpp
    * 5x5: 5.12 bpp
    * 6x5: 4.27 bpp
    * 6x6: 3.56 bpp
    * 8x5: 3.20 bpp
    * 8x6: 2.67 bpp
    * 10x5: 2.56 bpp
    * 10x6: 2.13 bpp
    * 8x8: 2.00 bpp
    * 10x8: 1.60 bpp
    * 10x10: 1.28 bpp
    * 12x10: 1.07 bpp
    * 12x12: 0.89 bpp

    Supported 3D block sizes are:

    * 3x3x3: 4.74 bpp
    * 4x3x3: 3.56 bpp
    * 4x4x3: 2.67 bpp
    * 4x4x4: 2.00 bpp
    * 5x4x4: 1.60 bpp
    * 5x5x4: 1.28 bpp
    * 5x5x5: 1.02 bpp
    * 6x5x5: 0.85 bpp
    * 6x6x5: 0.71 bpp
    * 6x6x6: 0.59 bpp

    The quality level sets the quality-performance tradeoff for the
    compressor. Higher quality levels make more complete searches of the
    search space to improve image quality, at the expense of compression
    time. The quality level can be set to any value between 0 (fastest) and
    100 (exhaustive), or to a pre-defined quality preset:

        -fastest        (quality 0)
        -fast           (quality 10)
        -medium         (quality 60)
        -thorough       (quality 98)
        -verythorough   (quality 99)
        -exhaustive     (quality 100)

    For compression of production content we recommend using the
    -thorough quality level. Higher quality levels are substantially
    slower with only minor improvements to image quality.

    There are a number of additional compressor options which are useful
    to consider, based on the type of image data being compressed.

    -decode_unorm8
        Set this for textures that use the decode_unorm8 extension
        decompression behavior at run-time. This improves image quality
        because the compression search will correctly round values to match
        the run-time behavior.

    -a <radius>
        Set this for color textures with an alpha translucency to scale
        per-texel weights by the alpha value. The alpha value used to scale
        the color weights is the average of all texels within the <radius>.
        Set <radius> to 0 to use only a texel's own alpha value, but note
        that this can cause issues when sampling with filtering.

        This improves image quality for color textures with alpha
        translucency because the compressor does not try to preserve color
        accuracy that would not be visible after quantization.

        Do not use this option when your textures store color values that
        have been premultiplied by the alpha value, as the equivalent
        scaling has already been applied.

    -normal
        Set this for normal map textures when the input contains
        three-component unit-length normals stored as (R=X, G=Y, B=Z). The
        output texture will be a two-component XY normal map stored as
        (RGB=X, A=Y). The Z component can be computed in shader code using
        the equation:

            nml.xy = texture(...).ga;                // Load in [0,1]
            nml.xy = nml.xy * 2.0 - 1.0;             // Unpack to [-1,1]
            nml.z = sqrt(1 - dot(nml.xy, nml.xy));   // Compute Z

        This improves image quality because the compressor only has to
        store two color endpoints per block.

        Alternative compressed component orders can be set by using the
        -esw and -dsw options after the -normal option on the command line.

    -rgbm <max>
        Set this for RGBM encoded textures, which store HDR values between
        0 and <max> in an LDR container format with a shared multiplier.
        The HDR value can be computed in shader code using the equation:

            vec3 hdr_value = tex.rgb * tex.a * max;

        This improves image quality because astcenc can compute error in
        the reconstructed HDR values when choosing the best encoding.

        To avoid severe block artifacts in the compressed image, we
        recommend that your HDR to RGBM conversion preprocess limits the
        value of M to keep it above a minimum threshold. We recommend
        using at least 16/255 as your lower limit.

    -perceptual
        Set this for LDR color textures to optimize for perceptual error
        instead of simple root-mean-square error.

        This improves perceived image quality in color data, but typically
        lowers standard metrics such as PSNR. Perceptual methods should not
        be used for non-color data.

    -zdim <zdim>
        Set this for 3D volumetric textures that are loaded from a set of
        <zdim> 2D input image files. The input filename given is decorated
        with the postfix "_<slice>" to find the file to load. For example,
        specifying the input file as input.png would load input_0.png,
        input_1.png, etc.

    -pp-normalize
        Set this for a normal map texture to run a preprocess over the
        image that forces normal vectors to be unit length. The
        preprocessing is applied before any codec encoding swizzle, so data
        must be in the RGB components of the input image.

    -pp-premultiply
        Set this for a color texture to run a preprocess over the image
        that scales RGB components by the alpha component value. The
        preprocessing is applied before any codec encoding swizzle, so
        color data must be in the RGB components of the input image and
        alpha data must be in the A component.

    -yflip
        Flip the image in the vertical axis prior to compression and after
        decompression. Note that using this option in test mode (-t*) will
        have no effect because the image will be flipped twice.

    -j <threads>
        Explicitly specify the number of threads to use in the codec. If
        not specified, the codec will use one thread per CPU detected in
        the system.)"
// This split in the literals is needed for Visual Studio; the compiler
// will concatenate these two strings together ...
R"(

ADVANCED COMPRESSION

    Options in this section provide extensive control over compressor
    settings and heuristics.

    Error weighting options
    -----------------------

    These options provide low-level control of the codec error metric
    computation, which is used to determine which compression trials
    produce the best result.

    -cw <red> <green> <blue> <alpha>
        Set this option to bias the error significance of each color
        component. Set a component error weight to zero to exclude it
        completely, which is useful if the component is unused.

        This improves the quality of the color components with a higher
        relative error weight, but can reduce the quality of the other
        components.

    -mpsnr <low> <high>
        Set this option for an HDR texture to set the low and high f-stop
        values for the mPSNR error metric.

        This option does not impact the compressed image quality, but
        changes how the mPSNR error metric is computed in test mode.

    Performance-quality tradeoff options
    ------------------------------------

    These options provide low-level control of the individual heuristics
    that determine the performance-quality tradeoff. Heuristic settings
    default to values that vary based on the quality level and the targeted
    bitrate. The values documented here are provided as an indication of
    recommended settings and are taken from the high bitrate set.

    -dblimit <number>
        Stop compression work on a block when the PSNR of the block,
        measured in dB, exceeds <number>. This option is ineffective for
        HDR textures. Preset defaults, where N is the number of texels in a
        block, are:

            -fastest      : MAX(63-19*log10(N),  85-35*log10(N))
            -fast         : MAX(63-19*log10(N),  85-35*log10(N))
            -medium       : MAX(70-19*log10(N),  95-35*log10(N))
            -thorough     : MAX(77-19*log10(N), 105-35*log10(N))
            -verythorough : 999
            -exhaustive   : 999

    -partitioncountlimit <number>
        Test up to <number> partitions for each block. Preset defaults are:

            -fastest      : 2
            -fast         : 3
            -medium       : 4
            -thorough     : 4
            -verythorough : 4
            -exhaustive   : 4

    -[2|3|4]partitionindexlimit <number>
        Estimate up to <number> partition indices for this partition count.
        Preset defaults are:

            -fastest      :   10 |   6 |   4
            -fast         :   18 |  10 |   8
            -medium       :   34 |  28 |  16
            -thorough     :   82 |  60 |  30
            -verythorough :  256 | 128 |  64
            -exhaustive   :  512 | 512 | 512

    -[2|3|4]partitioncandidatelimit <number>
        Test up to <number> partition indices for this partition count.
        Preset defaults are:

            -fastest      :   2 |  2 |  2
            -fast         :   2 |  2 |  2
            -medium       :   2 |  2 |  2
            -thorough     :   3 |  2 |  2
            -verythorough :  20 | 14 |  8
            -exhaustive   :  32 | 32 | 32

    -blockmodelimit <centile>
        Test up to <centile> block modes selected from a precomputed
        frequency distribution. Higher numbers give better quality, but
        increase compression time. This option is ineffective for 3D
        textures. Preset defaults are:

            -fastest      :  43
            -fast         :  55
            -medium       :  77
            -thorough     :  94
            -verythorough :  98
            -exhaustive   : 100

    -candidatelimit <number>
        Test up to <number> candidate encodings for each block mode:

            -fastest      : 2
            -fast         : 3
            -medium       : 3
            -thorough     : 4
            -verythorough : 6
            -exhaustive   : 8

    -refinementlimit <number>
        Perform up to <number> refinement iterations on colors and
        weights. Preset defaults are:

            -fastest      : 2
            -fast         : 3
            -medium       : 3
            -thorough     : 4
            -verythorough : 4
            -exhaustive   : 4

    -[2|3]partitionlimitfactor <factor>
        Stop compression before testing N+1 partitions unless the N
        partition error is lower than the error from encoding with N-1
        partitions by more than the specified factor. Preset defaults are:

            -fastest       : 1.00 | 1.00
            -fast          : 1.00 | 1.00
            -medium        : 1.10 | 1.05
            -thorough      : 1.35 | 1.15
            -verythorough  : 1.60 | 1.40
            -exhaustive    : 2.00 | 2.00

    -2planelimitcorrelation <factor>
        Stop compression after testing only one plane of weights, unless
        the minimum color correlation factor between any pair of color
        components is below this factor. This option is ineffective for
        normal maps. Preset defaults are:

            -fastest      : 0.50
            -fast         : 0.65
            -medium       : 0.85
            -thorough     : 0.95
            -verythorough : 0.98
            -exhaustive   : 0.99
)"
// This split in the literals is needed for Visual Studio; the compiler
// will concatenate these two strings together ...
R"(
    Other options
    -------------

    -esw <swizzle>
        Specify an encoding swizzle to reorder the color components before
        compression. The swizzle is specified using a four character
        string, which defines the compression component load order.

        The characters may be taken from the set [rgba01], selecting
        from the input color components or a constant 0 or 1. For example,
        to swap the RG components and replace alpha with a constant one
        the swizzle 'grb1' should be used.

        By default all 4 post-swizzle components are included in trial
        error metrics. When using -esw to map two component data to the L+A
        endpoint (e.g., -esw rrrg) the luminance data stored in the RGB
        components will be weighted 3 times more strongly than the alpha
        component. This can be corrected using -ssw to specify which
        components are sampled at run-time.

    -ssw <swizzle>
        Specify a sampling swizzle to adjust component error weights based
        on which components are actually read by your shader program. For
        example, using -ssw ra tells the compressor that the green and blue
        error can be ignored because the data is never read.

        The sampling swizzle is based on the channel ordering after the
        -esw transform has been applied. Using -ssw is equivalent to using
        -cw to set unused channel error weights to zero.

    -dsw <swizzle>
        Specify a decompression swizzle to reorder the color components
        after decompression. The swizzle is specified using a four
        character string, which defines the component store order used by
        the decompressor.

        The characters may be taken from the set [rgbaz01], selecting
        from the input color components, a reconstructed Z value, or a
        constant 0 or 1.

        The z character is used to specify that the compressed data
        stores an XY normal map, and that the Z component must be
        reconstructed from the two components stored in the data.

    -silent
        Suppress all non-essential diagnostic output from the codec.
)"
// This split in the literals is needed for Visual Studio; the compiler
// will concatenate these two strings together ...
R"(

DECOMPRESSION

    To decompress an image stored in the ASTC format you must specify the
    color profile, the input file name, and the output file name.

    The color profile is specified using the -dl (LDR linear), -ds (LDR
    sRGB), -dh (HDR RGB, LDR A), or -dH (HDR RGBA) decoder options.

    The input and output file paths and extensions must specify a valid
    file for the decompression operation. See the DECOMPRESSION FILE
    FORMATS section for the list of supported formats.

    The -dsw option documented in ADVANCED COMPRESSION is also relevant to
    decompression.

TEST

    To perform a compression test of an image using astcenc you must
    specify the color profile, the input file name, the output file name,
    the target block size, and the quality level. This mode will round-trip
    the image through compression and decompression, store the decompressed
    result, and emit image quality and performance metrics.

    The color profile is specified using the -tl (LDR linear), -ts (LDR
    sRGB), -th (HDR RGB, LDR A), or -tH (HDR RGBA) encoder options.

    The input file path and extension must specify a valid file for the
    compression operation. See the COMPRESSION FILE FORMATS section for the
    list of supported formats. The output file path and extension must
    specify a valid file for the decompression operation. See the
    DECOMPRESSION FILE FORMATS section for the list of supported formats.

COMPRESSION FILE FORMATS

    The following formats are supported as compression inputs:

    * LDR Formats:
        * BMP (*.bmp)
        * PNG (*.png)
        * Targa (*.tga)
        * JPEG (*.jpg)
    * HDR Formats:
        * OpenEXR (*.exr)
        * Radiance HDR (*.hdr)
    * Container Formats:
        * Khronos Texture KTX (*.ktx)
        * DirectDraw Surface DDS (*.dds)

    For the KTX and DDS formats only a subset of the features are
    supported:

    * Texel format must be R, RG, RGB, BGR, RGBA, BGRA, L, or LA.
    * Texture topology must be 2D, 2D-array, 3D, or cube-map.
    * Only the first mipmap level, cube face, or 2D array layer will be
      read.

    The following formats are supported as compression outputs:

    * ASTC (*.astc)
    * Khronos Texture KTX (*.ktx)

DECOMPRESSION FILE FORMATS

    The following formats are supported as decompression inputs:

    * ASTC (*.astc)
    * Khronos Texture KTX (*.ktx)

    The following formats are supported as decompression outputs:

    * LDR Formats:
        * BMP (*.bmp)
        * PNG (*.png)
        * Targa (*.tga)
    * HDR Formats:
        * OpenEXR (*.exr)
        * Radiance HDR (*.hdr)
    * Container Formats:
        * Khronos Texture KTX (*.ktx)
        * DirectDraw Surface DDS (*.dds)

QUICK REFERENCE

    To compress an image use:
        astcenc {-cl|-cs|-ch|-cH} <in> <out> <blocksize> <quality> [options]

    To decompress an image use:
        astcenc {-dl|-ds|-dh|-dH} <in> <out>

    To perform a quality test use:
        astcenc {-tl|-ts|-th|-tH} <in> <out> <blocksize> <quality> [options]

    Mode -*l = linear LDR, -*s = sRGB LDR, -*h = HDR RGB/LDR A, -*H = HDR.
    Quality = -fastest/-fast/-medium/-thorough/-verythorough/-exhaustive
              or a floating-point value between 0 and 100.
)";

/* See header for documentation. */
void astcenc_print_header()
{
#if (ASTCENC_AVX == 2)
	const char* simdtype = "avx2";
#elif (ASTCENC_SSE == 41)
	const char* simdtype = "sse4.1";
#elif (ASTCENC_SSE == 20)
	const char* simdtype = "sse2";
#elif (ASTCENC_SVE == 8)
	const char* simdtype = "sve.256b";
#elif (ASTCENC_SVE == 4)
	const char* simdtype = "sve.128b";
#elif (ASTCENC_NEON == 1)
	const char* simdtype = "neon";
#else
	const char* simdtype = "none";
#endif

#if (ASTCENC_POPCNT == 1)
	const char* pcnttype = "+popcnt";
#else
	const char* pcnttype = "";
#endif

#if (ASTCENC_F16C == 1)
	const char* f16ctype = "+f16c";
#else
	const char* f16ctype = "";
#endif

	unsigned int bits = static_cast<unsigned int>(sizeof(void*) * 8);
	printf(astcenc_copyright_string,
	       VERSION_STRING, bits, simdtype, pcnttype, f16ctype, YEAR_STRING);

    // If possible, print hint that 8-wide SVE could be used
#if ASTCENC_SVE == 4
    if (svcntw() == 8)
    {
        printf("Note: This CPU can support 256-bit SVE builds.\n");
    }
#endif
}

/* See header for documentation. */
void astcenc_print_shorthelp()
{
	astcenc_print_header();
	printf("%s", astcenc_short_help);
}

/* See header for documentation. */
void astcenc_print_longhelp()
{
	astcenc_print_header();
	printf("%s", astcenc_long_help);
}
