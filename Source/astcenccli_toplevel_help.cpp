// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2011-2020 Arm Limited
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

const char *astcenc_copyright_string =
R"(ASTC codec v2.0.alpha, %u-bit %s%s
Copyright 2011-2020 Arm Limited, all rights reserved
)";


const char *astcenc_short_help =
R"(
Basic usage:

To compress an image use the following command line:
    astcenc {-cl|-cs|-ch|-cH} <in> <out> <blockdim> <preset> [options]

For example, to compress to 8x6 blocks with the thorough preset use:
    astcenc -cl kodim01.png kodim01.astc 8x6 -thorough

To decompress an image use the following command line:
    astcenc {-dl|-ds|-dh|-dH} <in> <out>

For example, use:
    astcenc -dl kodim01.astc kodim01.png

To perform a compression test, writing back the decompressed output, use
the following command line:
    astcenc {-tl|-ts|-th|-tH} <in> <out> <blockdim> <preset> [options]

For example, use:
    astcenc -tl kodim01.png kodim01-test.png 8x6 -thorough

The -*l options are used to configure the codec to support only the linear
LDR profile, preventing use of the HDR encoding features.

The -*s options are used to configure the codec to support only
the sRGB LDR profile, preventing use of the HDR encoding features. Input
texture data must be encoded in the sRGB colorspace for this option to
provide correct output results.

The -*h/-*H options are used to configure the codec to support the HDR ASTC
color profile. Textures compressed with this profile may fail to decompress
correctly on GPU hardware without HDR profile support. The -*h options
configure the compressor for HDR RGB channels with an LDR alpha channel.
The -*H options configure the compressor for full HDR across all channels.

For full help documentation run 'astcenc -help'.
)";


const char *astcenc_long_help = R"(
NAME
       astcenc - compress or decompress images using the ASTC format

SYNOPSIS
       astcenc {-h|-help}
       astcenc {-v|-version}
       astcenc {-cl|-cs|-ch|-cH} <in> <out> <blocksize> <preset> [options]
       astcenc {-dl|-ds|-dh|-dH} <in> <out> <blocksize> <preset> [options]
       astcenc {-tl|-ts|-th|-tH} <in> <out> <blocksize> <preset> [options]

DESCRIPTION
       astcenc compresses image files into the Adaptive Scalable Texture
       Compression (ASTC) image format, a lossy compression format design
       for use in real-time graphics applications. It is a fully featured
       compressor implementation, supporting all of the compression
       profiles and block sizes specified by the ASTC format:

           All color profiles (LDR linear, LDR sRGB, and HDR)
           All 2D block sizes (4x4 though to 12x12)
           All 3D block sizes (3x3x3 through to 6x6x6)

       The compressor provides a number of pre-determined quality presets,
       which allow users to tradeoff compressed image quality against
       compression performance. For advanced users the compressor provides
       many additional control options.

       astcenc can also be used to decompress ASTC compressed images, and
       perform compression image quality analysis.

COMPRESSION
       To compress an image using the ASTC format you must specify the
       color profile, the input file name, the output file name, the
       target block size, and the quality preset.

       The color profile is specified using the -cl (LDR linear),
       -cs (LDR sRGB), -ch (HDR RGB, LDR A), or -cH (HDR RGBA) encoder
       options. Note that not all hardware implementations of ASTC support
       the HDR profile.

       The input file path must match a valid file format for compression,
       and the output file format must be a valid output for compression.
       See the FILE FORMATS section for the list of supported formats.

       The block size must be a valid ASTC block size. Every block
       compresses into 128 bits of compressed output, so the block size
       determines the compressed data bitrate.

       Supported 2D block sizes are:

             4x4: 8.00 bpp        10x5: 2.56 bpp
             5x4: 6.40 bpp        10x6: 2.13 bpp
             5x5: 5.12 bpp         8x8: 2.00 bpp
             6x5: 4.27 bpp        10x8: 1.60 bpp
             6x6: 3.56 bpp       10x10: 1.28 bpp
             8x5: 3.20 bpp       12x10: 1.07 bpp
             8x6: 2.67 bpp       12x12: 0.89 bpp

       Supported 3D block sizes are:

           3x3x3: 4.74 bpp       5x5x4: 1.28 bpp
           4x3x3: 3.56 bpp       5x5x5: 1.02 bpp
           4x4x3: 2.67 bpp       6x5x5: 0.85 bpp
           4x4x4: 2.00 bpp       6x6x5: 0.71 bpp
           5x4x4: 1.60 bpp       6x6x6: 0.59 bpp

       The quality preset configures the quality-performance tradeoff for
       the compressor; more complete searches of the search space improve
       image quality at the expense of compression time. The available
       presets are:

           -fast
           -medium
           -thorough
           -exhaustive

       Note that using -exhaustive significantly increases compression
       time, but typically only gives minor quality improvements over
       using -thorough.

       There are a number of additional compressor options which are
       useful to consider for common usage, based on the type of image
       data being compressed.

       -mask
           The input texture is a mask texture with unrelated data stored
           in the various color channels. This improves image quality by
           trying to minimize the effect of error cross-talk across the
           color channels.

       -normal
           The input texture is a three channel normal map, storing unit
           length normals as R=X, G=Y, and B=Z, optimized for angular PSNR.
           The compressor will compress this data as a two channel X+Y
           normal map with the following channel layout (RGB=X, A=Y). The
           Z component can be recovered programatically in shader core by
           using the equation:

               Z = sqrt(1 - X^2 - Y^2).

       -perceptual
           The codec should optimize perceptual error, instead of direct
           RMS error. This aims to improves perceived image quality, but
           typically lowers the measured PSNR score. Perceptual methods
           are currently only available for normal maps.

       -array <size>
           Loads an array of <size> 2D image slices to use as a 3D image.
           The input filename given is used is decorated with the postfix
           "_<slice>" to find the file to load. For example, an input
           named "input.png" would load as input_0.png, input_1.png, etc.

COMPRESSION TIPS & TRICKS
       ASTC is a block-based format that can be prone to block artifacts.
       If block artifacts are a problem when compressing a given texture,
       adding some or all of following command-line options may help:

           -b 1.8
           -v 2 1 1 0 25 0.1
           -va 1 1 0 25
           -dblimit 60

       The -b option is a general-purpose block-artifact reduction option.
       The -v and -va option settings will concentrate effort where smooth
       regions lie next to regions with high detail, which are particularly
       prone to block artifacts. Increasing the -dblimit option is
       sometimes also needed to force the compressor to keep searching for
       a better encoding, which can be needed in images with smooth
       gradients.

       If a texture exhibits severe block artifacts in only some of the
       color channels, which is a common problem for mask textures, then
       using the -cw option to raise the weighting of the affected color
       channel(s) may help. For example, if the green color channel is
       particularly badly encoded then try '-cw 1 6 1 1'.

ADVANCED COMPRESSION
       Error weighting options
       -----------------------

       These options provide low-level control of the codec error metric
       computation, used to determine what good compression looks like.

       -v <radius> <power> <base> <mean> <stdev> <mix>
           Compute the per-texel relative error weighting for the RGB color
           channels as follows:

           weight = 1 / (<base> + <mean> * mean^2 + <stdev> * stdev^2)

           The <radius> argument specifies the texel radius of the
           neighborhood over which the average and standard deviation are
           computed.

           The <mix> parameter is used to control the degree of mixing
           of the average and stddev error components across the color
           channels. Setting this parameter to 0 causes the computation
           to be done completely separately for each color channel; setting
           it to 1 causes the results from the RGB channels to be combined
           and applied to all three together. Intermediate values between
           these two extremes do a linear mix of the two error values.

           The <power> argument is a power used to raise the values of the
           input texels before computing average and standard deviation;
           e.g. a power of 0.5 causes the codec to take the square root
           of every input texel value.

       -va <power> <base> <mean> <stdev>
           Compute the per-texel relative error weighting for the alpha
           channel, when used in conjunction with -v. See documentation for
           -v for parameter documentation.

       -a <radius>
           For textures with alpha channel, scale per-texel weights by the
           alpha value. The alpha value chosen for scaling of any
           particular texel is taken as an average across a neighborhood
           of the texel defined by the <radius> argument. Setting <radius>
           to 0 causes only the texel's own alpha to be used.

       -cw <red> <green> <blue> <alpha>
           Assign an additional weight scaling to each color channel,
           allowing the channels to be treated differently in terms of
           error significance. Set values above 1 to increase a channel's
           significance, and values below 1 to decrease it. Set to 0
           to exclude a channel from error computation completely.

       -b <weight>
           Assign an additional weight scaling for texels at compression
           block edges and corners. Setting this to a value above 1
           increases the significance of texels closer to the edges of a
           block, and can help to reduce block artifacts.

       -mpsnr <low> <high>
           Set the low and high f-stop values for the mPSNR error metric.
           The mPSNR error metric only applies to HDR textures.

       Performance-quality tradeoff options
       ------------------------------------

       These options provide low-level control of the codec heuristics that
       drive the performance-quality trade off.

       -partitionlimit <number>
           Test only <number> block partitions. Higher numbers give
           better quality, however large values give diminishing returns
           especially for smaller block sizes. Preset defaults are:

               -fast       :    4
               -medium     :   25
               -thorough   :  100
               -exhaustive : 1024

       -blockmodelimit <number>
           Test only block modes below the <number> usage centile in an
           empirically determined distribution of block mode frequency.
           This option is ineffective for 3D textures. Preset defaults are:

               -fast       :  50
               -medium     :  75
               -thorough   :  95
               -exhaustive : 100

       -refinementlimit <value>
           Iterate only <value> refinement iterations on colors and
           weights. Minimum value is 1. Preset defaults are:

               -fast       : 1
               -medium     : 2
               -thorough   : 4
               -exhaustive : 4

       -dblimit <number>
           Stop compression work on a block as soon as the PSNR of the
           block, measured in dB, exceeds <number>. This option is
           ineffective for HDR textures. Preset defaults, where N is the
           number of texels in a block, are:

               -fast       : dblimit = MAX(63-19*log10(N),  85-35*log10(N))
               -medium     : dblimit = MAX(70-19*log10(N),  95-35*log10(N))
               -thorough   : dblimit = MAX(77-19*log10(N), 105-35*log10(N))
               -exhaustive : dblimit = 999

       -partitionearlylimit <factor>
           Stop compression work on a block after only testing blocks with
           up to two partions and one plane of weights, unless the two
           partition error term is lower than the error term from encoding
           with one partition by more than the specified factor. This
           option is ineffective for normal maps. Preset defaults are:

               -fast       :    1.0
               -medium     :    1.2
               -thorough   :    2.5
               -exhaustive : 1000.0

       -planecorlimit <factor>
           Stop compression after testing only one planes of weights,
           unless the minimum color correlation factor between any pair of
           color channels is below this factor. This option is ineffective
           for normal maps. Preset defaults are:

               -fast       : 0.50
               -medium     : 0.75
               -thorough   : 0.95
               -exhaustive : 0.99

       Other options
       -------------

       -esw <swizzle>
           Swizzle the color components before compression. The swizzle is
           specified using a 4-character string, which defines the output
           format ordering. The characters may be taken from the set
           [rgba01], selecting either input color channels or a literal
           zero or one. For example to swap the RG channels, and replace
           alpha with 1, the swizzle 'grb1' should be used.

           Note that the input swizzle is assumed to take place before any
           compression, and all error weighting applies to the post-swizzle
           channel ordering.

       -dsw <swizzle>
           Swizzle the color components after decompression. The swizzle is
           specified using the same method as the -esw option, with support
           for an additional "z" character. This is used to specify that
           the compressed data stores an X+Y normal map, and that the Z
           output channel should be reconstructed from the two channels
           stored in the data. For the typical ASTC normal encoding,
           which uses an 'rrrg' compression swizzle, you should specify an
           'raz1' swizzle for decompression.

       -yflip
           Flip the image in the vertical axis prior to compression and
           after decompression. Note that using this option in a test mode
           (-t*) will have no effect as the image will be flipped twice.

       -j <threads>
           Explicitly specify the number of compression/decompression
           theads to use in the codec. If not specified, the codec will
           use on thread per CPU detected in the system.

       -silent
           Suppresses all non-essential diagnostic output from the codec.
           Error messages will always be printed, as will mandatory outputs
           for the selected operation mode. For example, the test mode
           will always output image quality metrics and compression time
           but will suppress all other output.

DECOMPRESSION
       To decompress an image stored in the ASTC format you must specify
       the color profile, the input file name, and the output file name.

       The color profile is specified using the -dl (LDR linear), -ds
       (LDR sRGB), -dh (HDR RGB, LDR A), or -dH (HDR RGBA) decoder options.

       The input file path must match a valid file format for
       decompression, and the output file format must be a valid output for
       a decompressed image. Note that not all output formats that the
       coompression path can produce are supported for decompression. See
       the FILE FORMATS section for the list of supported formats.

       The -dsw options documented in ADVANCED COMPRESSION option documentation
	   are relevent to decompression.

TEST
       To perform a compression test which round-trips a single image
       through compression and decompression and stores the decompressed
       result back to file, you must specify same settings as COMPRESSION
       other than swapping the color profile to select test mode. Note that
       the compressed intermediate data is discarded in this mode.

       The color profile is specified using the -tl (LDR linear), -ts (LDR
       sRGB), -th (HDR RGB, LDR A), or -tH (HDR RGBA) encoder options.

       This operation mode will print error metrics suitable for either
       LDR and HDR images, allowing some assessment of the compression
       image quality.)"
// This split in the literals is needed for Visual Studio; the compiler
// will concatenate these two strings together ...
R"(

COMPRESSION FILE FORMATS
       The following formats are supported as compression inputs:

           LDR Formats:
               BMP (*.bmp)
               PNG (*.png)
               Targa (*.tga)
               JPEG (*.jpg)

           HDR Formats:
               OpenEXR (*.exr)
               Radiance HDR (*.hdr)

           Container Formats:
               Khronos Texture KTX (*.ktx)
               DirectDraw Surface DDS (*.dds)

       For the KTX and DDS formats only a subset of the features of the
       formats are supported:

           Texture topology must be 2D, 2D-array, 3D, or cube-map. Note
           that 2D-array textures are treated as 3D block input.

           Texel format must be R, RG, RGB, BGR, RGBA, BGRA, L, or LA.

           Only the first mipmap in the file will be read.

       The following formats are supported as compression outputs:

           ASTC (*.astc)
           Khronos Texture KTX (*.ktx)


DECOMPRESSION FILE FORMATS
       The following formats are supported as decompression inputs:

           ASTC (*.astc)
           Khronos Texture KTX (*.ktx)

       The following formats are supported as decompression outputs:

           LDR Formats:
               BMP (*.bmp)
               PNG (*.png)
               Targa (*.tga)

           HDR Formats:
               OpenEXR (*.exr)
               Radiance HDR (*.hdr)

           Container Formats:
               Khronos Texture KTX (*.ktx)
               DirectDraw Surface DDS (*.dds)
)";

// print version and basic build information
void astcenc_print_header()
{
#if (ASTCENC_AVX == 2)
	const char* simdtype = "avx2";
#elif (ASTCENC_SSE == 42)
	const char* simdtype = "sse4.2";
#else
	const char* simdtype = "sse2";
#endif

#if (ASTCENC_POPCNT == 1)
	const char* pcnttype = "+popcnt";
#else
	const char* pcnttype = "";
#endif

	unsigned int bits = (int)(sizeof(void*) * 8);
	printf(astcenc_copyright_string, bits, simdtype, pcnttype);
}

void astcenc_print_shorthelp() {
	astcenc_print_header();
	printf("%s", astcenc_short_help);
}

void astcenc_print_longhelp() {
	astcenc_print_header();
	printf("%s", astcenc_long_help);
}
