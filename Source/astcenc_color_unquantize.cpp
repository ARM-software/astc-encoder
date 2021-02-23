// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2011-2021 Arm Limited
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

#include <utility>

/**
 * @brief Functions for color unquantization.
 */

#include "astcenc_internal.h"

static ASTCENC_SIMD_INLINE vint4 unquant_color(
	int quant_level,
	vint4 inputq
) {
	const uint8_t* unq = color_unquant_tables[quant_level];
	return vint4(unq[inputq.lane<0>()], unq[inputq.lane<1>()],
	             unq[inputq.lane<2>()], unq[inputq.lane<3>()]);
}

static ASTCENC_SIMD_INLINE vint4 uncontract_color(
	vint4 input
) {
	vmask4 mask(true, true, false, false);
	vint4 bc0 = asr<1>(input + input.lane<2>());
	return select(input, bc0, mask);
}

static void rgba_delta_unpack(
	vint4 input0q,
	vint4 input1q,
	int quant_level,
	vint4& output0,
	vint4& output1
) {
	// Unquantize color endpoints
	vint4 input0 = unquant_color(quant_level, input0q);
	vint4 input1 = unquant_color(quant_level, input1q);

	// Perform bit-transfer
	input0 = input0 | lsl<1>(input1 & 0x80);
	input1 = input1 & 0x7F;
	vmask4 mask = (input1 & 0x40) != vint4::zero();
	input1 = select(input1, input1 - 0x80, mask);

	// Scale
	input0 = asr<1>(input0);
	input1 = asr<1>(input1);

	// Apply blue-uncontraction if needed
	int rgb_sum = hadd_rgb_s(input1);
	input1 = input1 + input0;
	if (rgb_sum < 0)
	{
		input0 = uncontract_color(input0);
		input1 = uncontract_color(input1);
		std::swap(input0, input1);
	}

	output0 = clamp(0, 255, input0);
	output1 = clamp(0, 255, input1);
}

static void rgb_delta_unpack(
	vint4 input0q,
	vint4 input1q,
	int quant_level,
	vint4& output0,
	vint4& output1
) {
	rgba_delta_unpack(input0q, input1q, quant_level, output0, output1);
	output0.set_lane<3>(255);
	output1.set_lane<3>(255);
}

static void rgba_unpack(
	vint4 input0q,
	vint4 input1q,
	int quant_level,
	vint4& output0,
	vint4& output1
) {
	// Unquantize color endpoints
	vint4 input0 = unquant_color(quant_level, input0q);
	vint4 input1 = unquant_color(quant_level, input1q);

	// Apply blue-uncontraction if needed
	if (hadd_rgb_s(input0) > hadd_rgb_s(input1))
	{
		input0 = uncontract_color(input0);
		input1 = uncontract_color(input1);
		std::swap(input0, input1);
	}

	output0 = input0;
	output1 = input1;
}

static void rgb_unpack(
	vint4 input0q,
	vint4 input1q,
	int quant_level,
	vint4& output0,
	vint4& output1
) {
	rgba_unpack(input0q, input1q, quant_level, output0, output1);
	output0.set_lane<3>(255);
	output1.set_lane<3>(255);
}

static void rgb_scale_alpha_unpack(
	vint4 input0q,
	int alpha1q,
	int scaleq,
	int quant_level,
	vint4& output0,
	vint4& output1
) {
	// Unquantize color endpoints
	vint4 input = unquant_color(quant_level, input0q);
	int alpha1 = color_unquant_tables[quant_level][alpha1q];
	int scale = color_unquant_tables[quant_level][scaleq];

	output1 = input;
	output1.set_lane<3>(alpha1);

	output0 = asr<8>(input * scale);
	output0.set_lane<3>(input.lane<3>());
}

static void rgb_scale_unpack(
	vint4 input0q,
	int scaleq,
	int quant_level,
	vint4& output0,
	vint4& output1
) {
	vint4 input = unquant_color(quant_level, input0q);
	int scale = color_unquant_tables[quant_level][scaleq];

	output1 = input;
	output1.set_lane<3>(255);

	output0 = asr<8>(input * scale);
	output0.set_lane<3>(255);
}

static void luminance_unpack(
	const int input[2],
	int quant_level,
	vint4* output0,
	vint4* output1
) {
	int lum0 = color_unquant_tables[quant_level][input[0]];
	int lum1 = color_unquant_tables[quant_level][input[1]];
	*output0 = vint4(lum0, lum0, lum0, 255);
	*output1 = vint4(lum1, lum1, lum1, 255);
}

static void luminance_delta_unpack(
	const int input[2],
	int quant_level,
	vint4* output0,
	vint4* output1
) {
	int v0 = color_unquant_tables[quant_level][input[0]];
	int v1 = color_unquant_tables[quant_level][input[1]];
	int l0 = (v0 >> 2) | (v1 & 0xC0);
	int l1 = l0 + (v1 & 0x3F);

	l1 = astc::min(l1, 255);

	*output0 = vint4(l0, l0, l0, 255);
	*output1 = vint4(l1, l1, l1, 255);
}

static void luminance_alpha_unpack(
	const int input[4],
	int quant_level,
	vint4* output0,
	vint4* output1
) {
	int lum0 = color_unquant_tables[quant_level][input[0]];
	int lum1 = color_unquant_tables[quant_level][input[1]];
	int alpha0 = color_unquant_tables[quant_level][input[2]];
	int alpha1 = color_unquant_tables[quant_level][input[3]];
	*output0 = vint4(lum0, lum0, lum0, alpha0);
	*output1 = vint4(lum1, lum1, lum1, alpha1);
}

static void luminance_alpha_delta_unpack(
	const int input[4],
	int quant_level,
	vint4* output0,
	vint4* output1
) {
	int lum0 = color_unquant_tables[quant_level][input[0]];
	int lum1 = color_unquant_tables[quant_level][input[1]];
	int alpha0 = color_unquant_tables[quant_level][input[2]];
	int alpha1 = color_unquant_tables[quant_level][input[3]];

	lum0 |= (lum1 & 0x80) << 1;
	alpha0 |= (alpha1 & 0x80) << 1;
	lum1 &= 0x7F;
	alpha1 &= 0x7F;
	if (lum1 & 0x40)
		lum1 -= 0x80;
	if (alpha1 & 0x40)
		alpha1 -= 0x80;

	lum0 >>= 1;
	lum1 >>= 1;
	alpha0 >>= 1;
	alpha1 >>= 1;
	lum1 += lum0;
	alpha1 += alpha0;

	lum1 = astc::clamp(lum1, 0, 255);
	alpha1 = astc::clamp(alpha1, 0, 255);

	*output0 = vint4(lum0, lum0, lum0, alpha0);
	*output1 = vint4(lum1, lum1, lum1, alpha1);
}

// RGB-offset format
static void hdr_rgbo_unpack3(
	const int input[4],
	int quant_level,
	vint4* output0,
	vint4* output1
) {
	int v0 = color_unquant_tables[quant_level][input[0]];
	int v1 = color_unquant_tables[quant_level][input[1]];
	int v2 = color_unquant_tables[quant_level][input[2]];
	int v3 = color_unquant_tables[quant_level][input[3]];

	int modeval = ((v0 & 0xC0) >> 6) | (((v1 & 0x80) >> 7) << 2) | (((v2 & 0x80) >> 7) << 3);

	int majcomp;
	int mode;
	if ((modeval & 0xC) != 0xC)
	{
		majcomp = modeval >> 2;
		mode = modeval & 3;
	}
	else if (modeval != 0xF)
	{
		majcomp = modeval & 3;
		mode = 4;
	}
	else
	{
		majcomp = 0;
		mode = 5;
	}

	int red = v0 & 0x3F;
	int green = v1 & 0x1F;
	int blue = v2 & 0x1F;
	int scale = v3 & 0x1F;

	int bit0 = (v1 >> 6) & 1;
	int bit1 = (v1 >> 5) & 1;
	int bit2 = (v2 >> 6) & 1;
	int bit3 = (v2 >> 5) & 1;
	int bit4 = (v3 >> 7) & 1;
	int bit5 = (v3 >> 6) & 1;
	int bit6 = (v3 >> 5) & 1;

	int ohcomp = 1 << mode;

	if (ohcomp & 0x30)
		green |= bit0 << 6;
	if (ohcomp & 0x3A)
		green |= bit1 << 5;
	if (ohcomp & 0x30)
		blue |= bit2 << 6;
	if (ohcomp & 0x3A)
		blue |= bit3 << 5;

	if (ohcomp & 0x3D)
		scale |= bit6 << 5;
	if (ohcomp & 0x2D)
		scale |= bit5 << 6;
	if (ohcomp & 0x04)
		scale |= bit4 << 7;

	if (ohcomp & 0x3B)
		red |= bit4 << 6;
	if (ohcomp & 0x04)
		red |= bit3 << 6;

	if (ohcomp & 0x10)
		red |= bit5 << 7;
	if (ohcomp & 0x0F)
		red |= bit2 << 7;

	if (ohcomp & 0x05)
		red |= bit1 << 8;
	if (ohcomp & 0x0A)
		red |= bit0 << 8;

	if (ohcomp & 0x05)
		red |= bit0 << 9;
	if (ohcomp & 0x02)
		red |= bit6 << 9;

	if (ohcomp & 0x01)
		red |= bit3 << 10;
	if (ohcomp & 0x02)
		red |= bit5 << 10;

	// expand to 12 bits.
	static const int shamts[6] { 1, 1, 2, 3, 4, 5 };
	int shamt = shamts[mode];
	red <<= shamt;
	green <<= shamt;
	blue <<= shamt;
	scale <<= shamt;

	// on modes 0 to 4, the values stored for "green" and "blue" are differentials,
	// not absolute values.
	if (mode != 5)
	{
		green = red - green;
		blue = red - blue;
	}

	// switch around components.
	int temp;
	switch (majcomp)
	{
	case 1:
		temp = red;
		red = green;
		green = temp;
		break;
	case 2:
		temp = red;
		red = blue;
		blue = temp;
		break;
	default:
		break;
	}

	int red0 = red - scale;
	int green0 = green - scale;
	int blue0 = blue - scale;

	// clamp to [0,0xFFF].
	if (red < 0)
		red = 0;
	if (green < 0)
		green = 0;
	if (blue < 0)
		blue = 0;

	if (red0 < 0)
		red0 = 0;
	if (green0 < 0)
		green0 = 0;
	if (blue0 < 0)
		blue0 = 0;

	*output0 = vint4(red0 << 4, green0 << 4, blue0 << 4, 0x7800);
	*output1 = vint4(red << 4, green << 4, blue << 4, 0x7800);
}

static void hdr_rgb_unpack3(
	const int input[6],
	int quant_level,
	vint4* output0,
	vint4* output1
) {

	int v0 = color_unquant_tables[quant_level][input[0]];
	int v1 = color_unquant_tables[quant_level][input[1]];
	int v2 = color_unquant_tables[quant_level][input[2]];
	int v3 = color_unquant_tables[quant_level][input[3]];
	int v4 = color_unquant_tables[quant_level][input[4]];
	int v5 = color_unquant_tables[quant_level][input[5]];

	// extract all the fixed-placement bitfields
	int modeval = ((v1 & 0x80) >> 7) | (((v2 & 0x80) >> 7) << 1) | (((v3 & 0x80) >> 7) << 2);

	int majcomp = ((v4 & 0x80) >> 7) | (((v5 & 0x80) >> 7) << 1);

	if (majcomp == 3)
	{
		*output0 = vint4(v0 << 8, v2 << 8, (v4 & 0x7F) << 9, 0x7800);
		*output1 = vint4(v1 << 8, v3 << 8, (v5 & 0x7F) << 9, 0x7800);
		return;
	}

	int a = v0 | ((v1 & 0x40) << 2);
	int b0 = v2 & 0x3f;
	int b1 = v3 & 0x3f;
	int c = v1 & 0x3f;
	int d0 = v4 & 0x7f;
	int d1 = v5 & 0x7f;

	// get hold of the number of bits in 'd0' and 'd1'
	static const int dbits_tab[8] { 7, 6, 7, 6, 5, 6, 5, 6 };
	int dbits = dbits_tab[modeval];

	// extract six variable-placement bits
	int bit0 = (v2 >> 6) & 1;
	int bit1 = (v3 >> 6) & 1;
	int bit2 = (v4 >> 6) & 1;
	int bit3 = (v5 >> 6) & 1;
	int bit4 = (v4 >> 5) & 1;
	int bit5 = (v5 >> 5) & 1;

	// and prepend the variable-placement bits depending on mode.
	int ohmod = 1 << modeval;	// one-hot-mode
	if (ohmod & 0xA4)
		a |= bit0 << 9;
	if (ohmod & 0x8)
		a |= bit2 << 9;
	if (ohmod & 0x50)
		a |= bit4 << 9;

	if (ohmod & 0x50)
		a |= bit5 << 10;
	if (ohmod & 0xA0)
		a |= bit1 << 10;

	if (ohmod & 0xC0)
		a |= bit2 << 11;

	if (ohmod & 0x4)
		c |= bit1 << 6;
	if (ohmod & 0xE8)
		c |= bit3 << 6;

	if (ohmod & 0x20)
		c |= bit2 << 7;

	if (ohmod & 0x5B)
	{
		b0 |= bit0 << 6;
		b1 |= bit1 << 6;
	}

	if (ohmod & 0x12)
	{
		b0 |= bit2 << 7;
		b1 |= bit3 << 7;
	}

	if (ohmod & 0xAF)
	{
		d0 |= bit4 << 5;
		d1 |= bit5 << 5;
	}

	if (ohmod & 0x5)
	{
		d0 |= bit2 << 6;
		d1 |= bit3 << 6;
	}

	// sign-extend 'd0' and 'd1'
	// note: this code assumes that signed right-shift actually sign-fills, not zero-fills.
	int32_t d0x = d0;
	int32_t d1x = d1;
	int sx_shamt = 32 - dbits;
	d0x <<= sx_shamt;
	d0x >>= sx_shamt;
	d1x <<= sx_shamt;
	d1x >>= sx_shamt;
	d0 = d0x;
	d1 = d1x;

	// expand all values to 12 bits, with left-shift as needed.
	int val_shamt = (modeval >> 1) ^ 3;
	a <<= val_shamt;
	b0 <<= val_shamt;
	b1 <<= val_shamt;
	c <<= val_shamt;
	d0 <<= val_shamt;
	d1 <<= val_shamt;

	// then compute the actual color values.
	int red1 = a;
	int green1 = a - b0;
	int blue1 = a - b1;
	int red0 = a - c;
	int green0 = a - b0 - c - d0;
	int blue0 = a - b1 - c - d1;

	// clamp the color components to [0,2^12 - 1]
	red0 = astc::clamp(red0, 0, 4095);
	green0 = astc::clamp(green0, 0, 4095);
	blue0 = astc::clamp(blue0, 0, 4095);

	red1 = astc::clamp(red1, 0, 4095);
	green1 = astc::clamp(green1, 0, 4095);
	blue1 = astc::clamp(blue1, 0, 4095);

	// switch around the color components
	int temp0, temp1;
	switch (majcomp)
	{
	case 1:					// switch around red and green
		temp0 = red0;
		temp1 = red1;
		red0 = green0;
		red1 = green1;
		green0 = temp0;
		green1 = temp1;
		break;
	case 2:					// switch around red and blue
		temp0 = red0;
		temp1 = red1;
		red0 = blue0;
		red1 = blue1;
		blue0 = temp0;
		blue1 = temp1;
		break;
	case 0:					// no switch
		break;
	}

	*output0 = vint4(red0 << 4, green0 << 4, blue0 << 4, 0x7800);
	*output1 = vint4(red1 << 4, green1 << 4, blue1 << 4, 0x7800);
}

static void hdr_rgb_ldr_alpha_unpack3(
	const int input[8],
	int quant_level,
	vint4* output0,
	vint4* output1
) {
	hdr_rgb_unpack3(input, quant_level, output0, output1);

	int v6 = color_unquant_tables[quant_level][input[6]];
	int v7 = color_unquant_tables[quant_level][input[7]];
	output0->set_lane<3>(v6);
	output1->set_lane<3>(v7);
}

static void hdr_luminance_small_range_unpack(
	const int input[2],
	int quant_level,
	vint4* output0,
	vint4* output1
) {
	int v0 = color_unquant_tables[quant_level][input[0]];
	int v1 = color_unquant_tables[quant_level][input[1]];

	int y0, y1;
	if (v0 & 0x80)
	{
		y0 = ((v1 & 0xE0) << 4) | ((v0 & 0x7F) << 2);
		y1 = (v1 & 0x1F) << 2;
	}
	else
	{
		y0 = ((v1 & 0xF0) << 4) | ((v0 & 0x7F) << 1);
		y1 = (v1 & 0xF) << 1;
	}

	y1 += y0;
	if (y1 > 0xFFF)
		y1 = 0xFFF;

	*output0 = vint4(y0 << 4, y0 << 4, y0 << 4, 0x7800);
	*output1 = vint4(y1 << 4, y1 << 4, y1 << 4, 0x7800);
}

static void hdr_luminance_large_range_unpack(
	const int input[2],
	int quant_level,
	vint4* output0,
	vint4* output1
) {
	int v0 = color_unquant_tables[quant_level][input[0]];
	int v1 = color_unquant_tables[quant_level][input[1]];

	int y0, y1;
	if (v1 >= v0)
	{
		y0 = v0 << 4;
		y1 = v1 << 4;
	}
	else
	{
		y0 = (v1 << 4) + 8;
		y1 = (v0 << 4) - 8;
	}
	*output0 = vint4(y0 << 4, y0 << 4, y0 << 4, 0x7800);
	*output1 = vint4(y1 << 4, y1 << 4, y1 << 4, 0x7800);
}

static void hdr_alpha_unpack(
	const int input[2],
	int quant_level,
	int* output0,
	int* output1
) {

	int v6 = color_unquant_tables[quant_level][input[0]];
	int v7 = color_unquant_tables[quant_level][input[1]];

	int selector = ((v6 >> 7) & 1) | ((v7 >> 6) & 2);
	v6 &= 0x7F;
	v7 &= 0x7F;
	if (selector == 3)
	{
		*output0 = v6 << 5;
		*output1 = v7 << 5;
	}
	else
	{
		v6 |= (v7 << (selector + 1)) & 0x780;
		v7 &= (0x3f >> selector);
		v7 ^= 32 >> selector;
		v7 -= 32 >> selector;
		v6 <<= (4 - selector);
		v7 <<= (4 - selector);
		v7 += v6;

		if (v7 < 0)
			v7 = 0;
		else if (v7 > 0xFFF)
			v7 = 0xFFF;

		*output0 = v6;
		*output1 = v7;
	}

	*output0 <<= 4;
	*output1 <<= 4;
}

static void hdr_rgb_hdr_alpha_unpack3(
	const int input[8],
	int quant_level,
	vint4* output0,
	vint4* output1
) {
	hdr_rgb_unpack3(input, quant_level, output0, output1);

	int alpha0, alpha1;
	hdr_alpha_unpack(input + 6, quant_level, &alpha0, &alpha1);

	output0->set_lane<3>(alpha0);
	output1->set_lane<3>(alpha1);
}

void unpack_color_endpoints(
	astcenc_profile decode_mode,
	int format,
	int quant_level,
	const int* input,
	int* rgb_hdr,
	int* alpha_hdr,
	int* nan_endpoint,
	vint4* output0,
	vint4* output1
) {
	// TODO: Make these bools ...

	// Assume no NaNs and LDR endpoints

	// TODO: Review use of NaN endpoint. It's never set for HDR images ...
	*nan_endpoint = 0;
	*rgb_hdr = 0;
	*alpha_hdr = 0;


	switch (format)
	{
	case FMT_LUMINANCE:
		luminance_unpack(input, quant_level, output0, output1);
		break;

	case FMT_LUMINANCE_DELTA:
		luminance_delta_unpack(input, quant_level, output0, output1);
		break;

	case FMT_HDR_LUMINANCE_SMALL_RANGE:
		*rgb_hdr = 1;
		*alpha_hdr = -1;
		hdr_luminance_small_range_unpack(input, quant_level, output0, output1);
		break;

	case FMT_HDR_LUMINANCE_LARGE_RANGE:
		*rgb_hdr = 1;
		*alpha_hdr = -1;
		hdr_luminance_large_range_unpack(input, quant_level, output0, output1);
		break;

	case FMT_LUMINANCE_ALPHA:
		luminance_alpha_unpack(input, quant_level, output0, output1);
		break;

	case FMT_LUMINANCE_ALPHA_DELTA:
		luminance_alpha_delta_unpack(input, quant_level, output0, output1);
		break;

	case FMT_RGB_SCALE:
		{
			vint4 input0q(input[0], input[1], input[2], 0);
			int scale = input[3];
			rgb_scale_unpack(input0q, scale, quant_level, *output0, *output1);
		}
		break;

	case FMT_RGB_SCALE_ALPHA:
		{
			vint4 input0q(input[0], input[1], input[2], input[4]);
			int alpha1q = input[5];
			int scaleq = input[3];
			rgb_scale_alpha_unpack(input0q, alpha1q, scaleq, quant_level, *output0, *output1);
		}
		break;

	case FMT_HDR_RGB_SCALE:
		*rgb_hdr = 1;
		*alpha_hdr = -1;
		hdr_rgbo_unpack3(input, quant_level, output0, output1);
		break;

	case FMT_RGB:
		{
			vint4 input0q(input[0], input[2], input[4], 0);
			vint4 input1q(input[1], input[3], input[5], 0);
			rgb_unpack(input0q, input1q, quant_level, *output0, *output1);
		}
		break;

	case FMT_RGB_DELTA:
		{
			vint4 input0q(input[0], input[2], input[4], 0);
			vint4 input1q(input[1], input[3], input[5], 0);
			rgb_delta_unpack(input0q, input1q, quant_level, *output0, *output1);
		}
		break;

	case FMT_HDR_RGB:
		*rgb_hdr = 1;
		*alpha_hdr = -1;
		hdr_rgb_unpack3(input, quant_level, output0, output1);
		break;

	case FMT_RGBA:
		{
			vint4 input0q(input[0], input[2], input[4], input[6]);
			vint4 input1q(input[1], input[3], input[5], input[7]);
			rgba_unpack(input0q, input1q, quant_level, *output0, *output1);
		}
		break;

	case FMT_RGBA_DELTA:
		{
			vint4 input0q(input[0], input[2], input[4], input[6]);
			vint4 input1q(input[1], input[3], input[5], input[7]);
			rgba_delta_unpack(input0q, input1q, quant_level, *output0, *output1);
		}
		break;

	case FMT_HDR_RGB_LDR_ALPHA:
		*rgb_hdr = 1;
		hdr_rgb_ldr_alpha_unpack3(input, quant_level, output0, output1);
		break;

	case FMT_HDR_RGBA:
		*rgb_hdr = 1;
		*alpha_hdr = 1;
		hdr_rgb_hdr_alpha_unpack3(input, quant_level, output0, output1);
		break;
	}

	// Assign a correct default alpha
	if (*alpha_hdr == -1)
	{
		if (decode_mode == ASTCENC_PRF_HDR)
		{
			output0->set_lane<3>(0x7800);
			output1->set_lane<3>(0x7800);
			*alpha_hdr = 1;
		}
		else
		{
			output0->set_lane<3>(0x00FF);
			output1->set_lane<3>(0x00FF);
			*alpha_hdr = 0;
		}
	}

	vint4 ldr_scale(257);
	vint4 hdr_scale(1);
	vint4 output_scale = ldr_scale;

	// An LDR profile image
	if ((decode_mode == ASTCENC_PRF_LDR_SRGB) ||
	    (decode_mode == ASTCENC_PRF_LDR_SRGB))
	{
		// Also matches HDR alpha, as cannot have HDR alpha without HDR RGB
		if (*rgb_hdr == 1)
		{
			*output0 = vint4(0xFF00, 0x0000, 0xFF00, 0xFF00);
			*output1 = vint4(0xFF00, 0x0000, 0xFF00, 0xFF00);
			output_scale = hdr_scale;

			*rgb_hdr = 0;
			*alpha_hdr = 0;
		}
	}
	// An HDR profile image
	else
	{
		bool hrgb = *rgb_hdr == 1;
		bool ha = *alpha_hdr == 1;
		vmask4 hdr_lanes(hrgb, hrgb, hrgb, ha);
		output_scale = select(ldr_scale, hdr_scale, hdr_lanes);
	}

	*output0 = *output0 * output_scale;
	*output1 = *output1 * output_scale;
}
