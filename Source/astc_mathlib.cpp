// ----------------------------------------------------------------------------
//  This confidential and proprietary software may be used only as authorised
//  by a licensing agreement from Arm Limited.
//      (C) COPYRIGHT 2011-2020 Arm Limited, ALL RIGHTS RESERVED
//  The entire notice above must be reproduced on all authorised copies and
//  copies may only be made to the extent permitted by a licensing agreement
//  from Arm Limited.
// ----------------------------------------------------------------------------

#include "astc_mathlib.h"

float3 cross(float3 p, float3 q)
{
	return float3(p.y * q.z - p.z * q.y,
	              p.z * q.x - p.x * q.z,
	              p.x * q.y - p.y * q.x);
}

float determinant(mat2 p)
{
	float2 v = float2(p.v[0].x * p.v[1].y, p.v[0].y * p.v[1].x);
	return v.x - v.y;
}

float2 transform(mat2 p, float2 q)
{
	return float2(dot(p.v[0], q), dot(p.v[1], q));
}

float4 transform(mat4 p, float4 q)
{
	return float4(dot(p.v[0], q), dot(p.v[1], q), dot(p.v[2], q), dot(p.v[3], q));
}

mat2 invert(mat2 p)
{
	float rdet = 1.0f / determinant(p);
	mat2 res;
	res.v[0] = float2(p.v[1].y, -p.v[0].y) * rdet;
	res.v[1] = float2(-p.v[1].x, p.v[0].x) * rdet;
	return res;
}

mat4 invert(mat4 p)
{
	// cross products between the bottom two rows
	float3 bpc0 = cross(float3(p.v[2].y, p.v[2].z, p.v[2].w), float3(p.v[3].y, p.v[3].z, p.v[3].w));
	float3 bpc1 = cross(float3(p.v[2].x, p.v[2].z, p.v[2].w), float3(p.v[3].x, p.v[3].z, p.v[3].w));
	float3 bpc2 = cross(float3(p.v[2].x, p.v[2].y, p.v[2].w), float3(p.v[3].x, p.v[3].y, p.v[3].w));
	float3 bpc3 = cross(float3(p.v[2].x, p.v[2].y, p.v[2].z), float3(p.v[3].x, p.v[3].y, p.v[3].z));

	// dot-products for the top rows
	float4 row1 = float4(dot(bpc0, float3(p.v[1].y, p.v[1].z, p.v[1].w)),
	                    -dot(bpc1, float3(p.v[1].x, p.v[1].z, p.v[1].w)),
	                     dot(bpc2, float3(p.v[1].x, p.v[1].y, p.v[1].w)),
	                    -dot(bpc3, float3(p.v[1].x, p.v[1].y, p.v[1].z)));

	float det = dot(p.v[0], row1);
	float rdet = 1.0f / det;

	mat4 res;

	float3 tpc0 = cross(float3(p.v[0].y, p.v[0].z, p.v[0].w), float3(p.v[1].y, p.v[1].z, p.v[1].w));
	res.v[0] = float4(row1.x,
	                 -dot(bpc0, float3(p.v[0].y, p.v[0].z, p.v[0].w)),
	                  dot(tpc0, float3(p.v[3].y, p.v[3].z, p.v[3].w)),
	                 -dot(tpc0, float3(p.v[2].y, p.v[2].z, p.v[2].w))) * rdet;

	float3 tpc1 = cross(float3(p.v[0].x, p.v[0].z, p.v[0].w), float3(p.v[1].x, p.v[1].z, p.v[1].w));
	res.v[1] = float4(row1.y,
	                  dot(bpc1, float3(p.v[0].x, p.v[0].z, p.v[0].w)),
	                 -dot(tpc1, float3(p.v[3].x, p.v[3].z, p.v[3].w)),
	                  dot(tpc1, float3(p.v[2].x, p.v[2].z, p.v[2].w))) * rdet;

	float3 tpc2 = cross(float3(p.v[0].x, p.v[0].y, p.v[0].w), float3(p.v[1].x, p.v[1].y, p.v[1].w));
	res.v[2] = float4(row1.z,
	                 -dot(bpc2, float3(p.v[0].x, p.v[0].y, p.v[0].w)),
	                  dot(tpc2, float3(p.v[3].x, p.v[3].y, p.v[3].w)),
	                 -dot(tpc2, float3(p.v[2].x, p.v[2].y, p.v[2].w))) * rdet;

	float3 tpc3 = cross(float3(p.v[0].x, p.v[0].z, p.v[0].z), float3(p.v[1].x, p.v[1].y, p.v[1].z));
	res.v[3] = float4(row1.w,
	                  dot(bpc3, float3(p.v[0].x, p.v[0].y, p.v[0].z)),
	                 -dot(tpc3, float3(p.v[3].x, p.v[3].y, p.v[3].z)),
	                  dot(tpc3, float3(p.v[2].x, p.v[2].y, p.v[2].z))) * rdet;

	return res;
}

/* Public function, see header file for detailed documentation */
float astc::log2(float val)
{
	if32 p;
	p.f = val;
	if( p.s < 0x800000 )
		p.s = 0x800000; // negative, 0, denormal get clamped to non-denormal.

	// normalize mantissa to range [0.66, 1.33] and extract an exponent
	// in such a way that 1.0 returns 0.
	p.s -= 0x3f2aaaab;
	int expo = p.s >> 23;
	p.s &= 0x7fffff;
	p.s += 0x3f2aaaab;

	float x = p.f - 1.0f;

	// taylor polynomial that, with horner's-rule style evaluation,
	// gives sufficient precision for our use
	// (relative error of about 1 in 10^6)

	float res = (float)expo
	          + x * ( 1.442695040888963f
	          + x * (-0.721347520444482f
	          + x * ( 0.480898346962988f
	          + x * (-0.360673760222241f
	          + x * ( 0.288539008177793f
	          + x * (-0.240449173481494f
	          + x * ( 0.206099291555566f
	          + x * (-0.180336880111120f
	          + x * ( 0.160299448987663f
	          )))))))));
	return res;
}

/* Public function, see header file for detailed documentation */
float astc::atan2(float y, float x)
{
	const float PI = (float)M_PI;
	const float PI_2 = PI / 2.f;

	// Handle the discontinuity at x == 0
	if (x == 0.0f)
	{
		if (y > 0.0f)
		{
			return PI_2;
		}
		else if (y == 0.0f)
		{
			return 0.0f;
		}
		return -PI_2;
	}

	float z = y / x;
	float z2 = z * z;
	if (std::fabs(z) < 1.0f)
	{
		float atan = z / (1.0f + (0.28f * z2));
		if (x < 0.0f)
		{
			if (y < 0.0f)
			{
				return atan - PI;
			}
			else
			{
				return atan + PI;
			}
		}
		return atan;
	}
	else
	{
		float atan = PI_2 - (z / (z2 + 0.28f));
		if (y < 0.0f)
		{
			return atan - PI;
		}
		else
		{
			return atan;
		}
	}
}
