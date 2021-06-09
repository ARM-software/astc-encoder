#include <iostream>

#define ASTCENC_SSE 20
#include "astcenc_mathlib.h"

/**
 * @brief Convert unorm16 [0, 65535] to float16 in range [0, 1].
 */
static ASTCENC_SIMD_INLINE vint4 unorm16_to_sf16_bugged(vint4 p)
{
	vint4 fp16_one = vint4(0x3C00);
	vint4 fp16_small = lsl<8>(p);

	vmask4 is_one = p == vint4(0xFFFF);
	vmask4 is_small = p < vint4(4);

	// Manually inline clz() on Visual Studio to avoid release build codegen bug
#if 0 && !defined(__clang__) && defined(_MSC_VER)
	vint4 a = (~lsr<8>(p)) & p;
	a = float_as_int(int_to_float(a));
	a = vint4(127 + 31) - lsr<23>(a);
	vint4 lz = clamp(0, 32, a) - 16;
#else
	vint4 lz = clz(p) - 16;
#endif
	// The value of p is corrupted after calling clz()

	p = p * two_to_the_n(lz + 1);
	p = p & vint4(0xFFFF);

	p = lsr<6>(p);

	p = p | lsl<10>(vint4(14) - lz);

	vint4 r = select(p, fp16_one, is_one);
	r = select(r, fp16_small, is_small);
	return r;
}


int main()
{
	vint4 value(65519);

	// This function inlines vint4 clz() as a workaround for the issue, which 
	// masks the problem and gives the correct result.
	vint4 result_good = unorm16_to_sf16(value);

	// This function uses the original code, calling clz() as a function, 
	// which corrupts the value of p in the caller in Release builds.
	vint4 result_bad = unorm16_to_sf16_bugged(value);

	print(result_good);
	print(result_bad);

	if (any(result_good != result_bad))
	{
		puts("Failed ...\n");
		return 1;
	}
	
	puts("Success ...\n");
	return 0;
}
