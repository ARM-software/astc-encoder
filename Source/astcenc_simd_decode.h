// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2022 Arm Limited
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

// ---------------------------------------------------------------------------
// ASTC SIMD Decode library.
//
// This library is provided as a single-file header in a style similar to
// "stb_image.h". In a program that uses the library, it can be included
// from multiple source code files, however, exactly ONE of them must
// define the macro ASTC_SIMD_DECODE_IMPLEMENTATION before including it.
// This inclusion will provide the library; all other inclusions will
// include only the public interface.
//
// The implementation must be compiled with support for one of the
// following SIMD ISAs:
//   * Intel SSSE3
//   * Intel AVX2
//   * ARMv8 NEON
// On Intel, it will also use BMI2 if available, for the PEXT and PDEP
// instructions.
//
// This library provides partial support for multithreading. This
// support takes the form of allowing different invocations of the
// main workhorse function 'astc_simd_decode_row_iterate()' function
// to be run in parallel. (The actual work of setting up a thread-pool
// that can launch these invocations in parallel is left to calling code.)
// ---------------------------------------------------------------------------



// Initial run: for astc 8x8 on a 4GHz haswell, LDR decode speed
// is about 730 MTex/s or about 6.0 cycles per texel,
// corresponding to about 385 cycles per block.
//
// The final loop makes up only about 1.0 cycles per texel.
// Index infill makes up about 1.37 cycles per texel.
// ISE-unpacking takes about 1.17 cycles per texel.
// Color endpoint unpacking takes about 0.41 cycles per texel.
//
// For other block sizes, measured single-thread performance is
// (in MTexels/s on 4GHz haswell, using gcc -O3 on a 3072x3072
// composite of the 24 kodim images):
//            LDR   HDR
//  *  4x4:    277   177
//  *  5x5:    406   227
//  *  6x6:    523   282
//  *  8x8:    730   379
//  * 10x10:   775   358
//  * 12x12:  1030   412




#if !defined(ASTC_SIMD_DECODE_H_INCLUDED) && !defined(ASTCSD_RECURSIVE_INCLUDE)
#define ASTC_SIMD_DECODE_H_INCLUDED

// *****************************************************************************
// *****************************************************************************
// *****                                                                   *****
// *****       Public interface of the ASTC SIMD Decode Library            *****
// *****                                                                   *****
// *****************************************************************************
// *****************************************************************************

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

// -------------------
//   Data structures
// -------------------
struct astc_simd_decode_processed_params_t {
	uint32_t block_xdim;  // 4 to 12 for ASTC-2D, 3 to 6 for ASTC-3D
	uint32_t block_ydim;  // 4 to 12 for ASTC-2D, 3 to 6 for ASTC-3D
	uint32_t block_zdim;  // 3 to 6 for ASTC-3D, not used for ASTC-2D

	uint32_t min_block_ydim_8; // block-Y-dimension or 8, whichever is smaller.
	uint32_t block_xpydim;     // block_xdim + block_ydim

	uint32_t block_xdim_m3_x16;  // (block_xdim-3)*16
	uint32_t block_ydim_m3_x16;  // (block_ydim-3)*16

	// partition-data pointers : astcsd_partition_data
	// 2D: astd_partition_data + (large_block ? 16 : 0):
	//   pointer0 : +0
	//   pointer1 : +64
	//   pointer2 : +65600
	//   pointer3 : +131136
	// 3D: astc_partition_data + (large_block ? 48 : 0) + 196736
	//   pointer0 : +0
	//   pointer1 : +128
	//   pointer2 : +131200
	//   pointer3 : +262272

	const uint8_t *partptrs[4];

	intptr_t src_row_stride;
	intptr_t src_layer_stride;

	intptr_t dst_row_stride;
	intptr_t dst_layer_stride;

	uint32_t texture_xdim;
	uint32_t texture_ydim;
	uint32_t texture_zdim;

	uint32_t blocks_per_row;
	uint32_t shingled_blocks_per_row;

	uint32_t block_rows_per_layer;
	uint32_t block_layers;

	uint8_t *dst_start_addr;
	const uint8_t *src_start_addr;

	void (*decode_row_func)(
		uint8_t *dst,
		const uint8_t *src,
		const struct astc_simd_decode_processed_params_t *bd,
		uint32_t dst_rows,
		uint32_t dst_layers);

	uint64_t hdr_error_color;
};

// Decode parameter data structure for ASTC decode. Must be prepared by callers.

struct astc_simd_decode_params_t
{
	// Pointers to start of data areas
	void* dst_start_ptr;
	const void *src_start_ptr;

	// Data strides, in bytes
	intptr_t src_row_stride; // stride for each block row of source ASTC data
	intptr_t src_layer_stride; // stride for each block layer of source ASTC data (3D textures only; may be set to 0 for 2D textures)
	intptr_t dst_row_stride; // stride for each pixel row of output texel data
	intptr_t dst_layer_stride; // stride for each pixel layer of output texel data (3D textures only; may be set to 0 for 2D textures)

	// Pixel resolution of result texture as a whole
	uint32_t xres;
	uint32_t yres;
	uint32_t zres; // must be 1 (not 0) for 2D textures

	// ASTC block dimensions
	uint8_t block_xdim;
	uint8_t block_ydim;
	uint8_t block_zdim; // must be 1 (not 0) for ASTC-2D

	// Decode flags:
	//    bit 0: HDR enable
	//       0 = decode to RGBA 8888 (unorm8 components, 4 bytes per pixel)
	//       1 = decode to RGBA 16:16:16:16 (FP16 components, 8 bytes per pixel)
	//    bit 1: Error color for HDR decode
	//       0 = 0xFFFF (NaN) in each component
	//       1 = Magenta (R,G,B,A=1,0,1,1)
	uint8_t decode_flags;
};


// ---------------------------------------------------------------------
// ASTC SIMD Decoder precomputation function. Must be run exactly
// once before any ASTC decoding is done. 'astc_simd_decode_prepare()'
// and 'astc_simd_decode_full()' will run it if it hasn't previously
// been run.
//
// Not MT-safe.
// ---------------------------------------------------------------------
void astc_simd_decode_precompute(void);

// --------------------------------------------------------------------
// Function to decode a complete ASTC texture in one go.
//
// MT-safe if 'astc_simd_decode_precompute()' has previously been run.
//
// This function will perform the complete decode within a single
// thread. If multithreaded decode of a single texture is desired,
// please use the below 'astc_simd_decode_prepare()' and
// 'astc_simd_decode_row_iterate()' functions.
// --------------------------------------------------------------------
void astc_simd_decode_full(const astc_simd_decode_params_t* inp_params);

// --------------------------------------------------------------------
// Function to preprocess ASTC decoding parameters for a given
// texture into the opaque data structure above to assist fast
// decoding. Note that it doesn't do any internal allocations -
// it relies on the caller to allocate a memory block large enough
// to hold a full 'astc_simd_decode_processed_params_t' data structure.
//
// Its return value is an integer specifying how many ASTC
// block rows to process. This gives the number of times the
// 'astc_simd_decode_row_iterate()' function must be called in
// order to decode the complete texture.
//
// MT-safe if 'astc_simd_decode_precompute()' has previously been run.
// ---------------------------------------------------------------------
uint32_t astc_simd_decode_prepare(
	astc_simd_decode_processed_params_t* processed_params,
	const astc_simd_decode_params_t* inp_params);

// ------------------------------------------------------------------
// Row iteration function for ASTC SIMD decode.
// This function must be called once per row for the decoding of
// a full ASTC texture.
//
// MT-safe: different threads may invoke the function for different
// rows of the same texture in parallel.
// ------------------------------------------------------------------
void astc_simd_decode_row_iterate(
	const astc_simd_decode_processed_params_t* args,
	uint32_t row_index);

// -----------------------------------------------------------------
// Per-block decode functions. All of these are MT-safe.
//
// Each of these functions will decode exactly 1 ASTC block.
// Note that when called as standalone functions, they will
// write somewhat more data per scanline than what the ASTC block
// size indicates that they should. As such, the 'dst_row_stride'
// argument must be set as follows:
//
//  * 2d_8bit:   at least 48 bytes
//  * 2d_16bit:  at least 96 bytes
//  * 3d_8bit:   at least 32 bytes
//  * 3d_16bit:  at least 64 bytes
// -----------------------------------------------------------------
void astc_simd_decode_block_2d_8bit(
	const uint8_t* blk,                   // Input ASTC block
	const astc_simd_decode_processed_params_t* bd, // ASTC-decode parameter block
	uint32_t dst_rows,                    // Number of rows to actually output
	intptr_t dst_row_stride,              // Byte stride for each row to output
	void* dst                             // Memory address to write block to
	);

void astc_simd_decode_block_2d_16bit(
	const uint8_t* blk,                   // Input ASTC block
	const astc_simd_decode_processed_params_t* bd, // ASTC-decode parameter block
	uint32_t dst_rows,                    // Number of rows to actually output
	intptr_t dst_row_stride,              // Byte stride for each row to output
	void* dst                             // Memory address to write block to
	);

void astc_simd_decode_block_3d_8bit(
	const uint8_t* blk,                   // Input ASTC block
	const astc_simd_decode_processed_params_t* bd, // ASTC-decode parameter block
	uint32_t dst_rows,                    // Number of rows to actually output
	intptr_t dst_row_stride,              // Byte stride for each row to output
	uint32_t dst_layers,                  // Number of layers to actually output
	intptr_t dst_layer_stride,            // Byte stride for each layer to output
	void* dst                             // Memory address to write block to
	);



void astc_simd_decode_block_3d_16bit(
	const uint8_t* blk,                   // Input ASTC block
	const astc_simd_decode_processed_params_t* bd, // ASTC-decode parameter block
	uint32_t dst_rows,                    // Number of rows to actually output
	intptr_t dst_row_stride,              // Byte stride for each row to output
	uint32_t dst_layers,                  // Number of layers to actually output
	intptr_t dst_layer_stride,            // Byte stride for each layer to output
	void* dst                             // Memory address to write block to
	);

#ifdef __cplusplus
}
#endif

// ***********************************
// ***********************************
// ***                             ***
// ***   End of public interface   ***
// ***                             ***
// ***********************************
// ***********************************
#endif // matches: #if !defined(ASTC_SIMD_DECODE_H_INCLUDED) && !defined(ASTCSD_RECURSIVE_INCLUDE)


#if defined(ASTC_SIMD_DECODE_IMPLEMENTATION) && !defined(ASTCSD_RECURSIVE_INCLUDE)
#include <string.h>  // for memset()

// ---------------------------------------------------------
//   various macros/decorators to assist optimization and
//   notify the compiler of unusual-but-intended behaviors.
// ---------------------------------------------------------

// Define a switch fallthrough macro. Depending on compiler and
// language version, there exists at least 4 syntaxes for this.

#if (__cplusplus >= 201703L)
	#define ASTCSD_FALLTHROUGH [[fallthrough]]
#elif defined(__clang__) && defined(__has_cpp_attribute)
	#if (__has_cpp_attribute(clang::fallthrough))
		#define ASTCSD_FALLTHROUGH [[clang::fallthrough]]
	#else
		#define ASTCSD_FALLTHROUGH ((void)0)
	#endif
#elif defined(__GNUC__) && __GNUC__ >= 7
	#define ASTCSD_FALLTHROUGH __attribute__((fallthrough))
#else
	#define ASTCSD_FALLTHROUGH ((void)0)
#endif

// Function attributes to indicate that the function
// is "pure" (depends only on input arguments and globals)
// or "const" (depends only on input arguments).
#if defined(__GNUC__)
	#define ASTCSD_ATTR_CONST __attribute__((const))
	#define ASTCSD_ATTR_PURE  __attribute__((pure))
#else
	#define ASTCSD_ATTR_CONST
	#define ASTCSD_ATTR_PURE
#endif

// Function attributes to indicate that the function performs
// deliberate unaligned memory accesses or too-large shift amounts,
// to reduce false positives from UBSAN under gcc/clang (-fsanitize=undefined)
#if defined(__clang__) || (defined(__GNUC__) && __GNUC__ >= 8)
	#define ASTCSD_ALLOW_UNALIGNED __attribute__ ((no_sanitize("alignment")))
	#define ASTCSD_UNCHECKED_SHIFT __attribute__ ((no_sanitize("shift")))
#else
	#define ASTCSD_ALLOW_UNALIGNED
	#define ASTCSD_UNCHECKED_SHIFT
#endif

#define ASTCSD_RESTRICT __restrict__
#define ASTCSD_UNUSED(x) (void)(x)

// *****************************************************************************
// *****************************************************************************
// *****                                                                   *****
// *****         ASTC SIMD Decode Library: platform-specific               *****
// *****         wrappers and macros around SIMD intrinsics                *****
// *****                                                                   *****
// *****************************************************************************
// *****************************************************************************

/*
The following SIMD intrinsic-wrappers are defined:
  ----------
  data types
  ----------

  vec128_t, vec256_t, vec256x2_t : vector data types

  --------------------------------
  memory load/store, data movement
  --------------------------------

  v128_load_unaligned(ptr)  : unaligned loads and stores
  v256_load_unaligned(ptr)
  v128_store_unaligned(ptr, datum)
  v256_store_unaligned(ptr, datum)


  v128_load_32zx(ptr)     : load a 32-bit value from memory, and zero-extend it to 128 bits.
  v128_load_64zx(ptr)     : load a 64-bit value from memory, and zero-extend it to 128 bits.

  v128_load_8rep(ptr )    : load an 8-bit value from memory, and replicate it (16x) to form a 128-bit vector
  v128_load_32rep(ptr)    : load a 32-bit value from memory, and replicate it (4x) to form a 128-bit vector
  v128_load_64rep(ptr)    : load a 64-bit value from memory, and replicate it (2x) to form a 128-bit vector
  v256_load_32rep(ptr)    : load a 32-bit value from memory, and replicate it (8x) to form a 256-bit vector
  v256_load_64rep(ptr)    : load a 64-bit value from memory, and replicate it (4x) to form a 256-bit vector

  v128_store_64low(ptr, datum)  : store 64-bit datum to memory from low half of 128-bit vector
  v128_store_64high(ptr, datum) : store 64-bit datum to memory from high half of 128-bit vector

  v128_set_64zx(src)      : load a 64-bit value from general-purpose register into vector register, and zero-extend it to 128 bits.

  v128_cat_64(high, low)  : create a 128-bit vector by concatenating two 64-bit data items
  v256_cat_128(high, low) : create a 256-bit vector by concatenating two 128-bit vectors

  v128_rep_low_half_of(src)    : replicate the bottom 64 bits of a 128-bit vector to form a new 128-bit vector

  v128_from_low_half_of(src)   : get low 128-bit half of a 256-bit vector
  v128_from_high_half_of(src)  : get high 128-bit half of a 256-bit vector
  u64_from_low_half_of(src)    : get low 64-bit half of a 128-bit vector, as a general-purpose register.

  --------------------
  data reorganizations
  --------------------

  v128_byterev_64(src)   : perform byte-order reversal within each 64-bit lane of a 128-bit vector
  v256_byterev_64(src)   : perform byte-order reversal within each 64-bit lane of a 256-bit vector

  v128_interleave_low_u16(hi, lo)  : Interleave low bytes of 16-bit lanes:
									   for each 16-bit lane, set lane to LO + (HI << 8).
									   LO and HI are both assumed to be 255 or less; undefined behavior otherwise.
									   Maps to the ARM NEON 'TRN1' instruction; on SSE, it maps to a 1-byte left-shift and an OR.

  v128_interleave_low_u32(hi, lo) : Interleave low words of 32-bit lanes:
									  for each 16-bit lane, set lane to LO + (HI << 16).
									  LO and HI are both assumed to be 65535 or less; undefined behavior otherwise.
									  Maps to the ARM NEON 'TRN1' instruction; on SSE, it maps to a 2-byte left-shift and an OR.

  v128_interleave_bits28_13_to_u32(hi, lo) : for each 32-bit lane of the inputs, extract bits 28:13 as a 16-bit value, then
											   concatenate the two values into one result 32-bit value. (This turns out to have
											   very different optimal sequences in SSSE3, AVX2 and NEON.)

  v256_pair_interleave32_and_store(dst, hi, lo) : Interleave 32-bit items from two 256-bit vectors to form a 512-bit vector
													of eight 64-bit elements. This has similar sequences under SSSE3 and NEON
													but a rather different sequence under AVX2.

  v128_interleave2_x816_low(hi, lo)  : Interleave bytes from two 128-bit vectors, return bottom 16 bytes of result.
  v128_interleave2_x816_high(hi, lo) : Interleave bytes from two 128-bit vectors, return top 16 bytes of result.

  v256x2_interleave2_8x32(hi, lo) : Interleave bytes from two 256-bit vectors, returning result as a pair of 32-byte vectors

  -----------------------
  data-dependent permutes
  -----------------------

  v128_permute_8x16(datum, perm) : perform a 16-byte permute, where 'datum' is the data to permute, and 'perm' is the permutation
									 that, for each byte lane, determines which byte lane datum to pick data from. Lane index goes from
									 0 to 15; a large index of 128 or greater results in the lane being filled with 0. Lane indexes
									 in the range 16..127 produce undefined behavior. This corresponds to SSE's PSHUFB or NEON's TBL1.

  v256_permute_8x16_pair(datum, perm) : performs two 16-byte permutes, one for each 128-bit lane of the input vectors.
										  This corresponds to AVX2's VPSHUFB.

  v128_permutez_8x16(datum, perm)      : 16-byte permute, but if the per-byte index is in the range 16..127, then it is taken modulo 16
  v256_permutez_8x16_pair(datum, perm) : and used as byte-index. This requires an extra masking before NEON's TBL1, but not with SSE's PSHUFB.

  ----------
  arithmetic
  ----------

  v128_mul_add_pair_u16(a, b) : interpret a and b as vectors of 8-bit unsigned integers, then perform pairwise multiply-add of each
  v256_mul_add_pair_u16(a, b) : pair of lanes. Values in 'b' are assumed to be in the range 0..127, with elements in a pair assumed to
								  have a sum of 127 or less. This corresponds to PMADDUBSW in SSE/AVX and a 3-instruction sequence in NEON.

  v128_mul_add_pair_s32(a, b) : interpret a and b as vectors of 16-bit signed integers, then perform pairwise multiply-add of each
  v256_mul_add_pair_s32(a, b) : pair of lanes. Values in 'b' are assumed to be in the range 0..127, with elements in a pair assumed to
								  have a sum of 127 or less. This corresponds to PMADDWD in SSE/AVX and a 3-instruction sequence in NEON.

  v128_add_8(a,b), v128_sub_8(a,b)     : Addition, subtraction with 8 bits per component
  v256_add_8(a,b), v256_sub_8(a,b)

  v128_satsub_u8(a, b), v128_satsub_u16(a, b)  : subtraction with unsigned saturation

  v128_avg_u16(a, b)  : 16-bit averaging with rounding :  (a+b+1)>>1
  v256_avg_u16(a, b)

  v128_mul2add1_8(a)  :  For each 8-bit lane, compute 2*a+1
  v256_mul2add1_8(a)  :  For each 8-bit lane, compute 2*a+1

  v128_min_u8(a, b), v128_min_s16(a, b)
  v128_max_u8(a, b), v128_max_u16(a, b)   : integer minmax

  v128_cmpgt_s8(a,b)   : compare greater-than, with 8-bit signed component. Returns 0xFF for each 'true' lane and 0 for each 'false' lane.
  v256_cmpgt_s8(a,b)
  v128_cmpgt_s16(a,b)

  ------------------
  Bit-ops and shifts
  ------------------

  v128_and(a,b), v128_or(a,b), v128_xor(a,b)  : 128-bit bit-operations
  v256_and(a,b), v256_or(a,b), v256_xor(a,b)  : 256-bit bit-operations

  v128_andnot(a, b)  :  bit-operation: A & ~B  (PANDN on x86, BIC on ARM)
  v256_andnot(a, b)

  v128_shlb(datum, shamt)  : shift 128-bit vector left; shift-amount given in bytes
  v128_shrb(datum, shamt)  : shift 128-bit vector right; shift-amount given in bytes

  v128_shri_u16(datum, shamt) : shift right immediate, with 16-bit lanes; shift-amount given in bits
  v256_shri_u16(datum, shamt) : shift right immediate, with 16-bit lanes; shift-amount given in bits
  v128_shli_u16(datum, shamt) : shift left immediate, with 16-bit lanes; shift-amount given in bits
  v256_shli_u16(datum, shamt) : shift left immediate, with 16-bit lanes; shift-amount given in bits

  v128_shl_u16(datum, shamt) : shift left of 16-bit vector-component, with shift-amount from scalar register.

  v128_shr2_u8(d), v256_shr2_u8(d) : 2-bit right-shift, with 8-bit lanes
  v128_shr4_u8(d), v256_shr4_u8(d) : 4-bit right-shift, with 8-bit lanes

  v128_shr6_round_u16(d), v256_shr6_round_u16(f) : 6-bit right-shift with rounding, with 16-bit lanes.
  v128_shr6_round_u32(d), v256_shr6_round_u32(f) : 6-bit right-shift with rounding, with 32-bit lanes.

  -------------------
  other/uncategorized
  -------------------

  v128_zero()   : 128-bit vector of all-0s (generated without memory accesses, using VPXOR on x86 and MOVI on NEON).
  v256_zero()   : 256-bit vector or all-0s

  v128_all1s()  : 128-bit vector of all-1s (no memory access: VCMPEQD on x86, MOVI on NEON)
  v256_all1s()  : 256-bit vector or all-1s

  v256_select8(a, b, c)  : for each byte lane: select b if top bit of a is set, else select c.

  v128_unorm16_to_fp16z(a),  :  Convert UNORM16 to FP16, with round-to-zero rounding.
  v256_unorm16_to_fp16z(a)      Doesn't need to work for values 1,2,3. This is sufficient
                                  for values produced by lerp.

  v128_unorm16_to_fp16z_full(a) : Convert UNORM16 to FP16, with round-to-zero rounding.
                                    Is required to work for values 1,2,3. This is only required
                                    for handling of LDR void-extent blocks. (This requires
                                    extra instructions on SSE/AVX but not on NEON)
*/

// ---------------------------------
//   Common for existing platforms
// ---------------------------------
static inline ASTCSD_ALLOW_UNALIGNED uint32_t u32_load_unaligned(const void *src) { return *(const uint32_t *)src; }
static inline ASTCSD_ALLOW_UNALIGNED uint64_t u64_load_unaligned(const void *src) { return *(const uint64_t *)src; }
static inline ASTCSD_ALLOW_UNALIGNED void u32_store_unaligned(void* dst, uint32_t datum) { *(uint32_t *)dst = datum; }
static inline ASTCSD_ALLOW_UNALIGNED void u64_store_unaligned(void* dst, uint64_t datum) { *(uint64_t *)dst = datum; }

static inline ASTCSD_UNCHECKED_SHIFT uint64_t u64_left_rotate(uint64_t datum, int shamt) { return (datum << shamt) | (datum >> (64-shamt)); }

// --------------------------------------------
//    ARMv8 NEON:  16-byte vector operations
// --------------------------------------------

#if defined(__ARM_NEON__)

	#define ASTCSD_NATIVE_16_BYTE_VECTORS
	#include <arm_neon.h>

	typedef uint8x16_t vec128_t;

	// Functions that map to well-defined AArch64 instructions but do not appear
	// to be easily reachable through NEON intrinsics. These use inline assembly.

	static inline uint64_t u64_byterev(uint64_t p)  { uint64_t res; asm("rev %x0, %x1" : "=r"(res) : "r"(p)); return res; }

	static inline ASTCSD_ATTR_PURE  vec128_t v128_load_32zx(const uint32_t *v)  { vec128_t res; asm("ldr %s0, %1": "=w"(res) : "m"(*v)); return res; }
	static inline ASTCSD_ATTR_PURE  vec128_t v128_load_64zx(const uint64_t *v)  { vec128_t res; asm("ldr %d0, %1": "=w"(res) : "m"(*v)); return res; }
	static inline ASTCSD_ATTR_CONST vec128_t v128_set_64zx(const uint64_t v)    { vec128_t res; asm("fmov %d0, %1": "=w"(res) : "r"(v)); return res; }
	static inline ASTCSD_ATTR_CONST vec128_t v128_set_32zx(const uint32_t v)    { vec128_t res; asm("fmov %s0, %w1": "=w"(res) : "r"(v)); return res; }
	static inline ASTCSD_ATTR_CONST vec128_t v128_zero(void)                      { vec128_t res; asm("movi %0.4s, #0": "=w"(res)); return res; }
	static inline ASTCSD_ATTR_CONST  vec128_t v128_all1s(void)                    { vec128_t res; asm("movi %0.16b, #255": "=w"(res)); return res; }
	static inline ASTCSD_ATTR_CONST  vec128_t v128_val1_8rep(void)    { vec128_t res; asm("movi %0.16b, #1": "=w"(res)); return res; }
	static inline ASTCSD_ATTR_CONST  vec128_t v128_val3_8rep(void)    { vec128_t res; asm("movi %0.16b, #3": "=w"(res)); return res; }
	static inline ASTCSD_ATTR_CONST  vec128_t v128_val15_8rep(void)   { vec128_t res; asm("movi %0.16b, #15": "=w"(res)); return res; }
	static inline ASTCSD_ATTR_CONST  vec128_t v128_val143_8rep(void)  { vec128_t res; asm("movi %0.16b, #143": "=w"(res)); return res; }

	// Functions that map to NEON intrinsics
	static inline vec128_t v128_cat_64(uint64_t high, uint64_t low) { return vreinterpretq_u64_u8(vsetq_lane_u64(high, vreinterpretq_u8_u64(v128_set_64zx(low)), 1)); }
	static inline vec128_t v128_rep_low_half_of(vec128_t src)       { uint64x2_t v = vreinterpretq_u8_u64(src); return vreinterpretq_u64_u8(vzip1q_u64(v, v)); }

	static inline vec128_t v128_cmpgt_s8(vec128_t a, vec128_t b)    { return vreinterpretq_s8_u8(vcgtq_s8(vreinterpretq_u8_s8(a), vreinterpretq_u8_s8(b))); }
	static inline vec128_t v128_cmpgt_s16(vec128_t a, vec128_t b)   { return vreinterpretq_s16_u8(vcgtq_s16(vreinterpretq_u8_s16(a), vreinterpretq_u8_s16(b))); }
	static inline vec128_t v128_cmpeq_16(vec128_t a, vec128_t b)    { return vreinterpretq_u16_u8(vceqq_u16(vreinterpretq_u8_u16(a), vreinterpretq_u8_u16(b))); }
	static inline vec128_t v128_add_16(vec128_t a, vec128_t b)      { return vreinterpretq_u16_u8(vaddq_u16(vreinterpretq_u8_u16(a), vreinterpretq_u8_u16(b))); }
	static inline vec128_t v128_sub_16(vec128_t a, vec128_t b)      { return vreinterpretq_u16_u8(vsubq_u16(vreinterpretq_u8_u16(a), vreinterpretq_u8_u16(b))); }
	static inline vec128_t v128_satsub_u16(vec128_t a, vec128_t b)  { return vreinterpretq_u16_u8(vqsubq_u16(vreinterpretq_u8_u16(a), vreinterpretq_u8_u16(b))); }
	static inline vec128_t v128_min_s16(vec128_t a, vec128_t b)     { return vreinterpretq_s16_u8(vminq_s16(vreinterpretq_u8_s16(a), vreinterpretq_u8_s16(b))); }
	static inline vec128_t v128_max_s16(vec128_t a, vec128_t b)     { return vreinterpretq_s16_u8(vmaxq_s16(vreinterpretq_u8_s16(a), vreinterpretq_u8_s16(b))); }
	static inline vec128_t v128_avg_u16(vec128_t a, vec128_t b)     { return vreinterpretq_u16_u8(vrhaddq_u16(vreinterpretq_u8_u16(a), vreinterpretq_u8_u16(b))); }

	static inline vec128_t v128_mul_add_pair_u16(vec128_t a, vec128_t b)
	{
		uint16x8_t prod_low = vmull_u8(vget_low_u8(a), vget_low_u8(b));
		uint16x8_t prod_high = vmull_high_u8(a, b);
		return vreinterpretq_u16_u8(vpaddq_u16(prod_low, prod_high));
	}

	static inline vec128_t v128_mul_add_pair_s32(vec128_t a, vec128_t b)
	{
		int32x4_t prod_low = vmull_s16(vget_low_s16(vreinterpretq_u8_s16(a)), vget_low_s16(vreinterpretq_u8_s16(b)));
		int32x4_t prod_high = vmull_high_s16(vreinterpretq_u8_s16(a), vreinterpretq_u8_s16(b));
		return vreinterpretq_s32_u8(vpaddq_s32(prod_low, prod_high));
	}

	static inline vec128_t v128_shl_u16(vec128_t v, int shamt)
	{
		return vreinterpretq_u16_u8(vshlq_u16(vreinterpretq_u8_u16(v), vdupq_n_u16(shamt)));
	}

	#define v128_load_8rep(ptr)   vld1q_dup_u8(ptr)
	#define v128_load_32rep(ptr)  vreinterpretq_u32_u8(vld1q_dup_u32(ptr))
	#define v128_load_64rep(ptr)  vreinterpretq_u64_u8(vld1q_dup_u64(ptr))

	#define v128_store_64low(ptr, datum)   vst1q_lane_u64((uint64_t *)(ptr), vreinterpretq_u8_u64(datum), 0)
	#define v128_store_64high(ptr, datum)  vst1q_lane_u64((uint64_t *)(ptr), vreinterpretq_u8_u64(datum), 1)

	#define v128_and(a, b)     vandq_u8(a, b)
	#define v128_or(a, b)      vorrq_u8(a, b)
	#define v128_xor(a, b)     veorq_u8(a, b)
	#define v128_andnot(a, b)  vbicq_u8(a, b)

	#define v128_add_8(a, b)  vaddq_u8(a, b)
	#define v128_sub_8(a, b)  vsubq_u8(a, b)

	#define v128_satsub_u8(a, b)  vqsubq_u8(a, b)

	#define v128_min_u8(a, b) vminq_u8(a, b)
	#define v128_max_u8(a, b) vmaxq_u8(a, b)

	#define v128_permute_8x16(a, b)   vqtbl1q_u8(a, b)
	#define v128_permutez_8x16(a, b)  vqtbl1q_u8(a, vandq_u8(b, v128_val143_8rep()))

	#define v128_byterev_64(a)       vreinterpretq_u64_u8(vrev64q_u8(vreinterpretq_u8_u64(a)))
	#define v128_convert_s32_f32(a)  vreinterpretq_f32_u8(vcvtq_f32_s32(vreinterpretq_u8_s32(a)))

	#define v128_mul2add1_8(a)  vsliq_n_u8(v128_val1_8rep(), (a), 1)

	#define v128_shlb(a, shamt)  vextq_u8(v128_zero(), (a), 16-(shamt))
	#define v128_shrb(a, shamt)  vextq_u8((a), v128_zero(), (shamt))

	#define v128_shri_u16(a, shamt)  vreinterpretq_u16_u8(vshrq_n_u16(vreinterpretq_u8_u16(a), (shamt)))
	#define v128_shri_u32(a, shamt)  vreinterpretq_u32_u8(vshrq_n_u32(vreinterpretq_u8_u32(a), (shamt)))
	#define v128_shli_u16(a, shamt)  vreinterpretq_u16_u8(vshlq_n_u16(vreinterpretq_u8_u16(a), (shamt)))
	#define v128_shli_u32(a, shamt)  vreinterpretq_u32_u8(vshlq_n_u32(vreinterpretq_u8_u32(a), (shamt)))

	#define v128_interleave_low_u16(hi, lo)  vtrn1q_u8(lo, hi)
	#define v128_interleave_low_u32(hi, lo)  vreinterpretq_u16_u8(vtrn1q_u16(vreinterpretq_u8_u16(lo), vreinterpretq_u8_u16(hi)))

	#define v128_shr2_u8(a)  vshrq_n_u8(a, 2)
	#define v128_shr4_u8(a)  vshrq_n_u8(a, 4)

	#define v128_and3_8(a)   v128_and(a, v128_val3_8rep())
	#define v128_and15_8(a)  v128_and(a, v128_val15_8rep())

	#define v128_shr6_round_u16(a)  vreinterpretq_u16_u8(vrshrq_n_u16(vreinterpretq_u8_u16(a), 6))
	#define v128_shr6_round_u32(a)  vreinterpretq_u32_u8(vrshrq_n_u32(vreinterpretq_u8_u32(a), 6))

	#define v128_load_unaligned(ptr)  vld1q_u8((const uint8_t *) (ptr))
	#define v128_store_unaligned(ptr, datum)  vst1q_u8((uint8_t *) (ptr), (datum))

	#define u64_from_low_half_of(vl) vgetq_lane_u64(vreinterpretq_u8_u64(vl), 0)

	static inline vec128_t v128_interleave_bits28_13_to_u32(vec128_t hi, vec128_t lo)
	{
		uint16x8_t lowv = vreinterpretq_u32_u16(vshrq_n_u32(vreinterpretq_u8_u32(lo), 13));
		uint16x8_t highv = vreinterpretq_u32_u16(vshrq_n_u32(vreinterpretq_u8_u32(hi), 13));

		return vreinterpretq_u16_u8(vtrn1q_u16(lowv, highv));
	}

	#define v128_interleave32_low(hi, lo)   vreinterpretq_u32_u8(vzip1q_u32(vreinterpretq_u8_u32(lo), vreinterpretq_u8_u32(hi)))
	#define v128_interleave32_high(hi, lo)  vreinterpretq_u32_u8(vzip2q_u32(vreinterpretq_u8_u32(lo), vreinterpretq_u8_u32(hi)))

	#define v128_interleave2_8x16_low(hi, lo) vzip1q_u8(lo, hi)
	#define v128_interleave2_8x16_high(hi, lo) vzip2q_u8(lo, hi)

	static inline vec128_t v128_unorm16_to_fp16z(vec128_t inp)
	{
		// compute X AND NOT (X>>11) in order to construct a value
		// that will be rounded downwards when converting to FP16.
		uint16x8_t v0 = vreinterpretq_u8_u16(inp);
		uint16x8_t v1 = vbicq_u16(v0, vshrq_n_u16(v0, 11));
		uint16x8_t v2 = vreinterpretq_f16_u16(vcvtq_n_f16_u16(v1, 16));
		// special-case adjustment for 0xFFFF input value, which must
		// be converted to 1.0 .
		uint16x8_t v3 = vceqq_u16(v0, vreinterpretq_u8_u16(v128_all1s()));
		uint16x8_t v4 = vsubq_u16(v2, v3);
		return vreinterpretq_u16_u8(v4);
	}

	#define v128_unorm16_to_fp16z_full(x) v128_unorm16_to_fp16z(x)

#endif

// ---------------------------------------------
//    Intel SSSE3:  16-byte vector operations
// ---------------------------------------------

#if defined(__SSSE3__)

	#define ASTCSD_NATIVE_16_BYTE_VECTORS
	#include <immintrin.h>

	typedef __m128i vec128_t;

	static inline uint64_t u64_byterev(uint64_t p)  { asm("bswapq %q0" : "+r"(p)); return p; }
	static inline vec128_t v128_load_32zx(const uint32_t *src) { return _mm_castps_si128(_mm_load_ss((const float*)src)); }
	static inline vec128_t v128_load_64zx(const uint64_t *src) { return _mm_castpd_si128(_mm_load_sd((const double*)src)); }

	#define v128_set_64zx(src)  _mm_cvtsi64_si128((int64_t)(src))
	#define v128_set_32zx(src)  _mm_cvtsi32_si128((int)(src))

	// for v128_store_64low, use the Intel intrinsic on platforms that
	// have it (clang 8.0+, gcc 9.0+) and inline assembly otherwise
	#if (defined(__clang__) && (__clang_major__ >= 8)) || ((defined(__GNUC__) && (__GNUC__ >= 9) && !defined(__clang__)))
		#define v128_store_64low(addr, datum)   _mm_storeu_si64((addr), (datum))
	#else
		void v128_store_64low(void *addr, vec128_t datum)
		{
			asm("movq %1,%0": "=m"(*(uint64_t *)addr) : "x"(datum));
		}
	#endif

	#define v128_store_64high(addr, datum)  _mm_storeh_pd((double *)(addr), _mm_castsi128_pd(datum))

	#define v128_zero()  _mm_setzero_si128()

	#define v128_and(a, b)    _mm_and_si128(a, b)
	#define v128_or( a, b)    _mm_or_si128( a, b)
	#define v128_xor(a, b)    _mm_xor_si128(a, b)
	#define v128_andnot(a, b) _mm_andnot_si128(b, a)

	#define v128_min_u8(a, b)   _mm_min_epu8(a, b)
	#define v128_min_s16(a, b)  _mm_min_epi16(a, b)
	#define v128_max_u8(a, b)   _mm_max_epu8(a, b)
	#define v128_max_s16(a, b)  _mm_max_epi16(a, b)
	#define v128_satsub_u8(a, b)   _mm_subs_epu8(a, b)
	#define v128_satsub_u16(a, b)  _mm_subs_epu16(a, b)

	#define v128_avg_u16(a, b) _mm_avg_epu16(a, b)

	#define v128_add_8(a, b)   _mm_add_epi8(a, b)
	#define v128_add_16(a, b)  _mm_add_epi16(a, b)
	#define v128_sub_8(a, b)   _mm_sub_epi8(a, b)
	#define v128_sub_16(a, b)  _mm_sub_epi16(a, b)
	#define v128_cmpgt_s8(a, b)  _mm_cmpgt_epi8(a, b)
	#define v128_cmpgt_s16(a, b) _mm_cmpgt_epi16(a, b)
	#define v128_cmpeq_16(a, b)  _mm_cmpeq_epi16(a, b)

	#define v128_cat_64(high, low) _mm_set_epi64x(high, low)

	#define v128_rep_low_half_of(v) _mm_shuffle_epi32(v, 0x44)
	#define u64_from_low_half_of(vl)  _mm_cvtsi128_si64(vl)

	#define v128_permute_8x16(a, b) _mm_shuffle_epi8(a, b)
	#define v128_permutez_8x16(a, b) _mm_shuffle_epi8(a, b)

	#define v128_shlb(a, b) _mm_slli_si128(a, b)
	#define v128_shrb(a, b) _mm_srli_si128(a, b)

	#define v128_shri_u16(a, b) _mm_srli_epi16(a, b)
	#define v128_shri_u32(a, b) _mm_srli_epi32(a, b)
	#define v128_shli_u16(a, b) _mm_slli_epi16(a, b)
	#define v128_shli_u32(a, b) _mm_slli_epi32(a, b)

	#define v128_mul_add_pair_u16(a, b) _mm_maddubs_epi16(a, b)
	#define v128_mul_add_pair_s32(a, b) _mm_madd_epi16(a, b)

	#define v128_load_unaligned(addr)         _mm_loadu_si128((const __m128i *) (addr))
	#define v128_store_unaligned(addr, datum) _mm_storeu_si128((__m128i *) (addr) , datum)

	#define v128_convert_s32_f32(a)  _mm_castps_si128(_mm_cvtepi32_ps(a))

	static inline ASTCSD_ATTR_CONST vec128_t v128_all1s(void)
	{
		return _mm_set1_epi64x(-1);
	}


	static inline vec128_t v128_byterev_64(vec128_t a)
	{
		static const union { uint8_t s[16]; vec128_t v; } bswap_vec = {{ 7,6,5,4,3,2,1,0, 15,14,13,12,11,10,9,8 }};
		return _mm_shuffle_epi8(a, bswap_vec.v);
	}

	static inline vec128_t v128_shr2_u8(vec128_t a)
	{
		static const union { uint8_t s[16]; vec128_t v; } mask_vec = {{ 63,63,63,63,63,63,63,63, 63,63,63,63,63,63,63,63 }};
		return _mm_and_si128(_mm_srli_epi16(a, 2), mask_vec.v);
	}

	static inline vec128_t v128_shr4_u8(vec128_t a)
	{
		static const union { uint8_t s[16]; vec128_t v; } mask_vec = {{ 15,15,15,15,15,15,15,15, 15,15,15,15,15,15,15,15 }};
		return _mm_and_si128(_mm_srli_epi16(a, 4), mask_vec.v);
	}

	static inline vec128_t v128_mul2add1_8(vec128_t p)
	{
		static const union { uint8_t s[16]; vec128_t v; } addend_vec = {{ 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1 }};
		return _mm_add_epi8(_mm_add_epi8(p, p), addend_vec.v);
	}

	static inline vec128_t v128_and3_8(vec128_t p)
	{
		static const union { uint8_t s[16]; vec128_t v; } val3_vec = {{ 3,3,3,3,3,3,3,3, 3,3,3,3,3,3,3,3 }};
		return _mm_and_si128(p, val3_vec.v);
	}

	static inline vec128_t v128_and15_8(vec128_t p)
	{
		static const union { uint8_t s[16]; vec128_t v; } val15_vec = {{ 15,15,15,15,15,15,15,15, 15,15,15,15,15,15,15,15 }};
		return _mm_and_si128(p, val15_vec.v);
	}

	static inline vec128_t v128_shl_u16(vec128_t p, int shamt)
	{
		vec128_t v0 = v128_set_32zx(shamt);
		return _mm_sll_epi16(p, v0);
	}

	static inline vec128_t v128_shr6_round_u16(vec128_t v)
	{
		static const union { uint16_t s[8]; vec128_t v; } val32_vec = {{ 32,32,32,32,32,32,32,32 }};
		return _mm_srli_epi16(_mm_add_epi16(v, val32_vec.v), 6);
	}

	static inline vec128_t v128_shr6_round_u32(vec128_t v)
	{
		static const union { uint32_t s[4]; vec128_t v; } val32_vec = {{ 32,32,32,32 }};
		return _mm_srli_epi32(_mm_add_epi32(v, val32_vec.v), 6);
	}


	static inline vec128_t v128_interleave_low_u16(vec128_t hi, vec128_t lo)
	{
		return _mm_or_si128(_mm_slli_si128(hi, 1), lo);
	}

	static inline vec128_t v128_interleave_low_u32(vec128_t hi, vec128_t lo)
	{
		static const union { uint32_t s[4]; vec128_t v; } val32_vec = {{ 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF }};

		return _mm_or_si128(_mm_slli_epi32(hi, 16), _mm_and_si128(lo, val32_vec.v));
	}

	#if defined(__GNUC__) && defined(__AVX2__)
		// broadcast instructions are AVX-only.
		static inline vec128_t v128_load_8rep( const uint8_t *src)   { vec128_t res; asm("vpbroadcastb %1,%0" :"=x"(res) : "m"(*src)); return res; }
		static inline vec128_t v128_load_32rep(const uint32_t *src)  { return _mm_castps_si128(_mm_broadcast_ss((const float *)src)); }
	#else
		static inline vec128_t v128_load_8rep( const uint8_t *src)   { return _mm_shuffle_epi32(v128_set_32zx(*src * 0x01010101u), 0); }
		static inline vec128_t v128_load_32rep(const uint32_t *src)  { return _mm_shuffle_epi32(v128_load_32zx(src), 0); }
	#endif

	static inline vec128_t v128_load_64rep(const uint64_t *src)  { return _mm_shuffle_epi32(v128_load_64zx(src), 0x44); }

	static inline vec128_t v128_interleave_bits28_13_to_u32(vec128_t hi, vec128_t lo)
	{
		static const union { uint32_t s[4]; vec128_t v; } mask = {{ 0xffff, 0xffff, 0xffff, 0xffff }};

		lo = _mm_srli_epi32(lo, 13);
		hi = _mm_slli_epi32(hi, 3);

		lo = _mm_and_si128(lo, mask.v);
		hi = _mm_andnot_si128(mask.v, hi);

		return _mm_or_si128(lo, hi);
	}


	#define v128_interleave32_low(hi, lo)   _mm_unpacklo_epi32(lo, hi)
	#define v128_interleave32_high(hi, lo)  _mm_unpackhi_epi32(lo, hi)

	#define v128_interleave2_8x16_low(hi, lo) _mm_unpacklo_epi8(lo, hi)
	#define v128_interleave2_8x16_high(hi, lo) _mm_unpackhi_epi8(lo, hi)

	static inline vec128_t v128_unorm16_to_fp16z(vec128_t inp)
	{
		static const union { int16_t s[8]; vec128_t v; } lowmask = {{ -1, 0,-1, 0, -1, 0,-1, 0 }};
		static const union { uint32_t s[16]; vec128_t v; } fpvmask = {{ 0x07ffe000, 0x07ffe000, 0x07ffe000, 0x07ffe000 }};
		vec128_t high_lanes = _mm_srli_epi32(inp, 16);
		vec128_t low_lanes = _mm_and_si128(inp, lowmask.v);

		// Convert to FP32. This will gives us the bits we want in bits 28:13.
		high_lanes = _mm_castps_si128(_mm_cvtepi32_ps(high_lanes));
		low_lanes  = _mm_castps_si128(_mm_cvtepi32_ps(low_lanes));
		high_lanes = _mm_slli_epi32(_mm_and_si128(high_lanes, fpvmask.v), 3);
		low_lanes  = _mm_srli_epi32(_mm_and_si128(low_lanes, fpvmask.v), 13);
		high_lanes = _mm_sub_epi16(high_lanes, _mm_cmpeq_epi16(inp, v128_all1s()));

		return _mm_add_epi16(high_lanes, low_lanes);
	}

	static inline vec128_t v128_unorm16_to_fp16z_full(vec128_t inp)
	{
		vec128_t temp = v128_unorm16_to_fp16z(inp);
		vec128_t mask = _mm_cmpeq_epi16(_mm_srli_epi16(inp, 2), _mm_setzero_si128());
		return _mm_or_si128(_mm_andnot_si128(mask, temp), _mm_and_si128(mask, _mm_slli_epi16(inp, 8)));
	}

#endif

#if !defined(ASTCSD_NATIVE_16_BYTE_VECTORS)
	// If no usable SIMD vector extensions are enabled, then
	// emit an error message and terminate the compilation.
	#error "Cannot compile without SIMD support (SSSE3, AVX2 or ARMv8 NEON)"
	// For compilers that do not by default stop on an #error directive
	// (e.g. gcc and clang), include the name of a file that cannot exist
	// in order to forcibly terminate the compilation at this point and
	// thereby avoid thousands of lines of useless error messages.
	#include <\\stop::here/>
#endif

// ------------------------------------------
//   Intel AVX2: native 32-byte operations
// ------------------------------------------
#if defined(__AVX2__)

	#define ASTCSD_NATIVE_32_BYTE_VECTORS

	typedef __m256i vec256_t;

	typedef struct _vec256x2_t { vec256_t low; vec256_t high; } vec256x2_t;

	// data assignments
	// for v256_cat_128, use the "_mm256_set_m128i" intrinsic on platforms that support it (clang 6.0+, gcc 8.0+)
	#if (defined(__clang__) && (__clang_major__ >= 6)) || (defined(__GNUC__) && (__GNUC__ >= 8) && !defined(__clang__))
		#define v256_cat_128(high, low)     _mm256_set_m128i(high, low)
	#else
		#define v256_cat_128(high, low)     _mm256_inserti128_si256(_mm256_castsi128_si256(low),(high),1)
	#endif
	#define v256_set_128rep(src)        _mm256_broadcastsi128_si256(src)
	#define v256_zero()                   _mm256_setzero_si256()
	#define v128_from_low_half_of(vl)   _mm256_castsi256_si128((vl))
	#define v128_from_high_half_of(vl)  _mm256_extracti128_si256((vl), 1)
	static inline ASTCSD_ATTR_CONST vec256_t v256_all1s(void) { return _mm256_set1_epi64x(-1); } // variant that placates valgrind

	// data loads and stores
	#define v256_load_unaligned(addr)          _mm256_loadu_si256((const __m256i *) (addr))
	#define v256_store_unaligned(addr, datum)  _mm256_storeu_si256((__m256i *) (addr) , datum)

	#if defined(__GNUC__)
		static inline vec256_t v256_load_32rep(const uint32_t *src) { vec256_t res; asm("vpbroadcastd %1,%0" :"=x"(res) : "m"(*src)); return res; }
		static inline vec256_t v256_load_64rep(const uint64_t *src) { vec256_t res; asm("vpbroadcastq %1,%0" :"=x"(res) : "m"(*src)); return res; }
	#else
		static inline vec256_t v256_load_32rep(const uint32_t *src) { return _mm256_castps_si256(_mm256_broadcast_ss((const float *)src)); }
		static inline vec256_t v256_load_64rep(const uint64_t *src) { return _mm256_castpd_si256(_mm256_broadcast_sd((const double *)src)); }
	#endif

	// data rearrangements
	#define v256_permute_8x16_pair(a, b)  _mm256_shuffle_epi8(a, b)
	#define v256_permutez_8x16_pair(a, b) _mm256_shuffle_epi8(a, b)
	static inline vec256_t v256_interleave_low_u16(vec256_t hi, vec256_t lo) { return _mm256_or_si256(_mm256_slli_si256(hi, 1), lo); }
	static inline vec256_t v256_interleave_low_u32(vec256_t hi, vec256_t lo)
	{
		static const union { uint32_t s[8]; vec256_t v; }
		val32_vec = {{ 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF }};
		return _mm256_or_si256(_mm256_slli_epi32(hi, 16), _mm256_and_si256(lo, val32_vec.v));
	}

	static inline vec256_t v256_byterev_64(vec256_t p)
	{
		static const union { uint8_t s[32]; vec256_t v; } bswap_vec = {{ 7,6,5,4,3,2,1,0, 15,14,13,12,11,10,9,8,  7,6,5,4,3,2,1,0, 15,14,13,12,11,10,9,8 }};
		return _mm256_shuffle_epi8(p, bswap_vec.v);
	}

	// shift operations
	#define v256_shri_u16(a, b)  _mm256_srli_epi16(a, b)
	#define v256_shri_u32(a, b)  _mm256_srli_epi32(a, b)
	#define v256_shli_u16(a, b)  _mm256_slli_epi16(a, b)
	#define v256_shli_u32(a, b)  _mm256_slli_epi32(a, b)
	#define v256_shlb_pair(a, b) _mm256_slli_si256(a, b)
	#define v256_shrb_pair(a, b) _mm256_srli_si256(a, b)

	// bitwise operations
	#define v256_and(a, b)    _mm256_and_si256(a, b)
	#define v256_or( a, b)    _mm256_or_si256( a, b)
	#define v256_xor(a, b)    _mm256_xor_si256(a, b)
	#define v256_andnot(a, b) _mm256_andnot_si256(b, a)

	// add/subtract
	#define v256_add_8(a, b)       _mm256_add_epi8(a, b)
	#define v256_add_16(a, b)      _mm256_add_epi16(a, b)
	#define v256_sub_8(a, b)       _mm256_sub_epi8(a, b)
	#define v256_sub_16(a, b)      _mm256_sub_epi16(a, b)
	#define v256_satsub_u8(a, b)   _mm256_subs_epu8(a, b)
	#define v256_satsub_u16(a, b)  _mm256_subs_epu16(a, b)
	#define v256_avg_u16(a, b)     _mm256_avg_epu16(a, b)

	// compare and minmax
	#define v256_min_u8(a, b)    _mm256_min_epu8(a, b)
	#define v256_min_s16(a, b)   _mm256_min_epi16(a, b)
	#define v256_max_u8(a, b)    _mm256_max_epu8(a, b)
	#define v256_max_s16(a, b)   _mm256_max_epi16(a, b)
	#define v256_cmpgt_s8(a, b)  _mm256_cmpgt_epi8(a, b)
	#define v256_cmpeq_16(a, b)  _mm256_cmpeq_epi16(a, b)

	// other important operations
	#define v256_mul_add_pair_u16(a, b) _mm256_maddubs_epi16(a, b)
	#define v256_mul_add_pair_s32(a, b) _mm256_madd_epi16(a, b)
	#define v256_convert_s32_f32(a)     _mm256_castps_si256(_mm256_cvtepi32_ps(a))

	#define v256_select8(a, b, c) _mm256_blendv_epi8(c, b, a)

	// specializations
	static inline vec256_t v256_shr2_u8(vec256_t a)
	{
		static const union { uint8_t s[32]; vec256_t v; } mask_vec = {{ 63,63,63,63,63,63,63,63, 63,63,63,63,63,63,63,63, 63,63,63,63,63,63,63,63, 63,63,63,63,63,63,63,63 }};
		return _mm256_and_si256(_mm256_srli_epi16(a, 2), mask_vec.v);
	}

	static inline vec256_t v256_shr4_u8(vec256_t a)
	{
		static const union { uint8_t s[32]; vec256_t v; } mask_vec = {{ 15,15,15,15,15,15,15,15, 15,15,15,15,15,15,15,15, 15,15,15,15,15,15,15,15, 15,15,15,15,15,15,15,15 }};
		return _mm256_and_si256(_mm256_srli_epi16(a, 4), mask_vec.v);
	}

	static inline vec256_t v256_shr6_round_u16(vec256_t v)
	{
		static const union { uint16_t s[16]; vec256_t v; } val32_vec = {{ 32,32,32,32,32,32,32,32, 32,32,32,32,32,32,32,32 }};
		return _mm256_srli_epi16(_mm256_add_epi16(v, val32_vec.v), 6);
	}

	static inline vec256_t v256_shr6_round_u32(vec256_t v)
	{
		static const union { uint32_t s[8]; vec256_t v; } val32_vec = {{ 32,32,32,32,32,32,32,32 }};
		return _mm256_srli_epi32(_mm256_add_epi32(v, val32_vec.v), 6);
	}

	static inline vec256_t v256_mul2add1_8(vec256_t p)
	{
		static const union { uint8_t s[32]; vec256_t v; } addend_vec = {{ 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1 }};
		return _mm256_add_epi8(_mm256_add_epi8(p, p), addend_vec.v);
	}

	static inline vec256_t v256_and3_8(vec256_t p)
	{
		static const union { uint8_t s[32]; vec256_t v; } val3_vec = {{ 3,3,3,3,3,3,3,3, 3,3,3,3,3,3,3,3, 3,3,3,3,3,3,3,3, 3,3,3,3,3,3,3,3 }};
		return _mm256_and_si256(p, val3_vec.v);
	}

	static inline vec256_t v256_and15_8(vec256_t p)
	{
		static const union { uint8_t s[32]; vec256_t v; } val15_vec = {{ 15,15,15,15,15,15,15,15, 15,15,15,15,15,15,15,15,  15,15,15,15,15,15,15,15, 15,15,15,15,15,15,15,15 }};
		return _mm256_and_si256(p, val15_vec.v);
	}

	static inline vec256_t v256_interleave_bits28_13_to_u32(vec256_t hi, vec256_t lo)
	{
		lo = _mm256_srli_epi32(lo, 13);
		hi = _mm256_slli_epi32(hi, 3);
		return _mm256_blend_epi16(lo, hi, 0xaa);
	}

	// interleave 32-bit items from two 256-bit vectors to form a 512-bit vector of 64-bit elements.
	static inline void v256_pair_interleave32_and_store(void* dst, vec256_t hi, vec256_t lo)
	{
		vec256_t il0 = _mm256_unpacklo_epi32(lo, hi);
		vec256_t il1 = _mm256_unpackhi_epi32(lo, hi);

		vec256_t r0 = _mm256_inserti128_si256(il0, _mm256_castsi256_si128(il1), 1);
		vec256_t r1 = _mm256_blend_epi32(il1, _mm256_castsi128_si256(_mm256_extracti128_si256(il0, 1)), 0xf);

		vec256_t *resp = (vec256_t *)dst;
		v256_store_unaligned(resp, r0);
		v256_store_unaligned(resp+1, r1);
	}

	// interleave 8-bit items from two 256-bit vectors to form a pair of 256-bit vectors
	static inline vec256x2_t v256x2_interleave2_8x32(vec256_t high, vec256_t low)
	{
		vec256x2_t res;
		vec256_t p = _mm256_unpacklo_epi8(low, high); // bits 15:0, 47:32
		vec256_t q = _mm256_unpackhi_epi8(low, high); // bits 31:16, 63:48
		vec256_t r = _mm256_castsi128_si256(_mm256_extracti128_si256(p, 1)); // bytes 47:32
		res.low = _mm256_inserti128_si256(p, _mm256_castsi256_si128(q), 1);
		res.high = _mm256_blend_epi32(r, q, 0xf0);
		return res;
	}


	// Convert a vector fo UNORM16 values to FP16, with round-to-zero rounding.
	// Does not necessarily convert values 1,2,3 correctly, however since these
	// values do not occur as result of lerp, this is not considered a problem.
	static inline vec256_t v256_unorm16_to_fp16z(vec256_t inp)
	{

		static const union { int16_t s[16]; vec256_t v; } lowmask  = {{ -1, 0,-1, 0, -1, 0,-1, 0, -1, 0,-1, 0, -1, 0,-1, 0 }};
		static const union { uint32_t s[16]; vec256_t v; } fpvmask = {{ 0x07ffe000, 0x07ffe000, 0x07ffe000, 0x07ffe000, 0x07ffe000, 0x07ffe000, 0x07ffe000, 0x07ffe000 }};

		vec256_t high_lanes = _mm256_srli_epi32(inp, 16);
		vec256_t low_lanes = _mm256_and_si256(inp, lowmask.v);


		// Convert to FP32. This will gives us the bits we want in bits 28:13.
		high_lanes = _mm256_castps_si256(_mm256_cvtepi32_ps(high_lanes));
		low_lanes  = _mm256_castps_si256(_mm256_cvtepi32_ps(low_lanes));
		high_lanes = _mm256_slli_epi32(_mm256_and_si256(high_lanes, fpvmask.v), 3);
		low_lanes  = _mm256_srli_epi32(_mm256_and_si256(low_lanes, fpvmask.v), 13);

		high_lanes = _mm256_sub_epi16(high_lanes, _mm256_cmpeq_epi16(inp, v256_all1s()));

		return _mm256_add_epi16(high_lanes, low_lanes);
	}

#endif

// -----------------------------------------------------------------------
//   If primitives for native 256-bit vectors have not been provided,
//   then construct them from pairs of native 128-bit vector operations.
//   This applies to SSSE3 and NEON, neither of which have 256-bit
//   native vector operations.
// -----------------------------------------------------------------------

#if !defined(ASTCSD_NATIVE_32_BYTE_VECTORS)

	typedef struct vec256_t_   { vec128_t low; vec128_t high; } vec256_t;
	typedef struct vec256x2_t_ { vec256_t low; vec256_t high; } vec256x2_t;

	// data assignments
	static inline vec256_t v256_set_128rep(vec128_t v)               { vec256_t res; res.low = v;   res.high = v;    return res; }
	static inline vec256_t v256_cat_128(vec128_t high, vec128_t low) { vec256_t res; res.low = low; res.high = high; return res; }
	static inline vec256_t v256_zero(void)                           { return v256_set_128rep(v128_zero()); }
	static inline vec256_t v256_all1s(void)                          { return v256_set_128rep(v128_all1s()); }
	#define v128_from_low_half_of(datum)   ((datum).low)
	#define v128_from_high_half_of(datum)  ((datum).high)

	// data loads
	static inline vec256_t v256_load_unaligned(const void *ptr)
	{
		const uint8_t *p = (const uint8_t *)ptr;
		vec256_t res;
		res.low  = v128_load_unaligned((const vec128_t *)p);
		res.high = v128_load_unaligned((const vec128_t *)(p+16));
		return res;
	}

	static inline void v256_store_unaligned(void *ptr, vec256_t datum)
	{
		uint8_t *p = (uint8_t *)ptr;
		v128_store_unaligned((vec128_t *)p, datum.low);
		v128_store_unaligned((vec128_t *)(p+16), datum.high);
	}

	static inline vec256_t v256_load_32rep(const uint32_t *v)        { return v256_set_128rep(v128_load_32rep(v)); }
	static inline vec256_t v256_load_64rep(const uint64_t *v)        { return v256_set_128rep(v128_load_64rep(v)); }

	// data rearrangments
	static inline vec256_t v256_permute_8x16_pair(vec256_t a, vec256_t b)  { vec256_t res; res.low = v128_permute_8x16(a.low, b.low);  res.high = v128_permute_8x16(a.high, b.high);  return res; }
	static inline vec256_t v256_permutez_8x16_pair(vec256_t a, vec256_t b) { vec256_t res; res.low = v128_permutez_8x16(a.low, b.low); res.high = v128_permutez_8x16(a.high, b.high); return res; }
	static inline vec256_t v256_interleave_low_u16(vec256_t a, vec256_t b) { vec256_t res; res.low = v128_interleave_low_u16(a.low, b.low); res.high = v128_interleave_low_u16(a.high, b.high); return res; }
	static inline vec256_t v256_interleave_low_u32(vec256_t a, vec256_t b) { vec256_t res; res.low = v128_interleave_low_u32(a.low, b.low); res.high = v128_interleave_low_u32(a.high, b.high); return res; }
	static inline vec256_t v256_byterev_64(vec256_t a)                     { vec256_t res; res.low = v128_byterev_64(a.low); res.high = v128_byterev_64(a.high); return res; }

	// shift operations
	#define v256_shri_u16(a, shamt) (vec256_t) { v128_shri_u16((a).low, shamt), v128_shri_u16((a).high, shamt) }
	#define v256_shri_u32(a, shamt) (vec256_t) { v128_shri_u32((a).low, shamt), v128_shri_u32((a).high, shamt) }
	#define v256_shli_u16(a, shamt) (vec256_t) { v128_shli_u16((a).low, shamt), v128_shli_u16((a).high, shamt) }
	#define v256_shli_u32(a, shamt) (vec256_t) { v128_shli_u32((a).low, shamt), v128_shli_u32((a).high, shamt) }
	#define v256_shlb_pair(a, shamt) (vec256_t) { v128_shlb((a).low, shamt), v128_shlb((a).high, shamt) }
	#define v256_shrb_pair(a, shamt) (vec256_t) { v128_shrb((a).low, shamt), v128_shrb((a).high, shamt) }

	// bitwise operations
	static inline vec256_t v256_and(vec256_t a, vec256_t b)     { vec256_t res; res.low = v128_and(a.low, b.low); res.high = v128_and(a.high, b.high); return res; }
	static inline vec256_t v256_or( vec256_t a, vec256_t b)     { vec256_t res; res.low = v128_or(a.low, b.low); res.high = v128_or(a.high, b.high); return res; }
	static inline vec256_t v256_xor(vec256_t a, vec256_t b)     { vec256_t res; res.low = v128_xor(a.low, b.low); res.high = v128_xor(a.high, b.high); return res; }
	static inline vec256_t v256_andnot(vec256_t a, vec256_t b)  { vec256_t res; res.low = v128_andnot(a.low, b.low); res.high = v128_andnot(a.high, b.high); return res; }

	// add/subtract
	static inline vec256_t v256_add_8(vec256_t a, vec256_t b)       { vec256_t res; res.low = v128_add_8(a.low, b.low); res.high = v128_add_8(a.high, b.high); return res; }
	static inline vec256_t v256_add_16(vec256_t a, vec256_t b)      { vec256_t res; res.low = v128_add_16(a.low, b.low); res.high = v128_add_16(a.high, b.high); return res; }
	static inline vec256_t v256_sub_8(vec256_t a, vec256_t b)       { vec256_t res; res.low = v128_sub_8(a.low, b.low); res.high = v128_sub_8(a.high, b.high); return res; }
	static inline vec256_t v256_sub_16(vec256_t a, vec256_t b)      { vec256_t res; res.low = v128_sub_16(a.low, b.low); res.high = v128_sub_16(a.high, b.high); return res; }
	static inline vec256_t v256_satsub_u8(vec256_t a, vec256_t b)   { vec256_t res; res.low = v128_satsub_u8(a.low, b.low); res.high = v128_satsub_u8(a.high, b.high); return res; }
	static inline vec256_t v256_satsub_u16(vec256_t a, vec256_t b)  { vec256_t res; res.low = v128_satsub_u16(a.low, b.low); res.high = v128_satsub_u16(a.high, b.high); return res; }
	static inline vec256_t v256_avg_u16(vec256_t a, vec256_t b)     { vec256_t res; res.low = v128_avg_u16(a.low, b.low); res.high = v128_avg_u16(a.high, b.high); return res; }

	// compare and minmax
	static inline vec256_t v256_min_u8(vec256_t a, vec256_t b)    { vec256_t res; res.low = v128_min_u8( a.low, b.low); res.high = v128_min_u8( a.high, b.high); return res; }
	static inline vec256_t v256_min_s16(vec256_t a, vec256_t b)   { vec256_t res; res.low = v128_min_s16(a.low, b.low); res.high = v128_min_s16(a.high, b.high); return res; }
	static inline vec256_t v256_max_u8(vec256_t a, vec256_t b)    { vec256_t res; res.low = v128_max_u8( a.low, b.low); res.high = v128_max_u8( a.high, b.high); return res; }
	static inline vec256_t v256_max_s16(vec256_t a, vec256_t b)   { vec256_t res; res.low = v128_max_s16(a.low, b.low); res.high = v128_max_s16(a.high, b.high); return res; }
	static inline vec256_t v256_cmpgt_s8(vec256_t a, vec256_t b)  { vec256_t res; res.low = v128_cmpgt_s8(a.low, b.low); res.high = v128_cmpgt_s8(a.high, b.high); return res; }
	static inline vec256_t v256_cmpeq_16(vec256_t a, vec256_t b)  { vec256_t res; res.low = v128_cmpeq_16(a.low, b.low); res.high = v128_cmpeq_16(a.high, b.high); return res; }

	// other important operations
	static inline vec256_t v256_convert_s32_f32(vec256_t a)               { vec256_t res; res.low = v128_convert_s32_f32(a.low); res.high = v128_convert_s32_f32(a.high); return res; }
	static inline vec256_t v256_mul_add_pair_u16(vec256_t a, vec256_t b)  { vec256_t res; res.low = v128_mul_add_pair_u16(a.low, b.low); res.high = v128_mul_add_pair_u16(a.high, b.high); return res; }
	static inline vec256_t v256_mul_add_pair_s32(vec256_t a, vec256_t b)  { vec256_t res; res.low = v128_mul_add_pair_s32(a.low, b.low); res.high = v128_mul_add_pair_s32(a.high, b.high); return res; }


	// specializations
	static inline vec256_t v256_shr2_u8(vec256_t a)        { vec256_t res; res.low = v128_shr2_u8(a.low); res.high = v128_shr2_u8(a.high); return res; }
	static inline vec256_t v256_shr4_u8(vec256_t a)        { vec256_t res; res.low = v128_shr4_u8(a.low); res.high = v128_shr4_u8(a.high); return res; }
	static inline vec256_t v256_shr6_round_u16(vec256_t a) { vec256_t res; res.low = v128_shr6_round_u16(a.low); res.high = v128_shr6_round_u16(a.high); return res; }
	static inline vec256_t v256_shr6_round_u32(vec256_t a) { vec256_t res; res.low = v128_shr6_round_u32(a.low); res.high = v128_shr6_round_u32(a.high); return res; }
	static inline vec256_t v256_mul2add1_8(vec256_t a)     { vec256_t res; res.low = v128_mul2add1_8(a.low); res.high = v128_mul2add1_8(a.high); return res; }

	static inline vec256_t v256_and3_8(vec256_t a)  { vec256_t res; res.low = v128_and3_8(a.low);  res.high = v128_and3_8(a.high); return res;}
	static inline vec256_t v256_and15_8(vec256_t a) { vec256_t res; res.low = v128_and15_8(a.low); res.high = v128_and15_8(a.high); return res;}

	static inline vec256_t v256_interleave_bits28_13_to_u32(vec256_t a, vec256_t b)  { vec256_t res; res.low = v128_interleave_bits28_13_to_u32(a.low, b.low); res.high = v128_interleave_bits28_13_to_u32(a.high, b.high); return res; }

	static inline void v256_pair_interleave32_and_store(void* dst, vec256_t hi, vec256_t lo)
	{
		vec128_t r0 = v128_interleave32_low( hi.low, lo.low);
		vec128_t r1 = v128_interleave32_high(hi.low, lo.low);
		vec128_t r2 = v128_interleave32_low( hi.high, lo.high);
		vec128_t r3 = v128_interleave32_high(hi.high, lo.high);
		vec128_t *res = (vec128_t *)dst;
		v128_store_unaligned(res,   r0);
		v128_store_unaligned(res+1, r1);
		v128_store_unaligned(res+2, r2);
		v128_store_unaligned(res+3, r3);
	}

	static inline vec256x2_t v256x2_interleave2_8x32(vec256_t hi, vec256_t lo)
	{
		vec256x2_t res;
		res.low.low   = v128_interleave2_8x16_low( hi.low, lo.low);
		res.low.high  = v128_interleave2_8x16_high(hi.low, lo.low);
		res.high.low  = v128_interleave2_8x16_low( hi.high, lo.high);
		res.high.high = v128_interleave2_8x16_high(hi.high, lo.high);

		return res;
	}

	static inline vec256_t v256_unorm16_to_fp16z(vec256_t inp)
	{
		vec256_t res;
		res.low = v128_unorm16_to_fp16z(inp.low);
		res.high = v128_unorm16_to_fp16z(inp.high);
		return res;
	}

	static inline vec256_t v256_select8(vec256_t mask, vec256_t b, vec256_t c)
	{
		mask = v256_cmpgt_s8(v256_zero(), mask); // 0xFF if top bit is 1, 0 otherwise.
		c = v256_andnot(c, mask);
		b = v256_and(b, mask);
		return v256_or(b, c);
	}

#endif

// ------------------------------------------------------------------
//   Bit-extract and bit-deposit functions, and supporting code
// ------------------------------------------------------------------

// For bit-extract and bit-deposit, we want to use the x86 PEXT and PDEP
// instructions if available - unless we're running on an AMD Zen1/2, in
// which case these instructions are extremely slow and need to be avoided.
// (These instructions are fast on Zen3, though.)

typedef struct astcsd_u64_pair_t_ { uint64_t v0; uint64_t v1; } astcsd_u64_pair_t;

#if defined(__GNUC__) && (defined(__znver1__) || defined(__znver2__))
	#define ASTCSD_AVOID_SLOW_PEXT_PDEP
#endif

#if defined(__BMI2__) && !defined(ASTCSD_AVOID_SLOW_PEXT_PDEP)
	#define ASTCSD_FAST_PEXT_PDEP
	#define u12_bit_extract(datum, mask)  _pext_u32(datum, mask)
	#define u64_bit_extract(datum, mask)  _pext_u64(datum, mask)
	#define u64_bit_deposit(datum, mask)  _pdep_u64(datum, mask)

	static inline astcsd_u64_pair_t u64_bit_extract_pair( uint64_t datum, uint64_t mask0, uint64_t mask1)
	{
		astcsd_u64_pair_t res;
		res.v0 = _pext_u64(datum, mask0);
		res.v1 = _pext_u64(datum, mask1);
		return res;
	}

	static int astcsd_precompute_bit_extract_and_deposit_tables() { return 0; }
#else
	static uint8_t astcsd_bit_count_tab[256];
	static uint8_t astcsd_bit_deposit_tab[65536];
	static uint8_t astcsd_bit_extract_tab[4096];

	static int astcsd_precompute_bit_extract_and_deposit_tables()
	{
		uint32_t i,j;

		for (i=0;i<256;i++)
		{
			uint32_t ix = i;
			uint32_t counter = 0;
			while(ix) { counter++; ix &= ix-1; }
			astcsd_bit_count_tab[i] = counter;
		}

		for (i=0;i<256;i++)
		{
			for (j=0;j<256;j++)
			{
				// perform 8-bit manual PDEP
				uint32_t k;
				uint32_t pdep_res = 0;
				uint32_t shamt = 0;
				uint32_t mask = 1;
				for (k=0;k<8;k++)
				{
					if (i & mask)
						pdep_res |= (j << shamt) & mask;
					else
						shamt++;
					mask <<= 1;
				}
				astcsd_bit_deposit_tab[256*i+j] = pdep_res & 0xFF;
			}
		}

		for (i=0;i<64;i++)
		{
			for (j=0;j<64;j++)
			{
				uint32_t k;
				uint32_t pext_res = 0;
				uint32_t shamt = 0;
				uint32_t mask = 1;
				for (k=0;k<8;k++)
				{
					if (i & mask)
						pext_res |= (j & mask) >> shamt;
					else
						shamt++;
					mask <<= 1;
				}
				astcsd_bit_extract_tab[i+64*j] = pext_res & 0x3F;
			}
		}
		return 0;
	}

	static uint32_t u12_bit_extract(uint64_t p, uint64_t q)
	{
		uint32_t i;
		uint64_t res = 0;
		uint32_t res_pos = 0;

		p = (p << 6) | (p >> 58); // left-rotate by 6 bits before we start.
		for (i=0;i<2;i++)
		{
			uint64_t mask_bits = q & 0x3F;
			q >>= 6;
			uint64_t data_bits = p & 0xFC0;
			p = (p >> 6) | (p << 58); // right-rotate data bits by 6 bits.
			uint64_t xd_data = astcsd_bit_extract_tab[mask_bits + data_bits];
			res |= xd_data << res_pos;
			res_pos += astcsd_bit_count_tab[mask_bits];
		}
		return res;
	}

	static uint64_t u64_bit_extract(uint64_t datum, uint64_t mask)
	{
		uint32_t i;
		uint64_t res = 0;
		uint32_t res_pos = 0;

		datum = (datum << 6) | (datum >> 58); // left-rotate by 6 bits before we start.
		for (i=0;i<11;i++)
		{
			uint64_t mask_bits = mask & 0x3F;
			mask >>= 6;
			uint64_t data_bits = datum & 0xFC0;
			datum = (datum >> 6) | (datum << 58); // right-rotate data bits by 6 bits.
			uint64_t xd_data = astcsd_bit_extract_tab[mask_bits + data_bits];
			res |= xd_data << res_pos;
			res_pos += astcsd_bit_count_tab[mask_bits];
		}
		return res;
	}

	static astcsd_u64_pair_t u64_bit_extract_pair(uint64_t datum, uint64_t mask0, uint64_t mask1)
	{
		uint32_t i;
		uint64_t res0 = 0, res1 = 0;
		uint32_t res0_pos = 0, res1_pos = 0;

		datum = (datum << 6) | (datum >> 58); // left-rotate by 6 bits before we start.
		for (i=0;i<11;i++)
		{
			uint64_t mask0_bits = mask0 & 0x3F;
			uint64_t mask1_bits = mask1 & 0x3F;
			mask0 >>= 6;
			mask1 >>= 6;
			uint64_t data_bits = datum & 0xFC0;
			datum = (datum >> 6) | (datum << 58); // right-rotate data bits by 6 bits.
			uint64_t xd0_data = astcsd_bit_extract_tab[mask0_bits + data_bits];
			uint64_t xd1_data = astcsd_bit_extract_tab[mask1_bits + data_bits];
			res0 |= xd0_data << res0_pos;
			res1 |= xd1_data << res1_pos;
			res0_pos += astcsd_bit_count_tab[mask0_bits];
			res1_pos += astcsd_bit_count_tab[mask1_bits];
		}
		astcsd_u64_pair_t res;
		res.v0 = res0;
		res.v1 = res1;
		return res;
	}

	// Not fully general version. Only works when mask is a vector of eight equal bytes.
	static uint64_t u64_bit_deposit(uint64_t datum, uint64_t mask)
	{
		uint32_t i;
		uint32_t bitcount_per_byte = astcsd_bit_count_tab[mask & 0xFF];
		const uint8_t *deposit_tab_ptr = astcsd_bit_deposit_tab + ((mask & 0xFF) << 8);

		uint64_t res = 0;
		for (i=0;i<8;i++)
		{
			uint64_t src_bits = datum & 0xFF;
			datum >>= bitcount_per_byte;
			uint64_t deposit_datum = deposit_tab_ptr[src_bits];
			res <<= 8;
			res |= deposit_datum;
		}
		return u64_byterev(res);
	}
#endif

// *****************************************************************************
// *****************************************************************************
// *****                                                                   *****
// *****          End of wrappers/macros for SIMD intrinsics               *****
// *****          and other platform-specific functions                    *****
// *****                                                                   *****
// *****************************************************************************
// *****************************************************************************

// macro to define a value to use for don't-care table elements.
#define _ 0
// ********************************************
//     Lookup tables and constants
// ********************************************


// Table of color-ISE lengths (minus 1, div 2), indexed by (CEM<<2) + (partcnt-1)
// A "length" of 32 is used to indicate an invalid encoding.
static const uint8_t astcsd_cise_count_tbl[256] = {
	 0, 1, 2, 3, 0, 1, 2, 3, 0, 3, 5, 7, 0, 5, 8,32,
	 0, 1, 2, 3, 0, 2, 3, 4, 0, 4, 6, 8, 0, 6,32,32,
	 0, 1, 2, 3, 0, 2, 3, 4, 0, 4, 6, 8, 0, 6,32,32,
	 0, 1, 2, 3, 0, 3, 4, 5, 0, 5, 7,32, 0, 7,32,32,
	 1, 3, 5, 7, 1, 1, 3, 4, 1, 3, 6, 8, 1, 5,32,32,
	 1, 3, 5, 7, 1, 2, 4, 5, 1, 4, 7,32, 1, 6,32,32,
	 1, 3, 5, 7, 1, 2, 4, 5, 1, 4, 7,32, 1, 6,32,32,
	 1, 3, 5, 7, 1, 3, 5, 6, 1, 5, 8,32, 1, 7,32,32,
	 2, 5, 8,32, 2, 1, 2, 4, 2, 3, 5, 8, 2, 5, 8,32,
	 2, 5, 8,32, 2, 2, 3, 5, 2, 4, 6,32, 2, 6,32,32,
	 2, 5, 8,32, 2, 2, 3, 5, 2, 4, 6,32, 2, 6,32,32,
	 2, 5, 8,32, 2, 3, 4, 6, 2, 5, 7,32, 2, 7,32,32,
	 3, 7,32,32, 3, 1, 3, 5, 3, 3, 6,32, 3, 5,32,32,
	 3, 7,32,32, 3, 2, 4, 6, 3, 4, 7,32, 3, 6,32,32,
	 3, 7,32,32, 3, 2, 4, 6, 3, 4, 7,32, 3, 6,32,32,
	 3, 7,32,32, 3, 3, 5, 7, 3, 5, 8,32, 3, 7,32,32,
};

// ASTC trit/quint unpacking tables

static const struct astcsd_trit_quint_table_t
{
	uint32_t quint_tab[128];
	uint64_t trit_tab[256];
} astcsd_trit_quint_table = {


// ASTC quint-block unpacking table.
//  * When indexed with a 7-bit ASTC quint-block, it will
//    return a word containing 3 unpacked quints, in
//    bits 0,5,4 of each byte.
//  * When indexed with a bit-reversed 7-bit ASTC quint-block,
//    it will return a word containing 3 unpacked quints, in
//    bits 3:1 of each byte.
//static const union {uint32_t bx; float v; } ix_quint_tab[128] =
	{
		0x000000,0x040010,0x020020,0x060030,0x000401,0x040500,0x020501,0x070501,
		0x001200,0x041210,0x021220,0x061230,0x001601,0x040710,0x120701,0x070701,
		0x002008,0x042018,0x022028,0x062038,0x002409,0x040528,0x220509,0x070509,
		0x003208,0x043218,0x023228,0x063238,0x003609,0x040738,0x320709,0x070709,
		0x100004,0x140014,0x120024,0x160034,0x100405,0x140504,0x030405,0x070504,
		0x101204,0x141214,0x121224,0x161234,0x101605,0x140714,0x031605,0x070714,
		0x102808,0x182014,0x182028,0x182030,0x142809,0x180524,0x092409,0x090520,
		0x123808,0x183214,0x183228,0x183230,0x163809,0x180734,0x093609,0x090730,
		0x200002,0x240012,0x220022,0x260032,0x200403,0x240502,0x030422,0x070432,
		0x201202,0x241212,0x221222,0x261232,0x201603,0x240712,0x031622,0x071632,
		0x202800,0x242810,0x222820,0x262830,0x202805,0x240924,0x032824,0x072834,
		0x203802,0x243812,0x223822,0x263832,0x203807,0x240936,0x033826,0x073836,
		0x300006,0x340016,0x320026,0x360036,0x300407,0x340506,0x030406,0x070416,
		0x301206,0x341216,0x321226,0x361236,0x301607,0x340716,0x031606,0x071616,
		0x382808,0x382016,0x382820,0x382032,0x382809,0x380526,0x092804,0x092412,
		0x383808,0x383216,0x383822,0x383232,0x383809,0x380736,0x093806,0x093612,
	},

// ASTC trit-block unpacking table.
//  * When indexed with an 8-bit ASTC trit-block, it will return
//    a word containing 5 unpacked trits, in bits 5:4 of each byte.
//  * When indexed with a bit-reversed 8-bit ASTC trit-block, it
//    will return a word containing 5 unpacked trits, in
//    bits 2:1 of each byte.
//static const union { uint64_t bx; double v; } ix_trit_tab[256] =
	{
		0x0000000000ull,0x0200000010ull,0x0004000020ull,0x0204200000ull,
		0x0002001000ull,0x0202001010ull,0x0400001020ull,0x0402200010ull,
		0x0000022000ull,0x0200022010ull,0x0004022020ull,0x0204220020ull,
		0x0002222000ull,0x0202222010ull,0x0400222020ull,0x0402220020ull,
		0x0000100400ull,0x0200100410ull,0x0004100420ull,0x0204201400ull,
		0x0002101400ull,0x0202101410ull,0x0400101420ull,0x0402201410ull,
		0x0000122400ull,0x0200122410ull,0x0004122420ull,0x0204221420ull,
		0x2022020400ull,0x2222020410ull,0x2420020420ull,0x2422220400ull,
		0x0010000200ull,0x0210000210ull,0x0014000220ull,0x0214200200ull,
		0x0012001200ull,0x0212001210ull,0x0410001220ull,0x0412200210ull,
		0x0010022200ull,0x0210022210ull,0x0014022220ull,0x0214220220ull,
		0x0012222200ull,0x0212222210ull,0x0410222220ull,0x0412220220ull,
		0x0010140400ull,0x0210140410ull,0x0014140420ull,0x0214241400ull,
		0x0012141400ull,0x0212141410ull,0x0410141420ull,0x0412241410ull,
		0x0414102000ull,0x0414122010ull,0x0414102420ull,0x0414221420ull,
		0x2424001200ull,0x2424021210ull,0x2424041420ull,0x2424240410ull,
		0x0020000004ull,0x0220000014ull,0x0024000024ull,0x0224200004ull,
		0x0022001004ull,0x0222001014ull,0x0420001024ull,0x0422200014ull,
		0x0020022004ull,0x0220022014ull,0x0024022024ull,0x0224220024ull,
		0x0022222004ull,0x0222222014ull,0x0420222024ull,0x0422220024ull,
		0x0020100404ull,0x0220100414ull,0x0024100424ull,0x0224201404ull,
		0x0022101404ull,0x0222101414ull,0x0420101424ull,0x0422201414ull,
		0x0020122404ull,0x0220122414ull,0x0024122424ull,0x0224221424ull,
		0x2022022404ull,0x2222022414ull,0x2420022424ull,0x2422220424ull,
		0x2000000204ull,0x2200000214ull,0x2004000224ull,0x2204200204ull,
		0x2002001204ull,0x2202001214ull,0x2400001224ull,0x2402200214ull,
		0x2000022204ull,0x2200022214ull,0x2004022224ull,0x2204220224ull,
		0x2002222204ull,0x2202222214ull,0x2400222224ull,0x2402220224ull,
		0x2000140404ull,0x2200140414ull,0x2004140424ull,0x2204241404ull,
		0x2002141404ull,0x2202141414ull,0x2400141424ull,0x2402241414ull,
		0x2404102004ull,0x2404122014ull,0x2404102424ull,0x2404221424ull,
		0x2424202204ull,0x2424222214ull,0x2424242424ull,0x2424240424ull,
		0x1000000002ull,0x1200000012ull,0x1004000022ull,0x1204200002ull,
		0x1002001002ull,0x1202001012ull,0x1400001022ull,0x1402200012ull,
		0x1000022002ull,0x1200022012ull,0x1004022022ull,0x1204220022ull,
		0x1002222002ull,0x1202222012ull,0x1400222022ull,0x1402220022ull,
		0x1000100402ull,0x1200100412ull,0x1004100422ull,0x1204201402ull,
		0x1002101402ull,0x1202101412ull,0x1400101422ull,0x1402201412ull,
		0x1000122402ull,0x1200122412ull,0x1004122422ull,0x1204221422ull,
		0x2022120402ull,0x2222120412ull,0x2420120422ull,0x2422221402ull,
		0x1010000202ull,0x1210000212ull,0x1014000222ull,0x1214200202ull,
		0x1012001202ull,0x1212001212ull,0x1410001222ull,0x1412200212ull,
		0x1010022202ull,0x1210022212ull,0x1014022222ull,0x1214220222ull,
		0x1012222202ull,0x1212222212ull,0x1410222222ull,0x1412220222ull,
		0x1010140402ull,0x1210140412ull,0x1014140422ull,0x1214241402ull,
		0x1012141402ull,0x1212141412ull,0x1410141422ull,0x1412241412ull,
		0x1414102002ull,0x1414122012ull,0x1414102422ull,0x1414221422ull,
		0x2424101202ull,0x2424121212ull,0x2424141422ull,0x2424241412ull,
		0x1020040000ull,0x1220040010ull,0x1024040020ull,0x1224240000ull,
		0x1022041000ull,0x1222041010ull,0x1420041020ull,0x1422240010ull,
		0x1020042200ull,0x1220042210ull,0x1024042220ull,0x1224240220ull,
		0x1022242200ull,0x1222242210ull,0x1420242220ull,0x1422240220ull,
		0x1020140004ull,0x1220140014ull,0x1024140024ull,0x1224241004ull,
		0x1022141004ull,0x1222141014ull,0x1420141024ull,0x1422241014ull,
		0x1020142204ull,0x1220142214ull,0x1024142224ull,0x1224241224ull,
		0x2022142204ull,0x2222142214ull,0x2420142224ull,0x2422241224ull,
		0x2010040002ull,0x2210040012ull,0x2014040022ull,0x2214240002ull,
		0x2012041002ull,0x2212041012ull,0x2410041022ull,0x2412240012ull,
		0x2010042202ull,0x2210042212ull,0x2014042222ull,0x2214240222ull,
		0x2012242202ull,0x2212242212ull,0x2410242222ull,0x2412240222ull,
		0x2010140004ull,0x2210140014ull,0x2014140024ull,0x2214241004ull,
		0x2012141004ull,0x2212141014ull,0x2410141024ull,0x2412241014ull,
		0x2414142000ull,0x2414142210ull,0x2414142024ull,0x2414241224ull,
		0x2424242002ull,0x2424242212ull,0x2424242024ull,0x2424241224ull,
	},
};

// ****************************************************
//      Initialization data structures and routines
// ****************************************************

// Array of ASTC unpacked block-mode-descriptors. The
// array is indexed by the 11 lowest bits of an ASTC-block,
// plus 2048 for ASTC-3D.

typedef struct astcsd_block_mode_desc_t_
{
	uint8_t block_type;  // 0 = regular, 1 = error, 2 = LDR void-extent, 3 = HDR void-extent
	uint8_t grid_xdim;
	uint8_t grid_ydim;//:4;
	uint8_t grid_zdim;
	uint8_t dual_plane_en;
	uint8_t index_quant_level;
	uint8_t index_count;
	uint8_t index_ise_bits;

	uint8_t grid_xstride; // grid_xdim << dual_plane_en
	uint8_t grid_xpydim;  // grid_xdim + grid_ydim. Used for quick detection of 1:1 mapping for index-grid

	uint8_t padding[6];  // to ensure a 16-byte data structure.

} astcsd_block_mode_desc_t;

static astcsd_block_mode_desc_t astcsd_block_mode_descs[4096]; // 2048 descriptors for 2D, 2048 descriptors for 3D.

// Lookup table used to determine ASTC-color-ISE quantization level.
// It is indexed by [a+9*b], where
//  'a' is the number of integer pairs in the color-ISE minus 1
//  'b' is the number of bits in the color-ISE + 7.
// The result value is in the range 4 to 20, corresponding to various
// ASTC color quantization levels; for cases where the color-ISE has too few
// bits to represent the required number of integers with at least 2.6 bits
// per integer, the result value is 32.
static uint8_t astcsd_cise_quant_tab[896];

// Lookup table used to hold ASTC-partition data. This table uses
// 64 bytes for each partition-pattern for ASTC-2D, and 128 bytes
// for ASTC-3D.
static vec128_t astcsd_partition_data[36880];

// helper function to compute the number of bits in an ISE, given
// an item-count and a quantization level. A too-high quantization
// level causes it to return a very large value.
static uint32_t astcsd_compute_ise_bitcount(uint32_t items, uint32_t quant_level)
{
	// bits per item, multiplied by 15.
	static const uint8_t multipliers[21] =
		{ 15,24,30,35,39,45,50,54,60,65,69,75,80,84,90,95,99,105,110,114,120 };
	// multiply, add 15, divide by 15.
	if (quant_level > 20) return 100000;
	return ((items * multipliers[quant_level] + 14u)*8739) >> 17;
}

// routine to initialize the 'astcsd_cise_quat_tab' lookup table above.
static void astcsd_precompute_cise_quant_tab(void)
{
	uint32_t i,j;
	memset(astcsd_cise_quant_tab + 832, 32, 64); // need to cover elements from 855 to 879, but more doesn't harm.

	for (i=0;i<=8;i++)
	{
		uint32_t items = 2*(i+1);
		uint32_t curr_level = 32;
		uint32_t next_level = 4;
		uint32_t next_pos = astcsd_compute_ise_bitcount(items, next_level) + 7;

		j = 0;
		uint8_t *res = astcsd_cise_quant_tab + i;
		while(1)
		{
			while(j == next_pos)
			{
				if (j == 95) goto END_OF_LOOP; // two-level break
				curr_level = next_level++;
				next_pos = astcsd_compute_ise_bitcount(items, next_level) + 7;
				if (next_pos > 95)
					next_pos = 95;
			}
			*res = curr_level;
			res += 9;
			j++;
		}
		END_OF_LOOP: ;
	}
}

// helper function to initialize the block-mode-descriptor array described above.
static void astcsd_precompute_block_mode_descrs(void)
{
	uint32_t i;

	for (i=0;i<4096;i++)
	{
		uint32_t block_mode = i & 2047;

		uint32_t N=1, M=1, Q=1;
		uint32_t is_invalid = 0;

		uint32_t base_quant_mode = (block_mode >> 4) & 1;
		uint32_t H = (block_mode >> 9) & 1;
		uint32_t D = (block_mode >> 10) & 1;
		uint32_t A = (block_mode >> 5) & 0x3;

		if ((block_mode & 0xF) == 0) is_invalid = 1;

		uint32_t bm2 = (block_mode & 3) == 0 ? (block_mode >> 2) : block_mode;
		base_quant_mode |= (bm2 & 3) << 1;
		uint32_t B = (bm2 >> 7) & 3;

		if (i < 2048)
		{
			// ASTC-2D block modes
			if ((block_mode & 3) != 0)
			{
				switch ((block_mode >> 2) & 3)
				{
					case 0:  N = B + 4; M = A + 2; break;
					case 1:  N = B + 8; M = A + 2; break;
					case 2:  N = A + 2; M = B + 8; break;
					case 3:
						B &= 1;
						if (block_mode & 0x100) { N = B + 2; M = A + 2; }
						else                   { N = A + 2; M = B + 6; }
						break;
				}
			}
			else
			{
				switch ((block_mode >> 7) & 3)
				{
					case 0:  N = 12;    M = A + 2; break;
					case 1:  N = A + 2; M = 12;    break;
					case 2:  N = A + 6; M = B + 6; D = 0; H = 0; break;
					case 3:
						switch ((block_mode >> 5) & 3)
							{
							case 0:  N = 6;  M = 10; break;
							case 1:  N = 10; M = 6;  break;
							case 2:
							case 3:  is_invalid = 1; break;  // includes void-extent
						}
						break;
				}
			}
		}
		else
		{
			// ASTC-3D block modes
			if ((block_mode & 3) != 0)
			{
				uint32_t C = (block_mode >> 2) & 0x3;
				N = A + 2;
				M = B + 2;
				Q = C + 2;
			}
			else
			{
				if (((block_mode >> 7) & 3) != 3)  { D = 0; H = 0; }
				switch ((block_mode >> 7) & 3)
				{
					case 0:   N = 6;      M = B + 2;  Q = A + 2; break;
					case 1:   N = A + 2;  M = 6;      Q = B + 2; break;
					case 2:   N = A + 2;  M = B + 2;  Q = 6;     break;
					case 3:   N = 2;      M = 2;      Q = 2;
						switch ((block_mode >> 5) & 3)
						{
							case 0: N = 6; break;
							case 1: M = 6; break;
							case 2: Q = 6; break;
							case 3: is_invalid = 1;   // includes void-extent
						}
					break;
				}
			}
		}

		uint32_t weight_count = N * M * Q* (D + 1);
		uint32_t quant_level = (base_quant_mode - 2) + 6 * H;

		if (weight_count > 64)
		{
			weight_count = 64;
			is_invalid = 1;
		}

		uint32_t weightbits = astcsd_compute_ise_bitcount(weight_count, quant_level);

		if (weightbits < 24)
		{
			weightbits = 24;
			is_invalid = 1;
		}

		if (weightbits > 96)
		{
			weightbits = 96;
			is_invalid = 1;
		}

		astcsd_block_mode_desc_t *bmd = &(astcsd_block_mode_descs[i]);
		bmd->grid_xdim = N;
		bmd->grid_ydim = M;
		bmd->grid_zdim = Q;
		bmd->dual_plane_en = D;
		bmd->grid_xstride = N << D;
		// since grid_xpydim is used only for quick detection of 1:1 index-grid mapping
		// and the index-grid 1:1 mapping fastpath only works for index grids with
		// 8 indexes or less, grid_xpydim is set to a known-unmatchable value
		// when then grid is wider than 8.
		bmd->grid_xpydim = N > 8 ? 99 : N+M;
		bmd->index_quant_level = quant_level;
		bmd->index_count = weight_count;
		bmd->index_ise_bits = weightbits;

		bmd->block_type = is_invalid;
	}

	astcsd_block_mode_descs[1532].block_type = 2;  // 2D LDR void-extent
	astcsd_block_mode_descs[2044].block_type = 3;  // 2D HDR void-extent

	astcsd_block_mode_descs[ 508+2048].block_type = 2; // 3D LDR void-extent
	astcsd_block_mode_descs[1532+2048].block_type = 2; // 3D LDR void-extent
	astcsd_block_mode_descs[1020+2048].block_type = 3; // 3D HDR void-extent
	astcsd_block_mode_descs[2044+2048].block_type = 3; // 3D HDR void-extent
}

/*****************************************************
  Data layout of precomputed packed-partition-data:

  The main block for most of the layout is a block of 8x2 texels stored in 32 bits.
	  (x=0,y=0) -> bits 1:0
	  (x=1,y=0) -> bits 9:8
	  (x=2,y=0) -> bits 17:16
	  (x=3,y=0) -> bits 25:24
	  (x=4..7,y=0) -> bits 3:2, 11:10, 19:18, 27:26
	  (x=0,y=1) xor (x=0,y=0) -> bits 5:4
	  (x=1,y=1) xor (x=1,y=1) -> bits 13:12  and so on
	 basically, each texel starts at an offset of  ((x & 3)*8) + (((x>>2) & 1) * 2) + (y*4); also, y=1 texels are XORed with y0-texels

  There is an alternate block layout that stores 4x4 texels in 32 bits:
	  (x=0, y=0) -> bits 1:0
	  (x=1, y=0) -> bits 9:8
	  (x=2, y=0) -> bits 17:16
	  (x=3, y=0) -> bits 25:14
	  (x=0, y=1) -> bits 3:2
	  (x=0, y=2) xor (x=0, y=0) -> bits 5:4
	  (x=0, y=3) xor (x=0, y=1) -> bits 7:6
	 basically, each texel starts at an offset of  (x*8) + (y*2); y=2/y=3 texels are XORed with y=0/y=1 texels, respectively.

  The data for a partitioning is then laid out as:
	* 2D:
		[0]:one 128-bit block for small partitions
		 * bits 31:0  : 8x2-block, start at x=0, y=0
		 * bits 63:32 : 8x2-block, start at x=0, y=2
		 * bits 95:64 : 8x2-block, start at x=0, y=4
		[1]: one 128-bit block for the upper-left 8x8 block of large partitions:
		 * bits 31:0  : 8x2-block, start at x=0, y=0
		 * bits 63:32 : 8x2-block, start at x=0, y=2
		 * bits 95:64 : 8x2-block, start at x=0, y=4
		 * bits 127:96 : 8x2-block, start at x=0, y=6
		[2] one 128-bit block for the lower-left 8x4 block:
		 * bits 31:0  : 8x2-block, start at x=0, y=8
		 * bits 63:32 : 8x2-block, start at x=0, y=10
		[3] one 128-bit block for the right 12x8 block:
		 * bits 31:0   : 4x4-block: start at x=8, y=0
		 * bits 63:32  : 4x4-block: start at x=8, y=4
		 * bits 95:64  : 4x4-block: start at x=8, y=8

	* 3D: first, flat array of 3 blocks for small partitions
		  then,  flat array of 6 blocks for large partitions
		   * each block:
			 * bits 31:0  : 8x2-block, start at x=0, y=0
			 * bits 63:32 : 8x2-block, start at x=0, y=2
			 * bits 95:64 : 8x2-block, start at x=0, y=4
		  Each block is 96 bits; they are stored densely-packed;
		  this allows both arrays to fit into 128 bytes.

*********************************************************************/

static inline void astcsd_precompute_one_packed_partition(
	vec128_t packed2d[4],  // packed representations of 2D partitions
	vec128_t packed3d[8],  // packed representations of 3D partitions
	int seed,
	int partitioncount
) {
	// hash52
	uint32_t v = seed + ((partitioncount-1) << 10);
	v *= 0xEEDE0891;
	v ^= v >> 5;
	v += v << 16;
	v ^= v >> 7;
	v ^= v >> 3;
	v ^= v << 6;
	v ^= v >> 17;


	uint8_t seed1  = (uint8_t) (v & 0xF);
	uint8_t seed2  = (uint8_t) ((v >> 4) & 0xF);
	uint8_t seed3  = (uint8_t) ((v >> 8) & 0xF);
	uint8_t seed4  = (uint8_t) ((v >> 12) & 0xF);
	uint8_t seed5  = (uint8_t) ((v >> 16) & 0xF);
	uint8_t seed6  = (uint8_t) ((v >> 20) & 0xF);
	uint8_t seed7  = (uint8_t) ((v >> 24) & 0xF);
	uint8_t seed8  = (uint8_t) ((v >> 28) & 0xF);
	uint8_t seed9  = (uint8_t) ((v >> 18) & 0xF);
	uint8_t seed10 = (uint8_t) ((v >> 22) & 0xF);
	uint8_t seed11 = (uint8_t) ((v >> 26) & 0xF);
	uint8_t seed12 = (uint8_t) (((v >> 30) | (v << 2)) & 0xF);

	static const uint8_t sh1_arr[16] = { 5,5,5,4, 5,5,5,4, 6,5,6,4, 5,5,5,4 };
	static const uint8_t sh2_arr[16] = { 5,5,4,5, 5,5,4,5, 5,6,4,6, 5,5,4,5 };

	uint32_t sh_idx = (seed & 3) + 4*(partitioncount-1);
	uint32_t sh1 = sh1_arr[sh_idx];
	uint32_t sh2 = sh2_arr[sh_idx];
	uint32_t sh3 = (seed & 0x10) ? sh1 : sh2;


	// squaring all the seeds in order to bias their distribution
	// towards lower values.
	seed1 = (seed1*seed1) >> sh1;
	seed2 = (seed2*seed2) >> sh2;
	seed3 = (seed3*seed3) >> sh1;
	seed4 = (seed4*seed4) >> sh2;
	seed5 = (seed5*seed5) >> sh1;
	seed6 = (seed6*seed6) >> sh2;
	seed7 = (seed7*seed7) >> sh1;
	seed8 = (seed8*seed8) >> sh2;
	seed9 = (seed9*seed9) >> sh3;
	seed10 = (seed10*seed10) >> sh3;
	seed11 = (seed11*seed11) >> sh3;
	seed12 = (seed12*seed12) >> sh3;

	seed2 <<= 2;
	seed4 <<= 2;
	seed6 <<= 2;
	seed8 <<= 2;

	seed9 <<= 2;
	seed10 <<= 2;
	seed11 <<= 2;
	seed12 <<= 2;

	// parameter adjustment for ASTC-3D pattern generation; the point of this
	// subtraction is to be able to use the seed-derived values as addends from
	// one plane to the next, reducing the number of temporaries needed.
	seed11 -= 6*seed2;
	seed12 -= 6*seed4;
	seed9  -= 6*seed6;
	seed10 -= 6*seed8;

	uint8_t b_mask = 0xFC;
	uint8_t c_mask = 0xFC;
	uint8_t d_mask = 0xFC;

	if (partitioncount <= 1) { b_mask = 0; seed3 = 0; seed4 = 0; seed12 = 0; }
	if (partitioncount <= 2) { c_mask = 0; seed5 = 0; seed6 = 0; seed9  = 0; }
	if (partitioncount <= 3) { d_mask = 0; seed7 = 0; seed8 = 0; seed10 = 0; }

	uint8_t vx[4];
	vx[0] = (v >> 12) & 0xFC;
	vx[1] = (v >> 8) & b_mask;
	vx[2] = (v >> 4) & c_mask;
	vx[3] = (v >> 0) & d_mask;

	static const union { uint8_t s[16]; vec128_t v; } vec_1s       = {{ 1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1 }};
	static const union { uint8_t s[16]; vec128_t v; } part_2d_mask = {{ 0xff,0xff,0xff,0xff, 0xff,0xff,0xff,0xff, 0xff,0xff,0xff,0xff, 0,0,0,0 }};
	static const union { uint8_t s[16]; vec128_t v; } part_3d_mask = {{ 0xff,0xff,0xff,0xff, 0xff,0xff,0, 0, 0xff,0xff,0xff,0xff, 0xff,0xff, 0,0 }};
	static const union { uint8_t s[16]; vec128_t v; } high_mask    = {{ 0, 0, 0, 0, 0, 0, 0, 0, 0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff }};

	static const union { int8_t s[16]; vec128_t v; } decim_pat_2d = {{ 0,2,4,6, 8,10,12,14, -1,-1,-1,-1, -1,-1,-1,-1 }};
	static const union { int8_t s[16]; vec128_t v; } decim_pat_3d = {{ 0,2,4,-1, -1,-1,-1,-1, -1,-1,-1,-1, -1,-1,-1,-1 }};

	uint32_t y,z;

	// **************************
	//   Generate 2D partitions
	// **************************

	// addend-table: adtbl[i].s[j] is equal to (i*j*4 + 128) mod 256.
	// the +128 is needed because SSE's _mm_cmpgt_epi8 only supports signed compare.
	static const union { uint8_t s[16]; vec128_t v; } adtbl[32] = {
		{{128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,128,}},
		{{128,132,136,140,144,148,152,156,160,164,168,172,176,180,184,188,}},
		{{128,136,144,152,160,168,176,184,192,200,208,216,224,232,240,248,}},
		{{128,140,152,164,176,188,200,212,224,236,248,  4, 16, 28, 40, 52,}},
		{{128,144,160,176,192,208,224,240,  0, 16, 32, 48, 64, 80, 96,112,}},
		{{128,148,168,188,208,228,248, 12, 32, 52, 72, 92,112,132,152,172,}},
		{{128,152,176,200,224,248, 16, 40, 64, 88,112,136,160,184,208,232,}},
		{{128,156,184,212,240, 12, 40, 68, 96,124,152,180,208,236,  8, 36,}},
		{{128,160,192,224,  0, 32, 64, 96,128,160,192,224,  0, 32, 64, 96,}},
		{{128,164,200,236, 16, 52, 88,124,160,196,232, 12, 48, 84,120,156,}},
		{{128,168,208,248, 32, 72,112,152,192,232, 16, 56, 96,136,176,216,}},
		{{128,172,216,  4, 48, 92,136,180,224, 12, 56,100,144,188,232, 20,}},
		{{128,176,224, 16, 64,112,160,208,  0, 48, 96,144,192,240, 32, 80,}},
		{{128,180,232, 28, 80,132,184,236, 32, 84,136,188,240, 36, 88,140,}},
		{{128,184,240, 40, 96,152,208,  8, 64,120,176,232, 32, 88,144,200,}},
		{{128,188,248, 52,112,172,232, 36, 96,156,216, 20, 80,140,200,  4,}},
		{{128,192,  0, 64,128,192,  0, 64,128,192,  0, 64,128,192,  0, 64,}},
		{{128,196,  8, 76,144,212, 24, 92,160,228, 40,108,176,244, 56,124,}},
		{{128,200, 16, 88,160,232, 48,120,192,  8, 80,152,224, 40,112,184,}},
		{{128,204, 24,100,176,252, 72,148,224, 44,120,196, 16, 92,168,244,}},
		{{128,208, 32,112,192, 16, 96,176,  0, 80,160,240, 64,144,224, 48,}},
		{{128,212, 40,124,208, 36,120,204, 32,116,200, 28,112,196, 24,108,}},
		{{128,216, 48,136,224, 56,144,232, 64,152,240, 72,160,248, 80,168,}},
		{{128,220, 56,148,240, 76,168,  4, 96,188, 24,116,208, 44,136,228,}},
		{{128,224, 64,160,  0, 96,192, 32,128,224, 64,160,  0, 96,192, 32,}},
		{{128,228, 72,172, 16,116,216, 60,160,  4,104,204, 48,148,248, 92,}},
		{{128,232, 80,184, 32,136,240, 88,192, 40,144,248, 96,200, 48,152,}},
		{{128,236, 88,196, 48,156,  8,116,224, 76,184, 36,144,252,104,212,}},
		{{128,240, 96,208, 64,176, 32,144,  0,112,224, 80,192, 48,160, 16,}},
		{{128,244,104,220, 80,196, 56,172, 32,148,  8,124,240,100,216, 76,}},
		{{128,248,112,232, 96,216, 80,200, 64,184, 48,168, 32,152, 16,136,}},
		{{128,252,120,244,112,236,104,228, 96,220, 88,212, 80,204, 72,196,}},
	};

	vec128_t a2v = v128_load_8rep(vx + 0);
	vec128_t b2v = v128_load_8rep(vx + 1);
	vec128_t c2v = v128_load_8rep(vx + 2);
	vec128_t d2v = v128_load_8rep(vx + 3);

	vec128_t a2a = v128_load_8rep(&seed2);
	vec128_t b2a = v128_load_8rep(&seed4);
	vec128_t c2a = v128_load_8rep(&seed6);
	vec128_t d2a = v128_load_8rep(&seed8);

	a2v = v128_add_8(a2v, adtbl[seed1].v);
	b2v = v128_add_8(b2v, adtbl[seed3].v);
	c2v = v128_add_8(c2v, adtbl[seed5].v);
	d2v = v128_add_8(d2v, adtbl[seed7].v);

	// stash away some data items that will be useful for 3d-partitioning later.
	vec128_t a3v = v128_rep_low_half_of(a2v);
	vec128_t b3v = v128_rep_low_half_of(b2v);
	vec128_t c3v = v128_rep_low_half_of(c2v);
	vec128_t d3v = v128_rep_low_half_of(d2v);

	union { uint8_t s[192]; vec128_t v[12]; } res2d_large;
	union { uint8_t s[96];  vec128_t v[6]; } res2d_small;

	// ------------------------------------------------------------
	//   main partition-pattern generation loop for 2d partitions
	// ------------------------------------------------------------
	for (y=0;y<12;y++)
	{
		vec128_t v0 = v128_or(v128_or(v128_cmpgt_s8(b2v, a2v), v128_cmpgt_s8(c2v,a2v)), v128_cmpgt_s8(d2v,a2v));
		vec128_t v1 = v128_or(v128_cmpgt_s8(d2v,b2v), v128_cmpgt_s8(c2v,b2v));
		vec128_t v2 = v128_cmpgt_s8(d2v,c2v);
		vec128_t v3 = v128_and(v1, v2);
		vec128_t partition = v128_sub_8(v128_sub_8(vec_1s.v , v1), v3);
		res2d_large.v[y] = v128_and(partition, v128_and(v0, part_2d_mask.v));

		a2v = v128_add_8(a2v, a2a);
		b2v = v128_add_8(b2v, b2a);
		c2v = v128_add_8(c2v, c2a);
		d2v = v128_add_8(d2v, d2a);
	}

	// generate small partition from the large partition by 2x decimation on each axis
	for (y=0;y<6;y++)
	{
		res2d_small.v[y] = v128_permute_8x16(res2d_large.v[2*y], decim_pat_2d.v);
	}

	// ------------------------------
	//   pack the 2D partition data
	//
	// Result after backing: four 128-bit words.
	//  * Result word 0 contains 8x8 data, for the small partitioning
	//  * Result word 1 contains 8x8 data, for the upper-left part of the block
	//  * Result word 2 contains 8x4 data, for the lower-left part of the block
	//  * Result word 3 contains 16x4 data, for the right part of the block
	// ------------------------------
	uint32_t *dstp = (uint32_t *)packed2d;
	const uint32_t *srcpl = (uint32_t *)res2d_large.s;
	const uint32_t *srcps = (uint32_t *)res2d_small.s;

	// for each loop iteration: process 2 rows of small-partitioning
	// and 4 rows of large-partitioning.
	for (y=0;y<3;y++)
	{
		uint32_t vl0, vl1;
		vl0 = srcps[0] + (srcps[1] << 2);  vl1 = srcps[ 4] + (srcps[ 5] << 2);   dstp[y]       = ((vl0^vl1) << 4) + vl0;  // two rows of small-partition
		vl0 = srcpl[0] + (srcpl[1] << 2);  vl1 = srcpl[ 4] + (srcpl[ 5] << 2);   dstp[2*y + 4] = ((vl0^vl1) << 4) + vl0;  // rows 0,1 (modulo 4; large-partition)
		vl0 = srcpl[8] + (srcpl[9] << 2);  vl1 = srcpl[12] + (srcpl[13] << 2);   dstp[2*y + 5] = ((vl0^vl1) << 4) + vl0;  // rows 2,3 (modulo 4; large-partition)
		vl0 = srcpl[2] + (srcpl[6] << 2);  vl1 = srcpl[10] + (srcpl[14] << 2);   dstp[y+12]    = ((vl0^vl1) << 4) + vl0;  // rows 0,1,2,3 (modulo 4)  (right column of large-partition)
		srcpl += 16;
		srcps += 8;
	}

	dstp[3]  = dstp[10] = dstp[11] = dstp[15] = 0;

	// **************************
	//   Generate 3D partitions
	// **************************
	vec128_t a3b = v128_load_8rep(&seed11);
	vec128_t b3b = v128_load_8rep(&seed12);
	vec128_t c3b = v128_load_8rep(&seed9 );
	vec128_t d3b = v128_load_8rep(&seed10);

	a3v = v128_add_8(a3v, v128_and(a2a, high_mask.v));
	b3v = v128_add_8(b3v, v128_and(b2a, high_mask.v));
	c3v = v128_add_8(c3v, v128_and(c2a, high_mask.v));
	d3v = v128_add_8(d3v, v128_and(d2a, high_mask.v));

	a2a = v128_add_8(a2a, a2a);
	b2a = v128_add_8(b2a, b2a);
	c2a = v128_add_8(c2a, c2a);
	d2a = v128_add_8(d2a, d2a);

	// buffer to hold 3d partitionings:
	//  v[0..8] : 3 layers for small partitioning
	//  v[9..27] : 6 layer for large partitioning
	union { uint32_t d[108]; uint64_t q[54]; vec128_t v[27]; } res3d;

	// ------------------------------------------------------------
	//   main partition-pattern generation loop for 3d partitions
	// ------------------------------------------------------------
	for (z=0;z<6;z++)
	{
		// innerloop processes 2 rows of 1 layer per iteration.
		for (y=0;y<3;y++)
		{
			vec128_t v0 = v128_or(v128_or(v128_cmpgt_s8(b3v, a3v), v128_cmpgt_s8(c3v,a3v)), v128_cmpgt_s8(d3v,a3v));
			vec128_t v1 = v128_or(v128_cmpgt_s8(d3v,b3v), v128_cmpgt_s8(c3v,b3v));
			vec128_t v2 = v128_cmpgt_s8(d3v,c3v);
			vec128_t v3 = v128_and(v1, v2);
			vec128_t partition = v128_sub_8(v128_sub_8(vec_1s.v , v1), v3);
			res3d.v[y + 3*z + 9] = v128_and(partition, v128_and(v0, part_3d_mask.v));

			a3v = v128_add_8(a3v, a2a);
			b3v = v128_add_8(b3v, b2a);
			c3v = v128_add_8(c3v, c2a);
			d3v = v128_add_8(d3v, d2a);
		}

		a3v = v128_add_8(a3v, a3b);
		b3v = v128_add_8(b3v, b3b);
		c3v = v128_add_8(c3v, c3b);
		d3v = v128_add_8(d3v, d3b);
	}

	// generate small partition from the large partition by 2x decimation on each axis
	for (z=0;z<3;z++)
	{
		for (y=0;y<3;y++)
		{
			res3d.q[6*z+y] = u64_from_low_half_of(v128_permute_8x16(res3d.v[6*z+y + 9], decim_pat_3d.v));
			res3d.q[6*z+y+3] = 0;
		}
	}

	// ------------------------------
	//   Pack the 3D partition-data
	// ------------------------------
	dstp = (uint32_t *)packed3d;
	const uint32_t *srcp = res3d.d;
	for (z=0;z<9;z++)
	{
		for (y=0;y<3;y++)
		{
			uint32_t vl0 = srcp[0] + (srcp[1] << 2);
			uint32_t vl1 = srcp[2] + (srcp[3] << 2);
			dstp[y] = ((vl0 ^ vl1) << 4) + vl0;

			srcp += 4;
		}
		dstp += 3;
	}

	dstp[0] = dstp[1] = dstp[2] = dstp[3] = dstp[4] = 0; // pad end of block
}

static void astcsd_precompute_packed_partitions(
	vec128_t dst[36880] // should be at least 128-byte aligned.
) {
	vec128_t *dst_2d = dst;
	vec128_t *dst_3d = dst + 12296;

	memset(dst_2d,     0,  64); // write one 2d partitioning as all-0s.
	memset(dst_3d - 4, 0, 192); // write one 3d partitioning as all-0s; also, write the padding area between 2d and 3d partitionings as all-0s.

	uint32_t seed, partcnt;
	for (partcnt=2; partcnt <=4; partcnt++)
	{
		for (seed=0; seed<1024; seed++)
		{
			dst_2d += 4;
			dst_3d += 8;
			astcsd_precompute_one_packed_partition(dst_2d, dst_3d, seed, partcnt);
		}
	}
}

static int astcsd_precompute_completed = 0;

void astc_simd_decode_precompute(void)
{
	astcsd_precompute_bit_extract_and_deposit_tables();
	astcsd_precompute_cise_quant_tab();
	astcsd_precompute_block_mode_descrs();
	astcsd_precompute_packed_partitions(astcsd_partition_data);
	astcsd_precompute_completed = 1;
}

// *****************************************************************************
// *****************************************************************************
// *****                                                                   *****
// *****              Subroutines for per-block ASTC processing            *****
// *****                                                                   *****
// *****************************************************************************
// *****************************************************************************

static void astcsd_index_ise_unpack(
	const uint8_t* blk,

	// ix_quant_level:
	//   0: 0..1 (1 bit)
	//   1: 0..2 (1 trit)
	//   2: 0..3 (2 bits)
	//   3: 0..4 (1 quint)
	//   4: 0..5 (1 bit + 1 trit)
	//   5: 0..7 (3 bits)
	//   6: 0..9 (1 bit + 1 quint)
	//   7: 0..11 (2 bits + 1 trit)
	//   8: 0..15 (4 bits)
	//   9: 0..19 (2 bits + 1 quint)
	//  10: 0..23 (3 bits + 1 trit)
	//  11: 0..31 (5 bits)
	uint32_t ix_quant_level,

	// 4 to 64
	uint32_t index_count,

	vec256_t *vpx
) {
	// -------------------------------------------------------
	//   AVX2 constant-value vectors that are used unindexed
	// -------------------------------------------------------
	static const union { uint8_t b[32]; vec256_t v; } val32_vec = {{ 32,32,32,32, 32,32,32,32, 32,32,32,32, 32,32,32,32, 32,32,32,32, 32,32,32,32, 32,32,32,32, 32,32,32,32 }};
	static const union { uint8_t b[32]; vec256_t v; } val63_vec = {{ 63,63,63,63, 63,63,63,63, 63,63,63,63, 63,63,63,63, 63,63,63,63, 63,63,63,63, 63,63,63,63, 63,63,63,63 }};

	// For trit and quint formats, we will perform expansion of the 'bc'-bits of the index
	// into an integer that is added to the rest of the index. This is done as follows:
	//  * format 7: the 'b'-bit is mapped into bit 2 of each lane : indexes 0 and 4 are mapped to values 0 and 69.
	//  * format 9: the 'b'-bit is mapped into bit 0 of each lane : indexes 0 and 1 are mapped to values 0 and 66.
	//  * format 10: the 'bc'-bits are mapped into bits 1:0 : indexes 0,1,2,3 are mapped to values 0,66,33,99.
	// For all other index-formats, bc-bits do not exist; for these formats, index is set to 0, resulting in value 0.
	// Index values 5-15 are not used.
	static const union { uint8_t b[32]; vec256_t v; } ix_tq_addends_tab = {{ 0,66,33,99, 69,_,_,_, _,_,_,_, _,_,_,_,   0,66,33,99, 69, _,_,_, _,_,_,_, _,_,_,_, }};

	// ------------------------------------------------------
	//   Indexable lookup tables. These are folded into one
	//   big data structure in order to limit the number
	//   of LEA instructions needed.
	// ------------------------------------------------------
	static const struct {
		uint64_t deposit_masks[12];
		uint64_t tq_masks_high[12];
		uint64_t tq_masks_low[12];
		uint64_t b_masks_low[12];
		uint64_t tq_high_shifts[12];
		union { uint8_t s[64]; vec256_t v[2]; } b_rrottab;
		union { uint8_t s[16]; vec128_t v; } xlat_tabs[12];
		uint32_t tq_shift_and_case_arr[112]; } ixd = {


		// Bit-deposit masks, applied to bits in the ISE.
		// These are used as arguments to PDEP, to help expand these bits
		// into byte-lanes.
		//static const uint64_t ix_deposit_masks[12] =
		{
			0x0101010101010101ull, // format 0: 1-bit
			0,                     // format 1: trit
			0x0303030303030303ull, // format 2: 2-bit
			0,                     // format 3: quint
			0x4040404040404040ull, // format 4: trit + 1 bit
			0x0707070707070707ull, // format 5: 3-bit
			0x4040404040404040ull, // format 6: quint + 1 bit
			0x4444444444444444ull, // format 7: trit + 2 bits; see note on 'ix_tw_addends_tab'
			0x0f0f0f0f0f0f0f0full, // format 8: 4-bit
			0x4141414141414141ull, // format 9: quint + 2 bits; see note on 'ix_tw_addends_tab'
			0x4343434343434343ull, // format 10: trit + 3 bits; see note on 'ix_tw_addends_tab'
			0x2f2f2f2f2f2f2f2full, // format 11: 5-bit
		},



		// bit-extraction masks:
		//  for trit/quint formats, extract the top 64 bits of the trit/quint-data

		//static const uint64_t ix_tq_masks_high[12] =
		{
			0,
			0xffffffffffffffffull,  // format 1: 1 trit: all-1s
			0,
			0xffffffffffffffffull,  // format 3: 1 quint: all-1s
			0x6d6b6b5b5adad6d6ull,  // format 4: 1 bit + 1 trit: bit-repeating pattern 0110110101101
			0,
			0x76ddb76ddb76ddb7ull,  // format 6: 1 bit + 1 quint: bit-repeating pattern 0111011011
			0x33264cc9933264ccull,  // format 7: 2 bits + 1 trit: bit-repeating pattern 001100110010011001
			0,
			0x3999ccce66733399ull,  // format 9: 2 bits + 1 quint: bit-repeating pattern 0011100110011
			0x18c4623188c46311ull,  // format 10: 3 bits + 1 trit: bit-repeating pattern 00011000110001000110001
			0
		},

	//static const uint64_t ix_tq_masks_low[12] =
		{
			0,0,0,0,
			0xb6b5b5ad0000001full,  // format 4
			0,
			0x6ddb76d000000000ull,  // format 6
			0x9933264c003fffffull,  // format 7
			0,
			0x9ccce66700000fffull,  // format 9
			0x88c623117fffffffull,  // format 10
			0
		},

	// bit-extraction masks:
	//  for trit/quint formats, extract the bit-data; this is at most 60 bits.

	//static const uint64_t ix_b_masks_low[12] =
		{
			0,0,0,0,
			0x494a4a5201ffffffull,  // format 4
			0,
			0x92248923ffffffffull,  // format 6
			0x66ccd9b3000000ffull,  // format 7
			0,
			0x633319980003ffffull,  // format 9
			0x7739dce000000003ull,  // format 10
			0
		},

		//static const uint8_t ix_tq_high_shifts[12] =
		{ 0,0,0,0,25,0,19,36,0,30,42,0 },


		// rotation-table for use with the bitwise-ISE unpacking.
		//static const union { uint8_t s[64]; vec256_t v[2]; } ix_b_rrottab =
		{{
			0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
			0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
			0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
			0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
		}},

	// translation table for bits, trits and quints

	// For the bit-based formats with 1-4 bits, this contains the final index value.
	// For the 5-bit format, the value contains an index value that can be derived
	// from the top 4 bits; if the bottom bit 'e' is set, then 4*e must be added
	// to the value from this table to form the final index value.
	// The bit-based index values here have a 2*N+1 transform baked in.

	// For the trit/quint formats, trits and quints are provided from lookup tables
	// and translated into an addend-value. The lookup table item provides the
	// trit or quint in bits 3:1, with bit 0 not necessarily 0. 2*N+1 is not
	// baked in for these entries.

	//static const union { uint8_t s[16]; vec128_t v; } ix_xlat_tabs[12] = {
		{
			{{  1,129,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _  }}, // 1 bit
			{{  0,  0,128,128,255,255,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _  }}, // 1 trit
			{{  1, 87, 43,129,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _  }}, // 2 bits
			{{  0,  0, 64, 64,128,128,191,191,255,255,  _,  _,  _,  _,  _,  _  }}, // 1 quint
			{{  0,  0, 50, 50,100,100,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _  }}, // 1 bit + 1 trit
			{{  1, 75, 37,111, 19, 93, 55,129,  _,  _,  _,  _,  _,  _,  _,  _  }}, // 3 bits
			{{  0,  0, 28, 28, 56, 56, 84, 84,112,112,  _,  _,  _,  _,  _,  _  }}, // 1 bit + 1 quint
			{{  0,  0, 23, 23, 46, 46,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _  }}, // 2 bits + 1 trit
			{{  1, 71, 35,105, 17, 87, 51,121,  9, 79, 43,113, 25, 95, 59,129  }}, // 4-bit
			{{  0,  0, 13, 13, 26, 26, 39, 39, 52, 52,  _,  _,  _,  _,  _,  _  }}, // 2 bits + 1 quint
			{{  0,  0, 11, 11, 22, 22,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _  }}, // 3 bits + 1 trit
			{{  1, 69, 33,101, 17, 85, 49,117,  9, 77, 41,109, 25, 93, 57,125  }}, // 5 bits
		},


	// For index ISEs with trit/quint encodings and a number of indexes that
	// are not divisible by the number of trits/quints in one block, we need
	// to zero out the low-order bits of the trit/quint block; this is done
	// by taking a mask of all-1s, left-shift by the appropriate amount, then
	// ANDing it with the packed trit/quint array. The left-shifts are given by
	// the arrays below; they're indexed by index-count and the trit-vs-quint
	// distinction. (For index-ISEs with trits, the maximum number indexes we
	// can have is 60; for index-ISEs with quints, the maximum number is 42.

	// The index that is used to select shift-amount is also used to select
	// switch-case; as such, we will look up a 32-bit item that contains both
	// the switch-case and the shift-amounts.
	//   bits 4:0 : switch-case
	//   bits 15:8 : shift-amount to use for the high part of the trit/quint-mask
	//   bits 23:16: shift-amount to use for the low part of the trit/quint-mask

	//static const uint32_t ix_tq_shift_and_case_arr[112] = {
		{
			0x3f3f00,0x3f3e0e,0x3f3c0e,0x3f3b0e,0x3f390e,0x3f380e,0x3f360f,0x3f340f,0x3f330f,0x3f310f,0x3f300f,0x3f2e10,0x3f2c10,0x3f2b10,0x3f2910,0x3f2810,
			0x3f2611,0x3f2411,0x3f2311,0x3f2111,0x3f2011,0x3f1e12,0x3f1c12,0x3f1b12,0x3f1912,0x3f1812,0x3f1613,0x3f1413,0x3f1313,0x3f1113,0x3f1013,0x3f0e14,
			0x3f0c14,0x3f0b14,0x3f0914,0x3f0814,0x3f0615,0x3f0415,0x3f0315,0x3f0115,0x3f0015,0x3e0016,0x3c0016,0x3b0016,0x390016,0x380016,0x360017,0x340017,
			0x330017,0x310017,0x300017,0x2e0018,0x2c0018,0x2b0018,0x290018,0x280018,0x260019,0x240019,0x230019,0x210019,0x200019,0x000000,0x000000,0x000000,
			0x3f3f00,0x3f3d00,0x3f3b00,0x3f3900,0x3f3601,0x3f3401,0x3f3201,0x3f2f02,0x3f2d02,0x3f2b02,0x3f2803,0x3f2603,0x3f2403,0x3f2104,0x3f1f04,0x3f1d04,
			0x3f1a05,0x3f1805,0x3f1605,0x3f1306,0x3f1106,0x3f0f06,0x3f0c07,0x3f0a07,0x3f0807,0x3f0508,0x3f0308,0x3f0108,0x3e0009,0x3c0009,0x3a0009,0x37000a,
			0x35000a,0x33000a,0x30000b,0x2e000b,0x2c000b,0x29000c,0x27000c,0x25000c,0x22000d,0x20000d,0x1e000d,0x000000,0x000000,0x000000,0x000000,0x000000,
		} };

	// ------------------------
	//   End of lookup tables
	// ------------------------

	// table lookups common for both bit and trit/quint index-ISE unpacking.
	uint64_t ix_deposit_mask = ixd.deposit_masks[ix_quant_level];
	vec256_t ix_xlat_tab = v256_set_128rep(ixd.xlat_tabs[ix_quant_level].v); // translation table for multiplying trits/quints.

	uint32_t ix_quant_level_x5 = ix_quant_level * 5;

	// Clear bytes 32..79 of result buffer. These bytes may not get written
	// otherwise, in which case they would be left uninitialized - this could
	// then cause uninitialized values to be used in infill.
	vpx[1] = v256_zero();
	*(vec128_t *)(&vpx[2]) = v128_zero();

	if ((2341 >> ix_quant_level) & 1)
	{
		// **************************************************************
		//   Index ISE unpacking and unquantization for bit-only formats
		// **************************************************************
		uint32_t bits_per_index = (ix_quant_level_x5+25u) >> 4;  // 0->1, 2->2, 5->3, 8->4, 11->5
		uint32_t blk_count = (index_count+7u) >> 3;  // number of indexes divided by 8, rounded up.

		// we will use VPSHUFB to repeatedly byte-rotate the block in order to
		// ensure that the bits that serve as input to PDEP are at the bottom of the block;
		// this requires an initial rotation-mask to do the first rotate and a second rotation
		// mask for the per-iteration rotations. These can be obtained by performing two unaligned
		// accesses into an array that only contains consecutive integers.
		vec128_t init_rot = v128_load_unaligned(ixd.b_rrottab.s + 32u - (bits_per_index * blk_count));
		vec128_t rrot_vec = v128_load_unaligned(ixd.b_rrottab.s + bits_per_index);

		// load the block, and apply initial rotation
		vec128_t blkv = v128_permute_8x16(v128_load_unaligned(blk), init_rot);

		vec128_t vhi = v128_zero();
		uint64_t b0,b1 = 0;

		// main loop
		switch(blk_count & 7u)
		{
			case 0:
				b1 = u64_bit_deposit(u64_from_low_half_of(blkv), ix_deposit_mask);
				blkv = v128_permute_8x16(blkv, rrot_vec);
				ASTCSD_FALLTHROUGH;
			case 7:
				b0 = u64_bit_deposit(u64_from_low_half_of(blkv), ix_deposit_mask);
				blkv = v128_permute_8x16(blkv, rrot_vec);
				vhi = v128_cat_64(b1, b0);
				ASTCSD_FALLTHROUGH;
			case 6:
				b1 = u64_bit_deposit(u64_from_low_half_of(blkv), ix_deposit_mask);
				blkv = v128_permute_8x16(blkv, rrot_vec);
				blk += bits_per_index;
				ASTCSD_FALLTHROUGH;
			case 5:
				{
					b0 = u64_bit_deposit(u64_from_low_half_of(blkv), ix_deposit_mask);
					blkv = v128_permute_8x16(blkv, rrot_vec);
					vec128_t vlo = v128_cat_64(b1, b0);

					// Processing of top 32 indexes. This consists of just a couple of shuffles.
					// For 5-bit indexes, more arithmetic is needed, but this can only
					// affect the bottom 32 indexes, we so will not do that arithmetic here.
					vec256_t v1 = v256_cat_128(vhi, vlo);
					vec256_t v3 = v256_permute_8x16_pair(ix_xlat_tab, v1);
					vpx[1] = v256_byterev_64(v3);
				}
				ASTCSD_FALLTHROUGH;

			case 4:
				b1 = u64_bit_deposit(u64_from_low_half_of(blkv), ix_deposit_mask);
				blkv = v128_permute_8x16(blkv, rrot_vec);
				ASTCSD_FALLTHROUGH;
			case 3:
				b0 = u64_bit_deposit(u64_from_low_half_of(blkv), ix_deposit_mask);
				blkv = v128_permute_8x16(blkv, rrot_vec);
				vhi = v128_cat_64(b1, b0);
				ASTCSD_FALLTHROUGH;
			case 2:
				b1 = u64_bit_deposit(u64_from_low_half_of(blkv), ix_deposit_mask);
				blkv = v128_permute_8x16(blkv, rrot_vec);
				ASTCSD_FALLTHROUGH;
			case 1:
				{
					b0 = u64_bit_deposit(u64_from_low_half_of(blkv), ix_deposit_mask);
					vec128_t vlo = v128_cat_64(b1, b0);

					// Processing of bottom 32 indexes. This contains the same swizzles
					// as the top 32 indexes, but also a few steps to help handle
					// the bottom bit of the 5-bit indexes correctly.
					vec256_t v1 = v256_cat_128(vhi, vlo);
					vec256_t v2 = v256_and(v1, val32_vec.v);
					vec256_t v3 = v256_permutez_8x16_pair(ix_xlat_tab, v1);
					vec256_t v4 = v256_shri_u16(v2, 3);
					vpx[0] = v256_byterev_64(v256_add_8(v3, v4));
				}
		}
	}
	else
	{
		// *****************************************************************
		//   Index ISE unpacking and unquantization for trit/quint formats
		// *****************************************************************

		// look up a couple of shift-amounts and a switch-case in one go.
		uintptr_t tqia = (37376u >> ix_quant_level) & 64u; // 64 for quint-formats, 0 for anything else
		uintptr_t tq_mask_idx = index_count + tqia;
		uint32_t tq_shamts_and_case = ixd.tq_shift_and_case_arr[tq_mask_idx];

		uint64_t blk_low = *(const uint64_t *)(blk);
		uint64_t blk_high = *(const uint64_t *)(blk+8);

		uint8_t ix_tq_high_shift = ixd.tq_high_shifts[ix_quant_level];
		uint64_t ix_tq_mask_high = ixd.tq_masks_high[ix_quant_level];

		// get hold of the trit/quint bits. also, get hold of the b-bits
		astcsd_u64_pair_t blk_high_split = u64_bit_extract_pair(blk_high, ix_tq_mask_high, ~ix_tq_mask_high);
		astcsd_u64_pair_t blk_low_split  = u64_bit_extract_pair(blk_low, ixd.tq_masks_low[ix_quant_level], ixd.b_masks_low[ix_quant_level]);

		uint64_t tq_bits_high = (blk_high_split.v0 << ix_tq_high_shift) + blk_low_split.v0;
		uint64_t b_bits = (blk_high_split.v1 << (-ix_tq_high_shift & 0x3f)) + (blk_low_split.v1 << 2);

		uint32_t tql_shift = (1792u >> ix_quant_level) & 28u;  // 1->0, 3->0, 6->28, other formats don't-care.
		uint64_t tq_bits_low = blk_low << (tql_shift & 0x3f);

		// prepare B-bit unpacking
		uint32_t b_rotamt_x2 = ix_quant_level_x5 & 48u;
		uint32_t b_rotamt = b_rotamt_x2 >> 1;  // 4->8, 6->8, 7->16, 9->16, 10->24, other formats don't-care.

		// ---------------------------------------------------------------
		//   preparation items steps before the main index-decode-switch
		// ---------------------------------------------------------------
		vec128_t t;
		vec128_t acc0 = v128_zero();
		vec128_t acc1 = acc0;
		const struct astcsd_trit_quint_table_t *tqt = &astcsd_trit_quint_table;

		// perform masking of the trits and quints
		tq_bits_high &= (-1ull) << ((tq_shamts_and_case >> 8) & 0x3f);
		tq_bits_low  &= (-1ull) << ((tq_shamts_and_case >> 16) & 0x3f);

		// Force operands to be simultaneously present in registers.
		// This is observed to improve register spilling decisions and
		// thereby slightly improve performance under AVX2.
		#if defined(__GNUC__) && defined(__AVX2__)
			asm("" : "+x"(acc1), "+x"(acc0), "+r"(tqt) :);
		#endif

		// ---------------------------------------------------------------------------
		// main index-decode-switch: this is an unrolled for-loop implemented
		//    as a switch-statement where most entries have deliberate fallthrough;
		//    this is done to fetch and unpack the trit/quint blocks so that each
		//    trit/quint ends up in 1 byte and so that we don't fetch
		//    any unneeded blocks.
		//
		// interspersed in this loop are a series of processing steps that
		// are used to postprocess the indexes after a full 32-byte vector's worth
		// of indexes has been gathered.
		// ---------------------------------------------------------------------------
		switch(tq_shamts_and_case & 0x1f)
		{
			case 13: // decode quint-block #13
				t = v128_load_32zx(tqt->quint_tab + ((tq_bits_low >> 30) & 0x7f));
				acc0 = v128_shlb(t, 7);
				ASTCSD_FALLTHROUGH;
			case 12: // decode quint-block #12
				t = v128_load_32zx(tqt->quint_tab + ((tq_bits_low >> 37) & 0x7f));
				acc0 = v128_or(acc0, v128_shlb(t, 4));
				ASTCSD_FALLTHROUGH;
			case 11: // decode quint-block #11
				t = v128_load_32zx(tqt->quint_tab + ((tq_bits_low >> 44) & 0x7f));
				acc0 = v128_or(acc0, v128_shlb(t, 1));
				ASTCSD_FALLTHROUGH;
			case 10: // decode quint-block #10
				{
					t = v128_load_32zx(tqt->quint_tab + ((tq_bits_low >> 51) & 0x7f));
					acc0 = v128_or(acc0, v128_shrb(t, 2));

					// process and emit a block of 16 indexes.
					// The only quint format that supports more than 32 indexes is format #3, in which we
					// have only a quint and no bits; as such, bit-processing can be skipped.
					vec128_t vl = acc0;
					vl = v128_permutez_8x16(v128_from_low_half_of(ix_xlat_tab), vl);
					vl = v128_shr2_u8(vl);
					vl = v128_sub_8(vl, v128_cmpgt_s8(vl, v128_from_low_half_of(val32_vec.v)));

					vl = v128_mul2add1_8(vl);
					*(vec128_t *)(&vpx[1]) = vl;

					acc1 = v128_shlb(t, 14);
				}
				ASTCSD_FALLTHROUGH;
			case 9: // decode quint-block #9
				t = v128_load_32zx(tqt->quint_tab + (((tq_bits_low >> 58) | (tq_bits_high << 6)) & 0x7f));
				acc1 = v128_or(acc1, v128_shlb(t, 11));
				ASTCSD_FALLTHROUGH;
			case 8: // decode quint-block #8
				t = v128_load_32zx(tqt->quint_tab + ((tq_bits_high >> 1) & 0x7f));
				acc1 = v128_or(acc1, v128_shlb(t, 8));
				ASTCSD_FALLTHROUGH;
			case 7: // decode quint-block #7
				t = v128_load_32zx(tqt->quint_tab + ((tq_bits_high >> 8) & 0x7f));
				acc1 = v128_or(acc1, v128_shlb(t, 5));
				ASTCSD_FALLTHROUGH;
			case 6: // decode quint-block #6
				t = v128_load_32zx(tqt->quint_tab + ((tq_bits_high >> 15) & 0x7f));
				acc1 = v128_or(acc1, v128_shlb(t, 2));
				ASTCSD_FALLTHROUGH;
			case 5: // decode quint-block #5
				t = v128_load_32zx(tqt->quint_tab + ((tq_bits_high >> 22) & 0x7f));
				acc1 = v128_or(acc1, v128_shrb(t, 1));
				acc0 = v128_shlb(t, 15);
				ASTCSD_FALLTHROUGH;
			case 4: // decode quint-block #4
				t = v128_load_32zx(tqt->quint_tab + ((tq_bits_high >> 29) & 0x7f));
				acc0 = v128_or(acc0, v128_shlb(t, 12));
				ASTCSD_FALLTHROUGH;
			case 3: // decode quint-block #3
				t = v128_load_32zx(tqt->quint_tab + ((tq_bits_high >> 36) & 0x7f));
				acc0 = v128_or(acc0, v128_shlb(t, 9));
				ASTCSD_FALLTHROUGH;
			case 2: // decode quint-block #2
				t = v128_load_32zx(tqt->quint_tab + ((tq_bits_high >> 43) & 0x7f));
				acc0 = v128_or(acc0, v128_shlb(t, 6));
				ASTCSD_FALLTHROUGH;
			case 1: // decode quint-block #1
				t = v128_load_32zx(tqt->quint_tab + ((tq_bits_high >> 50) & 0x7f));
				acc0 = v128_or(acc0, v128_shlb(t, 3));
				ASTCSD_FALLTHROUGH;
			case 0: // decode quint-block #0
				t = v128_load_32zx(tqt->quint_tab + ((tq_bits_high >> 57) & 0x7f));
				break;
			case 26:
			case 27:
			case 28:
			case 29:
			case 30:
			case 31:
			case 25: // decode trit-block #11
				t = v128_load_64zx(tqt->trit_tab + ((tq_bits_low >> 32) & 0xff));
				acc1 = v128_shlb(t, 7);
				ASTCSD_FALLTHROUGH;
			case 24: // decode trit-block #10
				t = v128_load_64zx(tqt->trit_tab + ((tq_bits_low >> 40) & 0xff));
				acc1 = v128_or(acc1, v128_shlb(t, 2));
				ASTCSD_FALLTHROUGH;
			case 23: // decode trit-block #9
				t = v128_load_64zx(tqt->trit_tab + ((tq_bits_low >> 48) & 0xff));
				acc1 = v128_or(acc1, v128_shrb(t, 3));
				acc0 = v128_shlb(t, 13);
				ASTCSD_FALLTHROUGH;
			case 22: // decode trit-block #8
				t = v128_load_64zx(tqt->trit_tab + ((tq_bits_low >> 56) & 0xff));
				acc0 = v128_or(acc0, v128_shlb(t, 8));
				ASTCSD_FALLTHROUGH;
			case 21: // decode trit-block #7
				t = v128_load_64zx(tqt->trit_tab + ((tq_bits_high >> 0) & 0xff));
				acc0 = v128_or(acc0, v128_shlb(t, 3));
				ASTCSD_FALLTHROUGH;
			case 20: // decode trit-block #6
				{
					t = v128_load_64zx(tqt->trit_tab + ((tq_bits_high >> 8) & 0xff));
					acc0 = v128_or(acc0, v128_shrb(t, 2));
					vec256_t vl = v256_cat_128(acc1, acc0);

					// Process and emit a block of 32 indexes.
					// For format #4, we can have up to 4 indexes that have a bit in addition to just a trit;
					// as such, for the remaining indexes, we do trit-processing only.
					vec128_t bx0 = v128_set_64zx(u64_bit_deposit(b_bits >> (b_rotamt + b_rotamt_x2), ix_deposit_mask));

					vl = v256_permutez_8x16_pair(ix_xlat_tab, vl);
					vec256_t vbx = v256_cat_128(v128_zero(), v128_byterev_64(bx0));

					vl = v256_add_8(vl, v256_permute_8x16_pair(ix_tq_addends_tab.v, vbx));
					vl = v256_xor(vl, v256_cmpgt_s8(vbx, val63_vec.v));
					vl = v256_shr2_u8(vl);
					vl = v256_sub_8(vl, v256_cmpgt_s8(vl, val32_vec.v));  // for each index, add 1 if index value is greater than 32.

					vl = v256_mul2add1_8(vl);
					vpx[1] = vl;

					acc1 = v128_shlb(t, 14);
				}
				ASTCSD_FALLTHROUGH;
			case 19: // decode trit-block #5
				t = v128_load_64zx(tqt->trit_tab + ((tq_bits_high >> 16) & 0xff));
				acc1 = v128_or(acc1, v128_shlb(t, 9));
				ASTCSD_FALLTHROUGH;
			case 18: // decode trit-block #4
				t = v128_load_64zx(tqt->trit_tab + ((tq_bits_high >> 24) & 0xff));
				acc1 = v128_or(acc1, v128_shlb(t, 4));
				ASTCSD_FALLTHROUGH;
			case 17: // decode trit-block #3
				t = v128_load_64zx(tqt->trit_tab + ((tq_bits_high >> 32) & 0xff));
				acc1 = v128_or(acc1, v128_shrb(t, 1));
				acc0 = v128_shlb(t, 15);
				ASTCSD_FALLTHROUGH;
			case 16: // decode trit-block #2
				t = v128_load_64zx(tqt->trit_tab + ((tq_bits_high >> 40) & 0xff));
				acc0 = v128_or(acc0, v128_shlb(t, 10));
				ASTCSD_FALLTHROUGH;
			case 15: // decode trit-block #1
				t = v128_load_64zx(tqt->trit_tab + ((tq_bits_high >> 48) & 0xff));
				acc0 = v128_or(acc0, v128_shlb(t, 5));
				ASTCSD_FALLTHROUGH;
			case 14: // decode trit-block #0
				t = v128_load_64zx(tqt->trit_tab + ((tq_bits_high >> 56) & 0xff));
				break;
		}

		// shared code between the end of quint-processing and end of trit-processing.
		acc0 = v128_or(acc0, t);
		vec256_t vl = v256_cat_128(acc1, acc0);

		// unpack the B-bits for the low 32-bit indexes
		uint64_t bbl = u64_left_rotate(b_bits, b_rotamt_x2);

		vec128_t bx0 = v128_cat_64(u64_bit_deposit(bbl, ix_deposit_mask), u64_bit_deposit(bbl >> b_rotamt, ix_deposit_mask));

		bbl = u64_left_rotate(bbl, b_rotamt_x2);
		vec128_t bx1 = v128_cat_64(u64_bit_deposit(bbl, ix_deposit_mask), u64_bit_deposit(bbl >> b_rotamt, ix_deposit_mask));
		vec256_t vbx = v256_byterev_64(v256_cat_128(bx1, bx0));

		// process and emit a block of 32 indexees
		vl = v256_permutez_8x16_pair(ix_xlat_tab, vl);
		vl = v256_add_8(vl, v256_permutez_8x16_pair(ix_tq_addends_tab.v, vbx));
		vl = v256_xor(vl, v256_cmpgt_s8(vbx, val63_vec.v));
		vl = v256_shr2_u8(vl);
		vl = v256_sub_8(vl, v256_cmpgt_s8(vl, val32_vec.v));  // for each index, add 1 if index value is greater than 32.

		vl = v256_mul2add1_8(vl);

		vpx[0] = vl;
	}
}

/*****************************************************************************************

Perform ASTC color-ISE decoding.

The maximum length of the color-ISE is 75 bits; the maximum number of items is 18.

This corresponds to a maximum of 6 quint-blocks or 4 trit-blocks. The maximum
number of bits in any of the trit/quint encodings appears to be about 50 or so.

Unlike the index-ISE, there is no bit-reversal present.

The trit/quint formats are:

index :   range :   content         :inpbits : B-layout  : mul : desired-layout
	4 :    0..5 :  1 trit + 1 bit   : a      : 000000000 : 204 :  _aTT|____
	6 :    0..9 :  1 quint + 1 bit  : a      : 000000000 : 113 :  _aQQ|___Q
	7 :   0..11 :  1 trit + 2 bits  : ba     : b000b0bb0 :  93 :  baTT|____
	9 :   0..19 :  1 quint + 2 bits : ba     : b0000bb00 :  54 :  baQQ|___Q
   10 :   0..23 :  1 trit + 3 bits  : cba    : cb000cbcb :  44 :  cbTT|a___
   12 :   0..39 :  1 quint + 3 bits : cba    : cb0000cbc :  26 :  cbQQ|a___
   13 :   0..47 :  1 trit + 4 bits  : dcba   : dcb000dcb :  22 :  dcTT|ba__
   15 :   0..79 :  1 quint + 4 bits : dcba   : dcb0000dc :  13 :  dcQQ|ba Q
   16 :   0..95 :  1 trit + 5 bits  : edcba  : edcb0000e :  11 :  edTT|cba_
   18 :  0..159 :  1 quint + 5 bits : edcba  : edcb0000e :   6 :  edQQ|cbaQ
   19 :  0..191 :  1 trit + 6 bits  : fedcba : fedcb000f :   5 :  feTT|dcba

When unpacking trits and quints, it appears to be useful
to put the bottom 2 bits of the TQ-value at the top, and
the top bit of the quint at the bottom. The TQ-table needs to be
prepared accordingly.

From the "desired-layout", we use the top 4 bits and the bottom 4 bits
for two distinct VPSHUFB lookups, and then add them together.

For bit-formats, we the "desired-layout" is set to be a plain left-justified layout;
same flow should work for bit-formats as for the trit/quint-formats.
*****************************************************************************************/

static vec256_t astcsd_color_ise_unpack(
	const uint8_t* blk,

	uint32_t quant_level,
	uint32_t index_pair_count  // 1..9
) {

	// -------------------------------------------------------
	//   AVX2 constant-value vectors that are used unindexed
	// -------------------------------------------------------
	static const union { uint8_t b[32]; vec256_t v; } val241_vec = {{ 241,241,241,241, 241,241,241,241, 241,241,241,241, 241,241,241,241, 241,241,241,241, 241,241,241,241, 241,241,241,241, 241,241,241,241}};

	// --------------------------------------------------------------
	//   Tables for use with the ISE decode/unquantization routines
	// --------------------------------------------------------------
	static const struct {
		uint32_t xt_bitcounts[21];
		uint32_t xt_masks_high[21];
		uint64_t xt_masks_low[21];
		uint64_t deposit_masks[21];
		union { uint8_t s[16]; vec128_t v; } hightabs[21];
		union { uint8_t s[16]; vec128_t v; } lowtabs[21];
		uint32_t amask_tab[21];
		uint8_t tq_lengths[48];
		uint8_t switch_item_tab[48];
		} cxd = {

	// bits 7:0 : number of bits set for each entry in the low 64 bits of the bit-extraction masks.
	// bits 15:8 : number of bits per integer (excluding trit/quint/top-3-bits) times 8.

	//static const uint32_t xt_bitcounts[21] = {
		{
		_,_,_,_,//0,0,0,0,
		0x0827, 0x003c, 0x082d, 0x101c,
		0x0830, 0x1022, 0x1816, 0x1026,
		0x181c, 0x2012, 0x181f, 0x2018,
		0x280f, 0x201b, 0x2813, 0x300d,
		0x2818 },

	// bit-extraction masks that are used to capture trit-blocks bits and
	// quint-block bits for any trit/quint format, and top 3 bits for any bit-format.
	// bit-extraction masks for the remaining bits are obtained by negating these masks.

	// bits 75:64 of the masks
	//static const uint32_t xt_masks_high[21] =
		{
		_,_,_,_, // modes 0,1,2,3: not used
		0x56d, // mode 4: 0..5 : 1 trit + 1 bit
		0x0  , // mode 5: 0..7 : 3 bits + 0 bits
		0x3b6, // mode 6: 0..9 : 1 quint + 1 bit    <-- last format that's guaranteed to fit into 64 bits for an 18-integer sequence
		0x499, // mode 7: 0..11 : 1 trit + 2 bits
		0x6ee, // mode 8: 0..15 : 3 bits + 1 bit    <-- last format that can produce 18 integers
		0x339, // mode 9: 0..19 : 1 quint + 2 bits
		0x311, // mode 10: 0..23 : 1 trit + 3 bits
		0x739, // mode 11: 0..31 : 3 bits + 2 bits
		0x638, // mode 12: 0..39 : 1 quint + 3 bits
		0x10c, // mode 13: 0..47 : 1 trit + 4 bits
		0x0e3, // mode 14: 0..63 : 3 bits + 3 bits
		0x430, // mode 15: 0..79 : 1 quint + 4 bits
		0x182, // mode 16: 0..95 : 1 trit + 5 bits
		0x438, // mode 17: 0..127: 3 bits + 4 bits
		0x383, // mode 18: 0..159 : 1 quint + 5 bits
		0x018, // mode 19: 0..191 : 1 trit + 6 bits
		0x0e0, // mode 20: 0..255 : 3 bits + 5 bits
		},

	// bits 63:0 of the masks
	//static const uint64_t xt_masks_low[21] =
		{
		_,_,_,_,               // modes 0,1,2,3: not used
		0x6b6b5b5adad6d6b6ull, // mode 4: 0..5 : 1 trit + 1 bit
		0x0fffffffffffffffull, // mode 5: 0..7 : 3 bits + 0 bits
		0xedbb6edbb6edbb6eull, // mode 6: 0..9 : 1 quint + 1 bit    <-- last format that's guaranteed to fit into 64 bits for an 18-integer sequence
		0x33264cc9933264ccull, // mode 7: 0..11 : 1 trit + 2 bits
		0xeeeeeeeeeeeeeeeeull, // mode 8: 0..15 : 3 bits + 1 bit    <-- last format that can produce 18 integers
		0x99ccce667333999cull, // mode 9: 0..19 : 1 quint + 2 bits
		0x88c623118c462318ull, // mode 10: 0..23 : 1 trit + 3 bits
		0xce739ce739ce739cull, // mode 11: 0..31 : 3 bits + 2 bits
		0xc638c638c638c638ull, // mode 12: 0..39 : 1 quint + 3 bits
		0x308610c308610c30ull, // mode 13: 0..47 : 1 trit + 4 bits
		0x8e38e38e38e38e38ull, // mode 14: 0..63 : 3 bits + 3 bits
		0xe1861c30c3861870ull, // mode 15: 0..79 : 1 quint + 4 bits
		0x0c1060c106083060ull, // mode 16: 0..95 : 1 trit + 5 bits
		0x70e1c3870e1c3870ull, // mode 17: 0..127: 3 bits + 4 bits
		0x060e0c18383060e0ull, // mode 18: 0..159 : 1 quint + 5 bits
		0x103030206040c0c0ull, // mode 19: 0..191 : 1 trit + 6 bits
		0xe0e0e0e0e0e0e0e0ull, // mode 20: 0..255 : 3 bits + 5 bits
		},


	// deposit-masks to use for the low-order bits.
	//  - for bit-formats, fill down 000xxxxx
	//  - for trit/quint-formats, generally fill down xx00xxxx
	//    except that 'a'-bit must be kept out
	//    of the top bit position.
	//static const uint64_t deposit_masks[21] =
		{
		_,_,_,_,                // modes 0,1,2,3: not used
		0x4040404040404040ull,  // 1 trit + 1 bit
		0,                      // 3 bits + 0 bits
		0x4040404040404040ull,  // 1 quint + 1 bit
		0xc0c0c0c0c0c0c0c0ull,  // 1 trit + 2 bits
		0x1010101010101010ull,  // 3 bits + 1 bit
		0xc0c0c0c0c0c0c0c0ull,  // 1 quint + 2 bits
		0xc8c8c8c8c8c8c8c8ull,  // 1 trit + 3 bits
		0x1818181818181818ull,  // 3 bits + 2 bits
		0xc8c8c8c8c8c8c8c8ull,  // 1 quint + 3 bits
		0xccccccccccccccccull,  // 1 trit + 4 bits
		0x1c1c1c1c1c1c1c1cull,  // 3 bits + 3 bits
		0xccccccccccccccccull,  // 1 quint + 4 bits
		0xcecececececececeull,  // 1 trit + 5 bits
		0x1e1e1e1e1e1e1e1eull,  // 3 bits + 4 bits
		0xcecececececececeull,  // 1 quint + 5 bits
		0xcfcfcfcfcfcfcfcfull,  // 1 trit + 6 bits
		0x1f1f1f1f1f1f1f1full   // 3 bits + 5 bits
		},

	// The unquantization process is done by collecting bits into the 'desired'-layout described above;
	// after this collection is done, each 8-bit value can be split into two 4-bit values; the two
	// 4-bit values can then be used to look up values from two 16-entry lookup tables, which can
	// then be added together in order to obtain the final quantization result (except for
	// the part where all bits are supposed to be XORed with the 'a'-bit).
	//static const union { uint8_t s[16]; vec128_t v; } hightabs[21] =
		{
		{{  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _ }},
		{{  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _ }},
		{{  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _ }},
		{{  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _ }},
		{{  0, 51,102,  _,  0, 51,102,  _,  _,  _,  _,  _,  _,  _,  _,  _ }},
		{{  0,  _, 36,  _, 73,  _,109,  _,146,  _,182,  _,219,  _,255,  _ }},
		{{  0, 28, 56, 84,  0, 28, 56, 84,  _,  _,  _,  _,  _,  _,  _,  _ }},
		{{  0, 23, 46,  _,  0, 23, 46,  _, 69, 92,116,  _, 69, 92,116,  _ }},
		{{  0, 17, 34, 51, 68, 85,102,119,136,153,170,187,204,221,238,255 }},
		{{  0, 13, 27, 40,  0, 13, 27, 40, 67, 80, 94,107, 67, 80, 94,107 }},
		{{  0, 11, 22,  _, 33, 44, 55,  _, 66, 77, 88,  _, 99,110,121,  _ }},
		{{  0, 16, 33, 49, 66, 82, 99,115,132,148,165,181,198,214,231,247 }},
		{{  0,  6, 13, 19, 32, 39, 45, 52, 65, 71, 78, 84, 97,104,110,117 }},
		{{  0,  5, 11,  _, 32, 38, 43,  _, 65, 70, 76,  _, 97,103,108,  _ }},
		{{  0, 16, 32, 48, 65, 81, 97,113,130,146,162,178,195,211,227,243 }},
		{{  0,  3,  6,  9, 32, 35, 38, 42, 64, 67, 71, 74, 96,100,103,106 }},
		{{  0,  2,  5,  _, 32, 35, 37,  _, 64, 67, 70,  _, 96, 99,102,  _ }},
		{{  0, 16, 32, 48, 64, 80, 96,112,129,145,161,177,193,209,225,241 }},
		{{  0,  1,  3,  4, 32, 33, 35, 36, 64, 65, 67, 68, 96, 97, 99,100 }},
		{{  0,  1,  2,  _, 32, 33, 34,  _, 64, 65, 66,  _, 96, 97, 98,  _ }},
		{{  0, 16, 32, 48, 64, 80, 96,112,128,144,160,176,192,208,224,240 }},
		},

	//static const union { uint8_t s[16]; vec128_t v; } lowtabs[21] =
		{
		{{  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _ }},
		{{  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _ }},
		{{  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _ }},
		{{  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _ }},
		{{  0,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _ }},
		{{  0,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _ }},
		{{  0,113,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _ }},
		{{  0,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _ }},
		{{  0,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _ }},
		{{  0, 54,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _,  _ }},
		{{  0,  _,  _,  _,  _,  _,  _,  _,  0,  _,  _,  _,  _,  _,  _,  _ }},
		{{  0,  _,  _,  _,  _,  _,  _,  _,  8,  _,  _,  _,  _,  _,  _,  _ }},
		{{  0, 26,  _,  _,  _,  _,  _,  _,  0, 26,  _,  _,  _,  _,  _,  _ }},
		{{  0,  _,  _,  _,  0,  _,  _,  _, 16,  _,  _,  _, 16,  _,  _,  _ }},
		{{  0,  _,  _,  _,  4,  _,  _,  _,  8,  _,  _,  _, 12,  _,  _,  _ }},
		{{  0, 13,  _,  _,  0, 13,  _,  _, 16, 29,  _,  _, 16, 29,  _,  _ }},
		{{  0,  _,  0,  _,  8,  _,  8,  _, 16,  _, 16,  _, 24,  _, 24,  _ }},
		{{  0,  _,  2,  _,  4,  _,  6,  _,  8,  _, 10,  _, 12,  _, 14,  _ }},
		{{  0,  6,  0,  6,  8, 14,  8, 14, 16, 22, 16, 22, 24, 30, 24, 30 }},
		{{  0,  0,  4,  4,  8,  8, 12, 12, 16, 16, 20, 20, 24, 24, 28, 28 }},
		{{  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15 }},
		},



	// bit-mask for location of the 'a'-bit. If this bit is set in 'dlvec',
	// then all bits of the unquantized value must be negated; this only
	// really applies to formats with trits and quints in them.
	//static const union {uint32_t bx; float v; } amask_tab[21] =
		{
		_,_,_,_,
		0x40404040, // format #4
		0,
		0x40404040, // format #6
		0x40404040, // format #7
		0,
		0x40404040, // format #9
		0x08080808, // format #10
		0,
		0x08080808, // format #12
		0x04040404, // format #13
		0,
		0x04040404, // format #15
		0x02020202, // format #16
		0,
		0x02020202, // format #18
		0x01010101, // format #19
		0
		},

	// Number of valid bits in the 'high-bits' bit-vector. This is indexed
	// by number of index-pairs and whether a trit/quint encoding is used.
	//static const uint8_t tq_lengths[48] =
		{
		0, 6,12,18,24,30,36,42,48,54, _,_,_,_,_,_,  // bit-formats
		0, 4, 7,10,13,16,20,23,26,29, _,_,_,_,_,_,  // trit-formats
		0, 5,10,14,19,24,28,33,38,42, _,_,_,_,_,_,   // quint-formats
		},

	// Switch-item table. The main part of the ISE decode process
	// takes the form of a 16-entry switch-statement containing
	// an unrolled fetch-trit/quint loop. This table specifies which
	// entry points in the table to use; it is indexed by number
	// of index-pairs and whether a trit/quint encoding is used.
	//static const uint8_t switch_item_tab[48] =
		{
		0, 0, 0, 0, 0, 1, 1, 1, 1, 2, _,_,_,_,_,_, // bit-formats
		3, 3, 3, 4, 4, 5, 6, 6, 7, 8, _,_,_,_,_,_, // trit-formats
		9, 9,10,10,11,12,12,13,14,15, _,_,_,_,_,_, // quint-formats
		},
	};


	// -----------------
	//   end of tables
	// -----------------
	uint64_t blk_low = *(uint64_t *)(blk);
	uint64_t blk_high = *(uint64_t *)(blk+8);

	uint32_t blk_shamt = 0x1d1d1d11u >> ((blk_low >> 8) & 0x18); // 17 if bits 12:11 are 0, 29 otherwise.

	uint64_t ise_low = (blk_low >> (blk_shamt & 0x3f)) | (blk_high << (-blk_shamt & 0x3f));  // bits 63:0 of ISE
	uint32_t ise_high = blk_high >> 29; // bits 75:64 of ISE if 2 or more partitions, don't-care if 1 partition.

	// compute an index that is needed for the switch-selection
	// and the high-bit masking below.
	uint32_t tqm_offset = (044444400 >> quant_level) & 48; // 32 for quint-formats, 16 for trit-format, 0 for bit-formats.
	uint32_t tqm_mask_idx = tqm_offset + index_pair_count;

	uint32_t xt_bitcount = cxd.xt_bitcounts[quant_level];

	// high_bits <- trit/quint-bits, or the top 3 bits of each integer
	// low_bits <- remaining bits of each integer
	uint64_t xt_mask_low = cxd.xt_masks_low[quant_level];
	uint32_t xt_mask_high = cxd.xt_masks_high[quant_level];

	uint64_t high_bits = u12_bit_extract(ise_high, xt_mask_high);
	uint64_t low_bits = u12_bit_extract(ise_high, xt_mask_high ^ 0xfff);

	astcsd_u64_pair_t ise_low_split = u64_bit_extract_pair(ise_low, xt_mask_low, ~xt_mask_low);

	high_bits = (high_bits << (xt_bitcount  & 0x3f)) + ise_low_split.v0;
	low_bits  = (low_bits  << (-xt_bitcount & 0x3f)) + ise_low_split.v1;

	uint32_t low_bitcount_x8 = xt_bitcount >> 8;

	// Zero out 6 bits of the 'high_bits' vector above the valid bits.
	// This is needed for correct decode of trit/quint blocks; when the
	// number of integers in the ISE is not divisible by the number
	// of trits/quints in 1 block (5 and 3 respectively), then only
	// a partial block is encoded, and decoding must progress as if the
	// remainder of the block is all-0s. This can require the zeroing
	// of as much as 6 bits, in order to ensure that the top 6 bits
	// of a trit-block are all zeroed out properly.
	high_bits &= ~(63ull << cxd.tq_lengths[tqm_mask_idx]);

	uint64_t tr_deposit_mask = cxd.deposit_masks[quant_level];

	uint32_t switch_item = cxd.switch_item_tab[tqm_mask_idx];

	const struct astcsd_trit_quint_table_t *tqt = &astcsd_trit_quint_table;

	// temporaries assigned within the switch-statement, and used afterwards.
	vec128_t t;
	vec128_t bhv0 = v128_zero();
	vec128_t bhv1 = v128_zero();
	vec128_t blv1 = v128_zero();

	uint64_t bh1=0, bl1=0;

	// ----------------------------------------------
	//   main switch-statement for color-ISE decode
	// ----------------------------------------------
	switch(switch_item & 0xFU)
	{
		case 2: // bit-format, 18 indexes: decode indexes 17-16
			blv1 = v128_set_64zx(u64_bit_deposit(low_bits >> (low_bitcount_x8*2), tr_deposit_mask));
			bhv1 = v128_set_64zx(((high_bits >> 43) & 0xe0) | ((high_bits >> 38) & 0xe000));
			ASTCSD_FALLTHROUGH;
		case 1: // bit-format, more than 8 indexes: decode indexes 15-8
			bh1 = u64_bit_deposit(high_bits >> 24, 0xe0e0e0e0e0e0e0e0);
			bl1 = u64_bit_deposit(low_bits >> low_bitcount_x8, tr_deposit_mask);
			ASTCSD_FALLTHROUGH;
		case 0: // bit-format final: decode indexes 7-0.
			bhv0 = v128_cat_64(bh1, u64_bit_deposit(high_bits, 0xe0e0e0e0e0e0e0e0));
			break;

		case 8: // trit-format, 18 indexes; decode low-bits for indexes 17-16
			blv1 = v128_set_64zx(u64_bit_deposit(low_bits >> (low_bitcount_x8*2), tr_deposit_mask));
			ASTCSD_FALLTHROUGH;
		case 7: // trit-format, 16 indexes; decode trit-block #3
			t = v128_load_64zx(tqt->trit_tab + ((high_bits >> 24) & 0xff));
			bhv0 = v128_shlb(t, 15);
			bhv1 = v128_shrb(t, 1);
			ASTCSD_FALLTHROUGH;
		case 6: // trit-format, 12 or more indexes; decode trit-block #2
			t = v128_load_64zx(tqt->trit_tab + ((high_bits >> 16) & 0xff));
			bhv0 = v128_or(bhv0, v128_shlb(t, 10));
			ASTCSD_FALLTHROUGH;
		case 5: // trit-format, 10 or more indexes; decode low-bits for indexes 15-8
			bl1 = u64_bit_deposit(low_bits >> low_bitcount_x8, tr_deposit_mask);
			ASTCSD_FALLTHROUGH;
		case 4: // trit-format, 6 or more indexes; decode trit-block #1
			t = v128_load_64zx(tqt->trit_tab + ((high_bits >> 8) & 0xff));
			bhv0 = v128_or(bhv0, v128_shlb(t, 5));
			ASTCSD_FALLTHROUGH;
		case 3: // trit-format final: decode trit-block #0
			t = v128_load_64zx(tqt->trit_tab + (high_bits & 0xFF));
			bhv0 = v128_or(bhv0, t);
			break;

		case 15: // quint-format, 18 indexes; decode low-bits for indexes 17-16
			blv1 = v128_set_64zx(u64_bit_deposit(low_bits >> (low_bitcount_x8*2), tr_deposit_mask));
			ASTCSD_FALLTHROUGH;
		case 14: // quint-format, 16 or more indexes; decode quint-block #5
			t = v128_load_32zx(tqt->quint_tab + ((high_bits >> 35) & 0x7f));
			bhv0 = v128_shlb(t, 15);
			bhv1 = v128_shrb(t, 1);
			ASTCSD_FALLTHROUGH;
		case 13: // quint-format, 14 or more indexes; decode quint-block #4
			t = v128_load_32zx(tqt->quint_tab + ((high_bits >> 28) & 0x7f));
			bhv0 = v128_or(bhv0, v128_shlb(t, 12));
			ASTCSD_FALLTHROUGH;
		case 12: // quint-format, 10 or more indexes; deocde quint-block #3 and low-bits for indexes 15-8
			t = v128_load_32zx(tqt->quint_tab + ((high_bits >> 21) & 0x7f));

			bhv0 = v128_or(bhv0, v128_shlb(t, 9));
			bl1 = u64_bit_deposit(low_bits >> low_bitcount_x8, tr_deposit_mask);
			ASTCSD_FALLTHROUGH;
		case 11: // quint-format, 8 or more indexes; decode quint-block #2
			t = v128_load_32zx(tqt->quint_tab + ((high_bits >> 14) & 0x7f));
			bhv0 = v128_or(bhv0, v128_shlb(t, 6));
			ASTCSD_FALLTHROUGH;
		case 10: // quint-format, 4 or more indexes: decode quint-block #1
			t = v128_load_32zx(tqt->quint_tab + ((high_bits >> 7) & 0x7f));
			bhv0 = v128_or(bhv0, v128_shlb(t, 3));
			ASTCSD_FALLTHROUGH;
		case 9: // quint-format final: decode quint-block #0
			t = v128_load_32zx(tqt->quint_tab + (high_bits & 0x7f));
			bhv0 = v128_or(bhv0, t);
			break;
	}

	// decode low-bits for indexes 7-0.
	uint64_t bl0 = u64_bit_deposit(low_bits, tr_deposit_mask);

	// prepare lookup tables for postprocessing.
	vec256_t hightab = v256_set_128rep(cxd.hightabs[quant_level].v);
	vec256_t lowtab  = v256_set_128rep(cxd.lowtabs[quant_level].v);

	vec256_t amask = v256_load_32rep(cxd.amask_tab + quant_level);


	// construct the 'desired-layout'-vector.
	vec256_t bhv = v256_cat_128(bhv1, bhv0); // "high-bits": top 3 bits of index and/or trit/quint bits
	vec256_t blv = v256_cat_128(blv1, v128_cat_64(bl1, bl0)); // "low-bits": remaining bits of index

	// bitwise-AND for trit/quint values to get rid of low-order bit detritus
	// that is present in the trit/quint tables. (This detritus is necessary
	// in order to be able to reuse the same tables for index-ISE and color-ISE decoding.)
	bhv = v256_and(bhv,  val241_vec.v);

	vec256_t dlvec = v256_or(bhv, blv); // 'desired-layout' vector.

	// perform unquantization of the 'desired-layout'
	vec256_t lowvec  = v256_and15_8(dlvec);
	vec256_t highvec = v256_shr4_u8(dlvec);

	vec256_t xorvec = v256_cmpgt_s8(v256_and(dlvec, amask), v256_zero());
	vec256_t sum    = v256_add_8(v256_permute_8x16_pair(hightab, highvec), v256_permute_8x16_pair(lowtab, lowvec));

	// conditionally bit-negate based on the 'a'-bit.
	return v256_xor(sum, xorvec);
}

/*******************************************************************************
  ASTC 2D index-grid infill routine.

  What we take as input are:
	* a decoded index-ISE (1 byte per index)
	* block-size (Xdim, Ydim)
	* grid-size (Xdim, Ydim, dual-plane-flag)
  What we want to produce as output is a set of vectors such that:
	* For each pixel, we have 4 bytes laid out as:
		byte 3: <second-plane-weight>
		byte 2: <64 minus second-plane-weight>
		byte 1: <first-plaen-weight>
		byte 0: <64 minus first-plane-weight>
	  bytes 3:2 don't matter for things that aren't dual-index-planes
	* We want the output to be 2x12 vectors, each 32 bytes.
	  The first set of 12 vectors contains items for the
	  first 8 texels in a scanline. The second set of 12
	  vectors contains items for the remaining texels
	  in the scanline.
*******************************************************************************/

static void astcsd_index_infill_2d(
	const astc_simd_decode_processed_params_t* bd,
	const astcsd_block_mode_desc_t *bmd,
	const uint8_t *indexes,  // unquantized indexes, modified by 2*X+1
	vec256_t *res
) {
	uint32_t y;

	// -------------------------------------------------------
	//   AVX2 constant-value vectors that are used unindexed
	// -------------------------------------------------------
	static const union { uint8_t  s[32]; vec256_t v; } byterep4_swz = {{ 0,0,0,0,  1,1,1,1,  2,2,2,2,  3,3,3,3,  4,4,4,4,  5,5,5,5,  6,6,6,6,  7,7,7,7 }};
	static const union { uint16_t s[16]; vec256_t v; } vh175        = {{ 175,175,175,175,  175,175,175,175,  175,175,175,175,  175,175,175,175 }};
	static const union { uint16_t s[16]; vec256_t v; } vh2          = {{ 2,2,2,2, 2,2,2,2, 2,2,2,2, 2,2,2,2 }};
	static const union { uint8_t  s[32]; vec256_t v; } xn_mask      = {{ 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0 }};
	static const union { uint16_t s[16]; vec256_t v; } xft_mask     = {{ 0xf000, 0xf000, 0xf000, 0xf000, 0xf000, 0xf000, 0xf000, 0xf000, 0xf000, 0xf000, 0xf000, 0xf000, 0xf000, 0xf000, 0xf000, 0xf000 }};
	static const union { uint8_t  s[32]; vec256_t v; } xf_rsubv     = {{ 64, 0, 64, 0,  64, 0, 64, 0,  64, 0, 64, 0,  64, 0,  64, 0,   64, 0, 64, 0,  64, 0, 64, 0,  64, 0, 64, 0, 64, 0, 64, 0}};
	static const union { uint8_t  s[32]; vec256_t v; } xh_mask      = {{ 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff }};
	static const union { uint8_t  s[32]; vec256_t v; } vsmask2      = {{ 64,0xff, 64,0xff, 64,0xff, 64,0xff,  64,0xff, 64,0xff, 64,0xff, 64,0xff,  64,0xff, 64,0xff, 64,0xff, 64,0xff,  64,0xff, 64,0xff, 64,0xff, 64,0xff }};
	static const union { uint8_t  s[32]; vec256_t v; } vsmask1      = {{ 63,0xff, 63,0xff, 63,0xff, 63,0xff,  63,0xff, 63,0xff, 63,0xff, 63,0xff,  63,0xff, 63,0xff, 63,0xff, 63,0xff,  63,0xff, 63,0xff, 63,0xff, 63,0xff }};
	static const union { uint16_t s[16]; vec256_t v; } or_lowv      = {{ 0xffff, 0, 0xffff, 0, 0xffff, 0, 0xffff, 0, 0xffff, 0, 0xffff, 0, 0xffff, 0, 0xffff, 0 }};
	static const union { uint16_t s[16]; vec256_t v; } or_highv     = {{ 0, 0xffff, 0, 0xffff, 0, 0xffff, 0, 0xffff, 0, 0xffff, 0, 0xffff, 0, 0xffff, 0, 0xffff }};

	// ------------------------------------------------------
	//   Indexable lookup tables. These are folded into one
	//   big data structure in order to limit the number
	//   of LEA instructions needed.
	// ------------------------------------------------------
	static const struct astcsd_infill2d_constants {
		uint32_t tosv_addends[2];
		uint64_t fvals_tab[160];
		uint64_t rsv_tab[160];
		uint64_t shs_tab[320];
		} icp = {

	// texel-index addends for each byte of a 4-byte lane:
	//    * no-dual-planes: 0,1,0,1
	//    * dual-planes:    0,2,1,3
	//static const uint32_t tosv_addends[2] =
		{ 0x01000100, 0x03010200 },

	// -----------------------------------------------------------
	//   constants and tables for use with the 2d infill-routine
	// -----------------------------------------------------------

	// Table of fractional-offsets for each texel in a row of texels.
	// For the first 8 texels in a row, the fractional-offset of texel N is found at bits [8*N+3:8*N]
	// For the remaining texels in a row with more than 8 texels, the fractional-offset
	// is found at bits [8*(N-8)+7 : 8*(N-8)+4] and [8*(N-8)+7+32 : 8*(N-8)+4+32]
	// (this duplication is helpful to avoid some duplication instructions when preparing
	// the 2-rows-at-a-time loop for handling more than 8 texels)
	// Indexing into the table is done with [grid_ydim + 16*(block_ydim-3)]
	//static const uint64_t fvals_tab[160] = {
		{
		0, 0, 0x0000000000000800ull, 0x0000000000000000ull, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0x00000000000B0500ull, 0x0000000000050B00ull, 0x0000000000000000ull, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0x000000000C080400ull, 0x0000000008000800ull, 0x0000000004080C00ull, 0x0000000000000000ull, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0x0000000D0A060300ull, 0x0000000A030D0600ull, 0x000000060D030A00ull, 0x00000003060A0D00ull, 0x0000000000000000ull, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0x00000D0B08050300ull, 0x00000B05000B0500ull, 0x0000080008000800ull, 0x0000050B00050B00ull, 0x00000305080B0D00ull, 0x0000000000000000ull, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0x000E0B0907050200ull, 0x000B07020E090500ull, 0x0009020B050E0700ull, 0x00070E050B020900ull, 0x0004090E02070B00ull, 0x00020407090B0E00ull, 0x0000000000000000ull, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0x0E0C0A0806040200ull, 0x0C0804000C080400ull, 0x0A040E08020C0600ull, 0x0800080008000800ull, 0x060C02080E040A00ull, 0x04080C0004080C00ull, 0x020406080A0C0E00ull, 0x0000000000000000ull, 0, 0, 0, 0, 0, 0,
		0, 0, 0x0C0B09E7050402E0ull, 0x090502DE0B0704D0ull, 0x05000BB5000B05B0ull, 0x020B049D050E0790ull, 0x0E050D740B020970ull, 0x0B00056B00050B60ull, 0x070B0E4205090C40ull, 0x040607290B0D0E20ull, 0x0000000000000000ull, 0, 0, 0, 0, 0,
		0, 0, 0x0B0AE8D60503E2D0ull, 0x0603D0AD0A06D3A0ull, 0x010DB8630E0AB560ull, 0x0D06903A030D9630ull, 0x0800880008008800ull, 0x030960D60D036AD0ull, 0x0E03489D01064B90ull, 0x090D3063060A3D60ull, 0x04F618390BFD1E30ull, 0x00F0F00000F0F000ull, 0, 0, 0, 0,
		0, 0, 0x0AF9D7C604F3D1C0ull, 0x04D1AF7C09D6A370ull, 0x0FCA76310DC97430ull, 0x09A34DF701AC46F0ull, 0x039C14AD069F17A0ull, 0x0D74EC630A71E960ull, 0x076DC3190F64CA10ull, 0x01469ADF03479CD0ull, 0x0C3E6194073A6D90ull, 0x0617394A0C1D3F40ull, 0x0000000000000000ull, 0, 0
		},

	// Table of offset-vectors for t-axis.
	// For destination row N of texels, bits[4*N+6 : 4*N+3] of entries from this table contain the index of the first source row.
	// Indexing into the table is done with [grid_ydim + 16*(block_ydim-3)]
	//static const uint64_t rsv_tab[160] = {
		{
		0, 0, 0x0000000000000800ull, 0x0000000000001080ull, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0x0000000000008000ull, 0x0000000000010800ull, 0x0000000000019080ull, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0x0000000000080000ull, 0x0000000000108800ull, 0x0000000000190800ull, 0x0000000000219080ull, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0x0000000000800000ull, 0x0000000001088000ull, 0x0000000001908800ull, 0x0000000002190800ull, 0x0000000002A19080ull, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0x0000000008000000ull, 0x0000000010888000ull, 0x0000000019108800ull, 0x0000000021910800ull, 0x000000002A190800ull, 0x0000000032A19080ull, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0x0000000080000000ull, 0x0000000108880000ull, 0x0000000191088000ull, 0x0000000219108800ull, 0x00000002A1910800ull, 0x000000032A190800ull, 0x00000003B2A19080ull, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0x0000000800000000ull, 0x0000001088880000ull, 0x0000001910888000ull, 0x0000002199108800ull, 0x0000002A19908800ull, 0x00000032A1990800ull, 0x0000003B2A190800ull, 0x00000043B2A19080ull, 0, 0, 0, 0, 0, 0,
		0, 0, 0x0000008000000000ull, 0x0000010888800000ull, 0x0000019110888000ull, 0x0000021991088000ull, 0x000002A199108800ull, 0x0000032A21910800ull, 0x000003B2A1990800ull, 0x0000043B2A190800ull, 0x000004C3B2A19080ull, 0, 0, 0, 0, 0,
		0, 0, 0x0000080000000000ull, 0x0000108888800000ull, 0x0000191108880000ull, 0x0000219911088000ull, 0x00002A2199108800ull, 0x000032A219908800ull, 0x00003B2A21910800ull, 0x000043B2A2190800ull, 0x0000443B2A190800ull, 0x00004C43B2A19080ull, 0, 0, 0, 0,
		0, 0, 0x0000800000000000ull, 0x0001088888000000ull, 0x0001911088880000ull, 0x0002199110888000ull, 0x0002A21991088000ull, 0x00032A2199108800ull, 0x0003B2AA19908800ull, 0x00043B2AA1910800ull, 0x0004C3B2A2190800ull, 0x00054C3B2A190800ull, 0x0005D4C3B2A19080ull, 0, 0, 0
		},

	// Table of integer texel offsets for each texel in a row of texels.
	// For the first 8 texels in a row, the fractional-offset of texel N is found at bits [8*N+3:8*N]
	// For the remaining texels in a row with more than 8 texels, the fractional-offset
	// is found at bits [8*(N-8)+7 : 8*(N-8)+4] and [8*(N-8)+7+32 : 8*(N-8)+4+32]
	// Indexing into the table is done with [2*(grid_ydim + 16*(block_ydim-3)) + dual_index_planes]
	//static const uint64_t shs_tab[320] = {
		{
						 0ull,                  0ull,                  0ull,                  0ull, 0x0000000000010000ull, 0x0000000000020000ull, 0x0000000000020100ull, 0x0000000000040200ull,
						 0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,
						 0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,
						 0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,

						 0ull,                  0ull,                  0ull,                  0ull, 0x0000000001000000ull, 0x0000000002000000ull, 0x0000000002010000ull, 0x0000000004020000ull,
		0x0000000003020100ull, 0x0000000006040200ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,
						 0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,
						 0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,

						 0ull,                  0ull,                  0ull,                  0ull, 0x0000000100000000ull, 0x0000000200000000ull, 0x0000000201010000ull, 0x0000000402020000ull,
		0x0000000302010000ull, 0x0000000604020000ull, 0x0000000403020100ull, 0x0000000806040200ull,                  0ull,                  0ull,                  0ull,                  0ull,
						 0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,
						 0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,

						 0ull,                  0ull,                  0ull,                  0ull, 0x0000010000000000ull, 0x0000020000000000ull, 0x0000020101000000ull, 0x0000040202000000ull,
		0x0000030201010000ull, 0x0000060402020000ull, 0x0000040302010000ull, 0x0000080604020000ull, 0x0000050403020100ull, 0x00000A0806040200ull,                  0ull,                  0ull,
						 0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,
						 0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,

						 0ull,                  0ull,                  0ull,                  0ull, 0x0001000000000000ull, 0x0002000000000000ull, 0x0002010101000000ull, 0x0004020202000000ull,
		0x0003020201010000ull, 0x0006040402020000ull, 0x0004030202010000ull, 0x0008060404020000ull, 0x0005040302010000ull, 0x000A080604020000ull, 0x0006050403020100ull, 0x000C0A0806040200ull,
						 0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,
						 0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,

						 0ull,                  0ull,                  0ull,                  0ull, 0x0100000000000000ull, 0x0200000000000000ull, 0x0201010100000000ull, 0x0402020200000000ull,
		0x0302020101000000ull, 0x0604040202000000ull, 0x0403020201010000ull, 0x0806040402020000ull, 0x0504030202010000ull, 0x0A08060404020000ull, 0x0605040302010000ull, 0x0C0A080604020000ull,
		0x0706050403020100ull, 0x0E0C0A0806040200ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,
						 0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,

						 0ull,                  0ull,                  0ull,                  0ull, 0x0000001000000010ull, 0x0000002000000020ull, 0x0101012100000020ull, 0x0202024200000040ull,
		0x0202013101000030ull, 0x0404026202000060ull, 0x0303024201010040ull, 0x0606048402020080ull, 0x0403035201010050ull, 0x080606A4020200A0ull, 0x0504036302010060ull, 0x0A0806C6040200C0ull,
		0x0605047302010070ull, 0x0C0A08E6040200E0ull, 0x0706058403020180ull, 0x0E0C0A8806040280ull,                  0ull,                  0ull,                  0ull,                  0ull,
						 0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,

						 0ull,                  0ull,                  0ull,                  0ull, 0x0000100000001000ull, 0x0000200000002000ull, 0x0101211000002010ull, 0x0202422000004020ull,
		0x0202312101003020ull, 0x0404624202006040ull, 0x0302423101004030ull, 0x0604846202008060ull, 0x0303524201015040ull, 0x0606A4840202A080ull, 0x0404635202016050ull, 0x0808C6A40402C0A0ull,
		0x0504736302017060ull, 0x0A08E6C60402E0C0ull, 0x0605847302018070ull, 0x0C0A886604028060ull, 0x0706958403029180ull, 0x0E0CAA880604A280ull,                  0ull,                  0ull,
						 0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,

						 0ull,                  0ull,                  0ull,                  0ull, 0x0010000000100000ull, 0x0020000000200000ull, 0x0121111000201010ull, 0x0242222000402020ull,
		0x0231212100302020ull, 0x0462424200604040ull, 0x0242323101403030ull, 0x0484646202806060ull, 0x0353424201514040ull, 0x06A6848402A28080ull, 0x0463534201615040ull, 0x08C6A68402C2A080ull,
		0x0474635202716050ull, 0x08E8C6A404E2C0A0ull, 0x0584746302817060ull, 0x0A88684604826040ull, 0x0685847302818070ull, 0x0C8A886604828060ull, 0x0796858403928180ull, 0x0EAC8A8806A48280ull,
						 0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,

						 0ull,                  0ull,                  0ull,                  0ull, 0x1000000010000000ull, 0x2000000020000000ull, 0x2111101020101010ull, 0x4222202040202020ull,
		0x3121212130202020ull, 0x6242424260404040ull, 0x4232312141303020ull, 0x8464624282606040ull, 0x5342423151404030ull, 0xA6848462A2808060ull, 0x6353424261514040ull, 0xC6A68484C2A28080ull,
		0x7463535271615050ull, 0xE8C6A6A4E2C2A0A0ull, 0x8574635282716050ull, 0x8A68462484624020ull, 0x9584746392817060ull, 0xAA886846A4826040ull, 0xA6958473A2918070ull, 0xCCAA8866C4A28060ull,
		0xB7A69584B3A29180ull, 0xEECCAA88E6C4A280ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,                  0ull,
		}
	};

	// ---------------------
	//   end of constants
	// ---------------------

	const struct astcsd_infill2d_constants *ic = &icp;

	uint32_t block_xdim = bd->block_xdim;
	uint32_t block_ydim = bd->block_ydim;

	uint32_t grid_xdim = bmd->grid_xdim;
	uint32_t grid_ydim = bmd->grid_ydim & 0xF;
	uint32_t dual_index_planes = bmd->dual_plane_en;

	// -------------------------------------------------------------------------
	//   step 1: unpack the index-array into a series of rows
	//
	// Each row of input indexes is unpacked into a fixed-size row of 16 bytes.
	// For cases where a row of input indexes is wider than 16 indexes,
	// only the first 16 indexes are actually copied. (The remaining indexes
	// are copied in step 4 in this case.)
	// -------------------------------------------------------------------------

	// buffer to contain rows of the index-grid.
	vec128_t grid_rows[13];

	uintptr_t xstride = bmd->grid_xstride; // grid_xdim << dual_plane_en

	switch(grid_ydim)
	{
		case 15:
		case 14:
		case 13:
		case 12: grid_rows[11] = v128_load_unaligned(indexes + xstride*11u); ASTCSD_FALLTHROUGH;
		case 11: grid_rows[10] = v128_load_unaligned(indexes + xstride*10u); ASTCSD_FALLTHROUGH;
		case 10: grid_rows[ 9] = v128_load_unaligned(indexes + xstride* 9u); ASTCSD_FALLTHROUGH;
		case  9: grid_rows[ 8] = v128_load_unaligned(indexes + xstride* 8u); ASTCSD_FALLTHROUGH;
		case  8: grid_rows[ 7] = v128_load_unaligned(indexes + xstride* 7u); ASTCSD_FALLTHROUGH;
		case  7: grid_rows[ 6] = v128_load_unaligned(indexes + xstride* 6u); ASTCSD_FALLTHROUGH;
		case  6: grid_rows[ 5] = v128_load_unaligned(indexes + xstride* 5u); ASTCSD_FALLTHROUGH;
		case  5: grid_rows[ 4] = v128_load_unaligned(indexes + xstride* 4u); ASTCSD_FALLTHROUGH;
		case  4: grid_rows[ 3] = v128_load_unaligned(indexes + xstride* 3u); ASTCSD_FALLTHROUGH;
		case  3: grid_rows[ 2] = v128_load_unaligned(indexes + xstride* 2u); ASTCSD_FALLTHROUGH;
		case  2: grid_rows[ 1] = v128_load_unaligned(indexes + xstride* 1u); ASTCSD_FALLTHROUGH;
		case  1: grid_rows[ 0] = v128_load_unaligned(indexes);               ASTCSD_FALLTHROUGH;
		case  0: break;
	}

	// the main loop can overread one row, so we will write zeroes
	// to an extra row in order to ensure that no uninitialized data are read.
	grid_rows[grid_ydim] = v128_zero();

	// Do a simplified infill-procedure if we have a 1:1 texel-index mapping.
	// This gives an apparent decode speedup of about 9% or so for
	// the 4x4 block size, diminishing with larger sizes.
	if (bd->block_xpydim == bmd->grid_xpydim)
	{
		vec256_t *resptr = res;
		uint32_t i;
		static const int8_t add_vec0[4] = { 126,-1,126,-1 };

		static const union { uint8_t s[32]; vec256_t v; } repvecs[2] = {
			{{ 0,0,0,0,  1,1,1,1,  2,2,2,2,  3,3,3,3,  4,4,4,4,  5,5,5,5,  6,6,6,6,  7,7,7,7 }},
			{{ 0,0,1,1,  2,2,3,3,  4,4,5,5,  6,6,7,7,  8,8,9,9,  10,10,11,11,  12,12,13,13,  14,14,15,15 }} };

		vec256_t perm_vec = repvecs[dual_index_planes].v;
		vec256_t add_vec = v256_load_32rep((const uint32_t *)add_vec0);

		for (i=0;i<block_ydim;i++)
		{
			vec256_t v0 = v256_set_128rep(grid_rows[i]);
			v0 = v256_permute_8x16_pair(v0, perm_vec);
			v0 = v256_add_8(v0, add_vec);
			v0 = v256_xor(v0, xn_mask.v);
			v0 = v256_shri_u16(v0, 1);
			*resptr ++ = v0;
		}
		return;
	}

	// -----------------------------------------------------------------------------------
	//   step 2: perform table lookup for per-axis data, and unpack to a representation
	//           that is maximally efficient for subsequent calculation steps.
	// -----------------------------------------------------------------------------------

	uint32_t so_idx = grid_xdim + bd->block_xdim_m3_x16;
	uint32_t to_idx = grid_ydim + bd->block_ydim_m3_x16;
	uint32_t so2_idx = 2*so_idx + dual_index_planes;
	uint64_t rsv = ic->rsv_tab[ to_idx ];   // texel-row-offset-vector for t-axis

	vec256_t tosv_addends_ext = v256_load_32rep(ic->tosv_addends + dual_index_planes);

	vec256_t ft_vecbase = v256_load_64rep(ic->fvals_tab + to_idx); // packed fractional-coordinates for t-axis
	vec256_t fs_vecbase = v256_load_64rep(ic->fvals_tab + so_idx); // packed fractional-coordinates for s-axis
	vec256_t shuf_vecbase = v256_load_64rep(ic->shs_tab + so2_idx); // packed texel-offsets-in-row vector for s-axis

	ft_vecbase = v256_permute_8x16_pair(ft_vecbase, byterep4_swz.v);
	fs_vecbase = v256_permute_8x16_pair(fs_vecbase, byterep4_swz.v);
	shuf_vecbase = v256_permute_8x16_pair(shuf_vecbase, byterep4_swz.v);

	// ft-values, for interpolation weights between rows.
	// these values are written as vectors but read as scalars - 4 bytes per row.
	// For each row, the 4 bytes have the following layout:
	//  - byte lanes 0 & 2 : lane <- 2
	//  - byte lanes 1 & 3 : lane <- (ft << 4)
	// rows 12 to 15 are don't-care; they are a duplicate of rows 8-11.
	union { uint32_t s[16]; vec256_t v[2]; } ftvec;  // processed-ft-values vector
	union { uint32_t s[16]; vec256_t v[2]; } ftvec2;  // processed-ft-values vector, shifted down by 2 bytes
	vec256_t ftv0 = v256_or(v256_and(v256_shli_u16(ft_vecbase, 4), xft_mask.v), vh2.v);
	vec256_t ftv1 = v256_or(v256_and(ft_vecbase,                     xft_mask.v), vh2.v);

	ftvec.v[0] = ftv0;
	ftvec.v[1] = ftv1;
	ftvec2.v[0] = v256_shrb_pair(ftv0, 2);
	ftvec2.v[1] = v256_shrb_pair(ftv1, 2);

	// fs-values, for interpolation weights between two texels in a row
	//  - fsvec0: fs-values for first 8 texels in a rows
	//  - fsvec1: shuffle-patterns for last 4 texels in a row, duplicated 2x
	// within each 4-byte lane:
	//  - byte lanes 0 & 2 :  lane <- (16-fs) | 64
	//  - byte lanes 1 & 3 :  lane <- fs
	vec256_t fsvec0 = v256_and15_8(fs_vecbase);
	vec256_t fsvec1 = v256_and15_8(v256_shri_u16(fs_vecbase, 4));
	fsvec0 = v256_xor(v256_add_8(fsvec0, vh175.v), xn_mask.v);
	fsvec1 = v256_xor(v256_add_8(fsvec1, vh175.v), xn_mask.v);

	// Shuffle-patterns, to pick indexes to use for texels in a row:
	//  - shuf_pat0: shuffle-patterns for first 8 texels in a row.
	//  - shuf_pat1: shuffle-patterns for last 4 texels in a row, duplicated 2x
	vec256_t shuf_pat0 = v256_and15_8(v256_add_8(shuf_vecbase,                     tosv_addends_ext));
	vec256_t shuf_pat1 = v256_and15_8(v256_add_8(v256_shri_u16(shuf_vecbase, 4), tosv_addends_ext));

	// ------------------------------------------------------------------------------------------------
	//   step 3: process each row in order to compute final infilled indexes (first 8 texels per row)
	// ------------------------------------------------------------------------------------------------

	vec256_t vs1 = v256_and(fsvec0, vsmask1.v); // fs-value vector #1
	vec256_t vs2 = v256_and(fsvec0, vsmask2.v); // fs-value vector #2

	uint64_t rsv2 = rsv;

	const uint8_t *gridptr = (const uint8_t *)grid_rows;
	vec256_t *resptr = res;

	if (dual_index_planes)
	{
		// for dual index planes, we process data one row at a time.
		// In this case, the low 16 bits of each 32-bit lane are used to process
		// indexes for the first index plane, and the high 16 bits are used
		// to process the second index plane.

		for (y=0;y<block_ydim;y++)
		{
			uintptr_t rowsel = rsv2 & 0x78;  // dword-offset of the grid-row that we wish to read.
			rsv2 >>= 4;

			vec256_t vt = v256_load_32rep(ftvec.s + y);  // broadcast of 'ft'-value for the current row.

			// fetch 2 rows of indexes; these are the rows that we will interpolate between.
			vec256_t row1 = v256_set_128rep(*(const vec128_t *)(gridptr + 2*rowsel + 16));
			vec256_t row0 = v256_set_128rep(*(const vec128_t *)(gridptr + 2*rowsel)      );

			// perform shuffles to line up indexes within each row
			row0 = v256_permute_8x16_pair(row0, shuf_pat0);
			row1 = v256_permute_8x16_pair(row1, shuf_pat0);

			// prepare index-infill-weights
			vec256_t a = v256_mul_add_pair_u16(vt, vs2); // 16-bit lane: bits[15:8] <- ((fs * ft) + 8) >> 4  = w11
			vec256_t c = v256_shri_u16(vt, 12);          // 16-bit lane: bits[7:0] <- ft
			vec256_t b = v256_shri_u16(a, 8);            // bits[7:0] <- w11
			c = v256_sub_16(c, b);                       // bits[7:0] <- ft - w11
			a = v256_and(a, xh_mask.v);
			a = v256_or(a, c);                           // bits[15:8] <- w11, bits[7:0] <- ft - w11
			b = v256_sub_8(vs1, a);

			// perform the multiplications for index infill, and add together the results; add-8 can be pushed out to an add-0.5 on indexes before infill.
			a = v256_mul_add_pair_u16(row1, a);
			b = v256_mul_add_pair_u16(row0, b);
			a = v256_add_16(a, b);

			// right-shift multiply-result by 4 bits, put it in upper 8-bit sublane of 16-bit lane, then compute 64-x in the lower lane.
			b = v256_shri_u16(a, 5);
			a = v256_sub_8(xf_rsubv.v, b);
			a = v256_interleave_low_u16(b, a);

			*resptr++ = a;
		}
	}
	else
	{
		// if we don't have dual-index-planes, then we will process rows 2 at a time in order to get
		// as much utilization of the AVX2 vector units as possible; the calculation is interleaved so that
		// the low 16 bits of each 32-bit lane is the first row, and the high 16 bits are the second row.
		// This results in a loop that takes 38 instructions (5 loads, 2 stores, 21 avx-arith, 10 int-arith)
		// while the 1-row-at-a-time loop takes 27 instructions. This loop will write 1 redundant row for
		// odd-height blocks.

		// use bit-ORing to flag lanes of the shuffle-pattern that should return 0
		// rather than an acutal byte-lane of the source argument.
		//  * shuf_pat_low: shuffle-pattern to select texels for low 2 bytes of 4-byte lane (first row)
		//  * shuf_pat_high: shuffle-pattern to select texels for high 2 bytes of 4-byte lane (second row)
		vec256_t shuf_pat_low  = v256_or(shuf_pat0, or_highv.v);
		vec256_t shuf_pat_high = v256_or(shuf_pat0, or_lowv.v);

		for (y=0;y<block_ydim;y+=2)
		{
			uintptr_t rowsel_l = rsv2 & 0x78;  // dword-offset of the grid-row that we wish to read.
			uintptr_t rowsel_h = (rsv2>>4) & 0x78;  // dword-offset of the grid-row that we wish to read.
			rsv2 >>= 8;

			// broadcast of 'ft'-values for the current pair of rows
			// 'ftvec2' is prepared so that the 4-byte lane contains ft-values for 2 rows.
			vec256_t vt = v256_load_32rep(ftvec2.s + y);

			// fetch 2x2 rows of indexes
			vec256_t row1l = v256_set_128rep(*(const vec128_t *)(gridptr + 2*rowsel_l + 16));
			vec256_t row0l = v256_set_128rep(*(const vec128_t *)(gridptr + 2*rowsel_l)      );

			vec256_t row1h = v256_set_128rep(*(const vec128_t *)(gridptr + 2*rowsel_h + 16));
			vec256_t row0h = v256_set_128rep(*(const vec128_t *)(gridptr + 2*rowsel_h)      );

			// perform shuffles to line up indexes within each row
			vec256_t row1 = v256_or(v256_permute_8x16_pair(row1l, shuf_pat_low),  v256_permute_8x16_pair(row1h, shuf_pat_high));
			vec256_t row0 = v256_or(v256_permute_8x16_pair(row0l, shuf_pat_low),  v256_permute_8x16_pair(row0h, shuf_pat_high));

			// prepare index-infill-weights
			vec256_t a = v256_mul_add_pair_u16(vt, vs2); // 16-bit lane: bits[15:8] <- ((fs * ft) + 8) >> 4  = w11
			vec256_t c = v256_shri_u16(vt, 12);          // 16-bit lane: bits[7:0] <- ft
			vec256_t b = v256_shri_u16(a, 8);            // bits[7:0] <- w11
			c = v256_sub_16(c, b);                       // bits[7:0] <- ft - w11
			a = v256_and(a, xh_mask.v);                  // bits[7:0] <- 0
			a = v256_or(a, c);                           // bits[15:8] <- w11, bits[7:0] <- ft - w11
			b = v256_sub_8(vs1, a);

			// perform the multiplications for index infill, and add together the results; add-8 can be pushed out to an add-0.5 on indexes before infill.
			a = v256_mul_add_pair_u16(row1, a);
			b = v256_mul_add_pair_u16(row0, b);
			a = v256_add_16(a, b);

			// right-shift multiply-result by 4 bits, put it in upper 8-bit sublane of 16-bit lane, then compute 64-x in the lower lane.
			b = v256_shri_u16(a, 5);
			a = v256_sub_8(xf_rsubv.v, b);
			a = v256_interleave_low_u16(b, a);

			resptr[0] = a;
			a = v256_shri_u32(a, 16);
			resptr[1] = a;
			resptr += 2;
		}
	}

	// -----------------------------------------------------------------------------------------------
	//   step 4: process each row in order to compute final infilled indexes (last 4 texels per row)
	// -----------------------------------------------------------------------------------------------

	// for block dimensions greater than 8, we need additional processing
	// in order to cover the remaining texels. We do this for two rows at a time, so that we
	// can use full-width 256-bit vectors for the main calculations. This does 'lock in' the
	// routine, in that it cannot easily be extended to block sizes wider than 12 texels;
	// this is not considered a problem.
	if (block_xdim > 8)
	{

		// if we have more than 16 indexes per row, then, for the first three rows,
		// copy indexes 23:16 of the row into the 16-byte per-row buffer. This is expected
		// to happen only rarely; it can only happen with dual-index-planes with
		// a grid dimension greater than 8.
		// (The reason we need copying for only 3 rows is that ASTC limits the number
		// of indexes per block to 64, so it's not possible to have more than 16 indexes
		// per row *and* more than 3 rows at the same time.)
		if (xstride > 16)
		{
			grid_rows[0] = v128_load_unaligned(indexes + 8);
			grid_rows[1] = v128_load_unaligned(indexes + 8 + xstride);
			grid_rows[2] = v128_load_unaligned(indexes + 8 + 2*xstride);
		}


		vs1 = v256_and(fsvec1, vsmask1.v);
		vs2 = v256_and(fsvec1, vsmask2.v);

		for (y=0;y<block_ydim;y+=2)
		{
			uintptr_t rowsel0 = rsv & 0x78;        // dword-offsets of the grid-rows that we wish to read
			uintptr_t rowsel1 = (rsv >> 4) & 0x78;
			rsv >>= 8;

			// fetch and expand ft-values for the two rows
			vec256_t vt = v256_cat_128(
				v128_load_32rep(ftvec.s + y + 1),
				v128_load_32rep(ftvec.s + y));

			// fetch row-data for processing of two rows
			vec256_t row1 = v256_cat_128(*(const vec128_t *)(gridptr + 2*rowsel1 + 16), *(const vec128_t *)(gridptr + 2*rowsel0 + 16));
			vec256_t row0 = v256_cat_128(*(const vec128_t *)(gridptr + 2*rowsel1 + 0),  *(const vec128_t *)(gridptr + 2*rowsel0 + 0) );

			// -- beginning of calculation that's a copy-paste of the above loop --
			row0 = v256_permute_8x16_pair(row0, shuf_pat1);
			row1 = v256_permute_8x16_pair(row1, shuf_pat1);

			vec256_t a = v256_mul_add_pair_u16(vt, vs2); // 16-bit lane: bits[15:8] <- ((fs * ft) + 8) >> 4  = w11
			vec256_t c = v256_shri_u16(vt, 12);          // 16-bit lane: bits[7:0] <- ft
			vec256_t b = v256_shri_u16(a, 8);            // bits[7:0] <- w11
			c = v256_sub_16(c, b);                       // bits[7:0] <- ft - w11
			a = v256_and(a, xh_mask.v);
			a = v256_or(a, c);                          // bits[15:8] <- w11, bits[7:0] <- ft - w11

			b = v256_sub_8(vs1, a);


			// perform the multiplications for index infill, and add together the results - and add 8.
			a = v256_mul_add_pair_u16(row1, a);
			b = v256_mul_add_pair_u16(row0, b);
			a = v256_add_16(a, b);

			// right-shift multiply-result by 4 bits, put it in upper 8-bit sublane of 16-bit lane, then compute 64-x in the lower lane.
			b = v256_shri_u16(a, 5);
			a = v256_sub_8(xf_rsubv.v, b);
			a = v256_interleave_low_u16(b, a);

			// -- copy-paste end --
			// write back the result in two parts.
			// Use 32-byte writes.
			res[(y>>1)+12] = a;
		}
	}
}

/**************************************************************************
  ASTC-3d index infill routine

  What we take as input are:
	* a decoded index-ISE (1 byte per index)
	* block-size (Xdim, Ydim)
	* grid-size (Xdim, Ydim, dual-plane-flag)
  What we want to produce as output is a set of vectors such that:
	* For each pixel, we have 4 bytes laid out as:
		byte 3: <second-plane-weight>
		byte 2: <64 minus second-plane-weight>
		byte 1: <first-plaen-weight>
		byte 0: <64 minus first-plane-weight>
	  bytes 3:2 don't matter for things that aren't dual-index-planes
	* We want the output to be 36 vectors, of 32 bytes each
**************************************************************************/

// Under gcc, unpacking-loop can be sped up with gcc's computed-goto extensions;
// This is not useful for clang or for any compiler that doesn't implement this extension.
#if defined(__GNUC__) && !defined(__clang__)
	#define ASTCSD_INFILL3D_USE_COMPUTED_GOTO
#endif

static void astcsd_index_infill_3d(
	const astc_simd_decode_processed_params_t* bd,
	const astcsd_block_mode_desc_t *bmd,
	const uint8_t *indexes,
	vec256_t *res
) {
	uint32_t y,z;

	uint32_t block_xdim = bd->block_xdim;
	uint32_t block_ydim = bd->block_ydim;
	uint32_t block_zdim = bd->block_zdim;

	uint32_t grid_xdim = bmd->grid_xdim;
	uint32_t grid_ydim = bmd->grid_ydim;
	uint32_t grid_zdim = bmd->grid_zdim;
	uint32_t dual_index_planes = bmd->dual_plane_en;

	// -------------------------------------------------------
	//   AVX2 constant-value vectors that are used unindexed
	// -------------------------------------------------------
	static const union { uint8_t s[32]; vec256_t v; } byterep4_swz = {{ 0,0,0,0,  1,1,1,1,  2,2,2,2,  3,3,3,3,  4,4,4,4,  5,5,5,5,  6,6,6,6,  7,7,7,7 }};
	static const union { uint16_t s[16]; vec256_t v; } vh64        = {{ 64,64,64,64, 64,64,64,64, 64,64,64,64, 64,64,64,64 }};
	static const union { uint8_t s[32]; vec256_t v; } v239         = {{ 239,239,239,239, 239,239,239,239,  239,239,239,239, 239,239,239,239,   239,239,239,239, 239,239,239,239,  239,239,239,239, 239,239,239,239 }};
	static const union { uint32_t s[8]; vec256_t v; } vmask16   = {{ 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff }};
	static const union { uint32_t s[8]; vec256_t v; } vmask_b0  = {{       0xff,       0xff,       0xff,       0xff,        0xff,       0xff,       0xff,       0xff }};
	static const union { uint32_t s[8]; vec256_t v; } vmask_b1  = {{     0xff00,     0xff00,     0xff00,     0xff00,      0xff00,     0xff00,     0xff00,     0xff00 }};
	static const union { uint32_t s[8]; vec256_t v; } vmask_b2  = {{   0xff0000,   0xff0000,   0xff0000,   0xff0000,    0xff0000,   0xff0000,   0xff0000,   0xff0000 }};
	static const union { uint32_t s[8]; vec256_t v; } vmask_b3  = {{ 0xff000000, 0xff000000, 0xff000000, 0xff000000,  0xff000000, 0xff000000, 0xff000000, 0xff000000 }};
	static const union { uint32_t s[8]; vec256_t v; } vmask_b12 = {{   0xffff00,   0xffff00,   0xffff00,   0xffff00,    0xffff00,   0xffff00,   0xffff00,   0xffff00 }};
	static const union { uint8_t s[32]; vec256_t v; } w01_swz   = {{ 1,3,1,3, 5,7,5,7, 9,11,9,11, 13,15,13,15,  1,3,1,3, 5,7,5,7, 9,11,9,11, 13,15,13,15 }};
	static const union { uint8_t s[32]; vec256_t v; } w23_swz   = {{ 0,3,0,3, 4,7,4,7, 8,11,8,11, 12,15,12,15,  0,3,0,3, 4,7,4,7, 8,11,8,11, 12,15,12,15 }};
	static const union { uint8_t s[32]; vec256_t v; } w45_swz   = {{ 1,2,1,2, 5,6,5,6, 9,10,9,10, 13,14,13,14,  1,2,1,2, 5,6,5,6, 9,10,9,10, 13,14,13,14 }};
	static const union { uint8_t s[32]; vec256_t v; } w67_swz   = {{ 2,0,2,0, 6,4,6,4, 10,8,10,8, 14,12,14,12,  2,0,2,0, 6,4,6,4, 10,8,10,8, 14,12,14,12 }};
	static const union { uint8_t s[32]; vec256_t v; } permvec   = {{ 0,1,2,3, 0,1,2,3, 0,1,2,3, 0,1,2,3,  4,5,6,7, 4,5,6,7, 4,5,6,7, 4,5,6,7 }};

	// ------------------------------------------------------
	//   Indexable lookup tables. These are folded into one
	//   big data structure in order to limit the number
	//   of LEA instructions needed.
	// ------------------------------------------------------
	static const struct infill3d_constants_t {
		uint32_t tosv_addends[2];
		uint32_t jtro_arr[20];
		uint64_t fo_arr[20];
		uint64_t tosv_arr[40];
		} icp = {


	// texel-index addends for each byte of a 4-byte lane:
	//    * no-dual-planes: 0,1,0,1
	//    * dual-planes:    0,2,1,3
	//static const uint32_t tosv_addends[2] =
		{ 0x01000100, 0x03010200 },

	// texel-offset vectors for T and R axes
	//static const uint32_t jtro_arr[20] =
		{
		0x02402400, 0x02002000, 0x02490000, 0x02480000,
		0x06886880, 0x04404400, 0x06d22400, 0x04912000,
		0xffffffff, 0x06886880, 0x0b1b4400, 0x08da2400,
		0xffffffff, 0xffffffff, 0x0fac6880, 0x0b234400,
		0xffffffff, 0xffffffff, 0xffffffff, 0x0fac6880,
		},

	// texel fractional-offset vectors
	//static const uint64_t fo_arr[20] =
		{
		0x0800080008000800ull, 0x000b0500000b0500ull, 0x0c0804000c080400ull, 0x0603000d0a060300ull,
		0x0000000000000000ull, 0x00050b0000050b00ull, 0x0800080008000800ull, 0x0d06000a030d0600ull,
		0xffffffffffffffffull, 0x0000000000000000ull, 0x04080c0004080c00ull, 0x030a00060d030a00ull,
		0xffffffffffffffffull, 0xffffffffffffffffull, 0x0000000000000000ull, 0x0a0d0003060a0d00ull,
		0xffffffffffffffffull, 0xffffffffffffffffull, 0xffffffffffffffffull, 0x0000000000000000ull,
		},

	// texel-offset vectors for S axis
	//static const uint64_t tosv_arr[40] =
		{
		0x0101000001010000ull, 0x0202000002020000ull, 0x0100000001000000ull, 0x0200000002000000ull, 0x0101010100000000ull, 0x0202020200000000ull, 0x0101010000000000ull, 0x0202020000000000ull,
		0x0302010003020100ull, 0x0604020006040200ull, 0x0201000002010000ull, 0x0402000004020000ull, 0x0303020201010000ull, 0x0606040402020000ull, 0x0202020101000000ull, 0x0404040202000000ull,
		0xffffffffffffffffull, 0xffffffffffffffffull, 0x0302010003020100ull, 0x0604020006040200ull, 0x0504030302010000ull, 0x0a08060604020000ull, 0x0403030201010000ull, 0x0806060402020000ull,
		0xffffffffffffffffull, 0xffffffffffffffffull, 0xffffffffffffffffull, 0xffffffffffffffffull, 0x0706050403020100ull, 0x0e0c0a0806040200ull, 0x0504040302010000ull, 0x0a08080604020000ull,
		0xffffffffffffffffull, 0xffffffffffffffffull, 0xffffffffffffffffull, 0xffffffffffffffffull, 0xffffffffffffffffull, 0xffffffffffffffffull, 0x0706050403020100ull, 0x0e0c0a0806040200ull,
		},
	};

	// --------------------
	//   end of constants
	// --------------------

	const struct infill3d_constants_t *ic = &icp;

	// prepare a collection of values that will be helpful during the index infill process.
	uint32_t so_idx = block_xdim-3 + 4*(grid_xdim-2);
	uint32_t to_idx = block_ydim-3 + 4*(grid_ydim-2);
	uint32_t ro_idx = block_zdim-3 + 4*(grid_zdim-2);

	uint32_t jto = ic->jtro_arr[ to_idx ];
	uint32_t jro = ic->jtro_arr[ ro_idx ];

	// --------------------------------------------------------------------------------
	//   prepare vectors of useful per-axis data for the processing of the S,T,R axes
	// --------------------------------------------------------------------------------
	vec256_t tosv;  // texel-offsets vector for S-axis; used to pick individual weights within a row
	vec256_t fsv;  // fractional-coordinates vector for S-axis: array of 8 items, each consisting of byte replicated 4 times
	vec256_t ftv;  // fractional-coordinates vector for T-axis: array of 8 items, each consisting of byte replicated 4 times
	vec256_t frv;  // fractional-coordinates vector for R-axis: packed array of 8 bytes, replicated 4 times

	vec256_t tosv_addends_ext = v256_load_32rep(ic->tosv_addends + dual_index_planes);

	tosv = v256_load_64rep(ic->tosv_arr + 2*so_idx + dual_index_planes);  // 64->256 broadcast
	tosv = v256_permute_8x16_pair(tosv, byterep4_swz.v); // replicate each byte 4 times
	tosv = v256_add_8(tosv, tosv_addends_ext);

	fsv = v256_load_64rep(ic->fo_arr + so_idx);  // 64->256 broadcast
	fsv = v256_permute_8x16_pair(fsv, byterep4_swz.v); // replicate each byte 4 times

	ftv = v256_load_64rep(ic->fo_arr + to_idx);  // 64->256 broadcast
	ftv = v256_permute_8x16_pair(ftv, byterep4_swz.v); // replicate each byte 4 times

	frv = v256_load_64rep(ic->fo_arr + ro_idx);  // 64->256 broadcast

	// --------------------------------------------------------------
	//   Unpack index-grid into constant-strided rows and layers,
	//   (row-stride: 16 bytes, layer-stride: 96 bytes)
	//   in order to minimize the number of misaligned accesses done
	//   and to simplify address calculations later on.
	// --------------------------------------------------------------

	// buffer to contain rows of the index-grid.
	// 6+1 layers, each with 6 rows of 12 bytes (+4 bytes padding).
	vec256_t grid_rows[22];

	vec256_t *gd_destx = grid_rows;

	uintptr_t xstride = grid_xdim << dual_index_planes;
	uintptr_t ystride = xstride * grid_ydim;

	vec256_t zeroes = v256_zero();
	const uint8_t *ix_ptr = indexes;

	#ifdef ASTCSD_INFILL3D_USE_COMPUTED_GOTO
		// The unpacking code runs a switch-statement in a loop; this
		// statement will look up a jump-target-address for each iteration.
		// The target-address is the same for every loop iteration; as such,
		// we can hoist the target-address-lookup out of the innerloop
		// using gcc's computed-goto feature. (This is only useful
		// under gcc; while clang supports computed-goto, it will actually
		// perform the desired target-address-lookup hoisting optimization
		// on its own, making the use of computed-goto redundant.)
		static const void *labels[8] = { &&LBL2, &&LBL2, &&LBL2, &&LBL3, &&LBL4, &&LBL5, &&LBL6, &&LBL7 };
		const void *label = labels[grid_ydim & 7];
	#endif

	for (z=0;z<grid_zdim;z++)
	{
		vec128_t *gd_dest = (vec128_t *)gd_destx;
		#ifdef ASTCSD_INFILL3D_USE_COMPUTED_GOTO
			goto *label;
			LBL7:
			LBL6: gd_dest[5] = v128_load_unaligned(ix_ptr + 5u*xstride);
			LBL5: gd_dest[4] = v128_load_unaligned(ix_ptr + 4u*xstride);
			LBL4: gd_dest[3] = v128_load_unaligned(ix_ptr + 3u*xstride);
			LBL3: gd_dest[2] = v128_load_unaligned(ix_ptr + 2u*xstride);
			LBL2: gd_dest[1] = v128_load_unaligned(ix_ptr + 1u*xstride);
				  gd_dest[0] = v128_load_unaligned(ix_ptr);
		#else
			switch(grid_ydim & 7)
			{
				case 7: ASTCSD_FALLTHROUGH;
				case 6: gd_dest[5] = v128_load_unaligned(ix_ptr + 5u*xstride); ASTCSD_FALLTHROUGH;
				case 5: gd_dest[4] = v128_load_unaligned(ix_ptr + 4u*xstride); ASTCSD_FALLTHROUGH;
				case 4: gd_dest[3] = v128_load_unaligned(ix_ptr + 3u*xstride); ASTCSD_FALLTHROUGH;
				case 3: gd_dest[2] = v128_load_unaligned(ix_ptr + 2u*xstride); ASTCSD_FALLTHROUGH;
				case 2: gd_dest[1] = v128_load_unaligned(ix_ptr + 1u*xstride); ASTCSD_FALLTHROUGH;
				case 1: gd_dest[0] = v128_load_unaligned(ix_ptr );             ASTCSD_FALLTHROUGH;
				case 0: ;
			}
		#endif

		// Pad a row with zeroes so that later stages will not read uninitialized memory.
		gd_dest[grid_ydim] = v128_zero();
		gd_destx += 3;
		ix_ptr += ystride;
	}


	// Pad an additional layer with zeroes so that later stages will not read uninitialized memory.
	gd_destx[0] = zeroes;
	gd_destx[1] = zeroes;
	gd_destx[2] = zeroes;
	gd_destx[3] = zeroes;


	// -----------------------------------------------------------------------
	// For ASTC-3D, the weight interpolation for weight infill can be phrased
	// as an interpolation between 8 items, taken as a grid of 2x2x2 items.
	// The weight calculation can then be reformulated from the ton of
	// conditionals listed in the official specification to:
	//   w0 = 16 - MAX(fs, MAX(ft,fr));
	//   w1 = MAX(fs-MAX(ft,fr), 0);
	//   w2 = MAX(ft-MAX(fs,fr), 0);
	//   w3 = MAX(MIN(fs,ft)-fr, 0);
	//   w4 = MAX(fr-MAX(fs,ft), 0);
	//   w5 = MAX(MIN(fs,fr)-ft, 0);
	//   w6 = MAX(MIN(ft,fr)-fs, 0);
	//   w7 = MIN(fs, MIN(ft,fr));
	// This can be reformulated further to:
	//   w0 = 16 - MAX(fs, MAX(ft,fr));
	//   w1 = MIN(MAX(fs-ft,0), MAX(fs-fr,0));
	//   w2 = MIN(MAX(ft-fs,0), MAX(ft-fr,0));
	//   w3 = MIN(MAX(fs-fr,0), MAX(ft-fr,0));
	//   w4 = MIN(MAX(fr-fs,0), MAX(fr-ft,0));
	//   w5 = MIN(MAX(fs-ft,0), MAX(fr-ft,0));
	//   w6 = MIN(MAX(ft-fs,0), MAX(fr-fs,0));
	//   w7 = MIN(fs, MIN(ft,fr));
	// which, at least for weights w1 to w6, allows a very efficient
	// vectorized calculation using saturating-unsigned-subtract
	// and minimum-value instructions.
	// -----------------------------------------------------------------------

	for (z=0;z<block_zdim;z++)
	{
		// prepare two vectors, where each 4-byte lane contains items for one row:
		// "ftrx[0]":
		//   * byte 0: ft
		//   * byte 1: fr
		//   * byte 2: MIN(ft,fr)
		//   * byte 3: MAX(ft,fr)
		// "ftrx[1]":
		//   * bytes 0 and 3: MAX(ft-fr, 0)
		//   * bytes 1 and 2: MAX(fr-ft, 0)
		union {uint32_t s[16]; vec256_t v[2]; } ftrx;

		// copy the bottom byte of 'frx' to the bottom byte of every 4-byte lane,
		// then right-shift 'frx' 1 byte to prepare it for the next pass.
		 vec256_t frx = v256_permute_8x16_pair(frv, v256_zero());
		frv = v256_shrb_pair(frv, 1);


		vec256_t ftmfr = v256_andnot(v256_satsub_u8(ftv, frx), vmask_b12.v);
		vec256_t frmft = v256_and(   vmask_b12.v, v256_satsub_u8(frx, ftv));

		vec256_t ftrmin = v256_and(  vmask_b2.v, v256_min_u8(ftv, frx));
		vec256_t ftrmax = v256_and(  vmask_b3.v, v256_max_u8(frx, ftv));

		ftrx.v[0] = v256_or(
			v256_or(v256_and(ftv, vmask_b0.v), v256_and(frx, vmask_b1.v)),
			v256_or(ftrmin, ftrmax));

		ftrx.v[1] = v256_or(ftmfr, frmft);

		// pointer to the first of the two index-layers to use
		// for infill for the current result-layer
		uint32_t jr = jro & 0x70;
		jro >>= 3;
		const uint8_t *layer_ptr = (const uint8_t *)grid_rows + 6 * jr;

		// load a vector of offsets to use for rows.
		uint32_t jtol = jto;

		if (block_xdim <= 4)
		{

			// if block-x-dimension is 4 or less, then we will process 2 rows at a time.
			// This appears to take about 62 instructions per loop iteration, while the 1-row-at-a-time version takes
			// about 50 instructions per iteration.
			for (y=0;y<block_ydim;y+=2)
			{
				// pointer to the first of the index-rows to use
				// for infill for the current result-layer
				const uint8_t *row_ptr1 = layer_ptr + (jtol & 0x70);
				jtol >>= 3;
				const uint8_t *row_ptr2 = layer_ptr + (jtol & 0x70);
				jtol >>= 3;

				// read and unpack index-grid contents
				vec256_t row00 = v256_cat_128(*(const vec128_t *)(row_ptr2),     *(const vec128_t *)(row_ptr1));
				vec256_t row01 = v256_cat_128(*(const vec128_t *)(row_ptr2+16),  *(const vec128_t *)(row_ptr1+16));
				vec256_t row10 = v256_cat_128(*(const vec128_t *)(row_ptr2+96),  *(const vec128_t *)(row_ptr1+96));
				vec256_t row11 = v256_cat_128(*(const vec128_t *)(row_ptr2+112), *(const vec128_t *)(row_ptr1+112));

				vec256_t ftrxv  = v256_load_64rep((uint64_t *)(ftrx.s + y));
				vec256_t ftrmmv = v256_load_64rep((uint64_t *)(ftrx.s + y + 8));

				ftrxv = v256_permute_8x16_pair(ftrxv, permvec.v);
				ftrmmv = v256_permute_8x16_pair(ftrmmv, permvec.v);

				vec256_t *res_ptr = res + 6*z + y;

				// ---- compute weights 2,3,4,5 -----
				vec256_t dv0 = v256_satsub_u8(ftrxv, fsv); // saturating-ussigned-subtracts
				vec256_t dv1 = v256_satsub_u8(fsv, ftrxv);

				// w2453: bytes 0-3 in each 4-byte lane: weights 2,4,5,3
				vec256_t w2453 = v256_min_u8(
					ftrmmv,
					v256_or(
						v256_shli_u32(dv1, 16),
						v256_and(dv0, vmask16.v)));

				// ---- compute weights 0,1,6,7 ----
				vec256_t dv4 = v256_andnot(vmask_b3.v, v256_add_8(v256_max_u8(fsv, ftrxv), v239.v)); // w0: use byte 3 of each 4-byte lane
				vec256_t dv5 = v256_and(v256_min_u8(fsv, ftrxv), vmask_b2.v);  // w7: use byte 2 of each 4-byte lane

				dv1 = v256_and(dv1, vmask_b3.v); // w1: byte 3
				dv0 = v256_and(dv0, vmask_b2.v); // w6: byte 2

				// w7061: bytes 0-3 in each 4-byte lane: weights 7,0,6,1
				vec256_t w7061 = v256_or(
					v256_or(dv0, dv1),
					v256_shri_u32(v256_or(dv4, dv5), 16));

				// ---- perform the interpolation ----
				vec256_t px0 = v256_mul_add_pair_u16(v256_permute_8x16_pair(row00, tosv), v256_permute_8x16_pair(w7061, w01_swz.v));  // apply weights 0,1
				vec256_t px1 = v256_mul_add_pair_u16(v256_permute_8x16_pair(row01, tosv), v256_permute_8x16_pair(w2453, w23_swz.v));  // apply weights 2,3
				vec256_t px2 = v256_mul_add_pair_u16(v256_permute_8x16_pair(row10, tosv), v256_permute_8x16_pair(w2453, w45_swz.v));  // apply weights 4,5
				vec256_t px3 = v256_mul_add_pair_u16(v256_permute_8x16_pair(row11, tosv), v256_permute_8x16_pair(w7061, w67_swz.v));  // apply weights 6,7

				vec256_t summa = v256_add_16(v256_add_16(px0, px1), v256_add_16(px2, px3));

				// right-shift by 5, then, for each 2-byte lane, put final weight
				// in upper byte and 64-minus-weight in the lower byte.
				summa = v256_shri_u16(summa, 5);
				vec256_t sl = v256_sub_16(vh64.v, summa);
				sl = v256_interleave_low_u16(summa, sl);

				// perform two 32-byte writes. We don't care about the contents of the top 16 bytes.
				res_ptr[0] = sl;
				res_ptr[1] = v256_set_128rep(v128_from_high_half_of(sl));
			}
		}
		else
		{
			for (y=0;y<block_ydim;y++)
			{
				// pointer to the first of the index-rows to use
				// for infill for the current result-layer
				const uint8_t *row_ptr = layer_ptr + (jtol & 0x70);
				jtol >>= 3;

				// read and unpack index-grid contents
				vec256_t row00 = v256_set_128rep(*(const vec128_t *)(row_ptr));
				vec256_t row01 = v256_set_128rep(*(const vec128_t *)(row_ptr + 16));
				vec256_t row10 = v256_set_128rep(*(const vec128_t *)(row_ptr + 96));
				vec256_t row11 = v256_set_128rep(*(const vec128_t *)(row_ptr + 112));

				vec256_t ftrxv  = v256_load_32rep(ftrx.s + y);
				vec256_t ftrmmv = v256_load_32rep(ftrx.s + y + 8);

				vec256_t *res_ptr = res + 6*z + y;

				// ---- compute weights 2,3,4,5 -----
				vec256_t dv0 = v256_satsub_u8(ftrxv, fsv); // saturating-ussigned-subtracts
				vec256_t dv1 = v256_satsub_u8(fsv, ftrxv);

				// w2453: bytes 0-3 in each 4-byte lane: weights 2,4,5,3
				vec256_t w2453 = v256_min_u8(
					ftrmmv,
					v256_or(
						v256_shli_u32(dv1, 16),
						v256_and(dv0, vmask16.v)));

				// ---- compute weights 0,1,6,7 ----
				vec256_t dv4 = v256_andnot(vmask_b3.v, v256_add_8(v256_max_u8(fsv, ftrxv), v239.v)); // w0: use byte 3 of each 4-byte lane
				vec256_t dv5 = v256_and(v256_min_u8(fsv, ftrxv), vmask_b2.v);  // w7: use byte 2 of each 4-byte lane

				dv1 = v256_and(dv1, vmask_b3.v); // w1: byte 3
				dv0 = v256_and(dv0, vmask_b2.v); // w6: byte 2

				// w7061: bytes 0-3 in each 4-byte lane: weights 7,0,6,1
				vec256_t w7061 = v256_or(
					v256_or(dv0, dv1),
					v256_shri_u32(v256_or(dv4, dv5), 16));

				// ---- perform the interpolation ----
				vec256_t px0 = v256_mul_add_pair_u16(v256_permute_8x16_pair(row00, tosv), v256_permute_8x16_pair(w7061, w01_swz.v));  // apply weights 0,1
				vec256_t px1 = v256_mul_add_pair_u16(v256_permute_8x16_pair(row01, tosv), v256_permute_8x16_pair(w2453, w23_swz.v));  // apply weights 2,3
				vec256_t px2 = v256_mul_add_pair_u16(v256_permute_8x16_pair(row10, tosv), v256_permute_8x16_pair(w2453, w45_swz.v));  // apply weights 4,5
				vec256_t px3 = v256_mul_add_pair_u16(v256_permute_8x16_pair(row11, tosv), v256_permute_8x16_pair(w7061, w67_swz.v));  // apply weights 6,7

				vec256_t summa = v256_add_16(v256_add_16(px0, px1), v256_add_16(px2, px3));

				// right-shift by 5, then, for each 2-byte lane, put final weight
				// in upper byte and 64-minus-weight in the lower byte.
				summa = v256_shri_u16(summa, 5);
				vec256_t sl = v256_sub_16(vh64.v, summa);
				*res_ptr = v256_interleave_low_u16(summa, sl);
			}
		}
	}
}

/*****************************************************************************
  Optimized color endpoint decode functions for ASTC.

  * 'astcsd_unpack_color_endpoints8' unpacks to a packed 8*8-bit representation.
  * 'astcsd_unpack_color_endpoints16' unpacks to a packed 8*16-bit representation.

  Input is a 64-bit integer containing a packed array of 8 consecutive 8-bit
  integers from the ASTC color-ISE representation. (For formats that require
  less than 8 integers, the high-order bytes of this 64-bit integer are
  treated as don't-care.)

  The result of 'astcsd_unpack_color_endpoints8' is interleaved:
	 bits  7: 0 : endpoint 0 red
	 bits 15: 8 : endpoint 1 red
	 bits 23:16 : endpoint 0 green
	 bits 31:24 : endpoint 1 green
	 bits 39:32 : endpoint 0 blue
	 bits 47:40 : endpoint 1 blue
	 bits 55:48 : endpoint 0 alpha
	 bits 63:56 : endpoint 1 alpha

  The result of 'astcsd_unpack_color_endpoints16' is not interleaved:
	 bits  15:  0 : endpoint 0 red
	 bits  31: 16 : endpoint 0 green
	 bits  47: 32 : endpoint 0 blue
	 bits  63: 48 : endpoint 0 alpha
	 bits  79: 64 : endpoint 1 red
	 bits  95: 80 : endpoint 1 green
	 bits 111: 96 : endpoint 1 blue
	 bits 127:112 : endpoint 1 alpha

  In the case of LDR endpoints with 16-bit decode, they're decoded
  in a special way:
	cccccccc -> 0ccc_cccc_c100_0001
  This way, the lowest bit can be used to detect that this is an LDR
  color endpoint, that needs to be postprocessed as LDR. (The HDR
  endpoint modes all set this bit to 0.)
*****************************************************************************/
static uint64_t astcsd_unpack_color_endpoints8(
	uint64_t inp,
	uint32_t fmt
) {
	switch(fmt & 0xf)
	{
		case 2:
		case 3:
		case 7:
		case 11:
		case 14:
		case 15:
			return 0xffffffff0000ffffull; // magenta.

		case 0: // LDR luminance
		{
			uint64_t v0 = inp & 0xFFFF;
			return (v0 * 0x100010001ull) + 0xFFFF000000000000ull;
		}

		case 1: // LDR luminance-delta
		{
			uint32_t l0 = ((inp >> 2) & 0x3F) | ((inp >> 8) & 0xC0);
			uint32_t l1 = l0 + ((inp >> 8) & 0x3F);
			if (l1 > 255) l1 = 255;
			uint64_t v2 = (l1 << 8) | l0;
			return (v2 * 0x100010001ull) + 0xFFFF000000000000ull;
		}

		case 4: // LDR luminance-and-alpha
		{
			uint64_t v0 = inp & 0xFFFF;
			return v0 + (v0 << 16) + (inp << 32);
		}

		case 5: // LDR luminance-and-alpha with delta
		{
			uint32_t la0 = ((inp >> 1) & 0x7F007F) | ((inp >> 8) & 0x800080);
			uint32_t la1 = ((inp >> 9) & 0x3F003F) ^ 0x200020;
			la1 += la0;
			// do a quick check to see if clamping is needed; skip it if it isn't.
			uint32_t la2 = la1 + 0x1e001e0;
			if (la2 & 0x1000100)
			{
				uint32_t lum1   = la1 & 0x1FF;
				uint32_t alpha1 = (la1 >> 16) & 0x1FF;
				if (lum1   < 32) lum1   = 32; else if (lum1   > 287) lum1   = 287;
				if (alpha1 < 32) alpha1 = 32; else if (alpha1 > 287) alpha1 = 287;
				la1 = lum1 | (alpha1 << 16);
			}
			la1 -= 0x200020;
			uint64_t v0 = la0 + (la1 << 8);
			uint64_t v1 = v0 & 0xFFFF;
			return (v0 << 32) | (v1 << 16) | v1;
		}

		case 6: // LDR RGB-and-scale
			inp |= 0xffff00000000ull;
			ASTCSD_FALLTHROUGH;

		case 10: // LDR RGBA-and-scale
		{
			#ifdef ASTCSD_FAST_PEXT_PDEP
				uint64_t v0 = u64_bit_deposit(inp, 0xff00ff00ff00ull);
			#else
				uint64_t v0 = ((inp & 0xFF) << 8) | ((inp & 0xFF00) << 16) | ((inp & 0xFF0000) << 24);
			#endif
			uint64_t scale = (uint32_t)inp >> 24;
			uint64_t v1 = v0 * scale;
			return ((v1 & 0xFF00FF00FF0000ull) >> 16) + v0 + ((inp << 16) & 0xffff000000000000ull);
		}

		case 8: // LDR RGB
			inp |= 0xffff000000000000ull;
			ASTCSD_FALLTHROUGH;

		case 12: // LDR RGBA
		{
			uint64_t mask = 0xFF00FF00FF;
			uint64_t multval = 0x0001000100010001;

			uint64_t vl0 = inp & mask;
			uint64_t vl1 = (inp >> 8) & mask;

			uint64_t sum0 = (vl0 * multval) & 0xffff000000000000ull;
			uint64_t sum1 = (vl1 * multval);
			if (sum0 > sum1)
			{
				// blue-contraction
				vl0 += (vl0 >> 32) * multval;
				vl1 += (vl1 >> 32) * multval;
				vl0 >>= 1;
				vl1 >>= 1;
				vl0 &= mask;
				vl1 &= mask;
				mask = 0x00FF000000000000ull;
				vl0 |= inp & mask;
				vl1 |= (inp >> 8) & mask;
				return (vl0 << 8) | vl1;
			}
			return inp;
		}

		case 13:  // LDR RGBA-delta
		{
			uint64_t a0 = ((inp >> 49) & 0x7F) | ((inp >> 56) & 0x80);
			uint64_t a1 = (((inp >> 57) & 0x3F) ^ 0x20) - 0x20;
			a1 += a0;
			if (a1 & 0x100) // check if clamping is needed; if so, do it.
			{
				if (a1 & 0x200) a1 = 0; else a1 = 255;
			}

			inp = (inp & 0xffffffffffffull) | (a0 << 48) | (a1 << 56);
		}
		goto CASE9_13_COMMON;

		case 9: // LDR RGB-delta
			inp |= 0xffff000000000000ull;
		CASE9_13_COMMON:
			{
			uint64_t rgb0 = ((inp >> 1) & 0x7F007F007Full) | ((inp >> 8) & 0x8000800080ull);
			uint64_t rgb1 = (inp ^ 0x400040004000ull) & 0x7e007e007e00ull;
			uint64_t rgbsum = (rgb1 * (0x0001000100010000ull>>9)) >> 48; // top 16 bits as actual sum
			rgb1 = rgb0 + (rgb1 >> 9);
			int shamt = 0;
			if (rgbsum < 0x60)
			{
				// perform blue-contraction
				static const volatile uint64_t mulv = 0x100010001ull;
				uint64_t mulv2 = mulv;
				rgb0 += (rgb0 >> 32) * mulv2;
				rgb0 >>= 1;
				rgb0 &= 0x1FF01FF01FFull;

				rgb1 += (rgb1 >> 32) * mulv2;
				rgb1 >>= 1;
				rgb1 &= 0x1FF01FF01FF;
				shamt = 8;
			}

			// do a quick check to see if any of the components of 'r1' need
			// to be clamped. This clamping can usually be avoided.
			uint64_t rgb1p = rgb1 + 0x1e001e001e0ull;
			if (rgb1p & 0x10001000100ull)
			{
				uint64_t r1 = rgb1 & 0x1FF;
				uint64_t g1 = (rgb1 >> 16) & 0x1FF;
				uint64_t b1 = (rgb1 >> 32) & 0x1FF;
				if (r1 < 32) r1 = 32; else if (r1 > 287) r1 = 287;
				if (g1 < 32) g1 = 32; else if (g1 > 287) g1 = 287;
				if (b1 < 32) b1 = 32; else if (b1 > 287) b1 = 287;
				rgb1 = r1 | (g1 << 16) | (b1 << 32);
			}

			rgb1 -= 0x2000200020ull;
			rgb0 += inp & 0x00ff000000000000ull;
			rgb1 += (inp >> 8) & 0x00ff000000000000ull;

			return (rgb0 << shamt) + (rgb1 << (8-shamt));
		}
	}
	return 0;
}

static vec128_t astcsd_unpack_color_endpoints16(
	uint64_t inp,
	uint32_t fmt
) {
	// early-out for the LDR endpoint formats.
	if ((1<<fmt) & 0x3773)
	{
		static const union { int8_t s[16]; vec128_t v; } ldr16_swizzle = {{ 0,-1, 2,-1, 4,-1, 6,-1,  1,-1, 3,-1, 5,-1, 7,-1 }};
		static const union { uint16_t s[8]; vec128_t v; } ldr_addend = {{ 65,65,65,65, 65,65,65,65 }};
		uint64_t endpts8 = astcsd_unpack_color_endpoints8(inp, fmt);

		vec128_t vl1 = v128_permute_8x16(v128_set_64zx(endpts8), ldr16_swizzle.v);
		vl1 = v128_shl_u16(vl1, 7);
		vl1 = v128_or(vl1, ldr_addend.v);

		return vl1;
	}

	// constants and variables used by more than one path
	static const union { uint16_t s[8]; vec128_t v; } hdr_alpha1 = {{ 0,0,0, 0x7800, 0,0,0, 0x7800 }};
	static const union { uint16_t s[8]; vec128_t v; } hdr_magenta = {{ 0xffff,0, 0xffff, 0xffff, 0xffff,0,0xffff, 0xffff }};
	vec128_t hdr_alpha;

	switch(fmt & 0xF)
	{
		case 0:
		case 1:
		case 4:
		case 5:
		case 6:
		case 8:
		case 9:
		case 10:
		case 12:
		case 13:
			return hdr_magenta.v;

		case 2: // HDR luminance, wide range
		{
			// --------------------
			//   useful constants
			// --------------------
			static const struct mode2_constants_t {
				union { int8_t s[16]; vec128_t v; } swizzles[2];
				union { int16_t s[8]; vec128_t v; } addends[2]; }
			m2c = {
				{
				{{ -1, 0, -1, 0, -1, 0, -1,-1,  -1, 1, -1, 1, -1, 1, -1, -1}},
				{{ -1, 1, -1, 1, -1, 1, -1,-1,  -1, 0, -1, 0, -1, 0, -1, -1}},
				},
				{
				{{ 0,0,0,0x7800,0,0,0,0x7800 }},
				{{ 128,128,128,0x7800, -128,-128,-128,0x7800}},
				} };
			// --------------------
			//   end of constants
			// --------------------
			uint32_t v0 = inp & 0xFF;
			uint32_t v1 = (inp >> 8) & 0xFF;
			uint32_t ofs = ((v1-v0) >> 8) & 0x10;
			const struct mode2_constants_t *m2v = &m2c;

			vec128_t swizzle = *(const vec128_t *) ((const uint8_t *)m2v->swizzles + ofs);
			vec128_t addend = *(const vec128_t *) ((const uint8_t *)m2v->addends + ofs);

			return v128_add_16(v128_permute_8x16(v128_set_32zx(inp), swizzle), addend);
		}

		case 3: // HDR luminance, narrow range
		{

			static const union { int8_t s[16]; vec128_t v; } m3_swizzle = {{ 2,3,2,3,2,3,4,5, 0,1,0,1,0,1,4,5 }};


			uint32_t v1 = (inp >> 8) & 0xFF;
			uint32_t shamt = ((inp >> 7) & 1) + 5;
			uint32_t v0m = ((inp >> 3) & 0x10) | 0xF;

			uint32_t y0 = ((v1 & ~v0m) << 8) | ((inp & 0x7F) << shamt);
			uint32_t y1 = (v1 & v0m) << shamt;
			y1 += y0;

			uint64_t vlx = (uint64_t)(y0 | 0x78000000) << 16;

			if (y1 > 0xFFF0)
				y1 = 0xFFF0;

			vlx |= y1;

			return v128_permute_8x16(v128_set_64zx(vlx), m3_swizzle.v);
		}


		case 7: // HDR RGB-and-scale
		{

			static const struct mode7_constants_t {
				uint32_t rha_tab[16];
				uint32_t rhb_tab[16];
				uint64_t inp_mask_tab[16];
				union { int8_t s[16]; vec128_t v; } minuend_swizzles[4];
				union { int8_t s[16]; vec128_t v; } subtrahend_swizzles[4]; }
			m7c = {

			// bits 10:6 : selection bits for red ('ra')
			// bits 5:4 : major color component (0,1,2); use value 3to indicate mode #5.
			//static const uint32_t rha_tab[16] =
				{ 0x780,0x080,0x3c0,0x080,0x790,0x090,0x3d0,0x090,0x7a0,0x0a0,0x3e0,0x0a0,0x000,0x010,0x020,0x030 },
			// bits 10:6 : selection bits for red ('rb')
			// bits 3:0 : shift-amount.
			//static const uint32_t rhb_tab[16] =
				{ 0x045,0x745,0x006,0x147,0x045,0x745,0x006,0x147,0x045,0x745,0x006,0x147,0x0c8,0x0c8,0x0c8,0x049 },

			//static const uint64_t inp_mask_tab[16] =
				{
				0x007f001f001f003full,
				0x001f003f003f003full,
				0x00ff001f001f003full,
				0x007f003f003f003full,
				0x007f001f001f003full,
				0x001f003f003f003full,
				0x00ff001f001f003full,
				0x007f003f003f003full,
				0x007f001f001f003full,
				0x001f003f003f003full,
				0x00ff001f001f003full,
				0x007f003f003f003full,
				0x003f007f007f003full,
				0x003f007f007f003full,
				0x003f007f007f003full,
				0x007f007f007f003full },

			//static const union { int8_t s[16]; vec128_t v; } minuend_swizzles[4] = {
				{
				{ { 0, 1, 0, 1, 0, 1,-1,-1,  0, 1, 0, 1, 0, 1,-1,-1} },
				{ { 0, 1, 0, 1, 0, 1,-1,-1,  0, 1, 0, 1, 0, 1,-1,-1} },
				{ { 0, 1, 0, 1, 0, 1,-1,-1,  0, 1, 0, 1, 0, 1,-1,-1} },
				{ { 0, 1, 4, 5, 2, 3,-1,-1,  0, 1, 4, 5, 2, 3,-1,-1} } },

			//static const union { int8_t s[16]; vec128_t v; } subtrahend_swizzles[4] = {
				{
				{ {-1,-1, 4, 5, 2, 3,-1,-1, -1,-1, 4, 5, 2, 3,-1,-1} },
				{ { 4, 5,-1,-1, 2, 3,-1,-1,  4, 5,-1,-1, 2, 3,-1,-1} },
				{ { 2, 3, 4, 5,-1,-1,-1,-1,  2, 3, 4, 5,-1,-1,-1,-1} },
				{ {-1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1} } }
			};

			static const union { int8_t s[16]; vec128_t v; } scale_swizzle = {{6,7,6,7,6,7,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1}};

			const struct mode7_constants_t *m7p = &m7c;

			// red: bits 6:0
			// blue: bits 23:16
			// green: bits 39:32
			// scale: bits 55:48
			uint64_t inpx = ((uint32_t)inp) | ((inp & 0xffffff00) << 24);

			// extract color components, including variable-placement bits
			#ifdef ASTCSD_FAST_PEXT_PDEP
				uintptr_t modeval = u64_bit_extract(inp, 0x8080c0);
				uint32_t ra = u64_bit_extract(inpx, 0x206000600000) << 6;
				uint32_t rb = ((inp >> 25) & 0x40) | (u64_bit_extract(inpx, 0x60004040000000) << 7);
			#else
				uintptr_t modeval = ((inp & 0xc0) >> 6) | ((inp & 0x8000) >> 13) | ((inp & 0x800000) >> 20);
				uint32_t ra =
						  ((inp >> 15) & 0xc0)      // bits 21,22
						| ((inp >> 5) & 0x300)      // bits 14,13
						| ((inp >> 11) & 0x400);    // bit 21
				uint32_t rb = ((inp >> 25) & 0x40)  // bit 31
						| ((inp >> 23) & 0x80)      // bit 30
						| ((inp >> 6) & 0x100)      // bit 14
						| ((inp >> 20) & 0x600);    // bits 30,29
			#endif

			uint32_t rha_val = m7p->rha_tab[modeval];
			uint32_t rhb_val = m7p->rhb_tab[modeval];

			inpx &= m7p->inp_mask_tab[modeval];
			inpx |= (ra & rha_val) | (rb & rhb_val); // apply variable-placement for red-component
			inpx <<= rhb_val & 0xF;

			vec128_t inpy = v128_set_64zx(inpx);

			// For modes 0-4, color components 1 and 2 are deltas from component 0, so perform
			// a subtraction. Additionally, bake in a major-color-component reordering.

			uint32_t majcomp = rha_val & 0x30;
			vec128_t minuend_swizzle    = *(const vec128_t *)((const uint8_t *)m7p->minuend_swizzles + majcomp);
			vec128_t subtrahend_swizzle = *(const vec128_t *)((const uint8_t *)m7p->subtrahend_swizzles + majcomp);

			vec128_t subtrahend =  v128_permute_8x16(inpy, subtrahend_swizzle);
			vec128_t minuend = v128_or(v128_permute_8x16(inpy, minuend_swizzle), hdr_alpha1.v);

			// prepare subtraction by scale-value
			vec128_t scalev = v128_permute_8x16(inpy, scale_swizzle.v);

			// perform subtractions and clamping
			return v128_satsub_u16(v128_satsub_u16(minuend, subtrahend), scalev);
		}


		case 15: // HDR RGB + HDR Alpha
		{
			uint8_t mode = ((inp >> 55) & 1) | ((inp >> 62) & 2);

			uint32_t v6 = (inp >> 48) & 0x7f;
			uint32_t v7 = (inp >> 55) & 0xfe;

			if (mode == 3)
			{
				v6 <<= 25;
				v7 <<= 8;
			}
			else
			{
				v6 |= (v7 << mode) & 0x780;
				v6 += v6;
				v7 &= 0x7f >> mode;
				v7 ^= 0x40 >> mode;
				v7 -= 0x40 >> mode;
				v7 += v6;

				v7 <<= (mode ^ 7);
				v6 <<= (mode ^ 7);
				if (v7 & 0x10000)
					v7 = (v7 & 0x20000) ? 0 : 0xfff0;
				v6 <<= 16;
			}
			uint32_t va = v6 | v7;
			static const union { int8_t s[16]; vec128_t v; } m15_alpha_shuffle = {{ -1,-1,-1,-1,-1,-1, 2,3, -1,-1,-1,-1,-1,-1, 0,1 }};
			hdr_alpha = v128_permute_8x16(v128_set_32zx(va), m15_alpha_shuffle.v);
			goto MODE_11_14_15_COMMON;
		}

		case 14: // HDR RGB + LDR Alpha
		{
			static const union { int8_t s[16]; vec128_t v; } m14_alpha_shuffle = {{ -1,-1,-1,-1,-1,-1, 6,-1,  -1,-1,-1,-1,-1,-1, 7,-1 }};
			static const union { uint16_t s[16]; vec128_t v; } m14_alpha_addend = {{ 0,0,0,65, 0,0,0,65 }};
			hdr_alpha = v128_permute_8x16(v128_set_64zx(inp), m14_alpha_shuffle.v);
			hdr_alpha = v128_shl_u16(hdr_alpha, 7);
			hdr_alpha = v128_or(hdr_alpha, m14_alpha_addend.v);

			goto MODE_11_14_15_COMMON;
		}

		case 11: // HDR RGB
			hdr_alpha = hdr_alpha1.v;

			MODE_11_14_15_COMMON:
			{

			// -------------------------------------------------------
			//   constants that will be useful for decoding mode #11
			// -------------------------------------------------------

			static const struct mode11_constants_t {
				uint8_t shift_tabs[8][4];
				union { int8_t s[16]; vec128_t v; } xov_shuffles[3];
				union { int8_t s[16]; vec128_t v; } b_shuffles[3];
				union { int8_t s[16]; vec128_t v; } d_shuffles[3];
				uint64_t exov_tab[8];
				uint64_t inp_mask_tab[8];  }
			m11c = {
			// shift-values to apply to each of the four masked-bits to move it into the
			// correct position. Bits that are to be thrown away get a shift-amount of 63.
			//static const uint8_t shift_tabs[8][4] =
				{
				{ 63,63,63,63 },
				{ 63,25,63,19 },
				{ 49,63,12,63 },
				{ 63, 8,63,45 },
				{ 46,25,45,19 },
				{ 49, 8,48,11 },
				{ 46, 8,45,47 },
				{ 49, 8,48,47 },
				},

			//static const union { int8_t s[16]; vec128_t v; } xov_shuffles[3] =
				{
				{{ -1,-1, 7,-1, 7,-1, -1,-1,  -1,-1,-1,-1,-1,-1,-1,-1 }},
				{{  7,-1,-1,-1, 7,-1, -1,-1,  -1,-1,-1,-1,-1,-1,-1,-1 }},
				{{  7,-1, 7,-1,-1,-1, -1,-1,  -1,-1,-1,-1,-1,-1,-1,-1 }} },
			//static const union { int8_t s[16]; vec128_t v; } b_shuffles[3]   =
				{
				{{ -1,-1, 2,-1, 3,-1, -1,-1,  -1,-1, 2,-1, 3,-1,-1,-1 }},
				{{  2,-1,-1,-1, 3,-1, -1,-1,   2,-1,-1,-1, 3,-1,-1,-1 }},
				{{  3,-1, 2,-1,-1,-1, -1,-1,   3,-1, 2,-1,-1,-1,-1,-1 }} },
			//static const union { int8_t s[16]; vec128_t v; } d_shuffles[3]   =
				{
				{{ -1,-1, 4,-1, 5,-1, -1,-1,  -1,-1,-1,-1,-1,-1,-1,-1 }},
				{{  4,-1,-1,-1, 5,-1, -1,-1,  -1,-1,-1,-1,-1,-1,-1,-1 }},
				{{  5,-1, 4,-1,-1,-1, -1,-1,  -1,-1,-1,-1,-1,-1,-1,-1 }} },

			//static const uint64_t exov_tab[8] =
				{
				0x4000404000000000ull,
				0x2000202000000000ull,
				0x4000404000000000ull,
				0x2000202000000000ull,
				0x1000101000000000ull,
				0x2000202000000000ull,
				0x1000101000000000ull,
				0x2000202000000000ull },

			//static const uint64_t inp_mask_tab[8] =
				{
				0x7f7f7f7f3fffull,
				0x3f3f7f7f3fffull,
				0x7f7f3f3f3fffull,
				0x3f3f7f7f3fffull,
				0x1f1f7f7f3fffull,
				0x3f3f3f3f3fffull,
				0x1f1f7f7f3fffull,
				0x3f3f3f3f3fffull } };


			static const union { int8_t s[16];  vec128_t v; } a_shuffle = {{  0, 6, 0, 6, 0, 6, -1,-1,   0, 6, 0, 6, 0, 6,-1,-1 }};
			static const union { int8_t s[16];  vec128_t v; } c_shuffle = {{  1,-1, 1,-1, 1,-1, -1,-1,  -1,-1,-1,-1,-1,-1,-1,-1 }};
			static const union { uint16_t s[8]; vec128_t v; } hclamp    = {{ 0xfff, 0xfff, 0xfff, 0, 0xfff, 0xfff, 0xfff, 0 }};

			static const union { int8_t s[16];  vec128_t v; } maj3_swizzle = {{ -1, 0, -1, 2, -1, 4, -1, -1,  -1, 1,-1, 3,-1, 5, -1, -1 }};
			static const union { uint16_t s[8]; vec128_t v; } maj3_mask    = {{ 0,0,0xffff, 0, 0,0,0, 0 }};

			// --------------------
			//   end of constants
			// --------------------


			// extract the mode-value bitfield.
			// this is bits [47,39,31,23,15] of the input vector
			// (the top bit of each byte, except the first byte)

			#ifdef ASTCSD_FAST_PEXT_PDEP
				uint8_t modeval = (uint8_t)u64_bit_extract(inp, 0x808080808000);
			#else
				uint8_t modeval = (uint8_t)(((inp & 0x808080808000ull) * 0x10204081ull) >> 43);
			#endif

			// majority-component bitfield (bits [47-39] of input) should usually be 0, 1 or 2; if
			// it is 3, then a special case is flagged. This case corresponds to the above 'modeval'
			// taking a value of 24 or greater.
			if (modeval >= 24)
			{
				vec128_t inpy = v128_set_64zx(inp);
				inpy = v128_add_8(inpy, v128_and(inpy, maj3_mask.v)); // left-shift byte lanes 4 and 5 by 1 bit
				inpy = v128_permute_8x16(inpy, maj3_swizzle.v); // pad each byte out to 16 bits
				return v128_or(inpy, hdr_alpha);
			}

			const struct mode11_constants_t *m11t = &m11c;


			uintptr_t majcomp = modeval & 0x18;
			modeval &= 7;


			// mask inputs, and shift a[8] into position as well.
			uint64_t inpx =
				  (inp & m11t->inp_mask_tab[modeval])
				| ((inp & 0x4000) << 34);

			// toggle the top bits of d0 and d1, also set bits 63:56 to an
			// addend-value that's needed to mimic the effects of the
			// sign-extension of d0 and d1 that the specification requires.
			inpx ^= m11t->exov_tab[modeval];

			// --------------------------------------------------------------------------
			//   Variable-position bit handling.
			//
			//   Extract from input:  bits: 46,45,38,37,30,29,22
			//   This corresponds to: bits: X3,X5,X2,X4,X1,__,X0 in the specification
			//   Then put these bits in the proper places in the modified-input vector.
			//   Only bits 6:2 and bit 0 of the 'inpp' result actually matter; other
			//   bits will be masked off later anyway.
			// --------------------------------------------------------------------------

			#ifdef ASTCSD_FAST_PEXT_PDEP
				uint64_t inp_vpb = u64_bit_extract(inp, 0x0000606060400000);
			#else
				uint64_t inp_vpb = ((inp & 0x0000606060400000) * 0x41041) >> 40;
			#endif

			// Each byte of this big constant is used to select two of the bits in 'inp_vpb'
			// that will be selected and shifted below; the two other bits of 'inp_vpb'
			// that will be masked and shifted have fixed positions.
			uint64_t mask_val = 0x0528052800050000ull >> (modeval << 3);
			mask_val &= inp_vpb;

			const uint8_t *shift_tab = m11t->shift_tabs[modeval];

			// do the actual masking, and shifting of these bits.
			//  ..clang will actually auto-vectorize this into 256-bit vector operations (!)
			inpx |=
				  ((mask_val & 9)    << shift_tab[0])
				| ((inp_vpb & 0x40)  << shift_tab[1])
				| ((mask_val & 36)   << shift_tab[2])
				| ((inp_vpb & 0x10)  << shift_tab[3]);

			// result-vector for bit-shuffles:
			//  bits  7: 0 : 'a', bits 7:0
			//  bits 15: 8 : 'c'
			//  bits 23:16 : 'b0'
			//  bits 31:24 : 'b1'
			//  bits 39:32 : 'd0', with sign toggled
			//  bits 47:32 : 'd1', with sign toggled
			//  bits 55:48 : 'a', bits 15:8
			//  bits 63:56 : 'xov': correction term for the missing sign-extension of 'd0' and 'd1'.
			vec128_t inpy = v128_set_64zx(inpx);

			// Shuffles to apply to the 'b', 'd' and 'xov' values.
			// We implicitly implement the majority-component reordering here,
			// by selecting between three shuffle-values.
			vec128_t xov_shuffle = *(const vec128_t *)((const uint8_t *)m11t->xov_shuffles + 2*majcomp);
			vec128_t b_shuffle   = *(const vec128_t *)((const uint8_t *)m11t->b_shuffles + 2*majcomp);
			vec128_t d_shuffle   = *(const vec128_t *)((const uint8_t *)m11t->d_shuffles + 2*majcomp);

			vec128_t subtrahend = v128_add_16(v128_permute_8x16(inpy, a_shuffle.v), v128_permute_8x16(inpy, xov_shuffle));
			vec128_t minuend = v128_add_16(v128_add_16(v128_permute_8x16(inpy, b_shuffle), v128_permute_8x16(inpy, c_shuffle.v)), v128_permute_8x16(inpy, d_shuffle));

			vec128_t sum = v128_satsub_u16(subtrahend, minuend); // saturating subtract, to apply low-clamp.

			sum = v128_shl_u16(sum, (modeval>>1) ^ 3);
			sum = v128_min_s16(sum, hclamp.v);   // minimum value, to apply high-clamp. This clamp must be applied before the final shift

			sum = v128_shli_u16(sum, 4); // final shift
			return v128_or(sum, hdr_alpha);  // set the alpha-channel items
		}

	}

	return hdr_magenta.v; // should not be reachable.
}

// -------------------------------------------------------------------
//   Function to convert ASTC-LNS to FP16. This is needed for
//   decompression of the HDR endpoint color formats for
//   ASTC-HDR.
// -------------------------------------------------------------------
static inline vec256_t astcsd_lns_to_fp16_vec(
	vec256_t inp
) {
	static const union { int16_t  s[16]; vec256_t v; } v449  = {{ 449,449,449,449, 449,449,449,449, 449,449,449,449, 449,449,449,449 }};
	static const union { int16_t  s[16]; vec256_t v; } v511  = {{ 511,511,511,511, 511,511,511,511, 511,511,511,511, 511,511,511,511 }};
	static const union { int16_t  s[16]; vec256_t v; } v1537 = {{ 1537,1537,1537,1537, 1537,1537,1537,1537, 1537,1537,1537,1537, 1537,1537,1537,1537 }};
	static const union { uint16_t s[16]; vec256_t v; } v2047 = {{ 2047,2047,2047,2047, 2047,2047,2047,2047, 2047,2047,2047,2047, 2047,2047,2047,2047 }};
	static const union { uint16_t s[16]; vec256_t v; } v7bff = {{ 0x7bff,0x7bff,0x7bff,0x7bff, 0x7bff,0x7bff,0x7bff,0x7bff, 0x7bff,0x7bff,0x7bff,0x7bff, 0x7bff,0x7bff,0x7bff,0x7bff }};

	vec256_t p0 = v256_and(inp, v2047.v);
	vec256_t pa = v256_satsub_u16(v511.v, p0);  // constant: 511 - 1
	vec256_t pb = v256_satsub_u16(p0, v1537.v); // constant: 1536 + 1
	vec256_t pc = v256_avg_u16(pa, pb);         // built-in add-1; don't-care if pa and pb are both zero, but we need then adjusted down by 1 otherwise.
												  // at this point, pc = mt>>1.
	vec256_t pd = v256_avg_u16(pc, v1537.v);    // right-shift and add 1. Also add 4*v, a constant we cancel out later. v < 1024.
	vec256_t pe = v256_avg_u16(pd, inp);
	vec256_t pf = v256_sub_16(pe, v449.v);      // subtract: 1 for +1s from AVG, 64 for a bias for the LNS, and the v above.
	return v256_min_s16(pf, v7bff.v);
}

// ---------------------------------------------------------
//   2D block copy; designed to make full use of vector
//   loads/stores without doing overread/overwrite.
// ---------------------------------------------------------

static void astcsd_2d_block_copy(
	void* dst,
	const void *src,
	uint32_t num_rows,   // must be greater than 0.
	uint32_t width,      // 1 unit = 4 bytes
	intptr_t src_stride,  // 1 unit = 1 byte
	intptr_t dst_stride   // 1 unit = 1 byte
) {
	uint8_t *dstp = (uint8_t *)dst;
	const uint8_t *srcp = (const uint8_t *)src;

	uint32_t rem_width = width;

	intptr_t dv0 = dst_stride * num_rows;
	dstp += dv0;
	dv0 = -dv0;

	// copy first 8 units with 256-bit moves
	while(rem_width >= 8u)
	{
		uint8_t * ASTCSD_RESTRICT d2 = dstp;
		const uint8_t * ASTCSD_RESTRICT s2 = srcp;
		intptr_t dv = dv0;
		do {
			v256_store_unaligned(d2 + dv, v256_load_unaligned(s2));
			s2 += src_stride;
			dv += dst_stride;
		} while(dv);

		srcp += 32;
		dstp += 32;
		rem_width -= 8;
	}

	// if input width is 8 or more, and there are more than 4 items left,
	// then copy remaining data with 256-bit moves.
	if (rem_width > 4u && width >= 8)
	{
		intptr_t rem_adj = 8 - rem_width;
		srcp -= rem_adj << 2;
		dstp -= rem_adj << 2;
		uint8_t * ASTCSD_RESTRICT d2 = dstp;
		const uint8_t * ASTCSD_RESTRICT s2 = srcp;
		intptr_t dv = dv0;
		do {
			v256_store_unaligned(d2 + dv, v256_load_unaligned(s2));
			s2 += src_stride;
			dv += dst_stride;
		} while(dv);
		return;
	}

	// if remaining-width is 4 or more, do 128-bit copies.
	if (rem_width >= 4u)
	{
		uint8_t * ASTCSD_RESTRICT d2 = dstp;
		const uint8_t * ASTCSD_RESTRICT s2 = srcp;
		intptr_t dv = dv0;
		do {
			v128_store_unaligned(d2 + dv, v128_load_unaligned(s2));
			s2 += src_stride;
			dv += dst_stride;
		} while(dv);
		srcp += 16;
		dstp += 16;
		rem_width -= 4;
	}

	if (rem_width & 2)
	{
		uint8_t * ASTCSD_RESTRICT d2 = dstp;
		const uint8_t * ASTCSD_RESTRICT s2 = srcp;
		intptr_t dv = dv0;

		do {
			uint64_t vl = *(const uint64_t *)s2;
			u64_store_unaligned(d2 + dv, vl);
			s2 += src_stride;
			dv += dst_stride;
		} while(dv);
		srcp += 8;
		dstp += 8;
	}

	if (rem_width & 1)
	{
		uint8_t * ASTCSD_RESTRICT d2 = dstp;
		const uint8_t * ASTCSD_RESTRICT s2 = srcp;
		intptr_t dv = dv0;

		do {
			uint32_t vl = *(const uint32_t *)s2;
			u32_store_unaligned(d2 + dv, vl);
			s2 += src_stride;
			dv += dst_stride;
		} while(dv);
	}
}

#endif // matches: #if defined(ASTC_SIMD_DECODE_IMPLEMENTATION) && !defined(ASTCSD_RECURSIVE_INCLUDE)

#if defined(ASTC_SIMD_DECODE_IMPLEMENTATION) && defined(ASTCSD_RECURSIVE_INCLUDE)

// -------------------------------------------------------------------
//   Functions that exist in specialized versions for:
//     * ASTC-2D-LDR
//     * ASTC-2D-HDR
//     * ASTC-3D-LDR
//     * ASTC-3D-HDR
//   In order to generate the specialized versions, the
//   functions are written with a bunch of preprocessor
//   directives to control the aspects of their operation
//   that differ between the specializations. The specializations
//   are then generated by having this file recursively include
//   itself, once for each of the desired specializations.
//
//   The functions that exist in such specialized versions
//   are the main block decode function and a row-processing
//   function.
// -------------------------------------------------------------------

//  ****************************************************************************
//  ****************************************************************************
//  *****                                                                  *****
//  *****               Main ASTC block decode function                    *****
//  *****                                                                  *****
//  ****************************************************************************
//  ****************************************************************************

void ASTCSD_BLOCK_DECODE_FUNC_NAME
	(
	const uint8_t* blk,                  // Input ASTC block
	const astc_simd_decode_processed_params_t* bd, // ASTC-decode parameter block
	uint32_t dst_rows,                   // Number of rows to actually output
	intptr_t dst_row_stride,             // Byte stride for each row to output
#ifdef ASTCSD_DEFINE_DECODE_3D
	uint32_t dst_layers,                 // Number of layers to actually output
	intptr_t dst_layer_stride,           // Byte stride for each layer to output
#endif
	void* dst                            // Memory address to write block to
) {
	uintptr_t i;
	#ifdef ASTCSD_DEFINE_DECODE_3D
		uintptr_t j;
	#endif
	#ifdef ASTCSD_DEFINE_DECODE_HDR
		uintptr_t k;
	#endif

	uint32_t block_hdr = *(const uint32_t *)blk;

	#ifdef ASTCSD_DEFINE_DECODE_3D
		const astcsd_block_mode_desc_t *bmd = &(astcsd_block_mode_descs[block_hdr & 0x7FF]) + 2048;
	#else
		const astcsd_block_mode_desc_t *bmd = &(astcsd_block_mode_descs[block_hdr & 0x7FF]);
	#endif
	uint32_t block_type = bmd->block_type; // block-type: 0=normal, 1=error, 2=LDR-void-extent, 3=HDR-void-extent

	// -------------------------------------------------
	//  ASTC void-extent/constant-color block handling
	// -------------------------------------------------
	if (block_type != 0)
	{
		if (block_type >= 2)
		{
			// possible void-extent; check whether it is actually valid.
			uint64_t block_low = *(const uint64_t *)blk;
			if ((block_low & 0xfffffffffffffc00ull) != 0xfffffffffffffc00ull)
			{
				#ifdef ASTCSD_DEFINE_DECODE_3D
					uint32_t vxt_min_s = (block_low >> 10) & 0x1ff;
					uint32_t vxt_max_s = (block_low >> 19) & 0x1ff;
					uint32_t vxt_min_t = (block_low >> 28) & 0x1ff;
					uint32_t vxt_max_t = (block_low >> 37) & 0x1ff;
					uint32_t vxt_min_r = (block_low >> 46) & 0x1ff;
					uint32_t vxt_max_r = (block_low >> 55) & 0x1ff;

					if (vxt_min_s >= vxt_max_s || vxt_min_t >= vxt_max_t || vxt_min_r >= vxt_max_r)
						block_type = 1; // invalid void-extent: error-block
				#else
					uint32_t vxt_min_s = (block_low >> 12) & 0x1fff;
					uint32_t vxt_max_s = (block_low >> 25) & 0x1fff;

					uint32_t vxt_min_t = (block_low >> 38) & 0x1fff;
					uint32_t vxt_max_t = (block_low >> 51) & 0x1fff;

					// invalid void-extent gets decoded as error-block.
					if ((block_low & 0xc00) != 0xc00 || vxt_min_s >= vxt_max_s || vxt_min_t >= vxt_max_t)
						block_type = 1; // invalid void-extent: error-block
				#endif
			}
		}

		OUTPUT_CONST_COLOR: ;
		// constant-color block: error-color or void-extent

		#ifdef ASTCSD_DEFINE_DECODE_HDR
			// prepare a 128-bit vector of FP16 constant-color-block data
			static const uint32_t pixels_in_256bits = 4;
			const uint64_t *ccolor_ptr = block_type >= 2 ? (const uint64_t *)(blk+8) : &(bd->hdr_error_color);
			vec128_t ccolor_x = v128_load_64rep(ccolor_ptr);
			if (block_type == 2)
			{
				// error block or LDR constant color block: convert UNORM16 to FP16.
				ccolor_x = v128_unorm16_to_fp16z_full(ccolor_x);
			}
			else
			{
				// perform NaN-quietening if we have a HDR constant-color block.
				static const union { uint16_t s[8]; vec128_t v; } nsmask   = {{ 0x7fff,0x7fff,0x7fff,0x7fff, 0x7fff,0x7fff,0x7fff,0x7fff }};
				static const union { uint16_t s[8]; vec128_t v; } fp16infs = {{ 0x7c00,0x7c00,0x7c00,0x7c00, 0x7c00,0x7c00,0x7c00,0x7c00 }};
				static const union { uint16_t s[8]; vec128_t v; } fp16_mantmsb = {{ 0x200,0x200,0x200,0x200, 0x200,0x200,0x200,0x200 }};
				vec128_t v0 = v128_cmpgt_s16(v128_and(ccolor_x, nsmask.v), fp16infs.v); // 0xffff for NaNs, 0 otherwise.
				ccolor_x = v128_or(ccolor_x, v128_and(v0, fp16_mantmsb.v));
			}
		#else
			// prepare a 128-bit vector of UNORM8 constant-color-block data
			static const uint64_t cc_errorcolor = 0xffffffff0000ffffull;
			static const uint32_t pixels_in_256bits = 8;
			const uint64_t *ccolor_ptr = block_type == 2 ? (const uint64_t *)(blk+8) : &cc_errorcolor;
			static const union { uint8_t s[16]; vec128_t v; } cc_perm = {{ 1,3,5,7, 1,3,5,7, 1,3,5,7, 1,3,5,7 }};
			vec128_t ccolor_x = v128_permute_8x16(v128_load_64zx(ccolor_ptr), cc_perm.v);
		#endif

		// expand the 128-bit vector to 256 bits.
		vec256_t ccolor_y = v256_set_128rep(ccolor_x);
		uint8_t *dst8 = (uint8_t *)dst;


		// constant-color blocks for ASTC-3D: iterate over the layers, performing one
		// wide write per row (or two writes, in case of HDR with blocks wider than 4 texels)
		#ifdef ASTCSD_DEFINE_DECODE_3D
			ASTCSD_UNUSED(pixels_in_256bits);
			for (j=0; j<dst_layers; j++)
			{
				uint8_t *dstb8 = dst8;
				#ifdef ASTCSD_DEFINE_DECODE_HDR
				if (bd->block_xdim > 4)
				{
					for (i=0; i<dst_rows; i++)
					{
						v256_store_unaligned(dstb8, ccolor_y);
						v128_store_unaligned(dstb8+32, ccolor_x);
						dstb8 += dst_row_stride;
					}
				}
				else
				#endif
				for (i=0; i<dst_rows; i++)
				{
					v256_store_unaligned(dstb8, ccolor_y);
					dstb8 += dst_row_stride;
				}
				dst8 += dst_layer_stride;
			}
		#else
		// Constant-color blocks for ASTC-2D: perform 1, 2 or 3 writes per row.
		// The case for 3 writes occurs only for HDR. For LDR but not HDR, write #2
		// can be limited to 128 bits because max block size is 12.
			#ifdef ASTCSD_DEFINE_DECODE_HDR
			if (bd->block_xdim > 8)
			{
				for (i=0; i<dst_rows; i++)
				{
					v256_store_unaligned(dst8, ccolor_y);
					v256_store_unaligned(dst8+32, ccolor_y);
					v256_store_unaligned(dst8+64, ccolor_y);
					dst8 += dst_row_stride;
				}
			}
			else
			#endif
			if (bd->block_xdim > pixels_in_256bits)
			{
				for (i=0; i<dst_rows; i++)
				{
					v256_store_unaligned(dst8, ccolor_y);
					#ifdef ASTCSD_DEFINE_DECODE_HDR
						v256_store_unaligned(dst8+32, ccolor_y);
					#else
						v128_store_unaligned(dst8+32, ccolor_x);
					#endif
					dst8 += dst_row_stride;
				}
			}
			else
			{
				for (i=0; i<dst_rows; i++)
				{
					v256_store_unaligned(dst8, ccolor_y);
					dst8 += dst_row_stride;
				}
			}
		#endif
		// ---------------------------------------------
		//   End of ASTC constant-color block handling
		// ---------------------------------------------
		return;
		}

	// ---------------------------------------------------
	//   Apparent regular block, analyze its mode fields
	// ---------------------------------------------------
	uint32_t partcnt_m1 = (block_hdr >> 11) & 3; // partition-count minus 1.

	// Compute number of color ISE bits plus 7.
	// The number of color-ISE bits after everything else is subtracted
	// may be -7 to +87; as such, we will add 7 so that we can use the value as an array index.

	// Also extract color endpoint modes and second-index-plane color component,
	// as these are byproducts of the same calculation.
	uint32_t color_ise_bits_p7 = (111+7+9) - (uint32_t)bmd->index_ise_bits - ((uint32_t)bmd->dual_plane_en << 1);
	uint32_t spi_pos = 128 - 2 - (uint32_t)bmd->index_ise_bits;

	static const uint32_t cem_pos_tab[4] = { 9,21,21,21 };

	uint32_t cpt = cem_pos_tab[partcnt_m1];

	uint32_t cem = (block_hdr >> cpt) & 0xfc;  // CEM-bits, left-shifted by 2
	color_ise_bits_p7 -= cpt;

	// put partcnt in bottom 2 bits of the cem-left-shift-by-2 value
	cem += partcnt_m1;

	static const uint8_t cxs_tab[16] = { 0,0,0,0, 0,2,5,8, 0,2,5,8, 0,2,5,8 };
	uint32_t cem_ext_size = cxs_tab[ cem & 0xF ];

	color_ise_bits_p7 -= cem_ext_size;
	spi_pos -= cem_ext_size;

	// compute the address of the partition table data to use
	#ifdef ASTCSD_DEFINE_DECODE_3D
		static const uint32_t partidx_mask[4] = { 0, (0x3ff<<7), (0x3ff<<7), (0x3ff<<7) };
		uint32_t partidx = ((block_hdr >> (13-7)) & partidx_mask[partcnt_m1]);
	#else
		static const uint32_t partidx_mask[4] = { 0, (0x3ff<<6), (0x3ff<<6), (0x3ff<<6) };
		uint32_t partidx = ((block_hdr >> (13-6)) & partidx_mask[partcnt_m1]);
	#endif

	const vec128_t *part_ptr = (const vec128_t *)(bd->partptrs[partcnt_m1] + partidx);

	#ifdef __GNUC__
		// prefetch partition data, as we will need it later.
		__builtin_prefetch(part_ptr);
	#endif

	// fetch an item containing second-index-plane-color-component and/or additional CEM-bits.
	uint32_t spi_item = u32_load_unaligned(blk + (spi_pos >> 3)); // potentially unaligned access, confined to the input block; assumes little-endian.
	spi_item >>= (spi_pos & 7);
	uint32_t second_plane = (spi_item & 3) + ((uint32_t)bmd->dual_plane_en << 2);

	uint32_t color_ise_paircount_m1 = astcsd_cise_count_tbl[ cem ]; // number of integers in color-ISE minus 2, divided by 2 : range 0..8. 32 if error.
	uint32_t color_ise_quant_level = astcsd_cise_quant_tab[ color_ise_bits_p7 * 9 + color_ise_paircount_m1]; // quantization level for color-ISE. 32 if error.

	cem >>= 2;
	cem |= (spi_item & 0x3fc) << 4;

	// perform a few checks for invalid blocks
	if ((partcnt_m1 + bmd->dual_plane_en) & 4) goto OUTPUT_CONST_COLOR; // check for dual-planes with 4 partitions.
	if (bd->block_xdim < bmd->grid_xdim) goto OUTPUT_CONST_COLOR;
	if (bd->block_ydim < bmd->grid_ydim) goto OUTPUT_CONST_COLOR;
	#ifdef ASTCSD_DEFINE_DECODE_3D
		if (bd->block_zdim < bmd->grid_zdim) goto OUTPUT_CONST_COLOR;
	#endif
	if ((color_ise_paircount_m1 | color_ise_quant_level) & 32) goto OUTPUT_CONST_COLOR; // check for too many color-ISE-integers or too few bits for the color-ISE.

	// ------------------------------------------------------------------
	//   If we have reached this point, then we know that the block is
	//   valid and a 'regular' block (not error-block, not void-extent)
	// ------------------------------------------------------------------

	// ---------------------------------------
	//   Unpack and unquantize the two ISEs.
	// ---------------------------------------
	union { uint8_t s[80]; vec128_t x[5]; vec256_t v[2]; } index_ise;
	union { uint8_t s[32]; vec256_t v; } color_ise;

	astcsd_index_ise_unpack(
		blk,
		bmd->index_quant_level,
		bmd->index_count,
		index_ise.v);

	color_ise.v = astcsd_color_ise_unpack(
		blk,
		color_ise_quant_level,
		color_ise_paircount_m1 + 1);


	// -------------------------------------------------------
	//   obtain the color endpoint format items from the CEM
	// -------------------------------------------------------
	uint32_t endpt_fmts = ((cem >> 2) & 0xF) * 0x01010101;

	if (cem & 3)
	{
		uint32_t cemh =  ((cem & 3) - 1);
		uint32_t cem_m = cem >> (partcnt_m1 + 3);
		#ifdef ASTCSD_FAST_PEXT_PDEP
			endpt_fmts = (uint32_t) ((cemh*0x04040404) + u64_bit_deposit(cem>>2, 0x04040404) + u64_bit_deposit(cem_m, 0x03030303));
		#else
			endpt_fmts = (cemh * 0x04040404) + (((cem & 0x3c) * 0x204081) & 0x04040404) + ((((cem_m & 0x3f) * 0x1041) + (cem_m << 18)) & 0x03030303);
		#endif
	}

	// ---------------------------------
	//   perform color endpoint decode
	// ---------------------------------

	uint32_t ce_offset = 0;

#ifdef ASTCSD_DEFINE_DECODE_HDR
	// a table that, for each of the 16 endpoint modes specifies:
	//  bit 0: set if any of the endpoint color components contain LDR
	//  bit 1: set if any of the endpoint color components contain HDR
	static const uint8_t lhdr_flags[16] = { 1,1,2,2,1,1,1,2,1,1,1,2,1,1,3,2 };
	// "lhdr_mode": used to pick what type of postprocessing we
	// will need to do for decoding to FP16.
	//  * bit 0 set: LDR color components present
	//  * bit 1 set: HDR color components present
	// If both bits are set, then we need to perform a selection
	// based on per-color-component metadata. If only one bit is set,
	// then all color components in the block use the same
	// postprocessing.
	uint8_t lhdr_mode = 0;


	union { uint64_t s[8]; vec128_t x[4]; vec256_t v[2]; } decoded_endpts;
	decoded_endpts.v[0] = v256_zero();
	decoded_endpts.v[1] = v256_zero();
	for (i=0; i<=partcnt_m1; i++)
	{
		static const union { int8_t s[16]; vec128_t v; } endpt_swz   = {{ 0,8,2,10,4,12,6,14, 1,9,3,11,5,13,7,15 }};
		static const union { uint16_t s[8]; vec128_t v; } signbits16 = {{ 0x8000,0x8000,0x8000,0x8000,0x8000,0x8000,0x8000,0x8000 }};

		uint32_t endpt_fmt = endpt_fmts & 0xFF;
		lhdr_mode |= lhdr_flags[endpt_fmt];
		endpt_fmts >>= 8;
		uint64_t endpt_inp = u64_load_unaligned(color_ise.s + 2*ce_offset);  // unaligned access
		vec128_t decoded_endpt = astcsd_unpack_color_endpoints16(endpt_inp, endpt_fmt);
		decoded_endpt = v128_xor(decoded_endpt, signbits16.v); // flip signs because lerp uses signed multiply.

		decoded_endpt = v128_permute_8x16(decoded_endpt, endpt_swz.v);

		v128_store_64low(&(decoded_endpts.s[i]), decoded_endpt);
		v128_store_64high(&(decoded_endpts.s[i+4]), decoded_endpt);

		ce_offset += (endpt_fmt >> 2) + 1;
	}

#else
	union { uint64_t s[4]; vec128_t x[2]; vec256_t v; } decoded_endpts;

	decoded_endpts.v = v256_zero();
	for (i=0; i<=partcnt_m1; i++)
	{
		uint32_t endpt_fmt = endpt_fmts & 0xFF;
		endpt_fmts >>= 8;
		uint64_t endpt_inp = u64_load_unaligned(color_ise.s + 2*ce_offset);  // unaligned access
		decoded_endpts.s[i] = astcsd_unpack_color_endpoints8(endpt_inp, endpt_fmt);
		ce_offset += (endpt_fmt >> 2) + 1;
	}
#endif

	// ----------------------------------
	//   Perform ASTC index-grid infill
	// ----------------------------------

#ifdef ASTCSD_DEFINE_DECODE_3D
	vec256_t infilled_indexes[36];
	astcsd_index_infill_3d(
#else
	vec256_t infilled_indexes[18];
	astcsd_index_infill_2d(
#endif
		bd,
		bmd,
		index_ise.s,
		infilled_indexes);

	// -------------------------------------------------------------------
	//   Perform data preparation before the main pixel-processing loops
	// -------------------------------------------------------------------

	// sort the bytes of the color endpoints into a form suitable for use with the final lerp routine:
	//  byte 0: partition 0, endpt 0, red
	//  byte 1: partition 1, endpt 0, red
	//  byte 2: partition 2, endpt 0, red
	//  byte 3: partition 3, endpt 0, red
	//  byte 4: partition 0, endpt 1, red
	// byte 5-7: partition 1-3, endpt1, red
	// byte 8-11: partition 0-3, endpt0, blue
	// byte 12-15: partition 0-3, endpt1, blue

	static const union { int8_t s[16]; vec128_t v; } rb01sel = {{ 0, 8,-1,-1, 1, 9,-1,-1,  4,12,-1,-1, 5,13,-1,-1 }};
	static const union { int8_t s[16]; vec128_t v; } rb23sel = {{-1,-1, 0, 8,-1,-1, 1, 9, -1,-1, 4,12,-1,-1, 5,13 }};

	static const union { int8_t s[16]; vec128_t v; } ga01sel = {{ 2,10,-1,-1, 3,11,-1,-1,  6,14,-1,-1, 7,15,-1,-1 }};
	static const union { int8_t s[16]; vec128_t v; } ga23sel = {{-1,-1, 2,10,-1,-1, 3,11, -1,-1, 6,14,-1,-1, 7,15 }};

	vec256_t epv0 = v256_set_128rep(decoded_endpts.x[0]);  // color endpoints 0 and 1
	vec256_t endpts_rb = v256_permute_8x16_pair(epv0, v256_set_128rep(rb01sel.v));
	vec256_t endpts_ga = v256_permute_8x16_pair(epv0, v256_set_128rep(ga01sel.v));
	#ifdef ASTCSD_DEFINE_DECODE_HDR
		epv0 = v256_set_128rep(decoded_endpts.x[2]);  // color endpoints 0 and 1, high byte
		vec256_t endpts_rbh = v256_permute_8x16_pair(epv0, v256_set_128rep(rb01sel.v));
		vec256_t endpts_gah = v256_permute_8x16_pair(epv0, v256_set_128rep(ga01sel.v));
	#endif

	if (partcnt_m1 >= 2)
	{
		vec256_t epv1 = v256_set_128rep(decoded_endpts.x[1]);  // color endpoints 2 and 3
		endpts_rb = v256_or(endpts_rb, v256_permute_8x16_pair(epv1, v256_set_128rep(rb23sel.v)));
		endpts_ga = v256_or(endpts_ga, v256_permute_8x16_pair(epv1, v256_set_128rep(ga23sel.v)));
		#ifdef ASTCSD_DEFINE_DECODE_HDR
			epv1 = v256_set_128rep(decoded_endpts.x[3]);  // color endpoints 2 and 3, high byte
			endpts_rbh = v256_or(endpts_rbh, v256_permute_8x16_pair(epv1, v256_set_128rep(rb23sel.v)));
			endpts_gah = v256_or(endpts_gah, v256_permute_8x16_pair(epv1, v256_set_128rep(ga23sel.v)));
		#endif
	}

	// generate index-select swizzles based on which color component, if any,
	// is assigned to the second index plane; these are applied to the result of index-infill.
	static const uint32_t sw_sels[9] = { 0,0,0,0,0, 0x202, 0, 0x2020000, 0 };
	static const union { uint8_t s[16]; vec128_t v; } ixs_orv = {{ 0,1,0,1, 4,5,4,5, 8,9,8,9, 12,13,12,13 }};

	vec256_t ixs_rb =  v256_load_32rep(sw_sels + second_plane+1);
	vec256_t ixs_ga =  v256_load_32rep(sw_sels + second_plane);
	ixs_rb = v256_or(ixs_rb, v256_set_128rep(ixs_orv.v));
	ixs_ga = v256_or(ixs_ga, v256_set_128rep(ixs_orv.v));

	// item that's ORed onto the data items that are read out of the partition-tables, to form indexes used to select into endpts_rb/ga.
	static const uint8_t pd_orvs[4] = { 0,4,8,12 };
	vec256_t pd_orv = v256_load_32rep((const uint32_t *)pd_orvs);

	static const union { uint8_t b[16]; vec128_t v; } row_swizzles[8] = {
		{{  0, 0, 0, 0,  1, 1, 1, 1,  2, 2, 2, 2,  3, 3, 3, 3 }},
		{{  0, 0, 0, 0,  1, 1, 1, 1,  2, 2, 2, 2,  3, 3, 3, 3 }},
		{{  4, 4, 4, 4,  5, 5, 5, 5,  6, 6, 6, 6,  7, 7, 7, 7 }},
		{{  4, 4, 4, 4,  5, 5, 5, 5,  6, 6, 6, 6,  7, 7, 7, 7 }},

		{{  8, 8, 8, 8,  9, 9, 9, 9, 10,10,10,10, 11,11,11,11 }},
		{{  8, 8, 8, 8,  9, 9, 9, 9, 10,10,10,10, 11,11,11,11 }},
		{{ 12,12,12,12, 13,13,13,13, 14,14,14,14, 15,15,15,15 }},
		{{ 12,12,12,12, 13,13,13,13, 14,14,14,14, 15,15,15,15 }},
		};


	// -------------------------------
	//   Main pixel processing loops
	// -------------------------------

#ifdef ASTCSD_DEFINE_DECODE_HDR
	// Metadata stash for use in case of mixed-mode HDR postprocess.
	// Since mixed-mode isn't common, our main priority is to just ensure
	// that data get into the stash with as little overhead as possible.
	// We can sort out the data later if needed.

	// layout for ASTC-2D: (4 entries per row, last 2 entries used for even rows only)
	//   * entry 0: row 0, RB (first 8 pixels)
	//   * entry 1: row 0, GA (first 8 pixels)
	//   * entry 2: row 0/1: RB (up to 8 pixels)
	//   * entry 3: row 0/1: GA (up to 8 pixels)

	// layout for ASTC-3D: (2 entries per row)
	//   * entry 0: row 0, RB (first 8 pixels)
	//   * entry 1: row 0, GA (first 8 pixels)

	#ifdef ASTCSD_DEFINE_DECODE_3D
		vec256_t mdata[2*6*6];
		vec256_t cdata[2*6*6];
	#else
		vec256_t mdata[4*12];
		vec256_t cdata[3*12];
	#endif
#endif

#ifdef ASTCSD_DEFINE_DECODE_3D
	// for ASTC-3D, we run a loop that's basically the same as the upper 8x8 block
	// of the ASTC-2D loop; this loop is run once for each layer.
	#ifndef ASTCSD_DEFINE_DECODE_HDR
		uint8_t *dst8 = (uint8_t *)dst;
	#endif

	const uint8_t *part_ptrx = (const uint8_t *)part_ptr;

	for (j=0; j<dst_layers; j++)
	{
		#ifdef ASTCSD_DEFINE_DECODE_HDR
			vec256_t *dst16b = cdata + 12*j;
			vec256_t *mdp = mdata + 12*j;
		#else
			uint8_t *dst8b = dst8;
		#endif
		vec128_t v0;
		vec256_t v1;
		vec256_t row_vl;
		vec256_t row_dif;
		v0 = v128_load_unaligned(part_ptrx);
		part_ptrx += 12;

		v1 = v256_cat_128(v128_shri_u32(v0, 2), v0);  // duplicate 128-bit item into two 128-bit lanes, with upper-lane left-shifted by 2.
		row_dif = v256_and3_8(v256_shri_u16(v1, 4));  // for odd rows: right-shift by 4, XOR with even-row item
		row_vl = v256_and3_8(v1);                       // for even-rows: keep bottom 2 bits
		// perform decode for the first 8x8 part of the block.
		for (i=0; i<dst_rows; i++)
		{
			v1 = v256_permute_8x16_pair(row_vl, v256_set_128rep(row_swizzles[i].v)); // partition-data: apply swizzle for current row
			row_vl = v256_xor(row_vl, row_dif);                     // partition-data: toggle between even and odd rows

			v1 = v256_or(v1, pd_orv);

			// for each of 8 pixels, pick color endpoints to use
			vec256_t rb_endpts = v256_permute_8x16_pair(endpts_rb, v1);
			vec256_t ga_endpts = v256_permute_8x16_pair(endpts_ga, v1);
			#ifdef ASTCSD_DEFINE_DECODE_HDR
				// for HDR, pick low and high bytes separately, then interleave.
				// also stash aside the low bytes for later possible postprocessing.
				mdp[2*i] = rb_endpts;     // byte order: ..BBRRBBRR
				mdp[2*i + 1] = ga_endpts; // byte order: ..AAGGAAGG

				vec256_t rbh_endpts = v256_permute_8x16_pair(endpts_rbh, v1);
				vec256_t gah_endpts = v256_permute_8x16_pair(endpts_gah, v1);
				vec256x2_t rb_endpts_ext = v256x2_interleave2_8x32(rbh_endpts, rb_endpts);
				vec256x2_t ga_endpts_ext = v256x2_interleave2_8x32(gah_endpts, ga_endpts);
			#endif

			// for each color component of 8 pixels, pick weights to use
			vec256_t weights = infilled_indexes[6*j + i];
			vec256_t rb_weights = v256_permute_8x16_pair(weights, ixs_rb);
			vec256_t ga_weights = v256_permute_8x16_pair(weights, ixs_ga);

			#ifdef ASTCSD_DEFINE_DECODE_HDR
				// expand weights from 8 bits to 16 bits
				vec256x2_t rb_weights_ext = v256x2_interleave2_8x32(v256_zero(), rb_weights);
				vec256x2_t ga_weights_ext = v256x2_interleave2_8x32(v256_zero(), ga_weights);
				// perform lerp and rounding
				vec256_t rb_prods0 = v256_shr6_round_u32(v256_mul_add_pair_s32(rb_endpts_ext.low,  rb_weights_ext.low ));
				vec256_t rb_prods1 = v256_shr6_round_u32(v256_mul_add_pair_s32(rb_endpts_ext.high, rb_weights_ext.high));
				vec256_t ga_prods0 = v256_shr6_round_u32(v256_mul_add_pair_s32(ga_endpts_ext.low,  ga_weights_ext.low ));
				vec256_t ga_prods1 = v256_shr6_round_u32(v256_mul_add_pair_s32(ga_endpts_ext.high, ga_weights_ext.high));
				// line up data for memory write
				vec256_t res0 = v256_interleave_low_u32(ga_prods0, rb_prods0);
				vec256_t res1 = v256_interleave_low_u32(ga_prods1, rb_prods1);
				dst16b[0] = res0;
				dst16b[1] = res1;
				dst16b += 2;
			#else
				// perform the lerp
				vec256_t rb_prods = v256_mul_add_pair_u16(rb_endpts, rb_weights);
				vec256_t ga_prods = v256_mul_add_pair_u16(ga_endpts, ga_weights);

				// perform rounding, and line up data for memory write
				rb_prods = v256_shr6_round_u16(rb_prods);
				ga_prods = v256_shr6_round_u16(ga_prods);
				vec256_t res = v256_interleave_low_u16(ga_prods, rb_prods);

				v256_store_unaligned(dst8b, res);
				dst8b += dst_row_stride;
			#endif
		}
		#ifndef ASTCSD_DEFINE_DECODE_HDR
			dst8 += dst_layer_stride;
		#endif
	}

#else

	// --------------------------------------------------------------------
	// The per-pixel processing loops for ASTC-2D.
	//
	// There are three loops:
	//   * The first one covers the upper-left 8x8 texels
	//     of each block, and is used for all block sizes.
	//   * The second one covers the right 4x12 texels of
	//     each block, and is used for block sizes with
	//     a width of more than 8 texels.
	//   * The third one covers the bottom-left 8x4 texels
	//     of each block, and is used for block sizes with a height
	//     of more than 8 texels.
	//
	// For each of the three loops, the code will fetch one 128-bit
	// block of precomputed packed-partition-data; this block
	// is prepared differently for each of the three loops.
	// --------------------------------------------------------------------

	// prepare partition-data for the first 8x8 part of the block
	vec128_t v0;
	vec256_t v1;
	vec256_t row_vl;
	vec256_t row_dif;
	v0 = part_ptr[0];
	v1 = v256_cat_128(v128_shri_u32(v0, 2), v0);  // duplicate 128-bit item into two 128-bit lanes, with upper-lane left-shifted by 2.
	row_dif = v256_and3_8(v256_shri_u16(v1, 4));  // for odd rows: right-shift by 4, XOR with even-row item
	row_vl = v256_and3_8(v1);                       // for even-rows: keep bottom 2 bits

	#ifdef ASTCSD_DEFINE_DECODE_HDR
		vec256_t *dst16 = cdata;
	#else
		uint8_t *dst8 = (uint8_t *)dst;
	#endif

	uint32_t pass1_limit = (dst_rows < bd->min_block_ydim_8) ? dst_rows : bd->min_block_ydim_8;

	// perform decode for the first 8x8 part of the block.
	for (i=0; i<pass1_limit; i++)
	{
		v1 = v256_permute_8x16_pair(row_vl, v256_set_128rep(row_swizzles[i].v)); // partition-data: apply swizzle for current row
		row_vl = v256_xor(row_vl, row_dif);                     // partition-data: toggle between even and odd rows

		v1 = v256_or(v1, pd_orv);

		// for each of 8 pixels, pick color endpoints to use
		vec256_t rb_endpts = v256_permute_8x16_pair(endpts_rb, v1);
		vec256_t ga_endpts = v256_permute_8x16_pair(endpts_ga, v1);
		#ifdef ASTCSD_DEFINE_DECODE_HDR
			// for HDR, pick low and high bytes separately, then interleave.
			// also stash aside the low bytes for later possible postprocessing.
			mdata[4*i] = rb_endpts;     // byte order: ..BBRRBBRR
			mdata[4*i + 1] = ga_endpts; // byte order: ..AAGGAAGG

			vec256_t rbh_endpts = v256_permute_8x16_pair(endpts_rbh, v1);
			vec256_t gah_endpts = v256_permute_8x16_pair(endpts_gah, v1);
			vec256x2_t rb_endpts_ext = v256x2_interleave2_8x32(rbh_endpts, rb_endpts);
			vec256x2_t ga_endpts_ext = v256x2_interleave2_8x32(gah_endpts, ga_endpts);
		#endif

		// for each color component of 8 pixels, pick weights to use
		vec256_t weights = infilled_indexes[i];
		vec256_t rb_weights = v256_permute_8x16_pair(weights, ixs_rb);
		vec256_t ga_weights = v256_permute_8x16_pair(weights, ixs_ga);

		#ifdef ASTCSD_DEFINE_DECODE_HDR
			// expand weights from 8 bits to 16 bits
			vec256x2_t rb_weights_ext = v256x2_interleave2_8x32(v256_zero(), rb_weights);
			vec256x2_t ga_weights_ext = v256x2_interleave2_8x32(v256_zero(), ga_weights);
			// perform lerp and rounding
			vec256_t rb_prods0 = v256_shr6_round_u32(v256_mul_add_pair_s32(rb_endpts_ext.low,  rb_weights_ext.low ));
			vec256_t rb_prods1 = v256_shr6_round_u32(v256_mul_add_pair_s32(rb_endpts_ext.high, rb_weights_ext.high));
			vec256_t ga_prods0 = v256_shr6_round_u32(v256_mul_add_pair_s32(ga_endpts_ext.low,  ga_weights_ext.low ));
			vec256_t ga_prods1 = v256_shr6_round_u32(v256_mul_add_pair_s32(ga_endpts_ext.high, ga_weights_ext.high));
			// line up data for memory write
			vec256_t res0 = v256_interleave_low_u32(ga_prods0, rb_prods0);
			vec256_t res1 = v256_interleave_low_u32(ga_prods1, rb_prods1);
			dst16[0] = res0;
			dst16[1] = res1;
			dst16 += 3;
		#else
			// perform the lerp
			vec256_t rb_prods = v256_mul_add_pair_u16(rb_endpts, rb_weights);
			vec256_t ga_prods = v256_mul_add_pair_u16(ga_endpts, ga_weights);

			// perform rounding, and line up data for memory write
			rb_prods = v256_shr6_round_u16(rb_prods);
			ga_prods = v256_shr6_round_u16(ga_prods);
			vec256_t res = v256_interleave_low_u16(ga_prods, rb_prods);

			v256_store_unaligned(dst8, res);
			dst8 += dst_row_stride;
		#endif
	}

	if (bd->block_xdim > 8)
	{
		// prepare partition-data for the 12x4 right-part of the block
		v0 = part_ptr[2];
		v1 = v256_cat_128(v128_shri_u32(v0, 2), v0);  // duplicate 128-bit item into two 128-bit lanes, with upper-lane left-shifted by 2.
		row_dif = v256_and3_8(v256_shri_u16(v1, 4));  // for odd rows: right-shift by 4, XOR with even-row item
		row_vl = v256_and3_8(v1);                       // for even-rows: keep bottom 2 bits

		// perform decode for the 12x4 right-part of the block
		i=0;
		uint32_t r=dst_rows;
		#ifdef ASTCSD_DEFINE_DECODE_HDR
			dst16 = cdata + 2;
		#else
			dst8 = (uint8_t *)dst;
			dst8 += 32;
		#endif

		while(1)
		{
			v1 = v256_permute_8x16_pair(row_vl, v256_set_128rep(row_swizzles[i].v)); // partition-data: apply swizzle for current row
			row_vl = v256_xor(row_vl, row_dif);                     // partition-data: toggle between even and odd rows

			v1 = v256_or(v1, pd_orv);

			// for each of 8 pixels, pick color endpoints to use
			vec256_t rb_endpts = v256_permute_8x16_pair(endpts_rb, v1);
			vec256_t ga_endpts = v256_permute_8x16_pair(endpts_ga, v1);
			#ifdef ASTCSD_DEFINE_DECODE_HDR
				// for HDR, pick low and high bytes separately, then interleave.
				// also stash aside the low bytes for later possible postprocessing.
				mdata[8*i+2] = rb_endpts;
				mdata[8*i+3] = ga_endpts;

				vec256_t rbh_endpts = v256_permute_8x16_pair(endpts_rbh, v1);
				vec256_t gah_endpts = v256_permute_8x16_pair(endpts_gah, v1);
				vec256x2_t rb_endpts_ext = v256x2_interleave2_8x32(rbh_endpts, rb_endpts);
				vec256x2_t ga_endpts_ext = v256x2_interleave2_8x32(gah_endpts, ga_endpts);
			#endif

			// for each color component of 8 pixels, pick weights to use
			vec256_t weights = infilled_indexes[i+12];
			vec256_t rb_weights = v256_permute_8x16_pair(weights, ixs_rb);
			vec256_t ga_weights = v256_permute_8x16_pair(weights, ixs_ga);

			#ifdef ASTCSD_DEFINE_DECODE_HDR
				// expand weights from 8 bits to 16 bits
				vec256x2_t rb_weights_ext = v256x2_interleave2_8x32(v256_zero(), rb_weights);
				vec256x2_t ga_weights_ext = v256x2_interleave2_8x32(v256_zero(), ga_weights);
				// perform lerp and rounding
				vec256_t rb_prods0 = v256_shr6_round_u32(v256_mul_add_pair_s32(rb_endpts_ext.low,  rb_weights_ext.low ));
				vec256_t rb_prods1 = v256_shr6_round_u32(v256_mul_add_pair_s32(rb_endpts_ext.high, rb_weights_ext.high));
				vec256_t ga_prods0 = v256_shr6_round_u32(v256_mul_add_pair_s32(ga_endpts_ext.low,  ga_weights_ext.low ));
				vec256_t ga_prods1 = v256_shr6_round_u32(v256_mul_add_pair_s32(ga_endpts_ext.high, ga_weights_ext.high));
				// line up data for memory write
				vec256_t res0 = v256_interleave_low_u32(ga_prods0, rb_prods0);
				vec256_t res1 = v256_interleave_low_u32(ga_prods1, rb_prods1);
				dst16[0] = res0;
				if (!--r) break;
				dst16[3] = res1;
				if (!--r) break;
				dst16 += 6;
				i++;
			#else
				// perform the lerp
				vec256_t rb_prods = v256_mul_add_pair_u16(rb_endpts, rb_weights);
				vec256_t ga_prods = v256_mul_add_pair_u16(ga_endpts, ga_weights);

				// perform rounding, and line up data for memory write
				rb_prods = v256_shr6_round_u16(rb_prods);
				ga_prods = v256_shr6_round_u16(ga_prods);
				vec256_t res = v256_interleave_low_u16(ga_prods, rb_prods);

				v128_store_unaligned(dst8, v128_from_low_half_of(res));
				if (!--r) break;
				dst8 += dst_row_stride;
				v128_store_unaligned(dst8, v128_from_high_half_of(res));
				if (!--r) break;
				dst8 += dst_row_stride;
				i++;
			#endif
		}
	}

	if (dst_rows > 8)
	{
		// perform decode for the bottom part of the block
		#ifdef ASTCSD_DEFINE_DECODE_HDR
		dst16 = cdata + 24;
		#else
		dst8 = (uint8_t *)dst;
		dst8 += dst_row_stride*8;
		#endif
		v0 = part_ptr[1];
		v1 = v256_cat_128(v128_shri_u32(v0, 2), v0); // duplicate 128-bit item into two 128-bit lanes, with upper-lane left-shifted by 2.
		row_dif = v256_and3_8(v256_shri_u16(v1, 4)); // for odd rows: right-shift by 4, XOR with even-row item
		row_vl = v256_and3_8(v1);                      // for even-rows: keep bottom 2 bits

		for (i=8; i<dst_rows; i++)
		{

			v1 = v256_permute_8x16_pair(row_vl, v256_set_128rep(row_swizzles[i-8].v)); // partition-data: apply swizzle for current row
			row_vl = v256_xor(row_vl, row_dif);                       // partition-data: toggle between even and odd rows

			v1 = v256_or(v1, pd_orv);

			// for each of 8 pixels, pick color endpoints to use
			vec256_t rb_endpts = v256_permute_8x16_pair(endpts_rb, v1);
			vec256_t ga_endpts = v256_permute_8x16_pair(endpts_ga, v1);
			#ifdef ASTCSD_DEFINE_DECODE_HDR
				// for HDR, pick low and high bytes separately, then interleave.
				// also stash aside the low bytes for later possible postprocessing.
				mdata[4*i] = rb_endpts;
				mdata[4*i+1] = ga_endpts;

				vec256_t rbh_endpts = v256_permute_8x16_pair(endpts_rbh, v1);
				vec256_t gah_endpts = v256_permute_8x16_pair(endpts_gah, v1);
				vec256x2_t rb_endpts_ext = v256x2_interleave2_8x32(rbh_endpts, rb_endpts);
				vec256x2_t ga_endpts_ext = v256x2_interleave2_8x32(gah_endpts, ga_endpts);
			#endif

			// for each color component of 8 pixels, pick weights to use
			vec256_t weights = infilled_indexes[i];
			vec256_t rb_weights = v256_permute_8x16_pair(weights, ixs_rb);
			vec256_t ga_weights = v256_permute_8x16_pair(weights, ixs_ga);

			#ifdef ASTCSD_DEFINE_DECODE_HDR
				// expand weights from 8 bits to 16 bits
				vec256x2_t rb_weights_ext = v256x2_interleave2_8x32(v256_zero(), rb_weights);
				vec256x2_t ga_weights_ext = v256x2_interleave2_8x32(v256_zero(), ga_weights);
				// perform lerp and rounding
				vec256_t rb_prods0 = v256_shr6_round_u32(v256_mul_add_pair_s32(rb_endpts_ext.low,  rb_weights_ext.low ));
				vec256_t rb_prods1 = v256_shr6_round_u32(v256_mul_add_pair_s32(rb_endpts_ext.high, rb_weights_ext.high));
				vec256_t ga_prods0 = v256_shr6_round_u32(v256_mul_add_pair_s32(ga_endpts_ext.low,  ga_weights_ext.low ));
				vec256_t ga_prods1 = v256_shr6_round_u32(v256_mul_add_pair_s32(ga_endpts_ext.high, ga_weights_ext.high));
				// line up data for memory write
				vec256_t res0 = v256_interleave_low_u32(ga_prods0, rb_prods0);
				vec256_t res1 = v256_interleave_low_u32(ga_prods1, rb_prods1);
				dst16[0] = res0;
				dst16[1] = res1;
				dst16 += 3;
			#else
				// perform the lerp
				vec256_t rb_prods = v256_mul_add_pair_u16(rb_endpts, rb_weights);
				vec256_t ga_prods = v256_mul_add_pair_u16(ga_endpts, ga_weights);

				// perform rounding, and line up data for memory write
				rb_prods = v256_shr6_round_u16(rb_prods);
				ga_prods = v256_shr6_round_u16(ga_prods);
				vec256_t res = v256_interleave_low_u16(ga_prods, rb_prods);

				v256_store_unaligned(dst8, res);
				dst8 += dst_row_stride;
			#endif
		}
	}
#endif

#ifdef ASTCSD_DEFINE_DECODE_HDR

	// -------------------------------------------------------------------
	//  After the main pixel processing loops, we are done for ASTC-LDR.
	//  However, for ASTC-HDR, we need to perform postprocessing. This
	//  postprocessing consists of converting data either from
	//  UNORM16 to FP16 or from ASTC-LNS to FP16. Depending on the color
	//  endpoint modes, the choice may be made on a per-color-component
	//  basis or a per-whole-block basis.
	//
	//  If it turns out that we need the former, then the data that
	//  were stashed aside for postprocessing become necessary.
	// -------------------------------------------------------------------

	uint32_t vstrips = (bd->block_xdim + 3) >> 2;

	static const union { uint16_t s[16]; vec256_t v; } v182 =
		{{ 0x182, 0x182, 0x182, 0x182,  0x182, 0x182, 0x182, 0x182,  0x182, 0x182, 0x182, 0x182,  0x182, 0x182, 0x182, 0x182 }};
	static const union { uint16_t s[16]; vec256_t v; } vmsb16 =
		{{ 0x8000, 0x8000, 0x8000, 0x8000,  0x8000, 0x8000, 0x8000, 0x8000,  0x8000, 0x8000, 0x8000, 0x8000,  0x8000, 0x8000, 0x8000, 0x8000 }};
	static const union { uint8_t s[32]; vec256_t v; } mx_swz =
		{{ 0,2,1,3, 4,6,5,7, 8,10,9,11, 12,14,13,15,  0,2,1,3, 4,6,5,7, 8,10,9,11, 12,14,13,15 }};

	#ifdef ASTCSD_DEFINE_DECODE_3D
		uintptr_t src16_addend = 2;
	#else
		uintptr_t src16_addend = 3;
	#endif

	switch(lhdr_mode)
	{
		case 1: // LDR: UNORM16 to FP16
			for (k=0;k<vstrips;k++)
			{
				#ifdef ASTCSD_DEFINE_DECODE_3D
				for (j=0;j<dst_layers;j++)
				{
					uint8_t *dst8 = (uint8_t *)dst + j*dst_layer_stride;
					vec256_t *src16 = cdata + j*12;
				#else
				{
					uint8_t *dst8 = (uint8_t *)dst;
					vec256_t *src16 = cdata;
				#endif
					dst8 += 32*k;
					src16 += k;
					for (i=0; i<dst_rows; i++)
					{
						vec256_t v = *src16;
						vec256_t v2 = v256_shli_u16(v, 1);
						vec256_t v3 = v256_shri_u16(v, 7);
						vec256_t v4 = v256_sub_16(v2, v182.v);
						vec256_t v5 = v256_add_16(v3, v4);
						v256_store_unaligned(dst8, v256_unorm16_to_fp16z(v5));
						dst8 += dst_row_stride;
						src16 += src16_addend;
					}
				}
			}
			break;

		case 2: // HDR: ASTC-LNS to FP16
			for (k=0;k<vstrips;k++)
			{
				#ifdef ASTCSD_DEFINE_DECODE_3D
				for (j=0;j<dst_layers;j++)
				{
					uint8_t *dst8 = (uint8_t *)dst + j*dst_layer_stride;
					vec256_t *src16 = cdata + j*12;
				#else
				{
					uint8_t *dst8 = (uint8_t *)dst;
					vec256_t *src16 = cdata;
				#endif
					dst8 += 32*k;
					src16 += k;
					for (i=0; i<dst_rows; i++)
					{
						vec256_t v = *src16;
						vec256_t v2 = v256_xor(v, vmsb16.v);
						v256_store_unaligned(dst8, astcsd_lns_to_fp16_vec(v2));
						dst8 += dst_row_stride;
						src16 += src16_addend;
					}
				}
			}
			break;

		case 3: // Mix of LDR and HDR
			{
			// ------------------------------------------------
			//   reinterleave metadata for use with selection
			// ------------------------------------------------

			#ifdef ASTCSD_DEFINE_DECODE_3D
				for (j=0; j<dst_layers; j++)
				{
					vec256_t *mdp = mdata + 12*j;
					for (i=0;i<dst_rows;i++)
					{
						vec256x2_t md2 = v256x2_interleave2_8x32(mdp[2*i+1], mdp[2*i]);
						mdp[2*i]   = v256_shli_u16(v256_permute_8x16_pair(md2.low, mx_swz.v), 7);
						mdp[2*i+1] = v256_shli_u16(v256_permute_8x16_pair(md2.high, mx_swz.v), 7);
					}
				}
			#else
				for (i=0; i<dst_rows; i++)
				{
					vec256x2_t md2 = v256x2_interleave2_8x32(mdata[4*i+1], mdata[4*i]); // byte order: ..ABABGRGR_ABABGRGR
					mdata[4*i]   = v256_shli_u16(v256_permute_8x16_pair(md2.low, mx_swz.v), 7);  // byte order: AABBGGRR_AABBGGRR
					mdata[4*i+1] = v256_shli_u16(v256_permute_8x16_pair(md2.high, mx_swz.v), 7);
				}

				if (vstrips > 2)
				{
					for (i=0;i<dst_rows; i+=2)
					{
						vec256x2_t md2 = v256x2_interleave2_8x32(mdata[4*i+3], mdata[4*i+2]);
						mdata[4*i+2] = v256_shli_u16(v256_permute_8x16_pair(md2.low, mx_swz.v), 7);
						mdata[4*i+6] = v256_shli_u16(v256_permute_8x16_pair(md2.high, mx_swz.v), 7);
					}
				}
			#endif

			// ---------------------------------------------------
			//   perform postprocessing for both LDR and HDR,
			//   and use metadata to select which one to return.
			// ---------------------------------------------------
			for (k=0;k<vstrips;k++)
			{
				#ifdef ASTCSD_DEFINE_DECODE_3D
				for (j=0; j<dst_layers; j++)
				{
					uint8_t *dst8 = (uint8_t *)dst + j*dst_layer_stride;
					vec256_t *src16 = cdata + j*12;
					vec256_t *mdp = mdata + 12*j;
				#else
				{
					uint8_t *dst8 = (uint8_t *)dst;
					vec256_t *src16 = cdata;
				#endif
					dst8 += 32*k;
					src16 += k;
					for (i=0; i<dst_rows; i++)
					{
						vec256_t v = *src16;

						vec256_t u2 = v256_xor(v, vmsb16.v);
						vec256_t u3 = astcsd_lns_to_fp16_vec(u2); // HDR result

						vec256_t v2 = v256_shli_u16(v, 1);
						vec256_t v3 = v256_shri_u16(v, 7);
						vec256_t v4 = v256_sub_16(v2, v182.v);
						vec256_t v5 = v256_add_16(v3, v4);

						vec256_t v6 = v256_unorm16_to_fp16z(v5); // LDR result

						#ifdef ASTCSD_DEFINE_DECODE_3D
							vec256_t res = v256_select8(mdp[2*i+k], v6, u3);
						#else
							vec256_t res = v256_select8(mdata[4*i+k], v6, u3);
						#endif

						v256_store_unaligned(dst8, res);
						dst8 += dst_row_stride;
						src16 += src16_addend;
					}
				}
			}
			break;
		}
	}
#endif // defined:ASTCSD_DEFINE_DECODE_HDR
}

// -----------------------------------------------------------------------------
//   Function to decode a row of ASTC blocks into a row of decoded content.
//
//   It takes into account the fact that the main ASTC block
//   decode function above might write more pixels than the block
//   size indicates - it arranges for the block decode function's
//   output to be shingled where safely possible and copied through
//   an intermediate buffer where shingling by itself would result
//   in out-of-bounds writes.
// -----------------------------------------------------------------------------
static void ASTCSD_ROW_DECODE_FUNC_NAME
	(
	uint8_t *dst,
	const uint8_t *src,
	const astc_simd_decode_processed_params_t* bd,
	uint32_t dst_rows,        // number of pixel rows to actually output
	uint32_t dst_layers       // number of pixel layers to actually output. Ignored for 2D.
) {
	uintptr_t i;

	#ifndef ASTCSD_DEFINE_DECODE_3D
		ASTCSD_UNUSED(dst_layers);
	#endif

	#ifdef ASTCSD_DEFINE_DECODE_HDR
		static const uintptr_t pixel_size = 8;
		static const uintptr_t pixel_words = 2;
		#ifdef ASTCSD_DEFINE_DECODE_3D
			static const uintptr_t tempbuf_row_stride = 8*8;
			static const uintptr_t tempbuf_layer_stride = 6*8*8;
			uint8_t temp_buf[6*6*8*8];
		#else // ASTC-2D
			static const uintptr_t tempbuf_row_stride = 12*8;
			uint8_t temp_buf[12*12*8];
		#endif
	#else // ASTC-LDR
		static const uintptr_t pixel_size = 4;
		static const uintptr_t pixel_words = 1;
		#ifdef ASTCSD_DEFINE_DECODE_3D
			static const uintptr_t tempbuf_row_stride = 8*4;
			static const uintptr_t tempbuf_layer_stride = 6*8*4;
			uint8_t temp_buf[6*6*8*4];
		#else // ASTC-2D
			static const uintptr_t tempbuf_row_stride = 12*4;
			uint8_t temp_buf[12*12*4];
		#endif
	#endif


	// Use shingled decoding for as many blocks as safely possible.
	// What this does is that the block is decoded into a width of 8 or 12
	// texels along the X axis (4 or 6 texels for 3D); the excess texels
	// are then overwritten by the next block.

	for (i=0; i<bd->shingled_blocks_per_row; i++)
	{
		const uint8_t* blk = src + 16*i;
		uint8_t *bdst = dst + pixel_size*i*bd->block_xdim;
		ASTCSD_BLOCK_DECODE_FUNC_NAME (
			blk,
			bd,
			dst_rows,
			bd->dst_row_stride,
			#ifdef ASTCSD_DEFINE_DECODE_3D
				dst_layers,
				bd->dst_layer_stride,
			#endif
			bdst);
	}

	// At the end of a texel-block row, it is not safe to do shingled
	// decoding, as this will cause a buffer overrun; as such, do
	// shingled decode into a small temp-buffer that can then be
	// copied into the main buffer.

	for (i = bd->shingled_blocks_per_row; i < bd->blocks_per_row; i++)
	{
		memset(temp_buf, 0, sizeof(temp_buf));
		const uint8_t* blk = src + 16*i;
		uint8_t *bdst = dst + (i * pixel_size * bd->block_xdim);

		ASTCSD_BLOCK_DECODE_FUNC_NAME (
			blk,
			bd,
			dst_rows,
			tempbuf_row_stride,
			#ifdef ASTCSD_DEFINE_DECODE_3D
				dst_layers,
				tempbuf_layer_stride,
			#endif
			temp_buf);

		// number of X-axis texels to copy for block
		uintptr_t rem_texels = bd->texture_xdim - i*bd->block_xdim;
		if (rem_texels > bd->block_xdim) rem_texels = bd->block_xdim;

		#ifdef ASTCSD_DEFINE_DECODE_3D
		// 3D texture: do copying for each 2D layer
		uintptr_t j;
		for (j = 0; j < dst_layers; j++)
		{
			astcsd_2d_block_copy(
				bdst     + j * bd->dst_layer_stride,
				temp_buf + j * tempbuf_layer_stride,
				dst_rows,                 // num-rows
				pixel_words * rem_texels, // width of each row
				tempbuf_row_stride,
				bd->dst_row_stride);
		}
		#else
		// 2D texture: copy a single layer
		astcsd_2d_block_copy(
			bdst,
			temp_buf,
			dst_rows,                 // num-rows
			pixel_words * rem_texels, // width of each row
			tempbuf_row_stride,
			bd->dst_row_stride);
		#endif
	}
}

#endif // matches: defined(ASTC_SIMD_DECODE_IMPLEMENTATION) && defined(ASTCSD_RECURSIVE_INCLUDE)

// --------------------------------------------------------------
//   In order to generate all the versions of the main ASTC
//   block decode function, we #include ourselves recursively,
//   with a different set of macros for each inclusion in
//   order to activate the features needed.
// --------------------------------------------------------------
#if defined(ASTC_SIMD_DECODE_IMPLEMENTATION) && !defined(ASTCSD_RECURSIVE_INCLUDE)

	#define ASTCSD_RECURSIVE_INCLUDE

	// define the ASTC-2D LDR decode routines
	#define ASTCSD_BLOCK_DECODE_FUNC_NAME astc_simd_decode_block_2d_8bit
	#define ASTCSD_ROW_DECODE_FUNC_NAME astc_simd_decode_row_2d_8bit
	#include __FILE__
	#undef ASTCSD_ROW_DECODE_FUNC_NAME
	#undef ASTCSD_BLOCK_DECODE_FUNC_NAME

	// define the ASTC-3D LDR decode routines
	#define ASTCSD_DEFINE_DECODE_3D
	#define ASTCSD_BLOCK_DECODE_FUNC_NAME astc_simd_decode_block_3d_8bit
	#define ASTCSD_ROW_DECODE_FUNC_NAME astc_simd_decode_row_3d_8bit
	#include __FILE__
	#undef ASTCSD_ROW_DECODE_FUNC_NAME
	#undef ASTCSD_BLOCK_DECODE_FUNC_NAME
	#undef ASTCSD_DEFINE_DECODE_3D

	// define the ASTC-2D HDR decode routines
	#define ASTCSD_DEFINE_DECODE_HDR
	#define ASTCSD_BLOCK_DECODE_FUNC_NAME astc_simd_decode_block_2d_16bit
	#define ASTCSD_ROW_DECODE_FUNC_NAME astc_simd_decode_row_2d_16bit
	#include __FILE__
	#undef ASTCSD_ROW_DECODE_FUNC_NAME
	#undef ASTCSD_BLOCK_DECODE_FUNC_NAME
	#undef ASTCSD_DEFINE_DECODE_HDR

	// define the ASTC-3D HDR decode routines
	#define ASTCSD_DEFINE_DECODE_3D
	#define ASTCSD_DEFINE_DECODE_HDR
	#define ASTCSD_BLOCK_DECODE_FUNC_NAME astc_simd_decode_block_3d_16bit
	#define ASTCSD_ROW_DECODE_FUNC_NAME astc_simd_decode_row_3d_16bit
	#include __FILE__
	#undef ASTCSD_ROW_DECODE_FUNC_NAME
	#undef ASTCSD_BLOCK_DECODE_FUNC_NAME
	#undef ASTCSD_DEFINE_DECODE_HDR
	#undef ASTCSD_DEFINE_DECODE_3D

	#undef ASTCSD_RECURSIVE_INCLUDE

// --------------------------------------------------------------
// Given:
//  * a texture X-dimension
//  * an ASTC block X-dimension
// compute the number of blocks that we can decode per row
// in a shingled manner. Any remaining blocks in a row must
// be decoded into a separate buffer, then copied from that
// separate buffer into the main result buffer.
// --------------------------------------------------------------

static void astcsd_compute_shingled_blocks_per_row(
	uint32_t texture_xdim,
	uint32_t astc_block_xdim,
	uint32_t *shingled_blocks,
	uint32_t *total_blocks
) {
	uint32_t num_blocks = (texture_xdim + astc_block_xdim - 1) / astc_block_xdim;
	*total_blocks = num_blocks;

	uint32_t shingle_size = astc_block_xdim > 8 ? 12 : 8;
	while(num_blocks > 0)
	{
		uint32_t shingled_width = (num_blocks - 1)*astc_block_xdim + shingle_size;
		if (shingled_width <= texture_xdim) break;
		num_blocks--;
	}
	*shingled_blocks = num_blocks;
}


// ----------------------------------------------
//   Prepare a decoder-arguments data structure
// ----------------------------------------------

uint32_t astc_simd_decode_prepare(
	astc_simd_decode_processed_params_t *processed_params,
	const astc_simd_decode_params_t *inp_params
) {
	if (astcsd_precompute_completed == 0)
		astc_simd_decode_precompute();

	const astc_simd_decode_params_t *inp = inp_params;
	astc_simd_decode_processed_params_t *res = processed_params;

	res->block_xdim = inp->block_xdim;
	res->block_ydim = inp->block_ydim;
	res->block_zdim = inp->block_zdim;

	res->min_block_ydim_8  = inp->block_ydim < 8 ? inp->block_ydim : 8;
	res->block_xpydim      = inp->block_xdim + inp->block_ydim;
	res->block_xdim_m3_x16 = (inp->block_xdim-3)*16;
	res->block_ydim_m3_x16 = (inp->block_ydim-3)*16;

	int is_large_block = (inp->block_xdim * inp->block_ydim * inp->block_zdim) > 31;

	const uint8_t *part_ptr0 = (const uint8_t *)astcsd_partition_data;
	if (inp->block_zdim > 1)
	{
		// ASTC-3D partition table pointers
		part_ptr0 += 196736;
		if (is_large_block) part_ptr0 += 36;
		res->partptrs[0] = part_ptr0;
		res->partptrs[1] = part_ptr0 + 128;
		res->partptrs[2] = part_ptr0 + 131200;
		res->partptrs[3] = part_ptr0 + 262272;
	}
	else
	{
		// ASTC-2D partition table pointers
		if (is_large_block) part_ptr0 += 16;
		res->partptrs[0] = part_ptr0;
		res->partptrs[1] = part_ptr0 + 64;
		res->partptrs[2] = part_ptr0 + 65600;
		res->partptrs[3] = part_ptr0 + 131136;
	}

	astcsd_compute_shingled_blocks_per_row(
		inp->xres,
		inp->block_xdim,
		&(res->shingled_blocks_per_row),
		&(res->blocks_per_row));


	res->src_start_addr = (const uint8_t *)(inp->src_start_ptr);
	res->dst_start_addr = (uint8_t *)(inp->dst_start_ptr);

	res->src_row_stride   = inp->src_row_stride;
	res->src_layer_stride = inp->src_layer_stride;

	res->dst_row_stride   = inp->dst_row_stride;
	res->dst_layer_stride = inp->dst_layer_stride;

	res->texture_xdim = inp->xres;
	res->texture_ydim = inp->yres;
	res->texture_zdim = inp->zres;

	res->decode_row_func = (inp->block_zdim > 1)
		? ((inp->decode_flags & 1) ? astc_simd_decode_row_3d_16bit : astc_simd_decode_row_3d_8bit)
		: ((inp->decode_flags & 1) ? astc_simd_decode_row_2d_16bit : astc_simd_decode_row_2d_8bit);

	res->hdr_error_color = (inp->decode_flags & 2) ? 0x3C003C0000003C00ULL : 0xFFFFFFFFFFFFFFFFULL;

	res->block_rows_per_layer = (inp->yres + inp->block_ydim - 1) / inp->block_ydim;
	res->block_layers = (inp->zres + inp->block_zdim - 1) / inp->block_zdim;
	return res->block_layers * res->block_rows_per_layer;
}

void astc_simd_decode_row_iterate(
	const astc_simd_decode_processed_params_t *args,
	uint32_t row_index
) {

	uintptr_t layer_index = 0;
	intptr_t src_layer_offset = 0;
	intptr_t dst_layer_offset = 0;

	if (row_index >= args->block_rows_per_layer)
	{
		layer_index  = row_index / args->block_rows_per_layer;
		row_index    = row_index % args->block_rows_per_layer;
		src_layer_offset = layer_index * args->src_layer_stride;
		dst_layer_offset = layer_index * args->block_zdim * args->dst_layer_stride;
	}

	intptr_t src_row_offset = row_index * args->src_row_stride;
	intptr_t dst_row_offset = row_index * args->block_ydim * args->dst_row_stride;

	uint32_t dst_rows =
		row_index == (args->block_rows_per_layer - 1)
			? (args->texture_ydim - args->block_ydim*row_index)
			: args->block_ydim;

	uint32_t dst_layers =
		layer_index == (args->block_layers - 1)
			? (args->texture_zdim - args->block_zdim*layer_index)
			: args->block_zdim;

	args->decode_row_func(
		args->dst_start_addr + dst_row_offset + dst_layer_offset,
		args->src_start_addr + src_row_offset + src_layer_offset,
		args,
		dst_rows,
		dst_layers);
}

// --------------------------------
//   Full-texture decode function
// --------------------------------
void astc_simd_decode_full(
	const astc_simd_decode_params_t *inp_params
) {
	astc_simd_decode_processed_params_t processed_params;
	uint32_t num_rows = astc_simd_decode_prepare(&processed_params, inp_params);
	for (uint32_t i = 0; i < num_rows; i++)
	{
		astc_simd_decode_row_iterate(&processed_params, i);
	}
}

#endif // matches: #if defined(ASTC_SIMD_DECODE_IMPLEMENTATION) && !defined(ASTCSD_RECURSIVE_INCLUDE)
