/*
* Copyright (c) Hanyu Jiang 2015
*
* Self-defined PTX function for signed 16-bits SIMD
*/

#ifndef SIMD_FUNS__
#define SIMD_FUNS__

/* extract the most significant 16-bits from one 32-bits register, and
*  put it into a 32-bits integer
*/
static __device__ __forceinline__ int extr_int16(int a)
{
	int r;
	asm("{                            \n\t"
		".reg .b32 r;                 \n\t"
		"mov.b32   r, %1;             \n\t" // r = 0xXXXXXXXX
		"shr.s32   %0, r, 16;		  \n\t" // %0 = (r) = 0xffffXXXX / 0x0000XXXX (must shift signed!)
		"}"
		: "=r"(r) : "r"(a));      /* r is %0, a is %1 */
	return r;
}

/* duplicate one 16-bits value to TWO that is wrapped in one 32-bits
*  integer (one physical register). For further SIMD operation.
*/
static __device__ __forceinline__ int dup_int16(int a)
{
	int r;
	asm("{                            \n\t"
		".reg .b32 r, t;              \n\t"
		"mov.b32   r, %1;             \n\t" // r = 0xffffXXXX / 0x0000XXXX
		"shl.b32   t, r, 16;		  \n\t" // t = 0xXXXX0000
		"and.b32   r, r, 0x0000ffff;  \n\t" // r = 0x0000XXXX
		"or.b32    %0, t, r;          \n\t" // %0 = 0xXXXXXXXX
		"}"
		: "=r"(r) : "r"(a));      /* r is %0, a is %1 */
	return r;
}

/* This function is used to get
* value between
* 2 16-bits signed integers which is packed
* as one 32-bits integer. Maximum value will be
* broadcast to all of them. (only sm_30 or higher)
*
*/
static __device__ __forceinline__ int max_int16(int a)
{
	int r;

	/* [1] [2] vs. [2] [1] */
	asm("{								\n\t"
		".reg .s32 r, a, b, s;			\n\t"
		"mov.b32   r, %1;				\n\t" // r = 0xXXXX XXXX ([1] [2])
		"and.b32   s, s, 0x00000000;	\n\t" // s = 0x00000000
		"shl.b32   a, r, 16;			\n\t" // a = 0xXXXX0000
		"shr.u32   b, r, 16;			\n\t" // b = 0x0000XXXX / shift usigned so no 0xffffXXXX case though its negative
		"or.b32	   b, b, a;				\n\t" // merge, rotate 1 16-bits value, b = 0xXXXX XXXX ([2] [1])
		"vmax2.s32.s32.s32 %0, r, b, s; \n\t"
		"}"
		: "=r"(r) : "r"(a));

	return r;
}

/* Integrate "extract_uint16" and "reorder_unit16" */
static __device__ __forceinline__ int reorder_int16(int a)
{
	int r;
	asm("{                                       \n\t" // volatile or NOT
		".reg .s32 r, t;                         \n\t"
		".reg .pred p;							 \n\t" // predicate

		"mov.b32   r, %1;                        \n\t" // r = 0xXXXXXXXX
		"mov.b32   t, %1;                        \n\t" // t = 0xXXXXXXXX
		"shr.u32   r, r, 16;                     \n\t" // to be 0x0000XXXX (like max_16int, ignore signed here)

		"shfl.idx.b32 r, r, %laneid + 31, 0x1f;  \n\t" // hard-coded "width" = 32, so "c"=0x1f 正好测试一下用imm代替reg作为输入是否合法: 合法！
		"shl.b32   t, t, 16;                     \n\t" // t = 0xXXXX0000 from [0x0001 0000] to [0x0000 0000]

		"setp.eq.u32 p, %laneid + 0, 0;			 \n\t" //
		"@p	mov .u32 r, 0x00008000;				 \n\t" // -infinity is -32768 for VIT

		"or.b32    %0, t, r;                     \n\t" //merge t = 0x0000XXXX with r = 0xXXXX0000
		"}"
		: "=r"(r) : "r"(a));      /* r is %0, a is %1 */
	return r;
}


/* extract the most significant 8-bits from one 32-bits register, and
*  put it into a 32-bits integer
*/
static __device__ __forceinline__ unsigned int extr_uint8(unsigned int a)
{
	unsigned int r;
	asm("{                            \n\t"
		".reg .u32 r;                 \n\t" //why there always be .u32 ?
		"mov.b32   r, %1;             \n\t" //r = 0xXXXXXXXX
		"and.b32   %0, r, 0x000000ff; \n\t" // extract least 8-bits
		"}"
		: "=r"(r) : "r"(a));      /* r is %0, a is %1 */
	return r;
}

/* duplicate one 8-bits value to four that is wrapped in one 32-bits
*  integer (one physical register). For further SIMD operation.
*/
static __device__ __forceinline__ unsigned int dup_uint8(unsigned int a)
{
	unsigned int r;
	asm("{                            \n\t"
		".reg .u32 r, t;              \n\t"
		"mov.b32   r, %1;             \n\t" // r = 0x000000XX
		"shl.b32   t, r, 8;			  \n\t" // t = 0x0000XX00
		"or.b32    r, r, t;           \n\t" // r = 0x0000XXXX
		"shl.b32   t, r, 16;          \n\t" // t = 0xXXXX0000
		"or.b32    %0, t, r;          \n\t" // %0 = 0xXXXXXXXX
		"}"
		: "=r"(r) : "r"(a));      /* r is %0, a is %1 */
	return r;
}

/* This function is used to get
value between
four 8-bits unsigned integers which is packed
as one 32-bits integer. Maximum value will be
broadcast to all of them. (only sm_30 or higher)
*/

static __device__ __forceinline__ unsigned int max_uint8(unsigned int a)
{
	unsigned int r;

	/* 1234 vs 4123 = 4234 */
	/* 4234 vs 4423 = 4434 */
	/* 4434 vs 4443 = 4444 */
	asm("{								\n\t"
		".reg .u32 r, a, b, s;			\n\t"
		"mov.b32   r, %1;				\n\t"
		"and.b32   s, s, 0x00000000;	\n\t"
		"shl.b32   a, r, 24;			\n\t" // a = 0xXX000000
		"shr.u32   b, r, 8;				\n\t" // b = 0x00XXXXXX
		"or.b32	   b, b, a;				\n\t" // merge, rotate 1 8-bits value, b = 0xXXXXXXXX
		"vmax4.u32.u32.u32 r, r, b, s;	\n\t"
		"and.b32   s, s, 0x00000000;	\n\t"
		"shl.b32   a, r, 24;			\n\t" // a = 0xXX000000
		"shr.u32   b, r, 8;				\n\t" // b = 0x00XXXXXX
		"or.b32	   b, b, a;				\n\t" // merge, rotate 1 8-bits value, b = 0xXXXXXXXX
		"vmax4.u32.u32.u32 r, r, b, s;	\n\t"
		"and.b32   s, s, 0x00000000;	\n\t"
		"shl.b32   a, r, 24;			\n\t" // a = 0xXX000000
		"shr.u32   b, r, 8;				\n\t" // b = 0x00XXXXXX
		"or.b32	   b, b, a;				\n\t" // merge, rotate 1 8-bits value, b = 0xXXXXXXXX
		"vmax4.u32.u32.u32 %0, r, b, s;	\n\t"
		"}"
		: "=r"(r) : "r"(a));

	return r;
}

/* Integrate "extract_uint8" and "reorder_unit8" */
static __device__ __forceinline__ unsigned int reorder_uint8(unsigned int a)
{
	unsigned int r;
	asm("{                                       \n\t" // volatile or NOT
		".reg .u32    r, t;                      \n\t"
		".reg .pred   p;						 \n\t"
		"mov.b32      r, %1;                     \n\t" // r = 0xXXXXXXXX
		"mov.b32      t, %1;                     \n\t" // t = 0xXXXXXXXX
		"shr.u32      r, r, 24;                  \n\t" // to be 0x000000XX [i.e: r = 0000 0000 0000 0000 0000 0000 0000 0011 = 0x00 00 00 03]
		// from [0x03 02 01 00] to [0x00 00 00 03]
		"shfl.idx.b32 r, r, %laneid + 31, 0x1f;  \n\t" // hard-coded "width" = 32, so "c"=0x1f
		"shl.b32      t, t, 8;                   \n\t" // t = 0xXXXXXX00 from [0x03 02 01 00] to [0x02 01 00 00]

		"setp.eq.u32  p, %laneid + 0, 0;		 \n\t" // if (threadIdx.x == 0) there is a branch
		"@p   mov.u32 r, 0x00000000;			 \n\t" // then r = 0x00 00 00 00 (avoid full filled hmm size)
		// to MSV, -infinity = 0; to VIT, -infinity = -32768
		"or.b32      %0, t, r;                   \n\t" // merge t = 0xXXXXXX00 with r = 0x00000000 to be [0x02 01 00 ff]
		"}"
		: "=r"(r) : "r"(a));      /* r is %0, a is %1 */
	return r;
}

/* Integrate "extract_uint8" and "reorder_unit8" */
static __device__ __forceinline__ unsigned int reorder_uint8_ssv(unsigned int a)
{
	unsigned int r;
	asm("{                                       \n\t" // volatile or NOT
		".reg .u32    r, t;                      \n\t"
		".reg .pred   p;						 \n\t"
		"mov.b32      r, %1;                     \n\t" // r = 0xXXXXXXXX
		"mov.b32      t, %1;                     \n\t" // t = 0xXXXXXXXX
		"shr.u32      r, r, 24;                  \n\t" // to be 0x000000XX [i.e: r = 0000 0000 0000 0000 0000 0000 0000 0011 = 0x00 00 00 03]
		// from [0x03 02 01 00] to [0x00 00 00 03]
		"shfl.idx.b32 r, r, %laneid + 31, 0x1f;  \n\t" // hard-coded "width" = 32, so "c"=0x1f
		"shl.b32      t, t, 8;                   \n\t" // t = 0xXXXXXX00 from [0x03 02 01 00] to [0x02 01 00 00]

		"setp.eq.u32  p, %laneid + 0, 0;		 \n\t" // if (threadIdx.x == 0) there is a branch
		"@p   mov.u32 r, 0x00000080;			 \n\t" // then r = 0x00 00 00 80 (for ssv only:0x80=10000000)
		"or.b32      %0, t, r;                   \n\t" // merge t = 0xXXXXXX00 with r = 0x00000080 for next row
		"}"
		: "=r"(r) : "r"(a));      /* r is %0, a is %1 */
	return r;
}

#endif /* SIMD_FUNS__ */