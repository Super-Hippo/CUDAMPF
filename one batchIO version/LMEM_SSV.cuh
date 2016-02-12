/*
 * local memory version
 * MMX is cached in local memory
 */

 
#ifndef _LOCALMEMORY_KERNEL_SSV_
#define _LOCALMEMORY_KERNEL_SSV_

																					
extern "C" __global__ void KERNEL(unsigned int* dseq, unsigned int total, unsigned int* doffset,													
double* dsc, unsigned int* dL, unsigned int* dL_6r, unsigned int* dmat,										
unsigned int dbase, unsigned int dbias, unsigned int dtbm, unsigned int dtec, 								
float dscale)																		
{																											
	volatile __shared__ unsigned int cache[RIB][32];														
	unsigned int MMX[SIZE];																					
	const int row = blockIdx.y * blockDim.y + threadIdx.y;													
	unsigned int xE, xJ, sv, mmx, i, j, q, res, tjb, off, Len;											
	unsigned int count = 0;																					
	unsigned int index = RIB * gridDim.y * count;																																					
	while (row + index < total)																				
	{																										
		Len = dL_6r[row + index];																			
		off = doffset[row + index];																			
		i = 0;																								
		mmx = 0;																					
		xJ = 0;																								
		xE = 0x80808080;																					
		sv = 0x80808080;																								
		tjb = -1.0f * roundf(dscale * logf(3.0f / (float)(dL[row + index] + 3)));							
		tjb = (tjb > 255.) ? 255 : (unsigned int)tjb;
		if (tjb + dtbm + dtec + dbias >= 127) 
		{
		 	dsc[row + index] = (double)0;
		 	goto Endseq;
		}																																
		for (q = 0; q < Q; q++)																				
			MMX[q] = 0x80808080;																						
		while (i < Len)																						
		{																									
 			cache[threadIdx.y][threadIdx.x] = dseq[off + i + threadIdx.x];									
 			for (j = 0; j < 32; j++)																		
 			{																								
 				res = (cache[threadIdx.y][j] & 0x000000ff);													
 				if (res == 31) break;	
 				res = res * Q * 32 + threadIdx.x;																	
 				mmx = reorder_uint8_ssv(sv);
#pragma unroll Q																	
 				for (q = 0; q < Q; q++)																		
 				{																																						
 					sv = vsubss4(mmx, __ldg(&dmat[res + q * 32]));					
 					xE = vmaxu4(xE, sv);																	
 					mmx = MMX[q];																			
 					MMX[q] = sv;																			
 				}																							
 				xE = vmaxu4(xE, __shfl_xor(xE, 16));														
 				xE = vmaxu4(xE, __shfl_xor(xE, 8));															
 				xE = vmaxu4(xE, __shfl_xor(xE, 4));															
 				xE = vmaxu4(xE, __shfl_xor(xE, 2));															
 				xE = vmaxu4(xE, __shfl_xor(xE, 1));															
 				xE = max_uint8(xE);																																				
 				res = ((cache[threadIdx.y][j] >> 8) & 0x000000ff);											
 				if (res == 31) break;																																						
 				res = res * Q * 32 + threadIdx.x;																	
 				mmx = reorder_uint8_ssv(sv);
#pragma unroll Q																	
 				for (q = 0; q < Q; q++)																		
 				{																																							
 					sv = vsubss4(mmx, __ldg(&dmat[res + q * 32]));					
 					xE = vmaxu4(xE, sv);																	
 					mmx = MMX[q];																			
 					MMX[q] = sv;																			
 				}																							
 				xE = vmaxu4(xE, __shfl_xor(xE, 16));														
 				xE = vmaxu4(xE, __shfl_xor(xE, 8));															
 				xE = vmaxu4(xE, __shfl_xor(xE, 4));															
 				xE = vmaxu4(xE, __shfl_xor(xE, 2));															
 				xE = vmaxu4(xE, __shfl_xor(xE, 1));															
 				xE = max_uint8(xE);																																				
 				res = ((cache[threadIdx.y][j] >> 16) & 0x000000ff);											
 				if (res == 31) break;																		
 				res = res * Q * 32 + threadIdx.x;																	
 				mmx = reorder_uint8_ssv(sv);
#pragma unroll Q																	
 				for (q = 0; q < Q; q++)																		
 				{																																							
 					sv = vsubss4(mmx, __ldg(&dmat[res + q * 32]));					
 					xE = vmaxu4(xE, sv);																	
 					mmx = MMX[q];																			
 					MMX[q] = sv;																			
 				}																							
 				xE = vmaxu4(xE, __shfl_xor(xE, 16));														
 				xE = vmaxu4(xE, __shfl_xor(xE, 8));															
 				xE = vmaxu4(xE, __shfl_xor(xE, 4));															
 				xE = vmaxu4(xE, __shfl_xor(xE, 2));															
 				xE = vmaxu4(xE, __shfl_xor(xE, 1));															
 				xE = max_uint8(xE);																																				
 				res = ((cache[threadIdx.y][j] >> 24) & 0x000000ff);											
 				if (res == 31) break;																		
 				res = res * Q * 32 + threadIdx.x;																	
 				mmx = reorder_uint8_ssv(sv);
#pragma unroll Q																	
 				for (q = 0; q < Q; q++)																		
 				{																																						
 					sv = vsubss4(mmx, __ldg(&dmat[res + q * 32]));					
 					xE = vmaxu4(xE, sv);																	
 					mmx = MMX[q];																			
 					MMX[q] = sv;																			
 				}																							
 				xE = vmaxu4(xE, __shfl_xor(xE, 16));														
 				xE = vmaxu4(xE, __shfl_xor(xE, 8));															
 				xE = vmaxu4(xE, __shfl_xor(xE, 4));															
 				xE = vmaxu4(xE, __shfl_xor(xE, 2));															
 				xE = vmaxu4(xE, __shfl_xor(xE, 1));															
 				xE = max_uint8(xE);																																			
 			}																								
 			i += 32;																						
		}
	 	xJ = extr_uint8(xE);
	 	if (xJ >= 255 - dbias) 
	 	{
	 		dsc[row + index] = (double)999999;
	 		if (dbase - tjb - dtbm < 128)
	 		{
	 			dsc[row + index] = (double)0;
	 			goto Endseq;
	 		}
	 		goto Endseq;
	 	}
	 	xJ += dbase - tjb - dtbm;
	 	xJ -= 128;
	 	if (xJ >= 255 - dbias)
	 	{
	 		dsc[row + index] = (double)999999;
	 		goto Endseq;
	 	}
	 	xJ = xJ - dtec;
	 	if (xJ > dbase)
	 	{
	 		dsc[row + index] = (double)0;
	 		goto Endseq;
	 	}
	 	dsc[row + index] = (double)(((float)(xJ - tjb) - (float)dbase) / dscale - (float)3.0);
	Endseq:
		count++;
		index = RIB * gridDim.y * count;
	}																										
}																											


#endif /* _LOCALMEMORY_KERNEL_SSV_ */
