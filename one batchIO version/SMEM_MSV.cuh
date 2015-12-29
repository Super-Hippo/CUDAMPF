/*
 * shared memory version
 * MMX is cached in shared memory
 */

 
#ifndef _SHAREDMEMORY_KERNEL_MSV_
#define _SHAREDMEMORY_KERNEL_MSV_
																				
extern "C" __global__ 																					
void KERNEL(unsigned int* dseq, unsigned int total, unsigned int* doffset,													
double* dsc, unsigned int* dL, unsigned int* dL_6r, unsigned int* dmat,										
unsigned int dbase, unsigned int dbias, unsigned int dtbm, unsigned int dtec, 								
float dscale, double mu, double lambda)																		
{																											
	volatile __shared__ unsigned int cache[RIB][32];														
	volatile __shared__ unsigned int MMX[RIB][Q * 32];														
	const int row = blockIdx.y * blockDim.y + threadIdx.y;													
	unsigned int xE, xJ, xB, sv, mmx, i, j, q, res, tjb, off, Len;											
	unsigned int count = 0;																					
	unsigned int index = RIB * gridDim.y * count;															
	float Score, Nullsc;																					
	double y, ey;																							
	while (row + index < total)																				
	{																										
		Len = dL_6r[row + index];																			
		off = doffset[row + index];																			
		i = 0;																								
		mmx = 0;																							
		xJ = 0;																								
		xE = 0;																								
		j = 0;																								
		res = 0;																							
		y = 0;																								
		ey = 0;																								
		Score = 0;																							
		sv = 0;																								
		tjb = -1.0f * roundf(dscale * logf(3.0f / (float)(dL[row + index] + 3)));							
		tjb = (tjb > 255.) ? 255 : (unsigned int)tjb;														
		tjb = dup_uint8(tjb);																				
		xB = vsubus4(dbase, tjb);																			
		Nullsc = (float)dL[row + index] * logf((float)dL[row + index] / ((float)dL[row + index] + (float)1)) + logf((float)1 - (float)dL[row + index] / ((float)dL[row + index] + (float)1));	
		for (q = 0; q < Q; q++)																				
			MMX[threadIdx.y][q * 32 + threadIdx.x] = 0;														
		while (i < Len)																						
		{																									
			cache[threadIdx.y][threadIdx.x] = dseq[off + i + threadIdx.x];									
			for (j = 0; j < 32; j++)																		
			{																								
				res = (cache[threadIdx.y][j] & 0x000000ff);													
				if (res == 31) break;																		
				xE = 0;																						
				xB = vsubus4(xB, dtbm);																		
				mmx = reorder_uint8(sv);																	
				for (q = 0; q < Q; q++)																		
				{																							
					sv = vmaxu4(mmx, xB);																	
					sv = vaddus4(sv, dbias);																
					sv = vsubus4(sv, __ldg(&dmat[res * Q * 32 + q * 32 + threadIdx.x]));					
					xE = vmaxu4(xE, sv);																	
					mmx = MMX[threadIdx.y][q * 32 + threadIdx.x];											
					MMX[threadIdx.y][q * 32 + threadIdx.x] = sv;											
				}																							
				xE = vmaxu4(xE, __shfl_xor(xE, 16));														
				xE = vmaxu4(xE, __shfl_xor(xE, 8));															
				xE = vmaxu4(xE, __shfl_xor(xE, 4));															
				xE = vmaxu4(xE, __shfl_xor(xE, 2));															
				xE = vmaxu4(xE, __shfl_xor(xE, 1));															
				xE = max_uint8(xE);																			
				if (vaddus4(xE, dbias) == 0xffffffff)														
				{																							
					dsc[row + index] = (double)999999;														
					break;																					
				}																							
				xE = vsubus4(xE, dtec);																		
				xJ = vmaxu4(xJ, xE);																		
				xB = vmaxu4(dbase, xJ);																		
				xB = vsubus4(xB, tjb);																		
				res = ((cache[threadIdx.y][j] >> 8) & 0x000000ff);											
				if (res == 31) break;																		
				xE = 0;																						
				xB = vsubus4(xB, dtbm);																		
				mmx = reorder_uint8(sv);																	
				for (q = 0; q < Q; q++)																		
				{																							
					sv = vmaxu4(mmx, xB);																	
					sv = vaddus4(sv, dbias);																
					sv = vsubus4(sv, __ldg(&dmat[res * Q * 32 + q * 32 + threadIdx.x]));					
					xE = vmaxu4(xE, sv);																	
					mmx = MMX[threadIdx.y][q * 32 + threadIdx.x];											
					MMX[threadIdx.y][q * 32 + threadIdx.x] = sv;											
				}																							
				xE = vmaxu4(xE, __shfl_xor(xE, 16));														
				xE = vmaxu4(xE, __shfl_xor(xE, 8));															
				xE = vmaxu4(xE, __shfl_xor(xE, 4));															
				xE = vmaxu4(xE, __shfl_xor(xE, 2));															
				xE = vmaxu4(xE, __shfl_xor(xE, 1));															
				xE = max_uint8(xE);																			
				if (vaddus4(xE, dbias) == 0xffffffff)														
				{																							
					dsc[row + index] = (double)999999;														
					break;																					
				}																							
				xE = vsubus4(xE, dtec);																		
				xJ = vmaxu4(xJ, xE);																		
				xB = vmaxu4(dbase, xJ);																		
				xB = vsubus4(xB, tjb);																		
				res = ((cache[threadIdx.y][j] >> 16) & 0x000000ff);											
				if (res == 31) break;																		
				xE = 0;																						
				xB = vsubus4(xB, dtbm);																		
				mmx = reorder_uint8(sv);																	
				for (q = 0; q < Q; q++)																		
				{																							
					sv = vmaxu4(mmx, xB);																	
					sv = vaddus4(sv, dbias);																
					sv = vsubus4(sv, __ldg(&dmat[res * Q * 32 + q * 32 + threadIdx.x]));					
					xE = vmaxu4(xE, sv);																	
					mmx = MMX[threadIdx.y][q * 32 + threadIdx.x];											
					MMX[threadIdx.y][q * 32 + threadIdx.x] = sv;											
				}																							
				xE = vmaxu4(xE, __shfl_xor(xE, 16));														
				xE = vmaxu4(xE, __shfl_xor(xE, 8));															
				xE = vmaxu4(xE, __shfl_xor(xE, 4));															
				xE = vmaxu4(xE, __shfl_xor(xE, 2));															
				xE = vmaxu4(xE, __shfl_xor(xE, 1));															
				xE = max_uint8(xE);																			
				if (vaddus4(xE, dbias) == 0xffffffff)														
				{																							
					dsc[row + index] = (double)999999;														
					break;																					
				}																							
				xE = vsubus4(xE, dtec);																		
				xJ = vmaxu4(xJ, xE);																		
				xB = vmaxu4(dbase, xJ);																		
				xB = vsubus4(xB, tjb);																		
				res = ((cache[threadIdx.y][j] >> 24) & 0x000000ff);											
				if (res == 31) break;																		
				xE = 0;																						
				xB = vsubus4(xB, dtbm);																		
				mmx = reorder_uint8(sv);																	
				for (q = 0; q < Q; q++)																		
				{																							
					sv = vmaxu4(mmx, xB);																	
					sv = vaddus4(sv, dbias);																
					sv = vsubus4(sv, __ldg(&dmat[res * Q * 32 + q * 32 + threadIdx.x]));					
					xE = vmaxu4(xE, sv);																	
					mmx = MMX[threadIdx.y][q * 32 + threadIdx.x];											
					MMX[threadIdx.y][q * 32 + threadIdx.x] = sv;											
				}																							
				xE = vmaxu4(xE, __shfl_xor(xE, 16));														
				xE = vmaxu4(xE, __shfl_xor(xE, 8));															
				xE = vmaxu4(xE, __shfl_xor(xE, 4));															
				xE = vmaxu4(xE, __shfl_xor(xE, 2));															
				xE = vmaxu4(xE, __shfl_xor(xE, 1));															
				xE = max_uint8(xE);																			
				if (vaddus4(xE, dbias) == 0xffffffff)														
				{																							
					dsc[row + index] = (double)999999;														
					break;																					
				}																							
				xE = vsubus4(xE, dtec);																		
				xJ = vmaxu4(xJ, xE);																		
				xB = vmaxu4(dbase, xJ);																		
				xB = vsubus4(xB, tjb);																		
			}																								
			i += 32;																						
		}																									
		if (fabs(dsc[row + index] - (double)999999.0f) < 1e-6)												
		{																									
			y = lambda * (dsc[row + index] - mu);			
			ey = -exp(-y);									
			if (fabs(ey) < (double)5e-9)					
			{												
				dsc[row + index] = -ey;						
			}												
			else {											
				dsc[row + index] = (double)1 - exp(ey);		
			}												
		}													
		else {												
			Score = ((float)(extr_uint8(xJ) - extr_uint8(tjb)) - (float)extr_uint8(dbase)) / dscale - (float)3.0;	
			Score = (Score - Nullsc) / (float)0.69314718055994529;													
			y = lambda * ((double)Score - mu);																		
			ey = -exp(-y);																							
			if (fabs(ey) < (double)5e-9)																			
			{																										
				dsc[row + index] = -ey;																				
			}																										
			else {																									
				dsc[row + index] = (double)1 - exp(ey);																
			}																										
		}																									
		count++;																							
		index = RIB * gridDim.y * count;																	
	}																										
}																											


#endif /* _SHAREDMEMORY_KERNEL_MSV_ */
