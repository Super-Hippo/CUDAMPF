/*
 * Local memory version
 * MMX/IMX/DMX are stored in local memory
*/


#ifndef _LOCALMEMORY_KERNEL_VIT_
#define _LOCALMEMORY_KERNEL_VIT_

extern "C" __global__
void KERNEL(unsigned int* seq, unsigned int total, unsigned int* offset,
double* sc, int* L, unsigned int* L_6r, int* mat, int* tran,
int base, int e_lm, int ddbound, float scale, int QV,
double mu, double lambda)									
{
	volatile __shared__ unsigned int cache[RIB][32];
	int MMX[SIZE];
	int IMX[SIZE];
	int DMX[SIZE];
	const int row = blockIdx.y * blockDim.y + threadIdx.y;
	int xE, xJ, xB, xN, xC;
	int mmx, imx, dmx;
	int NCJ_MOVE;
	int sv, dcv, Dmax;
	unsigned int i, j, q, z;
	unsigned int LEN, OFF, res, res_s;
	unsigned int count = 0;
	unsigned int index = 0;
	xE = 0x80008000;
	xJ = 0x80008000;
	xB = 0x80008000;
	xN = 0x80008000;
	xC = 0x80008000;
	mmx = 0x80008000;
	imx = 0x80008000;
	dmx = 0x80008000;
	NCJ_MOVE = 0x80008000;
	sv = 0x80008000;
	dcv = 0x80008000;
	Dmax = 0x80008000;
	i = 0;
	j = 0;
	q = 0;
	z = 0;
	LEN = 0;
	OFF = 0;
	res = 0;
	res_s = 0;
	float SCORE;
	float Nullsc;
	double y, ey;
	while (row + index < total)
	{
		LEN = L_6r[row + index];
		OFF = offset[row + index];
		NCJ_MOVE = rintf(scale * logf(3.0f / (float)(L[row + index] + 3.0f)));
		if (NCJ_MOVE >= 32767)
		{
			NCJ_MOVE = 0xffffffff;
		}
		else if (NCJ_MOVE <= -32768)
		{
			NCJ_MOVE = 0x80008000;
		}
		else 
		{
			NCJ_MOVE = dup_int16(NCJ_MOVE);
		}
		i = 0;
		sv = 0x80008000;
		xJ = 0x80008000;
		xC = 0x80008000;
		xN = base;
		xB = vaddss2(xN, NCJ_MOVE);
		Nullsc = (float)L[row + index] * logf((float)L[row + index] / ((float)L[row + index] + (float)1)) + logf((float)1 - (float)L[row + index] / ((float)L[row + index] + (float)1));
		for (q = 0; q < QV; q++)
		{
			MMX[q] = 0x80008000;
			IMX[q] = 0x80008000;
			DMX[q] = 0x80008000;
		}
		for (i = 0; i < LEN; i += 32)
		{
			cache[threadIdx.y][threadIdx.x] = seq[OFF + i + threadIdx.x];
			for (j = 0; j < 32; j++)
			{
				res = cache[threadIdx.y][j];	
				if ((res & 0x000000ff) == 31) goto Endseq;

				for (z = 0; z < 4; z++)
				{
					res_s = ((res >> (8 * z)) & 0x000000ff);
					if (res_s == 31) goto Endseq;
					res_s *= QV * 32;
					xE = 0x80008000;
					Dmax = 0x80008000;
					dcv = 0x80008000;
					mmx = MMX[QV - 1];
					mmx = reorder_int16(mmx);
					imx = IMX[QV - 1];
					imx = reorder_int16(imx);
					dmx = DMX[QV - 1];
					dmx = reorder_int16(dmx);
					for (q = 0; q < QV; q++)
					{
						sv = vaddss2(xB, __ldg(&tran[q * 224 + 0 * 32 + threadIdx.x]));	/* 0 * 32 = B_M * 32 */
						sv = vmaxs2(sv, vaddss2(mmx, __ldg(&tran[q * 224 + 1 * 32 + threadIdx.x])));	/* 1 * 32 = M_M * 32 */
						sv = vmaxs2(sv, vaddss2(imx, __ldg(&tran[q * 224 + 2 * 32 + threadIdx.x])));	/* 2 * 32 = I_M * 32 */
						sv = vmaxs2(sv, vaddss2(dmx, __ldg(&tran[q * 224 + 3 * 32 + threadIdx.x])));	/* 3 * 32 = D_M * 32 */
						sv = vaddss2(sv, __ldg(&mat[res_s + q * 32 + threadIdx.x]));
						xE = vmaxs2(sv, xE);
						mmx = MMX[q];
						imx = IMX[q];
						dmx = DMX[q];
						MMX[q] = sv;
						DMX[q] = dcv;
						dcv = vaddss2(sv, __ldg(&tran[q * 224 + 4 * 32 + threadIdx.x]));	/* 4 * 32 = M_D * 32 */
						Dmax = vmaxs2(dcv, Dmax);
						sv = vaddss2(mmx, __ldg(&tran[q * 224 + 5 * 32 + threadIdx.x]));	/* 5 * 32 = M_I * 32 */
						sv = vmaxs2(sv, vaddss2(imx, __ldg(&tran[q * 224 + 6 * 32 + threadIdx.x])));	/* 6 * 32 I_I * 32 */
						IMX[q] = sv;
					}
					xE = vmaxs2(xE, __shfl_xor(xE, 16));
					xE = vmaxs2(xE, __shfl_xor(xE, 8));
					xE = vmaxs2(xE, __shfl_xor(xE, 4));
					xE = vmaxs2(xE, __shfl_xor(xE, 2));
					xE = vmaxs2(xE, __shfl_xor(xE, 1));
					xE = max_int16(xE);
					if (extr_int16(xE) >= 32767)
					{
						sc[row + index] = (double)999999;
						goto Endseq;
					}
					xC = vmaxs2(xC, vaddss2(xE, e_lm));
					xJ = vmaxs2(xJ, vaddss2(xE, e_lm));
					xB = vmaxs2(vaddss2(xJ, NCJ_MOVE), vaddss2(xN, NCJ_MOVE));
					Dmax = vmaxs2(Dmax, __shfl_xor(Dmax, 16));
					Dmax = vmaxs2(Dmax, __shfl_xor(Dmax, 8));
					Dmax = vmaxs2(Dmax, __shfl_xor(Dmax, 4));
					Dmax = vmaxs2(Dmax, __shfl_xor(Dmax, 2));
					Dmax = vmaxs2(Dmax, __shfl_xor(Dmax, 1));
					Dmax = max_int16(Dmax);
					if ((extr_int16(Dmax) + extr_int16(ddbound)) > extr_int16(xB))
					{
						dcv = reorder_int16(dcv);
						for (q = 0; q < QV; q++)
						{
							DMX[q] = vmaxs2(dcv, DMX[q]);
							dcv = vaddss2(DMX[q], __ldg(&tran[QV * 224 + q * 32 + threadIdx.x]));
						}
						do {
							dcv = reorder_int16(dcv);
							for (q = 0; q < QV; q++)
							{
								if (__any(vsetgts2(dcv, DMX[q])) == 0)
									break;
								DMX[q] = vmaxs2(dcv, DMX[q]);
								dcv = vaddss2(DMX[q], __ldg(&tran[QV * 224 + q * 32 + threadIdx.x]));
							}
						} while (q == QV);
					}
					else {
						DMX[0] = reorder_int16(dcv);
					}
				}
			}
		}
		Endseq:;
		if (fabs(sc[row + index] - (double)999999) < 1e-6)
		{
			y = lambda * (sc[row + index] - mu);
			ey = -exp(-y);
			if (fabs(ey) < (double)5e-9)
			{
				sc[row + index] = -ey;
			}
			else {
				sc[row + index] = (double)1 - exp(ey);
			}
		}
		else if (extr_int16(xC) > -32768)
		{
			SCORE = ((float)extr_int16(xC) + (float)extr_int16(NCJ_MOVE) - (float)extr_int16(base)) / scale - (float)3.0;
			SCORE = (SCORE - Nullsc) / (float)0.69314718055994529;
			y = lambda * ((double)SCORE - mu);
			ey = -exp(-y);
			if (fabs(ey) < (double)5e-9)
			{
				sc[row + index] = -ey;
			}
			else {
				sc[row + index] = (double)1 - exp(ey);
			}
		}
		else
		{
			sc[row + index] = (double)-999999;
		}
		count++;
		index = RIB * gridDim.y * count;
	}
}


#endif /* _LOCALMEMORY_KERNEL_VIT_ */
