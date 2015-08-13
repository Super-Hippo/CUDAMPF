#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "header_def.h"


/**/
/* Transfer seq into int value for optimization in Kernel
1. padding to 128x
2. unsigned char for each residue
3. union {unsigned char; unsigned char [4];}
4. store as 2-D array for padding seq
5. store length of padding seq as 1-D array
*/
unsigned int seq_Padding(unsigned int** iSeq, unsigned int* iLen, char** seq, unsigned int* length, int num)
{
	unsigned int i, j;
	unsigned int a, b, c;
	unsigned int sum = 0;

	union { unsigned int vec; unsigned char res[4]; } tmp = { 0 };

	for (i = 0; i < num; i++) {

		a = length[i] / 4;
		b = length[i] % 4;
		c = (b == 0) ? a : (a + 1);		/* 4 residue (CHARS) in 1 (INT) */

		a = c / 32;						/* keep reading the address of 128B */
		b = c % 32;						/* 32 in 1 for global coalesce */
		c = (b == 0) ? a : (a + 1);		/* times of 32 ints */

		iLen[i] = c * 32;				/* length of padding seq */

		iSeq[i] = (unsigned int*)malloc(c * 32 * sizeof(unsigned int));
		memset(iSeq[i], 0, c * 32 * sizeof(unsigned int));          	  	 /* initial to 0 */

		/* begin fill */
		for (j = 0; j < c * 32; j++) {
			tmp.res[0] = (j * 4 + 0 < length[i]) ? Amino_Offset(seq[i][j * 4 + 0]) : 31;
			tmp.res[1] = (j * 4 + 1 < length[i]) ? Amino_Offset(seq[i][j * 4 + 1]) : 31;
			tmp.res[2] = (j * 4 + 2 < length[i]) ? Amino_Offset(seq[i][j * 4 + 2]) : 31;
			tmp.res[3] = (j * 4 + 3 < length[i]) ? Amino_Offset(seq[i][j * 4 + 3]) : 31;
			iSeq[i][j] = tmp.vec;
			//iSeq[sum + j] = tmp.vec;	/* assign one INT */
		}

		sum += c * 32;					/* for allocating 1D iSeq */
	}

	return sum;		// for 1D array alloction later
}

/*****************************************************************
* 2. Standard iid null model ("null1")
*****************************************************************/
int NullOne(unsigned int L, float *ret_sc)
{
	float p1 = 0.0f;

	p1 = (float)L / (float)(L + 1);

	*ret_sc = (float)L * logf(p1) + logf(1. - p1);

	return fileOK;
}

/* Function:  esl_gumbel_surv() */
double esl_gumbel_surv(double x, double mu, double lambda)
{
	double y = lambda*(x - mu);
	double ey = -exp(-y);

	/* Use 1-e^x ~ -x approximation here when e^-y is small. */
	if (fabs(ey) < eslSMALLX1) return -ey;
	else                       return 1 - exp(ey);
}


/* For protein models, default iid background frequencies */
/* Function:  p7_AminoFrequencies() */
int p7_AminoFrequencies(float *f)
{
	f[0] = 0.0787945;   /* A */
	f[1] = 0.0151600;   /* C */
	f[2] = 0.0535222;   /* D */
	f[3] = 0.0668298;   /* E */
	f[4] = 0.0397062;   /* F */
	f[5] = 0.0695071;   /* G */
	f[6] = 0.0229198;   /* H */
	f[7] = 0.0590092;   /* I */
	f[8] = 0.0594422;   /* K */
	f[9] = 0.0963728;   /* L */
	f[10] = 0.0237718;   /* M */
	f[11] = 0.0414386;   /* N */
	f[12] = 0.0482904;   /* P */
	f[13] = 0.0395639;   /* Q */
	f[14] = 0.0540978;   /* R */
	f[15] = 0.0683364;   /* S */
	f[16] = 0.0540687;   /* T */
	f[17] = 0.0673417;   /* V */
	f[18] = 0.0114135;   /* W */
	f[19] = 0.0304133;   /* Y */

	return fileOK;
}


/* Implement roundf() function
* Example: roundf(2.1) = 2; roundf(2.6) = 3
*      roundf(-2.1) = -2; roundf(-2.6) = -3
* Using existing functions: floor(); ceil();
* Example: floor(2.6) = 2; ceil(2.1) = 3
*/
int round_DIY(float input)
{
	int r;
	r = (input - floor(input) < 0.5) ? floor(input) : ceil(input);
	return r;
}

/* FREE mem of HMM model */
void freeHMM(HMMER_PROFILE* hmm)
{
	free(hmm->f);         /* 20 standard residues */
	free(hmm->MU);        /* 3 floats */
	free(hmm->LAMBDA);    /* 3 floats */

	free(hmm);            /* free 'om', not include above extra space */
	/* only include their pointer address: 3 float pointes: MU,LAMBDA,f */
}


/* get optimal grid for MSV with SMEM version */
/* "Q" is the times of inner loop
 * "smem" is the maximum shared memory per SMX
 * "ws" is warp size
 */ 
RIB_BLOCK get_opt_MSV(int Q, int smem, int ws)
{
	RIB_BLOCK rb;
	rb.warps = 0;
	rb.res_blocks = 0;

	int i, j;
	int old_value = 0;
	int new_value = 0;

	for(i = 1; i <= 32; i++)	/* "i" is the row in blocks (warp in block) */
	{
		for(j = 1; j <= 16; j++)	/* "j" is the resident blocks per SMX */
		{
			new_value = (Q + 1) * ws * 4 * i * j;

			if( (new_value <= smem) && (32 * i * j <= 2048) ) 	/* hard-coded 32,2048 */
			{
				if(new_value > old_value)
				{
					rb.warps = i;
					rb.res_blocks = j;
					old_value = new_value;
				} 
				else if ( (old_value == new_value) && (j < rb.res_blocks) )
				{
					rb.warps = i;
					rb.res_blocks = j;
					old_value = new_value;
				}
			}
		}
	}

	return rb;
}

/* get optimal grid for VIT with SMEM version */
/* "Q" is the times of inner loop
 * "smem" is the maximum shared memory per SMX
 * "ws" is warp size
 */ 
RIB_BLOCK get_opt_VIT(int Q, int smem, int ws)
{
	RIB_BLOCK rb;
	rb.warps = 0;
	rb.res_blocks = 0;

	int i, j;
	int old_value = 0;
	int new_value = 0;

	for(i = 1; i <= 32; i++)	/* "i" is the row in blocks (warp in block) */
	{
		for(j = 1; j <= 16; j++)	/* "j" is the resident blocks per SMX */
		{
			new_value = (3 * Q + 1) * ws * 4 * i * j;

			if( (new_value <= smem) && (32 * i * j <= 2048) ) 	/* hard-coded 32,2048 */
			{
				if(new_value > old_value)
				{
					rb.warps = i;
					rb.res_blocks = j;
					old_value = new_value;
				} 
				else if ( (old_value == new_value) && (j < rb.res_blocks) )
				{
					rb.warps = i;
					rb.res_blocks = j;
					old_value = new_value;
				}
			}
		}
	}

	return rb;
}