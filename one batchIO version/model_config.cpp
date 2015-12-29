/* Configure hmm model
* includes special nodes transition
* expon to prob from score
* entry score
* log-odd to score from prob
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>

#include "header_def.h"


/*****************************************************************
* 1. The HMMER_PROFILE structure: a score profile.
*****************************************************************/
HMMER_PROFILE* hmmer_profile_Create(int seq_size, int hmm_size)
{
	HMMER_PROFILE *om = NULL;                            /* claim here for return at last */
	om = (HMMER_PROFILE*)malloc(sizeof(HMMER_PROFILE));  /* this is an important step: after this,
														  * we can do any assignment works for all
														  * elements of om. Otherwise, we cannot do
														  * none of below assignment works
														  */

	/* ================================================== */
	/* Initialization for raw data and necessary variance */
	/* ================================================== */

	om->M = hmm_size;                                      	/* number of nodes in the model */
	om->MU = (double*)malloc(3 * sizeof(double));     		/* respectively, for MSV, VIT and FWD */
	om->LAMBDA = (double*)malloc(3 * sizeof(double));			/* respectively, for MSV, VIT and FWD */
	om->f = (float*)malloc(20 * sizeof(float));     		/* above 3 pointer get their own space which is not included in om's space */
	if (p7_AminoFrequencies(om->f) != 1) printf("error\n");	/* For protein models, default iid background frequencies */

	unsigned int i, j;

	/* Dynamic allocation 2D array (temperary) */
	om->mat_32bits = (float**)malloc((om->M + 1) * sizeof(float*));
	for (i = 0; i < (om->M + 1); i++)
	{
		om->mat_32bits[i] = (float*)malloc(PROTEIN_TYPE * sizeof(float));
		for (j = 0; j < PROTEIN_TYPE; j++)
			om->mat_32bits[i][j] = 0.0f;
	}

	om->ins_32bits = (float**)malloc((om->M + 1) * sizeof(float*));
	for (i = 0; i < (om->M + 1); i++)
	{
		om->ins_32bits[i] = (float*)malloc(PROTEIN_TYPE * sizeof(float));
		for (j = 0; j < PROTEIN_TYPE; j++)
			om->ins_32bits[i][j] = 0.0f;
	}

	om->tran_32bits = (float**)malloc((om->M + 1) * sizeof(float*));
	for (i = 0; i < (om->M + 1); i++)
	{
		om->tran_32bits[i] = (float*)malloc(7 * sizeof(float));
		for (j = 0; j < 7; j++)
			om->tran_32bits[i][j] = 0.0f;
	}

	om->log_tran_32bits = (float**)malloc((om->M) * sizeof(float*));
	for (i = 0; i < (om->M); i++)
	{
		om->log_tran_32bits[i] = (float*)malloc(TRANS_TYPE * sizeof(float));
		for (j = 0; j < TRANS_TYPE; j++)
			om->log_tran_32bits[i][j] = 0.0f;
	}

	/* ========================================================================= */
	/* MSVFilter.SSVFilter uses scaled, biased uchars: 16x unsigned byte vectors */
	/* ========================================================================= */

	om->msvQ = NQB(hmm_size);				printf("msv:%d\n", om->msvQ);										/* 128x parallel */
	om->msv_vec = (__32uint__*)malloc(om->msvQ * PROTEIN_TYPE * sizeof(__32uint__));
	om->ssv_vec = (__32uint__*)malloc(om->msvQ * PROTEIN_TYPE * sizeof(__32uint__));		/*ssv*/

	for (int j = 0; j < om->msvQ * PROTEIN_TYPE; j++) {
		memset(om->msv_vec[j], 0, sizeof(om->msv_vec[j]));  	/* initialization: each "msv_vec[j]" is an array consist of 32 int values */
		memset(om->ssv_vec[j], 0, sizeof(om->ssv_vec[j]));		/*ssv*/
	}

	om->tbm_b = 0;
	om->tec_b = 0;
	om->tjb_b = 0;
	om->base_b = 0;
	om->bias_b = 0;
	om->scale_b = 0.0f;
	om->nj = 1.0f;          /* H3 only uses multihit local alignment ( --fs mode) */

	om->tbm_s = 0;
	om->tec_s = 0;
	om->tjb_s = 0;
	om->base_s = 0;
	om->bias_s = 0;

	/* ================================================================== */
	/* ViterbiFilter uses scaled swords: 8x signed 16-bit integer vectors */
	/* ================================================================== */

	om->vitQ = NQW(hmm_size);					printf("vit %d\n",om->vitQ);									/* 128x parallel */
	om->vit_vec = (__32int__*)malloc(om->vitQ * PROTEIN_TYPE * sizeof(__32int__));	/* allocate (vitQ * PROTEIN_TYPE) arrays */
	om->trans_vec = (__32int__*)malloc(om->vitQ * TRANS_TYPE * sizeof(__32int__));	/* allocate (vitQ * TRANS_TYPE) arrays */

	om->ncj_move = NULL;	/* getting ncj_move out of VIT kernel (Temporary: we may move it into kernel again) */

	for (int j = 0; j < om->vitQ * PROTEIN_TYPE; j++)
		memset(om->vit_vec[j], 0, sizeof(om->vit_vec[j]));		/* initialize those arrays */

	for (int j = 0; j < om->vitQ * TRANS_TYPE; j++)
		memset(om->trans_vec[j], 0, sizeof(om->trans_vec[j]));	/* initialize those arrays */

	om->scale_w = 0.0f;
	om->base_w = 0;
	om->ddbound_w = 0;
	om->ncj_roundoff = 0.0f;

	om->base_vs = 0;
	om->E_lm = 0;
	om->ddbound_vs = 0;

	return om;
}

/* Get the log-odd score */
int log_Odd_score(HMMER_PROFILE *hmm)
{
	int i, j;

	/* 1. Match Emission */
	for (j = 0; j < PROTEIN_TYPE; j++)
		hmm->mat_32bits[0][j] = minusInfinity;      /* For initialize Edge values. what ever the degen is, we set the Edge -infinity */

	for (i = 1; i <= hmm->M; i++) {
		for (j = 0; j < 20; j++)
			hmm->mat_32bits[i][j] = logf(hmm->mat_32bits[i][j] / hmm->f[j]);
	}

	return fileOK;
}

/* Function:  p7_hmm_CalculateOccupancy() */
int p7_hmm_CalculateOccupancy(HMMER_PROFILE *hmm, float *mocc)
{
	int k;

	mocc[0] = 0.;                                                       /* no M_0 state */
	mocc[1] = hmm->tran_32bits[0][MI] + hmm->tran_32bits[0][MM];        /* initialize w/ 1 - B->D_1 */
	for (k = 2; k <= hmm->M; k++) {
		mocc[k] = mocc[k - 1] * (hmm->tran_32bits[k - 1][MM] + hmm->tran_32bits[k - 1][MI]) + (1.0 - mocc[k - 1]) * hmm->tran_32bits[k - 1][DM];
	}

	return fileOK;
}

int get_entryScore(HMMER_PROFILE *hmm)
{
	int k;
	float Z = 0.;
	float* occ = NULL;

	occ = (float*)malloc((hmm->M + 1) * sizeof(float));
	if (p7_hmm_CalculateOccupancy(hmm, occ) != fileOK) printf("error in get Occupancy!\n");

	for (k = 1; k <= hmm->M; k++){              /* since occ[0] is 0, so we start from k = 1 */
		Z += occ[k] * (float)(hmm->M - k + 1);
	}

	for (k = 1; k <= hmm->M; k++)
	{
		hmm->log_tran_32bits[k - 1][0] = logf(occ[k] / Z);      /* note off-by-one: entry at Mk stored as [k-1][BM] */
	}                                                       	/* '0' represents the first column, which is the "BM" */

	free(occ);

	return fileOK;
}

/* Get the log of main transition */
int log_Trans(HMMER_PROFILE *hmm)
{
	int i;

	/* only node 1 to M-1 have meaningful values */
	for (i = 1; i < hmm->M; i++) {
		hmm->log_tran_32bits[i][M_M] = logf(hmm->tran_32bits[i][MM]); 
		hmm->log_tran_32bits[i][I_M] = logf(hmm->tran_32bits[i][IM]);
		hmm->log_tran_32bits[i][D_M] = logf(hmm->tran_32bits[i][DM]);
	}

	/* off-by-one */
	for (i = 0; i < (hmm->M - 1); i++) {
		hmm->log_tran_32bits[i][M_D] = logf(hmm->tran_32bits[i + 1][MD]);
		hmm->log_tran_32bits[i][M_I] = logf(hmm->tran_32bits[i + 1][MI]);
		hmm->log_tran_32bits[i][D_D] = logf(hmm->tran_32bits[i + 1][DD]);
		hmm->log_tran_32bits[i][I_I] = logf(hmm->tran_32bits[i + 1][II]);
	}

	/* SET "-infinity" to Node 0 for MM,IM,DM */
	hmm->log_tran_32bits[0][M_M] = minusInfinity;
	hmm->log_tran_32bits[0][I_M] = minusInfinity;
	hmm->log_tran_32bits[0][D_M] = minusInfinity;

	/* MD, MI, II, DD */
	hmm->log_tran_32bits[hmm->M - 1][M_D] = minusInfinity;
	hmm->log_tran_32bits[hmm->M - 1][M_I] = minusInfinity;
	hmm->log_tran_32bits[hmm->M - 1][I_I] = minusInfinity;
	hmm->log_tran_32bits[hmm->M - 1][D_D] = minusInfinity;

	return fileOK;
}

/* Get special node transition (1-D arrary with 8 values)*/
int xTrans(HMMER_PROFILE *hmm)
{
	hmm->Xtran_32bits[E * XTRANS_TYPE + LOOP] = -eslCONST_LOG2;      /* H3 is --fs mode, which is multihit local alignment */
	hmm->Xtran_32bits[E * XTRANS_TYPE + MOVE] = -eslCONST_LOG2;      /* too */

	hmm->Xtran_32bits[N * XTRANS_TYPE + LOOP] = 0;                   /* Others need reconfigure based on different 'L', here just simply initialize to 0 */
	hmm->Xtran_32bits[N * XTRANS_TYPE + MOVE] = 0;
	hmm->Xtran_32bits[J * XTRANS_TYPE + LOOP] = 0;
	hmm->Xtran_32bits[J * XTRANS_TYPE + MOVE] = 0;
	hmm->Xtran_32bits[C * XTRANS_TYPE + LOOP] = 0;
	hmm->Xtran_32bits[C * XTRANS_TYPE + MOVE] = 0;

	return fileOK;
}

/* ====================================================== Configuration for MSV ======================================================================= */

/* mf_conversion() */
int mf_conversion(HMMER_PROFILE *hmm)
{
	float   max_rs = 0.0;                           		/* maximum residue score: used for unsigned emission score bias */
	int i, j, q, z;                                     	/* index */
	union { __32uint__ v; unsigned char i[128]; } tmp;  	/* use union for following assignment to each unsigned char value */
	union { __32uint__ v; unsigned char i[128]; } tmp_ssv;	/* for SSV */

	/* Here, we get the 'max_rs' values */
	for (i = 0; i <= hmm->M; i++) 
	{
		for (j = 0; j < 20; j++)
		{
			max_rs = _MAX(max_rs, hmm->mat_32bits[i][j]);
		}
	}

	hmm->scale_b = 3.0 / eslCONST_LOG2;                               /* scores in units of third-bits */
	hmm->base_b = 190;
	hmm->bias_b = unbiased_byteify(hmm, -1.0 * max_rs);              /* here, we use 'max_rs' */

	/* only for SSV */
	char val_1 = (char)(hmm->bias_b + 127);
	char val_2 = 127;

	/* V3: Striped Match emission cost */
	for (i = 0; i < PROTEIN_TYPE; i++)
	{
		for (q = 0, j = 1; q < hmm->msvQ; q++, j++)
		{
			for (z = 0; z < 128; z++) {
				tmp.i[z] = ((j + z * hmm->msvQ <= hmm->M) ? biased_byteify(hmm, hmm->mat_32bits[j + z * hmm->msvQ][i]) : 255);
				tmp_ssv.i[z] = ((((unsigned char)val_1 - (unsigned char)tmp.i[z]) < 0) ? 0 : ((unsigned char)val_1 - (unsigned char)tmp.i[z])) ^ val_2;
			}

			for (int l = 0; l < 32; l++) {
				hmm->msv_vec[i * hmm->msvQ + q][l] = tmp.v[l];		/* array copy */
				hmm->ssv_vec[i * hmm->msvQ + q][l] = tmp_ssv.v[l];	/*ssv*/
				//printf("%d ", tmp_ssv.v[l]);
			}
		}
	}

	/* transition costs */
	hmm->tbm_b = unbiased_byteify(hmm, logf(2.0f / ((float)hmm->M * (float)(hmm->M + 1))));   /* constant B->Mk penalty        */
	hmm->tec_b = unbiased_byteify(hmm, logf(0.5f));                                           /* constant multihit E->C = E->J */

	/* SIMD-lize base, bias, tbm, tec */
	hmm->base_s = simdlize_MSV(hmm->base_b);
	hmm->bias_s = simdlize_MSV(hmm->bias_b);
	hmm->tec_s = simdlize_MSV(hmm->tec_b);
	hmm->tbm_s = simdlize_MSV(hmm->tbm_b);

	return fileOK;
}

/* SIMD-lize function
* only for MSV (8-bits unchar)
*/
unsigned int simdlize_MSV(unsigned int a)
{
	union{ unsigned int v; unsigned char i[4]; } tmp;
	int p;
	for (p = 0; p < 4; p++)
		tmp.i[p] = a & 0x000000ff;
	return tmp.v;
}

/* unbiased_byteify() */
unsigned char unbiased_byteify(HMMER_PROFILE *hmm, float sc)
{
	unsigned char b;

	sc = -1.0f * round_DIY(hmm->scale_b * sc);       /* ugh. sc is now an integer cost represented in a float...    */
	b = (sc > 255.) ? 255 : (unsigned char)sc;     	/* and now we cast and saturate it to an unsigned char cost... */
	return b;
}

/* biased_byteify() */
unsigned char biased_byteify(HMMER_PROFILE *hmm, float sc)
{
	unsigned char b;

	sc = -1.0f * round_DIY(hmm->scale_b * sc);                            	/* ugh. sc is now an integer cost represented in a float...           */
	b = (sc > 255 - hmm->bias_b) ? 255 : (unsigned char)sc + hmm->bias_b;   /* and now we cast, saturate, and bias it to an unsigned char cost... */
	return b;
}

/* ====================================================== Configuration for Viterbi ===================================================================== */

/* wordify() */
short wordify(HMMER_PROFILE* hmm, float sc)
{
	sc = round_DIY(hmm->scale_w * sc);
	if (sc >= 32767.0) return  32767;
	else if (sc <= -32768.0) return -32768;
	else return (short)sc;
}

/* SIMD-lize function
* only for vit (16-bits word)
* not sure whether it works for negative values...
*/
int simdlize_VIT(int a)
{
	union{ int v; short i[2]; } tmp;
	int p;
	for (p = 0; p < 2; p++)
		tmp.i[p] = (short)a;
	return tmp.v;
}

/* vf_conversion() */
int vf_conversion(HMMER_PROFILE* hmm)
{
	int     ddtmp;    				   		   /* used in finding worst DD transition bound */
	int i, j, q, z, l, k, t;           		   /* index */
	int maxval;						   		   /* constraint */
	short val;						   		   /* temp value */
	union { __32int__ v; short i[64]; } tmp;   /* use union for following assignment to each (short) value */

	hmm->scale_w = 500.0 / eslCONST_LOG2;
	hmm->base_w = 12000;

	/* V3: Striped Match emission cost */
	for (i = 0; i < PROTEIN_TYPE; i++)
	{
		for (q = 0, j = 1; q < hmm->vitQ; q++, j++)
		{
			for (z = 0; z < 64; z++)
				tmp.i[z] = ((j + z * hmm->vitQ <= hmm->M) ? wordify(hmm, hmm->mat_32bits[j + z * hmm->vitQ][i]) : 32768);

			for (l = 0; l < 32; l++)
				hmm->vit_vec[i * hmm->vitQ + q][l] = tmp.v[l];
		}
	}

	/* Transition costs, all but the DD's. */
	for (k = 0, q = 0; q < hmm->vitQ; q++, k++)
	{
		for (t = B_M; t <= I_I; t++) 
		{
			switch (t)
			{
			case B_M: maxval = 0; break; /* gm has tBMk stored off by one! start from k=0 not 1   */
			case M_M: maxval = 0; break; /* MM, DM, IM vectors are rotated by -1, start from k=0  */
			case I_M: maxval = 0; break;
			case D_M: maxval = 0; break;
			case M_D: maxval = 0; break; /* the remaining ones are straight up  */
			case M_I: maxval = 0; break;
			case I_I: maxval = -1; break;
			}

			for (z = 0; z < 64; z++)
			{
				val = ((k + z * hmm->vitQ < hmm->M) ? wordify(hmm, hmm->log_tran_32bits[k + z * hmm->vitQ][t]) : -32768);
				tmp.i[z] = (val <= maxval) ? val : maxval;   /* do not allow an II transition cost of 0, or hell may occur. */
			}

			for (l = 0; l < 32; l++)
				hmm->trans_vec[q * 7 + t][l] = tmp.v[l];	// 7 is hard-coded since we have BM,MM,IM,DM,MD,MI,II in this loop...
		}
	}

	/* Finally the DD's, which are at the end of the optimized tsc vector; (j is already sitting there) */
	for (k = 0, q = 0; q < hmm->vitQ; q++, k++)				// k = 0 is only for our case
	{
		for (z = 0; z < 64; z++)
			tmp.i[z] = ((k + z * hmm->vitQ < hmm->M) ? wordify(hmm, hmm->log_tran_32bits[k + z * hmm->vitQ][D_D]) : -32768);

		for (l = 0; l < 32; l++)
			hmm->trans_vec[7 * hmm->vitQ + q][l] = tmp.v[l];		// since all others are done, DDs are placed at last, #8. so indexing it by "q"..
	}

	/* Transition score bound for "lazy F" DD path evaluation (xref J2/52) */
	hmm->ddbound_w = -32768;

	for (int i = 1; i < (hmm->M - 1) - 1; i++)
	{                                                                 
		ddtmp = (int)wordify(hmm, hmm->log_tran_32bits[i][D_D]);          /* DD */
		ddtmp += (int)wordify(hmm, hmm->log_tran_32bits[i + 2][D_M]);     /* DM */ 
		ddtmp -= (int)wordify(hmm, hmm->log_tran_32bits[i + 2][B_M]);     /* BM */ 
		hmm->ddbound_w = _MAX(hmm->ddbound_w, ddtmp);
	}

	/* Special nodes and parameters for kernel */
	hmm->E_lm = simdlize_VIT((int)wordify(hmm, hmm->Xtran_32bits[E * XTRANS_TYPE + LOOP]));		/* E_LOOP,MOVE same */
	hmm->base_vs = simdlize_VIT(hmm->base_w);													/* base_w */
	hmm->ddbound_vs = simdlize_VIT(hmm->ddbound_w);												/* ddbound */

	/* **************************************************************************************************** */
	/* Until now, we have finished MATCH emission, main node transition, and partly special node transition */
	/* **************************************************************************************************** */

	return fileOK;
}
