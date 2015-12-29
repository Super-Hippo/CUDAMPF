/*
 * Copyright 2014-2015 Hanyu Jiang. All rights reserved.
 *
 * This is main entrance file by g++
 *
 */
#include "header_def.h"

int main(int argc, char* argv[])
{
	/* Input Parameters */
	char *seqdb_filename;
	char *model_filename;
	model_filename = argv[2];
	seqdb_filename = argv[3];

	/* ************************************************ */
	/* 		            0. Preparation		            */
	/* ************************************************ */

	/* Get IMPORTANT parameters */
	int number = 0;								/* the total number of seq in this database */
	int model_size = 0; 						/* get hmm model size at first */
	number = get_Seqnumber(seqdb_filename);
	model_size = get_hmm_size(model_filename);

	printf("In this database, we have %d seqs.\n", number);
	printf("The HMM size is %d.\n", model_size);

	/* Initiate hmm model */
	HMMER_PROFILE *hmm = NULL;							/* claim hmm         */
	hmm = hmmer_profile_Create(number, model_size);		/* alloc mem for hmm */

	/* ************************************************ */
	/* 		1. Parameters, match, insert emission 		*/
	/* ************************************************ */
	StopWatchInterface *timer;
    sdkCreateTimer(&timer);
    sdkStartTimer(&timer);

	/* get parameters (M, MU[], LAMBDA[]) */
	if (get_Parameters(hmm, model_filename) != 1) printf("error!\n");

	/* Get raw match & insert emission */
	if (get_Emission(hmm, model_filename) != 1) printf("error!\n");

	/* Transfer to prob from raw data (need do it before get degen)*/
	if (log_Odd_score(hmm) != 1) printf("error!\n");

	/* Get degenerate residues */
	if (mat_degen(hmm) != 1) printf("error in getiing 'degenerate residue'! \n");

	/* Rescaling and rounding */
	if (mf_conversion(hmm) != 1) printf("error!\n");

	//Until now, all preparation works of MSV has been done!


	/* ************************************************************************ */
	/* 		2. Transition probability only for Viterbi filter 		*/
	/* ************************************************************************ */

	/**/
	if (get_transition(hmm, model_filename) != 1) printf("error!\n");

	/**/
	if (get_entryScore(hmm) != 1) printf("error!\n");

	/**/
	if (log_Trans(hmm) != 1) printf("error!\n");

	/**/
	if (xTrans(hmm) != 1) printf("error!\n");

	/**/
	if (vf_conversion(hmm) != 1) printf("error!\n");

	//Until now, all preparation works of Viterbi has been done!

	sdkStopTimer(&timer);
    printf("model process time: %f (ms)\n", sdkGetTimerValue(&timer));
    sdkDeleteTimer(&timer);

	/* **************************************************************************** */
	/* 		3. Protein sequence database read and pre-process               		*/
	/* **************************************************************************** */

	sdkCreateTimer(&timer);
    sdkStartTimer(&timer);

	/* This part is the simplest implementation without any decent speedup optimization */
	char** seq = (char**)malloc(number * sizeof(char*));     							/* dynamic memory for address value of each sequence: seq[i]->an address value */
	unsigned int* seq_len = (unsigned int*)malloc(number * sizeof(unsigned int));       /* for cache length of each sequence */
	if (alloc_Eachseq(seq, seq_len, number, seqdb_filename) != 1) printf("error!\n");
	if (fill_Eachseq(seq, seq_len, number, seqdb_filename) != 1)  printf("error!\n");

	/* Padding and digitalize */
	unsigned int sum = 0;
	unsigned int** iSeq = (unsigned int**)malloc(number * sizeof(unsigned int*)); /* 4 uchar in 1 */
	unsigned int* iLen = (unsigned int*)malloc(number * sizeof(unsigned int));		/* length of 4 in 1*/
	sum = seq_Padding(iSeq, iLen, seq, seq_len, number);

	/* imm check sum */
	if ((sum % 32 == 0) && (sum % 4 == 0)) {
		printf("the seq database has been aligned!!\n");
	}
	else {
		printf("padding sequences unaligned!\n");
		getchar();
		exit(0);
	}

	/* transfer iSeq from 2D to 1D */
	unsigned int* seq_1D = NULL;
	unsigned int* offset = NULL;
	seq_1D = (unsigned int*)malloc(sum * sizeof(unsigned int));
	offset = (unsigned int*)malloc(number * sizeof(unsigned int));
	memset(seq_1D, 0, sum * sizeof(unsigned int));
	memset(offset, 0, number * sizeof(unsigned int));

	/* fill the array (temporary solution without any optimization, this is painful) */
	unsigned int addr = 0;
	unsigned int i, j;
	for (i = 0; i < number; i++)
	{
		offset[i] = addr;				/* one addr is 1 "unsigned int" value consist of 4 uchars */

		for (j = 0; j < iLen[i]; j++)
		{
			seq_1D[addr] = iSeq[i][j];
			addr += 1;					/* addr is just a index */
		}
	}

	/* Until now, we have 1D array: iSeq and the offset for it */

	sdkStopTimer(&timer);
    printf("whole seq database process time: %f (ms)\n", sdkGetTimerValue(&timer));
    sdkDeleteTimer(&timer);

	/* ******************************************************************** */
	/* 			4. ALL COPY AND LAUNCH BY USING CUDA DRIVE API				*/
	/*						NEW THINGS BEGIN HERE  							*/
	/* ******************************************************************** */

	/* Get the Device Property */
	CUdevice cuDevice;
	checkCudaErrors(cuInit(0));
	checkCudaErrors(cuDeviceGet(&cuDevice, 0));

	/* Device Property: fixed based on Device */
	int WARP_SIZE;
	cuDeviceGetAttribute(&WARP_SIZE, CU_DEVICE_ATTRIBUTE_WARP_SIZE, cuDevice);
	int SMX;
	cuDeviceGetAttribute(&SMX, CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT, cuDevice);
	int MAX_SMEM;
	cuDeviceGetAttribute(&MAX_SMEM, CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_MULTIPROCESSOR, cuDevice);

	printf("WARP_SIZE = %d\n", WARP_SIZE);
	printf("SMX = %d\n", SMX);
	printf("MAX_SMEM = %d\n", MAX_SMEM);

    char* handle = NULL;	/* pointer to kernel file */
    RIB_BLOCK rb;			/* get optimal RIB for SMEM */
	int opt_Reg = 0;		/* get maximum available regs for each threads based on resident registers per SMX on current context */

    double *pValue;		/* score got from device for each seq */
	pValue = (double*)malloc(number * sizeof(double));
	memset(pValue, 0, number * sizeof(double));
        
    int defect_RIB = 0;		/*This is for making up the defect of nvrtc compilor to avoid weird OUT OF RESOURCE */

    /* switchable grid & block */
    dim3 GRID;
    dim3 BLOCK;

	/* Switchable Launch */
	switch(atoi(argv[1]))
	{
		case 1:	/* LMEM MSV only */
				handle = read_kernel("LMEM_MSV.cuh");	/* read LMEM */
				opt_Reg = 64;

				GRID = dim3(1, SMX, 1);
				BLOCK = dim3(WARP_SIZE, 32, 1);		/* [FIXED] one block per SMX, 1024 threads per block */

				RTC_MSV(number, handle, hmm,
					    seq_1D, offset, seq_len,
					    iLen, sum, pValue,
					    32, opt_Reg, GRID, BLOCK);	/* 32 warps per block; 64 reg per thread; */
				break;

		case 2:	/* SMEM MSV only */
    			rb = get_opt_MSV(hmm->msvQ, MAX_SMEM, WARP_SIZE);			/* scalable for different device */
    			printf("SMEM MSV::warps per block: %d\n", rb.warps);
    			printf("SMEM MSV::resident blocks per SMX: %d\n", rb.res_blocks);

    			handle = read_kernel("SMEM_MSV.cuh");	/* read SMEM */

    			/* 65536 is hard-coded (need change to device property) */
    			opt_Reg = (65536/(WARP_SIZE * rb.warps * rb.res_blocks) > 255) ? 255 : (65536/(WARP_SIZE * rb.warps * rb.res_blocks));
    			printf("Available registers per thread is: %d\n", opt_Reg);

    			GRID = dim3(1, rb.res_blocks * SMX, 1);
      			BLOCK = dim3(WARP_SIZE, rb.warps, 1);		/* RIB & Resident blocks always change depending size of model for optimal */

				RTC_MSV(number, handle, hmm,
					    seq_1D, offset, seq_len,
					    iLen, sum, pValue,
					    rb.warps, opt_Reg, GRID, BLOCK);			/* optimal warps per block; 32 reg per thread (based on test); */		
				break;

		case 3: /* LMEM VIT only */
				handle = read_kernel("LMEM_VIT.cuh");	/* read LMEM */
				opt_Reg = 64;

				GRID = dim3(1, SMX, 1);
				BLOCK = dim3(WARP_SIZE, 32, 1);	/* [FIXED] one block per SMX, 1024 threads per block */

				RTC_VIT(number, handle, hmm,
					    seq_1D, offset, seq_len,
					    iLen, sum, pValue,
					    32, opt_Reg, GRID, BLOCK);	/* 32 warps per block; 64 reg per thread; */

				break;

		case 4: /* SMEM VIT only */
				rb = get_opt_VIT(hmm->vitQ, MAX_SMEM, WARP_SIZE);			/* scalable for different device */
    			printf("SMEM VIT::warps per block: %d\n", rb.warps);
    			printf("SMEM VIT::resident blocks per SMX: %d\n", rb.res_blocks);

    			handle = read_kernel("SMEM_VIT.cuh");	/* read SMEM */

    			/* 65536 is hard-coded (need change to device property) */
    			opt_Reg = (65536/(WARP_SIZE * rb.warps * rb.res_blocks) > 255) ? 255 : (65536/(WARP_SIZE * rb.warps * rb.res_blocks));
    			printf("Available registers per thread is: %d\n", opt_Reg);

    			GRID = dim3(1, rb.res_blocks * SMX, 1);
      			BLOCK = dim3(WARP_SIZE, rb.warps, 1);		/* RIB & Resident blocks always change depending size of model for optimal */

				RTC_VIT(number, handle, hmm,
					    seq_1D, offset, seq_len,
					    iLen, sum, pValue,
					    rb.warps, opt_Reg, GRID, BLOCK);			/* optimal warps per block; 32 reg per thread (based on test); */	
				break;

		case 5: /* LMEM SSV only */
				handle = read_kernel("LMEM_SSV.cuh");	/* read LMEM */
				opt_Reg = 64;

				if ((hmm->msvQ > 15)&&(hmm->msvQ < 19))
				  defect_RIB = 28;
				else
				  defect_RIB = 32;
				GRID = dim3(1, SMX, 1);
				BLOCK = dim3(WARP_SIZE, defect_RIB, 1);		/* [FIXED] one block per SMX, 1024 threads per block */

				RTC_SSV(number, handle, hmm,
					    seq_1D, offset, seq_len,
					    iLen, sum, pValue,
					    defect_RIB, opt_Reg, GRID, BLOCK);	/* 32 warps per block; 64 reg per thread; */
				break;

		case 6:	/* SMEM SSV only */
    			rb = get_opt_MSV(hmm->msvQ, MAX_SMEM, WARP_SIZE);			/* scalable for different device */
    			printf("SMEM SSV::warps per block: %d\n", rb.warps);
    			printf("SMEM SSV::resident blocks per SMX: %d\n", rb.res_blocks);

    			handle = read_kernel("SMEM_SSV.cuh");	/* read SMEM */

    			/* 65536 is hard-coded (need change to device property) */
    			opt_Reg = (65536/(WARP_SIZE * rb.warps * rb.res_blocks) > 255) ? 255 : (65536/(WARP_SIZE * rb.warps * rb.res_blocks));
    			printf("Available registers per thread is: %d\n", opt_Reg);

    			GRID = dim3(1, rb.res_blocks * SMX, 1);
      			BLOCK = dim3(WARP_SIZE, rb.warps, 1);		/* RIB & Resident blocks always change depending size of model for optimal */

				RTC_SSV(number, handle, hmm,
					    seq_1D, offset, seq_len,
					    iLen, sum, pValue,
					    rb.warps, opt_Reg, GRID, BLOCK);			/* optimal warps per block; 32 reg per thread (based on test); */		
				break;
		case 7:
			if(hmm->M <= 300)
			{
				rb = get_opt_MSV(hmm->msvQ, MAX_SMEM, WARP_SIZE);
    				printf("SMEM MSV::warps per block: %d\n", rb.warps);
    				printf("SMEM MSV::resident blocks per SMX: %d\n", rb.res_blocks);
    				handle = read_kernel("SMEM_MSV.cuh");
    				opt_Reg = (65536/(WARP_SIZE * rb.warps * rb.res_blocks) > 255) ? 255 : (65536/(WARP_SIZE * rb.warps * rb.res_blocks));
    				printf("Available registers per thread is: %d\n", opt_Reg);
    				GRID = dim3(1, rb.res_blocks * SMX, 1);
      				BLOCK = dim3(WARP_SIZE, rb.warps, 1);
				RTC_MSV(number, handle, hmm,
					seq_1D, offset, seq_len,
					iLen, sum, pValue,
					rb.warps, opt_Reg, GRID, BLOCK);
			} else {
				handle = read_kernel("LMEM_MSV.cuh");
				opt_Reg = 64;
				GRID = dim3(1, SMX, 1);
				BLOCK = dim3(WARP_SIZE, 32, 1);
				RTC_MSV(number, handle, hmm,
					seq_1D, offset, seq_len,
					iLen, sum, pValue,
					32, opt_Reg, GRID, BLOCK);
			}
			break;
		case 8:
			if(hmm->M <= 200)
			{
				rb = get_opt_VIT(hmm->vitQ, MAX_SMEM, WARP_SIZE);
	    			printf("SMEM VIT::warps per block: %d\n", rb.warps);
	    			printf("SMEM VIT::resident blocks per SMX: %d\n", rb.res_blocks);
	    			handle = read_kernel("SMEM_VIT.cuh");
	    			opt_Reg = (65536/(WARP_SIZE * rb.warps * rb.res_blocks) > 255) ? 255 : (65536/(WARP_SIZE * rb.warps * rb.res_blocks));
	    			printf("Available registers per thread is: %d\n", opt_Reg);
	    			GRID = dim3(1, rb.res_blocks * SMX, 1);
	      			BLOCK = dim3(WARP_SIZE, rb.warps, 1);
				RTC_VIT(number, handle, hmm,
					seq_1D, offset, seq_len,
					iLen, sum, pValue,
					rb.warps, opt_Reg, GRID, BLOCK);
			} else {
				handle = read_kernel("LMEM_VIT.cuh");
				opt_Reg = 64;
				GRID = dim3(1, SMX, 1);
				BLOCK = dim3(WARP_SIZE, 32, 1);
				RTC_VIT(number, handle, hmm,
					seq_1D, offset, seq_len,
					iLen, sum, pValue,
					32, opt_Reg, GRID, BLOCK);
			}
			break;
		case 9:
			if(hmm->M <= 500)
			{
   				rb = get_opt_MSV(hmm->msvQ, MAX_SMEM, WARP_SIZE);			/* scalable for different device */
    			printf("SMEM SSV::warps per block: %d\n", rb.warps);
    			printf("SMEM SSV::resident blocks per SMX: %d\n", rb.res_blocks);

    			handle = read_kernel("SMEM_SSV.cuh");	/* read SMEM */

    			/* 65536 is hard-coded (need change to device property) */
    			opt_Reg = (65536/(WARP_SIZE * rb.warps * rb.res_blocks) > 255) ? 255 : (65536/(WARP_SIZE * rb.warps * rb.res_blocks));
    			printf("Available registers per thread is: %d\n", opt_Reg);

    			GRID = dim3(1, rb.res_blocks * SMX, 1);
      			BLOCK = dim3(WARP_SIZE, rb.warps, 1);		/* RIB & Resident blocks always change depending size of model for optimal */

				RTC_SSV(number, handle, hmm,
					    seq_1D, offset, seq_len,
					    iLen, sum, pValue,
					    rb.warps, opt_Reg, GRID, BLOCK);	/* optimal warps per block; 32 reg per thread (based on test); */	
			} else {
				handle = read_kernel("LMEM_SSV.cuh");		/* read LMEM */
				opt_Reg = 64;

				if ((hmm->msvQ > 15)&&(hmm->msvQ < 19))
				  defect_RIB = 28;
				else
				  defect_RIB = 32;
				GRID = dim3(1, SMX, 1);
				BLOCK = dim3(WARP_SIZE, defect_RIB, 1);		/* [FIXED] one block per SMX, 1024 threads per block */

				RTC_SSV(number, handle, hmm,
					    seq_1D, offset, seq_len,
					    iLen, sum, pValue,
					    defect_RIB, opt_Reg, GRID, BLOCK);	/* 32 warps per block; 64 reg per thread; */
			}
			break;
		default:
				break;
	}

	/* ***************************************** */
	/* 		  6. Host & Device memory release  	 */
	/* ***************************************** */
#if 0
	for (int i = 0; i < number; i++) {
		free(seq[i]);
	}
	free(seq);				/* PLEASE REMIND: free child first, then parent */
	printf("Host free: seq...\n");

	for (int i = 0; i < number; i++) {
		free(iSeq[i]);
	}
	free(iSeq);				/* PLEASE REMIND: free child first, then parent */
	printf("Host free: iSeq...\n");

	free(iLen);
	free(seq_len);
	free(seq_1D);
	free(offset);
	free(sc);
	free(nullsc);

	/* HMM model */
	freeHMM(hmm);
#endif

	return 0;
}
