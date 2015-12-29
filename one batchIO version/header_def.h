/* #### this is a .cuh file */
#ifndef HMMER_DEFS
#define HMMER_DEFS

/* This is a header file
* It includes all necessary marco value, enum numbers and function prototype
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <vector>

/* all necessary headers */
#include <cuda.h>
#include <nvrtc.h>
#include <helper_functions.h>
#include <cuda_runtime.h>
#include <helper_cuda_drvapi.h>

/* Define */
#define PROTEIN_TYPE 29           /* different types of protein, 20 basic + 1 degenerate (need change if there are more than 1 degen)*/
#define XTRANS_NODE 4             /* N,E,C,J*/
#define XTRANS_TYPE 2             /* Move, Loop*/
#define force_local_size 40		  /* DONT MODIFY! this is for making sure local memory allocation */

// Macro to catch CUDA errors in CUDA runtime calls
#define NVRTC_SAFE_CALL(Name, x)                                             \
  do {                                                                       \
    nvrtcResult result = x;                                                  \
    if (result != NVRTC_SUCCESS) {                                           \
      std::cerr << "\nerror: " << Name << " failed with error " <<           \
                                               nvrtcGetErrorString(result);  \
      exit(1);                                                               \
    }                                                                        \
  } while(0)


#define fileERROR 0
#define fileOK 1

#define _MAX(a,b) (((a) >= (b)) ? (a) : (b))
#define _MIN(a,b) (((a) <= (b)) ? (a) : (b))

#define INFINITY       999999999                    //infinity
#define minusInfinity  logf(0.0f)                   //-infinity
#define eslCONST_LOG2  0.69314718055994529          //ln2
#define F1             0.02                         //default threshold for MSV
#define F2             0.001						//default threshold for VIT

/* log(1+x) ~ x and  1-e^x = -x approximation.
* Same threshold appears to be optimal for float or double x. xref STL9/138.
*/
#define eslSMALLX1    5e-9

/* INDEX ONLY for reading .hmm file */
enum r_trans {
	MM = 0,
	MI = 1,
	MD = 2,
	IM = 3,
	II = 4,
	DM = 5,
	DD = 6
};
#define p7P_NTRANS 7    /* MM,MI,MD,IM,II,DM,DD (order)*/

/* INDEX for configuring transiton matrix */
enum c_trans {
	B_M = 0,
	M_M = 1,
	I_M = 2,
	D_M = 3,
	M_D = 4,
	M_I = 5,
	I_I = 6,
	D_D = 7
};
#define TRANS_TYPE 8     /* BM,MM,DM,IM,MD,MI,II,DD (order)*/

/* Indices for special state types in the length model, gm->xsc[x][]
*/
enum p7p_xstates_e {
	E = 0,
	N = 1,
	J = 2,
	C = 3
};
#define p7P_NXSTATES 4

/* Indices for transitions from the length modeling scores gm->xsc[][x]
*/
enum p7p_xtransitions_e {
	LOOP = 0,
	MOVE = 1
};
#define p7P_NXTRANS 2


/* For the striped implementation */
typedef unsigned int __32uint__[32]; 			//a data type consist of 32 unsigned int values...
typedef int __32int__[32];						/* signed */

#define NQB(X)  ( _MAX(2, ((((X)-1) / 128) + 1)) )   /* 128 uchars (for MSV) */
#define NQW(X)	( _MAX(2, ((((X)-1) / 64 ) + 1)) )    /* 64 uwords (for VIT) */

/**/
typedef struct hmmer_profile_s {

	/* ============================= */
	/* Raw data read from .hmm files */
	/* ============================= */

	float** mat_32bits;									/* match emission score */
	float** ins_32bits;									/* insert emisssion score */
	float** tran_32bits;								/* main mode transition. WHY '7'? Since there is not 'BM' */
	float** log_tran_32bits;							/* log-trans, ncludes 'BM' */
	float   Xtran_32bits[XTRANS_NODE * XTRANS_TYPE];   	/* special node transition: 4 x 2 = 8 */

	/* =============================================================== */
	/* MSVFilter uses scaled, biased uchars: 16x unsigned byte vectors */
	/* =============================================================== */
	int   msvQ;      								/* number of vectors */
	__32uint__ *msv_vec;                            /* striped pointer: each instance is a array of 32 int values */
	__32uint__ *ssv_vec;		/* SSV */

	/* Single Values */
	unsigned int   tbm_b;                              /* constant B->Mk cost:    scaled log 2/M(M+1)       */
	unsigned int   tec_b;                              /* constant E->C  cost:    scaled log 0.5            */
	unsigned int   tjb_b;                              /* constant NCJ move cost: scaled log 3/(L+3)        */
	unsigned int   base_b;                             /* typically +190: offset of uchar scores            */
	unsigned int   bias_b;                             /* positive bias to emission scores, make them >=0   */
	float scale_b;                                     /* typically 3 / log2: scores scale to 1/3 bits      */

	/* Striped Values */
	unsigned int tbm_s;
	unsigned int tec_s;
	unsigned int tjb_s;
	unsigned int base_s;
	unsigned int bias_s;

	/* ================================================================== */
	/* ViterbiFilter uses scaled swords: 8x signed 16-bit integer vectors */
	/* ================================================================== */

	/* This is for SIMD version kernel */
	int vitQ;										/* number of vectors */
	__32int__ *vit_vec;                             /* match emission */
	__32int__ *trans_vec;							/* match transition */

	/* Single Values */
	float   scale_w;                                 /* score units: typically 500 / log(2), 1/500 bits   */
	int     base_w;                                  /* offset of sword scores: typically +12000          */
	int     ddbound_w;                               /* threshold precalculated for lazy DD evaluation    */
	float   ncj_roundoff;                            /* missing precision on NN,CC,JJ after rounding      */

	/* Striped Values */
	int base_vs;									/* offset of sword scores: typically +12000 */
	int E_lm;										/* ONLY FOR H3: E node's loop and move: wordify(-eslCONST_LOG2) */
	int ddbound_vs; 								/* threshold precalculated for lazy DD evaluation */
	int* ncj_move;									/* ONLY FOR H3: related to length of each seq */
	/* ONLY FOR H3: NCJ loop = 0 */

	/* ========================================================= */
	/* Information about current configuration, size, allocation */
	/* ========================================================= */

	int    M;            /* model length                                      */
	double  *MU;          /* parameters for getting P-value                    */
	double  *LAMBDA;      /* parameters for getting P-value                    */
	float  *f;           /* amino acid background frequencies                 */
	float  nj;           /* expected # of J's: 0 or 1, uni vs. multihit       */

} HMMER_PROFILE;

/**/
typedef struct rib_block_per_smx {
	int warps;
	int res_blocks;
} RIB_BLOCK;

/* ------------------------------------------ Function Prototypes -------------------------------------------------------- */

/* degenerate.cpp */
extern int mat_degen(HMMER_PROFILE*);
extern void degen_enable_matrix(int**);
extern int Amino_Offset(char);

/* other_functions.cpp */
extern int NullOne(unsigned int, float*);
extern double esl_gumbel_surv(double, double, double);
extern int p7_AminoFrequencies(float*);
extern void get_block_thread(int*, int*, int, int, int, int);
extern void freeHMM(HMMER_PROFILE*);
extern int round_DIY(float);
extern void check_striped_data(HMMER_PROFILE*);
extern void check_sequential_data(HMMER_PROFILE*);
extern unsigned int seq_Padding(unsigned int**, unsigned int*, char**, unsigned int*, int);
extern void check_1D_array(unsigned int*, unsigned int*, unsigned int*, int*, int);
extern void check_striped_data_vit(HMMER_PROFILE*);
extern void check_sequential_data_vit(HMMER_PROFILE*);
extern void check_striped_trans_vit(HMMER_PROFILE*);
extern void check_sequential_trans_vit(HMMER_PROFILE*);
extern void check_dup_NCJ_move(int*, int);
extern RIB_BLOCK get_opt_MSV(int, int, int);
extern RIB_BLOCK get_opt_VIT(int, int, int);

/* read_files.cpp */
extern int get_Parameters(HMMER_PROFILE*, char*);
extern void nextLine(FILE*, int);
extern void moveCursor(FILE*, int);
extern int get_Emission(HMMER_PROFILE*, char*);
extern int get_transition(HMMER_PROFILE*, char*);
extern int get_Seqnumber(char*);
extern int get_seqID(char*, int*);                          /* Not always suitable in different database */
extern int alloc_Eachseq(char**, unsigned int*, int, char*);
extern int fill_Eachseq(char**, unsigned int*, int, char*);
extern int get_hmm_size(char*);
extern char* read_kernel(char*);

/* model_config.cpp */
extern HMMER_PROFILE* hmmer_profile_Create(int, int);
extern int mf_conversion(HMMER_PROFILE*);
extern unsigned char unbiased_byteify(HMMER_PROFILE*, float);
extern unsigned char biased_byteify(HMMER_PROFILE*, float);
extern int log_Odd_score(HMMER_PROFILE*);
extern int p7_hmm_CalculateOccupancy(HMMER_PROFILE*, float*);
extern int get_entryScore(HMMER_PROFILE*);
extern int log_Trans(HMMER_PROFILE*);
extern int xTrans(HMMER_PROFILE*);
extern short wordify(HMMER_PROFILE*, float);
extern int vf_conversion(HMMER_PROFILE*);
extern unsigned int simdlize_MSV(unsigned int);
extern void simd_NCJ_MOVE(int*, int, HMMER_PROFILE*, unsigned int);

/* launch kernels */
extern void RTC_SSV(unsigned int, const char*, HMMER_PROFILE*,
	     unsigned int*, unsigned int*, unsigned int*,
	     unsigned int*, unsigned int, double*,
	     int, int, dim3, dim3);
extern void RTC_MSV(unsigned int, const char*, HMMER_PROFILE*,
	     unsigned int*, unsigned int*, unsigned int*,
	     unsigned int*, unsigned int, double*,
	     int, int, dim3, dim3);
extern void RTC_VIT(unsigned int, const char*, HMMER_PROFILE*,
	     unsigned int*, unsigned int*, unsigned int*,
	     unsigned int*, unsigned int, double*,
	     int, int, dim3, dim3);


#endif /*HMMER_DEFS*/
