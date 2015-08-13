/* For reading useful data from .hmm files
* Now only focus on the MSV filter
*
*
*
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sstream>
#include <iostream>
#include <fstream>

#include "header_def.h"

static unsigned long total_residue = 0;
FILE* FilePointer = NULL;

/* read length (temprary solution) */
int get_hmm_size(char* hmm_Path)
{
	char* LEN = "LENG";
	FilePointer = fopen(hmm_Path, "r");
	char* tag_len = (char*)malloc(5 * sizeof(char));    //read width of 5 bytes
	char* lenTok = (char*)malloc(100 * sizeof(char));   //enough width of 100 bytes to ensure get accurate 'length'
	int length = 0;

	if (FilePointer == NULL) {
		printf("Fatal error: Cannot open or find <.hmm> file\n");
		exit(-1);
	}
	else {
		printf("Got hmm length\n");
	}

	/* get the size of hmm */
	do{
		fgets(tag_len, 5, FilePointer);
		if (strcmp(tag_len, LEN))
			nextLine(FilePointer, 1);
	} while (strcmp(tag_len, LEN) && !feof(FilePointer));

	fgets(lenTok, 100, FilePointer);
	length = (int)atof(lenTok);       //reading raw .hmm file and get the length of it

	free(tag_len);
	free(lenTok);
	fclose(FilePointer);

	return length;
}

/* Function:  get_Parameters()
* Synopsis:  Read .hmm file and get specified data.
*
* Purpose:  We got several important parameters here:
*           1. hmm->M
*           2. hmm->MU[0,1,2]
*           3. hmm->LAMBDA[0,1,2]
* Returns:
*
*/
int get_Parameters(HMMER_PROFILE *hmm, char* hmm_Path)
{
	//int BITS;
	//printf("%d", sizeof(BITS));
	//getchar();

	char* LEN = "LENG";
	char* PARA_MSV = "STATS LOCAL MSV";
	char* PARA_VIT = "STATS LOCAL VITERBI";   //Ã»×ö
	char* PARA_FWD = "STATS LOCAL FORWARD";

	FilePointer = fopen(hmm_Path, "r");
	if (FilePointer == NULL) {
		printf("Fatal error: Cannot open or find <.hmm> file\n");
		getchar();
		return fileERROR;
	}
	else {
		printf(".hmm file open well\n");
	}

	/* 2. get MSV parameters */
	char* tag_MSV = (char*)malloc(16 * sizeof(char)); //CAUTION: 16 is a dead data
	char* MSVTok = (char*)malloc(100 * sizeof(char));
	char *p_msv;
	do{
		fgets(tag_MSV, 16, FilePointer);
		if (strcmp(tag_MSV, PARA_MSV))
			nextLine(FilePointer, 1);
	} while (strcmp(tag_MSV, PARA_MSV) && !feof(FilePointer));
	fgets(MSVTok, 100, FilePointer);

	p_msv = strtok(MSVTok, "  ");
	hmm->MU[0] = (float)atof(p_msv);
	p_msv = strtok(NULL, "  ");
	hmm->LAMBDA[0] = (float)atof(p_msv);

	/* demo for extracting values from a string
	* ==============================================
	* char* test1 = (char*)malloc(100*sizeof(char));
	* fgets(test1, 100, FilePointer);
	* char *p;
	* p = strtok(test1, "  ");
	* while(p)
	* {
	*  printf("%s\n", p);
	*  p = strtok(NULL, "  ");
	* }
	* ==============================================
	*/

	/* 3. get Viterbi parameters */
	char* tag_VIT = (char*)malloc(20 * sizeof(char));
	char* VITTok = (char*)malloc(100 * sizeof(char));
	char *p_vit;
	do{
		fgets(tag_VIT, 20, FilePointer);
		if (strcmp(tag_VIT, PARA_VIT))
			nextLine(FilePointer, 1);
	} while (strcmp(tag_VIT, PARA_VIT) && !feof(FilePointer));
	fgets(VITTok, 100, FilePointer);

	p_vit = strtok(VITTok, "  ");
	hmm->MU[1] = (float)atof(p_vit);
	p_vit = strtok(NULL, "  ");
	hmm->LAMBDA[1] = (float)atof(p_vit);

	/* 4. get Forward parameters */
	char* tag_FWD = (char*)malloc(20 * sizeof(char));
	char* FWDTok = (char*)malloc(100 * sizeof(char));
	char *p_fwd;
	do{
		fgets(tag_FWD, 20, FilePointer);
		if (strcmp(tag_FWD, PARA_FWD))
			nextLine(FilePointer, 1);
	} while (strcmp(tag_FWD, PARA_FWD) && !feof(FilePointer));
	fgets(FWDTok, 100, FilePointer);

	p_fwd = strtok(FWDTok, "  ");
	hmm->MU[2] = (float)atof(p_fwd);
	p_fwd = strtok(NULL, "  ");
	hmm->LAMBDA[2] = (float)atof(p_fwd);

	printf("MSV MU: %f, LAMBDA: %f\n", hmm->MU[0], hmm->LAMBDA[0]);
	printf("viterbi MU: %f, LAMBDA: %f\n", hmm->MU[1], hmm->LAMBDA[1]);
	printf("Forward MU: %f, LAMBDA: %f\n", hmm->MU[2], hmm->LAMBDA[2]);

	/* 5. free */
	fclose(FilePointer);
	free(tag_MSV);
	free(MSVTok);
	free(tag_FWD);
	free(FWDTok);

	return fileOK;
}


/* always move build-in pointer before FIRST character in next 'times'line */
void nextLine(FILE* pointer_2, int times)
{
	int moveChar;
	for (int i = 0; i < times; i++)
	{
		moveChar = 0;
		while ((char)moveChar != /*CAUTION*/ '\n')       //condition is based on Windows. It will be different in Linux and Mac OS
		{
			moveChar = fgetc(pointer_2);
		}
	}
}

/* move file curser back with 'times' chars */
void moveCursor(FILE* pointer_2, int times)
{
	for (int i = 0; i < times; i++)
	{
		fgetc(pointer_2);
	}
}

/* Function:  get_MatchEmission()
* Synopsis:  Read .hmm file and get specified data.
*
*/
int get_Emission(HMMER_PROFILE *hmm, char* hmm_Path)
{
	char* BEGIN = "  COMPO   ";
	int i = 0;
	int j = 0;
	int q = 0;

	FilePointer = fopen(hmm_Path, "r");
	if (FilePointer == NULL) {
		printf("Fatal error: Cannot open or find <.hmm> file\n");
		getchar();
		return fileERROR;
	}
	else {
		printf(".hmm file open well\n");
	}

	/* 1. locate the 'COMPO' */
	char* locate = (char*)malloc(11 * sizeof(char));
	do{
		fgets(locate, 11, FilePointer);
		if (strcmp(locate, BEGIN))
			nextLine(FilePointer, 1);
	} while (strcmp(locate, BEGIN) && !feof(FilePointer));

	/* 2. Process node 0 for I (no M here) */
	char* m_Temp = (char*)malloc(179 * sizeof(char));
	char* match;

	char* i_Temp = (char*)malloc(179 * sizeof(char));
	char* insert;

	nextLine(FilePointer, 1);
	moveCursor(FilePointer, 10);

	fgets(i_Temp, 179, FilePointer);
	insert = strtok(i_Temp, "  ");
	while (insert) {
		//printf("%s\n", insert);
		hmm->ins_32bits[0][i] = expf(-1.0 * (float)atof(insert));
		i++;
		insert = strtok(NULL, "  ");
	}

	nextLine(FilePointer, 2);
	moveCursor(FilePointer, 10);

	/* 3. Process node 1 to M (both I and M) */
	j = j + 1;

	do{
		//Match Emission
		i = 0;
		fgets(m_Temp, 179, FilePointer);
		match = strtok(m_Temp, "  ");
		while (match) {
			hmm->mat_32bits[j][i] = expf(-1.0 * (float)atof(match));
			i++;
			match = strtok(NULL, "  ");
		}

		i = 0;
		nextLine(FilePointer, 1);
		moveCursor(FilePointer, 10);

		fgets(i_Temp, 179, FilePointer);
		insert = strtok(i_Temp, "  ");
		while (insert) {
			hmm->ins_32bits[j][i] = expf(-1.0 * (float)atof(insert));
			i++;
			insert = strtok(NULL, "  ");
		}

		//Finish one node..
		nextLine(FilePointer, 2);
		moveCursor(FilePointer, 10);
		j++;

	} while (j <= hmm->M);      //ensure that we can fill-up Mat/Ins[HMM_SIZE + 1][PROTEIN_TYPE]


	/* 4. free */
	fclose(FilePointer);
	free(locate);
	free(m_Temp);
	free(i_Temp);

	return fileOK;
}

int get_transition(HMMER_PROFILE *hmm, char* hmm_Path)
{
	int i;

	char* BEGIN = "  COMPO   ";
	FilePointer = fopen(hmm_Path, "r");
	if (FilePointer == NULL) {
		printf("Fatal error: Cannot open or find <.hmm> file\n");
		getchar();
		return fileERROR;
	}
	else {
		printf(".hmm file open well_trans\n");
	}

	/* 1. locate the line of 1st stage */
	char* locate = (char*)malloc(11 * sizeof(char));
	do{
		fgets(locate, 11, FilePointer);
		if (strcmp(locate, BEGIN))
			nextLine(FilePointer, 1);
	} while (strcmp(locate, BEGIN) && !feof(FilePointer));

	nextLine(FilePointer, 2);

	char* t_Temp = (char*)malloc(72 * sizeof(char));
	char* p;

	/* 2. Process node 0 */
	fgets(t_Temp, 72, FilePointer);

	p = strtok(t_Temp, "  ");
	hmm->tran_32bits[0][MM] = expf(-1.0 * (float)atof(p));
	p = strtok(NULL, "  ");
	hmm->tran_32bits[0][MI] = expf(-1.0 * (float)atof(p));
	p = strtok(NULL, "  ");
	hmm->tran_32bits[0][MD] = expf(-1.0 * (float)atof(p));
	p = strtok(NULL, "  ");
	hmm->tran_32bits[0][IM] = expf(-1.0 * (float)atof(p));
	p = strtok(NULL, "  ");
	hmm->tran_32bits[0][II] = expf(-1.0 * (float)atof(p));
	p = strtok(NULL, "  ");
	hmm->tran_32bits[0][DM] = 1.0f;

	hmm->tran_32bits[0][DD] = 0;

	nextLine(FilePointer, 3);             //move to next 'transition' line... 

	/* 3. Process node 1 to M-1 */
	for (i = 1; i < hmm->M; i++) {
		fgets(t_Temp, 72, FilePointer);

		p = strtok(t_Temp, "  ");
		hmm->tran_32bits[i][MM] = expf(-1.0 * (float)atof(p));
		p = strtok(NULL, "  ");
		hmm->tran_32bits[i][MI] = expf(-1.0 * (float)atof(p));
		p = strtok(NULL, "  ");
		hmm->tran_32bits[i][MD] = expf(-1.0 * (float)atof(p));
		p = strtok(NULL, "  ");
		hmm->tran_32bits[i][IM] = expf(-1.0 * (float)atof(p));
		p = strtok(NULL, "  ");
		hmm->tran_32bits[i][II] = expf(-1.0 * (float)atof(p));
		p = strtok(NULL, "  ");
		hmm->tran_32bits[i][DM] = expf(-1.0 * (float)atof(p));
		p = strtok(NULL, "  ");
		hmm->tran_32bits[i][DD] = expf(-1.0 * (float)atof(p));

		nextLine(FilePointer, 3);
	}

	/* 4. Process node M */
	fgets(t_Temp, 72, FilePointer);

	p = strtok(t_Temp, "  ");
	hmm->tran_32bits[hmm->M][MM] = expf(-1.0 * (float)atof(p));
	p = strtok(NULL, "  ");
	hmm->tran_32bits[hmm->M][MI] = expf(-1.0 * (float)atof(p));

	hmm->tran_32bits[hmm->M][MD] = expf(-1.0 * (float)atof(p));

	p = strtok(NULL, "        *  ");
	hmm->tran_32bits[hmm->M][IM] = expf(-1.0 * (float)atof(p));
	p = strtok(NULL, "  ");
	hmm->tran_32bits[hmm->M][II] = expf(-1.0 * (float)atof(p));
	p = strtok(NULL, "  ");

	hmm->tran_32bits[hmm->M][DM] = 1.0f;

	hmm->tran_32bits[hmm->M][DD] = 0;

	/* 3. free memory */
	fclose(FilePointer);
	free(locate);
	free(t_Temp);

	return fileOK;
}


/* Function:  get_Seqnumber() */
int get_Seqnumber(char* seq_Path)
{
	int fileChar = 0;
	int times = 0;

	FilePointer = fopen(seq_Path, "r");
	if (FilePointer == NULL) {
		printf("Fatal error: Cannot open or find <.fsa> file\n");
		getchar();
		return fileERROR;
	}
	else {
		printf(".fsa file open well\n");
	}

	fileChar = fgetc(FilePointer);
	while (!feof(FilePointer)) {
		if ((char)fileChar == '>') {
			times++;
			nextLine(FilePointer, 1);
		}
		fileChar = fgetc(FilePointer);
	}

	fclose(FilePointer);

	return times;
}

/* dynamic memory allocation for each sequence on host */
int alloc_Eachseq(char** addr, unsigned int* len_, int num, char* seq_Path)
{
	int fileChar = 0;
	int index = 0;

	FilePointer = fopen(seq_Path, "r");
	if (FilePointer == NULL) {
		printf("Fatal error: Cannot open or find <.fsa> file\n");
		getchar();
		return fileERROR;
	}
	else {
		printf(".fsa file open well\n");
	}

	fileChar = fgetc(FilePointer);
	while (!feof(FilePointer)) {
		if ((char)fileChar == '>') {
			nextLine(FilePointer, 1);
			int count = 0;
			int fileChar_in = 0;

			while ((char)fileChar_in != '>' && fileChar_in != EOF) {
				fileChar_in = fgetc(FilePointer);
				if (fileChar_in >= 65 && fileChar_in <= 90) {
					count++;
				}
			}

			if (index < num) {
				addr[index] = (char*)malloc(count*sizeof(char));
				len_[index] = count;
				total_residue += len_[index];
				index++;
			}
		}
	}

	fclose(FilePointer);
	printf("\n Total residue = %ld\n", total_residue);
	
	return fileOK;
}

/* cache each sequence into corresponding position */
int fill_Eachseq(char** addr, unsigned int* len_, int num, char* seq_Path)
{
	int fileChar = 0;
	int index = 0;

	FilePointer = fopen(seq_Path, "r");
	if (FilePointer == NULL) {
		printf("Fatal error: Cannot open or find <.fsa> file\n");
		getchar();
		return fileERROR;
	}
	else {
		printf(".fsa file open well\n");
	}

	fileChar = fgetc(FilePointer);
	while (!feof(FilePointer)) {
		if ((char)fileChar == '>') {
			nextLine(FilePointer, 1);
			int i = 0;
			int fileChar_in = 0;

			while ((char)fileChar_in != '>' && fileChar_in != EOF) {
				fileChar_in = fgetc(FilePointer);
				if (fileChar_in >= 65 && fileChar_in <= 90 && i < len_[index]) {
					addr[index][i] = (char)fileChar_in;
					i++;
				}
			}

			if (index < num) {
				index++;
			}
		}
	}

	fclose(FilePointer);

	return fileOK;
}

// read the .cu file and convert kernel content to char* array
char* read_kernel(char *filename)
{
	StopWatchInterface *timer;
    sdkCreateTimer(&timer);
    sdkStartTimer(&timer);

    std::ifstream inputFile(filename, std::ios::in | std::ios::binary |
                                std::ios::ate);

    if (!inputFile.is_open()) 
    {
        std::cerr << "\nerror: unable to open " << filename << " for reading!\n";
        exit(1);
    }

    std::streampos pos = inputFile.tellg();
    size_t inputSize = (size_t)pos;
    char * memBlock = new char [inputSize + 1];

    inputFile.seekg (0, std::ios::beg);
    inputFile.read (memBlock, inputSize);
    inputFile.close();
    memBlock[inputSize] = '\x0';

    sdkStopTimer(&timer);
    printf("read .cuh file time: %f (ms)\n", sdkGetTimerValue(&timer));
    sdkDeleteTimer(&timer);

    return memBlock;
}




