#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "header_def.h"

/* Add degenerate residue */
int mat_degen(HMMER_PROFILE *hmm)
{
	/* For each node, we calculate extra degenerate emissions,
	* and we have 6 degens with value, and other 3 no-value
	* degens.
	* According to different database(amino, DNA, RNA), those
	* degens represents different meaning, and certainly, their
	* emission score are completely different.
	* Here, we only consider about the 'amino' case
	*/

	/* Basic:              A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y                */
	int degen[6][20] = {
		{ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },  /* B => ND */
		{ 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },  /* J => IL */
		{ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 },  /* Z => QE */
		{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },  /* O => K  */
		{ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },  /* U => C  */
		{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }   /* X =>    */
	};
	/* Ö»ÓÐÔÚÊý×é¶¨ÒåµÄÊ±ºò²Å¿ÉÒÔÈç´Ë³õÊ¼»¯£¬ÏÈ¶¨Òå£¬ºó³õÊ¼»¯£¬Ôò¸Ã·¨²»¿ÉÐÐ*/

	float result = 0.;
	float denom = 0.;
	int p = 20;                                 /* 20 is beginning of degen: 'B' */

	/* Matrix */
	for (int i = 1; i <= hmm->M; i++)            /* 1 to 1901 */
	{
		result = 0.;
		denom = 0.;
		p = 20;

		for (int q = 0; q < 6; q++)             /*   0 to 5  */
		{
			result = 0.;
			denom = 0.;

			for (int j = 0; j < 20; j++)      /*  0 to 19  */
			{
				if (degen[q][j] == 1)
				{
					result += hmm->mat_32bits[i][j] * hmm->f[j];
					denom += hmm->f[j];
				}
			}
			result /= denom;
			hmm->mat_32bits[i][p] = result;
			p = p + 1;                        /* move to next degen */
		}

		/* For those 3 non-valued residue */
		hmm->mat_32bits[i][26] = minusInfinity;     /* Gap */
		hmm->mat_32bits[i][27] = minusInfinity;     /* nonresidue */
		hmm->mat_32bits[i][28] = minusInfinity;     /* missing data */
	}

	return 1;
}


/* ------------------------------------------- */
/*         Complete degenerate transfer        */
/* ------------------------------------------- */
/* include 20 standard CAPITAL residues.       */
/* include 6 degenerate residues with meaning. */
/* include 3 no meaningful residue.            */
/* allow low-case letter.                      */
/* totally 29 different residues               */
/* ------------------------------------------- */
int Amino_Offset(char amino)
{
	switch ((int)amino)
	{
	case 65:  return 0;  /* A */
	case 67:  return 1;  /* C */
	case 68:  return 2;  /* D */
	case 69:  return 3;  /* E */
	case 70:  return 4;  /* F */
	case 71:  return 5;  /* G */
	case 72:  return 6;  /* H */
	case 73:  return 7;  /* I */
	case 75:  return 8;  /* K */
	case 76:  return 9;  /* L */
	case 77:  return 10; /* M */
	case 78:  return 11; /* N */
	case 80:  return 12; /* P */
	case 81:  return 13; /* Q */
	case 82:  return 14; /* R */
	case 83:  return 15; /* S */
	case 84:  return 16; /* T */
	case 86:  return 17; /* V */
	case 87:  return 18; /* W */
	case 89:  return 19; /* Y */

		/*  6 degenerate residues  */
	case 66:  return 20; /* B */
	case 74:  return 21; /* J */
	case 90:  return 22; /* Z */
	case 79:  return 23; /* O */
	case 85:  return 24; /* U */
	case 88:  return 25; /* X */

		/* Gap, nonresidue, missing data */
	case 45:  return 26; /* - */
	case 42:  return 27; /* * */
	case 126: return 28; /* ~ */

		/* low-case & mutated gap */
	case 97:  return 0;  /* a */
	case 99:  return 1;  /* c */
	case 100: return 2;  /* d */
	case 101: return 3;  /* e */
	case 102: return 4;  /* f */
	case 103: return 5;  /* g */
	case 104: return 6;  /* h */
	case 105: return 7;  /* i */
	case 107: return 8;  /* k */
	case 108: return 9;  /* l */
	case 109: return 10; /* m */
	case 110: return 11; /* n */
	case 112: return 12; /* p */
	case 113: return 13; /* q */
	case 114: return 14; /* r */
	case 115: return 15; /* s */
	case 116: return 16; /* t */
	case 118: return 17; /* v */
	case 119: return 18; /* w */
	case 121: return 19; /* y */
	case 98:  return 20; /* b */
	case 106: return 21; /* j */
	case 122: return 22; /* z */
	case 111: return 23; /* o */
	case 117: return 24; /* u */
	case 120: return 25; /* x */

	case 95:  return 26; /* _ */
	case 46:  return 26; /* . */

	default:
		printf("CAUTION: This is not a valid residue!\n");
		exit(-1);
		return 9999999;   /* ERROR */
	}
}