/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   test_vardb.c
 * Author: eli
 *
 * Created on April 9, 2019, 12:09 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <wordexp.h>
#include <stdint.h>
#include <math.h>
#include "../cpu_vardb/vardb.h"

#define BLEN 16
#define NUMLEN 3
#define NUM_ROWS 4
char * fnt = "~/tmp/vardb_samples.txt";
static char line[4096];


/*
 * 
 */
int main(int argc, char** argv) {
	wordexp_t exp_result;
	wordexp(fnt, &exp_result, 0);
	printf("Opening %s.\n", exp_result.we_wordv[0]);
	FILE * fh;
	fh = fopen(exp_result.we_wordv[0], "rt");
	if (fh == NULL) {
		printf("Fail to open file.\n");
		exit(1);
	}
	wordfree(&exp_result);

	int iel = 0;
	int iline = 0;
    while (fgets(line, sizeof(line), fh) != NULL) {
        /* note that fgets don't strip the terminating \n, checking its
           presence would allow to handle lines longer that sizeof(line) */
        printf("%s", line); 
		bool bquery = true;
		if (line[0] == '0') {
			bquery = false;
		}
		char snum[NUMLEN+1];
		memcpy(snum, line+1, NUMLEN);
		snum[NUMLEN] = '\0';
		int numw = atoi(snum);
		printf("Got %d words.\n", numw);
		int pos = NUMLEN + 1;
		u16 data[NUM_ROWS * numw];
		u16 var_deltas[numw]; u16 hds[numw];
		for (int iw=0; iw<numw; iw++) {
			char bbuf[BLEN+1];
			memcpy(bbuf, &line[pos], BLEN);
			bbuf[BLEN] = '\0';
			pos += BLEN;
			for (int iu16 = 0; iu16 < NUM_ROWS; iu16++) {
				char * p = &(bbuf[iu16*4]);
				u16 v = 0;
				for (int inib = 0; inib < 4; inib++, p++) {
					char hbuf[2];
					hbuf[0] = *p; hbuf[1] = '\0';
					v += (u16)strtol(hbuf, NULL, 16) * (u16)pow((double)16, (3-inib));
				}
				data[(iu16 * numw) + iw] = v;
			}
			char intbuf[3];
			memcpy(intbuf, &(line[pos]), 2); intbuf[2] = '\0';
			hds[iw] = atoi(intbuf);
			pos += 2;
			memcpy(intbuf, &(line[pos]), 2); intbuf[2] = '\0';
			var_deltas[iw] = atoi(intbuf);
			pos += 2;
			
		}
		printf("About to call library.\n");
		if (!bquery) {
			vardb_add_rec(iel, iline, data, numw, NUM_ROWS);
			iel += numw;
		}
		else {
			u16 ret_buf[iline+1];
			u16 dist_buf[iline+1];
			u16 k = 3;
//			vardb_find_seq(ret_buf,iline+1,data, numw, NUM_ROWS);
//			vardb_find_match(ret_buf,iline+1,data, hds, var_deltas, numw, NUM_ROWS);
			vardb_find_dist(ret_buf, dist_buf, iline+1, k, data, numw, NUM_ROWS);
		}
		iline++;
    }

	fclose(fh);
	return (EXIT_SUCCESS);
}

