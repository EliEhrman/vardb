/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <stdio.h>
#include <stdint.h>
#include <memory.h>
#include "vardb.h"

#define c_num_rows 13
#define c_num_cols (1 << 15)
#define c_num_markers 5
u16 g_mem[c_num_rows][c_num_cols]; // in c this naming is normally reversed
bool g_markers[c_num_markers][c_num_cols];
int c_vr_idxs = 4;
int c_wstarts_imrk = 0;
int g_end_col = 0;

static u16 popcount(u16 a)
{
	u16 count = 0;
	while (a) {
		count += (u16)(a & (u16)0x1);
		a >>= 1;
	}
	return count;
}

int vardb_init() {
	
}
int vardb_add_rec(int col, int idx, u16 * data, int num_els, int num_rows) {
	printf("vardb_add_rec for cpu vardb.\n");
	if (col != g_end_col) {
		printf("Error! Record does not follow previous. Gaps are not implemented yet.\n");
		return 1;
	}
	for (int irow = 0; irow<num_rows; irow++) {
		memcpy(&(g_mem[irow][col]), &(data[(irow*num_els)]), sizeof(u16) * num_els);
	}
	for (int iel = 0; iel < num_els; iel++) {
		g_mem[c_vr_idxs][col+iel] = idx;
	}
	g_markers[c_wstarts_imrk][col] = true;
	if (col >= g_end_col) {
		g_end_col = col + num_els;
	}
	return 0;
}

int vardb_find_seq(u16 * ret_buf, int ret_buf_size, u16 * query, int num_els, int num_rows) {
	int num_found = 0;
	for (int icol = 0; icol < c_num_cols; icol++) {
		if (icol == g_end_col) {
			ret_buf[num_found] = 0xffff;
			return 0;
		}
		if (!g_markers[c_wstarts_imrk][icol]) continue;
		int start_col = icol;
		for (; icol < c_num_cols; icol++) {
			if (icol == g_end_col) {
				icol--;
				break;
			}
			int iq = icol - start_col;
			if (iq == num_els) {
				if (g_markers[c_wstarts_imrk][icol]) {
					ret_buf[num_found++] = g_mem[c_vr_idxs][icol] - 1;
					icol--;
				}				
				break;
			}
			if (iq != 0 && g_markers[c_wstarts_imrk][icol]) {
				icol--;
				break;
			}
			bool bmatched = true;
			for (int irow=0; irow < num_rows; irow++) {
				if (g_mem[irow][icol] != query[(irow*num_els) + iq]) {
					bmatched = false;
					break;
				}
			}
			if (!bmatched ) {
				break;
			}
		} 
	}
	printf("Error! Should not get here.");
	return 1;
}

int vardb_find_dist(u16 * ret_buf, u16 * dist_buf, int ret_buf_size, u16 k, u16 * query, int num_els, int num_rows) {
	u16 num_found = 0;
	u16 max_dist = 0;
	u16 iworst = 0;
	for (int icol = 0; icol < c_num_cols; icol++) {
		if (icol == g_end_col) {
			printf("vardb_find_dist completed with %d found.\n", num_found);
			ret_buf[num_found] = 0xffff;
			return 0;
		}
		if (!g_markers[c_wstarts_imrk][icol]) continue;
		int start_col = icol;
		u16 diff = 0;
		for (; icol < c_num_cols; icol++) {
			if (icol == g_end_col) {
				icol--;
				break;
			}
			int iq = icol - start_col;
			if (iq == num_els) {
				if (g_markers[c_wstarts_imrk][icol]) {
					if (num_found < k) {
						if (diff >= max_dist) {
							iworst = num_found;
							max_dist = diff;
						}
						dist_buf[num_found] = diff;
						ret_buf[num_found++] = g_mem[c_vr_idxs][icol] - 1;
						printf("vardb_find_dist placed idx %d at loc %d.\n", g_mem[c_vr_idxs][icol] - 1, num_found-1);
					}
					else {
						if (diff >= max_dist) {
							icol--;
							break;
						}
						dist_buf[iworst] = diff;
						ret_buf[iworst] = g_mem[c_vr_idxs][icol] - 1;
						printf("vardb_find_dist placed idx %d at loc %d.\n", g_mem[c_vr_idxs][icol] - 1, iworst);
						max_dist = 0;
						for (int ifound = 0; ifound < k; ifound++) {
							if (dist_buf[ifound] >= max_dist) {
								iworst = ifound;
								max_dist = dist_buf[ifound];
							}
						}						
					}
					icol--;
				}				
				break;
			}
			if (iq != 0 && g_markers[c_wstarts_imrk][icol]) {
				icol--;
				break;
			}
			for (int irow=0; irow < num_rows; irow++) {
				u16 qval = query[(irow*num_els) + iq];
				diff += popcount(g_mem[irow][icol] ^ qval);
			}
		} 
	}
	printf("Error! Should not get here.");
	return 1;
}


int vardb_find_match(u16 * ret_buf, int ret_buf_size, u16 * query, u16 * hds, u16 * var_delts, int num_els, int num_rows) {
	int num_found = 0;
	for (int icol = 0; icol < c_num_cols; icol++) {
		if (icol == g_end_col) {
			ret_buf[num_found] = 0xffff;
			return 0;
		}
		if (!g_markers[c_wstarts_imrk][icol]) continue;
		int start_col = icol;
		for (; icol < c_num_cols; icol++) {
			if (icol == g_end_col) {
				icol--;
				break;
			}
			int iq = icol - start_col;
			if (iq == num_els) {
				if (g_markers[c_wstarts_imrk][icol]) {
					ret_buf[num_found++] = g_mem[c_vr_idxs][icol] - 1;
					icol--;
				}				
				break;
			}
			if (iq != 0 && g_markers[c_wstarts_imrk][icol]) {
				icol--;
				break;
			}
			bool bmatched = true;
			u16 diff = 0;
			for (int irow=0; irow < num_rows; irow++) {
				u16 qval;
				if (var_delts[iq] > 0) {
					qval = g_mem[irow][icol - var_delts[iq]];
				}
				else {
					qval = query[(irow*num_els) + iq];
				}
				diff += popcount(g_mem[irow][icol] ^ qval);
				if (diff > hds[iq]) {
					bmatched = false;
					break;
				}
			}
			if (!bmatched ) {
				break;
			}
		} 
	}
	printf("Error! Should not get here.");
	return 1;
}
