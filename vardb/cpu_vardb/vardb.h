/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   vardb.h
 * Author: eli
 *
 * Created on April 16, 2019, 9:15 AM
 */

#ifndef VARDB_H
#define VARDB_H

#ifdef __cplusplus
extern "C" {
#endif

typedef uint16_t u16;
typedef _Bool bool;
#define true (bool)1
#define false (bool)0
    
int vardb_add_rec(int col, int idx, u16 * data, int num_els, int num_rows);
int vardb_find_seq(u16 * ret_buf, int ret_buf_size, u16 * query, int num_els, int num_rows);
int vardb_find_match(u16 * ret_buf, int ret_buf_size, u16 * query, u16 * hds, u16 * var_delts, int num_els, int num_rows);
int vardb_find_dist(u16 * ret_buf, u16 * dist_buf, int ret_buf_size, u16 k, u16 * query, int num_els, int num_rows);



#ifdef __cplusplus
}
#endif

#endif /* VARDB_H */

