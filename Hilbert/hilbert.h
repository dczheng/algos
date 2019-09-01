/*
    dczheng
    created 2019-08-30
    Reference: https://www.cs.dal.ca/research/techreports/cs-2006-07
*/

#include "stdio.h"

typedef unsigned int u_int;

extern u_int BITS, MASK;

u_int gc( u_int a );
u_int gc_inv( u_int a );
u_int dcgc( u_int a );
u_int rrot( u_int a, u_int n );
u_int lrot( u_int a, u_int n );
u_int GD( u_int a );
u_int dd( u_int a );
u_int ep( u_int a );
u_int Ttrans( u_int a, u_int e, u_int d );
u_int Ttrans_inv( u_int a, u_int e, u_int d );
void  print_bin( u_int n, int nn );

void hilbert_index( u_int *idx, u_int *p, u_int dim, u_int bits );
void hilbert_index_inv( u_int *idx, u_int *p, u_int dim, u_int bits );
