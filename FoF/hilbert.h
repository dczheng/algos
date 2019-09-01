/*
    dczheng
    created 2019-08-30
    Reference: ~/doc/hilbert.doc copyed from https://www.cs.dal.ca/research/techreports/cs-2006-07
    http://pdebuyl.be/blog/2015/hilbert-curve.html
*/

#include "stdio.h"

typedef unsigned int myint;

extern myint BITS, MASK;

myint gc( myint a );
myint gc_inv( myint a );
myint dcgc( myint a );
myint rrot( myint a, myint n );
myint lrot( myint a, myint n );
myint GD( myint a );
myint dd( myint a );
myint ep( myint a );
myint Ttrans( myint a, myint e, myint d );
myint Ttrans_inv( myint a, myint e, myint d );
void  print_bin( myint n, int nn );

void hilbert_index( myint *idx, myint *p, myint dim, myint bits );
void hilbert_index_inv( myint *idx, myint *p, myint dim, myint bits );
