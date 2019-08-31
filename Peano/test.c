/*
    dczheng
    created 2019-08-30
    Reference: https://www.cs.dal.ca/research/techreports/cs-2006-07
*/

#include "peano.h"
#include "stdio.h"

void print_bin( unsigned int n, int nn ) {
    int i, t, b;

    b = sizeof(n) * 8;

    if ( nn <= 0 )
        nn = b;

    for( i=b-nn; i<b; i++ ) {
        t = n & ( 1<<(b-1-i) ) ? 1 : 0;
        printf( "%i", t );
    }
    printf( "\n" );

}

void main() {

    unsigned int i, j, ii;
    for( i=0; i<32; i++ ) {
        j = BRGC(i);
        printf( "%4i: ", i );
        print_bin( j, 10 );
        ii = BRGC_inv( j );
        printf( "%4i: \n", ii );
        printf( "%4i: \n\n", DCGC(i) );

    }
}
