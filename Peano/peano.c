/*
    dczheng
    created 2019-08-30
    Reference: https://www.cs.dal.ca/research/techreports/cs-2006-07
*/

#include "stdio.h"

unsigned int BRGC( unsigned int n )  {
    return n ^ ( n>> 1 );
}

unsigned int BRGC_inv( unsigned int n ) {
    unsigned int i, j, r, nn;
    int b, t;
    nn = sizeof(n)*8;
    for(b=nn-1; b>=0; b-- ) {
        if ( (n>>b) & 1 ) 
            break;
    }
    b++;

    for( i=0, r=0; i<nn; i++ ) {
        for( j=i, t=0; j<b; j++ )
            t += ((n>>j) & 1);
        t %= 2;
        r |= t<<i;
    }

    return r;
}

unsigned int DCGC( unsigned int n )  {
    unsigned int a, i;
    a = BRGC( n ) ^ BRGC( n+1 );
    i=0;
    while(a) {
        i++;
        a >>= 1;
    }
    return i-1;
}

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
