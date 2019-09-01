/*
    dczheng
    created 2019-08-30
    Reference: ~/doc/hilbert.doc copyed from https://www.cs.dal.ca/research/techreports/cs-2006-07
    http://pdebuyl.be/blog/2015/hilbert-curve.html
*/

#include "stdio.h"
#include "stdlib.h"
#include "hilbert.h"

void test_gc() {

    u_int i, j, ii;
    for( i=0; i<8; i++ ) {
        j = gc(i);
        printf( "%4u: ", i );
        print_bin( j, 10 );
        ii = gc_inv( j );
        printf( "%4u \n", ii );
        printf( "dcgc: %4u \n\n", dcgc(i) );

    }


}

void test_Rot() {
    u_int i;

    BITS = 8;
    i = 13;
    print_bin( i, 8 );
    print_bin( i, -1 );
    i = rrot( i, 3 );
    print_bin( i, 8 );
    print_bin( i, -1 );

    i = lrot( i, 4 );
    print_bin( i, 8 );
    print_bin( i, -1 );

    i = lrot( i, 0 );
    print_bin( i, 8 );
    print_bin( i, -1 );
    printf( "\n" );
}

void test_Ttrans() {

    u_int a, g, aa, e1, d1, e2, d2, e, d;

    e1 = 10;
    d1 = 3;
    e2 = 11;
    d2 = 4;

    print_bin( Ttrans( e1, e1, d ), -1 );
    printf( "\n" );

    a = 13;
    g = Ttrans( a, e1, d1 );
    g = Ttrans( g, e2, d2 );
    print_bin( g, -1 );

    e = e1 ^ ( lrot( e2, d1+1 ) );
    d = d1 + d2 + 1;
    g = Ttrans( a, e, d );
    print_bin( g, -1 );
    printf( "\n" );

    g = Ttrans( a, e1, d1 );
    aa = Ttrans_inv( g, e1, d1 );
    print_bin( a, -1 );
    print_bin( aa, -1 );



}

void test_hilbert() {
    u_int p[2], idx, pp[2], i, j, dim, bits;
    dim = 2;
    bits = 3;
    for( i=0; i<(1<<bits); i++ )
        for( j=0; j<(1<<bits); j++ ) {
            p[0] = i;
            p[1] = j;
            hilbert_index( &idx, p, dim, bits );
            hilbert_index_inv( &idx, pp, dim, bits );
            printf( "(%3u,%3u) -> %4u -> (%3u,%3u) \n",
            p[0], p[1], idx, pp[0], pp[1] );
        }
}

int hidx_compare2( const void *a, const void *b ) {
    return ( ((u_int*)a)[2] > ((u_int*)b)[2] ) ? 1 : -1;
}

void test_hilbert2() {
    u_int *p_idx, i, j, dim, bits, N, idx;
    FILE *fd;
    dim = 2;
    bits = 4;
    N = 1<<bits;
    printf( "N: %u\n", N );
    p_idx = malloc( sizeof(u_int) * N * N  * 3 );
    for( i=0; i<N; i++ ) 
        for( j=0; j<N; j++ ) {
            idx = (i*N + j)*3;
            p_idx[idx] = i;
            p_idx[idx+1] = j;
            hilbert_index( p_idx+idx+2, p_idx+idx, dim, bits );
        }

    qsort( p_idx, N*N, sizeof(u_int)*3, hidx_compare2 );

    fd = fopen( "hidx2.dat", "w" );
    for( i=0; i<N; i++ ) 
        for( j=0; j<N; j++ ) {
            idx = (i*N + j)*3;
            fprintf( fd, "%u %u %u\n", p_idx[idx+2], p_idx[idx], p_idx[idx+1] );
        }
    fclose(fd);
    free( p_idx );

}

int hidx_compare3( const void *a, const void *b ) {
    return ( ((u_int*)a)[3] > ((u_int*)b)[3] ) ? 1 : -1;
}

void test_hilbert3() {
    u_int *p_idx, i, j, k, dim, bits, N, idx;
    FILE *fd;
    dim = 3;
    bits = 3;
    N = 1<<bits;
    printf( "N: %u\n", N );
    p_idx = malloc( sizeof(u_int) * N * N * N  * 4 );
    for( i=0; i<N; i++ ) 
        for( j=0; j<N; j++ )
            for( k=0; k<N; k++ ) {
            idx = ((i*N + j) * N + k) * 4;
            p_idx[idx] = i;
            p_idx[idx+1] = j;
            p_idx[idx+2] = k;
            //p_idx[idx+3] = i+j+k;
            hilbert_index( p_idx+idx+3, p_idx+idx, dim, bits );
        }

    qsort( p_idx, N*N*N, sizeof(u_int)*4, hidx_compare3 );

    fd = fopen( "hidx3.dat", "w" );
    for( i=0; i<N; i++ ) 
        for( j=0; j<N; j++ )
            for( k=0; k<N; k++ ) {
            idx = ( (i*N + j) * N + k)*4;
            fprintf( fd, "%u %u %u %u\n", p_idx[idx+3], p_idx[idx], p_idx[idx+1], p_idx[idx+2] );
        }
    fclose(fd);
    free( p_idx );

}

void main() {

    printf( "test:\n" );
    //test_gc();
    //test_Rot();
    //test_Ttrans();
    //test_hilbert();
    test_hilbert2();
    test_hilbert3();

}
