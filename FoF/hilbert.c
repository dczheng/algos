/*
    dczheng
    created 2019-08-30
    Reference: ~/doc/hilbert.doc copyed from https://www.cs.dal.ca/research/techreports/cs-2006-07
    http://pdebuyl.be/blog/2015/hilbert-curve.html
*/

#include "hilbert.h"

myint BITS, MASK, N_U_INT;

myint gc( myint a ) {
    a &= MASK;
    return ( a ^ ( a>>1 ) );
}

myint gc_inv( myint a ) {
    myint b, i;
    short n=BITS-1;
    a &= MASK;
    while( n>=0 ) {
        if ( (a>>n)&1 )
            break;
        n--;
    }
    n++;
    b = a;
    i = 1;
    while( i<n ){
        b ^= ( a>>i );
        i++;
    }
    return b;
}

myint dcgc( myint a ) {
    myint b, t;
    a &= MASK;
    b = 1;
    t = a & 1;
    while( t ) {
        t = ( a>>b ) & 1;
        b++;
    }
    return b-1;
}

myint rrot( myint a, myint n ) {
    a &= MASK;
    n %= BITS;
    while( n>0 ){
        a = ( a>>1 ) + ( (a&1)<<(BITS-1) );
        n--;
    } 
    return a;
}

myint lrot( myint a, myint n ) {
    n %= BITS;
    while( n>0 ){
        a = ( a<<1 ) + ( (a&(1<<(BITS-1)))>>(BITS-1) );
        a &= MASK;
        n--;
    } 
    return a;
}

myint Ttrans( myint a, myint e, myint d ) {
    return rrot( a^e, d+1 );
}

myint Ttrans_inv( myint a, myint e, myint d ) {
    return lrot( a, d+1 ) ^ e;
}

myint ep( myint a ) {
    a &= MASK;
    if ( a == 0 )
        return 0;
    return gc( 2 * ( (a-1)/2 ) );
}


myint gd( myint a ) {
    return dcgc( a );
}

myint dd( myint a ) {
    a &= MASK;
    if ( a==0 )
        return 0;
    if ( a%2 == 0 )
        return gd( a-1 ) % BITS;

    return gd(a) % BITS;
}

void print_bin( myint a, int n ) {
    short b, t;

    if ( n == 0 )
        n= BITS;
    if ( n < 0 )
        n = sizeof(myint)*8;
    n--;

    while( n>=0 ) {
        t = a & ( 1<<(n) ) ? 1 : 0;
        printf( "%i", t );
        n--;
    }
    printf( "\n" );

}

void bits_init( myint bits ) {
    int i;
    BITS = bits;
    MASK = 0;
    i = 0;
    while(i<BITS) {
        MASK |= 1<<i;
        i++;
    }
}

void hilbert_index( myint *idx, myint *p, myint dim, myint bits ) {
    myint e, l, w, d, h;
    int i, j;

    bits_init( dim );
    h = e = 0;
    d = 1;

    for ( i=bits-1; i>=0; i-- ) {
        //printf( "i:%i ", i );
        l = 0;
        for( j=0; j<dim; j++ ) {
            l |= ((p[j]>>i)&1)<<j;
        }
        //printf( "l:%u ", l );
        l = Ttrans( l, e, d );
        //printf( "tl:%u ", l );
        w = gc_inv( l ) ;
        //printf( "w: %u ", w );
        //printf( "e(w): %u, d(w): %u ", ep(w), dd(w) );
        e = e ^ lrot( ep(w), d+1 );
        d = ( d + dd(w) + 1 ) % dim;
        //printf( "e: %u, d: %u ", e, d );
        h = ( h<<dim ) | w;
        //printf( "h: %u\n\n", h );
    }

    *idx = h;
}

void hilbert_index_inv( myint *idx, myint *p, myint dim, myint bits ) {

    int i, j;
    myint e, d, w, h, t, l;

    h = *idx;
    bits_init( dim );
    for( i=0; i<dim; i++ )
        p[i] = 0;
    e = 0;
    d = 1;

    for( i=bits-1; i>=0; i-- ) {
        for( j=0,w=0; j<dim; j++ ) {
            t = i*dim+j;
            w |= ((h&(1<<t))>>t)<<j;
        }
        l = gc( w );
        l = Ttrans_inv( l, e, d );
        for (j=0; j<dim; j++)
            p[j] |= ((l&(1<<j))>>j)<<i;

        e = e ^ lrot(ep(w), d+1);
        d = (d + dd(w) + 1) % dim;

    }
}

