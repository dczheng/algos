/*
    dczheng
    created 2019-08-14
*/

#include "fof.h"

// a implementation of 2D fof.

int *Next, *Head, *Tail, *Len, Width, Height, FoFN, *FoFMap,
        Xmin, Ymin, Xmax, Ymax;

//#define FOF_SINGLE_DEBUG
void fof_single_finder( int p0 ) {

    int  i, j, x0, y0;
    int p1, l, s, ss, p_next;

    if ( FoFMap[p0] == 0 )
        return;

    x0 = p0 % Width;
    y0 = p0 / Width;
    p_next = -1;

#ifdef FOF_SINGLE_DEBUG
    printf( "start\n" );
    fflush( stdout );
#endif
    for( i=-1; i<2; i++ )
        for( j=-1; j<2; j++ ) {
            p1 = (y0+i) * Width + (x0+j);
#ifdef FOF_SINGLE_DEBUG
            printf( "p0: %li, p1:  %li, NPixs: %li (%i*%i)\n",
            p0, p1, FoFN, Width, Height );
            fflush( stdout );
#endif
            if ( p1 == p0     ||
                 p1 < 0       ||
                 p1 > FoFN-1 ||
                 FoFMap[p1]  == 0 )
                continue;

#ifdef FOF_SINGLE_DEBUG
            printf( "1\n" );
            fflush( stdout );
#endif
            if ( Head[p1] == Head[p0] )
                continue;

            if ( Len[Head[p1]] > Len[Head[p0]] ) {
                l = p1;
                s = p0;
             }
             else {
                l = p0;
                s = p1;
             }

#ifdef FOF_SINGLE_DEBUG
            printf( "2\n" );
            fflush( stdout );
#endif
            Next[Tail[Head[l]]] = Head[s];
            Tail[Head[l]] = Tail[Head[s]];
            Len[Head[l]] += Len[Head[s]];
            Len[Head[s]] = 1;

            ss = Head[s];
            do{
                /*
                printf( "p0: %li, p1: %li l: %li, s: %li, ss: %li\n"
                    "Head[l]: %li, Head[s]: %li\n",
                p0, p1, l, s, ss, Head[l], Head[s] );
                */
                Head[ss] = Head[l];
            } while( (ss=Next[ss]) >= 0 );

#ifdef FOF_SINGLE_DEBUG
            printf( "3\n" );
            fflush( stdout );
#endif
            p_next = p1;

    }

#ifdef FOF_SINGLE_DEBUG
    printf( "end\n\n" );
    fflush( stdout );
#endif

    if ( p_next >= 0 )
        fof_single_finder( p_next );

}

int hidx_compare( const void *a, const void *b ) {
    return ( *((int*)a) > *((int*)b) ) ? 1 : -1;
}

void fof_init( int *map, int W, int H ) {

    int i, *nidx_all, nidx_tot, p[2], *hidx_local, *hidx, *idx_dis,
        bits, nidx, n;

#define FOF_INIT_DEBUG

#ifdef FOF_INIT_DEBUG
    FILE *fd;
    char buf[100];
#endif

    Width = W;
    Height = H;
    FoFN = W * H;
    FoFMap = map;

    Next = malloc( sizeof(int) * FoFN );
    Head = malloc( sizeof(int) * FoFN );
    Len = malloc( sizeof(int) * FoFN );
    Tail = malloc( sizeof(int) * FoFN );
    hidx_local = malloc( sizeof(int) * FoFN );

    for( i=0; i<FoFN; i++ ) {
        Next[i] = -1;
        Head[i] = i;
        Tail[i] = i;
        Len[i] = 1;
        hidx_local[i] = -1;
    }

    for( bits=1; (1<<bits)<W; bits++ );
    for( ; (1<<bits)<H; bits++ );
    //printf( "bits: %i [%i]\n", bits, 1<<bits );

    for ( i=0, nidx=0; i<FoFN; i++ ) {
        if ( i % NTask != ThisTask ||
            FoFMap[i] == 0
            )
            continue;
        p[0] = i % W; 
        p[1] = i / W; 
        hilbert_index( hidx_local+nidx, p, 2, bits );
        nidx++;
    }
#ifdef FOF_INIT_DEBUG
    printf( "[%i] nidx: %i\n", nidx );
#endif

    nidx_all = malloc( sizeof(int)*NTask );
    idx_dis = malloc( sizeof(int)*NTask );

    MPI_Allgather( &nidx, 1, MPI_INT, 
                   nidx_all, 1, MPI_INT,
                   MPI_COMM_WORLD );
    nidx_tot = 0;
    for( i=0; i<NTask; i++ ){
        nidx_tot += nidx_all[i];
        if ( i==0 )
            idx_dis[i] = 0;
        else
            idx_dis[i] = idx_dis[i-1] + nidx_all[i-1];
    }

#ifdef FOF_INIT_DEBUG
    printf( "[%i] ndix_tot: %i\n", ThisTask, nidx_tot );
    if ( ThisTask==0 ) {
        printf( "nidx_all: ", ThisTask );
        for( i=0; i<NTask; i++ )
            printf( "%i ", nidx_all[i] );
        printf( "\nidx_dis: " );
        for( i=0; i<NTask; i++ )
            printf( "%i ", idx_dis[i] );
        printf( "\n" );
    }
#endif


    hidx = malloc( sizeof(int) * FoFN );
    MPI_Gatherv( hidx_local, nidx, MPI_INT, 
                 hidx, nidx_all, idx_dis, MPI_INT,
                 0, MPI_COMM_WORLD );
    free( hidx_local );
    free( nidx_all );

    if ( ThisTask == 0 ) {
        qsort( hidx, nidx_tot, sizeof(int), hidx_compare );

#ifdef FOF_INIT_DEBUG
    for( i=0; i<10; i++ )
        printf( "hidx[%i]: %i\n", i, hidx[i] );
#endif
    }
    MPI_Bcast( hidx, nidx_tot, MPI_INT, 0, MPI_COMM_WORLD );

    Ymin = Xmin = INT_MAX;
    Ymax = Xmax = -INT_MAX;
    n = nidx_tot / NTask;
#ifdef FOF_INIT_DEBUG
    sprintf( buf, "p_%03i.dat", ThisTask );
    fd = fopen( buf, "w" );
#endif
    for( i=ThisTask*n; i<(ThisTask+1)*n; i++ ) {
        if ( i > nidx_tot-1 )
            break;
        hilbert_index_inv( hidx+i, p, 2, bits );
        Xmin = ( p[0]<Xmin ) ? p[0] : Xmin;
        Xmax = ( p[0]>Xmax ) ? p[0] : Xmax;
        Ymin = ( p[1]<Ymin ) ? p[1] : Ymin;
        Ymax = ( p[1]>Ymax ) ? p[1] : Ymax;
#ifdef FOF_INIT_DEBUG
        fprintf( fd, "%i %i\n", p[0], p[1] );
#endif
    }
#ifdef FOF_INIT_DEBUG
    fclose( fd );
    printf( "[%i] X: (%i, %i), Y: (%i, %i)\n",
            ThisTask,
            Xmin, Xmax, Ymin, Ymax );
#endif
    free( hidx );

}

void fof_reset() {

    int i;
    for( i=0; i<FoFN; i++ ) {
        Next[i] = -1;
        Head[i] = i;
        Tail[i] = i;
        Len[i] = 1;
    }
}

void fof_free() {

    free( Next );
    free( Head );
    free( Tail );
    free( Len );

}

int fof_compare_len( const void *a, const void *b ) {
    return ( (*(int*)a) < *((int*)b) ) ? 1 : -1;
}

void fof_sort () {

    int *tmp, p;
    tmp = malloc( sizeof(int) * FoFN * 3 );
    for ( p=0; p<FoFN; p++ ) {
        tmp[ 3*p ] = Len[p];
        tmp[ 3*p+1 ] = Head[p];
        tmp[ 3*p+2 ] = Tail[p];
    }
    
    qsort( tmp, FoFN, sizeof(int)*3, fof_compare_len );

    for ( p=0; p<FoFN; p++ ) {
        Len[p] = tmp[ 3*p ];
        Head[p] = tmp[ 3*p+1 ];
        Tail[p] = tmp[ 3*p+2 ];
    }
    free( tmp );
}

void fof() {
    int p, i, j;
    for( i=Ymin; i<Ymax; i++ )
        for ( j=Xmin; j<Xmax; j++ ){
            p = i * Width + j; 
            fof_single_finder(  p );
        }
    fof_sort();
}
