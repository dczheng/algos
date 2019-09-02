/*
    dczheng
    created 2019-08-14
*/

#include "fof.h"

// a implementation of 2D fof.

int *Next, *Head, *Tail, *Len, Width, Height, FoFN, *FoFMap, *Flag,
#ifdef HILBERT_INDEX
     nidx_tot, bits, nidx_per_task,
#endif
     Xmin, Ymin, Xmax, Ymax, DX, DY;

//#define FOF_SINGLE_DEBUG
#define FOF_SINGLE_DEBUG2

#ifdef FOF_SINGLE_DEBUG2 
    int PN;
#endif

void fof_single_finder( int p0, int mode ) {

    int  i, j, x0, y0;
    int p1, l, s, ss, p_next;

    if ( FoFMap[p0] == 0 )
        return;

    Flag[p0] = 1;

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

            if ( mode )
                if ( y0+i<Ymin ||
                     y0+i>Ymax-1 ||
                     x0+j<Xmin ||
                     x0+j>Xmax-1 )
                    continue;

#ifdef FOF_SINGLE_DEBUG
            printf( "1\n" );
            fflush( stdout );
#endif
            if ( Head[p1] == Head[p0] )
                continue;

            Flag[p1] = 1;

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

#ifdef FOF_SINGLE_DEBUG2
        PN++;
#endif
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
        fof_single_finder( p_next, mode );

}

int hidx_compare( const void *a, const void *b ) {
    return ( *((int*)a) > *((int*)b) ) ? 1 : -1;
}

void fof_init( int *map, int W, int H ) {

    int i;

//#define FOF_INIT_DEBUG

#ifdef FOF_INIT_DEBUG
    FILE *fd;
    char buf[100];
    int p[2];
#endif

    Width = W;
    Height = H;
    FoFN = W * H;
    FoFMap = map;

    Next = malloc( sizeof(int) * FoFN );
    Head = malloc( sizeof(int) * FoFN );
    Len = malloc( sizeof(int) * FoFN );
    Tail = malloc( sizeof(int) * FoFN );
    Flag = malloc( sizeof(int) * FoFN );

    for( i=0; i<FoFN; i++ ) {
        Next[i] = -1;
        Head[i] = i;
        Tail[i] = i;
        Len[i]  = 1;
        Flag[i] = 0;
    }

    if ( NTask == 1 ) 
        return;

#ifdef HILBERT_INDEX
    for( bits=1; (1<<bits)<sqrt(NTask); bits++ );
    nidx_tot = 1<<(bits*2);
    nidx_per_task = nidx_tot / NTask;
    if ( nidx_tot % NTask != 0 )
        nidx_per_task ++;
    DX = Width  / (1<<bits);
    DY = Height / (1<<bits);

#ifdef FOF_INIT_DEBUG
    printf( "[%i] bits: %i [%i], nidx_tot: %i, nidx_per_task: %i, "
            "DXY: (%i, %i)\n",
    ThisTask,
    bits, 1<<bits, nidx_tot, nidx_per_task, DX, DY );
    sprintf( buf, "p_%03i.dat", ThisTask );
    fd = fopen( buf, "w" );
    for( i=0; i<nidx_tot; i++ ) {
        if ( i / nidx_per_task != ThisTask )
            continue;
        hilbert_index_inv( &i, p, 2, bits );
        fprintf( fd, "%i %i\n", p[0], p[1] );
    }
#endif
#ifdef FOF_INIT_DEBUG
    fclose( fd );
#endif

#else
    DX = Width / NTask;
    DY = Height / NTask;
#endif

}

void fof_reset() {

    int i;
    for( i=0; i<FoFN; i++ ) {
        Next[i] = -1;
        Head[i] = i;
        Tail[i] = i;
        Len[i] = 1;
        Flag[i] = 0;
    }
}

void fof_free() {

    free( Next );
    free( Head );
    free( Tail );
    free( Len );
    free( Flag );

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

typedef struct FoFData{
    int Next, Len, Head, Tail, index;
}FoFData;

void fof() {
    int p, k, i, j, n_local, *n_all, *n_dis, n_tot, n_pad;
    FoFData *FDat, *FDat_tot;
    n_pad = 0;
#ifdef HILBERT_INDEX
    int pp[2];
#endif

#define FOF_DEBUG
#ifdef FOF_SINGLE_DEBUG2
    PN = 0;
#endif
    if ( NTask == 1 ) {
        for( p=0; p<FoFN; p++ )
            fof_single_finder( p, 0 );
#ifdef FOF_SINGLE_DEBUG2
    printf( "[%i] PN: %i\n", ThisTask, PN );
#endif
        return;
    }

#ifdef HILBERT_INDEX
    for( k=0, n_local=0; k<nidx_tot; k++ ) {
        if ( k / nidx_per_task != ThisTask )
            continue;
        hilbert_index_inv( &k, pp, 2, bits );
        Xmin = pp[0] * DX;
        Xmax = Xmin + DX-n_pad; 
        Ymin = pp[1] * DY;
        Ymax = Ymin + DY-n_pad; 
#else
    for( k=0, n_local=0; k<NTask; k++ ){
        if ( k != ThisTask )
            continue;
        Xmin = DX * k;
        Xmax = Xmin + DX-n_pad;
        Ymin = 0;
        Ymax = Height;
#endif
        n_local += ( Xmax-Xmin ) * ( Ymax-Ymin );
#ifdef FOF_DEBUG
        printf( "[%i] X:(%i, %i), Y:(%i, %i)\n", ThisTask,
                Xmin, Xmax, Ymin, Ymax );
#endif
        for( i=Ymin; i<Ymax; i++ )
            for ( j=Xmin; j<Xmax; j++ ){
                p = i * Width + j; 
                fof_single_finder( p, 1 );
            }
    }

#ifdef FOF_SINGLE_DEBUG2
    printf( "[%i] PN: %i\n", ThisTask, PN );
#endif
#ifdef FOF_DEBUG
    printf( "[%i] n_local: %i\n", ThisTask, n_local );
#endif

    FDat = malloc( sizeof(FoFData) * n_local );
    if ( NULL == FDat ) {
        printf( "malloc error [%li]\n", sizeof(FoFData)*n_local );
        exit(0);
    }

#ifdef HILBERT_INDEX
    for( k=0, n_local=0; k<nidx_tot; k++ ) {
        if ( k / nidx_per_task != ThisTask )
            continue;
        hilbert_index_inv( &k, pp, 2, bits );
        Xmin = pp[0] * DX;
        Xmax = Xmin + DX-n_pad; 
        Ymin = pp[1] * DY;
        Ymax = Ymin + DY-n_pad; 
#else
    for( k=0, n_local=0; k<NTask; k++ ){
        if ( k != ThisTask )
            continue;
        Xmin = DX * k;
        Xmax = Xmin + DX-n_pad;
        Ymin = 0;
        Ymax = Height;
#endif
        for( i=Ymin; i<Ymax; i++ )
            for ( j=Xmin; j<Xmax; j++ ){
                p = i * Width + j; 
                if ( p < 0       ||
                     p > FoFN-1 )
                    continue;

                if ( Flag[p] ){
#ifndef HILBERT_INDEX
                    if ( (p%Width) / DX != ThisTask ) {
                        printf( "[%i] can't be, [%i]\n",
                        ThisTask, (p%Width)/DX );
                        exit(1);
                    }
#endif
                    FDat[n_local].Next = Next[p];
                    FDat[n_local].Head = Head[p];
                    FDat[n_local].Tail = Tail[p];
                    FDat[n_local].Len  = Len[p];
                    FDat[n_local].index  = p;
                    n_local ++;
                }
            }

    }

    n_all = malloc( sizeof(int) * NTask );
    n_dis = malloc( sizeof(int) * NTask );
    MPI_Gather( &n_local, 1, MPI_INT, n_all, 1, MPI_INT, 0, MPI_COMM_WORLD );
    if ( ThisTask == 0 ){
        for( i=0, n_tot=0; i<NTask; i++ ) {
            n_tot += n_all[i];
            if (i==0 )
                n_dis[i] = 0;
            else
                n_dis[i] = n_dis[i-1] + n_all[i-1];
        }

#ifdef FOF_DEBUG
        printf( "n_tot: %i\n", n_tot );
#endif
        FDat_tot = malloc( sizeof(FoFData) * n_tot );
    }

    MPI_Gatherv( FDat, n_local, MPI_INT, 
                 FDat_tot, n_all, n_dis, MPI_INT, 0, MPI_COMM_WORLD );
    free( n_all );
    free( n_dis );
    free( FDat );

    if ( ThisTask )
        return;

    for( i=0; i<n_tot; i++ ) {
        p = FDat_tot[i].index;
#ifdef HILBERT_INDEX
#else
        if ( (p%Width)/DX == 0 )
            continue;
#endif
        if ( Flag[p] == 1 ) {
            printf( "can't be\n" );
            exit(1);
        }
        Next[p] = FDat_tot[i].Next;
        Head[p] = FDat_tot[i].Head;
        Tail[p] = FDat_tot[i].Tail;
        Len[p] = FDat_tot[i].Len;

    }

    //fof_sort();
}
