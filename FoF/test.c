#include "stdio.h"
#include "time.h"
#include "fof.h"
#include "string.h"

int NTask, ThisTask;
int *map, W, H;

int main( int argc, char *argv[] ) {

    int Ncir, i, j, k, x, y, Ngroup;
    clock_t t0, t1;
    double dt, dt_min, dt_max, dt_tot, r;
    FILE *fd;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &ThisTask );
    MPI_Comm_size( MPI_COMM_WORLD, &NTask );

    Ncir = 3000;
    H = W = 10000;

    map = malloc( sizeof(int) * W*H );

    srand( 0 );
    memset( map, 0, sizeof(int) * W*H );

    for( k=0; k<Ncir; k++ ) {

        x = (int)( (double)rand() / RAND_MAX * ( W-10 ) + 5 );
        y = (int)( (double)rand() / RAND_MAX * ( H-10 ) + 5 );
        r = (double)rand() / RAND_MAX * ( W/200 ) + 5;
    
        for( i=y-r; i<y+r; i++ ) {
            for(j=x-r; j<x+r; j++ ) {
                if ( i<0 || i>H-1 
                  || j<0 || j>W-1 )
                    continue;
                    if ( (i-y)*(i-y ) + (j-x)*(j-x) < r*r )
                        map[ i*W + j ] = 1;
            }
        }
    }
    
    if ( ThisTask == -1 ) {
        fd = fopen( "map.dat", "w" );
        for( i=0; i<H; i++ ) {
            for( j=0; j<W; j++ )
                fprintf( fd, "%i ", map[i*W+j] );
            fprintf( fd, "\n" );
        }
        fclose( fd );
    }

    MPI_Barrier(MPI_COMM_WORLD);

    t0 = t1 = clock();
    fof_init( map, W, H );
    dt = (double)(clock()-t1) / CLOCKS_PER_SEC;
    MPI_Reduce( &dt, &dt_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Reduce( &dt, &dt_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce( &dt, &dt_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    if ( ThisTask == 0 )
        printf( "[fof_init] dt_min: %g, dt_max: %g, dt_mean: %g\n",
                dt_min, dt_max, dt_tot/NTask );

    //MPI_Barrier( MPI_COMM_WORLD );

    t1 = clock();
    fof();
    dt = (double)(clock()-t1) / CLOCKS_PER_SEC;
    MPI_Reduce( &dt, &dt_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
    MPI_Reduce( &dt, &dt_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    MPI_Reduce( &dt, &dt_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    if ( ThisTask == 0 )
        printf( "[fof] dt_min: %g, dt_max: %g, dt_mean: %g\n",
                dt_min, dt_max, dt_tot/NTask );

    if ( ThisTask == 0 ) {
        for( i=0, Ngroup=0; i<W*H; i++ )
            if ( Len[i] > 1 )
                Ngroup++ ;

        printf( "Maxgroup: %i, Mingroup: %i, Ngroup: %i \n",
            Len[0], Len[Ngroup-1], Ngroup );

        fd = fopen( "groups.dat", "w" );
        for( k=0; k<Ngroup; k++ ) {
            fprintf( fd, "%i ", k+1 );
            i = Head[k];
            while(i>0) {
                fprintf( fd, "%i %i ", i%W, i/W );
                i = Next[i];
            }
            fprintf( fd, "\n" );
        }

        fclose( fd );
    }

    fof_free();
    free( map );

    MPI_Finalize();

    return 0;

}
