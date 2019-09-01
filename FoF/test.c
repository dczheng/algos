#include "stdio.h"
#include "time.h"
#include "fof.h"
#include "string.h"

int NTask, ThisTask;
int *map, W, H;

void main( int argc, char *argv[] ) {

    int Ncir, i, j, k, x, y, Ngroup;
    time_t time1, time2;
    struct tm *tb;
    double r;
    FILE *fd;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &ThisTask );
    MPI_Comm_size( MPI_COMM_WORLD, &NTask );

    Ncir = 200;
    H = W = 4000;

    map = malloc( sizeof(int) * W*H );
    if ( ThisTask == 0 ) {

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
    
        fd = fopen( "map.dat", "w" );
        for( i=0; i<H; i++ ) {
            for( j=0; j<W; j++ )
                fprintf( fd, "%i ", map[i*W+j] );
            fprintf( fd, "\n" );
        }
        fclose( fd );
    }

    MPI_Bcast( map, sizeof(int)*W*H, MPI_BYTE, 0, MPI_COMM_WORLD );

    time1 = time( NULL );
    fof_init( map, W, H );
    time2 = time( NULL );

    printf( "[%i] init time: %li sec.\n", ThisTask, (int)difftime( time2, time1 ) );


    MPI_Barrier( MPI_COMM_WORLD );

    time1 = time( NULL );
    fof();
    time2 = time( NULL );
    printf( "[%i] fof time: %li sec.\n", ThisTask, (int)difftime( time2, time1 ) );
    MPI_Barrier( MPI_COMM_WORLD );
    exit(0);

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

}
