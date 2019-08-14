#include "stdio.h"
#include "fof.h"

void main() {

    int *map, N, i;
    N = 300;
    map = malloc( sizeof(int) * N * N );

    for( i=0; i<N*N; i++ )
        map[i] = 0;
    fof_init( map, N, N );
    fof();
    fof_free();

    free( map );

}
