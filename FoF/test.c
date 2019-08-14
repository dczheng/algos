#include "stdio.h"
#include "fof.h"

void main() {

    int *map, N;
    fof_init( N );
    fof( map );
    fof_free();

}
