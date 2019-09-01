#include "stdio.h"
#include "stdlib.h"
#include "mpi.h"
#include "hilbert.h"
#include "limits.h"

#define PARALLELFOF

#ifdef PARALLELFOF
    extern int NTask, ThisTask;
#endif

extern int *Next, *Head, *Tail, *Len;

void fof_init( int *map, int W, int H );
void fof_free();
void fof_reset();
void fof();
