#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

typedef double (*MYFUNC)( double, double );

struct CN_Data{
    double *x, t_step, *y0, *y, xmin, xmax, t_start;
    int xN, type;
    MYFUNC f_a;
    MYFUNC f_b;
    MYFUNC f_c;
    MYFUNC f_d;
    void (*boundary)( double *left, double *right, double *x, int N, double t );
    void (*condition)( double *y, double *x, int N );

};

void CN_init( struct CN_Data p );
void CN_forward( double t_forward );
void CN_test();
void CN_free();
void TridMat( double *a, double *b, double *c, double *d, double *x, int N, int flag );
