#include "stdio.h"
#include "Crank-Nicolson.h"
#include "math.h"

#define CUBE(X) ((X)*(X)*(X))
#define SQR(X) ((X)*(X))

void TridMat_test(){

    int i;

    double a[4] = { 0, 3, 1, 3 };
    double b[4] = { 10, 10, 7, 4 };
    double c[4] = { 2, 4, 5, 0 };
    double d[4] = { 3, 4, 5, 6 };
    double x[4];

    TridMat( a, b, c, d, x, 4, 1 );
    for ( i=0; i<4; i++ )
        printf( "%g ", x[i] );
    printf( "\n" );

    TridMat( a, b, c, d, x, 4, 0 );
    for ( i=0; i<4; i++ )
        printf( "%g ", x[i] );
    printf( "\n" );


}

double f_zeros( double x, double t ){
    return 0;
}

double f_ones( double x, double t ){
    return 1;
}

void boundary_dirichlet( double *left, double *right, double *x, int N, double t ){

    left[0] = 1;
    left[1] = 0;
    left[2] = 0;
    right[0] = 0;
    right[1] = 1;
    right[2] = 0;
}

void boundary_newmann( double *left, double *right, double *x, int N, double t ){

    left[0] = 1;
    left[1] = -1;
    left[2] = 0;
    right[0] = 1;
    right[1] = -1;
    right[2] = 0;
}


MYFUNC f_a, f_b, f_c, f_d, f_0, f_th;
double xmax, xmin, dt, dtt, t0, t1;
int xN, type;
char FileName[100];
void (*boundary)( double *left, double *right, double *x, int N, double t );

#define myprint( fd, s, a, N ) {\
    int i;\
    fprintf( fd, s );\
    for( i=0; i<N; i++ )\
        fprintf( fd, "%g ", a[i] );\
    fprintf( fd, "\n" );\
}

void Crank_Nicolson_test() {

    struct CN_Data p;
    double dx, *y0, *y, t, *x;
    int i;
    char buf[100];
    FILE *fd;

    if ( f_a == NULL )
        f_a = f_zeros;

    if ( f_b == NULL )
        f_b = f_zeros;

    if ( f_c == NULL )
        f_c = f_zeros;

    if ( f_d == NULL )
        f_d = f_zeros;

    if ( f_0 == NULL )
        f_0 = f_zeros;

    if ( boundary == NULL )
        boundary = boundary_dirichlet;

    sprintf( buf, "%s.dat", FileName );
    fd = fopen( buf, "w" );

    y0 = malloc( sizeof( double ) * xN );
    y = malloc( sizeof( double ) * xN );
    x = malloc( sizeof( double ) * xN );

    switch ( type ) {
        case 0:
            dx =  (xmax-xmin) / ( xN-1  );
            for( i=0; i<xN; i++ )
                x[i] = xmin + i * dx;
            break;
        case 1:
            dx =  log10(xmax/xmin) / ( xN-1  );
            for( i=0; i<xN; i++ )
                x[i] = pow( 10, log10(xmin) + i * dx);
            break;
    }

    for( i=0; i<xN; i++ ) {
        y0[i] = (*f_0)(x[i], 0);
    }

    myprint( fd, "0 ", x, xN );
    myprint( fd, "0 ", y0, xN );

    p.y0 = y0;
    p.y = y;
    p.xN = xN;
    p.t_step = dt;
    p.f_a = f_a;
    p.f_b = f_b;
    p.f_c = f_c;
    p.f_d = f_d;
    p.t_start = t0;
    p.boundary = boundary;
    p.xmax = xmax;
    p.xmin = xmin;
    p.type = type;
    p.x = x;
    p.condition = NULL;

    CN_init( p );

    for( t=t0; t<t1; t+=dtt ) {

        CN_forward( dtt );

        sprintf( buf, "%g ", t+dtt );
        myprint( fd, buf, y, xN );

    }


    CN_free();
    free( y );
    fclose( fd );

    sprintf( buf, "%s_th.dat", FileName );
    fd = fopen( buf, "w" );

    fprintf( fd, "0 " );
    for( i=0; i<xN; i++ ) {
        fprintf( fd, "%g ", x[i] );
    }
    fprintf( fd, "\n" );

    fprintf( fd, "0 " );
    for( i=0; i<xN; i++ ) {
        fprintf( fd, "%g ", y0[i] );
    }
    free( y0 );
    fprintf( fd, "\n" );

    if ( f_th == NULL ) {
        free( x );
        fclose( fd );
        return;
    }

    for( t=t0; t<t1; t+=dtt ) {
        fprintf( fd, "%g ", t+dtt );
        for( i=0; i<xN; i++ ) {
            fprintf( fd, "%g ", (*f_th)( xmin+i*dx, t+dtt ) );
        }
        fprintf( fd, "\n" );
    }

    free( x );
    fclose( fd );

}

/***************************************************/
double f_b1_log( double x, double t ){
    return 1/exp(x);
}

double f_c1_log( double x, double t ){
    return -2;
}

void bound1_log( double *left, double *right, double *x, int N, double t ){

    double xx[4];
    int i;
    for( i=0; i<2; i++ ) {
        xx[i] = exp(x[i]);
        xx[3-i] = exp(x[N-1-i]);
    }

    left[2] = right[2] = 0;

    left[0] = -SQR(xx[0]) / (xx[1]-xx[0]) + 1 - x[0];
    left[1] = SQR(xx[0]) / (xx[1]-xx[0]);

    right[0] = SQR(xx[3]) / (xx[2]-xx[3]);
    right[1] = -SQR(xx[3]) / (xx[2]-xx[3]) + 1 - xx[3];

}

double f_a1( double x, double t ){
    return x*x;
}

double f_b1( double x, double t ){
    return x+1;
}

double f_c1( double x, double t ){
    return -2;
}

void bound1( double *left, double *right, double *x, int N, double t ){

    left[2] = right[2] = 0;

    left[0] = -SQR(x[0]) / (x[1]-x[0]) + 1 - x[0];
    left[1] = SQR(x[0]) / (x[1]-x[0]);

    right[0] = SQR(x[N-1]) / (x[N-2]-x[N-1]);
    right[1] = -SQR(x[N-1]) / (x[N-2]-x[N-1]) + 1 - x[N-1];

}

/***************************************************/

double f_b2_log( double x, double t ){
    return -1/exp(x);
}

double f_c2_log( double x, double t ){
    return -2;
}

void bound2_log( double *left, double *right, double *x, int N, double t ){

    double xx[4];
    int i;
    for( i=0; i<2; i++ ) {
        xx[i] = exp(x[i]);
        xx[3-i] = exp(x[N-1-i]);
    }

    left[2] = right[2] = 0;

    left[0] = -SQR(xx[0]) / (xx[1]-xx[0]) - 1 - x[0];
    left[1] = SQR(xx[0]) / (xx[1]-xx[0]);

    right[0] = SQR(xx[3]) / (xx[2]-xx[3]);
    right[1] = -SQR(xx[3]) / (xx[2]-xx[3]) - 1 - xx[3];

}

double f_a2( double x, double t ){
    return x*x;
}

double f_b2( double x, double t ){
    return x-1;
}

double f_c2( double x, double t ){
    return -2;
}

void bound2( double *left, double *right, double *x, int N, double t ){

    left[2] = right[2] = 0;

    left[0] = -SQR(x[0]) / (x[1]-x[0]) - 1 - x[0];
    left[1] = SQR(x[0]) / (x[1]-x[0]);

    right[0] = SQR(x[N-1]) / (x[N-2]-x[N-1]);
    right[1] = -SQR(x[N-1]) / (x[N-2]-x[N-1]) - 1 - x[N-1];

}

/***************************************************/

double f_c3_log( double x, double t ){
    return -(1+1/exp(x));
}

void bound3_log( double *left, double *right, double *x, int N, double t ){

    double xx[4];
    int i;
    for( i=0; i<2; i++ ) {
        xx[i] = exp(x[i]);
        xx[3-i] = exp(x[N-1-i]);
    }

    left[2] = right[2] = 0;

    left[0] = -SQR(xx[0]) / (xx[1]-xx[0]) - x[0];
    left[1] = SQR(xx[0]) / (xx[1]-xx[0]);

    right[0] = SQR(xx[3]) / (xx[2]-xx[3]);
    right[1] = -SQR(xx[3]) / (xx[2]-xx[3]) - xx[3];

}

double f_a3( double x, double t ){
    return x*x;
}

double f_b3( double x, double t ){
    return x;
}

double f_c3( double x, double t ){
    return -(1+1/x);
}

void bound3( double *left, double *right, double *x, int N, double t ){

    left[2] = right[2] = 0;

    left[0] = -SQR(x[0]) / (x[1]-x[0]) - x[0];
    left[1] = SQR(x[0]) / (x[1]-x[0]);

    right[0] = SQR(x[N-1]) / (x[N-2]-x[N-1]);
    right[1] = -SQR(x[N-1]) / (x[N-2]-x[N-1]) - x[N-1];

}

/***************************************************/
double f_a4_log( double x, double t ){
    return exp(x);
}

double f_b4_log( double x, double t ){
    return exp(x);
}

double f_c4_log( double x, double t ){
    return -(2*exp(x)+1);
}

void bound4_log( double *left, double *right, double *x, int N, double t ){

    double xx[4];
    int i;
    for( i=0; i<2; i++ ) {
        xx[i] = exp(x[i]);
        xx[3-i] = exp(x[N-1-i]);
    }

    left[2] = right[2] = 0;

    left[0] = -CUBE(xx[0]) / (xx[1]-xx[0]) - SQR(x[0]);
    left[1] = CUBE(xx[0]) / (xx[1]-xx[0]);

    right[0] = CUBE(xx[3]) / (xx[2]-xx[3]);
    right[1] = -CUBE(xx[3]) / (xx[2]-xx[3]) - SQR(xx[3]);
    //printf( "%g %g %g %g\n", left[0], left[1], right[0], right[1] );

}

double f_a4( double x, double t ){
    return x*x*x;
}

double f_b4( double x, double t ){
    return 2*x*x;
}

double f_c4( double x, double t ){
    return -2*x-1;
}

void bound4( double *left, double *right, double *x, int N, double t ){

    left[2] = right[2] = 0;

    left[0] = -CUBE(x[0]) / (x[1]-x[0]) - SQR(x[0]);
    left[1] = CUBE(x[0]) / (x[1]-x[0]);

    right[0] = CUBE(x[N-1]) / (x[N-2]-x[N-1]);
    right[1] = -CUBE(x[N-1]) / (x[N-2]-x[N-1]) - SQR(x[N-1]);
    //printf( "%g %g %g %g\n", left[0], left[1], right[0], right[1] );
}
/***************************************************/

double f_delta_log( double x, double t ){
    double a, x0, r;
    x = exp(x);
    a=0.01;
    x0 = 0.1;
    //printf( "[in d] %g, %g, %g, %g\n", x, r, (x-x0)/a, pow( (x-x0)/a, 2 ) );
    r =  1 / (a*sqrt(M_PI)) * exp( -pow( (x-x0)/a, 2 ) );
    return r;
}

double f_delta( double x, double t ){
    double a, x0, r;
    a=0.01;
    x0 = 0.1;
    //printf( "[in d] %g, %g, %g, %g\n", x, r, (x-x0)/a, pow( (x-x0)/a, 2 ) );
    r =  1 / (a*sqrt(M_PI)) * exp( -pow( (x-x0)/a, 2 ) );
    return r;
}

void main() {

    int i, j, model;

    dt = 1e-3;
    t0 = 0;
    t1 = 5;
    dtt = 1;
    f_th = NULL;
    xN = 1000;

    for( i=1; i<5; i++)
        for( j=0; j<3; j++ ) {
            /*
            if ( i!=4 )
                continue;
                */

            type = j;


            model = i*10 + j;
            printf( "Model_%i%i ...\n", i, j );
            sprintf( FileName, "Model_%i%i", i, j );
            switch ( model ){
                case 10:
                    f_a = f_ones;
                    f_b = f_b1_log;
                    f_c = f_c1_log;
                    f_d = f_delta_log;
                    f_0 = NULL;
                    xmin = -10;
                    xmax = 11;
                    boundary = boundary_newmann;
                    break;

                case 11:
                    f_a = f_a1;
                    f_b = f_b1;
                    f_c = f_c1;
                    f_d = f_delta;
                    xmin = 1e-5;
                    xmax = 1e5;
                    f_0 = NULL;
                    boundary = boundary_newmann;
                    break;

                case 20:
                    f_a = f_ones;
                    f_b = f_b2_log;
                    f_c = f_c2_log;
                    f_d = f_delta_log;
                    xmin = -10;
                    xmax = 11;
                    f_0 = NULL;
                    boundary = boundary_dirichlet;
                    break;

                case 21:
                    f_a = f_a2;
                    f_b = f_b2;
                    f_c = f_c2;
                    f_d = f_delta;
                    xmin = 1e-5;
                    xmax = 1e5;
                    f_0 = NULL;
                    boundary = boundary_dirichlet;
                    break;

                case 30:
                    f_a = f_ones;
                    f_b = NULL;
                    f_c = f_c3_log;
                    f_d = f_delta_log;
                    xmin = -10;
                    xmax = 11;
                    f_0 = NULL;
                    xN = 300;
                    boundary = boundary_newmann;
                    break;

                case 31:
                    f_a = f_a3;
                    f_b = f_b3;
                    f_c = f_c3;
                    f_d = f_delta;
                    xmin = 1e-5;
                    xmax = 1e5;
                    xN = 300;
                    f_0 = NULL;
                    boundary = boundary_newmann;
                    break;

                case 40:
                    f_a = f_a4_log;
                    f_b = f_b4_log;
                    f_c = f_c4_log;
                    f_d = NULL;
                    f_0 = f_delta_log;
                    xmin = -10;
                    xmax = 11;
                    boundary = boundary_newmann;
                    break;

                case 41:
                    f_a = f_a4;
                    f_b = f_b4;
                    f_c = f_c4;
                    f_d = NULL;
                    f_0 = f_delta;
                    xmin = 1e-5;
                    xmax = 1e5;
                    boundary = boundary_newmann;
                    break;

                }

                Crank_Nicolson_test();
        }

}


