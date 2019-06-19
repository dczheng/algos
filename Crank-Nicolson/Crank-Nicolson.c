#include "Crank-Nicolson.h"

struct CN_Data cn;

static double *y_cur, *y_later, *dx1, *dx2, t_cur, *A, *B, *C, *D;
static double left[3], right[3];

void CN_init( struct CN_Data p) {

    double dx;
    int i;
    cn = p;

    t_cur = cn.t_start;
    y_cur   = malloc( sizeof( double ) * cn.xN );
    y_later = malloc( sizeof( double ) * cn.xN );
    memcpy( y_cur, cn.y0, sizeof(double) * cn.xN );

    switch ( cn.type ) {

        case 0:
            printf( "use equal step\n" );
            dx = ( cn.xmax-cn.xmin ) / ( cn.xN-1 );
            cn.x = malloc( sizeof( double ) * cn.xN );
            for( i=0; i<cn.xN; i++ )
                cn.x[i] = cn.xmin + i * dx;
            break;

        case 1:
            printf( "use log step\n" );
            dx = log10(cn.xmax/cn.xmin) / ( cn.xN-1 );
            cn.x = malloc( sizeof( double ) * cn.xN );
            dx1 = malloc( sizeof( double ) * cn.xN );
            for( i=0; i<cn.xN; i++ )
                cn.x[i] = cn.xmin * pow( 10, i * dx);
            for( i=0; i<cn.xN-1; i++ )
                dx1[i] = cn.x[i+1] - cn.x[i];
            break;

        case 2:
            dx1 = malloc( sizeof( double ) * cn.xN );
            dx2 = malloc( sizeof( double ) * cn.xN );
            for( i=0; i<cn.xN-2; i++ ) {
                dx1[i] = cn.x[i+1] - cn.x[i];
                dx2[i] = cn.x[i+2] - cn.x[i];
            }
            dx1[cn.xN-2] = cn.x[cn.xN-1] - cn.x[cn.xN-2];
            break;

        default:
            printf( "ERROR" );
            exit(0);

    }

    A = malloc( sizeof(double)*cn.xN );
    B = malloc( sizeof(double)*cn.xN );
    C = malloc( sizeof(double)*cn.xN );
    D = malloc( sizeof(double)*cn.xN );


}

void print_ABCD() {

    int i;

    printf( "A: " );
    for( i=0; i<cn.xN; i++ )
        printf( "%g ", A[i]);
    printf( "\n" );

    printf( "B: " );
    for( i=0; i<cn.xN; i++ )
        printf( "%g ", B[i]);
    printf( "\n" );

    printf( "C: " );
    for( i=0; i<cn.xN; i++ )
        printf( "%g ", C[i]);
    printf( "\n" );

    printf( "D: " );
    for( i=0; i<cn.xN; i++ )
        printf( "%g ", D[i]);
    printf( "\n" );
}

void CN_single_forward() {

    int i;
    double a, b, c, d, tmid, h, k, h2, t1, t2;

    tmid = t_cur + cn.t_step / 2.0;

    switch ( cn.type ){

        case 0:

            k = cn.t_step;
            h = cn.x[1]-cn.x[0];
            h2 = h*h;

            for( i=1; i<cn.xN-1; i++ ) {

                a = (*cn.f_a)( cn.x[i], tmid );
                b = (*cn.f_b)( cn.x[i], tmid );
                c = (*cn.f_c)( cn.x[i], tmid );
                d = (*cn.f_d)( cn.x[i], tmid );

                A[i] = -( 2*k*a - k*h*b );
                B[i] = 4*h2 + 4*k*a - 2*h2*k*c;
                C[i] = -(2*k*a+k*h*b);

                D[i] = ( 2*k*a-k*h*b ) * y_cur[i-1]
                 + ( 4*h2-4*k*a+2*h2*k*c ) * y_cur[i]
                 + ( 2*k*a + k*h*b ) * y_cur[i+1]
                 + 4*h2*k*d;

            }

            break;

        case 1:
            for( i=1; i<cn.xN-1; i++ ) {

                k = cn.t_step;
                a = (*cn.f_a)( cn.x[i], tmid );
                b = (*cn.f_b)( cn.x[i], tmid );
                c = (*cn.f_c)( cn.x[i], tmid );
                d = (*cn.f_d)( cn.x[i], tmid );

                t1 = dx1[i] * dx1[i];
                t2 = dx1[i] * dx1[i-1];
                A[i] = -0.5 * a/t2;
                B[i] =  0.5 * (a/t1 + a/t2  + b / dx1[i] );
                C[i] =  0.5 * ( 2/k -  a/t1 - b/dx1[i] - c );
                D[i] =  0.5 * ( a/t2 * y_cur[i-1]
                                - (a/t1 + a/t2 + b/dx1[i]) * y_cur[i]
                                + ( a/t1 + b/dx1[i] + c + 2/k ) * y_cur[i+1] )
                        + d;
                if ( D[i] > 100 ) {
                    printf(
                            "i: %i\n"
                        "x: %g %g %g\n"
                        "a:%g, b:%g, c:%g, d:%g\n"
                        "t1:%g, t2:%g, dx1:%g, k:%g\n"
                        "y_cur: %g %g %g\n"
                        "D: %g\n",
                        i,
                        cn.x[i-1], cn.x[i], cn.x[i+1],
                        a, b, c, d,
                        t1, t2, dx1[i], k,
                        y_cur[i-1], y_cur[i], y_cur[i+1],
                        D[i]
                        );
                    exit(0);
                }

            }

            break;

        case 2:
            for( i=1; i<cn.xN-1; i++ ) {

                k = cn.t_step;
                a = (*cn.f_a)( cn.x[i], tmid );
                b = (*cn.f_b)( cn.x[i], tmid );
                c = (*cn.f_c)( cn.x[i], tmid );
                d = (*cn.f_d)( cn.x[i], tmid );

                A[i] = -a / (dx1[i-1]*dx2[i-1]) + b / (2*dx2[i-1]);
                B[i] = 1/k + a / (dx1[i]*dx1[i-1]) - c/2;
                C[i] = -a / ( dx1[i]*dx2[i-1]) - b / (2*dx2[i-1]);

                D[i] = ( a / (dx1[i-1]*dx2[i-1]) - b / (2*dx2[i-1]) ) * y_cur[i-1]
                    + ( 1/k - a / (dx1[i]*dx1[i-1]) + c / 2 ) * y_cur[i]
                    + ( a / (dx1[i]*dx2[i-1]) + b / (2*dx2[i-1]) ) * y_cur[i+1]
                    + d;

            }

            break;

    }

    (*cn.boundary)( left, right, cn.x, tmid );
    A[0] = C[cn.xN-1] = 0;
    B[0] = left[0];
    C[0] = left[1];
    D[0] = left[2];
    A[cn.xN-1] = right[0];
    B[cn.xN-1] = right[1];
    D[cn.xN-1] = right[2];
    TridMat( A, B, C, D, y_later, cn.xN, 1 );

    if ( cn.condition != NULL )
        (*cn.condition)( y_later, cn.xN );

    /*
    double flag;
    flag = 0;
    for ( i=0; i<cn.xN; i++ )
        if ( isnan(y_later[i]) ) {
            flag = 1;
            break;
        }
    if ( flag ) {
        print_ABCD();
        exit(0);
    }
    */
}

void CN_forward( double t_forward ) {

    //int i;
    double dt_tmp, *y_tmp, t_end;
    t_end = t_cur + t_forward;
    printf( "[%g] --> [%g]\n", t_cur, t_end );

    while( t_cur+ cn.t_step<=t_end ){
        CN_single_forward();;
        t_cur += cn.t_step;

        y_tmp = y_cur;
        y_cur = y_later;
        y_later = y_tmp;
    }

    if ( t_cur < t_end ){

        dt_tmp = cn.t_step;
        cn.t_step = t_end - t_cur;
        CN_single_forward();
        cn.t_step = dt_tmp;

        y_tmp = y_cur;
        y_cur = y_later;
        y_later = y_tmp;
    }

    t_cur = t_end;
    memcpy( cn.y, y_cur, sizeof( double ) * cn.xN );


}

void CN_free() {

    free( y_cur );
    free( y_later );

    switch ( cn.type ) {
        case 0:
            free(cn.x);
            break;
        case 1:
            free(cn.x);
            free(dx1);
            break;
        case 2:
            free( dx1 );
            free( dx2 );
            break;
    }

    free( A );
    free( B );
    free( C );
    free( D );

}

void TridMat( double *a, double *b, double *c, double *d, double *x, int N, int flag ) {

    int i;
    double *cc, *dd, w;

    if ( flag==1 ) {

        cc = malloc( sizeof( double ) * N );
        dd = malloc( sizeof( double ) * N );

        cc[0] = c[0] / b[0];
        dd[0]  = d[0] / b[0];

        for( i=1; i<N; i++ ) {

            cc[i] = c[i] / ( b[i] - a[i] * cc[i-1] );
            dd[i] = ( d[i] - a[i] * dd[i-1]  ) / ( b[i] - a[i] * cc[i-1] );

        }

        x[N-1] = dd[N-1];
        for ( i=N-2; i>=0; i-- ) {
            x[i] = dd[i] - cc[i] * x[i+1];
        }

        free( cc );
        free( dd );

    }
    else {

        for ( i=1; i<N; i++ ) {
            w = a[i] / b[i-1];
            //printf( "[%i] b: %g, d: %g,", i, b[i], d[i] );
            b[i] = b[i] - w * c[i-1];
            d[i] = d[i] - w * d[i-1];
            //printf( "w: %g, b: %g, d: %g\n", w, b[i], d[i]  );
        }

        x[N-1] = d[N-1] / b[N-1];
        for ( i=N-2; i>=0; i-- )
            x[i] = ( d[i] - c[i] * x[i+1]) / b[i];

    }

}
