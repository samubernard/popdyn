
/* include */
#include <stdlib.h>
#include <math.h>
#include "gaussrand.h"

/* function implementation */

double rand01()
{
    return (double)rand()/RAND_MAX; 
}

double randN(double mu, double sigma)
{
    double u, v, x, y, q;
    do {
        u = rand01();
        v = 1.7156*( rand01()-0.5 );
        x = u - 0.449871;
        y = fabs(v) + 0.386595;
        q = x*x + y*(0.19600*y-0.25472*x);
    } while ( q > 0.27597 && (q > 0.27846 || v*v > -4.0*log(u)*u*u) );
    return (mu + sigma*v/u);
}
