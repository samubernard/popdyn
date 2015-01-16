/* coupled_oscillators.c Kuramoto model of coupled phase oscillators */
/* 
 * Dynamique Populations Cellulaires
 * Université Lyon 1, 2014
 *
 * Fichiers:
 *  coupled_oscillators.c
 *  makefile
 *
 * Ce code utilise la librairie gsl (GNU Scientific Library)
 *
 * Pour compiler sur mac ou linux:
 *
 * ouvrir un terminal 
 * >> make
 * >> ./coupled_oscillators
 * 
 * Deux fichiers sont générés en sortie:
 *  osc.txt : series temporelles
 *  order.txt: paramètre d'ordre
 *
 * Le fichier osc.txt contient est structuré en colonnes::
 * 0   w1  w2 ...      wNOSC      
 * t1  1   2   3   ...  NOSC
 * 
 * La premiere ligne est composee des frequences des oscillateurs. Les
 * oscillateurs sont places en ordre croissant de frequence.
 * t est le temps, r le module du paramètre d'ordre, psi l'argument du paramètre
 * d'ordre, 1...NOSC les phases des oscillateur.
 * Les phases des oscillateurs sont dans R (pas dans [0, 2*pi]).
 *
 * Pour visualiser les résultats, on peut utiliser gnuplot ou Matlab
 *
 * GNUPLOT:
 * plot 'osc.txt' u 1:2:(sin($3)) nonuniform matrix with image
 * 
 * Pour tracer le paramètre d'ordre
 *
 * pour r:
 * plot 'order.txt' using 1:2 w lines
 *
 * pour psi:
 * plot 'order.txt' using 1:3 w lines
 *
 * Pour Matlab, utiliser le script coupled_oscillators.m
 * 
 */

/* ================================================================= */
/*                              Libraries                            */
/* ================================================================= */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_odeiv2.h>

struct par {
    uint32_t NOSC;
    double K;
    double *w;
};

int ode_rhs(double t, const double y[], double f[], void *params);
int compare(const void *, const void *);

int main ( int argc, char *argv[] )
{
    uint32_t NOSC=50; /* nbr oscillators */
    double K=.5;  /* coupling stength */
    double sigma=0.5; /* std w distribution */
    double *y;     /* phase of oscillators */
    struct par mu; /* parameters */
    FILE *file_order,*file_out;
    uint32_t i;
    double t,tfinal;
    double h,hmin;
    double w_mean;
    double r,psi, rx, ry;

	/* get command line input parameters */
	if (argc == 4)
    {
        NOSC = strtoul(argv[1],NULL,10);	
    	K = strtod(argv[2],NULL);
        tfinal = strtod(argv[3],NULL);
    }	
	else
    {
    	fprintf(stderr, "Error: call to coupled_oscillators must have the following three mandatory arguments: unsigned long NOSC, double K, and double tfinal\n");
    	exit(EXIT_FAILURE);
    }

	/* init random number generator */
	const gsl_rng_type * rngT;
	gsl_rng * rand_nbr_gen;
	gsl_rng_env_setup();
	rngT = gsl_rng_default;
	rand_nbr_gen = gsl_rng_alloc(rngT);

    /* gsl ode */
    const gsl_odeiv2_step_type * odeT;
    gsl_odeiv2_step * s;
    gsl_odeiv2_control * c;
    gsl_odeiv2_evolve * e;
    gsl_odeiv2_system sys; 
    int status;
    
    mu.NOSC = NOSC;
    mu.K = K;
    mu.w = malloc(NOSC*sizeof(double));
    y = malloc(NOSC*sizeof(double));
    w_mean = 0.0;
    for (i = 0; i < NOSC; i++)
    {
        y[i] = 2*M_PI*gsl_rng_uniform(rand_nbr_gen)-M_PI;
        mu.w[i] = gsl_ran_gaussian(rand_nbr_gen,sigma*sigma);
        w_mean += mu.w[i];
    }
    w_mean /= NOSC;
    
    /* sort the frequencies from slowest to fastest */
    qsort(mu.w, NOSC, sizeof(double), compare);

    /* remove the mean to mu.w */
    file_out = fopen("osc.txt","w");
    fprintf(file_out,"%.5e ",0.0);
    for ( i = 0; i < NOSC; i++)
    {
        mu.w[i] -= w_mean;
        fprintf(file_out,"%.5e ",mu.w[i]);
    }
    fprintf(file_out,"\n");
    /* DO NOT TOUCH mu.w BELOW */

    /* initialise the ode solver */
    /* the solver used is a Runge-Kutta-Fehlberg (4, 5) method */
    odeT = gsl_odeiv2_step_rkf45;
    s = gsl_odeiv2_step_alloc(odeT,NOSC);
    c = gsl_odeiv2_control_y_new(1e-6,0.0);
    e = gsl_odeiv2_evolve_alloc(NOSC);
    
    h = 1e-2;
    hmin = 1e-5;
    t = 0.0;

    file_order = fopen("order.txt","w");
    while (t < tfinal )
    {

        /* apply one step of the ode solver */
        sys = (gsl_odeiv2_system) {ode_rhs, NULL, NOSC, &mu};
        status = gsl_odeiv2_evolve_apply(e,c,s,&sys,&t,tfinal,&h,y);
        if (status != GSL_SUCCESS)
            break;

        rx = 0;
        ry = 0;
        for ( i = 0; i < NOSC; i++)
        {
            rx += cos(y[i]);
            ry += sin(y[i]);
            
        }
        rx /= NOSC;
        ry /= NOSC;
        r = hypot(rx,ry);
        psi = atan2(ry,rx);
        
        if (h < hmin) h = hmin;
        printf(".\n");

        fprintf(file_order,"%.5e %.5e %.5e\n",t,r,psi); 
        fprintf(file_out,"%.5e ",t);

        for (i = 0; i < NOSC; i++)
        {
            fprintf (file_out,"%.5e ",y[i]); 
        }
        fprintf(file_out,"\n");
        

    }

    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);

    fclose(file_out);
    fclose(file_order);

    free(y);
    free(mu.w);

	return 0;
}


int ode_rhs(double t, const double y[], double f[], void *params)
{

    struct par mu = *(struct par*)params;
    uint32_t NOSC = mu.NOSC;
    double K  = mu.K;
    double *w = mu.w;
    uint32_t i,j;
    double *c;

    /* printf("call to ode_rhs, NOSC=%d\n",NOSC); */
    c = malloc(NOSC*sizeof(double));
    /* compute coupling */
    for ( i = 0; i < NOSC; i++)
    {
        c[i] = 0.0;
        for ( j = 0; j < NOSC; j++)
        {
            c[i] += sin(y[j]-y[i]);
        }

        f[i] = w[i] + K/NOSC*c[i];
        /* printf("f[%d]=%f\n",i,f[i]); */

    }

    free(c);

    return GSL_SUCCESS;

}

int compare(const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;

    return (*da > *db) - (*da < *db);
}
