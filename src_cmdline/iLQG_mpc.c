// MATLAB Mex function wrapper for iLQG algorithm
// Copyright (c) 2016 Jens Geisler


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "iLQG.h"
#include "printMat.h"
#include "matMult.h"

void init_params(tOptSet *o)
{
    // Set model parameters
    double g = 9.81;
    double L = 0.257;

    double h = 0.02;
    double m = 2.35;
    double b = 0.14328;
    double a = L-b;
    double G_f = m*g*b/L;
    double G_r = m*g*a/L;

    double c_x = 116;
    double c_a = 197;
    double Iz = 0.025;
    double mu = 1.31;
    double mu_s = 0.5;

    double cu[] = {1e-2, 1e-2};
    double cdu[] = {1e-1, 15e-1};

    double cf[] = {10, 10, 1, 0.1, 0.1, 0.1};
    double pf[] = {0.01, 0.01, 0.1, 0.1, 0.1, 0.1};

    double cx[] = {5e-2, 5e-2, 4e-2};
    double cdx[] = {5e-1, 5e-2, 2e-2};
    double px[] = {0.01, 0.01, 0.1};

    double cdrift = -0.001;

    double k_pos = 0.1;
    double k_vel = 0;
    double d_thres = 0.3;


    double limThr[] = {0, 4};
    double limSteer[] = {-0.68, 0.76};

    double xDes[] = {3, 0, 0, 3, 0, 0};
    double Obs[] = {1, 0};

    o->p= malloc(n_params*sizeof(double *));
    o->p[0] = &G_f;
    o->p[1] = &G_r;
    o->p[2] = &Iz;
    o->p[3] = Obs;
    o->p[4] = &a;
    o->p[5] = &b;
    o->p[6] = &c_a;
    o->p[7] = &c_x;
    o->p[8] = &cdrift;
    o->p[9] = cdu;
    o->p[10] = cdx;
    o->p[11] = cf;
    o->p[12] = cu;
    o->p[13] = cx;
    o->p[14] = &d_thres;
    o->p[15] = &h;
    o->p[16] = &k_pos;
    o->p[17] = &k_vel;
    o->p[18] = limSteer;
    o->p[19] = limThr;
    o->p[20] = &m;
    o->p[21] = &mu;
    o->p[22] = &mu_s;
    o->p[23] = pf;
    o->p[24] = px;
    o->p[25] = xDes;
}

int main()
{
    //TODO make o manually

    // dims
    int N, n, m, m_, n_, si, i, k;
    // inputs
    tOptSet o= INIT_OPTSET;
    int dims[3];

    //outputs
    double *x_nom, *u_nom; //, *l, *L;

    // aux
    char *err_msg, *fname;
    clock_t begin, end;

    double x0[] = {0,0,0, 3,0,0, 3,0,0,0};
    

    n= sizeof(x0)/sizeof(x0[0]);  // length of state vector
    m= 2; // number of inputs
    N= 51; // T+1 TODO how to make this variable?
    double u0[m*(N-1)]; // TODO row-first
    srand(time(NULL));
    for(int i=0; i<N-1; i++) {
        u0[i*m] = ((double)rand()/(double)(RAND_MAX)) * 0.5 + 3;
        u0[i*m+1] = ((double)rand()/(double)(RAND_MAX)) * 0.2 + 0.3;
    }

    // inputs
    o.x0= x0; //double *
    u_nom= u0;  // double **
    o.n_hor= N-1;

    standard_parameters(&o);
    // Set optimization parameters
    fname = "max_iter";
    double max_iter = 100;

    err_msg = setOptParam(&o, fname, &max_iter, 1);
    if(err_msg) {
        printf("MATLAB:dimagree, Error setting optimization parameter '%s': %s.\n", fname, err_msg);
    }

    // Set model and problem parameters
    init_params(&o);

    // outputs
    double success[1];
    double new_cost[1];

    double x_new[n*N];
    double u_new[m*(N-1)];

    // aux
    for(i= 0; i<NUMBER_OF_THREADS+1; i++)
        o.trajectories[i].t= malloc(sizeof(trajEl_t)*(N-1));
    o.multipliers.t= malloc(sizeof(multipliersEl_t)*N);

    printf("Set const vars\n");
    if(!init_opt(&o)) {
        success[0]= 0;
        new_cost[0]= o.cost;
    } else {
        printf("Initializing trajectory\n");
        for(k= 0; k<N-1; k++)
            for(i= 0; i<N_U; i++)
                o.nominal->t[k].u[i]= u_nom[MAT_IDX(i, k, N_U)];
        if(!forward_pass(o.candidates[0], &o, 0.0, &o.cost, 0)) {
            printf("forward_pass failed\n");
            success[0]= 0;
            new_cost[0]= o.cost;
        } else {
            makeCandidateNominal(&o, 0);

            printf("Starting iLQG\n");
            begin = clock();
            success[0]= iLQG(&o);
            end = clock();
            printf("Time for iLQG: %f seconds\n", (double)(end - begin) / CLOCKS_PER_SEC);
            for(k= 0; k<N-1; k++)
                for(i= 0; i<N_X; i++)
                    x_new[MAT_IDX(i, k, N_X)]= o.nominal->t[k].x[i];
            for(i= 0; i<N_X; i++)
                x_new[MAT_IDX(i, N-1, N_X)]= o.nominal->f.x[i];

            for(k= 0; k<N-1; k++)
                for(i= 0; i<N_U; i++)
                    u_new[MAT_IDX(i, k, N_U)]= o.nominal->t[k].u[i];

            new_cost[0]= o.cost;
        }
    }

    free(o.p);
    for(i= 0; i<NUMBER_OF_THREADS+1; i++)
        free(o.trajectories[i].t);

    FILE *U = fopen("U.csv", "w");
    for(int i=0; i<N-1; i++) {
        fprintf(U, "%f, %f\n", u_new[i*m], u_new[i*m+1]);
    }
    fclose(U);
}
