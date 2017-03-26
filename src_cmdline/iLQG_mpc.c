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

double* assignParam(double* values, int numVal) {
    double* temp = malloc(numVal*sizeof(double));
    memcpy(temp, values, numVal*sizeof(double));

    return temp;
}

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

    double c_x = 116.0;
    double c_a = 197.0;
    double Iz = 0.025;
    double mu = 1.31;
    double mu_s = 0.5;

    double cu[2] = {1e-2, 1e-2};
    double cdu[2] = {1e-1, 15e-1};

    double cf[6] = {10, 10, 1, 0.1, 0.1, 0.1};
    double pf[6] = {0.01, 0.01, 0.1, 0.1, 0.1, 0.1};

    double cx[3] = {5e-2, 5e-2, 4e-2};
    double cdx[3] = {5e-1, 5e-2, 2e-2};
    double px[3] = {0.01, 0.01, 0.1};

    double cdrift = -0.001;

    double k_pos = 0.5;
    double k_vel = 0.1;
    double d_thres = 0.3;


    double limThr[2] = {0.0, 4.0};
    double limSteer[2] = {-0.68, 0.76};

    double xDes[6] = {3.0, 0.0, 0.0, 3.0, 0.0, 0.0};

    double Obs[2] = {1.0, 0.0};

    o->p= malloc(n_params*sizeof(double *));
    o->p[0] = assignParam(&G_f,1);
    o->p[1] = assignParam(&G_r,1);;
    o->p[2] = assignParam(&Iz,1);;
    o->p[3] = assignParam(Obs,2);
    o->p[4] = assignParam(&a,1);
    o->p[5] = assignParam(&b,1);
    o->p[6] = assignParam(&c_a,1);
    o->p[7] = assignParam(&c_x,1);
    o->p[8] = assignParam(&cdrift,1);
    o->p[9] = assignParam(cdu,2);
    o->p[10] = assignParam(cdx,3);
    o->p[11] = assignParam(cf,6);
    o->p[12] = assignParam(cu,2);
    o->p[13] = assignParam(cx,3);
    o->p[14] = assignParam(&d_thres,1);
    o->p[15] = assignParam(&h,1);
    o->p[16] = assignParam(&k_pos,1);
    o->p[17] = assignParam(&k_vel,1);
    o->p[18] = assignParam(limSteer,2);
    o->p[19] = assignParam(limThr,2);
    o->p[20] = assignParam(&m,1);
    o->p[21] = assignParam(&mu,1);
    o->p[22] = assignParam(&mu_s,1);
    o->p[23] = assignParam(pf,6);
    o->p[24] = assignParam(px,3);
    o->p[25] = assignParam(xDes,6);

    // printf("%f", o->p[3][0]);
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

    double x0[] = {0,0,0, 3,0,0,3,0,0,0};
    

    n= sizeof(x0)/sizeof(x0[0]);  // length of state vector
    m= 2; // number of inputs
    N= 51; // T+1 TODO how to make this variable?
    double u0[m*(N-1)]; // TODO row-first
    srand(time(NULL));
    for(int i=0; i<N-1; i++) {
        u0[i*m] = ((double)rand()/(double)(RAND_MAX)) * 0.5 + 3;
        u0[i*m+1] = ((double)rand()/(double)(RAND_MAX)) * 0.2 + 0.1;
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

    for(int i=0; i<n_params; i++) {
        free(o.p[i]);
    }

    free(o.p);
    for(i= 0; i<NUMBER_OF_THREADS+1; i++)
        free(o.trajectories[i].t);

    FILE *U = fopen("U.csv", "w");
    for(int i=0; i<N-1; i++) {
        fprintf(U, "%f, %f\n", u_new[i*m], u_new[i*m+1]);
    }
    fclose(U);

    FILE *X = fopen("X.csv", "w");
    for(int i=0; i<N; i++) {
        for(int j=0; j<n; j++) {
            fprintf(X, "%f,", x_new[i*n+j]);
        }
        fprintf(X,"\n");
    }
    fclose(X);
}
