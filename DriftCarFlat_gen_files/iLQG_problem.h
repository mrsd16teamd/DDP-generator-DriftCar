/* File generated form template iLQG_problem.tem on 2017-03-22 11:44:54-04:00. Do not edit! */

#ifndef ILQG_PROBLEM_H
#define ILQG_PROBLEM_H

#include <math.h>
#include "mex.h"
#ifndef  HAVE_OCTAVE
#include "matrix.h"
#endif

#define isNANorINF(v) (mxIsNaN(v) || mxIsInf(v))
#define INF mxGetInf()
// #define isNANorINF(v) (isnan(v) || !finite(v))
// #define INF HUGE_VAL

#define N_X 10 
#define N_U 2 

#define sizeofQxx 55 
#define sizeofQuu 3 
#define sizeofQxu 20 

typedef struct {
    double x[N_X];
    double u[N_U];
    double lower[N_U];
    double upper[N_U];
    double lower_sign[N_U];
    double upper_sign[N_U];
    double lower_hx[N_X*N_U];
    double upper_hx[N_X*N_U];
    
    double l[N_U];
    double L[N_U*N_X]; 
    double c;
    double cx[N_X];
    double cxx[sizeofQxx];
    double cu[N_U];
    double cuu[sizeofQuu];
    double cxu[sizeofQxu];
    double fx[N_X*N_X];
    double fu[N_X*N_U];
#if FULL_DDP
    double fxx[N_X*sizeofQxx];
    double fuu[N_X*sizeofQuu];
    double fxu[N_X*sizeofQxu];
#endif
    double dy;
    double dx;
    double dr;
    double dUy;
    double dUx;

    double diff_dy_Ux;
    double diff_2dy_Ux_Ux;
    double diff_2dy_Ux_Uy;
    double diff_2dy_Ux_phi;
    double diff_dy_Uy;
    double diff_2dy_Uy_Uy;
    double diff_2dy_Uy_phi;
    double diff_dy_phi;
    double diff_2dy_phi_phi;
    double diff_dx_Ux;
    double diff_2dx_Ux_Ux;
    double diff_2dx_Ux_Uy;
    double diff_2dx_Ux_phi;
    double diff_dx_Uy;
    double diff_2dx_Uy_Uy;
    double diff_2dx_Uy_phi;
    double diff_dx_phi;
    double diff_2dx_phi_phi;
    double diff_dr_Ux;
    double diff_2dr_Ux_Ux;
    double diff_2dr_Ux_Uy;
    double diff_2dr_Ux_r;
    double diff_2dr_Ux_steer;
    double diff_2dr_Ux_thr;
    double diff_dr_Uy;
    double diff_2dr_Uy_Uy;
    double diff_2dr_Uy_r;
    double diff_2dr_Uy_steer;
    double diff_2dr_Uy_thr;
    double diff_dr_r;
    double diff_2dr_r_r;
    double diff_2dr_r_steer;
    double diff_2dr_r_thr;
    double diff_dr_steer;
    double diff_2dr_steer_steer;
    double diff_dr_thr;
    double diff_2dr_thr_thr;
    double diff_dUy_Ux;
    double diff_2dUy_Ux_Ux;
    double diff_2dUy_Ux_Uy;
    double diff_2dUy_Ux_r;
    double diff_2dUy_Ux_steer;
    double diff_2dUy_Ux_thr;
    double diff_dUy_Uy;
    double diff_2dUy_Uy_Uy;
    double diff_2dUy_Uy_r;
    double diff_2dUy_Uy_steer;
    double diff_2dUy_Uy_thr;
    double diff_dUy_r;
    double diff_2dUy_r_r;
    double diff_2dUy_r_steer;
    double diff_2dUy_r_thr;
    double diff_dUy_steer;
    double diff_2dUy_steer_steer;
    double diff_dUy_thr;
    double diff_2dUy_thr_thr;
    double diff_dUx_Ux;
    double diff_2dUx_Ux_Ux;
    double diff_2dUx_Ux_Uy;
    double diff_2dUx_Ux_r;
    double diff_2dUx_Ux_steer;
    double diff_2dUx_Ux_thr;
    double diff_dUx_Uy;
    double diff_2dUx_Uy_Uy;
    double diff_2dUx_Uy_r;
    double diff_2dUx_Uy_steer;
    double diff_2dUx_Uy_thr;
    double diff_dUx_r;
    double diff_2dUx_r_r;
    double diff_2dUx_r_steer;
    double diff_2dUx_r_thr;
    double diff_dUx_steer;
    double diff_2dUx_steer_steer;
    double diff_dUx_thr;
    double diff_2dUx_thr_thr;
} trajEl_t;

typedef struct {
    double x[N_X];

    double c;
    double cx[N_X];
    double cxx[sizeofQxx];


} trajFin_t;

typedef struct {
    trajEl_t* t;
    trajFin_t f;
} traj_t;

typedef struct {

} multipliersEl_t;

typedef struct {

} multipliersFin_t;

typedef struct {
    multipliersEl_t* t;
    multipliersFin_t f;
} multipliers_t;

#endif // ILQG_PROBLEM_H
