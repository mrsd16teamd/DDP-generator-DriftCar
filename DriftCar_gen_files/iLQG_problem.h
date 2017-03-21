/* File generated form template iLQG_problem.tem on 2017-03-21 12:19:40-04:00. Do not edit! */

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
    double s6;
    double s5;
    double s4;
    double s3;
    double s2;
    double s1;

    double diff_s6_Ux;
    double diff_2s6_Ux_Ux;
    double diff_2s6_Ux_Uy;
    double diff_2s6_Ux_r;
    double diff_2s6_Ux_steer;
    double diff_2s6_Ux_thr;
    double diff_s6_Uy;
    double diff_2s6_Uy_Uy;
    double diff_2s6_Uy_r;
    double diff_2s6_Uy_steer;
    double diff_2s6_Uy_thr;
    double diff_s6_r;
    double diff_2s6_r_r;
    double diff_2s6_r_steer;
    double diff_2s6_r_thr;
    double diff_s6_steer;
    double diff_2s6_steer_steer;
    double diff_s6_thr;
    double diff_2s6_thr_thr;
    double diff_s5_Ux;
    double diff_2s5_Ux_Ux;
    double diff_2s5_Ux_Uy;
    double diff_2s5_Ux_r;
    double diff_2s5_Ux_steer;
    double diff_2s5_Ux_thr;
    double diff_s5_Uy;
    double diff_2s5_Uy_Uy;
    double diff_2s5_Uy_r;
    double diff_2s5_Uy_steer;
    double diff_2s5_Uy_thr;
    double diff_s5_r;
    double diff_2s5_r_r;
    double diff_2s5_r_steer;
    double diff_2s5_r_thr;
    double diff_s5_steer;
    double diff_2s5_steer_steer;
    double diff_s5_thr;
    double diff_2s5_thr_thr;
    double diff_s4_Ux;
    double diff_2s4_Ux_Ux;
    double diff_2s4_Ux_Uy;
    double diff_2s4_Ux_r;
    double diff_2s4_Ux_steer;
    double diff_2s4_Ux_thr;
    double diff_s4_Uy;
    double diff_2s4_Uy_Uy;
    double diff_2s4_Uy_r;
    double diff_2s4_Uy_steer;
    double diff_2s4_Uy_thr;
    double diff_s4_r;
    double diff_2s4_r_r;
    double diff_2s4_r_steer;
    double diff_2s4_r_thr;
    double diff_s4_steer;
    double diff_2s4_steer_steer;
    double diff_s4_thr;
    double diff_2s4_thr_thr;
    double diff_s2_Ux;
    double diff_2s2_Ux_Ux;
    double diff_2s2_Ux_Uy;
    double diff_2s2_Ux_phi;
    double diff_s2_Uy;
    double diff_2s2_Uy_Uy;
    double diff_2s2_Uy_phi;
    double diff_s2_phi;
    double diff_2s2_phi_phi;
    double diff_s1_Ux;
    double diff_2s1_Ux_Ux;
    double diff_2s1_Ux_Uy;
    double diff_2s1_Ux_phi;
    double diff_s1_Uy;
    double diff_2s1_Uy_Uy;
    double diff_2s1_Uy_phi;
    double diff_s1_phi;
    double diff_2s1_phi_phi;
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
