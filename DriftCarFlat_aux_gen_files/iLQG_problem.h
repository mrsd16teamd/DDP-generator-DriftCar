/* File generated form template iLQG_problem.tem on 2017-03-22 11:58:35-04:00. Do not edit! */

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
    double beta_;
    double dy;
    double U;
    double dx;
    double alpha_R;
    double gamma_R;
    double F_r;
    double Fyr;
    double dr;
    double K_;
    double K;
    double alpha;
    double alpha_F;
    double Fyf;
    double gamma_F;
    double F_f;
    double dUy;
    double dUx;
    double rev;
    double Fxr;

    double diff_beta__Ux;
    double diff_2beta__Ux_Ux;
    double diff_2beta__Ux_Uy;
    double diff_beta__Uy;
    double diff_2beta__Uy_Uy;
    double diff_dy_Ux;
    double diff_2dy_Ux_Ux;
    double diff_2dy_Ux_Uy;
    double diff_2dy_Ux_phi;
    double diff_dy_Uy;
    double diff_2dy_Uy_Uy;
    double diff_2dy_Uy_phi;
    double diff_dy_phi;
    double diff_2dy_phi_phi;
    double diff_U_Ux;
    double diff_2U_Ux_Ux;
    double diff_2U_Ux_Uy;
    double diff_U_Uy;
    double diff_2U_Uy_Uy;
    double diff_dx_Ux;
    double diff_2dx_Ux_Ux;
    double diff_2dx_Ux_Uy;
    double diff_2dx_Ux_phi;
    double diff_dx_Uy;
    double diff_2dx_Uy_Uy;
    double diff_2dx_Uy_phi;
    double diff_dx_phi;
    double diff_2dx_phi_phi;
    double diff_alpha_R_Ux;
    double diff_2alpha_R_Ux_Ux;
    double diff_2alpha_R_Ux_Uy;
    double diff_2alpha_R_Ux_r;
    double diff_alpha_R_Uy;
    double diff_2alpha_R_Uy_Uy;
    double diff_2alpha_R_Uy_r;
    double diff_alpha_R_r;
    double diff_2alpha_R_r_r;
    double diff_gamma_R_Ux;
    double diff_2gamma_R_Ux_Ux;
    double diff_2gamma_R_Ux_Uy;
    double diff_2gamma_R_Ux_r;
    double diff_2gamma_R_Ux_thr;
    double diff_gamma_R_Uy;
    double diff_2gamma_R_Uy_Uy;
    double diff_2gamma_R_Uy_r;
    double diff_2gamma_R_Uy_thr;
    double diff_gamma_R_r;
    double diff_2gamma_R_r_r;
    double diff_2gamma_R_r_thr;
    double diff_gamma_R_thr;
    double diff_2gamma_R_thr_thr;
    double diff_F_r_Ux;
    double diff_2F_r_Ux_Ux;
    double diff_2F_r_Ux_Uy;
    double diff_2F_r_Ux_r;
    double diff_2F_r_Ux_thr;
    double diff_F_r_Uy;
    double diff_2F_r_Uy_Uy;
    double diff_2F_r_Uy_r;
    double diff_2F_r_Uy_thr;
    double diff_F_r_r;
    double diff_2F_r_r_r;
    double diff_2F_r_r_thr;
    double diff_F_r_thr;
    double diff_2F_r_thr_thr;
    double diff_Fyr_Ux;
    double diff_2Fyr_Ux_Ux;
    double diff_2Fyr_Ux_Uy;
    double diff_2Fyr_Ux_r;
    double diff_2Fyr_Ux_thr;
    double diff_Fyr_Uy;
    double diff_2Fyr_Uy_Uy;
    double diff_2Fyr_Uy_r;
    double diff_2Fyr_Uy_thr;
    double diff_Fyr_r;
    double diff_2Fyr_r_r;
    double diff_2Fyr_r_thr;
    double diff_Fyr_thr;
    double diff_2Fyr_thr_thr;
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
    double diff_K__Ux;
    double diff_2K__Ux_Ux;
    double diff_2K__Ux_thr;
    double diff_K__thr;
    double diff_K_Ux;
    double diff_2K_Ux_Ux;
    double diff_2K_Ux_thr;
    double diff_K_thr;
    double diff_2K_thr_thr;
    double diff_alpha_Ux;
    double diff_2alpha_Ux_Ux;
    double diff_2alpha_Ux_Uy;
    double diff_2alpha_Ux_r;
    double diff_2alpha_Ux_steer;
    double diff_alpha_Uy;
    double diff_2alpha_Uy_Uy;
    double diff_2alpha_Uy_r;
    double diff_alpha_r;
    double diff_2alpha_r_r;
    double diff_alpha_steer;
    double diff_alpha_F_Ux;
    double diff_2alpha_F_Ux_Ux;
    double diff_2alpha_F_Ux_Uy;
    double diff_2alpha_F_Ux_r;
    double diff_2alpha_F_Ux_steer;
    double diff_alpha_F_Uy;
    double diff_2alpha_F_Uy_Uy;
    double diff_2alpha_F_Uy_r;
    double diff_2alpha_F_Uy_steer;
    double diff_alpha_F_r;
    double diff_2alpha_F_r_r;
    double diff_2alpha_F_r_steer;
    double diff_alpha_F_steer;
    double diff_2alpha_F_steer_steer;
    double diff_Fyf_Ux;
    double diff_2Fyf_Ux_Ux;
    double diff_2Fyf_Ux_Uy;
    double diff_2Fyf_Ux_r;
    double diff_2Fyf_Ux_steer;
    double diff_Fyf_Uy;
    double diff_2Fyf_Uy_Uy;
    double diff_2Fyf_Uy_r;
    double diff_2Fyf_Uy_steer;
    double diff_Fyf_r;
    double diff_2Fyf_r_r;
    double diff_2Fyf_r_steer;
    double diff_Fyf_steer;
    double diff_2Fyf_steer_steer;
    double diff_gamma_F_Ux;
    double diff_2gamma_F_Ux_Ux;
    double diff_2gamma_F_Ux_Uy;
    double diff_2gamma_F_Ux_r;
    double diff_gamma_F_Uy;
    double diff_2gamma_F_Uy_Uy;
    double diff_2gamma_F_Uy_r;
    double diff_gamma_F_r;
    double diff_2gamma_F_r_r;
    double diff_F_f_Ux;
    double diff_2F_f_Ux_Ux;
    double diff_2F_f_Ux_Uy;
    double diff_2F_f_Ux_r;
    double diff_F_f_Uy;
    double diff_2F_f_Uy_Uy;
    double diff_2F_f_Uy_r;
    double diff_F_f_r;
    double diff_2F_f_r_r;
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
    double diff_rev_Ux;
    double diff_2rev_Ux_Ux;
    double diff_2rev_Ux_thr;
    double diff_rev_thr;
    double diff_2rev_thr_thr;
    double diff_Fxr_Ux;
    double diff_2Fxr_Ux_Ux;
    double diff_2Fxr_Ux_Uy;
    double diff_2Fxr_Ux_r;
    double diff_2Fxr_Ux_thr;
    double diff_Fxr_Uy;
    double diff_2Fxr_Uy_Uy;
    double diff_2Fxr_Uy_r;
    double diff_2Fxr_Uy_thr;
    double diff_Fxr_r;
    double diff_2Fxr_r_r;
    double diff_2Fxr_r_thr;
    double diff_Fxr_thr;
    double diff_2Fxr_thr_thr;
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
