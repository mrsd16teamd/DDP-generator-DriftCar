/* File generated form template iLQG_func.tem on 2017-03-22 11:58:35-04:00. Do not edit! */

#include "iLQG.h"
#include "matMult.h"

#define mcond(cond, a, dummy, b) ((cond)? a: b)
#define abs fabs
#define sec(x) (1/cos(x))


int n_params= 26;
    
tParamDesc p_name1= {"G_f", 1, 0}; // -> p[0][0]
tParamDesc p_name2= {"G_r", 1, 0}; // -> p[1][0]
tParamDesc p_name3= {"Iz", 1, 0}; // -> p[2][0]
tParamDesc p_name4= {"Obs", 2, 0}; // -> p[3][0]
tParamDesc p_name5= {"a", 1, 0}; // -> p[4][0]
tParamDesc p_name6= {"b", 1, 0}; // -> p[5][0]
tParamDesc p_name7= {"c_a", 1, 0}; // -> p[6][0]
tParamDesc p_name8= {"c_x", 1, 0}; // -> p[7][0]
tParamDesc p_name9= {"cdrift", 1, 0}; // -> p[8][0]
tParamDesc p_name10= {"cdu", 2, 0}; // -> p[9][0]
tParamDesc p_name11= {"cdx", 3, 0}; // -> p[10][0]
tParamDesc p_name12= {"cf", 6, 0}; // -> p[11][0]
tParamDesc p_name13= {"cu", 2, 0}; // -> p[12][0]
tParamDesc p_name14= {"cx", 3, 0}; // -> p[13][0]
tParamDesc p_name15= {"d_thres", 1, 0}; // -> p[14][0]
tParamDesc p_name16= {"h", 1, 0}; // -> p[15][0]
tParamDesc p_name17= {"k_pos", 1, 0}; // -> p[16][0]
tParamDesc p_name18= {"k_vel", 1, 0}; // -> p[17][0]
tParamDesc p_name19= {"limSteer", 2, 0}; // -> p[18][0]
tParamDesc p_name20= {"limThr", 2, 0}; // -> p[19][0]
tParamDesc p_name21= {"m", 1, 0}; // -> p[20][0]
tParamDesc p_name22= {"mu", 1, 0}; // -> p[21][0]
tParamDesc p_name23= {"mu_s", 1, 0}; // -> p[22][0]
tParamDesc p_name24= {"pf", 6, 0}; // -> p[23][0]
tParamDesc p_name25= {"px", 3, 0}; // -> p[24][0]
tParamDesc p_name26= {"xDes", 6, 0}; // -> p[25][0]
int n_vars= 0;
    
tParamDesc *paramdesc[]= {&p_name1, &p_name2, &p_name3, &p_name4, &p_name5, &p_name6, &p_name7, &p_name8, &p_name9, &p_name10, &p_name11, &p_name12, &p_name13, &p_name14, &p_name15, &p_name16, &p_name17, &p_name18, &p_name19, &p_name20, &p_name21, &p_name22, &p_name23, &p_name24, &p_name25, &p_name26};

#define aux_beta_ t->beta_
#define aux_dy t->dy
#define aux_U t->U
#define aux_dx t->dx
#define aux_alpha_R t->alpha_R
#define aux_gamma_R t->gamma_R
#define aux_F_r t->F_r
#define aux_Fyr t->Fyr
#define aux_dr t->dr
#define aux_K_ t->K_
#define aux_K t->K
#define aux_alpha t->alpha
#define aux_alpha_F t->alpha_F
#define aux_Fyf t->Fyf
#define aux_gamma_F t->gamma_F
#define aux_F_f t->F_f
#define aux_dUy t->dUy
#define aux_dUx t->dUx
#define aux_rev t->rev
#define aux_Fxr t->Fxr

#define daux_diff_beta__Ux t->diff_beta__Ux
#define daux_diff_2beta__Ux_Ux t->diff_2beta__Ux_Ux
#define daux_diff_2beta__Ux_Uy t->diff_2beta__Ux_Uy
#define daux_diff_beta__Uy t->diff_beta__Uy
#define daux_diff_2beta__Uy_Uy t->diff_2beta__Uy_Uy
#define daux_diff_dy_Ux t->diff_dy_Ux
#define daux_diff_2dy_Ux_Ux t->diff_2dy_Ux_Ux
#define daux_diff_2dy_Ux_Uy t->diff_2dy_Ux_Uy
#define daux_diff_2dy_Ux_phi t->diff_2dy_Ux_phi
#define daux_diff_dy_Uy t->diff_dy_Uy
#define daux_diff_2dy_Uy_Uy t->diff_2dy_Uy_Uy
#define daux_diff_2dy_Uy_phi t->diff_2dy_Uy_phi
#define daux_diff_dy_phi t->diff_dy_phi
#define daux_diff_2dy_phi_phi t->diff_2dy_phi_phi
#define daux_diff_U_Ux t->diff_U_Ux
#define daux_diff_2U_Ux_Ux t->diff_2U_Ux_Ux
#define daux_diff_2U_Ux_Uy t->diff_2U_Ux_Uy
#define daux_diff_U_Uy t->diff_U_Uy
#define daux_diff_2U_Uy_Uy t->diff_2U_Uy_Uy
#define daux_diff_dx_Ux t->diff_dx_Ux
#define daux_diff_2dx_Ux_Ux t->diff_2dx_Ux_Ux
#define daux_diff_2dx_Ux_Uy t->diff_2dx_Ux_Uy
#define daux_diff_2dx_Ux_phi t->diff_2dx_Ux_phi
#define daux_diff_dx_Uy t->diff_dx_Uy
#define daux_diff_2dx_Uy_Uy t->diff_2dx_Uy_Uy
#define daux_diff_2dx_Uy_phi t->diff_2dx_Uy_phi
#define daux_diff_dx_phi t->diff_dx_phi
#define daux_diff_2dx_phi_phi t->diff_2dx_phi_phi
#define daux_diff_alpha_R_Ux t->diff_alpha_R_Ux
#define daux_diff_2alpha_R_Ux_Ux t->diff_2alpha_R_Ux_Ux
#define daux_diff_2alpha_R_Ux_Uy t->diff_2alpha_R_Ux_Uy
#define daux_diff_2alpha_R_Ux_r t->diff_2alpha_R_Ux_r
#define daux_diff_alpha_R_Uy t->diff_alpha_R_Uy
#define daux_diff_2alpha_R_Uy_Uy t->diff_2alpha_R_Uy_Uy
#define daux_diff_2alpha_R_Uy_r t->diff_2alpha_R_Uy_r
#define daux_diff_alpha_R_r t->diff_alpha_R_r
#define daux_diff_2alpha_R_r_r t->diff_2alpha_R_r_r
#define daux_diff_gamma_R_Ux t->diff_gamma_R_Ux
#define daux_diff_2gamma_R_Ux_Ux t->diff_2gamma_R_Ux_Ux
#define daux_diff_2gamma_R_Ux_Uy t->diff_2gamma_R_Ux_Uy
#define daux_diff_2gamma_R_Ux_r t->diff_2gamma_R_Ux_r
#define daux_diff_2gamma_R_Ux_thr t->diff_2gamma_R_Ux_thr
#define daux_diff_gamma_R_Uy t->diff_gamma_R_Uy
#define daux_diff_2gamma_R_Uy_Uy t->diff_2gamma_R_Uy_Uy
#define daux_diff_2gamma_R_Uy_r t->diff_2gamma_R_Uy_r
#define daux_diff_2gamma_R_Uy_thr t->diff_2gamma_R_Uy_thr
#define daux_diff_gamma_R_r t->diff_gamma_R_r
#define daux_diff_2gamma_R_r_r t->diff_2gamma_R_r_r
#define daux_diff_2gamma_R_r_thr t->diff_2gamma_R_r_thr
#define daux_diff_gamma_R_thr t->diff_gamma_R_thr
#define daux_diff_2gamma_R_thr_thr t->diff_2gamma_R_thr_thr
#define daux_diff_F_r_Ux t->diff_F_r_Ux
#define daux_diff_2F_r_Ux_Ux t->diff_2F_r_Ux_Ux
#define daux_diff_2F_r_Ux_Uy t->diff_2F_r_Ux_Uy
#define daux_diff_2F_r_Ux_r t->diff_2F_r_Ux_r
#define daux_diff_2F_r_Ux_thr t->diff_2F_r_Ux_thr
#define daux_diff_F_r_Uy t->diff_F_r_Uy
#define daux_diff_2F_r_Uy_Uy t->diff_2F_r_Uy_Uy
#define daux_diff_2F_r_Uy_r t->diff_2F_r_Uy_r
#define daux_diff_2F_r_Uy_thr t->diff_2F_r_Uy_thr
#define daux_diff_F_r_r t->diff_F_r_r
#define daux_diff_2F_r_r_r t->diff_2F_r_r_r
#define daux_diff_2F_r_r_thr t->diff_2F_r_r_thr
#define daux_diff_F_r_thr t->diff_F_r_thr
#define daux_diff_2F_r_thr_thr t->diff_2F_r_thr_thr
#define daux_diff_Fyr_Ux t->diff_Fyr_Ux
#define daux_diff_2Fyr_Ux_Ux t->diff_2Fyr_Ux_Ux
#define daux_diff_2Fyr_Ux_Uy t->diff_2Fyr_Ux_Uy
#define daux_diff_2Fyr_Ux_r t->diff_2Fyr_Ux_r
#define daux_diff_2Fyr_Ux_thr t->diff_2Fyr_Ux_thr
#define daux_diff_Fyr_Uy t->diff_Fyr_Uy
#define daux_diff_2Fyr_Uy_Uy t->diff_2Fyr_Uy_Uy
#define daux_diff_2Fyr_Uy_r t->diff_2Fyr_Uy_r
#define daux_diff_2Fyr_Uy_thr t->diff_2Fyr_Uy_thr
#define daux_diff_Fyr_r t->diff_Fyr_r
#define daux_diff_2Fyr_r_r t->diff_2Fyr_r_r
#define daux_diff_2Fyr_r_thr t->diff_2Fyr_r_thr
#define daux_diff_Fyr_thr t->diff_Fyr_thr
#define daux_diff_2Fyr_thr_thr t->diff_2Fyr_thr_thr
#define daux_diff_dr_Ux t->diff_dr_Ux
#define daux_diff_2dr_Ux_Ux t->diff_2dr_Ux_Ux
#define daux_diff_2dr_Ux_Uy t->diff_2dr_Ux_Uy
#define daux_diff_2dr_Ux_r t->diff_2dr_Ux_r
#define daux_diff_2dr_Ux_steer t->diff_2dr_Ux_steer
#define daux_diff_2dr_Ux_thr t->diff_2dr_Ux_thr
#define daux_diff_dr_Uy t->diff_dr_Uy
#define daux_diff_2dr_Uy_Uy t->diff_2dr_Uy_Uy
#define daux_diff_2dr_Uy_r t->diff_2dr_Uy_r
#define daux_diff_2dr_Uy_steer t->diff_2dr_Uy_steer
#define daux_diff_2dr_Uy_thr t->diff_2dr_Uy_thr
#define daux_diff_dr_r t->diff_dr_r
#define daux_diff_2dr_r_r t->diff_2dr_r_r
#define daux_diff_2dr_r_steer t->diff_2dr_r_steer
#define daux_diff_2dr_r_thr t->diff_2dr_r_thr
#define daux_diff_dr_steer t->diff_dr_steer
#define daux_diff_2dr_steer_steer t->diff_2dr_steer_steer
#define daux_diff_dr_thr t->diff_dr_thr
#define daux_diff_2dr_thr_thr t->diff_2dr_thr_thr
#define daux_diff_K__Ux t->diff_K__Ux
#define daux_diff_2K__Ux_Ux t->diff_2K__Ux_Ux
#define daux_diff_2K__Ux_thr t->diff_2K__Ux_thr
#define daux_diff_K__thr t->diff_K__thr
#define daux_diff_K_Ux t->diff_K_Ux
#define daux_diff_2K_Ux_Ux t->diff_2K_Ux_Ux
#define daux_diff_2K_Ux_thr t->diff_2K_Ux_thr
#define daux_diff_K_thr t->diff_K_thr
#define daux_diff_2K_thr_thr t->diff_2K_thr_thr
#define daux_diff_alpha_Ux t->diff_alpha_Ux
#define daux_diff_2alpha_Ux_Ux t->diff_2alpha_Ux_Ux
#define daux_diff_2alpha_Ux_Uy t->diff_2alpha_Ux_Uy
#define daux_diff_2alpha_Ux_r t->diff_2alpha_Ux_r
#define daux_diff_2alpha_Ux_steer t->diff_2alpha_Ux_steer
#define daux_diff_alpha_Uy t->diff_alpha_Uy
#define daux_diff_2alpha_Uy_Uy t->diff_2alpha_Uy_Uy
#define daux_diff_2alpha_Uy_r t->diff_2alpha_Uy_r
#define daux_diff_alpha_r t->diff_alpha_r
#define daux_diff_2alpha_r_r t->diff_2alpha_r_r
#define daux_diff_alpha_steer t->diff_alpha_steer
#define daux_diff_alpha_F_Ux t->diff_alpha_F_Ux
#define daux_diff_2alpha_F_Ux_Ux t->diff_2alpha_F_Ux_Ux
#define daux_diff_2alpha_F_Ux_Uy t->diff_2alpha_F_Ux_Uy
#define daux_diff_2alpha_F_Ux_r t->diff_2alpha_F_Ux_r
#define daux_diff_2alpha_F_Ux_steer t->diff_2alpha_F_Ux_steer
#define daux_diff_alpha_F_Uy t->diff_alpha_F_Uy
#define daux_diff_2alpha_F_Uy_Uy t->diff_2alpha_F_Uy_Uy
#define daux_diff_2alpha_F_Uy_r t->diff_2alpha_F_Uy_r
#define daux_diff_2alpha_F_Uy_steer t->diff_2alpha_F_Uy_steer
#define daux_diff_alpha_F_r t->diff_alpha_F_r
#define daux_diff_2alpha_F_r_r t->diff_2alpha_F_r_r
#define daux_diff_2alpha_F_r_steer t->diff_2alpha_F_r_steer
#define daux_diff_alpha_F_steer t->diff_alpha_F_steer
#define daux_diff_2alpha_F_steer_steer t->diff_2alpha_F_steer_steer
#define daux_diff_Fyf_Ux t->diff_Fyf_Ux
#define daux_diff_2Fyf_Ux_Ux t->diff_2Fyf_Ux_Ux
#define daux_diff_2Fyf_Ux_Uy t->diff_2Fyf_Ux_Uy
#define daux_diff_2Fyf_Ux_r t->diff_2Fyf_Ux_r
#define daux_diff_2Fyf_Ux_steer t->diff_2Fyf_Ux_steer
#define daux_diff_Fyf_Uy t->diff_Fyf_Uy
#define daux_diff_2Fyf_Uy_Uy t->diff_2Fyf_Uy_Uy
#define daux_diff_2Fyf_Uy_r t->diff_2Fyf_Uy_r
#define daux_diff_2Fyf_Uy_steer t->diff_2Fyf_Uy_steer
#define daux_diff_Fyf_r t->diff_Fyf_r
#define daux_diff_2Fyf_r_r t->diff_2Fyf_r_r
#define daux_diff_2Fyf_r_steer t->diff_2Fyf_r_steer
#define daux_diff_Fyf_steer t->diff_Fyf_steer
#define daux_diff_2Fyf_steer_steer t->diff_2Fyf_steer_steer
#define daux_diff_gamma_F_Ux t->diff_gamma_F_Ux
#define daux_diff_2gamma_F_Ux_Ux t->diff_2gamma_F_Ux_Ux
#define daux_diff_2gamma_F_Ux_Uy t->diff_2gamma_F_Ux_Uy
#define daux_diff_2gamma_F_Ux_r t->diff_2gamma_F_Ux_r
#define daux_diff_gamma_F_Uy t->diff_gamma_F_Uy
#define daux_diff_2gamma_F_Uy_Uy t->diff_2gamma_F_Uy_Uy
#define daux_diff_2gamma_F_Uy_r t->diff_2gamma_F_Uy_r
#define daux_diff_gamma_F_r t->diff_gamma_F_r
#define daux_diff_2gamma_F_r_r t->diff_2gamma_F_r_r
#define daux_diff_F_f_Ux t->diff_F_f_Ux
#define daux_diff_2F_f_Ux_Ux t->diff_2F_f_Ux_Ux
#define daux_diff_2F_f_Ux_Uy t->diff_2F_f_Ux_Uy
#define daux_diff_2F_f_Ux_r t->diff_2F_f_Ux_r
#define daux_diff_F_f_Uy t->diff_F_f_Uy
#define daux_diff_2F_f_Uy_Uy t->diff_2F_f_Uy_Uy
#define daux_diff_2F_f_Uy_r t->diff_2F_f_Uy_r
#define daux_diff_F_f_r t->diff_F_f_r
#define daux_diff_2F_f_r_r t->diff_2F_f_r_r
#define daux_diff_dUy_Ux t->diff_dUy_Ux
#define daux_diff_2dUy_Ux_Ux t->diff_2dUy_Ux_Ux
#define daux_diff_2dUy_Ux_Uy t->diff_2dUy_Ux_Uy
#define daux_diff_2dUy_Ux_r t->diff_2dUy_Ux_r
#define daux_diff_2dUy_Ux_steer t->diff_2dUy_Ux_steer
#define daux_diff_2dUy_Ux_thr t->diff_2dUy_Ux_thr
#define daux_diff_dUy_Uy t->diff_dUy_Uy
#define daux_diff_2dUy_Uy_Uy t->diff_2dUy_Uy_Uy
#define daux_diff_2dUy_Uy_r t->diff_2dUy_Uy_r
#define daux_diff_2dUy_Uy_steer t->diff_2dUy_Uy_steer
#define daux_diff_2dUy_Uy_thr t->diff_2dUy_Uy_thr
#define daux_diff_dUy_r t->diff_dUy_r
#define daux_diff_2dUy_r_r t->diff_2dUy_r_r
#define daux_diff_2dUy_r_steer t->diff_2dUy_r_steer
#define daux_diff_2dUy_r_thr t->diff_2dUy_r_thr
#define daux_diff_dUy_steer t->diff_dUy_steer
#define daux_diff_2dUy_steer_steer t->diff_2dUy_steer_steer
#define daux_diff_dUy_thr t->diff_dUy_thr
#define daux_diff_2dUy_thr_thr t->diff_2dUy_thr_thr
#define daux_diff_dUx_Ux t->diff_dUx_Ux
#define daux_diff_2dUx_Ux_Ux t->diff_2dUx_Ux_Ux
#define daux_diff_2dUx_Ux_Uy t->diff_2dUx_Ux_Uy
#define daux_diff_2dUx_Ux_r t->diff_2dUx_Ux_r
#define daux_diff_2dUx_Ux_steer t->diff_2dUx_Ux_steer
#define daux_diff_2dUx_Ux_thr t->diff_2dUx_Ux_thr
#define daux_diff_dUx_Uy t->diff_dUx_Uy
#define daux_diff_2dUx_Uy_Uy t->diff_2dUx_Uy_Uy
#define daux_diff_2dUx_Uy_r t->diff_2dUx_Uy_r
#define daux_diff_2dUx_Uy_steer t->diff_2dUx_Uy_steer
#define daux_diff_2dUx_Uy_thr t->diff_2dUx_Uy_thr
#define daux_diff_dUx_r t->diff_dUx_r
#define daux_diff_2dUx_r_r t->diff_2dUx_r_r
#define daux_diff_2dUx_r_steer t->diff_2dUx_r_steer
#define daux_diff_2dUx_r_thr t->diff_2dUx_r_thr
#define daux_diff_dUx_steer t->diff_dUx_steer
#define daux_diff_2dUx_steer_steer t->diff_2dUx_steer_steer
#define daux_diff_dUx_thr t->diff_dUx_thr
#define daux_diff_2dUx_thr_thr t->diff_2dUx_thr_thr
#define daux_diff_rev_Ux t->diff_rev_Ux
#define daux_diff_2rev_Ux_Ux t->diff_2rev_Ux_Ux
#define daux_diff_2rev_Ux_thr t->diff_2rev_Ux_thr
#define daux_diff_rev_thr t->diff_rev_thr
#define daux_diff_2rev_thr_thr t->diff_2rev_thr_thr
#define daux_diff_Fxr_Ux t->diff_Fxr_Ux
#define daux_diff_2Fxr_Ux_Ux t->diff_2Fxr_Ux_Ux
#define daux_diff_2Fxr_Ux_Uy t->diff_2Fxr_Ux_Uy
#define daux_diff_2Fxr_Ux_r t->diff_2Fxr_Ux_r
#define daux_diff_2Fxr_Ux_thr t->diff_2Fxr_Ux_thr
#define daux_diff_Fxr_Uy t->diff_Fxr_Uy
#define daux_diff_2Fxr_Uy_Uy t->diff_2Fxr_Uy_Uy
#define daux_diff_2Fxr_Uy_r t->diff_2Fxr_Uy_r
#define daux_diff_2Fxr_Uy_thr t->diff_2Fxr_Uy_thr
#define daux_diff_Fxr_r t->diff_Fxr_r
#define daux_diff_2Fxr_r_r t->diff_2Fxr_r_r
#define daux_diff_2Fxr_r_thr t->diff_2Fxr_r_thr
#define daux_diff_Fxr_thr t->diff_Fxr_thr
#define daux_diff_2Fxr_thr_thr t->diff_2Fxr_thr_thr





static int calcXVariableAux(trajEl_t *t, multipliersEl_t *m, int k, tOptSet *o);
static int calcXUVariableAux(trajEl_t *t, multipliersEl_t *m, int k, tOptSet *o);
static int calcFVariableAux(trajFin_t *t, multipliersFin_t *m, tOptSet *o);
static int calcLAuxDeriv(trajEl_t *t, multipliersEl_t *m, int k, tOptSet *o);
static int calcFAuxDeriv(trajFin_t *t, multipliersFin_t *m, tOptSet *o);
static int bp_derivsL(trajEl_t *t, int k, double **p);
static int bp_derivsF(trajFin_t *t, int k, double **p);

static int ddpL(trajEl_t *t, int k, tOptSet *o) {
    const double *x= t->x;
    const double *u= t->u;
    double **p= o->p;
    
    t->c=(-1.2+sqrt(1.0+pow(x[3],2.0)))*p[8][0]+pow(x[8],2.0)*p[9][0]+pow(x[9],2.0)*p[9][1]+pow(u[0],2.0)*p[12][0]+pow(u[1],2.0)*p[12][1]+p[13][0]*(-p[24][0]+sqrt(pow(p[24][0],2.0)+pow(x[0]-p[25][0],2.0)))+p[13][1]*(-p[24][1]+sqrt(pow(p[24][1],2.0)+pow(x[1]-p[25][1],2.0)))+p[13][2]*(-p[24][2]+sqrt(pow(p[24][2],2.0)+pow(x[2]-p[25][2],2.0)))+p[10][0]*(-p[24][0]+sqrt(pow(p[24][0],2.0)+pow(x[3]-p[25][3],2.0)))+p[10][1]*(-p[24][1]+sqrt(pow(p[24][1],2.0)+pow(x[4]-p[25][4],2.0)))+p[10][2]*(-p[24][2]+sqrt(pow(p[24][2],2.0)+pow(x[5]-p[25][5],2.0)))+p[16][0]*mcond(sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))>p[14][0],0.0,1,pow(1/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))-1/p[14][0],2.0))+p[17][0]*mcond(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4]<0.0,0.0,1,(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4])/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0)));
    if (isNANorINF(t->c))
        {
            PRNT("    @k %d: t->c in line %d is nan or inf: %g\n", k, __LINE__-3,t->c);
            {
                PRNT("        p[10,0]= %g\n",p[10][0]);
                PRNT("        p[24,0]= %g\n",p[24][0]);
                PRNT("        p[25,3]= %g\n",p[25][3]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[10,1]= %g\n",p[10][1]);
                PRNT("        p[24,1]= %g\n",p[24][1]);
                PRNT("        p[25,4]= %g\n",p[25][4]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[8,0]= %g\n",p[8][0]);
                PRNT("        p[9,1]= %g\n",p[9][1]);
                PRNT("        x[9]= %g\n",x[9]);
                PRNT("        p[9,0]= %g\n",p[9][0]);
                PRNT("        x[8]= %g\n",x[8]);
                PRNT("        p[13,2]= %g\n",p[13][2]);
                PRNT("        p[24,2]= %g\n",p[24][2]);
                PRNT("        p[25,2]= %g\n",p[25][2]);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        p[10,2]= %g\n",p[10][2]);
                PRNT("        p[25,5]= %g\n",p[25][5]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[12,1]= %g\n",p[12][1]);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        p[12,0]= %g\n",p[12][0]);
                PRNT("        u[0]= %g\n",u[0]);
                PRNT("        p[13,0]= %g\n",p[13][0]);
                PRNT("        p[25,0]= %g\n",p[25][0]);
                PRNT("        x[0]= %g\n",x[0]);
                PRNT("        p[13,1]= %g\n",p[13][1]);
                PRNT("        p[25,1]= %g\n",p[25][1]);
                PRNT("        x[1]= %g\n",x[1]);
                PRNT("        p[17,0]= %g\n",p[17][0]);
                PRNT("        p[3,0]= %g\n",p[3][0]);
                PRNT("        p[3,1]= %g\n",p[3][1]);
                PRNT("        p[16,0]= %g\n",p[16][0]);
                PRNT("        p[14,0]= %g\n",p[14][0]);
            }
            return(0);
        }
    
    return 1;
}

static int ddpF(trajFin_t *t, tOptSet *o) {
    const double *x= t->x;
    const int k= o->n_hor;
    double **p= o->p;
    
    t->c=p[11][0]*(-p[23][0]+sqrt(pow(p[23][0],2.0)+pow(x[0]-p[25][0],2.0)))+p[11][1]*(-p[23][1]+sqrt(pow(p[23][1],2.0)+pow(x[1]-p[25][1],2.0)))+p[11][2]*(-p[23][2]+sqrt(pow(p[23][2],2.0)+pow(x[2]-p[25][2],2.0)))+p[11][3]*(-p[23][3]+sqrt(pow(p[23][3],2.0)+pow(x[3]-p[25][3],2.0)))+p[11][4]*(-p[23][4]+sqrt(pow(p[23][4],2.0)+pow(x[4]-p[25][4],2.0)))+p[11][5]*(-p[23][5]+sqrt(pow(p[23][5],2.0)+pow(x[5]-p[25][5],2.0)));
    if (isNANorINF(t->c))
        {
            PRNT("    @k %d: t->c in line %d is nan or inf: %g\n", k, __LINE__-3,t->c);
            {
                PRNT("        p[11,3]= %g\n",p[11][3]);
                PRNT("        p[23,3]= %g\n",p[23][3]);
                PRNT("        p[25,3]= %g\n",p[25][3]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[11,4]= %g\n",p[11][4]);
                PRNT("        p[23,4]= %g\n",p[23][4]);
                PRNT("        p[25,4]= %g\n",p[25][4]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[11,2]= %g\n",p[11][2]);
                PRNT("        p[23,2]= %g\n",p[23][2]);
                PRNT("        p[25,2]= %g\n",p[25][2]);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        p[11,5]= %g\n",p[11][5]);
                PRNT("        p[23,5]= %g\n",p[23][5]);
                PRNT("        p[25,5]= %g\n",p[25][5]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[11,0]= %g\n",p[11][0]);
                PRNT("        p[23,0]= %g\n",p[23][0]);
                PRNT("        p[25,0]= %g\n",p[25][0]);
                PRNT("        x[0]= %g\n",x[0]);
                PRNT("        p[11,1]= %g\n",p[11][1]);
                PRNT("        p[23,1]= %g\n",p[23][1]);
                PRNT("        p[25,1]= %g\n",p[25][1]);
                PRNT("        x[1]= %g\n",x[1]);
            }
            return(0);
        }
    
    return 1;
}

static int ddpf(double x_next[], trajEl_t *t, int k, double **p, int N) {
    const double *x= t->x;
    const double *u= t->u;
    
    x_next[0]=x[0]+p[15][0]*aux_dx;
    if (isNANorINF(x_next[0]))
        {
            PRNT("    @k %d: x_next[0] in line %d is nan or inf: %g\n", k, __LINE__-3,x_next[0]);
            {
                PRNT("        aux_dx= %g\n",aux_dx);
                PRNT("        p[15,0]= %g\n",p[15][0]);
                PRNT("        x[0]= %g\n",x[0]);
            }
            return(0);
        }
    x_next[1]=x[1]+p[15][0]*aux_dy;
    if (isNANorINF(x_next[1]))
        {
            PRNT("    @k %d: x_next[1] in line %d is nan or inf: %g\n", k, __LINE__-3,x_next[1]);
            {
                PRNT("        aux_dy= %g\n",aux_dy);
                PRNT("        p[15,0]= %g\n",p[15][0]);
                PRNT("        x[1]= %g\n",x[1]);
            }
            return(0);
        }
    x_next[2]=x[2]+x[5]*p[15][0];
    if (isNANorINF(x_next[2]))
        {
            PRNT("    @k %d: x_next[2] in line %d is nan or inf: %g\n", k, __LINE__-3,x_next[2]);
            {
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        p[15,0]= %g\n",p[15][0]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
    x_next[3]=x[3]+p[15][0]*aux_dUx;
    if (isNANorINF(x_next[3]))
        {
            PRNT("    @k %d: x_next[3] in line %d is nan or inf: %g\n", k, __LINE__-3,x_next[3]);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        aux_dUx= %g\n",aux_dUx);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    x_next[4]=x[4]+p[15][0]*aux_dUy;
    if (isNANorINF(x_next[4]))
        {
            PRNT("    @k %d: x_next[4] in line %d is nan or inf: %g\n", k, __LINE__-3,x_next[4]);
            {
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        aux_dUy= %g\n",aux_dUy);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    x_next[5]=x[5]+p[15][0]*aux_dr;
    if (isNANorINF(x_next[5]))
        {
            PRNT("    @k %d: x_next[5] in line %d is nan or inf: %g\n", k, __LINE__-3,x_next[5]);
            {
                PRNT("        aux_dr= %g\n",aux_dr);
                PRNT("        p[15,0]= %g\n",p[15][0]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
    x_next[6]=u[0];
    if (isNANorINF(x_next[6]))
        {
            PRNT("    @k %d: x_next[6] in line %d is nan or inf: %g\n", k, __LINE__-3,x_next[6]);
            {
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
    x_next[7]=u[1];
    if (isNANorINF(x_next[7]))
        {
            PRNT("    @k %d: x_next[7] in line %d is nan or inf: %g\n", k, __LINE__-3,x_next[7]);
            {
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
    x_next[8]=u[0]-x[6];
    if (isNANorINF(x_next[8]))
        {
            PRNT("    @k %d: x_next[8] in line %d is nan or inf: %g\n", k, __LINE__-3,x_next[8]);
            {
                PRNT("        x[6]= %g\n",x[6]);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
    x_next[9]=u[1]-x[7];
    if (isNANorINF(x_next[9]))
        {
            PRNT("    @k %d: x_next[9] in line %d is nan or inf: %g\n", k, __LINE__-3,x_next[9]);
            {
                PRNT("        x[7]= %g\n",x[7]);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
    return 1;
}

void clampU(double *u, trajEl_t *t, int k, double **p, int N) {
    double limit;
    const double *x= t->x;

// constraint h[1]= limThr[0]-thr
    limit=p[19][0];
    if (u[0]<limit)
        u[0]=limit;

// constraint h[2]= -limThr[1]+thr
    limit=p[19][1];
    if (u[0]>limit)
        u[0]=limit;

// constraint h[3]= limSteer[0]-steer
    limit=p[18][0];
    if (u[1]<limit)
        u[1]=limit;

// constraint h[4]= -limSteer[1]+steer
    limit=p[18][1];
    if (u[1]>limit)
        u[1]=limit;

}

static void limitsU(trajEl_t *t, int k, double **p, int N) {
    int i, j;
    int lower_idx[N_U], upper_idx[N_U], *idx_;
    double limit;
    const double *x= t->x;
    double *hx_, *h_sign;

    for(i= 0; i<N_U; i++) {
        lower_idx[i]= -1;
        upper_idx[i]= -1;
        t->lower[i]= -INF;
        t->upper[i]= INF;
    }        
    
// constraint h[1]= limThr[0]-thr
    limit=p[19][0];
    if (t->lower[0]<limit)
        {
            t->lower[0]=limit;
            lower_idx[0]=0;
        }

// constraint h[2]= -limThr[1]+thr
    limit=p[19][1];
    if (t->upper[0]>limit)
        {
            t->upper[0]=limit;
            upper_idx[0]=1;
        }

// constraint h[3]= limSteer[0]-steer
    limit=p[18][0];
    if (t->lower[1]<limit)
        {
            t->lower[1]=limit;
            lower_idx[1]=2;
        }

// constraint h[4]= -limSteer[1]+steer
    limit=p[18][1];
    if (t->upper[1]>limit)
        {
            t->upper[1]=limit;
            upper_idx[1]=3;
        }


    for(i= 0; i<N_U; i++) {
        t->lower[i]-= t->u[i];
        t->upper[i]-= t->u[i]; 
    }

    for(j= 0; j<2; j++) {
        if(j==0) {
            idx_= lower_idx;
            hx_= t->lower_hx;
            h_sign= t->lower_sign;
        } else {
            idx_= upper_idx;
            hx_= t->upper_hx;
            h_sign= t->upper_sign;
        }
        for(i= 0; i<N_U; i++, hx_+= N_X, h_sign++) {
            switch(idx_[i]) {
                case -1:
                    h_sign[0]= 0.0;
                    break;
                case 0:
// constraint h[1]= limThr[0]-thr
                    hx_[0]=0.0;
                    hx_[1]=0.0;
                    hx_[2]=0.0;
                    hx_[3]=0.0;
                    hx_[4]=0.0;
                    hx_[5]=0.0;
                    hx_[6]=0.0;
                    hx_[7]=0.0;
                    hx_[8]=0.0;
                    hx_[9]=0.0;
                    h_sign[0]=-1.0;
                    break;
                case 1:
// constraint h[2]= -limThr[1]+thr
                    hx_[0]=0.0;
                    hx_[1]=0.0;
                    hx_[2]=0.0;
                    hx_[3]=0.0;
                    hx_[4]=0.0;
                    hx_[5]=0.0;
                    hx_[6]=0.0;
                    hx_[7]=0.0;
                    hx_[8]=0.0;
                    hx_[9]=0.0;
                    h_sign[0]=1.0;
                    break;
                case 2:
// constraint h[3]= limSteer[0]-steer
                    hx_[0]=0.0;
                    hx_[1]=0.0;
                    hx_[2]=0.0;
                    hx_[3]=0.0;
                    hx_[4]=0.0;
                    hx_[5]=0.0;
                    hx_[6]=0.0;
                    hx_[7]=0.0;
                    hx_[8]=0.0;
                    hx_[9]=0.0;
                    h_sign[0]=-1.0;
                    break;
                case 3:
// constraint h[4]= -limSteer[1]+steer
                    hx_[0]=0.0;
                    hx_[1]=0.0;
                    hx_[2]=0.0;
                    hx_[3]=0.0;
                    hx_[4]=0.0;
                    hx_[5]=0.0;
                    hx_[6]=0.0;
                    hx_[7]=0.0;
                    hx_[8]=0.0;
                    hx_[9]=0.0;
                    h_sign[0]=1.0;
                    break;
            }
        }
    }
}

int forward_pass(traj_t *c, tOptSet *o, double alpha, double *csum, int cost_only) {
    int i, k, j;
    double dx;
    double *x0= o->x0;
    int N= o->n_hor;
    double **params= o->p;
    
    trajEl_t *t= o->nominal->t;
    trajFin_t *f= &o->nominal->f;
    trajEl_t *ct= c->t;
    trajFin_t *cf= &c->f;
    
    multipliersEl_t *m= o->multipliers.t;
    multipliersFin_t *mf= &o->multipliers.f;
    
    double *x_next;
    
    csum[0]= 0.0;

    if(!cost_only)
        for(i= 0; i<N_X; i++) ct->x[i]= x0[i]; // ic

    for(k= 0; k<N; k++, t++, ct++, m++) {
        if(!cost_only) {
            if(alpha) {
                for(j= 0; j<N_U; j++)
                    ct->u[j]= t->u[j] + t->l[j]*alpha;
                for(i= 0; i<N_X; i++) {
                    dx= ct->x[i] - t->x[i];
                    
                    for(j= 0; j<N_U; j++) {
                        ct->u[j]+= t->L[MAT_IDX(j, i, N_U)]*dx;
                    }
                }
            } else {
                for(j= 0; j<N_U; j++)
                    ct->u[j]= t->u[j];
            }
        }        
        if(!calcXVariableAux(ct, m, k, o)) return 0;
        
        if(!cost_only)
            clampU(ct->u, ct, k, params, N);
        if(!calcXUVariableAux(ct, m, k, o)) return 0;
        
        if(!cost_only) {
            if(k>=N-1)
                x_next= cf->x;
            else
                x_next= (ct+1)->x;
                
            if(!ddpf(x_next, ct, k, params, N)) return 0;
        }
        
        if(!ddpL(ct, k, o)) return 0;
        csum[0]+= ct->c;        
    }
    
    if(!calcFVariableAux(cf, mf, o)) return 0;
        
    if(!ddpF(cf, o))  return 0;
    csum[0]+= cf->c;
        
    return 1;
}

int calc_derivs(tOptSet *o) {
    int k, i_;
    int N= o->n_hor;

    trajEl_t *t= o->nominal->t + N -1;
    trajFin_t *f= &o->nominal->f;
    
    multipliersEl_t *m= o->multipliers.t + N - 1;
    multipliersFin_t *mf= &o->multipliers.f;

    if(!calcFAuxDeriv(f, mf, o)) return 0;
    if(!bp_derivsF(f, N, o->p)) return 0;
    
#if MULTI_THREADED   
    pthread_mutex_lock(&step_mutex);
    step_calc_done= N;
    pthread_cond_signal(&next_step_condition);
    pthread_mutex_unlock(&step_mutex);
#endif

    for(k= N-1; k>=0; k--, t--, m--) {
        if(!calcLAuxDeriv(t, m, k, o)) return 0;
        if(!bp_derivsL(t, k, o->p)) return 0;
        
        limitsU(t, k, o->p, N);

#if MULTI_THREADED   
        pthread_mutex_lock(&step_mutex);
        step_calc_done= k;
        pthread_cond_signal(&next_step_condition);
        pthread_mutex_unlock(&step_mutex);
#endif
    }
    return 1;
}

static int calcXVariableAux(trajEl_t *t, multipliersEl_t *m, int k, tOptSet *o) {
    const double *x= t->x;
    double **p= o->p;
    const double w_pen= o->w_pen_l;
    
// aux_def[1]
    aux_beta_=mcond(x[3]<0.0&&x[4]>0.0,3.141592653589793-atan(x[4]/(0.001+abs(x[3]))),1,mcond(x[3]<0.0&&x[4]<0.0,-3.141592653589793-atan(x[4]/(0.001+abs(x[3]))),1,atan(x[4]/(0.001+abs(x[3])))));
    if (isNANorINF(aux_beta_))
        {
            PRNT("    @k %d: aux_beta_ in line %d is nan or inf: %g\n", k, __LINE__-3,aux_beta_);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
            }
            return(0);
        }
// aux_def[2]
    aux_dy=aux_U*sin(x[2]+aux_beta_);
    if (isNANorINF(aux_dy))
        {
            PRNT("    @k %d: aux_dy in line %d is nan or inf: %g\n", k, __LINE__-3,aux_dy);
            {
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
            }
            return(0);
        }
// aux_def[3]
    aux_U=sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0));
    if (isNANorINF(aux_U))
        {
            PRNT("    @k %d: aux_U in line %d is nan or inf: %g\n", k, __LINE__-3,aux_U);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
            }
            return(0);
        }
// aux_def[4]
    aux_dx=aux_U*cos(x[2]+aux_beta_);
    if (isNANorINF(aux_dx))
        {
            PRNT("    @k %d: aux_dx in line %d is nan or inf: %g\n", k, __LINE__-3,aux_dx);
            {
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
            }
            return(0);
        }
// aux_def[5]
    aux_alpha_R=atan((x[4]-x[5]*p[5][0])/(0.001+abs(x[3])));
    if (isNANorINF(aux_alpha_R))
        {
            PRNT("    @k %d: aux_alpha_R in line %d is nan or inf: %g\n", k, __LINE__-3,aux_alpha_R);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_def[7]
    aux_F_r=mcond(aux_gamma_R<=3.0*p[1][0]*p[21][0],aux_gamma_R-0.3333333333333333*(2.0-p[22][0]/p[21][0])*pow(aux_gamma_R,2.0)/p[1][0]/p[21][0]+0.1111111111111111*(1.0-0.6666666666666666*p[22][0]/p[21][0])*pow(aux_gamma_R,3.0)/pow(p[1][0],2.0)/pow(p[21][0],2.0),1,p[1][0]*p[22][0]);
    if (isNANorINF(aux_F_r))
        {
            PRNT("    @k %d: aux_F_r in line %d is nan or inf: %g\n", k, __LINE__-3,aux_F_r);
            {
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        p[1,0]= %g\n",p[1][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        p[22,0]= %g\n",p[22][0]);
            }
            return(0);
        }
// aux_def[8]
    aux_Fyr=-p[6][0]*aux_F_r*tan(aux_alpha_R)/aux_gamma_R;
    if (isNANorINF(aux_Fyr))
        {
            PRNT("    @k %d: aux_Fyr in line %d is nan or inf: %g\n", k, __LINE__-3,aux_Fyr);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
            }
            return(0);
        }
// aux_def[15]
    aux_gamma_F=abs(p[6][0])*abs(tan(aux_alpha_R));
    if (isNANorINF(aux_gamma_F))
        {
            PRNT("    @k %d: aux_gamma_F in line %d is nan or inf: %g\n", k, __LINE__-3,aux_gamma_F);
            {
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
            }
            return(0);
        }
// aux_def[16]
    aux_F_f=mcond(aux_gamma_F<=3.0*p[0][0]*p[21][0],aux_gamma_F-0.3333333333333333*(2.0-p[22][0]/p[21][0])*pow(aux_gamma_F,2.0)/p[0][0]/p[21][0]+0.1111111111111111*(1.0-0.6666666666666666*p[22][0]/p[21][0])*pow(aux_gamma_F,3.0)/pow(p[0][0],2.0)/pow(p[21][0],2.0),1,p[0][0]*p[22][0]);
    if (isNANorINF(aux_F_f))
        {
            PRNT("    @k %d: aux_F_f in line %d is nan or inf: %g\n", k, __LINE__-3,aux_F_f);
            {
                PRNT("        aux_gamma_F= %g\n",aux_gamma_F);
                PRNT("        p[0,0]= %g\n",p[0][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        p[22,0]= %g\n",p[22][0]);
            }
            return(0);
        }
    return 1;
} 

static int calcXUVariableAux(trajEl_t *t, multipliersEl_t *m, int k, tOptSet *o) {
    const double *x= t->x;
    const double *u= t->u;
    double **p= o->p;
    const double w_pen= o->w_pen_l;
    
// aux_def[6]
    aux_gamma_R=sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0));
    if (isNANorINF(aux_gamma_R))
        {
            PRNT("    @k %d: aux_gamma_R in line %d is nan or inf: %g\n", k, __LINE__-3,aux_gamma_R);
            {
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[7,0]= %g\n",p[7][0]);
            }
            return(0);
        }
// aux_def[9]
    aux_dr=(p[4][0]*cos(u[1])*aux_Fyf-p[5][0]*aux_Fyr)/p[2][0];
    if (isNANorINF(aux_dr))
        {
            PRNT("    @k %d: aux_dr in line %d is nan or inf: %g\n", k, __LINE__-3,aux_dr);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        aux_Fyr= %g\n",aux_Fyr);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_def[10]
    aux_K_=(u[0]-x[3])/(0.001+abs(x[3]));
    if (isNANorINF(aux_K_))
        {
            PRNT("    @k %d: aux_K_ in line %d is nan or inf: %g\n", k, __LINE__-3,aux_K_);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
// aux_def[11]
    aux_K=abs(aux_K_);
    if (isNANorINF(aux_K))
        {
            PRNT("    @k %d: aux_K in line %d is nan or inf: %g\n", k, __LINE__-3,aux_K);
            {
                PRNT("        aux_K_= %g\n",aux_K_);
            }
            return(0);
        }
// aux_def[12]
    aux_alpha=-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])));
    if (isNANorINF(aux_alpha))
        {
            PRNT("    @k %d: aux_alpha in line %d is nan or inf: %g\n", k, __LINE__-3,aux_alpha);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_def[13]
    aux_alpha_F=mcond(abs(aux_alpha)>1.5707963267948966,aux_alpha*(3.141592653589793-abs(aux_alpha))/abs(aux_alpha),1,aux_alpha);
    if (isNANorINF(aux_alpha_F))
        {
            PRNT("    @k %d: aux_alpha_F in line %d is nan or inf: %g\n", k, __LINE__-3,aux_alpha_F);
            {
                PRNT("        aux_alpha= %g\n",aux_alpha);
            }
            return(0);
        }
// aux_def[14]
    aux_Fyf=-p[6][0]*aux_F_f*tan(aux_alpha_F)/aux_gamma_F;
    if (isNANorINF(aux_Fyf))
        {
            PRNT("    @k %d: aux_Fyf in line %d is nan or inf: %g\n", k, __LINE__-3,aux_Fyf);
            {
                PRNT("        aux_F_f= %g\n",aux_F_f);
                PRNT("        aux_alpha_F= %g\n",aux_alpha_F);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        aux_gamma_F= %g\n",aux_gamma_F);
            }
            return(0);
        }
// aux_def[17]
    aux_dUy=-x[3]*x[5]+(cos(u[1])*aux_Fyf+aux_Fyr)/p[20][0];
    if (isNANorINF(aux_dUy))
        {
            PRNT("    @k %d: aux_dUy in line %d is nan or inf: %g\n", k, __LINE__-3,aux_dUy);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        aux_Fyr= %g\n",aux_Fyr);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_def[18]
    aux_dUx=x[4]*x[5]+(aux_Fxr-sin(u[1])*aux_Fyf)/p[20][0];
    if (isNANorINF(aux_dUx))
        {
            PRNT("    @k %d: aux_dUx in line %d is nan or inf: %g\n", k, __LINE__-3,aux_dUx);
            {
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        aux_Fxr= %g\n",aux_Fxr);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_def[19]
    aux_rev=mcond((u[0]-x[3])/(0.001+abs(x[3]))<0.0,-1.0,1,1.0);
    if (isNANorINF(aux_rev))
        {
            PRNT("    @k %d: aux_rev in line %d is nan or inf: %g\n", k, __LINE__-3,aux_rev);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
// aux_def[20]
    aux_Fxr=p[7][0]*aux_F_r*aux_K*aux_rev/(1.0+aux_K)/aux_gamma_R;
    if (isNANorINF(aux_Fxr))
        {
            PRNT("    @k %d: aux_Fxr in line %d is nan or inf: %g\n", k, __LINE__-3,aux_Fxr);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        aux_rev= %g\n",aux_rev);
            }
            return(0);
        }
    return 1;
} 

static int calcFVariableAux(trajFin_t *t, multipliersFin_t *m, tOptSet *o) {
    const double *x= t->x;
    double **p= o->p;
    const double w_pen= o->w_pen_f;
    const int k= o->n_hor;
    
    return 1;
} 

static int calcLAuxDeriv(trajEl_t *t, multipliersEl_t *m, int k, tOptSet *o) {
    const double *x= t->x;
    const double *u= t->u;
    const double w_pen= o->w_pen_l;
    double **p= o->p;
    
// aux_deriv[1]
    daux_diff_beta__Ux=mcond(x[3]<0.0&&x[4]>0.0,x[4]*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0)),1,mcond(x[3]<0.0&&x[4]<0.0,x[4]*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0)),1,-x[4]*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0))));
    if (isNANorINF(daux_diff_beta__Ux))
        {
            PRNT("    @k %d: daux_diff_beta__Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_beta__Ux);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
            }
            return(0);
        }
// aux_deriv[2]
    daux_diff_2beta__Ux_Ux=mcond(x[3]<0.0&&x[4]>0.0,x[4]*mcond(x[3]>0.0,0.0,1,0.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0))+2.0*pow(x[4],3.0)*pow(mcond(x[3]>0.0,1.0,1,-1.0),2.0)/pow(0.001+abs(x[3]),5.0)/pow(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0),2.0)-2.0*x[4]*pow(mcond(x[3]>0.0,1.0,1,-1.0),2.0)/pow(0.001+abs(x[3]),3.0)/(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0)),1,mcond(x[3]<0.0&&x[4]<0.0,x[4]*mcond(x[3]>0.0,0.0,1,0.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0))+2.0*pow(x[4],3.0)*pow(mcond(x[3]>0.0,1.0,1,-1.0),2.0)/pow(0.001+abs(x[3]),5.0)/pow(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0),2.0)-2.0*x[4]*pow(mcond(x[3]>0.0,1.0,1,-1.0),2.0)/pow(0.001+abs(x[3]),3.0)/(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0)),1,-x[4]*mcond(x[3]>0.0,0.0,1,0.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0))-2.0*pow(x[4],3.0)*pow(mcond(x[3]>0.0,1.0,1,-1.0),2.0)/pow(0.001+abs(x[3]),5.0)/pow(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),
     2.0),2.0)+2.0*x[4]*pow(mcond(x[3]>0.0,1.0,1,-1.0),2.0)/pow(0.001+abs(x[3]),3.0)/(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0))));
    if (isNANorINF(daux_diff_2beta__Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2beta__Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2beta__Ux_Ux);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
            }
            return(0);
        }
// aux_deriv[3]
    daux_diff_2beta__Ux_Uy=mcond(x[3]<0.0&&x[4]>0.0,-2.0*pow(x[4],2.0)*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),4.0)/pow(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0),2.0)+mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0)),1,mcond(x[3]<0.0&&x[4]<0.0,-2.0*pow(x[4],2.0)*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),4.0)/pow(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0),2.0)+mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0)),1,2.0*pow(x[4],2.0)*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),4.0)/pow(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0),2.0)-mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0))));
    if (isNANorINF(daux_diff_2beta__Ux_Uy))
        {
            PRNT("    @k %d: daux_diff_2beta__Ux_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2beta__Ux_Uy);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
            }
            return(0);
        }
// aux_deriv[4]
    daux_diff_beta__Uy=mcond(x[3]<0.0&&x[4]>0.0,-1/(0.001+abs(x[3]))/(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0)),1,mcond(x[3]<0.0&&x[4]<0.0,-1/(0.001+abs(x[3]))/(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0)),1,1/(0.001+abs(x[3]))/(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0))));
    if (isNANorINF(daux_diff_beta__Uy))
        {
            PRNT("    @k %d: daux_diff_beta__Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_beta__Uy);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
            }
            return(0);
        }
// aux_deriv[5]
    daux_diff_2beta__Uy_Uy=mcond(x[3]<0.0&&x[4]>0.0,2.0*x[4]/pow(0.001+abs(x[3]),3.0)/pow(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0),2.0),1,mcond(x[3]<0.0&&x[4]<0.0,2.0*x[4]/pow(0.001+abs(x[3]),3.0)/pow(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0),2.0),1,-2.0*x[4]/pow(0.001+abs(x[3]),3.0)/pow(1.0+pow(x[4],2.0)/pow(0.001+abs(x[3]),2.0),2.0)));
    if (isNANorINF(daux_diff_2beta__Uy_Uy))
        {
            PRNT("    @k %d: daux_diff_2beta__Uy_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2beta__Uy_Uy);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
            }
            return(0);
        }
// aux_deriv[6]
    daux_diff_dy_Ux=derivative(aux_U,x[3],1.0)*sin(x[2]+aux_beta_)+aux_U*cos(x[2]+aux_beta_)*daux_diff_beta__Ux;
    if (isNANorINF(daux_diff_dy_Ux))
        {
            PRNT("    @k %d: daux_diff_dy_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dy_Ux);
            {
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        daux_diff_beta__Ux= %g\n",daux_diff_beta__Ux);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        x[3]= %g\n",x[3]);
            }
            return(0);
        }
// aux_deriv[7]
    daux_diff_2dy_Ux_Ux=derivative(aux_U,x[3],2.0)*sin(x[2]+aux_beta_)+aux_U*cos(x[2]+aux_beta_)*daux_diff_2beta__Ux_Ux+2.0*derivative(aux_U,x[3],1.0)*cos(x[2]+aux_beta_)*daux_diff_beta__Ux-aux_U*sin(x[2]+aux_beta_)*pow(daux_diff_beta__Ux,2.0);
    if (isNANorINF(daux_diff_2dy_Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2dy_Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dy_Ux_Ux);
            {
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        daux_diff_2beta__Ux_Ux= %g\n",daux_diff_2beta__Ux_Ux);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        daux_diff_beta__Ux= %g\n",daux_diff_beta__Ux);
            }
            return(0);
        }
// aux_deriv[8]
    daux_diff_2dy_Ux_Uy=derivative(aux_U,x[3],1.0,x[4],1.0)*sin(x[2]+aux_beta_)+aux_U*cos(x[2]+aux_beta_)*daux_diff_2beta__Ux_Uy+derivative(aux_U,x[4],1.0)*cos(x[2]+aux_beta_)*daux_diff_beta__Ux+derivative(aux_U,x[3],1.0)*cos(x[2]+aux_beta_)*daux_diff_beta__Uy-aux_U*sin(x[2]+aux_beta_)*daux_diff_beta__Ux*daux_diff_beta__Uy;
    if (isNANorINF(daux_diff_2dy_Ux_Uy))
        {
            PRNT("    @k %d: daux_diff_2dy_Ux_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dy_Ux_Uy);
            {
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        daux_diff_2beta__Ux_Uy= %g\n",daux_diff_2beta__Ux_Uy);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        daux_diff_beta__Ux= %g\n",daux_diff_beta__Ux);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        daux_diff_beta__Uy= %g\n",daux_diff_beta__Uy);
            }
            return(0);
        }
// aux_deriv[9]
    daux_diff_2dy_Ux_phi=derivative(aux_U,x[3],1.0)*cos(x[2]+aux_beta_)-aux_U*sin(x[2]+aux_beta_)*daux_diff_beta__Ux;
    if (isNANorINF(daux_diff_2dy_Ux_phi))
        {
            PRNT("    @k %d: daux_diff_2dy_Ux_phi in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dy_Ux_phi);
            {
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        daux_diff_beta__Ux= %g\n",daux_diff_beta__Ux);
            }
            return(0);
        }
// aux_deriv[10]
    daux_diff_dy_Uy=derivative(aux_U,x[4],1.0)*sin(x[2]+aux_beta_)+aux_U*cos(x[2]+aux_beta_)*daux_diff_beta__Uy;
    if (isNANorINF(daux_diff_dy_Uy))
        {
            PRNT("    @k %d: daux_diff_dy_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dy_Uy);
            {
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        daux_diff_beta__Uy= %g\n",daux_diff_beta__Uy);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        x[4]= %g\n",x[4]);
            }
            return(0);
        }
// aux_deriv[11]
    daux_diff_2dy_Uy_Uy=derivative(aux_U,x[4],2.0)*sin(x[2]+aux_beta_)+aux_U*cos(x[2]+aux_beta_)*daux_diff_2beta__Uy_Uy+2.0*derivative(aux_U,x[4],1.0)*cos(x[2]+aux_beta_)*daux_diff_beta__Uy-aux_U*sin(x[2]+aux_beta_)*pow(daux_diff_beta__Uy,2.0);
    if (isNANorINF(daux_diff_2dy_Uy_Uy))
        {
            PRNT("    @k %d: daux_diff_2dy_Uy_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dy_Uy_Uy);
            {
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        daux_diff_2beta__Uy_Uy= %g\n",daux_diff_2beta__Uy_Uy);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        daux_diff_beta__Uy= %g\n",daux_diff_beta__Uy);
            }
            return(0);
        }
// aux_deriv[12]
    daux_diff_2dy_Uy_phi=derivative(aux_U,x[4],1.0)*cos(x[2]+aux_beta_)-aux_U*sin(x[2]+aux_beta_)*daux_diff_beta__Uy;
    if (isNANorINF(daux_diff_2dy_Uy_phi))
        {
            PRNT("    @k %d: daux_diff_2dy_Uy_phi in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dy_Uy_phi);
            {
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        daux_diff_beta__Uy= %g\n",daux_diff_beta__Uy);
            }
            return(0);
        }
// aux_deriv[13]
    daux_diff_dy_phi=aux_U*cos(x[2]+aux_beta_);
    if (isNANorINF(daux_diff_dy_phi))
        {
            PRNT("    @k %d: daux_diff_dy_phi in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dy_phi);
            {
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
            }
            return(0);
        }
// aux_deriv[14]
    daux_diff_2dy_phi_phi=-aux_U*sin(x[2]+aux_beta_);
    if (isNANorINF(daux_diff_2dy_phi_phi))
        {
            PRNT("    @k %d: daux_diff_2dy_phi_phi in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dy_phi_phi);
            {
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
            }
            return(0);
        }
// aux_deriv[15]
    daux_diff_U_Ux=x[3]/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0));
    if (isNANorINF(daux_diff_U_Ux))
        {
            PRNT("    @k %d: daux_diff_U_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_U_Ux);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
            }
            return(0);
        }
// aux_deriv[16]
    daux_diff_2U_Ux_Ux=-pow(x[3],2.0)/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),1.5)+1/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0));
    if (isNANorINF(daux_diff_2U_Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2U_Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2U_Ux_Ux);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
            }
            return(0);
        }
// aux_deriv[17]
    daux_diff_2U_Ux_Uy=-x[3]*x[4]/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),1.5);
    if (isNANorINF(daux_diff_2U_Ux_Uy))
        {
            PRNT("    @k %d: daux_diff_2U_Ux_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2U_Ux_Uy);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
            }
            return(0);
        }
// aux_deriv[18]
    daux_diff_U_Uy=x[4]/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0));
    if (isNANorINF(daux_diff_U_Uy))
        {
            PRNT("    @k %d: daux_diff_U_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_U_Uy);
            {
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        x[3]= %g\n",x[3]);
            }
            return(0);
        }
// aux_deriv[19]
    daux_diff_2U_Uy_Uy=-pow(x[4],2.0)/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),1.5)+1/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0));
    if (isNANorINF(daux_diff_2U_Uy_Uy))
        {
            PRNT("    @k %d: daux_diff_2U_Uy_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2U_Uy_Uy);
            {
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        x[3]= %g\n",x[3]);
            }
            return(0);
        }
// aux_deriv[20]
    daux_diff_dx_Ux=cos(x[2]+aux_beta_)*daux_diff_U_Ux-aux_U*sin(x[2]+aux_beta_)*daux_diff_beta__Ux;
    if (isNANorINF(daux_diff_dx_Ux))
        {
            PRNT("    @k %d: daux_diff_dx_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dx_Ux);
            {
                PRNT("        daux_diff_U_Ux= %g\n",daux_diff_U_Ux);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        daux_diff_beta__Ux= %g\n",daux_diff_beta__Ux);
            }
            return(0);
        }
// aux_deriv[21]
    daux_diff_2dx_Ux_Ux=cos(x[2]+aux_beta_)*daux_diff_2U_Ux_Ux-aux_U*sin(x[2]+aux_beta_)*daux_diff_2beta__Ux_Ux-2.0*sin(x[2]+aux_beta_)*daux_diff_U_Ux*daux_diff_beta__Ux-aux_U*cos(x[2]+aux_beta_)*pow(daux_diff_beta__Ux,2.0);
    if (isNANorINF(daux_diff_2dx_Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2dx_Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dx_Ux_Ux);
            {
                PRNT("        daux_diff_2U_Ux_Ux= %g\n",daux_diff_2U_Ux_Ux);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        daux_diff_beta__Ux= %g\n",daux_diff_beta__Ux);
                PRNT("        daux_diff_2beta__Ux_Ux= %g\n",daux_diff_2beta__Ux_Ux);
                PRNT("        daux_diff_U_Ux= %g\n",daux_diff_U_Ux);
            }
            return(0);
        }
// aux_deriv[22]
    daux_diff_2dx_Ux_Uy=cos(x[2]+aux_beta_)*daux_diff_2U_Ux_Uy-aux_U*sin(x[2]+aux_beta_)*daux_diff_2beta__Ux_Uy-sin(x[2]+aux_beta_)*daux_diff_U_Uy*daux_diff_beta__Ux-sin(x[2]+aux_beta_)*daux_diff_U_Ux*daux_diff_beta__Uy-aux_U*cos(x[2]+aux_beta_)*daux_diff_beta__Ux*daux_diff_beta__Uy;
    if (isNANorINF(daux_diff_2dx_Ux_Uy))
        {
            PRNT("    @k %d: daux_diff_2dx_Ux_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dx_Ux_Uy);
            {
                PRNT("        daux_diff_2U_Ux_Uy= %g\n",daux_diff_2U_Ux_Uy);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        daux_diff_beta__Ux= %g\n",daux_diff_beta__Ux);
                PRNT("        daux_diff_beta__Uy= %g\n",daux_diff_beta__Uy);
                PRNT("        daux_diff_2beta__Ux_Uy= %g\n",daux_diff_2beta__Ux_Uy);
                PRNT("        daux_diff_U_Uy= %g\n",daux_diff_U_Uy);
                PRNT("        daux_diff_U_Ux= %g\n",daux_diff_U_Ux);
            }
            return(0);
        }
// aux_deriv[23]
    daux_diff_2dx_Ux_phi=-sin(x[2]+aux_beta_)*daux_diff_U_Ux-aux_U*cos(x[2]+aux_beta_)*daux_diff_beta__Ux;
    if (isNANorINF(daux_diff_2dx_Ux_phi))
        {
            PRNT("    @k %d: daux_diff_2dx_Ux_phi in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dx_Ux_phi);
            {
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        daux_diff_beta__Ux= %g\n",daux_diff_beta__Ux);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        daux_diff_U_Ux= %g\n",daux_diff_U_Ux);
            }
            return(0);
        }
// aux_deriv[24]
    daux_diff_dx_Uy=cos(x[2]+aux_beta_)*daux_diff_U_Uy-aux_U*sin(x[2]+aux_beta_)*daux_diff_beta__Uy;
    if (isNANorINF(daux_diff_dx_Uy))
        {
            PRNT("    @k %d: daux_diff_dx_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dx_Uy);
            {
                PRNT("        daux_diff_U_Uy= %g\n",daux_diff_U_Uy);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        daux_diff_beta__Uy= %g\n",daux_diff_beta__Uy);
            }
            return(0);
        }
// aux_deriv[25]
    daux_diff_2dx_Uy_Uy=cos(x[2]+aux_beta_)*daux_diff_2U_Uy_Uy-aux_U*sin(x[2]+aux_beta_)*daux_diff_2beta__Uy_Uy-2.0*sin(x[2]+aux_beta_)*daux_diff_U_Uy*daux_diff_beta__Uy-aux_U*cos(x[2]+aux_beta_)*pow(daux_diff_beta__Uy,2.0);
    if (isNANorINF(daux_diff_2dx_Uy_Uy))
        {
            PRNT("    @k %d: daux_diff_2dx_Uy_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dx_Uy_Uy);
            {
                PRNT("        daux_diff_2U_Uy_Uy= %g\n",daux_diff_2U_Uy_Uy);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        daux_diff_beta__Uy= %g\n",daux_diff_beta__Uy);
                PRNT("        daux_diff_2beta__Uy_Uy= %g\n",daux_diff_2beta__Uy_Uy);
                PRNT("        daux_diff_U_Uy= %g\n",daux_diff_U_Uy);
            }
            return(0);
        }
// aux_deriv[26]
    daux_diff_2dx_Uy_phi=-sin(x[2]+aux_beta_)*daux_diff_U_Uy-aux_U*cos(x[2]+aux_beta_)*daux_diff_beta__Uy;
    if (isNANorINF(daux_diff_2dx_Uy_phi))
        {
            PRNT("    @k %d: daux_diff_2dx_Uy_phi in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dx_Uy_phi);
            {
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        daux_diff_beta__Uy= %g\n",daux_diff_beta__Uy);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        daux_diff_U_Uy= %g\n",daux_diff_U_Uy);
            }
            return(0);
        }
// aux_deriv[27]
    daux_diff_dx_phi=-aux_U*sin(x[2]+aux_beta_);
    if (isNANorINF(daux_diff_dx_phi))
        {
            PRNT("    @k %d: daux_diff_dx_phi in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dx_phi);
            {
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
            }
            return(0);
        }
// aux_deriv[28]
    daux_diff_2dx_phi_phi=-aux_U*cos(x[2]+aux_beta_);
    if (isNANorINF(daux_diff_2dx_phi_phi))
        {
            PRNT("    @k %d: daux_diff_2dx_phi_phi in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dx_phi_phi);
            {
                PRNT("        aux_U= %g\n",aux_U);
                PRNT("        aux_beta_= %g\n",aux_beta_);
                PRNT("        x[2]= %g\n",x[2]);
            }
            return(0);
        }
// aux_deriv[29]
    daux_diff_alpha_R_Ux=-(x[4]-x[5]*p[5][0])*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0));
    if (isNANorINF(daux_diff_alpha_R_Ux))
        {
            PRNT("    @k %d: daux_diff_alpha_R_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_alpha_R_Ux);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[30]
    daux_diff_2alpha_R_Ux_Ux=-(x[4]-x[5]*p[5][0])*mcond(x[3]>0.0,0.0,1,0.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))-2.0*pow(x[4]-x[5]*p[5][0],3.0)*pow(mcond(x[3]>0.0,1.0,1,-1.0),2.0)/pow(0.001+abs(x[3]),5.0)/pow(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0),2.0)+2.0*(x[4]-x[5]*p[5][0])*pow(mcond(x[3]>0.0,1.0,1,-1.0),2.0)/pow(0.001+abs(x[3]),3.0)/(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0));
    if (isNANorINF(daux_diff_2alpha_R_Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2alpha_R_Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_R_Ux_Ux);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[31]
    daux_diff_2alpha_R_Ux_Uy=2.0*pow(x[4]-x[5]*p[5][0],2.0)*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),4.0)/pow(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0),2.0)-mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0));
    if (isNANorINF(daux_diff_2alpha_R_Ux_Uy))
        {
            PRNT("    @k %d: daux_diff_2alpha_R_Ux_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_R_Ux_Uy);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[32]
    daux_diff_2alpha_R_Ux_r=-2.0*p[5][0]*pow(x[4]-x[5]*p[5][0],2.0)*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),4.0)/pow(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0),2.0)+p[5][0]*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0));
    if (isNANorINF(daux_diff_2alpha_R_Ux_r))
        {
            PRNT("    @k %d: daux_diff_2alpha_R_Ux_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_R_Ux_r);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[33]
    daux_diff_alpha_R_Uy=1/(0.001+abs(x[3]))/(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0));
    if (isNANorINF(daux_diff_alpha_R_Uy))
        {
            PRNT("    @k %d: daux_diff_alpha_R_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_alpha_R_Uy);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[34]
    daux_diff_2alpha_R_Uy_Uy=-2.0*(x[4]-x[5]*p[5][0])/pow(0.001+abs(x[3]),3.0)/pow(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0),2.0);
    if (isNANorINF(daux_diff_2alpha_R_Uy_Uy))
        {
            PRNT("    @k %d: daux_diff_2alpha_R_Uy_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_R_Uy_Uy);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[35]
    daux_diff_2alpha_R_Uy_r=2.0*p[5][0]*(x[4]-x[5]*p[5][0])/pow(0.001+abs(x[3]),3.0)/pow(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0),2.0);
    if (isNANorINF(daux_diff_2alpha_R_Uy_r))
        {
            PRNT("    @k %d: daux_diff_2alpha_R_Uy_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_R_Uy_r);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[36]
    daux_diff_alpha_R_r=-p[5][0]/(0.001+abs(x[3]))/(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0));
    if (isNANorINF(daux_diff_alpha_R_r))
        {
            PRNT("    @k %d: daux_diff_alpha_R_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_alpha_R_r);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[37]
    daux_diff_2alpha_R_r_r=-2.0*pow(p[5][0],2.0)*(x[4]-x[5]*p[5][0])/pow(0.001+abs(x[3]),3.0)/pow(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0),2.0);
    if (isNANorINF(daux_diff_2alpha_R_r_r))
        {
            PRNT("    @k %d: daux_diff_2alpha_R_r_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_R_r_r);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[38]
    daux_diff_gamma_R_Ux=0.5*(-2.0*pow(p[7][0],2.0)*pow(aux_K,2.0)*derivative(aux_K,x[3],1.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[7][0],2.0)*aux_K*derivative(aux_K,x[3],1.0)/pow(1.0+aux_K,2.0)-2.0*pow(p[6][0],2.0)*derivative(aux_K,x[3],1.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_alpha_R_Ux/pow(1.0+aux_K,2.0))/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0));
    if (isNANorINF(daux_diff_gamma_R_Ux))
        {
            PRNT("    @k %d: daux_diff_gamma_R_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_gamma_R_Ux);
            {
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        daux_diff_alpha_R_Ux= %g\n",daux_diff_alpha_R_Ux);
            }
            return(0);
        }
// aux_deriv[39]
    daux_diff_2gamma_R_Ux_Ux=0.5*(6.0*pow(p[7][0],2.0)*pow(aux_K,2.0)*pow(derivative(aux_K,x[3],1.0),2.0)/pow(1.0+aux_K,4.0)-8.0*pow(p[7][0],2.0)*aux_K*pow(derivative(aux_K,x[3],1.0),2.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[7][0],2.0)*pow(derivative(aux_K,x[3],1.0),2.0)/pow(1.0+aux_K,2.0)-2.0*pow(p[7][0],2.0)*pow(aux_K,2.0)*derivative(aux_K,x[3],2.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[7][0],2.0)*aux_K*derivative(aux_K,x[3],2.0)/pow(1.0+aux_K,2.0)+6.0*pow(p[6][0],2.0)*pow(derivative(aux_K,x[3],1.0),2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,4.0)-2.0*pow(p[6][0],2.0)*derivative(aux_K,x[3],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_2alpha_R_Ux_Ux/pow(1.0+aux_K,2.0)-8.0*pow(p[6][0],2.0)*derivative(aux_K,x[3],1.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_alpha_R_Ux/pow(1.0+aux_K,3.0)+2.0*pow(p[6][0],2.0)*pow(sec(aux_alpha_R),4.0)*pow(daux_diff_alpha_R_Ux,2.0)/pow(1.0+aux_K,2.0)+4.0*pow(p[6][0],2.0)*pow(sec(
     aux_alpha_R),2.0)*pow(tan(aux_alpha_R),2.0)*pow(daux_diff_alpha_R_Ux,2.0)/pow(1.0+aux_K,2.0))/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))-0.25*pow(-2.0*pow(p[7][0],2.0)*pow(aux_K,2.0)*derivative(aux_K,x[3],1.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[7][0],2.0)*aux_K*derivative(aux_K,x[3],1.0)/pow(1.0+aux_K,2.0)-2.0*pow(p[6][0],2.0)*derivative(aux_K,x[3],1.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_alpha_R_Ux/pow(1.0+aux_K,2.0),2.0)/pow(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0),1.5);
    if (isNANorINF(daux_diff_2gamma_R_Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2gamma_R_Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2gamma_R_Ux_Ux);
            {
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        daux_diff_2alpha_R_Ux_Ux= %g\n",daux_diff_2alpha_R_Ux_Ux);
                PRNT("        daux_diff_alpha_R_Ux= %g\n",daux_diff_alpha_R_Ux);
            }
            return(0);
        }
// aux_deriv[40]
    daux_diff_2gamma_R_Ux_Uy=-0.5*pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*(-2.0*pow(p[7][0],2.0)*pow(aux_K,2.0)*derivative(aux_K,x[3],1.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[7][0],2.0)*aux_K*derivative(aux_K,x[3],1.0)/pow(1.0+aux_K,2.0)-2.0*pow(p[6][0],2.0)*derivative(aux_K,x[3],1.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_alpha_R_Ux/pow(1.0+aux_K,2.0))*daux_diff_alpha_R_Uy/pow(1.0+aux_K,2.0)/pow(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0),1.5)+0.5*(2.0*pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_2alpha_R_Ux_Uy/pow(1.0+aux_K,2.0)-4.0*pow(p[6][0],2.0)*derivative(aux_K,x[3],1.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_alpha_R_Uy/pow(1.0+aux_K,3.0)+2.0*pow(p[6][0],2.0)*pow(sec(aux_alpha_R),4.0)*daux_diff_alpha_R_Ux*daux_diff_alpha_R_Uy/pow(1.0+aux_K,2.0)+4.0*pow(p[6][0],2.0)*pow(sec(
     aux_alpha_R),2.0)*pow(tan(aux_alpha_R),2.0)*daux_diff_alpha_R_Ux*daux_diff_alpha_R_Uy/pow(1.0+aux_K,2.0))/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0));
    if (isNANorINF(daux_diff_2gamma_R_Ux_Uy))
        {
            PRNT("    @k %d: daux_diff_2gamma_R_Ux_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2gamma_R_Ux_Uy);
            {
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        daux_diff_alpha_R_Ux= %g\n",daux_diff_alpha_R_Ux);
                PRNT("        daux_diff_alpha_R_Uy= %g\n",daux_diff_alpha_R_Uy);
                PRNT("        daux_diff_2alpha_R_Ux_Uy= %g\n",daux_diff_2alpha_R_Ux_Uy);
            }
            return(0);
        }
// aux_deriv[41]
    daux_diff_2gamma_R_Ux_r=-0.5*pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*(-2.0*pow(p[7][0],2.0)*pow(aux_K,2.0)*derivative(aux_K,x[3],1.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[7][0],2.0)*aux_K*derivative(aux_K,x[3],1.0)/pow(1.0+aux_K,2.0)-2.0*pow(p[6][0],2.0)*derivative(aux_K,x[3],1.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_alpha_R_Ux/pow(1.0+aux_K,2.0))*daux_diff_alpha_R_r/pow(1.0+aux_K,2.0)/pow(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0),1.5)+0.5*(2.0*pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_2alpha_R_Ux_r/pow(1.0+aux_K,2.0)-4.0*pow(p[6][0],2.0)*derivative(aux_K,x[3],1.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_alpha_R_r/pow(1.0+aux_K,3.0)+2.0*pow(p[6][0],2.0)*pow(sec(aux_alpha_R),4.0)*daux_diff_alpha_R_Ux*daux_diff_alpha_R_r/pow(1.0+aux_K,2.0)+4.0*pow(p[6][0],2.0)*pow(sec(aux_alpha_R),
     2.0)*pow(tan(aux_alpha_R),2.0)*daux_diff_alpha_R_Ux*daux_diff_alpha_R_r/pow(1.0+aux_K,2.0))/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0));
    if (isNANorINF(daux_diff_2gamma_R_Ux_r))
        {
            PRNT("    @k %d: daux_diff_2gamma_R_Ux_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2gamma_R_Ux_r);
            {
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        daux_diff_alpha_R_Ux= %g\n",daux_diff_alpha_R_Ux);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
                PRNT("        daux_diff_2alpha_R_Ux_r= %g\n",daux_diff_2alpha_R_Ux_r);
            }
            return(0);
        }
// aux_deriv[42]
    daux_diff_2gamma_R_Ux_thr=-0.25*(-2.0*pow(p[7][0],2.0)*pow(aux_K,2.0)*derivative(aux_K,u[0],1.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[7][0],2.0)*aux_K*derivative(aux_K,u[0],1.0)/pow(1.0+aux_K,2.0)-2.0*pow(p[6][0],2.0)*derivative(aux_K,u[0],1.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,3.0))*(-2.0*pow(p[7][0],2.0)*pow(aux_K,2.0)*derivative(aux_K,x[3],1.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[7][0],2.0)*aux_K*derivative(aux_K,x[3],1.0)/pow(1.0+aux_K,2.0)-2.0*pow(p[6][0],2.0)*derivative(aux_K,x[3],1.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_alpha_R_Ux/pow(1.0+aux_K,2.0))/pow(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0),1.5)+0.5*(-2.0*pow(p[7][0],2.0)*pow(aux_K,2.0)*derivative(aux_K,u[0],1.0,x[3],1.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[7][0],2.0)*aux_K*derivative(aux_K,u[0],1.0,x[3],1.0)/pow(1.0+aux_K,2.0)+6.0*pow(p[7][0],2.0)*pow(aux_K,2.0)*derivative(aux_K,u[0]
     ,1.0)*derivative(aux_K,x[3],1.0)/pow(1.0+aux_K,4.0)-8.0*pow(p[7][0],2.0)*aux_K*derivative(aux_K,u[0],1.0)*derivative(aux_K,x[3],1.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[7][0],2.0)*derivative(aux_K,u[0],1.0)*derivative(aux_K,x[3],1.0)/pow(1.0+aux_K,2.0)-2.0*pow(p[6][0],2.0)*derivative(aux_K,u[0],1.0,x[3],1.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,3.0)+6.0*pow(p[6][0],2.0)*derivative(aux_K,u[0],1.0)*derivative(aux_K,x[3],1.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,4.0)-4.0*pow(p[6][0],2.0)*derivative(aux_K,u[0],1.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_alpha_R_Ux/pow(1.0+aux_K,3.0))/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0));
    if (isNANorINF(daux_diff_2gamma_R_Ux_thr))
        {
            PRNT("    @k %d: daux_diff_2gamma_R_Ux_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2gamma_R_Ux_thr);
            {
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        u[0]= %g\n",u[0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        daux_diff_alpha_R_Ux= %g\n",daux_diff_alpha_R_Ux);
            }
            return(0);
        }
// aux_deriv[43]
    daux_diff_gamma_R_Uy=pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_alpha_R_Uy/pow(1.0+aux_K,2.0)/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0));
    if (isNANorINF(daux_diff_gamma_R_Uy))
        {
            PRNT("    @k %d: daux_diff_gamma_R_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_gamma_R_Uy);
            {
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        daux_diff_alpha_R_Uy= %g\n",daux_diff_alpha_R_Uy);
            }
            return(0);
        }
// aux_deriv[44]
    daux_diff_2gamma_R_Uy_Uy=pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_2alpha_R_Uy_Uy/pow(1.0+aux_K,2.0)/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))-pow(p[6][0],4.0)*pow(sec(aux_alpha_R),4.0)*pow(tan(aux_alpha_R),2.0)*pow(daux_diff_alpha_R_Uy,2.0)/pow(1.0+aux_K,4.0)/pow(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0),1.5)+pow(p[6][0],2.0)*pow(sec(aux_alpha_R),4.0)*pow(daux_diff_alpha_R_Uy,2.0)/pow(1.0+aux_K,2.0)/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))+2.0*pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*pow(tan(aux_alpha_R),2.0)*pow(daux_diff_alpha_R_Uy,2.0)/pow(1.0+aux_K,2.0)/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0));
    if (isNANorINF(daux_diff_2gamma_R_Uy_Uy))
        {
            PRNT("    @k %d: daux_diff_2gamma_R_Uy_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2gamma_R_Uy_Uy);
            {
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        daux_diff_2alpha_R_Uy_Uy= %g\n",daux_diff_2alpha_R_Uy_Uy);
                PRNT("        daux_diff_alpha_R_Uy= %g\n",daux_diff_alpha_R_Uy);
            }
            return(0);
        }
// aux_deriv[45]
    daux_diff_2gamma_R_Uy_r=pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_2alpha_R_Uy_r/pow(1.0+aux_K,2.0)/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))-pow(p[6][0],4.0)*pow(sec(aux_alpha_R),4.0)*pow(tan(aux_alpha_R),2.0)*daux_diff_alpha_R_Uy*daux_diff_alpha_R_r/pow(1.0+aux_K,4.0)/pow(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0),1.5)+pow(p[6][0],2.0)*pow(sec(aux_alpha_R),4.0)*daux_diff_alpha_R_Uy*daux_diff_alpha_R_r/pow(1.0+aux_K,2.0)/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))+2.0*pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*pow(tan(aux_alpha_R),2.0)*daux_diff_alpha_R_Uy*daux_diff_alpha_R_r/pow(1.0+aux_K,2.0)/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0));
    if (isNANorINF(daux_diff_2gamma_R_Uy_r))
        {
            PRNT("    @k %d: daux_diff_2gamma_R_Uy_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2gamma_R_Uy_r);
            {
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        daux_diff_2alpha_R_Uy_r= %g\n",daux_diff_2alpha_R_Uy_r);
                PRNT("        daux_diff_alpha_R_Uy= %g\n",daux_diff_alpha_R_Uy);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
            }
            return(0);
        }
// aux_deriv[46]
    daux_diff_2gamma_R_Uy_thr=-0.5*pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*(-2.0*pow(p[7][0],2.0)*pow(aux_K,2.0)*derivative(aux_K,u[0],1.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[7][0],2.0)*aux_K*derivative(aux_K,u[0],1.0)/pow(1.0+aux_K,2.0)-2.0*pow(p[6][0],2.0)*derivative(aux_K,u[0],1.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,3.0))*daux_diff_alpha_R_Uy/pow(1.0+aux_K,2.0)/pow(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0),1.5)-2.0*pow(p[6][0],2.0)*derivative(aux_K,u[0],1.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_alpha_R_Uy/pow(1.0+aux_K,3.0)/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0));
    if (isNANorINF(daux_diff_2gamma_R_Uy_thr))
        {
            PRNT("    @k %d: daux_diff_2gamma_R_Uy_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2gamma_R_Uy_thr);
            {
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        u[0]= %g\n",u[0]);
                PRNT("        daux_diff_alpha_R_Uy= %g\n",daux_diff_alpha_R_Uy);
            }
            return(0);
        }
// aux_deriv[47]
    daux_diff_gamma_R_r=pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_alpha_R_r/pow(1.0+aux_K,2.0)/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0));
    if (isNANorINF(daux_diff_gamma_R_r))
        {
            PRNT("    @k %d: daux_diff_gamma_R_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_gamma_R_r);
            {
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
            }
            return(0);
        }
// aux_deriv[48]
    daux_diff_2gamma_R_r_r=pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_2alpha_R_r_r/pow(1.0+aux_K,2.0)/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))-pow(p[6][0],4.0)*pow(sec(aux_alpha_R),4.0)*pow(tan(aux_alpha_R),2.0)*pow(daux_diff_alpha_R_r,2.0)/pow(1.0+aux_K,4.0)/pow(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0),1.5)+pow(p[6][0],2.0)*pow(sec(aux_alpha_R),4.0)*pow(daux_diff_alpha_R_r,2.0)/pow(1.0+aux_K,2.0)/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))+2.0*pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*pow(tan(aux_alpha_R),2.0)*pow(daux_diff_alpha_R_r,2.0)/pow(1.0+aux_K,2.0)/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0));
    if (isNANorINF(daux_diff_2gamma_R_r_r))
        {
            PRNT("    @k %d: daux_diff_2gamma_R_r_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2gamma_R_r_r);
            {
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        daux_diff_2alpha_R_r_r= %g\n",daux_diff_2alpha_R_r_r);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
            }
            return(0);
        }
// aux_deriv[49]
    daux_diff_2gamma_R_r_thr=-0.5*pow(p[6][0],2.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*(-2.0*pow(p[7][0],2.0)*pow(aux_K,2.0)*derivative(aux_K,u[0],1.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[7][0],2.0)*aux_K*derivative(aux_K,u[0],1.0)/pow(1.0+aux_K,2.0)-2.0*pow(p[6][0],2.0)*derivative(aux_K,u[0],1.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,3.0))*daux_diff_alpha_R_r/pow(1.0+aux_K,2.0)/pow(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0),1.5)-2.0*pow(p[6][0],2.0)*derivative(aux_K,u[0],1.0)*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_alpha_R_r/pow(1.0+aux_K,3.0)/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0));
    if (isNANorINF(daux_diff_2gamma_R_r_thr))
        {
            PRNT("    @k %d: daux_diff_2gamma_R_r_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2gamma_R_r_thr);
            {
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        u[0]= %g\n",u[0]);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
            }
            return(0);
        }
// aux_deriv[50]
    daux_diff_gamma_R_thr=0.5*(-2.0*pow(p[7][0],2.0)*pow(aux_K,2.0)*derivative(aux_K,u[0],1.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[7][0],2.0)*aux_K*derivative(aux_K,u[0],1.0)/pow(1.0+aux_K,2.0)-2.0*pow(p[6][0],2.0)*derivative(aux_K,u[0],1.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,3.0))/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0));
    if (isNANorINF(daux_diff_gamma_R_thr))
        {
            PRNT("    @k %d: daux_diff_gamma_R_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_gamma_R_thr);
            {
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
// aux_deriv[51]
    daux_diff_2gamma_R_thr_thr=0.5*(6.0*pow(p[7][0],2.0)*pow(aux_K,2.0)*pow(derivative(aux_K,u[0],1.0),2.0)/pow(1.0+aux_K,4.0)-8.0*pow(p[7][0],2.0)*aux_K*pow(derivative(aux_K,u[0],1.0),2.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[7][0],2.0)*pow(derivative(aux_K,u[0],1.0),2.0)/pow(1.0+aux_K,2.0)-2.0*pow(p[7][0],2.0)*pow(aux_K,2.0)*derivative(aux_K,u[0],2.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[7][0],2.0)*aux_K*derivative(aux_K,u[0],2.0)/pow(1.0+aux_K,2.0)+6.0*pow(p[6][0],2.0)*pow(derivative(aux_K,u[0],1.0),2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,4.0)-2.0*pow(p[6][0],2.0)*derivative(aux_K,u[0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,3.0))/sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))-0.25*pow(-2.0*pow(p[7][0],2.0)*pow(aux_K,2.0)*derivative(aux_K,u[0],1.0)/pow(1.0+aux_K,3.0)+2.0*pow(p[7][0],2.0)*aux_K*derivative(aux_K,u[0],1.0)/pow(1.0+aux_K,2.0)-2.0*pow(p[6][0],2.0)*derivative(aux_K,u[0],1.0)*pow(tan(aux_alpha_R),2.0)/pow(
     1.0+aux_K,3.0),2.0)/pow(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0),1.5);
    if (isNANorINF(daux_diff_2gamma_R_thr_thr))
        {
            PRNT("    @k %d: daux_diff_2gamma_R_thr_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2gamma_R_thr_thr);
            {
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
// aux_deriv[52]
    daux_diff_F_r_Ux=mcond(sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))<=3.0*p[1][0]*p[21][0],daux_diff_gamma_R_Ux+0.3333333333333333*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))*daux_diff_gamma_R_Ux/pow(p[1][0],2.0)/pow(p[21][0],2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))*daux_diff_gamma_R_Ux/p[1][0]/p[21][0],1,0.0);
    if (isNANorINF(daux_diff_F_r_Ux))
        {
            PRNT("    @k %d: daux_diff_F_r_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_F_r_Ux);
            {
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        p[1,0]= %g\n",p[1][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_gamma_R_Ux= %g\n",daux_diff_gamma_R_Ux);
                PRNT("        p[22,0]= %g\n",p[22][0]);
            }
            return(0);
        }
// aux_deriv[53]
    daux_diff_2F_r_Ux_Ux=mcond(sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))<=3.0*p[1][0]*p[21][0],daux_diff_2gamma_R_Ux_Ux+0.3333333333333333*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_Ux_Ux/pow(p[1][0],2.0)/pow(p[21][0],2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_Ux_Ux/p[1][0]/p[21][0]+0.3333333333333333*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(-2.0*pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)*derivative(abs(aux_K_),x[3],1.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),3.0)-2.0*pow(p[7][0],2.0)*pow(aux_K_,2.0)*derivative(
     abs(aux_K_),x[3],1.0)/pow(1.0+abs(aux_K_),3.0)+2.0*pow(p[7][0],2.0)*abs(aux_K_)*derivative(abs(aux_K_),x[3],1.0)/pow(1.0+abs(aux_K_),2.0)+2.0*(x[4]-x[5]*p[5][0])*pow(p[6][0],2.0)*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*daux_diff_alpha_R_Ux/(0.001+abs(x[3]))/pow(1.0+abs(aux_K_),2.0))*daux_diff_gamma_R_Ux/pow(p[1][0],2.0)/pow(p[21][0],2.0)-0.3333333333333333*(2.0-p[22][0]/p[21][0])*(-2.0*pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)*derivative(abs(aux_K_),x[3],1.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),3.0)-2.0*pow(p[7][0],2.0)*pow(aux_K_,2.0)*derivative(abs(aux_K_),x[3],1.0)/pow(1.0+abs(aux_K_),3.0)+2.0*pow(p[7][0],2.0)*abs(aux_K_)*derivative(abs(aux_K_),x[3],1.0)/pow(1.0+abs(aux_K_),2.0)+2.0*(x[4]-x[5]*p[5][0])*pow(p[6][0],2.0)*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*daux_diff_alpha_R_Ux/(0.001+abs(x[3]))/pow(1.0+abs(aux_K_),2.0))*daux_diff_gamma_R_Ux/p[1][0]/p[21][0]/sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)
     /pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0)),1,0.0);
    if (isNANorINF(daux_diff_2F_r_Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2F_r_Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2F_r_Ux_Ux);
            {
                PRNT("        aux_K_= %g\n",aux_K_);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[1,0]= %g\n",p[1][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_2gamma_R_Ux_Ux= %g\n",daux_diff_2gamma_R_Ux_Ux);
                PRNT("        p[22,0]= %g\n",p[22][0]);
                PRNT("        daux_diff_gamma_R_Ux= %g\n",daux_diff_gamma_R_Ux);
                PRNT("        daux_diff_alpha_R_Ux= %g\n",daux_diff_alpha_R_Ux);
            }
            return(0);
        }
// aux_deriv[54]
    daux_diff_2F_r_Ux_Uy=mcond(sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))<=3.0*p[1][0]*p[21][0],daux_diff_2gamma_R_Ux_Uy+0.3333333333333333*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_Ux_Uy/pow(p[1][0],2.0)/pow(p[21][0],2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_Ux_Uy/p[1][0]/p[21][0]+0.6666666666666666*(x[4]-x[5]*p[5][0])*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*daux_diff_alpha_R_Uy*daux_diff_gamma_R_Ux/pow(p[1][0],2.0)/pow(p[21][0],2.0)/(0.001+abs(
     x[3]))/pow(1.0+abs(aux_K_),2.0)-0.6666666666666666*(x[4]-x[5]*p[5][0])*pow(p[6][0],2.0)*(2.0-p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*daux_diff_alpha_R_Uy*daux_diff_gamma_R_Ux/p[1][0]/p[21][0]/(0.001+abs(x[3]))/pow(1.0+abs(aux_K_),2.0)/sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0)),1,0.0);
    if (isNANorINF(daux_diff_2F_r_Ux_Uy))
        {
            PRNT("    @k %d: daux_diff_2F_r_Ux_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2F_r_Ux_Uy);
            {
                PRNT("        aux_K_= %g\n",aux_K_);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[1,0]= %g\n",p[1][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_2gamma_R_Ux_Uy= %g\n",daux_diff_2gamma_R_Ux_Uy);
                PRNT("        daux_diff_alpha_R_Uy= %g\n",daux_diff_alpha_R_Uy);
                PRNT("        daux_diff_gamma_R_Ux= %g\n",daux_diff_gamma_R_Ux);
                PRNT("        p[22,0]= %g\n",p[22][0]);
            }
            return(0);
        }
// aux_deriv[55]
    daux_diff_2F_r_Ux_r=mcond(sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))<=3.0*p[1][0]*p[21][0],daux_diff_2gamma_R_Ux_r+0.3333333333333333*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_Ux_r/pow(p[1][0],2.0)/pow(p[21][0],2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_Ux_r/p[1][0]/p[21][0]+0.6666666666666666*(x[4]-x[5]*p[5][0])*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*daux_diff_alpha_R_r*daux_diff_gamma_R_Ux/pow(p[1][0],2.0)/pow(p[21][0],2.0)/(0.001+abs(x[3])
     )/pow(1.0+abs(aux_K_),2.0)-0.6666666666666666*(x[4]-x[5]*p[5][0])*pow(p[6][0],2.0)*(2.0-p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*daux_diff_alpha_R_r*daux_diff_gamma_R_Ux/p[1][0]/p[21][0]/(0.001+abs(x[3]))/pow(1.0+abs(aux_K_),2.0)/sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0)),1,0.0);
    if (isNANorINF(daux_diff_2F_r_Ux_r))
        {
            PRNT("    @k %d: daux_diff_2F_r_Ux_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2F_r_Ux_r);
            {
                PRNT("        aux_K_= %g\n",aux_K_);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[1,0]= %g\n",p[1][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_2gamma_R_Ux_r= %g\n",daux_diff_2gamma_R_Ux_r);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
                PRNT("        daux_diff_gamma_R_Ux= %g\n",daux_diff_gamma_R_Ux);
                PRNT("        p[22,0]= %g\n",p[22][0]);
            }
            return(0);
        }
// aux_deriv[56]
    daux_diff_2F_r_Ux_thr=mcond(sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))<=3.0*p[1][0]*p[21][0],daux_diff_2gamma_R_Ux_thr+0.3333333333333333*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_Ux_thr/pow(p[1][0],2.0)/pow(p[21][0],2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_Ux_thr/p[1][0]/p[21][0]+0.3333333333333333*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(-2.0*pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)*derivative(abs(aux_K_),u[0],1.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),3.0)-2.0*pow(p[7][0],2.0)*pow(aux_K_,2.0)*
     derivative(abs(aux_K_),u[0],1.0)/pow(1.0+abs(aux_K_),3.0)+2.0*pow(p[7][0],2.0)*abs(aux_K_)*derivative(abs(aux_K_),u[0],1.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_gamma_R_Ux/pow(p[1][0],2.0)/pow(p[21][0],2.0)-0.3333333333333333*(2.0-p[22][0]/p[21][0])*(-2.0*pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)*derivative(abs(aux_K_),u[0],1.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),3.0)-2.0*pow(p[7][0],2.0)*pow(aux_K_,2.0)*derivative(abs(aux_K_),u[0],1.0)/pow(1.0+abs(aux_K_),3.0)+2.0*pow(p[7][0],2.0)*abs(aux_K_)*derivative(abs(aux_K_),u[0],1.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_gamma_R_Ux/p[1][0]/p[21][0]/sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0)),1,0.0);
    if (isNANorINF(daux_diff_2F_r_Ux_thr))
        {
            PRNT("    @k %d: daux_diff_2F_r_Ux_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2F_r_Ux_thr);
            {
                PRNT("        aux_K_= %g\n",aux_K_);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[1,0]= %g\n",p[1][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_2gamma_R_Ux_thr= %g\n",daux_diff_2gamma_R_Ux_thr);
                PRNT("        p[22,0]= %g\n",p[22][0]);
                PRNT("        daux_diff_gamma_R_Ux= %g\n",daux_diff_gamma_R_Ux);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
// aux_deriv[57]
    daux_diff_F_r_Uy=mcond(sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))<=3.0*p[1][0]*p[21][0],daux_diff_gamma_R_Uy+0.3333333333333333*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))*daux_diff_gamma_R_Uy/pow(p[1][0],2.0)/pow(p[21][0],2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))*daux_diff_gamma_R_Uy/p[1][0]/p[21][0],1,0.0);
    if (isNANorINF(daux_diff_F_r_Uy))
        {
            PRNT("    @k %d: daux_diff_F_r_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_F_r_Uy);
            {
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        p[1,0]= %g\n",p[1][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_gamma_R_Uy= %g\n",daux_diff_gamma_R_Uy);
                PRNT("        p[22,0]= %g\n",p[22][0]);
            }
            return(0);
        }
// aux_deriv[58]
    daux_diff_2F_r_Uy_Uy=mcond(sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))<=3.0*p[1][0]*p[21][0],daux_diff_2gamma_R_Uy_Uy+0.3333333333333333*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_Uy_Uy/pow(p[1][0],2.0)/pow(p[21][0],2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_Uy_Uy/p[1][0]/p[21][0]+0.6666666666666666*(x[4]-x[5]*p[5][0])*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*daux_diff_alpha_R_Uy*daux_diff_gamma_R_Uy/pow(p[1][0],2.0)/pow(p[21][0],2.0)/(0.001+abs(
     x[3]))/pow(1.0+abs(aux_K_),2.0)-0.6666666666666666*(x[4]-x[5]*p[5][0])*pow(p[6][0],2.0)*(2.0-p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*daux_diff_alpha_R_Uy*daux_diff_gamma_R_Uy/p[1][0]/p[21][0]/(0.001+abs(x[3]))/pow(1.0+abs(aux_K_),2.0)/sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0)),1,0.0);
    if (isNANorINF(daux_diff_2F_r_Uy_Uy))
        {
            PRNT("    @k %d: daux_diff_2F_r_Uy_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2F_r_Uy_Uy);
            {
                PRNT("        aux_K_= %g\n",aux_K_);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[1,0]= %g\n",p[1][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_2gamma_R_Uy_Uy= %g\n",daux_diff_2gamma_R_Uy_Uy);
                PRNT("        daux_diff_alpha_R_Uy= %g\n",daux_diff_alpha_R_Uy);
                PRNT("        daux_diff_gamma_R_Uy= %g\n",daux_diff_gamma_R_Uy);
                PRNT("        p[22,0]= %g\n",p[22][0]);
            }
            return(0);
        }
// aux_deriv[59]
    daux_diff_2F_r_Uy_r=mcond(sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))<=3.0*p[1][0]*p[21][0],daux_diff_2gamma_R_Uy_r+0.3333333333333333*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_Uy_r/pow(p[1][0],2.0)/pow(p[21][0],2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_Uy_r/p[1][0]/p[21][0]+0.6666666666666666*(x[4]-x[5]*p[5][0])*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*daux_diff_alpha_R_r*daux_diff_gamma_R_Uy/pow(p[1][0],2.0)/pow(p[21][0],2.0)/(0.001+abs(x[3])
     )/pow(1.0+abs(aux_K_),2.0)-0.6666666666666666*(x[4]-x[5]*p[5][0])*pow(p[6][0],2.0)*(2.0-p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*daux_diff_alpha_R_r*daux_diff_gamma_R_Uy/p[1][0]/p[21][0]/(0.001+abs(x[3]))/pow(1.0+abs(aux_K_),2.0)/sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0)),1,0.0);
    if (isNANorINF(daux_diff_2F_r_Uy_r))
        {
            PRNT("    @k %d: daux_diff_2F_r_Uy_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2F_r_Uy_r);
            {
                PRNT("        aux_K_= %g\n",aux_K_);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[1,0]= %g\n",p[1][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_2gamma_R_Uy_r= %g\n",daux_diff_2gamma_R_Uy_r);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
                PRNT("        daux_diff_gamma_R_Uy= %g\n",daux_diff_gamma_R_Uy);
                PRNT("        p[22,0]= %g\n",p[22][0]);
            }
            return(0);
        }
// aux_deriv[60]
    daux_diff_2F_r_Uy_thr=mcond(sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))<=3.0*p[1][0]*p[21][0],daux_diff_2gamma_R_Uy_thr+0.3333333333333333*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_Uy_thr/pow(p[1][0],2.0)/pow(p[21][0],2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_Uy_thr/p[1][0]/p[21][0]+0.3333333333333333*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(-2.0*pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)*derivative(abs(aux_K_),u[0],1.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),3.0)-2.0*pow(p[7][0],2.0)*pow(aux_K_,2.0)*
     derivative(abs(aux_K_),u[0],1.0)/pow(1.0+abs(aux_K_),3.0)+2.0*pow(p[7][0],2.0)*abs(aux_K_)*derivative(abs(aux_K_),u[0],1.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_gamma_R_Uy/pow(p[1][0],2.0)/pow(p[21][0],2.0)-0.3333333333333333*(2.0-p[22][0]/p[21][0])*(-2.0*pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)*derivative(abs(aux_K_),u[0],1.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),3.0)-2.0*pow(p[7][0],2.0)*pow(aux_K_,2.0)*derivative(abs(aux_K_),u[0],1.0)/pow(1.0+abs(aux_K_),3.0)+2.0*pow(p[7][0],2.0)*abs(aux_K_)*derivative(abs(aux_K_),u[0],1.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_gamma_R_Uy/p[1][0]/p[21][0]/sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0)),1,0.0);
    if (isNANorINF(daux_diff_2F_r_Uy_thr))
        {
            PRNT("    @k %d: daux_diff_2F_r_Uy_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2F_r_Uy_thr);
            {
                PRNT("        aux_K_= %g\n",aux_K_);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[1,0]= %g\n",p[1][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_2gamma_R_Uy_thr= %g\n",daux_diff_2gamma_R_Uy_thr);
                PRNT("        p[22,0]= %g\n",p[22][0]);
                PRNT("        daux_diff_gamma_R_Uy= %g\n",daux_diff_gamma_R_Uy);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
// aux_deriv[61]
    daux_diff_F_r_r=mcond(sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))<=3.0*p[1][0]*p[21][0],daux_diff_gamma_R_r+0.3333333333333333*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))*daux_diff_gamma_R_r/pow(p[1][0],2.0)/pow(p[21][0],2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))*daux_diff_gamma_R_r/p[1][0]/p[21][0],1,0.0);
    if (isNANorINF(daux_diff_F_r_r))
        {
            PRNT("    @k %d: daux_diff_F_r_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_F_r_r);
            {
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        p[1,0]= %g\n",p[1][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_gamma_R_r= %g\n",daux_diff_gamma_R_r);
                PRNT("        p[22,0]= %g\n",p[22][0]);
            }
            return(0);
        }
// aux_deriv[62]
    daux_diff_2F_r_r_r=mcond(sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))<=3.0*p[1][0]*p[21][0],daux_diff_2gamma_R_r_r+0.3333333333333333*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_r_r/pow(p[1][0],2.0)/pow(p[21][0],2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_r_r/p[1][0]/p[21][0]+0.6666666666666666*(x[4]-x[5]*p[5][0])*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*daux_diff_alpha_R_r*daux_diff_gamma_R_r/pow(p[1][0],2.0)/pow(p[21][0],2.0)/(0.001+abs(x[3]))/pow
     (1.0+abs(aux_K_),2.0)-0.6666666666666666*(x[4]-x[5]*p[5][0])*pow(p[6][0],2.0)*(2.0-p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*daux_diff_alpha_R_r*daux_diff_gamma_R_r/p[1][0]/p[21][0]/(0.001+abs(x[3]))/pow(1.0+abs(aux_K_),2.0)/sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0)),1,0.0);
    if (isNANorINF(daux_diff_2F_r_r_r))
        {
            PRNT("    @k %d: daux_diff_2F_r_r_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2F_r_r_r);
            {
                PRNT("        aux_K_= %g\n",aux_K_);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[1,0]= %g\n",p[1][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_2gamma_R_r_r= %g\n",daux_diff_2gamma_R_r_r);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
                PRNT("        daux_diff_gamma_R_r= %g\n",daux_diff_gamma_R_r);
                PRNT("        p[22,0]= %g\n",p[22][0]);
            }
            return(0);
        }
// aux_deriv[63]
    daux_diff_2F_r_r_thr=mcond(sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))<=3.0*p[1][0]*p[21][0],daux_diff_2gamma_R_r_thr+0.3333333333333333*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_r_thr/pow(p[1][0],2.0)/pow(p[21][0],2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_r_thr/p[1][0]/p[21][0]+0.3333333333333333*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(-2.0*pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)*derivative(abs(aux_K_),u[0],1.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),3.0)-2.0*pow(p[7][0],2.0)*pow(aux_K_,2.0)*derivative(
     abs(aux_K_),u[0],1.0)/pow(1.0+abs(aux_K_),3.0)+2.0*pow(p[7][0],2.0)*abs(aux_K_)*derivative(abs(aux_K_),u[0],1.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_gamma_R_r/pow(p[1][0],2.0)/pow(p[21][0],2.0)-0.3333333333333333*(2.0-p[22][0]/p[21][0])*(-2.0*pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)*derivative(abs(aux_K_),u[0],1.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),3.0)-2.0*pow(p[7][0],2.0)*pow(aux_K_,2.0)*derivative(abs(aux_K_),u[0],1.0)/pow(1.0+abs(aux_K_),3.0)+2.0*pow(p[7][0],2.0)*abs(aux_K_)*derivative(abs(aux_K_),u[0],1.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_gamma_R_r/p[1][0]/p[21][0]/sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0)),1,0.0);
    if (isNANorINF(daux_diff_2F_r_r_thr))
        {
            PRNT("    @k %d: daux_diff_2F_r_r_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2F_r_r_thr);
            {
                PRNT("        aux_K_= %g\n",aux_K_);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[1,0]= %g\n",p[1][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_2gamma_R_r_thr= %g\n",daux_diff_2gamma_R_r_thr);
                PRNT("        p[22,0]= %g\n",p[22][0]);
                PRNT("        daux_diff_gamma_R_r= %g\n",daux_diff_gamma_R_r);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
// aux_deriv[64]
    daux_diff_F_r_thr=mcond(sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))<=3.0*p[1][0]*p[21][0],daux_diff_gamma_R_thr+0.3333333333333333*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))*daux_diff_gamma_R_thr/pow(p[1][0],2.0)/pow(p[21][0],2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*sqrt(pow(p[7][0],2.0)*pow(aux_K,2.0)/pow(1.0+aux_K,2.0)+pow(p[6][0],2.0)*pow(tan(aux_alpha_R),2.0)/pow(1.0+aux_K,2.0))*daux_diff_gamma_R_thr/p[1][0]/p[21][0],1,0.0);
    if (isNANorINF(daux_diff_F_r_thr))
        {
            PRNT("    @k %d: daux_diff_F_r_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_F_r_thr);
            {
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        p[1,0]= %g\n",p[1][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_gamma_R_thr= %g\n",daux_diff_gamma_R_thr);
                PRNT("        p[22,0]= %g\n",p[22][0]);
            }
            return(0);
        }
// aux_deriv[65]
    daux_diff_2F_r_thr_thr=mcond(sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))<=3.0*p[1][0]*p[21][0],daux_diff_2gamma_R_thr_thr+0.3333333333333333*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_thr_thr/pow(p[1][0],2.0)/pow(p[21][0],2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_2gamma_R_thr_thr/p[1][0]/p[21][0]+0.3333333333333333*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(-2.0*pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)*derivative(abs(aux_K_),u[0],1.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),3.0)-2.0*pow(p[7][0],2.0)*pow(aux_K_,2.0)*
     derivative(abs(aux_K_),u[0],1.0)/pow(1.0+abs(aux_K_),3.0)+2.0*pow(p[7][0],2.0)*abs(aux_K_)*derivative(abs(aux_K_),u[0],1.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_gamma_R_thr/pow(p[1][0],2.0)/pow(p[21][0],2.0)-0.3333333333333333*(2.0-p[22][0]/p[21][0])*(-2.0*pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)*derivative(abs(aux_K_),u[0],1.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),3.0)-2.0*pow(p[7][0],2.0)*pow(aux_K_,2.0)*derivative(abs(aux_K_),u[0],1.0)/pow(1.0+abs(aux_K_),3.0)+2.0*pow(p[7][0],2.0)*abs(aux_K_)*derivative(abs(aux_K_),u[0],1.0)/pow(1.0+abs(aux_K_),2.0))*daux_diff_gamma_R_thr/p[1][0]/p[21][0]/sqrt(pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)/pow(0.001+abs(x[3]),2.0)/pow(1.0+abs(aux_K_),2.0)+pow(p[7][0],2.0)*pow(aux_K_,2.0)/pow(1.0+abs(aux_K_),2.0)),1,0.0);
    if (isNANorINF(daux_diff_2F_r_thr_thr))
        {
            PRNT("    @k %d: daux_diff_2F_r_thr_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2F_r_thr_thr);
            {
                PRNT("        aux_K_= %g\n",aux_K_);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[1,0]= %g\n",p[1][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_2gamma_R_thr_thr= %g\n",daux_diff_2gamma_R_thr_thr);
                PRNT("        p[22,0]= %g\n",p[22][0]);
                PRNT("        daux_diff_gamma_R_thr= %g\n",daux_diff_gamma_R_thr);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
// aux_deriv[66]
    daux_diff_Fyr_Ux=-p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_Ux/aux_gamma_R-p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_Ux/aux_gamma_R+p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_gamma_R_Ux/pow(aux_gamma_R,2.0);
    if (isNANorINF(daux_diff_Fyr_Ux))
        {
            PRNT("    @k %d: daux_diff_Fyr_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_Fyr_Ux);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_gamma_R_Ux= %g\n",daux_diff_gamma_R_Ux);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_F_r_Ux= %g\n",daux_diff_F_r_Ux);
                PRNT("        daux_diff_alpha_R_Ux= %g\n",daux_diff_alpha_R_Ux);
            }
            return(0);
        }
// aux_deriv[67]
    daux_diff_2Fyr_Ux_Ux=-p[6][0]*tan(aux_alpha_R)*daux_diff_2F_r_Ux_Ux/aux_gamma_R-p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_2alpha_R_Ux_Ux/aux_gamma_R+p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_2gamma_R_Ux_Ux/pow(aux_gamma_R,2.0)-2.0*p[6][0]*pow(sec(aux_alpha_R),2.0)*daux_diff_F_r_Ux*daux_diff_alpha_R_Ux/aux_gamma_R-2.0*p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*pow(daux_diff_alpha_R_Ux,2.0)/aux_gamma_R+2.0*p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_Ux*daux_diff_gamma_R_Ux/pow(aux_gamma_R,2.0)+2.0*p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_Ux*daux_diff_gamma_R_Ux/pow(aux_gamma_R,2.0)-2.0*p[6][0]*aux_F_r*tan(aux_alpha_R)*pow(daux_diff_gamma_R_Ux,2.0)/pow(aux_gamma_R,3.0);
    if (isNANorINF(daux_diff_2Fyr_Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2Fyr_Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyr_Ux_Ux);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_gamma_R_Ux= %g\n",daux_diff_gamma_R_Ux);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_2gamma_R_Ux_Ux= %g\n",daux_diff_2gamma_R_Ux_Ux);
                PRNT("        daux_diff_F_r_Ux= %g\n",daux_diff_F_r_Ux);
                PRNT("        daux_diff_alpha_R_Ux= %g\n",daux_diff_alpha_R_Ux);
                PRNT("        daux_diff_2F_r_Ux_Ux= %g\n",daux_diff_2F_r_Ux_Ux);
                PRNT("        daux_diff_2alpha_R_Ux_Ux= %g\n",daux_diff_2alpha_R_Ux_Ux);
            }
            return(0);
        }
// aux_deriv[68]
    daux_diff_2Fyr_Ux_Uy=-p[6][0]*tan(aux_alpha_R)*daux_diff_2F_r_Ux_Uy/aux_gamma_R-p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_2alpha_R_Ux_Uy/aux_gamma_R+p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_2gamma_R_Ux_Uy/pow(aux_gamma_R,2.0)-p[6][0]*pow(sec(aux_alpha_R),2.0)*daux_diff_F_r_Uy*daux_diff_alpha_R_Ux/aux_gamma_R-p[6][0]*pow(sec(aux_alpha_R),2.0)*daux_diff_F_r_Ux*daux_diff_alpha_R_Uy/aux_gamma_R-2.0*p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_alpha_R_Ux*daux_diff_alpha_R_Uy/aux_gamma_R+p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_Uy*daux_diff_gamma_R_Ux/pow(aux_gamma_R,2.0)+p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_Uy*daux_diff_gamma_R_Ux/pow(aux_gamma_R,2.0)+p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_Ux*daux_diff_gamma_R_Uy/pow(aux_gamma_R,2.0)+p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_Ux*daux_diff_gamma_R_Uy/pow(aux_gamma_R,2.0)-2.0*p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_gamma_R_Ux*daux_diff_gamma_R_Uy/pow(aux_gamma_R,
     3.0);
    if (isNANorINF(daux_diff_2Fyr_Ux_Uy))
        {
            PRNT("    @k %d: daux_diff_2Fyr_Ux_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyr_Ux_Uy);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_gamma_R_Ux= %g\n",daux_diff_gamma_R_Ux);
                PRNT("        daux_diff_gamma_R_Uy= %g\n",daux_diff_gamma_R_Uy);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_2gamma_R_Ux_Uy= %g\n",daux_diff_2gamma_R_Ux_Uy);
                PRNT("        daux_diff_F_r_Uy= %g\n",daux_diff_F_r_Uy);
                PRNT("        daux_diff_alpha_R_Uy= %g\n",daux_diff_alpha_R_Uy);
                PRNT("        daux_diff_F_r_Ux= %g\n",daux_diff_F_r_Ux);
                PRNT("        daux_diff_alpha_R_Ux= %g\n",daux_diff_alpha_R_Ux);
                PRNT("        daux_diff_2F_r_Ux_Uy= %g\n",daux_diff_2F_r_Ux_Uy);
                PRNT("        daux_diff_2alpha_R_Ux_Uy= %g\n",daux_diff_2alpha_R_Ux_Uy);
            }
            return(0);
        }
// aux_deriv[69]
    daux_diff_2Fyr_Ux_r=-p[6][0]*tan(aux_alpha_R)*daux_diff_2F_r_Ux_r/aux_gamma_R-p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_2alpha_R_Ux_r/aux_gamma_R+p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_2gamma_R_Ux_r/pow(aux_gamma_R,2.0)-p[6][0]*pow(sec(aux_alpha_R),2.0)*daux_diff_F_r_r*daux_diff_alpha_R_Ux/aux_gamma_R-p[6][0]*pow(sec(aux_alpha_R),2.0)*daux_diff_F_r_Ux*daux_diff_alpha_R_r/aux_gamma_R-2.0*p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_alpha_R_Ux*daux_diff_alpha_R_r/aux_gamma_R+p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_r*daux_diff_gamma_R_Ux/pow(aux_gamma_R,2.0)+p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_r*daux_diff_gamma_R_Ux/pow(aux_gamma_R,2.0)+p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_Ux*daux_diff_gamma_R_r/pow(aux_gamma_R,2.0)+p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_Ux*daux_diff_gamma_R_r/pow(aux_gamma_R,2.0)-2.0*p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_gamma_R_Ux*daux_diff_gamma_R_r/pow(aux_gamma_R,3.0);
    if (isNANorINF(daux_diff_2Fyr_Ux_r))
        {
            PRNT("    @k %d: daux_diff_2Fyr_Ux_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyr_Ux_r);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_gamma_R_Ux= %g\n",daux_diff_gamma_R_Ux);
                PRNT("        daux_diff_gamma_R_r= %g\n",daux_diff_gamma_R_r);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_2gamma_R_Ux_r= %g\n",daux_diff_2gamma_R_Ux_r);
                PRNT("        daux_diff_F_r_r= %g\n",daux_diff_F_r_r);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
                PRNT("        daux_diff_F_r_Ux= %g\n",daux_diff_F_r_Ux);
                PRNT("        daux_diff_alpha_R_Ux= %g\n",daux_diff_alpha_R_Ux);
                PRNT("        daux_diff_2F_r_Ux_r= %g\n",daux_diff_2F_r_Ux_r);
                PRNT("        daux_diff_2alpha_R_Ux_r= %g\n",daux_diff_2alpha_R_Ux_r);
            }
            return(0);
        }
// aux_deriv[70]
    daux_diff_2Fyr_Ux_thr=-p[6][0]*tan(aux_alpha_R)*daux_diff_2F_r_Ux_thr/aux_gamma_R+p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_2gamma_R_Ux_thr/pow(aux_gamma_R,2.0)-p[6][0]*pow(sec(aux_alpha_R),2.0)*daux_diff_F_r_thr*daux_diff_alpha_R_Ux/aux_gamma_R+p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_thr*daux_diff_gamma_R_Ux/pow(aux_gamma_R,2.0)+p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_Ux*daux_diff_gamma_R_thr/pow(aux_gamma_R,2.0)+p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_Ux*daux_diff_gamma_R_thr/pow(aux_gamma_R,2.0)-2.0*p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_gamma_R_Ux*daux_diff_gamma_R_thr/pow(aux_gamma_R,3.0);
    if (isNANorINF(daux_diff_2Fyr_Ux_thr))
        {
            PRNT("    @k %d: daux_diff_2Fyr_Ux_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyr_Ux_thr);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_gamma_R_Ux= %g\n",daux_diff_gamma_R_Ux);
                PRNT("        daux_diff_gamma_R_thr= %g\n",daux_diff_gamma_R_thr);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_2gamma_R_Ux_thr= %g\n",daux_diff_2gamma_R_Ux_thr);
                PRNT("        daux_diff_F_r_thr= %g\n",daux_diff_F_r_thr);
                PRNT("        daux_diff_F_r_Ux= %g\n",daux_diff_F_r_Ux);
                PRNT("        daux_diff_alpha_R_Ux= %g\n",daux_diff_alpha_R_Ux);
                PRNT("        daux_diff_2F_r_Ux_thr= %g\n",daux_diff_2F_r_Ux_thr);
            }
            return(0);
        }
// aux_deriv[71]
    daux_diff_Fyr_Uy=-p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_Uy/aux_gamma_R-p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_Uy/aux_gamma_R+p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_gamma_R_Uy/pow(aux_gamma_R,2.0);
    if (isNANorINF(daux_diff_Fyr_Uy))
        {
            PRNT("    @k %d: daux_diff_Fyr_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_Fyr_Uy);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_gamma_R_Uy= %g\n",daux_diff_gamma_R_Uy);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_F_r_Uy= %g\n",daux_diff_F_r_Uy);
                PRNT("        daux_diff_alpha_R_Uy= %g\n",daux_diff_alpha_R_Uy);
            }
            return(0);
        }
// aux_deriv[72]
    daux_diff_2Fyr_Uy_Uy=-p[6][0]*tan(aux_alpha_R)*daux_diff_2F_r_Uy_Uy/aux_gamma_R-p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_2alpha_R_Uy_Uy/aux_gamma_R+p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_2gamma_R_Uy_Uy/pow(aux_gamma_R,2.0)-2.0*p[6][0]*pow(sec(aux_alpha_R),2.0)*daux_diff_F_r_Uy*daux_diff_alpha_R_Uy/aux_gamma_R-2.0*p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*pow(daux_diff_alpha_R_Uy,2.0)/aux_gamma_R+2.0*p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_Uy*daux_diff_gamma_R_Uy/pow(aux_gamma_R,2.0)+2.0*p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_Uy*daux_diff_gamma_R_Uy/pow(aux_gamma_R,2.0)-2.0*p[6][0]*aux_F_r*tan(aux_alpha_R)*pow(daux_diff_gamma_R_Uy,2.0)/pow(aux_gamma_R,3.0);
    if (isNANorINF(daux_diff_2Fyr_Uy_Uy))
        {
            PRNT("    @k %d: daux_diff_2Fyr_Uy_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyr_Uy_Uy);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_gamma_R_Uy= %g\n",daux_diff_gamma_R_Uy);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_2gamma_R_Uy_Uy= %g\n",daux_diff_2gamma_R_Uy_Uy);
                PRNT("        daux_diff_F_r_Uy= %g\n",daux_diff_F_r_Uy);
                PRNT("        daux_diff_alpha_R_Uy= %g\n",daux_diff_alpha_R_Uy);
                PRNT("        daux_diff_2F_r_Uy_Uy= %g\n",daux_diff_2F_r_Uy_Uy);
                PRNT("        daux_diff_2alpha_R_Uy_Uy= %g\n",daux_diff_2alpha_R_Uy_Uy);
            }
            return(0);
        }
// aux_deriv[73]
    daux_diff_2Fyr_Uy_r=-p[6][0]*tan(aux_alpha_R)*daux_diff_2F_r_Uy_r/aux_gamma_R-p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_2alpha_R_Uy_r/aux_gamma_R+p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_2gamma_R_Uy_r/pow(aux_gamma_R,2.0)-p[6][0]*pow(sec(aux_alpha_R),2.0)*daux_diff_F_r_r*daux_diff_alpha_R_Uy/aux_gamma_R-p[6][0]*pow(sec(aux_alpha_R),2.0)*daux_diff_F_r_Uy*daux_diff_alpha_R_r/aux_gamma_R-2.0*p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*daux_diff_alpha_R_Uy*daux_diff_alpha_R_r/aux_gamma_R+p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_r*daux_diff_gamma_R_Uy/pow(aux_gamma_R,2.0)+p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_r*daux_diff_gamma_R_Uy/pow(aux_gamma_R,2.0)+p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_Uy*daux_diff_gamma_R_r/pow(aux_gamma_R,2.0)+p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_Uy*daux_diff_gamma_R_r/pow(aux_gamma_R,2.0)-2.0*p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_gamma_R_Uy*daux_diff_gamma_R_r/pow(aux_gamma_R,3.0);
    if (isNANorINF(daux_diff_2Fyr_Uy_r))
        {
            PRNT("    @k %d: daux_diff_2Fyr_Uy_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyr_Uy_r);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_gamma_R_Uy= %g\n",daux_diff_gamma_R_Uy);
                PRNT("        daux_diff_gamma_R_r= %g\n",daux_diff_gamma_R_r);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_2gamma_R_Uy_r= %g\n",daux_diff_2gamma_R_Uy_r);
                PRNT("        daux_diff_F_r_r= %g\n",daux_diff_F_r_r);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
                PRNT("        daux_diff_F_r_Uy= %g\n",daux_diff_F_r_Uy);
                PRNT("        daux_diff_alpha_R_Uy= %g\n",daux_diff_alpha_R_Uy);
                PRNT("        daux_diff_2F_r_Uy_r= %g\n",daux_diff_2F_r_Uy_r);
                PRNT("        daux_diff_2alpha_R_Uy_r= %g\n",daux_diff_2alpha_R_Uy_r);
            }
            return(0);
        }
// aux_deriv[74]
    daux_diff_2Fyr_Uy_thr=-p[6][0]*tan(aux_alpha_R)*daux_diff_2F_r_Uy_thr/aux_gamma_R+p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_2gamma_R_Uy_thr/pow(aux_gamma_R,2.0)-p[6][0]*pow(sec(aux_alpha_R),2.0)*daux_diff_F_r_thr*daux_diff_alpha_R_Uy/aux_gamma_R+p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_thr*daux_diff_gamma_R_Uy/pow(aux_gamma_R,2.0)+p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_Uy*daux_diff_gamma_R_thr/pow(aux_gamma_R,2.0)+p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_Uy*daux_diff_gamma_R_thr/pow(aux_gamma_R,2.0)-2.0*p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_gamma_R_Uy*daux_diff_gamma_R_thr/pow(aux_gamma_R,3.0);
    if (isNANorINF(daux_diff_2Fyr_Uy_thr))
        {
            PRNT("    @k %d: daux_diff_2Fyr_Uy_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyr_Uy_thr);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_gamma_R_Uy= %g\n",daux_diff_gamma_R_Uy);
                PRNT("        daux_diff_gamma_R_thr= %g\n",daux_diff_gamma_R_thr);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_2gamma_R_Uy_thr= %g\n",daux_diff_2gamma_R_Uy_thr);
                PRNT("        daux_diff_F_r_thr= %g\n",daux_diff_F_r_thr);
                PRNT("        daux_diff_F_r_Uy= %g\n",daux_diff_F_r_Uy);
                PRNT("        daux_diff_alpha_R_Uy= %g\n",daux_diff_alpha_R_Uy);
                PRNT("        daux_diff_2F_r_Uy_thr= %g\n",daux_diff_2F_r_Uy_thr);
            }
            return(0);
        }
// aux_deriv[75]
    daux_diff_Fyr_r=-p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_r/aux_gamma_R-p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_r/aux_gamma_R+p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_gamma_R_r/pow(aux_gamma_R,2.0);
    if (isNANorINF(daux_diff_Fyr_r))
        {
            PRNT("    @k %d: daux_diff_Fyr_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_Fyr_r);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_gamma_R_r= %g\n",daux_diff_gamma_R_r);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_F_r_r= %g\n",daux_diff_F_r_r);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
            }
            return(0);
        }
// aux_deriv[76]
    daux_diff_2Fyr_r_r=-p[6][0]*tan(aux_alpha_R)*daux_diff_2F_r_r_r/aux_gamma_R-p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_2alpha_R_r_r/aux_gamma_R+p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_2gamma_R_r_r/pow(aux_gamma_R,2.0)-2.0*p[6][0]*pow(sec(aux_alpha_R),2.0)*daux_diff_F_r_r*daux_diff_alpha_R_r/aux_gamma_R-2.0*p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*pow(daux_diff_alpha_R_r,2.0)/aux_gamma_R+2.0*p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_r*daux_diff_gamma_R_r/pow(aux_gamma_R,2.0)+2.0*p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_r*daux_diff_gamma_R_r/pow(aux_gamma_R,2.0)-2.0*p[6][0]*aux_F_r*tan(aux_alpha_R)*pow(daux_diff_gamma_R_r,2.0)/pow(aux_gamma_R,3.0);
    if (isNANorINF(daux_diff_2Fyr_r_r))
        {
            PRNT("    @k %d: daux_diff_2Fyr_r_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyr_r_r);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_gamma_R_r= %g\n",daux_diff_gamma_R_r);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_2gamma_R_r_r= %g\n",daux_diff_2gamma_R_r_r);
                PRNT("        daux_diff_F_r_r= %g\n",daux_diff_F_r_r);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
                PRNT("        daux_diff_2F_r_r_r= %g\n",daux_diff_2F_r_r_r);
                PRNT("        daux_diff_2alpha_R_r_r= %g\n",daux_diff_2alpha_R_r_r);
            }
            return(0);
        }
// aux_deriv[77]
    daux_diff_2Fyr_r_thr=-p[6][0]*tan(aux_alpha_R)*daux_diff_2F_r_r_thr/aux_gamma_R+p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_2gamma_R_r_thr/pow(aux_gamma_R,2.0)-p[6][0]*pow(sec(aux_alpha_R),2.0)*daux_diff_F_r_thr*daux_diff_alpha_R_r/aux_gamma_R+p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_thr*daux_diff_gamma_R_r/pow(aux_gamma_R,2.0)+p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_r*daux_diff_gamma_R_thr/pow(aux_gamma_R,2.0)+p[6][0]*aux_F_r*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_r*daux_diff_gamma_R_thr/pow(aux_gamma_R,2.0)-2.0*p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_gamma_R_r*daux_diff_gamma_R_thr/pow(aux_gamma_R,3.0);
    if (isNANorINF(daux_diff_2Fyr_r_thr))
        {
            PRNT("    @k %d: daux_diff_2Fyr_r_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyr_r_thr);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_gamma_R_r= %g\n",daux_diff_gamma_R_r);
                PRNT("        daux_diff_gamma_R_thr= %g\n",daux_diff_gamma_R_thr);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_2gamma_R_r_thr= %g\n",daux_diff_2gamma_R_r_thr);
                PRNT("        daux_diff_F_r_thr= %g\n",daux_diff_F_r_thr);
                PRNT("        daux_diff_F_r_r= %g\n",daux_diff_F_r_r);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
                PRNT("        daux_diff_2F_r_r_thr= %g\n",daux_diff_2F_r_r_thr);
            }
            return(0);
        }
// aux_deriv[78]
    daux_diff_Fyr_thr=-p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_thr/aux_gamma_R+p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_gamma_R_thr/pow(aux_gamma_R,2.0);
    if (isNANorINF(daux_diff_Fyr_thr))
        {
            PRNT("    @k %d: daux_diff_Fyr_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_Fyr_thr);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_gamma_R_thr= %g\n",daux_diff_gamma_R_thr);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_F_r_thr= %g\n",daux_diff_F_r_thr);
            }
            return(0);
        }
// aux_deriv[79]
    daux_diff_2Fyr_thr_thr=-p[6][0]*tan(aux_alpha_R)*daux_diff_2F_r_thr_thr/aux_gamma_R+p[6][0]*aux_F_r*tan(aux_alpha_R)*daux_diff_2gamma_R_thr_thr/pow(aux_gamma_R,2.0)+2.0*p[6][0]*tan(aux_alpha_R)*daux_diff_F_r_thr*daux_diff_gamma_R_thr/pow(aux_gamma_R,2.0)-2.0*p[6][0]*aux_F_r*tan(aux_alpha_R)*pow(daux_diff_gamma_R_thr,2.0)/pow(aux_gamma_R,3.0);
    if (isNANorINF(daux_diff_2Fyr_thr_thr))
        {
            PRNT("    @k %d: daux_diff_2Fyr_thr_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyr_thr_thr);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_gamma_R_thr= %g\n",daux_diff_gamma_R_thr);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_2gamma_R_thr_thr= %g\n",daux_diff_2gamma_R_thr_thr);
                PRNT("        daux_diff_F_r_thr= %g\n",daux_diff_F_r_thr);
                PRNT("        daux_diff_2F_r_thr_thr= %g\n",daux_diff_2F_r_thr_thr);
            }
            return(0);
        }
// aux_deriv[80]
    daux_diff_dr_Ux=(p[4][0]*cos(u[1])*derivative(aux_Fyf,x[3],1.0)-p[5][0]*daux_diff_Fyr_Ux)/p[2][0];
    if (isNANorINF(daux_diff_dr_Ux))
        {
            PRNT("    @k %d: daux_diff_dr_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dr_Ux);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        daux_diff_Fyr_Ux= %g\n",daux_diff_Fyr_Ux);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[81]
    daux_diff_2dr_Ux_Ux=(p[4][0]*cos(u[1])*derivative(aux_Fyf,x[3],2.0)-p[5][0]*daux_diff_2Fyr_Ux_Ux)/p[2][0];
    if (isNANorINF(daux_diff_2dr_Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2dr_Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dr_Ux_Ux);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        daux_diff_2Fyr_Ux_Ux= %g\n",daux_diff_2Fyr_Ux_Ux);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[82]
    daux_diff_2dr_Ux_Uy=(p[4][0]*cos(u[1])*derivative(aux_Fyf,x[3],1.0,x[4],1.0)-p[5][0]*daux_diff_2Fyr_Ux_Uy)/p[2][0];
    if (isNANorINF(daux_diff_2dr_Ux_Uy))
        {
            PRNT("    @k %d: daux_diff_2dr_Ux_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dr_Ux_Uy);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        daux_diff_2Fyr_Ux_Uy= %g\n",daux_diff_2Fyr_Ux_Uy);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[83]
    daux_diff_2dr_Ux_r=(p[4][0]*cos(u[1])*derivative(aux_Fyf,x[3],1.0,x[5],1.0)-p[5][0]*daux_diff_2Fyr_Ux_r)/p[2][0];
    if (isNANorINF(daux_diff_2dr_Ux_r))
        {
            PRNT("    @k %d: daux_diff_2dr_Ux_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dr_Ux_r);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        daux_diff_2Fyr_Ux_r= %g\n",daux_diff_2Fyr_Ux_r);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[84]
    daux_diff_2dr_Ux_steer=(p[4][0]*cos(u[1])*derivative(aux_Fyf,u[1],1.0,x[3],1.0)-p[4][0]*sin(u[1])*derivative(aux_Fyf,x[3],1.0))/p[2][0];
    if (isNANorINF(daux_diff_2dr_Ux_steer))
        {
            PRNT("    @k %d: daux_diff_2dr_Ux_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dr_Ux_steer);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
            }
            return(0);
        }
// aux_deriv[85]
    daux_diff_2dr_Ux_thr=-p[5][0]*daux_diff_2Fyr_Ux_thr/p[2][0];
    if (isNANorINF(daux_diff_2dr_Ux_thr))
        {
            PRNT("    @k %d: daux_diff_2dr_Ux_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dr_Ux_thr);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        daux_diff_2Fyr_Ux_thr= %g\n",daux_diff_2Fyr_Ux_thr);
            }
            return(0);
        }
// aux_deriv[86]
    daux_diff_dr_Uy=(p[4][0]*cos(u[1])*derivative(aux_Fyf,x[4],1.0)-p[5][0]*daux_diff_Fyr_Uy)/p[2][0];
    if (isNANorINF(daux_diff_dr_Uy))
        {
            PRNT("    @k %d: daux_diff_dr_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dr_Uy);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        daux_diff_Fyr_Uy= %g\n",daux_diff_Fyr_Uy);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[87]
    daux_diff_2dr_Uy_Uy=(p[4][0]*cos(u[1])*derivative(aux_Fyf,x[4],2.0)-p[5][0]*daux_diff_2Fyr_Uy_Uy)/p[2][0];
    if (isNANorINF(daux_diff_2dr_Uy_Uy))
        {
            PRNT("    @k %d: daux_diff_2dr_Uy_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dr_Uy_Uy);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        daux_diff_2Fyr_Uy_Uy= %g\n",daux_diff_2Fyr_Uy_Uy);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[88]
    daux_diff_2dr_Uy_r=(p[4][0]*cos(u[1])*derivative(aux_Fyf,x[4],1.0,x[5],1.0)-p[5][0]*daux_diff_2Fyr_Uy_r)/p[2][0];
    if (isNANorINF(daux_diff_2dr_Uy_r))
        {
            PRNT("    @k %d: daux_diff_2dr_Uy_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dr_Uy_r);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        daux_diff_2Fyr_Uy_r= %g\n",daux_diff_2Fyr_Uy_r);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[89]
    daux_diff_2dr_Uy_steer=(p[4][0]*cos(u[1])*derivative(aux_Fyf,u[1],1.0,x[4],1.0)-p[4][0]*sin(u[1])*derivative(aux_Fyf,x[4],1.0))/p[2][0];
    if (isNANorINF(daux_diff_2dr_Uy_steer))
        {
            PRNT("    @k %d: daux_diff_2dr_Uy_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dr_Uy_steer);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
            }
            return(0);
        }
// aux_deriv[90]
    daux_diff_2dr_Uy_thr=-p[5][0]*daux_diff_2Fyr_Uy_thr/p[2][0];
    if (isNANorINF(daux_diff_2dr_Uy_thr))
        {
            PRNT("    @k %d: daux_diff_2dr_Uy_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dr_Uy_thr);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        daux_diff_2Fyr_Uy_thr= %g\n",daux_diff_2Fyr_Uy_thr);
            }
            return(0);
        }
// aux_deriv[91]
    daux_diff_dr_r=(p[4][0]*cos(u[1])*derivative(aux_Fyf,x[5],1.0)-p[5][0]*daux_diff_Fyr_r)/p[2][0];
    if (isNANorINF(daux_diff_dr_r))
        {
            PRNT("    @k %d: daux_diff_dr_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dr_r);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        daux_diff_Fyr_r= %g\n",daux_diff_Fyr_r);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[92]
    daux_diff_2dr_r_r=(p[4][0]*cos(u[1])*derivative(aux_Fyf,x[5],2.0)-p[5][0]*daux_diff_2Fyr_r_r)/p[2][0];
    if (isNANorINF(daux_diff_2dr_r_r))
        {
            PRNT("    @k %d: daux_diff_2dr_r_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dr_r_r);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        daux_diff_2Fyr_r_r= %g\n",daux_diff_2Fyr_r_r);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[93]
    daux_diff_2dr_r_steer=(p[4][0]*cos(u[1])*derivative(aux_Fyf,u[1],1.0,x[5],1.0)-p[4][0]*sin(u[1])*derivative(aux_Fyf,x[5],1.0))/p[2][0];
    if (isNANorINF(daux_diff_2dr_r_steer))
        {
            PRNT("    @k %d: daux_diff_2dr_r_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dr_r_steer);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
            }
            return(0);
        }
// aux_deriv[94]
    daux_diff_2dr_r_thr=-p[5][0]*daux_diff_2Fyr_r_thr/p[2][0];
    if (isNANorINF(daux_diff_2dr_r_thr))
        {
            PRNT("    @k %d: daux_diff_2dr_r_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dr_r_thr);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        daux_diff_2Fyr_r_thr= %g\n",daux_diff_2Fyr_r_thr);
            }
            return(0);
        }
// aux_deriv[95]
    daux_diff_dr_steer=(-p[4][0]*sin(u[1])*aux_Fyf+p[4][0]*cos(u[1])*derivative(aux_Fyf,u[1],1.0))/p[2][0];
    if (isNANorINF(daux_diff_dr_steer))
        {
            PRNT("    @k %d: daux_diff_dr_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dr_steer);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
            }
            return(0);
        }
// aux_deriv[96]
    daux_diff_2dr_steer_steer=(-p[4][0]*cos(u[1])*aux_Fyf-2.0*p[4][0]*sin(u[1])*derivative(aux_Fyf,u[1],1.0)+p[4][0]*cos(u[1])*derivative(aux_Fyf,u[1],2.0))/p[2][0];
    if (isNANorINF(daux_diff_2dr_steer_steer))
        {
            PRNT("    @k %d: daux_diff_2dr_steer_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dr_steer_steer);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[97]
    daux_diff_dr_thr=-p[5][0]*daux_diff_Fyr_thr/p[2][0];
    if (isNANorINF(daux_diff_dr_thr))
        {
            PRNT("    @k %d: daux_diff_dr_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dr_thr);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        daux_diff_Fyr_thr= %g\n",daux_diff_Fyr_thr);
            }
            return(0);
        }
// aux_deriv[98]
    daux_diff_2dr_thr_thr=-p[5][0]*daux_diff_2Fyr_thr_thr/p[2][0];
    if (isNANorINF(daux_diff_2dr_thr_thr))
        {
            PRNT("    @k %d: daux_diff_2dr_thr_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dr_thr_thr);
            {
                PRNT("        p[2,0]= %g\n",p[2][0]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        daux_diff_2Fyr_thr_thr= %g\n",daux_diff_2Fyr_thr_thr);
            }
            return(0);
        }
// aux_deriv[99]
    daux_diff_K__Ux=-1/(0.001+abs(x[3]))-(u[0]-x[3])*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0);
    if (isNANorINF(daux_diff_K__Ux))
        {
            PRNT("    @k %d: daux_diff_K__Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_K__Ux);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
// aux_deriv[100]
    daux_diff_2K__Ux_Ux=-(u[0]-x[3])*mcond(x[3]>0.0,0.0,1,0.0)/pow(0.001+abs(x[3]),2.0)+2.0*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)+2.0*(u[0]-x[3])*pow(mcond(x[3]>0.0,1.0,1,-1.0),2.0)/pow(0.001+abs(x[3]),3.0);
    if (isNANorINF(daux_diff_2K__Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2K__Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2K__Ux_Ux);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
// aux_deriv[101]
    daux_diff_2K__Ux_thr=-mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0);
    if (isNANorINF(daux_diff_2K__Ux_thr))
        {
            PRNT("    @k %d: daux_diff_2K__Ux_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2K__Ux_thr);
            {
                PRNT("        x[3]= %g\n",x[3]);
            }
            return(0);
        }
// aux_deriv[102]
    daux_diff_K__thr=1/(0.001+abs(x[3]));
    if (isNANorINF(daux_diff_K__thr))
        {
            PRNT("    @k %d: daux_diff_K__thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_K__thr);
            {
                PRNT("        x[3]= %g\n",x[3]);
            }
            return(0);
        }
// aux_deriv[103]
    daux_diff_K_Ux=mcond(aux_K_>0.0,1.0,1,-1.0)*daux_diff_K__Ux;
    if (isNANorINF(daux_diff_K_Ux))
        {
            PRNT("    @k %d: daux_diff_K_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_K_Ux);
            {
                PRNT("        aux_K_= %g\n",aux_K_);
                PRNT("        daux_diff_K__Ux= %g\n",daux_diff_K__Ux);
            }
            return(0);
        }
// aux_deriv[104]
    daux_diff_2K_Ux_Ux=mcond(aux_K_>0.0,1.0,1,-1.0)*daux_diff_2K__Ux_Ux+mcond((u[0]-x[3])/(0.001+abs(x[3]))>0.0,0.0,1,0.0)*daux_diff_K__Ux;
    if (isNANorINF(daux_diff_2K_Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2K_Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2K_Ux_Ux);
            {
                PRNT("        aux_K_= %g\n",aux_K_);
                PRNT("        daux_diff_2K__Ux_Ux= %g\n",daux_diff_2K__Ux_Ux);
                PRNT("        daux_diff_K__Ux= %g\n",daux_diff_K__Ux);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
// aux_deriv[105]
    daux_diff_2K_Ux_thr=mcond(aux_K_>0.0,1.0,1,-1.0)*daux_diff_2K__Ux_thr+mcond((u[0]-x[3])/(0.001+abs(x[3]))>0.0,0.0,1,0.0)*daux_diff_K__Ux;
    if (isNANorINF(daux_diff_2K_Ux_thr))
        {
            PRNT("    @k %d: daux_diff_2K_Ux_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2K_Ux_thr);
            {
                PRNT("        aux_K_= %g\n",aux_K_);
                PRNT("        daux_diff_2K__Ux_thr= %g\n",daux_diff_2K__Ux_thr);
                PRNT("        daux_diff_K__Ux= %g\n",daux_diff_K__Ux);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
// aux_deriv[106]
    daux_diff_K_thr=mcond(aux_K_>0.0,1.0,1,-1.0)*daux_diff_K__thr;
    if (isNANorINF(daux_diff_K_thr))
        {
            PRNT("    @k %d: daux_diff_K_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_K_thr);
            {
                PRNT("        aux_K_= %g\n",aux_K_);
                PRNT("        daux_diff_K__thr= %g\n",daux_diff_K__thr);
            }
            return(0);
        }
// aux_deriv[107]
    daux_diff_2K_thr_thr=mcond((u[0]-x[3])/(0.001+abs(x[3]))>0.0,0.0,1,0.0)*daux_diff_K__thr;
    if (isNANorINF(daux_diff_2K_thr_thr))
        {
            PRNT("    @k %d: daux_diff_2K_thr_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2K_thr_thr);
            {
                PRNT("        daux_diff_K__thr= %g\n",daux_diff_K__thr);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
// aux_deriv[108]
    daux_diff_alpha_Ux=-u[1]/(0.001+abs(x[3]))+u[1]*x[3]*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)-(x[4]+p[4][0]*x[5])*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0));
    if (isNANorINF(daux_diff_alpha_Ux))
        {
            PRNT("    @k %d: daux_diff_alpha_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_alpha_Ux);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[109]
    daux_diff_2alpha_Ux_Ux=u[1]*x[3]*mcond(x[3]>0.0,0.0,1,0.0)/pow(0.001+abs(x[3]),2.0)-(x[4]+p[4][0]*x[5])*mcond(x[3]>0.0,0.0,1,0.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))+2.0*u[1]*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)-2.0*u[1]*x[3]*pow(mcond(x[3]>0.0,1.0,1,-1.0),2.0)/pow(0.001+abs(x[3]),3.0)-2.0*pow(x[4]+p[4][0]*x[5],3.0)*pow(mcond(x[3]>0.0,1.0,1,-1.0),2.0)/pow(0.001+abs(x[3]),5.0)/pow(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0),2.0)+2.0*(x[4]+p[4][0]*x[5])*pow(mcond(x[3]>0.0,1.0,1,-1.0),2.0)/pow(0.001+abs(x[3]),3.0)/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0));
    if (isNANorINF(daux_diff_2alpha_Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2alpha_Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_Ux_Ux);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[110]
    daux_diff_2alpha_Ux_Uy=2.0*pow(x[4]+p[4][0]*x[5],2.0)*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),4.0)/pow(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0),2.0)-mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0));
    if (isNANorINF(daux_diff_2alpha_Ux_Uy))
        {
            PRNT("    @k %d: daux_diff_2alpha_Ux_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_Ux_Uy);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[111]
    daux_diff_2alpha_Ux_r=2.0*p[4][0]*pow(x[4]+p[4][0]*x[5],2.0)*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),4.0)/pow(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0),2.0)-p[4][0]*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0));
    if (isNANorINF(daux_diff_2alpha_Ux_r))
        {
            PRNT("    @k %d: daux_diff_2alpha_Ux_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_Ux_r);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[112]
    daux_diff_2alpha_Ux_steer=-1/(0.001+abs(x[3]))+x[3]*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0);
    if (isNANorINF(daux_diff_2alpha_Ux_steer))
        {
            PRNT("    @k %d: daux_diff_2alpha_Ux_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_Ux_steer);
            {
                PRNT("        x[3]= %g\n",x[3]);
            }
            return(0);
        }
// aux_deriv[113]
    daux_diff_alpha_Uy=1/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0));
    if (isNANorINF(daux_diff_alpha_Uy))
        {
            PRNT("    @k %d: daux_diff_alpha_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_alpha_Uy);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[114]
    daux_diff_2alpha_Uy_Uy=-2.0*(x[4]+p[4][0]*x[5])/pow(0.001+abs(x[3]),3.0)/pow(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0),2.0);
    if (isNANorINF(daux_diff_2alpha_Uy_Uy))
        {
            PRNT("    @k %d: daux_diff_2alpha_Uy_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_Uy_Uy);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[115]
    daux_diff_2alpha_Uy_r=-2.0*p[4][0]*(x[4]+p[4][0]*x[5])/pow(0.001+abs(x[3]),3.0)/pow(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0),2.0);
    if (isNANorINF(daux_diff_2alpha_Uy_r))
        {
            PRNT("    @k %d: daux_diff_2alpha_Uy_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_Uy_r);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[116]
    daux_diff_alpha_r=p[4][0]/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0));
    if (isNANorINF(daux_diff_alpha_r))
        {
            PRNT("    @k %d: daux_diff_alpha_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_alpha_r);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[117]
    daux_diff_2alpha_r_r=-2.0*pow(p[4][0],2.0)*(x[4]+p[4][0]*x[5])/pow(0.001+abs(x[3]),3.0)/pow(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0),2.0);
    if (isNANorINF(daux_diff_2alpha_r_r))
        {
            PRNT("    @k %d: daux_diff_2alpha_r_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_r_r);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[118]
    daux_diff_alpha_steer=-x[3]/(0.001+abs(x[3]));
    if (isNANorINF(daux_diff_alpha_steer))
        {
            PRNT("    @k %d: daux_diff_alpha_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_alpha_steer);
            {
                PRNT("        x[3]= %g\n",x[3]);
            }
            return(0);
        }
// aux_deriv[119]
    daux_diff_alpha_F_Ux=mcond(abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))>1.5707963267948966,(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*daux_diff_alpha_Ux/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))),1,daux_diff_alpha_Ux);
    if (isNANorINF(daux_diff_alpha_F_Ux))
        {
            PRNT("    @k %d: daux_diff_alpha_F_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_alpha_F_Ux);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_alpha_Ux= %g\n",daux_diff_alpha_Ux);
            }
            return(0);
        }
// aux_deriv[120]
    daux_diff_2alpha_F_Ux_Ux=mcond(abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))>1.5707963267948966,(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*daux_diff_2alpha_Ux_Ux/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_2alpha_Ux_Ux/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_2alpha_Ux_Ux/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(u[1]/(0.001+abs(x[3]))-u[1]*x[3]*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)+(x[4]+p[4][0]*x[5])*
     mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0)))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(u[1]/(0.001+abs(x[3]))-u[1]*x[3]*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)+(x[4]+p[4][0]*x[5])*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0)))*(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,0.0,1,0.0)*
     daux_diff_alpha_Ux/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,0.0,1,0.0)*daux_diff_alpha_Ux/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(-u[1]/(0.001+abs(x[3]))+u[1]*x[3]*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)-(x[4]+p[4][0]*x[5])*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0)))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+(-u[1]/(0.001+abs(x[3]))+u[1]*x[3]*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)-(x[4]+p[4][0]*x[5])*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+
     abs(x[3]),2.0)))*(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/pow(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)+(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*(u[1]/(0.001+abs(x[3]))-u[1]*x[3]*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)+(x[4]+p[4][0]*x[5])*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0)))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)+(u[1]/(0.001+abs(x[3]))-u[1]*x[3]*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)+(x[4]+p[4][0]*x[5
     ])*mcond(x[3]>0.0,1.0,1,-1.0)/pow(0.001+abs(x[3]),2.0)/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0)))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))),1,daux_diff_2alpha_Ux_Ux);
    if (isNANorINF(daux_diff_2alpha_F_Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2alpha_F_Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_F_Ux_Ux);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_2alpha_Ux_Ux= %g\n",daux_diff_2alpha_Ux_Ux);
                PRNT("        daux_diff_alpha_Ux= %g\n",daux_diff_alpha_Ux);
            }
            return(0);
        }
// aux_deriv[121]
    daux_diff_2alpha_F_Ux_Uy=mcond(abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))>1.5707963267948966,(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*daux_diff_2alpha_Ux_Uy/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_2alpha_Ux_Uy/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_2alpha_Ux_Uy/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*
     daux_diff_alpha_Ux/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,0.0,1,0.0)*daux_diff_alpha_Ux/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,0.0,1,
     0.0)*daux_diff_alpha_Ux/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/pow(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,
     -1.0)*daux_diff_alpha_Ux/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))),1,daux_diff_2alpha_Ux_Uy);
    if (isNANorINF(daux_diff_2alpha_F_Ux_Uy))
        {
            PRNT("    @k %d: daux_diff_2alpha_F_Ux_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_F_Ux_Uy);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_2alpha_Ux_Uy= %g\n",daux_diff_2alpha_Ux_Uy);
                PRNT("        daux_diff_alpha_Ux= %g\n",daux_diff_alpha_Ux);
            }
            return(0);
        }
// aux_deriv[122]
    daux_diff_2alpha_F_Ux_r=mcond(abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))>1.5707963267948966,(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*daux_diff_2alpha_Ux_r/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_2alpha_Ux_r/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_2alpha_Ux_r/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+p[4][0]*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*
     daux_diff_alpha_Ux/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+p[4][0]*(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,0.0,1,0.0)*daux_diff_alpha_Ux/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0
     ,0.0,1,0.0)*daux_diff_alpha_Ux/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-p[4][0]*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+p[4][0]*(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/pow(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-p[4][0]*(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5]
     )/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-p[4][0]*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))),1,daux_diff_2alpha_Ux_r);
    if (isNANorINF(daux_diff_2alpha_F_Ux_r))
        {
            PRNT("    @k %d: daux_diff_2alpha_F_Ux_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_F_Ux_r);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_2alpha_Ux_r= %g\n",daux_diff_2alpha_Ux_r);
                PRNT("        daux_diff_alpha_Ux= %g\n",daux_diff_alpha_Ux);
            }
            return(0);
        }
// aux_deriv[123]
    daux_diff_2alpha_F_Ux_steer=mcond(abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))>1.5707963267948966,(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*daux_diff_2alpha_Ux_steer/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_2alpha_Ux_steer/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_2alpha_Ux_steer/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-x[3]*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,
     -1.0)*daux_diff_alpha_Ux/(0.001+abs(x[3]))/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-x[3]*(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/(0.001+abs(x[3]))/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,0.0,1,0.0)*daux_diff_alpha_Ux/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,0.0,1,0.0)*daux_diff_alpha_Ux/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+x[3]*
     mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/(0.001+abs(x[3]))/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-x[3]*(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/(0.001+abs(x[3]))/pow(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)+x[3]*(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/(0.001+abs(x[3]))/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)+x[3]*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))
     )>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Ux/(0.001+abs(x[3]))/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))),1,daux_diff_2alpha_Ux_steer);
    if (isNANorINF(daux_diff_2alpha_F_Ux_steer))
        {
            PRNT("    @k %d: daux_diff_2alpha_F_Ux_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_F_Ux_steer);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_2alpha_Ux_steer= %g\n",daux_diff_2alpha_Ux_steer);
                PRNT("        daux_diff_alpha_Ux= %g\n",daux_diff_alpha_Ux);
            }
            return(0);
        }
// aux_deriv[124]
    daux_diff_alpha_F_Uy=mcond(abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))>1.5707963267948966,(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*daux_diff_alpha_Uy/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Uy/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Uy/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))),1,daux_diff_alpha_Uy);
    if (isNANorINF(daux_diff_alpha_F_Uy))
        {
            PRNT("    @k %d: daux_diff_alpha_F_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_alpha_F_Uy);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_alpha_Uy= %g\n",daux_diff_alpha_Uy);
            }
            return(0);
        }
// aux_deriv[125]
    daux_diff_2alpha_F_Uy_Uy=mcond(abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))>1.5707963267948966,(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*daux_diff_2alpha_Uy_Uy/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_2alpha_Uy_Uy/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_2alpha_Uy_Uy/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*
     daux_diff_alpha_Uy/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Uy/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,0.0,1,0.0)*daux_diff_alpha_Uy/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,0.0,1,
     0.0)*daux_diff_alpha_Uy/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Uy/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Uy/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/pow(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,
     -1.0)*daux_diff_alpha_Uy/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Uy/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))),1,daux_diff_2alpha_Uy_Uy);
    if (isNANorINF(daux_diff_2alpha_F_Uy_Uy))
        {
            PRNT("    @k %d: daux_diff_2alpha_F_Uy_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_F_Uy_Uy);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_2alpha_Uy_Uy= %g\n",daux_diff_2alpha_Uy_Uy);
                PRNT("        daux_diff_alpha_Uy= %g\n",daux_diff_alpha_Uy);
            }
            return(0);
        }
// aux_deriv[126]
    daux_diff_2alpha_F_Uy_r=mcond(abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))>1.5707963267948966,(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*daux_diff_2alpha_Uy_r/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_2alpha_Uy_r/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_2alpha_Uy_r/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+p[4][0]*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*
     daux_diff_alpha_Uy/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+p[4][0]*(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Uy/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,0.0,1,0.0)*daux_diff_alpha_Uy/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0
     ,0.0,1,0.0)*daux_diff_alpha_Uy/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-p[4][0]*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Uy/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+p[4][0]*(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Uy/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/pow(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-p[4][0]*(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5]
     )/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Uy/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-p[4][0]*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Uy/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))),1,daux_diff_2alpha_Uy_r);
    if (isNANorINF(daux_diff_2alpha_F_Uy_r))
        {
            PRNT("    @k %d: daux_diff_2alpha_F_Uy_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_F_Uy_r);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_2alpha_Uy_r= %g\n",daux_diff_2alpha_Uy_r);
                PRNT("        daux_diff_alpha_Uy= %g\n",daux_diff_alpha_Uy);
            }
            return(0);
        }
// aux_deriv[127]
    daux_diff_2alpha_F_Uy_steer=mcond(abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))>1.5707963267948966,-x[3]*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Uy/(0.001+abs(x[3]))/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-x[3]*(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Uy/(0.001+abs(x[3]))/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,0.0,1,0.0)*daux_diff_alpha_Uy/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0
     ]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,0.0,1,0.0)*daux_diff_alpha_Uy/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+x[3]*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Uy/(0.001+abs(x[3]))/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-x[3]*(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Uy/(0.001+abs(x[3]))/pow(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)+x[3]*(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x
     [3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Uy/(0.001+abs(x[3]))/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)+x[3]*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_Uy/(0.001+abs(x[3]))/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))),1,0.0);
    if (isNANorINF(daux_diff_2alpha_F_Uy_steer))
        {
            PRNT("    @k %d: daux_diff_2alpha_F_Uy_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_F_Uy_steer);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_alpha_Uy= %g\n",daux_diff_alpha_Uy);
            }
            return(0);
        }
// aux_deriv[128]
    daux_diff_alpha_F_r=mcond(abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))>1.5707963267948966,(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*daux_diff_alpha_r/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_r/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_r/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))),1,daux_diff_alpha_r);
    if (isNANorINF(daux_diff_alpha_F_r))
        {
            PRNT("    @k %d: daux_diff_alpha_F_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_alpha_F_r);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_alpha_r= %g\n",daux_diff_alpha_r);
            }
            return(0);
        }
// aux_deriv[129]
    daux_diff_2alpha_F_r_r=mcond(abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))>1.5707963267948966,(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*daux_diff_2alpha_r_r/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_2alpha_r_r/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_2alpha_r_r/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+p[4][0]*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*
     daux_diff_alpha_r/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+p[4][0]*(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_r/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,0.0,1,0.0)*daux_diff_alpha_r/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,
     0.0,1,0.0)*daux_diff_alpha_r/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-p[4][0]*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_r/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+p[4][0]*(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_r/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/pow(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-p[4][0]*(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(
     0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_r/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-p[4][0]*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_r/(0.001+abs(x[3]))/(1.0+pow(x[4]+p[4][0]*x[5],2.0)/pow(0.001+abs(x[3]),2.0))/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))),1,daux_diff_2alpha_r_r);
    if (isNANorINF(daux_diff_2alpha_F_r_r))
        {
            PRNT("    @k %d: daux_diff_2alpha_F_r_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_F_r_r);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_2alpha_r_r= %g\n",daux_diff_2alpha_r_r);
                PRNT("        daux_diff_alpha_r= %g\n",daux_diff_alpha_r);
            }
            return(0);
        }
// aux_deriv[130]
    daux_diff_2alpha_F_r_steer=mcond(abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))>1.5707963267948966,-x[3]*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_r/(0.001+abs(x[3]))/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-x[3]*(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_r/(0.001+abs(x[3]))/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,0.0,1,0.0)*daux_diff_alpha_r/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[
     5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,0.0,1,0.0)*daux_diff_alpha_r/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+x[3]*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_r/(0.001+abs(x[3]))/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-x[3]*(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_r/(0.001+abs(x[3]))/pow(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)+x[3]*(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>
     0.0,1.0,1,-1.0)*daux_diff_alpha_r/(0.001+abs(x[3]))/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)+x[3]*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_r/(0.001+abs(x[3]))/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))),1,0.0);
    if (isNANorINF(daux_diff_2alpha_F_r_steer))
        {
            PRNT("    @k %d: daux_diff_2alpha_F_r_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_F_r_steer);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_alpha_r= %g\n",daux_diff_alpha_r);
            }
            return(0);
        }
// aux_deriv[131]
    daux_diff_alpha_F_steer=mcond(abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))>1.5707963267948966,(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*daux_diff_alpha_steer/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_steer/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_steer/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))),1,daux_diff_alpha_steer);
    if (isNANorINF(daux_diff_alpha_F_steer))
        {
            PRNT("    @k %d: daux_diff_alpha_F_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_alpha_F_steer);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_alpha_steer= %g\n",daux_diff_alpha_steer);
            }
            return(0);
        }
// aux_deriv[132]
    daux_diff_2alpha_F_steer_steer=mcond(abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))>1.5707963267948966,-x[3]*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_steer/(0.001+abs(x[3]))/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-x[3]*(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_steer/(0.001+abs(x[3]))/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)-(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,0.0,1,0.0)*daux_diff_alpha_steer/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan(
     (x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,0.0,1,0.0)*daux_diff_alpha_steer/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))+x[3]*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_steer/(0.001+abs(x[3]))/abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))-x[3]*(3.141592653589793-abs(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))))*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_steer/(0.001+abs(x[3]))/pow(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)+x[3]*(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))))*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][
     0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_steer/(0.001+abs(x[3]))/pow(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3]))),2.0)+x[3]*mcond(u[1]*x[3]/(0.001+abs(x[3]))-atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*mcond(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))>0.0,1.0,1,-1.0)*daux_diff_alpha_steer/(0.001+abs(x[3]))/(-u[1]*x[3]/(0.001+abs(x[3]))+atan((x[4]+p[4][0]*x[5])/(0.001+abs(x[3])))),1,0.0);
    if (isNANorINF(daux_diff_2alpha_F_steer_steer))
        {
            PRNT("    @k %d: daux_diff_2alpha_F_steer_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2alpha_F_steer_steer);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[4,0]= %g\n",p[4][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_alpha_steer= %g\n",daux_diff_alpha_steer);
            }
            return(0);
        }
// aux_deriv[133]
    daux_diff_Fyf_Ux=-p[6][0]*derivative(aux_F_f,x[3],1.0)*tan(aux_alpha_F)/aux_gamma_F+p[6][0]*aux_F_f*tan(aux_alpha_F)*derivative(aux_gamma_F,x[3],1.0)/pow(aux_gamma_F,2.0)-p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*daux_diff_alpha_F_Ux/aux_gamma_F;
    if (isNANorINF(daux_diff_Fyf_Ux))
        {
            PRNT("    @k %d: daux_diff_Fyf_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_Fyf_Ux);
            {
                PRNT("        aux_F_f= %g\n",aux_F_f);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        aux_alpha_F= %g\n",aux_alpha_F);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        aux_gamma_F= %g\n",aux_gamma_F);
                PRNT("        daux_diff_alpha_F_Ux= %g\n",daux_diff_alpha_F_Ux);
            }
            return(0);
        }
// aux_deriv[134]
    daux_diff_2Fyf_Ux_Ux=-p[6][0]*derivative(aux_F_f,x[3],2.0)*tan(aux_alpha_F)/aux_gamma_F+2.0*p[6][0]*derivative(aux_F_f,x[3],1.0)*tan(aux_alpha_F)*derivative(aux_gamma_F,x[3],1.0)/pow(aux_gamma_F,2.0)-2.0*p[6][0]*aux_F_f*tan(aux_alpha_F)*pow(derivative(aux_gamma_F,x[3],1.0),2.0)/pow(aux_gamma_F,3.0)+p[6][0]*aux_F_f*tan(aux_alpha_F)*derivative(aux_gamma_F,x[3],2.0)/pow(aux_gamma_F,2.0)-p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*daux_diff_2alpha_F_Ux_Ux/aux_gamma_F-2.0*p[6][0]*derivative(aux_F_f,x[3],1.0)*pow(sec(aux_alpha_F),2.0)*daux_diff_alpha_F_Ux/aux_gamma_F+2.0*p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*derivative(aux_gamma_F,x[3],1.0)*daux_diff_alpha_F_Ux/pow(aux_gamma_F,2.0)-2.0*p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*tan(aux_alpha_F)*pow(daux_diff_alpha_F_Ux,2.0)/aux_gamma_F;
    if (isNANorINF(daux_diff_2Fyf_Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2Fyf_Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyf_Ux_Ux);
            {
                PRNT("        aux_F_f= %g\n",aux_F_f);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        aux_alpha_F= %g\n",aux_alpha_F);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        aux_gamma_F= %g\n",aux_gamma_F);
                PRNT("        daux_diff_2alpha_F_Ux_Ux= %g\n",daux_diff_2alpha_F_Ux_Ux);
                PRNT("        daux_diff_alpha_F_Ux= %g\n",daux_diff_alpha_F_Ux);
            }
            return(0);
        }
// aux_deriv[135]
    daux_diff_2Fyf_Ux_Uy=-p[6][0]*derivative(aux_F_f,x[3],1.0,x[4],1.0)*tan(aux_alpha_F)/aux_gamma_F+p[6][0]*derivative(aux_F_f,x[4],1.0)*tan(aux_alpha_F)*derivative(aux_gamma_F,x[3],1.0)/pow(aux_gamma_F,2.0)+p[6][0]*aux_F_f*tan(aux_alpha_F)*derivative(aux_gamma_F,x[3],1.0,x[4],1.0)/pow(aux_gamma_F,2.0)+p[6][0]*derivative(aux_F_f,x[3],1.0)*tan(aux_alpha_F)*derivative(aux_gamma_F,x[4],1.0)/pow(aux_gamma_F,2.0)-2.0*p[6][0]*aux_F_f*tan(aux_alpha_F)*derivative(aux_gamma_F,x[3],1.0)*derivative(aux_gamma_F,x[4],1.0)/pow(aux_gamma_F,3.0)-p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*daux_diff_2alpha_F_Ux_Uy/aux_gamma_F-p[6][0]*derivative(aux_F_f,x[4],1.0)*pow(sec(aux_alpha_F),2.0)*daux_diff_alpha_F_Ux/aux_gamma_F+p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*derivative(aux_gamma_F,x[4],1.0)*daux_diff_alpha_F_Ux/pow(aux_gamma_F,2.0)-p[6][0]*derivative(aux_F_f,x[3],1.0)*pow(sec(aux_alpha_F),2.0)*daux_diff_alpha_F_Uy/aux_gamma_F+p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*derivative(aux_gamma_F,x[3],1.0)*
     daux_diff_alpha_F_Uy/pow(aux_gamma_F,2.0)-2.0*p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*tan(aux_alpha_F)*daux_diff_alpha_F_Ux*daux_diff_alpha_F_Uy/aux_gamma_F;
    if (isNANorINF(daux_diff_2Fyf_Ux_Uy))
        {
            PRNT("    @k %d: daux_diff_2Fyf_Ux_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyf_Ux_Uy);
            {
                PRNT("        aux_F_f= %g\n",aux_F_f);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        aux_alpha_F= %g\n",aux_alpha_F);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        aux_gamma_F= %g\n",aux_gamma_F);
                PRNT("        daux_diff_2alpha_F_Ux_Uy= %g\n",daux_diff_2alpha_F_Ux_Uy);
                PRNT("        daux_diff_alpha_F_Ux= %g\n",daux_diff_alpha_F_Ux);
                PRNT("        daux_diff_alpha_F_Uy= %g\n",daux_diff_alpha_F_Uy);
            }
            return(0);
        }
// aux_deriv[136]
    daux_diff_2Fyf_Ux_r=-p[6][0]*derivative(aux_F_f,x[3],1.0,x[5],1.0)*tan(aux_alpha_F)/aux_gamma_F+p[6][0]*derivative(aux_F_f,x[5],1.0)*tan(aux_alpha_F)*derivative(aux_gamma_F,x[3],1.0)/pow(aux_gamma_F,2.0)+p[6][0]*aux_F_f*tan(aux_alpha_F)*derivative(aux_gamma_F,x[3],1.0,x[5],1.0)/pow(aux_gamma_F,2.0)+p[6][0]*derivative(aux_F_f,x[3],1.0)*tan(aux_alpha_F)*derivative(aux_gamma_F,x[5],1.0)/pow(aux_gamma_F,2.0)-2.0*p[6][0]*aux_F_f*tan(aux_alpha_F)*derivative(aux_gamma_F,x[3],1.0)*derivative(aux_gamma_F,x[5],1.0)/pow(aux_gamma_F,3.0)-p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*daux_diff_2alpha_F_Ux_r/aux_gamma_F-p[6][0]*derivative(aux_F_f,x[5],1.0)*pow(sec(aux_alpha_F),2.0)*daux_diff_alpha_F_Ux/aux_gamma_F+p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*derivative(aux_gamma_F,x[5],1.0)*daux_diff_alpha_F_Ux/pow(aux_gamma_F,2.0)-p[6][0]*derivative(aux_F_f,x[3],1.0)*pow(sec(aux_alpha_F),2.0)*daux_diff_alpha_F_r/aux_gamma_F+p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*derivative(aux_gamma_F,x[3],1.0)*
     daux_diff_alpha_F_r/pow(aux_gamma_F,2.0)-2.0*p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*tan(aux_alpha_F)*daux_diff_alpha_F_Ux*daux_diff_alpha_F_r/aux_gamma_F;
    if (isNANorINF(daux_diff_2Fyf_Ux_r))
        {
            PRNT("    @k %d: daux_diff_2Fyf_Ux_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyf_Ux_r);
            {
                PRNT("        aux_F_f= %g\n",aux_F_f);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        aux_alpha_F= %g\n",aux_alpha_F);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        aux_gamma_F= %g\n",aux_gamma_F);
                PRNT("        daux_diff_2alpha_F_Ux_r= %g\n",daux_diff_2alpha_F_Ux_r);
                PRNT("        daux_diff_alpha_F_Ux= %g\n",daux_diff_alpha_F_Ux);
                PRNT("        daux_diff_alpha_F_r= %g\n",daux_diff_alpha_F_r);
            }
            return(0);
        }
// aux_deriv[137]
    daux_diff_2Fyf_Ux_steer=-p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*daux_diff_2alpha_F_Ux_steer/aux_gamma_F-p[6][0]*derivative(aux_F_f,x[3],1.0)*pow(sec(aux_alpha_F),2.0)*daux_diff_alpha_F_steer/aux_gamma_F+p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*derivative(aux_gamma_F,x[3],1.0)*daux_diff_alpha_F_steer/pow(aux_gamma_F,2.0)-2.0*p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*tan(aux_alpha_F)*daux_diff_alpha_F_Ux*daux_diff_alpha_F_steer/aux_gamma_F;
    if (isNANorINF(daux_diff_2Fyf_Ux_steer))
        {
            PRNT("    @k %d: daux_diff_2Fyf_Ux_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyf_Ux_steer);
            {
                PRNT("        aux_F_f= %g\n",aux_F_f);
                PRNT("        aux_alpha_F= %g\n",aux_alpha_F);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_2alpha_F_Ux_steer= %g\n",daux_diff_2alpha_F_Ux_steer);
                PRNT("        aux_gamma_F= %g\n",aux_gamma_F);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        daux_diff_alpha_F_steer= %g\n",daux_diff_alpha_F_steer);
                PRNT("        daux_diff_alpha_F_Ux= %g\n",daux_diff_alpha_F_Ux);
            }
            return(0);
        }
// aux_deriv[138]
    daux_diff_Fyf_Uy=-p[6][0]*derivative(aux_F_f,x[4],1.0)*tan(aux_alpha_F)/aux_gamma_F+p[6][0]*aux_F_f*tan(aux_alpha_F)*derivative(aux_gamma_F,x[4],1.0)/pow(aux_gamma_F,2.0)-p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*daux_diff_alpha_F_Uy/aux_gamma_F;
    if (isNANorINF(daux_diff_Fyf_Uy))
        {
            PRNT("    @k %d: daux_diff_Fyf_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_Fyf_Uy);
            {
                PRNT("        aux_F_f= %g\n",aux_F_f);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        aux_alpha_F= %g\n",aux_alpha_F);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        aux_gamma_F= %g\n",aux_gamma_F);
                PRNT("        daux_diff_alpha_F_Uy= %g\n",daux_diff_alpha_F_Uy);
            }
            return(0);
        }
// aux_deriv[139]
    daux_diff_2Fyf_Uy_Uy=-p[6][0]*derivative(aux_F_f,x[4],2.0)*tan(aux_alpha_F)/aux_gamma_F+2.0*p[6][0]*derivative(aux_F_f,x[4],1.0)*tan(aux_alpha_F)*derivative(aux_gamma_F,x[4],1.0)/pow(aux_gamma_F,2.0)-2.0*p[6][0]*aux_F_f*tan(aux_alpha_F)*pow(derivative(aux_gamma_F,x[4],1.0),2.0)/pow(aux_gamma_F,3.0)+p[6][0]*aux_F_f*tan(aux_alpha_F)*derivative(aux_gamma_F,x[4],2.0)/pow(aux_gamma_F,2.0)-p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*daux_diff_2alpha_F_Uy_Uy/aux_gamma_F-2.0*p[6][0]*derivative(aux_F_f,x[4],1.0)*pow(sec(aux_alpha_F),2.0)*daux_diff_alpha_F_Uy/aux_gamma_F+2.0*p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*derivative(aux_gamma_F,x[4],1.0)*daux_diff_alpha_F_Uy/pow(aux_gamma_F,2.0)-2.0*p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*tan(aux_alpha_F)*pow(daux_diff_alpha_F_Uy,2.0)/aux_gamma_F;
    if (isNANorINF(daux_diff_2Fyf_Uy_Uy))
        {
            PRNT("    @k %d: daux_diff_2Fyf_Uy_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyf_Uy_Uy);
            {
                PRNT("        aux_F_f= %g\n",aux_F_f);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        aux_alpha_F= %g\n",aux_alpha_F);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        aux_gamma_F= %g\n",aux_gamma_F);
                PRNT("        daux_diff_2alpha_F_Uy_Uy= %g\n",daux_diff_2alpha_F_Uy_Uy);
                PRNT("        daux_diff_alpha_F_Uy= %g\n",daux_diff_alpha_F_Uy);
            }
            return(0);
        }
// aux_deriv[140]
    daux_diff_2Fyf_Uy_r=-p[6][0]*derivative(aux_F_f,x[4],1.0,x[5],1.0)*tan(aux_alpha_F)/aux_gamma_F+p[6][0]*derivative(aux_F_f,x[5],1.0)*tan(aux_alpha_F)*derivative(aux_gamma_F,x[4],1.0)/pow(aux_gamma_F,2.0)+p[6][0]*aux_F_f*tan(aux_alpha_F)*derivative(aux_gamma_F,x[4],1.0,x[5],1.0)/pow(aux_gamma_F,2.0)+p[6][0]*derivative(aux_F_f,x[4],1.0)*tan(aux_alpha_F)*derivative(aux_gamma_F,x[5],1.0)/pow(aux_gamma_F,2.0)-2.0*p[6][0]*aux_F_f*tan(aux_alpha_F)*derivative(aux_gamma_F,x[4],1.0)*derivative(aux_gamma_F,x[5],1.0)/pow(aux_gamma_F,3.0)-p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*daux_diff_2alpha_F_Uy_r/aux_gamma_F-p[6][0]*derivative(aux_F_f,x[5],1.0)*pow(sec(aux_alpha_F),2.0)*daux_diff_alpha_F_Uy/aux_gamma_F+p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*derivative(aux_gamma_F,x[5],1.0)*daux_diff_alpha_F_Uy/pow(aux_gamma_F,2.0)-p[6][0]*derivative(aux_F_f,x[4],1.0)*pow(sec(aux_alpha_F),2.0)*daux_diff_alpha_F_r/aux_gamma_F+p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*derivative(aux_gamma_F,x[4],1.0)*
     daux_diff_alpha_F_r/pow(aux_gamma_F,2.0)-2.0*p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*tan(aux_alpha_F)*daux_diff_alpha_F_Uy*daux_diff_alpha_F_r/aux_gamma_F;
    if (isNANorINF(daux_diff_2Fyf_Uy_r))
        {
            PRNT("    @k %d: daux_diff_2Fyf_Uy_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyf_Uy_r);
            {
                PRNT("        aux_F_f= %g\n",aux_F_f);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        aux_alpha_F= %g\n",aux_alpha_F);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        aux_gamma_F= %g\n",aux_gamma_F);
                PRNT("        daux_diff_2alpha_F_Uy_r= %g\n",daux_diff_2alpha_F_Uy_r);
                PRNT("        daux_diff_alpha_F_Uy= %g\n",daux_diff_alpha_F_Uy);
                PRNT("        daux_diff_alpha_F_r= %g\n",daux_diff_alpha_F_r);
            }
            return(0);
        }
// aux_deriv[141]
    daux_diff_2Fyf_Uy_steer=-p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*daux_diff_2alpha_F_Uy_steer/aux_gamma_F-p[6][0]*derivative(aux_F_f,x[4],1.0)*pow(sec(aux_alpha_F),2.0)*daux_diff_alpha_F_steer/aux_gamma_F+p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*derivative(aux_gamma_F,x[4],1.0)*daux_diff_alpha_F_steer/pow(aux_gamma_F,2.0)-2.0*p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*tan(aux_alpha_F)*daux_diff_alpha_F_Uy*daux_diff_alpha_F_steer/aux_gamma_F;
    if (isNANorINF(daux_diff_2Fyf_Uy_steer))
        {
            PRNT("    @k %d: daux_diff_2Fyf_Uy_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyf_Uy_steer);
            {
                PRNT("        aux_F_f= %g\n",aux_F_f);
                PRNT("        aux_alpha_F= %g\n",aux_alpha_F);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_2alpha_F_Uy_steer= %g\n",daux_diff_2alpha_F_Uy_steer);
                PRNT("        aux_gamma_F= %g\n",aux_gamma_F);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        daux_diff_alpha_F_steer= %g\n",daux_diff_alpha_F_steer);
                PRNT("        daux_diff_alpha_F_Uy= %g\n",daux_diff_alpha_F_Uy);
            }
            return(0);
        }
// aux_deriv[142]
    daux_diff_Fyf_r=-p[6][0]*derivative(aux_F_f,x[5],1.0)*tan(aux_alpha_F)/aux_gamma_F+p[6][0]*aux_F_f*tan(aux_alpha_F)*derivative(aux_gamma_F,x[5],1.0)/pow(aux_gamma_F,2.0)-p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*daux_diff_alpha_F_r/aux_gamma_F;
    if (isNANorINF(daux_diff_Fyf_r))
        {
            PRNT("    @k %d: daux_diff_Fyf_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_Fyf_r);
            {
                PRNT("        aux_F_f= %g\n",aux_F_f);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        aux_alpha_F= %g\n",aux_alpha_F);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        aux_gamma_F= %g\n",aux_gamma_F);
                PRNT("        daux_diff_alpha_F_r= %g\n",daux_diff_alpha_F_r);
            }
            return(0);
        }
// aux_deriv[143]
    daux_diff_2Fyf_r_r=-p[6][0]*derivative(aux_F_f,x[5],2.0)*tan(aux_alpha_F)/aux_gamma_F+2.0*p[6][0]*derivative(aux_F_f,x[5],1.0)*tan(aux_alpha_F)*derivative(aux_gamma_F,x[5],1.0)/pow(aux_gamma_F,2.0)-2.0*p[6][0]*aux_F_f*tan(aux_alpha_F)*pow(derivative(aux_gamma_F,x[5],1.0),2.0)/pow(aux_gamma_F,3.0)+p[6][0]*aux_F_f*tan(aux_alpha_F)*derivative(aux_gamma_F,x[5],2.0)/pow(aux_gamma_F,2.0)-p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*daux_diff_2alpha_F_r_r/aux_gamma_F-2.0*p[6][0]*derivative(aux_F_f,x[5],1.0)*pow(sec(aux_alpha_F),2.0)*daux_diff_alpha_F_r/aux_gamma_F+2.0*p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*derivative(aux_gamma_F,x[5],1.0)*daux_diff_alpha_F_r/pow(aux_gamma_F,2.0)-2.0*p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*tan(aux_alpha_F)*pow(daux_diff_alpha_F_r,2.0)/aux_gamma_F;
    if (isNANorINF(daux_diff_2Fyf_r_r))
        {
            PRNT("    @k %d: daux_diff_2Fyf_r_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyf_r_r);
            {
                PRNT("        aux_F_f= %g\n",aux_F_f);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        aux_alpha_F= %g\n",aux_alpha_F);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        aux_gamma_F= %g\n",aux_gamma_F);
                PRNT("        daux_diff_2alpha_F_r_r= %g\n",daux_diff_2alpha_F_r_r);
                PRNT("        daux_diff_alpha_F_r= %g\n",daux_diff_alpha_F_r);
            }
            return(0);
        }
// aux_deriv[144]
    daux_diff_2Fyf_r_steer=-p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*daux_diff_2alpha_F_r_steer/aux_gamma_F-p[6][0]*derivative(aux_F_f,x[5],1.0)*pow(sec(aux_alpha_F),2.0)*daux_diff_alpha_F_steer/aux_gamma_F+p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*derivative(aux_gamma_F,x[5],1.0)*daux_diff_alpha_F_steer/pow(aux_gamma_F,2.0)-2.0*p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*tan(aux_alpha_F)*daux_diff_alpha_F_r*daux_diff_alpha_F_steer/aux_gamma_F;
    if (isNANorINF(daux_diff_2Fyf_r_steer))
        {
            PRNT("    @k %d: daux_diff_2Fyf_r_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyf_r_steer);
            {
                PRNT("        aux_F_f= %g\n",aux_F_f);
                PRNT("        aux_alpha_F= %g\n",aux_alpha_F);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_2alpha_F_r_steer= %g\n",daux_diff_2alpha_F_r_steer);
                PRNT("        aux_gamma_F= %g\n",aux_gamma_F);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        daux_diff_alpha_F_steer= %g\n",daux_diff_alpha_F_steer);
                PRNT("        daux_diff_alpha_F_r= %g\n",daux_diff_alpha_F_r);
            }
            return(0);
        }
// aux_deriv[145]
    daux_diff_Fyf_steer=-p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*daux_diff_alpha_F_steer/aux_gamma_F;
    if (isNANorINF(daux_diff_Fyf_steer))
        {
            PRNT("    @k %d: daux_diff_Fyf_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_Fyf_steer);
            {
                PRNT("        aux_F_f= %g\n",aux_F_f);
                PRNT("        aux_alpha_F= %g\n",aux_alpha_F);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_alpha_F_steer= %g\n",daux_diff_alpha_F_steer);
                PRNT("        aux_gamma_F= %g\n",aux_gamma_F);
            }
            return(0);
        }
// aux_deriv[146]
    daux_diff_2Fyf_steer_steer=-p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*daux_diff_2alpha_F_steer_steer/aux_gamma_F-2.0*p[6][0]*aux_F_f*pow(sec(aux_alpha_F),2.0)*tan(aux_alpha_F)*pow(daux_diff_alpha_F_steer,2.0)/aux_gamma_F;
    if (isNANorINF(daux_diff_2Fyf_steer_steer))
        {
            PRNT("    @k %d: daux_diff_2Fyf_steer_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fyf_steer_steer);
            {
                PRNT("        aux_F_f= %g\n",aux_F_f);
                PRNT("        aux_alpha_F= %g\n",aux_alpha_F);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_2alpha_F_steer_steer= %g\n",daux_diff_2alpha_F_steer_steer);
                PRNT("        aux_gamma_F= %g\n",aux_gamma_F);
                PRNT("        daux_diff_alpha_F_steer= %g\n",daux_diff_alpha_F_steer);
            }
            return(0);
        }
// aux_deriv[147]
    daux_diff_gamma_F_Ux=abs(p[6][0])*pow(sec(aux_alpha_R),2.0)*mcond(tan(aux_alpha_R)>0.0,1.0,1,-1.0)*daux_diff_alpha_R_Ux;
    if (isNANorINF(daux_diff_gamma_F_Ux))
        {
            PRNT("    @k %d: daux_diff_gamma_F_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_gamma_F_Ux);
            {
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_alpha_R_Ux= %g\n",daux_diff_alpha_R_Ux);
            }
            return(0);
        }
// aux_deriv[148]
    daux_diff_2gamma_F_Ux_Ux=abs(p[6][0])*pow(sec(aux_alpha_R),2.0)*mcond(tan(aux_alpha_R)>0.0,1.0,1,-1.0)*daux_diff_2alpha_R_Ux_Ux+abs(p[6][0])*mcond((x[4]-x[5]*p[5][0])/(0.001+abs(x[3]))>0.0,0.0,1,0.0)*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_Ux+2.0*abs(p[6][0])*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*mcond(tan(aux_alpha_R)>0.0,1.0,1,-1.0)*pow(daux_diff_alpha_R_Ux,2.0);
    if (isNANorINF(daux_diff_2gamma_F_Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2gamma_F_Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2gamma_F_Ux_Ux);
            {
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_2alpha_R_Ux_Ux= %g\n",daux_diff_2alpha_R_Ux_Ux);
                PRNT("        daux_diff_alpha_R_Ux= %g\n",daux_diff_alpha_R_Ux);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[149]
    daux_diff_2gamma_F_Ux_Uy=abs(p[6][0])*pow(sec(aux_alpha_R),2.0)*mcond(tan(aux_alpha_R)>0.0,1.0,1,-1.0)*daux_diff_2alpha_R_Ux_Uy+abs(p[6][0])*mcond((x[4]-x[5]*p[5][0])/(0.001+abs(x[3]))>0.0,0.0,1,0.0)*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_Ux+2.0*abs(p[6][0])*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*mcond(tan(aux_alpha_R)>0.0,1.0,1,-1.0)*daux_diff_alpha_R_Ux*daux_diff_alpha_R_Uy;
    if (isNANorINF(daux_diff_2gamma_F_Ux_Uy))
        {
            PRNT("    @k %d: daux_diff_2gamma_F_Ux_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2gamma_F_Ux_Uy);
            {
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_2alpha_R_Ux_Uy= %g\n",daux_diff_2alpha_R_Ux_Uy);
                PRNT("        daux_diff_alpha_R_Ux= %g\n",daux_diff_alpha_R_Ux);
                PRNT("        daux_diff_alpha_R_Uy= %g\n",daux_diff_alpha_R_Uy);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[150]
    daux_diff_2gamma_F_Ux_r=abs(p[6][0])*pow(sec(aux_alpha_R),2.0)*mcond(tan(aux_alpha_R)>0.0,1.0,1,-1.0)*daux_diff_2alpha_R_Ux_r+abs(p[6][0])*mcond((x[4]-x[5]*p[5][0])/(0.001+abs(x[3]))>0.0,0.0,1,0.0)*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_Ux+2.0*abs(p[6][0])*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*mcond(tan(aux_alpha_R)>0.0,1.0,1,-1.0)*daux_diff_alpha_R_Ux*daux_diff_alpha_R_r;
    if (isNANorINF(daux_diff_2gamma_F_Ux_r))
        {
            PRNT("    @k %d: daux_diff_2gamma_F_Ux_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2gamma_F_Ux_r);
            {
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_2alpha_R_Ux_r= %g\n",daux_diff_2alpha_R_Ux_r);
                PRNT("        daux_diff_alpha_R_Ux= %g\n",daux_diff_alpha_R_Ux);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[151]
    daux_diff_gamma_F_Uy=abs(p[6][0])*pow(sec(aux_alpha_R),2.0)*mcond(tan(aux_alpha_R)>0.0,1.0,1,-1.0)*daux_diff_alpha_R_Uy;
    if (isNANorINF(daux_diff_gamma_F_Uy))
        {
            PRNT("    @k %d: daux_diff_gamma_F_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_gamma_F_Uy);
            {
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_alpha_R_Uy= %g\n",daux_diff_alpha_R_Uy);
            }
            return(0);
        }
// aux_deriv[152]
    daux_diff_2gamma_F_Uy_Uy=abs(p[6][0])*pow(sec(aux_alpha_R),2.0)*mcond(tan(aux_alpha_R)>0.0,1.0,1,-1.0)*daux_diff_2alpha_R_Uy_Uy+abs(p[6][0])*mcond((x[4]-x[5]*p[5][0])/(0.001+abs(x[3]))>0.0,0.0,1,0.0)*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_Uy+2.0*abs(p[6][0])*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*mcond(tan(aux_alpha_R)>0.0,1.0,1,-1.0)*pow(daux_diff_alpha_R_Uy,2.0);
    if (isNANorINF(daux_diff_2gamma_F_Uy_Uy))
        {
            PRNT("    @k %d: daux_diff_2gamma_F_Uy_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2gamma_F_Uy_Uy);
            {
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_2alpha_R_Uy_Uy= %g\n",daux_diff_2alpha_R_Uy_Uy);
                PRNT("        daux_diff_alpha_R_Uy= %g\n",daux_diff_alpha_R_Uy);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[153]
    daux_diff_2gamma_F_Uy_r=abs(p[6][0])*pow(sec(aux_alpha_R),2.0)*mcond(tan(aux_alpha_R)>0.0,1.0,1,-1.0)*daux_diff_2alpha_R_Uy_r+abs(p[6][0])*mcond((x[4]-x[5]*p[5][0])/(0.001+abs(x[3]))>0.0,0.0,1,0.0)*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_Uy+2.0*abs(p[6][0])*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*mcond(tan(aux_alpha_R)>0.0,1.0,1,-1.0)*daux_diff_alpha_R_Uy*daux_diff_alpha_R_r;
    if (isNANorINF(daux_diff_2gamma_F_Uy_r))
        {
            PRNT("    @k %d: daux_diff_2gamma_F_Uy_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2gamma_F_Uy_r);
            {
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_2alpha_R_Uy_r= %g\n",daux_diff_2alpha_R_Uy_r);
                PRNT("        daux_diff_alpha_R_Uy= %g\n",daux_diff_alpha_R_Uy);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[154]
    daux_diff_gamma_F_r=abs(p[6][0])*pow(sec(aux_alpha_R),2.0)*mcond(tan(aux_alpha_R)>0.0,1.0,1,-1.0)*daux_diff_alpha_R_r;
    if (isNANorINF(daux_diff_gamma_F_r))
        {
            PRNT("    @k %d: daux_diff_gamma_F_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_gamma_F_r);
            {
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
            }
            return(0);
        }
// aux_deriv[155]
    daux_diff_2gamma_F_r_r=abs(p[6][0])*pow(sec(aux_alpha_R),2.0)*mcond(tan(aux_alpha_R)>0.0,1.0,1,-1.0)*daux_diff_2alpha_R_r_r+abs(p[6][0])*mcond((x[4]-x[5]*p[5][0])/(0.001+abs(x[3]))>0.0,0.0,1,0.0)*pow(sec(aux_alpha_R),2.0)*daux_diff_alpha_R_r+2.0*abs(p[6][0])*pow(sec(aux_alpha_R),2.0)*tan(aux_alpha_R)*mcond(tan(aux_alpha_R)>0.0,1.0,1,-1.0)*pow(daux_diff_alpha_R_r,2.0);
    if (isNANorINF(daux_diff_2gamma_F_r_r))
        {
            PRNT("    @k %d: daux_diff_2gamma_F_r_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2gamma_F_r_r);
            {
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        daux_diff_2alpha_R_r_r= %g\n",daux_diff_2alpha_R_r_r);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
            }
            return(0);
        }
// aux_deriv[156]
    daux_diff_F_f_Ux=mcond(abs(p[6][0])*abs(tan(aux_alpha_R))<=3.0*p[0][0]*p[21][0],daux_diff_gamma_F_Ux+0.3333333333333333*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*pow(tan(aux_alpha_R),2.0)*daux_diff_gamma_F_Ux/pow(p[0][0],2.0)/pow(p[21][0],2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*abs(p[6][0])*abs(tan(aux_alpha_R))*daux_diff_gamma_F_Ux/p[0][0]/p[21][0],1,0.0);
    if (isNANorINF(daux_diff_F_f_Ux))
        {
            PRNT("    @k %d: daux_diff_F_f_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_F_f_Ux);
            {
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[0,0]= %g\n",p[0][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_gamma_F_Ux= %g\n",daux_diff_gamma_F_Ux);
                PRNT("        p[22,0]= %g\n",p[22][0]);
            }
            return(0);
        }
// aux_deriv[157]
    daux_diff_2F_f_Ux_Ux=mcond(abs(-x[4]+x[5]*p[5][0])*abs(p[6][0])/(0.001+abs(x[3]))<=3.0*p[0][0]*p[21][0],daux_diff_2gamma_F_Ux_Ux+0.3333333333333333*pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*daux_diff_2gamma_F_Ux_Ux/pow(p[0][0],2.0)/pow(p[21][0],2.0)/pow(0.001+abs(x[3]),2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*abs(-x[4]+x[5]*p[5][0])*abs(p[6][0])*daux_diff_2gamma_F_Ux_Ux/p[0][0]/p[21][0]/(0.001+abs(x[3]))+0.6666666666666666*(x[4]-x[5]*p[5][0])*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*daux_diff_alpha_R_Ux*daux_diff_gamma_F_Ux/pow(p[0][0],2.0)/pow(p[21][0],2.0)/(0.001+abs(x[3]))-0.6666666666666666*(2.0-p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*abs(p[6][0])*mcond((x[4]-x[5]*p[5][0])/(0.001+abs(x[3]))>0.0,1.0,1,-1.0)*daux_diff_alpha_R_Ux*daux_diff_gamma_F_Ux/p[0][0]/p[21][0],1,0.0);
    if (isNANorINF(daux_diff_2F_f_Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2F_f_Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2F_f_Ux_Ux);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[0,0]= %g\n",p[0][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_2gamma_F_Ux_Ux= %g\n",daux_diff_2gamma_F_Ux_Ux);
                PRNT("        p[22,0]= %g\n",p[22][0]);
                PRNT("        daux_diff_alpha_R_Ux= %g\n",daux_diff_alpha_R_Ux);
                PRNT("        daux_diff_gamma_F_Ux= %g\n",daux_diff_gamma_F_Ux);
            }
            return(0);
        }
// aux_deriv[158]
    daux_diff_2F_f_Ux_Uy=mcond(abs(-x[4]+x[5]*p[5][0])*abs(p[6][0])/(0.001+abs(x[3]))<=3.0*p[0][0]*p[21][0],daux_diff_2gamma_F_Ux_Uy+0.3333333333333333*pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*daux_diff_2gamma_F_Ux_Uy/pow(p[0][0],2.0)/pow(p[21][0],2.0)/pow(0.001+abs(x[3]),2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*abs(-x[4]+x[5]*p[5][0])*abs(p[6][0])*daux_diff_2gamma_F_Ux_Uy/p[0][0]/p[21][0]/(0.001+abs(x[3]))+0.6666666666666666*(x[4]-x[5]*p[5][0])*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*daux_diff_alpha_R_Uy*daux_diff_gamma_F_Ux/pow(p[0][0],2.0)/pow(p[21][0],2.0)/(0.001+abs(x[3]))-0.6666666666666666*(2.0-p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*abs(p[6][0])*mcond((x[4]-x[5]*p[5][0])/(0.001+abs(x[3]))>0.0,1.0,1,-1.0)*daux_diff_alpha_R_Uy*daux_diff_gamma_F_Ux/p[0][0]/p[21][0],1,0.0);
    if (isNANorINF(daux_diff_2F_f_Ux_Uy))
        {
            PRNT("    @k %d: daux_diff_2F_f_Ux_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2F_f_Ux_Uy);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[0,0]= %g\n",p[0][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_2gamma_F_Ux_Uy= %g\n",daux_diff_2gamma_F_Ux_Uy);
                PRNT("        p[22,0]= %g\n",p[22][0]);
                PRNT("        daux_diff_alpha_R_Uy= %g\n",daux_diff_alpha_R_Uy);
                PRNT("        daux_diff_gamma_F_Ux= %g\n",daux_diff_gamma_F_Ux);
            }
            return(0);
        }
// aux_deriv[159]
    daux_diff_2F_f_Ux_r=mcond(abs(-x[4]+x[5]*p[5][0])*abs(p[6][0])/(0.001+abs(x[3]))<=3.0*p[0][0]*p[21][0],daux_diff_2gamma_F_Ux_r+0.3333333333333333*pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*daux_diff_2gamma_F_Ux_r/pow(p[0][0],2.0)/pow(p[21][0],2.0)/pow(0.001+abs(x[3]),2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*abs(-x[4]+x[5]*p[5][0])*abs(p[6][0])*daux_diff_2gamma_F_Ux_r/p[0][0]/p[21][0]/(0.001+abs(x[3]))+0.6666666666666666*(x[4]-x[5]*p[5][0])*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*daux_diff_alpha_R_r*daux_diff_gamma_F_Ux/pow(p[0][0],2.0)/pow(p[21][0],2.0)/(0.001+abs(x[3]))-0.6666666666666666*(2.0-p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*abs(p[6][0])*mcond((x[4]-x[5]*p[5][0])/(0.001+abs(x[3]))>0.0,1.0,1,-1.0)*daux_diff_alpha_R_r*daux_diff_gamma_F_Ux/p[0][0]/p[21][0],1,0.0);
    if (isNANorINF(daux_diff_2F_f_Ux_r))
        {
            PRNT("    @k %d: daux_diff_2F_f_Ux_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2F_f_Ux_r);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[0,0]= %g\n",p[0][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_2gamma_F_Ux_r= %g\n",daux_diff_2gamma_F_Ux_r);
                PRNT("        p[22,0]= %g\n",p[22][0]);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
                PRNT("        daux_diff_gamma_F_Ux= %g\n",daux_diff_gamma_F_Ux);
            }
            return(0);
        }
// aux_deriv[160]
    daux_diff_F_f_Uy=mcond(abs(p[6][0])*abs(tan(aux_alpha_R))<=3.0*p[0][0]*p[21][0],daux_diff_gamma_F_Uy+0.3333333333333333*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*pow(tan(aux_alpha_R),2.0)*daux_diff_gamma_F_Uy/pow(p[0][0],2.0)/pow(p[21][0],2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*abs(p[6][0])*abs(tan(aux_alpha_R))*daux_diff_gamma_F_Uy/p[0][0]/p[21][0],1,0.0);
    if (isNANorINF(daux_diff_F_f_Uy))
        {
            PRNT("    @k %d: daux_diff_F_f_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_F_f_Uy);
            {
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[0,0]= %g\n",p[0][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_gamma_F_Uy= %g\n",daux_diff_gamma_F_Uy);
                PRNT("        p[22,0]= %g\n",p[22][0]);
            }
            return(0);
        }
// aux_deriv[161]
    daux_diff_2F_f_Uy_Uy=mcond(abs(-x[4]+x[5]*p[5][0])*abs(p[6][0])/(0.001+abs(x[3]))<=3.0*p[0][0]*p[21][0],daux_diff_2gamma_F_Uy_Uy+0.3333333333333333*pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*daux_diff_2gamma_F_Uy_Uy/pow(p[0][0],2.0)/pow(p[21][0],2.0)/pow(0.001+abs(x[3]),2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*abs(-x[4]+x[5]*p[5][0])*abs(p[6][0])*daux_diff_2gamma_F_Uy_Uy/p[0][0]/p[21][0]/(0.001+abs(x[3]))+0.6666666666666666*(x[4]-x[5]*p[5][0])*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*daux_diff_alpha_R_Uy*daux_diff_gamma_F_Uy/pow(p[0][0],2.0)/pow(p[21][0],2.0)/(0.001+abs(x[3]))-0.6666666666666666*(2.0-p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*abs(p[6][0])*mcond((x[4]-x[5]*p[5][0])/(0.001+abs(x[3]))>0.0,1.0,1,-1.0)*daux_diff_alpha_R_Uy*daux_diff_gamma_F_Uy/p[0][0]/p[21][0],1,0.0);
    if (isNANorINF(daux_diff_2F_f_Uy_Uy))
        {
            PRNT("    @k %d: daux_diff_2F_f_Uy_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2F_f_Uy_Uy);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[0,0]= %g\n",p[0][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_2gamma_F_Uy_Uy= %g\n",daux_diff_2gamma_F_Uy_Uy);
                PRNT("        p[22,0]= %g\n",p[22][0]);
                PRNT("        daux_diff_alpha_R_Uy= %g\n",daux_diff_alpha_R_Uy);
                PRNT("        daux_diff_gamma_F_Uy= %g\n",daux_diff_gamma_F_Uy);
            }
            return(0);
        }
// aux_deriv[162]
    daux_diff_2F_f_Uy_r=mcond(abs(-x[4]+x[5]*p[5][0])*abs(p[6][0])/(0.001+abs(x[3]))<=3.0*p[0][0]*p[21][0],daux_diff_2gamma_F_Uy_r+0.3333333333333333*pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*daux_diff_2gamma_F_Uy_r/pow(p[0][0],2.0)/pow(p[21][0],2.0)/pow(0.001+abs(x[3]),2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*abs(-x[4]+x[5]*p[5][0])*abs(p[6][0])*daux_diff_2gamma_F_Uy_r/p[0][0]/p[21][0]/(0.001+abs(x[3]))+0.6666666666666666*(x[4]-x[5]*p[5][0])*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*daux_diff_alpha_R_r*daux_diff_gamma_F_Uy/pow(p[0][0],2.0)/pow(p[21][0],2.0)/(0.001+abs(x[3]))-0.6666666666666666*(2.0-p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*abs(p[6][0])*mcond((x[4]-x[5]*p[5][0])/(0.001+abs(x[3]))>0.0,1.0,1,-1.0)*daux_diff_alpha_R_r*daux_diff_gamma_F_Uy/p[0][0]/p[21][0],1,0.0);
    if (isNANorINF(daux_diff_2F_f_Uy_r))
        {
            PRNT("    @k %d: daux_diff_2F_f_Uy_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2F_f_Uy_r);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[0,0]= %g\n",p[0][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_2gamma_F_Uy_r= %g\n",daux_diff_2gamma_F_Uy_r);
                PRNT("        p[22,0]= %g\n",p[22][0]);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
                PRNT("        daux_diff_gamma_F_Uy= %g\n",daux_diff_gamma_F_Uy);
            }
            return(0);
        }
// aux_deriv[163]
    daux_diff_F_f_r=mcond(abs(p[6][0])*abs(tan(aux_alpha_R))<=3.0*p[0][0]*p[21][0],daux_diff_gamma_F_r+0.3333333333333333*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*pow(tan(aux_alpha_R),2.0)*daux_diff_gamma_F_r/pow(p[0][0],2.0)/pow(p[21][0],2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*abs(p[6][0])*abs(tan(aux_alpha_R))*daux_diff_gamma_F_r/p[0][0]/p[21][0],1,0.0);
    if (isNANorINF(daux_diff_F_f_r))
        {
            PRNT("    @k %d: daux_diff_F_f_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_F_f_r);
            {
                PRNT("        aux_alpha_R= %g\n",aux_alpha_R);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        p[0,0]= %g\n",p[0][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_gamma_F_r= %g\n",daux_diff_gamma_F_r);
                PRNT("        p[22,0]= %g\n",p[22][0]);
            }
            return(0);
        }
// aux_deriv[164]
    daux_diff_2F_f_r_r=mcond(abs(-x[4]+x[5]*p[5][0])*abs(p[6][0])/(0.001+abs(x[3]))<=3.0*p[0][0]*p[21][0],daux_diff_2gamma_F_r_r+0.3333333333333333*pow(x[4]-x[5]*p[5][0],2.0)*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*daux_diff_2gamma_F_r_r/pow(p[0][0],2.0)/pow(p[21][0],2.0)/pow(0.001+abs(x[3]),2.0)-0.6666666666666666*(2.0-p[22][0]/p[21][0])*abs(-x[4]+x[5]*p[5][0])*abs(p[6][0])*daux_diff_2gamma_F_r_r/p[0][0]/p[21][0]/(0.001+abs(x[3]))+0.6666666666666666*(x[4]-x[5]*p[5][0])*pow(p[6][0],2.0)*(1.0-0.6666666666666666*p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*daux_diff_alpha_R_r*daux_diff_gamma_F_r/pow(p[0][0],2.0)/pow(p[21][0],2.0)/(0.001+abs(x[3]))-0.6666666666666666*(2.0-p[22][0]/p[21][0])*(1.0+pow(x[4]-x[5]*p[5][0],2.0)/pow(0.001+abs(x[3]),2.0))*abs(p[6][0])*mcond((x[4]-x[5]*p[5][0])/(0.001+abs(x[3]))>0.0,1.0,1,-1.0)*daux_diff_alpha_R_r*daux_diff_gamma_F_r/p[0][0]/p[21][0],1,0.0);
    if (isNANorINF(daux_diff_2F_f_r_r))
        {
            PRNT("    @k %d: daux_diff_2F_f_r_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2F_f_r_r);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[6,0]= %g\n",p[6][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[5,0]= %g\n",p[5][0]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[0,0]= %g\n",p[0][0]);
                PRNT("        p[21,0]= %g\n",p[21][0]);
                PRNT("        daux_diff_2gamma_F_r_r= %g\n",daux_diff_2gamma_F_r_r);
                PRNT("        p[22,0]= %g\n",p[22][0]);
                PRNT("        daux_diff_alpha_R_r= %g\n",daux_diff_alpha_R_r);
                PRNT("        daux_diff_gamma_F_r= %g\n",daux_diff_gamma_F_r);
            }
            return(0);
        }
// aux_deriv[165]
    daux_diff_dUy_Ux=-x[5]+(cos(u[1])*daux_diff_Fyf_Ux+daux_diff_Fyr_Ux)/p[20][0];
    if (isNANorINF(daux_diff_dUy_Ux))
        {
            PRNT("    @k %d: daux_diff_dUy_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dUy_Ux);
            {
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        daux_diff_Fyr_Ux= %g\n",daux_diff_Fyr_Ux);
                PRNT("        daux_diff_Fyf_Ux= %g\n",daux_diff_Fyf_Ux);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[166]
    daux_diff_2dUy_Ux_Ux=(cos(u[1])*daux_diff_2Fyf_Ux_Ux+daux_diff_2Fyr_Ux_Ux)/p[20][0];
    if (isNANorINF(daux_diff_2dUy_Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2dUy_Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUy_Ux_Ux);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        daux_diff_2Fyr_Ux_Ux= %g\n",daux_diff_2Fyr_Ux_Ux);
                PRNT("        daux_diff_2Fyf_Ux_Ux= %g\n",daux_diff_2Fyf_Ux_Ux);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[167]
    daux_diff_2dUy_Ux_Uy=(cos(u[1])*daux_diff_2Fyf_Ux_Uy+daux_diff_2Fyr_Ux_Uy)/p[20][0];
    if (isNANorINF(daux_diff_2dUy_Ux_Uy))
        {
            PRNT("    @k %d: daux_diff_2dUy_Ux_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUy_Ux_Uy);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        daux_diff_2Fyr_Ux_Uy= %g\n",daux_diff_2Fyr_Ux_Uy);
                PRNT("        daux_diff_2Fyf_Ux_Uy= %g\n",daux_diff_2Fyf_Ux_Uy);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[168]
    daux_diff_2dUy_Ux_r=-1.0+(cos(u[1])*daux_diff_2Fyf_Ux_r+daux_diff_2Fyr_Ux_r)/p[20][0];
    if (isNANorINF(daux_diff_2dUy_Ux_r))
        {
            PRNT("    @k %d: daux_diff_2dUy_Ux_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUy_Ux_r);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        daux_diff_2Fyr_Ux_r= %g\n",daux_diff_2Fyr_Ux_r);
                PRNT("        daux_diff_2Fyf_Ux_r= %g\n",daux_diff_2Fyf_Ux_r);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[169]
    daux_diff_2dUy_Ux_steer=(cos(u[1])*daux_diff_2Fyf_Ux_steer-sin(u[1])*daux_diff_Fyf_Ux)/p[20][0];
    if (isNANorINF(daux_diff_2dUy_Ux_steer))
        {
            PRNT("    @k %d: daux_diff_2dUy_Ux_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUy_Ux_steer);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        daux_diff_2Fyf_Ux_steer= %g\n",daux_diff_2Fyf_Ux_steer);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_Fyf_Ux= %g\n",daux_diff_Fyf_Ux);
            }
            return(0);
        }
// aux_deriv[170]
    daux_diff_2dUy_Ux_thr=daux_diff_2Fyr_Ux_thr/p[20][0];
    if (isNANorINF(daux_diff_2dUy_Ux_thr))
        {
            PRNT("    @k %d: daux_diff_2dUy_Ux_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUy_Ux_thr);
            {
                PRNT("        daux_diff_2Fyr_Ux_thr= %g\n",daux_diff_2Fyr_Ux_thr);
                PRNT("        p[20,0]= %g\n",p[20][0]);
            }
            return(0);
        }
// aux_deriv[171]
    daux_diff_dUy_Uy=(cos(u[1])*daux_diff_Fyf_Uy+daux_diff_Fyr_Uy)/p[20][0];
    if (isNANorINF(daux_diff_dUy_Uy))
        {
            PRNT("    @k %d: daux_diff_dUy_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dUy_Uy);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        daux_diff_Fyr_Uy= %g\n",daux_diff_Fyr_Uy);
                PRNT("        daux_diff_Fyf_Uy= %g\n",daux_diff_Fyf_Uy);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[172]
    daux_diff_2dUy_Uy_Uy=(cos(u[1])*daux_diff_2Fyf_Uy_Uy+daux_diff_2Fyr_Uy_Uy)/p[20][0];
    if (isNANorINF(daux_diff_2dUy_Uy_Uy))
        {
            PRNT("    @k %d: daux_diff_2dUy_Uy_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUy_Uy_Uy);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        daux_diff_2Fyr_Uy_Uy= %g\n",daux_diff_2Fyr_Uy_Uy);
                PRNT("        daux_diff_2Fyf_Uy_Uy= %g\n",daux_diff_2Fyf_Uy_Uy);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[173]
    daux_diff_2dUy_Uy_r=(cos(u[1])*daux_diff_2Fyf_Uy_r+daux_diff_2Fyr_Uy_r)/p[20][0];
    if (isNANorINF(daux_diff_2dUy_Uy_r))
        {
            PRNT("    @k %d: daux_diff_2dUy_Uy_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUy_Uy_r);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        daux_diff_2Fyr_Uy_r= %g\n",daux_diff_2Fyr_Uy_r);
                PRNT("        daux_diff_2Fyf_Uy_r= %g\n",daux_diff_2Fyf_Uy_r);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[174]
    daux_diff_2dUy_Uy_steer=(cos(u[1])*daux_diff_2Fyf_Uy_steer-sin(u[1])*daux_diff_Fyf_Uy)/p[20][0];
    if (isNANorINF(daux_diff_2dUy_Uy_steer))
        {
            PRNT("    @k %d: daux_diff_2dUy_Uy_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUy_Uy_steer);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        daux_diff_2Fyf_Uy_steer= %g\n",daux_diff_2Fyf_Uy_steer);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_Fyf_Uy= %g\n",daux_diff_Fyf_Uy);
            }
            return(0);
        }
// aux_deriv[175]
    daux_diff_2dUy_Uy_thr=daux_diff_2Fyr_Uy_thr/p[20][0];
    if (isNANorINF(daux_diff_2dUy_Uy_thr))
        {
            PRNT("    @k %d: daux_diff_2dUy_Uy_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUy_Uy_thr);
            {
                PRNT("        daux_diff_2Fyr_Uy_thr= %g\n",daux_diff_2Fyr_Uy_thr);
                PRNT("        p[20,0]= %g\n",p[20][0]);
            }
            return(0);
        }
// aux_deriv[176]
    daux_diff_dUy_r=-x[3]+(cos(u[1])*daux_diff_Fyf_r+daux_diff_Fyr_r)/p[20][0];
    if (isNANorINF(daux_diff_dUy_r))
        {
            PRNT("    @k %d: daux_diff_dUy_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dUy_r);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        daux_diff_Fyr_r= %g\n",daux_diff_Fyr_r);
                PRNT("        daux_diff_Fyf_r= %g\n",daux_diff_Fyf_r);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[177]
    daux_diff_2dUy_r_r=(cos(u[1])*daux_diff_2Fyf_r_r+daux_diff_2Fyr_r_r)/p[20][0];
    if (isNANorINF(daux_diff_2dUy_r_r))
        {
            PRNT("    @k %d: daux_diff_2dUy_r_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUy_r_r);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        daux_diff_2Fyr_r_r= %g\n",daux_diff_2Fyr_r_r);
                PRNT("        daux_diff_2Fyf_r_r= %g\n",daux_diff_2Fyf_r_r);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[178]
    daux_diff_2dUy_r_steer=(cos(u[1])*daux_diff_2Fyf_r_steer-sin(u[1])*daux_diff_Fyf_r)/p[20][0];
    if (isNANorINF(daux_diff_2dUy_r_steer))
        {
            PRNT("    @k %d: daux_diff_2dUy_r_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUy_r_steer);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        daux_diff_2Fyf_r_steer= %g\n",daux_diff_2Fyf_r_steer);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_Fyf_r= %g\n",daux_diff_Fyf_r);
            }
            return(0);
        }
// aux_deriv[179]
    daux_diff_2dUy_r_thr=daux_diff_2Fyr_r_thr/p[20][0];
    if (isNANorINF(daux_diff_2dUy_r_thr))
        {
            PRNT("    @k %d: daux_diff_2dUy_r_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUy_r_thr);
            {
                PRNT("        daux_diff_2Fyr_r_thr= %g\n",daux_diff_2Fyr_r_thr);
                PRNT("        p[20,0]= %g\n",p[20][0]);
            }
            return(0);
        }
// aux_deriv[180]
    daux_diff_dUy_steer=(-sin(u[1])*aux_Fyf+cos(u[1])*daux_diff_Fyf_steer)/p[20][0];
    if (isNANorINF(daux_diff_dUy_steer))
        {
            PRNT("    @k %d: daux_diff_dUy_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dUy_steer);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        daux_diff_Fyf_steer= %g\n",daux_diff_Fyf_steer);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
            }
            return(0);
        }
// aux_deriv[181]
    daux_diff_2dUy_steer_steer=(-cos(u[1])*aux_Fyf+cos(u[1])*daux_diff_2Fyf_steer_steer-2.0*sin(u[1])*daux_diff_Fyf_steer)/p[20][0];
    if (isNANorINF(daux_diff_2dUy_steer_steer))
        {
            PRNT("    @k %d: daux_diff_2dUy_steer_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUy_steer_steer);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_2Fyf_steer_steer= %g\n",daux_diff_2Fyf_steer_steer);
                PRNT("        daux_diff_Fyf_steer= %g\n",daux_diff_Fyf_steer);
            }
            return(0);
        }
// aux_deriv[182]
    daux_diff_dUy_thr=daux_diff_Fyr_thr/p[20][0];
    if (isNANorINF(daux_diff_dUy_thr))
        {
            PRNT("    @k %d: daux_diff_dUy_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dUy_thr);
            {
                PRNT("        daux_diff_Fyr_thr= %g\n",daux_diff_Fyr_thr);
                PRNT("        p[20,0]= %g\n",p[20][0]);
            }
            return(0);
        }
// aux_deriv[183]
    daux_diff_2dUy_thr_thr=daux_diff_2Fyr_thr_thr/p[20][0];
    if (isNANorINF(daux_diff_2dUy_thr_thr))
        {
            PRNT("    @k %d: daux_diff_2dUy_thr_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUy_thr_thr);
            {
                PRNT("        daux_diff_2Fyr_thr_thr= %g\n",daux_diff_2Fyr_thr_thr);
                PRNT("        p[20,0]= %g\n",p[20][0]);
            }
            return(0);
        }
// aux_deriv[184]
    daux_diff_dUx_Ux=(derivative(aux_Fxr,x[3],1.0)-sin(u[1])*daux_diff_Fyf_Ux)/p[20][0];
    if (isNANorINF(daux_diff_dUx_Ux))
        {
            PRNT("    @k %d: daux_diff_dUx_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dUx_Ux);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        aux_Fxr= %g\n",aux_Fxr);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        daux_diff_Fyf_Ux= %g\n",daux_diff_Fyf_Ux);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[185]
    daux_diff_2dUx_Ux_Ux=(derivative(aux_Fxr,x[3],2.0)-sin(u[1])*daux_diff_2Fyf_Ux_Ux)/p[20][0];
    if (isNANorINF(daux_diff_2dUx_Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2dUx_Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUx_Ux_Ux);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        aux_Fxr= %g\n",aux_Fxr);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        daux_diff_2Fyf_Ux_Ux= %g\n",daux_diff_2Fyf_Ux_Ux);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[186]
    daux_diff_2dUx_Ux_Uy=(derivative(aux_Fxr,x[3],1.0,x[4],1.0)-sin(u[1])*daux_diff_2Fyf_Ux_Uy)/p[20][0];
    if (isNANorINF(daux_diff_2dUx_Ux_Uy))
        {
            PRNT("    @k %d: daux_diff_2dUx_Ux_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUx_Ux_Uy);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        aux_Fxr= %g\n",aux_Fxr);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        daux_diff_2Fyf_Ux_Uy= %g\n",daux_diff_2Fyf_Ux_Uy);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[187]
    daux_diff_2dUx_Ux_r=(derivative(aux_Fxr,x[3],1.0,x[5],1.0)-sin(u[1])*daux_diff_2Fyf_Ux_r)/p[20][0];
    if (isNANorINF(daux_diff_2dUx_Ux_r))
        {
            PRNT("    @k %d: daux_diff_2dUx_Ux_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUx_Ux_r);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        aux_Fxr= %g\n",aux_Fxr);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        daux_diff_2Fyf_Ux_r= %g\n",daux_diff_2Fyf_Ux_r);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[188]
    daux_diff_2dUx_Ux_steer=(-sin(u[1])*daux_diff_2Fyf_Ux_steer-cos(u[1])*daux_diff_Fyf_Ux)/p[20][0];
    if (isNANorINF(daux_diff_2dUx_Ux_steer))
        {
            PRNT("    @k %d: daux_diff_2dUx_Ux_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUx_Ux_steer);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        daux_diff_Fyf_Ux= %g\n",daux_diff_Fyf_Ux);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_2Fyf_Ux_steer= %g\n",daux_diff_2Fyf_Ux_steer);
            }
            return(0);
        }
// aux_deriv[189]
    daux_diff_2dUx_Ux_thr=derivative(aux_Fxr,u[0],1.0,x[3],1.0)/p[20][0];
    if (isNANorINF(daux_diff_2dUx_Ux_thr))
        {
            PRNT("    @k %d: daux_diff_2dUx_Ux_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUx_Ux_thr);
            {
                PRNT("        aux_Fxr= %g\n",aux_Fxr);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        u[0]= %g\n",u[0]);
                PRNT("        p[20,0]= %g\n",p[20][0]);
            }
            return(0);
        }
// aux_deriv[190]
    daux_diff_dUx_Uy=x[5]+(derivative(aux_Fxr,x[4],1.0)-sin(u[1])*daux_diff_Fyf_Uy)/p[20][0];
    if (isNANorINF(daux_diff_dUx_Uy))
        {
            PRNT("    @k %d: daux_diff_dUx_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dUx_Uy);
            {
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        aux_Fxr= %g\n",aux_Fxr);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        daux_diff_Fyf_Uy= %g\n",daux_diff_Fyf_Uy);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[191]
    daux_diff_2dUx_Uy_Uy=(derivative(aux_Fxr,x[4],2.0)-sin(u[1])*daux_diff_2Fyf_Uy_Uy)/p[20][0];
    if (isNANorINF(daux_diff_2dUx_Uy_Uy))
        {
            PRNT("    @k %d: daux_diff_2dUx_Uy_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUx_Uy_Uy);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        aux_Fxr= %g\n",aux_Fxr);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        daux_diff_2Fyf_Uy_Uy= %g\n",daux_diff_2Fyf_Uy_Uy);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[192]
    daux_diff_2dUx_Uy_r=1.0+(derivative(aux_Fxr,x[4],1.0,x[5],1.0)-sin(u[1])*daux_diff_2Fyf_Uy_r)/p[20][0];
    if (isNANorINF(daux_diff_2dUx_Uy_r))
        {
            PRNT("    @k %d: daux_diff_2dUx_Uy_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUx_Uy_r);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        aux_Fxr= %g\n",aux_Fxr);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        daux_diff_2Fyf_Uy_r= %g\n",daux_diff_2Fyf_Uy_r);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[193]
    daux_diff_2dUx_Uy_steer=(-sin(u[1])*daux_diff_2Fyf_Uy_steer-cos(u[1])*daux_diff_Fyf_Uy)/p[20][0];
    if (isNANorINF(daux_diff_2dUx_Uy_steer))
        {
            PRNT("    @k %d: daux_diff_2dUx_Uy_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUx_Uy_steer);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        daux_diff_Fyf_Uy= %g\n",daux_diff_Fyf_Uy);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_2Fyf_Uy_steer= %g\n",daux_diff_2Fyf_Uy_steer);
            }
            return(0);
        }
// aux_deriv[194]
    daux_diff_2dUx_Uy_thr=derivative(aux_Fxr,u[0],1.0,x[4],1.0)/p[20][0];
    if (isNANorINF(daux_diff_2dUx_Uy_thr))
        {
            PRNT("    @k %d: daux_diff_2dUx_Uy_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUx_Uy_thr);
            {
                PRNT("        aux_Fxr= %g\n",aux_Fxr);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        u[0]= %g\n",u[0]);
                PRNT("        p[20,0]= %g\n",p[20][0]);
            }
            return(0);
        }
// aux_deriv[195]
    daux_diff_dUx_r=x[4]+(derivative(aux_Fxr,x[5],1.0)-sin(u[1])*daux_diff_Fyf_r)/p[20][0];
    if (isNANorINF(daux_diff_dUx_r))
        {
            PRNT("    @k %d: daux_diff_dUx_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dUx_r);
            {
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        aux_Fxr= %g\n",aux_Fxr);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        daux_diff_Fyf_r= %g\n",daux_diff_Fyf_r);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[196]
    daux_diff_2dUx_r_r=(derivative(aux_Fxr,x[5],2.0)-sin(u[1])*daux_diff_2Fyf_r_r)/p[20][0];
    if (isNANorINF(daux_diff_2dUx_r_r))
        {
            PRNT("    @k %d: daux_diff_2dUx_r_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUx_r_r);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        aux_Fxr= %g\n",aux_Fxr);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        daux_diff_2Fyf_r_r= %g\n",daux_diff_2Fyf_r_r);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }
// aux_deriv[197]
    daux_diff_2dUx_r_steer=(-sin(u[1])*daux_diff_2Fyf_r_steer-cos(u[1])*daux_diff_Fyf_r)/p[20][0];
    if (isNANorINF(daux_diff_2dUx_r_steer))
        {
            PRNT("    @k %d: daux_diff_2dUx_r_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUx_r_steer);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        daux_diff_Fyf_r= %g\n",daux_diff_Fyf_r);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_2Fyf_r_steer= %g\n",daux_diff_2Fyf_r_steer);
            }
            return(0);
        }
// aux_deriv[198]
    daux_diff_2dUx_r_thr=derivative(aux_Fxr,u[0],1.0,x[5],1.0)/p[20][0];
    if (isNANorINF(daux_diff_2dUx_r_thr))
        {
            PRNT("    @k %d: daux_diff_2dUx_r_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUx_r_thr);
            {
                PRNT("        aux_Fxr= %g\n",aux_Fxr);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        u[0]= %g\n",u[0]);
                PRNT("        p[20,0]= %g\n",p[20][0]);
            }
            return(0);
        }
// aux_deriv[199]
    daux_diff_dUx_steer=(-cos(u[1])*aux_Fyf-sin(u[1])*daux_diff_Fyf_steer)/p[20][0];
    if (isNANorINF(daux_diff_dUx_steer))
        {
            PRNT("    @k %d: daux_diff_dUx_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dUx_steer);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        daux_diff_Fyf_steer= %g\n",daux_diff_Fyf_steer);
            }
            return(0);
        }
// aux_deriv[200]
    daux_diff_2dUx_steer_steer=(sin(u[1])*aux_Fyf-sin(u[1])*daux_diff_2Fyf_steer_steer-2.0*cos(u[1])*daux_diff_Fyf_steer)/p[20][0];
    if (isNANorINF(daux_diff_2dUx_steer_steer))
        {
            PRNT("    @k %d: daux_diff_2dUx_steer_steer in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUx_steer_steer);
            {
                PRNT("        p[20,0]= %g\n",p[20][0]);
                PRNT("        daux_diff_Fyf_steer= %g\n",daux_diff_Fyf_steer);
                PRNT("        u[1]= %g\n",u[1]);
                PRNT("        aux_Fyf= %g\n",aux_Fyf);
                PRNT("        daux_diff_2Fyf_steer_steer= %g\n",daux_diff_2Fyf_steer_steer);
            }
            return(0);
        }
// aux_deriv[201]
    daux_diff_dUx_thr=derivative(aux_Fxr,u[0],1.0)/p[20][0];
    if (isNANorINF(daux_diff_dUx_thr))
        {
            PRNT("    @k %d: daux_diff_dUx_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_dUx_thr);
            {
                PRNT("        aux_Fxr= %g\n",aux_Fxr);
                PRNT("        u[0]= %g\n",u[0]);
                PRNT("        p[20,0]= %g\n",p[20][0]);
            }
            return(0);
        }
// aux_deriv[202]
    daux_diff_2dUx_thr_thr=derivative(aux_Fxr,u[0],2.0)/p[20][0];
    if (isNANorINF(daux_diff_2dUx_thr_thr))
        {
            PRNT("    @k %d: daux_diff_2dUx_thr_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2dUx_thr_thr);
            {
                PRNT("        aux_Fxr= %g\n",aux_Fxr);
                PRNT("        u[0]= %g\n",u[0]);
                PRNT("        p[20,0]= %g\n",p[20][0]);
            }
            return(0);
        }
// aux_deriv[203]
    daux_diff_rev_Ux=mcond((u[0]-x[3])/(0.001+abs(x[3]))<0.0,0.0,1,0.0);
    if (isNANorINF(daux_diff_rev_Ux))
        {
            PRNT("    @k %d: daux_diff_rev_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_rev_Ux);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
// aux_deriv[204]
    daux_diff_2rev_Ux_Ux=mcond((u[0]-x[3])/(0.001+abs(x[3]))<0.0,0.0,1,0.0);
    if (isNANorINF(daux_diff_2rev_Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2rev_Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2rev_Ux_Ux);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
// aux_deriv[205]
    daux_diff_2rev_Ux_thr=mcond((u[0]-x[3])/(0.001+abs(x[3]))<0.0,0.0,1,0.0);
    if (isNANorINF(daux_diff_2rev_Ux_thr))
        {
            PRNT("    @k %d: daux_diff_2rev_Ux_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2rev_Ux_thr);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
// aux_deriv[206]
    daux_diff_rev_thr=mcond((u[0]-x[3])/(0.001+abs(x[3]))<0.0,0.0,1,0.0);
    if (isNANorINF(daux_diff_rev_thr))
        {
            PRNT("    @k %d: daux_diff_rev_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_rev_thr);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
// aux_deriv[207]
    daux_diff_2rev_thr_thr=mcond((u[0]-x[3])/(0.001+abs(x[3]))<0.0,0.0,1,0.0);
    if (isNANorINF(daux_diff_2rev_thr_thr))
        {
            PRNT("    @k %d: daux_diff_2rev_thr_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2rev_thr_thr);
            {
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
// aux_deriv[208]
    daux_diff_Fxr_Ux=p[7][0]*aux_K*aux_rev*daux_diff_F_r_Ux/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_K_Ux/pow(1.0+aux_K,2.0)/aux_gamma_R+p[7][0]*aux_F_r*aux_rev*daux_diff_K_Ux/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_gamma_R_Ux/(1.0+aux_K)/pow(aux_gamma_R,2.0)+p[7][0]*aux_F_r*aux_K*daux_diff_rev_Ux/(1.0+aux_K)/aux_gamma_R;
    if (isNANorINF(daux_diff_Fxr_Ux))
        {
            PRNT("    @k %d: daux_diff_Fxr_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_Fxr_Ux);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        daux_diff_rev_Ux= %g\n",daux_diff_rev_Ux);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_gamma_R_Ux= %g\n",daux_diff_gamma_R_Ux);
                PRNT("        aux_rev= %g\n",aux_rev);
                PRNT("        daux_diff_F_r_Ux= %g\n",daux_diff_F_r_Ux);
                PRNT("        daux_diff_K_Ux= %g\n",daux_diff_K_Ux);
            }
            return(0);
        }
// aux_deriv[209]
    daux_diff_2Fxr_Ux_Ux=p[7][0]*aux_K*aux_rev*daux_diff_2F_r_Ux_Ux/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_2K_Ux_Ux/pow(1.0+aux_K,2.0)/aux_gamma_R+p[7][0]*aux_F_r*aux_rev*daux_diff_2K_Ux_Ux/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_2gamma_R_Ux_Ux/(1.0+aux_K)/pow(aux_gamma_R,2.0)+p[7][0]*aux_F_r*aux_K*daux_diff_2rev_Ux_Ux/(1.0+aux_K)/aux_gamma_R-2.0*p[7][0]*aux_K*aux_rev*daux_diff_F_r_Ux*daux_diff_K_Ux/pow(1.0+aux_K,2.0)/aux_gamma_R+2.0*p[7][0]*aux_rev*daux_diff_F_r_Ux*daux_diff_K_Ux/(1.0+aux_K)/aux_gamma_R+2.0*p[7][0]*aux_F_r*aux_K*aux_rev*pow(daux_diff_K_Ux,2.0)/pow(1.0+aux_K,3.0)/aux_gamma_R-2.0*p[7][0]*aux_F_r*aux_rev*pow(daux_diff_K_Ux,2.0)/pow(1.0+aux_K,2.0)/aux_gamma_R-2.0*p[7][0]*aux_K*aux_rev*daux_diff_F_r_Ux*daux_diff_gamma_R_Ux/(1.0+aux_K)/pow(aux_gamma_R,2.0)+2.0*p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_K_Ux*daux_diff_gamma_R_Ux/pow(1.0+aux_K,2.0)/pow(aux_gamma_R,2.0)-2.0*p[7][0]*aux_F_r*aux_rev*daux_diff_K_Ux*
     daux_diff_gamma_R_Ux/(1.0+aux_K)/pow(aux_gamma_R,2.0)+2.0*p[7][0]*aux_F_r*aux_K*aux_rev*pow(daux_diff_gamma_R_Ux,2.0)/(1.0+aux_K)/pow(aux_gamma_R,3.0)+2.0*p[7][0]*aux_K*daux_diff_F_r_Ux*daux_diff_rev_Ux/(1.0+aux_K)/aux_gamma_R-2.0*p[7][0]*aux_F_r*aux_K*daux_diff_K_Ux*daux_diff_rev_Ux/pow(1.0+aux_K,2.0)/aux_gamma_R+2.0*p[7][0]*aux_F_r*daux_diff_K_Ux*daux_diff_rev_Ux/(1.0+aux_K)/aux_gamma_R-2.0*p[7][0]*aux_F_r*aux_K*daux_diff_gamma_R_Ux*daux_diff_rev_Ux/(1.0+aux_K)/pow(aux_gamma_R,2.0);
    if (isNANorINF(daux_diff_2Fxr_Ux_Ux))
        {
            PRNT("    @k %d: daux_diff_2Fxr_Ux_Ux in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fxr_Ux_Ux);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        daux_diff_gamma_R_Ux= %g\n",daux_diff_gamma_R_Ux);
                PRNT("        daux_diff_rev_Ux= %g\n",daux_diff_rev_Ux);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_2rev_Ux_Ux= %g\n",daux_diff_2rev_Ux_Ux);
                PRNT("        daux_diff_F_r_Ux= %g\n",daux_diff_F_r_Ux);
                PRNT("        daux_diff_K_Ux= %g\n",daux_diff_K_Ux);
                PRNT("        aux_rev= %g\n",aux_rev);
                PRNT("        daux_diff_2gamma_R_Ux_Ux= %g\n",daux_diff_2gamma_R_Ux_Ux);
                PRNT("        daux_diff_2F_r_Ux_Ux= %g\n",daux_diff_2F_r_Ux_Ux);
                PRNT("        daux_diff_2K_Ux_Ux= %g\n",daux_diff_2K_Ux_Ux);
            }
            return(0);
        }
// aux_deriv[210]
    daux_diff_2Fxr_Ux_Uy=p[7][0]*aux_K*aux_rev*daux_diff_2F_r_Ux_Uy/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_2gamma_R_Ux_Uy/(1.0+aux_K)/pow(aux_gamma_R,2.0)-p[7][0]*aux_K*aux_rev*daux_diff_F_r_Uy*daux_diff_K_Ux/pow(1.0+aux_K,2.0)/aux_gamma_R+p[7][0]*aux_rev*daux_diff_F_r_Uy*daux_diff_K_Ux/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_K*aux_rev*daux_diff_F_r_Uy*daux_diff_gamma_R_Ux/(1.0+aux_K)/pow(aux_gamma_R,2.0)-p[7][0]*aux_K*aux_rev*daux_diff_F_r_Ux*daux_diff_gamma_R_Uy/(1.0+aux_K)/pow(aux_gamma_R,2.0)+p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_K_Ux*daux_diff_gamma_R_Uy/pow(1.0+aux_K,2.0)/pow(aux_gamma_R,2.0)-p[7][0]*aux_F_r*aux_rev*daux_diff_K_Ux*daux_diff_gamma_R_Uy/(1.0+aux_K)/pow(aux_gamma_R,2.0)+2.0*p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_gamma_R_Ux*daux_diff_gamma_R_Uy/(1.0+aux_K)/pow(aux_gamma_R,3.0)+p[7][0]*aux_K*daux_diff_F_r_Uy*daux_diff_rev_Ux/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*daux_diff_gamma_R_Uy*daux_diff_rev_Ux/(1.0+aux_K)/pow(aux_gamma_R,2.0);
    if (isNANorINF(daux_diff_2Fxr_Ux_Uy))
        {
            PRNT("    @k %d: daux_diff_2Fxr_Ux_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fxr_Ux_Uy);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        daux_diff_gamma_R_Uy= %g\n",daux_diff_gamma_R_Uy);
                PRNT("        daux_diff_rev_Ux= %g\n",daux_diff_rev_Ux);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_F_r_Uy= %g\n",daux_diff_F_r_Uy);
                PRNT("        daux_diff_gamma_R_Ux= %g\n",daux_diff_gamma_R_Ux);
                PRNT("        aux_rev= %g\n",aux_rev);
                PRNT("        daux_diff_2gamma_R_Ux_Uy= %g\n",daux_diff_2gamma_R_Ux_Uy);
                PRNT("        daux_diff_F_r_Ux= %g\n",daux_diff_F_r_Ux);
                PRNT("        daux_diff_K_Ux= %g\n",daux_diff_K_Ux);
                PRNT("        daux_diff_2F_r_Ux_Uy= %g\n",daux_diff_2F_r_Ux_Uy);
            }
            return(0);
        }
// aux_deriv[211]
    daux_diff_2Fxr_Ux_r=p[7][0]*aux_K*aux_rev*daux_diff_2F_r_Ux_r/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_2gamma_R_Ux_r/(1.0+aux_K)/pow(aux_gamma_R,2.0)-p[7][0]*aux_K*aux_rev*daux_diff_F_r_r*daux_diff_K_Ux/pow(1.0+aux_K,2.0)/aux_gamma_R+p[7][0]*aux_rev*daux_diff_F_r_r*daux_diff_K_Ux/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_K*aux_rev*daux_diff_F_r_r*daux_diff_gamma_R_Ux/(1.0+aux_K)/pow(aux_gamma_R,2.0)-p[7][0]*aux_K*aux_rev*daux_diff_F_r_Ux*daux_diff_gamma_R_r/(1.0+aux_K)/pow(aux_gamma_R,2.0)+p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_K_Ux*daux_diff_gamma_R_r/pow(1.0+aux_K,2.0)/pow(aux_gamma_R,2.0)-p[7][0]*aux_F_r*aux_rev*daux_diff_K_Ux*daux_diff_gamma_R_r/(1.0+aux_K)/pow(aux_gamma_R,2.0)+2.0*p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_gamma_R_Ux*daux_diff_gamma_R_r/(1.0+aux_K)/pow(aux_gamma_R,3.0)+p[7][0]*aux_K*daux_diff_F_r_r*daux_diff_rev_Ux/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*daux_diff_gamma_R_r*daux_diff_rev_Ux/(1.0+aux_K)/pow(aux_gamma_R,2.0);
    if (isNANorINF(daux_diff_2Fxr_Ux_r))
        {
            PRNT("    @k %d: daux_diff_2Fxr_Ux_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fxr_Ux_r);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        daux_diff_gamma_R_r= %g\n",daux_diff_gamma_R_r);
                PRNT("        daux_diff_rev_Ux= %g\n",daux_diff_rev_Ux);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_F_r_r= %g\n",daux_diff_F_r_r);
                PRNT("        daux_diff_gamma_R_Ux= %g\n",daux_diff_gamma_R_Ux);
                PRNT("        aux_rev= %g\n",aux_rev);
                PRNT("        daux_diff_2gamma_R_Ux_r= %g\n",daux_diff_2gamma_R_Ux_r);
                PRNT("        daux_diff_F_r_Ux= %g\n",daux_diff_F_r_Ux);
                PRNT("        daux_diff_K_Ux= %g\n",daux_diff_K_Ux);
                PRNT("        daux_diff_2F_r_Ux_r= %g\n",daux_diff_2F_r_Ux_r);
            }
            return(0);
        }
// aux_deriv[212]
    daux_diff_2Fxr_Ux_thr=p[7][0]*aux_K*aux_rev*daux_diff_2F_r_Ux_thr/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_2K_Ux_thr/pow(1.0+aux_K,2.0)/aux_gamma_R+p[7][0]*aux_F_r*aux_rev*daux_diff_2K_Ux_thr/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_2gamma_R_Ux_thr/(1.0+aux_K)/pow(aux_gamma_R,2.0)+p[7][0]*aux_F_r*aux_K*daux_diff_2rev_Ux_thr/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_K*aux_rev*daux_diff_F_r_thr*daux_diff_K_Ux/pow(1.0+aux_K,2.0)/aux_gamma_R+p[7][0]*aux_rev*daux_diff_F_r_thr*daux_diff_K_Ux/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_K*aux_rev*daux_diff_F_r_Ux*daux_diff_K_thr/pow(1.0+aux_K,2.0)/aux_gamma_R+p[7][0]*aux_rev*daux_diff_F_r_Ux*daux_diff_K_thr/(1.0+aux_K)/aux_gamma_R+2.0*p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_K_Ux*daux_diff_K_thr/pow(1.0+aux_K,3.0)/aux_gamma_R-2.0*p[7][0]*aux_F_r*aux_rev*daux_diff_K_Ux*daux_diff_K_thr/pow(1.0+aux_K,2.0)/aux_gamma_R-p[7][0]*aux_K*aux_rev*daux_diff_F_r_thr*daux_diff_gamma_R_Ux/(1.0+aux_K)/pow(aux_gamma_R,2.0)+p[
     7][0]*aux_F_r*aux_K*aux_rev*daux_diff_K_thr*daux_diff_gamma_R_Ux/pow(1.0+aux_K,2.0)/pow(aux_gamma_R,2.0)-p[7][0]*aux_F_r*aux_rev*daux_diff_K_thr*daux_diff_gamma_R_Ux/(1.0+aux_K)/pow(aux_gamma_R,2.0)-p[7][0]*aux_K*aux_rev*daux_diff_F_r_Ux*daux_diff_gamma_R_thr/(1.0+aux_K)/pow(aux_gamma_R,2.0)+p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_K_Ux*daux_diff_gamma_R_thr/pow(1.0+aux_K,2.0)/pow(aux_gamma_R,2.0)-p[7][0]*aux_F_r*aux_rev*daux_diff_K_Ux*daux_diff_gamma_R_thr/(1.0+aux_K)/pow(aux_gamma_R,2.0)+2.0*p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_gamma_R_Ux*daux_diff_gamma_R_thr/(1.0+aux_K)/pow(aux_gamma_R,3.0)+p[7][0]*aux_K*daux_diff_F_r_thr*daux_diff_rev_Ux/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*daux_diff_K_thr*daux_diff_rev_Ux/pow(1.0+aux_K,2.0)/aux_gamma_R+p[7][0]*aux_F_r*daux_diff_K_thr*daux_diff_rev_Ux/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*daux_diff_gamma_R_thr*daux_diff_rev_Ux/(1.0+aux_K)/pow(aux_gamma_R,2.0)+p[7][0]*aux_K*daux_diff_F_r_Ux*daux_diff_rev_thr/(1.0+aux_K)/
     aux_gamma_R-p[7][0]*aux_F_r*aux_K*daux_diff_K_Ux*daux_diff_rev_thr/pow(1.0+aux_K,2.0)/aux_gamma_R+p[7][0]*aux_F_r*daux_diff_K_Ux*daux_diff_rev_thr/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*daux_diff_gamma_R_Ux*daux_diff_rev_thr/(1.0+aux_K)/pow(aux_gamma_R,2.0);
    if (isNANorINF(daux_diff_2Fxr_Ux_thr))
        {
            PRNT("    @k %d: daux_diff_2Fxr_Ux_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fxr_Ux_thr);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        daux_diff_gamma_R_thr= %g\n",daux_diff_gamma_R_thr);
                PRNT("        daux_diff_rev_Ux= %g\n",daux_diff_rev_Ux);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_gamma_R_Ux= %g\n",daux_diff_gamma_R_Ux);
                PRNT("        daux_diff_rev_thr= %g\n",daux_diff_rev_thr);
                PRNT("        daux_diff_2rev_Ux_thr= %g\n",daux_diff_2rev_Ux_thr);
                PRNT("        daux_diff_F_r_thr= %g\n",daux_diff_F_r_thr);
                PRNT("        daux_diff_K_thr= %g\n",daux_diff_K_thr);
                PRNT("        daux_diff_F_r_Ux= %g\n",daux_diff_F_r_Ux);
                PRNT("        daux_diff_K_Ux= %g\n",daux_diff_K_Ux);
                PRNT("        aux_rev= %g\n",aux_rev);
                PRNT("        daux_diff_2gamma_R_Ux_thr= %g\n",daux_diff_2gamma_R_Ux_thr);
                PRNT("        daux_diff_2F_r_Ux_thr= %g\n",daux_diff_2F_r_Ux_thr);
                PRNT("        daux_diff_2K_Ux_thr= %g\n",daux_diff_2K_Ux_thr);
            }
            return(0);
        }
// aux_deriv[213]
    daux_diff_Fxr_Uy=p[7][0]*aux_K*aux_rev*daux_diff_F_r_Uy/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_gamma_R_Uy/(1.0+aux_K)/pow(aux_gamma_R,2.0);
    if (isNANorINF(daux_diff_Fxr_Uy))
        {
            PRNT("    @k %d: daux_diff_Fxr_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_Fxr_Uy);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        daux_diff_gamma_R_Uy= %g\n",daux_diff_gamma_R_Uy);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        aux_rev= %g\n",aux_rev);
                PRNT("        daux_diff_F_r_Uy= %g\n",daux_diff_F_r_Uy);
            }
            return(0);
        }
// aux_deriv[214]
    daux_diff_2Fxr_Uy_Uy=p[7][0]*aux_K*aux_rev*daux_diff_2F_r_Uy_Uy/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_2gamma_R_Uy_Uy/(1.0+aux_K)/pow(aux_gamma_R,2.0)-2.0*p[7][0]*aux_K*aux_rev*daux_diff_F_r_Uy*daux_diff_gamma_R_Uy/(1.0+aux_K)/pow(aux_gamma_R,2.0)+2.0*p[7][0]*aux_F_r*aux_K*aux_rev*pow(daux_diff_gamma_R_Uy,2.0)/(1.0+aux_K)/pow(aux_gamma_R,3.0);
    if (isNANorINF(daux_diff_2Fxr_Uy_Uy))
        {
            PRNT("    @k %d: daux_diff_2Fxr_Uy_Uy in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fxr_Uy_Uy);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        daux_diff_gamma_R_Uy= %g\n",daux_diff_gamma_R_Uy);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        aux_rev= %g\n",aux_rev);
                PRNT("        daux_diff_2gamma_R_Uy_Uy= %g\n",daux_diff_2gamma_R_Uy_Uy);
                PRNT("        daux_diff_F_r_Uy= %g\n",daux_diff_F_r_Uy);
                PRNT("        daux_diff_2F_r_Uy_Uy= %g\n",daux_diff_2F_r_Uy_Uy);
            }
            return(0);
        }
// aux_deriv[215]
    daux_diff_2Fxr_Uy_r=p[7][0]*aux_K*aux_rev*daux_diff_2F_r_Uy_r/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_2gamma_R_Uy_r/(1.0+aux_K)/pow(aux_gamma_R,2.0)-p[7][0]*aux_K*aux_rev*daux_diff_F_r_r*daux_diff_gamma_R_Uy/(1.0+aux_K)/pow(aux_gamma_R,2.0)-p[7][0]*aux_K*aux_rev*daux_diff_F_r_Uy*daux_diff_gamma_R_r/(1.0+aux_K)/pow(aux_gamma_R,2.0)+2.0*p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_gamma_R_Uy*daux_diff_gamma_R_r/(1.0+aux_K)/pow(aux_gamma_R,3.0);
    if (isNANorINF(daux_diff_2Fxr_Uy_r))
        {
            PRNT("    @k %d: daux_diff_2Fxr_Uy_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fxr_Uy_r);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        daux_diff_gamma_R_Uy= %g\n",daux_diff_gamma_R_Uy);
                PRNT("        daux_diff_gamma_R_r= %g\n",daux_diff_gamma_R_r);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        aux_rev= %g\n",aux_rev);
                PRNT("        daux_diff_2gamma_R_Uy_r= %g\n",daux_diff_2gamma_R_Uy_r);
                PRNT("        daux_diff_F_r_r= %g\n",daux_diff_F_r_r);
                PRNT("        daux_diff_F_r_Uy= %g\n",daux_diff_F_r_Uy);
                PRNT("        daux_diff_2F_r_Uy_r= %g\n",daux_diff_2F_r_Uy_r);
            }
            return(0);
        }
// aux_deriv[216]
    daux_diff_2Fxr_Uy_thr=p[7][0]*aux_K*aux_rev*daux_diff_2F_r_Uy_thr/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_2gamma_R_Uy_thr/(1.0+aux_K)/pow(aux_gamma_R,2.0)-p[7][0]*aux_K*aux_rev*daux_diff_F_r_Uy*daux_diff_K_thr/pow(1.0+aux_K,2.0)/aux_gamma_R+p[7][0]*aux_rev*daux_diff_F_r_Uy*daux_diff_K_thr/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_K*aux_rev*daux_diff_F_r_thr*daux_diff_gamma_R_Uy/(1.0+aux_K)/pow(aux_gamma_R,2.0)+p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_K_thr*daux_diff_gamma_R_Uy/pow(1.0+aux_K,2.0)/pow(aux_gamma_R,2.0)-p[7][0]*aux_F_r*aux_rev*daux_diff_K_thr*daux_diff_gamma_R_Uy/(1.0+aux_K)/pow(aux_gamma_R,2.0)-p[7][0]*aux_K*aux_rev*daux_diff_F_r_Uy*daux_diff_gamma_R_thr/(1.0+aux_K)/pow(aux_gamma_R,2.0)+2.0*p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_gamma_R_Uy*daux_diff_gamma_R_thr/(1.0+aux_K)/pow(aux_gamma_R,3.0)+p[7][0]*aux_K*daux_diff_F_r_Uy*daux_diff_rev_thr/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*daux_diff_gamma_R_Uy*daux_diff_rev_thr/(1.0+aux_K)/pow(
     aux_gamma_R,2.0);
    if (isNANorINF(daux_diff_2Fxr_Uy_thr))
        {
            PRNT("    @k %d: daux_diff_2Fxr_Uy_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fxr_Uy_thr);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        daux_diff_gamma_R_Uy= %g\n",daux_diff_gamma_R_Uy);
                PRNT("        daux_diff_rev_thr= %g\n",daux_diff_rev_thr);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_F_r_Uy= %g\n",daux_diff_F_r_Uy);
                PRNT("        daux_diff_gamma_R_thr= %g\n",daux_diff_gamma_R_thr);
                PRNT("        aux_rev= %g\n",aux_rev);
                PRNT("        daux_diff_2gamma_R_Uy_thr= %g\n",daux_diff_2gamma_R_Uy_thr);
                PRNT("        daux_diff_F_r_thr= %g\n",daux_diff_F_r_thr);
                PRNT("        daux_diff_K_thr= %g\n",daux_diff_K_thr);
                PRNT("        daux_diff_2F_r_Uy_thr= %g\n",daux_diff_2F_r_Uy_thr);
            }
            return(0);
        }
// aux_deriv[217]
    daux_diff_Fxr_r=p[7][0]*aux_K*aux_rev*daux_diff_F_r_r/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_gamma_R_r/(1.0+aux_K)/pow(aux_gamma_R,2.0);
    if (isNANorINF(daux_diff_Fxr_r))
        {
            PRNT("    @k %d: daux_diff_Fxr_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_Fxr_r);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        daux_diff_gamma_R_r= %g\n",daux_diff_gamma_R_r);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        aux_rev= %g\n",aux_rev);
                PRNT("        daux_diff_F_r_r= %g\n",daux_diff_F_r_r);
            }
            return(0);
        }
// aux_deriv[218]
    daux_diff_2Fxr_r_r=p[7][0]*aux_K*aux_rev*daux_diff_2F_r_r_r/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_2gamma_R_r_r/(1.0+aux_K)/pow(aux_gamma_R,2.0)-2.0*p[7][0]*aux_K*aux_rev*daux_diff_F_r_r*daux_diff_gamma_R_r/(1.0+aux_K)/pow(aux_gamma_R,2.0)+2.0*p[7][0]*aux_F_r*aux_K*aux_rev*pow(daux_diff_gamma_R_r,2.0)/(1.0+aux_K)/pow(aux_gamma_R,3.0);
    if (isNANorINF(daux_diff_2Fxr_r_r))
        {
            PRNT("    @k %d: daux_diff_2Fxr_r_r in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fxr_r_r);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        daux_diff_gamma_R_r= %g\n",daux_diff_gamma_R_r);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        aux_rev= %g\n",aux_rev);
                PRNT("        daux_diff_2gamma_R_r_r= %g\n",daux_diff_2gamma_R_r_r);
                PRNT("        daux_diff_F_r_r= %g\n",daux_diff_F_r_r);
                PRNT("        daux_diff_2F_r_r_r= %g\n",daux_diff_2F_r_r_r);
            }
            return(0);
        }
// aux_deriv[219]
    daux_diff_2Fxr_r_thr=p[7][0]*aux_K*aux_rev*daux_diff_2F_r_r_thr/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_2gamma_R_r_thr/(1.0+aux_K)/pow(aux_gamma_R,2.0)-p[7][0]*aux_K*aux_rev*daux_diff_F_r_r*daux_diff_K_thr/pow(1.0+aux_K,2.0)/aux_gamma_R+p[7][0]*aux_rev*daux_diff_F_r_r*daux_diff_K_thr/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_K*aux_rev*daux_diff_F_r_thr*daux_diff_gamma_R_r/(1.0+aux_K)/pow(aux_gamma_R,2.0)+p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_K_thr*daux_diff_gamma_R_r/pow(1.0+aux_K,2.0)/pow(aux_gamma_R,2.0)-p[7][0]*aux_F_r*aux_rev*daux_diff_K_thr*daux_diff_gamma_R_r/(1.0+aux_K)/pow(aux_gamma_R,2.0)-p[7][0]*aux_K*aux_rev*daux_diff_F_r_r*daux_diff_gamma_R_thr/(1.0+aux_K)/pow(aux_gamma_R,2.0)+2.0*p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_gamma_R_r*daux_diff_gamma_R_thr/(1.0+aux_K)/pow(aux_gamma_R,3.0)+p[7][0]*aux_K*daux_diff_F_r_r*daux_diff_rev_thr/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*daux_diff_gamma_R_r*daux_diff_rev_thr/(1.0+aux_K)/pow(aux_gamma_R,2.0);
    if (isNANorINF(daux_diff_2Fxr_r_thr))
        {
            PRNT("    @k %d: daux_diff_2Fxr_r_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fxr_r_thr);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        daux_diff_gamma_R_r= %g\n",daux_diff_gamma_R_r);
                PRNT("        daux_diff_rev_thr= %g\n",daux_diff_rev_thr);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_F_r_r= %g\n",daux_diff_F_r_r);
                PRNT("        daux_diff_gamma_R_thr= %g\n",daux_diff_gamma_R_thr);
                PRNT("        aux_rev= %g\n",aux_rev);
                PRNT("        daux_diff_2gamma_R_r_thr= %g\n",daux_diff_2gamma_R_r_thr);
                PRNT("        daux_diff_F_r_thr= %g\n",daux_diff_F_r_thr);
                PRNT("        daux_diff_K_thr= %g\n",daux_diff_K_thr);
                PRNT("        daux_diff_2F_r_r_thr= %g\n",daux_diff_2F_r_r_thr);
            }
            return(0);
        }
// aux_deriv[220]
    daux_diff_Fxr_thr=p[7][0]*aux_K*aux_rev*daux_diff_F_r_thr/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_K_thr/pow(1.0+aux_K,2.0)/aux_gamma_R+p[7][0]*aux_F_r*aux_rev*daux_diff_K_thr/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_gamma_R_thr/(1.0+aux_K)/pow(aux_gamma_R,2.0)+p[7][0]*aux_F_r*aux_K*daux_diff_rev_thr/(1.0+aux_K)/aux_gamma_R;
    if (isNANorINF(daux_diff_Fxr_thr))
        {
            PRNT("    @k %d: daux_diff_Fxr_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_Fxr_thr);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        daux_diff_rev_thr= %g\n",daux_diff_rev_thr);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_gamma_R_thr= %g\n",daux_diff_gamma_R_thr);
                PRNT("        aux_rev= %g\n",aux_rev);
                PRNT("        daux_diff_F_r_thr= %g\n",daux_diff_F_r_thr);
                PRNT("        daux_diff_K_thr= %g\n",daux_diff_K_thr);
            }
            return(0);
        }
// aux_deriv[221]
    daux_diff_2Fxr_thr_thr=p[7][0]*aux_K*aux_rev*daux_diff_2F_r_thr_thr/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_2K_thr_thr/pow(1.0+aux_K,2.0)/aux_gamma_R+p[7][0]*aux_F_r*aux_rev*daux_diff_2K_thr_thr/(1.0+aux_K)/aux_gamma_R-p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_2gamma_R_thr_thr/(1.0+aux_K)/pow(aux_gamma_R,2.0)+p[7][0]*aux_F_r*aux_K*daux_diff_2rev_thr_thr/(1.0+aux_K)/aux_gamma_R-2.0*p[7][0]*aux_K*aux_rev*daux_diff_F_r_thr*daux_diff_K_thr/pow(1.0+aux_K,2.0)/aux_gamma_R+2.0*p[7][0]*aux_rev*daux_diff_F_r_thr*daux_diff_K_thr/(1.0+aux_K)/aux_gamma_R+2.0*p[7][0]*aux_F_r*aux_K*aux_rev*pow(daux_diff_K_thr,2.0)/pow(1.0+aux_K,3.0)/aux_gamma_R-2.0*p[7][0]*aux_F_r*aux_rev*pow(daux_diff_K_thr,2.0)/pow(1.0+aux_K,2.0)/aux_gamma_R-2.0*p[7][0]*aux_K*aux_rev*daux_diff_F_r_thr*daux_diff_gamma_R_thr/(1.0+aux_K)/pow(aux_gamma_R,2.0)+2.0*p[7][0]*aux_F_r*aux_K*aux_rev*daux_diff_K_thr*daux_diff_gamma_R_thr/pow(1.0+aux_K,2.0)/pow(aux_gamma_R,2.0)-2.0*p[7][0]*aux_F_r*aux_rev*
     daux_diff_K_thr*daux_diff_gamma_R_thr/(1.0+aux_K)/pow(aux_gamma_R,2.0)+2.0*p[7][0]*aux_F_r*aux_K*aux_rev*pow(daux_diff_gamma_R_thr,2.0)/(1.0+aux_K)/pow(aux_gamma_R,3.0)+2.0*p[7][0]*aux_K*daux_diff_F_r_thr*daux_diff_rev_thr/(1.0+aux_K)/aux_gamma_R-2.0*p[7][0]*aux_F_r*aux_K*daux_diff_K_thr*daux_diff_rev_thr/pow(1.0+aux_K,2.0)/aux_gamma_R+2.0*p[7][0]*aux_F_r*daux_diff_K_thr*daux_diff_rev_thr/(1.0+aux_K)/aux_gamma_R-2.0*p[7][0]*aux_F_r*aux_K*daux_diff_gamma_R_thr*daux_diff_rev_thr/(1.0+aux_K)/pow(aux_gamma_R,2.0);
    if (isNANorINF(daux_diff_2Fxr_thr_thr))
        {
            PRNT("    @k %d: daux_diff_2Fxr_thr_thr in line %d is nan or inf: %g\n", k, __LINE__-3,daux_diff_2Fxr_thr_thr);
            {
                PRNT("        aux_F_r= %g\n",aux_F_r);
                PRNT("        aux_K= %g\n",aux_K);
                PRNT("        p[7,0]= %g\n",p[7][0]);
                PRNT("        daux_diff_gamma_R_thr= %g\n",daux_diff_gamma_R_thr);
                PRNT("        daux_diff_rev_thr= %g\n",daux_diff_rev_thr);
                PRNT("        aux_gamma_R= %g\n",aux_gamma_R);
                PRNT("        daux_diff_2rev_thr_thr= %g\n",daux_diff_2rev_thr_thr);
                PRNT("        daux_diff_F_r_thr= %g\n",daux_diff_F_r_thr);
                PRNT("        daux_diff_K_thr= %g\n",daux_diff_K_thr);
                PRNT("        aux_rev= %g\n",aux_rev);
                PRNT("        daux_diff_2gamma_R_thr_thr= %g\n",daux_diff_2gamma_R_thr_thr);
                PRNT("        daux_diff_2F_r_thr_thr= %g\n",daux_diff_2F_r_thr_thr);
                PRNT("        daux_diff_2K_thr_thr= %g\n",daux_diff_2K_thr_thr);
            }
            return(0);
        }
    return 1;
}

static int bp_derivsL(trajEl_t *t, int k, double **p) {
    const double *x= t->x;
    const double *u= t->u;
    
// derivatives of f
// df[i]/d x[0]
// df[i]/d x[1]
// df[i]/d x[2]
    t->fx[20]=p[15][0]*daux_diff_dx_phi;
    if (isNANorINF(t->fx[20]))
        {
            PRNT("    @k %d: t->fx[20] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fx[20]);
            {
                PRNT("        daux_diff_dx_phi= %g\n",daux_diff_dx_phi);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fx[21]=p[15][0]*daux_diff_dy_phi;
    if (isNANorINF(t->fx[21]))
        {
            PRNT("    @k %d: t->fx[21] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fx[21]);
            {
                PRNT("        daux_diff_dy_phi= %g\n",daux_diff_dy_phi);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// df[i]/d x[3]
    t->fx[30]=p[15][0]*daux_diff_dx_Ux;
    if (isNANorINF(t->fx[30]))
        {
            PRNT("    @k %d: t->fx[30] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fx[30]);
            {
                PRNT("        daux_diff_dx_Ux= %g\n",daux_diff_dx_Ux);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fx[31]=p[15][0]*daux_diff_dy_Ux;
    if (isNANorINF(t->fx[31]))
        {
            PRNT("    @k %d: t->fx[31] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fx[31]);
            {
                PRNT("        daux_diff_dy_Ux= %g\n",daux_diff_dy_Ux);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fx[33]=1.0+p[15][0]*daux_diff_dUx_Ux;
    if (isNANorINF(t->fx[33]))
        {
            PRNT("    @k %d: t->fx[33] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fx[33]);
            {
                PRNT("        daux_diff_dUx_Ux= %g\n",daux_diff_dUx_Ux);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fx[34]=p[15][0]*daux_diff_dUy_Ux;
    if (isNANorINF(t->fx[34]))
        {
            PRNT("    @k %d: t->fx[34] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fx[34]);
            {
                PRNT("        daux_diff_dUy_Ux= %g\n",daux_diff_dUy_Ux);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fx[35]=p[15][0]*daux_diff_dr_Ux;
    if (isNANorINF(t->fx[35]))
        {
            PRNT("    @k %d: t->fx[35] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fx[35]);
            {
                PRNT("        daux_diff_dr_Ux= %g\n",daux_diff_dr_Ux);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// df[i]/d x[4]
    t->fx[40]=p[15][0]*daux_diff_dx_Uy;
    if (isNANorINF(t->fx[40]))
        {
            PRNT("    @k %d: t->fx[40] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fx[40]);
            {
                PRNT("        daux_diff_dx_Uy= %g\n",daux_diff_dx_Uy);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fx[41]=p[15][0]*daux_diff_dy_Uy;
    if (isNANorINF(t->fx[41]))
        {
            PRNT("    @k %d: t->fx[41] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fx[41]);
            {
                PRNT("        daux_diff_dy_Uy= %g\n",daux_diff_dy_Uy);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fx[43]=p[15][0]*daux_diff_dUx_Uy;
    if (isNANorINF(t->fx[43]))
        {
            PRNT("    @k %d: t->fx[43] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fx[43]);
            {
                PRNT("        daux_diff_dUx_Uy= %g\n",daux_diff_dUx_Uy);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fx[44]=1.0+p[15][0]*daux_diff_dUy_Uy;
    if (isNANorINF(t->fx[44]))
        {
            PRNT("    @k %d: t->fx[44] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fx[44]);
            {
                PRNT("        daux_diff_dUy_Uy= %g\n",daux_diff_dUy_Uy);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fx[45]=p[15][0]*daux_diff_dr_Uy;
    if (isNANorINF(t->fx[45]))
        {
            PRNT("    @k %d: t->fx[45] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fx[45]);
            {
                PRNT("        daux_diff_dr_Uy= %g\n",daux_diff_dr_Uy);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// df[i]/d x[5]
    t->fx[53]=p[15][0]*daux_diff_dUx_r;
    if (isNANorINF(t->fx[53]))
        {
            PRNT("    @k %d: t->fx[53] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fx[53]);
            {
                PRNT("        daux_diff_dUx_r= %g\n",daux_diff_dUx_r);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fx[54]=p[15][0]*daux_diff_dUy_r;
    if (isNANorINF(t->fx[54]))
        {
            PRNT("    @k %d: t->fx[54] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fx[54]);
            {
                PRNT("        daux_diff_dUy_r= %g\n",daux_diff_dUy_r);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fx[55]=1.0+p[15][0]*daux_diff_dr_r;
    if (isNANorINF(t->fx[55]))
        {
            PRNT("    @k %d: t->fx[55] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fx[55]);
            {
                PRNT("        daux_diff_dr_r= %g\n",daux_diff_dr_r);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// df[i]/d x[6]
// df[i]/d x[7]
// df[i]/d x[8]
// df[i]/d x[9]

// df[i]/d u[0]
    t->fu[3]=p[15][0]*daux_diff_dUx_thr;
    if (isNANorINF(t->fu[3]))
        {
            PRNT("    @k %d: t->fu[3] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fu[3]);
            {
                PRNT("        daux_diff_dUx_thr= %g\n",daux_diff_dUx_thr);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fu[4]=p[15][0]*daux_diff_dUy_thr;
    if (isNANorINF(t->fu[4]))
        {
            PRNT("    @k %d: t->fu[4] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fu[4]);
            {
                PRNT("        daux_diff_dUy_thr= %g\n",daux_diff_dUy_thr);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fu[5]=p[15][0]*daux_diff_dr_thr;
    if (isNANorINF(t->fu[5]))
        {
            PRNT("    @k %d: t->fu[5] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fu[5]);
            {
                PRNT("        daux_diff_dr_thr= %g\n",daux_diff_dr_thr);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// df[i]/d u[1]
    t->fu[13]=p[15][0]*daux_diff_dUx_steer;
    if (isNANorINF(t->fu[13]))
        {
            PRNT("    @k %d: t->fu[13] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fu[13]);
            {
                PRNT("        daux_diff_dUx_steer= %g\n",daux_diff_dUx_steer);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fu[14]=p[15][0]*daux_diff_dUy_steer;
    if (isNANorINF(t->fu[14]))
        {
            PRNT("    @k %d: t->fu[14] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fu[14]);
            {
                PRNT("        daux_diff_dUy_steer= %g\n",daux_diff_dUy_steer);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fu[15]=p[15][0]*daux_diff_dr_steer;
    if (isNANorINF(t->fu[15]))
        {
            PRNT("    @k %d: t->fu[15] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fu[15]);
            {
                PRNT("        daux_diff_dr_steer= %g\n",daux_diff_dr_steer);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }

#if FULL_DDP
// d^2f[0]/dx[i] dx[j]
// j= 0
// j= 1
// j= 2
    t->fxx[5]=p[15][0]*daux_diff_2dx_phi_phi;
    if (isNANorINF(t->fxx[5]))
        {
            PRNT("    @k %d: t->fxx[5] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[5]);
            {
                PRNT("        daux_diff_2dx_phi_phi= %g\n",daux_diff_2dx_phi_phi);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 3
    t->fxx[9]=p[15][0]*daux_diff_2dx_Ux_Ux;
    if (isNANorINF(t->fxx[9]))
        {
            PRNT("    @k %d: t->fxx[9] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[9]);
            {
                PRNT("        daux_diff_2dx_Ux_Ux= %g\n",daux_diff_2dx_Ux_Ux);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 4
    t->fxx[13]=p[15][0]*daux_diff_2dx_Ux_Uy;
    if (isNANorINF(t->fxx[13]))
        {
            PRNT("    @k %d: t->fxx[13] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[13]);
            {
                PRNT("        daux_diff_2dx_Ux_Uy= %g\n",daux_diff_2dx_Ux_Uy);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxx[14]=p[15][0]*daux_diff_2dx_Uy_Uy;
    if (isNANorINF(t->fxx[14]))
        {
            PRNT("    @k %d: t->fxx[14] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[14]);
            {
                PRNT("        daux_diff_2dx_Uy_Uy= %g\n",daux_diff_2dx_Uy_Uy);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 5
// j= 6
// j= 7
// j= 8
// j= 9
// d^2f[1]/dx[i] dx[j]
// j= 0
// j= 1
// j= 2
    t->fxx[60]=p[15][0]*daux_diff_2dy_phi_phi;
    if (isNANorINF(t->fxx[60]))
        {
            PRNT("    @k %d: t->fxx[60] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[60]);
            {
                PRNT("        daux_diff_2dy_phi_phi= %g\n",daux_diff_2dy_phi_phi);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 3
    t->fxx[64]=p[15][0]*daux_diff_2dy_Ux_Ux;
    if (isNANorINF(t->fxx[64]))
        {
            PRNT("    @k %d: t->fxx[64] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[64]);
            {
                PRNT("        daux_diff_2dy_Ux_Ux= %g\n",daux_diff_2dy_Ux_Ux);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 4
    t->fxx[68]=p[15][0]*daux_diff_2dy_Ux_Uy;
    if (isNANorINF(t->fxx[68]))
        {
            PRNT("    @k %d: t->fxx[68] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[68]);
            {
                PRNT("        daux_diff_2dy_Ux_Uy= %g\n",daux_diff_2dy_Ux_Uy);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxx[69]=p[15][0]*daux_diff_2dy_Uy_Uy;
    if (isNANorINF(t->fxx[69]))
        {
            PRNT("    @k %d: t->fxx[69] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[69]);
            {
                PRNT("        daux_diff_2dy_Uy_Uy= %g\n",daux_diff_2dy_Uy_Uy);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 5
// j= 6
// j= 7
// j= 8
// j= 9
// d^2f[2]/dx[i] dx[j]
// j= 0
// j= 1
// j= 2
// j= 3
// j= 4
// j= 5
// j= 6
// j= 7
// j= 8
// j= 9
// d^2f[3]/dx[i] dx[j]
// j= 0
// j= 1
// j= 2
// j= 3
    t->fxx[174]=p[15][0]*daux_diff_2dUx_Ux_Ux;
    if (isNANorINF(t->fxx[174]))
        {
            PRNT("    @k %d: t->fxx[174] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[174]);
            {
                PRNT("        daux_diff_2dUx_Ux_Ux= %g\n",daux_diff_2dUx_Ux_Ux);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 4
    t->fxx[178]=p[15][0]*daux_diff_2dUx_Ux_Uy;
    if (isNANorINF(t->fxx[178]))
        {
            PRNT("    @k %d: t->fxx[178] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[178]);
            {
                PRNT("        daux_diff_2dUx_Ux_Uy= %g\n",daux_diff_2dUx_Ux_Uy);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxx[179]=p[15][0]*daux_diff_2dUx_Uy_Uy;
    if (isNANorINF(t->fxx[179]))
        {
            PRNT("    @k %d: t->fxx[179] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[179]);
            {
                PRNT("        daux_diff_2dUx_Uy_Uy= %g\n",daux_diff_2dUx_Uy_Uy);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 5
    t->fxx[183]=p[15][0]*daux_diff_2dUx_Ux_r;
    if (isNANorINF(t->fxx[183]))
        {
            PRNT("    @k %d: t->fxx[183] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[183]);
            {
                PRNT("        daux_diff_2dUx_Ux_r= %g\n",daux_diff_2dUx_Ux_r);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxx[184]=p[15][0]*daux_diff_2dUx_Uy_r;
    if (isNANorINF(t->fxx[184]))
        {
            PRNT("    @k %d: t->fxx[184] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[184]);
            {
                PRNT("        daux_diff_2dUx_Uy_r= %g\n",daux_diff_2dUx_Uy_r);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxx[185]=p[15][0]*daux_diff_2dUx_r_r;
    if (isNANorINF(t->fxx[185]))
        {
            PRNT("    @k %d: t->fxx[185] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[185]);
            {
                PRNT("        daux_diff_2dUx_r_r= %g\n",daux_diff_2dUx_r_r);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 6
// j= 7
// j= 8
// j= 9
// d^2f[4]/dx[i] dx[j]
// j= 0
// j= 1
// j= 2
// j= 3
    t->fxx[229]=p[15][0]*daux_diff_2dUy_Ux_Ux;
    if (isNANorINF(t->fxx[229]))
        {
            PRNT("    @k %d: t->fxx[229] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[229]);
            {
                PRNT("        daux_diff_2dUy_Ux_Ux= %g\n",daux_diff_2dUy_Ux_Ux);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 4
    t->fxx[233]=p[15][0]*daux_diff_2dUy_Ux_Uy;
    if (isNANorINF(t->fxx[233]))
        {
            PRNT("    @k %d: t->fxx[233] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[233]);
            {
                PRNT("        daux_diff_2dUy_Ux_Uy= %g\n",daux_diff_2dUy_Ux_Uy);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxx[234]=p[15][0]*daux_diff_2dUy_Uy_Uy;
    if (isNANorINF(t->fxx[234]))
        {
            PRNT("    @k %d: t->fxx[234] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[234]);
            {
                PRNT("        daux_diff_2dUy_Uy_Uy= %g\n",daux_diff_2dUy_Uy_Uy);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 5
    t->fxx[238]=p[15][0]*daux_diff_2dUy_Ux_r;
    if (isNANorINF(t->fxx[238]))
        {
            PRNT("    @k %d: t->fxx[238] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[238]);
            {
                PRNT("        daux_diff_2dUy_Ux_r= %g\n",daux_diff_2dUy_Ux_r);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxx[239]=p[15][0]*daux_diff_2dUy_Uy_r;
    if (isNANorINF(t->fxx[239]))
        {
            PRNT("    @k %d: t->fxx[239] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[239]);
            {
                PRNT("        daux_diff_2dUy_Uy_r= %g\n",daux_diff_2dUy_Uy_r);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxx[240]=p[15][0]*daux_diff_2dUy_r_r;
    if (isNANorINF(t->fxx[240]))
        {
            PRNT("    @k %d: t->fxx[240] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[240]);
            {
                PRNT("        daux_diff_2dUy_r_r= %g\n",daux_diff_2dUy_r_r);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 6
// j= 7
// j= 8
// j= 9
// d^2f[5]/dx[i] dx[j]
// j= 0
// j= 1
// j= 2
// j= 3
    t->fxx[284]=p[15][0]*daux_diff_2dr_Ux_Ux;
    if (isNANorINF(t->fxx[284]))
        {
            PRNT("    @k %d: t->fxx[284] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[284]);
            {
                PRNT("        daux_diff_2dr_Ux_Ux= %g\n",daux_diff_2dr_Ux_Ux);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 4
    t->fxx[288]=p[15][0]*daux_diff_2dr_Ux_Uy;
    if (isNANorINF(t->fxx[288]))
        {
            PRNT("    @k %d: t->fxx[288] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[288]);
            {
                PRNT("        daux_diff_2dr_Ux_Uy= %g\n",daux_diff_2dr_Ux_Uy);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxx[289]=p[15][0]*daux_diff_2dr_Uy_Uy;
    if (isNANorINF(t->fxx[289]))
        {
            PRNT("    @k %d: t->fxx[289] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[289]);
            {
                PRNT("        daux_diff_2dr_Uy_Uy= %g\n",daux_diff_2dr_Uy_Uy);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 5
    t->fxx[293]=p[15][0]*daux_diff_2dr_Ux_r;
    if (isNANorINF(t->fxx[293]))
        {
            PRNT("    @k %d: t->fxx[293] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[293]);
            {
                PRNT("        daux_diff_2dr_Ux_r= %g\n",daux_diff_2dr_Ux_r);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxx[294]=p[15][0]*daux_diff_2dr_Uy_r;
    if (isNANorINF(t->fxx[294]))
        {
            PRNT("    @k %d: t->fxx[294] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[294]);
            {
                PRNT("        daux_diff_2dr_Uy_r= %g\n",daux_diff_2dr_Uy_r);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxx[295]=p[15][0]*daux_diff_2dr_r_r;
    if (isNANorINF(t->fxx[295]))
        {
            PRNT("    @k %d: t->fxx[295] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxx[295]);
            {
                PRNT("        daux_diff_2dr_r_r= %g\n",daux_diff_2dr_r_r);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 6
// j= 7
// j= 8
// j= 9
// d^2f[6]/dx[i] dx[j]
// j= 0
// j= 1
// j= 2
// j= 3
// j= 4
// j= 5
// j= 6
// j= 7
// j= 8
// j= 9
// d^2f[7]/dx[i] dx[j]
// j= 0
// j= 1
// j= 2
// j= 3
// j= 4
// j= 5
// j= 6
// j= 7
// j= 8
// j= 9
// d^2f[8]/dx[i] dx[j]
// j= 0
// j= 1
// j= 2
// j= 3
// j= 4
// j= 5
// j= 6
// j= 7
// j= 8
// j= 9
// d^2f[9]/dx[i] dx[j]
// j= 0
// j= 1
// j= 2
// j= 3
// j= 4
// j= 5
// j= 6
// j= 7
// j= 8
// j= 9

// d^2f[0]/du[i] du[j]
// j= 0
// j= 1
// d^2f[1]/du[i] du[j]
// j= 0
// j= 1
// d^2f[2]/du[i] du[j]
// j= 0
// j= 1
// d^2f[3]/du[i] du[j]
// j= 0
    t->fuu[9]=p[15][0]*daux_diff_2dUx_thr_thr;
    if (isNANorINF(t->fuu[9]))
        {
            PRNT("    @k %d: t->fuu[9] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fuu[9]);
            {
                PRNT("        daux_diff_2dUx_thr_thr= %g\n",daux_diff_2dUx_thr_thr);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 1
    t->fuu[11]=p[15][0]*daux_diff_2dUx_steer_steer;
    if (isNANorINF(t->fuu[11]))
        {
            PRNT("    @k %d: t->fuu[11] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fuu[11]);
            {
                PRNT("        daux_diff_2dUx_steer_steer= %g\n",daux_diff_2dUx_steer_steer);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// d^2f[4]/du[i] du[j]
// j= 0
    t->fuu[12]=p[15][0]*daux_diff_2dUy_thr_thr;
    if (isNANorINF(t->fuu[12]))
        {
            PRNT("    @k %d: t->fuu[12] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fuu[12]);
            {
                PRNT("        daux_diff_2dUy_thr_thr= %g\n",daux_diff_2dUy_thr_thr);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 1
    t->fuu[14]=p[15][0]*daux_diff_2dUy_steer_steer;
    if (isNANorINF(t->fuu[14]))
        {
            PRNT("    @k %d: t->fuu[14] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fuu[14]);
            {
                PRNT("        daux_diff_2dUy_steer_steer= %g\n",daux_diff_2dUy_steer_steer);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// d^2f[5]/du[i] du[j]
// j= 0
    t->fuu[15]=p[15][0]*daux_diff_2dr_thr_thr;
    if (isNANorINF(t->fuu[15]))
        {
            PRNT("    @k %d: t->fuu[15] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fuu[15]);
            {
                PRNT("        daux_diff_2dr_thr_thr= %g\n",daux_diff_2dr_thr_thr);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 1
    t->fuu[17]=p[15][0]*daux_diff_2dr_steer_steer;
    if (isNANorINF(t->fuu[17]))
        {
            PRNT("    @k %d: t->fuu[17] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fuu[17]);
            {
                PRNT("        daux_diff_2dr_steer_steer= %g\n",daux_diff_2dr_steer_steer);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// d^2f[6]/du[i] du[j]
// j= 0
// j= 1
// d^2f[7]/du[i] du[j]
// j= 0
// j= 1
// d^2f[8]/du[i] du[j]
// j= 0
// j= 1
// d^2f[9]/du[i] du[j]
// j= 0
// j= 1

// d^2f[0]/dx[i] du[j]
// j= 0
// j= 1
// d^2f[1]/dx[i] du[j]
// j= 0
// j= 1
// d^2f[2]/dx[i] du[j]
// j= 0
// j= 1
// d^2f[3]/dx[i] du[j]
// j= 0
    t->fxu[63]=p[15][0]*daux_diff_2dUx_Ux_thr;
    if (isNANorINF(t->fxu[63]))
        {
            PRNT("    @k %d: t->fxu[63] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxu[63]);
            {
                PRNT("        daux_diff_2dUx_Ux_thr= %g\n",daux_diff_2dUx_Ux_thr);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxu[64]=p[15][0]*daux_diff_2dUx_Uy_thr;
    if (isNANorINF(t->fxu[64]))
        {
            PRNT("    @k %d: t->fxu[64] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxu[64]);
            {
                PRNT("        daux_diff_2dUx_Uy_thr= %g\n",daux_diff_2dUx_Uy_thr);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxu[65]=p[15][0]*daux_diff_2dUx_r_thr;
    if (isNANorINF(t->fxu[65]))
        {
            PRNT("    @k %d: t->fxu[65] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxu[65]);
            {
                PRNT("        daux_diff_2dUx_r_thr= %g\n",daux_diff_2dUx_r_thr);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 1
    t->fxu[73]=p[15][0]*daux_diff_2dUx_Ux_steer;
    if (isNANorINF(t->fxu[73]))
        {
            PRNT("    @k %d: t->fxu[73] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxu[73]);
            {
                PRNT("        daux_diff_2dUx_Ux_steer= %g\n",daux_diff_2dUx_Ux_steer);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxu[74]=p[15][0]*daux_diff_2dUx_Uy_steer;
    if (isNANorINF(t->fxu[74]))
        {
            PRNT("    @k %d: t->fxu[74] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxu[74]);
            {
                PRNT("        daux_diff_2dUx_Uy_steer= %g\n",daux_diff_2dUx_Uy_steer);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxu[75]=p[15][0]*daux_diff_2dUx_r_steer;
    if (isNANorINF(t->fxu[75]))
        {
            PRNT("    @k %d: t->fxu[75] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxu[75]);
            {
                PRNT("        daux_diff_2dUx_r_steer= %g\n",daux_diff_2dUx_r_steer);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// d^2f[4]/dx[i] du[j]
// j= 0
    t->fxu[83]=p[15][0]*daux_diff_2dUy_Ux_thr;
    if (isNANorINF(t->fxu[83]))
        {
            PRNT("    @k %d: t->fxu[83] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxu[83]);
            {
                PRNT("        daux_diff_2dUy_Ux_thr= %g\n",daux_diff_2dUy_Ux_thr);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxu[84]=p[15][0]*daux_diff_2dUy_Uy_thr;
    if (isNANorINF(t->fxu[84]))
        {
            PRNT("    @k %d: t->fxu[84] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxu[84]);
            {
                PRNT("        daux_diff_2dUy_Uy_thr= %g\n",daux_diff_2dUy_Uy_thr);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxu[85]=p[15][0]*daux_diff_2dUy_r_thr;
    if (isNANorINF(t->fxu[85]))
        {
            PRNT("    @k %d: t->fxu[85] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxu[85]);
            {
                PRNT("        daux_diff_2dUy_r_thr= %g\n",daux_diff_2dUy_r_thr);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 1
    t->fxu[93]=p[15][0]*daux_diff_2dUy_Ux_steer;
    if (isNANorINF(t->fxu[93]))
        {
            PRNT("    @k %d: t->fxu[93] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxu[93]);
            {
                PRNT("        daux_diff_2dUy_Ux_steer= %g\n",daux_diff_2dUy_Ux_steer);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxu[94]=p[15][0]*daux_diff_2dUy_Uy_steer;
    if (isNANorINF(t->fxu[94]))
        {
            PRNT("    @k %d: t->fxu[94] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxu[94]);
            {
                PRNT("        daux_diff_2dUy_Uy_steer= %g\n",daux_diff_2dUy_Uy_steer);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxu[95]=p[15][0]*daux_diff_2dUy_r_steer;
    if (isNANorINF(t->fxu[95]))
        {
            PRNT("    @k %d: t->fxu[95] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxu[95]);
            {
                PRNT("        daux_diff_2dUy_r_steer= %g\n",daux_diff_2dUy_r_steer);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// d^2f[5]/dx[i] du[j]
// j= 0
    t->fxu[103]=p[15][0]*daux_diff_2dr_Ux_thr;
    if (isNANorINF(t->fxu[103]))
        {
            PRNT("    @k %d: t->fxu[103] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxu[103]);
            {
                PRNT("        daux_diff_2dr_Ux_thr= %g\n",daux_diff_2dr_Ux_thr);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxu[104]=p[15][0]*daux_diff_2dr_Uy_thr;
    if (isNANorINF(t->fxu[104]))
        {
            PRNT("    @k %d: t->fxu[104] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxu[104]);
            {
                PRNT("        daux_diff_2dr_Uy_thr= %g\n",daux_diff_2dr_Uy_thr);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxu[105]=p[15][0]*daux_diff_2dr_r_thr;
    if (isNANorINF(t->fxu[105]))
        {
            PRNT("    @k %d: t->fxu[105] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxu[105]);
            {
                PRNT("        daux_diff_2dr_r_thr= %g\n",daux_diff_2dr_r_thr);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// j= 1
    t->fxu[113]=p[15][0]*daux_diff_2dr_Ux_steer;
    if (isNANorINF(t->fxu[113]))
        {
            PRNT("    @k %d: t->fxu[113] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxu[113]);
            {
                PRNT("        daux_diff_2dr_Ux_steer= %g\n",daux_diff_2dr_Ux_steer);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxu[114]=p[15][0]*daux_diff_2dr_Uy_steer;
    if (isNANorINF(t->fxu[114]))
        {
            PRNT("    @k %d: t->fxu[114] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxu[114]);
            {
                PRNT("        daux_diff_2dr_Uy_steer= %g\n",daux_diff_2dr_Uy_steer);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
    t->fxu[115]=p[15][0]*daux_diff_2dr_r_steer;
    if (isNANorINF(t->fxu[115]))
        {
            PRNT("    @k %d: t->fxu[115] in line %d is nan or inf: %g\n", k, __LINE__-3,t->fxu[115]);
            {
                PRNT("        daux_diff_2dr_r_steer= %g\n",daux_diff_2dr_r_steer);
                PRNT("        p[15,0]= %g\n",p[15][0]);
            }
            return(0);
        }
// d^2f[6]/dx[i] du[j]
// j= 0
// j= 1
// d^2f[7]/dx[i] du[j]
// j= 0
// j= 1
// d^2f[8]/dx[i] du[j]
// j= 0
// j= 1
// d^2f[9]/dx[i] du[j]
// j= 0
// j= 1
#endif

// derivatives of L
    t->cx[0]=p[13][0]*(x[0]-p[25][0])/sqrt(pow(p[24][0],2.0)+pow(x[0]-p[25][0],2.0))+p[16][0]*mcond(sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))>p[14][0],0.0,1,-2.0*(x[0]-p[3][0])*(1/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))-1/p[14][0])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5))+p[17][0]*mcond(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4]<0.0,0.0,1,x[3]/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0))-(x[0]-p[3][0])*(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5)/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0)));
    if (isNANorINF(t->cx[0]))
        {
            PRNT("    @k %d: t->cx[0] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cx[0]);
            {
                PRNT("        p[13,0]= %g\n",p[13][0]);
                PRNT("        p[25,0]= %g\n",p[25][0]);
                PRNT("        x[0]= %g\n",x[0]);
                PRNT("        p[24,0]= %g\n",p[24][0]);
                PRNT("        p[17,0]= %g\n",p[17][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[3,0]= %g\n",p[3][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[3,1]= %g\n",p[3][1]);
                PRNT("        x[1]= %g\n",x[1]);
                PRNT("        p[16,0]= %g\n",p[16][0]);
                PRNT("        p[14,0]= %g\n",p[14][0]);
            }
            return(0);
        }
    t->cx[1]=p[13][1]*(x[1]-p[25][1])/sqrt(pow(p[24][1],2.0)+pow(x[1]-p[25][1],2.0))+p[16][0]*mcond(sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))>p[14][0],0.0,1,-2.0*(x[1]-p[3][1])*(1/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))-1/p[14][0])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5))+p[17][0]*mcond(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4]<0.0,0.0,1,x[4]/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0))-(x[1]-p[3][1])*(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5)/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0)));
    if (isNANorINF(t->cx[1]))
        {
            PRNT("    @k %d: t->cx[1] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cx[1]);
            {
                PRNT("        p[13,1]= %g\n",p[13][1]);
                PRNT("        p[25,1]= %g\n",p[25][1]);
                PRNT("        x[1]= %g\n",x[1]);
                PRNT("        p[24,1]= %g\n",p[24][1]);
                PRNT("        p[17,0]= %g\n",p[17][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[3,0]= %g\n",p[3][0]);
                PRNT("        x[0]= %g\n",x[0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[3,1]= %g\n",p[3][1]);
                PRNT("        p[16,0]= %g\n",p[16][0]);
                PRNT("        p[14,0]= %g\n",p[14][0]);
            }
            return(0);
        }
    t->cx[2]=p[13][2]*(x[2]-p[25][2])/sqrt(pow(p[24][2],2.0)+pow(x[2]-p[25][2],2.0));
    if (isNANorINF(t->cx[2]))
        {
            PRNT("    @k %d: t->cx[2] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cx[2]);
            {
                PRNT("        p[13,2]= %g\n",p[13][2]);
                PRNT("        p[25,2]= %g\n",p[25][2]);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        p[24,2]= %g\n",p[24][2]);
            }
            return(0);
        }
    t->cx[3]=x[3]*p[8][0]/sqrt(1.0+pow(x[3],2.0))+p[10][0]*(x[3]-p[25][3])/sqrt(pow(p[24][0],2.0)+pow(x[3]-p[25][3],2.0))+p[17][0]*mcond(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4]<0.0,0.0,1,-x[3]*(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4])/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),1.5)+(x[0]-p[3][0])/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0)));
    if (isNANorINF(t->cx[3]))
        {
            PRNT("    @k %d: t->cx[3] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cx[3]);
            {
                PRNT("        p[10,0]= %g\n",p[10][0]);
                PRNT("        p[25,3]= %g\n",p[25][3]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[24,0]= %g\n",p[24][0]);
                PRNT("        p[8,0]= %g\n",p[8][0]);
                PRNT("        p[17,0]= %g\n",p[17][0]);
                PRNT("        p[3,0]= %g\n",p[3][0]);
                PRNT("        x[0]= %g\n",x[0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[3,1]= %g\n",p[3][1]);
                PRNT("        x[1]= %g\n",x[1]);
            }
            return(0);
        }
    t->cx[4]=p[10][1]*(x[4]-p[25][4])/sqrt(pow(p[24][1],2.0)+pow(x[4]-p[25][4],2.0))+p[17][0]*mcond(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4]<0.0,0.0,1,-x[4]*(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4])/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),1.5)+(x[1]-p[3][1])/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0)));
    if (isNANorINF(t->cx[4]))
        {
            PRNT("    @k %d: t->cx[4] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cx[4]);
            {
                PRNT("        p[10,1]= %g\n",p[10][1]);
                PRNT("        p[25,4]= %g\n",p[25][4]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[24,1]= %g\n",p[24][1]);
                PRNT("        p[17,0]= %g\n",p[17][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[3,0]= %g\n",p[3][0]);
                PRNT("        x[0]= %g\n",x[0]);
                PRNT("        p[3,1]= %g\n",p[3][1]);
                PRNT("        x[1]= %g\n",x[1]);
            }
            return(0);
        }
    t->cx[5]=p[10][2]*(x[5]-p[25][5])/sqrt(pow(p[24][2],2.0)+pow(x[5]-p[25][5],2.0));
    if (isNANorINF(t->cx[5]))
        {
            PRNT("    @k %d: t->cx[5] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cx[5]);
            {
                PRNT("        p[10,2]= %g\n",p[10][2]);
                PRNT("        p[25,5]= %g\n",p[25][5]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[24,2]= %g\n",p[24][2]);
            }
            return(0);
        }
    t->cx[6]=0.0;
    t->cx[7]=0.0;
    t->cx[8]=2.0*x[8]*p[9][0];
    if (isNANorINF(t->cx[8]))
        {
            PRNT("    @k %d: t->cx[8] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cx[8]);
            {
                PRNT("        p[9,0]= %g\n",p[9][0]);
                PRNT("        x[8]= %g\n",x[8]);
            }
            return(0);
        }
    t->cx[9]=2.0*x[9]*p[9][1];
    if (isNANorINF(t->cx[9]))
        {
            PRNT("    @k %d: t->cx[9] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cx[9]);
            {
                PRNT("        p[9,1]= %g\n",p[9][1]);
                PRNT("        x[9]= %g\n",x[9]);
            }
            return(0);
        }

    t->cxx[0]=p[13][0]/sqrt(pow(p[24][0],2.0)+pow(x[0]-p[25][0],2.0))-p[13][0]*pow(x[0]-p[25][0],2.0)/pow(pow(p[24][0],2.0)+pow(x[0]-p[25][0],2.0),1.5)+p[16][0]*mcond(sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))>p[14][0],0.0,1,2.0*pow(x[0]-p[3][0],2.0)/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),3.0)+6.0*pow(x[0]-p[3][0],2.0)*(1/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))-1/p[14][0])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),2.5)-2.0*(1/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))-1/p[14][0])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5))+p[17][0]*mcond(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4]<0.0,0.0,1,-2.0*x[3]*(x[0]-p[3][0])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5)/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0))+3.0*pow(x[0]-p[3][0],2.0)*(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),2.5)/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0))-(x[3]*(x[0]-p[3][0])+(x[1]-p[3
     ][1])*x[4])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5)/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0)));
    if (isNANorINF(t->cxx[0]))
        {
            PRNT("    @k %d: t->cxx[0] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[0]);
            {
                PRNT("        p[13,0]= %g\n",p[13][0]);
                PRNT("        p[25,0]= %g\n",p[25][0]);
                PRNT("        x[0]= %g\n",x[0]);
                PRNT("        p[24,0]= %g\n",p[24][0]);
                PRNT("        p[17,0]= %g\n",p[17][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[3,0]= %g\n",p[3][0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[3,1]= %g\n",p[3][1]);
                PRNT("        x[1]= %g\n",x[1]);
                PRNT("        p[16,0]= %g\n",p[16][0]);
                PRNT("        p[14,0]= %g\n",p[14][0]);
            }
            return(0);
        }
    t->cxx[1]=p[16][0]*mcond(sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))>p[14][0],0.0,1,2.0*(x[0]-p[3][0])*(x[1]-p[3][1])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),3.0)+6.0*(x[0]-p[3][0])*(x[1]-p[3][1])*(1/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))-1/p[14][0])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),2.5))+p[17][0]*mcond(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4]<0.0,0.0,1,-x[3]*(x[1]-p[3][1])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5)/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0))-(x[0]-p[3][0])*x[4]/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5)/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0))+3.0*(x[0]-p[3][0])*(x[1]-p[3][1])*(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),2.5)/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0)));
    if (isNANorINF(t->cxx[1]))
        {
            PRNT("    @k %d: t->cxx[1] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[1]);
            {
                PRNT("        p[17,0]= %g\n",p[17][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[3,0]= %g\n",p[3][0]);
                PRNT("        x[0]= %g\n",x[0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[3,1]= %g\n",p[3][1]);
                PRNT("        x[1]= %g\n",x[1]);
                PRNT("        p[16,0]= %g\n",p[16][0]);
                PRNT("        p[14,0]= %g\n",p[14][0]);
            }
            return(0);
        }
    t->cxx[2]=p[13][1]/sqrt(pow(p[24][1],2.0)+pow(x[1]-p[25][1],2.0))-p[13][1]*pow(x[1]-p[25][1],2.0)/pow(pow(p[24][1],2.0)+pow(x[1]-p[25][1],2.0),1.5)+p[16][0]*mcond(sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))>p[14][0],0.0,1,2.0*pow(x[1]-p[3][1],2.0)/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),3.0)-2.0*(1/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))-1/p[14][0])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5)+6.0*pow(x[1]-p[3][1],2.0)*(1/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))-1/p[14][0])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),2.5))+p[17][0]*mcond(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4]<0.0,0.0,1,-2.0*(x[1]-p[3][1])*x[4]/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5)/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0))-(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5)/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0))+3.0*pow(x[1]-p[3][1],2.0)*(x[3]*(x[0]-p[3][0])+(x[1]-p[3
     ][1])*x[4])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),2.5)/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0)));
    if (isNANorINF(t->cxx[2]))
        {
            PRNT("    @k %d: t->cxx[2] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[2]);
            {
                PRNT("        p[13,1]= %g\n",p[13][1]);
                PRNT("        p[25,1]= %g\n",p[25][1]);
                PRNT("        x[1]= %g\n",x[1]);
                PRNT("        p[24,1]= %g\n",p[24][1]);
                PRNT("        p[17,0]= %g\n",p[17][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[3,0]= %g\n",p[3][0]);
                PRNT("        x[0]= %g\n",x[0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[3,1]= %g\n",p[3][1]);
                PRNT("        p[16,0]= %g\n",p[16][0]);
                PRNT("        p[14,0]= %g\n",p[14][0]);
            }
            return(0);
        }
    t->cxx[3]=0.0;
    t->cxx[4]=0.0;
    t->cxx[5]=p[13][2]/sqrt(pow(p[24][2],2.0)+pow(x[2]-p[25][2],2.0))-p[13][2]*pow(x[2]-p[25][2],2.0)/pow(pow(p[24][2],2.0)+pow(x[2]-p[25][2],2.0),1.5);
    if (isNANorINF(t->cxx[5]))
        {
            PRNT("    @k %d: t->cxx[5] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[5]);
            {
                PRNT("        p[13,2]= %g\n",p[13][2]);
                PRNT("        p[25,2]= %g\n",p[25][2]);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        p[24,2]= %g\n",p[24][2]);
            }
            return(0);
        }
    t->cxx[6]=p[17][0]*mcond(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4]<0.0,0.0,1,-pow(x[3],2.0)/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),1.5)+x[3]*(x[0]-p[3][0])*(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5)/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),1.5)-pow(x[0]-p[3][0],2.0)/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5)/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0))+1/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0)));
    if (isNANorINF(t->cxx[6]))
        {
            PRNT("    @k %d: t->cxx[6] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[6]);
            {
                PRNT("        p[17,0]= %g\n",p[17][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[3,0]= %g\n",p[3][0]);
                PRNT("        x[0]= %g\n",x[0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[3,1]= %g\n",p[3][1]);
                PRNT("        x[1]= %g\n",x[1]);
            }
            return(0);
        }
    t->cxx[7]=p[17][0]*mcond(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4]<0.0,0.0,1,-x[3]*x[4]/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),1.5)+x[3]*(x[1]-p[3][1])*(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5)/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),1.5)-(x[0]-p[3][0])*(x[1]-p[3][1])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5)/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0)));
    if (isNANorINF(t->cxx[7]))
        {
            PRNT("    @k %d: t->cxx[7] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[7]);
            {
                PRNT("        p[17,0]= %g\n",p[17][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[3,0]= %g\n",p[3][0]);
                PRNT("        x[0]= %g\n",x[0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[3,1]= %g\n",p[3][1]);
                PRNT("        x[1]= %g\n",x[1]);
            }
            return(0);
        }
    t->cxx[8]=0.0;
    t->cxx[9]=-pow(x[3],2.0)*p[8][0]/pow(1.0+pow(x[3],2.0),1.5)+p[8][0]/sqrt(1.0+pow(x[3],2.0))+p[10][0]/sqrt(pow(p[24][0],2.0)+pow(x[3]-p[25][3],2.0))-p[10][0]*pow(x[3]-p[25][3],2.0)/pow(pow(p[24][0],2.0)+pow(x[3]-p[25][3],2.0),1.5)+p[17][0]*mcond(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4]<0.0,0.0,1,3.0*pow(x[3],2.0)*(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4])/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),2.5)-2.0*x[3]*(x[0]-p[3][0])/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),1.5)-(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4])/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),1.5));
    if (isNANorINF(t->cxx[9]))
        {
            PRNT("    @k %d: t->cxx[9] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[9]);
            {
                PRNT("        p[10,0]= %g\n",p[10][0]);
                PRNT("        p[25,3]= %g\n",p[25][3]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[24,0]= %g\n",p[24][0]);
                PRNT("        p[8,0]= %g\n",p[8][0]);
                PRNT("        p[17,0]= %g\n",p[17][0]);
                PRNT("        p[3,0]= %g\n",p[3][0]);
                PRNT("        x[0]= %g\n",x[0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[3,1]= %g\n",p[3][1]);
                PRNT("        x[1]= %g\n",x[1]);
            }
            return(0);
        }
    t->cxx[10]=p[17][0]*mcond(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4]<0.0,0.0,1,-x[3]*x[4]/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),1.5)+(x[0]-p[3][0])*x[4]*(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5)/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),1.5)-(x[0]-p[3][0])*(x[1]-p[3][1])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5)/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0)));
    if (isNANorINF(t->cxx[10]))
        {
            PRNT("    @k %d: t->cxx[10] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[10]);
            {
                PRNT("        p[17,0]= %g\n",p[17][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[3,0]= %g\n",p[3][0]);
                PRNT("        x[0]= %g\n",x[0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[3,1]= %g\n",p[3][1]);
                PRNT("        x[1]= %g\n",x[1]);
            }
            return(0);
        }
    t->cxx[11]=p[17][0]*mcond(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4]<0.0,0.0,1,(x[1]-p[3][1])*x[4]*(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4])/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5)/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),1.5)-pow(x[4],2.0)/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),1.5)+1/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0))-pow(x[1]-p[3][1],2.0)/pow(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0),1.5)/sqrt(1.e-6+pow(x[3],2.0)+pow(x[4],2.0)));
    if (isNANorINF(t->cxx[11]))
        {
            PRNT("    @k %d: t->cxx[11] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[11]);
            {
                PRNT("        p[17,0]= %g\n",p[17][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[3,0]= %g\n",p[3][0]);
                PRNT("        x[0]= %g\n",x[0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[3,1]= %g\n",p[3][1]);
                PRNT("        x[1]= %g\n",x[1]);
            }
            return(0);
        }
    t->cxx[12]=0.0;
    t->cxx[13]=p[17][0]*mcond(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4]<0.0,0.0,1,3.0*x[3]*x[4]*(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4])/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),2.5)-x[3]*(x[1]-p[3][1])/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),1.5)-(x[0]-p[3][0])*x[4]/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),1.5));
    if (isNANorINF(t->cxx[13]))
        {
            PRNT("    @k %d: t->cxx[13] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[13]);
            {
                PRNT("        p[17,0]= %g\n",p[17][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[3,0]= %g\n",p[3][0]);
                PRNT("        x[0]= %g\n",x[0]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[3,1]= %g\n",p[3][1]);
                PRNT("        x[1]= %g\n",x[1]);
            }
            return(0);
        }
    t->cxx[14]=p[10][1]/sqrt(pow(p[24][1],2.0)+pow(x[4]-p[25][4],2.0))-p[10][1]*pow(x[4]-p[25][4],2.0)/pow(pow(p[24][1],2.0)+pow(x[4]-p[25][4],2.0),1.5)+p[17][0]*mcond(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4]<0.0,0.0,1,3.0*pow(x[4],2.0)*(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4])/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),2.5)-2.0*(x[1]-p[3][1])*x[4]/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),1.5)-(x[3]*(x[0]-p[3][0])+(x[1]-p[3][1])*x[4])/sqrt(1.e-6+pow(x[0]-p[3][0],2.0)+pow(x[1]-p[3][1],2.0))/pow(1.e-6+pow(x[3],2.0)+pow(x[4],2.0),1.5));
    if (isNANorINF(t->cxx[14]))
        {
            PRNT("    @k %d: t->cxx[14] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[14]);
            {
                PRNT("        p[10,1]= %g\n",p[10][1]);
                PRNT("        p[25,4]= %g\n",p[25][4]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[24,1]= %g\n",p[24][1]);
                PRNT("        p[17,0]= %g\n",p[17][0]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[3,0]= %g\n",p[3][0]);
                PRNT("        x[0]= %g\n",x[0]);
                PRNT("        p[3,1]= %g\n",p[3][1]);
                PRNT("        x[1]= %g\n",x[1]);
            }
            return(0);
        }
    t->cxx[15]=0.0;
    t->cxx[16]=0.0;
    t->cxx[17]=0.0;
    t->cxx[18]=0.0;
    t->cxx[19]=0.0;
    t->cxx[20]=p[10][2]/sqrt(pow(p[24][2],2.0)+pow(x[5]-p[25][5],2.0))-p[10][2]*pow(x[5]-p[25][5],2.0)/pow(pow(p[24][2],2.0)+pow(x[5]-p[25][5],2.0),1.5);
    if (isNANorINF(t->cxx[20]))
        {
            PRNT("    @k %d: t->cxx[20] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[20]);
            {
                PRNT("        p[10,2]= %g\n",p[10][2]);
                PRNT("        p[25,5]= %g\n",p[25][5]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[24,2]= %g\n",p[24][2]);
            }
            return(0);
        }
    t->cxx[21]=0.0;
    t->cxx[22]=0.0;
    t->cxx[23]=0.0;
    t->cxx[24]=0.0;
    t->cxx[25]=0.0;
    t->cxx[26]=0.0;
    t->cxx[27]=0.0;
    t->cxx[28]=0.0;
    t->cxx[29]=0.0;
    t->cxx[30]=0.0;
    t->cxx[31]=0.0;
    t->cxx[32]=0.0;
    t->cxx[33]=0.0;
    t->cxx[34]=0.0;
    t->cxx[35]=0.0;
    t->cxx[36]=0.0;
    t->cxx[37]=0.0;
    t->cxx[38]=0.0;
    t->cxx[39]=0.0;
    t->cxx[40]=0.0;
    t->cxx[41]=0.0;
    t->cxx[42]=0.0;
    t->cxx[43]=0.0;
    t->cxx[44]=2.0*p[9][0];
    if (isNANorINF(t->cxx[44]))
        {
            PRNT("    @k %d: t->cxx[44] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[44]);
            {
                PRNT("        p[9,0]= %g\n",p[9][0]);
            }
            return(0);
        }
    t->cxx[45]=0.0;
    t->cxx[46]=0.0;
    t->cxx[47]=0.0;
    t->cxx[48]=0.0;
    t->cxx[49]=0.0;
    t->cxx[50]=0.0;
    t->cxx[51]=0.0;
    t->cxx[52]=0.0;
    t->cxx[53]=0.0;
    t->cxx[54]=2.0*p[9][1];
    if (isNANorINF(t->cxx[54]))
        {
            PRNT("    @k %d: t->cxx[54] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[54]);
            {
                PRNT("        p[9,1]= %g\n",p[9][1]);
            }
            return(0);
        }

    t->cu[0]=2.0*u[0]*p[12][0];
    if (isNANorINF(t->cu[0]))
        {
            PRNT("    @k %d: t->cu[0] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cu[0]);
            {
                PRNT("        p[12,0]= %g\n",p[12][0]);
                PRNT("        u[0]= %g\n",u[0]);
            }
            return(0);
        }
    t->cu[1]=2.0*u[1]*p[12][1];
    if (isNANorINF(t->cu[1]))
        {
            PRNT("    @k %d: t->cu[1] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cu[1]);
            {
                PRNT("        p[12,1]= %g\n",p[12][1]);
                PRNT("        u[1]= %g\n",u[1]);
            }
            return(0);
        }

    t->cuu[0]=2.0*p[12][0];
    if (isNANorINF(t->cuu[0]))
        {
            PRNT("    @k %d: t->cuu[0] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cuu[0]);
            {
                PRNT("        p[12,0]= %g\n",p[12][0]);
            }
            return(0);
        }
    t->cuu[1]=0.0;
    t->cuu[2]=2.0*p[12][1];
    if (isNANorINF(t->cuu[2]))
        {
            PRNT("    @k %d: t->cuu[2] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cuu[2]);
            {
                PRNT("        p[12,1]= %g\n",p[12][1]);
            }
            return(0);
        }

    t->cxu[0]=0.0;
    t->cxu[1]=0.0;
    t->cxu[2]=0.0;
    t->cxu[3]=0.0;
    t->cxu[4]=0.0;
    t->cxu[5]=0.0;
    t->cxu[6]=0.0;
    t->cxu[7]=0.0;
    t->cxu[8]=0.0;
    t->cxu[9]=0.0;
    t->cxu[10]=0.0;
    t->cxu[11]=0.0;
    t->cxu[12]=0.0;
    t->cxu[13]=0.0;
    t->cxu[14]=0.0;
    t->cxu[15]=0.0;
    t->cxu[16]=0.0;
    t->cxu[17]=0.0;
    t->cxu[18]=0.0;
    t->cxu[19]=0.0;

    return 1;
}

static int calcFAuxDeriv(trajFin_t *t, multipliersFin_t *m, tOptSet *o) {
    const double *x= t->x;
    const double w_pen= o->w_pen_f;
    double **p= o->p;
    const int k= o->n_hor;
    
    return 1;
}

static int bp_derivsF(trajFin_t *t, int k, double **p) {
    const double *x= t->x;
    
    t->cx[0]=p[11][0]*(x[0]-p[25][0])/sqrt(pow(p[23][0],2.0)+pow(x[0]-p[25][0],2.0));
    if (isNANorINF(t->cx[0]))
        {
            PRNT("    @k %d: t->cx[0] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cx[0]);
            {
                PRNT("        p[11,0]= %g\n",p[11][0]);
                PRNT("        p[25,0]= %g\n",p[25][0]);
                PRNT("        x[0]= %g\n",x[0]);
                PRNT("        p[23,0]= %g\n",p[23][0]);
            }
            return(0);
        }
    t->cx[1]=p[11][1]*(x[1]-p[25][1])/sqrt(pow(p[23][1],2.0)+pow(x[1]-p[25][1],2.0));
    if (isNANorINF(t->cx[1]))
        {
            PRNT("    @k %d: t->cx[1] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cx[1]);
            {
                PRNT("        p[11,1]= %g\n",p[11][1]);
                PRNT("        p[25,1]= %g\n",p[25][1]);
                PRNT("        x[1]= %g\n",x[1]);
                PRNT("        p[23,1]= %g\n",p[23][1]);
            }
            return(0);
        }
    t->cx[2]=p[11][2]*(x[2]-p[25][2])/sqrt(pow(p[23][2],2.0)+pow(x[2]-p[25][2],2.0));
    if (isNANorINF(t->cx[2]))
        {
            PRNT("    @k %d: t->cx[2] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cx[2]);
            {
                PRNT("        p[11,2]= %g\n",p[11][2]);
                PRNT("        p[25,2]= %g\n",p[25][2]);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        p[23,2]= %g\n",p[23][2]);
            }
            return(0);
        }
    t->cx[3]=p[11][3]*(x[3]-p[25][3])/sqrt(pow(p[23][3],2.0)+pow(x[3]-p[25][3],2.0));
    if (isNANorINF(t->cx[3]))
        {
            PRNT("    @k %d: t->cx[3] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cx[3]);
            {
                PRNT("        p[11,3]= %g\n",p[11][3]);
                PRNT("        p[25,3]= %g\n",p[25][3]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[23,3]= %g\n",p[23][3]);
            }
            return(0);
        }
    t->cx[4]=p[11][4]*(x[4]-p[25][4])/sqrt(pow(p[23][4],2.0)+pow(x[4]-p[25][4],2.0));
    if (isNANorINF(t->cx[4]))
        {
            PRNT("    @k %d: t->cx[4] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cx[4]);
            {
                PRNT("        p[11,4]= %g\n",p[11][4]);
                PRNT("        p[25,4]= %g\n",p[25][4]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[23,4]= %g\n",p[23][4]);
            }
            return(0);
        }
    t->cx[5]=p[11][5]*(x[5]-p[25][5])/sqrt(pow(p[23][5],2.0)+pow(x[5]-p[25][5],2.0));
    if (isNANorINF(t->cx[5]))
        {
            PRNT("    @k %d: t->cx[5] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cx[5]);
            {
                PRNT("        p[11,5]= %g\n",p[11][5]);
                PRNT("        p[25,5]= %g\n",p[25][5]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[23,5]= %g\n",p[23][5]);
            }
            return(0);
        }
    t->cx[6]=0.0;
    t->cx[7]=0.0;
    t->cx[8]=0.0;
    t->cx[9]=0.0;

    t->cxx[0]=p[11][0]/sqrt(pow(p[23][0],2.0)+pow(x[0]-p[25][0],2.0))-p[11][0]*pow(x[0]-p[25][0],2.0)/pow(pow(p[23][0],2.0)+pow(x[0]-p[25][0],2.0),1.5);
    if (isNANorINF(t->cxx[0]))
        {
            PRNT("    @k %d: t->cxx[0] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[0]);
            {
                PRNT("        p[11,0]= %g\n",p[11][0]);
                PRNT("        p[25,0]= %g\n",p[25][0]);
                PRNT("        x[0]= %g\n",x[0]);
                PRNT("        p[23,0]= %g\n",p[23][0]);
            }
            return(0);
        }
    t->cxx[1]=0.0;
    t->cxx[2]=p[11][1]/sqrt(pow(p[23][1],2.0)+pow(x[1]-p[25][1],2.0))-p[11][1]*pow(x[1]-p[25][1],2.0)/pow(pow(p[23][1],2.0)+pow(x[1]-p[25][1],2.0),1.5);
    if (isNANorINF(t->cxx[2]))
        {
            PRNT("    @k %d: t->cxx[2] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[2]);
            {
                PRNT("        p[11,1]= %g\n",p[11][1]);
                PRNT("        p[25,1]= %g\n",p[25][1]);
                PRNT("        x[1]= %g\n",x[1]);
                PRNT("        p[23,1]= %g\n",p[23][1]);
            }
            return(0);
        }
    t->cxx[3]=0.0;
    t->cxx[4]=0.0;
    t->cxx[5]=p[11][2]/sqrt(pow(p[23][2],2.0)+pow(x[2]-p[25][2],2.0))-p[11][2]*pow(x[2]-p[25][2],2.0)/pow(pow(p[23][2],2.0)+pow(x[2]-p[25][2],2.0),1.5);
    if (isNANorINF(t->cxx[5]))
        {
            PRNT("    @k %d: t->cxx[5] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[5]);
            {
                PRNT("        p[11,2]= %g\n",p[11][2]);
                PRNT("        p[25,2]= %g\n",p[25][2]);
                PRNT("        x[2]= %g\n",x[2]);
                PRNT("        p[23,2]= %g\n",p[23][2]);
            }
            return(0);
        }
    t->cxx[6]=0.0;
    t->cxx[7]=0.0;
    t->cxx[8]=0.0;
    t->cxx[9]=p[11][3]/sqrt(pow(p[23][3],2.0)+pow(x[3]-p[25][3],2.0))-p[11][3]*pow(x[3]-p[25][3],2.0)/pow(pow(p[23][3],2.0)+pow(x[3]-p[25][3],2.0),1.5);
    if (isNANorINF(t->cxx[9]))
        {
            PRNT("    @k %d: t->cxx[9] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[9]);
            {
                PRNT("        p[11,3]= %g\n",p[11][3]);
                PRNT("        p[25,3]= %g\n",p[25][3]);
                PRNT("        x[3]= %g\n",x[3]);
                PRNT("        p[23,3]= %g\n",p[23][3]);
            }
            return(0);
        }
    t->cxx[10]=0.0;
    t->cxx[11]=0.0;
    t->cxx[12]=0.0;
    t->cxx[13]=0.0;
    t->cxx[14]=p[11][4]/sqrt(pow(p[23][4],2.0)+pow(x[4]-p[25][4],2.0))-p[11][4]*pow(x[4]-p[25][4],2.0)/pow(pow(p[23][4],2.0)+pow(x[4]-p[25][4],2.0),1.5);
    if (isNANorINF(t->cxx[14]))
        {
            PRNT("    @k %d: t->cxx[14] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[14]);
            {
                PRNT("        p[11,4]= %g\n",p[11][4]);
                PRNT("        p[25,4]= %g\n",p[25][4]);
                PRNT("        x[4]= %g\n",x[4]);
                PRNT("        p[23,4]= %g\n",p[23][4]);
            }
            return(0);
        }
    t->cxx[15]=0.0;
    t->cxx[16]=0.0;
    t->cxx[17]=0.0;
    t->cxx[18]=0.0;
    t->cxx[19]=0.0;
    t->cxx[20]=p[11][5]/sqrt(pow(p[23][5],2.0)+pow(x[5]-p[25][5],2.0))-p[11][5]*pow(x[5]-p[25][5],2.0)/pow(pow(p[23][5],2.0)+pow(x[5]-p[25][5],2.0),1.5);
    if (isNANorINF(t->cxx[20]))
        {
            PRNT("    @k %d: t->cxx[20] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[20]);
            {
                PRNT("        p[11,5]= %g\n",p[11][5]);
                PRNT("        p[25,5]= %g\n",p[25][5]);
                PRNT("        x[5]= %g\n",x[5]);
                PRNT("        p[23,5]= %g\n",p[23][5]);
            }
            return(0);
        }
    t->cxx[21]=0.0;
    t->cxx[22]=0.0;
    t->cxx[23]=0.0;
    t->cxx[24]=0.0;
    t->cxx[25]=0.0;
    t->cxx[26]=0.0;
    t->cxx[27]=0.0;
    t->cxx[28]=0.0;
    t->cxx[29]=0.0;
    t->cxx[30]=0.0;
    t->cxx[31]=0.0;
    t->cxx[32]=0.0;
    t->cxx[33]=0.0;
    t->cxx[34]=0.0;
    t->cxx[35]=0.0;
    t->cxx[36]=0.0;
    t->cxx[37]=0.0;
    t->cxx[38]=0.0;
    t->cxx[39]=0.0;
    t->cxx[40]=0.0;
    t->cxx[41]=0.0;
    t->cxx[42]=0.0;
    t->cxx[43]=0.0;
    t->cxx[44]=0.0;
    t->cxx[45]=0.0;
    t->cxx[46]=0.0;
    t->cxx[47]=0.0;
    t->cxx[48]=0.0;
    t->cxx[49]=0.0;
    t->cxx[50]=0.0;
    t->cxx[51]=0.0;
    t->cxx[52]=0.0;
    t->cxx[53]=0.0;
    t->cxx[54]=0.0;
    return 1;
}

static int init_running(trajEl_t *t, tOptSet *o) {
    int k;
    double **p= o->p;
    
    for(k= 0; k<o->n_hor; k++, t++) {


// derivatives of L
    t->cx[6]=0.0;
    t->cx[7]=0.0;

    t->cxx[3]=0.0;
    t->cxx[4]=0.0;
    t->cxx[8]=0.0;
    t->cxx[12]=0.0;
    t->cxx[15]=0.0;
    t->cxx[16]=0.0;
    t->cxx[17]=0.0;
    t->cxx[18]=0.0;
    t->cxx[19]=0.0;
    t->cxx[21]=0.0;
    t->cxx[22]=0.0;
    t->cxx[23]=0.0;
    t->cxx[24]=0.0;
    t->cxx[25]=0.0;
    t->cxx[26]=0.0;
    t->cxx[27]=0.0;
    t->cxx[28]=0.0;
    t->cxx[29]=0.0;
    t->cxx[30]=0.0;
    t->cxx[31]=0.0;
    t->cxx[32]=0.0;
    t->cxx[33]=0.0;
    t->cxx[34]=0.0;
    t->cxx[35]=0.0;
    t->cxx[36]=0.0;
    t->cxx[37]=0.0;
    t->cxx[38]=0.0;
    t->cxx[39]=0.0;
    t->cxx[40]=0.0;
    t->cxx[41]=0.0;
    t->cxx[42]=0.0;
    t->cxx[43]=0.0;
    t->cxx[44]=2.0*p[9][0];
    if (isNANorINF(t->cxx[44]))
        {
            PRNT("    @k %d: t->cxx[44] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[44]);
            {
                PRNT("        p[9,0]= %g\n",p[9][0]);
            }
            return(0);
        }
    t->cxx[45]=0.0;
    t->cxx[46]=0.0;
    t->cxx[47]=0.0;
    t->cxx[48]=0.0;
    t->cxx[49]=0.0;
    t->cxx[50]=0.0;
    t->cxx[51]=0.0;
    t->cxx[52]=0.0;
    t->cxx[53]=0.0;
    t->cxx[54]=2.0*p[9][1];
    if (isNANorINF(t->cxx[54]))
        {
            PRNT("    @k %d: t->cxx[54] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cxx[54]);
            {
                PRNT("        p[9,1]= %g\n",p[9][1]);
            }
            return(0);
        }


    t->cuu[0]=2.0*p[12][0];
    if (isNANorINF(t->cuu[0]))
        {
            PRNT("    @k %d: t->cuu[0] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cuu[0]);
            {
                PRNT("        p[12,0]= %g\n",p[12][0]);
            }
            return(0);
        }
    t->cuu[1]=0.0;
    t->cuu[2]=2.0*p[12][1];
    if (isNANorINF(t->cuu[2]))
        {
            PRNT("    @k %d: t->cuu[2] in line %d is nan or inf: %g\n", k, __LINE__-3,t->cuu[2]);
            {
                PRNT("        p[12,1]= %g\n",p[12][1]);
            }
            return(0);
        }

    t->cxu[0]=0.0;
    t->cxu[1]=0.0;
    t->cxu[2]=0.0;
    t->cxu[3]=0.0;
    t->cxu[4]=0.0;
    t->cxu[5]=0.0;
    t->cxu[6]=0.0;
    t->cxu[7]=0.0;
    t->cxu[8]=0.0;
    t->cxu[9]=0.0;
    t->cxu[10]=0.0;
    t->cxu[11]=0.0;
    t->cxu[12]=0.0;
    t->cxu[13]=0.0;
    t->cxu[14]=0.0;
    t->cxu[15]=0.0;
    t->cxu[16]=0.0;
    t->cxu[17]=0.0;
    t->cxu[18]=0.0;
    t->cxu[19]=0.0;

// derivatives of f
// df[i]/d x[0]
    t->fx[0]=1.0;
    t->fx[1]=0.0;
    t->fx[2]=0.0;
    t->fx[3]=0.0;
    t->fx[4]=0.0;
    t->fx[5]=0.0;
    t->fx[6]=0.0;
    t->fx[7]=0.0;
    t->fx[8]=0.0;
    t->fx[9]=0.0;
// df[i]/d x[1]
    t->fx[10]=0.0;
    t->fx[11]=1.0;
    t->fx[12]=0.0;
    t->fx[13]=0.0;
    t->fx[14]=0.0;
    t->fx[15]=0.0;
    t->fx[16]=0.0;
    t->fx[17]=0.0;
    t->fx[18]=0.0;
    t->fx[19]=0.0;
// df[i]/d x[2]
    t->fx[22]=1.0;
    t->fx[23]=0.0;
    t->fx[24]=0.0;
    t->fx[25]=0.0;
    t->fx[26]=0.0;
    t->fx[27]=0.0;
    t->fx[28]=0.0;
    t->fx[29]=0.0;
// df[i]/d x[3]
    t->fx[32]=0.0;
    t->fx[36]=0.0;
    t->fx[37]=0.0;
    t->fx[38]=0.0;
    t->fx[39]=0.0;
// df[i]/d x[4]
    t->fx[42]=0.0;
    t->fx[46]=0.0;
    t->fx[47]=0.0;
    t->fx[48]=0.0;
    t->fx[49]=0.0;
// df[i]/d x[5]
    t->fx[50]=0.0;
    t->fx[51]=0.0;
    t->fx[52]=p[15][0];
    t->fx[56]=0.0;
    t->fx[57]=0.0;
    t->fx[58]=0.0;
    t->fx[59]=0.0;
// df[i]/d x[6]
    t->fx[60]=0.0;
    t->fx[61]=0.0;
    t->fx[62]=0.0;
    t->fx[63]=0.0;
    t->fx[64]=0.0;
    t->fx[65]=0.0;
    t->fx[66]=0.0;
    t->fx[67]=0.0;
    t->fx[68]=-1.0;
    t->fx[69]=0.0;
// df[i]/d x[7]
    t->fx[70]=0.0;
    t->fx[71]=0.0;
    t->fx[72]=0.0;
    t->fx[73]=0.0;
    t->fx[74]=0.0;
    t->fx[75]=0.0;
    t->fx[76]=0.0;
    t->fx[77]=0.0;
    t->fx[78]=0.0;
    t->fx[79]=-1.0;
// df[i]/d x[8]
    t->fx[80]=0.0;
    t->fx[81]=0.0;
    t->fx[82]=0.0;
    t->fx[83]=0.0;
    t->fx[84]=0.0;
    t->fx[85]=0.0;
    t->fx[86]=0.0;
    t->fx[87]=0.0;
    t->fx[88]=0.0;
    t->fx[89]=0.0;
// df[i]/d x[9]
    t->fx[90]=0.0;
    t->fx[91]=0.0;
    t->fx[92]=0.0;
    t->fx[93]=0.0;
    t->fx[94]=0.0;
    t->fx[95]=0.0;
    t->fx[96]=0.0;
    t->fx[97]=0.0;
    t->fx[98]=0.0;
    t->fx[99]=0.0;

// df[i]/d u[0]
    t->fu[0]=0.0;
    t->fu[1]=0.0;
    t->fu[2]=0.0;
    t->fu[6]=1.0;
    t->fu[7]=0.0;
    t->fu[8]=1.0;
    t->fu[9]=0.0;
// df[i]/d u[1]
    t->fu[10]=0.0;
    t->fu[11]=0.0;
    t->fu[12]=0.0;
    t->fu[16]=0.0;
    t->fu[17]=1.0;
    t->fu[18]=0.0;
    t->fu[19]=1.0;

#if FULL_DDP
// d^2f[0]/dx[i] dx[j]
// j= 0
    t->fxx[0]=0.0;
// j= 1
    t->fxx[1]=0.0;
    t->fxx[2]=0.0;
// j= 2
    t->fxx[3]=0.0;
    t->fxx[4]=0.0;
// j= 3
    t->fxx[6]=0.0;
    t->fxx[7]=0.0;
    t->fxx[8]=0.0;
// j= 4
    t->fxx[10]=0.0;
    t->fxx[11]=0.0;
    t->fxx[12]=0.0;
// j= 5
    t->fxx[15]=0.0;
    t->fxx[16]=0.0;
    t->fxx[17]=0.0;
    t->fxx[18]=0.0;
    t->fxx[19]=0.0;
    t->fxx[20]=0.0;
// j= 6
    t->fxx[21]=0.0;
    t->fxx[22]=0.0;
    t->fxx[23]=0.0;
    t->fxx[24]=0.0;
    t->fxx[25]=0.0;
    t->fxx[26]=0.0;
    t->fxx[27]=0.0;
// j= 7
    t->fxx[28]=0.0;
    t->fxx[29]=0.0;
    t->fxx[30]=0.0;
    t->fxx[31]=0.0;
    t->fxx[32]=0.0;
    t->fxx[33]=0.0;
    t->fxx[34]=0.0;
    t->fxx[35]=0.0;
// j= 8
    t->fxx[36]=0.0;
    t->fxx[37]=0.0;
    t->fxx[38]=0.0;
    t->fxx[39]=0.0;
    t->fxx[40]=0.0;
    t->fxx[41]=0.0;
    t->fxx[42]=0.0;
    t->fxx[43]=0.0;
    t->fxx[44]=0.0;
// j= 9
    t->fxx[45]=0.0;
    t->fxx[46]=0.0;
    t->fxx[47]=0.0;
    t->fxx[48]=0.0;
    t->fxx[49]=0.0;
    t->fxx[50]=0.0;
    t->fxx[51]=0.0;
    t->fxx[52]=0.0;
    t->fxx[53]=0.0;
    t->fxx[54]=0.0;
// d^2f[1]/dx[i] dx[j]
// j= 0
    t->fxx[55]=0.0;
// j= 1
    t->fxx[56]=0.0;
    t->fxx[57]=0.0;
// j= 2
    t->fxx[58]=0.0;
    t->fxx[59]=0.0;
// j= 3
    t->fxx[61]=0.0;
    t->fxx[62]=0.0;
    t->fxx[63]=0.0;
// j= 4
    t->fxx[65]=0.0;
    t->fxx[66]=0.0;
    t->fxx[67]=0.0;
// j= 5
    t->fxx[70]=0.0;
    t->fxx[71]=0.0;
    t->fxx[72]=0.0;
    t->fxx[73]=0.0;
    t->fxx[74]=0.0;
    t->fxx[75]=0.0;
// j= 6
    t->fxx[76]=0.0;
    t->fxx[77]=0.0;
    t->fxx[78]=0.0;
    t->fxx[79]=0.0;
    t->fxx[80]=0.0;
    t->fxx[81]=0.0;
    t->fxx[82]=0.0;
// j= 7
    t->fxx[83]=0.0;
    t->fxx[84]=0.0;
    t->fxx[85]=0.0;
    t->fxx[86]=0.0;
    t->fxx[87]=0.0;
    t->fxx[88]=0.0;
    t->fxx[89]=0.0;
    t->fxx[90]=0.0;
// j= 8
    t->fxx[91]=0.0;
    t->fxx[92]=0.0;
    t->fxx[93]=0.0;
    t->fxx[94]=0.0;
    t->fxx[95]=0.0;
    t->fxx[96]=0.0;
    t->fxx[97]=0.0;
    t->fxx[98]=0.0;
    t->fxx[99]=0.0;
// j= 9
    t->fxx[100]=0.0;
    t->fxx[101]=0.0;
    t->fxx[102]=0.0;
    t->fxx[103]=0.0;
    t->fxx[104]=0.0;
    t->fxx[105]=0.0;
    t->fxx[106]=0.0;
    t->fxx[107]=0.0;
    t->fxx[108]=0.0;
    t->fxx[109]=0.0;
// d^2f[2]/dx[i] dx[j]
// j= 0
    t->fxx[110]=0.0;
// j= 1
    t->fxx[111]=0.0;
    t->fxx[112]=0.0;
// j= 2
    t->fxx[113]=0.0;
    t->fxx[114]=0.0;
    t->fxx[115]=0.0;
// j= 3
    t->fxx[116]=0.0;
    t->fxx[117]=0.0;
    t->fxx[118]=0.0;
    t->fxx[119]=0.0;
// j= 4
    t->fxx[120]=0.0;
    t->fxx[121]=0.0;
    t->fxx[122]=0.0;
    t->fxx[123]=0.0;
    t->fxx[124]=0.0;
// j= 5
    t->fxx[125]=0.0;
    t->fxx[126]=0.0;
    t->fxx[127]=0.0;
    t->fxx[128]=0.0;
    t->fxx[129]=0.0;
    t->fxx[130]=0.0;
// j= 6
    t->fxx[131]=0.0;
    t->fxx[132]=0.0;
    t->fxx[133]=0.0;
    t->fxx[134]=0.0;
    t->fxx[135]=0.0;
    t->fxx[136]=0.0;
    t->fxx[137]=0.0;
// j= 7
    t->fxx[138]=0.0;
    t->fxx[139]=0.0;
    t->fxx[140]=0.0;
    t->fxx[141]=0.0;
    t->fxx[142]=0.0;
    t->fxx[143]=0.0;
    t->fxx[144]=0.0;
    t->fxx[145]=0.0;
// j= 8
    t->fxx[146]=0.0;
    t->fxx[147]=0.0;
    t->fxx[148]=0.0;
    t->fxx[149]=0.0;
    t->fxx[150]=0.0;
    t->fxx[151]=0.0;
    t->fxx[152]=0.0;
    t->fxx[153]=0.0;
    t->fxx[154]=0.0;
// j= 9
    t->fxx[155]=0.0;
    t->fxx[156]=0.0;
    t->fxx[157]=0.0;
    t->fxx[158]=0.0;
    t->fxx[159]=0.0;
    t->fxx[160]=0.0;
    t->fxx[161]=0.0;
    t->fxx[162]=0.0;
    t->fxx[163]=0.0;
    t->fxx[164]=0.0;
// d^2f[3]/dx[i] dx[j]
// j= 0
    t->fxx[165]=0.0;
// j= 1
    t->fxx[166]=0.0;
    t->fxx[167]=0.0;
// j= 2
    t->fxx[168]=0.0;
    t->fxx[169]=0.0;
    t->fxx[170]=0.0;
// j= 3
    t->fxx[171]=0.0;
    t->fxx[172]=0.0;
    t->fxx[173]=0.0;
// j= 4
    t->fxx[175]=0.0;
    t->fxx[176]=0.0;
    t->fxx[177]=0.0;
// j= 5
    t->fxx[180]=0.0;
    t->fxx[181]=0.0;
    t->fxx[182]=0.0;
// j= 6
    t->fxx[186]=0.0;
    t->fxx[187]=0.0;
    t->fxx[188]=0.0;
    t->fxx[189]=0.0;
    t->fxx[190]=0.0;
    t->fxx[191]=0.0;
    t->fxx[192]=0.0;
// j= 7
    t->fxx[193]=0.0;
    t->fxx[194]=0.0;
    t->fxx[195]=0.0;
    t->fxx[196]=0.0;
    t->fxx[197]=0.0;
    t->fxx[198]=0.0;
    t->fxx[199]=0.0;
    t->fxx[200]=0.0;
// j= 8
    t->fxx[201]=0.0;
    t->fxx[202]=0.0;
    t->fxx[203]=0.0;
    t->fxx[204]=0.0;
    t->fxx[205]=0.0;
    t->fxx[206]=0.0;
    t->fxx[207]=0.0;
    t->fxx[208]=0.0;
    t->fxx[209]=0.0;
// j= 9
    t->fxx[210]=0.0;
    t->fxx[211]=0.0;
    t->fxx[212]=0.0;
    t->fxx[213]=0.0;
    t->fxx[214]=0.0;
    t->fxx[215]=0.0;
    t->fxx[216]=0.0;
    t->fxx[217]=0.0;
    t->fxx[218]=0.0;
    t->fxx[219]=0.0;
// d^2f[4]/dx[i] dx[j]
// j= 0
    t->fxx[220]=0.0;
// j= 1
    t->fxx[221]=0.0;
    t->fxx[222]=0.0;
// j= 2
    t->fxx[223]=0.0;
    t->fxx[224]=0.0;
    t->fxx[225]=0.0;
// j= 3
    t->fxx[226]=0.0;
    t->fxx[227]=0.0;
    t->fxx[228]=0.0;
// j= 4
    t->fxx[230]=0.0;
    t->fxx[231]=0.0;
    t->fxx[232]=0.0;
// j= 5
    t->fxx[235]=0.0;
    t->fxx[236]=0.0;
    t->fxx[237]=0.0;
// j= 6
    t->fxx[241]=0.0;
    t->fxx[242]=0.0;
    t->fxx[243]=0.0;
    t->fxx[244]=0.0;
    t->fxx[245]=0.0;
    t->fxx[246]=0.0;
    t->fxx[247]=0.0;
// j= 7
    t->fxx[248]=0.0;
    t->fxx[249]=0.0;
    t->fxx[250]=0.0;
    t->fxx[251]=0.0;
    t->fxx[252]=0.0;
    t->fxx[253]=0.0;
    t->fxx[254]=0.0;
    t->fxx[255]=0.0;
// j= 8
    t->fxx[256]=0.0;
    t->fxx[257]=0.0;
    t->fxx[258]=0.0;
    t->fxx[259]=0.0;
    t->fxx[260]=0.0;
    t->fxx[261]=0.0;
    t->fxx[262]=0.0;
    t->fxx[263]=0.0;
    t->fxx[264]=0.0;
// j= 9
    t->fxx[265]=0.0;
    t->fxx[266]=0.0;
    t->fxx[267]=0.0;
    t->fxx[268]=0.0;
    t->fxx[269]=0.0;
    t->fxx[270]=0.0;
    t->fxx[271]=0.0;
    t->fxx[272]=0.0;
    t->fxx[273]=0.0;
    t->fxx[274]=0.0;
// d^2f[5]/dx[i] dx[j]
// j= 0
    t->fxx[275]=0.0;
// j= 1
    t->fxx[276]=0.0;
    t->fxx[277]=0.0;
// j= 2
    t->fxx[278]=0.0;
    t->fxx[279]=0.0;
    t->fxx[280]=0.0;
// j= 3
    t->fxx[281]=0.0;
    t->fxx[282]=0.0;
    t->fxx[283]=0.0;
// j= 4
    t->fxx[285]=0.0;
    t->fxx[286]=0.0;
    t->fxx[287]=0.0;
// j= 5
    t->fxx[290]=0.0;
    t->fxx[291]=0.0;
    t->fxx[292]=0.0;
// j= 6
    t->fxx[296]=0.0;
    t->fxx[297]=0.0;
    t->fxx[298]=0.0;
    t->fxx[299]=0.0;
    t->fxx[300]=0.0;
    t->fxx[301]=0.0;
    t->fxx[302]=0.0;
// j= 7
    t->fxx[303]=0.0;
    t->fxx[304]=0.0;
    t->fxx[305]=0.0;
    t->fxx[306]=0.0;
    t->fxx[307]=0.0;
    t->fxx[308]=0.0;
    t->fxx[309]=0.0;
    t->fxx[310]=0.0;
// j= 8
    t->fxx[311]=0.0;
    t->fxx[312]=0.0;
    t->fxx[313]=0.0;
    t->fxx[314]=0.0;
    t->fxx[315]=0.0;
    t->fxx[316]=0.0;
    t->fxx[317]=0.0;
    t->fxx[318]=0.0;
    t->fxx[319]=0.0;
// j= 9
    t->fxx[320]=0.0;
    t->fxx[321]=0.0;
    t->fxx[322]=0.0;
    t->fxx[323]=0.0;
    t->fxx[324]=0.0;
    t->fxx[325]=0.0;
    t->fxx[326]=0.0;
    t->fxx[327]=0.0;
    t->fxx[328]=0.0;
    t->fxx[329]=0.0;
// d^2f[6]/dx[i] dx[j]
// j= 0
    t->fxx[330]=0.0;
// j= 1
    t->fxx[331]=0.0;
    t->fxx[332]=0.0;
// j= 2
    t->fxx[333]=0.0;
    t->fxx[334]=0.0;
    t->fxx[335]=0.0;
// j= 3
    t->fxx[336]=0.0;
    t->fxx[337]=0.0;
    t->fxx[338]=0.0;
    t->fxx[339]=0.0;
// j= 4
    t->fxx[340]=0.0;
    t->fxx[341]=0.0;
    t->fxx[342]=0.0;
    t->fxx[343]=0.0;
    t->fxx[344]=0.0;
// j= 5
    t->fxx[345]=0.0;
    t->fxx[346]=0.0;
    t->fxx[347]=0.0;
    t->fxx[348]=0.0;
    t->fxx[349]=0.0;
    t->fxx[350]=0.0;
// j= 6
    t->fxx[351]=0.0;
    t->fxx[352]=0.0;
    t->fxx[353]=0.0;
    t->fxx[354]=0.0;
    t->fxx[355]=0.0;
    t->fxx[356]=0.0;
    t->fxx[357]=0.0;
// j= 7
    t->fxx[358]=0.0;
    t->fxx[359]=0.0;
    t->fxx[360]=0.0;
    t->fxx[361]=0.0;
    t->fxx[362]=0.0;
    t->fxx[363]=0.0;
    t->fxx[364]=0.0;
    t->fxx[365]=0.0;
// j= 8
    t->fxx[366]=0.0;
    t->fxx[367]=0.0;
    t->fxx[368]=0.0;
    t->fxx[369]=0.0;
    t->fxx[370]=0.0;
    t->fxx[371]=0.0;
    t->fxx[372]=0.0;
    t->fxx[373]=0.0;
    t->fxx[374]=0.0;
// j= 9
    t->fxx[375]=0.0;
    t->fxx[376]=0.0;
    t->fxx[377]=0.0;
    t->fxx[378]=0.0;
    t->fxx[379]=0.0;
    t->fxx[380]=0.0;
    t->fxx[381]=0.0;
    t->fxx[382]=0.0;
    t->fxx[383]=0.0;
    t->fxx[384]=0.0;
// d^2f[7]/dx[i] dx[j]
// j= 0
    t->fxx[385]=0.0;
// j= 1
    t->fxx[386]=0.0;
    t->fxx[387]=0.0;
// j= 2
    t->fxx[388]=0.0;
    t->fxx[389]=0.0;
    t->fxx[390]=0.0;
// j= 3
    t->fxx[391]=0.0;
    t->fxx[392]=0.0;
    t->fxx[393]=0.0;
    t->fxx[394]=0.0;
// j= 4
    t->fxx[395]=0.0;
    t->fxx[396]=0.0;
    t->fxx[397]=0.0;
    t->fxx[398]=0.0;
    t->fxx[399]=0.0;
// j= 5
    t->fxx[400]=0.0;
    t->fxx[401]=0.0;
    t->fxx[402]=0.0;
    t->fxx[403]=0.0;
    t->fxx[404]=0.0;
    t->fxx[405]=0.0;
// j= 6
    t->fxx[406]=0.0;
    t->fxx[407]=0.0;
    t->fxx[408]=0.0;
    t->fxx[409]=0.0;
    t->fxx[410]=0.0;
    t->fxx[411]=0.0;
    t->fxx[412]=0.0;
// j= 7
    t->fxx[413]=0.0;
    t->fxx[414]=0.0;
    t->fxx[415]=0.0;
    t->fxx[416]=0.0;
    t->fxx[417]=0.0;
    t->fxx[418]=0.0;
    t->fxx[419]=0.0;
    t->fxx[420]=0.0;
// j= 8
    t->fxx[421]=0.0;
    t->fxx[422]=0.0;
    t->fxx[423]=0.0;
    t->fxx[424]=0.0;
    t->fxx[425]=0.0;
    t->fxx[426]=0.0;
    t->fxx[427]=0.0;
    t->fxx[428]=0.0;
    t->fxx[429]=0.0;
// j= 9
    t->fxx[430]=0.0;
    t->fxx[431]=0.0;
    t->fxx[432]=0.0;
    t->fxx[433]=0.0;
    t->fxx[434]=0.0;
    t->fxx[435]=0.0;
    t->fxx[436]=0.0;
    t->fxx[437]=0.0;
    t->fxx[438]=0.0;
    t->fxx[439]=0.0;
// d^2f[8]/dx[i] dx[j]
// j= 0
    t->fxx[440]=0.0;
// j= 1
    t->fxx[441]=0.0;
    t->fxx[442]=0.0;
// j= 2
    t->fxx[443]=0.0;
    t->fxx[444]=0.0;
    t->fxx[445]=0.0;
// j= 3
    t->fxx[446]=0.0;
    t->fxx[447]=0.0;
    t->fxx[448]=0.0;
    t->fxx[449]=0.0;
// j= 4
    t->fxx[450]=0.0;
    t->fxx[451]=0.0;
    t->fxx[452]=0.0;
    t->fxx[453]=0.0;
    t->fxx[454]=0.0;
// j= 5
    t->fxx[455]=0.0;
    t->fxx[456]=0.0;
    t->fxx[457]=0.0;
    t->fxx[458]=0.0;
    t->fxx[459]=0.0;
    t->fxx[460]=0.0;
// j= 6
    t->fxx[461]=0.0;
    t->fxx[462]=0.0;
    t->fxx[463]=0.0;
    t->fxx[464]=0.0;
    t->fxx[465]=0.0;
    t->fxx[466]=0.0;
    t->fxx[467]=0.0;
// j= 7
    t->fxx[468]=0.0;
    t->fxx[469]=0.0;
    t->fxx[470]=0.0;
    t->fxx[471]=0.0;
    t->fxx[472]=0.0;
    t->fxx[473]=0.0;
    t->fxx[474]=0.0;
    t->fxx[475]=0.0;
// j= 8
    t->fxx[476]=0.0;
    t->fxx[477]=0.0;
    t->fxx[478]=0.0;
    t->fxx[479]=0.0;
    t->fxx[480]=0.0;
    t->fxx[481]=0.0;
    t->fxx[482]=0.0;
    t->fxx[483]=0.0;
    t->fxx[484]=0.0;
// j= 9
    t->fxx[485]=0.0;
    t->fxx[486]=0.0;
    t->fxx[487]=0.0;
    t->fxx[488]=0.0;
    t->fxx[489]=0.0;
    t->fxx[490]=0.0;
    t->fxx[491]=0.0;
    t->fxx[492]=0.0;
    t->fxx[493]=0.0;
    t->fxx[494]=0.0;
// d^2f[9]/dx[i] dx[j]
// j= 0
    t->fxx[495]=0.0;
// j= 1
    t->fxx[496]=0.0;
    t->fxx[497]=0.0;
// j= 2
    t->fxx[498]=0.0;
    t->fxx[499]=0.0;
    t->fxx[500]=0.0;
// j= 3
    t->fxx[501]=0.0;
    t->fxx[502]=0.0;
    t->fxx[503]=0.0;
    t->fxx[504]=0.0;
// j= 4
    t->fxx[505]=0.0;
    t->fxx[506]=0.0;
    t->fxx[507]=0.0;
    t->fxx[508]=0.0;
    t->fxx[509]=0.0;
// j= 5
    t->fxx[510]=0.0;
    t->fxx[511]=0.0;
    t->fxx[512]=0.0;
    t->fxx[513]=0.0;
    t->fxx[514]=0.0;
    t->fxx[515]=0.0;
// j= 6
    t->fxx[516]=0.0;
    t->fxx[517]=0.0;
    t->fxx[518]=0.0;
    t->fxx[519]=0.0;
    t->fxx[520]=0.0;
    t->fxx[521]=0.0;
    t->fxx[522]=0.0;
// j= 7
    t->fxx[523]=0.0;
    t->fxx[524]=0.0;
    t->fxx[525]=0.0;
    t->fxx[526]=0.0;
    t->fxx[527]=0.0;
    t->fxx[528]=0.0;
    t->fxx[529]=0.0;
    t->fxx[530]=0.0;
// j= 8
    t->fxx[531]=0.0;
    t->fxx[532]=0.0;
    t->fxx[533]=0.0;
    t->fxx[534]=0.0;
    t->fxx[535]=0.0;
    t->fxx[536]=0.0;
    t->fxx[537]=0.0;
    t->fxx[538]=0.0;
    t->fxx[539]=0.0;
// j= 9
    t->fxx[540]=0.0;
    t->fxx[541]=0.0;
    t->fxx[542]=0.0;
    t->fxx[543]=0.0;
    t->fxx[544]=0.0;
    t->fxx[545]=0.0;
    t->fxx[546]=0.0;
    t->fxx[547]=0.0;
    t->fxx[548]=0.0;
    t->fxx[549]=0.0;

// d^2f[0]/du[i] du[j]
// j= 0
    t->fuu[0]=0.0;
// j= 1
    t->fuu[1]=0.0;
    t->fuu[2]=0.0;
// d^2f[1]/du[i] du[j]
// j= 0
    t->fuu[3]=0.0;
// j= 1
    t->fuu[4]=0.0;
    t->fuu[5]=0.0;
// d^2f[2]/du[i] du[j]
// j= 0
    t->fuu[6]=0.0;
// j= 1
    t->fuu[7]=0.0;
    t->fuu[8]=0.0;
// d^2f[3]/du[i] du[j]
// j= 0
// j= 1
    t->fuu[10]=0.0;
// d^2f[4]/du[i] du[j]
// j= 0
// j= 1
    t->fuu[13]=0.0;
// d^2f[5]/du[i] du[j]
// j= 0
// j= 1
    t->fuu[16]=0.0;
// d^2f[6]/du[i] du[j]
// j= 0
    t->fuu[18]=0.0;
// j= 1
    t->fuu[19]=0.0;
    t->fuu[20]=0.0;
// d^2f[7]/du[i] du[j]
// j= 0
    t->fuu[21]=0.0;
// j= 1
    t->fuu[22]=0.0;
    t->fuu[23]=0.0;
// d^2f[8]/du[i] du[j]
// j= 0
    t->fuu[24]=0.0;
// j= 1
    t->fuu[25]=0.0;
    t->fuu[26]=0.0;
// d^2f[9]/du[i] du[j]
// j= 0
    t->fuu[27]=0.0;
// j= 1
    t->fuu[28]=0.0;
    t->fuu[29]=0.0;

// d^2f[0]/dx[i] du[j]
// j= 0
    t->fxu[0]=0.0;
    t->fxu[1]=0.0;
    t->fxu[2]=0.0;
    t->fxu[3]=0.0;
    t->fxu[4]=0.0;
    t->fxu[5]=0.0;
    t->fxu[6]=0.0;
    t->fxu[7]=0.0;
    t->fxu[8]=0.0;
    t->fxu[9]=0.0;
// j= 1
    t->fxu[10]=0.0;
    t->fxu[11]=0.0;
    t->fxu[12]=0.0;
    t->fxu[13]=0.0;
    t->fxu[14]=0.0;
    t->fxu[15]=0.0;
    t->fxu[16]=0.0;
    t->fxu[17]=0.0;
    t->fxu[18]=0.0;
    t->fxu[19]=0.0;
// d^2f[1]/dx[i] du[j]
// j= 0
    t->fxu[20]=0.0;
    t->fxu[21]=0.0;
    t->fxu[22]=0.0;
    t->fxu[23]=0.0;
    t->fxu[24]=0.0;
    t->fxu[25]=0.0;
    t->fxu[26]=0.0;
    t->fxu[27]=0.0;
    t->fxu[28]=0.0;
    t->fxu[29]=0.0;
// j= 1
    t->fxu[30]=0.0;
    t->fxu[31]=0.0;
    t->fxu[32]=0.0;
    t->fxu[33]=0.0;
    t->fxu[34]=0.0;
    t->fxu[35]=0.0;
    t->fxu[36]=0.0;
    t->fxu[37]=0.0;
    t->fxu[38]=0.0;
    t->fxu[39]=0.0;
// d^2f[2]/dx[i] du[j]
// j= 0
    t->fxu[40]=0.0;
    t->fxu[41]=0.0;
    t->fxu[42]=0.0;
    t->fxu[43]=0.0;
    t->fxu[44]=0.0;
    t->fxu[45]=0.0;
    t->fxu[46]=0.0;
    t->fxu[47]=0.0;
    t->fxu[48]=0.0;
    t->fxu[49]=0.0;
// j= 1
    t->fxu[50]=0.0;
    t->fxu[51]=0.0;
    t->fxu[52]=0.0;
    t->fxu[53]=0.0;
    t->fxu[54]=0.0;
    t->fxu[55]=0.0;
    t->fxu[56]=0.0;
    t->fxu[57]=0.0;
    t->fxu[58]=0.0;
    t->fxu[59]=0.0;
// d^2f[3]/dx[i] du[j]
// j= 0
    t->fxu[60]=0.0;
    t->fxu[61]=0.0;
    t->fxu[62]=0.0;
    t->fxu[66]=0.0;
    t->fxu[67]=0.0;
    t->fxu[68]=0.0;
    t->fxu[69]=0.0;
// j= 1
    t->fxu[70]=0.0;
    t->fxu[71]=0.0;
    t->fxu[72]=0.0;
    t->fxu[76]=0.0;
    t->fxu[77]=0.0;
    t->fxu[78]=0.0;
    t->fxu[79]=0.0;
// d^2f[4]/dx[i] du[j]
// j= 0
    t->fxu[80]=0.0;
    t->fxu[81]=0.0;
    t->fxu[82]=0.0;
    t->fxu[86]=0.0;
    t->fxu[87]=0.0;
    t->fxu[88]=0.0;
    t->fxu[89]=0.0;
// j= 1
    t->fxu[90]=0.0;
    t->fxu[91]=0.0;
    t->fxu[92]=0.0;
    t->fxu[96]=0.0;
    t->fxu[97]=0.0;
    t->fxu[98]=0.0;
    t->fxu[99]=0.0;
// d^2f[5]/dx[i] du[j]
// j= 0
    t->fxu[100]=0.0;
    t->fxu[101]=0.0;
    t->fxu[102]=0.0;
    t->fxu[106]=0.0;
    t->fxu[107]=0.0;
    t->fxu[108]=0.0;
    t->fxu[109]=0.0;
// j= 1
    t->fxu[110]=0.0;
    t->fxu[111]=0.0;
    t->fxu[112]=0.0;
    t->fxu[116]=0.0;
    t->fxu[117]=0.0;
    t->fxu[118]=0.0;
    t->fxu[119]=0.0;
// d^2f[6]/dx[i] du[j]
// j= 0
    t->fxu[120]=0.0;
    t->fxu[121]=0.0;
    t->fxu[122]=0.0;
    t->fxu[123]=0.0;
    t->fxu[124]=0.0;
    t->fxu[125]=0.0;
    t->fxu[126]=0.0;
    t->fxu[127]=0.0;
    t->fxu[128]=0.0;
    t->fxu[129]=0.0;
// j= 1
    t->fxu[130]=0.0;
    t->fxu[131]=0.0;
    t->fxu[132]=0.0;
    t->fxu[133]=0.0;
    t->fxu[134]=0.0;
    t->fxu[135]=0.0;
    t->fxu[136]=0.0;
    t->fxu[137]=0.0;
    t->fxu[138]=0.0;
    t->fxu[139]=0.0;
// d^2f[7]/dx[i] du[j]
// j= 0
    t->fxu[140]=0.0;
    t->fxu[141]=0.0;
    t->fxu[142]=0.0;
    t->fxu[143]=0.0;
    t->fxu[144]=0.0;
    t->fxu[145]=0.0;
    t->fxu[146]=0.0;
    t->fxu[147]=0.0;
    t->fxu[148]=0.0;
    t->fxu[149]=0.0;
// j= 1
    t->fxu[150]=0.0;
    t->fxu[151]=0.0;
    t->fxu[152]=0.0;
    t->fxu[153]=0.0;
    t->fxu[154]=0.0;
    t->fxu[155]=0.0;
    t->fxu[156]=0.0;
    t->fxu[157]=0.0;
    t->fxu[158]=0.0;
    t->fxu[159]=0.0;
// d^2f[8]/dx[i] du[j]
// j= 0
    t->fxu[160]=0.0;
    t->fxu[161]=0.0;
    t->fxu[162]=0.0;
    t->fxu[163]=0.0;
    t->fxu[164]=0.0;
    t->fxu[165]=0.0;
    t->fxu[166]=0.0;
    t->fxu[167]=0.0;
    t->fxu[168]=0.0;
    t->fxu[169]=0.0;
// j= 1
    t->fxu[170]=0.0;
    t->fxu[171]=0.0;
    t->fxu[172]=0.0;
    t->fxu[173]=0.0;
    t->fxu[174]=0.0;
    t->fxu[175]=0.0;
    t->fxu[176]=0.0;
    t->fxu[177]=0.0;
    t->fxu[178]=0.0;
    t->fxu[179]=0.0;
// d^2f[9]/dx[i] du[j]
// j= 0
    t->fxu[180]=0.0;
    t->fxu[181]=0.0;
    t->fxu[182]=0.0;
    t->fxu[183]=0.0;
    t->fxu[184]=0.0;
    t->fxu[185]=0.0;
    t->fxu[186]=0.0;
    t->fxu[187]=0.0;
    t->fxu[188]=0.0;
    t->fxu[189]=0.0;
// j= 1
    t->fxu[190]=0.0;
    t->fxu[191]=0.0;
    t->fxu[192]=0.0;
    t->fxu[193]=0.0;
    t->fxu[194]=0.0;
    t->fxu[195]=0.0;
    t->fxu[196]=0.0;
    t->fxu[197]=0.0;
    t->fxu[198]=0.0;
    t->fxu[199]=0.0;
#endif
    }
    
    return 1;
}

static int init_final(trajFin_t *t, tOptSet *o) {
    double **p= o->p;
    const int k= o->n_hor;



    t->cx[6]=0.0;
    t->cx[7]=0.0;
    t->cx[8]=0.0;
    t->cx[9]=0.0;

    t->cxx[1]=0.0;
    t->cxx[3]=0.0;
    t->cxx[4]=0.0;
    t->cxx[6]=0.0;
    t->cxx[7]=0.0;
    t->cxx[8]=0.0;
    t->cxx[10]=0.0;
    t->cxx[11]=0.0;
    t->cxx[12]=0.0;
    t->cxx[13]=0.0;
    t->cxx[15]=0.0;
    t->cxx[16]=0.0;
    t->cxx[17]=0.0;
    t->cxx[18]=0.0;
    t->cxx[19]=0.0;
    t->cxx[21]=0.0;
    t->cxx[22]=0.0;
    t->cxx[23]=0.0;
    t->cxx[24]=0.0;
    t->cxx[25]=0.0;
    t->cxx[26]=0.0;
    t->cxx[27]=0.0;
    t->cxx[28]=0.0;
    t->cxx[29]=0.0;
    t->cxx[30]=0.0;
    t->cxx[31]=0.0;
    t->cxx[32]=0.0;
    t->cxx[33]=0.0;
    t->cxx[34]=0.0;
    t->cxx[35]=0.0;
    t->cxx[36]=0.0;
    t->cxx[37]=0.0;
    t->cxx[38]=0.0;
    t->cxx[39]=0.0;
    t->cxx[40]=0.0;
    t->cxx[41]=0.0;
    t->cxx[42]=0.0;
    t->cxx[43]=0.0;
    t->cxx[44]=0.0;
    t->cxx[45]=0.0;
    t->cxx[46]=0.0;
    t->cxx[47]=0.0;
    t->cxx[48]=0.0;
    t->cxx[49]=0.0;
    t->cxx[50]=0.0;
    t->cxx[51]=0.0;
    t->cxx[52]=0.0;
    t->cxx[53]=0.0;
    t->cxx[54]=0.0;

    return 1;
}

int init_trajectory(traj_t *t, tOptSet *o) {
    if(!init_running(t->t, o)) return 0;
    if(!init_final(&t->f, o)) return 0;
    
    return 1;
}

static int init_multipliers_running(tOptSet *o) {
    multipliersEl_t *m= o->multipliers.t;
    int k, i;

    for(k= 0; k<o->n_hor; k++, m++) {

    }
    
    return 1;
}

static int init_multipliers_final(tOptSet *o) {
    multipliersFin_t *m= &o->multipliers.f;
    int i;


    
    return 1;
}

int init_multipliers(tOptSet *o) {
    if(!init_multipliers_running(o)) return 0;
    if(!init_multipliers_final(o)) return 0;
    
    return 1;
}

int init_opt(tOptSet *o) {
    int i;
    
    for(i= 0; i<NUMBER_OF_THREADS+1; i++)
        if(!init_trajectory(&o->trajectories[i], o)) return 0;

    o->nominal= &o->trajectories[0];
    for(i= 1; i<NUMBER_OF_THREADS+1; i++)
    o->candidates[i-1]= &o->trajectories[i];
    
    if(!init_multipliers(o)) return 0;
    
    return 1;
}

static int update_multipliers_running(tOptSet *o, int init) {
    trajEl_t *t= o->nominal->t;
    multipliersEl_t *m= o->multipliers.t;
    const double *x= t->x;
    const double *u= t->u;
    const double w_pen= o->w_pen_l;
    double **p= o->p;
    int increase_pen= 0;
    int k, i;
    
    for(k= 0; k<o->n_hor; k++, m++, t++) {

// TODO: maybe check test for sufficient complementarity reduction

    if(init) return 1;

// inequality constraints according to D. Ruxton: Differential dynamic programming applied to continuous optimal control problems with state variable inequality constraints
    }
    
    if(!init && increase_pen)
        o->w_pen_l= min(o->w_pen_max_l, o->w_pen_l*o->w_pen_fact1);
        
    return 1;
}

static int update_multipliers_final(tOptSet *o, int init) {
    trajFin_t *t= &o->nominal->f;
    multipliersFin_t *m= &o->multipliers.f;
    const double *x= t->x;
    const double w_pen= o->w_pen_f;
    double **p= o->p;
    int increase_pen= 0;
    int k= o->n_hor;


// TODO: maybe check test for sufficient complementarity reduction

    if(!init && increase_pen)
        o->w_pen_f= min(o->w_pen_max_f, o->w_pen_f*o->w_pen_fact1);

    if(init) return 1;

// inequality constraints according to D. Ruxton: Differential dynamic programming applied to continuous optimal control problems with state variable inequality constraints

    return 1;
}

int update_multipliers(tOptSet *o, int init) {
    if(!update_multipliers_running(o, init)) return 0;
    if(!update_multipliers_final(o, init)) return 0;

    return 1;
}

int get_g_size() {
    return(0);
}

int calcG(double g[], trajEl_t *t, int k, double **p) {
    const double *x= t->x;
    const double *u= t->u;
    
    return(1);
} 
