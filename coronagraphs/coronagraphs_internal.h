/**
 * @file    coronagraphs_internal.h
 * @brief   Internal state, constants, and cross-file prototypes for coffeecoronagraphs
 */

#ifndef CORONAGRAPHS_INTERNAL_H
#define CORONAGRAPHS_INTERNAL_H

#include <stdio.h>

#include "milk_runnable.h"
#include "coronagraphs.h"

#define CORO_INTERNAL __attribute__((visibility("hidden")))

#define CORO_PI               3.14159265358979323846264338328
#define CORONAGRAPHSDATALOCAL "/data/tmp/coronagraphs"
#define CORONAGRAPHSDATADIR   ""
#define CORONAGRAPHS_PSCALE   0.006
#define CORONAGRAPHS_TDIAM    2.4
#define CORONAGRAPHS_LAMBDA   0.0000005

#define BL8MASK_NBSTEP 100000
#define BL8MASK_STEP   0.01

#define CPAAPO_NBPOINTS 1930
#define CPAAPO_FNAME    "pup_10.411256_2.0"
#define CPAFPMASKRAD    0.0
#define CPAPPMASKRAD    1.0

#define PIAAPPMASKRAD        1.0
#define PIAACPPMASKRAD1      1.0
#define PIAA_HYBRID_CST      0.00000001
#define piaaconfNPUPFILESIZE 200000
#define piaaconfNBr1fr0_pts  20000
#define PIAAAIC_FIELDLIMIT   15.0

extern double CORONAGRAPHS_PIXSCALE;
extern double COROTMP1;
extern double MASKSIZELD;

typedef struct coro_state
{
    int     useDFT;
    double  DFTZFACTOR;
    int     FPMASKSIZE_ERROR;
    double  FPMASK_FACTOR;
    double  FPMASK_FACTOR1;
    double  FPMASK_FACTOR2;
    double  FPM_TRANSM1;
    double  FPM_TRANSM2;
    double  BL8MASK_m;
    double  BL8MASK_l;
    double  BL8MASK_eps;
    int     BL8MODE;
    double  BL4MASK_eps;
    double  ODC_GAUSS;
    double  ODCMASK_eps;
    double  SHEAR4_SHEAR;
    double  RRPM_RADIUS;
    double  RRPM_P2;
    double  RRPM_P3;
    double  RRPM_P4;
    double  RRPM_P5;
    int     RRPM_PIAA;
    double  PIAACENTOBS;
    long    PIAAAPO_NBPOINTS;
    double *PIAAAPO;
    double *PIAA_HYBRID_CPAAPO;
    double  PIAAAPO2D[CORONAGRAPHS_ARRAYSIZE * CORONAGRAPHS_ARRAYSIZE];
    char    PIAAAPO_FNAME[400];
    char    PIAAAPODIZE_2DAPOFNAME[400];
    int     PIAAFPMASK;
    int     PIAALOWFS;
    double  PIAAFPMASKRAD;
    double *piaaconfpup_amp_profile;
    double *piaaconfr0;
    double *piaaconfr1;
    long    piaaconfNBpoints;
    double *piaaconfr0fr1;
    double *piaaconfr1fr0;
    double *piaaconfM0;
    double *piaaconfM1;
    int     piaaconfdirection;
    double  PIAAOVERSIZE;
    int     initPIAA;
    int     AUTOPIAACMASK;
    double  APLC_CentOBS0;
    double  APLC_CentOBS1;
    double  PIAAextfactor0;
    double  PIAAextfactor1;
    int     APLC_PIAA;
    long    NB_APLC_STEP;
    double  APLC_FPMASKsize;
    int     APLC_FLIP;
    int     APLC_PMASK;
    double  FPMASK_transm_error;
    double  FPMASK_size_error;
    double  STRIPCOFFSET;
    long    OVC_CHARGE;
    long    aporawN;
    double *aporaw_r;
    double *aporaw_v;
    int     fitapoINIT;
    long    fitapoN;
    double *fitapo_a;
    double *fitapo_b;
    double *fitapo_c;
    double *fitapo_c1;
    long    LOOPCNT;
    double  optval0;
    long    IDprol_init;
    long    IDprol_ffrac;
    long    IDprol_transm;
    long    IDprol_peak;
    long    IDprol_fitapo_a;
    long    IDprol_fitapo_b;
    long    IDprol_fitapo_c;
    long    IDprol_fitfit;
    double  fitapo_minc;
    double  APLCapo_CO_START;
    double  APLCapo_CO_END;
    double  APLCapo_CO_STEP;
    double  APLCapo_FPMRAD_START;
    double  APLCapo_FPMRAD_END;
    double  APLCapo_FPMRAD_STEP;
} coro_state_t;

extern CORO_INTERNAL coro_state_t g_coro;

typedef struct coro_ref_pair
{
    const char *im_re;
    const char *file_re;
    const char *im_im;
    const char *file_im;
    const char *label;
    imageID     id_re;
    imageID     id_im;
    int         needs_init;
} coro_ref_pair_t;

typedef struct coro_scan_ctx
{
    double  obs0input;
    double  obs0;
    long    kmax;
    long    size;
    double  pixscale;
    double  starprof_step;
    double  starprof_maxrad;
    long    starprofNBstep;
    double  star_radius;
    double *separray;
    double *transmarray;
    double *starprofrad;
    double *starprofarray;
    double *starprofarraycnt;
    float  *obs1array;
    long    obs1arraysize;
    long    obs1_index_start;
    float  *fpmaskradarray;
    long    fpmaskradarraysize;
    int     compute_data;
    int     cent_obs_flag;
} coro_scan_ctx_t;

/* Shared internal utilities in coronagraphs.c */
CORO_INTERNAL double coro_pupil_radius_pix(void);

CORO_INTERNAL errno_t coro_check_fscanf(
    FILE *fp,
    int   fscanfcnt,
    int   expected);

CORO_INTERNAL errno_t coro_fft2d_centered(
    const char *in_name,
    const char *out_name);

CORO_INTERNAL void coro_setup_pupil_tiptilt(
    double      xld,
    double      yld,
    int         subpix,
    double      oversize,
    const char *pha_name,
    double     *out_total);

CORO_INTERNAL void coro_complex_to_psf(
    const char *cname,
    const char *fname,
    const char *aname,
    const char *pname,
    const char *psfname,
    double      norm_denom);

CORO_INTERNAL void coro_load_ref_pair(coro_ref_pair_t *ref);

CORO_INTERNAL void coro_apply_or_save_ref_pair(
    coro_ref_pair_t *ref,
    imageID          id_c);

/* Internal PIAA, APLC, and scan helpers */
CORO_INTERNAL errno_t coro_piaa_build_lookup_and_shapes(long nb_points);

CORO_INTERNAL void coro_piaa_inverse_unwrap(
    const char *in_amp,
    const char *in_pha,
    const char *out_amp,
    const char *out_pha);

CORO_INTERNAL void coro_aplc_propagate_steps(
    int    write_files,
    double mask_ca_mult);

CORO_INTERNAL void coro_aplc_finish_psf(
    const char *psfname,
    double      total0,
    int         write_files);

CORO_INTERNAL errno_t coro_scan_compute_grid(coro_scan_ctx_t *ctx);

/* Module functions exported across coffeecoronagraphs */
double coronagraph_make_2Dprolate_DFT(
    double      fpmradpix,
    double      beamradpix,
    double      centralObs,
    const char *outname,
    long        size,
    const char *pupmask_name,
    double     *transmval);

errno_t coronagraph_apofit_eval_r(
    double  r,
    double *outval);

errno_t coronagraph_apofit_eval(double *outval);

double coronagraph_apofit(const char *fnameout);

int coronagraphs_free_PIAA(void);

int coronagraph_telescope_pupil_Subaru_inside1(
    double x,
    double y);

int coronagraph_PIAAperturbation(
    double *zarray,
    long   *zindex,
    long    NBzern,
    double  ratio);

int coronagraph_init_CPA(void);
int coronagraph_init_ODC(void);
int coronagraph_init_BL8(void);
int coronagraph_init_BL4(void);
int coronagraph_init_RRPM(void);
int coronagraph_init_OVC(long charge);

int coronagraph_simul_SHEAR4(
    double      xld,
    double      yld,
    const char *psfname);

int coronagraph_simul_DICC(
    double      xld,
    double      yld,
    const char *psfname);

int coronagraph_simul_ODC(
    double      xld,
    double      yld,
    const char *psfname);

int coronagraph_simul_BL4(
    double      xld,
    double      yld,
    const char *psfname);

int coronagraph_simul_PPA(
    double      xld,
    double      yld,
    const char *psfname);

int coronagraph_simul_NOCORO(
    double      xld,
    double      yld,
    const char *psfname);

int coronagraph_simul_STRIPC(
    double      xld,
    double      yld,
    const char *psfname);

int coronagraph_simul_SIMXY(
    double      xld,
    double      yld,
    const char *psfname);

int coronagraph_RRPM_optimize(void);

int coronagraph_init_EXTERNAL_OCCULTER(
    double D,
    double l,
    double lambda,
    long   FACTOR);

int coronagraph_simul_EXTERNAL_OCCULTER(
    double      xld,
    double      yld,
    const char *psfname);

#endif /* CORONAGRAPHS_INTERNAL_H */
