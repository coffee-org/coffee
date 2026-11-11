/**
 * @file    AOsystSim_internal.h
 * @brief   Internal shared state, context structures, and subsystem helpers for AOsystSim
 */

#ifndef AOSYSTSIM_INTERNAL_H
#define AOSYSTSIM_INTERNAL_H

#include <stdint.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include "milk_runnable.h"
#include "coffee_compat.h"
#include "OptSystProp/OptSystProp.h"

#define AOSIM_INTERNAL __attribute__((visibility("hidden")))

/* --- Shared Module State & Working Directory Helper (AOsystSim.c) --- */

extern OPTSYST *optsystsim;

/* Create AOsystSim_wdir working directory on first invocation. */
AOSIM_INTERNAL void aosim_ensure_wdir(void);

/* --- Atmospheric Wavefront Generator Context (AOsystSim_wf_conf.c / AOsystSim_wf.c) --- */

/**
 * @brief Internal runtime state for multi-layer atmospheric wavefront streaming.
 */
typedef struct
{
    char    WFDIR[STRINGMAXLEN_DIRNAME];
    float   LAMBDA;
    char    conf_atmWFconf_fname[STRINGMAXLEN_FULLFILENAME];
    int     OUT0FITSFILE;
    char    OUT0FITSFILENAMEOPD[STRINGMAXLEN_FULLFILENAME];
    char    OUT0FITSFILENAMEAMP[STRINGMAXLEN_FULLFILENAME];
    char    OUTPHYSTIME[STRINGMAXLEN_STREAMNAME];
    int     OUT0STREAM;
    char    OUT0STREAMNAMEOPD[STRINGMAXLEN_STREAMNAME];
    char    OUT0STREAMNAMEAMP[STRINGMAXLEN_STREAMNAME];
    int     TRIGGERMODE;
    char    TRIGGERFILE[STRINGMAXLEN_FULLFILENAME];
    float   TRIGGERDT;
    char    TRIGGER0STREAM[STRINGMAXLEN_STREAMNAME];
    int     TRIGGER0SEM;
    char    TRIGGER1STREAM[STRINGMAXLEN_STREAMNAME];
    int     TRIGGER1SEM;
    float   DT;
    float   OPDMFACT;
    float   AMPMFACT;
    long    ARRAYSIZE;
    int     PIXSCALEMODE;
    int     PIXBINFACTOR;
    float   PUPDIAM;
    char    PUPILFILE[STRINGMAXLEN_FULLFILENAME];
    int     DM0MODE;
    char    DM0NAME[STRINGMAXLEN_STREAMNAME];
    char    OUT1FITSFILENAMEOPD[STRINGMAXLEN_FULLFILENAME];
    char    OUT1FITSFILENAMEAMP[STRINGMAXLEN_FULLFILENAME];
    int     OUT1STREAM;
    char    OUT1STREAMNAMEOPD[STRINGMAXLEN_STREAMNAME];
    char    OUT1STREAMNAMEAMP[STRINGMAXLEN_STREAMNAME];
    float   wfin_TIME_STEP;
    long    wfin_NB_TSPAN;
    float   wfin_TIME_SPAN;
    char    wfin_PREFIX[200];
    float   wfin_PUPIL_SCALE;
    long    wfin_WFsize;
    long    NBframes;
    long    kmax;
    int     AMPfile;
    float   pupscale;
    long    iioffset;
    long    jjoffset;
    long    ii1start0;
    long    jj1start0;
    long    ii1end0;
    long    jj1end0;
    double  phystime;
    imageID IDphystime;
    imageID IDopd0;
    imageID IDamp0;
    imageID IDopd1;
    imageID IDamp1;
    imageID IDdm0opd;
    imageID IDampmask;
    imageID IDbin_re;
    imageID IDbin_im;
    imageID IDbin_amp;
    imageID IDbin_opd;
} AOSIM_MKWF_CTX;

/* Parse atmospheric wavefront configuration files into context structure. */
AOSIM_INTERNAL void aosim_mkwf_read_conf(
    const char     *conf_fname,
    AOSIM_MKWF_CTX *ctx);

/* Create output shared-memory wavefront streams and connect optional DM0 stream. */
AOSIM_INTERNAL void aosim_mkwf_create_streams(AOSIM_MKWF_CTX *ctx);

/* Scan input WF FITS cubes and build telescope pupil amplitude mask. */
AOSIM_INTERNAL void aosim_mkwf_init_pupil(AOSIM_MKWF_CTX *ctx);

/* --- Coronagraph & Lyot LOWFS Context (AOsystSim_lowfs_init.c / AOsystSim_lowfs.c) --- */

/**
 * @brief Internal runtime state for coronagraph and Lyot LOWFS simulation loop.
 */
typedef struct
{
    int      INMODE;
    char     INAMPSTREAMNAME[200];
    char     INOPDSTREAMNAME[200];
    char     INAMPFITSFILENAME[STRINGMAXLEN_FULLFILENAME];
    char     INOPDFITSFILENAME[STRINGMAXLEN_FULLFILENAME];
    char     INTRIGSTREAMNAME[200];
    int      INTRIGSEMCHAN;
    char     INTRIGGERFILE[200];
    char     INPHYSTIME[200];
    int      OUTMODE;
    uint32_t OUTLOWFSARRAYSIZE;
    uint64_t OUTLOWFSARRAYSIZE2;
    char     OUTLOWFSSTREAMNAME[200];
    char     OUTLOWFSFITSFILENAME[STRINGMAXLEN_FULLFILENAME];
    char     OUTTRIGGERFILE[200];
    char     COROFPMAMP[200];
    char     COROFPMPHA[200];
    char     LYOTSTOPREFLAMP[200];
    char     LYOTSTOPREFLPHA[200];
    char     LYOTSTOPTRANSMAMP[200];
    char     LYOTSTOPTRANSMPHA[200];
    char     LOWFSOPDMAP[200];
    float    LOWFSCAMETIME;
    float    LOWFScamstarttime;
    long     LOWFScamcnt;
    uint32_t ARRAYSIZE;
    uint64_t ARRAYSIZE2;
    double   LAMBDA;
    uint32_t xsizein;
    uint32_t ysizein;
    long     iioffset;
    long     jjoffset;
    imageID  IDoutLOWFS;
    imageID  IDinTRIG;
    imageID  IDinOPD;
    imageID  IDinAMP;
    imageID  IDphystime;
    imageID  IDfpmamp;
    imageID  IDfpmpha;
    imageID  ID_LS_ramp;
    imageID  ID_LS_rpha;
    imageID  ID_LS_tamp;
    imageID  ID_LS_tpha;
    imageID  ID_lowfsopd;
    imageID  IDwfa;
    imageID  IDwfp;
    imageID  IDfoc1a;
    imageID  IDfoc1p;
    imageID  IDpup1ta;
    imageID  IDpup1tp;
    imageID  IDpup1ra;
    imageID  IDpup1rp;
    imageID  IDimlowfs;
    imageID  IDimcamlowfstmp;
} AOSIM_LOWFS_CTX;

/* Read coronagraph and LOWFS daemon configuration parameters. */
AOSIM_INTERNAL void aosim_lowfs_read_conf(
    const char      *conf_fname,
    AOSIM_LOWFS_CTX *ctx);

/* Attach input wavefront streams and validate preloaded optical mask images. */
AOSIM_INTERNAL void aosim_lowfs_connect_inputs(AOSIM_LOWFS_CTX *ctx);

/* Allocate output LOWFS streams and internal complex-plane working images. */
AOSIM_INTERNAL void aosim_lowfs_init_buffers(AOSIM_LOWFS_CTX *ctx);

/* --- Focal-Plane Wavefront Sensing State & Helpers (AOsystSim_fpwfs_*.c) --- */

extern long        NBprobesG;
extern int         CENTERprobe;
extern long        NBoptVar;
extern long double probe_re[100];
extern long double probe_im[100];
extern long double nprobe_re[100];
extern long double nprobe_im[100];
extern long double probe_tflux[100];
extern long double Cflux;
extern long double probe_nmflux[100];
extern long double probe_nmnoise[100];
extern double      tmpvalue1;
extern double      tmpvalue2;

/* Evaluate GSL simplex chi-squared objective for FPWFS complex-amplitude estimation. */
double f_eval(
    const gsl_vector *v,
    void             *params);

/**
 * @brief Internal state for FPWFS Monte Carlo sensitivity analysis.
 */
typedef struct
{
    int                      mapmode;
    int                      mode;
    int                      optmode;
    int                      NBprobes;
    int                      execmode;
    int                      imSimulMode;
    int                      NewWF;
    int                      initmap;
    double                   ProbeNoise;
    double                   probe_noise_prop;
    double                   probeamp;
    double                   Iflux;
    double                   probeampl;
    double                   WFerr;
    double                   contrast;
    double                   totFlux;
    double                   FLUXph;
    double                   RON;
    double                   CnoiseFloor;
    double                   mapampl;
    long                     mapsize;
    long                     mapxsize;
    long                     mapysize;
    long                     mapzsize;
    long                     kxtest;
    long                     kytest;
    long                     iterMax;
    imageID                  IDpsfC;
    imageID                  IDpsfCnoise;
    imageID                  IDprobampC;
    imageID                  IDmap;
    imageID                  IDmap_ptre;
    imageID                  IDmap_ptim;
    imageID                  IDmap_ptre_in;
    imageID                  IDmap_ptim_in;
    imageID                  IDmap_Iflux;
    imageID                  IDmap_are;
    imageID                  IDmap_aim;
    imageID                  IDmap_e;
    gsl_multimin_fminimizer *s;
    gsl_vector              *ss;
    gsl_vector              *xvect;
    gsl_multimin_function    feval_func;
} AOSIM_FPWFS_SENS_CTX;

/* Solve FPWFS electric-field estimation across all pixels of one Monte Carlo slice. */
AOSIM_INTERNAL void aosim_fpwfs_solve_slice(
    AOSIM_FPWFS_SENS_CTX *ctx,
    long                  mapz);

#endif /* AOSYSTSIM_INTERNAL_H */
