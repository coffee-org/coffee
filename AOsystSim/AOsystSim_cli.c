/**
 * @file    AOsystSim_cli.c
 * @brief   MILK_RUNNABLE command wrappers and module registry for coffeeAOsystSim
 */

#define MODULE_SHORTNAME_DEFAULT "aosystsim"
#define MODULE_DESCRIPTION       "Adaptive Optics system simulator"

#include "milk_runnable.h"
#include "OptSystProp/OptSystProp.h"
#include "AOsystSim.h"

MODULE_DEPS("milkinfo", "milkimagegen", "milkimagefilter", "milkstatistic",
            "milkfft", "milkOpticsMaterials", "milkWFpropagate", "coffeeOptSystProp");

/* 1. AOsimfilt */
static char *filt_inname;
static char *filt_outname;

static CLICMDARGDEF farg_aosimfilt[] = {
    { CLIARG_STR_NOT_IMG, ".inwf", "input wavefront image name", "wfin",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &filt_inname, NULL },
    { CLIARG_STR_NOT_IMG, ".outwf", "output wavefront image name", "wfout",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &filt_outname, NULL },
};

/**
 * cli_aosim_filter - Run Fourier-domain open-loop AO spatial/temporal filter
 *
 * Return: RETURN_SUCCESS on completion.
 */
static errno_t cli_aosim_filter(void)
{
    AOsystSim_simpleAOfilter(filt_inname, filt_outname);
    return RETURN_SUCCESS;
}

MILK_RUNNABLE_RT_DEFINE(AOsimfilt);
static const MILK_RUNNABLE runnable_aosim_filt = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "AOsimfilt", .description = "AO simple Fourier filtering" },
    .srcfile  = __FILE__,
    .farg     = farg_aosimfilt,
    .nbarg    = (int) (sizeof(farg_aosimfilt) / sizeof(farg_aosimfilt[0])),
    .cli_func = cli_aosim_filter,
    .rt       = &AOsimfilt_rt,
};

/* 2. AOsystsfitpup */
static char *fitpup_dmim;
static char *fitpup_outname;

static CLICMDARGDEF farg_fitpup[] = {
    { CLIARG_IMG, ".dmresp", "measured DM response image", "dmresp",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &fitpup_dmim, NULL },
    { CLIARG_STR_NOT_IMG, ".telpup", "output fitted telescope pupil", "telpup",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &fitpup_outname, NULL },
};

/**
 * cli_aosim_fitpup - Fit telescope pupil geometry to DM response image
 *
 * Return: RETURN_SUCCESS on completion.
 */
static errno_t cli_aosim_fitpup(void)
{
    AOsystSim_fitTelPup(fitpup_dmim, fitpup_outname);
    return RETURN_SUCCESS;
}

MILK_RUNNABLE_RT_DEFINE(AOsystsfitpup);
static const MILK_RUNNABLE runnable_aosim_fitpup = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "AOsystsfitpup", .description = "fit telescope pupil" },
    .srcfile  = __FILE__,
    .farg     = farg_fitpup,
    .nbarg    = (int) (sizeof(farg_fitpup) / sizeof(farg_fitpup[0])),
    .cli_func = cli_aosim_fitpup,
    .rt       = &AOsystsfitpup_rt,
};

/* 3. AOsimmkWF */
static char *mkwf_conf;

static CLICMDARGDEF farg_mkwf[] = {
    { CLIARG_STR, ".conffile", "WF configuration filename", "WF.conf",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &mkwf_conf, NULL },
};

/**
 * cli_aosim_mkWF - Generate atmospheric wavefront series from configuration file
 *
 * Return: RETURN_SUCCESS on completion.
 */
static errno_t cli_aosim_mkWF(void)
{
    AOsystSim_mkWF(mkwf_conf);
    return RETURN_SUCCESS;
}

MILK_RUNNABLE_RT_DEFINE(AOsimmkWF);
static const MILK_RUNNABLE runnable_aosim_mkwf = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "AOsimmkWF", .description = "make WF series for AO simulation" },
    .srcfile  = __FILE__,
    .farg     = farg_mkwf,
    .nbarg    = (int) (sizeof(farg_mkwf) / sizeof(farg_mkwf[0])),
    .cli_func = cli_aosim_mkWF,
    .rt       = &AOsimmkWF_rt,
};

/* 4. AOsimPyrWFS */
static char *pyrwfs_conf;

static CLICMDARGDEF farg_pyrwfs[] = {
    { CLIARG_STR, ".conffile", "Pyramid WFS configuration filename", "PyrWFS.conf",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &pyrwfs_conf, NULL },
};

/**
 * cli_aosim_PyrWFS - Run Pyramid WFS simulation from configuration file
 *
 * Return: RETURN_SUCCESS on completion.
 */
static errno_t cli_aosim_PyrWFS(void)
{
    AOsystSim_PyrWFS(pyrwfs_conf);
    return RETURN_SUCCESS;
}

MILK_RUNNABLE_RT_DEFINE(AOsimPyrWFS);
static const MILK_RUNNABLE runnable_aosim_pyrwfs = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "AOsimPyrWFS", .description = "run pyramid WFS simulation" },
    .srcfile  = __FILE__,
    .farg     = farg_pyrwfs,
    .nbarg    = (int) (sizeof(farg_pyrwfs) / sizeof(farg_pyrwfs[0])),
    .cli_func = cli_aosim_PyrWFS,
    .rt       = &AOsimPyrWFS_rt,
};

/* 5. AOsimDM */
static char *dm_conf;

static CLICMDARGDEF farg_dm[] = {
    { CLIARG_STR, ".conffile", "DM configuration filename", "DM.conf",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &dm_conf, NULL },
};

/**
 * cli_aosim_DM - Run deformable mirror simulation from configuration file
 *
 * Return: RETURN_SUCCESS on completion.
 */
static errno_t cli_aosim_DM(void)
{
    AOsystSim_DM(dm_conf);
    return RETURN_SUCCESS;
}

MILK_RUNNABLE_RT_DEFINE(AOsimDM);
static const MILK_RUNNABLE runnable_aosim_dm = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "AOsimDM", .description = "run DM for AO simulation" },
    .srcfile  = __FILE__,
    .farg     = farg_dm,
    .nbarg    = (int) (sizeof(farg_dm) / sizeof(farg_dm[0])),
    .cli_func = cli_aosim_DM,
    .rt       = &AOsimDM_rt,
};

/* 6. AOsimcoroLOWFS */
static char *lowfs_conf;

static CLICMDARGDEF farg_lowfs[] = {
    { CLIARG_STR, ".conffile", "Coronagraph/LOWFS configuration filename", "LOWFS.conf",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &lowfs_conf, NULL },
};

/**
 * cli_aosim_coroLOWFS - Run coronagraph and LOWFS simulation from configuration file
 *
 * Return: RETURN_SUCCESS on completion.
 */
static errno_t cli_aosim_coroLOWFS(void)
{
    AOsystSim_coroLOWFS(lowfs_conf);
    return RETURN_SUCCESS;
}

MILK_RUNNABLE_RT_DEFINE(AOsimcoroLOWFS);
static const MILK_RUNNABLE runnable_aosim_lowfs = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "AOsimcoroLOWFS", .description = "run coro and LOWFS simulation" },
    .srcfile  = __FILE__,
    .farg     = farg_lowfs,
    .nbarg    = (int) (sizeof(farg_lowfs) / sizeof(farg_lowfs[0])),
    .cli_func = cli_aosim_coroLOWFS,
    .rt       = &AOsimcoroLOWFS_rt,
};

/* 7. AOsystsim */
static long *sim_syncmode;
static long *sim_dmindex;
static long *sim_delayus;

static CLICMDARGDEF farg_aosystsim[] = {
    { CLIARG_INT64, ".syncmode", "synchronization mode (0:turb, 1:DM, 2:both, >2:delay)", "0",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &sim_syncmode, NULL },
    { CLIARG_INT64, ".dmindex", "DM stream index", "0",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &sim_dmindex, NULL },
    { CLIARG_INT64, ".delayus", "loop delay [us] when syncmode > 2", "1000",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &sim_delayus, NULL },
};

/**
 * cli_aosim_run - Run coupled DM-to-WFS AO simulation loop
 *
 * Return: RETURN_SUCCESS on completion.
 */
static errno_t cli_aosim_run(void)
{
    AOsystSim_run((int) *sim_syncmode, *sim_dmindex, *sim_delayus);
    return RETURN_SUCCESS;
}

MILK_RUNNABLE_RT_DEFINE(AOsystsim);
static const MILK_RUNNABLE runnable_aosim_run = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "AOsystsim", .description = "run coupled AO system simulator" },
    .srcfile  = __FILE__,
    .farg     = farg_aosystsim,
    .nbarg    = (int) (sizeof(farg_aosystsim) / sizeof(farg_aosystsim[0])),
    .cli_func = cli_aosim_run,
    .rt       = &AOsystsim_rt,
};

/* 8. AOsystexaosim */
MILK_RUNNABLE_RT_DEFINE(AOsystexaosim);
static const MILK_RUNNABLE runnable_aosim_exao = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "AOsystexaosim", .description = "run analytical extreme-AO analysis" },
    .srcfile  = __FILE__,
    .cli_func = AOsystSim_extremeAO_contrast_sim,
    .rt       = &AOsystexaosim_rt,
};

/* 9. AOsystmkABprobes */
static char   *prb_nameA;
static char   *prb_nameB;
static long   *prb_dmx;
static long   *prb_dmy;
static double *prb_cpamax;
static double *prb_cparmin;
static double *prb_cparmax;
static double *prb_rmsampl;
static long   *prb_geom;

static CLICMDARGDEF farg_mkprobes[] = {
    { CLIARG_STR_NOT_IMG, ".probeA", "output probe A image name", "prA",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &prb_nameA, NULL },
    { CLIARG_STR_NOT_IMG, ".probeB", "output probe B image name", "prB",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &prb_nameB, NULL },
    { CLIARG_INT64, ".dmxsize", "DM X actuator count", "50",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &prb_dmx, NULL },
    { CLIARG_INT64, ".dmysize", "DM Y actuator count", "50",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &prb_dmy, NULL },
    { CLIARG_FLOAT64, ".CPAmax", "maximum spatial frequency [CPA]", "0.7",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &prb_cpamax, NULL },
    { CLIARG_FLOAT64, ".CPArmin", "inner spatial frequency [CPA]", "0.1",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &prb_cparmin, NULL },
    { CLIARG_FLOAT64, ".CPArmax", "outer spatial frequency [CPA]", "0.7",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &prb_cparmax, NULL },
    { CLIARG_FLOAT64, ".RMSampl", "probe RMS amplitude", "0.1",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &prb_rmsampl, NULL },
    { CLIARG_INT64, ".modegeom", "probe geometry mode", "1",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &prb_geom, NULL },
};

/**
 * cli_aosim_mkprobes - Generate pair of focal-plane wavefront sensing DM probes
 *
 * Return: RETURN_SUCCESS on completion.
 */
static errno_t cli_aosim_mkprobes(void)
{
    AOsystSim_FPWFS_mkprobes(
        prb_nameA,
        prb_nameB,
        *prb_dmx,
        *prb_dmy,
        *prb_cpamax,
        *prb_cparmin,
        *prb_cparmax,
        *prb_rmsampl,
        *prb_geom);
    return RETURN_SUCCESS;
}

MILK_RUNNABLE_RT_DEFINE(AOsystmkABprobes);
static const MILK_RUNNABLE runnable_aosim_mkprobes = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "AOsystmkABprobes", .description = "make AB focal-plane WFS probes" },
    .srcfile  = __FILE__,
    .farg     = farg_mkprobes,
    .nbarg    = (int) (sizeof(farg_mkprobes) / sizeof(farg_mkprobes[0])),
    .cli_func = cli_aosim_mkprobes,
    .rt       = &AOsystmkABprobes_rt,
};

/* 10. AOsystFPWFSan */
static long *fpwfs_mapmode;
static long *fpwfs_mode;
static long *fpwfs_optmode;
static long *fpwfs_nbprobes;

static CLICMDARGDEF farg_fpwfsan[] = {
    { CLIARG_INT64, ".mapmode", "sensitivity map mode", "0",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &fpwfs_mapmode, NULL },
    { CLIARG_INT64, ".mode", "sensing mode", "0",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &fpwfs_mode, NULL },
    { CLIARG_INT64, ".optmode", "optimization mode", "0",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &fpwfs_optmode, NULL },
    { CLIARG_INT64, ".NBprobes", "number of probes", "3",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &fpwfs_nbprobes, NULL },
};

/**
 * cli_aosim_fpwfsan - Run focal-plane WFS coherence and sensitivity analysis
 *
 * Return: RETURN_SUCCESS on completion.
 */
static errno_t cli_aosim_fpwfsan(void)
{
    AOsystSim_FPWFS_sensitivityAnalysis(
        (int) *fpwfs_mapmode,
        (int) *fpwfs_mode,
        (int) *fpwfs_optmode,
        (int) *fpwfs_nbprobes);
    return RETURN_SUCCESS;
}

MILK_RUNNABLE_RT_DEFINE(AOsystFPWFSan);
static const MILK_RUNNABLE runnable_aosim_fpwfsan = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "AOsystFPWFSan", .description = "focal-plane WFS sensitivity" },
    .srcfile  = __FILE__,
    .farg     = farg_fpwfsan,
    .nbarg    = (int) (sizeof(farg_fpwfsan) / sizeof(farg_fpwfsan[0])),
    .cli_func = cli_aosim_fpwfsan,
    .rt       = &AOsystFPWFSan_rt,
};

static const MILK_RUNNABLE *const runnables[] = {
    &runnable_aosim_filt,
    &runnable_aosim_fitpup,
    &runnable_aosim_mkwf,
    &runnable_aosim_pyrwfs,
    &runnable_aosim_dm,
    &runnable_aosim_lowfs,
    &runnable_aosim_run,
    &runnable_aosim_exao,
    &runnable_aosim_mkprobes,
    &runnable_aosim_fpwfsan,
};

MILK_MODULE_RUNNABLES(coffeeAOsystSim, _module_deps, runnables);
