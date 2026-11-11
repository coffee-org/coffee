/**
 * @file    coronagraphs_cli.c
 * @brief   MILK_RUNNABLE command wrappers and module registry for coffeecoronagraphs
 */

#define MODULE_SHORTNAME_DEFAULT "coro"
#define MODULE_DESCRIPTION       "Coronagraph simulation tools"

#include "milk_runnable.h"
#include "coronagraphs.h"

MODULE_DEPS("milkimagegen", "milkfft", "milkimagebasic",
            "milkinfo", "milkstatistic", "milkZernikePolyn",
            "milkOpticsMaterials", "milkWFpropagate");

/* 1. cormk2Dprolate */
static double *mk2dp_masksizepix;
static double *mk2dp_beamradpix;
static double *mk2dp_centobs;
static char   *mk2dp_outname;
static long   *mk2dp_size;

static CLICMDARGDEF farg_mk2dprolate[] = {
    { CLIARG_FLOAT64, ".masksizepix", "focal plane mask radius [pix]", "10.0",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &mk2dp_masksizepix, NULL },
    { CLIARG_FLOAT64, ".beamradpix", "pupil beam radius [pix]", "200.0",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &mk2dp_beamradpix, NULL },
    { CLIARG_FLOAT64, ".centobs", "central obstruction ratio", "0.2",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &mk2dp_centobs, NULL },
    { CLIARG_STR_NOT_IMG, ".outname", "output prolate image name", "prolim",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &mk2dp_outname, NULL },
    { CLIARG_INT64, ".size", "array size [pix]", "1024",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &mk2dp_size, NULL },
};

/**
 * cli_coronagraph_make_2Dprolate - Invoke 2D prolate apodizer computation (pixel mask radius)
 *
 * Return: errno_t status code from coronagraph_make_2Dprolate.
 */
static errno_t cli_coronagraph_make_2Dprolate(void)
{
    return coronagraph_make_2Dprolate(
        *mk2dp_masksizepix,
        *mk2dp_beamradpix,
        *mk2dp_centobs,
        mk2dp_outname,
        *mk2dp_size,
        "NULLim",
        NULL);
}

MILK_RUNNABLE_RT_DEFINE(cormk2Dprolate);
static const MILK_RUNNABLE runnable_coro_mk2Dprolate = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "cormk2Dprolate", .description = "make 2D prolate apodizer" },
    .srcfile  = __FILE__,
    .farg     = farg_mk2dprolate,
    .nbarg    = (int) (sizeof(farg_mk2dprolate) / sizeof(farg_mk2dprolate[0])),
    .cli_func = cli_coronagraph_make_2Dprolate,
    .rt       = &cormk2Dprolate_rt,
};

/* 2. cormk2Dprolateld */
static double *mk2dpld_masksizeld;
static double *mk2dpld_beamradpix;
static double *mk2dpld_centobs;
static char   *mk2dpld_outname;
static long   *mk2dpld_size;

static CLICMDARGDEF farg_mk2dprolateld[] = {
    { CLIARG_FLOAT64, ".masksizeld", "focal plane mask radius [l/D]", "1.4",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &mk2dpld_masksizeld, NULL },
    { CLIARG_FLOAT64, ".beamradpix", "pupil beam radius [pix]", "200.0",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &mk2dpld_beamradpix, NULL },
    { CLIARG_FLOAT64, ".centobs", "central obstruction ratio", "0.2",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &mk2dpld_centobs, NULL },
    { CLIARG_STR_NOT_IMG, ".outname", "output prolate image name", "prolim",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &mk2dpld_outname, NULL },
    { CLIARG_INT64, ".size", "array size [pix]", "1024",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &mk2dpld_size, NULL },
};

/**
 * cli_coronagraph_make_2Dprolateld - Invoke 2D prolate apodizer computation (l/D mask radius)
 *
 * Return: errno_t status code from coronagraph_make_2Dprolateld.
 */
static errno_t cli_coronagraph_make_2Dprolateld(void)
{
    return coronagraph_make_2Dprolateld(
        *mk2dpld_masksizeld,
        *mk2dpld_beamradpix,
        *mk2dpld_centobs,
        mk2dpld_outname,
        *mk2dpld_size,
        "NULL");
}

MILK_RUNNABLE_RT_DEFINE(cormk2Dprolateld);
static const MILK_RUNNABLE runnable_coro_mk2Dprolateld = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "cormk2Dprolateld", .description = "make 2D prolate (l/D mask size)" },
    .srcfile  = __FILE__,
    .farg     = farg_mk2dprolateld,
    .nbarg    = (int) (sizeof(farg_mk2dprolateld) / sizeof(farg_mk2dprolateld[0])),
    .cli_func = cli_coronagraph_make_2Dprolateld,
    .rt       = &cormk2Dprolateld_rt,
};

/* 3. corup2Dprolate */
static double *up2dp_masksizeld;
static double *up2dp_beamradpix;
static double *up2dp_centobs;
static double *up2dp_zfactor;

static CLICMDARGDEF farg_up2dprolate[] = {
    { CLIARG_FLOAT64, ".masksizeld", "focal plane mask radius [l/D]", "1.4",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &up2dp_masksizeld, NULL },
    { CLIARG_FLOAT64, ".beamradpix", "pupil beam radius [pix]", "200.0",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &up2dp_beamradpix, NULL },
    { CLIARG_FLOAT64, ".centobs", "central obstruction ratio", "0.2",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &up2dp_centobs, NULL },
    { CLIARG_FLOAT64, ".zfactor", "DFT zoom factor", "4.0",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &up2dp_zfactor, NULL },
};

/**
 * cli_coronagraph_update_2Dprolate - Update existing 2D prolate solution iteratively
 *
 * Return: errno_t status code from coronagraph_update_2Dprolate.
 */
static errno_t cli_coronagraph_update_2Dprolate(void)
{
    return coronagraph_update_2Dprolate(
        *up2dp_masksizeld,
        *up2dp_beamradpix,
        *up2dp_centobs,
        *up2dp_zfactor);
}

MILK_RUNNABLE_RT_DEFINE(corup2Dprolate);
static const MILK_RUNNABLE runnable_coro_up2Dprolate = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "corup2Dprolate", .description = "update 2D prolate apodizer" },
    .srcfile  = __FILE__,
    .farg     = farg_up2dprolate,
    .nbarg    = (int) (sizeof(farg_up2dprolate) / sizeof(farg_up2dprolate[0])),
    .cli_func = cli_coronagraph_update_2Dprolate,
    .rt       = &corup2Dprolate_rt,
};

/* 4. corprolcomp */
MILK_RUNNABLE_RT_DEFINE(corprolcomp);
static const MILK_RUNNABLE runnable_coro_prolcomp = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "corprolcomp", .description = "compile APLC prolate functions data" },
    .srcfile  = __FILE__,
    .cli_func = coronagraph_APLCapo_compile,
    .rt       = &corprolcomp_rt,
};

/* 5. corsimpsf */
static double *simpsf_xld;
static double *simpsf_yld;
static char   *simpsf_name;
static long   *simpsf_type;

static CLICMDARGDEF farg_simpsf[] = {
    { CLIARG_FLOAT64, ".xld", "source X offset [l/D]", "0.0",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &simpsf_xld, NULL },
    { CLIARG_FLOAT64, ".yld", "source Y offset [l/D]", "0.0",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &simpsf_yld, NULL },
    { CLIARG_STR_NOT_IMG, ".psfname", "output PSF image name", "corpsf",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &simpsf_name, NULL },
    { CLIARG_INT64, ".corotype", "coronagraph type index", "50",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &simpsf_type, NULL },
};

/**
 * cli_coronagraph_simulPSF - Simulate off-axis or on-axis coronagraph PSF
 *
 * Return: errno_t status code from coronagraph_simulPSF.
 */
static errno_t cli_coronagraph_simulPSF(void)
{
    return coronagraph_simulPSF(
        *simpsf_xld,
        *simpsf_yld,
        simpsf_name,
        *simpsf_type,
        "");
}

MILK_RUNNABLE_RT_DEFINE(corsimpsf);
static const MILK_RUNNABLE runnable_coro_simpsf = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "corsimpsf", .description = "simulate coronagraph PSF" },
    .srcfile  = __FILE__,
    .farg     = farg_simpsf,
    .nbarg    = (int) (sizeof(farg_simpsf) / sizeof(farg_simpsf[0])),
    .cli_func = cli_coronagraph_simulPSF,
    .rt       = &corsimpsf_rt,
};

/* 6. coroscanpiaacmcperf */
static double *scan_obs0;

static CLICMDARGDEF farg_scanpiaacmc[] = {
    { CLIARG_FLOAT64, ".obs0", "input central obstruction ratio", "0.2",
      (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT), (void **) &scan_obs0, NULL },
};

/**
 * cli_scanPIAACMC_centObs_perf - Scan idealized PIAACMC performance for central obstruction
 *
 * Return: errno_t status code from CORONAGRAPHS_scanPIAACMC_centObs_perf.
 */
static errno_t cli_scanPIAACMC_centObs_perf(void)
{
    return CORONAGRAPHS_scanPIAACMC_centObs_perf(*scan_obs0);
}

MILK_RUNNABLE_RT_DEFINE(coroscanpiaacmcperf);
static const MILK_RUNNABLE runnable_coro_scanpiaacmcperf = {
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "coroscanpiaacmcperf", .description = "scan PIAACMC performance" },
    .srcfile  = __FILE__,
    .farg     = farg_scanpiaacmc,
    .nbarg    = (int) (sizeof(farg_scanpiaacmc) / sizeof(farg_scanpiaacmc[0])),
    .cli_func = cli_scanPIAACMC_centObs_perf,
    .rt       = &coroscanpiaacmcperf_rt,
};

static const MILK_RUNNABLE *const runnables[] = {
    &runnable_coro_mk2Dprolate,
    &runnable_coro_mk2Dprolateld,
    &runnable_coro_up2Dprolate,
    &runnable_coro_prolcomp,
    &runnable_coro_simpsf,
    &runnable_coro_scanpiaacmcperf,
};

MILK_MODULE_RUNNABLES(coffeecoronagraphs, _module_deps, runnables);
