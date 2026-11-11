/**
 * @file    PIAACMCsimul.c
 * @brief   PIAA-type coronagraph design module descriptor and global state
 *
 * Supports design and Fresnel/polychromatic simulation of both APLCMC and
 * PIAACMC coronagraphs.
 */

#define _GNU_SOURCE

#define MODULE_SHORTNAME_DEFAULT "coffeepiaacmcsim"
#define MODULE_DESCRIPTION       "PIAACMC simulation"

#include <stdlib.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul.h"
#include "PIAACMCsimul_run.h"

#include "FocalPlaneMask/FPM_process.h"
#include "FocalPlaneMask/FPMresp_resample.h"
#include "FocalPlaneMask/FPMresp_rmzones.h"
#include "FocalPlaneMask/rings2sectors.h"

#include "LyotStop/geomProp.h"

PIAACMCSIMUL_PARAMS piaacmcparams =
{
    .optsystinit       = 0,
    .LAMBDASTART       = 0.5e-6,
    .LAMBDAEND         = 0.605e-6,
    .focmMode          = -1,
    .FPMSCALEFACTOR    = 0.9,
    .CnormFactor       = 1.0,
    .SCORINGTOTAL      = 1.0,
    .MODampl           = 1.0e-6,
    .PIAACMC_save      = 1,
    .WRITE_OK          = 1,
};

OPTPIAACMCDESIGN piaacmcopticaldesign;
OPTSYST          piaacmcopticalsystem;

static int piaacmc_atexit_registered = 0;

MODULE_DEPS(
    "milkinfo",
    "milklinoptimtools",
    "milkimagebasic",
    "milkimagefilter",
    "milkimagegen",
    "milkfft",
    "milkstatistic",
    "milkZernikePolyn",
    "milkOpticsMaterials",
    "milkWFpropagate",
    "coffeeOptSystProp",
    "coffeecoronagraphs");

/**
 * init_PIAACMCsimul - Initialize PIAACMCsimul module defaults and exit handler
 *
 * Return: RETURN_SUCCESS.
 */
errno_t init_PIAACMCsimul(void)
{
    if (!piaacmc_atexit_registered)
    {
        atexit(PIAACMCsimul_free);
        piaacmc_atexit_registered = 1;
    }
    return RETURN_SUCCESS;
}

/**
 * libinit_PIAACMCsimul - Shared library constructor registering cleanup hook
 */
void __attribute__((constructor)) libinit_PIAACMCsimul(void)
{
    init_PIAACMCsimul();
}

static const MILK_RUNNABLE *const runnables[] =
{
    &runnable_piaacmc_ring2sect,
    &runnable_piaacmc_run,
    &runnable_piaacmc_fpmresprm,
    &runnable_piaacmc_fpmresprs,
    &runnable_piaacmc_fpmprocess,
    &runnable_piaacmc_geomprop,
};

MILK_MODULE_RUNNABLES(coffeePIAACMCsimul, _module_deps, runnables);
