/**
 * @file    PIAAACMCsimul.c
 * @brief   PIAA-type coronagraph design
 *
 * Can design both APLCMC and PIAACMC coronagraphs
 */

#define _GNU_SOURCE

#define MODULE_SHORTNAME_DEFAULT "coffeepiaacmcsim"
#define MODULE_DESCRIPTION       "PIAACMC simulation"

#include "CommandLineInterface/CLIcore.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul.h"

#include "PIAACMCsimul_run.h"

#include "FocalPlaneMask/FPM_process.h"
#include "FocalPlaneMask/FPMresp_resample.h"
#include "FocalPlaneMask/FPMresp_rmzones.h"
#include "FocalPlaneMask/rings2sectors.h"

#include "LyotStop/geomProp.h"

#ifdef HAVE_LIBGOMP
#include <omp.h>
#define OMP_NELEMENT_LIMIT 1000000
#endif

// static int INITSTATUS_PIAACMCsimul = 0;

/// optical system description
//OPTSYST *optsyst;

//OPTPIAACMCDESIGN *piaacmcoptdesign;

PIAACMCSIMUL_PARAMS piaacmcparams;
OPTPIAACMCDESIGN    piaacmcopticaldesign;
OPTSYST             piaacmcopticalsystem;

INIT_MODULE_LIB(coffee_PIAACMCsimul)

/** @name MODULE INITIALIZATION
 * Registers CLI commands
*/

static errno_t init_module_CLI()
{

    CLIADDCMD_PIAACMCsimul__ring2sectors();

    CLIADDCMD_PIAACMCsimul__run();

    CLIADDCMD_PIAACMCsimul__FPMresp_rmzones();

    CLIADDCMD_PIAACMCsimul__FPMresp_resample();

    CLIADDCMD_PIAACMCsimul__FPM_process();

    CLIADDCMD_PIAACMCsimul__LyotStop__PIAACMCsimul_geomProp();

    piaacmcparams.optsystinit = 0;

    // this makes 20% bandwidth, from 0.55/1.1 to 0.55*1.1
    piaacmcparams.LAMBDASTART = 0.5e-6;
    piaacmcparams.LAMBDAEND   = 0.605e-6;

    piaacmcparams.FORCE_CREATE_Cmodes = 0;
    piaacmcparams.CREATE_Cmodes       = 0;
    piaacmcparams.FORCE_CREATE_Fmodes = 0;
    piaacmcparams.CREATE_Fmodes       = 0;

    piaacmcparams.FORCE_CREATE_fpmzmap = 0;
    piaacmcparams.CREATE_fpmzmap       = 0;
    piaacmcparams.FORCE_CREATE_fpmzt   = 0;
    piaacmcparams.CREATE_fpmzt         = 0;

    piaacmcparams.FORCE_CREATE_fpmza = 0;
    piaacmcparams.CREATE_fpmza       = 0;

    piaacmcparams.FORCE_MAKE_PIAA0shape = 0;
    piaacmcparams.MAKE_PIAA0shape       = 0;
    piaacmcparams.FORCE_MAKE_PIAA1shape = 0;
    piaacmcparams.MAKE_PIAA1shape       = 0;

    // if != -1, compute only impulse response to corresponding zone
    piaacmcparams.focmMode           = -1;
    piaacmcparams.PIAACMC_FPMsectors = 0;

    piaacmcparams.FPMSCALEFACTOR    = 0.9;
    piaacmcparams.PIAACMC_MASKRADLD = 0.0;

    piaacmcparams.LOOPCNT = 0;

    piaacmcparams.CnormFactor = 1.0;

    piaacmcparams.computePSF_FAST_FPMresp = 0;
    piaacmcparams.computePSF_ResolvedTarget =
        0; // source size = 1e-{0.1*computePSF_ResolvedTarget}
    piaacmcparams.computePSF_ResolvedTarget_mode =
        0; // 0: source is simulated as 3 points, 1: source is simulated as 6 points
    piaacmcparams.PIAACMC_FPM_FASTDERIVATIVES = 0;

    piaacmcparams.SCORINGTOTAL    = 1.0;
    piaacmcparams.MODampl         = 1.0e-6;
    piaacmcparams.SCORINGMASKTYPE = 0;
    piaacmcparams.PIAACMC_save    = 1;
    //piaacmcparams.PIAACMC_MASKregcoeff = 1.0;
    piaacmcparams.PIAACMC_fpmtype = 0;

    piaacmcparams.WRITE_OK = 1;

    piaacmcparams.LINOPT = 0;

    // add atexit functions here
    atexit(PIAACMCsimul_free);

    return RETURN_SUCCESS;
}
