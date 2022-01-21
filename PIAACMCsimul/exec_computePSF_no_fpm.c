/**
 * @file    PIAAACMCsimul_computePSF_no_fpm.c
 * @brief   Compute on-axis PSF with no focal plane mask
 *
 *
 */

// System includes
#include <stdio.h>
#include <stdlib.h>

// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_memory/COREMOD_memory.h"

#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"

#include "PIAACMCsimul_computePSF.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"
#include "init_piaacmcopticaldesign.h"

#include "PIAAshape/makePIAAshapes.h"

/**
 * ---
 *
 * ## Mode 3: Calibrate, no focal plane mask (not currently used)
 *
 * Compute PSF and contrast with no focal plane mask with the current design.
*
* Provides the denominator for the contrast estimate
*
* Saved by PIAACMCsimul_computePSF as fits file "psfi0"

*/

errno_t exec_computePSF_no_fpm(double *outval)
{
    DEBUG_TRACE_FSTART();

    double fpmradld = 0.95; // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;
    double val;

    printf(
        "=================================== mode 003 "
        "===================================\n");

    // load some more cli variables
    {
        variableID IDv;

        if ((IDv = variable_ID("PIAACMC_centobs0")) != -1)
            centobs0 = data.variable[IDv].value.f;
        if ((IDv = variable_ID("PIAACMC_centobs1")) != -1)
            centobs1 = data.variable[IDv].value.f;
        if ((IDv = variable_ID("PIAACMC_fpmradld")) != -1)
        {
            fpmradld = data.variable[IDv].value.f;
            printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
        }
    }

    // init as in mode 0
    {
        uint64_t initflag = INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__READCONF;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF;

        FUNC_CHECK_RETURN(init_piaacmcopticaldesign(fpmradld,
                                                    centobs0,
                                                    centobs1,
                                                    initflag,
                                                    NULL));
    }

    FUNC_CHECK_RETURN(makePIAAshapes());

    piaacmcopticalsystem.FOCMASKarray[0].mode = 1; // use 1-fpm

    piaacmcparams.linopt_paramrefval[0] = piaacmcopticaldesign.fpmaskamptransm;

    piaacmcopticaldesign.fpmaskamptransm = -1.0; // Remove focal plane mask
    piaacmcparams.FORCE_CREATE_fpmza     = 1;

    FUNC_CHECK_RETURN(
        init_piaacmcopticaldesign(fpmradld,
                                  centobs0,
                                  centobs1,
                                  INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT,
                                  NULL));

    // compute the PSF for an on-axis source, all optical elements
    FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(0.0,
                                              0.0,
                                              0,
                                              piaacmcopticalsystem.NBelem,
                                              0,
                                              0,
                                              0,
                                              0,
                                              &val));

    // restore original configuration
    piaacmcopticaldesign.fpmaskamptransm = piaacmcparams.linopt_paramrefval[0];

    {
        uint64_t initflag = INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT;
        FUNC_CHECK_RETURN(init_piaacmcopticaldesign(fpmradld,
                                                    centobs0,
                                                    centobs1,
                                                    initflag,
                                                    NULL));
    }

    FUNC_CHECK_RETURN(
        PIAACMCsimul_savepiaacmcconf(piaacmcparams.piaacmcconfdir));

    piaacmcparams.FORCE_CREATE_fpmza = 0;

    if (outval == NULL)
    {
        DEBUG_TRACEPOINT("FOUT %f -> NULL", val);
    }
    else
    {
        *outval = val;
        DEBUG_TRACEPOINT("FOUT %f -> return ptr", val);
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
