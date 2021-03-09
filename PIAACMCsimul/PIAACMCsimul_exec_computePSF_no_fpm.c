/**
 * @file    PIAAACMCsimul_computePSF_no_fpm.c
 * @brief   Compute on-axis PSF with no focal plane mask
 * 
 *  
 * @author  O. Guyon
 * @date    25 nov 2017
 *
 * 
 * @bug No known bugs.
 * 
 */



// System includes
#include <stdio.h>
#include <stdlib.h>




// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_memory/COREMOD_memory.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"



extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTSYST *optsyst;

extern OPTPIAACMCDESIGN *piaacmc;






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

double PIAACMCsimul_exec_computePSF_no_fpm()
{
    long IDv;
    double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;
    double val;

    printf("=================================== mode 003 ===================================\n");

    // load some more cli variables
    if( (IDv = variable_ID("PIAACMC_centobs0")) != -1)
        centobs0 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_centobs1")) != -1)
        centobs1 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }


    // init as in mode 0
    PIAACMCsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 1);
    PIAACMCsimul_makePIAAshapes(piaacmc, 0);
    optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm

    piaacmcsimul_var.linopt_paramrefval[0] = piaacmc[0].fpmaskamptransm;

    piaacmc[0].fpmaskamptransm = -1.0;  // Remove focal plane mask
    piaacmcsimul_var.FORCE_CREATE_fpmza = 1;
    PIAACMCsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 0);
    // compute the PSF for an on-axis source, all optical elements
    val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);

    // restore original configuration
    piaacmc[0].fpmaskamptransm = piaacmcsimul_var.linopt_paramrefval[0];
    PIAACMCsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 0);
    PIAACMCsimul_savepiaacmcconf( piaacmcsimul_var.piaacmcconfdir);
    piaacmcsimul_var.FORCE_CREATE_fpmza = 0;

    return val;
}

