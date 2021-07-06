/**
 * @file    PIAAACMCsimul_optimize_PIAA_shapes.c
 * @brief   Optimize PIAA optics shapes
 *
 *
 */



// System includes
#include <stdio.h>
#include <stdlib.h>




// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_memory/COREMOD_memory.h"

#include "PIAACMCsimul.h"
#include "init_piaacmcopticaldesign.h"






/**
 * ---
 *
 * ## Mode 4: Optimize PIAA optics shapes, cosine modes only (not currently used, replaced by mode 40. skipping)
 *
 */
errno_t exec_optimize_PIAA_shapes()
{
    DEBUG_TRACE_FSTART();

    imageID IDv;
    double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;
    //long NBiter = 1000;
    long kmax;


    printf("=================================== mode 004 ===================================\n");
    // load some more cli variables
    if((IDv = variable_ID("PIAACMC_centobs0")) != -1)
    {
        centobs0 = data.variable[IDv].value.f;
    }
    if((IDv = variable_ID("PIAACMC_centobs1")) != -1)
    {
        centobs1 = data.variable[IDv].value.f;
    }
    if((IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }

    FUNC_CHECK_RETURN(init_piaacmcopticaldesign(0, fpmradld, centobs0, centobs1, 0, 1));

    piaacmcparams.LINOPT = 1; // perform linear optimization
    /*if((IDv = variable_ID("PIAACMC_nbiter")) != -1)
    {
        NBiter = (long) data.variable[IDv].value.f + 0.01;
    }
    else
    {
        NBiter = 1000;
    }*/

    kmax = data.image[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0];
    if((IDv = variable_ID("PIAACMC_maxoptCterm")) != -1)
    {
        kmax = (long) data.variable[IDv].value.f + 0.01;
    }

    if(kmax > data.image[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0])
    {
        kmax = data.image[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0];
    }

    piaacmcparams.linopt_number_param = 0;
    for(long k = 0; k < kmax; k++)
    {
        piaacmcparams.linopt_paramtype[piaacmcparams.linopt_number_param] =
            _DATATYPE_FLOAT;
        piaacmcparams.linopt_paramvalf[piaacmcparams.linopt_number_param] =
            &data.image[piaacmcopticaldesign.piaa0CmodesID].array.F[k];
        piaacmcparams.linopt_paramdelta[piaacmcparams.linopt_number_param] =
            1.0e-9;
        piaacmcparams.linopt_parammaxstep[piaacmcparams.linopt_number_param] =
            1.0e-8;
        piaacmcparams.linopt_parammin[piaacmcparams.linopt_number_param] =
            -1.0e-5;
        piaacmcparams.linopt_parammax[piaacmcparams.linopt_number_param] = 1.0e-5;
        piaacmcparams.linopt_number_param++;
    }

    for(long k = 0; k < kmax; k++)
    {
        piaacmcparams.linopt_paramtype[piaacmcparams.linopt_number_param] =
            _DATATYPE_FLOAT;
        piaacmcparams.linopt_paramvalf[piaacmcparams.linopt_number_param] =
            &data.image[piaacmcopticaldesign.piaa1CmodesID].array.F[k];
        piaacmcparams.linopt_paramdelta[piaacmcparams.linopt_number_param] =
            1.0e-9;
        piaacmcparams.linopt_parammaxstep[piaacmcparams.linopt_number_param] =
            1.0e-8;
        piaacmcparams.linopt_parammin[piaacmcparams.linopt_number_param] =
            -1.0e-5;
        piaacmcparams.linopt_parammax[piaacmcparams.linopt_number_param] = 1.0e-5;
        piaacmcparams.linopt_number_param++;
    }
    piaacmcparams.FORCE_MAKE_PIAA0shape = 1;
    piaacmcparams.FORCE_MAKE_PIAA1shape = 1;


    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

