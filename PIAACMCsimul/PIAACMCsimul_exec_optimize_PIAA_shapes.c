/**
 * @file    PIAAACMCsimul_optimize_PIAA_shapes.c
 * @brief   Optimize PIAA optics shapes
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

#include "PIAACMCsimul/PIAACMCsimul.h"


 
 
extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTPIAACMCDESIGN *piaacmc;





/**
 * ---
 *
 * ## Mode 4: Optimize PIAA optics shapes, cosine modes only (not currently used, replaced by mode 40. skipping)
 *
 */
int PIAACMCsimul_exec_optimize_PIAA_shapes()
{
    long IDv;
    double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;
    long NBiter = 1000;
    long kmax;
    long k;


    printf("=================================== mode 004 ===================================\n");
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

    PIAACMCsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 1);
    piaacmcsimul_var.LINOPT = 1; // perform linear optimization
    if((IDv=variable_ID("PIAACMC_nbiter"))!=-1)
        NBiter = (long) data.variable[IDv].value.f+0.01;
    else
        NBiter = 1000;

    kmax = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];
    if((IDv=variable_ID("PIAACMC_maxoptCterm"))!=-1)
        kmax = (long) data.variable[IDv].value.f+0.01;

    if(kmax>data.image[piaacmc[0].piaa0CmodesID].md[0].size[0])
        kmax = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];

    piaacmcsimul_var.linopt_number_param = 0;
    for(k=0; k<kmax; k++)
    {
        piaacmcsimul_var.linopt_paramtype[piaacmcsimul_var.linopt_number_param] = _DATATYPE_FLOAT;
        piaacmcsimul_var.linopt_paramvalf[piaacmcsimul_var.linopt_number_param] = &data.image[piaacmc[0].piaa0CmodesID].array.F[k];
        piaacmcsimul_var.linopt_paramdelta[piaacmcsimul_var.linopt_number_param] = 1.0e-9;
        piaacmcsimul_var.linopt_parammaxstep[piaacmcsimul_var.linopt_number_param] = 1.0e-8;
        piaacmcsimul_var.linopt_parammin[piaacmcsimul_var.linopt_number_param] = -1.0e-5;
        piaacmcsimul_var.linopt_parammax[piaacmcsimul_var.linopt_number_param] = 1.0e-5;
        piaacmcsimul_var.linopt_number_param++;
    }

    for(k=0; k<kmax; k++)
    {
        piaacmcsimul_var.linopt_paramtype[piaacmcsimul_var.linopt_number_param] = _DATATYPE_FLOAT;
        piaacmcsimul_var.linopt_paramvalf[piaacmcsimul_var.linopt_number_param] = &data.image[piaacmc[0].piaa1CmodesID].array.F[k];
        piaacmcsimul_var.linopt_paramdelta[piaacmcsimul_var.linopt_number_param] = 1.0e-9;
        piaacmcsimul_var.linopt_parammaxstep[piaacmcsimul_var.linopt_number_param] = 1.0e-8;
        piaacmcsimul_var.linopt_parammin[piaacmcsimul_var.linopt_number_param] = -1.0e-5;
        piaacmcsimul_var.linopt_parammax[piaacmcsimul_var.linopt_number_param] = 1.0e-5;
        piaacmcsimul_var.linopt_number_param++;
    }
    piaacmcsimul_var.FORCE_MAKE_PIAA0shape = 1;
    piaacmcsimul_var.FORCE_MAKE_PIAA1shape = 1;


    return 0;
}

