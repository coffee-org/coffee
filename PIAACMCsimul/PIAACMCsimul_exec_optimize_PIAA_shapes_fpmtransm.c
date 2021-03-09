/**
 * @file    PIAAACMCsimul_exec_optimize_PIAA_shapes_fpmtransm.c
 * @brief   Optimize PIAA shapes and focal plane mask transmission
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
#include "COREMOD_iofits/COREMOD_iofits.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"



extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTSYST *optsyst;

extern OPTPIAACMCDESIGN *piaacmc;






/**
 * ---
 *
 * ## Mode 40: Optimize PIAA optics shapes (and focal plane mask transmission for idealized PIAACMC)
 *
 *
 */

int PIAACMCsimul_exec_optimize_PIAA_shapes_fpmtransm()
{
    long IDv;
    double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;
    long NBiter = 1000;


    long mz, k;
    long ID_CPAfreq;
    long kmaxC, kmaxF;





    printf("=================================== mode 040 ===================================\n");
    //		piaacmcsimul_var.FORCE_CREATE_fpmza = 1;

    if( (IDv = variable_ID("PIAACMC_centobs0")) != -1)
        centobs0 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_centobs1")) != -1)
        centobs1 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }



    piaacmcsimul_var.PIAACMC_fpmtype = 0; // idealized (default)
    if((IDv=variable_ID("PIAACMC_fpmtype"))!=-1)
        piaacmcsimul_var.PIAACMC_fpmtype = (int) (data.variable[IDv].value.f + 0.1);
    printf("PIAACMC_fpmtype = %d\n", piaacmcsimul_var.PIAACMC_fpmtype);


    piaacmcsimul_var.FORCE_CREATE_fpmza = 1;
    PIAACMCsimul_initpiaacmcconf(piaacmcsimul_var.PIAACMC_fpmtype, fpmradld, centobs0, centobs1, 0, 1);
    //       printf("data.image[piaacmc[0].zoneaID].array.D[0] = %lf\n", data.image[piaacmc[0].zoneaID].array.D[0]);
    //        sleep(10);



    PIAACMCsimul_makePIAAshapes(piaacmc, 0);
    optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm


    if(0) // TEST
    {
        double valref;

        valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 1, 0, 0, 1);
        printf("valref = %g\n", valref);
        printf("EXEC CASE 40 COMPLETED\n");
    }
    else
    {
        piaacmcsimul_var.LINOPT = 1; // perform linear optimization

        if((IDv=variable_ID("PIAACMC_nbiter"))!=-1)
            piaacmcsimul_var.linopt_NBiter = (long) data.variable[IDv].value.f+0.01;
        else
            piaacmcsimul_var.linopt_NBiter = 1000;


        kmaxC = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];
        if((IDv=variable_ID("PIAACMC_maxoptCterm"))!=-1)
            kmaxC = (long) data.variable[IDv].value.f+0.01;
        if(kmaxC>data.image[piaacmc[0].piaa0CmodesID].md[0].size[0])
            kmaxC = data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];


        kmaxF = data.image[piaacmc[0].piaa0FmodesID].md[0].size[0];
        if((IDv=variable_ID("PIAACMC_maxoptFterm"))!=-1)
            kmaxF = (long) data.variable[IDv].value.f+0.01;
        if(kmaxF>data.image[piaacmc[0].piaa0FmodesID].md[0].size[0])
            kmaxF = data.image[piaacmc[0].piaa0FmodesID].md[0].size[0];




        // PIAA shapes regularization

        piaacmcsimul_var.linopt_REGPIAASHAPES = 0; // default
        if((IDv=variable_ID("REGPIAASHAPES"))!=-1)
            piaacmcsimul_var.linopt_REGPIAASHAPES = (long) data.variable[IDv].value.f+0.01;

        piaacmcsimul_var.linopt_piaa0C_regcoeff = 0.0e-7; // regularization coeff
        piaacmcsimul_var.linopt_piaa1C_regcoeff = 0.0e-7; // regularization coeff
        if((IDv=variable_ID("REGPIAA_C_COEFF"))!=-1)
        {
            piaacmcsimul_var.linopt_piaa0C_regcoeff = data.variable[IDv].value.f;
            piaacmcsimul_var.linopt_piaa1C_regcoeff = data.variable[IDv].value.f;
        }

        piaacmcsimul_var.linopt_piaa0C_regcoeff_alpha = 1.0; // regularization coeff power
        piaacmcsimul_var.linopt_piaa1C_regcoeff_alpha = 1.0; // regularization coeff power
        if((IDv=variable_ID("REGPIAA_C_ALPHA"))!=-1)
        {
            piaacmcsimul_var.linopt_piaa0C_regcoeff_alpha = data.variable[IDv].value.f;
            piaacmcsimul_var.linopt_piaa1C_regcoeff_alpha = data.variable[IDv].value.f;
        }

        piaacmcsimul_var.linopt_piaa0F_regcoeff = 0.0e-7; // regularization coeff
        piaacmcsimul_var.linopt_piaa1F_regcoeff = 0.0e-7; // regularization coeff
        if((IDv=variable_ID("REGPIAA_F_COEFF"))!=-1)
        {
            piaacmcsimul_var.linopt_piaa0F_regcoeff = data.variable[IDv].value.f;
            piaacmcsimul_var.linopt_piaa1F_regcoeff = data.variable[IDv].value.f;
        }

        piaacmcsimul_var.linopt_piaa0F_regcoeff_alpha = 1.0; // regularization coeff power
        piaacmcsimul_var.linopt_piaa1F_regcoeff_alpha = 1.0; // regularization coeff power
        if((IDv=variable_ID("REGPIAA_F_ALPHA"))!=-1)
        {
            piaacmcsimul_var.linopt_piaa0F_regcoeff_alpha = data.variable[IDv].value.f;
            piaacmcsimul_var.linopt_piaa1F_regcoeff_alpha = data.variable[IDv].value.f;
        }

        if(piaacmcsimul_var.linopt_REGPIAASHAPES==1)
        {
            printf("loading CPA modes frequ\n");
            fflush(stdout);
            ID_CPAfreq = image_ID("cpamodesfreq");
            if(ID_CPAfreq == -1)
                ID_CPAfreq = load_fits("cpamodesfreq.fits", "cpamodesfreq", 2);
        }


        // FPM SAG regularization

        piaacmcsimul_var.linopt_REGFPMSAG = 1; // default
        if((IDv=variable_ID("REGFPMSAG"))!=-1)
            piaacmcsimul_var.linopt_REGFPMSAG = (long) data.variable[IDv].value.f+0.01;




        piaacmcsimul_var.linopt_number_param = 0;

        if(piaacmcsimul_var.PIAACMC_fpmtype==0) // ideal mask
        {
            piaacmcsimul_var.linopt_paramtype[piaacmcsimul_var.linopt_number_param] = _DATATYPE_DOUBLE;
            piaacmcsimul_var.linopt_paramval[piaacmcsimul_var.linopt_number_param] = &data.image[piaacmc[0].zoneaID].array.D[0];
            piaacmcsimul_var.linopt_paramdelta[piaacmcsimul_var.linopt_number_param] = 1.0e-3;
            piaacmcsimul_var.linopt_parammaxstep[piaacmcsimul_var.linopt_number_param] = 2.0e-1;
            piaacmcsimul_var.linopt_parammin[piaacmcsimul_var.linopt_number_param] = -1.0;
            piaacmcsimul_var.linopt_parammax[piaacmcsimul_var.linopt_number_param] = 1.0;
            piaacmcsimul_var.linopt_number_param++;
        }
        else // real physical mask
        {
            if(variable_ID("PIAACMC_mzOPT")!=-1) // optimize zones
            {
                for(mz=0; mz<data.image[piaacmc[0].zonezID].md[0].size[0]; mz++)
                {
                    piaacmcsimul_var.linopt_paramtype[piaacmcsimul_var.linopt_number_param] = _DATATYPE_DOUBLE;
                    piaacmcsimul_var.linopt_paramval[piaacmcsimul_var.linopt_number_param] = &data.image[piaacmc[0].zonezID].array.D[mz];
                    piaacmcsimul_var.linopt_paramdelta[piaacmcsimul_var.linopt_number_param] = 1.0e-9;
                    piaacmcsimul_var.linopt_parammaxstep[piaacmcsimul_var.linopt_number_param] = 5.0e-8;
                    piaacmcsimul_var.linopt_parammin[piaacmcsimul_var.linopt_number_param] = -2.0e-6;
                    piaacmcsimul_var.linopt_parammax[piaacmcsimul_var.linopt_number_param] = 2.0e-6;
                    piaacmcsimul_var.linopt_number_param++;
                }
            }
        }


        for(k=0; k<kmaxC; k++)
        {
            piaacmcsimul_var.linopt_paramtype[piaacmcsimul_var.linopt_number_param] = _DATATYPE_FLOAT;
            piaacmcsimul_var.linopt_paramvalf[piaacmcsimul_var.linopt_number_param] = &data.image[piaacmc[0].piaa0CmodesID].array.F[k];
            piaacmcsimul_var.linopt_paramdelta[piaacmcsimul_var.linopt_number_param] = 1.0e-10;
            piaacmcsimul_var.linopt_parammaxstep[piaacmcsimul_var.linopt_number_param] = 1.0e-7;
            piaacmcsimul_var.linopt_parammin[piaacmcsimul_var.linopt_number_param] = -1.0e-3;
            piaacmcsimul_var.linopt_parammax[piaacmcsimul_var.linopt_number_param] = 1.0e-3;
            piaacmcsimul_var.linopt_number_param++;
        }

        for(k=0; k<kmaxC; k++)
        {
            piaacmcsimul_var.linopt_paramtype[piaacmcsimul_var.linopt_number_param] = _DATATYPE_FLOAT;
            piaacmcsimul_var.linopt_paramvalf[piaacmcsimul_var.linopt_number_param] = &data.image[piaacmc[0].piaa1CmodesID].array.F[k];
            piaacmcsimul_var.linopt_paramdelta[piaacmcsimul_var.linopt_number_param] = 1.0e-10;
            piaacmcsimul_var.linopt_parammaxstep[piaacmcsimul_var.linopt_number_param] = 1.0e-7;
            piaacmcsimul_var.linopt_parammin[piaacmcsimul_var.linopt_number_param] = -1.0e-3;
            piaacmcsimul_var.linopt_parammax[piaacmcsimul_var.linopt_number_param] = 1.0e-3;
            piaacmcsimul_var.linopt_number_param++;
        }

        for(k=0; k<kmaxF; k++)
        {
            piaacmcsimul_var.linopt_paramtype[piaacmcsimul_var.linopt_number_param] = _DATATYPE_FLOAT;
            piaacmcsimul_var.linopt_paramvalf[piaacmcsimul_var.linopt_number_param] = &data.image[piaacmc[0].piaa0FmodesID].array.F[k];
            piaacmcsimul_var.linopt_paramdelta[piaacmcsimul_var.linopt_number_param] = 1.0e-10;
            piaacmcsimul_var.linopt_parammaxstep[piaacmcsimul_var.linopt_number_param] = 1.0e-7;
            piaacmcsimul_var.linopt_parammin[piaacmcsimul_var.linopt_number_param] = -1.0e-3;
            piaacmcsimul_var.linopt_parammax[piaacmcsimul_var.linopt_number_param] = 1.0e-3;
            piaacmcsimul_var.linopt_number_param++;
        }

        for(k=0; k<kmaxF; k++)
        {
            piaacmcsimul_var.linopt_paramtype[piaacmcsimul_var.linopt_number_param] = _DATATYPE_FLOAT;
            piaacmcsimul_var.linopt_paramvalf[piaacmcsimul_var.linopt_number_param] = &data.image[piaacmc[0].piaa1FmodesID].array.F[k];
            piaacmcsimul_var.linopt_paramdelta[piaacmcsimul_var.linopt_number_param] = 1.0e-10;
            piaacmcsimul_var.linopt_parammaxstep[piaacmcsimul_var.linopt_number_param] = 1.0e-7;
            piaacmcsimul_var.linopt_parammin[piaacmcsimul_var.linopt_number_param] = -1.0e-3;
            piaacmcsimul_var.linopt_parammax[piaacmcsimul_var.linopt_number_param] = 1.0e-3;
            piaacmcsimul_var.linopt_number_param++;
        }

        piaacmcsimul_var.FORCE_MAKE_PIAA0shape = 1;
        piaacmcsimul_var.FORCE_MAKE_PIAA1shape = 1;
    }

    return 0;
}




