/**
 * @file    PIAAACMCsimul_exec_optimize_fpm_zones.c
 * @brief   Optimize focal plane mask zone values
 *
 *
 *
 */



#include <stdio.h>
#include <stdlib.h>



// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "statistic/statistic.h"

#include "OpticsMaterials/OpticsMaterials.h"
#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"

#include "PIAACMCsimul_init.h"
#include "PIAACMCsimul_initpiaacmcconf.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"

#include "PIAAshape/makePIAAshapes.h"



extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTSYST *optsyst;

extern OPTPIAACMCDESIGN *piaacmc;





/**
 * ---
 *
 *  ## Mode 13: Optimize focal plane mask zones only
 *
 * Uses "fast" mode:
 *
 * After mode 11, we can use the (complex) light propagated from each zone to compute the impact of
 * any thickness (sag) of that zone: the zone thickness induces a phase rotation for that zone,
 * which is applied to the unobstructed light from that zone as a complex rotation.
 *
 * The search is via steepest descent from random starting points.
 *
 * This mode only sets up the optimization that actually happens after exiting the switch statement if piaacmcsimul_var.LINOPT = 1 (as does mode 40)
 *
 */


errno_t exec_optimize_fpm_zones()
{
    DEBUG_TRACE_FSTART();

    imageID IDv;
    double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;

    char stopfile[STRINGMAXLEN_FULLFILENAME];

    imageID IDfpmresp;

    //imageID IDstatus;
    char fname[STRINGMAXLEN_FULLFILENAME];

    FILE *fp;


    int tmpd1;



    printf("=================================== mode 013 ===================================\n");


    // set the name of the stopfile
    WRITE_FULLFILENAME(stopfile, "%s/stoploop13.txt", piaacmcsimul_var.piaacmcconfdir);



    // get cli variables

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


    // FPM sag regularization control flag, do if == 1
    piaacmcsimul_var.linopt_REGFPMSAG = 1; // default
    if((IDv = variable_ID("REGFPMSAG")) != -1)
    {
        piaacmcsimul_var.linopt_REGFPMSAG = (long) data.variable[IDv].value.f + 0.01;
    }




    // set current state for statistical tracking
    //    data.image[IDstatus].array.UI16[0] = 0;

    // usual initialization
    FUNC_CHECK_RETURN(
        PIAACMCsimul_initpiaacmcconf(1, fpmradld, centobs0, centobs1, 0, 1)
    );

    FUNC_CHECK_RETURN(makePIAAshapes(piaacmc));

    FUNC_CHECK_RETURN(PIAACMCsimul_init(piaacmc, 0, 0.0, 0.0));

    // set current state for statistical tracking
    //data.image[IDstatus].array.UI16[0] = 1;

    // tracking diagnostic, giving the total flux in each plane
    WRITE_FULLFILENAME(fname, "%s/tmp_flux.txt", piaacmcsimul_var.piaacmcconfdir);
    fp = fopen(fname, "r");
    for(long elem = 0; elem < optsyst[0].NBelem; elem++)
    {
        //ret = fscanf(fp, "%lf %lf  %d\n", &tmplf1, &tmplf2, &tmpd1);

        int fscanfcnt;
        double tmplf1, tmplf2;
        fscanfcnt = fscanf(fp, "%lf %lf  %d\n", &tmplf1, &tmplf2, &tmpd1);
        if(fscanfcnt == EOF)
        {
            if(ferror(fp))
            {
                perror("fscanf");
            }
            else
            {
                fprintf(stderr,
                        "Error: fscanf reached end of file, no matching characters, no matching failure\n");
            }
            FUNC_RETURN_FAILURE("Call to fscanf failed");
        }
        else if(fscanfcnt != 3)
        {
            fprintf(stderr,
                    "Error: fscanf successfully matched and assigned %i input items, 3 expected\n",
                    fscanfcnt);
            FUNC_RETURN_FAILURE("Call to fscanf failed");
        }

        // scale flux to current number of lambda
        optsyst[0].flux[elem] = tmplf1 / tmpd1 * optsyst[0].nblambda;
    }
    fclose(fp);


    piaacmcsimul_var.LINOPT =
        1; // perform linear optimization after the switch exits
    // get the number of iterations
    /*
    if((IDv = variable_ID("PIAACMC_nbiter")) != -1)
    {
        NBiter = (long) data.variable[IDv].value.f + 0.01;
    }
    else
    {
        NBiter = 50;    // default number of iterations
    }
    */
    // set current state for statistical tracking
    //data.image[IDstatus].array.UI16[0] = 2;

    // get the FPMresp array computed in mode 11
    PIAACMCsimul_update_fnamedescr_conf();
    WRITE_FULLFILENAME(fname, "%s/FPMresp%d.%s.fits", piaacmcsimul_var.piaacmcconfdir,
                       piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.fnamedescr_conf);



    FUNC_CHECK_RETURN(
        load_fits(fname, "FPMresp", 1, &IDfpmresp)
    );

    piaacmcsimul_var.vsize =
        data.image[IDfpmresp].md[0].size[0]; // number of eval pts x2
    // make an array that holds the resulting light for evaluation point given the FPM solution, for each wavelenth
    //ID =
    FUNC_CHECK_RETURN(
        create_2Dimage_ID("imvect1", piaacmcsimul_var.vsize, piaacmc[0].nblambda, NULL)
    );

    // allocate arrays for fast routine
    // define convenient array variables
    piaacmcsimul_var.fpmresp_array = data.image[IDfpmresp].array.D;
    piaacmcsimul_var.zonez_array = data.image[piaacmc[0].zonezID].array.D;
    // allocate derivative of phase against thickness array
    piaacmcsimul_var.dphadz_array = (double *) malloc(sizeof(
                                        double) * piaacmc[0].nblambda);
    // compute this derivative
    for(long k = 0; k < piaacmc[0].nblambda; k++)
    {
        // OPTICSMATERIALS_pha_lambda computes change in phase per unit thickness at specified wavelength
        // second arg is thickness, so 1.0 meters determines result per meter thickness
        piaacmcsimul_var.dphadz_array[k] = OpticsMaterials_pha_lambda(
                                               piaacmc[0].fpmmaterial_code, 1.0, optsyst[0].lambdaarray[k]);
        printf("%ld  %g %g\n", k, optsyst[0].lambdaarray[k],
               piaacmcsimul_var.dphadz_array[k]);
    }
    piaacmcsimul_var.outtmp_array = (double *) malloc(sizeof(double) *
                                    (piaacmcsimul_var.vsize * piaacmc[0].nblambda +
                                     data.image[piaacmc[0].zonezID].md[0].size[0]));

    // do the fast optimization using the results of mode 11
    piaacmcsimul_var.computePSF_FAST_FPMresp = 1;

    // set current state for statistical tracking
    //data.image[IDstatus].array.UI16[0] = 3;

    // read the contrast normalization factor into CnormFactor
    WRITE_FULLFILENAME(fname, "%s/CnormFactor.txt", piaacmcsimul_var.piaacmcconfdir);
    fp = fopen(fname, "r");
    {
        int fscanfcnt;

        fscanfcnt = fscanf(fp, "%lf", &piaacmcsimul_var.CnormFactor);
        if(fscanfcnt == EOF)
        {
            if(ferror(fp))
            {
                perror("fscanf");
            }
            else
            {
                fprintf(stderr,
                        "Error: fscanf reached end of file, no matching characters, no matching failure\n");
            }
            FUNC_RETURN_FAILURE("Call to fscanf failed");
        }
        else if(fscanfcnt != 1)
        {
            fprintf(stderr,
                    "Error: fscanf successfully matched and assigned %i input items, 1 expected\n",
                    fscanfcnt);
            FUNC_RETURN_FAILURE("Call to fscanf failed");
        }
    }
    //ret = fscanf(fp, "%lf", &piaacmcsimul_var.CnormFactor);
    fclose(fp);
    // for each zone, add a random offset in range +- MODampl
    // this randomizes the starting point for each zone
    // data.image[piaacmc[0].zonezID].array.D[k] is set in PIAACMCsimul_run()
    for(long k = 0; k < data.image[piaacmc[0].zonezID].md[0].size[0]; k++)
    {
        data.image[piaacmc[0].zonezID].array.D[k] += piaacmcsimul_var.MODampl *
                (1.0 - 2.0 * ran1());
    }

    // set up optimization parameters for each zone
    // uses abstract specification of optimization parameters called paramval, ...
    piaacmcsimul_var.linopt_number_param = 0;
    for(long mz = 0; mz < data.image[piaacmc[0].zonezID].md[0].size[0]; mz++)
    {
        // parameter type
        piaacmcsimul_var.linopt_paramtype[piaacmcsimul_var.linopt_number_param] =
            _DATATYPE_DOUBLE;
        // value: sag of each zone
        piaacmcsimul_var.linopt_paramval[piaacmcsimul_var.linopt_number_param] =
            &data.image[piaacmc[0].zonezID].array.D[mz];
        // derivative step size
        piaacmcsimul_var.linopt_paramdelta[piaacmcsimul_var.linopt_number_param] =
            3.0e-9;
        // max parameter step size
        piaacmcsimul_var.linopt_parammaxstep[piaacmcsimul_var.linopt_number_param] =
            1.0e-6;
        // max and min allowed values for parameter
        piaacmcsimul_var.linopt_parammin[piaacmcsimul_var.linopt_number_param] =
            piaacmc[0].fpmminsag;
        piaacmcsimul_var.linopt_parammax[piaacmcsimul_var.linopt_number_param] =
            piaacmc[0].fpmmaxsag;
        // move on to next parameter
        piaacmcsimul_var.linopt_number_param++;
    }
    piaacmcsimul_var.PIAACMC_FPM_FASTDERIVATIVES =
        1; // for fast execution using analytic derivatives

    // set current state for statistical tracking
    //data.image[IDstatus].array.UI16[0] = 4;
    // Now on to the actual optimization, after exit from the switch statement
    // I hope you have a lot of time...

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}









