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
#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "statistic/statistic.h"

#include "OptSystProp/OptSystProp.h"
#include "OpticsMaterials/OpticsMaterials.h"

#include "PIAACMCsimul.h"

#include "PIAACMCsimul_loadsavepiaacmcconf.h"
#include "init_piaacmcopticaldesign.h"
#include "init_piaacmcopticalsystem.h"

#include "PIAAshape/makePIAAshapes.h"

/**
 * read_cnorm_factor - Read contrast normalization factor from CnormFactor.txt
 * @fname: Path to CnormFactor.txt.
 * @val:   Output pointer receiving the parsed normalization factor.
 *
 * Skips leading blank lines and '#' header comments written by
 * PIAACMCsimul_computePSF().
 *
 * Return: RETURN_SUCCESS on valid parse, RETURN_FAILURE otherwise.
 */
static errno_t read_cnorm_factor(
    const char *restrict fname,
    double     *restrict val)
{
    FILE *fp = fopen(fname, "r");
    if (fp == NULL)
    {
        FUNC_RETURN_FAILURE("Cannot open %s", fname);
    }

    char line[256];
    int  parsed = 0;
    while (fgets(line, (int) sizeof(line), fp) != NULL)
    {
        const char *p = line;
        while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n')
        {
            p++;
        }
        if (*p == '\0' || *p == '#')
        {
            continue;
        }
        if (sscanf(p, "%lf", val) == 1)
        {
            parsed = 1;
        }
        break;
    }

    fclose(fp);
    if (parsed != 1)
    {
        FUNC_RETURN_FAILURE("Failed to parse CnormFactor from %s", fname);
    }

    return RETURN_SUCCESS;
}

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
 * This mode only sets up the optimization that actually happens after exiting the switch statement if piaacmcparams.LINOPT = 1 (as does mode 40)
 *
 */

errno_t exec_optimize_fpm_zones()
{
    DEBUG_TRACE_FSTART();

    imageID IDv;
    double  fpmradld = 0.95; // default
    double  centobs0 = 0.3;
    double  centobs1 = 0.2;

    char stopfile[STRINGMAXLEN_FULLFILENAME];

    imageID IDfpmresp;

    //imageID IDstatus;
    char fname[STRINGMAXLEN_FULLFILENAME];

    FILE *fp;

    int tmpd1;

    printf(
        "=================================== mode 013 "
        "===================================\n");

    // set the name of the stopfile
    WRITE_FULLFILENAME(stopfile,
                       "%s/stoploop13.txt",
                       piaacmcparams.piaacmcconfdir);

    // get cli variables

    if((IDv = variable_ID("PIAACMC_centobs0")) != -1)
    {
        centobs0 = coffee_dcvar_f(IDv);
    }
    if((IDv = variable_ID("PIAACMC_centobs1")) != -1)
    {
        centobs1 = coffee_dcvar_f(IDv);
    }
    if((IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = coffee_dcvar_f(IDv);
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }

    // FPM sag regularization control flag, do if == 1
    piaacmcparams.linopt_REGFPMSAG = 1; // default
    if((IDv = variable_ID("REGFPMSAG")) != -1)
    {
        piaacmcparams.linopt_REGFPMSAG =
            (long) coffee_dcvar_f(IDv) + 0.01;
    }

    // set current state for statistical tracking
    //    dcimg[IDstatus].array.UI16[0] = 0;

    // usual initialization
    {
        uint64_t initflag = INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__READCONF;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL;

        FUNC_CHECK_RETURN(init_piaacmcopticaldesign(fpmradld,
                          centobs0,
                          centobs1,
                          initflag,
                          NULL));
    }

    FUNC_CHECK_RETURN(makePIAAshapes());

    FUNC_CHECK_RETURN(init_piaacmcopticalsystem(0.0, 0.0));

    // set current state for statistical tracking
    //dcimg[IDstatus].array.UI16[0] = 1;

    // tracking diagnostic, giving the total flux in each plane
    WRITE_FULLFILENAME(fname, "%s/tmp_flux.txt", piaacmcparams.piaacmcconfdir);
    fp = fopen(fname, "r");
    if(fp == NULL)
    {
        FUNC_RETURN_FAILURE("Cannot open tmp_flux.txt");
    }
    for(long elem = 0; elem < piaacmcopticalsystem.NBelem; elem++)
    {
        //ret = fscanf(fp, "%lf %lf  %d\n", &tmplf1, &tmplf2, &tmpd1);
        double tmplf1, tmplf2;
        int fscanfcnt = fscanf(fp, "%lf %lf  %d\n", &tmplf1, &tmplf2, &tmpd1);
        if(fscanfcnt == EOF)
        {
            if(ferror(fp))
            {
                perror("fscanf");
            }
            else
            {
                fprintf(stderr, "Error: fscanf reached end of file\n");
            }
            fclose(fp);
            FUNC_RETURN_FAILURE("Call to fscanf failed");
        }
        else if(fscanfcnt != 3)
        {
            fprintf(stderr, "Error: fscanf matched %i items, 3 expected\n", fscanfcnt);
            fclose(fp);
            FUNC_RETURN_FAILURE("Call to fscanf failed");
        }

        // scale flux to current number of lambda
        piaacmcopticalsystem.flux[elem] = tmplf1 / tmpd1 * piaacmcopticalsystem.nblambda;
    }
    fclose(fp);

    piaacmcparams.LINOPT =
        1; // perform linear optimization after the switch exits
    // get the number of iterations
    /*
    if((IDv = variable_ID("PIAACMC_nbiter")) != -1)
    {
        NBiter = (long) coffee_dcvar_f(IDv) + 0.01;
    }
    else
    {
        NBiter = 50;    // default number of iterations
    }
    */
    // set current state for statistical tracking
    //dcimg[IDstatus].array.UI16[0] = 2;

    // get the FPMresp array computed in mode 11
    PIAACMCsimul_update_fnamedescr_conf();
    SNPRINTF_CHECK(fname, STRINGMAXLEN_FULLFILENAME, "%s/FPMresp%d.%s.fits",
                   piaacmcparams.piaacmcconfdir, piaacmcparams.SCORINGMASKTYPE,
                   piaacmcparams.fnamedescr_conf);

    FUNC_CHECK_RETURN(load_fits(fname, "FPMresp", 1, &IDfpmresp));

    piaacmcparams.vsize =
        dcimg[IDfpmresp].md[0].size[0]; // number of eval pts x2
    // make an array that holds the resulting light for evaluation point given the FPM solution, for each wavelenth
    //ID =
    FUNC_CHECK_RETURN(create_2Dimage_ID("imvect1",
                                        piaacmcparams.vsize,
                                        piaacmcopticaldesign.nblambda,
                                        NULL));

    // allocate arrays for fast routine
    // define convenient array variables
    piaacmcparams.fpmresp_array = dcimg[IDfpmresp].array.D;
    piaacmcparams.zonez_array =
        dcimg[piaacmcopticaldesign.zonezID].array.D;
    // allocate derivative of phase against thickness array
    piaacmcparams.dphadz_array =
        (double *) malloc(sizeof(double) * piaacmcopticaldesign.nblambda);
    // compute this derivative
    for(long k = 0; k < piaacmcopticaldesign.nblambda; k++)
    {
        // OPTICSMATERIALS_pha_lambda computes change in phase per unit thickness at specified wavelength
        // second arg is thickness, so 1.0 meters determines result per meter thickness
        piaacmcparams.dphadz_array[k] =
            OpticsMaterials_pha_lambda(piaacmcopticaldesign.fpmmaterial_code,
                                       1.0,
                                       piaacmcopticalsystem.lambdaarray[k]);
        printf("%ld  %g %g\n",
               k,
               piaacmcopticalsystem.lambdaarray[k],
               piaacmcparams.dphadz_array[k]);
    }
    piaacmcparams.outtmp_array = (double *) malloc(
                                     sizeof(double) *
                                     (piaacmcparams.vsize * piaacmcopticaldesign.nblambda +
                                      dcimg[piaacmcopticaldesign.zonezID].md[0].size[0]));

    // do the fast optimization using the results of mode 11
    piaacmcparams.computePSF_FAST_FPMresp = 1;

    // set current state for statistical tracking
    //dcimg[IDstatus].array.UI16[0] = 3;

    // read the contrast normalization factor into CnormFactor
    WRITE_FULLFILENAME(fname, "%s/CnormFactor.txt", piaacmcparams.piaacmcconfdir);
    FUNC_CHECK_RETURN(read_cnorm_factor(fname, &piaacmcparams.CnormFactor));
    // for each zone, add a random offset in range +- MODampl
    // this randomizes the starting point for each zone
    // dcimg[piaacmcopticaldesign.zonezID].array.D[k] is set in PIAACMCsimul_run()
    for(long k = 0; k < dcimg[piaacmcopticaldesign.zonezID].md[0].size[0];
            k++)
    {
        dcimg[piaacmcopticaldesign.zonezID].array.D[k] +=
            piaacmcparams.MODampl * (1.0 - 2.0 * ran1());
    }

    // set up optimization parameters for each zone
    // uses abstract specification of optimization parameters called paramval, ...
    piaacmcparams.linopt_number_param = 0;
    for(long mz = 0;
            mz < dcimg[piaacmcopticaldesign.zonezID].md[0].size[0];
            mz++)
    {
        // parameter type
        piaacmcparams.linopt_paramtype[piaacmcparams.linopt_number_param] =
            _DATATYPE_DOUBLE;
        // value: sag of each zone
        piaacmcparams.linopt_paramval[piaacmcparams.linopt_number_param] =
            &dcimg[piaacmcopticaldesign.zonezID].array.D[mz];
        // derivative step size
        piaacmcparams.linopt_paramdelta[piaacmcparams.linopt_number_param] =
            3.0e-9;
        // max parameter step size
        piaacmcparams.linopt_parammaxstep[piaacmcparams.linopt_number_param] =
            1.0e-6;
        // max and min allowed values for parameter
        piaacmcparams.linopt_parammin[piaacmcparams.linopt_number_param] =
            piaacmcopticaldesign.fpmminsag;
        piaacmcparams.linopt_parammax[piaacmcparams.linopt_number_param] =
            piaacmcopticaldesign.fpmmaxsag;
        // move on to next parameter
        piaacmcparams.linopt_number_param++;
    }
    piaacmcparams.PIAACMC_FPM_FASTDERIVATIVES =
        1; // for fast execution using analytic derivatives

    // set current state for statistical tracking
    //dcimg[IDstatus].array.UI16[0] = 4;
    // Now on to the actual optimization, after exit from the switch statement
    // I hope you have a lot of time...

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
