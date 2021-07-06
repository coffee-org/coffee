/**
 * @file    PIAAACMCsimul_exec_optimize_lyot_stop_position.c
 * @brief   Optimize Lyot stop(s) conjugations
 *
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

#include "PIAACMCsimul.h"

#include "PIAACMCsimul_computePSF.h"
#include "init_piaacmcopticaldesign.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"

#include "PIAAshape/makePIAAshapes.h"




/**
 * ---
 *
 * ## Mode 2: Optimize focal plane mask transmission for monochromatic idealized PIAACMC

For monochromatic, idealized PIAACMC, find the scalar transimssion of the uniform focal plane mask
that provides best contrast in the evaluation zone

Very similar to the Lyot stop search in mode 1: iterative refined marching, changing the
the transmission value piaacmcopticaldesign.fpmaskamptransm, which
is between 0 and 1

Uses single on-axis light source
**/
errno_t exec_optimize_fpmtransmission()
{
    DEBUG_TRACE_FSTART();

    imageID IDv;
    double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;

    double range, stepsize;
    long NBiter = 1000;
    double parambest[10000]; // for scanning
    double paramref[10000];
    char fnamelog[STRINGMAXLEN_FULLFILENAME];

    printf("=================================== mode 002 ===================================\n");

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




    /// ### Initialize as in mode 0

    FUNC_CHECK_RETURN(
        init_piaacmcopticaldesign(0, fpmradld, centobs0, centobs1, 0, 1)
    );

    FUNC_CHECK_RETURN(makePIAAshapes());

    piaacmcopticalsystem.FOCMASKarray[0].mode = 1; // use 1-fpm


    /// ### Initialize search range and step
    range = 0.3;
    stepsize = range / 3.0;
    // initialized in PIAACMCsimul_initpiaacmcconf
    paramref[0] = piaacmcopticaldesign.fpmaskamptransm;
    NBiter = 6;


    /// ### Scan parameter value
    WRITE_FULLFILENAME(fnamelog, "%s/result_fpmt.log", piaacmcparams.piaacmcconfdir);
    {
        FILE * fp = fopen(fnamelog, "w");
        fclose(fp);
    }

    for(long iter = 0; iter < NBiter; iter++)
    {
        // starting point of march
        piaacmcopticaldesign.fpmaskamptransm = paramref[0] - range;

        // store current value as best
        parambest[0] = piaacmcopticaldesign.fpmaskamptransm;

        int loopOK = 1;
        double valbest = 1.0;

        /// While within the search loop :
        while(loopOK == 1)
        {
            double val;

            printf("\n\n\n");

            piaacmcparams.FORCE_CREATE_fpmza =
                1; // forces creation of new focal plane mask in the next two routines

            /// - Call PIAACMCsimul_initpiaacmcconf()
            FUNC_CHECK_RETURN(init_piaacmcopticaldesign(0, fpmradld, centobs0, centobs1, 0, 0));

            // compute on-axis PSF of all optical elements returning contrast in evaluation zone
            // ************************* need to do all piaacmcopticalsystem.NBelem?
            /// - call PIAACMCsimul_computePSF() to evaluate design
            FUNC_CHECK_RETURN(
                PIAACMCsimul_computePSF(
                    0.0, 0.0, 0, piaacmcopticalsystem.NBelem, 0, 0, 0, 0, &val)
            );

            if(val < valbest)
            {
                // we have a better contrast!  Store it
                parambest[0] = piaacmcopticaldesign.fpmaskamptransm;
                valbest = val;
            }

            /// - write entry to output log file
            {
                FILE * fp = fopen(fnamelog, "a");
                fprintf(fp, " %+011.8lf", piaacmcopticaldesign.fpmaskamptransm);
                fprintf(fp, " %12g  %8ld %12g %12g\n", val, iter, range, stepsize);
                fclose(fp);
            }

            /// -  increment parameter
            piaacmcopticaldesign.fpmaskamptransm += stepsize;

            // if we've reached the end of the range stop the loop
            if(piaacmcopticaldesign.fpmaskamptransm > paramref[0] + range + 0.001 * stepsize)
            {
                loopOK = 0;
            }
        }


        printf("BEST SOLUTION :  ");
        // store best solution
        paramref[0] = parambest[0];
        printf(" %lf", parambest[0]);

        printf(" %g\n", valbest);


        {
            FILE * fp = fopen(fnamelog, "a");
            fprintf(fp, "\n");
            fclose(fp);
        }
        // refine range and stepsize
        range *= 0.3;
        stepsize = range / 3.0;
    }
    // save final result
    piaacmcopticaldesign.fpmaskamptransm = parambest[0];

// why? **************************
    FUNC_CHECK_RETURN(init_piaacmcopticaldesign(0, fpmradld, centobs0, centobs1, 0, 0));

    // save final result to disk
    FUNC_CHECK_RETURN(PIAACMCsimul_savepiaacmcconf(piaacmcparams.piaacmcconfdir));


    piaacmcparams.FORCE_CREATE_fpmza = 0; // turning off to be good citizens

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
