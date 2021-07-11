/**
 * @file    PIAAACMCsimul_exec_optimize_lyot_stop_position.c
 * @brief   Optimize Lyot stop(s) conjugations
 *
 *
 */



// System includes
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>



// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"

#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"

#include "PIAACMCsimul_computePSF.h"
#include "init_piaacmcopticalsystem.h"
#include "init_piaacmcopticaldesign.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"

#include "PIAAshape/makePIAAshapes.h"




/**
* ---
*
* ## Mode 1: Optimize Lyot stop positions

Lyot stop positions are encoded as piaacmc[0].LyotStop_zpos

there can be multiple LyotStop_zpos

Vary these zpos, looking for the best contrast returned by PIAACMCsimul_computePSF

Search is performed by iterative refined marching
**/
errno_t exec_optimize_lyot_stop_position()
{
    DEBUG_TRACE_FSTART();

    imageID IDv;
    double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;

    long NBiter = 1000;

    double parambest[10000]; // for scanning
    double paramref[10000];

    char fnamelog[STRINGMAXLEN_FULLFILENAME];

    printf("=================================== mode 001 ===================================\n");
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



    // init as in mode 0
    {
        uint64_t initflag = INIT_PIAACMCOPTICALDESIGN_MODE__READCONF;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF;
        FUNC_CHECK_RETURN(
            init_piaacmcopticaldesign(
                fpmradld,
                centobs0,
                centobs1,
                initflag,
                NULL
            )
        );
    }

    FUNC_CHECK_RETURN(
        makePIAAshapes()
    );

    piaacmcopticalsystem.FOCMASKarray[0].mode = 1; // use 1-fpm

    // initialization
    // necessary to initialize optical design
    FUNC_CHECK_RETURN(
        init_piaacmcopticalsystem(0.0, 0.0)
    );

    // set initial lyot stop marching range (current position +- range)
    double range, stepsize;
    {
        variableID IDv;
        if((IDv = variable_ID("PIAACMC_lsoptrange")) != -1)
        {
            range = data.variable[IDv].value.f;    // from cli
        }
        else
        {
            range = 3.0;    // default, in meters
        }
        stepsize = range / 3.0; // initial march stepsize
    }


    // store initial Lyot stop positions
    // NBLyotStop = length(LyotStop_zpos)
    for(long ls = 0; ls < piaacmcopticaldesign.NBLyotStop; ls++)
    {
        paramref[ls] = piaacmcopticaldesign.LyotStop_zpos[ls];
    }
    NBiter = 4; // number of iterations

    // start up a log
    WRITE_FULLFILENAME(fnamelog, "%s/result_LMpos.log", piaacmcparams.piaacmcconfdir);
    {
        FILE * fp = fopen(fnamelog, "w");
        fclose(fp);
    }


    // pick another initial march stepsize
    stepsize = range / 5.0;
    // start the iterative refined march
    for(long iter = 0; iter < NBiter; iter++)
    {
        // for each Lyot stop, find its best position
        for(long ls = 0; ls < piaacmcopticaldesign.NBLyotStop; ls++)
        {
            // start position for march.  paramref is current best value
            piaacmcopticaldesign.LyotStop_zpos[ls] = paramref[ls] - range;
            // current best position
            parambest[ls] = piaacmcopticaldesign.LyotStop_zpos[ls];

            // loopOK = 1;
            double valbest = 1.0;

            // march to the other other end of range
            while(piaacmcopticaldesign.LyotStop_zpos[ls] < paramref[ls] + range)
            {
                long elem;
                double val;

                long elem0 = 6; // elem0 is the starting point of the PSF propagation.  This is a staring default

                // look for the element called "Lyot mask 0" as the actual starting point
                elem0 = -1;
                printf("Number of elements = %ld\n", piaacmcopticalsystem.NBelem);
                assert(piaacmcopticalsystem.NBelem > 0);

                for(elem = 0; elem < piaacmcopticalsystem.NBelem; elem++)
                {
                    printf("elem %ld :  %s\n", elem, piaacmcopticalsystem.name[elem]);
                    if(strcmp("Lyot mask 0", piaacmcopticalsystem.name[elem]) == 0)
                    {
                        elem0 = elem;
                    }
                }
                assert(elem0 != -1);  // throw a message if this was not found


                piaacmcopticalsystem.keepMem[elem0] = 1; // save this element and reuse

                // compute the PSF for this Lyot stop position, returning contrast in the evaluation zone
                FUNC_CHECK_RETURN(
                    PIAACMCsimul_computePSF(
                        0.0, 0.0, elem0, piaacmcopticalsystem.NBelem, 0, 0, 0, 0, &val)
                );

                // if this is the best contrast for this stop, save it for it and the position of this stop
                if(val < valbest)
                {
                    parambest[ls] = piaacmcopticaldesign.LyotStop_zpos[ls];
                    valbest = val;
                }

                // say what's happening
                {
                    FILE * fp = fopen(fnamelog, "a");
                    for(long ls1 = 0; ls1 < piaacmcopticaldesign.NBLyotStop; ls1++)
                    {
                        fprintf(fp, " %lf", piaacmcopticaldesign.LyotStop_zpos[ls1]);
                    }
                    fprintf(fp, " %g\n", val);
                    fclose(fp);
                }

                // march along by the step size
                piaacmcopticaldesign.LyotStop_zpos[ls] += stepsize;
            }
            printf("BEST SOLUTION :  ");
            paramref[ls] = parambest[ls]; // update best position for this stop
            piaacmcopticaldesign.LyotStop_zpos[ls] =
                paramref[ls]; // store in case this is last iteration
            printf(" %lf", parambest[ls]);
            printf(" %g\n", valbest);
        }

        {
            FILE * fp = fopen(fnamelog, "a");
            fprintf(fp, "\n");
            fclose(fp);
        }

        // reduce the range and stepsize, refining the march
        range *= 0.3;
        stepsize = range / 3.0;
    }
    // store all best positions  Done!!
    for(long ls = 0; ls < piaacmcopticaldesign.NBLyotStop; ls++)
    {
        piaacmcopticaldesign.LyotStop_zpos[ls] = parambest[ls];
    }
    PIAACMCsimul_savepiaacmcconf(piaacmcparams.piaacmcconfdir);


    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

