/**
 * @file    PIAAACMCsimul_measure_transm_curve.c
 * @brief   Measure transmission as a function of angular separation
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
#include "PIAACMCsimul_initpiaacmcconf.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"

#include "PIAAshape/makePIAAshapes.h"


// externs
extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTSYST *optsyst;

extern OPTPIAACMCDESIGN *piaacmc;





/**
 * ---
 *
 * ## Mode 101: Measure transmission as a function of angular separation
 *
 */
errno_t PIAACMCsimul_measure_transm_curve()
{
    DEBUG_TRACE_FSTART();

    imageID IDv;
    double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;
    //double valref;
//    char fname[1500];
//    char fnametransm[1500];

    printf("=================================== mode 101 ===================================\n");

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

    piaacmcsimul_var.PIAACMC_fpmtype = 0; // idealized (default)
    if((IDv = variable_ID("PIAACMC_fpmtype")) != -1)
    {
        piaacmcsimul_var.PIAACMC_fpmtype = (int)(data.variable[IDv].value.f + 0.1);
    }

    piaacmcsimul_var.FORCE_CREATE_fpmza = 1;
    FUNC_CHECK_RETURN(
        PIAACMCsimul_initpiaacmcconf(piaacmcsimul_var.PIAACMC_fpmtype, fpmradld,
                                     centobs0, centobs1, 0, 1)
    );

    FUNC_CHECK_RETURN(makePIAAshapes(piaacmc));
    optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm


    {
        double cval = 0.0;
        FUNC_CHECK_RETURN(
            PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 1, 0, 0, 0, &cval)
        );
    }

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/psfi0test_x00_y00.fits", piaacmcsimul_var.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("psfi0", fname));
    }


    char fnametransm[STRINGMAXLEN_FULLFILENAME];
    PIAACMCsimul_update_fnamedescr();
    WRITE_FULLFILENAME(fnametransm, "%s/transmCurve_sm%d.%s.txt",
                       piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.SCORINGMASKTYPE,
                       piaacmcsimul_var.fnamedescr);

    FILE *fpt;
    fpt = fopen(fnametransm, "w");
    fclose(fpt);

    double stepld = 0.001;
    double xld;
    for(xld = 0.0; xld < 10.0; xld += stepld)
    {
        double val;

        //valref =
        {
            double cval = 0.0;
            FUNC_CHECK_RETURN(
                PIAACMCsimul_computePSF(xld, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0, &cval)
            );
        }


        imageID ID = image_ID("psfi0");
        //            sprintf(fname, "psfi0transm_%04.1f.fits", xld);
        //           save_fits("psfi0", fname);
        printf("ID = %ld\n", ID);
        uint32_t xsize = data.image[ID].md[0].size[0];
        uint32_t ysize = data.image[ID].md[0].size[1];
        uint32_t zsize = data.image[ID].md[0].size[2];
        printf("image size = %u %u %u\n", xsize, ysize, zsize);
        val = 0.0;

        for(uint32_t kk = 0; kk < zsize; kk++)
        {
            for(uint32_t ii = 0; ii < xsize; ii++)
                for(uint32_t jj = 0; jj < ysize; jj++)
                {
                    double dx, dy;

                    dx = 1.0 * ii - 0.5 * xsize;
                    dy = 1.0 * jj - 0.5 * ysize;
                    if((dx * dx + dy * dy) < 30.0 * 30.0)
                    {
                        val += data.image[ID].array.F[kk * xsize * ysize + jj * ysize + ii];
                    }
                }
        }
        val /= zsize;

        fpt = fopen(fnametransm, "a");
        fprintf(fpt, "%10f %.18f\n", xld, val);
        fclose(fpt);

        FUNC_CHECK_RETURN(
            delete_image_ID("psfi0", DELETE_IMAGE_ERRMODE_WARNING)
        );

        stepld = 0.001;
        stepld += 0.1 * xld;
        if(stepld > 0.2)
        {
            stepld = 0.2;
        }
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}



