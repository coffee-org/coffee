/**
 * @file    PIAAACMCsimul_measure_transm_curve.c
 * @brief   Measure transmission as a function of angular separation
 * 
 *  
 * @author  O. Guyon
 * @date    25 nov 2017
 *
 * 
 * @bug No known bugs.
 * 
 */


/* =============================================================================================== */
/* =============================================================================================== */
/*                                        HEADER FILES                                             */
/* =============================================================================================== */
/* =============================================================================================== */

// System includes
#include <stdio.h>
#include <stdlib.h>



// milk includes
#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"



/* =============================================================================================== */
/* =============================================================================================== */
/*                                      DEFINES, MACROS                                            */
/* =============================================================================================== */
/* =============================================================================================== */



/* =============================================================================================== */
/* =============================================================================================== */
/*                                  GLOBAL DATA DECLARATION                                        */
/* =============================================================================================== */
/* =============================================================================================== */

// externs
extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTSYST *optsyst;

extern OPTPIAACMCDESIGN *piaacmc;





/* =============================================================================================== */
/* =============================================================================================== */
/*                                    FUNCTIONS SOURCE CODE                                        */
/* =============================================================================================== */
/* =============================================================================================== */

/**
 * ---
 *
 * ## Mode 101: Measure transmission as a function of angular separation
 *
 */
int PIAACMCsimul_measure_transm_curve()
{
    long IDv;
    double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;
    double ldoffset;
    double valref;
    char fname[1500];
    char fnametransm[1500];

    printf("=================================== mode 101 ===================================\n");

    if( (IDv = variable_ID("PIAACMC_centobs0")) != -1)
        centobs0 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_centobs1")) != -1)
        centobs1 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }


    printf("101: transm as a function of angular separation  ldoffset = %f\n", ldoffset);

    piaacmcsimul_var.PIAACMC_fpmtype = 0; // idealized (default)
    if((IDv=variable_ID("PIAACMC_fpmtype"))!=-1)
        piaacmcsimul_var.PIAACMC_fpmtype = (int) (data.variable[IDv].value.f + 0.1);

    piaacmcsimul_var.FORCE_CREATE_fpmza = 1;
    PIAACMCsimul_initpiaacmcconf(piaacmcsimul_var.PIAACMC_fpmtype, fpmradld, centobs0, centobs1, 0, 1);

    PIAACMCsimul_makePIAAshapes(piaacmc, 0);
    optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm


    valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 1, 0, 0, 0);
    sprintf(fname,"!%s/psfi0test_x00_y00.fits", piaacmcsimul_var.piaacmcconfdir);
    save_fits("psfi0", fname);

	PIAACMCsimul_update_fnamedescr();
    sprintf(fnametransm, "%s/transmCurve_sm%d.%s.txt", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.fnamedescr);

    FILE *fpt;
    fpt = fopen(fnametransm, "w");
    fclose(fpt);

    double stepld = 0.001;
    double xld;
    for(xld=0.0; xld<10.0; xld+=stepld)
    {
        long ID;
        long xsize, ysize, zsize;
        double val;

        valref = PIAACMCsimul_computePSF(xld, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);
        ID = image_ID("psfi0");
        //            sprintf(fname, "!psfi0transm_%04.1f.fits", xld);
        //           save_fits("psfi0", fname);
        printf("ID = %ld\n", ID);
        xsize = data.image[ID].md[0].size[0];
        ysize = data.image[ID].md[0].size[1];
        zsize = data.image[ID].md[0].size[2];
        printf("image size = %ld %ld %ld\n", xsize, ysize, zsize);
        val = 0.0;

        long kk;
        for(kk=0; kk<zsize; kk++)
        {
            long ii, jj;
            for(ii=0; ii<xsize; ii++)
                for(jj=0; jj<ysize; jj++)
                {
                    double dx, dy;

                    dx = 1.0*ii-0.5*xsize;
                    dy = 1.0*jj-0.5*ysize;
                    if((dx*dx+dy*dy)<30.0*30.0)
                        val += data.image[ID].array.F[kk*xsize*ysize+jj*ysize+ii];
                }
        }
        val /= zsize;

        fpt = fopen(fnametransm, "a");
        fprintf(fpt, "%10f %.18f\n", xld, val);
        fclose(fpt);
        delete_image_ID("psfi0");

        stepld = 0.001;
        stepld += 0.1*xld;
        if(stepld>0.2)
            stepld = 0.2;
    }

    return 0;
}



