/**
 * @file    PIAACMCsimul_mkLyotMask.c
 * @brief   PIAA-type coronagraph design
 * 
 * Can design both APLCMC and PIAACMC coronagraphs
 *  
 * 
 * ## Change log
 * - 20180323  Guyon   minor syntax cleanup   
 * 
 * @author  O. Guyon
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
#include <stdlib.h>
#include <stdio.h>


// milk includes
#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "info/info.h"


#include "PIAACMCsimul/PIAACMCsimul.h"


/* =============================================================================================== */
/* =============================================================================================== */
/*                                  GLOBAL DATA DECLARATION                                        */
/* =============================================================================================== */
/* =============================================================================================== */

extern OPTPIAACMCDESIGN *piaacmc;




/* =============================================================================================== */
/* =============================================================================================== */
/*                                    FUNCTIONS SOURCE CODE                                        */
/* =============================================================================================== */
/* =============================================================================================== */

/// Make Lyot stop geometry
/// param[in] IDincoh_name   Incoherent Lyot pupil intensity response to off-axis sources
/// parampin] IDmc_name      Intensity Lyot pupil image for on-axis source
///
/// explores two thresholding methods applied together :
/// (1) keeps pixels for which offaxisLight / onaxisLight > rsl
/// (2) keeps pixels for which onaxisLight < v0
/// selects the mask that achieves the strongest on-axis rejection while satifying the throughput constraint

long PIAACMCsimul_mkLyotMask(
	const char *IDincoh_name, 
	const char *IDmc_name, 
	const char *IDzone_name, 
	double throughput, 
	const char *IDout_name
	)
{
    long ID, ID1;
    long IDmc, IDincoh, IDzone;
    double val, val1, v, v0, bestval, v_best, rsl_best;
    double rsl, rsl0;
    long iter, NBiter;
    long ii;
    long xsize = 0;
    long ysize = 0;
    long IDout;
    float sigma = 4.0;
    int filter_size = 10;

	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif


    NBiter = 100;

    sigma = 0.01*piaacmc[0].beamrad/piaacmc[0].pixscale;
    filter_size = (long) (sigma*2.0);

    printf("IDincoh_name : %s   %ld\n", IDincoh_name, image_ID(IDincoh_name));
    printf("IDmc_name    : %s   %ld\n", IDmc_name, image_ID(IDmc_name));
    printf("IDzone_name  : %s   %ld\n", IDzone_name, image_ID(IDzone_name));

    //IDincoh = gauss_filter(IDincoh_name, "incohg", sigma, filter_size);
    IDincoh = image_ID(IDincoh_name);

    IDmc = image_ID(IDmc_name);
    //	IDmc = gauss_filter(IDmc_name, "mcg", sigma, filter_size);

    IDzone = image_ID(IDzone_name);
    xsize = data.image[IDmc].md[0].size[0];
    ysize = data.image[IDmc].md[0].size[1];

    IDout = create_2Dimage_ID(IDout_name, xsize, ysize);

    // normalize both images to 1.0
    val = 0.0;
    for(ii=0; ii<xsize*ysize; ii++)
        val += data.image[IDmc].array.F[ii];
    for(ii=0; ii<xsize*ysize; ii++)
        data.image[IDmc].array.F[ii] /= val;

    val = 0.0;
    for(ii=0; ii<xsize*ysize; ii++)
        val += data.image[IDincoh].array.F[ii];
    for(ii=0; ii<xsize*ysize; ii++)
        data.image[IDincoh].array.F[ii] /= val;



    // estimate iteratively rsl0, the threshold in offaxis/onaxis starlight that will achieve the goal throughput
    rsl = 1.0;
    for(iter=0; iter<NBiter; iter++)
    {
        val = 0.0;
        val1 = 0.0;

        for(ii=0; ii<xsize*ysize; ii++)
        {
            if((data.image[IDzone].array.F[ii]>-1)&&(data.image[IDincoh].array.F[ii]/data.image[IDmc].array.F[ii]>rsl))
            {
                val += data.image[IDincoh].array.F[ii];
                val1 += data.image[IDmc].array.F[ii];
                data.image[IDout].array.F[ii] = 1.0;
            }
            else
                data.image[IDout].array.F[ii] = 0.0;
        }
        printf("rsl = %f  ->  %f %f   (%f)\n", rsl, val, val1, throughput);
        if(val>throughput) // too much light came through
            rsl *= 1.1;
        else
            rsl *= 0.9;
    }
    rsl0 = rsl;

    // v0 = img_percentile("mcg", 0.99);
    v0 = img_percentile(IDmc_name, 0.99);
    printf("v0 = %lf\n", v0);


    

    bestval = 1.0; // to be minized: total starlight transmitted
    for(rsl=0.0*rsl0; rsl< 2.0*rsl0; rsl+=0.02*rsl0)
        for(v=0.00000001*v0; v<50.0*v0; v*=1.2)
        {
            val = 0.0;
            val1 = 0.0;

            for(ii=0; ii<xsize*ysize; ii++)
            {
                if((data.image[IDzone].array.F[ii]>-1)&&(data.image[IDincoh].array.F[ii]/data.image[IDmc].array.F[ii]>rsl)&&(data.image[IDmc].array.F[ii]<v))
                {
                    val += data.image[IDincoh].array.F[ii];
                    val1 += data.image[IDmc].array.F[ii];
                }
            }

            if(val>throughput)
            {
                if(val1<bestval)
                {
                    bestval = val1;
                    rsl_best = rsl;
                    v_best = v;
                    printf("BEST SOLUTION: %.12lf / %.12lf    %.12lf / %.12lf  -> %.12lf  %.12lf\n", rsl_best, rsl0, v_best, v0, val, bestval);
                }

            }
        }

    for(ii=0; ii<xsize*ysize; ii++)
    {
        if((data.image[IDzone].array.F[ii]>-1)&&(data.image[IDincoh].array.F[ii]/data.image[IDmc].array.F[ii]>rsl_best)&&(data.image[IDmc].array.F[ii]<v_best))
            data.image[IDout].array.F[ii] = 1.0;
        else
            data.image[IDout].array.F[ii] = 0.0;
    }

    if(0)
    {
        ID1 = create_2Dimage_ID("postLMim", xsize, ysize);
        for(ii=0; ii<xsize*ysize; ii++)
            data.image[ID1].array.F[ii] = data.image[IDmc].array.F[ii]*data.image[IDout].array.F[ii];
        save_fits("postLMim", "!postLMim.fits");
        delete_image_ID("postLMim");
    }


    return(IDout);
}





