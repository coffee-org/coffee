/**
 * @file    PIAACMCsimul_mkLyotMask.c
 * @brief   PIAA-type coronagraph design
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>

// milk includes
#include "CLIcore.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "info/info.h"

#include "PIAACMCsimul/PIAACMCsimul.h"

/**
 * @brief Make Lyot stop geometry
 *
 * Explores two thresholding methods applied together
 * (1) keeps pixels for which offaxisLight / onaxisLight > rsl
 * (2) keeps pixels for which onaxisLight < v0
 * Selects the mask that achieves the strongest on-axis rejection while satifying the throughput constraint
 *
 * @param[in]  IDincoh_name  Incoherent Lyot pupil intensity response to off-axis sources
 * @param[in]  IDmc_name     Intensity Lyot pupil image for on-axis source
 * @param[in]  IDzone_name
 * @param[in]  throughput
 * @param[out] IDout_name
 * @param[out] outID
 * @return errno_t
 */
errno_t mkLyotMask(const char *__restrict__ IDincoh_name,
                   const char *__restrict__ IDmc_name,
                   const char *__restrict__ IDzone_name,
                   double throughput,
                   const char *__restrict__ IDout_name,
                   imageID *outID)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %s %s %s %lf %s",
                     IDincoh_name,
                     IDmc_name,
                     IDzone_name,
                     throughput,
                     IDout_name);

    imageID ID1;
    imageID IDmc, IDincoh, IDzone;
    long    iter, NBiter;
    imageID IDout;
    //float sigma = 4.0;
    //int filter_size = 10;

    NBiter = 100;

    //filter_size = (long)(sigma * 2.0);

    printf("IDincoh_name : %s   %ld\n", IDincoh_name, image_ID(IDincoh_name, data.core.image, data.core.NB_MAX_IMAGE));
    printf("IDmc_name    : %s   %ld\n", IDmc_name, image_ID(IDmc_name, data.core.image, data.core.NB_MAX_IMAGE));
    printf("IDzone_name  : %s   %ld\n", IDzone_name, image_ID(IDzone_name, data.core.image, data.core.NB_MAX_IMAGE));

    //IDincoh = gauss_filter(IDincoh_name, "incohg", sigma, filter_size);
    IDincoh = image_ID(IDincoh_name, data.core.image, data.core.NB_MAX_IMAGE);

    IDmc = image_ID(IDmc_name, data.core.image, data.core.NB_MAX_IMAGE);
    //	IDmc = gauss_filter(IDmc_name, "mcg", sigma, filter_size);

    IDzone          = image_ID(IDzone_name, data.core.image, data.core.NB_MAX_IMAGE);
    uint32_t xsize  = data.core.image[IDmc].md[0].size[0];
    uint32_t ysize  = data.core.image[IDmc].md[0].size[1];
    uint64_t xysize = xsize;
    xysize *= ysize;

    FUNC_CHECK_RETURN(create_2Dimage_ID(IDout_name, xsize, ysize, &IDout));

    // normalize both images to 1.0
    {
        double val = 0.0;
        for(uint64_t ii = 0; ii < xysize; ii++)
        {
            val += data.core.image[IDmc].array.F[ii];
        }
        for(uint64_t ii = 0; ii < xysize; ii++)
        {
            data.core.image[IDmc].array.F[ii] /= val;
        }
    }

    {
        double val = 0.0;
        for(uint64_t ii = 0; ii < xysize; ii++)
        {
            val += data.core.image[IDincoh].array.F[ii];
        }
        for(uint64_t ii = 0; ii < xysize; ii++)
        {
            data.core.image[IDincoh].array.F[ii] /= val;
        }
    }

    // estimate iteratively rsl0, the threshold in offaxis/onaxis starlight that will achieve the goal throughput
    double rsl = 1.0;
    for(iter = 0; iter < NBiter; iter++)
    {
        double val  = 0.0;
        double val1 = 0.0;

        for(uint64_t ii = 0; ii < xysize; ii++)
        {
            if((data.core.image[IDzone].array.F[ii] > -1) &&
                    (data.core.image[IDincoh].array.F[ii] /
                     data.core.image[IDmc].array.F[ii] >
                     rsl))
            {
                val += data.core.image[IDincoh].array.F[ii];
                val1 += data.core.image[IDmc].array.F[ii];
                data.core.image[IDout].array.F[ii] = 1.0;
            }
            else
            {
                data.core.image[IDout].array.F[ii] = 0.0;
            }
        }
        printf("rsl = %f  ->  %f %f   (%f)\n", rsl, val, val1, throughput);
        if(val > throughput)  // too much light came through
        {
            rsl *= 1.1;
        }
        else
        {
            rsl *= 0.9;
        }
    }
    double rsl0 = rsl;

    // v0 = img_percentile("mcg", 0.99);
    double v0 = img_percentile(IDmc_name, 0.99);
    printf("v0 = %lf\n", v0);

    double bestval  = 1.0; // to be minized: total starlight transmitted
    double rsl_best = 0.0;
    double v_best   = 0.0;

    for(rsl = 0.0 * rsl0; rsl < 2.0 * rsl0; rsl += 0.02 * rsl0)
        for(double v = 0.00000001 * v0; v < 50.0 * v0; v *= 1.2)
        {
            double val  = 0.0;
            double val1 = 0.0;

            for(uint64_t ii = 0; ii < xysize; ii++)
            {
                if((data.core.image[IDzone].array.F[ii] > -1) &&
                        (data.core.image[IDincoh].array.F[ii] /
                         data.core.image[IDmc].array.F[ii] >
                         rsl) &&
                        (data.core.image[IDmc].array.F[ii] < v))
                {
                    val += data.core.image[IDincoh].array.F[ii];
                    val1 += data.core.image[IDmc].array.F[ii];
                }
            }

            if(val > throughput)
            {
                if(val1 < bestval)
                {
                    bestval  = val1;
                    rsl_best = rsl;
                    v_best   = v;
                    printf(
                        "BEST SOLUTION: %.12lf / %.12lf    %.12lf / %.12lf  -> "
                        "%.12lf  %.12lf\n",
                        rsl_best,
                        rsl0,
                        v_best,
                        v0,
                        val,
                        bestval);
                }
            }
        }

    for(uint64_t ii = 0; ii < xysize; ii++)
    {
        if((data.core.image[IDzone].array.F[ii] > -1) &&
                (data.core.image[IDincoh].array.F[ii] / data.core.image[IDmc].array.F[ii] >
                 rsl_best) &&
                (data.core.image[IDmc].array.F[ii] < v_best))
        {
            data.core.image[IDout].array.F[ii] = 1.0;
        }
        else
        {
            data.core.image[IDout].array.F[ii] = 0.0;
        }
    }

    if(0)
    {
        FUNC_CHECK_RETURN(create_2Dimage_ID("postLMim", xsize, ysize, &ID1));
        for(uint64_t ii = 0; ii < xysize; ii++)
        {
            data.core.image[ID1].array.F[ii] =
                data.core.image[IDmc].array.F[ii] * data.core.image[IDout].array.F[ii];
        }

        FUNC_CHECK_RETURN(save_fits("postLMim", "postLMim.fits"));
        FUNC_CHECK_RETURN(
            delete_image_ID("postLMim", DELETE_IMAGE_ERRMODE_WARNING));
    }

    if(outID != NULL)
    {
        *outID = IDout;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
