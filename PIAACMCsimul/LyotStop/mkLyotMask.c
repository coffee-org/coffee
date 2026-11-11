/**
 * @file    mkLyotMask.c
 * @brief   Synthesize binary Lyot stop masks from off-axis/on-axis intensity maps
 */

#include <stdint.h>
#include <stdio.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_memory/COREMOD_memory.h"
#include "info/info.h"

#include "PIAACMCsimul/LyotStop/mkLyotMask.h"
#include "PIAACMCsimul/PIAACMCsimul.h"

struct lyot_mask_inputs
{
    imageID  IDincoh;
    imageID  IDmc;
    imageID  IDzone;
    imageID  IDout;
    uint64_t xysize;
    double   throughput;
};

/**
 * normalize_image_sum_to_unity - Normalize a 2D float image in-place so its pixel sum equals 1.
 * @ID:     Target milk image ID
 * @xysize: Total number of pixels in the 2D image
 */
static void normalize_image_sum_to_unity(
    imageID  ID,
    uint64_t xysize)
{
    double val = 0.0;
    for (uint64_t ii = 0; ii < xysize; ii++)
    {
        val += dcimg[ID].array.F[ii];
    }
    for (uint64_t ii = 0; ii < xysize; ii++)
    {
        dcimg[ID].array.F[ii] /= val;
    }
}

/**
 * estimate_initial_ratio_threshold - Bracket off-axis/on-axis ratio for target throughput.
 * @in: Normalized pupil intensity maps, zone mask, and target throughput
 *
 * Return: Estimated threshold `rsl0` achieving approximately `in->throughput`.
 */
static double estimate_initial_ratio_threshold(const struct lyot_mask_inputs *in)
{
    long   NBiter = 100;
    double rsl    = 1.0;

    for (long iter = 0; iter < NBiter; iter++)
    {
        double val  = 0.0;
        double val1 = 0.0;

        for (uint64_t ii = 0; ii < in->xysize; ii++)
        {
            float incoh = dcimg[in->IDincoh].array.F[ii];
            float mc    = dcimg[in->IDmc].array.F[ii];
            if ((dcimg[in->IDzone].array.F[ii] > -1) && (incoh / mc > rsl))
            {
                val += incoh;
                val1 += mc;
                dcimg[in->IDout].array.F[ii] = 1.0;
            }
            else
            {
                dcimg[in->IDout].array.F[ii] = 0.0;
            }
        }

        printf("rsl = %f  ->  %f %f   (%f)\n", rsl, val, val1, in->throughput);
        if (val > in->throughput)
        {
            rsl *= 1.1;
        }
        else
        {
            rsl *= 0.9;
        }
    }

    return rsl;
}

/**
 * grid_search_dual_thresholds - Jointly scan ratio and on-axis flux cap to minimize stellar leak.
 * @in:           Normalized pupil intensity maps, zone mask, and target throughput
 * @rsl0:         Reference off-axis/on-axis ratio threshold
 * @v0:           99th-percentile on-axis intensity scale
 * @rsl_best_out: Output optimal off-axis/on-axis ratio threshold
 * @v_best_out:   Output optimal on-axis intensity ceiling
 */
static void grid_search_dual_thresholds(
    const struct lyot_mask_inputs *in,
    double                         rsl0,
    double                         v0,
    double                        *rsl_best_out,
    double                        *v_best_out)
{
    double bestval  = 1.0;
    double rsl_best = 0.0;
    double v_best   = 0.0;

    for (double rsl = 0.0 * rsl0; rsl < 2.0 * rsl0; rsl += 0.02 * rsl0)
    {
        for (double v = 0.00000001 * v0; v < 50.0 * v0; v *= 1.2)
        {
            double val  = 0.0;
            double val1 = 0.0;

            for (uint64_t ii = 0; ii < in->xysize; ii++)
            {
                float incoh = dcimg[in->IDincoh].array.F[ii];
                float mc    = dcimg[in->IDmc].array.F[ii];
                if ((dcimg[in->IDzone].array.F[ii] > -1) && (incoh / mc > rsl) && (mc < v))
                {
                    val += incoh;
                    val1 += mc;
                }
            }

            if ((val > in->throughput) && (val1 < bestval))
            {
                bestval  = val1;
                rsl_best = rsl;
                v_best   = v;
                printf("BEST SOLUTION: %.12lf / %.12lf    %.12lf / %.12lf  -> %.12lf  %.12lf\n",
                       rsl_best,
                       rsl0,
                       v_best,
                       v0,
                       val,
                       bestval);
            }
        }
    }

    *rsl_best_out = rsl_best;
    *v_best_out   = v_best;
}

/**
 * mkLyotMask - Synthesize optimal binary Lyot stop mask satisfying a geometric throughput target.
 * @IDincoh_name: Incoherent Lyot pupil intensity response to off-axis sources
 * @IDmc_name:    On-axis stellar residual Lyot pupil intensity image
 * @IDzone_name:  Annular zone map (`> -1` inside active pupil zones)
 * @throughput:   Minimum required off-axis geometric throughput (`0..1`)
 * @IDout_name:   Name of output binary Lyot mask image
 * @outID:        Optional output image ID pointer (may be NULL)
 *
 * Return: RETURN_SUCCESS on completion, or milk image creation error status.
 */
errno_t mkLyotMask(
    const char *restrict IDincoh_name,
    const char *restrict IDmc_name,
    const char *restrict IDzone_name,
    double               throughput,
    const char *restrict IDout_name,
    imageID             *outID)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %s %s %s %lf %s",
                     IDincoh_name,
                     IDmc_name,
                     IDzone_name,
                     throughput,
                     IDout_name);

    printf("IDincoh_name : %s   %ld\n", IDincoh_name, image_ID(IDincoh_name, dcimg, dcnimg));
    printf("IDmc_name    : %s   %ld\n", IDmc_name, image_ID(IDmc_name, dcimg, dcnimg));
    printf("IDzone_name  : %s   %ld\n", IDzone_name, image_ID(IDzone_name, dcimg, dcnimg));

    struct lyot_mask_inputs in;
    in.IDincoh    = image_ID(IDincoh_name, dcimg, dcnimg);
    in.IDmc       = image_ID(IDmc_name, dcimg, dcnimg);
    in.IDzone     = image_ID(IDzone_name, dcimg, dcnimg);
    in.throughput = throughput;

    uint32_t xsize = dcimg[in.IDmc].md[0].size[0];
    uint32_t ysize = dcimg[in.IDmc].md[0].size[1];
    in.xysize      = (uint64_t) xsize * ysize;

    FUNC_CHECK_RETURN(create_2Dimage_ID(IDout_name, xsize, ysize, &in.IDout));

    normalize_image_sum_to_unity(in.IDmc, in.xysize);
    normalize_image_sum_to_unity(in.IDincoh, in.xysize);

    double rsl0 = estimate_initial_ratio_threshold(&in);
    double v0   = img_percentile(IDmc_name, 0.99);
    printf("v0 = %lf\n", v0);

    double rsl_best = 0.0;
    double v_best   = 0.0;
    grid_search_dual_thresholds(&in, rsl0, v0, &rsl_best, &v_best);

    for (uint64_t ii = 0; ii < in.xysize; ii++)
    {
        float incoh = dcimg[in.IDincoh].array.F[ii];
        float mc    = dcimg[in.IDmc].array.F[ii];
        if ((dcimg[in.IDzone].array.F[ii] > -1) && (incoh / mc > rsl_best) && (mc < v_best))
        {
            dcimg[in.IDout].array.F[ii] = 1.0;
        }
        else
        {
            dcimg[in.IDout].array.F[ii] = 0.0;
        }
    }

    if (outID != NULL)
    {
        *outID = in.IDout;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
