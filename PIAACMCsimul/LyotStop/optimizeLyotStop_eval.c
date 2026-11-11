/**
 * @file    optimizeLyotStop_eval.c
 * @brief   On-axis Fresnel smoothing and axial plane selection for Lyot stop optimization
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "image_filter/image_filter.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"

#include "optimizeLyotStop_internal.h"

struct smooth_work_ids
{
    imageID IDre;
    imageID IDim;
    imageID IDintg;
    imageID ID_LMintC;
    float   sigma;
    int     filter_size;
};

/**
 * accumulate_smoothed_plane - Accumulate Gaussian-smoothed polychromatic intensity at plane l.
 * @ctx:  Lyot stop optimization workspace context
 * @work: Pre-allocated 2D/3D scratch image IDs and filter parameters
 * @l:    Axial propagation plane index
 * @IDa1: Propagated amplitude cube ID
 * @IDp1: Propagated phase cube ID
 *
 * Return: RETURN_SUCCESS or milk image cleanup error code.
 */
static errno_t accumulate_smoothed_plane(
    const struct lyot_opt_ctx    *ctx,
    const struct smooth_work_ids *work,
    long                          l,
    imageID                       IDa1,
    imageID                       IDp1)
{
    uint64_t xysize = ctx->xysize;

    for (long k = 0; k < ctx->nblambda; k++)
    {
        for (uint64_t ii = 0; ii < xysize; ii++)
        {
            double amp = dcimg[IDa1].array.F[k * xysize + ii];
            double pha = dcimg[IDp1].array.F[k * xysize + ii];
            dcimg[work->IDre].array.F[ii] = amp * cos(pha);
            dcimg[work->IDim].array.F[ii] = amp * sin(pha);
        }

        imageID IDreg = gauss_filter("retmpim", "retmpimg", work->sigma, work->filter_size);
        imageID IDimg = gauss_filter("imtmpim", "imtmpimg", work->sigma, work->filter_size);

        for (uint64_t ii = 0; ii < xysize; ii++)
        {
            double re                   = dcimg[IDreg].array.F[ii];
            double im                   = dcimg[IDimg].array.F[ii];
            dcimg[work->IDintg].array.F[ii] = re * re + im * im;
        }

        imageID IDintgg =
            gauss_filter("tmpintg", "tmpintgg", 2.0 * work->sigma, work->filter_size);

        for (uint64_t ii = 0; ii < xysize; ii++)
        {
            dcimg[work->ID_LMintC].array.F[l * xysize + ii] += dcimg[IDintgg].array.F[ii];
        }

        FUNC_CHECK_RETURN(delete_image_ID("retmpimg", DELETE_IMAGE_ERRMODE_WARNING));
        FUNC_CHECK_RETURN(delete_image_ID("imtmpimg", DELETE_IMAGE_ERRMODE_WARNING));
        FUNC_CHECK_RETURN(delete_image_ID("tmpintgg", DELETE_IMAGE_ERRMODE_WARNING));
    }

    return RETURN_SUCCESS;
}

/**
 * lyot_opt_propagate_onaxis_cube - Propagate and smooth on-axis wavefront across NBz planes.
 * @ctx:        Lyot stop optimization workspace context
 * @IDamp_name: Input 2D/3D amplitude image name
 * @IDpha_name: Input 2D/3D phase image name
 *
 * Return: RETURN_SUCCESS or propagation/I/O error status.
 */
errno_t lyot_opt_propagate_onaxis_cube(
    const struct lyot_opt_ctx *ctx,
    const char       *restrict IDamp_name,
    const char       *restrict IDpha_name)
{
    DEBUG_TRACE_FSTART();

    struct smooth_work_ids work;
    FUNC_CHECK_RETURN(create_2Dimage_ID("retmpim", ctx->xsize, ctx->ysize, &work.IDre));
    FUNC_CHECK_RETURN(create_2Dimage_ID("imtmpim", ctx->xsize, ctx->ysize, &work.IDim));
    FUNC_CHECK_RETURN(
        create_3Dimage_ID("LMintC", ctx->xsize, ctx->ysize, ctx->NBz, &work.ID_LMintC));
    FUNC_CHECK_RETURN(create_2Dimage_ID("tmpintg", ctx->xsize, ctx->ysize, &work.IDintg));

    work.sigma = 0.01 * piaacmcopticaldesign.beamrad / piaacmcopticaldesign.pixscale;
    work.filter_size = (long)(work.sigma * 2.0);

    for (long l = 0; l < ctx->NBz; l++)
    {
        char nameamp[STRINGMAXLEN_IMGNAME];
        WRITE_IMAGENAME(nameamp, "LMPamp%02ld", l);

        char namepha[STRINGMAXLEN_IMGNAME];
        WRITE_IMAGENAME(namepha, "LMPpha%02ld", l);

        double zprop = ctx->zarray[l];
        FUNC_CHECK_RETURN(OptSystProp_propagateCube(&piaacmcopticalsystem,
                                                    0,
                                                    IDamp_name,
                                                    IDpha_name,
                                                    nameamp,
                                                    namepha,
                                                    zprop,
                                                    0));

        imageID IDa1 = image_ID(nameamp, dcimg, dcnimg);
        imageID IDp1 = image_ID(namepha, dcimg, dcnimg);

        FUNC_CHECK_RETURN(accumulate_smoothed_plane(ctx, &work, l, IDa1, IDp1));
        FUNC_CHECK_RETURN(delete_image_ID(nameamp, DELETE_IMAGE_ERRMODE_WARNING));
        FUNC_CHECK_RETURN(delete_image_ID(namepha, DELETE_IMAGE_ERRMODE_WARNING));
    } // for l

    FUNC_CHECK_RETURN(delete_image_ID("retmpim", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("imtmpim", DELETE_IMAGE_ERRMODE_WARNING));

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/LMintC.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("LMintC", fname));
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * compute_incoherent_zone_norms - Accumulate L1 and L_alpha norms of incoherent cube per zone.
 * @ctx:   Lyot stop optimization workspace context
 * @alpha: Concentration norm exponent
 */
static void compute_incoherent_zone_norms(
    const struct lyot_opt_ctx *ctx,
    double                     alpha)
{
    for (long l = 0; l < ctx->NBz; l++)
    {
        for (long m = 0; m < ctx->NBmasks; m++)
        {
            ctx->totarray[l * ctx->NBmasks + m]  = 0.0;
            ctx->tot2array[l * ctx->NBmasks + m] = 0.0;
        }
    }

    double rmin_obs = 0.9 * piaacmcopticaldesign.centObs1;
    for (long l = 0; l < ctx->NBz; l++)
    {
        for (uint64_t ii = 0; ii < ctx->xysize; ii++)
        {
            long m = (long)(dcimg[ctx->IDzone].array.F[ii] + 0.1);

            if ((m > -1) && (m < ctx->NBmasks) && (ctx->rarray[ii] < 1.0) &&
                (ctx->rarray[ii] > rmin_obs))
            {
                double val = dcimg[ctx->IDincohc].array.F[l * ctx->xysize + ii];
                ctx->totarray[l * ctx->NBmasks + m] += val;
                ctx->tot2array[l * ctx->NBmasks + m] += pow(val, alpha);
            }
        }
    }
}

/**
 * select_best_plane_for_mask - Find optimal conjugation plane for mask m and copy zone intensity.
 * @ctx:   Lyot stop optimization workspace context
 * @m:     Annular mask zone index
 * @alpha: Concentration norm exponent
 * @IDint: Smoothed on-axis 3D cube ID
 * @IDmc:  Combined on-axis 2D pupil image ID
 * @IDmc1: Combined off-axis incoherent 2D pupil image ID
 */
static void select_best_plane_for_mask(
    const struct lyot_opt_ctx *ctx,
    long                       m,
    double                     alpha,
    imageID                    IDint,
    imageID                    IDmc,
    imageID                    IDmc1)
{
    double valbest = 0.0;
    long   lbest   = 0;
    double zbest   = 0.0;

    for (long l = 0; l < ctx->NBz; l++)
    {
        long   idx = l * ctx->NBmasks + m;
        double val = ctx->tot2array[idx] / pow(ctx->totarray[idx], alpha);
        printf("MASK %ld   z(%ld)= %f  ->  %g   ( %g %g) \n",
               m,
               l,
               ctx->zarray[l],
               val,
               ctx->tot2array[idx],
               ctx->totarray[idx]);
        if (val > valbest)
        {
            valbest = val;
            zbest   = ctx->zarray[l];
            lbest   = l;
        }
    }

    printf(" ==========  MASK %ld   BEST CONJUGATION : %ld %f (%g)\n", m, lbest, zbest, valbest);
    piaacmcopticaldesign.LyotStop_zpos[m] = zbest;

    {
        char fname1[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname1, "%s/LyotMasks_zpos.txt", piaacmcparams.piaacmcconfdir);
        FILE *fp = fopen(fname1, "w");
        if (fp != NULL)
        {
            fprintf(fp, "%02ld %f\n", lbest, zbest);
            fclose(fp);
        }
    }

    for (uint64_t ii = 0; ii < ctx->xysize; ii++)
    {
        if (m == dcimg[ctx->IDzone].array.F[ii])
        {
            dcimg[IDmc].array.F[ii]  = dcimg[IDint].array.F[lbest * ctx->xysize + ii];
            dcimg[IDmc1].array.F[ii] = dcimg[ctx->IDincohc].array.F[lbest * ctx->xysize + ii];
        }
    }
}

/**
 * lyot_opt_select_conjugations - Evaluate zone norms and assemble optimal Lcomb/LcombOA maps.
 * @ctx:   Lyot stop optimization workspace context
 * @alpha: Concentration norm exponent
 *
 * Return: RETURN_SUCCESS or FITS/image error status.
 */
errno_t lyot_opt_select_conjugations(
    const struct lyot_opt_ctx *ctx,
    double                     alpha)
{
    DEBUG_TRACE_FSTART();

    compute_incoherent_zone_norms(ctx, alpha);

    imageID IDmc;
    imageID IDmc1;
    FUNC_CHECK_RETURN(create_2Dimage_ID("Lcomb", ctx->xsize, ctx->ysize, &IDmc));
    FUNC_CHECK_RETURN(create_2Dimage_ID("LcombOA", ctx->xsize, ctx->ysize, &IDmc1));

    for (uint64_t ii = 0; ii < ctx->xysize; ii++)
    {
        dcimg[IDmc1].array.F[ii] = 0.0;
    }

    imageID IDint = image_ID("LMintC", dcimg, dcnimg);
    for (long m = 0; m < ctx->NBmasks; m++)
    {
        select_best_plane_for_mask(ctx, m, alpha, IDint, IDmc, IDmc1);
    }

    if (piaacmcparams.PIAACMC_save == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];

        WRITE_FULLFILENAME(fname, "%s/Lcomb.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("Lcomb", fname));

        WRITE_FULLFILENAME(fname, "%s/LcombOA.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("LcombOA", fname));
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
