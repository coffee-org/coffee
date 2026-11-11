/**
 * @file    PIAACMCsimul_eval_poly_sens.c
 * @brief   Pointing jitter, OPD error sensitivity, and coherent/incoherent contrast evaluation
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_tools/COREMOD_tools.h"

#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"
#include "PIAACMCsimul_computePSF.h"
#include "PIAACMCsimul_eval_poly_internal.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"
#include "init_piaacmcopticalsystem.h"

/**
 * @brief Radially binned total, coherent, and incoherent contrast curves for extended-source eval.
 */
struct extsrc_curves
{
    double *tot;
    double *coh;
    double *inc;
    double *cnt;
    long    sepNBpt;
    double  aveC_mean;
};

/**
 * accum_psf_to_starim - Accumulate latest psfi0, psfre0, and psfim0 cubes into starim totals.
 * @ctx: Polychromatic evaluation context
 */
static void accum_psf_to_starim(struct eval_poly_ctx *ctx)
{
    imageID  ID    = image_ID("psfi0", dcimg, dcnimg);
    imageID  IDre  = image_ID("psfre0", dcimg, dcnimg);
    imageID  IDim  = image_ID("psfim0", dcimg, dcnimg);
    uint64_t nelem = ctx->xysize * ctx->zsize;

    for (uint64_t ii = 0; ii < nelem; ii++)
    {
        dcimg[ctx->IDps].array.F[ii]    += dcimg[ID].array.F[ii];
        dcimg[ctx->IDps_re].array.F[ii] += dcimg[IDre].array.F[ii];
        dcimg[ctx->IDps_im].array.F[ii] += dcimg[IDim].array.F[ii];
    }

    ctx->NBpt++;
    delete_image_ID("psfre0", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("psfim0", DELETE_IMAGE_ERRMODE_WARNING);
}

/**
 * eval_pointing_offset - Propagate one tip/tilt offset, save PSF FITS, and accumulate field.
 * @ctx:     Polychromatic evaluation context
 * @xld:     Tip offset [lambda/D]
 * @yld:     Tilt offset [lambda/D]
 * @relfits: Output FITS filename relative to piaacmcconfdir
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t eval_pointing_offset(
    struct eval_poly_ctx *ctx,
    double                xld,
    double                yld,
    const char           *relfits)
{
    double cval = 0.0;
    FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(
        xld, yld, 0, piaacmcopticalsystem.NBelem, 0, 0, 0, 0, &cval));
    ctx->valref += 0.25 * cval;

    char fname[STRINGMAXLEN_FULLFILENAME];
    WRITE_FULLFILENAME(fname, "%s/%s", piaacmcparams.piaacmcconfdir, relfits);
    save_fits("psfi0", fname);

    accum_psf_to_starim(ctx);
    return RETURN_SUCCESS;
}

/**
 * eval_opd_error_modes - Propagate each OPDerrC mode, save PSF cube, and accumulate field.
 * @ctx: Polychromatic evaluation context
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t eval_opd_error_modes(struct eval_poly_ctx *ctx)
{
    for (long OPDmode = 0; OPDmode < ctx->nbOPDerr; OPDmode++)
    {
        long    size = dcimg[ctx->IDopderrC].md[0].size[0];
        imageID IDopderr;
        FUNC_CHECK_RETURN(create_2Dimage_ID("opderr", size, size, &IDopderr));

        for (long ii = 0; ii < size * size; ii++)
        {
            dcimg[IDopderr].array.F[ii] = dcimg[ctx->IDopderrC].array.F[size * size * OPDmode + ii];
        }

        FUNC_CHECK_RETURN(init_piaacmcopticalsystem(0.0, 0.0));
        FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(
            0.0, 0.0, 0, piaacmcopticalsystem.NBelem, 0, 0, 0, 0, NULL));

        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(
            fname, "%s/psfi0_opderr%02ld.fits", piaacmcparams.piaacmcconfdir, OPDmode);
        save_fits("psfi0", fname);

        FUNC_CHECK_RETURN(delete_image_ID("opderr", DELETE_IMAGE_ERRMODE_WARNING));
        accum_psf_to_starim(ctx);
    }

    return RETURN_SUCCESS;
}

/**
 * build_and_save_coh_inc_cubes - Compute coherent/incoherent cubes and write extsrc FITS files.
 * @ctx: Polychromatic evaluation context
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t build_and_save_coh_inc_cubes(struct eval_poly_ctx *ctx)
{
    uint64_t nelem = ctx->xysize * ctx->zsize;

    for (uint64_t ii = 0; ii < nelem; ii++)
    {
        dcimg[ctx->IDps].array.F[ii]    /= ctx->NBpt;
        dcimg[ctx->IDps_re].array.F[ii] /= ctx->NBpt;
        dcimg[ctx->IDps_im].array.F[ii] /= ctx->NBpt;

        float re = dcimg[ctx->IDps_re].array.F[ii];
        float im = dcimg[ctx->IDps_im].array.F[ii];

        dcimg[ctx->IDps_COH].array.F[ii] = re * re + im * im;
        dcimg[ctx->IDps_INC].array.F[ii] =
            dcimg[ctx->IDps].array.F[ii] - dcimg[ctx->IDps_COH].array.F[ii];
    }

    long extsrc_code = (long) (-log10(ctx->ldoffset) * 10.0 + 0.1);
    char fname[STRINGMAXLEN_FULLFILENAME];

    PIAACMCsimul_update_fnamedescr();
    SNPRINTF_CHECK(
        fname,
        STRINGMAXLEN_FULLFILENAME,
        "%s/psfi0_extsrc%2ld_sm%d.%s.fits",
        piaacmcparams.piaacmcconfdir,
        extsrc_code,
        piaacmcparams.SCORINGMASKTYPE,
        piaacmcparams.fnamedescr);
    FUNC_CHECK_RETURN(save_fits("starim", fname));

    PIAACMCsimul_update_fnamedescr();
    SNPRINTF_CHECK(
        fname,
        STRINGMAXLEN_FULLFILENAME,
        "%s/psfi0INC_extsrc%2ld_sm%d.%s.fits",
        piaacmcparams.piaacmcconfdir,
        extsrc_code,
        piaacmcparams.SCORINGMASKTYPE,
        piaacmcparams.fnamedescr);
    FUNC_CHECK_RETURN(save_fits("starimINC", fname));

    PIAACMCsimul_update_fnamedescr();
    SNPRINTF_CHECK(
        fname,
        STRINGMAXLEN_FULLFILENAME,
        "%s/psfi0COH_extsrc%2ld_sm%d.%s.fits",
        piaacmcparams.piaacmcconfdir,
        extsrc_code,
        piaacmcparams.SCORINGMASKTYPE,
        piaacmcparams.fnamedescr);
    FUNC_CHECK_RETURN(save_fits("starimCOH", fname));

    return RETURN_SUCCESS;
}

/**
 * write_percentile_contrast_curve - Write radial p50/p20 percentile curve ContrastCurveP_extsrc.
 * @ctx:  Polychromatic evaluation context
 * @IDrc: 3D radial separation coordinate cube [lambda/D]
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t write_percentile_contrast_curve(
    const struct eval_poly_ctx *ctx,
    imageID                     IDrc)
{
    uint64_t nelem   = ctx->xysize * ctx->zsize;
    float   *rcarray = (float *) malloc(sizeof(float) * nelem);
    if (rcarray == NULL)
    {
        FUNC_RETURN_FAILURE("malloc error");
    }

    char fname[STRINGMAXLEN_FULLFILENAME];
    PIAACMCsimul_update_fnamedescr();
    SNPRINTF_CHECK(
        fname,
        STRINGMAXLEN_FULLFILENAME,
        "%s/ContrastCurveP_extsrc%2ld_sm%d.%s.txt",
        piaacmcparams.piaacmcconfdir,
        (long) (-log10(ctx->ldoffset) * 10.0 + 0.1),
        piaacmcparams.SCORINGMASKTYPE,
        piaacmcparams.fnamedescr);

    FILE *fp = fopen(fname, "w");
    if (fp == NULL)
    {
        free(rcarray);
        FUNC_RETURN_FAILURE("Cannot open ContrastCurveP file");
    }

    float drc      = 0.2f;
    long  rcbinmax = 100;
    for (long rcbin = 0; rcbin < rcbinmax; rcbin++)
    {
        float rcmin  = drc * rcbin;
        float rcmax  = drc * (rcbin + 1);
        long  rc_cnt = 0;

        for (uint64_t idx = 0; idx < nelem; idx++)
        {
            double rc = dcimg[IDrc].array.F[idx];
            if ((rc > rcmin) && (rc < rcmax))
            {
                rcarray[rc_cnt] = dcimg[ctx->IDps].array.F[idx] / ctx->avpeak;
                rc_cnt++;
            }
        }

        float p50 = 0.0f;
        float p20 = 0.0f;
        if (rc_cnt > 0)
        {
            quick_sort_float(rcarray, rc_cnt);
            p50 = rcarray[(long) (0.5 * rc_cnt)];
            p20 = rcarray[(long) (0.2 * rc_cnt)];
        }
        fprintf(fp, "%5.2f %12g %12g\n", drc * (rcbin + 0.5), p50, p20);
    } // for rcbin

    fclose(fp);
    free(rcarray);
    return RETURN_SUCCESS;
}

/**
 * bin_extsrc_contrast_cubes - Bin total/COH/INC radial contrast and write percentile curve.
 * @ctx:    Polychromatic evaluation context
 * @curves: Output radial contrast curves and 2..6 l/D mean
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t bin_extsrc_contrast_cubes(
    const struct eval_poly_ctx *ctx,
    struct extsrc_curves       *curves)
{
    double focscale = (2.0 * piaacmcopticaldesign.beamrad / piaacmcopticaldesign.pixscale) /
                      piaacmcopticaldesign.size;
    printf("focscale = %f\n", focscale);

    imageID IDrc;
    FUNC_CHECK_RETURN(create_3Dimage_ID("_tmp_rc", ctx->xsize, ctx->ysize, ctx->zsize, &IDrc));

    double aveC    = 0.0;
    long   aveCcnt = 0;

    for (uint32_t kk = 0; kk < ctx->zsize; kk++)
    {
        for (uint32_t ii = 0; ii < ctx->xsize; ii++)
        {
            for (uint32_t jj = 0; jj < ctx->ysize; jj++)
            {
                uint64_t pix = kk * ctx->xysize + jj * ctx->xsize + ii;
                double   xc  = (1.0 * ii - 0.5 * ctx->xsize) * focscale;
                double   yc  = (1.0 * jj - 0.5 * ctx->ysize) * focscale;
                double   rc  = sqrt(xc * xc + yc * yc);

                dcimg[IDrc].array.F[pix] = rc;

                long ri = (long) (rc / ctx->eval_sepstepld - 0.5);
                if (ri < 0)
                {
                    ri = 0;
                }
                if (ri < curves->sepNBpt)
                {
                    curves->tot[ri] += dcimg[ctx->IDps].array.F[pix] / ctx->avpeak;
                    curves->coh[ri] += dcimg[ctx->IDps_COH].array.F[pix] / ctx->avpeak;
                    curves->inc[ri] += dcimg[ctx->IDps_INC].array.F[pix] / ctx->avpeak;
                    curves->cnt[ri] += 1.0;
                }
                if ((rc > 2.0) && (rc < 6.0))
                {
                    aveC += dcimg[ctx->IDps].array.F[pix] / ctx->avpeak;
                    aveCcnt++;
                }
            }
        }
    } // for kk

    curves->aveC_mean = (aveCcnt > 0) ? (aveC / aveCcnt) : 0.0;

    errno_t err = write_percentile_contrast_curve(ctx, IDrc);
    FUNC_CHECK_RETURN(delete_image_ID("_tmp_rc", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(err);

    FUNC_CHECK_RETURN(delete_image_ID("starimINC", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("starimCOH", DELETE_IMAGE_ERRMODE_WARNING));

    return RETURN_SUCCESS;
}

/**
 * write_extsrc_summary_files - Write ContrastCurve_extsrc and ContrastVal_extsrc ASCII tables.
 * @ctx:    Polychromatic evaluation context
 * @curves: Radially binned total/COH/INC contrast curves
 */
static void write_extsrc_summary_files(
    const struct eval_poly_ctx *ctx,
    struct extsrc_curves       *curves)
{
    long extsrc_code = (long) (-log10(ctx->ldoffset) * 10.0 + 0.1);
    char fname[STRINGMAXLEN_FULLFILENAME];

    PIAACMCsimul_update_fnamedescr();
    SNPRINTF_CHECK(
        fname,
        STRINGMAXLEN_FULLFILENAME,
        "%s/ContrastCurve_extsrc%2ld_sm%d.%s.txt",
        piaacmcparams.piaacmcconfdir,
        extsrc_code,
        piaacmcparams.SCORINGMASKTYPE,
        piaacmcparams.fnamedescr);

    FILE *fp = fopen(fname, "w");
    if (fp != NULL)
    {
        fprintf(fp, "# Contrast curves\n");
        fprintf(fp, "# col 1 : separation [l/D]\n");
        fprintf(fp, "# col 2 : contrast, total intensity, radially averaged\n");
        fprintf(fp, "# col 3 : contrast, Coherent component, radially averaged\n");
        fprintf(fp, "# col 4 : contrast, incoherent component, radially averaged\n\n");
        for (long ri = 0; ri < curves->sepNBpt; ri++)
        {
            double denom     = curves->cnt[ri] + 0.000001;
            curves->tot[ri] /= denom;
            curves->coh[ri] /= denom;
            curves->inc[ri] /= denom;
            fprintf(
                fp,
                "%10f %10g %10g %10g %10g\n",
                ctx->eval_sepstepld * ri,
                curves->tot[ri],
                curves->cnt[ri],
                curves->coh[ri],
                curves->inc[ri]);
        }
        fclose(fp);
    }

    PIAACMCsimul_update_fnamedescr();
    SNPRINTF_CHECK(
        fname,
        STRINGMAXLEN_FULLFILENAME,
        "%s/ContrastVal_extsrc%2ld_sm%d.%s.txt",
        piaacmcparams.piaacmcconfdir,
        extsrc_code,
        piaacmcparams.SCORINGMASKTYPE,
        piaacmcparams.fnamedescr);

    fp = fopen(fname, "w");
    if (fp != NULL)
    {
        fprintf(fp, "# Contrast value\n");
        fprintf(
            fp,
            "# valref, aveC/aveCcnt, piaacmcopticaldesign.fpmaskradld, "
            "piaacmcparams.PIAACMC_MASKRADLD, piaacmcopticaldesign.centObs0, "
            "piaacmcopticaldesign.centObs1, (long) (piaacmcopticaldesign.lambda*1e9), "
            "(long) (piaacmcopticaldesign.lambdaB+0.1), piaacmcparams.PIAACMC_FPMsectors, "
            "piaacmcopticaldesign.NBrings, piaacmcopticaldesign.focmNBzone, "
            "piaacmcopticaldesign.nblambda, (long) (1000.0*ldoffset), "
            "piaacmcopticaldesign.fpmsagreg_coeff, piaacmcparams.computePSF_ResolvedTarget, "
            "piaacmcparams.computePSF_ResolvedTarget_mode\n#\n");
        fprintf(
            fp,
            "%10g %10g %4.2f %4.2f %4.2f %4.2f %04ld %02ld %d %03ld "
            "%03ld %02d %03ld %7.3f %02d %d\n",
            ctx->valref,
            curves->aveC_mean,
            piaacmcopticaldesign.fpmaskradld,
            piaacmcparams.PIAACMC_MASKRADLD,
            piaacmcopticaldesign.centObs0,
            piaacmcopticaldesign.centObs1,
            (long) (piaacmcopticaldesign.lambda * 1e9),
            (long) (piaacmcopticaldesign.lambdaB + 0.1),
            piaacmcparams.PIAACMC_FPMsectors,
            piaacmcopticaldesign.NBrings,
            piaacmcopticaldesign.focmNBzone,
            piaacmcopticaldesign.nblambda,
            (long) (1000.0 * ctx->ldoffset),
            piaacmcopticaldesign.fpmsagreg_coeff,
            piaacmcparams.computePSF_ResolvedTarget,
            piaacmcparams.computePSF_ResolvedTarget_mode);
        fclose(fp);
    }
}

/**
 * PIAACMCsimul_eval_poly_sens - Run pointing/OPD sensitivity and write extended-source curves.
 * @ctx: Polychromatic evaluation context
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
errno_t PIAACMCsimul_eval_poly_sens(struct eval_poly_ctx *ctx)
{
    long sz = piaacmcopticaldesign.size;

    FUNC_CHECK_RETURN(create_3Dimage_ID("starim", sz, sz, ctx->zsize, &ctx->IDps));
    FUNC_CHECK_RETURN(create_3Dimage_ID("starim_re", sz, sz, ctx->zsize, &ctx->IDps_re));
    FUNC_CHECK_RETURN(create_3Dimage_ID("starim_im", sz, sz, ctx->zsize, &ctx->IDps_im));
    FUNC_CHECK_RETURN(create_3Dimage_ID("starimCOH", sz, sz, ctx->zsize, &ctx->IDps_COH));
    FUNC_CHECK_RETURN(create_3Dimage_ID("starimINC", sz, sz, ctx->zsize, &ctx->IDps_INC));

    ctx->NBpt   = 0;
    ctx->valref = 0.0;

    FUNC_CHECK_RETURN(eval_pointing_offset(ctx, ctx->ldoffset, 0.0, "psfi0_p0.fits"));
    FUNC_CHECK_RETURN(eval_pointing_offset(ctx, -ctx->ldoffset, 0.0, "psfi0_m0.fits"));
    FUNC_CHECK_RETURN(eval_pointing_offset(ctx, 0.0, ctx->ldoffset, "psfi0_0p.fits"));
    FUNC_CHECK_RETURN(eval_pointing_offset(ctx, 0.0, -ctx->ldoffset, "psfi0_0m.fits"));

    FUNC_CHECK_RETURN(eval_opd_error_modes(ctx));
    FUNC_CHECK_RETURN(build_and_save_coh_inc_cubes(ctx));

    struct extsrc_curves curves =
    {
        .sepNBpt = (long) (ctx->eval_sepmaxld / ctx->eval_sepstepld),
    };
    curves.tot = (double *) calloc((size_t) curves.sepNBpt, sizeof(double));
    curves.coh = (double *) calloc((size_t) curves.sepNBpt, sizeof(double));
    curves.inc = (double *) calloc((size_t) curves.sepNBpt, sizeof(double));
    curves.cnt = (double *) calloc((size_t) curves.sepNBpt, sizeof(double));
    if ((curves.tot == NULL) || (curves.coh == NULL) ||
        (curves.inc == NULL) || (curves.cnt == NULL))
    {
        free(curves.tot);
        free(curves.coh);
        free(curves.inc);
        free(curves.cnt);
        FUNC_RETURN_FAILURE("malloc error");
    }

    errno_t err = bin_extsrc_contrast_cubes(ctx, &curves);
    if (err == RETURN_SUCCESS)
    {
        write_extsrc_summary_files(ctx, &curves);
    }

    free(curves.tot);
    free(curves.coh);
    free(curves.inc);
    free(curves.cnt);

    FUNC_CHECK_RETURN(err);
    return RETURN_SUCCESS;
}
