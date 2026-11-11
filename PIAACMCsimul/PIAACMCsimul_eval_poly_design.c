/**
 * @file    PIAACMCsimul_eval_poly_design.c
 * @brief   Polychromatic coronagraph point-source and sensitivity evaluation driver (Mode 100)
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"
#include "PIAACMCsimul_computePSF.h"
#include "PIAACMCsimul_eval_poly_design.h"
#include "PIAACMCsimul_eval_poly_internal.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"
#include "init_piaacmcopticaldesign.h"

#include "PIAAshape/makePIAAshapes.h"

/**
 * eval_poly_init_design - Load Mode 100 CLI overrides, OPDerrC cube, and optical design.
 * @ctx: Polychromatic evaluation context
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t eval_poly_init_design(struct eval_poly_ctx *ctx)
{
    variableID IDv;

    if ((IDv = variable_ID("PIAACMC_centobs0")) != -1)
    {
        ctx->centobs0 = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_centobs1")) != -1)
    {
        ctx->centobs1 = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        ctx->fpmradld = coffee_dcvar_f(IDv);
        printf("MASK RADIUS = %lf lambda/D\n", ctx->fpmradld);
    }

    ctx->IDopderrC = image_ID("OPDerrC", dcimg, dcnimg);
    if (ctx->IDopderrC == -1)
    {
        FUNC_CHECK_RETURN(load_fits("OPDerrC.fits", "OPDerrC", 0, &ctx->IDopderrC));
    }
    ctx->nbOPDerr = (ctx->IDopderrC != -1) ? dcimg[ctx->IDopderrC].md[0].size[2] : 0;
    printf("Will add optional OPD error modes (%ld modes)\n", ctx->nbOPDerr);
    fflush(stdout);

    piaacmcparams.PIAACMC_fpmtype = 0;
    if ((IDv = variable_ID("PIAACMC_fpmtype")) != -1)
    {
        piaacmcparams.PIAACMC_fpmtype = (int) (coffee_dcvar_f(IDv) + 0.1);
    }

    piaacmcparams.FORCE_CREATE_fpmza = 1;
    {
        uint64_t initflag = INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__READCONF;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF;
        if (piaacmcparams.PIAACMC_fpmtype == 1)
        {
            initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL;
        }

        FUNC_CHECK_RETURN(init_piaacmcopticaldesign(
            ctx->fpmradld, ctx->centobs0, ctx->centobs1, initflag, NULL));
    }

    FUNC_CHECK_RETURN(makePIAAshapes());
    piaacmcopticalsystem.FOCMASKarray[0].mode = 1;

    if ((IDv = variable_ID("PIAACMC_ldoffset")) != -1)
    {
        ctx->ldoffset = coffee_dcvar_f(IDv);
    }
    printf("ldoffset = %f\n", ctx->ldoffset);

    return RETURN_SUCCESS;
}

/**
 * eval_poly_measure_offaxis_ref - Compute 5 l/D off-axis PSF peak normalization avpeak.
 * @ctx: Polychromatic evaluation context (populates dimensions, valref, and avpeak)
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t eval_poly_measure_offaxis_ref(struct eval_poly_ctx *ctx)
{
    FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(
        5.0, 0.0, 0, piaacmcopticalsystem.NBelem, 1, 0, 0, 0, &ctx->valref));

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/psfi0_x50_y00.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("psfi0", fname));
    }

    imageID ID  = image_ID("psfi0", dcimg, dcnimg);
    ctx->xsize  = dcimg[ID].md[0].size[0];
    ctx->ysize  = dcimg[ID].md[0].size[1];
    ctx->xysize = (uint64_t) ctx->xsize * ctx->ysize;
    ctx->zsize  = dcimg[ID].md[0].size[2];

    float *peakarray = (float *) malloc(sizeof(float) * ctx->zsize);
    if (peakarray == NULL)
    {
        FUNC_RETURN_FAILURE("malloc error");
    }

    for (uint32_t kk = 0; kk < ctx->zsize; kk++)
    {
        peakarray[kk] = 0.0f;
        for (uint64_t ii = 0; ii < ctx->xysize; ii++)
        {
            float val = dcimg[ID].array.F[kk * ctx->xysize + ii];
            if (val > peakarray[kk])
            {
                peakarray[kk] = val;
            }
        }
    }

    ctx->avpeak = 0.0f;
    for (uint32_t kk = 0; kk < ctx->zsize; kk++)
    {
        printf("peak %02ld  %10lf\n", (long) kk, peakarray[kk]);
        ctx->avpeak += peakarray[kk];
    }
    ctx->avpeak /= ctx->zsize;
    free(peakarray);

    FUNC_CHECK_RETURN(delete_image_ID("psfi0", DELETE_IMAGE_ERRMODE_WARNING));
    return RETURN_SUCCESS;
}

/**
 * bin_onaxis_psf_contrast - Radially bin on-axis PSF contrast and compute 2..6 l/D mean.
 * @ctx:     Polychromatic evaluation context
 * @curve:   Radial contrast accumulator array (length sepNBpt)
 * @cnt:     Radial sample count array (length sepNBpt)
 * @sepNBpt: Number of radial separation bins
 *
 * Return: Mean contrast within 2..6 lambda/D.
 */
static double bin_onaxis_psf_contrast(
    const struct eval_poly_ctx *ctx,
    double                     *curve,
    double                     *cnt,
    long                        sepNBpt)
{
    imageID ID       = image_ID("psfi0", dcimg, dcnimg);
    double  focscale = (2.0 * piaacmcopticaldesign.beamrad / piaacmcopticaldesign.pixscale) /
                       piaacmcopticaldesign.size;
    printf("focscale = %f\n", focscale);

    double aveC    = 0.0;
    long   aveCcnt = 0;

    for (uint32_t kk = 0; kk < ctx->zsize; kk++)
    {
        for (uint32_t ii = 0; ii < ctx->xsize; ii++)
        {
            for (uint32_t jj = 0; jj < ctx->ysize; jj++)
            {
                double xc = (1.0 * ii - 0.5 * ctx->xsize) * focscale;
                double yc = (1.0 * jj - 0.5 * ctx->ysize) * focscale;
                double rc = sqrt(xc * xc + yc * yc);
                long   ri = (long) (rc / ctx->eval_sepstepld - 0.5);
                if (ri < 0)
                {
                    ri = 0;
                }

                double cval = dcimg[ID].array.F[kk * ctx->xysize + jj * ctx->xsize + ii] /
                              ctx->avpeak;
                if (ri < sepNBpt)
                {
                    curve[ri] += cval;
                    cnt[ri]   += 1.0;
                }
                if ((rc > 2.0) && (rc < 6.0))
                {
                    aveC += cval;
                    aveCcnt++;
                }
            }
        }
    } // for kk

    return (aveCcnt > 0) ? (aveC / aveCcnt) : 0.0;
}

/**
 * eval_poly_eval_onaxis - Compute on-axis point-source PSF and write tt000 contrast files.
 * @ctx: Polychromatic evaluation context
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t eval_poly_eval_onaxis(struct eval_poly_ctx *ctx)
{
    FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(
        0.0, 0.0, 0, piaacmcopticalsystem.NBelem, 1, 0, 0, 1, &ctx->valref));

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/psfi0_x00_y00.fits", piaacmcparams.piaacmcconfdir);
        save_fits("psfi0", fname);
    }

    long    sepNBpt = (long) (ctx->eval_sepmaxld / ctx->eval_sepstepld);
    double *curve   = (double *) calloc((size_t) sepNBpt, sizeof(double));
    double *cnt     = (double *) calloc((size_t) sepNBpt, sizeof(double));
    if ((curve == NULL) || (cnt == NULL))
    {
        free(curve);
        free(cnt);
        FUNC_RETURN_FAILURE("malloc error");
    }

    double aveC_mean = bin_onaxis_psf_contrast(ctx, curve, cnt, sepNBpt);

    char fname[STRINGMAXLEN_FULLFILENAME];
    PIAACMCsimul_update_fnamedescr();
    SNPRINTF_CHECK(
        fname,
        STRINGMAXLEN_FULLFILENAME,
        "%s/ContrastCurve_tt000.%s.txt",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.fnamedescr);

    FILE *fp = fopen(fname, "w");
    if (fp != NULL)
    {
        fprintf(fp, "# Contrast curve\n");
        fprintf(fp, "# col 1: Angular separation [l/D]\n");
        fprintf(fp, "# col 2: Contrast value\n");
        fprintf(fp, "# col 3: Number of points used for contrast measurement\n\n");
        for (long ri = 0; ri < sepNBpt; ri++)
        {
            curve[ri] /= cnt[ri] + 0.000001;
            fprintf(fp, "%10f %10g %10g\n", ctx->eval_sepstepld * ri, curve[ri], cnt[ri]);
        }
        fclose(fp);
    }

    free(curve);
    free(cnt);

    PIAACMCsimul_update_fnamedescr();
    SNPRINTF_CHECK(
        fname,
        STRINGMAXLEN_FULLFILENAME,
        "%s/ContrastVal_tt000.%s.txt",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.fnamedescr);

    fp = fopen(fname, "w");
    if (fp != NULL)
    {
        fprintf(
            fp,
            "%10g %10g %4.2f %4.2f %4.2f %4.2f %04ld %02ld %d %03ld %03ld %02d %03ld %02d %d\n",
            ctx->valref,
            aveC_mean,
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
            (long) 0,
            piaacmcparams.computePSF_ResolvedTarget,
            piaacmcparams.computePSF_ResolvedTarget_mode);
        fclose(fp);
    }

    return RETURN_SUCCESS;
}

/**
 * PIAACMCsimul_eval_poly_design - Evaluate polychromatic contrast and sensitivity (Mode 100).
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
errno_t PIAACMCsimul_eval_poly_design(void)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FUNC");

    struct eval_poly_ctx ctx =
    {
        .fpmradld       = 0.95,
        .centobs0       = 0.3,
        .centobs1       = 0.2,
        .ldoffset       = 0.01,
        .eval_sepstepld = 0.2,
        .eval_sepmaxld  = 20.0,
        .IDopderrC      = -1,
    };

    printf("=================================== mode 100 ===================================\n");

    FUNC_CHECK_RETURN(eval_poly_init_design(&ctx));
    FUNC_CHECK_RETURN(eval_poly_measure_offaxis_ref(&ctx));
    FUNC_CHECK_RETURN(eval_poly_eval_onaxis(&ctx));
    FUNC_CHECK_RETURN(PIAACMCsimul_eval_poly_sens(&ctx));

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
