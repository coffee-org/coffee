/**
 * @file    PIAACMCsimul_linopt_scan.c
 * @brief   Multi-cutoff SVD search-direction synthesis and clamped gain line-search
 */

#include <math.h>
#include <stdio.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "linopt_imtools/linopt_imtools.h"

#include "PIAACMCsimul.h"
#include "PIAACMCsimul_computePSF.h"
#include "PIAACMCsimul_linopt.h"
#include "PIAACMCsimul_regularization.h"

/**
 * PIAACMCsimul_linopt_eval_objective - Compute PSF contrast plus active regularizers.
 * @val_total:    Output total regularized objective (`contrast + val0 + val1`)
 * @val_contrast: Optional output raw dark-hole contrast term (may be NULL)
 * @val_reg_piaa: Optional output PIAA shape regularization penalty (may be NULL)
 * @val_reg_fpm:  Optional output FPM sag regularization penalty (may be NULL)
 *
 * Return: RETURN_SUCCESS or PSF propagation error status.
 */
errno_t PIAACMCsimul_linopt_eval_objective(
    double *val_total,
    double *val_contrast,
    double *val_reg_piaa,
    double *val_reg_fpm)
{
    DEBUG_TRACE_FSTART();

    double cval = 0.0;
    double reg0 = 0.0;
    double reg1 = 0.0;

    FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(0.0,
                                              0.0,
                                              0,
                                              piaacmcopticalsystem.NBelem,
                                              0,
                                              piaacmcparams.computePSF_ResolvedTarget,
                                              piaacmcparams.computePSF_ResolvedTarget_mode,
                                              0,
                                              &cval));

    double total = cval;
    if (piaacmcparams.linopt_REGPIAASHAPES == 1)
    {
        reg0 = PIAACMCsimul_regularization_PIAAshapes_value();
        total += reg0;
    }
    if (piaacmcparams.linopt_REGFPMSAG == 1)
    {
        reg1 = PIAACMCsimul_regularization_fpmsag_value();
        total += reg1;
    }

    *val_total = total;
    if (val_contrast != NULL)
    {
        *val_contrast = cval;
    }
    if (val_reg_piaa != NULL)
    {
        *val_reg_piaa = reg0;
    }
    if (val_reg_fpm != NULL)
    {
        *val_reg_fpm = reg1;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * apply_clamped_step - Step active optimization parameters along `optcoeff` with clamping.
 * @IDdir:        Search direction image ID (`optcoeff`)
 * @scangain:     Scalar step gain along `-optcoeff`
 * @limflag_out:  Set to 1 if any parameter hit its per-step `linopt_parammaxstep` bound
 */
static void apply_clamped_step(
    imageID IDdir,
    double  scangain,
    int    *limflag_out)
{
    *limflag_out = 0;

    for (long i = 0; i < piaacmcparams.linopt_number_param; i++)
    {
        double step = -scangain * dcimg[IDdir].array.F[i] * piaacmcparams.linopt_paramdelta[i];
        double maxs = piaacmcparams.linopt_parammaxstep[i];

        if (step < -maxs)
        {
            printf("MIN LIMIT [%3ld   %20g]   %20g -> ",
                   i,
                   piaacmcparams.linopt_paramdelta[i],
                   step);
            step = -maxs;
            printf(" %20g\n", step);
            *limflag_out = 1;
        }
        if (step > maxs)
        {
            printf("MAX LIMIT [%3ld   %20g]   %20g -> ",
                   i,
                   piaacmcparams.linopt_paramdelta[i],
                   step);
            step = maxs;
            printf(" %20g\n", step);
            *limflag_out = 1;
        }

        if (piaacmcparams.linopt_paramtype[i] == _DATATYPE_FLOAT)
        {
            float cur = *(piaacmcparams.linopt_paramvalf[i]);
            if (cur + (float) step > piaacmcparams.linopt_parammax[i])
            {
                step = piaacmcparams.linopt_parammax[i] - cur;
            }
            if (cur + (float) step < piaacmcparams.linopt_parammin[i])
            {
                step = piaacmcparams.linopt_parammin[i] - cur;
            }
            piaacmcparams.linopt_paramdeltaval[i]   = step;
            *(piaacmcparams.linopt_paramvalf[i]) += (float) step;
        }
        else
        {
            double cur = *(piaacmcparams.linopt_paramval[i]);
            if (cur + step > piaacmcparams.linopt_parammax[i])
            {
                step = piaacmcparams.linopt_parammax[i] - cur;
            }
            if (cur + step < piaacmcparams.linopt_parammin[i])
            {
                step = piaacmcparams.linopt_parammin[i] - cur;
            }
            piaacmcparams.linopt_paramdeltaval[i] = step;
            *(piaacmcparams.linopt_paramval[i]) += step;
        }
    }
}

/**
 * revert_clamped_step - Restore optimization parameters to pre-step values.
 */
static void revert_clamped_step(void)
{
    for (long i = 0; i < piaacmcparams.linopt_number_param; i++)
    {
        if (piaacmcparams.linopt_paramtype[i] == _DATATYPE_FLOAT)
        {
            *(piaacmcparams.linopt_paramvalf[i]) -= (float) piaacmcparams.linopt_paramdeltaval[i];
        }
        else
        {
            *(piaacmcparams.linopt_paramval[i]) -= piaacmcparams.linopt_paramdeltaval[i];
        }
    }
}

/**
 * record_and_log_linescan_eval - Log candidate evaluation and snapshot best parameter vector.
 * @ctx:         Linear optimization iteration context
 * @alphareg:    SVD interpolation parameter in `[0, 1]`
 * @scangain:    Evaluated step gain
 * @contrastval: Total regularized objective value
 * @valContrast: Raw dark-hole contrast component
 * @val0:        PIAA shape regularization penalty
 * @val1:        FPM zone sag regularization penalty
 * @limflag:     Non-zero if step hit `parammaxstep`
 * @IDoptvec:    Image ID storing best parameter vector
 * @initbestval: In/out flag tracking whether `bestval` has been seeded
 */
static void record_and_log_linescan_eval(
    PIAACMC_linopt_ctx *ctx,
    double              alphareg,
    double              scangain,
    double              contrastval,
    double              valContrast,
    double              val0,
    double              val1,
    int                 limflag,
    imageID             IDoptvec,
    int                *initbestval)
{
    char fname[STRINGMAXLEN_FULLFILENAME];

    WRITE_FULLFILENAME(fname, "%s/linoptval.txt", piaacmcparams.piaacmcconfdir);
    FILE *fp = fopen(fname, "a");
    fprintf(fp,
            "##  [ %5ld / %5ld ]   %5.3f  %12lf  %12g   "
            "(reg = %12g [%1d] %12g [%1d]   contrast = %20g)       [%d] [%ld]",
            ctx->iter,
            piaacmcparams.linopt_NBiter,
            alphareg,
            scangain,
            contrastval,
            val0,
            piaacmcparams.linopt_REGPIAASHAPES,
            val1,
            piaacmcparams.linopt_REGFPMSAG,
            valContrast,
            limflag,
            piaacmcparams.linopt_number_param);

    if ((contrastval < ctx->bestval) || (*initbestval == 0))
    {
        for (long i = 0; i < piaacmcparams.linopt_number_param; i++)
        {
            if (piaacmcparams.linopt_paramtype[i] == _DATATYPE_FLOAT)
            {
                dcimg[IDoptvec].array.F[i] = *(piaacmcparams.linopt_paramvalf[i]);
            }
            else
            {
                dcimg[IDoptvec].array.F[i] = (float) *(piaacmcparams.linopt_paramval[i]);
            }
        }
        ctx->bestval = contrastval;
        fprintf(fp,
                (*initbestval == 0) ? " ===== START POINT =====\n" : "  -> BEST VECTOR =======\n");
        ctx->bestgain = scangain;
        *initbestval  = 1;
    }
    else
    {
        fprintf(fp, " bestval = %12g\n", ctx->bestval);
    }
    fclose(fp);
}

/**
 * build_interpolated_svd_direction - Blend `optcoeff0..2` into `optcoeff` for `alphareg`.
 * @IDstatus: Shared-memory status telemetry image ID
 * @alphareg: Piecewise-linear SVD conditioning parameter in `[0, 1]`
 *
 * Return: RETURN_SUCCESS or FITS/image error status.
 */
static errno_t build_interpolated_svd_direction(
    imageID IDstatus,
    double  alphareg)
{
    DEBUG_TRACE_FSTART();

    double acoeff0 = (1.0 - 2.0 * alphareg > 0.0) ? (1.0 - 2.0 * alphareg) : 0.0;
    double acoeff1 = 1.0 - fabs(2.0 * (alphareg - 0.5));
    double acoeff2 = (2.0 * alphareg - 1.0 > 0.0) ? (2.0 * alphareg - 1.0) : 0.0;

    dcimg[IDstatus].array.UI16[0] = 20;

    arith_image_cstmult("optcoeff0", acoeff0, "optcoeff0m");
    arith_image_cstmult("optcoeff1", acoeff1, "optcoeff1m");
    arith_image_cstmult("optcoeff2", acoeff2, "optcoeff2m");

    FUNC_CHECK_RETURN(save_fl_fits("optcoeff0", "optcoeff0.fits"));
    FUNC_CHECK_RETURN(save_fl_fits("optcoeff1", "optcoeff1.fits"));
    FUNC_CHECK_RETURN(save_fl_fits("optcoeff2", "optcoeff2.fits"));

    arith_image_add("optcoeff0m", "optcoeff1m", "optcoeff01m");
    arith_image_add("optcoeff01m", "optcoeff2m", "optcoeff");

    FUNC_CHECK_RETURN(delete_image_ID("optcoeff0m", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("optcoeff1m", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("optcoeff2m", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("optcoeff01m", DELETE_IMAGE_ERRMODE_WARNING));

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * scan_single_alphareg - Perform geometric gain line-search along blended `optcoeff`.
 * @ctx:          Linear optimization iteration context
 * @alphareg:     SVD blending parameter in `[0, 1]`
 * @IDoptvec:     Image ID storing best parameter vector
 * @initbestval:  In/out flag tracking whether `bestval` has been initialized
 * @NBlinoptgain: In/out maximum step count taken across `alphareg` scans
 *
 * Return: RETURN_SUCCESS or PSF evaluation error status.
 */
static errno_t scan_single_alphareg(
    PIAACMC_linopt_ctx *ctx,
    double              alphareg,
    imageID             IDoptvec,
    int                *initbestval,
    long               *NBlinoptgain)
{
    DEBUG_TRACE_FSTART();

    const float scangainfact = 1.2f;
    FUNC_CHECK_RETURN(build_interpolated_svd_direction(ctx->IDstatus, alphareg));

    dcimg[ctx->IDstatus].array.UI16[0] = 21;

    int    linscanOK   = 1;
    double scangain    = 0.001;
    double contrastval = 100000000000.0;
    int    k           = 0;

    while (linscanOK == 1)
    {
        int     limflag = 0;
        imageID IDdir   = image_ID("optcoeff", dcimg, dcnimg);

        dcimg[ctx->IDstatus].array.UI16[0] = 22;
        apply_clamped_step(IDdir, scangain, &limflag);

        double valold      = contrastval;
        double valContrast = 0.0;
        double val0        = 0.0;
        double val1        = 0.0;

        dcimg[ctx->IDstatus].array.UI16[0] = 23;
        FUNC_CHECK_RETURN(
            PIAACMCsimul_linopt_eval_objective(&contrastval, &valContrast, &val0, &val1));

        dcimg[ctx->IDstatus].array.UI16[0] = 24;
        record_and_log_linescan_eval(ctx,
                                     alphareg,
                                     scangain,
                                     contrastval,
                                     valContrast,
                                     val0,
                                     val1,
                                     limflag,
                                     IDoptvec,
                                     initbestval);

        revert_clamped_step();
        dcimg[ctx->IDstatus].array.UI16[0] = 25;

        k++;
        linscanOK = (contrastval < valold && k <= 90) ? 1 : 0;

        scangain += 0.001;
        scangain *= scangainfact;
    }

    if (k > *NBlinoptgain)
    {
        *NBlinoptgain = k;
    }

    delete_image_ID("optcoeff", DELETE_IMAGE_ERRMODE_WARNING);
    dcimg[ctx->IDstatus].array.UI16[0] = 26;

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * PIAACMCsimul_linopt_svd_linescan - Compute 3 SVD pseudo-inverses and scan `alphareg`.
 * @ctx:              Linear optimization iteration context
 * @IDoptvec_out:     Output image ID holding best parameter vector (`optvec`)
 * @NBlinoptgain_out: Output maximum step count across `alphareg` directions
 *
 * Return: RETURN_SUCCESS or error status.
 */
errno_t PIAACMCsimul_linopt_svd_linescan(
    PIAACMC_linopt_ctx *ctx,
    imageID            *IDoptvec_out,
    long               *NBlinoptgain_out)
{
    DEBUG_TRACE_FSTART();

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/linoptval.txt", piaacmcparams.piaacmcconfdir);
        FILE *fp = fopen(fname, "a");
        fprintf(fp, "### scanning gain \n");
        fprintf(fp, "### <alphareg>  <gain>  <contrast>\n");
        fclose(fp);
    }
    dcimg[ctx->IDstatus].array.UI16[0] = 17;

    printf("ref = 0.1   -- ");
    fflush(stdout);
    linopt_imtools_image_fitModes("vecDHref1D", "DHmodes", "DHmask", 0.1, "optcoeff0", 0, NULL);
    printf("- DONE\n");
    fflush(stdout);

    printf("ref = 0.01  -- ");
    fflush(stdout);
    linopt_imtools_image_fitModes("vecDHref1D", "DHmodes", "DHmask", 0.01, "optcoeff1", 0, NULL);
    printf("- DONE\n");
    fflush(stdout);

    printf("ref = 0.001 -- ");
    fflush(stdout);
    linopt_imtools_image_fitModes("vecDHref1D", "DHmodes", "DHmask", 0.001, "optcoeff2", 0, NULL);
    printf("- DONE\n");
    fflush(stdout);

    dcimg[ctx->IDstatus].array.UI16[0] = 18;

    arith_image_cstmult("optcoeff0", 0.0, "optvec");
    imageID IDoptvec     = image_ID("optvec", dcimg, dcnimg);
    int     initbestval  = 0;
    long    NBlinoptgain = 0;

    ctx->bestval                       = ctx->valref;
    ctx->bestgain                      = 0.0;
    dcimg[ctx->IDstatus].array.UI16[0] = 19;

    for (double alphareg = 0.0; alphareg < 1.01; alphareg += 0.2)
    {
        FUNC_CHECK_RETURN(
            scan_single_alphareg(ctx, alphareg, IDoptvec, &initbestval, &NBlinoptgain));
    }

    FUNC_CHECK_RETURN(delete_image_ID("optcoeff0", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("optcoeff1", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("DHmodes", DELETE_IMAGE_ERRMODE_WARNING));

    *IDoptvec_out     = IDoptvec;
    *NBlinoptgain_out = NBlinoptgain;

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
