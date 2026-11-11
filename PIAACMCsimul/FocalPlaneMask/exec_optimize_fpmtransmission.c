/**
 * @file    exec_optimize_fpmtransmission.c
 * @brief   Optimize idealized monochromatic focal plane mask transmission
 */

#include <stdint.h>
#include <stdio.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_memory/COREMOD_memory.h"

#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"
#include "PIAACMCsimul_computePSF.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"
#include "init_piaacmcopticaldesign.h"

#include "FocalPlaneMask/exec_optimize_fpmtransmission.h"
#include "PIAAshape/makePIAAshapes.h"

struct fpmtransm_opt_ctx
{
    double fpmradld;
    double centobs0;
    double centobs1;
    double range;
    double stepsize;
    double paramref;
    double parambest;
    char   fnamelog[STRINGMAXLEN_FULLFILENAME];
};

/**
 * fpmtransm_load_cli_vars - Populate optical geometry overrides from CLI variables
 * @ctx: Optimization context receiving central obstruction and mask radius values.
 */
static void fpmtransm_load_cli_vars(struct fpmtransm_opt_ctx *ctx)
{
    imageID IDv;

    ctx->fpmradld = 0.95;
    ctx->centobs0 = 0.3;
    ctx->centobs1 = 0.2;

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
}

/**
 * fpmtransm_eval_step - Evaluate on-axis PSF contrast at current FPM transmission
 * @ctx:     Active scan context holding geometry and log file path.
 * @iter:    Outer refinement iteration index for logging.
 * @val_out: Output pointer receiving measured dark-hole contrast.
 *
 * Return: RETURN_SUCCESS on valid evaluation and log append, error code otherwise.
 */
static errno_t fpmtransm_eval_step(
    const struct fpmtransm_opt_ctx *ctx,
    long                            iter,
    double                         *val_out)
{
    printf("\n\n\n");

    piaacmcparams.FORCE_CREATE_fpmza = 1;
    DEBUG_TRACEPOINT(
        "computePSF piaacmcopticaldesign.fpmaskamptransm %lf",
        piaacmcopticaldesign.fpmaskamptransm);

    FUNC_CHECK_RETURN(init_piaacmcopticaldesign(
        ctx->fpmradld,
        ctx->centobs0,
        ctx->centobs1,
        INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT,
        NULL));

    FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(
        0.0,
        0.0,
        0,
        piaacmcopticalsystem.NBelem,
        0,
        0,
        0,
        0,
        val_out));

    FILE *fp = fopen(ctx->fnamelog, "a");
    if (fp == NULL)
    {
        FUNC_RETURN_FAILURE("Cannot append to file %s", ctx->fnamelog);
    }
    fprintf(fp, " %+011.8lf", piaacmcopticaldesign.fpmaskamptransm);
    fprintf(fp, " %12g  %8ld %12g %12g\n", *val_out, iter, ctx->range, ctx->stepsize);
    fclose(fp);

    return RETURN_SUCCESS;
}

/**
 * fpmtransm_scan_iter - Run one 1D grid-search refinement pass over FPM transmission
 * @ctx:  Optimization context updated in-place with best parameter and narrowed range.
 * @iter: Current refinement pass index.
 *
 * Return: RETURN_SUCCESS on pass completion, error code otherwise.
 */
static errno_t fpmtransm_scan_iter(
    struct fpmtransm_opt_ctx *ctx,
    long                      iter)
{
    piaacmcopticaldesign.fpmaskamptransm = ctx->paramref - ctx->range;
    ctx->parambest                       = piaacmcopticaldesign.fpmaskamptransm;

    int    loopOK  = 1;
    double valbest = 1.0;

    while (loopOK == 1)
    {
        double val = 1.0;

        FUNC_CHECK_RETURN(fpmtransm_eval_step(ctx, iter, &val));
        if (val < valbest)
        {
            ctx->parambest = piaacmcopticaldesign.fpmaskamptransm;
            valbest        = val;
        }

        piaacmcopticaldesign.fpmaskamptransm += ctx->stepsize;
        double maxval = ctx->paramref + ctx->range + 0.001 * ctx->stepsize;
        if (piaacmcopticaldesign.fpmaskamptransm > maxval)
        {
            loopOK = 0;
        }
    }

    ctx->paramref = ctx->parambest;
    printf("BEST SOLUTION :   %lf %g\n", ctx->parambest, valbest);

    FILE *fp = fopen(ctx->fnamelog, "a");
    if (fp == NULL)
    {
        FUNC_RETURN_FAILURE("Cannot append to file %s", ctx->fnamelog);
    }
    fprintf(fp, "\n");
    fclose(fp);

    ctx->range *= 0.3;
    ctx->stepsize = ctx->range / 3.0;

    return RETURN_SUCCESS;
}

/**
 * exec_optimize_fpmtransmission - Optimize idealized monochromatic FPM transmission (Mode 2)
 *
 * Performs iterative multi-scale 1D parameter marching over
 * piaacmcopticaldesign.fpmaskamptransm to minimize evaluation-zone contrast
 * for an on-axis point source, then persists the optimal configuration.
 *
 * Return: RETURN_SUCCESS on completion, error code otherwise.
 */
errno_t exec_optimize_fpmtransmission(void)
{
    DEBUG_TRACE_FSTART();

    struct fpmtransm_opt_ctx ctx;
    const long               NBiter = 6;

    printf("=================================== mode 002 ===================================\n");
    fpmtransm_load_cli_vars(&ctx);

    {
        uint64_t initflag = INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__READCONF;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF;
        FUNC_CHECK_RETURN(init_piaacmcopticaldesign(
            ctx.fpmradld,
            ctx.centobs0,
            ctx.centobs1,
            initflag,
            NULL));
    }

    FUNC_CHECK_RETURN(makePIAAshapes());
    piaacmcopticalsystem.FOCMASKarray[0].mode = 1;

    ctx.range     = 0.3;
    ctx.stepsize  = ctx.range / 3.0;
    ctx.paramref  = piaacmcopticaldesign.fpmaskamptransm;
    ctx.parambest = ctx.paramref;

    WRITE_FULLFILENAME(ctx.fnamelog, "%s/result_fpmt.log", piaacmcparams.piaacmcconfdir);
    {
        FILE *fp = fopen(ctx.fnamelog, "w");
        if (fp == NULL)
        {
            FUNC_RETURN_FAILURE("Cannot create file %s", ctx.fnamelog);
        }
        fclose(fp);
    }

    for (long iter = 0; iter < NBiter; iter++)
    {
        DEBUG_TRACEPOINT("iter %ld/%ld %lf %lf", iter, NBiter, ctx.paramref, ctx.range);
        FUNC_CHECK_RETURN(fpmtransm_scan_iter(&ctx, iter));
    }

    piaacmcopticaldesign.fpmaskamptransm = ctx.parambest;
    FUNC_CHECK_RETURN(init_piaacmcopticaldesign(
        ctx.fpmradld,
        ctx.centobs0,
        ctx.centobs1,
        INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT,
        NULL));
    FUNC_CHECK_RETURN(PIAACMCsimul_savepiaacmcconf(piaacmcparams.piaacmcconfdir));

    piaacmcparams.FORCE_CREATE_fpmza = 0;

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
