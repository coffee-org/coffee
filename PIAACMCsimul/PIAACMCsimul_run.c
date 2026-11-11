/**
 * @file    PIAACMCsimul_run.c
 * @brief   CLI entry point and Mode 13 multi-pass FPM optimization driver
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/time.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "OptSystProp/OptSystProp.h"
#include "OpticsMaterials/OpticsMaterials.h"

#include "PIAACMCsimul.h"
#include "PIAACMCsimul_exec.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"
#include "PIAACMCsimul_run.h"
#include "PIAACMCsimul_run_internal.h"

static char    *confindex;
static int64_t *mode;

static CLICMDARGDEF piaacmc_run_farg[] =
{
    {
        CLIARG_STR,
        ".confindex",
        "configuration index",
        "000",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &confindex,
        NULL
    },
    {
        CLIARG_INT64,
        ".mode",
        "mode",
        "0",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &mode,
        NULL
    }
};

MILK_RUNNABLE_RT_DEFINE(piaacmc_run);

/**
 * cli_piaacmc_run - CLI entry wrapper invoking PIAACMCsimul_run.
 *
 * Return: errno_t status from PIAACMCsimul_run.
 */
static errno_t cli_piaacmc_run(void)
{
    DEBUG_TRACE_FSTART();
    FUNC_CHECK_RETURN(PIAACMCsimul_run(confindex, (long) *mode));
    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

const MILK_RUNNABLE runnable_piaacmc_run =
{
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "run", .description = "Simulate PIAACMC" },
    .srcfile  = __FILE__,
    .farg     = piaacmc_run_farg,
    .nbarg    = (int) (sizeof(piaacmc_run_farg) / sizeof(piaacmc_run_farg[0])),
    .cli_func = cli_piaacmc_run,
    .rt       = &piaacmc_run_rt,
};

/**
 * mode13_sync_best_fits - Update in-memory fpmbestsol and copy FPM FITS products to .best.fits.
 * @ctx: Mutable Mode 13 outer loop context
 */
static void mode13_sync_best_fits(struct mode13_loop_ctx *ctx)
{
    char fname1[STRINGMAXLEN_FULLFILENAME];
    char fnamebestsol[STRINGMAXLEN_FULLFILENAME];

    PIAACMCsimul_update_fnamedescr();
    SNPRINTF_CHECK(
        fnamebestsol,
        STRINGMAXLEN_FULLFILENAME,
        "%s/fpm_zonez.%s.best.fits",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.fnamedescr);

    if (ctx->IDbestsol == -1)
    {
        load_fits(fnamebestsol, "fpmbestsol", LOADFITS_ERRMODE_IGNORE, &ctx->IDbestsol);
    }
    else
    {
        imageID IDbestsoltmp = -1;
        load_fits(fnamebestsol, "fpmbestsoltmp", LOADFITS_ERRMODE_IGNORE, &IDbestsoltmp);
        if (IDbestsoltmp != -1)
        {
            for (uint32_t k = 0; k < dcimg[ctx->IDbestsol].md[0].size[0]; k++)
            {
                dcimg[ctx->IDbestsol].array.D[k] = dcimg[IDbestsoltmp].array.D[k];
            }
            delete_image_ID("fpmbestsoltmp", DELETE_IMAGE_ERRMODE_WARNING);
        }
    }

    PIAACMCsimul_update_fnamedescr();
    SNPRINTF_CHECK(
        fname1,
        STRINGMAXLEN_FULLFILENAME,
        "%s/fpm_zonez.%s.fits",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.fnamedescr);
    SNPRINTF_CHECK(
        fnamebestsol,
        STRINGMAXLEN_FULLFILENAME,
        "%s/fpm_zonez.%s.best.fits",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.fnamedescr);

    EXECUTE_SYSTEM_COMMAND("cp %s %s", fname1, fnamebestsol);
    EXECUTE_SYSTEM_COMMAND("echo \"cp %s %s\" > cmdlogtest.txt", fname1, fnamebestsol);

    SNPRINTF_CHECK(
        fname1,
        STRINGMAXLEN_FULLFILENAME,
        "%s/fpm_sagmapHR.%s.fits",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.fnamedescr);
    SNPRINTF_CHECK(
        fnamebestsol,
        STRINGMAXLEN_FULLFILENAME,
        "%s/fpm_sagmapHR.%s.best.fits",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.fnamedescr);

    EXECUTE_SYSTEM_COMMAND("cp %s %s", fname1, fnamebestsol);
    EXECUTE_SYSTEM_COMMAND("echo \"cp %s %s\" > cmdlogtest1.txt", fname1, fnamebestsol);
}

/**
 * mode13_write_bestval_file - Record updated best contrast merit and touch newbestsol.txt.
 * @fnamebestval: Path to mode13.<descr>.bestval.txt
 * @bestval:      Updated best objective merit value
 */
static void mode13_write_bestval_file(
    const char *fnamebestval,
    double      bestval)
{
    FILE *fp = fopen(fnamebestval, "w");
    if (fp != NULL)
    {
        fprintf(
            fp,
            "%30g %d %04ld %02ld %03ld %5.2f %02d %d %s %02d %07.3f\n",
            bestval,
            piaacmcparams.PIAACMC_FPMsectors,
            (long) (1.0e9 * piaacmcopticaldesign.lambda + 0.1),
            (long) (1.0 * piaacmcopticaldesign.lambdaB + 0.1),
            piaacmcopticaldesign.NBrings,
            piaacmcparams.PIAACMC_MASKRADLD,
            piaacmcparams.computePSF_ResolvedTarget,
            piaacmcparams.computePSF_ResolvedTarget_mode,
            piaacmcopticaldesign.fpmmaterial_name,
            piaacmcopticaldesign.nblambda,
            piaacmcopticaldesign.fpmsagreg_coeff);
        fclose(fp);
    }

    EXECUTE_SYSTEM_COMMAND("touch %s/newbestsol.txt", piaacmcparams.piaacmcconfdir);
}

/**
 * mode13_log_opt_entry - Append per-iteration optimization summary line to mode13.*.opt.txt.
 * @ctx: Mode 13 outer loop context
 * @bOK: Non-zero if this iteration achieved a new global best
 */
static void mode13_log_opt_entry(
    const struct mode13_loop_ctx *ctx,
    int                           bOK)
{
    char fname[STRINGMAXLEN_FULLFILENAME];

    PIAACMCsimul_update_fnamedescr();
    SNPRINTF_CHECK(
        fname,
        STRINGMAXLEN_FULLFILENAME,
        "%s/mode13.%s.opt.txt",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.fnamedescr);

    printf("Add current solution (possibly not best) to file: %s\n", fname);

    FILE *fp = fopen(fname, "a");
    if (fp == NULL)
    {
        return;
    }

    fprintf(
        fp,
        "%10d %20.5g   %16.5g -> %16.5g   (%16.5g) %d "
        "%5ld/%5ld [%12g %2d %12g %12g  %12g]%s\n",
        ctx->loopiter,
        piaacmcparams.MODampl,
        piaacmcparams.PIAACMCSIMUL_VALREF,
        piaacmcparams.PIAACMCSIMUL_VAL,
        ctx->bestval,
        ctx->zeroST,
        ctx->cnt0,
        ctx->cnt00,
        piaacmcparams.CnormFactor,
        piaacmcopticaldesign.nblambda,
        piaacmcopticalsystem.flux[0],
        piaacmcparams.SCORINGTOTAL,
        piaacmcparams.PIAACMCSIMUL_VAL0,
        (bOK == 1) ? " BEST" : "");
    fclose(fp);
}

/**
 * mode13_should_continue - Check stoploop13.txt and wall-clock searchtime budget.
 * @ctx: Mode 13 outer loop context
 *
 * Return: 1 to continue outer search loop, 0 to terminate.
 */
static int mode13_should_continue(const struct mode13_loop_ctx *ctx)
{
    int loopOK = 1;

    if (file_exists(ctx->stopfile) == 1)
    {
        printf("FILE \"%s\" found -> stopping\n", ctx->stopfile);
        loopOK = 0;
        EXECUTE_SYSTEM_COMMAND("touch stop.stopfile.txt");
        EXECUTE_SYSTEM_COMMAND("rm %s", ctx->stopfile);
    }
    else
    {
        printf("File \"%s\" not found -> continuing\n", ctx->stopfile);
    }

    struct timeval end;
    gettimeofday(&end, NULL);
    printf("start: %ld secs, %ld usecs\n", (long) ctx->start.tv_sec, (long) ctx->start.tv_usec);
    printf("end: %ld secs, %ld usecs\n", (long) end.tv_sec, (long) end.tv_usec);

    long secs_used   = end.tv_sec - ctx->start.tv_sec;
    long micros_used = (secs_used * 1000000) + end.tv_usec - ctx->start.tv_usec;

    FILE *fp = fopen("timeused.txt", "w");
    if (fp != NULL)
    {
        fprintf(fp, "# Time used vs. time limit\n");
        fprintf(fp, "# Time limit can be written in file searchtime.txt\n\n");
        fprintf(fp, "%12.3f    %12.3f\n", 1.0e-6 * micros_used, ctx->searchtime);
        fclose(fp);
    }

    if (micros_used > 1000000.0 * ctx->searchtime)
    {
        loopOK = 0;
        EXECUTE_SYSTEM_COMMAND("touch stop.time_elapsed.ttxt");
    }

    printf("End of loop\n");
    return loopOK;
}

/**
 * mode13_finalize_summary - Print final Mode 13 configuration summary and restore best FPM sags.
 */
static void mode13_finalize_summary(void)
{
    printf("piaacmcconfdir              : %s\n", piaacmcparams.piaacmcconfdir);
    printf("computePSF_ResolvedTarget   : %d\n", piaacmcparams.computePSF_ResolvedTarget);
    printf("computePSF_ResolvedTarget_mode   : %d\n", piaacmcparams.computePSF_ResolvedTarget_mode);
    printf("PIAACMC_FPMsectors   : %d\n", piaacmcparams.PIAACMC_FPMsectors);
    printf(
        "(long) (10.0*PIAACMC_MASKRADLD+0.1)   : %ld\n",
        (long) (10.0 * piaacmcparams.PIAACMC_MASKRADLD + 0.1));
    printf("piaacmcopticaldesign.NBrings   : %ld\n", piaacmcopticaldesign.NBrings);
    printf("piaacmcopticaldesign.nblambda   : %d\n", piaacmcopticaldesign.nblambda);
    fflush(stdout);

    char fname[STRINGMAXLEN_FULLFILENAME];
    char fnamebestsol[STRINGMAXLEN_FULLFILENAME];

    PIAACMCsimul_update_fnamedescr();
    SNPRINTF_CHECK(
        fname,
        STRINGMAXLEN_FULLFILENAME,
        "%s/fpm_zonez.%s.fits",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.fnamedescr);
    SNPRINTF_CHECK(
        fnamebestsol,
        STRINGMAXLEN_FULLFILENAME,
        "%s/fpm_zonez.%s.best.fits",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.fnamedescr);

    EXECUTE_SYSTEM_COMMAND("cp %s %s", fnamebestsol, fname);
}

/**
 * PIAACMCsimul_run_mode13 - Execute randomized multi-start FPM zone sag optimization (Mode 13).
 * @confindex: Configuration directory name/index
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t PIAACMCsimul_run_mode13(const char *confindex)
{
    struct mode13_loop_ctx ctx;

    ctx.bestval    = 1.0;
    ctx.loopin     = 0;
    ctx.loopiter   = 0;
    ctx.IDbestsol  = -1;
    ctx.sag0       = 1.0;
    ctx.searchtime = PIAACMCsimul_run_read_searchtime();

    int loopOK = (ctx.searchtime < 0.1) ? 0 : 1;
    printf("loopOK = %d\n", loopOK);
    fflush(stdout);

    gettimeofday(&ctx.start, NULL);
    SNPRINTF_CHECK(
        ctx.stopfile,
        STRINGMAXLEN_FULLFILENAME,
        "%s/stoploop13.txt",
        piaacmcparams.piaacmcconfdir);

    while ((loopOK == 1) && (ctx.loopiter < 1000000))
    {
        printf("LOOP start\n");
        fflush(stdout);

        ctx.searchtime = PIAACMCsimul_run_read_searchtime();
        ctx.loopin     = 1;
        PIAACMCsimul_run_prep_mode13_zones(&ctx);

        printf("Execute optimization\n");
        FUNC_CHECK_RETURN(PIAACMCsimul_exec(confindex, 13));

        double pha0 = OpticsMaterials_pha_lambda(
            piaacmcopticaldesign.fpmmaterial_code, ctx.sag0, piaacmcopticaldesign.lambda);
        printf("%g m  -> %g rad\n", ctx.sag0, pha0);
        fflush(stdout);
        ctx.sag0 = ctx.sag0 / (pha0 / 2.0 / M_PI);
        printf(
            "======================= sag0 = %g m  -> %g rad\n",
            ctx.sag0,
            (double) OpticsMaterials_pha_lambda(
                piaacmcopticaldesign.fpmmaterial_code, ctx.sag0, piaacmcopticaldesign.lambda));

        if (ctx.IDbestsol == -1)
        {
            char fnamebestsol[STRINGMAXLEN_FULLFILENAME];
            PIAACMCsimul_update_fnamedescr();
            SNPRINTF_CHECK(
                fnamebestsol,
                STRINGMAXLEN_FULLFILENAME,
                "%s/fpm_zonez.%s.best.fits",
                piaacmcparams.piaacmcconfdir,
                piaacmcparams.fnamedescr);
            printf("LOADING \"%s\"...\n", fnamebestsol);
            fflush(stdout);
            FUNC_CHECK_RETURN(
                load_fits(fnamebestsol, "fpmbestsol", LOADFITS_ERRMODE_IGNORE, &ctx.IDbestsol));
        }

        char fnamebestval[STRINGMAXLEN_FULLFILENAME];
        PIAACMCsimul_update_fnamedescr();
        SNPRINTF_CHECK(
            fnamebestval,
            STRINGMAXLEN_FULLFILENAME,
            "%s/mode13.%s.bestval.txt",
            piaacmcparams.piaacmcconfdir,
            piaacmcparams.fnamedescr);

        if (ctx.loopiter == 0)
        {
            printf("READING FILE \"%s\"\n", fnamebestval);
            fflush(stdout);
            FILE *fp = fopen(fnamebestval, "r");
            if (fp != NULL)
            {
                int r = fscanf(fp, "%lf", &ctx.bestval);
                (void) r;
                fclose(fp);
            }
        }

        printf("\n\n\n\n======= val = %g [%g]\n", piaacmcparams.PIAACMCSIMUL_VAL, ctx.bestval);
        fflush(stdout);

        int bOK = 0;
        if (piaacmcparams.PIAACMCSIMUL_VAL < ctx.bestval)
        {
            bOK         = 1;
            ctx.bestval = piaacmcparams.PIAACMCSIMUL_VAL;
            printf("========================================= ");
            printf("SAVING BEST MASK SOLUTION -> fpm_zonez.[name].best.fits\n");
            mode13_sync_best_fits(&ctx);
            mode13_write_bestval_file(fnamebestval, ctx.bestval);
        }

        mode13_log_opt_entry(&ctx, bOK);
        ctx.loopiter++;
        loopOK = mode13_should_continue(&ctx);
    } // while loopOK

    if (ctx.loopin == 1)
    {
        mode13_finalize_summary();
    }

    return RETURN_SUCCESS;
}

/**
 * PIAACMCsimul_run - Primary entry point for PIAACMC design and evaluation modes.
 * @confindex: Configuration directory name/index
 * @mode:      Operation mode index
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
errno_t PIAACMCsimul_run(
    const char *confindex,
    long        mode)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %s mode %ld", confindex, mode);

#ifdef PIAASIMUL_LOGFUNC0
    PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
#endif

    FUNC_CHECK_RETURN(PIAACMCsimul_run_setparam_variables(confindex, mode));

    if (mode == 13)
    {
        FUNC_CHECK_RETURN(PIAACMCsimul_run_mode13(confindex));
    }
    else
    {
        DEBUG_TRACEPOINT("running exec function");
        FUNC_CHECK_RETURN(PIAACMCsimul_exec(confindex, mode));
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
