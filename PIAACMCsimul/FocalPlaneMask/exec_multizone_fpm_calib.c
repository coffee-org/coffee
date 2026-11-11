/**
 * @file    exec_multizone_fpm_calib.c
 * @brief   Setup and calibrate response of multizone focal plane mask (Mode 11)
 */

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul.h"
#include "PIAACMCsimul_computePSF.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"
#include "init_piaacmcopticaldesign.h"
#include "init_piaacmcopticalsystem.h"

#include "FocalPlaneMask/exec_multizone_fpm_calib.h"
#include "FocalPlaneMask/exec_multizone_fpm_calib_internal.h"
#include "PIAAshape/makePIAAshapes.h"

/**
 * fpm_calib_load_cli_vars - Read CLI optical geometry overrides into calibration context
 * @ctx: Calibration context initialized with default and CLI-overridden parameters.
 */
static void fpm_calib_load_cli_vars(struct fpm_calib_ctx *ctx)
{
    variableID IDv;

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
 * fpm_calib_print_summary - Print active Mode 11 optical design and bandpass summary
 */
static void fpm_calib_print_summary(void)
{
    printf("piaacmcconfdir     : %s\n", piaacmcparams.piaacmcconfdir);
    printf("SCORINGMASKTYPE    : %d\n", piaacmcparams.SCORINGMASKTYPE);
    printf("PIAACMC_FPMsectors : %d\n", piaacmcparams.PIAACMC_FPMsectors);
    printf("lamda              : %ld nm\n", (long) (1.0e9 * piaacmcopticaldesign.lambda + 0.1));
    printf("lamdaB             : %ld \n", (long) (1.0 * piaacmcopticaldesign.lambdaB + 0.1));
    printf("piaacmcopticaldesign.NBrings : %ld\n", piaacmcopticaldesign.NBrings);
    printf("mask rad           : %ld\n", (long) (100.0 * piaacmcparams.PIAACMC_MASKRADLD + 0.1));
    printf("computePSF_ResolvedTarget : %d\n", piaacmcparams.computePSF_ResolvedTarget);
    printf("computePSF_ResolvedTarget_mode : %d\n", piaacmcparams.computePSF_ResolvedTarget_mode);
    printf("piaacmcopticaldesign.fpmmaterial_name : %s\n", piaacmcopticaldesign.fpmmaterial_name);
    printf("piaacmcopticaldesign.nblambda         : %d\n", piaacmcopticaldesign.nblambda);
}

/**
 * fpm_calib_init_mp_mode - Configure multi-process role flags and reinitialize physical FPM
 * @ctx: Calibration context holding pupil and mask geometry parameters.
 *
 * Return: RETURN_SUCCESS on design initialization, error code otherwise.
 */
static errno_t fpm_calib_init_mp_mode(const struct fpm_calib_ctx *ctx)
{
    variableID IDv;

    piaacmcparams.PIAACMC_FPMresp_mp = 1;
    if ((IDv = variable_ID("PIAACMC_FPMresp_mp")) != -1)
    {
        piaacmcparams.PIAACMC_FPMresp_mp = (long) coffee_dcvar_f(IDv) + 0.01;
    }
    printf("PIAACMC_FPMresp_mp = %ld\n", piaacmcparams.PIAACMC_FPMresp_mp);

    printf("------------------------------------- STEP02\n");
    printf(
        "piaacmcopticaldesign.focmNBzone  =  %ld   (%ld)\n",
        piaacmcopticaldesign.focmNBzone,
        piaacmcopticaldesign.NBrings);
    fflush(stdout);

    piaacmcparams.PIAACMC_FPMresp_thread = 0;
    if ((IDv = variable_ID("PIAACMC_FPMresp_thread")) != -1)
    {
        piaacmcparams.PIAACMC_FPMresp_thread = (long) coffee_dcvar_f(IDv) + 0.01;
    }
    printf("PIAACMC_FPMresp_thread = %ld\n", piaacmcparams.PIAACMC_FPMresp_thread);

    if (piaacmcparams.PIAACMC_FPMresp_mp == 1
        || piaacmcparams.PIAACMC_FPMresp_thread > piaacmcparams.PIAACMC_FPMresp_mp - 1)
    {
        piaacmcparams.FORCE_CREATE_fpmzmap = 1;
        piaacmcparams.FORCE_CREATE_fpmzt   = 1;
        piaacmcparams.FORCE_CREATE_fpmza   = 1;
    }
    else
    {
        printf("NO initOK file created\n");
        piaacmcparams.FORCE_CREATE_fpmzmap = 0;
        piaacmcparams.FORCE_CREATE_fpmzt   = 0;
        piaacmcparams.FORCE_CREATE_fpmza   = 0;
    }

    uint64_t initflag = INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT
                        | INIT_PIAACMCOPTICALDESIGN_MODE__READCONF
                        | INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF
                        | INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL;
    return init_piaacmcopticaldesign(ctx->fpmradld, ctx->centobs0, ctx->centobs1, initflag, NULL);
}

/**
 * fpm_calib_resolve_mp_files - Resolve combined or worker FPMresp FITS paths and probe cache
 * @ctx: Calibration context updated with output paths, zone stride, and cached image IDs.
 *
 * Return: RETURN_SUCCESS on path formatting, error code otherwise.
 */
static errno_t fpm_calib_resolve_mp_files(struct fpm_calib_ctx *ctx)
{
    char fname[STRINGMAXLEN_FULLFILENAME];

    PIAACMCsimul_update_fnamedescr_conf();

    if (piaacmcparams.PIAACMC_FPMresp_mp == 1
        || piaacmcparams.PIAACMC_FPMresp_thread > piaacmcparams.PIAACMC_FPMresp_mp - 1)
    {
        SNPRINTF_CHECK(
            fname,
            STRINGMAXLEN_FULLFILENAME,
            "%s/FPMresp%d.%s.fits",
            piaacmcparams.piaacmcconfdir,
            piaacmcparams.SCORINGMASKTYPE,
            piaacmcparams.fnamedescr_conf);
        SNPRINTF_CHECK(ctx->fname1, STRINGMAXLEN_FULLFILENAME, "%s.tmp", fname);
        SNPRINTF_CHECK(ctx->fname2, STRINGMAXLEN_FULLFILENAME, "%s", fname);

        ctx->mzoffset = 0;
        ctx->mzstep   = 1;

        load_fits(fname, "FPMresp", 1, &ctx->ID_FPMresp);
        ctx->IDcomb = ctx->ID_FPMresp;
        return RETURN_SUCCESS;
    }

    char fnamecomb[STRINGMAXLEN_FULLFILENAME];
    SNPRINTF_CHECK(
        fnamecomb,
        STRINGMAXLEN_FULLFILENAME,
        "%s/FPMresp%d.%s.fits",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.SCORINGMASKTYPE,
        piaacmcparams.fnamedescr_conf);

    PIAACMCsimul_update_fnamedescr_conf();
    SNPRINTF_CHECK(
        fname,
        STRINGMAXLEN_FULLFILENAME,
        "%s/FPMresp%d.%s.mp%02ld_thread%02ld.fits",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.SCORINGMASKTYPE,
        piaacmcparams.fnamedescr_conf,
        piaacmcparams.PIAACMC_FPMresp_mp,
        piaacmcparams.PIAACMC_FPMresp_thread);

    SNPRINTF_CHECK(ctx->fname1, STRINGMAXLEN_FULLFILENAME, "%s.tmp", fname);
    SNPRINTF_CHECK(ctx->fname2, STRINGMAXLEN_FULLFILENAME, "%s", fname);
    ctx->mzoffset = piaacmcparams.PIAACMC_FPMresp_thread;
    ctx->mzstep   = piaacmcparams.PIAACMC_FPMresp_mp;

    load_fits(fname, "FPMresp", LOADFITS_ERRMODE_WARNING, &ctx->ID_FPMresp);
    load_fits(fnamecomb, "FPMresp", LOADFITS_ERRMODE_WARNING, &ctx->IDcomb);

    return RETURN_SUCCESS;
}

/**
 * fpm_calib_create_resp_cube - Propagate unobstructed PSF and initialize FPMresp zone 0
 * @ctx:          Active calibration context.
 * @outIDfpmresp: Output pointer receiving allocated FPMresp 3D double image identifier.
 *
 * Return: RETURN_SUCCESS on PSF propagation and cube allocation, error code otherwise.
 */
static errno_t fpm_calib_create_resp_cube(
    const struct fpm_calib_ctx *ctx,
    imageID                    *outIDfpmresp)
{
    piaacmcopticalsystem.FOCMASKarray[0].mode = 1;
    piaacmcopticaldesign.fpmRad =
        0.5 * (piaacmcparams.LAMBDASTART + piaacmcparams.LAMBDAEND) * piaacmcopticaldesign.Fratio
        * piaacmcparams.PIAACMC_MASKRADLD;

    uint64_t initflag = INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT
                        | INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL;
    FUNC_CHECK_RETURN(
        init_piaacmcopticaldesign(ctx->fpmradld, ctx->centobs0, ctx->centobs1, initflag, NULL));
    FUNC_CHECK_RETURN(makePIAAshapes());
    FUNC_CHECK_RETURN(init_piaacmcopticalsystem(0.0, 0.0));
    FUNC_CHECK_RETURN(makePIAAshapes());

    piaacmcparams.focmMode = dcimg[piaacmcopticaldesign.zonezID].md[0].size[0] + 10;
    piaacmcopticalsystem.FOCMASKarray[0].mode = 1;

    FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(
        0.0,
        0.0,
        0,
        piaacmcopticalsystem.NBelem,
        0,
        piaacmcparams.computePSF_ResolvedTarget,
        piaacmcparams.computePSF_ResolvedTarget_mode,
        0,
        NULL));

    imageID  IDimvect = image_ID("imvect", dcimg, dcnimg);
    uint32_t neval    = dcimg[IDimvect].md[0].size[0];
    uint32_t nzones1  = (uint32_t) (piaacmcopticaldesign.focmNBzone + 1);
    imageID  IDfpmresp;

    FUNC_CHECK_RETURN(create_3Dimage_ID_double(
        "FPMresp",
        neval,
        nzones1,
        piaacmcopticaldesign.nblambda,
        &IDfpmresp));

    for (int k = 0; k < piaacmcopticaldesign.nblambda; k++)
    {
        size_t kbase = (size_t) k * nzones1 * neval;
        for (uint32_t ii = 0; ii < neval; ii++)
        {
            dcimg[IDfpmresp].array.D[kbase + ii] = dcimg[IDimvect].array.F[(size_t) k * neval + ii];
        }
    }

    *outIDfpmresp = IDfpmresp;
    return RETURN_SUCCESS;
}

/**
 * exec_multizone_fpm_calib - Calibrate polychromatic multi-zone FPM response cube (Mode 11)
 *
 * Computes or assembles the 3D complex electric-field impulse-response matrix (FPMresp)
 * across evaluation dark-hole points, individual FPM zones (plus exterior light in zone 0),
 * and bandpass wavelengths, supporting both single-process and multi-worker execution.
 *
 * Return: RETURN_SUCCESS on calibration completion, error code otherwise.
 */
errno_t exec_multizone_fpm_calib(void)
{
    DEBUG_TRACE_FSTART();

    struct fpm_calib_ctx ctx;
    memset(&ctx, 0, sizeof(ctx));

    printf("=================================== mode 011 ===================================\n");
    fpm_calib_load_cli_vars(&ctx);

    {
        uint64_t initflag = INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT
                            | INIT_PIAACMCOPTICALDESIGN_MODE__READCONF
                            | INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF
                            | INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL;
        FUNC_CHECK_RETURN(
            init_piaacmcopticaldesign(ctx.fpmradld, ctx.centobs0, ctx.centobs1, initflag, NULL));
    }
    fpm_calib_print_summary();

    {
        char fname[STRINGMAXLEN_FULLFILENAME];

        PIAACMCsimul_update_fnamedescr_conf();
        SNPRINTF_CHECK(
            fname,
            STRINGMAXLEN_FULLFILENAME,
            "%s/FPMresp%d.%s.fits",
            piaacmcparams.piaacmcconfdir,
            piaacmcparams.SCORINGMASKTYPE,
            piaacmcparams.fnamedescr_conf);
        load_fits(fname, "FPMresp", 1, &ctx.ID_FPMresp);
    }

    if (ctx.ID_FPMresp == -1)
    {
        FUNC_CHECK_RETURN(fpm_calib_init_mp_mode(&ctx));
        FUNC_CHECK_RETURN(fpm_calib_resolve_mp_files(&ctx));

        if (ctx.IDcomb == -1 && ctx.ID_FPMresp == -1)
        {
            imageID IDfpmresp = -1;

            FUNC_CHECK_RETURN(fpm_calib_create_resp_cube(&ctx, &IDfpmresp));
            if (piaacmcparams.PIAACMC_FPMresp_thread > piaacmcparams.PIAACMC_FPMresp_mp - 1)
            {
                FUNC_CHECK_RETURN(fpm_calib_combine_threads(IDfpmresp));
            }
            else
            {
                FUNC_CHECK_RETURN(fpm_calib_compute_worker_zones(&ctx, IDfpmresp));
            }
        }
        else
        {
            printf("File FPMresp exists\n");
        }
    }

    piaacmcparams.focmMode = -1;

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
