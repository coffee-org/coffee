/**
 * @file    PIAACMCsimul_linopt.c
 * @brief   Reference initialization, Jacobian assembly, and outer SVD optimization loop
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_tools/COREMOD_tools.h"

#include "PIAACMCsimul.h"
#include "PIAACMCsimul_achromFPMsol_eval_zonezderivative.h"
#include "PIAACMCsimul_computePSF.h"
#include "PIAACMCsimul_linopt.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"
#include "PIAACMCsimul_regularization.h"
#include "PIAAshape/makePIAAshapes.h"

/**
 * linopt_stage_reference_files - Snapshot initial state into `<confdir>_linopt` and load refs.
 *
 * Return: RETURN_SUCCESS or FITS/config I/O error status.
 */
static errno_t linopt_stage_reference_files(void)
{
    DEBUG_TRACE_FSTART();

    char dirname[STRINGMAXLEN_DIRNAME];

    SNPRINTF_CHECK(dirname, STRINGMAXLEN_DIRNAME, "%s_linopt", piaacmcparams.piaacmcconfdir);
    FUNC_CHECK_RETURN(PIAACMCsimul_savepiaacmcconf(dirname));

    EXECUTE_SYSTEM_COMMAND("rsync -au --progress %s/* ./%s/",
                           dirname,
                           piaacmcparams.piaacmcconfdir);
    EXECUTE_SYSTEM_COMMAND("cd \"%s\" && { test -f piaacmcparams.ref.conf || "
                           "cp piaacmcparams.conf piaacmcparams.ref.conf; }",
                           dirname);

    if (piaacmcopticaldesign.PIAAmode == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];

        EXECUTE_SYSTEM_COMMAND("cd \"%s\" && { test -f piaa0Cmodes.ref.fits || "
                               "cp piaa0Cmodes.fits piaa0Cmodes.ref.fits; }",
                               dirname);
        EXECUTE_SYSTEM_COMMAND("cd \"%s\" && { test -f piaa0Fmodes.ref.fits || "
                               "cp piaa0Fmodes.fits piaa0Fmodes.ref.fits; }",
                               dirname);
        EXECUTE_SYSTEM_COMMAND("cd \"%s\" && { test -f piaa1Cmodes.ref.fits || "
                               "cp piaa1Cmodes.fits piaa1Cmodes.ref.fits; }",
                               dirname);
        EXECUTE_SYSTEM_COMMAND("cd \"%s\" && { test -f piaa1Fmodes.ref.fits || "
                               "cp piaa1Fmodes.fits piaa1Fmodes.ref.fits; }",
                               dirname);

        WRITE_FULLFILENAME(fname, "%s/piaa0Cmodes.ref.fits", dirname);
        FUNC_CHECK_RETURN(load_fits(fname, "piaa0Cmref", 1, NULL));
        WRITE_FULLFILENAME(fname, "%s/piaa1Cmodes.ref.fits", dirname);
        FUNC_CHECK_RETURN(load_fits(fname, "piaa1Cmref", 1, NULL));
        WRITE_FULLFILENAME(fname, "%s/piaa0Fmodes.ref.fits", dirname);
        FUNC_CHECK_RETURN(load_fits(fname, "piaa0Fmref", 1, NULL));
        WRITE_FULLFILENAME(fname, "%s/piaa1Fmodes.ref.fits", dirname);
        FUNC_CHECK_RETURN(load_fits(fname, "piaa1Fmref", 1, NULL));
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * linopt_init_reference_state - Compute initial reference contrast and build `vecDHref1D`.
 * @ctx: Linear optimization context to initialize
 *
 * Return: RETURN_SUCCESS or propagation/memory error status.
 */
static errno_t linopt_init_reference_state(PIAACMC_linopt_ctx *ctx)
{
    DEBUG_TRACE_FSTART();

    double valref = 0.0;

    dcimg[ctx->IDstatus].array.UI16[0] = 5;
    FUNC_CHECK_RETURN(makePIAAshapes());

    piaacmcopticalsystem.FOCMASKarray[0].mode = 1;
    FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(0.0,
                                              0.0,
                                              0,
                                              piaacmcopticalsystem.NBelem,
                                              0,
                                              piaacmcparams.computePSF_ResolvedTarget,
                                              piaacmcparams.computePSF_ResolvedTarget_mode,
                                              0,
                                              &valref));

    FUNC_CHECK_RETURN(linopt_stage_reference_files());

    if (piaacmcparams.linopt_REGPIAASHAPES == 1)
    {
        double val0 = PIAACMCsimul_regularization_PIAAshapes_value();
        printf("VALREF = %g + %g -> %g\n", valref, val0, valref + val0);
        valref += val0;
    }
    if (piaacmcparams.linopt_REGFPMSAG == 1)
    {
        valref += PIAACMCsimul_regularization_fpmsag_value();
    }

    dcimg[ctx->IDstatus].array.UI16[0] = 6;
    printf("================================ Reference = %g\n", valref);

    chname_image_ID("imvect", "vecDHref");
    imageID ID = image_ID("vecDHref", dcimg, dcnimg);

    dcimg[ctx->IDstatus].array.UI16[0] = 7;
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/vecDMref.fits", piaacmcparams.piaacmcconfdir);
        save_fits("vecDHref", fname);
    }

    long size1Dvec0 = dcimg[ID].md[0].nelement;
    long size1Dvec  = size1Dvec0;
    if (piaacmcparams.linopt_REGPIAASHAPES == 1)
    {
        size1Dvec += dcimg[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0];
        size1Dvec += dcimg[piaacmcopticaldesign.piaa1CmodesID].md[0].size[0];
        size1Dvec += dcimg[piaacmcopticaldesign.piaa0FmodesID].md[0].size[0];
        size1Dvec += dcimg[piaacmcopticaldesign.piaa1FmodesID].md[0].size[0];
    }
    if (piaacmcparams.linopt_REGFPMSAG == 1)
    {
        size1Dvec += dcimg[piaacmcopticaldesign.zonezID].md[0].size[0];
    }

    imageID IDm;
    imageID ID1Dref;
    FUNC_CHECK_RETURN(create_2Dimage_ID("DHmask", size1Dvec, 1, &IDm));
    FUNC_CHECK_RETURN(create_2Dimage_ID("vecDHref1D", size1Dvec, 1, &ID1Dref));

    ID = image_ID("vecDHref", dcimg, dcnimg);
    /* Index `ii` advances across dark-hole field entries and optional regularization tails */
    long ii;
    for (ii = 0; ii < (long) dcimg[ID].md[0].nelement; ii++)
    {
        dcimg[ID1Dref].array.F[ii] = dcimg[ID].array.F[ii];
        dcimg[IDm].array.F[ii]     = 1.0;
    }

    dcimg[ctx->IDstatus].array.UI16[0] = 8;
    if (piaacmcparams.linopt_REGPIAASHAPES == 1)
    {
        ii = PIAACMCsimul_regularization_PIAAshapes_add1Dvector(ID1Dref, ii);
    }
    if (piaacmcparams.linopt_REGFPMSAG == 1)
    {
        ii = PIAACMCsimul_regularization_fpmsag_add1Dvector(ID1Dref, IDm, ii);
    }

    FUNC_CHECK_RETURN(delete_image_ID("vecDHref", DELETE_IMAGE_ERRMODE_WARNING));

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/linoptval.txt", piaacmcparams.piaacmcconfdir);
        FILE *fp = fopen(fname, "w");
        fclose(fp);
    }

    ctx->valref     = valref;
    ctx->size1Dvec  = size1Dvec;
    ctx->size1Dvec0 = size1Dvec0;

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * linopt_compute_fast_fpm_jacobian - Assemble analytical FPM zone-sag Jacobian `DHmodes`.
 * @ctx:       Linear optimization context
 * @IDmodes:   Output 3D Jacobian cube image ID (`DHmodes`)
 * @IDmodes2D: Diagnostic 2D Jacobian image ID (`DHmodes2D`)
 *
 * Return: RETURN_SUCCESS or derivative/FITS error status.
 */
static errno_t linopt_compute_fast_fpm_jacobian(
    const PIAACMC_linopt_ctx *ctx,
    imageID                   IDmodes,
    imageID                   IDmodes2D)
{
    DEBUG_TRACE_FSTART();

    imageID  IDz     = piaacmcopticaldesign.zonezID;
    uint32_t nbzones = dcimg[IDz].md[0].size[0];

    dcimg[ctx->IDstatus].array.UI16[0]        = 12;
    piaacmcopticalsystem.FOCMASKarray[0].mode = 1;

    printf("Computing %ld derivatives ", (long) nbzones);
    fflush(stdout);

    for (uint32_t mz = 0; mz < nbzones; mz++)
    {
        printf(" %ld", (long) mz);
        FUNC_CHECK_RETURN(
            PIAACMCsimul_achromFPMsol_eval_zonezderivative(mz,
                                                           piaacmcparams.fpmresp_array,
                                                           piaacmcparams.zonez_array,
                                                           piaacmcparams.dphadz_array,
                                                           piaacmcparams.outtmp_array,
                                                           piaacmcparams.vsize,
                                                           nbzones,
                                                           piaacmcopticaldesign.nblambda));
        for (long ii = 0; ii < ctx->size1Dvec0; ii++)
        {
            dcimg[IDmodes].array.F[mz * ctx->size1Dvec + ii] =
                piaacmcparams.outtmp_array[ii] * piaacmcparams.linopt_paramdelta[mz];
        }
    }

    if (piaacmcparams.linopt_REGFPMSAG == 1)
    {
        for (uint32_t mz = 0; mz < nbzones; mz++)
        {
            dcimg[IDmodes].array.F[mz * ctx->size1Dvec + (ctx->size1Dvec0 + mz)] =
                (piaacmcopticaldesign.fpmsagreg_alpha / piaacmcopticaldesign.fpmsagreg_coeff) *
                pow(dcimg[IDz].array.D[mz] / piaacmcopticaldesign.fpmsagreg_coeff,
                    piaacmcopticaldesign.fpmsagreg_alpha - 1.0) *
                piaacmcparams.linopt_paramdelta[mz];
        }
    }

    memcpy(dcimg[IDmodes2D].array.F,
           dcimg[IDmodes].array.F,
           sizeof(float) * ctx->size1Dvec * piaacmcparams.linopt_number_param);
    FUNC_CHECK_RETURN(save_fl_fits("DHmodes2D", "test_DHmodes2D.fits"));

    dcimg[ctx->IDstatus].array.UI16[0] = 13;
    printf("Done computing derivatives (FAST MODE)\n");
    fflush(stdout);

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * linopt_compute_numeric_jacobian - Assemble finite-difference PIAA Jacobian `DHmodes`.
 * @ctx:     Linear optimization context
 * @IDmodes: Output 3D Jacobian cube image ID (`DHmodes`)
 *
 * Return: RETURN_SUCCESS or PSF propagation/FITS error status.
 */
static errno_t linopt_compute_numeric_jacobian(
    const PIAACMC_linopt_ctx *ctx,
    imageID                   IDmodes)
{
    DEBUG_TRACE_FSTART();

    imageID ID1Dref = image_ID("vecDHref1D", dcimg, dcnimg);

    dcimg[ctx->IDstatus].array.UI16[0] = 14;
    for (int i = 0; i < piaacmcparams.linopt_number_param; i++)
    {
        double contrastval = 0.0;

        piaacmcopticalsystem.FOCMASKarray[0].mode = 1;
        if (piaacmcparams.linopt_paramtype[i] == _DATATYPE_FLOAT)
        {
            *(piaacmcparams.linopt_paramvalf[i]) += (float) piaacmcparams.linopt_paramdelta[i];
        }
        else
        {
            *(piaacmcparams.linopt_paramval[i]) += piaacmcparams.linopt_paramdelta[i];
        }

        FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(0.0,
                                                  0.0,
                                                  0,
                                                  piaacmcopticalsystem.NBelem,
                                                  0,
                                                  piaacmcparams.computePSF_ResolvedTarget,
                                                  piaacmcparams.computePSF_ResolvedTarget_mode,
                                                  0,
                                                  &contrastval));

        imageID ID = image_ID("imvect", dcimg, dcnimg);
        {
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(fname, "%s/linoptval.txt", piaacmcparams.piaacmcconfdir);
            FILE *fp = fopen(fname, "a");
            fprintf(fp,
                    "# %5ld/%5ld %5d/%5ld %6ld %6ld    %20.15g %20.15g %20.15g      %20.15g\n",
                    ctx->iter,
                    piaacmcparams.linopt_NBiter,
                    i,
                    piaacmcparams.linopt_number_param,
                    dcimg[ID].md[0].nelement,
                    ctx->size1Dvec,
                    piaacmcparams.linopt_paramdelta[i],
                    contrastval,
                    ctx->valref,
                    ctx->bestval);
            fclose(fp);
        }

        imageID ID1D;
        FUNC_CHECK_RETURN(create_2Dimage_ID("imvect1D", ctx->size1Dvec, 1, &ID1D));
        long ii;
        for (ii = 0; ii < (long) dcimg[ID].md[0].nelement; ii++)
        {
            dcimg[ID1D].array.F[ii] = dcimg[ID].array.F[ii];
        }
        if (piaacmcparams.linopt_REGPIAASHAPES == 1)
        {
            PIAACMCsimul_regularization_PIAAshapes_add1Dvector(ID1D, ii);
        }

        FUNC_CHECK_RETURN(delete_image_ID("imvect", DELETE_IMAGE_ERRMODE_WARNING));

        if (piaacmcparams.linopt_paramtype[i] == _DATATYPE_FLOAT)
        {
            *(piaacmcparams.linopt_paramvalf[i]) -= (float) piaacmcparams.linopt_paramdelta[i];
        }
        else
        {
            *(piaacmcparams.linopt_paramval[i]) -= piaacmcparams.linopt_paramdelta[i];
        }

        for (uint64_t kk = 0; kk < dcimg[ID1D].md[0].nelement; kk++)
        {
            dcimg[IDmodes].array.F[i * dcimg[ID1D].md[0].nelement + kk] =
                dcimg[ID1D].array.F[kk] - dcimg[ID1Dref].array.F[kk];
        }

        imageID IDdiag;
        FUNC_CHECK_RETURN(create_2Dimage_ID("DHmodes2D",
                                            ctx->size1Dvec,
                                            piaacmcparams.linopt_number_param,
                                            &IDdiag));
        for (uint64_t kk = 0; kk < dcimg[IDmodes].md[0].nelement; kk++)
        {
            dcimg[IDdiag].array.F[kk] = dcimg[IDmodes].array.F[kk];
        }

        {
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(fname, "%s/DMmodes.fits", piaacmcparams.piaacmcconfdir);
            FUNC_CHECK_RETURN(save_fits("DHmodes2D", fname));
        }
        FUNC_CHECK_RETURN(delete_image_ID("DHmodes2D", DELETE_IMAGE_ERRMODE_WARNING));
    }

    dcimg[ctx->IDstatus].array.UI16[0] = 15;
    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * linopt_commit_iteration - Apply `optvec`, refresh `vecDHref1D`, and checkpoint `_linopt`.
 * @ctx:             Linear optimization context
 * @IDoptvec:        Best parameter vector image ID from line-search
 * @contrastval_out: Output updated regularized objective value
 *
 * Return: RETURN_SUCCESS or propagation/FITS error status.
 */
static errno_t linopt_commit_iteration(
    const PIAACMC_linopt_ctx *ctx,
    imageID                   IDoptvec,
    double                   *contrastval_out)
{
    DEBUG_TRACE_FSTART();

    for (long i = 0; i < piaacmcparams.linopt_number_param; i++)
    {
        if (piaacmcparams.linopt_paramtype[i] == _DATATYPE_FLOAT)
        {
            *(piaacmcparams.linopt_paramvalf[i]) = dcimg[IDoptvec].array.F[i];
        }
        else
        {
            *(piaacmcparams.linopt_paramval[i]) = (double) dcimg[IDoptvec].array.F[i];
        }
    }

    double contrastval = 0.0;
    FUNC_CHECK_RETURN(PIAACMCsimul_linopt_eval_objective(&contrastval, NULL, NULL, NULL));

    printf("gain: %lf -> val = %20g\n", ctx->bestgain, contrastval);
    dcimg[ctx->IDstatus].array.UI16[0] = 27;

    imageID ID1Dref = image_ID("vecDHref1D", dcimg, dcnimg);
    imageID ID      = image_ID("imvect", dcimg, dcnimg);
    long    ii;
    for (ii = 0; ii < (long) dcimg[ID].md[0].nelement; ii++)
    {
        dcimg[ID1Dref].array.F[ii] = dcimg[ID].array.F[ii];
    }
    dcimg[ctx->IDstatus].array.UI16[0] = 28;

    if (piaacmcparams.linopt_REGPIAASHAPES == 1)
    {
        ii = PIAACMCsimul_regularization_PIAAshapes_add1Dvector(ID1Dref, ii);
    }
    if (piaacmcparams.linopt_REGFPMSAG == 1)
    {
        PIAACMCsimul_regularization_fpmsag_add1Dvector(ID1Dref, -1, ii);
    }

    FUNC_CHECK_RETURN(delete_image_ID("imvect", DELETE_IMAGE_ERRMODE_WARNING));
    dcimg[ctx->IDstatus].array.UI16[0] = 29;

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/linoptval.txt", piaacmcparams.piaacmcconfdir);
        FILE *fp = fopen(fname, "a");
        if (fp == NULL)
        {
            printf("ERROR: cannot open file \"%s\"\n", fname);
            exit(0);
        }
        fprintf(fp, "-> %5ld    %20g <- %20g \n", ctx->iter, contrastval, ctx->valref);
        printf("%5ld %20g %20g \n", ctx->iter, contrastval, ctx->valref);
        fflush(stdout);
        fclose(fp);
    }

    piaacmcparams.PIAACMCSIMUL_VAL    = contrastval;
    piaacmcparams.PIAACMCSIMUL_VALREF = ctx->valref;
    if (piaacmcparams.PIAACMC_fpmtype == 0)
    {
        piaacmcopticaldesign.fpmaskamptransm = dcimg[piaacmcopticaldesign.zoneaID].array.D[0];
    }

    dcimg[ctx->IDstatus].array.UI16[0] = 30;
    {
        char dirname[STRINGMAXLEN_DIRNAME];
        SNPRINTF_CHECK(dirname, STRINGMAXLEN_DIRNAME, "%s_linopt", piaacmcparams.piaacmcconfdir);
        PIAACMCsimul_savepiaacmcconf(dirname);
        EXECUTE_SYSTEM_COMMAND("rsync -au --progress %s/* ./%s/",
                               dirname,
                               piaacmcparams.piaacmcconfdir);
        if (piaacmcopticaldesign.PIAAmode == 1)
        {
            EXECUTE_SYSTEM_COMMAND("cp %s/piaa0Cmodes.fits %s/piaa0Cmodes.%04ld.fits",
                                   dirname,
                                   dirname,
                                   ctx->iter);
            EXECUTE_SYSTEM_COMMAND("cp %s/piaa0Fmodes.fits %s/piaa0Fmodes.%04ld.fits",
                                   dirname,
                                   dirname,
                                   ctx->iter);
            EXECUTE_SYSTEM_COMMAND("cp %s/piaa1Cmodes.fits %s/piaa1Cmodes.%04ld.fits",
                                   dirname,
                                   dirname,
                                   ctx->iter);
            EXECUTE_SYSTEM_COMMAND("cp %s/piaa1Fmodes.fits %s/piaa1Fmodes.%04ld.fits",
                                   dirname,
                                   dirname,
                                   ctx->iter);
        }
        EXECUTE_SYSTEM_COMMAND("cp %s/piaacmcparams.conf %s/piaacmcparams.%04ld.conf",
                               dirname,
                               dirname,
                               ctx->iter);
    }

    *contrastval_out = contrastval;
    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * PIAACMCsimul_run_linopt - Run regularized Jacobian + SVD line-search optimization loop.
 * @IDstatus: Shared-memory status telemetry image ID
 * @stopfile: Path to stop-signal sentinel file (`stopmode<mode>.txt`)
 *
 * Return: RETURN_SUCCESS on convergence or stop signal, or error code on failure.
 */
errno_t PIAACMCsimul_run_linopt(
    imageID     IDstatus,
    const char *stopfile)
{
    DEBUG_TRACE_FSTART();

    PIAACMC_linopt_ctx ctx = {
        .IDstatus   = IDstatus,
        .iter       = 0,
        .size1Dvec  = 0,
        .size1Dvec0 = 0,
        .valref     = 0.0,
        .bestval    = 0.0,
        .bestgain   = 0.0,
    };

    FUNC_CHECK_RETURN(linopt_init_reference_state(&ctx));

    int    iterOK      = 1;
    double oldval      = 1.0;
    dcimg[IDstatus].array.UI16[0] = 9;

    while (iterOK == 1)
    {
        imageID IDmodes;
        imageID IDmodes2D;
        imageID IDoptvec     = -1;
        long    NBlinoptgain = 0;
        double  contrastval  = 0.0;

        dcimg[IDstatus].array.UI16[0] = 10;
        printf("Iteration %ld/%ld\n", ctx.iter, piaacmcparams.linopt_NBiter);
        fflush(stdout);

        FUNC_CHECK_RETURN(create_3Dimage_ID("DHmodes",
                                            ctx.size1Dvec,
                                            1,
                                            piaacmcparams.linopt_number_param,
                                            &IDmodes));
        FUNC_CHECK_RETURN(create_2Dimage_ID("DHmodes2D",
                                            ctx.size1Dvec,
                                            piaacmcparams.linopt_number_param,
                                            &IDmodes2D));

        {
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(fname, "%s/linoptval.txt", piaacmcparams.piaacmcconfdir);
            FILE *fp = fopen(fname, "a");
            fprintf(fp,
                    "### PIAACMC_FPM_FASTDERIVATIVES = %d\n",
                    piaacmcparams.PIAACMC_FPM_FASTDERIVATIVES);
            fclose(fp);
        }
        dcimg[IDstatus].array.UI16[0] = 11;

        printf("Compute local derivatives of output vector against input "
               "focal plane mask zones (sags)\n");
        fflush(stdout);

        if (piaacmcparams.PIAACMC_FPM_FASTDERIVATIVES == 1)
        {
            FUNC_CHECK_RETURN(linopt_compute_fast_fpm_jacobian(&ctx, IDmodes, IDmodes2D));
        }
        else
        {
            FUNC_CHECK_RETURN(linopt_compute_numeric_jacobian(&ctx, IDmodes));
        }
        dcimg[IDstatus].array.UI16[0] = 16;

        FUNC_CHECK_RETURN(PIAACMCsimul_linopt_svd_linescan(&ctx, &IDoptvec, &NBlinoptgain));
        FUNC_CHECK_RETURN(linopt_commit_iteration(&ctx, IDoptvec, &contrastval));

        if (file_exists(stopfile) == 1 ||
            ctx.iter == piaacmcparams.linopt_NBiter ||
            (ctx.iter > 2 && contrastval > 0.98 * oldval) ||
            NBlinoptgain < 3)
        {
            iterOK = 0;
        }

        oldval = contrastval;
        ctx.iter++;

        printf("END OF LOOP ITERATION\n");
        fflush(stdout);
        dcimg[IDstatus].array.UI16[0] = 31;
    }

    printf(" ============ END OF OPTIMIZATION LOOP ======= \n");
    dcimg[IDstatus].array.UI16[0] = 32;

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
