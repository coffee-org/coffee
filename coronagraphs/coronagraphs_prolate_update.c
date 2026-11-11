/**
 * @file    coronagraphs_prolate_update.c
 * @brief   Iterative zoomed-DFT sweep and refinement of precomputed 2D prolate apodizers
 */

#include <stdio.h>

#include "coffee_compat.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "info/info.h"

#include "coronagraphs_internal.h"

typedef struct prolate_upd_ctx
{
    int    co_flag;
    double beamradpix;
    double centralObs;
    long   size;
    long   NBiter1;
    long   ID0;
    long   ID1;
    long   ID2;
} prolate_upd_ctx_t;

/**
 * prolate_format_paths - Format APLCapo reference FITS, info, and profile paths
 * @co_flag:    1 if central obstruction is active, 0 otherwise
 * @maskld:     Focal plane mask radius [lambda/D]
 * @centralObs: Central obstruction ratio
 * @size:       Array size [pix]
 * @fname_raw:  Output buffer for raw .ref.gz path
 * @fname_info: Output buffer for .ref.info path (may be NULL)
 * @fname_prof: Output buffer for .ref.prof path (may be NULL)
 */
static void prolate_format_paths(
    int     co_flag,
    double  maskld,
    double  centralObs,
    long    size,
    char   *fname_raw,
    char   *fname_info,
    char   *fname_prof)
{
    if (co_flag == 0)
    {
        sprintf(
            fname_raw,
            "%s/APLCapo/raw/APLCapo_%.3f.%ld.ref.gz",
            CORONAGRAPHSDATALOCAL,
            maskld,
            size);
        if (fname_info != NULL)
        {
            sprintf(
                fname_info,
                "%s/APLCapo/APLCapo_%.3f.%ld.ref.info",
                CORONAGRAPHSDATALOCAL,
                maskld,
                size);
        }
        if (fname_prof != NULL)
        {
            sprintf(
                fname_prof,
                "%s/APLCapo/raw/APLCapo_%.3f.%ld.ref.prof",
                CORONAGRAPHSDATALOCAL,
                maskld,
                size);
        }
    }
    else
    {
        sprintf(
            fname_raw,
            "%s/APLCapo/raw/APLCapo_%.3f.%.3f.%ld.ref.gz",
            CORONAGRAPHSDATALOCAL,
            maskld,
            centralObs,
            size);
        if (fname_info != NULL)
        {
            sprintf(
                fname_info,
                "%s/APLCapo/APLCapo_%.3f.%.3f.%ld.ref.info",
                CORONAGRAPHSDATALOCAL,
                maskld,
                centralObs,
                size);
        }
        if (fname_prof != NULL)
        {
            sprintf(
                fname_prof,
                "%s/APLCapo/raw/APLCapo_%.3f.%.3f.%ld.ref.prof",
                CORONAGRAPHSDATALOCAL,
                maskld,
                centralObs,
                size);
        }
    }
}

/**
 * prolate_update_inner_step - Run one DFT optimization pass for the active MASKSIZELD
 * @ctx:   Update sweep context
 * @iter:  Sub-iteration index (0 or 1)
 * @iter1: Outer iteration counter
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE on error.
 */
static errno_t prolate_update_inner_step(
    const prolate_upd_ctx_t *ctx,
    long                     iter,
    long                     iter1)
{
    char   fname[500];
    char   fnameinfo[500];
    char   fnameprof[500];
    double transm;
    double err = 0.0;
    long   ID;

    printf(" ---- %f lambda/d  [%ld / %ld] ----\n", MASKSIZELD, iter1, ctx->NBiter1);
    prolate_format_paths(
        ctx->co_flag,
        MASKSIZELD,
        ctx->centralObs,
        ctx->size,
        fname,
        fnameinfo,
        fnameprof);

    load_fits(fname, "apostart", 1, NULL);
    printf(
        ".... STARTING OPTIMIZATION [%s] (%f %f %f)... \n",
        fname,
        MASKSIZELD,
        CORONAGRAPHS_PIXSCALE,
        ctx->centralObs);

    if (coronagraph_make_2Dprolate_DFT(
            MASKSIZELD / CORONAGRAPHS_PIXSCALE,
            ctx->beamradpix,
            ctx->centralObs,
            "out",
            ctx->size,
            "NULLim",
            &transm) != RETURN_SUCCESS)
    {
        FUNC_RETURN_FAILURE("Call to coronagraph_make_2Dprolate_DFT failed");
    }

    save_fl_fits("out", fname);
    ID = image_ID("out", dcimg, dcnimg);
    for (long ii = 0; ii < ctx->size * ctx->size; ii++)
    {
        dcimg[ctx->ID1].array.F[ii] = dcimg[ctx->ID2].array.F[ii];
        dcimg[ctx->ID2].array.F[ii] = dcimg[ID].array.F[ii];
        double tmp1                 = dcimg[ctx->ID1].array.F[ii] - dcimg[ctx->ID2].array.F[ii];
        err += tmp1 * tmp1;
    }
    if (iter1 == 0)
    {
        err = 0.0;
    }

    profile("out", "out.prof", ctx->size / 2, ctx->size / 2, 1.0, ctx->size / 2);
    EXECUTE_SYSTEM_COMMAND("mv out.prof %s", fnameprof);
    EXECUTE_SYSTEM_COMMAND("cp APLCapo.%.3f.%.3f.info %s", MASKSIZELD, ctx->centralObs, fnameinfo);

    delete_image_ID("apostart", DELETE_IMAGE_ERRMODE_WARNING);

    if (iter == 1)
    {
        FILE *fp = fopen("prol_conv.log", "a");
        fprintf(fp, "%ld %.18f %.18g\n", iter1, transm, err);
        fclose(fp);
    }

    return RETURN_SUCCESS;
}

/**
 * prolate_update_masksize - Update prolate files for the current MASKSIZELD step
 * @ctx:          Update sweep context
 * @masksizestep: Mask size increment [lambda/D]
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE on error.
 */
static errno_t prolate_update_masksize(
    prolate_upd_ctx_t *ctx,
    float              masksizestep)
{
    char   fname[500];
    char   fnameprev[500];
    char   fnameinfo[500];
    char   fnameinfoprev[500];
    char   fnamenew[500];
    double FFcoeff = 0.0;

    printf("--------- MASK SIZE = %f --------\n", MASKSIZELD);
    prolate_format_paths(
        ctx->co_flag, MASKSIZELD, ctx->centralObs, ctx->size, fname, fnameinfo, NULL);
    prolate_format_paths(
        ctx->co_flag,
        MASKSIZELD - masksizestep,
        ctx->centralObs,
        ctx->size,
        fnameprev,
        fnameinfoprev,
        NULL);
    sprintf(
        fnamenew,
        "%s/APLCapo/APLCapo_%.3f.%.3f.%ld.ref.new1",
        CORONAGRAPHSDATALOCAL,
        MASKSIZELD,
        ctx->centralObs,
        ctx->size);

    if (file_exists(fnamenew) != 0)
    {
        return RETURN_SUCCESS;
    }

    printf("Mask size = %f l/D\n", MASKSIZELD);
    if ((file_exists(fnameprev) == 1) && (file_exists(fnameinfoprev) == 1) &&
        ((file_exists(fname) == 0) || (file_exists(fnameinfo) == 0)))
    {
        printf("COPY FROM PREVIOUS MASK SIZE\n");
        EXECUTE_SYSTEM_COMMAND("cp %s %s", fnameprev, fname);
        EXECUTE_SYSTEM_COMMAND("cp %s %s", fnameprev, fname);
        EXECUTE_SYSTEM_COMMAND("cp %s %s", fnameinfoprev, fnameinfo);
    }

    create_2Dimage_ID("apo", ctx->size, ctx->size, &ctx->ID0);
    create_2Dimage_ID("apo1", ctx->size, ctx->size, &ctx->ID1);
    create_2Dimage_ID("apo2", ctx->size, ctx->size, &ctx->ID2);
    for (long ii = 0; ii < ctx->size * ctx->size; ii++)
    {
        dcimg[ctx->ID2].array.F[ii] = 0.0;
    }

    {
        FILE *fp = fopen("prol_conv.log", "w");
        fclose(fp);
    }

    long iter1 = 0;
    while (iter1 < ctx->NBiter1)
    {
        iter1++;
        long iter = 0;
        while ((iter != 2) && (iter1 < ctx->NBiter1 + 1))
        {
            if (prolate_update_inner_step(ctx, iter, iter1) != RETURN_SUCCESS)
            {
                return RETURN_FAILURE;
            }
            iter++;
            iter1++;
        }

        for (long ii = 0; ii < ctx->size * ctx->size; ii++)
        {
            dcimg[ctx->ID0].array.F[ii] =
                dcimg[ctx->ID2].array.F[ii] +
                FFcoeff * (dcimg[ctx->ID2].array.F[ii] - dcimg[ctx->ID1].array.F[ii]);
        }
        save_fl_fits("apo", fname);
    }

    delete_image_ID("apo", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("apo1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("apo2", DELETE_IMAGE_ERRMODE_WARNING);
    EXECUTE_SYSTEM_COMMAND("touch %s", fnamenew);

    return RETURN_SUCCESS;
}

/**
 * coronagraph_update_2Dprolate - Iteratively refine 2D prolate apodizers across mask sizes
 * @masksizeld: Initial mask radius [lambda/D]
 * @beamradpix: Entrance pupil beam radius [pix]
 * @centralObs: Central obstruction ratio
 * @zfactor:    Zoomed DFT oversampling factor
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE on error.
 */
errno_t coronagraph_update_2Dprolate(
    double masksizeld,
    double beamradpix,
    double centralObs,
    double zfactor)
{
    DEBUG_TRACE_FSTART();

    prolate_upd_ctx_t ctx = {
        .co_flag    = 0,
        .beamradpix = beamradpix,
        .centralObs = centralObs,
        .size       = CORONAGRAPHS_ARRAYSIZE,
        .NBiter1    = 1,
    };
    float masksizestep = 0.01f;
    long  ID;

    g_coro.DFTZFACTOR = zfactor;
    if ((ID = variable_ID("DFTZFACTOR")) != -1)
    {
        g_coro.DFTZFACTOR = coffee_dcvar_f(ID);
    }

    if (centralObs > 0.001)
    {
        ctx.co_flag = 1;
        printf("Central obstruction = %f\n", centralObs);
    }

    MASKSIZELD = masksizeld;
    for (MASKSIZELD = 0.8; MASKSIZELD < 5.01; MASKSIZELD += masksizestep)
    {
        if (prolate_update_masksize(&ctx, masksizestep) != RETURN_SUCCESS)
        {
            return RETURN_FAILURE;
        }
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
