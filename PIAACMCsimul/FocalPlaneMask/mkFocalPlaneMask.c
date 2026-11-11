/**
 * @file    mkFocalPlaneMask.c
 * @brief   PIAA-type coronagraph design, make focal plane mask
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
#include "OpticsMaterials/OpticsMaterials.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"

#include "FocalPlaneMask/mkFocalPlaneMask.h"
#include "FocalPlaneMask/mkFocalPlaneMask_internal.h"

/**
 * mkFocalPlaneMask_alloc_zones - Allocate and precompute polychromatic zone phasor tables
 * @ctx: FPM synthesis context receiving allocated thickness, amplitude, and phase arrays.
 */
static void mkFocalPlaneMask_alloc_zones(struct fpm_mask_ctx *ctx)
{
    printf("Allocate memory\n");
    fflush(stdout);

    size_t nelem     = (size_t) piaacmcopticaldesign.focmNBzone * (size_t) ctx->nblambda;
    ctx->tarray      = (double *) malloc(sizeof(double) * nelem);
    ctx->aarray      = (double *) malloc(sizeof(double) * nelem);
    ctx->phaarray    = (double *) malloc(sizeof(double) * nelem);
    ctx->cosphaarray = (double *) malloc(sizeof(double) * nelem);
    ctx->sinphaarray = (double *) malloc(sizeof(double) * nelem);

    if (ctx->tarray == NULL || ctx->aarray == NULL || ctx->phaarray == NULL
        || ctx->cosphaarray == NULL || ctx->sinphaarray == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    printf(
        "Precompute zones phase shifts  %ld zones, %d wavelengths\n",
        piaacmcopticaldesign.focmNBzone,
        ctx->nblambda);
    fflush(stdout);

    for (int k = 0; k < ctx->nblambda; k++)
    {
        for (long zi = 0; zi < piaacmcopticaldesign.focmNBzone; zi++)
        {
            long idx = piaacmcopticaldesign.focmNBzone * k + zi;

            ctx->tarray[idx]   = dcimg[piaacmcopticaldesign.zonezID].array.D[zi];
            ctx->aarray[idx]   = dcimg[piaacmcopticaldesign.zoneaID].array.D[zi];
            ctx->phaarray[idx] = OpticsMaterials_pha_lambda(
                piaacmcopticaldesign.fpmmaterial_code,
                ctx->tarray[idx],
                piaacmcopticalsystem.lambdaarray[k]);
            ctx->cosphaarray[idx] = (double) cosf((float) ctx->phaarray[idx]);
            ctx->sinphaarray[idx] = (double) sinf((float) ctx->phaarray[idx]);
        }
    }
}

/**
 * mkFocalPlaneMask_free_zones - Release precomputed zone phasor working buffers
 * @ctx: FPM synthesis context whose arrays are freed.
 */
static void mkFocalPlaneMask_free_zones(struct fpm_mask_ctx *ctx)
{
    free(ctx->tarray);
    free(ctx->aarray);
    free(ctx->phaarray);
    free(ctx->cosphaarray);
    free(ctx->sinphaarray);
}

/**
 * mkFocalPlaneMask_save_diagnostics - Export 2D/3D FPM sag, zone, amplitude, and phase FITS maps
 * @ctx: Completed FPM synthesis context.
 *
 * Return: RETURN_SUCCESS on all FITS writes and cleanup, error code otherwise.
 */
static errno_t mkFocalPlaneMask_save_diagnostics(const struct fpm_mask_ctx *ctx)
{
    char fname[STRINGMAXLEN_FULLFILENAME];

    PIAACMCsimul_update_fnamedescr();

    SNPRINTF_CHECK(
        fname,
        STRINGMAXLEN_FULLFILENAME,
        "%s/fpm_sagmap2D.%s.fits.gz",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.fnamedescr);
    FUNC_CHECK_RETURN(save_fits("fpmsag", fname));

    SNPRINTF_CHECK(
        fname,
        STRINGMAXLEN_FULLFILENAME,
        "%s/fpm_zonemap2D.%s.fits.gz",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.fnamedescr);
    FUNC_CHECK_RETURN(save_fits("fpmzone", fname));

    imageID IDm;
    FUNC_CHECK_RETURN(create_3DCimage_ID("fpmCA", ctx->size, ctx->size, ctx->nblambda, &IDm));

    for (int k = 0; k < ctx->nblambda; k++)
    {
        for (uint32_t ii = 0; ii < ctx->size; ii++)
        {
            for (uint32_t jj = 0; jj < ctx->size; jj++)
            {
                uint64_t pixidx = (uint64_t) k * ctx->size2 + (uint64_t) jj * ctx->size + ii;

                dcimg[IDm].array.CF[pixidx].re = 1.0f - dcimg[ctx->ID].array.CF[pixidx].re;
                dcimg[IDm].array.CF[pixidx].im = -dcimg[ctx->ID].array.CF[pixidx].im;
            }
        }
    }

    mk_amph_from_complex("fpmCA", "tfpma", "tfpmp", 0);
    FUNC_CHECK_RETURN(delete_image_ID("fpmCA", DELETE_IMAGE_ERRMODE_WARNING));

    SNPRINTF_CHECK(
        fname,
        STRINGMAXLEN_FULLFILENAME,
        "%s/fpm_CAampmap2D.%s.fits.gz",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.fnamedescr);
    FUNC_CHECK_RETURN(save_fits("tfpma", fname));

    SNPRINTF_CHECK(
        fname,
        STRINGMAXLEN_FULLFILENAME,
        "%s/fpm_CAphamap2D.%s.fits.gz",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.fnamedescr);
    FUNC_CHECK_RETURN(save_fits("tfpmp", fname));

    FUNC_CHECK_RETURN(delete_image_ID("tfpma", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("tfpmp", DELETE_IMAGE_ERRMODE_WARNING));

    return RETURN_SUCCESS;
}

/**
 * mkFocalPlaneMask - Synthesize 3D polychromatic complex-amplitude focal plane mask cube
 * @IDzonemap_name: Input 2D UINT16 high-resolution zone index map image name.
 * @ID_name:        Output 3D complex-float focal plane mask cube image name.
 * @mode:           Synthesis mode (-1 or > focmNBzone for 1-FPM; 1..focmNBzone for zone impulse).
 * @saveMask:       Non-zero to write diagnostic sag, zone, amplitude, and phase FITS files.
 * @outID:          Optional output pointer receiving created cube image identifier.
 *
 * Return: RETURN_SUCCESS on completion, error code otherwise.
 */
errno_t mkFocalPlaneMask(
    const char *IDzonemap_name,
    const char *ID_name,
    int         mode,
    int         saveMask,
    imageID    *outID)
{
    DEBUG_TRACE_FSTART();

    struct fpm_mask_ctx ctx;
    memset(&ctx, 0, sizeof(ctx));

    ctx.size     = piaacmcopticalsystem.size;
    ctx.size2    = (uint64_t) ctx.size * ctx.size;
    ctx.nblambda = piaacmcopticalsystem.nblambda;
    ctx.mode     = mode;
    ctx.NBsubPix = 64;
    ctx.IDz      = image_ID(IDzonemap_name, dcimg, dcnimg);

    FUNC_CHECK_RETURN(create_3DCimage_ID(ID_name, ctx.size, ctx.size, ctx.nblambda, &ctx.ID));
    FUNC_CHECK_RETURN(create_3Dimage_ID("fpmsag", ctx.size, ctx.size, ctx.nblambda, &ctx.IDsag));
    FUNC_CHECK_RETURN(create_3Dimage_ID("fpmzone", ctx.size, ctx.size, ctx.nblambda, &ctx.IDzone));

    if (piaacmcopticaldesign.NBrings > 2)
    {
        ctx.CentCone  = 1;
        ctx.OuterCone = 1;
    }
    if (fabs(piaacmcopticaldesign.fpmOuterConeZ) < 1.0e-12)
    {
        ctx.OuterCone = 0;
    }
    if (mode == -1 || mode > piaacmcopticaldesign.focmNBzone)
    {
        ctx.FPMmode = 1;
    }

    printf(
        "===================== Make focal plane mask  %s %s   mode=%d/%ld   [%d %d]\n",
        IDzonemap_name,
        ID_name,
        mode,
        piaacmcopticaldesign.focmNBzone,
        (int) ctx.CentCone,
        (int) ctx.OuterCone);

    double fpscale0 = (2.0 * piaacmcopticaldesign.beamrad / piaacmcopticaldesign.pixscale)
                      / piaacmcopticaldesign.size / piaacmcopticaldesign.fpzfactor
                      * piaacmcopticalsystem.lambdaarray[0] * piaacmcopticaldesign.Fratio;
    printf(
        "piaacmcopticaldesign.fpmRad = %g m    fpscale[0] = %g m/pix   mode = %d\n",
        piaacmcopticaldesign.fpmRad,
        fpscale0,
        mode);

    mkFocalPlaneMask_alloc_zones(&ctx);
    mkFocalPlaneMask_compute_cube(&ctx);

    printf("===================== focal plane mask  : Done\n");
    fflush(stdout);

    if (saveMask == 1)
    {
        errno_t status = mkFocalPlaneMask_save_diagnostics(&ctx);
        if (status != RETURN_SUCCESS)
        {
            mkFocalPlaneMask_free_zones(&ctx);
            return status;
        }
    }

    errno_t del_status = delete_image_ID("fpmsag", DELETE_IMAGE_ERRMODE_WARNING);
    mkFocalPlaneMask_free_zones(&ctx);
    FUNC_CHECK_RETURN(del_status);

    if (outID != NULL)
    {
        *outID = ctx.ID;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
