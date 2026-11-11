/**
 * @file    optimizeLyotStop.c
 * @brief   Optimize multi-plane Lyot stop axial conjugations across Fresnel cube
 *
 * Can design both APLCMC and PIAACMC coronagraphs.
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "PIAACMCsimul/LyotStop/mkLyotMask.h"
#include "PIAACMCsimul/LyotStop/optimizeLyotStop.h"
#include "PIAACMCsimul/PIAACMCsimul.h"

#include "optimizeLyotStop_internal.h"

/**
 * free_lyot_opt_ctx - Release heap buffers held by a Lyot stop optimization context.
 * @ctx: Context whose member arrays should be freed
 */
static void free_lyot_opt_ctx(struct lyot_opt_ctx *ctx)
{
    free(ctx->totarray);
    free(ctx->tot2array);
    free(ctx->rinarray);
    free(ctx->routarray);
    free(ctx->zarray);
    free(ctx->rarray);
}

/**
 * init_lyot_opt_ctx - Allocate optimization arrays and compute annular zone radii and z planes.
 * @ctx:           Output context to initialize
 * @IDamp_name:    Input 2D/3D amplitude image name
 * @IDincohc_name: Input 3D incoherent off-axis intensity cube name
 * @zmin:          Minimum axial propagation distance [m]
 * @zmax:          Maximum axial propagation distance [m]
 */
static void init_lyot_opt_ctx(
    struct lyot_opt_ctx *ctx,
    const char *restrict IDamp_name,
    const char *restrict IDincohc_name,
    float                zmin,
    float                zmax)
{
    long NBz     = ctx->NBz;
    long NBmasks = ctx->NBmasks;

    ctx->zarray    = (float *) malloc(sizeof(float) * NBz);
    ctx->rinarray  = (float *) malloc(sizeof(float) * NBmasks);
    ctx->routarray = (float *) malloc(sizeof(float) * NBmasks);
    ctx->totarray  = (double *) malloc(sizeof(double) * NBmasks * NBz);
    ctx->tot2array = (double *) malloc(sizeof(double) * NBmasks * NBz);
    if ((ctx->zarray == NULL) || (ctx->rinarray == NULL) || (ctx->routarray == NULL) ||
        (ctx->totarray == NULL) || (ctx->tot2array == NULL))
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    ctx->routarray[0] = 1.0;
    ctx->rinarray[0]  = 1.0 - 1.0 / NBmasks;
    for (long m = 1; m < NBmasks; m++)
    {
        ctx->routarray[m] = ctx->rinarray[m - 1];
        ctx->rinarray[m] =
            ctx->routarray[m] - (1.0 - piaacmcopticaldesign.centObs1) / NBmasks;
    }
    ctx->rinarray[NBmasks - 1] = 0.0;

    for (long m = 0; m < NBmasks; m++)
    {
        DEBUG_TRACEPOINT("annulus %ld : %f - %f\n", m, ctx->routarray[m], ctx->rinarray[m]);
    }

    imageID IDa  = image_ID(IDamp_name, dcimg, dcnimg);
    ctx->xsize   = dcimg[IDa].md[0].size[0];
    ctx->ysize   = dcimg[IDa].md[0].size[1];
    ctx->xysize  = (uint64_t) ctx->xsize * ctx->ysize;
    ctx->nblambda = (dcimg[IDa].md[0].naxis == 3) ? dcimg[IDa].md[0].size[2] : 1;
    ctx->IDincohc = image_ID(IDincohc_name, dcimg, dcnimg);

    ctx->rarray = (float *) malloc(sizeof(float) * ctx->xysize);
    if (ctx->rarray == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    for (long l = 0; l < NBz; l++)
    {
        ctx->zarray[l] = zmin + (zmax - zmin) * l / (NBz - 1);
    }
}

/**
 * build_lyot_zone_map - Populate 2D radial coordinate buffer and annular LMzonemap image.
 * @ctx: Initialized Lyot optimization context
 *
 * Return: RETURN_SUCCESS or image/FITS error status.
 */
static errno_t build_lyot_zone_map(struct lyot_opt_ctx *ctx)
{
    FUNC_CHECK_RETURN(create_2Dimage_ID("LMzonemap", ctx->xsize, ctx->ysize, &ctx->IDzone));

    double rad_pix = piaacmcopticaldesign.beamrad / piaacmcopticaldesign.pixscale;
    for (uint32_t ii = 0; ii < ctx->xsize; ii++)
    {
        for (uint32_t jj = 0; jj < ctx->ysize; jj++)
        {
            uint64_t idx = (uint64_t) jj * ctx->xsize + ii;
            dcimg[ctx->IDzone].array.F[idx] = -2;

            double x         = (1.0 * ii - 0.5 * ctx->xsize) / rad_pix;
            double y         = (1.0 * jj - 0.5 * ctx->xsize) / rad_pix;
            double r         = sqrt(x * x + y * y);
            ctx->rarray[idx] = r;

            for (long m = 0; m < ctx->NBmasks; m++)
            {
                if ((r > ctx->rinarray[m] - 0.0001) && (r < ctx->routarray[m] + 0.0001))
                {
                    dcimg[ctx->IDzone].array.F[idx] = m;
                }
            }
        }
    }

    if (piaacmcparams.PIAACMC_save == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/LMzonemap.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("LMzonemap", fname));
    }

    return RETURN_SUCCESS;
}

/**
 * fill_single_annular_mask - Populate individual Lyot mask m from synthesized full LMask.
 * @ctx:        Lyot optimization context
 * @m:          Mask index
 * @IDm:        Output 2D mask image ID
 * @IDlyotmask: Full synthesized binary Lyot mask ID
 * @dr:         Radial overlap margin around zone boundaries
 */
static void fill_single_annular_mask(
    const struct lyot_opt_ctx *ctx,
    long                       m,
    imageID                    IDm,
    imageID                    IDlyotmask,
    float                      dr)
{
    double rad_pix = piaacmcopticaldesign.beamrad / piaacmcopticaldesign.pixscale;

    for (uint32_t ii = 0; ii < ctx->xsize; ii++)
    {
        for (uint32_t jj = 0; jj < ctx->ysize; jj++)
        {
            uint64_t idx = (uint64_t) jj * ctx->xsize + ii;
            double   x   = (1.0 * ii - 0.5 * ctx->xsize) / rad_pix;
            double   y   = (1.0 * jj - 0.5 * ctx->xsize) / rad_pix;
            double   r   = sqrt(x * x + y * y);

            if ((r > ctx->rinarray[m] - dr) && (r < ctx->routarray[m] + dr))
            {
                dcimg[IDm].array.F[idx] = dcimg[IDlyotmask].array.F[idx];
            }
            else
            {
                dcimg[IDm].array.F[idx] = 1.0;
            }

            if ((m == 0) && (r > 1.0))
            {
                dcimg[IDm].array.F[idx] = 0.0;
            }
        }
    }
}

/**
 * export_annular_lyot_masks - Partition LMask into cumulative per-plane optLM%02ld.fits files.
 * @ctx:        Lyot optimization context
 * @IDlyotmask: Synthesized binary Lyot mask image ID
 * @dr:         Radial overlap margin around zone boundaries
 *
 * Return: RETURN_SUCCESS or image/FITS error status.
 */
static errno_t export_annular_lyot_masks(
    const struct lyot_opt_ctx *ctx,
    imageID                    IDlyotmask,
    float                      dr)
{
    imageID IDlscumul;
    FUNC_CHECK_RETURN(create_2Dimage_ID("LMcumul", ctx->xsize, ctx->ysize, &IDlscumul));

    for (uint64_t ii = 0; ii < ctx->xysize; ii++)
    {
        dcimg[IDlscumul].array.F[ii] = 1.0;
    }

    for (long m = 0; m < ctx->NBmasks; m++)
    {
        char name[STRINGMAXLEN_IMGNAME];
        WRITE_IMAGENAME(name, "optLM%02ld", m);

        imageID IDm;
        FUNC_CHECK_RETURN(create_2Dimage_ID(name, ctx->xsize, ctx->ysize, &IDm));

        fill_single_annular_mask(ctx, m, IDm, IDlyotmask, dr);

        for (uint64_t ii = 0; ii < ctx->xysize; ii++)
        {
            dcimg[IDm].array.F[ii] *= dcimg[IDlscumul].array.F[ii];
            dcimg[IDlscumul].array.F[ii] = dcimg[IDm].array.F[ii];
        }

        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/optLM%02ld.fits", piaacmcparams.piaacmcconfdir, m);
        FUNC_CHECK_RETURN(save_fits(name, fname));
    }

    FUNC_CHECK_RETURN(delete_image_ID("LMcumul", DELETE_IMAGE_ERRMODE_WARNING));
    return RETURN_SUCCESS;
}

/**
 * optimizeLyotStop - Optimize Lyot stop axial conjugations and masks from zmin to zmax.
 * @IDamp_name:    Input 2D/3D amplitude image name
 * @IDpha_name:    Input 2D/3D phase image name
 * @IDincohc_name: Input 3D incoherent off-axis intensity cube name
 * @zmin:          Minimum propagation distance [m]
 * @zmax:          Maximum propagation distance [m]
 * @throughput:    Target geometric throughput of Lyot stop(s)
 * @NBz:           Number of discrete propagation planes between zmin and zmax
 * @NBmasks:       Number of Lyot stop planes
 * @outratioval:   Optional output metric placeholder (may be NULL)
 *
 * Return: RETURN_SUCCESS on completion, or milk error code on failure.
 */
errno_t optimizeLyotStop(
    const char *restrict IDamp_name,
    const char *restrict IDpha_name,
    const char *restrict IDincohc_name,
    float                zmin,
    float                zmax,
    double               throughput,
    long                 NBz,
    long                 NBmasks,
    double              *outratioval)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %s %s %s %f %f %lf %ld %ld",
                     IDamp_name,
                     IDpha_name,
                     IDincohc_name,
                     zmin,
                     zmax,
                     throughput,
                     NBz,
                     NBmasks);

    double ratio = 1.0;
    float  dr    = 0.02;
    double alpha = 1.01;

    struct lyot_opt_ctx ctx = {0};
    ctx.NBz                 = NBz;
    ctx.NBmasks             = NBmasks;

    init_lyot_opt_ctx(&ctx, IDamp_name, IDincohc_name, zmin, zmax);

    errno_t status = build_lyot_zone_map(&ctx);
    if (status == RETURN_SUCCESS)
    {
        status = lyot_opt_propagate_onaxis_cube(&ctx, IDamp_name, IDpha_name);
    }
    if (status == RETURN_SUCCESS)
    {
        status = lyot_opt_select_conjugations(&ctx, alpha);
    }

    imageID IDlyotmask = -1;
    if (status == RETURN_SUCCESS)
    {
        status = mkLyotMask("LcombOA", "Lcomb", "LMzonemap", throughput, "LMask", &IDlyotmask);
    }
    if ((status == RETURN_SUCCESS) && (piaacmcparams.PIAACMC_save == 1))
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/LMask.fits", piaacmcparams.piaacmcconfdir);
        status = save_fits("LMask", fname);
    }
    if (status == RETURN_SUCCESS)
    {
        status = delete_image_ID("Lcomb", DELETE_IMAGE_ERRMODE_WARNING);
    }
    if (status == RETURN_SUCCESS)
    {
        status = export_annular_lyot_masks(&ctx, IDlyotmask, dr);
    }

    free_lyot_opt_ctx(&ctx);
    if (status != RETURN_SUCCESS)
    {
        return status;
    }

    FUNC_CHECK_RETURN(delete_image_ID("LMzonemap", DELETE_IMAGE_ERRMODE_WARNING));

    if (outratioval != NULL)
    {
        *outratioval = ratio;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
