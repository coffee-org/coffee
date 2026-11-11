/**
 * @file    coronagraphs_prolate_dft.c
 * @brief   High-definition zoomed-DFT 2D prolate spheroidal apodizer computation
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "coffee_compat.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "fft/fft.h"
#include "image_gen/image_gen.h"

#include "coronagraphs_internal.h"

typedef struct prolate_dft_ctx
{
    long        size;
    uint64_t    size2;
    double      fpmradpix;
    double      beamradpix;
    double      centralObs;
    imageID     IDfpmask;
    imageID     IDfpmaskc;
    imageID     IDfpmaskz;
    imageID     IDpupa0;
    imageID     IDpupp0;
    imageID     IDpupa0m;
    imageID     IDpupp0m;
    imageID     IDprol;
    imageID     IDprolr;
    imageID     IDproli;
    imageID     IDprolp;
    long double total0;
    double      transm;
    double      peak;
} prolate_dft_ctx_t;

/**
 * prolate_dft_create_fpm - Build zoomed perturbed focal plane mask images
 * @ctx: DFT prolate solver context
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE on allocation/image failure.
 */
static errno_t prolate_dft_create_fpm(prolate_dft_ctx_t *ctx)
{
    long    fpmshape_n = 20;
    double *fpmshape_ra = (double *) malloc(sizeof(double) * fpmshape_n);
    double *fpmshape_ka = (double *) malloc(sizeof(double) * fpmshape_n);
    double *fpmshape_pa = (double *) malloc(sizeof(double) * fpmshape_n);

    if ((fpmshape_ra == NULL) || (fpmshape_ka == NULL) || (fpmshape_pa == NULL))
    {
        free(fpmshape_ra);
        free(fpmshape_ka);
        free(fpmshape_pa);
        PRINT_ERROR("malloc returns NULL pointer");
        return RETURN_FAILURE;
    }

    for (long k = 0; k < fpmshape_n; k++)
    {
        char    vname[STRINGMAXLEN_IMGNAME];
        imageID ID;

        fpmshape_ra[k] = 0.0;
        fpmshape_ka[k] = 0.0;
        fpmshape_pa[k] = 0.0;

        WRITE_IMAGENAME(vname, "FPMSHAPE_%ld_r", k);
        if ((ID = variable_ID(vname)) != -1)
        {
            fpmshape_ra[k] = coffee_dcvar_f(ID);
        }

        WRITE_IMAGENAME(vname, "FPMSHAPE_%ld_k", k);
        if ((ID = variable_ID(vname)) != -1)
        {
            fpmshape_ka[k] = coffee_dcvar_f(ID);
        }

        WRITE_IMAGENAME(vname, "FPMSHAPE_%ld_p", k);
        if ((ID = variable_ID(vname)) != -1)
        {
            fpmshape_pa[k] = coffee_dcvar_f(ID);
        }

        printf(
            "FPM SHAPE TERM %02ld : %f %f %f\n",
            k,
            fpmshape_ra[k],
            fpmshape_ka[k],
            fpmshape_pa[k]);
    }

    ctx->IDfpmask  = make_subpixdisk(
        "FPmask",
        ctx->size,
        ctx->size,
        ctx->size / 2,
        ctx->size / 2,
        ctx->fpmradpix);
    ctx->IDfpmaskz = make_subpixdisk_perturb(
        "FPmaskz",
        ctx->size,
        ctx->size,
        ctx->size / 2,
        ctx->size / 2,
        ctx->fpmradpix * g_coro.DFTZFACTOR,
        fpmshape_n,
        fpmshape_ra,
        fpmshape_ka,
        fpmshape_pa);

    free(fpmshape_ra);
    free(fpmshape_ka);
    free(fpmshape_pa);

    FUNC_CHECK_RETURN(create_2DCimage_ID("_fpmz", ctx->size, ctx->size, &ctx->IDfpmaskc));

    for (uint64_t ii = 0; ii < ctx->size2; ii++)
    {
        dcimg[ctx->IDfpmask].array.F[ii]     = dcimg[ctx->IDfpmaskz].array.F[ii];
        dcimg[ctx->IDfpmaskc].array.CF[ii].re = dcimg[ctx->IDfpmaskz].array.F[ii];
        dcimg[ctx->IDfpmaskc].array.CF[ii].im = 0.0;
    }

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/FPmask.tmp.fits", dcsavedir);
        save_fl_fits("FPmask", fname);
    }

    return RETURN_SUCCESS;
}

/**
 * prolate_dft_init_pupil_and_apo - Prepare entrance pupil and seed prolate profile
 * @ctx:          DFT prolate solver context
 * @outname:      Output prolate image name
 * @pupmask_name: Optional extra pupil mask image name
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE on error.
 */
static errno_t prolate_dft_init_pupil_and_apo(
    prolate_dft_ctx_t *ctx,
    const char        *outname,
    const char        *pupmask_name)
{
    long size    = ctx->size;
    int  co_flag = (ctx->centralObs > 0.001) ? 1 : 0;

    FUNC_CHECK_RETURN(create_2Dimage_ID("pupp0", size, size, &ctx->IDpupp0));

    ctx->IDpupa0 = image_ID("pupa0", dcimg, dcnimg);
    if (ctx->IDpupa0 != -1)
    {
        if (dcimg[ctx->IDpupa0].md[0].size[0] != (uint32_t) size)
        {
            printf("ERROR: pupa0 should be %ld x %ld\n", size, size);
            return RETURN_FAILURE;
        }
    }
    else
    {
        ctx->IDpupa0 = make_disk("pupa0", size, size, size / 2, size / 2, ctx->beamradpix);
        if (co_flag == 1)
        {
            imageID IDco = make_disk(
                "pupa0co",
                size,
                size,
                size / 2,
                size / 2,
                ctx->beamradpix * ctx->centralObs);
            for (uint64_t ii = 0; ii < ctx->size2; ii++)
            {
                dcimg[ctx->IDpupa0].array.F[ii] -= dcimg[IDco].array.F[ii];
            }
            FUNC_CHECK_RETURN(delete_image_ID("pupa0co", DELETE_IMAGE_ERRMODE_WARNING));
        }

        imageID IDpupmask = image_ID(pupmask_name, dcimg, dcnimg);
        if (IDpupmask != -1)
        {
            for (uint64_t ii = 0; ii < ctx->size2; ii++)
            {
                dcimg[ctx->IDpupa0].array.F[ii] *= dcimg[IDpupmask].array.F[ii];
            }
        }
    }

    ctx->total0 = 0.0;
    for (uint64_t ii = 0; ii < ctx->size2; ii++)
    {
        long double v1 = dcimg[ctx->IDpupa0].array.F[ii];
        ctx->total0 += v1 * v1;
    }

    FUNC_CHECK_RETURN(create_2Dimage_ID("pupa0m", size, size, &ctx->IDpupa0m));
    FUNC_CHECK_RETURN(create_2Dimage_ID("pupp0m", size, size, &ctx->IDpupp0m));
    FUNC_CHECK_RETURN(create_2Dimage_ID(outname, size, size, &ctx->IDprol));

    for (uint32_t ii = 0; ii < (uint32_t) size; ii++)
    {
        for (uint32_t jj = 0; jj < (uint32_t) size; jj++)
        {
            uint64_t idx = (uint64_t) jj * size + ii;
            if (dcimg[ctx->IDpupa0].array.F[idx] > 0.0001)
            {
                double x                  = 1.0 * ii - size / 2;
                double y                  = 1.0 * jj - size / 2;
                double r2                 = (x * x + y * y) / ctx->beamradpix / ctx->beamradpix;
                double val                = exp(-r2 * 0.73 / 5.0 * ctx->fpmradpix);
                dcimg[ctx->IDprolr].array.F[idx] = val;
                dcimg[ctx->IDproli].array.F[idx] = 0.0;
                dcimg[ctx->IDprol].array.F[idx]  = val;
            }
            else
            {
                dcimg[ctx->IDprolr].array.F[idx] = 0.0;
                dcimg[ctx->IDproli].array.F[idx] = 0.0;
                dcimg[ctx->IDprol].array.F[idx]  = 0.0;
            }
        }
    }

    imageID ID = image_ID("apostart", dcimg, dcnimg);
    if (ID != -1)
    {
        for (uint64_t ii = 0; ii < ctx->size2; ii++)
        {
            dcimg[ctx->IDprolr].array.F[ii] = dcimg[ID].array.F[ii];
            dcimg[ctx->IDprol].array.F[ii]  = dcimg[ID].array.F[ii];
        }
    }
    FUNC_CHECK_RETURN(delete_image_ID("apostart", DELETE_IMAGE_ERRMODE_WARNING));

    return RETURN_SUCCESS;
}

/**
 * prolate_dft_write_and_normalize - Save APLCapo info file and normalize prolate iterate
 * @ctx:      DFT prolate solver context
 * @fpm_thru: Fraction of focal energy inside focal plane mask
 */
static void prolate_dft_write_and_normalize(
    prolate_dft_ctx_t *ctx,
    double             fpm_thru)
{
    uint64_t size2 = ctx->size2;
    imageID  IDr   = image_ID("pr2", dcimg, dcnimg);
    imageID  IDi   = image_ID("pi2", dcimg, dcnimg);
    char     fname[STRINGMAXLEN_FULLFILENAME];
    FILE    *fp;

    ctx->peak = 0.0;
    for (uint64_t ii = 0; ii < size2; ii++)
    {
        if ((dcimg[ctx->IDpupa0].array.F[ii] > 0.0001) &&
            (fabs(dcimg[IDr].array.F[ii]) > fabs(ctx->peak)))
        {
            ctx->peak = dcimg[IDr].array.F[ii];
        }
    }
    printf("%f\n", ctx->peak);
    printf("------- MASK SIZE = %f l/D ------\n", MASKSIZELD);

    WRITE_FULLFILENAME(fname, "%s/APLCapo.%.3f.%.3f.info", dcsavedir, MASKSIZELD, ctx->centralObs);
    fp = fopen(fname, "w");
    fprintf(fp, "%f %.18f %.18f %.18f\n", MASKSIZELD, fpm_thru, ctx->transm, ctx->peak);
    fprintf(fp, "# Created by %s line %d\n", __FILE__, __LINE__);
    fprintf(fp, "# col 1: Mask size in SYSTEM l/D\n");
    fprintf(fp, "# col 2: Fraction of light going through FPmask\n");
    fprintf(fp, "# col 3: Prolate throughput if not done by PIAA\n");
    fprintf(fp, "# col 4: Peak\n");
    fprintf(
        fp,
        "# ideal complex transmission for FPM is MASKCAMULT = %.18f\n",
        -(1.0 - ctx->peak) / ctx->peak);
    fclose(fp);
    printf("IDEAL COMPLEX TRANSMISSION = %g\n", -(1.0 - ctx->peak) / ctx->peak);

    create_variable_ID("APLCmaskCtransm", -(1.0 - ctx->peak) / ctx->peak);

    for (uint64_t ii = 0; ii < size2; ii++)
    {
        if (dcimg[ctx->IDpupa0].array.F[ii] > 0.0001)
        {
            dcimg[ctx->IDprolr].array.F[ii] = dcimg[IDr].array.F[ii];
            dcimg[ctx->IDproli].array.F[ii] = dcimg[IDi].array.F[ii];
        }
    }

    for (uint64_t ii = 0; ii < size2; ii++)
    {
        double re                       = dcimg[ctx->IDprolr].array.F[ii];
        double im                       = dcimg[ctx->IDproli].array.F[ii];
        dcimg[ctx->IDprol].array.F[ii] = sqrt(re * re + im * im);
    }

    ctx->peak = 0.0;
    for (uint64_t ii = 0; ii < size2; ii++)
    {
        if ((dcimg[ctx->IDpupa0].array.F[ii] > 0.0001) &&
            (fabs(dcimg[ctx->IDprol].array.F[ii]) > fabs(ctx->peak)))
        {
            ctx->peak = dcimg[ctx->IDprol].array.F[ii];
        }
    }

    for (uint64_t ii = 0; ii < size2; ii++)
    {
        if (dcimg[ctx->IDpupa0].array.F[ii] > 0.001)
        {
            dcimg[ctx->IDprol].array.F[ii] /= ctx->peak;
        }
        else
        {
            dcimg[ctx->IDprol].array.F[ii] = 0.0;
        }
    }
}

/**
 * prolate_dft_step - Execute one zoomed-DFT prolate power-method iteration
 * @ctx:    DFT prolate solver context
 * @iter:   Current iteration index
 * @NBiter: Total number of iterations
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE on error.
 */
static errno_t prolate_dft_step(
    prolate_dft_ctx_t *ctx,
    long               iter,
    long               NBiter)
{
    uint64_t    size2  = ctx->size2;
    long double total  = 0.0;
    long double total2 = 0.0;
    imageID     ID;

    printf("[%ld/%ld - %lu] ", iter, NBiter, (unsigned long) size2);
    for (uint64_t ii = 0; ii < size2; ii++)
    {
        dcimg[ctx->IDpupa0m].array.F[ii] =
            dcimg[ctx->IDpupa0].array.F[ii] * dcimg[ctx->IDprol].array.F[ii];
        dcimg[ctx->IDpupp0m].array.F[ii] = dcimg[ctx->IDpupp0].array.F[ii];
    }

    FUNC_CHECK_RETURN(mk_complex_from_amph("pupa0m", "pupp0", "pc1", 0));
    FUNC_CHECK_RETURN(fft_DFTinsertFPM("pc1", "_fpmz", g_coro.DFTZFACTOR, "pc3", NULL));

    for (uint64_t ii = 0; ii < size2; ii++)
    {
        total += dcimg[ctx->IDprol].array.F[ii] * dcimg[ctx->IDprol].array.F[ii];
    }
    ctx->transm = (double) (total / ctx->total0);

    total = dcdoublearr[0] / (1.0 * size2 * ctx->total0);

    ID = image_ID("pc1", dcimg, dcnimg);
    for (uint64_t ii = 0; ii < size2; ii++)
    {
        dcimg[ID].array.CF[ii].re *= ctx->size;
        dcimg[ID].array.CF[ii].im *= ctx->size;
        double re = dcimg[ID].array.CF[ii].re;
        double im = dcimg[ID].array.CF[ii].im;
        total2 += re * re + im * im;
    }
    total2 /= 1.0 * size2 * ctx->total0;

    printf(
        "%ld/%ld   FPmask throughput = %g    prolate throughput = %.18g\n",
        iter,
        NBiter,
        (double) (total / total2),
        ctx->transm);

    FUNC_CHECK_RETURN(mk_reim_from_complex("pc3", "pr2", "pi2", 0));
    prolate_dft_write_and_normalize(ctx, (double) (total / total2));

    FUNC_CHECK_RETURN(delete_image_ID("pc1", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("pc3", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("pr2", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("pi2", DELETE_IMAGE_ERRMODE_WARNING));
    return RETURN_SUCCESS;
}

/**
 * coronagraph_make_2Dprolate_DFT - Compute 2D prolate apodizer via zoomed DFTs
 * @fpmradpix:    Focal plane mask radius [pix]
 * @beamradpix:   Entrance pupil beam radius [pix]
 * @centralObs:   Central obstruction ratio
 * @outname:      Output prolate apodizer image name
 * @size:         2D square grid size [pix]
 * @pupmask_name: Optional extra pupil mask image name
 * @transmval:    Optional output pointer receiving prolate intensity throughput
 *
 * Return: RETURN_SUCCESS (0.0) on success, RETURN_FAILURE on error.
 */
double coronagraph_make_2Dprolate_DFT(
    double      fpmradpix,
    double      beamradpix,
    double      centralObs,
    const char *outname,
    long        size,
    const char *pupmask_name,
    double     *transmval)
{
    DEBUG_TRACE_FSTART();

    prolate_dft_ctx_t ctx = {
        .size       = size,
        .size2      = (uint64_t) size * size,
        .fpmradpix  = fpmradpix,
        .beamradpix = beamradpix,
        .centralObs = centralObs,
    };
    long    NBiter = 10;
    imageID ID;

    if ((ID = variable_ID("DFTZFACTOR")) != -1)
    {
        g_coro.DFTZFACTOR = coffee_dcvar_f(ID);
    }
    if ((ID = variable_ID("PNBITER")) != -1)
    {
        NBiter = (long) (1.0 * coffee_dcvar_f(ID) + 0.01);
    }

    MASKSIZELD = (2.0 * beamradpix / size) * fpmradpix;
    printf("MASK SIZE = %f pixel = %f l/D\n", fpmradpix, MASKSIZELD);
    printf("PUPIL RADIUS = %f\n", beamradpix);
    printf("SIZE = %ld\n", size);

    FUNC_CHECK_RETURN(create_2Dimage_ID("proli", size, size, &ctx.IDproli));
    FUNC_CHECK_RETURN(create_2Dimage_ID("prolr", size, size, &ctx.IDprolr));
    FUNC_CHECK_RETURN(create_2Dimage_ID("prolp", size, size, &ctx.IDprolp));

    FUNC_CHECK_RETURN(prolate_dft_create_fpm(&ctx));
    FUNC_CHECK_RETURN(prolate_dft_init_pupil_and_apo(&ctx, outname, pupmask_name));

    for (long iter = 0; iter < NBiter; iter++)
    {
        FUNC_CHECK_RETURN(prolate_dft_step(&ctx, iter, NBiter));
    }

    for (uint64_t ii = 0; ii < ctx.size2; ii++)
    {
        double mult                     = 1.0 + (1.0 - ctx.peak) / ctx.peak;
        dcimg[ctx.IDfpmask].array.F[ii]  = 1.0 - dcimg[ctx.IDfpmask].array.F[ii] * mult;
        dcimg[ctx.IDfpmaskz].array.F[ii] = 1.0 - dcimg[ctx.IDfpmaskz].array.F[ii] * mult;
    }

    FUNC_CHECK_RETURN(delete_image_ID("FPmaskz", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("FPmask", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("_fpmz", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("pupp0", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("pupa0", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("pupa0m", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("proli", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("prolr", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("prolp", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("pupp0m", DELETE_IMAGE_ERRMODE_WARNING));

    if (transmval != NULL)
    {
        *transmval = ctx.transm;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
