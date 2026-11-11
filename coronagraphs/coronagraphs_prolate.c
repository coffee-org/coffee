/**
 * @file    coronagraphs_prolate.c
 * @brief   FFT-based 2D prolate spheroidal apodizer computation and entry wrappers
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "coffee_compat.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "image_gen/image_gen.h"

#include "coronagraphs_internal.h"

typedef struct prolate_fft_ctx
{
    long        size;
    long        size2;
    double      fpmradpix;
    double      beamradpix;
    double      centralObs;
    long        IDfpmask;
    long        IDpupa0;
    long        IDpupp0;
    long        IDpupa0m;
    long        IDpupp0m;
    long        IDprol;
    long        IDprolr;
    long        IDproli;
    long        IDprolp;
    long double total0;
    double      transm;
    double      peak;
} prolate_fft_ctx_t;

/**
 * prolate_fft_init_pupil - Initialize or validate entrance pupil amplitude "pupa0"
 * @ctx:          FFT prolate solver context
 * @pupmask_name: Optional external pupil mask image name
 */
static void prolate_fft_init_pupil(
    prolate_fft_ctx_t *ctx,
    const char        *pupmask_name)
{
    long size    = ctx->size;
    int  co_flag = (ctx->centralObs > 0.001) ? 1 : 0;

    create_2Dimage_ID("pupp0", size, size, &ctx->IDpupp0);

    ctx->IDpupa0 = image_ID("pupa0", dcimg, dcnimg);
    if (ctx->IDpupa0 != -1)
    {
        if (dcimg[ctx->IDpupa0].md[0].size[0] != (uint32_t) size)
        {
            printf("ERROR: pupa0 should be %ld x %ld\n", size, size);
            exit(0);
        }
        save_fits("pupa0", "test__pupa0.fits");
    }
    else
    {
        ctx->IDpupa0 = make_disk("pupa0", size, size, size / 2, size / 2, ctx->beamradpix);
        if (co_flag == 1)
        {
            long IDco = make_disk(
                "pupa0co",
                size,
                size,
                size / 2,
                size / 2,
                ctx->beamradpix * ctx->centralObs);
            for (long ii = 0; ii < ctx->size2; ii++)
            {
                dcimg[ctx->IDpupa0].array.F[ii] -= dcimg[IDco].array.F[ii];
            }
            delete_image_ID("pupa0co", DELETE_IMAGE_ERRMODE_WARNING);
        }

        long IDpupmask = image_ID(pupmask_name, dcimg, dcnimg);
        if (IDpupmask != -1)
        {
            for (long ii = 0; ii < ctx->size2; ii++)
            {
                dcimg[ctx->IDpupa0].array.F[ii] *= dcimg[IDpupmask].array.F[ii];
            }
        }
        save_fits("pupa0", "test_pupa0.fits");
    }

    ctx->total0 = 0.0;
    for (long ii = 0; ii < ctx->size2; ii++)
    {
        double v1 = dcimg[ctx->IDpupa0].array.F[ii];
        ctx->total0 += v1 * v1;
    }
}

/**
 * prolate_fft_init_apo - Seed initial Gaussian or "apostart" prolate estimate
 * @ctx:     FFT prolate solver context
 * @outname: Output prolate image name
 */
static void prolate_fft_init_apo(
    prolate_fft_ctx_t *ctx,
    const char        *outname)
{
    long size = ctx->size;
    long ID;

    create_2Dimage_ID("pupa0m", size, size, &ctx->IDpupa0m);
    create_2Dimage_ID("pupp0m", size, size, &ctx->IDpupp0m);
    create_2Dimage_ID(outname, size, size, &ctx->IDprol);

    for (long ii = 0; ii < size; ii++)
    {
        for (long jj = 0; jj < size; jj++)
        {
            long idx = jj * size + ii;
            if (dcimg[ctx->IDpupa0].array.F[idx] > 0.01)
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

    ID = image_ID("apostart", dcimg, dcnimg);
    if (ID != -1)
    {
        for (long ii = 0; ii < ctx->size2; ii++)
        {
            dcimg[ctx->IDprolr].array.F[ii] = dcimg[ID].array.F[ii];
            dcimg[ctx->IDprol].array.F[ii]  = dcimg[ID].array.F[ii];
        }
    }
    delete_image_ID("apostart", DELETE_IMAGE_ERRMODE_WARNING);
}

/**
 * prolate_fft_flip_norm - Invert coordinates and normalize 2D FFT output plane
 * @name:   Target float image name
 * @cpname: Temporary copy image name
 * @size:   Grid dimension [pix]
 * @out_id: Output image ID of target
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE.
 */
static errno_t prolate_fft_flip_norm(
    const char *name,
    const char *cpname,
    long        size,
    long       *out_id)
{
    long size2 = size * size;
    long id_im = image_ID(name, dcimg, dcnimg);
    long id_cp;

    if (id_im == -1)
    {
        return RETURN_FAILURE;
    }
    copy_image_ID(name, cpname, 0);
    id_cp = image_ID(cpname, dcimg, dcnimg);
    if (id_cp == -1)
    {
        return RETURN_FAILURE;
    }

    for (long ii = 0; ii < size2; ii++)
    {
        dcimg[id_im].array.F[ii] = 0.0;
    }
    for (long ii = 1; ii < size; ii++)
    {
        for (long jj = 1; jj < size; jj++)
        {
            dcimg[id_im].array.F[jj * size + ii] =
                dcimg[id_cp].array.F[(size - jj) * size + (size - ii)] / size2;
        }
    }

    delete_image_ID(cpname, DELETE_IMAGE_ERRMODE_WARNING);
    *out_id = id_im;
    return RETURN_SUCCESS;
}

/**
 * prolate_fft_write_and_normalize - Log eigenmode metrics and normalize next iterate
 * @ctx:      FFT prolate solver context
 * @IDr:      Real pupil image ID
 * @IDi:      Imaginary pupil image ID
 * @fpm_thru: Fraction of focal energy inside focal plane mask
 */
static void prolate_fft_write_and_normalize(
    prolate_fft_ctx_t *ctx,
    long               IDr,
    long               IDi,
    double             fpm_thru)
{
    long  size  = ctx->size;
    long  size2 = ctx->size2;
    char  fname[1500];
    FILE *fp;

    ctx->peak = dcimg[IDr].array.F[(size / 2) * size + size / 2];
    printf("PEAK =  %f ", ctx->peak);
    ctx->peak = 0.0;
    for (long ii = 0; ii < size2; ii++)
    {
        if ((dcimg[ctx->IDpupa0].array.F[ii] > 0.01) &&
            (fabs(dcimg[IDr].array.F[ii]) > fabs(ctx->peak)))
        {
            ctx->peak = dcimg[IDr].array.F[ii];
        }
    }
    printf(" %f\n", ctx->peak);

    printf("------- MASK SIZE = %f l/D ------\n", MASKSIZELD);
    sprintf(fname, "%s/APLCapo.%.3f.%.3f.info", dcsavedir, MASKSIZELD, ctx->centralObs);
    fp = fopen(fname, "w");
    fprintf(
        fp,
        "%f %.18f %.18f %.18f %.18f\n",
        MASKSIZELD,
        fpm_thru,
        ctx->transm,
        ctx->peak,
        ctx->centralObs);
    fprintf(fp, "# col 1 : Mask size in SYSTEM l/D\n");
    fprintf(fp, "# col 2 : Fraction of light going through FPmask\n");
    fprintf(fp, "# col 3 : Prolate throughput if not done by PIAA\n");
    fprintf(fp, "# col 4 : Peak\n");
    fprintf(fp, "# col 5 : Central obstruction");
    fclose(fp);
    printf("IDEAL COMPLEX TRANSMISSION = %g\n", -(1.0 - ctx->peak) / ctx->peak);

    for (long ii = 0; ii < size2; ii++)
    {
        if (dcimg[ctx->IDpupa0].array.F[ii] > 0.01)
        {
            dcimg[ctx->IDprolr].array.F[ii] = dcimg[IDr].array.F[ii];
            dcimg[ctx->IDproli].array.F[ii] = dcimg[IDi].array.F[ii];
        }
    }

    for (long ii = 1; ii < size - 1; ii++)
    {
        for (long jj = 1; jj < size - 1; jj++)
        {
            long   idx = jj * size + ii;
            double re  = dcimg[ctx->IDprolr].array.F[idx];
            double im  = dcimg[ctx->IDproli].array.F[idx];
            dcimg[ctx->IDprol].array.F[idx]  = sqrt(re * re + im * im);
            dcimg[ctx->IDprolp].array.F[idx] = -atan2(im, re);
        }
    }

    ctx->peak = 0.0;
    for (long ii = 0; ii < size2; ii++)
    {
        if ((dcimg[ctx->IDpupa0].array.F[ii] > 0.01) &&
            (fabs(dcimg[ctx->IDprol].array.F[ii]) > fabs(ctx->peak)))
        {
            ctx->peak = dcimg[ctx->IDprol].array.F[ii];
        }
    }

    for (long ii = 0; ii < size2; ii++)
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
 * prolate_fft_step - Execute one FFT prolate power-method iteration
 * @ctx:    FFT prolate solver context
 * @iter:   Current iteration index
 * @NBiter: Total number of iterations
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE.
 */
static errno_t prolate_fft_step(
    prolate_fft_ctx_t *ctx,
    long               iter,
    long               NBiter)
{
    long        size2 = ctx->size2;
    long double total = 0.0;
    long double total2 = 0.0;
    long        IDfa1m, IDfa1, IDr, IDi;

    for (long ii = 0; ii < size2; ii++)
    {
        dcimg[ctx->IDpupa0m].array.F[ii] =
            dcimg[ctx->IDpupa0].array.F[ii] * dcimg[ctx->IDprol].array.F[ii];
        dcimg[ctx->IDpupp0m].array.F[ii] = dcimg[ctx->IDpupp0].array.F[ii];
    }

    if ((mk_complex_from_amph("pupa0m", "pupp0", "pc1", 0) != RETURN_SUCCESS) ||
        (image_ID("pc1", dcimg, dcnimg) == -1) ||
        (coro_fft2d_centered("pc1", "fc1") != RETURN_SUCCESS) ||
        (mk_amph_from_complex("fc1", "fa1", "fp1", 0) != RETURN_SUCCESS))
    {
        PRINT_ERROR("forward FFT step pupa0m->fa1,fp1 failed");
        return RETURN_FAILURE;
    }

    arith_image_mult("fa1", "FPmask", "fa1m");

    for (long ii = 0; ii < size2; ii++)
    {
        total += dcimg[ctx->IDprol].array.F[ii] * dcimg[ctx->IDprol].array.F[ii];
    }
    ctx->transm = (double) (total / ctx->total0);

    IDfa1m = image_ID("fa1m", dcimg, dcnimg);
    IDfa1  = image_ID("fa1", dcimg, dcnimg);
    if ((IDfa1m == -1) || (IDfa1 == -1) || (image_ID("fp1", dcimg, dcnimg) == -1))
    {
        return RETURN_FAILURE;
    }

    total = 0.0;
    for (long ii = 0; ii < size2; ii++)
    {
        total += dcimg[IDfa1m].array.F[ii] * dcimg[IDfa1m].array.F[ii];
        total2 += dcimg[IDfa1].array.F[ii] * dcimg[IDfa1].array.F[ii];
    }
    total /= 1.0 * size2 * ctx->total0;
    total2 /= 1.0 * size2 * ctx->total0;

    printf(
        "%ld/%ld   FPmask throughput = %g    prolate throughput = %.18g\n",
        iter,
        NBiter,
        (double) (total / total2),
        ctx->transm);

    if ((mk_complex_from_amph("fa1m", "fp1", "fc2", 0) != RETURN_SUCCESS) ||
        (image_ID("fc2", dcimg, dcnimg) == -1) ||
        (coro_fft2d_centered("fc2", "pc2") != RETURN_SUCCESS))
    {
        PRINT_ERROR("inverse FFT step fa1m->pc2 failed");
        return RETURN_FAILURE;
    }

    mk_reim_from_complex("pc2", "pr2", "pi2", 0);
    if ((prolate_fft_flip_norm("pr2", "pr2cp", ctx->size, &IDr) != RETURN_SUCCESS) ||
        (prolate_fft_flip_norm("pi2", "pi2cp", ctx->size, &IDi) != RETURN_SUCCESS))
    {
        return RETURN_FAILURE;
    }

    prolate_fft_write_and_normalize(ctx, IDr, IDi, (double) (total / total2));

    delete_image_ID("pc1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fc1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fa1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fp1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fa1m", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pr2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pi2", DELETE_IMAGE_ERRMODE_WARNING);
    return RETURN_SUCCESS;
}

/**
 * coronagraph_make_2Dprolate - Compute 2D prolate spheroidal pupil apodizer via FFTs
 * @fpmradpix:    Focal plane mask radius [pix]
 * @beamradpix:   Entrance pupil beam radius [pix]
 * @centralObs:   Central obstruction ratio (0..1)
 * @outname:      Output prolate apodizer image name
 * @size:         2D square grid size [pix]
 * @pupmask_name: Optional extra pupil mask image name
 * @transmval:    Optional output pointer receiving prolate intensity throughput
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE on error.
 */
errno_t coronagraph_make_2Dprolate(
    double      fpmradpix,
    double      beamradpix,
    double      centralObs,
    const char *outname,
    long        size,
    const char *pupmask_name,
    double     *transmval)
{
    DEBUG_TRACE_FSTART();

    prolate_fft_ctx_t ctx = {
        .size       = size,
        .size2      = size * size,
        .fpmradpix  = fpmradpix,
        .beamradpix = beamradpix,
        .centralObs = centralObs,
    };
    long NBiter = 40;
    long ID     = variable_ID("PNBITER");

    if (ID != -1)
    {
        NBiter = (long) (1.0 * coffee_dcvar_f(ID) + 0.01);
    }

    MASKSIZELD = (2.0 * beamradpix / size) * fpmradpix;
    printf("MASK RADIUS = %f pix = %f l/D\n", fpmradpix, MASKSIZELD);
    printf("PUPIL RADIUS = %f\n", beamradpix);
    printf("SIZE = %ld\n", size);

    create_2Dimage_ID("proli", size, size, &ctx.IDproli);
    create_2Dimage_ID("prolr", size, size, &ctx.IDprolr);
    create_2Dimage_ID("prolp", size, size, &ctx.IDprolp);

    ctx.IDfpmask = make_subpixdisk("FPmask", size, size, size / 2, size / 2, fpmradpix);
    prolate_fft_init_pupil(&ctx, pupmask_name);
    prolate_fft_init_apo(&ctx, outname);

    for (long iter = 0; iter < NBiter; iter++)
    {
        if (prolate_fft_step(&ctx, iter, NBiter) != RETURN_SUCCESS)
        {
            return RETURN_FAILURE;
        }
    }

    for (long ii = 0; ii < ctx.size2; ii++)
    {
        dcimg[ctx.IDfpmask].array.F[ii] =
            1.0 - dcimg[ctx.IDfpmask].array.F[ii] * (1.0 + (1.0 - ctx.peak) / ctx.peak);
    }

    delete_image_ID("FPmask", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pupp0", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pupa0", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pupa0m", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("proli", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("prolr", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("prolp", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pupp0m", DELETE_IMAGE_ERRMODE_WARNING);

    if (transmval != NULL)
    {
        *transmval = ctx.transm;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * coronagraph_make_2Dprolateld - Compute 2D prolate apodizer for mask radius in lambda/D
 * @masksizeld:   Focal plane mask radius [lambda/D]
 * @beamradpix:   Entrance pupil beam radius [pix]
 * @centralObs:   Central obstruction ratio
 * @outname:      Output prolate apodizer image name
 * @size:         2D square grid size [pix]
 * @pupmask_name: Optional extra pupil mask image name
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE on error.
 */
errno_t coronagraph_make_2Dprolateld(
    double      masksizeld,
    double      beamradpix,
    double      centralObs,
    const char *outname,
    long        size,
    const char *pupmask_name)
{
    DEBUG_TRACE_FSTART();

    double fpmradpix = masksizeld / (2.0 * beamradpix / size);

    if (g_coro.useDFT == 0)
    {
        if (coronagraph_make_2Dprolate(
                fpmradpix,
                beamradpix,
                centralObs,
                outname,
                size,
                pupmask_name,
                NULL) != RETURN_SUCCESS)
        {
            FUNC_RETURN_FAILURE("Call to coronagraph_make_2Dprolate failed");
        }
    }
    else
    {
        if (coronagraph_make_2Dprolate_DFT(
                fpmradpix,
                beamradpix,
                centralObs,
                outname,
                size,
                pupmask_name,
                NULL) != RETURN_SUCCESS)
        {
            FUNC_RETURN_FAILURE("Call to coronagraph_make_2Dprolate_DFT failed");
        }
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * coronagraph_make_2Dprolate_CS - Compute 2D prolate apodizer from custom pupil "pupaCS"
 * @masksize:   Focal plane mask radius [pix]
 * @beamradpix: Pupil beam radius [pix]
 * @outname:    Unused legacy output name (result stored in "outapo")
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE on error.
 */
errno_t coronagraph_make_2Dprolate_CS(
    double      masksize,
    double      beamradpix,
    const char *outname)
{
    DEBUG_TRACE_FSTART();

    long size = CORONAGRAPHS_ARRAYSIZE;

    (void) outname;
    copy_image_ID("pupaCS", "pupa0", 0);

    if (coronagraph_make_2Dprolate_DFT(
            masksize,
            beamradpix,
            0.0,
            "outapo",
            size,
            "NULLim",
            NULL) != RETURN_SUCCESS)
    {
        FUNC_RETURN_FAILURE("Call to coronagraph_make_2Dprolate_DFT failed");
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
