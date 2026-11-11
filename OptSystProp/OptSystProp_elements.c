/**
 * @file    OptSystProp_elements.c
 * @brief   Optical element phase, amplitude, and zoomed FPM DFT application kernels
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "OpticsMaterials/OpticsMaterials.h"
#include "fft/fft.h"

#include "OptSystProp/OptSystProp.h"
#include "OptSystProp_internal.h"

/**
 * OptSystProp_apply_ampmask - Multiply wavefront amplitude cube by an amplitude mask.
 * @ctx: Per-element propagation context
 *
 * Supports both 2D achromatic masks (broadcast across spectral channels) and
 * 3D chromatic amplitude cubes matching `nblambda`.
 *
 * Return: RETURN_SUCCESS on success.
 */
errno_t OptSystProp_apply_ampmask(const OptSystProp_elem_ctx *ctx)
{
    OPTSYST *optsyst  = ctx->optsyst;
    long     index    = ctx->index;
    long     elem     = ctx->elem;
    long     nblambda = optsyst[index].nblambda;
    uint64_t size2    = (uint64_t) optsyst[index].size * optsyst[index].size;
    imageID  ID       = optsyst[index].elemarrayindex[elem];

    printf("============= elem %ld:  Opaque mask (%s) =================\n",
           elem,
           dcimg[ID].name);
    fflush(stdout);

    if (ID == -1)
    {
        printf("ERROR: ID = -1, missing mask image\n");
        exit(0);
    }

    if ((dcimg[ID].md[0].naxis == 2) || (dcimg[ID].md[0].size[2] != nblambda))
    {
#ifdef HAVE_LIBGOMP
#pragma omp parallel default(shared)
        {
#pragma omp for
#endif
            for (long kl = 0; kl < nblambda; kl++)
            {
                for (uint64_t ii = 0; ii < size2; ii++)
                {
                    dcimg[ctx->IDa].array.F[size2 * kl + ii] *= dcimg[ID].array.F[ii];
                }
            }
#ifdef HAVE_LIBGOMP
        }
#endif
    }
    else
    {
#ifdef HAVE_LIBGOMP
#pragma omp parallel
        {
#pragma omp for
#endif
            for (uint64_t ii = 0; ii < size2 * nblambda; ii++)
            {
                dcimg[ctx->IDa].array.F[ii] *= dcimg[ID].array.F[ii];
            }
#ifdef HAVE_LIBGOMP
        }
#endif
    }

    return RETURN_SUCCESS;
}

/**
 * OptSystProp_apply_mirror - Apply reflective mirror surface sag phase across wavelengths.
 * @ctx: Per-element propagation context
 *
 * Converts 2D or 3D mirror surface sag `z` [m] into reflected phase `-4 * pi * z / lambda`.
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE if surface image is missing.
 */
errno_t OptSystProp_apply_mirror(const OptSystProp_elem_ctx *ctx)
{
    OPTSYST *optsyst  = ctx->optsyst;
    long     index    = ctx->index;
    long     elem     = ctx->elem;
    long     nblambda = optsyst[index].nblambda;
    uint64_t size2    = (uint64_t) optsyst[index].size * optsyst[index].size;
    int      arridx   = optsyst[index].elemarrayindex[elem];
    imageID  ID       = optsyst[index].ASPHSURFMarray[arridx].surfID;

    printf("============= Mirror surface =======================\n");
    fflush(stdout);
    printf("%d surface ID = %ld\n", arridx, ID);

    if (ID == -1)
    {
        printf("ERROR: Surface ID does not exist\n");
        printf("   %s %d\n", __FILE__, __LINE__);
        printf("   Function %s return FAILURE\n", __func__);
        return RETURN_FAILURE;
    }

    if (dcimg[ID].md[0].naxis == 2)
    {
#ifdef HAVE_LIBGOMP
#pragma omp parallel default(shared)
        {
#pragma omp for
#endif
            for (long kl = 0; kl < nblambda; kl++)
            {
                for (uint64_t ii = 0; ii < size2; ii++)
                {
                    dcimg[ctx->IDp].array.F[size2 * kl + ii] -=
                        4.0 * M_PI * dcimg[ID].array.F[ii] / optsyst[index].lambdaarray[kl];
                }
            }
#ifdef HAVE_LIBGOMP
        }
#endif
    }
    else
    {
#ifdef HAVE_LIBGOMP
#pragma omp parallel default(shared)
        {
#pragma omp for
#endif
            for (long kl = 0; kl < nblambda; kl++)
            {
                for (uint64_t ii = 0; ii < size2; ii++)
                {
                    dcimg[ctx->IDp].array.F[size2 * kl + ii] -=
                        4.0 * M_PI * dcimg[ID].array.F[size2 * kl + ii] /
                        optsyst[index].lambdaarray[kl];
                }
            }
#ifdef HAVE_LIBGOMP
        }
#endif
    }

    return RETURN_SUCCESS;
}

/**
 * OptSystProp_init_refractive_ncoeff - Populate dispersive sag-to-phase coefficients.
 * @optsyst: Optical system array
 * @index:   Active optical system index
 * @elem:    Current element index
 */
static void OptSystProp_init_refractive_ncoeff(
    OPTSYST *optsyst,
    long     index,
    long     elem)
{
    int        arridx = optsyst[index].elemarrayindex[elem];
    ASPHSURFR *surfr  = &optsyst[index].ASPHSURFRarray[arridx];
    long       nblam  = optsyst[index].nblambda;

    for (long kl = 0; kl < nblam; kl++)
    {
        double lam = optsyst[index].lambdaarray[kl];
        double n0  = OpticsMaterials_n(surfr->mat0, lam);
        double n1  = OpticsMaterials_n(surfr->mat1, lam);

        surfr->ncoeff[kl] = 2.0 * M_PI * (n0 - n1) / lam;
        DEBUG_TRACEPOINT("elem %ld ASPHSURFRarray %f %f -> ncoeff %ld = %lf",
                         elem,
                         n0,
                         n1,
                         kl,
                         surfr->ncoeff[kl]);
    }
    surfr->init = 1;
}

/**
 * OptSystProp_apply_refractive - Apply refractive surface dispersive OPD phase cube.
 * @ctx: Per-element propagation context
 *
 * Lazily initializes per-wavelength material refractive index coefficients via
 * `OpticsMaterials_n` and accumulates phase across all spectral slices.
 *
 * Return: RETURN_SUCCESS or error code if surface image is missing.
 */
errno_t OptSystProp_apply_refractive(const OptSystProp_elem_ctx *ctx)
{
    DEBUG_TRACE_FSTART();

    OPTSYST   *optsyst  = ctx->optsyst;
    long       index    = ctx->index;
    long       elem     = ctx->elem;
    long       nblambda = optsyst[index].nblambda;
    uint64_t   size2    = (uint64_t) optsyst[index].size * optsyst[index].size;
    int        arridx   = optsyst[index].elemarrayindex[elem];
    ASPHSURFR *surfr    = &optsyst[index].ASPHSURFRarray[arridx];
    imageID    ID       = surfr->surfID;

    printf("============= [%ld] Refractive surface =======================\n", elem);
    fflush(stdout);

    DEBUG_TRACEPOINT("index %ld  elem %ld  %d surface ID : %ld", index, elem, arridx, ID);
    DEBUG_TRACEPOINT("refractive surface ID %ld %s", (long) ID, dcimg[ID].md[0].name);

    if (ID == -1)
    {
        list_image_ID();
        FUNC_RETURN_FAILURE("image ID not found");
    }

    if (surfr->init != 1)
    {
        OptSystProp_init_refractive_ncoeff(optsyst, index, elem);
    }

    if (dcimg[ID].md[0].naxis == 2)
    {
#ifdef HAVE_LIBGOMP
#pragma omp parallel default(shared)
        {
#pragma omp for
#endif
            for (long kl = 0; kl < nblambda; kl++)
            {
                for (uint64_t ii = 0; ii < size2; ii++)
                {
                    dcimg[ctx->IDp].array.F[size2 * kl + ii] +=
                        dcimg[ID].array.F[ii] * surfr->ncoeff[kl];
                }
            }
#ifdef HAVE_LIBGOMP
        }
#endif
    }
    else
    {
#ifdef HAVE_LIBGOMP
#pragma omp parallel default(shared)
        {
#pragma omp for
#endif
            for (long kl = 0; kl < nblambda; kl++)
            {
                for (uint64_t ii = 0; ii < size2; ii++)
                {
                    dcimg[ctx->IDp].array.F[size2 * kl + ii] +=
                        dcimg[ID].array.F[size2 * kl + ii] * surfr->ncoeff[kl];
                }
            }
#ifdef HAVE_LIBGOMP
        }
#endif
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * OptSystProp_focmask_bin_sparse - Bin complex pupil wavefront onto sparse supergrid.
 * @size:     2D wavefront linear dimension [pix]
 * @nblambda: Number of wavelength channels
 * @pad:      Sub-sampling half-stride (`DFTgridpad`)
 */
static void OptSystProp_focmask_bin_sparse(
    long size,
    long nblambda,
    long pad)
{
    uint64_t size2  = (uint64_t) size * size;
    long     gsize  = 2 * pad + 1;
    long     offset = pad;
    imageID  IDre;
    imageID  IDim;

    create_3Dimage_ID("dftgridre", size, size, nblambda, &IDre);
    create_3Dimage_ID("dftgridim", size, size, nblambda, &IDim);

    imageID ID = image_ID("_WFctmp", dcimg, dcnimg);
    for (long kl = 0; kl < nblambda; kl++)
    {
        for (uint32_t ii = 0; ii < (uint32_t) size; ii++)
        {
            for (uint32_t jj = 0; jj < (uint32_t) size; jj++)
            {
                float re  = dcimg[ID].array.CF[size2 * kl + jj * size + ii].re;
                float im  = dcimg[ID].array.CF[size2 * kl + jj * size + ii].im;
                long  ii1 = offset + ((long) (ii / gsize)) * gsize;
                long  jj1 = offset + ((long) (jj / gsize)) * gsize;

                if ((ii1 < size) && (jj1 < size))
                {
                    dcimg[IDre].array.F[size2 * kl + jj1 * size + ii1] += re;
                    dcimg[IDim].array.F[size2 * kl + jj1 * size + ii1] += im;
                }
            }
        }
    }

    mk_complex_from_reim("dftgridre", "dftgridim", "_WFctmpc", 0);
    delete_image_ID("dftgridre", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("dftgridim", DELETE_IMAGE_ERRMODE_WARNING);
}

/**
 * OptSystProp_focmask_interp_sparse - Bilinearly interpolate sparse DFT output to full grid.
 * @size:     2D wavefront linear dimension [pix]
 * @nblambda: Number of wavelength channels
 * @pad:      Sub-sampling half-stride (`DFTgridpad`)
 */
static void OptSystProp_focmask_interp_sparse(
    long size,
    long nblambda,
    long pad)
{
    uint64_t size2  = (uint64_t) size * size;
    long     gsize  = 2 * pad + 1;
    long     offset = pad;
    long     kwidth = 2 * gsize + 1;
    float   *convkern;

    convkern = (float *) malloc(sizeof(float) * kwidth * kwidth);
    if (convkern == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    double tot = 0.0;
    for (long i = 0; i < kwidth; i++)
    {
        for (long j = 0; j < kwidth; j++)
        {
            float u              = fabs(1.0 * (i - gsize) / gsize);
            float t              = fabs(1.0 * (j - gsize) / gsize);
            float val            = (1.0 - u) * (1.0 - t);
            convkern[j * kwidth + i] = val;
            tot += val;
        }
    }
    for (long i = 0; i < kwidth * kwidth; i++)
    {
        convkern[i] *= gsize * gsize / tot;
    }

    imageID ID = image_ID("_WFcout", dcimg, dcnimg);
    imageID IDre1;
    imageID IDim1;
    create_3Dimage_ID("dftgridre1", size, size, nblambda, &IDre1);
    create_3Dimage_ID("dftgridim1", size, size, nblambda, &IDim1);

    for (long kl = 0; kl < nblambda; kl++)
    {
        for (long ii1 = offset + gsize; ii1 < size - gsize; ii1 += gsize)
        {
            for (long jj1 = offset + gsize; jj1 < size - gsize; jj1 += gsize)
            {
                float re = dcimg[ID].array.CF[size2 * kl + jj1 * size + ii1].re;
                float im = dcimg[ID].array.CF[size2 * kl + jj1 * size + ii1].im;

                for (long i = 0; i < kwidth; i++)
                {
                    for (long j = 0; j < kwidth; j++)
                    {
                        long  ii = ii1 + (i - gsize);
                        long  jj = jj1 + (j - gsize);
                        float kw = convkern[j * kwidth + i];
                        dcimg[IDre1].array.F[size2 * kl + jj * size + ii] += re * kw;
                        dcimg[IDim1].array.F[size2 * kl + jj * size + ii] += im * kw;
                    }
                }
            }
        }
    }

    free(convkern);
    delete_image_ID("_WFcout", DELETE_IMAGE_ERRMODE_WARNING);
    mk_complex_from_reim("dftgridre1", "dftgridim1", "_WFcout", 0);
    delete_image_ID("dftgridre1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("dftgridim1", DELETE_IMAGE_ERRMODE_WARNING);
}

/**
 * OptSystProp_focmask_sparse_dft - Execute sub-sampled pupil-to-FPM-to-pupil DFT pass.
 * @ctx: Per-element propagation context
 */
static void OptSystProp_focmask_sparse_dft(const OptSystProp_elem_ctx *ctx)
{
    OPTSYST *optsyst   = ctx->optsyst;
    long     index     = ctx->index;
    long     elemindex = optsyst[index].elemarrayindex[ctx->elem];
    imageID  ID        = optsyst[index].FOCMASKarray[elemindex].fpmID;

    OptSystProp_focmask_bin_sparse(optsyst[index].size,
                                   optsyst[index].nblambda,
                                   optsyst[index].DFTgridpad);

    printf("focm # %ld: %s\n", elemindex, dcimg[ID].name);

    delete_image_ID("fpmatest", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fpmptest", DELETE_IMAGE_ERRMODE_WARNING);
    mk_amph_from_complex("piaacmcfpm", "fpmatest", "fpmptest", 0);

    if (optsyst[index].SAVE == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];

        mk_amph_from_complex(dcimg[ID].name, "fpma", "fpmp", 0);
        WRITE_FULLFILENAME(fname, "%s/fpm__ampl.fits", ctx->savedir);
        save_fits("fpma", fname);
        WRITE_FULLFILENAME(fname, "%s/fpm__pha.fits", ctx->savedir);
        save_fits("fpmp", fname);
        delete_image_ID("fpma", DELETE_IMAGE_ERRMODE_WARNING);
        delete_image_ID("fpmp", DELETE_IMAGE_ERRMODE_WARNING);
    }

    fft_DFTinsertFPM("_WFctmpc",
                     dcimg[ID].name,
                     optsyst[index].FOCMASKarray[elemindex].zfactor,
                     "_WFcout",
                     NULL);
    delete_image_ID("_WFctmpc", DELETE_IMAGE_ERRMODE_WARNING);

    OptSystProp_focmask_interp_sparse(optsyst[index].size,
                                      optsyst[index].nblambda,
                                      optsyst[index].DFTgridpad);
}

/**
 * OptSystProp_apply_focmask - Apply focal plane mask via zoomed 2D DFT (direct or Babinet).
 * @ctx: Per-element propagation context
 *
 * Converts pupil amplitude/phase to complex `_WFctmp`, evaluates `fft_DFTinsertFPM`
 * (optionally on a sub-sampled grid when `DFTgridpad > 0`), and reconstructs output
 * amplitude/phase cubes using Babinet subtraction (`mode == 1`) or direct (`mode == 0`).
 *
 * Return: RETURN_SUCCESS on success.
 */
errno_t OptSystProp_apply_focmask(const OptSystProp_elem_ctx *ctx)
{
    OPTSYST *optsyst   = ctx->optsyst;
    long     index     = ctx->index;
    long     elem      = ctx->elem;
    long     elemindex = optsyst[index].elemarrayindex[elem];

    printf("============= [%ld] Focal Plane Mask ==============\n", elem);
    fflush(stdout);

    mk_complex_from_amph(ctx->imnameamp_out, ctx->imnamepha_out, "_WFctmp", 0);
    if (ctx->sharedmem == 0)
    {
        delete_image_ID(ctx->imnameamp_out, DELETE_IMAGE_ERRMODE_WARNING);
        delete_image_ID(ctx->imnamepha_out, DELETE_IMAGE_ERRMODE_WARNING);
    }

    printf("optsyst[index].DFTgridpad = %ld\n", optsyst[index].DFTgridpad);

    if (optsyst[index].DFTgridpad > 0)
    {
        OptSystProp_focmask_sparse_dft(ctx);
    }
    else
    {
        imageID ID = optsyst[index].FOCMASKarray[elemindex].fpmID;
        printf("focm : %s\n", dcimg[ID].name);
        fflush(stdout);

        fft_DFTinsertFPM("_WFctmp",
                         dcimg[ID].name,
                         optsyst[index].FOCMASKarray[elemindex].zfactor,
                         "_WFcout",
                         NULL);
    }

    printf("optsyst[index].FOCMASKarray[elemindex].mode = %d\n",
           optsyst[index].FOCMASKarray[elemindex].mode);
    if (optsyst[index].FOCMASKarray[elemindex].mode == 1)
    {
        arith_image_sub_inplace("_WFctmp", "_WFcout");
        mk_amph_from_complex("_WFctmp", ctx->imnameamp_out, ctx->imnamepha_out, 0);
    }
    else
    {
        mk_amph_from_complex("_WFcout", ctx->imnameamp_out, ctx->imnamepha_out, 0);
    }

    delete_image_ID("_WFctmp", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("_WFcout", DELETE_IMAGE_ERRMODE_WARNING);

    return RETURN_SUCCESS;
}

/**
 * OptSystProp_apply_element - Dispatch optical element transfer function by elemtype.
 * @ctx: Per-element propagation context
 *
 * Return: RETURN_SUCCESS or element-specific failure status.
 */
errno_t OptSystProp_apply_element(const OptSystProp_elem_ctx *ctx)
{
    int etype = ctx->optsyst[ctx->index].elemtype[ctx->elem];

    switch (etype)
    {
    case 1:
        return OptSystProp_apply_ampmask(ctx);
    case 3:
        return OptSystProp_apply_mirror(ctx);
    case 4:
        return OptSystProp_apply_refractive(ctx);
    case 5:
        return OptSystProp_apply_focmask(ctx);
    default:
        return RETURN_SUCCESS;
    }
}
