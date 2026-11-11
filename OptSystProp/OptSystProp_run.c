/**
 * @file    OptSystProp_run.c
 * @brief   Sequential multi-element Fresnel optical propagation and final PSF orchestrator
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
#include "fft/fft.h"

#include "OptSystProp/OptSystProp.h"
#include "OptSystProp_internal.h"

/**
 * OptSystProp_init_base_cubes - Allocate and initialize entrance wavefront amp/phase cubes.
 * @optsyst:   Optical system array
 * @index:     Active optical system index
 * @sharedmem: Shared-memory allocation flag
 */
static void OptSystProp_init_base_cubes(
    OPTSYST *optsyst,
    long     index,
    int      sharedmem)
{
    long      size     = optsyst[index].size;
    long      nblambda = optsyst[index].nblambda;
    uint64_t  size2    = (uint64_t) size * size;
    uint32_t *imsizearray;

    imsizearray = (uint32_t *) malloc(sizeof(uint32_t) * 3);
    if (imsizearray == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    imsizearray[0] = size;
    imsizearray[1] = size;
    imsizearray[2] = nblambda;

    {
        char imname[STRINGMAXLEN_IMGNAME];
        WRITE_IMAGENAME(imname, "WFamp%ld", index);
        imageID IDa = image_ID(imname, dcimg, dcnimg);
        if (IDa == -1)
        {
            create_image_ID(imname, 3, imsizearray, _DATATYPE_FLOAT, sharedmem, 0, 0, &IDa);
        }

        for (uint64_t ii = 0; ii < size2; ii++)
        {
            for (long kl = 0; kl < nblambda; kl++)
            {
                dcimg[IDa].array.F[size2 * kl + ii] = 1.0;
            }
        }
    }

    {
        char imname[STRINGMAXLEN_IMGNAME];
        WRITE_IMAGENAME(imname, "WFpha%ld", index);
        imageID IDp = image_ID(imname, dcimg, dcnimg);
        if (IDp == -1)
        {
            create_image_ID(imname, 3, imsizearray, _DATATYPE_FLOAT, sharedmem, 0, 0, NULL);
        }
    }

    free(imsizearray);
}

/**
 * OptSystProp_find_start_elem - Locate earliest cached element up to requested elemstart.
 * @index:     Active optical system index
 * @elemstart: Requested starting element index
 *
 * Return: Effective starting element index `elemstart1`.
 */
static long OptSystProp_find_start_elem(
    long index,
    long elemstart)
{
    char imnameamp_in[STRINGMAXLEN_IMGNAME];
    char imnamepha_in[STRINGMAXLEN_IMGNAME];
    long elemstart1 = 0;
    int  elemOK     = 1;

    while (elemOK == 1)
    {
        long ID1 = -1;
        long ID2 = -1;

        if (elemstart1 == 0)
        {
            WRITE_IMAGENAME(imnameamp_in, "WFamp%ld", index);
            WRITE_IMAGENAME(imnamepha_in, "WFpha%ld", index);
        }
        else
        {
            WRITE_IMAGENAME(imnameamp_in, "WFamp%ld_%03ld", index, elemstart1 - 1);
            WRITE_IMAGENAME(imnamepha_in, "WFpha%ld_%03ld", index, elemstart1 - 1);
        }

        ID1 = image_ID(imnameamp_in, dcimg, dcnimg);
        if ((ID1 != -1) &&
            ((ID2 = image_ID(imnamepha_in, dcimg, dcnimg)) != -1) &&
            (elemstart1 < elemstart + 1))
        {
            elemstart1++;
            elemOK = 1;
        }
        else
        {
            elemOK = 0;
        }

        printf("%ld/%ld %d    %s %ld   %s %ld\n",
               elemstart1,
               elemstart,
               elemOK,
               imnameamp_in,
               ID1,
               imnamepha_in,
               ID2);
    }

    elemstart1--;
    printf("STARTING AT ELEMENT %ld\n", elemstart1);
    return elemstart1;
}

/**
 * OptSystProp_propagate_to_elem - Propagate or copy wavefront from previous plane to elem.
 * @optsyst:       Optical system array
 * @index:         Active optical system index
 * @elem:          Destination element index
 * @sharedmem:     Shared-memory retention flag
 * @imnameamp_out: Output buffer receiving destination amplitude cube name
 * @imnamepha_out: Output buffer receiving destination phase cube name
 */
static void OptSystProp_propagate_to_elem(
    OPTSYST *optsyst,
    long     index,
    long     elem,
    int      sharedmem,
    char    *imnameamp_out,
    char    *imnamepha_out)
{
    const double proplim = 1.0e-4;
    char         imnameamp_in[STRINGMAXLEN_IMGNAME];
    char         imnamepha_in[STRINGMAXLEN_IMGNAME];
    double       propdist;

    if (elem == 0)
    {
        WRITE_IMAGENAME(imnameamp_in, "WFamp%ld", index);
        WRITE_IMAGENAME(imnamepha_in, "WFpha%ld", index);
        WRITE_IMAGENAME(imnameamp_out, "WFamp%ld_000", index);
        WRITE_IMAGENAME(imnamepha_out, "WFpha%ld_000", index);
        propdist = optsyst[index].elemZpos[0];
    }
    else
    {
        WRITE_IMAGENAME(imnameamp_in, "WFamp%ld_%03ld", index, elem - 1);
        WRITE_IMAGENAME(imnamepha_in, "WFpha%ld_%03ld", index, elem - 1);
        WRITE_IMAGENAME(imnameamp_out, "WFamp%ld_%03ld", index, elem);
        WRITE_IMAGENAME(imnamepha_out, "WFpha%ld_%03ld", index, elem);
        propdist = optsyst[index].elemZpos[elem] - optsyst[index].elemZpos[elem - 1];
    }

    if ((image_ID(imnameamp_out, dcimg, dcnimg) != -1) && (sharedmem == 0))
    {
        delete_image_ID(imnameamp_out, DELETE_IMAGE_ERRMODE_WARNING);
    }
    if ((image_ID(imnamepha_out, dcimg, dcnimg) != -1) && (sharedmem == 0))
    {
        delete_image_ID(imnamepha_out, DELETE_IMAGE_ERRMODE_WARNING);
    }

    if (fabs(propdist) > proplim)
    {
        printf("Propagating to element %ld  (%lf m)\n", elem, propdist);
        OptSystProp_propagateCube(optsyst,
                                  0,
                                  imnameamp_in,
                                  imnamepha_in,
                                  imnameamp_out,
                                  imnamepha_out,
                                  propdist,
                                  sharedmem);
    }
    else
    {
        copy_image_ID(imnameamp_in, imnameamp_out, sharedmem);
        copy_image_ID(imnamepha_in, imnamepha_out, sharedmem);
    }

    printf("*********** %ld  -> %d\n", elem - 1, optsyst[index].keepMem[elem - 1]);
    if ((optsyst[index].keepMem[elem - 1] == 0) && (sharedmem == 0))
    {
        printf("********** Deleting element %ld      %s %s\n",
               elem - 1,
               imnameamp_in,
               imnamepha_in);
        delete_image_ID(imnameamp_in, DELETE_IMAGE_ERRMODE_WARNING);
        delete_image_ID(imnamepha_in, DELETE_IMAGE_ERRMODE_WARNING);
    }
}

/**
 * OptSystProp_record_elem_plane - Compute post-element integrated flux and optional FITS dump.
 * @ctx: Per-element propagation context
 */
static void OptSystProp_record_elem_plane(const OptSystProp_elem_ctx *ctx)
{
    OPTSYST *optsyst  = ctx->optsyst;
    long     index    = ctx->index;
    long     elem     = ctx->elem;
    long     nblambda = optsyst[index].nblambda;
    uint64_t size2    = (uint64_t) optsyst[index].size * optsyst[index].size;
    imageID  IDa      = image_ID(ctx->imnameamp_out, dcimg, dcnimg);

    optsyst[index].flux[elem] = 0.0;
    for (long kl = 0; kl < nblambda; kl++)
    {
        for (uint64_t ii = 0; ii < size2; ii++)
        {
            float amp = dcimg[IDa].array.F[kl * size2 + ii];
            optsyst[index].flux[elem] += amp * amp;
        }
    }

    printf("Element %ld  [%ld %ld]  Flux = %lf\n",
           elem,
           nblambda,
           (long) size2,
           optsyst[index].flux[elem] / nblambda);

    if (isnan(optsyst[index].flux[elem]) != 0)
    {
        exit(0);
    }

    if (optsyst[index].SAVE == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];

        printf("Saving intermediate plane [%ld] ... ", elem);
        WRITE_FULLFILENAME(fname, "./%s/WFamp%ld_%03ld.fits", ctx->savedir, index, elem);
        save_fits(ctx->imnameamp_out, fname);
        WRITE_FULLFILENAME(fname, "./%s/WFpha%ld_%03ld.fits", ctx->savedir, index, elem);
        save_fits(ctx->imnamepha_out, fname);
        printf("done\n");
        fflush(stdout);
    }
}

/**
 * OptSystProp_compute_final_psf - Compute normalized focal-plane FFT PSF from exit pupil.
 * @optsyst:       Optical system array
 * @index:         Active optical system index
 * @imnameamp_out: Exit-pupil amplitude cube image name
 * @imnamepha_out: Exit-pupil phase cube image name
 * @savedir:       Diagnostic output directory
 * @sharedmem:     Shared-memory retention flag
 */
static void OptSystProp_compute_final_psf(
    OPTSYST    *optsyst,
    long        index,
    const char *imnameamp_out,
    const char *imnamepha_out,
    const char *savedir,
    int         sharedmem)
{
    long     nblambda = optsyst[index].nblambda;
    uint64_t size2    = (uint64_t) optsyst[index].size * optsyst[index].size;
    char     imname[STRINGMAXLEN_IMGNAME];
    char     imnameamp[STRINGMAXLEN_IMGNAME];
    char     imnamepha[STRINGMAXLEN_IMGNAME];
    char     imnamere[STRINGMAXLEN_IMGNAME];
    char     imnameim[STRINGMAXLEN_IMGNAME];

    printf("COMPUTING FINAL IMAGE AS FFT OF %ld\n", optsyst[index].NBelem - 1);
    mk_complex_from_amph(imnameamp_out, imnamepha_out, "_WFctmp", 0);
    permut("_WFctmp");

    WRITE_IMAGENAME(imname, "psfc%ld", index);
    do2dfft("_WFctmp", imname);
    delete_image_ID("_WFctmp", DELETE_IMAGE_ERRMODE_WARNING);
    permut(imname);

    WRITE_IMAGENAME(imnameamp, "psfa%ld", index);
    WRITE_IMAGENAME(imnamepha, "psfp%ld", index);
    WRITE_IMAGENAME(imnamere, "psfre%ld", index);
    WRITE_IMAGENAME(imnameim, "psfim%ld", index);

    mk_reim_from_complex(imname, imnamere, imnameim, sharedmem);
    mk_amph_from_complex(imname, imnameamp, imnamepha, sharedmem);

    if (optsyst[index].SAVE == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];

        WRITE_FULLFILENAME(fname, "%s/psfa%ld.fits", savedir, index);
        save_fits(imnameamp, fname);
        WRITE_FULLFILENAME(fname, "%s/psfp%ld.fits", savedir, index);
        save_fits(imnamepha, fname);
        WRITE_FULLFILENAME(fname, "%s/psfre%ld.fits", savedir, index);
        save_fits(imnamere, fname);
        WRITE_FULLFILENAME(fname, "%s/psfim%ld.fits", savedir, index);
        save_fits(imnameim, fname);
    }

    {
        imageID ID   = image_ID(imnameamp, dcimg, dcnimg);
        imageID IDre = image_ID(imnamere, dcimg, dcnimg);
        imageID IDim = image_ID(imnameim, dcimg, dcnimg);
        double  norm = sqrt(size2 * optsyst[index].flux[0] / nblambda);

        for (uint64_t ii = 0; ii < size2 * nblambda; ii++)
        {
            dcimg[ID].array.F[ii] /= norm;
            dcimg[IDre].array.F[ii] /= norm;
            dcimg[IDim].array.F[ii] /= norm;
        }
    }

    WRITE_IMAGENAME(imname, "psfi%ld", index);
    arith_image_mult(imnameamp, imnameamp, imname);

    double total = arith_image_total(imname) / nblambda;
    printf("TOTAL = %lf\n", total);

    if (optsyst[index].SAVE == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];

        WRITE_FULLFILENAME(fname, "%s/psfi%ld.fits", savedir, index);
        save_fits(imname, fname);
    }
}

/**
 * OptSystProp_run - Propagate complex wavefront sequentially across optical elements.
 * @optsyst:   Optical system model array
 * @index:     Optical system index (nominally 0)
 * @elemstart: First optical element index to evaluate
 * @elemend:   One-past-last optical element index to evaluate
 * @savedir:   Output directory for diagnostic FITS files when `SAVE == 1`
 * @sharedmem: Retain intermediate/output wavefront cubes in shared memory if non-zero
 *
 * Return: RETURN_SUCCESS on completion, or element error status on failure.
 */
errno_t OptSystProp_run(
    OPTSYST    *optsyst,
    long        index,
    long        elemstart,
    long        elemend,
    const char *savedir,
    int         sharedmem)
{
    DEBUG_TRACE_FSTART();

    char imnameamp_out[STRINGMAXLEN_IMGNAME] = "";
    char imnamepha_out[STRINGMAXLEN_IMGNAME] = "";

    OptSystProp_init_base_cubes(optsyst, index, sharedmem);

    long elemstart1 = OptSystProp_find_start_elem(index, elemstart);
    long emax       = elemend;
    if (emax > optsyst[index].NBelem)
    {
        emax = optsyst[index].NBelem;
    }

    /* Used after loop to check whether final optical element was reached */
    long elem;
    for (elem = elemstart1; elem < emax; elem++)
    {
        OptSystProp_propagate_to_elem(optsyst,
                                      index,
                                      elem,
                                      sharedmem,
                                      imnameamp_out,
                                      imnamepha_out);

        printf("Applying element %ld\n", elem);
        fflush(stdout);

        OptSystProp_elem_ctx ctx = {
            .optsyst       = optsyst,
            .index         = index,
            .elem          = elem,
            .IDa           = image_ID(imnameamp_out, dcimg, dcnimg),
            .IDp           = image_ID(imnamepha_out, dcimg, dcnimg),
            .imnameamp_out = imnameamp_out,
            .imnamepha_out = imnamepha_out,
            .savedir       = savedir,
            .sharedmem     = sharedmem,
        };

        errno_t status = OptSystProp_apply_element(&ctx);
        if (status != RETURN_SUCCESS)
        {
            DEBUG_TRACE_FEXIT();
            return status;
        }

        OptSystProp_record_elem_plane(&ctx);
    }

    if ((elem == optsyst[index].NBelem) && (optsyst[index].endmode == 0))
    {
        OptSystProp_compute_final_psf(optsyst,
                                      index,
                                      imnameamp_out,
                                      imnamepha_out,
                                      savedir,
                                      sharedmem);
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
