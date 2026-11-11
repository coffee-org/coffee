/**
 * @file    PIAAshape_modes.c
 * @brief   Pupil coordinate grids and radial Cosine / 2D Fourier modal bases for PIAA optics
 */

#define DEBUGLOG

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "linopt_imtools/linopt_imtools.h"

#include "PIAACMCsimul.h"
#include "PIAAshape/PIAAshape_modes.h"

/**
 * make_x_y_r_PA_images - Pre-compute normalized pupil x, y, r, and position-angle maps.
 * @size:    2D square grid linear size [pix]
 * @beamrad: Beam radius in pixel units
 *
 * Return: RETURN_SUCCESS or image allocation error status.
 */
errno_t make_x_y_r_PA_images(
    long  size,
    float beamrad)
{
    DEBUG_TRACE_FSTART();

    imageID IDx;
    imageID IDy;
    imageID IDr;
    imageID IDPA;

    FUNC_CHECK_RETURN(create_2Dimage_ID("xcoord", size, size, &IDx));
    FUNC_CHECK_RETURN(create_2Dimage_ID("ycoord", size, size, &IDy));
    FUNC_CHECK_RETURN(create_2Dimage_ID("rcoord", size, size, &IDr));
    FUNC_CHECK_RETURN(create_2Dimage_ID("PAcoord", size, size, &IDPA));

    printf("pre-computing x, y, r, and PA\n");
    fflush(stdout);

    for (uint32_t ii = 0; ii < (uint32_t) size; ii++)
    {
        for (uint32_t jj = 0; jj < (uint32_t) size; jj++)
        {
            float x                             = (1.0 * ii - 0.5 * size) / beamrad;
            float y                             = (1.0 * jj - 0.5 * size) / beamrad;
            dcimg[IDx].array.F[jj * size + ii]  = x;
            dcimg[IDy].array.F[jj * size + ii]  = y;
            dcimg[IDr].array.F[jj * size + ii]  = sqrt(x * x + y * y);
            dcimg[IDPA].array.F[jj * size + ii] = atan2(y, x);
        }
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * load_or_create_cmodes - Load cached radial cosine basis `Cmodes` or build and save it.
 * @size:    2D wavefront grid size [pix]
 * @beamrad: Beam radius in pixel units
 */
static void load_or_create_cmodes(
    long  size,
    float beamrad)
{
    char fname[STRINGMAXLEN_FULLFILENAME];

    piaacmcparams.CREATE_Cmodes = 0;
    WRITE_FULLFILENAME(fname, "ref/Cmodes_%ld.fits", piaacmcopticaldesign.size);

    DEBUG_TRACEPOINT("piaacmcparams.FORCE_CREATE_Cmodes = %d",
                     piaacmcparams.FORCE_CREATE_Cmodes);

    if (piaacmcparams.FORCE_CREATE_Cmodes == 0)
    {
        piaacmcopticaldesign.CmodesID = image_ID("Cmodes", dcimg, dcnimg);
        if (piaacmcopticaldesign.CmodesID == -1)
        {
            load_fits(fname, "Cmodes", 0, &(piaacmcopticaldesign.CmodesID));
        }
        if (piaacmcopticaldesign.CmodesID == -1)
        {
            piaacmcparams.CREATE_Cmodes = 1;
        }
    }
    else
    {
        piaacmcparams.CREATE_Cmodes = 1;
    }

    if (piaacmcparams.CREATE_Cmodes == 1)
    {
        if (piaacmcopticaldesign.CmodesID != -1)
        {
            delete_image_ID("Cmodes", DELETE_IMAGE_ERRMODE_WARNING);
        }
        long Cmsize = (long) (beamrad * 4);
        if (Cmsize > size)
        {
            Cmsize = size;
        }
        printf("beamradpix = %f -> Cmsize = %ld\n", beamrad, Cmsize);
        if (Cmsize % 2 == 1)
        {
            Cmsize++;
        }
        piaacmcopticaldesign.Cmsize = Cmsize;
        DEBUG_TRACEPOINT("run linopt_imtools_makeCosRadModes");
        linopt_imtools_makeCosRadModes("Cmodes",
                                       Cmsize,
                                       piaacmcopticaldesign.piaaNBCmodesmax,
                                       ApoFitCosFact * beamrad,
                                       2.0,
                                       NULL);
        piaacmcopticaldesign.CmodesID = image_ID("Cmodes", dcimg, dcnimg);
        save_fits("Cmodes", fname);

        EXECUTE_SYSTEM_COMMAND("mv ModesExpr_CosRad.txt %s/", piaacmcparams.piaacmcconfdir);
    }

    piaacmcopticaldesign.NBCmodes = dcimg[piaacmcopticaldesign.CmodesID].md[0].size[2];
    piaacmcopticaldesign.Cmsize   = dcimg[piaacmcopticaldesign.CmodesID].md[0].size[0];
}

/**
 * load_or_create_fmodes - Load cached 2D Fourier basis `Fmodes` or build and save it.
 * @size:    2D wavefront grid size [pix]
 * @beamrad: Beam radius in pixel units
 *
 * Return: RETURN_SUCCESS or FITS/mode-generation error status.
 */
static errno_t load_or_create_fmodes(
    long  size,
    float beamrad)
{
    DEBUG_TRACE_FSTART();

    char fname[STRINGMAXLEN_FULLFILENAME];

    piaacmcparams.CREATE_Fmodes = 0;
    WRITE_FULLFILENAME(fname, "ref/Fmodes_%ld.fits", piaacmcopticaldesign.size);

    DEBUG_TRACEPOINT("piaacmcparams.FORCE_CREATE_Fmodes = %d",
                     piaacmcparams.FORCE_CREATE_Fmodes);

    if (piaacmcparams.FORCE_CREATE_Fmodes == 0)
    {
        piaacmcopticaldesign.FmodesID = image_ID("Fmodes", dcimg, dcnimg);
        if (piaacmcopticaldesign.FmodesID == -1)
        {
            load_fits(fname, "Fmodes", LOADFITS_ERRMODE_IGNORE, &(piaacmcopticaldesign.FmodesID));
        }
        if (piaacmcopticaldesign.FmodesID == -1)
        {
            piaacmcparams.CREATE_Fmodes = 1;
        }
    }
    else
    {
        piaacmcparams.CREATE_Fmodes = 1;
    }

    if (piaacmcparams.CREATE_Fmodes == 1)
    {
        long Fmsize = (long) (beamrad * 4);
        if (Fmsize > size)
        {
            Fmsize = size;
        }
        piaacmcopticaldesign.Fmsize = Fmsize;
        DEBUG_TRACEPOINT("run linopt_imtools_makeCPAmodes");
        IMGID imgFmodes = imgid_make_from_name("Fmodes");
        IMGID imgNull   = imgid_make_from_name("NULL");
        FUNC_CHECK_RETURN(linopt_imtools_makeCPAmodes(&imgFmodes,
                                                      Fmsize,
                                                      Fmsize,
                                                      0.5 * Fmsize,
                                                      0.5 * Fmsize,
                                                      0.0,
                                                      piaacmcopticaldesign.piaaCPAmax,
                                                      piaacmcopticaldesign.piaaCPAmax,
                                                      0.8,
                                                      beamrad,
                                                      2.0,
                                                      0.0,
                                                      0.0,
                                                      0.0,
                                                      1,
                                                      NULL,
                                                      imgNull,
                                                      0.0,
                                                      0.0));
        imgid_free(&imgNull);
        piaacmcopticaldesign.FmodesID = image_ID("Fmodes", dcimg, dcnimg);

        FUNC_CHECK_RETURN(save_fits("Fmodes", fname));
        FUNC_CHECK_RETURN(save_fits("cpamodesfreq", "cpamodesfreq.fits"));
        EXECUTE_SYSTEM_COMMAND("mv ModesExpr_CPA.txt %s/", piaacmcparams.piaacmcconfdir);
    }

    piaacmcopticaldesign.NBFmodes = dcimg[piaacmcopticaldesign.FmodesID].md[0].size[2];
    piaacmcopticaldesign.Fmsize   = dcimg[piaacmcopticaldesign.FmodesID].md[0].size[0];

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * make_C_F_modes - Initialize radial cosine (`Cmodes`) and 2D Fourier (`Fmodes`) bases.
 * @size:    2D wavefront grid size [pix]
 * @beamrad: Beam radius in pixel units
 *
 * Return: RETURN_SUCCESS or mode-generation error status.
 */
errno_t make_C_F_modes(
    long  size,
    float beamrad)
{
    DEBUG_TRACE_FSTART();

    printf("Creating / loading Cmodes and Fmodes ...\n");
    fflush(stdout);

    load_or_create_cmodes(size, beamrad);
    FUNC_CHECK_RETURN(load_or_create_fmodes(size, beamrad));

    printf("DONE Creating / loading Cmodes and Fmodes\n");
    fflush(stdout);

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
