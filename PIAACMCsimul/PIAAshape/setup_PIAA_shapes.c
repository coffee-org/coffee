/**
 * @file    setup_PIAA_shapes.c
 * @brief   Initial 2D prolate apodization, radial PIAA sag solver, and modal fit pipeline
 */

#define DEBUGLOG

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "coronagraphs/coronagraphs.h"
#include "image_basic/image_basic.h"
#include "image_gen/image_gen.h"
#include "linopt_imtools/linopt_imtools.h"

#include "PIAACMCsimul.h"
#include "PIAAshape/init_geomPIAA_rad.h"
#include "PIAAshape/load2DRadialApodization.h"
#include "PIAAshape/mkPIAAMshapes_from_RadSag.h"
#include "PIAAshape/setup_PIAA_shapes.h"

/**
 * try_load_piaa_mode_coeffs - Load existing PIAA modal coefficient vectors from piaaref/.
 *
 * Return: RETURN_SUCCESS or failure if `APLCmaskCtransm.txt` exists but is unreadable.
 */
static errno_t try_load_piaa_mode_coeffs(void)
{
    DEBUG_TRACE_FSTART();

    piaacmcopticaldesign.piaa0CmodesID = image_ID("piaa0Cmodescoeff", dcimg, dcnimg);
    piaacmcopticaldesign.piaa0FmodesID = image_ID("piaa0Fmodescoeff", dcimg, dcnimg);
    piaacmcopticaldesign.piaa1CmodesID = image_ID("piaa1Cmodescoeff", dcimg, dcnimg);
    piaacmcopticaldesign.piaa1FmodesID = image_ID("piaa1Fmodescoeff", dcimg, dcnimg);

    EXECUTE_SYSTEM_COMMAND("mkdir -p %s/piaaref/", piaacmcparams.piaacmcconfdir);

    if ((piaacmcopticaldesign.piaa0CmodesID != -1) &&
        (piaacmcopticaldesign.piaa0FmodesID != -1) &&
        (piaacmcopticaldesign.piaa1CmodesID != -1) &&
        (piaacmcopticaldesign.piaa1FmodesID != -1))
    {
        DEBUG_TRACE_FEXIT();
        return RETURN_SUCCESS;
    }

    char fname[STRINGMAXLEN_FULLFILENAME];
    DEBUG_TRACEPOINT("loading mode coeffs");

    WRITE_FULLFILENAME(fname, "%s/piaaref/piaa0Cmodes.fits", piaacmcparams.piaacmcconfdir);
    FUNC_CHECK_RETURN(load_fits(fname,
                                "piaa0Cmodescoeff",
                                LOADFITS_ERRMODE_WARNING,
                                &(piaacmcopticaldesign.piaa0CmodesID)));

    WRITE_FULLFILENAME(fname, "%s/piaaref/piaa0Fmodes.fits", piaacmcparams.piaacmcconfdir);
    FUNC_CHECK_RETURN(load_fits(fname,
                                "piaa0Fmodescoeff",
                                LOADFITS_ERRMODE_WARNING,
                                &(piaacmcopticaldesign.piaa0FmodesID)));

    WRITE_FULLFILENAME(fname, "%s/piaaref/piaa1Cmodes.fits", piaacmcparams.piaacmcconfdir);
    FUNC_CHECK_RETURN(load_fits(fname,
                                "piaa1Cmodescoeff",
                                LOADFITS_ERRMODE_WARNING,
                                &(piaacmcopticaldesign.piaa1CmodesID)));

    WRITE_FULLFILENAME(fname, "%s/piaaref/piaa1Fmodes.fits", piaacmcparams.piaacmcconfdir);
    FUNC_CHECK_RETURN(load_fits(fname,
                                "piaa1Fmodescoeff",
                                LOADFITS_ERRMODE_WARNING,
                                &(piaacmcopticaldesign.piaa1FmodesID)));

    WRITE_FULLFILENAME(fname, "%s/piaaref/APLCmaskCtransm.txt", piaacmcparams.piaacmcconfdir);
    FILE *fp = fopen(fname, "r");
    if (fp != NULL)
    {
        DEBUG_TRACEPOINT("loading %s", fname);
        float tmpf = 0.0;
        if (fscanf(fp, "%f", &tmpf) == 1)
        {
            piaacmcopticaldesign.fpmaskamptransm = tmpf;
        }
        else
        {
            fclose(fp);
            FUNC_RETURN_FAILURE("Cannot read value from file %s", fname);
        }
        fclose(fp);
    }
    else
    {
        DEBUG_TRACEPOINT("file %s not found", fname);
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * compute_2d_prolate_apodization - Solve multi-grid 2D prolate spheroidal pupil apodizer.
 * @piaacmctype: 0 for idealized monochromatic FPM, non-zero for physical FPM
 * @size:        2D grid size [pix]
 * @beamrad:     Beam radius [pix]
 * @loaded:      Non-zero if configuration was already loaded from disk
 * @saveconf:    Set to 1 when default idealized mask transmission is newly calibrated
 *
 * Return: RETURN_SUCCESS or prolate solver error status.
 */
static errno_t compute_2d_prolate_apodization(
    long     piaacmctype,
    uint32_t size,
    float    beamrad,
    int      loaded,
    int     *saveconf)
{
    DEBUG_TRACE_FSTART();

    char fname[STRINGMAXLEN_FULLFILENAME];

    DEBUG_TRACEPOINT("Creating 2D apodization for idealized circular monochromatic PIAACMC");

    create_variable_ID("DFTZFACTOR", 2);
    create_variable_ID("PNBITER", 15);
    FUNC_CHECK_RETURN(coronagraph_make_2Dprolateld(piaacmcopticaldesign.fpmaskradld,
                                                   beamrad * 0.5,
                                                   piaacmcopticaldesign.centObs1,
                                                   "apotmp1",
                                                   size / 2,
                                                   "NULLim"));

    basic_resizeim("apotmp1", "apostart", size, size);
    FUNC_CHECK_RETURN(delete_image_ID("apotmp1", DELETE_IMAGE_ERRMODE_WARNING));

    create_variable_ID("DFTZFACTOR", 4);
    create_variable_ID("PNBITER", 5);
    FUNC_CHECK_RETURN(coronagraph_make_2Dprolateld(piaacmcopticaldesign.fpmaskradld,
                                                   beamrad,
                                                   piaacmcopticaldesign.centObs1,
                                                   "apo",
                                                   size,
                                                   "pupmaskim"));

    chname_image_ID("apo", "apostart");
    create_variable_ID("DFTZFACTOR", 8);
    create_variable_ID("PNBITER", 5);
    FUNC_CHECK_RETURN(coronagraph_make_2Dprolateld(piaacmcopticaldesign.fpmaskradld,
                                                   beamrad,
                                                   piaacmcopticaldesign.centObs1,
                                                   "apo",
                                                   size,
                                                   "pupmaskim"));

    chname_image_ID("apo", "apostart");
    create_variable_ID("DFTZFACTOR", 16);
    create_variable_ID("PNBITER", 10);
    FUNC_CHECK_RETURN(coronagraph_make_2Dprolateld(piaacmcopticaldesign.fpmaskradld,
                                                   beamrad,
                                                   piaacmcopticaldesign.centObs1,
                                                   "apo",
                                                   size,
                                                   "pupmaskim"));

    chname_image_ID("apo", "apo2Drad");
    WRITE_FULLFILENAME(fname, "%s/apo2Drad.fits", piaacmcparams.piaacmcconfdir);
    FUNC_CHECK_RETURN(save_fits("apo2Drad", fname));

    if (piaacmcopticaldesign.PIAAmode == 1)
    {
        WRITE_FULLFILENAME(fname, "%s/piaaref/apo2Drad.fits", piaacmcparams.piaacmcconfdir);
        save_fits("apo2Drad", fname);
    }

    if ((piaacmctype == 0) && (loaded == 0))
    {
        DEBUG_TRACEPOINT("Idealized focal plane mask");
        piaacmcopticaldesign.fpmaskamptransm = -coffee_dcvar_f(variable_ID("APLCmaskCtransm"));
        printf("FOCAL PLANE MASK TRANSM = %f\n", piaacmcopticaldesign.fpmaskamptransm);
        printf("Saving default configuration\n");
        fflush(stdout);
        *saveconf = 1;

        WRITE_FULLFILENAME(fname, "%s/piaaref/APLCmaskCtransm.txt", piaacmcparams.piaacmcconfdir);
        FILE *fp = fopen(fname, "w");
        fprintf(fp, "%.20f\n", piaacmcopticaldesign.fpmaskamptransm);
        fclose(fp);
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * ensure_apo2drad_and_split - Load or compute `apo2Drad` and split into PIAA and CPA parts.
 * @piaacmctype: 0 for idealized monochromatic FPM, non-zero for physical FPM
 * @size:        2D grid size [pix]
 * @beamrad:     Beam radius [pix]
 * @loaded:      Non-zero if configuration was already loaded from disk
 * @saveconf:    Output flag set if configuration must be saved
 *
 * Return: RETURN_SUCCESS or error status.
 */
static errno_t ensure_apo2drad_and_split(
    long     piaacmctype,
    uint32_t size,
    float    beamrad,
    int      loaded,
    int     *saveconf)
{
    DEBUG_TRACE_FSTART();

    char    fname[STRINGMAXLEN_FULLFILENAME];
    imageID IDtmp1 = -1;

    WRITE_FULLFILENAME(fname, "%s/apo2Drad.fits", piaacmcparams.piaacmcconfdir);
    FUNC_CHECK_RETURN(load_fits(fname, "apo2Drad", LOADFITS_ERRMODE_WARNING, &IDtmp1));

    if (IDtmp1 == -1)
    {
        DEBUG_TRACEPOINT("creating apodization");
        EXECUTE_SYSTEM_COMMAND("cd \"%s\" && { test ! -f piaaref/apo2Drad.fits || "
                               "cp piaaref/apo2Drad.fits apo2Drad.fits; }",
                               piaacmcparams.piaacmcconfdir);

        imageID IDtmp2 = -1;
        FUNC_CHECK_RETURN(load_fits(fname, "apo2Drad", LOADFITS_ERRMODE_WARNING, &IDtmp2));
        if (IDtmp2 == -1)
        {
            FUNC_CHECK_RETURN(
                compute_2d_prolate_apodization(piaacmctype, size, beamrad, loaded, saveconf));
        }
    }

    imageID  IDapo  = image_ID("apo2Drad", dcimg, dcnimg);
    uint32_t xsize  = dcimg[IDapo].md[0].size[0];
    uint32_t ysize  = dcimg[IDapo].md[0].size[1];
    uint64_t xysize = (uint64_t) xsize * ysize;

    imageID IDapo_PIAA;
    imageID IDapo_CPA;
    FUNC_CHECK_RETURN(create_2Dimage_ID("apo2Drad_PIAA", xsize, ysize, &IDapo_PIAA));
    FUNC_CHECK_RETURN(create_2Dimage_ID("apo2Drad_CPA", xsize, ysize, &IDapo_CPA));

    if (piaacmcopticaldesign.PIAAmode == 0)
    {
        DEBUG_TRACEPOINT("everything goes to the conventional apodizer");
        for (uint64_t ii = 0; ii < xysize; ii++)
        {
            dcimg[IDapo_PIAA].array.F[ii] = 1.0;
            dcimg[IDapo_CPA].array.F[ii]  = dcimg[IDapo].array.F[ii];
        }
    }
    else
    {
        DEBUG_TRACEPOINT("sharing podization, coeff = %f", piaacmcopticaldesign.PIAAcoeff);
        for (uint64_t ii = 0; ii < xysize; ii++)
        {
            double coeff                  = piaacmcopticaldesign.PIAAcoeff;
            dcimg[IDapo_PIAA].array.F[ii] = pow(dcimg[IDapo].array.F[ii], coeff);
            dcimg[IDapo_CPA].array.F[ii]  = pow(dcimg[IDapo].array.F[ii], 1.0 - coeff);
        }
    }

    copy_image_ID("apo2Drad_CPA", "prePIAA0mask", 0);
    FUNC_CHECK_RETURN(save_fits("prePIAA0mask", "prePIAA0mask.fits"));

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * build_2d_piaa_sags_from_apo - Solve 1D radial PIAA ODEs and interpolate 2D sag maps.
 * @beamrad: Beam radius [pix]
 *
 * Return: RETURN_SUCCESS or solver error status.
 */
static errno_t build_2d_piaa_sags_from_apo(float beamrad)
{
    DEBUG_TRACE_FSTART();

    char fname[STRINGMAXLEN_FULLFILENAME];

    FUNC_CHECK_RETURN(load2DRadialApodization("apo2Drad_PIAA", beamrad, "outApofit"));
    FUNC_CHECK_RETURN(init_geomPIAA_rad("outApofit"));

    DEBUG_TRACEPOINT("Make 2D sag maps");
    WRITE_FULLFILENAME(fname, "%s/PIAA_Mshapes.txt", piaacmcparams.piaacmcconfdir);

    double beamradpix = piaacmcopticaldesign.beamrad / piaacmcopticaldesign.pixscale;
    FUNC_CHECK_RETURN(mkPIAAMshapes_from_RadSag(fname,
                                                piaacmcopticaldesign.PIAAsep,
                                                piaacmcopticaldesign.beamrad,
                                                piaacmcopticaldesign.r0lim,
                                                piaacmcopticaldesign.r1lim,
                                                piaacmcopticaldesign.size,
                                                beamradpix,
                                                piaacmcopticaldesign.NBradpts,
                                                "piaam0z",
                                                "piaam1z"));

    DEBUG_TRACEPOINT("piaacmcparams.PIAACMC_save %d", piaacmcparams.PIAACMC_save);
    if (piaacmcparams.PIAACMC_save == 1)
    {
        DEBUG_TRACEPOINT("saving files");
        WRITE_FULLFILENAME(fname, "%s/piaam0z.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("piaam0z", fname));
        WRITE_FULLFILENAME(fname, "%s/piaam1z.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("piaam1z", fname));
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * fit_piaa_cosine_modes - Crop 2D PIAA sags, fit radial `Cmodes`, and reconstruct `piaa*Cz`.
 * @size:    Full 2D wavefront grid size [pix]
 * @beamrad: Beam radius [pix]
 *
 * Return: RETURN_SUCCESS or modal fit error status.
 */
static errno_t fit_piaa_cosine_modes(
    uint32_t size,
    float    beamrad)
{
    DEBUG_TRACE_FSTART();

    char     fname[STRINGMAXLEN_FULLFILENAME];
    imageID  ID0   = image_ID("Cmodes", dcimg, dcnimg);
    uint32_t size0 = dcimg[ID0].md[0].size[0];
    uint32_t off   = (size - size0) / 2;

    imageID IDcrop0;
    imageID IDcrop1;
    FUNC_CHECK_RETURN(create_2Dimage_ID("piaa0zcrop", size0, size0, &IDcrop0));
    FUNC_CHECK_RETURN(create_2Dimage_ID("piaa1zcrop", size0, size0, &IDcrop1));

    imageID IDm0 = image_ID("piaam0z", dcimg, dcnimg);
    imageID IDm1 = image_ID("piaam1z", dcimg, dcnimg);
    for (uint32_t ii = 0; ii < size0; ii++)
    {
        for (uint32_t jj = 0; jj < size0; jj++)
        {
            uint32_t src_idx                    = (jj + off) * size + (ii + off);
            dcimg[IDcrop0].array.F[jj * size0 + ii] = dcimg[IDm0].array.F[src_idx];
            dcimg[IDcrop1].array.F[jj * size0 + ii] = dcimg[IDm1].array.F[src_idx];
        }
    }

    make_disk("maskd", size0, size0, 0.5 * size0, 0.5 * size0, beamrad);
    make_2Dgridpix("gridpix", size0, size0, 1, 1, 0, 0);
    arith_image_mult("maskd", "gridpix", "maskfit");

    printf("--------- FITTING COSINE MODES ---------\n");
    fflush(stdout);

    FUNC_CHECK_RETURN(linopt_imtools_image_fitModes("piaa0zcrop",
                                                    "Cmodes",
                                                    "maskfit",
                                                    1.0e-6,
                                                    "piaa0Cmodescoeff",
                                                    0,
                                                    NULL));
    FUNC_CHECK_RETURN(linopt_imtools_image_fitModes("piaa1zcrop",
                                                    "Cmodes",
                                                    "maskfit",
                                                    1.0e-6,
                                                    "piaa1Cmodescoeff",
                                                    0,
                                                    NULL));

    FUNC_CHECK_RETURN(
        linopt_imtools_image_construct("Cmodes", "piaa0Cmodescoeff", "piaa0Cz", NULL));
    FUNC_CHECK_RETURN(
        linopt_imtools_image_construct("Cmodes", "piaa1Cmodescoeff", "piaa1Cz", NULL));

    if (piaacmcparams.PIAACMC_save == 1)
    {
        WRITE_FULLFILENAME(fname, "%s/initpiaacmc_piaa0Cz.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("piaa0Cz", fname));
        WRITE_FULLFILENAME(fname,
                           "%s/initpiaacmc_piaa0Cmodescoeff.fits",
                           piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("piaa0Cmodescoeff", fname));
        WRITE_FULLFILENAME(fname, "%s/initpiaacmc_Cmodes.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("Cmodes", fname));
        WRITE_FULLFILENAME(fname, "%s/initpiaacmc_piaa1Cz.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("piaa1Cz", fname));
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * fit_piaa_fourier_residuals_and_save - Fit `Fmodes` to cosine residuals and persist `piaaref/`.
 * @size: Full 2D wavefront grid size [pix]
 *
 * Return: RETURN_SUCCESS or modal fit / FITS error status.
 */
static errno_t fit_piaa_fourier_residuals_and_save(uint32_t size)
{
    DEBUG_TRACE_FSTART();

    char     fname[STRINGMAXLEN_FULLFILENAME];
    imageID  ID0cz = image_ID("piaa0Cz", dcimg, dcnimg);
    imageID  ID1cz = image_ID("piaa1Cz", dcimg, dcnimg);
    uint32_t size0 = dcimg[ID0cz].md[0].size[0];
    uint32_t off   = (size - size0) / 2;
    imageID  IDm0  = image_ID("piaam0z", dcimg, dcnimg);
    imageID  IDm1  = image_ID("piaam1z", dcimg, dcnimg);
    imageID  IDr0;
    imageID  IDr1;

    FUNC_CHECK_RETURN(create_2Dimage_ID("piaa0Cres", size0, size0, &IDr0));
    FUNC_CHECK_RETURN(create_2Dimage_ID("piaa1Cres", size0, size0, &IDr1));

    for (uint32_t ii = 0; ii < size0; ii++)
    {
        for (uint32_t jj = 0; jj < size0; jj++)
        {
            uint32_t src_idx = (jj + off) * size + (ii + off);
            uint32_t dst_idx = jj * size0 + ii;
            dcimg[IDr0].array.F[dst_idx] =
                dcimg[IDm0].array.F[src_idx] - dcimg[ID0cz].array.F[dst_idx];
            dcimg[IDr1].array.F[dst_idx] =
                dcimg[IDm1].array.F[src_idx] - dcimg[ID1cz].array.F[dst_idx];
        }
    }

    if (piaacmcparams.PIAACMC_save == 1)
    {
        WRITE_FULLFILENAME(fname, "%s/initpiaacmc_piaa0Cres.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("piaa0Cres", fname));
        WRITE_FULLFILENAME(fname, "%s/initpiaacmc_piaa1Cres.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("piaa1Cres", fname));
    }

    printf("--------- FITTING FOURIER MODES ---------\n");
    fflush(stdout);

    FUNC_CHECK_RETURN(linopt_imtools_image_fitModes("piaa0Cres",
                                                    "Fmodes",
                                                    "maskfit",
                                                    0.01,
                                                    "piaa0Fmodescoeff",
                                                    0,
                                                    NULL));
    FUNC_CHECK_RETURN(linopt_imtools_image_fitModes("piaa1Cres",
                                                    "Fmodes",
                                                    "maskfit",
                                                    0.01,
                                                    "piaa1Fmodescoeff",
                                                    0,
                                                    NULL));

    FUNC_CHECK_RETURN(delete_image_ID("piaa0zcrop", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("piaa1zcrop", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("maskfit", DELETE_IMAGE_ERRMODE_WARNING));

    piaacmcopticaldesign.piaa0CmodesID = image_ID("piaa0Cmodescoeff", dcimg, dcnimg);
    piaacmcopticaldesign.piaa0FmodesID = image_ID("piaa0Fmodescoeff", dcimg, dcnimg);
    piaacmcopticaldesign.piaa1CmodesID = image_ID("piaa1Cmodescoeff", dcimg, dcnimg);
    piaacmcopticaldesign.piaa1FmodesID = image_ID("piaa1Fmodescoeff", dcimg, dcnimg);

    WRITE_FULLFILENAME(fname, "%s/piaaref/piaa0Cmodes.fits", piaacmcparams.piaacmcconfdir);
    FUNC_CHECK_RETURN(save_fits(dcimg[piaacmcopticaldesign.piaa0CmodesID].name, fname));
    WRITE_FULLFILENAME(fname, "%s/piaaref/piaa0Fmodes.fits", piaacmcparams.piaacmcconfdir);
    FUNC_CHECK_RETURN(save_fits(dcimg[piaacmcopticaldesign.piaa0FmodesID].name, fname));
    WRITE_FULLFILENAME(fname, "%s/piaaref/piaa1Cmodes.fits", piaacmcparams.piaacmcconfdir);
    FUNC_CHECK_RETURN(save_fits(dcimg[piaacmcopticaldesign.piaa1CmodesID].name, fname));
    WRITE_FULLFILENAME(fname, "%s/piaaref/piaa1Fmodes.fits", piaacmcparams.piaacmcconfdir);
    FUNC_CHECK_RETURN(save_fits(dcimg[piaacmcopticaldesign.piaa1FmodesID].name, fname));

    EXECUTE_SYSTEM_COMMAND("cp %s/piaaref/* %s/",
                           piaacmcparams.piaacmcconfdir,
                           piaacmcparams.piaacmcconfdir);

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * setupPIAAshapes - Load or synthesize PIAA optical shapes and modal decompositions.
 * @piaacmctype: 0 for idealized monochromatic FPM, non-zero for physical FPM
 * @size:        2D wavefront grid size [pix]
 * @beamrad:     Beam radius [pix]
 * @loaded:      Non-zero if configuration was loaded from disk
 * @saveconf:    Output flag set to 1 if configuration must be persisted
 *
 * Return: RETURN_SUCCESS or error code on failure.
 */
errno_t setupPIAAshapes(
    long     piaacmctype,
    uint32_t size,
    float    beamrad,
    int      loaded,
    int     *saveconf)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %ld %u %f", piaacmctype, size, beamrad);

    EXECUTE_SYSTEM_COMMAND("mkdir -p %s/piaaref/", piaacmcparams.piaacmcconfdir);

    if (piaacmcopticaldesign.PIAAmode == 1)
    {
        FUNC_CHECK_RETURN(try_load_piaa_mode_coeffs());
    }

    if ((piaacmcopticaldesign.piaa0CmodesID == -1) ||
        (piaacmcopticaldesign.piaa0FmodesID == -1) ||
        (piaacmcopticaldesign.piaa1CmodesID == -1) ||
        (piaacmcopticaldesign.piaa1FmodesID == -1))
    {
        DEBUG_TRACEPOINT("creating piaa modes files");
        FUNC_CHECK_RETURN(
            ensure_apo2drad_and_split(piaacmctype, size, beamrad, loaded, saveconf));
        FUNC_CHECK_RETURN(build_2d_piaa_sags_from_apo(beamrad));
        FUNC_CHECK_RETURN(fit_piaa_cosine_modes(size, beamrad));
        FUNC_CHECK_RETURN(fit_piaa_fourier_residuals_and_save(size));
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
