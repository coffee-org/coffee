/**
 * @file    init_piaacmcopticalsystem.c
 * @brief   Initialize OptSystProp OPTSYST optical element train for PIAACMC/APLCMC simulation
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"
#include "init_piaacmcopticalsystem.h"
#include "init_piaacmcopticalsystem_internal.h"

/**
 * optsyst_setup_spectral_grid - Configure wavelength channels, grid geometry, and lambdalist.txt.
 * @beamradpix_out: Output beam radius in pixels
 */
static void optsyst_setup_spectral_grid(double *beamradpix_out)
{
    piaacmcopticalsystem.nblambda = piaacmcopticaldesign.nblambda;

    printf("lambda = %g\n", piaacmcopticaldesign.lambda);
    printf("LAMBDASTART = %g\n", piaacmcparams.LAMBDASTART);
    printf("LAMBDAEND = %g\n", piaacmcparams.LAMBDAEND);

    for (int k = 0; k < piaacmcopticalsystem.nblambda; k++)
    {
        piaacmcopticalsystem.lambdaarray[k] =
            piaacmcparams.LAMBDASTART +
            (0.5 + k) * (piaacmcparams.LAMBDAEND - piaacmcparams.LAMBDASTART) /
                piaacmcopticalsystem.nblambda;
    }

    if (piaacmcparams.PIAACMC_save == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/lambdalist.txt", piaacmcparams.piaacmcconfdir);

        FILE *fp = fopen(fname, "w");
        if (fp == NULL)
        {
            PRINT_ERROR("Cannot create file %s", fname);
            abort();
        }
        for (int k = 0; k < piaacmcopticalsystem.nblambda; k++)
        {
            fprintf(fp, "%02d %20g\n", k, piaacmcopticalsystem.lambdaarray[k]);
        }
        fclose(fp);
    }

    piaacmcopticalsystem.beamrad    = piaacmcopticaldesign.beamrad;
    piaacmcopticalsystem.size       = piaacmcopticaldesign.size;
    piaacmcopticalsystem.pixscale   = piaacmcopticaldesign.pixscale;
    piaacmcopticalsystem.DFTgridpad = 0;

    variableID IDv = variable_ID("PIAACMC_dftgrid");
    if (IDv != -1)
    {
        piaacmcopticalsystem.DFTgridpad = (long) (coffee_dcvar_f(IDv) + 0.001);
    }

    piaacmcopticalsystem.NB_asphsurfm = 2 + piaacmcopticaldesign.nbDM;
    piaacmcopticalsystem.NB_asphsurfr = 0;
    piaacmcopticalsystem.NBelem       = 100;

    *beamradpix_out = piaacmcopticalsystem.beamrad / piaacmcopticalsystem.pixscale;
}

/**
 * fill_default_pupil_cube - Populate 3D pupa0 cube from telpup or centrally obscured annulus.
 * @IDa:      Destination 3D pupa0 image ID
 * @size2:    Pixels per spectral slice (size * size)
 * @nblambda: Number of spectral slices
 */
static void fill_default_pupil_cube(
    imageID  IDa,
    uint64_t size2,
    long     nblambda)
{
    imageID IDr = image_ID("rcoord", dcimg, dcnimg);
    imageID ID  = image_ID("telpup", dcimg, dcnimg);
    if ((ID == -1) && (file_exists("telpup.fits") == 1))
    {
        load_fits("telpup.fits", "telpup", 1, &ID);
    }

    if (ID == -1)
    {
        for (long k = 0; k < nblambda; k++)
        {
            for (uint64_t ii = 0; ii < size2; ii++)
            {
                float rval = dcimg[IDr].array.F[ii];
                dcimg[IDa].array.F[k * size2 + ii] =
                    ((rval > piaacmcopticaldesign.centObs0) && (rval < 1.0f)) ? 1.0f : 0.0f;
            }
        }
    }
    else
    {
        for (long k = 0; k < nblambda; k++)
        {
            for (uint64_t ii = 0; ii < size2; ii++)
            {
                dcimg[IDa].array.F[k * size2 + ii] = (dcimg[ID].array.F[ii] > 0.5f) ? 1.0f : 0.0f;
            }
        }
    }
}

/**
 * optsyst_add_input_pupil - Load or synthesize pupa0 cube and register entrance pupil element.
 * @elem:     Current element index cursor
 * @size:     2D wavefront array dimension [pix]
 * @nblambda: Number of spectral slices
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t optsyst_add_input_pupil(
    long *elem,
    long  size,
    long  nblambda)
{
    long idx = *elem;

    snprintf(piaacmcopticalsystem.name[idx], STRINGMAXLEN_OPTSYST_ELEMNAME, "input pupil");
    piaacmcopticalsystem.elemtype[idx] = 1;

    char fname_pupa0[STRINGMAXLEN_FULLFILENAME];
    WRITE_FULLFILENAME(fname_pupa0, "%s/pupa0_%ld.fits", piaacmcparams.piaacmcconfdir, size);
    if (file_exists(fname_pupa0) == 1)
    {
        load_fits(fname_pupa0, "pupa0", LOADFITS_ERRMODE_ERROR, NULL);
    }

    imageID IDa = image_ID("pupa0", dcimg, dcnimg);
    if (IDa == -1)
    {
        printf("CREATING INPUT PUPIL\n");
        FUNC_CHECK_RETURN(create_3Dimage_ID("pupa0", size, size, nblambda, &IDa));
        fill_default_pupil_cube(IDa, (uint64_t) size * size, nblambda);
        WRITE_FULLFILENAME(fname_pupa0, "%s/pupa0_%ld.fits", piaacmcparams.piaacmcconfdir, size);
        FUNC_CHECK_RETURN(save_fl_fits("pupa0", fname_pupa0));
    }

    piaacmcopticalsystem.elemarrayindex[idx] = IDa;
    piaacmcopticalsystem.elemZpos[idx]       = 0.0;
    (*elem)++;

    return RETURN_SUCCESS;
}

/**
 * optsyst_add_ttm_and_dms - Create entrance tip/tilt + OPD mirror and append active DMs.
 * @elem:       Current element index cursor
 * @TTxld:      Source x pointing offset [lambda/D]
 * @TTyld:      Source y pointing offset [lambda/D]
 * @beamradpix: Nominal beam radius [pix]
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t optsyst_add_ttm_and_dms(
    long  *elem,
    double TTxld,
    double TTyld,
    double beamradpix)
{
    long    size = piaacmcopticalsystem.size;
    imageID ID;
    FUNC_CHECK_RETURN(create_2Dimage_ID("TTm", size, size, &ID));

    double lam_mid = (piaacmcparams.LAMBDAEND + piaacmcparams.LAMBDASTART) * 0.5;
    for (long ii = 0; ii < size; ii++)
    {
        for (long jj = 0; jj < size; jj++)
        {
            double x = (1.0 * ii - 0.5 * size) / beamradpix;
            double y = (1.0 * jj - 0.5 * size) / beamradpix;

            dcimg[ID].array.F[jj * size + ii] = 0.25 * (TTxld * x + TTyld * y) * lam_mid;
        }
    }

    imageID IDopderr = image_ID("opderr", dcimg, dcnimg);
    if (IDopderr != -1)
    {
        for (long ii = 0; ii < size; ii++)
        {
            for (long jj = 0; jj < size; jj++)
            {
                dcimg[ID].array.F[jj * size + ii] += dcimg[IDopderr].array.F[jj * size + ii] * 0.5;
            }
        }
    }

    long idx = *elem;
    snprintf(piaacmcopticalsystem.name[idx], STRINGMAXLEN_OPTSYST_ELEMNAME, "TT mirror");
    piaacmcopticalsystem.elemtype[idx]            = 3;
    piaacmcopticalsystem.elemarrayindex[idx]      = 0;
    piaacmcopticalsystem.ASPHSURFMarray[0].surfID = ID;
    piaacmcopticalsystem.elemZpos[idx]            = 0.0;
    (*elem)++;

    for (int iDM = 0; iDM < piaacmcopticaldesign.nbDM; iDM++)
    {
        idx = *elem;
        snprintf(piaacmcopticalsystem.name[idx], STRINGMAXLEN_OPTSYST_ELEMNAME, "DM %d", iDM);
        piaacmcopticalsystem.elemtype[idx]       = 3;
        piaacmcopticalsystem.elemarrayindex[idx] = 3 + iDM;
        piaacmcopticalsystem.ASPHSURFMarray[3 + iDM].surfID = piaacmcopticaldesign.ID_DM[iDM];
        piaacmcopticalsystem.elemZpos[idx]                  = piaacmcopticaldesign.DMpos[iDM];
        (*elem)++;
    }

    return RETURN_SUCCESS;
}

/**
 * optsyst_save_conjugations - Write optical element conjugations to conjugations.txt when enabled.
 */
static void optsyst_save_conjugations(void)
{
    if (piaacmcparams.PIAACMC_save != 1)
    {
        return;
    }

    char fname[STRINGMAXLEN_FULLFILENAME];
    WRITE_FULLFILENAME(fname, "%s/conjugations.txt", piaacmcparams.piaacmcconfdir);

    FILE *fp = fopen(fname, "w");
    if (fp == NULL)
    {
        PRINT_ERROR("Cannot create %s", fname);
        abort();
    }

    for (long elem = 0; elem < piaacmcopticalsystem.NBelem; elem++)
    {
        fprintf(
            fp,
            "%02ld  %f    %s\n",
            elem,
            piaacmcopticalsystem.elemZpos[elem],
            piaacmcopticalsystem.name[elem]);
    }
    fclose(fp);
}

/**
 * init_piaacmcopticalsystem - Assemble OPTSYST plane sequence for source offset (TTxld, TTyld).
 * @TTxld: Off-axis source tip coordinate [lambda/D]
 * @TTyld: Off-axis source tilt coordinate [lambda/D]
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
errno_t init_piaacmcopticalsystem(
    double TTxld,
    double TTyld)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %f %f", TTxld, TTyld);

    double beamradpix = 0.0;
    optsyst_setup_spectral_grid(&beamradpix);

    long size     = piaacmcopticalsystem.size;
    long nblambda = piaacmcopticalsystem.nblambda;
    long elem     = 0;

    FUNC_CHECK_RETURN(optsyst_add_input_pupil(&elem, size, nblambda));
    FUNC_CHECK_RETURN(optsyst_add_ttm_and_dms(&elem, TTxld, TTyld, beamradpix));
    FUNC_CHECK_RETURN(piaacmc_optsyst_add_piaa_stage(&elem, size, beamradpix));
    FUNC_CHECK_RETURN(piaacmc_optsyst_add_fpm(&elem));
    piaacmc_optsyst_add_lyot_and_backend(&elem, size);

    piaacmcopticalsystem.NBelem  = elem;
    piaacmcopticalsystem.endmode = 0;

    optsyst_save_conjugations();
    piaacmcparams.optsystinit = 1;

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
