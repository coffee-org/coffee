/**
 * @file    PIAACMCsimul_regularization.c
 * @brief   PIAA modal surface and FPM zone sag regularization penalties and 1D vectors
 */

#include <math.h>
#include <stdint.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_memory/COREMOD_memory.h"

#include "PIAACMCsimul.h"
#include "PIAACMCsimul_regularization.h"

/**
 * ensure_piaa_ref_vector - Look up or zero-initialize a 1D PIAA modal reference vector.
 * @refname:   Reference image name (`piaa0Cmref`, `piaa1Cmref`, `piaa0Fmref`, `piaa1Fmref`)
 * @alloc_len: Vector length to allocate if not yet present in memory
 *
 * Return: Image ID of the reference modal vector.
 */
static imageID ensure_piaa_ref_vector(
    const char *refname,
    uint32_t    alloc_len)
{
    imageID IDref = image_ID(refname, dcimg, dcnimg);

    if (IDref == -1)
    {
        create_2Dimage_ID(refname, alloc_len, 1, &IDref);
        for (uint32_t jj = 0; jj < alloc_len; jj++)
        {
            dcimg[IDref].array.F[jj] = 0.0;
        }
    }

    return IDref;
}

/**
 * PIAACMCsimul_regularization_PIAAshapes_value - Compute PIAA shape regularization penalty.
 *
 * Sums squared frequency-weighted modal deviations from `piaa[01][CF]mref` across
 * radial cosine (`Cmodes`) and 2D Fourier (`Fmodes`) bases for PIAA optics 0 and 1.
 *
 * Return: Scalar sum-of-squares regularization penalty.
 */
double PIAACMCsimul_regularization_PIAAshapes_value(void)
{
    DEBUG_TRACE_FSTART();

    double   value      = 0.0;
    uint32_t nc0        = dcimg[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0];
    uint32_t nc1        = dcimg[piaacmcopticaldesign.piaa1CmodesID].md[0].size[0];
    uint32_t nf0        = dcimg[piaacmcopticaldesign.piaa0FmodesID].md[0].size[0];
    uint32_t nf1        = dcimg[piaacmcopticaldesign.piaa1FmodesID].md[0].size[0];
    imageID  ID_CPAfreq = image_ID("cpamodesfreq", dcimg, dcnimg);

    imageID ID    = piaacmcopticaldesign.piaa0CmodesID;
    imageID IDref = ensure_piaa_ref_vector("piaa0Cmref", nc0);
    for (uint32_t jj = 0; jj < nc0; jj++)
    {
        double tmp = piaacmcparams.linopt_piaa0C_regcoeff *
                     (dcimg[ID].array.F[jj] - dcimg[IDref].array.F[jj]) *
                     pow(1.0 * jj, piaacmcparams.linopt_piaa0C_regcoeff_alpha);
        value += tmp * tmp;
    }

    ID    = piaacmcopticaldesign.piaa1CmodesID;
    IDref = ensure_piaa_ref_vector("piaa1Cmref", nc0);
    for (uint32_t jj = 0; jj < nc1; jj++)
    {
        double tmp = piaacmcparams.linopt_piaa1C_regcoeff *
                     (dcimg[ID].array.F[jj] - dcimg[IDref].array.F[jj]) *
                     pow(1.0 * jj, piaacmcparams.linopt_piaa1C_regcoeff_alpha);
        value += tmp * tmp;
    }

    ID    = piaacmcopticaldesign.piaa0FmodesID;
    IDref = ensure_piaa_ref_vector("piaa0Fmref", nf0);
    for (uint32_t jj = 0; jj < nf0; jj++)
    {
        double tmp = piaacmcparams.linopt_piaa0F_regcoeff *
                     (dcimg[ID].array.F[jj] - dcimg[IDref].array.F[jj]) *
                     pow(1.0 * dcimg[ID_CPAfreq].array.F[jj],
                         piaacmcparams.linopt_piaa0F_regcoeff_alpha);
        value += tmp * tmp;
    }

    ID    = piaacmcopticaldesign.piaa1FmodesID;
    IDref = ensure_piaa_ref_vector("piaa1Fmref", nf1);
    for (uint32_t jj = 0; jj < nf1; jj++)
    {
        double tmp = piaacmcparams.linopt_piaa1F_regcoeff *
                     (dcimg[ID].array.F[jj] - dcimg[IDref].array.F[jj]) *
                     pow(1.0 * dcimg[ID_CPAfreq].array.F[jj],
                         piaacmcparams.linopt_piaa1F_regcoeff_alpha);
        value += tmp * tmp;
    }

    DEBUG_TRACE_FEXIT();
    return value;
}

/**
 * PIAACMCsimul_regularization_fpmsag_value - Compute FPM zone sag regularization penalty.
 *
 * Return: Sum over zones of `((sag / fpmsagreg_coeff) ^ fpmsagreg_alpha) ^ 2`.
 */
double PIAACMCsimul_regularization_fpmsag_value(void)
{
    DEBUG_TRACE_FSTART();

    double  regvalue = 0.0;
    imageID IDzonez  = piaacmcopticaldesign.zonezID;

    for (long zoneindex = 0; zoneindex < dcimg[IDzonez].md[0].size[0]; zoneindex++)
    {
        double tmp = pow(dcimg[IDzonez].array.D[zoneindex] / piaacmcopticaldesign.fpmsagreg_coeff,
                         piaacmcopticaldesign.fpmsagreg_alpha);
        regvalue += tmp * tmp;
    }

    DEBUG_TRACE_FEXIT();
    return regvalue;
}

/**
 * PIAACMCsimul_regularization_PIAAshapes_add1Dvector - Append PIAA regularization residuals.
 * @ID1D:   Destination 1D float evaluation vector image ID
 * @index0: Starting element index in `ID1D`
 *
 * Return: Next write index in `ID1D` after appending all 4 PIAA modal penalty blocks.
 */
long PIAACMCsimul_regularization_PIAAshapes_add1Dvector(
    imageID ID1D,
    long    index0)
{
    DEBUG_TRACE_FSTART();

    long     vindex     = index0;
    uint32_t nc0        = dcimg[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0];
    uint32_t nc1        = dcimg[piaacmcopticaldesign.piaa1CmodesID].md[0].size[0];
    uint32_t nf0        = dcimg[piaacmcopticaldesign.piaa0FmodesID].md[0].size[0];
    uint32_t nf1        = dcimg[piaacmcopticaldesign.piaa1FmodesID].md[0].size[0];
    imageID  ID_CPAfreq = image_ID("cpamodesfreq", dcimg, dcnimg);

    imageID ID    = piaacmcopticaldesign.piaa0CmodesID;
    imageID IDref = ensure_piaa_ref_vector("piaa0Cmref", nc0);
    for (uint32_t ii = 0; ii < nc0; ii++)
    {
        dcimg[ID1D].array.F[vindex] = piaacmcparams.linopt_piaa0C_regcoeff *
                                      (dcimg[ID].array.F[ii] - dcimg[IDref].array.F[ii]) *
                                      pow(1.0 * ii, piaacmcparams.linopt_piaa0C_regcoeff_alpha);
        vindex++;
    }

    ID    = piaacmcopticaldesign.piaa1CmodesID;
    IDref = ensure_piaa_ref_vector("piaa1Cmref", nc0);
    for (uint32_t ii = 0; ii < nc1; ii++)
    {
        dcimg[ID1D].array.F[vindex] = piaacmcparams.linopt_piaa1C_regcoeff *
                                      (dcimg[ID].array.F[ii] - dcimg[IDref].array.F[ii]) *
                                      pow(1.0 * ii, piaacmcparams.linopt_piaa1C_regcoeff_alpha);
        vindex++;
    }

    ID    = piaacmcopticaldesign.piaa0FmodesID;
    IDref = ensure_piaa_ref_vector("piaa0Fmref", nf0);
    for (uint32_t ii = 0; ii < nf0; ii++)
    {
        dcimg[ID1D].array.F[vindex] = piaacmcparams.linopt_piaa0F_regcoeff *
                                      (dcimg[ID].array.F[ii] - dcimg[IDref].array.F[ii]) *
                                      pow(1.0 * dcimg[ID_CPAfreq].array.F[ii],
                                          piaacmcparams.linopt_piaa0F_regcoeff_alpha);
        vindex++;
    }

    ID    = piaacmcopticaldesign.piaa1FmodesID;
    IDref = ensure_piaa_ref_vector("piaa1Fmref", nf1);
    for (uint32_t ii = 0; ii < nf1; ii++)
    {
        dcimg[ID1D].array.F[vindex] = piaacmcparams.linopt_piaa1F_regcoeff *
                                      (dcimg[ID].array.F[ii] - dcimg[IDref].array.F[ii]) *
                                      pow(1.0 * dcimg[ID_CPAfreq].array.F[ii],
                                          piaacmcparams.linopt_piaa1F_regcoeff_alpha);
        vindex++;
    }

    DEBUG_TRACE_FEXIT();
    return vindex;
}

/**
 * PIAACMCsimul_regularization_fpmsag_add1Dvector - Append FPM sag regularization residuals.
 * @ID1D:   Destination 1D float evaluation vector image ID
 * @IDmask: Optional mask image ID to fill with `1.0` (`-1` to skip mask write)
 * @index0: Starting element index in `ID1D`
 *
 * Return: Next write index in `ID1D` after appending all zone sag terms.
 */
long PIAACMCsimul_regularization_fpmsag_add1Dvector(
    imageID ID1D,
    imageID IDmask,
    long    index0)
{
    long    vindex = index0;
    imageID ID     = piaacmcopticaldesign.zonezID;

    for (uint32_t jj = 0; jj < dcimg[ID].md[0].size[0]; jj++)
    {
        dcimg[ID1D].array.F[vindex] =
            pow(dcimg[ID].array.D[jj] / piaacmcopticaldesign.fpmsagreg_coeff,
                piaacmcopticaldesign.fpmsagreg_alpha);
        if (IDmask != -1)
        {
            dcimg[IDmask].array.F[vindex] = 1.0;
        }
        vindex++;
    }

    return vindex;
}
