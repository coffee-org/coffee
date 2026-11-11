/**
 * @file    makePIAAshapes.c
 * @brief   PIAA-type coronagraph design, make PIAA shapes
 *
 * Can design both APLCMC and PIAACMC coronagraphs.
 */

#include <stdint.h>
#include <stdio.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "OptSystProp/OptSystProp.h"
#include "OpticsMaterials/OpticsMaterials.h"
#include "linopt_imtools/linopt_imtools.h"

#include "PIAACMCsimul.h"
#include "PIAAshape/makePIAAshapes.h"

/**
 * add_centered_2d_image - Add a smaller square mode image into the center of a square map.
 * @dst_id:   Destination square image ID (`dst_size x dst_size`)
 * @dst_size: Linear dimension [pix] of destination image
 * @src_id:   Source square mode image ID (`src_size x src_size`)
 */
static void add_centered_2d_image(
    imageID  dst_id,
    uint32_t dst_size,
    imageID  src_id)
{
    uint32_t src_size = dcimg[src_id].md[0].size[0];
    uint32_t offset   = (dst_size - src_size) / 2;

    for (uint32_t ii = 0; ii < src_size; ii++)
    {
        for (uint32_t jj = 0; jj < src_size; jj++)
        {
            uint64_t dst_idx = (uint64_t)(jj + offset) * dst_size + (ii + offset);
            dcimg[dst_id].array.F[dst_idx] += dcimg[src_id].array.F[(uint64_t) jj * src_size + ii];
        }
    }
}

/**
 * assemble_piaa_modal_surface - Reconstruct radial cosine + 2D Fourier modes into a PIAA sag map.
 * @c_coeff:  Input radial cosine coefficient image name
 * @c_tmp:    Scratch radial cosine 2D surface image name
 * @f_coeff:  Input 2D Fourier coefficient image name
 * @f_tmp:    Scratch 2D Fourier surface image name
 * @out_name: Output combined 2D PIAA surface image name
 * @out_id:   Output milk image ID for `out_name`
 *
 * Return: RETURN_SUCCESS or image creation failure status.
 */
static errno_t assemble_piaa_modal_surface(
    const char *restrict c_coeff,
    const char *restrict c_tmp,
    const char *restrict f_coeff,
    const char *restrict f_tmp,
    const char *restrict out_name,
    imageID             *out_id)
{
    imageID ID0;
    linopt_imtools_image_construct("Cmodes", c_coeff, c_tmp, &ID0);

    imageID ID1;
    linopt_imtools_image_construct("Fmodes", f_coeff, f_tmp, &ID1);

    uint32_t dst_size = piaacmcopticaldesign.size;
    imageID  ID       = image_ID(out_name, dcimg, dcnimg);
    if (ID == -1)
    {
        FUNC_CHECK_RETURN(create_2Dimage_ID(out_name, dst_size, dst_size, &ID));
    }

    uint64_t npix = (uint64_t) dst_size * dst_size;
    for (uint64_t ii = 0; ii < npix; ii++)
    {
        dcimg[ID].array.F[ii] = 0.0;
    }

    add_centered_2d_image(ID, dst_size, ID0);
    add_centered_2d_image(ID, dst_size, ID1);

    *out_id = ID;
    return RETURN_SUCCESS;
}

/**
 * write_refractive_index_table - Export wavelength dispersion table `ri_array.txt` for lenses.
 * @sag2opd_coeff0_out: Output sag-to-OPD conversion factor `(n0 - 1) / 2` at central wavelength
 *
 * Return: RETURN_SUCCESS or failure if `ri_array.txt` cannot be opened.
 */
static errno_t write_refractive_index_table(double *sag2opd_coeff0_out)
{
    char fname[STRINGMAXLEN_FULLFILENAME];
    WRITE_FULLFILENAME(fname, "%s/ri_array.txt", piaacmcparams.piaacmcconfdir);

    FILE *fpri = fopen(fname, "w");
    if (fpri == NULL)
    {
        FUNC_RETURN_FAILURE("Cannot create file %s", fname);
    }

    double ri0            = OpticsMaterials_n(piaacmcopticaldesign.PIAAmaterial_code,
                                   piaacmcopticaldesign.lambda);
    double sag2opd_coeff0 = (ri0 - 1.0) / 2.0;

    for (long k = 0; k < piaacmcopticaldesign.nblambda; k++)
    {
        double ri            = OpticsMaterials_n(piaacmcopticaldesign.PIAAmaterial_code,
                                      piaacmcopticaldesign.lambdaarray[k]);
        double sag2opd_coeff = (ri - 1.0) / 2.0;
        fprintf(fpri,
                "%g %.16f %.16f %.16f %.16f\n",
                piaacmcopticaldesign.lambdaarray[k],
                ri,
                ri0,
                sag2opd_coeff,
                sag2opd_coeff / sag2opd_coeff0);
    }
    fclose(fpri);

    *sag2opd_coeff0_out = sag2opd_coeff0;
    return RETURN_SUCCESS;
}

/**
 * build_refractive_piaa_sag - Convert effective mirror sag map into physical refractive lens sag.
 * @IDmirror: Input mirror-equivalent 2D sag image ID (`piaam0z` or `piaam1z`)
 * @sag_name: Output refractive lens sag image name (`piaar0zsag` or `piaar1zsag`)
 *
 * Return: RETURN_SUCCESS or image/FITS error status.
 */
static errno_t build_refractive_piaa_sag(
    imageID              IDmirror,
    const char *restrict sag_name)
{
    uint32_t dst_size = piaacmcopticaldesign.size;
    imageID  IDsag    = image_ID(sag_name, dcimg, dcnimg);
    int      mksag    = 0;

    if (IDsag == -1)
    {
        mksag = 1;
    }
    else if ((dcimg[IDsag].md[0].size[0] != dst_size) || (dcimg[IDsag].md[0].size[1] != dst_size))
    {
        FUNC_CHECK_RETURN(delete_image_ID(sag_name, DELETE_IMAGE_ERRMODE_IGNORE));
        mksag = 1;
    }

    if (mksag == 1)
    {
        FUNC_CHECK_RETURN(create_2Dimage_ID(sag_name, dst_size, dst_size, &IDsag));
    }

    double sag2opd_coeff0 = 1.0;
    FUNC_CHECK_RETURN(write_refractive_index_table(&sag2opd_coeff0));

    uint64_t npix = (uint64_t) dst_size * dst_size;
    for (uint64_t ii = 0; ii < npix; ii++)
    {
        dcimg[IDsag].array.F[ii] = dcimg[IDmirror].array.F[ii] / sag2opd_coeff0;
    }

    if (piaacmcparams.PIAACMC_save == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/%s.fits", piaacmcparams.piaacmcconfdir, sag_name);
        FUNC_CHECK_RETURN(save_fl_fits(sag_name, fname));
        DEBUG_TRACEPOINT("Saved %s to %s", sag_name, fname);
    }

    return RETURN_SUCCESS;
}

/**
 * build_piaa0_surface - Reconstruct PIAA element 0 mirror and optional refractive sag maps.
 *
 * Return: RETURN_SUCCESS or image/FITS error status.
 */
static errno_t build_piaa0_surface(void)
{
    imageID IDpiaam0z;
    FUNC_CHECK_RETURN(assemble_piaa_modal_surface("piaa0Cmodescoeff",
                                                  "piaa0Cz",
                                                  "piaa0Fmodescoeff",
                                                  "piaa0Fz",
                                                  "piaam0z",
                                                  &IDpiaam0z));

    DEBUG_TRACEPOINT("piaacmcparams.PIAACMC_save %d", piaacmcparams.PIAACMC_save);
    if (piaacmcparams.PIAACMC_save == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];

        WRITE_FULLFILENAME(fname, "%s/mkPIAAshapes_Cmodes.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("Cmodes", fname));

        WRITE_FULLFILENAME(fname,
                           "%s/mkPIAAshapes_piaa0Cmodescoeff.fits",
                           piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("piaa0Cmodescoeff", fname));

        WRITE_FULLFILENAME(fname, "%s/mkPIAAshapes_piaa0Cz.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("piaa0Cz", fname));

        WRITE_FULLFILENAME(fname, "%s/mkPIAAshapes_piaa0Fz.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("piaa0Fz", fname));

        WRITE_FULLFILENAME(fname, "%s/mkPIAAshapes_piaam0z.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("piaam0z", fname));
    }

    FUNC_CHECK_RETURN(delete_image_ID("piaa0Cz", DELETE_IMAGE_ERRMODE_IGNORE));
    FUNC_CHECK_RETURN(delete_image_ID("piaa0Fz", DELETE_IMAGE_ERRMODE_IGNORE));

    DEBUG_TRACEPOINT("PIAA material code = %d", piaacmcopticaldesign.PIAAmaterial_code);
    if (piaacmcopticaldesign.PIAAmaterial_code != 0)
    {
        FUNC_CHECK_RETURN(build_refractive_piaa_sag(IDpiaam0z, "piaar0zsag"));
    }

    return RETURN_SUCCESS;
}

/**
 * build_piaa1_surface - Reconstruct PIAA element 1 mirror and optional refractive sag maps.
 *
 * Return: RETURN_SUCCESS or image/FITS error status.
 */
static errno_t build_piaa1_surface(void)
{
    imageID IDpiaam1z;
    FUNC_CHECK_RETURN(assemble_piaa_modal_surface("piaa1Cmodescoeff",
                                                  "piaa1Cz",
                                                  "piaa1Fmodescoeff",
                                                  "piaa1Fz",
                                                  "piaam1z",
                                                  &IDpiaam1z));

    if (piaacmcparams.PIAACMC_save == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];

        WRITE_FULLFILENAME(fname, "%s/piaa1Cz.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("piaa1Cz", fname));

        WRITE_FULLFILENAME(fname, "%s/piaa1Fz.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("piaa1Fz", fname));

        WRITE_FULLFILENAME(fname, "%s/piaam1z.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("piaam1z", fname));
    }

    FUNC_CHECK_RETURN(delete_image_ID("piaa1Cz", DELETE_IMAGE_ERRMODE_IGNORE));
    FUNC_CHECK_RETURN(delete_image_ID("piaa1Fz", DELETE_IMAGE_ERRMODE_IGNORE));

    if (piaacmcopticaldesign.PIAAmaterial_code != 0)
    {
        FUNC_CHECK_RETURN(build_refractive_piaa_sag(IDpiaam1z, "piaar1zsag"));
    }

    return RETURN_SUCCESS;
}

/**
 * makePIAAshapes - Construct 2D PIAA mirror (`piaam0z`, `piaam1z`) and lens sag maps from modes.
 *
 * Evaluates linear combinations of radial cosine (`Cmodes`) and 2D Fourier (`Fmodes`) basis cubes
 * whenever PIAA mode is active and surface maps are missing or forced to rebuild.
 *
 * Return: RETURN_SUCCESS on completion, or milk error code on failure.
 */
errno_t makePIAAshapes(void)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("PIAA material code = %d", piaacmcopticaldesign.PIAAmaterial_code);
    DEBUG_TRACEPOINT("PIAAmode %d", piaacmcopticaldesign.PIAAmode);

    if (piaacmcopticaldesign.PIAAmode != 1)
    {
        DEBUG_TRACE_FEXIT();
        return RETURN_SUCCESS;
    }

    piaacmcparams.MAKE_PIAA0shape = 0;
    if ((piaacmcparams.FORCE_MAKE_PIAA0shape == 1) || (image_ID("piaam0z", dcimg, dcnimg) == -1) ||
        ((piaacmcopticaldesign.PIAAmaterial_code != 0) &&
         (image_ID("piaar0zsag", dcimg, dcnimg) == -1)))
    {
        piaacmcparams.MAKE_PIAA0shape = 1;
    }

    piaacmcparams.MAKE_PIAA1shape = 0;
    if ((piaacmcparams.FORCE_MAKE_PIAA1shape == 1) || (image_ID("piaam1z", dcimg, dcnimg) == -1) ||
        ((piaacmcopticaldesign.PIAAmaterial_code != 0) &&
         (image_ID("piaar1zsag", dcimg, dcnimg) == -1)))
    {
        piaacmcparams.MAKE_PIAA1shape = 1;
    }

    DEBUG_TRACEPOINT("piaacmcparams.MAKE_PIAA0shape %d", piaacmcparams.MAKE_PIAA0shape);
    if (piaacmcparams.MAKE_PIAA0shape == 1)
    {
        FUNC_CHECK_RETURN(build_piaa0_surface());
    }

    if (piaacmcparams.MAKE_PIAA1shape == 1)
    {
        FUNC_CHECK_RETURN(build_piaa1_surface());
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
