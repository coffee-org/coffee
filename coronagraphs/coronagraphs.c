/**
 * @file    coronagraphs.c
 * @brief   Shared module state and common pupil/FFT/PSF helpers for coffeecoronagraphs
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "coffee_compat.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "fft/fft.h"
#include "image_gen/image_gen.h"

#include "coronagraphs_internal.h"

double CORONAGRAPHS_PIXSCALE =
    CORONAGRAPHS_TDIAM / CORONAGRAPHS_PSCALE / CORONAGRAPHS_ARRAYSIZE;
double COROTMP1;
double MASKSIZELD;

coro_state_t g_coro = {
    .useDFT                 = 1,
    .DFTZFACTOR             = 8.0,
    .FPMASKSIZE_ERROR       = 0,
    .FPMASK_FACTOR          = 1.0,
    .FPMASK_FACTOR1         = 1.0,
    .FPMASK_FACTOR2         = 1.1,
    .FPM_TRANSM1            = 0.0,
    .FPM_TRANSM2            = 1.0,
    .BL8MASK_m              = 1.0,
    .BL8MASK_l              = 3.0,
    .BL8MASK_eps            = 0.6,
    .BL8MODE                = 1,
    .BL4MASK_eps            = 0.21,
    .ODC_GAUSS              = 10.0,
    .ODCMASK_eps            = 0.85,
    .SHEAR4_SHEAR           = 0.1,
    .RRPM_RADIUS            = 0.53773,
    .RRPM_P2                = -0.3207,
    .RRPM_P3                = 0.00952,
    .RRPM_P4                = 0.000658,
    .RRPM_P5                = 0.001127,
    .RRPM_PIAA              = 1,
    .PIAACENTOBS            = 0.0,
    .PIAAAPO_NBPOINTS       = 2042,
    .PIAAAPO_FNAME          = "pup_10.411256_2.0",
    .PIAAAPODIZE_2DAPOFNAME = "APLCapo_4.200.1024.ref.gz",
    .PIAAFPMASK             = 0,
    .PIAALOWFS              = 1,
    .PIAAFPMASKRAD          = 5.5,
    .piaaconfdirection      = 0,
    .PIAAOVERSIZE           = 1.025,
    .initPIAA               = 0,
    .AUTOPIAACMASK          = 1,
    .APLC_CentOBS0          = 0.0,
    .APLC_CentOBS1          = 0.0,
    .PIAAextfactor0         = 1.1,
    .PIAAextfactor1         = 1.1,
    .APLC_PIAA              = 0,
    .NB_APLC_STEP           = 0,
    .APLC_FPMASKsize        = 4.2,
    .APLC_FLIP              = 0,
    .APLC_PMASK             = 0,
    .FPMASK_transm_error    = 0.0,
    .FPMASK_size_error      = 0.0,
    .STRIPCOFFSET           = 0.4,
    .OVC_CHARGE             = 2,
    .fitapoINIT             = 0,
    .fitapoN                = 5,
    .LOOPCNT                = 0,
    .fitapo_minc            = 0.1,
    .APLCapo_CO_START       = 0.0,
    .APLCapo_CO_END         = 0.5,
    .APLCapo_CO_STEP        = 0.001,
    .APLCapo_FPMRAD_START   = 0.8,
    .APLCapo_FPMRAD_END     = 5.0,
    .APLCapo_FPMRAD_STEP    = 0.001,
};

/**
 * coro_pupil_radius_pix - Return nominal telescope pupil radius in pixels
 *
 * Return: Pupil radius [pix].
 */
double coro_pupil_radius_pix(void)
{
    return CORONAGRAPHS_TDIAM / CORONAGRAPHS_PSCALE / 2.0;
}

/**
 * coro_check_fscanf - Validate fscanf item count and report I/O errors
 * @fp:        Open file stream used in fscanf
 * @fscanfcnt: Return value from fscanf
 * @expected:  Expected number of matched items
 *
 * Return: RETURN_SUCCESS when fscanfcnt == expected, RETURN_FAILURE otherwise.
 */
errno_t coro_check_fscanf(
    FILE *fp,
    int   fscanfcnt,
    int   expected)
{
    if (fscanfcnt == EOF)
    {
        if (ferror(fp))
        {
            perror("fscanf");
        }
        else
        {
            fprintf(stderr, "Error: fscanf reached end of file\n");
        }
        return RETURN_FAILURE;
    }

    if (fscanfcnt != expected)
    {
        fprintf(
            stderr,
            "Error: fscanf matched %i input items, %i expected\n",
            fscanfcnt,
            expected);
        return RETURN_FAILURE;
    }

    return RETURN_SUCCESS;
}

/**
 * coro_fft2d_centered - Execute centered 2D FFT (permut -> do2dfft -> permut)
 * @in_name:  Input complex image name (permuted in-place)
 * @out_name: Output complex image name
 *
 * Return: RETURN_SUCCESS on success, RETURN_FAILURE on error.
 */
errno_t coro_fft2d_centered(
    const char *in_name,
    const char *out_name)
{
    permut(in_name);
    if (do2dfft(in_name, out_name) == -1)
    {
        return RETURN_FAILURE;
    }
    if (image_ID(out_name, dcimg, dcnimg) == -1)
    {
        return RETURN_FAILURE;
    }
    permut(out_name);
    return RETURN_SUCCESS;
}

/**
 * coro_setup_pupil_tiptilt - Build entrance pupil "pa1" and tip-tilt phase "pp1"
 * @xld:       Source X offset [lambda/D]
 * @yld:       Source Y offset [lambda/D]
 * @subpix:    1 for subpixel disk edge, 0 for binary disk
 * @oversize:  Multiplicative factor on pupil radius
 * @pha_name:  Optional phase aberration image name (e.g. "corphase")
 * @out_total: Output unapodized pupil flux integral
 */
void coro_setup_pupil_tiptilt(
    double      xld,
    double      yld,
    int         subpix,
    double      oversize,
    const char *pha_name,
    double     *out_total)
{
    uint32_t size     = CORONAGRAPHS_ARRAYSIZE;
    uint64_t size2    = (uint64_t) size * size;
    double   trad_pix = coro_pupil_radius_pix();
    imageID  IDp1;

    if (subpix != 0)
    {
        make_subpixdisk("pa1", size, size, 0.5 * size, 0.5 * size, trad_pix * oversize);
    }
    else
    {
        make_disk("pa1", size, size, 0.5 * size, 0.5 * size, trad_pix * oversize);
    }

    *out_total = arith_image_total("pa1") / (oversize * oversize);
    IDp1       = make_slopexy(
        "pp1",
        size,
        size,
        CORO_PI * xld / trad_pix,
        CORO_PI * yld / trad_pix);

    if (pha_name != NULL)
    {
        imageID id_pha = image_ID(pha_name, dcimg, dcnimg);
        if (id_pha != -1)
        {
            for (uint64_t ii = 0; ii < size2; ii++)
            {
                dcimg[IDp1].array.F[ii] += dcimg[id_pha].array.F[ii];
            }
        }
    }
}

/**
 * coro_complex_to_psf - Propagate exit-pupil complex field to normalized focal PSF
 * @cname:      Input pupil complex image name
 * @fname:      Temporary focal complex image name
 * @aname:      Temporary focal amplitude image name
 * @pname:      Temporary focal phase image name
 * @psfname:    Output PSF intensity image name
 * @norm_denom: Normalization divisor applied to every PSF pixel
 */
void coro_complex_to_psf(
    const char *cname,
    const char *fname,
    const char *aname,
    const char *pname,
    const char *psfname,
    double      norm_denom)
{
    uint32_t size  = CORONAGRAPHS_ARRAYSIZE;
    uint64_t size2 = (uint64_t) size * size;
    imageID  ID;

    permut(cname);
    do2dfft(cname, fname);
    permut(fname);
    mk_amph_from_complex(fname, aname, pname, 0);
    arith_image_mult(aname, aname, psfname);

    ID = image_ID(psfname, dcimg, dcnimg);
    delete_image_ID(fname, DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID(aname, DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID(pname, DELETE_IMAGE_ERRMODE_WARNING);

    for (uint64_t ii = 0; ii < size2; ii++)
    {
        dcimg[ID].array.F[ii] /= norm_denom;
    }
}

/**
 * coro_load_ref_pair - Load cached real/imaginary coronagraph reference maps if present
 * @ref: Reference pair descriptor (sets id_re, id_im, and needs_init)
 */
void coro_load_ref_pair(coro_ref_pair_t *ref)
{
    ref->needs_init = 0;

    ref->id_re = image_ID(ref->im_re, dcimg, dcnimg);
    if (ref->id_re == -1)
    {
        if (file_exists(ref->file_re) == 1)
        {
            load_fits(ref->file_re, ref->im_re, 1, &ref->id_re);
        }
        else
        {
            ref->needs_init = 1;
        }
    }

    ref->id_im = image_ID(ref->im_im, dcimg, dcnimg);
    if (ref->id_im == -1)
    {
        if (file_exists(ref->file_im) == 1)
        {
            load_fits(ref->file_im, ref->im_im, 1, &ref->id_im);
        }
        else
        {
            ref->needs_init = 1;
        }
    }
}

/**
 * coro_apply_or_save_ref_pair - Save reference field on first run or subtract it
 * @ref:  Reference pair descriptor
 * @id_c: Pupil complex image ID ("pc3")
 */
void coro_apply_or_save_ref_pair(
    coro_ref_pair_t *ref,
    imageID          id_c)
{
    uint32_t size  = CORONAGRAPHS_ARRAYSIZE;
    uint64_t size2 = (uint64_t) size * size;

    if (ref->needs_init == 1)
    {
        create_2Dimage_ID(ref->im_re, size, size, &ref->id_re);
        create_2Dimage_ID(ref->im_im, size, size, &ref->id_im);
        for (uint64_t ii = 0; ii < size2; ii++)
        {
            dcimg[ref->id_re].array.F[ii] = dcimg[id_c].array.CF[ii].re;
            dcimg[ref->id_im].array.F[ii] = dcimg[id_c].array.CF[ii].im;
        }
        save_fl_fits(ref->im_re, ref->file_re);
        save_fl_fits(ref->im_im, ref->file_im);
        printf("%s initialized, please rerun program.\n", ref->label);
        exit(0);
    }

    for (uint64_t ii = 0; ii < size2; ii++)
    {
        dcimg[id_c].array.CF[ii].re -= dcimg[ref->id_re].array.F[ii];
        dcimg[id_c].array.CF[ii].im -= dcimg[ref->id_im].array.F[ii];
    }
}
