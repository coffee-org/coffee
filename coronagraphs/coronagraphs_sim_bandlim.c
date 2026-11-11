/**
 * @file    coronagraphs_sim_bandlim.c
 * @brief   Band-limited and optical differentiation Lyot coronagraph simulators (ODC, BL8, BL4)
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "fft/fft.h"

#include "coronagraphs_internal.h"

/**
 * bandlim_load_or_init_mask - Load cached focal mask image or invoke initializer callback
 * @im_name:   Shared-memory mask image name
 * @file_path: FITS file path on disk
 * @init_fn:   Mask initialization callback if file is absent
 *
 * Return: Image identifier of loaded/created focal mask.
 */
static imageID bandlim_load_or_init_mask(
    const char *im_name,
    const char *file_path,
    int       (*init_fn)(void))
{
    imageID id_mask = image_ID(im_name, dcimg, dcnimg);
    if (id_mask == -1)
    {
        if (file_exists(file_path) == 1)
        {
            load_fits(file_path, im_name, 1, &id_mask);
        }
        else
        {
            init_fn();
        }
    }
    return image_ID(im_name, dcimg, dcnimg);
}

/**
 * bandlim_propagate_to_lyot - Propagate tilted pupil through focal mask to Lyot plane "pc3"
 * @xld:       Source X offset [lambda/D]
 * @yld:       Source Y offset [lambda/D]
 * @id_mask:   Focal plane amplitude mask image ID
 * @out_total: Output unapodized pupil flux integral
 *
 * Return: Image identifier of pre-stop Lyot complex field "pc3".
 */
static imageID bandlim_propagate_to_lyot(
    double  xld,
    double  yld,
    imageID id_mask,
    double *out_total)
{
    const uint64_t size2 = (uint64_t) CORONAGRAPHS_ARRAYSIZE * CORONAGRAPHS_ARRAYSIZE;

    coro_setup_pupil_tiptilt(xld, yld, 1, 1.0, "corphase", out_total);
    mk_complex_from_amph("pa1", "pp1", "pc2", 0);
    coro_fft2d_centered("pc2", "fc2");

    imageID id_fc2 = image_ID("fc2", dcimg, dcnimg);
    for (uint64_t ii = 0; ii < size2; ii++)
    {
        float m = dcimg[id_mask].array.F[ii];
        dcimg[id_fc2].array.CF[ii].re *= m;
        dcimg[id_fc2].array.CF[ii].im *= m;
    }

    coro_fft2d_centered("fc2", "pc3");
    return image_ID("pc3", dcimg, dcnimg);
}

/**
 * bandlim_finish_psf - Subtract/save Lyot reference, FFT to PSF, and delete scratch images
 * @ref:     Reference field descriptor
 * @id_pc3:  Lyot-stop complex field image ID
 * @psfname: Output PSF image name
 * @total:   Unapodized pupil flux integral
 */
static void bandlim_finish_psf(
    coro_ref_pair_t *ref,
    imageID          id_pc3,
    const char      *psfname,
    double           total)
{
    const double size2 = (double) CORONAGRAPHS_ARRAYSIZE * CORONAGRAPHS_ARRAYSIZE;

    coro_apply_or_save_ref_pair(ref, id_pc3);
    coro_complex_to_psf("pc3", "fc4", "fa4", "fp4", psfname, total * size2 * size2 * size2);

    delete_image_ID("pc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pc3", DELETE_IMAGE_ERRMODE_WARNING);
}

/**
 * coronagraph_simul_ODC - Simulate Optical Differentiation Coronagraph PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output PSF image name
 *
 * Return: 0 on completion.
 */
int coronagraph_simul_ODC(
    double      xld,
    double      yld,
    const char *psfname)
{
    const uint32_t size     = CORONAGRAPHS_ARRAYSIZE;
    const double   trad_pix = coro_pupil_radius_pix();
    char           file_re[200];
    char           file_im[200];
    char           refname[200];

    snprintf(file_re, sizeof(file_re), "%s/odcre.ref.%ul", CORONAGRAPHSDATALOCAL, size);
    snprintf(file_im, sizeof(file_im), "%s/odcim.ref.%ul", CORONAGRAPHSDATALOCAL, size);
    coro_ref_pair_t ref = {
        .im_re   = "refodcre",
        .file_re = file_re,
        .im_im   = "refodcim",
        .file_im = file_im,
        .label   = "ODC",
    };
    coro_load_ref_pair(&ref);
    if (ref.needs_init == 1)
    {
        printf("initialize ODC coronagraph\n");
    }

    snprintf(refname, sizeof(refname), "%s/odcm.ref.%ul", CORONAGRAPHSDATALOCAL, size);
    imageID id_mask = bandlim_load_or_init_mask("odcmref", refname, coronagraph_init_ODC);

    double  total  = 0.0;
    imageID id_pc3 = bandlim_propagate_to_lyot(xld, yld, id_mask, &total);

    for (uint32_t ii = 0; ii < size; ii++)
    {
        for (uint32_t jj = 0; jj < size; jj++)
        {
            double dx = 1.0 * ii - size / 2;
            double dy = 1.0 * jj - size / 2;
            if (sqrt(dx * dx + dy * dy) / trad_pix > g_coro.ODCMASK_eps)
            {
                dcimg[id_pc3].array.CF[jj * size + ii].re = 0.0f;
                dcimg[id_pc3].array.CF[jj * size + ii].im = 0.0f;
            }
        }
    }

    bandlim_finish_psf(&ref, id_pc3, psfname, total);
    return 0;
}

/**
 * coronagraph_simul_BL8 - Simulate 8th-order Band-Limited Lyot coronagraph PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output PSF image name
 *
 * Return: RETURN_SUCCESS on completion.
 */
errno_t coronagraph_simul_BL8(
    double      xld,
    double      yld,
    const char *psfname)
{
    const uint32_t size     = CORONAGRAPHS_ARRAYSIZE;
    const double   trad_pix = coro_pupil_radius_pix();
    char           file_re[200];
    char           file_im[200];
    char           refname[200];

    const char *im_re = (g_coro.BL8MODE == 0) ? "refbl8re" : "refbl8lre";
    const char *im_im = (g_coro.BL8MODE == 0) ? "refbl8im" : "refbl8lim";
    const char *tag_r = (g_coro.BL8MODE == 0) ? "bl8re" : "bl8lre";
    const char *tag_i = (g_coro.BL8MODE == 0) ? "bl8im" : "bl8lim";
    const char *im_m  = (g_coro.BL8MODE == 0) ? "bl8mref" : "bl8mrefl";
    const char *tag_m = (g_coro.BL8MODE == 0) ? "bl8m" : "bl8ml";

    snprintf(file_re, sizeof(file_re), "%s/%s.ref.%ul", CORONAGRAPHSDATALOCAL, tag_r, size);
    snprintf(file_im, sizeof(file_im), "%s/%s.ref.%ul", CORONAGRAPHSDATALOCAL, tag_i, size);
    coro_ref_pair_t ref = {
        .im_re   = im_re,
        .file_re = file_re,
        .im_im   = im_im,
        .file_im = file_im,
        .label   = "BL8",
    };
    coro_load_ref_pair(&ref);
    if (ref.needs_init == 1)
    {
        printf("initialize BL8 coronagraph\n");
    }

    snprintf(refname, sizeof(refname), "%s/%s.ref.%ul", CORONAGRAPHSDATALOCAL, tag_m, size);
    imageID id_mask = bandlim_load_or_init_mask(im_m, refname, coronagraph_init_BL8);

    double  total    = 0.0;
    imageID id_pc3   = bandlim_propagate_to_lyot(xld, yld, id_mask, &total);
    long    iioffset = (long) (g_coro.BL8MASK_eps * trad_pix);

    for (uint32_t ii = 0; ii < size; ii++)
    {
        for (uint32_t jj = 0; jj < size; jj++)
        {
            double dy  = 1.0 * jj - size / 2;
            int    out = 0;
            if (g_coro.BL8MODE == 0)
            {
                double dx = 1.0 * ii - size / 2;
                out       = (sqrt(dx * dx + dy * dy) / trad_pix > 1.0 - g_coro.BL8MASK_eps);
            }
            else
            {
                double dx1 = (1.0 * ii - size / 2) - iioffset;
                double dx2 = (1.0 * ii - size / 2) + iioffset;
                double r1  = sqrt(dx1 * dx1 + dy * dy) / trad_pix;
                double r2  = sqrt(dx2 * dx2 + dy * dy) / trad_pix;
                out        = ((r1 > 1.0) || (r2 > 1.0));
            }
            if (out)
            {
                dcimg[id_pc3].array.CF[jj * size + ii].re = 0.0f;
                dcimg[id_pc3].array.CF[jj * size + ii].im = 0.0f;
            }
        }
    }

    bandlim_finish_psf(&ref, id_pc3, psfname, total);
    return RETURN_SUCCESS;
}

/**
 * coronagraph_simul_BL4 - Simulate 4th-order Band-Limited Lyot coronagraph PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output PSF image name
 *
 * Return: 0 on completion.
 */
int coronagraph_simul_BL4(
    double      xld,
    double      yld,
    const char *psfname)
{
    const uint32_t size     = CORONAGRAPHS_ARRAYSIZE;
    const double   trad_pix = coro_pupil_radius_pix();
    char           file_re[200];
    char           file_im[200];
    char           refname[200];

    snprintf(file_re, sizeof(file_re), "%s/bl4re.ref.%ul", CORONAGRAPHSDATALOCAL, size);
    snprintf(file_im, sizeof(file_im), "%s/bl4im.ref.%ul", CORONAGRAPHSDATALOCAL, size);
    coro_ref_pair_t ref = {
        .im_re   = "refbl4re",
        .file_re = file_re,
        .im_im   = "refbl4im",
        .file_im = file_im,
        .label   = "BL4",
    };
    coro_load_ref_pair(&ref);
    if (ref.needs_init == 1)
    {
        printf("initialize BL4 coronagraph\n");
    }

    snprintf(refname, sizeof(refname), "%s/bl4m.ref.%ul", CORONAGRAPHSDATALOCAL, size);
    imageID id_mask = bandlim_load_or_init_mask("bl4mref", refname, coronagraph_init_BL4);

    double  total    = 0.0;
    imageID id_pc3   = bandlim_propagate_to_lyot(xld, yld, id_mask, &total);
    long    iioffset = (long) (g_coro.BL4MASK_eps * trad_pix + 0.2) + 1;

    for (uint32_t ii = 0; ii < size; ii++)
    {
        for (uint32_t jj = 0; jj < size; jj++)
        {
            double dx1 = (1.0 * ii - size / 2) - iioffset;
            double dx2 = (1.0 * ii - size / 2) + iioffset;
            double dy  = 1.0 * jj - size / 2;
            double r1  = (sqrt(dx1 * dx1 + dy * dy) + 1.0) / trad_pix;
            double r2  = (sqrt(dx2 * dx2 + dy * dy) + 1.0) / trad_pix;
            if ((r1 > 1.0) || (r2 > 1.0))
            {
                dcimg[id_pc3].array.CF[jj * size + ii].re = 0.0f;
                dcimg[id_pc3].array.CF[jj * size + ii].im = 0.0f;
            }
        }
    }

    bandlim_finish_psf(&ref, id_pc3, psfname, total);
    return 0;
}
