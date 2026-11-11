/**
 * @file    coronagraphs_sim_piaac.c
 * @brief   PIAA, PIAAC, and hybrid AIC-PIAAC coronagraph PSF simulators
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "coffee_compat.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "fft/fft.h"

#include "coronagraphs_internal.h"

/**
 * piaa_make_lowfs_psfs - Generate plus/minus defocused LOWFS PSFs ("psfWFSp", "psfWFSm")
 * @trad_pix:   Pupil radius in pixels
 * @norm_denom: PSF intensity normalization divisor
 */
static void piaa_make_lowfs_psfs(
    double trad_pix,
    double norm_denom)
{
    const uint32_t size = CORONAGRAPHS_ARRAYSIZE;

    copy_image_ID("pc3", "pc3p", 0);
    copy_image_ID("pc3", "pc3m", 0);
    imageID id_p = image_ID("pc3p", dcimg, dcnimg);
    imageID id_m = image_ID("pc3m", dcimg, dcnimg);
    imageID id_c = image_ID("pc3", dcimg, dcnimg);

    for (uint32_t ii = 0; ii < size; ii++)
    {
        for (uint32_t jj = 0; jj < size; jj++)
        {
            double dx = 1.0 * ii - size / 2;
            double dy = 1.0 * jj - size / 2;
            double r  = sqrt(dx * dx + dy * dy) / trad_pix;
            if (r < PIAAPPMASKRAD)
            {
                uint32_t idx  = jj * size + ii;
                double   re   = dcimg[id_c].array.CF[idx].re;
                double   im   = dcimg[id_c].array.CF[idx].im;
                double   amp  = sqrt(re * re + im * im);
                double   pha  = atan2(im, re);
                double   phap = pha + r * r * 6.0;
                double   pham = pha - r * r * 6.0;

                dcimg[id_p].array.CF[idx].re = (float) (amp * cos(phap));
                dcimg[id_p].array.CF[idx].im = (float) (amp * sin(phap));
                dcimg[id_m].array.CF[idx].re = (float) (amp * cos(pham));
                dcimg[id_m].array.CF[idx].im = (float) (amp * sin(pham));
            }
        }
    }

    coro_complex_to_psf("pc3p", "fc4p", "fa4p", "fp4p", "psfWFSp", norm_denom);
    delete_image_ID("pc3p", DELETE_IMAGE_ERRMODE_WARNING);

    coro_complex_to_psf("pc3m", "fc4m", "fa4m", "fp4m", "psfWFSm", norm_denom);
    delete_image_ID("pc3m", DELETE_IMAGE_ERRMODE_WARNING);
}

/**
 * coronagraph_simul_PIAA - Simulate lossless PIAA coronagraph PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output PSF image name
 *
 * Return: RETURN_SUCCESS on completion.
 */
errno_t coronagraph_simul_PIAA(
    double      xld,
    double      yld,
    const char *psfname)
{
    const uint32_t size     = CORONAGRAPHS_ARRAYSIZE;
    const double   size2    = (double) size * size;
    const double   trad_pix = coro_pupil_radius_pix();
    double         total    = 0.0;

    if (g_coro.initPIAA == 0)
    {
        coronagraph_init_PIAA(NULL);
        g_coro.initPIAA = 1;
    }

    coro_setup_pupil_tiptilt(xld, yld, 1, g_coro.PIAAOVERSIZE, "corphase", &total);
    coronagraphs_PIAA_apodize_beam("pa1", "pp1", "pa2", "pp2");

    mk_complex_from_amph("pa2", "pp2", "pc2", 0);
    coro_fft2d_centered("pc2", "fc2");
    imageID id_fc2 = image_ID("fc2", dcimg, dcnimg);

    if ((g_coro.PIAAFPMASK == 1) || (g_coro.PIAALOWFS == 1))
    {
        for (uint32_t ii = 0; ii < size; ii++)
        {
            for (uint32_t jj = 0; jj < size; jj++)
            {
                double dx   = 1.0 * ii - size / 2;
                double dy   = 1.0 * jj - size / 2;
                double r    = sqrt(dx * dx + dy * dy) * CORONAGRAPHS_PIXSCALE;
                int    zero = 0;
                if ((g_coro.PIAAFPMASK == 1) && (r < g_coro.PIAAFPMASKRAD))
                {
                    zero = 1;
                }
                if ((g_coro.PIAALOWFS == 1) &&
                    ((r > g_coro.PIAAFPMASKRAD * 1.2) || (r < 1.2 * g_coro.PIAAFPMASKRAD / 4.0)))
                {
                    zero = 1;
                }
                if (zero)
                {
                    dcimg[id_fc2].array.CF[jj * size + ii].re = 0.0f;
                    dcimg[id_fc2].array.CF[jj * size + ii].im = 0.0f;
                }
            }
        }
    }

    coro_fft2d_centered("fc2", "pc3");
    imageID id_pc3 = image_ID("pc3", dcimg, dcnimg);
    for (uint32_t ii = 0; ii < size; ii++)
    {
        for (uint32_t jj = 0; jj < size; jj++)
        {
            double dx = 1.0 * ii - size / 2;
            double dy = 1.0 * jj - size / 2;
            if (sqrt(dx * dx + dy * dy) / trad_pix > PIAAPPMASKRAD)
            {
                dcimg[id_pc3].array.CF[jj * size + ii].re = 0.0f;
                dcimg[id_pc3].array.CF[jj * size + ii].im = 0.0f;
            }
        }
    }

    double norm_denom = total * size2 * size2 * size2;
    if (g_coro.PIAALOWFS == 1)
    {
        piaa_make_lowfs_psfs(trad_pix, norm_denom);
    }

    coro_complex_to_psf("pc3", "fc4", "fa4", "fp4", psfname, norm_denom);

    delete_image_ID("pc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pc3", DELETE_IMAGE_ERRMODE_WARNING);
    return RETURN_SUCCESS;
}

/**
 * piaac_cheat_subtract - Subtract focal-mask-reflected complex beam from entrance pupil
 * @trad_pix: Pupil radius in pixels
 */
static void piaac_cheat_subtract(double trad_pix)
{
    const uint32_t size = CORONAGRAPHS_ARRAYSIZE;
    imageID        id1  = image_ID("pa1", dcimg, dcnimg);
    imageID        id2  = image_ID("pp1", dcimg, dcnimg);
    imageID        id3  = image_ID("pa4", dcimg, dcnimg);
    imageID        id4  = image_ID("pp4", dcimg, dcnimg);

    printf("------------ MASK SUBTRACTION ----------\n");
    for (uint32_t ii = 10; ii < size - 10; ii++)
    {
        for (uint32_t jj = 10; jj < size - 10; jj++)
        {
            long   idx1 = jj * size + ii;
            long   idx2 = (size - jj) * size + (size - ii);
            double re1  = dcimg[id1].array.F[idx2] * cos(dcimg[id2].array.F[idx2]);
            double im1  = dcimg[id1].array.F[idx2] * sin(dcimg[id2].array.F[idx2]);
            double re2  = dcimg[id3].array.F[idx1] * cos(dcimg[id4].array.F[idx1]);
            double im2  = dcimg[id3].array.F[idx1] * sin(dcimg[id4].array.F[idx1]);
            double re   = re1 * size * size - re2;
            double im   = im1 * size * size - im2;

            double dx = 1.0 * ii - size / 2;
            double dy = 1.0 * jj - size / 2;
            if (sqrt(dx * dx + dy * dy) / trad_pix > PIAAPPMASKRAD)
            {
                re = 0.0;
                im = 0.0;
            }

            dcimg[id3].array.F[idx1] = (float) sqrt(re * re + im * im);
            dcimg[id4].array.F[idx1] = (float) atan2(im, re);
        }
    }
}

/**
 * coronagraph_simul_PIAAC - Simulate PIAA + inverse-PIAA coronagraph PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output PSF image name
 *
 * Return: RETURN_SUCCESS on completion.
 */
errno_t coronagraph_simul_PIAAC(
    double      xld,
    double      yld,
    const char *psfname)
{
    const uint32_t size     = CORONAGRAPHS_ARRAYSIZE;
    const double   size2    = (double) size * size;
    const double   trad_pix = coro_pupil_radius_pix();
    const int      cheat    = 0;
    double         total    = 0.0;

    if (g_coro.initPIAA == 0)
    {
        coronagraph_init_PIAA(NULL);
        g_coro.initPIAA = 1;
    }

    coro_setup_pupil_tiptilt(xld, yld, 1, g_coro.PIAAOVERSIZE, "corphase", &total);
    g_coro.piaaconfdirection = 0;
    coronagraphs_PIAA_apodize_beam("pa1", "pp1", "pa2", "pp2");

    mk_complex_from_amph("pa2", "pp2", "pc2", 0);
    coro_fft2d_centered("pc2", "fc2");
    imageID id_fc2 = image_ID("fc2", dcimg, dcnimg);

    for (uint32_t ii = 0; ii < size; ii++)
    {
        for (uint32_t jj = 0; jj < size; jj++)
        {
            double dx   = 1.0 * ii - size / 2;
            double dy   = 1.0 * jj - size / 2;
            double r    = sqrt(dx * dx + dy * dy) * CORONAGRAPHS_PIXSCALE;
            int    mask = (cheat == 1) ? (r > g_coro.PIAAFPMASKRAD) : (r < g_coro.PIAAFPMASKRAD);
            if (mask)
            {
                dcimg[id_fc2].array.CF[jj * size + ii].re = 0.0f;
                dcimg[id_fc2].array.CF[jj * size + ii].im = 0.0f;
            }
        }
    }

    coro_fft2d_centered("fc2", "pc3");
    imageID id_pc3 = image_ID("pc3", dcimg, dcnimg);
    for (uint32_t ii = 0; ii < size; ii++)
    {
        for (uint32_t jj = 0; jj < size; jj++)
        {
            double dx = 1.0 * ii - size / 2;
            double dy = 1.0 * jj - size / 2;
            if (sqrt(dx * dx + dy * dy) / trad_pix > PIAAPPMASKRAD)
            {
                dcimg[id_pc3].array.CF[jj * size + ii].re = 0.0f;
                dcimg[id_pc3].array.CF[jj * size + ii].im = 0.0f;
            }
        }
    }

    mk_amph_from_complex("pc3", "pa3", "pp3", 0);
    coro_piaa_inverse_unwrap("pa3", "pp3", "pa4", "pp4");
    if (cheat == 1)
    {
        piaac_cheat_subtract(trad_pix);
    }

    mk_complex_from_amph("pa4", "pp4", "pc4", 0);
    coro_complex_to_psf("pc4", "fc4", "fa4", "fp4", psfname, total * size2 * size2 * size2);

    delete_image_ID("pc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa3", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa4", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp3", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp4", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pc3", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pc4", DELETE_IMAGE_ERRMODE_WARNING);
    return RETURN_SUCCESS;
}

/**
 * aic_piaac_combine_output - Compute inner AIC channels and tile with outer PIAAC PSF
 * @id_aicfield: Inner focal complex field ("fc2_aic")
 * @total:       Unapodized pupil flux integral
 * @psfname:     Output composite PSF image name
 */
static void aic_piaac_combine_output(
    imageID     id_aicfield,
    double      total,
    const char *psfname)
{
    const uint32_t size  = CORONAGRAPHS_ARRAYSIZE;
    const double   size2 = (double) size * size;
    imageID        id_aicc;
    imageID        id_aicd;
    imageID        id_out;

    create_2Dimage_ID("psfAICc", size, size, &id_aicc);
    create_2Dimage_ID("psfAICd", size, size, &id_aicd);
    for (uint32_t ii = 1; ii < size; ii++)
    {
        for (uint32_t jj = 1; jj < size; jj++)
        {
            uint32_t idx1 = jj * size + ii;
            uint32_t idx2 = (size - jj) * size + (size - ii);
            double   r11  = dcimg[id_aicfield].array.CF[idx1].re;
            double   r12  = dcimg[id_aicfield].array.CF[idx2].re;
            double   i11  = dcimg[id_aicfield].array.CF[idx1].im;
            double   i12  = dcimg[id_aicfield].array.CF[idx2].im;

            dcimg[id_aicd].array.F[idx1] =
                (float) (0.25 * ((r11 - r12) * (r11 - r12) + (i11 - i12) * (i11 - i12)) /
                         (total * size2));
            dcimg[id_aicc].array.F[idx1] =
                (float) (0.25 * ((r11 + r12) * (r11 + r12) + (i11 + i12) * (i11 + i12)) /
                         (total * size2));
        }
    }
    delete_image_ID("fc2_aic", DELETE_IMAGE_ERRMODE_WARNING);

    imageID id_piaac = image_ID("psfPIAAC", dcimg, dcnimg);
    create_2Dimage_ID(psfname, size, size, &id_out);

    for (uint32_t ii = 0; ii < size; ii++)
    {
        for (uint32_t jj = 0; jj < size / 2; jj++)
        {
            dcimg[id_out].array.F[jj * size + ii] =
                dcimg[id_piaac].array.F[(jj + size / 4) * size + ii];
        }
    }
    for (uint32_t ii = 0; ii < size / 2; ii++)
    {
        for (uint32_t jj = size / 2; jj < size; jj++)
        {
            dcimg[id_out].array.F[jj * size + ii] =
                dcimg[id_aicd].array.F[(jj - size / 4) * size + (ii + size / 4)];
        }
    }
    for (uint32_t ii = size / 2; ii < size; ii++)
    {
        for (uint32_t jj = size / 2; jj < size; jj++)
        {
            dcimg[id_out].array.F[jj * size + ii] =
                dcimg[id_aicc].array.F[(jj - size / 4) * size + (ii - size / 4)];
        }
    }

    delete_image_ID("psfPIAAC", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("psfAICc", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("psfAICd", DELETE_IMAGE_ERRMODE_WARNING);
}

/**
 * coronagraph_simul_AIC_PIAAC - Simulate hybrid inner-AIC + outer-PIAAC coronagraph PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output composite PSF image name
 *
 * Return: RETURN_SUCCESS on completion.
 */
errno_t coronagraph_simul_AIC_PIAAC(
    double      xld,
    double      yld,
    const char *psfname)
{
    const uint32_t size     = CORONAGRAPHS_ARRAYSIZE;
    const double   size2    = (double) size * size;
    const double   trad_pix = coro_pupil_radius_pix();
    double         total    = 0.0;
    imageID        id_aicfield;

    if (g_coro.initPIAA == 0)
    {
        coronagraph_init_PIAA(NULL);
        g_coro.initPIAA = 1;
    }

    coro_setup_pupil_tiptilt(xld, yld, 1, g_coro.PIAAOVERSIZE, "corphase", &total);
    g_coro.piaaconfdirection = 0;
    coronagraphs_PIAA_apodize_beam("pa1", "pp1", "pa2", "pp2");

    mk_complex_from_amph("pa2", "pp2", "pc2", 0);
    coro_fft2d_centered("pc2", "fc2");

    create_2DCimage_ID("fc2_aic", size, size, &id_aicfield);
    imageID id_fc2 = image_ID("fc2", dcimg, dcnimg);
    for (uint32_t ii = 0; ii < size; ii++)
    {
        for (uint32_t jj = 0; jj < size; jj++)
        {
            double dx = 1.0 * ii - size / 2;
            double dy = 1.0 * jj - size / 2;
            if (sqrt(dx * dx + dy * dy) * CORONAGRAPHS_PIXSCALE < PIAAAIC_FIELDLIMIT)
            {
                uint32_t idx                       = jj * size + ii;
                dcimg[id_aicfield].array.CF[idx]   = dcimg[id_fc2].array.CF[idx];
                dcimg[id_fc2].array.CF[idx].re     = 0.0f;
                dcimg[id_fc2].array.CF[idx].im     = 0.0f;
            }
        }
    }

    coro_fft2d_centered("fc2", "pc3");
    imageID id_pc3 = image_ID("pc3", dcimg, dcnimg);
    for (uint32_t ii = 0; ii < size; ii++)
    {
        for (uint32_t jj = 0; jj < size; jj++)
        {
            double dx = 1.0 * ii - size / 2;
            double dy = 1.0 * jj - size / 2;
            if (sqrt(dx * dx + dy * dy) / trad_pix > PIAAPPMASKRAD)
            {
                dcimg[id_pc3].array.CF[jj * size + ii].re = 0.0f;
                dcimg[id_pc3].array.CF[jj * size + ii].im = 0.0f;
            }
        }
    }

    mk_amph_from_complex("pc3", "pa3", "pp3", 0);
    g_coro.piaaconfdirection = 1;
    coronagraphs_PIAA_apodize_beam("pa3", "pp3", "pa4", "pp4");
    g_coro.piaaconfdirection = 0;
    mk_complex_from_amph("pa4", "pp4", "pc4", 0);

    coro_complex_to_psf("pc4", "fc4", "fa4", "fp4", "psfPIAAC", total * size2 * size2 * size2);

    delete_image_ID("pc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa3", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa4", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp3", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp4", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pc3", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pc4", DELETE_IMAGE_ERRMODE_WARNING);

    aic_piaac_combine_output(id_aicfield, total, psfname);
    return RETURN_SUCCESS;
}
