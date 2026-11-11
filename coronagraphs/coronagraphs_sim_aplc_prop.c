/**
 * @file    coronagraphs_sim_aplc_prop.c
 * @brief   Multi-step APLC/PIAACMC focal-plane propagation steps and final PSF formation
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "coffee_compat.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "fft/fft.h"

#include "coronagraphs_internal.h"

/**
 * aplc_step_fft - Propagate "pc1" across one APLC stage via full-array 2D FFT and "LyotMask"
 * @write_files: Non-zero to save intermediate FITS files
 */
static void aplc_step_fft(int write_files)
{
    const uint32_t size   = CORONAGRAPHS_ARRAYSIZE;
    const uint64_t size2  = (uint64_t) size * size;
    imageID        id_fpm = image_ID("aplcfpm", dcimg, dcnimg);

    permut("pc1");
    do2dfft("pc1", "fc1");
    delete_image_ID("pc1", DELETE_IMAGE_ERRMODE_WARNING);
    permut("fc1");

    imageID id_fc1 = image_ID("fc1", dcimg, dcnimg);
    for (uint64_t ii = 0; ii < size2; ii++)
    {
        float m = dcimg[id_fpm].array.F[ii];
        dcimg[id_fc1].array.CF[ii].re *= m;
        dcimg[id_fc1].array.CF[ii].im *= m;
    }

    mk_amph_from_complex("fc1", "fc1_a", "fc1_p", 0);
    if (write_files == 1)
    {
        save_fl_fits("aplcfpm", "aplcfpm.fits");
        save_fl_fits("fc1_a", "fc1_a.fits");
        save_fl_fits("fc1_p", "fc1_p.fits");
    }
    delete_image_ID("fc1_a", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fc1_p", DELETE_IMAGE_ERRMODE_WARNING);

    permut("fc1");
    do2dfft("fc1", "pc1");
    delete_image_ID("fc1", DELETE_IMAGE_ERRMODE_WARNING);
    permut("pc1");

    if (write_files == 1)
    {
        mk_amph_from_complex("pc1", "pc1a", "pc1p", 0);
        save_fl_fits("pc1a", "pc1a.fits");
        save_fl_fits("pc1p", "pc1p.fits");
        delete_image_ID("pc1a", DELETE_IMAGE_ERRMODE_WARNING);
        delete_image_ID("pc1p", DELETE_IMAGE_ERRMODE_WARNING);
    }

    imageID id_pc1 = image_ID("pc1", dcimg, dcnimg);
    imageID id_lm  = image_ID("LyotMask", dcimg, dcnimg);
    for (uint32_t ii = 0; ii < size; ii++)
    {
        for (uint32_t jj = 0; jj < size; jj++)
        {
            uint32_t idx = jj * size + ii;
            if ((ii < 5) || (jj < 5) || (ii > size - 5) || (jj > size - 5))
            {
                dcimg[id_pc1].array.CF[idx].re = 0.0f;
                dcimg[id_pc1].array.CF[idx].im = 0.0f;
            }
            else
            {
                float scale = dcimg[id_lm].array.F[(size - jj) * size + (size - ii)] / size2;
                dcimg[id_pc1].array.CF[idx].re *= scale;
                dcimg[id_pc1].array.CF[idx].im *= scale;
            }
        }
    }

    if (write_files == 1)
    {
        mk_amph_from_complex("pc1", "pc1am", "pc1pm", 0);
        save_fl_fits("pc1am", "pc1am.fits");
        save_fl_fits("pc1pm", "pc1pm.fits");
        delete_image_ID("pc1am", DELETE_IMAGE_ERRMODE_WARNING);
        delete_image_ID("pc1pm", DELETE_IMAGE_ERRMODE_WARNING);
    }
}

/**
 * aplc_step_dft - Propagate "pc1" across one APLC stage via zoomed focal-plane DFT
 * @write_files:  Non-zero to save intermediate FITS files
 * @mask_ca_mult: Focal-plane complex amplitude multiplier
 */
static void aplc_step_dft(
    int    write_files,
    double mask_ca_mult)
{
    const uint32_t size  = CORONAGRAPHS_ARRAYSIZE;
    const uint64_t size2 = (uint64_t) size * size;

    printf("STARTING DFT  [%f %f]\n", mask_ca_mult, -(1.0 - mask_ca_mult));
    fflush(stdout);

    if (write_files == 1)
    {
        save_fl_fits("aplcfpm", "TEST_aplcfpm.fits");
    }

    imageID id_fpm = image_ID("aplcfpm", dcimg, dcnimg);
    imageID id_fpm1;
    create_2Dimage_ID("aplcfpm1", size, size, &id_fpm1);
    for (uint64_t ii = 0; ii < size2; ii++)
    {
        dcimg[id_fpm1].array.F[ii] = 1.0f - dcimg[id_fpm].array.F[ii];
    }
    if (write_files == 1)
    {
        save_fl_fits("aplcfpm1", "TEST_aplcfpm1.fits");
    }

    imageID id_fpma = image_ID("aplcfpm_a", dcimg, dcnimg);
    imageID id_fpmp = image_ID("aplcfpm_p", dcimg, dcnimg);
    if ((id_fpma != -1) && (id_fpmp != -1))
    {
        printf("============= USING CUSTOM COMPLEX FOCAL PLANE MASK ============\n");
        mk_complex_from_amph("aplcfpm_a", "aplcfpm_p", "aplcfpm_c", 0);
        imageID id_c = image_ID("aplcfpm_c", dcimg, dcnimg);
        for (uint64_t ii = 0; ii < size2; ii++)
        {
            dcimg[id_c].array.CF[ii].re = 1.0f - dcimg[id_c].array.CF[ii].re;
            dcimg[id_c].array.CF[ii].im = -dcimg[id_c].array.CF[ii].im;
        }
        fft_DFTinsertFPM("pc1", "aplcfpm_c", g_coro.DFTZFACTOR, "pc2", NULL);
        delete_image_ID("aplcfpm_c", DELETE_IMAGE_ERRMODE_WARNING);
    }
    else
    {
        fft_DFTinsertFPM_re("pc1", "aplcfpm1", g_coro.DFTZFACTOR, "pc2", NULL);
    }
    delete_image_ID("aplcfpm1", DELETE_IMAGE_ERRMODE_WARNING);

    imageID id1 = image_ID("pc1", dcimg, dcnimg);
    imageID id2 = image_ID("pc2", dcimg, dcnimg);
    for (uint64_t ii = 0; ii < size2; ii++)
    {
        dcimg[id1].array.CF[ii].re -= dcimg[id2].array.CF[ii].re;
        dcimg[id1].array.CF[ii].im -= dcimg[id2].array.CF[ii].im;
    }
    delete_image_ID("pc2", DELETE_IMAGE_ERRMODE_WARNING);

    if (write_files == 1)
    {
        mk_amph_from_complex("pc1", "pc1a", "pc1p", 0);
        save_fl_fits("pc1a", "pc1a_01.fits");
        save_fl_fits("pc1p", "pc1p_01.fits");
        delete_image_ID("pc1a", DELETE_IMAGE_ERRMODE_WARNING);
        delete_image_ID("pc1p", DELETE_IMAGE_ERRMODE_WARNING);
    }
}

/**
 * coro_aplc_propagate_steps - Run NB_APLC_STEP pupil-to-pupil coronagraphic mask stages
 * @write_files:  Non-zero to write intermediate diagnostic FITS files
 * @mask_ca_mult: Focal-plane complex amplitude multiplier
 */
void coro_aplc_propagate_steps(
    int    write_files,
    double mask_ca_mult)
{
    for (long iter = 0; iter < g_coro.NB_APLC_STEP; iter++)
    {
        if (g_coro.useDFT == 0)
        {
            aplc_step_fft(write_files);
        }
        else
        {
            aplc_step_dft(write_files, mask_ca_mult);
        }
    }
}

/**
 * aplc_symmetrize_focal_flip - Symmetrize complex focal plane field under 180-deg inversion
 * @im_name: Complex focal field image name
 */
static void aplc_symmetrize_focal_flip(const char *im_name)
{
    const uint32_t size = CORONAGRAPHS_ARRAYSIZE;
    imageID        id   = image_ID(im_name, dcimg, dcnimg);

    for (uint32_t ii = 1; ii < size; ii++)
    {
        for (uint32_t jj = 1; jj < size / 2 + 1; jj++)
        {
            uint32_t idx1 = jj * size + ii;
            uint32_t idx2 = (size - jj) * size + (size - ii);
            float    re   = 0.5f * (dcimg[id].array.CF[idx1].re + dcimg[id].array.CF[idx2].re);
            float    im   = 0.5f * (dcimg[id].array.CF[idx1].im + dcimg[id].array.CF[idx2].im);

            dcimg[id].array.CF[idx1].re = re;
            dcimg[id].array.CF[idx1].im = im;
            dcimg[id].array.CF[idx2].re = re;
            dcimg[id].array.CF[idx2].im = im;
        }
    }
}

/**
 * aplc_finish_inverse_piaa - Apply inverse PIAA unwrapping, stop, and final FFT to PSF
 * @psfname:     Output PSF image name
 * @total0:      Unapodized entrance pupil energy
 * @write_files: Non-zero to write intermediate diagnostic FITS files
 */
static void aplc_finish_inverse_piaa(
    const char *psfname,
    double      total0,
    int         write_files)
{
    const uint32_t size     = CORONAGRAPHS_ARRAYSIZE;
    const uint64_t size2    = (uint64_t) size * size;
    const double   trad_pix = coro_pupil_radius_pix();

    mk_amph_from_complex("pc1", "pa3", "pp3", 0);
    if (write_files == 1)
    {
        save_fl_fits("pa3", "pa3.fits");
    }

    coro_piaa_inverse_unwrap("pa3", "pp3", "pa4", "pp4");
    printf(" CENT OBS 0 = %f x %f -------------------\n", g_coro.APLC_CentOBS0, trad_pix);

    imageID id3 = image_ID("pa4", dcimg, dcnimg);
    imageID id4 = image_ID("pp4", dcimg, dcnimg);
    for (uint32_t ii = 0; ii < size; ii++)
    {
        for (uint32_t jj = 0; jj < size; jj++)
        {
            long   idx = jj * size + ii;
            double re  = dcimg[id3].array.F[idx] * cos(dcimg[id4].array.F[idx]);
            double im  = dcimg[id3].array.F[idx] * sin(dcimg[id4].array.F[idx]);
            double dx  = 1.0 * ii - size / 2;
            double dy  = 1.0 * jj - size / 2;
            double r   = sqrt(dx * dx + dy * dy);

            if ((r > PIAAPPMASKRAD * trad_pix) || (r < g_coro.APLC_CentOBS0 * trad_pix + 3.0))
            {
                re = 0.0;
                im = 0.0;
            }
            dcimg[id3].array.F[idx] = (float) sqrt(re * re + im * im);
            dcimg[id4].array.F[idx] = (float) atan2(im, re);
        }
    }

    mk_complex_from_amph("pa4", "pp4", "pc4", 0);
    if (write_files == 1)
    {
        save_fl_fits("pa4", "pa4.fits");
        save_fl_fits("pp4", "pp4.fits");
    }

    coro_fft2d_centered("pc4", "fc4");
    if (g_coro.APLC_FLIP == 1)
    {
        printf("-------- FLIP --------\n");
        aplc_symmetrize_focal_flip("fc4");
    }
    mk_amph_from_complex("fc4", "fa4", "fp4", 0);

    arith_image_mult("fa4", "fa4", psfname);
    imageID id_psf = image_ID(psfname, dcimg, dcnimg);
    for (uint64_t ii = 0; ii < size2; ii++)
    {
        dcimg[id_psf].array.F[ii] /= (float) (size2 * total0);
    }

    delete_image_ID("pa3", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp3", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pc4", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa4", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp4", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fc4", DELETE_IMAGE_ERRMODE_WARNING);
    if (write_files == 1)
    {
        save_fl_fits("fa4", "fa4.fits");
        save_fl_fits("fp4", "fp4.fits");
    }
    delete_image_ID("fa4", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fp4", DELETE_IMAGE_ERRMODE_WARNING);
}

/**
 * coro_aplc_finish_psf - Convert post-Lyot complex field "pc1" into final APLC/PIAACMC PSF
 * @psfname:     Output PSF image name
 * @total0:      Unapodized entrance pupil energy
 * @write_files: Non-zero to write intermediate diagnostic FITS files
 */
void coro_aplc_finish_psf(
    const char *psfname,
    double      total0,
    int         write_files)
{
    if (g_coro.APLC_PIAA == 2)
    {
        aplc_finish_inverse_piaa(psfname, total0, write_files);
        return;
    }

    const uint32_t size  = CORONAGRAPHS_ARRAYSIZE;
    const uint64_t size2 = (uint64_t) size * size;

    do2dfft("pc1", "fc2");
    permut("fc2");
    if (g_coro.APLC_FLIP == 1)
    {
        aplc_symmetrize_focal_flip("fc2");
    }

    mk_amph_from_complex("fc2", "fa2", "fp2", 0);
    imageID id_fa2 = image_ID("fa2", dcimg, dcnimg);
    delete_image_ID("fp2", DELETE_IMAGE_ERRMODE_WARNING);

    imageID id_out;
    create_2Dimage_ID(psfname, size, size, &id_out);
    for (uint64_t ii = 0; ii < size2; ii++)
    {
        float a                   = dcimg[id_fa2].array.F[ii];
        dcimg[id_out].array.F[ii] = (float) (a / size2 / total0 * a);
    }

    delete_image_ID("fc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fa2", DELETE_IMAGE_ERRMODE_WARNING);
}
