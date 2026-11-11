/**
 * @file    coronagraphs_sim_vortex_rrpm.c
 * @brief   Roddier & Roddier Phase Mask (RRPM) and Optical Vortex Coronagraph (OVC) simulators
 */

#include <math.h>
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
#include "statistic/statistic.h"

#include "coronagraphs_internal.h"

/**
 * rrpm_load_masks - Load or initialize RRPM focal phase disk and pupil apodizer maps
 * @out_fm: Output focal phase mask image ID
 * @out_pm: Output pupil amplitude mask image ID
 */
static void rrpm_load_masks(
    imageID *out_fm,
    imageID *out_pm)
{
    const uint32_t size = CORONAGRAPHS_ARRAYSIZE;
    char           refname[200];

    snprintf(g_coro.PIAAAPO_FNAME, sizeof(g_coro.PIAAAPO_FNAME), "rrpmref_pm.ref.%ul.prof", size);
    snprintf(
        g_coro.PIAAAPODIZE_2DAPOFNAME,
        sizeof(g_coro.PIAAAPODIZE_2DAPOFNAME),
        "rrpmref_pm.ref.%ul",
        size);

    snprintf(refname, sizeof(refname), "%s/rrpmref_fm.ref.%ul", CORONAGRAPHSDATALOCAL, size);
    imageID id_fm = image_ID("rrpmreffm", dcimg, dcnimg);
    if (id_fm == -1)
    {
        if (file_exists(refname) == 1)
        {
            load_fits(refname, "rrpmreffm", 1, &id_fm);
        }
        else
        {
            coronagraph_init_RRPM();
        }
    }
    *out_fm = image_ID("rrpmreffm", dcimg, dcnimg);

    snprintf(refname, sizeof(refname), "%s/rrpmref_pm.ref.%ul", CORONAGRAPHSDATALOCAL, size);
    imageID id_pm = image_ID("rrpmrefpm", dcimg, dcnimg);
    if (id_pm == -1)
    {
        if (file_exists(refname) == 1)
        {
            load_fits(refname, "rrpmrefpm", 1, &id_pm);
        }
        else
        {
            coronagraph_init_RRPM();
        }
    }
    *out_pm = image_ID("rrpmrefpm", dcimg, dcnimg);
}

/**
 * rrpm_apodize_pupil - Apply conventional or lossless PIAA pupil apodization for RRPM
 * @id_a1: Entrance pupil amplitude image ID
 * @id_p1: Entrance pupil phase image ID
 * @id_pm: RRPM polynomial pupil apodizer image ID
 *
 * Return: Flux throughput scaling factor (1.0 for conventional, >1.0 for PIAA).
 */
static double rrpm_apodize_pupil(
    imageID id_a1,
    imageID id_p1,
    imageID id_pm)
{
    const uint32_t size   = CORONAGRAPHS_ARRAYSIZE;
    const uint64_t size2  = (uint64_t) size * size;
    double         factor = 1.0;
    imageID        id_a2;

    if (g_coro.RRPM_PIAA == 0)
    {
        create_2Dimage_ID("pa2", size, size, &id_a2);
        for (uint64_t ii = 0; ii < size2; ii++)
        {
            dcimg[id_a2].array.F[ii] = dcimg[id_a1].array.F[ii] * dcimg[id_pm].array.F[ii];
        }
        return 1.0;
    }

    g_coro.PIAAAPO_NBPOINTS = 80;
    if (g_coro.initPIAA == 0)
    {
        coronagraph_init_PIAA(NULL);
        g_coro.initPIAA = 1;
    }

    coronagraphs_PIAA_apodize_beam("pa1", "pp1", "pa2", "pp2");
    id_a2         = image_ID("pa2", dcimg, dcnimg);
    imageID id_p2 = image_ID("pp2", dcimg, dcnimg);
    for (uint64_t ii = 0; ii < size2; ii++)
    {
        dcimg[id_p1].array.F[ii] = dcimg[id_p2].array.F[ii];
    }

    double tmp1 = 0.0;
    double tmp2 = 0.0;
    for (uint64_t ii = 0; ii < size2; ii++)
    {
        double v1 = dcimg[id_a1].array.F[ii];
        double v2 = dcimg[id_pm].array.F[ii];
        tmp1 += v1 * v1;
        if (v1 > 0.5)
        {
            tmp2 += v2 * v2;
        }
        dcimg[id_a2].array.F[ii] = (float) (v1 * v2);
    }
    factor = tmp1 / tmp2;
    delete_image_ID("pp2", DELETE_IMAGE_ERRMODE_WARNING);
    return factor;
}

/**
 * coronagraph_simul_RRPM - Simulate Roddier & Roddier Phase Mask coronagraph PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output PSF image name
 *
 * Return: RETURN_SUCCESS on completion.
 */
errno_t coronagraph_simul_RRPM(
    double      xld,
    double      yld,
    const char *psfname)
{
    const uint32_t size  = CORONAGRAPHS_ARRAYSIZE;
    const uint64_t size2 = (uint64_t) size * size;
    imageID        id_fm;
    imageID        id_pm;
    char           file_re[200];
    char           file_im[200];

    rrpm_load_masks(&id_fm, &id_pm);

    snprintf(file_re, sizeof(file_re), "%s/rrpmre.ref.%ul", CORONAGRAPHSDATALOCAL, size);
    snprintf(file_im, sizeof(file_im), "%s/rrpmim.ref.%ul", CORONAGRAPHSDATALOCAL, size);
    coro_ref_pair_t ref = {
        .im_re   = "refrrpmre",
        .file_re = file_re,
        .im_im   = "refrrpmim",
        .file_im = file_im,
        .label   = "RRPM",
    };
    coro_load_ref_pair(&ref);

    double x = xld;
    double y = yld;
    if (ref.needs_init == 1)
    {
        printf("initialize RRPM coronagraph\n");
        x = 0.0;
        y = 0.0;
    }

    double total = 0.0;
    coro_setup_pupil_tiptilt(x, y, 1, 1.0, "corphase", &total);
    imageID id_a1  = image_ID("pa1", dcimg, dcnimg);
    imageID id_p1  = image_ID("pp1", dcimg, dcnimg);
    double  factor = rrpm_apodize_pupil(id_a1, id_p1, id_pm);

    mk_complex_from_amph("pa2", "pp1", "pc2", 0);
    coro_fft2d_centered("pc2", "fc2");
    imageID id_fc2 = image_ID("fc2", dcimg, dcnimg);
    for (uint64_t ii = 0; ii < size2; ii++)
    {
        float fm = dcimg[id_fm].array.F[ii];
        dcimg[id_fc2].array.CF[ii].re *= fm;
        dcimg[id_fc2].array.CF[ii].im *= fm;
    }

    coro_fft2d_centered("fc2", "pc3");
    imageID id_pc3 = image_ID("pc3", dcimg, dcnimg);
    for (uint64_t ii = 0; ii < size2; ii++)
    {
        if (dcimg[id_a1].array.F[ii] < 0.999f)
        {
            dcimg[id_pc3].array.CF[ii].re = 0.0f;
            dcimg[id_pc3].array.CF[ii].im = 0.0f;
        }
    }

    coro_apply_or_save_ref_pair(&ref, id_pc3);
    double s2 = (double) size2;
    coro_complex_to_psf("pc3", "fc4", "fa4", "fp4", psfname, total * s2 * s2 * s2 / factor);

    delete_image_ID("pc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pc3", DELETE_IMAGE_ERRMODE_WARNING);
    return RETURN_SUCCESS;
}

/**
 * coronagraph_RRPM_optimize - Random-walk optimization of RRPM mask radius and polynomial
 *
 * Return: 0 (runs unbounded random walk).
 */
int coronagraph_RRPM_optimize(void)
{
    double best_rad = g_coro.RRPM_RADIUS;
    double best_p2  = g_coro.RRPM_P2;
    double best_p3  = g_coro.RRPM_P3;
    double best_p4  = g_coro.RRPM_P4;
    double best_p5  = g_coro.RRPM_P5;

    EXECUTE_SYSTEM_COMMAND(
        "rm %s/rrpmref_fm.ref %s/rrpmref_pm.ref",
        CORONAGRAPHSDATALOCAL,
        CORONAGRAPHSDATALOCAL);
    delete_image_ID("rrpmref_fm.ref", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("rrpmref_pm.ref", DELETE_IMAGE_ERRMODE_WARNING);
    coronagraph_simul_RRPM(0.0, 0.0, "psf");
    double best_value = arith_image_total("psf");
    delete_image_ID("psf", DELETE_IMAGE_ERRMODE_WARNING);

    for (;;)
    {
        g_coro.RRPM_RADIUS += 0.01 * pow((0.5 - ran1()), 3.0);
        g_coro.RRPM_P2 += 0.0005 * pow((0.5 - ran1()), 3.0);
        g_coro.RRPM_P3 += 0.0001 * pow((0.5 - ran1()), 3.0);
        g_coro.RRPM_P4 += 0.0001 * pow((0.5 - ran1()), 3.0);
        g_coro.RRPM_P5 += 0.0001 * pow((0.5 - ran1()), 3.0);

        EXECUTE_SYSTEM_COMMAND(
            "rm %s/rrpmref_fm.ref %s/rrpmref_pm.ref",
            CORONAGRAPHSDATALOCAL,
            CORONAGRAPHSDATALOCAL);
        delete_image_ID("rrpmref_fm.ref", DELETE_IMAGE_ERRMODE_WARNING);
        delete_image_ID("rrpmref_pm.ref", DELETE_IMAGE_ERRMODE_WARNING);
        coronagraph_simul_RRPM(0.0, 0.0, "psf");
        double value = arith_image_total("psf");
        delete_image_ID("psf", DELETE_IMAGE_ERRMODE_WARNING);

        printf(
            "[%g %g %g %g %g] %g     [%g %g %g %g %g] %g\n",
            g_coro.RRPM_RADIUS,
            g_coro.RRPM_P2,
            g_coro.RRPM_P3,
            g_coro.RRPM_P4,
            g_coro.RRPM_P5,
            value,
            best_rad,
            best_p2,
            best_p3,
            best_p4,
            best_p5,
            best_value);

        if (value < best_value)
        {
            best_rad   = g_coro.RRPM_RADIUS;
            best_p2    = g_coro.RRPM_P2;
            best_p3    = g_coro.RRPM_P3;
            best_p4    = g_coro.RRPM_P4;
            best_p5    = g_coro.RRPM_P5;
            best_value = value;
        }
        else
        {
            g_coro.RRPM_RADIUS = best_rad;
            g_coro.RRPM_P2     = best_p2;
            g_coro.RRPM_P3     = best_p3;
            g_coro.RRPM_P4     = best_p4;
            g_coro.RRPM_P5     = best_p5;
        }
    }
}

/**
 * ovc_load_masks - Load or initialize OVC focal phase and amplitude reference maps
 * @out_fmp: Output vortex phase mask image ID
 * @out_fma: Output vortex amplitude mask image ID
 */
static void ovc_load_masks(
    imageID *out_fmp,
    imageID *out_fma)
{
    const long size   = CORONAGRAPHS_ARRAYSIZE;
    const long charge = g_coro.OVC_CHARGE;
    char       refname[200];

    snprintf(
        refname,
        sizeof(refname),
        "%s/ovcref_fmp%ld.ref.%ld",
        CORONAGRAPHSDATALOCAL,
        charge,
        size);
    imageID id_fmp = image_ID("ovcreffmp", dcimg, dcnimg);
    if (id_fmp == -1)
    {
        if (file_exists(refname) == 0)
        {
            coronagraph_init_OVC(charge);
        }
        load_fits(refname, "ovcreffmp", 1, &id_fmp);
    }
    *out_fmp = image_ID("ovcreffmp", dcimg, dcnimg);

    snprintf(
        refname,
        sizeof(refname),
        "%s/ovcref_fma%ld.ref.%ld",
        CORONAGRAPHSDATALOCAL,
        charge,
        size);
    imageID id_fma = image_ID("ovcreffma", dcimg, dcnimg);
    if (id_fma == -1)
    {
        if (file_exists(refname) == 0)
        {
            coronagraph_init_OVC(charge);
        }
        load_fits(refname, "ovcreffma", 1, &id_fma);
    }
    *out_fma = image_ID("ovcreffma", dcimg, dcnimg);
}

/**
 * ovc_propagate_to_lyot - Propagate entrance pupil through OVC vortex mask to Lyot plane "pc3"
 * @x:         Source X offset [lambda/D]
 * @y:         Source Y offset [lambda/D]
 * @id_fma:    Vortex amplitude mask image ID
 * @id_fmp:    Vortex phase mask image ID
 * @out_total: Output unapodized pupil flux integral
 *
 * Return: Image identifier of pre-stop Lyot complex field "pc3".
 */
static imageID ovc_propagate_to_lyot(
    double  x,
    double  y,
    imageID id_fma,
    imageID id_fmp,
    double *out_total)
{
    const long   size     = CORONAGRAPHS_ARRAYSIZE;
    const long   size2    = size * size;
    const double trad_pix = coro_pupil_radius_pix();

    imageID id_a1 = make_subpixdisk("pa1", size, size, 0.5 * size, 0.5 * size, trad_pix);
    *out_total    = arith_image_total("pa1");

    imageID id_amp = image_ID("coramp", dcimg, dcnimg);
    if (id_amp != -1)
    {
        printf("---------- USING AMPLITUDE FILE ----------\n");
        for (long ii = 0; ii < size2; ii++)
        {
            dcimg[id_a1].array.F[ii] = dcimg[id_amp].array.F[ii];
        }
    }

    imageID id_p1 =
        make_slopexy("pp1", size, size, CORO_PI * x / trad_pix, CORO_PI * y / trad_pix);
    imageID id_pha = image_ID("corpha", dcimg, dcnimg);
    if (id_pha != -1)
    {
        printf("---------- USING PHASE FILE ----------\n");
        for (long ii = 0; ii < size2; ii++)
        {
            dcimg[id_p1].array.F[ii] += dcimg[id_pha].array.F[ii];
        }
    }

    mk_complex_from_amph("pa1", "pp1", "pc2", 0);
    save_fl_fits("pa1", "pa1");
    save_fl_fits("pp1", "pp1");

    coro_fft2d_centered("pc2", "fc2");
    imageID id_fc2 = image_ID("fc2", dcimg, dcnimg);

    mk_amph_from_complex("fc2", "fa2", "fp2", 0);
    save_fl_fits("fa2", "fa2");
    save_fl_fits("fp2", "fp2");

    for (long ii = 0; ii < size2; ii++)
    {
        double re  = dcimg[id_fc2].array.CF[ii].re;
        double im  = dcimg[id_fc2].array.CF[ii].im;
        double amp = sqrt(re * re + im * im) * dcimg[id_fma].array.F[ii];
        double pha = atan2(im, re) + dcimg[id_fmp].array.F[ii];
        dcimg[id_fc2].array.CF[ii].re = (float) (amp * cos(pha));
        dcimg[id_fc2].array.CF[ii].im = (float) (amp * sin(pha));
    }

    coro_fft2d_centered("fc2", "pc3");
    mk_amph_from_complex("pc3", "pc3a", "pc3p", 0);
    save_fl_fits("pc3a", "pc3a");
    save_fl_fits("pc3p", "pc3p");
    delete_image_ID("pc3a", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pc3p", DELETE_IMAGE_ERRMODE_WARNING);

    return image_ID("pc3", dcimg, dcnimg);
}

/**
 * coronagraph_simul_OVC - Simulate Optical Vortex Coronagraph PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output PSF image name
 *
 * Return: RETURN_SUCCESS on completion.
 */
errno_t coronagraph_simul_OVC(
    double      xld,
    double      yld,
    const char *psfname)
{
    const long   size     = CORONAGRAPHS_ARRAYSIZE;
    const long   size2    = size * size;
    const double trad_pix = coro_pupil_radius_pix();
    imageID      id_fmp;
    imageID      id_fma;
    char         file_re[200];
    char         file_im[200];
    char         label[64];

    ovc_load_masks(&id_fmp, &id_fma);

    snprintf(
        file_re,
        sizeof(file_re),
        "%s/ovcre%ld.ref.%ld",
        CORONAGRAPHSDATALOCAL,
        g_coro.OVC_CHARGE,
        size);
    snprintf(
        file_im,
        sizeof(file_im),
        "%s/ovcim%ld.ref.%ld",
        CORONAGRAPHSDATALOCAL,
        g_coro.OVC_CHARGE,
        size);
    snprintf(label, sizeof(label), "OVC charge %ld", g_coro.OVC_CHARGE);

    coro_ref_pair_t ref = {
        .im_re   = "refovcre",
        .file_re = file_re,
        .im_im   = "refovcim",
        .file_im = file_im,
        .label   = label,
    };
    coro_load_ref_pair(&ref);

    double x = xld;
    double y = yld;
    if (ref.needs_init == 1)
    {
        printf("initialize OVC coronagraph charge %ld\n", g_coro.OVC_CHARGE);
        x = 0.0;
        y = 0.0;
    }

    double  total  = 0.0;
    imageID id_pc3 = ovc_propagate_to_lyot(x, y, id_fma, id_fmp, &total);

    coro_apply_or_save_ref_pair(&ref, id_pc3);

    for (long ii = 0; ii < size; ii++)
    {
        for (long jj = 0; jj < size; jj++)
        {
            double dx = 1.0 * ii - size / 2;
            double dy = 1.0 * jj - size / 2;
            if (sqrt(dx * dx + dy * dy) / trad_pix > 0.96)
            {
                dcimg[id_pc3].array.CF[jj * size + ii].re = 0.0f;
                dcimg[id_pc3].array.CF[jj * size + ii].im = 0.0f;
            }
        }
    }

    double s2 = (double) size2;
    coro_complex_to_psf("pc3", "fc4", "fa4", "fp4", psfname, total * s2 * s2 * s2);

    delete_image_ID("pc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pc3", DELETE_IMAGE_ERRMODE_WARNING);
    return RETURN_SUCCESS;
}
