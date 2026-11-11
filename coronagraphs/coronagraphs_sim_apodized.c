/**
 * @file    coronagraphs_sim_apodized.c
 * @brief   Pupil-apodization and pupil-interferometry simulators (CPA, PPA, NOCORO, STRIPC, SIMXY)
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "fft/fft.h"

#include "coronagraphs_internal.h"

typedef struct quad_interf_ctx
{
    imageID id_a1;
    imageID id_p1;
    imageID id_c;
} quad_interf_ctx_t;

/**
 * coronagraph_simul_CPA - Simulate Conventional Pupil Apodization coronagraph PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output PSF image name
 *
 * Return: RETURN_SUCCESS on completion.
 */
errno_t coronagraph_simul_CPA(
    double      xld,
    double      yld,
    const char *psfname)
{
    const long size  = CORONAGRAPHS_ARRAYSIZE;
    const long size2 = size * size;
    char       refname[200];

    snprintf(refname, sizeof(refname), "%s/cpapup.ref.%ld", CORONAGRAPHSDATALOCAL, size);
    imageID id_ref = image_ID("cpapupref", dcimg, dcnimg);
    if (id_ref == -1)
    {
        if (file_exists(refname) == 1)
        {
            load_fits(refname, "cpapupref", 1, &id_ref);
        }
        else
        {
            coronagraph_init_CPA();
        }
    }
    id_ref = image_ID("cpapupref", dcimg, dcnimg);

    double total = 0.0;
    coro_setup_pupil_tiptilt(xld, yld, 1, 1.0, "corphase", &total);
    imageID id_a1 = image_ID("pa1", dcimg, dcnimg);
    imageID id_a2;
    create_2Dimage_ID("pa2", size, size, &id_a2);

    for (long ii = 0; ii < size2; ii++)
    {
        dcimg[id_a2].array.F[ii] = dcimg[id_a1].array.F[ii] * dcimg[id_ref].array.F[ii];
    }

    mk_complex_from_amph("pa2", "pp1", "pc2", 0);
    coro_fft2d_centered("pc2", "fc2");

    imageID id_fc2 = image_ID("fc2", dcimg, dcnimg);
    for (long ii = 0; ii < size; ii++)
    {
        for (long jj = 0; jj < size; jj++)
        {
            double dx = 1.0 * ii - size / 2;
            double dy = 1.0 * jj - size / 2;
            if (sqrt(dx * dx + dy * dy) * CORONAGRAPHS_PIXSCALE < CPAFPMASKRAD)
            {
                dcimg[id_fc2].array.CF[jj * size + ii].re = 0.0f;
                dcimg[id_fc2].array.CF[jj * size + ii].im = 0.0f;
            }
        }
    }

    coro_fft2d_centered("fc2", "pc3");
    imageID id_pc3 = image_ID("pc3", dcimg, dcnimg);
    for (long ii = 0; ii < size2; ii++)
    {
        if (dcimg[id_a1].array.F[ii] > CPAPPMASKRAD)
        {
            dcimg[id_pc3].array.CF[ii].re = 0.0f;
            dcimg[id_pc3].array.CF[ii].im = 0.0f;
        }
    }

    double s2 = (double) size2;
    coro_complex_to_psf("pc3", "fc4", "fa4", "fp4", psfname, total * s2 * s2 * s2);

    delete_image_ID("pc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pc3", DELETE_IMAGE_ERRMODE_WARNING);
    return RETURN_SUCCESS;
}

/**
 * ppa_init_ref_maps - Build Pupil Phase Apodization amplitude/phase reference maps
 * @trad_pix: Pupil radius in pixels
 * @out_refa: Output amplitude reference image ID ("pparefa")
 * @out_refp: Output phase reference image ID ("pparefp")
 */
static void ppa_init_ref_maps(
    double   trad_pix,
    imageID *out_refa,
    imageID *out_refp)
{
    const long   size      = CORONAGRAPHS_ARRAYSIZE;
    const long   size2     = size * size;
    const long   nb_sample = 10000;
    const double a         = 2.0;
    const double epsilon   = 0.001;
    double      *pha_array = (double *) malloc((size_t) nb_sample * sizeof(double));
    long        *idx_array = (long *) malloc((size_t) nb_sample * sizeof(long));
    imageID      id_refr;
    imageID      id_refi;
    imageID      id_cnt;

    for (long k = 0; k < nb_sample; k++)
    {
        double x     = 2.0 * k / nb_sample - 1.0;
        pha_array[k] = a * log(((1.0 + epsilon) + x) / ((1.0 + epsilon) - x));
        idx_array[k] = (long) (x * (trad_pix - 1.0) / sqrt(2.0) + size / 2);
    }

    create_2Dimage_ID("pparefr", size, size, &id_refr);
    create_2Dimage_ID("pparefi", size, size, &id_refi);
    create_2Dimage_ID("pparefa", size, size, out_refa);
    create_2Dimage_ID("pparefp", size, size, out_refp);
    create_2Dimage_ID("cntref", size, size, &id_cnt);

    for (long kx = 0; kx < nb_sample; kx++)
    {
        double phax = pha_array[kx];
        long   ii   = idx_array[kx];
        for (long ky = 0; ky < nb_sample; ky++)
        {
            double pha   = phax + pha_array[ky];
            long   index = idx_array[ky] * size + ii;
            dcimg[id_refr].array.F[index] += (float) cos(pha);
            dcimg[id_refi].array.F[index] += (float) sin(pha);
            dcimg[id_cnt].array.F[index] += 1.0f;
        }
    }

    for (long ii = 0; ii < size2; ii++)
    {
        if (dcimg[id_cnt].array.F[ii] > 0.5f)
        {
            double re                     = dcimg[id_refr].array.F[ii] / dcimg[id_cnt].array.F[ii];
            double im                     = dcimg[id_refi].array.F[ii] / dcimg[id_cnt].array.F[ii];
            dcimg[*out_refa].array.F[ii] = (float) sqrt(re * re + im * im);
            dcimg[*out_refp].array.F[ii] = (float) atan2(im, re);
        }
    }

    delete_image_ID("pparefr", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pparefi", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("cntref", DELETE_IMAGE_ERRMODE_WARNING);
    free(pha_array);
    free(idx_array);
}

/**
 * coronagraph_simul_PPA - Simulate Pupil Phase Apodization coronagraph PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output PSF image name
 *
 * Return: 0 on completion.
 */
int coronagraph_simul_PPA(
    double      xld,
    double      yld,
    const char *psfname)
{
    const long   size     = CORONAGRAPHS_ARRAYSIZE;
    const long   size2    = size * size;
    const double trad_pix = coro_pupil_radius_pix();
    double       total    = 0.0;

    coro_setup_pupil_tiptilt(xld, yld, 1, 1.0, "corphase", &total);
    imageID id_a1 = image_ID("pa1", dcimg, dcnimg);
    imageID id_p1 = image_ID("pp1", dcimg, dcnimg);

    imageID id_refa = image_ID("pparefa", dcimg, dcnimg);
    imageID id_refp = image_ID("pparefp", dcimg, dcnimg);
    if ((id_refa == -1) || (id_refp == -1))
    {
        ppa_init_ref_maps(trad_pix, &id_refa, &id_refp);
    }

    imageID id_a2;
    create_2Dimage_ID("pa2", size, size, &id_a2);
    for (long ii = 0; ii < size2; ii++)
    {
        dcimg[id_p1].array.F[ii] += dcimg[id_refp].array.F[ii];
        dcimg[id_a2].array.F[ii] = dcimg[id_a1].array.F[ii] * dcimg[id_refa].array.F[ii];
    }

    mk_complex_from_amph("pa2", "pp1", "pc2", 0);
    coro_complex_to_psf("pc2", "fc2", "fa2", "fp2", psfname, total * size2);

    delete_image_ID("pc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp1", DELETE_IMAGE_ERRMODE_WARNING);
    return 0;
}

/**
 * coronagraph_simul_NOCORO - Simulate unocculted telescope Airy PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output PSF image name
 *
 * Return: 0 on completion.
 */
int coronagraph_simul_NOCORO(
    double      xld,
    double      yld,
    const char *psfname)
{
    const uint64_t size2 = (uint64_t) CORONAGRAPHS_ARRAYSIZE * CORONAGRAPHS_ARRAYSIZE;
    double         total = 0.0;

    coro_setup_pupil_tiptilt(xld, yld, 1, 1.0, "corphase", &total);
    mk_complex_from_amph("pa1", "pp1", "pc2", 0);
    coro_complex_to_psf("pc2", "fc2", "fa2", "fp2", psfname, total * size2);

    delete_image_ID("pc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp1", DELETE_IMAGE_ERRMODE_WARNING);
    return 0;
}

/**
 * apod_combine_quad_interf - Combine 4 symmetric pupil points with nulling signs into "pc2"
 * @ctx:  Pupil amplitude/phase and complex output image IDs
 * @idx1: Quadrant 1 pixel index
 * @idx2: Quadrant 2 pixel index
 * @idx3: Quadrant 3 pixel index
 * @idx4: Quadrant 4 pixel index
 */
static void apod_combine_quad_interf(
    const quad_interf_ctx_t *ctx,
    long                     idx1,
    long                     idx2,
    long                     idx3,
    long                     idx4)
{
    double a1 = dcimg[ctx->id_a1].array.F[idx1];
    double p1 = dcimg[ctx->id_p1].array.F[idx1];
    double a2 = dcimg[ctx->id_a1].array.F[idx2];
    double p2 = dcimg[ctx->id_p1].array.F[idx2];
    double p3 = dcimg[ctx->id_p1].array.F[idx3];
    double p4 = dcimg[ctx->id_p1].array.F[idx4];

    double r1 = a1 * cos(p1);
    double i1 = a1 * sin(p1);
    double r2 = a2 * cos(p2);
    double i2 = a2 * sin(p2);
    double r3 = a1 * cos(p3);
    double i3 = a1 * sin(p3);
    double r4 = a2 * cos(p4);
    double i4 = a2 * sin(p4);

    float dr = (float) (0.25 * ((r1 - r4) - (r2 - r3)));
    float di = (float) (0.25 * ((i1 - i4) - (i2 - i3)));

    dcimg[ctx->id_c].array.CF[idx1].re = dr;
    dcimg[ctx->id_c].array.CF[idx1].im = di;
    dcimg[ctx->id_c].array.CF[idx2].re = -dr;
    dcimg[ctx->id_c].array.CF[idx2].im = -di;
    dcimg[ctx->id_c].array.CF[idx3].re = dr;
    dcimg[ctx->id_c].array.CF[idx3].im = di;
    dcimg[ctx->id_c].array.CF[idx4].re = -dr;
    dcimg[ctx->id_c].array.CF[idx4].im = -di;
}

/**
 * coronagraph_simul_STRIPC - Simulate strip-offset 4-beam interferometric coronagraph PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output PSF image name
 *
 * Return: 0 on completion.
 */
int coronagraph_simul_STRIPC(
    double      xld,
    double      yld,
    const char *psfname)
{
    const uint32_t size     = CORONAGRAPHS_ARRAYSIZE;
    const uint64_t size2    = (uint64_t) size * size;
    const double   trad_pix = coro_pupil_radius_pix();
    const long     offset   = (long) (trad_pix * g_coro.STRIPCOFFSET);
    double         total    = 0.0;
    imageID        id_c;

    coro_setup_pupil_tiptilt(xld, yld, 0, 1.0, "corphase", &total);
    create_2DCimage_ID("pc2", size, size, &id_c);
    quad_interf_ctx_t ctx = {
        .id_a1 = image_ID("pa1", dcimg, dcnimg),
        .id_p1 = image_ID("pp1", dcimg, dcnimg),
        .id_c  = id_c,
    };

    for (long ii = 0; ii < (long) trad_pix; ii++)
    {
        for (long jj = 0; jj < (long) trad_pix; jj++)
        {
            double dx   = 2.0 * offset - ii;
            double dy   = 2.0 * offset - jj;
            double rad1 = sqrt(1.0 * ii * ii + 1.0 * jj * jj);
            double rad2 = sqrt(dx * dx + dy * dy);
            double rad3 = sqrt(dx * dx + 1.0 * jj * jj);
            double rad4 = sqrt(1.0 * ii * ii + dy * dy);

            if ((rad1 < trad_pix) && (rad2 < trad_pix) && (rad3 < trad_pix) && (rad4 < trad_pix) &&
                (ii < 2.0 * offset) && (jj < 2.0 * offset))
            {
                long idx1 = (size / 2 + jj) * size + (size / 2 + ii);
                long idx2 = (size / 2 - 2 * offset + jj) * size + (size / 2 + ii);
                long idx3 = (size / 2 - 2 * offset + jj) * size + (size / 2 - 2 * offset + ii);
                long idx4 = (size / 2 + jj) * size + (size / 2 - 2 * offset + ii);
                apod_combine_quad_interf(&ctx, idx1, idx2, idx3, idx4);
            }
        }
    }

    coro_complex_to_psf("pc2", "fc2", "fa2", "fp2", psfname, total * size2);
    delete_image_ID("pc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp1", DELETE_IMAGE_ERRMODE_WARNING);
    return 0;
}

/**
 * coronagraph_simul_SIMXY - Simulate symmetric XY-folding interferometric coronagraph PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output PSF image name
 *
 * Return: 0 on completion.
 */
int coronagraph_simul_SIMXY(
    double      xld,
    double      yld,
    const char *psfname)
{
    const uint32_t size     = CORONAGRAPHS_ARRAYSIZE;
    const uint64_t size2    = (uint64_t) size * size;
    const double   trad_pix = coro_pupil_radius_pix();
    double         total    = 0.0;
    imageID        id_c;

    coro_setup_pupil_tiptilt(xld, yld, 0, 1.0, "corphase", &total);
    create_2DCimage_ID("pc2", size, size, &id_c);
    quad_interf_ctx_t ctx = {
        .id_a1 = image_ID("pa1", dcimg, dcnimg),
        .id_p1 = image_ID("pp1", dcimg, dcnimg),
        .id_c  = id_c,
    };

    for (long ii = 0; ii < (long) trad_pix; ii++)
    {
        for (long jj = 0; jj < (long) trad_pix; jj++)
        {
            if (sqrt(1.0 * ii * ii + 1.0 * jj * jj) < trad_pix)
            {
                long idx1 = (size / 2 + jj) * size + (size / 2 + ii);
                long idx2 = (size / 2 - jj) * size + (size / 2 + ii);
                long idx3 = (size / 2 - jj) * size + (size / 2 - ii);
                long idx4 = (size / 2 + jj) * size + (size / 2 - ii);
                apod_combine_quad_interf(&ctx, idx1, idx2, idx3, idx4);
            }
        }
    }

    coro_complex_to_psf("pc2", "fc2", "fa2", "fp2", psfname, total * size2);
    delete_image_ID("pc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp1", DELETE_IMAGE_ERRMODE_WARNING);
    return 0;
}
