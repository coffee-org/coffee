/**
 * @file    coronagraphs_sim_phase.c
 * @brief   Interferometric and phase-mask coronagraph simulators (SHEAR4, DICC, AIC, 4QPM)
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "coffee_compat.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "fft/fft.h"
#include "image_gen/image_gen.h"

#include "coronagraphs_internal.h"

/**
 * coronagraph_simul_SHEAR4 - Simulate 4-way lateral shearing interfero-coronagraph PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output PSF image name
 *
 * Return: 0 on completion.
 */
int coronagraph_simul_SHEAR4(
    double      xld,
    double      yld,
    const char *psfname)
{
    const long   size     = CORONAGRAPHS_ARRAYSIZE;
    const long   size2    = size * size;
    const double epsilon  = 1.0e-10;
    const double trad_pix = coro_pupil_radius_pix();
    const long   shearpix = (long) (g_coro.SHEAR4_SHEAR * trad_pix);
    double       total    = 0.0;
    imageID      IDc;

    coro_setup_pupil_tiptilt(xld, yld, 0, 1.0, "corphase", &total);
    imageID IDa1 = image_ID("pa1", dcimg, dcnimg);
    imageID IDp1 = image_ID("pp1", dcimg, dcnimg);
    create_2DCimage_ID("pc2", size, size, &IDc);

    long iistart = (long) (0.5 * size - trad_pix - shearpix);
    long iiend   = (long) (0.5 * size + trad_pix + shearpix);
    if (iistart < 1)
    {
        iistart = 1;
    }
    if (iiend > size - 1)
    {
        iiend = size - 1;
    }

    for (long jj = iistart; jj < iiend; jj++)
    {
        for (long ii = iistart; ii < iiend; ii++)
        {
            long idx1 = (jj + shearpix) * size + (ii - shearpix);
            long idx2 = (jj + shearpix) * size + (ii + shearpix);
            long idx3 = (jj - shearpix) * size + (ii - shearpix);
            long idx4 = (jj - shearpix) * size + (ii + shearpix);

            double a1 = dcimg[IDa1].array.F[idx1];
            double a2 = dcimg[IDa1].array.F[idx2];
            double a3 = dcimg[IDa1].array.F[idx3];
            double a4 = dcimg[IDa1].array.F[idx4];

            if (a1 * a2 * a3 * a4 > epsilon)
            {
                double p1 = dcimg[IDp1].array.F[idx1];
                double p2 = dcimg[IDp1].array.F[idx2];
                double p3 = dcimg[IDp1].array.F[idx3];
                double p4 = dcimg[IDp1].array.F[idx4];

                double re = a1 * cos(p1) - a2 * cos(p2) - a3 * cos(p3) + a4 * cos(p4);
                double im = a1 * sin(p1) - a2 * sin(p2) - a3 * sin(p3) + a4 * sin(p4);

                dcimg[IDc].array.CF[jj * size + ii].re = (float) (0.25 * re);
                dcimg[IDc].array.CF[jj * size + ii].im = (float) (0.25 * im);
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
 * dicc_init_modes - Create and orthonormalize 6 low-order pupil modes for DICC
 * @id_aref: Reference pupil disk image ID
 * @total:   Total pupil pixel area
 */
static void dicc_init_modes(
    imageID id_aref,
    double  total)
{
    const long size  = CORONAGRAPHS_ARRAYSIZE;
    const long size2 = size * size;
    imageID    id_1;
    imageID    id_x;
    imageID    id_y;
    imageID    id_x2;
    imageID    id_y2;
    imageID    id_xy;

    printf("Initializing modes...");
    fflush(stdout);
    create_2Dimage_ID("dicc1", size, size, &id_1);
    create_2Dimage_ID("diccx", size, size, &id_x);
    create_2Dimage_ID("diccy", size, size, &id_y);
    create_2Dimage_ID("diccx2", size, size, &id_x2);
    create_2Dimage_ID("diccy2", size, size, &id_y2);
    create_2Dimage_ID("diccxy", size, size, &id_xy);

    double totalx2s = 0.0;
    double totaly2s = 0.0;
    for (long ii = 0; ii < size; ii++)
    {
        for (long jj = 0; jj < size; jj++)
        {
            if (dcimg[id_aref].array.F[jj * size + ii] > 0.5f)
            {
                double x = 1.0 * (ii - size / 2);
                double y = 1.0 * (jj - size / 2);
                totalx2s += x * x;
                totaly2s += y * y;
            }
        }
    }
    totalx2s /= total;
    totaly2s /= total;

    double t1  = 0.0;
    double tx  = 0.0;
    double ty  = 0.0;
    double tx2 = 0.0;
    double ty2 = 0.0;
    double txy = 0.0;

    for (long ii = 0; ii < size; ii++)
    {
        for (long jj = 0; jj < size; jj++)
        {
            long idx = jj * size + ii;
            if (dcimg[id_aref].array.F[idx] > 0.5f)
            {
                double x  = 1.0 * (ii - size / 2);
                double y  = 1.0 * (jj - size / 2);
                double x2 = x * x - totalx2s;
                double y2 = y * y - totaly2s;

                dcimg[id_1].array.F[idx]  = 1.0f;
                dcimg[id_x].array.F[idx]  = (float) x;
                dcimg[id_y].array.F[idx]  = (float) y;
                dcimg[id_x2].array.F[idx] = (float) x2;
                dcimg[id_y2].array.F[idx] = (float) y2;
                dcimg[id_xy].array.F[idx] = (float) (x * y);

                t1 += 1.0;
                tx += x * x;
                ty += y * y;
                tx2 += x2 * x2;
                ty2 += y2 * y2;
                txy += x * y * x * y;
            }
        }
    }

    double tx2y2 = 0.0;
    for (long idx = 0; idx < size2; idx++)
    {
        if (dcimg[id_aref].array.F[idx] > 0.5f)
        {
            dcimg[id_1].array.F[idx] /= (float) sqrt(t1);
            dcimg[id_x].array.F[idx] /= (float) sqrt(tx);
            dcimg[id_y].array.F[idx] /= (float) sqrt(ty);
            dcimg[id_x2].array.F[idx] /= (float) sqrt(tx2);
            dcimg[id_y2].array.F[idx] /= (float) sqrt(ty2);
            dcimg[id_xy].array.F[idx] /= (float) sqrt(txy);
            tx2y2 += dcimg[id_x2].array.F[idx] * dcimg[id_y2].array.F[idx];
        }
    }

    ty2 = 0.0;
    for (long idx = 0; idx < size2; idx++)
    {
        if (dcimg[id_aref].array.F[idx] > 0.5f)
        {
            dcimg[id_y2].array.F[idx] -= (float) (tx2y2 * dcimg[id_x2].array.F[idx]);
            ty2 += dcimg[id_y2].array.F[idx] * dcimg[id_y2].array.F[idx];
        }
    }

    for (long idx = 0; idx < size2; idx++)
    {
        if (dcimg[id_aref].array.F[idx] > 0.5f)
        {
            dcimg[id_y2].array.F[idx] /= (float) sqrt(ty2);
        }
    }

    printf(" done\n");
    fflush(stdout);
}

/**
 * dicc_project_out_modes - Project out 6 low-order modes from pupil amplitude/phase
 * @id_aref: Reference pupil mask ID
 * @id_a1:   Pupil amplitude image ID
 * @id_p1:   Pupil phase image ID
 */
static void dicc_project_out_modes(
    imageID id_aref,
    imageID id_a1,
    imageID id_p1)
{
    const long   size2 = (long) CORONAGRAPHS_ARRAYSIZE * CORONAGRAPHS_ARRAYSIZE;
    const char  *names[6] = {"dicc1", "diccx", "diccy", "diccx2", "diccy2", "diccxy"};
    imageID      modes[6];

    for (int m = 0; m < 6; m++)
    {
        modes[m] = image_ID(names[m], dcimg, dcnimg);
    }

    for (long iter = 0; iter < 2; iter++)
    {
        double vr[6] = {0};
        double vi[6] = {0};

        for (long idx = 0; idx < size2; idx++)
        {
            if (dcimg[id_aref].array.F[idx] > 0.5f)
            {
                double amp = dcimg[id_a1].array.F[idx];
                double pha = dcimg[id_p1].array.F[idx];
                double re  = amp * cos(pha);
                double im  = amp * sin(pha);
                for (int m = 0; m < 6; m++)
                {
                    double mv = dcimg[modes[m]].array.F[idx];
                    vr[m] += re * mv;
                    vi[m] += im * mv;
                }
            }
        }

        for (long idx = 0; idx < size2; idx++)
        {
            if (dcimg[id_aref].array.F[idx] > 0.5f)
            {
                double amp = dcimg[id_a1].array.F[idx];
                double pha = dcimg[id_p1].array.F[idx];
                double re  = amp * cos(pha);
                double im  = amp * sin(pha);
                for (int m = 0; m < 6; m++)
                {
                    double mv = dcimg[modes[m]].array.F[idx];
                    re -= vr[m] * mv;
                    im -= vi[m] * mv;
                }
                dcimg[id_a1].array.F[idx] = (float) sqrt(re * re + im * im);
                dcimg[id_p1].array.F[idx] = (float) atan2(im, re);
            }
        }
    }
}

/**
 * coronagraph_simul_DICC - Simulate Differential Interferometric Coronagraph PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output PSF image name
 *
 * Return: 0 on completion.
 */
int coronagraph_simul_DICC(
    double      xld,
    double      yld,
    const char *psfname)
{
    const long   size     = CORONAGRAPHS_ARRAYSIZE;
    const long   size2    = size * size;
    const double trad_pix = coro_pupil_radius_pix();
    double       total    = 0.0;
    imageID      IDc;

    imageID id_aref = make_disk("paref", size, size, 0.5 * size, 0.5 * size, trad_pix);
    coro_setup_pupil_tiptilt(xld, yld, 0, 1.0, "corphase", &total);
    imageID id_a1 = image_ID("pa1", dcimg, dcnimg);
    imageID id_p1 = image_ID("pp1", dcimg, dcnimg);
    create_2DCimage_ID("pc2", size, size, &IDc);

    if (image_ID("dicc1", dcimg, dcnimg) == -1)
    {
        dicc_init_modes(id_aref, total);
    }

    dicc_project_out_modes(id_aref, id_a1, id_p1);
    delete_image_ID("paref", DELETE_IMAGE_ERRMODE_WARNING);

    for (long idx = 0; idx < size2; idx++)
    {
        double amp                    = dcimg[id_a1].array.F[idx];
        double pha                    = dcimg[id_p1].array.F[idx];
        dcimg[IDc].array.CF[idx].re = (float) (amp * cos(pha));
        dcimg[IDc].array.CF[idx].im = (float) (amp * sin(pha));
    }

    coro_complex_to_psf("pc2", "fc2", "fa2", "fp2", psfname, total * size2);
    delete_image_ID("pc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp1", DELETE_IMAGE_ERRMODE_WARNING);
    return 0;
}

/**
 * coronagraph_simul_AIC - Simulate Achromatic Interfero-Coronagraph PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output PSF image name
 *
 * Return: RETURN_SUCCESS on completion.
 */
errno_t coronagraph_simul_AIC(
    double      xld,
    double      yld,
    const char *psfname)
{
    const long   size     = CORONAGRAPHS_ARRAYSIZE;
    const long   size2    = size * size;
    const double epsilon  = 1.0e-10;
    const double trad_pix = coro_pupil_radius_pix();
    double       total    = 0.0;
    imageID      IDc;

    coro_setup_pupil_tiptilt(xld, yld, 0, 1.0, "corphase", &total);
    imageID IDa1 = image_ID("pa1", dcimg, dcnimg);
    imageID IDp1 = image_ID("pp1", dcimg, dcnimg);
    create_2DCimage_ID("pc2", size, size, &IDc);

    long iistart = (long) (0.5 * size - trad_pix - 2.0);
    long iiend   = (long) (0.5 * size + trad_pix + 2.0);
    if (iistart < 1)
    {
        iistart = 1;
    }
    if (iiend > size - 1)
    {
        iiend = size - 1;
    }

    for (long jj = iistart; jj < iiend; jj++)
    {
        long index1a = jj * size;
        long index2a = (size - jj) * size + size;
        for (long ii = iistart; ii < iiend; ii++)
        {
            long   index1 = index1a + ii;
            long   index2 = index2a - ii;
            double a11    = dcimg[IDa1].array.F[index1];
            if (a11 > epsilon)
            {
                double p11                     = dcimg[IDp1].array.F[index1];
                double a12                     = dcimg[IDa1].array.F[index2];
                double p12                     = dcimg[IDp1].array.F[index2];
                dcimg[IDc].array.CF[index1].re = (float) (0.5 * (a11 * cos(p11) - a12 * cos(p12)));
                dcimg[IDc].array.CF[index1].im = (float) (0.5 * (a11 * sin(p11) - a12 * sin(p12)));
            }
        }
    }

    coro_complex_to_psf("pc2", "fc2", "fa2", "fp2", psfname, total * size2);
    delete_image_ID("pc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp1", DELETE_IMAGE_ERRMODE_WARNING);
    return RETURN_SUCCESS;
}

/**
 * coronagraph_simul_4QPM - Simulate Four-Quadrant Phase Mask coronagraph PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output PSF image name
 *
 * Return: RETURN_SUCCESS on completion.
 */
errno_t coronagraph_simul_4QPM(
    double      xld,
    double      yld,
    const char *psfname)
{
    const long size  = CORONAGRAPHS_ARRAYSIZE;
    const long size2 = size * size;
    char       file_re[200];
    char       file_im[200];

    snprintf(file_re, sizeof(file_re), "%s/4qre.ref.%ld", CORONAGRAPHSDATALOCAL, size);
    snprintf(file_im, sizeof(file_im), "%s/4qim.ref.%ld", CORONAGRAPHSDATALOCAL, size);

    coro_ref_pair_t ref = {
        .im_re   = "ref4qre",
        .file_re = file_re,
        .im_im   = "ref4qim",
        .file_im = file_im,
        .label   = "4QPM",
    };
    coro_load_ref_pair(&ref);

    double x = xld;
    double y = yld;
    if (ref.needs_init == 1)
    {
        printf("initialize 4 quadrant coronagraph\n");
        x = 0.0;
        y = 0.0;
    }

    double total = 0.0;
    coro_setup_pupil_tiptilt(x, y, 1, 1.0, "corphase", &total);
    imageID IDa1 = image_ID("pa1", dcimg, dcnimg);

    mk_complex_from_amph("pa1", "pp1", "pc2", 0);
    coro_fft2d_centered("pc2", "fc2");

    imageID ID = image_ID("fc2", dcimg, dcnimg);
    for (long ii = 0; ii < size; ii++)
    {
        for (long jj = 0; jj < size; jj++)
        {
            if (((ii > size / 2) && (jj > size / 2)) || ((ii < size / 2) && (jj < size / 2)))
            {
                dcimg[ID].array.CF[jj * size + ii].re *= -1.0f;
                dcimg[ID].array.CF[jj * size + ii].im *= -1.0f;
            }
        }
    }
    for (long ii = 0; ii < size; ii++)
    {
        dcimg[ID].array.CF[(size / 2) * size + ii].re = 0.0f;
        dcimg[ID].array.CF[(size / 2) * size + ii].im = 0.0f;
        dcimg[ID].array.CF[ii * size + (size / 2)].re = 0.0f;
        dcimg[ID].array.CF[ii * size + (size / 2)].im = 0.0f;
    }

    coro_fft2d_centered("fc2", "pc3");
    ID = image_ID("pc3", dcimg, dcnimg);
    for (long ii = 0; ii < size2; ii++)
    {
        if (dcimg[IDa1].array.F[ii] < 0.99999f)
        {
            dcimg[ID].array.CF[ii].re = 0.0f;
            dcimg[ID].array.CF[ii].im = 0.0f;
        }
    }

    coro_apply_or_save_ref_pair(&ref, ID);
    coro_complex_to_psf(
        "pc3",
        "fc4",
        "fa4",
        "fp4",
        psfname,
        total * (double) size2 * (double) size2 * (double) size2);

    delete_image_ID("pc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pc3", DELETE_IMAGE_ERRMODE_WARNING);
    return RETURN_SUCCESS;
}
