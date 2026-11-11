/**
 * @file    coronagraphs_masks_init.c
 * @brief   Initialization of CPA, ODC, BL8, BL4, RRPM, and OVC coronagraphic masks
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "image_gen/image_gen.h"

#include "coronagraphs_internal.h"

static double CPAAPO[CPAAPO_NBPOINTS];
static double BL8MASK[BL8MASK_NBSTEP];

/**
 * coronagraph_init_CPA - Build oversampled 2D Conventional Pupil Apodizer reference map
 *
 * Return: 0 on success, RETURN_FAILURE on file read error.
 */
int coronagraph_init_CPA(void)
{
    const uint32_t size     = CORONAGRAPHS_ARRAYSIZE;
    const long     factor   = 9;
    const double   trad_pix = coro_pupil_radius_pix();
    char           fname[200];
    char           fname1[200];
    imageID        ID;

    snprintf(fname, sizeof(fname), "%s/%s", CORONAGRAPHSDATALOCAL, CPAAPO_FNAME);
    FILE *fp = fopen(fname, "r");
    if (fp == NULL)
    {
        printf("ERROR : cannot open file \"%s\"\n", fname);
        exit(0);
    }

    for (long i = 0; i < CPAAPO_NBPOINTS; i++)
    {
        long   tmpl      = 0;
        double tmpf      = 0.0;
        int    fscanfcnt = fscanf(fp, "%ld %lf\n", &tmpl, &tmpf);
        if (coro_check_fscanf(fp, fscanfcnt, 2) != RETURN_SUCCESS)
        {
            fclose(fp);
            return RETURN_FAILURE;
        }
        CPAAPO[i] = tmpf;
    }
    fclose(fp);

    for (long i = 1; i < CPAAPO_NBPOINTS; i++)
    {
        CPAAPO[i] /= CPAAPO[0];
    }
    CPAAPO[0] = 1.0;

    create_2Dimage_ID("cpapupref", size, size, &ID);
    for (uint32_t ii = 0; ii < (uint32_t) factor * size; ii++)
    {
        for (uint32_t jj = 0; jj < (uint32_t) factor * size; jj++)
        {
            long   ii1   = (long) ((1.0 * ii - 0.5) / factor + 0.5);
            long   jj1   = (long) ((1.0 * jj - 0.5) / factor + 0.5);
            double x     = (1.0 / factor) * (1.0 * ii - 0.5 * factor * size);
            double y     = (1.0 / factor) * (1.0 * jj - 0.5 * factor * size);
            double r     = sqrt(x * x + y * y) / trad_pix;
            double ri    = r * CPAAPO_NBPOINTS;
            long   ri1   = (long) ri;
            double alpha = ri - ri1;
            double value = (ri < CPAAPO_NBPOINTS - 1)
                               ? ((1.0 - alpha) * CPAAPO[ri1] + alpha * CPAAPO[ri1 + 1])
                               : CPAAPO[CPAAPO_NBPOINTS - 1];

            if ((ii1 > 0) && (jj1 > 0) && (ii1 < (long) size) && (jj1 < (long) size))
            {
                dcimg[ID].array.F[jj1 * size + ii1] += (float) (value / factor / factor);
            }
        }
    }

    snprintf(fname1, sizeof(fname1), "%s/cpapup.ref.%ul", CORONAGRAPHSDATALOCAL, size);
    save_fl_fits("cpapupref", fname1);
    return 0;
}

/**
 * coronagraph_init_ODC - Build oversampled Optical Differentiation Coronagraph focal mask
 *
 * Return: 0 on completion.
 */
int coronagraph_init_ODC(void)
{
    const long size     = CORONAGRAPHS_ARRAYSIZE;
    const int  factor   = 4;
    const long size_odc = (size / 2) * factor;
    double    *odc1m    = (double *) malloc(sizeof(double) * size_odc);
    double     max      = 0.0;
    imageID    ID;
    char       fname1[200];

    for (long ii = 0; ii < size_odc; ii++)
    {
        double x  = CORONAGRAPHS_PIXSCALE * ii;
        odc1m[ii] = x * exp(-x * x / g_coro.ODC_GAUSS / g_coro.ODC_GAUSS);
        if (odc1m[ii] > max)
        {
            max = odc1m[ii];
            printf("x = %f\n", x);
        }
    }

    for (long ii = 0; ii < size_odc; ii++)
    {
        odc1m[ii] /= max;
    }

    create_2Dimage_ID("odcmref", size, size, &ID);
    for (long ii = 0; ii < factor * size; ii++)
    {
        long   ii1   = (long) ((1.0 * ii - 0.5) / factor + 0.5);
        long   index = (ii > factor * size / 2) ? (ii - factor * size / 2)
                                                : (factor * size / 2 - ii);
        double tmp   = (ii > factor * size / 2) ? odc1m[index] : -odc1m[index];

        if ((ii1 > 0) && (ii1 < size))
        {
            for (long jj = 0; jj < size; jj++)
            {
                dcimg[ID].array.F[jj * size + ii1] += (float) (tmp / factor);
            }
        }
    }

    free(odc1m);
    snprintf(fname1, sizeof(fname1), "%s/odcm.ref.%ld", CORONAGRAPHSDATALOCAL, size);
    save_fl_fits("odcmref", fname1);
    return 0;
}

/**
 * bl8_sample_oversampled_mask - Rasterize circular or linear 8th-order band-limited mask
 * @im_name:     Output shared-memory image name
 * @file_prefix: Output file basename prefix under CORONAGRAPHSDATALOCAL
 * @circular:    1 for radial mask, 0 for 1D linear mask
 */
static void bl8_sample_oversampled_mask(
    const char *im_name,
    const char *file_prefix,
    int         circular)
{
    const long size   = CORONAGRAPHS_ARRAYSIZE;
    const long factor = 9;
    imageID    ID;
    char       fname1[200];

    create_2Dimage_ID(im_name, size, size, &ID);
    for (long ii = 0; ii < factor * size; ii++)
    {
        for (long jj = 0; jj < factor * size; jj++)
        {
            long   ii1 = (long) ((1.0 * ii - 0.5) / factor + 0.5);
            long   jj1 = (long) ((1.0 * jj - 0.5) / factor + 0.5);
            double x   = (1.0 / factor) * (ii - factor * size / 2);
            double y   = (1.0 / factor) * (jj - factor * size / 2);
            double r   = circular ? sqrt(x * x + y * y) : fabs(x);
            long   idx = (long) ((r * CORONAGRAPHS_PIXSCALE) / BL8MASK_STEP);

            if ((ii1 > 0) && (jj1 > 0) && (ii1 < size) && (jj1 < size))
            {
                double val = (idx >= BL8MASK_NBSTEP) ? 1.0 : BL8MASK[idx];
                dcimg[ID].array.F[jj1 * size + ii1] += (float) (val / factor / factor);
            }
        }
    }

    snprintf(fname1, sizeof(fname1), "%s/%s.ref.%ld", CORONAGRAPHSDATALOCAL, file_prefix, size);
    save_fl_fits(im_name, fname1);
}

/**
 * coronagraph_init_BL8 - Build circular and linear 8th-order band-limited focal masks
 *
 * Return: 0 on completion.
 */
int coronagraph_init_BL8(void)
{
    double tmp1 = (g_coro.BL8MASK_l - g_coro.BL8MASK_m) / g_coro.BL8MASK_l;
    double max  = 0.0;

    BL8MASK[0] = 0.0;
    for (long i = 1; i < BL8MASK_NBSTEP; i++)
    {
        double x    = BL8MASK_STEP * i;
        double x2   = CORO_PI * x * g_coro.BL8MASK_eps / g_coro.BL8MASK_l;
        double tmp2 = pow(sin(x2) / x2, g_coro.BL8MASK_l);
        double x3   = CORO_PI * x * g_coro.BL8MASK_eps / g_coro.BL8MASK_m;
        double tmp3 = (g_coro.BL8MASK_m / g_coro.BL8MASK_l) * pow(sin(x3) / x3, g_coro.BL8MASK_m);
        BL8MASK[i]  = tmp1 - tmp2 + tmp3;
        if (BL8MASK[i] > max)
        {
            max = BL8MASK[i];
        }
    }

    for (long i = 0; i < BL8MASK_NBSTEP; i++)
    {
        BL8MASK[i] /= max;
    }

    bl8_sample_oversampled_mask("bl8mref", "bl8m", 1);
    bl8_sample_oversampled_mask("bl8mrefl", "bl8ml", 0);
    return 0;
}

/**
 * coronagraph_init_BL4 - Build oversampled 4th-order band-limited focal mask
 *
 * Return: 0 on completion.
 */
int coronagraph_init_BL4(void)
{
    const long size   = CORONAGRAPHS_ARRAYSIZE;
    const long factor = 9;
    imageID    ID;
    char       fname1[200];

    create_2Dimage_ID("bl4mref", size, size, &ID);
    for (long ii = 0; ii < factor * size; ii++)
    {
        for (long jj = 0; jj < factor * size; jj++)
        {
            long   ii1  = (long) ((1.0 * ii - 0.5) / factor + 0.5);
            long   jj1  = (long) ((1.0 * jj - 0.5) / factor + 0.5);
            double x    = (1.0 / factor) * (ii - factor * size / 2);
            double y    = (1.0 / factor) * (jj - factor * size / 2);
            double rx   = fabs(x) * CORONAGRAPHS_PIXSCALE;
            double r    = sqrt(x * x + y * y) * CORONAGRAPHS_PIXSCALE;
            double tmp1 = sin(rx * CORO_PI * g_coro.BL4MASK_eps / 2.0);
            double tmp2 = exp(-(r / 50.0) * (r / 50.0) / 2.0);

            if ((ii1 > 0) && (jj1 > 0) && (ii1 < size) && (jj1 < size))
            {
                dcimg[ID].array.F[jj1 * size + ii1] +=
                    (float) (tmp1 * tmp1 * tmp2 / factor / factor);
            }
        }
    }

    snprintf(fname1, sizeof(fname1), "%s/bl4m.ref.%ld", CORONAGRAPHSDATALOCAL, size);
    save_fl_fits("bl4mref", fname1);
    return 0;
}

/**
 * coronagraph_init_RRPM - Build Roddier & Roddier focal phase disk and polynomial pupil mask
 *
 * Return: 0 on completion.
 */
int coronagraph_init_RRPM(void)
{
    const long   size     = CORONAGRAPHS_ARRAYSIZE;
    const long   size2    = size * size;
    const double trad_pix = coro_pupil_radius_pix();
    imageID      ID;
    char         fname1[200];

    ID = make_subpixdisk(
        "rrpmreffm",
        size,
        size,
        size / 2,
        size / 2,
        g_coro.RRPM_RADIUS / CORONAGRAPHS_PIXSCALE);
    for (long ii = 0; ii < size2; ii++)
    {
        dcimg[ID].array.F[ii] = 1.0f - 2.0f * dcimg[ID].array.F[ii];
    }
    snprintf(fname1, sizeof(fname1), "%s/rrpmref_fm.ref.%ld", CORONAGRAPHSDATALOCAL, size);
    save_fl_fits("rrpmreffm", fname1);

    create_2Dimage_ID("rrpmrefpm", size, size, &ID);
    for (long ii = 0; ii < size; ii++)
    {
        for (long jj = 0; jj < size; jj++)
        {
            double x   = 1.0 * (ii - size / 2);
            double y   = 1.0 * (jj - size / 2);
            double r   = sqrt(x * x + y * y) / trad_pix;
            double val = 1.0 + g_coro.RRPM_P2 * r * r + g_coro.RRPM_P3 * pow(r, 3.0) +
                         g_coro.RRPM_P4 * pow(r, 4.0) + g_coro.RRPM_P5 * pow(r, 5.0);
            if (val > 1.0)
            {
                val = 1.0;
            }
            dcimg[ID].array.F[jj * size + ii] = (float) val;
        }
    }
    snprintf(fname1, sizeof(fname1), "%s/rrpmref_pm.ref.%ld", CORONAGRAPHSDATALOCAL, size);
    save_fl_fits("rrpmrefpm", fname1);
    return 0;
}

/**
 * coronagraph_init_OVC - Build subpixel-averaged Optical Vortex Coronagraph focal mask
 * @charge: Topological vortex charge
 *
 * Return: 0 on completion.
 */
int coronagraph_init_OVC(long charge)
{
    const long   size  = CORONAGRAPHS_ARRAYSIZE;
    const long   ssamp = 3;
    const double eps   = 1.0e-8;
    imageID      ID1;
    imageID      ID2;
    char         fname1[200];
    char         fnamea[200];
    char         fnamep[200];

    printf("Initializing OVC charge %ld...\n", charge);
    fflush(stdout);

    snprintf(fnamea, sizeof(fnamea), "ovcreffma%ld", charge);
    snprintf(fnamep, sizeof(fnamep), "ovcreffmp%ld", charge);
    create_2Dimage_ID(fnamep, size, size, &ID1);
    create_2Dimage_ID(fnamea, size, size, &ID2);

    for (long ii = 0; ii < size; ii++)
    {
        for (long jj = 0; jj < size; jj++)
        {
            double re  = 0.0;
            double im  = 0.0;
            long   cnt = 0;
            for (long i = -ssamp; i < ssamp + 1; i++)
            {
                for (long j = -ssamp; j < ssamp + 1; j++)
                {
                    double x     = 1.0 * (ii - size / 2) + 0.5 * i / (0.5 + ssamp);
                    double y     = 1.0 * (jj - size / 2) + 0.5 * j / (0.5 + ssamp);
                    double theta = 1.0 * charge * atan2(y, x);
                    if ((x * x + y * y) > eps)
                    {
                        re += cos(theta);
                        im += sin(theta);
                    }
                    cnt++;
                }
            }
            re /= cnt;
            im /= cnt;
            dcimg[ID2].array.F[jj * size + ii] = (float) sqrt(re * re + im * im);
            dcimg[ID1].array.F[jj * size + ii] = (float) atan2(im, re);
        }
    }

    snprintf(
        fname1, sizeof(fname1), "%s/ovcref_fmp%ld.ref.%ld", CORONAGRAPHSDATALOCAL, charge, size);
    save_fl_fits(fnamep, fname1);
    snprintf(
        fname1, sizeof(fname1), "%s/ovcref_fma%ld.ref.%ld", CORONAGRAPHSDATALOCAL, charge, size);
    save_fl_fits(fnamea, fname1);

    printf("done\n");
    fflush(stdout);
    return 0;
}
