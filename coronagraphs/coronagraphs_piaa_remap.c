/**
 * @file    coronagraphs_piaa_remap.c
 * @brief   PIAA ray-remapping lookup tables, mirror sag shapes, and 2D beam apodization
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "coronagraphs_internal.h"

typedef struct piaa_remap_ctx
{
    long    size;
    double  trad_pix;
    int     use_2dprol;
    imageID id_a1;
    imageID id_o1;
    imageID id_a2;
    imageID id_o2;
    imageID id_prol2d;
} piaa_remap_ctx_t;

/**
 * piaa_resample_lookup - Bin and linearly fill uniform ray-mapping lookup table
 * @x_src:      Source independent radial samples
 * @y_src:      Source dependent radial samples
 * @nb_points:  Number of source radial samples
 * @ext_factor: Radial extent normalization factor for x_src
 * @out_fname:  Optional ASCII file name to write lookup table
 * @out_table:  Output allocated lookup array (piaaconfNBr1fr0_pts elements)
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE on allocation failure.
 */
static errno_t piaa_resample_lookup(
    const double *x_src,
    const double *y_src,
    long          nb_points,
    double        ext_factor,
    const char   *out_fname,
    double      **out_table)
{
    const long npts     = piaaconfNBr1fr0_pts;
    double    *table    = (double *) calloc((size_t) npts, sizeof(double));
    long      *cntarray = (long *) calloc((size_t) npts, sizeof(long));

    if ((table == NULL) || (cntarray == NULL))
    {
        PRINT_ERROR("calloc returns NULL pointer");
        free(table);
        free(cntarray);
        return RETURN_FAILURE;
    }

    for (long ii = 0; ii < nb_points; ii++)
    {
        long jj = (long) ((x_src[ii] / ext_factor) * npts + 0.5);
        if ((jj > -1) && (jj < npts))
        {
            table[jj] += y_src[ii];
            cntarray[jj]++;
        }
    }

    for (long ii = 0; ii < npts; ii++)
    {
        table[ii] /= cntarray[ii];
    }

    for (long ii = 0; ii < npts; ii++)
    {
        if (cntarray[ii] == 0)
        {
            uint32_t ii0 = (uint32_t) ii;
            uint32_t ii1 = (uint32_t) ii;
            while (cntarray[ii0] == 0)
            {
                ii0--;
            }
            while (cntarray[ii1] == 0)
            {
                ii1++;
            }
            double x  = (1.0 * ii - ii0) / (ii1 - ii0);
            table[ii] = (1.0 - x) * table[ii0] + x * table[ii1];
        }
    }

    if (out_fname != NULL)
    {
        FILE *fpout = fopen(out_fname, "w");
        if (fpout != NULL)
        {
            for (long ii = 0; ii < npts; ii++)
            {
                fprintf(fpout, "%18.16f %18.16f\n", ext_factor * ii / npts, table[ii]);
            }
            fclose(fpout);
        }
    }

    free(cntarray);
    *out_table = table;
    return RETURN_SUCCESS;
}

/**
 * piaa_compute_mirror_shapes - Integrate radial PIAA M0/M1 mirror sag profiles
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE on allocation failure.
 */
static errno_t piaa_compute_mirror_shapes(void)
{
    const long   npts        = piaaconfNBr1fr0_pts;
    const double beam_radius = 0.01;
    const double mirror_dist = 1.0;

    printf("Compute PIAA mirror shapes\n");
    g_coro.piaaconfM0 = (double *) malloc(sizeof(double) * npts);
    g_coro.piaaconfM1 = (double *) malloc(sizeof(double) * npts);
    if ((g_coro.piaaconfM0 == NULL) || (g_coro.piaaconfM1 == NULL))
    {
        PRINT_ERROR("malloc returns NULL pointer");
        return RETURN_FAILURE;
    }

    g_coro.piaaconfM0[0] = 0.0;
    g_coro.piaaconfM1[0] = mirror_dist;

    for (long i = 0; i < npts - 1; i++)
    {
        double r0c   = g_coro.PIAAextfactor0 * i / npts;
        double r1c   = g_coro.piaaconfr1fr0[i];
        double dx    = (r0c - r1c) * beam_radius;
        double dz    = g_coro.piaaconfM1[i] - g_coro.piaaconfM0[i];
        double dist  = sqrt(dx * dx + dz * dz);
        double y3    = dist - dz;
        double slope = (fabs(dx) > 1.0e-9) ? (y3 / dx) : 0.0;
        double r0n   = g_coro.PIAAextfactor0 * (i + 1) / npts;

        g_coro.piaaconfM0[i + 1] = g_coro.piaaconfM0[i] + slope * (r0n - r0c) * beam_radius;
        g_coro.piaaconfM1[i + 1] =
            g_coro.piaaconfM1[i] + slope * (g_coro.piaaconfr1fr0[i + 1] - r1c) * beam_radius;
    }

    FILE *fp = fopen("PIAA_Mshapes.txt", "w");
    if (fp != NULL)
    {
        for (long ii = 0; ii < npts; ii++)
        {
            fprintf(
                fp,
                "%18.16f %18.16f %18.16f %18.16f\n",
                g_coro.PIAAextfactor0 * ii / npts * beam_radius,
                g_coro.piaaconfM0[ii],
                g_coro.piaaconfr1fr0[ii] * beam_radius,
                g_coro.piaaconfM1[ii]);
        }
        fclose(fp);
    }

    return RETURN_SUCCESS;
}

/**
 * coro_piaa_build_lookup_and_shapes - Build forward/inverse PIAA lookup tables and mirror sags
 * @nb_points: Number of adaptive radial samples in piaaconfr0/piaaconfr1
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE on error.
 */
errno_t coro_piaa_build_lookup_and_shapes(long nb_points)
{
    printf("compute r1 as a function of r0\n");
    fflush(stdout);
    if (piaa_resample_lookup(
            g_coro.piaaconfr0,
            g_coro.piaaconfr1,
            nb_points,
            g_coro.PIAAextfactor0,
            "PIAA_r0r1.txt",
            &g_coro.piaaconfr1fr0) != RETURN_SUCCESS)
    {
        return RETURN_FAILURE;
    }

    if (piaa_compute_mirror_shapes() != RETURN_SUCCESS)
    {
        return RETURN_FAILURE;
    }

    printf("compute r0 as a function of r1\n");
    fflush(stdout);
    return piaa_resample_lookup(
        g_coro.piaaconfr1,
        g_coro.piaaconfr0,
        nb_points,
        g_coro.PIAAextfactor1,
        "PIAA_r1r0.txt",
        &g_coro.piaaconfr0fr1);
}

/**
 * piaa_load_norm_prol2d - Load and RMS-normalize "apo2Dimt.fits" into image "prol2d"
 * @size2: Total pixel count in 2D image
 *
 * Return: Image identifier for "prol2d".
 */
static imageID piaa_load_norm_prol2d(long size2)
{
    imageID id_prol = image_ID("prol2d", dcimg, dcnimg);
    if (id_prol == -1)
    {
        load_fits("apo2Dimt.fits", "prol2d", 1, &id_prol);
    }

    double total_prol = 0.0;
    long   prol_cnt   = 0;
    for (long ii = 0; ii < size2; ii++)
    {
        double val = dcimg[id_prol].array.F[ii];
        if (val > 0.00001)
        {
            total_prol += val * val;
            prol_cnt++;
        }
    }

    double norm = sqrt(total_prol / prol_cnt);
    for (long ii = 0; ii < size2; ii++)
    {
        dcimg[id_prol].array.F[ii] /= (float) norm;
    }

    return id_prol;
}

/**
 * piaa_remap_pixel - Remap amplitude and bilinearly interpolated OPD at pixel (ii2, jj2)
 * @ctx: PIAA 2D remapping context
 * @ii2: Output pixel x index
 * @jj2: Output pixel y index
 */
static void piaa_remap_pixel(
    const piaa_remap_ctx_t *ctx,
    long                    ii2,
    long                    jj2)
{
    const double  epsilon = 1.0e-10;
    const long    npts    = piaaconfNBr1fr0_pts;
    const long    size    = ctx->size;
    const long    out_idx = jj2 * size + ii2;
    const double  x       = 1.0 * (ii2 - size / 2) / ctx->trad_pix;
    const double  y       = 1.0 * (jj2 - size / 2) / ctx->trad_pix;
    double        r2      = sqrt(x * x + y * y);
    const double  r       = r2;
    const double  pa      = atan2(y, x);

    if (r2 >= 1.0)
    {
        return;
    }

    double        ext   = (g_coro.piaaconfdirection == 0) ? g_coro.PIAAextfactor1
                                                          : g_coro.PIAAextfactor0;
    const double *table = (g_coro.piaaconfdirection == 0) ? g_coro.piaaconfr0fr1
                                                          : g_coro.piaaconfr1fr0;
    long          i;

    if (r2 < epsilon)
    {
        i  = 1;
        r2 = 1.0 * i / npts * ext;
    }
    else
    {
        i = (long) (r2 * npts / ext + 0.5);
    }
    if (i > npts - 2)
    {
        i = npts - 2;
    }

    double r1  = table[i];
    double fx  = 0.5 * size + r1 * ctx->trad_pix * cos(pa);
    double fy  = 0.5 * size + r1 * ctx->trad_pix * sin(pa);
    long   ii  = (long) fx;
    long   jj  = (long) fy;
    double t   = fx - ii;
    double u   = fy - jj;
    long   iin = (long) (fx + 0.5);
    long   jjn = (long) (fy + 0.5);

    if ((ii <= -1) || (jj <= -1) || (ii >= size - 1) || (jj >= size - 1))
    {
        dcimg[ctx->id_o2].array.F[out_idx] = 0.0f;
        return;
    }

    long ia = (long) (r * g_coro.PIAAAPO_NBPOINTS);
    if (ia >= g_coro.PIAAAPO_NBPOINTS)
    {
        return;
    }

    if (g_coro.piaaconfdirection == 0)
    {
        if (ctx->use_2dprol == 1)
        {
            long in_idx                        = jjn * size + iin;
            dcimg[ctx->id_a2].array.F[out_idx] = dcimg[ctx->id_prol2d].array.F[out_idx] *
                                                 dcimg[ctx->id_a1].array.F[in_idx];
        }
        else
        {
            dcimg[ctx->id_a2].array.F[out_idx] =
                (float) (dcimg[ctx->id_a1].array.F[jj * size + ii] * g_coro.PIAAAPO2D[out_idx] *
                         g_coro.PIAA_HYBRID_CPAAPO[ia]);
        }
    }
    else
    {
        dcimg[ctx->id_a2].array.F[out_idx] =
            (float) (dcimg[ctx->id_a1].array.F[jj * size + ii] /
                     (g_coro.PIAAAPO2D[jj * size + ii] + epsilon));
    }

    if (r > PIAACPPMASKRAD1)
    {
        dcimg[ctx->id_a2].array.F[out_idx] = 0.0f;
    }

    const float *o1 = dcimg[ctx->id_o1].array.F;
    double       interp =
        (1.0 - t) * (1.0 - u) * o1[jj * size + ii] + t * (1.0 - u) * o1[jj * size + ii + 1] +
        t * u * o1[(jj + 1) * size + ii + 1] + (1.0 - t) * u * o1[(jj + 1) * size + ii];
    dcimg[ctx->id_o2].array.F[out_idx] = (float) interp;
}

/**
 * coronagraphs_PIAA_apodize_beam - Remap pupil amplitude and OPD maps across PIAA optics
 * @ampl1: Input amplitude image name
 * @opd1:  Input OPD/phase image name
 * @ampl2: Output remapped amplitude image name
 * @opd2:  Output remapped OPD/phase image name
 *
 * Return: RETURN_SUCCESS on completion.
 */
errno_t coronagraphs_PIAA_apodize_beam(
    const char *ampl1,
    const char *opd1,
    const char *ampl2,
    const char *opd2)
{
    piaa_remap_ctx_t ctx = {
        .size       = CORONAGRAPHS_ARRAYSIZE,
        .trad_pix   = coro_pupil_radius_pix(),
        .use_2dprol = 1,
        .id_prol2d  = -1,
    };

    if (ctx.use_2dprol == 1)
    {
        ctx.id_prol2d = piaa_load_norm_prol2d(ctx.size * ctx.size);
    }

    ctx.id_a1 = image_ID(ampl1, dcimg, dcnimg);
    ctx.id_o1 = image_ID(opd1, dcimg, dcnimg);
    create_2Dimage_ID(ampl2, ctx.size, ctx.size, &ctx.id_a2);
    create_2Dimage_ID(opd2, ctx.size, ctx.size, &ctx.id_o2);

    for (long ii2 = 0; ii2 < ctx.size; ii2++)
    {
        for (long jj2 = 0; jj2 < ctx.size; jj2++)
        {
            piaa_remap_pixel(&ctx, ii2, jj2);
        }
    }

    return RETURN_SUCCESS;
}

/**
 * coro_piaa_inverse_unwrap - Run 3-phase-shifted inverse PIAA remap to avoid wrap artifacts
 * @in_amp:  Input post-Lyot amplitude image name
 * @in_pha:  Input post-Lyot phase image name
 * @out_amp: Output inverse-remapped amplitude image name
 * @out_pha: Output inverse-remapped unwrapped phase image name
 */
void coro_piaa_inverse_unwrap(
    const char *in_amp,
    const char *in_pha,
    const char *out_amp,
    const char *out_pha)
{
    const uint64_t size2 = (uint64_t) CORONAGRAPHS_ARRAYSIZE * CORONAGRAPHS_ARRAYSIZE;
    char           in_pha_b[64];
    char           in_pha_c[64];
    char           out_amp_b[64];
    char           out_pha_b[64];
    char           out_amp_c[64];
    char           out_pha_c[64];

    snprintf(in_pha_b, sizeof(in_pha_b), "%sb", in_pha);
    snprintf(in_pha_c, sizeof(in_pha_c), "%sc", in_pha);
    snprintf(out_amp_b, sizeof(out_amp_b), "%sb", out_amp);
    snprintf(out_pha_b, sizeof(out_pha_b), "%sb", out_pha);
    snprintf(out_amp_c, sizeof(out_amp_c), "%sc", out_amp);
    snprintf(out_pha_c, sizeof(out_pha_c), "%sc", out_pha);

    copy_image_ID(in_pha, in_pha_b, 0);
    copy_image_ID(in_pha, in_pha_c, 0);
    imageID id_in_b = image_ID(in_pha_b, dcimg, dcnimg);
    imageID id_in_c = image_ID(in_pha_c, dcimg, dcnimg);

    for (uint64_t ii = 0; ii < size2; ii++)
    {
        double vb = dcimg[id_in_b].array.F[ii] + 2.0 * CORO_PI / 3.0;
        if (vb > CORO_PI)
        {
            vb -= 2.0 * CORO_PI;
        }
        dcimg[id_in_b].array.F[ii] = (float) vb;

        double vc = dcimg[id_in_c].array.F[ii] - 2.0 * CORO_PI / 3.0;
        if (vc < -CORO_PI)
        {
            vc += 2.0 * CORO_PI;
        }
        dcimg[id_in_c].array.F[ii] = (float) vc;
    }

    g_coro.piaaconfdirection = 1;
    coronagraphs_PIAA_apodize_beam(in_amp, in_pha, out_amp, out_pha);
    coronagraphs_PIAA_apodize_beam(in_amp, in_pha_b, out_amp_b, out_pha_b);
    coronagraphs_PIAA_apodize_beam(in_amp, in_pha_c, out_amp_c, out_pha_c);
    g_coro.piaaconfdirection = 0;

    imageID id_out   = image_ID(out_pha, dcimg, dcnimg);
    imageID id_out_b = image_ID(out_pha_b, dcimg, dcnimg);
    imageID id_out_c = image_ID(out_pha_c, dcimg, dcnimg);

    for (uint64_t ii = 0; ii < size2; ii++)
    {
        double tmp2 = dcimg[id_out_b].array.F[ii];
        double tmp3 = dcimg[id_out_c].array.F[ii];
        if ((fabs(tmp2 - (2.0 * CORO_PI / 3.0 - CORO_PI)) < 0.3 * CORO_PI) &&
            (fabs(tmp3 + (2.0 * CORO_PI / 3.0 - CORO_PI)) < 0.3 * CORO_PI))
        {
            double tmp1 = CORO_PI + 0.5 * (tmp2 + tmp3);
            if (tmp1 > CORO_PI)
            {
                tmp1 -= 2.0 * CORO_PI;
            }
            dcimg[id_out].array.F[ii] = (float) tmp1;
        }
    }

    delete_image_ID(in_pha_b, DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID(in_pha_c, DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID(out_amp_b, DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID(out_pha_b, DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID(out_amp_c, DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID(out_pha_c, DELETE_IMAGE_ERRMODE_WARNING);
}
