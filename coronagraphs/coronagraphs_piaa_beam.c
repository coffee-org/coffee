/**
 * @file    coronagraphs_piaa_beam.c
 * @brief   PIAA state cleanup, Subaru pupil geometry, and Zernike PIAA sag perturbations
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "ZernikePolyn/ZernikePolyn.h"
#include "ZernikePolyn/zernike_value.h"

#include "coronagraphs_internal.h"

typedef struct piaa_pert_ctx
{
    const double *zarray;
    const long   *zindex;
    long          nb_zern;
    double        ratio;
    double        eps;
    double        epst;
} piaa_pert_ctx_t;

typedef struct piaa_vertex_res
{
    double x;
    double y;
    double v0;
    double dx;
    double dy;
} piaa_vertex_res_t;

/**
 * coronagraphs_free_PIAA - Free analytical apodizer fit coefficient buffers
 *
 * Return: 0 on completion.
 */
int coronagraphs_free_PIAA(void)
{
    free(g_coro.fitapo_a);
    free(g_coro.fitapo_b);
    free(g_coro.fitapo_c);
    free(g_coro.fitapo_c1);
    return 0;
}

/**
 * coronagraph_telescope_pupil_Subaru_inside1 - Test whether normalized point lies in Subaru pupil
 * @x: Normalized pupil x coordinate (unit pupil radius)
 * @y: Normalized pupil y coordinate (unit pupil radius)
 *
 * Return: 1 if inside clear aperture, 0 if obstructed.
 */
int coronagraph_telescope_pupil_Subaru_inside1(
    double x,
    double y)
{
    const double centobs       = 0.285;
    const double spider_offset = -0.11;
    const double spider_width  = 0.0275;
    const long   nbarms        = 4;
    const double eps           = 1.0e-20;
    double       x1[4];
    double       y1[4];
    double       x2[4];
    double       y2[4];

    int    inside = 1;
    double r      = sqrt(x * x + y * y);
    if ((r > 1.0) || (r < centobs))
    {
        inside = 0;
    }

    double pa = 0.775;
    x1[0]     = -sin(pa) * spider_offset;
    x2[0]     = cos(pa);
    y1[0]     = cos(pa) * spider_offset;
    y2[0]     = sin(pa);

    x1[1] = x1[0];
    x2[1] = x2[0];
    y1[1] = -y1[0];
    y2[1] = -y2[0];

    x1[2] = -x1[0];
    x2[2] = -x2[0];
    y1[2] = y1[0];
    y2[2] = y2[0];

    x1[3] = -x1[0];
    x2[3] = -x2[0];
    y1[3] = -y1[0];
    y2[3] = -y2[0];

    pa = 0.0;
    for (long arm = 0; arm < nbarms; arm++)
    {
        double tmpx = x1[arm] * cos(pa) + y1[arm] * sin(pa);
        double tmpy = -x1[arm] * sin(pa) + y1[arm] * cos(pa);
        x1[arm]     = tmpx;
        y1[arm]     = tmpy;
        tmpx        = x2[arm] * cos(pa) + y2[arm] * sin(pa);
        tmpy        = -x2[arm] * sin(pa) + y2[arm] * cos(pa);
        x2[arm]     = tmpx;
        y2[arm]     = tmpy;
    }

    for (long arm = 0; arm < nbarms; arm++)
    {
        double x12  = x2[arm] - x1[arm];
        double y12  = y2[arm] - y1[arm];
        double asq  = x12 * x12 + y12 * y12;
        double x13  = x - x1[arm];
        double y13  = y - y1[arm];
        double bsq  = x13 * x13 + y13 * y13;
        double dot  = x12 * x13 + y12 * y13;
        double absq = dot * dot;
        double dist = (bsq - absq / asq < eps) ? 0.0 : sqrt(bsq - absq / asq);

        if ((dot > 0.0) && (dist < spider_width))
        {
            inside = 0;
        }
    }

    return inside;
}

/**
 * coronagraphs_make_SUBARU_pupil - Create 4096x4096 Subaru pupil geometry image "subpup"
 *
 * Return: RETURN_SUCCESS on completion.
 */
errno_t coronagraphs_make_SUBARU_pupil(void)
{
    const long   ssize = 4096;
    const double rad   = 201.0;
    imageID      ID;

    g_coro.PIAACENTOBS = 0.278;
    create_2Dimage_ID("subpup", ssize, ssize, &ID);

    for (uint32_t ii = 0; ii < (uint32_t) ssize; ii++)
    {
        for (uint32_t jj = 0; jj < (uint32_t) ssize; jj++)
        {
            double x = (1.0 * ii - ssize / 2) / rad;
            double y = (1.0 * jj - ssize / 2) / rad;
            dcimg[ID].array.F[jj * ssize + ii] =
                (float) coronagraph_telescope_pupil_Subaru_inside1(x, y);
        }
    }

    return RETURN_SUCCESS;
}

/**
 * piaa_pert_eval_vertex - Compute Zernike sag and 1st-order ray deflection at one triangle vertex
 * @ctx:   Zernike perturbation parameters
 * @x0:    Center ray x coordinate [pupil radii]
 * @y0:    Center ray y coordinate [pupil radii]
 * @angle: Vertex azimuth offset [rad]
 * @res:   Output vertex coordinates, sag, and deflections
 */
static void piaa_pert_eval_vertex(
    const piaa_pert_ctx_t *ctx,
    double                 x0,
    double                 y0,
    double                 angle,
    piaa_vertex_res_t     *res)
{
    double px0 = x0 + ctx->epst * cos(angle);
    double py0 = y0 + ctx->epst * sin(angle);
    double pr0 = sqrt(px0 * px0 + py0 * py0);
    double pa0 = atan2(py0, px0);

    double px1 = px0 + ctx->eps;
    double py1 = py0;
    double pr1 = sqrt(px1 * px1 + py1 * py1);
    double pa1 = atan2(py1, px1);

    double px2 = px0;
    double py2 = py0 + ctx->eps;
    double pr2 = sqrt(px2 * px2 + py2 * py2);
    double pa2 = atan2(py2, px2);

    double v0 = 0.0;
    double v1 = 0.0;
    double v2 = 0.0;
    for (long k = 0; k < ctx->nb_zern; k++)
    {
        v0 += ctx->zarray[k] * Zernike_value(ctx->zindex[k], pr0, pa0);
        v1 += ctx->zarray[k] * Zernike_value(ctx->zindex[k], pr1, pa1);
        v2 += ctx->zarray[k] * Zernike_value(ctx->zindex[k], pr2, pa2);
    }

    res->x  = px0;
    res->y  = py0;
    res->v0 = v0;
    res->dx = -2.0 * ctx->ratio * (v1 - v0) / ctx->eps;
    res->dy = -2.0 * ctx->ratio * (v2 - v0) / ctx->eps;
}

/**
 * coronagraph_PIAAperturbation - Compute 1st-order PIAA remapping perturbation from Zernikes
 * @zarray: Zernike sag amplitude array [pupil radii]
 * @zindex: Noll/Zernike mode index array
 * @NBzern: Number of Zernike terms
 * @ratio:  Separation between PIAA optics [pupil radii]
 *
 * Return: 0 on completion.
 */
int coronagraph_PIAAperturbation(
    double *zarray,
    long   *zindex,
    long    NBzern,
    double  ratio)
{
    const uint32_t  size     = CORONAGRAPHS_ARRAYSIZE;
    const double    trad_pix = coro_pupil_radius_pix();
    piaa_pert_ctx_t ctx      = {
             .zarray  = zarray,
             .zindex  = zindex,
             .nb_zern = NBzern,
             .ratio   = ratio,
             .eps     = 1.0e-4,
             .epst    = 1.0e-3,
    };
    imageID IDzmap;
    imageID IDdx;
    imageID IDdy;
    imageID IDi;

    printf("RADIUS = %f\n", trad_pix);
    create_2Dimage_ID("PIAApert_z", size, size, &IDzmap);
    create_2Dimage_ID("PIAApert_dx", size, size, &IDdx);
    create_2Dimage_ID("PIAApert_dy", size, size, &IDdy);
    create_2Dimage_ID("PIAApert_i", size, size, &IDi);

    zernike_init();

    for (uint32_t ii = 0; ii < size; ii++)
    {
        for (uint32_t jj = 0; jj < size; jj++)
        {
            double x0 = (1.0 * ii - size / 2) / trad_pix;
            double y0 = (1.0 * jj - size / 2) / trad_pix;
            double r0 = sqrt(x0 * x0 + y0 * y0);
            if (r0 >= 1.1)
            {
                continue;
            }

            piaa_vertex_res_t p0;
            piaa_vertex_res_t p1;
            piaa_vertex_res_t p2;
            piaa_pert_eval_vertex(&ctx, x0, y0, 0.0, &p0);
            piaa_pert_eval_vertex(&ctx, x0, y0, 2.0 * CORO_PI / 3.0, &p1);
            piaa_pert_eval_vertex(&ctx, x0, y0, 4.0 * CORO_PI / 3.0, &p2);

            uint32_t idx                = jj * size + ii;
            dcimg[IDzmap].array.F[idx]  = (float) ((p0.v0 + p1.v0 + p2.v0) / 3.0);
            dcimg[IDdx].array.F[idx]    = (float) ((p0.dx + p1.dx + p2.dx) / 3.0);
            dcimg[IDdy].array.F[idx]    = (float) ((p0.dy + p1.dy + p2.dy) / 3.0);

            double p0xp = p0.x + p0.dx;
            double p0yp = p0.y + p0.dy;
            double p1xp = p1.x + p1.dx;
            double p1yp = p1.y + p1.dy;
            double p2xp = p2.x + p2.dx;
            double p2yp = p2.y + p2.dy;

            double tpre =
                0.5 * fabs((p0.x - p2.x) * (p1.y - p0.y) - (p0.x - p1.x) * (p2.y - p0.y));
            double tpost =
                0.5 * fabs((p0xp - p2xp) * (p1yp - p0yp) - (p0xp - p1xp) * (p2yp - p0yp));
            dcimg[IDi].array.F[idx] = (float) (tpre / tpost);
        }
    }

    save_fl_fits("PIAApert_z", "PIAApert_z.fits");
    save_fl_fits("PIAApert_dx", "PIAApert_dx.fits");
    save_fl_fits("PIAApert_dy", "PIAApert_dy.fits");
    save_fl_fits("PIAApert_i", "PIAApert_i.fits");

    return 0;
}
