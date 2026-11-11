/**
 * @file    coronagraphs_piaa_init.c
 * @brief   Initialization of radial PIAA apodization profiles and adaptive ray coordinates
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "coronagraphs_internal.h"

typedef struct piaa_init_ctx
{
    double  normfactor;
    double  total_ir;
    double  total_or;
    double  flux_in;
    long    nb_istep;
    double *innerprof_cumul;
    long    nb_points;
} piaa_init_ctx_t;

/**
 * piaa_eval_2d_profile - Sample 2D analytical apodization profile and flux fractions
 * @ctx: PIAA initialization context
 */
static void piaa_eval_2d_profile(piaa_init_ctx_t *ctx)
{
    const uint32_t size     = CORONAGRAPHS_ARRAYSIZE;
    const long     factor   = 3;
    const double   trad_pix = coro_pupil_radius_pix();
    imageID        IDt;
    imageID        IDtt;

    printf("factor = %ld\n", factor);
    create_2Dimage_ID("apo2Dim", size, size, &IDt);
    create_2Dimage_ID("apo2Dimt", size, size, &IDtt);

    double total    = 0.0;
    double total1   = 0.0;
    double total_or = 0.0;
    double total_ir = 0.0;

    for (long ii = 0; ii < factor * (long) size; ii++)
    {
        long ii1 = (long) ((1.0 * ii - 0.5) / factor + 0.5);
        for (long jj = 0; jj < factor * (long) size; jj++)
        {
            long   jj1   = (long) ((1.0 * jj - 0.5) / factor + 0.5);
            double x     = (1.0 / factor) * (ii - factor * (long) size / 2);
            double y     = (1.0 / factor) * (jj - factor * (long) size / 2);
            double r     = sqrt(x * x + y * y) / trad_pix;
            double value = 0.0;

            if (r < g_coro.PIAAextfactor1)
            {
                coronagraph_apofit_eval_r(r, &value);
                if (r > 1.0)
                {
                    total_or += value * value;
                }
                else if (r < g_coro.APLC_CentOBS1)
                {
                    total_ir += value * value;
                }
                else
                {
                    total += value * value;
                    total1 += 1.0;
                }
            }

            if ((ii1 > 0) && (jj1 > 0) && (ii1 < (long) size) && (jj1 < (long) size))
            {
                long idx = jj1 * (long) size + ii1;
                g_coro.PIAAAPO2D[idx] += value / factor / factor;
                dcimg[IDt].array.F[idx] = (float) g_coro.PIAAAPO2D[idx];
                dcimg[IDtt].array.F[idx] =
                    ((r > g_coro.APLC_CentOBS1) && (r < 1.0)) ? (float) g_coro.PIAAAPO2D[idx]
                                                              : 0.0f;
            }
        }
    }

    ctx->normfactor = sqrt(total1 / total);
    printf("NORMALIZATION FACTOR = %f\n", ctx->normfactor);
    fflush(stdout);

    for (uint64_t ii = 0; ii < (uint64_t) size * size; ii++)
    {
        dcimg[IDt].array.F[ii] *= (float) ctx->normfactor;
    }

    printf("%g %g %g\n", total_ir / total, total / total, total_or / total);
    ctx->total_ir = total_ir / total;
    ctx->total_or = total_or / total;
    ctx->flux_in  = ctx->total_ir;

    printf("Light within central obstruction (0->%f) = %g\n", g_coro.APLC_CentOBS1, ctx->total_ir);
    printf("Light in outer ring (1->%f) = %g\n", g_coro.PIAAextfactor1, ctx->total_or);

    save_fl_fits("apo2Dim", "apo2Dim.fits");
    save_fl_fits("apo2Dimt", "apo2Dimt.fits");

    double co0_sq = g_coro.APLC_CentOBS0 * g_coro.APLC_CentOBS0;
    printf("GOAL = %g\n", ctx->total_ir / (co0_sq / (1.0 - co0_sq)));
}

/**
 * piaa_solve_inner_profile - Solve smooth central-obstruction cumulative flux extension
 * @ctx: PIAA initialization context
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE on allocation failure.
 */
static errno_t piaa_solve_inner_profile(piaa_init_ctx_t *ctx)
{
    double b      = 0.5;
    double bstep  = 0.1;
    double verr   = 1.0;
    double dir    = 1.0;
    double t0     = 0.0;
    double t0cnt  = 0.0;
    double co0_sq = g_coro.APLC_CentOBS0 * g_coro.APLC_CentOBS0;
    double target = ctx->total_ir / (co0_sq / (1.0 - co0_sq));

    ctx->nb_istep        = 1000000;
    ctx->innerprof_cumul = (double *) malloc(sizeof(double) * ctx->nb_istep);
    if (ctx->innerprof_cumul == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        return RETURN_FAILURE;
    }

    while (fabs(verr) > 1.0e-9)
    {
        t0    = 0.0;
        t0cnt = 0.0;
        for (long i = 0; i < ctx->nb_istep; i++)
        {
            const double a    = 0.5;
            const double eps1 = 1.0e-8;
            double       x    = 1.0 * i / ctx->nb_istep;
            if (x < eps1)
            {
                x = eps1;
            }
            if (x > 1.0 - eps1)
            {
                x = 1.0 - eps1;
            }
            double xm1 = x - 1.0;
            double val = b + (1.0 - b) * (0.5 + atan(-a / b / (x * x) + a / (xm1 * xm1)) / CORO_PI);
            t0 += x * val * val;
            ctx->innerprof_cumul[i] = t0;
            t0cnt += x;
        }

        verr        = t0 / t0cnt - target;
        double odir = dir;
        if (verr > 0.0)
        {
            b /= (1.0 + bstep);
            dir = -1.0;
        }
        else
        {
            b *= (1.0 + bstep);
            dir = 1.0;
        }
        if (odir * dir < 0.0)
        {
            bstep *= 0.1;
        }
        printf(".");
        fflush(stdout);
    }

    printf("\nTOTAL = %f -> %g (%g %g)\n", b, t0 / t0cnt, bstep, verr);
    double last_val = ctx->innerprof_cumul[ctx->nb_istep - 1];
    for (long i = 0; i < ctx->nb_istep; i++)
    {
        ctx->innerprof_cumul[i] *= ctx->flux_in / last_val;
    }

    return RETURN_SUCCESS;
}

/**
 * piaa_fill_pup0sum - Evaluate normalized cumulative input flux at radial array r0_arr
 * @ctx:     PIAA initialization context
 * @r0_arr:  Input radial coordinate array
 * @npts:    Number of radial samples
 * @pup0sum: Output normalized cumulative flux array
 */
static void piaa_fill_pup0sum(
    const piaa_init_ctx_t *ctx,
    const double          *r0_arr,
    long                   npts,
    double                *pup0sum)
{
    double co_sq = g_coro.PIAACENTOBS * g_coro.PIAACENTOBS;
    double norm  = 1.0 + ctx->total_ir + ctx->total_or;

    for (long i = 0; i < npts; i++)
    {
        if (r0_arr[i] < g_coro.APLC_CentOBS0)
        {
            long i1 = (long) ((r0_arr[i] / g_coro.APLC_CentOBS0) * ctx->nb_istep);
            if (i1 > ctx->nb_istep - 1)
            {
                i1 = ctx->nb_istep - 1;
            }
            pup0sum[i] = ctx->innerprof_cumul[i1];
        }
        else
        {
            pup0sum[i] = ctx->flux_in + (r0_arr[i] * r0_arr[i] - co_sq) / (1.0 - co_sq);
        }
        pup0sum[i] /= norm;
    }
}

/**
 * piaa_compute_uniform_r0_r1 - Compute initial uniform-r0 to r1 cumulative flux mapping
 * @ctx:        PIAA initialization context
 * @out_r0:     Output uniform r0 array (piaaconfNPUPFILESIZE elements)
 * @out_r1:     Output mapped r1 array (piaaconfNPUPFILESIZE elements)
 * @out_pupsum: Output working array for output cumulative flux
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE on allocation failure.
 */
static errno_t piaa_compute_uniform_r0_r1(
    const piaa_init_ctx_t *ctx,
    double               **out_r0,
    double               **out_r1,
    double               **out_pupsum)
{
    const long nb0     = piaaconfNPUPFILESIZE;
    double    *r0      = (double *) malloc(sizeof(double) * nb0);
    double    *r1      = (double *) malloc(sizeof(double) * nb0);
    double    *pupsum  = (double *) malloc(sizeof(double) * nb0);
    double    *pup0sum = (double *) malloc(sizeof(double) * nb0);

    if ((r0 == NULL) || (r1 == NULL) || (pupsum == NULL) || (pup0sum == NULL))
    {
        PRINT_ERROR("malloc returns NULL pointer");
        free(r0);
        free(r1);
        free(pupsum);
        free(pup0sum);
        return RETURN_FAILURE;
    }

    for (long i = 0; i < nb0; i++)
    {
        r0[i] = g_coro.PIAAextfactor0 * i / (nb0 - 1);
    }
    if (g_coro.APLC_CentOBS0 > 0.01)
    {
        g_coro.PIAACENTOBS = g_coro.APLC_CentOBS0;
    }

    piaa_fill_pup0sum(ctx, r0, nb0, pup0sum);

    double total = 0.0;
    for (long ii = 0; ii < nb0; ii++)
    {
        double amp = g_coro.piaaconfpup_amp_profile[ii];
        total += amp * amp * ii;
        pupsum[ii] = total;
    }

    long ii = 0;
    r0[0]   = 0.0;
    r1[0]   = 0.0;
    for (long i = 1; i < nb0; i++)
    {
        double f0 = pup0sum[i];
        while (((pupsum[ii] / pupsum[nb0 - 1]) < f0) && (ii < nb0))
        {
            ii++;
        }
        double f1       = pupsum[ii - 1] / pupsum[nb0 - 1];
        double f2       = pupsum[ii] / pupsum[nb0 - 1];
        double fluxdens = (fabs(f2 - f1) > 1.0e-10) ? ((f2 - f1) / (2.0 * ii - 1.0)) : 1.0e-10;
        double x        = sqrt((f0 - f1) / fluxdens + (1.0 * ii * ii - 2.0 * ii + 1.0)) + 1.0 - ii;
        r1[i]           = g_coro.PIAAextfactor1 * (1.0 * ii - 1.0 + x) / (nb0 - 1);
    }

    free(pup0sum);
    *out_r0     = r0;
    *out_r1     = r1;
    *out_pupsum = pupsum;
    return RETURN_SUCCESS;
}

/**
 * piaa_resample_adaptive_r0_r1 - Build adaptively sampled piaaconfr0 and piaaconfr1 arrays
 * @ctx:    PIAA initialization context
 * @r0:     Uniform r0 array (freed before return)
 * @r1:     Uniform r1 array (freed before return)
 * @pupsum: Working cumulative output flux array (freed before return)
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE on allocation failure.
 */
static errno_t piaa_resample_adaptive_r0_r1(
    piaa_init_ctx_t *ctx,
    double          *r0,
    double          *r1,
    double          *pupsum)
{
    const long   nb0     = piaaconfNPUPFILESIZE;
    const double epsilon = 1.0e-12;

    printf("Computing number of points [%ld]...", nb0);
    fflush(stdout);

    long   i   = 0;
    long   ii  = 0;
    double tmp = 0.0;
    while (i < nb0 - 1)
    {
        double dr0  = r0[i + 1] - r0[i];
        double dr1  = (r1[i + 1] - r1[i] < epsilon) ? epsilon : (r1[i + 1] - r1[i]);
        double ndr0 = (dr0 < dr0 * dr0 / dr1) ? dr0 : (dr0 * dr0 / dr1);
        tmp += ndr0;
        ii++;
        if (tmp > r0[i])
        {
            i++;
        }
    }

    long nb_points   = ii;
    ctx->nb_points   = nb_points;
    printf("\n%ld Points\n", nb_points);
    fflush(stdout);

    g_coro.piaaconfr0 = (double *) malloc(sizeof(double) * nb_points);
    g_coro.piaaconfr1 = (double *) malloc(sizeof(double) * nb_points);
    if ((g_coro.piaaconfr0 == NULL) || (g_coro.piaaconfr1 == NULL))
    {
        printf("COULD NOT ALLOCATE MEMORY FOR piaaconfr0/piaaconfr1\n");
        exit(0);
    }

    printf("Filling up nr0 ...");
    fflush(stdout);
    i                    = 0;
    ii                   = 0;
    g_coro.piaaconfr0[0] = 0.0;
    while ((i < nb0 - 1) && (ii < nb_points - 1))
    {
        double dr0  = r0[i + 1] - r0[i];
        double dr1  = r1[i + 1] - r1[i];
        double ndr0 = (dr0 < dr0 * dr0 / dr1) ? dr0 : (dr0 * dr0 / dr1);
        if (ndr0 < 0.0)
        {
            ndr0 = 0.0;
        }
        g_coro.piaaconfr0[ii + 1] = g_coro.piaaconfr0[ii] + ndr0;
        ii++;
        if (g_coro.piaaconfr0[ii] > r0[i])
        {
            i++;
        }
    }
    free(r0);
    free(r1);
    printf("----- %ld    %ld %ld ----\n", ii, nb_points, nb0);

    double *pup0sum = (double *) malloc(sizeof(double) * nb_points);
    if (pup0sum == NULL)
    {
        free(pupsum);
        return RETURN_FAILURE;
    }
    piaa_fill_pup0sum(ctx, g_coro.piaaconfr0, nb_points, pup0sum);
    free(ctx->innerprof_cumul);
    ctx->innerprof_cumul = NULL;

    double total = 0.0;
    for (ii = 0; ii < nb0; ii++)
    {
        double amp = g_coro.piaaconfpup_amp_profile[ii];
        total += 2.0 * CORO_PI * amp * amp * ii;
        pupsum[ii] = total;
    }

    ii                   = 0;
    g_coro.piaaconfr1[0] = 0.0;
    for (i = 1; i < nb_points; i++)
    {
        double f0 = pup0sum[i];
        while (((pupsum[ii] / pupsum[nb0 - 1]) < f0) && (ii < nb0))
        {
            ii++;
        }
        double f1       = (ii != 0) ? (pupsum[ii - 1] / pupsum[nb0 - 1]) : 0.0;
        double f2       = pupsum[ii] / pupsum[nb0 - 1];
        double fluxdens = (f2 - f1) / (2.0 * CORO_PI * ii);
        double psr0_sq  = 1.0 * ii * ii - 1.0 * ii;
        g_coro.piaaconfr1[i] =
            g_coro.PIAAextfactor1 * sqrt((f0 - f1) / fluxdens / CORO_PI + psr0_sq) / (nb0 - 1);
    }

    g_coro.piaaconfNBpoints = nb_points;
    free(pup0sum);
    free(pupsum);
    return RETURN_SUCCESS;
}

/**
 * coronagraph_init_PIAA - Initialize radial PIAA apodization and ray-remapping tables
 * @outNBpoints: Optional output receiving number of adaptive radial samples
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE on error.
 */
errno_t coronagraph_init_PIAA(long *outNBpoints)
{
    DEBUG_TRACE_FSTART();

    piaa_init_ctx_t ctx = {0};

    printf("PIAA init    %ld points\n", g_coro.PIAAAPO_NBPOINTS);
    printf("piaaconfNPUPFILESIZE = %ld\n", (long) piaaconfNPUPFILESIZE);
    printf("NBpoints0 = %ld\n", (long) piaaconfNPUPFILESIZE);
    printf("APLC_CentOBS0 = %f\n", g_coro.APLC_CentOBS0);
    printf("APLC_CentOBS1 = %f\n", g_coro.APLC_CentOBS1);
    fflush(stdout);

    g_coro.PIAAextfactor1     = 2.0;
    g_coro.PIAAAPO            = (double *) malloc(sizeof(double) * g_coro.PIAAAPO_NBPOINTS);
    g_coro.PIAA_HYBRID_CPAAPO = (double *) malloc(sizeof(double) * g_coro.PIAAAPO_NBPOINTS);
    if ((g_coro.PIAAAPO == NULL) || (g_coro.PIAA_HYBRID_CPAAPO == NULL))
    {
        PRINT_ERROR("malloc returns NULL pointer");
        return RETURN_FAILURE;
    }

    piaa_eval_2d_profile(&ctx);
    if (piaa_solve_inner_profile(&ctx) != RETURN_SUCCESS)
    {
        return RETURN_FAILURE;
    }

    g_coro.PIAAextfactor0 = sqrt((ctx.total_ir + 1.0 + ctx.total_or) / (1.0 + ctx.total_ir));
    printf("PIAAextfactor0 = %f  \n", g_coro.PIAAextfactor0);

    g_coro.piaaconfpup_amp_profile = (double *) malloc(sizeof(double) * piaaconfNPUPFILESIZE);
    if (g_coro.piaaconfpup_amp_profile == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        free(ctx.innerprof_cumul);
        return RETURN_FAILURE;
    }

    for (long ii = 0; ii < piaaconfNPUPFILESIZE; ii++)
    {
        double v = 0.0;
        double r = g_coro.PIAAextfactor1 * ii / piaaconfNPUPFILESIZE;
        coronagraph_apofit_eval_r(r, &v);
        g_coro.piaaconfpup_amp_profile[ii] = ctx.normfactor * v;
    }

    free(g_coro.piaaconfr0);
    free(g_coro.piaaconfr1);
    g_coro.piaaconfr0 = NULL;
    g_coro.piaaconfr1 = NULL;

    double *r0     = NULL;
    double *r1     = NULL;
    double *pupsum = NULL;
    if (piaa_compute_uniform_r0_r1(&ctx, &r0, &r1, &pupsum) != RETURN_SUCCESS)
    {
        free(ctx.innerprof_cumul);
        return RETURN_FAILURE;
    }

    if (piaa_resample_adaptive_r0_r1(&ctx, r0, r1, pupsum) != RETURN_SUCCESS)
    {
        return RETURN_FAILURE;
    }

    if (coro_piaa_build_lookup_and_shapes(ctx.nb_points) != RETURN_SUCCESS)
    {
        return RETURN_FAILURE;
    }

    printf("init PIAA is done\n");
    fflush(stdout);
    if (outNBpoints != NULL)
    {
        *outNBpoints = ctx.nb_points;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
