/**
 * @file    init_geomPIAA_rad.c
 * @brief   Solve geometric ray-mapping ODEs for 1D radial PIAA mirror sags
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_memory/COREMOD_memory.h"

#include "PIAACMCsimul.h"
#include "PIAAshape/init_geomPIAA_rad.h"

#include "init_geomPIAA_rad_internal.h"

/**
 * build_output_apodization_flux - Reconstruct cosine-fit output amplitude `pup1` and enclosed flux.
 * @ctx:           Radial PIAA profile context to populate (`pup1`, `flux1cumul`, `flux1_in/out`)
 * @IDapofit_name: Name of 1D radial cosine fit coefficient image
 */
static void build_output_apodization_flux(
    struct geom_piaa_rad_ctx *ctx,
    const char      *restrict IDapofit_name)
{
    double coeffa  = 3.0;
    double coeffa1 = 0.5;

    imageID IDcoeff = image_ID(IDapofit_name, dcimg, dcnimg);
    long    nbcoeff = dcimg[IDcoeff].md->size[0];
    printf("%ld coefficients\n", nbcoeff);

    double total     = 0.0;
    double flux1_in  = 0.0;
    double flux1_out = 0.0;

    for (long ii = 0; ii < ctx->NBradpts; ii++)
    {
        ctx->pup1[ii] = 0.0;
        double r      = 1.0 * ii / ctx->NBradpts * piaacmcopticaldesign.r1lim;
        double r1;

        if (r < 1.0)
        {
            r1 = r;
        }
        else
        {
            r1 = 1.0 + (r - 1) / pow((1.0 + pow(1.0 / coeffa1 * (r - 1), coeffa)), 1.0 / coeffa);
        }

        for (long k = 0; k < nbcoeff; k++)
        {
            ctx->pup1[ii] += dcimg[IDcoeff].array.F[k] * cos(r1 * k * M_PI / ApoFitCosFact);
        }

        if (r < piaacmcopticaldesign.centObs1)
        {
            flux1_in += ctx->pup1[ii] * ctx->pup1[ii] * r;
        }
        if (r > 1.0)
        {
            flux1_out += ctx->pup1[ii] * ctx->pup1[ii] * r;
        }

        total += ctx->pup1[ii] * ctx->pup1[ii] * r;
        ctx->flux1cumul[ii] = total;
    }

    double normcoeff = 1.0 / (total - flux1_in - flux1_out);
    flux1_in *= normcoeff;
    flux1_out *= normcoeff;
    for (long ii = 0; ii < ctx->NBradpts; ii++)
    {
        ctx->flux1cumul[ii] *= normcoeff;
    }

    printf("outer fluxes 1: %lf %lf\n", flux1_in, flux1_out);
    ctx->flux1_in  = flux1_in;
    ctx->flux1_out = flux1_out;
}

/**
 * init_input_flat_flux - Initialize uniform input beam `pup0` and annular flux normalization.
 * @ctx: Radial PIAA profile context (`pup0`, `flux0cumul`, `normcoeff0`)
 */
static void init_input_flat_flux(struct geom_piaa_rad_ctx *ctx)
{
    double total     = 0.0;
    double flux0_in  = 0.0;
    double flux0_out = 0.0;

    for (long ii = 0; ii < ctx->NBradpts; ii++)
    {
        double r      = 1.0 * ii / ctx->NBradpts * piaacmcopticaldesign.r0lim;
        ctx->pup0[ii] = 1.0;

        if (r < piaacmcopticaldesign.centObs0)
        {
            flux0_in += ctx->pup0[ii] * ctx->pup0[ii] * r;
        }
        if (r > 1.0)
        {
            flux0_out += ctx->pup0[ii] * ctx->pup0[ii] * r;
        }

        total += ctx->pup0[ii] * ctx->pup0[ii] * r;
        ctx->flux0cumul[ii] = total;
    }

    double normcoeff = 1.0 / (total - flux0_in - flux0_out);
    flux0_in *= normcoeff;
    flux0_out *= normcoeff;

    printf("outer fluxes 0: %lf (%lf)    %lf\n", flux0_in, ctx->flux1_in, flux0_out);
    ctx->normcoeff0 = normcoeff;
}

/**
 * solve_inner_pseudo_profile - Match input central-obstruction pseudo-flux to output `flux1_in`.
 * @ctx: Radial PIAA profile context (`pup0`, `flux0cumul`)
 */
static void solve_inner_pseudo_profile(struct geom_piaa_rad_ctx *ctx)
{
    double a       = 0.5;
    double b       = 0.5;
    double bstep   = 0.1;
    double verr    = 1.0;
    double eps1    = 1e-8;
    long   NBistep =
        piaacmcopticaldesign.centObs0 * ctx->NBradpts / piaacmcopticaldesign.r0lim;
    int    dir     = 1;
    double t0      = 0.0;

    while (fabs(verr) > 1.0e-9)
    {
        t0 = 0.0;
        for (long ii = 0; ii < NBistep; ii++)
        {
            double x = 1.0 * ii / NBistep;
            double r = 1.0 * ii / ctx->NBradpts * piaacmcopticaldesign.r0lim;

            if (x < eps1)
            {
                x = eps1;
            }
            if (x > 1.0 - eps1)
            {
                x = 1.0 - eps1;
            }

            ctx->pup0[ii] =
                b + (1.0 - b) * (0.5 + atan(-a / b / (x * x) + a / pow(x - 1.0, 2)) / M_PI);
            t0 += r * ctx->pup0[ii] * ctx->pup0[ii];
            ctx->flux0cumul[ii] = t0;
        }

        verr     = t0 * ctx->normcoeff0 - ctx->flux1_in;
        int odir = dir;
        if (verr > 0.0)
        {
            b /= (1.0 + bstep);
            dir = -1;
        }
        else
        {
            b *= (1.0 + bstep);
            dir = 1;
        }
        if (odir * dir < 0)
        {
            bstep *= 0.1;
        }
        printf(".");
        fflush(stdout);
    }
    printf("\n");
    printf("TOTAL = %f -> %g (%g %g)\n", b, t0 * ctx->normcoeff0, bstep, verr);
}

/**
 * solve_outer_pseudo_profile - Match input outer-extrapolation pseudo-flux to output `flux1_out`.
 * @ctx: Radial PIAA profile context (`pup0`, `flux0cumul`)
 */
static void solve_outer_pseudo_profile(struct geom_piaa_rad_ctx *ctx)
{
    double a        = 0.5;
    double b        = 0.5;
    double bstep    = 0.1;
    double verr     = 1.0;
    double eps1     = 1e-8;
    long   iioffset = (long)(1.0 * ctx->NBradpts / piaacmcopticaldesign.r0lim);
    long   NBistep  = ctx->NBradpts - iioffset;
    int    dir      = 1;
    double t0       = 0.0;

    while (fabs(verr) > 1.0e-9)
    {
        t0 = 0.0;
        for (long ii = 0; ii < NBistep; ii++)
        {
            double x = 1.0 - 1.0 * ii / NBistep;
            double r = 1.0 + 1.0 * ii / ctx->NBradpts * piaacmcopticaldesign.r0lim;

            if (x < eps1)
            {
                x = eps1;
            }
            if (x > 1.0 - eps1)
            {
                x = 1.0 - eps1;
            }

            long idx       = ii + iioffset;
            ctx->pup0[idx] =
                b + (1.0 - b) * (0.5 + atan(-a / b / (x * x) + a / pow(x - 1.0, 2)) / M_PI);
            t0 += r * ctx->pup0[idx] * ctx->pup0[idx];
            ctx->flux0cumul[idx] = t0;
        }

        verr     = t0 * ctx->normcoeff0 - ctx->flux1_out;
        int odir = dir;
        if (verr > 0.0)
        {
            b /= (1.0 + bstep);
            dir = -1;
        }
        else
        {
            b *= (1.0 + bstep);
            dir = 1;
        }
        if (odir * dir < 0)
        {
            bstep *= 0.1;
        }
        printf(".");
        fflush(stdout);
    }
    printf("\n");
    printf("TOTAL = %f -> %g (%g %g)\n", b, t0 * ctx->normcoeff0, bstep, verr);
}

/**
 * finalize_and_save_pup01_prof - Normalize input cumulative flux and export `pup01.prof`.
 * @ctx: Radial PIAA profile context (`pup0`, `pup1`, `flux0cumul`, `flux1cumul`)
 */
static void finalize_and_save_pup01_prof(struct geom_piaa_rad_ctx *ctx)
{
    double total     = 0.0;
    double flux0_in  = 0.0;
    double flux0_out = 0.0;

    for (long ii = 0; ii < ctx->NBradpts; ii++)
    {
        double r = 1.0 * ii / ctx->NBradpts * piaacmcopticaldesign.r0lim;
        if (r < piaacmcopticaldesign.centObs0)
        {
            flux0_in += ctx->pup0[ii] * ctx->pup0[ii] * r;
        }
        if (r > 1.0)
        {
            flux0_out += ctx->pup0[ii] * ctx->pup0[ii] * r;
        }

        total += ctx->pup0[ii] * ctx->pup0[ii] * r;
        ctx->flux0cumul[ii] = total * ctx->normcoeff0;
    }
    flux0_in *= ctx->normcoeff0;
    flux0_out *= ctx->normcoeff0;

    printf("outer fluxes 0: %lf (%lf)    %lf\n", flux0_in, ctx->flux1_in, flux0_out);

    char fname[STRINGMAXLEN_FULLFILENAME];
    WRITE_FULLFILENAME(fname, "%s/pup01.prof", piaacmcparams.piaacmcconfdir);
    FILE *fp = fopen(fname, "w");
    if (fp == NULL)
    {
        PRINT_ERROR("Cannot create file %s", fname);
        abort();
    }
    for (long ii = 0; ii < ctx->NBradpts; ii++)
    {
        double r0 = 1.0 * ii / ctx->NBradpts * piaacmcopticaldesign.r0lim;
        double r1 = 1.0 * ii / ctx->NBradpts * piaacmcopticaldesign.r1lim;
        fprintf(fp,
                "%f %f %g %g %g %g\n",
                r0,
                r1,
                ctx->pup0[ii],
                ctx->pup1[ii],
                ctx->flux0cumul[ii],
                ctx->flux1cumul[ii]);
    }
    fclose(fp);
}

/**
 * init_geomPIAA_rad - Compute lossless 1D radial PIAA mirror profiles from cosine apodization fit.
 * @IDapofit_name: Input 1D cosine-fit radial amplitude apodization image name
 *
 * Only works for circular PIAA geometries. Writes `pup01.prof` and `PIAA_Mshapes.txt`.
 *
 * Return: RETURN_SUCCESS on completion, or milk error code if output file creation fails.
 */
errno_t init_geomPIAA_rad(const char *restrict IDapofit_name)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %s", IDapofit_name);

    struct geom_piaa_rad_ctx ctx = {0};
    ctx.NBradpts                 = piaacmcopticaldesign.NBradpts;

    ctx.pup0       = (double *) malloc(sizeof(double) * ctx.NBradpts);
    ctx.pup1       = (double *) malloc(sizeof(double) * ctx.NBradpts);
    ctx.flux0cumul = (double *) malloc(sizeof(double) * ctx.NBradpts);
    ctx.flux1cumul = (double *) malloc(sizeof(double) * ctx.NBradpts);
    if ((ctx.pup0 == NULL) || (ctx.pup1 == NULL) || (ctx.flux0cumul == NULL) ||
        (ctx.flux1cumul == NULL))
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    build_output_apodization_flux(&ctx, IDapofit_name);
    init_input_flat_flux(&ctx);
    solve_inner_pseudo_profile(&ctx);
    solve_outer_pseudo_profile(&ctx);
    finalize_and_save_pup01_prof(&ctx);

    errno_t status = geom_piaa_rad_march_and_save(&ctx);

    free(ctx.flux0cumul);
    free(ctx.flux1cumul);
    free(ctx.pup0);
    free(ctx.pup1);

    if (status != RETURN_SUCCESS)
    {
        return status;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
