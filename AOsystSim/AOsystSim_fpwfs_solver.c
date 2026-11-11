/**
 * @file    AOsystSim_fpwfs_solver.c
 * @brief   Per-pixel Nelder-Mead simplex electric-field solver for FPWFS sensitivity analysis
 */

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "statistic/statistic.h"

#include "AOsystSim_internal.h"

/**
 * @brief Solved per-pixel focal-plane electric-field and incoherent flux estimate.
 */
typedef struct
{
    double ptre_in;
    double ptim_in;
    double ptre;
    double ptim;
    double Iflux;
    double are;
    double aim;
    double e;
    double bestvalue;
    double optsize;
    long   iter;
} FPWFS_PIXEL_SOL;

/**
 * fpwfs_sample_pixel_probes - Populate ideal or image-simulated noisy probe intensities
 * @ctx:     Sensitivity analysis context
 * @kx:      Pixel X index
 * @ky:      Pixel Y index
 * @ptre:    True real electric-field amplitude
 * @ptim:    True imaginary electric-field amplitude
 * @is_test: Non-zero when pixel matches diagnostic test coordinates
 */
static void fpwfs_sample_pixel_probes(
    const AOSIM_FPWFS_SENS_CTX *ctx,
    long                        kx,
    long                        ky,
    double                      ptre,
    double                      ptim,
    int                         is_test)
{
    long idx2d = ky * ctx->mapxsize + kx;
    long plane = ctx->mapxsize * ctx->mapysize;

    for (int pr = 0; pr < ctx->NBprobes; pr++)
    {
        double re    = (double) probe_re[pr] - ptre;
        double im    = (double) probe_im[pr] - ptim;
        double mflux = ctx->Iflux + (double) Cflux * (re * re + im * im);

        double nre  = (double) nprobe_re[pr] - ptre;
        double nim  = (double) nprobe_im[pr] - ptim;
        double mean = ((ctx->Iflux + (double) Cflux * (nre * nre + nim * nim)) * ctx->FLUXph) /
                      ctx->NBprobes;
        double val  = fast_poisson(mean) + ctx->RON * gauss();

        probe_nmflux[pr]  = val / (ctx->FLUXph / ctx->NBprobes);
        probe_nmnoise[pr] = sqrt((val < 1.0) ? 1.0 : val);

        if ((ctx->mapmode == 0) && is_test)
        {
            printf(
                "M0    NO NOISE:  probe %3d  normalized flux = %8.5lf  (%g ph)\n",
                pr,
                mflux / (double) Cflux,
                mflux);
            printf(
                "M0  WITH NOISE:             normalized flux = %8.5lf  (%g ph)\n",
                (double) (probe_nmflux[pr] / Cflux),
                (double) probe_nmflux[pr]);
            fflush(stdout);
        }

        if (ctx->imSimulMode == 1)
        {
            float pamp        = dcimg[ctx->IDprobampC].array.F[idx2d];
            probe_nmflux[pr]  = dcimg[ctx->IDpsfC].array.F[pr * plane + idx2d] / pamp;
            probe_nmnoise[pr] = dcimg[ctx->IDpsfCnoise].array.F[pr * plane + idx2d] / pamp;
            if (is_test)
            {
                printf("PROBE %d   -> %g\n", pr, (double) probe_nmflux[pr]);
            }
        }
    }
}

/**
 * fpwfs_extract_solution - Unpack and canonicalize simplex parameter vector
 * @x:   GSL simplex state vector
 * @sol: Output parameter struct
 */
static void fpwfs_extract_solution(
    const gsl_vector *x,
    FPWFS_PIXEL_SOL  *sol)
{
    sol->ptre  = gsl_vector_get(x, 0);
    sol->ptim  = gsl_vector_get(x, 1);
    sol->Iflux = gsl_vector_get(x, 2);
    if (NBoptVar > 3)
    {
        sol->are = gsl_vector_get(x, 3);
    }
    if (NBoptVar > 4)
    {
        sol->aim = gsl_vector_get(x, 4);
        sol->e   = gsl_vector_get(x, 5);
        if (sol->e > 1.0)
        {
            double tmp = sol->are;
            sol->e     = 1.0 / sol->e;
            sol->are   = sol->aim;
            sol->aim   = -tmp;
        }
    }
}

/**
 * fpwfs_print_solution - Log converged FPWFS simplex solution for diagnostic pixels
 * @kx:      Pixel X index
 * @ky:      Pixel Y index
 * @iterMax: Maximum simplex iteration budget
 * @sol:     Solved pixel parameters
 */
static void fpwfs_print_solution(
    long                   kx,
    long                   ky,
    long                   iterMax,
    const FPWFS_PIXEL_SOL *sol)
{
    printf("\n\n");
    printf(
        "[%3ld %3ld] OPTIMAL SOLUTION  [ %6ld / %6ld  %g ]: \n",
        kx,
        ky,
        sol->iter,
        iterMax,
        sol->optsize);
    printf("       ptre = %.18f\n", sol->ptre);
    printf("       ptim = %.18f\n", sol->ptim);
    printf("      Iflux = %.18f\n", sol->Iflux);
    if (NBoptVar > 3)
    {
        printf("        are = %.18f\n", sol->are);
    }
    if (NBoptVar > 4)
    {
        printf("        aim = %.18f\n", sol->aim);
        printf("          e = %.18f\n", sol->e);
    }
    printf("OPT VALUES :  %g  %g  %g\n\n\n", tmpvalue1, tmpvalue2, sol->bestvalue);
}

/**
 * fpwfs_run_simplex_pixel - Run GSL Nelder-Mead simplex fit on a single focal-plane pixel
 * @ctx:     Sensitivity analysis context
 * @ptre:    Initial guess for real field component
 * @ptim:    Initial guess for imaginary field component
 * @is_test: Non-zero when diagnostic trace output is enabled
 * @sol:     Output solved parameters
 */
static void fpwfs_run_simplex_pixel(
    AOSIM_FPWFS_SENS_CTX *ctx,
    double                ptre,
    double                ptim,
    int                   is_test,
    FPWFS_PIXEL_SOL      *sol)
{
    for (long optvar = 0; optvar < NBoptVar; optvar++)
    {
        gsl_vector_set(ctx->xvect, (size_t) optvar, 0.0);
    }
    gsl_vector_set(ctx->xvect, 0, ptre);
    gsl_vector_set(ctx->xvect, 1, ptim);
    if (NBoptVar > 3)
    {
        gsl_vector_set(ctx->xvect, 3, 1.0);
    }
    if (NBoptVar > 4)
    {
        gsl_vector_set(ctx->xvect, 4, 0.0);
        gsl_vector_set(ctx->xvect, 5, 1.0);
    }
    gsl_vector_set_all(ctx->ss, 0.1);
    gsl_multimin_fminimizer_set(ctx->s, &ctx->feval_func, ctx->xvect, ctx->ss);

    sol->iter      = 0;
    sol->bestvalue = 100000.0;
    sol->optsize   = 0.0;

    int status;
    do
    {
        sol->iter++;
        status = gsl_multimin_fminimizer_iterate(ctx->s);
        if (status)
        {
            break;
        }
        sol->optsize = gsl_multimin_fminimizer_size(ctx->s);
        double tol   = 1.0 / pow(10.0, 10.0 - 4.0 * (double) sol->iter / (double) ctx->iterMax);
        status       = gsl_multimin_test_size(sol->optsize, tol);
        if (status == GSL_SUCCESS)
        {
            if (is_test)
            {
                printf("  ->  %e [%5ld]  (%e)\n", ctx->s->fval, sol->iter, sol->bestvalue);
            }
            fpwfs_extract_solution(ctx->s->x, sol);
            if (ctx->s->fval < sol->bestvalue)
            {
                sol->bestvalue = ctx->s->fval;
            }
        }
    }
    while ((status == GSL_CONTINUE) && (sol->iter < ctx->iterMax));

    if (sol->iter > ctx->iterMax - 1)
    {
        printf("Max number of iteration reached (optsize = %g)\n", sol->optsize);
    }
    fpwfs_extract_solution(ctx->s->x, sol);
}

/**
 * fpwfs_store_pixel_solution - Record solved parameters into 3D Monte Carlo output cubes
 * @ctx:  Sensitivity analysis context
 * @mapz: Monte Carlo realization slice index
 * @kx:   Pixel X index
 * @ky:   Pixel Y index
 * @sol:  Solved pixel parameters
 */
static void fpwfs_store_pixel_solution(
    const AOSIM_FPWFS_SENS_CTX *ctx,
    long                        mapz,
    long                        kx,
    long                        ky,
    const FPWFS_PIXEL_SOL      *sol)
{
    long   idx2d  = ky * ctx->mapxsize + kx;
    long   idx3d  = mapz * ctx->mapxsize * ctx->mapysize + idx2d;
    double pamp   = dcimg[ctx->IDprobampC].array.F[idx2d];
    double spamp  = sqrt(pamp);

    dcimg[ctx->IDmap].array.F[idx3d]         = (float) (sol->Iflux * spamp);
    dcimg[ctx->IDmap_ptre].array.F[idx3d]    = (float) (sol->ptre * spamp);
    dcimg[ctx->IDmap_ptim].array.F[idx3d]    = (float) (sol->ptim * spamp);
    dcimg[ctx->IDmap_ptre_in].array.F[idx3d] = (float) (sol->ptre_in * spamp);
    dcimg[ctx->IDmap_ptim_in].array.F[idx3d] = (float) (sol->ptim_in * spamp);
    dcimg[ctx->IDmap_Iflux].array.F[idx3d]   = (float) (sol->Iflux * pamp);
    if (NBoptVar > 3)
    {
        dcimg[ctx->IDmap_are].array.F[idx3d] = (float) (sol->are * spamp);
    }
    if (NBoptVar > 4)
    {
        dcimg[ctx->IDmap_aim].array.F[idx3d] = (float) (sol->aim * spamp);
        dcimg[ctx->IDmap_e].array.F[idx3d]   = (float) sol->e;
    }
}

/**
 * aosim_fpwfs_solve_slice - Solve FPWFS electric-field estimation across one 2D slice
 * @ctx:  Sensitivity analysis context
 * @mapz: Monte Carlo realization slice index
 */
void aosim_fpwfs_solve_slice(
    AOSIM_FPWFS_SENS_CTX *ctx,
    long                  mapz)
{
    for (long kx = 0; kx < ctx->mapxsize; kx++)
    {
        for (long ky = 0; ky < ctx->mapysize; ky++)
        {
            printf("\r slice %ld    pixel %5ld %5ld       ", mapz, kx, ky);
            fflush(stdout);

            double ptre;
            double ptim;
            if (ctx->mapmode == 0)
            {
                ptre = 2.0 * ran1() - 1.0;
                ptim = 2.0 * ran1() - 1.0;
            }
            else if (ctx->mapmode == 2)
            {
                ptre = 0.0;
                ptim = 0.0;
            }
            else
            {
                ptre = ctx->mapampl * (2.0 * (double) kx / (ctx->mapxsize - 1) - 1.0);
                ptim = ctx->mapampl * (2.0 * (double) ky / (ctx->mapysize - 1) - 1.0);
            }

            int is_test = ((kx == ctx->kxtest) && (ky == ctx->kytest)) ||
                          ((kx == ctx->kxtest + 1) && (ky == ctx->kytest));
            int verbose = (ctx->mapmode == 0) || is_test;
            if (verbose)
            {
                printf("\n\n ptre, ptim = %g, %g\n", ptre, ptim);
            }

            fpwfs_sample_pixel_probes(ctx, kx, ky, ptre, ptim, is_test);
            if (((ctx->mapmode == 0) && (kx == ctx->kxtest) && (ky == ctx->kytest)) ||
                ((kx == ctx->kxtest + 1) && (ky == ctx->kytest)))
            {
                fflush(stdout);
            }

            FPWFS_PIXEL_SOL sol = {.ptre_in = ptre, .ptim_in = ptim};
            fpwfs_run_simplex_pixel(ctx, ptre, ptim, is_test, &sol);
            if (verbose)
            {
                fpwfs_print_solution(kx, ky, ctx->iterMax, &sol);
            }
            if (ctx->mapmode > 0)
            {
                fpwfs_store_pixel_solution(ctx, mapz, kx, ky, &sol);
            }
        }
    }
}
