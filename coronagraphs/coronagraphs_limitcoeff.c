/**
 * @file    coronagraphs_limitcoeff.c
 * @brief   Fundamental coronagraphic inner-working-angle Taylor-series limit calculation
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "statistic/statistic.h"

#include "coronagraphs_internal.h"

typedef struct limitcoeff_ctx
{
    long    csize;
    long    nb_coeff;
    long    nb_pts;
    double  xmax;
    double  xeval;
    double *coeffmax;
    double *coeff;
    double *bestcoeff;
    double *xpow;
    double *xval;
    double *pval;
    double *pvalmax;
} limitcoeff_ctx_t;

/**
 * limitcoeff_free - Release all working buffers in limitcoeff_ctx_t
 * @ctx: Context structure to clean up
 */
static void limitcoeff_free(limitcoeff_ctx_t *ctx)
{
    free(ctx->coeffmax);
    free(ctx->coeff);
    free(ctx->bestcoeff);
    free(ctx->xpow);
    free(ctx->xval);
    free(ctx->pval);
    free(ctx->pvalmax);
}

/**
 * limitcoeff_alloc - Allocate Taylor-series bound tables and random-search working buffers
 * @ctx: Context structure to populate
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE.
 */
static errno_t limitcoeff_alloc(limitcoeff_ctx_t *ctx)
{
    ctx->coeffmax  = (double *) malloc(sizeof(double) * ctx->nb_coeff);
    ctx->coeff     = (double *) malloc(sizeof(double) * ctx->nb_coeff);
    ctx->bestcoeff = (double *) malloc(sizeof(double) * ctx->nb_coeff);
    ctx->xval      = (double *) malloc(sizeof(double) * ctx->nb_pts);
    ctx->pval      = (double *) malloc(sizeof(double) * ctx->nb_pts);
    ctx->xpow      = (double *) malloc(sizeof(double) * ctx->nb_coeff * ctx->nb_pts);
    ctx->pvalmax   = (double *) malloc(sizeof(double) * ctx->nb_pts);

    if ((ctx->coeffmax == NULL) || (ctx->coeff == NULL) || (ctx->bestcoeff == NULL) ||
        (ctx->xval == NULL) || (ctx->pval == NULL) || (ctx->xpow == NULL) || (ctx->pvalmax == NULL))
    {
        limitcoeff_free(ctx);
        PRINT_ERROR("malloc returns NULL pointer");
        return RETURN_FAILURE;
    }

    return RETURN_SUCCESS;
}

/**
 * limitcoeff_compute_bounds - Compute moment upper bounds and print sum-of-terms envelope
 * @ctx: Initialized limitcoeff context
 */
static void limitcoeff_compute_bounds(limitcoeff_ctx_t *ctx)
{
    double fact = 1.0;
    for (long k = 0; k < ctx->nb_coeff; k++)
    {
        if (k != 0)
        {
            fact *= k;
        }
        double value = 0.0;
        long   cnt   = 0;
        for (long ii = 0; ii < ctx->csize; ii++)
        {
            for (long jj = 0; jj < ctx->csize; jj++)
            {
                double x = 1.0 * ii / ctx->csize;
                double y = 1.0 * jj / ctx->csize;
                if (sqrt(x * x + y * y) < 1.0)
                {
                    value += pow(x, 2 * k);
                    cnt++;
                }
            }
        }
        value /= cnt;
        ctx->coeffmax[k] = sqrt(value) * pow(CORO_PI, k) / fact;
    }

    for (long i = 0; i < ctx->nb_pts; i++)
    {
        ctx->xval[i]    = 1.0 * i * ctx->xmax / ctx->nb_pts;
        ctx->pvalmax[i] = 0.0;
        for (long k = 0; k < ctx->nb_coeff; k++)
        {
            ctx->xpow[i * ctx->nb_coeff + k] = ctx->coeffmax[k] * pow(ctx->xval[i], k);
        }
    }

    for (long i = 0; i < ctx->nb_pts; i++)
    {
        double tmp = 0.0;
        for (long k = 1; k < ctx->nb_coeff; k++)
        {
            tmp += ctx->xpow[i * ctx->nb_coeff + k];
        }
        printf("%f %g\n", ctx->xval[i], tmp);
    }
}

/**
 * limitcoeff_random_search - Run Monte Carlo coefficient search for peak throughput at xeval
 * @ctx: Initialized limitcoeff context
 */
static void limitcoeff_random_search(limitcoeff_ctx_t *ctx)
{
    for (long k = 0; k < ctx->nb_coeff; k++)
    {
        ctx->bestcoeff[k] = 0.0;
    }

    const long i0 = (long) (ctx->xeval / ctx->xmax * ctx->nb_pts);
    while (1)
    {
        double amp = pow(ran1(), 5.0);
        for (long k = 0; k < ctx->nb_coeff; k++)
        {
            double c = ctx->bestcoeff[k] + amp * (ran1() - 0.5);
            if (c > 1.0)
            {
                c = 1.0;
            }
            if (c < -1.0)
            {
                c = -1.0;
            }
            ctx->coeff[k] = c;
        }
        ctx->coeff[0] = 0.0;

        double maxpval = 0.0;
        for (long i = 0; i < ctx->nb_pts; i++)
        {
            double sum = 0.0;
            for (long k = 0; k < ctx->nb_coeff; k++)
            {
                sum += ctx->xpow[i * ctx->nb_coeff + k] * ctx->coeff[k];
            }
            ctx->pval[i] = sum;
            if (fabs(sum) > fabs(maxpval))
            {
                maxpval = sum;
            }
        }
        for (long i = 0; i < ctx->nb_pts; i++)
        {
            ctx->pval[i] /= maxpval;
        }

        if (ctx->pval[i0] > ctx->pvalmax[i0])
        {
            for (long i = 0; i < ctx->nb_pts; i++)
            {
                ctx->pvalmax[i] = ctx->pval[i];
            }
            for (long k = 0; k < ctx->nb_coeff; k++)
            {
                ctx->bestcoeff[k] = ctx->coeff[k];
            }

            printf("%g\n", ctx->pvalmax[i0] * ctx->pvalmax[i0]);
            fflush(stdout);
            FILE *fp = fopen("transm.txt", "w");
            for (long i = 0; i < ctx->nb_pts; i++)
            {
                fprintf(fp, "%f %f\n", ctx->xval[i], ctx->pvalmax[i] * ctx->pvalmax[i]);
            }
            fclose(fp);

            fp = fopen("bestcoeff.txt", "w");
            for (long k = 0; k < ctx->nb_coeff; k++)
            {
                fprintf(fp, "%ld %f\n", k, fabs(ctx->coeff[k]));
            }
            fclose(fp);
        }
    }
}

/**
 * coronagraph_compute_limitcoeff - Compute fundamental coronagraphic contrast/throughput limits
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE.
 */
errno_t coronagraph_compute_limitcoeff(void)
{
    limitcoeff_ctx_t ctx = {
        .csize    = 200,
        .nb_coeff = 50,
        .nb_pts   = 100,
        .xmax     = 5.0,
        .xeval    = 0.5,
    };

    if (limitcoeff_alloc(&ctx) != RETURN_SUCCESS)
    {
        return RETURN_FAILURE;
    }

    limitcoeff_compute_bounds(&ctx);
    exit(0);

    limitcoeff_random_search(&ctx);
    limitcoeff_free(&ctx);
    return RETURN_SUCCESS;
}
