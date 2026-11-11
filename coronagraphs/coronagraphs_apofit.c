/**
 * @file    coronagraphs_apofit.c
 * @brief   GSL Nelder-Mead simplex fitting of analytical apodization profiles
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_multimin.h>

#include "statistic/statistic.h"

#include "coronagraphs_internal.h"

double coronagraph_apofit_f_evalmask(
    const gsl_vector *v,
    void             *params);

typedef struct apofit_ctx
{
    long                     nb_opt_var;
    double                  *a_best;
    double                  *b_best;
    double                  *c_best;
    double                   bestvalue;
    double                   keep_limit;
    double                   keep_limit1;
    int                      keep;
    long                     keep_cnt;
    int                      ok;
    gsl_vector              *x;
    gsl_vector              *ss;
    gsl_multimin_fminimizer *s;
    gsl_multimin_function    func;
} apofit_ctx_t;

/**
 * coronagraph_apofit_f_evalmask - GSL objective callback for apodizer profile fitting
 * @v:      Candidate parameter vector (triplets a, b, c)
 * @params: Unused callback context pointer
 *
 * Return: Regularized RMS fit residual.
 */
double coronagraph_apofit_f_evalmask(
    const gsl_vector *v,
    void             *params)
{
    (void) params;
    double value = 0.0;

    for (long n = 0; n < g_coro.fitapoN; n++)
    {
        g_coro.fitapo_a[n] = gsl_vector_get(v, 3 * n);
        g_coro.fitapo_b[n] = gsl_vector_get(v, 3 * n + 1);
        g_coro.fitapo_c[n] = gsl_vector_get(v, 3 * n + 2);
    }

    coronagraph_apofit_eval(&value);

    if (g_coro.LOOPCNT == 0)
    {
        g_coro.optval0 = value;
    }
    g_coro.LOOPCNT++;

    return value;
}

/**
 * apofit_load_initial_best - Load prior best fit from file or set standard defaults
 * @ctx:      Simplex fitting context
 * @fnameout: Path to saved fit parameter file
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE on malformed file.
 */
static errno_t apofit_load_initial_best(
    apofit_ctx_t *ctx,
    const char   *fnameout)
{
    printf("Reading previous best solution [%s] ... ", fnameout);
    fflush(stdout);

    FILE *fpbest = fopen(fnameout, "r");
    int   use_prev = 0;

    if (fpbest != NULL)
    {
        use_prev      = 1;
        int fscanfcnt = fscanf(fpbest, "%lf", &ctx->bestvalue);
        if (coro_check_fscanf(fpbest, fscanfcnt, 1) != RETURN_SUCCESS)
        {
            fclose(fpbest);
            return RETURN_FAILURE;
        }

        printf("PREVIOUS BEST SOLUTION = %g\n", ctx->bestvalue);
        fflush(stdout);

        for (long n = 0; n < g_coro.fitapoN; n++)
        {
            fscanfcnt = fscanf(
                fpbest,
                " %lf %lf %lf",
                &ctx->a_best[n],
                &ctx->b_best[n],
                &ctx->c_best[n]);
            if (coro_check_fscanf(fpbest, fscanfcnt, 3) != RETURN_SUCCESS)
            {
                fclose(fpbest);
                return RETURN_FAILURE;
            }
            g_coro.fitapo_a[n] = ctx->a_best[n];
            g_coro.fitapo_b[n] = ctx->b_best[n];
            g_coro.fitapo_c[n] = ctx->c_best[n];
            if (g_coro.fitapo_c[n] < 0.0)
            {
                use_prev = 0;
            }
        }
        fclose(fpbest);
    }

    if (use_prev == 0)
    {
        printf(" Using std starting point\n");
        fflush(stdout);
        ctx->bestvalue = 100.0;
        ctx->a_best[0] = 1.0;
        ctx->b_best[0] = -1.0;
        ctx->c_best[0] = 1.0;
        for (long n = 1; n < g_coro.fitapoN; n++)
        {
            ctx->a_best[n] = 0.0;
            ctx->b_best[n] = -1.0;
            ctx->c_best[n] = 1.0;
        }
    }

    return RETURN_SUCCESS;
}

/**
 * apofit_select_start_point - Populate simplex starting vector x for outer pass iter1
 * @ctx:   Simplex fitting context
 * @iter1: Outer restart index
 */
static void apofit_select_start_point(
    apofit_ctx_t *ctx,
    long          iter1)
{
    const double a_range = 10.0;
    const double b_range = 5.0;
    const double c_range = 5.0;

    if ((iter1 == 0) || (ctx->ok == 1))
    {
        ctx->keep_cnt = 0;
        printf("BEST  ");
        fflush(stdout);
        for (long n = 0; n < g_coro.fitapoN; n++)
        {
            gsl_vector_set(ctx->x, 3 * n, ctx->a_best[n]);
            gsl_vector_set(ctx->x, 3 * n + 1, ctx->b_best[n]);
            gsl_vector_set(ctx->x, 3 * n + 2, ctx->c_best[n]);
        }
    }
    else if (ctx->keep == 0)
    {
        ctx->keep_cnt = 0;
        printf("RAND  ");
        fflush(stdout);
        for (long n = 0; n < g_coro.fitapoN; n++)
        {
            double v1 = (1.0 - 2.0 * ran1()) * a_range;
            double v2 = (1.0 - 2.0 * ran1()) * b_range;
            double v3 = 0.5 + (2.0 * ran1()) * c_range;
            gsl_vector_set(ctx->x, 3 * n, v1);
            gsl_vector_set(ctx->x, 3 * n + 1, v2);
            gsl_vector_set(ctx->x, 3 * n + 2, v3);
            printf("%lf %lf %lf\n", v1, v2, v3);
        }
    }
    else
    {
        printf("KEEP  ");
        fflush(stdout);
        ctx->keep_cnt++;
        ctx->keep_limit =
            (0.94 * ctx->keep_limit + 0.05 * ctx->bestvalue) * 0.7 + 0.3 * ctx->s->fval;
        for (long n = 0; n < g_coro.fitapoN; n++)
        {
            gsl_vector_set(ctx->x, 3 * n, g_coro.fitapo_a[n]);
            gsl_vector_set(ctx->x, 3 * n + 1, g_coro.fitapo_b[n]);
            gsl_vector_set(ctx->x, 3 * n + 2, g_coro.fitapo_c[n]);
        }
    }

    for (long n = 0; n < g_coro.fitapoN; n++)
    {
        g_coro.fitapo_a[n] = gsl_vector_get(ctx->x, 3 * n);
        g_coro.fitapo_b[n] = gsl_vector_get(ctx->x, 3 * n + 1);
        g_coro.fitapo_c[n] = gsl_vector_get(ctx->x, 3 * n + 2);
    }
}

/**
 * apofit_save_best - Write improved best parameters and gnuplot expression to fnameout
 * @ctx:      Simplex fitting context
 * @fnameout: Output parameter file path
 */
static void apofit_save_best(
    apofit_ctx_t *ctx,
    const char   *fnameout)
{
    ctx->ok        = 1;
    ctx->bestvalue = ctx->s->fval;

    for (long n = 0; n < g_coro.fitapoN; n++)
    {
        ctx->a_best[n] = g_coro.fitapo_a[n];
        ctx->b_best[n] = g_coro.fitapo_b[n];
        ctx->c_best[n] = g_coro.fitapo_c[n];
    }

    FILE *fpbest = fopen(fnameout, "w");
    if (fpbest == NULL)
    {
        return;
    }

    fprintf(fpbest, "%.20g", ctx->bestvalue);
    for (long n = 0; n < g_coro.fitapoN; n++)
    {
        fprintf(fpbest, " %15.10f %15.10f %15.10f", ctx->a_best[n], ctx->b_best[n], ctx->c_best[n]);
    }
    fprintf(fpbest, "\n");
    fprintf(
        fpbest,
        "f(x) = %15.10f * exp(%15.10f * x**%15.10f)",
        ctx->a_best[0],
        ctx->b_best[0],
        ctx->c_best[0]);
    for (long n = 1; n < g_coro.fitapoN; n++)
    {
        fprintf(
            fpbest,
            "+ %15.10f * exp(%15.10f * x**%15.10f)",
            ctx->a_best[n],
            ctx->b_best[n],
            ctx->c_best[n]);
    }
    fprintf(fpbest, "\n");
    fclose(fpbest);
}

/**
 * apofit_run_simplex_pass - Execute one inner Nelder-Mead minimization pass
 * @ctx:      Simplex fitting context
 * @fnameout: Output parameter file path
 */
static void apofit_run_simplex_pass(
    apofit_ctx_t *ctx,
    const char   *fnameout)
{
    const long   iter_max = 50000;
    const double eps1     = 0.99999999;
    long         iter     = 0;
    int          status   = GSL_CONTINUE;

    gsl_vector_set_all(ctx->ss, 1.0e-9);
    printf("I");
    fflush(stdout);
    gsl_multimin_fminimizer_set(ctx->s, &ctx->func, ctx->x, ctx->ss);
    printf("I");
    fflush(stdout);

    ctx->ok        = 0;
    g_coro.LOOPCNT = 0;

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(ctx->s);
        if (status)
        {
            break;
        }

        double size = gsl_multimin_fminimizer_size(ctx->s);
        status      = gsl_multimin_test_size(size, 1e-11);

        if ((status == GSL_SUCCESS) || (ctx->s->fval < ctx->bestvalue * eps1))
        {
            if (status == GSL_SUCCESS)
            {
                printf(
                    " %e ->  %e [%5ld]  (%e)%s\n",
                    g_coro.optval0,
                    ctx->s->fval,
                    iter,
                    ctx->bestvalue,
                    (ctx->ok == 1) ? "  *" : "");
            }

            if (ctx->keep == 0)
            {
                double w         = (ctx->s->fval < ctx->keep_limit1) ? 0.2 : 0.02;
                ctx->keep_limit1 = (1.0 - w) * ctx->keep_limit1 + w * ctx->s->fval;
            }

            if (ctx->s->fval < ctx->bestvalue * eps1)
            {
                apofit_save_best(ctx, fnameout);
            }

            if (ctx->keep == 0)
            {
                ctx->keep       = (ctx->s->fval < ctx->keep_limit1) ? 1 : 0;
                ctx->keep_limit = ctx->keep_limit1;
            }
            else
            {
                ctx->keep = (ctx->s->fval < ctx->keep_limit) ? 1 : 0;
                if (ctx->keep_cnt > 50)
                {
                    ctx->keep = 0;
                }
            }
        }
    }
    while (status == GSL_CONTINUE && iter < iter_max);

    if (iter > iter_max - 1)
    {
        printf("Max number of iteration reached... starting over\n");
        ctx->keep = 0;
    }
}

/**
 * apofit_write_profile - Save raw vs analytical apodizer profile table to "apofitsol.txt"
 */
static void apofit_write_profile(void)
{
    FILE *fp = fopen("apofitsol.txt", "w");
    if (fp == NULL)
    {
        return;
    }

    for (long i = 0; i < g_coro.aporawN; i++)
    {
        double r = g_coro.aporaw_r[i];
        double v = 0.0;
        for (long n = 0; n < g_coro.fitapoN; n++)
        {
            v += g_coro.fitapo_a[n] * exp(g_coro.fitapo_b[n] * pow(r, g_coro.fitapo_c[n]));
        }
        fprintf(fp, "%11.8f %11.8f %11.8f\n", r, g_coro.aporaw_v[i], v);
    }
    fclose(fp);
}

/**
 * coronagraph_apofit - Fit current raw apodization profile with sum of stretched exponentials
 * @fnameout: File path to read/write best parameter set
 *
 * Return: Best achieved RMS objective value, or RETURN_FAILURE on error.
 */
double coronagraph_apofit(const char *fnameout)
{
    const long                          outer_max = 200;
    const gsl_multimin_fminimizer_type *min_type  = gsl_multimin_fminimizer_nmsimplex2rand;
    apofit_ctx_t                        ctx       = {0};

    printf("Fitting apodization profile with analytical function\n");
    fflush(stdout);
    coronagraph_apofit_eval(NULL);

    printf("Preparing optimization ... \n");
    fflush(stdout);
    ctx.nb_opt_var  = 3 * g_coro.fitapoN;
    ctx.x           = gsl_vector_alloc(ctx.nb_opt_var);
    ctx.ss          = gsl_vector_alloc(ctx.nb_opt_var);
    ctx.func.n      = ctx.nb_opt_var;
    ctx.func.f      = &coronagraph_apofit_f_evalmask;
    ctx.func.params = NULL;
    ctx.s           = gsl_multimin_fminimizer_alloc(min_type, ctx.nb_opt_var);

    ctx.a_best = (double *) malloc(sizeof(double) * g_coro.fitapoN);
    ctx.b_best = (double *) malloc(sizeof(double) * g_coro.fitapoN);
    ctx.c_best = (double *) malloc(sizeof(double) * g_coro.fitapoN);
    if ((ctx.a_best == NULL) || (ctx.b_best == NULL) || (ctx.c_best == NULL))
    {
        PRINT_ERROR("malloc returns NULL pointer");
        free(ctx.a_best);
        free(ctx.b_best);
        free(ctx.c_best);
        gsl_vector_free(ctx.x);
        gsl_vector_free(ctx.ss);
        gsl_multimin_fminimizer_free(ctx.s);
        return RETURN_FAILURE;
    }

    if (apofit_load_initial_best(&ctx, fnameout) != RETURN_SUCCESS)
    {
        free(ctx.a_best);
        free(ctx.b_best);
        free(ctx.c_best);
        gsl_vector_free(ctx.x);
        gsl_vector_free(ctx.ss);
        gsl_multimin_fminimizer_free(ctx.s);
        return RETURN_FAILURE;
    }

    printf("----------- STARTING DOWNHILL SIMPLEX ------------- [%ld]\n", outer_max);
    ctx.ok          = 0;
    ctx.keep_limit  = ctx.bestvalue;
    ctx.keep_limit1 = 3.0 * ctx.bestvalue;
    ctx.keep        = 0;
    ctx.keep_cnt    = 0;

    for (long iter1 = 0; iter1 < outer_max; iter1++)
    {
        printf("%05ld [KL %e %e]  ", iter1, ctx.keep_limit, ctx.keep_limit1);
        fflush(stdout);
        apofit_select_start_point(&ctx, iter1);
        apofit_run_simplex_pass(&ctx, fnameout);
    }

    gsl_vector_free(ctx.x);
    gsl_vector_free(ctx.ss);
    gsl_multimin_fminimizer_free(ctx.s);

    for (long n = 0; n < g_coro.fitapoN; n++)
    {
        g_coro.fitapo_a[n] = ctx.a_best[n];
        g_coro.fitapo_b[n] = ctx.b_best[n];
        g_coro.fitapo_c[n] = ctx.c_best[n];
    }

    free(ctx.a_best);
    free(ctx.b_best);
    free(ctx.c_best);

    apofit_write_profile();
    return ctx.bestvalue;
}
