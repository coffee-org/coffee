/**
 * @file    coronagraphs_apofit_eval.c
 * @brief   Analytical multi-Gaussian/exponential apodization profile evaluation
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "coronagraphs_internal.h"

/**
 * coronagraph_apofit_eval_r - Evaluate analytical apodization fit at normalized radius r
 * @r:      Normalized pupil radius
 * @outval: Output apodization amplitude (may be NULL)
 *
 * Return: RETURN_SUCCESS on completion.
 */
errno_t coronagraph_apofit_eval_r(
    double  r,
    double *outval)
{
    DEBUG_TRACE_FSTART();

    const double eps  = 1.0e-8;
    const double eps1 = 1.0e-4;
    const double a    = 8.0;
    double       rp;

    if (r > 1.0 + eps)
    {
        rp = r - 1.0;
        rp = rp / (1.0 + a * rp * exp(-1.0 / (a * rp)));
        rp += 1.0;
    }
    else
    {
        rp = r;
    }

    if (rp < eps1)
    {
        rp = eps1;
    }

    double value = 0.0;
    for (long n = 0; n < g_coro.fitapoN; n++)
    {
        value += g_coro.fitapo_a[n] * exp(g_coro.fitapo_b[n] * pow(rp, g_coro.fitapo_c1[n]));
    }

    if (isfinite(value) == 0)
    {
        printf("ERROR : infinite value in coronagraph_apofit_eval_r\n");
        value = 0.0;
        for (long n = 0; n < g_coro.fitapoN; n++)
        {
            printf(
                "        :::: %ld  %f %f %f :::::\n",
                n,
                g_coro.fitapo_a[n],
                g_coro.fitapo_b[n],
                g_coro.fitapo_c[n]);
        }

        for (long n = 0; n < g_coro.fitapoN; n++)
        {
            double pval = pow(rp, g_coro.fitapo_c[n]);
            double eval = exp(g_coro.fitapo_b[n] * pval);
            value += g_coro.fitapo_a[n] * eval;
            printf(
                "---- %ld %f   %f %f %f -> %f %f -> %f\n",
                n,
                rp,
                g_coro.fitapo_a[n],
                g_coro.fitapo_b[n],
                g_coro.fitapo_c[n],
                pval,
                eval,
                value);
            fflush(stdout);
        }
    }

    if (outval != NULL)
    {
        *outval = value;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * coronagraph_apofit_eval - Evaluate regularized RMS residual of current apodizer fit
 * @outval: Output residual cost value (may be NULL)
 *
 * Return: RETURN_SUCCESS on completion.
 */
errno_t coronagraph_apofit_eval(double *outval)
{
    DEBUG_TRACE_FSTART();

    double value = 0.0;

    for (long n = 0; n < g_coro.fitapoN; n++)
    {
        double c_pen =
            0.5 * (atan(-(g_coro.fitapo_c[n] - 0.5) * 10.0 / (CORO_PI / 2.0)) / (CORO_PI / 2.0) +
                   1.0);
        value += 100.0 * pow(c_pen, 8.0);

        g_coro.fitapo_c1[n] = g_coro.fitapo_c[n];
        if (g_coro.fitapo_c1[n] < g_coro.fitapo_minc)
        {
            g_coro.fitapo_c1[n] = g_coro.fitapo_minc;
        }

        value += exp(-1.0 / (g_coro.fitapo_a[n] * g_coro.fitapo_a[n] / 200.0 / 200.0));
        value += exp(-1.0 / (g_coro.fitapo_b[n] * g_coro.fitapo_b[n] / 200.0 / 200.0));
    }

    for (long i = 1; i < g_coro.aporawN - 1; i++)
    {
        double v    = 0.0;
        double r    = g_coro.aporaw_r[i];
        coronagraph_apofit_eval_r(r, &v);
        double tmpv = g_coro.aporaw_v[i] - v;
        value += tmpv * tmpv;
    }
    value = sqrt(value / g_coro.aporawN);

    if (isfinite(value) == 0)
    {
        printf("ERROR: INFINITE VALUE in coronagraph_apofit_eval ----\n");
        fflush(stdout);

        for (long n = 0; n < g_coro.fitapoN; n++)
        {
            printf(
                "   :::: %ld  %f %f %f :::::\n",
                n,
                g_coro.fitapo_a[n],
                g_coro.fitapo_b[n],
                g_coro.fitapo_c[n]);
        }

        printf("aporawN = %ld\n", g_coro.aporawN);
        value = 0.0;
        for (long i = 0; i < g_coro.aporawN; i++)
        {
            double v    = 0.0;
            double r    = g_coro.aporaw_r[i];
            coronagraph_apofit_eval_r(r, &v);
            double tmpv = g_coro.aporaw_v[i] - v;
            value += tmpv * tmpv;
            printf("%ld %f %f %f %f\n", i, r, v, tmpv, value);
            fflush(stdout);
        }
        value = sqrt(value / g_coro.aporawN);
        exit(0);
    }

    if (outval != NULL)
    {
        *outval = value;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
