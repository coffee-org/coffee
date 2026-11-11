/**
 * @file    coronagraphs_scan.c
 * @brief   PIAACMC central-obscuration contrast, IWA, and throughput design-space scanner
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "coffee_compat.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "coronagraphs_internal.h"

typedef struct scan_out_maps
{
    FILE   *fpresult;
    imageID id_transm;
    imageID id_iwa;
    imageID id_contrastmax;
    imageID id_fpmasktransm;
} scan_out_maps_t;

/**
 * scan_ctx_free - Free all dynamically allocated arrays inside coro_scan_ctx_t
 * @ctx: Scan context to release
 */
static void scan_ctx_free(coro_scan_ctx_t *ctx)
{
    free(ctx->separray);
    free(ctx->transmarray);
    free(ctx->starprofrad);
    free(ctx->starprofarray);
    free(ctx->starprofarraycnt);
    free(ctx->obs1array);
    free(ctx->fpmaskradarray);
}

/**
 * scan_ctx_init - Initialize grid axes, separation samples, and buffers for PIAACMC scan
 * @ctx:       Output scan context
 * @obs0input: Signed input central obscuration parameter
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE.
 */
static errno_t scan_ctx_init(
    coro_scan_ctx_t *ctx,
    double           obs0input)
{
    static const double sep_init[18] = {0.0, 0.005, 0.01, 0.02, 0.1, 0.2, 0.5, 1.0, 1.5,
                                        2.0, 2.5,   3.0,  3.5,  4.0, 5.0, 6.0, 7.0, 8.0};
    const float obs1min       = 0.0f;
    const float obs1max       = 0.4f;
    const float obs1step      = 0.05f;
    const float fpmaskradmin  = 0.8f;
    const float fpmaskradmax  = 3.0f;
    const float fpmaskradstep = 0.1f;

    memset(ctx, 0, sizeof(*ctx));
    ctx->obs0input          = obs0input;
    ctx->obs0               = fabs(obs0input);
    ctx->kmax               = 18;
    ctx->size               = 256;
    ctx->pixscale           = 10.25;
    ctx->starprof_step      = 0.1;
    ctx->starprof_maxrad    = 10.0;
    ctx->starprofNBstep     = (long) (ctx->starprof_maxrad / ctx->starprof_step);
    ctx->star_radius        = 0.01;
    ctx->obs1arraysize      = (long) ((obs1max - obs1min) / obs1step);
    ctx->fpmaskradarraysize = (long) ((fpmaskradmax - fpmaskradmin) / fpmaskradstep);
    ctx->compute_data       = 1;

    ctx->starprofrad      = (double *) malloc(sizeof(double) * ctx->starprofNBstep);
    ctx->starprofarray    = (double *) malloc(sizeof(double) * ctx->starprofNBstep);
    ctx->starprofarraycnt = (double *) malloc(sizeof(double) * ctx->starprofNBstep);
    ctx->transmarray      = (double *) malloc(sizeof(double) * ctx->kmax);
    ctx->separray         = (double *) malloc(sizeof(double) * ctx->kmax);
    ctx->obs1array        = (float *) malloc(sizeof(float) * ctx->obs1arraysize);
    ctx->fpmaskradarray   = (float *) malloc(sizeof(float) * ctx->fpmaskradarraysize);

    if ((ctx->starprofrad == NULL) || (ctx->starprofarray == NULL) ||
        (ctx->starprofarraycnt == NULL) || (ctx->transmarray == NULL) || (ctx->separray == NULL) ||
        (ctx->obs1array == NULL) || (ctx->fpmaskradarray == NULL))
    {
        scan_ctx_free(ctx);
        PRINT_ERROR("malloc returns NULL pointer");
        return RETURN_FAILURE;
    }

    memcpy(ctx->separray, sep_init, sizeof(sep_init));
    for (long ii = 0; ii < ctx->obs1arraysize; ii++)
    {
        ctx->obs1array[ii] = obs1min + obs1step * ii;
    }
    for (long ii = 0; ii < ctx->fpmaskradarraysize; ii++)
    {
        ctx->fpmaskradarray[ii] = fpmaskradmin + fpmaskradstep * ii;
    }

    FILE *fp = fopen("obs1indexstart.txt", "r");
    if (fp != NULL)
    {
        int fscanfcnt = fscanf(fp, "%ld", &ctx->obs1_index_start);
        if (coro_check_fscanf(fp, fscanfcnt, 1) != RETURN_SUCCESS)
        {
            fclose(fp);
            scan_ctx_free(ctx);
            return RETURN_FAILURE;
        }
        fclose(fp);
    }

    return RETURN_SUCCESS;
}

/**
 * scan_load_prol_info - Read apodizer .info file and record focal-plane mask transmission
 * @ctx:              Scan context
 * @maps:             Output summary images
 * @obs1_idx:         Output-obscuration grid index
 * @fpm_idx:          Focal-plane mask radius grid index
 * @inout_fpm_transm: Updated focal-plane mask complex-amplitude transmission if file exists
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE on malformed file.
 */
static errno_t scan_load_prol_info(
    const coro_scan_ctx_t *ctx,
    const scan_out_maps_t *maps,
    long                   obs1_idx,
    long                   fpm_idx,
    double                *inout_fpm_transm)
{
    double      obs1      = ctx->obs1array[obs1_idx];
    double      fpmaskrad = ctx->fpmaskradarray[fpm_idx];
    struct stat buffer;
    char        fname[500];

    if (ctx->cent_obs_flag == 0)
    {
        snprintf(fname,
                 sizeof(fname),
                 "%s/APLCapo/APLCapo_%5.3f.4096.ref.info",
                 CORONAGRAPHSDATALOCAL,
                 fpmaskrad);
    }
    else
    {
        snprintf(fname,
                 sizeof(fname),
                 "%s/APLCapo/APLCapo_%5.3f.%5.3f.4096.ref.info",
                 CORONAGRAPHSDATALOCAL,
                 fpmaskrad,
                 obs1);
    }

    if (stat(fname, &buffer) != 0)
    {
        return RETURN_SUCCESS;
    }

    FILE *fp = fopen(fname, "r");
    float t0 = 0.0f;
    float t1 = 0.0f;
    float t2 = 0.0f;
    float t3 = 0.0f;
    int   rc = fscanf(fp, "%f %f %f %f\n", &t0, &t1, &t2, &t3);
    if (coro_check_fscanf(fp, rc, 4) != RETURN_SUCCESS)
    {
        fclose(fp);
        return RETURN_FAILURE;
    }
    fclose(fp);

    *inout_fpm_transm = -(1.0 - t3) / t3;
    dcimg[maps->id_fpmasktransm].array.F[fpm_idx * ctx->obs1arraysize + obs1_idx] =
        (float) (*inout_fpm_transm);
    return RETURN_SUCCESS;
}

/**
 * scan_load_star_contrast - Read stellar profile and record worst contrast in [2.5, 2.8] l/D
 * @ctx:        Scan context
 * @maps:       Output summary images
 * @obs1_idx:   Output-obscuration grid index
 * @fpm_idx:    Focal-plane mask radius grid index
 * @out_maxval: Output maximum contrast on evaluation annulus
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE on malformed file.
 */
static errno_t scan_load_star_contrast(
    const coro_scan_ctx_t *ctx,
    const scan_out_maps_t *maps,
    long                   obs1_idx,
    long                   fpm_idx,
    double                *out_maxval)
{
    const double eval_start = 2.5;
    const double eval_end   = 2.8;
    double       obs1       = ctx->obs1array[obs1_idx];
    double       fpmaskrad  = ctx->fpmaskradarray[fpm_idx];
    double       maxval     = 0.0;
    struct stat  buffer;
    char         profname[500];

    snprintf(profname,
             sizeof(profname),
             "./PROF/star_%5.3f_%5.3f_%5.3f.prof",
             fpmaskrad,
             ctx->obs0,
             obs1);
    if (stat(profname, &buffer) == 0)
    {
        printf("Reading file %s\n", profname);
        FILE *fp = fopen(profname, "r");
        for (long i = 0; i < ctx->starprofNBstep; i++)
        {
            float t0 = 0.0f;
            float t1 = 0.0f;
            int   rc = fscanf(fp, "%f %f\n", &t0, &t1);
            if (coro_check_fscanf(fp, rc, 2) != RETURN_SUCCESS)
            {
                fclose(fp);
                return RETURN_FAILURE;
            }
            if ((t0 > eval_start) && (t0 < eval_end) && (t1 > maxval))
            {
                maxval = t1;
            }
        }
        fclose(fp);
    }

    printf("== %s ==   %f  %f  %f          %e   [%ld %ld]\n",
           profname,
           fpmaskrad,
           ctx->obs0,
           obs1,
           maxval,
           fpm_idx,
           obs1_idx);
    dcimg[maps->id_contrastmax].array.F[fpm_idx * ctx->obs1arraysize + obs1_idx] = (float) maxval;
    *out_maxval                                                                  = maxval;
    return RETURN_SUCCESS;
}

/**
 * scan_load_transm_iwa - Read separation throughput profile, interpolate IWA50, and update cube
 * @ctx:      Scan context
 * @maps:     Output summary images
 * @obs1_idx: Output-obscuration grid index
 * @fpm_idx:  Focal-plane mask radius grid index
 * @out_iwa:  Output interpolated 50%-throughput inner working angle [lambda/D]
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE on malformed file.
 */
static errno_t scan_load_transm_iwa(
    coro_scan_ctx_t       *ctx,
    const scan_out_maps_t *maps,
    long                   obs1_idx,
    long                   fpm_idx,
    float                 *out_iwa)
{
    double      obs1      = ctx->obs1array[obs1_idx];
    double      fpmaskrad = ctx->fpmaskradarray[fpm_idx];
    struct stat buffer;
    char        transmname[500];

    snprintf(transmname,
             sizeof(transmname),
             "./TRANSM/transm_%5.3f_%5.3f_%5.3f.prof",
             fpmaskrad,
             ctx->obs0,
             obs1);
    if (stat(transmname, &buffer) != 0)
    {
        *out_iwa = 0.0f;
        return RETURN_SUCCESS;
    }

    printf("Reading file %s\n", transmname);
    FILE *fp   = fopen(transmname, "r");
    long  cell = fpm_idx * ctx->obs1arraysize + obs1_idx;
    long   area = ctx->fpmaskradarraysize * ctx->obs1arraysize;
    for (long k = 0; k < ctx->kmax; k++)
    {
        float t0 = 0.0f;
        float t1 = 0.0f;
        int   rc = fscanf(fp, "%f %f\n", &t0, &t1);
        if (coro_check_fscanf(fp, rc, 2) != RETURN_SUCCESS)
        {
            fclose(fp);
            return RETURN_FAILURE;
        }
        dcimg[maps->id_transm].array.F[k * area + cell] = t1;
        ctx->transmarray[k]                             = t1;
    }
    fclose(fp);

    long  k     = 1;
    float tmpf0 = (float) ctx->transmarray[k - 1];
    float tmpf1 = (float) ctx->transmarray[k];
    while ((tmpf1 < 0.5f) && (k < ctx->kmax - 1))
    {
        tmpf0 = (float) ctx->transmarray[k - 1];
        tmpf1 = (float) ctx->transmarray[k];
        k++;
    }
    k--;
    float alpha =
        (float) ((0.5 - ctx->transmarray[k - 1]) / (ctx->transmarray[k] - ctx->transmarray[k - 1]));
    float iwa = (float) (ctx->separray[k - 1] + alpha * (ctx->separray[k] - ctx->separray[k - 1]));
    dcimg[maps->id_iwa].array.F[cell] = iwa;
    printf("%ld %f %f   %ld %f %f  -> [%f %f] %f %f\n",
           k - 1,
           ctx->separray[k - 1],
           tmpf0,
           k,
           ctx->separray[k],
           tmpf1,
           0.5 - ctx->transmarray[k - 1],
           ctx->transmarray[k] - ctx->transmarray[k - 1],
           alpha,
           iwa);
    *out_iwa = iwa;
    return RETURN_SUCCESS;
}

/**
 * scan_save_summary_fits - Save 3D transmission cube and 2D performance maps to FITS
 * @obs0input: Input central obscuration parameter used in output filenames
 */
static void scan_save_summary_fits(double obs0input)
{
    char fname[500];

    snprintf(fname, sizeof(fname), "PIAACMCtransm_%5.3f.fits", obs0input);
    save_fl_fits("PIAACMCtransm", fname);
    snprintf(fname, sizeof(fname), "PIAACMCiwa_%5.3f.fits", obs0input);
    save_fl_fits("PIAACMCiwa", fname);
    snprintf(fname, sizeof(fname), "PIAACMCcontrastmax_%5.3f.fits", obs0input);
    save_fl_fits("PIAACMCcontrastmax", fname);
    snprintf(fname, sizeof(fname), "PIAACMCfpmasktransm_%5.3f.fits", obs0input);
    save_fl_fits("PIAACMCfpmasktransm", fname);
}

/**
 * CORONAGRAPHS_scanPIAACMC_centObs_perf - Scan PIAACMC performance vs. central obscuration
 * @obs0input: Input pupil central obscuration fraction
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE.
 */
errno_t CORONAGRAPHS_scanPIAACMC_centObs_perf(double obs0input)
{
    coro_scan_ctx_t ctx;

    printf("RUNNING PIAACMC DESIGN SCAN FOR CENTRAL OBSTRUCTION %f\n", obs0input);
    fflush(stdout);

    if (scan_ctx_init(&ctx, obs0input) != RETURN_SUCCESS)
    {
        return RETURN_FAILURE;
    }

    if (ctx.compute_data == 1)
    {
        coro_scan_compute_grid(&ctx);
    }

    scan_out_maps_t maps = {0};
    maps.fpresult        = fopen("result.log", "w");
    create_3Dimage_ID(
        "PIAACMCtransm", ctx.obs1arraysize, ctx.fpmaskradarraysize, ctx.kmax, &maps.id_transm);
    create_2Dimage_ID("PIAACMCiwa", ctx.obs1arraysize, ctx.fpmaskradarraysize, &maps.id_iwa);
    create_2Dimage_ID(
        "PIAACMCcontrastmax", ctx.obs1arraysize, ctx.fpmaskradarraysize, &maps.id_contrastmax);
    create_2Dimage_ID(
        "PIAACMCfpmasktransm", ctx.obs1arraysize, ctx.fpmaskradarraysize, &maps.id_fpmasktransm);

    double fpmasktransm = 0.0;
    for (long obs1_idx = 0; obs1_idx < ctx.obs1arraysize; obs1_idx++)
    {
        double obs1 = ctx.obs1array[obs1_idx];
        if (obs1 > 0.001)
        {
            ctx.cent_obs_flag = 1;
            printf("Central obstruction = %f\n", obs1);
        }

        for (long fpm_idx = 0; fpm_idx < ctx.fpmaskradarraysize; fpm_idx++)
        {
            double maxval = 0.0;
            float  iwa    = 0.0f;
            if ((scan_load_prol_info(&ctx, &maps, obs1_idx, fpm_idx, &fpmasktransm) !=
                 RETURN_SUCCESS) ||
                (scan_load_star_contrast(&ctx, &maps, obs1_idx, fpm_idx, &maxval) !=
                 RETURN_SUCCESS) ||
                (scan_load_transm_iwa(&ctx, &maps, obs1_idx, fpm_idx, &iwa) != RETURN_SUCCESS))
            {
                fclose(maps.fpresult);
                scan_ctx_free(&ctx);
                return RETURN_FAILURE;
            }

            if (iwa > 0.5f)
            {
                fprintf(maps.fpresult,
                        "%5.3f %5.3f %g %g %g\n",
                        obs1,
                        ctx.fpmaskradarray[fpm_idx],
                        iwa,
                        maxval,
                        fpmasktransm);
            }
        }
    }
    fclose(maps.fpresult);

    scan_save_summary_fits(obs0input);
    scan_ctx_free(&ctx);
    return RETURN_SUCCESS;
}
