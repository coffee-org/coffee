/**
 * @file    coronagraphs_scan_compute.c
 * @brief   Grid PSF generation and radial/throughput profile synthesis for PIAACMC scan
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "coffee_compat.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "coronagraphs_internal.h"

/**
 * scan_write_aplc_conf - Write APLCPIAA.conf for a single (fpmaskrad, obs0, obs1) design point
 * @fpmaskrad: Focal-plane mask radius [lambda/D]
 * @obs0:      Input pupil central obscuration ratio
 * @obs1:      Output pupil central obscuration ratio
 */
static void scan_write_aplc_conf(
    double fpmaskrad,
    double obs0,
    double obs1)
{
    FILE *fpconf = fopen("APLCPIAA.conf", "w");
    fprintf(fpconf, "1 %5.3f-%5.3f-%5.3f 2 1\n", fpmaskrad, obs0, obs1);
    printf("1 %5.3f-%5.3f-%5.3f 2 1\n", fpmaskrad, obs0, obs1);
    fclose(fpconf);
}

/**
 * scan_ensure_cell_psfs - Generate any missing separation PSFs and load psf00..psf17 into memory
 * @ctx:       Scan context
 * @fpmaskrad: Focal-plane mask radius [lambda/D]
 * @obs1:      Output central obscuration ratio
 *
 * Return: 1 if all separation PSFs are loaded, 0 if skipped.
 */
static int scan_ensure_cell_psfs(
    const coro_scan_ctx_t *ctx,
    double                 fpmaskrad,
    double                 obs1)
{
    const int   skip_psf_comp = 0;
    int         psfsok        = 1;
    struct stat buffer;
    char        fname[500];
    char        fname1[64];

    for (long k = 0; k < ctx->kmax; k++)
    {
        if (stat("STOP", &buffer) == 0)
        {
            EXECUTE_SYSTEM_COMMAND("rm STOP");
            exit(0);
        }

        snprintf(fname,
                 sizeof(fname),
                 "./PSF/PSF_%5.3f_%5.3f_%5.3f_%5.3f.fits",
                 fpmaskrad,
                 ctx->obs0,
                 obs1,
                 ctx->separray[k]);

        if (stat(fname, &buffer) != 0)
        {
            if (skip_psf_comp == 1)
            {
                psfsok = 0;
            }
            else
            {
                printf("============= CREATING %s ===========\n", fname);
                scan_write_aplc_conf(fpmaskrad, ctx->obs0, obs1);
                EXECUTE_SYSTEM_COMMAND("./mkpsf %f %f", ctx->separray[k], 0.0);
                EXECUTE_SYSTEM_COMMAND("mv psf.fits %s", fname);
            }
        }

        if (psfsok == 1)
        {
            snprintf(fname1, sizeof(fname1), "psf%02ld", k);
            load_fits(fname, fname1, 1, NULL);
        }
    }

    return psfsok;
}

/**
 * scan_write_star_profile - Compute azimuthally binned finite-star contrast profile and save ASCII
 * @ctx:       Scan context
 * @fpmaskrad: Focal-plane mask radius [lambda/D]
 * @obs1:      Output central obscuration ratio
 */
static void scan_write_star_profile(
    coro_scan_ctx_t *ctx,
    double           fpmaskrad,
    double           obs1)
{
    for (long i = 0; i < ctx->starprofNBstep; i++)
    {
        ctx->starprofarraycnt[i] = 0.0;
        ctx->starprofarray[i]    = 0.0;
    }

    const long k_star = 2;
    char       fname1[64];
    snprintf(fname1, sizeof(fname1), "psf%02ld", k_star);
    imageID id = image_ID(fname1, dcimg, dcnimg);

    for (long ii = 0; ii < ctx->size; ii++)
    {
        for (long jj = 0; jj < ctx->size; jj++)
        {
            double x = (1.0 * ii - ctx->size / 2) / ctx->pixscale;
            double y = (1.0 * jj - ctx->size / 2) / ctx->pixscale;
            double r = sqrt(x * x + y * y);
            long   i = (long) (r / ctx->starprof_step + 0.5);
            if (i < ctx->starprofNBstep)
            {
                ctx->starprofarray[i] += dcimg[id].array.F[jj * ctx->size + ii];
                ctx->starprofarraycnt[i] += 1.0;
            }
        }
    }

    snprintf(fname1, sizeof(fname1), "psf%02d", (int) (ctx->kmax - 1));
    double peak = arith_image_max(fname1);
    printf("PEAK = %f\n", peak);

    char profname[500];
    snprintf(profname,
             sizeof(profname),
             "./PROF/star_%5.3f_%5.3f_%5.3f.prof",
             fpmaskrad,
             ctx->obs0,
             obs1);
    FILE *fp          = fopen(profname, "w");
    double scale_mult = 0.5 * pow(ctx->star_radius / ctx->separray[k_star], 2.0) / peak;
    for (long i = 0; i < ctx->starprofNBstep; i++)
    {
        ctx->starprofarray[i] /= (ctx->starprofarraycnt[i] + 0.000001);
        ctx->starprofarray[i] *= scale_mult;
        fprintf(fp, "%5.3f %g\n", ctx->starprof_step * i, ctx->starprofarray[i]);
    }
    fclose(fp);
}

/**
 * scan_write_transm_and_cleanup - Save normalized separation transmission curve and free PSF images
 * @ctx:       Scan context
 * @fpmaskrad: Focal-plane mask radius [lambda/D]
 * @obs1:      Output central obscuration ratio
 */
static void scan_write_transm_and_cleanup(
    coro_scan_ctx_t *ctx,
    double           fpmaskrad,
    double           obs1)
{
    struct stat buffer;
    char        transmname[500];
    char        fname1[64];

    snprintf(transmname,
             sizeof(transmname),
             "./TRANSM/transm_%5.3f_%5.3f_%5.3f.prof",
             fpmaskrad,
             ctx->obs0,
             obs1);
    if (stat(transmname, &buffer) != 0)
    {
        float tmax = 0.0f;
        for (long k = 0; k < ctx->kmax; k++)
        {
            snprintf(fname1, sizeof(fname1), "psf%02ld", k);
            ctx->transmarray[k] = arith_image_total(fname1) / (1.0 - ctx->obs0 * ctx->obs0);
            if (ctx->transmarray[k] > tmax)
            {
                tmax = (float) ctx->transmarray[k];
            }
        }
        for (long k = 0; k < ctx->kmax; k++)
        {
            ctx->transmarray[k] /= tmax;
        }

        FILE *fp = fopen(transmname, "w");
        for (long k = 0; k < ctx->kmax; k++)
        {
            fprintf(fp, "%5.3f %15.13f\n", ctx->separray[k], ctx->transmarray[k]);
        }
        fclose(fp);
    }

    for (long k = 0; k < ctx->kmax; k++)
    {
        snprintf(fname1, sizeof(fname1), "psf%02ld", k);
        delete_image_ID(fname1, DELETE_IMAGE_ERRMODE_WARNING);
    }
}

/**
 * coro_scan_compute_grid - Populate PSF, stellar halo, and throughput profiles across scan grid
 * @ctx: Initialized PIAACMC scan context
 *
 * Return: RETURN_SUCCESS.
 */
errno_t coro_scan_compute_grid(coro_scan_ctx_t *ctx)
{
    struct stat buffer;
    char        fnameprol[500];

    for (long fpm_idx = 0; fpm_idx < ctx->fpmaskradarraysize; fpm_idx++)
    {
        for (long obs1_idx = ctx->obs1_index_start; obs1_idx < ctx->obs1arraysize; obs1_idx++)
        {
            double obs1 = ctx->obs1array[obs1_idx];
            if (obs1 > 0.001)
            {
                ctx->cent_obs_flag = 1;
            }
            double fpmaskrad = ctx->fpmaskradarray[fpm_idx];

            if ((fpmaskrad <= 0.8) || (fpmaskrad >= 2.0) || (obs1 <= 0.0) || (obs1 >= 0.4))
            {
                continue;
            }

            scan_write_aplc_conf(fpmaskrad, ctx->obs0, obs1);
            if (ctx->cent_obs_flag == 0)
            {
                snprintf(fnameprol,
                         sizeof(fnameprol),
                         "%s/APLCapo/raw/APLCapo_%5.3f.4096.ref.gz",
                         CORONAGRAPHSDATALOCAL,
                         fpmaskrad);
            }
            else
            {
                snprintf(fnameprol,
                         sizeof(fnameprol),
                         "%s/APLCapo/raw/APLCapo_%5.3f.%5.3f.4096.ref.gz",
                         CORONAGRAPHSDATALOCAL,
                         fpmaskrad,
                         obs1);
            }

            if (stat(fnameprol, &buffer) != 0)
            {
                printf("--------- FILE %s MISSING -> NO PROCESSING ----\n", fnameprol);
                continue;
            }

            printf("FILE %s EXISTS... MOVING AHEAD\n", fnameprol);
            if (scan_ensure_cell_psfs(ctx, fpmaskrad, obs1) == 1)
            {
                scan_write_star_profile(ctx, fpmaskrad, obs1);
                scan_write_transm_and_cleanup(ctx, fpmaskrad, obs1);
            }
        }
    }

    return RETURN_SUCCESS;
}
