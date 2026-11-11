/**
 * @file    AOsystSim_filter.c
 * @brief   Open-loop spatial-frequency and time-lag AO transfer filter stream simulator
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "fft/fft.h"

#include "AOsystSim/AOsystSim.h"
#include "AOsystSim_internal.h"

/**
 * @brief Internal state for the Fourier-domain AO spatial/temporal filter loop.
 */
struct aosim_filt_ctx
{
    uint32_t sizearray[2];
    uint64_t size2;
    imageID  IDin;
    imageID  IDin1;
    imageID  IDmask;
    imageID  IDout;
    imageID  IDdm;
    imageID  IDwfe;
    imageID  IDaosf_noise;
    imageID  IDaosf_mult;
    imageID  IDaosf_gain;
    imageID  IDwfcbuff;
    long     wfcbuff_size;
    double  *wfcbuff_time;
    long    *wfcbuff_status;
    float    loopgain;
    float    multr0;
    float    wfsetime;
    float    timedelay;
    double   dmmovetime;
    long     dmmoveNBpt;
    double   tnowdouble;
    double   time_wfse0;
    long     wfsecnt;
    uint64_t cnt0;
};

/**
 * aosim_filt_init_streams - Attach input turbulence stream and allocate delay ring buffer
 * @ctx:       Filter loop context
 * @IDin_name: Input wavefront shared-memory stream name
 */
static void aosim_filt_init_streams(
    struct aosim_filt_ctx *ctx,
    const char            *IDin_name)
{
    ctx->IDin         = read_sharedmem_image(IDin_name, dcimg, dcnimg);
    ctx->sizearray[0] = dcimg[ctx->IDin].md[0].size[0];
    ctx->sizearray[1] = dcimg[ctx->IDin].md[0].size[1];
    ctx->size2        = (uint64_t) ctx->sizearray[0] * ctx->sizearray[1];

    create_3Dimage_ID("wfcbuff",
                      ctx->sizearray[0],
                      ctx->sizearray[1],
                      ctx->wfcbuff_size,
                      &ctx->IDwfcbuff);

    ctx->wfcbuff_time = (double *) malloc(sizeof(double) * ctx->wfcbuff_size);
    if (ctx->wfcbuff_time == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    ctx->wfcbuff_status = (long *) malloc(sizeof(long) * ctx->wfcbuff_size);
    if (ctx->wfcbuff_status == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    for (long kk = 0; kk < ctx->wfcbuff_size; kk++)
    {
        ctx->wfcbuff_status[kk] = 0;
    }

    create_2Dimage_ID("aosf_noise", ctx->sizearray[0], ctx->sizearray[1], &ctx->IDaosf_noise);
    create_2Dimage_ID("aosf_mult", ctx->sizearray[0], ctx->sizearray[1], &ctx->IDaosf_mult);
    create_2Dimage_ID("aosf_gain", ctx->sizearray[0], ctx->sizearray[1], &ctx->IDaosf_gain);
    create_2Dimage_ID("aosf_mask", ctx->sizearray[0], ctx->sizearray[1], &ctx->IDmask);
}

/**
 * aosim_filt_init_masks - Build pupil apodization mask, Fourier transfer filters, and output stream
 * @ctx:        Filter loop context
 * @IDout_name: Output filtered wavefront shared-memory stream name
 */
static void aosim_filt_init_masks(
    struct aosim_filt_ctx *ctx,
    const char            *IDout_name)
{
    uint32_t nx    = ctx->sizearray[0];
    uint32_t ny    = ctx->sizearray[1];
    double   rmask = 0.4 * nx;

    aosim_ensure_wdir();

    for (uint32_t ii = 0; ii < nx; ii++)
    {
        for (uint32_t jj = 0; jj < ny; jj++)
        {
            double   x   = 1.0 * ii - 0.5 * nx;
            double   y   = 1.0 * jj - 0.5 * ny;
            double   r   = sqrt(x * x + y * y);
            double   r1  = (r - rmask) / (0.5 * nx - rmask);
            uint64_t idx = (uint64_t) jj * nx + ii;

            if (r1 < 0.0)
            {
                dcimg[ctx->IDmask].array.F[idx] = 1.0f;
            }
            else if (r1 > 1.0)
            {
                dcimg[ctx->IDmask].array.F[idx] = 0.0f;
            }
            else
            {
                dcimg[ctx->IDmask].array.F[idx] = (float) (0.5 * (cos(r1 * M_PI) + 1.0));
            }

            dcimg[ctx->IDaosf_gain].array.F[idx] = ctx->loopgain;
            dcimg[ctx->IDaosf_mult].array.F[idx] =
                (float) exp(-pow(r / (ctx->multr0 * nx), 8.0));
            if (r > ctx->multr0 * nx)
            {
                dcimg[ctx->IDaosf_mult].array.F[idx] = 0.0f;
            }
            dcimg[ctx->IDaosf_noise].array.F[idx] = 0.0f;
        }
    }

    save_fits("aosf_noise", "AOsystSim_wdir/aosf_noise.fits");
    save_fits("aosf_mult", "AOsystSim_wdir/aosf_mult.fits");
    save_fits("aosf_gain", "AOsystSim_wdir/aosf_gain.fits");

    permut("aosf_mult");
    permut("aosf_noise");
    permut("aosf_gain");

    create_image_ID("aofiltdm", 2, ctx->sizearray, _DATATYPE_FLOAT, 1, 0, 0, &ctx->IDdm);
    create_image_ID("aofiltwfe", 2, ctx->sizearray, _DATATYPE_FLOAT, 1, 0, 0, &ctx->IDwfe);
    create_image_ID(IDout_name, 2, ctx->sizearray, _DATATYPE_FLOAT, 1, 0, 0, &ctx->IDout);

    strcpy(dcimg[ctx->IDout].kw[0].name, "TIME");
    dcimg[ctx->IDout].kw[0].type       = 'D';
    dcimg[ctx->IDout].kw[0].value.numf = 0.0;
    strcpy(dcimg[ctx->IDout].kw[0].comment, "Physical time [sec]");

    create_image_ID("aofiltin", 2, ctx->sizearray, _DATATYPE_FLOAT, 1, 0, 0, &ctx->IDin1);
}

/**
 * aosim_filt_measure_wfe - Subtract DM OPD, measure residual RMS, and filter WFE in Fourier space
 * @ctx:        Filter loop context
 * @IDout_name: Output stream name used for FFT input
 */
static void aosim_filt_measure_wfe(
    struct aosim_filt_ctx *ctx,
    const char            *IDout_name)
{
    double  WFrms0   = 0.0;
    double  WFrms    = 0.0;
    long    WFrmscnt = 0;
    imageID ID;

    for (uint64_t ii = 0; ii < ctx->size2; ii++)
    {
        dcimg[ctx->IDin1].array.F[ii] =
            dcimg[ctx->IDin].array.F[ii] * dcimg[ctx->IDmask].array.F[ii];
    }

    for (uint64_t ii = 0; ii < ctx->size2; ii++)
    {
        float in_val  = dcimg[ctx->IDin1].array.F[ii];
        float out_val = in_val - dcimg[ctx->IDdm].array.F[ii];

        dcimg[ctx->IDout].array.F[ii] = out_val;
        if (dcimg[ctx->IDmask].array.F[ii] > 0.99f)
        {
            WFrms0 += (double) in_val * in_val;
            WFrms += (double) out_val * out_val;
            WFrmscnt++;
        }
    }
    WFrms0                             = sqrt(WFrms0 / WFrmscnt);
    WFrms                              = sqrt(WFrms / WFrmscnt);
    dcimg[ctx->IDout].kw[0].value.numf = ctx->tnowdouble;

    ctx->tnowdouble = dcimg[ctx->IDin].kw[0].value.numf;
    printf("\r Time : %.6lf    WFSexposure time:  %.6f   [%5ld]    %10.5lf -> %10.5lf ",
           ctx->tnowdouble,
           ctx->tnowdouble - ctx->time_wfse0,
           ctx->wfsecnt,
           WFrms0,
           WFrms);
    fflush(stdout);
    ctx->cnt0 = dcimg[ctx->IDin].md[0].cnt0;

    do2drfft(IDout_name, "aosf_tmpfft");
    ID = image_ID("aosf_tmpfft", dcimg, dcnimg);
    for (uint64_t ii = 0; ii < ctx->size2; ii++)
    {
        float scale = dcimg[ctx->IDaosf_mult].array.F[ii] / (float) ctx->size2;

        dcimg[ID].array.CF[ii].re *= scale;
        dcimg[ID].array.CF[ii].im *= scale;
    }

    do2dffti("aosf_tmpfft", "testo");
    delete_image_ID("aosf_tmpfft", DELETE_IMAGE_ERRMODE_WARNING);
    ID = image_ID("testo", dcimg, dcnimg);

    for (uint64_t ii = 0; ii < ctx->size2; ii++)
    {
        dcimg[ctx->IDwfe].array.F[ii] += dcimg[ID].array.CF[ii].re;
    }
    delete_image_ID("testo", DELETE_IMAGE_ERRMODE_WARNING);
    ctx->wfsecnt++;
}

/**
 * aosim_filt_update_dm - Enqueue completed WFS exposure into delay buffer and apply DM step
 * @ctx: Filter loop context
 */
static void aosim_filt_update_dm(struct aosim_filt_ctx *ctx)
{
    if ((ctx->tnowdouble - ctx->time_wfse0) > ctx->wfsetime)
    {
        if (ctx->wfsecnt > 0)
        {
            for (uint64_t ii = 0; ii < ctx->size2; ii++)
            {
                dcimg[ctx->IDwfe].array.F[ii] /= (float) ctx->wfsecnt;
            }
        }

        for (long k1 = 0; k1 < ctx->dmmoveNBpt; k1++)
        {
            long k0 = 0;

            while ((k0 < ctx->wfcbuff_size) && (ctx->wfcbuff_status[k0] == 1))
            {
                k0++;
            }
            if (k0 > ctx->wfcbuff_size - 1)
            {
                printf("\n WFC buffer full [%ld]\n", k0);
                exit(0);
            }

            for (uint64_t ii = 0; ii < ctx->size2; ii++)
            {
                dcimg[ctx->IDwfcbuff].array.F[(uint64_t) k0 * ctx->size2 + ii] =
                    dcimg[ctx->IDwfe].array.F[ii] / (float) ctx->dmmoveNBpt;
            }
            ctx->wfcbuff_time[k0] = ctx->tnowdouble + ctx->timedelay +
                                    1.0 * k1 / (ctx->dmmoveNBpt - 1) * ctx->dmmovetime;
            ctx->wfcbuff_status[k0] = 1;
        }

        for (uint64_t ii = 0; ii < ctx->size2; ii++)
        {
            dcimg[ctx->IDwfe].array.F[ii] = 0.0f;
        }
        ctx->time_wfse0 += ctx->wfsetime;
        ctx->wfsecnt = 0;
    }

    for (long kk = 0; kk < ctx->wfcbuff_size; kk++)
    {
        if ((ctx->wfcbuff_status[kk] == 1) && (ctx->wfcbuff_time[kk] < ctx->tnowdouble))
        {
            for (uint64_t ii = 0; ii < ctx->size2; ii++)
            {
                dcimg[ctx->IDdm].array.F[ii] +=
                    dcimg[ctx->IDaosf_gain].array.F[ii] *
                    dcimg[ctx->IDwfcbuff].array.F[(uint64_t) kk * ctx->size2 + ii];
            }
            ctx->wfcbuff_status[kk] = 0;
        }
    }
    dcimg[ctx->IDout].md[0].cnt0 = ctx->cnt0;
}

/**
 * AOsystSim_simpleAOfilter - Apply open-loop spatial-frequency and time-delay AO transfer filter
 * @IDin_name:  Input atmospheric wavefront stream name
 * @IDout_name: Output AO-filtered wavefront stream name
 *
 * Return: 0 on completion.
 */
int AOsystSim_simpleAOfilter(
    const char *IDin_name,
    const char *IDout_name)
{
    struct aosim_filt_ctx ctx = { 0 };

    ctx.loopgain     = 0.001f;
    ctx.multr0       = 0.2f;
    ctx.wfsetime     = 0.001f;
    ctx.timedelay    = 0.0005f;
    ctx.wfcbuff_size = 10000;
    ctx.dmmovetime   = 0.001;
    ctx.dmmoveNBpt   = 10;

    aosim_filt_init_streams(&ctx, IDin_name);
    aosim_filt_init_masks(&ctx, IDout_name);

    printf("%s -> %s\n", IDin_name, IDout_name);
    ctx.cnt0       = (uint64_t) -1;
    ctx.time_wfse0 = dcimg[ctx.IDin].kw[0].value.numf;
    ctx.wfsecnt    = 0;

    for (;;)
    {
        usleep(10);
        if (dcimg[ctx.IDin].md[0].cnt0 != ctx.cnt0)
        {
            aosim_filt_measure_wfe(&ctx, IDout_name);
            aosim_filt_update_dm(&ctx);
        }
    }

    delete_image_ID("circbuff", DELETE_IMAGE_ERRMODE_WARNING);
    free(ctx.wfcbuff_time);
    free(ctx.wfcbuff_status);

    return 0;
}
