/**
 * @file    AOsystSim_lowfs.c
 * @brief   Real-time coronagraph and Lyot low-order wavefront sensor (LOWFS) propagation loop
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "fft/fft.h"

#include "AOsystSim/AOsystSim.h"
#include "AOsystSim_internal.h"

/**
 * lowfs_fft_amph - Perform centered 2D FFT propagation between amplitude/phase image pairs
 * @in_amp:  Input amplitude image name
 * @in_pha:  Input phase image name
 * @tmp_c1:  Scratch input complex image name
 * @tmp_c2:  Scratch transformed complex image name
 * @out_amp: Output shared-memory amplitude image name
 * @out_pha: Output shared-memory phase image name
 */
static void lowfs_fft_amph(
    const char *in_amp,
    const char *in_pha,
    const char *tmp_c1,
    const char *tmp_c2,
    const char *out_amp,
    const char *out_pha)
{
    mk_complex_from_amph(in_amp, in_pha, tmp_c1, 0);
    permut(tmp_c1);
    do2dfft(tmp_c1, tmp_c2);
    permut(tmp_c2);
    mk_amph_from_complex(tmp_c2, out_amp, out_pha, 1);
    delete_image_ID(tmp_c1, DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID(tmp_c2, DELETE_IMAGE_ERRMODE_WARNING);
}

/**
 * lowfs_publish_pair - Increment counters and post semaphores on an updated amp/phase stream pair
 * @IDa: Amplitude stream image ID
 * @IDp: Phase stream image ID
 */
static void lowfs_publish_pair(
    imageID IDa,
    imageID IDp)
{
    COREMOD_MEMORY_image_set_sempost_byID(IDa, -1);
    COREMOD_MEMORY_image_set_sempost_byID(IDp, -1);
    dcimg[IDa].md[0].cnt0++;
    dcimg[IDp].md[0].cnt0++;
    dcimg[IDa].md[0].write = 0;
    dcimg[IDp].md[0].write = 0;
}

/**
 * lowfs_propagate_coro - Propagate entrance WF through FPM and split Lyot stop arms
 * @ctx: Coronagraph/LOWFS runtime context
 */
static void lowfs_propagate_coro(AOSIM_LOWFS_CTX *ctx)
{
    imageID IDfoc0a;
    imageID IDfoc0p;
    imageID IDpup1a;
    imageID IDpup1p;

    for (uint32_t ii = 0; ii < ctx->xsizein; ii++)
    {
        for (uint32_t jj = 0; jj < ctx->ysizein; jj++)
        {
            long dst = ((long) jj + ctx->jjoffset) * ctx->ARRAYSIZE + ((long) ii + ctx->iioffset);
            long src = (long) jj * ctx->xsizein + ii;

            dcimg[ctx->IDwfa].array.F[dst] = dcimg[ctx->IDinAMP].array.F[src];
            dcimg[ctx->IDwfp].array.F[dst] =
                (float) (M_PI * dcimg[ctx->IDinOPD].array.F[src] * 1.0e-6 / ctx->LAMBDA);
        }
    }

    lowfs_fft_amph("aosim_wfa",
                   "aosim_wfp",
                   "loc_aosim_wfc",
                   "loc_aosim_fc0",
                   "aosim_foc0_amp",
                   "aosim_foc0_pha");

    IDfoc0a                        = image_ID("aosim_foc0_amp", dcimg, dcnimg);
    IDfoc0p                        = image_ID("aosim_foc0_pha", dcimg, dcnimg);
    dcimg[ctx->IDfoc1a].md[0].write = 1;
    dcimg[ctx->IDfoc1p].md[0].write = 1;
    for (uint64_t ii = 0; ii < ctx->ARRAYSIZE2; ii++)
    {
        dcimg[ctx->IDfoc1a].array.F[ii] =
            dcimg[IDfoc0a].array.F[ii] * dcimg[ctx->IDfpmamp].array.F[ii];
        dcimg[ctx->IDfoc1p].array.F[ii] =
            dcimg[IDfoc0p].array.F[ii] + dcimg[ctx->IDfpmpha].array.F[ii];
    }
    lowfs_publish_pair(ctx->IDfoc1a, ctx->IDfoc1p);

    lowfs_fft_amph("aosim_foc1_amp",
                   "aosim_foc1_pha",
                   "loc_aosim_foc1_c",
                   "loc_aosim_pup1_c",
                   "aosim_pup1_amp",
                   "aosim_pup1_pha");

    IDpup1a                         = image_ID("aosim_pup1_amp", dcimg, dcnimg);
    IDpup1p                         = image_ID("aosim_pup1_pha", dcimg, dcnimg);
    dcimg[ctx->IDpup1ta].md[0].write = 1;
    dcimg[ctx->IDpup1tp].md[0].write = 1;
    for (uint64_t ii = 0; ii < ctx->ARRAYSIZE2; ii++)
    {
        dcimg[ctx->IDpup1ta].array.F[ii] =
            dcimg[IDpup1a].array.F[ii] * dcimg[ctx->ID_LS_tamp].array.F[ii];
        dcimg[ctx->IDpup1tp].array.F[ii] =
            dcimg[IDpup1p].array.F[ii] + dcimg[ctx->ID_LS_tpha].array.F[ii];
    }
    lowfs_publish_pair(ctx->IDpup1ta, ctx->IDpup1tp);

    lowfs_fft_amph("aosim_pup1t_amp",
                   "aosim_pup1t_pha",
                   "loc_aosim_pup1t_c",
                   "loc_aosim_foc2_c",
                   "aosim_foc2_amp",
                   "aosim_foc2_pha");

    dcimg[ctx->IDpup1ra].md[0].write = 1;
    dcimg[ctx->IDpup1rp].md[0].write = 1;
    for (uint64_t ii = 0; ii < ctx->ARRAYSIZE2; ii++)
    {
        dcimg[ctx->IDpup1ra].array.F[ii] =
            dcimg[IDpup1a].array.F[ii] * dcimg[ctx->ID_LS_ramp].array.F[ii];
        dcimg[ctx->IDpup1rp].array.F[ii] =
            dcimg[IDpup1p].array.F[ii] + dcimg[ctx->ID_LS_rpha].array.F[ii] +
            (float) (2.0 * M_PI * dcimg[ctx->ID_lowfsopd].array.F[ii] * 1.0e-6 / ctx->LAMBDA);
    }
    lowfs_publish_pair(ctx->IDpup1ra, ctx->IDpup1rp);

    lowfs_fft_amph("aosim_pup1r_amp",
                   "aosim_pup1r_pha",
                   "loc_aosim_pup1r_c",
                   "loc_aosim_foc1r_c",
                   "aosim_foclowfs_amp",
                   "aosim_foclowfs_pha");
}

/**
 * lowfs_accumulate_and_publish - Crop reflected focal intensity and integrate LOWFS camera frames
 * @ctx:  Coronagraph/LOWFS runtime context
 * @iter: Loop iteration counter
 */
static void lowfs_accumulate_and_publish(
    AOSIM_LOWFS_CTX *ctx,
    long             iter)
{
    imageID  IDfoclowfsa = image_ID("aosim_foclowfs_amp", dcimg, dcnimg);
    uint32_t outsz       = ctx->OUTLOWFSARRAYSIZE;
    long     crop_off    = ((long) ctx->ARRAYSIZE - (long) outsz) / 2;
    float    tnow;
    char     fname[STRINGMAXLEN_FULLFILENAME];

    dcimg[ctx->IDimlowfs].md[0].write = 1;
    for (uint32_t ii = 0; ii < outsz; ii++)
    {
        for (uint32_t jj = 0; jj < outsz; jj++)
        {
            long  ii1 = (long) ii + crop_off;
            long  jj1 = (long) jj + crop_off;
            float amp = dcimg[IDfoclowfsa].array.F[jj1 * ctx->ARRAYSIZE + ii1];
            float val = amp * amp;

            dcimg[ctx->IDimlowfs].array.F[outsz * jj + ii] = val;
            dcimg[ctx->IDimcamlowfstmp].array.F[outsz * jj + ii] += val;
        }
    }
    COREMOD_MEMORY_image_set_sempost_byID(ctx->IDimlowfs, -1);
    dcimg[ctx->IDimlowfs].md[0].cnt0++;
    dcimg[ctx->IDimlowfs].md[0].write = 0;
    ctx->LOWFScamcnt++;

    tnow = dcimg[ctx->IDphystime].array.F[0];
    printf("%8ld  TIME: %16f  %16f   %16f/%16f\n",
           iter,
           tnow,
           ctx->LOWFScamstarttime,
           tnow - ctx->LOWFScamstarttime,
           ctx->LOWFSCAMETIME);
    fflush(stdout);

    if ((tnow - ctx->LOWFScamstarttime) > ctx->LOWFSCAMETIME)
    {
        dcimg[ctx->IDoutLOWFS].md[0].write      = 1;
        dcimg[ctx->IDimcamlowfstmp].md[0].write = 1;
        for (uint64_t ii = 0; ii < ctx->OUTLOWFSARRAYSIZE2; ii++)
        {
            dcimg[ctx->IDoutLOWFS].array.F[ii] =
                dcimg[ctx->IDimcamlowfstmp].array.F[ii] / (float) ctx->LOWFScamcnt;
            dcimg[ctx->IDimcamlowfstmp].array.F[ii] = 0.0f;
        }
        COREMOD_MEMORY_image_set_sempost_byID(ctx->IDoutLOWFS, -1);
        COREMOD_MEMORY_image_set_sempost_byID(ctx->IDimcamlowfstmp, -1);
        dcimg[ctx->IDoutLOWFS].md[0].cnt0++;
        dcimg[ctx->IDimcamlowfstmp].md[0].cnt0++;
        dcimg[ctx->IDoutLOWFS].md[0].write      = 0;
        dcimg[ctx->IDimcamlowfstmp].md[0].write = 0;
        ctx->LOWFScamcnt                        = 0;
        ctx->LOWFScamstarttime                  = tnow;
    }

    if (ctx->OUTMODE == 1)
    {
        EXECUTE_SYSTEM_COMMAND("rm %s", ctx->OUTLOWFSFITSFILENAME);
        SNPRINTF_CHECK(fname, STRINGMAXLEN_FULLFILENAME, "%s", ctx->OUTLOWFSFITSFILENAME);
        save_fits(ctx->OUTLOWFSSTREAMNAME, fname);
        EXECUTE_SYSTEM_COMMAND("touch %s", ctx->OUTTRIGGERFILE);
    }
}

/**
 * lowfs_wait_trigger - Wait for next wavefront OPD trigger semaphore or trigger file
 * @ctx: Coronagraph/LOWFS runtime context
 */
static void lowfs_wait_trigger(const AOSIM_LOWFS_CTX *ctx)
{
    int OKf;

    switch (ctx->INMODE)
    {
        case 0:
            COREMOD_MEMORY_image_set_semwait(ctx->INTRIGSTREAMNAME, ctx->INTRIGSEMCHAN);
            break;

        case 1:
            OKf = 0;
            while (OKf == 0)
            {
                usleep(10000);
                if (file_exists(ctx->INTRIGGERFILE) == 1)
                {
                    EXECUTE_SYSTEM_COMMAND("rm %s", ctx->INTRIGGERFILE);
                    OKf = 1;
                }
            }
            break;

        default:
            printf("INMODE value %d not valid\n", ctx->INMODE);
            exit(0);
            break;
    }
}

/**
 * AOsystSim_coroLOWFS - Run real-time shared-memory coronagraph + Lyot LOWFS simulation loop
 * @CONF_FNAME: Coronagraph/LOWFS daemon configuration file path
 *
 * Return: 0 on completion.
 */
int AOsystSim_coroLOWFS(const char *CONF_FNAME)
{
    AOSIM_LOWFS_CTX ctx  = { 0 };
    long            kmax = 100000000;

    printf("AOsystSim coro LOWFS...\n");

    aosim_lowfs_read_conf(CONF_FNAME, &ctx);
    aosim_lowfs_connect_inputs(&ctx);
    aosim_lowfs_init_buffers(&ctx);

    for (long k = 0; k < kmax; k++)
    {
        lowfs_propagate_coro(&ctx);
        lowfs_accumulate_and_publish(&ctx, k);
        lowfs_wait_trigger(&ctx);
    }

    return 0;
}
