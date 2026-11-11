/**
 * @file    AOsystSim_wf.c
 * @brief   Real-time atmospheric wavefront interpolation, binning, and stream publishing loop
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

#include "AOsystSim/AOsystSim.h"
#include "AOsystSim_internal.h"

/**
 * @brief Per-iteration temporal interpolation and spatial binning state.
 */
struct mkwf_interp_state
{
    imageID ID0;
    imageID ID1;
    imageID ID0amp;
    imageID ID1amp;
    long    kk0;
    long    kk1;
    float   alpha;
    long    csize;
};

/**
 * mkwf_ensure_cube - Lazy-load phase and amplitude WF FITS cubes for span index idx
 * @ctx:   Wavefront generator runtime context
 * @idx:   Span file index
 * @IDpha: Output phase cube image ID
 * @IDamp: Output amplitude cube image ID
 */
static void mkwf_ensure_cube(
    AOSIM_MKWF_CTX *ctx,
    long            idx,
    imageID        *IDpha,
    imageID        *IDamp)
{
    char wfimname_pha[STRINGMAXLEN_STREAMNAME];
    char wfimname_amp[STRINGMAXLEN_STREAMNAME];
    char wf_fname[STRINGMAXLEN_FULLFILENAME];
    long lam_pm = (long) (1.0e12 * ctx->LAMBDA + 0.5);

    sprintf(wfimname_pha, "wf%08ld_pha", idx);
    sprintf(wfimname_amp, "wf%08ld_amp", idx);

    if (image_ID(wfimname_pha, dcimg, dcnimg) == -1)
    {
        SNPRINTF_CHECK(wf_fname,
                       STRINGMAXLEN_FULLFILENAME,
                       "%s/%s%08ld.%09ld.pha.fits",
                       ctx->WFDIR,
                       ctx->wfin_PREFIX,
                       idx,
                       lam_pm);
        printf("Loading WF file name : %s\n", wf_fname);
        load_fits(wf_fname, wfimname_pha, 1, IDpha);

        SNPRINTF_CHECK(wf_fname,
                       STRINGMAXLEN_FULLFILENAME,
                       "%s/%s%08ld.%09ld.amp.fits",
                       ctx->WFDIR,
                       ctx->wfin_PREFIX,
                       idx,
                       lam_pm);
        printf("Loading WF file name : %s\n", wf_fname);
        load_fits(wf_fname, wfimname_amp, 1, IDamp);
        if (*IDamp == -1)
        {
            ctx->AMPfile = 0;
        }
    }

    *IDpha = image_ID(wfimname_pha, dcimg, dcnimg);
    if (ctx->AMPfile == 1)
    {
        *IDamp = image_ID(wfimname_amp, dcimg, dcnimg);
    }
}

/**
 * mkwf_resolve_cubes - Advance span indices for current physical time and load bracketing WF cubes
 * @ctx: Wavefront generator runtime context
 * @k:   Current circular WF file index (updated in place)
 * @k1:  Cumulative span counter (updated in place)
 * @st:  Output interpolation state
 */
static void mkwf_resolve_cubes(
    AOSIM_MKWF_CTX           *ctx,
    long                     *k,
    long                     *k1,
    struct mkwf_interp_state *st)
{
    float   frame_f = (float) ((ctx->phystime - ctx->wfin_TIME_SPAN * (*k1)) / ctx->wfin_TIME_STEP);
    long    frame_n = (long) frame_f;
    long    knext;
    imageID IDwf0    = -1;
    imageID IDwf1    = -1;
    imageID IDwf0amp = -1;
    imageID IDwf1amp = -1;

    st->alpha = frame_f - frame_n;
    if (frame_n > (ctx->NBframes - 1))
    {
        (*k)++;
        if (*k == ctx->kmax)
        {
            *k = 0;
        }
        (*k1)++;
        frame_n -= ctx->NBframes;
    }

    mkwf_ensure_cube(ctx, *k, &IDwf0, &IDwf0amp);
    knext = (*k == ctx->kmax) ? 0 : (*k + 1);
    mkwf_ensure_cube(ctx, knext, &IDwf1, &IDwf1amp);

    st->ID0    = IDwf0;
    st->ID0amp = IDwf0amp;
    st->kk0    = frame_n;
    if (frame_n + 1 == ctx->NBframes)
    {
        st->kk1    = 0;
        st->ID1    = IDwf1;
        st->ID1amp = IDwf1amp;
    }
    else
    {
        st->kk1    = frame_n + 1;
        st->ID1    = IDwf0;
        st->ID1amp = IDwf0amp;
    }
}

/**
 * mkwf_bin_step - Perform coherent complex-phasor spatial binning over active pupil subpixels
 * @ctx: Wavefront generator runtime context
 * @st:  Per-step interpolation state
 */
static void mkwf_bin_step(
    AOSIM_MKWF_CTX           *ctx,
    struct mkwf_interp_state *st)
{
    long csize    = ctx->wfin_WFsize / ctx->PIXBINFACTOR;
    long arr      = ctx->ARRAYSIZE;
    long wfsize   = ctx->wfin_WFsize;
    long wf2      = wfsize * wfsize;
    long off1     = (csize - arr) / 2;
    long ii1start = (off1 + ctx->ii1start0 < 0) ? 0 : (off1 + ctx->ii1start0);
    long ii1end   = (off1 + ctx->ii1end0 > csize) ? csize : (off1 + ctx->ii1end0);
    long jj1start = (off1 + ctx->jj1start0 < 0) ? 0 : (off1 + ctx->jj1start0);
    long jj1end   = (off1 + ctx->jj1end0 > csize) ? csize : (off1 + ctx->jj1end0);
    int  bf       = ctx->PIXBINFACTOR;

    st->csize = csize;
    if ((ctx->PIXSCALEMODE != 2) || (bf == 1))
    {
        return;
    }

    if (image_ID("tmpre", dcimg, dcnimg) == -1)
    {
        create_2Dimage_ID("tmpre", csize, csize, &ctx->IDbin_re);
        create_2Dimage_ID("tmpim", csize, csize, &ctx->IDbin_im);
        create_2Dimage_ID("tmpamp", csize, csize, &ctx->IDbin_amp);
        create_2Dimage_ID("tmpopd", csize, csize, &ctx->IDbin_opd);
    }

    for (long ii1 = ii1start; ii1 < ii1end; ii1++)
    {
        for (long jj1 = jj1start; jj1 < jj1end; jj1++)
        {
            if (dcimg[ctx->IDampmask].array.F[(jj1 - off1) * arr + (ii1 - off1)] <= 1.0e-6f)
            {
                continue;
            }
            dcimg[ctx->IDbin_re].array.F[jj1 * csize + ii1]  = 0.0f;
            dcimg[ctx->IDbin_im].array.F[jj1 * csize + ii1]  = 0.0f;
            dcimg[ctx->IDbin_opd].array.F[jj1 * csize + ii1] = 0.0f;

            for (int i = 0; i < bf; i++)
            {
                for (int j = 0; j < bf; j++)
                {
                    long   pix = (bf * jj1 + j) * wfsize + (bf * ii1 + i);
                    double pha = (1.0 - st->alpha) * dcimg[st->ID0].array.F[st->kk0 * wf2 + pix] +
                                 st->alpha * dcimg[st->ID1].array.F[st->kk1 * wf2 + pix];
                    float  opd = (float) (ctx->OPDMFACT * pha / (2.0 * M_PI) * ctx->LAMBDA * 1.0e6);
                    float  amp = 1.0f;
                    double phi = 1.0e-6 * opd / ctx->LAMBDA * 2.0 * M_PI;

                    if (ctx->AMPfile == 1)
                    {
                        amp = (1.0f - st->alpha) * dcimg[st->ID0amp].array.F[st->kk0 * wf2 + pix] +
                              st->alpha * dcimg[st->ID1amp].array.F[st->kk1 * wf2 + pix];
                    }
                    dcimg[ctx->IDbin_re].array.F[jj1 * csize + ii1] += amp * (float) cos(phi);
                    dcimg[ctx->IDbin_im].array.F[jj1 * csize + ii1] += amp * (float) sin(phi);
                    dcimg[ctx->IDbin_opd].array.F[jj1 * csize + ii1] += opd;
                }
            }
            {
                float re    = dcimg[ctx->IDbin_re].array.F[jj1 * csize + ii1];
                float im    = dcimg[ctx->IDbin_im].array.F[jj1 * csize + ii1];
                float norm2 = (float) (bf * bf);

                dcimg[ctx->IDbin_amp].array.F[jj1 * csize + ii1] = sqrtf(re * re + im * im) / norm2;
                dcimg[ctx->IDbin_opd].array.F[jj1 * csize + ii1] /= norm2;
            }
        }
    } // for (long ii1)
}

/**
 * mkwf_write_out0 - Populate OUT0 OPD/AMP wavefront streams and optional FITS snapshots
 * @ctx: Wavefront generator runtime context
 * @st:  Per-step interpolation state
 */
static void mkwf_write_out0(
    AOSIM_MKWF_CTX                 *ctx,
    const struct mkwf_interp_state *st)
{
    long arr     = ctx->ARRAYSIZE;
    long wfsize  = ctx->wfin_WFsize;
    long wf2     = wfsize * wfsize;
    long iisize  = wfsize / ctx->PIXBINFACTOR;
    long iistart = (arr / 2 - iisize / 2 < 0) ? 0 : (arr / 2 - iisize / 2);
    long iiend   = (arr / 2 + iisize / 2 > arr) ? arr : (arr / 2 + iisize / 2);
    char fname[STRINGMAXLEN_FULLFILENAME];

    dcimg[ctx->IDopd0].md[0].write = 1;
    for (long ii = iistart; ii < iiend; ii++)
    {
        for (long jj = iistart; jj < iiend; jj++)
        {
            if (ctx->PIXBINFACTOR == 1)
            {
                long pix = (jj + ctx->jjoffset) * wfsize + (ii + ctx->iioffset);

                dcimg[ctx->IDopd0].array.F[jj * arr + ii] =
                    ctx->OPDMFACT *
                    ((1.0f - st->alpha) * dcimg[st->ID0].array.F[st->kk0 * wf2 + pix] +
                     st->alpha * dcimg[st->ID1].array.F[st->kk1 * wf2 + pix]);
            }
            else
            {
                long ii1 = (ii - iistart) + ctx->iioffset;
                long jj1 = (jj - iistart) + ctx->jjoffset;

                dcimg[ctx->IDopd0].array.F[jj * arr + ii] =
                    dcimg[ctx->IDbin_opd].array.F[jj1 * st->csize + ii1];
            }
        }
    }
    dcimg[ctx->IDopd0].md[0].cnt0++;
    dcimg[ctx->IDopd0].md[0].write = 0;
    COREMOD_MEMORY_image_set_sempost(ctx->OUT0STREAMNAMEOPD, -1);

    if (ctx->OUT0STREAM == 2)
    {
        dcimg[ctx->IDamp0].md[0].write = 1;
        for (long ii = iistart; ii < iiend; ii++)
        {
            for (long jj = iistart; jj < iiend; jj++)
            {
                float val = ctx->AMPMFACT * dcimg[ctx->IDampmask].array.F[jj * arr + ii];

                if (ctx->PIXBINFACTOR == 1)
                {
                    if (ctx->AMPfile == 1)
                    {
                        long pix = (jj + ctx->jjoffset) * wfsize + (ii + ctx->iioffset);

                        val *= (1.0f - st->alpha) * dcimg[st->ID0amp].array.F[st->kk0 * wf2 + pix] +
                               st->alpha * dcimg[st->ID1amp].array.F[st->kk1 * wf2 + pix];
                    }
                }
                else if (ctx->AMPfile == 1)
                {
                    long ii1 = ii + ctx->iioffset;
                    long jj1 = jj + ctx->jjoffset;
                    long pix = jj1 * wfsize + ii1;

                    val *= dcimg[ctx->IDbin_amp].array.F[jj1 * st->csize + ii1] *
                           ((1.0f - st->alpha) * dcimg[st->ID0amp].array.F[st->kk0 * wf2 + pix] +
                            st->alpha * dcimg[st->ID1amp].array.F[st->kk1 * wf2 + pix]);
                }
                else
                {
                    long ii1 = (ii - iistart) + ctx->iioffset;
                    long jj1 = (jj - iistart) + ctx->jjoffset;

                    val *= dcimg[ctx->IDbin_amp].array.F[jj1 * st->csize + ii1];
                }
                dcimg[ctx->IDamp0].array.F[jj * arr + ii] = val;
            }
        }
        dcimg[ctx->IDamp0].md[0].cnt0++;
        dcimg[ctx->IDamp0].md[0].write = 0;
        COREMOD_MEMORY_image_set_sempost(ctx->OUT0STREAMNAMEAMP, -1);
    } // if (ctx->OUT0STREAM == 2)

    if (ctx->OUT0FITSFILE > 0)
    {
        EXECUTE_SYSTEM_COMMAND("rm %s", ctx->OUT0FITSFILENAMEOPD);
        SNPRINTF_CHECK(fname, STRINGMAXLEN_FULLFILENAME, "%s", ctx->OUT0FITSFILENAMEOPD);
        save_fits(ctx->OUT0STREAMNAMEOPD, fname);
    }
    if ((ctx->OUT0FITSFILE == 2) && (ctx->OUT0STREAM == 2))
    {
        EXECUTE_SYSTEM_COMMAND("rm %s", ctx->OUT0FITSFILENAMEAMP);
        SNPRINTF_CHECK(fname, STRINGMAXLEN_FULLFILENAME, "%s", ctx->OUT0FITSFILENAMEAMP);
        save_fits(ctx->OUT0STREAMNAMEAMP, fname);
    }
}

/**
 * mkwf_postproc_and_wait - Apply DM0 reflection OPD, update physical clock, and wait on trigger
 * @ctx: Wavefront generator runtime context
 */
static void mkwf_postproc_and_wait(AOSIM_MKWF_CTX *ctx)
{
    uint64_t npix = (uint64_t) ctx->ARRAYSIZE * ctx->ARRAYSIZE;
    imageID  ID;
    int      OKf;

    if (ctx->DM0MODE > 0)
    {
        dcimg[ctx->IDopd1].md[0].write = 1;
        dcimg[ctx->IDamp1].md[0].write = 1;
        for (uint64_t ii = 0; ii < npix; ii++)
        {
            dcimg[ctx->IDopd1].array.F[ii] =
                dcimg[ctx->IDopd0].array.F[ii] - 2.0f * dcimg[ctx->IDdm0opd].array.F[ii];
        }
        memcpy(dcimg[ctx->IDamp1].array.F, dcimg[ctx->IDamp0].array.F, sizeof(float) * npix);
        dcimg[ctx->IDopd1].md[0].cnt0++;
        dcimg[ctx->IDamp1].md[0].cnt0++;
        dcimg[ctx->IDopd1].md[0].write = 0;
        dcimg[ctx->IDamp1].md[0].write = 0;
        COREMOD_MEMORY_image_set_sempost(ctx->OUT1STREAMNAMEOPD, -1);
        COREMOD_MEMORY_image_set_sempost(ctx->OUT1STREAMNAMEAMP, -1);
    }

    ctx->phystime += ctx->DT;
    dcimg[ctx->IDphystime].md[0].write = 0;
    dcimg[ctx->IDphystime].array.F[0]  = (float) ctx->phystime;
    COREMOD_MEMORY_image_set_sempost_byID(ctx->IDphystime, -1);
    dcimg[ctx->IDphystime].md[0].cnt0++;
    dcimg[ctx->IDphystime].md[0].write = 0;

    switch (ctx->TRIGGERMODE)
    {
        case 0:
            OKf = 0;
            while (OKf == 0)
            {
                usleep((long) (1.0e6 * ctx->TRIGGERDT));
                if (file_exists(ctx->TRIGGERFILE) == 1)
                {
                    EXECUTE_SYSTEM_COMMAND("rm %s", ctx->TRIGGERFILE);
                    OKf = 1;
                }
            }
            break;

        case 1:
            ID = image_ID(ctx->TRIGGER0STREAM, dcimg, dcnimg);
            while (ID == -1)
            {
                usleep(1000000);
                ID = read_sharedmem_image(ctx->TRIGGER0STREAM, dcimg, dcnimg);
            }
            COREMOD_MEMORY_image_set_semwait(ctx->TRIGGER0STREAM, ctx->TRIGGER0SEM);

            ID = image_ID(ctx->TRIGGER1STREAM, dcimg, dcnimg);
            while (ID == -1)
            {
                usleep(1000000);
                ID = read_sharedmem_image(ctx->TRIGGER1STREAM, dcimg, dcnimg);
            }
            COREMOD_MEMORY_image_set_semwait(ctx->TRIGGER1STREAM, ctx->TRIGGER1SEM);
            break;

        case 2:
            usleep((long) (1.0e6 * ctx->TRIGGERDT));
            break;

        default:
            printf("TRIGGERMODE = %d not valid\n", ctx->TRIGGERMODE);
            exit(0);
            break;
    } // switch (ctx->TRIGGERMODE)
}

/**
 * AOsystSim_mkWF - Stream multi-layer Fresnel-propagated atmospheric wavefronts from config file
 * @CONF_FNAME: Atmospheric wavefront daemon configuration file path
 *
 * Return: 0 on completion.
 */
int AOsystSim_mkWF(const char *CONF_FNAME)
{
    AOSIM_MKWF_CTX ctx = { 0 };
    long           k   = 0;
    long           k1  = 0;

    printf("AOsystSim mkWFS...\n");

    aosim_mkwf_read_conf(CONF_FNAME, &ctx);
    aosim_mkwf_create_streams(&ctx);
    aosim_mkwf_init_pupil(&ctx);

    for (;;)
    {
        struct mkwf_interp_state st = { 0 };

        mkwf_resolve_cubes(&ctx, &k, &k1, &st);
        mkwf_bin_step(&ctx, &st);
        mkwf_write_out0(&ctx, &st);
        mkwf_postproc_and_wait(&ctx);

        if (k > ctx.kmax)
        {
            k = 0;
        }
    }

    list_image_ID();
    return 0;
}
