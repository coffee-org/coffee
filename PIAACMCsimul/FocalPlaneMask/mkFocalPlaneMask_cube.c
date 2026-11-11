/**
 * @file    mkFocalPlaneMask_cube.c
 * @brief   Subpixel-integrated polychromatic FPM complex-amplitude cube evaluation
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "OpticsMaterials/OpticsMaterials.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul.h"

#include "FocalPlaneMask/mkFocalPlaneMask_internal.h"

struct fpm_subpix_state
{
    double retmp;
    double imtmp;
    double ttmp;
    double zonetmp;
    double t;
    float  cospha;
    float  sinpha;
    double amp;
    long   zi;
};

/**
 * mkFocalPlaneMask_lookup_zone - Resolve central cone and discrete zone properties at (ii1, jj1)
 * @ctx: Active FPM synthesis context.
 * @k:   Wavelength slice index.
 * @r:   Radial subpixel distance from mask center [m].
 * @ii1: Horizontal pixel index in high-resolution zone map.
 * @jj1: Vertical pixel index in high-resolution zone map.
 * @st:  Subpixel state updated with sag, amplitude, and phase phasor.
 */
static void mkFocalPlaneMask_lookup_zone(
    const struct fpm_mask_ctx *ctx,
    int                        k,
    double                     r,
    long                       ii1,
    long                       jj1,
    struct fpm_subpix_state   *st)
{
    if (ctx->CentCone == 1 && r < 0.99 * piaacmcopticaldesign.fpmRad
        && r < piaacmcopticaldesign.fpmCentConeRad)
    {
        double sag_mid = 0.5 * (piaacmcopticaldesign.fpmminsag + piaacmcopticaldesign.fpmmaxsag);
        st->t          = piaacmcopticaldesign.fpmCentConeZ
                         + (r / piaacmcopticaldesign.fpmCentConeRad)
                               * (sag_mid - piaacmcopticaldesign.fpmCentConeZ);
        double pha = OpticsMaterials_pha_lambda(
            piaacmcopticaldesign.fpmmaterial_code,
            st->t,
            piaacmcopticalsystem.lambdaarray[k]);
        st->cospha = cosf((float) pha);
        st->sinpha = sinf((float) pha);
        st->amp    = 1.0;
    }

    long arrsize = piaacmcopticaldesign.fpmarraysize;
    st->zi       = (long) dcimg[ctx->IDz].array.UI16[jj1 * arrsize + ii1];
    long maxzone = (long) dcimg[piaacmcopticaldesign.zonezID].md[0].size[0];

    if (st->zi - 1 > maxzone - 1)
    {
        printf(
            "ERROR: Zone %d does not exist (image %s has size %ld %ld)   pix %ld %ld   %ld\n",
            (int) st->zi,
            dcimg[piaacmcopticaldesign.zonezID].md[0].name,
            maxzone,
            (long) dcimg[piaacmcopticaldesign.zonezID].md[0].size[1],
            jj1,
            jj1,
            arrsize);
        exit(0);
    }

    if (st->zi > 0)
    {
        long zidx  = piaacmcopticaldesign.focmNBzone * k + st->zi - 1;
        st->t      = dcimg[piaacmcopticaldesign.zonezID].array.D[st->zi - 1];
        st->amp    = dcimg[piaacmcopticaldesign.zoneaID].array.D[st->zi - 1];
        st->cospha = (float) ctx->cosphaarray[zidx];
        st->sinpha = (float) ctx->sinphaarray[zidx];
    }
}

/**
 * mkFocalPlaneMask_eval_subpix - Accumulate one subpixel sample into pixel complex amplitude
 * @ctx: Active FPM synthesis context.
 * @k:   Wavelength slice index.
 * @x:   Subpixel horizontal coordinate in focal plane [m].
 * @y:   Subpixel vertical coordinate in focal plane [m].
 * @st:  Running subpixel accumulator for current pixel.
 */
static void mkFocalPlaneMask_eval_subpix(
    const struct fpm_mask_ctx *ctx,
    int                        k,
    double                     x,
    double                     y,
    struct fpm_subpix_state   *st)
{
    double r = sqrt(x * x + y * y);

    st->zi     = 0;
    st->cospha = 1.0f;
    st->sinpha = 0.0f;
    st->amp    = 1.0;

    if (ctx->OuterCone == 1 && r > 0.9 * piaacmcopticaldesign.fpmRad
        && r < piaacmcopticaldesign.fpmOuterConeRad)
    {
        st->t = piaacmcopticaldesign.fpmOuterConeZ * (piaacmcopticaldesign.fpmOuterConeRad - r)
                / (piaacmcopticaldesign.fpmOuterConeRad - piaacmcopticaldesign.fpmRad);
        long zidx  = piaacmcopticaldesign.focmNBzone * k + st->zi - 1;
        st->cospha = (float) ctx->cosphaarray[zidx];
        st->sinpha = (float) ctx->sinphaarray[zidx];
        st->amp    = 1.0;
    }

    long   arrsize = piaacmcopticaldesign.fpmarraysize;
    double scale   = 0.5 / piaacmcopticaldesign.fpmRad * piaacmcparams.FPMSCALEFACTOR;
    long   ii1     = (long) ((0.5 + x * scale) * arrsize + 0.5);
    long   jj1     = (long) ((0.5 + y * scale) * arrsize + 0.5);

    if (ii1 > -1 && ii1 < arrsize && jj1 > -1 && jj1 < arrsize)
    {
        mkFocalPlaneMask_lookup_zone(ctx, k, r, ii1, jj1, st);
    }

    if (ctx->FPMmode == 1)
    {
        st->retmp += 1.0 - st->amp * st->cospha;
        st->imtmp += -st->amp * st->sinpha;
    }
    else if (ctx->mode == st->zi)
    {
        st->amp = 1.0;
        st->retmp += st->amp;
    }

    st->ttmp += st->t;
    st->zonetmp += 1.0 * st->zi;
}

/**
 * mkFocalPlaneMask_eval_pixel - Compute FPM complex amplitude and sag for pixel (ii, jj)
 * @ctx:     Active FPM synthesis context.
 * @k:       Wavelength slice index.
 * @ii:      Pixel column index.
 * @jj:      Pixel row index.
 * @fpscale: Focal plane spatial scale at wavelength k [m/pix].
 */
static void mkFocalPlaneMask_eval_pixel(
    const struct fpm_mask_ctx *ctx,
    int                        k,
    uint32_t                   ii,
    uint32_t                   jj,
    double                     fpscale)
{
    double   x      = (1.0 * ii - ctx->size / 2) * fpscale;
    double   y      = (1.0 * jj - ctx->size / 2) * fpscale;
    double   r      = sqrt(x * x + y * y);
    float    cospha = 1.0f;
    float    sinpha = 0.0f;
    double   amp    = 1.0;
    uint64_t pixidx = (uint64_t) k * ctx->size2 + (uint64_t) jj * ctx->size + ii;

    if (ctx->OuterCone == 1 && r > 0.9 * piaacmcopticaldesign.fpmRad
        && r < piaacmcopticaldesign.fpmOuterConeRad)
    {
        double t_cone = piaacmcopticaldesign.fpmOuterConeZ
                        * (piaacmcopticaldesign.fpmOuterConeRad - r)
                        / (piaacmcopticaldesign.fpmOuterConeRad - piaacmcopticaldesign.fpmRad);
        double pha = OpticsMaterials_pha_lambda(
            piaacmcopticaldesign.fpmmaterial_code,
            t_cone,
            piaacmcopticalsystem.lambdaarray[k]);
        cospha = cosf((float) pha);
        sinpha = sinf((float) pha);
    }

    dcimg[ctx->IDzone].array.F[pixidx] = 0.0f;

    if (r >= 1.1 * piaacmcopticaldesign.fpmRad)
    {
        if (ctx->FPMmode == 1)
        {
            dcimg[ctx->ID].array.CF[pixidx].re = (float) (1.0 - amp * cospha);
            dcimg[ctx->ID].array.CF[pixidx].im = (float) (-amp * sinpha);
        }
        else
        {
            dcimg[ctx->ID].array.CF[pixidx].re = 0.0f;
            dcimg[ctx->ID].array.CF[pixidx].im = 0.0f;
        }
        dcimg[ctx->IDsag].array.F[pixidx]  = 0.0f;
        dcimg[ctx->IDzone].array.F[pixidx] = 0.0f;
        return;
    }

    struct fpm_subpix_state st;
    memset(&st, 0, sizeof(st));
    st.cospha = 1.0f;
    st.amp    = 1.0;

    for (long iii = 0; iii < ctx->NBsubPix; iii++)
    {
        double sx = (1.0 * ii - ctx->size / 2 + 1.0 * (0.5 + iii) / ctx->NBsubPix - 0.5) * fpscale;
        for (long jjj = 0; jjj < ctx->NBsubPix; jjj++)
        {
            double sy = (1.0 * jj - ctx->size / 2 + 1.0 * (0.5 + jjj) / ctx->NBsubPix - 0.5)
                        * fpscale;
            mkFocalPlaneMask_eval_subpix(ctx, k, sx, sy, &st);
        }
    }

    double norm                        = (double) (ctx->NBsubPix * ctx->NBsubPix);
    dcimg[ctx->ID].array.CF[pixidx].re = (float) (st.retmp / norm);
    dcimg[ctx->ID].array.CF[pixidx].im = (float) (st.imtmp / norm);
    dcimg[ctx->IDsag].array.F[pixidx]  = (float) (st.ttmp / norm);
    dcimg[ctx->IDzone].array.F[pixidx] = (float) (st.zonetmp / norm);
}

/**
 * mkFocalPlaneMask_compute_slice - Synthesize one monochromatic FPM slice at index k
 * @ctx: Active FPM synthesis context.
 * @k:   Wavelength slice index.
 */
static void mkFocalPlaneMask_compute_slice(
    const struct fpm_mask_ctx *ctx,
    int                        k)
{
    double fpscale = (2.0 * piaacmcopticaldesign.beamrad / piaacmcopticaldesign.pixscale)
                     / piaacmcopticaldesign.size / piaacmcopticaldesign.fpzfactor
                     * piaacmcopticalsystem.lambdaarray[k] * piaacmcopticaldesign.Fratio;

    printf(
        "LAMBDA %3d / %3d = %10.5g m    SCALE = %10.5g m/pix   size=%4ul  rad=%g\n",
        k,
        ctx->nblambda,
        piaacmcopticalsystem.lambdaarray[k],
        fpscale,
        ctx->size,
        piaacmcopticaldesign.fpmRad);
    printf(
        "Zone 0 amplitude [%ld]: %lf\n",
        piaacmcopticaldesign.zoneaID,
        dcimg[piaacmcopticaldesign.zoneaID].array.D[0]);
    printf("Zone 0 thickness: %g\n", dcimg[piaacmcopticaldesign.zonezID].array.D[0]);
    printf("Number of zones: %ld\n", piaacmcopticaldesign.focmNBzone);
    printf("piaacmcopticaldesign.fpmRad = %g m\n", piaacmcopticaldesign.fpmRad);
    printf("piaacmcopticaldesign.fpmCentConeRad = %g m\n", piaacmcopticaldesign.fpmCentConeRad);
    printf(
        "piaacmcopticaldesign.fpmOuterConeRad [%d] = %g m\n",
        (int) ctx->OuterCone,
        piaacmcopticaldesign.fpmOuterConeRad);

    for (uint32_t ii = 0; ii < ctx->size; ii++)
    {
        for (uint32_t jj = 0; jj < ctx->size; jj++)
        {
            mkFocalPlaneMask_eval_pixel(ctx, k, ii, jj, fpscale);
        }
    }
}

/**
 * mkFocalPlaneMask_compute_cube - Evaluate all polychromatic FPM slices across wavelengths
 * @ctx: Initialized FPM synthesis context with precomputed zone phasors.
 */
void mkFocalPlaneMask_compute_cube(const struct fpm_mask_ctx *ctx)
{
#ifdef HAVE_LIBGOMP
#pragma omp parallel for default(shared)
#endif
    for (int k = 0; k < ctx->nblambda; k++)
    {
        mkFocalPlaneMask_compute_slice(ctx, k);
    }
}
