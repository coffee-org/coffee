/**
 * @file    mkFocalPlaneMask_internal.h
 * @brief   Internal shared context and slice kernel declarations for mkFocalPlaneMask
 */

#ifndef PIAACMCSIMUL_FOCALPLANEMASK_MKFOCALPLANEMASK_INTERNAL_H
#define PIAACMCSIMUL_FOCALPLANEMASK_MKFOCALPLANEMASK_INTERNAL_H

#include <stdint.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

struct fpm_mask_ctx
{
    imageID      ID;
    imageID      IDz;
    imageID      IDsag;
    imageID      IDzone;
    uint32_t     size;
    uint64_t     size2;
    int          nblambda;
    int          mode;
    int          FPMmode;
    long         NBsubPix;
    uint_fast8_t CentCone;
    uint_fast8_t OuterCone;
    double      *tarray;
    double      *aarray;
    double      *phaarray;
    double      *cosphaarray;
    double      *sinphaarray;
};

/* Evaluate polychromatic complex-amplitude FPM cube slices across wavelengths. */
void mkFocalPlaneMask_compute_cube(const struct fpm_mask_ctx *ctx);

#endif /* PIAACMCSIMUL_FOCALPLANEMASK_MKFOCALPLANEMASK_INTERNAL_H */
