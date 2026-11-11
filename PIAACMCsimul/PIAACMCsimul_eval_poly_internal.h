/**
 * @file    PIAACMCsimul_eval_poly_internal.h
 * @brief   Shared evaluation state and sensitivity driver for Mode 100 polychromatic evaluation
 */

#ifndef PIAACMCSIMUL_EVAL_POLY_INTERNAL_H
#define PIAACMCSIMUL_EVAL_POLY_INTERNAL_H

#include <stdint.h>

#include "milk_runnable.h"

/**
 * @brief Shared dimensions, calibration constants, and cube IDs across Mode 100 evaluation.
 */
struct eval_poly_ctx
{
    double   fpmradld;
    double   centobs0;
    double   centobs1;
    double   ldoffset;
    double   eval_sepstepld;
    double   eval_sepmaxld;
    double   valref;
    float    avpeak;
    uint32_t xsize;
    uint32_t ysize;
    uint32_t zsize;
    uint64_t xysize;
    imageID  IDopderrC;
    long     nbOPDerr;
    imageID  IDps;
    imageID  IDps_re;
    imageID  IDps_im;
    imageID  IDps_COH;
    imageID  IDps_INC;
    long     NBpt;
};

/* Evaluate pointing offsets, OPD error modes, and coherent/incoherent contrast curves. */
errno_t PIAACMCsimul_eval_poly_sens(struct eval_poly_ctx *ctx);

#endif /* PIAACMCSIMUL_EVAL_POLY_INTERNAL_H */
