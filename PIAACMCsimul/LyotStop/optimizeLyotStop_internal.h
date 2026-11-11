/**
 * @file    optimizeLyotStop_internal.h
 * @brief   Internal context and evaluation helpers for multi-plane Lyot stop optimization
 */

#ifndef PIAACMCSIMUL_LYOTSTOP_OPTIMIZELYOTSTOP_INTERNAL_H
#define PIAACMCSIMUL_LYOTSTOP_OPTIMIZELYOTSTOP_INTERNAL_H

#include <stdint.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

struct lyot_opt_ctx
{
    uint32_t xsize;
    uint32_t ysize;
    uint64_t xysize;
    long     nblambda;
    long     NBz;
    long     NBmasks;
    imageID  IDincohc;
    imageID  IDzone;
    float   *zarray;
    float   *rinarray;
    float   *routarray;
    float   *rarray;
    double  *totarray;
    double  *tot2array;
};

/* Propagate and spatially smooth the on-axis complex wavefront across NBz axial planes. */
errno_t lyot_opt_propagate_onaxis_cube(
    const struct lyot_opt_ctx *ctx,
    const char       *restrict IDamp_name,
    const char       *restrict IDpha_name);

/* Score incoherent pupil concentration per zone and select optimal axial conjugations. */
errno_t lyot_opt_select_conjugations(
    const struct lyot_opt_ctx *ctx,
    double                     alpha);

#endif /* PIAACMCSIMUL_LYOTSTOP_OPTIMIZELYOTSTOP_INTERNAL_H */
