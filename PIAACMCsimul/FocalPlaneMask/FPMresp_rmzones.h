/**
 * @file    FPMresp_rmzones.h
 * @brief   Outer-zone truncation utility for FPM complex response matrices
 */

#ifndef PIAACMCSIMUL_FOCALPLANEMASK_FPMRESP_RMZONES_H
#define PIAACMCSIMUL_FOCALPLANEMASK_FPMRESP_RMZONES_H

#include "milk_runnable.h"

extern const MILK_RUNNABLE runnable_piaacmc_fpmresprm;

/* Collapse zones beyond NBzones into the unmodulated reference row (y = 0). */
errno_t FPMresp_rmzones(
    const char *restrict FPMresp_in_name,
    const char *restrict FPMresp_out_name,
    long                 NBzones,
    imageID             *outID);

#endif /* PIAACMCSIMUL_FOCALPLANEMASK_FPMRESP_RMZONES_H */
