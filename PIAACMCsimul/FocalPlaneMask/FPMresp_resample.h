/**
 * @file    FPMresp_resample.h
 * @brief   Spectral interpolation and dark-hole decimation of FPM response matrices
 */

#ifndef PIAACMCSIMUL_FOCALPLANEMASK_FPMRESP_RESAMPLE_H
#define PIAACMCSIMUL_FOCALPLANEMASK_FPMRESP_RESAMPLE_H

#include "milk_runnable.h"

extern const MILK_RUNNABLE runnable_piaacmc_fpmresprs;

/* Resample 3D FPM complex response cube to NBlambda channels and PTstep stride. */
errno_t PIAACMC_FPMresp_resample(
    const char *restrict FPMresp_in_name,
    const char *restrict FPMresp_out_name,
    long                 NBlambda,
    long                 PTstep,
    imageID             *outID);

#endif /* PIAACMCSIMUL_FOCALPLANEMASK_FPMRESP_RESAMPLE_H */
