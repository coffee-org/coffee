/**
 * @file    FPM_process.h
 * @brief   Discrete lithographic sag quantization of multi-zone focal plane masks
 */

#ifndef PIAACMCSIMUL_FOCALPLANEMASK_FPM_PROCESS_H
#define PIAACMCSIMUL_FOCALPLANEMASK_FPM_PROCESS_H

#include "milk_runnable.h"

extern const MILK_RUNNABLE runnable_piaacmc_fpmprocess;

/* Quantize continuous FPM zone sags [m] into 2^NBexp discrete etching levels. */
errno_t PIAACMC_FPM_process(
    const char *restrict FPMsag_name,
    const char *restrict zonescoord_name,
    long                 NBexp,
    const char *restrict outname);

#endif /* PIAACMCSIMUL_FOCALPLANEMASK_FPM_PROCESS_H */
