#ifndef PIAACMCSIMUL_FOCALPLANEMASK_FPM_PROCESS_H
#define PIAACMCSIMUL_FOCALPLANEMASK_FPM_PROCESS_H

#include "milk_runnable.h"

extern const MILK_RUNNABLE runnable_piaacmc_fpmprocess;

/* Quantize continuous FPM zone sags into discrete lithographic etching levels. */
errno_t PIAACMC_FPM_process(
    const char *restrict FPMsag_name,
    const char *restrict zonescoord_name,
    long                 NBexp,
    const char *restrict outname);

#endif /* PIAACMCSIMUL_FOCALPLANEMASK_FPM_PROCESS_H */
