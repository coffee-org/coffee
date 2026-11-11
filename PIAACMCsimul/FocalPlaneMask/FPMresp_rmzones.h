#ifndef PIAACMCSIMUL_FOCALPLANEMASK_FPMRESP_RMZONES_H
#define PIAACMCSIMUL_FOCALPLANEMASK_FPMRESP_RMZONES_H

#include "milk_runnable.h"

extern const MILK_RUNNABLE runnable_piaacmc_fpmresprm;

/* Collapse outer zones of an FPM complex response matrix into the background slice. */
errno_t FPMresp_rmzones(
    const char *restrict FPMresp_in_name,
    const char *restrict FPMresp_out_name,
    long                 NBzones,
    imageID             *outID);

#endif /* PIAACMCSIMUL_FOCALPLANEMASK_FPMRESP_RMZONES_H */
