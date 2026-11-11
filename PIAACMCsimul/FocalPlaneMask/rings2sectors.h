#ifndef PIAACMCSIMUL_FOCALPLANEMASK_RING2SECTORS_H
#define PIAACMCSIMUL_FOCALPLANEMASK_RING2SECTORS_H

#include "milk_runnable.h"

extern const MILK_RUNNABLE runnable_piaacmc_ring2sect;

/* Expand radial ring FPM sag values onto azimuthal sector zones. */
errno_t rings2sectors(
    const char *IDin_name,
    const char *sectfname,
    const char *IDout_name,
    imageID    *outID);

#endif /* PIAACMCSIMUL_FOCALPLANEMASK_RING2SECTORS_H */
