/**
 * @file    rings2sectors.h
 * @brief   Mapping of radial ring FPM sag solutions onto 2D hexagonal sector zones
 */

#ifndef PIAACMCSIMUL_FOCALPLANEMASK_RING2SECTORS_H
#define PIAACMCSIMUL_FOCALPLANEMASK_RING2SECTORS_H

#include "milk_runnable.h"

extern const MILK_RUNNABLE runnable_piaacmc_ring2sect;

/* Expand radial ring FPM sag image IDin_name onto sector table sectfname -> IDout_name. */
errno_t rings2sectors(
    const char *IDin_name,
    const char *sectfname,
    const char *IDout_name,
    imageID    *outID);

#endif /* PIAACMCSIMUL_FOCALPLANEMASK_RING2SECTORS_H */
