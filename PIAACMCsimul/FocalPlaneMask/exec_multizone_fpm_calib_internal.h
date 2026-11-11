/**
 * @file    exec_multizone_fpm_calib_internal.h
 * @brief   Internal shared context and worker/combiner prototypes for FPM calibration
 */

#ifndef PIAACMCSIMUL_FOCALPLANEMASK_EXEC_MULTIZONE_FPM_CALIB_INTERNAL_H
#define PIAACMCSIMUL_FOCALPLANEMASK_EXEC_MULTIZONE_FPM_CALIB_INTERNAL_H

#include "milk_runnable.h"
#include "coffee_compat.h"

struct fpm_calib_ctx
{
    double  fpmradld;
    double  centobs0;
    double  centobs1;
    long    mzoffset;
    long    mzstep;
    imageID ID_FPMresp;
    imageID IDcomb;
    char    fname1[STRINGMAXLEN_FULLFILENAME];
    char    fname2[STRINGMAXLEN_FULLFILENAME];
};

/* Merge per-thread FPMresp FITS files into combined FPMresp cube (combiner thread). */
errno_t fpm_calib_combine_threads(imageID IDfpmresp);

/* Evaluate FPM zone impulse responses assigned to current worker or single process. */
errno_t fpm_calib_compute_worker_zones(
    const struct fpm_calib_ctx *ctx,
    imageID                     IDfpmresp);

#endif /* PIAACMCSIMUL_FOCALPLANEMASK_EXEC_MULTIZONE_FPM_CALIB_INTERNAL_H */
