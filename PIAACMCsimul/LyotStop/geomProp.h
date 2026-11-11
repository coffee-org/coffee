#ifndef PIAACMCSIMUL_LYOTSTOP_GEOMPROP_H
#define PIAACMCSIMUL_LYOTSTOP_GEOMPROP_H

#include "milk_runnable.h"

extern const MILK_RUNNABLE runnable_piaacmc_geomprop;

/* Propagate 2D pupil intensity geometrically using local surface sag slopes. */
errno_t PIAACMCsimul_geomProp(
    const char *restrict IDin_name,
    const char *restrict IDsag_name,
    const char *restrict IDout_name,
    const char *restrict IDoutcnt_name,
    double               drindex,
    double               pscale,
    double               zprop,
    double               krad,
    double               kstep,
    double               rlim,
    imageID             *outID);

#endif /* PIAACMCSIMUL_LYOTSTOP_GEOMPROP_H */
