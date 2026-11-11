/**
 * @file    PIAACMCsimul_run_internal.h
 * @brief   Internal context and parameter helpers for PIAACMCsimul_run and Mode 13 loop
 */

#ifndef PIAACMCSIMUL_RUN_INTERNAL_H
#define PIAACMCSIMUL_RUN_INTERNAL_H

#include <sys/time.h>

#include "milk_runnable.h"

/**
 * @brief Mutable state across outer iterations of the Mode 13 multi-zone FPM search.
 */
struct mode13_loop_ctx
{
    struct timeval start;
    char           stopfile[STRINGMAXLEN_FULLFILENAME];
    double         searchtime;
    double         bestval;
    double         sag0;
    imageID        IDbestsol;
    long           cnt00;
    long           cnt0;
    int            zeroST;
    int            loopiter;
    int            loopin;
};

/* Load CLI variable overrides into piaacmcparams and initialize piaacmcopticaldesign. */
errno_t PIAACMCsimul_run_setparam_variables(
    const char *confindex,
    long        mode);

/* Read Mode 13 search time budget [s] from searchtime.txt (default 36000 s). */
double PIAACMCsimul_run_read_searchtime(void);

/* Select starting sag vector, apply 2*pi wraps, and write sag0.txt / fpsagtest.txt. */
void PIAACMCsimul_run_prep_mode13_zones(struct mode13_loop_ctx *ctx);

#endif /* PIAACMCSIMUL_RUN_INTERNAL_H */
