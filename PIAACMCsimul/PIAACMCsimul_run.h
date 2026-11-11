/**
 * @file    PIAACMCsimul_run.h
 * @brief   CLI entry point and multi-pass Mode 13 optimization driver
 */

#ifndef PIAACMCSIMUL_RUN_H
#define PIAACMCSIMUL_RUN_H

#include "milk_runnable.h"

extern const MILK_RUNNABLE runnable_piaacmc_run;

/* Read mc_var overrides and run PIAACMC design/simulation mode on confindex. */
errno_t PIAACMCsimul_run(
    const char *confindex,
    long        mode);

#endif /* PIAACMCSIMUL_RUN_H */
