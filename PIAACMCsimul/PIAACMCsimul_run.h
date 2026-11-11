#ifndef PIAACMCSIMUL_RUN_H
#define PIAACMCSIMUL_RUN_H

#include "milk_runnable.h"

extern const MILK_RUNNABLE runnable_piaacmc_run;

/* Execute PIAACMC design/simulation step for the given configuration and mode. */
errno_t PIAACMCsimul_run(
    const char *confindex,
    long        mode);

#endif /* PIAACMCSIMUL_RUN_H */
