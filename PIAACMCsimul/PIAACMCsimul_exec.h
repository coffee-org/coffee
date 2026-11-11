/**
 * @file    PIAACMCsimul_exec.h
 * @brief   Primary mode dispatcher for PIAACMC / APLCMC design and evaluation steps
 */

#ifndef PIAACMCSIMUL_EXEC_H
#define PIAACMCSIMUL_EXEC_H

/* Execute PIAACMC design/evaluation operation 'mode' in piaacmcconf_i<confindex>. */
errno_t PIAACMCsimul_exec(
    const char *confindex,
    long        mode);

#endif /* PIAACMCSIMUL_EXEC_H */
