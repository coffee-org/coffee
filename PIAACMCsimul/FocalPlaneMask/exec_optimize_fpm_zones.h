/**
 * @file    exec_optimize_fpm_zones.h
 * @brief   Iterative regularized SVD / gradient optimization of FPM zone sags (Mode 13)
 */

#ifndef PIAACMCSIMUL_FOCALPLANEMASK_EXEC_OPTIMIZE_FPM_ZONES_H
#define PIAACMCSIMUL_FOCALPLANEMASK_EXEC_OPTIMIZE_FPM_ZONES_H

/* Optimize multi-zone FPM sags [m] using analytical Jacobian SVD line search (Mode 13). */
errno_t exec_optimize_fpm_zones(void);

#endif /* PIAACMCSIMUL_FOCALPLANEMASK_EXEC_OPTIMIZE_FPM_ZONES_H */
