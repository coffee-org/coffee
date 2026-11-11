/**
 * @file    exec_computePSF_no_fpm.h
 * @brief   Unocculted reference PSF calibration without focal plane mask (Mode 3)
 */

#ifndef EXEC_COMPUTEPSF_NO_FPM_H
#define EXEC_COMPUTEPSF_NO_FPM_H

/* Propagate on-axis star with FPM removed to calibrate CnormFactor.txt. */
errno_t exec_computePSF_no_fpm(double *outval);

#endif /* EXEC_COMPUTEPSF_NO_FPM_H */
