/**
 * @file    PIAACMC_f_evalmask.h
 * @brief   GSL Nelder-Mead / annealing merit callback for FPM zone sag search
 */

#ifndef PIAACMCSIMUL_F_EVALMASK_H
#define PIAACMCSIMUL_F_EVALMASK_H

#include <gsl/gsl_vector.h>

/* Evaluate regularized polychromatic dark-hole contrast for candidate FPM vector v. */
errno_t f_evalmask(
    const gsl_vector *v,
    void             *params,
    double           *outval);

#endif /* PIAACMCSIMUL_F_EVALMASK_H */
