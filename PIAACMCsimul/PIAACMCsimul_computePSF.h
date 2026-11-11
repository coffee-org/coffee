/**
 * @file    PIAACMCsimul_computePSF.h
 * @brief   Polychromatic coronagraphic PSF and dark-hole contrast evaluation
 */

#ifndef PIAACMCSIMUL_COMPUTEPSF_H
#define PIAACMCSIMUL_COMPUTEPSF_H

/* Propagate point or extended source at (xld, yld) [l/D] and return mean scoring contrast. */
errno_t PIAACMCsimul_computePSF(
    float   xld,
    float   yld,
    long    startelem,
    long    endelem,
    int     savepsf,
    int     sourcesize,
    int     extmode,
    int     outsave,
    double *contrastval);

#endif /* PIAACMCSIMUL_COMPUTEPSF_H */
