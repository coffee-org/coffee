/**
 * @file    PIAACMCsimul_regularization.h
 * @brief   PIAA modal surface and FPM zone sag regularization penalties and 1D vectors
 */

#ifndef PIAACMCSIMUL_REGULARIZATION_H
#define PIAACMCSIMUL_REGULARIZATION_H

#include "COREMOD_memory/COREMOD_memory.h"

/* Compute scalar L2 regularization penalty for PIAA element 0/1 Cosine and Fourier modes. */
double PIAACMCsimul_regularization_PIAAshapes_value(void);

/* Compute scalar L2 regularization penalty for physical focal plane mask zone sags. */
double PIAACMCsimul_regularization_fpmsag_value(void);

/* Append weighted PIAA modal deviations to 1D optimization vector starting at index0. */
long PIAACMCsimul_regularization_PIAAshapes_add1Dvector(
    imageID ID1D,
    long    index0);

/* Append weighted FPM zone sag terms to 1D optimization vector (and optional mask). */
long PIAACMCsimul_regularization_fpmsag_add1Dvector(
    imageID ID1D,
    imageID IDmask,
    long    index0);

#endif /* PIAACMCSIMUL_REGULARIZATION_H */
