/**
 * @file    init_piaacmcopticalsystem_internal.h
 * @brief   Internal element-assembly helpers for init_piaacmcopticalsystem
 */

#ifndef PIAACMCSIMUL_INIT_PIAACMCOPTICALSYSTEM_INTERNAL_H
#define PIAACMCSIMUL_INIT_PIAACMCOPTICALSYSTEM_INTERNAL_H

#include "milk_runnable.h"

/* Append pre/post PIAA masks and PIAA aspheric elements 0 and 1 to piaacmcopticalsystem. */
errno_t piaacmc_optsyst_add_piaa_stage(
    long  *elem,
    long   size,
    double beamradpix);

/* Build and append the focal plane mask element unless PIAACMC_NOFPM is defined. */
errno_t piaacmc_optsyst_add_fpm(long *elem);

/* Append inverse PIAA optics, Lyot stop planes, and output pupil stop to piaacmcopticalsystem. */
void piaacmc_optsyst_add_lyot_and_backend(
    long *elem,
    long  size);

#endif /* PIAACMCSIMUL_INIT_PIAACMCOPTICALSYSTEM_INTERNAL_H */
