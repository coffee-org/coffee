/**
 * @file    init_piaacmcopticalsystem.h
 * @brief   Population of OptSystProp OPTSYST element sequence from OPTPIAACMCDESIGN
 */

#ifndef PIAACMCSIMUL_INIT_PIAACMCOPTICALSYSTEM_H
#define PIAACMCSIMUL_INIT_PIAACMCOPTICALSYSTEM_H

/* Assemble optical element train in piaacmcopticalsystem for source offset (xld, yld) [l/D]. */
errno_t init_piaacmcopticalsystem(
    double xld,
    double yld);

#endif /* PIAACMCSIMUL_INIT_PIAACMCOPTICALSYSTEM_H */
