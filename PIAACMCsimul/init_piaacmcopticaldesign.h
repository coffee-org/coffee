/**
 * @file    init_piaacmcopticaldesign.h
 * @brief   Initialization of PIAACMC optical prescription, modal bases, PIAA sags, and masks
 */

#ifndef PIAACMCSIMUL_INIT_PIAACMCOPTICALDESIGN_H
#define PIAACMCSIMUL_INIT_PIAACMCOPTICALDESIGN_H

#include <stdint.h>

#define INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT         0x0000000000000000ULL
#define INIT_PIAACMCOPTICALDESIGN_MODE__READCONF        0x0000000000000001ULL
#define INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF 0x0000000000000002ULL
#define INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL     0x0000000000000004ULL
#define INIT_PIAACMCOPTICALDESIGN_MODE__WSCMODE         0x0000000000000008ULL

/* Initialize or load optical design parameters, PIAA modes/shapes, FPM, and Lyot stops. */
errno_t init_piaacmcopticaldesign(
    double    fpmradld,
    double    centobs0,
    double    centobs1,
    uint64_t  flags,
    uint64_t *status);

#endif /* PIAACMCSIMUL_INIT_PIAACMCOPTICALDESIGN_H */
