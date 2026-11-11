/**
 * @file    setup_PIAA_shapes.h
 * @brief   Initial 2D prolate apodization, radial PIAA sag solver, and modal fit pipeline
 */

#ifndef PIAACMCSIMUL_SETUP_PIAA_SHAPES_H
#define PIAACMCSIMUL_SETUP_PIAA_SHAPES_H

#include <stdint.h>

/* Load or construct PIAA 2D apodization, radial sag tables, and Cmodes/Fmodes fits. */
errno_t setupPIAAshapes(
    long     piaacmctype,
    uint32_t size,
    float    beamrad,
    int      loaded,
    int     *saveconf);

#endif /* PIAACMCSIMUL_SETUP_PIAA_SHAPES_H */
