/**
 * @file    PIAAshape_modes.h
 * @brief   Pupil coordinate grids and radial Cosine / 2D Fourier modal bases for PIAA optics
 */

#ifndef PIAACMCSIMUL_PIAASHAPE_MODES_H
#define PIAACMCSIMUL_PIAASHAPE_MODES_H

/* Pre-compute normalized pupil coordinate maps (xcoord, ycoord, rcoord, PAcoord). */
errno_t make_x_y_r_PA_images(
    long  size,
    float beamrad);

/* Load or generate 2D radial cosine (Cmodes) and Fourier (Fmodes) basis cubes. */
errno_t make_C_F_modes(
    long  size,
    float beamrad);

#endif /* PIAACMCSIMUL_PIAASHAPE_MODES_H */
