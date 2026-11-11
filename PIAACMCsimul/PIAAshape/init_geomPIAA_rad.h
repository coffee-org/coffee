/**
 * @file    init_geomPIAA_rad.h
 * @brief   1D geometric ray-trace solver for lossless radial PIAA mirror profiles
 */

#ifndef INIT_GEOMPIAA_RAD_H
#define INIT_GEOMPIAA_RAD_H

/* Integrate radial PIAA ray-mapping equations from cosine fit IDapofit_name -> PIAA_Mshapes.txt. */
errno_t init_geomPIAA_rad(const char *restrict IDapofit_name);

#endif /* INIT_GEOMPIAA_RAD_H */
