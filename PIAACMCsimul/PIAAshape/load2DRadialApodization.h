/**
 * @file    load2DRadialApodization.h
 * @brief   Modal projection of 2D radial pupil apodizer onto radial cosine basis
 */

#ifndef LOAD2DRADIALAPODIZATION_H
#define LOAD2DRADIALAPODIZATION_H

/* Fit 2D apodizer IDapo_name (beam radius beamradpix [pix]) onto cosine vector IDapofit_name. */
errno_t load2DRadialApodization(
    const char *restrict IDapo_name,
    float                beamradpix,
    const char *restrict IDapofit_name);

#endif /* LOAD2DRADIALAPODIZATION_H */
