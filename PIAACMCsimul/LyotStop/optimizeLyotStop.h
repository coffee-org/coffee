/**
 * @file    optimizeLyotStop.h
 * @brief   Multi-plane Lyot stop conjugation search across 3D Fresnel intensity cubes
 */

#ifndef OPTIMIZELYOT_H
#define OPTIMIZELYOT_H

/* Search NBz planes in [zmin, zmax] [m] for optimal NBmasks Lyot stop conjugations. */
errno_t optimizeLyotStop(
    const char *restrict IDamp_name,
    const char *restrict IDpha_name,
    const char *restrict Dincohc_name,
    float                zmin,
    float                zmax,
    double               throughput,
    long                 NBz,
    long                 NBmasks,
    double              *outratioval);

#endif /* OPTIMIZELYOT_H */
