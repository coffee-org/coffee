/**
 * @file    mkPIAAMshapes_from_RadSag.h
 * @brief   2D radial interpolation of 1D geometric PIAA mirror sag profiles
 */

#ifndef MKPIAASHAPES_FROM_RADSAG_H
#define MKPIAASHAPES_FROM_RADSAG_H

#include <stdint.h>

/* Interpolate 1D table piaa1Dsagfname onto 2D sag maps ID_PIAAM0_name and ID_PIAAM1_name [m]. */
errno_t mkPIAAMshapes_from_RadSag(
    const char *restrict piaa1Dsagfname,
    double               PIAAsep,
    double               beamrad,
    double               r0limfact,
    double               r1limfact,
    uint32_t             size,
    double               beamradpix,
    long                 NBradpts,
    const char *restrict ID_PIAAM0_name,
    const char *restrict ID_PIAAM1_name);

#endif /* MKPIAASHAPES_FROM_RADSAG_H */
