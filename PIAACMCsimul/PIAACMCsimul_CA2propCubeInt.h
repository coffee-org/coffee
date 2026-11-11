/**
 * @file    PIAACMCsimul_CA2propCubeInt.h
 * @brief   Axial Fresnel intensity cube generation from complex pupil field
 */

#ifndef PIAACMCSIMUL_CA2PROPCUBEINT_H
#define PIAACMCSIMUL_CA2PROPCUBEINT_H

/* Propagate (IDamp_name, IDpha_name) over NBz planes in [zmin, zmax] [m] into 3D intensity cube. */
errno_t PIAACMCsimul_CA2propCubeInt(
    const char *restrict IDamp_name,
    const char *restrict IDpha_name,
    float                zmin,
    float                zmax,
    long                 NBz,
    const char *restrict IDout_name,
    imageID             *outID);

#endif /* PIAACMCSIMUL_CA2PROPCUBEINT_H */
