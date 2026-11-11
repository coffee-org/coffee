/**
 * @file    init_geomPIAA_rad_internal.h
 * @brief   Internal context and radial ODE marcher declaration for 1D PIAA sag synthesis
 */

#ifndef PIAACMCSIMUL_PIAASHAPE_INIT_GEOMPIAA_RAD_INTERNAL_H
#define PIAACMCSIMUL_PIAASHAPE_INIT_GEOMPIAA_RAD_INTERNAL_H

#include "milk_runnable.h"
#include "coffee_compat.h"

struct geom_piaa_rad_ctx
{
    long    NBradpts;
    double *pup0;
    double *pup1;
    double *flux0cumul;
    double *flux1cumul;
    double  flux1_in;
    double  flux1_out;
    double  normcoeff0;
};

/* Invert cumulative flux maps, integrate PIAA mirror slopes, and write PIAA_Mshapes.txt. */
errno_t geom_piaa_rad_march_and_save(const struct geom_piaa_rad_ctx *ctx);

#endif /* PIAACMCSIMUL_PIAASHAPE_INIT_GEOMPIAA_RAD_INTERNAL_H */
