/**
 * @file    PIAACMCsimul_achromFPMsol_eval_zonezderivative.h
 * @brief   Analytical complex electric-field derivative with respect to FPM zone sag
 */

#ifndef PIAACMCSIMUL_ACHROMFPMSOL_EVAL_ZONEZDERIVATIVE_H
#define PIAACMCSIMUL_ACHROMFPMSOL_EVAL_ZONEZDERIVATIVE_H

/* Compute complex focal-plane field derivative wrt sag [m] of FPM index 'zone'. */
errno_t PIAACMCsimul_achromFPMsol_eval_zonezderivative(
    long    zone,
    double *fpmresp_array,
    double *zonez_array,
    double *dphadz_array,
    double *outtmp_array,
    long    vsize,
    long    nbz,
    long    nbl);

#endif /* PIAACMCSIMUL_ACHROMFPMSOL_EVAL_ZONEZDERIVATIVE_H */
