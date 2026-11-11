/**
 * @file    PIAACMCsimul_achromFPMsol_eval.h
 * @brief   Fast linear-superposition evaluation of polychromatic FPM dark-hole contrast
 */

#ifndef PIAACMCSIMUL_ACHROMFPMSOL_EVAL_H
#define PIAACMCSIMUL_ACHROMFPMSOL_EVAL_H

/* Evaluate mean dark-hole contrast outval from precomputed FPM zone response cube. */
errno_t PIAACMCsimul_achromFPMsol_eval(
    double *restrict fpmresp_array,
    double *restrict zonez_array,
    double *restrict dphadz_array,
    double *restrict outtmp_array,
    long             vsize,
    long             nbz,
    long             nbl,
    double          *outval);

#endif /* PIAACMCSIMUL_ACHROMFPMSOL_EVAL_H */
