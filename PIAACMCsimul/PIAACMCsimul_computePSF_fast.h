/**
 * @file    PIAACMCsimul_computePSF_fast.h
 * @brief   Focal-plane scoring mask, fast FPMresp evaluator, and PSF flux/Cnorm telemetry
 */

#ifndef PIAACMCSIMUL_COMPUTEPSF_FAST_H
#define PIAACMCSIMUL_COMPUTEPSF_FAST_H

#include <stdint.h>

/* Ensure focal-plane dark-hole `scoringmask` and pixel lookup tables (`pixindex`) exist. */
errno_t PIAACMCsimul_init_scoringmask(
    uint32_t size,
    float    focscale);

/* Evaluate dark-hole electric field `imvect` and mean contrast from pre-computed `FPMresp`. */
errno_t PIAACMCsimul_eval_fast_fpmresp(
    uint32_t size,
    float    focscale,
    double  *avContrast_out);

/* Print per-element fluxes and write optional flux_*.txt and CnormFactor.txt files. */
void PIAACMCsimul_record_psf_telemetry(
    uint32_t size,
    float    focscale,
    int      sourcesize,
    int      outsave,
    int      is_ext);

#endif /* PIAACMCSIMUL_COMPUTEPSF_FAST_H */
