#ifndef PIAACMCSIMUL_ACHROMFPMSOL_EVAL_ZONEZDERIVATIVE_H
#define PIAACMCSIMUL_ACHROMFPMSOL_EVAL_ZONEZDERIVATIVE_H

errno_t PIAACMCsimul_achromFPMsol_eval_zonezderivative(long    zone,
                                                       double *fpmresp_array,
                                                       double *zonez_array,
                                                       double *dphadz_array,
                                                       double *outtmp_array,
                                                       long    vsize,
                                                       long    nbz,
                                                       long    nbl);

#endif
