/**
 * @file    PIAACMCsimul_f_evalmask.c
 * @brief   PIAA-type coronagraph design, run
 *
 * Can design both APLCMC and PIAACMC coronagraphs
 *
 *
 */



#include <stdlib.h>
#include <stdio.h>


#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>



// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "PIAACMCsimul.h"
#include "PIAACMCsimul_achromFPMsol_eval.h"




errno_t f_evalmask(
    const gsl_vector *v,
    void *params,
    double *outval
)
{
    DEBUG_TRACE_FSTART();

    double *p = (double *)params;
    double value;
    long k;

    (void) p;


    for(k = 0; k < data.image[piaacmcopticaldesign.zonezID].md[0].size[0]; k++)
    {
        piaacmcparams.zonez_array[k] = gsl_vector_get(v, k);
    }

    FUNC_CHECK_RETURN(
        PIAACMCsimul_achromFPMsol_eval(
            piaacmcparams.fpmresp_array,
            piaacmcparams.zonez_array,
            piaacmcparams.dphadz_array,
            piaacmcparams.outtmp_array,
            piaacmcparams.vsize,
            data.image[piaacmcopticaldesign.zonezID].md[0].size[0],
            piaacmcopticaldesign.nblambda,
            &value
        );
    );
    value /= piaacmcparams.CnormFactor * piaacmcparams.SCORINGTOTAL *
             piaacmcopticaldesign.nblambda;

    piaacmcparams.LOOPCNT++;

    if(outval != NULL)
    {
        *outval = value;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}



