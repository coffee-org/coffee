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





extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTPIAACMCDESIGN *piaacmc;



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


    for(k = 0; k < data.image[piaacmc[0].zonezID].md[0].size[0]; k++)
    {
        piaacmcsimul_var.zonez_array[k] = gsl_vector_get(v, k);
    }

    FUNC_CHECK_RETURN(
        PIAACMCsimul_achromFPMsol_eval(
            piaacmcsimul_var.fpmresp_array,
            piaacmcsimul_var.zonez_array,
            piaacmcsimul_var.dphadz_array,
            piaacmcsimul_var.outtmp_array,
            piaacmcsimul_var.vsize,
            data.image[piaacmc[0].zonezID].md[0].size[0],
            piaacmc[0].nblambda,
            &value
        );
    );
    value /= piaacmcsimul_var.CnormFactor * piaacmcsimul_var.SCORINGTOTAL *
             piaacmc[0].nblambda;

    piaacmcsimul_var.LOOPCNT++;

    if(outval != NULL)
    {
        *outval = value;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}



