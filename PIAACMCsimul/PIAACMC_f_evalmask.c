/**
 * @file    PIAACMCsimul_f_evalmask.c
 * @brief   PIAA-type coronagraph design, run
 * 
 * Can design both APLCMC and PIAACMC coronagraphs
 *  
 * @author  O. Guyon
 * @date    21 nov 2017
 *
 * 
 * @bug No known bugs.
 * 
 */



#include <stdlib.h>
#include <stdio.h>


#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>



// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "PIAACMCsimul/PIAACMCsimul.h"






extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTPIAACMCDESIGN *piaacmc;




static double f_evalmask (const gsl_vector *v, void *params)
{
    double *p = (double *)params;
    double value;
    long k;


	#ifdef PIAASIMUL_LOGFUNC1
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif


    for(k=0; k<data.image[piaacmc[0].zonezID].md[0].size[0]; k++)
        piaacmcsimul_var.zonez_array[k] = gsl_vector_get(v, k);


    value = PIAACMCsimul_achromFPMsol_eval(piaacmcsimul_var.fpmresp_array, piaacmcsimul_var.zonez_array, piaacmcsimul_var.dphadz_array, piaacmcsimul_var.outtmp_array, piaacmcsimul_var.vsize, data.image[piaacmc[0].zonezID].md[0].size[0], piaacmc[0].nblambda);
    value /= piaacmcsimul_var.CnormFactor*piaacmcsimul_var.SCORINGTOTAL*piaacmc[0].nblambda;

    piaacmcsimul_var.LOOPCNT++;

    return (value);

}



