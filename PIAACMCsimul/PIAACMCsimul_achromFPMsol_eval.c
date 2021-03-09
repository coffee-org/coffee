/**
 * @file    PIAACMCsimul_achromFPMsol_eval.c
 * @brief   PIAA-type coronagraph design, run
 * 
 * Can design both APLCMC and PIAACMC coronagraphs
 *  
 * @author  O. Guyon
 * @date    2017-12-23
 *
 * 
 * @bug No known bugs.
 * 
 */


#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"




extern OPTSYST *optsyst;

extern PIAACMCsimul_varType piaacmcsimul_var;





///
/// solves for focal plane mask solution using pre-computed zone responses
///
///
/// written to be fast, no checking of array sizes
/// all arrays pre-allocated outside this function
///
double PIAACMCsimul_achromFPMsol_eval(
		double *restrict fpmresp_array, 	/// @param[in] fpmresp_array   Mask zones responses, double array
		double *restrict zonez_array, 		/// @param[in] zonez_array     Zone thicknesses, double array
		double *restrict dphadz_array,		/// @param[in] dphadz_array    For each lambda, pha = thickness x dphadt_array[lambdaindex]
		double *restrict outtmp_array, 		/// @param[out] outtmp_array   Output temp array
		long vsize, 
		long nbz, 
		long nbl)
{
	
	double evalval;
	long evalk;

//	long evali;
//	long evalk, evalki, evalki1, evalmz, evalii, evalii1, evalii2, evalkv;
//	double evalcosp, evalsinp, evalre, evalim, evalre1, evalim1, evalpha;
//	double evalv1;

    // axis 0: eval pts (ii)   size = data.image[IDfpmresp].md[0].size[0] -> vsize
    // axis 1: zones (mz)      size = data.image[piaacmc[0].zonezID].md[0].size[0]+1 = nbz+1
    // axis 2: lambda (k)      size = piaacmc[0].nblambda -> nbl
    //
    // indexing :  k*(data.image[piaacmc[0].zonezID].md[0].size[0]+1)*vsize + mz*vsize + ii


	#ifdef PIAASIMUL_LOGFUNC1
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif


    for(evalk=0; evalk<nbl; evalk++) // wavelength index
    {
		long evalki;
		
        evalki = evalk*(nbz+1)*vsize;

        if(optsyst[0].FOCMASKarray[0].mode == 1) // include outer zone
        {
            // outer zone
            long evalii;
            for(evalii=0; evalii<vsize/2; evalii++)
            {
                outtmp_array[evalk*vsize+2*evalii] = fpmresp_array[evalk*(nbz+1)*vsize + 2*evalii]; // mz=0 -> mz*vsize not included in index
                outtmp_array[evalk*vsize+2*evalii+1] = fpmresp_array[evalk*(nbz+1)*vsize + 2*evalii + 1];
            }

			// mask zones
			long evalmz, evalki1, evalkv;
			double evalpha, evalcosp, evalsinp;
            for(evalmz=0; evalmz<nbz; evalmz++)
            {
                evalpha = zonez_array[evalmz]*dphadz_array[evalk];   // CHANGED sign to + on 2017-12-23 to adopt new sign convention
                evalcosp = cos(evalpha);
                evalsinp = sin(evalpha);
                evalki1 = evalki + (evalmz+1)*vsize;
                evalkv = evalk*vsize;

				long evalii1, evalii2;
				double evalre, evalim, evalre1, evalim1;
                for(evalii=0; evalii<vsize/2; evalii++)
                {
                    evalii1 = 2*evalii;
                    evalii2 = 2*evalii+1;
                    evalre = fpmresp_array[evalki1 + evalii1];
                    evalim = fpmresp_array[evalki1 + evalii2];
                    evalre1 = evalre*evalcosp - evalim*evalsinp;
                    evalim1 = evalre*evalsinp + evalim*evalcosp;
                    outtmp_array[evalkv + evalii1] += evalre1;
                    outtmp_array[evalkv + evalii2] += evalim1;
                }
            }

        }
        else  // single zone impulse
        {
			long evalmz, evalki1, evalkv;
			double evalpha, evalcosp, evalsinp;
			
            evalmz = piaacmcsimul_var.focmMode - 1;
            evalpha = zonez_array[evalmz]*dphadz_array[evalk];
            evalcosp = 1.0; //cos(evalpha);
            evalsinp = 0.0; //sin(evalpha);
            evalki1 = evalki + (evalmz+1)*vsize;
            evalkv = evalk*vsize;
            
            long evalii, evalii1, evalii2;
            double evalre, evalim, evalre1, evalim1;
            for(evalii=0; evalii<vsize/2; evalii++)
            {
                evalii1 = 2*evalii;
                evalii2 = 2*evalii+1;
                evalre = fpmresp_array[evalki1 + evalii1];
                evalim = fpmresp_array[evalki1 + evalii2];
                evalre1 = evalre*evalcosp - evalim*evalsinp;
                evalim1 = evalre*evalsinp + evalim*evalcosp;
                outtmp_array[evalkv + evalii1] = evalre1;
                outtmp_array[evalkv + evalii2] = evalim1;
            }

        }
    }


    //	for(evalmz=0; evalmz<nbz; evalmz++)
    //	outtmp_array[nbl*vsize + evalmz] = piaacmcsimul_var.PIAACMC_MASKregcoeff*zonez_array[evalmz]*sqrt(vsize*nbl/nbz);


    evalval = 0.0;
    long evalii;
    double evalv1;
    for(evalii=0; evalii<vsize*nbl; evalii++)
    {
        evalv1 = outtmp_array[evalii];
        evalval += evalv1*evalv1;
    }
    //  evalval /= vsize*nbl;

    // note that evalval is prop to bumber of spectral channels x number of evaluation pixels 
    return evalval;
}

