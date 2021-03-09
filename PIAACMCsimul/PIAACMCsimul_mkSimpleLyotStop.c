/**
 * @file    PIAACMCsimul_mkSimpleLyotStop.c
 * @brief   PIAA-type coronagraph design, run
 * 
 * Can design both APLCMC and PIAACMC coronagraphs
 *  
 * 
 * ## Change log
 * - 20180323  Guyon   minor syntax cleanup   
 * 
 * 
 * @author  O. Guyon
 *
 * 
 * @bug No known bugs.
 * 
 */




/* =============================================================================================== */
/* =============================================================================================== */
/*                                        HEADER FILES                                             */
/* =============================================================================================== */
/* =============================================================================================== */

#include <stdlib.h>
#include <stdio.h>



// milk includes
#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "PIAACMCsimul/PIAACMCsimul.h"



/* =============================================================================================== */
/* =============================================================================================== */
/*                                  GLOBAL DATA DECLARATION                                        */
/* =============================================================================================== */
/* =============================================================================================== */


extern OPTPIAACMCDESIGN *piaacmc;
extern PIAACMCsimul_varType piaacmcsimul_var;





/* =============================================================================================== */
/* =============================================================================================== */
/*                                    FUNCTION(S) SOURCE CODE                                      */
/* =============================================================================================== */
/* =============================================================================================== */

// transmits between rin and rout
long PIAACMCsimul_mkSimpleLyotStop(
	const char *ID_name, 
	float rin, 
	float rout
	)
{
    long size;
    long size2;
    long ii, k;
    long ID, IDr;
    float r;

	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif


    size = piaacmc[0].size;
    size2 = size*size;

	IDr = image_ID("rcoord");
	
    ID = create_3Dimage_ID(ID_name, size, size, piaacmc[0].nblambda);
    for(k=0; k<piaacmc[0].nblambda; k++)
        for(ii=0; ii<size2; ii++)
        {
            if((data.image[IDr].array.F[ii]<rout)&&(data.image[IDr].array.F[ii]>rin))
                data.image[ID].array.F[k*size2+ii] = 1.0;
            else
                data.image[ID].array.F[k*size2+ii] = 0.0;
        }

    return(ID);
}




