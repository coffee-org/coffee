/**
 * @file    PIAAACMCsimul_free.c
 * @brief   PIAA-type coronagraph design, initialize
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



// System includes
#include <stdio.h>
#include <stdint.h>
#include <malloc.h>
#include <math.h>

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"



extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTSYST *optsyst;



/**
 * Frees memory for module
 */
void PIAACMCsimul_free( void )
{
	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "\n");
	#endif

	
    if( piaacmcsimul_var.optsystinit == 1 )
    {
        free(optsyst);
    }
}
