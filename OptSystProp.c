/**
 * @file    OptSystProp.c
 * @brief   Optical system propagation
 * 
 * Propagate beam through optical systsm
 *  
 * @author  O. Guyon
 * @date    7 Jul 2017
 *
 * 
 * @bug No known bugs.
 * 
 */



#include <string.h>


#include "CommandLineInterface/CLIcore.h"

#include "OptSystProp/OptSystProp.h"




#define SBUFFERSIZE 2000

static int INITSTATUS_OptSystProp = 0;



// CLI commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string
// 4: existing image
//





void __attribute__ ((constructor)) libinit_OptSystProp()
{
	if ( INITSTATUS_OptSystProp == 0 )
	{
		init_OptSystProp();
		RegisterModule(__FILE__, "coffee", "Optical propagation through system");
		INITSTATUS_OptSystProp = 1; 
	}
}






int_fast8_t init_OptSystProp()
{
    // add atexit functions here
   

    return 0;

}










