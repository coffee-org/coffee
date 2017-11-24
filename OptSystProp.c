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



extern DATA data;

#define SBUFFERSIZE 2000





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
	init_OptSystProp();
//	printf(" ...... Loading module %s\n", __FILE__);
}






int_fast8_t init_OptSystProp()
{
    strcpy(data.module[data.NBmodule].name, __FILE__);
    strcpy(data.module[data.NBmodule].package, "coffee"); 
    strcpy(data.module[data.NBmodule].info, "Optical propagation through system");
    data.NBmodule++;


    // add atexit functions here
   

    return 0;

}










