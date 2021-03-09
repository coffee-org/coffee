/**
 * @file    PIAACMCsimul_PIAACMC_FPMresp_rmzones.c
 * @brief   PIAA-type coronagraph design
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


// milk includes
#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"










// remove outer zones to FPMresp
long PIAACMC_FPMresp_rmzones(const char *FPMresp_in_name, const char *FPMresp_out_name, long NBzones)
{
    long ID, IDout;
    long ii, jj, kk;
    long xsize = 0;
    long ysize = 0;
    long zsize = 0;
    long ysize1 = 0;

	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif


    ID = image_ID(FPMresp_in_name);
    xsize = data.image[ID].md[0].size[0];
    ysize = data.image[ID].md[0].size[1];
    zsize = data.image[ID].md[0].size[2];


    ysize1 = data.image[ID].md[0].size[1]-NBzones;

    IDout = create_3Dimage_ID_double(FPMresp_out_name, xsize, ysize1, zsize);

    for(ii=0; ii<xsize; ii++)
        for(kk=0; kk<zsize; kk++)
        {
            for(jj=0; jj<ysize1; jj++)
                data.image[IDout].array.D[kk*xsize*ysize1 + jj*xsize + ii] = data.image[ID].array.D[kk*xsize*ysize + jj*xsize + ii];
            for(jj=ysize1; jj<ysize; jj++)
                data.image[IDout].array.D[kk*xsize*ysize1 + ii] += data.image[ID].array.D[kk*xsize*ysize + jj*xsize + ii];
        }

    return(IDout);
}


