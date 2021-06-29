/**
 * @file    PIAACMCsimul_PIAACMC_FPMresp_rmzones.c
 * @brief   PIAA-type coronagraph design
 *
 * Can design both APLCMC and PIAACMC coronagraphs
 *
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
    DEBUG_TRACE_FSTART();

    long ID, IDout;

#ifdef PIAASIMUL_LOGFUNC0
    PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
#endif


    ID = image_ID(FPMresp_in_name);
    uint32_t xsize = data.image[ID].md[0].size[0];
    uint32_t ysize = data.image[ID].md[0].size[1];
    uint32_t zsize = data.image[ID].md[0].size[2];


    long ysize1 = data.image[ID].md[0].size[1]-NBzones;

    create_3Dimage_ID_double(FPMresp_out_name, xsize, ysize1, zsize, &IDout);

    for(uint32_t ii=0; ii<xsize; ii++)
        for(uint32_t kk=0; kk<zsize; kk++)
        {
            for(long jj=0; jj<ysize1; jj++)
                data.image[IDout].array.D[kk*xsize*ysize1 + jj*xsize + ii] = data.image[ID].array.D[kk*xsize*ysize + jj*xsize + ii];
            for(long jj=ysize1; jj<ysize; jj++)
                data.image[IDout].array.D[kk*xsize*ysize1 + ii] += data.image[ID].array.D[kk*xsize*ysize + jj*xsize + ii];
        }

    DEBUG_TRACE_FEXIT();
    return(IDout);
}


