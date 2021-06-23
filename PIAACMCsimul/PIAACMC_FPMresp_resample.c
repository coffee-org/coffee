/**
 * @file    PIAACMCsimul_FPMresp_resample.c
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







// compress FPMresp to smaller number of lambda and points
long PIAACMC_FPMresp_resample(
    const char *FPMresp_in_name,
    const char *FPMresp_out_name,
    long NBlambda,
    long PTstep
)
{
    DEBUG_TRACE_FSTART();

    imageID ID = -1;
    imageID IDout = -1;

#ifdef PIAASIMUL_LOGFUNC0
    PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__,
                                 "");
#endif


    ID = image_ID(FPMresp_in_name);
    uint32_t xsize = data.image[ID].md[0].size[0];
    uint32_t ysize = data.image[ID].md[0].size[1];
    uint32_t zsize = data.image[ID].md[0].size[2];


    long xsize1 = (long)(xsize / PTstep);
    long zsize1 = NBlambda;


    IDout = create_3Dimage_ID_double(FPMresp_out_name, xsize1, ysize, zsize1);
    /*	for(kk1=0;kk1<zsize1;kk1++)
    		for(ii=0;ii<xsize;ii++)
    			data.image[IDout].array.D[kk1*xsize1*ysize + ii] = 1;
    	*/



    for(long kk1 = 0; kk1 < zsize1; kk1++)
    {
        double kk1x = 1.0 * kk1 / (zsize1 - 1);
        long kk1xi = (long) (kk1x * (zsize - 1));
        double alpha = kk1x * (zsize - 1) - kk1xi;

        if(kk1xi == zsize - 1)
        {
            kk1xi = zsize - 2;
            alpha = 1.0;
        }

        printf("lambda index %ld  (%lf)   :  %ld %lf (%lf)\n", kk1, kk1x, kk1xi, alpha,
               (kk1xi + alpha) / (zsize - 1));
        long ii1 = 0;



        for(uint32_t ii = 0; ii < xsize; ii += 2 * PTstep)
        {

            for(uint32_t jj = 0; jj < ysize; jj++)
            {
                double re0 = data.image[ID].array.D[kk1xi * xsize * ysize + jj * xsize + ii];
                double im0 = data.image[ID].array.D[kk1xi * xsize * ysize + jj * xsize + ii + 1];

                double re1 = data.image[ID].array.D[(kk1xi + 1) * xsize * ysize + jj * xsize + ii];
                double im1 = data.image[ID].array.D[(kk1xi + 1) * xsize * ysize + jj * xsize + ii + 1];

                double re = (1.0 - alpha) * re0 + alpha * re1;
                double im = (1.0 - alpha) * im0 + alpha * im1;


                data.image[IDout].array.D[kk1 * xsize1 * ysize + jj * xsize1 + ii1] =
                    re; //*sqrt(PTstep);
                data.image[IDout].array.D[kk1 * xsize1 * ysize + jj * xsize1 + ii1 + 1] =
                    im; //*sqrt(PTstep);
            }

            ii1 += 2;
        }
    }


    /*		for(kk=0;kk<zsize;kk++)
    			{
    				for(jj=0;jj<ysize1;jj++)
    					data.image[IDout].array.D[kk*xsize*ysize1 + jj*xsize + ii] = data.image[ID].array.D[kk*xsize*ysize + jj*xsize + ii];
    				for(jj=ysize1;jj<ysize;jj++)
    					data.image[IDout].array.D[kk*xsize*ysize1 + ii] += data.image[ID].array.D[kk*xsize*ysize + jj*xsize + ii];
    			}
    */

    DEBUG_TRACE_FEXIT();
    return(IDout);
}
