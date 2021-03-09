/**
 * @file    PIAACMCsimul_load2DRadialApodization.c
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
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_arith/COREMOD_arith.h"

#include "info/info.h"
#include "linopt_imtools/linopt_imtools.h"

#include "PIAACMCsimul/PIAACMCsimul.h"




extern PIAACMCsimul_varType piaacmcsimul_var;





//
// load and fit radial apodization profile
// modal basis is mk(r) : cos(r*k*M_PI/1.3)
//
uint_fast8_t PIAACMCsimul_load2DRadialApodization(
		const char *IDapo_name, 
		float beamradpix, 
		const char *IDapofit_name
		)
{
    long NBpts;
    long IDm;
    long sizem;
    long kmax = 10;
    long ID, IDmask, IDin;
    long ii, jj;
    long offset;
    long sizein;
    float eps = 1.0e-4;
    char fname[500];
    int ret;
    char command[1000];
    int debug = 0;

	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif


    sizem = (long) (beamradpix*2);

    // CREATE MODES IF THEY DO NOT EXIST
    if((IDm=image_ID("APOmodesCos"))==-1)
    {
        IDm = linopt_imtools_makeCosRadModes("APOmodesCos", sizem, kmax, ApoFitCosFact*beamradpix, 1.0);
        sprintf(fname, "!%s/APOmodesCos.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("APOmodesCos", fname);
    }

    // CREATE MASK AND CROP INPUT
    IDmask = create_2Dimage_ID("fitmaskapo", sizem, sizem);

    IDin = image_ID(IDapo_name);
    sizein = data.image[IDin].md[0].size[0];
    ID = create_2Dimage_ID("_apoincrop", sizem, sizem);
    offset = (sizein-sizem)/2;
    for(ii=0; ii<sizem; ii++)
        for(jj=0; jj<sizem; jj++)
        {
            data.image[ID].array.F[jj*sizem+ii] = data.image[IDin].array.F[(jj+offset)*sizein+(ii+offset)];
            if((data.image[ID].array.F[jj*sizem+ii]>eps)&&(ii%1==0)&&(jj%1==0))
                data.image[IDmask].array.F[jj*sizem+ii] = 1.0;
        }

    if(debug==1)
    {
        sprintf(fname, "!%s/_apoincrop.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("_apoincrop", fname);

        sprintf(fname, "!%s/fitmaskapo.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("fitmaskapo", fname);
    }

    linopt_imtools_image_fitModes("_apoincrop", "APOmodesCos", "fitmaskapo", 1.0e-8, IDapofit_name, 0);
    sprintf(command, "mv %s/eigenv.dat %s/eigenv_APOmodesCos.dat", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.piaacmcconfdir);
    ret = system(command);

    if(debug==1) // test fit quality
    {
        linopt_imtools_image_construct("APOmodesCos", IDapofit_name, "testapofitsol");

        sprintf(fname, "!%s/testapofitsol.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("testapofitsol", fname);

        arith_image_sub("_apoincrop", "testapofitsol", "apofitres");
        arith_image_mult("apofitres", "fitmaskapo", "apofitresm");

        sprintf(fname, "!%s/apofitres.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("apofitres", fname);

        sprintf(fname, "!%s/apofitresm.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("apofitresm", fname);

        // linopt_imtools_image_fitModes("apofitres", "APOmodesCos", "fitmaskapo", 1.0e-5, "test2c", 0);
        info_image_stats("apofitresm", "");
    }

    delete_image_ID("_apoincrop");
    delete_image_ID("fitmaskapo");


    return 0;
}







