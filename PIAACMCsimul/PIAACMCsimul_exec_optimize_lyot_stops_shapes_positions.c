/**
 * @file    PIAAACMCsimul_exec_optimize_lyot_stops_shapes_positions.c
 * @brief   PIAA-type coronagraph design, execute compute image
 * 
 *  
 * @author  O. Guyon
 * @date    25 nov 2017
 *
 * 
 * @bug No known bugs.
 * 
 */



// System includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"




extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTSYST *optsyst;

extern OPTPIAACMCDESIGN *piaacmc;











/** @brief Make Lyot stops using off-axis light minimums
 * 
 *  Finds minumum flux level in 3D intensity data cube
 * 
 * @param[in]  IDincohc_name             image: Input 3D intensity 
 * @param[out] test_oals_val.fits        FITS : Minimum 2D image
 * @param[out] test_oals_index.fits      FITS : z-index for each pixel of the 2D output minimum
 * 
*/

long PIAACMCsimul_optimizeLyotStop_offaxis_min(
		const char *IDincohc_name
		)
{
    long IDincohc;
    
    long IDindex;
    long IDminflux;
    long ii, jj, kk;
    long xsize = 0;
    long ysize = 0;
    long xysize = 0;
    double minv;
    long minindex;
    double tmpv;
    long NBz;
    

	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif

    
    IDincohc = image_ID(IDincohc_name);

    xsize = data.image[IDincohc].md[0].size[0];
    ysize = data.image[IDincohc].md[0].size[1];
    xysize = xsize*ysize;
    NBz = data.image[IDincohc].md[0].size[2];
    
    IDindex = create_2Dimage_ID("oals_index", xsize, ysize);
    IDminflux = create_2Dimage_ID("oals_val", xsize, ysize);
    
    for(ii=0;ii<xsize;ii++)
        for(jj=0;jj<ysize;jj++)
            {
                minv = data.image[IDincohc].array.F[jj*xsize+ii];
                minindex = 0;

                for(kk=1;kk<NBz;kk++)
                {
                    tmpv = data.image[IDincohc].array.F[kk*xysize+jj*xsize+ii];
                    if(tmpv<minv)
                        {
                            minv = tmpv;
                            minindex = kk;
                        }
                }
                data.image[IDminflux].array.F[jj*xsize+ii] = minv;
                data.image[IDindex].array.F[jj*xsize+ii] = minindex;
            }    
    
    save_fits("oals_index", "!test_oals_index.fits");
    save_fits("oals_val", "!test_oals_val.fits");
    
    return 0;
}
















/**
  * ---
  *
  * ## Mode 5: Optimize Lyot stops shapes and positions
  *
  *
  *
  */
int PIAACMCsimul_exec_optimize_lyot_stops_shapes_positions()
{
    long IDv;
    double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;
    long elem, elem0;
    long ls, NBpropstep, k;
    double oaoffset; // off-axis offset for Lyot stop optimization
    double zmin, zmax;
    char fnamea[500];
    char fnamep[500];
    double lstransm;
    long NBincpt, k1;
    long ii;
    long xsize = 0;
    long ysize = 0;
    long NBkr, kr;
    long ID1, IDa, IDc;
    char command[1000];
    char fname[1000];
    long cnt;
    char fptestname[1000];
    FILE *fptest;
    int r;

    printf("=================================== mode 005 ===================================\n");
    // load some more cli variables
    if( (IDv = variable_ID("PIAACMC_centobs0")) != -1)
        centobs0 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_centobs1")) != -1)
        centobs1 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }

    /// ### Initialize as in mode 0
    PIAACMCsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 1);
    PIAACMCsimul_makePIAAshapes(piaacmc, 0);
    optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm



    /// ### Load CLI variables as appropriate

    /// - <- **PIAACMC_nbpropstep** : # of propagation steps along the beam
    NBpropstep = 150;
    if((IDv=variable_ID("PIAACMC_nbpropstep"))!=-1)
        NBpropstep = (long) data.variable[IDv].value.f+0.01;

    /// - <- **PIAACMC_lstransm** : desired Lyot stop transmission
    lstransm = 0.85;
    if((IDv=variable_ID("PIAACMC_lstransm"))!=-1)
        lstransm = (double) data.variable[IDv].value.f;
    printf("lstransm  = %f\n", lstransm);





    /// ### Identify post focal plane pupil plane (first pupil after focal plane mask)

    /// Provides reference complex amplitude plane for downstream analysis
    PIAACMCsimul_init(piaacmc, 0, 0.0, 0.0);

    printf("=========== %ld elements ======================================================\n", optsyst[0].NBelem);
    // find the ID of the "post focal plane mask pupil" element
    for(elem=0; elem<optsyst[0].NBelem; elem++)
    {
        if( strcmp("post focal plane mask pupil", optsyst[0].name[elem]) == 0 )
        {
            elem0 = elem;
            printf("post focal plane mask pupil = %ld\n", elem);
        }
        else
            printf("elem %ld : %s\n", elem, optsyst[0].name[elem]);
    }
    optsyst[0].keepMem[elem0] = 1; // keep it for future use





    /// ### Compute incoherent 3D illumination near pupil plane

    /// Multiple off-axis sources are propagated and the corresponding intensities added
    /// - -> **OAincohc**  Output incoherent image (3D)

    oaoffset = 20.0; // off axis amplitude
    // compute the reference on-axis PSF
    PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);

    // filenames of the complex amplitude and phase in the post FPM pupil plane indexed by elem0
    sprintf(fnamea, "WFamp0_%03ld", elem0);
    sprintf(fnamep, "WFpha0_%03ld", elem0);

    printf("elem0 = %ld\n", elem0);

    // args for the PIAACMCsimul_CA2propCubeInt function

    // set range of propagation
    zmin = piaacmc[0].LyotZmin;
    zmax = piaacmc[0].LyotZmax;

    // args that determine "extended" off-axis source
    // number of off-axis sources on each circle
    NBincpt = 15;
    // number of circle radii
    NBkr = 5;

    // propagate complex amplitude in a range from zmin to zmax, where 0 is elem0
    // computes the diffracted light from the on-axis source
    ID1 = PIAACMCsimul_CA2propCubeInt(fnamea, fnamep, zmin, zmax, NBpropstep, "iproptmp");
    // complex amplitude at elem0, only used to determine image size
    IDa = image_ID(fnamea);

    xsize = data.image[IDa].md[0].size[0];
    ysize = data.image[IDa].md[0].size[1];

    /// OAincohc is the summed light "all" from off-axis sources in the pupil,
    /// including the on-axis source(!),
    /// giving intensity contribution of all off-axis sources
    /// in order to preserve the intensity of the off-axis in the design.
    /// load OAincohc if exist, maybe we've been here before
    sprintf(fname, "%s/OAincohc.fits", piaacmcsimul_var.piaacmcconfdir);
    IDc = load_fits(fname, "OAincohc", 1);


    if(IDc==-1) // OAincohc does not exist so we have to make it
    {
        // create image to receive sum
        IDc = create_3Dimage_ID("OAincohc", xsize, ysize, NBpropstep);

        cnt = 0; // initialize counter so later we can normalize by number of sources
        // loop over radii
        for(kr=0; kr<NBkr; kr++)
        {
            // loop over points at current radius
            for(k1=0; k1<NBincpt; k1++)
            {
                // compute PSF for a point at this angle with scaled offset
                // PIAACMCsimul_computePSF changes fnamea and fnamep (in call to OptSystProp_run)!
                PIAACMCsimul_computePSF(oaoffset*(1.0+kr)/NBkr*cos(2.0*M_PI*k1/NBincpt), oaoffset*(1.0+kr)/NBkr*sin(2.0*M_PI*k1/NBincpt), 0, optsyst[0].NBelem, 0, 0, 0, 0);
                // propagate that elem0 from zmin to zmax with new PSF
                ID1 = PIAACMCsimul_CA2propCubeInt(fnamea, fnamep, zmin, zmax, NBpropstep, "iproptmp");
                for(ii=0; ii<xsize*ysize; ii++) // ii is indexing x-y plane
                    for(k=0; k<NBpropstep; k++) // k is indexing z-direction. adding to IDc which looks right
                        data.image[IDc].array.F[k*xsize*ysize+ii] += data.image[ID1].array.F[k*xsize*ysize+ii];
                delete_image_ID("iproptmp");
                cnt ++;
            }
        }
        // scale by the number of sources to give average
        for(ii=0; ii<xsize*ysize; ii++)
            for(k=0; k<NBpropstep; k++)
                data.image[IDc].array.F[k*xsize*ysize+ii] /= cnt;


        sprintf(fname, "!%s/OAincohc.fits", piaacmcsimul_var.piaacmcconfdir);
        // save the final result
        save_fits("OAincohc", fname);
    }



    /// ### Compute on-axis PSF 3D intensity to define light to reject
    PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);
    // propagate it into the optical system, with result in image named "iprop00"

    /// call function PIAACMCsimul_CA2propCubeInt() to compute 3D intensity cube
    /// @param[out] iprop00 3D intensity image
    PIAACMCsimul_CA2propCubeInt(fnamea, fnamep, zmin, zmax, NBpropstep, "iprop00");
    //  save_fits("iprop00", "!test_iprop00.fits");

    /// ### Compute image that has the min along z of OAincohc at each x,y
    /// Function PIAACMCsimul_optimizeLyotStop_offaxis_min() computes minimal intensity image.
    ///
    PIAACMCsimul_optimizeLyotStop_offaxis_min("OAincohc");

    /// ### Optimize Lyot stops

    /// The actual Lyot stop shape and location optimization is done by PIAACMCsimul_optimizeLyotStop(), producing optimal Lyot stops in optLM*.fits
    /// and position relative to elem0 in piaacmc[0].LyotStop_zpos
    PIAACMCsimul_optimizeLyotStop(fnamea, fnamep, "OAincohc", zmin, zmax, lstransm, NBpropstep, piaacmc[0].NBLyotStop);

    sprintf(fptestname, "conj_test.txt");
    fptest = fopen(fptestname, "w");
    fprintf(fptest, "# Optimal Lyot stop conjugations\n");
    fprintf(fptest, "# \n");
    fprintf(fptest, "# DO NOT EDIT THIS FILE\n");
    fprintf(fptest, "# Written by %s in %s\n", __FUNCTION__, __FILE__);
    fprintf(fptest, "# \n");
    fprintf(fptest, "# Lyot stop index   zmin   zmax    LyotStop_zpos    elemZpos[elem0]\n");
    for(ls=0; ls<piaacmc[0].NBLyotStop; ls++)
        fprintf(fptest, "%5ld  %f  %f     %f  %f\n", ls, zmin, zmax, piaacmc[0].LyotStop_zpos[ls], optsyst[0].elemZpos[elem0]);
    fclose(fptest);

    // convert Lyot stop position from relative to elem0 to absolute
    for(ls=0; ls<piaacmc[0].NBLyotStop; ls++)
        piaacmc[0].LyotStop_zpos[ls] += optsyst[0].elemZpos[elem0];

    // and we're done!  save.
    PIAACMCsimul_savepiaacmcconf( piaacmcsimul_var.piaacmcconfdir);
    // copy to the final Lyot stop file for this mode
    for(ls=0; ls<piaacmc[0].NBLyotStop; ls++)
    {
        sprintf(command, "cp ./%s/optLM%02ld.fits ./%s/LyotStop%ld.fits", piaacmcsimul_var.piaacmcconfdir, ls, piaacmcsimul_var.piaacmcconfdir, ls);
        r = system(command);
    }

    return 0;
}





