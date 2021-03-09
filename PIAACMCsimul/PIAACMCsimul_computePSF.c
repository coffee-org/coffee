/**
 * @file    PIAACMCsimul_computePSF.c
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
#include <math.h>


// milk includes
#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"


#include "linopt_imtools/linopt_imtools.h"
#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"




extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTSYST *optsyst;

extern OPTPIAACMCDESIGN *piaacmc;






/**
 * @brief Compute PSF
 * 
 * @return Average contrast in evaluation zone
 * 
 * 
 *  
 * Source is defined by parameters sourcesize and extmode :
 * - source size = 1e-{sourcesize*0.1}, except if sourcesize = 0 (point source)
 * - sourcesize is a 2-digit number ( 10 = 0.1 l/D, 20 = 0.01 l/D etc..)
 * - extmode = 0 : 1 point (point source)
 * - extmode = 1 : 3 point sources, 120 apart on circle radius = source size
 * - extmode = 2 : 6 point sources. 3 as above on circle radius 1/sqrt(2.5) + 3 on outer circle, radius 2/sqrt(2.5), 120 apart, clockled 60 deg off inner points
 * 
 * @note If opderrcube exists, include each slice as a WF mode  
 * 
 * PSF is held in shared memory by default
 * 
 * ---
 * 
 * ### Output
 * 
 * | name                              |  type      | Description                            |
 * |-----------------------------------|------------|----------------------------------------|
 * | scoringmask                       | 2D image   | focal plane points used for evaluation |
 * | <piaacmcdir>/scoringmask<N>.fits  | 2D FITS    | focal plane points used for evaluation |
 * | imvec                             | 1D image   | output vector                          |
 * | psfi0                             | 3D image   | output PSF                             |
 * 
 * 
 * ---
 * ---
 * 
 * 
 */ 



double PIAACMCsimul_computePSF(
		float xld, 			/// @param[in]   xld         float: Source X position [l/D]
		float yld, 			/// @param[in]   yld         float: Source Y position [l/D]
		long startelem, 	/// @param[in]   startelem   long : First element in propagation
		long endelem, 		/// @param[in]   endelem     long : Last element in propagation
		int savepsf, 		/// @param[in]   savepsf     int  : Save PSF flag
		int sourcesize, 	/// @param[in]   sourcezise  int  : Source size (10x log10)
		int extmode, 		/// @param[in]   extmode     int  : Source extended type
		int outsave			/// @param[in]   outsave     int  : Save output flag
		)
{
    FILE *fp;
    FILE *fpflux;
    double x, y;
    long IDa, IDp;
    long size;
    long nblambda;
    long size2;
    long ii, jj, k;
    long IDpiaa1z, IDpiaa2z;
    long elem;
    long kl;

    char fname_piaa1z[500];
    char fname_piaa2z[500];
    char fname_pupa0[500];
    char fname_pupp0[500];
    char fname[1500];

    long ID;
    long index;

    double proplim = 1.0e-4;
    double total;

    long size0, size1;
    long Cmsize, Fmsize;


    // how to measure quality
    float focscale; // l/D per pix
    float scoringIWA = 1.5; 
    float scoringOWA = 20.0;
    float scoringOWAhr = 8.0;
    float scoringIWAx = -20.5;
    long IDsm;
    float r;

    double value;
    double avContrast;
    double peakcontrast;
    double tmpv;
    double value1;


    double dld;
    long nbelem;
    long ID1;
    long offset;
    long offset1;
    double normcoeff = 1.0;
    double pha;

    int ret;
    char command[1000];
    double rad1, rad2;

    float val;
    long IDv;

	long nbOPDerr = 0;
	long OPDmode;
	long IDopderrC = -1;
	long IDopderr;
	
	long imindex = 0;
	long NBimindex = 0;
	char imname[200];
	
	long naxis;
	

	#ifdef PIAASIMUL_LOGFUNC1
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif


    // size of one side of each image array
    size = piaacmc[0].size;
    
    // number of pixels in each image
    size2 = size*size;

	//printf("Loading (optional) OPDerr file\n");
	//fflush(stdout);

	// load an error if it exists
	IDopderrC = image_ID("OPDerrC");
	if(IDopderrC == -1)
		IDopderrC = load_fits("OPDerrC.fits", "OPDerrC", 0);
		
	
		

	if(IDopderrC != -1)
		{
			naxis = data.image[IDopderrC].md[0].naxis;
			if(naxis==2)
				nbOPDerr = data.image[IDopderrC].md[0].size[2];  // number of error arrays
			else
				nbOPDerr = 0;
			printf("INCLUDING %ld OPD ERROR MODES\n", nbOPDerr);
			fflush(stdout);
		}
	else
		{
			//printf("NO OPD ERROR MODES\n");
			//fflush(stdout);
			nbOPDerr = 0;
		}
    // focal plane plate scale in lambda/D per pixel
    focscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size;


    /// ### Create scoring mask if it doesn't exist
    /// The scoring mask is the array of evaluation points on the focal plane
    if((IDsm=image_ID("scoringmask"))==-1)
    {
		printf("CREATING SCORING MASK\n");
        printf("FOCAL PLANE SCALE = %f l/d per pix\n", focscale);
		fflush(stdout);
        IDsm = create_2Dimage_ID("scoringmask", size, size);

        if(piaacmcsimul_var.SCORINGMASKTYPE==0) // high density, wide
        {
            // draw an array of points, clip to desired subregions
            for(ii=0; ii<size; ii++)
                for(jj=0; jj<size; jj++)
                {   // create regular array and the raduis of each point
                    x = (1.0*ii-0.5*size)*focscale;
                    y = (1.0*jj-0.5*size)*focscale;
                    r = sqrt(x*x+y*y);

                    // clip the regular array to the desired annulus inside high-resolution
                    // region defined by scoringOWAhr
                    // use every other point
                    // and clip to an x > scoringIWAx part of the annulus if desired
                    if((r>scoringIWA)&&(r<scoringOWAhr)&&(x>scoringIWAx)&&((ii+jj)%2==0))
                        data.image[IDsm].array.F[jj*size+ii] = 1.0;
                    // pick every other row and column between scoringOWAhr and scoringOWA
                    if((r>scoringOWAhr)&&(r<scoringOWA)&&(x>scoringIWAx)&&(ii%2==0)&&(jj%2==0))
                        data.image[IDsm].array.F[jj*size+ii] = 1.0;
                    // draw a single radial line of points out to IWA = 70 (every other point)
                    if((x>scoringOWA)&&(fabs(y)<scoringIWA*0.25)&&(r<70.0)&&((ii+jj)%2==0)) // single line
                        data.image[IDsm].array.F[jj*size+ii] = 1.0;
                }
        }
        else // focused on central pixels, fewer pixels for faster convergence - used for FPMresp based optimization
        {
            for(ii=0; ii<size; ii++)
                for(jj=0; jj<size; jj++)
                {
                    x = (1.0*ii-0.5*size)*focscale;
                    y = (1.0*jj-0.5*size)*focscale;
                    r = sqrt(x*x+y*y);
                    // clip from scoringIWA to scoringOWAhr only using every other column and row
                    if((r>scoringIWA)&&(r<scoringOWAhr)&&(x>scoringIWAx)&&(ii%2==0)&&(jj%2==0))
                        data.image[IDsm].array.F[jj*size+ii] = 1.0;
                }
        }
        if( piaacmcsimul_var.PIAACMC_save==1)
        {
            // save a disgnostic image
            sprintf(fname, "!%s/scoringmask%d.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.SCORINGMASKTYPE);
            save_fits("scoringmask", fname);
        }
        // a pixtable is a list of non-zero pixels with their coordinates
        linopt_imtools_mask_to_pixtable("scoringmask", "pixindex", "pixmult");
        // sums the image, giving the total number of pixels in the scoring mask
        piaacmcsimul_var.SCORINGTOTAL = arith_image_total("scoringmask");

        //exit(0);

    }


	
	
	
	/// ## Fast PSF computattion (if piaacmcsimul_var.computePSF_FAST_FPMresp = 1)
	
    if(piaacmcsimul_var.computePSF_FAST_FPMresp==1) /// @note Only possible if mode 11 has already been executed
    {
        /// Compute the PSF as the complex amplitude for the evaluation points on the focal plane
        /// for a given FPM zone thickness based on the FPMresp array computed in mode 11
        value1 = PIAACMCsimul_achromFPMsol_eval(piaacmcsimul_var.fpmresp_array, piaacmcsimul_var.zonez_array, piaacmcsimul_var.dphadz_array, piaacmcsimul_var.outtmp_array, piaacmcsimul_var.vsize, data.image[piaacmc[0].zonezID].md[0].size[0], optsyst[0].nblambda);

//		printf("FAST FPMresp CALLED\n");
//		sleep(1000000);//TEST
		
	
        ///
        /// PSF result is stored in outtmp_array

        value = 0.0;
        peakcontrast = 0.0;

        ID = image_ID("imvect"); /// - Use \c imvect for storage if it exists, or create it
        if(ID==-1)
            ID = create_2Dimage_ID("imvect", piaacmcsimul_var.vsize*optsyst[0].nblambda, 1);
        /// - Write the result into \c imvect 
        for(ii=0; ii< piaacmcsimul_var.vsize*optsyst[0].nblambda; ii++) // for each wavelength
        {
            data.image[ID].array.F[ii] = piaacmcsimul_var.outtmp_array[ii];
            // square to give intensity
            tmpv = piaacmcsimul_var.outtmp_array[ii]*piaacmcsimul_var.outtmp_array[ii];
            // total intensity = sum(intensity_i) = sum(Re_i^2 + Im_i^2)
            value += tmpv;
        }
        /// - Total flux in the output vector is stored in \c piaacmcsimul_var.PIAACMCSIMUL_VAL0 as total flux
        piaacmcsimul_var.PIAACMCSIMUL_VAL0 = value;

        /// set \c value to average value per area normalized to flux
        value = value/size/size/optsyst[0].flux[0]; // flux[0] is proportional to the number of lambda channels, so this normalization makes value independant of number of spectral channels
        // here value is the total light (averaged across spectral channels) in the measurement points, normalized to the input flux
        // actual average contrast, averaging over # of pixels and physical area
        avContrast = value/(piaacmcsimul_var.SCORINGTOTAL*focscale*focscale);


        //        printf("Peak constrast (rough estimate)= %g\n", peakcontrast/size/size/optsyst[0].flux[0]/focscale/focscale*3.0);
        //        printf("value1 = %g\n", value1);
        //		printf("Total light in scoring field = %g  -> Average contrast = %g   (%g)\n", value, value/(arith_image_total("scoringmask")*focscale*focscale), value1/piaacmcsimul_var.CnormFactor/optsyst[0].nblambda);

    }
    else /// ## Full/Slow PSF computation (if piaacmcsimul_var.computePSF_FAST_FPMresp = 0)
    {
        /// The PSF for an extended source is approximated as
        /// a collection of point sources.
        /// Sourcesize determines the separation of the point sources
        if( sourcesize != 0 ) // sourcesize > 0 only if in linear optimization (step >= 100)
        {
            printf("COMPUTING RESOLVED SOURCE PSF / ADDING OPD MODES\n");	
            fflush(stdout);

            // dld is the radius of the circle containing the point sources in lamba/D
            dld = 1.0/pow(10.0, 0.1*sourcesize); // nominal pointing offset [l/D]
            // extmode controls how many point sources we use to model the extended source
			// extmode == 0 => single point (unresolved source)
            // extmode == 1 => three point sources
            // extmode == 2 => six point sources
            if (extmode==2)
            {
                // if I have six sources, put them in two rings of three each
                rad1 = dld/sqrt(2.5);
                rad2 = 2.0*dld/sqrt(2.5);
            }
            else
            {
                // if I have three sources keep them at the same radii
                rad1 = dld;
                rad2 = dld;
            }
            
            // we will collect propagation results in a set of vectors called imvectp<1,2,...>
            // which will ultimately be collected into a single long imvect

            // image index, counts the number of PSFs we make, one for each point source
			imindex = 0;
            // name of the image vector to hold the propagation result
			sprintf(imname, "imvectp%02ld", imindex);
			
            // xld and yld are the input positions of the input source in lambda/D
			// initialize first point source, which sets optsyst
			if (extmode==0)
				PIAACMCsimul_init(piaacmc, 0, xld, yld);
			else
				PIAACMCsimul_init(piaacmc, 0, xld+rad1, yld);	
			
			PIAACMCsimul_makePIAAshapes(piaacmc, 0);
			
            // propagate it (optsyst is a global), output in psfc0 (complex amlitude)
            // and psfi0 (intensity)
            OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcsimul_var.piaacmcconfdir, 0);


            // convert the image to the vector
            linopt_imtools_Image_to_vec("psfc0", "pixindex", "pixmult", imname);
            // save the intensity of the first point
            copy_image_ID("psfi0", "psfi0ext", 0);
			//sprintf(fname, "!%s/psfi0_pt%03ld.fits", piaacmcsimul_var.piaacmcconfdir, imindex);
			//save_fits("psfi0", fname);
			imindex++;            
            
            if (extmode>0)
            {
            // do the same for the second point
			sprintf(imname, "imvectp%02ld", imindex);
            pha = 2.0*M_PI/3.0; // 1/3 around the circle
            PIAACMCsimul_init(piaacmc, 0, xld+rad1*cos(pha), yld+rad1*sin(pha)); 
            PIAACMCsimul_makePIAAshapes(piaacmc, 0);
            OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcsimul_var.piaacmcconfdir, 0);
            linopt_imtools_Image_to_vec("psfc0", "pixindex", "pixmult", imname);
            // add the intensity to build up PSF for extended source
            arith_image_add_inplace("psfi0ext", "psfi0");
			//sprintf(fname, "!%s/psfi0_pt%03ld.fits", piaacmcsimul_var.piaacmcconfdir, imindex);
			//save_fits("psfi0", fname);
			imindex++;
						
            // do the same for the third point
			sprintf(imname, "imvectp%02ld", imindex);
            pha = 4.0*M_PI/3.0; // 2/3 around the circle
            PIAACMCsimul_init(piaacmc, 0, xld+rad1*cos(pha), yld+rad1*sin(pha));
            PIAACMCsimul_makePIAAshapes(piaacmc, 0);
            OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcsimul_var.piaacmcconfdir, 0);
            linopt_imtools_Image_to_vec("psfc0", "pixindex", "pixmult", imname);
            // add the intensity to build up PSF for extended source
            arith_image_add_inplace("psfi0ext", "psfi0");
       		//sprintf(fname, "!%s/psfi0_pt%03ld.fits", piaacmcsimul_var.piaacmcconfdir, imindex);
			//save_fits("psfi0", fname);
			imindex++;
			}
            
            
            if (extmode==2)
            { // keep going for the other three points if desired, on the outer radius
				sprintf(imname, "imvectp%02ld", imindex);
                pha = M_PI/3.0;
                PIAACMCsimul_init(piaacmc, 0, xld+rad2*cos(pha), yld+rad2*sin(pha));
                PIAACMCsimul_makePIAAshapes(piaacmc, 0);
                OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcsimul_var.piaacmcconfdir, 0);
                linopt_imtools_Image_to_vec("psfc0", "pixindex", "pixmult", imname);
                arith_image_add_inplace("psfi0ext","psfi0");
				//sprintf(fname, "!%s/psfi0_pt%03ld.fits", piaacmcsimul_var.piaacmcconfdir, imindex);
				//save_fits("psfi0", fname);
				imindex++;
			

				sprintf(imname, "imvectp%02ld", imindex);
                pha = 2.0*M_PI/3.0 + M_PI/3.0;
                PIAACMCsimul_init(piaacmc, 0, xld+rad2*cos(pha), yld+rad2*sin(pha));
                PIAACMCsimul_makePIAAshapes(piaacmc, 0);
                OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcsimul_var.piaacmcconfdir, 0);
                linopt_imtools_Image_to_vec("psfc0", "pixindex", "pixmult", imname);
                arith_image_add_inplace("psfi0ext", "psfi0");
				//sprintf(fname, "!%s/psfi0_pt%03ld.fits", piaacmcsimul_var.piaacmcconfdir, imindex);
				//save_fits("psfi0", fname);
				imindex++;

				sprintf(imname, "imvectp%02ld", imindex);
                pha = 4.0*M_PI/3.0 + M_PI/3.0;
                PIAACMCsimul_init(piaacmc, 0, xld+rad2*cos(pha), yld+rad2*sin(pha));
                PIAACMCsimul_makePIAAshapes(piaacmc, 0);
                OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcsimul_var.piaacmcconfdir, 0);
                linopt_imtools_Image_to_vec("psfc0", "pixindex", "pixmult", imname);
                arith_image_add_inplace("psfi0ext", "psfi0");
				//sprintf(fname, "!%s/psfi0_pt%03ld.fits", piaacmcsimul_var.piaacmcconfdir, imindex);
				//save_fits("psfi0", fname);
				imindex++;
			
                // but multiply by 0.5 'cause we have twice as many points
				//    arith_image_cstmult_inplace("psfi0ext", 0.5);
            }
            
            
			/// ### OPTIONAL: Add OPD error to list of modes
			
			printf("Adding optional OPD error modes (%ld modes)\n", nbOPDerr);
			fflush(stdout);
			// add error modes if any
            // add new evaluation points for the error to imvect so we minimize the error
            // as well as the non-error
			for(OPDmode=0; OPDmode < nbOPDerr; OPDmode++)
			{
				sprintf(imname, "imvectp%02ld", imindex);

				IDopderr = create_2Dimage_ID("opderr", size, size);
                // "opderr" is a standard name read by PIAACMCsimul_init
				for(ii=0;ii<size*size;ii++)
					data.image[IDopderr].array.F[ii] = data.image[IDopderrC].array.F[size*size*OPDmode + ii];
				PIAACMCsimul_init(piaacmc, 0, 0.0, 0.0); // add error to the data
				PIAACMCsimul_makePIAAshapes(piaacmc, 0);
                OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcsimul_var.piaacmcconfdir, 0);
                linopt_imtools_Image_to_vec("psfc0", "pixindex", "pixmult", imname);
                arith_image_add_inplace("psfi0ext", "psfi0");
	       		//sprintf(fname, "!%s/psfi0_pt%03ld.fits", piaacmcsimul_var.piaacmcconfdir, imindex); //TEST
				//save_fits("psfi0", fname); //TEST
				delete_image_ID("opderr");
		
				imindex++;
			}
			
            /// - Average over all the PSFs we've created to simulate this extended source
			NBimindex = imindex;
			arith_image_cstmult_inplace("psfi0ext", 1.0/NBimindex);

			/// - If \c outsave = 1, save PSF to FITS file
            if(outsave==1)
            {
				PIAACMCsimul_update_fnamedescr();
                sprintf(fname, "!%s/psfi0_exsrc%03d_sm%d.%s.fits", piaacmcsimul_var.piaacmcconfdir, sourcesize, piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.fnamedescr);                    
                save_fits("psfi0ext", fname);
            }
            
      

            // get the number of elements in a single-PSF vector
            ID = image_ID("imvectp00");
            nbelem = data.image[ID].md[0].nelement;
            // make big vector to collect the complex amplitudes of all the above PSFs
            ID = image_ID("imvect");
            if(ID!=-1)
                delete_image_ID("imvect");

            offset = nbelem/piaacmc[0].nblambda; // number of pixels per lambda x2 (re, im)
            printf("offset = %ld\n", offset);

            // number of pixels per lambda x 2 times the number of PSFs
            offset1 = NBimindex*offset;
            normcoeff = 1.0/sqrt(NBimindex);

            // make an imvect for each lambda
            ID = create_2Dimage_ID("imvect", offset1, piaacmc[0].nblambda);


            // fill in with each imvectp* created above
			for(imindex=0; imindex<NBimindex; imindex++)
			{
				sprintf(imname, "imvectp%02ld", imindex);
				ID1 = image_ID(imname);
				for(kl=0; kl<piaacmc[0].nblambda; kl++)
					for(ii=0; ii<offset; ii++)
						data.image[ID].array.F[kl*offset1 + imindex*offset + ii] = data.image[ID1].array.F[kl*offset+ii]*normcoeff;
				delete_image_ID(imname);
			}


            //linopt_imtools_Image_to_vec("psfc0", "pixindex", "pixmult", "imvect");
            
            // measure the contrast for all aimplitudes in imvect
            value = 0.0;
            peakcontrast = 0.0;
            ID = image_ID("imvect");
            for(ii=0; ii<data.image[ID].md[0].nelement; ii++)
            {
                tmpv = data.image[ID].array.F[ii]*data.image[ID].array.F[ii];
                value += tmpv;
                if(tmpv>peakcontrast)
                    peakcontrast = tmpv;
            }

            for(elem=0; elem<optsyst[0].NBelem; elem++)
                printf("    FLUX %3ld   %12.4lf %8.6lf\n", elem, optsyst[0].flux[elem], optsyst[0].flux[elem]/optsyst[0].flux[0]);


			/// - If \c outsave = 1, save flux to txt file
            if(outsave==1)
            {				
                PIAACMCsimul_update_fnamedescr();
                sprintf(fname, "!%s/flux_exsrc%03d_sm%d.%s.txt", piaacmcsimul_var.piaacmcconfdir, sourcesize, piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.fnamedescr);                
                fpflux = fopen(fname, "w");
                for(elem=0; elem<optsyst[0].NBelem; elem++)
                    fprintf(fpflux, "%18.16lf %18.16lf  %d\n", optsyst[0].flux[elem], optsyst[0].flux[elem]/optsyst[0].flux[0], optsyst[0].nblambda);
                fprintf(fpflux, "W0  %d\n", optsyst[0].nblambda);
                fclose(fpflux);
                
                sprintf(command, "cp %s %s/tmp_flux.txt", fname, piaacmcsimul_var.piaacmcconfdir);
                ret = system(command);
            }



            // compute average contrast
            value = value/size/size/optsyst[0].flux[0];
            avContrast = value/(piaacmcsimul_var.SCORINGTOTAL*focscale*focscale);

            //           piaacmcsimul_var.CnormFactor = size*size*optsyst[0].flux[0]*arith_image_total("scoringmask")*focscale*focscale; // /optsyst[0].nblambda;
            piaacmcsimul_var.CnormFactor = focscale*focscale*size*size*optsyst[0].flux[0]/optsyst[0].nblambda;

            if(piaacmcsimul_var.WRITE_OK==1)
            {
                sprintf(fname, "%s/CnormFactor.txt", piaacmcsimul_var.piaacmcconfdir);
                fp = fopen(fname, "w");
                fprintf(fp, "# Written by function PIAACMCsimul_computePSF()\n");
                fprintf(fp, "# Cnormfactor = focscale*focscale*size*size*optsyst[0].flux[0]/optsyst[0].nblambda\n");
                fprintf(fp, "# 0    focscale  size  flux[0]  nblambda\n");
                fprintf(fp,"\n");

                fprintf(fp, "%g\n", piaacmcsimul_var.CnormFactor);
                fprintf(fp, "0      %g %ld %g %d\n", focscale, size, optsyst[0].flux[0], optsyst[0].nblambda);
                fclose(fp);
            }

            printf("COMPUTING RESOLVED SOURCE PSF\n");
            printf("SCORINGTOTAL = %f  %f\n", piaacmcsimul_var.SCORINGTOTAL, arith_image_total("scoringmask"));
            printf("Peak constrast (rough estimate)= %g\n", peakcontrast/size/size/optsyst[0].flux[0]/focscale/focscale/normcoeff/normcoeff);
            avContrast = value/(arith_image_total("scoringmask")*focscale*focscale);
            printf("Total light in scoring field = %g, peak PSF = %g   -> Average contrast = %g\n", value, piaacmc[0].peakPSF, avContrast);
            

            if(outsave==1)
            {
				PIAACMCsimul_update_fnamedescr();				                                          
                sprintf(fname, "%s/contrast_exsrc%03d_sm%d.%s.txt", piaacmcsimul_var.piaacmcconfdir, sourcesize, piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.fnamedescr);                    
                fp = fopen(fname, "w");
                fprintf(fp, "%g", avContrast);
                fclose(fp);
            }

        }
        else // called for step 0 through 15.  Does not use OPDerr
        { // compute the PSF for a single point source at offset xld, yld
            printf("COMPUTING UNRESOLVED SOURCE PSF [%f x %f]\n", xld, yld);


            // ========== initializes optical system to piaacmc design ===========
            // xld and yld are the input positions of the input source in lamba/D
            // initialize first point source, which sets optsyst
            
            /// calls PIAACMCsimul_init() 
            PIAACMCsimul_init(piaacmc, 0, xld, yld);
            
            /// calls PIAACMCsimul_makePIAAshapes()
            PIAACMCsimul_makePIAAshapes(piaacmc, 0);



            // ============ perform propagations ================
            // propagate it (optsyst is a global), output in psfc0 (complex amlitude)
            // and psfi0 (intensity)
            OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcsimul_var.piaacmcconfdir, 0);

            if(outsave==1)
            {
				PIAACMCsimul_update_fnamedescr();
				sprintf(fname, "!%s/psfi0_ptsrc_sm%d.%s.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.fnamedescr);    
				save_fits("psfi0", fname);
            }


 //           list_image_ID();
            // linearize the result into imvect
            linopt_imtools_Image_to_vec("psfc0", "pixindex", "pixmult", "imvect");
           
           
            save_fits("imvect", "!./testdir/test_imvect.fits");
            
            
            // extract amplitude and phase for diagnostics
            //mk_amph_from_complex("psfc0", "psfc0a", "psfc0p", 0);
            //save_fits("psfc0a", "!test_psfc0a.fits");
            //list_image_ID();            
            //delete_image_ID("psfc0a");
            //delete_image_ID("psfc0p");
            //printf("saved -> test_psfc0a.fits\n");
            //fflush(stdout);
            
           // if(savepsf==2)
			//	{	
				//	mk_reim_from_complex("psfc0", "psfc0re", "psfc0im", 0);		
/*					save_fits("psfi0", "!test_psfi0.fits");
					save_fits("psfc0re", "!test_psfc0re.fits");
					save_fits("psfc0im", "!test_psfc0im.fits");
					sleep(100000);*/

			//	}
            
            // compute average contrast
            value = 0.0;
            peakcontrast = 0.0;
            ID = image_ID("imvect");
            for(ii=0; ii<data.image[ID].md[0].nelement; ii+=2)
            {
                // intensity as Re^2 + Im^2
                tmpv = data.image[ID].array.F[ii]*data.image[ID].array.F[ii] + data.image[ID].array.F[ii+1]*data.image[ID].array.F[ii+1];
                value += tmpv;
                if(tmpv>peakcontrast)
                    peakcontrast = tmpv;
            }
            // report the contrast
            for(elem=0; elem<optsyst[0].NBelem; elem++)
                printf("    FLUX %3ld   %12.4lf %8.6lf\n", elem, optsyst[0].flux[elem], optsyst[0].flux[elem]/optsyst[0].flux[0]);
//            value = value/size/size/optsyst[0].flux[0];

            if(outsave==1)
            {
				PIAACMCsimul_update_fnamedescr();
                sprintf(fname, "%s/flux_ptsrc_sm%d.%s.txt", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.fnamedescr); 
                fpflux = fopen(fname, "w");
                for(elem=0; elem<optsyst[0].NBelem; elem++)
                    fprintf(fpflux, "%18.16lf %18.16lf  %d\n", optsyst[0].flux[elem], optsyst[0].flux[elem]/optsyst[0].flux[0], optsyst[0].nblambda);
                fprintf(fpflux, "W1\n");
                fclose(fpflux);
        
                sprintf(command, "cp %s %s/tmp_flux.txt", fname, piaacmcsimul_var.piaacmcconfdir);
                ret = system(command);

           }

  
                
            //         piaacmcsimul_var.CnormFactor = size*size*optsyst[0].flux[0]*arith_image_total("scoringmask")*focscale*focscale; // /optsyst[0].nblambda;
            piaacmcsimul_var.CnormFactor = focscale*focscale*size*size*optsyst[0].flux[0]/optsyst[0].nblambda;
            sprintf(fname, "%s/CnormFactor.txt", piaacmcsimul_var.piaacmcconfdir);
            fp = fopen(fname, "w");
            fprintf(fp, "# Written by function PIAACMCsimul_computePSF()\n");
            fprintf(fp, "# Cnormfactor = focscale*focscale*size*size*optsyst[0].flux[0]/optsyst[0].nblambda\n");
            fprintf(fp, "# 1    focscale  size  flux[0]  nblambda\n");
            fprintf(fp,"\n");
            fprintf(fp, "%g\n", piaacmcsimul_var.CnormFactor);
            fprintf(fp, "1      %g %ld %g %d\n", focscale, size, optsyst[0].flux[0], optsyst[0].nblambda);
            fclose(fp);
            
            // here we're essentially done!


            printf("COMPUTING UNRESOLVED SOURCE PSF -*- [%f x %f]\n", xld, yld);
            printf("SCORINGTOTAL = %f  %f\n", piaacmcsimul_var.SCORINGTOTAL, arith_image_total("scoringmask"));
   
            if((IDv=variable_ID("PIAACMC_NOFPM"))!=-1)
                {
                    ID = image_ID("psfc0");
                    piaacmc[0].peakPSF = 0.0;
                    for(ii=0;ii<size*size;ii++)
                        {
                            val = data.image[ID].array.CF[ii].re*data.image[ID].array.CF[ii].re + data.image[ID].array.CF[ii].im*data.image[ID].array.CF[ii].im;
                            if(val>piaacmc[0].peakPSF)
                                piaacmc[0].peakPSF = val;
                        }

                    fp = fopen("conf/conf_peakPSF.txt", "w");
                    fprintf(fp, "%g\n", piaacmc[0].peakPSF);
                    fclose(fp);
                }
            if(piaacmc[0].peakPSF<1.0)
                {
                    printf("Peak constrast (rough estimate)= %g -> %g\n", peakcontrast, peakcontrast/(optsyst[0].flux[0]*optsyst[0].flux[0]));
                    printf("optsyst[0].flux[0]  = %g\n", optsyst[0].flux[0]);
                    printf("SCORINGMASKTYPE = %d\n", piaacmcsimul_var.SCORINGMASKTYPE);
                    avContrast = value/(optsyst[0].flux[0]*optsyst[0].flux[0])/piaacmcsimul_var.SCORINGTOTAL;
                    printf("[0] Total light in scoring field = %g, peak PSF = %g, SCOTINGTOTAL = %g   -> Average contrast = %g\n", value, piaacmc[0].peakPSF,  piaacmcsimul_var.SCORINGTOTAL, avContrast);
                }
            else
                {
                    printf("Peak constrast = %g -> %g\n", peakcontrast, peakcontrast/piaacmc[0].peakPSF);
                    printf("SCORINGMASKTYPE = %d\n", piaacmcsimul_var.SCORINGMASKTYPE);
                    avContrast = value/piaacmc[0].peakPSF/piaacmcsimul_var.SCORINGTOTAL;
                    printf("[1] Total light in scoring field = %g, peak PSF = %g, SCOTINGTOTAL = %g  -> Average contrast = %g\n", value, piaacmc[0].peakPSF, piaacmcsimul_var.SCORINGTOTAL, avContrast);
                    
                    if((fp=fopen("PSFcontrastval.txt","w"))!=NULL)
                    {
                        fprintf(fp, "%g %g\n", peakcontrast/piaacmc[0].peakPSF, value/piaacmc[0].peakPSF/piaacmcsimul_var.SCORINGTOTAL);
                        fclose(fp);
                    }
                }

            if(outsave==1)
            {    
				PIAACMCsimul_update_fnamedescr();
                sprintf(fname, "%s/contrast_ptsrc_sm%d.%s.txt", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.fnamedescr); 
                printf("saving contrast value [%g] -> %s\n", avContrast, fname);
                
                fp = fopen(fname, "w");
                fprintf(fp, "%g", avContrast);
                fclose(fp);
            }
        }
    }

    return(avContrast);
}

