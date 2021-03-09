/**
 * @file    PIAAACMCsimul_init.c
 * @brief   Initializes the optsyst structure to simulate PIAACMC system
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
#include <assert.h>


// milk includes

#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "image_gen/image_gen.h"


#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"


extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTSYST *optsyst;

extern OPTPIAACMCDESIGN *piaacmc;



 





/**
 * 
 * @brief Initializes the optsyst structure to simulate PIAACMC system
 * 
 * Fills in an OPTSYST global optsyst (see OptSysProp.h) which describes
 * the optical system as a series of planes based on the input design structure
 * 
 * TTxld and TTyld are tip/tilt x-y coordinates specifying the location of the source relative
 * to the optical axis in units of lambda/D
 * 
 * @note Index allows multiple configurations, but it's always 0.  Nonzero values are untested
 * 
 */
void PIAACMCsimul_init(
		OPTPIAACMCDESIGN *design, 
		long index, 
		double TTxld, 
		double TTyld 
		)
{
    FILE *fp;
    FILE *fpri;
    long k, i;
    long size;
    double x, y, PA;
    long ii, jj;
    long nblambda;
    long size2;
    double beamradpix;
    long kx, ky, kxy;
    long IDpiaaz0, IDpiaaz1;
    long surf;
    long IDa;
    char fname_pupa0[500];
    long ID;
    long IDr;
    long elem;
    char fname[500];
    long IDv;
    long IDopderr;
    long iDM; // DM index

    int savefpm;

	long ID_DFTmask00;
	double r;

    //    double ri, ri0, sag2opd_coeff;
    //    long IDpiaar0zsag, IDpiaar1zsag;
    //    int mkpiaar0zsag, mkpiaar1zsag;
    //    double sag2opd_coeff0;
    int IDpiaam0z, IDpiaam1z;

    int ret;
    char command[1000];


	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif


    assert(index == 0);   // test that index is always 0

    optsyst[0].nblambda = design[index].nblambda;
    nblambda = optsyst[0].nblambda;

    if(piaacmcsimul_var.PIAACMC_save==1)
    {
        sprintf(fname, "%s/lambdalist.txt", piaacmcsimul_var.piaacmcconfdir);
        fp = fopen(fname, "w");
    }


    printf("lambda = %g\n", design[index].lambda);
    printf("LAMBDASTART = %g\n", piaacmcsimul_var.LAMBDASTART);
    printf("LAMBDAEND = %g\n", piaacmcsimul_var.LAMBDAEND);

    // sets up the wavelengths over specifed bandwidth
    for(k=0; k<optsyst[0].nblambda; k++)
    {
        optsyst[0].lambdaarray[k] = piaacmcsimul_var.LAMBDASTART + (0.5+k)*( piaacmcsimul_var.LAMBDAEND - piaacmcsimul_var.LAMBDASTART)/optsyst[0].nblambda;
        if(piaacmcsimul_var.PIAACMC_save==1)
            fprintf(fp, "%02ld %20g\n", k, optsyst[0].lambdaarray[k]);
    }
    if(piaacmcsimul_var.PIAACMC_save==1)
        fclose(fp);

    // the physical pupil size in meters
    optsyst[0].beamrad = design[index].beamrad;
    // the number of pixels in each side of the square complex amplitude arrays at each plane
    optsyst[0].size = design[index].size;
    size = optsyst[0].size;
    size2 = size*size; // area
    optsyst[0].pixscale = design[index].pixscale;
    optsyst[0].DFTgridpad = 0; // 0 for full DFT sampling, >0 for faster execution
    // beam radius in pixels
    beamradpix = optsyst[0].beamrad/optsyst[0].pixscale;

	// Create "_DFTmask00" : force pixels within 10% of nominal pupil radius to be part of the DFT
/*	ID_DFTmask00 = create_2Dimage_ID("_DFTmask00", size, size);
	for(ii=0;ii<size;ii++)
		for(jj=0;jj<size;jj++)
			{
				x = (1.0*ii-0.5*size)/beamradpix;
				y = (1.0*jj-0.5*size)/beamradpix;
				r = sqrt(x*x+y*y);
				if(r<1.1)
					data.image[ID_DFTmask00].array.F[jj*size+ii] = 1.0;
				else
					data.image[ID_DFTmask00].array.F[jj*size+ii] = 0.0;					
			}
*/



    // printf("BEAM RADIUS = %f / %f  = %f pix,   piaacmc[0].beamrad = %f\n", optsyst[0].beamrad, optsyst[0].pixscale, beamradpix, piaacmc[0].beamrad );
    // sleep(10);

    // parameter that determines sampling of DFTs for progation onto the FPM
    // 0 => full sampling
    // 1 => every other pixel in each dimension
    // 2 => every third pixel in each dimension
    // etc: n => every (n+1)th pixel in each dimension
    // allows subsampling to speed up the DFT computation
    if((IDv=variable_ID("PIAACMC_dftgrid"))!=-1)
        optsyst[0].DFTgridpad = (long) (data.variable[IDv].value.f+0.001);


    // define optical elements and locations
    // have at least two aspheric mirrors in addition to the DMs
    // tyically design[index].nbDM = 0
    optsyst[0].NB_asphsurfm = 2+design[index].nbDM;
    // no aspheric lenses
    optsyst[0].NB_asphsurfr = 0;

    optsyst[0].NBelem = 100; // to be updated later

    if(piaacmcsimul_var.PIAACMC_save==1)
    {
        sprintf(fname, "%s/conjugations.txt", piaacmcsimul_var.piaacmcconfdir);
        fp = fopen(fname, "w");
    }





    elem = 0;
    // ------------------- elem 0: input pupil -----------------------
    sprintf(optsyst[0].name[elem], "input pupil");
    optsyst[0].elemtype[elem] = 1; // pupil mask
    // input pupil from file - will always exist
    sprintf(fname_pupa0, "%s/pupa0_%ld.fits", piaacmcsimul_var.piaacmcconfdir, size);

    if(file_exists(fname_pupa0)==1)
        load_fits(fname_pupa0, "pupa0", 1);

    IDa = image_ID("pupa0");
    if(IDa==-1) // if pupil does not exist, use circular one (this occurs in initial design steps)
    { 
        printf("CREATING INPUT PUPIL\n");
        if(IDa!=-1)
            delete_image_ID("pupa0");
        IDa = create_3Dimage_ID("pupa0", size, size, nblambda);

		IDr = image_ID("rcoord");

        ID = image_ID("telpup");
        if(ID==-1)
            if(file_exists("telpup.fits")==1)
                ID = load_fits("telpup.fits", "telpup", 1);


        if(ID==-1)
        {
            for(k=0; k<nblambda; k++)
                for(ii=0; ii<size2; ii++)
                {
                    if((data.image[IDr].array.F[ii]>design[index].centObs0)&&(data.image[IDr].array.F[ii]<1.0))
                        data.image[IDa].array.F[k*size2+ii] = 1.0;
                    else
                        data.image[IDa].array.F[k*size2+ii] = 0.0;
                }
        }
        else
            for(k=0; k<nblambda; k++)
                for(ii=0; ii<size2; ii++)
                {
                    if(data.image[ID].array.F[ii]>0.5)
                        data.image[IDa].array.F[k*size2+ii] = 1.0;
                    else
                        data.image[IDa].array.F[k*size2+ii] = 0.0;
                }

        sprintf(fname_pupa0, "!%s/pupa0_%ld.fits", piaacmcsimul_var.piaacmcconfdir, size);
        save_fl_fits("pupa0", fname_pupa0);
    }
    optsyst[0].elemarrayindex[elem] = IDa;
    optsyst[0].elemZpos[elem] = 0.0; // pupil is at z = 0
    if(piaacmcsimul_var.PIAACMC_save==1)
        fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
    elem++;







    // set up the shape of a mirror to insert tip/tilt and optical error
    
    // initialize this mirror by setting pointing (simulated as mirror shape), defining off-axis source
    ID = create_2Dimage_ID("TTm", size, size);

    for(ii=0; ii<size; ii++)
        for(jj=0; jj<size; jj++)
        {
            x = (1.0*ii-0.5*size)/beamradpix; // position of pixel (ii,jj) in pupil radius units
            y = (1.0*jj-0.5*size)/beamradpix;
            // set the mirror shape as a linear tilt of pixel position reflecting the tilt
            data.image[ID].array.F[jj*size+ii] = 0.25*(TTxld*x+TTyld*y)*(piaacmcsimul_var.LAMBDAEND + piaacmcsimul_var.LAMBDASTART)*0.5; // xld -> half-OPD
        }


    // add OPD error on TTM if it exists
    IDopderr = image_ID("opderr");
    if(IDopderr != -1)
    {
        for(ii=0; ii<size; ii++)
            for(jj=0; jj<size; jj++)
            {
                x = (1.0*ii-0.5*size)/beamradpix;
                y = (1.0*jj-0.5*size)/beamradpix;
                // add the error shape to the mirror shape
                data.image[ID].array.F[jj*size+ii] += data.image[IDopderr].array.F[jj*size+ii]*0.5;
            }
    }

    // sprintf(fname, "!%s/TTm.fits", piaacmcsimul_var.piaacmcconfdir);
    // save_fits("TTm", fname);

    // finish the definition of the TT mirror specifying various properties
    sprintf(optsyst[0].name[elem], "TT mirror");
    optsyst[0].elemtype[elem] = 3; // reflective mirror
    optsyst[0].elemarrayindex[elem] = 0; // not used because this field is only relevant for
                                        // DM or aspheric mirrors, which this mirror is not
                                        // this mirror is "flat" except for possible injected OPD error
    optsyst[0].ASPHSURFMarray[0].surfID = ID; // store array ID
    optsyst[0].elemZpos[elem] = 0.0; // put it at the entrance pupil
    if(piaacmcsimul_var.PIAACMC_save==1)
        fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
    //        fprintf(fp,"%02ld  %f    Fold mirror used to induce pointing offsets\n", elem, optsyst[0].elemZpos[elem]);
    elem++;




    // set up the deformable mirrors
    // tyically design[index].nbDM = 0 so we will skip this
    for(iDM=0; iDM<design[index].nbDM; iDM++)
    {
        // ----------------- DM (s) -----------------------------------------------
        sprintf(optsyst[0].name[elem], "DM %ld", iDM);
        optsyst[0].elemtype[elem] = 3; // reflective element
        optsyst[0].elemarrayindex[elem] = 3+iDM; // index
        optsyst[0].ASPHSURFMarray[optsyst[0].elemarrayindex[elem]].surfID = design[index].ID_DM[iDM];
        optsyst[0].elemZpos[elem] = design[index].DMpos[iDM];
        if(piaacmcsimul_var.PIAACMC_save==1)
            fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
        //          fprintf(fp,"%02ld  %f    DM %ld\n", elem, optsyst[0].elemZpos[elem], iDM);
        elem++;
    }




    // shape/sag for the first aspheric mirror
    IDpiaam0z = image_ID("piaam0z");  // nominal sag (mirror equivalent)
    // shape/sag for the second aspheric mirror
    IDpiaam1z = image_ID("piaam1z");  //







    // ------------------- [OPTIONAL] pre-apodizer  -----------------------
    // typically not present for PIAACMC
    ID = image_ID("prePIAA0mask");
    if(ID==-1)
        ID = load_fits("prePIAA0mask.fits", "prePIAA0mask", 1);
    if(ID!=-1)
    {
        // tell the design that this element exists (was found on disk)
        design[index].prePIAA0mask = 1;
        sprintf(optsyst[0].name[elem], "pupil plane apodizer");
        optsyst[0].elemtype[elem] = 1; // opaque mask
        optsyst[0].elemarrayindex[elem] = ID;
        optsyst[0].elemZpos[elem] = design[index].prePIAA0maskpos;

        if(piaacmcsimul_var.PIAACMC_save==1)
            fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
        elem++;
    }
    else
        design[index].prePIAA0mask = 0;

    // PIAAmode = 1 => this is a PIAA system
    // PIAAmode = 0 => this is not a PIAA system
    if(piaacmc[0].PIAAmode == 1)
    {
        // ------------------- elem 2:  PIAA M/L 0  -----------------------
        // (M/L is "mirror or lens" - in our case it's a mirror)
        sprintf(optsyst[0].name[elem], "PIAA optics 0");
        optsyst[0].elemtype[elem] = 3; // reflective PIAA M0
        // this is an aspheric mirror, so we need actual sag shapes
        optsyst[0].elemarrayindex[elem] = 1; // index = 1 implied aspheric
        optsyst[0].elemZpos[elem] = design[index].PIAA0pos;  // location of this element relative to pupil

        printf("============ (2) PIAA0pos = %f ==================\n", optsyst[0].elemZpos[elem]);
//        sleep(5);
        // set up the element properties
        if(design[index].PIAAmaterial_code == 0) // mirror
            // specify the sag array, put in data global by routine named something like createPIAA_mirror_shapes
            optsyst[0].ASPHSURFMarray[optsyst[0].elemarrayindex[elem]].surfID = IDpiaam0z;
        else // lens
        {
            optsyst[0].elemtype[elem] = 4;
            optsyst[0].ASPHSURFRarray[optsyst[0].elemarrayindex[elem]].surfID = image_ID("piaar0zsag"); //IDpiaar0zsag;
            optsyst[0].ASPHSURFRarray[optsyst[0].elemarrayindex[elem]].mat0 = 100;
            optsyst[0].ASPHSURFRarray[optsyst[0].elemarrayindex[elem]].mat1 = design[0].PIAAmaterial_code; // vacuum
        }
        // make sure the above did something
        if(optsyst[0].ASPHSURFMarray[1].surfID==-1)
        {
            printf("ERROR: surface 0 not identified\n");
            list_image_ID();
            exit(0);
        }
        // print this element to tracking file if desired
        if(piaacmcsimul_var.PIAACMC_save==1)
            fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
        elem++;
    }





    // ------------------- [OPTIONAL] opaque mask after last elem -----------------------
    // get opaque mask from the file, with a standard filename for the first mask
    // we don't have one in the nominal design
    ID = load_fits("postPIAA0mask.fits", "postPIAA0mask", 1);
    if(ID!=-1)
    {
        // tell the design that this element exists (was found on disk)
        design[index].postPIAA0mask = 1;
        sprintf(optsyst[0].name[elem], "opaque mask after PIAA element 0");
        optsyst[0].elemtype[elem] = 1; // opaque mask
        optsyst[0].elemarrayindex[elem] = ID;
        optsyst[0].elemZpos[elem] = design[index].postPIAA0maskpos; // get position from design input

        if(piaacmcsimul_var.PIAACMC_save==1)
            fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
        elem++;
    }
    else
        design[index].postPIAA0mask = 0;



    // if we're a PIAA
    if(piaacmc[0].PIAAmode == 1)
    {
        // add one more mirror and mask
        // ------------------- elem 3: reflective PIAA M1  -----------------------
        sprintf(optsyst[0].name[elem], "PIAA optics 1");
        optsyst[0].elemtype[elem] = 3; // reflective PIAA M1
        optsyst[0].elemarrayindex[elem] = 2;
        optsyst[0].elemZpos[elem] = design[index].PIAA0pos + design[index].PIAAsep;

        if(design[index].PIAAmaterial_code == 0) // mirror
            optsyst[0].ASPHSURFMarray[2].surfID = IDpiaam1z;
        else // lens
        {
            optsyst[0].elemtype[elem] = 4;
            optsyst[0].ASPHSURFRarray[2].surfID = image_ID("piaar1zsag"); //IDpiaar0zsag;
            optsyst[0].ASPHSURFRarray[optsyst[0].elemarrayindex[elem]].mat0 = 100; // vacuum
            optsyst[0].ASPHSURFRarray[optsyst[0].elemarrayindex[elem]].mat1 = design[0].PIAAmaterial_code;
        }


        if(optsyst[0].ASPHSURFMarray[2].surfID==-1)
        {
            printf("ERROR: surface 1 not identified\n");
            list_image_ID();
            exit(0);
        }

        if(piaacmcsimul_var.PIAACMC_save==1)
            fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
        //       fprintf(fp,"%02ld  %f    PIAAM1\n", elem, optsyst[0].elemZpos[elem]);
        elem++;





        // ------------------- elem 4 opaque mask at reflective PIAA M1  -----------------------
        sprintf(optsyst[0].name[elem], "opaque mask at PIAA elem 1");
        optsyst[0].elemtype[elem] = 1; // opaque mask
        ID = load_fits("piaa1mask.fits", "piaa1mask", 1);
        if(ID==-1)
            ID = make_disk("piaa1mask", size, size, 0.5*size, 0.5*size, design[index].r1lim*beamradpix);
        optsyst[0].elemarrayindex[elem] = ID;
        optsyst[0].elemZpos[elem] = optsyst[0].elemZpos[elem-1];

        if(piaacmcsimul_var.PIAACMC_save==1)
            fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
        //        fprintf(fp,"%02ld  %f    PIAAM1 edge opaque mask\n", elem, optsyst[0].elemZpos[elem]);
        elem++;
    }



    // --------------------  elem 5: focal plane mask ------------------------
    if((IDv=variable_ID("PIAACMC_NOFPM"))==-1)
    {
        sprintf(optsyst[0].name[elem], "post focal plane mask pupil");
        optsyst[0].elemtype[elem] = 5; // focal plane mask
        optsyst[0].elemarrayindex[elem] = 0;

        printf("=========== MAKE FOCAL PLANE MASK ===========\n");
        savefpm = 0;
        if((IDv=variable_ID("PIAACMC_SAVE_fpm"))!=-1)
            savefpm = (int) (data.variable[IDv].value.f+0.001);

        /// call PIAACMCsimul_mkFocalPlaneMask() to make the focal plane mask
        optsyst[0].FOCMASKarray[0].fpmID = PIAACMCsimul_mkFocalPlaneMask("fpmzmap", "piaacmcfpm", piaacmcsimul_var.focmMode, savefpm); // if -1, this is 1-fpm; otherwise, this is impulse response from single zone

		// TEST
		
		mk_reim_from_complex("piaacmcfpm", "piaacmcfpm_re", "piaacmcfpm_im", 0);
        save_fits("piaacmcfpm_re", "!./testdir/test_piaacmcfpm_re.fits");
        save_fits("piaacmcfpm_im", "!./testdir/test_piaacmcfpm_im.fits");
        delete_image_ID("piaacmcfpm_re");
        delete_image_ID("piaacmcfpm_im");
		


        // zfactor is the zoom factor for the DFT, driving sample resolution in the z direction
        // to allow faster DFTs.  Similar to DFTgridpad.
        optsyst[0].FOCMASKarray[0].zfactor = design[index].fpzfactor;
        // set the position of the pupil from which the DFT propagates to the FPM
        // NOT the position along the beam of the FPM.  That's OK and intended.
        // For this element, this defines the conjugation of the pupil from which we are computing the DFT
        optsyst[0].elemZpos[elem] = optsyst[0].elemZpos[elem-1]; // plane from which FT is done
        if(piaacmcsimul_var.PIAACMC_save==1)
            fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
        //      fprintf(fp,"%02ld  %f    post-focal plane mask pupil\n", elem, optsyst[0].elemZpos[elem]);
		printf("=========== FOCAL PLANE MASK : DONE ===========\n");
        elem++;
    }



    // if we're a PIAA
    if(piaacmc[0].PIAAmode == 1)
    {
        // if the the inverse PIAA is prior to the Lyot stop
        // (invPIAAmode = 0 has no inverse PIAA, = 1 has Lyot stop prior to inverse PIAA)
        // there is no inverse PIAA in the WFIRST design
        if(design[index].invPIAAmode == 2) // inv PIAA -> Lyot stops
        {
            // --------------------  elem 8: inv PIAA1 ------------------------
            sprintf(optsyst[0].name[elem], "invPIAA optics 1");

            if(design[index].PIAAmaterial_code == 0) // mirror
                optsyst[0].elemtype[elem] = 3; // reflective PIAA M/L 1
            else
                optsyst[0].elemtype[elem] = 4; // refractive PIAA M/L 1

            optsyst[0].elemarrayindex[elem] = 2;
            // put an element at z=0 in conjugation space (conjugate to the pupil)
            optsyst[0].elemZpos[elem] = 0.0;
            if(piaacmcsimul_var.PIAACMC_save==1)
                fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
            //          fprintf(fp,"%02ld  %f    invPIAA1\n", elem, optsyst[0].elemZpos[elem]);
            elem++;

            // --------------------  elem 9: inv PIAA0 ------------------------
            sprintf(optsyst[0].name[elem], "invPIAA optics 0");

            if(design[index].PIAAmaterial_code == 0) //  mirror
                optsyst[0].elemtype[elem] = 3; // reflective PIAA M/L 0
            else
                optsyst[0].elemtype[elem] = 4; // refractive PIAA M/L 0

            optsyst[0].elemarrayindex[elem] = 1;
            // previous element + PIAAsep
            optsyst[0].elemZpos[elem] = design[index].PIAAsep;
            if(piaacmcsimul_var.PIAACMC_save==1)
                fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
            //         fprintf(fp,"%02ld  %f    invPIAA0\n", elem, optsyst[0].elemZpos[elem]);
            elem++;
        }
    }


    // --------------------  Lyot masks  ------------------------
    // add Lyot masks as specified in the design
    for(i=0; i<design[index].NBLyotStop; i++)
    {
        sprintf(optsyst[0].name[elem], "Lyot mask %ld", i);
        optsyst[0].elemtype[elem] = 1; // Lyot mask
        optsyst[0].elemarrayindex[elem] = design[index].IDLyotStop[i];
        printf("elem %ld  Lyot mask %ld : %ld\n", elem, i, design[index].IDLyotStop[i]);
        optsyst[0].elemZpos[elem] =  design[index].LyotStop_zpos[i];
        if(piaacmcsimul_var.PIAACMC_save==1)
            fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
        //          fprintf(fp,"%02ld  %f  Lyot Stop %ld\n", elem, optsyst[0].elemZpos[elem], i);
        elem++;
    }

    // if we're a PIAA
    if(piaacmc[0].PIAAmode == 1)
    {
        // add stops for the inverse PIAA
        // not in WFIRST design, skipping...
        if(design[index].invPIAAmode == 1) // Lyot masks -> inv PIAA
        {
            // --------------------  elem 8: inv PIAA1 ------------------------
            sprintf(optsyst[0].name[elem], "invPIAA optics 1");
            if(design[index].PIAAmaterial_code == 0) // mirror
                optsyst[0].elemtype[elem] = 3; // reflective PIAA M/L 1
            else
                optsyst[0].elemtype[elem] = 4; // refractive PIAA M/L 1
            optsyst[0].elemarrayindex[elem] = 2;
            optsyst[0].elemZpos[elem] = 0.0;
            if(piaacmcsimul_var.PIAACMC_save==1)
                fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
            //           fprintf(fp,"%02ld  %f    invPIAA1\n", elem, optsyst[0].elemZpos[elem]);
            elem++;

            // --------------------  elem 9: inv PIAA0 ------------------------
            sprintf(optsyst[0].name[elem], "invPIAA optics 0");
            if(design[index].PIAAmaterial_code == 0) //  mirror
                optsyst[0].elemtype[elem] = 3; // reflective PIAA M/L 0
            else
                optsyst[0].elemtype[elem] = 4; // refractive PIAA M/L 0
            optsyst[0].elemarrayindex[elem] = 1;
            optsyst[0].elemZpos[elem] = design[index].PIAAsep;
            if(piaacmcsimul_var.PIAACMC_save==1)
                fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
            //           fprintf(fp,"%02ld  %f    invPIAA0\n", elem, optsyst[0].elemZpos[elem]);
            elem++;
        }
    }


    // if we're a PIAA
    if(piaacmc[0].PIAAmode == 1)
    {
        // --------------------  elem 9: back end mask  ------------------------
        // not in WFIRST design, skipping, but it looks very straightforward

        sprintf(optsyst[0].name[elem], "back end pupil stop  (rad = %f)", design[index].pupoutmaskrad);

        optsyst[0].elemtype[elem] = 1;
        ID = make_disk("pupoutmask", size, size, 0.5*size, 0.5*size, design[index].pupoutmaskrad*design[index].beamrad/design[index].pixscale);
        optsyst[0].elemarrayindex[elem] = ID;
        optsyst[0].elemZpos[elem] =  optsyst[0].elemZpos[elem-1];
        if(piaacmcsimul_var.PIAACMC_save==1)
            fprintf(fp,"%02ld  %f    %s\n", elem, optsyst[0].elemZpos[elem], optsyst[0].name[elem]);
        //     fprintf(fp,"%02ld  %f   back end mask\n", elem, optsyst[0].elemZpos[elem]);
        elem++;
    }



    if(piaacmcsimul_var.PIAACMC_save==1)
        fclose(fp);

    optsyst[0].NBelem = elem;
    optsyst[0].endmode = 0;


	


    piaacmcsimul_var.optsystinit = 1;
}

