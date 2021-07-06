/**
 * @file    PIAAACMCsimul_init.c
 * @brief   Initializes the optsyst structure to simulate PIAACMC system
 *
 * Can design both APLCMC and PIAACMC coronagraphs
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

#include "PIAACMCsimul.h"
#include "FocalPlaneMask/mkFocalPlaneMask.h"



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
 *
 */
errno_t init_piaacmcopticalsystem(
    double TTxld,
    double TTyld
)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %f %f", TTxld, TTyld);

    long size;
    long nblambda;
    uint64_t size2;
    double beamradpix;
    long elem;

    int savefpm;


    piaacmcopticalsystem.nblambda = piaacmcopticaldesign.nblambda;
    nblambda = piaacmcopticalsystem.nblambda;

    printf("lambda = %g\n", piaacmcopticaldesign.lambda);
    printf("LAMBDASTART = %g\n", piaacmcparams.LAMBDASTART);
    printf("LAMBDAEND = %g\n", piaacmcparams.LAMBDAEND);





    // sets up the wavelengths over specifed bandwidth
    for(int k = 0; k < piaacmcopticalsystem.nblambda; k++)
    {
        piaacmcopticalsystem.lambdaarray[k] =
            piaacmcparams.LAMBDASTART + (0.5 + k) *
            (piaacmcparams.LAMBDAEND - piaacmcparams.LAMBDASTART) /
            piaacmcopticalsystem.nblambda;
    }

    if(piaacmcparams.PIAACMC_save == 1)
    {   // write wavelength list
        char fname[STRINGMAXLEN_FULLFILENAME];
        FILE *fp;

        WRITE_FULLFILENAME(
            fname,
            "%s/lambdalist.txt",
            piaacmcparams.piaacmcconfdir
        );

        fp = fopen(fname, "w");
        if(fp == NULL)
        {
            PRINT_ERROR("Cannot create file %s", fname);
            abort();
        }

        for(int k = 0; k < piaacmcopticalsystem.nblambda; k++)
        {
            fprintf(fp, "%02d %20g\n", k, piaacmcopticalsystem.lambdaarray[k]);
        }

        fclose(fp);
    }



    // the physical pupil size in meters
    piaacmcopticalsystem.beamrad = piaacmcopticaldesign.beamrad;

    // the number of pixels in each side of the square complex amplitude arrays at each plane
    piaacmcopticalsystem.size = piaacmcopticaldesign.size;
    size = piaacmcopticalsystem.size;
    size2 = size * size; // area

    piaacmcopticalsystem.pixscale = piaacmcopticaldesign.pixscale;
    piaacmcopticalsystem.DFTgridpad = 0; // 0 for full DFT sampling, >0 for faster execution

    // beam radius in pixels
    beamradpix = piaacmcopticalsystem.beamrad / piaacmcopticalsystem.pixscale;

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



    // printf("BEAM RADIUS = %f / %f  = %f pix,   piaacmcopticaldesign.beamrad = %f\n", piaacmcopticalsystem.beamrad, piaacmcopticalsystem.pixscale, beamradpix, piaacmcopticaldesign.beamrad );
    // sleep(10);

    // parameter that determines sampling of DFTs for progation onto the FPM
    // 0 => full sampling
    // 1 => every other pixel in each dimension
    // 2 => every third pixel in each dimension
    // etc: n => every (n+1)th pixel in each dimension
    // allows subsampling to speed up the DFT computation
    {
        variableID IDv;
        if((IDv = variable_ID("PIAACMC_dftgrid")) != -1)
        {
            piaacmcopticalsystem.DFTgridpad = (long)(data.variable[IDv].value.f + 0.001);
        }
    }


    // define optical elements and locations
    // have at least two aspheric mirrors in addition to the DMs
    // tyically design[index].nbDM = 0
    piaacmcopticalsystem.NB_asphsurfm = 2 + piaacmcopticaldesign.nbDM;
    // no aspheric lenses
    piaacmcopticalsystem.NB_asphsurfr = 0;

    piaacmcopticalsystem.NBelem = 100; // to be updated later







    elem = 0;
    // ------------------- elem 0: input pupil -----------------------
    snprintf(piaacmcopticalsystem.name[elem], STRINGMAXLEN_OPTSYST_ELEMNAME, "input pupil");
    piaacmcopticalsystem.elemtype[elem] = 1; // pupil mask
    // input pupil from file - will always exist
    char fname_pupa0[STRINGMAXLEN_FULLFILENAME];
    WRITE_FULLFILENAME(
        fname_pupa0,
        "%s/pupa0_%ld.fits",
        piaacmcparams.piaacmcconfdir,
        size
    );

    if(file_exists(fname_pupa0) == 1)
    {
        load_fits(
            fname_pupa0,
            "pupa0",
            LOADFITS_ERRMODE_ERROR,
            NULL
        );
    }

    imageID IDa = image_ID("pupa0");
    if(IDa == -1) // if pupil does not exist, use circular one (this occurs in initial design steps)
    {
        printf("CREATING INPUT PUPIL\n");
        create_3Dimage_ID("pupa0", size, size, nblambda, &IDa);

        imageID IDr = image_ID("rcoord");

        imageID ID = image_ID("telpup");
        if(ID == -1)
            if(file_exists("telpup.fits") == 1)
            {
                load_fits("telpup.fits", "telpup", 1, &ID);
            }


        if(ID == -1)
        {
            for(long k = 0; k < nblambda; k++)
                for(uint64_t ii = 0; ii < size2; ii++)
                {
                    if((data.image[IDr].array.F[ii] > piaacmcopticaldesign.centObs0)
                            && (data.image[IDr].array.F[ii] < 1.0))
                    {
                        data.image[IDa].array.F[k * size2 + ii] = 1.0;
                    }
                    else
                    {
                        data.image[IDa].array.F[k * size2 + ii] = 0.0;
                    }
                }
        }
        else
            for(long k = 0; k < nblambda; k++)
                for(uint64_t ii = 0; ii < size2; ii++)
                {
                    if(data.image[ID].array.F[ii] > 0.5)
                    {
                        data.image[IDa].array.F[k * size2 + ii] = 1.0;
                    }
                    else
                    {
                        data.image[IDa].array.F[k * size2 + ii] = 0.0;
                    }
                }

        WRITE_FULLFILENAME(fname_pupa0,
                           "%s/pupa0_%ld.fits",
                           piaacmcparams.piaacmcconfdir,
                           size);

        FUNC_CHECK_RETURN(
            save_fl_fits("pupa0", fname_pupa0)
        );
    }
    piaacmcopticalsystem.elemarrayindex[elem] = IDa;
    piaacmcopticalsystem.elemZpos[elem] = 0.0; // pupil is at z = 0
    elem++;







    // set up the shape of a mirror to insert tip/tilt and optical error

    // initialize this mirror by setting pointing (simulated as mirror shape), defining off-axis source
    {
        imageID ID;
        FUNC_CHECK_RETURN(
            create_2Dimage_ID("TTm", size, size, &ID)
        );

        for(uint32_t ii = 0; ii < size; ii++)
            for(uint32_t jj = 0; jj < size; jj++)
            {

                // position of pixel (ii,jj) in pupil radius units
                double x = (1.0 * ii - 0.5 * size) / beamradpix;
                double y = (1.0 * jj - 0.5 * size) / beamradpix;

                // set the mirror shape as a linear tilt of pixel position reflecting the tilt
                data.image[ID].array.F[jj * size + ii] = 0.25 * (TTxld * x + TTyld * y) *
                        (piaacmcparams.LAMBDAEND + piaacmcparams.LAMBDASTART) *
                        0.5; // xld -> half-OPD
            }


        // add OPD error on TTM if it exists
        imageID IDopderr = image_ID("opderr");
        if(IDopderr != -1)
        {
            for(uint32_t ii = 0; ii < size; ii++)
                for(uint32_t jj = 0; jj < size; jj++)
                {
                    // add the error shape to the mirror shape
                    data.image[ID].array.F[jj * size + ii] +=
                        data.image[IDopderr].array.F[jj * size + ii] * 0.5;
                }
        }

        // sprintf(fname, "%s/TTm.fits", piaacmcparams.piaacmcconfdir);
        // save_fits("TTm", fname);

        // finish the definition of the TT mirror specifying various properties
        snprintf(piaacmcopticalsystem.name[elem], STRINGMAXLEN_OPTSYST_ELEMNAME, "TT mirror");
        piaacmcopticalsystem.elemtype[elem] = 3; // reflective mirror

        piaacmcopticalsystem.elemarrayindex[elem] = 0;
        // not used because this field is only relevant for
        // DM or aspheric mirrors, which this mirror is not
        // this mirror is "flat" except for possible injected OPD error

        // store array ID
        piaacmcopticalsystem.ASPHSURFMarray[0].surfID = ID;
    }

    // put it at the entrance pupil
    piaacmcopticalsystem.elemZpos[elem] = 0.0;
    //        fprintf(fp,"%02ld  %f    Fold mirror used to induce pointing offsets\n", elem, piaacmcopticalsystem.elemZpos[elem]);
    elem++;




    // set up the deformable mirrors
    // tyically design[index].nbDM = 0 so we will skip this
    for(int iDM = 0; iDM < piaacmcopticaldesign.nbDM; iDM++)
    {
        // ----------------- DM (s) -----------------------------------------------
        snprintf(piaacmcopticalsystem.name[elem], STRINGMAXLEN_OPTSYST_ELEMNAME, "DM %d", iDM);
        piaacmcopticalsystem.elemtype[elem] = 3; // reflective element
        piaacmcopticalsystem.elemarrayindex[elem] = 3 + iDM; // index

        piaacmcopticalsystem.ASPHSURFMarray[piaacmcopticalsystem.elemarrayindex[elem]].surfID =
            piaacmcopticaldesign.ID_DM[iDM];

        piaacmcopticalsystem.elemZpos[elem] = piaacmcopticaldesign.DMpos[iDM];
        //          fprintf(fp,"%02ld  %f    DM %ld\n", elem, piaacmcopticalsystem.elemZpos[elem], iDM);
        elem++;
    }










    // ------------------- [OPTIONAL] pre-apodizer  -----------------------
    // typically not present for PIAACMC
    {
        imageID ID = image_ID("prePIAA0mask");
        if(ID == -1)
        {
            FUNC_CHECK_RETURN(
                load_fits("prePIAA0mask.fits", "prePIAA0mask", LOADFITS_ERRMODE_WARNING, &ID)
            );
        }
        if(ID != -1)
        {
            // tell the design that this element exists (was found on disk)
            piaacmcopticaldesign.prePIAA0mask = 1;
            snprintf(piaacmcopticalsystem.name[elem], STRINGMAXLEN_OPTSYST_ELEMNAME, "pupil plane apodizer");
            piaacmcopticalsystem.elemtype[elem] = 1; // opaque mask
            piaacmcopticalsystem.elemarrayindex[elem] = ID;
            piaacmcopticalsystem.elemZpos[elem] = piaacmcopticaldesign.prePIAA0maskpos;
            elem++;
        }
        else
        {
            piaacmcopticaldesign.prePIAA0mask = 0;
        }
    }

    // PIAAmode = 1 => this is a PIAA system
    // PIAAmode = 0 => this is not a PIAA system
    if(piaacmcopticaldesign.PIAAmode == 1)
    {
        // shape/sag for the first aspheric mirror
        imageID IDpiaam0z = image_ID("piaam0z");  // nominal sag (mirror equivalent)


        // ------------------- elem 2:  PIAA M/L 0  -----------------------
        // (M/L is "mirror or lens" )
        snprintf(piaacmcopticalsystem.name[elem], STRINGMAXLEN_OPTSYST_ELEMNAME, "PIAA optics 0");
        piaacmcopticalsystem.elemtype[elem] = 3; // reflective PIAA M0
        // this is an aspheric mirror, so we need actual sag shapes
        piaacmcopticalsystem.elemarrayindex[elem] = 1; // index = 1 implied aspheric
        piaacmcopticalsystem.elemZpos[elem] =
            piaacmcopticaldesign.PIAA0pos;  // location of this element relative to pupil

        printf(
            "============ (2) PIAA0pos = %f ==================\n",
            piaacmcopticalsystem.elemZpos[elem]
        );

        // set up the element properties
        if(piaacmcopticaldesign.PIAAmaterial_code == 0) // mirror
            // specify the sag array, put in data global by routine named something like createPIAA_mirror_shapes
        {
            DEBUG_TRACEPOINT("Element %ld is mirror", elem);
            piaacmcopticalsystem.ASPHSURFMarray[piaacmcopticalsystem.elemarrayindex[elem]].surfID = IDpiaam0z;
            if(piaacmcopticalsystem.ASPHSURFMarray[piaacmcopticalsystem.elemarrayindex[elem]].surfID == -1)
            {
                FUNC_RETURN_FAILURE("PIAA mirror surface 0 not identified");
            }
        }
        else // lens
        {
            piaacmcopticalsystem.elemtype[elem] = 4;

            DEBUG_TRACEPOINT("Element %ld is lens", elem);
            piaacmcopticalsystem.ASPHSURFRarray[piaacmcopticalsystem.elemarrayindex[elem]].surfID = image_ID("piaar0zsag");

            if(piaacmcopticalsystem.ASPHSURFRarray[piaacmcopticalsystem.elemarrayindex[elem]].surfID == -1)
            {
                FUNC_RETURN_FAILURE("PIAA lens surface 0 not identified");
            }

            // material, vacuum
            piaacmcopticalsystem.ASPHSURFRarray[piaacmcopticalsystem.elemarrayindex[elem]].mat0 = 100;
            piaacmcopticalsystem.ASPHSURFRarray[piaacmcopticalsystem.elemarrayindex[elem]].mat1 = piaacmcopticaldesign.PIAAmaterial_code;
        }

        elem++;
    }





    // ------------------- [OPTIONAL] opaque mask after last elem -----------------------
    // get opaque mask from the file, with a standard filename for the first mask
    // we don't have one in the nominal design
    {
        imageID ID;
        load_fits("postPIAA0mask.fits", "postPIAA0mask", LOADFITS_ERRMODE_IGNORE, &ID);
        if(ID != -1)
        {
            // tell the design that this element exists (was found on disk)
            piaacmcopticaldesign.postPIAA0mask = 1;
            snprintf(piaacmcopticalsystem.name[elem], STRINGMAXLEN_OPTSYST_ELEMNAME, "opaque mask after PIAA element 0");
            piaacmcopticalsystem.elemtype[elem] = 1; // opaque mask
            piaacmcopticalsystem.elemarrayindex[elem] = ID;
            piaacmcopticalsystem.elemZpos[elem] =
                piaacmcopticaldesign.postPIAA0maskpos; // get position from design input
            elem++;
        }
        else
        {
            piaacmcopticaldesign.postPIAA0mask = 0;
        }
    }



    // if we're a PIAA
    if(piaacmcopticaldesign.PIAAmode == 1)
    {
        // shape/sag for the second aspheric mirror
        imageID IDpiaam1z = image_ID("piaam1z");

        // add one more mirror and mask
        // ------------------- elem 3: reflective PIAA M1  -----------------------
        snprintf(piaacmcopticalsystem.name[elem], STRINGMAXLEN_OPTSYST_ELEMNAME, "PIAA optics 1");
        piaacmcopticalsystem.elemtype[elem] = 3; // reflective PIAA M1
        piaacmcopticalsystem.elemarrayindex[elem] = 2;
        piaacmcopticalsystem.elemZpos[elem] = piaacmcopticaldesign.PIAA0pos + piaacmcopticaldesign.PIAAsep;

        if(piaacmcopticaldesign.PIAAmaterial_code == 0) // mirror
        {
            DEBUG_TRACEPOINT("Element %ld is mirror", elem);
            piaacmcopticalsystem.ASPHSURFMarray[piaacmcopticalsystem.elemarrayindex[elem]].surfID = IDpiaam1z;
            if(piaacmcopticalsystem.ASPHSURFMarray[piaacmcopticalsystem.elemarrayindex[elem]].surfID == -1)
            {
                FUNC_RETURN_FAILURE("PIAA mirror surface 1 not identified");
            }
        }
        else // lens
        {
            piaacmcopticalsystem.elemtype[elem] = 4;

            DEBUG_TRACEPOINT("Element %ld is lens", elem);
            piaacmcopticalsystem.ASPHSURFRarray[piaacmcopticalsystem.elemarrayindex[elem]].surfID = image_ID("piaar1zsag");

            if(piaacmcopticalsystem.ASPHSURFRarray[piaacmcopticalsystem.elemarrayindex[elem]].surfID == -1)
            {
                FUNC_RETURN_FAILURE("PIAA lens surface 1 not identified");
            }

            piaacmcopticalsystem.ASPHSURFRarray[piaacmcopticalsystem.elemarrayindex[elem]].mat0 = 100;
            piaacmcopticalsystem.ASPHSURFRarray[piaacmcopticalsystem.elemarrayindex[elem]].mat1 = piaacmcopticaldesign.PIAAmaterial_code;
        }
        //       fprintf(fp,"%02ld  %f    PIAAM1\n", elem, piaacmcopticalsystem.elemZpos[elem]);
        elem++;





        // ------------------- elem 4 opaque mask at reflective PIAA M1  -----------------------
        snprintf(piaacmcopticalsystem.name[elem], STRINGMAXLEN_OPTSYST_ELEMNAME, "opaque mask at PIAA elem 1");
        piaacmcopticalsystem.elemtype[elem] = 1; // opaque mask
        {
            imageID ID;
            load_fits("piaa1mask.fits", "piaa1mask", LOADFITS_ERRMODE_IGNORE, &ID);
            if(ID == -1)
            {
                ID = make_disk("piaa1mask", size, size, 0.5 * size, 0.5 * size,
                               piaacmcopticaldesign.r1lim * beamradpix);
            }
            piaacmcopticalsystem.elemarrayindex[elem] = ID;
        }
        piaacmcopticalsystem.elemZpos[elem] = piaacmcopticalsystem.elemZpos[elem - 1];
        //        fprintf(fp,"%02ld  %f    PIAAM1 edge opaque mask\n", elem, piaacmcopticalsystem.elemZpos[elem]);
        elem++;
    }



    // --------------------  elem 5: focal plane mask ------------------------
    if((variable_ID("PIAACMC_NOFPM")) == -1)
    {
        snprintf(piaacmcopticalsystem.name[elem], STRINGMAXLEN_OPTSYST_ELEMNAME, "post focal plane mask pupil");
        piaacmcopticalsystem.elemtype[elem] = 5; // focal plane mask
        piaacmcopticalsystem.elemarrayindex[elem] = 0;

        printf("=========== MAKE FOCAL PLANE MASK ===========\n");
        savefpm = 0;

        {
            variableID IDv;
            if((IDv = variable_ID("PIAACMC_SAVE_fpm")) != -1)
            {
                savefpm = (int)(data.variable[IDv].value.f + 0.001);
            }
        }

        /// call PIAACMCsimul_mkFocalPlaneMask() to make the focal plane mask
        FUNC_CHECK_RETURN(
            mkFocalPlaneMask(
                "fpmzmap",
                "piaacmcfpm",
                piaacmcparams.focmMode,
                savefpm,
                &(piaacmcopticalsystem.FOCMASKarray[0].fpmID)
            )
        );

        // if -1, this is 1-fpm; otherwise, this is impulse response from single zone

        // TEST

        mk_reim_from_complex("piaacmcfpm", "piaacmcfpm_re", "piaacmcfpm_im", 0);

        FUNC_CHECK_RETURN(
            save_fits("piaacmcfpm_re", "./testdir/test_piaacmcfpm_re.fits")
        );

        FUNC_CHECK_RETURN(
            save_fits("piaacmcfpm_im", "./testdir/test_piaacmcfpm_im.fits")
        );

        FUNC_CHECK_RETURN(
            delete_image_ID("piaacmcfpm_re", DELETE_IMAGE_ERRMODE_WARNING)
        );

        FUNC_CHECK_RETURN(
            delete_image_ID("piaacmcfpm_im", DELETE_IMAGE_ERRMODE_WARNING)
        );



        // zfactor is the zoom factor for the DFT, driving sample resolution in the z direction
        // to allow faster DFTs.  Similar to DFTgridpad.
        piaacmcopticalsystem.FOCMASKarray[0].zfactor = piaacmcopticaldesign.fpzfactor;
        // set the position of the pupil from which the DFT propagates to the FPM
        // NOT the position along the beam of the FPM.  That's OK and intended.
        // For this element, this defines the conjugation of the pupil from which we are computing the DFT
        piaacmcopticalsystem.elemZpos[elem] = piaacmcopticalsystem.elemZpos[elem -
                                              1]; // plane from which FT is done
        //      fprintf(fp,"%02ld  %f    post-focal plane mask pupil\n", elem, piaacmcopticalsystem.elemZpos[elem]);
        printf("=========== FOCAL PLANE MASK : DONE ===========\n");
        elem++;
    }



    // if we're a PIAA
    if(piaacmcopticaldesign.PIAAmode == 1)
    {
        // if the the inverse PIAA is prior to the Lyot stop
        // (invPIAAmode = 0 has no inverse PIAA, = 1 has Lyot stop prior to inverse PIAA)
        // there is no inverse PIAA in the WFIRST design
        if(piaacmcopticaldesign.invPIAAmode == 2) // inv PIAA -> Lyot stops
        {
            // --------------------  elem 8: inv PIAA1 ------------------------
            snprintf(piaacmcopticalsystem.name[elem], STRINGMAXLEN_OPTSYST_ELEMNAME, "invPIAA optics 1");

            if(piaacmcopticaldesign.PIAAmaterial_code == 0) // mirror
            {
                piaacmcopticalsystem.elemtype[elem] = 3;    // reflective PIAA M/L 1
            }
            else
            {
                piaacmcopticalsystem.elemtype[elem] = 4;    // refractive PIAA M/L 1
            }

            piaacmcopticalsystem.elemarrayindex[elem] = 2;
            // put an element at z=0 in conjugation space (conjugate to the pupil)
            piaacmcopticalsystem.elemZpos[elem] = 0.0;
            //          fprintf(fp,"%02ld  %f    invPIAA1\n", elem, piaacmcopticalsystem.elemZpos[elem]);
            elem++;

            // --------------------  elem 9: inv PIAA0 ------------------------
            snprintf(piaacmcopticalsystem.name[elem], STRINGMAXLEN_OPTSYST_ELEMNAME, "invPIAA optics 0");

            if(piaacmcopticaldesign.PIAAmaterial_code == 0) //  mirror
            {
                piaacmcopticalsystem.elemtype[elem] = 3;    // reflective PIAA M/L 0
            }
            else
            {
                piaacmcopticalsystem.elemtype[elem] = 4;    // refractive PIAA M/L 0
            }

            piaacmcopticalsystem.elemarrayindex[elem] = 1;
            // previous element + PIAAsep
            piaacmcopticalsystem.elemZpos[elem] = piaacmcopticaldesign.PIAAsep;
            //         fprintf(fp,"%02ld  %f    invPIAA0\n", elem, piaacmcopticalsystem.elemZpos[elem]);
            elem++;
        }
    }


    // --------------------  Lyot masks  ------------------------
    // add Lyot masks as specified in the design
    for(long i = 0; i < piaacmcopticaldesign.NBLyotStop; i++)
    {
        snprintf(piaacmcopticalsystem.name[elem], STRINGMAXLEN_OPTSYST_ELEMNAME, "Lyot mask %ld", i);
        piaacmcopticalsystem.elemtype[elem] = 1; // Lyot mask
        piaacmcopticalsystem.elemarrayindex[elem] = piaacmcopticaldesign.IDLyotStop[i];
        printf("elem %ld  Lyot mask %ld : %ld\n", elem, i, piaacmcopticaldesign.IDLyotStop[i]);
        piaacmcopticalsystem.elemZpos[elem] =  piaacmcopticaldesign.LyotStop_zpos[i];
        //          fprintf(fp,"%02ld  %f  Lyot Stop %ld\n", elem, piaacmcopticalsystem.elemZpos[elem], i);
        elem++;
    }

    // if we're a PIAA
    if(piaacmcopticaldesign.PIAAmode == 1)
    {
        // add stops for the inverse PIAA
        // not in WFIRST design, skipping...
        if(piaacmcopticaldesign.invPIAAmode == 1) // Lyot masks -> inv PIAA
        {
            // --------------------  elem 8: inv PIAA1 ------------------------
            snprintf(piaacmcopticalsystem.name[elem], STRINGMAXLEN_OPTSYST_ELEMNAME, "invPIAA optics 1");
            if(piaacmcopticaldesign.PIAAmaterial_code == 0) // mirror
            {
                piaacmcopticalsystem.elemtype[elem] = 3;    // reflective PIAA M/L 1
            }
            else
            {
                piaacmcopticalsystem.elemtype[elem] = 4;    // refractive PIAA M/L 1
            }
            piaacmcopticalsystem.elemarrayindex[elem] = 2;
            piaacmcopticalsystem.elemZpos[elem] = 0.0;
            //           fprintf(fp,"%02ld  %f    invPIAA1\n", elem, piaacmcopticalsystem.elemZpos[elem]);
            elem++;

            // --------------------  elem 9: inv PIAA0 ------------------------
            snprintf(piaacmcopticalsystem.name[elem], STRINGMAXLEN_OPTSYST_ELEMNAME, "invPIAA optics 0");
            if(piaacmcopticaldesign.PIAAmaterial_code == 0) //  mirror
            {
                piaacmcopticalsystem.elemtype[elem] = 3;    // reflective PIAA M/L 0
            }
            else
            {
                piaacmcopticalsystem.elemtype[elem] = 4;    // refractive PIAA M/L 0
            }
            piaacmcopticalsystem.elemarrayindex[elem] = 1;
            piaacmcopticalsystem.elemZpos[elem] = piaacmcopticaldesign.PIAAsep;
            //           fprintf(fp,"%02ld  %f    invPIAA0\n", elem, piaacmcopticalsystem.elemZpos[elem]);
            elem++;
        }
    }


    // if we're a PIAA
    if(piaacmcopticaldesign.PIAAmode == 1)
    {
        // --------------------  elem 9: back end mask  ------------------------
        // not in WFIRST design, skipping, but it looks very straightforward

        snprintf(piaacmcopticalsystem.name[elem], STRINGMAXLEN_OPTSYST_ELEMNAME, "back end pupil stop  (rad = %f)",
                 piaacmcopticaldesign.pupoutmaskrad);

        piaacmcopticalsystem.elemtype[elem] = 1;
        {
            imageID ID;
            ID = make_disk("pupoutmask", size, size, 0.5 * size, 0.5 * size,
                           piaacmcopticaldesign.pupoutmaskrad * piaacmcopticaldesign.beamrad / piaacmcopticaldesign.pixscale);
            piaacmcopticalsystem.elemarrayindex[elem] = ID;
        }
        piaacmcopticalsystem.elemZpos[elem] =  piaacmcopticalsystem.elemZpos[elem - 1];
        //     fprintf(fp,"%02ld  %f   back end mask\n", elem, piaacmcopticalsystem.elemZpos[elem]);
        elem++;
    }





    piaacmcopticalsystem.NBelem = elem;
    piaacmcopticalsystem.endmode = 0;


    if(piaacmcparams.PIAACMC_save == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        FILE *fp;

        WRITE_FULLFILENAME(
            fname,
            "%s/conjugations.txt",
            piaacmcparams.piaacmcconfdir
        );
        fp = fopen(fname, "w");
        if(fp == NULL)
        {
            PRINT_ERROR("Cannot create %s", fname);
            abort();
        }

        for(elem=0; elem<piaacmcopticalsystem.NBelem; elem++)
        {
            fprintf(
                fp,
                "%02ld  %f    %s\n",
                elem,
                piaacmcopticalsystem.elemZpos[elem],
                piaacmcopticalsystem.name[elem]
            );
        }

        fclose(fp);
    }

    piaacmcparams.optsystinit = 1;

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

