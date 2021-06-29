/**
 * @file    PIAAACMCsimul_exec_multizone_fpm_calib.c
 * @brief   Setup and calibrate response of multizone focal plane mask
 *
 *
 */



// System includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"



extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTSYST *optsyst;

extern OPTPIAACMCDESIGN *piaacmc;






/**
 * ---
 *
 * ## Mode 11: Setup multizone ring mask and Compute polychromatic response to zones, store result in FPMresp

   here we compute how the light propagates from each individual mask zone to the focal plane
   (where each mask zone is completely tranparent)
   *
   */
errno_t PIAACMCsimul_exec_multizone_fpm_calib()
{
    DEBUG_TRACE_FSTART();

    double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;

    // spreading computations over multiple processes for resp matrix



    printf("=================================== mode 011 ===================================\n");

    // get cli variables
    {
        variableID IDv;

        if((IDv = variable_ID("PIAACMC_centobs0")) != -1)
        {
            centobs0 = data.variable[IDv].value.f;
        }
        if((IDv = variable_ID("PIAACMC_centobs1")) != -1)
        {
            centobs1 = data.variable[IDv].value.f;
        }
        if((IDv = variable_ID("PIAACMC_fpmradld")) != -1)
        {
            fpmradld = data.variable[IDv].value.f;
            printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
        }
    }
    /*
        if((IDv = variable_ID("PIAACMC_nblambda")) != -1)
        {
            tmpnblambda = data.variable[IDv].value.f;
        }

        if((IDv = variable_ID("PIAACMC_NBrings")) != -1)
        {
            tmpNBrings = data.variable[IDv].value.f;
        }
    */
    // initialize
    if(PIAACMCsimul_initpiaacmcconf(1, fpmradld, centobs0, centobs1, 0, 1) != RETURN_SUCCESS)
    {
        FUNC_RETURN_FAILURE("Call to PIAACMCsimul_initpiaacmcconf failed");
    }

    printf("piaacmcconfdir     : %s\n", piaacmcsimul_var.piaacmcconfdir);
    printf("SCORINGMASKTYPE    : %d\n", piaacmcsimul_var.SCORINGMASKTYPE);
    printf("PIAACMC_FPMsectors : %d\n", piaacmcsimul_var.PIAACMC_FPMsectors);
    printf("lamda              : %ld nm\n",
           (long)(1.0e9 * piaacmc[0].lambda + 0.1));
    printf("lamdaB             : %ld \n", (long)(1.0 * piaacmc[0].lambdaB + 0.1));
    printf("piaacmc[0].NBrings : %ld\n", piaacmc[0].NBrings);
    printf("mask rad           : %ld\n",
           (long)(100.0 * piaacmcsimul_var.PIAACMC_MASKRADLD + 0.1));
    printf("computePSF_ResolvedTarget : %d\n",
           piaacmcsimul_var.computePSF_ResolvedTarget);
    printf("computePSF_ResolvedTarget_mode : %d\n",
           piaacmcsimul_var.computePSF_ResolvedTarget_mode);
    printf("piaacmc[0].fpmmaterial_name : %s\n", piaacmc[0].fpmmaterial_name);
    printf("piaacmc[0].nblambda         : %d\n", piaacmc[0].nblambda);




    imageID ID_FPMresp;
    {   // get the combined focal plane mask response
        // set output filename of the combined focal plane mask response file
        char fname[STRINGMAXLEN_FULLFILENAME];

        PIAACMCsimul_update_fnamedescr_conf();

        WRITE_FULLFILENAME(
            fname,
            "%s/FPMresp%d.%s.fits",
            piaacmcsimul_var.piaacmcconfdir,
            piaacmcsimul_var.SCORINGMASKTYPE,
            piaacmcsimul_var.fnamedescr_conf
        );

        load_fits(fname, "FPMresp", 1, &ID_FPMresp);
    }

    // if it did not exist, create it
    if(ID_FPMresp == -1)
    {
        long mzoffset;
        long mzstep;

        // get the number of tmux threads from cli

        piaacmcsimul_var.PIAACMC_FPMresp_mp = 1;
        // 1: all computations on a single thread

        {
            variableID IDv;
            if((IDv = variable_ID("PIAACMC_FPMresp_mp")) != -1) // multi threaded
            {
                piaacmcsimul_var.PIAACMC_FPMresp_mp = (long) data.variable[IDv].value.f + 0.01;
            }
        }
        printf("PIAACMC_FPMresp_mp = %ld\n", piaacmcsimul_var.PIAACMC_FPMresp_mp);

        printf("------------------------------------- STEP02\n");
        printf("piaacmc[0].focmNBzone  =  %ld   (%ld)\n", piaacmc[0].focmNBzone,
               piaacmc[0].NBrings);
        fflush(stdout);

        // get our tmux thread number in [0 PIAACMC_FPMresp_mp]
        // where the master thread has PIAACMC_FPMresp_thread == PIAACMC_FPMresp_mp
        piaacmcsimul_var.PIAACMC_FPMresp_thread = 0;
        {
            variableID IDv;
            if((IDv = variable_ID("PIAACMC_FPMresp_thread")) != -1) // multi threaded
            {
                piaacmcsimul_var.PIAACMC_FPMresp_thread = (long) data.variable[IDv].value.f +
                        0.01;
            }
        }
        printf("PIAACMC_FPMresp_thread = %ld\n",
               piaacmcsimul_var.PIAACMC_FPMresp_thread);



        //index = 0;
        if((piaacmcsimul_var.PIAACMC_FPMresp_mp == 1)
                || (piaacmcsimul_var.PIAACMC_FPMresp_thread >
                    piaacmcsimul_var.PIAACMC_FPMresp_mp - 1)) // main or combine process
            // why not test PIAACMC_FPMresp_thread == PIAACMC_FPMresp_mp?
        {
            // we're the parent set up the FPM zone map
            piaacmcsimul_var.FORCE_CREATE_fpmzmap = 1;
            piaacmcsimul_var.FORCE_CREATE_fpmzt = 1;
            piaacmcsimul_var.FORCE_CREATE_fpmza = 1;

            if(PIAACMCsimul_initpiaacmcconf(1, fpmradld, centobs0, centobs1, 0, 1) != RETURN_SUCCESS)
            {
                FUNC_RETURN_FAILURE("Call to PIAACMCsimul_initpiaacmcconf failed");
            }
        }
        else
        {
            printf("NO initOK file created\n");
            // we're a child tmux thread, do not set up the FPM zone map, get it from parent via file
            piaacmcsimul_var.FORCE_CREATE_fpmzmap = 0;
            piaacmcsimul_var.FORCE_CREATE_fpmzt = 0;
            piaacmcsimul_var.FORCE_CREATE_fpmza = 0;

            if(PIAACMCsimul_initpiaacmcconf(1, fpmradld, centobs0, centobs1, 0, 1) != RETURN_SUCCESS)
            {
                FUNC_RETURN_FAILURE("Call to PIAACMCsimul_initpiaacmcconf failed");
            }
        }



        //char fname[STRINGMAXLEN_FULLFILENAME];


        char fname1[STRINGMAXLEN_FULLFILENAME];
        char fname2[STRINGMAXLEN_FULLFILENAME];
        imageID IDcomb;

        // if we're the parent load
        if((piaacmcsimul_var.PIAACMC_FPMresp_mp == 1)
                || (piaacmcsimul_var.PIAACMC_FPMresp_thread >
                    piaacmcsimul_var.PIAACMC_FPMresp_mp - 1))
        {
            PIAACMCsimul_update_fnamedescr_conf();
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(
                fname,
                "%s/FPMresp%d.%s.fits",
                piaacmcsimul_var.piaacmcconfdir,
                piaacmcsimul_var.SCORINGMASKTYPE,
                piaacmcsimul_var.fnamedescr_conf
            );

            WRITE_FILENAME(fname1, "%s.tmp", fname);

            WRITE_FILENAME(fname2, "%s", fname);

            mzoffset = 0;
            mzstep = 1;

            // this will always fail in the current state (see line 6606) ************************
            load_fits(fname, "FPMresp", 1, &ID_FPMresp);
            IDcomb = ID_FPMresp;
        }
        else // we're a child tmux thread.
        {
            // combined FPMresp file
            PIAACMCsimul_update_fnamedescr_conf();
            char fnamecomb[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(
                fnamecomb,
                "%s/FPMresp%d.%s.fits",
                piaacmcsimul_var.piaacmcconfdir,
                piaacmcsimul_var.SCORINGMASKTYPE,
                piaacmcsimul_var.fnamedescr_conf
            );

            // partial FPMresp file
            PIAACMCsimul_update_fnamedescr_conf();
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(
                fname,
                "%s/FPMresp%d.%s.mp%02ld_thread%02ld.fits",
                piaacmcsimul_var.piaacmcconfdir,
                piaacmcsimul_var.SCORINGMASKTYPE,
                piaacmcsimul_var.fnamedescr_conf,
                piaacmcsimul_var.PIAACMC_FPMresp_mp,
                piaacmcsimul_var.PIAACMC_FPMresp_thread
            );

            // stash the filename of the partial file for later
            WRITE_FILENAME(fname1, "%s.tmp", fname);
            WRITE_FILENAME(fname2, "%s", fname);
            // set region of the partial file that this child computes
            mzoffset = piaacmcsimul_var.PIAACMC_FPMresp_thread;
            mzstep = piaacmcsimul_var.PIAACMC_FPMresp_mp;

            // may exist from a previous execution with restart
            load_fits(fname, "FPMresp", LOADFITS_ERRMODE_WARNING, &ID_FPMresp);

            // will always fail
            load_fits(fnamecomb, "FPMresp", LOADFITS_ERRMODE_WARNING, &IDcomb);
        }
        // at this point IDcomb==-1, and in the parent ID==-1 always, and in the child ID==-1 if this is not a restart
        // actually create the FPMresp file either as a part by a child or combined by the parent

        if((IDcomb == -1)
                && (ID_FPMresp == -1))
        {
            // this will always fire for the parent thread,
            // and will always fire for children in a fresh run

            //                printf("--------------------------------------------------------STEP 0005 File \"%s\" does not exist: creating\n", fname);
            //   printf("piaacmc[0].focmNBzone  =  %ld   (%ld)\n", piaacmc[0].focmNBzone, piaacmc[0].NBrings);
            //  sleep(3);

            // usual initialzation
            optsyst[0].FOCMASKarray[0].mode =
                1; // use 1-fpm computation in optical propagation
            //piaacmc[0].fpmaskamptransm = 1.0;
            // set the physical size of the FPM as mean(lambda/D)*mask radius in units of lambda/D
            piaacmc[0].fpmRad = 0.5 * (piaacmcsimul_var.LAMBDASTART +
                                       piaacmcsimul_var.LAMBDAEND) * piaacmc[0].Fratio *
                                piaacmcsimul_var.PIAACMC_MASKRADLD; // piaacmcsimul_var.PIAACMC_MASKRADLD l/D radius at central lambda
            // initialize the optical system
            if(PIAACMCsimul_initpiaacmcconf(1, fpmradld, centobs0, centobs1, 0, 0) != RETURN_SUCCESS)
            {
                FUNC_RETURN_FAILURE("Call to PIAACMCsimul_initpiaacmcconf failed");
            }


            //     printf("-------------------------- STEP 0005a  piaacmc[0].focmNBzone  =  %ld   (%ld)\n", piaacmc[0].focmNBzone, piaacmc[0].NBrings);
            //            sleep(3);

            // computes or loads the piaa optics from the piaacmc structure
            if(PIAACMCsimul_makePIAAshapes(piaacmc, 0) != RETURN_SUCCESS)
            {
                FUNC_RETURN_FAILURE("Call to PIAACMCsimul_makePIAAshapes failed");
            }


            //   printf("-------------------------- STEP 0005b  piaacmc[0].focmNBzone  =  %ld   (%ld)\n", piaacmc[0].focmNBzone, piaacmc[0].NBrings);
            //          sleep(3);

            // initialize the optical system to be on axis
            if(PIAACMCsimul_init(piaacmc, 0, 0.0, 0.0) != RETURN_SUCCESS)
            {
                FUNC_RETURN_FAILURE("Call to PIAACMCsimul_init failed");
            }
            // printf("-------------------------- STEP 0005c  piaacmc[0].focmNBzone  =  %ld   (%ld)\n", piaacmc[0].focmNBzone, piaacmc[0].NBrings);
            //        sleep(3);

            // make the shapes again why?  *****************************************
            if(PIAACMCsimul_makePIAAshapes(piaacmc, 0) != RETURN_SUCCESS)
            {
                FUNC_RETURN_FAILURE("Call to PIAACMCsimul_makePIAAshapes failed");
            }

            /*               printf("piaacmc[0].NBrings =                             %ld\n", piaacmc[0].NBrings);
                            printf("piaacmc[0].focmNBzone =                          %ld\n", piaacmc[0].focmNBzone);
                            printf("piaacmc[0].nblambda =                            %d\n", piaacmc[0].nblambda);

                            fflush(stdout);
            sleep(5);*/



            // focmMode controls which part of the FPM is propagated
            // if focmMode is a legal zone index, that zone index is propagated
            // here set focmMode beyond a legal zone index, so all zones are transparent and all
            // light including that which misses the FPM is propagated.
            // Later, we will subtract off the actual zone contributions, which will leave only
            // the light that misses the FPM.
            piaacmcsimul_var.focmMode = data.image[piaacmc[0].zonezID].md[0].size[0] +
                                        10;  // response for no focal plane mask
            optsyst[0].FOCMASKarray[0].mode =
                1; // use 1-fpm computation in optical propagation
            // To compute the on-axis PSF component for light around the FPM, we first compute the full PSF (no focal plane mask)
            // We will, later on, subtract each zone response to it to get the response to light outside the focal plane mask
            //


            {
                double val;
                errno_t fret =
                    PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0,
                                            piaacmcsimul_var.computePSF_ResolvedTarget,
                                            piaacmcsimul_var.computePSF_ResolvedTarget_mode, 0, &val);
                if( fret != RETURN_SUCCESS)
                {
                    FUNC_RETURN_FAILURE("Call to PIAACMCsimul_computePSF failed");
                }

                printf("val = %g\n", val);
            }

            imageID IDfpmresp;
            {
                imageID IDimvect = image_ID("imvect");
                // FPMresp geometry:
                // first dimension (size[0]) is twice the number of evaluation points in the focal plane, giving Re and Im
                //      of the field at that evaluation point
                // second dimension (size[1]) is nbzones+1 zone indices, where nbzones is the number of mask zones (hexagons)
                // +1 because the first zone index stores the response to light that misses the FPM
                // WARNING: FPMresp size[1] is nbzones+1, as first vector stored is the response for light outside the mask
                // third dimension (size[2]) is wavelength


                // axis 0: eval pts (ii) - size = data.image[ID].md[0].size[0]
                // axis 1: zones (mz) - size = data.image[piaacmc[0].zonezID].md[0].size[0]+1
                // axis 3: lambda (k) - size = piaacmc[0].nblambda

                // allocate the combined FPMresp 3D array
                // ID is the "imvect" array created by PIAACMCsimul_computePSF and contains the pixels in the
                // evaluation set as a 1D vector (0th index) per wavelength
                create_3Dimage_ID_double("FPMresp",
                                         data.image[IDimvect].md[0].size[0],
                                         piaacmc[0].focmNBzone + 1,
                                         piaacmc[0].nblambda,
                                         &IDfpmresp);
                //     list_image_ID();
                //    sleep(100);

                // light outside mask, initially set as full PSF
                for(int k = 0; k < piaacmc[0].nblambda; k++) // loop over wavelengths
                {
                    for(uint32_t ii = 0; ii < data.image[IDimvect].md[0].size[0]; ii++)
                    {
                        // loop over evaluation points
                        // set the 0th zone to be the light from the above on-axis PSF computation with a
                        // black FPM on the evaluation pixels in "imvect"

                        data.image[IDfpmresp].array.D[k * (piaacmc[0].focmNBzone + 1)
                                                      *data.image[IDimvect].md[0].size[0] + ii] = data.image[IDimvect].array.F[k *
                                                              data.image[IDimvect].md[0].size[0] + ii];
                    }
                }
            }


            // if we're the parent, combine files
            if(piaacmcsimul_var.PIAACMC_FPMresp_thread > piaacmcsimul_var.PIAACMC_FPMresp_mp - 1)
            {
                long index = 0;
                {
                    variableID IDv;
                    if((IDv = variable_ID("PID")) != -1)
                    {
                        index = (long) data.variable[IDv].value.f + 0.01;
                    }
                }
                // this file is looked for in the bash script, which waits for this
                // file to spawn the tmux child processes
                EXECUTE_SYSTEM_COMMAND("touch initOK_%ld", index);
                //printf("EXECUTING : %s\n", command);

                // now the tmux children have been kicked off by the bash script and are
                // computing their partial FPMresp files

                // begin combining the partial files into the final FPMresp file when ready
                printf("COMBINING FILES\n");
                fflush(stdout);


                char fnamet[STRINGMAXLEN_FULLFILENAME];
                WRITE_FULLFILENAME(fnamet,
                                   "%s/FPMthreadstatus.txt",
                                   piaacmcsimul_var.piaacmcconfdir
                                  );
                FILE * fpt = fopen(fnamet, "w");
                fclose(fpt);


                // we now wait for the children to complete their partial FPM resp files, then combine them
                for(long thr = 0; thr < piaacmcsimul_var.PIAACMC_FPMresp_mp; thr++)
                {
                    printf("thr = %ld\n", thr);
                    fflush(stdout);

                    // each child creates an FPMresp...thread*.fits.tmp file, which is moved to
                    // FPMresp...thread*.fits (set as fname in the next sprintf) when the child is done
                    // signaling to the parent process that this part is ready to ingest.

                    imageID ID1 = -1;
                    while(ID1 == -1) // wait for the partial FPMresp file from each child
                    {
                        // name of final child partial FPMresp file
                        PIAACMCsimul_update_fnamedescr_conf();
                        char fname[STRINGMAXLEN_FULLFILENAME];
                        WRITE_FULLFILENAME(fname,
                                           "%s/FPMresp%d.%s.mp%02ld_thread%02ld.fits",
                                           piaacmcsimul_var.piaacmcconfdir,
                                           piaacmcsimul_var.SCORINGMASKTYPE,
                                           piaacmcsimul_var.fnamedescr_conf,
                                           piaacmcsimul_var.PIAACMC_FPMresp_mp,
                                           thr);

                        printf("Waiting for file \"%s\" ...\n", fname);
                        fflush(stdout);

                        {   // update thread status file
                            FILE * fpt = fopen(fnamet, "a");
                            fprintf(fpt,
                                    "Process %ld (thread %ld) --- Waiting for file \"%s\" ...\n",
                                    (long) getpid(), thr, fname);
                            fclose(fpt);
                            sleep(1.0);
                        }

                        // safely remove image with this name
                        delete_image_ID("tmpFPMresp", DELETE_IMAGE_ERRMODE_WARNING);
                        //  list_image_ID();
                        // try to load the final child's partial FPMresp file
                        load_fits(fname, "tmpFPMresp", 1, &ID1);
                        //  list_image_ID();
                    }
                    // we found this child's partial FPMresp file!
                    // now insert it into our combined FPMresp file

                    {
                        FILE * fpt = fopen(fnamet, "a");
                        fprintf(fpt, "READING %s\n", fnamet);
                        fclose(fpt);
                    }

                    /*     list_image_ID();

                         printf("piaacmc[0].NBrings =                             %ld\n", piaacmc[0].NBrings);
                         printf("piaacmc[0].focmNBzone =                          %ld\n", piaacmc[0].focmNBzone);
                         printf("piaacmc[0].nblambda =                            %d\n", piaacmc[0].nblambda);
                         printf("data.image[ID].md[0].size[0] =                   %ld\n", data.image[ID].md[0].size[0]);
                         printf("data.image[piaacmc[0].zonezID].md[0].size[0] =   %ld\n", data.image[piaacmc[0].zonezID].md[0].size[0]);
                         fflush(stdout);
                         sleep(100);
                    */


                    // ID1 is now != -1

                    mzstep = piaacmcsimul_var.PIAACMC_FPMresp_mp; // total number of tmux threads
                    mzoffset = thr; // the thread number of the child that just delivered its result
                    // insert the partial result of child thr into the combined FPMresp array
                    // be sure to skip the first line 'cause we already set it to be the light the went around the FPM
                    for(long mz = 1 + mzoffset; mz < piaacmc[0].focmNBzone + 1;
                            mz += mzstep) // loop over zone, do every PIAACMC_FPMresp_mp line
                    {

                        printf("mz = %ld    %ld %ld\n", mz, IDfpmresp, ID1);
                        fflush(stdout);
                        for(int k = 0; k < piaacmc[0].nblambda; k++)
                        {   // for each wavelenth
                            for(uint32_t ii = 0; ii < data.image[ID_FPMresp].md[0].size[0]; ii++)
                            {   // for each evaluation point
                                // index of this evaluation point and wavelength and zone
                                // tmpl1 = k*(nzones+1)*nEvaluationPoints) + zoneIndex*nEvaluationPoints + evaluationPoint
                                long tmpl1 =
                                    k * (data.image[piaacmc[0].zonezID].md[0].size[0] + 1) * data.image[ID_FPMresp].md[0].size[0]
                                    + mz * data.image[ID_FPMresp].md[0].size[0] + ii;

                                // set the combined array value from the partial file (both are same shape and size, of course)
                                data.image[IDfpmresp].array.D[tmpl1] = data.image[ID1].array.D[tmpl1];

                                // subtract the current zone value from the first zone line, which contained all light
                                // (with no mask).  Eventually this will contain only light that misses the FPM.
                                data.image[IDfpmresp].array.D[k * (data.image[piaacmc[0].zonezID].md[0].size[0]
                                                                   + 1)*data.image[ID_FPMresp].md[0].size[0] + ii] -= data.image[ID1].array.D[tmpl1];
                            }
                        }
                    }

                    // we're done with the partial array, so delete it
                    delete_image_ID("tmpFPMresp", DELETE_IMAGE_ERRMODE_WARNING);

                }
                // write out the current state of the combined FPMresp file
                {
                    PIAACMCsimul_update_fnamedescr_conf();
                    char fname[STRINGMAXLEN_FULLFILENAME];
                    WRITE_FULLFILENAME(
                        fname,
                        "%s/FPMresp%d.%s.fits",
                        piaacmcsimul_var.piaacmcconfdir,
                        piaacmcsimul_var.SCORINGMASKTYPE,
                        piaacmcsimul_var.fnamedescr_conf
                    );

                    save_fits("FPMresp", fname);
                }
                // remove the child's .tmp file just in case (it should no longer exist 'cause we renamed, not copied, the .tmp file)
                EXECUTE_SYSTEM_COMMAND("rm %s/FPMresp*.fits.tmp",
                                       piaacmcsimul_var.piaacmcconfdir);
            }
            else // we're a tmux child, so compute the response for our portion
            {
                // name of the child partial FPMresp file (to become .tmp)
                PIAACMCsimul_update_fnamedescr_conf();
                char fname[STRINGMAXLEN_FULLFILENAME];
                WRITE_FULLFILENAME(
                    fname,
                    "%s/FPMresp%d.%s.mp%02ld_thread%02ld.fits",
                    piaacmcsimul_var.piaacmcconfdir,
                    piaacmcsimul_var.SCORINGMASKTYPE,
                    piaacmcsimul_var.fnamedescr_conf,
                    piaacmcsimul_var.PIAACMC_FPMresp_mp,
                    piaacmcsimul_var.PIAACMC_FPMresp_thread
                );

                // diagnostic file to make sure the child is working with the right zones
                {
                    char fnametmp[STRINGMAXLEN_FULLFILENAME];
                    WRITE_FULLFILENAME(
                        fnametmp,
                        "%s/fpmzmap_thread%02ld.fits",
                        piaacmcsimul_var.piaacmcconfdir,
                        piaacmcsimul_var.PIAACMC_FPMresp_thread
                    );
                    save_fits("fpmzmap", fnametmp);
                }

                printf("Making component %ld / %ld\n",
                       piaacmcsimul_var.PIAACMC_FPMresp_thread,
                       piaacmcsimul_var.PIAACMC_FPMresp_mp);
                fflush(stdout);
                piaacmcsimul_var.WRITE_OK = 0;
                // for each FPM zone, compute the response
                // skip the first one 'cause it is not computed by the children
                for(long mz = 1 + mzoffset; mz < piaacmc[0].focmNBzone + 1; mz += mzstep)
                {
                    piaacmcsimul_var.focmMode =
                        mz;  // focmMode can be a zone index, in which case operations are on that zone
                    optsyst[0].FOCMASKarray[0].mode = 0; // direct focal plane mask response

                    // default the reference to the 4th element
                    // but look for the element called "opaque mask at PIAA elem 1"
                    long elem0 = 4;
                    for(long elem = 0; elem < optsyst[0].NBelem; elem++)
                    {
                        if(strcmp("opaque mask at PIAA elem 1", optsyst[0].name[elem]) == 0)
                        {
                            elem0 = elem;
                            printf("opaque mask at PIAA elem 1 = %ld\n", elem);
                        }
                        // raise an alarm if "opaque mask at PIAA elem 1" is not found *************************************
                    }

                    optsyst[0].keepMem[elem0] = 1; // keep it in memory

                    printf("piaacmc[0].NBrings =                             %ld\n",
                           piaacmc[0].NBrings);
                    printf("piaacmc[0].focmNBzone =                          %ld\n",
                           piaacmc[0].focmNBzone);
                    printf("piaacmc[0].nblambda =                            %d\n",
                           piaacmc[0].nblambda);
                    printf("data.image[ID].md[0].size[0] =                   %ld\n",
                           (long) data.image[ID_FPMresp].md[0].size[0]);
                    printf("data.image[piaacmc[0].zonezID].md[0].size[0] =   %ld\n",
                           (long) data.image[piaacmc[0].zonezID].md[0].size[0]);
                    fflush(stdout);

                    // compute the on-axis PSF
                    {
                        double val;

                        errno_t fret = PIAACMCsimul_computePSF(0.0, 0.0, elem0, optsyst[0].NBelem, 0,
                                                               piaacmcsimul_var.computePSF_ResolvedTarget,
                                                               piaacmcsimul_var.computePSF_ResolvedTarget_mode, 0, &val);
                        if( fret != RETURN_SUCCESS)
                        {
                            FUNC_RETURN_FAILURE("Call to PIAACMCsimul_computePSF failed");
                        }
                    }
                    // The PSF result for the evaluation points is put in array "imvect" which previously was
                    // assigned to another PSF result.
                    // Should put another ID = image_ID("imvect") here *************************************

                    // set the response of this zone from the PSF result
                    for(int k = 0; k < piaacmc[0].nblambda; k++)
                    {   // loop over wavelength
                        for(uint32_t ii = 0; ii < data.image[ID_FPMresp].md[0].size[0]; ii++)
                        {   // loop over evaluation points
                            // see previous example for explanation of indexing
                            // save response, which is just the value of the on-axis PSF at each evaluation point
                            data.image[IDfpmresp].array.D[k * (data.image[piaacmc[0].zonezID].md[0].size[0]
                                                               + 1)*data.image[ID_FPMresp].md[0].size[0] + mz * data.image[ID_FPMresp].md[0].size[0] + ii] =
                                                                   data.image[ID_FPMresp].array.F[k * data.image[ID_FPMresp].md[0].size[0] + ii];

                            if(piaacmcsimul_var.PIAACMC_FPMresp_mp == 1)
                            {   // if we're single threaded (no children)
                                // subtract the current zone value from the first zone line, which contained all light
                                // (with no mask).  Eventually this will contain only light that misses the FPM.

                                data.image[IDfpmresp].array.D[k * (data.image[piaacmc[0].zonezID].md[0].size[0]
                                                                   + 1)*data.image[ID_FPMresp].md[0].size[0] + ii] -= data.image[ID_FPMresp].array.F[k *
                                                                           data.image[ID_FPMresp].md[0].size[0] + ii];
                            }
                        }
                    }


                    printf("Saving FPMresp (ID = %ld) as \"%s\" ...", image_ID("FPMresp"), fname1);
                    fflush(stdout);
                    // fname1 is the .tmp name
                    save_fits("FPMresp", fname1);
                    printf("Done \n");
                    fflush(stdout);
                }

                //*************************** make up our minds about single threading, in the meantime say we don't support it

                // partial file complete!  move it to the final file name so parent can see it
                printf("Saving FPMresp (ID = %ld) as \"%s\" ...", image_ID("FPMresp"), fname2);
                fflush(stdout);
                // fname2 is the final name
                save_fits("FPMresp", fname2);
                printf("Done \n");
                fflush(stdout);
                piaacmcsimul_var.WRITE_OK = 0;

                //   EXECUTE_SYSTEM_COMMAND("mv %s %s", fname1, fname2);
            }
        }
        else
        {
            printf("File FPMresp exists\n");
            //printf("File \"%s\" or \"%s\" exists\n", fname, fnamecomb);
        }
    }
    piaacmcsimul_var.focmMode = -1;


    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}




