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


#include "PIAACMCsimul.h"

#include "PIAACMCsimul_computePSF.h"
#include "init_piaacmcopticalsystem.h"
#include "init_piaacmcopticaldesign.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"

#include "PIAAshape/makePIAAshapes.h"




/**
 * ---
 *
 * ## Mode 11: Setup multizone ring mask and Compute polychromatic response to zones, store result in FPMresp

   here we compute how the light propagates from each individual mask zone to the focal plane
   (where each mask zone is completely tranparent)
   *
   */
errno_t exec_multizone_fpm_calib()
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


    // initialize
    FUNC_CHECK_RETURN(
        init_piaacmcopticaldesign(1, fpmradld, centobs0, centobs1, 0, 1)
    );

    printf("piaacmcconfdir     : %s\n", piaacmcparams.piaacmcconfdir);
    printf("SCORINGMASKTYPE    : %d\n", piaacmcparams.SCORINGMASKTYPE);
    printf("PIAACMC_FPMsectors : %d\n", piaacmcparams.PIAACMC_FPMsectors);
    printf("lamda              : %ld nm\n",
           (long)(1.0e9 * piaacmcopticaldesign.lambda + 0.1));
    printf("lamdaB             : %ld \n", (long)(1.0 * piaacmcopticaldesign.lambdaB + 0.1));
    printf("piaacmcopticaldesign.NBrings : %ld\n", piaacmcopticaldesign.NBrings);
    printf("mask rad           : %ld\n",
           (long)(100.0 * piaacmcparams.PIAACMC_MASKRADLD + 0.1));
    printf("computePSF_ResolvedTarget : %d\n",
           piaacmcparams.computePSF_ResolvedTarget);
    printf("computePSF_ResolvedTarget_mode : %d\n",
           piaacmcparams.computePSF_ResolvedTarget_mode);
    printf("piaacmcopticaldesign.fpmmaterial_name : %s\n", piaacmcopticaldesign.fpmmaterial_name);
    printf("piaacmcopticaldesign.nblambda         : %d\n", piaacmcopticaldesign.nblambda);




    imageID ID_FPMresp;
    {   // get the combined focal plane mask response
        // set output filename of the combined focal plane mask response file
        char fname[STRINGMAXLEN_FULLFILENAME];

        PIAACMCsimul_update_fnamedescr_conf();

        WRITE_FULLFILENAME(
            fname,
            "%s/FPMresp%d.%s.fits",
            piaacmcparams.piaacmcconfdir,
            piaacmcparams.SCORINGMASKTYPE,
            piaacmcparams.fnamedescr_conf
        );

        load_fits(fname, "FPMresp", 1, &ID_FPMresp);
    }

    // if it did not exist, create it
    if(ID_FPMresp == -1)
    {
        long mzoffset;
        long mzstep;

        // get the number of tmux threads from cli

        piaacmcparams.PIAACMC_FPMresp_mp = 1;
        // 1: all computations on a single thread

        {
            variableID IDv;
            if((IDv = variable_ID("PIAACMC_FPMresp_mp")) != -1) // multi threaded
            {
                piaacmcparams.PIAACMC_FPMresp_mp = (long) data.variable[IDv].value.f + 0.01;
            }
        }
        printf("PIAACMC_FPMresp_mp = %ld\n", piaacmcparams.PIAACMC_FPMresp_mp);

        printf("------------------------------------- STEP02\n");
        printf("piaacmcopticaldesign.focmNBzone  =  %ld   (%ld)\n", piaacmcopticaldesign.focmNBzone,
               piaacmcopticaldesign.NBrings);
        fflush(stdout);

        // get our tmux thread number in [0 PIAACMC_FPMresp_mp]
        // where the master thread has PIAACMC_FPMresp_thread == PIAACMC_FPMresp_mp
        piaacmcparams.PIAACMC_FPMresp_thread = 0;
        {
            variableID IDv;
            if((IDv = variable_ID("PIAACMC_FPMresp_thread")) != -1) // multi threaded
            {
                piaacmcparams.PIAACMC_FPMresp_thread = (long) data.variable[IDv].value.f +
                                                       0.01;
            }
        }
        printf("PIAACMC_FPMresp_thread = %ld\n",
               piaacmcparams.PIAACMC_FPMresp_thread);



        //index = 0;
        if((piaacmcparams.PIAACMC_FPMresp_mp == 1)
                || (piaacmcparams.PIAACMC_FPMresp_thread >
                    piaacmcparams.PIAACMC_FPMresp_mp - 1)) // main or combine process
            // why not test PIAACMC_FPMresp_thread == PIAACMC_FPMresp_mp?
        {
            // we're the parent set up the FPM zone map
            piaacmcparams.FORCE_CREATE_fpmzmap = 1;
            piaacmcparams.FORCE_CREATE_fpmzt = 1;
            piaacmcparams.FORCE_CREATE_fpmza = 1;

            FUNC_CHECK_RETURN(
                init_piaacmcopticaldesign(1, fpmradld, centobs0, centobs1, 0, 1)
            );
        }
        else
        {
            printf("NO initOK file created\n");
            // we're a child tmux thread, do not set up the FPM zone map, get it from parent via file
            piaacmcparams.FORCE_CREATE_fpmzmap = 0;
            piaacmcparams.FORCE_CREATE_fpmzt = 0;
            piaacmcparams.FORCE_CREATE_fpmza = 0;

            FUNC_CHECK_RETURN(
                init_piaacmcopticaldesign(1, fpmradld, centobs0, centobs1, 0, 1)
            );
        }



        //char fname[STRINGMAXLEN_FULLFILENAME];


        char fname1[STRINGMAXLEN_FULLFILENAME];
        char fname2[STRINGMAXLEN_FULLFILENAME];
        imageID IDcomb;

        // if we're the parent load
        if((piaacmcparams.PIAACMC_FPMresp_mp == 1)
                || (piaacmcparams.PIAACMC_FPMresp_thread >
                    piaacmcparams.PIAACMC_FPMresp_mp - 1))
        {
            PIAACMCsimul_update_fnamedescr_conf();
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(
                fname,
                "%s/FPMresp%d.%s.fits",
                piaacmcparams.piaacmcconfdir,
                piaacmcparams.SCORINGMASKTYPE,
                piaacmcparams.fnamedescr_conf
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
                piaacmcparams.piaacmcconfdir,
                piaacmcparams.SCORINGMASKTYPE,
                piaacmcparams.fnamedescr_conf
            );

            // partial FPMresp file
            PIAACMCsimul_update_fnamedescr_conf();
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(
                fname,
                "%s/FPMresp%d.%s.mp%02ld_thread%02ld.fits",
                piaacmcparams.piaacmcconfdir,
                piaacmcparams.SCORINGMASKTYPE,
                piaacmcparams.fnamedescr_conf,
                piaacmcparams.PIAACMC_FPMresp_mp,
                piaacmcparams.PIAACMC_FPMresp_thread
            );

            // stash the filename of the partial file for later
            WRITE_FILENAME(fname1, "%s.tmp", fname);
            WRITE_FILENAME(fname2, "%s", fname);
            // set region of the partial file that this child computes
            mzoffset = piaacmcparams.PIAACMC_FPMresp_thread;
            mzstep = piaacmcparams.PIAACMC_FPMresp_mp;

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
            //   printf("piaacmcopticaldesign.focmNBzone  =  %ld   (%ld)\n", piaacmcopticaldesign.focmNBzone, piaacmcopticaldesign.NBrings);
            //  sleep(3);

            // usual initialzation
            piaacmcopticalsystem.FOCMASKarray[0].mode =
                1; // use 1-fpm computation in optical propagation
            //piaacmcopticaldesign.fpmaskamptransm = 1.0;
            // set the physical size of the FPM as mean(lambda/D)*mask radius in units of lambda/D
            piaacmcopticaldesign.fpmRad = 0.5 * (piaacmcparams.LAMBDASTART +
                                                 piaacmcparams.LAMBDAEND) * piaacmcopticaldesign.Fratio *
                                          piaacmcparams.PIAACMC_MASKRADLD; // piaacmcparams.PIAACMC_MASKRADLD l/D radius at central lambda
            // initialize the optical system
            FUNC_CHECK_RETURN(
                init_piaacmcopticaldesign(1, fpmradld, centobs0, centobs1, 0, 0)
            );


            //     printf("-------------------------- STEP 0005a  piaacmcopticaldesign.focmNBzone  =  %ld   (%ld)\n", piaacmcopticaldesign.focmNBzone, piaacmcopticaldesign.NBrings);
            //            sleep(3);

            // computes or loads the piaa optics from the piaacmc structure
            FUNC_CHECK_RETURN(makePIAAshapes());


            //   printf("-------------------------- STEP 0005b  piaacmcopticaldesign.focmNBzone  =  %ld   (%ld)\n", piaacmcopticaldesign.focmNBzone, piaacmcopticaldesign.NBrings);
            //          sleep(3);

            // initialize the optical system to be on axis
            FUNC_CHECK_RETURN(init_piaacmcopticalsystem(0.0, 0.0));

            // printf("-------------------------- STEP 0005c  piaacmcopticaldesign.focmNBzone  =  %ld   (%ld)\n", piaacmcopticaldesign.focmNBzone, piaacmcopticaldesign.NBrings);
            //        sleep(3);

            // make the shapes again why?  *****************************************
            FUNC_CHECK_RETURN(makePIAAshapes());

            /*               printf("piaacmcopticaldesign.NBrings =                             %ld\n", piaacmcopticaldesign.NBrings);
                            printf("piaacmcopticaldesign.focmNBzone =                          %ld\n", piaacmcopticaldesign.focmNBzone);
                            printf("piaacmcopticaldesign.nblambda =                            %d\n", piaacmcopticaldesign.nblambda);

                            fflush(stdout);
            sleep(5);*/



            // focmMode controls which part of the FPM is propagated
            // if focmMode is a legal zone index, that zone index is propagated
            // here set focmMode beyond a legal zone index, so all zones are transparent and all
            // light including that which misses the FPM is propagated.
            // Later, we will subtract off the actual zone contributions, which will leave only
            // the light that misses the FPM.
            piaacmcparams.focmMode = data.image[piaacmcopticaldesign.zonezID].md[0].size[0] +
                                     10;  // response for no focal plane mask
            piaacmcopticalsystem.FOCMASKarray[0].mode = 1;
            // use 1-fpm computation in optical propagation
            // To compute the on-axis PSF component for light around the FPM, we first compute the full PSF (no focal plane mask)
            // We will, later on, subtract each zone response to it to get the response to light outside the focal plane mask
            //


            FUNC_CHECK_RETURN(
                PIAACMCsimul_computePSF(
                    0.0, 0.0, 0, piaacmcopticalsystem.NBelem, 0,
                    piaacmcparams.computePSF_ResolvedTarget,
                    piaacmcparams.computePSF_ResolvedTarget_mode, 0, NULL)
            );

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
                // axis 1: zones (mz) - size = data.image[piaacmcopticaldesign.zonezID].md[0].size[0]+1
                // axis 3: lambda (k) - size = piaacmcopticaldesign.nblambda

                // allocate the combined FPMresp 3D array
                // ID is the "imvect" array created by PIAACMCsimul_computePSF and contains the pixels in the
                // evaluation set as a 1D vector (0th index) per wavelength
                FUNC_CHECK_RETURN(
                    create_3Dimage_ID_double("FPMresp",
                                             data.image[IDimvect].md[0].size[0],
                                             piaacmcopticaldesign.focmNBzone + 1,
                                             piaacmcopticaldesign.nblambda,
                                             &IDfpmresp)
                );

                // light outside mask, initially set as full PSF
                for(int k = 0; k < piaacmcopticaldesign.nblambda; k++) // loop over wavelengths
                {
                    for(uint32_t ii = 0; ii < data.image[IDimvect].md[0].size[0]; ii++)
                    {
                        // loop over evaluation points
                        // set the 0th zone to be the light from the above on-axis PSF computation with a
                        // black FPM on the evaluation pixels in "imvect"

                        data.image[IDfpmresp].array.D[k * (piaacmcopticaldesign.focmNBzone + 1)
                                                      *data.image[IDimvect].md[0].size[0] + ii] = data.image[IDimvect].array.F[k *
                                                              data.image[IDimvect].md[0].size[0] + ii];
                    }
                }
            }


            // if we're the parent, combine files
            if(piaacmcparams.PIAACMC_FPMresp_thread > piaacmcparams.PIAACMC_FPMresp_mp - 1)
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
                                   piaacmcparams.piaacmcconfdir
                                  );
                FILE * fpt = fopen(fnamet, "w");
                fclose(fpt);


                // we now wait for the children to complete their partial FPM resp files, then combine them
                for(long thr = 0; thr < piaacmcparams.PIAACMC_FPMresp_mp; thr++)
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
                                           piaacmcparams.piaacmcconfdir,
                                           piaacmcparams.SCORINGMASKTYPE,
                                           piaacmcparams.fnamedescr_conf,
                                           piaacmcparams.PIAACMC_FPMresp_mp,
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
                        FUNC_CHECK_RETURN(
                            delete_image_ID("tmpFPMresp", DELETE_IMAGE_ERRMODE_WARNING)
                        );

                        // try to load the final child's partial FPMresp file
                        FUNC_CHECK_RETURN(
                            load_fits(fname, "tmpFPMresp", LOADFITS_ERRMODE_WARNING, &ID1)
                        );
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

                         printf("piaacmcopticaldesign.NBrings =                             %ld\n", piaacmcopticaldesign.NBrings);
                         printf("piaacmcopticaldesign.focmNBzone =                          %ld\n", piaacmcopticaldesign.focmNBzone);
                         printf("piaacmcopticaldesign.nblambda =                            %d\n", piaacmcopticaldesign.nblambda);
                         printf("data.image[ID].md[0].size[0] =                   %ld\n", data.image[ID].md[0].size[0]);
                         printf("data.image[piaacmcopticaldesign.zonezID].md[0].size[0] =   %ld\n", data.image[piaacmcopticaldesign.zonezID].md[0].size[0]);
                         fflush(stdout);
                         sleep(100);
                    */


                    // ID1 is now != -1

                    mzstep = piaacmcparams.PIAACMC_FPMresp_mp; // total number of tmux threads
                    mzoffset = thr; // the thread number of the child that just delivered its result
                    // insert the partial result of child thr into the combined FPMresp array
                    // be sure to skip the first line 'cause we already set it to be the light the went around the FPM
                    for(long mz = 1 + mzoffset; mz < piaacmcopticaldesign.focmNBzone + 1;
                            mz += mzstep) // loop over zone, do every PIAACMC_FPMresp_mp line
                    {

                        printf("mz = %ld    %ld %ld\n", mz, IDfpmresp, ID1);
                        fflush(stdout);
                        for(int k = 0; k < piaacmcopticaldesign.nblambda; k++)
                        {   // for each wavelenth
                            for(uint32_t ii = 0; ii < data.image[ID_FPMresp].md[0].size[0]; ii++)
                            {   // for each evaluation point
                                // index of this evaluation point and wavelength and zone
                                // tmpl1 = k*(nzones+1)*nEvaluationPoints) + zoneIndex*nEvaluationPoints + evaluationPoint
                                long tmpl1 =
                                    k * (data.image[piaacmcopticaldesign.zonezID].md[0].size[0] + 1) * data.image[ID_FPMresp].md[0].size[0]
                                    + mz * data.image[ID_FPMresp].md[0].size[0] + ii;

                                // set the combined array value from the partial file (both are same shape and size, of course)
                                data.image[IDfpmresp].array.D[tmpl1] = data.image[ID1].array.D[tmpl1];

                                // subtract the current zone value from the first zone line, which contained all light
                                // (with no mask).  Eventually this will contain only light that misses the FPM.
                                data.image[IDfpmresp].array.D[k * (data.image[piaacmcopticaldesign.zonezID].md[0].size[0]
                                                                   + 1)*data.image[ID_FPMresp].md[0].size[0] + ii] -= data.image[ID1].array.D[tmpl1];
                            }
                        }
                    }

                    // we're done with the partial array, so delete it
                    FUNC_CHECK_RETURN(
                        delete_image_ID("tmpFPMresp", DELETE_IMAGE_ERRMODE_WARNING)
                    );

                }
                // write out the current state of the combined FPMresp file
                {
                    PIAACMCsimul_update_fnamedescr_conf();
                    char fname[STRINGMAXLEN_FULLFILENAME];
                    WRITE_FULLFILENAME(
                        fname,
                        "%s/FPMresp%d.%s.fits",
                        piaacmcparams.piaacmcconfdir,
                        piaacmcparams.SCORINGMASKTYPE,
                        piaacmcparams.fnamedescr_conf
                    );

                    save_fits("FPMresp", fname);
                }
                // remove the child's .tmp file just in case (it should no longer exist 'cause we renamed, not copied, the .tmp file)
                EXECUTE_SYSTEM_COMMAND("rm %s/FPMresp*.fits.tmp",
                                       piaacmcparams.piaacmcconfdir);
            }
            else // we're a tmux child, so compute the response for our portion
            {
                // name of the child partial FPMresp file (to become .tmp)
                PIAACMCsimul_update_fnamedescr_conf();
                char fname[STRINGMAXLEN_FULLFILENAME];
                WRITE_FULLFILENAME(
                    fname,
                    "%s/FPMresp%d.%s.mp%02ld_thread%02ld.fits",
                    piaacmcparams.piaacmcconfdir,
                    piaacmcparams.SCORINGMASKTYPE,
                    piaacmcparams.fnamedescr_conf,
                    piaacmcparams.PIAACMC_FPMresp_mp,
                    piaacmcparams.PIAACMC_FPMresp_thread
                );

                // diagnostic file to make sure the child is working with the right zones
                {
                    char fnametmp[STRINGMAXLEN_FULLFILENAME];
                    WRITE_FULLFILENAME(
                        fnametmp,
                        "%s/fpmzmap_thread%02ld.fits",
                        piaacmcparams.piaacmcconfdir,
                        piaacmcparams.PIAACMC_FPMresp_thread
                    );
                    FUNC_CHECK_RETURN(save_fits("fpmzmap", fnametmp));
                }

                printf("Making component %ld / %ld\n",
                       piaacmcparams.PIAACMC_FPMresp_thread,
                       piaacmcparams.PIAACMC_FPMresp_mp);
                fflush(stdout);
                piaacmcparams.WRITE_OK = 0;
                // for each FPM zone, compute the response
                // skip the first one 'cause it is not computed by the children
                for(long mz = 1 + mzoffset; mz < piaacmcopticaldesign.focmNBzone + 1; mz += mzstep)
                {
                    piaacmcparams.focmMode =
                        mz;  // focmMode can be a zone index, in which case operations are on that zone
                    piaacmcopticalsystem.FOCMASKarray[0].mode = 0; // direct focal plane mask response

                    // default the reference to the 4th element
                    // but look for the element called "opaque mask at PIAA elem 1"
                    long elem0 = 4;
                    for(long elem = 0; elem < piaacmcopticalsystem.NBelem; elem++)
                    {
                        if(strcmp("opaque mask at PIAA elem 1", piaacmcopticalsystem.name[elem]) == 0)
                        {
                            elem0 = elem;
                            printf("opaque mask at PIAA elem 1 = %ld\n", elem);
                        }
                        // raise an alarm if "opaque mask at PIAA elem 1" is not found *************************************
                    }

                    piaacmcopticalsystem.keepMem[elem0] = 1; // keep it in memory

                    printf("piaacmcopticaldesign.NBrings =                             %ld\n",
                           piaacmcopticaldesign.NBrings);
                    printf("piaacmcopticaldesign.focmNBzone =                          %ld\n",
                           piaacmcopticaldesign.focmNBzone);
                    printf("piaacmcopticaldesign.nblambda =                            %d\n",
                           piaacmcopticaldesign.nblambda);
                    printf("data.image[ID].md[0].size[0] =                   %ld\n",
                           (long) data.image[ID_FPMresp].md[0].size[0]);
                    printf("data.image[piaacmcopticaldesign.zonezID].md[0].size[0] =   %ld\n",
                           (long) data.image[piaacmcopticaldesign.zonezID].md[0].size[0]);
                    fflush(stdout);

                    // compute the on-axis PSF
                    FUNC_CHECK_RETURN(
                        PIAACMCsimul_computePSF(
                            0.0, 0.0, elem0, piaacmcopticalsystem.NBelem, 0,
                            piaacmcparams.computePSF_ResolvedTarget,
                            piaacmcparams.computePSF_ResolvedTarget_mode, 0, NULL)
                    );

                    // The PSF result for the evaluation points is put in array "imvect" which previously was
                    // assigned to another PSF result.
                    // Should put another ID = image_ID("imvect") here *************************************

                    // set the response of this zone from the PSF result
                    for(int k = 0; k < piaacmcopticaldesign.nblambda; k++)
                    {   // loop over wavelength
                        for(uint32_t ii = 0; ii < data.image[ID_FPMresp].md[0].size[0]; ii++)
                        {   // loop over evaluation points
                            // see previous example for explanation of indexing
                            // save response, which is just the value of the on-axis PSF at each evaluation point
                            data.image[IDfpmresp].array.D[k * (data.image[piaacmcopticaldesign.zonezID].md[0].size[0]
                                                               + 1)*data.image[ID_FPMresp].md[0].size[0] + mz * data.image[ID_FPMresp].md[0].size[0] + ii] =
                                                                   data.image[ID_FPMresp].array.F[k * data.image[ID_FPMresp].md[0].size[0] + ii];

                            if(piaacmcparams.PIAACMC_FPMresp_mp == 1)
                            {   // if we're single threaded (no children)
                                // subtract the current zone value from the first zone line, which contained all light
                                // (with no mask).  Eventually this will contain only light that misses the FPM.

                                data.image[IDfpmresp].array.D[k * (data.image[piaacmcopticaldesign.zonezID].md[0].size[0]
                                                                   + 1)*data.image[ID_FPMresp].md[0].size[0] + ii] -= data.image[ID_FPMresp].array.F[k *
                                                                           data.image[ID_FPMresp].md[0].size[0] + ii];
                            }
                        }
                    }


                    printf("Saving FPMresp (ID = %ld) as \"%s\" ...", image_ID("FPMresp"), fname1);
                    fflush(stdout);
                    // fname1 is the .tmp name
                    FUNC_CHECK_RETURN(save_fits("FPMresp", fname1));
                    printf("Done \n");
                    fflush(stdout);
                }

                //*************************** make up our minds about single threading, in the meantime say we don't support it

                // partial file complete!  move it to the final file name so parent can see it
                printf("Saving FPMresp (ID = %ld) as \"%s\" ...", image_ID("FPMresp"), fname2);
                fflush(stdout);
                // fname2 is the final name
                FUNC_CHECK_RETURN(save_fits("FPMresp", fname2));
                printf("Done \n");
                fflush(stdout);
                piaacmcparams.WRITE_OK = 0;

                //   EXECUTE_SYSTEM_COMMAND("mv %s %s", fname1, fname2);
            }
        }
        else
        {
            printf("File FPMresp exists\n");
            //printf("File \"%s\" or \"%s\" exists\n", fname, fnamecomb);
        }
    }
    piaacmcparams.focmMode = -1;


    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}




