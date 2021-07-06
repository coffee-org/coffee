/**
 * @file    OptSystProp_run.c
 * @brief   Optical propagation execution
 *
 *
 *
 *
 */

#include <stdlib.h>
#include <math.h>


#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "fft/fft.h"

#include "OpticsMaterials/OpticsMaterials.h"
#include "OptSystProp/OptSystProp.h"





/**
 *
 * @brief Optical propagation execution
 *
 * @param[in]   index	system index (usually 0)
 * @param[in]   elemstart	starting element index
 * @param[in]	elemend		ending element index
 * @param[in]   savedir  directory to which image results are saved
 * @param[in]   sharedmem  1 if WF* arrays should be kept in memory after use
 *
 * optsyst.elemkeepmem	1 if element complex amplitude should be kept in memory after use
 */
errno_t OptSystProp_run(OPTSYST    *optsyst,
                        long        index,
                        long        elemstart,
                        long        elemend,
                        const char *savedir,
                        int         sharedmem
                       )
{
    DEBUG_TRACE_FSTART();

    long size;
    long nblambda;
    uint64_t size2;

    double proplim = 1.0e-4;

    char imnameamp_in[STRINGMAXLEN_IMGNAME];
    char imnamepha_in[STRINGMAXLEN_IMGNAME];
    char imnameamp_out[STRINGMAXLEN_IMGNAME];
    char imnamepha_out[STRINGMAXLEN_IMGNAME];

    long emax;

    long elemstart1 = 0;
    int elemOK;


    uint32_t *imsizearray;

    // number of pixels in one side of each square data array
    size = optsyst[index].size;
    size2 = size * size; // area = total number of pixels
    nblambda = optsyst[index].nblambda;

    imsizearray = (uint32_t *) malloc(sizeof(uint32_t) * 3);
    if(imsizearray == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    // create base complex amplitude output arrays

    // set base complex amplitude array sizes
    imsizearray[0] = size;
    imsizearray[1] = size;
    imsizearray[2] = nblambda;

    // WFamp is the standard output of the wave front complex amplitude
    {
        char imname[STRINGMAXLEN_IMGNAME];
        WRITE_IMAGENAME(imname, "WFamp%ld", index);
        imageID IDa = image_ID(imname);
        if(IDa == -1)
        {
            create_image_ID(imname, 3, imsizearray, _DATATYPE_FLOAT, sharedmem, 0, 0, &IDa);
            //    create_3Dimage_ID(imname, size, size, nblambda);
        }

        // initialize wavefront amplitude to 1
        for(uint64_t ii = 0; ii < size2; ii++)
            for(long kl = 0; kl < nblambda; kl++)
            {
                data.image[IDa].array.F[size2 * kl + ii] = 1.0;
            }
    }

    // WFpha is the standard output of the wave front complex phase

    {
        char imname[STRINGMAXLEN_IMGNAME];
        WRITE_IMAGENAME(imname, "WFpha%ld", index);
        imageID IDp = image_ID(imname);
        if(IDp == -1)
        {
            create_image_ID(imname, 3, imsizearray, _DATATYPE_FLOAT, sharedmem, 0, 0, NULL);
            //create_3Dimage_ID(imname, size, size, nblambda);
        }
    }



    // start propagation at the first defined element or the specified input start point

    // set up the data name stubs for later definition of file and data names
    // one for each defined element
    elemstart1 = 0;
    elemOK = 1;
    // looking for the first defined element prior to input elemstart
    // or elemstart, whichever comes first
    while(elemOK == 1)
    {
        long ID1, ID2;

        if(elemstart1 == 0)
        {
            WRITE_IMAGENAME(imnameamp_in, "WFamp%ld", index);
            WRITE_IMAGENAME(imnamepha_in, "WFpha%ld", index);
        }
        else
        {
            WRITE_IMAGENAME(imnameamp_in, "WFamp%ld_%03ld", index, elemstart1 - 1);
            WRITE_IMAGENAME(imnamepha_in, "WFpha%ld_%03ld", index, elemstart1 - 1);
        }
        // skip elements that are not defined and and prior to the input elemstart
        if(((ID1 = image_ID(imnameamp_in)) != -1)
                && ((ID2 = image_ID(imnamepha_in)) != -1) && (elemstart1 < elemstart + 1))
        {
            elemstart1++;
            elemOK = 1;
        }
        else
        {
            elemOK = 0;
        }

        printf("%ld/%ld %d    %s %ld   %s %ld\n", elemstart1, elemstart, elemOK,
               imnameamp_in, ID1, imnamepha_in, ID2);
    }
    // back up so we propagate from the current element to the next element
    elemstart1--;

    printf("STARTING AT ELEMENT %ld\n", elemstart1);

    emax = elemend;
    if(emax > optsyst[index].NBelem)
    {
        emax = optsyst[index].NBelem;    // make sure we're not doing too many elements, so user does not need
    }
    // to know how many elements there are

    // loop over the elements, propagating from elem to elem + 1
    long elem;
    for(elem = elemstart1; elem < emax; elem++)
    {
        double propdist;

        // naming convention: WFamp0/pha0 is the running buffer for propagation
        // WFamp0_element#/pha0_element# persisted saved results for each element
        if(elem == 0)
        {
            // initialize the running buffer
            WRITE_IMAGENAME(imnameamp_in, "WFamp%ld", index);
            WRITE_IMAGENAME(imnamepha_in, "WFpha%ld", index);
            WRITE_IMAGENAME(imnameamp_out, "WFamp%ld_000", index);
            WRITE_IMAGENAME(imnamepha_out, "WFpha%ld_000", index);
            // propagation distance between elements
            propdist = optsyst[index].elemZpos[0];
        }
        else
        {
            // now there is a previous element
            // this is how we move down the elements
            WRITE_IMAGENAME(imnameamp_in, "WFamp%ld_%03ld", index, elem - 1); // already exists
            WRITE_IMAGENAME(imnamepha_in, "WFpha%ld_%03ld", index, elem - 1);
            WRITE_IMAGENAME(imnameamp_out, "WFamp%ld_%03ld", index, elem);
            WRITE_IMAGENAME(imnamepha_out, "WFpha%ld_%03ld", index, elem);
            // propagation distance between elements
            propdist = optsyst[index].elemZpos[elem] - optsyst[index].elemZpos[elem - 1];
        }

        // delete the element data if it exists
        if((image_ID(imnameamp_out) != -1) && (sharedmem == 0))
        {
            delete_image_ID(imnameamp_out, DELETE_IMAGE_ERRMODE_WARNING);
        }

        if((image_ID(imnamepha_out) != -1) && (sharedmem == 0))
        {
            delete_image_ID(imnamepha_out, DELETE_IMAGE_ERRMODE_WARNING);
        }

        // if our propagation distance big enough to require propagation
        // (if they're close enough together then the propagation is a trivial copy)
        if(fabs(propdist) > proplim)
        {

            printf("Propagating to element %ld  (%lf m)\n", elem,  propdist);
            // do the real propagation via Fresnel propagation
            OptSystProp_propagateCube(optsyst, 0, imnameamp_in, imnamepha_in, imnameamp_out,
                                      imnamepha_out, propdist, sharedmem);
        }
        else // do the trivial identity propagation
        {
            copy_image_ID(imnameamp_in, imnameamp_out, sharedmem);
            copy_image_ID(imnamepha_in, imnamepha_out, sharedmem);
        }
        imageID IDa = image_ID(imnameamp_out);
        imageID IDp = image_ID(imnamepha_out);

        /// discard element memory after used
        printf("*********** %ld  -> %d\n", elem - 1, optsyst[index].keepMem[elem - 1]);
        if((optsyst[index].keepMem[elem - 1] == 0) && (sharedmem == 0))
        {
            printf("********** Deleting element %ld      %s %s\n", elem - 1, imnameamp_in,
                   imnamepha_in);
            delete_image_ID(imnameamp_in, DELETE_IMAGE_ERRMODE_WARNING);
            delete_image_ID(imnamepha_in, DELETE_IMAGE_ERRMODE_WARNING);
        }

        printf("Applying element %ld\n", elem);
        fflush(stdout);






        // multiply each pixel by an amplitude
        // this would be 0 or 1 for an opaque mask
        // or between 0 and 1 for an apodization
        if(optsyst[index].elemtype[elem] == 1) // AMPLITUDE MASK
        {
            imageID ID = optsyst[index].elemarrayindex[elem];
            printf("============= elem %ld:  Opaque mask (%s) =================\n", elem,
                   data.image[ID].name);
            fflush(stdout);
            //	list_image_ID();

            if(ID == -1)
            {
                printf("ERROR: ID = -1, missing mask image\n");
                exit(0);
            }

            //	save_fits(data.image[ID].name, "opmask.fits"); //TEST
            //save_fits(data.image[IDa].name, "opmask1.fits"); //TEST

            printf("ID = %ld\n", ID);
            fflush(stdout);

            // achromatic (no wavelength dimension) vs chromatic (has wavelength dimension)
            if((data.image[ID].md[0].naxis == 2)
                    || (data.image[ID].md[0].size[2] != nblambda))
            {
                // chromatic case: apply mask to each wavelength
                //		printf("single dim %ld %ld\n", data.image[ID].md[0].size[2], nblambda);
                //	fflush(stdout);
# ifdef HAVE_LIBGOMP
                #pragma omp parallel default(shared) private(ii)
                {
                    #pragma omp for
# endif
                    for(long kl = 0; kl < nblambda; kl++) // loop over wavelengths
                        for(uint64_t ii = 0; ii < size2; ii++)
                            // actually apply the mask
                        {
                            data.image[IDa].array.F[size2 * kl + ii] *= data.image[ID].array.F[ii];
                        }
# ifdef HAVE_LIBGOMP
                }
# endif
            }
            else
            {
                //	printf("multi dim %ld %ld\n", data.image[ID].md[0].size[2], nblambda);
                //	fflush(stdout);
# ifdef HAVE_LIBGOMP
                #pragma omp parallel
                {
                    #pragma omp for
# endif
                    // don't loop over wavelength
                    for(uint64_t ii = 0; ii < size2 *nblambda; ii++)
                        // actually apply the mask
                    {
                        data.image[IDa].array.F[ii] *= data.image[ID].array.F[ii];
                    }
# ifdef HAVE_LIBGOMP
                }
# endif
            }

            //	save_fits(data.image[IDa].name, "opmask2.fits"); //TEST
            //	printf("POINT 1.1\n");

        }










        // apply a change in phase, which depends on the mirror shape
        // same chromatic vs. achromatic choice as above, but this time the phase data is chromatic
        // (phase is driven by wavelength!!)

        if(optsyst[index].elemtype[elem] == 3)
        {
            // MIRROR SURFACE - STORED AS OPD MAP AS A SINGLE MAP (ACHROMATIC) OR A CUBE (CHROMATIC)
            printf("============= Mirror surface =======================\n");
            fflush(stdout);
            imageID ID = optsyst[index].ASPHSURFMarray[optsyst[index].elemarrayindex[elem]].surfID;
            printf("%d surface ID = %ld\n", optsyst[index].elemarrayindex[elem], ID);

            if(ID == -1)
            {
                printf("ERROR: Surface ID does not exist\n");
                printf("   %s %d\n", __FILE__, __LINE__);
                printf("   Function %s return FAILURE\n", __func__);
                DEBUG_TRACE_FEXIT();
                return RETURN_FAILURE;
            }


            if(data.image[ID].md[0].naxis == 2)
            {
# ifdef HAVE_LIBGOMP
                #pragma omp parallel default(shared) private(ii)
                {
                    #pragma omp for
# endif
                    // in the "achromatic mirror" case, the phase is still chromatic so we
                    // have to loop over wavelength
                    for(long kl = 0; kl < nblambda; kl++)
                        for(uint64_t ii = 0; ii < size2; ii++)
                        {   // compute the change in phase
                            data.image[IDp].array.F[size2 * kl + ii] -= 4.0 * M_PI *
                            data.image[ID].array.F[ii] / optsyst[index].lambdaarray[kl];
                        }
# ifdef HAVE_LIBGOMP
                }
# endif
            }
            else // chromatic "mirror"
            {
# ifdef HAVE_LIBGOMP
                #pragma omp parallel default(shared) private(ii)
                {
                    #pragma omp for
# endif
                    for(long kl = 0; kl < nblambda; kl++)
                        for(uint64_t ii = 0; ii < size2; ii++)
                        {   // compute the change in phase
                            data.image[IDp].array.F[size2 * kl + ii] -= 4.0 * M_PI *
                            data.image[ID].array.F[size2 * kl + ii] / optsyst[index].lambdaarray[kl];
                        }
# ifdef HAVE_LIBGOMP
                }
# endif
            }
        }








        // apply a change in phase
        if(optsyst[index].elemtype[elem] == 4)
        {
            // REFRACTIVE SURFACE - STORED AS SAG MAP AS A SINGLE MAP (ACHROMATIC) OR A CUBE (CHROMATIC)
            printf("============= [%ld] Refractive surface =======================\n", elem);
            fflush(stdout);

            imageID ID = optsyst[index].ASPHSURFRarray[optsyst[index].elemarrayindex[elem]].surfID;
            DEBUG_TRACEPOINT("index %ld  elem %ld  %d surface ID : %ld",
                             index, elem,
                             optsyst[index].elemarrayindex[elem], ID);

            DEBUG_TRACEPOINT("refractive surface ID %ld %s", (long) ID, data.image[ID].md[0].name);

            if(ID == -1)
            {
                list_image_ID();
                FUNC_RETURN_FAILURE("image ID not found");
            }

            list_image_ID();

            if(optsyst[index].ASPHSURFRarray[optsyst[index].elemarrayindex[elem]].init != 1)
            {
                for(long kl = 0; kl < nblambda; kl++)
                {
                    // get the indices of refraction, n0 for ambient index, n1 for the lens index etc.

                    double n0 = OpticsMaterials_n(
                                    optsyst[index].ASPHSURFRarray[optsyst[index].elemarrayindex[elem]].mat0,
                                    optsyst[index].lambdaarray[kl]);

                    double n1 = OpticsMaterials_n(
                                    optsyst[index].ASPHSURFRarray[optsyst[index].elemarrayindex[elem]].mat1,
                                    optsyst[index].lambdaarray[kl]);

                    // set the resulting wavelength-dependent phase coefficient
                    optsyst[index].ASPHSURFRarray[optsyst[index].elemarrayindex[elem]].ncoeff[kl] =
                        2.0 * M_PI * (n0 - n1) / optsyst[index].lambdaarray[kl];

                    DEBUG_TRACEPOINT("elem %ld ASPHSURFRarray ncoeff %ld = %lf",
                                     elem, kl,
                                     optsyst[index].ASPHSURFRarray[optsyst[index].elemarrayindex[elem]].ncoeff[kl]);
                }
                optsyst[index].ASPHSURFRarray[optsyst[index].elemarrayindex[elem]].init = 1;
            }

            //TEST
            {
                char fnamepre[STRINGMAXLEN_FILENAME];
                WRITE_FILENAME(fnamepre, "test_refract_elem%ld_prepha.fits", elem);
                save_fl_fits(data.image[IDp].md[0].name, fnamepre);

                char fnameoptpha[STRINGMAXLEN_FILENAME];
                WRITE_FILENAME(fnameoptpha, "test_refract_elem%ld_optpha.fits", elem);
                save_fl_fits(data.image[ID].md[0].name, fnameoptpha);
            }

            if(data.image[ID].md[0].naxis == 2)
            {
# ifdef HAVE_LIBGOMP
                #pragma omp parallel default(shared) private(ii)
                {
                    #pragma omp for
# endif
                    for(long kl = 0; kl < nblambda; kl++)
                        for(uint64_t ii = 0; ii < size2; ii++)
                        {   // apply change in phase
                            data.image[IDp].array.F[size2 * kl + ii] +=
                            data.image[ID].array.F[ii] * optsyst[index].ASPHSURFRarray[optsyst[index].elemarrayindex[elem]].ncoeff[kl];
                        }
# ifdef HAVE_LIBGOMP
                }
# endif
            }
            else
            {
# ifdef HAVE_LIBGOMP
                #pragma omp parallel default(shared) private(ii)
                {
                    #pragma omp for
# endif
                    for(long kl = 0; kl < nblambda; kl++)
                        for(uint64_t ii = 0; ii < size2; ii++)
                        {   // apply change in phase
                            data.image[IDp].array.F[size2 * kl + ii] +=
                            data.image[ID].array.F[size2 * kl + ii] * optsyst[index].ASPHSURFRarray[optsyst[index].elemarrayindex[elem]].ncoeff[kl];
                        }
# ifdef HAVE_LIBGOMP
                }
# endif
            }
            //TEST
            {
                if(elem==5)
                {
                    for(long kl = 0; kl < nblambda; kl++)
                    for(uint64_t ii = 0; ii < size2; ii++)
                    {
                        data.image[IDp].array.F[size2 * kl + ii] = 0.0;
                    }
                }


                char fnamepost[STRINGMAXLEN_FILENAME];
                WRITE_FILENAME(fnamepost, "test_refract_elem%ld_postpha.fits", elem);
                save_fl_fits(data.image[IDp].md[0].name, fnamepost);
            }
        }









        // covers both transmission and phase FPM
        // apply a change in phase or amplitude
        if(optsyst[index].elemtype[elem] == 5)
        {   // FOCAL PLANE MASK - MASK INPUT IS 1-MASK FOR EFFICIENT DFT
            printf("============= [%ld] Focal Plane Mask ==============\n", elem);
            fflush(stdout);
            // uses 1-fpm

            /* {
                 // TEST: INPUT
                 char fname[STRINGMAXLEN_FULLFILENAME];
                 WRITE_FULLFILENAME(fname, "%s/test_preFPMamp_%02ld.fits", savedir, elem);
                 save_fits(imnameamp_out, fname);
                 WRITE_FULLFILENAME(fname, "%s/test_preFPMpha_%02ld.fits", savedir, elem);
                 save_fits(imnamepha_out, fname);
             }*/


            // convert input pupil amplitude and phase to Re and Im
            imageID ID = mk_complex_from_amph(imnameamp_out, imnamepha_out, "_WFctmp", 0);
            if(sharedmem == 0)
            {
                delete_image_ID(imnameamp_out, DELETE_IMAGE_ERRMODE_WARNING);
                delete_image_ID(imnamepha_out, DELETE_IMAGE_ERRMODE_WARNING);
            }


            printf("optsyst[index].DFTgridpad = %ld\n", optsyst[index].DFTgridpad);

            // if we're subsampling in pixel space for faster DFTs
            if(optsyst[index].DFTgridpad > 0)
            {
                // RESAMPLE ON A SPARSE GRID TO SPEED UP DFT
                // we're really reducing the resolution of the grid (not subsampling)
                // so we need to accumulate the values in the neighborhood of each
                // subsampled point onto that subsampled point
                imageID IDre;
                create_3Dimage_ID("dftgridre", size, size, nblambda, &IDre);

                imageID IDim;
                create_3Dimage_ID("dftgridim", size, size, nblambda, &IDim);

                // grid size, odd number - this is the space between subsampled pixels
                long gsize = 2 * optsyst[index].DFTgridpad + 1;

                // offset from box edge to active pixel
                long offset = optsyst[index].DFTgridpad;

                ID = image_ID("_WFctmp");
                for(long kl = 0; kl < nblambda; kl++)
                    for(uint32_t ii = 0; ii < size; ii++)
                        for(uint32_t jj = 0; jj < size; jj++)
                        {   // for each pixel

                            float re = data.image[ID].array.CF[size2 * kl + jj * size + ii].re;
                            float im = data.image[ID].array.CF[size2 * kl + jj * size + ii].im;

                            // find the nearest subsampled point to ii,jj
                            long ii1 = offset + ((long)(ii / gsize)) * gsize;
                            long jj1 = offset + ((long)(jj / gsize)) * gsize;
                            if((ii1 < size) && (jj1 < size))
                            {
                                // add the value at ii,jj to the nearest subsampled point
                                data.image[IDre].array.F[size2 * kl + jj1 * size + ii1] += re;
                                data.image[IDim].array.F[size2 * kl + jj1 * size + ii1] += im;
                            }
                        }

                // save_fits("dftgridre", "dftgridre.fits");
                // save_fits("dftgridim", "dftgridim.fits");

                // combine separate Re and Im files into single Re/Im complex array
                mk_complex_from_reim("dftgridre", "dftgridim", "_WFctmpc", 0);
                delete_image_ID("dftgridre", DELETE_IMAGE_ERRMODE_WARNING);
                delete_image_ID("dftgridim", DELETE_IMAGE_ERRMODE_WARNING);



                long elemindex = optsyst[index].elemarrayindex[elem];
                ID = optsyst[index].FOCMASKarray[elemindex].fpmID;
                printf("focm # %ld: %s\n", elemindex, data.image[ID].name);


                //      printf("Saving to testfpm.fits\n");
                //      fflush(stdout);


                delete_image_ID("fpmatest", DELETE_IMAGE_ERRMODE_WARNING);
                delete_image_ID("fpmptest", DELETE_IMAGE_ERRMODE_WARNING);
                mk_amph_from_complex("piaacmcfpm", "fpmatest", "fpmptest", 0);

                if(optsyst[index].SAVE == 1)
                {
                    // make amp and phase files from complex
                    mk_amph_from_complex(data.image[ID].name, "fpma", "fpmp", 0);

                    char fname[STRINGMAXLEN_FULLFILENAME];

                    WRITE_FULLFILENAME(fname, "%s/fpm__ampl.fits", savedir);
                    save_fits("fpma", fname);

                    WRITE_FULLFILENAME(fname, "%s/fpm__pha.fits", savedir);
                    save_fits("fpmp", fname);

                    delete_image_ID("fpma", DELETE_IMAGE_ERRMODE_WARNING);
                    delete_image_ID("fpmp", DELETE_IMAGE_ERRMODE_WARNING);
                }


                //	      exit(0);
                /*     list_image_ID();
                     printf("fft_DFTinsertFPM  args :  %s %f\n", data.image[ID].name, optsyst[index].FOCMASKarray[i].zfactor);
                     sleep(10); // TEST*/



                // do DFT from pupil to pupil with an FPM in the middle
                // so this does pupil1 -> DFT -> FPM -> DFT -> pupil2
                // output in last argument
                // the focal plane mask is the second argument, which contains properties of the FPM
                // that are applied as part of fft_DFTinsertFPM
                fft_DFTinsertFPM("_WFctmpc", data.image[ID].name,
                                 optsyst[index].FOCMASKarray[elemindex].zfactor, "_WFcout", NULL);

                delete_image_ID("_WFctmpc", DELETE_IMAGE_ERRMODE_WARNING);




                /*sprintf(command, "mv _DFT_foca %s/_DFT_foca_%02ld.fits", savedir, elem);
                r = system(command);
                sprintf(command, "mv _DFT_focp %s/_DFT_focp_%02ld.fits", savedir, elem);
                r = system(command);
                */


                // TEST
                /*mk_reim_from_complex("_WFcout", "_twfre", "_twfim");
                sprintf(fname, "%s/test_twfre.fits", savedir);
                save_fits("_twfre", fname);
                sprintf(fname, "%s/test_twfim.fits", savedir);
                save_fits("_twfim", fname);
                delete_image_ID("_twfre", DELETE_IMAGE_ERRMODE_WARNING);
                delete_image_ID("_twfim", DELETE_IMAGE_ERRMODE_WARNING);
                */

                //
                // INTERPOLATE SPARSE RESULT ON CONTINUOUS GRID
                //
                // go back to the normal pixel grid via a convolution with a bilinear tent kernel
                // gsize is the grid spacing between subsampled pixels
                // so the kernel spans one subsampled area
                float *convkern;
                convkern = (float *) malloc(sizeof(float) * (2 * gsize + 1) * (2 * gsize + 1));
                if(convkern == NULL) {
                    PRINT_ERROR("malloc returns NULL pointer");
                    abort();
                }
                double tot = 0.0;
                // set up the convolution kernel
                for(long i = 0; i < 2 * gsize + 1; i++)
                    for(long j = 0; j < 2 * gsize + 1; j++)
                    {
                        float u = fabs(1.0 * (i - gsize) / gsize);
                        float t = fabs(1.0 * (j - gsize) / gsize);
                        float val = (1.0 - u) * (1.0 - t);
                        convkern[j * (2 * gsize + 1) + i] = val;
                        //printf("   %d %d %f\n", i, j, val);
                        tot += val;
                    }
                for(long i = 0; i < (2 * gsize + 1) * (2 * gsize + 1); i++)
                    // normalize the kernel to unit integral over area
                {
                    convkern[i] *= gsize * gsize / tot;
                }

                // apply the kernel
                ID = image_ID("_WFcout");

                imageID IDre1;
                create_3Dimage_ID("dftgridre1", size, size, nblambda, &IDre1);

                imageID IDim1;
                create_3Dimage_ID("dftgridim1", size, size, nblambda, &IDim1);

                for(long kl = 0; kl < nblambda; kl++) // for each wavelength
                    for(long ii1 = offset + gsize; ii1 < size - gsize;
                            ii1 += gsize) // for each subsampled pixel
                        for(long jj1 = offset + gsize; jj1 < size - gsize; jj1 += gsize)
                        {
                            // values of the DFT output on the subsampled pixels
                            float re = data.image[ID].array.CF[size2 * kl + jj1 * size + ii1].re;
                            float im = data.image[ID].array.CF[size2 * kl + jj1 * size + ii1].im;

                            for(long i = 0; i < 2 * gsize + 1; i++)
                                for(long j = 0; j < 2 * gsize + 1; j++)
                                {
                                    long ii = ii1 + (i - gsize); // indices of the original pixels to interpolate onto
                                    long jj = jj1 + (j - gsize);
                                    // perform the convolution, saving the results in separate Re and Im arrays
                                    data.image[IDre1].array.F[size2 * kl + jj * size + ii] += re * convkern[j *
                                            (2 * gsize + 1) + i];
                                    data.image[IDim1].array.F[size2 * kl + jj * size + ii] += im * convkern[j *
                                            (2 * gsize + 1) + i];
                                }
                        }

                // TEST

                /*         sprintf(fname, "%s/test_dftgridre1_elem%ld.fits", savedir, elem);
                         save_fits("dftgridre1", fname);
                         sprintf(fname, "%s/test_dftgridim1_elem%ld.fits", savedir, elem);
                         save_fits("dftgridim1", fname);
                  */

                free(convkern);
                delete_image_ID("_WFcout", DELETE_IMAGE_ERRMODE_WARNING);
                // convert the interpolated full-resolution DFT result arrays to a single complex array
                mk_complex_from_reim("dftgridre1", "dftgridim1", "_WFcout", 0);
                delete_image_ID("dftgridre1", DELETE_IMAGE_ERRMODE_WARNING);
                delete_image_ID("dftgridim1", DELETE_IMAGE_ERRMODE_WARNING);

            }
            else /// If not subsampled :
            {
                long i = optsyst[index].elemarrayindex[elem];
                ID = optsyst[index].FOCMASKarray[i].fpmID;
                printf("focm : %s\n", data.image[ID].name);
                fflush(stdout);


                /*{   // TEST: write fpm to disk
                    // make amp and phase files from complex
                    mk_amph_from_complex(data.image[ID].name, "fpma", "fpmp", 0);

                    char fname[STRINGMAXLEN_FULLFILENAME];

                    WRITE_FULLFILENAME(fname, "%s/fpm__ampl.fits", savedir);
                    save_fits("fpma", fname);

                    WRITE_FULLFILENAME(fname, "%s/fpm__pha.fits", savedir);
                    save_fits("fpmp", fname);

                    delete_image_ID("fpma", DELETE_IMAGE_ERRMODE_WARNING);
                    delete_image_ID("fpmp", DELETE_IMAGE_ERRMODE_WARNING);
                }*/


                /// - do DFT from pupil to pupil with a FPM in the middle, so this does pupil1 -> DFT -> FPM -> DFT -> pupil2
                ///
                /// - output in last argument
                ///
                /// - the focal plane mask is the second argument, which contains properties of the FPM that are applied as part of fft_DFTinsertFPM()
                fft_DFTinsertFPM("_WFctmp", data.image[ID].name,
                                 optsyst[index].FOCMASKarray[i].zfactor, "_WFcout", NULL);

                /* {   // TEST: OUPUT 1
                     mk_amph_from_complex("_WFcout", "postFPMamp", "postFPMpha", 0);

                     char fname[STRINGMAXLEN_FULLFILENAME];

                     WRITE_FULLFILENAME(fname, "%s/test_postFPMamp0_%02ld.fits", savedir, elem);
                     save_fits("postFPMamp", fname);

                     WRITE_FULLFILENAME(fname, "%s/test_postFPMpha0_%02ld.fits", savedir, elem);
                     save_fits("postFPMpha", fname);
                 }*/

                // save diagnostics
                /* sprintf(command, "mv _DFT_foca %s/_DFT_foca_%02ld.fits", savedir, elem);
                 if(system(command) != 0)
                	printERROR(__FILE__,__func__,__LINE__, "system() returns non-zero value");

                 sprintf(command, "mv _DFT_focp %s/_DFT_focp_%02ld.fits", savedir, elem);
                 if(system(command) != 0)
                	printERROR(__FILE__,__func__,__LINE__, "system() returns non-zero value");
                	*/

            }
            long elemindex = optsyst[index].elemarrayindex[elem];

            list_image_ID();

            printf("optsyst[index].FOCMASKarray[elemindex].mode = %d\n", optsyst[index].FOCMASKarray[elemindex].mode);
            if(optsyst[index].FOCMASKarray[elemindex].mode == 1)
            {
                // we are computing using the 1 - FPM trick, so subtract the DFT result from the input light
                // this makes the larger computation more efficient
                arith_image_sub_inplace("_WFctmp", "_WFcout");
                mk_amph_from_complex("_WFctmp", imnameamp_out, imnamepha_out, 0);
            }
            else
            {
                mk_amph_from_complex("_WFcout", imnameamp_out, imnamepha_out, 0);
            }

            /* {
                 // TEST: OUTPUT
                 char fname[STRINGMAXLEN_FULLFILENAME];
                 WRITE_FULLFILENAME(fname, "%s/test_postFPMamp1_%02ld.fits", savedir, elem);
                 save_fits(imnameamp_out, fname);
                 WRITE_FULLFILENAME(fname, "%s/test_postFPMpha1_%02ld.fits", savedir, elem);
                 save_fits(imnamepha_out, fname);
             }*/

            delete_image_ID("_WFctmp", DELETE_IMAGE_ERRMODE_WARNING);
            delete_image_ID("_WFcout", DELETE_IMAGE_ERRMODE_WARNING);
            //  delete_image_ID("dftgrid", DELETE_IMAGE_ERRMODE_WARNING);
        }





        // computes the total flux at this point
        IDa = image_ID(imnameamp_out); // output of the current element
        optsyst[index].flux[elem] = 0.0;
        for(long kl = 0; kl < nblambda; kl++)
            for(uint64_t ii = 0; ii < size2; ii++)
            {
                optsyst[index].flux[elem] += data.image[IDa].array.F[kl * size2 + ii] *
                                             data.image[IDa].array.F[kl * size2 + ii];
            }

        printf("Element %ld  [%ld %ld]  Flux = %lf\n",
               elem,
               nblambda,
               (long) size2,
               optsyst[index].flux[elem] / nblambda);

        if(isnan(optsyst[index].flux[elem]) != 0)
        {
            exit(0);
        }

        if(optsyst[index].SAVE == 1)
        {
            printf("Saving intermediate plane [%ld] ... ", elem);

            char fname[STRINGMAXLEN_FULLFILENAME];

            WRITE_FULLFILENAME(fname, "./%s/WFamp%ld_%03ld.fits", savedir, index, elem);
            save_fits(imnameamp_out, fname);

            WRITE_FULLFILENAME(fname, "./%s/WFpha%ld_%03ld.fits", savedir, index, elem);
            save_fits(imnamepha_out, fname);

            printf("done\n");
            fflush(stdout);
        }
    }

    // we're at the last element so compute and save the final images psfi0 and psfc0
    //
    if((elem == optsyst[index].NBelem)
            && (optsyst[index].endmode == 0)) // Compute final focal plane image
    {
        char imnameamp[STRINGMAXLEN_IMGNAME];
        char imnamepha[STRINGMAXLEN_IMGNAME];
        char imnamere[STRINGMAXLEN_IMGNAME];
        char imnameim[STRINGMAXLEN_IMGNAME];


        printf("COMPUTING FINAL IMAGE AS FFT OF %ld\n", elem - 1);
        mk_complex_from_amph(imnameamp_out, imnamepha_out, "_WFctmp", 0);
        permut("_WFctmp"); // permute as needed for FFTW

        {
            char imname[STRINGMAXLEN_IMGNAME];
            WRITE_IMAGENAME(imname, "psfc%ld", index);
            // propagate to the final image plane via a direct 2d fft
            do2dfft("_WFctmp", imname);

            delete_image_ID("_WFctmp", DELETE_IMAGE_ERRMODE_WARNING);
            permut(imname);
            WRITE_IMAGENAME(imnameamp, "psfa%ld", index);
            WRITE_IMAGENAME(imnamepha, "psfp%ld", index);
            WRITE_IMAGENAME(imnamere, "psfre%ld", index);
            WRITE_IMAGENAME(imnameim, "psfim%ld", index);

            mk_reim_from_complex(imname, imnamere, imnameim, sharedmem);
            mk_amph_from_complex(imname, imnameamp, imnamepha, sharedmem);
        }


        if(optsyst[index].SAVE == 1)
        {
            char fname[STRINGMAXLEN_FULLFILENAME];

            // save PSF in both Re/Im and amp/phase representations if desired
            WRITE_FULLFILENAME(fname, "%s/psfa%ld.fits", savedir, index);
            save_fits(imnameamp, fname);

            WRITE_FULLFILENAME(fname, "%s/psfp%ld.fits", savedir, index);
            save_fits(imnamepha, fname);

            WRITE_FULLFILENAME(fname, "%s/psfre%ld.fits", savedir, index);
            save_fits(imnamere, fname);

            WRITE_FULLFILENAME(fname, "%s/psfim%ld.fits", savedir, index);
            save_fits(imnameim, fname);
        }


        {
            imageID ID = image_ID(imnameamp);
            imageID IDre = image_ID(imnamere);
            imageID IDim = image_ID(imnameim);
            // normalize so the intensity sums to 1
            for(uint64_t ii = 0; ii < size2 * nblambda; ii++)
            {
                data.image[ID].array.F[ii] /= sqrt(size2 * optsyst[index].flux[0] / nblambda);
                data.image[IDre].array.F[ii] /= sqrt(size2 * optsyst[index].flux[0] / nblambda);
                data.image[IDim].array.F[ii] /= sqrt(size2 * optsyst[index].flux[0] / nblambda);
            }
        }

        // compute and print total flux
        {
            char imname[STRINGMAXLEN_IMGNAME];

            WRITE_IMAGENAME(imname, "psfi%ld", index);
            arith_image_mult(imnameamp, imnameamp, imname); // intensity is amp^2

            double total = arith_image_total(imname) /
                           nblambda; // total flux "averaged" over wavelength
            printf("TOTAL = %lf\n", total);


            if(optsyst[index].SAVE == 1)
            {
                // save PSF intensity if desired

                char fname[STRINGMAXLEN_FULLFILENAME];

                WRITE_FULLFILENAME(fname, "%s/psfi%ld.fits", savedir, index);
                save_fits(imname, fname);
            }
        }
    }


    free(imsizearray);

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

