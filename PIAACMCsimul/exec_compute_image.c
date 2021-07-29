/**
 * @file    PIAAACMCsimul_compute_image.c
 * @brief   PIAA-type coronagraph design, execute compute image
 *
 */



// System includes
#include <stdio.h>
#include <stdlib.h>




// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"

#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"

#include "PIAACMCsimul_computePSF.h"
#include "init_piaacmcopticaldesign.h"

#include "PIAAshape/makePIAAshapes.h"







/**
 *
 * @brief Compute PSF or image scene
 *
 */
errno_t exec_compute_image()
{
    DEBUG_TRACE_FSTART();

    double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;
    int PIAACMC_WFCmode = 0;
    double xpos, ypos, fval;
    imageID IDscene;
    imageID ID;
    imageID IDpsfi0;
    double valref;


    // Run existing config for on-axis point source. If new, create centrally obscured idealized PIAACMC
    // compatible with wavefront control
    printf("=================================== mode 000 ===================================\n");
    // Either load a set of point sources from "scene.txt" or use a single on-axis point source,
    // and create the image for these sources by computing and adding their PSFs

    {   // load some more cli variables
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

        piaacmcparams.PIAACMC_fpmtype = 0; // idealized (default)
        if((IDv = variable_ID("PIAACMC_fpmtype")) != -1)
        {
            piaacmcparams.PIAACMC_fpmtype = (int)(data.variable[IDv].value.f + 0.1);
        }
        printf("PIAACMC_fpmtype = %d\n", piaacmcparams.PIAACMC_fpmtype);

        PIAACMC_WFCmode = 0; // number of DMs
        if((IDv = variable_ID("PIAACMC_WFCmode")) != -1)
        {
            PIAACMC_WFCmode = (int)(data.variable[IDv].value.f + 0.1);
        }
        printf("PIAACMC_WFCmode = %d\n", PIAACMC_WFCmode);
    }



    // force creation of the FPM zone amplitudes by called functions
    piaacmcparams.FORCE_CREATE_fpmza = 1;

    // main initialization function to set up the piaacmc structure
    {
        uint64_t initflag = INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__READCONF;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF;

        if(piaacmcparams.PIAACMC_fpmtype == 1)
        {
            initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL;
        }

        if(PIAACMC_WFCmode == 1)
        {
            initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__WSCMODE;
        }
        FUNC_CHECK_RETURN(
            init_piaacmcopticaldesign(
                fpmradld,
                centobs0,
                centobs1,
                initflag,
                NULL
            )
        );
    }

    // make the mirror or lenses shapes
    FUNC_CHECK_RETURN(
        makePIAAshapes()
    );

    // use 1-fpm normalization for efficiency
    piaacmcopticalsystem.FOCMASKarray[0].mode = 1;



    // if file "LOOPMODE" exists, run PSF computation as a loop, waiting on OPDerrC to change
    FILE * fp = fopen("LOOPMODE.txt", "r");
    if(fp != NULL)
    {
        printf("RUNNING PSF LOOP COMPUTATION\n");

        uint32_t *sizearray;
        sizearray = (uint32_t *) malloc(sizeof(uint32_t) * 2);
        if(sizearray == NULL)
        {
            FUNC_RETURN_FAILURE("malloc error");
        }
        sizearray[0] = piaacmcopticaldesign.size;
        sizearray[1] = piaacmcopticaldesign.size;

        imageID IDopderrC;
        FUNC_CHECK_RETURN(
            create_image_ID("opderr", 2, sizearray, _DATATYPE_FLOAT, 1, 0, 0, &IDopderrC)
        );
        COREMOD_MEMORY_image_set_createsem("opderr", 10);
        free(sizearray);

        long sizecrop = piaacmcopticaldesign.size / 16;
        sizearray = (uint32_t *) malloc(sizeof(uint32_t) * 3);
        if(sizearray == NULL)
        {
            FUNC_RETURN_FAILURE("malloc error");
        }
        sizearray[0] = sizecrop;
        sizearray[1] = sizecrop;
        sizearray[2] = piaacmcopticaldesign.nblambda;
        FUNC_CHECK_RETURN(
            create_image_ID("psfiout0", 3, sizearray, _DATATYPE_FLOAT, 1, 0, 0, &IDpsfi0)
        );
        free(sizearray);

        long iter = 0;
        while(iter < 10)
        {
            {
                double cval = 0.0;
                FUNC_CHECK_RETURN(
                    PIAACMCsimul_computePSF(xpos, ypos, 0, piaacmcopticalsystem.NBelem, 1, 0, 0, 1, &cval)
                );
            }

            ID = image_ID("psfi0");

            // copy results to IDpsfi0
            data.image[IDpsfi0].md[0].write = 1;

            for(long k = 0; k < piaacmcopticaldesign.nblambda; k++)
                for(long ii1 = 0; ii1 < sizecrop; ii1++)
                    for(long jj1 = 0; jj1 < sizecrop; jj1++)
                    {
                        long ii = ii1 + (piaacmcopticaldesign.size - sizecrop) / 2;
                        long jj = jj1 + (piaacmcopticaldesign.size - sizecrop) / 2;
                        data.image[IDpsfi0].array.F[k * sizecrop * sizecrop + jj1 * sizecrop + ii1] =
                            data.image[ID].array.F[k * piaacmcopticaldesign.size * piaacmcopticaldesign.size + jj *
                                                   piaacmcopticaldesign.size + ii];
                    }
            COREMOD_MEMORY_image_set_sempost_byID(IDpsfi0, -1);
            data.image[IDpsfi0].md[0].cnt0 ++;
            data.image[IDpsfi0].md[0].write = 0;

            COREMOD_MEMORY_image_set_semwait("opderr", 0);
            // drive semaphore #1 to zero
            while(sem_trywait(data.image[IDopderrC].semptr[0]) == 0) {}
            //iter++;
        }
    }
    else
    {
        // if file "scene.txt" exists, compute series of PSFs and sum
        FILE * fpscene = fopen("SCENE.txt", "r");
        if(fpscene != NULL)
        {
            int initscene = 0;
            // for each source in the scene, read position and flux
            while(fscanf(fpscene, "%lf %lf %lf\n", &xpos, &ypos, &fval) == 3)
            {
                printf("COMPUTING PSF AT POSITION %lf %lf, flux  = %g\n", xpos, ypos, fval);
                // make the actual PSF
                FUNC_CHECK_RETURN(
                    PIAACMCsimul_computePSF(
                        xpos, ypos, 0, piaacmcopticalsystem.NBelem, 1, 0, 0, 1, NULL)
                );


                // get the image "psfi0" index, which was created in PIAACMCsimul_computePSF
                ID = image_ID("psfi0");
                // get image size.  3rd dimension is wavelength
                uint32_t xsize = data.image[ID].md[0].size[0];
                uint32_t ysize = data.image[ID].md[0].size[1];
                uint32_t zsize = data.image[ID].md[0].size[2];

                if(initscene == 0)
                {
                    initscene = 1;
                    // create 3D image to sum the PSFs into
                    create_3Dimage_ID("scene", xsize, ysize, zsize, &IDscene);
                }
                ID = image_ID("psfi0");
                // sum the current PSF into the image: summed image is IDscene, source is ID
                for(uint64_t ii = 0; ii < xsize * ysize * zsize; ii++)
                {
                    data.image[IDscene].array.F[ii] += fval * data.image[ID].array.F[ii];
                }


            }
            fclose(fpscene);
            // we're done!  Save it, overwriting previous scene.fits file
            save_fits("scene", "scene.fits");
        }
        else // scene.txt does not exist, just do an on-axis source
        {
            {
                errno_t fret =
                    PIAACMCsimul_computePSF(
                        0.0,
                        0.0,
                        0,
                        piaacmcopticalsystem.NBelem,
                        1,
                        0,
                        0,
                        1,
                        &valref
                    );

                if( fret != RETURN_SUCCESS)
                {
                    FUNC_RETURN_FAILURE("Call to PIAACMCsimul_computePSF failed");
                }
            }
            printf("valref = %g\n", valref);
        }
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}


