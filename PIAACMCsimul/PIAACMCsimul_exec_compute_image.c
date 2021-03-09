/**
 * @file    PIAAACMCsimul_compute_image.c
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
 *
 * @brief Compute PSF or image scene
 *
 */

int PIAACMCsimul_exec_compute_image()
{
    long IDv;
    double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;
    int PIAACMC_WFCmode = 0;
    uint32_t *sizearray;
    FILE *fp;
    long IDopderrC;
    long sizecrop;
    long k, ii, jj, ii1, jj1;
    double xpos, ypos, fval;
    int initscene;
    long IDscene;
    long xsize, ysize, zsize;
    long ID;
    long iter;
    long IDpsfi0;
    double valref;


    // Run existing config for on-axis point source. If new, create centrally obscured idealized PIAACMC
    // compatible with wavefront control
    printf("=================================== mode 000 ===================================\n");
    // Either load a set of point sources from "scene.txt" or use a single on-axis point source,
    // and create the image for these sources by computing and adding their PSFs

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

    piaacmcsimul_var.PIAACMC_fpmtype = 0; // idealized (default)
    if((IDv=variable_ID("PIAACMC_fpmtype"))!=-1)
        piaacmcsimul_var.PIAACMC_fpmtype = (int) (data.variable[IDv].value.f + 0.1);
    printf("PIAACMC_fpmtype = %d\n", piaacmcsimul_var.PIAACMC_fpmtype);

    PIAACMC_WFCmode = 0; // number of DMs
    if((IDv=variable_ID("PIAACMC_WFCmode"))!=-1)
        PIAACMC_WFCmode = (int) (data.variable[IDv].value.f + 0.1);
    printf("PIAACMC_WFCmode = %d\n", PIAACMC_WFCmode);




    // force creation of the FPM zone amplitudes by called functions
    piaacmcsimul_var.FORCE_CREATE_fpmza = 1;

    // main initialization function to set up the piaacmc structure
    PIAACMCsimul_initpiaacmcconf(piaacmcsimul_var.PIAACMC_fpmtype, fpmradld, centobs0, centobs1, PIAACMC_WFCmode, 1);
    // make the mirror or lenses shapes 
    PIAACMCsimul_makePIAAshapes(piaacmc, 0);
    optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm normalization for efficiency



    // if file "LOOPMODE" exists, run PSF computation as a loop, waiting on OPDerrC to change
    fp = fopen("LOOPMODE.txt", "r");
    if(fp != NULL)
    {
        printf("RUNNING PSF LOOP COMPUTATION\n");
        sizearray = (uint32_t*) malloc(sizeof(uint32_t)*2);
        sizearray[0] = piaacmc[0].size;
        sizearray[1] = piaacmc[0].size;

        IDopderrC = create_image_ID("opderr", 2, sizearray, _DATATYPE_FLOAT, 1, 0);
        COREMOD_MEMORY_image_set_createsem("opderr", 10);
        free(sizearray);

        sizecrop = piaacmc[0].size/16;
        sizearray = (uint32_t*) malloc(sizeof(uint32_t)*3);
        sizearray[0] = sizecrop;
        sizearray[1] = sizecrop;
        sizearray[2] = piaacmc[0].nblambda;
        IDpsfi0 = create_image_ID("psfiout0", 3, sizearray, _DATATYPE_FLOAT, 1, 0);
        free(sizearray);

        iter = 0;
        while(iter<10)
        {
            PIAACMCsimul_computePSF(xpos, ypos, 0, optsyst[0].NBelem, 1, 0, 0, 1);
            ID = image_ID("psfi0");

            // copy results to IDpsfi0
            data.image[IDpsfi0].md[0].write = 1;

            for(k=0; k<piaacmc[0].nblambda; k++)
                for(ii1=0; ii1<sizecrop; ii1++)
                    for(jj1=0; jj1<sizecrop; jj1++)
                    {
                        ii = ii1 + (piaacmc[0].size - sizecrop)/2;
                        jj = jj1 + (piaacmc[0].size - sizecrop)/2;
                        data.image[IDpsfi0].array.F[k*sizecrop*sizecrop + jj1*sizecrop + ii1] = data.image[ID].array.F[k*piaacmc[0].size*piaacmc[0].size + jj*piaacmc[0].size + ii];
                    }
            COREMOD_MEMORY_image_set_sempost_byID(IDpsfi0, -1);
            data.image[IDpsfi0].md[0].cnt0 ++;
            data.image[IDpsfi0].md[0].write = 0;

            COREMOD_MEMORY_image_set_semwait("opderr", 0);
            // drive semaphore #1 to zero
            while(sem_trywait(data.image[IDopderrC].semptr[0])==0) {}
            //iter++;
        }
    }
    else
    {
        // if file "scene.txt" exists, compute series of PSFs and sum
        fp = fopen("SCENE.txt", "r");
        if(fp!=NULL)
        {
            initscene = 0;
            // for each source in the scene, read position and flux
            while(fscanf(fp, "%lf %lf %lf\n", &xpos, &ypos, &fval) == 3)
            {
                printf("COMPUTING PSF AT POSITION %lf %lf, flux  = %g\n", xpos, ypos, fval);
                // make the actual PSF
                PIAACMCsimul_computePSF(xpos, ypos, 0, optsyst[0].NBelem, 1, 0, 0, 1);
                // get the image "psfi0" index, which was created in PIAACMCsimul_computePSF
                ID = image_ID("psfi0");
                // get image size.  3rd dimension is wavelength
                xsize = data.image[ID].md[0].size[0];
                ysize = data.image[ID].md[0].size[1];
                zsize = data.image[ID].md[0].size[2];

                if(initscene==0)
                {
                    initscene = 1;
                    // create 3D image to sum the PSFs into
                    IDscene = create_3Dimage_ID("scene", xsize, ysize, zsize);
                }
                ID = image_ID("psfi0");
                // sum the current PSF into the image: summed image is IDscene, source is ID
                for(ii=0; ii<xsize*ysize*zsize; ii++)
                    data.image[IDscene].array.F[ii] += fval*data.image[ID].array.F[ii];


            }
            fclose(fp);
            // we're done!  Save it, overwriting previous scene.fits file
            save_fits("scene", "!scene.fits");
        }
        else // scene.txt does not exist, just do an on-axis source
        {
            valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 1, 0, 0, 1);
            printf("valref = %g\n", valref);
        }
    }

    return 0;

}


