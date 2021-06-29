/**
 * @file    PIAACMCsimul_optimizeLyotStop.c
 * @brief   PIAA-type coronagraph design
 *
 * Can design both APLCMC and PIAACMC coronagraphs
 *
 *
 */



/* =============================================================================================== */
/* =============================================================================================== */
/*                                        HEADER FILES                                             */
/* =============================================================================================== */
/* =============================================================================================== */

// System includes
#include <stdlib.h>
#include <stdio.h>
#include <math.h>



// milk includes
#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "image_filter/image_filter.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"




extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTSYST *optsyst;

extern OPTPIAACMCDESIGN *piaacmc;





/* =============================================================================================== */
/* =============================================================================================== */
/*                                    FUNCTION(S) SOURCE CODE                                      */
/* =============================================================================================== */
/* =============================================================================================== */

/**
 * @brief Lyot stops positions from zmin to zmax relative to current, working back (light goes from 0 to zmax)
 *
 * @param[in]	IDamp_name    image : 2D amplitude
 * @param[in]	IDpha_name    image : 2D phase
 * @param[in]   IDincohc_name image : 3D incoherent intensity
 * @param[in]   zmin          float : minimum propagation value
 * @param[in]   zmax          float : maximum propagation value
 * @param[in]   throughput    double: Geometric throughput of Lyot stop(s)
 * @param[in]   NBz           long  : Number of discrete propagation planes between zmin and zmax
 * @param[in]   NBmasks       long  : Number of Lyot stop(s)
 *
*/
double PIAACMCsimul_optimizeLyotStop(
    const char *IDamp_name,
    const char *IDpha_name,
    const char *IDincohc_name,
    float zmin,
    float zmax,
    double throughput,
    long NBz,
    long NBmasks
)
{
    DEBUG_TRACE_FSTART();

    // initial guess places Lyot stops regularly from zmin to zmax
    // light propagates from zmin to zmax
    // we start with a single mask in zmax, and work back
    //
    double ratio = 1.0; // output metric... not used yet, currently a placeholder

    long nblambda; // number of wavelengths, read from input cube
    float *zarray;

    float *rinarray;
    float *routarray;
    float dr = 0.02;


    double *totarray;
    double *tot2array;

    imageID IDincohc, IDint, IDmc, IDmc1;

    float *rarray;

    double alpha = 1.01; // norm alpha used to identify best plane


#ifdef PIAASIMUL_LOGFUNC0
    PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__,
                                 "");
#endif



    zarray = (float *) malloc(sizeof(float) * NBz);
    if(zarray == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }

    rinarray = (float *) malloc(sizeof(float) * NBmasks);
    if(rinarray == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }

    routarray = (float *) malloc(sizeof(float) * NBmasks);
    if(routarray == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }

    totarray = (double *) malloc(sizeof(double) * NBmasks * NBz);
    if(totarray == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }

    tot2array = (double *) malloc(sizeof(double) * NBmasks * NBz);
    if(tot2array == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }

    routarray[0] = 1.0;
    rinarray[0] = 1.0 - 1.0 / NBmasks;
    for(long m = 1; m < NBmasks; m++)
    {
        routarray[m] = rinarray[m - 1];
        rinarray[m] = routarray[m] - (1.0 - piaacmc[0].centObs1) / NBmasks;
    }
    rinarray[NBmasks - 1] = 0.0;

    for(long m = 0; m < NBmasks; m++)
    {
        printf("annulus %ld : %f - %f\n", m, routarray[m], rinarray[m]);
    }


    imageID IDa = image_ID(IDamp_name);
//    imageID IDp = image_ID(IDpha_name);
    uint32_t xsize = data.image[IDa].md[0].size[0];
    uint32_t ysize = data.image[IDa].md[0].size[1];
    uint64_t xysize = xsize;
    xysize *= ysize;


    rarray = (float *) malloc(sizeof(float) * xsize * ysize);
    if(rarray == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }

    IDincohc = image_ID(IDincohc_name);



    if(data.image[IDa].md[0].naxis == 3)
    {
        nblambda = data.image[IDa].md[0].size[2];
    }
    else
    {
        nblambda = 1;
    }

    imageID IDzone;
    create_2Dimage_ID("LMzonemap", xsize, ysize, &IDzone);

    for(uint32_t ii = 0; ii < xsize; ii++)
        for(uint32_t jj = 0; jj < ysize; jj++)
        {
            data.image[IDzone].array.F[jj * xsize + ii] = -2;
            double x = (1.0 * ii - 0.5 * xsize) / (piaacmc[0].beamrad / piaacmc[0].pixscale);
            double y = (1.0 * jj - 0.5 * xsize) / (piaacmc[0].beamrad / piaacmc[0].pixscale);
            double r = sqrt(x * x + y * y);
            rarray[jj * xsize + ii] = r;
            for(long m = 0; m < NBmasks; m++)
                if((r > rinarray[m] - 0.0001) && (r < routarray[m] + 0.0001))
                {
                    data.image[IDzone].array.F[jj * xsize + ii] = m;
                }
        }
    if(piaacmcsimul_var.PIAACMC_save == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/LMzonemap.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("LMzonemap", fname);
    }
    // initialize zarray
    for(long l = 0; l < NBz; l++)
    {
        zarray[l] = zmin + (zmax - zmin) * l / (NBz - 1);
    }

    //  save_fits(nameint, fname);


    imageID IDre;
    create_2Dimage_ID("retmpim", xsize, ysize, &IDre);

    imageID IDim;
    create_2Dimage_ID("imtmpim", xsize, ysize, &IDim);

    imageID ID_LMintC;
    create_3Dimage_ID("LMintC", xsize, ysize, NBz, &ID_LMintC);

    imageID IDintg;
    create_2Dimage_ID("tmpintg", xsize, ysize, &IDintg);

    float sigma = 0.01 * piaacmc[0].beamrad / piaacmc[0].pixscale;
    int filter_size = (long)(sigma * 2.0);

    for(long l = 0; l < NBz; l++)
    {
        char nameamp[STRINGMAXLEN_IMGNAME];
        WRITE_IMAGENAME(nameamp, "LMPamp%02ld", l);

        char namepha[STRINGMAXLEN_IMGNAME];
        WRITE_IMAGENAME(namepha, "LMPpha%02ld", l);

        double zprop = zarray[l];
        OptSystProp_propagateCube(optsyst, 0, IDamp_name, IDpha_name, nameamp, namepha,
                                  zprop, 0);


        imageID IDa1 = image_ID(nameamp);
        imageID IDp1 = image_ID(namepha);

        for(long k = 0; k < nblambda; k++)
        {
            for(uint64_t ii = 0; ii < xysize; ii++)
            {
                double amp = data.image[IDa1].array.F[k * xsize * ysize + ii];
                double pha = data.image[IDp1].array.F[k * xsize * ysize + ii];
                data.image[IDre].array.F[ii] = amp * cos(pha);
                data.image[IDim].array.F[ii] = amp * sin(pha);
            }
            imageID IDreg = gauss_filter("retmpim", "retmpimg", sigma, filter_size);
            imageID IDimg = gauss_filter("imtmpim", "imtmpimg", sigma, filter_size);

            for(uint64_t ii = 0; ii < xysize; ii++)
            {
                double re = data.image[IDreg].array.F[ii];
                double im = data.image[IDimg].array.F[ii];
                data.image[IDintg].array.F[ii] = re * re + im * im;
            }
            imageID IDintgg = gauss_filter("tmpintg", "tmpintgg", 2.0 * sigma, filter_size);


            for(uint64_t ii = 0; ii < xysize; ii++)
            {
                data.image[ID_LMintC].array.F[l * xsize * ysize + ii] +=
                    data.image[IDintgg].array.F[ii];
            }
            //data.image[IDa1].array.F[ii]*data.image[IDa1].array.F[ii]; //data.image[IDintgg].array.F[ii];

            delete_image_ID("retmpimg", DELETE_IMAGE_ERRMODE_WARNING);
            delete_image_ID("imtmpimg", DELETE_IMAGE_ERRMODE_WARNING);
            delete_image_ID("tmpintgg", DELETE_IMAGE_ERRMODE_WARNING);
        }

        /*    for(ii=0; ii<xsize*ysize; ii++)
            {
                m = (long) (data.image[IDzone].array.F[ii]+0.1);


                if((m>-1)&&(m<NBmasks)&&(rarray[ii]<1.0)&&(rarray[ii]>0.9*piaacmc[0].centObs1))
                {
                    totarray[l*NBmasks+m] += data.image[ID].array.F[l*xsize*ysize+ii];
                    tot2array[l*NBmasks+m] += pow(data.image[ID].array.F[l*xsize*ysize+ii], alpha);
                }
            }*/
        delete_image_ID(nameamp, DELETE_IMAGE_ERRMODE_WARNING);
        delete_image_ID(namepha, DELETE_IMAGE_ERRMODE_WARNING);
    }
    delete_image_ID("retmpim", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("imtmpim", DELETE_IMAGE_ERRMODE_WARNING);

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname,  "%s/LMintC.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("LMintC", fname);
    }



    for(long l = 0; l < NBz; l++)
        for(long m = 0; m < NBmasks; m++)
        {
            totarray[l * NBmasks + m] = 0.0;
            tot2array[l * NBmasks + m] = 0.0;
        }

    for(long l = 0; l < NBz; l++)
    {
        for(uint64_t ii = 0; ii < xysize; ii++)
        {
            long m = (long) (data.image[IDzone].array.F[ii] + 0.1);

            if((m > -1) && (m < NBmasks) && (rarray[ii] < 1.0)
                    && (rarray[ii] > 0.9 * piaacmc[0].centObs1))
            {
                totarray[l * NBmasks + m] += data.image[IDincohc].array.F[l * xsize * ysize +
                                             ii];
                tot2array[l * NBmasks + m] += pow(data.image[IDincohc].array.F[l * xsize * ysize
                                                  + ii], alpha);
            }
        }
    }

    create_2Dimage_ID("Lcomb", xsize, ysize, &IDmc);
    create_2Dimage_ID("LcombOA", xsize, ysize, &IDmc1);
    for(uint64_t ii = 0; ii < xysize; ii++)
    {
        data.image[IDmc1].array.F[ii] = 0.0;
    }

    {
        char fname1[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname1, "%s/LyotMasks_zpos.txt", piaacmcsimul_var.piaacmcconfdir);
        FILE * fp = fopen(fname1, "w");
        IDint = image_ID("LMintC");
        for(long m = 0; m < NBmasks; m++)
        {
            double valbest = 0.0;
            long lbest = 0;
            double zbest = 0.0;
            for(long l = 0; l < NBz; l++)
            {
                double val =  tot2array[l * NBmasks + m] / pow(totarray[l * NBmasks + m], alpha);
                printf("MASK %ld   z(%ld)= %f  ->  %g   ( %g %g) \n", m, l, zarray[l], val,
                       tot2array[l * NBmasks + m], totarray[l * NBmasks + m]);
                if(val > valbest)
                {
                    valbest = val;
                    zbest = zarray[l];
                    lbest = l;
                }
            }
            printf(" ==========  MASK %ld   BEST CONJUGATION : %ld %f (%g)\n", m, lbest,
                   zbest, valbest);
            piaacmc[0].LyotStop_zpos[m] = zbest; // relative to starting plane
            fprintf(fp, "%02ld %f\n", lbest, zbest);

            for(uint64_t ii = 0; ii < xysize; ii++)
                if(m == data.image[IDzone].array.F[ii])
                {
                    data.image[IDmc].array.F[ii] = data.image[IDint].array.F[lbest * xsize * ysize +
                                                   ii];
                    data.image[IDmc1].array.F[ii] = data.image[IDincohc].array.F[lbest * xsize *
                                                    ysize + ii];
                }
        }
        fclose(fp);
    }

    if(piaacmcsimul_var.PIAACMC_save == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];

        WRITE_FULLFILENAME(fname, "%s/Lcomb.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("Lcomb", fname);

        WRITE_FULLFILENAME(fname, "%s/LcombOA.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("LcombOA", fname);
    }


    /// call PIAACMCsimul_mkLyotMask()
    imageID IDlyotmask = PIAACMCsimul_mkLyotMask("LcombOA", "Lcomb", "LMzonemap", throughput,
                         "LMask");

    if(piaacmcsimul_var.PIAACMC_save == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/LMask.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("LMask", fname);
    }

    delete_image_ID("Lcomb", DELETE_IMAGE_ERRMODE_WARNING);

    imageID IDlscumul;
    create_2Dimage_ID("LMcumul", xsize, ysize, &IDlscumul);

    for(uint64_t ii = 0; ii < xysize; ii++)
    {
        data.image[IDlscumul].array.F[ii] = 1.0;
    }

    for(long m = 0; m < NBmasks; m++)
    {
        char name[STRINGMAXLEN_IMGNAME];
        WRITE_IMAGENAME(name, "optLM%02ld", m);

        imageID IDm;
        create_2Dimage_ID(name, xsize, ysize, &IDm);

        for(uint32_t ii = 0; ii < xsize; ii++)
            for(uint32_t jj = 0; jj < ysize; jj++)
            {
                double x = (1.0 * ii - 0.5 * xsize) / (piaacmc[0].beamrad / piaacmc[0].pixscale);
                double y = (1.0 * jj - 0.5 * xsize) / (piaacmc[0].beamrad / piaacmc[0].pixscale);
                double r = sqrt(x * x + y * y);

                if((r > rinarray[m] - dr) && (r < routarray[m] + dr))
                {
                    data.image[IDm].array.F[jj * xsize + ii] = data.image[IDlyotmask].array.F[jj * xsize +
                            ii];
                }
                else
                {
                    data.image[IDm].array.F[jj * xsize + ii] = 1.0;
                }
            }
        if(m == 0)
            for(uint32_t ii = 0; ii < xsize; ii++)
                for(uint32_t jj = 0; jj < ysize; jj++)
                {
                    double x = (1.0 * ii - 0.5 * xsize) / (piaacmc[0].beamrad / piaacmc[0].pixscale);
                    double y = (1.0 * jj - 0.5 * xsize) / (piaacmc[0].beamrad / piaacmc[0].pixscale);
                    double r = sqrt(x * x + y * y);
                    if(r > 1.0)
                    {
                        data.image[IDm].array.F[jj * xsize + ii] = 0.0;
                    }
                }


        for(uint64_t ii = 0; ii < xysize; ii++)
        {
            data.image[IDm].array.F[ii] *= data.image[IDlscumul].array.F[ii];
            data.image[IDlscumul].array.F[ii] = data.image[IDm].array.F[ii];
        }


        {
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(fname, "%s/optLM%02ld.fits", piaacmcsimul_var.piaacmcconfdir, m);
            save_fits(name, fname);
        }
    }


    delete_image_ID("LMcumul", DELETE_IMAGE_ERRMODE_WARNING);

    free(totarray);
    free(tot2array);
    free(rinarray);
    free(routarray);
    free(zarray);
    free(rarray);

    delete_image_ID("LMzonemap", DELETE_IMAGE_ERRMODE_WARNING);

    DEBUG_TRACE_FEXIT();
    return(ratio);
}



