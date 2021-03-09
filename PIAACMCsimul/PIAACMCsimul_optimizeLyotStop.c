/**
 * @file    PIAACMCsimul_optimizeLyotStop.c
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
    // initial guess places Lyot stops regularly from zmin to zmax
    // light propagates from zmin to zmax
    // we start with a single mask in zmax, and work back
    //
    double ratio = 1.0; // output metric... not used yet, currently a placeholder

    long ID, IDa, IDp;
    long nblambda; // number of wavelengths, read from input cube
    float *zarray;
    long l;
    double zprop;

    char nameamp[500];
    char namepha[500];
    char nameint[500];
    char fname[500];
    char fname1[500];
    long xsize = 0;
    long ysize = 0;
    long ii, jj, k, m;

    float *rinarray;
    float *routarray;
    float dr = 0.02;
    double tot;


    double *totarray;
    double *tot2array;
    long IDzone;
    double x, y, r;

    double zbest, valbest, val;
    long lbest;

    long IDincohc, IDint, IDmc, IDmc1;
    float rsl;
    long iter;
    long NBiter = 100;
    double val1;
    long IDm;
    char name[500];

    long IDre, IDim, IDreg, IDimg;
    double amp, pha, re, im;
    float sigma, sigma1;
    int filter_size;

    long IDlscumul;
    long IDintg, IDintgg;

    float *rarray;


    FILE *fp;
    double alpha = 1.01; // norm alpha used to identify best plane

	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif



    sigma = 0.01*piaacmc[0].beamrad/piaacmc[0].pixscale;
    filter_size = (long) (sigma*2.0);

    zarray = (float*) malloc(sizeof(float)*NBz);

    rinarray = (float*) malloc(sizeof(float)*NBmasks);
    routarray = (float*) malloc(sizeof(float)*NBmasks);

    totarray = (double*) malloc(sizeof(double)*NBmasks*NBz);
    tot2array = (double*) malloc(sizeof(double)*NBmasks*NBz);


    routarray[0] = 1.0;
    rinarray[0] = 1.0 - 1.0/NBmasks;
    for(m=1; m<NBmasks; m++)
    {
        routarray[m] = rinarray[m-1];
        rinarray[m] = routarray[m] - (1.0-piaacmc[0].centObs1)/NBmasks;
    }
    rinarray[NBmasks-1] = 0.0;

    for(m=0; m<NBmasks; m++)
        printf("annulus %ld : %f - %f\n", m, routarray[m], rinarray[m]);


    IDa = image_ID(IDamp_name);
    IDp = image_ID(IDpha_name);
    xsize = data.image[IDa].md[0].size[0];
    ysize = data.image[IDa].md[0].size[1];

    
    rarray = (float*) malloc(sizeof(float)*xsize*ysize);
    IDincohc = image_ID(IDincohc_name);



    if(data.image[IDa].md[0].naxis==3)
        nblambda = data.image[IDa].md[0].size[2];
    else
        nblambda = 1;

    IDzone = create_2Dimage_ID("LMzonemap", xsize, ysize);
    for(ii=0; ii<xsize; ii++)
        for(jj=0; jj<ysize; jj++)
        {
            data.image[IDzone].array.F[jj*xsize+ii] = -2;
            x = (1.0*ii-0.5*xsize)/(piaacmc[0].beamrad/piaacmc[0].pixscale);
            y = (1.0*jj-0.5*xsize)/(piaacmc[0].beamrad/piaacmc[0].pixscale);
            r = sqrt(x*x+y*y);
            rarray[jj*xsize+ii] = r;
            for(m=0; m<NBmasks; m++)
                if((r>rinarray[m]-0.0001)&&(r<routarray[m]+0.0001))
                    data.image[IDzone].array.F[jj*xsize+ii] = m;
        }
    if( piaacmcsimul_var.PIAACMC_save==1)
    {
        sprintf(fname, "!%s/LMzonemap.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("LMzonemap", fname);
    }
    // initialize zarray
    for(l=0; l<NBz; l++)
        zarray[l] = zmin + (zmax-zmin)*l/(NBz-1);

    //  save_fits(nameint, fname);


    IDre = create_2Dimage_ID("retmpim", xsize, ysize);
    IDim = create_2Dimage_ID("imtmpim", xsize, ysize);

    ID = create_3Dimage_ID("LMintC", xsize, ysize, NBz);
    IDintg = create_2Dimage_ID("tmpintg", xsize, ysize);

    for(l=0; l<NBz; l++)
    {
        sprintf(nameamp, "LMPamp%02ld", l);
        sprintf(namepha, "LMPpha%02ld", l);
        zprop = zarray[l];
        OptSystProp_propagateCube(optsyst, 0, IDamp_name, IDpha_name, nameamp, namepha, zprop, 0);

     
        IDa = image_ID(nameamp);
        IDp = image_ID(namepha);
     
        for(k=0; k<nblambda; k++)
        {
        for(ii=0; ii<xsize*ysize; ii++)
            {
                amp = data.image[IDa].array.F[k*xsize*ysize+ii];
                pha = data.image[IDp].array.F[k*xsize*ysize+ii];
                data.image[IDre].array.F[ii] = amp*cos(pha);
                data.image[IDim].array.F[ii] = amp*sin(pha);
            }
            IDreg = gauss_filter("retmpim", "retmpimg", sigma, filter_size);
            IDimg = gauss_filter("imtmpim", "imtmpimg", sigma, filter_size);

            for(ii=0; ii<xsize*ysize; ii++)
            {
                re = data.image[IDreg].array.F[ii];
                im = data.image[IDimg].array.F[ii];
                data.image[IDintg].array.F[ii] = re*re + im*im;
            }
            IDintgg = gauss_filter("tmpintg", "tmpintgg", 2.0*sigma, filter_size);


            for(ii=0; ii<xsize*ysize; ii++)
                data.image[ID].array.F[l*xsize*ysize+ii] += data.image[IDintgg].array.F[ii];
                //data.image[IDa].array.F[ii]*data.image[IDa].array.F[ii]; //data.image[IDintgg].array.F[ii];

            delete_image_ID("retmpimg");
            delete_image_ID("imtmpimg");
            delete_image_ID("tmpintgg");
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
        delete_image_ID(nameamp);
        delete_image_ID(namepha);
    }
    delete_image_ID("retmpim");
   delete_image_ID("imtmpim");


    sprintf(fname,  "!%s/LMintC.fits", piaacmcsimul_var.piaacmcconfdir);
    save_fits("LMintC", fname);




    for(l=0; l<NBz; l++)
        for(m=0;m<NBmasks;m++)
            {
                totarray[l*NBmasks+m] = 0.0;
                tot2array[l*NBmasks+m] = 0.0;
            }

    for(l=0; l<NBz; l++)
    {
        for(ii=0; ii<xsize*ysize; ii++)
        {
            m = (long) (data.image[IDzone].array.F[ii]+0.1);
    
            if((m>-1)&&(m<NBmasks)&&(rarray[ii]<1.0)&&(rarray[ii]>0.9*piaacmc[0].centObs1))
            {
                totarray[l*NBmasks+m] += data.image[IDincohc].array.F[l*xsize*ysize+ii];
                tot2array[l*NBmasks+m] += pow(data.image[IDincohc].array.F[l*xsize*ysize+ii], alpha);
            }
        }
    }

    IDmc = create_2Dimage_ID("Lcomb", xsize, ysize);
    IDmc1 = create_2Dimage_ID("LcombOA", xsize, ysize);
    for(ii=0; ii<xsize*ysize; ii++)
        data.image[IDmc1].array.F[ii] = 0.0;

    sprintf(fname1, "%s/LyotMasks_zpos.txt", piaacmcsimul_var.piaacmcconfdir);
    fp = fopen(fname1, "w");
    IDint = image_ID("LMintC");
    for(m=0; m<NBmasks; m++)
    {
        valbest = 0.0;
        lbest = 0;
        zbest = 0.0;
        for(l=0; l<NBz; l++)
        {
            val =  tot2array[l*NBmasks+m]/pow(totarray[l*NBmasks+m], alpha);
            printf("MASK %ld   z(%ld)= %f  ->  %g   ( %g %g) \n", m, l, zarray[l], val, tot2array[l*NBmasks+m], totarray[l*NBmasks+m]);
            if(val>valbest)
            {
                valbest = val;
                zbest = zarray[l];
                lbest = l;
            }
        }
        printf(" ==========  MASK %ld   BEST CONJUGATION : %ld %f (%g)\n", m, lbest, zbest, valbest);
        piaacmc[0].LyotStop_zpos[m] = zbest; // relative to starting plane
        fprintf(fp, "%02ld %f\n", lbest, zbest);

        for(ii=0; ii<xsize*ysize; ii++)
            if(m==data.image[IDzone].array.F[ii])
                {
                    data.image[IDmc].array.F[ii] = data.image[IDint].array.F[lbest*xsize*ysize+ii];
                    data.image[IDmc1].array.F[ii] = data.image[IDincohc].array.F[lbest*xsize*ysize+ii];
                }
    }
    fclose(fp);

    if( piaacmcsimul_var.PIAACMC_save==1)
    {
        sprintf(fname, "!%s/Lcomb.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("Lcomb", fname);
        sprintf(fname, "!%s/LcombOA.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("LcombOA", fname);
    }
    
    
    /// call PIAACMCsimul_mkLyotMask() 
    ID = PIAACMCsimul_mkLyotMask("LcombOA", "Lcomb", "LMzonemap", throughput, "LMask");
    
    if( piaacmcsimul_var.PIAACMC_save==1)
    {
        sprintf(fname, "!%s/LMask.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("LMask", fname);
    }
    delete_image_ID("Lcomb");

    IDlscumul = create_2Dimage_ID("LMcumul", xsize, ysize);
    for(ii=0; ii<xsize*ysize; ii++)
        data.image[IDlscumul].array.F[ii] = 1.0;

    for(m=0; m<NBmasks; m++)
    {
        sprintf(name, "optLM%02ld", m);
        IDm = create_2Dimage_ID(name, xsize, ysize);
        for(ii=0; ii<xsize; ii++)
            for(jj=0; jj<ysize; jj++)
            {
                x = (1.0*ii-0.5*xsize)/(piaacmc[0].beamrad/piaacmc[0].pixscale);
                y = (1.0*jj-0.5*xsize)/(piaacmc[0].beamrad/piaacmc[0].pixscale);
                r = sqrt(x*x+y*y);

                if((r>rinarray[m]-dr)&&(r<routarray[m]+dr))
                    data.image[IDm].array.F[jj*xsize+ii] = data.image[ID].array.F[jj*xsize+ii];
                else
                    data.image[IDm].array.F[jj*xsize+ii] = 1.0;
            }
        if(m==0)
            for(ii=0; ii<xsize; ii++)
                for(jj=0; jj<ysize; jj++)
                {
                    x = (1.0*ii-0.5*xsize)/(piaacmc[0].beamrad/piaacmc[0].pixscale);
                    y = (1.0*jj-0.5*xsize)/(piaacmc[0].beamrad/piaacmc[0].pixscale);
                    r = sqrt(x*x+y*y);
                    if(r>1.0)
                        data.image[IDm].array.F[jj*xsize+ii] = 0.0;
                }


        for(ii=0; ii<xsize*ysize; ii++)
        {
            data.image[IDm].array.F[ii] *= data.image[IDlscumul].array.F[ii];
            data.image[IDlscumul].array.F[ii] = data.image[IDm].array.F[ii];
        }


        sprintf(fname, "!%s/optLM%02ld.fits", piaacmcsimul_var.piaacmcconfdir, m);
        save_fits(name, fname);
    }


    delete_image_ID("LMcumul");

    free(totarray);
    free(tot2array);
    free(rinarray);
    free(routarray);
    free(zarray);
    free(rarray);

    delete_image_ID("LMzonemap");

    return(ratio);
}



