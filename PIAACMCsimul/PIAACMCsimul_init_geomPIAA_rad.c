/**
 * @file    PIAACMCsimul_init_geomPIAA_rad.c
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

#include "PIAACMCsimul/PIAACMCsimul.h"


extern OPTPIAACMCDESIGN *piaacmc;
extern PIAACMCsimul_varType piaacmcsimul_var;





/**
 * computes radial PIAA optics sag
 *
 * this function only works for circular PIAA
 * uses radial PIAACMC design to initialize PIAA optics shapes and focal plane mask
 */
int PIAACMCsimul_init_geomPIAA_rad(
		const char *IDapofit_name
		)
{
    long i, ii, k;
    double *pup0;
    double *pup1;
    double *flux0cumul;
    double *flux1cumul;

    long IDcoeff;
    long nbcoeff;
    double r;
    FILE *fp;
    double total;

    // to convert r ro r1 (assymptotic outer radius on pup1)
    double coeffa = 3.0; // convergence rate from linear to assymptotic value
    double coeffa1 = 0.5; // convergence radius limit (added to 1.0)

    double r0, r1;
    double FLUX0_in = 0.0; // inside central obstruction
    double FLUX0_out = 0.0; // outside beam edge
    double FLUX1_in = 0.0; // inside central obstruction
    double FLUX1_out = 0.0; // outside beam edge
    double normcoeff;

    // inner profile adjustment
    double a, b, verr, bstep, x, eps1;
    double t0, t0cnt, value;
    int dir, odir;
    long NBistep;
    long iioffset;

    double fluxdens, F0, F1, F2;
    double dr0, dr1, ndr0, ndr1;

    long NBpoints;
    double *piaar00;
    double *piaar11;
    double *piaar01;
    double *piaar10;
    long cnt;
    double tmp;
    double epsilon = 0.000000000001;

    double *piaaM0z;
    double *piaaM1z;
    double r0c, r1c, dx, dy, dist, y3, r0n, slope, dz;

    char fname[500];

	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif


    pup0 = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
    pup1 = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
    flux0cumul = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
    flux1cumul = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);


    // STEP 1: CREATE AMPLITUDE AND CUMULATIVE INTENSITY PROFILES


    // CREATE OUTPUT AMPLITUDE APODIZATION PROFILE AND ITS CUMUL

    IDcoeff = image_ID(IDapofit_name);
    nbcoeff = data.image[IDcoeff].md[0].size[0];
    printf("%ld coefficients\n", nbcoeff);

    total = 0.0;
    for(ii=0; ii<piaacmc[0].NBradpts; ii++)
    {
        pup1[ii] = 0.0;
        r = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r1lim;
        if(r<1.0)
            r1 = r;
        else
            r1 = 1.0 + (r-1) / pow((1.0 + pow(1.0/coeffa1 * (r-1),coeffa)), 1.0/coeffa);

        for(k=0; k<nbcoeff; k++)
            pup1[ii] += data.image[IDcoeff].array.F[k]*cos(r1*k*M_PI/ApoFitCosFact);
        if(r<piaacmc[0].centObs1)
            FLUX1_in += pup1[ii]*pup1[ii]*r;
        if(r>1.0)
            FLUX1_out += pup1[ii]*pup1[ii]*r;
        total += pup1[ii]*pup1[ii]*r;
        flux1cumul[ii] = total;
    }


    normcoeff = 1.0/(total-FLUX1_in-FLUX1_out);

    FLUX1_in *= normcoeff;
    FLUX1_out *= normcoeff;
    for(ii=0; ii<piaacmc[0].NBradpts; ii++)
    {
        r = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r1lim;
        flux1cumul[ii] *= normcoeff;
    }


    printf("outer fluxes 1: %lf %lf\n", FLUX1_in, FLUX1_out);




    // CREATE FLUX0

    total = 0.0;
    for(ii=0; ii<piaacmc[0].NBradpts; ii++)
    {
        r = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r0lim;
        pup0[ii] = 1.0;

        if(r<piaacmc[0].centObs0)
            FLUX0_in += pup0[ii]*pup0[ii]*r;
        if(r>1.0)
            FLUX0_out += pup0[ii]*pup0[ii]*r;

        total += pup0[ii]*pup0[ii]*r;
        flux0cumul[ii] = total;
    }
    normcoeff = 1.0/(total-FLUX0_in-FLUX0_out);

    FLUX0_in *= normcoeff;
    FLUX0_out *= normcoeff;

    printf("outer fluxes 0: %lf (%lf)    %lf\n", FLUX0_in, FLUX1_in, FLUX0_out);


    //
    // Compute inner pseudo profile
    //
    b = 0.5;
    bstep = 0.1;
    verr = 1.0;
    NBistep = piaacmc[0].centObs0*piaacmc[0].NBradpts/piaacmc[0].r0lim;
    //  innerprof_cumul = (double*) malloc(sizeof(double)*NBistep);
	dir = 1.0; // initial direction
    while(fabs(verr)>1.0e-9)
    {
        t0 = 0.0;
        t0cnt = 0.0;
        for(ii=0; ii<NBistep; ii++)
        {
            x = 1.0*ii/NBistep;
            r = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r0lim;
            a = 0.5;
            eps1 = 1e-8;
            if(x<eps1)
                x = eps1;
            if(x>1.0-eps1)
                x = 1.0-eps1;
            pup0[ii] =  b + (1.0-b)*(0.5+atan(-a/b/(x*x) + a/pow(x-1.0,2))/M_PI);
            t0 += r*pup0[ii]* pup0[ii];
            flux0cumul[ii] = t0;
        }

        verr = t0*normcoeff - FLUX1_in;

        odir = dir;
        if(verr>0.0) // too much light
        {
            b /= (1.0+bstep);
            dir = -1.0;
        }
        else
        {
            b *= (1.0+bstep);
            dir = 1.0;
        }
        if(odir*dir<0.0)
            bstep *= 0.1;
        printf(".");
        fflush(stdout);
    }
    printf("\n");
    printf("TOTAL = %f -> %g (%g %g)\n", b, t0*normcoeff, bstep, verr);


    // outer region
    b = 0.5;
    bstep = 0.1;
    verr = 1.0;
    NBistep = piaacmc[0].NBradpts*(piaacmc[0].r0lim-1.0)/piaacmc[0].r0lim;
    //  innerprof_cumul = (double*) malloc(sizeof(double)*NBistep);
    iioffset = (long) (1.0*piaacmc[0].NBradpts/piaacmc[0].r0lim);
    NBistep = piaacmc[0].NBradpts-iioffset;
    while(fabs(verr)>1.0e-9)
    {
        t0 = 0.0;
        t0cnt = 0.0;
        for(ii=0; ii<NBistep; ii++)
        {
            x = 1.0-1.0*ii/NBistep;
            r = 1.0+1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r0lim;
            a = 0.5;
            eps1 = 1e-8;
            if(x<eps1)
                x = eps1;
            if(x>1.0-eps1)
                x = 1.0-eps1;
            pup0[ii+iioffset] =  b + (1.0-b)*(0.5+atan(-a/b/(x*x) + a/pow(x-1.0,2))/M_PI);
            t0 += r*pup0[ii+iioffset]* pup0[ii+iioffset];
            flux0cumul[ii+iioffset] = t0;
        }

        verr = t0*normcoeff - FLUX1_out;

        odir = dir;
        if(verr>0.0) // too much light
        {
            b /= (1.0+bstep);
            dir = -1.0;
        }
        else
        {
            b *= (1.0+bstep);
            dir = 1.0;
        }
        if(odir*dir<0.0)
            bstep *= 0.1;
        printf(".");
        fflush(stdout);
    }
    printf("\n");
    printf("TOTAL = %f -> %g (%g %g)\n", b, t0*normcoeff, bstep, verr);



    total = 0.0;
    FLUX0_in = 0.0;
    FLUX0_out = 0.0;
    for(ii=0; ii<piaacmc[0].NBradpts; ii++)
    {
        r = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r0lim;
        if(r<piaacmc[0].centObs0)
            FLUX0_in += pup0[ii]*pup0[ii]*r;
        if(r>1.0)
            FLUX0_out += pup0[ii]*pup0[ii]*r;

        total += pup0[ii]*pup0[ii]*r;
        flux0cumul[ii] = total;
        flux0cumul[ii] *= normcoeff;;
    }
    FLUX0_in *= normcoeff;
    FLUX0_out *= normcoeff;

    printf("outer fluxes 0: %lf (%lf)    %lf\n", FLUX0_in, FLUX1_in, FLUX0_out);




    sprintf(fname, "%s/pup01.prof", piaacmcsimul_var.piaacmcconfdir);
    fp = fopen(fname, "w");
    for(ii=0; ii<piaacmc[0].NBradpts; ii++)
    {
        r0 = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r0lim;
        r1 = 1.0*ii/piaacmc[0].NBradpts*piaacmc[0].r1lim;
        fprintf(fp, "%f %f %g %g %g %g\n", r0, r1, pup0[ii], pup1[ii], flux0cumul[ii], flux1cumul[ii]);
    }
    fclose(fp);




    // STEP 2: COMPUTE r0 - r1 CORRESPONDANCE

    piaar00 = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts); // r0 as a function of r0 index
    piaar11 = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts); // r1 as a function of r1 index
    piaar10 = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts); // r1 as a function of r0 index
    piaar01 = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts); // r0 as a function of r1 index

    /* computing r0 and r1 */
    /* r0 and r1 are dimensionless */

    /* first, r0 is evenly distributed on the first optic */
    for(i=0; i<piaacmc[0].NBradpts; i++)
    {
        piaar00[i] = piaacmc[0].r0lim*i/piaacmc[0].NBradpts;
        piaar11[i] = piaacmc[0].r1lim*i/piaacmc[0].NBradpts;
    }

    i=0;
    ii=0;
    cnt = 0;
    piaar00[0] = 0.0;
    piaar10[0] = 0.0;
    //  fp = fopen("test0.txt", "w");
    for(i=1; i<piaacmc[0].NBradpts; i++)
    {
        F0 = flux0cumul[i];
        while((flux1cumul[ii]<flux0cumul[i])&&(ii<piaacmc[0].NBradpts))
            ii++;
        F1 = flux1cumul[ii-1];
        F2 = flux1cumul[ii];

        /* F0 = F1 + ( (F2-F1)/(ii^2-(ii-1)^2) * ((ii-1+x)^2-(ii-1)^2) ) */
        if(fabs(F2-F1)>0.0000000001)
            fluxdens = (F2-F1)/(2.0*ii-1.0);
        else
            fluxdens = 0.0000000001;
        x = sqrt((F0-F1)/fluxdens+(1.0*ii*ii-2.0*ii+1.0)) + 1.0 - 1.0*ii;

        piaar10[i] = piaacmc[0].r1lim*(1.0*ii-1.0+x)/piaacmc[0].NBradpts;
        //  fprintf(fp, "%lf %lf %lf\n", piaar00[i], piaar10[i], F0);
    }
    //  fclose(fp);



    i=0;
    ii=0;
    cnt = 0;
    piaar01[0] = 0.0;
    piaar11[0] = 0.0;
    //  fp = fopen("test1.txt", "w");
    for(i=1; i<piaacmc[0].NBradpts; i++)
    {
        F0 = flux1cumul[i];
        while((flux0cumul[ii]<flux1cumul[i])&&(ii<piaacmc[0].NBradpts))
            ii++;
        F1 = flux0cumul[ii-1];
        F2 = flux0cumul[ii];

        /* F0 = F1 + ( (F2-F1)/(ii^2-(ii-1)^2) * ((ii-1+x)^2-(ii-1)^2) ) */
        if(fabs(F2-F1)>0.0000000001)
            fluxdens = (F2-F1)/(2.0*ii-1.0);
        else
            fluxdens = 0.0000000001;
        x = sqrt((F0-F1)/fluxdens+(1.0*ii*ii-2.0*ii+1.0)) + 1.0 - 1.0*ii;

        piaar01[i] = piaacmc[0].r0lim*(1.0*ii-1.0+x)/piaacmc[0].NBradpts;
        //  fprintf(fp, "%lf %lf %lf\n", piaar11[i], piaar01[i], F0);
    }
    //  fclose(fp);





    printf("======== Compute PIAA optics shapes ============\n");
    piaaM0z = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
    piaaM1z = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);

    piaaM0z[0] = 0.0;
    piaaM1z[0] = piaacmc[0].PIAAsep;


    for(i=0; i<piaacmc[0].NBradpts-1; i++)
    {
        r0c = piaar00[i];
        r1c = piaar10[i];
        dx = (r0c-r1c)*piaacmc[0].beamrad;
        dz = piaaM1z[i]-piaaM0z[i];
   //     dist = sqrt(dx*dx+dz*dz);
        dist = dz * sqrt(1. + (dx/dz)*(dx/dz)); // preserve sign of dz
        y3 = dist - dz;
        if(fabs(dx)>0.000000001)
            slope = y3/dx;
        else
            slope = 0.0;
        r0n = piaacmc[0].r0lim*(i+1)/piaacmc[0].NBradpts;
        piaaM0z[i+1] = piaaM0z[i] + slope*(r0n-r0c)*piaacmc[0].beamrad;

        if(fabs(dx)>0.000000001)
            slope = y3/dx;
        else
            slope = 0.0;
        piaaM1z[i+1] = piaaM1z[i] + slope*(piaar10[i+1]-r1c)*piaacmc[0].beamrad;
    }

    sprintf(fname, "%s/PIAA_Mshapes.txt", piaacmcsimul_var.piaacmcconfdir);
    fp = fopen(fname, "w");
    for(ii=0; ii<piaacmc[0].NBradpts; ii++)
        fprintf(fp, "%18.16f %18.16f %18.16f %18.16f\n", piaar00[ii]*piaacmc[0].beamrad, piaaM0z[ii], piaar10[ii]*piaacmc[0].beamrad, piaaM1z[ii]);
    fclose(fp);




    free(flux0cumul);
    free(flux1cumul);
    free(pup0);
    free(pup1);

    free(piaar00);
    free(piaar10);
    free(piaar01);
    free(piaar11);


    return(0);
}




