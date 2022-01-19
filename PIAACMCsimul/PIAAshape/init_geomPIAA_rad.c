/**
 * @file    PIAACMCsimul_init_geomPIAA_rad.c
 * @brief   PIAA-type coronagraph design
 *
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_memory/COREMOD_memory.h"

#include "PIAACMCsimul.h"

/**
 * @brief Computes radial PIAA mirror shapes
 *
 * Only works for circular PIAA.
 *
 * @param[in] IDapofit_name
 *      Radial amplitude apodization, expressed as linear sum of cosines
 *
 * @return errno_t
 */
errno_t init_geomPIAA_rad(const char *__restrict__ IDapofit_name)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT_LOG("FARG %s", IDapofit_name);

    long nbcoeff;
    double total;

    // to convert r ro r1 (assymptotic outer radius on pup1)
    double coeffa = 3.0;  // convergence rate from linear to assymptotic value
    double coeffa1 = 0.5; // convergence radius limit (added to 1.0)

    double FLUX0_in = 0.0;  // inside central obstruction
    double FLUX0_out = 0.0; // outside beam edge
    double FLUX1_in = 0.0;  // inside central obstruction
    double FLUX1_out = 0.0; // outside beam edge
    double normcoeff;

    //double *piaar11;
    double *piaar01;
    double *piaar10;

    double *piaaM0z;
    double *piaaM1z;

    double *pup0 = (double *)malloc(sizeof(double) * piaacmcopticaldesign.NBradpts);
    if (pup0 == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }

    double *pup1 = (double *)malloc(sizeof(double) * piaacmcopticaldesign.NBradpts);
    if (pup1 == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }

    double *flux0cumul = (double *)malloc(sizeof(double) * piaacmcopticaldesign.NBradpts);
    if (flux0cumul == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }

    double *flux1cumul = (double *)malloc(sizeof(double) * piaacmcopticaldesign.NBradpts);
    if (flux1cumul == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }

    // STEP 1: CREATE AMPLITUDE AND CUMULATIVE INTENSITY PROFILES

    // CREATE OUTPUT AMPLITUDE APODIZATION PROFILE AND ITS CUMUL

    imageID IDcoeff = image_ID(IDapofit_name);
    nbcoeff = data.image[IDcoeff].md->size[0];
    printf("%ld coefficients\n", nbcoeff);

    total = 0.0;
    for (long ii = 0; ii < piaacmcopticaldesign.NBradpts; ii++)
    {
        double r1;

        pup1[ii] = 0.0;
        double r = 1.0 * ii / piaacmcopticaldesign.NBradpts * piaacmcopticaldesign.r1lim;
        if (r < 1.0)
            r1 = r;
        else
            r1 = 1.0 + (r - 1) / pow((1.0 + pow(1.0 / coeffa1 * (r - 1), coeffa)), 1.0 / coeffa);

        // reconstruct apodization profile from cosine coefficients
        for (long k = 0; k < nbcoeff; k++)
        {
            pup1[ii] += data.image[IDcoeff].array.F[k] * cos(r1 * k * M_PI / ApoFitCosFact);
        }

        // sum up flux inside central obstruction
        if (r < piaacmcopticaldesign.centObs1)
        {
            FLUX1_in += pup1[ii] * pup1[ii] * r;
        }

        // sum up flux outside pupil edge
        if (r > 1.0)
        {
            FLUX1_out += pup1[ii] * pup1[ii] * r;
        }

        total += pup1[ii] * pup1[ii] * r;
        flux1cumul[ii] = total;
    }

    normcoeff = 1.0 / (total - FLUX1_in - FLUX1_out);

    FLUX1_in *= normcoeff;
    FLUX1_out *= normcoeff;
    for (long ii = 0; ii < piaacmcopticaldesign.NBradpts; ii++)
    {
        // r = 1.0*ii/piaacmc->NBradpts*piaacmc->r1lim;
        flux1cumul[ii] *= normcoeff;
    }

    printf("outer fluxes 1: %lf %lf\n", FLUX1_in, FLUX1_out);

    // CREATE FLUX0

    total = 0.0;
    for (long ii = 0; ii < piaacmcopticaldesign.NBradpts; ii++)
    {
        double r = 1.0 * ii / piaacmcopticaldesign.NBradpts * piaacmcopticaldesign.r0lim;
        pup0[ii] = 1.0;

        if (r < piaacmcopticaldesign.centObs0)
            FLUX0_in += pup0[ii] * pup0[ii] * r;
        if (r > 1.0)
            FLUX0_out += pup0[ii] * pup0[ii] * r;

        total += pup0[ii] * pup0[ii] * r;
        flux0cumul[ii] = total;
    }
    normcoeff = 1.0 / (total - FLUX0_in - FLUX0_out);

    FLUX0_in *= normcoeff;
    FLUX0_out *= normcoeff;

    printf("outer fluxes 0: %lf (%lf)    %lf\n", FLUX0_in, FLUX1_in, FLUX0_out);

    //
    // Compute inner pseudo profile
    //
    {
        double a = 0.5;
        double b = 0.5;
        double bstep = 0.1;
        double verr = 1.0;
        double eps1 = 1e-8;
        long NBistep = piaacmcopticaldesign.centObs0 * piaacmcopticaldesign.NBradpts / piaacmcopticaldesign.r0lim;
        //  innerprof_cumul = (double*) malloc(sizeof(double)*NBistep);

        int dir = 1; // initial direction

        double t0;
        while (fabs(verr) > 1.0e-9)
        {
            int odir;
            t0 = 0.0;
            //t0cnt = 0.0;
            for (long ii = 0; ii < NBistep; ii++)
            {
                double x = 1.0 * ii / NBistep;
                double r = 1.0 * ii / piaacmcopticaldesign.NBradpts * piaacmcopticaldesign.r0lim;

                if (x < eps1)
                    x = eps1;
                if (x > 1.0 - eps1)
                    x = 1.0 - eps1;

                pup0[ii] = b + (1.0 - b) * (0.5 + atan(-a / b / (x * x) + a / pow(x - 1.0, 2)) / M_PI);
                t0 += r * pup0[ii] * pup0[ii];
                flux0cumul[ii] = t0;
            }

            verr = t0 * normcoeff - FLUX1_in;

            odir = dir;
            if (verr > 0.0) // too much light
            {
                b /= (1.0 + bstep);
                dir = -1;
            }
            else
            {
                b *= (1.0 + bstep);
                dir = 1;
            }
            if (odir * dir < 0)
                bstep *= 0.1;
            printf(".");
            fflush(stdout);
        }
        printf("\n");
        printf("TOTAL = %f -> %g (%g %g)\n", b, t0 * normcoeff, bstep, verr);
    }

    // outer region
    {
        double a = 0.5;
        double b = 0.5;
        double bstep = 0.1;
        double verr = 1.0;
        double eps1 = 1e-8;

        long NBistep = piaacmcopticaldesign.NBradpts * (piaacmcopticaldesign.r0lim - 1.0) / piaacmcopticaldesign.r0lim;
        //  innerprof_cumul = (double*) malloc(sizeof(double)*NBistep);
        long iioffset = (long)(1.0 * piaacmcopticaldesign.NBradpts / piaacmcopticaldesign.r0lim);
        NBistep = piaacmcopticaldesign.NBradpts - iioffset;

        int dir = 1; // initial direction
        int odir;

        double t0;
        while (fabs(verr) > 1.0e-9)
        {
            t0 = 0.0;
            //t0cnt = 0.0;
            for (long ii = 0; ii < NBistep; ii++)
            {
                double x = 1.0 - 1.0 * ii / NBistep;
                double r = 1.0 + 1.0 * ii / piaacmcopticaldesign.NBradpts * piaacmcopticaldesign.r0lim;

                if (x < eps1)
                    x = eps1;
                if (x > 1.0 - eps1)
                    x = 1.0 - eps1;
                pup0[ii + iioffset] = b + (1.0 - b) * (0.5 + atan(-a / b / (x * x) + a / pow(x - 1.0, 2)) / M_PI);
                t0 += r * pup0[ii + iioffset] * pup0[ii + iioffset];
                flux0cumul[ii + iioffset] = t0;
            }

            verr = t0 * normcoeff - FLUX1_out;

            odir = dir;
            if (verr > 0.0) // too much light
            {
                b /= (1.0 + bstep);
                dir = -1;
            }
            else
            {
                b *= (1.0 + bstep);
                dir = 1;
            }
            if (odir * dir < 0)
                bstep *= 0.1;
            printf(".");
            fflush(stdout);
        }
        printf("\n");
        printf("TOTAL = %f -> %g (%g %g)\n", b, t0 * normcoeff, bstep, verr);
    }

    total = 0.0;
    FLUX0_in = 0.0;
    FLUX0_out = 0.0;
    for (long ii = 0; ii < piaacmcopticaldesign.NBradpts; ii++)
    {
        double r = 1.0 * ii / piaacmcopticaldesign.NBradpts * piaacmcopticaldesign.r0lim;
        if (r < piaacmcopticaldesign.centObs0)
            FLUX0_in += pup0[ii] * pup0[ii] * r;
        if (r > 1.0)
            FLUX0_out += pup0[ii] * pup0[ii] * r;

        total += pup0[ii] * pup0[ii] * r;
        flux0cumul[ii] = total;
        flux0cumul[ii] *= normcoeff;
        ;
    }
    FLUX0_in *= normcoeff;
    FLUX0_out *= normcoeff;

    printf("outer fluxes 0: %lf (%lf)    %lf\n", FLUX0_in, FLUX1_in, FLUX0_out);

    { // write beam amplitude profiles to file
        //
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/pup01.prof", piaacmcparams.piaacmcconfdir);
        FILE *fp = fopen(fname, "w");
        if (fp == NULL)
        {
            PRINT_ERROR("Cannot create file %s", fname);
            abort();
        }
        for (long ii = 0; ii < piaacmcopticaldesign.NBradpts; ii++)
        {
            double r0 = 1.0 * ii / piaacmcopticaldesign.NBradpts * piaacmcopticaldesign.r0lim;
            double r1 = 1.0 * ii / piaacmcopticaldesign.NBradpts * piaacmcopticaldesign.r1lim;

            fprintf(fp, "%f %f %g %g %g %g\n", r0, r1, pup0[ii], pup1[ii], flux0cumul[ii], flux1cumul[ii]);
        }
        fclose(fp);
    }

    // STEP 2: COMPUTE r0 - r1 CORRESPONDANCE

    double *piaar00 = (double *)malloc(sizeof(double) * piaacmcopticaldesign.NBradpts); // r0 as a function of r0 index
    if (piaar00 == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }
    /*
        piaar11 = (double*) malloc(sizeof(double)*piaacmc->NBradpts); // r1 as a function of r1 index
        if(piaar11 == NULL) {
            PRINT_ERROR("malloc returns NULL pointer");
            abort(); // or handle error in other ways
        }
    */
    piaar10 = (double *)malloc(sizeof(double) * piaacmcopticaldesign.NBradpts); // r1 as a function of r0 index
    if (piaar10 == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }

    piaar01 = (double *)malloc(sizeof(double) * piaacmcopticaldesign.NBradpts); // r0 as a function of r1 index
    if (piaar01 == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }
    /* computing r0 and r1 */
    /* r0 and r1 are dimensionless */

    /* first, r0 is evenly distributed on the first optic */
    for (long i = 0; i < piaacmcopticaldesign.NBradpts; i++)
    {
        piaar00[i] = piaacmcopticaldesign.r0lim * i / piaacmcopticaldesign.NBradpts;
        //piaar11[i] = piaacmc->r1lim*i/piaacmc->NBradpts;
    }

    piaar00[0] = 0.0;
    piaar10[0] = 0.0;
    //  fp = fopen("test0.txt", "w");
    {
        long ii = 0;
        for (long i = 1; i < piaacmcopticaldesign.NBradpts; i++)
        {
            double F0 = flux0cumul[i];

            while ((ii < piaacmcopticaldesign.NBradpts) && (flux1cumul[ii] < flux0cumul[i]))
            {
                ii++;
            }

            double F1 = flux1cumul[ii - 1];
            double F2 = flux1cumul[ii];

            /* F0 = F1 + ( (F2-F1)/(ii^2-(ii-1)^2) * ((ii-1+x)^2-(ii-1)^2) ) */
            double fluxdens;
            if (fabs(F2 - F1) > 0.0000000001)
                fluxdens = (F2 - F1) / (2.0 * ii - 1.0);
            else
                fluxdens = 0.0000000001;

            double x = sqrt((F0 - F1) / fluxdens + (1.0 * ii * ii - 2.0 * ii + 1.0)) + 1.0 - 1.0 * ii;

            piaar10[i] = piaacmcopticaldesign.r1lim * (1.0 * ii - 1.0 + x) / piaacmcopticaldesign.NBradpts;
            //  fprintf(fp, "%lf %lf %lf\n", piaar00[i], piaar10[i], F0);
        }
    }
    //  fclose(fp);

    piaar01[0] = 0.0;
    //piaar11[0] = 0.0;
    //  fp = fopen("test1.txt", "w");
    {
        long ii = 0;
        for (long i = 1; i < piaacmcopticaldesign.NBradpts; i++)
        {
            double F0 = flux1cumul[i];

            while ((ii < piaacmcopticaldesign.NBradpts) && (flux0cumul[ii] < flux1cumul[i]))
            {
                ii++;
            }

            double F1 = flux0cumul[ii - 1];
            double F2 = flux0cumul[ii];

            /* F0 = F1 + ( (F2-F1)/(ii^2-(ii-1)^2) * ((ii-1+x)^2-(ii-1)^2) ) */
            double fluxdens;
            if (fabs(F2 - F1) > 0.0000000001)
                fluxdens = (F2 - F1) / (2.0 * ii - 1.0);
            else
                fluxdens = 0.0000000001;

            double x = sqrt((F0 - F1) / fluxdens + (1.0 * ii * ii - 2.0 * ii + 1.0)) + 1.0 - 1.0 * ii;

            piaar01[i] = piaacmcopticaldesign.r0lim * (1.0 * ii - 1.0 + x) / piaacmcopticaldesign.NBradpts;
            //  fprintf(fp, "%lf %lf %lf\n", piaar11[i], piaar01[i], F0);
        }
    }
    //  fclose(fp);

    printf("======== Compute PIAA optics shapes ============\n");

    piaaM0z = (double *)malloc(sizeof(double) * piaacmcopticaldesign.NBradpts);
    if (piaaM0z == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }

    piaaM1z = (double *)malloc(sizeof(double) * piaacmcopticaldesign.NBradpts);
    if (piaaM1z == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }

    piaaM0z[0] = 0.0;
    piaaM1z[0] = piaacmcopticaldesign.PIAAsep;

    for (long i = 0; i < piaacmcopticaldesign.NBradpts - 1; i++)
    {
        double r0c = piaar00[i];
        double r1c = piaar10[i];
        double dx = (r0c - r1c) * piaacmcopticaldesign.beamrad;
        double dz = piaaM1z[i] - piaaM0z[i];
        //     dist = sqrt(dx*dx+dz*dz);
        double dist = dz * sqrt(1. + (dx / dz) * (dx / dz)); // preserve sign of dz
        double y3 = dist - dz;

        double slope;
        if (fabs(dx) > 0.000000001)
            slope = y3 / dx;
        else
            slope = 0.0;

        double r0n = piaacmcopticaldesign.r0lim * (i + 1) / piaacmcopticaldesign.NBradpts;

        piaaM0z[i + 1] = piaaM0z[i] + slope * (r0n - r0c) * piaacmcopticaldesign.beamrad;

        if (fabs(dx) > 0.000000001)
            slope = y3 / dx;
        else
            slope = 0.0;

        piaaM1z[i + 1] = piaaM1z[i] + slope * (piaar10[i + 1] - r1c) * piaacmcopticaldesign.beamrad;
    }

    { // write PIAA optics radial shapes to file

        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/PIAA_Mshapes.txt", piaacmcparams.piaacmcconfdir);
        FILE *fp = fopen(fname, "w");
        if (fp != NULL)
        {
            for (long ii = 0; ii < piaacmcopticaldesign.NBradpts; ii++)
            {
                fprintf(fp, "%18.16f %18.16f %18.16f %18.16f\n", piaar00[ii] * piaacmcopticaldesign.beamrad,
                        piaaM0z[ii], piaar10[ii] * piaacmcopticaldesign.beamrad, piaaM1z[ii]);
            }
            fclose(fp);
        }
        else
        {
            FUNC_RETURN_FAILURE("cannot write file %s", fname);
        }
    }

    free(piaaM0z);
    free(piaaM1z);

    free(flux0cumul);
    free(flux1cumul);
    free(pup0);
    free(pup1);

    free(piaar00);
    free(piaar10);
    free(piaar01);
    //free(piaar11);

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
