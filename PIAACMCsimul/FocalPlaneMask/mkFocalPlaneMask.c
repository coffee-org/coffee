/**
 * @file    PIAAACMCsimul_mkFocalPlaneMask.c
 * @brief   PIAA-type coronagraph design, make focal plane mask
 *
 *
 */

// System includes
#include <malloc.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "OpticsMaterials/OpticsMaterials.h"

#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"

/**
 * @brief Make complex amplitude focal plane mask
 *
 @param[in]  IDzonemap_name	zones
 @param[in]  ID_name
 @param[in]  mode       if mode = -1, make whole 1-fpm, if mode = zone, make only 1 zone with CA = (1.0, 0.0)
 @param[in]  saveMask   1 if mask saved to file system

if mode is invalid number, no focal plane mask, AND assume 1-fpm is computed

 zone numbering starts here from 1 (zone 1 = outermost ring)
*/
errno_t mkFocalPlaneMask(const char *IDzonemap_name,
                         const char *ID_name,
                         int         mode,
                         int         saveMask,
                         imageID    *outID)
{
    DEBUG_TRACE_FSTART();

    uint32_t size;

    //  double fpscale; // [m/pix]
    uint64_t size2;

    uint_fast8_t CentCone  = 0;
    uint_fast8_t OuterCone = 0;

    double *tarray;
    double *aarray;
    double *phaarray;
    double *cosphaarray;
    double *sinphaarray;

    long NBsubPix = 64;
    int  FPMmode  = 0; // 1 if using 1-fpm

    size         = piaacmcopticalsystem.size;
    size2        = size * size;
    int nblambda = piaacmcopticalsystem.nblambda;

    imageID IDz = image_ID(IDzonemap_name);

    imageID ID;
    FUNC_CHECK_RETURN(create_3DCimage_ID(ID_name, size, size, nblambda, &ID));

    imageID IDsag;
    FUNC_CHECK_RETURN(
        create_3Dimage_ID("fpmsag", size, size, nblambda, &IDsag));

    imageID IDzone;
    FUNC_CHECK_RETURN(
        create_3Dimage_ID("fpmzone", size, size, nblambda, &IDzone));

    if(piaacmcopticaldesign.NBrings > 2)
    {
        CentCone  = 1;
        OuterCone = 1;
    }
    if(fabs(piaacmcopticaldesign.fpmOuterConeZ) < 1.0e-12)
    {
        OuterCone = 0;
    }

    printf(
        "===================== Make focal plane mask  %s %s   mode=%d/%ld   "
        "[%d %d]\n",
        IDzonemap_name,
        ID_name,
        mode,
        piaacmcopticaldesign.focmNBzone,
        CentCone,
        OuterCone);

    if(mode == -1)
    {
        FPMmode = 1;
    }
    if(mode > piaacmcopticaldesign.focmNBzone)
    {
        FPMmode = 1;
    }

    // pixel scale [m/pix] at first wavelength in array
    double fpscale =
        (2.0 * piaacmcopticaldesign.beamrad / piaacmcopticaldesign.pixscale) /
        piaacmcopticaldesign.size / piaacmcopticaldesign.fpzfactor *
        piaacmcopticalsystem.lambdaarray[0] * piaacmcopticaldesign.Fratio;
    printf(
        "piaacmcopticaldesign.fpmRad = %g m    fpscale[0] = %g m/pix   mode = "
        "%d\n",
        piaacmcopticaldesign.fpmRad,
        fpscale,
        mode);

    printf("Allocate memory\n");
    fflush(stdout);
    tarray = (double *) malloc(sizeof(double) *
                               piaacmcopticaldesign.focmNBzone * nblambda);
    if(tarray == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    aarray = (double *) malloc(sizeof(double) *
                               piaacmcopticaldesign.focmNBzone * nblambda);
    if(aarray == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    phaarray = (double *) malloc(sizeof(double) *
                                 piaacmcopticaldesign.focmNBzone * nblambda);
    if(phaarray == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    cosphaarray = (double *) malloc(sizeof(double) *
                                    piaacmcopticaldesign.focmNBzone * nblambda);
    if(cosphaarray == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    sinphaarray = (double *) malloc(sizeof(double) *
                                    piaacmcopticaldesign.focmNBzone * nblambda);
    if(sinphaarray == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    // precompute zones phase shifts
    printf("Precompute zones phase shifts  %ld zones, %d wavelengths\n",
           piaacmcopticaldesign.focmNBzone,
           nblambda);
    fflush(stdout);

    for(int k = 0; k < nblambda; k++)
        for(long zi = 0; zi < piaacmcopticaldesign.focmNBzone; zi++)
        {
            //printf("lamdba %3ld  zone %4ld   ", k, zi);
            //fflush(stdout);

            //printf("  ID=%3ld  zi=%3ld -> %4ld/%4ld :\n", piaacmcopticaldesign.zonezID, zi, piaacmcopticaldesign.focmNBzone*k+zi, piaacmcopticaldesign.focmNBzone*nblambda);
            //fflush(stdout);

            // print material thickness
            tarray[piaacmcopticaldesign.focmNBzone * k + zi] =
                data.image[piaacmcopticaldesign.zonezID].array.D[zi];
            //printf("     thickness = %8.4g m\n", tarray[piaacmcopticaldesign.focmNBzone*k+zi]);
            //fflush(stdout);

            //printf("     (ID = %3ld  zi = %3ld -> pix = %4ld/%4ld)\n", piaacmcopticaldesign.zoneaID, zi, piaacmcopticaldesign.focmNBzone*k+zi, piaacmcopticaldesign.focmNBzone*nblambda);
            //fflush(stdout);
            aarray[piaacmcopticaldesign.focmNBzone * k + zi] =
                data.image[piaacmcopticaldesign.zoneaID].array.D[zi];
            //printf("     amp = %8.4f\n", aarray[piaacmcopticaldesign.focmNBzone*k+zi]);
            //fflush(stdout);

            //
            // compute phase from thickness
            // phase sign is positive for outgoing beam element ahead of main beam
            //
            phaarray[piaacmcopticaldesign.focmNBzone * k + zi] =
                OpticsMaterials_pha_lambda(
                    piaacmcopticaldesign.fpmmaterial_code,
                    tarray[piaacmcopticaldesign.focmNBzone * k + zi],
                    piaacmcopticalsystem.lambdaarray[k]);
            //printf("     pha = %8.4f rad\n", phaarray[piaacmcopticaldesign.focmNBzone*k+zi]);
            //fflush(stdout);

            cosphaarray[piaacmcopticaldesign.focmNBzone * k + zi] =
                cosf(phaarray[piaacmcopticaldesign.focmNBzone * k + zi]);
            sinphaarray[piaacmcopticaldesign.focmNBzone * k + zi] =
                sinf(phaarray[piaacmcopticaldesign.focmNBzone * k + zi]);
        }

    //printf("Entering parallel loop\n");
    //fflush(stdout);

#ifdef HAVE_LIBGOMP
    #pragma omp parallel default(shared) private(ii,                               \
    jj,                               \
    x,                                \
    y,                                \
    r,                                \
    retmp,                            \
    imtmp,                            \
    iii,                              \
    jjj,                              \
    ii1,                              \
    jj1,                              \
    zi,                               \
    t,                                \
    a,                                \
    fpscale,                          \
    amp,                              \
    pha,                              \
    cospha,                           \
    sinpha,                           \
    ttmp,                             \
    zonetmp)
    {
        #pragma omp for
#endif
        for(int k = 0; k < nblambda; k++)
        {
            fpscale = (2.0 * piaacmcopticaldesign.beamrad /
                       piaacmcopticaldesign.pixscale) /
                      piaacmcopticaldesign.size /
                      piaacmcopticaldesign.fpzfactor *
                      piaacmcopticalsystem.lambdaarray[k] *
                      piaacmcopticaldesign.Fratio;
            printf(
                "LAMBDA %3d / %3d = %10.5g m    SCALE = %10.5g m/pix   "
                "size=%4ul  rad=%g\n",
                k,
                nblambda,
                piaacmcopticalsystem.lambdaarray[k],
                fpscale,
                size,
                piaacmcopticaldesign.fpmRad);
            printf("Zone 0 amplitude [%ld]: %lf\n",
                   piaacmcopticaldesign.zoneaID,
                   data.image[piaacmcopticaldesign.zoneaID].array.D[0]);
            printf("Zone 0 thickness: %g\n",
                   data.image[piaacmcopticaldesign.zonezID].array.D[0]);
            printf("Number of zones: %ld\n", piaacmcopticaldesign.focmNBzone);
            printf("piaacmcopticaldesign.fpmRad = %g m\n",
                   piaacmcopticaldesign.fpmRad);
            printf("piaacmcopticaldesign.fpmCentConeRad = %g m\n",
                   piaacmcopticaldesign.fpmCentConeRad);
            printf("piaacmcopticaldesign.fpmOuterConeRad [%d] = %g m\n",
                   OuterCone,
                   piaacmcopticaldesign.fpmOuterConeRad);

            for(uint32_t ii = 0; ii < size; ii++)
                for(uint32_t jj = 0; jj < size; jj++)
                {
                    //printf("[ %4ld %4ld ] ", ii, jj);

                    double x = (1.0 * ii - size / 2) * fpscale; // [m]
                    double y = (1.0 * jj - size / 2) * fpscale; // [m]
                    double r = sqrt(x * x + y * y);             // [m]

                    // default
                    double t = 0.0;
                    double a = 1.0;
                    //                  float pha = 0.0;
                    float  cospha = 1.0;
                    float  sinpha = 0.0;
                    double amp    = 1.0;

                    if(OuterCone == 1)
                    {
                        if((r > 0.9 * piaacmcopticaldesign.fpmRad) &&
                                (r < piaacmcopticaldesign
                                 .fpmOuterConeRad)) // outer cone
                        {
                            double t =
                                piaacmcopticaldesign.fpmOuterConeZ *
                                (piaacmcopticaldesign.fpmOuterConeRad - r) /
                                (piaacmcopticaldesign.fpmOuterConeRad -
                                 piaacmcopticaldesign.fpmRad);
                            double pha = OpticsMaterials_pha_lambda(
                                             piaacmcopticaldesign.fpmmaterial_code,
                                             t,
                                             piaacmcopticalsystem.lambdaarray[k]);
                            cospha = cosf(pha);
                            sinpha = sinf(pha);
                            a      = 1.0;
                            amp    = a;
                        }
                    }

                    data.image[IDzone].array.F[k * size2 + jj * size + ii] = 0;

                    if(r <
                            1.1 * piaacmcopticaldesign.fpmRad) //     fine sampling
                    {
                        //       printf("pix outercone ...");
                        //      fflush(stdout);

                        double retmp   = 0.0;
                        double imtmp   = 0.0;
                        double ttmp    = 0.0;
                        double zonetmp = 0.0;
                        for(long iii = 0; iii < NBsubPix; iii++)
                        {
                            for(long jjj = 0; jjj < NBsubPix; jjj++)
                            {
                                // physical coordinates on mask
                                // x and y in [m]
                                double x =
                                    (1.0 * ii - size / 2 +
                                     1.0 * (0.5 + iii) / NBsubPix - 0.5) *
                                    fpscale;
                                double y =
                                    (1.0 * jj - size / 2 +
                                     1.0 * (0.5 + jjj) / NBsubPix - 0.5) *
                                    fpscale;
                                double r = sqrt(x * x + y * y); // [m]

                                long zi = 0; // default
                                cospha  = 1.0;
                                sinpha  = 0.0;
                                amp     = 1.0;

                                if(OuterCone == 1)
                                    if((r >
                                            0.9 * piaacmcopticaldesign.fpmRad) &&
                                            (r <
                                             piaacmcopticaldesign
                                             .fpmOuterConeRad)) // outer cone
                                    {
                                        t = piaacmcopticaldesign.fpmOuterConeZ *
                                            (piaacmcopticaldesign
                                             .fpmOuterConeRad -
                                             r) /
                                            (piaacmcopticaldesign
                                             .fpmOuterConeRad -
                                             piaacmcopticaldesign.fpmRad);
                                        a = 1.0;
                                        cospha =
                                            cosphaarray[piaacmcopticaldesign
                                                        .focmNBzone *
                                                        k +
                                                        zi - 1];
                                        sinpha =
                                            sinphaarray[piaacmcopticaldesign
                                                        .focmNBzone *
                                                        k +
                                                        zi - 1];
                                        amp = a;
                                    }

                                long ii1 =
                                    (long)((0.5 +
                                            0.5 * x /
                                            piaacmcopticaldesign.fpmRad *
                                            piaacmcparams.FPMSCALEFACTOR) *
                                           piaacmcopticaldesign
                                           .fpmarraysize +
                                           0.5);
                                long jj1 =
                                    (long)((0.5 +
                                            0.5 * y /
                                            piaacmcopticaldesign.fpmRad *
                                            piaacmcparams.FPMSCALEFACTOR) *
                                           piaacmcopticaldesign
                                           .fpmarraysize +
                                           0.5);

                                if((ii1 > -1) &&
                                        (ii1 < piaacmcopticaldesign.fpmarraysize) &&
                                        (jj1 > -1) &&
                                        (jj1 < piaacmcopticaldesign.fpmarraysize))
                                {
                                    if(CentCone == 1)
                                    {
                                        // central cone
                                        if((r <
                                                0.99 *
                                                piaacmcopticaldesign.fpmRad) &&
                                                (r < piaacmcopticaldesign
                                                 .fpmCentConeRad))
                                        {
                                            t = piaacmcopticaldesign
                                                .fpmCentConeZ +
                                                (r / piaacmcopticaldesign
                                                 .fpmCentConeRad) *
                                                (0.5 * (piaacmcopticaldesign
                                                        .fpmminsag +
                                                        piaacmcopticaldesign
                                                        .fpmmaxsag) -
                                                 piaacmcopticaldesign
                                                 .fpmCentConeZ);
                                            // piaacmcopticaldesign.fpmCentConeZ*(piaacmcopticaldesign.fpmCentConeRad-r)/(piaacmcopticaldesign.fpmCentConeRad); //piaacmcopticaldesign.fpmCentConeZ
                                            a = 1.0;
                                            double pha =
                                                OpticsMaterials_pha_lambda(
                                                    piaacmcopticaldesign
                                                    .fpmmaterial_code,
                                                    t,
                                                    piaacmcopticalsystem
                                                    .lambdaarray[k]);
                                            cospha = cosf(pha);
                                            sinpha = sinf(pha);
                                            amp    = a;
                                        }
                                    }

                                    // Zone number
                                    zi = (long)(data.image[IDz].array.UI16
                                                [jj1 * piaacmcopticaldesign
                                                     .fpmarraysize +
                                                     ii1]);
                                    if(zi - 1 >
                                            data.image[piaacmcopticaldesign.zonezID]
                                            .md[0]
                                            .size[0] -
                                            1)
                                    {
                                        printf(
                                            "ERROR: Zone %d does not exist "
                                            "(image %s has size %ld %ld)   pix "
                                            "%ld "
                                            "%ld   %ld\n",
                                            (int) zi,
                                            data.image[piaacmcopticaldesign
                                                       .zonezID]
                                            .md[0]
                                            .name,
                                            (long) data
                                            .image[piaacmcopticaldesign
                                                   .zonezID]
                                            .md[0]
                                            .size[0],
                                            (long) data
                                            .image[piaacmcopticaldesign
                                                   .zonezID]
                                            .md[0]
                                            .size[1],
                                            (long) jj1,
                                            (long) jj1,
                                            piaacmcopticaldesign.fpmarraysize);
                                        exit(0);
                                    }
                                    if(zi > 0)
                                    {
                                        t = data.image[piaacmcopticaldesign
                                                       .zonezID]
                                            .array.D[zi - 1]; // thickness
                                        a = data.image[piaacmcopticaldesign
                                                       .zoneaID]
                                            .array
                                            .D[zi -
                                                  1]; // amplitude transmission
                                        cospha =
                                            cosphaarray[piaacmcopticaldesign
                                                        .focmNBzone *
                                                        k +
                                                        zi - 1];
                                        sinpha =
                                            sinphaarray[piaacmcopticaldesign
                                                        .focmNBzone *
                                                        k +
                                                        zi - 1];
                                        amp = a;
                                    }
                                }

                                if(FPMmode == 1)  // make 1-fpm
                                {
                                    /*                                    amp = a;
                                                                        pha = OpticsMaterials_pha_lambda(piaacmcopticaldesign.fpmmaterial_code, t, piaacmcopticalsystem.lambdaarray[k]);
                                                                        cospha = cosf(pha);
                                                                        sinpha = sinf(pha);*/

                                    //									cospha = cosphaarray[piaacmcopticaldesign.focmNBzone*k+zi-1];
                                    //									sinpha = sinphaarray[piaacmcopticaldesign.focmNBzone*k+zi-1];

                                    retmp += 1.0 - amp * cospha;
                                    imtmp += -amp * sinpha;
                                }
                                else // impulse response from single zone
                                {
                                    if(mode == zi)
                                    {
                                        amp = 1.0;
                                        //                                      pha = OPTICSMATERIALS_pha_lambda(piaacmcopticaldesign.fpmmaterial_code, t, piaacmcopticalsystem.lambdaarray[k]);
                                        //										cospha = cosf(pha);
                                        //										sinpha = sinf(pha);
                                        retmp += amp;
                                        imtmp += 0.0;
                                    }
                                }

                                //bad location
                                ttmp += t;
                                zonetmp += 1.0 * zi;
                            }
                        }

                        data.image[ID].array.CF[k * size2 + jj * size + ii].re =
                            retmp / (NBsubPix * NBsubPix);
                        data.image[ID].array.CF[k * size2 + jj * size + ii].im =
                            imtmp / (NBsubPix * NBsubPix);
                        data.image[IDsag].array.F[k * size2 + jj * size + ii] =
                            ttmp / (NBsubPix * NBsubPix);
                        data.image[IDzone].array.F[k * size2 + jj * size + ii] =
                            zonetmp / (NBsubPix * NBsubPix);
                    }
                    else // coarse sampling, outside zones
                    {
                        if(FPMmode == 1)  // make 1-fpm
                        {
                            data.image[ID]
                            .array.CF[k * size2 + jj * size + ii]
                            .re = 1.0 - amp * cospha;
                            data.image[ID]
                            .array.CF[k * size2 + jj * size + ii]
                            .im = -amp * sinpha;
                        }
                        else // single zone
                        {
                            data.image[ID]
                            .array.CF[k * size2 + jj * size + ii]
                            .re = 0.0;
                            data.image[ID]
                            .array.CF[k * size2 + jj * size + ii]
                            .im = 0.0;
                        }

                        data.image[IDsag].array.F[k * size2 + jj * size + ii] =
                            t;
                        data.image[IDzone].array.F[k * size2 + jj * size + ii] =
                            0.0;
                    }
                }
        }
#ifdef HAVE_LIBGOMP
    }
#endif

    printf("===================== focal plane mask  : Done\n");
    fflush(stdout);

    if(saveMask == 1)
{
    /* save mask sag */
    //save_fits("fpmsag", "tmp_fpmsag.fits");

    PIAACMCsimul_update_fnamedescr();

        char fname[STRINGMAXLEN_FULLFILENAME];

        WRITE_FULLFILENAME(fname,
                           "%s/fpm_sagmap2D.%s.fits.gz",
                           piaacmcparams.piaacmcconfdir,
                           piaacmcparams.fnamedescr);
        FUNC_CHECK_RETURN(save_fits("fpmsag", fname));

        /* save zones */
        WRITE_FULLFILENAME(fname,
                           "%s/fpm_zonemap2D.%s.fits.gz",
                           piaacmcparams.piaacmcconfdir,
                           piaacmcparams.fnamedescr);
        FUNC_CHECK_RETURN(save_fits("fpmzone", fname));

        /* save mask transmission */
        imageID IDm;
        FUNC_CHECK_RETURN(
            create_3DCimage_ID("fpmCA", size, size, nblambda, &IDm));

        for(int k = 0; k < nblambda; k++)
            for(uint32_t ii = 0; ii < size; ii++)
                for(uint32_t jj = 0; jj < size; jj++)
                {
                    data.image[IDm].array.CF[k * size2 + jj * size + ii].re =
                        1.0 -
                        data.image[ID].array.CF[k * size2 + jj * size + ii].re;
                    data.image[IDm].array.CF[k * size2 + jj * size + ii].im =
                        -data.image[ID].array.CF[k * size2 + jj * size + ii].im;
                    // [re,im] = 1-fpm  -> fpm = [(1.0-re), -im]
                }

        mk_amph_from_complex("fpmCA", "tfpma", "tfpmp", 0);
        FUNC_CHECK_RETURN(
            delete_image_ID("fpmCA", DELETE_IMAGE_ERRMODE_WARNING));

        WRITE_FULLFILENAME(fname,
                           "%s/fpm_CAampmap2D.%s.fits.gz",
                           piaacmcparams.piaacmcconfdir,
                           piaacmcparams.fnamedescr);
        FUNC_CHECK_RETURN(save_fits("tfpma", fname));

        WRITE_FULLFILENAME(fname,
                           "%s/fpm_CAphamap2D.%s.fits.gz",
                           piaacmcparams.piaacmcconfdir,
                           piaacmcparams.fnamedescr);
        FUNC_CHECK_RETURN(save_fits("tfpmp", fname));

        FUNC_CHECK_RETURN(
            delete_image_ID("tfpma", DELETE_IMAGE_ERRMODE_WARNING));
        FUNC_CHECK_RETURN(
            delete_image_ID("tfpmp", DELETE_IMAGE_ERRMODE_WARNING));
    }

    FUNC_CHECK_RETURN(delete_image_ID("fpmsag", DELETE_IMAGE_ERRMODE_WARNING));

    free(tarray);
    free(aarray);
    free(phaarray);
    free(cosphaarray);
    free(sinphaarray);

    if(outID != NULL)
{
    *outID = ID;
}

DEBUG_TRACE_FEXIT();
return RETURN_SUCCESS;
}
