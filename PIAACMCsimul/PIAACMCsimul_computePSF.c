/**
 * @file    PIAACMCsimul_computePSF.c
 * @brief   PIAA-type coronagraph design
 *
 * Can design both APLCMC and PIAACMC coronagraphs
 *
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "OptSystProp/OptSystProp.h"
#include "linopt_imtools/linopt_imtools.h"

#include "PIAACMCsimul.h"

#include "PIAACMCsimul_achromFPMsol_eval.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"
#include "init_piaacmcopticalsystem.h"

#include "PIAAshape/makePIAAshapes.h"

/**
 * @brief Compute PSF
 *
 * @return Average contrast in evaluation zone
 *
 *
 *
 * Source is defined by parameters sourcesize and extmode :
 * - source size = 1e-{sourcesize*0.1}, except if sourcesize = 0 (point source)
 * - sourcesize is a 2-digit number ( 10 = 0.1 l/D, 20 = 0.01 l/D etc..)
 * - extmode = 0 : 1 point (point source)
 * - extmode = 1 : 3 point sources, 120 apart on circle radius = source size
 * - extmode = 2 : 6 point sources. 3 as above on circle radius 1/sqrt(2.5) + 3 on outer circle, radius 2/sqrt(2.5), 120 apart, clockled 60 deg off inner points
 *
 * @note If opderrcube exists, include each slice as a WF mode
 *
 * PSF is held in shared memory by default
 *
 * ---
 *
 * ### Output
 *
 * | name                              |  type      | Description                            |
 * |-----------------------------------|------------|----------------------------------------|
 * | scoringmask                       | 2D image   | focal plane points used for evaluation |
 * | <piaacmcdir>/scoringmask<N>.fits  | 2D FITS    | focal plane points used for evaluation |
 * | imvec                             | 1D image   | output vector                          |
 * | psfi0                             | 3D image   | output PSF                             |
 *
 *
 * ---
 * ---
 *
 *
 */
errno_t PIAACMCsimul_computePSF(
    float xld, /// @param[in]   xld         float: Source X position [l/D]
    float yld, /// @param[in]   yld         float: Source Y position [l/D]
    long
        startelem, /// @param[in]   startelem   long : First element in propagation
    long
        endelem, /// @param[in]   endelem     long : Last element in propagation
    int savepsf, /// @param[in]   savepsf     int  : Save PSF flag
    int sourcesize, /// @param[in]   sourcezise  int  : Source size (10x log10)
    int extmode,    /// @param[in]   extmode     int  : Source extended type
    int outsave,    /// @param[in]   outsave     int  : Save output flag
    double *contrastval /// @param[out]  contrastval *double : contrast value
)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %f %f", xld, yld);

    // how to measure quality
    float focscale; // l/D per pix
    float scoringIWA   = 1.5;
    float scoringOWA   = 20.0;
    float scoringOWAhr = 8.0;
    float scoringIWAx  = -20.5;

    double avContrast;
    double peakcontrast;

    double normcoeff = 1.0;

    (void) savepsf;
    (void) endelem;

    // size of one side of each image array
    uint32_t size = piaacmcopticaldesign.size;

    // load an error if it exists
    imageID IDopderrC = image_ID("OPDerrC");
    long    nbOPDerr  = 0;
    if (IDopderrC == -1)
    {
        load_fits("OPDerrC.fits", "OPDerrC", 0, &IDopderrC);
    }

    if (IDopderrC != -1)
    {
        uint8_t naxis = data.image[IDopderrC].md[0].naxis;
        if (naxis == 2)
        {
            nbOPDerr =
                data.image[IDopderrC].md[0].size[2]; // number of error arrays
        }
        printf("INCLUDING %ld OPD ERROR MODES\n", nbOPDerr);
        fflush(stdout);
    }

    // focal plane plate scale in lambda/D per pixel
    focscale =
        (2.0 * piaacmcopticaldesign.beamrad / piaacmcopticaldesign.pixscale) /
        piaacmcopticaldesign.size;

    /// ### Create scoring mask if it doesn't exist
    /// The scoring mask is the array of evaluation points on the focal plane
    {
        if (image_ID("scoringmask") == -1)
        {
            printf("CREATING SCORING MASK\n");
            printf("FOCAL PLANE SCALE = %f l/d per pix\n", focscale);
            fflush(stdout);
            imageID IDsm;
            create_2Dimage_ID("scoringmask", size, size, &IDsm);

            if (piaacmcparams.SCORINGMASKTYPE == 0) // high density, wide
            {
                // draw an array of points, clip to desired subregions
                for (uint32_t ii = 0; ii < size; ii++)
                    for (uint32_t jj = 0; jj < size; jj++)
                    {
                        // create regular array and the raduis of each point
                        double x = (1.0 * ii - 0.5 * size) * focscale;
                        double y = (1.0 * jj - 0.5 * size) * focscale;
                        double r = sqrt(x * x + y * y);

                        // clip the regular array to the desired annulus inside high-resolution
                        // region defined by scoringOWAhr
                        // use every other point
                        // and clip to an x > scoringIWAx part of the annulus if desired
                        if ((r > scoringIWA) && (r < scoringOWAhr) &&
                            (x > scoringIWAx) && ((ii + jj) % 2 == 0))
                        {
                            data.image[IDsm].array.F[jj * size + ii] = 1.0;
                        }
                        // pick every other row and column between scoringOWAhr and scoringOWA
                        if ((r > scoringOWAhr) && (r < scoringOWA) &&
                            (x > scoringIWAx) && (ii % 2 == 0) && (jj % 2 == 0))
                        {
                            data.image[IDsm].array.F[jj * size + ii] = 1.0;
                        }
                        // draw a single radial line of points out to IWA = 70 (every other point)
                        if ((x > scoringOWA) && (fabs(y) < scoringIWA * 0.25) &&
                            (r < 70.0) && ((ii + jj) % 2 == 0)) // single line
                        {
                            data.image[IDsm].array.F[jj * size + ii] = 1.0;
                        }
                    }
            }
            else // focused on central pixels, fewer pixels for faster convergence - used for FPMresp based optimization
            {
                for (uint32_t ii = 0; ii < size; ii++)
                    for (uint32_t jj = 0; jj < size; jj++)
                    {
                        double x = (1.0 * ii - 0.5 * size) * focscale;
                        double y = (1.0 * jj - 0.5 * size) * focscale;
                        double r = sqrt(x * x + y * y);
                        // clip from scoringIWA to scoringOWAhr only using every other column and row
                        if ((r > scoringIWA) && (r < scoringOWAhr) &&
                            (x > scoringIWAx) && (ii % 2 == 0) && (jj % 2 == 0))
                        {
                            data.image[IDsm].array.F[jj * size + ii] = 1.0;
                        }
                    }
            }
            if (piaacmcparams.PIAACMC_save == 1)
            {
                // save a disgnostic image

                char fname[STRINGMAXLEN_FULLFILENAME];

                WRITE_FULLFILENAME(fname,
                                   "%s/scoringmask%d.fits",
                                   piaacmcparams.piaacmcconfdir,
                                   piaacmcparams.SCORINGMASKTYPE);

                save_fits("scoringmask", fname);
            }
            // a pixtable is a list of non-zero pixels with their coordinates
            FUNC_CHECK_RETURN(linopt_imtools_mask_to_pixtable("scoringmask",
                                                              "pixindex",
                                                              "pixmult",
                                                              NULL));

            // sums the image, giving the total number of pixels in the scoring mask
            piaacmcparams.SCORINGTOTAL = arith_image_total("scoringmask");

            //exit(0);
        }
    }

    /// ## Fast PSF computattion (if piaacmcparams.computePSF_FAST_FPMresp = 1)

    /// @note Only possible if mode 11 has already been executed
    if (piaacmcparams.computePSF_FAST_FPMresp == 1)
    {
        /// Compute the PSF as the complex amplitude for the evaluation points on the focal plane
        /// for a given FPM zone thickness based on the FPMresp array computed in mode 11
        PIAACMCsimul_achromFPMsol_eval(
            piaacmcparams.fpmresp_array,
            piaacmcparams.zonez_array,
            piaacmcparams.dphadz_array,
            piaacmcparams.outtmp_array,
            piaacmcparams.vsize,
            data.image[piaacmcopticaldesign.zonezID].md[0].size[0],
            piaacmcopticalsystem.nblambda,
            NULL);

        //		printf("FAST FPMresp CALLED\n");
        //		sleep(1000000);//TEST

        ///
        /// PSF result is stored in outtmp_array

        //peakcontrast = 0.0;

        {
            imageID ID = image_ID(
                "imvect"); /// - Use \c imvect for storage if it exists, or create it
            if (ID == -1)
            {
                create_2Dimage_ID("imvect",
                                  piaacmcparams.vsize *
                                      piaacmcopticalsystem.nblambda,
                                  1,
                                  &ID);
            }
            /// - Write the result into \c imvect
            double value = 0.0;
            for (uint32_t ii = 0;
                 ii < piaacmcparams.vsize * piaacmcopticalsystem.nblambda;
                 ii++) // for each wavelength
            {
                data.image[ID].array.F[ii] = piaacmcparams.outtmp_array[ii];
                // square to give intensity
                double tmpv = piaacmcparams.outtmp_array[ii] *
                              piaacmcparams.outtmp_array[ii];
                // total intensity = sum(intensity_i) = sum(Re_i^2 + Im_i^2)
                value += tmpv;
            }

            /// - Total flux in the output vector is stored in \c piaacmcparams.PIAACMCSIMUL_VAL0 as total flux
            piaacmcparams.PIAACMCSIMUL_VAL0 = value;

            /// set \c value to average value per area normalized to flux
            value =
                value / size / size /
                piaacmcopticalsystem.flux
                    [0]; // flux[0] is proportional to the number of lambda channels, so this normalization makes value independant of number of spectral channels
            // here value is the total light (averaged across spectral channels) in the measurement points, normalized to the input flux
            // actual average contrast, averaging over # of pixels and physical area
            avContrast =
                value / (piaacmcparams.SCORINGTOTAL * focscale * focscale);

            //        printf("Peak constrast (rough estimate)= %g\n", peakcontrast/size/size/piaacmcopticalsystem.flux[0]/focscale/focscale*3.0);
            //        printf("value1 = %g\n", value1);
            //		printf("Total light in scoring field = %g  -> Average contrast = %g   (%g)\n", value, value/(arith_image_total("scoringmask")*focscale*focscale), value1/piaacmcparams.CnormFactor/piaacmcopticalsystem.nblambda);
        }
    }
    else /// ## Full/Slow PSF computation (if piaacmcparams.computePSF_FAST_FPMresp = 0)
    {
        /// The PSF for an extended source is approximated as
        /// a collection of point sources.
        /// Sourcesize determines the separation of the point sources
        // sourcesize > 0 only if in linear optimization (step >= 100)
        if (sourcesize != 0)
        {
            printf("COMPUTING RESOLVED SOURCE PSF / ADDING OPD MODES\n");
            fflush(stdout);

            double rad1 = 0.0;
            double rad2 = 0.0;

            // dld is the radius of the circle containing the point sources in lamba/D
            double dld =
                1.0 /
                pow(10.0, 0.1 * sourcesize); // nominal pointing offset [l/D]
            // extmode controls how many point sources we use to model the extended source
            // extmode == 0 => single point (unresolved source)
            // extmode == 1 => three point sources
            // extmode == 2 => six point sources
            if (extmode == 2)
            {
                // if I have six sources, put them in two rings of three each
                rad1 = dld / sqrt(2.5);
                rad2 = 2.0 * dld / sqrt(2.5);
            }
            else
            {
                // if I have three sources keep them at the same radii
                rad1 = dld;
                rad2 = dld;
            }

            // we will collect propagation results in a set of vectors called imvectp<1,2,...>
            // which will ultimately be collected into a single long imvect

            // image index, counts the number of PSFs we make, one for each point source
            long imindex = 0;

            // xld and yld are the input positions of the input source in lambda/D
            // initialize first point source, which sets optsyst
            if (extmode == 0)
            {
                FUNC_CHECK_RETURN(init_piaacmcopticalsystem(xld, yld));
            }
            else
            {
                FUNC_CHECK_RETURN(init_piaacmcopticalsystem(xld + rad1, yld));
            }

            makePIAAshapes();

            // propagate it (optsyst is a global), output in psfc0 (complex amlitude)
            // and psfi0 (intensity)
            {
                FUNC_CHECK_RETURN(OptSystProp_run(&piaacmcopticalsystem,
                                                  0,
                                                  startelem,
                                                  piaacmcopticalsystem.NBelem,
                                                  piaacmcparams.piaacmcconfdir,
                                                  0));
            }

            {
                // name of the image vector to hold the propagation result
                char imname[STRINGMAXLEN_IMGNAME];
                WRITE_IMAGENAME(imname, "imvectp%02ld", imindex);

                // convert the image to the vector
                FUNC_CHECK_RETURN(linopt_imtools_image_to_vec("psfc0",
                                                              "pixindex",
                                                              "pixmult",
                                                              imname,
                                                              NULL));
            }

            // save the intensity of the first point
            copy_image_ID("psfi0", "psfi0ext", 0);
            //sprintf(fname, "%s/psfi0_pt%03ld.fits", piaacmcparams.piaacmcconfdir, imindex);
            //save_fits("psfi0", fname);
            imindex++;

            if (extmode > 0)
            {
                // do the same for the second point

                double pha = 2.0 * M_PI / 3.0; // 1/3 around the circle

                {
                    FUNC_CHECK_RETURN(
                        init_piaacmcopticalsystem(xld + rad1 * cos(pha),
                                                  yld + rad1 * sin(pha)));
                }

                FUNC_CHECK_RETURN(makePIAAshapes());

                {
                    errno_t fret = OptSystProp_run(&piaacmcopticalsystem,
                                                   0,
                                                   startelem,
                                                   piaacmcopticalsystem.NBelem,
                                                   piaacmcparams.piaacmcconfdir,
                                                   0);
                    if (fret != RETURN_SUCCESS)
                    {
                        FUNC_RETURN_FAILURE("Call to OptSystProp_run failed");
                    }
                }

                {
                    char imname[STRINGMAXLEN_IMGNAME];
                    WRITE_IMAGENAME(imname, "imvectp%02ld", imindex);

                    FUNC_CHECK_RETURN(linopt_imtools_image_to_vec("psfc0",
                                                                  "pixindex",
                                                                  "pixmult",
                                                                  imname,
                                                                  NULL));
                }

                // add the intensity to build up PSF for extended source
                arith_image_add_inplace("psfi0ext", "psfi0");
                //sprintf(fname, "%s/psfi0_pt%03ld.fits", piaacmcparams.piaacmcconfdir, imindex);
                //save_fits("psfi0", fname);
                imindex++;

                DEBUG_TRACEPOINT(" ");
                // do the same for the third point

                pha = 4.0 * M_PI / 3.0; // 2/3 around the circle
                FUNC_CHECK_RETURN(
                    init_piaacmcopticalsystem(xld + rad1 * cos(pha),
                                              yld + rad1 * sin(pha)));

                FUNC_CHECK_RETURN(makePIAAshapes());

                {
                    FUNC_CHECK_RETURN(
                        OptSystProp_run(&piaacmcopticalsystem,
                                        0,
                                        startelem,
                                        piaacmcopticalsystem.NBelem,
                                        piaacmcparams.piaacmcconfdir,
                                        0));
                }

                {
                    char imname[STRINGMAXLEN_IMGNAME];
                    WRITE_IMAGENAME(imname, "imvectp%02ld", imindex);
                    FUNC_CHECK_RETURN(linopt_imtools_image_to_vec("psfc0",
                                                                  "pixindex",
                                                                  "pixmult",
                                                                  imname,
                                                                  NULL));
                }

                // add the intensity to build up PSF for extended source
                arith_image_add_inplace("psfi0ext", "psfi0");
                //sprintf(fname, "%s/psfi0_pt%03ld.fits", piaacmcparams.piaacmcconfdir, imindex);
                //save_fits("psfi0", fname);
                imindex++;
            }

            if (extmode == 2)
            {
                // keep going for the other three points if desired, on the outer radius
                double pha = M_PI / 3.0;
                FUNC_CHECK_RETURN(
                    init_piaacmcopticalsystem(xld + rad2 * cos(pha),
                                              yld + rad2 * sin(pha)));

                FUNC_CHECK_RETURN(makePIAAshapes());

                FUNC_CHECK_RETURN(OptSystProp_run(&piaacmcopticalsystem,
                                                  0,
                                                  startelem,
                                                  piaacmcopticalsystem.NBelem,
                                                  piaacmcparams.piaacmcconfdir,
                                                  0));

                {
                    char imname[STRINGMAXLEN_IMGNAME];
                    WRITE_IMAGENAME(imname, "imvectp%02ld", imindex);
                    FUNC_CHECK_RETURN(linopt_imtools_image_to_vec("psfc0",
                                                                  "pixindex",
                                                                  "pixmult",
                                                                  imname,
                                                                  NULL));
                }

                arith_image_add_inplace("psfi0ext", "psfi0");
                //sprintf(fname, "%s/psfi0_pt%03ld.fits", piaacmcparams.piaacmcconfdir, imindex);
                //save_fits("psfi0", fname);
                imindex++;

                pha = 2.0 * M_PI / 3.0 + M_PI / 3.0;
                FUNC_CHECK_RETURN(
                    init_piaacmcopticalsystem(xld + rad2 * cos(pha),
                                              yld + rad2 * sin(pha)));

                FUNC_CHECK_RETURN(makePIAAshapes());

                FUNC_CHECK_RETURN(OptSystProp_run(&piaacmcopticalsystem,
                                                  0,
                                                  startelem,
                                                  piaacmcopticalsystem.NBelem,
                                                  piaacmcparams.piaacmcconfdir,
                                                  0));

                {
                    char imname[STRINGMAXLEN_IMGNAME];
                    WRITE_IMAGENAME(imname, "imvectp%02ld", imindex);
                    FUNC_CHECK_RETURN(linopt_imtools_image_to_vec("psfc0",
                                                                  "pixindex",
                                                                  "pixmult",
                                                                  imname,
                                                                  NULL));
                }

                arith_image_add_inplace("psfi0ext", "psfi0");
                //sprintf(fname, "%s/psfi0_pt%03ld.fits", piaacmcparams.piaacmcconfdir, imindex);
                //save_fits("psfi0", fname);
                imindex++;

                pha = 4.0 * M_PI / 3.0 + M_PI / 3.0;

                FUNC_CHECK_RETURN(
                    init_piaacmcopticalsystem(xld + rad2 * cos(pha),
                                              yld + rad2 * sin(pha)));

                FUNC_CHECK_RETURN(makePIAAshapes());

                FUNC_CHECK_RETURN(OptSystProp_run(&piaacmcopticalsystem,
                                                  0,
                                                  startelem,
                                                  piaacmcopticalsystem.NBelem,
                                                  piaacmcparams.piaacmcconfdir,
                                                  0));

                {
                    char imname[STRINGMAXLEN_IMGNAME];
                    WRITE_IMAGENAME(imname, "imvectp%02ld", imindex);
                    FUNC_CHECK_RETURN(linopt_imtools_image_to_vec("psfc0",
                                                                  "pixindex",
                                                                  "pixmult",
                                                                  imname,
                                                                  NULL));
                }

                arith_image_add_inplace("psfi0ext", "psfi0");
                //sprintf(fname, "%s/psfi0_pt%03ld.fits", piaacmcparams.piaacmcconfdir, imindex);
                //save_fits("psfi0", fname);
                imindex++;

                // but multiply by 0.5 'cause we have twice as many points
                //    arith_image_cstmult_inplace("psfi0ext", 0.5);
            }

            /// ### OPTIONAL: Add OPD error to list of modes

            printf("Adding optional OPD error modes (%ld modes)\n", nbOPDerr);
            fflush(stdout);
            // add error modes if any
            // add new evaluation points for the error to imvect so we minimize the error
            // as well as the non-error
            for (long OPDmode = 0; OPDmode < nbOPDerr; OPDmode++)
            {
                imageID IDopderr;
                FUNC_CHECK_RETURN(
                    create_2Dimage_ID("opderr", size, size, &IDopderr));

                // "opderr" is a standard name read by PIAACMCsimul_init
                for (uint64_t ii = 0; ii < size * size; ii++)
                {
                    data.image[IDopderr].array.F[ii] =
                        data.image[IDopderrC]
                            .array.F[size * size * OPDmode + ii];
                }

                FUNC_CHECK_RETURN(
                    init_piaacmcopticalsystem(0.0,
                                              0.0)); // add error to the data

                FUNC_CHECK_RETURN(makePIAAshapes());

                FUNC_CHECK_RETURN(OptSystProp_run(&piaacmcopticalsystem,
                                                  0,
                                                  startelem,
                                                  piaacmcopticalsystem.NBelem,
                                                  piaacmcparams.piaacmcconfdir,
                                                  0));

                {
                    char imname[STRINGMAXLEN_IMGNAME];
                    WRITE_IMAGENAME(imname, "imvectp%02ld", imindex);
                    FUNC_CHECK_RETURN(linopt_imtools_image_to_vec("psfc0",
                                                                  "pixindex",
                                                                  "pixmult",
                                                                  imname,
                                                                  NULL));
                }

                arith_image_add_inplace("psfi0ext", "psfi0");
                //sprintf(fname, "%s/psfi0_pt%03ld.fits", piaacmcparams.piaacmcconfdir, imindex); //TEST
                //save_fits("psfi0", fname); //TEST
                delete_image_ID("opderr", DELETE_IMAGE_ERRMODE_WARNING);

                imindex++;
            }

            /// - Average over all the PSFs we've created to simulate this extended source
            long NBimindex = imindex;
            arith_image_cstmult_inplace("psfi0ext", 1.0 / NBimindex);

            /// - If \c outsave = 1, save PSF to FITS file
            if (outsave == 1)
            {
                PIAACMCsimul_update_fnamedescr();

                char fname[STRINGMAXLEN_FULLFILENAME];
                WRITE_FULLFILENAME(fname,
                                   "%s/psfi0_exsrc%03d_sm%d.%s.fits",
                                   piaacmcparams.piaacmcconfdir,
                                   sourcesize,
                                   piaacmcparams.SCORINGMASKTYPE,
                                   piaacmcparams.fnamedescr);

                save_fits("psfi0ext", fname);
            }

            // get the number of elements in a single-PSF vector
            {
                imageID ID     = image_ID("imvectp00");
                long    nbelem = data.image[ID].md[0].nelement;

                // make big vector to collect the complex amplitudes of all the above PSFs
                ID = image_ID("imvect");
                if (ID != -1)
                {
                    delete_image_ID("imvect", DELETE_IMAGE_ERRMODE_WARNING);
                }

                {
                    long offset =
                        nbelem /
                        piaacmcopticaldesign
                            .nblambda; // number of pixels per lambda x2 (re, im)
                    printf("offset = %ld\n", offset);

                    // number of pixels per lambda x 2 times the number of PSFs
                    long   offset1   = NBimindex * offset;
                    double normcoeff = 1.0 / sqrt(NBimindex);

                    // make an imvect for each lambda
                    imageID ID;
                    create_2Dimage_ID("imvect",
                                      offset1,
                                      piaacmcopticaldesign.nblambda,
                                      &ID);

                    // fill in with each imvectp* created above
                    for (imindex = 0; imindex < NBimindex; imindex++)
                    {
                        char imname[STRINGMAXLEN_IMGNAME];
                        WRITE_IMAGENAME(imname, "imvectp%02ld", imindex);
                        imageID ID1 = image_ID(imname);
                        for (long kl = 0; kl < piaacmcopticaldesign.nblambda;
                             kl++)
                            for (long ii = 0; ii < offset; ii++)
                            {
                                data.image[ID].array.F[kl * offset1 +
                                                       imindex * offset + ii] =
                                    data.image[ID1].array.F[kl * offset + ii] *
                                    normcoeff;
                            }
                        delete_image_ID(imname, DELETE_IMAGE_ERRMODE_ERROR);
                    }
                }
            }

            //linopt_imtools_image_to_vec("psfc0", "pixindex", "pixmult", "imvect");

            // measure the contrast for all aimplitudes in imvect
            double value = 0.0;

            double peakcontrast = 0.0;
            {
                imageID ID = image_ID("imvect");
                for (uint64_t ii = 0; ii < data.image[ID].md[0].nelement; ii++)
                {
                    double tmpv =
                        data.image[ID].array.F[ii] * data.image[ID].array.F[ii];
                    value += tmpv;
                    if (tmpv > peakcontrast)
                    {
                        peakcontrast = tmpv;
                    }
                }
            }

            for (long elem = 0; elem < piaacmcopticalsystem.NBelem; elem++)
            {
                printf("    FLUX %3ld   %12.4lf %8.6lf\n",
                       elem,
                       piaacmcopticalsystem.flux[elem],
                       piaacmcopticalsystem.flux[elem] /
                           piaacmcopticalsystem.flux[0]);
            }

            /// - If \c outsave = 1, save flux to txt file
            if (outsave == 1)
            {
                PIAACMCsimul_update_fnamedescr();

                char fname[STRINGMAXLEN_FULLFILENAME];
                WRITE_FULLFILENAME(fname,
                                   "%s/flux_exsrc%03d_sm%d.%s.txt",
                                   piaacmcparams.piaacmcconfdir,
                                   sourcesize,
                                   piaacmcparams.SCORINGMASKTYPE,
                                   piaacmcparams.fnamedescr);

                FILE *fpflux = fopen(fname, "w");
                for (long elem = 0; elem < piaacmcopticalsystem.NBelem; elem++)
                {
                    fprintf(fpflux,
                            "%18.16lf %18.16lf  %d\n",
                            piaacmcopticalsystem.flux[elem],
                            piaacmcopticalsystem.flux[elem] /
                                piaacmcopticalsystem.flux[0],
                            piaacmcopticalsystem.nblambda);
                }
                fprintf(fpflux, "W0  %d\n", piaacmcopticalsystem.nblambda);
                fclose(fpflux);

                EXECUTE_SYSTEM_COMMAND("cp %s %s/tmp_flux.txt",
                                       fname,
                                       piaacmcparams.piaacmcconfdir);
            }

            // compute average contrast
            value = value / size / size / piaacmcopticalsystem.flux[0];
            avContrast =
                value / (piaacmcparams.SCORINGTOTAL * focscale * focscale);

            //           piaacmcparams.CnormFactor = size*size*piaacmcopticalsystem.flux[0]*arith_image_total("scoringmask")*focscale*focscale; // /piaacmcopticalsystem.nblambda;
            piaacmcparams.CnormFactor = focscale * focscale * size * size *
                                        piaacmcopticalsystem.flux[0] /
                                        piaacmcopticalsystem.nblambda;

            if (piaacmcparams.WRITE_OK == 1)
            {
                char fname[STRINGMAXLEN_FULLFILENAME];
                WRITE_FULLFILENAME(fname,
                                   "%s/CnormFactor.txt",
                                   piaacmcparams.piaacmcconfdir);

                FILE *fp = fopen(fname, "w");
                fprintf(fp,
                        "# Written by function PIAACMCsimul_computePSF()\n");
                fprintf(fp,
                        "# Cnormfactor = "
                        "focscale*focscale*size*size*piaacmcopticalsystem.flux["
                        "0]/piaacmcopticalsystem.nblambda\n");
                fprintf(fp, "# 0    focscale  size  flux[0]  nblambda\n");
                fprintf(fp, "\n");

                fprintf(fp, "%g\n", piaacmcparams.CnormFactor);
                fprintf(fp,
                        "0      %g %ul %g %d\n",
                        focscale,
                        size,
                        piaacmcopticalsystem.flux[0],
                        piaacmcopticalsystem.nblambda);
                fclose(fp);
            }

            printf("COMPUTING RESOLVED SOURCE PSF\n");
            printf("SCORINGTOTAL = %f  %f\n",
                   piaacmcparams.SCORINGTOTAL,
                   arith_image_total("scoringmask"));
            printf("Peak constrast (rough estimate)= %g\n",
                   peakcontrast / size / size / piaacmcopticalsystem.flux[0] /
                       focscale / focscale / normcoeff / normcoeff);

            avContrast = value / (arith_image_total("scoringmask") * focscale *
                                  focscale);
            printf(
                "Total light in scoring field = %g, peak PSF = %g   -> Average "
                "contrast = %g\n",
                value,
                piaacmcopticaldesign.peakPSF,
                avContrast);

            if (outsave == 1)
            {
                PIAACMCsimul_update_fnamedescr();

                char fname[STRINGMAXLEN_FULLFILENAME];
                WRITE_FULLFILENAME(fname,
                                   "%s/contrast_exsrc%03d_sm%d.%s.txt",
                                   piaacmcparams.piaacmcconfdir,
                                   sourcesize,
                                   piaacmcparams.SCORINGMASKTYPE,
                                   piaacmcparams.fnamedescr);

                FILE *fp = fopen(fname, "w");
                fprintf(fp, "%g", avContrast);
                fclose(fp);
            }
        }
        else // called for step 0 through 15.  Does not use OPDerr
        {
            // compute the PSF for a single point source at offset xld, yld
            printf("COMPUTING UNRESOLVED SOURCE PSF [%f x %f]\n", xld, yld);

            // ========== initializes optical system to piaacmc design ===========
            // xld and yld are the input positions of the input source in lamba/D
            // initialize first point source, which sets optsyst

            /// calls PIAACMCsimul_init()
            FUNC_CHECK_RETURN(init_piaacmcopticalsystem(xld, yld));

            FUNC_CHECK_RETURN(makePIAAshapes());

            // ============ perform propagations ================
            // propagate it (optsyst is a global), output in psfc0 (complex amlitude)
            // and psfi0 (intensity)

            FUNC_CHECK_RETURN(OptSystProp_run(&piaacmcopticalsystem,
                                              0,
                                              startelem,
                                              piaacmcopticalsystem.NBelem,
                                              piaacmcparams.piaacmcconfdir,
                                              0));

            if (outsave == 1)
            {
                PIAACMCsimul_update_fnamedescr();
                char fname[STRINGMAXLEN_FULLFILENAME];
                WRITE_FULLFILENAME(fname,
                                   "%s/psfi0_ptsrc_sm%d.%s.fits",
                                   piaacmcparams.piaacmcconfdir,
                                   piaacmcparams.SCORINGMASKTYPE,
                                   piaacmcparams.fnamedescr);

                FUNC_CHECK_RETURN(save_fits("psfi0", fname));
            }

            // linearize the result into imvect
            FUNC_CHECK_RETURN(linopt_imtools_image_to_vec("psfc0",
                                                          "pixindex",
                                                          "pixmult",
                                                          "imvect",
                                                          NULL));

            FUNC_CHECK_RETURN(
                save_fits("imvect", "./testdir/test_imvect.fits"));

            // extract amplitude and phase for diagnostics
            //mk_amph_from_complex("psfc0", "psfc0a", "psfc0p", 0);
            //save_fits("psfc0a", "test_psfc0a.fits");
            //delete_image_ID("psfc0a");
            //delete_image_ID("psfc0p");
            //printf("saved -> test_psfc0a.fits\n");
            //fflush(stdout);

            // if(savepsf==2)
            //	{
            //	mk_reim_from_complex("psfc0", "psfc0re", "psfc0im", 0);
            /*					save_fits("psfi0", "test_psfi0.fits");
            					save_fits("psfc0re", "test_psfc0re.fits");
            					save_fits("psfc0im", "test_psfc0im.fits");
            					sleep(100000);*/

            //	}

            {
                DEBUG_TRACEPOINT("compute average contrast");
                double value = 0.0;
                peakcontrast = 0.0;
                {
                    imageID ID = image_ID("imvect");
                    for (uint64_t ii = 0; ii < data.image[ID].md[0].nelement;
                         ii += 2)
                    {
                        // intensity as Re^2 + Im^2
                        double tmpv = data.image[ID].array.F[ii] *
                                          data.image[ID].array.F[ii] +
                                      data.image[ID].array.F[ii + 1] *
                                          data.image[ID].array.F[ii + 1];
                        value += tmpv;
                        if (tmpv > peakcontrast)
                        {
                            peakcontrast = tmpv;
                        }
                    }
                }

                DEBUG_TRACEPOINT("report the contrast");
                for (long elem = 0; elem < piaacmcopticalsystem.NBelem; elem++)
                {
                    printf("    FLUX %3ld   %12.4lf %8.6lf\n",
                           elem,
                           piaacmcopticalsystem.flux[elem],
                           piaacmcopticalsystem.flux[elem] /
                               piaacmcopticalsystem.flux[0]);
                }
                //            value = value/size/size/piaacmcopticalsystem.flux[0];

                if (outsave == 1)
                {
                    DEBUG_TRACEPOINT("writing flux file");

                    PIAACMCsimul_update_fnamedescr();

                    char fname[STRINGMAXLEN_FULLFILENAME];
                    WRITE_FULLFILENAME(fname,
                                       "%s/flux_ptsrc_sm%d.%s.txt",
                                       piaacmcparams.piaacmcconfdir,
                                       piaacmcparams.SCORINGMASKTYPE,
                                       piaacmcparams.fnamedescr);

                    FILE *fpflux = fopen(fname, "w");
                    if (fpflux == NULL)
                    {
                        PRINT_ERROR("cannot create file %s", fname);
                        abort();
                    }
                    for (long elem = 0; elem < piaacmcopticalsystem.NBelem;
                         elem++)
                    {
                        fprintf(fpflux,
                                "%18.16lf %18.16lf  %d\n",
                                piaacmcopticalsystem.flux[elem],
                                piaacmcopticalsystem.flux[elem] /
                                    piaacmcopticalsystem.flux[0],
                                piaacmcopticalsystem.nblambda);
                    }
                    fprintf(fpflux, "W1\n");
                    fclose(fpflux);

                    EXECUTE_SYSTEM_COMMAND("cp %s %s/tmp_flux.txt",
                                           fname,
                                           piaacmcparams.piaacmcconfdir);
                }

                //         piaacmcparams.CnormFactor = size*size*piaacmcopticalsystem.flux[0]*arith_image_total("scoringmask")*focscale*focscale; // /piaacmcopticalsystem.nblambda;
                piaacmcparams.CnormFactor = focscale * focscale * size * size *
                                            piaacmcopticalsystem.flux[0] /
                                            piaacmcopticalsystem.nblambda;

                {
                    DEBUG_TRACEPOINT("writing CnormFactor file");

                    char fname[STRINGMAXLEN_FULLFILENAME];
                    WRITE_FULLFILENAME(fname,
                                       "%s/CnormFactor.txt",
                                       piaacmcparams.piaacmcconfdir);

                    FILE *fp = fopen(fname, "w");
                    if (fp == NULL)
                    {
                        PRINT_ERROR("cannot create file %s", fname);
                        abort();
                    }
                    fprintf(
                        fp,
                        "# Written by function PIAACMCsimul_computePSF()\n");
                    fprintf(fp,
                            "# Cnormfactor = "
                            "focscale*focscale*size*size*piaacmcopticalsystem."
                            "flux[0]/piaacmcopticalsystem.nblambda\n");
                    fprintf(fp, "# 1    focscale  size  flux[0]  nblambda\n");
                    fprintf(fp, "\n");
                    fprintf(fp, "%g\n", piaacmcparams.CnormFactor);
                    fprintf(fp,
                            "1      %g %ul %g %d\n",
                            focscale,
                            size,
                            piaacmcopticalsystem.flux[0],
                            piaacmcopticalsystem.nblambda);
                    fclose(fp);
                }
                // here we're essentially done!

                printf("COMPUTING UNRESOLVED SOURCE PSF -*- [%f x %f]\n",
                       xld,
                       yld);
                printf("SCORINGTOTAL = %f  %f\n",
                       piaacmcparams.SCORINGTOTAL,
                       arith_image_total("scoringmask"));

                if ((variable_ID("PIAACMC_NOFPM")) != -1)
                {
                    imageID ID                   = image_ID("psfc0");
                    piaacmcopticaldesign.peakPSF = 0.0;
                    for (uint64_t ii = 0; ii < size * size; ii++)
                    {
                        double val = data.image[ID].array.CF[ii].re *
                                         data.image[ID].array.CF[ii].re +
                                     data.image[ID].array.CF[ii].im *
                                         data.image[ID].array.CF[ii].im;
                        if (val > piaacmcopticaldesign.peakPSF)
                        {
                            piaacmcopticaldesign.peakPSF = val;
                        }
                    }

                    FILE *fp = fopen("conf/conf_peakPSF.txt", "w");
                    fprintf(fp, "%g\n", piaacmcopticaldesign.peakPSF);
                    fclose(fp);
                }

                if (piaacmcopticaldesign.peakPSF < 1.0)
                {
                    printf("Peak constrast (rough estimate)= %g -> %g\n",
                           peakcontrast,
                           peakcontrast / (piaacmcopticalsystem.flux[0] *
                                           piaacmcopticalsystem.flux[0]));
                    printf("piaacmcopticalsystem.flux[0]  = %g\n",
                           piaacmcopticalsystem.flux[0]);
                    printf("SCORINGMASKTYPE = %d\n",
                           piaacmcparams.SCORINGMASKTYPE);
                    avContrast = value /
                                 (piaacmcopticalsystem.flux[0] *
                                  piaacmcopticalsystem.flux[0]) /
                                 piaacmcparams.SCORINGTOTAL;
                    printf(
                        "[0] Total light in scoring field = %g, peak PSF = %g, "
                        "SCOTINGTOTAL = %g   -> Average "
                        "contrast = %g\n",
                        value,
                        piaacmcopticaldesign.peakPSF,
                        piaacmcparams.SCORINGTOTAL,
                        avContrast);
                }
                else
                {
                    printf("Peak constrast = %g -> %g\n",
                           peakcontrast,
                           peakcontrast / piaacmcopticaldesign.peakPSF);
                    printf("SCORINGMASKTYPE = %d\n",
                           piaacmcparams.SCORINGMASKTYPE);
                    avContrast = value / piaacmcopticaldesign.peakPSF /
                                 piaacmcparams.SCORINGTOTAL;
                    printf(
                        "[1] Total light in scoring field = %g, peak PSF = %g, "
                        "SCOTINGTOTAL = %g  -> Average "
                        "contrast = %g\n",
                        value,
                        piaacmcopticaldesign.peakPSF,
                        piaacmcparams.SCORINGTOTAL,
                        avContrast);

                    DEBUG_TRACEPOINT("writing contrast value to file");
                    FILE *fp;
                    if ((fp = fopen("PSFcontrastval.txt", "w")) != NULL)
                    {
                        fprintf(fp,
                                "%g %g\n",
                                peakcontrast / piaacmcopticaldesign.peakPSF,
                                value / piaacmcopticaldesign.peakPSF /
                                    piaacmcparams.SCORINGTOTAL);
                        fclose(fp);
                    }
                }

                DEBUG_TRACEPOINT("outsave = %d", outsave);
                if (outsave == 1)
                {
                    PIAACMCsimul_update_fnamedescr();
                    char fname[STRINGMAXLEN_FULLFILENAME];
                    WRITE_FULLFILENAME(fname,
                                       "%s/contrast_ptsrc_sm%d.%s.txt",
                                       piaacmcparams.piaacmcconfdir,
                                       piaacmcparams.SCORINGMASKTYPE,
                                       piaacmcparams.fnamedescr);

                    printf("saving contrast value [%g] -> %s\n",
                           avContrast,
                           fname);

                    FILE *fp = fopen(fname, "w");
                    fprintf(fp, "%g", avContrast);
                    fclose(fp);
                }
            }
        }
    }

    if (contrastval != NULL)
    {
        *contrastval = avContrast;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
