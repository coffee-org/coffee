/**
 * @file    PIAAACMCsimul_exec_optimize_lyot_stops_shapes_positions.c
 * @brief   PIAA-type coronagraph design, execute compute image
 *
 */

// log all debug trace points to file
#define DEBUGLOG

// System includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



// milk includes
#include "CommandLineInterface/CLIcore.h"


#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"

#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"

#include "PIAACMCsimul_CA2propCubeInt.h"
#include "PIAACMCsimul_computePSF.h"
#include "init_piaacmcopticalsystem.h"
#include "init_piaacmcopticaldesign.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"

#include "PIAAshape/makePIAAshapes.h"
#include "LyotStop/optimizeLyotStop.h"







/** @brief Make Lyot stops using off-axis light minimums
 *
 *  Finds minimum flux level in 3D intensity data cube
 *
 * @param[in]  IDincohc_name             image: Input 3D intensity
 * @param[out] test_oals_val.fits        FITS : Minimum 2D image
 * @param[out] test_oals_index.fits      FITS : z-index for each pixel of the 2D output minimum
 *
*/
errno_t optimizeLyotStop_offaxis_min(
    const char *__restrict__ IDincohc_name
)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %s", IDincohc_name);


    imageID IDincohc = image_ID(IDincohc_name);
    if(IDincohc == -1)
    {
        FUNC_RETURN_FAILURE("Image %s not found in memory", IDincohc_name);
    }

    uint32_t xsize = data.image[IDincohc].md[0].size[0];
    uint32_t ysize = data.image[IDincohc].md[0].size[1];
    uint64_t xysize = xsize;
    xysize *= ysize;

    uint32_t NBz = data.image[IDincohc].md[0].size[2];

    imageID IDindex;
    FUNC_CHECK_RETURN(
        create_2Dimage_ID("oals_index", xsize, ysize, &IDindex)
    );

    imageID IDminflux;
    FUNC_CHECK_RETURN(
        create_2Dimage_ID("oals_val", xsize, ysize, &IDminflux)
    );

    DEBUG_TRACEPOINT("scanning for minimum");
    for(uint64_t ii = 0; ii < xysize; ii++)
    {
        float minv = data.image[IDincohc].array.F[ii];
        uint32_t minindex = 0;

        for(uint32_t kk=1; kk < NBz; kk++)
        {
            float tmpv = data.image[IDincohc].array.F[xysize*kk + ii];
            if(tmpv < minv)
            {
                minv = tmpv;
                minindex = kk;
            }
        }
        data.image[IDminflux].array.F[ii] = minv;
        data.image[IDindex].array.F[ii] = (float) minindex;
    }

    DEBUG_TRACEPOINT("Saving minimum map to filesystem");

    FUNC_CHECK_RETURN(
        save_fits("oals_index", "test_oals_index.fits")
    );

    FUNC_CHECK_RETURN(
        save_fits("oals_val", "test_oals_val.fits")
    );

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}








/**
  * ---
  *
  * ## Mode 5: Optimize Lyot stops shapes and positions
  *
  *
  *
  */
errno_t exec_optimize_lyot_stops_shapes_positions()
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG");

    imageID IDv;
    double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;


    imageID IDa, IDc;

    printf("=================================== mode 005 ===================================\n");
    // load some more cli variables
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

    /// ### Initialize as in mode 0
    {
        uint64_t initflag = INIT_PIAACMCOPTICALDESIGN_MODE__READCONF;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF;
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

    FUNC_CHECK_RETURN(
        makePIAAshapes()
    );

    piaacmcopticalsystem.FOCMASKarray[0].mode = 1; // use 1-fpm



    /// ### Load CLI variables as appropriate

    /// - <- **PIAACMC_nbpropstep** : # of propagation steps along the beam
    long NBpropstep = 150;
    if((IDv = variable_ID("PIAACMC_nbpropstep")) != -1)
    {
        NBpropstep = (long) data.variable[IDv].value.f + 0.01;
    }

    /// - <- **PIAACMC_lstransm** : desired Lyot stop transmission
    double lstransm = 0.85;
    if((IDv = variable_ID("PIAACMC_lstransm")) != -1)
    {
        lstransm = (double) data.variable[IDv].value.f;
    }
    printf("lstransm  = %f\n", lstransm);



    /// ### Identify post focal plane pupil plane (first pupil after focal plane mask)

    /// Provides reference complex amplitude plane for downstream analysis
    FUNC_CHECK_RETURN(init_piaacmcopticalsystem(0.0, 0.0));

    printf("=========== %ld elements ======================================================\n",
           piaacmcopticalsystem.NBelem);
    // find the ID of the "post focal plane mask pupil" element
    long elem0 = 0;
    for(long elem = 0; elem < piaacmcopticalsystem.NBelem; elem++)
    {
        if(strcmp("post focal plane mask pupil", piaacmcopticalsystem.name[elem]) == 0)
        {
            elem0 = elem;
            DEBUG_TRACEPOINT("post focal plane mask pupil = %ld", elem);
        }
        else
        {
            DEBUG_TRACEPOINT("elem %ld : %s", elem, piaacmcopticalsystem.name[elem]);
        }
    }
    piaacmcopticalsystem.keepMem[elem0] = 1; // keep it for future use





    /// ### Compute incoherent 3D illumination near pupil plane

    /// Multiple off-axis sources are propagated and the corresponding intensities added
    /// - -> **OAincohc**  Output incoherent image (3D)

    DEBUG_TRACEPOINT("compute the reference on-axis PSF");
    {
        double cval = 0.0;
        FUNC_CHECK_RETURN(
            PIAACMCsimul_computePSF(
                0.0, 0.0, 0, piaacmcopticalsystem.NBelem, 0, 0, 0, 0, &cval)
        );
    }

    // filenames of the complex amplitude and phase in the post FPM pupil plane indexed by elem0
    char imnamea[STRINGMAXLEN_IMGNAME];
    WRITE_IMAGENAME(imnamea, "WFamp0_%03ld", elem0);

    char imnamep[STRINGMAXLEN_IMGNAME];
    WRITE_IMAGENAME(imnamep, "WFpha0_%03ld", elem0);


    // args for the PIAACMCsimul_CA2propCubeInt function

    // set range of propagation
    float zmin = piaacmcopticaldesign.LyotZmin;
    float zmax = piaacmcopticaldesign.LyotZmax;



    DEBUG_TRACEPOINT("propagate complex amplitude in a range from %f to %f", zmin, zmax);
    // computes the diffracted light from the on-axis source
    FUNC_CHECK_RETURN(
        PIAACMCsimul_CA2propCubeInt(imnamea, imnamep, zmin, zmax, NBpropstep,
                                    "iproptmp", NULL)
    );
    // complex amplitude at elem0, only used to determine image size
    IDa = image_ID(imnamea);

    uint32_t xsize = data.image[IDa].md[0].size[0];
    uint32_t ysize = data.image[IDa].md[0].size[1];
    uint64_t xysize = xsize;
    xysize *= ysize;

    /// OAincohc is the summed light "all" from off-axis sources in the pupil,
    /// including the on-axis source(!),
    /// giving intensity contribution of all off-axis sources
    /// in order to preserve the intensity of the off-axis in the design.
    /// load OAincohc if exist, maybe we've been here before
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/OAincohc.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(
            load_fits(fname, "OAincohc", LOADFITS_ERRMODE_WARNING, &IDc)
        );
    }

    DEBUG_TRACEPOINT("IDc = %ld", IDc);

    if(IDc == -1)
    {   // OAincohc does not exist so we have to make it

        DEBUG_TRACEPOINT("Creating off-axis light 3D volumetric cube");
        // create image to receive sum
        FUNC_CHECK_RETURN(
            create_3Dimage_ID("OAincohc", xsize, ysize, NBpropstep, &IDc)
        );

        long cnt = 0; // initialize counter so later we can normalize by number of sources


        // number of circle radii
        long NBkr = 2;  //5

        // loop over radii
        for(long kr = 0; kr < NBkr; kr++)
        {
            // number of off-axis sources on each circle
            long NBincpt = 9; //15;

            // loop over points at current radius
            for(long k1 = 0; k1 < NBincpt; k1++)
            {
                DEBUG_TRACEPOINT("OAincohc point %ld/%ld %ld/%ld", kr, NBkr, k1, NBincpt);
                // compute PSF for a point at this angle with scaled offset
                // PIAACMCsimul_computePSF changes fnamea and fnamep (in call to OptSystProp_run)!
                {
                    double cval = 0.0;
                    double oaoffset = 20.0; // off axis amplitude

                    FUNC_CHECK_RETURN(
                        PIAACMCsimul_computePSF(
                            oaoffset * (1.0 + kr) / NBkr * cos(2.0 * M_PI * k1 / NBincpt),
                            oaoffset * (1.0 + kr) / NBkr * sin(2.0 * M_PI * k1 / NBincpt),
                            0,
                            piaacmcopticalsystem.NBelem,
                            0,
                            0,
                            0,
                            0,
                            &cval
                        )
                    );
                }

                // propagate that elem0 from zmin to zmax with new PSF
                imageID ID1;
                FUNC_CHECK_RETURN(
                    PIAACMCsimul_CA2propCubeInt(imnamea, imnamep, zmin, zmax, NBpropstep,
                                                "iproptmp", &ID1)
                );
                for(uint64_t ii = 0; ii < xysize; ii++)
                {   // ii is indexing x-y plane
                    for(long k = 0; k < NBpropstep; k++)
                    {   // k is indexing z-direction. adding to IDc
                        data.image[IDc].array.F[xysize*k + ii] +=
                            data.image[ID1].array.F[xysize*k + ii];
                    }
                }
                FUNC_CHECK_RETURN(
                    delete_image_ID("iproptmp", DELETE_IMAGE_ERRMODE_WARNING)
                );
                cnt ++;
            }
        }
        // scale by the number of sources to give average
        for(uint64_t ii = 0; ii < xysize; ii++)
        {
            for(long k = 0; k < NBpropstep; k++)
            {
                data.image[IDc].array.F[k * xsize * ysize + ii] /= cnt;
            }
        }


        {   // save the final result
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(fname, "%s/OAincohc.fits", piaacmcparams.piaacmcconfdir);
            FUNC_CHECK_RETURN(save_fits("OAincohc", fname));
        }
    }



    DEBUG_TRACEPOINT("Compute on-axis PSF 3D intensity to define light to reject");
    FUNC_CHECK_RETURN(
        PIAACMCsimul_computePSF(
            0.0, 0.0, 0,
            piaacmcopticalsystem.NBelem,
            0, 0, 0, 0, NULL)
    );
    // propagate it into the optical system, with result in image named "iprop00"

    DEBUG_TRACEPOINT("compute on-axis 3D intensity cube");
    FUNC_CHECK_RETURN(
        PIAACMCsimul_CA2propCubeInt(
            imnamea,
            imnamep,
            zmin,
            zmax,
            NBpropstep,
            "iprop00",
            NULL
        )
    );
    //  save_fits("iprop00", "test_iprop00.fits");

    /// ### Compute image that has the min along z of OAincohc at each x,y
    /// Function optimizeLyotStop_offaxis_min() computes minimal intensity image.
    ///
    FUNC_CHECK_RETURN(
        optimizeLyotStop_offaxis_min(
            "OAincohc"
        )
    );

    /// ### Optimize Lyot stops

    // The actual Lyot stop shape and location optimization is done by PIAACMCsimul_optimizeLyotStop(),
    // producing optimal Lyot stops in optLM*.fits
    // and position relative to elem0 in piaacmcopticaldesign.LyotStop_zpos
    DEBUG_TRACEPOINT("Optimize Lyot stop geometry");
    FUNC_CHECK_RETURN(
        optimizeLyotStop(
            imnamea,
            imnamep,
            "OAincohc",
            zmin,
            zmax,
            lstransm,
            NBpropstep,
            piaacmcopticaldesign.NBLyotStop,
            NULL
        )
    );

    {
        char fptestname[STRINGMAXLEN_FILENAME];
        WRITE_FILENAME(fptestname, "conj_test.txt");

        DEBUG_TRACEPOINT("write %s", fptestname);

        FILE * fptest = fopen(fptestname, "w");
        if(fptest == NULL)
        {
            PRINT_ERROR("cannot create file %s", fptestname);
            abort();
        }
        fprintf(fptest, "# Optimal Lyot stop conjugations\n");
        fprintf(fptest, "# \n");
        fprintf(fptest, "# DO NOT EDIT THIS FILE\n");
        fprintf(fptest, "# Written by %s in %s\n", __FUNCTION__, __FILE__);
        fprintf(fptest, "# \n");
        fprintf(fptest,
                "# Lyot stop index   zmin   zmax    LyotStop_zpos    elemZpos[elem0]\n");
        for(long ls = 0; ls < piaacmcopticaldesign.NBLyotStop; ls++)
        {
            fprintf(fptest, "%5ld  %f  %f     %f  %f\n", ls, zmin, zmax,
                    piaacmcopticaldesign.LyotStop_zpos[ls], piaacmcopticalsystem.elemZpos[elem0]);
        }
        fclose(fptest);
    }

    // convert Lyot stop position from relative to elem0 to absolute
    for(long ls = 0; ls < piaacmcopticaldesign.NBLyotStop; ls++)
    {
        piaacmcopticaldesign.LyotStop_zpos[ls] += piaacmcopticalsystem.elemZpos[elem0];
    }

    // and we're done!  save.
    FUNC_CHECK_RETURN(
        PIAACMCsimul_savepiaacmcconf(piaacmcparams.piaacmcconfdir)
    );

    // copy to the final Lyot stop file for this mode
    for(long ls = 0; ls < piaacmcopticaldesign.NBLyotStop; ls++)
    {
        EXECUTE_SYSTEM_COMMAND("cp ./%s/optLM%02ld.fits ./%s/LyotStop%ld.fits",
                               piaacmcparams.piaacmcconfdir, ls, piaacmcparams.piaacmcconfdir, ls);
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}





