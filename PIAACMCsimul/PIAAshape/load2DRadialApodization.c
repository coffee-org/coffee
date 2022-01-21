/**
 * @file    PIAACMCsimul_load2DRadialApodization.c
 * @brief   PIAA-type coronagraph design
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>

// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "info/info.h"
#include "linopt_imtools/linopt_imtools.h"

#include "PIAACMCsimul.h"

/**
 * @brief Load and fit radial apodization profile
 *
 * Purpose
 * -------
 *
 * Input is 2D apodization map image.
 * Output is fit.
 * modal basis is mk(r) : cos(r*k*M_PI/1.3)
 *
 * Arguments
 * -------
 * @param[in] IDapo_name
 * @param[in] beamradpix
 * @param[put] IDapofit_name
 * @return errno_t
 */
errno_t load2DRadialApodization(const char *__restrict__ IDapo_name,
                                float beamradpix,
                                const char *__restrict__ IDapofit_name)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %s %f %s", IDapo_name, beamradpix, IDapofit_name);

    long  kmax  = 10;
    float eps   = 1.0e-4;
    int   debug = 1;

    uint32_t sizem = (long) (beamradpix * 2);

    { // Create 2D radial cos modes
        if (image_ID("APOmodesCos") == -1)
        {
            char fname[STRINGMAXLEN_FULLFILENAME];
            FUNC_CHECK_RETURN(
                linopt_imtools_makeCosRadModes("APOmodesCos",
                                               sizem,
                                               kmax,
                                               ApoFitCosFact * beamradpix,
                                               1.0,
                                               NULL));

            WRITE_FULLFILENAME(fname,
                               "%s/APOmodesCos.fits",
                               piaacmcparams.piaacmcconfdir);
            FUNC_CHECK_RETURN(save_fits("APOmodesCos", fname));
        }
    }

    { // CREATE MASK AND CROP INPUT
        imageID IDmask;
        FUNC_CHECK_RETURN(
            create_2Dimage_ID("fitmaskapo", sizem, sizem, &IDmask));

        imageID  IDin   = image_ID(IDapo_name);
        uint32_t sizein = data.image[IDin].md[0].size[0];

        imageID ID;
        FUNC_CHECK_RETURN(create_2Dimage_ID("_apoincrop", sizem, sizem, &ID));

        long offset = (sizein - sizem) / 2;
        for (uint32_t ii = 0; ii < sizem; ii++)
            for (uint32_t jj = 0; jj < sizem; jj++)
            {
                data.image[ID].array.F[jj * sizem + ii] =
                    data.image[IDin]
                        .array.F[(jj + offset) * sizein + (ii + offset)];
                if ((data.image[ID].array.F[jj * sizem + ii] > eps) &&
                    (ii % 1 == 0) && (jj % 1 == 0))
                {
                    data.image[IDmask].array.F[jj * sizem + ii] = 1.0;
                }
            }
    }

    { // for debugging, write files to filesystem
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname,
                           "%s/_apoincrop.fits",
                           piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("_apoincrop", fname));

        WRITE_FULLFILENAME(fname,
                           "%s/fitmaskapo.fits",
                           piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("fitmaskapo", fname));
    }

    // Perform linear fit
    FUNC_CHECK_RETURN(linopt_imtools_image_fitModes("_apoincrop",
                                                    "APOmodesCos",
                                                    "fitmaskapo",
                                                    1.0e-8,
                                                    IDapofit_name,
                                                    0,
                                                    NULL));

    EXECUTE_SYSTEM_COMMAND("mv %s/eigenv.dat %s/eigenv_APOmodesCos.dat",
                           piaacmcparams.piaacmcconfdir,
                           piaacmcparams.piaacmcconfdir);

    if (debug == 1) // test fit quality
    {
        char fname[STRINGMAXLEN_FULLFILENAME];

        FUNC_CHECK_RETURN(linopt_imtools_image_construct("APOmodesCos",
                                                         IDapofit_name,
                                                         "testapofitsol",
                                                         NULL));

        WRITE_FULLFILENAME(fname,
                           "%s/testapofitsol.fits",
                           piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("testapofitsol", fname));

        arith_image_sub("_apoincrop", "testapofitsol", "apofitres");
        arith_image_mult("apofitres", "fitmaskapo", "apofitresm");

        WRITE_FULLFILENAME(fname,
                           "%s/_2Dapofitres.fits",
                           piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("apofitres", fname));

        WRITE_FULLFILENAME(fname,
                           "%s/_2Dapofitresm.fits",
                           piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("apofitresm", fname));

        // linopt_imtools_image_fitModes("apofitres", "APOmodesCos", "fitmaskapo", 1.0e-5, "test2c", 0);
        info_image_stats("apofitresm", "");
    }

    FUNC_CHECK_RETURN(
        delete_image_ID("_apoincrop", DELETE_IMAGE_ERRMODE_WARNING));

    FUNC_CHECK_RETURN(
        delete_image_ID("fitmaskapo", DELETE_IMAGE_ERRMODE_WARNING));

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
