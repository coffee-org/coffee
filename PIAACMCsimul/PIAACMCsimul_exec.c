/**
 * @file    PIAAACMCsimul_init.c
 * @brief   PIAA-type coronagraph design, execute
 *
 * Can design both APLCMC and PIAACMC coronagraphs
 *
 */

// System includes
#include <assert.h>
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_tools/COREMOD_tools.h"

#include "image_gen/image_gen.h"
#include "statistic/statistic.h"

#include "OptSystProp/OptSystProp.h"
#include "OpticsMaterials/OpticsMaterials.h"
#include "linopt_imtools/linopt_imtools.h"

#include "PIAACMCsimul.h"

#include "PIAACMCsimul_achromFPMsol_eval_zonezderivative.h"
#include "PIAACMCsimul_computePSF.h"
#include "PIAACMCsimul_eval_poly_design.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"
#include "PIAACMCsimul_measure_transm_curve.h"
#include "exec_computePSF_no_fpm.h"
#include "exec_compute_image.h"
#include "exec_optimize_PIAA_shapes_fpmtransm.h"

#include "FocalPlaneMask/exec_multizone_fpm_calib.h"
#include "FocalPlaneMask/exec_optimize_fpm_zones.h"
#include "FocalPlaneMask/exec_optimize_fpmtransmission.h"

#include "LyotStop/exec_optimize_lyot_stop_position.h"
#include "LyotStop/exec_optimize_lyot_stops_shapes_positions.h"

#include "PIAAshape/exec_optimize_PIAA_shapes.h"
#include "PIAAshape/makePIAAshapes.h"

static double PIAACMCsimul_regularization_PIAAshapes_value()
{
    DEBUG_TRACE_FSTART();

    double value = 0.0;
    imageID IDref;
    imageID ID;
    imageID ID_CPAfreq;

    // add regularization component of the evaluation metrix
    // first, compute and add PIAA shape regularization value (value) if applicable
    // this is the same as the previous computation of val0

    // index of the PIAA element 0 (first mirror) shapes via cosine modes
    ID = piaacmcopticaldesign.piaa0CmodesID;

    // index of PIAA shapes reference image
    IDref = image_ID("piaa0Cmref");
    if (IDref == -1)
    {
        // error message if we get here?  ***************************************
        // if the reference image doesn't exist, create it
        create_2Dimage_ID("piaa0Cmref", data.image[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0], 1, &IDref);

        // initialize to zero shape
        for (uint32_t jj = 0; jj < data.image[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0]; jj++)
        {
            data.image[IDref].array.F[jj] = 0.0;
        }
    }

    // This section of code does not actually fill in the regularization terms in the output vector
    // filling in is done later.  Here we are only computing the initial reference scalar objective value

    // For each cosine mode set the optimization parameter = cosine mode modified by regularization
    for (uint32_t jj = 0; jj < data.image[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0]; jj++)
    {
        // compute square of C*(deviation from reference)*(mode index)^(alpha)
        // so higher-index zones (higher spatial frequency) are more
        // heavily penalized

        double tmp;
        tmp = piaacmcparams.linopt_piaa0C_regcoeff * (data.image[ID].array.F[jj] - data.image[IDref].array.F[jj]) *
              pow(1.0 * jj, piaacmcparams.linopt_piaa0C_regcoeff_alpha);
        value += tmp * tmp;
    }

    // do the same for PIAA element 1
    ID = piaacmcopticaldesign.piaa1CmodesID;
    IDref = image_ID("piaa1Cmref");
    if (IDref == -1)
    {
        create_2Dimage_ID("piaa1Cmref", data.image[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0], 1, &IDref);

        for (uint32_t jj = 0; jj < data.image[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0]; jj++)
        {
            data.image[IDref].array.F[jj] = 0.0;
        }
    }
    for (uint32_t jj = 0; jj < data.image[piaacmcopticaldesign.piaa1CmodesID].md[0].size[0]; jj++)
    {
        double tmp;
        tmp = piaacmcparams.linopt_piaa1C_regcoeff * (data.image[ID].array.F[jj] - data.image[IDref].array.F[jj]) *
              pow(1.0 * jj, piaacmcparams.linopt_piaa1C_regcoeff_alpha);
        value += tmp * tmp;
    }

    // get spatial frequency of each mode in cycles/aperture
    ID_CPAfreq = image_ID("cpamodesfreq");

    // do the same for PIAA element 0 and 1 for the Fourier modes
    // this time use the actual spatial frequency rather than mode index as proxy for frequency
    ID = piaacmcopticaldesign.piaa0FmodesID;
    IDref = image_ID("piaa0Fmref");
    if (IDref == -1)
    {
        create_2Dimage_ID("piaa0Fmref", data.image[piaacmcopticaldesign.piaa0FmodesID].md[0].size[0], 1, &IDref);

        for (uint32_t jj = 0; jj < data.image[piaacmcopticaldesign.piaa0FmodesID].md[0].size[0]; jj++)
        {
            data.image[IDref].array.F[jj] = 0.0;
        }
    }
    for (uint32_t jj = 0; jj < data.image[piaacmcopticaldesign.piaa0FmodesID].md[0].size[0]; jj++)
    {
        double tmp;
        tmp = piaacmcparams.linopt_piaa0F_regcoeff * (data.image[ID].array.F[jj] - data.image[IDref].array.F[jj]) *
              pow(1.0 * data.image[ID_CPAfreq].array.F[jj], piaacmcparams.linopt_piaa0F_regcoeff_alpha);
        value += tmp * tmp;
    }

    ID = piaacmcopticaldesign.piaa1FmodesID;
    IDref = image_ID("piaa1Fmref");
    if (IDref == -1)
    {
        create_2Dimage_ID("piaa1Fmref", data.image[piaacmcopticaldesign.piaa1FmodesID].md[0].size[0], 1, &IDref);

        for (uint32_t jj = 0; jj < data.image[piaacmcopticaldesign.piaa1FmodesID].md[0].size[0]; jj++)
        {
            data.image[IDref].array.F[jj] = 0.0;
        }
    }
    for (uint32_t jj = 0; jj < data.image[piaacmcopticaldesign.piaa1FmodesID].md[0].size[0]; jj++)
    {
        double tmp;
        tmp = piaacmcparams.linopt_piaa1F_regcoeff * (data.image[ID].array.F[jj] - data.image[IDref].array.F[jj]) *
              pow(1.0 * data.image[ID_CPAfreq].array.F[jj], piaacmcparams.linopt_piaa1F_regcoeff_alpha);
        value += tmp * tmp;
    }

    DEBUG_TRACE_FEXIT();
    return value;
}

/**
 * @brief Compute regularization term for focal plane mask sag values
 *
 */

static double PIAACMCsimul_regularization_fpmsag_value()
{
    DEBUG_TRACE_FSTART();

    double regvalue; // output regularization value
    imageID IDzonez; // image: zone sag values

    // regvalue is the regularization value for the focal plane mask sag values
    // same as above, but not dependent on position
    regvalue = 0.0;
    IDzonez = piaacmcopticaldesign.zonezID;

    for (long zoneindex = 0; zoneindex < data.image[IDzonez].md[0].size[0]; zoneindex++)
    {
        // compute the square of (sag/coeff)^alpha
        double tmp;
        tmp = pow(data.image[IDzonez].array.D[zoneindex] / piaacmcopticaldesign.fpmsagreg_coeff,
                  piaacmcopticaldesign.fpmsagreg_alpha);
        regvalue += tmp * tmp;
    }

    DEBUG_TRACE_FEXIT();
    return regvalue;
}

static long PIAACMCsimul_regularization_PIAAshapes_add1Dvector(imageID ID1D, long index0)
{
    DEBUG_TRACE_FSTART();

    long vindex; // index in output vector

    imageID IDref;
    imageID ID_CPAfreq;
    imageID ID;

    // fill in the shape regularization value's response to piaacmcparams.linopt_paramdelta using the
    // same formulas as before with the delta PSF as input
    ID = piaacmcopticaldesign.piaa0CmodesID;
    vindex = index0;

    IDref = image_ID("piaa0Cmref");
    if (IDref == -1)
    {
        create_2Dimage_ID("piaa0Cmref", data.image[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0], 1, &IDref);

        for (uint32_t ii = 0; ii < data.image[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0]; ii++)
        {
            data.image[IDref].array.F[ii] = 0.0;
        }
    }

    for (uint32_t ii = 0; ii < data.image[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0]; ii++)
    {
        data.image[ID1D].array.F[vindex] = piaacmcparams.linopt_piaa0C_regcoeff *
                                           (data.image[ID].array.F[ii] - data.image[IDref].array.F[ii]) *
                                           pow(1.0 * ii, piaacmcparams.linopt_piaa0C_regcoeff_alpha);
        vindex++;
    }

    ID = piaacmcopticaldesign.piaa1CmodesID;
    IDref = image_ID("piaa1Cmref");
    if (IDref == -1)
    {
        create_2Dimage_ID("piaa1Cmref", data.image[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0], 1, &IDref);

        for (uint32_t ii = 0; ii < data.image[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0]; ii++)
        {
            data.image[IDref].array.F[ii] = 0.0;
        }
    }
    for (uint32_t ii = 0; ii < data.image[piaacmcopticaldesign.piaa1CmodesID].md[0].size[0]; ii++)
    {
        data.image[ID1D].array.F[vindex] = piaacmcparams.linopt_piaa1C_regcoeff *
                                           (data.image[ID].array.F[ii] - data.image[IDref].array.F[ii]) *
                                           pow(1.0 * ii, piaacmcparams.linopt_piaa1C_regcoeff_alpha);
        vindex++;
    }

    ID_CPAfreq = image_ID("cpamodesfreq");

    ID = piaacmcopticaldesign.piaa0FmodesID;
    IDref = image_ID("piaa0Fmref");
    if (IDref == -1)
    {
        create_2Dimage_ID("piaa0Fmref", data.image[piaacmcopticaldesign.piaa0FmodesID].md[0].size[0], 1, &IDref);

        for (uint32_t ii = 0; ii < data.image[piaacmcopticaldesign.piaa0FmodesID].md[0].size[0]; ii++)
        {
            data.image[IDref].array.F[ii] = 0.0;
        }
    }
    for (uint32_t ii = 0; ii < data.image[piaacmcopticaldesign.piaa0FmodesID].md[0].size[0]; ii++)
    {
        data.image[ID1D].array.F[vindex] =
            piaacmcparams.linopt_piaa0F_regcoeff * (data.image[ID].array.F[ii] - data.image[IDref].array.F[ii]) *
            pow(1.0 * data.image[ID_CPAfreq].array.F[ii], piaacmcparams.linopt_piaa0F_regcoeff_alpha);
        vindex++;
    }

    ID = piaacmcopticaldesign.piaa1FmodesID;
    IDref = image_ID("piaa1Fmref");
    if (IDref == -1)
    {
        create_2Dimage_ID("piaa1Fmref", data.image[piaacmcopticaldesign.piaa1FmodesID].md[0].size[0], 1, &IDref);

        for (uint32_t ii = 0; ii < data.image[piaacmcopticaldesign.piaa1FmodesID].md[0].size[0]; ii++)
        {
            data.image[IDref].array.F[ii] = 0.0;
        }
    }
    for (uint32_t ii = 0; ii < data.image[piaacmcopticaldesign.piaa1FmodesID].md[0].size[0]; ii++)
    {
        data.image[ID1D].array.F[vindex] =
            piaacmcparams.linopt_piaa1F_regcoeff * (data.image[ID].array.F[ii] - data.image[IDref].array.F[ii]) *
            pow(1.0 * data.image[ID_CPAfreq].array.F[ii], piaacmcparams.linopt_piaa1F_regcoeff_alpha);
        vindex++;
    }

    DEBUG_TRACE_FEXIT();
    return vindex;
}

/**
 *
 * @brief Main simulation routine
 *
 *
 * @param[in] confindex  PIAACMC configuration index pointing to the input/output directory number
 * @param[in] mode       Type of operation to be performed
 *
 */
errno_t PIAACMCsimul_exec(const char *confindex, long mode)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %s %ld", confindex, mode);

    double valref;

    imageID IDv, ID;
    imageID IDmodes, IDmodes2D;

    imageID IDm, ID1D, ID1Dref;

    double contrastval;

    char stopfile[STRINGMAXLEN_FULLFILENAME];

    piaacmcparams.linopt_REGPIAASHAPES = 0;

    piaacmcparams.linopt_piaa0C_regcoeff = 0.0e-7;
    piaacmcparams.linopt_piaa1C_regcoeff = 0.0e-7;
    piaacmcparams.linopt_piaa0C_regcoeff_alpha = 1.0;
    piaacmcparams.linopt_piaa1C_regcoeff_alpha = 1.0;

    piaacmcparams.linopt_piaa0F_regcoeff = 0.0e-7;
    piaacmcparams.linopt_piaa1F_regcoeff = 0.0e-7;
    piaacmcparams.linopt_piaa0F_regcoeff_alpha = 1.0;
    piaacmcparams.linopt_piaa1F_regcoeff_alpha = 1.0;

    piaacmcparams.linopt_REGFPMSAG = 0;

    imageID IDstatus = -1;
    { // Create status shared variable
        // this allows realtime monitoring of the code by other processes
        // sets status at different points in the code
        // has 4 pixels :
        // - pixel 0    major step
        // - pixel 1    minor step
        // - pixel 2    source code file
        // - pixel 3    source code line number
        IDstatus = image_ID("stat_PIAACMCsimulexec");
        if (IDstatus == -1)
        {
            printf("Looking for stat_PIAACMCsimulexec\n");
            fflush(stdout);
            IDstatus = read_sharedmem_image("stat_PIAACMCsimulexec");
            printf("ID for stat_PIAACMCsimulexec: %ld\n", IDstatus);
            fflush(stdout);
        }
        if (IDstatus == -1)
        {
            uint32_t *sizearray;
            sizearray = (uint32_t *)malloc(sizeof(uint32_t) * 2);
            if (sizearray == NULL)
            {
                FUNC_RETURN_FAILURE("malloc returns NULL pointer");
            }
            sizearray[0] = 4;
            sizearray[1] = 1;
            create_image_ID("stat_PIAACMCsimulexec", 2, sizearray, _DATATYPE_UINT16, 1, 0, 0, &IDstatus);
            free(sizearray);
        }
    }

    for (int elem = 0; elem < 100; elem++)
    {
        piaacmcopticalsystem.keepMem[elem] = 0; // flag that says save this element for reuse
    }

    // set the result directories
    snprintf(piaacmcparams.piaacmcconfdir, STRINGMAXLEN_DIRNAME, "%s", confindex);
    snprintf(data.SAVEDIR, STRINGMAXLEN_DIRNAME, "%s", piaacmcparams.piaacmcconfdir);

    piaacmcparams.linopt_NBiter = 1000;
    piaacmcparams.linopt_number_param = 0;

    piaacmcopticalsystem.SAVE = piaacmcparams.PIAACMC_save;

    // get variables from command line, possibly sets globals
    /*    if( (IDv = variable_ID("PIAACMC_centobs0")) != -1)
            centobs0 = data.variable[IDv].value.f;
        if( (IDv = variable_ID("PIAACMC_centobs1")) != -1)
            centobs1 = data.variable[IDv].value.f;
            */
    double fpmradld = 0.95;
    if ((IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }

    // start a log of code mode entry/exit times
    EXECUTE_SYSTEM_COMMAND("echo \"%03ld     $(date)\" >> ./log/PIAACMC_mode_log.txt", mode);

    // set the name of the stopfile
    WRITE_FULLFILENAME(stopfile, "%s/stopmode%ld.txt", piaacmcparams.piaacmcconfdir, mode);

    switch (mode)
    {
    case 0:
        FUNC_CHECK_RETURN(exec_compute_image());
        printf("EXEC CASE 0 COMPLETED\n");
        fflush(stdout);
        break;

    case 1:
        FUNC_CHECK_RETURN(exec_optimize_lyot_stop_position());
        break;

    case 2:
        FUNC_CHECK_RETURN(exec_optimize_fpmtransmission());
        break;

    case 3:
        FUNC_CHECK_RETURN(exec_computePSF_no_fpm(NULL));
        break;

    case 4:
        FUNC_CHECK_RETURN(exec_optimize_PIAA_shapes());
        break;

    case 5:
        FUNC_CHECK_RETURN(exec_optimize_lyot_stops_shapes_positions());
        break;

    case 11:
        FUNC_CHECK_RETURN(exec_multizone_fpm_calib());
        break;

    case 13:
        FUNC_CHECK_RETURN(exec_optimize_fpm_zones());
        break;

    case 40:
        FUNC_CHECK_RETURN(exec_optimize_PIAA_shapes_fpmtransm());
        break;

    case 100: // evaluate current design: polychromatic contrast, pointing sensitivity
        FUNC_CHECK_RETURN(PIAACMCsimul_eval_poly_design());
        break;

    case 101:
        FUNC_CHECK_RETURN(PIAACMCsimul_measure_transm_curve());
        break;

    case 302: // restore configuration settings
        printf("=================================== mode 302 ===================================\n");

        EXECUTE_SYSTEM_COMMAND("cp %s/saveconf/conf_*.txt %s/", piaacmcparams.piaacmcconfdir,
                               piaacmcparams.piaacmcconfdir);
        break;

    default:
        PRINT_ERROR("mode not recognized");
        break;
    }

    // linear optimization set up in modes 13 and 40
    //
    // the output parameters start as the evaluation zone values in "imvect"
    // If we are regularizing, we supplement the output parameters by adding
    // penalty terms as additional output parameters
    if (piaacmcparams.LINOPT == 1) // linear optimization
    {
        // for state tracking and statistics
        data.image[IDstatus].array.UI16[0] = 5;

        // Compute Reference on-axis performance contrast (valref)
        FUNC_CHECK_RETURN(makePIAAshapes());

        piaacmcopticalsystem.FOCMASKarray[0].mode = 1; // use 1-fpm
        {
            FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(0.0, 0.0, 0, piaacmcopticalsystem.NBelem, 0,
                                                      piaacmcparams.computePSF_ResolvedTarget,
                                                      piaacmcparams.computePSF_ResolvedTarget_mode, 0, &valref));
        }

        {
            // save the current configuration to the _linopt directory
            char dirname[STRINGMAXLEN_DIRNAME];

            WRITE_DIRNAME(dirname, "%s_linopt", piaacmcparams.piaacmcconfdir);

            FUNC_CHECK_RETURN(PIAACMCsimul_savepiaacmcconf(dirname));

            // import configuration from _linopt directory
            EXECUTE_SYSTEM_COMMAND("rsync -au --progress %s/* ./%s/", dirname, piaacmcparams.piaacmcconfdir);

            // "cp -n" will only copy the file if destination does not exist
            // this will ensure that the .ref.fits files are the PIAA shapes before any linear optimization
            // these will be the "reference" shapes used to regularize PIAA shapes: the deviation from this reference will be kept small by the regularization

            EXECUTE_SYSTEM_COMMAND("cp -n %s/piaa0Cmodes.fits %s/piaa0Cmodes.ref.fits", dirname, dirname);
            EXECUTE_SYSTEM_COMMAND("cp -n %s/piaa0Fmodes.fits %s/piaa0Fmodes.ref.fits", dirname, dirname);

            EXECUTE_SYSTEM_COMMAND("cp -n %s/piaa1Cmodes.fits %s/piaa1Cmodes.ref.fits", dirname, dirname);
            EXECUTE_SYSTEM_COMMAND("cp -n %s/piaa1Fmodes.fits %s/piaa1Fmodes.ref.fits", dirname, dirname);

            EXECUTE_SYSTEM_COMMAND("cp -n %s/piaacmcparams.conf %s/piaacmcparams.ref.conf", dirname, dirname);

            // load PIAA reference shapes (used for regularization)

            char fname[STRINGMAXLEN_FULLFILENAME];

            WRITE_FULLFILENAME(fname, "%s/piaa0Cmodes.ref.fits", dirname);
            FUNC_CHECK_RETURN(load_fits(fname, "piaa0Cmref", 1, NULL));

            WRITE_FULLFILENAME(fname, "%s/piaa1Cmodes.ref.fits", dirname);
            FUNC_CHECK_RETURN(load_fits(fname, "piaa1Cmref", 1, NULL));

            WRITE_FULLFILENAME(fname, "%s/piaa0Fmodes.ref.fits", dirname);
            FUNC_CHECK_RETURN(load_fits(fname, "piaa0Fmref", 1, NULL));

            WRITE_FULLFILENAME(fname, "%s/piaa1Fmodes.ref.fits", dirname);
            FUNC_CHECK_RETURN(load_fits(fname, "piaa1Fmref", 1, NULL));
        }

        // we have now saved the starting point of the optimization for future comparison
        // in the <piaacmcconfdir>_linopt directory

        // here we compute regularization value of piaashapes and store it in the val0 variable
        // if regularization is of PIAA shapes is ON, then val0 will be computed and added to the overal performance metric valref

        // regularize the piaashapes via a penalty added to the reference contrast valref
        // The optimization minimizes the summed contrast + val0 + val1.
        // Regularization is via adding a constant val0 + val1 to the contrast we're minimizing
        // note that here we're setting output parameters.
        if (piaacmcparams.linopt_REGPIAASHAPES == 1)
        {
            double val0 = PIAACMCsimul_regularization_PIAAshapes_value();
            printf("VALREF = %g + %g -> %g\n", valref, val0, valref + val0);
            valref += val0;
        }

        // val1 is the regularization value for the focal plane mask sag values
        // same as above, but not dependent on position
        //        double val1;
        if (piaacmcparams.linopt_REGFPMSAG == 1)
        {
            double val1 = PIAACMCsimul_regularization_fpmsag_value();
            valref += val1;
        }
        /*      else
              {
                  val1 = 0.0;
              }*/

        // At this point all we've done is compute the overall performance metric including
        // regularization in valref.

        // for state tracking and statistics
        data.image[IDstatus].array.UI16[0] = 6;
        printf("================================ Reference = %g\n", valref);

        // copy imvect to vecDHref "vector dark hole reference"
        // vecDHref is the dark hole complex amplitude state at the beginning of the linear optimization
        // (this dark hole is nominally a full annulus from 1.5 to ~8 lambda/D, created at the top
        // of PIAACMCsimul_computePSF with size controlled by scoringIWA and scoringOWA
        // the corresponding performance metric is valref
        chname_image_ID("imvect",
                        "vecDHref"); // note: imvect was computed by PIAACMCsimul_computePSF called ~150 lines above
        ID = image_ID("vecDHref");   // ID changed identity
        //xsize = data.image[ID].md[0].size[0];
        //ysize = data.image[ID].md[0].size[1];

        // for state tracking and statistics
        data.image[IDstatus].array.UI16[0] = 7;

        // now we will just determine the size of the size of the
        // optimization vectors that we will actually fill in later
        // save vecDHref initial state as a reference
        {
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(fname, "%s/vecDMref.fits", piaacmcparams.piaacmcconfdir);
            save_fits("vecDHref", fname);
        }
        // get and copy size of vecDHref, 'cause we're manipulating size1Dvec
        long size1Dvec = data.image[ID].md[0].nelement;
        long size1Dvec0 = size1Dvec;

        // PIAA shapes regularization
        // if regularization is turned on, the size of the evaluation vector is increased to include PIAA shape coefficients, in addition to complex amplitude in focal plane
        // the optimization code will then simultaneously minimize the light in the focal plane AND the PIAA shape deviations from the nominal shape
        if (piaacmcparams.linopt_REGPIAASHAPES == 1)
        {
            // there are 4 groups of PIAA shape parameters: 2 sets of cosine modes and 2 sets of Fourier modes
            size1Dvec += data.image[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0];
            size1Dvec += data.image[piaacmcopticaldesign.piaa1CmodesID].md[0].size[0];
            size1Dvec += data.image[piaacmcopticaldesign.piaa0FmodesID].md[0].size[0];
            size1Dvec += data.image[piaacmcopticaldesign.piaa1FmodesID].md[0].size[0];
        }

        // The same approach is used for regularization of the focal plane mask sag values
        // the sag values are appended to the evaluation vector
        if (piaacmcparams.linopt_REGFPMSAG == 1)
        {
            size1Dvec += data.image[piaacmcopticaldesign.zonezID].md[0].size[0];
        }

        // re-package vector into 1D array and add regularization terms
        // the resulting vector is stored in vecDHref1D
        // we also create a mask image which is intended to mask out some of the pixels from the evaluation
        // the mask is currently not used, so we will write 1.0 in all of its pixels
        // DHmask and vecDHref1D contain all the optimization parameters as set above
        // "ID of evaluation mode mask"
        FUNC_CHECK_RETURN(create_2Dimage_ID("DHmask", size1Dvec, 1, &IDm));

        // "ID of 1D dark zone reference"
        FUNC_CHECK_RETURN(create_2Dimage_ID("vecDHref1D", size1Dvec, 1, &ID1Dref));

        // we first write 1.0 into the focal plane complex amplitudes in the vector
        ID = image_ID("vecDHref");
        {
            uint64_t ii;
            for (ii = 0; ii < data.image[ID].md[0].nelement; ii++)
            {
                // imbed vecDHref into the evaluation zone part of of the full parameter vector vecDHref1D
                data.image[ID1Dref].array.F[ii] = data.image[ID].array.F[ii];
                // sets the evaluation zone part of of the full parameter vector vecDHref1D to 1
                // 1 means the evaluation zone is on
                data.image[IDm].array.F[ii] = 1.0;
            }
            // !!!!!! WARNING !!!!!!!
            // the state of ii at this point drives the code below and will evolve until the comment
            // that says we're done with ii

            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 8;

            // Now actually fill in the regularized output vector.
            // If we are not regularizing, the output evaluation zone values are filled in by the
            // PSF calls in the optimization loop
            // and then append regularization vectors in the main evaluation vector
            // Note: index ii is incremented as we add "sub-vectors" into the main evaluation vector
            if (piaacmcparams.linopt_REGPIAASHAPES == 1)
            {
                // initialize by filling in the regularization terms of the output,
                // !! starting at the current value of ii !!
                // This means we're actually writing the output vector regularization terms.
                // Otherwise this is the same as the if(REGPIAASHAPES == 1) code block above
                ii = PIAACMCsimul_regularization_PIAAshapes_add1Dvector(ID1Dref, ii);
            }

            // same for the sags, starting at the current value of ii
            if (piaacmcparams.linopt_REGFPMSAG == 1)
            {
                ID = piaacmcopticaldesign.zonezID;
                for (uint32_t jj = 0; jj < data.image[ID].md[0].size[0]; jj++)
                {
                    data.image[ID1Dref].array.F[ii] =
                        pow(data.image[ID].array.D[jj] / piaacmcopticaldesign.fpmsagreg_coeff,
                            piaacmcopticaldesign.fpmsagreg_alpha);
                    data.image[IDm].array.F[ii] = 1.0;
                    ii++;
                }
            }
            // !!!!!!
            // we're done with ii
        }

        // vecDHref has beem embedded into vecDHref1D
        FUNC_CHECK_RETURN(delete_image_ID("vecDHref", DELETE_IMAGE_ERRMODE_WARNING));

        // at this point, we have completed the initialization, and the optimization loop starts

        // file that will track optimization loop progress
        {
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(fname, "%s/linoptval.txt", piaacmcparams.piaacmcconfdir);

            FILE *fp;
            fp = fopen(fname, "w");
            fclose(fp);
        }
        //	list_image_ID();
        //
        // LINEAR OPTIMIZATION AROUND CURRENT POINT
        //

        double bestval = 0.0;

        int iterOK = 1;
        long iter = 0;
        double oldval = 1.0;
        data.image[IDstatus].array.UI16[0] = 9;

        // while # of iterations < piaacmcparams.linopt_NBiter
        //  and the ojective changes by more than 2% after the second iteration
        //  and something about NBlinoptgain ???????
        while (iterOK == 1)
        {
            double bestgain = 0.0;

            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 10;
            printf("Iteration %ld/%ld\n", iter, piaacmcparams.linopt_NBiter);
            fflush(stdout);

            // array for collecting dark hole mode derivatives
            // stores derivative of output vector against input parameters
            FUNC_CHECK_RETURN(create_3Dimage_ID("DHmodes", size1Dvec, 1, piaacmcparams.linopt_number_param, &IDmodes));
            // 2D array for diagnostic display
            FUNC_CHECK_RETURN(create_2Dimage_ID("DHmodes2D", size1Dvec, piaacmcparams.linopt_number_param, &IDmodes2D));

            // get ready to update optimization tracking file
            {
                FILE *fp;
                char fname[STRINGMAXLEN_FULLFILENAME];

                WRITE_FULLFILENAME(fname, "%s/linoptval.txt", piaacmcparams.piaacmcconfdir);
                fp = fopen(fname, "a");
                fprintf(fp, "### PIAACMC_FPM_FASTDERIVATIVES = %d\n", piaacmcparams.PIAACMC_FPM_FASTDERIVATIVES);
                fclose(fp);
            }
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 11;

            printf("Compute local derivatives of output vector against input focal plane mask zones (sags)\n");
            fflush(stdout);

            if (piaacmcparams.PIAACMC_FPM_FASTDERIVATIVES == 1) // TO BE USED ONLY FOR FOCAL PLANE MASK OPTIMIZATION
            {
                // this only happens in mode 13
                // the fast derivative mode only works for focal plane mask optimization, for which derivatives against sag values can be comptuted by simple rotation of pre-computed vectors from mode 11
                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 12;
                piaacmcopticalsystem.FOCMASKarray[0].mode = 1; // use 1-fpm
                //				ID = create_2Dimage_ID("DHmodes2Dtest", size1Dvec, piaacmcparams.linopt_number_param);

                printf("Computing %ld derivatives ", (long)data.image[piaacmcopticaldesign.zonezID].md[0].size[0]);
                fflush(stdout);
                for (uint32_t mz = 0; mz < data.image[piaacmcopticaldesign.zonezID].md[0].size[0];
                     mz++) // loop over mask zones
                {
                    printf(" %ld", (long)mz);
                    // actually compute the derivative
                    // fpmresp_array is results from mode 11
                    // from mode 13 above:
                    //      fpmresp_array = data.image[IDfpmresp].array.D;
                    //      zonez_array = data.image[piaacmcopticaldesign.zonezID].array.D;
                    // dphadz_array was computed in mode 13 shortly afterwards
                    // outtmp_array is output
                    FUNC_CHECK_RETURN(PIAACMCsimul_achromFPMsol_eval_zonezderivative(
                        mz, piaacmcparams.fpmresp_array, piaacmcparams.zonez_array, piaacmcparams.dphadz_array,
                        piaacmcparams.outtmp_array, piaacmcparams.vsize,
                        data.image[piaacmcopticaldesign.zonezID].md[0].size[0], piaacmcopticaldesign.nblambda));
                    for (long ii = 0; ii < size1Dvec0; ii++)
                    {
                        data.image[IDmodes].array.F[mz * size1Dvec + ii] =
                            piaacmcparams.outtmp_array[ii] * piaacmcparams.linopt_paramdelta[mz];
                    }
                }
                // derivatives of regularization values against sag values can also be computed analytically without requiring diffraction propagation
                if (piaacmcparams.linopt_REGFPMSAG == 1)
                {
                    ID = piaacmcopticaldesign.zonezID;
                    // following should be derivative of (sag/coeff)^alpha
                    // w.r.t. sag
                    for (uint32_t mz = 0; mz < data.image[ID].md[0].size[0]; mz++)
                    {
                        data.image[IDmodes].array.F[mz * size1Dvec + (size1Dvec0 + mz)] =
                            (piaacmcopticaldesign.fpmsagreg_alpha / piaacmcopticaldesign.fpmsagreg_coeff) *
                            pow(data.image[ID].array.D[mz] / piaacmcopticaldesign.fpmsagreg_coeff,
                                piaacmcopticaldesign.fpmsagreg_alpha - 1.0) *
                            piaacmcparams.linopt_paramdelta[mz];
                    }
                }

                // TEST diagnostic
                memcpy(data.image[IDmodes2D].array.F, data.image[IDmodes].array.F,
                       sizeof(float) * size1Dvec * piaacmcparams.linopt_number_param);

                FUNC_CHECK_RETURN(save_fl_fits("DHmodes2D", "test_DHmodes2D.fits"));

                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 13;

                printf("Done computing derivatives (FAST MODE)\n");
                fflush(stdout);
            }
            else // ONLY FOR PIAA SHAPES OPTIMIZATION
            {
                // derivatives against PIAA shapes must be computed numerically

                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 14;
                for (int i = 0; i < piaacmcparams.linopt_number_param; i++)
                {
                    piaacmcopticalsystem.FOCMASKarray[0].mode = 1; // use 1-fpm
                    // get delta on-axis PSF response to the change in piaacmcparams.linopt_paramdelta
                    // to later give derivative w.r.t. piaacmcparams.linopt_paramdelta
                    if (piaacmcparams.linopt_paramtype[i] == _DATATYPE_FLOAT)
                    {
                        *(piaacmcparams.linopt_paramvalf[i]) += (float)piaacmcparams.linopt_paramdelta[i];
                        FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(
                            0.0, 0.0, 0, piaacmcopticalsystem.NBelem, 0, piaacmcparams.computePSF_ResolvedTarget,
                            piaacmcparams.computePSF_ResolvedTarget_mode, 0, &contrastval));
                    }
                    else
                    {
                        *(piaacmcparams.linopt_paramval[i]) += piaacmcparams.linopt_paramdelta[i];
                        FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(
                            0.0, 0.0, 0, piaacmcopticalsystem.NBelem, 0, piaacmcparams.computePSF_ResolvedTarget,
                            piaacmcparams.computePSF_ResolvedTarget_mode, 0, &contrastval));
                    }

                    //      sprintf(fname,"%s/imvect_%02ld.fits", piaacmcparams.piaacmcconfdir, i);
                    //       save_fits("imvect", fname);
                    ID = image_ID("imvect");

                    {
                        char fname[STRINGMAXLEN_FULLFILENAME];
                        WRITE_FULLFILENAME(fname, "%s/linoptval.txt", piaacmcparams.piaacmcconfdir);
                        FILE *fp = fopen(fname, "a");
                        fprintf(fp, "# %5ld/%5ld %5d/%5ld %6ld %6ld    %20.15g %20.15g %20.15g      %20.15g\n", iter,
                                piaacmcparams.linopt_NBiter, i, piaacmcparams.linopt_number_param,
                                data.image[ID].md[0].nelement, size1Dvec, piaacmcparams.linopt_paramdelta[i],
                                contrastval, valref, bestval);
                        fclose(fp);
                    }
                    // re-package vector into 1D array and add regularization terms
                    // evaluation vector is "imvect1D", ID = ID1D
                    // similar to vecDHref before
                    FUNC_CHECK_RETURN(create_2Dimage_ID("imvect1D", size1Dvec, 1, &ID1D));
                    // fill in the evaluation point portion
                    {
                        uint64_t ii;
                        for (ii = 0; ii < data.image[ID].md[0].nelement; ii++)
                        {
                            data.image[ID1D].array.F[ii] = data.image[ID].array.F[ii];
                        }

                        if (piaacmcparams.linopt_REGPIAASHAPES == 1)
                        {
                            // fill in the shape regularization value's response to piaacmcparams.linopt_paramdelta using the
                            // same formulas as before with the delta PSF as input
                            ii = PIAACMCsimul_regularization_PIAAshapes_add1Dvector(ID1D, ii);
                        }
                    }

                    FUNC_CHECK_RETURN(
                        delete_image_ID("imvect", DELETE_IMAGE_ERRMODE_WARNING)); // has been imbedded into imvect1D

                    // restore original state (return to original staring point)
                    if (piaacmcparams.linopt_paramtype[i] == _DATATYPE_FLOAT)
                    {
                        *(piaacmcparams.linopt_paramvalf[i]) -= (float)piaacmcparams.linopt_paramdelta[i];
                    }
                    else
                    {
                        *(piaacmcparams.linopt_paramval[i]) -= piaacmcparams.linopt_paramdelta[i];
                    }

                    // compute actual derivative as first difference from reference
                    // this is the starting derivative
                    for (uint64_t ii = 0; ii < data.image[ID1D].md[0].nelement; ii++)
                    {
                        data.image[IDmodes].array.F[i * data.image[ID1D].md[0].nelement + ii] =
                            (data.image[ID1D].array.F[ii] - data.image[ID1Dref].array.F[ii]);
                    }

                    //    printf("%3ld %g %g\n", i, val, valref);

                    // create diagnostic image
                    FUNC_CHECK_RETURN(
                        create_2Dimage_ID("DHmodes2D", size1Dvec, piaacmcparams.linopt_number_param, &ID));

                    for (uint64_t ii = 0; ii < data.image[IDmodes].md[0].nelement; ii++)
                    {
                        data.image[ID].array.F[ii] = data.image[IDmodes].array.F[ii];
                    }

                    {
                        char fname[STRINGMAXLEN_FULLFILENAME];
                        WRITE_FULLFILENAME(fname, "%s/DMmodes.fits", piaacmcparams.piaacmcconfdir);
                        FUNC_CHECK_RETURN(save_fits("DHmodes2D", fname));
                    }

                    FUNC_CHECK_RETURN(delete_image_ID("DHmodes2D", DELETE_IMAGE_ERRMODE_WARNING));
                }
                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 15;
            }

            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 16;

            { // print the results to file for human tracking
                char fname[STRINGMAXLEN_FULLFILENAME];
                WRITE_FULLFILENAME(fname, "%s/linoptval.txt", piaacmcparams.piaacmcconfdir);
                FILE *fp = fopen(fname, "a");
                fprintf(fp, "### scanning gain \n");
                fprintf(fp, "### <alphareg>  <gain>  <contrast>\n");
                fclose(fp);
            }
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 17;

            // first three arguments are names of the input arrays
            // vecDHref1D is the input data
            // DHmodes is the basis of modes to expand vecDHref1D into
            // DHmask weights elements of vecDHref1D (nominally all weights = 1)
            // 4th arg is the pseudoinverse eigenvalue (via eigenvalue decomposition)
            // this decomposes vecDHref1D into the DHmodes
            // 5th arg is the output: optcoeff0*DHmodes = vecDHref1D
            // computed via pseudoinverse of DHmodes
            // This decomposition facilitates the cancellation of vecDHref1D by
            // searching in the DHmodes basis
            //
            // use three cutoff values to give three options for future evaluation
            // smallest cutoff values produce the largest changes (are least well conditioned)
            printf("ref = 0.1   -- ");
            fflush(stdout);
            linopt_imtools_image_fitModes("vecDHref1D", "DHmodes", "DHmask", 0.1, "optcoeff0", 0, NULL);
            printf("- DONE\n");
            fflush(stdout);

            printf("ref = 0.01  -- ");
            fflush(stdout);
            linopt_imtools_image_fitModes("vecDHref1D", "DHmodes", "DHmask", 0.01, "optcoeff1", 0, NULL);
            printf("- DONE\n");
            fflush(stdout);

            printf("ref = 0.001 -- ");
            fflush(stdout);
            linopt_imtools_image_fitModes("vecDHref1D", "DHmodes", "DHmask", 0.001, "optcoeff2", 0, NULL);
            printf("- DONE\n");
            fflush(stdout);

            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 18;

            // initialize zero "optimal" vector optvec giving direction to move in the search for the min
            arith_image_cstmult("optcoeff0", 0.0, "optvec"); // create optimal vector
            imageID IDoptvec = image_ID("optvec");
            // initialize the objective value
            int initbestval = 0;
            bestval = valref;

            double alphareg = 1.0; // has no effect (see next loop)

            // say something here ???
            long NBlinoptgain = 0;

            float scangainfact = 1.2;
            //alphascaninit = 0;
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 19;
            // alphareg controls linear combinations of the directions optcoeff0,1,2
            // alphareg = 0 => moving along optcoeff0
            // alphareg = 0.5 => moving along optcoeff1
            // alphareg = 1 => moving along optcoeff2
            // with interpolation
            // look in 5 steps if alphareg += 0.2
            for (alphareg = 0.0; alphareg < 1.01; alphareg += 0.2)
            {
                // produce piecewise linear interoplation coefficients
                double acoeff0 = 1.0 - 2.0 * alphareg; // ranges from -1 to 1
                if (acoeff0 < 0.0)
                {
                    acoeff0 = 0.0; // clip to 0 to 1, = 0 if alphareg > 0.5
                }

                double acoeff1 = 1.0 - fabs(2.0 * (alphareg - 0.5)); // two lines: from 0,1 to .5,0 to 1,1

                double acoeff2 = 2.0 * alphareg - 1.0; // ranges from -1 to 1
                if (acoeff2 < 0.0)
                {
                    acoeff2 = 0.0; // clip to 0 to 1, = 0 if alphareg < 0.5
                }

                // sum of acoeff0,1,2 = 1 at all values of alphareg.

                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 20;

                // optcoeff0m = acoeff0*optcoeff0 etc.
                arith_image_cstmult("optcoeff0", acoeff0, "optcoeff0m");
                arith_image_cstmult("optcoeff1", acoeff1, "optcoeff1m");
                arith_image_cstmult("optcoeff2", acoeff2, "optcoeff2m");

                // diagnostic
                FUNC_CHECK_RETURN(save_fl_fits("optcoeff0", "optcoeff0.fits")); //TEST
                FUNC_CHECK_RETURN(save_fl_fits("optcoeff1", "optcoeff1.fits"));
                FUNC_CHECK_RETURN(save_fl_fits("optcoeff2", "optcoeff2.fits"));

                // optcoeff01m = acoeff0*optcoeff0 + acoeff1*optcoeff1
                arith_image_add("optcoeff0m", "optcoeff1m", "optcoeff01m");
                // optcoeff = acoeff0*optcoeff0 + acoeff1*optcoeff1 + acoeff2*optcoeff2
                arith_image_add("optcoeff01m", "optcoeff2m", "optcoeff");
                // optcoeff now has our search direction
                FUNC_CHECK_RETURN(delete_image_ID("optcoeff0m", DELETE_IMAGE_ERRMODE_WARNING));
                FUNC_CHECK_RETURN(delete_image_ID("optcoeff1m", DELETE_IMAGE_ERRMODE_WARNING));
                FUNC_CHECK_RETURN(delete_image_ID("optcoeff2m", DELETE_IMAGE_ERRMODE_WARNING));
                FUNC_CHECK_RETURN(delete_image_ID("optcoeff01m", DELETE_IMAGE_ERRMODE_WARNING));

                ID = image_ID("optcoeff");
                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 21;

                // do linear scan along the direction optcoeff from current parameter location
                int linscanOK = 1;
                // size of step in this direction
                double scangain = 0.0;        //scanstepgain; overriden below
                contrastval = 100000000000.0; // initialize minimization objective to a big number
                //double bestgain = 0.0;
                // iteration counter of steps in the current direction optcoeff
                int k = 0;

                // if(alphascaninit==1)
                //		scangain = bestgain/scangainfact/scangainfact/scangainfact/scangainfact/scangainfact;
                //	alphascaninit = 1;

                // scangain is our location along the direction optcoeff
                scangain = 0.001;
                //              scanstepgain = 0.000001; // TEST
                //              scangainfact = 1.00001; // TEST

                int linoptlimflagarray[100];

                // while objective value < previous value and we've taken no more than than 90 steps
                while (linscanOK == 1)
                {
                    // for state tracking and statistics
                    data.image[IDstatus].array.UI16[0] = 22;

                    // compute offsets
                    ID = image_ID("optcoeff"); // direction vector
                    linoptlimflagarray[k] = 0;
                    // step each parameter by optcoeff
                    for (long i = 0; i < piaacmcparams.linopt_number_param; i++) // looping over parameters
                    {
                        // compute step delta for this parameter
                        // image[ID] = optcoeff is a derivative w.r.t. piaacmcparams.linopt_paramdelta, so has dimension
                        // (parameter dimension)/(piaacmcparams.linopt_paramdelta dimension) so we have to mulitply by
                        // piaacmcparams.linopt_paramdelta to put our step in physical parameter units
                        // negative because we want to cancel the value from the delta PSF
                        piaacmcparams.linopt_paramdeltaval[i] =
                            -scangain * data.image[ID].array.F[i] * piaacmcparams.linopt_paramdelta[i];
                        if (piaacmcparams.linopt_paramdeltaval[i] <
                            -piaacmcparams.linopt_parammaxstep[i]) // if the step is too large in the negative direction
                        {
                            printf("MIN LIMIT [%3ld   %20g]   %20g -> ", i, piaacmcparams.linopt_paramdelta[i],
                                   piaacmcparams.linopt_paramdeltaval[i]); //TEST
                            piaacmcparams.linopt_paramdeltaval[i] =
                                -piaacmcparams.linopt_parammaxstep[i]; // set it to the negative largest allowed step
                            printf(" %20g\n", piaacmcparams.linopt_paramdeltaval[i]); //TEST
                            linoptlimflagarray[k] = 1;
                        }
                        if (piaacmcparams.linopt_paramdeltaval[i] >
                            piaacmcparams.linopt_parammaxstep[i]) // if the step is too large in the positive direction
                        {
                            printf("MAX LIMIT [%3ld   %20g]   %20g -> ", i, piaacmcparams.linopt_paramdelta[i],
                                   piaacmcparams.linopt_paramdeltaval[i]); //TEST
                            piaacmcparams.linopt_paramdeltaval[i] =
                                piaacmcparams.linopt_parammaxstep[i]; // set it to the positive largest allowed step
                            printf(" %20g\n", piaacmcparams.linopt_paramdeltaval[i]); //TEST
                            linoptlimflagarray[k] = 1;
                        }

                        // apply offsets to the global data object via the pointers piaacmcparams.linopt_paramvalf, which
                        // point into the (hopefully) coorect locations of each parameter's value in the
                        // data object
                        if (piaacmcparams.linopt_paramtype[i] == _DATATYPE_FLOAT)
                        {
                            if (*(piaacmcparams.linopt_paramvalf[i]) + (float)piaacmcparams.linopt_paramdeltaval[i] >
                                piaacmcparams.linopt_parammax[i])
                            // if we're about to step too far, set the step to the
                            // limit - current parameter value, so we step to the limit
                            {
                                piaacmcparams.linopt_paramdeltaval[i] =
                                    piaacmcparams.linopt_parammax[i] - *(piaacmcparams.linopt_paramvalf[i]);
                            }

                            if (*(piaacmcparams.linopt_paramvalf[i]) + (float)piaacmcparams.linopt_paramdeltaval[i] <
                                piaacmcparams.linopt_parammin[i])
                            // if we're about to step too far, set the step to the
                            // limit - current parameter value, so we step to the limit
                            {
                                piaacmcparams.linopt_paramdeltaval[i] =
                                    piaacmcparams.linopt_parammin[i] - *(piaacmcparams.linopt_paramvalf[i]);
                            }
                            // take the actual step (piaacmcparams.linopt_paramvalf is a 1D array of pointers)
                            *(piaacmcparams.linopt_paramvalf[i]) += (float)piaacmcparams.linopt_paramdeltaval[i];
                        }
                        else // same for the double case
                        {
                            if (*(piaacmcparams.linopt_paramval[i]) + piaacmcparams.linopt_paramdeltaval[i] >
                                piaacmcparams.linopt_parammax[i])
                            {
                                piaacmcparams.linopt_paramdeltaval[i] =
                                    piaacmcparams.linopt_parammax[i] - *(piaacmcparams.linopt_paramval[i]);
                            }

                            if (*(piaacmcparams.linopt_paramval[i]) + piaacmcparams.linopt_paramdeltaval[i] <
                                piaacmcparams.linopt_parammin[i])
                            {
                                piaacmcparams.linopt_paramdeltaval[i] =
                                    piaacmcparams.linopt_parammin[i] - *(piaacmcparams.linopt_paramval[i]);
                            }

                            *(piaacmcparams.linopt_paramval[i]) += piaacmcparams.linopt_paramdeltaval[i];
                        }
                    }
                    // store the current objective value for later comparison
                    double valold = contrastval;
                    // for state tracking and statistics
                    data.image[IDstatus].array.UI16[0] = 23;

                    // compute new state and compute assossiated evaluation metric
                    // using the modified global data object
                    FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(
                        0.0, 0.0, 0, piaacmcopticalsystem.NBelem, 0, piaacmcparams.computePSF_ResolvedTarget,
                        piaacmcparams.computePSF_ResolvedTarget_mode, 0, &contrastval));
                    double valContrast = contrastval; // contrast component of the evaluation metric
                    // we've now only done the light portion

                    // add regularization component of the evaluation metrix
                    // first, compute and add PIAA shape regularization value (val0) if applicable
                    // this is the same as the previous computation of val0
                    double val0 = 0.0;
                    if (piaacmcparams.linopt_REGPIAASHAPES == 1)
                    {
                        val0 = PIAACMCsimul_regularization_PIAAshapes_value();
                        contrastval += val0;
                    }

                    // add sag regularization (val1) if applicable, as before
                    double val1 = 0.0;
                    if (piaacmcparams.linopt_REGFPMSAG == 1)
                    {
                        val1 = PIAACMCsimul_regularization_fpmsag_value();
                        ;
                        contrastval += val1;
                    }
                    // val is now our complete objective!! Yay!!

                    // for state tracking and statistics
                    data.image[IDstatus].array.UI16[0] = 24;

                    { // print it for monitoring
                        char fname[STRINGMAXLEN_FULLFILENAME];
                        WRITE_FULLFILENAME(fname, "%s/linoptval.txt", piaacmcparams.piaacmcconfdir);
                        FILE *fp = fopen(fname, "a");
                        // printf the first part of the line reporting current values
                        fprintf(fp,
                                "##  [ %5ld / %5ld ]   %5.3f  %12lf  %12g   (reg = %12g [%1d] %12g [%1d]   contrast = "
                                "%20g)       [%d] [%ld]",
                                iter, piaacmcparams.linopt_NBiter, alphareg, scangain, contrastval, val0,
                                piaacmcparams.linopt_REGPIAASHAPES, val1, piaacmcparams.linopt_REGFPMSAG, valContrast,
                                linoptlimflagarray[k], piaacmcparams.linopt_number_param);

                        // now add text indicating status and complete line
                        // and store all parameters for the current best solution
                        if ((contrastval < bestval) || (initbestval == 0))
                        {
                            for (long i = 0; i < piaacmcparams.linopt_number_param; i++)
                                if (piaacmcparams.linopt_paramtype[i] == _DATATYPE_FLOAT)
                                {
                                    data.image[IDoptvec].array.F[i] = *(piaacmcparams.linopt_paramvalf[i]);
                                }
                                else
                                {
                                    data.image[IDoptvec].array.F[i] = (float)*(piaacmcparams.linopt_paramval[i]);
                                }
                            bestval = contrastval;
                            if (initbestval == 0)
                            {
                                fprintf(fp, " ===== START POINT =====\n");
                            }
                            else
                            {
                                fprintf(fp, "  -> BEST VECTOR =======\n");
                            }
                            bestgain = scangain;
                            initbestval = 1;
                        }
                        else
                        {
                            fprintf(fp, " bestval = %12g\n", bestval);
                        }
                        fclose(fp);
                    }

                    // remove offsets returning the global data object to its original state
                    for (long i = 0; i < piaacmcparams.linopt_number_param; i++)
                    {
                        if (piaacmcparams.linopt_paramtype[i] == _DATATYPE_FLOAT)
                        {
                            *(piaacmcparams.linopt_paramvalf[i]) -= (float)piaacmcparams.linopt_paramdeltaval[i];
                        }
                        else
                        {
                            *(piaacmcparams.linopt_paramval[i]) -= piaacmcparams.linopt_paramdeltaval[i];
                        }
                    }
                    // for state tracking and statistics
                    data.image[IDstatus].array.UI16[0] = 25;

                    // store the current position and value
                    //linoptgainarray[k] = scangain;
                    //linoptvalarray[k] = val;
                    k++; // next step

                    // test to see if we're no longer getting better
                    if (contrastval < valold)
                    {
                        linscanOK = 1; // if we're getting better keep going
                        // bestgain = scangain;
                        // scangain += scanstepgain;
                    }
                    else // otherwise stop stepping
                    {
                        linscanOK = 0;
                    }

                    if (k > 90) // stop if we've taken too many steps
                    {
                        linscanOK = 0;
                    }

                    // increment our location along the line
                    double scanstepgain = 0.001;
                    scangain += scanstepgain; // scanstepgain is an initilizaed function local (currently 0.001)
                    scangain *= scangainfact; // causes later steps to be larger
                    // (implicit scangainfact^n for the nth step)
                }
                // NBlinoptgain is counting the largest number of steps needed in this inner loop
                // stepping in the current direction.  When this is small (< 3) we declare victory
                // and stop the outer linear optimization
                if (k > NBlinoptgain)
                {
                    NBlinoptgain = k;
                }

                delete_image_ID("optcoeff", DELETE_IMAGE_ERRMODE_WARNING); // delete the current direction
                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 26;
            }
            // best solution after this linear linescan is stored in IDoptvec
            FUNC_CHECK_RETURN(delete_image_ID("optcoeff0", DELETE_IMAGE_ERRMODE_WARNING));
            FUNC_CHECK_RETURN(delete_image_ID("optcoeff1", DELETE_IMAGE_ERRMODE_WARNING));
            FUNC_CHECK_RETURN(delete_image_ID("DHmodes", DELETE_IMAGE_ERRMODE_WARNING));

            // we've now found the minimum using the three directions from the
            // alternative decompositions of the parameter space with the DHmodes basis
            // (with different conditioning).
            // Now we check the result by recomputing from scratch the objective with the current
            // (hopefully) optimal parameters.  Belts and suspenders

            // At this point the global data object has been restored to its original
            // (non-optimal) state.

            // (re-)compute best solution identified in previous linescan
            // update state to best solution (IDoptvec), setting the global data object
            // to the optimal state
            for (long i = 0; i < piaacmcparams.linopt_number_param; i++)
            {
                if (piaacmcparams.linopt_paramtype[i] == _DATATYPE_FLOAT)
                {
                    *(piaacmcparams.linopt_paramvalf[i]) = data.image[IDoptvec].array.F[i];
                }
                else
                {
                    *(piaacmcparams.linopt_paramval[i]) = (double)data.image[IDoptvec].array.F[i];
                }
            }
            //double valold = contrastval;
            // compute contrast metric component -> val using the data object in the latest optimal state
            FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(0.0, 0.0, 0, piaacmcopticalsystem.NBelem, 0,
                                                      piaacmcparams.computePSF_ResolvedTarget,
                                                      piaacmcparams.computePSF_ResolvedTarget_mode, 0, &contrastval));
            // add PIAA shape regularization component (val0) if applicable
            // same val0 and val1 computations are before, using the latest optimal state
            // double val0 = 0.0;
            if (piaacmcparams.linopt_REGPIAASHAPES == 1)
            {
                double val0 = PIAACMCsimul_regularization_PIAAshapes_value();
                contrastval += val0;
            }

            // add sag regularization component if applicable
            if (piaacmcparams.linopt_REGFPMSAG == 1)
            {
                double val1 = PIAACMCsimul_regularization_fpmsag_value();
                ;
                contrastval += val1;
            }

            // now val is the objective including any desired regularization terms using the
            // latest optimal solution

            printf("gain: %lf -> val = %20g\n", bestgain, contrastval);
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 27;

            // update reference state evaluation vector of optimal results including evaluation zone values
            // and desired regularization terms
            // which sets up starting the next iteration at the best solution

            ID1Dref = image_ID("vecDHref1D");
            ID = image_ID("imvect");
            // first fill in evaluation zone (complex) values
            {
                uint64_t ii;
                for (ii = 0; ii < data.image[ID].md[0].nelement; ii++)
                {
                    data.image[ID1Dref].array.F[ii] = data.image[ID].array.F[ii];
                }
                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 28;

                // now fill in the regularization terms if desired
                // same code as before.
                if (piaacmcparams.linopt_REGPIAASHAPES == 1)
                {
                    ii = PIAACMCsimul_regularization_PIAAshapes_add1Dvector(ID1Dref, ii);
                }

                if (piaacmcparams.linopt_REGFPMSAG == 1)
                {
                    ID = piaacmcopticaldesign.zonezID;
                    for (uint32_t jj = 0; jj < data.image[ID].md[0].size[0]; jj++)
                    {
                        data.image[ID1Dref].array.F[ii] =
                            pow(data.image[ID].array.D[jj] / piaacmcopticaldesign.fpmsagreg_coeff,
                                piaacmcopticaldesign.fpmsagreg_alpha);
                        ii++;
                    }
                }
            }

            FUNC_CHECK_RETURN(delete_image_ID("imvect", DELETE_IMAGE_ERRMODE_WARNING));
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 29;

            { // print out current best value for tracking
                char fname[STRINGMAXLEN_FULLFILENAME];
                WRITE_FULLFILENAME(fname, "%s/linoptval.txt", piaacmcparams.piaacmcconfdir);
                FILE *fp = fopen(fname, "a");
                if (fp == NULL)
                {
                    printf("ERROR: cannot open file \"%s\"\n", fname);
                    exit(0);
                }
                fprintf(fp, "-> %5ld    %20g <- %20g \n", iter, contrastval, valref);
                printf("%5ld %20g %20g \n", iter, contrastval, valref);
                fflush(stdout);
                fclose(fp);
            }

            // save current best value and reference value in globals
            piaacmcparams.PIAACMCSIMUL_VAL = contrastval;
            piaacmcparams.PIAACMCSIMUL_VALREF = valref;

            // Nominally if we're in this linear
            // optimization piaacmcparams.PIAACMC_fpmtype = 1, so the next line is not executed
            if (piaacmcparams.PIAACMC_fpmtype == 0) // in the idealized PIAACMC case
            {
                piaacmcopticaldesign.fpmaskamptransm =
                    data.image[piaacmcopticaldesign.zoneaID].array.D
                        [0]; // required to ensure that the new optimal focal plane mask transmission is written to disk
            }

            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 30;

            { // tracking diagnostics giving behavior of the modes by iteration
                char dirname[STRINGMAXLEN_DIRNAME];
                WRITE_DIRNAME(dirname, "%s_linopt", piaacmcparams.piaacmcconfdir);

                PIAACMCsimul_savepiaacmcconf(dirname); // staging area
                EXECUTE_SYSTEM_COMMAND("rsync -au --progress %s/* ./%s/", dirname, piaacmcparams.piaacmcconfdir);

                EXECUTE_SYSTEM_COMMAND("cp %s/piaa0Cmodes.fits %s/piaa0Cmodes.%04ld.fits", dirname, dirname, iter);
                EXECUTE_SYSTEM_COMMAND("cp %s/piaa0Fmodes.fits %s/piaa0Fmodes.%04ld.fits", dirname, dirname, iter);

                EXECUTE_SYSTEM_COMMAND("cp %s/piaa1Cmodes.fits %s/piaa1Cmodes.%04ld.fits", dirname, dirname, iter);
                EXECUTE_SYSTEM_COMMAND("cp %s/piaa1Fmodes.fits %s/piaa1Fmodes.%04ld.fits", dirname, dirname, iter);
                EXECUTE_SYSTEM_COMMAND("cp %s/piaacmcparams.conf %s/piaacmcparams.%04ld.conf", dirname, dirname, iter);
            }

            if (file_exists(stopfile) == 1)
            {
                iterOK = 0;
            }

            // Figure out if current loop should continue optimization
            // if optimization ends, then iterOK set to 0
            // if we've reached the allowed number of iterations
            if (iter == piaacmcparams.linopt_NBiter)
            {
                iterOK = 0;
            }
            if (iter > 2)
            {
                if (contrastval > 0.98 * oldval) // if after second iteration and we'be improved by less than 10%
                {
                    iterOK = 0;
                }
            }

            if (NBlinoptgain < 3) // if we've stopped moving much
            {
                iterOK = 0;
            }

            // set up for next iteration
            oldval = contrastval;
            iter++;

            printf("END OF LOOP ITERATION\n");
            fflush(stdout);
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 31;
        }
        printf(" ============ END OF OPTIMIZATION LOOP ======= \n");
        // for state tracking and statistics
        data.image[IDstatus].array.UI16[0] = 32;
    } // end of if (piaacmcparams.LINOPT==1): done with the linear optimization

    piaacmcparams.LINOPT = 0;

    //  PIAACMCsimul_savepiaacmcconf("piaacmc0");
    //  PIAACMCsimul_loadpiaacmcconf("piaacmc0");
    // PIAACMCsimul_savepiaacmcconf("piaacmc1");
    //exit(0);

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
