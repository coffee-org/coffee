/**
 * @file    PIAAACMCsimul_initpiaacmcconf.c
 * @brief   PIAA-type coronagraph design, initialize
 *
 * Can design both APLCMC and PIAACMC coronagraphs
 *
 */

// log all debug trace points to file
#define DEBUGLOG

// System includes

#include <assert.h>
#include <ctype.h>
#include <malloc.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>

#include <time.h>

// External libraries
#include <fitsio.h>

// milk includes
#include "CommandLineInterface/CLIcore.h"
//   core modules
#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_tools/COREMOD_tools.h"

//   other modules
#include "OpticsMaterials/OpticsMaterials.h"
#include "WFpropagate/WFpropagate.h"
#include "coronagraphs/coronagraphs.h"
#include "fft/fft.h"
#include "image_basic/image_basic.h"
#include "image_filter/image_filter.h"
#include "image_gen/image_gen.h"
#include "info/info.h"
#include "linopt_imtools/linopt_imtools.h"
#include "statistic/statistic.h"

#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"

#include "PIAACMCsimul_loadsavepiaacmcconf.h"

#include "FocalPlaneMask/mkFPM_zonemap.h"

#include "LyotStop/mkSimpleLyotStop.h"

#include "PIAAshape/init_geomPIAA_rad.h"
#include "PIAAshape/load2DRadialApodization.h"
#include "PIAAshape/mkPIAAMshapes_from_RadSag.h"

#include "init_piaacmcopticaldesign.h"

#ifdef HAVE_LIBGOMP
#include <omp.h>
#define OMP_NELEMENT_LIMIT 1000000
#endif

// has the configuration been loaded ?
static int loaded = 0;

// if 1, save conf at end of this function
static int saveconf = 0;

/**
 *
 * # List of Configuration Parameters {#page_PIAACMCsimul_Configuration}
 *
 * Each parameter is stored on disk as conf/conf_<param>.txt
 *
 *
 *
 *
 * PIAA OPTICS DESIGN:
 *
 * | Parameter              | Meaning                                                  |
 * | ---------------------- | -------------------------------------------------------- |
 * | PIAAmode               | PIAA mode (0: classical apodization, 1: PIAA)            |
 * | fpmradld               | Focal plane mask radius [l/D]                            |
 * | PIAAcoeff              | Fraction of PIAA apodization                             |
 * | coin                   | Central obstruction at input beam [beam radius]          |
 * | coout                  | Central obstruction at output beam [beam radius]
 * | PIAAmaterial           | PIAA optics material
 * | PIAAcirc               | FLAG: 1 if PIAA shapes are circular (no Fourier modes)
 * | REGPIAACCOEFF          | PIAA regularization amplitude, cosine modes
 * | REGPIAACALPHA          | PIAA regularization power law, cosine modes
 * | REGPIAAFCOEFF          | PIAA regularization amplitude, Fourier modes
 * | REGPIAAFALPHA          | PIAA regularization power law, Fourier modes
 *
 *
 * LYOT STOP(S) DESIGN:
 *
 * | Parameter              | Meaning                                                  |
 * | ---------------------- | -------------------------------------------------------- |
 * | LStransm               | Lyot stop transmission
 * | NBls                   | Number of Lyot stops
 * | lambda                 | Wavelength for monochromatic design [nm]
 *
 *
 *
 *
 * FOCAL PLANE MASK DESIGN:
 *
 *
 * | Parameter              | Meaning                                                  |
 * | ---------------------- | -------------------------------------------------------- |
 * | fpmmaterial            | focal plane mask material
 * | FPMsectors             | mask geometry: 0=disk, 1=sectors, 2=hexagonal tiling
 * | NBrings                | number of rings in focal plane mask
 * | maskradld              | mask outer radius at central wavelength [l/D]
 * | fpmminsag              | min focal plane mask sag
 * | fpmmaxsag              | max focal plane mask sag
 * | fpmregsag_coeff        | sag regularization coefficient
 * | fpmregsag_alpha        | sag regularization coefficient exponent
 * | fpmccnbr               | how many central rings replaced by cone (set to 0 if no central cone
 * | fpmccz                 | sag at cone center (sag at cone edge will be midpoint between minsag and maxsag)
 * | fpmocradld             | outer cone outer radius [l/D]
 * | fpmocz                 | sag at inner edge of outer cone (sag = 0 at outer edge), set to 0 if no outer cone
 *
 *
 *
 * OPTIMIZATION PARAMETERS:
 *
 * | Parameter              | Meaning                                                  |
 * | ---------------------- | -------------------------------------------------------- |
 * | mlambda                | central wavelength for polychromatic design [nm]
 * | mlambdaB               | spectral bandwidth [%]
 * | nblambda               | Number of wavelength values
 * | ssize                  | source angular size for mask optimization (20: 0.01 l/D; 10: 0.1 l/D)
 * | extmode                | source extent mode (0: 1 point, 1: 3 points; 2: 6 points)
 *
 *
 *
 * OPTICAL DESIGN:
 *
 *
 * | Parameter              | Meaning                                                  |
 * | ---------------------- | -------------------------------------------------------- |
 * | size                   | array size
 * | beamrad                | beam radius [mm]
 * | pscale                 | pixel scale in pupil [m/pix]
 * | Fratio                 | F ratio at focal plane mask
 * | PIAAr0lim              | outer edge of PIAA optic 0 [beam radius unit]
 * | PIAAr1lim              | outer edge of PIAA optic 1 [beam radius unit]
 * | PIAAsep                | distance between PIAA optics [m]
 * | PIAA0pos               | PIAA optic 0 distance from pupil plane [m]
 * | invPIAAmode            | 0: no inv PIAA, 1: inv PIAA after Lyot stops, 2: inv PIAA before Lyot stops
 * | prePIAA0maskpos        | pre-PIAA optic 0 mask distance from pupil plane [m] (if prePIAA0mask.fits exists)
 * | postPIAA0maskpos       | post-PIAA optic 0 mask distance from pupil plane [m] (if postPIAA0mask.fits exits)
 * | piaaNBCmodesmax        | maximum number of radial cosine modes for PIAA optics
 * | piaaCPAmax             | maximum spatial frequency (CPA) for PIAA optics
 * | LyotZmin               | minimum value for Lyot stop(s) conjugation range [m] - relative to element named "post focal plane mask pupil"
 * | LyotZmax               | maximum value for Lyot stop(s) conjugation range [m] - relative to element named "post focal plane mask pupil"
 * | pupoutmaskrad          | output pupil mask radius (scaled to pupil radius)
 *
 *
 */

static errno_t PIAACMCsimul_initpiaacmcconf_readconfparams(long   piaacmctype,
                                                           double fpmradld,
                                                           double centobs0,
                                                           double centobs1,
                                                           int    WFCmode)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %ld %lf %lf %lf %d",
                     piaacmctype,
                     fpmradld,
                     centobs0,
                     centobs1,
                     WFCmode);

    // Default Values for PIAACMC (will adopt them unless configuration file exists)
    piaacmcopticaldesign.nblambda = 8;

    // high resolution
    //piaacmcopticaldesign.size = 4096;
    //piaacmcopticaldesign.pixscale = 0.000055;

    // mid resolution
    piaacmcopticaldesign.size     = 2048;
    piaacmcopticaldesign.pixscale = 0.000055;

    // low resolution
    //piaacmcopticaldesign.size = 1024;
    //piaacmcopticaldesign.pixscale = 0.00011;

    // very low resolution
    //    piaacmcopticaldesign.size = 512;
    // piaacmcopticaldesign.pixscale = 0.00022;

    piaacmcopticaldesign.beamrad   = 0.01; // beam physical radius
    piaacmcopticaldesign.PIAA0pos  = 1.0;  // piaa 0 position [m]
    piaacmcopticaldesign.PIAAsep   = 1.00; // [m]
    piaacmcopticaldesign.fpzfactor = 8.0;
    piaacmcopticaldesign.Fratio    = 80.0;                    // default
    strcpy(piaacmcopticaldesign.PIAAmaterial_name, "Mirror"); // mirrors
    piaacmcopticaldesign.prePIAA0mask     = 0;
    piaacmcopticaldesign.prePIAA0maskpos  = 0.0;
    piaacmcopticaldesign.postPIAA0mask    = 0;
    piaacmcopticaldesign.postPIAA0maskpos = 0.0;
    piaacmcopticaldesign.piaaNBCmodesmax  = 40;
    piaacmcopticaldesign.piaaCPAmax       = 10.0;

    piaacmcopticaldesign.centObs0 = centobs0; // input central obstruction
    piaacmcopticaldesign.centObs1 = centobs1; // output central obstruction
    piaacmcopticaldesign.NBradpts = 50000;
    piaacmcopticaldesign.r0lim =
        1.15; // outer radius after extrapolation, piaa optics 0
    piaacmcopticaldesign.r1lim =
        1.5; // outer radius after extrapolation, piaa optics 1

    /// Wavefront control
    piaacmcopticaldesign.nbDM =
        WFCmode; // number of deformable mirrors (10 max)
    for (long iDM = 0; iDM < piaacmcopticaldesign.nbDM; iDM++)
    {
        piaacmcopticaldesign.DMpos[iDM] =
            0.0 + 0.6 * iDM /
                      (0.01 + piaacmcopticaldesign.nbDM -
                       1.0); // DM conjugation in collimated space
        piaacmcopticaldesign.ID_DM[iDM] =
            -1; // DM image identifier - to be updated later
    }

    piaacmcopticaldesign.NBLyotStop = 2;
    for (int i = 0; i < 10; i++)
    {
        piaacmcopticaldesign.LyotStop_zpos[i] = 0.0;
        piaacmcopticaldesign.IDLyotStop[i]    = -1;
    }

    piaacmcopticaldesign.fpmaskradld =
        fpmradld; // to compute prolate spheroidal function
    piaacmcopticaldesign.fpmarraysize = 2048;

    piaacmcopticaldesign.fpmRad    = 100.0e-6; // focal plane radius [m]
    piaacmcopticaldesign.NBrings   = 4; // number of rings in focal plane mask
    piaacmcopticaldesign.fpmminsag = -1.0e-5;
    piaacmcopticaldesign.fpmmaxsag = 1.0e-5;
    piaacmcopticaldesign.fpmsagreg_coeff   = 1.0;
    piaacmcopticaldesign.fpmsagreg_alpha   = 1.0;
    piaacmcopticaldesign.NBringCentCone    = 0; // central cone
    piaacmcopticaldesign.fpmCentConeZ      = piaacmcopticaldesign.fpmmaxsag;
    piaacmcopticaldesign.fpmOuterConeZ     = piaacmcopticaldesign.fpmmaxsag;
    piaacmcopticaldesign.fpmOuterConeRadld = 80.0;
    piaacmcopticaldesign.fpmmaterial_code  = 0; // 0: mirror
    piaacmcopticaldesign.fpmaskamptransm   = 1.0;

    {
        FILE *fp;
        if ((fp = fopen("conf/conf_peakPSF.txt", "r")) != NULL)
        {
            DEBUG_TRACEPOINT("File conf/conf_peakPSF.txt open for reading");
            if (fscanf(fp, "%f", &piaacmcopticaldesign.peakPSF) != 1)
            {
                FUNC_RETURN_FAILURE(
                    "Cannot read value from file conf/conf_peakPSF.txt");
            }
            fclose(fp);
        }
        else
        {
            piaacmcopticaldesign.peakPSF = -1.0;
        }
        DEBUG_TRACEPOINT("piaacmcopticaldesign.peakPSF = %f",
                         piaacmcopticaldesign.peakPSF);
    }

    piaacmcopticaldesign.PIAAmode = 1;
    {
        variableID IDv = variable_ID("PIAACMC_PIAAmode");
        if (IDv != -1)
        {
            piaacmcopticaldesign.PIAAmode =
                (int) (data.variable[IDv].value.f + 0.01);
        }
    }
    DEBUG_TRACEPOINT("piaacmcopticaldesign.PIAAmode = %d",
                     piaacmcopticaldesign.PIAAmode);

    piaacmcopticaldesign.PIAAcoeff = 1.0;
    {
        variableID IDv = variable_ID("PIAACMC_PIAAcoeff");
        if (IDv != -1)
        {
            piaacmcopticaldesign.PIAAcoeff = data.variable[IDv].value.f;
        }
    }
    DEBUG_TRACEPOINT("piaacmcopticaldesign.PIAAcoeff = %f",
                     piaacmcopticaldesign.PIAAcoeff);

    if (piaacmcopticaldesign.PIAAmode == 0)
    {
        piaacmcopticaldesign.invPIAAmode = 0;
    }
    else
    {
        piaacmcopticaldesign.invPIAAmode = 1;
        {
            variableID IDv = variable_ID("PIAACMC_invPIAAmode");
            if (IDv != -1)
            {
                piaacmcopticaldesign.invPIAAmode =
                    (long) (data.variable[IDv].value.f + 0.001);
            }
        }
    }
    DEBUG_TRACEPOINT("piaacmcopticaldesign.invPIAAmode = %d",
                     piaacmcopticaldesign.invPIAAmode);

    // FOCAL PLANE MATERIAL

    { // read focal plane mask material
        // -> piaacmcopticaldesign.fpmmaterial_name
        FILE *fp;
        char  fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname,
                           "%s/conf_fpmmaterial_name.txt",
                           piaacmcparams.piaacmcconfdir);
        if ((fp = fopen(fname, "r")) != NULL)
        {
            char name[200];
            if (fscanf(fp, "%199s", name) == 1)
            {
                strcpy(piaacmcopticaldesign.fpmmaterial_name, name);
                printf(
                    "Reading %s   piaacmcopticaldesign.fpmmaterial_name : %s\n",
                    fname,
                    piaacmcopticaldesign.fpmmaterial_name);
            }
            else
            {
                FUNC_RETURN_FAILURE("Cannot read value from file %s", fname);
            }
            fclose(fp);
        }
        else
        {
            snprintf(piaacmcopticaldesign.fpmmaterial_name,
                     STRINGMAXLEN_PIAACMCSIMUL_MATERIALNAME,
                     "Mirror");
            WRITE_FULLFILENAME(fname,
                               "%s/conf_fpmmaterial_name.txt",
                               piaacmcparams.piaacmcconfdir);
            printf("Writing %s   piaacmcopticaldesign.fpmmaterial_name : %s\n",
                   fname,
                   piaacmcopticaldesign.fpmmaterial_name);
            if ((fp = fopen(fname, "w")) != NULL)
            {
                fprintf(fp, "%s\n", piaacmcopticaldesign.fpmmaterial_name);
                fclose(fp);
            }
            else
            {
                FUNC_RETURN_FAILURE("Cannot create file \"%s\"", fname);
            }
        }
    }

    printf("piaacmcopticaldesign.fpmmaterial_name : %s\n",
           piaacmcopticaldesign.fpmmaterial_name);
    piaacmcopticaldesign.fpmmaterial_code =
        OpticsMaterials_code(piaacmcopticaldesign.fpmmaterial_name);

    { // write focal plane mask material
        FILE *fp;
        char  fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname,
                           "%s/conf_fpmmaterial_code.txt",
                           piaacmcparams.piaacmcconfdir);
        if ((fp = fopen(fname, "w")) != NULL)
        {
            fprintf(fp, "%d\n", piaacmcopticaldesign.fpmmaterial_code);
            fclose(fp);
        }
        else
        {
            FUNC_RETURN_FAILURE("Cannot create file \"%s\"", fname);
        }
    }
    DEBUG_TRACEPOINT("piaacmcopticaldesign.fpmmaterial_code = %d",
                     piaacmcopticaldesign.fpmmaterial_code);

    { // read input parameters

        variableID IDv;
        if ((IDv = variable_ID("PIAACMC_beamrad")) != -1)
        {
            piaacmcopticaldesign.beamrad =
                data.variable[IDv].value.f; // beam physical radius
        }

        if ((IDv = variable_ID("PIAACMC_Fratio")) != -1)
        {
            piaacmcopticaldesign.Fratio =
                data.variable[IDv].value.f; // Focal ratio
        }
        if ((IDv = variable_ID("PIAACMC_r0lim")) != -1)
        {
            piaacmcopticaldesign.r0lim = data.variable[IDv].value.f;
        }
        if ((IDv = variable_ID("PIAACMC_r1lim")) != -1)
        {
            piaacmcopticaldesign.r1lim = data.variable[IDv].value.f;
        }

        if ((IDv = variable_ID("PIAACMC_PIAAsep")) != -1)
        {
            piaacmcopticaldesign.PIAAsep =
                data.variable[IDv].value.f; // piaa separation
        }
        if ((IDv = variable_ID("PIAACMC_PIAA0pos")) != -1)
        {
            piaacmcopticaldesign.PIAA0pos =
                data.variable[IDv].value.f; // piaa elem 0 position
        }

        if ((IDv = variable_ID("PIAACMC_prePIAA0maskpos")) != -1)
        {
            piaacmcopticaldesign.prePIAA0maskpos =
                data.variable[IDv].value.f; // pre piaa elem 0 mask position
        }
        if ((IDv = variable_ID("PIAACMC_postPIAA0maskpos")) != -1)
        {
            piaacmcopticaldesign.postPIAA0maskpos =
                data.variable[IDv].value.f; // post piaa elem 0 mask position
        }

        piaacmcopticaldesign.LyotZmin = -3.0;
        if ((IDv = variable_ID("PIAACMC_LyotZmin")) != -1)
        {
            piaacmcopticaldesign.LyotZmin = data.variable[IDv].value.f;
        }
        piaacmcopticaldesign.LyotZmax = 3.0;
        if ((IDv = variable_ID("PIAACMC_LyotZmax")) != -1)
        {
            piaacmcopticaldesign.LyotZmax = data.variable[IDv].value.f;
        }

        piaacmcopticaldesign.pupoutmaskrad = 0.95;
        if ((IDv = variable_ID("PIAACMC_pupoutmaskrad")) != -1)
        {
            piaacmcopticaldesign.pupoutmaskrad = data.variable[IDv].value.f;
        }

        if ((IDv = variable_ID("PIAACMC_piaaNBCmodesmax")) != -1)
        {
            piaacmcopticaldesign.piaaNBCmodesmax =
                (long) (data.variable[IDv].value.f +
                        0.01); // max number of Cosine terms
        }
        if ((IDv = variable_ID("PIAACMC_piaaCPAmax")) != -1)
        {
            piaacmcopticaldesign.piaaCPAmax =
                data.variable[IDv].value.f; // max CPA for PIAA shapes tuning
        }

        piaacmcopticaldesign.NBLyotStop = 1;
        if (piaacmcopticaldesign.PIAAmode == 0)
        {
            piaacmcopticaldesign.NBLyotStop = 1;
        }
        else
        {
            if ((IDv = variable_ID("PIAACMC_nblstop")) != -1)
            {
                piaacmcopticaldesign.NBLyotStop =
                    (long) data.variable[IDv].value.f + 0.01;
            }
        }

        if ((IDv = variable_ID("PIAACMC_lambda")) != -1)
        {
            piaacmcopticaldesign.lambda =
                1.0e-9 * data.variable[IDv].value.f; // central wavelength [m]
        }
        //             printf("lambda = %g\n", piaacmcopticaldesign.lambda);

        if ((IDv = variable_ID("PIAACMC_lambdaB")) != -1)
        {
            piaacmcopticaldesign.lambdaB =
                data.variable[IDv].value.f; // spectral bandwidth [%]
        }

        piaacmcparams.LAMBDASTART =
            piaacmcopticaldesign.lambda *
            (1.0 - 0.005 * piaacmcopticaldesign.lambdaB);
        piaacmcparams.LAMBDAEND = piaacmcopticaldesign.lambda *
                                  (1.0 + 0.005 * piaacmcopticaldesign.lambdaB);

        if ((IDv = variable_ID("PIAACMC_nblambda")) != -1)
        {
            piaacmcopticaldesign.nblambda = data.variable[IDv].value.f;
        }

        if ((IDv = variable_ID("PIAACMC_NBrings")) != -1)
        {
            piaacmcopticaldesign.NBrings = data.variable[IDv].value.f;
        }

        if ((IDv = variable_ID("PIAACMC_fpmminsag")) != -1)
        {
            piaacmcopticaldesign.fpmminsag = data.variable[IDv].value.f;
        }
        if ((IDv = variable_ID("PIAACMC_fpmmaxsag")) != -1)
        {
            piaacmcopticaldesign.fpmmaxsag = data.variable[IDv].value.f;
        }

        if ((IDv = variable_ID("PIAACMC_fpmsagreg_coeff")) != -1)
        {
            piaacmcopticaldesign.fpmsagreg_coeff = data.variable[IDv].value.f;
        }
        if ((IDv = variable_ID("PIAACMC_fpmsagreg_alpha")) != -1)
        {
            piaacmcopticaldesign.fpmsagreg_alpha = data.variable[IDv].value.f;
        }

        if ((IDv = variable_ID("PIAACMC_NBringCentCone")) != -1)
        {
            piaacmcopticaldesign.NBringCentCone = data.variable[IDv].value.f;
        }

        if ((IDv = variable_ID("PIAACMC_fpmCentConeZ")) != -1)
        {
            piaacmcopticaldesign.fpmCentConeZ = data.variable[IDv].value.f;
        }
        if ((IDv = variable_ID("PIAACMC_fpmOuterConeZ")) != -1)
        {
            piaacmcopticaldesign.fpmOuterConeZ = data.variable[IDv].value.f;
        }
        if ((IDv = variable_ID("PIAACMC_fpmOuterConeRadld")) != -1)
        {
            piaacmcopticaldesign.fpmOuterConeRadld = data.variable[IDv].value.f;
        }
        piaacmcopticaldesign.fpmOuterConeRad =
            0.5 * (piaacmcparams.LAMBDASTART + piaacmcparams.LAMBDAEND) *
            piaacmcopticaldesign.Fratio *
            piaacmcopticaldesign.fpmOuterConeRadld; // [l/D] radius

        if ((IDv = variable_ID("PIAACMC_size")) != -1)
        {
            piaacmcopticaldesign.size =
                (long) (data.variable[IDv].value.f + 0.01);
        }

        if ((IDv = variable_ID("PIAACMC_pixscale")) != -1)
        {
            piaacmcopticaldesign.pixscale = data.variable[IDv].value.f;
        }
    }

    DEBUG_TRACEPOINT("piaacmctype = %ld", piaacmctype);

    if (piaacmctype == 0) // idealized focal plane mask
    {
        piaacmcparams.FORCE_CREATE_fpmzt =
            1; // force making the focal plane mask
        piaacmcopticaldesign.NBrings         = 1;
        piaacmcopticaldesign.NBringCentCone  = 0;
        piaacmcopticaldesign.fpmOuterConeZ   = 0.0;
        piaacmcopticaldesign.fpmminsag       = -1e-5;
        piaacmcopticaldesign.fpmmaxsag       = 1e-5;
        piaacmcopticaldesign.fpmsagreg_coeff = 1.0;
        piaacmcopticaldesign.fpmsagreg_alpha = 1.0;
        piaacmcopticaldesign.fpmRad =
            0.5 * (piaacmcparams.LAMBDASTART + piaacmcparams.LAMBDAEND) *
            piaacmcopticaldesign.Fratio * fpmradld; // [l/D] radius
        printf(
            "Idealized focal plane mask  radius = %f l/D  = %g m    [lambda = "
            "%g (%g - %g)] [Fratio = %f]\n",
            fpmradld,
            piaacmcopticaldesign.fpmRad,
            0.5 * (piaacmcparams.LAMBDASTART + piaacmcparams.LAMBDAEND),
            piaacmcparams.LAMBDASTART,
            piaacmcparams.LAMBDAEND,
            piaacmcopticaldesign.Fratio);
    }
    else
    {
        if (piaacmcparams.PIAACMC_MASKRADLD < 0.2) // not initialized
        {
            // 1.2x nominal radius, rounded to nearest 0.1 l/D
            piaacmcparams.PIAACMC_MASKRADLD =
                0.1 * ((long) (10.0 * 1.2 * fpmradld));
        }

        piaacmcopticaldesign.fpmRad =
            0.5 * (piaacmcparams.LAMBDASTART + piaacmcparams.LAMBDAEND) *
            piaacmcopticaldesign.Fratio * piaacmcparams.PIAACMC_MASKRADLD;
        printf(
            "Physical focal plane mask - rad = %f l/D -> %g m   [lambda = %g "
            "(%g - %g)] [Fratio = %f]\n",
            piaacmcparams.PIAACMC_MASKRADLD,
            piaacmcopticaldesign.fpmRad,
            0.5 * (piaacmcparams.LAMBDASTART + piaacmcparams.LAMBDAEND),
            piaacmcparams.LAMBDASTART,
            piaacmcparams.LAMBDAEND,
            piaacmcopticaldesign.Fratio);
    }

    piaacmcopticaldesign.CmodesID      = -1; // Cosine radial mode
    piaacmcopticaldesign.FmodesID      = -1; // Fourier 2D modes
    piaacmcopticaldesign.piaa0CmodesID = -1;
    piaacmcopticaldesign.piaa0FmodesID = -1;
    piaacmcopticaldesign.piaa1CmodesID = -1;
    piaacmcopticaldesign.piaa1FmodesID = -1;

    // focm zone material thickness, double precision image
    piaacmcopticaldesign.zonezID = -1;

    // focm zone amplitude transmission, double precision image
    piaacmcopticaldesign.zoneaID = -1;

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

static errno_t make_x_y_r_PA_images(long size, float beamrad)
{
    DEBUG_TRACE_FSTART();

    imageID IDx;
    FUNC_CHECK_RETURN(create_2Dimage_ID("xcoord", size, size, &IDx));

    imageID IDy;
    FUNC_CHECK_RETURN(create_2Dimage_ID("ycoord", size, size, &IDy));

    imageID IDr;
    FUNC_CHECK_RETURN(create_2Dimage_ID("rcoord", size, size, &IDr));

    imageID IDPA;
    FUNC_CHECK_RETURN(create_2Dimage_ID("PAcoord", size, size, &IDPA));

    printf("pre-computing x, y, r, and PA\n");
    fflush(stdout);
    // list_image_ID();

    for (uint32_t ii = 0; ii < size; ii++)
    {
        for (uint32_t jj = 0; jj < size; jj++)
        {
            float x = (1.0 * ii - 0.5 * size) / beamrad;
            float y = (1.0 * jj - 0.5 * size) / beamrad;
            data.image[IDx].array.F[jj * size + ii]  = x;
            data.image[IDy].array.F[jj * size + ii]  = y;
            data.image[IDr].array.F[jj * size + ii]  = sqrt(x * x + y * y);
            data.image[IDPA].array.F[jj * size + ii] = atan2(y, x);
        }
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

static errno_t make_C_F_modes(long size, float beamrad)
{
    DEBUG_TRACE_FSTART();

    printf("Creating / loading Cmodes and Fmodes ...\n");
    fflush(stdout);

    { // Cmodes
        piaacmcparams.CREATE_Cmodes = 0;
        //   sprintf(fname, "%s/Cmodes.fits", piaacmcparams.piaacmcconfdir);

        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname,
                           "ref/Cmodes_%ld.fits",
                           piaacmcopticaldesign.size);

        DEBUG_TRACEPOINT("piaacmcparams.FORCE_CREATE_Cmodes = %d",
                         piaacmcparams.FORCE_CREATE_Cmodes);

        if (piaacmcparams.FORCE_CREATE_Cmodes == 0)
        {
            piaacmcopticaldesign.CmodesID = image_ID("Cmodes");
            if (piaacmcopticaldesign.CmodesID == -1)
            {
                load_fits(fname, "Cmodes", 0, &(piaacmcopticaldesign.CmodesID));
            }
            if (piaacmcopticaldesign.CmodesID == -1)
            {
                piaacmcparams.CREATE_Cmodes = 1;
            }
        }
        else
        {
            piaacmcparams.CREATE_Cmodes = 1;
        }

        DEBUG_TRACEPOINT("piaacmcparams.FORCE_CREATE_Cmodes = %d",
                         piaacmcparams.FORCE_CREATE_Cmodes);

        if (piaacmcparams.CREATE_Cmodes == 1)
        {
            if (piaacmcopticaldesign.CmodesID != -1)
            {
                delete_image_ID("Cmodes", DELETE_IMAGE_ERRMODE_WARNING);
            }
            long Cmsize = (long) (beamrad * 4);
            if (Cmsize > size)
            {
                Cmsize = size;
            }
            printf("beamradpix = %f -> Cmsize = %ld\n", beamrad, Cmsize);
            // make sure Cmsize if even
            if (Cmsize % 2 == 1)
            {
                Cmsize++;
            }
            piaacmcopticaldesign.Cmsize = Cmsize;
            DEBUG_TRACEPOINT("run linopt_imtools_makeCosRadModes");
            linopt_imtools_makeCosRadModes("Cmodes",
                                           Cmsize,
                                           piaacmcopticaldesign.piaaNBCmodesmax,
                                           ApoFitCosFact * beamrad,
                                           2.0,
                                           NULL);
            piaacmcopticaldesign.CmodesID = image_ID("Cmodes");
            save_fits("Cmodes", fname);

            EXECUTE_SYSTEM_COMMAND("mv ModesExpr_CosRad.txt %s/",
                                   piaacmcparams.piaacmcconfdir);
        }
        piaacmcopticaldesign.NBCmodes =
            data.image[piaacmcopticaldesign.CmodesID].md[0].size[2];
        piaacmcopticaldesign.Cmsize =
            data.image[piaacmcopticaldesign.CmodesID].md[0].size[0];
    }

    { // Fmodes
        piaacmcparams.CREATE_Fmodes = 0;
        //    sprintf(fname, "%s/Fmodes.fits", piaacmcparams.piaacmcconfdir);

        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname,
                           "ref/Fmodes_%ld.fits",
                           piaacmcopticaldesign.size);

        DEBUG_TRACEPOINT("piaacmcparams.FORCE_CREATE_Fmodes = %d",
                         piaacmcparams.FORCE_CREATE_Fmodes);

        if (piaacmcparams.FORCE_CREATE_Fmodes == 0)
        {
            piaacmcopticaldesign.FmodesID = image_ID("Fmodes");
            if (piaacmcopticaldesign.FmodesID == -1)
            {
                load_fits(fname,
                          "Fmodes",
                          LOADFITS_ERRMODE_IGNORE,
                          &(piaacmcopticaldesign.FmodesID));
            }
            if (piaacmcopticaldesign.FmodesID == -1)
            {
                piaacmcparams.CREATE_Fmodes = 1;
            }
        }
        else
        {
            piaacmcparams.CREATE_Fmodes = 1;
        }

        DEBUG_TRACEPOINT("piaacmcparams.FORCE_CREATE_Fmodes = %d",
                         piaacmcparams.FORCE_CREATE_Fmodes);

        if (piaacmcparams.CREATE_Fmodes == 1)
        {
            long Fmsize = (long) (beamrad * 4);
            if (Fmsize > size)
            {
                Fmsize = size;
            }
            piaacmcopticaldesign.Fmsize = Fmsize;
            DEBUG_TRACEPOINT("run linopt_imtools_makeCPAmodes");
            FUNC_CHECK_RETURN(
                linopt_imtools_makeCPAmodes("Fmodes",
                                            Fmsize,
                                            piaacmcopticaldesign.piaaCPAmax,
                                            0.8,
                                            beamrad,
                                            2.0,
                                            1,
                                            NULL));
            piaacmcopticaldesign.FmodesID = image_ID("Fmodes");

            FUNC_CHECK_RETURN(save_fits("Fmodes", fname));

            FUNC_CHECK_RETURN(save_fits("cpamodesfreq", "cpamodesfreq.fits"));

            EXECUTE_SYSTEM_COMMAND("mv ModesExpr_CPA.txt %s/",
                                   piaacmcparams.piaacmcconfdir);
        }
        piaacmcopticaldesign.NBFmodes =
            data.image[piaacmcopticaldesign.FmodesID].md[0].size[2];
        piaacmcopticaldesign.Fmsize =
            data.image[piaacmcopticaldesign.FmodesID].md[0].size[0];
    }

    printf("DONE Creating / loading Cmodes and Fmodes\n");
    fflush(stdout);

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * @brief Construct PIAA optics shape, or load if exists
 *
 * Purpose
 * -------
 * PIAA apodization is first specificed as a 2D amplitude image, named apo2Drad.fits on filesystem.
 *
 *
 * Arguments
 * ---------
 * @param[in]
 * piaacmctype  long
 *
 * @param[in]
 * size         uint32_t
 *
 * @param[in]
 * beamrad      float
 *              Optical beam radius in pixel unit
 *
 * @return errno_t
 */
static errno_t setupPIAAshapes(long piaacmctype, uint32_t size, float beamrad)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %ld %u %f", piaacmctype, size, beamrad);

    EXECUTE_SYSTEM_COMMAND("mkdir -p %s/piaaref/",
                           piaacmcparams.piaacmcconfdir);

    if (piaacmcopticaldesign.PIAAmode == 1)
    { // if coronagraph has PIAA optics

        // Do the files already exist in local memory ?
        piaacmcopticaldesign.piaa0CmodesID = image_ID("piaa0Cmodescoeff");
        piaacmcopticaldesign.piaa0FmodesID = image_ID("piaa0Fmodescoeff");
        piaacmcopticaldesign.piaa1CmodesID = image_ID("piaa1Cmodescoeff");
        piaacmcopticaldesign.piaa1FmodesID = image_ID("piaa1Fmodescoeff");

        // reference files will be placed in piaaref directory
        EXECUTE_SYSTEM_COMMAND("mkdir -p %s/piaaref/",
                               piaacmcparams.piaacmcconfdir);

        if ((piaacmcopticaldesign.piaa0CmodesID == -1) ||
            (piaacmcopticaldesign.piaa0FmodesID == -1) ||
            (piaacmcopticaldesign.piaa1CmodesID == -1) ||
            (piaacmcopticaldesign.piaa1FmodesID == -1))
        { // if any of the 4 files does not exist, all will be (re)computed

            // trying to load the 4 files from filesystem
            char fname[STRINGMAXLEN_FULLFILENAME];
            DEBUG_TRACEPOINT("loading mode coeffs");

            WRITE_FULLFILENAME(fname,
                               "%s/piaaref/piaa0Cmodes.fits",
                               piaacmcparams.piaacmcconfdir);
            FUNC_CHECK_RETURN(load_fits(fname,
                                        "piaa0Cmodescoeff",
                                        LOADFITS_ERRMODE_WARNING,
                                        &(piaacmcopticaldesign.piaa0CmodesID)));

            WRITE_FULLFILENAME(fname,
                               "%s/piaaref/piaa0Fmodes.fits",
                               piaacmcparams.piaacmcconfdir);
            FUNC_CHECK_RETURN(load_fits(fname,
                                        "piaa0Fmodescoeff",
                                        LOADFITS_ERRMODE_WARNING,
                                        &(piaacmcopticaldesign.piaa0FmodesID)));

            WRITE_FULLFILENAME(fname,
                               "%s/piaaref/piaa1Cmodes.fits",
                               piaacmcparams.piaacmcconfdir);
            FUNC_CHECK_RETURN(load_fits(fname,
                                        "piaa1Cmodescoeff",
                                        LOADFITS_ERRMODE_WARNING,
                                        &(piaacmcopticaldesign.piaa1CmodesID)));

            WRITE_FULLFILENAME(fname,
                               "%s/piaaref/piaa1Fmodes.fits",
                               piaacmcparams.piaacmcconfdir);
            FUNC_CHECK_RETURN(load_fits(fname,
                                        "piaa1Fmodescoeff",
                                        LOADFITS_ERRMODE_WARNING,
                                        &(piaacmcopticaldesign.piaa1FmodesID)));

            { // Loading APLC mask transmission file if it exists
                WRITE_FULLFILENAME(fname,
                                   "%s/piaaref/APLCmaskCtransm.txt",
                                   piaacmcparams.piaacmcconfdir);

                FILE *fp = fopen(fname, "r");
                if (fp != NULL)
                {
                    // if file not found, this is not a failure - it will be created below
                    DEBUG_TRACEPOINT("loading %s", fname);

                    float tmpf = 0.0;
                    if (fscanf(fp, "%f", &tmpf) == 1)
                    {
                        piaacmcopticaldesign.fpmaskamptransm = tmpf;
                    }
                    else
                    {
                        // but if file exists, value should be readable
                        FUNC_RETURN_FAILURE("Cannot read value from file %s",
                                            fname);
                    }
                    fclose(fp);
                }
                else
                {
                    DEBUG_TRACEPOINT("file %s not found", fname);
                }
            }
        }
    }

    // If above files could not be loaded, create them
    //
    if ((piaacmcopticaldesign.piaa0CmodesID == -1) ||
        (piaacmcopticaldesign.piaa0FmodesID == -1) ||
        (piaacmcopticaldesign.piaa1CmodesID == -1) ||
        (piaacmcopticaldesign.piaa1FmodesID == -1))
    {
        DEBUG_TRACEPOINT("creating piaa modes files");

        // We start from 2D apodization map
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname,
                           "%s/apo2Drad.fits",
                           piaacmcparams.piaacmcconfdir);

        imageID IDtmp1 = -1;
        FUNC_CHECK_RETURN(
            load_fits(fname, "apo2Drad", LOADFITS_ERRMODE_WARNING, &IDtmp1));

        // if the 2D apodization map cannot be loaded, then we create it
        if (IDtmp1 == -1) // CREATE APODIZATION
        {
            DEBUG_TRACEPOINT("creating apodization");

            EXECUTE_SYSTEM_COMMAND(
                "cp %s/piaaref/apo2Drad.fits %s/apo2Drad.fits",
                piaacmcparams.piaacmcconfdir,
                piaacmcparams.piaacmcconfdir);

            //sprintf(fname, "%s/apo2Drad.fits", piaacmcparams.piaacmcconfdir);
            imageID IDtmp2 = -1;
            FUNC_CHECK_RETURN(load_fits(fname,
                                        "apo2Drad",
                                        LOADFITS_ERRMODE_WARNING,
                                        &IDtmp2));

            if (IDtmp2 == -1)
            {
                DEBUG_TRACEPOINT(
                    "Creating 2D apodization for idealized circular "
                    "monochromatic PIAACMC");

                // 2D apodization is a prolate spheroidal function, iteratively computed
                // The function is first computed at low resolution (faster), and then expanded
                // to become the starting point of a higher resolution set of iterations

                // first iteration: half size image, 2x zoom, without pupil mask
                create_variable_ID("DFTZFACTOR", 2);
                create_variable_ID("PNBITER", 15);

                FUNC_CHECK_RETURN(coronagraph_make_2Dprolateld(
                    piaacmcopticaldesign.fpmaskradld,
                    beamrad * 0.5,
                    piaacmcopticaldesign.centObs1,
                    "apotmp1",
                    size / 2,
                    "NULLim"));

                // expand solution to full size
                basic_resizeim("apotmp1", "apostart", size, size);
                FUNC_CHECK_RETURN(
                    delete_image_ID("apotmp1", DELETE_IMAGE_ERRMODE_WARNING));

                // full size, 4x zoom
                create_variable_ID("DFTZFACTOR", 4);
                create_variable_ID("PNBITER", 5);
                FUNC_CHECK_RETURN(coronagraph_make_2Dprolateld(
                    piaacmcopticaldesign.fpmaskradld,
                    beamrad,
                    piaacmcopticaldesign.centObs1,
                    "apo",
                    size,
                    "pupmaskim"));

                // full size, 8x zoom
                chname_image_ID("apo", "apostart");
                create_variable_ID("DFTZFACTOR", 8);
                create_variable_ID("PNBITER", 5);
                FUNC_CHECK_RETURN(coronagraph_make_2Dprolateld(
                    piaacmcopticaldesign.fpmaskradld,
                    beamrad,
                    piaacmcopticaldesign.centObs1,
                    "apo",
                    size,
                    "pupmaskim"));

                // full size, 16x zoom
                chname_image_ID("apo", "apostart");
                create_variable_ID("DFTZFACTOR", 16);
                create_variable_ID("PNBITER", 10);
                FUNC_CHECK_RETURN(coronagraph_make_2Dprolateld(
                    piaacmcopticaldesign.fpmaskradld,
                    beamrad,
                    piaacmcopticaldesign.centObs1,
                    "apo",
                    size,
                    "pupmaskim"));

                chname_image_ID("apo", "apo2Drad");
                WRITE_FULLFILENAME(fname,
                                   "%s/apo2Drad.fits",
                                   piaacmcparams.piaacmcconfdir);
                FUNC_CHECK_RETURN(save_fits("apo2Drad", fname));

                if (piaacmcopticaldesign.PIAAmode == 1)
                {
                    WRITE_FULLFILENAME(fname,
                                       "%s/piaaref/apo2Drad.fits",
                                       piaacmcparams.piaacmcconfdir);
                    save_fits("apo2Drad", fname);
                }

                // When 2D optimization is complete, the focal plane mask transmission has converged.
                // Its value is saved to filesystem

                if ((piaacmctype == 0) &&
                    (loaded == 0)) // idealized focal plane mask
                {
                    DEBUG_TRACEPOINT("Idealized focal plane mask");

                    piaacmcopticaldesign.fpmaskamptransm =
                        -data.variable[variable_ID("APLCmaskCtransm")].value.f;
                    printf("FOCAL PLANE MASK TRANSM = %f\n",
                           piaacmcopticaldesign.fpmaskamptransm);
                    printf("Saving default configuration\n");
                    fflush(stdout);
                    saveconf = 1;

                    WRITE_FULLFILENAME(fname,
                                       "%s/piaaref/APLCmaskCtransm.txt",
                                       piaacmcparams.piaacmcconfdir);
                    {
                        FILE *fp = fopen(fname, "w");
                        fprintf(fp,
                                "%.20f\n",
                                piaacmcopticaldesign.fpmaskamptransm);
                        fclose(fp);
                    }
                }
            }
        }

        // split apodization in conventional pupil apodizer (apoCPA) and PIAA apodization (apo2Drad_PIAA)
        imageID  IDapo  = image_ID("apo2Drad");
        uint64_t xysize = 1;
        uint32_t xsize  = data.image[IDapo].md[0].size[0];
        xysize *= xsize;
        uint32_t ysize = data.image[IDapo].md[0].size[1];
        xysize *= ysize;

        imageID IDapo_PIAA;
        FUNC_CHECK_RETURN(
            create_2Dimage_ID("apo2Drad_PIAA", xsize, ysize, &IDapo_PIAA));

        imageID IDapo_CPA;
        FUNC_CHECK_RETURN(
            create_2Dimage_ID("apo2Drad_CPA", xsize, ysize, &IDapo_CPA));

        if (piaacmcopticaldesign.PIAAmode == 0)
        { // no PIAA optics

            DEBUG_TRACEPOINT("everything goes to the conventional apodizer");

            for (uint64_t ii = 0; ii < xysize; ii++)
            {
                data.image[IDapo_PIAA].array.F[ii] = 1.0;
                data.image[IDapo_CPA].array.F[ii] =
                    data.image[IDapo].array.F[ii];
            }
        }
        else
        { // PIAA optics

            DEBUG_TRACEPOINT("sharing podization, coeff = %f",
                             piaacmcopticaldesign.PIAAcoeff);
            for (uint64_t ii = 0; ii < xysize; ii++)
            {
                // fraction of apodization done by PIAA - between 0 and 1
                double coeff = piaacmcopticaldesign.PIAAcoeff;

                data.image[IDapo_PIAA].array.F[ii] =
                    pow(data.image[IDapo].array.F[ii], coeff);

                data.image[IDapo_CPA].array.F[ii] =
                    pow(data.image[IDapo].array.F[ii], 1.0 - coeff);
            }
        }

        copy_image_ID("apo2Drad_CPA", "prePIAA0mask", 0);

        FUNC_CHECK_RETURN(save_fits("prePIAA0mask", "prePIAA0mask.fits"));

        // To go from 2D apodization map to PIAA shapes, the apodization is first approximated
        // as a radial function, represented as a linear sum of cosines.
        // PIAA shapes are computed for this cosine apodization

        // load PIAA apodization profile and fit it a series of radial cosines
        // Output image "outApofit" is list of cosine coefficients
        FUNC_CHECK_RETURN(
            load2DRadialApodization("apo2Drad_PIAA", beamrad, "outApofit"));

        // compute radial PIAA mirror sag corresponding to the radial cosines
        // Input is radial cosine coefficients
        // output :
        // -> <piaacmcconfdir>/PIAA_Mshapes.txt
        FUNC_CHECK_RETURN(init_geomPIAA_rad("outApofit"));

        // make 2D sag maps from 1D radial sag profiles
        // output :
        // piaam0z, piaam1z
        //
        {
            DEBUG_TRACEPOINT("Make 2D sag maps");
            WRITE_FULLFILENAME(fname,
                               "%s/PIAA_Mshapes.txt",
                               piaacmcparams.piaacmcconfdir);

            double beamradpix =
                piaacmcopticaldesign.beamrad / piaacmcopticaldesign.pixscale;
            FUNC_CHECK_RETURN(
                mkPIAAMshapes_from_RadSag(fname,
                                          piaacmcopticaldesign.PIAAsep,
                                          piaacmcopticaldesign.beamrad,
                                          piaacmcopticaldesign.r0lim,
                                          piaacmcopticaldesign.r1lim,
                                          piaacmcopticaldesign.size,
                                          beamradpix,
                                          piaacmcopticaldesign.NBradpts,
                                          "piaam0z",
                                          "piaam1z"));
        }

        DEBUG_TRACEPOINT("piaacmcparams.PIAACMC_save %d",
                         piaacmcparams.PIAACMC_save);
        if (piaacmcparams.PIAACMC_save == 1)
        {
            DEBUG_TRACEPOINT("saving files");

            WRITE_FULLFILENAME(fname,
                               "%s/piaam0z.fits",
                               piaacmcparams.piaacmcconfdir);
            FUNC_CHECK_RETURN(save_fits("piaam0z", fname));

            WRITE_FULLFILENAME(fname,
                               "%s/piaam1z.fits",
                               piaacmcparams.piaacmcconfdir);
            FUNC_CHECK_RETURN(save_fits("piaam1z", fname));
        }

        // crop piaam0z and piaam1z to Cmodes size
        {
            imageID  ID0   = image_ID("Cmodes");
            uint32_t size0 = data.image[ID0].md[0].size[0];

            { // crop piaam0z -> piaa0zcrop
                imageID ID1;
                FUNC_CHECK_RETURN(
                    create_2Dimage_ID("piaa0zcrop", size0, size0, &ID1));

                imageID ID = image_ID("piaam0z");
                for (uint32_t ii = 0; ii < size0; ii++)
                    for (uint32_t jj = 0; jj < size0; jj++)
                    {
                        data.image[ID1].array.F[jj * size0 + ii] =
                            data.image[ID]
                                .array.F[(jj + (size - size0) / 2) * size +
                                         (ii + (size - size0) / 2)];
                    }
            }

            { // crop piaam1z -> piaa1zcrop

                imageID ID1;
                FUNC_CHECK_RETURN(
                    create_2Dimage_ID("piaa1zcrop", size0, size0, &ID1));

                imageID ID = image_ID("piaam1z");
                for (uint32_t ii = 0; ii < size0; ii++)
                    for (uint32_t jj = 0; jj < size0; jj++)
                    {
                        data.image[ID1].array.F[jj * size0 + ii] =
                            data.image[ID]
                                .array.F[(jj + (size - size0) / 2) * size +
                                         (ii + (size - size0) / 2)];
                    }
            }

            make_disk("maskd", size0, size0, 0.5 * size0, 0.5 * size0, beamrad);
            make_2Dgridpix("gridpix", size0, size0, 1, 1, 0, 0);
        }
        arith_image_mult("maskd", "gridpix", "maskfit");

        //sprintf(fname, "%s/maskfit.fits", piaacmcparams.piaacmcconfdir);
        //save_fits("maskfit", fname);

        // The 2D sag maps (piaam0z and piaam1z) are now fitted as a linear sum of cosines
        // Results (coefficient values) are stored in piaa0Cmodescoeff and piaa1Cmodescoeff
        //
        printf("--------- FITTING COSINE MODES ---------\n");
        fflush(stdout);

        FUNC_CHECK_RETURN(linopt_imtools_image_fitModes("piaa0zcrop",
                                                        "Cmodes",
                                                        "maskfit",
                                                        1.0e-6,
                                                        "piaa0Cmodescoeff",
                                                        0,
                                                        NULL));

        EXECUTE_SYSTEM_COMMAND("mv %s/eigenv.dat %s/eigenv_piaa0Cmodes.dat",
                               piaacmcparams.piaacmcconfdir,
                               piaacmcparams.piaacmcconfdir);

        //      sprintf(fname, "%s/piaa0Cmodescoeff.fits", piaacmcparams.piaacmcconfdir);
        //      save_fits("piaa0Cmodescoeff", fname);

        FUNC_CHECK_RETURN(linopt_imtools_image_fitModes("piaa1zcrop",
                                                        "Cmodes",
                                                        "maskfit",
                                                        1.0e-6,
                                                        "piaa1Cmodescoeff",
                                                        0,
                                                        NULL));

        EXECUTE_SYSTEM_COMMAND("mv eigenv.dat %s/eigenv_piaa1Cmodes.dat",
                               piaacmcparams.piaacmcconfdir);

        // The coefficients are re-expanded to 2D maps :
        // -> piaa0Cz
        // -> piaa1Cz
        //
        FUNC_CHECK_RETURN(linopt_imtools_image_construct("Cmodes",
                                                         "piaa0Cmodescoeff",
                                                         "piaa0Cz",
                                                         NULL));

        if (piaacmcparams.PIAACMC_save == 1)
        {
            WRITE_FULLFILENAME(fname,
                               "%s/initpiaacmc_piaa0Cz.fits",
                               piaacmcparams.piaacmcconfdir);
            DEBUG_TRACEPOINT("saving %s -> %s", "piaa0Cz", fname);
            FUNC_CHECK_RETURN(save_fits("piaa0Cz", fname));

            WRITE_FULLFILENAME(fname,
                               "%s/initpiaacmc_piaa0Cmodescoeff.fits",
                               piaacmcparams.piaacmcconfdir);
            FUNC_CHECK_RETURN(save_fits("piaa0Cmodescoeff", fname));

            WRITE_FULLFILENAME(fname,
                               "%s/initpiaacmc_Cmodes.fits",
                               piaacmcparams.piaacmcconfdir);
            FUNC_CHECK_RETURN(save_fits("Cmodes", fname));
        }

        FUNC_CHECK_RETURN(linopt_imtools_image_construct("Cmodes",
                                                         "piaa1Cmodescoeff",
                                                         "piaa1Cz",
                                                         NULL));

        if (piaacmcparams.PIAACMC_save == 1)
        {
            WRITE_FULLFILENAME(fname,
                               "%s/initpiaacmc_piaa1Cz.fits",
                               piaacmcparams.piaacmcconfdir);
            DEBUG_TRACEPOINT("saving %s -> %s", "piaa1Cz", fname);
            FUNC_CHECK_RETURN(save_fits("piaa1Cz", fname));
        }

        // Compute residual difference between original 2D shape and fit
        //
        { // piaam0z -> piaa0Cres
            imageID  ID0   = image_ID("piaa0Cz");
            uint32_t size0 = data.image[ID0].md[0].size[0];
            imageID  ID1   = image_ID("piaam0z");

            imageID ID;
            FUNC_CHECK_RETURN(
                create_2Dimage_ID("piaa0Cres", size0, size0, &ID));

            for (uint32_t ii = 0; ii < size0; ii++)
                for (uint32_t jj = 0; jj < size0; jj++)
                {
                    data.image[ID].array.F[jj * size0 + ii] =
                        data.image[ID1]
                            .array.F[(jj + (size - size0) / 2) * size +
                                     (ii + (size - size0) / 2)] -
                        data.image[ID0].array.F[jj * size0 + ii];
                }
        }
        if (piaacmcparams.PIAACMC_save == 1)
        {
            WRITE_FULLFILENAME(fname,
                               "%s/initpiaacmc_piaa0Cres.fits",
                               piaacmcparams.piaacmcconfdir);
            DEBUG_TRACEPOINT("saving %s -> %s", "piaa0Cres", fname);
            FUNC_CHECK_RETURN(save_fits("piaa0Cres", fname));
        }

        { // piaam1z -> piaa1Cres
            imageID  ID0   = image_ID("piaa1Cz");
            uint32_t size0 = data.image[ID0].md[0].size[0];
            imageID  ID1   = image_ID("piaam1z");

            imageID ID;
            FUNC_CHECK_RETURN(
                create_2Dimage_ID("piaa1Cres", size0, size0, &ID));

            for (uint32_t ii = 0; ii < size0; ii++)
                for (uint32_t jj = 0; jj < size0; jj++)
                {
                    data.image[ID].array.F[jj * size0 + ii] =
                        data.image[ID1]
                            .array.F[(jj + (size - size0) / 2) * size +
                                     (ii + (size - size0) / 2)] -
                        data.image[ID0].array.F[jj * size0 + ii];
                }
        }
        if (piaacmcparams.PIAACMC_save == 1)
        {
            WRITE_FULLFILENAME(fname,
                               "%s/initpiaacmc_piaa1Cres.fits",
                               piaacmcparams.piaacmcconfdir);
            DEBUG_TRACEPOINT("saving %s -> %s", "piaa1Cres", fname);
            FUNC_CHECK_RETURN(save_fits("piaa1Cres", fname));
        }

        // The residuals are fitted as Fourier modes
        // This should be small
        //
        printf("--------- FITTING FOURIER MODES ---------\n");
        fflush(stdout);

        FUNC_CHECK_RETURN(linopt_imtools_image_fitModes("piaa0Cres",
                                                        "Fmodes",
                                                        "maskfit",
                                                        0.01,
                                                        "piaa0Fmodescoeff",
                                                        0,
                                                        NULL));

        EXECUTE_SYSTEM_COMMAND("mv %s/eigenv.dat %s/eigenv_piaa0Fmodes.dat",
                               piaacmcparams.piaacmcconfdir,
                               piaacmcparams.piaacmcconfdir);

        FUNC_CHECK_RETURN(linopt_imtools_image_fitModes("piaa1Cres",
                                                        "Fmodes",
                                                        "maskfit",
                                                        0.01,
                                                        "piaa1Fmodescoeff",
                                                        0,
                                                        NULL));

        EXECUTE_SYSTEM_COMMAND("mv %s/eigenv.dat %s/eigenv_piaa1Fmodes.dat",
                               piaacmcparams.piaacmcconfdir,
                               piaacmcparams.piaacmcconfdir);

        // save_fits("piaa1Fmodescoeff", "piaa1Fmodescoeff.fits");

        //linopt_imtools_image_construct("Fmodes", "piaa0Fmodescoeff", "piaa0Fz");
        //   save_fits("piaa0Fz", "piaa0Fz.fits");
        //arith_image_sub("piaa0Cres", "piaa0Fz", "piaa0CFres");
        //save_fits("piaa0CFres", "piaa0CFres.fits");
        FUNC_CHECK_RETURN(
            delete_image_ID("piaa0zcrop", DELETE_IMAGE_ERRMODE_WARNING));

        //linopt_imtools_image_construct("Fmodes", "piaa1Fmodescoeff", "piaa1Fz");
        //save_fits("piaa1Fz", "piaa1Fz.fits");
        //arith_image_sub("piaa1Cres", "piaa1Fz", "piaa1CFres");
        //save_fits("piaa1CFres", "piaa1CFres.fits");
        FUNC_CHECK_RETURN(
            delete_image_ID("piaa1zcrop", DELETE_IMAGE_ERRMODE_WARNING));

        FUNC_CHECK_RETURN(
            delete_image_ID("maskfit", DELETE_IMAGE_ERRMODE_WARNING));

        piaacmcopticaldesign.piaa0CmodesID = image_ID("piaa0Cmodescoeff");
        piaacmcopticaldesign.piaa0FmodesID = image_ID("piaa0Fmodescoeff");
        piaacmcopticaldesign.piaa1CmodesID = image_ID("piaa1Cmodescoeff");
        piaacmcopticaldesign.piaa1FmodesID = image_ID("piaa1Fmodescoeff");

        // At this point, we have a modal representation of PIAA optics shapes
        //
        DEBUG_TRACEPOINT("saving piaa modes files");

        WRITE_FULLFILENAME(fname,
                           "%s/piaaref/piaa0Cmodes.fits",
                           piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(
            save_fits(data.image[piaacmcopticaldesign.piaa0CmodesID].name,
                      fname));

        WRITE_FULLFILENAME(fname,
                           "%s/piaaref/piaa0Fmodes.fits",
                           piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(
            save_fits(data.image[piaacmcopticaldesign.piaa0FmodesID].name,
                      fname));

        WRITE_FULLFILENAME(fname,
                           "%s/piaaref/piaa1Cmodes.fits",
                           piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(
            save_fits(data.image[piaacmcopticaldesign.piaa1CmodesID].name,
                      fname));

        WRITE_FULLFILENAME(fname,
                           "%s/piaaref/piaa1Fmodes.fits",
                           piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(
            save_fits(data.image[piaacmcopticaldesign.piaa1FmodesID].name,
                      fname));

        EXECUTE_SYSTEM_COMMAND("cp %s/piaaref/* %s/",
                               piaacmcparams.piaacmcconfdir,
                               piaacmcparams.piaacmcconfdir);
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * @brief Creates/initializes piaacmcconf structure and directory
 *
 * @param[in] fpmradld     Focal plane mask nominal radius
 * @param[in] centobs0     Input central obstruction
 * @param[in] centobs1     Output central obstruction
 * @param[in] flags        flags, see .h file
 *
 * piaacmctype:
 * - 0: if configuration does not exist, create Monochromatic idealized PIAACMC, otherwise, read configuration
 * - 1: physical mask
 *
 *
 */
errno_t init_piaacmcopticaldesign(double    fpmradld,
                                  double    centobs0,
                                  double    centobs1,
                                  uint64_t  flags,
                                  uint64_t *status)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %f %f %f %lu", fpmradld, centobs0, centobs1, flags);

    FILE *fp;
    float beamradpix;
    long  size;

    uint32_t xsize  = 0;
    uint32_t ysize  = 0;
    uint64_t xysize = 1;

    // Create required directories
    int ret;

    ret = mkdir("testdir", 0777); // test files, output from code
    (void) ret;

    ret = mkdir("conf", 0777); // configuration
    (void) ret;

    ret = mkdir("status", 0777); // status
    (void) ret;

    ret = mkdir("ref", 0777); // reference files
    (void) ret;

    ret = mkdir("log", 0777); // log
    (void) ret;

    (void) status;

    DEBUG_TRACEPOINT("piaacmcopticaldesign.fpmaskamptransm %lf",
                     piaacmcopticaldesign.fpmaskamptransm);

    if (flags & INIT_PIAACMCOPTICALDESIGN_MODE__READCONF)
    {
        FUNC_CHECK_RETURN(PIAACMCsimul_initpiaacmcconf_readconfparams(
            (int) (flags & INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL),
            fpmradld,
            centobs0,
            centobs1,
            (int) (flags & INIT_PIAACMCOPTICALDESIGN_MODE__WSCMODE)));
    }

    DEBUG_TRACEPOINT("piaacmcopticaldesign.fpmaskamptransm %lf",
                     piaacmcopticaldesign.fpmaskamptransm);

    piaacmcopticaldesign.fpmCentConeRad = piaacmcopticaldesign.fpmRad *
                                          piaacmcopticaldesign.NBringCentCone /
                                          piaacmcopticaldesign.NBrings;

    printf("fpmRad = %g m\n",
           0.5 * (piaacmcparams.LAMBDASTART + piaacmcparams.LAMBDAEND) *
               piaacmcopticaldesign.Fratio * fpmradld);

    printf("factor = %f   (%ld / %ld)\n",
           1.0 * piaacmcopticaldesign.NBringCentCone /
               piaacmcopticaldesign.NBrings,
           piaacmcopticaldesign.NBringCentCone,
           piaacmcopticaldesign.NBrings);

    printf("fpmCentConeRad =  %g\n",
           (0.5 * (piaacmcparams.LAMBDASTART + piaacmcparams.LAMBDAEND) *
            piaacmcopticaldesign.Fratio * fpmradld) *
               piaacmcopticaldesign.NBringCentCone /
               piaacmcopticaldesign.NBrings);

    if (flags & INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF)
    {
        DEBUG_TRACEPOINT("Loading PIAACMC configuration");
        EXECUTE_SYSTEM_COMMAND("mkdir -p %s", piaacmcparams.piaacmcconfdir);
        loaded = PIAACMCsimul_loadpiaacmcconf(piaacmcparams.piaacmcconfdir);
        if (loaded == 0)
        {
            DEBUG_TRACEPOINT("Saving default configuration");
            saveconf = 1;
        }
    }

    DEBUG_TRACEPOINT("piaacmcopticaldesign.fpmaskamptransm %lf",
                     piaacmcopticaldesign.fpmaskamptransm);

    { // read PIAA material
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname,
                           "%s/conf_PIAAmaterial_name.txt",
                           piaacmcparams.piaacmcconfdir);
        if ((fp = fopen(fname, "r")) != NULL)
        {
            char name[200];
            if (fscanf(fp, "%199s", name) == 1)
            {
                strcpy(piaacmcopticaldesign.PIAAmaterial_name, name);
            }
            else
            {
                FUNC_RETURN_FAILURE("Cannot read value from file %s", fname);
            }
            fclose(fp);
        }
        else
        {
            snprintf(piaacmcopticaldesign.PIAAmaterial_name,
                     STRINGMAXLEN_PIAACMCSIMUL_MATERIALNAME,
                     "Mirror");
            WRITE_FULLFILENAME(fname,
                               "%s/conf_PIAAmaterial_name.txt",
                               piaacmcparams.piaacmcconfdir);
            if ((fp = fopen(fname, "w")) != NULL)
            {
                fprintf(fp, "%s\n", piaacmcopticaldesign.PIAAmaterial_name);
                fclose(fp);
            }
            else
            {
                FUNC_RETURN_FAILURE("Cannot create file \"%s\"", fname);
            }
        }
    }

    DEBUG_TRACEPOINT("piaacmcopticaldesign.PIAAmaterial_name = %s",
                     piaacmcopticaldesign.PIAAmaterial_name);
    piaacmcopticaldesign.PIAAmaterial_code =
        OpticsMaterials_code(piaacmcopticaldesign.PIAAmaterial_name);
    DEBUG_TRACEPOINT("piaacmcopticaldesign.PIAAmaterial_code = %d",
                     piaacmcopticaldesign.PIAAmaterial_code);

    { // write PIAA material
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname,
                           "%s/conf_PIAAmaterial_code.txt",
                           piaacmcparams.piaacmcconfdir);
        if ((fp = fopen(fname, "w")) != NULL)
        {
            fprintf(fp, "%d\n", piaacmcopticaldesign.PIAAmaterial_code);
            fclose(fp);
        }
        else
        {
            FUNC_RETURN_FAILURE("cannot create file \"%s\"", fname);
        }
    }

    printf("lambda = %g\n", piaacmcopticaldesign.lambda);
    printf("LAMBDASTART = %g\n", piaacmcparams.LAMBDASTART);
    printf("LAMBDAEND = %g\n", piaacmcparams.LAMBDAEND);

    for (int k = 0; k < piaacmcopticaldesign.nblambda; k++)
    {
        piaacmcopticaldesign.lambdaarray[k] =
            piaacmcparams.LAMBDASTART +
            (0.5 + k) * (piaacmcparams.LAMBDAEND - piaacmcparams.LAMBDASTART) /
                piaacmcopticaldesign.nblambda;
    }

    // create modes for aspheric optical surfaces description
    beamradpix = piaacmcopticaldesign.beamrad / piaacmcopticaldesign.pixscale;
    size       = piaacmcopticaldesign.size;

    DEBUG_TRACEPOINT("BEAM RADIUS :  %f / %f =  %f pix, size = %ld",
                     piaacmcopticaldesign.beamrad,
                     piaacmcopticaldesign.pixscale,
                     beamradpix,
                     size);

    // x, y, r and PA coordinates in beam (for convenience & speed)
    make_x_y_r_PA_images(size, beamradpix);

    // ==================== CREATE DMs ===============
    DEBUG_TRACEPOINT("%d DM(s)", piaacmcopticaldesign.nbDM);
    for (long iDM = 0; iDM < piaacmcopticaldesign.nbDM; iDM++)
    {
        printf("DM # %ld\n", iDM);

        char imname[STRINGMAXLEN_IMGNAME];
        WRITE_IMAGENAME(imname, "wfcDM%ld", iDM);

        // DM image identifier - to be updated later
        piaacmcopticaldesign.ID_DM[iDM] = image_ID(imname);
        printf("ID = %ld", piaacmcopticaldesign.ID_DM[iDM]);

        if (piaacmcopticaldesign.ID_DM[iDM] == -1)
        {
            read_sharedmem_image(imname);
            piaacmcopticaldesign.ID_DM[iDM] = image_ID(imname);
        }

        if (piaacmcopticaldesign.ID_DM[iDM] == -1)
        {
            uint32_t *sizearray;
            sizearray = (uint32_t *) malloc(sizeof(uint32_t) * 2);
            if (sizearray == NULL)
            {
                FUNC_RETURN_FAILURE("malloc returns NULL pointer");
            }
            sizearray[0] = size;
            sizearray[1] = size;
            FUNC_CHECK_RETURN(
                create_image_ID(imname,
                                2,
                                sizearray,
                                _DATATYPE_FLOAT,
                                1,
                                0,
                                0,
                                &(piaacmcopticaldesign.ID_DM[iDM])));
            free(sizearray);
        }
    }

    // ==================== CREATE MODES USED TO FIT AND DESCRIBE PIAA SHAPES ===============
    make_C_F_modes(size, beamradpix);

    // =================== IMPORT / CREATE PIAA SHAPES =====================
    FUNC_CHECK_RETURN(setupPIAAshapes(
        (int) (flags & INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL),
        size,
        beamradpix));

    // ============ MAKE FOCAL PLANE MASK ===============

    PIAACMCsimul_update_fnamedescr();

    EXECUTE_SYSTEM_COMMAND("echo \"%s\" > fpm_name.txt",
                           piaacmcparams.fnamedescr);

    EXECUTE_SYSTEM_COMMAND("echo \"%s\" > fpm_name_conf.txt",
                           piaacmcparams.fnamedescr_conf);

    /*    piaacmcparams.CREATE_fpmzmap = 0;
        if( piaacmcparams.FORCE_CREATE_fpmzmap == 0 )
        {
            if(image_ID("fpmzmap")==-1)
                {
                    sprintf(fname, "%s/fpmzmap%d_%03ld_%03ld.fits", piaacmcparams.piaacmcconfdir, piaacmcparams.PIAACMC_FPMsectors, piaacmcopticaldesign.NBrings, piaacmcopticaldesign.focmNBzone);
                    load_fits(fname, "fpmzmap", 1);
                    if(image_ID("fpmzmap")==-1)
                        piaacmcparams.CREATE_fpmzmap = 1;
                }
        }
        else
            piaacmcparams.CREATE_fpmzmap = 1;

        if( piaacmcparams.CREATE_fpmzmap == 1 )
        {
            if(image_ID("fpmzmap")!=-1)
                delete_image_ID("fpmzmap", DELETE_IMAGE_ERRMODE_WARNING);
            PIAACMCsimul_mkFPM_zonemap("fpmzmap");
            sprintf(fname, "%s/fpmzmap%d_%03ld_%03ld.fits", piaacmcparams.piaacmcconfdir, piaacmcparams.PIAACMC_FPMsectors, piaacmcopticaldesign.NBrings, piaacmcopticaldesign.focmNBzone);
            save_fits("fpmzmap", fname);
        }
    */

    if (image_ID("fpmzmap") == -1)
    {
        printf("Make zonemap ...\n");
        fflush(stdout);
        FUNC_CHECK_RETURN(mkFPM_zonemap("fpmzmap", NULL));
    }
    else
    {
        printf("zonemap already exists\n");
        fflush(stdout);
    }

    //    sprintf(fname, "%s/fpmzmap.fits", piaacmcparams.piaacmcconfdir);
    //    save_fits("fpmzmap", fname);

    /* sprintf(fname, "%s/fpmzmap.fits", piaacmcparams.piaacmcconfdir);
     save_fits("fpmzmap", fname);
     exit(0);*/

    // zones thickness

    piaacmcparams.CREATE_fpmzt = 0;
    if (piaacmcparams.FORCE_CREATE_fpmzt == 0)
    {
        piaacmcopticaldesign.zonezID = image_ID("fpmzt");
        if (piaacmcopticaldesign.zonezID == -1)
        {
            PIAACMCsimul_update_fnamedescr();

            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(fname,
                               "%s/fpm_zonez.%s.fits",
                               piaacmcparams.piaacmcconfdir,
                               piaacmcparams.fnamedescr);

            printf("LOADING FILE NAME : \"%s\"  -  %d %d \n",
                   fname,
                   (int) (flags & INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL),
                   loaded);

            load_fits(fname,
                      "fpmzt",
                      LOADFITS_ERRMODE_WARNING,
                      &(piaacmcopticaldesign.zonezID));
            if (piaacmcopticaldesign.zonezID == -1)
            {
                piaacmcparams.CREATE_fpmzt = 1;
            }
        }
    }
    else
    {
        piaacmcparams.CREATE_fpmzt = 1;
    }
    DEBUG_TRACEPOINT("piaacmcparams.CREATE_fpmzt = %d",
                     piaacmcparams.CREATE_fpmzt);

    if (piaacmcparams.CREATE_fpmzt == 1)
    {
        printf("Creating fpmzt, saving as fpm_zonez.fits - %d %d\n",
               (int) (flags & INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL),
               loaded);
        fflush(stdout);
        piaacmcopticaldesign.zonezID = image_ID("fpmzt");
        if (piaacmcopticaldesign.zonezID != -1)
        {
            delete_image_ID("fpmzt", DELETE_IMAGE_ERRMODE_WARNING);
        }

        create_2Dimage_ID_double("fpmzt",
                                 piaacmcopticaldesign.focmNBzone,
                                 1,
                                 &(piaacmcopticaldesign.zonezID));
        double t = 1.0e-9;

        if (flags & INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL)
        { // physical focal plane mask
            DEBUG_TRACEPOINT(
                "CREATING EXAMPLE FOCAL PLANE MASK  %d %d",
                (int) (flags & INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL),
                loaded);
        }
        else
        { // idealized focal plane mask
            DEBUG_TRACEPOINT("IDEALIZED FOCAL PLANE MASK");

            // measure dpha/dt
            double t0   = 1.0e-8;
            double pha0 = OpticsMaterials_pha_lambda(
                piaacmcopticaldesign.fpmmaterial_code,
                t0,
                0.5 * (piaacmcparams.LAMBDASTART + piaacmcparams.LAMBDAEND));
            // set t to get PI phase
            t = (M_PI / pha0) * t0;
            printf("t = %g m (%lf %g) -> %g %g\n",
                   t,
                   pha0,
                   t0,
                   OpticsMaterials_pha_lambda(
                       piaacmcopticaldesign.fpmmaterial_code,
                       t,
                       0.5 * (piaacmcparams.LAMBDASTART +
                              piaacmcparams.LAMBDAEND)),
                   0.5 * (piaacmcparams.LAMBDASTART + piaacmcparams.LAMBDAEND));
            fflush(stdout);

            printf(" -- lambda = %g\n", piaacmcopticaldesign.lambda);
            printf(" -- lambdaB = %g\n", piaacmcopticaldesign.lambdaB);
            printf(" -- LAMBDASTART = %g\n", piaacmcparams.LAMBDASTART);
            printf(" -- LAMBDAEND = %g\n", piaacmcparams.LAMBDAEND);
        }

        for (long ii = 0; ii < piaacmcopticaldesign.focmNBzone; ii++)
        {
            data.image[piaacmcopticaldesign.zonezID].array.D[ii] = t;
        }

        {
            PIAACMCsimul_update_fnamedescr();
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(fname,
                               "%s/fpm_zonez.%s.fits",
                               piaacmcparams.piaacmcconfdir,
                               piaacmcparams.fnamedescr);

            printf("Writing %s\n", fname);
            save_fits("fpmzt", fname);
        }
    }

    // zones transmission amplitude

    printf("CREATE_fpmza = %d\n", piaacmcparams.CREATE_fpmza);
    fflush(stdout);
    if (piaacmcparams.FORCE_CREATE_fpmza == 0)
    {
        piaacmcopticaldesign.zoneaID = image_ID("fpmza");
        if (piaacmcopticaldesign.zoneaID == -1)
        {
            PIAACMCsimul_update_fnamedescr();
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(fname,
                               "%s/fpm_zonea.%s.fits",
                               piaacmcparams.piaacmcconfdir,
                               piaacmcparams.fnamedescr);

            printf("LOADING FILE NAME : \"%s\"\n", fname);
            load_fits(fname,
                      "fpmza",
                      LOADFITS_ERRMODE_IGNORE,
                      &(piaacmcopticaldesign.zoneaID));

            if (piaacmcopticaldesign.zoneaID == -1)
            {
                piaacmcparams.CREATE_fpmza = 1;
            }
        }
    }
    else
    {
        piaacmcparams.CREATE_fpmza = 1;
    }

    if (piaacmcparams.CREATE_fpmza == 1)
    {
        if (piaacmcopticaldesign.zoneaID != -1)
        {
            delete_image_ID("fpmza", DELETE_IMAGE_ERRMODE_WARNING);
        }
        create_2Dimage_ID_double("fpmza",
                                 piaacmcopticaldesign.focmNBzone,
                                 1,
                                 &(piaacmcopticaldesign.zoneaID));

        if (piaacmcparams.PIAACMC_MASKRADLD > 0.2) // physical mask
        {
            printf("PHYSICAL MASK ... %ld zones\n",
                   piaacmcopticaldesign.focmNBzone);
            fflush(stdout);

            for (long ii = 0; ii < piaacmcopticaldesign.focmNBzone; ii++)
            {
                data.image[piaacmcopticaldesign.zoneaID].array.D[ii] = 1.0;
            }
        }
        else // idealized mask
        {
            printf("IDEALIZED MASK ... %ld zones\n",
                   piaacmcopticaldesign.focmNBzone);
            for (long ii = 0; ii < piaacmcopticaldesign.focmNBzone; ii++)
            {
                data.image[piaacmcopticaldesign.zoneaID].array.D[ii] =
                    piaacmcopticaldesign.fpmaskamptransm;
            }
        }

        {
            PIAACMCsimul_update_fnamedescr();
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(fname,
                               "%s/fpm_zonea.%s.fits",
                               piaacmcparams.piaacmcconfdir,
                               piaacmcparams.fnamedescr);

            printf("Writing %s\n", fname);
            FUNC_CHECK_RETURN(save_fits("fpmza", fname));
        }
    }

    //   printf("%d piaacmcopticaldesign.fpmaskamptransm = %f       %lf\n", piaacmcparams.CREATE_fpmza, piaacmcopticaldesign.fpmaskamptransm, data.image[piaacmcopticaldesign.zoneaID].array.D[0]);
    //   sleep(10);

    // ============= MAKE LYOT STOPS =======================
    printf(
        "LOADING/CREATING LYOT MASK  - %ld masks  (PIAAmode = %d, %ld x %ld)\n",
        piaacmcopticaldesign.NBLyotStop,
        piaacmcopticaldesign.PIAAmode,
        (long) xsize,
        (long) ysize);
    // list_image_ID();
    //size2 = size * size;

    if (piaacmcopticaldesign.PIAAmode == 1)
    {
        for (long i = 0; i < piaacmcopticaldesign.NBLyotStop; i++)
        {
            printf("LYOT MASK %ld\n", i);
            fflush(stdout);

            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(fname,
                               "%s/LyotStop%ld.fits",
                               piaacmcparams.piaacmcconfdir,
                               i);

            char name[STRINGMAXLEN_IMGNAME];
            WRITE_IMAGENAME(name, "lyotstop%ld", i);

            piaacmcopticaldesign.IDLyotStop[i] = image_ID(name);
            if (piaacmcopticaldesign.IDLyotStop[i] == -1)
            {
                WRITE_FULLFILENAME(fname,
                                   "%s/LyotStop%ld.fits",
                                   piaacmcparams.piaacmcconfdir,
                                   i);
                if ((i == 0) && (piaacmcopticaldesign.NBLyotStop > 1))
                {
                    FUNC_CHECK_RETURN(mkSimpleLyotStop(
                        name,
                        -0.01,
                        0.98,
                        &(piaacmcopticaldesign.IDLyotStop[i])));
                }
                else if ((i == 1) && (piaacmcopticaldesign.NBLyotStop > 2))
                {
                    FUNC_CHECK_RETURN(mkSimpleLyotStop(
                        name,
                        piaacmcopticaldesign.centObs1 + 0.02,
                        1.2,
                        &(piaacmcopticaldesign.IDLyotStop[i])));
                }
                else
                {
                    FUNC_CHECK_RETURN(mkSimpleLyotStop(
                        name,
                        piaacmcopticaldesign.centObs1 + 0.02,
                        0.98,
                        &(piaacmcopticaldesign.IDLyotStop[i])));
                }

                save_fl_fits(name, fname);
            }
        }
    }
    else
    {
        long i = 0;

        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname,
                           "%s/LyotStop%ld.fits",
                           piaacmcparams.piaacmcconfdir,
                           i);

        char name[STRINGMAXLEN_IMGNAME];
        WRITE_IMAGENAME(name, "lyotstop%ld", i);

        piaacmcopticaldesign.IDLyotStop[i] = image_ID(name);
        if (piaacmcopticaldesign.IDLyotStop[i] == -1)
        {
            FUNC_CHECK_RETURN(
                create_2Dimage_ID(name,
                                  xsize,
                                  ysize,
                                  &(piaacmcopticaldesign.IDLyotStop[i])));
            imageID ID = image_ID("pupmaskim");
            for (uint64_t ii = 0; ii < xysize; ii++)
                if (data.image[ID].array.F[ii] < 0.99999999)
                {
                    data.image[piaacmcopticaldesign.IDLyotStop[i]].array.F[ii] =
                        0.0;
                }
                else
                {
                    data.image[piaacmcopticaldesign.IDLyotStop[i]].array.F[ii] =
                        1.0;
                }

            for (uint32_t ii = 0; ii < xsize; ii++)
                for (uint32_t jj = 0; jj < ysize; jj++)
                {
                    double x   = 1.0 * ii - 0.5 * xsize;
                    double y   = 1.0 * jj - 0.5 * ysize;
                    double rad = sqrt(x * x + y * y);
                    rad /= beamradpix;
                    if (rad <
                        (piaacmcopticaldesign.centObs1 + 0.5 / beamradpix))
                    {
                        data.image[piaacmcopticaldesign.IDLyotStop[i]]
                            .array.F[jj * xsize + ii] = 0.0;
                    }
                    if (rad > (1.0 - 0.5 / beamradpix))
                    {
                        data.image[piaacmcopticaldesign.IDLyotStop[i]]
                            .array.F[jj * xsize + ii] = 0.0;
                    }
                }
            FUNC_CHECK_RETURN(save_fl_fits(name, fname));
        }
    }

    if (saveconf == 1)
    {
        FUNC_CHECK_RETURN(
            PIAACMCsimul_savepiaacmcconf(piaacmcparams.piaacmcconfdir));
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
