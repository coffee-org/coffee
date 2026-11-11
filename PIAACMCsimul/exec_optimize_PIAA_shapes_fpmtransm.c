/**
 * @file    exec_optimize_PIAA_shapes_fpmtransm.c
 * @brief   Optimize PIAA shapes and focal plane mask transmission (Mode 40)
 */

#include <stdint.h>
#include <stdio.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"
#include "exec_optimize_PIAA_shapes_fpmtransm.h"
#include "init_piaacmcopticaldesign.h"

#include "PIAAshape/makePIAAshapes.h"

/**
 * opt_piaa_fpm_init_design - Read CLI overrides and initialize optical design and PIAA surfaces.
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t opt_piaa_fpm_init_design(void)
{
    variableID IDv;
    double     fpmradld = 0.95;
    double     centobs0 = 0.3;
    double     centobs1 = 0.2;

    if ((IDv = variable_ID("PIAACMC_centobs0")) != -1)
    {
        centobs0 = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_centobs1")) != -1)
    {
        centobs1 = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = coffee_dcvar_f(IDv);
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }

    piaacmcparams.PIAACMC_fpmtype = 0;
    if ((IDv = variable_ID("PIAACMC_fpmtype")) != -1)
    {
        piaacmcparams.PIAACMC_fpmtype = (int) (coffee_dcvar_f(IDv) + 0.1);
    }
    printf("PIAACMC_fpmtype = %d\n", piaacmcparams.PIAACMC_fpmtype);

    piaacmcparams.FORCE_CREATE_fpmza = 1;
    {
        uint64_t initflag = INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__READCONF;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF;
        if (piaacmcparams.PIAACMC_fpmtype == 1)
        {
            initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL;
        }

        FUNC_CHECK_RETURN(
            init_piaacmcopticaldesign(fpmradld, centobs0, centobs1, initflag, NULL));
    }

    FUNC_CHECK_RETURN(makePIAAshapes());
    piaacmcopticalsystem.FOCMASKarray[0].mode = 1;

    return RETURN_SUCCESS;
}

/**
 * opt_piaa_fpm_configure_regularization - Configure PIAA modal and FPM sag regularization weights.
 */
static void opt_piaa_fpm_configure_regularization(void)
{
    variableID IDv;

    piaacmcparams.linopt_REGPIAASHAPES = 0;
    if ((IDv = variable_ID("REGPIAASHAPES")) != -1)
    {
        piaacmcparams.linopt_REGPIAASHAPES = (long) coffee_dcvar_f(IDv) + 0.01;
    }

    piaacmcparams.linopt_piaa0C_regcoeff = 0.0e-7f;
    piaacmcparams.linopt_piaa1C_regcoeff = 0.0e-7f;
    if ((IDv = variable_ID("REGPIAA_C_COEFF")) != -1)
    {
        piaacmcparams.linopt_piaa0C_regcoeff = coffee_dcvar_f(IDv);
        piaacmcparams.linopt_piaa1C_regcoeff = coffee_dcvar_f(IDv);
    }

    piaacmcparams.linopt_piaa0C_regcoeff_alpha = 1.0f;
    piaacmcparams.linopt_piaa1C_regcoeff_alpha = 1.0f;
    if ((IDv = variable_ID("REGPIAA_C_ALPHA")) != -1)
    {
        piaacmcparams.linopt_piaa0C_regcoeff_alpha = coffee_dcvar_f(IDv);
        piaacmcparams.linopt_piaa1C_regcoeff_alpha = coffee_dcvar_f(IDv);
    }

    piaacmcparams.linopt_piaa0F_regcoeff = 0.0e-7f;
    piaacmcparams.linopt_piaa1F_regcoeff = 0.0e-7f;
    if ((IDv = variable_ID("REGPIAA_F_COEFF")) != -1)
    {
        piaacmcparams.linopt_piaa0F_regcoeff = coffee_dcvar_f(IDv);
        piaacmcparams.linopt_piaa1F_regcoeff = coffee_dcvar_f(IDv);
    }

    piaacmcparams.linopt_piaa0F_regcoeff_alpha = 1.0f;
    piaacmcparams.linopt_piaa1F_regcoeff_alpha = 1.0f;
    if ((IDv = variable_ID("REGPIAA_F_ALPHA")) != -1)
    {
        piaacmcparams.linopt_piaa0F_regcoeff_alpha = coffee_dcvar_f(IDv);
        piaacmcparams.linopt_piaa1F_regcoeff_alpha = coffee_dcvar_f(IDv);
    }

    if (piaacmcparams.linopt_REGPIAASHAPES == 1)
    {
        printf("loading CPA modes frequ\n");
        fflush(stdout);
        imageID ID_CPAfreq = image_ID("cpamodesfreq", dcimg, dcnimg);
        if (ID_CPAfreq == -1)
        {
            load_fits("cpamodesfreq.fits", "cpamodesfreq", LOADFITS_ERRMODE_EXIT, &ID_CPAfreq);
        }
    }

    piaacmcparams.linopt_REGFPMSAG = 1;
    if ((IDv = variable_ID("REGFPMSAG")) != -1)
    {
        piaacmcparams.linopt_REGFPMSAG = (long) coffee_dcvar_f(IDv) + 0.01;
    }
}

/**
 * opt_piaa_fpm_register_fpm_params - Register idealized FPM transmission or physical zone sags.
 */
static void opt_piaa_fpm_register_fpm_params(void)
{
    piaacmcparams.linopt_number_param = 0;

    if (piaacmcparams.PIAACMC_fpmtype == 0)
    {
        long idx = piaacmcparams.linopt_number_param;

        piaacmcparams.linopt_paramtype[idx]    = _DATATYPE_DOUBLE;
        piaacmcparams.linopt_paramval[idx]     = &dcimg[piaacmcopticaldesign.zoneaID].array.D[0];
        piaacmcparams.linopt_paramdelta[idx]   = 1.0e-3;
        piaacmcparams.linopt_parammaxstep[idx] = 2.0e-1;
        piaacmcparams.linopt_parammin[idx]     = -1.0;
        piaacmcparams.linopt_parammax[idx]     = 1.0;
        piaacmcparams.linopt_number_param++;
        return;
    }

    if (variable_ID("PIAACMC_mzOPT") != -1)
    {
        uint32_t nbzones = dcimg[piaacmcopticaldesign.zonezID].md[0].size[0];
        for (uint32_t mz = 0; mz < nbzones; mz++)
        {
            long idx = piaacmcparams.linopt_number_param;

            piaacmcparams.linopt_paramtype[idx]    = _DATATYPE_DOUBLE;
            piaacmcparams.linopt_paramval[idx]     =
                &dcimg[piaacmcopticaldesign.zonezID].array.D[mz];
            piaacmcparams.linopt_paramdelta[idx]   = 1.0e-9;
            piaacmcparams.linopt_parammaxstep[idx] = 5.0e-8;
            piaacmcparams.linopt_parammin[idx]     = -2.0e-6;
            piaacmcparams.linopt_parammax[idx]     = 2.0e-6;
            piaacmcparams.linopt_number_param++;
        }
    }
}

/**
 * opt_piaa_fpm_append_float_modes - Append float modal coefficients to linear optimizer state.
 * @modesID: Coefficient vector image ID
 * @kmax:    Number of leading modes to optimize
 */
static void opt_piaa_fpm_append_float_modes(
    imageID modesID,
    long    kmax)
{
    for (long k = 0; k < kmax; k++)
    {
        long idx = piaacmcparams.linopt_number_param;

        piaacmcparams.linopt_paramtype[idx]    = _DATATYPE_FLOAT;
        piaacmcparams.linopt_paramvalf[idx]    = &dcimg[modesID].array.F[k];
        piaacmcparams.linopt_paramdelta[idx]   = 1.0e-10;
        piaacmcparams.linopt_parammaxstep[idx] = 1.0e-7;
        piaacmcparams.linopt_parammin[idx]     = -1.0e-3;
        piaacmcparams.linopt_parammax[idx]     = 1.0e-3;
        piaacmcparams.linopt_number_param++;
    }
}

/**
 * exec_optimize_PIAA_shapes_fpmtransm - Set up linear optimization of PIAA modes and FPM (Mode 40).
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
errno_t exec_optimize_PIAA_shapes_fpmtransm(void)
{
    DEBUG_TRACE_FSTART();

    printf("=================================== mode 040 ===================================\n");

    FUNC_CHECK_RETURN(opt_piaa_fpm_init_design());

    piaacmcparams.LINOPT = 1;

    variableID IDv;
    if ((IDv = variable_ID("PIAACMC_nbiter")) != -1)
    {
        piaacmcparams.linopt_NBiter = (long) coffee_dcvar_f(IDv) + 0.01;
    }
    else
    {
        piaacmcparams.linopt_NBiter = 1000;
    }

    long maxCavail = dcimg[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0];
    long kmaxC     = maxCavail;
    if ((IDv = variable_ID("PIAACMC_maxoptCterm")) != -1)
    {
        kmaxC = (long) coffee_dcvar_f(IDv) + 0.01;
    }
    if (kmaxC > maxCavail)
    {
        kmaxC = maxCavail;
    }

    long maxFavail = dcimg[piaacmcopticaldesign.piaa0FmodesID].md[0].size[0];
    long kmaxF     = maxFavail;
    if ((IDv = variable_ID("PIAACMC_maxoptFterm")) != -1)
    {
        kmaxF = (long) coffee_dcvar_f(IDv) + 0.01;
    }
    if (kmaxF > maxFavail)
    {
        kmaxF = maxFavail;
    }

    opt_piaa_fpm_configure_regularization();
    opt_piaa_fpm_register_fpm_params();

    opt_piaa_fpm_append_float_modes(piaacmcopticaldesign.piaa0CmodesID, kmaxC);
    opt_piaa_fpm_append_float_modes(piaacmcopticaldesign.piaa1CmodesID, kmaxC);
    opt_piaa_fpm_append_float_modes(piaacmcopticaldesign.piaa0FmodesID, kmaxF);
    opt_piaa_fpm_append_float_modes(piaacmcopticaldesign.piaa1FmodesID, kmaxF);

    piaacmcparams.FORCE_MAKE_PIAA0shape = 1;
    piaacmcparams.FORCE_MAKE_PIAA1shape = 1;

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
