/**
 * @file    exec_optimize_fpm_zones.c
 * @brief   Optimize focal plane mask zone values (Mode 13)
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "statistic/statistic.h"

#include "OptSystProp/OptSystProp.h"
#include "OpticsMaterials/OpticsMaterials.h"

#include "PIAACMCsimul.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"
#include "init_piaacmcopticaldesign.h"
#include "init_piaacmcopticalsystem.h"

#include "FocalPlaneMask/exec_optimize_fpm_zones.h"
#include "PIAAshape/makePIAAshapes.h"

/**
 * read_cnorm_factor - Read contrast normalization factor from CnormFactor.txt
 * @fname: Path to CnormFactor.txt.
 * @val:   Output pointer receiving the parsed normalization factor.
 *
 * Skips leading blank lines and '#' header comments written by
 * PIAACMCsimul_computePSF().
 *
 * Return: RETURN_SUCCESS on valid parse, RETURN_FAILURE otherwise.
 */
static errno_t read_cnorm_factor(
    const char *restrict fname,
    double     *restrict val)
{
    FILE *fp = fopen(fname, "r");
    if (fp == NULL)
    {
        FUNC_RETURN_FAILURE("Cannot open %s", fname);
    }

    char line[256];
    int  parsed = 0;
    while (fgets(line, (int) sizeof(line), fp) != NULL)
    {
        const char *p = line;
        while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n')
        {
            p++;
        }
        if (*p == '\0' || *p == '#')
        {
            continue;
        }
        if (sscanf(p, "%lf", val) == 1)
        {
            parsed = 1;
        }
        break;
    }

    fclose(fp);
    if (parsed != 1)
    {
        FUNC_RETURN_FAILURE("Failed to parse CnormFactor from %s", fname);
    }

    return RETURN_SUCCESS;
}

/**
 * fpm_zones_init_design - Parse CLI overrides and initialize physical FPM optical system
 *
 * Return: RETURN_SUCCESS on initialization, error code otherwise.
 */
static errno_t fpm_zones_init_design(void)
{
    imageID IDv;
    double  fpmradld = 0.95;
    double  centobs0 = 0.3;
    double  centobs1 = 0.2;

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

    piaacmcparams.linopt_REGFPMSAG = 1;
    if ((IDv = variable_ID("REGFPMSAG")) != -1)
    {
        piaacmcparams.linopt_REGFPMSAG = (long) coffee_dcvar_f(IDv) + 0.01;
    }

    uint64_t initflag = INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT;
    initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__READCONF;
    initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF;
    initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL;

    FUNC_CHECK_RETURN(init_piaacmcopticaldesign(fpmradld, centobs0, centobs1, initflag, NULL));
    FUNC_CHECK_RETURN(makePIAAshapes());
    FUNC_CHECK_RETURN(init_piaacmcopticalsystem(0.0, 0.0));

    return RETURN_SUCCESS;
}

/**
 * fpm_zones_load_flux - Read per-element beam flux diagnostics from tmp_flux.txt
 *
 * Return: RETURN_SUCCESS on valid read across all elements, RETURN_FAILURE otherwise.
 */
static errno_t fpm_zones_load_flux(void)
{
    char fname[STRINGMAXLEN_FULLFILENAME];

    WRITE_FULLFILENAME(fname, "%s/tmp_flux.txt", piaacmcparams.piaacmcconfdir);
    FILE *fp = fopen(fname, "r");
    if (fp == NULL)
    {
        FUNC_RETURN_FAILURE("Cannot open tmp_flux.txt");
    }

    for (long elem = 0; elem < piaacmcopticalsystem.NBelem; elem++)
    {
        double tmplf1;
        double tmplf2;
        int    tmpd1;
        int    fscanfcnt = fscanf(fp, "%lf %lf  %d\n", &tmplf1, &tmplf2, &tmpd1);
        if (fscanfcnt == EOF)
        {
            if (ferror(fp))
            {
                perror("fscanf");
            }
            else
            {
                fprintf(stderr, "Error: fscanf reached end of file\n");
            }
            fclose(fp);
            FUNC_RETURN_FAILURE("Call to fscanf failed");
        }
        if (fscanfcnt != 3)
        {
            fprintf(stderr, "Error: fscanf matched %i items, 3 expected\n", fscanfcnt);
            fclose(fp);
            FUNC_RETURN_FAILURE("Call to fscanf failed");
        }

        piaacmcopticalsystem.flux[elem] = tmplf1 / tmpd1 * piaacmcopticalsystem.nblambda;
    }

    fclose(fp);
    return RETURN_SUCCESS;
}

/**
 * fpm_zones_setup_fast_arrays - Load FPMresp cube and allocate fast-evaluation buffers
 *
 * Return: RETURN_SUCCESS on successful load and allocation, error code otherwise.
 */
static errno_t fpm_zones_setup_fast_arrays(void)
{
    char    fname[STRINGMAXLEN_FULLFILENAME];
    imageID IDfpmresp;

    piaacmcparams.LINOPT = 1;

    PIAACMCsimul_update_fnamedescr_conf();
    SNPRINTF_CHECK(
        fname,
        STRINGMAXLEN_FULLFILENAME,
        "%s/FPMresp%d.%s.fits",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.SCORINGMASKTYPE,
        piaacmcparams.fnamedescr_conf);

    FUNC_CHECK_RETURN(load_fits(fname, "FPMresp", 1, &IDfpmresp));

    piaacmcparams.vsize = dcimg[IDfpmresp].md[0].size[0];
    FUNC_CHECK_RETURN(create_2Dimage_ID(
        "imvect1",
        piaacmcparams.vsize,
        piaacmcopticaldesign.nblambda,
        NULL));

    piaacmcparams.fpmresp_array = dcimg[IDfpmresp].array.D;
    piaacmcparams.zonez_array   = dcimg[piaacmcopticaldesign.zonezID].array.D;
    piaacmcparams.dphadz_array  = (double *) malloc(
        sizeof(double) * (size_t) piaacmcopticaldesign.nblambda);
    if (piaacmcparams.dphadz_array == NULL)
    {
        FUNC_RETURN_FAILURE("malloc failed for dphadz_array");
    }

    for (long k = 0; k < piaacmcopticaldesign.nblambda; k++)
    {
        piaacmcparams.dphadz_array[k] = OpticsMaterials_pha_lambda(
            piaacmcopticaldesign.fpmmaterial_code,
            1.0,
            piaacmcopticalsystem.lambdaarray[k]);
        printf(
            "%ld  %g %g\n",
            k,
            piaacmcopticalsystem.lambdaarray[k],
            piaacmcparams.dphadz_array[k]);
    }

    size_t out_nelem = (size_t) piaacmcparams.vsize * (size_t) piaacmcopticaldesign.nblambda
                       + (size_t) dcimg[piaacmcopticaldesign.zonezID].md[0].size[0];
    piaacmcparams.outtmp_array = (double *) malloc(sizeof(double) * out_nelem);
    if (piaacmcparams.outtmp_array == NULL)
    {
        FUNC_RETURN_FAILURE("malloc failed for outtmp_array");
    }

    piaacmcparams.computePSF_FAST_FPMresp = 1;
    return RETURN_SUCCESS;
}

/**
 * fpm_zones_register_linopt_params - Perturb zone sags and register linear optimizer state
 *
 * Return: RETURN_SUCCESS on completion, error code otherwise.
 */
static errno_t fpm_zones_register_linopt_params(void)
{
    char fname[STRINGMAXLEN_FULLFILENAME];

    WRITE_FULLFILENAME(fname, "%s/CnormFactor.txt", piaacmcparams.piaacmcconfdir);
    FUNC_CHECK_RETURN(read_cnorm_factor(fname, &piaacmcparams.CnormFactor));

    long nbzones = (long) dcimg[piaacmcopticaldesign.zonezID].md[0].size[0];
    for (long k = 0; k < nbzones; k++)
    {
        double offset = piaacmcparams.MODampl * (1.0 - 2.0 * ran1());
        dcimg[piaacmcopticaldesign.zonezID].array.D[k] += offset;
    }

    piaacmcparams.linopt_number_param = 0;
    for (long mz = 0; mz < nbzones; mz++)
    {
        long idx = piaacmcparams.linopt_number_param;

        piaacmcparams.linopt_paramtype[idx]    = _DATATYPE_DOUBLE;
        piaacmcparams.linopt_paramval[idx]     = &dcimg[piaacmcopticaldesign.zonezID].array.D[mz];
        piaacmcparams.linopt_paramdelta[idx]   = 3.0e-9;
        piaacmcparams.linopt_parammaxstep[idx] = 1.0e-6;
        piaacmcparams.linopt_parammin[idx]     = piaacmcopticaldesign.fpmminsag;
        piaacmcparams.linopt_parammax[idx]     = piaacmcopticaldesign.fpmmaxsag;
        piaacmcparams.linopt_number_param++;
    }

    piaacmcparams.PIAACMC_FPM_FASTDERIVATIVES = 1;
    return RETURN_SUCCESS;
}

/**
 * exec_optimize_fpm_zones - Configure fast linear optimization of multi-zone FPM sags (Mode 13)
 *
 * Loads the pre-calibrated polychromatic FPM zone impulse-response matrix (FPMresp),
 * precomputes analytic material phase-thickness derivatives, perturbs starting zone
 * sags, and registers per-zone optimization parameters for the outer LINOPT loop.
 *
 * Return: RETURN_SUCCESS on successful setup, error code otherwise.
 */
errno_t exec_optimize_fpm_zones(void)
{
    DEBUG_TRACE_FSTART();

    printf("=================================== mode 013 ===================================\n");

    FUNC_CHECK_RETURN(fpm_zones_init_design());
    FUNC_CHECK_RETURN(fpm_zones_load_flux());
    FUNC_CHECK_RETURN(fpm_zones_setup_fast_arrays());
    FUNC_CHECK_RETURN(fpm_zones_register_linopt_params());

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
