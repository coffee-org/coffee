/**
 * @file    PIAACMCsimul_conf_io.c
 * @brief   Configuration parameter file and CLI variable ingestion for PIAACMCsimul
 */

#define DEBUGLOG

#include <stdio.h>
#include <string.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_memory/COREMOD_memory.h"
#include "OpticsMaterials/OpticsMaterials.h"

#include "PIAACMCsimul.h"
#include "PIAACMCsimul_conf_io.h"

/**
 * init_design_defaults - Populate baseline optical design defaults and mode flags.
 * @fpmradld: Nominal idealized focal plane mask radius [lambda/D]
 * @centobs0: Entrance pupil central obscuration ratio [beamrad]
 * @centobs1: Output pupil central obscuration ratio [beamrad]
 * @WFCmode:  Number of active wavefront control deformable mirrors
 *
 * Return: RETURN_SUCCESS or failure if `conf/conf_peakPSF.txt` is malformed.
 */
static errno_t init_design_defaults(
    double fpmradld,
    double centobs0,
    double centobs1,
    int    WFCmode)
{
    piaacmcopticaldesign.nblambda = 8;
    piaacmcopticaldesign.size     = 2048;
    piaacmcopticaldesign.pixscale = 0.000055;

    piaacmcopticaldesign.beamrad   = 0.01;
    piaacmcopticaldesign.PIAA0pos  = 1.0;
    piaacmcopticaldesign.PIAAsep   = 1.00;
    piaacmcopticaldesign.fpzfactor = 8.0;
    piaacmcopticaldesign.Fratio    = 80.0;
    strcpy(piaacmcopticaldesign.PIAAmaterial_name, "Mirror");
    piaacmcopticaldesign.prePIAA0mask     = 0;
    piaacmcopticaldesign.prePIAA0maskpos  = 0.0;
    piaacmcopticaldesign.postPIAA0mask    = 0;
    piaacmcopticaldesign.postPIAA0maskpos = 0.0;
    piaacmcopticaldesign.piaaNBCmodesmax  = 40;
    piaacmcopticaldesign.piaaCPAmax       = 10.0;

    piaacmcopticaldesign.centObs0 = centobs0;
    piaacmcopticaldesign.centObs1 = centobs1;
    piaacmcopticaldesign.NBradpts = 50000;
    piaacmcopticaldesign.r0lim    = 1.15;
    piaacmcopticaldesign.r1lim    = 1.5;

    piaacmcopticaldesign.nbDM = WFCmode;
    for (long iDM = 0; iDM < piaacmcopticaldesign.nbDM; iDM++)
    {
        piaacmcopticaldesign.DMpos[iDM] =
            0.0 + 0.6 * iDM / (0.01 + piaacmcopticaldesign.nbDM - 1.0);
        piaacmcopticaldesign.ID_DM[iDM] = -1;
    }

    piaacmcopticaldesign.NBLyotStop = 2;
    for (int i = 0; i < 10; i++)
    {
        piaacmcopticaldesign.LyotStop_zpos[i] = 0.0;
        piaacmcopticaldesign.IDLyotStop[i]    = -1;
    }

    piaacmcopticaldesign.fpmaskradld       = fpmradld;
    piaacmcopticaldesign.fpmarraysize      = 2048;
    piaacmcopticaldesign.fpmRad            = 100.0e-6;
    piaacmcopticaldesign.NBrings           = 4;
    piaacmcopticaldesign.fpmminsag         = -1.0e-5;
    piaacmcopticaldesign.fpmmaxsag         = 1.0e-5;
    piaacmcopticaldesign.fpmsagreg_coeff   = 1.0;
    piaacmcopticaldesign.fpmsagreg_alpha   = 1.0;
    piaacmcopticaldesign.NBringCentCone    = 0;
    piaacmcopticaldesign.fpmCentConeZ      = piaacmcopticaldesign.fpmmaxsag;
    piaacmcopticaldesign.fpmOuterConeZ     = piaacmcopticaldesign.fpmmaxsag;
    piaacmcopticaldesign.fpmOuterConeRadld = 80.0;
    piaacmcopticaldesign.fpmmaterial_code  = 0;
    piaacmcopticaldesign.fpmaskamptransm   = 1.0;

    {
        FILE *fp = fopen("conf/conf_peakPSF.txt", "r");
        if (fp != NULL)
        {
            DEBUG_TRACEPOINT("File conf/conf_peakPSF.txt open for reading");
            if (fscanf(fp, "%f", &piaacmcopticaldesign.peakPSF) != 1)
            {
                fclose(fp);
                FUNC_RETURN_FAILURE("Cannot read value from file conf/conf_peakPSF.txt");
            }
            fclose(fp);
        }
        else
        {
            piaacmcopticaldesign.peakPSF = -1.0;
        }
        DEBUG_TRACEPOINT("piaacmcopticaldesign.peakPSF = %f", piaacmcopticaldesign.peakPSF);
    }

    piaacmcopticaldesign.PIAAmode = 1;
    {
        variableID IDv = variable_ID("PIAACMC_PIAAmode");
        if (IDv != -1)
        {
            piaacmcopticaldesign.PIAAmode = (int) (coffee_dcvar_f(IDv) + 0.01);
        }
    }
    DEBUG_TRACEPOINT("piaacmcopticaldesign.PIAAmode = %d", piaacmcopticaldesign.PIAAmode);

    piaacmcopticaldesign.PIAAcoeff = 1.0;
    {
        variableID IDv = variable_ID("PIAACMC_PIAAcoeff");
        if (IDv != -1)
        {
            piaacmcopticaldesign.PIAAcoeff = coffee_dcvar_f(IDv);
        }
    }
    DEBUG_TRACEPOINT("piaacmcopticaldesign.PIAAcoeff = %f", piaacmcopticaldesign.PIAAcoeff);

    if (piaacmcopticaldesign.PIAAmode == 0)
    {
        piaacmcopticaldesign.invPIAAmode = 0;
    }
    else
    {
        piaacmcopticaldesign.invPIAAmode = 1;
        variableID IDv                   = variable_ID("PIAACMC_invPIAAmode");
        if (IDv != -1)
        {
            piaacmcopticaldesign.invPIAAmode = (long) (coffee_dcvar_f(IDv) + 0.001);
        }
    }
    DEBUG_TRACEPOINT("piaacmcopticaldesign.invPIAAmode = %d", piaacmcopticaldesign.invPIAAmode);

    return RETURN_SUCCESS;
}

/**
 * read_fpm_material - Read or initialize FPM material name and write material code file.
 *
 * Return: RETURN_SUCCESS or failure if material files cannot be read/created.
 */
static errno_t read_fpm_material(void)
{
    FILE *fp;
    char  fname[STRINGMAXLEN_FULLFILENAME];

    WRITE_FULLFILENAME(fname, "%s/conf_fpmmaterial_name.txt", piaacmcparams.piaacmcconfdir);
    fp = fopen(fname, "r");
    if (fp != NULL)
    {
        char name[200];
        if (fscanf(fp, "%199s", name) == 1)
        {
            strcpy(piaacmcopticaldesign.fpmmaterial_name, name);
            printf("Reading %s   piaacmcopticaldesign.fpmmaterial_name : %s\n",
                   fname,
                   piaacmcopticaldesign.fpmmaterial_name);
        }
        else
        {
            fclose(fp);
            FUNC_RETURN_FAILURE("Cannot read value from file %s", fname);
        }
        fclose(fp);
    }
    else
    {
        snprintf(piaacmcopticaldesign.fpmmaterial_name,
                 STRINGMAXLEN_PIAACMCSIMUL_MATERIALNAME,
                 "Mirror");
        printf("Writing %s   piaacmcopticaldesign.fpmmaterial_name : %s\n",
               fname,
               piaacmcopticaldesign.fpmmaterial_name);
        fp = fopen(fname, "w");
        if (fp != NULL)
        {
            fprintf(fp, "%s\n", piaacmcopticaldesign.fpmmaterial_name);
            fclose(fp);
        }
        else
        {
            FUNC_RETURN_FAILURE("Cannot create file \"%s\"", fname);
        }
    }

    printf("piaacmcopticaldesign.fpmmaterial_name : %s\n", piaacmcopticaldesign.fpmmaterial_name);
    piaacmcopticaldesign.fpmmaterial_code =
        OpticsMaterials_code(piaacmcopticaldesign.fpmmaterial_name);

    WRITE_FULLFILENAME(fname, "%s/conf_fpmmaterial_code.txt", piaacmcparams.piaacmcconfdir);
    fp = fopen(fname, "w");
    if (fp != NULL)
    {
        fprintf(fp, "%d\n", piaacmcopticaldesign.fpmmaterial_code);
        fclose(fp);
    }
    else
    {
        FUNC_RETURN_FAILURE("Cannot create file \"%s\"", fname);
    }
    DEBUG_TRACEPOINT("piaacmcopticaldesign.fpmmaterial_code = %d",
                     piaacmcopticaldesign.fpmmaterial_code);

    return RETURN_SUCCESS;
}

/**
 * read_dcvar_optics_params - Ingest pupil, PIAA, and Lyot geometry CLI variables.
 */
static void read_dcvar_optics_params(void)
{
    variableID IDv;

    if ((IDv = variable_ID("PIAACMC_beamrad")) != -1)
    {
        piaacmcopticaldesign.beamrad = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_Fratio")) != -1)
    {
        piaacmcopticaldesign.Fratio = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_r0lim")) != -1)
    {
        piaacmcopticaldesign.r0lim = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_r1lim")) != -1)
    {
        piaacmcopticaldesign.r1lim = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_PIAAsep")) != -1)
    {
        piaacmcopticaldesign.PIAAsep = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_PIAA0pos")) != -1)
    {
        piaacmcopticaldesign.PIAA0pos = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_prePIAA0maskpos")) != -1)
    {
        piaacmcopticaldesign.prePIAA0maskpos = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_postPIAA0maskpos")) != -1)
    {
        piaacmcopticaldesign.postPIAA0maskpos = coffee_dcvar_f(IDv);
    }

    piaacmcopticaldesign.LyotZmin = -3.0;
    if ((IDv = variable_ID("PIAACMC_LyotZmin")) != -1)
    {
        piaacmcopticaldesign.LyotZmin = coffee_dcvar_f(IDv);
    }
    piaacmcopticaldesign.LyotZmax = 3.0;
    if ((IDv = variable_ID("PIAACMC_LyotZmax")) != -1)
    {
        piaacmcopticaldesign.LyotZmax = coffee_dcvar_f(IDv);
    }

    piaacmcopticaldesign.pupoutmaskrad = 0.95;
    if ((IDv = variable_ID("PIAACMC_pupoutmaskrad")) != -1)
    {
        piaacmcopticaldesign.pupoutmaskrad = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_piaaNBCmodesmax")) != -1)
    {
        piaacmcopticaldesign.piaaNBCmodesmax = (long) (coffee_dcvar_f(IDv) + 0.01);
    }
    if ((IDv = variable_ID("PIAACMC_piaaCPAmax")) != -1)
    {
        piaacmcopticaldesign.piaaCPAmax = coffee_dcvar_f(IDv);
    }

    piaacmcopticaldesign.NBLyotStop = 1;
    if (piaacmcopticaldesign.PIAAmode != 0 && (IDv = variable_ID("PIAACMC_nblstop")) != -1)
    {
        piaacmcopticaldesign.NBLyotStop = (long) coffee_dcvar_f(IDv) + 0.01;
    }

    if ((IDv = variable_ID("PIAACMC_size")) != -1)
    {
        piaacmcopticaldesign.size = (long) (coffee_dcvar_f(IDv) + 0.01);
    }
    if ((IDv = variable_ID("PIAACMC_pixscale")) != -1)
    {
        piaacmcopticaldesign.pixscale = coffee_dcvar_f(IDv);
    }
}

/**
 * read_dcvar_spectral_fpm_params - Ingest bandpass and focal plane mask CLI variables.
 * @piaacmctype: 0 for idealized monochromatic FPM, non-zero for physical FPM
 * @fpmradld:    Nominal idealized focal plane mask radius [lambda/D]
 */
static void read_dcvar_spectral_fpm_params(
    long   piaacmctype,
    double fpmradld)
{
    variableID IDv;

    if ((IDv = variable_ID("PIAACMC_lambda")) != -1)
    {
        piaacmcopticaldesign.lambda = 1.0e-9 * coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_lambdaB")) != -1)
    {
        piaacmcopticaldesign.lambdaB = coffee_dcvar_f(IDv);
    }

    piaacmcparams.LAMBDASTART =
        piaacmcopticaldesign.lambda * (1.0 - 0.005 * piaacmcopticaldesign.lambdaB);
    piaacmcparams.LAMBDAEND =
        piaacmcopticaldesign.lambda * (1.0 + 0.005 * piaacmcopticaldesign.lambdaB);

    if ((IDv = variable_ID("PIAACMC_nblambda")) != -1)
    {
        piaacmcopticaldesign.nblambda = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_NBrings")) != -1)
    {
        piaacmcopticaldesign.NBrings = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_fpmminsag")) != -1)
    {
        piaacmcopticaldesign.fpmminsag = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_fpmmaxsag")) != -1)
    {
        piaacmcopticaldesign.fpmmaxsag = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_fpmsagreg_coeff")) != -1)
    {
        piaacmcopticaldesign.fpmsagreg_coeff = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_fpmsagreg_alpha")) != -1)
    {
        piaacmcopticaldesign.fpmsagreg_alpha = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_NBringCentCone")) != -1)
    {
        piaacmcopticaldesign.NBringCentCone = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_fpmCentConeZ")) != -1)
    {
        piaacmcopticaldesign.fpmCentConeZ = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_fpmOuterConeZ")) != -1)
    {
        piaacmcopticaldesign.fpmOuterConeZ = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_fpmOuterConeRadld")) != -1)
    {
        piaacmcopticaldesign.fpmOuterConeRadld = coffee_dcvar_f(IDv);
    }

    double lam_mid = 0.5 * (piaacmcparams.LAMBDASTART + piaacmcparams.LAMBDAEND);
    piaacmcopticaldesign.fpmOuterConeRad =
        lam_mid * piaacmcopticaldesign.Fratio * piaacmcopticaldesign.fpmOuterConeRadld;

    DEBUG_TRACEPOINT("piaacmctype = %ld", piaacmctype);
    if (piaacmctype == 0)
    {
        piaacmcparams.FORCE_CREATE_fpmzt     = 1;
        piaacmcopticaldesign.NBrings         = 1;
        piaacmcopticaldesign.NBringCentCone  = 0;
        piaacmcopticaldesign.fpmOuterConeZ   = 0.0;
        piaacmcopticaldesign.fpmminsag       = -1e-5;
        piaacmcopticaldesign.fpmmaxsag       = 1e-5;
        piaacmcopticaldesign.fpmsagreg_coeff = 1.0;
        piaacmcopticaldesign.fpmsagreg_alpha = 1.0;
        piaacmcopticaldesign.fpmRad          = lam_mid * piaacmcopticaldesign.Fratio * fpmradld;
        printf("Idealized focal plane mask  radius = %f l/D  = %g m    "
               "[lambda = %g (%g - %g)] [Fratio = %f]\n",
               fpmradld,
               piaacmcopticaldesign.fpmRad,
               lam_mid,
               piaacmcparams.LAMBDASTART,
               piaacmcparams.LAMBDAEND,
               piaacmcopticaldesign.Fratio);
    }
    else
    {
        if (piaacmcparams.PIAACMC_MASKRADLD < 0.2)
        {
            piaacmcparams.PIAACMC_MASKRADLD = 0.1 * ((long) (10.0 * 1.2 * fpmradld));
        }
        piaacmcopticaldesign.fpmRad =
            lam_mid * piaacmcopticaldesign.Fratio * piaacmcparams.PIAACMC_MASKRADLD;
        printf("Physical focal plane mask - rad = %f l/D -> %g m   "
               "[lambda = %g (%g - %g)] [Fratio = %f]\n",
               piaacmcparams.PIAACMC_MASKRADLD,
               piaacmcopticaldesign.fpmRad,
               lam_mid,
               piaacmcparams.LAMBDASTART,
               piaacmcparams.LAMBDAEND,
               piaacmcopticaldesign.Fratio);
    }
}

/**
 * PIAACMCsimul_initpiaacmcconf_readconfparams - Read defaults, CLI vars, and conf files.
 * @piaacmctype: 0 for idealized monochromatic FPM, non-zero for physical FPM
 * @fpmradld:    Nominal idealized focal plane mask radius [lambda/D]
 * @centobs0:    Entrance pupil central obscuration ratio [beamrad]
 * @centobs1:    Output pupil central obscuration ratio [beamrad]
 * @WFCmode:     Number of active deformable mirrors
 *
 * Return: RETURN_SUCCESS on success, or error code on I/O failure.
 */
errno_t PIAACMCsimul_initpiaacmcconf_readconfparams(
    long   piaacmctype,
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

    FUNC_CHECK_RETURN(init_design_defaults(fpmradld, centobs0, centobs1, WFCmode));
    FUNC_CHECK_RETURN(read_fpm_material());
    read_dcvar_optics_params();
    read_dcvar_spectral_fpm_params(piaacmctype, fpmradld);

    piaacmcopticaldesign.CmodesID      = -1;
    piaacmcopticaldesign.FmodesID      = -1;
    piaacmcopticaldesign.piaa0CmodesID = -1;
    piaacmcopticaldesign.piaa0FmodesID = -1;
    piaacmcopticaldesign.piaa1CmodesID = -1;
    piaacmcopticaldesign.piaa1FmodesID = -1;
    piaacmcopticaldesign.zonezID       = -1;
    piaacmcopticaldesign.zoneaID       = -1;

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * PIAACMCsimul_read_piaa_material - Read or initialize PIAA material name and code files.
 *
 * Return: RETURN_SUCCESS or failure if material files cannot be read/created.
 */
errno_t PIAACMCsimul_read_piaa_material(void)
{
    DEBUG_TRACE_FSTART();

    FILE *fp;
    char  fname[STRINGMAXLEN_FULLFILENAME];

    WRITE_FULLFILENAME(fname, "%s/conf_PIAAmaterial_name.txt", piaacmcparams.piaacmcconfdir);
    fp = fopen(fname, "r");
    if (fp != NULL)
    {
        char name[200];
        if (fscanf(fp, "%199s", name) == 1)
        {
            strcpy(piaacmcopticaldesign.PIAAmaterial_name, name);
        }
        else
        {
            fclose(fp);
            FUNC_RETURN_FAILURE("Cannot read value from file %s", fname);
        }
        fclose(fp);
    }
    else
    {
        snprintf(piaacmcopticaldesign.PIAAmaterial_name,
                 STRINGMAXLEN_PIAACMCSIMUL_MATERIALNAME,
                 "Mirror");
        fp = fopen(fname, "w");
        if (fp != NULL)
        {
            fprintf(fp, "%s\n", piaacmcopticaldesign.PIAAmaterial_name);
            fclose(fp);
        }
        else
        {
            FUNC_RETURN_FAILURE("Cannot create file \"%s\"", fname);
        }
    }

    DEBUG_TRACEPOINT("piaacmcopticaldesign.PIAAmaterial_name = %s",
                     piaacmcopticaldesign.PIAAmaterial_name);
    piaacmcopticaldesign.PIAAmaterial_code =
        OpticsMaterials_code(piaacmcopticaldesign.PIAAmaterial_name);
    DEBUG_TRACEPOINT("piaacmcopticaldesign.PIAAmaterial_code = %d",
                     piaacmcopticaldesign.PIAAmaterial_code);

    WRITE_FULLFILENAME(fname, "%s/conf_PIAAmaterial_code.txt", piaacmcparams.piaacmcconfdir);
    fp = fopen(fname, "w");
    if (fp != NULL)
    {
        fprintf(fp, "%d\n", piaacmcopticaldesign.PIAAmaterial_code);
        fclose(fp);
    }
    else
    {
        FUNC_RETURN_FAILURE("cannot create file \"%s\"", fname);
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
