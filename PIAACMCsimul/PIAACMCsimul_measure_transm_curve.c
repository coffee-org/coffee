/**
 * @file    PIAACMCsimul_measure_transm_curve.c
 * @brief   Measure coronagraphic off-axis transmission as a function of angular separation
 */

#include <stdint.h>
#include <stdio.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"
#include "PIAACMCsimul_computePSF.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"
#include "PIAACMCsimul_measure_transm_curve.h"
#include "init_piaacmcopticaldesign.h"

#include "PIAAshape/makePIAAshapes.h"

/**
 * measure_transm_init_design - Load CLI overrides and initialize the PIAACMC optical design.
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t measure_transm_init_design(void)
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
 * measure_transm_sample - Propagate one off-axis separation and append core flux to curve file.
 * @xld:         Angular source separation along x [lambda/D]
 * @fnametransm: Output ASCII transmission curve path
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t measure_transm_sample(
    double      xld,
    const char *fnametransm)
{
    FUNC_CHECK_RETURN(
        PIAACMCsimul_computePSF(xld, 0.0, 0, piaacmcopticalsystem.NBelem, 0, 0, 0, 0, NULL));

    imageID  ID    = image_ID("psfi0", dcimg, dcnimg);
    uint32_t xsize = dcimg[ID].md[0].size[0];
    uint32_t ysize = dcimg[ID].md[0].size[1];
    uint32_t zsize = dcimg[ID].md[0].size[2];
    printf("ID = %ld\n", ID);
    printf("image size = %u %u %u\n", xsize, ysize, zsize);

    double val = 0.0;
    for (uint32_t kk = 0; kk < zsize; kk++)
    {
        for (uint32_t ii = 0; ii < xsize; ii++)
        {
            for (uint32_t jj = 0; jj < ysize; jj++)
            {
                double dx = 1.0 * ii - 0.5 * xsize;
                double dy = 1.0 * jj - 0.5 * ysize;
                if ((dx * dx + dy * dy) < 30.0 * 30.0)
                {
                    val += dcimg[ID].array.F[kk * xsize * ysize + jj * xsize + ii];
                }
            }
        }
    } // for kk
    val /= zsize;

    FILE *fpt = fopen(fnametransm, "a");
    if (fpt == NULL)
    {
        FUNC_RETURN_FAILURE("Cannot open transmission curve file");
    }
    fprintf(fpt, "%10f %.18f\n", xld, val);
    fclose(fpt);

    FUNC_CHECK_RETURN(delete_image_ID("psfi0", DELETE_IMAGE_ERRMODE_WARNING));
    return RETURN_SUCCESS;
}

/**
 * PIAACMCsimul_measure_transm_curve - Measure planetary throughput vs. separation (Mode 101).
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
errno_t PIAACMCsimul_measure_transm_curve(void)
{
    DEBUG_TRACE_FSTART();

    printf("=================================== mode 101 ===================================\n");

    FUNC_CHECK_RETURN(measure_transm_init_design());

    {
        double cval = 0.0;
        FUNC_CHECK_RETURN(
            PIAACMCsimul_computePSF(0.0, 0.0, 0, piaacmcopticalsystem.NBelem, 1, 0, 0, 0, &cval));
    }

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/psfi0test_x00_y00.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("psfi0", fname));
    }

    char fnametransm[STRINGMAXLEN_FULLFILENAME];
    PIAACMCsimul_update_fnamedescr();
    SNPRINTF_CHECK(
        fnametransm,
        STRINGMAXLEN_FULLFILENAME,
        "%s/transmCurve_sm%d.%s.txt",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.SCORINGMASKTYPE,
        piaacmcparams.fnamedescr);

    FILE *fpt = fopen(fnametransm, "w");
    if (fpt == NULL)
    {
        FUNC_RETURN_FAILURE("Cannot create transmission curve file");
    }
    fclose(fpt);

    double stepld = 0.001;
    for (double xld = 0.0; xld < 10.0; xld += stepld)
    {
        FUNC_CHECK_RETURN(measure_transm_sample(xld, fnametransm));

        stepld = 0.001 + 0.1 * xld;
        if (stepld > 0.2)
        {
            stepld = 0.2;
        }
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
