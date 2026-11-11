/**
 * @file    exec_optimize_lyot_stop_position.c
 * @brief   Optimize Lyot stop(s) conjugations (Mode 1)
 */

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_memory/COREMOD_memory.h"

#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"
#include "PIAACMCsimul_computePSF.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"
#include "init_piaacmcopticaldesign.h"
#include "init_piaacmcopticalsystem.h"

#include "LyotStop/exec_optimize_lyot_stop_position.h"
#include "PIAAshape/makePIAAshapes.h"

struct mode1_march_step
{
    long        ls;
    double      ref_pos;
    double      range;
    double      stepsize;
    const char *fnamelog;
};

/**
 * init_mode1_design_and_range - Load CLI options, initialize optical system, and read scan range.
 * @range_out: Output initial half-range [m] for Lyot stop marching scan
 *
 * Return: RETURN_SUCCESS or optical initialization error code.
 */
static errno_t init_mode1_design_and_range(double *range_out)
{
    imageID IDv;
    double  fpmradld = 0.95;
    double  centobs0 = 0.3;
    double  centobs1 = 0.2;

    printf("=================================== mode 001 ===================================\n");
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

    uint64_t initflag = INIT_PIAACMCOPTICALDESIGN_MODE__READCONF;
    initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF;
    FUNC_CHECK_RETURN(
        init_piaacmcopticaldesign(fpmradld, centobs0, centobs1, initflag, NULL));

    FUNC_CHECK_RETURN(makePIAAshapes());
    piaacmcopticalsystem.FOCMASKarray[0].mode = 1;
    FUNC_CHECK_RETURN(init_piaacmcopticalsystem(0.0, 0.0));

    variableID var_range = variable_ID("PIAACMC_lsoptrange");
    if (var_range != -1)
    {
        *range_out = coffee_dcvar_f(var_range);
    }
    else
    {
        *range_out = 3.0;
    }

    return RETURN_SUCCESS;
}

/**
 * find_lyot_mask0_elem - Locate `"Lyot mask 0"` in the optical system and mark it keepMem.
 *
 * Return: Element index of `"Lyot mask 0"`.
 */
static long find_lyot_mask0_elem(void)
{
    long elem0 = -1;

    printf("Number of elements = %ld\n", piaacmcopticalsystem.NBelem);
    assert(piaacmcopticalsystem.NBelem > 0);

    for (long elem = 0; elem < piaacmcopticalsystem.NBelem; elem++)
    {
        printf("elem %ld :  %s\n", elem, piaacmcopticalsystem.name[elem]);
        if (strcmp("Lyot mask 0", piaacmcopticalsystem.name[elem]) == 0)
        {
            elem0 = elem;
        }
    }
    assert(elem0 != -1);

    piaacmcopticalsystem.keepMem[elem0] = 1;
    return elem0;
}

/**
 * scan_single_lyot_stop_zpos - March one Lyot stop over `[ref - range, ref + range)` for PSF.
 * @step:         March parameters (stop index, reference z, half-range, step size, log path)
 * @best_pos_out: Output axial conjugation [m] achieving minimum evaluation-zone contrast
 *
 * Return: RETURN_SUCCESS or PSF evaluation error status.
 */
static errno_t scan_single_lyot_stop_zpos(
    const struct mode1_march_step *step,
    double                        *best_pos_out)
{
    long   ls       = step->ls;
    double best_pos = step->ref_pos - step->range;
    double valbest  = 1.0;

    piaacmcopticaldesign.LyotStop_zpos[ls] = best_pos;

    while (piaacmcopticaldesign.LyotStop_zpos[ls] < step->ref_pos + step->range)
    {
        long   elem0 = find_lyot_mask0_elem();
        double val   = 0.0;

        FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(0.0,
                                                  0.0,
                                                  elem0,
                                                  piaacmcopticalsystem.NBelem,
                                                  0,
                                                  0,
                                                  0,
                                                  0,
                                                  &val));

        if (val < valbest)
        {
            best_pos = piaacmcopticaldesign.LyotStop_zpos[ls];
            valbest  = val;
        }

        FILE *fp = fopen(step->fnamelog, "a");
        if (fp != NULL)
        {
            for (long ls1 = 0; ls1 < piaacmcopticaldesign.NBLyotStop; ls1++)
            {
                fprintf(fp, " %lf", piaacmcopticaldesign.LyotStop_zpos[ls1]);
            }
            fprintf(fp, " %g\n", val);
            fclose(fp);
        }

        piaacmcopticaldesign.LyotStop_zpos[ls] += step->stepsize;
    }

    printf("BEST SOLUTION :  ");
    piaacmcopticaldesign.LyotStop_zpos[ls] = best_pos;
    printf(" %lf", best_pos);
    printf(" %g\n", valbest);

    *best_pos_out = best_pos;
    return RETURN_SUCCESS;
}

/**
 * exec_optimize_lyot_stop_position - Mode 1 iterative refined marching search for Lyot z positions.
 *
 * Varies `piaacmcopticaldesign.LyotStop_zpos` across multi-scale axial scans to minimize residual
 * focal-plane contrast returned by `PIAACMCsimul_computePSF`.
 *
 * Return: RETURN_SUCCESS on completion, or milk error code on failure.
 */
errno_t exec_optimize_lyot_stop_position(void)
{
    DEBUG_TRACE_FSTART();

    double parambest[10000];
    double paramref[10000];
    char   fnamelog[STRINGMAXLEN_FULLFILENAME];
    double range = 3.0;

    FUNC_CHECK_RETURN(init_mode1_design_and_range(&range));

    for (long ls = 0; ls < piaacmcopticaldesign.NBLyotStop; ls++)
    {
        paramref[ls]  = piaacmcopticaldesign.LyotStop_zpos[ls];
        parambest[ls] = paramref[ls];
    }

    WRITE_FULLFILENAME(fnamelog, "%s/result_LMpos.log", piaacmcparams.piaacmcconfdir);
    {
        FILE *fp = fopen(fnamelog, "w");
        if (fp != NULL)
        {
            fclose(fp);
        }
    }

    long   NBiter   = 4;
    double stepsize = range / 5.0;

    for (long iter = 0; iter < NBiter; iter++)
    {
        for (long ls = 0; ls < piaacmcopticaldesign.NBLyotStop; ls++)
        {
            struct mode1_march_step step = {ls, paramref[ls], range, stepsize, fnamelog};
            FUNC_CHECK_RETURN(scan_single_lyot_stop_zpos(&step, &parambest[ls]));
            paramref[ls] = parambest[ls];
        }

        FILE *fp = fopen(fnamelog, "a");
        if (fp != NULL)
        {
            fprintf(fp, "\n");
            fclose(fp);
        }

        range *= 0.3;
        stepsize = range / 3.0;
    }

    for (long ls = 0; ls < piaacmcopticaldesign.NBLyotStop; ls++)
    {
        piaacmcopticaldesign.LyotStop_zpos[ls] = parambest[ls];
    }
    PIAACMCsimul_savepiaacmcconf(piaacmcparams.piaacmcconfdir);

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
