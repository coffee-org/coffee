/**
 * @file    PIAACMCsimul_exec.c
 * @brief   Primary mode dispatcher for PIAACMC and APLCMC design and evaluation steps
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_memory/COREMOD_memory.h"

#include "FocalPlaneMask/exec_multizone_fpm_calib.h"
#include "FocalPlaneMask/exec_optimize_fpm_zones.h"
#include "FocalPlaneMask/exec_optimize_fpmtransmission.h"
#include "LyotStop/exec_optimize_lyot_stop_position.h"
#include "LyotStop/exec_optimize_lyot_stops_shapes_positions.h"
#include "PIAACMCsimul.h"
#include "PIAACMCsimul_eval_poly_design.h"
#include "PIAACMCsimul_exec.h"
#include "PIAACMCsimul_linopt.h"
#include "PIAACMCsimul_measure_transm_curve.h"
#include "PIAAshape/exec_optimize_PIAA_shapes.h"
#include "exec_computePSF_no_fpm.h"
#include "exec_compute_image.h"
#include "exec_optimize_PIAA_shapes_fpmtransm.h"

/**
 * exec_attach_status_shm - Attach or create `stat_PIAACMCsimulexec` shared status image.
 * @IDstatus_out: Output image ID for status monitoring array
 *
 * Return: RETURN_SUCCESS or allocation error status.
 */
static errno_t exec_attach_status_shm(imageID *IDstatus_out)
{
    DEBUG_TRACE_FSTART();

    imageID IDstatus = image_ID("stat_PIAACMCsimulexec", dcimg, dcnimg);
    if (IDstatus == -1)
    {
        printf("Looking for stat_PIAACMCsimulexec\n");
        fflush(stdout);
        IDstatus = read_sharedmem_image("stat_PIAACMCsimulexec", dcimg, dcnimg);
        printf("ID for stat_PIAACMCsimulexec: %ld\n", IDstatus);
        fflush(stdout);
    }
    if (IDstatus == -1)
    {
        uint32_t *sizearray = (uint32_t *) malloc(sizeof(uint32_t) * 2);
        if (sizearray == NULL)
        {
            FUNC_RETURN_FAILURE("malloc returns NULL pointer");
        }
        sizearray[0] = 4;
        sizearray[1] = 1;
        create_image_ID("stat_PIAACMCsimulexec",
                        2,
                        sizearray,
                        _DATATYPE_UINT16,
                        1,
                        0,
                        0,
                        &IDstatus);
        free(sizearray);
    }

    *IDstatus_out = IDstatus;
    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * exec_dispatch_mode - Dispatch PIAACMC design/evaluation step handler for `mode`.
 * @mode: Operation code (`0..5`, `11`, `13`, `40`, `100`, `101`, `302`)
 *
 * Return: RETURN_SUCCESS or mode handler error status.
 */
static errno_t exec_dispatch_mode(long mode)
{
    DEBUG_TRACE_FSTART();

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
    case 100:
        FUNC_CHECK_RETURN(PIAACMCsimul_eval_poly_design());
        break;
    case 101:
        FUNC_CHECK_RETURN(PIAACMCsimul_measure_transm_curve());
        break;
    case 302:
        printf("=================================== mode 302 "
               "===================================\n");
        EXECUTE_SYSTEM_COMMAND("cp %s/saveconf/conf_*.txt %s/",
                               piaacmcparams.piaacmcconfdir,
                               piaacmcparams.piaacmcconfdir);
        break;
    default:
        PRINT_ERROR("mode not recognized");
        break;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * PIAACMCsimul_exec - Execute PIAACMC/APLCMC design or evaluation mode in `confindex`.
 * @confindex: Configuration directory path (e.g. `piaacmcconf_i000`)
 * @mode:      Operation mode code
 *
 * Return: RETURN_SUCCESS on completion, or error code on failure.
 */
errno_t PIAACMCsimul_exec(
    const char *confindex,
    long        mode)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %s %ld", confindex, mode);

    char       stopfile[STRINGMAXLEN_FULLFILENAME];
    imageID    IDstatus = -1;
    variableID IDv;

    piaacmcparams.linopt_REGPIAASHAPES         = 0;
    piaacmcparams.linopt_piaa0C_regcoeff       = 0.0e-7;
    piaacmcparams.linopt_piaa1C_regcoeff       = 0.0e-7;
    piaacmcparams.linopt_piaa0C_regcoeff_alpha = 1.0;
    piaacmcparams.linopt_piaa1C_regcoeff_alpha = 1.0;
    piaacmcparams.linopt_piaa0F_regcoeff       = 0.0e-7;
    piaacmcparams.linopt_piaa1F_regcoeff       = 0.0e-7;
    piaacmcparams.linopt_piaa0F_regcoeff_alpha = 1.0;
    piaacmcparams.linopt_piaa1F_regcoeff_alpha = 1.0;
    piaacmcparams.linopt_REGFPMSAG             = 0;

    FUNC_CHECK_RETURN(exec_attach_status_shm(&IDstatus));

    for (int elem = 0; elem < 100; elem++)
    {
        piaacmcopticalsystem.keepMem[elem] = 0;
    }

    snprintf(piaacmcparams.piaacmcconfdir, STRINGMAXLEN_DIRNAME, "%s", confindex);
    snprintf(dcsavedir, STRINGMAXLEN_DIRNAME, "%s", piaacmcparams.piaacmcconfdir);

    piaacmcparams.linopt_NBiter       = 1000;
    piaacmcparams.linopt_number_param = 0;
    piaacmcopticalsystem.SAVE         = piaacmcparams.PIAACMC_save;

    if ((IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        double fpmradld = coffee_dcvar_f(IDv);
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }

    EXECUTE_SYSTEM_COMMAND("echo \"%03ld     $(date)\" >> ./log/PIAACMC_mode_log.txt", mode);
    WRITE_FULLFILENAME(stopfile, "%s/stopmode%ld.txt", piaacmcparams.piaacmcconfdir, mode);

    FUNC_CHECK_RETURN(exec_dispatch_mode(mode));

    if (piaacmcparams.LINOPT == 1)
    {
        FUNC_CHECK_RETURN(PIAACMCsimul_run_linopt(IDstatus, stopfile));
    }

    piaacmcparams.LINOPT = 0;

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
