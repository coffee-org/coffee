/**
 * @file    PIAACMCsimul_run_params.c
 * @brief   CLI parameter loading and Mode 13 starting-point perturbation helpers
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_memory/COREMOD_memory.h"
#include "statistic/statistic.h"

#include "PIAACMCsimul.h"
#include "PIAACMCsimul_run_internal.h"
#include "init_piaacmcopticaldesign.h"

/**
 * read_run_cli_flags - Populate global piaacmcparams runtime flags from CLI variables.
 * @mode: Active execution mode index
 */
static void read_run_cli_flags(long mode)
{
    variableID IDv = -1;

    if ((IDv = variable_ID("PIAACMC_MASKRADLD")) != -1)
    {
        piaacmcparams.PIAACMC_MASKRADLD = coffee_dcvar_f(IDv);
    }

    if ((IDv = variable_ID("PIAACMC_FPMsectors")) != -1)
    {
        piaacmcparams.PIAACMC_FPMsectors = (long) coffee_dcvar_f(IDv) + 0.01;
    }
    printf("PIAACMC_FPMsectors = %d\n", piaacmcparams.PIAACMC_FPMsectors);

    if ((IDv = variable_ID("SCORINGMASKTYPE")) != -1)
    {
        piaacmcparams.SCORINGMASKTYPE = (long) coffee_dcvar_f(IDv) + 0.01;
    }
    printf("SCORINGMASKTYPE = %d\n", piaacmcparams.SCORINGMASKTYPE);

    if ((IDv = variable_ID("PIAACMC_save")) != -1)
    {
        piaacmcparams.PIAACMC_save = (long) coffee_dcvar_f(IDv) + 0.01;
    }
    printf("PIAACMC_save = %d\n", piaacmcparams.PIAACMC_save);

    if ((IDv = variable_ID("PIAACMC_resolved")) != -1)
    {
        piaacmcparams.computePSF_ResolvedTarget = (long) (coffee_dcvar_f(IDv) + 0.01);
    }

    if ((IDv = variable_ID("PIAACMC_extmode")) != -1)
    {
        piaacmcparams.computePSF_ResolvedTarget_mode = (long) (coffee_dcvar_f(IDv) + 0.01);
    }

    printf("mode = %ld\n", mode);
    fflush(stdout);
}

/**
 * PIAACMCsimul_run_setparam_variables - Read CLI variables and initialize piaacmcopticaldesign.
 * @confindex: Configuration directory name/index
 * @mode:      Execution mode index
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
errno_t PIAACMCsimul_run_setparam_variables(
    const char *confindex,
    long        mode)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %s %ld", confindex, mode);

    read_run_cli_flags(mode);

    double     fpmradld = 0.95;
    double     centobs0 = 0.3;
    double     centobs1 = 0.2;
    variableID IDv      = -1;

    WRITE_DIRNAME(piaacmcparams.piaacmcconfdir, "%s", confindex);

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

    {
        uint64_t initflags = INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT;
        initflags |= INIT_PIAACMCOPTICALDESIGN_MODE__READCONF;
        initflags |= INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF;
        if (piaacmcparams.PIAACMC_fpmtype == 1)
        {
            initflags |= INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL;
        }

        FUNC_CHECK_RETURN(
            init_piaacmcopticaldesign(fpmradld, centobs0, centobs1, initflags, NULL));
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * PIAACMCsimul_run_read_searchtime - Read optimization wall-clock limit [s] from searchtime.txt.
 *
 * Return: Search time budget in seconds (defaults to 36000 s when file is absent).
 */
double PIAACMCsimul_run_read_searchtime(void)
{
    double searchtime = 3600.0 * 10.0;
    char   fname[STRINGMAXLEN_FILENAME];

    WRITE_FILENAME(fname, "searchtime.txt");
    FILE *fp = fopen(fname, "r");
    if (fp != NULL)
    {
        int r = fscanf(fp, "%lf\n", &searchtime);
        (void) r;
        fclose(fp);
    }

    printf("searchtime = %f sec\n", searchtime);
    fflush(stdout);

    return searchtime;
}

/**
 * reset_zonez_from_source - Copy best sag solution into zonezID, or zero sags if IDsrc == -1.
 * @IDsrc: Source double array image ID, or -1 to fill with zeros
 */
static void reset_zonez_from_source(imageID IDsrc)
{
    uint32_t nbzones = dcimg[piaacmcopticaldesign.zonezID].md[0].size[0];

    if (IDsrc != -1)
    {
        for (uint32_t k = 0; k < nbzones; k++)
        {
            dcimg[piaacmcopticaldesign.zonezID].array.D[k] = dcimg[IDsrc].array.D[k];
        }
    }
    else
    {
        for (uint32_t k = 0; k < nbzones; k++)
        {
            dcimg[piaacmcopticaldesign.zonezID].array.D[k] = 0.0;
        }
    }
}

/**
 * select_mode13_seed_sags - Choose initial sag vector (previous, zero, or best) for Mode 13 pass.
 * @ctx: Mutable Mode 13 outer loop context
 */
static void select_mode13_seed_sags(struct mode13_loop_ctx *ctx)
{
    if (ctx->loopiter < 1)
    {
        piaacmcparams.MODampl = 0.0;
    }
    else
    {
        piaacmcparams.MODampl = 1.0e-6 * pow(ran1(), 8.0);
    }

    ctx->cnt00  = 0;
    ctx->cnt0   = 0;
    ctx->zeroST = 0;

    if ((ctx->loopiter > 1) && (ran1() > 0.5))
    {
        if ((ran1() > 0.5) && (ctx->IDbestsol != -1))
        {
            ctx->zeroST = 2;
            printf("[%d] Adopting best solution as starting point\n", __LINE__);
            fflush(stdout);
            reset_zonez_from_source(ctx->IDbestsol);
        }
        else
        {
            ctx->zeroST = 1;
            printf("[%d] Adopting zero as starting point\n", __LINE__);
            printf("piaacmcopticaldesign.zonezID = %ld\n", piaacmcopticaldesign.zonezID);
            fflush(stdout);
            reset_zonez_from_source(-1);
        }
    }

    if (ctx->loopiter == 3)
    {
        if (ctx->IDbestsol != -1)
        {
            printf("[%d] Adopting best solution as starting point\n", __LINE__);
            fflush(stdout);
            ctx->zeroST           = 3;
            piaacmcparams.MODampl = 0.0;
            reset_zonez_from_source(ctx->IDbestsol);
        }
        else
        {
            ctx->zeroST = 1;
            printf("[%d] Adopting zero as starting point\n", __LINE__);
            printf("piaacmcopticaldesign.zonezID = %ld\n", piaacmcopticaldesign.zonezID);
            fflush(stdout);
            reset_zonez_from_source(-1);
        }
    }
}

/**
 * PIAACMCsimul_run_prep_mode13_zones - Perturb Mode 13 starting zone sags and dump fpsagtest.txt.
 * @ctx: Mutable Mode 13 outer loop context
 */
void PIAACMCsimul_run_prep_mode13_zones(struct mode13_loop_ctx *ctx)
{
    select_mode13_seed_sags(ctx);

    if (ctx->loopiter <= 0)
    {
        return;
    }

    uint32_t nbzones = dcimg[piaacmcopticaldesign.zonezID].md[0].size[0];
    EXECUTE_SYSTEM_COMMAND("echo \"%g  %ld\" > sag0.txt", ctx->sag0, (long) nbzones);

    double prob1 = pow(ran1(), 8.0);
    for (uint32_t k = 0; k < nbzones; k++)
    {
        if (dcimg[piaacmcopticaldesign.zonezID].array.D[k] > ctx->sag0 / 2.0)
        {
            ctx->cnt00++;
            if (ran1() < prob1)
            {
                dcimg[piaacmcopticaldesign.zonezID].array.D[k] -= ctx->sag0;
                ctx->cnt0++;
            }
        }
        if (dcimg[piaacmcopticaldesign.zonezID].array.D[k] < -ctx->sag0 / 2.0)
        {
            ctx->cnt00++;
            if (ran1() < prob1)
            {
                dcimg[piaacmcopticaldesign.zonezID].array.D[k] += ctx->sag0;
                ctx->cnt0++;
            }
        }
    } // for k

    printf("Write sag values to file fpsagtest.txt\n");
    fflush(stdout);

    FILE *fp = fopen("fpsagtest.txt", "w");
    if (fp != NULL)
    {
        fprintf(fp, "# %9.6f\n", ctx->sag0 * 1.0e6);
        fprintf(fp, "#    %5ld    %5ld    %5ld\n", ctx->cnt0, ctx->cnt00, (long) nbzones);
        for (uint32_t k = 0; k < nbzones; k++)
        {
            fprintf(
                fp,
                "%5ld %9.6f\n",
                (long) k,
                dcimg[piaacmcopticaldesign.zonezID].array.D[k] * 1.0e6);
        }
        fclose(fp);
    }
}
