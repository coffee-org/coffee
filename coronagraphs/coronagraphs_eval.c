/**
 * @file    coronagraphs_eval.c
 * @brief   Top-level coronagraph PSF dispatcher and user evaluation entry point
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "coronagraphs_internal.h"

typedef struct eval_aplc_preset
{
    long   coro_type;
    int    set_piaa;
    int    aplc_piaa;
    long   nb_step;
    double fpm_size;
    int    set_centobs;
    double centobs;
    int    log_size_err;
} eval_aplc_preset_t;

/**
 * eval_parse_fpmserr - Parse optional -FPMSERR flag from simulation options string
 * @options: Command-line style options string
 */
static void eval_parse_fpmserr(const char *options)
{
    const char *pos = strstr(options, "-FPMSERR ");
    if (pos == NULL)
    {
        return;
    }

    pos += strlen("-FPMSERR ");
    char input[200];
    int  ii = 0;
    while ((pos[ii] != ' ') && (pos[ii] != '\n') && (pos[ii] != '\0'))
    {
        input[ii] = pos[ii];
        ii++;
    }
    input[ii]                = '\0';
    g_coro.FPMASKSIZE_ERROR  = 1;
    g_coro.FPMASK_size_error = atof(input);
    printf("-------------- MASK size error = %f\n", g_coro.FPMASK_size_error);
}

/**
 * eval_try_aplc_preset - Match and run a preset multi-step APLC/PIAACMC configuration
 * @coro_type: Coronagraph architecture code
 * @xld:       Source X offset [lambda/D]
 * @yld:       Source Y offset [lambda/D]
 * @psfname:   Output PSF image name
 *
 * Return: 1 if matched and executed, 0 otherwise.
 */
static int eval_try_aplc_preset(
    long        coro_type,
    double      xld,
    double      yld,
    const char *psfname)
{
    static const eval_aplc_preset_t presets[] = {
        {10, 0, 0, 1, 1.8, 0, 0.0, 0},   {11, 0, 0, 2, 1.4, 0, 0.0, 0},
        {12, 0, 0, 3, 1.2, 0, 0.0, 0},   {13, 0, 0, 4, 1.0, 0, 0.0, 0},
        {14, 0, 0, 5, 1.0, 0, 0.0, 0},   {15, 0, 0, 0, 4.2, 0, 0.0, 0},
        {31, 1, 1, 2, 1.6, 0, 0.0, 0},   {32, 1, 1, 3, 1.3, 0, 0.0, 0},
        {33, 1, 1, 4, 1.2, 0, 0.0, 0},   {34, 1, 1, 2, 2.0, 0, 0.0, 0},
        {35, 1, 1, 1, 2.4, 0, 0.0, 0},   {36, 1, 0, 0, 2.6, 0, 0.0, 0},
        {37, 1, 0, 1, 1.4, 0, 0.0, 0},   {38, 1, 0, 2, 1.0, 0, 0.0, 0},
        {39, 1, 1, 0, 2.6, 0, 0.0, 0},   {40, 1, 1, 1, 1.4, 0, 0.0, 0},
        {41, 1, 1, 2, 1.0, 0, 0.0, 0},   {42, 1, 1, 0, 2.6, 1, 0.275, 0},
        {43, 1, 1, 1, 1.4, 1, 0.275, 0}, {44, 1, 1, 2, 1.0, 1, 0.275, 0},
        {45, 1, 2, 1, 2.0, 0, 0.0, 1},
    };
    const size_t n_presets = sizeof(presets) / sizeof(presets[0]);

    for (size_t idx = 0; idx < n_presets; idx++)
    {
        if (presets[idx].coro_type != coro_type)
        {
            continue;
        }
        if (presets[idx].set_centobs)
        {
            g_coro.PIAACENTOBS = presets[idx].centobs;
        }
        if (presets[idx].set_piaa)
        {
            g_coro.APLC_PIAA = presets[idx].aplc_piaa;
        }
        g_coro.NB_APLC_STEP    = presets[idx].nb_step;
        g_coro.APLC_FPMASKsize = presets[idx].fpm_size;
        if (presets[idx].log_size_err)
        {
            fprintf(stderr, "MASK SIZE ERROR = %f\n", g_coro.FPMASK_size_error);
        }
        coronagraph_simul_MULTISTEP_APLC(xld, yld, psfname);
        return 1;
    }

    return 0;
}

/**
 * eval_try_param_coro - Dispatch parameterized PIAAC, BL8, SHEAR4, STRIPC, and OVC modes
 * @coro_type: Coronagraph architecture code
 * @xld:       Source X offset [lambda/D]
 * @yld:       Source Y offset [lambda/D]
 * @psfname:   Output PSF image name
 *
 * Return: 1 if matched and executed, 0 otherwise.
 */
static int eval_try_param_coro(
    long        coro_type,
    double      xld,
    double      yld,
    const char *psfname)
{
    if ((coro_type >= 25) && (coro_type <= 29))
    {
        static const long ovc_charges[5] = {2, 4, 5, 6, 8};
        g_coro.OVC_CHARGE                = ovc_charges[coro_type - 25];
        coronagraph_simul_OVC(xld, yld, psfname);
        return 1;
    }

    if ((coro_type == 5) || (coro_type == 21) || (coro_type == 24))
    {
        g_coro.PIAAFPMASKRAD = (coro_type == 5) ? 5.5 : ((coro_type == 21) ? 5.0 : 3.5);
        coronagraph_simul_PIAAC(xld, yld, psfname);
        return 1;
    }

    if ((coro_type == 6) || (coro_type == 7))
    {
        g_coro.BL8MODE = (coro_type == 7) ? 1 : 0;
        coronagraph_simul_BL8(xld, yld, psfname);
        return 1;
    }

    if ((coro_type == 19) || (coro_type == 22) || (coro_type == 23))
    {
        g_coro.STRIPCOFFSET = (coro_type == 19) ? 0.4 : ((coro_type == 22) ? 0.35 : 0.45);
        coronagraph_simul_STRIPC(xld, yld, psfname);
        return 1;
    }

    if (coro_type == 17)
    {
        g_coro.SHEAR4_SHEAR = 0.3;
        coronagraph_simul_SHEAR4(xld, yld, psfname);
        return 1;
    }

    return 0;
}

/**
 * eval_try_direct_coro - Dispatch non-parameterized coronagraph simulator functions
 * @coro_type: Coronagraph architecture code
 * @xld:       Source X offset [lambda/D]
 * @yld:       Source Y offset [lambda/D]
 * @psfname:   Output PSF image name
 *
 * Return: 1 if matched and executed, 0 otherwise.
 */
static int eval_try_direct_coro(
    long        coro_type,
    double      xld,
    double      yld,
    const char *psfname)
{
    typedef int (*coro_sim_fn_t)(double, double, const char *);
    typedef struct eval_direct_entry
    {
        long          coro_type;
        coro_sim_fn_t fn;
    } eval_direct_entry_t;

    static const eval_direct_entry_t entries[] = {
        {-1, coronagraph_simul_NOCORO},
        {0, coronagraph_simul_DICC},
        {1, coronagraph_simul_AIC},
        {2, coronagraph_simul_RRPM},
        {3, coronagraph_simul_4QPM},
        {4, coronagraph_simul_PIAA},
        {8, coronagraph_simul_AIC_PIAAC},
        {9, coronagraph_simul_CPA},
        {16, coronagraph_simul_ODC},
        {18, coronagraph_simul_BL4},
        {20, coronagraph_simul_SIMXY},
        {30, coronagraph_simul_PPA},
        {60, coronagraph_simul_EXTERNAL_OCCULTER},
    };
    const size_t n_entries = sizeof(entries) / sizeof(entries[0]);

    for (size_t idx = 0; idx < n_entries; idx++)
    {
        if (entries[idx].coro_type == coro_type)
        {
            entries[idx].fn(xld, yld, psfname);
            return 1;
        }
    }

    return 0;
}

/**
 * eval_run_aplc_conf - Load APLCPIAA.conf and optional MASKerr.conf, then simulate PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output PSF image name
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE.
 */
static errno_t eval_run_aplc_conf(
    double      xld,
    double      yld,
    const char *psfname)
{
    g_coro.APLC_PIAA     = 1;
    g_coro.APLC_PMASK    = 1;
    g_coro.APLC_CentOBS0 = 0.0;
    g_coro.APLC_CentOBS1 = 0.0;

    FILE *fp = fopen("APLCPIAA.conf", "r");
    if (fp == NULL)
    {
        PRINT_ERROR("Cannot read configuration file \"APLCPIAA.conf\"\n");
        exit(0);
    }

    char str1[200];
    int  fscanfcnt = fscanf(
        fp, "%ld %199s %d %d\n", &g_coro.NB_APLC_STEP, str1, &g_coro.APLC_PIAA, &g_coro.APLC_PMASK);
    if (coro_check_fscanf(fp, fscanfcnt, 4) != RETURN_SUCCESS)
    {
        fclose(fp);
        return RETURN_FAILURE;
    }
    fclose(fp);

    int has_obs = 0;
    for (size_t ii = 0; str1[ii] != '\0'; ii++)
    {
        if (str1[ii] == '-')
        {
            has_obs  = 1;
            str1[ii] = ' ';
        }
    }
    if (has_obs == 1)
    {
        double tmpf0 = 0.0;
        double tmpf1 = 0.0;
        sscanf(str1, "%lf %lf %lf", &g_coro.APLC_FPMASKsize, &tmpf0, &tmpf1);
        g_coro.APLC_CentOBS0 = tmpf0;
        g_coro.APLC_CentOBS1 = tmpf1;
    }
    else
    {
        sscanf(str1, "%lf", &g_coro.APLC_FPMASKsize);
    }

    fprintf(stderr, "NB_APLC_STEP = %ld\n", g_coro.NB_APLC_STEP);
    fprintf(stderr, "APLC_FPMASKsize = %f\n", g_coro.APLC_FPMASKsize);
    fprintf(stderr, "APLC_CentOBS0 = %f\n", g_coro.APLC_CentOBS0);
    fprintf(stderr, "APLC_CentOBS1 = %f\n", g_coro.APLC_CentOBS1);
    fprintf(stderr, "APLC_PIAA = %d\n", g_coro.APLC_PIAA);
    fprintf(stderr, "APLC_PMASK = %d\n", g_coro.APLC_PMASK);

    fp = fopen("MASKerr.conf", "r");
    if (fp != NULL)
    {
        fscanfcnt = fscanf(fp, "%lf %lf\n", &g_coro.FPMASK_transm_error, &g_coro.FPMASK_size_error);
        if (coro_check_fscanf(fp, fscanfcnt, 2) != RETURN_SUCCESS)
        {
            fclose(fp);
            return RETURN_FAILURE;
        }
        fclose(fp);
        fprintf(stderr, "MASK TRANSM ERROR = %f\n", g_coro.FPMASK_transm_error);
        fprintf(stderr, "MASK SIZE ERROR = %f\n", g_coro.FPMASK_size_error);
    }
    else
    {
        printf("No focal plane mask error file - using perfect focal plane mask\n");
        g_coro.FPMASK_transm_error = 0.0;
        g_coro.FPMASK_size_error   = 0.0;
    }

    coronagraph_simul_MULTISTEP_APLC(xld, yld, psfname);
    return RETURN_SUCCESS;
}

/**
 * coronagraph_simulPSF - Dispatch coronagraph PSF simulation by coronagraph_type code
 * @xld:              Source X offset [lambda/D]
 * @yld:              Source Y offset [lambda/D]
 * @psfname:          Output PSF image name
 * @coronagraph_type: Architecture selector index
 * @options:          Option flags string (e.g. "-FPMSERR <val>")
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE.
 */
errno_t coronagraph_simulPSF(
    double      xld,
    double      yld,
    const char *psfname,
    long        coronagraph_type,
    const char *options)
{
    printf("SCALE : %f l/D per pix\n", CORONAGRAPHS_PIXSCALE);
    printf("CORONAGRAPHSDATADIR = %s\n", CORONAGRAPHSDATADIR);
    printf("CORONAGRAPHSDATALOCAL = %s\n", CORONAGRAPHSDATALOCAL);
    printf("CORONAGRAPHS_ARRAYSIZE = %d\n", CORONAGRAPHS_ARRAYSIZE);
    printf("OPTIONS = \"%s\"\n", options);

    eval_parse_fpmserr(options);

    if (coronagraph_type == 50)
    {
        return eval_run_aplc_conf(xld, yld, psfname);
    }

    if (eval_try_aplc_preset(coronagraph_type, xld, yld, psfname))
    {
        return RETURN_SUCCESS;
    }
    if (eval_try_param_coro(coronagraph_type, xld, yld, psfname))
    {
        return RETURN_SUCCESS;
    }
    eval_try_direct_coro(coronagraph_type, xld, yld, psfname);

    return RETURN_SUCCESS;
}

/**
 * coronagraph_userfunc - Execute user-configured PIAA Zernike perturbation evaluation hook
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE.
 */
errno_t coronagraph_userfunc(void)
{
    long   nb_zern = 1;
    double ratio   = 10.0;

    fprintf(stderr, "CORONAGRAPH USER FUNCTION\n");

    double *zarray = (double *) malloc(sizeof(double) * nb_zern);
    if (zarray == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        return RETURN_FAILURE;
    }

    long *zindex = (long *) malloc(sizeof(long) * nb_zern);
    if (zindex == NULL)
    {
        free(zarray);
        PRINT_ERROR("malloc returns NULL pointer");
        return RETURN_FAILURE;
    }

    zindex[0] = 12;
    zarray[0] = 1.0e-5;
    coronagraph_PIAAperturbation(zarray, zindex, nb_zern, ratio);
    free(zarray);
    free(zindex);
    exit(0);

    return RETURN_SUCCESS;
}
