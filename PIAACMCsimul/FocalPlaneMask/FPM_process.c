/**
 * @file    FPM_process.c
 * @brief   Quantize focal plane mask zone sags into discrete lithographic steps
 */

#include <stdio.h>
#include <stdlib.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_memory/COREMOD_memory.h"
#include "FPM_process.h"

static char    *FPMsagimname;
static char    *zonescoordname;
static int64_t *NBexpval;
static char    *FPMsagoutimname;

static CLICMDARGDEF fpmprocess_farg[] =
{
    {
        CLIARG_IMG,
        ".FPMsag",
        "input FPM sag",
        "inFPMsag",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &FPMsagimname,
        NULL
    },
    {
        CLIARG_STR,
        ".zcoordname",
        "sectors ASCII file",
        "coord.txt",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &zonescoordname,
        NULL
    },
    {
        CLIARG_INT64,
        ".NBexp",
        "number of exposures",
        "4",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &NBexpval,
        NULL
    },
    {
        CLIARG_STR,
        ".outFPMsag",
        "output FPM sags",
        "outFPMsag",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &FPMsagoutimname,
        NULL
    }
};

MILK_RUNNABLE_RT_DEFINE(fpmprocess);

/**
 * PIAACMC_FPM_process - Quantize FPM zone sags into 2^NBexp lithographic levels
 * @FPMsag_name:     Input 1D FPM sag image name [m].
 * @zonescoord_name: Sector coordinates ASCII file (reserved).
 * @NBexp:           Number of binary etching exposures (2^NBexp levels).
 * @outname:         Output ASCII table path receiving quantized sags [um].
 *
 * Return: RETURN_SUCCESS on completion, RETURN_FAILURE on invalid input.
 */
errno_t PIAACMC_FPM_process(
    const char *restrict FPMsag_name,
    const char *restrict zonescoord_name,
    long                 NBexp,
    const char *restrict outname)
{
    DEBUG_TRACE_FSTART();
    (void) zonescoord_name;

    imageID IDin = image_ID(FPMsag_name, dcimg, dcnimg);
    if (IDin == -1)
    {
        PRINT_ERROR("Image \"%s\" not found", FPMsag_name);
        return RETURN_FAILURE;
    }

    long NBzones = (long) dcimg[IDin].md[0].size[0];
    int  atype   = dcimg[IDin].md[0].datatype;

    if (atype != _DATATYPE_DOUBLE && atype != _DATATYPE_FLOAT)
    {
        PRINT_ERROR("Unsupported datatype %d for \"%s\"", atype, FPMsag_name);
        return RETURN_FAILURE;
    }

    double *sagarray_in = (double *) malloc(sizeof(double) * (size_t) NBzones);
    if (sagarray_in == NULL)
    {
        return RETURN_FAILURE;
    }

    for (long zone = 0; zone < NBzones; zone++)
    {
        if (atype == _DATATYPE_FLOAT)
        {
            sagarray_in[zone] = (double) dcimg[IDin].array.F[zone];
        }
        else
        {
            sagarray_in[zone] = dcimg[IDin].array.D[zone];
        }
    }

    double sagmin = sagarray_in[0];
    double sagmax = sagarray_in[0];
    for (long zone = 1; zone < NBzones; zone++)
    {
        if (sagarray_in[zone] < sagmin)
        {
            sagmin = sagarray_in[zone];
        }
        if (sagarray_in[zone] > sagmax)
        {
            sagmax = sagarray_in[zone];
        }
    }

    printf("Sag range [um]  :   %10.8f  ->  %10.8f\n", sagmin * 1.0e6, sagmax * 1.0e6);
    long NBsagsteps = 2;
    for (long k = 1; k < NBexp; k++)
    {
        NBsagsteps *= 2;
    }
    printf("NBsagsteps = %ld\n", NBsagsteps);

    double *sagstepval = (double *) malloc(sizeof(double) * (size_t) NBsagsteps);
    FILE   *fp         = fopen("saglevels.dat", "w");
    if (sagstepval == NULL || fp == NULL)
    {
        if (fp != NULL)
        {
            fclose(fp);
        }
        free(sagstepval);
        free(sagarray_in);
        return RETURN_FAILURE;
    }

    for (long k = 0; k < NBsagsteps; k++)
    {
        sagstepval[k] = sagmin + (sagmax - sagmin) * k / NBsagsteps
                        + 0.5 * (sagmax - sagmin) / NBsagsteps;
        fprintf(fp, "%4ld     %10.8f\n", k, sagstepval[k] * 1.0e6);
    }
    fclose(fp);

    FILE *fpout = fopen(outname, "w");
    if (fpout == NULL)
    {
        free(sagstepval);
        free(sagarray_in);
        return RETURN_FAILURE;
    }

    for (long zone = 0; zone < NBzones; zone++)
    {
        long k = 0;
        if (sagmax > sagmin)
        {
            k = (long) ((sagarray_in[zone] - sagmin) / ((sagmax - sagmin) / NBsagsteps));
        }
        if (k < 0)
        {
            k = 0;
        }
        if (k >= NBsagsteps || sagarray_in[zone] > sagstepval[NBsagsteps - 1])
        {
            k = NBsagsteps - 1;
        }
        fprintf(
            fpout,
            "%4ld  %+10.8f  %4ld  %+10.8f\n",
            zone,
            sagarray_in[zone] * 1.0e6,
            k,
            sagstepval[k] * 1.0e6);
    }
    fclose(fpout);

    free(sagstepval);
    free(sagarray_in);

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * cli_fpm_process - CLI wrapper invoking PIAACMC_FPM_process
 *
 * Return: errno_t status code from PIAACMC_FPM_process.
 */
static errno_t cli_fpm_process(void)
{
    return PIAACMC_FPM_process(
        FPMsagimname,
        zonescoordname,
        (long) *NBexpval,
        FPMsagoutimname);
}

const MILK_RUNNABLE runnable_piaacmc_fpmprocess =
{
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "fpmprocess", .description = "Quantize FPM" },
    .srcfile  = __FILE__,
    .farg     = fpmprocess_farg,
    .nbarg    = (int) (sizeof(fpmprocess_farg) / sizeof(fpmprocess_farg[0])),
    .cli_func = cli_fpm_process,
    .rt       = &fpmprocess_rt,
};
