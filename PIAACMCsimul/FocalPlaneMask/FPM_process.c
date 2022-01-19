/**
 * @file    PIAACMCsimul_FPM_process.c
 * @brief   PIAA-type coronagraph design
 *
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>

// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_memory/COREMOD_memory.h"

// Local variables pointers
static char *FPMsagimname;
static char *zonescoordname;
static long *NBexpval;
static char *FPMsagoutimname;

static CLICMDARGDEF farg[] = {
    {CLIARG_IMG, ".FPMsag", "input FPM sag", "inFPMsag", CLIARG_VISIBLE_DEFAULT, (void **)&FPMsagimname, NULL},
    {CLIARG_STR, ".zcoordname", "sectors ASCII file", "coord.txt", CLIARG_VISIBLE_DEFAULT, (void **)&zonescoordname,
     NULL},
    {CLIARG_LONG, ".NBexp", "number of exposures", "4", CLIARG_VISIBLE_DEFAULT, (void **)&NBexpval, NULL},
    {CLIARG_STR, ".outFPMsag", "output FPM sags", "outFPMsag", CLIARG_VISIBLE_DEFAULT, (void **)&FPMsagoutimname,
     NULL}};

static CLICMDDATA CLIcmddata = {"fpmprocess", "Quantize FPM", CLICMD_FIELDS_DEFAULTS};

// detailed help
static errno_t help_function() { return RETURN_SUCCESS; }

errno_t PIAACMC_FPM_process(const char *__restrict__ FPMsag_name, const char *__restrict__ zonescoord_name, long NBexp,
                            const char *__restrict__ outname)
{
    DEBUG_TRACE_FSTART();

    imageID IDin;
    long NBzones;
    int atype;

    double *sagarray_in;
    double *sagarray_out;
    double sagmax, sagmin;
    long NBsagsteps;

    double *sagstepval;
    FILE *fp;
    FILE *fpout;

    (void)zonescoord_name;

    IDin = image_ID(FPMsag_name);
    NBzones = data.image[IDin].md[0].size[0];
    atype = data.image[IDin].md[0].datatype;

    switch (atype)
    {
    case _DATATYPE_DOUBLE:
        printf("atype = _DATATYPE_DOUBLE\n");
        break;
    case _DATATYPE_FLOAT:
        printf("atype = _DATATYPE_FLOAT\n");
        break;
    default:
        printf("ERROR: atype not supported\n");
        exit(0);
        break;
    }

    printf("%ld zones\n", NBzones);
    sagarray_in = (double *)malloc(sizeof(double) * NBzones);
    sagarray_out = (double *)malloc(sizeof(double) * NBzones);

    for (long zone = 0; zone < NBzones; zone++)
    {
        if (atype == _DATATYPE_FLOAT)
        {
            sagarray_in[zone] = (double)data.image[IDin].array.F[zone];
        }
        else
        {
            sagarray_in[zone] = data.image[IDin].array.D[zone];
        }
    }

    sagmin = sagarray_in[0];
    sagmax = sagarray_in[0];
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
    NBsagsteps = 2;
    for (long k = 1; k < NBexp; k++)
    {
        NBsagsteps *= 2;
    }
    printf("NBsagsteps = %ld\n", NBsagsteps);

    fp = fopen("saglevels.dat", "w");

    sagstepval = (double *)malloc(sizeof(double) * NBsagsteps);
    for (long k = 0; k < NBsagsteps; k++)
    {
        sagstepval[k] = sagmin + (sagmax - sagmin) * k / NBsagsteps + 0.5 * (sagmax - sagmin) / NBsagsteps;
        fprintf(fp, "%4ld     %10.8f\n", k, sagstepval[k] * 1.0e6);
    }
    fclose(fp);
    printf("\n");

    fpout = fopen(outname, "w");

    //for(zone=0; zone<NBzones; zone++)
    for (long zone = 0; zone < NBzones; zone++)
    {
        long k = (long)((sagarray_in[zone] - sagmin) / ((sagmax - sagmin) / NBsagsteps));
        if (sagarray_in[zone] > sagstepval[NBsagsteps - 1])
        {
            k = NBsagsteps - 1;
        }
        //		printf("zone %4ld   %10.8f   %10.8f  %4ld\n", zone, sagarray_in[zone]*1e6, (sagarray_in[zone]-sagmin)*1e6, k);
        fprintf(fpout, "%4ld  %+10.8f  %4ld  %+10.8f\n", zone, sagarray_in[zone] * 1e6, k, sagstepval[k] * 1e6);
    }
    fclose(fpout);

    free(sagstepval);
    free(sagarray_in);
    free(sagarray_out);

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

static errno_t compute_function()
{
    DEBUG_TRACE_FSTART();

    INSERT_STD_PROCINFO_COMPUTEFUNC_START

    PIAACMC_FPM_process(FPMsagimname, zonescoordname, *NBexpval, FPMsagoutimname);

    INSERT_STD_PROCINFO_COMPUTEFUNC_END

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

INSERT_STD_FPSCLIfunctions

    // Register function in CLI
    errno_t
    CLIADDCMD_PIAACMCsimul__FPM_process()
{
    INSERT_STD_CLIREGISTERFUNC
    return RETURN_SUCCESS;
}
