/**
 * @file    rings2sectors.c
 * @brief   Expand circular ring FPM sags onto azimuthal sector zones
 */

#include <stdio.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_memory/COREMOD_memory.h"
#include "rings2sectors.h"

static char *inimname;
static char *secfname;
static char *outimname;

static CLICMDARGDEF ring2sect_farg[] =
{
    {
        CLIARG_IMG,
        ".inimname",
        "input image: circular mask design",
        "imin",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &inimname,
        NULL
    },
    {
        CLIARG_STR,
        ".secfname",
        "text file specifying which zones belong to which rings",
        "sec.txt",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &secfname,
        NULL
    },
    {
        CLIARG_STR_NOT_IMG,
        ".outimname",
        "output sector mask design",
        "outim",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &outimname,
        NULL
    }
};

MILK_RUNNABLE_RT_DEFINE(ring2sect);

/**
 * rings2sectors - Map 1D ring sag solution onto 1D sector zone sag array
 * @IDin_name:  Input circular ring mask sag image name.
 * @sectfname:  ASCII mapping file listing (zone_index, ring_index) pairs.
 * @IDout_name: Output sector mask sag image name.
 * @outID:      Optional output imageID pointer.
 *
 * Return: RETURN_SUCCESS on completion, error code otherwise.
 */
errno_t rings2sectors(
    const char *IDin_name,
    const char *sectfname,
    const char *IDout_name,
    imageID    *outID)
{
    DEBUG_TRACE_FSTART();

    imageID IDin = image_ID(IDin_name, dcimg, dcnimg);
    if (IDin == -1)
    {
        PRINT_ERROR("Image \"%s\" not found", IDin_name);
        return RETURN_FAILURE;
    }

    FILE *fp = fopen(sectfname, "r");
    if (fp == NULL)
    {
        PRINT_ERROR("Cannot open sector file \"%s\"", sectfname);
        return RETURN_FAILURE;
    }

    long arrayring[5000] = { 0 };
    long nbzone          = 0;
    long nbring          = 0;
    long tmpl1           = 0;
    long tmpl2           = 0;

    while (fscanf(fp, "%ld %ld\n", &tmpl1, &tmpl2) == 2)
    {
        if (tmpl1 >= 0 && tmpl1 < 5000)
        {
            arrayring[tmpl1] = tmpl2;
        }
        if (tmpl2 > nbring)
        {
            nbring = tmpl2;
        }
        if (tmpl1 > nbzone)
        {
            nbzone = tmpl1;
        }
    }
    fclose(fp);
    nbring++;
    nbzone++;

    imageID IDout = -1;
    FUNC_CHECK_RETURN(create_2Dimage_ID_double(IDout_name, nbzone, 1, &IDout));

    for (long zone = 0; zone < nbzone; zone++)
    {
        dcimg[IDout].array.D[zone] = dcimg[IDin].array.D[arrayring[zone]];
    }

    printf("%ld zones in %ld rings\n", nbzone, nbring);

    if (outID != NULL)
    {
        *outID = IDout;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * cli_rings2sectors - CLI wrapper invoking rings2sectors
 *
 * Return: errno_t status code from rings2sectors.
 */
static errno_t cli_rings2sectors(void)
{
    return rings2sectors(inimname, secfname, outimname, NULL);
}

const MILK_RUNNABLE runnable_piaacmc_ring2sect =
{
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "ring2sect", .description = "turn ring fpm design into sectors" },
    .srcfile  = __FILE__,
    .farg     = ring2sect_farg,
    .nbarg    = (int) (sizeof(ring2sect_farg) / sizeof(ring2sect_farg[0])),
    .cli_func = cli_rings2sectors,
    .rt       = &ring2sect_rt,
};
