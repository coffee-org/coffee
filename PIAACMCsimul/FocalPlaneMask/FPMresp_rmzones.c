/**
 * @file    FPMresp_rmzones.c
 * @brief   Remove outer zones from focal plane mask response matrix
 */

#include <stdio.h>
#include <stdlib.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_memory/COREMOD_memory.h"
#include "FPMresp_rmzones.h"

static char    *FPMrespinimname;
static char    *FPMrespoutimname;
static int64_t *NBzrmval;

static CLICMDARGDEF fpmresprm_farg[] =
{
    {
        CLIARG_IMG,
        ".FPMrespin",
        "input FPM response image",
        "FPMresp",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &FPMrespinimname,
        NULL
    },
    {
        CLIARG_STR,
        ".FPMrespout",
        "output FPM response image",
        "FPMrespout",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &FPMrespoutimname,
        NULL
    },
    {
        CLIARG_INT64,
        ".NBzrm",
        "NBzone removed",
        "125",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &NBzrmval,
        NULL
    }
};

MILK_RUNNABLE_RT_DEFINE(fpmresprm);

/**
 * FPMresp_rmzones - Fold outer NBzones of an FPM response cube into zone 0
 * @FPMresp_in_name:  Input 3D FPM response cube name.
 * @FPMresp_out_name: Output reduced 3D FPM response cube name.
 * @NBzones:          Number of outer zones to fold into the reference slice.
 * @outID:            Optional output imageID pointer.
 *
 * Return: RETURN_SUCCESS on completion, error code otherwise.
 */
errno_t FPMresp_rmzones(
    const char *restrict FPMresp_in_name,
    const char *restrict FPMresp_out_name,
    long                 NBzones,
    imageID             *outID)
{
    DEBUG_TRACE_FSTART();

    imageID ID = image_ID(FPMresp_in_name, dcimg, dcnimg);
    if (ID == -1)
    {
        PRINT_ERROR("Image \"%s\" not found", FPMresp_in_name);
        return RETURN_FAILURE;
    }

    uint32_t xsize  = dcimg[ID].md[0].size[0];
    uint32_t ysize  = dcimg[ID].md[0].size[1];
    uint32_t zsize  = dcimg[ID].md[0].size[2];
    long     ysize1 = (long) ysize - NBzones;
    if (ysize1 < 1)
    {
        PRINT_ERROR("NBzones (%ld) >= ysize (%u)", NBzones, ysize);
        return RETURN_FAILURE;
    }

    imageID IDout = -1;
    FUNC_CHECK_RETURN(
        create_3Dimage_ID_double(FPMresp_out_name, xsize, ysize1, zsize, &IDout));

    for (uint32_t ii = 0; ii < xsize; ii++)
    {
        for (uint32_t kk = 0; kk < zsize; kk++)
        {
            for (long jj = 0; jj < ysize1; jj++)
            {
                dcimg[IDout].array.D[(size_t) kk * xsize * ysize1 + (size_t) jj * xsize + ii] =
                    dcimg[ID].array.D[(size_t) kk * xsize * ysize + (size_t) jj * xsize + ii];
            }
            for (long jj = ysize1; jj < (long) ysize; jj++)
            {
                dcimg[IDout].array.D[(size_t) kk * xsize * ysize1 + ii] +=
                    dcimg[ID].array.D[(size_t) kk * xsize * ysize + (size_t) jj * xsize + ii];
            }
        }
    }

    if (outID != NULL)
    {
        *outID = IDout;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * cli_fpmresp_rmzones - CLI wrapper invoking FPMresp_rmzones
 *
 * Return: errno_t status code from FPMresp_rmzones.
 */
static errno_t cli_fpmresp_rmzones(void)
{
    return FPMresp_rmzones(FPMrespinimname, FPMrespoutimname, (long) *NBzrmval, NULL);
}

const MILK_RUNNABLE runnable_piaacmc_fpmresprm =
{
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "fpmresprm", .description = "remove outer zones in FPM resp matrix" },
    .srcfile  = __FILE__,
    .farg     = fpmresprm_farg,
    .nbarg    = (int) (sizeof(fpmresprm_farg) / sizeof(fpmresprm_farg[0])),
    .cli_func = cli_fpmresp_rmzones,
    .rt       = &fpmresprm_rt,
};
