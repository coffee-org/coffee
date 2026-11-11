/**
 * @file    FPMresp_resample.c
 * @brief   Resample focal plane mask response cube in wavelength and evaluation points
 */

#include <stdio.h>
#include <stdlib.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_memory/COREMOD_memory.h"
#include "FPMresp_resample.h"

static char    *FPMrespinimname;
static char    *FPMrespoutimname;
static int64_t *NBlambdaval;
static int64_t *PTstepval;

static CLICMDARGDEF fpmresprs_farg[] =
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
        ".NBlambda",
        "number of output wavelength slices",
        "10",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &NBlambdaval,
        NULL
    },
    {
        CLIARG_INT64,
        ".PTstep",
        "EvalPts step",
        "2",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &PTstepval,
        NULL
    }
};

MILK_RUNNABLE_RT_DEFINE(fpmresprs);

/**
 * PIAACMC_FPMresp_resample - Compress FPM response cube in wavelength and point stride
 * @FPMresp_in_name:  Input 3D FPM response cube name.
 * @FPMresp_out_name: Output resampled 3D FPM response cube name.
 * @NBlambda:         Target number of wavelength planes.
 * @PTstep:           Evaluation-point decimation stride.
 * @outID:            Optional output imageID pointer.
 *
 * Return: RETURN_SUCCESS on completion, error code otherwise.
 */
errno_t PIAACMC_FPMresp_resample(
    const char *restrict FPMresp_in_name,
    const char *restrict FPMresp_out_name,
    long                 NBlambda,
    long                 PTstep,
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
    long     xsize1 = (long) (xsize / PTstep);
    long     zsize1 = NBlambda;

    imageID IDout = -1;
    FUNC_CHECK_RETURN(
        create_3Dimage_ID_double(FPMresp_out_name, xsize1, ysize, zsize1, &IDout));

    for (long kk1 = 0; kk1 < zsize1; kk1++)
    {
        double kk1x  = (zsize1 > 1) ? (1.0 * kk1 / (zsize1 - 1)) : 0.0;
        long   kk1xi = (long) (kk1x * (zsize - 1));
        double alpha = kk1x * (zsize - 1) - kk1xi;

        if (zsize > 1 && kk1xi >= (long) zsize - 1)
        {
            kk1xi = (long) zsize - 2;
            alpha = 1.0;
        }

        long   ii1      = 0;
        size_t plane_sz = (size_t) xsize * ysize;
        size_t oplane   = (size_t) kk1 * xsize1 * ysize;
        for (uint32_t ii = 0; ii < xsize; ii += (uint32_t) (2 * PTstep))
        {
            for (uint32_t jj = 0; jj < ysize; jj++)
            {
                size_t row  = (size_t) jj * xsize + ii;
                size_t idx0 = (size_t) kk1xi * plane_sz + row;
                size_t idx1 = (zsize > 1) ? (idx0 + plane_sz) : idx0;
                double re0  = dcimg[ID].array.D[idx0];
                double im0  = dcimg[ID].array.D[idx0 + 1];
                double re1  = dcimg[ID].array.D[idx1];
                double im1  = dcimg[ID].array.D[idx1 + 1];

                size_t oidx                    = oplane + (size_t) jj * xsize1 + ii1;
                dcimg[IDout].array.D[oidx]     = (1.0 - alpha) * re0 + alpha * re1;
                dcimg[IDout].array.D[oidx + 1] = (1.0 - alpha) * im0 + alpha * im1;
            }
            ii1 += 2;
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
 * cli_fpmresp_resample - CLI wrapper invoking PIAACMC_FPMresp_resample
 *
 * Return: errno_t status code from PIAACMC_FPMresp_resample.
 */
static errno_t cli_fpmresp_resample(void)
{
    return PIAACMC_FPMresp_resample(
        FPMrespinimname,
        FPMrespoutimname,
        (long) *NBlambdaval,
        (long) *PTstepval,
        NULL);
}

const MILK_RUNNABLE runnable_piaacmc_fpmresprs =
{
    .kind     = MILK_RUNNABLE_CLI,
    .app      = { .cmdkey = "fpmresprs", .description = "resample FPM resp matrix" },
    .srcfile  = __FILE__,
    .farg     = fpmresprs_farg,
    .nbarg    = (int) (sizeof(fpmresprs_farg) / sizeof(fpmresprs_farg[0])),
    .cli_func = cli_fpmresp_resample,
    .rt       = &fpmresprs_rt,
};
