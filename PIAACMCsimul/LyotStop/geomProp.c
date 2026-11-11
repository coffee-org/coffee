/**
 * @file    geomProp.c
 * @brief   Beam geometrical propagation based on local slopes
 *
 * Uses local slopes of sag map to propagate (un-corrected intensity only) a beam
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// milk includes
#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "image_filter/image_filter.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul.h"
#include "geomProp.h"

// Local variables pointers
static char   *inintensity_imname;
static char   *sag2D_imname;
static char   *outintensity_imname;
static char   *outcnt_imname;
static double *drindexval;
static double *pscaleval;
static double *zpropval;
static double *kradval;
static double *kstepval;
static double *rlimval;

static CLICMDARGDEF geomprop_farg[] =
{
    {
        CLIARG_IMG,
        ".inintim",
        "input intensity image",
        "pupin",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &inintensity_imname,
        NULL
    },
    {
        CLIARG_IMG,
        ".insagim",
        "input 2D sag image",
        "piaa0z",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &sag2D_imname,
        NULL
    },
    {
        CLIARG_STR,
        ".outintim",
        "output intensity image",
        "pupout",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &outintensity_imname,
        NULL
    },
    {
        CLIARG_STR,
        ".outcntim",
        "output intensity image",
        "cntout",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &outcnt_imname,
        NULL
    },
    {
        CLIARG_FLOAT64,
        ".drindex",
        "delta refractive index",
        "2.0",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &drindexval,
        NULL
    },
    {
        CLIARG_FLOAT64,
        ".pscale",
        "pixel scale",
        "0.00011",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &pscaleval,
        NULL
    },
    {
        CLIARG_FLOAT64,
        ".zprop",
        "propagation dist",
        "2.302606",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &zpropval,
        NULL
    },
    {
        CLIARG_FLOAT64,
        ".krad",
        "kernel radius",
        "3.0",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &kradval,
        NULL
    },
    {
        CLIARG_FLOAT64,
        ".kstep",
        "kernel step",
        "0.5",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &kstepval,
        NULL
    },
    {
        CLIARG_FLOAT64,
        ".rlim",
        "clear aperture",
        "200.0",
        (FPFLAG_DEFAULT_INPUT | FPFLAG_PRIMARY_CLI_INPUT),
        (void **) &rlimval,
        NULL
    }
};

MILK_RUNNABLE_RT_DEFINE(geomprop);

// detailed help

/**
 * PIAACMCsimul_geomProp - Propagate 2D intensity geometrically from surface sag slopes
 * @IDin_name:     Input 2D intensity image name.
 * @IDsag_name:    Input 2D surface sag image name [m].
 * @IDout_name:    Output 2D propagated intensity image name.
 * @IDoutcnt_name: Output 2D ray-count image name.
 * @drindex:       Effective refractive index delta (2.0 for reflective mirror).
 * @pscale:        Pupil pixel scale [m/pix].
 * @zprop:         Propagation distance [m].
 * @krad:          Cosine-tapered slope estimation kernel radius [pix].
 * @kstep:         Ray launch grid step in input pupil [pix].
 * @rlim:          Clear aperture radius around array center [pix].
 * @outID:         Optional output imageID pointer.
 *
 * Return: RETURN_SUCCESS on completion, error code otherwise.
 */
errno_t PIAACMCsimul_geomProp(
    const char *restrict IDin_name,
    const char *restrict IDsag_name,
    const char *restrict IDout_name,
    const char *restrict IDoutcnt_name,
    double               drindex,
    double               pscale,
    double               zprop,
    double               krad,
    double               kstep,
    double               rlim,
    imageID             *outID)
{
    DEBUG_TRACE_FSTART();

    imageID  IDin  = image_ID(IDin_name, dcimg, dcnimg);
    uint32_t xsize = dcimg[IDin].md[0].size[0];
    uint32_t ysize = dcimg[IDin].md[0].size[1];
    imageID  IDsag = image_ID(IDsag_name, dcimg, dcnimg);

    imageID IDout;
    FUNC_CHECK_RETURN(create_2Dimage_ID(IDout_name, xsize, ysize, &IDout));

    imageID IDoutcnt;
    FUNC_CHECK_RETURN(create_2Dimage_ID(IDoutcnt_name, xsize, ysize, &IDoutcnt));

    printf("kstep = %f\n", kstep);

    for (float x = 0.5f * (float) xsize - (float) rlim;
            x < 0.5f * (float) xsize + (float) rlim;
            x += (float) kstep)
    {
        for (float y = 0.5f * (float) ysize - (float) rlim;
                y < 0.5f * (float) ysize + (float) rlim;
                y += (float) kstep)
        {
            long   ii0    = (long) (x + 0.5f);
            long   jj0    = (long) (y + 0.5f);
            double sumvx  = 0.0;
            double sumvy  = 0.0;
            double sumvx0 = 0.0;
            double sumvy0 = 0.0;
            double sumv   = 0.0;
            double sumc   = 0.0;

            long iimin = ii0 - ((long) krad + 1);
            long iimax = ii0 + ((long) krad + 2);
            long jjmin = jj0 - ((long) krad + 1);
            long jjmax = jj0 + ((long) krad + 2);

            if (iimin < 0 || iimax > (long) xsize - 1 || jjmin < 0 || jjmax > (long) ysize - 1)
            {
                PRINT_ERROR("Ray kernel out of bounds at (%ld, %ld)", ii0, jj0);
                abort();
            }

            for (long ii = iimin; ii < iimax; ii++)
            {
                for (long jj = jjmin; jj < jjmax; jj++)
                {
                    float dx = (float) ii - x;
                    float dy = (float) jj - y;
                    float dr = sqrtf(dx * dx + dy * dy) / (float) krad;
                    if (dr > 1.0f)
                    {
                        dr = 1.0f;
                    }
                    float coeff = 0.5f + 0.5f * cosf(dr * (float) M_PI);
                    float sagv  = dcimg[IDsag].array.F[jj * xsize + ii];

                    sumv += coeff * sagv;
                    sumvx += coeff * dx * sagv;
                    sumvy += coeff * dy * sagv;
                    sumvx0 += coeff * dx;
                    sumvy0 += coeff * dy;
                    sumc += coeff;
                }
            }

            float slx = (float) (((sumvx / sumc) - (sumvx0 * (sumv / sumc) / sumc)) / pscale);
            float sly = (float) (((sumvy / sumc) - (sumvy0 * (sumv / sumc) / sumc)) / pscale);
            float dii = (float) ((slx * zprop * drindex) / pscale);
            float djj = (float) ((sly * zprop * drindex) / pscale);

            long ii1 = (long) (x + dii + 0.5f);
            long jj1 = (long) (y + djj + 0.5f);

            if ((ii1 > 0) && (ii1 < (long) xsize) && (jj1 > 0) && (jj1 < (long) ysize))
            {
                dcimg[IDout].array.F[jj1 * xsize + ii1] +=
                    dcimg[IDin].array.F[jj0 * xsize + ii0];
                dcimg[IDoutcnt].array.F[jj1 * xsize + ii1] += 1.0f;
            }
        }
    }

    for (uint32_t ii1 = 0; ii1 < xsize; ii1++)
    {
        for (uint32_t jj1 = 0; jj1 < ysize; jj1++)
        {
            if (dcimg[IDoutcnt].array.F[jj1 * xsize + ii1] > 0.1f)
            {
                dcimg[IDout].array.F[jj1 * xsize + ii1] /=
                    dcimg[IDoutcnt].array.F[jj1 * xsize + ii1];
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
 * cli_geomprop - CLI wrapper invoking PIAACMCsimul_geomProp
 *
 * Return: errno_t status code from PIAACMCsimul_geomProp.
 */
static errno_t cli_geomprop(void)
{
    return PIAACMCsimul_geomProp(
        inintensity_imname,
        sag2D_imname,
        outintensity_imname,
        outcnt_imname,
        *drindexval,
        *pscaleval,
        *zpropval,
        *kradval,
        *kstepval,
        *rlimval,
        NULL);
}

const MILK_RUNNABLE runnable_piaacmc_geomprop =
{
    .kind     = MILK_RUNNABLE_CLI,
    .app      = {
        .cmdkey      = "piaacmcgeomprop",
        .description = "Geometric propagation from surface"
    },
    .srcfile  = __FILE__,
    .farg     = geomprop_farg,
    .nbarg    = (int) (sizeof(geomprop_farg) / sizeof(geomprop_farg[0])),
    .cli_func = cli_geomprop,
    .rt       = &geomprop_rt,
};
