/**
 * @file    init_geomPIAA_rad_march.c
 * @brief   Flux-conserving radial ray-map inversion and ODE marching for 1D PIAA mirror sags
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "PIAACMCsimul.h"

#include "init_geomPIAA_rad_internal.h"

/**
 * invert_cumulative_flux_map - Map source radial grid to destination radii of equal flux.
 * @src_flux: Cumulative enclosed flux array on the source pupil (`NBradpts` elements)
 * @dst_flux: Cumulative enclosed flux array on the destination pupil (`NBradpts` elements)
 * @NBradpts: Number of 1D radial sample points
 * @rlim:     Outer dimensionless radius bound on the destination pupil
 * @out_r:    Output mapped dimensionless destination radii (`NBradpts` elements)
 */
static void invert_cumulative_flux_map(
    const double *restrict src_flux,
    const double *restrict dst_flux,
    long                   NBradpts,
    double                 rlim,
    double       *restrict out_r)
{
    out_r[0] = 0.0;
    long ii  = 0;

    for (long i = 1; i < NBradpts; i++)
    {
        double F0 = src_flux[i];

        while ((ii < NBradpts) && (dst_flux[ii] < F0))
        {
            ii++;
        }

        double F1 = dst_flux[ii - 1];
        double F2 = dst_flux[ii];

        /* F0 = F1 + ( (F2-F1)/(ii^2-(ii-1)^2) * ((ii-1+x)^2-(ii-1)^2) ) */
        double fluxdens;
        if (fabs(F2 - F1) > 0.0000000001)
        {
            fluxdens = (F2 - F1) / (2.0 * ii - 1.0);
        }
        else
        {
            fluxdens = 0.0000000001;
        }

        double x = sqrt((F0 - F1) / fluxdens + (1.0 * ii * ii - 2.0 * ii + 1.0)) + 1.0 - 1.0 * ii;
        out_r[i] = rlim * (1.0 * ii - 1.0 + x) / NBradpts;
    }
}

/**
 * march_piaa_mirror_sags - Integrate local surface slopes to obtain axial sags M0z and M1z.
 * @piaar00:  Uniformly spaced dimensionless radii on PIAA element 0
 * @piaar10:  Mapped dimensionless radii on PIAA element 1 for each element-0 ray
 * @NBradpts: Number of 1D radial sample points
 * @piaaM0z:  Output axial sag profile [m] for PIAA element 0
 * @piaaM1z:  Output axial sag profile [m] for PIAA element 1
 */
static void march_piaa_mirror_sags(
    const double *restrict piaar00,
    const double *restrict piaar10,
    long                   NBradpts,
    double       *restrict piaaM0z,
    double       *restrict piaaM1z)
{
    piaaM0z[0] = 0.0;
    piaaM1z[0] = piaacmcopticaldesign.PIAAsep;

    for (long i = 0; i < NBradpts - 1; i++)
    {
        double r0c  = piaar00[i];
        double r1c  = piaar10[i];
        double dx   = (r0c - r1c) * piaacmcopticaldesign.beamrad;
        double dz   = piaaM1z[i] - piaaM0z[i];
        double dist = dz * sqrt(1.0 + (dx / dz) * (dx / dz));
        double y3   = dist - dz;

        double slope;
        if (fabs(dx) > 0.000000001)
        {
            slope = y3 / dx;
        }
        else
        {
            slope = 0.0;
        }

        double r0n     = piaacmcopticaldesign.r0lim * (i + 1) / NBradpts;
        piaaM0z[i + 1] = piaaM0z[i] + slope * (r0n - r0c) * piaacmcopticaldesign.beamrad;
        piaaM1z[i + 1] = piaaM1z[i] + slope * (piaar10[i + 1] - r1c) * piaacmcopticaldesign.beamrad;
    }
}

/**
 * write_piaa_mshapes_table - Write 1D physical radii and axial sags to `PIAA_Mshapes.txt`.
 * @piaar00:  Dimensionless radii on PIAA element 0
 * @piaaM0z:  Axial sag profile [m] on PIAA element 0
 * @piaar10:  Mapped dimensionless radii on PIAA element 1
 * @piaaM1z:  Axial sag profile [m] on PIAA element 1
 * @NBradpts: Number of 1D radial sample points
 *
 * Return: RETURN_SUCCESS or failure if output file cannot be opened.
 */
static errno_t write_piaa_mshapes_table(
    const double *restrict piaar00,
    const double *restrict piaaM0z,
    const double *restrict piaar10,
    const double *restrict piaaM1z,
    long                   NBradpts)
{
    char fname[STRINGMAXLEN_FULLFILENAME];
    WRITE_FULLFILENAME(fname, "%s/PIAA_Mshapes.txt", piaacmcparams.piaacmcconfdir);

    FILE *fp = fopen(fname, "w");
    if (fp == NULL)
    {
        FUNC_RETURN_FAILURE("cannot write file %s", fname);
    }

    for (long ii = 0; ii < NBradpts; ii++)
    {
        fprintf(fp,
                "%18.16f %18.16f %18.16f %18.16f\n",
                piaar00[ii] * piaacmcopticaldesign.beamrad,
                piaaM0z[ii],
                piaar10[ii] * piaacmcopticaldesign.beamrad,
                piaaM1z[ii]);
    }
    fclose(fp);
    return RETURN_SUCCESS;
}

/**
 * geom_piaa_rad_march_and_save - Compute r0-r1 ray mapping, integrate sags, and save table.
 * @ctx: Initialized radial amplitude and cumulative enclosed flux profiles
 *
 * Return: RETURN_SUCCESS or file I/O failure status.
 */
errno_t geom_piaa_rad_march_and_save(const struct geom_piaa_rad_ctx *ctx)
{
    DEBUG_TRACE_FSTART();

    long    NBradpts = ctx->NBradpts;
    double *piaar00  = (double *) malloc(sizeof(double) * NBradpts);
    double *piaar10  = (double *) malloc(sizeof(double) * NBradpts);
    double *piaar01  = (double *) malloc(sizeof(double) * NBradpts);
    double *piaaM0z  = (double *) malloc(sizeof(double) * NBradpts);
    double *piaaM1z  = (double *) malloc(sizeof(double) * NBradpts);
    if ((piaar00 == NULL) || (piaar10 == NULL) || (piaar01 == NULL) || (piaaM0z == NULL) ||
        (piaaM1z == NULL))
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    for (long i = 0; i < NBradpts; i++)
    {
        piaar00[i] = piaacmcopticaldesign.r0lim * i / NBradpts;
    }
    piaar00[0] = 0.0;

    invert_cumulative_flux_map(ctx->flux0cumul,
                               ctx->flux1cumul,
                               NBradpts,
                               piaacmcopticaldesign.r1lim,
                               piaar10);
    invert_cumulative_flux_map(ctx->flux1cumul,
                               ctx->flux0cumul,
                               NBradpts,
                               piaacmcopticaldesign.r0lim,
                               piaar01);

    printf("======== Compute PIAA optics shapes ============\n");
    march_piaa_mirror_sags(piaar00, piaar10, NBradpts, piaaM0z, piaaM1z);

    errno_t status = write_piaa_mshapes_table(piaar00, piaaM0z, piaar10, piaaM1z, NBradpts);

    free(piaaM0z);
    free(piaaM1z);
    free(piaar00);
    free(piaar10);
    free(piaar01);

    if (status != RETURN_SUCCESS)
    {
        return status;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
