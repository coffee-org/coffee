/**
 * @file    mkPIAAMshapes_from_RadSag.c
 * @brief   PIAA-type coronagraph design, make PIAA shapes from radial sag
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "PIAAshape/mkPIAAMshapes_from_RadSag.h"

struct rad_sag_table
{
    long    NBradpts;
    double *r0;
    double *z0;
    double *r1;
    double *z1;
};

/**
 * free_radial_sag_table - Release 1D radial sag profile buffers.
 * @tbl: Table structure whose heap arrays should be freed
 */
static void free_radial_sag_table(struct rad_sag_table *tbl)
{
    free(tbl->r0);
    free(tbl->z0);
    free(tbl->r1);
    free(tbl->z1);
}

/**
 * load_radial_sag_table - Allocate and read 4-column radial PIAA sag table from ASCII file.
 * @fname:    Path to 1D radial sag file (`r0 z0 r1 z1` per line)
 * @PIAAsep:  Axial separation [m] between PIAA elements subtracted from `z1`
 * @NBradpts: Number of radial sample rows to read
 * @tbl:      Output table receiving allocated radial and sag arrays
 *
 * Return: RETURN_SUCCESS or failure if `fname` cannot be opened.
 */
static errno_t load_radial_sag_table(
    const char  *restrict fname,
    double                PIAAsep,
    long                  NBradpts,
    struct rad_sag_table *tbl)
{
    tbl->NBradpts = NBradpts;
    tbl->r0       = (double *) malloc(sizeof(double) * NBradpts);
    tbl->z0       = (double *) malloc(sizeof(double) * NBradpts);
    tbl->r1       = (double *) malloc(sizeof(double) * NBradpts);
    tbl->z1       = (double *) malloc(sizeof(double) * NBradpts);
    if ((tbl->r0 == NULL) || (tbl->z0 == NULL) || (tbl->r1 == NULL) || (tbl->z1 == NULL))
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    DEBUG_TRACEPOINT("read sag radial profile %s", fname);
    FILE *fp = fopen(fname, "r");
    if (fp == NULL)
    {
        free_radial_sag_table(tbl);
        FUNC_RETURN_FAILURE("Cannot read file %s", fname);
    }

    for (long k = 0; k < NBradpts; k++)
    {
        int ret =
            fscanf(fp, "%lf %lf %lf %lf\n", &tbl->r0[k], &tbl->z0[k], &tbl->r1[k], &tbl->z1[k]);
        (void) ret;
    }
    fclose(fp);

    for (long k = 0; k < NBradpts; k++)
    {
        tbl->z1[k] -= PIAAsep;
    }

    return RETURN_SUCCESS;
}

/**
 * interp_radial_sag - Linearly interpolate 1D sag value at physical radius `r`.
 * @rarr:     Monotonically increasing physical radius array [m]
 * @zarr:     Corresponding axial sag array [m]
 * @NBradpts: Number of sample points in `rarr` and `zarr`
 * @r:        Query physical radius [m]
 *
 * Return: Interpolated axial sag [m] at radius `r`.
 */
static inline double interp_radial_sag(
    const double *restrict rarr,
    const double *restrict zarr,
    long                   NBradpts,
    double                 r)
{
    long k = 1;
    while ((k < NBradpts - 2) && (rarr[k] < r))
    {
        k++;
    }

    double r00   = rarr[k - 1];
    double r01   = rarr[k];
    double alpha = (r - r00) / (r01 - r00);
    if (alpha > 1.0)
    {
        alpha = 1.0;
    }

    return (1.0 - alpha) * zarr[k - 1] + alpha * zarr[k];
}

/**
 * mkPIAAMshapes_from_RadSag - Interpolate 1D radial sag profile onto 2D PIAA OPD/sag screens.
 * @piaa1Dsagfname: Path to ASCII 1D radial sag table (`PIAA_Mshapes.txt`)
 * @PIAAsep:        Physical separation between PIAA surfaces [m]
 * @beamrad:        Physical beam radius [m]
 * @r0limfact:      Input beam extrapolation radius factor
 * @r1limfact:      Output beam extrapolation radius factor
 * @size:           Output 2D square image linear dimension [pix]
 * @beamradpix:     Beam radius [pix]
 * @NBradpts:       Number of radial sample points in `piaa1Dsagfname`
 * @ID_PIAAM0_name: Output 2D sag image name for PIAA surface 0
 * @ID_PIAAM1_name: Output 2D sag image name for PIAA surface 1
 *
 * Return: RETURN_SUCCESS on completion, or milk error code on file/image failure.
 */
errno_t mkPIAAMshapes_from_RadSag(
    const char *restrict piaa1Dsagfname,
    double               PIAAsep,
    double               beamrad,
    double               r0limfact,
    double               r1limfact,
    uint32_t             size,
    double               beamradpix,
    long                 NBradpts,
    const char *restrict ID_PIAAM0_name,
    const char *restrict ID_PIAAM1_name)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %s %u %f %ld %s %s",
                     piaa1Dsagfname,
                     size,
                     beamradpix,
                     NBradpts,
                     ID_PIAAM0_name,
                     ID_PIAAM1_name);

    struct rad_sag_table tbl = {0};
    FUNC_CHECK_RETURN(load_radial_sag_table(piaa1Dsagfname, PIAAsep, NBradpts, &tbl));

    imageID ID_PIAAM0;
    imageID ID_PIAAM1;
    errno_t status = create_2Dimage_ID(ID_PIAAM0_name, size, size, &ID_PIAAM0);
    if (status == RETURN_SUCCESS)
    {
        status = create_2Dimage_ID(ID_PIAAM1_name, size, size, &ID_PIAAM1);
    }
    if (status != RETURN_SUCCESS)
    {
        free_radial_sag_table(&tbl);
        return status;
    }

    printf("\n\n");

    long   size_l = (long) size;
    double r0max  = r0limfact * beamrad;
    double r1max  = r1limfact * beamrad;

#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp for
#endif
        for (long ii = 0; ii < size_l; ii++)
        {
            for (long jj = 0; jj < size_l; jj++)
            {
                double x   = (1.0 * ii - 0.5 * size) / beamradpix;
                double y   = (1.0 * jj - 0.5 * size) / beamradpix;
                double r   = sqrt(x * x + y * y) * beamrad;
                long   idx = jj * size_l + ii;

                if (r < r0max)
                {
                    dcimg[ID_PIAAM0].array.F[idx] = interp_radial_sag(tbl.r0, tbl.z0, NBradpts, r);
                }
                else
                {
                    dcimg[ID_PIAAM0].array.F[idx] = 0.0;
                }

                if (r < r1max)
                {
                    dcimg[ID_PIAAM1].array.F[idx] = -interp_radial_sag(tbl.r1, tbl.z1, NBradpts, r);
                }
                else
                {
                    dcimg[ID_PIAAM1].array.F[idx] = 0.0;
                }
            }
        }
#ifdef _OPENMP
    } // omp parallel
#endif

    char fnametest[STRINGMAXLEN_FILENAME];
    WRITE_FILENAME(fnametest, "test_piaam0z.fits");
    status = save_fl_fits(ID_PIAAM0_name, fnametest);

    printf("\n\n");
    free_radial_sag_table(&tbl);

    if (status != RETURN_SUCCESS)
    {
        return status;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
