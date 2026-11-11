/**
 * @file    mkFPM_zonemap.c
 * @brief   Generate 2D focal plane mask hexagonal ring/sector zone index map
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "image_gen/image_gen.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"

#include "FocalPlaneMask/mkFPM_zonemap.h"

#define FPM_ZONEMAP_MAX_HEX 10000

struct fpm_zonemap_ctx
{
    uint_fast32_t *nbsector;
    uint_fast32_t *nbsectorcumul;
    uint_fast32_t  nbzonescc;
    double         hex_x[FPM_ZONEMAP_MAX_HEX];
    double         hex_y[FPM_ZONEMAP_MAX_HEX];
    long           hex_ring[FPM_ZONEMAP_MAX_HEX];
    uint_fast32_t  hex_number[FPM_ZONEMAP_MAX_HEX];
    uint_fast32_t  hindexMax;
};

/**
 * zonemap_init_sectors - Build polar ring/sector counts and write zone coordinates file
 * @ctx: Zone-map context receiving per-ring and cumulative sector counts.
 *
 * Return: RETURN_SUCCESS on file creation, RETURN_FAILURE otherwise.
 */
static errno_t zonemap_init_sectors(struct fpm_zonemap_ctx *ctx)
{
    char fname[STRINGMAXLEN_FULLFILENAME];

    WRITE_FULLFILENAME(
        fname,
        "%s/fpm_zonescoord_%d_%03ld.txt",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.PIAACMC_FPMsectors,
        piaacmcopticaldesign.NBrings);

    FILE *fp = fopen(fname, "w");
    if (fp == NULL)
    {
        FUNC_RETURN_FAILURE("Cannot create file %s", fname);
    }

    uint_fast64_t cnt    = 0;
    ctx->nbsector[0]      = 1;
    ctx->nbsectorcumul[0] = 1;
    fprintf(fp, "0 0\n");

    for (long ring = 1; ring < piaacmcopticaldesign.NBrings; ring++)
    {
        ctx->nbsector[ring]      = (uint_fast32_t) (2 * (ring + 1));
        ctx->nbsectorcumul[ring] = ctx->nbsectorcumul[ring - 1] + ctx->nbsector[ring];
        for (uint_fast64_t cnt1 = 0; cnt1 < ctx->nbsector[ring]; cnt1++)
        {
            cnt++;
            fprintf(fp, "%ld %ld\n", (long) cnt, ring);
        }
    }
    fclose(fp);

    for (long ring = 0; ring < piaacmcopticaldesign.NBrings; ring++)
    {
        printf(
            "ring %ld : %ld %ld\n",
            ring,
            (long) ctx->nbsector[ring],
            (long) ctx->nbsectorcumul[ring]);
    }

    return RETURN_SUCCESS;
}

/**
 * zonemap_enumerate_hex_grid - Populate hexagonal lattice centers within NBrings
 * @ctx: Zone-map context receiving hexagon center coordinates and ring assignments.
 */
static void zonemap_enumerate_hex_grid(struct fpm_zonemap_ctx *ctx)
{
    const double hexstep = 1.0;

    ctx->nbsector[0]      = 1;
    ctx->nbsectorcumul[0] = 1;
    for (long ring = 1; ring < piaacmcopticaldesign.NBrings; ring++)
    {
        ctx->nbsector[ring]      = 0;
        ctx->nbsectorcumul[ring] = 0;
    }

    uint_fast32_t hindex = 0;
    long          ii1max = (long) (piaacmcopticaldesign.NBrings / 3 + 2);
    long          jj1max = (long) (piaacmcopticaldesign.NBrings / sqrt(3.0) + 2);

    for (long ii1 = -ii1max; ii1 < ii1max; ii1++)
    {
        for (long jj1 = -jj1max; jj1 < jj1max; jj1++)
        {
            double hx   = hexstep * ii1 * 3.0;
            double hy   = hexstep * sqrt(3.0) * jj1;
            long   ring = (long) sqrt(hx * hx + hy * hy);
            if (ring < piaacmcopticaldesign.NBrings)
            {
                ctx->nbsector[ring]++;
                ctx->hex_x[hindex]    = hx;
                ctx->hex_y[hindex]    = hy;
                ctx->hex_ring[hindex] = ring;
                hindex++;
            }

            hx += hexstep * 1.5;
            hy += hexstep * sqrt(3.0) / 2.0;
            ring = (long) sqrt(hx * hx + hy * hy);
            if (ring < piaacmcopticaldesign.NBrings)
            {
                ctx->nbsector[ring]++;
                ctx->hex_x[hindex]    = hx;
                ctx->hex_y[hindex]    = hy;
                ctx->hex_ring[hindex] = ring;
                hindex++;
            }
        }
    }

    ctx->hindexMax = hindex;
}

/**
 * zonemap_write_hex_coords - Number hexagonal rings outward and write coordinate table
 * @ctx: Zone-map context updated with sequential hexagon indices and cumulative counts.
 *
 * Return: RETURN_SUCCESS on file write, RETURN_FAILURE otherwise.
 */
static errno_t zonemap_write_hex_coords(struct fpm_zonemap_ctx *ctx)
{
    char fname[STRINGMAXLEN_FULLFILENAME];

    WRITE_FULLFILENAME(
        fname,
        "%s/fpm_zonescoord_%d_%03ld.txt",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.PIAACMC_FPMsectors,
        piaacmcopticaldesign.NBrings);

    FILE *fp = fopen(fname, "w");
    if (fp == NULL)
    {
        FUNC_RETURN_FAILURE("Cannot create file %s", fname);
    }

    fprintf(fp, "# focal plane mask zones geometry\n");
    fprintf(fp, "# hexagonal tiling\n");
    fprintf(fp, "# col 1: hexagon index\n");
    fprintf(fp, "# col 2: ring index\n");
    fprintf(fp, "# col 3: hexagon center x coordinate\n");
    fprintf(fp, "# col 4: hexagon center y coordinate\n");
    fprintf(fp, "# Note: unit = hexagon center to tip radius (circumradius)\n");
    fprintf(fp, "# \n");
    fprintf(fp, "%5ld %5ld  %11.6f %11.6f\n", 0L, 0L, 0.0, 0.0);

    uint_fast32_t hcnt = 1;
    for (long ring = 1; ring < piaacmcopticaldesign.NBrings; ring++)
    {
        for (uint_fast32_t hindex = 0; hindex < ctx->hindexMax; hindex++)
        {
            if (ctx->hex_ring[hindex] == ring)
            {
                ctx->hex_number[hindex] = hcnt;
                fprintf(
                    fp,
                    "%5ld %5ld  %11.6f %11.6f\n",
                    (long) hcnt,
                    ring,
                    ctx->hex_x[hindex],
                    ctx->hex_y[hindex]);
                hcnt++;
            }
        }
        ctx->nbsectorcumul[ring] = ctx->nbsectorcumul[ring - 1] + ctx->nbsector[ring];
    }
    fclose(fp);

    for (long ring = 0; ring < piaacmcopticaldesign.NBrings; ring++)
    {
        printf(
            "ring %ld : %ld %ld\n",
            ring,
            (long) ctx->nbsector[ring],
            (long) ctx->nbsectorcumul[ring]);
    }

    return RETURN_SUCCESS;
}

/**
 * zonemap_setup_rings - Initialize ring/sector/hexagonal geometry tables for active mode
 * @ctx: Zone-map context to populate.
 *
 * Return: RETURN_SUCCESS on valid setup, error code otherwise.
 */
static errno_t zonemap_setup_rings(struct fpm_zonemap_ctx *ctx)
{
    switch (piaacmcparams.PIAACMC_FPMsectors)
    {
    case 0:
        ctx->nbsectorcumul[0] = 1;
        for (long ring = 1; ring < piaacmcopticaldesign.NBrings; ring++)
        {
            ctx->nbsectorcumul[ring] = ctx->nbsectorcumul[ring - 1] + 1;
        }
        break;

    case 1:
        FUNC_CHECK_RETURN(zonemap_init_sectors(ctx));
        break;

    case 2:
        zonemap_enumerate_hex_grid(ctx);
        FUNC_CHECK_RETURN(zonemap_write_hex_coords(ctx));
        break;

    default:
        printf("ERROR: FPMsector mode (%d) not recognized\n", piaacmcparams.PIAACMC_FPMsectors);
        exit(0);
    }

    if (piaacmcopticaldesign.NBringCentCone > 0)
    {
        ctx->nbzonescc = ctx->nbsectorcumul[piaacmcopticaldesign.NBringCentCone - 1];
    }
    else
    {
        ctx->nbzonescc = 0;
    }

    return RETURN_SUCCESS;
}

/**
 * zonemap_compute_polar_zone - Evaluate ring or polar-sector zone index at normalized (x, y)
 * @ctx: Initialized zone-map geometry context.
 * @x:   Normalized mask horizontal coordinate.
 * @y:   Normalized mask vertical coordinate.
 *
 * Return: 1-based zone index inside active FPM rings, or 0 outside/central-cone.
 */
static uint16_t zonemap_compute_polar_zone(
    const struct fpm_zonemap_ctx *ctx,
    double                        x,
    double                        y)
{
    const double eps = 1.0e-6;
    double       r   = sqrt(x * x + y * y);
    double       PA  = atan2(y, x);
    double       PAf = 0.5 * ((PA / M_PI) + 1.0);

    if (PAf < eps)
    {
        PAf = eps;
    }
    if (PAf > 1.0 - eps)
    {
        PAf = 1.0 - eps;
    }

    long zi = (long) ceil((1.0 - r) * piaacmcopticaldesign.NBrings);
    if (zi <= 0)
    {
        zi = 0;
    }
    if (zi > piaacmcopticaldesign.NBrings)
    {
        zi = piaacmcopticaldesign.NBrings;
    }

    if (piaacmcparams.PIAACMC_FPMsectors == 0)
    {
        return (zi == 0) ? 0U : (uint16_t) (piaacmcopticaldesign.NBrings - zi + 1);
    }

    if (zi == 0)
    {
        return 0U;
    }

    long ring      = piaacmcopticaldesign.NBrings - zi;
    long zoneindex = 1;
    if (ring > 0)
    {
        zoneindex = (long) ctx->nbsectorcumul[ring - 1] + 1 + (long) (PAf * ctx->nbsector[ring]);
    }

    if (piaacmcopticaldesign.NBrings > 1)
    {
        if (ring < piaacmcopticaldesign.NBringCentCone)
        {
            zoneindex = 0;
        }
        else
        {
            zoneindex -= (long) ctx->nbzonescc;
        }
    }

    return (uint16_t) zoneindex;
}

/**
 * zonemap_fill_polar - Rasterize concentric ring or polar-sector zone indices onto ID
 * @ctx: Initialized zone-map geometry context.
 * @ID:  Destination 2D UINT16 zone-map image identifier.
 */
static void zonemap_fill_polar(
    const struct fpm_zonemap_ctx *ctx,
    imageID                       ID)
{
    long   arrsize = piaacmcopticaldesign.fpmarraysize;
    double scale   = (double) arrsize * piaacmcparams.FPMSCALEFACTOR;

    for (long ii = 0; ii < arrsize; ii++)
    {
        double x = (2.0 * ii - 1.0 * arrsize) / scale;
        for (long jj = 0; jj < arrsize; jj++)
        {
            double y = (2.0 * jj - 1.0 * arrsize) / scale;

            dcimg[ID].array.UI16[jj * arrsize + ii] = zonemap_compute_polar_zone(ctx, x, y);
        }
    }

    if (piaacmcopticaldesign.NBrings <= 1)
    {
        piaacmcopticaldesign.focmNBzone = 1;
    }
    else if (piaacmcparams.PIAACMC_FPMsectors == 0)
    {
        piaacmcopticaldesign.focmNBzone = piaacmcopticaldesign.NBrings - (long) ctx->nbzonescc;
    }
    else
    {
        piaacmcopticaldesign.focmNBzone =
            (long) ctx->nbsectorcumul[piaacmcopticaldesign.NBrings - 1] - (long) ctx->nbzonescc;
    }
}

/**
 * zonemap_fill_hex - Rasterize hexagonal FPM tiling zones onto ID
 * @ctx: Initialized hexagonal zone-map context.
 * @ID:  Destination 2D UINT16 zone-map image identifier.
 */
static void zonemap_fill_hex(
    const struct fpm_zonemap_ctx *ctx,
    imageID                       ID)
{
    const double hexgap     = -0.0001;
    long         arrsize    = piaacmcopticaldesign.fpmarraysize;
    long         npix       = arrsize * arrsize;
    double       hexsteppix = 0.5 * arrsize / piaacmcopticaldesign.NBrings
                              * piaacmcparams.FPMSCALEFACTOR;
    double       hexrad     = hexsteppix * (1.0 - hexgap) * (sqrt(3.0) / 2.0);

    for (uint_fast32_t hindex = 0; hindex < ctx->hindexMax; hindex++)
    {
        imageID ID1 = make_hexagon(
            "_TMPhex",
            arrsize,
            arrsize,
            0.5 * arrsize + ctx->hex_x[hindex] * hexsteppix,
            0.5 * arrsize + ctx->hex_y[hindex] * hexsteppix,
            hexrad);

        uint16_t zval = (uint16_t) (ctx->hex_number[hindex] + 1);
        for (long ii = 0; ii < npix; ii++)
        {
            if (dcimg[ID1].array.F[ii] > 0.5f)
            {
                dcimg[ID].array.UI16[ii] = zval;
            }
        }
        delete_image_ID("_TMPhex", DELETE_IMAGE_ERRMODE_WARNING);
    }

    if (piaacmcopticaldesign.NBrings > 1)
    {
        piaacmcopticaldesign.focmNBzone =
            (long) ctx->nbsectorcumul[piaacmcopticaldesign.NBrings - 1];
    }
    else
    {
        piaacmcopticaldesign.focmNBzone = 1;
    }
}

/**
 * mkFPM_zonemap - Generate 2D UINT16 focal plane mask zone index map
 * @IDname: Name of output 2D zone-map image to create in shared/local memory.
 * @outID:  Optional output pointer receiving created image identifier.
 *
 * Supports concentric rings (mode 0), polar ring sectors (mode 1), and hexagonal
 * tilings (mode 2), updating piaacmcopticaldesign.focmNBzone accordingly.
 *
 * Return: RETURN_SUCCESS on completion, error code otherwise.
 */
errno_t mkFPM_zonemap(
    const char *restrict IDname,
    imageID             *outID)
{
    DEBUG_TRACE_FSTART();

    printf("function PIAACMCsimul_mkFPM_zonemap\n");
    fflush(stdout);

    uint32_t sizearray[2];
    sizearray[0] = (uint32_t) piaacmcopticaldesign.fpmarraysize;
    sizearray[1] = (uint32_t) piaacmcopticaldesign.fpmarraysize;

    imageID ID;
    create_image_ID(IDname, 2, sizearray, _DATATYPE_UINT16, 0, 0, 0, &ID);

    struct fpm_zonemap_ctx ctx;
    memset(&ctx, 0, sizeof(ctx));

    size_t nrings     = (size_t) piaacmcopticaldesign.NBrings;
    ctx.nbsector      = (uint_fast32_t *) malloc(sizeof(uint_fast32_t) * nrings);
    ctx.nbsectorcumul = (uint_fast32_t *) malloc(sizeof(uint_fast32_t) * nrings);
    if (ctx.nbsector == NULL || ctx.nbsectorcumul == NULL)
    {
        free(ctx.nbsector);
        free(ctx.nbsectorcumul);
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    errno_t status = zonemap_setup_rings(&ctx);
    if (status != RETURN_SUCCESS)
    {
        free(ctx.nbsector);
        free(ctx.nbsectorcumul);
        return status;
    }

    if (piaacmcparams.PIAACMC_FPMsectors == 0 || piaacmcparams.PIAACMC_FPMsectors == 1)
    {
        zonemap_fill_polar(&ctx, ID);
    }
    else if (piaacmcparams.PIAACMC_FPMsectors == 2)
    {
        zonemap_fill_hex(&ctx, ID);
    }

    printf(
        "[%d] piaacmcopticaldesign.focmNBzone  =  %ld %ld    %ld %ld   ->  %ld   (%ld)\n",
        piaacmcparams.PIAACMC_FPMsectors,
        piaacmcopticaldesign.NBrings,
        (long) ctx.nbsectorcumul[piaacmcopticaldesign.NBrings - 1],
        piaacmcopticaldesign.NBringCentCone,
        (long) ctx.nbzonescc,
        piaacmcopticaldesign.focmNBzone,
        piaacmcopticaldesign.NBrings);

    if (piaacmcparams.PIAACMC_FPMsectors != 0)
    {
        printf("Saving %s ....\n", IDname);
        save_fits(IDname, "__test_zonemap_00.fits");
    }

    free(ctx.nbsector);
    free(ctx.nbsectorcumul);

    printf("NUMBER OF ZONES = %ld\n", piaacmcopticaldesign.focmNBzone);

    if (outID != NULL)
    {
        *outID = ID;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
