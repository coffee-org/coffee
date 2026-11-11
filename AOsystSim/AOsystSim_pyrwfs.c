/**
 * @file    AOsystSim_pyrwfs.c
 * @brief   Configuration generator and real-time shared-memory daemon for modulated Pyramid WFS
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "fft/fft.h"
#include "statistic/statistic.h"

#include "AOsystSim/AOsystSim.h"
#include "AOsystSim_internal.h"

/**
 * @brief Runtime configuration and stream state for the Pyramid WFS daemon.
 */
struct pyrwfs_daemon_ctx
{
    float   LAMBDA;
    float   WFSFLUX;
    float   WFSCAMRON;
    int     INMODE;
    char    INSTREAMNAMEOPD[200];
    char    INSTREAMNAMEAMP[200];
    int     INSEMCHAN;
    char    INFITSFILENAMEOPD[200];
    char    INFITSFILENAMEAMP[200];
    char    INTRIGGERFILE[200];
    int     OUTMODE;
    char    OUTSTREAMNAME[200];
    char    OUTINSTSTREAMNAME[200];
    char    OUTFITSFILENAME[200];
    char    OUTTRIGGERFILE[200];
    long    ARRAYSIZE;
    int     PYRMODMODE;
    float   PYRMODAMP;
    long    PYRMODNBPT;
    float   PYRAPERTURE;
    long    PUPPIXDIAM;
    float   PYRPUPSEP;
    long    OUTBINFACT;
    long    OUTARRAYSIZE;
    float   fpscale;
    long    wfinsize;
    long    offset;
    imageID IDinOPD;
    imageID IDinAMP;
    imageID IDout;
    imageID IDoutinst;
    imageID IDpyr_amp;
    imageID IDpyr_pha;
    imageID IDpupamp;
    imageID IDpuppha;
    imageID IDfocamp;
    imageID IDfocpha;
    imageID IDpyrpupi;
};

/**
 * AOsystSim_PyrWFS_mkCONF - Write template configuration file for real-time Pyramid WFS daemon
 * @fname: Output configuration file path
 *
 * Return: 0 on completion.
 */
int AOsystSim_PyrWFS_mkCONF(const char *fname)
{
    FILE *fp = fopen(fname, "w");

    fprintf(fp, "\n");
    fprintf(fp, "LAMBDANM                 800              # wavelength [nm]\n");
    fprintf(fp, "WFSFLUX              150000.0              # flux per frame [photon]\n");
    fprintf(fp, "WFSCAMRON                1.0               # WFS camera RON [e-]\n\n");
    fprintf(fp, "# ============== INPUT TYPE (OPD unit = um) ====================\n");
    fprintf(fp, "INMODE                   0                 "
                "# 0:stream, 1:file system (FITS)\n");
    fprintf(fp, "INSTREAMNAMEOPD          wf1opd            # input OPD stream\n");
    fprintf(fp, "INSTREAMNAMEAMP          wf1amp            # input AMP stream\n");
    fprintf(fp, "INSEMCHAN                5                 "
                "# input semaphore channel (using OPD input)\n");
    fprintf(fp, "INFITSFILENAMEOPD        inwfopd.fits      # input FITS file name (OPD)\n");
    fprintf(fp, "INFITSFILENAMEAMP        inwfamp.fits      # input FITS file name (AMP)\n");
    fprintf(fp, "INTRIGGERFILE            inwf.txt          # input trigger file\n\n");
    fprintf(fp, "# ============== OUTPUT TYPE ===================================\n");
    fprintf(fp, "OUTMODE                  0                 # 0: stream, 1: file system\n");
    fprintf(fp, "OUTSTREAMNAME            pWFSim            # output image stream name \n");
    fprintf(fp, "OUTINSTSTREAMNAME        WFSinst           "
                "# output instantaneous image stream name \n");
    fprintf(fp, "OUTFITSFILENAME          PyrWFS_im.fits    "
                "# FITS file name output (intensity)\n");
    fprintf(fp, "OUTTRIGGERFILE           outpyr.txt        # output trigger file\n\n");
    fprintf(fp, "# ============== FREQUENCY, GEOMETRY =============================\n");
    fprintf(fp, "ARRAYSIZE                256               # array size for computations\n");
    fprintf(fp, "PYRMODMODE               1                 "
                "# 0: no modulation, 1: modulation\n");
    fprintf(fp, "PYRMODAMP                4.0               # Modulation radius [l/D]\n");
    fprintf(fp, "PYRMODNBPT               8                 "
                "# number of pyramid modulation points\n");
    fprintf(fp, "PYRAPERTURE              50.0              "
                "# spatial filter aperture [l/D]\n");
    fprintf(fp, "PUPPIXDIAM               100.0             "
                "# pupil diameter [pix] - used to compute l/D scale\n");
    fprintf(fp, "PYRPUPSEP                1.2               "
                "# separation between pupil imagees (relative to pup diam)\n");
    fprintf(fp, "OUTBINFACT               2                 # output binning factor\n");
    fprintf(fp, "OUTARRAYSIZE             120               # output array size\n\n");
    fclose(fp);

    return 0;
}

/**
 * pyrwfs_read_conf - Read Pyramid WFS daemon configuration parameters from disk
 * @conf_fname: Configuration file path
 * @ctx:        Pyramid WFS daemon context
 */
static void pyrwfs_read_conf(
    const char               *conf_fname,
    struct pyrwfs_daemon_ctx *ctx)
{
    FILE *fp;
    char  fname[STRINGMAXLEN_FULLFILENAME];

    if ((fp = fopen(conf_fname, "r")) == NULL)
    {
        sprintf(fname, "%s.default", conf_fname);
        printf("configuration file %s not found. Creating default template as %s\n",
               conf_fname,
               fname);
        AOsystSim_PyrWFS_mkCONF(fname);
        exit(0);
    }
    fclose(fp);

    ctx->LAMBDA    = 1.0e-9f * read_config_parameter_float(conf_fname, "LAMBDANM");
    ctx->WFSFLUX   = read_config_parameter_float(conf_fname, "WFSFLUX");
    ctx->WFSCAMRON = read_config_parameter_float(conf_fname, "WFSCAMRON");

    ctx->INMODE = read_config_parameter_long(conf_fname, "INMODE");
    read_config_parameter(conf_fname, "INSTREAMNAMEOPD", ctx->INSTREAMNAMEOPD);
    read_config_parameter(conf_fname, "INSTREAMNAMEAMP", ctx->INSTREAMNAMEAMP);
    ctx->INSEMCHAN = read_config_parameter_long(conf_fname, "INSEMCHAN");
    read_config_parameter(conf_fname, "INFITSFILENAMEOPD", ctx->INFITSFILENAMEOPD);
    read_config_parameter(conf_fname, "INFITSFILENAMEAMP", ctx->INFITSFILENAMEAMP);
    read_config_parameter(conf_fname, "INTRIGGERFILE", ctx->INTRIGGERFILE);

    ctx->OUTMODE = read_config_parameter_long(conf_fname, "OUTMODE");
    read_config_parameter(conf_fname, "OUTSTREAMNAME", ctx->OUTSTREAMNAME);
    read_config_parameter(conf_fname, "OUTINSTSTREAMNAME", ctx->OUTINSTSTREAMNAME);
    read_config_parameter(conf_fname, "OUTFITSFILENAME", ctx->OUTFITSFILENAME);
    read_config_parameter(conf_fname, "OUTTRIGGERFILE", ctx->OUTTRIGGERFILE);

    ctx->ARRAYSIZE    = read_config_parameter_long(conf_fname, "ARRAYSIZE");
    ctx->PYRMODMODE   = read_config_parameter_long(conf_fname, "PYRMODMODE");
    ctx->PYRMODAMP    = read_config_parameter_float(conf_fname, "PYRMODAMP");
    ctx->PYRMODNBPT   = read_config_parameter_long(conf_fname, "PYRMODNBPT");
    ctx->PYRAPERTURE  = read_config_parameter_float(conf_fname, "PYRAPERTURE");
    ctx->PUPPIXDIAM   = (long) read_config_parameter_float(conf_fname, "PUPPIXDIAM");
    ctx->PYRPUPSEP    = read_config_parameter_float(conf_fname, "PYRPUPSEP");
    ctx->OUTBINFACT   = read_config_parameter_long(conf_fname, "OUTBINFACT");
    ctx->OUTARRAYSIZE = read_config_parameter_long(conf_fname, "OUTARRAYSIZE");
    ctx->fpscale      = 1.0f * ctx->ARRAYSIZE / ctx->PUPPIXDIAM;

    if (ctx->PYRMODMODE == 0)
    {
        ctx->PYRMODNBPT = 1;
        ctx->PYRMODAMP  = 0.0f;
    }
}

/**
 * pyrwfs_init_masks - Synthesize 3D focal-plane pyramid amplitude and phase modulation cubes
 * @ctx: Pyramid WFS daemon context
 */
static void pyrwfs_init_masks(struct pyrwfs_daemon_ctx *ctx)
{
    long arr = ctx->ARRAYSIZE;

    printf("PYRAPERTURE = %f  -> %f pix\n", ctx->PYRAPERTURE, ctx->PYRAPERTURE * ctx->fpscale);
    printf("PYTPUPSEP = %f  ->  %f\n", ctx->PYRPUPSEP, M_PI * ctx->PUPPIXDIAM * ctx->PYRPUPSEP);

    create_3Dimage_ID("pyramp", arr, arr, ctx->PYRMODNBPT, &ctx->IDpyr_amp);
    create_3Dimage_ID("pyrpha", arr, arr, ctx->PYRMODNBPT, &ctx->IDpyr_pha);

    for (long pmodpt = 0; pmodpt < ctx->PYRMODNBPT; pmodpt++)
    {
        float PA    = (float) (2.0 * M_PI * pmodpt / ctx->PYRMODNBPT);
        float xc    = ctx->PYRMODAMP * cosf(PA);
        float yc    = ctx->PYRMODAMP * sinf(PA);
        float slope = (float) (M_PI * ctx->PUPPIXDIAM * ctx->PYRPUPSEP / arr);
        float ap    = ctx->PYRAPERTURE * ctx->fpscale;

        for (long ii = 0; ii < arr; ii++)
        {
            for (long jj = 0; jj < arr; jj++)
            {
                float x   = 1.0f * (ii - arr / 2) - xc * ctx->fpscale;
                float y   = 1.0f * (jj - arr / 2) - yc * ctx->fpscale;
                long  idx = pmodpt * arr * arr + jj * arr + ii;

                dcimg[ctx->IDpyr_pha].array.F[idx] = slope * (fabsf(x) + fabsf(y));
                if ((fabsf(1.0f * ii - arr / 2) > ap) || (fabsf(1.0f * jj - arr / 2) > ap))
                {
                    dcimg[ctx->IDpyr_amp].array.F[idx] = 0.0f;
                }
                else
                {
                    dcimg[ctx->IDpyr_amp].array.F[idx] = 1.0f;
                }
            }
        }
    }
    save_fits("pyrpha", "./AOsystSim_wdir/pyrpha.fits");
    save_fits("pyramp", "./AOsystSim_wdir/pyramp.fits");
}

/**
 * pyrwfs_init_streams - Attach input wavefront streams and allocate output/working buffers
 * @ctx: Pyramid WFS daemon context
 */
static void pyrwfs_init_streams(struct pyrwfs_daemon_ctx *ctx)
{
    uint32_t sizearray[2];

    switch (ctx->INMODE)
    {
        case 0:
            ctx->IDinOPD = read_sharedmem_image(ctx->INSTREAMNAMEOPD, dcimg, dcnimg);
            while (ctx->IDinOPD == -1)
            {
                usleep(1000000);
                ctx->IDinOPD = read_sharedmem_image(ctx->INSTREAMNAMEOPD, dcimg, dcnimg);
            }
            ctx->IDinAMP = read_sharedmem_image(ctx->INSTREAMNAMEAMP, dcimg, dcnimg);
            while (ctx->IDinAMP == -1)
            {
                usleep(1000000);
                ctx->IDinAMP = read_sharedmem_image(ctx->INSTREAMNAMEAMP, dcimg, dcnimg);
            }
            break;

        case 1:
            load_fits(ctx->INFITSFILENAMEOPD, "inOPD", 1, &ctx->IDinOPD);
            load_fits(ctx->INFITSFILENAMEAMP, "inAMP", 1, &ctx->IDinAMP);
            break;

        default:
            printf("INMODE value %d not valid\n", ctx->INMODE);
            exit(0);
            break;
    }

    sizearray[0]  = (uint32_t) ctx->OUTARRAYSIZE;
    sizearray[1]  = (uint32_t) ctx->OUTARRAYSIZE;
    ctx->wfinsize = dcimg[ctx->IDinOPD].md[0].size[0];
    if (ctx->OUTMODE == 0)
    {
        create_image_ID(ctx->OUTSTREAMNAME, 2, sizearray, _DATATYPE_FLOAT, 1, 0, 0, &ctx->IDout);
        COREMOD_MEMORY_image_set_createsem(ctx->OUTSTREAMNAME, 10);
        sizearray[0] = (uint32_t) ctx->ARRAYSIZE;
        sizearray[1] = (uint32_t) ctx->ARRAYSIZE;
        create_image_ID(ctx->OUTINSTSTREAMNAME,
                        2,
                        sizearray,
                        _DATATYPE_FLOAT,
                        1,
                        0,
                        0,
                        &ctx->IDoutinst);
        COREMOD_MEMORY_image_set_createsem(ctx->OUTINSTSTREAMNAME, 10);
    }
    else
    {
        create_image_ID(ctx->OUTSTREAMNAME, 2, sizearray, _DATATYPE_FLOAT, 0, 0, 0, &ctx->IDout);
        create_image_ID(ctx->OUTINSTSTREAMNAME,
                        2,
                        sizearray,
                        _DATATYPE_FLOAT,
                        0,
                        0,
                        0,
                        &ctx->IDoutinst);
    }

    create_2Dimage_ID("pupamp", ctx->ARRAYSIZE, ctx->ARRAYSIZE, &ctx->IDpupamp);
    create_2Dimage_ID("puppha", ctx->ARRAYSIZE, ctx->ARRAYSIZE, &ctx->IDpuppha);
    create_2Dimage_ID("focamp", ctx->ARRAYSIZE, ctx->ARRAYSIZE, &ctx->IDfocamp);
    create_2Dimage_ID("focpha", ctx->ARRAYSIZE, ctx->ARRAYSIZE, &ctx->IDfocpha);
    create_2Dimage_ID("pyrpupi", ctx->ARRAYSIZE, ctx->ARRAYSIZE, &ctx->IDpyrpupi);
    ctx->offset = (ctx->ARRAYSIZE - ctx->wfinsize) / 2;
}

/**
 * pyrwfs_step_propagate - Propagate input WF through one modulation point and update WFSinst
 * @ctx:    Pyramid WFS daemon context
 * @pmodpt: Current modulation point index
 */
static void pyrwfs_step_propagate(
    struct pyrwfs_daemon_ctx *ctx,
    long                      pmodpt)
{
    long     arr  = ctx->ARRAYSIZE;
    uint64_t arr2 = (uint64_t) arr * arr;
    imageID  IDfoca;
    imageID  IDfocp;
    imageID  IDpupa;

    for (long ii = 0; ii < ctx->wfinsize; ii++)
    {
        for (long jj = 0; jj < ctx->wfinsize; jj++)
        {
            long dst = (jj + ctx->offset) * arr + (ii + ctx->offset);
            long src = jj * ctx->wfinsize + ii;

            dcimg[ctx->IDpuppha].array.F[dst] =
                (float) (2.0 * M_PI * (1.0e-6 * dcimg[ctx->IDinOPD].array.F[src]) / ctx->LAMBDA);
            dcimg[ctx->IDpupamp].array.F[dst] = dcimg[ctx->IDinAMP].array.F[src];
        }
    }

    mk_complex_from_amph("pupamp", "puppha", "wfc", 0);
    permut("wfc");
    do2dfft("wfc", "imc");
    delete_image_ID("wfc", DELETE_IMAGE_ERRMODE_WARNING);
    permut("imc");
    mk_amph_from_complex("imc", "ima", "imp", 0);
    delete_image_ID("imc", DELETE_IMAGE_ERRMODE_WARNING);
    IDfoca = image_ID("ima", dcimg, dcnimg);
    IDfocp = image_ID("imp", dcimg, dcnimg);

    for (uint64_t ii = 0; ii < arr2; ii++)
    {
        dcimg[ctx->IDfocamp].array.F[ii] =
            dcimg[IDfoca].array.F[ii] * dcimg[ctx->IDpyr_amp].array.F[pmodpt * arr2 + ii];
        dcimg[ctx->IDfocpha].array.F[ii] =
            dcimg[IDfocp].array.F[ii] + dcimg[ctx->IDpyr_pha].array.F[pmodpt * arr2 + ii];
    }

    mk_complex_from_amph("focamp", "focpha", "focc", 0);
    permut("focc");
    do2dfft("focc", "pupc");
    delete_image_ID("focc", DELETE_IMAGE_ERRMODE_WARNING);
    permut("pupc");
    mk_amph_from_complex("pupc", "pupa", "pupp", 0);
    IDpupa = image_ID("pupa", dcimg, dcnimg);
    delete_image_ID("pupc", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pupp", DELETE_IMAGE_ERRMODE_WARNING);

    dcimg[ctx->IDoutinst].md[0].write = 1;
    for (uint64_t ii = 0; ii < arr2; ii++)
    {
        dcimg[ctx->IDoutinst].array.F[ii] = dcimg[IDpupa].array.F[ii] * dcimg[IDpupa].array.F[ii];
    }
    dcimg[ctx->IDoutinst].md[0].cnt0++;
    dcimg[ctx->IDoutinst].md[0].write = 0;
    if (ctx->OUTMODE == 0)
    {
        COREMOD_MEMORY_image_set_sempost(ctx->OUTINSTSTREAMNAME, -1);
    }
    delete_image_ID("pupa", DELETE_IMAGE_ERRMODE_WARNING);

    for (uint64_t ii = 0; ii < arr2; ii++)
    {
        dcimg[ctx->IDpyrpupi].array.F[ii] += dcimg[ctx->IDoutinst].array.F[ii];
    }
}

/**
 * pyrwfs_write_readout - Bin modulated detector image, inject noise, and publish WFS frame
 * @ctx: Pyramid WFS daemon context
 */
static void pyrwfs_write_readout(struct pyrwfs_daemon_ctx *ctx)
{
    uint64_t arr2      = (uint64_t) ctx->ARRAYSIZE * ctx->ARRAYSIZE;
    double   TotalFlux = 0.0;
    long     outoffset = (ctx->ARRAYSIZE - ctx->OUTBINFACT * ctx->OUTARRAYSIZE) / 2;
    char     fname[STRINGMAXLEN_FULLFILENAME];

    for (uint64_t ii = 0; ii < arr2; ii++)
    {
        TotalFlux += dcimg[ctx->IDpyrpupi].array.F[ii];
    }
    for (uint64_t ii = 0; ii < arr2; ii++)
    {
        dcimg[ctx->IDpyrpupi].array.F[ii] *= (float) (ctx->WFSFLUX / TotalFlux);
    }

    dcimg[ctx->IDout].md[0].write = 1;
    for (long ii = 0; ii < ctx->OUTARRAYSIZE; ii++)
    {
        for (long jj = 0; jj < ctx->OUTARRAYSIZE; jj++)
        {
            double tmpval = 0.0;

            for (long i = 0; i < ctx->OUTBINFACT; i++)
            {
                for (long j = 0; j < ctx->OUTBINFACT; j++)
                {
                    long ii1 = ii * ctx->OUTBINFACT + i + outoffset;
                    long jj1 = jj * ctx->OUTBINFACT + j + outoffset;

                    tmpval += dcimg[ctx->IDpyrpupi].array.F[jj1 * ctx->ARRAYSIZE + ii1];
                }
            }
            dcimg[ctx->IDout].array.F[jj * ctx->OUTARRAYSIZE + ii] =
                (float) (fast_poisson(tmpval) + ctx->WFSCAMRON * gauss());
        }
    }
    dcimg[ctx->IDout].md[0].cnt0++;
    dcimg[ctx->IDout].md[0].write = 0;

    if (ctx->OUTMODE == 0)
    {
        COREMOD_MEMORY_image_set_sempost(ctx->OUTSTREAMNAME, -1);
    }
    if (ctx->OUTMODE == 1)
    {
        EXECUTE_SYSTEM_COMMAND("rm %s", ctx->OUTFITSFILENAME);
        WRITE_FILENAME(fname, "%s", ctx->OUTFITSFILENAME);
        save_fits(ctx->OUTSTREAMNAME, fname);
        EXECUTE_SYSTEM_COMMAND("touch %s", ctx->OUTTRIGGERFILE);
    }
    for (uint64_t ii = 0; ii < arr2; ii++)
    {
        dcimg[ctx->IDpyrpupi].array.F[ii] = 0.0f;
    }
}

/**
 * pyrwfs_wait_trigger - Wait for next input wavefront frame semaphore or trigger file
 * @ctx: Pyramid WFS daemon context
 */
static void pyrwfs_wait_trigger(const struct pyrwfs_daemon_ctx *ctx)
{
    int OKf;

    switch (ctx->INMODE)
    {
        case 0:
            COREMOD_MEMORY_image_set_semwait(ctx->INSTREAMNAMEOPD, ctx->INSEMCHAN);
            break;

        case 1:
            OKf = 0;
            while (OKf == 0)
            {
                usleep(10000);
                if (file_exists(ctx->INTRIGGERFILE) == 1)
                {
                    EXECUTE_SYSTEM_COMMAND("rm %s", ctx->INTRIGGERFILE);
                    OKf = 1;
                }
            }
            break;

        default:
            printf("INMODE value %d not valid\n", ctx->INMODE);
            exit(0);
            break;
    }
}

/**
 * AOsystSim_PyrWFS - Run real-time shared-memory Pyramid WFS simulation loop
 * @CONF_FNAME: Pyramid WFS daemon configuration file path
 *
 * Return: 0 on completion.
 */
int AOsystSim_PyrWFS(const char *CONF_FNAME)
{
    struct pyrwfs_daemon_ctx ctx    = { 0 };
    long                     pmodpt = 0;
    long                     k      = 0;
    long                     kmax   = 100000000;

    aosim_ensure_wdir();
    printf("AOsystSim PyrWFS...\n");

    pyrwfs_read_conf(CONF_FNAME, &ctx);
    pyrwfs_init_masks(&ctx);
    pyrwfs_init_streams(&ctx);

    while (k < kmax)
    {
        pyrwfs_step_propagate(&ctx, pmodpt);
        pmodpt++;
        if (pmodpt == ctx.PYRMODNBPT)
        {
            pmodpt = 0;
            pyrwfs_write_readout(&ctx);
        }
        pyrwfs_wait_trigger(&ctx);
    }

    return 0;
}
