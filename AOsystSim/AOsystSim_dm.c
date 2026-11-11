/**
 * @file    AOsystSim_dm.c
 * @brief   Real-time deformable mirror temporal-response and 2D OPD surface daemon
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "image_filter/image_filter.h"

#include "AOsystSim/AOsystSim.h"
#include "AOsystSim_internal.h"

/**
 * @brief Runtime state and sparse influence-function tables for the DM daemon.
 */
struct aosim_dm_ctx
{
    int     INMODE;
    char    INSTREAMNAMEDM[200];
    char    INFITSFILENAMEDM[200];
    char    INTRIGSTREAMNAME[200];
    int     INTRIGSEMCHAN;
    char    INTRIGGERFILE[200];
    int     OUTMODE;
    char    OUTSTREAMNAMEDM[200];
    char    OUTFITSFILENAMEDM[200];
    char    OUTTRIGGERFILE[200];
    char    OUTSTREAMNAMEOPD[200];
    char    OUTSTREAMNAMEAMP[200];
    char    OUTFITSFILENAMEOPD[200];
    char    OUTFITSFILENAMEAMP[200];
    long    ARRAYSIZE;
    float   DMRAD;
    float   DMDT;
    long    NBTSAMPLES;
    float   DMLAGSTART;
    float   DMTIMECST;
    long    DMsize;
    imageID IDinDM;
    imageID IDinTRIG;
    imageID IDout;
    imageID IDout_tmp;
    imageID IDdmdispC;
    imageID IDdmdispC1;
    long    NBdmifcarray;
    float  *dmifcarray_value;
    long   *dmifcarray_act;
    long   *dmifcarray_iijj;
};

/**
 * aosim_dm_read_conf - Read DM daemon configuration parameters from disk
 * @conf_fname: Configuration file path
 * @ctx:        DM daemon runtime context
 */
static void aosim_dm_read_conf(
    const char          *conf_fname,
    struct aosim_dm_ctx *ctx)
{
    FILE *fp;
    char  fname[STRINGMAXLEN_FULLFILENAME];

    if ((fp = fopen(conf_fname, "r")) == NULL)
    {
        sprintf(fname, "%s.default", conf_fname);
        printf("configuration file %s not found. Creating default template as %s\n",
               conf_fname,
               fname);
        AOsystSim_DM_mkCONF(fname);
        exit(0);
    }
    fclose(fp);

    aosim_ensure_wdir();

    ctx->INMODE = read_config_parameter_long(conf_fname, "INMODE");
    read_config_parameter(conf_fname, "INSTREAMNAMEDM", ctx->INSTREAMNAMEDM);
    read_config_parameter(conf_fname, "INFITSFILENAMEDM", ctx->INFITSFILENAMEDM);
    read_config_parameter(conf_fname, "INTRIGSTREAMNAME", ctx->INTRIGSTREAMNAME);
    ctx->INTRIGSEMCHAN = read_config_parameter_long(conf_fname, "INTRIGSEMCHAN");
    read_config_parameter(conf_fname, "INTRIGGERFILE", ctx->INTRIGGERFILE);

    ctx->OUTMODE = read_config_parameter_long(conf_fname, "OUTMODE");
    read_config_parameter(conf_fname, "OUTSTREAMNAMEDM", ctx->OUTSTREAMNAMEDM);
    read_config_parameter(conf_fname, "OUTFITSFILENAMEDM", ctx->OUTFITSFILENAMEDM);
    read_config_parameter(conf_fname, "OUTTRIGGERFILE", ctx->OUTTRIGGERFILE);
    read_config_parameter(conf_fname, "OUTSTREAMNAMEOPD", ctx->OUTSTREAMNAMEOPD);
    read_config_parameter(conf_fname, "OUTSTREAMNAMEAMP", ctx->OUTSTREAMNAMEAMP);
    read_config_parameter(conf_fname, "OUTFITSFILENAMEOPD", ctx->OUTFITSFILENAMEOPD);
    read_config_parameter(conf_fname, "OUTFITSFILENAMEAMP", ctx->OUTFITSFILENAMEAMP);

    ctx->ARRAYSIZE  = read_config_parameter_long(conf_fname, "ARRAYSIZE");
    ctx->DMRAD      = read_config_parameter_float(conf_fname, "DMRAD");
    ctx->DMDT       = read_config_parameter_float(conf_fname, "DMDT");
    ctx->NBTSAMPLES = read_config_parameter_long(conf_fname, "NBTSAMPLES");
    ctx->DMLAGSTART = read_config_parameter_float(conf_fname, "DMLAGSTART");
    ctx->DMTIMECST  = read_config_parameter_float(conf_fname, "DMTIMECST");
}

/**
 * aosim_dm_init_streams - Attach input command/trigger streams and allocate temporal cubes
 * @ctx: DM daemon runtime context
 */
static void aosim_dm_init_streams(struct aosim_dm_ctx *ctx)
{
    uint32_t sizearray[2] = { (uint32_t) ctx->ARRAYSIZE, (uint32_t) ctx->ARRAYSIZE };
    int      shared       = (ctx->OUTMODE == 0) ? 1 : 0;

    switch (ctx->INMODE)
    {
        case 0:
            ctx->IDinDM = read_sharedmem_image(ctx->INSTREAMNAMEDM, dcimg, dcnimg);
            while (ctx->IDinDM == -1)
            {
                usleep(100000);
                ctx->IDinDM = read_sharedmem_image(ctx->INSTREAMNAMEDM, dcimg, dcnimg);
            }

            printf("READING INPUT TRIGGER STREAM ...\n");
            ctx->IDinTRIG = read_sharedmem_image(ctx->INTRIGSTREAMNAME, dcimg, dcnimg);
            printf("  IDinTRIG = %ld\n", ctx->IDinTRIG);
            while (ctx->IDinTRIG == -1)
            {
                usleep(100000);
                ctx->IDinTRIG = read_sharedmem_image(ctx->INTRIGSTREAMNAME, dcimg, dcnimg);
                printf("  IDinTRIG = %ld\n", ctx->IDinTRIG);
            }
            break;

        case 1:
            break;

        default:
            printf("INMODE value %d not valid\n", ctx->INMODE);
            exit(0);
            break;
    }

    ctx->DMsize = dcimg[ctx->IDinDM].md[0].size[0];
    create_image_ID(ctx->OUTSTREAMNAMEDM, 2, sizearray, _DATATYPE_FLOAT, shared, 0, 0, &ctx->IDout);
    if (shared)
    {
        COREMOD_MEMORY_image_set_createsem(ctx->OUTSTREAMNAMEDM, 10);
    }
    create_image_ID("tmpDMshape", 2, sizearray, _DATATYPE_FLOAT, 0, 0, 0, &ctx->IDout_tmp);
    create_3Dimage_ID("dmdispC", ctx->DMsize, ctx->DMsize, ctx->NBTSAMPLES, &ctx->IDdmdispC);
    create_3Dimage_ID("dmdispC1", ctx->DMsize, ctx->DMsize, ctx->NBTSAMPLES, &ctx->IDdmdispC1);
}

/**
 * aosim_dm_pack_sparse - Compress 3D influence-function cube into sparse coordinate arrays
 * @ctx:        DM daemon runtime context
 * @IDifc:      3D influence-function cube image ID
 * @dmif_limit: Sparsity cutoff threshold
 */
static void aosim_dm_pack_sparse(
    struct aosim_dm_ctx *ctx,
    imageID              IDifc,
    float                dmif_limit)
{
    long arr = ctx->ARRAYSIZE;

    ctx->dmifcarray_value = (float *) malloc(sizeof(float) * ctx->NBdmifcarray);
    ctx->dmifcarray_act   = (long *) malloc(sizeof(long) * ctx->NBdmifcarray);
    ctx->dmifcarray_iijj  = (long *) malloc(sizeof(long) * ctx->NBdmifcarray);
    if ((ctx->dmifcarray_value == NULL) || (ctx->dmifcarray_act == NULL) ||
        (ctx->dmifcarray_iijj == NULL))
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    ctx->NBdmifcarray = 0;
    for (long mx = 0; mx < ctx->DMsize; mx++)
    {
        for (long my = 0; my < ctx->DMsize; my++)
        {
            long act = my * ctx->DMsize + mx;

            for (long ii = 0; ii < arr; ii++)
            {
                for (long jj = 0; jj < arr; jj++)
                {
                    float val = dcimg[IDifc].array.F[act * arr * arr + jj * arr + ii];

                    if (fabsf(val) > dmif_limit)
                    {
                        ctx->dmifcarray_value[ctx->NBdmifcarray] = val;
                        ctx->dmifcarray_act[ctx->NBdmifcarray]   = act;
                        ctx->dmifcarray_iijj[ctx->NBdmifcarray]  = jj * arr + ii;
                        ctx->NBdmifcarray++;
                    }
                }
            }
        }
    }
}

/**
 * aosim_dm_build_ifc - Construct Gaussian-smoothed DM influence cube and sparse index
 * @ctx: DM daemon runtime context
 */
static void aosim_dm_build_ifc(struct aosim_dm_ctx *ctx)
{
    long     arr        = ctx->ARRAYSIZE;
    long     DMsize     = ctx->DMsize;
    double   dmifscale  = 4.0;
    float    dmif_limit = 0.0001f;
    float    actpix     = (float) (dmifscale * (2.0 * ctx->DMRAD / DMsize));
    float    sig        = 0.5f * actpix;
    uint32_t imsize[3]  = { (uint32_t) arr, (uint32_t) arr, (uint32_t) (DMsize * DMsize) };
    imageID  IDif;
    imageID  IDifc;

    create_image_ID("dmif0", 2, imsize, _DATATYPE_FLOAT, 0, 0, 0, &IDif);
    printf("DM size = %ld x %ld actuators\n", DMsize, DMsize);
    printf("actuator pix size = %f pix\n", actpix);
    list_image_ID();

    for (long ii = 0; ii < arr; ii++)
    {
        for (long jj = 0; jj < arr; jj++)
        {
            float x = (1.0f * ii - 0.5f * arr) / actpix;
            float y = (1.0f * jj - 0.5f * arr) / actpix;

            dcimg[IDif].array.F[jj * arr + ii] =
                ((fabsf(x) < 0.5f) && (fabsf(y) < 0.5f)) ? 1.0f : 0.0f;
        }
    }

    printf("convolve\n");
    fflush(stdout);
    save_fits("dmif0", "AOsystSim_wdir/dmif0.fits");
    printf("gauss filter   %lf %ld\n", (double) sig, (long) (2.0 * sig));
    fflush(stdout);
    gauss_filter("dmif0", "dmif", sig, (long) (2.0 * sig));
    list_image_ID();
    delete_image_ID("dmif0", DELETE_IMAGE_ERRMODE_WARNING);
    IDif = image_ID("dmif", dcimg, dcnimg);
    save_fits("dmif", "AOsystSim_wdir/dmif.fits");

    create_image_ID("dmifc", 3, imsize, _DATATYPE_FLOAT, 0, 0, 0, &IDifc);
    printf("\ndmifc = %ld   %ld %ld %ld    %ld %ld\n",
           IDifc,
           (long) dcimg[IDifc].md[0].size[0],
           (long) dcimg[IDifc].md[0].size[1],
           (long) dcimg[IDifc].md[0].size[2],
           arr,
           arr);

    ctx->NBdmifcarray = 0;
    for (long mx = 0; mx < DMsize; mx++)
    {
        for (long my = 0; my < DMsize; my++)
        {
            float ii0f = 0.5f * arr + 2.0f * (mx - DMsize / 2) / DMsize * ctx->DMRAD;
            float jj0f = 0.5f * arr + 2.0f * (my - DMsize / 2) / DMsize * ctx->DMRAD;

            printf("\r actuator %2ld %2ld    ", mx, my);
            fflush(stdout);
            for (long ii = 0; ii < arr; ii++)
            {
                for (long jj = 0; jj < arr; jj++)
                {
                    float rxif = (1.0f * ii - ii0f) * (float) dmifscale + arr / 2;
                    float ryif = (1.0f * jj - jj0f) * (float) dmifscale + arr / 2;
                    long  rxi  = (long) rxif;
                    long  ryi  = (long) ryif;

                    if ((rxi > 0) && (rxi < arr - 1) && (ryi > 0) && (ryi < arr - 1))
                    {
                        float u   = rxif - rxi;
                        float t   = ryif - ryi;
                        float v00 = dcimg[IDif].array.F[ryi * arr + rxi];
                        float v01 = dcimg[IDif].array.F[ryi * arr + rxi + 1];
                        float v10 = dcimg[IDif].array.F[(ryi + 1) * arr + rxi];
                        float v11 = dcimg[IDif].array.F[(ryi + 1) * arr + rxi + 1];
                        float val = (1.0f - u) * (1.0f - t) * v00 + (1.0f - u) * t * v10 +
                                    u * (1.0f - t) * v01 + u * t * v11;

                        dcimg[IDifc].array.F[(my * DMsize + mx) * arr * arr + jj * arr + ii] = val;
                        if (fabsf(val) > dmif_limit)
                        {
                            ctx->NBdmifcarray++;
                        }
                    }
                }
            }
        }
    } // for (long mx)
    save_fits("dmifc", "AOsystSim_wdir/dmifc.fits");
    printf("\n");

    aosim_dm_pack_sparse(ctx, IDifc, dmif_limit);
}

/**
 * aosim_dm_step - Advance actuator temporal response, compute 2D DM OPD, and wait on trigger
 * @ctx: DM daemon runtime context
 */
static void aosim_dm_step(struct aosim_dm_ctx *ctx)
{
    uint64_t dm2      = (uint64_t) ctx->DMsize * ctx->DMsize;
    long     nelement = ctx->ARRAYSIZE * ctx->ARRAYSIZE;
    long     kkstart  = (long) (ctx->DMLAGSTART / ctx->DMDT);
    float    alpha    = (float) pow(0.5, ctx->DMDT / ctx->DMTIMECST);
    char     fname[STRINGMAXLEN_FULLFILENAME];
    int      OKf;

    memcpy(dcimg[ctx->IDdmdispC1].array.F,
           dcimg[ctx->IDdmdispC].array.F,
           sizeof(float) * dm2 * ctx->NBTSAMPLES);
    memcpy(dcimg[ctx->IDdmdispC].array.F,
           dcimg[ctx->IDdmdispC1].array.F + dm2,
           sizeof(float) * dm2 * (ctx->NBTSAMPLES - 1));

    for (long kk = kkstart; kk < ctx->NBTSAMPLES; kk++)
    {
        for (uint64_t ii = 0; ii < dm2; ii++)
        {
            float cmd = dcimg[ctx->IDinDM].array.F[ii];
            float cur = dcimg[ctx->IDdmdispC].array.F[(uint64_t) kk * dm2 + ii];

            dcimg[ctx->IDdmdispC].array.F[(uint64_t) kk * dm2 + ii] = cmd - alpha * (cmd - cur);
        }
    }

    memset(dcimg[ctx->IDout_tmp].array.F, '\0', sizeof(float) * nelement);
    for (long kk = 0; kk < ctx->NBdmifcarray; kk++)
    {
        dcimg[ctx->IDout_tmp].array.F[ctx->dmifcarray_iijj[kk]] +=
            dcimg[ctx->IDdmdispC].array.F[ctx->dmifcarray_act[kk]] * ctx->dmifcarray_value[kk];
    }

    dcimg[ctx->IDout].md[0].write = 1;
    memcpy(dcimg[ctx->IDout].array.F, dcimg[ctx->IDout_tmp].array.F, sizeof(float) * nelement);
    dcimg[ctx->IDout].md[0].cnt0++;
    dcimg[ctx->IDout].md[0].write = 0;

    if (ctx->OUTMODE == 0)
    {
        COREMOD_MEMORY_image_set_sempost(ctx->OUTSTREAMNAMEDM, -1);
    }
    if (ctx->OUTMODE == 1)
    {
        EXECUTE_SYSTEM_COMMAND("rm %s", ctx->OUTFITSFILENAMEDM);
        WRITE_FILENAME(fname, "%s", ctx->OUTFITSFILENAMEDM);
        save_fits(ctx->OUTSTREAMNAMEDM, fname);
        EXECUTE_SYSTEM_COMMAND("touch %s", ctx->OUTTRIGGERFILE);
    }

    switch (ctx->INMODE)
    {
        case 0:
            COREMOD_MEMORY_image_set_semwait(ctx->INTRIGSTREAMNAME, ctx->INTRIGSEMCHAN);
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
 * AOsystSim_DM - Run real-time shared-memory deformable mirror OPD surface daemon
 * @CONF_FNAME: Deformable mirror daemon configuration file path
 *
 * Return: 0 on completion.
 */
int AOsystSim_DM(const char *CONF_FNAME)
{
    struct aosim_dm_ctx ctx  = { 0 };
    long                kmax = 100000000;

    printf("AOsystSim DM...\n");

    aosim_dm_read_conf(CONF_FNAME, &ctx);
    aosim_dm_init_streams(&ctx);
    aosim_dm_build_ifc(&ctx);

    for (long k = 0; k < kmax; k++)
    {
        aosim_dm_step(&ctx);
    }

    free(ctx.dmifcarray_value);
    free(ctx.dmifcarray_act);
    free(ctx.dmifcarray_iijj);

    return 0;
}
