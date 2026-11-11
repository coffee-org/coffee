/**
 * @file    AOsystSim_lowfs_init.c
 * @brief   Configuration file generation, stream attachment, and buffer setup for coroLOWFS
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_tools/COREMOD_tools.h"

#include "AOsystSim/AOsystSim.h"
#include "AOsystSim_internal.h"

/**
 * AOsystSim_coroLOWFS_mkCONF - Write template configuration file for coronagraph + LOWFS daemon
 * @fname: Output configuration file path
 *
 * Return: 0 on completion.
 */
int AOsystSim_coroLOWFS_mkCONF(const char *fname)
{
    FILE *fp = fopen(fname, "w");

    fprintf(fp, "\n");
    fprintf(fp, "# ============== AOsim: coro LOWFS process =============================\n");
    fprintf(fp, "# LOWFS process is triggered by WF OPD\n\n");
    fprintf(fp, "# ============== INPUT & TRIGGER (OPD unit = um) ===============\n");
    fprintf(fp, "INMODE                   0                 "
                "# 0:stream, 1:file system (FITS)\n");
    fprintf(fp, "INAMPSTREAMNAME          wf1amp            # input AMP stream\n");
    fprintf(fp, "INOPDSTREAMNAME          wf1opd            # input OPD stream\n");
    fprintf(fp, "INAMPFITSFILENAME        wf1amp.fits       # input FITS file name\n");
    fprintf(fp, "INOPDFITSFILENAME        wf1opd.fits       # input FITS file name\n");
    fprintf(fp, "INTRIGSTREAMNAME         wf1opd            # input trigger stream\n");
    fprintf(fp, "INTRIGSEMCHAN            4                 "
                "# input semaphore channel (using OPD input)\n");
    fprintf(fp, "INTRIGGERFILE            inwf.txt          # input trigger file\n");
    fprintf(fp, "INPHYSTIME               aosim_phystime    # physical time\n\n");
    fprintf(fp, "# ============== OUTPUT TYPE ===================================\n");
    fprintf(fp, "OUTMODE                  0                       "
                "# 0: stream, 1: file system\n");
    fprintf(fp, "OUTLOWFSARRAYSIZE        16                      "
                "# output LOWFS array size, pix\n");
    fprintf(fp, "OUTLOWFSSTREAMNAME       aosim_imcamLOWFS        # output WF stream name \n");
    fprintf(fp, "OUTLOWFSFITSFILENAME     aosim_imcamLOWFS.fits   "
                "# FITS file name output [um]\n");
    fprintf(fp, "OUTTRIGGERFILE           outLOWFS.txt            # output trigger file\n\n");
    fprintf(fp, "# ============== GEOMETRY, OPTICAL DESIGN =============================\n");
    fprintf(fp, "LAMBDAum                 1.65              # wavelength [um]\n");
    fprintf(fp, "ARRAYSIZE                256               "
                "# array size, pix, used for internal computations\n");
    fprintf(fp, "COROFPMAMP               corofpmamp        "
                "# coronagraph focal plane mask amplitude\n");
    fprintf(fp, "COROFPMPHA               corofpmpha        "
                "# coronagraph focal plane mask amplitude\n");
    fprintf(fp, "LYOTSTOPREFLAMP          lyotstopreflamp   "
                "# LYOT stop reflectivity map (amp)\n");
    fprintf(fp, "LYOTSTOPREFLPHA          lyotstopreflpha   "
                "# LYOT stop reflectivity map (pha)\n");
    fprintf(fp, "LYOTSTOPTRANSMAMP        lyotstoptransmamp "
                "# LYOT stop transmission map (amp)\n");
    fprintf(fp, "LYOTSTOPTRANSMPHA        lyotstoptransmpha "
                "# LYOT stop transmission map (pha)\n\n");
    fprintf(fp, "# ============== LOWFS =============================\n");
    fprintf(fp, "LOWFSOPDMAP              lowfsopdmap       # LOWFS defocus map [um]\n");
    fprintf(fp, "LOWFSCAMETIME            0.01              "
                "# LOWFS camera exposure time [s]\n\n");
    fclose(fp);

    return 0;
}

/**
 * aosim_lowfs_read_conf - Read coronagraph and LOWFS daemon configuration parameters
 * @conf_fname: Configuration file path
 * @ctx:        Coronagraph/LOWFS runtime context
 */
void aosim_lowfs_read_conf(
    const char      *conf_fname,
    AOSIM_LOWFS_CTX *ctx)
{
    FILE *fp;
    char  fname[STRINGMAXLEN_FULLFILENAME];

    if ((fp = fopen(conf_fname, "r")) == NULL)
    {
        sprintf(fname, "%s.default", conf_fname);
        printf("configuration file %s not found. Creating default template as %s\n",
               conf_fname,
               fname);
        AOsystSim_coroLOWFS_mkCONF(fname);
        exit(0);
    }
    fclose(fp);

    ctx->INMODE = read_config_parameter_long(conf_fname, "INMODE");
    read_config_parameter(conf_fname, "INAMPSTREAMNAME", ctx->INAMPSTREAMNAME);
    read_config_parameter(conf_fname, "INOPDSTREAMNAME", ctx->INOPDSTREAMNAME);
    read_config_parameter(conf_fname, "INAMPFITSFILENAME", ctx->INAMPFITSFILENAME);
    read_config_parameter(conf_fname, "INOPDFITSFILENAME", ctx->INOPDFITSFILENAME);
    read_config_parameter(conf_fname, "INTRIGSTREAMNAME", ctx->INTRIGSTREAMNAME);
    ctx->INTRIGSEMCHAN = read_config_parameter_long(conf_fname, "INTRIGSEMCHAN");
    read_config_parameter(conf_fname, "INTRIGGERFILE", ctx->INTRIGGERFILE);
    read_config_parameter(conf_fname, "INPHYSTIME", ctx->INPHYSTIME);

    ctx->OUTMODE           = read_config_parameter_long(conf_fname, "OUTMODE");
    ctx->OUTLOWFSARRAYSIZE = (uint32_t) read_config_parameter_long(conf_fname, "OUTLOWFSARRAYSIZE");
    ctx->OUTLOWFSARRAYSIZE2 = (uint64_t) ctx->OUTLOWFSARRAYSIZE * ctx->OUTLOWFSARRAYSIZE;
    read_config_parameter(conf_fname, "OUTLOWFSSTREAMNAME", ctx->OUTLOWFSSTREAMNAME);
    read_config_parameter(conf_fname, "OUTLOWFSFITSFILENAME", ctx->OUTLOWFSFITSFILENAME);
    read_config_parameter(conf_fname, "OUTTRIGGERFILE", ctx->OUTTRIGGERFILE);

    ctx->LAMBDA     = 1.0e-6 * read_config_parameter_float(conf_fname, "LAMBDAum");
    ctx->ARRAYSIZE  = (uint32_t) read_config_parameter_long(conf_fname, "ARRAYSIZE");
    ctx->ARRAYSIZE2 = (uint64_t) ctx->ARRAYSIZE * ctx->ARRAYSIZE;
    read_config_parameter(conf_fname, "COROFPMAMP", ctx->COROFPMAMP);
    read_config_parameter(conf_fname, "COROFPMPHA", ctx->COROFPMPHA);
    read_config_parameter(conf_fname, "LYOTSTOPREFLAMP", ctx->LYOTSTOPREFLAMP);
    read_config_parameter(conf_fname, "LYOTSTOPREFLPHA", ctx->LYOTSTOPREFLPHA);
    read_config_parameter(conf_fname, "LYOTSTOPTRANSMAMP", ctx->LYOTSTOPTRANSMAMP);
    read_config_parameter(conf_fname, "LYOTSTOPTRANSMPHA", ctx->LYOTSTOPTRANSMPHA);

    read_config_parameter(conf_fname, "LOWFSOPDMAP", ctx->LOWFSOPDMAP);
    ctx->LOWFSCAMETIME = read_config_parameter_float(conf_fname, "LOWFSCAMETIME");
}

/**
 * lowfs_require_image - Lookup preloaded optical mask image or exit on error
 * @name: Required image name in memory
 *
 * Return: Valid image ID.
 */
static imageID lowfs_require_image(const char *name)
{
    imageID ID = image_ID(name, dcimg, dcnimg);

    if (ID == -1)
    {
        printf("ERROR: image %s not loaded\n", name);
        exit(0);
    }
    return ID;
}

/**
 * lowfs_wait_shm - Poll until a shared-memory stream becomes available
 * @label: Diagnostic stream label
 * @name:  Shared-memory stream name
 *
 * Return: Valid attached image ID.
 */
static imageID lowfs_wait_shm(
    const char *label,
    const char *name)
{
    imageID ID = read_sharedmem_image(name, dcimg, dcnimg);

    while (ID == -1)
    {
        usleep(100000);
        ID = read_sharedmem_image(name, dcimg, dcnimg);
        printf("  %s stream \"%s\" = %ld\n", label, name, ID);
    }
    return ID;
}

/**
 * aosim_lowfs_connect_inputs - Attach input wavefront streams and validate preloaded optical masks
 * @ctx: Coronagraph/LOWFS runtime context
 */
void aosim_lowfs_connect_inputs(AOSIM_LOWFS_CTX *ctx)
{
    printf("READING INPHYSTIME STREAM \"%s\" ...\n", ctx->INPHYSTIME);
    ctx->IDphystime = read_sharedmem_image(ctx->INPHYSTIME, dcimg, dcnimg);
    printf("  IDphystime = %ld\n", ctx->IDphystime);
    while (ctx->IDphystime == -1)
    {
        usleep(100000);
        ctx->IDphystime = read_sharedmem_image(ctx->INPHYSTIME, dcimg, dcnimg);
        printf("  IDphystime = %ld\n", ctx->IDphystime);
    }

    switch (ctx->INMODE)
    {
        case 0:
            ctx->IDinOPD  = lowfs_wait_shm("INOPDSTREAMNAME", ctx->INOPDSTREAMNAME);
            ctx->IDinAMP  = lowfs_wait_shm("INAMPSTREAMNAME", ctx->INAMPSTREAMNAME);
            ctx->IDinTRIG = lowfs_wait_shm("INTRIGSTREAMNAME", ctx->INTRIGSTREAMNAME);
            break;

        case 1:
            break;

        default:
            printf("INMODE value %d not valid\n", ctx->INMODE);
            exit(0);
            break;
    }

    ctx->IDfpmamp    = lowfs_require_image(ctx->COROFPMAMP);
    ctx->IDfpmpha    = lowfs_require_image(ctx->COROFPMPHA);
    ctx->ID_LS_ramp  = lowfs_require_image(ctx->LYOTSTOPREFLAMP);
    ctx->ID_LS_rpha  = lowfs_require_image(ctx->LYOTSTOPREFLPHA);
    ctx->ID_LS_tamp  = lowfs_require_image(ctx->LYOTSTOPTRANSMAMP);
    ctx->ID_LS_tpha  = lowfs_require_image(ctx->LYOTSTOPTRANSMPHA);
    ctx->ID_lowfsopd = lowfs_require_image(ctx->LOWFSOPDMAP);
}

/**
 * aosim_lowfs_init_buffers - Allocate output LOWFS streams and internal optical-plane buffers
 * @ctx: Coronagraph/LOWFS runtime context
 */
void aosim_lowfs_init_buffers(AOSIM_LOWFS_CTX *ctx)
{
    uint32_t sz_out[2]  = { ctx->OUTLOWFSARRAYSIZE, ctx->OUTLOWFSARRAYSIZE };
    uint32_t sz_grid[2] = { ctx->ARRAYSIZE, ctx->ARRAYSIZE };
    int      shared     = (ctx->OUTMODE == 0) ? 1 : 0;

    create_image_ID(ctx->OUTLOWFSSTREAMNAME,
                    2,
                    sz_out,
                    _DATATYPE_FLOAT,
                    shared,
                    0,
                    0,
                    &ctx->IDoutLOWFS);
    if (shared)
    {
        COREMOD_MEMORY_image_set_createsem(ctx->OUTLOWFSSTREAMNAME, 10);
    }

    create_image_ID("aosim_imlowfs", 2, sz_out, _DATATYPE_FLOAT, 1, 0, 0, &ctx->IDimlowfs);
    create_image_ID("aosim_foc1_amp", 2, sz_grid, _DATATYPE_FLOAT, 1, 0, 0, &ctx->IDfoc1a);
    create_image_ID("aosim_foc1_pha", 2, sz_grid, _DATATYPE_FLOAT, 1, 0, 0, &ctx->IDfoc1p);
    create_image_ID("aosim_pup1t_amp", 2, sz_grid, _DATATYPE_FLOAT, 1, 0, 0, &ctx->IDpup1ta);
    create_image_ID("aosim_pup1t_pha", 2, sz_grid, _DATATYPE_FLOAT, 1, 0, 0, &ctx->IDpup1tp);
    create_image_ID("aosim_pup1r_amp", 2, sz_grid, _DATATYPE_FLOAT, 1, 0, 0, &ctx->IDpup1ra);
    create_image_ID("aosim_pup1r_pha", 2, sz_grid, _DATATYPE_FLOAT, 1, 0, 0, &ctx->IDpup1rp);

    ctx->xsizein = dcimg[ctx->IDinOPD].md[0].size[0];
    ctx->ysizein = dcimg[ctx->IDinOPD].md[0].size[1];
    create_2Dimage_ID("aosim_wfa", ctx->ARRAYSIZE, ctx->ARRAYSIZE, &ctx->IDwfa);
    create_2Dimage_ID("aosim_wfp", ctx->ARRAYSIZE, ctx->ARRAYSIZE, &ctx->IDwfp);
    ctx->iioffset = ((long) ctx->ARRAYSIZE - (long) ctx->xsizein) / 2;
    ctx->jjoffset = ((long) ctx->ARRAYSIZE - (long) ctx->ysizein) / 2;

    create_2Dimage_ID("aosim_imcamlowfstmp",
                      ctx->ARRAYSIZE,
                      ctx->ARRAYSIZE,
                      &ctx->IDimcamlowfstmp);
    ctx->LOWFScamstarttime = dcimg[ctx->IDphystime].array.F[0];
    ctx->LOWFScamcnt       = 0;

    list_image_ID();
}
