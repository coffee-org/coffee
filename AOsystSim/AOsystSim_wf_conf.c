/**
 * @file    AOsystSim_wf_conf.c
 * @brief   Configuration file generation, parsing, and stream initialization for AOsystSim_mkWF
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "image_gen/image_gen.h"

#include "AOsystSim/AOsystSim.h"
#include "AOsystSim_internal.h"

/**
 * AOsystSim_mkWF_mkCONF - Write template configuration file for multi-layer atmospheric WF stream
 * @fname: Output configuration file path
 *
 * Return: 0 on completion.
 */
int AOsystSim_mkWF_mkCONF(const char *fname)
{
    FILE *fp = fopen(fname, "w");

    fprintf(fp, "# ======== AO simulation : creating atmospheric wavefronts stream "
                "=================\n");
    fprintf(fp, "# script relies on pre-computed physical atmospheric wavefronts\n\n\n");
    fprintf(fp, "MODE             0                "
                "# 0: read pre-computed WFs, 1: empty WFs\n");
    fprintf(fp, "WFDIR            ./atmwf          "
                "# atmospheric wavefronts directory\n\n");
    fprintf(fp, "LAMBDANM         1650             # wavelength [nm]\n\n");
    fprintf(fp, "# ============== OUTPUT TYPE (OPD unit = um) ====================\n");
    fprintf(fp, "OUT0FITSFILE         2                "
                "# 0: none, 1 if FITS file output OPD only, 2 if OPD+AMP\n");
    fprintf(fp, "OUT0FITSFILENAMEOPD  wf0opd.fits      # FITS file name output (OPD)\n");
    fprintf(fp, "OUT0FITSFILENAMEAMP  wf0amp.fits      # FITS file name output (AMP)\n");
    fprintf(fp, "OUTPHYSTIME          aosim_phystime   # physical time [s]\n");
    fprintf(fp, "OUT0STREAM           2                "
                "# 1 if shared memory stream OPD only, 2 if OPD+AMP \n");
    fprintf(fp, "OUT0STREAMNAMEOPD    wf0opd           # output WF stream name (OPD)\n");
    fprintf(fp, "OUT0STREAMNAMEAMP    wf0amp           # output WF stream name (AMP)\n\n");
    fprintf(fp, "# ============== PROCESS TRIGGER ==================================\n");
    fprintf(fp, "TRIGGERMODE           2                "
                "# 0: file, 1: semaphore from stream, 2: time interval\n");
    fprintf(fp, "TRIGGERFILE           wfup.txt         "
                "# update file name (TRIGGERMODE=0,1\n");
    fprintf(fp, "TRIGGERDT             0.05             "
                "# update interval (TRIGGERMODE=0,2)\n");
    fprintf(fp, "TRIGGER0STREAM        WFSinst          # trigger stream\n");
    fprintf(fp, "TRIGGER0SEM           5                "
                "# trigger semaphore in TRIGGERSTREAM\n");
    fprintf(fp, "TRIGGER1STREAM        dm05dispmap      # trigger stream\n");
    fprintf(fp, "TRIGGER1SEM           5                "
                "# trigger semaphore in TRIGGERSTREAM\n\n");
    fprintf(fp, "# ============== PARAMETERS ======================================\n");
    fprintf(fp, "DT                    0.001            "
                "# time interval between computed wavefronts\n");
    fprintf(fp, "OPDMFACT              1.0              # OPD multiplicative factor\n");
    fprintf(fp, "AMPMFACT              1.0              # AMP multiplicative factor\n");
    fprintf(fp, "ARRAYSIZE             128              # output size [pix]\n");
    fprintf(fp, "PIXSCALEMODE          2                "
                "# 0: adopt input WF pixel scale, 1: custom pixel scale, 2: bin\n");
    fprintf(fp, "PIXSCALECUSTOM        0.1              # custom pixel size\n");
    fprintf(fp, "PIXBINFACTOR          4                # bin factor\n");
    fprintf(fp, "PUPDIAMM              8                # Pupil diameter [m]\n");
    fprintf(fp, "PUPILFILE    aohardsim/tpup.fits.gz    "
                "# Pupil file (optional), full image size is pupil diam\n\n");
    fprintf(fp, "# ============== POST-PROCESSING =============================\n");
    fprintf(fp, "DM0MODE               1                # 0 if no DM0, 1 if OPD DMn\n");
    fprintf(fp, "DM0NAME               dm05dispmap      # DM displacement map stream\n");
    fprintf(fp, "OUT1FITSFILE          0                # 1 if FITS file output\n");
    fprintf(fp, "OUT1FITSFILENAMEOPD  wf1opd.fits       # FITS file name output (OPD)\n");
    fprintf(fp, "OUT1FITSFILENAMEAMP  wf1amp.fits       # FITS file name output (AMP)\n");
    fprintf(fp, "OUT1STREAM           1                 "
                "# 1 if shared memory stream OPD only, 2 if OPD+AMP \n");
    fprintf(fp, "OUT1STREAMNAMEOPD     wf1opd           # output WF stream name (OPD)\n");
    fprintf(fp, "OUT1STREAMNAMEAMP     wf1amp           # output WF stream name (AMP)\n\n\n");
    fclose(fp);

    return 0;
}

/**
 * aosim_mkwf_read_conf - Parse atmospheric wavefront configuration files into context structure
 * @conf_fname: Top-level WF daemon configuration file path
 * @ctx:        Wavefront generator runtime context
 */
void aosim_mkwf_read_conf(
    const char     *conf_fname,
    AOSIM_MKWF_CTX *ctx)
{
    FILE *fp;
    FILE *fp1;
    char  fname[STRINGMAXLEN_FULLFILENAME];

    if ((fp = fopen(conf_fname, "r")) == NULL)
    {
        sprintf(fname, "%s.default", conf_fname);
        printf("configuration file %s not found. Creating default template as %s\n",
               conf_fname,
               fname);
        AOsystSim_mkWF_mkCONF(fname);
        exit(0);
    }
    fclose(fp);

    read_config_parameter(conf_fname, "WFDIR", ctx->WFDIR);
    printf("WFDIR = %s\n", ctx->WFDIR);
    sprintf(ctx->conf_atmWFconf_fname, "%s/WFsim.conf", ctx->WFDIR);
    printf("WFsim.conf   : %s\n", ctx->conf_atmWFconf_fname);
    read_config_parameter(ctx->conf_atmWFconf_fname, "SWF_FILE_PREFIX", ctx->wfin_PREFIX);

    ctx->LAMBDA = 1.0e-9f * read_config_parameter_float(conf_fname, "LAMBDANM");
    printf("LAMBDA = %.3f nm\n", ctx->LAMBDA * 1.0e9);

    ctx->OUT0FITSFILE = read_config_parameter_long(conf_fname, "OUT0FITSFILE");
    read_config_parameter(conf_fname, "OUT0FITSFILENAMEOPD", ctx->OUT0FITSFILENAMEOPD);
    read_config_parameter(conf_fname, "OUT0FITSFILENAMEAMP", ctx->OUT0FITSFILENAMEAMP);
    read_config_parameter(conf_fname, "OUTPHYSTIME", ctx->OUTPHYSTIME);
    ctx->OUT0STREAM = read_config_parameter_long(conf_fname, "OUT0STREAM");
    read_config_parameter(conf_fname, "OUT0STREAMNAMEOPD", ctx->OUT0STREAMNAMEOPD);
    read_config_parameter(conf_fname, "OUT0STREAMNAMEAMP", ctx->OUT0STREAMNAMEAMP);

    ctx->TRIGGERMODE = read_config_parameter_long(conf_fname, "TRIGGERMODE");
    ctx->TRIGGERDT   = read_config_parameter_float(conf_fname, "TRIGGERDT");
    read_config_parameter(conf_fname, "TRIGGERFILE", ctx->TRIGGERFILE);
    read_config_parameter(conf_fname, "TRIGGER0STREAM", ctx->TRIGGER0STREAM);
    ctx->TRIGGER0SEM = read_config_parameter_long(conf_fname, "TRIGGER0SEM");
    read_config_parameter(conf_fname, "TRIGGER1STREAM", ctx->TRIGGER1STREAM);
    ctx->TRIGGER1SEM = read_config_parameter_long(conf_fname, "TRIGGER1SEM");

    ctx->DT = read_config_parameter_float(conf_fname, "DT");
    printf("DT = %.3f ms\n", ctx->DT * 1.0e3);
    ctx->OPDMFACT  = read_config_parameter_float(conf_fname, "OPDMFACT");
    ctx->AMPMFACT  = read_config_parameter_float(conf_fname, "AMPMFACT");
    ctx->ARRAYSIZE = read_config_parameter_long(conf_fname, "ARRAYSIZE");
    printf("ARRAYSIZE = %ld   (%s)\n", ctx->ARRAYSIZE, conf_fname);

    ctx->PIXSCALEMODE = read_config_parameter_long(conf_fname, "PIXSCALEMODE");
    ctx->PIXBINFACTOR = read_config_parameter_long(conf_fname, "PIXBINFACTOR");
    ctx->PUPDIAM      = read_config_parameter_float(conf_fname, "PUPDIAMM");
    read_config_parameter(conf_fname, "PUPILFILE", ctx->PUPILFILE);

    ctx->DM0MODE = read_config_parameter_long(conf_fname, "DM0MODE");
    read_config_parameter(conf_fname, "DM0NAME", ctx->DM0NAME);
    read_config_parameter(conf_fname, "OUT1FITSFILENAMEOPD", ctx->OUT1FITSFILENAMEOPD);
    read_config_parameter(conf_fname, "OUT1FITSFILENAMEAMP", ctx->OUT1FITSFILENAMEAMP);
    ctx->OUT1STREAM = read_config_parameter_long(conf_fname, "OUT1STREAM");
    read_config_parameter(conf_fname, "OUT1STREAMNAMEOPD", ctx->OUT1STREAMNAMEOPD);
    read_config_parameter(conf_fname, "OUT1STREAMNAMEAMP", ctx->OUT1STREAMNAMEAMP);

    ctx->wfin_TIME_STEP = read_config_parameter_float(ctx->conf_atmWFconf_fname, "WFTIME_STEP");
    printf("Input WF time step = %f\n", ctx->wfin_TIME_STEP);
    ctx->wfin_NB_TSPAN = read_config_parameter_long(ctx->conf_atmWFconf_fname, "NB_TSPAN");
    printf("Input WF NB_TSPAN = %ld\n", ctx->wfin_NB_TSPAN);
    ctx->wfin_TIME_SPAN = read_config_parameter_float(ctx->conf_atmWFconf_fname, "TIME_SPAN");
    printf("Input WF TIME_SPAN = %f\n", ctx->wfin_TIME_SPAN);

    ctx->NBframes = (long) ((ctx->wfin_TIME_SPAN + 0.1f * ctx->wfin_TIME_STEP) /
                            ctx->wfin_TIME_STEP);
    printf("NBframes = %ld\n", ctx->NBframes);

    ctx->wfin_PUPIL_SCALE = read_config_parameter_float(ctx->conf_atmWFconf_fname, "PUPIL_SCALE");
    printf("Input WF PUPIL_SCALE = %f\n", ctx->wfin_PUPIL_SCALE);
    aosim_ensure_wdir();
    fp1 = fopen("AOsystSim_wdir/pupscale.txt", "w");
    fprintf(fp1, "%f", ctx->wfin_PUPIL_SCALE);
    fclose(fp1);

    ctx->wfin_WFsize = read_config_parameter_long(ctx->conf_atmWFconf_fname, "WFsize");
    printf("Input WF WFsize = %ld\n", ctx->wfin_WFsize);
}

/**
 * aosim_mkwf_create_streams - Create output wavefront streams and attach optional DM0 stream
 * @ctx: Wavefront generator runtime context
 */
void aosim_mkwf_create_streams(AOSIM_MKWF_CTX *ctx)
{
    uint32_t sz_time[2] = { 1, 1 };
    uint32_t sz_wf[2];
    int      sh0;

    create_image_ID(ctx->OUTPHYSTIME, 2, sz_time, _DATATYPE_FLOAT, 1, 0, 0, &ctx->IDphystime);
    ctx->phystime                     = 0.0;
    dcimg[ctx->IDphystime].array.F[0] = 0.0f;

    if (ctx->PIXSCALEMODE != 2)
    {
        ctx->PIXBINFACTOR = 1;
    }

    ctx->iioffset = (ctx->wfin_WFsize / ctx->PIXBINFACTOR - ctx->ARRAYSIZE) / 2;
    if (ctx->iioffset < 0)
    {
        ctx->iioffset = 0;
    }
    ctx->jjoffset = ctx->iioffset;
    ctx->pupscale = ctx->wfin_PUPIL_SCALE * ctx->PIXBINFACTOR;

    sz_wf[0] = (uint32_t) ctx->ARRAYSIZE;
    sz_wf[1] = (uint32_t) ctx->ARRAYSIZE;
    sh0      = (ctx->OUT0STREAM > 0) ? 1 : 0;

    create_image_ID(ctx->OUT0STREAMNAMEOPD, 2, sz_wf, _DATATYPE_FLOAT, sh0, 0, 0, &ctx->IDopd0);
    create_image_ID(ctx->OUT0STREAMNAMEAMP, 2, sz_wf, _DATATYPE_FLOAT, sh0, 0, 0, &ctx->IDamp0);
    if (sh0)
    {
        COREMOD_MEMORY_image_set_createsem(ctx->OUT0STREAMNAMEOPD, 10);
        COREMOD_MEMORY_image_set_createsem(ctx->OUT0STREAMNAMEAMP, 10);
    }

    printf("DM0MODE = %d\n", ctx->DM0MODE);
    if (ctx->DM0MODE > 0)
    {
        int sh1 = (ctx->OUT1STREAM > 0) ? 1 : 0;

        create_image_ID(ctx->OUT1STREAMNAMEOPD, 2, sz_wf, _DATATYPE_FLOAT, sh1, 0, 0, &ctx->IDopd1);
        printf("CREATED %s stream  %ld %ld   [%d]\n",
               ctx->OUT1STREAMNAMEOPD,
               (long) sz_wf[0],
               (long) sz_wf[1],
               ctx->OUT1STREAM);
        create_image_ID(ctx->OUT1STREAMNAMEAMP, 2, sz_wf, _DATATYPE_FLOAT, sh1, 0, 0, &ctx->IDamp1);
        if (sh1)
        {
            COREMOD_MEMORY_image_set_createsem(ctx->OUT1STREAMNAMEOPD, 10);
            COREMOD_MEMORY_image_set_createsem(ctx->OUT1STREAMNAMEAMP, 10);
        }

        printf("READING INPUT TRIGGER STREAM ...\n");
        ctx->IDdm0opd = read_sharedmem_image(ctx->DM0NAME, dcimg, dcnimg);
        printf("  IDdm0opd = %ld\n", ctx->IDdm0opd);
        while (ctx->IDdm0opd == -1)
        {
            usleep(100000);
            ctx->IDdm0opd = read_sharedmem_image(ctx->DM0NAME, dcimg, dcnimg);
            printf("  IDdm0opd = %ld\n", ctx->IDdm0opd);
        }

        if ((dcimg[ctx->IDdm0opd].md[0].size[0] != (uint32_t) ctx->ARRAYSIZE) ||
            (dcimg[ctx->IDdm0opd].md[0].size[1] != (uint32_t) ctx->ARRAYSIZE))
        {
            printf("ERROR: stream %s has wrong size: is %ld x %ld, should be %ld x %ld\n",
                   ctx->DM0NAME,
                   (long) dcimg[ctx->IDdm0opd].md[0].size[0],
                   (long) dcimg[ctx->IDdm0opd].md[0].size[1],
                   ctx->ARRAYSIZE,
                   ctx->ARRAYSIZE);
            exit(0);
        }
    } // if (ctx->DM0MODE > 0)
}

/**
 * aosim_mkwf_init_pupil - Scan available WF FITS files and synthesize telescope pupil mask
 * @ctx: Wavefront generator runtime context
 */
void aosim_mkwf_init_pupil(AOSIM_MKWF_CTX *ctx)
{
    char    wf_fname[STRINGMAXLEN_FULLFILENAME];
    int     fOK = 1;
    imageID IDpupil;
    long    arr = ctx->ARRAYSIZE;

    ctx->AMPfile   = 1;
    ctx->kmax      = 0;
    ctx->IDbin_re  = -1;
    ctx->IDbin_im  = -1;
    ctx->IDbin_amp = -1;
    ctx->IDbin_opd = -1;

    while (fOK == 1)
    {
        SNPRINTF_CHECK(wf_fname,
                       STRINGMAXLEN_FULLFILENAME,
                       "%s/%s%08ld.%09ld.pha.fits",
                       ctx->WFDIR,
                       ctx->wfin_PREFIX,
                       ctx->kmax,
                       (long) (1.0e12 * ctx->LAMBDA + 0.5));
        fOK = file_exists(wf_fname);
        ctx->kmax++;
    }
    ctx->kmax -= 2;
    printf("kmax = %ld\n", ctx->kmax);

    ctx->IDampmask = make_disk("pupmask",
                               arr,
                               arr,
                               0.5 * arr,
                               0.5 * arr,
                               0.5 * ctx->PUPDIAM / ctx->pupscale);
    load_fits(ctx->PUPILFILE, "pupilgeom", 0, &IDpupil);

    if (IDpupil != -1)
    {
        long  pupilsize = dcimg[IDpupil].md[0].size[0];
        float rad_pix   = 0.5f * ctx->PUPDIAM / ctx->pupscale;

        printf("LOADED %s\n", ctx->PUPILFILE);
        for (long ii = 0; ii < arr; ii++)
        {
            for (long jj = 0; jj < arr; jj++)
            {
                float x   = (1.0f * ii - 0.5f * arr) / rad_pix;
                float y   = (1.0f * jj - 0.5f * arr) / rad_pix;
                long  ii1 = (long) (0.5f * x * pupilsize + 0.5f * pupilsize);
                long  jj1 = (long) (0.5f * y * pupilsize + 0.5f * pupilsize);

                if ((ii1 > -1) && (ii1 < pupilsize) && (jj1 > -1) && (jj1 < pupilsize))
                {
                    dcimg[ctx->IDampmask].array.F[jj * arr + ii] *=
                        dcimg[IDpupil].array.F[jj1 * pupilsize + ii1];
                }
            }
        }
    }
    else
    {
        printf("COULD NOT FIND %s\n", ctx->PUPILFILE);
    }

    for (long ii = 0; ii < arr * arr; ii++)
    {
        dcimg[ctx->IDamp0].array.F[ii] = dcimg[ctx->IDampmask].array.F[ii];
    }

    ctx->ii1start0 = arr;
    ctx->jj1start0 = arr;
    ctx->ii1end0   = 0;
    ctx->jj1end0   = 0;
    for (long ii1 = 0; ii1 < arr; ii1++)
    {
        for (long jj1 = 0; jj1 < arr; jj1++)
        {
            if (dcimg[ctx->IDampmask].array.F[jj1 * arr + ii1] > 1.0e-6f)
            {
                if (ii1 > ctx->ii1end0)
                {
                    ctx->ii1end0 = ii1;
                }
                if (ii1 < ctx->ii1start0)
                {
                    ctx->ii1start0 = ii1;
                }
                if (jj1 > ctx->jj1end0)
                {
                    ctx->jj1end0 = jj1;
                }
                if (jj1 < ctx->jj1start0)
                {
                    ctx->jj1start0 = jj1;
                }
            }
        }
    }
}
