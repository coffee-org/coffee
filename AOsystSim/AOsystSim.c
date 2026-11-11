/**
 * @file    AOsystSim.c
 * @brief   Coupled real-time Adaptive Optics loop simulator and shared module state
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
#include "OptSystProp/OptSystProp.h"
#include "fft/fft.h"
#include "image_filter/image_filter.h"
#include "image_gen/image_gen.h"
#include "statistic/statistic.h"

#include "AOsystSim/AOsystSim.h"
#include "AOsystSim_internal.h"

static int WDIR_INIT = 0;

OPTSYST *optsystsim;

/**
 * @brief Internal state for the coupled real-time AO simulation runner.
 */
struct aosim_run_ctx
{
    long     arraysize;
    double   puprad;
    double   dmrad;
    double   dmifscale;
    long     DMsize;
    float    pupradcoeff;
    float    dmradcoeff;
    float    iwald;
    int      COROmode;
    double   dftzoomfact;
    imageID  IDdmctrl;
    imageID  IDifc;
    imageID  IDturb;
    imageID  IDdm0shape;
    imageID  IDout;
    imageID  IDdhmask;
    imageID  IDarray[2];
    uint32_t dhsizearray[2];
    long     dhxsize;
    long     dhysize;
    long     dhxoffset;
    long     dhyoffset;
};

/**
 * aosim_ensure_wdir - Create the AOsystSim_wdir working directory on first invocation
 */
void aosim_ensure_wdir(void)
{
    if (WDIR_INIT == 0)
    {
        EXECUTE_SYSTEM_COMMAND("mkdir -p AOsystSim_wdir");
        WDIR_INIT = 1;
    }
}

/**
 * aosim_sync_wait - Synchronize AO simulation loop iteration on stream semaphores or delay
 * @syncmode: 0=turbulence stream, 1=DM stream, 2=either stream, >2=fixed microsecond delay
 * @IDarray:  Two-element array of turbulence and DM image IDs
 * @delayus:  Fallback sleep duration in microseconds when syncmode > 2
 */
static void aosim_sync_wait(
    int            syncmode,
    const imageID *IDarray,
    long           delayus)
{
    if (syncmode == 0)
    {
        ImageStreamIO_semwait(&dcimg[IDarray[0]], 0);
        return;
    }
    if (syncmode == 1)
    {
        ImageStreamIO_semwait(&dcimg[IDarray[1]], 0);
        return;
    }
    if (syncmode == 2)
    {
        while (!DCSIG_ANY_SET())
        {
            if (ImageStreamIO_semtrywait(&dcimg[IDarray[0]], 0) == 0 ||
                ImageStreamIO_semtrywait(&dcimg[IDarray[1]], 0) == 0)
            {
                return;
            }
            usleep(100);
        }
        return;
    }
    printf("WAITING %ld us\n", delayus);
    usleep(delayus);
}

/**
 * aosim_run_init_dmctrl - Attach or create the shared-memory DM actuator command stream
 * @ctx:     Simulation runner context
 * @DMindex: Numeric DM stream index
 */
static void aosim_run_init_dmctrl(
    struct aosim_run_ctx *ctx,
    long                  DMindex)
{
    char name[200];

    sprintf(name, "dm%02lddisp", DMindex);
    ctx->IDdmctrl = image_ID(name, dcimg, dcnimg);
    if (ctx->IDdmctrl == -1)
    {
        ctx->IDdmctrl = read_sharedmem_image(name, dcimg, dcnimg);
    }
    if (ctx->IDdmctrl == -1)
    {
        uint32_t dmsizearray[2];

        dmsizearray[0] = (uint32_t) ctx->DMsize;
        dmsizearray[1] = (uint32_t) ctx->DMsize;
        create_image_ID(name, 2, dmsizearray, _DATATYPE_FLOAT, 1, 0, 0, &ctx->IDdmctrl);
        COREMOD_MEMORY_image_set_createsem(name, 2);
    }
    else
    {
        ctx->DMsize = dcimg[ctx->IDdmctrl].md[0].size[0];
        COREMOD_MEMORY_image_set_createsem(name, 2);
    }

    for (long kk = 0; kk < ctx->DMsize * ctx->DMsize; kk++)
    {
        dcimg[ctx->IDdmctrl].array.F[kk] = ran1() * 2.0e-8;
    }
    ctx->dmifscale = 0.5 * ctx->DMsize;
}

/**
 * aosim_run_build_dmifc - Construct Gaussian-smoothed 3D DM influence-function cube
 * @ctx: Simulation runner context
 */
static void aosim_run_build_dmifc(struct aosim_run_ctx *ctx)
{
    long     arraysize = ctx->arraysize;
    long     DMsize    = ctx->DMsize;
    long     DMnbact   = DMsize * DMsize;
    uint32_t imsize[3];
    imageID  IDif;
    double   actpix = ctx->dmifscale * (2.0 * ctx->dmrad / DMsize);
    double   sig;

    imsize[0] = (uint32_t) arraysize;
    imsize[1] = (uint32_t) arraysize;
    imsize[2] = (uint32_t) DMnbact;
    create_image_ID("dmif0", 2, imsize, _DATATYPE_FLOAT, 0, 0, 0, &IDif);
    printf("DM size = %ld x %ld actuators\n", DMsize, DMsize);
    printf("actuator pix size = %f pix\n", actpix);

    list_image_ID();
    for (long ii = 0; ii < arraysize; ii++)
    {
        for (long jj = 0; jj < arraysize; jj++)
        {
            double x = (1.0 * ii - 0.5 * arraysize) / actpix;
            double y = (1.0 * jj - 0.5 * arraysize) / actpix;

            if ((fabs(x) < 0.5) && (fabs(y) < 0.5))
            {
                dcimg[IDif].array.F[jj * arraysize + ii] = 1.0e-6f;
            }
            else
            {
                dcimg[IDif].array.F[jj * arraysize + ii] = 0.0f;
            }
        }
    }

    printf("convolve\n");
    fflush(stdout);
    save_fits("dmif0", "AOsystSim_wdir/dmif0.fits");
    sig = 0.5 * actpix;
    printf("gauss filter   %lf %ld\n", sig, (long) (2.0 * sig));
    fflush(stdout);
    gauss_filter("dmif0", "dmif", sig, (long) (2.0 * sig));
    list_image_ID();
    delete_image_ID("dmif0", DELETE_IMAGE_ERRMODE_WARNING);
    IDif = image_ID("dmif", dcimg, dcnimg);
    list_image_ID();
    save_fits("dmif", "AOsystSim_wdir/dmif.fits");

    create_image_ID("dmifc", 3, imsize, _DATATYPE_FLOAT, 0, 0, 0, &ctx->IDifc);
    printf("\n");
    list_image_ID();
    printf("dmifc = %ld   %ld %ld %ld    %ld %ld\n",
           ctx->IDifc,
           (long) dcimg[ctx->IDifc].md[0].size[0],
           (long) dcimg[ctx->IDifc].md[0].size[1],
           (long) dcimg[ctx->IDifc].md[0].size[2],
           arraysize,
           arraysize);

    for (long mx = 0; mx < DMsize; mx++)
    {
        for (long my = 0; my < DMsize; my++)
        {
            double ii0f = 0.5 * arraysize + 2.0 * (mx - DMsize / 2) / DMsize * ctx->dmrad;
            double jj0f = 0.5 * arraysize + 2.0 * (my - DMsize / 2) / DMsize * ctx->dmrad;

            printf("\r actuator %2ld %2ld    ", mx, my);
            fflush(stdout);
            for (long ii = 0; ii < arraysize; ii++)
            {
                for (long jj = 0; jj < arraysize; jj++)
                {
                    double rxif = (1.0 * ii - ii0f) * ctx->dmifscale + arraysize / 2;
                    double ryif = (1.0 * jj - jj0f) * ctx->dmifscale + arraysize / 2;
                    long   rxi  = (long) rxif;
                    long   ryi  = (long) ryif;

                    if ((rxi > 0) && (rxi < arraysize - 1) && (ryi > 0) && (ryi < arraysize - 1))
                    {
                        double u   = rxif - rxi;
                        double t   = ryif - ryi;
                        double v00 = dcimg[IDif].array.F[ryi * arraysize + rxi];
                        double v01 = dcimg[IDif].array.F[ryi * arraysize + rxi + 1];
                        double v10 = dcimg[IDif].array.F[(ryi + 1) * arraysize + rxi];
                        double v11 = dcimg[IDif].array.F[(ryi + 1) * arraysize + rxi + 1];
                        long   idx = (my * DMsize + mx) * arraysize * arraysize +
                                     jj * arraysize + ii;

                        dcimg[ctx->IDifc].array.F[idx] =
                            (float) ((1.0 - u) * (1.0 - t) * v00 + (1.0 - u) * t * v10 +
                                     u * (1.0 - t) * v01 + u * t * v11);
                    }
                }
            }
        }
    } // for (long mx)
    save_fits("dmifc", "AOsystSim_wdir/dmifc.fits");
    printf("\n");
}

/**
 * aosim_run_init_coro_fpm - Add optional opaque focal-plane coronagraph mask element
 * @arraysize:   Grid size in pixels
 * @dftzoomfact: Zoomed DFT oversampling factor
 * @elem:        Optical element slot index in optsystsim[0]
 *
 * Return: Updated optical element count.
 */
static int aosim_run_init_coro_fpm(
    long   arraysize,
    double dftzoomfact,
    int    elem)
{
    imageID IDfocmask;

    create_2DCimage_ID("focpm", arraysize, arraysize, &IDfocmask);
    for (long ii = 0; ii < arraysize; ii++)
    {
        for (long jj = 0; jj < arraysize; jj++)
        {
            double x = 1.0 * ii - 0.5 * arraysize;
            double y = 1.0 * jj - 0.5 * arraysize;
            double r = sqrt(x * x + y * y);

            if (r < 20.0 * dftzoomfact)
            {
                dcimg[IDfocmask].array.CF[jj * arraysize + ii].re = 1.0f;
                dcimg[IDfocmask].array.CF[jj * arraysize + ii].im = 0.0f;
            }
            else
            {
                dcimg[IDfocmask].array.CF[jj * arraysize + ii].re = 0.0f;
                dcimg[IDfocmask].array.CF[jj * arraysize + ii].im = 0.0f;
            }
        }
    }
    optsystsim[0].elemtype[elem]          = 5;
    optsystsim[0].FOCMASKarray[0].fpmID   = IDfocmask;
    optsystsim[0].FOCMASKarray[0].zfactor = dftzoomfact;
    optsystsim[0].FOCMASKarray[0].mode    = 1;
    optsystsim[0].elemZpos[elem]          = optsystsim[0].elemZpos[elem - 1];
    return elem + 1;
}

/**
 * aosim_run_init_optsyst - Initialize turbulence stream, OPTSYST model, and PSF output stream
 * @ctx:     Simulation runner context
 * @DMindex: Numeric DM stream index
 */
static void aosim_run_init_optsyst(
    struct aosim_run_ctx *ctx,
    long                  DMindex)
{
    long     arraysize = ctx->arraysize;
    uint32_t imsize2[2];
    uint32_t imsize3[3];
    char     name[200];
    imageID  IDpupm;
    imageID  ID;
    int      elem = 0;

    imsize2[0] = (uint32_t) arraysize;
    imsize2[1] = (uint32_t) arraysize;
    create_image_ID("WFturb", 2, imsize2, _DATATYPE_FLOAT, 1, 0, 0, &ctx->IDturb);
    COREMOD_MEMORY_image_set_createsem("WFturb", 2);
    list_image_ID();

    sprintf(name, "dm%02lddisp", DMindex);
    AOsystSim_DMshape(name, "dmifc", "dm2Ddisp");
    ctx->IDdm0shape = image_ID("dm2Ddisp", dcimg, dcnimg);
    save_fits("dm2Ddisp", "AOsystSim_wdir/dm2Ddisp.fits");

    optsystsim = (OPTSYST *) malloc(sizeof(OPTSYST));
    if (optsystsim == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }
    optsystsim[0].nblambda       = 1;
    optsystsim[0].lambdaarray[0] = 1.6e-6;
    optsystsim[0].beamrad        = 0.008;
    optsystsim[0].size           = arraysize;
    optsystsim[0].pixscale       = optsystsim[0].beamrad / 50.0;
    optsystsim[0].DFTgridpad     = 0;
    optsystsim[0].NB_asphsurfm   = 2;
    optsystsim[0].NB_asphsurfr   = 0;
    optsystsim[0].NBelem         = 100;

    IDpupm                             = make_disk("pupmask", arraysize, arraysize,
                                                   0.5 * arraysize, 0.5 * arraysize, ctx->puprad);
    optsystsim[0].elemtype[elem]       = 1;
    optsystsim[0].elemarrayindex[elem] = IDpupm;
    optsystsim[0].elemZpos[elem]       = 0.0;
    elem++;

    optsystsim[0].elemtype[elem]           = 3;
    optsystsim[0].elemarrayindex[elem]     = 0;
    optsystsim[0].ASPHSURFMarray[0].surfID = ctx->IDturb;
    optsystsim[0].elemZpos[elem]           = 0.0;
    elem++;

    optsystsim[0].elemtype[elem]           = 3;
    optsystsim[0].elemarrayindex[elem]     = 1;
    optsystsim[0].ASPHSURFMarray[1].surfID = ctx->IDdm0shape;
    optsystsim[0].elemZpos[elem]           = 0.0;
    optsystsim[0].keepMem[elem]            = 1;
    elem++;

    if (ctx->COROmode == 1)
    {
        elem = aosim_run_init_coro_fpm(arraysize, ctx->dftzoomfact, elem);
    }
    optsystsim[0].NBelem = elem;
    optsystsim[0].SAVE   = 1;

    OptSystProp_run(optsystsim, 0, 0, optsystsim[0].NBelem, "./testconf/", 1);

    ID         = image_ID("psfi0", dcimg, dcnimg);
    imsize3[0] = dcimg[ID].md[0].size[0];
    imsize3[1] = dcimg[ID].md[0].size[1];
    imsize3[2] = dcimg[ID].md[0].size[2];
    create_image_ID("aosimpsfout", 3, imsize3, _DATATYPE_FLOAT, 1, 0, 0, &ctx->IDout);
    COREMOD_MEMORY_image_set_createsem("aosimpsfout", 2);
    dcimg[ctx->IDout].md[0].write = 1;
    memcpy(dcimg[ctx->IDout].array.F,
           dcimg[ID].array.F,
           sizeof(float) * imsize3[0] * imsize3[1] * imsize3[2]);
    dcimg[ctx->IDout].md[0].cnt0++;
    dcimg[ctx->IDout].md[0].write = 0;
    COREMOD_MEMORY_image_set_sempost("aosimpsfout", -1);
}

/**
 * aosim_run_init_dhmask - Initialize synchronization IDs and focal-plane dark-hole mask
 * @ctx:     Simulation runner context
 * @DMindex: Numeric DM stream index
 */
static void aosim_run_init_dhmask(
    struct aosim_run_ctx *ctx,
    long                  DMindex)
{
    char name[200];
    long dhsize;

    ctx->IDarray[0] = image_ID("WFturb", dcimg, dcnimg);
    sprintf(name, "dm%02lddisp", DMindex);
    ctx->IDarray[1] = image_ID(name, dcimg, dcnimg);

    dhsize = (long) (0.5 / ctx->pupradcoeff * ctx->DMsize *
                     (ctx->pupradcoeff / ctx->dmradcoeff) * 0.5);
    ctx->dhxsize        = dhsize;
    ctx->dhysize        = dhsize * 2;
    ctx->dhxoffset      = ctx->arraysize / 2;
    ctx->dhyoffset      = ctx->arraysize / 2 - dhsize;
    ctx->dhsizearray[0] = (uint32_t) (ctx->dhxsize * 2);
    ctx->dhsizearray[1] = (uint32_t) ctx->dhysize;

    create_2Dimage_ID("dhmask", ctx->dhxsize, ctx->dhysize, &ctx->IDdhmask);
    for (long ii = 0; ii < ctx->dhxsize; ii++)
    {
        for (long jj = 0; jj < ctx->dhysize; jj++)
        {
            long   ii1   = (ii + ctx->dhxoffset) - ctx->arraysize / 2;
            long   jj1   = (jj + ctx->dhyoffset) - ctx->arraysize / 2;
            double ii1ld = 2.0 * ii1 * ctx->pupradcoeff;
            double rld   = sqrt((double) (ii1 * ii1 + jj1 * jj1)) * ctx->pupradcoeff * 2.0;

            dcimg[ctx->IDdhmask].array.F[jj * ctx->dhxsize + ii] = 1.0f;
            if ((ii1ld < 0.5) || (rld < ctx->iwald))
            {
                dcimg[ctx->IDdhmask].array.F[jj * ctx->dhxsize + ii] = 0.0f;
            }
        }
    }
    save_fits("dhmask", "AOsystSim_wdir/dhmask.fits");
}

/**
 * aosim_run_step - Execute one iteration of DM update, beam propagation, WFS, and dark-hole output
 * @ctx:      Simulation runner context
 * @syncmode: Synchronization mode
 * @DMindex:  Numeric DM stream index
 * @delayus:  Fallback delay in microseconds
 */
static void aosim_run_step(
    struct aosim_run_ctx *ctx,
    int                   syncmode,
    long                  DMindex,
    long                  delayus)
{
    char    name[200];
    char    imnameamp[200];
    char    imnamepha[200];
    long    index = 2;
    imageID ID;
    imageID IDdh;
    imageID IDre;
    imageID IDim;

    sprintf(name, "dm%02lddisp", DMindex);
    printf("Compute DM shape ...\n");
    fflush(stdout);
    AOsystSim_DMshape(name, "dmifc", "dm2Ddisp");
    printf("done\n");
    fflush(stdout);

    printf("Computing propagation ...\n");
    fflush(stdout);
    OptSystProp_run(optsystsim, 0, 0, optsystsim[0].NBelem, "./testconf/", 1);
    printf("done\n");
    fflush(stdout);

    if (sprintf(imnameamp, "WFamp0_%03ld", index) < 1)
    {
        PRINT_ERROR("sprintf wrote <1 char");
    }
    if (sprintf(imnamepha, "WFpha0_%03ld", index) < 1)
    {
        PRINT_ERROR("sprintf wrote <1 char");
    }

    mk_complex_from_amph(imnameamp, imnamepha, "_tmpwfc", 0);
    AOsystSim_WFSsim_Pyramid("_tmpwfc", "aosimwfsim", 0.0, 1);
    delete_image_ID("_tmpwfc", DELETE_IMAGE_ERRMODE_WARNING);
    COREMOD_MEMORY_image_set_sempost("aosimwfsim", 0);

    ID                            = image_ID("psfi0", dcimg, dcnimg);
    dcimg[ctx->IDout].md[0].write = 1;
    memcpy(dcimg[ctx->IDout].array.F,
           dcimg[ID].array.F,
           sizeof(float) * dcimg[ID].md[0].size[0] * dcimg[ID].md[0].size[1] *
               dcimg[ID].md[0].size[2]);
    dcimg[ctx->IDout].md[0].cnt0++;
    dcimg[ctx->IDout].md[0].write = 0;
    COREMOD_MEMORY_image_set_sempost("aosimpsfout", -1);

    IDre = image_ID("psfre0", dcimg, dcnimg);
    IDim = image_ID("psfim0", dcimg, dcnimg);
    create_image_ID("dhfield", 2, ctx->dhsizearray, _DATATYPE_FLOAT, 1, 0, 0, &IDdh);
    dcimg[IDdh].md[0].write = 1;
    for (long ii = 0; ii < ctx->dhxsize; ii++)
    {
        for (long jj = 0; jj < ctx->dhysize; jj++)
        {
            long  src = (jj + ctx->dhyoffset) * ctx->arraysize + (ii + ctx->dhxoffset);
            float msk = dcimg[ctx->IDdhmask].array.F[jj * ctx->dhxsize + ii];

            dcimg[IDdh].array.F[jj * (2 * ctx->dhxsize) + ii] =
                dcimg[IDre].array.F[src] * msk;
            dcimg[IDdh].array.F[jj * (2 * ctx->dhxsize) + (ii + ctx->dhxsize)] =
                dcimg[IDim].array.F[src] * msk;
        }
    }
    dcimg[IDdh].md[0].cnt0++;
    dcimg[IDdh].md[0].write = 0;

    aosim_sync_wait(syncmode, ctx->IDarray, delayus);
    COREMOD_MEMORY_image_set_semflush_IDarray(ctx->IDarray, 2);
}

/**
 * AOsystSim_run - Run coupled real-time shared-memory AO loop simulator
 * @syncmode: Stream synchronization mode (0:turb, 1:DM, 2:either, >2:fixed delay)
 * @DMindex:  Deformable mirror stream index
 * @delayus:  Loop delay in microseconds when syncmode > 2
 *
 * Return: 0 on completion.
 */
int AOsystSim_run(
    int  syncmode,
    long DMindex,
    long delayus)
{
    struct aosim_run_ctx ctx = { 0 };

    ctx.arraysize   = 128;
    ctx.dmifscale   = 2.0;
    ctx.DMsize      = 50;
    ctx.pupradcoeff = 0.17f;
    ctx.dmradcoeff  = 0.20f;
    ctx.iwald       = 2.0f;
    ctx.COROmode    = 0;
    ctx.dftzoomfact = 2.0;
    ctx.puprad      = ctx.pupradcoeff * ctx.arraysize;
    ctx.dmrad       = ctx.dmradcoeff * ctx.arraysize;

    aosim_ensure_wdir();
    aosim_run_init_dmctrl(&ctx, DMindex);
    aosim_run_build_dmifc(&ctx);
    aosim_run_init_optsyst(&ctx, DMindex);
    aosim_run_init_dhmask(&ctx, DMindex);

    for (;;)
    {
        aosim_run_step(&ctx, syncmode, DMindex, delayus);
    }

    return 0;
}
