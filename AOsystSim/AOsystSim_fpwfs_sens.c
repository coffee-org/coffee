/**
 * @file    AOsystSim_fpwfs_sens.c
 * @brief   Monte Carlo focal-plane wavefront sensor sensitivity analysis orchestrator
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "statistic/statistic.h"

#include "AOsystSim/AOsystSim.h"
#include "AOsystSim_internal.h"

/**
 * fpwfs_read_opt_double - Read optional Milk CLI scalar variable into target double
 * @name: Variable name in Milk variable table
 * @fmt:  Optional printf format string (or NULL to suppress log)
 * @dst:  In/out parameter value (preserved if variable is absent)
 */
static void fpwfs_read_opt_double(
    const char *name,
    const char *fmt,
    double     *dst)
{
    imageID vID = variable_ID(name);
    if (vID != -1)
    {
        *dst = coffee_dcvar_f(vID);
        if (fmt != NULL)
        {
            printf(fmt, *dst);
        }
    }
}

/**
 * fpwfs_sens_read_vars - Initialize sensitivity parameters and read optional CLI overrides
 * @ctx: Sensitivity analysis context
 */
static void fpwfs_sens_read_vars(AOSIM_FPWFS_SENS_CTX *ctx)
{
    imageID vID;

    ctx->imSimulMode = (ctx->mapmode == 2) ? 1 : 0;
    ctx->Iflux       = 0.0;
    ctx->mapampl     = 2.0;
    ctx->mapsize     = 41;
    ctx->iterMax     = 100000;
    ctx->FLUXph      = 1.0e4;
    ctx->kxtest      = 21;
    ctx->kytest      = 11;
    ctx->IDprobampC  = -1;
    ctx->mapzsize    = 10;
    ctx->ProbeNoise  = 0.2;
    ctx->RON         = 1.0;
    ctx->CnoiseFloor = 1.0e-5;
    ctx->totFlux     = 1.0e12;
    CENTERprobe      = 1;

    if ((vID = variable_ID("mapzsize")) != -1)
    {
        ctx->mapzsize = (long) (0.1 + coffee_dcvar_f(vID));
        printf("mapzsize = %ld\n", ctx->mapzsize);
    }
    if ((vID = variable_ID("NewWF")) != -1)
    {
        ctx->NewWF = (int) (0.1 + coffee_dcvar_f(vID));
        printf("NewWF = %d\n", ctx->NewWF);
    }
    fpwfs_read_opt_double("ProbeNoise", "ProbeNoise = %g\n", &ctx->ProbeNoise);
    fpwfs_read_opt_double("RON", "RON = %f\n", &ctx->RON);
    fpwfs_read_opt_double("CnoiseFloor", "CnoiseFloor = %g\n", &ctx->CnoiseFloor);
    if ((vID = variable_ID("CENTERprobe")) != -1)
    {
        CENTERprobe = (int) (0.1 + coffee_dcvar_f(vID));
        printf("CENTERprobe = %d\n", CENTERprobe);
    }

    if (ctx->imSimulMode == 1)
    {
        ctx->probeampl = 0.0001;
        ctx->contrast  = 5.0e-8;
        ctx->WFerr     = 0.0;
        fpwfs_read_opt_double("probeampl", "probeamp = %f rad\n", &ctx->probeampl);
        fpwfs_read_opt_double("contrast", "contrast = %e\n", &ctx->contrast);
        fpwfs_read_opt_double("WFerr", "WFerr = %f rad\n", &ctx->WFerr);
        ctx->mapmode = 2;
    }

    if ((vID = variable_ID("totFlux")) != -1)
    {
        ctx->totFlux = coffee_dcvar_f(vID);
        ctx->FLUXph  = coffee_dcvar_f(vID);
        printf("totFlux = %f ph\n", ctx->totFlux);
    }
    if (ctx->mapmode == 2)
    {
        CENTERprobe = 1;
    }
}

/**
 * fpwfs_sens_init_probes - Populate ideal and noisy complex probe coordinates
 * @ctx: Sensitivity analysis context
 */
static void fpwfs_sens_init_probes(AOSIM_FPWFS_SENS_CTX *ctx)
{
    ctx->probeamp = 1.0;
    if (CENTERprobe == 1)
    {
        probe_re[0] = 0.0;
        probe_im[0] = 0.0;
        for (int pr = 1; pr < ctx->NBprobes; pr++)
        {
            double ang   = 2.0 * M_PI / (ctx->NBprobes - 1) * (pr - 1);
            probe_re[pr] = ctx->probeamp * cos(ang);
            probe_im[pr] = ctx->probeamp * sin(ang);
        }
    }
    else
    {
        for (int pr = 0; pr < ctx->NBprobes; pr++)
        {
            double ang   = 2.0 * M_PI / ctx->NBprobes * pr;
            probe_re[pr] = ctx->probeamp * cos(ang);
            probe_im[pr] = ctx->probeamp * sin(ang);
        }
    }

    ctx->execmode = 1;
    if (ctx->mapmode == 0)
    {
        ctx->mapsize = 1;
    }
    ctx->mapxsize = ctx->mapsize;
    ctx->mapysize = ctx->mapsize;

    if (ctx->mode == 1)
    {
        ctx->probe_noise_prop = 0.0;
        for (int pr = 0; pr < ctx->NBprobes; pr++)
        {
            nprobe_re[pr] = probe_re[pr];
            nprobe_im[pr] = probe_im[pr];
        }
    }
    else if (ctx->mode == 2)
    {
        ctx->probe_noise_prop = ctx->ProbeNoise;
        for (int pr = 0; pr < ctx->NBprobes; pr++)
        {
            nprobe_re[pr] =
                probe_re[pr] + ctx->probe_noise_prop * ctx->probeamp * (2.0 * ran1() - 1.0);
            nprobe_im[pr] =
                probe_im[pr] + ctx->probe_noise_prop * ctx->probeamp * (2.0 * ran1() - 1.0);
        }
    }
    else if (ctx->mode == 3)
    {
        ctx->probe_noise_prop = ctx->ProbeNoise;
        double re_re          = 1.0 + (2.0 * ran1() - 1.0) * ctx->probe_noise_prop;
        double re_im          = (2.0 * ran1() - 1.0) * ctx->probe_noise_prop;
        double im_re          = (2.0 * ran1() - 1.0) * ctx->probe_noise_prop;
        double im_im          = 1.0 + (2.0 * ran1() - 1.0) * ctx->probe_noise_prop;
        for (int pr = 0; pr < ctx->NBprobes; pr++)
        {
            nprobe_re[pr] = probe_re[pr] * re_re + probe_im[pr] * im_re;
            nprobe_im[pr] = probe_re[pr] * re_im + probe_im[pr] * im_im;
        }
    }
    else
    {
        printf("ERROR: mode not supported\n");
        ctx->execmode = 0;
    }
}

/**
 * fpwfs_sens_run_imsimul - Simulate probed PSF cube for current Monte Carlo realization
 * @ctx: Sensitivity analysis context
 */
static void fpwfs_sens_run_imsimul(AOSIM_FPWFS_SENS_CTX *ctx)
{
    ctx->IDpsfC = AOsystSim_FPWFS_imsimul(
        ctx->probeampl,
        30.0,
        30.0,
        ctx->contrast,
        ctx->WFerr,
        ctx->totFlux,
        ctx->probe_noise_prop,
        ctx->RON,
        ctx->CnoiseFloor);
    ctx->IDpsfCnoise = image_ID("psfCcropnCn", dcimg, dcnimg);
    ctx->IDprobampC  = image_ID("psfprobeampC", dcimg, dcnimg);
    save_fl_fits("psfC", "psfC.fits");

    double   ave    = 0.0;
    long     avecnt = 0;
    uint32_t sx     = dcimg[ctx->IDpsfC].md[0].size[0];
    uint32_t sy     = dcimg[ctx->IDpsfC].md[0].size[1];

    for (int pr = 1; pr < ctx->NBprobes; pr++)
    {
        for (uint32_t ii = 0; ii < sx * sy; ii++)
        {
            ave += dcimg[ctx->IDpsfC].array.F[(uint32_t) pr * sx * sy + ii];
            avecnt++;
        }
    }
    printf("ave = %g\n", ave / (double) avecnt);
    ctx->mapmode  = 2;
    ctx->mapxsize = sx;
    ctx->mapysize = sy;
}

/**
 * fpwfs_sens_init_maps - Allocate 3D Monte Carlo solution cubes on first realization
 * @ctx: Sensitivity analysis context
 */
static void fpwfs_sens_init_maps(AOSIM_FPWFS_SENS_CTX *ctx)
{
    if (ctx->mapmode > 0)
    {
        long sx = ctx->mapxsize;
        long sy = ctx->mapysize;
        long sz = ctx->mapzsize;

        create_3Dimage_ID("WFSerrmap", sx, sy, sz, &ctx->IDmap);
        create_3Dimage_ID("WFSsol_ptre", sx, sy, sz, &ctx->IDmap_ptre);
        create_3Dimage_ID("WFSsol_ptim", sx, sy, sz, &ctx->IDmap_ptim);
        create_3Dimage_ID("WFSsol_ptre_in", sx, sy, sz, &ctx->IDmap_ptre_in);
        create_3Dimage_ID("WFSsol_ptim_in", sx, sy, sz, &ctx->IDmap_ptim_in);
        create_3Dimage_ID("WFSsol_Iflux", sx, sy, sz, &ctx->IDmap_Iflux);
        create_3Dimage_ID("WFSsol_are", sx, sy, sz, &ctx->IDmap_are);
        create_3Dimage_ID("WFSsol_aim", sx, sy, sz, &ctx->IDmap_aim);
        create_3Dimage_ID("WFSsol_e", sx, sy, sz, &ctx->IDmap_e);
    }
    else
    {
        ctx->IDmap = -1;
    }
    ctx->initmap = 1;
}

/**
 * fpwfs_sens_save_and_cleanup - Save solution cubes and release per-slice scratch images
 * @ctx: Sensitivity analysis context
 */
static void fpwfs_sens_save_and_cleanup(const AOSIM_FPWFS_SENS_CTX *ctx)
{
    char imname[200];

    if (ctx->mapmode > 0)
    {
        save_fits("WFSerrmap", "AOsystSim_wdir/WFSerrmap.fits");
        save_fits("WFSsol_ptre", "AOsystSim_wdir/WFSsol_ptre.fits");
        save_fits("WFSsol_ptim", "AOsystSim_wdir/WFSsol_ptim.fits");
        save_fits("WFSsol_ptre_in", "AOsystSim_wdir/WFSsol_ptre_in.fits");
        save_fits("WFSsol_ptim_in", "AOsystSim_wdir/WFSsol_ptim_in.fits");
        save_fits("WFSsol_Iflux", "AOsystSim_wdir/WFSsol_Iflux.fits");
        if (NBoptVar > 3)
        {
            save_fits("WFSsol_are", "AOsystSim_wdir/WFSsol_are.fits");
        }
        if (NBoptVar > 4)
        {
            save_fits("WFSsol_aim", "AOsystSim_wdir/WFSsol_aim.fits");
            save_fits("WFSsol_e", "AOsystSim_wdir/WFSsol_e.fits");
        }
    }

    delete_image_ID("pupa", DELETE_IMAGE_ERRMODE_WARNING);
    if (ctx->NewWF == 1)
    {
        delete_image_ID("wf0", DELETE_IMAGE_ERRMODE_WARNING);
    }
    delete_image_ID("wfA", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("wfB", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("wf", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("psfC", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("wfc", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("psfCcrop", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("psfprobeampC", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("psfCcropn", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("psfCcropnC", DELETE_IMAGE_ERRMODE_WARNING);
    for (int pr = 0; pr < ctx->NBprobes; pr++)
    {
        snprintf(imname, sizeof(imname), "psfC_%03d", pr);
        delete_image_ID(imname, DELETE_IMAGE_ERRMODE_WARNING);
    }
    list_image_ID();
}

/**
 * fpwfs_ensure_2d_map - Return existing 2D image ID or create a new one
 * @name:  Image name
 * @xsize: X dimension
 * @ysize: Y dimension
 *
 * Return: Valid 2D image ID.
 */
static imageID fpwfs_ensure_2d_map(
    const char *name,
    long        xsize,
    long        ysize)
{
    imageID ID = image_ID(name, dcimg, dcnimg);
    if (ID == -1)
    {
        create_2Dimage_ID(name, xsize, ysize, &ID);
    }
    return ID;
}

/**
 * fpwfs_sens_stats_iflux - Compute and export incoherent flux mean and RMS error maps
 * @ctx:  Sensitivity analysis context
 * @mapz: Number of completed Monte Carlo slices
 */
static void fpwfs_sens_stats_iflux(
    const AOSIM_FPWFS_SENS_CTX *ctx,
    long                        mapz)
{
    long    sx    = ctx->mapxsize;
    long    sy    = ctx->mapysize;
    imageID IDave = fpwfs_ensure_2d_map("WFSsol_Iflux_ave", sx, sy);
    imageID IDrms = fpwfs_ensure_2d_map("WFSsol_Iflux_rms", sx, sy);

    for (long kx = 0; kx < sx; kx++)
    {
        for (long ky = 0; ky < sy; ky++)
        {
            double ave = 0.0;
            double rms = 0.0;
            for (long kz = 0; kz < mapz; kz++)
            {
                double v = dcimg[ctx->IDmap_Iflux].array.F[kz * sy * sx + ky * sx + kx];
                ave += v;
                rms += v * v;
            }
            ave /= (double) mapz;
            rms = sqrt(rms / (double) mapz - ave * ave);
            dcimg[IDave].array.F[ky * sx + kx] = (float) ave;
            dcimg[IDrms].array.F[ky * sx + kx] = (float) rms;
        }
    }
    save_fits("WFSsol_Iflux_ave", "AOsystSim_wdir/WFSsol_Iflux_ave.fits");
    save_fits("WFSsol_Iflux_rms", "AOsystSim_wdir/WFSsol_Iflux_rms.fits");

    if (ctx->imSimulMode != 0)
    {
        return;
    }

    imageID IDrmsn  = fpwfs_ensure_2d_map("WFSsol_Iflux_rmsn", sx, sy);
    imageID IDrmsn1 = fpwfs_ensure_2d_map("WFSsol_Iflux_rmsn1", sx, sy);
    double  inv_sn  = 1.0 / sqrt(1.0 / ctx->FLUXph);

    for (long kx = 0; kx < sx; kx++)
    {
        for (long ky = 0; ky < sy; ky++)
        {
            double ptre = ctx->mapampl * (2.0 * (double) kx / (double) sx - 1.0);
            double ptim = ctx->mapampl * (2.0 * (double) ky / (double) sy - 1.0);
            double rval = dcimg[IDrms].array.F[ky * sx + kx] * inv_sn;

            dcimg[IDrmsn].array.F[ky * sx + kx]  = (float) rval;
            dcimg[IDrmsn1].array.F[ky * sx + kx] = (float) (rval / sqrt(ptre * ptre + ptim * ptim));
        }
    }
    save_fits("WFSsol_Iflux_rmsn", "AOsystSim_wdir/WFSsol_Iflux_rmsn.fits");
    save_fits("WFSsol_Iflux_rmsn1", "AOsystSim_wdir/WFSsol_Iflux_rmsn1.fits");
}

/**
 * fpwfs_sens_stats_ca - Compute and export coherent complex-amplitude RMS error maps
 * @ctx:  Sensitivity analysis context
 * @mapz: Number of completed Monte Carlo slices
 */
static void fpwfs_sens_stats_ca(
    const AOSIM_FPWFS_SENS_CTX *ctx,
    long                        mapz)
{
    long    sx     = ctx->mapxsize;
    long    sy     = ctx->mapysize;
    imageID IDrms  = fpwfs_ensure_2d_map("WFSsol_CA_rms", sx, sy);
    imageID IDrmsn = fpwfs_ensure_2d_map("WFSsol_CA_rmsn", sx, sy);
    double  sflux  = sqrt(ctx->FLUXph);

    for (long kx = 0; kx < sx; kx++)
    {
        for (long ky = 0; ky < sy; ky++)
        {
            double ave = 0.0;
            for (long kz = 0; kz < mapz; kz++)
            {
                long   idx = kz * sy * sx + ky * sx + kx;
                double dx  = dcimg[ctx->IDmap_ptre].array.F[idx] -
                             dcimg[ctx->IDmap_ptre_in].array.F[idx];
                double dy  = dcimg[ctx->IDmap_ptim].array.F[idx] -
                             dcimg[ctx->IDmap_ptim_in].array.F[idx];
                ave += dx * dx + dy * dy;
            }
            double ca_rms                       = sqrt(ave / (double) mapz);
            dcimg[IDrms].array.F[ky * sx + kx]  = (float) ca_rms;
            dcimg[IDrmsn].array.F[ky * sx + kx] = (float) (ca_rms * sflux);
        }
    }
    save_fits("WFSsol_CA_rms", "AOsystSim_wdir/WFSsol_CA_rms.fits");
    save_fits("WFSsol_CA_rmsn", "AOsystSim_wdir/WFSsol_CA_rmsn.fits");
}

/**
 * AOsystSim_FPWFS_sensitivityAnalysis - Run Monte Carlo FPWFS sensitivity analysis
 * @mapmode:  Evaluation grid mode (0: single point, 1: 2D CA grid, 2: image cube)
 * @mode:     Probe calibration error model (1: ideal, 2: uncorrelated, 3: axis)
 * @optmode:  Free parameter count selector (0: 3 vars, 1: 4 vars, 2: 6 vars)
 * @NBprobes: Number of focal-plane probes per measurement cycle
 *
 * Return: 0 on completion.
 */
int AOsystSim_FPWFS_sensitivityAnalysis(
    int mapmode,
    int mode,
    int optmode,
    int NBprobes)
{
    AOSIM_FPWFS_SENS_CTX ctx = {
        .mapmode  = mapmode,
        .mode     = mode,
        .optmode  = optmode,
        .NBprobes = NBprobes,
    };

    aosim_ensure_wdir();
    NBprobesG = NBprobes;
    fpwfs_sens_read_vars(&ctx);
    fpwfs_sens_init_probes(&ctx);

    printf("noise mode  = %d (%f)\n", ctx.mode, ctx.probe_noise_prop);
    printf("mapmode     = %d\n", ctx.mapmode);
    printf("mapsize     = %ld   %ld   %ld\n", ctx.mapxsize, ctx.mapysize, ctx.mapzsize);
    printf("execmode    = %d\n", ctx.execmode);
    printf("mapzsize    = %ld\n", ctx.mapzsize);
    printf("imSimulMode = %d\n\n\n", ctx.imSimulMode);

    create_2Dimage_ID("psfprobeamp", ctx.mapxsize, ctx.mapysize, &ctx.IDprobampC);
    for (uint64_t ii = 0; ii < (uint64_t) (ctx.mapxsize * ctx.mapysize); ii++)
    {
        dcimg[ctx.IDprobampC].array.F[ii] = 1.0f;
    }

    NBoptVar              = (optmode == 2) ? 6 : ((optmode == 1) ? 4 : 3);
    ctx.xvect             = gsl_vector_alloc((size_t) NBoptVar);
    ctx.ss                = gsl_vector_alloc((size_t) NBoptVar);
    ctx.feval_func.n      = (size_t) NBoptVar;
    ctx.feval_func.f      = &f_eval;
    ctx.feval_func.params = NULL;
    ctx.s                 = gsl_multimin_fminimizer_alloc(
        gsl_multimin_fminimizer_nmsimplex2,
        (size_t) NBoptVar);

    for (long mapz = 0; mapz < ctx.mapzsize; mapz++)
    {
        if (ctx.imSimulMode == 1)
        {
            fpwfs_sens_run_imsimul(&ctx);
        }
        else if (ctx.mapmode == 2)
        {
            ctx.mapxsize = dcimg[ctx.IDpsfC].md[0].size[0];
            ctx.mapysize = dcimg[ctx.IDpsfC].md[0].size[1];
        }

        if (ctx.initmap == 0)
        {
            fpwfs_sens_init_maps(&ctx);
        }
        if (ctx.mapmode == 0)
        {
            printf("Preparing optimization ... \n");
            fflush(stdout);
        }
        if (ctx.execmode > 0)
        {
            aosim_fpwfs_solve_slice(&ctx, mapz);
        }

        fpwfs_sens_save_and_cleanup(&ctx);
        if ((mapz > 2) && (ctx.mapmode > 0))
        {
            fpwfs_sens_stats_iflux(&ctx, mapz);
            if (ctx.imSimulMode == 0)
            {
                fpwfs_sens_stats_ca(&ctx, mapz);
            }
        }
    }

    gsl_multimin_fminimizer_free(ctx.s);
    gsl_vector_free(ctx.ss);
    gsl_vector_free(ctx.xvect);
    return 0;
}
