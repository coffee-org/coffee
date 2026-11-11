/**
 * @file    exec_optimize_lyot_stops_shapes_positions.c
 * @brief   Jointly optimize multi-plane Lyot stop shapes and axial conjugations
 */

#define DEBUGLOG

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"
#include "PIAACMCsimul_CA2propCubeInt.h"
#include "PIAACMCsimul_computePSF.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"
#include "init_piaacmcopticaldesign.h"
#include "init_piaacmcopticalsystem.h"

#include "LyotStop/exec_optimize_lyot_stops_shapes_positions.h"
#include "LyotStop/optimizeLyotStop.h"
#include "PIAAshape/makePIAAshapes.h"

struct mode5_prop_cfg
{
    char     imnamea[STRINGMAXLEN_IMGNAME];
    char     imnamep[STRINGMAXLEN_IMGNAME];
    float    zmin;
    float    zmax;
    long     NBpropstep;
    uint32_t xsize;
    uint32_t ysize;
    uint64_t xysize;
};

/**
 * optimizeLyotStop_offaxis_min - Find minimum flux level along z across a 3D intensity cube.
 * @IDincohc_name: Input 3D incoherent intensity cube name
 *
 * Writes `test_oals_index.fits` and `test_oals_val.fits` to the working directory.
 *
 * Return: RETURN_SUCCESS or failure if input image is missing or FITS write fails.
 */
errno_t optimizeLyotStop_offaxis_min(const char *restrict IDincohc_name)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %s", IDincohc_name);

    imageID IDincohc = image_ID(IDincohc_name, dcimg, dcnimg);
    if (IDincohc == -1)
    {
        FUNC_RETURN_FAILURE("Image %s not found in memory", IDincohc_name);
    }

    uint32_t xsize  = dcimg[IDincohc].md[0].size[0];
    uint32_t ysize  = dcimg[IDincohc].md[0].size[1];
    uint64_t xysize = (uint64_t) xsize * ysize;
    uint32_t NBz    = dcimg[IDincohc].md[0].size[2];

    imageID IDindex;
    FUNC_CHECK_RETURN(create_2Dimage_ID("oals_index", xsize, ysize, &IDindex));

    imageID IDminflux;
    FUNC_CHECK_RETURN(create_2Dimage_ID("oals_val", xsize, ysize, &IDminflux));

    DEBUG_TRACEPOINT("scanning for minimum");
    for (uint64_t ii = 0; ii < xysize; ii++)
    {
        float    minv     = dcimg[IDincohc].array.F[ii];
        uint32_t minindex = 0;

        for (uint32_t kk = 1; kk < NBz; kk++)
        {
            float tmpv = dcimg[IDincohc].array.F[xysize * kk + ii];
            if (tmpv < minv)
            {
                minv     = tmpv;
                minindex = kk;
            }
        }
        dcimg[IDminflux].array.F[ii] = minv;
        dcimg[IDindex].array.F[ii]   = (float) minindex;
    }

    DEBUG_TRACEPOINT("Saving minimum map to filesystem");
    FUNC_CHECK_RETURN(save_fits("oals_index", "test_oals_index.fits"));
    FUNC_CHECK_RETURN(save_fits("oals_val", "test_oals_val.fits"));

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * load_mode5_cli_and_design - Read CLI parameters and initialize optical design and PIAA shapes.
 * @nbpropstep_out: Output number of axial Fresnel propagation planes
 * @lstransm_out:   Output target Lyot stop geometric transmission
 *
 * Return: RETURN_SUCCESS or initialization error status.
 */
static errno_t load_mode5_cli_and_design(
    long   *nbpropstep_out,
    double *lstransm_out)
{
    imageID IDv;
    double  fpmradld = 0.95;
    double  centobs0 = 0.3;
    double  centobs1 = 0.2;

    printf("=================================== mode 005 ===================================\n");
    if ((IDv = variable_ID("PIAACMC_centobs0")) != -1)
    {
        centobs0 = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_centobs1")) != -1)
    {
        centobs1 = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = coffee_dcvar_f(IDv);
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }

    uint64_t initflag = INIT_PIAACMCOPTICALDESIGN_MODE__READCONF;
    initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF;
    FUNC_CHECK_RETURN(
        init_piaacmcopticaldesign(fpmradld, centobs0, centobs1, initflag, NULL));

    FUNC_CHECK_RETURN(makePIAAshapes());
    piaacmcopticalsystem.FOCMASKarray[0].mode = 1;

    long NBpropstep = 150;
    if ((IDv = variable_ID("PIAACMC_nbpropstep")) != -1)
    {
        NBpropstep = (long) coffee_dcvar_f(IDv) + 0.01;
    }

    double lstransm = 0.85;
    if ((IDv = variable_ID("PIAACMC_lstransm")) != -1)
    {
        lstransm = (double) coffee_dcvar_f(IDv);
    }
    printf("lstransm  = %f\n", lstransm);

    FUNC_CHECK_RETURN(init_piaacmcopticalsystem(0.0, 0.0));

    *nbpropstep_out = NBpropstep;
    *lstransm_out   = lstransm;
    return RETURN_SUCCESS;
}

/**
 * find_post_fpm_pupil_elem - Locate the post focal plane mask pupil element and mark keepMem.
 *
 * Return: Index of the post-FPM pupil element in `piaacmcopticalsystem`.
 */
static long find_post_fpm_pupil_elem(void)
{
    printf("=========== %ld elements ======================================================\n",
           piaacmcopticalsystem.NBelem);

    long elem0 = 0;
    for (long elem = 0; elem < piaacmcopticalsystem.NBelem; elem++)
    {
        if (strcmp("post focal plane mask pupil", piaacmcopticalsystem.name[elem]) == 0)
        {
            elem0 = elem;
            DEBUG_TRACEPOINT("post focal plane mask pupil = %ld", elem);
        }
        else
        {
            DEBUG_TRACEPOINT("elem %ld : %s", elem, piaacmcopticalsystem.name[elem]);
        }
    }
    piaacmcopticalsystem.keepMem[elem0] = 1;
    return elem0;
}

/**
 * accumulate_offaxis_point - Propagate one off-axis point source and add its 3D pupil intensity.
 * @cfg: Propagation range and post-FPM wavefront image configuration
 * @xld: Source X angular offset [lambda/D]
 * @yld: Source Y angular offset [lambda/D]
 * @IDc: Accumulator 3D image ID (`OAincohc`)
 *
 * Return: RETURN_SUCCESS or PSF/propagation error status.
 */
static errno_t accumulate_offaxis_point(
    const struct mode5_prop_cfg *cfg,
    double                       xld,
    double                       yld,
    imageID                      IDc)
{
    double cval = 0.0;
    FUNC_CHECK_RETURN(
        PIAACMCsimul_computePSF(xld, yld, 0, piaacmcopticalsystem.NBelem, 0, 0, 0, 0, &cval));

    imageID ID1;
    FUNC_CHECK_RETURN(PIAACMCsimul_CA2propCubeInt(cfg->imnamea,
                                                  cfg->imnamep,
                                                  cfg->zmin,
                                                  cfg->zmax,
                                                  cfg->NBpropstep,
                                                  "iproptmp",
                                                  &ID1));

    for (uint64_t ii = 0; ii < cfg->xysize; ii++)
    {
        for (long k = 0; k < cfg->NBpropstep; k++)
        {
            dcimg[IDc].array.F[cfg->xysize * k + ii] +=
                dcimg[ID1].array.F[cfg->xysize * k + ii];
        }
    }

    FUNC_CHECK_RETURN(delete_image_ID("iproptmp", DELETE_IMAGE_ERRMODE_WARNING));
    return RETURN_SUCCESS;
}

/**
 * ensure_offaxis_incoh_cube - Load cached `OAincohc.fits` or synthesize it over off-axis rings.
 * @cfg: Propagation range and post-FPM wavefront image configuration
 *
 * Return: RETURN_SUCCESS or propagation/FITS error status.
 */
static errno_t ensure_offaxis_incoh_cube(const struct mode5_prop_cfg *cfg)
{
    imageID IDc;
    char    fname[STRINGMAXLEN_FULLFILENAME];

    WRITE_FULLFILENAME(fname, "%s/OAincohc.fits", piaacmcparams.piaacmcconfdir);
    FUNC_CHECK_RETURN(load_fits(fname, "OAincohc", LOADFITS_ERRMODE_WARNING, &IDc));
    DEBUG_TRACEPOINT("IDc = %ld", IDc);

    if (IDc != -1)
    {
        return RETURN_SUCCESS;
    }

    DEBUG_TRACEPOINT("Creating off-axis light 3D volumetric cube");
    FUNC_CHECK_RETURN(
        create_3Dimage_ID("OAincohc", cfg->xsize, cfg->ysize, cfg->NBpropstep, &IDc));

    long   cnt      = 0;
    long   NBkr     = 2;
    long   NBincpt  = 9;
    double oaoffset = 20.0;

    for (long kr = 0; kr < NBkr; kr++)
    {
        for (long k1 = 0; k1 < NBincpt; k1++)
        {
            DEBUG_TRACEPOINT("OAincohc point %ld/%ld %ld/%ld", kr, NBkr, k1, NBincpt);
            double rad   = oaoffset * (1.0 + kr) / NBkr;
            double angle = 2.0 * M_PI * k1 / NBincpt;
            FUNC_CHECK_RETURN(
                accumulate_offaxis_point(cfg, rad * cos(angle), rad * sin(angle), IDc));
            cnt++;
        }
    }

    for (uint64_t ii = 0; ii < cfg->xysize; ii++)
    {
        for (long k = 0; k < cfg->NBpropstep; k++)
        {
            dcimg[IDc].array.F[k * cfg->xysize + ii] /= cnt;
        }
    }

    WRITE_FULLFILENAME(fname, "%s/OAincohc.fits", piaacmcparams.piaacmcconfdir);
    FUNC_CHECK_RETURN(save_fits("OAincohc", fname));
    return RETURN_SUCCESS;
}

/**
 * save_and_install_lyot_stops - Log optimal conjugations, convert to absolute z, and install FITS.
 * @elem0:       Reference post-FPM pupil element index
 * @zmin:        Minimum propagation range [m]
 * @zmax:        Maximum propagation range [m]
 * @caller_func: Name of calling entry function for header provenance comment
 *
 * Return: RETURN_SUCCESS or configuration save error status.
 */
static errno_t save_and_install_lyot_stops(
    long        elem0,
    float       zmin,
    float       zmax,
    const char *caller_func)
{
    char fptestname[STRINGMAXLEN_FILENAME];
    WRITE_FILENAME(fptestname, "conj_test.txt");
    DEBUG_TRACEPOINT("write %s", fptestname);

    FILE *fptest = fopen(fptestname, "w");
    if (fptest == NULL)
    {
        PRINT_ERROR("cannot create file %s", fptestname);
        abort();
    }
    fprintf(fptest, "# Optimal Lyot stop conjugations\n");
    fprintf(fptest, "# \n");
    fprintf(fptest, "# DO NOT EDIT THIS FILE\n");
    fprintf(fptest, "# Written by %s in %s\n", caller_func, __FILE__);
    fprintf(fptest, "# \n");
    fprintf(fptest, "# Lyot stop index   zmin   zmax    LyotStop_zpos    elemZpos[elem0]\n");
    for (long ls = 0; ls < piaacmcopticaldesign.NBLyotStop; ls++)
    {
        fprintf(fptest,
                "%5ld  %f  %f     %f  %f\n",
                ls,
                zmin,
                zmax,
                piaacmcopticaldesign.LyotStop_zpos[ls],
                piaacmcopticalsystem.elemZpos[elem0]);
    }
    fclose(fptest);

    for (long ls = 0; ls < piaacmcopticaldesign.NBLyotStop; ls++)
    {
        piaacmcopticaldesign.LyotStop_zpos[ls] += piaacmcopticalsystem.elemZpos[elem0];
    }

    FUNC_CHECK_RETURN(PIAACMCsimul_savepiaacmcconf(piaacmcparams.piaacmcconfdir));

    for (long ls = 0; ls < piaacmcopticaldesign.NBLyotStop; ls++)
    {
        EXECUTE_SYSTEM_COMMAND("cp ./%s/optLM%02ld.fits ./%s/LyotStop%ld.fits",
                               piaacmcparams.piaacmcconfdir,
                               ls,
                               piaacmcparams.piaacmcconfdir,
                               ls);
    }

    return RETURN_SUCCESS;
}

/**
 * exec_optimize_lyot_stops_shapes_positions - Mode 5 joint Lyot stop shape and conjugation solver.
 *
 * Computes 3D on-axis and off-axis Fresnel intensity cubes downstream of the focal plane mask,
 * synthesizes multi-plane binary Lyot stops, and updates configuration files.
 *
 * Return: RETURN_SUCCESS on completion, or milk error code on failure.
 */
errno_t exec_optimize_lyot_stops_shapes_positions(void)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG");

    struct mode5_prop_cfg cfg;
    double                lstransm = 0.85;

    FUNC_CHECK_RETURN(load_mode5_cli_and_design(&cfg.NBpropstep, &lstransm));
    long elem0 = find_post_fpm_pupil_elem();

    DEBUG_TRACEPOINT("compute the reference on-axis PSF");
    {
        double cval = 0.0;
        FUNC_CHECK_RETURN(
            PIAACMCsimul_computePSF(0.0, 0.0, 0, piaacmcopticalsystem.NBelem, 0, 0, 0, 0, &cval));
    }

    WRITE_IMAGENAME(cfg.imnamea, "WFamp0_%03ld", elem0);
    WRITE_IMAGENAME(cfg.imnamep, "WFpha0_%03ld", elem0);
    cfg.zmin = piaacmcopticaldesign.LyotZmin;
    cfg.zmax = piaacmcopticaldesign.LyotZmax;

    DEBUG_TRACEPOINT("propagate complex amplitude in a range from %f to %f", cfg.zmin, cfg.zmax);
    FUNC_CHECK_RETURN(PIAACMCsimul_CA2propCubeInt(cfg.imnamea,
                                                  cfg.imnamep,
                                                  cfg.zmin,
                                                  cfg.zmax,
                                                  cfg.NBpropstep,
                                                  "iproptmp",
                                                  NULL));

    imageID IDa = image_ID(cfg.imnamea, dcimg, dcnimg);
    cfg.xsize   = dcimg[IDa].md[0].size[0];
    cfg.ysize   = dcimg[IDa].md[0].size[1];
    cfg.xysize  = (uint64_t) cfg.xsize * cfg.ysize;

    FUNC_CHECK_RETURN(ensure_offaxis_incoh_cube(&cfg));

    DEBUG_TRACEPOINT("Compute on-axis PSF 3D intensity to define light to reject");
    FUNC_CHECK_RETURN(
        PIAACMCsimul_computePSF(0.0, 0.0, 0, piaacmcopticalsystem.NBelem, 0, 0, 0, 0, NULL));

    DEBUG_TRACEPOINT("compute on-axis 3D intensity cube");
    FUNC_CHECK_RETURN(PIAACMCsimul_CA2propCubeInt(cfg.imnamea,
                                                  cfg.imnamep,
                                                  cfg.zmin,
                                                  cfg.zmax,
                                                  cfg.NBpropstep,
                                                  "iprop00",
                                                  NULL));

    FUNC_CHECK_RETURN(optimizeLyotStop_offaxis_min("OAincohc"));

    DEBUG_TRACEPOINT("Optimize Lyot stop geometry");
    FUNC_CHECK_RETURN(optimizeLyotStop(cfg.imnamea,
                                       cfg.imnamep,
                                       "OAincohc",
                                       cfg.zmin,
                                       cfg.zmax,
                                       lstransm,
                                       cfg.NBpropstep,
                                       piaacmcopticaldesign.NBLyotStop,
                                       NULL));

    FUNC_CHECK_RETURN(save_and_install_lyot_stops(elem0, cfg.zmin, cfg.zmax, __FUNCTION__));

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
