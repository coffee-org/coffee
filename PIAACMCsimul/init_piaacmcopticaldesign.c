/**
 * @file    init_piaacmcopticaldesign.c
 * @brief   Top-level initialization of PIAACMC and APLCMC optical design prescription
 */

#define DEBUGLOG

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "OpticsMaterials/OpticsMaterials.h"

#include "FocalPlaneMask/mkFPM_zonemap.h"
#include "LyotStop/mkSimpleLyotStop.h"
#include "PIAACMCsimul.h"
#include "PIAACMCsimul_conf_io.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"
#include "PIAAshape/PIAAshape_modes.h"
#include "PIAAshape/setup_PIAA_shapes.h"
#include "init_piaacmcopticaldesign.h"

static int loaded   = 0;
static int saveconf = 0;

/**
 * init_design_dm_surfaces - Attach or allocate 2D deformable mirror surface maps.
 * @size: 2D wavefront grid size [pix]
 *
 * Return: RETURN_SUCCESS or memory allocation error status.
 */
static errno_t init_design_dm_surfaces(long size)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("%d DM(s)", piaacmcopticaldesign.nbDM);

    for (long iDM = 0; iDM < piaacmcopticaldesign.nbDM; iDM++)
    {
        char imname[STRINGMAXLEN_IMGNAME];

        printf("DM # %ld\n", iDM);
        WRITE_IMAGENAME(imname, "wfcDM%ld", iDM);

        piaacmcopticaldesign.ID_DM[iDM] = image_ID(imname, dcimg, dcnimg);
        printf("ID = %ld", piaacmcopticaldesign.ID_DM[iDM]);

        if (piaacmcopticaldesign.ID_DM[iDM] == -1)
        {
            read_sharedmem_image(imname, dcimg, dcnimg);
            piaacmcopticaldesign.ID_DM[iDM] = image_ID(imname, dcimg, dcnimg);
        }

        if (piaacmcopticaldesign.ID_DM[iDM] == -1)
        {
            uint32_t *sizearray = (uint32_t *) malloc(sizeof(uint32_t) * 2);
            if (sizearray == NULL)
            {
                FUNC_RETURN_FAILURE("malloc returns NULL pointer");
            }
            sizearray[0] = size;
            sizearray[1] = size;
            FUNC_CHECK_RETURN(create_image_ID(imname,
                                              2,
                                              sizearray,
                                              _DATATYPE_FLOAT,
                                              1,
                                              0,
                                              0,
                                              &(piaacmcopticaldesign.ID_DM[iDM])));
            free(sizearray);
        }
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * init_design_fpm_zonemap_and_zt - Ensure FPM zone map `fpmzmap` and zone sags `fpmzt` exist.
 * @flags: Initialization mode bitmask
 *
 * Return: RETURN_SUCCESS or FPM zone map error status.
 */
static errno_t init_design_fpm_zonemap_and_zt(uint64_t flags)
{
    DEBUG_TRACE_FSTART();

    int is_phys = (int) (flags & INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL);

    PIAACMCsimul_update_fnamedescr();
    EXECUTE_SYSTEM_COMMAND("echo \"%s\" > fpm_name.txt", piaacmcparams.fnamedescr);
    EXECUTE_SYSTEM_COMMAND("echo \"%s\" > fpm_name_conf.txt", piaacmcparams.fnamedescr_conf);

    if (image_ID("fpmzmap", dcimg, dcnimg) == -1)
    {
        printf("Make zonemap ...\n");
        fflush(stdout);
        FUNC_CHECK_RETURN(mkFPM_zonemap("fpmzmap", NULL));
    }
    else
    {
        printf("zonemap already exists\n");
        fflush(stdout);
    }

    piaacmcparams.CREATE_fpmzt = 0;
    if (piaacmcparams.FORCE_CREATE_fpmzt == 0)
    {
        piaacmcopticaldesign.zonezID = image_ID("fpmzt", dcimg, dcnimg);
        if (piaacmcopticaldesign.zonezID == -1)
        {
            char fname[STRINGMAXLEN_FULLFILENAME];

            PIAACMCsimul_update_fnamedescr();
            SNPRINTF_CHECK(fname,
                           STRINGMAXLEN_FULLFILENAME,
                           "%s/fpm_zonez.%s.fits",
                           piaacmcparams.piaacmcconfdir,
                           piaacmcparams.fnamedescr);
            printf("LOADING FILE NAME : \"%s\"  -  %d %d \n", fname, is_phys, loaded);
            load_fits(fname, "fpmzt", LOADFITS_ERRMODE_WARNING, &(piaacmcopticaldesign.zonezID));
            if (piaacmcopticaldesign.zonezID == -1)
            {
                piaacmcparams.CREATE_fpmzt = 1;
            }
        }
    }
    else
    {
        piaacmcparams.CREATE_fpmzt = 1;
    }

    if (piaacmcparams.CREATE_fpmzt == 1)
    {
        printf("Creating fpmzt, saving as fpm_zonez.fits - %d %d\n", is_phys, loaded);
        fflush(stdout);
        piaacmcopticaldesign.zonezID = image_ID("fpmzt", dcimg, dcnimg);
        if (piaacmcopticaldesign.zonezID != -1)
        {
            delete_image_ID("fpmzt", DELETE_IMAGE_ERRMODE_WARNING);
        }

        create_2Dimage_ID_double("fpmzt",
                                 piaacmcopticaldesign.focmNBzone,
                                 1,
                                 &(piaacmcopticaldesign.zonezID));
        double t = 1.0e-9;
        if (!is_phys)
        {
            double lam_mid = 0.5 * (piaacmcparams.LAMBDASTART + piaacmcparams.LAMBDAEND);
            double t0      = 1.0e-8;
            double pha0    = OpticsMaterials_pha_lambda(piaacmcopticaldesign.fpmmaterial_code,
                                                     t0,
                                                     lam_mid);
            t              = (M_PI / pha0) * t0;
            printf("t = %g m (%lf %g) -> %g %g\n",
                   t,
                   pha0,
                   t0,
                   OpticsMaterials_pha_lambda(piaacmcopticaldesign.fpmmaterial_code, t, lam_mid),
                   lam_mid);
            fflush(stdout);
            printf(" -- lambda = %g\n", piaacmcopticaldesign.lambda);
            printf(" -- lambdaB = %g\n", piaacmcopticaldesign.lambdaB);
            printf(" -- LAMBDASTART = %g\n", piaacmcparams.LAMBDASTART);
            printf(" -- LAMBDAEND = %g\n", piaacmcparams.LAMBDAEND);
        }

        for (long ii = 0; ii < piaacmcopticaldesign.focmNBzone; ii++)
        {
            dcimg[piaacmcopticaldesign.zonezID].array.D[ii] = t;
        }

        char fname[STRINGMAXLEN_FULLFILENAME];
        PIAACMCsimul_update_fnamedescr();
        SNPRINTF_CHECK(fname,
                       STRINGMAXLEN_FULLFILENAME,
                       "%s/fpm_zonez.%s.fits",
                       piaacmcparams.piaacmcconfdir,
                       piaacmcparams.fnamedescr);
        printf("Writing %s\n", fname);
        save_fits("fpmzt", fname);
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * init_design_fpm_zonea - Load or initialize FPM zone amplitude vector `fpmza`.
 *
 * Return: RETURN_SUCCESS or FITS write error status.
 */
static errno_t init_design_fpm_zonea(void)
{
    DEBUG_TRACE_FSTART();

    printf("CREATE_fpmza = %d\n", piaacmcparams.CREATE_fpmza);
    fflush(stdout);

    if (piaacmcparams.FORCE_CREATE_fpmza == 0)
    {
        piaacmcopticaldesign.zoneaID = image_ID("fpmza", dcimg, dcnimg);
        if (piaacmcopticaldesign.zoneaID == -1)
        {
            char fname[STRINGMAXLEN_FULLFILENAME];

            PIAACMCsimul_update_fnamedescr();
            SNPRINTF_CHECK(fname,
                           STRINGMAXLEN_FULLFILENAME,
                           "%s/fpm_zonea.%s.fits",
                           piaacmcparams.piaacmcconfdir,
                           piaacmcparams.fnamedescr);
            printf("LOADING FILE NAME : \"%s\"\n", fname);
            load_fits(fname, "fpmza", LOADFITS_ERRMODE_IGNORE, &(piaacmcopticaldesign.zoneaID));
            if (piaacmcopticaldesign.zoneaID == -1)
            {
                piaacmcparams.CREATE_fpmza = 1;
            }
        }
    }
    else
    {
        piaacmcparams.CREATE_fpmza = 1;
    }

    if (piaacmcparams.CREATE_fpmza == 1)
    {
        if (piaacmcopticaldesign.zoneaID != -1)
        {
            delete_image_ID("fpmza", DELETE_IMAGE_ERRMODE_WARNING);
        }
        create_2Dimage_ID_double("fpmza",
                                 piaacmcopticaldesign.focmNBzone,
                                 1,
                                 &(piaacmcopticaldesign.zoneaID));

        double amp_val = 1.0;
        if (piaacmcparams.PIAACMC_MASKRADLD > 0.2)
        {
            printf("PHYSICAL MASK ... %ld zones\n", piaacmcopticaldesign.focmNBzone);
            fflush(stdout);
        }
        else
        {
            printf("IDEALIZED MASK ... %ld zones\n", piaacmcopticaldesign.focmNBzone);
            amp_val = piaacmcopticaldesign.fpmaskamptransm;
        }

        for (long ii = 0; ii < piaacmcopticaldesign.focmNBzone; ii++)
        {
            dcimg[piaacmcopticaldesign.zoneaID].array.D[ii] = amp_val;
        }

        char fname[STRINGMAXLEN_FULLFILENAME];
        PIAACMCsimul_update_fnamedescr();
        SNPRINTF_CHECK(fname,
                       STRINGMAXLEN_FULLFILENAME,
                       "%s/fpm_zonea.%s.fits",
                       piaacmcparams.piaacmcconfdir,
                       piaacmcparams.fnamedescr);
        printf("Writing %s\n", fname);
        FUNC_CHECK_RETURN(save_fits("fpmza", fname));
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * init_design_lyot_stops - Load or generate initial Lyot stop masks for PIAACMC or APLCMC.
 * @xsize:      2D grid X dimension [pix]
 * @ysize:      2D grid Y dimension [pix]
 * @beamradpix: Beam radius [pix]
 *
 * Return: RETURN_SUCCESS or image/FITS error status.
 */
static errno_t init_design_lyot_stops(
    uint32_t xsize,
    uint32_t ysize,
    float    beamradpix)
{
    DEBUG_TRACE_FSTART();

    uint64_t xysize = (uint64_t) xsize * ysize;

    printf("LOADING/CREATING LYOT MASK  - %ld masks  (PIAAmode = %d, %ld x %ld)\n",
           piaacmcopticaldesign.NBLyotStop,
           piaacmcopticaldesign.PIAAmode,
           (long) xsize,
           (long) ysize);

    if (piaacmcopticaldesign.PIAAmode == 1)
    {
        for (long i = 0; i < piaacmcopticaldesign.NBLyotStop; i++)
        {
            char fname[STRINGMAXLEN_FULLFILENAME];
            char name[STRINGMAXLEN_IMGNAME];

            printf("LYOT MASK %ld\n", i);
            fflush(stdout);
            WRITE_FULLFILENAME(fname, "%s/LyotStop%ld.fits", piaacmcparams.piaacmcconfdir, i);
            WRITE_IMAGENAME(name, "lyotstop%ld", i);

            piaacmcopticaldesign.IDLyotStop[i] = image_ID(name, dcimg, dcnimg);
            if (piaacmcopticaldesign.IDLyotStop[i] == -1)
            {
                double r_in  = piaacmcopticaldesign.centObs1 + 0.02;
                double r_out = 0.98;
                if ((i == 0) && (piaacmcopticaldesign.NBLyotStop > 1))
                {
                    r_in = -0.01;
                }
                else if ((i == 1) && (piaacmcopticaldesign.NBLyotStop > 2))
                {
                    r_out = 1.2;
                }
                FUNC_CHECK_RETURN(
                    mkSimpleLyotStop(name, r_in, r_out, &(piaacmcopticaldesign.IDLyotStop[i])));
                save_fl_fits(name, fname);
            }
        }
    }
    else
    {
        long i = 0;
        char fname[STRINGMAXLEN_FULLFILENAME];
        char name[STRINGMAXLEN_IMGNAME];

        WRITE_FULLFILENAME(fname, "%s/LyotStop%ld.fits", piaacmcparams.piaacmcconfdir, i);
        WRITE_IMAGENAME(name, "lyotstop%ld", i);

        piaacmcopticaldesign.IDLyotStop[i] = image_ID(name, dcimg, dcnimg);
        if (piaacmcopticaldesign.IDLyotStop[i] == -1 && access(fname, F_OK) == 0)
        {
            load_fits(fname,
                      name,
                      LOADFITS_ERRMODE_IGNORE,
                      &(piaacmcopticaldesign.IDLyotStop[i]));
        }
        if (piaacmcopticaldesign.IDLyotStop[i] == -1)
        {
            FUNC_CHECK_RETURN(
                create_2Dimage_ID(name, xsize, ysize, &(piaacmcopticaldesign.IDLyotStop[i])));
            imageID ID = image_ID("pupmaskim", dcimg, dcnimg);
            for (uint64_t ii = 0; ii < xysize; ii++)
            {
                float pval = (ID != -1) ? dcimg[ID].array.F[ii] : 1.0f;
                dcimg[piaacmcopticaldesign.IDLyotStop[i]].array.F[ii] =
                    (pval < 0.99999999f) ? 0.0f : 1.0f;
            }

            for (uint32_t ii = 0; ii < xsize; ii++)
            {
                for (uint32_t jj = 0; jj < ysize; jj++)
                {
                    double x   = 1.0 * ii - 0.5 * xsize;
                    double y   = 1.0 * jj - 0.5 * ysize;
                    double rad = sqrt(x * x + y * y) / beamradpix;
                    if (rad < (piaacmcopticaldesign.centObs1 + 0.5 / beamradpix) ||
                        rad > (1.0 - 0.5 / beamradpix))
                    {
                        dcimg[piaacmcopticaldesign.IDLyotStop[i]].array.F[jj * xsize + ii] = 0.0;
                    }
                }
            }
            FUNC_CHECK_RETURN(save_fl_fits(name, fname));
        }
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * init_piaacmcopticaldesign - Initialize or load PIAACMC/APLCMC optical design prescription.
 * @fpmradld: Focal plane mask nominal radius [lambda/D]
 * @centobs0: Input pupil central obscuration ratio [beamrad]
 * @centobs1: Output pupil central obscuration ratio [beamrad]
 * @flags:    Bitmask of `INIT_PIAACMCOPTICALDESIGN_MODE__*` flags
 * @status:   Optional status pointer (unused)
 *
 * Return: RETURN_SUCCESS on success, or error code on failure.
 */
errno_t init_piaacmcopticaldesign(
    double    fpmradld,
    double    centobs0,
    double    centobs1,
    uint64_t  flags,
    uint64_t *status)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %f %f %f %lu", fpmradld, centobs0, centobs1, flags);

    (void) status;
    (void) mkdir("testdir", 0777);
    (void) mkdir("conf", 0777);
    (void) mkdir("status", 0777);
    (void) mkdir("ref", 0777);
    (void) mkdir("log", 0777);

    if (flags & INIT_PIAACMCOPTICALDESIGN_MODE__READCONF)
    {
        FUNC_CHECK_RETURN(PIAACMCsimul_initpiaacmcconf_readconfparams(
            (int) (flags & INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL),
            fpmradld,
            centobs0,
            centobs1,
            (int) (flags & INIT_PIAACMCOPTICALDESIGN_MODE__WSCMODE)));
    }

    double lam_mid = 0.5 * (piaacmcparams.LAMBDASTART + piaacmcparams.LAMBDAEND);
    piaacmcopticaldesign.fpmCentConeRad = piaacmcopticaldesign.fpmRad *
                                          piaacmcopticaldesign.NBringCentCone /
                                          piaacmcopticaldesign.NBrings;

    printf("fpmRad = %g m\n", lam_mid * piaacmcopticaldesign.Fratio * fpmradld);
    printf("factor = %f   (%ld / %ld)\n",
           1.0 * piaacmcopticaldesign.NBringCentCone / piaacmcopticaldesign.NBrings,
           piaacmcopticaldesign.NBringCentCone,
           piaacmcopticaldesign.NBrings);
    printf("fpmCentConeRad =  %g\n",
           (lam_mid * piaacmcopticaldesign.Fratio * fpmradld) *
               piaacmcopticaldesign.NBringCentCone / piaacmcopticaldesign.NBrings);

    if (flags & INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF)
    {
        EXECUTE_SYSTEM_COMMAND("mkdir -p %s", piaacmcparams.piaacmcconfdir);
        loaded = PIAACMCsimul_loadpiaacmcconf(piaacmcparams.piaacmcconfdir);
        if (loaded == 0)
        {
            saveconf = 1;
        }
    }

    FUNC_CHECK_RETURN(PIAACMCsimul_read_piaa_material());

    printf("lambda = %g\n", piaacmcopticaldesign.lambda);
    printf("LAMBDASTART = %g\n", piaacmcparams.LAMBDASTART);
    printf("LAMBDAEND = %g\n", piaacmcparams.LAMBDAEND);

    for (int k = 0; k < piaacmcopticaldesign.nblambda; k++)
    {
        piaacmcopticaldesign.lambdaarray[k] =
            piaacmcparams.LAMBDASTART +
            (0.5 + k) * (piaacmcparams.LAMBDAEND - piaacmcparams.LAMBDASTART) /
                piaacmcopticaldesign.nblambda;
    }

    float beamradpix = piaacmcopticaldesign.beamrad / piaacmcopticaldesign.pixscale;
    long  size       = piaacmcopticaldesign.size;

    make_x_y_r_PA_images(size, beamradpix);
    FUNC_CHECK_RETURN(init_design_dm_surfaces(size));
    make_C_F_modes(size, beamradpix);
    FUNC_CHECK_RETURN(setupPIAAshapes((int) (flags & INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL),
                                      (uint32_t) size,
                                      beamradpix,
                                      loaded,
                                      &saveconf));

    FUNC_CHECK_RETURN(init_design_fpm_zonemap_and_zt(flags));
    FUNC_CHECK_RETURN(init_design_fpm_zonea());
    FUNC_CHECK_RETURN(init_design_lyot_stops((uint32_t) size, (uint32_t) size, beamradpix));

    if (saveconf == 1)
    {
        FUNC_CHECK_RETURN(PIAACMCsimul_savepiaacmcconf(piaacmcparams.piaacmcconfdir));
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
