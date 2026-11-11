/**
 * @file    exec_multizone_fpm_calib_zones.c
 * @brief   Worker zone propagation and multi-process FITS combination for Mode 11
 */

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul.h"
#include "PIAACMCsimul_computePSF.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"

#include "FocalPlaneMask/exec_multizone_fpm_calib_internal.h"

/**
 * fpm_calib_wait_child_file - Poll until worker thread thr completes its partial FPMresp FITS
 * @fnamet: Path to FPMthreadstatus.txt diagnostic log file.
 * @thr:    Worker thread index to await.
 * @outID1: Output pointer receiving loaded tmpFPMresp image identifier.
 *
 * Return: RETURN_SUCCESS once partial FITS is loaded, error code otherwise.
 */
static errno_t fpm_calib_wait_child_file(
    const char *fnamet,
    long        thr,
    imageID    *outID1)
{
    imageID ID1 = -1;

    while (ID1 == -1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];

        PIAACMCsimul_update_fnamedescr_conf();
        SNPRINTF_CHECK(
            fname,
            STRINGMAXLEN_FULLFILENAME,
            "%s/FPMresp%d.%s.mp%02ld_thread%02ld.fits",
            piaacmcparams.piaacmcconfdir,
            piaacmcparams.SCORINGMASKTYPE,
            piaacmcparams.fnamedescr_conf,
            piaacmcparams.PIAACMC_FPMresp_mp,
            thr);

        printf("Waiting for file \"%s\" ...\n", fname);
        fflush(stdout);

        FILE *fpt = fopen(fnamet, "a");
        if (fpt != NULL)
        {
            fprintf(
                fpt,
                "Process %ld (thread %ld) --- Waiting for file \"%s\" ...\n",
                (long) getpid(),
                thr,
                fname);
            fclose(fpt);
        }
        sleep(1);

        FUNC_CHECK_RETURN(delete_image_ID("tmpFPMresp", DELETE_IMAGE_ERRMODE_WARNING));
        FUNC_CHECK_RETURN(load_fits(fname, "tmpFPMresp", LOADFITS_ERRMODE_WARNING, &ID1));
    }

    *outID1 = ID1;
    return RETURN_SUCCESS;
}

/**
 * fpm_calib_merge_child_zones - Accumulate worker thread thr zones into combined FPMresp cube
 * @IDfpmresp: Destination combined FPMresp 3D double image identifier.
 * @ID1:       Loaded partial FPMresp image identifier from worker thr.
 * @thr:       Worker thread index defining zone stride offset.
 */
static void fpm_calib_merge_child_zones(
    imageID IDfpmresp,
    imageID ID1,
    long    thr)
{
    long     mzstep  = piaacmcparams.PIAACMC_FPMresp_mp;
    uint32_t neval   = dcimg[IDfpmresp].md[0].size[0];
    uint32_t nzones1 = dcimg[piaacmcopticaldesign.zonezID].md[0].size[0] + 1;

    for (long mz = 1 + thr; mz < piaacmcopticaldesign.focmNBzone + 1; mz += mzstep)
    {
        printf("mz = %ld    %ld %ld\n", mz, IDfpmresp, ID1);
        fflush(stdout);

        for (int k = 0; k < piaacmcopticaldesign.nblambda; k++)
        {
            size_t kbase = (size_t) k * nzones1 * neval;
            for (uint32_t ii = 0; ii < neval; ii++)
            {
                size_t tmpl1 = kbase + (size_t) mz * neval + ii;

                dcimg[IDfpmresp].array.D[tmpl1]      = dcimg[ID1].array.D[tmpl1];
                dcimg[IDfpmresp].array.D[kbase + ii] -= dcimg[ID1].array.D[tmpl1];
            }
        }
    }
}

/**
 * fpm_calib_combine_threads - Signal workers, gather partial FPMresp cubes, and write final FITS
 * @IDfpmresp: Combined FPMresp cube initialized with unobstructed PSF in zone 0.
 *
 * Return: RETURN_SUCCESS on full assembly and disk write, error code otherwise.
 */
errno_t fpm_calib_combine_threads(imageID IDfpmresp)
{
    long       index = 0;
    variableID IDv;

    if ((IDv = variable_ID("PID")) != -1)
    {
        index = (long) coffee_dcvar_f(IDv) + 0.01;
    }
    EXECUTE_SYSTEM_COMMAND("touch initOK_%ld", index);

    printf("COMBINING FILES\n");
    fflush(stdout);

    char fnamet[STRINGMAXLEN_FULLFILENAME];
    WRITE_FULLFILENAME(fnamet, "%s/FPMthreadstatus.txt", piaacmcparams.piaacmcconfdir);
    FILE *fpt = fopen(fnamet, "w");
    if (fpt != NULL)
    {
        fclose(fpt);
    }

    for (long thr = 0; thr < piaacmcparams.PIAACMC_FPMresp_mp; thr++)
    {
        imageID ID1 = -1;

        printf("thr = %ld\n", thr);
        fflush(stdout);

        FUNC_CHECK_RETURN(fpm_calib_wait_child_file(fnamet, thr, &ID1));

        fpt = fopen(fnamet, "a");
        if (fpt != NULL)
        {
            fprintf(fpt, "READING %s\n", fnamet);
            fclose(fpt);
        }

        fpm_calib_merge_child_zones(IDfpmresp, ID1, thr);
        FUNC_CHECK_RETURN(delete_image_ID("tmpFPMresp", DELETE_IMAGE_ERRMODE_WARNING));
    }

    char fname[STRINGMAXLEN_FULLFILENAME];
    PIAACMCsimul_update_fnamedescr_conf();
    SNPRINTF_CHECK(
        fname,
        STRINGMAXLEN_FULLFILENAME,
        "%s/FPMresp%d.%s.fits",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.SCORINGMASKTYPE,
        piaacmcparams.fnamedescr_conf);
    save_fits("FPMresp", fname);

    EXECUTE_SYSTEM_COMMAND("rm -f %s/FPMresp*.fits.tmp", piaacmcparams.piaacmcconfdir);
    return RETURN_SUCCESS;
}

/**
 * fpm_calib_find_pre_fpm_elem - Locate optical element index immediately preceding FPM pupil
 *
 * Return: Element index before "post focal plane mask pupil", or 0 if not found.
 */
static long fpm_calib_find_pre_fpm_elem(void)
{
    for (long elem = 1; elem < piaacmcopticalsystem.NBelem; elem++)
    {
        if (strcmp("post focal plane mask pupil", piaacmcopticalsystem.name[elem]) == 0)
        {
            return elem - 1;
        }
    }

    return 0;
}

/**
 * fpm_calib_eval_single_zone - Propagate FPM impulse response for zone mz and checkpoint to tmp
 * @IDfpmresp: Target FPMresp 3D double image identifier.
 * @mz:        1-based FPM zone index to evaluate.
 * @fname1:    Temporary FITS path for checkpoint writes.
 *
 * Return: RETURN_SUCCESS on PSF propagation and checkpoint save, error code otherwise.
 */
static errno_t fpm_calib_eval_single_zone(
    imageID     IDfpmresp,
    long        mz,
    const char *fname1)
{
    piaacmcparams.focmMode                    = mz;
    piaacmcopticalsystem.FOCMASKarray[0].mode = 0;

    long elem0                          = fpm_calib_find_pre_fpm_elem();
    piaacmcopticalsystem.keepMem[elem0] = 1;

    printf("piaacmcopticaldesign.NBrings =                             %ld\n",
           piaacmcopticaldesign.NBrings);
    printf("piaacmcopticaldesign.focmNBzone =                          %ld\n",
           piaacmcopticaldesign.focmNBzone);
    printf("piaacmcopticaldesign.nblambda =                            %d\n",
           piaacmcopticaldesign.nblambda);
    printf("dcimg[IDfpmresp].md[0].size[0] =            %ld\n",
           (long) dcimg[IDfpmresp].md[0].size[0]);
    printf("dcimg[piaacmcopticaldesign.zonezID].md[0].size[0] =   %ld\n",
           (long) dcimg[piaacmcopticaldesign.zonezID].md[0].size[0]);
    fflush(stdout);

    FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(
        0.0,
        0.0,
        elem0,
        piaacmcopticalsystem.NBelem,
        0,
        piaacmcparams.computePSF_ResolvedTarget,
        piaacmcparams.computePSF_ResolvedTarget_mode,
        0,
        NULL));

    imageID  IDimvect = image_ID("imvect", dcimg, dcnimg);
    uint32_t neval    = dcimg[IDimvect].md[0].size[0];
    uint32_t nzones1  = dcimg[piaacmcopticaldesign.zonezID].md[0].size[0] + 1;

    for (int k = 0; k < piaacmcopticaldesign.nblambda; k++)
    {
        size_t kbase = (size_t) k * nzones1 * neval;
        for (uint32_t ii = 0; ii < neval; ii++)
        {
            float val = dcimg[IDimvect].array.F[(size_t) k * neval + ii];

            dcimg[IDfpmresp].array.D[kbase + (size_t) mz * neval + ii] = val;
            if (piaacmcparams.PIAACMC_FPMresp_mp == 1)
            {
                dcimg[IDfpmresp].array.D[kbase + ii] -= val;
            }
        }
    }

    printf("Saving FPMresp (ID = %ld) as \"%s\" ...", image_ID("FPMresp", dcimg, dcnimg), fname1);
    fflush(stdout);
    FUNC_CHECK_RETURN(save_fits("FPMresp", fname1));
    printf("Done \n");
    fflush(stdout);

    return RETURN_SUCCESS;
}

/**
 * fpm_calib_compute_worker_zones - Evaluate assigned zone stride and atomically publish FITS
 * @ctx:       Calibration context specifying zone stride and output file paths.
 * @IDfpmresp: Target FPMresp 3D double image identifier.
 *
 * Return: RETURN_SUCCESS on completion and atomic rename, error code otherwise.
 */
errno_t fpm_calib_compute_worker_zones(
    const struct fpm_calib_ctx *ctx,
    imageID                     IDfpmresp)
{
    char fnametmp[STRINGMAXLEN_FULLFILENAME];

    WRITE_FULLFILENAME(
        fnametmp,
        "%s/fpmzmap_thread%02ld.fits",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.PIAACMC_FPMresp_thread);
    FUNC_CHECK_RETURN(save_fits("fpmzmap", fnametmp));

    printf(
        "Making component %ld / %ld\n",
        piaacmcparams.PIAACMC_FPMresp_thread,
        piaacmcparams.PIAACMC_FPMresp_mp);
    fflush(stdout);

    piaacmcparams.WRITE_OK = 0;
    for (long mz = 1 + ctx->mzoffset; mz < piaacmcopticaldesign.focmNBzone + 1; mz += ctx->mzstep)
    {
        FUNC_CHECK_RETURN(fpm_calib_eval_single_zone(IDfpmresp, mz, ctx->fname1));
    }

    printf(
        "Saving FPMresp (ID = %ld) as \"%s\" ...",
        image_ID("FPMresp", dcimg, dcnimg),
        ctx->fname2);
    fflush(stdout);
    FUNC_CHECK_RETURN(save_fits("FPMresp", ctx->fname1));
    if (rename(ctx->fname1, ctx->fname2) != 0)
    {
        perror("rename");
        return RETURN_FAILURE;
    }
    printf("Done \n");
    fflush(stdout);

    piaacmcparams.WRITE_OK = 0;
    return RETURN_SUCCESS;
}
