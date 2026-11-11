/**
 * @file    PIAACMCsimul_loadsavepiaacmcconf.c
 * @brief   Load/save PIAACMC design parameters and FPM filename descriptors
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "image_basic/image_basic.h"

#include "PIAACMCsimul.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"

#include "FocalPlaneMask/mkFPM_zonemap.h"

/**
 * check_fscanf_count - Verify fscanf item match count and abort on malformed config input.
 * @fp:       Active configuration file stream
 * @actual:   Return value from fscanf
 * @expected: Expected number of matched conversions
 */
static void check_fscanf_count(
    FILE *fp,
    int   actual,
    int   expected)
{
    if (actual == EOF)
    {
        if (ferror(fp))
        {
            perror("fscanf");
        }
        else
        {
            fprintf(
                stderr,
                "Error: fscanf reached end of file, no matching characters, "
                "no matching failure\n");
        }
        abort();
    }
    if (actual != expected)
    {
        fprintf(
            stderr,
            "Error: fscanf successfully matched and assigned %i input items, "
            "%i expected\n",
            actual,
            expected);
        abort();
    }
}

/**
 * load_lyot_stops_conf - Load active LyotStop FITS masks and read 10 Lyot axial positions.
 * @fp:    Open piaacmcparams.conf stream
 * @dname: Configuration directory path
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t load_lyot_stops_conf(
    FILE       *fp,
    const char *dname)
{
    for (int i = 0; i < 10; i++)
    {
        if (i < piaacmcopticaldesign.NBLyotStop)
        {
            char fname[STRINGMAXLEN_FULLFILENAME];
            char imname[STRINGMAXLEN_IMGNAME];

            WRITE_FULLFILENAME(fname, "%s/LyotStop%d.fits", dname, i);
            WRITE_IMAGENAME(imname, "lyotstop%d", i);
            printf("Loading \"%s\" as \"%s\"\n", fname, imname);
            FUNC_CHECK_RETURN(
                load_fits(fname, imname, LOADFITS_ERRMODE_WARNING,
                          &(piaacmcopticaldesign.IDLyotStop[i])));
        }

        double tmplf;
        long   tmpl;
        int    fscanfcnt = fscanf(fp, "%lf   LyotStop_zpos %ld\n", &tmplf, &tmpl);
        check_fscanf_count(fp, fscanfcnt, 2);

        piaacmcopticaldesign.LyotStop_zpos[i] = tmplf;
        printf("LYOT STOP %d POS : %lf\n", i, piaacmcopticaldesign.LyotStop_zpos[i]);
    } // for i

    return RETURN_SUCCESS;
}

/**
 * load_piaa_mode_coeffs - Load a PIAA modal coefficient vector with fallback to piaaref/.
 * @dname:    Configuration directory path
 * @fitsname: Base FITS filename (e.g. "piaa0Cmodes.fits")
 * @imname:   In-memory image name
 * @outID:    Destination image ID pointer
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t load_piaa_mode_coeffs(
    const char *dname,
    const char *fitsname,
    const char *imname,
    long       *outID)
{
    char fname[STRINGMAXLEN_FULLFILENAME];

    WRITE_FULLFILENAME(fname, "%s/%s", dname, fitsname);
    FUNC_CHECK_RETURN(load_fits(fname, imname, 1, outID));
    if (*outID == -1)
    {
        WRITE_FULLFILENAME(fname, "%s/piaaref/%s", dname, fitsname);
        FUNC_CHECK_RETURN(load_fits(fname, imname, 1, outID));
    }

    return RETURN_SUCCESS;
}

/**
 * PIAACMCsimul_loadpiaacmcconf - Read piaacmcparams.conf, Lyot stops, and PIAA mode tables.
 * @dname: Configuration directory path
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
errno_t PIAACMCsimul_loadpiaacmcconf(const char *dname)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG dname %s", dname);

    char fname[STRINGMAXLEN_FULLFILENAME];
    WRITE_FULLFILENAME(fname, "%s/piaacmcparams.conf", dname);
    printf("%s\n", fname);

    FILE *fp = fopen(fname, "r");
    if (fp == NULL)
    {
        printf(
            "Configuration file \"%s\" does not exist (yet), using previously "
            "set configuration\n",
            fname);
        fflush(stdout);
        DEBUG_TRACE_FEXIT();
        return RETURN_SUCCESS;
    }

    long tmpl;
    int  fscanfcnt = fscanf(fp, "%ld   NBradpts\n", &tmpl);
    check_fscanf_count(fp, fscanfcnt, 1);
    piaacmcopticaldesign.NBradpts = tmpl;

    errno_t err = load_lyot_stops_conf(fp, dname);
    if (err != RETURN_SUCCESS)
    {
        fclose(fp);
        return err;
    }

    err = load_piaa_mode_coeffs(
        dname, "piaa0Cmodes.fits", "piaa0Cmodescoeff", &piaacmcopticaldesign.piaa0CmodesID);
    if (err == RETURN_SUCCESS)
    {
        err = load_piaa_mode_coeffs(
            dname, "piaa0Fmodes.fits", "piaa0Fmodescoeff", &piaacmcopticaldesign.piaa0FmodesID);
    }
    if (err == RETURN_SUCCESS)
    {
        err = load_piaa_mode_coeffs(
            dname, "piaa1Cmodes.fits", "piaa1Cmodescoeff", &piaacmcopticaldesign.piaa1CmodesID);
    }
    if (err == RETURN_SUCCESS)
    {
        err = load_piaa_mode_coeffs(
            dname, "piaa1Fmodes.fits", "piaa1Fmodescoeff", &piaacmcopticaldesign.piaa1FmodesID);
    }
    if (err != RETURN_SUCCESS)
    {
        fclose(fp);
        return err;
    }

    float tmpf;
    fscanfcnt = fscanf(fp, "%f   fpmaskamptransm\n", &tmpf);
    check_fscanf_count(fp, fscanfcnt, 1);
    piaacmcopticaldesign.fpmaskamptransm = tmpf;

    fclose(fp);

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * PIAACMCsimul_update_fnamedescr_conf - Format geometry and bandpass FPM tag fnamedescr_conf.
 *
 * Return: RETURN_SUCCESS.
 */
errno_t PIAACMCsimul_update_fnamedescr_conf(void)
{
    DEBUG_TRACE_FSTART();

    WRITE_FILENAME(
        piaacmcparams.fnamedescr_conf,
        "s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_wb%02d",
        piaacmcparams.PIAACMC_FPMsectors,
        (long) (1.0e9 * piaacmcopticaldesign.lambda + 0.1),
        (long) (1.0 * piaacmcopticaldesign.lambdaB + 0.1),
        piaacmcopticaldesign.NBrings,
        (long) (100.0 * piaacmcparams.PIAACMC_MASKRADLD + 0.1),
        piaacmcparams.computePSF_ResolvedTarget,
        piaacmcparams.computePSF_ResolvedTarget_mode,
        piaacmcopticaldesign.nblambda);

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * PIAACMCsimul_update_fnamedescr - Format full FPM solution filename tag fnamedescr.
 *
 * Return: RETURN_SUCCESS.
 */
errno_t PIAACMCsimul_update_fnamedescr(void)
{
    DEBUG_TRACE_FSTART();

    PIAACMCsimul_update_fnamedescr_conf();

    SNPRINTF_CHECK(
        piaacmcparams.fnamedescr,
        STRINGMAXLEN_FILENAME,
        "%s.minsag%06ld_maxsag%06ld_fpmregc%08ld_fpmrega%06ld_%s",
        piaacmcparams.fnamedescr_conf,
        (long) (1.0e9 * piaacmcopticaldesign.fpmminsag - 0.1),
        (long) (1.0e9 * piaacmcopticaldesign.fpmmaxsag + 0.1),
        (long) (1.0e9 * piaacmcopticaldesign.fpmsagreg_coeff + 0.1),
        (long) (1000.0 * piaacmcopticaldesign.fpmsagreg_alpha + 0.1),
        piaacmcopticaldesign.fpmmaterial_name);

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * save_image_if_valid - Write an in-memory image to FITS when its image ID is valid.
 * @ID:    Image ID (skipped if -1)
 * @dname: Destination directory
 * @rel:   Relative FITS filename
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t save_image_if_valid(
    long        ID,
    const char *dname,
    const char *rel)
{
    if (ID == -1)
    {
        return RETURN_SUCCESS;
    }

    char fname[STRINGMAXLEN_FULLFILENAME];
    WRITE_FULLFILENAME(fname, "%s/%s", dname, rel);
    FUNC_CHECK_RETURN(save_fits(dcimg[ID].name, fname));

    return RETURN_SUCCESS;
}

/**
 * save_fpm_maps - Save FPM zone sag/amplitude vectors and high-resolution 2D sag map.
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t save_fpm_maps(void)
{
    char fname[STRINGMAXLEN_FULLFILENAME];

    PIAACMCsimul_update_fnamedescr();

    if (piaacmcopticaldesign.zonezID != -1)
    {
        SNPRINTF_CHECK(
            fname,
            STRINGMAXLEN_FULLFILENAME,
            "%s/fpm_zonez.%s.fits",
            piaacmcparams.piaacmcconfdir,
            piaacmcparams.fnamedescr);
        FUNC_CHECK_RETURN(save_fits(dcimg[piaacmcopticaldesign.zonezID].name, fname));
    }

    if (piaacmcopticaldesign.zoneaID != -1)
    {
        SNPRINTF_CHECK(
            fname,
            STRINGMAXLEN_FULLFILENAME,
            "%s/fpm_zonea.%s.fits",
            piaacmcparams.piaacmcconfdir,
            piaacmcparams.fnamedescr);
        FUNC_CHECK_RETURN(save_fits(dcimg[piaacmcopticaldesign.zoneaID].name, fname));
    }

    imageID IDfpmzmap1 = image_ID("fpmzmap1", dcimg, dcnimg);
    if (IDfpmzmap1 == -1)
    {
        printf("Creating fpmzmap1 ...\n");
        fflush(stdout);

        FUNC_CHECK_RETURN(mkFPM_zonemap("fpmzmap1", &IDfpmzmap1));
        uint64_t xysize = (uint64_t) dcimg[IDfpmzmap1].md[0].size[0] *
                          dcimg[IDfpmzmap1].md[0].size[1];
        for (uint64_t ii = 0; ii < xysize; ii++)
        {
            dcimg[IDfpmzmap1].array.UI16[ii] -= 1;
        }
    }

    image_basic_indexmap("fpmzmap1", dcimg[piaacmcopticaldesign.zonezID].name, "fpmsagmapHR");

    SNPRINTF_CHECK(
        fname,
        STRINGMAXLEN_FULLFILENAME,
        "%s/fpm_sagmapHR.%s.fits",
        piaacmcparams.piaacmcconfdir,
        piaacmcparams.fnamedescr);
    FUNC_CHECK_RETURN(save_fits("fpmsagmapHR", fname));
    FUNC_CHECK_RETURN(delete_image_ID("fpmsagmapHR", DELETE_IMAGE_ERRMODE_WARNING));

    return RETURN_SUCCESS;
}

/**
 * PIAACMCsimul_savepiaacmcconf - Write piaacmcparams.conf, Lyot stops, PIAA modes, and FPM maps.
 * @dname: Destination configuration directory path
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
errno_t PIAACMCsimul_savepiaacmcconf(const char *restrict dname)
{
    DEBUG_TRACE_FSTART();

    EXECUTE_SYSTEM_COMMAND("mkdir -p %s", dname);

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/piaacmcparams.conf", dname);
        FILE *fp = fopen(fname, "w");
        if (fp == NULL)
        {
            FUNC_RETURN_FAILURE("Cannot write piaacmcparams.conf");
        }

        fprintf(fp, "%10ld   NBradpts\n", piaacmcopticaldesign.NBradpts);
        for (int i = 0; i < 10; i++)
        {
            if ((i < piaacmcopticaldesign.NBLyotStop) &&
                (piaacmcopticaldesign.IDLyotStop[i] != -1))
            {
                WRITE_FULLFILENAME(fname, "%s/LyotStop%d.fits", dname, i);
                save_fits(dcimg[piaacmcopticaldesign.IDLyotStop[i]].name, fname);
            }
            fprintf(fp, "%10.6lf   LyotStop_zpos %d\n", piaacmcopticaldesign.LyotStop_zpos[i], i);
        }
        fprintf(fp, "%10.6f    fpmaskamptransm\n", piaacmcopticaldesign.fpmaskamptransm);
        fclose(fp);
    }

    FUNC_CHECK_RETURN(
        save_image_if_valid(piaacmcopticaldesign.piaa0CmodesID, dname, "piaa0Cmodes.fits"));
    FUNC_CHECK_RETURN(
        save_image_if_valid(piaacmcopticaldesign.piaa0FmodesID, dname, "piaa0Fmodes.fits"));
    FUNC_CHECK_RETURN(
        save_image_if_valid(piaacmcopticaldesign.piaa1CmodesID, dname, "piaa1Cmodes.fits"));
    FUNC_CHECK_RETURN(
        save_image_if_valid(piaacmcopticaldesign.piaa1FmodesID, dname, "piaa1Fmodes.fits"));

    FUNC_CHECK_RETURN(save_fpm_maps());

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
