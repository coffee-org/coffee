/**
 * @file    PIAAACMCsimul_loadsavepiaacmcconf.c
 * @brief   PIAA-type coronagraph design
 *
 * Can design both APLCMC and PIAACMC coronagraphs
 *
 *
 */

// System includes
#include <stdio.h>

// milk includes
#include "CommandLineInterface/CLIcore.h"
//   core modules
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "image_basic/image_basic.h"
//   other modules

#include "PIAACMCsimul.h"

#include "FocalPlaneMask/mkFPM_zonemap.h"

errno_t PIAACMCsimul_loadpiaacmcconf(const char *dname)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG dname %s", dname);

    FILE *fp;
    char  fname[STRINGMAXLEN_FULLFILENAME];

    WRITE_FULLFILENAME(fname, "%s/piaacmcparams.conf", dname);
    printf("%s\n", fname);
    fp = fopen(fname, "r");
    if (fp == NULL)
    {
        printf(
            "Configuration file \"%s\" does not exist (yet), using previously "
            "set configuration\n",
            fname);
        fflush(stdout);
    }
    else
    {

        { // Get number of rad points
            int  fscanfcnt;
            long tmpl;

            fscanfcnt = fscanf(fp, "%ld   NBradpts\n", &tmpl);
            if (fscanfcnt == EOF)
            {
                if (ferror(fp))
                {
                    perror("fscanf");
                }
                else
                {
                    fprintf(stderr,
                            "Error: fscanf reached end of file, no matching "
                            "characters, no matching failure\n");
                }
                abort();
            }
            else if (fscanfcnt != 1)
            {
                fprintf(stderr,
                        "Error: fscanf successfully matched and assigned %i "
                        "input items, 1 expected\n",
                        fscanfcnt);
                abort();
            }
            piaacmcopticaldesign.NBradpts = tmpl;
        }

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
                    load_fits(fname,
                              imname,
                              LOADFITS_ERRMODE_WARNING,
                              &(piaacmcopticaldesign.IDLyotStop[i])));

                {
                    int    fscanfcnt;
                    double tmplf;
                    long   tmpl;

                    fscanfcnt =
                        fscanf(fp, "%lf   LyotStop_zpos %ld\n", &tmplf, &tmpl);
                    if (fscanfcnt == EOF)
                    {
                        if (ferror(fp))
                        {
                            perror("fscanf");
                        }
                        else
                        {
                            fprintf(
                                stderr,
                                "Error: fscanf reached end of file, no "
                                "matching characters, no matching failure\n");
                        }
                        abort();
                    }
                    else if (fscanfcnt != 2)
                    {
                        fprintf(stderr,
                                "Error: fscanf successfully matched and "
                                "assigned %i input items, 2 expected\n",
                                fscanfcnt);
                        abort();
                    }
                    piaacmcopticaldesign.LyotStop_zpos[i] = tmplf;
                }
            }
            else
            {
                int    fscanfcnt;
                double tmplf;
                long   tmpl;

                fscanfcnt =
                    fscanf(fp, "%lf   LyotStop_zpos %ld\n", &tmplf, &tmpl);
                if (fscanfcnt == EOF)
                {
                    if (ferror(fp))
                    {
                        perror("fscanf");
                    }
                    else
                    {
                        fprintf(stderr,
                                "Error: fscanf reached end of file, no "
                                "matching characters, no matching failure\n");
                    }
                    abort();
                }
                else if (fscanfcnt != 2)
                {
                    fprintf(stderr,
                            "Error: fscanf successfully matched and assigned "
                            "%i input items, 2 expected\n",
                            fscanfcnt);
                    abort();
                }

                piaacmcopticaldesign.LyotStop_zpos[i] = tmplf;
            }
            printf("LYOT STOP %d POS : %lf\n",
                   i,
                   piaacmcopticaldesign.LyotStop_zpos[i]);
        }

        WRITE_FULLFILENAME(fname, "%s/piaa0Cmodes.fits", dname);
        FUNC_CHECK_RETURN(load_fits(fname,
                                    "piaa0Cmodescoeff",
                                    1,
                                    &(piaacmcopticaldesign.piaa0CmodesID)));
        if (piaacmcopticaldesign.piaa0CmodesID == -1)
        {
            WRITE_FULLFILENAME(fname, "%s/piaaref/piaa0Cmodes.fits", dname);
            FUNC_CHECK_RETURN(load_fits(fname,
                                        "piaa0Cmodescoeff",
                                        1,
                                        &(piaacmcopticaldesign.piaa0CmodesID)));
        }

        WRITE_FULLFILENAME(fname, "%s/piaa0Fmodes.fits", dname);
        FUNC_CHECK_RETURN(load_fits(fname,
                                    "piaa0Fmodescoeff",
                                    1,
                                    &(piaacmcopticaldesign.piaa0FmodesID)));
        if (piaacmcopticaldesign.piaa0FmodesID == -1)
        {
            WRITE_FULLFILENAME(fname, "%s/piaaref/piaa0Fmodes.fits", dname);
            FUNC_CHECK_RETURN(load_fits(fname,
                                        "piaa0Fmodescoeff",
                                        1,
                                        &(piaacmcopticaldesign.piaa0FmodesID)));
        }

        WRITE_FULLFILENAME(fname, "%s/piaa1Cmodes.fits", dname);
        FUNC_CHECK_RETURN(load_fits(fname,
                                    "piaa1Cmodescoeff",
                                    1,
                                    &(piaacmcopticaldesign.piaa1CmodesID)));
        if (piaacmcopticaldesign.piaa1CmodesID == -1)
        {
            WRITE_FULLFILENAME(fname, "%s/piaaref/piaa1Cmodes.fits", dname);
            FUNC_CHECK_RETURN(load_fits(fname,
                                        "piaa1Cmodescoeff",
                                        1,
                                        &(piaacmcopticaldesign.piaa1CmodesID)));
        }

        WRITE_FULLFILENAME(fname, "%s/piaa1Fmodes.fits", dname);
        FUNC_CHECK_RETURN(load_fits(fname,
                                    "piaa1Fmodescoeff",
                                    1,
                                    &(piaacmcopticaldesign.piaa1FmodesID)));
        if (piaacmcopticaldesign.piaa1FmodesID == -1)
        {
            imageID IDtmp = -1;
            WRITE_FULLFILENAME(fname, "%s/piaaref/piaa1Fmodes.fits", dname);
            FUNC_CHECK_RETURN(load_fits(fname, "piaa1Fmodescoeff", 1, &IDtmp));
            piaacmcopticaldesign.piaa1FmodesID = IDtmp;
        }

        { // read focal plane mask transmission
            int   fscanfcnt;
            float tmpf;

            fscanfcnt = fscanf(fp, "%f   fpmaskamptransm\n", &tmpf);
            if (fscanfcnt == EOF)
            {
                if (ferror(fp))
                {
                    perror("fscanf");
                }
                else
                {
                    fprintf(stderr,
                            "Error: fscanf reached end of file, no matching "
                            "characters, no matching failure\n");
                }
                abort();
            }
            else if (fscanfcnt != 1)
            {
                fprintf(stderr,
                        "Error: fscanf successfully matched and assigned %i "
                        "input items, 1 expected\n",
                        fscanfcnt);
                abort();
            }

            piaacmcopticaldesign.fpmaskamptransm = tmpf;
        }
        //r = 1;

        fclose(fp);
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/** @brief Assemble focal plane mask configuration name
 *
 */
errno_t PIAACMCsimul_update_fnamedescr_conf()
{
    DEBUG_TRACE_FSTART();

    //WRITE_FILENAME()
    WRITE_FILENAME(piaacmcparams.fnamedescr_conf,
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

/** @brief Assemble focal plane mask name
 *
 *
 */
errno_t PIAACMCsimul_update_fnamedescr()
{
    DEBUG_TRACE_FSTART();

    PIAACMCsimul_update_fnamedescr_conf();

    WRITE_FILENAME(piaacmcparams.fnamedescr,
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

//
// Save PIAACMC optical design
//

errno_t PIAACMCsimul_savepiaacmcconf(const char *__restrict dname)
{
    DEBUG_TRACE_FSTART();

    EXECUTE_SYSTEM_COMMAND("mkdir -p %s", dname);

    // piaacmcparam.conf
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/piaacmcparams.conf", dname);
        FILE *fp = fopen(fname, "w");

        fprintf(fp, "%10ld   NBradpts\n", piaacmcopticaldesign.NBradpts);
        for (int i = 0; i < 10; i++)
        {
            if (i < piaacmcopticaldesign.NBLyotStop)
            {
                WRITE_FULLFILENAME(fname, "%s/LyotStop%d.fits", dname, i);
                if (piaacmcopticaldesign.IDLyotStop[i] != -1)
                {
                    save_fits(
                        data.image[piaacmcopticaldesign.IDLyotStop[i]].name,
                        fname);
                }
                fprintf(fp,
                        "%10.6lf   LyotStop_zpos %d\n",
                        piaacmcopticaldesign.LyotStop_zpos[i],
                        i);
            }
            else
            {
                fprintf(fp,
                        "%10.6lf   LyotStop_zpos %d\n",
                        piaacmcopticaldesign.LyotStop_zpos[i],
                        i);
            }
        }
        fprintf(fp,
                "%10.6f    fpmaskamptransm\n",
                piaacmcopticaldesign.fpmaskamptransm);
        fclose(fp);
    }

    // PIAACMC optics
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/piaa0Cmodes.fits", dname);
        if (piaacmcopticaldesign.piaa0CmodesID != -1)
        {
            FUNC_CHECK_RETURN(
                save_fits(data.image[piaacmcopticaldesign.piaa0CmodesID].name,
                          fname));
        }
    }

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/piaa0Fmodes.fits", dname);
        if (piaacmcopticaldesign.piaa0FmodesID != -1)
        {
            FUNC_CHECK_RETURN(
                save_fits(data.image[piaacmcopticaldesign.piaa0FmodesID].name,
                          fname));
        }
    }

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/piaa1Cmodes.fits", dname);
        if (piaacmcopticaldesign.piaa1CmodesID != -1)
        {
            FUNC_CHECK_RETURN(
                save_fits(data.image[piaacmcopticaldesign.piaa1CmodesID].name,
                          fname));
        }
    }

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/piaa1Fmodes.fits", dname);
        if (piaacmcopticaldesign.piaa1FmodesID != -1)
        {
            FUNC_CHECK_RETURN(
                save_fits(data.image[piaacmcopticaldesign.piaa1FmodesID].name,
                          fname));
        }
    }

    // Focal plane mask
    PIAACMCsimul_update_fnamedescr();
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname,
                           "%s/fpm_zonez.%s.fits",
                           piaacmcparams.piaacmcconfdir,
                           piaacmcparams.fnamedescr);
        if (piaacmcopticaldesign.zonezID != -1)
        {
            FUNC_CHECK_RETURN(
                save_fits(data.image[piaacmcopticaldesign.zonezID].name,
                          fname));
        }
    }

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname,
                           "%s/fpm_zonea.%s.fits",
                           piaacmcparams.piaacmcconfdir,
                           piaacmcparams.fnamedescr);
        if (piaacmcopticaldesign.zoneaID != -1)
        {
            FUNC_CHECK_RETURN(
                save_fits(data.image[piaacmcopticaldesign.zoneaID].name,
                          fname));
        }
    }

    imageID IDfpmzmap1;
    IDfpmzmap1 = image_ID("fpmzmap1");
    if (IDfpmzmap1 == -1)
    {
        printf("Creating fpmzmap1 ...\n");
        fflush(stdout);

        FUNC_CHECK_RETURN(mkFPM_zonemap("fpmzmap1", &IDfpmzmap1););
        uint32_t xsize  = data.image[IDfpmzmap1].md[0].size[0];
        uint32_t ysize  = data.image[IDfpmzmap1].md[0].size[1];
        uint64_t xysize = xsize;
        xysize *= ysize;

        for (uint64_t ii = 0; ii < xysize; ii++)
        {
            data.image[IDfpmzmap1].array.UI16[ii] -= 1;
        }
    }
    //list_image_ID();
    //printf("data.image[piaacmcopticaldesign.zonezID].name = %s\n", data.image[piaacmcopticaldesign.zonezID].name);

    image_basic_indexmap("fpmzmap1",
                         data.image[piaacmcopticaldesign.zonezID].name,
                         "fpmsagmapHR");

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname,
                           "%s/fpm_sagmapHR.%s.fits",
                           piaacmcparams.piaacmcconfdir,
                           piaacmcparams.fnamedescr);

        FUNC_CHECK_RETURN(save_fits("fpmsagmapHR", fname));
    }

    FUNC_CHECK_RETURN(
        delete_image_ID("fpmsagmapHR", DELETE_IMAGE_ERRMODE_WARNING));

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
