/**
 * @file    PIAAACMCsimul_eval_poly_design.c
 * @brief   Evaluate polychromatic design
 *
 */

// System includes
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_tools/COREMOD_tools.h"

#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"

#include "PIAACMCsimul_computePSF.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"
#include "init_piaacmcopticaldesign.h"
#include "init_piaacmcopticalsystem.h"

#include "PIAAshape/makePIAAshapes.h"

errno_t PIAACMCsimul_eval_poly_design()
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FUNC");

    double fpmradld = 0.95; // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;
    double ldoffset = 0.01; // default;

    // off-axis point source reference flux
    double valref;

    uint32_t xsize, ysize, zsize;
    uint64_t xysize;

    double eval_sepstepld = 0.2; // in l/D
    double eval_sepmaxld = 20.0; // in l/D

    long size;

    printf("=================================== mode 100 ===================================\n");

    {
        variableID IDv;

        IDv = variable_ID("PIAACMC_centobs0");
        if (IDv != -1)
        {
            centobs0 = data.variable[IDv].value.f;
        }

        IDv = variable_ID("PIAACMC_centobs1");
        if (IDv != -1)
        {
            centobs1 = data.variable[IDv].value.f;
        }

        IDv = variable_ID("PIAACMC_fpmradld");
        if (IDv != -1)
        {
            fpmradld = data.variable[IDv].value.f;
            printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
        }
    }

    // measure sensitivity to errors
    //	printf("Loading (optional) OPDerr file\n");
    //	fflush(stdout);
    // load an error if it exists
    imageID IDopderrC = image_ID("OPDerrC");
    if (IDopderrC == -1)
    {
        FUNC_CHECK_RETURN(load_fits("OPDerrC.fits", "OPDerrC", 0, &IDopderrC));
    }

    long nbOPDerr = 0;
    if (IDopderrC != -1)
    {
        nbOPDerr = data.image[IDopderrC].md[0].size[2]; // number of error arrays
        //printf("INCLUDING %ld OPD ERROR MODES\n", nbOPDerr);
        //fflush(stdout);
    }

    printf("Will add optional OPD error modes (%ld modes)\n", nbOPDerr);
    fflush(stdout);

    piaacmcparams.PIAACMC_fpmtype = 0; // idealized (default)
    {
        variableID IDv;
        IDv = variable_ID("PIAACMC_fpmtype");
        if (IDv != -1)
        {
            piaacmcparams.PIAACMC_fpmtype = (int)(data.variable[IDv].value.f + 0.1);
        }
    }

    piaacmcparams.FORCE_CREATE_fpmza = 1;

    {
        uint64_t initflag = INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__READCONF;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF;

        if (piaacmcparams.PIAACMC_fpmtype == 1)
        {
            initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL;
        }
        FUNC_CHECK_RETURN(init_piaacmcopticaldesign(fpmradld, centobs0, centobs1, initflag, NULL));
    }

    FUNC_CHECK_RETURN(makePIAAshapes());

    piaacmcopticalsystem.FOCMASKarray[0].mode = 1; // use 1-fpm

    //        valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, piaacmcopticalsystem.NBelem, 1, 0);
    //       printf("valref = %g\n", valref);

    //exit(0);

    {
        variableID IDv = variable_ID("PIAACMC_ldoffset");
        if (IDv != -1)
        {
            ldoffset = data.variable[IDv].value.f;
        }
    }

    printf("ldoffset = %f\n", ldoffset);

    // compute off-axis POINT source reference
    FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(5.0, 0.0, 0, piaacmcopticalsystem.NBelem, 1, 0, 0, 0, &valref));

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/psfi0_x50_y00.fits", piaacmcparams.piaacmcconfdir);
        FUNC_CHECK_RETURN(save_fits("psfi0", fname));
    }
    //load_fits(fname, "psfi");

    float avpeak;
    {
        imageID ID = image_ID("psfi0");
        xsize = data.image[ID].md[0].size[0];
        ysize = data.image[ID].md[0].size[1];
        xysize = xsize;
        xysize *= ysize;
        zsize = data.image[ID].md[0].size[2];

        float *peakarray = (float *)malloc(sizeof(float) * zsize);
        if (peakarray == NULL)
        {
            FUNC_RETURN_FAILURE("malloc error");
        }

        for (uint32_t kk = 0; kk < zsize; kk++)
        {
            peakarray[kk] = 0.0;
            for (uint64_t ii = 0; ii < xysize; ii++)
            {
                double val = data.image[ID].array.F[kk * xsize * ysize + ii];
                if (val > peakarray[kk])
                {
                    peakarray[kk] = val;
                }
            }
        }

        avpeak = 0.0;
        for (uint32_t kk = 0; kk < zsize; kk++)
        {
            printf("peak %02ld  %10lf\n", (long)kk, peakarray[kk]);
            avpeak += peakarray[kk];
        }
        avpeak /= zsize;
        free(peakarray);
        FUNC_CHECK_RETURN(delete_image_ID("psfi0", DELETE_IMAGE_ERRMODE_WARNING));
    }

    // compute on-axis POINT source
    FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(0.0, 0.0, 0, piaacmcopticalsystem.NBelem, 1, 0, 0, 1, &valref));

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/psfi0_x00_y00.fits", piaacmcparams.piaacmcconfdir);
        save_fits("psfi0", fname);
    }
    //load_fits(fname, "psfi");

    { /// compute contrast curve
        imageID ID = image_ID("psfi0");

        double focscale =
            (2.0 * piaacmcopticaldesign.beamrad / piaacmcopticaldesign.pixscale) / piaacmcopticaldesign.size;
        printf("focscale = %f\n", focscale);
        long eval_sepNBpt = (long)(eval_sepmaxld / eval_sepstepld);

        double *eval_contrastCurve = (double *)malloc(sizeof(double) * eval_sepNBpt);
        if (eval_contrastCurve == NULL)
        {
            FUNC_RETURN_FAILURE("malloc error");
        }

        double *eval_contrastCurve_cnt = (double *)malloc(sizeof(double) * eval_sepNBpt);
        if (eval_contrastCurve_cnt == NULL)
        {
            FUNC_RETURN_FAILURE("malloc error");
        }

        for (long ri = 0; ri < eval_sepNBpt; ri++)
        {
            eval_contrastCurve[ri] = 0.0;
            eval_contrastCurve_cnt[ri] = 0.0;
        }

        double aveC = 0.0;
        long aveCcnt = 0;

        for (uint32_t kk = 0; kk < zsize; kk++)
            for (uint32_t ii = 0; ii < xsize; ii++)
                for (uint32_t jj = 0; jj < ysize; jj++)
                {
                    double xc = 1.0 * ii - 0.5 * xsize;
                    double yc = 1.0 * jj - 0.5 * ysize;
                    xc *= focscale;
                    yc *= focscale;
                    double rc = sqrt(xc * xc + yc * yc);
                    long ri = (long)(rc / eval_sepstepld - 0.5);
                    if (ri < 0)
                    {
                        ri = 0;
                    }
                    if (ri < eval_sepNBpt)
                    {
                        eval_contrastCurve[ri] += data.image[ID].array.F[kk * xsize * ysize + jj * xsize + ii] / avpeak;
                        eval_contrastCurve_cnt[ri] += 1.0;
                    }
                    if ((rc > 2.0) && (rc < 6.0))
                    {
                        aveC += data.image[ID].array.F[kk * xsize * ysize + jj * xsize + ii] / avpeak;

                        aveCcnt++;
                    }
                }

        PIAACMCsimul_update_fnamedescr();
        {
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(fname, "%s/ContrastCurve_tt000.%s.txt", piaacmcparams.piaacmcconfdir,
                               piaacmcparams.fnamedescr);
            FILE *fp = fopen(fname, "w");
            fprintf(fp, "# Contrast curve\n");
            fprintf(fp, "# col 1: Angular separation [l/D]\n");
            fprintf(fp, "# col 2: Contrast value\n");
            fprintf(fp, "# col 3: Number of points used for contrast measurement\n");
            fprintf(fp, "\n");
            for (long ri = 0; ri < eval_sepNBpt; ri++)
            {
                eval_contrastCurve[ri] /= eval_contrastCurve_cnt[ri] + 0.000001;
                fprintf(fp, "%10f %10g %10g\n", eval_sepstepld * ri, eval_contrastCurve[ri],
                        eval_contrastCurve_cnt[ri]);
            }
            fclose(fp);
        }

        free(eval_contrastCurve);
        free(eval_contrastCurve_cnt);

        PIAACMCsimul_update_fnamedescr();
        {
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(fname, "%s/ContrastVal_tt000.%s.txt", piaacmcparams.piaacmcconfdir,
                               piaacmcparams.fnamedescr);
            FILE *fp = fopen(fname, "w");
            fprintf(fp, "%10g %10g %4.2f %4.2f %4.2f %4.2f %04ld %02ld %d %03ld %03ld %02d %03ld %02d %d\n", valref,
                    aveC / aveCcnt, piaacmcopticaldesign.fpmaskradld, piaacmcparams.PIAACMC_MASKRADLD,
                    piaacmcopticaldesign.centObs0, piaacmcopticaldesign.centObs1,
                    (long)(piaacmcopticaldesign.lambda * 1e9), (long)(piaacmcopticaldesign.lambdaB + 0.1),
                    piaacmcparams.PIAACMC_FPMsectors, piaacmcopticaldesign.NBrings, piaacmcopticaldesign.focmNBzone,
                    piaacmcopticaldesign.nblambda, (long)0, piaacmcparams.computePSF_ResolvedTarget,
                    piaacmcparams.computePSF_ResolvedTarget_mode);
            fclose(fp);
        }
    }

    // measure pointing sensitivity
    imageID IDps;
    FUNC_CHECK_RETURN(create_3Dimage_ID("starim", piaacmcopticaldesign.size, piaacmcopticaldesign.size, zsize, &IDps));

    imageID IDps_re;
    FUNC_CHECK_RETURN(
        create_3Dimage_ID("starim_re", piaacmcopticaldesign.size, piaacmcopticaldesign.size, zsize, &IDps_re));

    imageID IDps_im;
    FUNC_CHECK_RETURN(
        create_3Dimage_ID("starim_im", piaacmcopticaldesign.size, piaacmcopticaldesign.size, zsize, &IDps_im));

    imageID IDps_COH;
    FUNC_CHECK_RETURN(
        create_3Dimage_ID("starimCOH", piaacmcopticaldesign.size, piaacmcopticaldesign.size, zsize, &IDps_COH));

    imageID IDps_INC;
    FUNC_CHECK_RETURN(
        create_3Dimage_ID("starimINC", piaacmcopticaldesign.size, piaacmcopticaldesign.size, zsize, &IDps_INC));

    long NBpt = 0;

    {
        double cval = 0.0;
        FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(ldoffset, 0.0, 0, piaacmcopticalsystem.NBelem, 0, 0, 0, 0, &cval));
        valref = 0.25 * cval;
    }

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/psfi0_p0.fits", piaacmcparams.piaacmcconfdir);
        save_fits("psfi0", fname);
    }

    {
        imageID ID = image_ID("psfi0");
        imageID IDre = image_ID("psfre0");
        imageID IDim = image_ID("psfim0");
        for (uint64_t ii = 0; ii < xysize * zsize; ii++)
        {
            data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];
            data.image[IDps_re].array.F[ii] += data.image[IDre].array.F[ii];
            data.image[IDps_im].array.F[ii] += data.image[IDim].array.F[ii];
        }
    }
    NBpt++;
    delete_image_ID("psfre0", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("psfim0", DELETE_IMAGE_ERRMODE_WARNING);

    {
        double cval = 0.0;
        FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(-ldoffset, 0.0, 0, piaacmcopticalsystem.NBelem, 0, 0, 0, 0, &cval));

        valref += 0.25 * cval;
    }

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/psfi0_m0.fits", piaacmcparams.piaacmcconfdir);
        save_fits("psfi0", fname);
    }

    {
        imageID ID = image_ID("psfi0");
        imageID IDre = image_ID("psfre0");
        imageID IDim = image_ID("psfim0");
        for (uint64_t ii = 0; ii < xysize * zsize; ii++)
        {
            data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];
            data.image[IDps_re].array.F[ii] += data.image[IDre].array.F[ii];
            data.image[IDps_im].array.F[ii] += data.image[IDim].array.F[ii];
        }
    }
    NBpt++;
    delete_image_ID("psfre0", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("psfim0", DELETE_IMAGE_ERRMODE_WARNING);

    {
        double cval = 0.0;

        FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(0.0, ldoffset, 0, piaacmcopticalsystem.NBelem, 0, 0, 0, 0, &cval));
        valref += 0.25 * cval;
    }

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/psfi0_0p.fits", piaacmcparams.piaacmcconfdir);
        save_fits("psfi0", fname);
    }

    {
        imageID ID = image_ID("psfi0");
        imageID IDre = image_ID("psfre0");
        imageID IDim = image_ID("psfim0");
        for (uint64_t ii = 0; ii < xysize * zsize; ii++)
        {
            data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];
            data.image[IDps_re].array.F[ii] += data.image[IDre].array.F[ii];
            data.image[IDps_im].array.F[ii] += data.image[IDim].array.F[ii];
        }
    }
    NBpt++;
    delete_image_ID("psfre0", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("psfim0", DELETE_IMAGE_ERRMODE_WARNING);

    {
        double cval = 0.0;
        FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(0.0, -ldoffset, 0, piaacmcopticalsystem.NBelem, 0, 0, 0, 0, &cval));

        valref += 0.25 * cval;
    }

    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/psfi0_0m.fits", piaacmcparams.piaacmcconfdir);
        save_fits("psfi0", fname);
    }

    {
        imageID ID = image_ID("psfi0");
        imageID IDre = image_ID("psfre0");
        imageID IDim = image_ID("psfim0");
        for (uint64_t ii = 0; ii < xysize * zsize; ii++)
        {
            data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];
            data.image[IDps_re].array.F[ii] += data.image[IDre].array.F[ii];
            data.image[IDps_im].array.F[ii] += data.image[IDim].array.F[ii];
        }
    }
    NBpt++;
    delete_image_ID("psfre0", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("psfim0", DELETE_IMAGE_ERRMODE_WARNING);

    // measure sensitivity to errors

    //	printf("Adding optional OPD error modes (%ld modes)\n", nbOPDerr);
    //	fflush(stdout);
    // add error modes if any
    for (long OPDmode = 0; OPDmode < nbOPDerr; OPDmode++)
    {
        size = data.image[IDopderrC].md[0].size[0];

        imageID IDopderr;
        create_2Dimage_ID("opderr", size, size, &IDopderr);

        // "opderr" is a standard name read by PIAACMCsimul_init
        for (long ii = 0; ii < size * size; ii++)
        {
            data.image[IDopderr].array.F[ii] = data.image[IDopderrC].array.F[size * size * OPDmode + ii];
        }

        FUNC_CHECK_RETURN(init_piaacmcopticalsystem(0.0, 0.0)); // add error to the data

        FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(0.0, 0.0, 0, piaacmcopticalsystem.NBelem, 0, 0, 0, 0, NULL));

        {
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(fname, "%s/psfi0_opderr%02ld.fits", piaacmcparams.piaacmcconfdir, OPDmode);
            save_fits("psfi0", fname);
        }

        FUNC_CHECK_RETURN(delete_image_ID("opderr", DELETE_IMAGE_ERRMODE_WARNING));

        {
            imageID ID = image_ID("psfi0");
            imageID IDre = image_ID("psfre0");
            imageID IDim = image_ID("psfim0");
            for (uint64_t ii = 0; ii < xysize * zsize; ii++)
            {
                data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];
                data.image[IDps_re].array.F[ii] += data.image[IDre].array.F[ii];
                data.image[IDps_im].array.F[ii] += data.image[IDim].array.F[ii];
            }
        }
        NBpt++;
        delete_image_ID("psfre0", DELETE_IMAGE_ERRMODE_WARNING);
        delete_image_ID("psfim0", DELETE_IMAGE_ERRMODE_WARNING);
    }

    for (uint64_t ii = 0; ii < xysize * zsize; ii++)
    {
        data.image[IDps].array.F[ii] /= NBpt;
        data.image[IDps_re].array.F[ii] /= NBpt; // average re
        data.image[IDps_im].array.F[ii] /= NBpt; // average im
    }

    for (uint64_t ii = 0; ii < xysize * zsize; ii++)
    {
        data.image[IDps_COH].array.F[ii] = data.image[IDps_re].array.F[ii] * data.image[IDps_re].array.F[ii] +
                                           data.image[IDps_im].array.F[ii] * data.image[IDps_im].array.F[ii];

        data.image[IDps_INC].array.F[ii] = data.image[IDps].array.F[ii] - data.image[IDps_COH].array.F[ii];
    }

    // same as psfi0_extsrc
    //    sprintf(fname, "%s/psfi0_starim.%s.fits", piaacmcparams.piaacmcconfdir, piaacmcparams.fnamedescr);
    //    save_fits("starim", fname);

    PIAACMCsimul_update_fnamedescr();
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/psfi0_extsrc%2ld_sm%d.%s.fits", piaacmcparams.piaacmcconfdir,
                           (long)(-log10(ldoffset) * 10.0 + 0.1), piaacmcparams.SCORINGMASKTYPE,
                           piaacmcparams.fnamedescr);
        FUNC_CHECK_RETURN(save_fits("starim", fname));
    }

    PIAACMCsimul_update_fnamedescr();
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/psfi0INC_extsrc%2ld_sm%d.%s.fits", piaacmcparams.piaacmcconfdir,
                           (long)(-log10(ldoffset) * 10.0 + 0.1), piaacmcparams.SCORINGMASKTYPE,
                           piaacmcparams.fnamedescr);
        FUNC_CHECK_RETURN(save_fits("starimINC", fname));
    }

    PIAACMCsimul_update_fnamedescr();
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname, "%s/psfi0COH_extsrc%2ld_sm%d.%s.fits", piaacmcparams.piaacmcconfdir,
                           (long)(-log10(ldoffset) * 10.0 + 0.1), piaacmcparams.SCORINGMASKTYPE,
                           piaacmcparams.fnamedescr);
        FUNC_CHECK_RETURN(save_fits("starimCOH", fname));
    }

    { /// compute contrast curve
        /// measure average contrast value, 2-6 lambda/D
        double focscale =
            (2.0 * piaacmcopticaldesign.beamrad / piaacmcopticaldesign.pixscale) / piaacmcopticaldesign.size;
        printf("focscale = %f\n", focscale);
        long eval_sepNBpt = (long)(eval_sepmaxld / eval_sepstepld);

        double *eval_contrastCurve = (double *)malloc(sizeof(double) * eval_sepNBpt);
        if (eval_contrastCurve == NULL)
        {
            FUNC_RETURN_FAILURE("malloc error");
        }

        double *eval_COHcontrastCurve = (double *)malloc(sizeof(double) * eval_sepNBpt);
        if (eval_COHcontrastCurve == NULL)
        {
            FUNC_RETURN_FAILURE("malloc error");
        }

        double *eval_INCcontrastCurve = (double *)malloc(sizeof(double) * eval_sepNBpt);
        if (eval_INCcontrastCurve == NULL)
        {
            FUNC_RETURN_FAILURE("malloc error");
        }

        double *eval_contrastCurve_cnt = (double *)malloc(sizeof(double) * eval_sepNBpt);
        if (eval_contrastCurve_cnt == NULL)
        {
            FUNC_RETURN_FAILURE("malloc error");
        }

        for (long ri = 0; ri < eval_sepNBpt; ri++)
        {
            eval_contrastCurve[ri] = 0.0;
            eval_COHcontrastCurve[ri] = 0.0;
            eval_INCcontrastCurve[ri] = 0.0;
            eval_contrastCurve_cnt[ri] = 0.0;
        }

        double aveC = 0.0;
        double aveC_INC = 0.0;
        double aveC_COH = 0.0;
        long aveCcnt = 0;

        imageID IDrc;
        FUNC_CHECK_RETURN(create_3Dimage_ID("_tmp_rc", xsize, ysize, zsize, &IDrc));

        for (uint32_t kk = 0; kk < zsize; kk++)
        {
            for (uint32_t ii = 0; ii < xsize; ii++)
            {
                for (uint32_t jj = 0; jj < ysize; jj++)
                {
                    double xc = 1.0 * ii - 0.5 * xsize;
                    double yc = 1.0 * jj - 0.5 * ysize;
                    xc *= focscale;
                    yc *= focscale;
                    double rc = sqrt(xc * xc + yc * yc);
                    data.image[IDrc].array.F[kk * xsize * ysize + jj * xsize + ii] = rc;

                    long ri = (long)(rc / eval_sepstepld - 0.5);
                    if (ri < 0)
                    {
                        ri = 0;
                    }
                    if (ri < eval_sepNBpt)
                    {
                        eval_contrastCurve[ri] +=
                            data.image[IDps].array.F[kk * xsize * ysize + jj * xsize + ii] / avpeak;

                        eval_COHcontrastCurve[ri] +=
                            data.image[IDps_COH].array.F[kk * xsize * ysize + jj * xsize + ii] / avpeak;

                        eval_INCcontrastCurve[ri] +=
                            data.image[IDps_INC].array.F[kk * xsize * ysize + jj * xsize + ii] / avpeak;

                        eval_contrastCurve_cnt[ri] += 1.0;
                    }
                    if ((rc > 2.0) && (rc < 6.0))
                    {
                        aveC += data.image[IDps].array.F[kk * xsize * ysize + jj * xsize + ii] / avpeak;

                        aveC_COH += data.image[IDps_COH].array.F[kk * xsize * ysize + jj * xsize + ii] / avpeak;

                        aveC_INC += data.image[IDps_INC].array.F[kk * xsize * ysize + jj * xsize + ii] / avpeak;
                        aveCcnt++;
                    }
                }
            }
        }

        PIAACMCsimul_update_fnamedescr();
        {
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(fname, "%s/ContrastCurveP_extsrc%2ld_sm%d.%s.txt", piaacmcparams.piaacmcconfdir,
                               (long)(-log10(ldoffset) * 10.0 + 0.1), piaacmcparams.SCORINGMASKTYPE,
                               piaacmcparams.fnamedescr);

            FILE *fp = fopen(fname, "w");

            float *rcarray = (float *)malloc(sizeof(float) * xsize * ysize * zsize);
            if (rcarray == NULL)
            {
                FUNC_RETURN_FAILURE("malloc error");
            }

            float drc = 0.2;
            long rcbinmax = 100;
            for (long rcbin = 0; rcbin < rcbinmax; rcbin++)
            {
                float rcmin = drc * rcbin;
                float rcmax = drc * (rcbin + 1);

                long rc_cnt = 0;

                for (uint32_t kk = 0; kk < zsize; kk++)
                    for (uint32_t ii = 0; ii < xsize; ii++)
                        for (uint32_t jj = 0; jj < ysize; jj++)
                        {
                            double rc = data.image[IDrc].array.F[kk * xsize * ysize + jj * xsize + ii];
                            if ((rc > rcmin) && (rc < rcmax))
                            {
                                rcarray[rc_cnt] =
                                    data.image[IDps].array.F[kk * xsize * ysize + jj * xsize + ii] / avpeak;

                                rc_cnt++;
                            }
                        }
                quick_sort_float(rcarray, rc_cnt);

                float p50 = rcarray[(long)(0.5 * rc_cnt)];
                float p20 = rcarray[(long)(0.2 * rc_cnt)];
                fprintf(fp, "%5.2f %12g %12g\n", drc * (rcbin + 0.5), p50, p20);
            }

            free(rcarray);
            FUNC_CHECK_RETURN(delete_image_ID("_tmp_rc", DELETE_IMAGE_ERRMODE_WARNING));
            fclose(fp);
        }

        FUNC_CHECK_RETURN(delete_image_ID("starimINC", DELETE_IMAGE_ERRMODE_WARNING));
        FUNC_CHECK_RETURN(delete_image_ID("starimCOH", DELETE_IMAGE_ERRMODE_WARNING));

        PIAACMCsimul_update_fnamedescr();
        {
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(fname, "%s/ContrastCurve_extsrc%2ld_sm%d.%s.txt", piaacmcparams.piaacmcconfdir,
                               (long)(-log10(ldoffset) * 10.0 + 0.1), piaacmcparams.SCORINGMASKTYPE,
                               piaacmcparams.fnamedescr);

            FILE *fp = fopen(fname, "w");
            fprintf(fp, "# Contrast curves\n");
            fprintf(fp, "# col 1 : separation [l/D]\n");
            fprintf(fp, "# col 2 : contrast, total intensity, radially averaged\n");
            fprintf(fp, "# col 3 : contrast, Coherent component, radially averaged\n");
            fprintf(fp, "# col 4 : contrast, incoherent component, radially averaged\n");
            fprintf(fp, "\n");
            for (long ri = 0; ri < eval_sepNBpt; ri++)
            {
                eval_contrastCurve[ri] /= eval_contrastCurve_cnt[ri] + 0.000001;
                eval_COHcontrastCurve[ri] /= eval_contrastCurve_cnt[ri] + 0.000001;
                eval_INCcontrastCurve[ri] /= eval_contrastCurve_cnt[ri] + 0.000001;
                fprintf(fp, "%10f %10g %10g %10g %10g\n", eval_sepstepld * ri, eval_contrastCurve[ri],
                        eval_contrastCurve_cnt[ri], eval_COHcontrastCurve[ri], eval_INCcontrastCurve[ri]);
            }
            fclose(fp);
        }

        free(eval_contrastCurve);
        free(eval_COHcontrastCurve);
        free(eval_INCcontrastCurve);
        free(eval_contrastCurve_cnt);

        PIAACMCsimul_update_fnamedescr();
        {
            char fname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(fname, "%s/ContrastVal_extsrc%2ld_sm%d.%s.txt", piaacmcparams.piaacmcconfdir,
                               (long)(-log10(ldoffset) * 10.0 + 0.1), piaacmcparams.SCORINGMASKTYPE,
                               piaacmcparams.fnamedescr);
            FILE *fp = fopen(fname, "w");
            fprintf(fp, "# Contrast value\n");
            fprintf(fp,
                    "# valref, aveC/aveCcnt, piaacmcopticaldesign.fpmaskradld, piaacmcparams.PIAACMC_MASKRADLD, "
                    "piaacmcopticaldesign.centObs0, piaacmcopticaldesign.centObs1, (long) "
                    "(piaacmcopticaldesign.lambda*1e9), (long) (piaacmcopticaldesign.lambdaB+0.1), "
                    "piaacmcparams.PIAACMC_FPMsectors, piaacmcopticaldesign.NBrings, piaacmcopticaldesign.focmNBzone, "
                    "piaacmcopticaldesign.nblambda, (long) (1000.0*ldoffset), piaacmcopticaldesign.fpmsagreg_coeff, "
                    "piaacmcparams.computePSF_ResolvedTarget, piaacmcparams.computePSF_ResolvedTarget_mode\n");
            fprintf(fp, "#\n");
            fprintf(fp, "%10g %10g %4.2f %4.2f %4.2f %4.2f %04ld %02ld %d %03ld %03ld %02d %03ld %7.3f %02d %d\n",
                    valref, aveC / aveCcnt, piaacmcopticaldesign.fpmaskradld, piaacmcparams.PIAACMC_MASKRADLD,
                    piaacmcopticaldesign.centObs0, piaacmcopticaldesign.centObs1,
                    (long)(piaacmcopticaldesign.lambda * 1e9), (long)(piaacmcopticaldesign.lambdaB + 0.1),
                    piaacmcparams.PIAACMC_FPMsectors, piaacmcopticaldesign.NBrings, piaacmcopticaldesign.focmNBzone,
                    piaacmcopticaldesign.nblambda, (long)(1000.0 * ldoffset), piaacmcopticaldesign.fpmsagreg_coeff,
                    piaacmcparams.computePSF_ResolvedTarget, piaacmcparams.computePSF_ResolvedTarget_mode);
            fclose(fp);
        }
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
