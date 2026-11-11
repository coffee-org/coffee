/**
 * @file    PIAACMCsimul_computePSF_fast.c
 * @brief   Focal-plane scoring mask, fast FPMresp evaluator, and PSF flux/Cnorm telemetry
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "linopt_imtools/linopt_imtools.h"

#include "PIAACMCsimul.h"
#include "PIAACMCsimul_achromFPMsol_eval.h"
#include "PIAACMCsimul_computePSF_fast.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"

/**
 * PIAACMCsimul_init_scoringmask - Build focal-plane dark-hole evaluation mask and pixel list.
 * @size:     2D focal-plane grid dimension [pix]
 * @focscale: Focal-plane sampling scale [lambda/D per pix]
 *
 * Creates `"scoringmask"` (if not yet present), converts active pixels into
 * `"pixindex"` and `"pixmult"`, and stores total active pixel count in
 * `piaacmcparams.SCORINGTOTAL`.
 *
 * Return: RETURN_SUCCESS or `linopt_imtools_mask_to_pixtable` error status.
 */
errno_t PIAACMCsimul_init_scoringmask(
    uint32_t size,
    float    focscale)
{
    DEBUG_TRACE_FSTART();

    const float scoringIWA   = 1.5f;
    const float scoringOWA   = 20.0f;
    const float scoringOWAhr = 8.0f;
    const float scoringIWAx  = -20.5f;

    if (image_ID("scoringmask", dcimg, dcnimg) != -1)
    {
        DEBUG_TRACE_FEXIT();
        return RETURN_SUCCESS;
    }

    printf("CREATING SCORING MASK\n");
    printf("FOCAL PLANE SCALE = %f l/d per pix\n", focscale);
    fflush(stdout);

    imageID IDsm;
    create_2Dimage_ID("scoringmask", size, size, &IDsm);

    if (piaacmcparams.SCORINGMASKTYPE == 0)
    {
        for (uint32_t ii = 0; ii < size; ii++)
        {
            for (uint32_t jj = 0; jj < size; jj++)
            {
                double x = (1.0 * ii - 0.5 * size) * focscale;
                double y = (1.0 * jj - 0.5 * size) * focscale;
                double r = sqrt(x * x + y * y);

                if ((r > scoringIWA) && (r < scoringOWAhr) && (x > scoringIWAx) &&
                    ((ii + jj) % 2 == 0))
                {
                    dcimg[IDsm].array.F[jj * size + ii] = 1.0;
                }
                if ((r > scoringOWAhr) && (r < scoringOWA) && (x > scoringIWAx) &&
                    (ii % 2 == 0) && (jj % 2 == 0))
                {
                    dcimg[IDsm].array.F[jj * size + ii] = 1.0;
                }
                if ((x > scoringOWA) && (fabs(y) < scoringIWA * 0.25) && (r < 70.0) &&
                    ((ii + jj) % 2 == 0))
                {
                    dcimg[IDsm].array.F[jj * size + ii] = 1.0;
                }
            }
        }
    }
    else
    {
        for (uint32_t ii = 0; ii < size; ii++)
        {
            for (uint32_t jj = 0; jj < size; jj++)
            {
                double x = (1.0 * ii - 0.5 * size) * focscale;
                double y = (1.0 * jj - 0.5 * size) * focscale;
                double r = sqrt(x * x + y * y);

                if ((r > scoringIWA) && (r < scoringOWAhr) && (x > scoringIWAx) &&
                    (ii % 2 == 0) && (jj % 2 == 0))
                {
                    dcimg[IDsm].array.F[jj * size + ii] = 1.0;
                }
            }
        }
    }

    if (piaacmcparams.PIAACMC_save == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        WRITE_FULLFILENAME(fname,
                           "%s/scoringmask%d.fits",
                           piaacmcparams.piaacmcconfdir,
                           piaacmcparams.SCORINGMASKTYPE);
        save_fits("scoringmask", fname);
    }

    FUNC_CHECK_RETURN(
        linopt_imtools_mask_to_pixtable("scoringmask", "pixindex", "pixmult", NULL));
    piaacmcparams.SCORINGTOTAL = arith_image_total("scoringmask");

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * PIAACMCsimul_eval_fast_fpmresp - Evaluate polychromatic contrast via pre-computed `FPMresp`.
 * @size:           2D wavefront grid dimension [pix]
 * @focscale:       Focal-plane sampling scale [lambda/D per pix]
 * @avContrast_out: Output average normalized contrast across scoring pixels
 *
 * Return: RETURN_SUCCESS on completion.
 */
errno_t PIAACMCsimul_eval_fast_fpmresp(
    uint32_t size,
    float    focscale,
    double  *avContrast_out)
{
    DEBUG_TRACE_FSTART();

    PIAACMCsimul_achromFPMsol_eval(piaacmcparams.fpmresp_array,
                                   piaacmcparams.zonez_array,
                                   piaacmcparams.dphadz_array,
                                   piaacmcparams.outtmp_array,
                                   piaacmcparams.vsize,
                                   dcimg[piaacmcopticaldesign.zonezID].md[0].size[0],
                                   piaacmcopticalsystem.nblambda,
                                   NULL);

    uint32_t ntot = (uint32_t) (piaacmcparams.vsize * piaacmcopticalsystem.nblambda);
    imageID  ID   = image_ID("imvect", dcimg, dcnimg);
    if (ID == -1)
    {
        create_2Dimage_ID("imvect", ntot, 1, &ID);
    }

    double value = 0.0;
    for (uint32_t ii = 0; ii < ntot; ii++)
    {
        double fval           = piaacmcparams.outtmp_array[ii];
        dcimg[ID].array.F[ii] = fval;
        value += fval * fval;
    }

    piaacmcparams.PIAACMCSIMUL_VAL0 = value;
    value                           = value / size / size / piaacmcopticalsystem.flux[0];
    *avContrast_out                 = value / (piaacmcparams.SCORINGTOTAL * focscale * focscale);

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * PIAACMCsimul_record_psf_telemetry - Log element fluxes and write flux / CnormFactor files.
 * @size:       2D wavefront grid dimension [pix]
 * @focscale:   Focal-plane sampling scale [lambda/D per pix]
 * @sourcesize: Extended stellar size code (`0` for point source)
 * @outsave:    Non-zero to write `flux_*.txt` and copy to `tmp_flux.txt`
 * @is_ext:     `1` for resolved extended-source format, `0` for point-source format
 */
void PIAACMCsimul_record_psf_telemetry(
    uint32_t size,
    float    focscale,
    int      sourcesize,
    int      outsave,
    int      is_ext)
{
    for (long elem = 0; elem < piaacmcopticalsystem.NBelem; elem++)
    {
        printf("    FLUX %3ld   %12.4lf %8.6lf\n",
               elem,
               piaacmcopticalsystem.flux[elem],
               piaacmcopticalsystem.flux[elem] / piaacmcopticalsystem.flux[0]);
    }

    if (outsave == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        PIAACMCsimul_update_fnamedescr();
        if (is_ext)
        {
            SNPRINTF_CHECK(fname,
                           STRINGMAXLEN_FULLFILENAME,
                           "%s/flux_exsrc%03d_sm%d.%s.txt",
                           piaacmcparams.piaacmcconfdir,
                           sourcesize,
                           piaacmcparams.SCORINGMASKTYPE,
                           piaacmcparams.fnamedescr);
        }
        else
        {
            SNPRINTF_CHECK(fname,
                           STRINGMAXLEN_FULLFILENAME,
                           "%s/flux_ptsrc_sm%d.%s.txt",
                           piaacmcparams.piaacmcconfdir,
                           piaacmcparams.SCORINGMASKTYPE,
                           piaacmcparams.fnamedescr);
        }

        FILE *fpflux = fopen(fname, "w");
        if (fpflux == NULL && !is_ext)
        {
            PRINT_ERROR("cannot create file %s", fname);
            abort();
        }
        if (fpflux != NULL)
        {
            for (long elem = 0; elem < piaacmcopticalsystem.NBelem; elem++)
            {
                fprintf(fpflux,
                        "%18.16lf %18.16lf  %d\n",
                        piaacmcopticalsystem.flux[elem],
                        piaacmcopticalsystem.flux[elem] / piaacmcopticalsystem.flux[0],
                        piaacmcopticalsystem.nblambda);
            }
            if (is_ext)
            {
                fprintf(fpflux, "W0  %d\n", piaacmcopticalsystem.nblambda);
            }
            else
            {
                fprintf(fpflux, "W1\n");
            }
            fclose(fpflux);
            EXECUTE_SYSTEM_COMMAND("cp %s %s/tmp_flux.txt", fname, piaacmcparams.piaacmcconfdir);
        }
    }

    piaacmcparams.CnormFactor =
        focscale * focscale * size * size * piaacmcopticalsystem.flux[0] /
        piaacmcopticalsystem.nblambda;

    if (!is_ext || piaacmcparams.WRITE_OK == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        int  tag = is_ext ? 0 : 1;
        WRITE_FULLFILENAME(fname, "%s/CnormFactor.txt", piaacmcparams.piaacmcconfdir);
        FILE *fp = fopen(fname, "w");
        if (fp == NULL && !is_ext)
        {
            PRINT_ERROR("cannot create file %s", fname);
            abort();
        }
        if (fp != NULL)
        {
            fprintf(fp, "# Written by function PIAACMCsimul_computePSF()\n");
            fprintf(fp,
                    "# Cnormfactor = focscale*focscale*size*size*"
                    "piaacmcopticalsystem.flux[0]/piaacmcopticalsystem.nblambda\n");
            fprintf(fp, "# %d    focscale  size  flux[0]  nblambda\n\n", tag);
            fprintf(fp, "%g\n", piaacmcparams.CnormFactor);
            fprintf(fp,
                    "%d      %g %ul %g %d\n",
                    tag,
                    focscale,
                    size,
                    piaacmcopticalsystem.flux[0],
                    piaacmcopticalsystem.nblambda);
            fclose(fp);
        }
    }
}
