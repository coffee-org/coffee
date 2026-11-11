/**
 * @file    PIAACMCsimul_computePSF.c
 * @brief   Full Fresnel polychromatic point-source and extended-star PSF evaluator
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "OptSystProp/OptSystProp.h"
#include "linopt_imtools/linopt_imtools.h"

#include "PIAACMCsimul.h"
#include "PIAACMCsimul_computePSF.h"
#include "PIAACMCsimul_computePSF_fast.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"
#include "PIAAshape/makePIAAshapes.h"
#include "init_piaacmcopticalsystem.h"

/**
 * propagate_quadrature_point - Propagate a single off-axis quadrature point into `imvectpXX`.
 * @x_pt:      Source X offset [lambda/D]
 * @y_pt:      Source Y offset [lambda/D]
 * @startelem: Starting optical element index
 * @imindex:   Quadrature point index (`0` initializes `psfi0ext`, `>0` accumulates)
 *
 * Return: RETURN_SUCCESS or propagation error status.
 */
static errno_t propagate_quadrature_point(
    float x_pt,
    float y_pt,
    long  startelem,
    long  imindex)
{
    DEBUG_TRACE_FSTART();

    char imname[STRINGMAXLEN_IMGNAME];

    FUNC_CHECK_RETURN(init_piaacmcopticalsystem(x_pt, y_pt));
    FUNC_CHECK_RETURN(makePIAAshapes());
    FUNC_CHECK_RETURN(OptSystProp_run(&piaacmcopticalsystem,
                                      0,
                                      startelem,
                                      piaacmcopticalsystem.NBelem,
                                      piaacmcparams.piaacmcconfdir,
                                      0));

    WRITE_IMAGENAME(imname, "imvectp%02ld", imindex);
    FUNC_CHECK_RETURN(
        linopt_imtools_image_to_vec("psfc0", "pixindex", "pixmult", imname, NULL));

    if (imindex == 0)
    {
        copy_image_ID("psfi0", "psfi0ext", 0);
    }
    else
    {
        arith_image_add_inplace("psfi0ext", "psfi0");
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * propagate_extended_quadrature - Evaluate stellar disk points and optional OPD error slices.
 * @xld:        Nominal source X offset [lambda/D]
 * @yld:        Nominal source Y offset [lambda/D]
 * @startelem:  Starting optical element index
 * @sourcesize: Stellar diameter code (`dld = 10^(-0.1 * sourcesize)`)
 * @extmode:    Quadrature geometry (`0`: 1-pt, `1`: 3-pt ring, `2`: 6-pt dual ring)
 * @IDopderrC:  Optional OPD error cube image ID (`-1` if absent)
 * @nbOPDerr:   Number of OPD error slices to append
 * @nb_pts_out: Output total quadrature + OPD mode evaluations accumulated
 *
 * Return: RETURN_SUCCESS or propagation error status.
 */
static errno_t propagate_extended_quadrature(
    float    xld,
    float    yld,
    long     startelem,
    int      sourcesize,
    int      extmode,
    imageID  IDopderrC,
    long     nbOPDerr,
    long    *nb_pts_out)
{
    DEBUG_TRACE_FSTART();

    uint32_t size    = piaacmcopticaldesign.size;
    double   dld     = 1.0 / pow(10.0, 0.1 * sourcesize);
    double   rad1    = (extmode == 2) ? (dld / sqrt(2.5)) : dld;
    double   rad2    = (extmode == 2) ? (2.0 * dld / sqrt(2.5)) : dld;
    long     imindex = 0;

    float x0 = (extmode == 0) ? xld : (float) (xld + rad1);
    FUNC_CHECK_RETURN(propagate_quadrature_point(x0, yld, startelem, imindex++));

    if (extmode > 0)
    {
        for (int k = 1; k <= 2; k++)
        {
            double pha = (2.0 * M_PI / 3.0) * k;
            FUNC_CHECK_RETURN(propagate_quadrature_point((float) (xld + rad1 * cos(pha)),
                                                         (float) (yld + rad1 * sin(pha)),
                                                         startelem,
                                                         imindex++));
        }
    }

    if (extmode == 2)
    {
        for (int k = 0; k < 3; k++)
        {
            double pha = M_PI / 3.0 + (2.0 * M_PI / 3.0) * k;
            FUNC_CHECK_RETURN(propagate_quadrature_point((float) (xld + rad2 * cos(pha)),
                                                         (float) (yld + rad2 * sin(pha)),
                                                         startelem,
                                                         imindex++));
        }
    }

    printf("Adding optional OPD error modes (%ld modes)\n", nbOPDerr);
    fflush(stdout);
    for (long OPDmode = 0; OPDmode < nbOPDerr; OPDmode++)
    {
        imageID IDopderr;
        FUNC_CHECK_RETURN(create_2Dimage_ID("opderr", size, size, &IDopderr));
        for (uint64_t ii = 0; ii < (uint64_t) size * size; ii++)
        {
            dcimg[IDopderr].array.F[ii] =
                dcimg[IDopderrC].array.F[(uint64_t) size * size * OPDmode + ii];
        }
        FUNC_CHECK_RETURN(propagate_quadrature_point(0.0f, 0.0f, startelem, imindex++));
        delete_image_ID("opderr", DELETE_IMAGE_ERRMODE_WARNING);
    }

    arith_image_cstmult_inplace("psfi0ext", 1.0 / imindex);
    *nb_pts_out = imindex;

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * assemble_extended_imvect - Pack per-point `imvectpXX` vectors into unified `imvect`.
 * @NBimindex: Number of quadrature/OPD sub-vectors to concatenate
 */
static void assemble_extended_imvect(long NBimindex)
{
    imageID ID0    = image_ID("imvectp00", dcimg, dcnimg);
    long    nbelem = dcimg[ID0].md[0].nelement;

    if (image_ID("imvect", dcimg, dcnimg) != -1)
    {
        delete_image_ID("imvect", DELETE_IMAGE_ERRMODE_WARNING);
    }

    long   offset    = nbelem / piaacmcopticaldesign.nblambda;
    long   offset1   = NBimindex * offset;
    double normcoeff = 1.0 / sqrt(NBimindex);

    printf("offset = %ld\n", offset);

    imageID ID;
    create_2Dimage_ID("imvect", offset1, piaacmcopticaldesign.nblambda, &ID);

    for (long imindex = 0; imindex < NBimindex; imindex++)
    {
        char imname[STRINGMAXLEN_IMGNAME];
        WRITE_IMAGENAME(imname, "imvectp%02ld", imindex);
        imageID ID1 = image_ID(imname, dcimg, dcnimg);
        for (long kl = 0; kl < piaacmcopticaldesign.nblambda; kl++)
        {
            for (long ii = 0; ii < offset; ii++)
            {
                dcimg[ID].array.F[kl * offset1 + imindex * offset + ii] =
                    dcimg[ID1].array.F[kl * offset + ii] * normcoeff;
            }
        }
        delete_image_ID(imname, DELETE_IMAGE_ERRMODE_ERROR);
    }
}

/**
 * finalize_extended_contrast - Compute contrast and write optional flux/contrast reports.
 * @size:           2D wavefront grid dimension [pix]
 * @focscale:       Focal-plane sampling scale [lambda/D per pix]
 * @sourcesize:     Extended stellar size code
 * @outsave:        Non-zero to write `flux_exsrc*` and `contrast_exsrc*` files
 * @avContrast_out: Output average normalized dark-hole contrast
 */
static void finalize_extended_contrast(
    uint32_t size,
    float    focscale,
    int      sourcesize,
    int      outsave,
    double  *avContrast_out)
{
    double  value        = 0.0;
    double  peakcontrast = 0.0;
    imageID ID           = image_ID("imvect", dcimg, dcnimg);

    for (uint64_t ii = 0; ii < dcimg[ID].md[0].nelement; ii++)
    {
        double tmpv = dcimg[ID].array.F[ii] * dcimg[ID].array.F[ii];
        value += tmpv;
        if (tmpv > peakcontrast)
        {
            peakcontrast = tmpv;
        }
    }

    PIAACMCsimul_record_psf_telemetry(size, focscale, sourcesize, outsave, 1);

    value             = value / size / size / piaacmcopticalsystem.flux[0];
    double sm_tot     = arith_image_total("scoringmask");
    double avContrast = value / (sm_tot * focscale * focscale);

    printf("COMPUTING RESOLVED SOURCE PSF\n");
    printf("SCORINGTOTAL = %f  %f\n", piaacmcparams.SCORINGTOTAL, sm_tot);
    printf("Peak constrast (rough estimate)= %g\n",
           peakcontrast / size / size / piaacmcopticalsystem.flux[0] / focscale / focscale);
    printf("Total light in scoring field = %g, peak PSF = %g   -> Average contrast = %g\n",
           value,
           piaacmcopticaldesign.peakPSF,
           avContrast);

    if (outsave == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        PIAACMCsimul_update_fnamedescr();
        SNPRINTF_CHECK(fname,
                       STRINGMAXLEN_FULLFILENAME,
                       "%s/contrast_exsrc%03d_sm%d.%s.txt",
                       piaacmcparams.piaacmcconfdir,
                       sourcesize,
                       piaacmcparams.SCORINGMASKTYPE,
                       piaacmcparams.fnamedescr);
        FILE *fp = fopen(fname, "w");
        fprintf(fp, "%g", avContrast);
        fclose(fp);
    }

    *avContrast_out = avContrast;
}

/**
 * compute_psf_point_source - Propagate unresolved point source and evaluate contrast.
 * @xld:            Source X offset [lambda/D]
 * @yld:            Source Y offset [lambda/D]
 * @startelem:      Starting optical element index
 * @size:           2D wavefront grid dimension [pix]
 * @focscale:       Focal-plane sampling scale [lambda/D per pix]
 * @outsave:        Non-zero to save FITS PSF, flux, and contrast outputs
 * @avContrast_out: Output average normalized dark-hole contrast
 *
 * Return: RETURN_SUCCESS or propagation/FITS error status.
 */
static errno_t compute_psf_point_source(
    float    xld,
    float    yld,
    long     startelem,
    uint32_t size,
    float    focscale,
    int      outsave,
    double  *avContrast_out)
{
    DEBUG_TRACE_FSTART();

    printf("COMPUTING UNRESOLVED SOURCE PSF [%f x %f]\n", xld, yld);

    FUNC_CHECK_RETURN(init_piaacmcopticalsystem(xld, yld));
    FUNC_CHECK_RETURN(makePIAAshapes());
    FUNC_CHECK_RETURN(OptSystProp_run(&piaacmcopticalsystem,
                                      0,
                                      startelem,
                                      piaacmcopticalsystem.NBelem,
                                      piaacmcparams.piaacmcconfdir,
                                      0));

    if (outsave == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        PIAACMCsimul_update_fnamedescr();
        SNPRINTF_CHECK(fname,
                       STRINGMAXLEN_FULLFILENAME,
                       "%s/psfi0_ptsrc_sm%d.%s.fits",
                       piaacmcparams.piaacmcconfdir,
                       piaacmcparams.SCORINGMASKTYPE,
                       piaacmcparams.fnamedescr);
        FUNC_CHECK_RETURN(save_fits("psfi0", fname));
    }

    FUNC_CHECK_RETURN(linopt_imtools_image_to_vec("psfc0", "pixindex", "pixmult", "imvect", NULL));
    FUNC_CHECK_RETURN(save_fits("imvect", "./testdir/test_imvect.fits"));

    double  value        = 0.0;
    double  peakcontrast = 0.0;
    imageID ID           = image_ID("imvect", dcimg, dcnimg);
    for (uint64_t ii = 0; ii < dcimg[ID].md[0].nelement; ii += 2)
    {
        double re   = dcimg[ID].array.F[ii];
        double im   = dcimg[ID].array.F[ii + 1];
        double tmpv = re * re + im * im;
        value += tmpv;
        if (tmpv > peakcontrast)
        {
            peakcontrast = tmpv;
        }
    }

    PIAACMCsimul_record_psf_telemetry(size, focscale, 0, outsave, 0);

    printf("COMPUTING UNRESOLVED SOURCE PSF -*- [%f x %f]\n", xld, yld);
    printf("SCORINGTOTAL = %f  %f\n", piaacmcparams.SCORINGTOTAL, arith_image_total("scoringmask"));

    if (variable_ID("PIAACMC_NOFPM") != -1)
    {
        imageID IDc                  = image_ID("psfc0", dcimg, dcnimg);
        piaacmcopticaldesign.peakPSF = 0.0;
        for (uint64_t ii = 0; ii < (uint64_t) size * size; ii++)
        {
            double re  = dcimg[IDc].array.CF[ii].re;
            double im  = dcimg[IDc].array.CF[ii].im;
            double val = re * re + im * im;
            if (val > piaacmcopticaldesign.peakPSF)
            {
                piaacmcopticaldesign.peakPSF = val;
            }
        }
        FILE *fp = fopen("conf/conf_peakPSF.txt", "w");
        fprintf(fp, "%g\n", piaacmcopticaldesign.peakPSF);
        fclose(fp);
    }

    double avContrast = 0.0;
    if (piaacmcopticaldesign.peakPSF < 1.0)
    {
        double f0sq = piaacmcopticalsystem.flux[0] * piaacmcopticalsystem.flux[0];
        avContrast  = value / f0sq / piaacmcparams.SCORINGTOTAL;
        printf("Peak constrast (rough estimate)= %g -> %g\n", peakcontrast, peakcontrast / f0sq);
        printf("piaacmcopticalsystem.flux[0]  = %g\n", piaacmcopticalsystem.flux[0]);
        printf("SCORINGMASKTYPE = %d\n", piaacmcparams.SCORINGMASKTYPE);
        printf("[0] Total light in scoring field = %g, peak PSF = %g, "
               "SCOTINGTOTAL = %g   -> Average contrast = %g\n",
               value,
               piaacmcopticaldesign.peakPSF,
               piaacmcparams.SCORINGTOTAL,
               avContrast);
    }
    else
    {
        avContrast = value / piaacmcopticaldesign.peakPSF / piaacmcparams.SCORINGTOTAL;
        printf("Peak constrast = %g -> %g\n",
               peakcontrast,
               peakcontrast / piaacmcopticaldesign.peakPSF);
        printf("SCORINGMASKTYPE = %d\n", piaacmcparams.SCORINGMASKTYPE);
        printf("[1] Total light in scoring field = %g, peak PSF = %g, "
               "SCOTINGTOTAL = %g  -> Average contrast = %g\n",
               value,
               piaacmcopticaldesign.peakPSF,
               piaacmcparams.SCORINGTOTAL,
               avContrast);

        FILE *fp = fopen("PSFcontrastval.txt", "w");
        if (fp != NULL)
        {
            fprintf(fp, "%g %g\n", peakcontrast / piaacmcopticaldesign.peakPSF, avContrast);
            fclose(fp);
        }
    }

    if (outsave == 1)
    {
        char fname[STRINGMAXLEN_FULLFILENAME];
        PIAACMCsimul_update_fnamedescr();
        SNPRINTF_CHECK(fname,
                       STRINGMAXLEN_FULLFILENAME,
                       "%s/contrast_ptsrc_sm%d.%s.txt",
                       piaacmcparams.piaacmcconfdir,
                       piaacmcparams.SCORINGMASKTYPE,
                       piaacmcparams.fnamedescr);
        printf("saving contrast value [%g] -> %s\n", avContrast, fname);
        FILE *fp = fopen(fname, "w");
        fprintf(fp, "%g", avContrast);
        fclose(fp);
    }

    *avContrast_out = avContrast;
    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

/**
 * PIAACMCsimul_computePSF - Evaluate polychromatic PSF and dark-hole mean contrast.
 * @xld:         Source X pointing offset [lambda/D]
 * @yld:         Source Y pointing offset [lambda/D]
 * @startelem:   First optical element index in propagation
 * @endelem:     Unused legacy ending element parameter
 * @savepsf:     Unused legacy save flag
 * @sourcesize:  Stellar angular size code (`0` = point source, `>0` = `10^(-0.1*code)`)
 * @extmode:     Extended-source quadrature (`0`: 1-pt, `1`: 3-pt, `2`: 6-pt)
 * @outsave:     Save output FITS/flux/contrast files when `1`
 * @contrastval: Optional output pointer receiving average scoring contrast
 *
 * Return: RETURN_SUCCESS on success, or propagation/I/O error status.
 */
errno_t PIAACMCsimul_computePSF(
    float   xld,
    float   yld,
    long    startelem,
    long    endelem,
    int     savepsf,
    int     sourcesize,
    int     extmode,
    int     outsave,
    double *contrastval)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %f %f", xld, yld);

    (void) savepsf;
    (void) endelem;

    uint32_t size      = piaacmcopticaldesign.size;
    imageID  IDopderrC = image_ID("OPDerrC", dcimg, dcnimg);
    long     nbOPDerr  = 0;

    if (IDopderrC == -1)
    {
        load_fits("OPDerrC.fits", "OPDerrC", 0, &IDopderrC);
    }
    if (IDopderrC != -1)
    {
        if (dcimg[IDopderrC].md[0].naxis == 2)
        {
            nbOPDerr = dcimg[IDopderrC].md[0].size[2];
        }
        printf("INCLUDING %ld OPD ERROR MODES\n", nbOPDerr);
        fflush(stdout);
    }

    float focscale =
        (2.0 * piaacmcopticaldesign.beamrad / piaacmcopticaldesign.pixscale) /
        piaacmcopticaldesign.size;

    FUNC_CHECK_RETURN(PIAACMCsimul_init_scoringmask(size, focscale));

    double avContrast = 0.0;
    if (piaacmcparams.computePSF_FAST_FPMresp == 1)
    {
        FUNC_CHECK_RETURN(PIAACMCsimul_eval_fast_fpmresp(size, focscale, &avContrast));
    }
    else if (sourcesize != 0)
    {
        long NBimindex = 0;

        printf("COMPUTING RESOLVED SOURCE PSF / ADDING OPD MODES\n");
        fflush(stdout);

        FUNC_CHECK_RETURN(propagate_extended_quadrature(xld,
                                                        yld,
                                                        startelem,
                                                        sourcesize,
                                                        extmode,
                                                        IDopderrC,
                                                        nbOPDerr,
                                                        &NBimindex));
        if (outsave == 1)
        {
            char fname[STRINGMAXLEN_FULLFILENAME];
            PIAACMCsimul_update_fnamedescr();
            SNPRINTF_CHECK(fname,
                           STRINGMAXLEN_FULLFILENAME,
                           "%s/psfi0_exsrc%03d_sm%d.%s.fits",
                           piaacmcparams.piaacmcconfdir,
                           sourcesize,
                           piaacmcparams.SCORINGMASKTYPE,
                           piaacmcparams.fnamedescr);
            save_fits("psfi0ext", fname);
        }
        assemble_extended_imvect(NBimindex);
        finalize_extended_contrast(size, focscale, sourcesize, outsave, &avContrast);
    }
    else
    {
        FUNC_CHECK_RETURN(
            compute_psf_point_source(xld, yld, startelem, size, focscale, outsave, &avContrast));
    }

    if (contrastval != NULL)
    {
        *contrastval = avContrast;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
