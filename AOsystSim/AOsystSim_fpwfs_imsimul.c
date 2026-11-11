/**
 * @file    AOsystSim_fpwfs_imsimul.c
 * @brief   GSL simplex objective and focal-plane WFS probed PSF cube simulation
 */

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "fft/fft.h"
#include "image_gen/image_gen.h"
#include "statistic/statistic.h"

#include "AOsystSim/AOsystSim.h"
#include "AOsystSim_internal.h"

/* --- Shared FPWFS Module Globals --- */

long        NBprobesG   = 3;
int         CENTERprobe = 1;
long        NBoptVar;
long double probe_re[100];
long double probe_im[100];
long double nprobe_re[100];
long double nprobe_im[100];
long double probe_tflux[100];
long double Cflux = 1.0;
long double probe_nmflux[100];
long double probe_nmnoise[100];
double      tmpvalue1;
double      tmpvalue2;

/**
 * @brief Local state for focal-plane WFS probed image cube simulation.
 */
typedef struct
{
    long    size;
    double  puprad;
    double  CPAxmin;
    double  CPAxmax;
    double  CPAymin;
    double  CPAymax;
    double  CPAstep;
    double  tot1;
    double  peak;
    imageID IDpupa;
    imageID IDwf0;
    imageID IDwfA;
    imageID IDwfB;
    imageID IDwf;
    imageID IDpsfC;
} FPWFS_IMSIMUL_CTX;

/**
 * f_eval - Evaluate regularized chi-squared objective for FPWFS complex-amplitude fit
 * @v:      Candidate parameter vector (ptre, ptim, Iflux, [are, [aim, e]])
 * @params: Unused GSL callback context pointer
 *
 * Return: Objective cost value.
 */
double f_eval(
    const gsl_vector *v,
    void             *params)
{
    (void) params;

    long double ptre_test  = gsl_vector_get(v, 0);
    long double ptim_test  = gsl_vector_get(v, 1);
    long double Iflux_test = gsl_vector_get(v, 2);
    long double are_test   = (NBoptVar > 3) ? gsl_vector_get(v, 3) : 1.0;
    long double aim_test   = (NBoptVar > 4) ? gsl_vector_get(v, 4) : 0.0;
    long double e_test     = (NBoptVar > 4) ? gsl_vector_get(v, 5) : 1.0;

    long double value = 0.0;
    for (long pr = 0; pr < NBprobesG; pr++)
    {
        long double re = probe_re[pr] - ptre_test;
        long double im = probe_im[pr] - ptim_test;
        if (NBoptVar > 3)
        {
            long double x   = re * are_test + im * aim_test;
            long double y   = (-re * aim_test + im * are_test) * e_test;
            probe_tflux[pr] = Iflux_test + Cflux * (x * x + y * y);
        }
        else
        {
            probe_tflux[pr] = Iflux_test + Cflux * (re * re + im * im);
        }

        long double tmp1 = (probe_tflux[pr] - probe_nmflux[pr]) / probe_nmnoise[pr];
        value += tmp1 * tmp1;
    }

    tmpvalue1 = (double) value;
    value += pow(0.01 * (double) (ptre_test * ptre_test), 4.0);
    value += pow(0.01 * (double) (ptim_test * ptim_test), 4.0);

    if (NBoptVar > 4)
    {
        long double ai = are_test * are_test + aim_test * aim_test;
        value += 0.0001 * pow((double) ((ai - 1.0) * (ai - 1.0)), 4.0);
        tmpvalue2 = (double) value;
        value += 0.01 * pow((double) ((e_test - 1.0) * (e_test - 1.0)), 8.0);
    }

    return (double) value;
}

/**
 * fpwfs_imsimul_init_probes - Build apodized pupil, static WF error, and DM probe maps
 * @c:         Simulation state structure
 * @probeamp:  DM probe amplitude in radians
 * @wferramp:  Static wavefront error amplitude in radians
 * @DMgainErr: Fractional random actuator gain error
 */
static void fpwfs_imsimul_init_probes(
    FPWFS_IMSIMUL_CTX *c,
    double             probeamp,
    double             wferramp,
    double             DMgainErr)
{
    const long size     = c->size;
    const double puprad = c->puprad;

    c->IDpupa = make_subpixdisk("pupa", size, size, 0.5 * size, 0.5 * size, puprad);
    c->IDwf0  = image_ID("wf0", dcimg, dcnimg);
    if (c->IDwf0 == -1)
    {
        create_2Dimage_ID("wf0", size, size, &c->IDwf0);
        for (long ii = 0; ii < size * size; ii++)
        {
            dcimg[c->IDwf0].array.F[ii] = (float) (wferramp * (1.0 - 2.0 * ran1()));
        }
        save_fl_fits("wf0", "AOsystSim_wdir/wf0.fits");
    }

    create_2Dimage_ID("wfA", size, size, &c->IDwfA);
    create_2Dimage_ID("wfB", size, size, &c->IDwfB);
    create_2Dimage_ID("wf", size, size, &c->IDwf);
    create_3Dimage_ID("psfC", size, size, NBprobesG, &c->IDpsfC);

    const double step2  = probeamp * c->CPAstep * c->CPAstep;
    const double b_offs = 0.5 * M_PI;
    for (long ii = 0; ii < size; ii++)
    {
        for (long jj = 0; jj < size; jj++)
        {
            double x = (1.0 * ii - 0.5 * size) / puprad;
            double y = (1.0 * jj - 0.5 * size) / puprad;
            for (double CPAx = c->CPAxmin; CPAx < c->CPAxmax; CPAx += c->CPAstep)
            {
                for (double CPAy = c->CPAymin; CPAy < c->CPAymax; CPAy += c->CPAstep)
                {
                    double arg = M_PI * (x * CPAx + y * CPAy);
                    dcimg[c->IDwfA].array.F[jj * size + ii] += (float) (step2 * cos(arg));
                    dcimg[c->IDwfB].array.F[jj * size + ii] +=
                        (float) (step2 * cos(arg + b_offs));
                }
            }
            double gain  = 1.0 + DMgainErr * (1.0 - 2.0 * ran1());
            float  pmask = dcimg[c->IDpupa].array.F[jj * size + ii];
            dcimg[c->IDwfA].array.F[jj * size + ii] *= (float) gain * pmask;
            dcimg[c->IDwfB].array.F[jj * size + ii] *= (float) gain * pmask;
            dcimg[c->IDpupa].array.F[jj * size + ii] *= (float) exp(-8.0 * (x * x + y * y));
        }
    }
    save_fl_fits("wfA", "AOsystSim_wdir/wfA.fits");
    save_fl_fits("wfB", "AOsystSim_wdir/wfB.fits");
    save_fl_fits("pupa", "AOsystSim_wdir/pupa.fits");
}

/**
 * fpwfs_imsimul_compute_psfs - Propagate all probed pupil wavefronts to focal-plane PSFs
 * @c: Simulation state structure
 */
static void fpwfs_imsimul_compute_psfs(FPWFS_IMSIMUL_CTX *c)
{
    const long size  = c->size;
    const long size2 = size * size;
    char       fname[200];

    c->tot1 = 0.0;
    for (long pr = 0; pr < NBprobesG; pr++)
    {
        double coeffA = 0.0;
        double coeffB = 0.0;
        if (pr != 0)
        {
            double ang = 2.0 * M_PI / (NBprobesG - 1) * (pr - 1);
            coeffA     = cos(ang);
            coeffB     = sin(ang);
        }

        for (long ii = 0; ii < size2; ii++)
        {
            dcimg[c->IDwf].array.F[ii] =
                dcimg[c->IDwf0].array.F[ii] +
                (float) coeffA * dcimg[c->IDwfA].array.F[ii] +
                (float) coeffB * dcimg[c->IDwfB].array.F[ii];
        }

        snprintf(fname, sizeof(fname), "AOsystSim_wdir/DMprobe%02ld.fits", pr);
        save_fl_fits("wf", fname);

        mk_complex_from_amph("pupa", "wf", "wfc", 0);
        permut("wfc");
        do2dfft("wfc", "imc");
        permut("imc");
        mk_amph_from_complex("imc", "ima", "imp", 0);
        delete_image_ID("imc", DELETE_IMAGE_ERRMODE_WARNING);
        delete_image_ID("imp", DELETE_IMAGE_ERRMODE_WARNING);

        imageID IDa = image_ID("ima", dcimg, dcnimg);
        for (long ii = 0; ii < size2; ii++)
        {
            float amp = dcimg[IDa].array.F[ii];
            float psf = amp * amp;
            dcimg[c->IDpsfC].array.F[pr * size2 + ii] = psf;
            c->tot1 += psf;
        }
        delete_image_ID("ima", DELETE_IMAGE_ERRMODE_WARNING);
    }
}

/**
 * fpwfs_imsimul_add_companions - Inject five off-axis companions and normalize total flux
 * @c:        Simulation state structure
 * @sepx:     Primary companion X offset in pixels
 * @sepy:     Primary companion Y offset in pixels
 * @contrast: Primary companion contrast ratio
 * @totFlux:  Total integrated flux across PSF cube in photons
 */
static void fpwfs_imsimul_add_companions(
    FPWFS_IMSIMUL_CTX *c,
    double             sepx,
    double             sepy,
    double             contrast,
    double             totFlux)
{
    static const long   dx_comp[5]  = {0, -10, 10, 10, -10};
    static const long   dy_comp[5]  = {0, 10, 10, -10, -10};
    static const double rel_comp[5] = {1.0, 0.3, 0.1, 0.03, 0.01};

    const long size  = c->size;
    const long size2 = size * size;
    imageID    IDtmp;

    printf("Adding companions\n");
    fflush(stdout);
    create_3Dimage_ID("tmp3dim", size, size, NBprobesG, &IDtmp);
    for (long ii = 0; ii < size2 * NBprobesG; ii++)
    {
        dcimg[IDtmp].array.F[ii] = dcimg[c->IDpsfC].array.F[ii];
    }

    for (long pr = 0; pr < NBprobesG; pr++)
    {
        for (long ii = 0; ii < size; ii++)
        {
            for (long jj = 0; jj < size; jj++)
            {
                float src = dcimg[IDtmp].array.F[pr * size2 + jj * size + ii];
                for (int cc = 0; cc < 5; cc++)
                {
                    long ii1 = ii + (long) sepx + dx_comp[cc];
                    long jj1 = jj + (long) sepy + dy_comp[cc];
                    if ((ii1 > 0) && (ii1 < size) && (jj1 > 0) && (jj1 < size))
                    {
                        dcimg[c->IDpsfC].array.F[pr * size2 + jj1 * size + ii1] +=
                            (float) (rel_comp[cc] * contrast) * src;
                    }
                }
            }
        }
    }
    delete_image_ID("tmp3dim", DELETE_IMAGE_ERRMODE_WARNING);

    float scale = (float) (totFlux / c->tot1);
    for (long ii = 0; ii < size2 * NBprobesG; ii++)
    {
        dcimg[c->IDpsfC].array.F[ii] *= scale;
    }

    c->peak = 0.0;
    for (long ii = 0; ii < size2; ii++)
    {
        if (dcimg[c->IDpsfC].array.F[ii] > c->peak)
        {
            c->peak = dcimg[c->IDpsfC].array.F[ii];
        }
    }
}

/**
 * fpwfs_imsimul_crop_probes - Crop dark-hole region and compute probe amplitude map
 * @c:         Simulation state structure
 * @xsize_out: Output cropped X dimension
 * @ysize_out: Output cropped Y dimension
 */
static void fpwfs_imsimul_crop_probes(
    const FPWFS_IMSIMUL_CTX *c,
    long                    *xsize_out,
    long                    *ysize_out)
{
    double  pixscaleld = c->size / (2.0 * c->puprad);
    long    xmin       = (long) (0.5 * c->size + pixscaleld * c->CPAxmin);
    long    xmax       = (long) (0.5 * c->size + pixscaleld * c->CPAxmax);
    long    ymin       = (long) (0.5 * c->size + pixscaleld * c->CPAymin);
    long    ymax       = (long) (0.5 * c->size + pixscaleld * c->CPAymax);
    long    xsize      = xmax - xmin;
    long    ysize      = ymax - ymin;
    imageID ID;
    imageID ID1;

    printf("Cropping and adding photon noise\n");
    fflush(stdout);

    create_3Dimage_ID("psfCcrop", xsize, ysize, NBprobesG, &ID);
    for (long pr = 0; pr < NBprobesG; pr++)
    {
        for (long ii1 = 0; ii1 < xsize; ii1++)
        {
            for (long jj1 = 0; jj1 < ysize; jj1++)
            {
                long ii = ii1 + xmin;
                long jj = jj1 + ymin;
                dcimg[ID].array.F[pr * xsize * ysize + jj1 * xsize + ii1] =
                    dcimg[c->IDpsfC].array.F[pr * c->size * c->size + jj * c->size + ii];
            }
        }
    }

    create_2Dimage_ID("psfprobeampC", xsize, ysize, &ID1);
    for (long ii1 = 0; ii1 < xsize; ii1++)
    {
        for (long jj1 = 0; jj1 < ysize; jj1++)
        {
            double tot1 = 0.0;
            for (long pr = CENTERprobe; pr < NBprobesG; pr++)
            {
                tot1 += dcimg[ID].array.F[pr * xsize * ysize + jj1 * xsize + ii1] / c->peak;
            }
            tot1 /= (NBprobesG - CENTERprobe);
            dcimg[ID1].array.F[jj1 * xsize + ii1] =
                (float) (tot1 - dcimg[ID].array.F[jj1 * xsize + ii1] / c->peak);
        }
    }
    save_fl_fits("psfprobeampC", "AOsystSim_wdir/psfprobeampC.fits");

    *xsize_out = xsize;
    *ysize_out = ysize;
}

/**
 * fpwfs_imsimul_apply_noise - Add photon/readout noise and export normalized PSF cubes
 * @c:           Simulation state structure
 * @xsize:       Cropped X dimension
 * @ysize:       Cropped Y dimension
 * @RON:         Detector readout noise in photo-electrons
 * @CnoiseFloor: Contrast noise floor
 *
 * Return: Image ID of normalized noisy PSF cube ("psfCcropnC").
 */
static imageID fpwfs_imsimul_apply_noise(
    const FPWFS_IMSIMUL_CTX *c,
    long                     xsize,
    long                     ysize,
    double                   RON,
    double                   CnoiseFloor)
{
    imageID IDnoise;
    imageID ID;
    imageID ID1;
    char    imname[200];
    char    fname[200];

    create_3Dimage_ID("psfCcropnCn", xsize, ysize, NBprobesG, &IDnoise);
    put_poisson_noise("psfCcrop", "psfCcropn");
    ID = image_ID("psfCcropn", dcimg, dcnimg);
    for (long ii1 = 0; ii1 < xsize; ii1++)
    {
        for (long jj1 = 0; jj1 < ysize; jj1++)
        {
            dcimg[ID].array.F[jj1 * xsize + ii1] += (float) (RON * gauss());
        }
    }
    save_fl_fits("psfCcrop", "AOsystSim_wdir/psfCcrop.fits");
    save_fl_fits("psfCcropn", "AOsystSim_wdir/psfCcropn.fits");

    create_3Dimage_ID("psfCcropnC", xsize, ysize, NBprobesG, &ID1);
    for (long pr = 0; pr < NBprobesG; pr++)
    {
        for (long ii1 = 0; ii1 < xsize; ii1++)
        {
            for (long jj1 = 0; jj1 < ysize; jj1++)
            {
                long   idx = pr * xsize * ysize + jj1 * xsize + ii1;
                double raw = dcimg[ID].array.F[idx];
                double val = (raw < 1.0) ? 1.0 : raw;
                double nval = sqrt(val + RON * RON) / c->peak;

                dcimg[ID1].array.F[idx]     = (float) (raw / c->peak);
                dcimg[IDnoise].array.F[idx] =
                    (float) ((nval < CnoiseFloor) ? CnoiseFloor : nval);
            }
        }
    }
    save_fl_fits("psfCcropnC", "AOsystSim_wdir/psfCcropnC.fits");
    printf("Saving psfCcropnCn\n");
    save_fl_fits("psfCcropnCn", "AOsystSim_wdir/psfCcropnCnoise.fits");

    for (long pr = 0; pr < NBprobesG; pr++)
    {
        imageID IDs;
        snprintf(imname, sizeof(imname), "psfC_%03ld", pr);
        create_2Dimage_ID(imname, c->size, c->size, &IDs);
        for (long ii = 0; ii < c->size * c->size; ii++)
        {
            dcimg[IDs].array.F[ii] =
                (float) (dcimg[c->IDpsfC].array.F[pr * c->size * c->size + ii] / c->peak);
        }
        snprintf(fname, sizeof(fname), "AOsystSim_wdir/psfC_%03ld.fits", pr);
        save_fl_fits(imname, fname);
    }

    return ID1;
}

/**
 * AOsystSim_FPWFS_imsimul - Simulate focal-plane WFS probed PSF cubes with companions and noise
 * @probeamp:    DM probe amplitude in radians
 * @sepx:        Primary companion X offset in pixels
 * @sepy:        Primary companion Y offset in pixels
 * @contrast:    Primary companion contrast ratio
 * @wferramp:    Static wavefront error amplitude in radians
 * @totFlux:     Total integrated flux in photons
 * @DMgainErr:   Fractional random DM gain error
 * @RON:         Detector readout noise in photo-electrons
 * @CnoiseFloor: Minimum contrast noise floor
 *
 * Return: Image ID of normalized noisy cropped PSF cube ("psfCcropnC").
 */
long AOsystSim_FPWFS_imsimul(
    double probeamp,
    double sepx,
    double sepy,
    double contrast,
    double wferramp,
    double totFlux,
    double DMgainErr,
    double RON,
    double CnoiseFloor)
{
    FPWFS_IMSIMUL_CTX c = {
        .size    = 256,
        .puprad  = 50.0,
        .CPAxmin = 4.0,
        .CPAxmax = 20.0,
        .CPAymin = 4.0,
        .CPAymax = 20.0,
        .CPAstep = 0.2,
    };
    long xsize = 0;
    long ysize = 0;

    aosim_ensure_wdir();
    printf("Creating images .... Contrast = %g  \n", contrast);
    fflush(stdout);

    fpwfs_imsimul_init_probes(&c, probeamp, wferramp, DMgainErr);
    fpwfs_imsimul_compute_psfs(&c);
    fpwfs_imsimul_add_companions(&c, sepx, sepy, contrast, totFlux);
    fpwfs_imsimul_crop_probes(&c, &xsize, &ysize);
    imageID ID1 = fpwfs_imsimul_apply_noise(&c, xsize, ysize, RON, CnoiseFloor);

    printf("PSFs creation is complete\n");
    fflush(stdout);

    return ID1;
}
