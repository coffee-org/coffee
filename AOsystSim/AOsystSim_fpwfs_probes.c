/**
 * @file    AOsystSim_fpwfs_probes.c
 * @brief   Deformable mirror speckle-nulling phase probe synthesis and test PSF export
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "fft/fft.h"
#include "image_gen/image_gen.h"

#include "AOsystSim/AOsystSim.h"
#include "AOsystSim_internal.h"

/**
 * fpwfs_probe_cpa_bounds - Resolve spatial-frequency domain bounds for requested geometry
 * @modegeom: Geometry selector (0: horiz, 1: vert, 2: quadrant 1, 3: quadrant 2)
 * @CPAmax:   Outer spatial frequency cutoff in cycles per aperture
 * @bounds:   Output 4-element array [CPAxmin, CPAxmax, CPAymin, CPAymax]
 */
static void fpwfs_probe_cpa_bounds(
    long   modegeom,
    double CPAmax,
    double bounds[4])
{
    switch (modegeom)
    {
    case 0:
        bounds[0] = 0.5;
        bounds[1] = CPAmax;
        bounds[2] = -CPAmax;
        bounds[3] = CPAmax;
        break;
    case 1:
        bounds[0] = -CPAmax;
        bounds[1] = CPAmax;
        bounds[2] = 0.5;
        bounds[3] = CPAmax;
        break;
    case 2:
        bounds[0] = 0.5;
        bounds[1] = CPAmax;
        bounds[2] = 0.5;
        bounds[3] = CPAmax;
        break;
    case 3:
        bounds[0] = 0.5;
        bounds[1] = CPAmax;
        bounds[2] = -0.5;
        bounds[3] = -CPAmax;
        break;
    default:
        printf("ERROR: mode not supported\n");
        exit(0);
    }
}

/**
 * fpwfs_save_probe_test_psf - Embed scaled DM probe into pupil phase and export test PSF
 * @IDsrc:   Source DM probe image ID
 * @scale:   Multiplicative scale (+1, -1, or 0) applied to probe phase
 * @dmxsize: DM X dimension
 * @dmysize: DM Y dimension
 * @imsize:  Zero-padded FFT grid size
 * @tag:     Filename suffix ("A", "B", "mA", "mB", "00")
 */
static void fpwfs_save_probe_test_psf(
    imageID     IDsrc,
    float       scale,
    long        dmxsize,
    long        dmysize,
    long        imsize,
    const char *tag)
{
    char    fname_pupp[200];
    char    fname_foca[200];
    imageID ID   = image_ID("pupp", dcimg, dcnimg);
    long    offx = (imsize - dmxsize) / 2;
    long    offy = (imsize - dmysize) / 2;

    for (long ii = 0; ii < dmxsize; ii++)
    {
        for (long jj = 0; jj < dmysize; jj++)
        {
            dcimg[ID].array.F[(jj + offy) * imsize + (ii + offx)] =
                scale * dcimg[IDsrc].array.F[jj * dmxsize + ii];
        }
    }

    mk_complex_from_amph("pupa", "pupp", "pupc", 0);
    permut("pupc");
    do2dfft("pupc", "focc");
    permut("focc");
    mk_amph_from_complex("focc", "foca", "focp", 0);

    snprintf(fname_pupp, sizeof(fname_pupp), "AOsystSim_wdir/test_pupp_%s.fits", tag);
    snprintf(fname_foca, sizeof(fname_foca), "AOsystSim_wdir/test_foca_%s.fits", tag);
    save_fits("pupp", fname_pupp);
    save_fits("foca", fname_foca);

    delete_image_ID("pupc", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("focc", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("foca", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("focp", DELETE_IMAGE_ERRMODE_WARNING);
}

/**
 * fpwfs_export_probe_tests - Build annular test pupil and export +/-A, +/-B, and 00 PSFs
 * @IDdmA:   Quadrature cosine DM probe image ID
 * @IDdmB:   Quadrature sine DM probe image ID
 * @dmxsize: DM X dimension
 * @dmysize: DM Y dimension
 * @dmsize:  Maximum of dmxsize and dmysize
 */
static void fpwfs_export_probe_tests(
    imageID IDdmA,
    imageID IDdmB,
    long    dmxsize,
    long    dmysize,
    long    dmsize)
{
    long    imsize = 5 * dmsize;
    double  rad    = 0.45 * dmsize;
    imageID IDpupa = make_disk("pupa", imsize, imsize, 0.5 * imsize, 0.5 * imsize, rad);
    imageID IDpupp;

    for (long ii = 0; ii < imsize; ii++)
    {
        for (long jj = 0; jj < imsize; jj++)
        {
            double x = (1.0 * ii - 0.5 * imsize) / rad;
            double y = (1.0 * jj - 0.5 * imsize) / rad;
            if (sqrt(x * x + y * y) < 0.3)
            {
                dcimg[IDpupa].array.F[jj * imsize + ii] = 0.0f;
            }
        }
    }
    save_fits("pupa", "AOsystSim_wdir/test_pupa.fits");
    create_2Dimage_ID("pupp", imsize, imsize, &IDpupp);

    fpwfs_save_probe_test_psf(IDdmA, 1.0f, dmxsize, dmysize, imsize, "A");
    fpwfs_save_probe_test_psf(IDdmB, 1.0f, dmxsize, dmysize, imsize, "B");
    fpwfs_save_probe_test_psf(IDdmA, -1.0f, dmxsize, dmysize, imsize, "mA");
    fpwfs_save_probe_test_psf(IDdmB, -1.0f, dmxsize, dmysize, imsize, "mB");
    fpwfs_save_probe_test_psf(IDdmA, 0.0f, dmxsize, dmysize, imsize, "00");
}

/**
 * AOsystSim_FPWFS_mkprobes - Generate quadrature DM probe maps over target CPA band
 * @IDprobeA_name: Output image name for cosine quadrature probe A
 * @IDprobeB_name: Output image name for sine quadrature probe B
 * @dmxsize:       Number of DM actuators along X
 * @dmysize:       Number of DM actuators along Y
 * @CPAmax:        Maximum spatial frequency box limit in cycles per aperture
 * @CPArmin:       Inner radial spatial frequency cutoff in CPA
 * @CPArmax:       Outer radial spatial frequency cutoff in CPA
 * @RMSampl:       Target RMS probe surface amplitude
 * @modegeom:      Dark-hole geometry mode (0: horiz, 1: vert, 2: quad1, 3: quad2)
 *
 * Return: 0 on completion.
 */
int AOsystSim_FPWFS_mkprobes(
    const char *IDprobeA_name,
    const char *IDprobeB_name,
    long        dmxsize,
    long        dmysize,
    double      CPAmax,
    double      CPArmin,
    double      CPArmax,
    double      RMSampl,
    long        modegeom)
{
    const double CPAstep = 0.1;
    double       bounds[4];
    imageID      IDdmA;
    imageID      IDdmB;

    aosim_ensure_wdir();
    create_2Dimage_ID(IDprobeA_name, dmxsize, dmysize, &IDdmA);
    create_2Dimage_ID(IDprobeB_name, dmxsize, dmysize, &IDdmB);
    fpwfs_probe_cpa_bounds(modegeom, CPAmax, bounds);

    long dmsize = (dmysize > dmxsize) ? dmysize : dmxsize;
    for (double CPAx = bounds[0]; CPAx < bounds[1]; CPAx += CPAstep)
    {
        for (double CPAy = bounds[2]; CPAy < bounds[3]; CPAy += CPAstep)
        {
            double CPAr = sqrt(CPAx * CPAx + CPAy * CPAy);
            if ((CPAr <= CPArmin) || (CPAr >= CPArmax))
            {
                continue;
            }
            double pha = -1.6 * CPAx;
            for (long ii = 0; ii < dmxsize; ii++)
            {
                for (long jj = 0; jj < dmysize; jj++)
                {
                    double x   = 2.0 * (1.0 * ii - 0.5 * dmxsize) / dmsize;
                    double y   = 2.0 * (1.0 * jj - 0.5 * dmysize) / dmsize;
                    double arg = M_PI * (x * CPAx + y * CPAy) + pha;
                    dcimg[IDdmA].array.F[jj * dmxsize + ii] += (float) cos(arg);
                    dcimg[IDdmB].array.F[jj * dmxsize + ii] += (float) cos(arg + 0.5 * M_PI);
                }
            }
        }
    }

    double   rms   = 0.0;
    uint64_t npix  = (uint64_t) (dmxsize * dmysize);
    for (uint64_t ii = 0; ii < npix; ii++)
    {
        rms += (double) dcimg[IDdmA].array.F[ii] * dcimg[IDdmA].array.F[ii];
    }
    float norm = (float) (RMSampl / sqrt(rms / (double) npix));
    for (uint64_t ii = 0; ii < npix; ii++)
    {
        dcimg[IDdmA].array.F[ii] *= norm;
        dcimg[IDdmB].array.F[ii] *= norm;
    }

    fpwfs_export_probe_tests(IDdmA, IDdmB, dmxsize, dmysize, dmsize);
    return 0;
}
