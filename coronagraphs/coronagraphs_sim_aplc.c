/**
 * @file    coronagraphs_sim_aplc.c
 * @brief   Multi-step APLC/PIAACMC setup and external occulter simulators
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "coffee_compat.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "WFpropagate/WFpropagate.h"
#include "fft/fft.h"
#include "image_gen/image_gen.h"

#include "coronagraphs_internal.h"

/**
 * aplc_load_table_if_missing - Load an APLC lookup FITS image into memory if not cached
 * @im_name:   Shared-memory image name
 * @fits_file: FITS filename relative to CORONAGRAPHSDATADIR
 * @slot:      In/out imageID slot
 */
static void aplc_load_table_if_missing(
    const char *im_name,
    const char *fits_file,
    long       *slot)
{
    *slot = image_ID(im_name, dcimg, dcnimg);
    if (*slot == -1)
    {
        char fname[200];
        snprintf(fname, sizeof(fname), "%s/%s", CORONAGRAPHSDATADIR, fits_file);
        load_fits(fname, im_name, 1, slot);
    }
}

/**
 * aplc_load_apo_fitparams - Load APLC lookup cubes and populate analytical fit coefficients
 *
 * Return: Ideal focal-plane complex amplitude multiplier (MASKCAMULT).
 */
static double aplc_load_apo_fitparams(void)
{
    aplc_load_table_if_missing("pinit", "APLCapo_init.fits.gz", &g_coro.IDprol_init);
    aplc_load_table_if_missing("pffrac", "APLCapo_ffrac.fits.gz", &g_coro.IDprol_ffrac);
    aplc_load_table_if_missing("ptransm", "APLCapo_transm.fits.gz", &g_coro.IDprol_transm);
    aplc_load_table_if_missing("ppeak", "APLCapo_peak.fits.gz", &g_coro.IDprol_peak);
    aplc_load_table_if_missing("fitapoa", "APLCapo_fitapo_a.fits.gz", &g_coro.IDprol_fitapo_a);
    aplc_load_table_if_missing("fitapob", "APLCapo_fitapo_b.fits.gz", &g_coro.IDprol_fitapo_b);
    aplc_load_table_if_missing("fitapoc", "APLCapo_fitapo_c.fits.gz", &g_coro.IDprol_fitapo_c);
    aplc_load_table_if_missing("fitfit", "APLCapo_fitfit.fits.gz", &g_coro.IDprol_fitfit);

    long iisize = dcimg[g_coro.IDprol_init].md[0].size[0];
    long jjsize = dcimg[g_coro.IDprol_init].md[0].size[1];
    long iico =
        (long) ((g_coro.APLC_CentOBS1 - g_coro.APLCapo_CO_START) / g_coro.APLCapo_CO_STEP + 0.1);
    long jjfpmrad = (long) ((g_coro.APLC_FPMASKsize - g_coro.APLCapo_FPMRAD_START) /
                                g_coro.APLCapo_FPMRAD_STEP +
                            0.1);

    printf("%ld %ld\n", iico, jjfpmrad);
    if ((iico < 0) || (iico > iisize - 1))
    {
        printf("ERROR: iico out of range : %ld / %ld\n", iico, iisize);
        exit(0);
    }
    if ((jjfpmrad < 0) || (jjfpmrad > jjsize - 1))
    {
        printf("ERROR: iico out of range : %ld / %ld\n", jjfpmrad, jjsize);
        exit(0);
    }

    long cell = jjfpmrad * iisize + iico;
    if (dcimg[g_coro.IDprol_init].array.F[cell] <= 0.1f)
    {
        printf("NO DATA FOR THIS POINT\n");
        exit(0);
    }

    double tmp4         = dcimg[g_coro.IDprol_peak].array.F[cell];
    double mask_ca_mult = -(1.0 - tmp4) / tmp4;
    printf("Focal plane mask CA = %f (peak = %f)\n", mask_ca_mult, tmp4);

    if (g_coro.fitapoINIT == 0)
    {
        g_coro.fitapoN    = dcimg[g_coro.IDprol_fitapo_a].md[0].size[2];
        g_coro.fitapo_a   = (double *) malloc(sizeof(double) * g_coro.fitapoN);
        g_coro.fitapo_b   = (double *) malloc(sizeof(double) * g_coro.fitapoN);
        g_coro.fitapo_c   = (double *) malloc(sizeof(double) * g_coro.fitapoN);
        g_coro.fitapo_c1  = (double *) malloc(sizeof(double) * g_coro.fitapoN);
        g_coro.fitapoINIT = 1;
    }

    for (long n = 0; n < g_coro.fitapoN; n++)
    {
        long idx            = n * iisize * jjsize + cell;
        g_coro.fitapo_a[n]  = dcimg[g_coro.IDprol_fitapo_a].array.F[idx];
        g_coro.fitapo_b[n]  = dcimg[g_coro.IDprol_fitapo_b].array.F[idx];
        g_coro.fitapo_c[n]  = dcimg[g_coro.IDprol_fitapo_c].array.F[idx];
        g_coro.fitapo_c1[n] = (g_coro.fitapo_c[n] < g_coro.fitapo_minc) ? g_coro.fitapo_minc
                                                                        : g_coro.fitapo_c[n];
    }

    printf(
        "f(x) = %15.10f * exp(%15.10f * x**%15.10f)",
        g_coro.fitapo_a[0],
        g_coro.fitapo_b[0],
        g_coro.fitapo_c1[0]);
    for (long n = 1; n < g_coro.fitapoN; n++)
    {
        printf(
            "+ %15.10f * exp(%15.10f * x**%15.10f)",
            g_coro.fitapo_a[n],
            g_coro.fitapo_b[n],
            g_coro.fitapo_c1[n]);
    }
    printf("\n");

    return mask_ca_mult;
}

/**
 * aplc_setup_fpm - Ensure focal plane mask image "aplcfpm" exists in shared memory
 * @mask_ca_mult: Ideal focal-plane complex amplitude multiplier
 */
static void aplc_setup_fpm(double mask_ca_mult)
{
    const uint32_t size   = CORONAGRAPHS_ARRAYSIZE;
    const uint64_t size2  = (uint64_t) size * size;
    imageID        id_fpm = image_ID("aplcfpm", dcimg, dcnimg);

    if (id_fpm != -1)
    {
        return;
    }

    if (g_coro.FPMASKSIZE_ERROR == 0)
    {
        double rad_pix = g_coro.DFTZFACTOR * g_coro.APLC_FPMASKsize / CORONAGRAPHS_PIXSCALE;
        printf("MASK RADIUS = %f pix\n", rad_pix);
        id_fpm      = make_subpixdisk("aplcfpm", size, size, size / 2, size / 2, rad_pix);
        double mult = 1.0 - mask_ca_mult * (1.0 + g_coro.FPMASK_transm_error);
        for (uint64_t ii = 0; ii < size2; ii++)
        {
            dcimg[id_fpm].array.F[ii] = (float) (1.0 - mult * dcimg[id_fpm].array.F[ii]);
        }
        save_fl_fits("aplcfpm", "FPmask.fits");
        return;
    }

    const int fpm_zones = 1;
    printf(
        "-------------- FPMASKSIZE_ERROR = 1, FPMASK_FACTOR = %f, "
        "FPMZONES = %d ----------- (%f %f)\n",
        g_coro.FPMASK_FACTOR,
        fpm_zones,
        g_coro.FPMASK_transm_error,
        g_coro.FPMASK_size_error);

    double rad_pix = g_coro.DFTZFACTOR * g_coro.APLC_FPMASKsize * g_coro.FPMASK_FACTOR /
                     CORONAGRAPHS_PIXSCALE * (1.0 + g_coro.FPMASK_size_error);
    make_subpixdisk("aplcfpm", size, size, size / 2, size / 2, rad_pix);
    save_fl_fits("aplcfpm", "fpm.test.fits");
}

/**
 * aplc_setup_entrance_pupil - Construct and apodize entrance pupil into complex field "pc1"
 * @xld:         Source X offset [lambda/D]
 * @yld:         Source Y offset [lambda/D]
 * @write_files: Non-zero to write intermediate diagnostic FITS files
 *
 * Return: Unapodized entrance pupil total energy (total0).
 */
static double aplc_setup_entrance_pupil(
    double xld,
    double yld,
    int    write_files)
{
    const uint32_t size     = CORONAGRAPHS_ARRAYSIZE;
    const uint64_t size2    = (uint64_t) size * size;
    const double   trad_pix = coro_pupil_radius_pix();

    printf("MAKE INPUT PUPIL\n");
    printf("RADIUS = %f\n", trad_pix);
    imageID id_a1   = make_disk("pa1", size, size, size / 2, size / 2, trad_pix + 2.0);
    imageID id_pa0  = make_disk("pa0", size, size, size / 2, size / 2, trad_pix);
    imageID id_p1   = make_slopexy(
        "pp1",
        size,
        size,
        CORO_PI * xld / trad_pix,
        CORO_PI * yld / trad_pix);

    imageID id_pha = image_ID("corpha", dcimg, dcnimg);
    if (id_pha != -1)
    {
        printf("Applying pre phase errors\n");
        for (uint64_t ii = 0; ii < size2; ii++)
        {
            dcimg[id_p1].array.F[ii] += dcimg[id_pha].array.F[ii];
        }
    }

    imageID id_amp = image_ID("coramp", dcimg, dcnimg);
    if (id_amp != -1)
    {
        printf("Applying pre amplitude errors\n");
        for (uint64_t ii = 0; ii < size2; ii++)
        {
            dcimg[id_a1].array.F[ii] *= dcimg[id_amp].array.F[ii];
        }
    }

    double total0 = 0.0;
    for (uint64_t ii = 0; ii < size2; ii++)
    {
        total0 += dcimg[id_a1].array.F[ii] * dcimg[id_a1].array.F[ii];
    }

    if (g_coro.APLC_PIAA == 0)
    {
        imageID id_apomask = image_ID("apomask", dcimg, dcnimg);
        if (id_apomask != -1)
        {
            printf("Applying amplitude apodization mask\n");
            for (uint64_t ii = 0; ii < size2; ii++)
            {
                dcimg[id_a1].array.F[ii] *= dcimg[id_apomask].array.F[ii];
            }
        }
        mk_complex_from_amph("pa1", "pp1", "pc1", 0);
        delete_image_ID("pa1", DELETE_IMAGE_ERRMODE_WARNING);
        delete_image_ID("pp1", DELETE_IMAGE_ERRMODE_WARNING);

        if (id_amp != -1)
        {
            for (long ii = 3; ii < (long) size - 3; ii++)
            {
                for (long jj = 3; jj < (long) size - 3; jj++)
                {
                    if (dcimg[id_amp].array.F[jj * size + ii] < 0.001f)
                    {
                        dcimg[id_pa0].array.F[jj * size + ii] = 0.0f;
                    }
                }
            }
        }
        return total0;
    }

    g_coro.PIAAAPO_NBPOINTS = (long) trad_pix;
    if (g_coro.initPIAA == 0)
    {
        coronagraph_init_PIAA(NULL);
        g_coro.initPIAA = 1;
    }

    coronagraphs_PIAA_apodize_beam("pa1", "pp1", "pa1a", "pp1a");
    if (write_files == 1)
    {
        save_fl_fits("pa1", "pa1.fits");
        save_fl_fits("pp1", "pp1.fits");
        save_fl_fits("pa1a", "pa1a.fits");
        save_fl_fits("pp1a", "pp1a.fits");
    }
    delete_image_ID("pa1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp1", DELETE_IMAGE_ERRMODE_WARNING);

    imageID id_pa1a = image_ID("pa1a", dcimg, dcnimg);
    imageID id_post = image_ID("postPIAAmask", dcimg, dcnimg);
    if (id_post != -1)
    {
        printf("Applying postPIAA amplitude mask\n");
        for (uint64_t ii = 0; ii < size2; ii++)
        {
            dcimg[id_pa1a].array.F[ii] *= dcimg[id_post].array.F[ii];
        }
        if (write_files == 1)
        {
            save_fl_fits("pa1a", "pa1a_post.fits");
            save_fl_fits("pp1a", "pp1a_post.fits");
        }
    }

    mk_complex_from_amph("pa1a", "pp1a", "pc1", 0);
    if (write_files == 1)
    {
        save_fl_fits("pa1a", "pa1a");
        save_fl_fits("pp1a", "pp1a");
    }
    delete_image_ID("pa1a", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp1a", DELETE_IMAGE_ERRMODE_WARNING);
    return total0;
}

/**
 * coronagraph_simul_MULTISTEP_APLC - Simulate multi-step APLC / PIAACMC coronagraph PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output PSF image name
 *
 * Return: RETURN_SUCCESS on completion.
 */
errno_t coronagraph_simul_MULTISTEP_APLC(
    double      xld,
    double      yld,
    const char *psfname)
{
    const uint32_t size        = CORONAGRAPHS_ARRAYSIZE;
    int            write_files = 1;
    imageID        id_var      = variable_ID("DFTZFACTOR");

    if (id_var != -1)
    {
        g_coro.DFTZFACTOR = coffee_dcvar_f(id_var);
    }
    if (g_coro.useDFT == 0)
    {
        g_coro.DFTZFACTOR = 1.0;
    }
    else
    {
        printf("DFTZFACTOR = %f\n", g_coro.DFTZFACTOR);
    }

    if (variable_ID("WRITEFILES") != -1)
    {
        write_files = 1;
    }

    printf("SIZE = %ul\n", size);
    fflush(stdout);

    double mask_ca_mult = aplc_load_apo_fitparams();
    aplc_setup_fpm(mask_ca_mult);
    double total0 = aplc_setup_entrance_pupil(xld, yld, write_files);

    coro_aplc_propagate_steps(write_files, mask_ca_mult);
    coro_aplc_finish_psf(psfname, total0, write_files);

    delete_image_ID("pc1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa0", DELETE_IMAGE_ERRMODE_WARNING);
    return RETURN_SUCCESS;
}

/**
 * coronagraph_init_EXTERNAL_OCCULTER - Fresnel-propagate hyper-Gaussian starshade to telescope
 * @D:      Telescope aperture diameter [m]
 * @l:      Starshade separation distance [m]
 * @lambda: Optical wavelength [m]
 * @FACTOR: Grid oversize factor relative to CORONAGRAPHS_ARRAYSIZE
 *
 * Return: Image identifier of propagated complex shadow field "extoccp".
 */
int coronagraph_init_EXTERNAL_OCCULTER(
    double D,
    double l,
    double lambda,
    long   FACTOR)
{
    const uint32_t size     = CORONAGRAPHS_ARRAYSIZE;
    const double   param_a  = 12.5;
    const double   param_b  = 12.5;
    const long     param_n  = 6;
    const double   trad_pix = coro_pupil_radius_pix();
    const double   pixscale = D / 2.0 / trad_pix;

    make_offsetHyperGaussian(
        size * FACTOR,
        trad_pix / (D / 2.0) * param_a,
        trad_pix / (D / 2.0) * param_b,
        param_n,
        "extocc");
    printf("pix:  %g  %g\n", trad_pix / (D / 2.0) * param_a, trad_pix / (D / 2.0) * param_b);
    save_fl_fits("extocc", "extocc");
    create_2Dimage_ID("zero", size * FACTOR, size * FACTOR, NULL);
    mk_complex_from_amph("extocc", "zero", "extocc_c", 0);
    printf("pixscale = %g    l = %g    lambda = %g\n", pixscale, l, lambda);
    Fresnel_propagate_wavefront("extocc_c", "extoccp", pixscale, l, lambda);
    imageID id_occ = image_ID("extoccp", dcimg, dcnimg);

    delete_image_ID("zero", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("extocc", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("extocc_c", DELETE_IMAGE_ERRMODE_WARNING);
    return (int) id_occ;
}

/**
 * coronagraph_simul_EXTERNAL_OCCULTER - Simulate external starshade occulter PSF
 * @xld:     Source X offset [lambda/D]
 * @yld:     Source Y offset [lambda/D]
 * @psfname: Output PSF image name
 *
 * Return: 0 on completion.
 */
int coronagraph_simul_EXTERNAL_OCCULTER(
    double      xld,
    double      yld,
    const char *psfname)
{
    const long   size     = CORONAGRAPHS_ARRAYSIZE;
    const long   size2    = size * size;
    const double D        = 6.5;
    const double l        = 50000000.0;
    const double lambda   = 0.0000005;
    const long   factor   = 2;
    const double trad_pix = coro_pupil_radius_pix();

    if (image_ID("extoccp", dcimg, dcnimg) == -1)
    {
        coronagraph_init_EXTERNAL_OCCULTER(D, l, lambda, factor);
    }

    mk_amph_from_complex("extoccp", "exta", "extp", 0);
    imageID id_occa = image_ID("exta", dcimg, dcnimg);
    imageID id_occp = image_ID("extp", dcimg, dcnimg);

    double pixscale        = D / 2.0 / trad_pix;
    double shadow_pixscale = pixscale / l;
    double xrad            = xld * lambda / D;
    double yrad            = yld * lambda / D;
    printf("rad = %g %g   [%g]\n", xrad, yrad, shadow_pixscale);

    long xtransl = (long) (xrad / shadow_pixscale);
    long ytransl = (long) (yrad / shadow_pixscale);
    printf("transl = %ld %ld\n", xtransl, ytransl);

    imageID id_a1;
    create_2Dimage_ID("pa1", size, size, &id_a1);
    imageID id_p1 = make_slopexy(
        "pp1",
        size,
        size,
        CORO_PI * xld / trad_pix,
        CORO_PI * yld / trad_pix);

    double total0 = 0.0;
    for (long ii = 0; ii < size; ii++)
    {
        for (long jj = 0; jj < size; jj++)
        {
            long ii1 = (ii - size / 2) + size * factor / 2 - xtransl;
            long jj1 = (jj - size / 2) + size * factor / 2 - ytransl;
            long idx = jj * size + ii;

            if ((ii1 > 0) && (ii1 < size * factor) && (jj1 > 0) && (jj1 < size * factor))
            {
                dcimg[id_a1].array.F[idx] = dcimg[id_occa].array.F[jj1 * size * factor + ii1];
                dcimg[id_p1].array.F[idx] += dcimg[id_occp].array.F[jj1 * size * factor + ii1];
            }
            else
            {
                dcimg[id_a1].array.F[idx] = 1.0f;
            }

            double dx = 1.0 * (ii - size / 2);
            double dy = 1.0 * (jj - size / 2);
            if (sqrt(dx * dx + dy * dy) > trad_pix)
            {
                dcimg[id_a1].array.F[idx] = 0.0f;
            }
            else
            {
                total0 += 1.0;
            }
        }
    }

    mk_complex_from_amph("pa1", "pp1", "pc1", 0);
    do2dfft("pc1", "fc2");
    permut("fc2");
    mk_amph_from_complex("fc2", "fa2", "fp2", 0);

    imageID id_fa2 = image_ID("fa2", dcimg, dcnimg);
    delete_image_ID("fp2", DELETE_IMAGE_ERRMODE_WARNING);
    imageID id_out;
    create_2Dimage_ID(psfname, size, size, &id_out);
    for (long ii = 0; ii < size2; ii++)
    {
        float a                   = dcimg[id_fa2].array.F[ii];
        dcimg[id_out].array.F[ii] = (float) (a / size2 / total0 * a);
    }

    delete_image_ID("exta", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("extp", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pa1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pp1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pc1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fc2", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fa2", DELETE_IMAGE_ERRMODE_WARNING);
    return 0;
}
