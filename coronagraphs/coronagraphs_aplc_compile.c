/**
 * @file    coronagraphs_aplc_compile.c
 * @brief   Compile precomputed 2D APLC prolate apodizer parameters into FITS lookup cubes
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "coffee_compat.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "statistic/statistic.h"

#include "coronagraphs_internal.h"

/**
 * aplc_read_fitparams_file - Read residual and parameter triplets from a .fitparams file
 * @fname: Path to .fitparams file
 * @fitv:  Output fit residual value
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE on I/O error.
 */
static errno_t aplc_read_fitparams_file(
    const char *fname,
    double     *fitv)
{
    FILE *fp = fopen(fname, "r");
    if (fp == NULL)
    {
        printf("ERROR: file \"%s\" does not exist\n", fname);
        exit(0);
    }

    int fscanfcnt = fscanf(fp, "%lf", fitv);
    if (coro_check_fscanf(fp, fscanfcnt, 1) != RETURN_SUCCESS)
    {
        fclose(fp);
        return RETURN_FAILURE;
    }

    for (long n = 0; n < g_coro.fitapoN; n++)
    {
        fscanfcnt = fscanf(
            fp,
            " %lf %lf %lf",
            &g_coro.fitapo_a[n],
            &g_coro.fitapo_b[n],
            &g_coro.fitapo_c[n]);
        if (coro_check_fscanf(fp, fscanfcnt, 3) != RETURN_SUCCESS)
        {
            fclose(fp);
            return RETURN_FAILURE;
        }
    }

    fclose(fp);
    return RETURN_SUCCESS;
}

/**
 * aplc_refit_profile - Load raw .prof radial profile and run analytical simplex fit
 * @fnameapoprof: Path to raw radial profile file
 * @fnameapofit:  Path to output .fitparams file
 * @fitv:         Output fit residual value
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE on error.
 */
static errno_t aplc_refit_profile(
    const char *fnameapoprof,
    const char *fnameapofit,
    double     *fitv)
{
    const double eps0 = 1.0e-8;

    EXECUTE_SYSTEM_COMMAND("awk '{print $1,$2}' %s > apotmp.prof", fnameapoprof);

    FILE *fp = fopen("apotmp.prof", "r");
    if (fp == NULL)
    {
        printf("ERROR : cannot open file \"apotmp.prof\"\n");
        exit(0);
    }
    printf("FILE \"apotmp.prof\" opened \n");
    fflush(stdout);

    g_coro.aporaw_r = (double *) malloc(sizeof(double) * g_coro.PIAAAPO_NBPOINTS);
    g_coro.aporaw_v = (double *) malloc(sizeof(double) * g_coro.PIAAAPO_NBPOINTS);
    if ((g_coro.aporaw_r == NULL) || (g_coro.aporaw_v == NULL))
    {
        PRINT_ERROR("malloc returns NULL pointer");
        fclose(fp);
        free(g_coro.aporaw_r);
        free(g_coro.aporaw_v);
        return RETURN_FAILURE;
    }

    long cnt = 0;
    for (long i = 0; i < g_coro.PIAAAPO_NBPOINTS; i++)
    {
        double tmpf0     = 0.0;
        double tmpf      = 0.0;
        int    fscanfcnt = fscanf(fp, "%lf %lf\n", &tmpf0, &tmpf);
        if (coro_check_fscanf(fp, fscanfcnt, 2) != RETURN_SUCCESS)
        {
            fclose(fp);
            free(g_coro.aporaw_r);
            free(g_coro.aporaw_v);
            return RETURN_FAILURE;
        }
        if (tmpf > eps0)
        {
            g_coro.aporaw_r[cnt] = tmpf0 / 200.0;
            g_coro.aporaw_v[cnt] = tmpf;
            cnt++;
        }
    }
    fclose(fp);
    g_coro.aporawN = cnt;

    *fitv = coronagraph_apofit(fnameapofit);
    printf("fit value = %g\n", *fitv);

    errno_t rc = aplc_read_fitparams_file(fnameapofit, fitv);
    printf("FIT SOLUTION RESIDUAL = %g\n", *fitv);
    fflush(stdout);

    free(g_coro.aporaw_r);
    free(g_coro.aporaw_v);
    return rc;
}

/**
 * aplc_compile_cell - Compile prolate metadata and analytical fit for grid cell (ii, jj)
 * @ii:     Central obstruction grid index
 * @jj:     Focal plane mask radius grid index
 * @iisize: Grid width (central obstruction axis)
 * @jjsize: Grid height (focal plane mask radius axis)
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE on error.
 */
static errno_t aplc_compile_cell(
    long ii,
    long jj,
    long iisize,
    long jjsize)
{
    const double fitv_limit = 2.0e-5;
    struct stat  buffer;
    char         fnameprolinfo[500];
    char         fnameapoprof[500];
    char         fnameapofit[500];

    double co     = g_coro.APLCapo_CO_START + g_coro.APLCapo_CO_STEP * ii;
    double fpmrad = g_coro.APLCapo_FPMRAD_START + g_coro.APLCapo_FPMRAD_STEP * jj;

    if (co < 0.0001)
    {
        snprintf(
            fnameprolinfo,
            sizeof(fnameprolinfo),
            "%s/APLCapo/APLCapo_%5.3f.4096.ref.info",
            CORONAGRAPHSDATALOCAL,
            fpmrad);
        snprintf(
            fnameapoprof,
            sizeof(fnameapoprof),
            "%s/APLCapo/raw/APLCapo_%.3f.4096.ref.prof",
            CORONAGRAPHSDATALOCAL,
            fpmrad);
    }
    else
    {
        snprintf(
            fnameprolinfo,
            sizeof(fnameprolinfo),
            "%s/APLCapo/APLCapo_%5.3f.%5.3f.4096.ref.info",
            CORONAGRAPHSDATALOCAL,
            fpmrad,
            co);
        snprintf(
            fnameapoprof,
            sizeof(fnameapoprof),
            "%s/APLCapo/raw/APLCapo_%.3f.%.3f.4096.ref.prof",
            CORONAGRAPHSDATALOCAL,
            fpmrad,
            co);
    }

    if ((stat(fnameprolinfo, &buffer) != 0) || (stat(fnameapoprof, &buffer) != 0))
    {
        dcimg[g_coro.IDprol_init].array.F[jj * iisize + ii] = 0.0f;
        return RETURN_SUCCESS;
    }

    printf("\r%5ld/%5ld   %5ld/%5ld  %f %f  ", ii, iisize, jj, jjsize, co, fpmrad);

    FILE  *fp = fopen(fnameprolinfo, "r");
    double tmpf0;
    double tmpf1;
    double tmpf2;
    double tmpf3;
    int    fscanfcnt = fscanf(fp, "%lf %lf %lf %lf\n", &tmpf0, &tmpf1, &tmpf2, &tmpf3);
    if (coro_check_fscanf(fp, fscanfcnt, 4) != RETURN_SUCCESS)
    {
        fclose(fp);
        return RETURN_FAILURE;
    }
    fclose(fp);

    dcimg[g_coro.IDprol_init].array.F[jj * iisize + ii]   = 1.0f;
    dcimg[g_coro.IDprol_ffrac].array.F[jj * iisize + ii]  = (float) tmpf1;
    dcimg[g_coro.IDprol_transm].array.F[jj * iisize + ii] = (float) tmpf2;
    dcimg[g_coro.IDprol_peak].array.F[jj * iisize + ii]   = (float) tmpf3;

    snprintf(
        fnameapofit,
        sizeof(fnameapofit),
        "%s/APLCapo/APLCapo_%.3f.%.3f.4096.ref.prof.fitparams",
        CORONAGRAPHSDATALOCAL,
        fpmrad,
        co);

    g_coro.fitapo_a  = (double *) malloc(sizeof(double) * g_coro.fitapoN);
    g_coro.fitapo_b  = (double *) malloc(sizeof(double) * g_coro.fitapoN);
    g_coro.fitapo_c  = (double *) malloc(sizeof(double) * g_coro.fitapoN);
    g_coro.fitapo_c1 = (double *) malloc(sizeof(double) * g_coro.fitapoN);
    if (!g_coro.fitapo_a || !g_coro.fitapo_b || !g_coro.fitapo_c || !g_coro.fitapo_c1)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        return RETURN_FAILURE;
    }

    for (long n = 0; n < g_coro.fitapoN; n++)
    {
        g_coro.fitapo_a[n] = 0.0;
        g_coro.fitapo_b[n] = -1.0;
        g_coro.fitapo_c[n] = 2.0 + n;
    }

    int    fit_ok = 0;
    double fitv   = 0.0;
    if (stat(fnameapofit, &buffer) == 0)
    {
        if (aplc_read_fitparams_file(fnameapofit, &fitv) != RETURN_SUCCESS)
        {
            return RETURN_FAILURE;
        }
        printf("FIT SOLUTION RESIDUAL = %g ", fitv);
        fflush(stdout);
        fit_ok = (fitv <= fitv_limit) ? 1 : 0;
    }
    else
    {
        printf("                          ");
        fflush(stdout);
    }

    if (fit_ok == 0)
    {
        if (aplc_refit_profile(fnameapoprof, fnameapofit, &fitv) != RETURN_SUCCESS)
        {
            return RETURN_FAILURE;
        }
    }

    quick_sort3_double(g_coro.fitapo_a, g_coro.fitapo_b, g_coro.fitapo_c, g_coro.fitapoN);
    for (long n = 0; n < g_coro.fitapoN; n++)
    {
        long idx                                   = n * iisize * jjsize + jj * iisize + ii;
        dcimg[g_coro.IDprol_fitapo_a].array.F[idx] = (float) g_coro.fitapo_a[n];
        dcimg[g_coro.IDprol_fitapo_b].array.F[idx] = (float) g_coro.fitapo_b[n];
        dcimg[g_coro.IDprol_fitapo_c].array.F[idx] = (float) g_coro.fitapo_c[n];
    }
    dcimg[g_coro.IDprol_fitfit].array.F[jj * iisize + ii] = (float) fitv;

    free(g_coro.fitapo_a);
    free(g_coro.fitapo_b);
    free(g_coro.fitapo_c);
    free(g_coro.fitapo_c1);
    return RETURN_SUCCESS;
}

/**
 * aplc_save_and_cleanup - Save compiled APLC parameter cubes to FITS and release images
 */
static void aplc_save_and_cleanup(void)
{
    char fnameout[500];

    snprintf(fnameout, sizeof(fnameout), "%s/APLCapo_init.fits.gz", CORONAGRAPHSDATADIR);
    save_fl_fits("pinit", fnameout);
    snprintf(fnameout, sizeof(fnameout), "%s/APLCapo_ffrac.fits.gz", CORONAGRAPHSDATADIR);
    save_fl_fits("pffrac", fnameout);
    snprintf(fnameout, sizeof(fnameout), "%s/APLCapo_transm.fits.gz", CORONAGRAPHSDATADIR);
    save_fl_fits("ptransm", fnameout);
    snprintf(fnameout, sizeof(fnameout), "%s/APLCapo_peak.fits.gz", CORONAGRAPHSDATADIR);
    save_fl_fits("ppeak", fnameout);
    snprintf(fnameout, sizeof(fnameout), "%s/APLCapo_fitapo_a.fits.gz", CORONAGRAPHSDATADIR);
    save_fl_fits("fitapoa", fnameout);
    snprintf(fnameout, sizeof(fnameout), "%s/APLCapo_fitapo_b.fits.gz", CORONAGRAPHSDATADIR);
    save_fl_fits("fitapob", fnameout);
    snprintf(fnameout, sizeof(fnameout), "%s/APLCapo_fitapo_c.fits.gz", CORONAGRAPHSDATADIR);
    save_fl_fits("fitapoc", fnameout);
    snprintf(fnameout, sizeof(fnameout), "%s/APLCapo_fitfit.fits.gz", CORONAGRAPHSDATADIR);
    save_fl_fits("fitfit", fnameout);

    delete_image_ID("pffrac", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("ptransm", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("ppeak", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fitfit", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fitapoa", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fitapob", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("fitapoc", DELETE_IMAGE_ERRMODE_WARNING);
    printf("\n");
}

/**
 * coronagraph_APLCapo_compile - Compile precomputed APLC prolate apodizer lookup tables
 *
 * Return: RETURN_SUCCESS or RETURN_FAILURE on error.
 */
errno_t coronagraph_APLCapo_compile(void)
{
    long iisize =
        (long) ((g_coro.APLCapo_CO_END - g_coro.APLCapo_CO_START) / g_coro.APLCapo_CO_STEP);
    long jjsize = (long) ((g_coro.APLCapo_FPMRAD_END - g_coro.APLCapo_FPMRAD_START) /
                          g_coro.APLCapo_FPMRAD_STEP);

    create_2Dimage_ID("pinit", iisize, jjsize, &g_coro.IDprol_init);
    create_2Dimage_ID("pffrac", iisize, jjsize, &g_coro.IDprol_ffrac);
    create_2Dimage_ID("ptransm", iisize, jjsize, &g_coro.IDprol_transm);
    create_2Dimage_ID("ppeak", iisize, jjsize, &g_coro.IDprol_peak);

    create_3Dimage_ID("fitapoa", iisize, jjsize, g_coro.fitapoN, &g_coro.IDprol_fitapo_a);
    create_3Dimage_ID("fitapob", iisize, jjsize, g_coro.fitapoN, &g_coro.IDprol_fitapo_b);
    create_3Dimage_ID("fitapoc", iisize, jjsize, g_coro.fitapoN, &g_coro.IDprol_fitapo_c);
    create_2Dimage_ID("fitfit", iisize, jjsize, &g_coro.IDprol_fitfit);

    printf("\n\n");

    for (long ii = 0; ii < iisize; ii++)
    {
        for (long jj = 0; jj < jjsize; jj++)
        {
            if (ran1() < 1.01)
            {
                if (aplc_compile_cell(ii, jj, iisize, jjsize) != RETURN_SUCCESS)
                {
                    return RETURN_FAILURE;
                }
            }
        }
    }

    aplc_save_and_cleanup();
    return RETURN_SUCCESS;
}
