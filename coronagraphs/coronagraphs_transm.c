/**
 * @file    coronagraphs_transm.c
 * @brief   Off-axis planetary throughput and IWA evaluation for coronagraph models
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "coffee_compat.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "ZernikePolyn/ZernikePolyn.h"

#include "coronagraphs_internal.h"

typedef struct transm_cfg
{
    long    coro_type;
    double  logcontrast;
    double  contrast_scale;
    double  sepstep;
    double  sepmax;
    double  trad_pix;
    double  tt_error;
    double  foc_error;
    double  z_error;
    long    z_index;
    int     tt_included;
    int     foc_included;
    int     z_included;
    int     radaver;
    int     nb_tt_pts;
    int     nb_foc_pts;
    long    nb_anglestep;
    double  max_angle;
    imageID id_psf0;
    double *arrayS;
    double *arrayS1;
    double *arrayS2;
} transm_cfg_t;

/**
 * transm_extract_opt_token - Extract whitespace-delimited token following an option flag
 * @options: Full option string
 * @flag:    Option prefix including trailing space (e.g. "-TT ")
 * @buf:     Destination buffer
 * @cap:     Destination buffer capacity
 *
 * Return: 1 if found and extracted, 0 otherwise.
 */
static int transm_extract_opt_token(
    const char *options,
    const char *flag,
    char       *buf,
    size_t      cap)
{
    const char *pos = strstr(options, flag);
    if (pos == NULL)
    {
        return 0;
    }

    pos += strlen(flag);
    size_t ii = 0;
    while ((pos[ii] != ' ') && (pos[ii] != '\n') && (pos[ii] != '\0') && (ii + 1 < cap))
    {
        buf[ii] = pos[ii];
        ii++;
    }
    buf[ii] = '\0';
    return 1;
}

/**
 * transm_parse_options - Populate throughput evaluation configuration from option flags
 * @options: Option flags string
 * @cfg:     In/out throughput configuration
 */
static void transm_parse_options(
    const char   *options,
    transm_cfg_t *cfg)
{
    const double epsilon = 1.0e-10;
    char         input[50];

    if (strstr(options, "-RADAVER") != NULL)
    {
        cfg->radaver = 1;
        printf("RADIAL AVERAGE\n");
    }

    if (transm_extract_opt_token(options, "-TT ", input, sizeof(input)))
    {
        cfg->tt_error    = atof(input);
        cfg->tt_included = (cfg->tt_error > epsilon) ? 1 : 0;
        printf("TT error is %f\n", cfg->tt_error);
    }

    if (transm_extract_opt_token(options, "-FPMSERR ", input, sizeof(input)))
    {
        g_coro.FPMASKSIZE_ERROR = 1;
        g_coro.FPMASK_FACTOR    = atof(input);
        printf("--------------- FPMASK_FACTOR = %f ----------------------------\n",
               g_coro.FPMASK_FACTOR);
    }

    if (transm_extract_opt_token(options, "-FOC ", input, sizeof(input)))
    {
        cfg->foc_error    = atof(input);
        cfg->foc_included = (cfg->foc_error > epsilon) ? 1 : 0;
        printf("FOC error is %f  (lambda at edge)\n", cfg->foc_error);
    }

    if (transm_extract_opt_token(options, "-Znb ", input, sizeof(input)))
    {
        cfg->z_index    = atol(input);
        cfg->z_included = 1;
        printf("Z index is %ld\n", cfg->z_index);
    }

    if (transm_extract_opt_token(options, "-Zerr ", input, sizeof(input)))
    {
        cfg->z_error = atof(input);
        if (cfg->z_error > epsilon)
        {
            cfg->z_included = 1;
        }
        printf("Zerr error is %f rad RMS\n", cfg->z_error);
    }
}

/**
 * transm_accum_tt - Build stellar leakage PSF "psf0" averaged over a uniform tip-tilt disk
 * @cfg:       Throughput configuration
 * @coro_type: Coronagraph architecture code
 *
 * Return: Image ID of "psf0".
 */
static imageID transm_accum_tt(
    const transm_cfg_t *cfg,
    long                coro_type)
{
    const long size  = CORONAGRAPHS_ARRAYSIZE;
    const long size2 = size * size;
    imageID    id0;
    long       cnt = 0;

    create_2Dimage_ID("psf0", size, size, &id0);
    for (int i = -cfg->nb_tt_pts; i < cfg->nb_tt_pts; i++)
    {
        for (int j = -cfg->nb_tt_pts; j < cfg->nb_tt_pts; j++)
        {
            printf("TT point %d %d\n", i, j);
            double x = 1.0 * i * cfg->tt_error / cfg->nb_tt_pts;
            double y = 1.0 * j * cfg->tt_error / cfg->nb_tt_pts;
            if (sqrt(x * x + y * y) < cfg->tt_error)
            {
                coronagraph_simulPSF(x, y, "psf1", coro_type, "");
                imageID id1 = image_ID("psf1", dcimg, dcnimg);
                for (long ii = 0; ii < size2; ii++)
                {
                    dcimg[id0].array.F[ii] += dcimg[id1].array.F[ii];
                }
                delete_image_ID("psf1", DELETE_IMAGE_ERRMODE_WARNING);
                printf(".");
                fflush(stdout);
                cnt++;
            }
        }
    }
    printf("\n");

    for (long ii = 0; ii < size2; ii++)
    {
        dcimg[id0].array.F[ii] /= cnt;
    }
    return id0;
}

/**
 * transm_accum_foc - Build stellar leakage PSF "psf0" averaged over defocus excursion range
 * @cfg:       Throughput configuration
 * @coro_type: Coronagraph architecture code
 *
 * Return: Image ID of "psf0".
 */
static imageID transm_accum_foc(
    const transm_cfg_t *cfg,
    long                coro_type)
{
    const long size  = CORONAGRAPHS_ARRAYSIZE;
    const long size2 = size * size;
    imageID    id0;
    imageID    id_pha;
    long       cnt = 0;

    create_2Dimage_ID("psf0", size, size, &id0);
    create_2Dimage_ID("corphase", size, size, &id_pha);
    for (int i = -cfg->nb_foc_pts; i < cfg->nb_foc_pts; i++)
    {
        double tmp = 1.0 * i * cfg->foc_error / cfg->nb_foc_pts;
        for (long ii = 0; ii < size; ii++)
        {
            for (long jj = 0; jj < size; jj++)
            {
                double x = 1.0 * (ii - size / 2) / cfg->trad_pix;
                double y = 1.0 * (jj - size / 2) / cfg->trad_pix;
                dcimg[id_pha].array.F[jj * size + ii] = (x * x + y * y) * tmp * 2.0 * CORO_PI;
            }
        }
        coronagraph_simulPSF(0.0, 0.0, "psf1", coro_type, "");
        imageID id1 = image_ID("psf1", dcimg, dcnimg);
        for (long ii = 0; ii < size2; ii++)
        {
            dcimg[id0].array.F[ii] += dcimg[id1].array.F[ii];
        }
        delete_image_ID("psf1", DELETE_IMAGE_ERRMODE_WARNING);
        printf(".");
        fflush(stdout);
        cnt++;
    }
    delete_image_ID("corphase", DELETE_IMAGE_ERRMODE_WARNING);
    printf("\n");

    for (long ii = 0; ii < size2; ii++)
    {
        dcimg[id0].array.F[ii] /= cnt;
    }
    return id0;
}

/**
 * transm_auto_piaac_mask - Automatically size PIAAC focal-plane mask radius for target contrast
 * @cfg: Throughput configuration
 */
static void transm_auto_piaac_mask(const transm_cfg_t *cfg)
{
    const long size = CORONAGRAPHS_ARRAYSIZE;
    imageID    id0;

    g_coro.PIAAFPMASKRAD = 0.0;
    if (cfg->tt_included == 1)
    {
        id0 = transm_accum_tt(cfg, 4);
    }
    else if (cfg->foc_included == 1)
    {
        id0 = transm_accum_foc(cfg, 4);
    }
    else
    {
        coronagraph_simulPSF(0.0, 0.0, "psf0", 4, "");
        id0 = image_ID("psf0", dcimg, dcnimg);
    }

    double peak = dcimg[id0].array.F[(size / 2) * size + size / 2];
    long   ii   = 0;
    printf("%g %g\n", dcimg[id0].array.F[(size / 2) * size + ii], peak);
    while (dcimg[id0].array.F[(size / 2) * size + ii] < 2.0 * peak / cfg->contrast_scale)
    {
        ii++;
    }
    g_coro.PIAAFPMASKRAD = CORONAGRAPHS_PIXSCALE * (size / 2 - (ii - 1));
    printf("PIAAC:  MASK RADIUS = %f lambda/d  (ii = %ld)\n", g_coro.PIAAFPMASKRAD, ii);
    delete_image_ID("psf0", DELETE_IMAGE_ERRMODE_WARNING);
}

/**
 * transm_build_psf0 - Compute on-axis stellar leakage PSF "psf0" with configured aberrations
 * @cfg: Throughput configuration
 *
 * Return: Image ID of "psf0".
 */
static imageID transm_build_psf0(const transm_cfg_t *cfg)
{
    const long size = CORONAGRAPHS_ARRAYSIZE;

    if (cfg->tt_included == 1)
    {
        return transm_accum_tt(cfg, cfg->coro_type);
    }
    if (cfg->foc_included == 1)
    {
        return transm_accum_foc(cfg, cfg->coro_type);
    }
    if (cfg->z_included == 1)
    {
        mk_zer("corphase", size, cfg->z_index, cfg->trad_pix + 1);
        imageID id_pha = image_ID("corphase", dcimg, dcnimg);
        for (long ii = 0; ii < size * size; ii++)
        {
            dcimg[id_pha].array.F[ii] *= cfg->z_error;
        }
        coronagraph_simulPSF(0.0, 0.0, "psf0", cfg->coro_type, "");
        save_fl_fits("corphase", "corphase");
        delete_image_ID("corphase", DELETE_IMAGE_ERRMODE_WARNING);
        return image_ID("psf0", dcimg, dcnimg);
    }

    coronagraph_simulPSF(0.0, 0.0, "psf0", cfg->coro_type, "");
    return image_ID("psf0", dcimg, dcnimg);
}

/**
 * transm_score_planet_psf - Evaluate planetary throughput metrics against scaled stellar halo
 * @cfg:        Throughput configuration and sort scratch buffers
 * @id_psf:     Off-axis planetary PSF image ID
 * @out_total2: Output noise-dominated aperture throughput
 * @out_total:  Output pixel-wise contrast-exceeding throughput
 */
static void transm_score_planet_psf(
    transm_cfg_t *cfg,
    imageID       id_psf,
    double       *out_total2,
    double       *out_total)
{
    const long   size  = CORONAGRAPHS_ARRAYSIZE;
    const float *psf0  = dcimg[cfg->id_psf0].array.F;
    const float *psf   = dcimg[id_psf].array.F;
    double       total = 0.0;
    double       peak  = 0.0;

    for (long ii = 0; ii < size; ii++)
    {
        for (long jj = 0; jj < size; jj++)
        {
            long   idx   = jj * size + ii;
            double flux0 = psf0[idx] * cfg->contrast_scale;
            double flux  = psf[idx];
            if (flux > peak)
            {
                peak = flux;
            }
            if (flux > flux0)
            {
                total += flux;
            }
        }
    }

    long kk = 0;
    for (long ii = 0; ii < size; ii++)
    {
        for (long jj = 0; jj < size; jj++)
        {
            long   idx  = jj * size + ii;
            double flux = psf[idx];
            if (flux > 0.000001 * peak)
            {
                cfg->arrayS1[kk] = psf0[idx] * cfg->contrast_scale;
                cfg->arrayS2[kk] = flux;
                cfg->arrayS[kk]  = cfg->arrayS1[kk] / cfg->arrayS2[kk];
                kk++;
            }
        }
    }

    quick_sort3(cfg->arrayS, cfg->arrayS1, cfg->arrayS2, kk);
    double total1 = cfg->arrayS1[0];
    double total2 = cfg->arrayS2[0];
    long   ii     = 1;
    while ((total1 < total2) && (ii < kk))
    {
        total1 += cfg->arrayS1[ii];
        total2 += cfg->arrayS2[ii];
        ii++;
    }

    *out_total2 = total2;
    *out_total  = total;
}

/**
 * transm_is_diagonal_coro - Return 1 if coronagraph type scans along the 45-deg diagonal
 * @coro_type: Coronagraph architecture code
 *
 * Return: 1 for diagonal scan modes, 0 for X-axis scan modes.
 */
static int transm_is_diagonal_coro(long coro_type)
{
    return (coro_type == 3) || (coro_type == 17) || (coro_type == 19) || (coro_type == 20) ||
           (coro_type == 22) || (coro_type == 23) || (coro_type == 30);
}

/**
 * transm_scan_single_axis - Sweep off-axis separation along single optimal azimuth and log IWA50
 * @cfg:   Throughput configuration
 * @fname: Output throughput log filename
 */
static void transm_scan_single_axis(
    transm_cfg_t *cfg,
    const char   *fname)
{
    int    iwa_ok    = 0;
    double sep       = 0.0;
    double sepold    = 0.0;
    double total2old = 0.0;

    while (sep < cfg->sepmax)
    {
        if (transm_is_diagonal_coro(cfg->coro_type))
        {
            double diag = (cfg->coro_type == 30) ? (-sep / sqrt(2.0)) : (sep / sqrt(2.0));
            coronagraph_simulPSF(diag, diag, "psf", cfg->coro_type, "");
        }
        else
        {
            coronagraph_simulPSF(sep, 0.0, "psf", cfg->coro_type, "");
        }

        double total2 = 0.0;
        double total  = 0.0;
        transm_score_planet_psf(cfg, image_ID("psf", dcimg, dcnimg), &total2, &total);

        COROTMP1 += sqrt(total2);
        FILE *fp = fopen(fname, "a");
        if (fp == NULL)
        {
            printf("ERROR: cannot create file \"%s\"\n", fname);
            exit(0);
        }
        fprintf(fp, "%g %g %g\n", sep, total2, total);
        fclose(fp);
        printf("%g %g %g\n", sep, total2, total);
        delete_image_ID("psf", DELETE_IMAGE_ERRMODE_WARNING);

        if ((total2 > 0.5) && (iwa_ok == 0))
        {
            double iwa = sep - (sep - sepold) * (total2 - 0.5) / (total2 - total2old);
            iwa_ok     = 1;
            fp         = fopen(fname, "a");
            fprintf(fp,
                    "#IWA50 %f %f %f %f\n",
                    cfg->tt_error,
                    cfg->logcontrast,
                    g_coro.APLC_FPMASKsize,
                    iwa);
            fclose(fp);
        }

        total2old = total2;
        sepold    = sep;
        sep += cfg->sepstep;
    }
}

/**
 * transm_scan_radaver - Sweep separation with azimuthal averaging over NBanglestep angles
 * @cfg:   Throughput configuration
 * @fname: Output throughput log filename
 */
static void transm_scan_radaver(
    transm_cfg_t *cfg,
    const char   *fname)
{
    double sep = 0.0;
    while (sep < cfg->sepmax)
    {
        double maxtransm = 0.0;
        double total3    = 0.0;
        double total4    = 0.0;

        for (long anglestep = 0; anglestep < cfg->nb_anglestep; anglestep++)
        {
            double angle = cfg->max_angle * anglestep / cfg->nb_anglestep;
            coronagraph_simulPSF(sep * cos(angle), sep * sin(angle), "psf", cfg->coro_type, "");

            double total2 = 0.0;
            double total  = 0.0;
            transm_score_planet_psf(cfg, image_ID("psf", dcimg, dcnimg), &total2, &total);

            total3 += total2;
            total4 += total;
            printf("%g (%ld/%ld) %g %g\n", sep, anglestep, cfg->nb_anglestep, total2, total);
            delete_image_ID("psf", DELETE_IMAGE_ERRMODE_WARNING);
            if (total2 > maxtransm)
            {
                maxtransm = total2;
            }
        }

        FILE *fp = fopen(fname, "a");
        if (fp == NULL)
        {
            printf("ERROR: cannot create file \"%s\"\n", fname);
            exit(0);
        }
        fprintf(fp,
                "%g %g %g %g\n",
                sep,
                total3 / cfg->nb_anglestep,
                total4 / cfg->nb_anglestep,
                maxtransm);
        fclose(fp);
        printf("%g %g %g\n", sep, total3 / cfg->nb_anglestep, total4 / cfg->nb_anglestep);
        sep += cfg->sepstep;
    }
}

/**
 * coronagraph_transm - Measure off-axis planetary throughput curve for target log10 contrast
 * @fname:            Output ASCII profile path
 * @coronagraph_type: Architecture selector index
 * @logcontrast:      Target star-to-planet log10 contrast
 * @options:          Option flags string (-RADAVER, -TT, -FPMSERR, -FOC, -Znb, -Zerr)
 *
 * Return: RETURN_SUCCESS.
 */
errno_t coronagraph_transm(
    const char *fname,
    long        coronagraph_type,
    double      logcontrast,
    const char *options)
{
    const long   size  = CORONAGRAPHS_ARRAYSIZE;
    const size_t size2 = (size_t) size * (size_t) size;

    transm_cfg_t cfg = {
        .coro_type      = coronagraph_type,
        .logcontrast    = logcontrast,
        .contrast_scale = pow(10.0, logcontrast),
        .sepstep        = 0.05,
        .sepmax         = 6.0,
        .trad_pix       = coro_pupil_radius_pix(),
        .tt_error       = 0.001,
        .foc_error      = 0.001,
        .nb_tt_pts      = 20,
        .nb_foc_pts     = 10,
        .nb_anglestep   = (coronagraph_type == 30) ? 64 : 16,
        .max_angle      = (coronagraph_type == 30) ? (2.0 * CORO_PI) : (0.5 * CORO_PI),
    };

    printf("CORONAGRAPH TYPE %ld    OPTION : \"%s\"\n", coronagraph_type, options);
    printf("LOGCONTRAST = %f\n", logcontrast);

    transm_parse_options(options, &cfg);

    if ((g_coro.AUTOPIAACMASK == 1) && (coronagraph_type == 5))
    {
        transm_auto_piaac_mask(&cfg);
    }

    cfg.id_psf0 = transm_build_psf0(&cfg);
    save_fl_fits("psf0", "psf0");

    cfg.arrayS  = (double *) malloc(sizeof(double) * size2);
    cfg.arrayS1 = (double *) malloc(sizeof(double) * size2);
    cfg.arrayS2 = (double *) malloc(sizeof(double) * size2);

    FILE *fp = fopen(fname, "w");
    if (fp == NULL)
    {
        printf("ERROR: cannot create file \"%s\"\n", fname);
        exit(0);
    }
    fclose(fp);

    COROTMP1 = 0.0;
    if (cfg.radaver == 0)
    {
        transm_scan_single_axis(&cfg, fname);
    }
    else
    {
        transm_scan_radaver(&cfg, fname);
    }

    delete_image_ID("psf0", DELETE_IMAGE_ERRMODE_WARNING);
    free(cfg.arrayS);
    free(cfg.arrayS1);
    free(cfg.arrayS2);

    return RETURN_SUCCESS;
}
