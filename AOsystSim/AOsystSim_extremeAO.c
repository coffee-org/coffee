/**
 * @file    AOsystSim_extremeAO.c
 * @brief   Analytical extreme-AO error-budget and raw contrast evaluation across spatial frequency
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "OpticsMaterials/OpticsMaterials.h"

#include "AOsystSim/AOsystSim.h"

/**
 * @brief Auxiliary atmospheric and photometric parameters for the analytical ExAO model.
 */
struct exao_sim_params
{
    double nwfs;
    double nsci;
    double sourcemag_wfs;
    double tobs;
    double WFStlim;
    double IWAld;
};

/**
 * exao_init_conf - Populate target, atmosphere, refractive shear, and flux parameters
 * @cfg: Analytical ExAO configuration struct
 * @p:   Auxiliary simulation parameters
 */
static void exao_init_conf(
    EXAOSIMCONF            *cfg,
    struct exao_sim_params *p)
{
    double lambda_I         = 0.797e-6;
    double zeropt_I         = 4.5825e10;
    double lambda_J         = 1.22e-6;
    double zeropt_J         = 1.9422e10;
    double sourcemag_sci    = 7.3;
    double lambdaBwfs       = 0.2;
    double lambdaBsci       = 0.2;
    double systemEfficiency = 0.2;
    double refract0         = 30.0;
    double drefract;

    cfg->D          = 30.0;
    cfg->lambdawfs  = lambda_I;
    cfg->lambdai    = lambda_J;
    cfg->r0         = 0.15;
    cfg->windspeed  = 10.0;
    p->sourcemag_wfs = 9.0;
    p->tobs          = 3600.0;
    p->WFStlim       = 0.001;
    p->IWAld         = 1.3;

    for (long layer = 0; layer < 20; layer++)
    {
        cfg->CN2layer_h[layer]          = 0.0;
        cfg->CN2layer_coeff[layer]      = 0.0;
        cfg->CN2layer_Chrom_disp[layer] = 0.0;
    }
    cfg->CN2layer_h[0]     = 500.0;
    cfg->CN2layer_coeff[0] = 0.2283;
    cfg->CN2layer_h[1]     = 1000.0;
    cfg->CN2layer_coeff[1] = 0.0883;
    cfg->CN2layer_h[2]     = 2000.0;
    cfg->CN2layer_coeff[2] = 0.0666;
    cfg->CN2layer_h[3]     = 4000.0;
    cfg->CN2layer_coeff[3] = 0.1458;
    cfg->CN2layer_h[4]     = 8000.0;
    cfg->CN2layer_coeff[4] = 0.3350;
    cfg->CN2layer_h[5]     = 16000.0;
    cfg->CN2layer_coeff[5] = 0.1360;

    cfg->lambda0        = 0.55e-6;
    cfg->betapWFS       = sqrt(2.0);
    cfg->betaaWFS       = sqrt(2.0);
    cfg->betapWFSsci    = 2.0;
    cfg->betaaWFSsci    = 2.0;
    cfg->framedelayMult = 2.0;

    p->nwfs = OpticsMaterials_n(OpticsMaterials_code("Air"), cfg->lambdawfs);
    p->nsci = OpticsMaterials_n(OpticsMaterials_code("Air"), cfg->lambdai);
    printf("n = %.12f %.12f\n", p->nwfs, p->nsci);

    drefract = refract0 * (p->nwfs - p->nsci) / (p->nwfs - 1.0);
    printf("differential refraction angle = %.3f arcsec\n", drefract);

    for (long layer = 0; layer < 6; layer++)
    {
        double h        = cfg->CN2layer_h[layer];
        double h0       = 8000.0;
        double ang_rad  = drefract / 3600.0 / 180.0 * M_PI;
        double disp_exp = h0 * ang_rad * (1.0 - exp(-h / h0));

        printf("Simple geom projection   %10.3f km ->  %.6f m\n", h * 0.001, h * ang_rad);
        printf("                                       ->  %.6f m\n", disp_exp);
        cfg->CN2layer_Chrom_disp[layer] = disp_exp;
    }

    cfg->Fwfs = zeropt_I * 1.0e6 * (cfg->lambdawfs * lambdaBwfs) *
                pow(100.0, -0.2 * p->sourcemag_wfs) * lambdaBwfs * systemEfficiency;
    cfg->Fsci = zeropt_J * 1.0e6 * (cfg->lambdai * lambdaBsci) *
                pow(100.0, -0.2 * sourcemag_sci) * lambdaBsci * systemEfficiency;
    cfg->alpha_arcsec = 0.25;
}

/**
 * exao_eval_fresnel - Compute multi-layer Talbot/Fresnel fractions and open-loop contrast
 * @cfg: Analytical ExAO configuration struct at current spatial frequency f
 */
static void exao_eval_fresnel(EXAOSIMCONF *cfg)
{
    double dsa = 0.15;

    cfg->betapWFS = 1.48 / (cfg->f * dsa) * sqrt(1.0 + (dsa * dsa / cfg->r0 / cfg->r0));
    cfg->hf       = 0.22 * cfg->lambda0 /
                    (pow(cfg->f, 11.0 / 6.0) * cfg->D * pow(cfg->r0, 5.0 / 6.0));

    cfg->X    = 0.0;
    cfg->Xwfs = 0.0;
    cfg->dX   = 0.0;
    cfg->dY   = 0.0;

    for (long layer = 0; layer < 20; layer++)
    {
        double h     = cfg->CN2layer_h[layer];
        double coeff = cfg->CN2layer_coeff[layer];
        double f2    = cfg->f * cfg->f;
        double tmpv1 = cos(M_PI * h * f2 * cfg->lambdai);
        double tmpv2 = cos(M_PI * h * f2 * cfg->lambdawfs);
        double tmpv3 = sin(M_PI * h * f2 * cfg->lambdai);
        double tmpv4 = sin(M_PI * h * f2 * cfg->lambdawfs);

        cfg->X += coeff * tmpv1 * tmpv1;
        cfg->Xwfs += coeff * tmpv2 * tmpv2;
        cfg->dX += coeff * (tmpv1 - tmpv2) * (tmpv1 - tmpv2);
        cfg->dY += coeff * (tmpv3 - tmpv4) * (tmpv3 - tmpv4);
    }

    cfg->Y       = 1.0 - cfg->X * cfg->X;
    cfg->Ywfs    = 1.0 - cfg->Xwfs * cfg->Xwfs;
    cfg->CP_UOPD = pow(M_PI * cfg->hf / cfg->lambdai, 2.0) * cfg->Xwfs;
    cfg->CP_UAMP = pow(M_PI * cfg->hf / cfg->lambdai, 2.0) * cfg->Ywfs;
}

/**
 * exao_eval_loop_and_chrom - Evaluate optimal WFS exposure, closed-loop lag, and chromaticity
 * @cfg: Analytical ExAO configuration struct
 * @p:   Auxiliary simulation parameters
 */
static void exao_eval_loop_and_chrom(
    EXAOSIMCONF                  *cfg,
    const struct exao_sim_params *p)
{
    double tmpA;
    double tmpB;
    double n_ratio = (p->nwfs - p->nsci) / (p->nwfs - 1.0);

    tmpA          = 2.0 * M_PI * cfg->hf * cfg->windspeed * cfg->f * cfg->framedelayMult;
    tmpB          = cfg->lambdawfs / M_PI * cfg->betapWFS / sqrt(cfg->Fwfs * M_PI) / cfg->D;
    cfg->twfs_opt = pow(0.5 * tmpB / tmpA * tmpB / tmpA, 1.0 / 3.0);
    cfg->twfs     = (cfg->twfs_opt < p->WFStlim) ? p->WFStlim : cfg->twfs_opt;
    cfg->hfca     = tmpA * cfg->twfs;
    cfg->hfcb     = tmpB / sqrt(cfg->twfs);
    cfg->hfc      = sqrt(cfg->hfca * cfg->hfca + cfg->hfcb * cfg->hfcb);
    cfg->CP_OPHN  = pow(M_PI * cfg->hfcb / cfg->lambdai, 2.0);
    cfg->CP_OTEM  = pow(M_PI * cfg->hfca / cfg->lambdai, 2.0);

    tmpA = 2.0 * M_PI * cfg->hf * cfg->Ywfs * cfg->windspeed * cfg->f * cfg->framedelayMult;
    tmpB = cfg->lambdawfs / M_PI * cfg->betaaWFS / sqrt(cfg->Fwfs * M_PI) / cfg->D;
    cfg->twfs_Aopt = pow(0.5 * tmpB / tmpA * tmpB / tmpA, 1.0 / 3.0);
    cfg->twfs_A    = (cfg->twfs_Aopt < p->WFStlim) ? p->WFStlim : cfg->twfs_Aopt;
    cfg->hfca      = tmpA * cfg->twfs_A;
    cfg->hfcb      = tmpB / sqrt(cfg->twfs_A);
    cfg->CP_APHN   = pow(M_PI * cfg->hfcb / cfg->lambdai, 2.0);
    cfg->CP_ATEM   = pow(M_PI * cfg->hfca / cfg->lambdai, 2.0);

    cfg->CC_OMUL = cfg->CP_UOPD * n_ratio * n_ratio;
    cfg->CC_AMUL = cfg->CP_UAMP * n_ratio * n_ratio;
    cfg->CC_OPRO = cfg->CP_UOPD * cfg->dX;
    cfg->CC_APRO = cfg->CP_UAMP * cfg->dY;

    cfg->CC_ORPA = 0.0;
    cfg->CC_ARPA = 0.0;
    for (long layer = 0; layer < 20; layer++)
    {
        double dpha  = 2.0 * M_PI * cfg->CN2layer_Chrom_disp[layer] * cfg->f;
        double tmpv1 = cfg->hf * sin(dpha / 2.0) * 2.0;
        double c_lay = pow(M_PI * tmpv1 / cfg->lambdai, 2.0) * cfg->CN2layer_coeff[layer];

        cfg->CC_ORPA += c_lay * cfg->Xwfs;
        cfg->CC_ARPA += c_lay * cfg->Ywfs;
    }

    cfg->twfs_opt_amp = cfg->twfs_opt * pow(cfg->X / cfg->Y, 1.0 / 3.0) *
                        pow(cfg->betaaWFS / cfg->betapWFS, 2.0 / 3.0);
}

/**
 * exao_print_summary - Print single-separation error budget diagnostics to stdout
 * @cfg: Analytical ExAO configuration struct
 * @p:   Auxiliary simulation parameters
 */
static void exao_print_summary(
    const EXAOSIMCONF            *cfg,
    const struct exao_sim_params *p)
{
    double crosstime = cfg->D / cfg->windspeed;
    double c_lag     = pow(M_PI * cfg->hfca / cfg->lambdai, 2.0);
    double c_noise   = pow(M_PI * cfg->hfcb / cfg->lambdai, 2.0);

    printf("alpha_ld = %.2lf l/D\n", cfg->alpha_ld);
    printf("f = %f  -> p = %g m\n", cfg->f, 1.0 / cfg->f);
    printf("X = %g\n", cfg->X);
    printf("Y = %g\n", cfg->Y);
    printf("dX = %g\n", cfg->dX);
    printf("dY = %g\n", cfg->dY);
    printf("sourcemag_wfs = %f\n", p->sourcemag_wfs);
    printf("OPTIMAL WFS exposure time = %g sec  -> %.3lf kHz\n",
           cfg->twfs_opt,
           0.001 / cfg->twfs_opt);
    printf("single frequ WF error: %g -> %g  (%g + %g)\n",
           cfg->hf,
           cfg->hfc,
           cfg->hfca,
           cfg->hfcb);
    printf("WFS total flux per frame = %lf ph\n",
           cfg->twfs * cfg->Fwfs * cfg->D * cfg->D / 4.0);
    printf("CP_OPD contrast = %20g\n", cfg->CP_UOPD);
    printf("CP_AMP contrast = %20g\n", cfg->CP_UAMP);
    printf("   WFS  lag = %20g    %20g\n", c_lag, c_lag / sqrt(p->tobs / crosstime) * 2.0);
    printf("   WFS noise= %20g    %20g\n", c_noise, c_noise / sqrt(p->tobs / cfg->twfs) * 2.0);
}

/**
 * exao_write_header - Write 24-column legend header to ContrastTerms.out.txt
 * @fp: Open output file stream
 */
static void exao_write_header(FILE *fp)
{
    fprintf(fp, "# col # 1 : arcsec\n");
    fprintf(fp, "# col # 2 : f [m-1]\n");
    fprintf(fp, "# col # 3 : CP_UOPD  uncorrected OPD\n");
    fprintf(fp, "# col # 4 : CP_UAMP  uncorrected amplitude\n");
    fprintf(fp, "# col # 5 : X\n");
    fprintf(fp, "# col # 6 : Y\n");
    fprintf(fp, "# col # 7 : Xwfs\n");
    fprintf(fp, "# col # 8 : Ywfs\n");
    fprintf(fp, "# col # 9 : dX\n");
    fprintf(fp, "# col #10 : dY\n");
    fprintf(fp, "# col #11 : betaaWFS\n");
    fprintf(fp, "# col #12 : betapWFS\n");
    fprintf(fp, "# col #13 : CR_OPHN (log10)\n");
    fprintf(fp, "# col #14 : CR_OTEM (log10)\n");
    fprintf(fp, "# col #15 : exaosimconf[0].twfs_opt Optimal effective WFS "
                "exposure time for OPD correction\n");
    fprintf(fp, "# col #16 : CR_APHN (log10)\n");
    fprintf(fp, "# col #17 : CR_ATEM (log10)\n");
    fprintf(fp, "# col #18 : exaosimconf[0].twfs_Aopt Optimal effective "
                "WFS exposure time for AMP correction\n");
    fprintf(fp, "# col #19 : CC_OMUL (log10)\n");
    fprintf(fp, "# col #20 : CC_AMUL (log10)\n");
    fprintf(fp, "# col #21 : CC_OPRO (log10)\n");
    fprintf(fp, "# col #22 : CC_APRO (log10)\n");
    fprintf(fp, "# col #23 : CC_ORPA (log10)\n");
    fprintf(fp, "# col #24 : CC_ARPA (log10)\n\n");
}

/**
 * exao_write_row - Append one separation row of contrast terms to ContrastTerms.out.txt
 * @fp:  Open output file stream
 * @cfg: Analytical ExAO configuration struct
 */
static void exao_write_row(
    FILE              *fp,
    const EXAOSIMCONF *cfg)
{
    fprintf(fp, "%10f  %10f    ", cfg->alpha_arcsec, cfg->f);
    fprintf(fp, " %15g %15g   ", log10(cfg->CP_UOPD), log10(cfg->CP_UAMP));
    fprintf(fp, "%10.8f %10.8f %10.8f %10.8f   ", cfg->X, cfg->Y, cfg->Xwfs, cfg->Ywfs);
    fprintf(fp, "%10.8f %10.8f    ", cfg->dX, cfg->dY);
    fprintf(fp, "%10.8f %10.8f    ", cfg->betaaWFS, cfg->betapWFS);
    fprintf(fp, "%15g %15g   ", log10(cfg->CP_OPHN), log10(cfg->CP_OTEM));
    fprintf(fp, "%10.8f   ", cfg->twfs_opt);
    fprintf(fp, "%15g %15g   ", log10(cfg->CP_APHN), log10(cfg->CP_ATEM));
    fprintf(fp, "%10.8g   ", cfg->twfs_Aopt);
    fprintf(fp, "%15g %15g   ", log10(cfg->CC_OMUL), log10(cfg->CC_AMUL));
    fprintf(fp, "%15g %15g   ", log10(cfg->CC_OPRO), log10(cfg->CC_APRO));
    fprintf(fp, "%15g %15g   \n", log10(cfg->CC_ORPA), log10(cfg->CC_ARPA));
}

/**
 * exao_write_gnuplot - Write helper gnuplot script printcmd for ContrastTerms.out.txt
 */
static void exao_write_gnuplot(void)
{
    FILE *fp = fopen("printcmd", "w");

    fprintf(fp, "# gnuplot script\n");
    fprintf(fp, "gnuplot << EOF\n");
    fprintf(fp, "set term pngcairo size 1024,768\n");
    fprintf(fp, "set out \"ContrastTerms.png\"\n");
    fprintf(fp, "set grid\n");
    fprintf(fp, "set size 1,1\n");
    fprintf(fp, "set xlabel \"Angular Separation [arcsec]\"\n");
    fprintf(fp, "set ylabel \"Raw Contrast [10-base log]\"\n");
    fprintf(fp, "plot [0:0.3][-11:0] \"ContrastTerms.out.txt\" u 1:13 lw 3 w l "
                "title \"Residual OPD correction error - WFS photon noise [CR_{OPHN}]\"");
    fprintf(fp, ", \"ContrastTerms.out.txt\" u 1:14 lw 3 w l "
                "title \"Residual OPD correction error - Temporal lag [CR_{OTEM}]\"");
    fprintf(fp, ", \"ContrastTerms.out.txt\" u 1:4 lw 3 w l "
                "title \"Residual Amplitude - Scintillation [CR_{UAMP}]\"");
    fprintf(fp, ", \"ContrastTerms.out.txt\" u 1:16 lw 3 w l "
                "title \"Residual Amplitude correction error - WFS photon noise [CR_{APHN}]\"");
    fprintf(fp, ", \"ContrastTerms.out.txt\" u 1:17 lw 3 w l "
                "title \"Residual Amplitude correction error - Temporal lag [CR_{ATEM}]\"");
    fprintf(fp, ", \"ContrastTerms.out.txt\" u 1:19 lw 3 w l "
                "title \"Chromatic OPD - Multiplicative refractivity [CC_{OMUL}]\"");
    fprintf(fp, ", \"ContrastTerms.out.txt\" u 1:20 lw 3 w l "
                "title \"Chromatic Amplitude - Multiplicative refractivity [CC_{AMUL}]\"");
    fprintf(fp, ", \"ContrastTerms.out.txt\" u 1:21 lw 3 w l "
                "title \"Chromatic OPD - Fresnel propagation [CC_{OPRO}]\"");
    fprintf(fp, ", \"ContrastTerms.out.txt\" u 1:22 lw 3 w l "
                "title \"Chromatic Amplitude - Fresnel propagation [CC_{APRO}]\"");
    fprintf(fp, ", \"ContrastTerms.out.txt\" u 1:23 lw 3 w l "
                "title \"Chromatic OPD - Refractive Light Path [CC_{OPRA}]\"");
    fprintf(fp, ", \"ContrastTerms.out.txt\" u 1:24 lw 3 w l "
                "title \"Chromatic Amplitude - Refractive Light Path [CC_{APRA}]\"\n");
    fprintf(fp, "exit\nEOF\n\n");
    fclose(fp);
}

/**
 * AOsystSim_extremeAO_contrast_sim - Evaluate analytical extreme-AO raw contrast budget
 *
 * Return: RETURN_SUCCESS on completion.
 */
errno_t AOsystSim_extremeAO_contrast_sim(void)
{
    EXAOSIMCONF           *cfg;
    struct exao_sim_params params = { 0 };
    FILE                  *fp;
    int                    initprint = 0;

    cfg = (EXAOSIMCONF *) malloc(sizeof(EXAOSIMCONF));
    if (cfg == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    exao_init_conf(cfg, &params);

    printf("Writing output file ... \n");
    fflush(stdout);
    fp = fopen("ContrastTerms.out.txt", "w");

    for (cfg->alpha_arcsec = 0.010; cfg->alpha_arcsec < 2.0; cfg->alpha_arcsec += 0.001)
    {
        int OK = 0;

        printf("\nSeparation = %f\n", cfg->alpha_arcsec);
        fflush(stdout);

        while (OK == 0)
        {
            cfg->alpha    = cfg->alpha_arcsec / 3600.0 / 180.0 * M_PI;
            cfg->alpha_ld = cfg->alpha / (cfg->lambdai / cfg->D);
            cfg->f        = cfg->alpha / cfg->lambdai;
            if (cfg->alpha_ld > params.IWAld)
            {
                OK = 1;
            }
            else
            {
                cfg->alpha_arcsec += 0.001;
            }
        }

        exao_eval_fresnel(cfg);
        exao_eval_loop_and_chrom(cfg, &params);
        exao_print_summary(cfg, &params);

        if (initprint == 0)
        {
            exao_write_header(fp);
            initprint = 1;
        }
        exao_write_row(fp, cfg);
    }

    fclose(fp);
    exao_write_gnuplot();
    free(cfg);

    return RETURN_SUCCESS;
}
