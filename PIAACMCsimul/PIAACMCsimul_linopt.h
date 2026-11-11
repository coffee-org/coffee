/**
 * @file    PIAACMCsimul_linopt.h
 * @brief   Regularized Jacobian perturbation and multi-alpha SVD line-search optimizer
 */

#ifndef PIAACMCSIMUL_LINOPT_H
#define PIAACMCSIMUL_LINOPT_H

#include "COREMOD_memory/COREMOD_memory.h"

/**
 * @brief Shared per-iteration context for the regularized SVD line-search optimizer.
 */
typedef struct
{
    imageID IDstatus;   /**< Shared-memory status telemetry image ID */
    long    iter;       /**< Current outer linear optimization iteration */
    long    size1Dvec;  /**< Total 1D residual + regularization vector length */
    long    size1Dvec0; /**< Dark-hole complex electric-field vector length */
    double  valref;     /**< Initial reference regularized objective value */
    double  bestval;    /**< Best regularized objective value in current scan */
    double  bestgain;   /**< Line-search step gain achieving `bestval` */
} PIAACMC_linopt_ctx;

/* Evaluate on-axis PSF contrast plus active PIAA/FPM regularization penalties. */
errno_t PIAACMCsimul_linopt_eval_objective(
    double *val_total,
    double *val_contrast,
    double *val_reg_piaa,
    double *val_reg_fpm);

/* Solve 3-cutoff SVD directions and scan step gains across alphareg in [0, 1]. */
errno_t PIAACMCsimul_linopt_svd_linescan(
    PIAACMC_linopt_ctx *ctx,
    imageID            *IDoptvec_out,
    long               *NBlinoptgain_out);

/* Execute the full regularized linear perturbation loop when piaacmcparams.LINOPT == 1. */
errno_t PIAACMCsimul_run_linopt(
    imageID     IDstatus,
    const char *stopfile);

#endif /* PIAACMCSIMUL_LINOPT_H */
