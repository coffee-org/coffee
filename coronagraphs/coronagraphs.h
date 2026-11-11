/**
 * @file    coronagraphs.h
 * @brief   Public API for prolate apodizers, APLC tables, PIAA remapping, and coronagraph PSFs
 */

#ifndef _CORONAGRAPHS_H
#define _CORONAGRAPHS_H

#define CORONAGRAPHS_ARRAYSIZE 4096

/* Compute 2D prolate spheroidal pupil apodizer for focal mask radius in pixels. */
errno_t coronagraph_make_2Dprolate(
    double      masksizepix,
    double      beamradpix,
    double      centralObs,
    const char *outname,
    long        size,
    const char *pupmask_name,
    double     *transmval);

/* Compute 2D prolate spheroidal pupil apodizer for focal mask radius in lambda/D. */
errno_t coronagraph_make_2Dprolateld(
    double      masksizeld,
    double      beamradpix,
    double      centralObs,
    const char *outname,
    long        size,
    const char *pupmask_name);

/* Iteratively refine an existing 2D prolate apodizer using zoomed DFTs. */
errno_t coronagraph_update_2Dprolate(
    double masksizeld,
    double beamradpix,
    double centralObs,
    double zfactor);

/* Compute circularly symmetric 1D/2D prolate apodizer profile. */
errno_t coronagraph_make_2Dprolate_CS(
    double      masksize,
    double      centralObs,
    const char *outname);

/* Compile precomputed APLC prolate apodizer lookup tables. */
errno_t coronagraph_APLCapo_compile(void);

/* Initialize radial PIAA apodization and ray-remapping tables. */
errno_t coronagraph_init_PIAA(long *outNBpoints);

/* Generate nominal Subaru telescope pupil geometry mask. */
errno_t coronagraphs_make_SUBARU_pupil(void);

/* Remap input pupil amplitude and OPD maps across radial PIAA optics. */
errno_t coronagraphs_PIAA_apodize_beam(
    const char *ampl1,
    const char *opd1,
    const char *ampl2,
    const char *opd2);

/* Simulate Achromatic Interfero-Coronagraph (AIC) PSF for source at (xld, yld) [l/D]. */
errno_t coronagraph_simul_AIC(
    double      xld,
    double      yld,
    const char *psfname);

/* Simulate Four-Quadrant Phase Mask (4QPM) coronagraph PSF at (xld, yld) [l/D]. */
errno_t coronagraph_simul_4QPM(
    double      xld,
    double      yld,
    const char *psfname);

/* Simulate 8th-order Band-Limited (BL8) Lyot coronagraph PSF at (xld, yld) [l/D]. */
errno_t coronagraph_simul_BL8(
    double      xld,
    double      yld,
    const char *psfname);

/* Simulate Roddier & Roddier Phase Mask (RRPM) coronagraph PSF at (xld, yld) [l/D]. */
errno_t coronagraph_simul_RRPM(
    double      xld,
    double      yld,
    const char *psfname);

/* Simulate Optical Vortex Coronagraph (OVC) PSF at (xld, yld) [l/D]. */
errno_t coronagraph_simul_OVC(
    double      xld,
    double      yld,
    const char *psfname);

/* Simulate Conventional Pupil Apodization (CPA) coronagraph PSF at (xld, yld) [l/D]. */
errno_t coronagraph_simul_CPA(
    double      xld,
    double      yld,
    const char *psfname);

/* Simulate lossless PIAA coronagraph PSF at (xld, yld) [l/D]. */
errno_t coronagraph_simul_PIAA(
    double      xld,
    double      yld,
    const char *psfname);

/* Simulate PIAA Complex Mask Coronagraph (PIAAC) PSF at (xld, yld) [l/D]. */
errno_t coronagraph_simul_PIAAC(
    double      xld,
    double      yld,
    const char *psfname);

/* Simulate hybrid AIC + PIAAC coronagraph PSF at (xld, yld) [l/D]. */
errno_t coronagraph_simul_AIC_PIAAC(
    double      xld,
    double      yld,
    const char *psfname);

/* Simulate multi-step Apodized Pupil Lyot Coronagraph (APLC) PSF at (xld, yld) [l/D]. */
errno_t coronagraph_simul_MULTISTEP_APLC(
    double      xld,
    double      yld,
    const char *psfname);

/* Dispatch coronagraph PSF simulation by coronagraph_type code at (xld, yld) [l/D]. */
errno_t coronagraph_simulPSF(
    double      xld,
    double      yld,
    const char *psfname,
    long        coronagraph_type,
    const char *options);

/* Measure off-axis planetary throughput curve for a given target log10 contrast. */
errno_t coronagraph_transm(
    const char *fname,
    long        coronagraph_type,
    double      logcontrast,
    const char *options);

/* Execute user-configured coronagraph evaluation hook. */
errno_t coronagraph_userfunc(void);

/* Compute fundamental coronagraphic inner working angle / contrast limit coefficients. */
errno_t coronagraph_compute_limitcoeff(void);

/* Scan idealized PIAACMC contrast performance as a function of input central obscuration. */
errno_t CORONAGRAPHS_scanPIAACMC_centObs_perf(double obs0input);

#endif /* _CORONAGRAPHS_H */
