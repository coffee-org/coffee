/**
 * @file    PIAACMCsimul.h
 * @brief   Core state structures and global context for PIAACMC and APLCMC coronagraph design
 */

#ifndef PIAACMCSIMUL_H
#define PIAACMCSIMUL_H

#define ApoFitCosFact 1.0

#include <stdint.h>

#include "milk_runnable.h"
#include "OptSystProp/OptSystProp.h"

/**
 * @brief Runtime execution flags, working buffers, and optimizer state for PIAACMCsimul.
 */
typedef struct
{
    char piaacmcconfdir[STRINGMAXLEN_DIRNAME]; /**< Active configuration directory path */
    int  optsystinit;                          /**< Non-zero once OPTSYST has been initialized */

    int FORCE_CREATE_Cmodes;
    int CREATE_Cmodes;
    int FORCE_CREATE_Fmodes;
    int CREATE_Fmodes;

    int FORCE_CREATE_fpmzmap;
    int CREATE_fpmzmap;
    int FORCE_CREATE_fpmzt;
    int CREATE_fpmzt;

    int FORCE_CREATE_fpmza;
    int CREATE_fpmza;

    int FORCE_MAKE_PIAA0shape;
    int MAKE_PIAA0shape;
    int FORCE_MAKE_PIAA1shape;
    int MAKE_PIAA1shape;

    int focmMode;           /**< Zone index for single-zone impulse response (-1 = full mask) */
    int PIAACMC_FPMsectors; /**< 0: concentric rings, 1: hex zones, 2: 60-deg sectors */

    double FPMSCALEFACTOR; /**< Mask array undersizing factor to prevent DFT edge clipping */

    double LAMBDASTART; /**< Short-wavelength bandpass edge [m] */
    double LAMBDAEND;   /**< Long-wavelength bandpass edge [m] */

    float PIAACMC_MASKRADLD; /**< Physical FPM outer radius override [lambda/D] */

    /* Pre-allocated minimization working buffers */
    double *fpmresp_array;
    double *zonez_array;
    double *zonez0_array;
    double *zonez1_array;
    double *zonezbest_array;
    double *dphadz_array;
    double *outtmp_array;
    long    LOOPCNT;

    long   vsize;       /**< Dark-hole evaluation vector length (2 * points) */
    double CnormFactor; /**< Unocculted PSF peak normalization factor */

    int computePSF_FAST_FPMresp;
    int computePSF_ResolvedTarget;      /**< Stellar radius code: 10^(-0.1 * code) [lambda/D] */
    int computePSF_ResolvedTarget_mode; /**< Extended source quadrature: 0 = 3-pt, 1 = 6-pt */
    int PIAACMC_FPM_FASTDERIVATIVES;

    double SCORINGTOTAL;
    double MODampl;
    int    SCORINGMASKTYPE;
    int    PIAACMC_save;
    int    PIAACMC_fpmtype; /**< 0: idealized monochromatic FPM, 1: physical multi-zone FPM */

    long PIAACMC_FPMresp_mp;     /**< Worker process count in polychromatic FPM calibration */
    long PIAACMC_FPMresp_thread; /**< Zero-based worker index for Mode 11 slice calibration */

    int WRITE_OK;

    double PIAACMCSIMUL_VAL;
    double PIAACMCSIMUL_VAL0;
    double PIAACMCSIMUL_VALREF;

    /* Linear perturbation optimizer (Modes 40 / 41) */
    int     LINOPT;
    long    linopt_number_param;
    int     linopt_paramtype[10000];
    float  *linopt_paramvalf[10000];
    double *linopt_paramval[10000];
    double  linopt_paramrefval[10000];

    double linopt_paramdelta[10000];
    double linopt_paramdeltaval[10000];
    double linopt_parammaxstep[10000];
    double linopt_parammin[10000];
    double linopt_parammax[10000];

    int   linopt_REGPIAASHAPES;
    float linopt_piaa0C_regcoeff;
    float linopt_piaa1C_regcoeff;
    float linopt_piaa0C_regcoeff_alpha;
    float linopt_piaa1C_regcoeff_alpha;

    float linopt_piaa0F_regcoeff;
    float linopt_piaa1F_regcoeff;
    float linopt_piaa0F_regcoeff_alpha;
    float linopt_piaa1F_regcoeff_alpha;

    int  linopt_REGFPMSAG;
    long linopt_NBiter;

    char fnamedescr[STRINGMAXLEN_FILENAME];      /**< Full FPM solution filename tag */
    char fnamedescr_conf[STRINGMAXLEN_FILENAME]; /**< FPM geometry/bandpass filename tag */
} PIAACMCSIMUL_PARAMS;

#define STRINGMAXLEN_PIAACMCSIMUL_MATERIALNAME 10

/**
 * @brief Physical optical prescription of a PIAACMC or APLCMC coronagraph design.
 *
 * Focal plane mask radial layout (from center outward):
 * 1. **Central cone** (optional, `NBringCentCone` inner rings): linear sag profile from
 *    peak `fpmCentConeZ` [m] at `r = 0` to `0` at `fpmCentConeRad` [m].
 * 2. **Multi-zone rings / sectors**: independent sag per zone bounded in
 *    `[fpmminsag, fpmmaxsag]` [m].
 * 3. **Outer cone** (optional, out to `fpmOuterConeRad` [m]): linear sag ramp from `0`
 *    at the outer ring boundary to `fpmOuterConeZ` [m] beyond `fpmOuterConeRad`.
 */
typedef struct
{
    /* Seed radial PIAA geometry */
    double centObs0; /**< Entrance pupil central obscuration ratio [beamrad] */
    double centObs1; /**< Remapped output pupil central obscuration ratio [beamrad] */
    double r0lim;    /**< Extrapolated outer radius of PIAA element 0 [beamrad] */
    double r1lim;    /**< Extrapolated outer radius of PIAA element 1 [beamrad] */
    long   NBradpts; /**< Radial grid sample count in 1D PIAA sag table */

    /* Spectral bandpass */
    int    nblambda;          /**< Number of discrete spectral channels */
    double lambda;            /**< Central design wavelength [m] */
    double lambdaB;           /**< Fractional spectral bandwidth [%] */
    double lambdaarray[2000]; /**< Discrete channel wavelengths [m] */

    /* Collimated pupil & PIAA conjugation geometry */
    float beamrad;         /**< Nominal collimated beam radius [m] */
    long  size;            /**< 2D numerical wavefront grid dimension [pix] */
    float pixscale;        /**< Collimated pupil pixel scale [m/pix] */
    int   PIAAmode;        /**< 0: APLCMC (no PIAA optics), 1: PIAACMC */
    float PIAA0pos;        /**< Axial conjugation of PIAA element 0 from pupil [m] */
    float PIAAsep;         /**< Axial separation between PIAA elements 0 and 1 [m] */
    int   prePIAA0mask;    /**< Enable opaque mask upstream of PIAA element 0 (0/1) */
    float prePIAA0maskpos; /**< Axial conjugation of pre-PIAA0 mask [m] */
    int   postPIAA0mask;   /**< Enable opaque mask downstream of PIAA element 0 (0/1) */
    float postPIAA0maskpos; /**< Axial conjugation of post-PIAA0 mask [m] */
    float PIAAcoeff;       /**< Fraction of pupil apodization carried by PIAA (0..1) */
    int   invPIAAmode;     /**< 0: none, 1: inverse PIAA after Lyot, 2: before Lyot */
    float LyotZmin;        /**< Minimum axial search bound for Lyot stops [m] */
    float LyotZmax;        /**< Maximum axial search bound for Lyot stops [m] */
    float pupoutmaskrad;   /**< Output pupil stop radius [beamrad] */

    /* Deformable mirrors (up to 10) */
    int    nbDM;      /**< Active deformable mirror count */
    double DMpos[10]; /**< DM axial conjugations in collimated space [m] */
    long   ID_DM[10]; /**< DM surface OPD image IDs */

    /* Multi-plane Lyot stops (up to 10) */
    long   NBLyotStop;       /**< Active Lyot stop count */
    long   IDLyotStop[10];   /**< Lyot stop transmission mask image IDs */
    double LyotStop_zpos[10]; /**< Lyot stop axial conjugations [m] */

    /* PIAA surface modal basis (radial cosine + 2D Fourier) */
    char  PIAAmaterial_name[10];
    int   PIAAmaterial_code;
    long  CmodesID;        /**< 2D radial cosine basis cube image ID */
    long  Cmsize;          /**< Cropped modal grid dimension [pix] */
    long  NBCmodes;        /**< Active radial cosine mode count */
    long  piaaNBCmodesmax; /**< Maximum radial cosine mode count */
    long  FmodesID;        /**< 2D Fourier basis cube image ID */
    long  Fmsize;          /**< Cropped Fourier modal grid dimension [pix] */
    long  NBFmodes;        /**< Active 2D Fourier mode count */
    float piaaCPAmax;      /**< Maximum 2D Fourier spatial frequency [CPA] */

    long piaa0CmodesID; /**< PIAA element 0 cosine coefficient vector ID */
    long piaa0FmodesID; /**< PIAA element 0 Fourier coefficient vector ID */
    long piaa1CmodesID; /**< PIAA element 1 cosine coefficient vector ID */
    long piaa1FmodesID; /**< PIAA element 1 Fourier coefficient vector ID */

    float peakPSF; /**< Reference unocculted PSF peak intensity */

    /* Focal plane mask prescription */
    double fpmaskradld;     /**< Idealized monochromatic FPM radius [lambda/D] */
    long   focmNBzone;      /**< Active independent FPM zone count */
    double Fratio;          /**< Beam focal ratio at focal plane mask */
    long   zonezID;         /**< Zone sag double image ID (`fpmzt` / `fpm_zonez.fits`) [m] */
    double fpmaskamptransm; /**< Idealized FPM amplitude transmission */
    long   zoneaID;         /**< Zone amplitude image ID (`fpmza` / `fpm_zonea.fits`) */
    double fpzfactor;       /**< Focal plane DFT oversampling factor */
    double fpmRad;          /**< Physical FPM outer radius [m] */
    long   NBrings;         /**< Concentric hexagonal ring count */
    double fpmminsag;       /**< Minimum allowed FPM zone sag [m] */
    double fpmmaxsag;       /**< Maximum allowed FPM zone sag [m] */
    double fpmsagreg_coeff; /**< FPM sag regularization weight */
    double fpmsagreg_alpha; /**< FPM sag regularization exponent */
    long   NBringCentCone;  /**< Number of inner rings spanned by central cone */
    double fpmCentConeRad;  /**< Outer radius of central cone [m] */
    double fpmCentConeZ;    /**< Peak sag at center of central cone [m] */
    double fpmOuterConeZ;   /**< Sag amplitude at outer edge of outer cone [m] */
    double fpmOuterConeRadld; /**< Outer cone outer radius [lambda/D] */
    double fpmOuterConeRad; /**< Outer cone outer radius [m] */
    long   fpmarraysize;    /**< High-resolution 2D FPM complex grid size [pix] */
    char   fpmmaterial_name[STRINGMAXLEN_PIAACMCSIMUL_MATERIALNAME];
    int    fpmmaterial_code;
} OPTPIAACMCDESIGN;

extern PIAACMCSIMUL_PARAMS piaacmcparams;
extern OPTPIAACMCDESIGN    piaacmcopticaldesign;
extern OPTSYST             piaacmcopticalsystem;

/* Free dynamically allocated PIAACMC optimization working arrays. */
void PIAACMCsimul_free(void);

#endif /* PIAACMCSIMUL_H */
