/**
 * @file    AOsystSim.h
 * @brief   Public structures and API for Adaptive Optics system, WFS, DM, and LOWFS simulation
 */

#ifndef _AOSYSTSIM_H
#define _AOSYSTSIM_H

/**
 * @brief Analytical extreme-AO error-budget and contrast model configuration.
 */
typedef struct
{
    double alpha;        /**< On-sky angular separation [rad] */
    double alpha_ld;     /**< Angular separation [lambda/D] */
    double alpha_arcsec; /**< Angular separation [arcsec] */
    double lambda0;      /**< Seeing reference wavelength [m] */
    double lambdai;      /**< Science imaging wavelength [m] */
    double lambdawfs;    /**< Wavefront sensor wavelength [m] */

    double CN2layer_h[20];          /**< Turbulence layer altitudes [m] */
    double CN2layer_coeff[20];      /**< Normalized Cn^2 weights (sum = 1) */
    double CN2layer_Chrom_disp[20]; /**< Atmospheric refractive chromatic shear [m] */

    double D;           /**< Telescope aperture diameter [m] */
    double r0;          /**< Fried parameter at lambda0 [m] */
    double windspeed;   /**< Effective wind speed [m/s] */
    double betapWFS;    /**< Phase WFS photon-noise sensitivity factor */
    double betaaWFS;    /**< Amplitude WFS photon-noise sensitivity factor */
    double betapWFSsci; /**< Science-band phase WFS sensitivity factor */
    double betaaWFSsci; /**< Science-band amplitude WFS sensitivity factor */
    double Fwfs;        /**< Guide star photon flux in WFS band [ph/s/m^2] */
    double Fsci;        /**< Guide star photon flux in science band [ph/s/m^2] */

    double framedelayMult; /**< Loop latency in WFS frame units */

    double f;  /**< Spatial frequency [m^-1] */
    double hf; /**< Open-loop Fourier mode RMS amplitude [m] */
    double X;
    double Xwfs;
    double Y;
    double Ywfs;
    double dX;
    double dY;

    double CP_UOPD; /**< Uncorrected OPD contrast contribution */
    double CP_OPHN; /**< WFS photon-noise OPD contrast contribution */
    double CP_OTEM; /**< Servo-lag OPD contrast contribution */

    double CP_UAMP; /**< Uncorrected scintillation contrast contribution */
    double CP_APHN; /**< WFS photon-noise amplitude contrast contribution */
    double CP_ATEM; /**< Servo-lag amplitude contrast contribution */

    double CC_OMUL;
    double CC_OPRO;
    double CC_ORPA;

    double CC_AMUL;
    double CC_APRO;
    double CC_ARPA;

    double twfs;     /**< WFS integration time for OPD control [s] */
    double twfs_opt; /**< Optimal OPD WFS integration time [s] */

    double twfs_A;    /**< WFS integration time for amplitude control [s] */
    double twfs_Aopt; /**< Optimal amplitude WFS integration time [s] */

    double twfssci;
    double hfc;          /**< Residual closed-loop Fourier mode amplitude [m] */
    double twfssci_opt;  /**< Optimal science-band phase WFS exposure time [s] */
    double twfs_opt_amp; /**< Optimal science-band amplitude WFS exposure time [s] */
    double hfca;         /**< Servo time-lag error term [m] */
    double hfcb;         /**< Photon-noise error term [m] */
    double C2;           /**< Combined servo-lag + photon-noise phase contrast */
    double C2_wfs;       /**< WFS-band phase contrast */
    double C3;           /**< Residual scintillation contrast */
    double C4;           /**< Chromatic OPD contrast */
    double C5;           /**< Chromatic amplitude contrast */
    double C6;           /**< Dispersive refractive-index contrast */
    double Csum;         /**< Total raw PSF contrast */
    double Csum_detection; /**< 5-sigma post-detection contrast limit */

    double RIC_hfca; /**< Refractive index chromaticity: servo lag term */
    double RIC_hfcb; /**< Refractive index chromaticity: photon noise term */
    double RIC_hfc;
    double TL_hfca; /**< Fast near-IR loop: servo lag term */
    double TL_hfcb; /**< Fast near-IR loop: photon noise term */
    double TL_hfc;
    double C7;      /**< Post-near-IR servo-lag contrast */
    double C8;      /**< Post-near-IR photon-noise contrast */
    double C9;
    double C10;
    double C11; /**< Post-near-IR chromatic refractive-index contrast */
    double Csum2;
    double Csum2ave;
} EXAOSIMCONF;

/* Run real-time shared-memory AO loop simulator (DM + WF stream -> WFS images). */
int AOsystSim_run(
    int  syncmode,
    long DMindex,
    long delayus);

/* --- 1. Analytical High-Contrast PSF & Spatial-Frequency Filtering --- */

/* Apply open-loop spatial-frequency AO transfer filter (lag, gain, WFS noise). */
int AOsystSim_simpleAOfilter(
    const char *IDin_name,
    const char *IDout_name);

/* Evaluate analytical extreme-AO raw contrast budget across spatial frequencies. */
errno_t AOsystSim_extremeAO_contrast_sim(void);

/* --- 2. Telescope Pupil Geometry & Fitting --- */

/* Synthesize parametric obstructed/spidered telescope pupil on DM actuator grid. */
long AOsystSim_mkTelPupDM(
    const char *ID_name,
    long        msize,
    double      xc,
    double      xy,
    double      rin,
    double      rout,
    double      pupPA,
    double      spiderPA,
    double      spideroffset,
    double      spiderthick,
    double      stretchx);

/* Fit pupil registration parameters (offset, scale, rotation, spider geometry) to DM map. */
long AOsystSim_fitTelPup(
    const char *ID_name,
    const char *IDtelpup_name);

/* --- 3. Atmospheric Wavefront Stream Generation --- */

/* Write template configuration file for multi-layer atmospheric WF simulation. */
int AOsystSim_mkWF_mkCONF(const char *fname);

/* Stream multi-layer Fresnel-propagated atmospheric wavefronts from configuration file. */
int AOsystSim_mkWF(const char *CONF_FNAME);

/* --- 4. Pyramid Wavefront Sensor Simulation --- */

/* Compute modulated Pyramid WFS intensity image from complex pupil wavefront cube. */
int AOsystSim_WFSsim_Pyramid(
    const char *inWFc_name,
    const char *outWFSim_name,
    double      modampl,
    long        modnbpts);

/* Execute configured WFS simulation instance and publish output frame. */
int AOsystSim_runWFS(
    long        index,
    const char *IDout_name);

/* Write template configuration file for real-time Pyramid WFS daemon. */
int AOsystSim_PyrWFS_mkCONF(const char *fname);

/* Run real-time shared-memory Pyramid WFS simulation loop. */
int AOsystSim_PyrWFS(const char *CONF_FNAME);

/* --- 5. Deformable Mirror Surface Simulation --- */

/* Convolve actuator command map with 3D influence-function cube into 2D OPD surface. */
int AOsystSim_DMshape(
    const char *IDdmctrl_name,
    const char *IDdmifc_name,
    const char *IDdm_name);

/* Write template configuration file for deformable mirror stream process. */
int AOsystSim_DM_mkCONF(const char *fname);

/* Run real-time shared-memory deformable mirror OPD surface daemon. */
int AOsystSim_DM(const char *CONF_FNAME);

/* --- 6. Coronagraph & Lyot Low-Order WFS (LOWFS) --- */

/* Write template configuration file for coronagraph + LOWFS stream process. */
int AOsystSim_coroLOWFS_mkCONF(const char *fname);

/* Run real-time shared-memory coronagraph + LOWFS simulation loop. */
int AOsystSim_coroLOWFS(const char *CONF_FNAME);

/* --- 7. Focal-Plane Wavefront Sensing & Speckle Nulling (FPWFS) --- */

/* Simulate focal-plane electric-field estimation images under DM phase probes. */
long AOsystSim_FPWFS_imsimul(
    double probeamp,
    double sepx,
    double sepy,
    double contrast,
    double wferramp,
    double totFlux,
    double DMgainErr,
    double RON,
    double CnoiseFloor);

/* Generate pair of DM speckle-nulling probe maps over target dark-hole CPA band. */
int AOsystSim_FPWFS_mkprobes(
    const char *IDprobeA_name,
    const char *IDprobeB_name,
    long        dmxsize,
    long        dmysize,
    double      CPAmax,
    double      CPArmin,
    double      CPArmax,
    double      RMSampl,
    long        modegeom);

/* Run Monte Carlo sensitivity analysis of focal-plane wavefront sensing accuracy. */
int AOsystSim_FPWFS_sensitivityAnalysis(
    int mapmode,
    int mode,
    int optmode,
    int NBprobes);

#endif /* _AOSYSTSIM_H */
