/**
 * @file    OptSystProp.h
 * @brief   Public structures and API for multi-wavelength Fresnel optical propagation
 */

#ifndef _OPTSYSTPROP_H
#define _OPTSYSTPROP_H

#include "milkDebugTools.h"

/**
 * @brief Square-grid deformable mirror simulation descriptor.
 */
typedef struct
{
    long   NBact1D;    /**< Linear actuator count across pupil */
    double pitch;      /**< Inter-actuator pitch [m] */
    double maxstroke;  /**< Actuator stroke half-range (-maxstroke..+maxstroke) [m] */
    long   dispID;     /**< Actuator command displacement matrix image ID */
    long   IF_ID;      /**< Actuator influence function map image ID */
    double IFpixscale; /**< Influence function pixel scale [m/pix] */
    double IFsize;     /**< Influence function map linear size [pix] */
    double pixscale;   /**< Output surface map pixel scale [m/pix] */
    long   dispmapID;  /**< Interpolated 2D surface displacement map image ID */
} DM_SIM;

/**
 * @brief Reflective aspheric mirror surface (sampled on nominal collimated grid).
 */
typedef struct
{
    long surfID; /**< Surface Z sag image ID [m] */
} ASPHSURFM;

/**
 * @brief Refractive aspheric optical surface (sampled on nominal collimated grid).
 */
typedef struct
{
    long   mat0;        /**< Incident medium material code (see OpticsMaterials) */
    long   mat1;        /**< Transmitted medium material code (see OpticsMaterials) */
    int    init;        /**< Non-zero once chromatic OPD scale factors are initialized */
    double ncoeff[100]; /**< Per-wavelength sag-to-OPD conversion factor (max 100 ch) */
    long   surfID;      /**< Surface Z sag image ID [m] */
} ASPHSURFR;

/**
 * @brief Focal plane mask element applied via zoomed pupil-to-focal DFT.
 */
typedef struct
{
    long   fpmID;   /**< Complex focal plane mask cube image ID (CA or 1 - CA) */
    double zfactor; /**< Focal-plane DFT oversampling zoom factor */
    int    mode;    /**< Babinet mode: 1 if fpmID stores (1 - CA), 0 if direct CA */
} FOCMASK;

#define OPTSYST_ELEM_MAXNB            100
#define STRINGMAXLEN_OPTSYST_ELEMNAME 100

/**
 * @brief Sequential multi-element optical system model in equivalent collimated space.
 *
 * Elements are ordered sequentially along the optical train with axial conjugation
 * @ref elemZpos referenced to the entrance pupil (`z = 0`). Supported @ref elemtype
 * codes:
 * - `1`: Static opaque/apodizing amplitude mask (`elemarrayindex` = image ID)
 * - `3`: Reflective aspheric mirror (`ASPHSURFMarray[elemarrayindex]`)
 * - `4`: Refractive aspheric lens/plate (`ASPHSURFRarray[elemarrayindex]`)
 * - `5`: Focal plane mask via zoomed 2D DFT (`FOCMASKarray[elemarrayindex]`)
 */
typedef struct
{
    int    nblambda;          /**< Number of discrete spectral channels */
    double lambdaarray[2000]; /**< Wavelength grid [m] */

    double beamrad;    /**< Entrance beam radius in collimated space [m] */
    double pixscale;   /**< Transverse pixel scale in collimated space [m/pix] */
    long   size;       /**< 2D square wavefront grid dimension [pix] */
    long   DFTgridpad; /**< Focal-plane DFT sub-sampling stride (0 = full resolution) */

    long NBelem; /**< Active optical element count */
    char name[OPTSYST_ELEM_MAXNB][STRINGMAXLEN_OPTSYST_ELEMNAME];

    long      NB_asphsurfm; /**< Registered reflective aspheric surface count */
    ASPHSURFM ASPHSURFMarray[OPTSYST_ELEM_MAXNB];

    long      NB_asphsurfr; /**< Registered refractive aspheric surface count */
    ASPHSURFR ASPHSURFRarray[OPTSYST_ELEM_MAXNB];

    long    NB_focmask; /**< Registered focal plane mask count */
    FOCMASK FOCMASKarray[OPTSYST_ELEM_MAXNB];

    int elemtype[OPTSYST_ELEM_MAXNB];       /**< Element type code (1, 3, 4, or 5) */
    int elemarrayindex[OPTSYST_ELEM_MAXNB]; /**< Sub-array index or mask image ID */

    double flux[OPTSYST_ELEM_MAXNB];     /**< Integrated beam flux after element */
    double elemZpos[OPTSYST_ELEM_MAXNB]; /**< Axial conjugation in collimated space [m] */
    int    keepMem[OPTSYST_ELEM_MAXNB];  /**< Retain intermediate WF cubes in memory (0/1) */

    long elem_amp_ID_array[OPTSYST_ELEM_MAXNB]; /**< Per-element multiplicative amplitude ID */
    long elem_pha_ID_array[OPTSYST_ELEM_MAXNB]; /**< Per-element additive phase cube ID [rad] */

    int endmode; /**< 0: compute final focal-plane FFT PSF, 1: stop at last pupil plane */
    int SAVE;    /**< 1: write per-plane diagnostic FITS cubes to disk, 0: in-memory only */
} OPTSYST;

/* Propagate a multi-wavelength complex wavefront cube across axial distance zprop [m]. */
errno_t OptSystProp_propagateCube(
    OPTSYST    *optsyst,
    long        index,
    const char *IDin_amp_name,
    const char *IDin_pha_name,
    const char *IDout_amp_name,
    const char *IDout_pha_name,
    double      zprop,
    int         sharedmem);

/* Execute element-by-element Fresnel propagation from elemstart to elemend. */
errno_t OptSystProp_run(
    OPTSYST    *optsyst,
    long        index,
    long        elemstart,
    long        elemend,
    const char *savedir,
    int         sharedmem);

#endif /* _OPTSYSTPROP_H */
