#ifndef _OPTSYSTPROP_H
#define _OPTSYSTPROP_H

#ifndef __STDC_LIB_EXT1__
typedef int errno_t;
#endif

// ************************************************************************
// ------------------- DEFINITION OF OPTICAL ELEMENTS ---------------------
// ************************************************************************

//
// -------- DM -----------
// square grid geom only
//
typedef struct
{
    long   NBact1D;
    double pitch; // [m]
    double
    maxstroke; // max deviation; actuator moves from -maxstroke to +maxstroke [m]
    long dispID; // points to displacement matrix

    long   IF_ID;      // points to influence function map
    double IFpixscale; // influence function pixel scale [m]
    double IFsize;     // influence function map size (linear)
    double pixscale;   // map pixel scale [m/pix]
    long   dispmapID;  // points to displacement map
} DM_SIM;

//
// --------- aspheric surface mirror --------------
// using same sampling as nominal beam
//
typedef struct
{
    long surfID; // surface Z sag
} ASPHSURFM;

//
// ----------- aspheric surface, refractive ------------
// using same sampling as nominal beam
//
typedef struct
{
    long   mat0; // material before surface - see codes in OpticsMaterial.c
    long   mat1; // material after surface
    int    init; // has refractive index been computed ?
    double ncoeff
    [100]; // for each wavelength (max 100), multiplicative coeff between sag and induced OPD
    long surfID; // surface Z sag
} ASPHSURFR;

// ------- Focal plane mask ------------
typedef struct
{
    long   fpmID;   // 1-focal plane mask complex amplitude cube
    double zfactor; // oversampling factor
    int    mode;    // 1 if 1-CA, 0 if CA
} FOCMASK;

#define OPTSYST_ELEM_MAXNB            100
#define STRINGMAXLEN_OPTSYST_ELEMNAME 100

//
// *****************************************************************************************************
// ------------------------------ structure defining an optical system ---------------------------------
// *****************************************************************************************************

// Fresnel propagation used for simulations
// All elements except focal plane mask are placed in equivalent collimated beam space
// Input pupil is adopted as reference
//
// optical elements are applied sequentially, and consist of amplitude and phase cubes (1 slize per lambda)
// types of elements :
// 1: static opaque mask, defined by image identifyer
// 2: -
// 3: reflective aspheric surface. Defined by OPD map (used for mirrors, including deformable mirrors)
// 4: refractive aspheric surface. Defined by OPD map, material type (index of refraction) before and after surface
// 5: focal plane mask
//
typedef struct
{

    int    nblambda;
    double lambdaarray[2000];

    double beamrad;    // beam radius at input in collimated space [m]
    double pixscale;   // pixel scale in collimated beam [m/pix]
    long   size;       // array size
    long   DFTgridpad; // 0 for full res DFT, >0 for faster coarser DFTs

    // =============== OPTICAL ELEMENTS ===================
    long NBelem; // number of optical elements
    char name[OPTSYST_ELEM_MAXNB][STRINGMAXLEN_OPTSYST_ELEMNAME];

    long      NB_asphsurfm; // max number of aspheric mirrors
    ASPHSURFM ASPHSURFMarray[OPTSYST_ELEM_MAXNB];

    long      NB_asphsurfr; // max number of aspheric refractive surfaces
    ASPHSURFR ASPHSURFRarray[OPTSYST_ELEM_MAXNB];

    long    NB_focmask; // max number of focal plane masks
    FOCMASK FOCMASKarray[OPTSYST_ELEM_MAXNB];

    int elemtype[OPTSYST_ELEM_MAXNB]; // element type

    // if element is DM or aspheric surface, this is the index in the corresponding array of elements, otherwise, this is the image index
    int elemarrayindex[OPTSYST_ELEM_MAXNB];

    double flux[OPTSYST_ELEM_MAXNB];     // total flux AFTER element
    double elemZpos[OPTSYST_ELEM_MAXNB]; // position along beam
    int    keepMem
    [OPTSYST_ELEM_MAXNB]; // set to 1 if memory should be kept, 0 otherwise
    // this is what is used for propagations, created from info above
    long elem_amp_ID_array
    [OPTSYST_ELEM_MAXNB]; // amplitude map identifyer, multiplicative
    long
    elem_pha_ID_array[OPTSYST_ELEM_MAXNB]; // phase map identifyer, additive

    int endmode; // 0: compute PSF at the end of the sequence, 1: no PSF

    int SAVE; // 1 if intermediate results are saved, 0 otherwise

} OPTSYST;

errno_t init_OptSystProp();

errno_t OptSystProp_propagateCube(OPTSYST    *optsyst,
                                  long        index,
                                  const char *IDin_amp_name,
                                  const char *IDin_pha_name,
                                  const char *IDout_amp_name,
                                  const char *IDout_pha_name,
                                  double      zprop,
                                  int         sharedmem);

errno_t OptSystProp_run(OPTSYST    *optsyst,
                        long        index,
                        long        elemstart,
                        long        elemend,
                        const char *savedir,
                        int         sharedmem);

#endif
