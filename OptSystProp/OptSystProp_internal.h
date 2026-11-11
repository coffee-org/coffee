/**
 * @file    OptSystProp_internal.h
 * @brief   Internal optical element application helpers for OptSystProp
 */

#ifndef OPTSYSTPROP_INTERNAL_H
#define OPTSYSTPROP_INTERNAL_H

#include "COREMOD_memory/COREMOD_memory.h"
#include "OptSystProp/OptSystProp.h"

/**
 * @brief Per-element execution context passed to internal element kernels.
 */
typedef struct
{
    OPTSYST    *optsyst;       /**< Optical system array */
    long        index;         /**< Active optical system index */
    long        elem;          /**< Current optical element index */
    imageID     IDa;           /**< Output amplitude cube image ID */
    imageID     IDp;           /**< Output phase cube image ID */
    const char *imnameamp_out; /**< Output amplitude cube image name */
    const char *imnamepha_out; /**< Output phase cube image name */
    const char *savedir;       /**< Diagnostic output directory */
    int         sharedmem;     /**< Shared-memory retention flag */
} OptSystProp_elem_ctx;

/* Multiply complex wavefront amplitude cube by static/chromatic mask element. */
errno_t OptSystProp_apply_ampmask(const OptSystProp_elem_ctx *ctx);

/* Apply reflective aspheric mirror surface OPD phase delay across wavelengths. */
errno_t OptSystProp_apply_mirror(const OptSystProp_elem_ctx *ctx);

/* Apply refractive aspheric surface dispersive OPD phase delay across wavelengths. */
errno_t OptSystProp_apply_refractive(const OptSystProp_elem_ctx *ctx);

/* Apply focal plane mask element via zoomed pupil-to-focal-to-pupil 2D DFT. */
errno_t OptSystProp_apply_focmask(const OptSystProp_elem_ctx *ctx);

/* Dispatch optical element application by elemtype (1, 3, 4, or 5). */
errno_t OptSystProp_apply_element(const OptSystProp_elem_ctx *ctx);

#endif /* OPTSYSTPROP_INTERNAL_H */
