/**
 * @file    mkFocalPlaneMask.h
 * @brief   Polychromatic complex-amplitude focal plane mask cube synthesis
 */

#ifndef PIAACMCSIMUL_MKFOCALPLANEMASK_H
#define PIAACMCSIMUL_MKFOCALPLANEMASK_H

/* Construct 3D polychromatic complex-amplitude FPM cube ID_name from zone map and sags. */
errno_t mkFocalPlaneMask(
    const char *IDzonemap_name,
    const char *ID_name,
    int         mode,
    int         saveMask,
    imageID    *outID);

#endif /* PIAACMCSIMUL_MKFOCALPLANEMASK_H */
