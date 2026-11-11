/**
 * @file    mkFPM_zonemap.h
 * @brief   2D hexagonal ring and sector zone-index map generator for physical FPMs
 */

#ifndef PIAACMCSIMUL_FOCALPLANEMASK_MKFPM_ZONEMAP_H
#define PIAACMCSIMUL_FOCALPLANEMASK_MKFPM_ZONEMAP_H

/* Build high-resolution 2D integer zone index map IDname for physical FPM geometry. */
errno_t mkFPM_zonemap(
    const char *restrict IDname,
    imageID             *outID);

#endif /* PIAACMCSIMUL_FOCALPLANEMASK_MKFPM_ZONEMAP_H */
