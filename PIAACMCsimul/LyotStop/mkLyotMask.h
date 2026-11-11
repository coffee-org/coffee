/**
 * @file    mkLyotMask.h
 * @brief   Optimal binary Lyot stop synthesis from off-axis / on-axis pupil ratio maps
 */

#ifndef PIAACMCSIMUL_LYOTSTOP_MKLYOTMASK_H
#define PIAACMCSIMUL_LYOTSTOP_MKLYOTMASK_H

/* Construct binary Lyot mask IDout_name achieving target geometric throughput (0..1). */
errno_t mkLyotMask(
    const char *restrict IDincoh_name,
    const char *restrict IDmc_name,
    const char *restrict IDzone_name,
    double               throughput,
    const char *restrict IDout_name,
    imageID             *outID);

#endif /* PIAACMCSIMUL_LYOTSTOP_MKLYOTMASK_H */
