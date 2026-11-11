/**
 * @file    mkSimpleLyotStop.h
 * @brief   Circular annular Lyot stop mask generator
 */

#ifndef MKSIMPLELYOTSTOP_H
#define MKSIMPLELYOTSTOP_H

/* Create annular pupil transmission mask ID_name transmitting rin <= r <= rout [beamrad]. */
errno_t mkSimpleLyotStop(
    const char *restrict ID_name,
    float                rin,
    float                rout,
    imageID             *outID);

#endif /* MKSIMPLELYOTSTOP_H */
