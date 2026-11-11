/**
 * @file    exec_optimize_lyot_stop_position.h
 * @brief   Axial conjugation scan for fixed Lyot stop masks (Mode 1)
 */

#ifndef PIAACMCSIMUL_LYOTSTOP_EXEC_OPTIMIZE_LYOT_STOP_POSITION_H
#define PIAACMCSIMUL_LYOTSTOP_EXEC_OPTIMIZE_LYOT_STOP_POSITION_H

/* Scan axial conjugation z [m] to minimize on-axis residual starlight leak (Mode 1). */
errno_t exec_optimize_lyot_stop_position(void);

#endif /* PIAACMCSIMUL_LYOTSTOP_EXEC_OPTIMIZE_LYOT_STOP_POSITION_H */
