/**
 * @file    PIAACMCsimul_loadsavepiaacmcconf.h
 * @brief   Serialization of piaacmcparams.conf and FPM solution filename tags
 */

#ifndef PIAACMCSIMUL_LOADPIAACMCCONF_H
#define PIAACMCSIMUL_LOADPIAACMCCONF_H

/* Load optical prescription from <dname>/piaacmcparams.conf into piaacmcopticaldesign. */
errno_t PIAACMCsimul_loadpiaacmcconf(const char *dname);

/* Format geometry/bandpass descriptor string piaacmcparams.fnamedescr_conf. */
errno_t PIAACMCsimul_update_fnamedescr_conf(void);

/* Format full FPM sag/regularization descriptor string piaacmcparams.fnamedescr. */
errno_t PIAACMCsimul_update_fnamedescr(void);

/* Write current piaacmcopticaldesign prescription to <dname>/piaacmcparams.conf. */
errno_t PIAACMCsimul_savepiaacmcconf(const char *restrict dname);

#endif /* PIAACMCSIMUL_LOADPIAACMCCONF_H */
