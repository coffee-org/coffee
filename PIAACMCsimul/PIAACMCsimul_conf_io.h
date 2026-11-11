/**
 * @file    PIAACMCsimul_conf_io.h
 * @brief   Configuration parameter file and CLI variable ingestion for PIAACMCsimul
 */

#ifndef PIAACMCSIMUL_CONF_IO_H
#define PIAACMCSIMUL_CONF_IO_H

/* Populate piaacmcopticaldesign from defaults, CLI variables, and conf_*.txt files. */
errno_t PIAACMCsimul_initpiaacmcconf_readconfparams(
    long   piaacmctype,
    double fpmradld,
    double centobs0,
    double centobs1,
    int    WFCmode);

/* Read or initialize conf_PIAAmaterial_name.txt and write conf_PIAAmaterial_code.txt. */
errno_t PIAACMCsimul_read_piaa_material(void);

#endif /* PIAACMCSIMUL_CONF_IO_H */
