/**
 * @file    PIAAACMCsimul.c
 * @brief   PIAA-type coronagraph design
 * 
 * Can design both APLCMC and PIAACMC coronagraphs
 *  
 * @author  O. Guyon
 * @date    21 nov 2017
 *
 * 
 * @bug No known bugs.
 * 
 */


// System include

#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>

#ifdef __MACH__
#include <mach/mach_time.h>
#define CLOCK_REALTIME 0
#define CLOCK_MONOTONIC 0
int clock_gettime(int clk_id, struct mach_timespec *t){
    mach_timebase_info_data_t timebase;
    mach_timebase_info(&timebase);
    uint64_t time;
    time = mach_absolute_time();
    double nseconds = ((double)time * (double)timebase.numer)/((double)timebase.denom);
    double seconds = ((double)time * (double)timebase.numer)/((double)timebase.denom * 1e9);
    t->tv_sec = seconds;
    t->tv_nsec = nseconds;
    return 0;
}
#else
#include <time.h>
#endif

// External libraries

#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <fitsio.h>


// milk includes
//   core modules
#include "CommandLineInterface/CLIcore.h"
#include "00CORE/00CORE.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
//   other modules
#include "info/info.h"
#include "fft/fft.h"
#include "image_gen/image_gen.h"
#include "WFpropagate/WFpropagate.h"
#include "statistic/statistic.h"
#include "linopt_imtools/linopt_imtools.h"
#include "OpticsMaterials/OpticsMaterials.h"
#include "image_filter/image_filter.h"
#include "image_basic/image_basic.h"
#include "coronagraphs/coronagraphs.h"
#include "PIAACMCsimul/PIAACMCsimul.h"
#include "OptSystProp/OptSystProp.h"



# ifdef HAVE_LIBGOMP
#include <omp.h>
#define OMP_NELEMENT_LIMIT 1000000
#endif



static int INITSTATUS_PIAACMCsimul = 0;



PIAACMCsimul_varType piaacmcsimul_var; 

/// optical system description
OPTSYST *optsyst;


OPTPIAACMCDESIGN *piaacmc;











// command line interface (CLI) commands
//
// function CLI_checkarg used to check arguments
// 1: float
// 2: long
// 3: string, not existing image
// 4: existing image
// 5: string 
//



// local function(s)
static double f_evalmask (const gsl_vector *v, void *params);

// for testing
static char flogcomment[200];
//#define PIAASIMUL_LOGFUNC0 // top level
//#define PIAASIMUL_LOGFUNC1 // lower level



/* =============================================================================================== */
/* =============================================================================================== */
/** @name  Command line interface (CLI)
 *  CLI commands
 */
///@{
/* =============================================================================================== */
/* =============================================================================================== */


/* =============================================================================================== */
/*  2. Focal plane mask construction                                                               */
/* =============================================================================================== */

int_fast8_t PIAACMCsimul_rings2sectors_cli(){
    if(CLI_checkarg(1,4)+CLI_checkarg(2,3)+CLI_checkarg(3,3)==0)    {
        PIAACMCsimul_rings2sectors(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.string);
        return 0;    }    else        return 1;
}


/* =============================================================================================== */
/*  4. Lyot Stop(s)                                                              */
/* =============================================================================================== */

int_fast8_t PIAACMCsimul_geomProp_cli(){
	if(CLI_checkarg(1,4)+CLI_checkarg(2,4)+CLI_checkarg(3,3)+CLI_checkarg(4,3)+CLI_checkarg(5,1)+CLI_checkarg(6,1)+CLI_checkarg(7,1)+CLI_checkarg(8,1)+CLI_checkarg(9,1)+CLI_checkarg(10,1)==0)    {
		PIAACMCsimul_geomProp(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.string, data.cmdargtoken[4].val.string, data.cmdargtoken[5].val.numf, data.cmdargtoken[6].val.numf, data.cmdargtoken[7].val.numf, data.cmdargtoken[8].val.numf, data.cmdargtoken[9].val.numf, data.cmdargtoken[10].val.numf);
	} else return 1;
		
}


/* =============================================================================================== */
/*  5. Focal plane mask optimization                                                               */
/* =============================================================================================== */


int_fast8_t PIAACMC_FPMresp_rmzones_cli(){
    if(CLI_checkarg(1,4)+CLI_checkarg(2,3)+CLI_checkarg(3,2)==0)    {
        PIAACMC_FPMresp_rmzones(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.numl);
        return 0;    }    else        return 1;
}

int_fast8_t PIAACMC_FPMresp_resample_cli(){
    if(CLI_checkarg(1,4)+CLI_checkarg(2,3)+CLI_checkarg(3,2)+CLI_checkarg(4,2)==0)    {
        PIAACMC_FPMresp_resample(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.numl, data.cmdargtoken[4].val.numl);
        return 0;    }    else        return 1;
}

/* =============================================================================================== */
/*  6. Focal plane processing                                                                      */
/* =============================================================================================== */

int_fast8_t PIAACMC_FPM_process_cli(){
	if(CLI_checkarg(1,4)+CLI_checkarg(2,5)+CLI_checkarg(3,2)+CLI_checkarg(4,3)==0)    {
        PIAACMC_FPM_process(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.string, data.cmdargtoken[3].val.numl, data.cmdargtoken[4].val.string);
        return 0;    }    else        return 1;
}


/* =============================================================================================== */
/*  7. High level routines                                                                         */
/* =============================================================================================== */

int_fast8_t PIAACMCsimul_run_cli(){
    if(CLI_checkarg(1,3)+CLI_checkarg(2,2)==0)    {
        PIAACMCsimul_run(data.cmdargtoken[1].val.string, data.cmdargtoken[2].val.numl);
        return 0;    }    else        return 1;
}

///@}



void __attribute__ ((constructor)) libinit_PIAACMCsimul()
{
	if ( INITSTATUS_PIAACMCsimul == 0 )
	{
		init_PIAACMCsimul();
		RegisterModule(__FILE__, "coffee", "PIAACMC system simulation");
		INITSTATUS_PIAACMCsimul = 1; 
	}
}





 /** @name MODULE INITIALIZATION
  * Registers CLI commands
 */
///@{
 
int_fast8_t init_PIAACMCsimul()
{
	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif
	
    strcpy(data.cmd[data.NBcmd].key,"piaacmcsimring2sect");
    strcpy(data.cmd[data.NBcmd].module, __FILE__);
    data.cmd[data.NBcmd].fp = PIAACMCsimul_rings2sectors_cli;
    strcpy(data.cmd[data.NBcmd].info,"turn ring fpm design into sectors");
    strcpy(data.cmd[data.NBcmd].syntax,"<input ring fpm> <zone-ring table> <output sector fpm>");
    strcpy(data.cmd[data.NBcmd].example,"piaacmcsimring2sect");
    strcpy(data.cmd[data.NBcmd].Ccall,"long PIAACMCsimul_rings2sectors(const char *IDin_name, const char *sectfname, const char *IDout_name)");
    data.NBcmd++;


    strcpy(data.cmd[data.NBcmd].key,"piaacmcsimrun");
    strcpy(data.cmd[data.NBcmd].module, __FILE__);
    data.cmd[data.NBcmd].fp = PIAACMCsimul_run_cli;
    strcpy(data.cmd[data.NBcmd].info,"Simulate PIAACMC");
    strcpy(data.cmd[data.NBcmd].syntax,"<configuration index [string]> <mode[int]>");
    strcpy(data.cmd[data.NBcmd].example,"piaacmcsimrun");
    strcpy(data.cmd[data.NBcmd].Ccall,"int PIAACMCsimul_run(const char *confindex, long mode)");
    data.NBcmd++;


    strcpy(data.cmd[data.NBcmd].key,"piaacmsimfpmresprm");
    strcpy(data.cmd[data.NBcmd].module, __FILE__);
    data.cmd[data.NBcmd].fp = PIAACMC_FPMresp_rmzones_cli;
    strcpy(data.cmd[data.NBcmd].info,"remove zones in FPM resp matrix");
    strcpy(data.cmd[data.NBcmd].syntax,"<input FPMresp> <output FPMresp> <NBzone removed>");
    strcpy(data.cmd[data.NBcmd].example,"piaacmsimfpmresprm FPMresp FPMrespout 125");
    strcpy(data.cmd[data.NBcmd].Ccall,"long PIAACMC_FPMresp_rmzones(const char *FPMresp_in_name, const char *FPMresp_out_name, long NBzones)");
    data.NBcmd++;

    strcpy(data.cmd[data.NBcmd].key,"piaacmsimfpmresprs");
    strcpy(data.cmd[data.NBcmd].module, __FILE__);
    data.cmd[data.NBcmd].fp = PIAACMC_FPMresp_resample_cli;
    strcpy(data.cmd[data.NBcmd].info,"resample FPM resp matrix");
    strcpy(data.cmd[data.NBcmd].syntax,"<input FPMresp> <output FPMresp> <NBlambda> <EvalPts step>");
    strcpy(data.cmd[data.NBcmd].example,"piaacmsimfpmresprs FPMresp FPMrespout 10 2");
    strcpy(data.cmd[data.NBcmd].Ccall,"long PIAACMC_FPMresp_resample(const char *FPMresp_in_name, const char *FPMresp_out_name, long NBlambda, long PTstep)");
    data.NBcmd++;

    strcpy(data.cmd[data.NBcmd].key,"piaacmcfpmprocess");
    strcpy(data.cmd[data.NBcmd].module, __FILE__);
    data.cmd[data.NBcmd].fp = PIAACMC_FPM_process_cli;
    strcpy(data.cmd[data.NBcmd].info,"Quantize FPM");
    strcpy(data.cmd[data.NBcmd].syntax,"<input FPM sags> <sectors ASCII file> <number of exposures> <output FPM sags>");
    strcpy(data.cmd[data.NBcmd].example,"piaacmcfpmprocess");
    strcpy(data.cmd[data.NBcmd].Ccall,"long PIAACMC_FPM_process(const char *FPMsag_name, const char *zonescoord_name, long NBexp, const char *outname)");
    data.NBcmd++;

    strcpy(data.cmd[data.NBcmd].key,"piaacmcgeomprop");
    strcpy(data.cmd[data.NBcmd].module, __FILE__);
    data.cmd[data.NBcmd].fp = PIAACMCsimul_geomProp_cli;
    strcpy(data.cmd[data.NBcmd].info,"Geometric propagation from surface");
    strcpy(data.cmd[data.NBcmd].syntax,"<input map> <input sag> <output map> <output map counter> <delta refractive index> <pixel scale> <propagation dist> <kernel radius> <kernel step> <clear aperture>");
    strcpy(data.cmd[data.NBcmd].example,"piaacmcgeomprop pupin piaa0z pupout cntout 2.0 0.00011 2.302606 3.0 0.5 200.0");
	strcpy(data.cmd[data.NBcmd].Ccall,"long PIAACMCsimul_geomProp(const char *IDin_name, const char *IDsag_name, const char *IDout_name, const char *IDoutcnt_name, float drindex, float pscale, float zprop, float krad, float kstep, float rlim)");
	data.NBcmd++;


	piaacmcsimul_var.optsystinit = 0;
	
	// this makes 20% bandwidth, from 0.55/1.1 to 0.55*1.1
	piaacmcsimul_var.LAMBDASTART = 0.5e-6;
	piaacmcsimul_var.LAMBDAEND = 0.605e-6;

	piaacmcsimul_var.FORCE_CREATE_Cmodes = 0;
	piaacmcsimul_var.CREATE_Cmodes = 0;
	piaacmcsimul_var.FORCE_CREATE_Fmodes = 0;
	piaacmcsimul_var.CREATE_Fmodes = 0;

	piaacmcsimul_var.FORCE_CREATE_fpmzmap = 0;
	piaacmcsimul_var.CREATE_fpmzmap = 0;
	piaacmcsimul_var.FORCE_CREATE_fpmzt = 0;
	piaacmcsimul_var.CREATE_fpmzt = 0;

	piaacmcsimul_var.FORCE_CREATE_fpmza = 0;
	piaacmcsimul_var.CREATE_fpmza;

	piaacmcsimul_var.FORCE_MAKE_PIAA0shape = 0;
	piaacmcsimul_var.MAKE_PIAA0shape = 0;
	piaacmcsimul_var.FORCE_MAKE_PIAA1shape = 0;
	piaacmcsimul_var.MAKE_PIAA1shape = 0;
	
	piaacmcsimul_var.focmMode = -1; // if != -1, compute only impulse response to corresponding zone
	piaacmcsimul_var.PIAACMC_FPMsectors = 0;

	piaacmcsimul_var.FPMSCALEFACTOR = 0.9;
	piaacmcsimul_var.PIAACMC_MASKRADLD = 0.0;
	
	piaacmcsimul_var.LOOPCNT = 0;
	
	piaacmcsimul_var.CnormFactor = 1.0;

	piaacmcsimul_var.computePSF_FAST_FPMresp = 0;
	piaacmcsimul_var.computePSF_ResolvedTarget = 0; // source size = 1e-{0.1*computePSF_ResolvedTarget}
	piaacmcsimul_var.computePSF_ResolvedTarget_mode = 0; // 0: source is simulated as 3 points, 1: source is simulated as 6 points
	piaacmcsimul_var.PIAACMC_FPM_FASTDERIVATIVES = 0;

	piaacmcsimul_var.SCORINGTOTAL = 1.0;
	piaacmcsimul_var.MODampl = 1.0e-6;
	piaacmcsimul_var.SCORINGMASKTYPE = 0;
	piaacmcsimul_var.PIAACMC_save = 1;
	//piaacmcsimul_var.PIAACMC_MASKregcoeff = 1.0;
	piaacmcsimul_var.PIAACMC_fpmtype = 0;
	
	piaacmcsimul_var.WRITE_OK = 1;
	
	piaacmcsimul_var.LINOPT = 0;
	
    // add atexit functions here
    atexit(PIAACMCsimul_free);

    return 0;

}


///@}


// first argument should be "PIAACMCsimul.fcall.log"
// second argument should be __FUNCTION__
// PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, "");
static void PIAACMCsimul_logFunctionCall(char *LogFileName, const char *FunctionName, long line, char *comments)
{
	FILE *fp;
	time_t tnow;
	struct tm *uttime;
	struct timespec timenow;

	char string[21];
	
	tnow = time(NULL);
    uttime = gmtime(&tnow);
	clock_gettime(CLOCK_REALTIME, &timenow);
	
	// add custom parameter
	if(piaacmc == NULL)
		sprintf(string, "NULL");
	else
		sprintf(string, "%20ld", piaacmc[0].focmNBzone);
	
	fp = fopen(LogFileName, "a");
	fprintf(fp, "%02d:%02d:%02ld.%09ld  %10d  %40s %6ld   %20s %s\n", uttime->tm_hour, uttime->tm_min, timenow.tv_sec % 60, timenow.tv_nsec, getpid(), FunctionName, line, string, comments);
	fclose(fp);
}





















































































































