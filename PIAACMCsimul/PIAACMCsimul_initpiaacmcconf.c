/**
 * @file    PIAAACMCsimul_initpiaacmcconf.c
 * @brief   PIAA-type coronagraph design, initialize
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


// System includes

#include <stdint.h>
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
#include <sys/stat.h>
#include <sys/types.h>

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
#include <fitsio.h>


// cfitsTK includes
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




/// All global images and variables 

extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTSYST *optsyst;

extern OPTPIAACMCDESIGN *piaacmc;










/**
 * 
 * # List of Configuration Parameters {#page_PIAACMCsimul_Configuration}
 * 
 * Each parameter is stored on disk as conf/conf_<param>.txt
 * 
 * 
 * 
 * 
 * PIAA OPTICS DESIGN: 
 * 
 * | Parameter              | Meaning                                                  |
 * | ---------------------- | -------------------------------------------------------- |
 * | PIAAmode               | PIAA mode (0: classical apodization, 1: PIAA)            |
 * | fpmradld               | Focal plane mask radius [l/D]                            |
 * | PIAAcoeff              | Fraction of PIAA apodization                             |
 * | coin                   | Central obstruction at input beam [beam radius]          |
 * | coout                  | Central obstruction at output beam [beam radius]
 * | PIAAmaterial           | PIAA optics material
 * | PIAAcirc               | FLAG: 1 if PIAA shapes are circular (no Fourier modes)
 * | REGPIAACCOEFF          | PIAA regularization amplitude, cosine modes
 * | REGPIAACALPHA          | PIAA regularization power law, cosine modes
 * | REGPIAAFCOEFF          | PIAA regularization amplitude, Fourier modes
 * | REGPIAAFALPHA          | PIAA regularization power law, Fourier modes
 *
 * 
 * LYOT STOP(S) DESIGN:          
 * 
 * | Parameter              | Meaning                                                  |
 * | ---------------------- | -------------------------------------------------------- |
 * | LStransm               | Lyot stop transmission
 * | NBls                   | Number of Lyot stops
 * | lambda                 | Wavelength for monochromatic design [nm]
 * 
 * 
 * 
 * 
 * FOCAL PLANE MASK DESIGN: 
 * 
 * 
 * | Parameter              | Meaning                                                  |
 * | ---------------------- | -------------------------------------------------------- |
 * | fpmmaterial            | focal plane mask material
 * | FPMsectors             | mask geometry: 0=disk, 1=sectors, 2=hexagonal tiling
 * | NBrings                | number of rings in focal plane mask
 * | maskradld              | mask outer radius at central wavelength [l/D]
 * | fpmminsag              | min focal plane mask sag
 * | fpmmaxsag              | max focal plane mask sag
 * | fpmregsag_coeff        | sag regularization coefficient
 * | fpmregsag_alpha        | sag regularization coefficient exponent
 * | fpmccnbr               | how many central rings replaced by cone (set to 0 if no central cone
 * | fpmccz                 | sag at cone center (sag at cone edge will be midpoint between minsag and maxsag)
 * | fpmocradld             | outer cone outer radius [l/D]
 * | fpmocz                 | sag at inner edge of outer cone (sag = 0 at outer edge), set to 0 if no outer cone
 * 
 * 
 * 
 * OPTIMIZATION PARAMETERS: 
 * 
 * | Parameter              | Meaning                                                  |
 * | ---------------------- | -------------------------------------------------------- |
 * | mlambda                | central wavelength for polychromatic design [nm]
 * | mlambdaB               | spectral bandwidth [%]
 * | nblambda               | Number of wavelength values
 * | ssize                  | source angular size for mask optimization (20: 0.01 l/D; 10: 0.1 l/D)
 * | extmode                | source extent mode (0: 1 point, 1: 3 points; 2: 6 points)
 * 
 * 
 * 
 * OPTICAL DESIGN: 
 * 
 * 
 * | Parameter              | Meaning                                                  |
 * | ---------------------- | -------------------------------------------------------- |
 * | size                   | array size
 * | beamrad                | beam radius [mm]
 * | pscale                 | pixel scale in pupil [m/pix]
 * | Fratio                 | F ratio at focal plane mask
 * | PIAAr0lim              | outer edge of PIAA optic 0 [beam radius unit]
 * | PIAAr1lim              | outer edge of PIAA optic 1 [beam radius unit]
 * | PIAAsep                | distance between PIAA optics [m]
 * | PIAA0pos               | PIAA optic 0 distance from pupil plane [m]
 * | invPIAAmode            | 0: no inv PIAA, 1: inv PIAA after Lyot stops, 2: inv PIAA before Lyot stops
 * | prePIAA0maskpos        | pre-PIAA optic 0 mask distance from pupil plane [m] (if prePIAA0mask.fits exists)
 * | postPIAA0maskpos       | post-PIAA optic 0 mask distance from pupil plane [m] (if postPIAA0mask.fits exits)
 * | piaaNBCmodesmax        | maximum number of radial cosine modes for PIAA optics
 * | piaaCPAmax             | maximum spatial frequency (CPA) for PIAA optics
 * | LyotZmin               | minimum value for Lyot stop(s) conjugation range [m] - relative to element named "post focal plane mask pupil"
 * | LyotZmax               | maximum value for Lyot stop(s) conjugation range [m] - relative to element named "post focal plane mask pupil"
 * | pupoutmaskrad          | output pupil mask radius (scaled to pupil radius)
 * 
 * 
 */
















/**
 * @brief Creates/initializes piaacmcconf structure and directory
 *
 * @param[in] piaacmctype  Type of system: 0=idealized mask, 1=physical mask
 * @param[in] fpmradld     Focal plane mask nominal radius
 * @param[in] centobs0     Input central obstruction
 * @param[in] centobs1     Output central obstruction
 * @param[in] WFSmode      Number of DMs (0: no WFC)
 * @param[in] load         if 1, attempt to load configuration from file
 * 
 * piaacmctype:
 * - 0: if configuration does not exist, create Monochromatic idealized PIAACMC, otherwise, read configuration
 * - 1: physical mask
 *
 * 
 */

int PIAACMCsimul_initpiaacmcconf(
		long piaacmctype, 
		double fpmradld, 
		double centobs0, 
		double centobs1, 
		int WFCmode, 
		int load
		)
{
    FILE *fp;
    float beamradpix;
    long NBpiaacmcdesign = 1;
    long ii, jj, k, i;
    double x, y;
    long size, size0;
    long Cmsize;
    long Fmsize;
    long ID, ID0, ID1;
    long size2;
    double rad;
    char command[1000];
    long IDv1, IDv2;
    char fname[1500];
    char name[500];
    float tmpf;
    double pha0, t0, t;
    int loaded = 0; // has the configuration been loaded ?

    int saveconf = 0; // if 1, save conf at end of this function
    long IDv;
    int ret;
    double tmplf;

    int IDlscumul;
    long iDM; // DM index
    uint32_t *sizearray;


    long IDapo;
    long xsize = 0;
    long ysize = 0;
    long IDapo_PIAA, IDapo_CPA;
    double coeff;

	long IDx, IDy, IDr, IDPA;


	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif


	// Create required directories
	ret = mkdir("testdir", 0777);  // test files, output from code
	ret = mkdir("conf",    0777);  // configuration
	ret = mkdir("status",  0777);  // status
	ret = mkdir("ref",     0777);  // reference files
	ret = mkdir("log",     0777);  // log




    if(piaacmc == NULL)
    {
        piaacmc = (OPTPIAACMCDESIGN*) malloc(sizeof(OPTPIAACMCDESIGN)*NBpiaacmcdesign);


        // Default Values for PIAACMC (will adopt them unless configuration file exists)
        piaacmc[0].nblambda = 8;


        // high resolution
        //piaacmc[0].size = 4096;
        //piaacmc[0].pixscale = 0.000055;

        // mid resolution
        piaacmc[0].size = 2048;
        piaacmc[0].pixscale = 0.000055;

        // low resolution
        //piaacmc[0].size = 1024;
        //piaacmc[0].pixscale = 0.00011;


        // very low resolution
        //    piaacmc[0].size = 512;
        // piaacmc[0].pixscale = 0.00022;


        piaacmc[0].beamrad = 0.01; // beam physical radius
        piaacmc[0].PIAA0pos = 1.0; // piaa 0 position [m]
        piaacmc[0].PIAAsep = 1.00; // [m]
        piaacmc[0].fpzfactor = 8.0;
        piaacmc[0].Fratio = 80.0; // default
        strcpy(piaacmc[0].PIAAmaterial_name, "Mirror");  // mirrors
        piaacmc[0].prePIAA0mask = 0;
        piaacmc[0].prePIAA0maskpos = 0.0;
        piaacmc[0].postPIAA0mask = 0;
        piaacmc[0].postPIAA0maskpos = 0.0;
        piaacmc[0].piaaNBCmodesmax =  40;
        piaacmc[0].piaaCPAmax = 10.0;

        piaacmc[0].centObs0 = centobs0; // input central obstruction
        piaacmc[0].centObs1 = centobs1; // output central obstruction
        piaacmc[0].NBradpts = 50000;
        piaacmc[0].r0lim = 1.15; // outer radius after extrapolation, piaa optics 0
        piaacmc[0].r1lim = 1.5; // outer radius after extrapolation, piaa optics 1


        /// Wavefront control
        piaacmc[0].nbDM = WFCmode; // number of deformable mirrors (10 max)
        for(iDM=0; iDM<piaacmc[0].nbDM; iDM++)
        {
            piaacmc[0].DMpos[iDM] = 0.0 + 0.6*iDM/(0.01+piaacmc[0].nbDM-1.0); // DM conjugation in collimated space
            piaacmc[0].ID_DM[iDM] = -1;  // DM image identifier - to be updated later
        }


        piaacmc[0].NBLyotStop = 2;
        for(i=0; i<10; i++)
        {
            piaacmc[0].LyotStop_zpos[i] = 0.0;
            piaacmc[0].IDLyotStop[i] = -1;
        }

        piaacmc[0].fpmaskradld = fpmradld; // to compute prolate spheroidal function
        piaacmc[0].fpmarraysize = 2048;

        piaacmc[0].fpmRad = 100.0e-6; // focal plane radius [m]
        piaacmc[0].NBrings = 4; // number of rings in focal plane mask
        piaacmc[0].fpmminsag = -1.0e-5;
        piaacmc[0].fpmmaxsag = 1.0e-5;
        piaacmc[0].fpmsagreg_coeff = 1.0;
        piaacmc[0].fpmsagreg_alpha = 1.0;
        piaacmc[0].NBringCentCone = 0; // central cone
        piaacmc[0].fpmCentConeZ = piaacmc[0].fpmmaxsag;
        piaacmc[0].fpmOuterConeZ = piaacmc[0].fpmmaxsag;
        piaacmc[0].fpmOuterConeRadld = 80.0;
        piaacmc[0].fpmmaterial_code = 0;  // 0: mirror
        piaacmc[0].fpmaskamptransm = 1.0;


        if((fp = fopen("conf/conf_peakPSF.txt", "r"))!=NULL)
        {
            ret = fscanf(fp, "%f", &piaacmc[0].peakPSF);
            fclose(fp);
        }
        else
            piaacmc[0].peakPSF = -1.0;



        piaacmc[0].PIAAmode = 1;
        if((IDv=variable_ID("PIAACMC_PIAAmode"))!=-1)
            piaacmc[0].PIAAmode = (int) (data.variable[IDv].value.f+0.01);

        piaacmc[0].PIAAcoeff = 1.0;
        if((IDv=variable_ID("PIAACMC_PIAAcoeff"))!=-1)
            piaacmc[0].PIAAcoeff = data.variable[IDv].value.f;



        if(piaacmc[0].PIAAmode == 0)
        {
            piaacmc[0].invPIAAmode = 0;
        }
        else
        {
            piaacmc[0].invPIAAmode = 1;
            if((IDv=variable_ID("PIAACMC_invPIAAmode"))!=-1)
                piaacmc[0].invPIAAmode = (long) (data.variable[IDv].value.f+0.001);
        }



        
                sprintf(fname, "%s/conf_fpmmaterial_name.txt", piaacmcsimul_var.piaacmcconfdir );
                if( (fp = fopen(fname, "r")) != NULL)
                {
                    ret = fscanf(fp, "%s", name);
                    strcpy(piaacmc[0].fpmmaterial_name, name);
                    printf("Reading %s   piaacmc[0].fpmmaterial_name : %s\n", fname, piaacmc[0].fpmmaterial_name);
                    fclose(fp);
                }
                else
                {
                    sprintf(piaacmc[0].fpmmaterial_name, "Mirror");
                    sprintf(fname, "%s/conf_fpmmaterial_name.txt", piaacmcsimul_var.piaacmcconfdir);
                    printf("Writing %s   piaacmc[0].fpmmaterial_name : %s\n", fname, piaacmc[0].fpmmaterial_name);
                    if((fp=fopen(fname,"w"))!=NULL)
                    {
                        fprintf(fp, "%s\n", piaacmc[0].fpmmaterial_name);
                        fclose(fp);
                    }
                    else
                    {
                        printf("ERROR: cannot create file \"%s\"\n", fname);
                        exit(0);
                    }
                }

                printf("piaacmc[0].fpmmaterial_name : %s\n", piaacmc[0].fpmmaterial_name);
                piaacmc[0].fpmmaterial_code = OPTICSMATERIALS_code(piaacmc[0].fpmmaterial_name);

                sprintf(fname, "%s/conf_fpmmaterial_code.txt", piaacmcsimul_var.piaacmcconfdir);
                if((fp=fopen(fname,"w"))!=NULL)
                {
                    fprintf(fp, "%d\n", piaacmc[0].fpmmaterial_code);
                    fclose(fp);
                }
                else
                {
                    printf("ERROR: cannot create file \"%s\"\n", fname);
                    exit(0);
                }
        



        if((IDv=variable_ID("PIAACMC_beamrad"))!=-1)
            piaacmc[0].beamrad = data.variable[IDv].value.f; // beam physical radius

        if((IDv=variable_ID("PIAACMC_Fratio"))!=-1)
            piaacmc[0].Fratio = data.variable[IDv].value.f; // Focal ratio
        if((IDv=variable_ID("PIAACMC_r0lim"))!=-1)
            piaacmc[0].r0lim = data.variable[IDv].value.f;
        if((IDv=variable_ID("PIAACMC_r1lim"))!=-1)
            piaacmc[0].r1lim = data.variable[IDv].value.f;

        if((IDv=variable_ID("PIAACMC_PIAAsep"))!=-1)
            piaacmc[0].PIAAsep = data.variable[IDv].value.f; // piaa separation
        if((IDv=variable_ID("PIAACMC_PIAA0pos"))!=-1)
            piaacmc[0].PIAA0pos = data.variable[IDv].value.f; // piaa elem 0 position



        if((IDv=variable_ID("PIAACMC_prePIAA0maskpos"))!=-1)
            piaacmc[0].prePIAA0maskpos = data.variable[IDv].value.f; // pre piaa elem 0 mask position
        if((IDv=variable_ID("PIAACMC_postPIAA0maskpos"))!=-1)
            piaacmc[0].postPIAA0maskpos = data.variable[IDv].value.f; // post piaa elem 0 mask position



        piaacmc[0].LyotZmin = -3.0;
        if((IDv=variable_ID("PIAACMC_LyotZmin"))!=-1)
            piaacmc[0].LyotZmin = data.variable[IDv].value.f;
        piaacmc[0].LyotZmax = 3.0;
        if((IDv=variable_ID("PIAACMC_LyotZmax"))!=-1)
            piaacmc[0].LyotZmax = data.variable[IDv].value.f;

        piaacmc[0].pupoutmaskrad = 0.95;
        if((IDv=variable_ID("PIAACMC_pupoutmaskrad"))!=-1)
            piaacmc[0].pupoutmaskrad = data.variable[IDv].value.f;




        if((IDv=variable_ID("PIAACMC_piaaNBCmodesmax"))!=-1)
            piaacmc[0].piaaNBCmodesmax = (long) (data.variable[IDv].value.f +0.01); // max number of Cosine terms
        if((IDv=variable_ID("PIAACMC_piaaCPAmax"))!=-1)
            piaacmc[0].piaaCPAmax = data.variable[IDv].value.f; // max CPA for PIAA shapes tuning


        piaacmc[0].NBLyotStop = 1;
        if(piaacmc[0].PIAAmode == 0)
        {
            piaacmc[0].NBLyotStop = 1;
        }
        else
        {
            if((IDv=variable_ID("PIAACMC_nblstop"))!=-1)
                piaacmc[0].NBLyotStop = (long) data.variable[IDv].value.f+0.01;
        }




        if((IDv=variable_ID("PIAACMC_lambda"))!=-1)
            piaacmc[0].lambda = 1.0e-9*data.variable[IDv].value.f; // central wavelength [m]
        //             printf("lambda = %g\n", piaacmc[0].lambda);

        if((IDv=variable_ID("PIAACMC_lambdaB"))!=-1)
            piaacmc[0].lambdaB = data.variable[IDv].value.f; // spectral bandwidth [%]

        piaacmcsimul_var.LAMBDASTART = piaacmc[0].lambda * (1.0 - 0.005*piaacmc[0].lambdaB);
        piaacmcsimul_var.LAMBDAEND = piaacmc[0].lambda * (1.0 + 0.005*piaacmc[0].lambdaB);


        if((IDv=variable_ID("PIAACMC_nblambda"))!=-1)
            piaacmc[0].nblambda = data.variable[IDv].value.f;

        if((IDv=variable_ID("PIAACMC_NBrings"))!=-1)
            piaacmc[0].NBrings = data.variable[IDv].value.f;

        if((IDv=variable_ID("PIAACMC_fpmminsag"))!=-1)
            piaacmc[0].fpmminsag = data.variable[IDv].value.f;
        if((IDv=variable_ID("PIAACMC_fpmmaxsag"))!=-1)
            piaacmc[0].fpmmaxsag = data.variable[IDv].value.f;

        if((IDv=variable_ID("PIAACMC_fpmsagreg_coeff"))!=-1)
            piaacmc[0].fpmsagreg_coeff = data.variable[IDv].value.f;
        if((IDv=variable_ID("PIAACMC_fpmsagreg_alpha"))!=-1)
            piaacmc[0].fpmsagreg_alpha = data.variable[IDv].value.f;
         
         
         
            
        if((IDv=variable_ID("PIAACMC_NBringCentCone"))!=-1)
            piaacmc[0].NBringCentCone = data.variable[IDv].value.f;

        if((IDv=variable_ID("PIAACMC_fpmCentConeZ"))!=-1)
            piaacmc[0].fpmCentConeZ = data.variable[IDv].value.f;
        if((IDv=variable_ID("PIAACMC_fpmOuterConeZ"))!=-1)
            piaacmc[0].fpmOuterConeZ = data.variable[IDv].value.f;
        if((IDv=variable_ID("PIAACMC_fpmOuterConeRadld"))!=-1)
            piaacmc[0].fpmOuterConeRadld = data.variable[IDv].value.f;
        piaacmc[0].fpmOuterConeRad = 0.5*(piaacmcsimul_var.LAMBDASTART + piaacmcsimul_var.LAMBDAEND)*piaacmc[0].Fratio*piaacmc[0].fpmOuterConeRadld;  // [l/D] radius

        if((IDv=variable_ID("PIAACMC_size"))!=-1)
            piaacmc[0].size = (long) (data.variable[IDv].value.f+0.01);

        if((IDv=variable_ID("PIAACMC_pixscale"))!=-1)
            piaacmc[0].pixscale = data.variable[IDv].value.f;


        if(piaacmctype==0) // idealized focal plane mask
        {
            piaacmcsimul_var.FORCE_CREATE_fpmzt = 1; // force making the focal plane mask
            piaacmc[0].NBrings = 1;
            piaacmc[0].NBringCentCone = 0;
            piaacmc[0].fpmOuterConeZ = 0.0;
            piaacmc[0].fpmminsag = -1e-5;
            piaacmc[0].fpmmaxsag = 1e-5;
            piaacmc[0].fpmsagreg_coeff = 1.0;
            piaacmc[0].fpmsagreg_alpha = 1.0;
            piaacmc[0].fpmRad = 0.5*( piaacmcsimul_var.LAMBDASTART + piaacmcsimul_var.LAMBDAEND)*piaacmc[0].Fratio*fpmradld;  // [l/D] radius
            printf("Idealized focal plane mask  radius = %f l/D  = %g m    [lambda = %g (%g - %g)] [Fratio = %f]\n", fpmradld, piaacmc[0].fpmRad, 0.5*(piaacmcsimul_var.LAMBDASTART + piaacmcsimul_var.LAMBDAEND), piaacmcsimul_var.LAMBDASTART, piaacmcsimul_var.LAMBDAEND, piaacmc[0].Fratio);
        }
        else
        {
            if( piaacmcsimul_var.PIAACMC_MASKRADLD < 0.2) // not initialized
                piaacmcsimul_var.PIAACMC_MASKRADLD = 0.1*( (long) (10.0*1.2*fpmradld)); // 1.2x nominal radius, rounded to nearest 0.1 l/D

            piaacmc[0].fpmRad = 0.5*(piaacmcsimul_var.LAMBDASTART + piaacmcsimul_var.LAMBDAEND)*piaacmc[0].Fratio * piaacmcsimul_var.PIAACMC_MASKRADLD;
            printf("Physical focal plane mask - rad = %f l/D -> %g m   [lambda = %g (%g - %g)] [Fratio = %f]\n", piaacmcsimul_var.PIAACMC_MASKRADLD, piaacmc[0].fpmRad, 0.5*(piaacmcsimul_var.LAMBDASTART + piaacmcsimul_var.LAMBDAEND), piaacmcsimul_var.LAMBDASTART, piaacmcsimul_var.LAMBDAEND, piaacmc[0].Fratio);
        }

        piaacmc[0].CmodesID = -1; // Cosine radial mode
        piaacmc[0].FmodesID = -1; // Fourier 2D modes
        piaacmc[0].piaa0CmodesID = -1;
        piaacmc[0].piaa0FmodesID = -1;
        piaacmc[0].piaa1CmodesID = -1;
        piaacmc[0].piaa1FmodesID = -1;
        piaacmc[0].zonezID = -1;  // focm zone material thickness, double precision image
        piaacmc[0].zoneaID = -1;  // focm zone amplitude transmission, double precision image
    }
    

    piaacmc[0].fpmCentConeRad = piaacmc[0].fpmRad*piaacmc[0].NBringCentCone/piaacmc[0].NBrings;

    printf("fpmRad = %g m\n", 0.5*( piaacmcsimul_var.LAMBDASTART + piaacmcsimul_var.LAMBDAEND)*piaacmc[0].Fratio*fpmradld);
    printf("factor = %f   (%ld / %ld)\n", 1.0*piaacmc[0].NBringCentCone/piaacmc[0].NBrings, piaacmc[0].NBringCentCone, piaacmc[0].NBrings);
    printf("fpmCentConeRad =  %g\n", (0.5*( piaacmcsimul_var.LAMBDASTART + piaacmcsimul_var.LAMBDAEND)*piaacmc[0].Fratio*fpmradld)*piaacmc[0].NBringCentCone/piaacmc[0].NBrings);




    if(load==1)
    {
        printf("Loading PIAACMC configuration\n");
        fflush(stdout);
        sprintf(command, "mkdir -p %s", piaacmcsimul_var.piaacmcconfdir);
        ret = system(command);
        loaded = PIAACMCsimul_loadpiaacmcconf(piaacmcsimul_var.piaacmcconfdir);
        if(loaded==0)
        {
            printf("Saving default configuration\n");
            fflush(stdout);
            saveconf = 1;
        }
    }





    sprintf(fname, "%s/conf_PIAAmaterial_name.txt", piaacmcsimul_var.piaacmcconfdir );
    if( (fp = fopen(fname, "r")) != NULL)
    {
        ret = fscanf(fp, "%s", name);
        strcpy(piaacmc[0].PIAAmaterial_name, name);
        fclose(fp);
    }
    else
    {
        sprintf(piaacmc[0].PIAAmaterial_name, "Mirror");
        sprintf(fname, "%s/conf_PIAAmaterial_name.txt", piaacmcsimul_var.piaacmcconfdir);
        if((fp=fopen(fname,"w"))!=NULL)
        {
            fprintf(fp, "%s\n", piaacmc[0].PIAAmaterial_name);
            fclose(fp);
        }
        else
        {
            printf("ERROR: cannot create file \"%s\"\n", fname);
            exit(0);
        }
    }

    printf("piaacmc[0].PIAAmaterial_name = %s\n", piaacmc[0].PIAAmaterial_name);
    piaacmc[0].PIAAmaterial_code = OPTICSMATERIALS_code(piaacmc[0].PIAAmaterial_name);

    sprintf(fname, "%s/conf_PIAAmaterial_code.txt", piaacmcsimul_var.piaacmcconfdir);
    if((fp=fopen(fname,"w"))!=NULL)
    {
        fprintf(fp, "%d\n", piaacmc[0].PIAAmaterial_code);
        fclose(fp);
    }
    else
    {
        printf("ERROR: cannot create file \"%s\"\n", fname);
        exit(0);
    }



    printf("lambda = %g\n", piaacmc[0].lambda);
    printf("LAMBDASTART = %g\n", piaacmcsimul_var.LAMBDASTART);
    printf("LAMBDAEND = %g\n", piaacmcsimul_var.LAMBDAEND);


    for(k=0; k<piaacmc[0].nblambda; k++)
        piaacmc[0].lambdaarray[k] = piaacmcsimul_var.LAMBDASTART + (0.5+k)*(piaacmcsimul_var.LAMBDAEND - piaacmcsimul_var.LAMBDASTART)/piaacmc[0].nblambda;


    // create modes for aspheric optical surfaces description
    beamradpix = piaacmc[0].beamrad/piaacmc[0].pixscale;
    size = piaacmc[0].size;

    printf("BEAM RADIUS :  %f / %f =  %f pix, size = %ld\n", piaacmc[0].beamrad, piaacmc[0].pixscale, beamradpix, size);
    fflush(stdout);

    // x, y, r and PA coordinates in beam (for convenience & speed)
    IDx = create_2Dimage_ID("xcoord", size, size);
    IDy = create_2Dimage_ID("ycoord", size, size);
    IDr = create_2Dimage_ID("rcoord", size, size);
    IDPA = create_2Dimage_ID("PAcoord", size, size);
    printf("pre-computing x, y, r, and PA\n");
    fflush(stdout);
   // list_image_ID();

    for(ii=0; ii<size; ii++)
    {
        for(jj=0; jj<size; jj++)
        {
            x = (1.0*ii-0.5*size)/beamradpix;
            y = (1.0*jj-0.5*size)/beamradpix;
            data.image[IDx].array.F[jj*size+ii] = x;
            data.image[IDy].array.F[jj*size+ii] = y;
            data.image[IDr].array.F[jj*size+ii] = sqrt(x*x+y*y);
            data.image[IDPA].array.F[jj*size+ii] = atan2(y,x);
        }
    }

    // ==================== CREATE DMs ===============
    printf("%d DM(s)\n", piaacmc[0].nbDM);
    fflush(stdout);
    for(iDM=0; iDM<piaacmc[0].nbDM; iDM++)
    {
        printf("DM # %ld\n",iDM);
        sprintf(fname, "wfcDM%ld", iDM);
        piaacmc[0].ID_DM[iDM] = image_ID(fname);  // DM image identifier - to be updated later
        printf("ID = %ld", piaacmc[0].ID_DM[iDM]);

        if(piaacmc[0].ID_DM[iDM] == -1)
        {
            read_sharedmem_image(fname);
            piaacmc[0].ID_DM[iDM] = image_ID(fname);
        }

        if(piaacmc[0].ID_DM[iDM] == -1)
        {
            sizearray = (uint32_t*) malloc(sizeof(uint32_t)*2);
            sizearray[0] = size;
            sizearray[1] = size;
            piaacmc[0].ID_DM[iDM] = create_image_ID(fname, 2, sizearray, _DATATYPE_FLOAT, 1, 0);
            free(sizearray);
        }
    }

    // ==================== CREATE MODES USED TO FIT AND DESCRIBE PIAA SHAPES ===============
    printf("Creating / loading Cmodes and Fmodes ...\n");
    fflush(stdout);

	
    piaacmcsimul_var.CREATE_Cmodes = 0;
    //   sprintf(fname, "%s/Cmodes.fits", piaacmcsimul_var.piaacmcconfdir);
    sprintf(fname, "ref/Cmodes_%ld.fits", piaacmc[0].size);
    if( piaacmcsimul_var.FORCE_CREATE_Cmodes == 0 )
    {
        piaacmc[0].CmodesID = image_ID("Cmodes");
        if(piaacmc[0].CmodesID==-1)
            piaacmc[0].CmodesID = load_fits(fname, "Cmodes", 0);
        if(piaacmc[0].CmodesID==-1)
            piaacmcsimul_var.CREATE_Cmodes = 1;
    }
    else
        piaacmcsimul_var.CREATE_Cmodes = 1;
    if(piaacmcsimul_var.CREATE_Cmodes == 1)
    {
        if(piaacmc[0].CmodesID!=-1)
            delete_image_ID("Cmodes");
        Cmsize = (long) (beamradpix*4);
        if(Cmsize>size)
			Cmsize = size;
        printf("beamradpix = %f -> Cmsize = %ld\n", beamradpix, Cmsize);
        // make sure Cmsize if even
        if (Cmsize%2 == 1)
            Cmsize++;
        piaacmc[0].Cmsize = Cmsize;
        linopt_imtools_makeCosRadModes("Cmodes", Cmsize, piaacmc[0].piaaNBCmodesmax, ApoFitCosFact*beamradpix, 2.0);
        piaacmc[0].CmodesID = image_ID("Cmodes");
        save_fits("Cmodes", fname);
        
        sprintf(command, "mv ModesExpr_CosRad.txt %s/", piaacmcsimul_var.piaacmcconfdir);
        ret = system(command);
    }
    piaacmc[0].NBCmodes = data.image[piaacmc[0].CmodesID].md[0].size[2];
    piaacmc[0].Cmsize = data.image[piaacmc[0].CmodesID].md[0].size[0];



    piaacmcsimul_var.CREATE_Fmodes = 0;
    //    sprintf(fname, "%s/Fmodes.fits", piaacmcsimul_var.piaacmcconfdir);
    sprintf(fname, "ref/Fmodes_%ld.fits", piaacmc[0].size);
    if( piaacmcsimul_var.FORCE_CREATE_Fmodes == 0 )
    {
        piaacmc[0].FmodesID = image_ID("Fmodes");
        if(piaacmc[0].FmodesID==-1)
            piaacmc[0].FmodesID = load_fits(fname, "Fmodes", 0);
        if(piaacmc[0].FmodesID==-1)
            piaacmcsimul_var.CREATE_Fmodes = 1;
    }
    else
        piaacmcsimul_var.CREATE_Fmodes = 1;
    if( piaacmcsimul_var.CREATE_Fmodes == 1 )
    {
        Fmsize = (long) (beamradpix*4);
         if(Fmsize>size)
			Fmsize = size;
        piaacmc[0].Fmsize = Fmsize;
        linopt_imtools_makeCPAmodes("Fmodes",  Fmsize, piaacmc[0].piaaCPAmax, 0.8, beamradpix, 2.0, 1);
        piaacmc[0].FmodesID = image_ID("Fmodes");
        save_fits("Fmodes", fname);
		save_fits("cpamodesfreq", "!cpamodesfreq.fits");
        sprintf(command, "mv ModesExpr_CPA.txt %s/", piaacmcsimul_var.piaacmcconfdir);
        ret = system(command);
    }
    piaacmc[0].NBFmodes = data.image[piaacmc[0].FmodesID].md[0].size[2];
    piaacmc[0].Fmsize = data.image[piaacmc[0].FmodesID].md[0].size[0];

    printf("DONE Creating / loading Cmodes and Fmodes\n");
    fflush(stdout);









    // =================== IMPORT / CREATE PIAA SHAPES =====================
    sprintf(command, "mkdir -p %s/piaaref/", piaacmcsimul_var.piaacmcconfdir);
    ret = system(command);

    if(piaacmc[0].PIAAmode == 1)
    {
        piaacmc[0].piaa0CmodesID = image_ID("piaa0Cmodescoeff");
        piaacmc[0].piaa0FmodesID = image_ID("piaa0Fmodescoeff");
        piaacmc[0].piaa1CmodesID = image_ID("piaa1Cmodescoeff");
        piaacmc[0].piaa1FmodesID = image_ID("piaa1Fmodescoeff");

        sprintf(command, "mkdir -p %s/piaaref/", piaacmcsimul_var.piaacmcconfdir);
        ret = system(command);

        if((piaacmc[0].piaa0CmodesID==-1)||( piaacmc[0].piaa0FmodesID==-1)||(piaacmc[0].piaa1CmodesID==-1)||( piaacmc[0].piaa1FmodesID==-1))
        {
            sprintf(fname, "%s/piaaref/piaa0Cmodes.fits", piaacmcsimul_var.piaacmcconfdir);
            piaacmc[0].piaa0CmodesID = load_fits(fname, "piaa0Cmodescoeff", 1);

            sprintf(fname, "%s/piaaref/piaa0Fmodes.fits", piaacmcsimul_var.piaacmcconfdir);
            piaacmc[0].piaa0FmodesID = load_fits(fname, "piaa0Fmodescoeff", 1);

            sprintf(fname, "%s/piaaref/piaa1Cmodes.fits", piaacmcsimul_var.piaacmcconfdir);
            piaacmc[0].piaa1CmodesID = load_fits(fname, "piaa1Cmodescoeff", 1);

            sprintf(fname, "%s/piaaref/piaa1Fmodes.fits", piaacmcsimul_var.piaacmcconfdir);
            piaacmc[0].piaa1FmodesID = load_fits(fname, "piaa1Fmodescoeff", 1);

            sprintf(fname, "%s/piaaref/APLCmaskCtransm.txt", piaacmcsimul_var.piaacmcconfdir);
            fp = fopen(fname, "r");
            if(fp!=NULL)
            {
                ret = fscanf(fp, "%f", &tmpf);
                piaacmc[0].fpmaskamptransm = tmpf;
                fclose(fp);
            }
        }
    }



    if((piaacmc[0].piaa0CmodesID==-1)||( piaacmc[0].piaa0FmodesID==-1)||(piaacmc[0].piaa1CmodesID==-1)||( piaacmc[0].piaa1FmodesID==-1))
    {
        sprintf(fname, "%s/apo2Drad.fits", piaacmcsimul_var.piaacmcconfdir);
        if(load_fits(fname, "apo2Drad", 1)==-1)  // CREATE APODIZATION
        {
            sprintf(command, "cp %s/piaaref/apo2Drad.fits %s/apo2Drad.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.piaacmcconfdir);
            ret = system(command);

            sprintf(fname, "%s/apo2Drad.fits", piaacmcsimul_var.piaacmcconfdir);
            if(load_fits(fname, "apo2Drad", 1)==-1)
            {

                printf("Creating 2D apodization for idealized circular monochromatic PIAACMC\n");
                fflush(stdout);

                // first iteration: half size image, 2x zoom, without pupil mask
                IDv1 = create_variable_ID("DFTZFACTOR", 2);
                IDv2 = create_variable_ID("PNBITER", 15);
                coronagraph_make_2Dprolateld(piaacmc[0].fpmaskradld, beamradpix*0.5, piaacmc[0].centObs1, "apotmp1", size/2, "NULLim");

                // expand solution to full size
                basic_resizeim("apotmp1", "apostart", size, size);
                delete_image_ID("apotmp1");

                // full size, 4x zoom
                IDv1 = create_variable_ID("DFTZFACTOR", 4);
                IDv2 = create_variable_ID("PNBITER", 5);
                coronagraph_make_2Dprolateld(piaacmc[0].fpmaskradld, beamradpix, piaacmc[0].centObs1, "apo", size, "pupmaskim");

                // full size, 8x zoom
                chname_image_ID("apo", "apostart");
                IDv1 = create_variable_ID("DFTZFACTOR", 8);
                IDv2 = create_variable_ID("PNBITER", 5);
                coronagraph_make_2Dprolateld(piaacmc[0].fpmaskradld, beamradpix, piaacmc[0].centObs1, "apo", size, "pupmaskim");


                // full size, 16x zoom
                chname_image_ID("apo", "apostart");
                IDv1 = create_variable_ID("DFTZFACTOR", 16);
                IDv2 = create_variable_ID("PNBITER", 10);
                coronagraph_make_2Dprolateld(piaacmc[0].fpmaskradld, beamradpix, piaacmc[0].centObs1, "apo", size, "pupmaskim");


                chname_image_ID("apo", "apo2Drad");
                sprintf(fname, "!%s/apo2Drad.fits", piaacmcsimul_var.piaacmcconfdir);
                save_fits("apo2Drad", fname);

                if(piaacmc[0].PIAAmode == 1)
                {
                    sprintf(fname, "!%s/piaaref/apo2Drad.fits", piaacmcsimul_var.piaacmcconfdir);
                    save_fits("apo2Drad", fname);
                }



                if((piaacmctype==0)&&(loaded==0)) // idealized focal plane mask
                {
                    piaacmc[0].fpmaskamptransm =  -data.variable[variable_ID("APLCmaskCtransm")].value.f;
                    printf("FOCAL PLANE MASK TRANSM = %f\n", piaacmc[0].fpmaskamptransm);
                    printf("Saving default configuration\n");
                    fflush(stdout);
                    saveconf = 1;

                    sprintf(fname, "%s/piaaref/APLCmaskCtransm.txt", piaacmcsimul_var.piaacmcconfdir);
                    fp = fopen(fname, "w");
                    fprintf(fp, "%.20f\n", piaacmc[0].fpmaskamptransm);
                    fclose(fp);
                }

            }
        }
        
        


        // split apodization in conventional pupil apodizer (apoCPA) and PIAA apodization (apo2Drad_PIAA)
        IDapo = image_ID("apo2Drad");
        xsize = data.image[IDapo].md[0].size[0];
        ysize = data.image[IDapo].md[0].size[1];
        IDapo_PIAA = create_2Dimage_ID("apo2Drad_PIAA", xsize, ysize);
        IDapo_CPA = create_2Dimage_ID("apo2Drad_CPA", xsize, ysize);

        if(piaacmc[0].PIAAmode==0)
        {
            for(ii=0; ii<xsize*ysize; ii++) // everything goes to the conventional apodizer
            {
                data.image[IDapo_PIAA].array.F[ii] = 1.0;
                data.image[IDapo_CPA].array.F[ii] = data.image[IDapo].array.F[ii];
            }
        }
        else
        {
            for(ii=0; ii<xsize*ysize; ii++)
            {
                coeff = piaacmc[0].PIAAcoeff; // fraction of apodization done by PIAA - between 0 and 1
                data.image[IDapo_PIAA].array.F[ii] = pow(data.image[IDapo].array.F[ii], coeff);
                data.image[IDapo_CPA].array.F[ii] = pow(data.image[IDapo].array.F[ii], 1.0-coeff);
            }
        }

        copy_image_ID("apo2Drad_CPA", "prePIAA0mask", 0);
        save_fits("prePIAA0mask", "!prePIAA0mask.fits");
        
        
        


        // load PIAA apodization profile and fit it a series of cosines
        PIAACMCsimul_load2DRadialApodization("apo2Drad_PIAA", beamradpix, "outApofit");

        // compute radial PIAA sag -> <piaacmcconfdir>/PIAA_Mshapes.txt
        PIAACMCsimul_init_geomPIAA_rad("outApofit");


        // make 2D sag maps
        sprintf(fname, "%s/PIAA_Mshapes.txt", piaacmcsimul_var.piaacmcconfdir);
        PIAACMCsimul_mkPIAAMshapes_from_RadSag(fname, "piaam0z", "piaam1z");


        if(piaacmcsimul_var.PIAACMC_save==1)
        {
            sprintf(fname, "!%s/piaam0z.fits", piaacmcsimul_var.piaacmcconfdir);
            save_fits("piaam0z", fname);

            sprintf(fname, "!%s/piaam1z.fits", piaacmcsimul_var.piaacmcconfdir);
            save_fits("piaam1z", fname);
        }

        // crop piaam0z and piaam1z to Cmodes size
        ID0 = image_ID("Cmodes");
        size0 = data.image[ID0].md[0].size[0];
        ID1 = create_2Dimage_ID("piaa0zcrop", size0, size0);
        ID = image_ID("piaam0z");
        for(ii=0; ii<size0; ii++)
            for(jj=0; jj<size0; jj++)
                data.image[ID1].array.F[jj*size0+ii] = data.image[ID].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)];

        ID1 = create_2Dimage_ID("piaa1zcrop", size0, size0);
        ID = image_ID("piaam1z");
        for(ii=0; ii<size0; ii++)
            for(jj=0; jj<size0; jj++)
                data.image[ID1].array.F[jj*size0+ii] = data.image[ID].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)];

        make_disk("maskd", size0, size0, 0.5*size0, 0.5*size0, beamradpix);
        make_2Dgridpix("gridpix", size0, size0, 1, 1, 0, 0);
        arith_image_mult("maskd", "gridpix", "maskfit");

        //sprintf(fname, "!%s/maskfit.fits", piaacmcsimul_var.piaacmcconfdir);
        //save_fits("maskfit", fname);

        printf("--------- FITTING COSINE MODES ---------\n");
        fflush(stdout);

        linopt_imtools_image_fitModes("piaa0zcrop", "Cmodes", "maskfit", 1.0e-6, "piaa0Cmodescoeff", 0);
        sprintf(command, "mv %s/eigenv.dat %s/eigenv_piaa0Cmodes.dat", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.piaacmcconfdir);
        ret = system(command);


        //      sprintf(fname, "!%s/piaa0Cmodescoeff.fits", piaacmcsimul_var.piaacmcconfdir);
        //      save_fits("piaa0Cmodescoeff", fname);

        linopt_imtools_image_fitModes("piaa1zcrop", "Cmodes", "maskfit", 1.0e-6, "piaa1Cmodescoeff", 0);
        sprintf(command, "mv %s/eigenv.dat %s/eigenv_piaa1Cmodes.dat", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.piaacmcconfdir);
        ret = system(command);

        //      sprintf(fname, "!%s/piaa1Cmodescoeff.fits", piaacmcsimul_var.piaacmcconfdir);
        //      save_fits("piaa1Cmodescoeff", fname);

        linopt_imtools_image_construct("Cmodes", "piaa0Cmodescoeff", "piaa0Cz");

        if(piaacmcsimul_var.PIAACMC_save==1)
        {
            sprintf(fname, "!%s/piaa0Cz.fits", piaacmcsimul_var.piaacmcconfdir);
            save_fits("piaa0Cz", fname);
        }

        linopt_imtools_image_construct("Cmodes", "piaa1Cmodescoeff", "piaa1Cz");

        if(piaacmcsimul_var.PIAACMC_save==1)
        {
            sprintf(fname, "!%s/piaa1Cz.fits", piaacmcsimul_var.piaacmcconfdir);
            save_fits("piaa1Cz", fname);
        }


        ID0 = image_ID("piaa0Cz");
        size0 = data.image[ID0].md[0].size[0];
        ID1 = image_ID("piaam0z");
        ID = create_2Dimage_ID("piaa0Cres", size0, size0);
        for(ii=0; ii<size0; ii++)
            for(jj=0; jj<size0; jj++)
                data.image[ID].array.F[jj*size0+ii] = data.image[ID1].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)]-data.image[ID0].array.F[jj*size0+ii];
        if(piaacmcsimul_var.PIAACMC_save==1)
        {
            sprintf(fname, "!%s/piaa0Cres.fits", piaacmcsimul_var.piaacmcconfdir);
            save_fits("piaa0Cres", fname);
        }




        ID0 = image_ID("piaa1Cz");
        size0 = data.image[ID0].md[0].size[0];
        ID1 = image_ID("piaam1z");
        ID = create_2Dimage_ID("piaa1Cres", size0, size0);
        for(ii=0; ii<size0; ii++)
            for(jj=0; jj<size0; jj++)
                data.image[ID].array.F[jj*size0+ii] = data.image[ID1].array.F[(jj+(size-size0)/2)*size+(ii+(size-size0)/2)]-data.image[ID0].array.F[jj*size0+ii];

        if(piaacmcsimul_var.PIAACMC_save==1)
        {   sprintf(fname, "!%s/piaa1Cres.fits", piaacmcsimul_var.piaacmcconfdir);
            save_fits("piaa1Cres", fname);
        }



        printf("--------- FITTING FOURIER MODES ---------\n");
        fflush(stdout);

        linopt_imtools_image_fitModes("piaa0Cres", "Fmodes", "maskfit", 0.01, "piaa0Fmodescoeff", 0);
        sprintf(command, "mv %s/eigenv.dat %s/eigenv_piaa0Fmodes.dat", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.piaacmcconfdir);
        ret = system(command);

        linopt_imtools_image_fitModes("piaa1Cres", "Fmodes", "maskfit", 0.01, "piaa1Fmodescoeff", 0);
        sprintf(command, "mv %s/eigenv.dat %s/eigenv_piaa1Fmodes.dat", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.piaacmcconfdir);
        ret = system(command);


        // save_fits("piaa1Fmodescoeff", "!piaa1Fmodescoeff.fits");

        //linopt_imtools_image_construct("Fmodes", "piaa0Fmodescoeff", "piaa0Fz");
        //   save_fits("piaa0Fz", "!piaa0Fz.fits");
        //arith_image_sub("piaa0Cres", "piaa0Fz", "piaa0CFres");
        //save_fits("piaa0CFres", "!piaa0CFres.fits");
        delete_image_ID("piaa0zcrop");

        //linopt_imtools_image_construct("Fmodes", "piaa1Fmodescoeff", "piaa1Fz");
        //save_fits("piaa1Fz", "!piaa1Fz.fits");
        //arith_image_sub("piaa1Cres", "piaa1Fz", "piaa1CFres");
        //save_fits("piaa1CFres", "!piaa1CFres.fits");
        delete_image_ID("piaa1zcrop");

        delete_image_ID("maskfit");



        piaacmc[0].piaa0CmodesID = image_ID("piaa0Cmodescoeff");
        piaacmc[0].piaa0FmodesID = image_ID("piaa0Fmodescoeff");
        piaacmc[0].piaa1CmodesID = image_ID("piaa1Cmodescoeff");
        piaacmc[0].piaa1FmodesID = image_ID("piaa1Fmodescoeff");



        sprintf(fname, "!%s/piaaref/piaa0Cmodes.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits(data.image[piaacmc[0].piaa0CmodesID].name, fname);

        sprintf(fname, "!%s/piaaref/piaa0Fmodes.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits(data.image[piaacmc[0].piaa0FmodesID].name, fname);

        sprintf(fname, "!%s/piaaref/piaa1Cmodes.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits(data.image[piaacmc[0].piaa1CmodesID].name, fname);

        sprintf(fname, "!%s/piaaref/piaa1Fmodes.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits(data.image[piaacmc[0].piaa1FmodesID].name, fname);

        sprintf(command, "cp %s/piaaref/* %s/", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.piaacmcconfdir);
        ret = system(command);
    }


    // ============ MAKE FOCAL PLANE MASK ===============

	PIAACMCsimul_update_fnamedescr();
	
	sprintf(command, "echo \"%s\" > fpm_name.txt", piaacmcsimul_var.fnamedescr);
	ret = system(command);
	
	sprintf(command, "echo \"%s\" > fpm_name_conf.txt", piaacmcsimul_var.fnamedescr_conf);
	ret = system(command);
	
	


    /*    piaacmcsimul_var.CREATE_fpmzmap = 0;
        if( piaacmcsimul_var.FORCE_CREATE_fpmzmap == 0 )
        {
            if(image_ID("fpmzmap")==-1)
                {
                    sprintf(fname, "%s/fpmzmap%d_%03ld_%03ld.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.PIAACMC_FPMsectors, piaacmc[0].NBrings, piaacmc[0].focmNBzone);
                    load_fits(fname, "fpmzmap", 1);
                    if(image_ID("fpmzmap")==-1)
                        piaacmcsimul_var.CREATE_fpmzmap = 1;
                }
        }
        else
            piaacmcsimul_var.CREATE_fpmzmap = 1;

        if( piaacmcsimul_var.CREATE_fpmzmap == 1 )
        {
            if(image_ID("fpmzmap")!=-1)
                delete_image_ID("fpmzmap");
            PIAACMCsimul_mkFPM_zonemap("fpmzmap");
            sprintf(fname, "!%s/fpmzmap%d_%03ld_%03ld.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.PIAACMC_FPMsectors, piaacmc[0].NBrings, piaacmc[0].focmNBzone);
            save_fits("fpmzmap", fname);
        }
    */

    if(image_ID("fpmzmap")==-1)
    {
        printf("Make zonemap ...\n");
        fflush(stdout);
        //debug from Justin
        //sleep(10);
        PIAACMCsimul_mkFPM_zonemap("fpmzmap");
    }
    else
    {
        printf("zonemap already exists\n");
        fflush(stdout);
    }

    //    sprintf(fname, "!%s/fpmzmap.fits", piaacmcsimul_var.piaacmcconfdir);
    //    save_fits("fpmzmap", fname);



    /* sprintf(fname, "!%s/fpmzmap.fits", piaacmcsimul_var.piaacmcconfdir);
     save_fits("fpmzmap", fname);
     exit(0);*/


    // zones thickness

    piaacmcsimul_var.CREATE_fpmzt = 0;
    if( piaacmcsimul_var.FORCE_CREATE_fpmzt == 0 )
    {
        piaacmc[0].zonezID = image_ID("fpmzt");
        if(piaacmc[0].zonezID == -1)
        {   
			PIAACMCsimul_update_fnamedescr();         
            sprintf(fname, "%s/fpm_zonez.%s.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.fnamedescr);

            printf("LOADING FILE NAME : \"%s\"  -  %ld %d \n", fname, piaacmctype, loaded);
			fflush(stdout);
			// debug from Justin
			sleep(10);
			
            piaacmc[0].zonezID = load_fits(fname, "fpmzt", 1);
            if(piaacmc[0].zonezID == -1)
                piaacmcsimul_var.CREATE_fpmzt = 1;
        }
    }
    else
        piaacmcsimul_var.CREATE_fpmzt = 1;



    if( piaacmcsimul_var.CREATE_fpmzt == 1 )
    {
        printf("Creating fpmzt, saving as fpm_zonez.fits - %ld %d\n", piaacmctype, loaded);
        fflush(stdout);
        piaacmc[0].zonezID = image_ID("fpmzt");
        if(piaacmc[0].zonezID!=-1)
            delete_image_ID("fpmzt");

        piaacmc[0].zonezID = create_2Dimage_ID_double("fpmzt", piaacmc[0].focmNBzone, 1);
        t = 1.0e-9;

        if(piaacmctype==0) // idealized focal plane mask
        {
            printf("IDEALIZED FOCAL PLANE MASK\n");
            fflush(stdout);
            // debug from Justin
            //sleep(10);
            //          exit(0);

            // measure dpha/dt
            t0 = 1.0e-8;
            pha0 = OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial_code, t0, 0.5*(piaacmcsimul_var.LAMBDASTART + piaacmcsimul_var.LAMBDAEND));
            // set t to get PI phase
            t = (M_PI/pha0)*t0;
            printf("t = %g m (%lf %g) -> %g %g\n", t, pha0, t0, OPTICSMATERIALS_pha_lambda(piaacmc[0].fpmmaterial_code, t, 0.5*(piaacmcsimul_var.LAMBDASTART + piaacmcsimul_var.LAMBDAEND)), 0.5*(piaacmcsimul_var.LAMBDASTART + piaacmcsimul_var.LAMBDAEND));
			fflush(stdout);
			// debug from Justin
			//sleep(10);
            printf(" -- lambda = %g\n", piaacmc[0].lambda);
            printf(" -- lambdaB = %g\n", piaacmc[0].lambdaB);
            printf(" -- LAMBDASTART = %g\n", piaacmcsimul_var.LAMBDASTART);
            printf(" -- LAMBDAEND = %g\n", piaacmcsimul_var.LAMBDAEND);
        }
        else
        {
            printf("CREATING EXAMPLE FOCAL PLANE MASK  %ld %d\n", piaacmctype, loaded);
            fflush(stdout);
            //                exit(0);
        }

        for(ii=0; ii<piaacmc[0].focmNBzone; ii++)
            data.image[piaacmc[0].zonezID].array.D[ii] = t;

		PIAACMCsimul_update_fnamedescr();
		sprintf(fname, "!%s/fpm_zonez.%s.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.fnamedescr);

        printf("Writing %s\n", fname);
        fflush(stdout);
        // debug from Justin
        sleep(10);
        save_fits("fpmzt", fname);
    }


    // zones transmission amplitude

    printf("CREATE_fpmza = %d\n", piaacmcsimul_var.CREATE_fpmza);
    fflush(stdout);
    if( piaacmcsimul_var.FORCE_CREATE_fpmza == 0 )
    {
        piaacmc[0].zoneaID = image_ID("fpmza");
        if(piaacmc[0].zoneaID == -1)
        {
			PIAACMCsimul_update_fnamedescr();
            sprintf(fname, "%s/fpm_zonea.%s.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.fnamedescr);

            printf("LOADING FILE NAME : \"%s\"\n", fname);
            piaacmc[0].zoneaID = load_fits(fname, "fpmza", 1);

            if(piaacmc[0].zoneaID == -1)
                piaacmcsimul_var.CREATE_fpmza = 1;
        }
    }
    else
        piaacmcsimul_var.CREATE_fpmza = 1;



    if( piaacmcsimul_var.CREATE_fpmza == 1 )
    {
        if(piaacmc[0].zoneaID != -1)
            delete_image_ID("fpmza");
        piaacmc[0].zoneaID = create_2Dimage_ID_double("fpmza", piaacmc[0].focmNBzone, 1);

        if( piaacmcsimul_var.PIAACMC_MASKRADLD > 0.2 ) // physical mask
        {
            printf("PHYSICAL MASK ... %ld zones\n", piaacmc[0].focmNBzone);
            fflush(stdout);
            
            for(ii=0; ii<piaacmc[0].focmNBzone; ii++)
                data.image[piaacmc[0].zoneaID].array.D[ii] = 1.0;
        }
        else // idealized mask
        {
            printf("IDEALIZED MASK ... %ld zones\n", piaacmc[0].focmNBzone);
            for(ii=0; ii<piaacmc[0].focmNBzone; ii++)
                data.image[piaacmc[0].zoneaID].array.D[ii] = piaacmc[0].fpmaskamptransm;
        }

		PIAACMCsimul_update_fnamedescr();
        sprintf(fname, "!%s/fpm_zonea.%s.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.fnamedescr);

        printf("Writing %s\n", fname);
        save_fits("fpmza", fname);
    }

    //   printf("%d piaacmc[0].fpmaskamptransm = %f       %lf\n", piaacmcsimul_var.CREATE_fpmza, piaacmc[0].fpmaskamptransm, data.image[piaacmc[0].zoneaID].array.D[0]);
    //   sleep(10);


    // ============= MAKE LYOT STOPS =======================
    printf("LOADING/CREATING LYOT MASK  - %ld masks  (PIAAmode = %d, %ld x %ld)\n", piaacmc[0].NBLyotStop, piaacmc[0].PIAAmode, xsize, ysize);
   // list_image_ID();
    size2 = size*size;

	
    if(piaacmc[0].PIAAmode == 1)
    {
        for(i=0; i<piaacmc[0].NBLyotStop; i++)
        {
            printf("LYOT MASK %ld\n", i);
            fflush(stdout);

            sprintf(fname, "%s/LyotStop%ld.fits", piaacmcsimul_var.piaacmcconfdir, i);
            sprintf(name, "lyotstop%ld", i);

            piaacmc[0].IDLyotStop[i] = image_ID(name);
            if(piaacmc[0].IDLyotStop[i]==-1)
            {
                sprintf(fname, "!%s/LyotStop%ld.fits", piaacmcsimul_var.piaacmcconfdir, i);
                if ((i == 0) && (piaacmc[0].NBLyotStop > 1))
                    piaacmc[0].IDLyotStop[i] = PIAACMCsimul_mkSimpleLyotStop(name, -0.01, 0.98);
                else if ((i == 1) && (piaacmc[0].NBLyotStop > 2))
                    piaacmc[0].IDLyotStop[i] = PIAACMCsimul_mkSimpleLyotStop(name, piaacmc[0].centObs1+0.02, 1.2);
                else
                    piaacmc[0].IDLyotStop[i] = PIAACMCsimul_mkSimpleLyotStop(name, piaacmc[0].centObs1+0.02, 0.98);

                save_fl_fits(name, fname);
            }
        }
    }
    else
    {
        i = 0;
        sprintf(fname, "!%s/LyotStop%ld.fits", piaacmcsimul_var.piaacmcconfdir, i);
        sprintf(name, "lyotstop%ld", i);


        
        piaacmc[0].IDLyotStop[i] = image_ID(name);
        if(piaacmc[0].IDLyotStop[i]==-1)
			{
				piaacmc[0].IDLyotStop[i] = create_2Dimage_ID(name, xsize, ysize);
				ID = image_ID("pupmaskim");
				for(ii=0; ii<xsize*ysize; ii++)
					if(data.image[ID].array.F[ii] < 0.99999999)
						data.image[piaacmc[0].IDLyotStop[i]].array.F[ii] = 0.0;
					else
						data.image[piaacmc[0].IDLyotStop[i]].array.F[ii] = 1.0;
						
				for(ii=0;ii<xsize;ii++)
					for(jj=0;jj<ysize;jj++)
						{
							x = 1.0*ii-0.5*xsize;
							y = 1.0*jj-0.5*ysize;
							rad = sqrt(x*x+y*y);
							rad /= beamradpix;
							if(rad<(piaacmc[0].centObs1+0.5/beamradpix))
								data.image[piaacmc[0].IDLyotStop[i]].array.F[jj*xsize+ii] = 0.0;
							if(rad>(1.0-0.5/beamradpix))
								data.image[piaacmc[0].IDLyotStop[i]].array.F[jj*xsize+ii] = 0.0;
						}
				save_fl_fits(name, fname);
			}
    }



    if(saveconf==1)
        PIAACMCsimul_savepiaacmcconf( piaacmcsimul_var.piaacmcconfdir);


    return(0);
}

