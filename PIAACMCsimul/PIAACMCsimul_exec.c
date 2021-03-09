/**
 * @file    PIAAACMCsimul_init.c
 * @brief   PIAA-type coronagraph design, execute
 * 
 * Can design both APLCMC and PIAACMC coronagraphs
 *  
 * @author  O. Guyon
 * @date    26 nov 2017
 *
 *
 * @bug No known bugs.
 * 
 */


// System includes
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <sys/stat.h>
#include <time.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>



// milk includes
#include "CommandLineInterface/CLIcore.h"
#include "00CORE/00CORE.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "image_gen/image_gen.h"
#include "statistic/statistic.h"

#include "linopt_imtools/linopt_imtools.h"
#include "OpticsMaterials/OpticsMaterials.h"
#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"



extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTSYST *optsyst;

extern OPTPIAACMCDESIGN *piaacmc;











double PIAACMCsimul_regularization_PIAAshapes_value()
{
    double value = 0.0;
    long IDref;
    long ID;
    long jj;
    long ID_CPAfreq;

    // add regularization component of the evaluation metrix
    // first, compute and add PIAA shape regularization value (value) if applicable
    // this is the same as the previous computation of val0

    // index of the PIAA element 0 (first mirror) shapes via cosine modes
    ID = piaacmc[0].piaa0CmodesID;

    // index of PIAA shapes reference image
    IDref = image_ID("piaa0Cmref");
    if(IDref==-1)
    {
        // error message if we get here?  ***************************************
        // if the reference image doesn't exist, create it
        IDref = create_2Dimage_ID("piaa0Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);

        // initialize to zero shape
        for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
            data.image[IDref].array.F[jj] = 0.0;
    }


    // This section of code does not actually fill in the regularization terms in the output vector
    // filling in is done later.  Here we are only computing the initial reference scalar objective value

    // For each cosine mode set the optimization parameter = cosine mode modified by regularization
    for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
    {
        // compute square of C*(deviation from reference)*(mode index)^(alpha)
        // so higher-index zones (higher spatial frequency) are more
        // heavily penalized

        double tmp;
        tmp = piaacmcsimul_var.linopt_piaa0C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj, piaacmcsimul_var.linopt_piaa0C_regcoeff_alpha);
        value += tmp*tmp;
    }


    // do the same for PIAA element 1
    ID = piaacmc[0].piaa1CmodesID;
    IDref = image_ID("piaa1Cmref");
    if(IDref==-1)
    {
        IDref = create_2Dimage_ID("piaa1Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
        for(jj=0; jj<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; jj++)
            data.image[IDref].array.F[jj] = 0.0;
    }
    for(jj=0; jj<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0]; jj++)
    {
        double tmp;
        tmp = piaacmcsimul_var.linopt_piaa1C_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*jj, piaacmcsimul_var.linopt_piaa1C_regcoeff_alpha);
        value += tmp*tmp;
    }


    // get spatial frequency of each mode in cycles/aperture
    ID_CPAfreq = image_ID("cpamodesfreq");

    // do the same for PIAA element 0 and 1 for the Fourier modes
    // this time use the actual spatial frequency rather than mode index as proxy for frequency
    ID = piaacmc[0].piaa0FmodesID;
    IDref = image_ID("piaa0Fmref");
    if(IDref==-1)
    {
        IDref = create_2Dimage_ID("piaa0Fmref", data.image[piaacmc[0].piaa0FmodesID].md[0].size[0], 1);
        for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
            data.image[IDref].array.F[jj] = 0.0;
    }
    for(jj=0; jj<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; jj++)
    {
        double tmp;
        tmp = piaacmcsimul_var.linopt_piaa0F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaacmcsimul_var.linopt_piaa0F_regcoeff_alpha);
        value += tmp*tmp;
    }

    ID = piaacmc[0].piaa1FmodesID;
    IDref = image_ID("piaa1Fmref");
    if(IDref==-1)
    {
        IDref = create_2Dimage_ID("piaa1Fmref", data.image[piaacmc[0].piaa1FmodesID].md[0].size[0], 1);
        for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
            data.image[IDref].array.F[jj] = 0.0;
    }
    for(jj=0; jj<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; jj++)
    {
        double tmp;
        tmp = piaacmcsimul_var.linopt_piaa1F_regcoeff * (data.image[ID].array.F[jj]-data.image[IDref].array.F[jj]) * pow(1.0*data.image[ID_CPAfreq].array.F[jj], piaacmcsimul_var.linopt_piaa1F_regcoeff_alpha);
        value += tmp*tmp;
    }


    return value;

}




/** 
 * @brief Compute regularization term for focal plane mask sag values
 * 
 */ 

double PIAACMCsimul_regularization_fpmsag_value()
{
    double regvalue;  // output regularization value
    long IDzonez;     // image: zone sag values

    // regvalue is the regularization value for the focal plane mask sag values
    // same as above, but not dependent on position
    regvalue = 0.0;
    IDzonez = piaacmc[0].zonezID;

    long zoneindex;
    for(zoneindex=0; zoneindex < data.image[IDzonez].md[0].size[0]; zoneindex++)
    {
        // compute the square of (sag/coeff)^alpha
        double tmp;
        tmp = pow(data.image[IDzonez].array.D[zoneindex]/piaacmc[0].fpmsagreg_coeff, piaacmc[0].fpmsagreg_alpha);
        regvalue += tmp*tmp;
    }
    
    return regvalue;
}









long PIAACMCsimul_regularization_PIAAshapes_add1Dvector(long ID1D, long index0)
{
    long vindex;  // index in output vector
    long ii;     // running (temporary) index

    long IDref;
    long ID_CPAfreq;
    long ID;


    // fill in the shape regularization value's response to piaacmcsimul_var.linopt_paramdelta using the
    // same formulas as before with the delta PSF as input
    ID = piaacmc[0].piaa0CmodesID;
    vindex = index0;

    IDref = image_ID("piaa0Cmref");
    if(IDref==-1)
    {
        IDref = create_2Dimage_ID("piaa0Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
        long ii;
        for(ii=0; ii<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; ii++)
            data.image[IDref].array.F[ii] = 0.0;
    }

    for(ii=0; ii<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; ii++)
    {
        data.image[ID1D].array.F[vindex] = piaacmcsimul_var.linopt_piaa0C_regcoeff * (data.image[ID].array.F[ii]-data.image[IDref].array.F[ii]) * pow(1.0*ii, piaacmcsimul_var.linopt_piaa0C_regcoeff_alpha);
        vindex++;
    }

    ID = piaacmc[0].piaa1CmodesID;
    IDref = image_ID("piaa1Cmref");
    if(IDref==-1)
    {
        IDref = create_2Dimage_ID("piaa1Cmref", data.image[piaacmc[0].piaa0CmodesID].md[0].size[0], 1);
        for(ii=0; ii<data.image[piaacmc[0].piaa0CmodesID].md[0].size[0]; ii++)
            data.image[IDref].array.F[ii] = 0.0;
    }
    for(ii=0; ii<data.image[piaacmc[0].piaa1CmodesID].md[0].size[0]; ii++)
    {
        data.image[ID1D].array.F[vindex] = piaacmcsimul_var.linopt_piaa1C_regcoeff * (data.image[ID].array.F[ii]-data.image[IDref].array.F[ii]) * pow(1.0*ii, piaacmcsimul_var.linopt_piaa1C_regcoeff_alpha);
        vindex++;
    }

    ID_CPAfreq = image_ID("cpamodesfreq");

    ID = piaacmc[0].piaa0FmodesID;
    IDref = image_ID("piaa0Fmref");
    if(IDref==-1)
    {
        IDref = create_2Dimage_ID("piaa0Fmref", data.image[piaacmc[0].piaa0FmodesID].md[0].size[0], 1);
        for(ii=0; ii<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; ii++)
            data.image[IDref].array.F[ii] = 0.0;
    }
    for(ii=0; ii<data.image[piaacmc[0].piaa0FmodesID].md[0].size[0]; ii++)
    {
        data.image[ID1D].array.F[vindex] = piaacmcsimul_var.linopt_piaa0F_regcoeff * (data.image[ID].array.F[ii]-data.image[IDref].array.F[ii]) * pow(1.0*data.image[ID_CPAfreq].array.F[ii], piaacmcsimul_var.linopt_piaa0F_regcoeff_alpha);
        vindex++;
    }

    ID = piaacmc[0].piaa1FmodesID;
    IDref = image_ID("piaa1Fmref");
    if(IDref==-1)
    {
        IDref = create_2Dimage_ID("piaa1Fmref", data.image[piaacmc[0].piaa1FmodesID].md[0].size[0], 1);
        for(ii=0; ii<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; ii++)
            data.image[IDref].array.F[ii] = 0.0;
    }
    for(ii=0; ii<data.image[piaacmc[0].piaa1FmodesID].md[0].size[0]; ii++)
    {
        data.image[ID1D].array.F[vindex] = piaacmcsimul_var.linopt_piaa1F_regcoeff * (data.image[ID].array.F[ii]-data.image[IDref].array.F[ii]) * pow(1.0*data.image[ID_CPAfreq].array.F[ii], piaacmcsimul_var.linopt_piaa1F_regcoeff_alpha);
        vindex++;
    }

    return vindex;
}

















/**
 *
 * @brief Main simulation routine
 *
 *
 * @param[in] confindex  PIAACMC configuration index pointing to the input/output directory number
 * @param[in] mode       Type of operation to be performed
 *
 */

int PIAACMCsimul_exec(
	const char *confindex, 
	long mode
	)
{
    FILE *fp;
    char command[1000];
    char dirname[500];

    double valref, val0;
    long i, jj, ii;

    long IDv, ID, IDref;
    long IDmodes, IDmodes2D;
    long xsize = 0;
    long ysize = 0;
    long k;

    long iter;

    char fname[800];
    long IDm, ID1D, ID1Dref;
    long size1Dvec, size1Dvec0;


    int r;
    double val;

    double fpmradld = 0.95;
    double centobs0 = 0.3;
    double centobs1 = 0.2;

    long mz;

    long ii1, ii2, ki, kv, ki1;

    double scangain;
    double scanstepgain = 0.001;
    int linscanOK;
    double valold, oldval;
    double bestgain;
    double linoptgainarray[100];
    double linoptvalarray[100];
    int linoptlimflagarray[100];
    long NBlinoptgain;
    long kmax;


    double valtest;
    int ret;
    int kmaxC, kmaxF;

    long st;
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2rand;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function fpmeval_func;
    int status;

    // simulated annealing loop(s)
    long NBITER_SA = 1000000;
    long ITERMAX;
    long iter1;
    double SAcoeff;
    double amp1;
    double val1;
    int SA_MV = 0;
    double tmp1;
    double bestvalue;
    FILE *fpbest;

    // downhill simplex
    long NBITER_DS = 1000;
    long iterMax;
    int OK;
    int KEEP;
    double KEEPlimit = 1.0;
    double KEEPlimit1 = 1.0;
    double eps1;
    double mmsize;

    long elem;


    int iterOK;
    char stopfile[500];

    double alphareg;
    double bestval = 0.0;
    long IDoptvec = -1;

    double acoeff0, acoeff1, acoeff2;
    float scangainfact;
    int alphascaninit = 0;

    uint32_t *sizearray;
    long IDstatus;

    double valContrast;
    double tmp;
    int initbestval = 0;

    long ID_CPAfreq;





#ifdef PIAASIMUL_LOGFUNC0
    sprintf(flogcomment, "%s %ld", confindex, mode);
    PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, flogcomment);
#endif




    piaacmcsimul_var.linopt_REGPIAASHAPES = 0;

    piaacmcsimul_var.linopt_piaa0C_regcoeff = 0.0e-7;
    piaacmcsimul_var.linopt_piaa1C_regcoeff = 0.0e-7;
    piaacmcsimul_var.linopt_piaa0C_regcoeff_alpha = 1.0;
    piaacmcsimul_var.linopt_piaa1C_regcoeff_alpha = 1.0;

    piaacmcsimul_var.linopt_piaa0F_regcoeff = 0.0e-7;
    piaacmcsimul_var.linopt_piaa1F_regcoeff = 0.0e-7;
    piaacmcsimul_var.linopt_piaa0F_regcoeff_alpha = 1.0;
    piaacmcsimul_var.linopt_piaa1F_regcoeff_alpha = 1.0;


    piaacmcsimul_var.linopt_REGFPMSAG = 0;





    // Create status shared variable
    // this allows realtime monitoring of the code by other processes
    // sets status at different points in the code
    IDstatus = image_ID("stat_PIAACMCsimulexec");
    if(IDstatus == -1)
        IDstatus = read_sharedmem_image("stat_PIAACMCsimulexec");

    sizearray = (uint32_t*) malloc(sizeof(uint32_t)*2);
    sizearray[0] = 1;
    sizearray[1] = 1;
    IDstatus = create_image_ID("stat_PIAACMCsimulexec", 2, sizearray, _DATATYPE_UINT16, 1, 0);
    free(sizearray);



    //piaacmc = NULL; // set the pointer to the piaacmc structure to null

    // if the optical system pointer is empty, create an empty version
    if(optsyst==NULL)
    {
        optsyst = (OPTSYST*) malloc(sizeof(OPTSYST));
        optsyst[0].SAVE = 1;
    }

    for(elem=0; elem<100; elem++)
        optsyst[0].keepMem[elem] = 0; // flag that says save this element for reuse


    // set the result directories
    sprintf( piaacmcsimul_var.piaacmcconfdir, "%s", confindex);
    sprintf( data.SAVEDIR, "%s", piaacmcsimul_var.piaacmcconfdir);


    piaacmcsimul_var.linopt_NBiter = 1000;
    piaacmcsimul_var.linopt_number_param = 0;



    optsyst[0].SAVE = piaacmcsimul_var.PIAACMC_save;

    // get variables from command line, possibly sets globals
    if( (IDv = variable_ID("PIAACMC_centobs0")) != -1)
        centobs0 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_centobs1")) != -1)
        centobs1 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }

    // start a log of code mode entry/exit times
    sprintf(command, "echo \"%03ld     $(date)\" >> ./log/PIAACMC_mode_log.txt", mode);
    ret = system(command);
    printf("command = %s\n", command);



    // set the name of the stopfile
    sprintf(stopfile, "%s/stopmode%ld.txt", piaacmcsimul_var.piaacmcconfdir, mode);

    switch (mode) {

    case 0 :
        PIAACMCsimul_exec_compute_image();
        printf("EXEC CASE 0 COMPLETED\n");
        fflush(stdout);
        break;

    case 1 :
        PIAACMCsimul_exec_optimize_lyot_stop_position();
        break;

    case 2 :
        PIAACMCsimul_exec_optimize_fpmtransmission();
        break;

    case 3 :
        val = PIAACMCsimul_exec_computePSF_no_fpm();
        break;

    case 4 :
        PIAACMCsimul_exec_optimize_PIAA_shapes();
        break;

    case 5 :
        PIAACMCsimul_exec_optimize_lyot_stops_shapes_positions();
        break;

    case 11 :
        PIAACMCsimul_exec_multizone_fpm_calib();
        break;

    case 13 :
        PIAACMCsimul_exec_optimize_fpm_zones();
        break;

    case 40 :
        PIAACMCsimul_exec_optimize_PIAA_shapes_fpmtransm();
        break;

    case 100 : // evaluate current design: polychromatic contrast, pointing sensitivity
        PIAACMCsimul_eval_poly_design();
        break;

    case 101 :
        PIAACMCsimul_measure_transm_curve();
        break;




    case 302 : // restore configuration settings
        printf("=================================== mode 302 ===================================\n");

        sprintf(command, "cp %s/saveconf/conf_*.txt %s/", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.piaacmcconfdir);
        r = system(command);
        break;


    default :
        printERROR(__FILE__,__func__,__LINE__, "mode not recognized");
        break;
    }



































    // linear optimization set up in modes 13 and 40
    //
    // the output parameters start as the evaluation zone values in "imvect"
    // If we are regularizing, we supplement the output parameters by adding
    // penalty terms as additional output parameters
    if(piaacmcsimul_var.LINOPT == 1) // linear optimization
    {
        // for state tracking and statistics
        data.image[IDstatus].array.UI16[0] = 5;

        // Compute Reference on-axis performance contrast (valref)
        PIAACMCsimul_makePIAAshapes(piaacmc, 0);
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm
        valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, 0);


        // save the current configuration to the _linopt directory
        sprintf(dirname, "%s_linopt", piaacmcsimul_var.piaacmcconfdir);
        PIAACMCsimul_savepiaacmcconf(dirname);

        // import configuration from _linopt directory
        sprintf(command, "rsync -au --progress %s/* ./%s/", dirname, piaacmcsimul_var.piaacmcconfdir);
        r = system(command);

        // "cp -n" will only copy the file if destination does not exist
        // this will ensure that the .ref.fits files are the PIAA shapes before any linear optimization
        // these will be the "reference" shapes used to regularize PIAA shapes: the deviation from this reference will be kept small by the regularization

        sprintf(command, "cp -n %s/piaa0Cmodes.fits %s/piaa0Cmodes.ref.fits", dirname, dirname);
        r = system(command);
        sprintf(command, "cp -n %s/piaa0Fmodes.fits %s/piaa0Fmodes.ref.fits", dirname, dirname);
        r = system(command);

        sprintf(command, "cp -n %s/piaa1Cmodes.fits %s/piaa1Cmodes.ref.fits", dirname, dirname);
        r = system(command);
        sprintf(command, "cp -n %s/piaa1Fmodes.fits %s/piaa1Fmodes.ref.fits", dirname, dirname);
        r = system(command);

        sprintf(command, "cp -n %s/piaacmcparams.conf %s/piaacmcparams.ref.conf", dirname, dirname);
        r = system(command);


        // load PIAA reference shapes (used for regularization)
        sprintf(fname, "%s/piaa0Cmodes.ref.fits", dirname);
        load_fits(fname, "piaa0Cmref", 1);

        sprintf(fname, "%s/piaa1Cmodes.ref.fits", dirname);
        load_fits(fname, "piaa1Cmref", 1);

        sprintf(fname, "%s/piaa0Fmodes.ref.fits", dirname);
        load_fits(fname, "piaa0Fmref", 1);

        sprintf(fname, "%s/piaa1Fmodes.ref.fits", dirname);
        load_fits(fname, "piaa1Fmref", 1);

        // we have now saved the starting point of the optimization for future comparison
        // in the <piaacmcconfdir>_linopt directory

        // here we compute regularization value of piaashapes and store it in the val0 variable
        // if regularization is of PIAA shapes is ON, then val0 will be computed and added to the overal performance metric valref

        // regularize the piaashapes via a penalty added to the reference contrast valref
        // The optimization minimizes the summed contrast + val0 + val1.
        // Regularization is via adding a constant val0 + val1 to the contrast we're minimizing
        // note that here we're setting output parameters.
        if(piaacmcsimul_var.linopt_REGPIAASHAPES==1)
        {
            val0 = PIAACMCsimul_regularization_PIAAshapes_value();
            printf("VALREF = %g + %g -> %g\n", valref, val0, valref+val0);
            valref += val0;
        }


        // val1 is the regularization value for the focal plane mask sag values
        // same as above, but not dependent on position
        val1 = 0.0;
        if(piaacmcsimul_var.linopt_REGFPMSAG == 1)
        {
            val1 = PIAACMCsimul_regularization_fpmsag_value();
            valref += val1;
        }
        // At this point all we've done is compute the overall performance metric including
        // regularization in valref.


        // for state tracking and statistics
        data.image[IDstatus].array.UI16[0] = 6;
        printf("================================ Reference = %g\n", valref);

        // copy imvect to vecDHref "vector dark hole reference"
        // vecDHref is the dark hole complex amplitude state at the beginning of the linear optimization
        // (this dark hole is nominally a full annulus from 1.5 to ~8 lambda/D, created at the top
        // of PIAACMCsimul_computePSF with size controlled by scoringIWA and scoringOWA
        // the corresponding performance metric is valref
        chname_image_ID("imvect", "vecDHref"); // note: imvect was computed by PIAACMCsimul_computePSF called ~150 lines above
        ID = image_ID("vecDHref"); // ID changed identity
        xsize = data.image[ID].md[0].size[0];
        ysize = data.image[ID].md[0].size[1];



        // for state tracking and statistics
        data.image[IDstatus].array.UI16[0] = 7;

        // now we will just determine the size of the size of the
        // optimization vectors that we will actually fill in later
        // save vecDHref initial state as a reference
        sprintf(fname, "!%s/vecDMref.fits", piaacmcsimul_var.piaacmcconfdir);
        save_fits("vecDHref", fname);
        // get and copy size of vecDHref, 'cause we're manipulating size1Dvec
        size1Dvec = data.image[ID].md[0].nelement;
        size1Dvec0 = size1Dvec;

        // PIAA shapes regularization
        // if regularization is turned on, the size of the evaluation vector is increased to include PIAA shape coefficients, in addition to complex amplitude in focal plane
        // the optimization code will then simultaneously minimize the light in the focal plane AND the PIAA shape deviations from the nominal shape
        if(piaacmcsimul_var.linopt_REGPIAASHAPES==1)
        {
            // there are 4 groups of PIAA shape parameters: 2 sets of cosine modes and 2 sets of Fourier modes
            size1Dvec += data.image[piaacmc[0].piaa0CmodesID].md[0].size[0];
            size1Dvec += data.image[piaacmc[0].piaa1CmodesID].md[0].size[0];
            size1Dvec += data.image[piaacmc[0].piaa0FmodesID].md[0].size[0];
            size1Dvec += data.image[piaacmc[0].piaa1FmodesID].md[0].size[0];
        }

        // The same approach is used for regularization of the focal plane mask sag values
        // the sag values are appended to the evaluation vector
        if(piaacmcsimul_var.linopt_REGFPMSAG==1)
        {
            size1Dvec += data.image[piaacmc[0].zonezID].md[0].size[0];
        }


        // re-package vector into 1D array and add regularization terms
        // the resulting vector is stored in vecDHref1D
        // we also create a mask image which is intended to mask out some of the pixels from the evaluation
        // the mask is currently not used, so we will write 1.0 in all of its pixels
        // DHmask and vecDHref1D contain all the optimization parameters as set above
        IDm = create_2Dimage_ID("DHmask", size1Dvec, 1); // "ID of evaluation mode mask"
        ID1Dref = create_2Dimage_ID("vecDHref1D", size1Dvec, 1); // "ID of 1D dark zone reference"

        // we first write 1.0 into the focal plane complex amplitudes in the vector
        ID = image_ID("vecDHref");
        for(ii=0; ii<data.image[ID].md[0].nelement; ii++)
        {
            // imbed vecDHref into the evaluation zone part of of the full parameter vector vecDHref1D
            data.image[ID1Dref].array.F[ii] = data.image[ID].array.F[ii];
            // sets the evaluation zone part of of the full parameter vector vecDHref1D to 1
            // 1 means the evaluation zone is on
            data.image[IDm].array.F[ii] = 1.0;
        }
        // !!!!!! WARNING !!!!!!!
        // the state of ii at this point drives the code below and will evolve until the comment
        // that says we're done with ii

        // for state tracking and statistics
        data.image[IDstatus].array.UI16[0] = 8;

        // Now actually fill in the regularized output vector.
        // If we are not regularizing, the output evaluation zone values are filled in by the
        // PSF calls in the optimization loop
        // and then append regularization vectors in the main evaluation vector
        // Note: index ii is incremented as we add "sub-vectors" into the main evaluation vector
        if(piaacmcsimul_var.linopt_REGPIAASHAPES == 1)
        {
            // initialize by filling in the regularization terms of the output,
            // !! starting at the current value of ii !!
            // This means we're actually writing the output vector regularization terms.
            // Otherwise this is the same as the if(REGPIAASHAPES == 1) code block above
            ii = PIAACMCsimul_regularization_PIAAshapes_add1Dvector(ID1Dref, ii);
        }




        // same for the sags, starting at the current value of ii
        if(piaacmcsimul_var.linopt_REGFPMSAG == 1)
        {
            ID = piaacmc[0].zonezID;
            for(jj=0; jj < data.image[ID].md[0].size[0]; jj++)
            {
                data.image[ID1Dref].array.F[ii] = pow( data.image[ID].array.D[jj]/piaacmc[0].fpmsagreg_coeff, piaacmc[0].fpmsagreg_alpha);
                data.image[IDm].array.F[ii] = 1.0;
                ii++;
            }
        }
        // !!!!!!
        // we're done with ii


        delete_image_ID("vecDHref"); // vecDHref has beem embedded into vecDHref1D

        // at this point, we have completed the initialization, and the optimization loop starts


        initbestval = 0;
        // file that will track optimization loop progress
        sprintf(fname, "%s/linoptval.txt", piaacmcsimul_var.piaacmcconfdir);
        fp = fopen(fname, "w");
        fclose(fp);

        //	list_image_ID();
        //
        // LINEAR OPTIMIZATION AROUND CURRENT POINT
        //

        iterOK=1;
        iter = 0;
        oldval = 1.0;
        data.image[IDstatus].array.UI16[0] = 9;

        // while # of iterations < piaacmcsimul_var.linopt_NBiter
        //  and the ojective changes by more than 2% after the second iteration
        //  and something about NBlinoptgain ???????
        while(iterOK==1)//        for(iter=0; iter<piaacmcsimul_var.linopt_NBiter; iter++)
        {
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 10;
            printf("Iteration %ld/%ld\n", iter, piaacmcsimul_var.linopt_NBiter);
            fflush(stdout);

            // array for collecting dark hole mode derivatives
            // stores derivative of output vector against input parameters
            IDmodes = create_3Dimage_ID("DHmodes", size1Dvec, 1, piaacmcsimul_var.linopt_number_param);
            // 2D array for diagnostic display
            IDmodes2D = create_2Dimage_ID("DHmodes2D", size1Dvec, piaacmcsimul_var.linopt_number_param); //TEST

            // get ready to update optimization tracking file
            sprintf(fname, "%s/linoptval.txt", piaacmcsimul_var.piaacmcconfdir);
            fp = fopen(fname, "a");
            fprintf(fp, "### PIAACMC_FPM_FASTDERIVATIVES = %d\n", piaacmcsimul_var.PIAACMC_FPM_FASTDERIVATIVES);
            fclose(fp);
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 11;

            printf("Compute local derivatives of output vector against input focal plane mask zones (sags)\n");
            fflush(stdout);

            if( piaacmcsimul_var.PIAACMC_FPM_FASTDERIVATIVES == 1) // TO BE USED ONLY FOR FOCAL PLANE MASK OPTIMIZATION
            {   // this only happens in mode 13
                // the fast derivative mode only works for focal plane mask optimization, for which derivatives against sag values can be comptuted by simple rotation of pre-computed vectors from mode 11
                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 12;
                optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm
                //				ID = create_2Dimage_ID("DHmodes2Dtest", size1Dvec, piaacmcsimul_var.linopt_number_param);

                printf("Computing %ld derivatives ", (long) data.image[piaacmc[0].zonezID].md[0].size[0]);
                fflush(stdout);
                for(mz=0; mz<data.image[piaacmc[0].zonezID].md[0].size[0]; mz++) // loop over mask zones
                {
                    printf(" %ld", mz);
                    fflush(stdout);
                    // actually compute the derivative
                    // fpmresp_array is results from mode 11
                    // from mode 13 above:
                    //      fpmresp_array = data.image[IDfpmresp].array.D;
                    //      zonez_array = data.image[piaacmc[0].zonezID].array.D;
                    // dphadz_array was computed in mode 13 shortly afterwards
                    // outtmp_array is output
                    PIAACMCsimul_achromFPMsol_eval_zonezderivative(mz, piaacmcsimul_var.fpmresp_array, piaacmcsimul_var.zonez_array, piaacmcsimul_var.dphadz_array, piaacmcsimul_var.outtmp_array, piaacmcsimul_var.vsize, data.image[piaacmc[0].zonezID].md[0].size[0], piaacmc[0].nblambda);
                    for(ii=0; ii<size1Dvec0; ii++)
                        data.image[IDmodes].array.F[mz*size1Dvec+ii] = piaacmcsimul_var.outtmp_array[ii]*piaacmcsimul_var.linopt_paramdelta[mz];
                }
                // derivatives of regularization values against sag values can also be computed analytically without requiring diffraction propagation
                if(piaacmcsimul_var.linopt_REGFPMSAG == 1)
                {
                    ID = piaacmc[0].zonezID;
                    // following should be derivative of (sag/coeff)^alpha
                    // w.r.t. sag
                    for(mz=0; mz < data.image[ID].md[0].size[0]; mz++)
                        data.image[IDmodes].array.F[mz*size1Dvec + (size1Dvec0+mz)] = (piaacmc[0].fpmsagreg_alpha/piaacmc[0].fpmsagreg_coeff) * pow( data.image[ID].array.D[mz]/piaacmc[0].fpmsagreg_coeff, piaacmc[0].fpmsagreg_alpha-1.0)*piaacmcsimul_var.linopt_paramdelta[mz];
                }


                // TEST diagnostic
                memcpy(data.image[IDmodes2D].array.F, data.image[IDmodes].array.F, sizeof(float)*size1Dvec*piaacmcsimul_var.linopt_number_param);
                save_fl_fits("DHmodes2D", "!test_DHmodes2D.fits");

                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 13;

                printf("Done computing derivatives (FAST MODE)\n");
                fflush(stdout);
            }
            else // ONLY FOR PIAA SHAPES OPTIMIZATION
            {
                // derivatives against PIAA shapes must be computed numerically

                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 14;
                for(i=0; i<piaacmcsimul_var.linopt_number_param; i++)
                {
                    optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm
                    // get delta on-axis PSF response to the change in piaacmcsimul_var.linopt_paramdelta
                    // to later give derivative w.r.t. piaacmcsimul_var.linopt_paramdelta
                    if(piaacmcsimul_var.linopt_paramtype[i]==_DATATYPE_FLOAT)
                    {
                        *(piaacmcsimul_var.linopt_paramvalf[i]) += (float) piaacmcsimul_var.linopt_paramdelta[i];
                        val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, 0);
                    }
                    else
                    {
                        *(piaacmcsimul_var.linopt_paramval[i]) += piaacmcsimul_var.linopt_paramdelta[i];
                        val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, 0);
                    }

                    //      sprintf(fname,"!%s/imvect_%02ld.fits", piaacmcsimul_var.piaacmcconfdir, i);
                    //       save_fits("imvect", fname);
                    ID = image_ID("imvect");


                    sprintf(fname, "%s/linoptval.txt", piaacmcsimul_var.piaacmcconfdir);
                    fp = fopen(fname, "a");
                    fprintf(fp, "# %5ld/%5ld %5ld/%5ld %6ld %6ld    %20.15g %20.15g %20.15g      %20.15g\n", iter, piaacmcsimul_var.linopt_NBiter, i, piaacmcsimul_var.linopt_number_param, data.image[ID].md[0].nelement, size1Dvec, piaacmcsimul_var.linopt_paramdelta[i], val, valref, bestval);
                    fclose(fp);

                    // re-package vector into 1D array and add regularization terms
                    // evaluation vector is "imvect1D", ID = ID1D
                    // similar to vecDHref before
                    ID1D = create_2Dimage_ID("imvect1D", size1Dvec, 1);
                    // fill in the evaluation point portion
                    for(ii=0; ii<data.image[ID].md[0].nelement; ii++)
                        data.image[ID1D].array.F[ii] = data.image[ID].array.F[ii];

                    if(piaacmcsimul_var.linopt_REGPIAASHAPES==1)
                    {   // fill in the shape regularization value's response to piaacmcsimul_var.linopt_paramdelta using the
                        // same formulas as before with the delta PSF as input
                        ii = PIAACMCsimul_regularization_PIAAshapes_add1Dvector(ID1D, ii);
                    }


                    delete_image_ID("imvect"); // has been imbedded into imvect1D

                    // restore original state (return to original staring point)
                    if(piaacmcsimul_var.linopt_paramtype[i]==_DATATYPE_FLOAT)
                        *(piaacmcsimul_var.linopt_paramvalf[i]) -= (float) piaacmcsimul_var.linopt_paramdelta[i];
                    else
                        *(piaacmcsimul_var.linopt_paramval[i]) -= piaacmcsimul_var.linopt_paramdelta[i];


                    // compute actual derivative as first difference from reference
                    // this is the starting derivative
                    for(ii=0; ii<data.image[ID1D].md[0].nelement; ii++)
                        data.image[IDmodes].array.F[i*data.image[ID1D].md[0].nelement+ii] = (data.image[ID1D].array.F[ii] - data.image[ID1Dref].array.F[ii]);


                    //    printf("%3ld %g %g\n", i, val, valref);

                    // create diagnostic image
                    ID = create_2Dimage_ID("DHmodes2D", size1Dvec, piaacmcsimul_var.linopt_number_param);
                    for(ii=0; ii<data.image[IDmodes].md[0].nelement; ii++)
                        data.image[ID].array.F[ii] = data.image[IDmodes].array.F[ii];

                    sprintf(fname, "!%s/DMmodes.fits", piaacmcsimul_var.piaacmcconfdir);
                    save_fits("DHmodes2D", fname);

                    delete_image_ID("DHmodes2D");
                }
                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 15;
            }


            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 16;
            // print the results to file for human tracking
            sprintf(fname, "%s/linoptval.txt", piaacmcsimul_var.piaacmcconfdir);
            fp = fopen(fname, "a");
            fprintf(fp, "### scanning gain \n");
            fprintf(fp, "### <alphareg>  <gain>  <contrast>\n");
            fclose(fp);
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 17;

            // first three arguments are names of the input arrays
            // vecDHref1D is the input data
            // DHmodes is the basis of modes to expand vecDHref1D into
            // DHmask weights elements of vecDHref1D (nominally all weights = 1)
            // 4th arg is the pseudoinverse eigenvalue (via eigenvalue decomposition)
            // this decomposes vecDHref1D into the DHmodes
            // 5th arg is the output: optcoeff0*DHmodes = vecDHref1D
            // computed via pseudoinverse of DHmodes
            // This decomposition facilitates the cancellation of vecDHref1D by
            // searching in the DHmodes basis
            //
            // use three cutoff values to give three options for future evaluation
            // smallest cutoff values produce the largest changes (are least well conditioned)
            printf("ref = 0.1   -- ");
            fflush(stdout);
            linopt_imtools_image_fitModes("vecDHref1D", "DHmodes", "DHmask", 0.1, "optcoeff0", 0);
            printf("- DONE\n");
            fflush(stdout);

            printf("ref = 0.01  -- ");
            fflush(stdout);
            linopt_imtools_image_fitModes("vecDHref1D", "DHmodes", "DHmask", 0.01, "optcoeff1", 0);
            printf("- DONE\n");
            fflush(stdout);

            printf("ref = 0.001 -- ");
            fflush(stdout);
            linopt_imtools_image_fitModes("vecDHref1D", "DHmodes", "DHmask", 0.001, "optcoeff2", 0);
            printf("- DONE\n");
            fflush(stdout);

            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 18;

            // initialize zero "optimal" vector optvec giving direction to move in the search for the min
            arith_image_cstmult("optcoeff0", 0.0, "optvec"); // create optimal vector
            IDoptvec = image_ID("optvec");
            // initialize the objective value
            initbestval = 0;
            bestval = valref;

            alphareg = 1.0; // has no effect (see next loop)

            // say something here ???
            NBlinoptgain = 0;


            scangainfact = 1.2;
            alphascaninit = 0;
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 19;
            // alphareg controls linear combinations of the directions optcoeff0,1,2
            // alphareg = 0 => moving along optcoeff0
            // alphareg = 0.5 => moving along optcoeff1
            // alphareg = 1 => moving along optcoeff2
            // with interpolation
            // look in 5 steps if alphareg += 0.2
            for(alphareg=0.0; alphareg<1.01; alphareg += 0.2)
            {
                // produce piecewise linear interoplation coefficients
                acoeff0 = 1.0 - 2.0*alphareg; // ranges from -1 to 1
                if(acoeff0<0.0)
                    acoeff0 = 0.0; // clip to 0 to 1, = 0 if alphareg > 0.5

                acoeff1 = 1.0 - fabs(2.0*(alphareg-0.5)); // two lines: from 0,1 to .5,0 to 1,1

                acoeff2 = 2.0*alphareg - 1.0; // ranges from -1 to 1
                if(acoeff2<0.0)
                    acoeff2 = 0.0; // clip to 0 to 1, = 0 if alphareg < 0.5

                // sum of acoeff0,1,2 = 1 at all values of alphareg.

                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 20;

                // optcoeff0m = acoeff0*optcoeff0 etc.
                arith_image_cstmult("optcoeff0", acoeff0, "optcoeff0m");
                arith_image_cstmult("optcoeff1", acoeff1, "optcoeff1m");
                arith_image_cstmult("optcoeff2", acoeff2, "optcoeff2m");

                // diagnostic
                save_fl_fits("optcoeff0", "!optcoeff0.fits");//TEST
                save_fl_fits("optcoeff1", "!optcoeff1.fits");
                save_fl_fits("optcoeff2", "!optcoeff2.fits");

                // optcoeff01m = acoeff0*optcoeff0 + acoeff1*optcoeff1
                arith_image_add("optcoeff0m", "optcoeff1m", "optcoeff01m");
                // optcoeff = acoeff0*optcoeff0 + acoeff1*optcoeff1 + acoeff2*optcoeff2
                arith_image_add("optcoeff01m", "optcoeff2m", "optcoeff");
                // optcoeff now has our search direction
                delete_image_ID("optcoeff0m");
                delete_image_ID("optcoeff1m");
                delete_image_ID("optcoeff2m");
                delete_image_ID("optcoeff01m");

                ID = image_ID("optcoeff");
                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 21;

                // do linear scan along the direction optcoeff from current parameter location
                linscanOK = 1;
                // size of step in this direction
                scangain = 0.0; //scanstepgain; overriden below
                val = 100000000000.0; // initialize minimization objective to a big number
                bestgain = 0.0;
                // iteration counter of steps in the current direction optcoeff
                k = 0;

                // if(alphascaninit==1)
                //		scangain = bestgain/scangainfact/scangainfact/scangainfact/scangainfact/scangainfact;
                //	alphascaninit = 1;

                // scangain is our location along the direction optcoeff
                scangain = 0.001;
                //              scanstepgain = 0.000001; // TEST
                //              scangainfact = 1.00001; // TEST

                // while objective value < previous value and we've taken no more than than 90 steps
                while(linscanOK==1)
                {
                    // for state tracking and statistics
                    data.image[IDstatus].array.UI16[0] = 22;

                    // compute offsets
                    ID = image_ID("optcoeff"); // direction vector
                    linoptlimflagarray[k] = 0;
                    // step each parameter by optcoeff
                    for(i=0; i<piaacmcsimul_var.linopt_number_param; i++) // looping over parameters
                    {
                        // compute step delta for this parameter
                        // image[ID] = optcoeff is a derivative w.r.t. piaacmcsimul_var.linopt_paramdelta, so has dimension
                        // (parameter dimension)/(piaacmcsimul_var.linopt_paramdelta dimension) so we have to mulitply by
                        // piaacmcsimul_var.linopt_paramdelta to put our step in physical parameter units
                        // negative because we want to cancel the value from the delta PSF
                        piaacmcsimul_var.linopt_paramdeltaval[i] = -scangain * data.image[ID].array.F[i] * piaacmcsimul_var.linopt_paramdelta[i];
                        if(piaacmcsimul_var.linopt_paramdeltaval[i]<-piaacmcsimul_var.linopt_parammaxstep[i]) // if the step is too large in the negative direction
                        {
                            printf("MIN LIMIT [%3ld   %20g]   %20g -> ", i, piaacmcsimul_var.linopt_paramdelta[i], piaacmcsimul_var.linopt_paramdeltaval[i]); //TEST
                            piaacmcsimul_var.linopt_paramdeltaval[i] = -piaacmcsimul_var.linopt_parammaxstep[i]; // set it to the negative largest allowed step
                            printf(" %20g\n", piaacmcsimul_var.linopt_paramdeltaval[i]); //TEST
                            linoptlimflagarray[k] = 1;
                        }
                        if(piaacmcsimul_var.linopt_paramdeltaval[i]>piaacmcsimul_var.linopt_parammaxstep[i])// if the step is too large in the positive direction
                        {
                            printf("MAX LIMIT [%3ld   %20g]   %20g -> ", i, piaacmcsimul_var.linopt_paramdelta[i], piaacmcsimul_var.linopt_paramdeltaval[i]); //TEST
                            piaacmcsimul_var.linopt_paramdeltaval[i] = piaacmcsimul_var.linopt_parammaxstep[i]; // set it to the positive largest allowed step
                            printf(" %20g\n", piaacmcsimul_var.linopt_paramdeltaval[i]); //TEST
                            linoptlimflagarray[k] = 1;
                        }

                        // apply offsets to the global data object via the pointers piaacmcsimul_var.linopt_paramvalf, which
                        // point into the (hopefully) coorect locations of each parameter's value in the
                        // data object
                        if(piaacmcsimul_var.linopt_paramtype[i]==_DATATYPE_FLOAT)
                        {
                            if(  *(piaacmcsimul_var.linopt_paramvalf[i]) + (float) piaacmcsimul_var.linopt_paramdeltaval[i]  > piaacmcsimul_var.linopt_parammax[i] )
                                // if we're about to step too far, set the step to the
                                // limit - current parameter value, so we step to the limit
                                piaacmcsimul_var.linopt_paramdeltaval[i] = piaacmcsimul_var.linopt_parammax[i] - *(piaacmcsimul_var.linopt_paramvalf[i]);

                            if(  *(piaacmcsimul_var.linopt_paramvalf[i]) + (float) piaacmcsimul_var.linopt_paramdeltaval[i]  < piaacmcsimul_var.linopt_parammin[i] )
                                // if we're about to step too far, set the step to the
                                // limit - current parameter value, so we step to the limit
                                piaacmcsimul_var.linopt_paramdeltaval[i] = piaacmcsimul_var.linopt_parammin[i] - *(piaacmcsimul_var.linopt_paramvalf[i]);
                            // take the actual step (piaacmcsimul_var.linopt_paramvalf is a 1D array of pointers)
                            *(piaacmcsimul_var.linopt_paramvalf[i]) += (float) piaacmcsimul_var.linopt_paramdeltaval[i];
                        }
                        else // same for the double case
                        {
                            if(  *(piaacmcsimul_var.linopt_paramval[i]) + piaacmcsimul_var.linopt_paramdeltaval[i]  > piaacmcsimul_var.linopt_parammax[i] )
                                piaacmcsimul_var.linopt_paramdeltaval[i] = piaacmcsimul_var.linopt_parammax[i] - *(piaacmcsimul_var.linopt_paramval[i]);

                            if(  *(piaacmcsimul_var.linopt_paramval[i]) + piaacmcsimul_var.linopt_paramdeltaval[i]  < piaacmcsimul_var.linopt_parammin[i] )
                                piaacmcsimul_var.linopt_paramdeltaval[i] = piaacmcsimul_var.linopt_parammin[i] - *(piaacmcsimul_var.linopt_paramval[i]);

                            *(piaacmcsimul_var.linopt_paramval[i]) += piaacmcsimul_var.linopt_paramdeltaval[i];
                        }
                    }
                    // store the current objective value for later comparison
                    valold = val;
                    // for state tracking and statistics
                    data.image[IDstatus].array.UI16[0] = 23;

                    // compute new state and compute assossiated evaluation metric
                    // using the modified global data object
                    val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, 0);
                    valContrast = val; // contrast component of the evaluation metric
                    // we've now only done the light portion

                    // add regularization component of the evaluation metrix
                    // first, compute and add PIAA shape regularization value (val0) if applicable
                    // this is the same as the previous computation of val0
                    val0 = 0.0;
                    if(piaacmcsimul_var.linopt_REGPIAASHAPES==1)
                    {
                        val0 = PIAACMCsimul_regularization_PIAAshapes_value();
                        val += val0;
                    }

                    // add sag regularization (val1) if applicable, as before
                    val1 = 0.0;
                    if(piaacmcsimul_var.linopt_REGFPMSAG == 1)
                    {
                        val1 = PIAACMCsimul_regularization_fpmsag_value();;
                        val += val1;
                    }
                    // val is now our complete objective!! Yay!!


                    // for state tracking and statistics
                    data.image[IDstatus].array.UI16[0] = 24;
                    // print it for monitoring
                    sprintf(fname, "%s/linoptval.txt", piaacmcsimul_var.piaacmcconfdir);
                    fp = fopen(fname, "a");
                    // printf the first part of the line reporting current values
                    fprintf(fp, "##  [ %5ld / %5ld ]   %5.3f  %12lf  %12g   (reg = %12g [%1d] %12g [%1d]   contrast = %20g)       [%d] [%ld]", iter, piaacmcsimul_var.linopt_NBiter, alphareg, scangain, val, val0, piaacmcsimul_var.linopt_REGPIAASHAPES, val1, piaacmcsimul_var.linopt_REGFPMSAG, valContrast, linoptlimflagarray[k], piaacmcsimul_var.linopt_number_param);

                    // now add text indicating status and complete line
                    // and store all parameters for the current best solution
                    if((val<bestval)||(initbestval==0))
                    {
                        for(i=0; i<piaacmcsimul_var.linopt_number_param; i++)
                            if(piaacmcsimul_var.linopt_paramtype[i]==_DATATYPE_FLOAT)
                                data.image[IDoptvec].array.F[i] = *(piaacmcsimul_var.linopt_paramvalf[i]);
                            else
                                data.image[IDoptvec].array.F[i] = (float) *(piaacmcsimul_var.linopt_paramval[i]);
                        bestval = val;
                        if(initbestval == 0)
                            fprintf(fp, " ===== START POINT =====\n");
                        else
                            fprintf(fp, "  -> BEST VECTOR =======\n");
                        bestgain = scangain;
                        initbestval = 1;
                    }
                    else
                    {
                        fprintf(fp, " bestval = %12g\n", bestval);
                    }
                    fclose(fp);

                    // remove offsets returning the global data object to its original state
                    for(i=0; i<piaacmcsimul_var.linopt_number_param; i++)
                    {
                        if(piaacmcsimul_var.linopt_paramtype[i]==_DATATYPE_FLOAT)
                            *(piaacmcsimul_var.linopt_paramvalf[i]) -= (float) piaacmcsimul_var.linopt_paramdeltaval[i];
                        else
                            *(piaacmcsimul_var.linopt_paramval[i]) -= piaacmcsimul_var.linopt_paramdeltaval[i];
                    }
                    // for state tracking and statistics
                    data.image[IDstatus].array.UI16[0] = 25;

                    // store the current position and value
                    linoptgainarray[k] = scangain;
                    linoptvalarray[k] = val;
                    k++; // next step

                    // test to see if we're no longer getting better
                    if(val<valold)
                    {
                        linscanOK = 1; // if we're getting better keep going
                        // bestgain = scangain;
                        // scangain += scanstepgain;
                    }
                    else // otherwise stop stepping
                        linscanOK = 0;



                    if(k>90)  // stop if we've taken too many steps
                        linscanOK = 0;

                    // increment our location along the line
                    scangain += scanstepgain; // scanstepgain is an initilizaed function local (currently 0.001)
                    scangain *= scangainfact; // causes later steps to be larger
                    // (implicit scangainfact^n for the nth step)
                }
                // NBlinoptgain is counting the largest number of steps needed in this inner loop
                // stepping in the current direction.  When this is small (< 3) we declare victory
                // and stop the outer linear optimization
                if(k>NBlinoptgain)
                    NBlinoptgain = k;

                delete_image_ID("optcoeff"); // delete the current direction
                // for state tracking and statistics
                data.image[IDstatus].array.UI16[0] = 26;
            }
            // best solution after this linear linescan is stored in IDoptvec
            delete_image_ID("optcoeff0");
            delete_image_ID("optcoeff1");
            delete_image_ID("DHmodes");

            // we've now found the minimum using the three directions from the
            // alternative decompositions of the parameter space with the DHmodes basis
            // (with different conditioning).
            // Now we check the result by recomputing from scratch the objective with the current
            // (hopefully) optimal parameters.  Belts and suspenders

            // At this point the global data object has been restored to its original
            // (non-optimal) state.

            // (re-)compute best solution identified in previous linescan
            // update state to best solution (IDoptvec), setting the global data object
            // to the optimal state
            for(i=0; i<piaacmcsimul_var.linopt_number_param; i++)
            {
                if(piaacmcsimul_var.linopt_paramtype[i]==_DATATYPE_FLOAT)
                    *(piaacmcsimul_var.linopt_paramvalf[i]) = data.image[IDoptvec].array.F[i];
                else
                    *(piaacmcsimul_var.linopt_paramval[i]) = (double) data.image[IDoptvec].array.F[i];
            }
            valold = val;
            // compute contrast metric component -> val using the data object in the latest optimal state
            val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, 0);

            // add PIAA shape regularization component (val0) if applicable
            // same val0 and val1 computations are before, using the latest optimal state
            val0 = 0.0;
            if(piaacmcsimul_var.linopt_REGPIAASHAPES==1)
            {
                val0 = PIAACMCsimul_regularization_PIAAshapes_value();
                val += val0;
            }

            // add sag regularization component if applicable
            val1 = 0.0;
            if(piaacmcsimul_var.linopt_REGFPMSAG == 1)
            {
                val1 = PIAACMCsimul_regularization_fpmsag_value();;
                val += val1;
            }

            // now val is the objective including any desired regularization terms using the
            // latest optimal solution


            printf("gain: %lf -> val = %20g\n", bestgain, val);
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 27;



            // update reference state evaluation vector of optimal results including evaluation zone values
            // and desired regularization terms
            // which sets up starting the next iteration at the best solution

            ID1Dref = image_ID("vecDHref1D");
            ID = image_ID("imvect");
            // first fill in evaluation zone (complex) values
            for(ii=0; ii<data.image[ID].md[0].nelement; ii++)
                data.image[ID1Dref].array.F[ii] = data.image[ID].array.F[ii];
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 28;

            // now fill in the regularization terms if desired
            // same code as before.
            if(piaacmcsimul_var.linopt_REGPIAASHAPES==1)
            {
                ii = PIAACMCsimul_regularization_PIAAshapes_add1Dvector(ID1Dref, ii);
            }


            if(piaacmcsimul_var.linopt_REGFPMSAG == 1)
            {
                ID = piaacmc[0].zonezID;
                for(jj=0; jj < data.image[ID].md[0].size[0]; jj++)
                {
                    data.image[ID1Dref].array.F[ii] = pow(data.image[ID].array.D[jj]/piaacmc[0].fpmsagreg_coeff, piaacmc[0].fpmsagreg_alpha);
                    ii++;
                }
            }

            delete_image_ID("imvect");
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 29;









            // print out current best value for tracking
            sprintf(fname, "%s/linoptval.txt", piaacmcsimul_var.piaacmcconfdir);
            fp = fopen(fname, "a");
            if(fp==NULL)
            {
                printf("ERROR: cannot open file \"%s\"\n", fname);
                exit(0);
            }
            fprintf(fp, "-> %5ld    %20g <- %20g \n", iter, val, valref);
            printf("%5ld %20g %20g \n", iter, val, valref);
            fflush(stdout);
            fclose(fp);

            // save current best value and reference value in globals
            piaacmcsimul_var.PIAACMCSIMUL_VAL = val;
            piaacmcsimul_var.PIAACMCSIMUL_VALREF = valref;



            // Nominally if we're in this linear
            // optimization piaacmcsimul_var.PIAACMC_fpmtype = 1, so the next line is not executed
            if(piaacmcsimul_var.PIAACMC_fpmtype==0) // in the idealized PIAACMC case
                piaacmc[0].fpmaskamptransm = data.image[piaacmc[0].zoneaID].array.D[0]; // required to ensure that the new optimal focal plane mask transmission is written to disk

            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 30;


            // tracking diagnostics giving behavior of the modes by iteration
            sprintf(dirname, "%s_linopt", piaacmcsimul_var.piaacmcconfdir);
            PIAACMCsimul_savepiaacmcconf(dirname); // staging area
            sprintf(command, "rsync -au --progress %s/* ./%s/", dirname, piaacmcsimul_var.piaacmcconfdir);
            r = system(command);

            sprintf(command, "cp %s/piaa0Cmodes.fits %s/piaa0Cmodes.%04ld.fits", dirname, dirname, iter);
            r = system(command);
            sprintf(command, "cp %s/piaa0Fmodes.fits %s/piaa0Fmodes.%04ld.fits", dirname, dirname, iter);
            r = system(command);

            sprintf(command, "cp %s/piaa1Cmodes.fits %s/piaa1Cmodes.%04ld.fits", dirname, dirname, iter);
            r = system(command);
            sprintf(command, "cp %s/piaa1Fmodes.fits %s/piaa1Fmodes.%04ld.fits", dirname, dirname, iter);
            r = system(command);

            sprintf(command, "cp %s/piaacmcparams.conf %s/piaacmcparams.%04ld.conf", dirname, dirname, iter);
            r = system(command);

            if(file_exists(stopfile)==1)
                iterOK = 0;

            // Figure out if current loop should continue optimization
            // if optimization ends, then iterOK set to 0
            // if we've reached the allowed number of iterations
            if(iter==piaacmcsimul_var.linopt_NBiter)
                iterOK = 0;
            if(iter>2)
            {
                if(val>0.98*oldval) // if after second iteration and we'be improved by less than 10%
                    iterOK = 0;
            }

            if(NBlinoptgain<3) // if we've stopped moving much
                iterOK = 0;


            // set up for next iteration
            oldval = val;
            iter++;

            printf("END OF LOOP ITERATION\n");
            fflush(stdout);
            // for state tracking and statistics
            data.image[IDstatus].array.UI16[0] = 31;
        }
        printf(" ============ END OF OPTIMIZATION LOOP ======= \n");
        // for state tracking and statistics
        data.image[IDstatus].array.UI16[0] = 32;
    } // end of if (piaacmcsimul_var.LINOPT==1): done with the linear optimization

    piaacmcsimul_var.LINOPT = 0;



    //  PIAACMCsimul_savepiaacmcconf("piaacmc0");
    //  PIAACMCsimul_loadpiaacmcconf("piaacmc0");
    // PIAACMCsimul_savepiaacmcconf("piaacmc1");
    //exit(0);

    return 0;
}





