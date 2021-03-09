/**
 * @file    PIAAACMCsimul_exec_optimize_lyot_stop_position.c
 * @brief   Optimize Lyot stop(s) conjugations
 * 
 *  
 * @author  O. Guyon
 * @date    25 nov 2017
 *
 * 
 * @bug No known bugs.
 * 
 */



// System includes
#include <stdio.h>
#include <stdlib.h>
//#include <assert.h>
//#include <string.h>



// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"



extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTSYST *optsyst;

extern OPTPIAACMCDESIGN *piaacmc;







		/** 
		 * ---
		 * 
		 * ## Mode 2: Optimize focal plane mask transmission for monochromatic idealized PIAACMC
        
        For monochromatic, idealized PIAACMC, find the scalar transimssion of the uniform focal plane mask
        that provides best contrast in the evaluation zone
        
        Very similar to the Lyot stop search in mode 1: iterative refined marching, changing the
        the transmission value piaacmc[0].fpmaskamptransm, which
        is between 0 and 1
        
        Uses single on-axis light source
        **/
int PIAACMCsimul_exec_optimize_fpmtransmission()
{   
	long IDv;
	double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;    
    double range, stepsize;
    FILE *fp;     
    int loopOK;    
	double valbest;
    long iter;
    long NBiter = 1000;
    double parambest[10000]; // for scanning
    double paramref[10000];       
    char fnamelog[500];	        
    char fname[800];
     char command[1000];   
    	
        printf("=================================== mode 002 ===================================\n");

        // load some more cli variables
    if( (IDv = variable_ID("PIAACMC_centobs0")) != -1)
        centobs0 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_centobs1")) != -1)
        centobs1 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }

	


        /// ### Initialize as in mode 0
        PIAACMCsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 1);
        PIAACMCsimul_makePIAAshapes(piaacmc, 0);
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm

        /// ### Initialize search range and step
        range = 0.3;
        stepsize = range/3.0;
        paramref[0] = piaacmc[0].fpmaskamptransm; // initialized in PIAACMCsimul_initpiaacmcconf
        NBiter = 6;


		/// ### Scan parameter value
        sprintf(fnamelog, "%s/result_fpmt.log", piaacmcsimul_var.piaacmcconfdir);
        fp = fopen(fnamelog, "w");
        fclose(fp);

        for(iter=0; iter<NBiter; iter++)
        {
            // starting point of march
            piaacmc[0].fpmaskamptransm = paramref[0]-range;
 
            // store current value as best
            parambest[0] = piaacmc[0].fpmaskamptransm;

            loopOK = 1;
            valbest = 1.0;

            /// While within the search loop :
            while(loopOK==1)
            {
				double val;
				
				printf("\n\n\n");
				
                piaacmcsimul_var.FORCE_CREATE_fpmza = 1; // forces creation of new focal plane mask in the next two routines

                /// - Call PIAACMCsimul_initpiaacmcconf() 
                PIAACMCsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 0);
                                
                // compute on-axis PSF of all optical elements returning contrast in evaluation zone
                // ************************* need to do all optsyst[0].NBelem?
				/// - call PIAACMCsimul_computePSF() to evaluate design
                val = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);

                if(val<valbest)
                {
                    // we have a better contrast!  Store it
                    parambest[0] = piaacmc[0].fpmaskamptransm;
                    valbest = val;
                }

				/// - write entry to output log file
				fp = fopen(fnamelog, "a");
                fprintf(fp," %+011.8lf", piaacmc[0].fpmaskamptransm);
                fprintf(fp, " %12g  %8ld %12g %12g\n", val, iter, range, stepsize);
                fclose(fp);
								
                /// -  increment parameter
                piaacmc[0].fpmaskamptransm += stepsize;
                
                // if we've reached the end of the range stop the loop
                if(piaacmc[0].fpmaskamptransm>paramref[0]+range+0.001*stepsize)
                    loopOK = 0;
            }


            printf("BEST SOLUTION :  ");
            // store best solution
            paramref[0] = parambest[0];
            printf(" %lf", parambest[0]);

            printf(" %g\n", valbest);


            fp = fopen(fnamelog, "a");
            fprintf(fp, "\n");
            fclose(fp);
            // refine range and stepsize
            range *= 0.3;
            stepsize = range/3.0;
        }
        // save final result
        piaacmc[0].fpmaskamptransm = parambest[0];
        PIAACMCsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 0); // why? **************************
        PIAACMCsimul_savepiaacmcconf( piaacmcsimul_var.piaacmcconfdir); // save final result to disk
        piaacmcsimul_var.FORCE_CREATE_fpmza = 0; // turning off to be good citizens

	return 0;
}
