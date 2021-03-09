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
#include <assert.h>
#include <string.h>



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
     * ## Mode 1: Optimize Lyot stop positions
    
    Lyot stop positions are encoded as piaacmc[0].LyotStop_zpos
    
    there can be multiple LyotStop_zpos
    
    Vary these zpos, looking for the best contrast returned by PIAACMCsimul_computePSF
    
    Search is performed by iterative refined marching
    **/ 
int PIAACMCsimul_exec_optimize_lyot_stop_position()
{   
	long IDv;
	double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;    
    double range, stepsize;
    FILE *fp;     
    long ls, ls1;
    long iter;
    long NBiter = 1000;    
    double parambest[10000]; // for scanning
    double paramref[10000];    
    int loopOK;    
	double valbest;
    char fnamelog[500];	    
    long elem0;
    
        printf("=================================== mode 001 ===================================\n");
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
        


        // init as in mode 0
        PIAACMCsimul_initpiaacmcconf(0, fpmradld, centobs0, centobs1, 0, 1);
        PIAACMCsimul_makePIAAshapes(piaacmc, 0);
        optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm

        // initialization
		PIAACMCsimul_init(piaacmc, 0, 0.0, 0.0); // necessary to initialize optical design
        // set initial lyot stop marching range (current position +- range)
        if((IDv=variable_ID("PIAACMC_lsoptrange"))!=-1)
            range = data.variable[IDv].value.f; // from cli
        else
            range = 3.0; // default, in meters
        stepsize = range/3.0; // initial march stepsize
        // store initial Lyot stop positions
        for(ls=0; ls<piaacmc[0].NBLyotStop; ls++) // NBLyotStop = length(LyotStop_zpos)
            paramref[ls] = piaacmc[0].LyotStop_zpos[ls];
        NBiter = 4; // number of iterations

        // start up a log
        sprintf(fnamelog, "%s/result_LMpos.log", piaacmcsimul_var.piaacmcconfdir);
        fp = fopen(fnamelog, "w");
        fclose(fp);



        // pick another initial march stepsize
        stepsize = range/5.0;
        // start the iterative refined march
        for(iter=0; iter<NBiter; iter++)
        {
            // for each Lyot stop, find its best position
            for(ls=0; ls<piaacmc[0].NBLyotStop; ls++)
            {
                // start position for march.  paramref is current best value
                piaacmc[0].LyotStop_zpos[ls] = paramref[ls]-range;
                // current best position
                parambest[ls] = piaacmc[0].LyotStop_zpos[ls];

                loopOK = 1;
                valbest = 1.0;
                
                // march to the other other end of range
                while(piaacmc[0].LyotStop_zpos[ls]<paramref[ls]+range)
                {
                    long elem;
                    double val;
                    
                    elem0 = 6; // elem0 is the starting point of the PSF propagation.  This is a staring default

                    // look for the element called "Lyot mask 0" as the actual starting point
					elem0 = -1;
					printf("Number of elements = %ld\n", optsyst[0].NBelem);
					assert ( optsyst[0].NBelem > 0 );
					
					for(elem=0; elem<optsyst[0].NBelem; elem++)
					{
						printf("elem %ld :  %s\n", elem, optsyst[0].name[elem]);
                        if(strcmp("Lyot mask 0", optsyst[0].name[elem])==0)
                            elem0 = elem;
					}
					assert( elem0 != -1); // throw a message if this was not found


                    optsyst[0].keepMem[elem0] = 1; // save this element and reuse

                    // compute the PSF for this Lyot stop position, returning contrast in the evaluation zone
                    val = PIAACMCsimul_computePSF(0.0, 0.0, elem0, optsyst[0].NBelem, 0, 0, 0, 0);

                    // if this is the best contrast for this stop, save it for it and the position of this stop
                    if(val<valbest)
                    {
                        parambest[ls] = piaacmc[0].LyotStop_zpos[ls];
                        valbest = val;
                    }

                    // say what's happening
                    fp = fopen(fnamelog, "a");
                    for(ls1=0; ls1<piaacmc[0].NBLyotStop; ls1++)
                        fprintf(fp," %lf", piaacmc[0].LyotStop_zpos[ls1]);
                    fprintf(fp, " %g\n", val);
                    fclose(fp);

                    // march along by the step size
                    piaacmc[0].LyotStop_zpos[ls] += stepsize;
                }
                printf("BEST SOLUTION :  ");
                paramref[ls] = parambest[ls]; // update best position for this stop
                piaacmc[0].LyotStop_zpos[ls] = paramref[ls]; // store in case this is last iteration
                printf(" %lf", parambest[ls]);
                printf(" %g\n", valbest);
            }

            fp = fopen(fnamelog, "a");
            fprintf(fp, "\n");
            fclose(fp);

            // reduce the range and stepsize, refining the march
            range *= 0.3;
            stepsize = range/3.0;
        }
        // store all best positions  Done!!
        for(ls=0; ls<piaacmc[0].NBLyotStop; ls++)
            piaacmc[0].LyotStop_zpos[ls] = parambest[ls];
        PIAACMCsimul_savepiaacmcconf( piaacmcsimul_var.piaacmcconfdir);

		return 0;
}

