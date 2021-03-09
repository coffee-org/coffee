/**
 * @file    PIAACMCsimul_FPM_process.c
 * @brief   PIAA-type coronagraph design
 * 
 * Can design both APLCMC and PIAACMC coronagraphs
 * 
 * @bug No known bugs.
 * 
 */



#include <stdlib.h>
#include <stdio.h>



// milk includes
#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"









long PIAACMC_FPM_process(const char *FPMsag_name, const char *zonescoord_name, long NBexp, const char *outname)
{
	long IDin;
	long NBzones;
	int atype;
	
	double *sagarray_in;
	double *sagarray_out;
	long zone;
	double sagmax, sagmin;
	long k;
	long NBsagsteps;
	
	double *sagstepval;
	FILE *fp;
	FILE *fpout;
	
	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif
	
	
	IDin = image_ID(FPMsag_name);
	NBzones = data.image[IDin].md[0].size[0];
	atype = data.image[IDin].md[0].datatype;
	
	switch (atype) {
		case _DATATYPE_DOUBLE:
		printf("atype = _DATATYPE_DOUBLE\n");
		break;
		case _DATATYPE_FLOAT:
		printf("atype = _DATATYPE_FLOAT\n");
		break;
		default :
		printf("ERROR: atype not supported\n");
		exit(0);
		break;
	}
	
	printf("%ld zones\n", NBzones);
	sagarray_in = (double*) malloc(sizeof(double)*NBzones);
	sagarray_out = (double*) malloc(sizeof(double)*NBzones);
	
	
	for(zone=0; zone<NBzones; zone++)
		{
			if(atype == _DATATYPE_FLOAT)
				sagarray_in[zone] = (double) data.image[IDin].array.F[zone];
			else
				sagarray_in[zone] = data.image[IDin].array.D[zone];
		}
	
	sagmin = sagarray_in[0];
	sagmax = sagarray_in[0];
	for(zone=1; zone<NBzones; zone++)
	{
		if(sagarray_in[zone]<sagmin)
			sagmin = sagarray_in[zone];
		if(sagarray_in[zone]>sagmax)
			sagmax = sagarray_in[zone];
	}
	
	printf("Sag range [um]  :   %10.8f  ->  %10.8f\n", sagmin*1.0e6, sagmax*1.0e6);
	NBsagsteps = 2;
	for(k=1; k<NBexp; k++)
		NBsagsteps *= 2;
	printf("NBsagsteps = %ld\n", NBsagsteps);
	
	
	fp = fopen("saglevels.dat", "w");
	
	sagstepval = (double*) malloc(sizeof(double)*NBsagsteps);
	for(k=0;k<NBsagsteps;k++)
		{
			sagstepval[k] = sagmin + (sagmax-sagmin)*k/NBsagsteps + 0.5*(sagmax-sagmin)/NBsagsteps;
			fprintf(fp, "%4ld     %10.8f\n", k, sagstepval[k]*1.0e6);
		}
	fclose(fp);
	printf("\n");
	
	
	fpout = fopen(outname, "w");
	
	
	//for(zone=0; zone<NBzones; zone++)
	for(zone=0; zone<NBzones; zone++)
	{
		k = (long) ((sagarray_in[zone]-sagmin)/((sagmax-sagmin)/NBsagsteps));
		if(sagarray_in[zone] > sagstepval[NBsagsteps-1])
			k = NBsagsteps-1;
//		printf("zone %4ld   %10.8f   %10.8f  %4ld\n", zone, sagarray_in[zone]*1e6, (sagarray_in[zone]-sagmin)*1e6, k);
		fprintf(fpout, "%4ld  %+10.8f  %4ld  %+10.8f\n", zone, sagarray_in[zone]*1e6, k, sagstepval[k]*1e6);
	}
	fclose(fpout);
	
	free(sagstepval);
	free(sagarray_in);
	free(sagarray_out);

	return 0;
}



