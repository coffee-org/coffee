/**
 * @file    PIAAACMCsimul_loadsavepiaacmcconf.c
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



// System includes
#include <stdio.h>



// milk includes
//   core modules
#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "image_basic/image_basic.h"
//   other modules



#include "PIAACMCsimul/PIAACMCsimul.h"





extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTPIAACMCDESIGN *piaacmc;









int PIAACMCsimul_loadpiaacmcconf(
		const char *dname
		)
{
    char command[1000];
    int r;
    FILE *fp;
    char fname[1500];
    char imname[1500];
    long i;

    int tmpi;
    long tmpl;
    float tmpf;
    double tmplf;

	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif


    sprintf(fname,"%s/piaacmcparams.conf", dname);
    printf("%s\n", fname);
    fp = fopen(fname, "r");
    if(fp==NULL)
    {
        printf("Configuration file \"%s\" does not exist (yet), using previously set configuration\n", fname);
        fflush(stdout);
        r = 0;
    }
    else
    {

        r = fscanf(fp, "%ld   NBradpts\n", &tmpl);
        piaacmc[0].NBradpts = tmpl;

         for(i=0; i<10; i++)
        {
            if(i<piaacmc[0].NBLyotStop)
            {
                sprintf(fname, "%s/LyotStop%ld.fits", dname, i);
                sprintf(imname, "lyotstop%ld", i);
                printf("Loading \"%s\" as \"%s\"\n", fname, imname);
                piaacmc[0].IDLyotStop[i] = load_fits(fname, imname, 1);
                r = fscanf(fp, "%lf   LyotStop_zpos %ld\n", &tmplf, &tmpl);
                piaacmc[0].LyotStop_zpos[i] = tmplf;
            }
            else
            {
                r = fscanf(fp, "%lf   LyotStop_zpos %ld\n", &tmplf, &tmpl);
                piaacmc[0].LyotStop_zpos[i] = tmplf;
            }
        printf("LYOT STOP %ld POS : %lf\n", i, tmplf);
            
        }


        sprintf(fname, "%s/piaa0Cmodes.fits", dname);
        piaacmc[0].piaa0CmodesID = load_fits(fname, "piaa0Cmodescoeff", 1);
        if(piaacmc[0].piaa0CmodesID==-1)
        {
            sprintf(fname, "%s/piaaref/piaa0Cmodes.fits", dname);
            piaacmc[0].piaa0CmodesID = load_fits(fname, "piaa0Cmodescoeff", 1);
        }


        sprintf(fname, "%s/piaa0Fmodes.fits", dname);
        piaacmc[0].piaa0FmodesID = load_fits(fname, "piaa0Fmodescoeff", 1);
        if(piaacmc[0].piaa0FmodesID==-1)
        {
            sprintf(fname, "%s/piaaref/piaa0Fmodes.fits", dname);
            piaacmc[0].piaa0FmodesID = load_fits(fname, "piaa0Fmodescoeff", 1);
        }

        sprintf(fname, "%s/piaa1Cmodes.fits", dname);
        piaacmc[0].piaa1CmodesID = load_fits(fname, "piaa1Cmodescoeff", 1);
        if(piaacmc[0].piaa1CmodesID==-1)
        {
            sprintf(fname, "%s/piaaref/piaa1Cmodes.fits", dname);
            piaacmc[0].piaa1CmodesID = load_fits(fname, "piaa1Cmodescoeff", 1);
        }

        sprintf(fname, "%s/piaa1Fmodes.fits", dname);
        piaacmc[0].piaa1FmodesID = load_fits(fname, "piaa1Fmodescoeff", 1);
        if(piaacmc[0].piaa1FmodesID==-1)
        {
            sprintf(fname, "%s/piaaref/piaa1Fmodes.fits", dname);
            piaacmc[0].piaa1FmodesID = load_fits(fname, "piaa1Fmodescoeff", 1);
        }


        r = fscanf(fp, "%f   fpmaskamptransm\n",    &tmpf);
        piaacmc[0].fpmaskamptransm = tmpf;

        r = 1;

        fclose(fp);
    }
    
    return(r);
}



/** @brief Assemble focal plane mask configuration name 
 * 
 */ 
int PIAACMCsimul_update_fnamedescr_conf()
{
	sprintf(piaacmcsimul_var.fnamedescr_conf, "s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_wb%02d", piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].nblambda);
	// debug from Justin
	// printf("s%d_l%04ld_sr%02ld_nbr%03ld_mr%03ld_ssr%02d_ssm%d_wb%02d\n", piaacmcsimul_var.PIAACMC_FPMsectors, (long) (1.0e9*piaacmc[0].lambda + 0.1), (long) (1.0*piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings, (long) (100.0*piaacmcsimul_var.PIAACMC_MASKRADLD+0.1), piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].nblambda);
	// fflush(stdout);
	// sleep(2);
	
	return EXIT_SUCCESS;
}


/** @brief Assemble focal plane mask name
 * 
 * 
 */ 
int PIAACMCsimul_update_fnamedescr()
{
	PIAACMCsimul_update_fnamedescr_conf();
	
	sprintf(piaacmcsimul_var.fnamedescr, "%s.minsag%06ld_maxsag%06ld_fpmregc%08ld_fpmrega%06ld_%s", piaacmcsimul_var.fnamedescr_conf, (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1.0e9*piaacmc[0].fpmsagreg_coeff+0.1),  (long) (1000.0*piaacmc[0].fpmsagreg_alpha+0.1), piaacmc[0].fpmmaterial_name);
	// debug from Justin
	// printf("Running PIAACMCsimul_update_fnamedescr: %s.minsag%06ld_maxsag%06ld_fpmregc%08ld_fpmrega%06ld_%s\n", piaacmcsimul_var.fnamedescr_conf, (long) (1.0e9*piaacmc[0].fpmminsag - 0.1), (long) (1.0e9*piaacmc[0].fpmmaxsag + 0.1), (long) (1.0e9*piaacmc[0].fpmsagreg_coeff+0.1),  (long) (1000.0*piaacmc[0].fpmsagreg_alpha+0.1), piaacmc[0].fpmmaterial_name);
	// fflush(stdout);
	// sleep(2);
		
	return EXIT_SUCCESS;
}





//
// Save PIAACMC optical design
//

int PIAACMCsimul_savepiaacmcconf(const char *dname)
{
    char command[1000];
    int r;
    FILE *fp;
    char fname[1500];
    long i;
    
	

	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif

    
    sprintf(command, "mkdir -p %s", dname);
    r = system(command);	  
    
    
    // piaacmcparam.conf
    
    sprintf(fname,"%s/piaacmcparams.conf", dname);

    fp = fopen(fname, "w");

    fprintf(fp, "%10ld   NBradpts\n", piaacmc[0].NBradpts);
    for(i=0; i<10; i++)
    {
        if(i<piaacmc[0].NBLyotStop)
        {
            sprintf(fname, "!%s/LyotStop%ld.fits", dname, i);
            if(piaacmc[0].IDLyotStop[i]!=-1)
                save_fits(data.image[piaacmc[0].IDLyotStop[i]].name, fname);
            fprintf(fp, "%10.6lf   LyotStop_zpos %ld\n", piaacmc[0].LyotStop_zpos[i], i);
        }
        else
            fprintf(fp, "%10.6lf   LyotStop_zpos %ld\n", piaacmc[0].LyotStop_zpos[i], i);
    }
    
    fprintf(fp, "%10.6f    fpmaskamptransm\n", piaacmc[0].fpmaskamptransm);

    fclose(fp);



	// PIAACMC optics

    sprintf(fname, "!%s/piaa0Cmodes.fits", dname);
    if(piaacmc[0].piaa0CmodesID!=-1)
        save_fits(data.image[piaacmc[0].piaa0CmodesID].name, fname);

    sprintf(fname, "!%s/piaa0Fmodes.fits", dname);
    if(piaacmc[0].piaa0FmodesID!=-1)
        save_fits(data.image[piaacmc[0].piaa0FmodesID].name, fname);

    sprintf(fname, "!%s/piaa1Cmodes.fits", dname);
    if(piaacmc[0].piaa1CmodesID!=-1)
        save_fits(data.image[piaacmc[0].piaa1CmodesID].name, fname);

    sprintf(fname, "!%s/piaa1Fmodes.fits", dname);
    if(piaacmc[0].piaa1FmodesID!=-1)
        save_fits(data.image[piaacmc[0].piaa1FmodesID].name, fname);





	// Focal plane mask
	PIAACMCsimul_update_fnamedescr();  
    
    sprintf(fname, "!%s/fpm_zonez.%s.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.fnamedescr);
    if(piaacmc[0].zonezID!=-1)
        save_fits(data.image[piaacmc[0].zonezID].name, fname);


    sprintf(fname, "!%s/fpm_zonea.%s.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.fnamedescr);
    if(piaacmc[0].zoneaID!=-1)
        save_fits(data.image[piaacmc[0].zoneaID].name, fname);

	
	long IDfpmzmap1;
	long ii;
	long xsize, ysize, xysize;

	IDfpmzmap1 = image_ID("fpmzmap1");
	if(IDfpmzmap1 == -1)
	{
		printf("Creating fpmzmap1 ...\n");
		fflush(stdout);
		IDfpmzmap1 = PIAACMCsimul_mkFPM_zonemap("fpmzmap1");
		xsize = data.image[IDfpmzmap1].md[0].size[0];
		ysize = data.image[IDfpmzmap1].md[0].size[1];
		xysize = xsize * ysize;
		
		for(ii=0;ii<xysize;ii++)
			data.image[IDfpmzmap1].array.UI16[ii] -= 1;
	}
	//list_image_ID();
	//printf("data.image[piaacmc[0].zonezID].name = %s\n", data.image[piaacmc[0].zonezID].name);

	
	image_basic_indexmap("fpmzmap1", data.image[piaacmc[0].zonezID].name, "fpmsagmapHR");	
    sprintf(fname, "!%s/fpm_sagmapHR.%s.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.fnamedescr);
    save_fits("fpmsagmapHR", fname);
	delete_image_ID("fpmsagmapHR");

    return(0);
}


