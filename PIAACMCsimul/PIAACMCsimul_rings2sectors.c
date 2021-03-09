/**
 * @file    PIAAACMCsimul_rings2sectors.c
 * @brief   PIAA-type coronagraph design, rings to sectors
 * 
 * Can design both APLCMC and PIAACMC coronagraphs
 *  
 * @author  O. Guyon
 * @date    24 nov 2017
 *
 * 
 * | date        |  Code Change      |
 * |-------------|-------------------|
 * | 2017-11-24  | documentation     |
 * 
 * 
 * @bug No known bugs.
 * 
 */



#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "PIAACMCsimul/PIAACMCsimul.h"






/**
 * @brief Rings to sectors
 * 
 * @param[in] IDin_name	input image: circular mask design
 * @param[in] sectfname	text file specifying which zones belong to which rings
 * @param[out] IDout_name	output sector mask design
 * 
 */ 
long PIAACMCsimul_rings2sectors(
		const char *IDin_name, 
		const char *sectfname, 
		const char *IDout_name
		)
{
    long IDin, IDout;
    FILE *fp;
    long nbring, nbzone;
    long tmpl1, tmpl2;
    long zone;
    long arrayring[5000];

	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif


    IDin = image_ID(IDin_name);
    nbring = data.image[IDin].md[0].size[0];

    nbzone = 0;
    nbring = 0;
    fp = fopen(sectfname,"r");
    while(fscanf(fp,"%ld %ld\n", &tmpl1, &tmpl2)==2)
    {
        arrayring[tmpl1] = tmpl2;
        if(tmpl2>nbring)
            nbring = tmpl2;
        if(tmpl1>nbzone)
            nbzone = tmpl1;
    }
    fclose(fp);
    nbring++;
    nbzone++;

    IDout = create_2Dimage_ID_double(IDout_name, nbzone, 1);
    for(zone=0; zone<nbzone; zone++)
        data.image[IDout].array.D[zone] = data.image[IDin].array.D[arrayring[zone]];

    printf("%ld zones in %ld rings\n", nbzone, nbring);

    return(IDout);
}



