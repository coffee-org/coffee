/**
 * @file    PIAAACMCsimul_rings2sectors.c
 * @brief   PIAA-type coronagraph design, rings to sectors
 *
 * Can design both APLCMC and PIAACMC coronagraphs
 *
 */



#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "PIAACMCsimul/PIAACMCsimul.h"






// Local variables pointers
static char *inimname;
static char *secfname;
static char *outimname;



static CLICMDARGDEF farg[] =
{
    {
        CLIARG_IMG, ".inimname", "input image: circular mask design", "imin",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &inimname
    },
    {
        CLIARG_STR, ".secfname", "text file specifying which zones belong to which rings", "sec.txt",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &secfname
    },
    {
        CLIARG_STR_NOT_IMG, ".outimname", "output sector mask design", "outim",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &outimname
    }
};

static CLICMDDATA CLIcmddata =
{
    "ring2sect",
    "turn ring fpm design into sectors",
    CLICMD_FIELDS_DEFAULTS
};


// detailed help
static errno_t help_function()
{
    return RETURN_SUCCESS;
}





/**
 * @brief Rings to sectors
 *
 * @param[in] IDin_name	input image: circular mask design
 * @param[in] sectfname	text file specifying which zones belong to which rings
 * @param[out] IDout_name	output sector mask design
 *
 */
errno_t rings2sectors(
    const char *IDin_name,
    const char *sectfname,
    const char *IDout_name,
    imageID    *outID
)
{
    DEBUG_TRACE_FSTART();

    imageID IDin, IDout;
    FILE *fp;
    long nbring, nbzone;
    long tmpl1, tmpl2;
    long zone;
    long arrayring[5000];


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

    FUNC_CHECK_RETURN(
        create_2Dimage_ID_double(IDout_name, nbzone, 1, &IDout);
    );

    for(zone=0; zone<nbzone; zone++)
        data.image[IDout].array.D[zone] = data.image[IDin].array.D[arrayring[zone]];

    printf("%ld zones in %ld rings\n", nbzone, nbring);

    if(outID != NULL)
    {
        *outID = IDout;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}





static errno_t compute_function()
{


    INSERT_STD_PROCINFO_COMPUTEFUNC_START

    rings2sectors(
        inimname,
        secfname,
        outimname,
        NULL
    );

    INSERT_STD_PROCINFO_COMPUTEFUNC_END

    return RETURN_SUCCESS;
}




INSERT_STD_FPSCLIfunctions

// Register function in CLI
errno_t CLIADDCMD_PIAACMCsimul__ring2sectors()
{
    INSERT_STD_CLIREGISTERFUNC
    return RETURN_SUCCESS;
}


