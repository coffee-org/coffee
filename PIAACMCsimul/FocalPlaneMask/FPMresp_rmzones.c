/**
 * @file    PIAACMCsimul_PIAACMC_FPMresp_rmzones.c
 * @brief   PIAA-type coronagraph design
 *
 * Can design both APLCMC and PIAACMC coronagraphs
 *
 *
 */



#include <stdlib.h>
#include <stdio.h>


// milk includes
#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"



// Local variables pointers
static char *FPMrespinimname;
static char *FPMrespoutimname;
static long *NBzrmval;


static CLICMDARGDEF farg[] =
{
    {
        CLIARG_IMG, ".FPMrespin", "input FPM response image", "FPMresp",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &FPMrespinimname, NULL
    },
    {
        CLIARG_IMG, ".FPMrespout", "output FPM response image", "FPMrespout",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &FPMrespoutimname, NULL
    },
    {
        CLIARG_LONG, ".NBzrm", "NBzone removed", "125",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &NBzrmval, NULL
    }
};


static CLICMDDATA CLIcmddata =
{
    "fpmresprm",
    "remove outer zones in FPM resp matrix",
    CLICMD_FIELDS_DEFAULTS
};


// detailed help
static errno_t help_function()
{
    return RETURN_SUCCESS;
}






// remove outer zones to FPMresp
errno_t FPMresp_rmzones(
    const char *__restrict__ FPMresp_in_name,
    const char *__restrict__ FPMresp_out_name,
    long NBzones,
    imageID *outID
)
{
    DEBUG_TRACE_FSTART();

    imageID ID, IDout;


    ID = image_ID(FPMresp_in_name);
    uint32_t xsize = data.image[ID].md[0].size[0];
    uint32_t ysize = data.image[ID].md[0].size[1];
    uint32_t zsize = data.image[ID].md[0].size[2];


    long ysize1 = data.image[ID].md[0].size[1]-NBzones;

    FUNC_CHECK_RETURN(
        create_3Dimage_ID_double(FPMresp_out_name, xsize, ysize1, zsize, &IDout)
    );

    for(uint32_t ii=0; ii<xsize; ii++)
        for(uint32_t kk=0; kk<zsize; kk++)
        {
            for(long jj=0; jj<ysize1; jj++)
                data.image[IDout].array.D[kk*xsize*ysize1 + jj*xsize + ii] = data.image[ID].array.D[kk*xsize*ysize + jj*xsize + ii];
            for(long jj=ysize1; jj<ysize; jj++)
                data.image[IDout].array.D[kk*xsize*ysize1 + ii] += data.image[ID].array.D[kk*xsize*ysize + jj*xsize + ii];
        }

    if(outID != NULL)
    {
        *outID = IDout;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}



static errno_t compute_function()
{
    DEBUG_TRACE_FSTART();


    INSERT_STD_PROCINFO_COMPUTEFUNC_START

    FPMresp_rmzones(
        FPMrespinimname,
        FPMrespoutimname,
        *NBzrmval,
        NULL
    );

    INSERT_STD_PROCINFO_COMPUTEFUNC_END

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}




INSERT_STD_FPSCLIfunctions

// Register function in CLI
errno_t CLIADDCMD_PIAACMCsimul__FPMresp_rmzones()
{
    INSERT_STD_CLIREGISTERFUNC
    return RETURN_SUCCESS;
}

