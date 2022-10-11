/**
 * @file    PIAACMCsimul_FPMresp_resample.c
 * @brief   PIAA-type coronagraph design
 *
 * Can design both APLCMC and PIAACMC coronagraphs
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>

// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_memory/COREMOD_memory.h"

// Local variables pointers
static char *FPMrespinimname;
static char *FPMrespoutimname;
static long *NBlambdaval;
static long *PTstepval;

static CLICMDARGDEF farg[] = {{
        CLIARG_IMG,
        ".FPMrespin",
        "input FPM response image",
        "FPMresp",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &FPMrespinimname,
        NULL
    },
    {
        CLIARG_IMG,
        ".FPMrespout",
        "output FPM response image",
        "FPMrespout",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &FPMrespoutimname,
        NULL
    },
    {
        CLIARG_LONG,
        ".NBlambda",
        "",
        "10",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &NBlambdaval,
        NULL
    },
    {
        CLIARG_LONG,
        ".PTstep",
        "EvalPts step",
        "2",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &PTstepval,
        NULL
    }
};

static CLICMDDATA CLIcmddata =
{
    "fpmresprs", "resample FPM resp matrix", CLICMD_FIELDS_DEFAULTS
};

// detailed help
static errno_t help_function()
{
    return RETURN_SUCCESS;
}

// compress FPMresp to smaller number of lambda and points
errno_t PIAACMC_FPMresp_resample(const char *__restrict__ FPMresp_in_name,
                                 const char *__restrict__ FPMresp_out_name,
                                 long     NBlambda,
                                 long     PTstep,
                                 imageID *outID)
{
    DEBUG_TRACE_FSTART();

    imageID ID    = -1;
    imageID IDout = -1;

    ID             = image_ID(FPMresp_in_name);
    uint32_t xsize = data.image[ID].md[0].size[0];
    uint32_t ysize = data.image[ID].md[0].size[1];
    uint32_t zsize = data.image[ID].md[0].size[2];

    long xsize1 = (long)(xsize / PTstep);
    long zsize1 = NBlambda;

    FUNC_CHECK_RETURN(create_3Dimage_ID_double(FPMresp_out_name,
                      xsize1,
                      ysize,
                      zsize1,
                      &IDout));
    /*	for(kk1=0;kk1<zsize1;kk1++)
    		for(ii=0;ii<xsize;ii++)
    			data.image[IDout].array.D[kk1*xsize1*ysize + ii] = 1;
    	*/

    for(long kk1 = 0; kk1 < zsize1; kk1++)
    {
        double kk1x  = 1.0 * kk1 / (zsize1 - 1);
        long   kk1xi = (long)(kk1x * (zsize - 1));
        double alpha = kk1x * (zsize - 1) - kk1xi;

        if(kk1xi == zsize - 1)
        {
            kk1xi = zsize - 2;
            alpha = 1.0;
        }

        printf("lambda index %ld  (%lf)   :  %ld %lf (%lf)\n",
               kk1,
               kk1x,
               kk1xi,
               alpha,
               (kk1xi + alpha) / (zsize - 1));
        long ii1 = 0;

        for(uint32_t ii = 0; ii < xsize; ii += 2 * PTstep)
        {

            for(uint32_t jj = 0; jj < ysize; jj++)
            {
                double re0 =
                    data.image[ID]
                    .array.D[kk1xi * xsize * ysize + jj * xsize + ii];
                double im0 =
                    data.image[ID]
                    .array.D[kk1xi * xsize * ysize + jj * xsize + ii + 1];

                double re1 =
                    data.image[ID]
                    .array.D[(kk1xi + 1) * xsize * ysize + jj * xsize + ii];
                double im1 =
                    data.image[ID].array.D[(kk1xi + 1) * xsize * ysize +
                                           jj * xsize + ii + 1];

                double re = (1.0 - alpha) * re0 + alpha * re1;
                double im = (1.0 - alpha) * im0 + alpha * im1;

                data.image[IDout]
                .array.D[kk1 * xsize1 * ysize + jj * xsize1 + ii1] =
                    re; //*sqrt(PTstep);
                data.image[IDout]
                .array.D[kk1 * xsize1 * ysize + jj * xsize1 + ii1 + 1] =
                    im; //*sqrt(PTstep);
            }

            ii1 += 2;
        }
    }

    /*		for(kk=0;kk<zsize;kk++)
    			{
    				for(jj=0;jj<ysize1;jj++)
    					data.image[IDout].array.D[kk*xsize*ysize1 + jj*xsize + ii] = data.image[ID].array.D[kk*xsize*ysize + jj*xsize + ii];
    				for(jj=ysize1;jj<ysize;jj++)
    					data.image[IDout].array.D[kk*xsize*ysize1 + ii] += data.image[ID].array.D[kk*xsize*ysize + jj*xsize + ii];
    			}
    */

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

    PIAACMC_FPMresp_resample(FPMrespinimname,
                             FPMrespoutimname,
                             *NBlambdaval,
                             *PTstepval,
                             NULL);

    INSERT_STD_PROCINFO_COMPUTEFUNC_END

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

INSERT_STD_FPSCLIfunctions

// Register function in CLI
errno_t
CLIADDCMD_PIAACMCsimul__FPMresp_resample()
{
    INSERT_STD_CLIREGISTERFUNC
    return RETURN_SUCCESS;
}
