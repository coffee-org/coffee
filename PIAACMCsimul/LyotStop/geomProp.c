/**
 * @file    PIAACMCsimul_geomProp.c
 * @brief   Beam geometrical propagation based on local slopes
 *
 * Uses local slopes of sag map to propagate (un-corrected intensity only) a beam
 *
 */



#include <stdlib.h>
#include <stdio.h>
#include <math.h>



// milk includes
#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "image_filter/image_filter.h"

#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"



// Local variables pointers
static char *inintensity_imname;
static char *sag2D_imname;
static char *outintensity_imname;
static char *outcnt_imname;
static double *drindexval;
static double *pscaleval;
static double *zpropval;
static double *kradval;
static double *kstepval;
static double *rlimval;




static CLICMDARGDEF farg[] =
{
    {
        CLIARG_IMG, ".inintim", "input intensity image", "pupin",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &inintensity_imname, NULL
    },
    {
        CLIARG_IMG, ".insagim", "input 2D sag image", "piaa0z",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &sag2D_imname, NULL
    },
    {
        CLIARG_STR, ".outintim", "output intensity image", "pupout",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &outintensity_imname, NULL
    },
    {
        CLIARG_STR, ".outcntim", "output intensity image", "cntout",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &outcnt_imname, NULL
    },
    {
        CLIARG_FLOAT, ".drindex", "delta refractive index", "2.0",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &drindexval, NULL
    },
    {
        CLIARG_FLOAT, ".pscale", "pixel scale", "0.00011",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &pscaleval, NULL
    },
    {
        CLIARG_FLOAT, ".zprop", "propagation dist", "2.302606",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &zpropval, NULL
    },
    {
        CLIARG_FLOAT, ".krad", "kernel radius", "3.0",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &kradval, NULL
    },
    {
        CLIARG_FLOAT, ".kstep", "kernel step", "0.5",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &kstepval, NULL
    },
    {
        CLIARG_FLOAT, ".rlim", "clear aperture", "200.0",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &rlimval, NULL
    }
};

static CLICMDDATA CLIcmddata =
{
    "piaacmcgeomprop",
    "Geometric propagation from surface",
    CLICMD_FIELDS_DEFAULTS
};


// detailed help
static errno_t help_function()
{
    return RETURN_SUCCESS;
}



/**
 * @brief Lyot stops positions from zmin to zmax relative to current, working back (light goes from 0 to zmax)
 *
 * @param[in]	IDin_name		image : Input 2D intensity
 * @param[in]	IDsag_name	    image : 2D sag
 * @param[out]  IDout_name		image : 2D propagated intensity
 * @param[out]  IDoutcnt_name	image : 2D rays counter
 * @param[in]   drindex         float : refractive index (2 for mirror)
 * @param[in]	pscale			float : pixel scale [m]
 * @param[in]   zprop          	float : propagation distance [m]
 * @param[in]	krad			float : kernel radius used to evaluate slope [pixel]
 * @param[in]	kstep			float : step size in input pupil [pixel]
 * @param[in]	rlim			float : clear aperture radius (don't compute outside this value) [pixel]
 * @param[out]  outID           output ID
*/
errno_t PIAACMCsimul_geomProp(
    const char *__restrict__ IDin_name,
    const char *__restrict__ IDsag_name,
    const char *__restrict__ IDout_name,
    const char *__restrict__ IDoutcnt_name,
    double drindex,
    double pscale,
    double zprop,
    double krad,
    double kstep,
    double rlim,
    imageID *outID
)
{
    DEBUG_TRACE_FSTART();


    imageID IDin = image_ID(IDin_name);
    uint32_t xsize = data.image[IDin].md[0].size[0];
    uint32_t ysize = data.image[IDin].md[0].size[1];

    imageID IDsag = image_ID(IDsag_name);

    imageID IDout;
    FUNC_CHECK_RETURN(
        create_2Dimage_ID(IDout_name, xsize, ysize, &IDout)
    );

    imageID IDoutcnt;
    FUNC_CHECK_RETURN(
        create_2Dimage_ID(IDoutcnt_name, xsize, ysize, &IDoutcnt)
    );

    printf("kstep = %f\n", kstep);

    for (float x = 0.5*xsize-rlim; x < 0.5*xsize+rlim; x += kstep )
        for (float y = 0.5*ysize-rlim; y < 0.5*ysize+rlim; y += kstep )
        {
            long ii0 = (long) (x+0.5);
            long jj0 = (long) (y+0.5);

            double sumvx = 0.0;
            double sumvy = 0.0;
            double sumvx0 = 0.0;
            double sumvy0 = 0.0;
            double sumv = 0.0;
            double sumc = 0.0;

            long iimin, iimax, jjmin, jjmax;
            iimin = ii0 - ((long) krad+1);
            iimax = ii0 + ((long) krad+2);
            jjmin = jj0 - ((long) krad+1);
            jjmax = jj0 + ((long) krad+2);

            if(iimin<0) {
                PRINT_ERROR("iimin = %ld < 0  ii0 = %ld", iimin, ii0);
                abort();
            }
            if(iimax>xsize-1) {
                PRINT_ERROR("iimax = %ld > %ld  ii0 = %ld", iimax, (long) (xsize-1), ii0);
                abort();
            }

            if(jjmin<0) {
                PRINT_ERROR("jjmin = %ld < 0  jj0 = %ld", jjmin, jj0);
                abort();
            }
            if(jjmax>ysize-1) {
                PRINT_ERROR("jjmax = %ld > %ld  jj0 = %ld", jjmax, (long) (ysize-1), jj0);
                abort();
            }

            for (long ii = iimin; ii < iimax; ii++)
                for (long jj = jjmin; jj < jjmax; jj++)
                {
                    float dx = 1.0*ii - x;
                    float dy = 1.0*jj - y;

                    float dr = sqrt(dx*dx+dy*dy)/krad;
                    if(dr>1.0)
                        dr = 1.0;
                    float coeff = 0.5 + 0.5*cos(dr*M_PI);

                    sumv += coeff * data.image[IDsag].array.F[jj*xsize+ii];
                    sumvx += coeff * dx * data.image[IDsag].array.F[jj*xsize+ii];
                    sumvy += coeff * dy * data.image[IDsag].array.F[jj*xsize+ii];
                    sumvx0 += coeff * dx;
                    sumvy0 += coeff * dy;
                    sumc += coeff;
                }

            // local slopes [unitless]
            float slx = ( (sumvx / sumc) - (sumvx0 * (sumv/sumc) / sumc) ) / pscale;
            float sly = ( (sumvy / sumc) - (sumvy0 * (sumv/sumc) / sumc) ) / pscale;

            // displacements [pix]
            float dii = (slx * zprop * drindex) / pscale;
            float djj = (sly * zprop * drindex) / pscale;

            long ii1 = (long) ( x + dii + 0.5 );
            long jj1 = (long) ( y + djj + 0.5 );

            if((ii1>0)&&(ii1<xsize)&&(jj1>0)&&(jj1<ysize))
            {
                data.image[IDout].array.F[jj1*xsize + ii1] += data.image[IDin].array.F[jj0*xsize + ii0];
                data.image[IDoutcnt].array.F[jj1*xsize + ii1] += 1.0;
            }
        }

    for(uint32_t ii1 = 0; ii1 < xsize; ii1++)
        for(uint32_t jj1 = 0; jj1 < ysize; jj1++)
        {
            if (data.image[IDoutcnt].array.F[jj1*xsize + ii1] > 0.1)
                data.image[IDout].array.F[jj1*xsize + ii1] /= data.image[IDoutcnt].array.F[jj1*xsize + ii1];
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

    PIAACMCsimul_geomProp(
        inintensity_imname,
        sag2D_imname,
        outintensity_imname,
        outcnt_imname,
        *drindexval,
        *pscaleval,
        *zpropval,
        *kradval,
        *kstepval,
        *rlimval,
        NULL
    );

    INSERT_STD_PROCINFO_COMPUTEFUNC_END

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}




INSERT_STD_FPSCLIfunctions

// Register function in CLI
errno_t CLIADDCMD_PIAACMCsimul__LyotStop__PIAACMCsimul_geomProp()
{
    INSERT_STD_CLIREGISTERFUNC
    return RETURN_SUCCESS;
}
