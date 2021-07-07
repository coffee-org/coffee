/**
 * @file    PIAACMCsimul_CA2propCubeInt.c
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



#include <stdlib.h>
#include <stdio.h>

#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"



/**
 * @brief Propagate complex amplitude image into intensity map cube
 *
 *
 */

errno_t PIAACMCsimul_CA2propCubeInt(
    const char *__restrict__ IDamp_name,
    const char *__restrict__ IDpha_name,
    float zmin,
    float zmax,
    long NBz,
    const char *__restrict__ IDout_name,
    imageID *outID
)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %s %s %f %f %ld %s",
                     IDamp_name, IDpha_name, zmin, zmax, NBz, IDout_name);

    imageID IDout;
    long nblambda;
    float *zarray;
    float zprop;

    uint32_t xsize;
    uint32_t ysize;
    uint64_t xysize;


    {
        imageID IDa = image_ID(IDamp_name);
        xsize = data.image[IDa].md[0].size[0];
        ysize = data.image[IDa].md[0].size[1];
        xysize = xsize;
        xysize *= ysize;

        create_2Dimage_ID("retmpim", xsize, ysize, NULL);
        create_2Dimage_ID("imtmpim", xsize, ysize, NULL);

        if(data.image[IDa].md[0].naxis == 3)
        {
            nblambda = data.image[IDa].md[0].size[2];
        }
        else
        {
            nblambda = 1;
        }
    }
    FUNC_CHECK_RETURN(
        create_3Dimage_ID(IDout_name, xsize, ysize, NBz, &IDout)
    );

    FUNC_CHECK_RETURN(
        create_2Dimage_ID("tmpintg", xsize, ysize, NULL)
    );


    // initialize zarray
    zarray = (float *) malloc(sizeof(float) * NBz);
    if(zarray == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    for(long l = 0; l < NBz; l++)
    {
        zarray[l] = zmin + (zmax - zmin) * l / (NBz - 1);
    }



    for(long l = 0; l < NBz; l++)
    {
        printf("l = %ld/%ld\n", l, NBz);
        fflush(stdout);

        zprop = zarray[l];
        FUNC_CHECK_RETURN(
            OptSystProp_propagateCube(
                &piaacmcopticalsystem,
                0,
                IDamp_name,
                IDpha_name,
                "_tmppropamp",
                "_tmpproppha",
                zprop,
                0
            )
        );

        imageID IDa = image_ID("_tmppropamp");
        // imageID IDp = image_ID("_tmpproppha");


        // write intensity
        for(long k = 0; k < nblambda; k++)
            for(long ii = 0; ii < xsize * ysize; ii++)
            {
                data.image[IDout].array.F[l*xysize + ii] +=
                    data.image[IDa].array.F[k*xysize+ii] * data.image[IDa].array.F[k*xysize+ii];
            }



        FUNC_CHECK_RETURN(
            delete_image_ID("_tmppropamp", DELETE_IMAGE_ERRMODE_WARNING)
        );

        FUNC_CHECK_RETURN(
            delete_image_ID("_tmpproppha", DELETE_IMAGE_ERRMODE_WARNING)
        );
    }

    free(zarray);
    FUNC_CHECK_RETURN(
        delete_image_ID("retmpim", DELETE_IMAGE_ERRMODE_WARNING)
    );

    FUNC_CHECK_RETURN(
        delete_image_ID("imtmpim", DELETE_IMAGE_ERRMODE_WARNING)
    );

    if(outID != NULL)
    {
        *outID = IDout;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

