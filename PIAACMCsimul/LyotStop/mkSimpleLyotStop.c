/**
 * @file    mkSimpleLyotStop.c
 * @brief   PIAA-type coronagraph design, run
 *
 */

#include <stdio.h>
#include <stdlib.h>

// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_memory/COREMOD_memory.h"

#include "PIAACMCsimul/PIAACMCsimul.h"

/**
 * @brief Make Lyot stop that transmits between rin and rout
 *
 * @param[out] ID_name  Output image name
 * @param[in]  rin      Inner Lyot mask radius
 * @param[in]  rout     Outer Lyot mask radius
 * @param[out] outID    Output image identifier
 *
 * @return errno_t
 */
errno_t mkSimpleLyotStop(const char *__restrict__ ID_name,
                         float    rin,
                         float    rout,
                         imageID *outID)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %s %f %f", ID_name, rin, rout);

    uint32_t size;
    uint64_t size2;
    imageID  ID, IDr;

    size  = piaacmcopticaldesign.size;
    size2 = size;
    size2 *= size;

    IDr = image_ID("rcoord");

    FUNC_CHECK_RETURN(create_3Dimage_ID(ID_name,
                                        size,
                                        size,
                                        piaacmcopticaldesign.nblambda,
                                        &ID));

    for (long k = 0; k < piaacmcopticaldesign.nblambda; k++)
        for (uint64_t ii = 0; ii < size2; ii++)
        {
            if ((data.image[IDr].array.F[ii] < rout) &&
                (data.image[IDr].array.F[ii] > rin))
            {
                data.image[ID].array.F[k * size2 + ii] = 1.0;
            }
            else
            {
                data.image[ID].array.F[k * size2 + ii] = 0.0;
            }
        }

    if (outID != NULL)
    {
        *outID = ID;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
