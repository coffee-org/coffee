/**
 * @file    OptSystProp_propagateCube.c
 * @brief   Optical propagation of 2D complex amplitude to 3D cube
 *
 *
 */

#include "CommandLineInterface/CLIcore.h"



#include <math.h>
#include "COREMOD_memory/COREMOD_memory.h"
#include "WFpropagate/WFpropagate.h"
#include "OptSystProp/OptSystProp.h"








errno_t OptSystProp_propagateCube(
    OPTSYST *optsyst,
    long index,
    const char *IDin_amp_name,
    const char *IDin_pha_name,
    const char *IDout_amp_name,
    const char *IDout_pha_name,
    double zprop,
    int sharedmem
)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG prop %s %s -> %s %s, dist = %lf", IDin_amp_name, IDin_pha_name, IDout_amp_name, IDout_pha_name, zprop);

    printf("propagating by %lf m\n", zprop);

    imageID IDin_amp = image_ID(IDin_amp_name);
    imageID IDin_pha = image_ID(IDin_pha_name);
    uint32_t size = data.image[IDin_amp].md[0].size[0];
    uint64_t size2 = size * size;
    imageID IDc_in = create_2DCimage_ID("tmppropCin", size, size);

    uint32_t * imsizearray = (uint32_t *) malloc(sizeof(uint32_t) * 3);
    if(imsizearray == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }
    imsizearray[0] = size;
    imsizearray[1] = size;
    imsizearray[2] = optsyst[index].nblambda;

    imageID IDout_amp = image_ID(IDout_amp_name);
    if(IDout_amp == -1)
    {
        create_image_ID(IDout_amp_name, 3, imsizearray, _DATATYPE_FLOAT,
                        sharedmem, 0, 0, &IDout_amp);
    }

    imageID IDout_pha = image_ID(IDout_pha_name);
    if(IDout_pha == -1)
    {
        create_image_ID(IDout_pha_name, 3, imsizearray, _DATATYPE_FLOAT,
                        sharedmem, 0, 0, &IDout_pha);
    }
    free(imsizearray);

    data.image[IDout_amp].md[0].write = 1;
    data.image[IDout_pha].md[0].write = 1;

    for(int kl = 0; kl < optsyst[index].nblambda; kl++)
    {
        long IDc_out;

        printf("kl = %d / %d  %g\n", kl, optsyst[index].nblambda,
               optsyst[index].lambdaarray[kl]);
        // convert from amp/phase to Re/Im
        for(uint64_t ii = 0; ii < size2; ii++)
        {
            double amp = data.image[IDin_amp].array.F[kl * size2 + ii];
            double pha = data.image[IDin_pha].array.F[kl * size2 + ii];
            data.image[IDc_in].array.CF[ii].re = amp * cos(pha);
            data.image[IDc_in].array.CF[ii].im = amp * sin(pha);
        }
        // do the actual propagation
        Fresnel_propagate_wavefront("tmppropCin", "tmppropCout",
                                    optsyst[index].pixscale, zprop, optsyst[index].lambdaarray[kl]);

        IDc_out = image_ID("tmppropCout");
        // convert back from Re/Im to amp/phase
        for(uint64_t ii = 0; ii < size2; ii++)
        {
            double re = data.image[IDc_out].array.CF[ii].re;
            double im = data.image[IDc_out].array.CF[ii].im;
            double amp = sqrt(re * re + im * im);
            double pha = atan2(im, re);
            data.image[IDout_amp].array.F[kl * size2 + ii] = amp;
            data.image[IDout_pha].array.F[kl * size2 + ii] = pha;
        }
        delete_image_ID("tmppropCout", DELETE_IMAGE_ERRMODE_WARNING);
    }

    data.image[IDout_amp].md[0].cnt0++;
    data.image[IDout_pha].md[0].cnt0++;

    data.image[IDout_amp].md[0].write = 0;
    data.image[IDout_pha].md[0].write = 0;

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

