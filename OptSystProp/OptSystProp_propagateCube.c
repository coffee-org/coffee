/**
 * @file    OptSystProp_propagateCube.c
 * @brief   Optical propagation of 2D complex amplitude to 3D cube
 *
 *
 * @author  O. Guyon
 * @date    24 nov 2017
 *
 *
 * | date        |  Code Change                                       |
 * |-------------|----------------------------------------------------|
 * | 2017-11-24  | File creation (broken off from original main .c)   |
 *
 *
 * @bug No known bugs.
 *
 */


#include <math.h>

#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "WFpropagate/WFpropagate.h"
#include "OptSystProp/OptSystProp.h"








int OptSystProp_propagateCube(
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
    int kl;
    long ii;
    long size;
    long size2;
    long IDin_amp, IDin_pha;
    long IDc_in;
    long IDout_amp, IDout_pha;
    uint32_t *imsizearray;
    double amp, pha, re, im;



    printf("propagating by %lf m\n", zprop);

    IDin_amp = image_ID(IDin_amp_name);
    IDin_pha = image_ID(IDin_pha_name);
    size = data.image[IDin_amp].md[0].size[0];
    size2 = size * size;
    IDc_in = create_2DCimage_ID("tmppropCin", size, size);



    imsizearray = (uint32_t *) malloc(sizeof(uint32_t) * 3);
    imsizearray[0] = size;
    imsizearray[1] = size;
    imsizearray[2] = optsyst[index].nblambda;

    IDout_amp = image_ID(IDout_amp_name);

    if(IDout_amp == -1)
    {
        IDout_amp = create_image_ID(IDout_amp_name, 3, imsizearray, _DATATYPE_FLOAT,
                                    sharedmem, 0);
    }


    IDout_pha = image_ID(IDout_pha_name);
    if(IDout_pha == -1)
    {
        IDout_pha = create_image_ID(IDout_pha_name, 3, imsizearray, _DATATYPE_FLOAT,
                                    sharedmem, 0);
    }


    data.image[IDout_amp].md[0].write = 1;
    data.image[IDout_pha].md[0].write = 1;

    for(kl = 0; kl < optsyst[index].nblambda; kl++)
    {
        long IDc_out;

        printf("kl = %d / %d  %g\n", kl, optsyst[index].nblambda,
               optsyst[index].lambdaarray[kl]);
        // convert from amp/phase to Re/Im
        for(ii = 0; ii < size2; ii++)
        {
            amp = data.image[IDin_amp].array.F[kl * size2 + ii];
            pha = data.image[IDin_pha].array.F[kl * size2 + ii];
            data.image[IDc_in].array.CF[ii].re = amp * cos(pha);
            data.image[IDc_in].array.CF[ii].im = amp * sin(pha);
        }
        // do the actual propagation
        Fresnel_propagate_wavefront("tmppropCin", "tmppropCout",
                                    optsyst[index].pixscale, zprop, optsyst[index].lambdaarray[kl]);

        IDc_out = image_ID("tmppropCout");
        // convert back from Re/Im to amp/phase
        for(ii = 0; ii < size2; ii++)
        {
            re = data.image[IDc_out].array.CF[ii].re;
            im = data.image[IDc_out].array.CF[ii].im;
            amp = sqrt(re * re + im * im);
            pha = atan2(im, re);
            data.image[IDout_amp].array.F[kl * size2 + ii] = amp;
            data.image[IDout_pha].array.F[kl * size2 + ii] = pha;
        }
        delete_image_ID("tmppropCout");
    }

    data.image[IDout_amp].md[0].cnt0++;
    data.image[IDout_pha].md[0].cnt0++;

    data.image[IDout_amp].md[0].write = 0;
    data.image[IDout_pha].md[0].write = 0;


    return 0;
}

