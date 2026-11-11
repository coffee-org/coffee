/**
 * @file    AOsystSim_wfs_pyramid.c
 * @brief   Modulated Pyramid WFS intensity image kernel and shared-memory instance runner
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "fft/fft.h"
#include "image_filter/image_filter.h"

#include "AOsystSim/AOsystSim.h"
#include "AOsystSim_internal.h"

/**
 * pyrwfs_sim_init_masks - Build and cache focal-plane pyramid amplitude and phase modulation masks
 * @arraysize: Square wavefront grid size [pix]
 * @modampl:   Modulation radius [pix]
 * @modnbpts:  Number of discrete modulation points around modulation circle
 */
static void pyrwfs_sim_init_masks(
    uint32_t arraysize,
    double   modampl,
    long     modnbpts)
{
    double   lenssize  = 0.4 * arraysize;
    double   pcoeff    = 1.4;
    uint32_t imsize[2] = { arraysize, arraysize };
    char     pnamea[200];
    char     pnamep[200];
    char     pfnamea[200];
    char     pfnamep[200];

    for (long pmodpt = 0; pmodpt < modnbpts; pmodpt++)
    {
        imageID IDpyramp;
        imageID IDpyrpha;
        double  PA;
        double  xc;
        double  yc;

        sprintf(pnamea, "pyramp_%03ld", pmodpt);
        sprintf(pnamep, "pyrpha_%03ld", pmodpt);
        IDpyramp = image_ID(pnamea, dcimg, dcnimg);
        IDpyrpha = image_ID(pnamep, dcimg, dcnimg);
        if ((IDpyramp != -1) && (IDpyrpha != -1))
        {
            continue;
        }

        create_image_ID(pnamea, 2, imsize, _DATATYPE_FLOAT, 0, 0, 0, &IDpyramp);
        create_image_ID("pyrpha0", 2, imsize, _DATATYPE_FLOAT, 0, 0, 0, &IDpyrpha);

        PA = 2.0 * M_PI * pmodpt / modnbpts;
        xc = modampl * cos(PA);
        yc = modampl * sin(PA);

        for (uint32_t ii = 0; ii < arraysize; ii++)
        {
            for (uint32_t jj = 0; jj < arraysize; jj++)
            {
                double x = 1.0 * ((long) ii - (long) (arraysize / 2)) - xc;
                double y = 1.0 * ((long) jj - (long) (arraysize / 2)) - yc;

                dcimg[IDpyrpha].array.F[jj * arraysize + ii] =
                    (float) (pcoeff * (fabs(x) + fabs(y)));
                if ((fabs(x) > lenssize) || (fabs(y) > lenssize))
                {
                    dcimg[IDpyramp].array.F[jj * arraysize + ii] = 0.0f;
                }
                else
                {
                    dcimg[IDpyramp].array.F[jj * arraysize + ii] = 1.0f;
                }
            }
        }

        gauss_filter("pyrpha0", pnamep, 1.0, 10);
        delete_image_ID("pyrpha0", DELETE_IMAGE_ERRMODE_WARNING);

        sprintf(pfnamea, "AOsystSim_wdir/pyramp_%03ld.fits", pmodpt);
        sprintf(pfnamep, "AOsystSim_wdir/pyrpha_%03ld.fits", pmodpt);
        printf("SAVING: %s -> %s\n", pnamea, pfnamea);
        save_fits(pnamea, pfnamea);
        save_fits(pnamep, pfnamep);
    }
}

/**
 * pyrwfs_sim_init_buffers - Ensure output and scratch Pyramid WFS images exist
 * @outWFSim_name: Output WFS intensity stream name
 * @arraysize:     Square wavefront grid size [pix]
 * @ID_out:        Output WFS intensity image ID
 * @ID_tmp:        Output scratch accumulation image ID
 * @ID_cp:         Output scratch complex pupil copy image ID
 */
static void pyrwfs_sim_init_buffers(
    const char *outWFSim_name,
    uint32_t    arraysize,
    imageID    *ID_out,
    imageID    *ID_tmp,
    imageID    *ID_cp)
{
    uint32_t imsize[2] = { arraysize, arraysize };

    *ID_out = image_ID(outWFSim_name, dcimg, dcnimg);
    if (*ID_out == -1)
    {
        create_image_ID(outWFSim_name, 2, imsize, _DATATYPE_FLOAT, 1, 0, 0, ID_out);
        COREMOD_MEMORY_image_set_createsem(outWFSim_name, 7);
    }

    *ID_tmp = image_ID("outpwfsimtmp", dcimg, dcnimg);
    if (*ID_tmp == -1)
    {
        create_image_ID("outpwfsimtmp", 2, imsize, _DATATYPE_FLOAT, 1, 0, 0, ID_tmp);
    }

    *ID_cp = image_ID("pyrwfcin", dcimg, dcnimg);
    if (*ID_cp == -1)
    {
        create_image_ID("pyrwfcin", 2, imsize, _DATATYPE_COMPLEX_FLOAT, 0, 0, 0, ID_cp);
    }
}

/**
 * pyrwfs_sim_mod_step - Propagate one modulation angle through the Pyramid focal-plane mask
 * @pmodpt:    Modulation point index
 * @modnbpts:  Total modulation point count
 * @arraysize: Square wavefront grid size [pix]
 * @ID_inWFc:  Input complex pupil wavefront image ID
 * @ID_cp:     Scratch complex pupil image ID
 * @ID_tmp:    Scratch intensity accumulator image ID
 */
static void pyrwfs_sim_mod_step(
    long     pmodpt,
    long     modnbpts,
    uint32_t arraysize,
    imageID  ID_inWFc,
    imageID  ID_cp,
    imageID  ID_tmp)
{
    uint64_t arraysize2 = (uint64_t) arraysize * arraysize;
    char     pnamea[200];
    char     pnamep[200];
    imageID  IDpyramp;
    imageID  IDpyrpha;
    imageID  IDa;
    imageID  IDp;

    memcpy(dcimg[ID_cp].array.CF, dcimg[ID_inWFc].array.CF, sizeof(complex_float) * arraysize2);

    permut("pyrwfcin");
    do2dfft("pyrwfcin", "pyrpsfcin");
    permut("pyrpsfcin");
    mk_amph_from_complex("pyrpsfcin", "pyrpsfa", "pyrpsfp", 0);
    delete_image_ID("pyrpsfcin", DELETE_IMAGE_ERRMODE_WARNING);

    sprintf(pnamea, "pyramp_%03ld", pmodpt);
    sprintf(pnamep, "pyrpha_%03ld", pmodpt);
    IDpyramp = image_ID(pnamea, dcimg, dcnimg);
    IDpyrpha = image_ID(pnamep, dcimg, dcnimg);
    IDa      = image_ID("pyrpsfa", dcimg, dcnimg);
    IDp      = image_ID("pyrpsfp", dcimg, dcnimg);

    for (uint64_t ii = 0; ii < arraysize2; ii++)
    {
        dcimg[IDa].array.F[ii] *= dcimg[IDpyramp].array.F[ii];
        dcimg[IDp].array.F[ii] += dcimg[IDpyrpha].array.F[ii];
    }

    mk_complex_from_amph("pyrpsfa", "pyrpsfp", "pyrpsfc", 0);
    delete_image_ID("pyrpsfa", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pyrpsfp", DELETE_IMAGE_ERRMODE_WARNING);

    permut("pyrpsfc");
    do2dfft("pyrpsfc", "pyrwfs_pupc");
    delete_image_ID("pyrpsfc", DELETE_IMAGE_ERRMODE_WARNING);
    permut("pyrwfs_pupc");
    mk_amph_from_complex("pyrwfs_pupc", "pyrwfs_pupa", "pyrwfs_pupp", 0);

    delete_image_ID("pyrwfs_pupp", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("pyrwfs_pupc", DELETE_IMAGE_ERRMODE_WARNING);

    IDa = image_ID("pyrwfs_pupa", dcimg, dcnimg);
    for (uint64_t ii = 0; ii < arraysize2; ii++)
    {
        dcimg[ID_tmp].array.F[ii] +=
            dcimg[IDa].array.F[ii] * dcimg[IDa].array.F[ii] / (float) modnbpts;
    }
    delete_image_ID("pyrwfs_pupa", DELETE_IMAGE_ERRMODE_WARNING);
}

/**
 * AOsystSim_WFSsim_Pyramid - Compute modulated Pyramid WFS intensity image from complex pupil WF
 * @inWFc_name:    Input complex pupil wavefront image name
 * @outWFSim_name: Output Pyramid WFS intensity image name
 * @modampl:       Modulation amplitude [pix]
 * @modnbpts:      Number of modulation points
 *
 * Return: 0 on completion.
 */
int AOsystSim_WFSsim_Pyramid(
    const char *inWFc_name,
    const char *outWFSim_name,
    double      modampl,
    long        modnbpts)
{
    imageID  ID_inWFc   = image_ID(inWFc_name, dcimg, dcnimg);
    uint32_t arraysize  = dcimg[ID_inWFc].md[0].size[0];
    uint64_t arraysize2 = (uint64_t) arraysize * arraysize;
    imageID  ID_outWFSim;
    imageID  ID_outWFSim_tmp;
    imageID  ID_inWFccp;

    aosim_ensure_wdir();
    pyrwfs_sim_init_masks(arraysize, modampl, modnbpts);
    pyrwfs_sim_init_buffers(outWFSim_name, arraysize, &ID_outWFSim, &ID_outWFSim_tmp, &ID_inWFccp);

    dcimg[ID_outWFSim].md[0].write = 1;
    for (uint64_t ii = 0; ii < arraysize2; ii++)
    {
        dcimg[ID_outWFSim_tmp].array.F[ii] = 0.0f;
    }

    for (long pmodpt = 0; pmodpt < modnbpts; pmodpt++)
    {
        pyrwfs_sim_mod_step(pmodpt, modnbpts, arraysize, ID_inWFc, ID_inWFccp, ID_outWFSim_tmp);
    }

    memcpy(dcimg[ID_outWFSim].array.F, dcimg[ID_outWFSim_tmp].array.F, sizeof(float) * arraysize2);
    dcimg[ID_outWFSim].md[0].cnt0++;
    dcimg[ID_outWFSim].md[0].write = 0;

    return 0;
}

/**
 * AOsystSim_runWFS - Poll optical-element WF amplitude/phase streams and publish Pyramid WFS image
 * @index:      Optical element plane index
 * @IDout_name: Output WFS intensity stream name
 *
 * Return: 0 on completion.
 */
int AOsystSim_runWFS(
    long        index,
    const char *IDout_name)
{
    uint64_t cnt0 = 0;
    imageID  IDinamp;
    char     imnameamp[STRINGMAXLEN_IMGNAME];
    char     imnamepha[STRINGMAXLEN_IMGNAME];

    WRITE_IMAGENAME(imnameamp, "WFamp0_%03ld", index);
    WRITE_IMAGENAME(imnamepha, "WFpha0_%03ld", index);
    IDinamp = image_ID(imnameamp, dcimg, dcnimg);

    while (1)
    {
        while (cnt0 == dcimg[IDinamp].md[0].cnt0)
        {
            usleep(50);
        }
        cnt0 = dcimg[IDinamp].md[0].cnt0;

        mk_complex_from_amph(imnameamp, imnamepha, "_tmpwfc", 0);
        AOsystSim_WFSsim_Pyramid("_tmpwfc", IDout_name, 0.0, 1);
        delete_image_ID("_tmpwfc", DELETE_IMAGE_ERRMODE_WARNING);
    }

    return 0;
}
