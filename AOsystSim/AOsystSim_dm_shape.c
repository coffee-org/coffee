/**
 * @file    AOsystSim_dm_shape.c
 * @brief   Sparse influence-function DM surface evaluation and DM daemon configuration template
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_memory/COREMOD_memory.h"

#include "AOsystSim/AOsystSim.h"

int    DMifpixarray_init = 0;
float *DMifpixarray_val;
long  *DMifpixarray_index;
long  *DMifpixarray_pixindex;
long   DMifpixarray_NBpix;
long   DMifpixarray_NBpix0;
long  *dmifpixactarray;
long  *dmifpixactarray_dmact;
long  *dmifpixactarray_ii;

/**
 * dmshape_init_cache - Build sparse lookup tables of non-zero pixels in DM influence-function cube
 * @IDdmifc: Influence-function 3D cube image ID
 * @npix:    Number of 2D surface pixels per slice
 * @DMnbact: Number of actuators (cube slices)
 */
static void dmshape_init_cache(
    imageID  IDdmifc,
    uint64_t npix,
    long     DMnbact)
{
    double eps = 1.0e-12;

    DMifpixarray_NBpix = 0;
    for (long dmact = 0; dmact < DMnbact; dmact++)
    {
        for (uint64_t ii = 0; ii < npix; ii++)
        {
            if (fabs(dcimg[IDdmifc].array.F[(uint64_t) dmact * npix + ii]) > eps)
            {
                DMifpixarray_NBpix++;
            }
        }
    }

    dmifpixactarray       = (long *) malloc(sizeof(long) * DMifpixarray_NBpix);
    dmifpixactarray_dmact = (long *) malloc(sizeof(long) * DMifpixarray_NBpix);
    dmifpixactarray_ii    = (long *) malloc(sizeof(long) * DMifpixarray_NBpix);
    if ((dmifpixactarray == NULL) || (dmifpixactarray_dmact == NULL) ||
        (dmifpixactarray_ii == NULL))
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    DMifpixarray_NBpix = 0;
    for (long dmact = 0; dmact < DMnbact; dmact++)
    {
        for (uint64_t ii = 0; ii < npix; ii++)
        {
            uint64_t idx = (uint64_t) dmact * npix + ii;

            if (fabs(dcimg[IDdmifc].array.F[idx]) > eps)
            {
                dmifpixactarray[DMifpixarray_NBpix]       = (long) idx;
                dmifpixactarray_dmact[DMifpixarray_NBpix] = dmact;
                dmifpixactarray_ii[DMifpixarray_NBpix]    = (long) ii;
                DMifpixarray_NBpix++;
            }
        }
    }

    if (DMifpixarray_val != NULL)
    {
        printf("WARNING: array DMifpixarray_val is not NULL");
        fflush(stdout);
    }
    DMifpixarray_val = (float *) malloc(sizeof(float) * DMifpixarray_NBpix);
    if (DMifpixarray_val == NULL)
    {
        printf("ERROR: could not allocate array DMifpixarray_val\n");
        exit(0);
    }

    DMifpixarray_index    = (long *) malloc(sizeof(long) * DMifpixarray_NBpix);
    DMifpixarray_pixindex = (long *) malloc(sizeof(long) * DMifpixarray_NBpix);
    if ((DMifpixarray_index == NULL) || (DMifpixarray_pixindex == NULL))
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    DMifpixarray_NBpix0 = DMifpixarray_NBpix;
    DMifpixarray_init   = 1;
}

/**
 * AOsystSim_DMshape - Convolve actuator command map with 3D influence-function cube into 2D OPD
 * @IDdmctrl_name: Input actuator command 2D image name
 * @IDdmifc_name:  Input 3D influence-function cube image name
 * @IDdm_name:     Output 2D DM OPD surface image name
 *
 * Return: 0 on completion.
 */
int AOsystSim_DMshape(
    const char *IDdmctrl_name,
    const char *IDdmifc_name,
    const char *IDdm_name)
{
    imageID  IDdmctrl = image_ID(IDdmctrl_name, dcimg, dcnimg);
    imageID  IDdmifc  = image_ID(IDdmifc_name, dcimg, dcnimg);
    uint32_t dmsizex  = dcimg[IDdmifc].md[0].size[0];
    uint32_t dmsizey  = dcimg[IDdmifc].md[0].size[1];
    long     DMnbact  = dcimg[IDdmifc].md[0].size[2];
    uint64_t npix     = (uint64_t) dmsizex * dmsizey;
    imageID  IDdm;

    if (DMifpixarray_init == 0)
    {
        dmshape_init_cache(IDdmifc, npix, DMnbact);
    }

    for (long kk = 0; kk < DMifpixarray_NBpix0; kk++)
    {
        DMifpixarray_val[kk]      = dcimg[IDdmifc].array.F[dmifpixactarray[kk]];
        DMifpixarray_index[kk]    = dmifpixactarray_dmact[kk];
        DMifpixarray_pixindex[kk] = dmifpixactarray_ii[kk];
    }

    IDdm = image_ID(IDdm_name, dcimg, dcnimg);
    if (IDdm == -1)
    {
        create_2Dimage_ID(IDdm_name, dmsizex, dmsizey, &IDdm);
    }

    for (uint64_t ii = 0; ii < npix; ii++)
    {
        dcimg[IDdm].array.F[ii] = 0.0f;
    }

    for (long kk = 0; kk < DMifpixarray_NBpix; kk++)
    {
        dcimg[IDdm].array.F[DMifpixarray_pixindex[kk]] +=
            dcimg[IDdmctrl].array.F[DMifpixarray_index[kk]] * DMifpixarray_val[kk];
    }

    return 0;
}

/**
 * AOsystSim_DM_mkCONF - Write template configuration file for deformable mirror stream daemon
 * @fname: Output configuration file path
 *
 * Return: 0 on completion.
 */
int AOsystSim_DM_mkCONF(const char *fname)
{
    FILE *fp = fopen(fname, "w");

    fprintf(fp, "\n");
    fprintf(fp, "# ============== AOsim: DM process =============================\n");
    fprintf(fp, "# DM process is triggered by WF OPD, not by DM command\n\n");
    fprintf(fp, "# ============== INPUT & TRIGGER (OPD unit = um) ===============\n");
    fprintf(fp, "INMODE                   0                 "
                "# 0:stream, 1:file system (FITS)\n");
    fprintf(fp, "INSTREAMNAMEDM           dm05disp          # input DM command stream\n");
    fprintf(fp, "INFITSFILENAMEDM         dm05disp.fits     "
                "# input FITS file name (DM command)\n");
    fprintf(fp, "INTRIGSTREAMNAME         wf1opd            # input trigger stream\n");
    fprintf(fp, "INTRIGSEMCHAN            6                 "
                "# input semaphore channel (using OPD input)\n");
    fprintf(fp, "INTRIGGERFILE            inwf.txt          # input trigger file\n\n");
    fprintf(fp, "# ============== OUTPUT TYPE ===================================\n");
    fprintf(fp, "OUTMODE                  0                 # 0: stream, 1: file system\n");
    fprintf(fp, "OUTSTREAMNAMEDM          dm05dispmap       # output WF stream name \n");
    fprintf(fp, "OUTFITSFILENAMEDM        dm05dispmap.fits  # FITS file name output [um]\n");
    fprintf(fp, "OUTTRIGGERFILE           outdm.txt         # output trigger file\n\n");
    fprintf(fp, "# ============== GEOMETRY, TIME =============================\n");
    fprintf(fp, "ARRAYSIZE                128               # array size, pix\n");
    fprintf(fp, "DMRAD                    52                # DM radius on array\n");
    fprintf(fp, "DMDT                     0.0003            # DM time sampling\n");
    fprintf(fp, "NBTSAMPLES               100               # number of time samples\n");
    fprintf(fp, "DMLAGSTART               0.0003            # time lag start\n");
    fprintf(fp, "DMTIMECST                0.0001            # DM time constant\n\n");
    fclose(fp);

    return 0;
}
