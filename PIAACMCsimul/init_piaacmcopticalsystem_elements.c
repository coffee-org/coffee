/**
 * @file    init_piaacmcopticalsystem_elements.c
 * @brief   PIAA optics, focal plane mask, and Lyot stop element builders for OPTSYST
 */

#include <stdio.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "image_gen/image_gen.h"

#include "OptSystProp/OptSystProp.h"

#include "FocalPlaneMask/mkFocalPlaneMask.h"
#include "PIAACMCsimul.h"
#include "init_piaacmcopticalsystem_internal.h"

/**
 * add_pre_piaa0_mask - Append optional upstream pupil-plane apodizer mask if present.
 * @elem: Current element index cursor (incremented when element is added)
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t add_pre_piaa0_mask(long *elem)
{
    imageID ID = image_ID("prePIAA0mask", dcimg, dcnimg);
    if (ID == -1)
    {
        FUNC_CHECK_RETURN(
            load_fits("prePIAA0mask.fits", "prePIAA0mask", LOADFITS_ERRMODE_WARNING, &ID));
    }

    if (ID != -1)
    {
        long idx = *elem;

        piaacmcopticaldesign.prePIAA0mask = 1;
        snprintf(
            piaacmcopticalsystem.name[idx],
            STRINGMAXLEN_OPTSYST_ELEMNAME,
            "pupil plane apodizer");
        piaacmcopticalsystem.elemtype[idx]       = 1;
        piaacmcopticalsystem.elemarrayindex[idx] = ID;
        piaacmcopticalsystem.elemZpos[idx]       = piaacmcopticaldesign.prePIAA0maskpos;
        (*elem)++;
    }
    else
    {
        piaacmcopticaldesign.prePIAA0mask = 0;
    }

    return RETURN_SUCCESS;
}

/**
 * add_post_piaa0_mask - Append optional opaque mask downstream of PIAA element 0 if on disk.
 * @elem: Current element index cursor (incremented when element is added)
 */
static void add_post_piaa0_mask(long *elem)
{
    imageID ID = -1;
    load_fits("postPIAA0mask.fits", "postPIAA0mask", LOADFITS_ERRMODE_IGNORE, &ID);

    if (ID != -1)
    {
        long idx = *elem;

        piaacmcopticaldesign.postPIAA0mask = 1;
        snprintf(
            piaacmcopticalsystem.name[idx],
            STRINGMAXLEN_OPTSYST_ELEMNAME,
            "opaque mask after PIAA element 0");
        piaacmcopticalsystem.elemtype[idx]       = 1;
        piaacmcopticalsystem.elemarrayindex[idx] = ID;
        piaacmcopticalsystem.elemZpos[idx]       = piaacmcopticaldesign.postPIAA0maskpos;
        (*elem)++;
    }
    else
    {
        piaacmcopticaldesign.postPIAA0mask = 0;
    }
}

/**
 * bind_piaa_surface - Bind reflective or refractive sag map to an OPTSYST PIAA element.
 * @elem:     OPTSYST element index
 * @surfslot: Slot in ASPHSURFMarray or ASPHSURFRarray (1 or 2)
 * @mname:    Reflective mirror sag image name
 * @rname:    Refractive lens sag image name
 * @piaa_idx: PIAA element index (0 or 1) for error reporting
 *
 * Return: RETURN_SUCCESS on success, error code if sag surface image is missing.
 */
static errno_t bind_piaa_surface(
    long        elem,
    long        surfslot,
    const char *mname,
    const char *rname,
    int         piaa_idx)
{
    if (piaacmcopticaldesign.PIAAmaterial_code == 0)
    {
        DEBUG_TRACEPOINT("Element %ld is mirror", elem);
        imageID IDm = image_ID(mname, dcimg, dcnimg);
        piaacmcopticalsystem.ASPHSURFMarray[surfslot].surfID = IDm;
        if (IDm == -1)
        {
            if (piaa_idx == 0)
            {
                FUNC_RETURN_FAILURE("PIAA mirror surface 0 not identified");
            }
            FUNC_RETURN_FAILURE("PIAA mirror surface 1 not identified");
        }
        return RETURN_SUCCESS;
    }

    piaacmcopticalsystem.elemtype[elem] = 4;
    DEBUG_TRACEPOINT("Element %ld is lens", elem);
    imageID IDr = image_ID(rname, dcimg, dcnimg);
    piaacmcopticalsystem.ASPHSURFRarray[surfslot].surfID = IDr;
    if (IDr == -1)
    {
        if (piaa_idx == 0)
        {
            FUNC_RETURN_FAILURE("PIAA lens surface 0 not identified");
        }
        FUNC_RETURN_FAILURE("PIAA lens surface 1 not identified");
    }

    piaacmcopticalsystem.ASPHSURFRarray[surfslot].mat0 = 100;
    piaacmcopticalsystem.ASPHSURFRarray[surfslot].mat1 = piaacmcopticaldesign.PIAAmaterial_code;

    return RETURN_SUCCESS;
}

/**
 * piaacmc_optsyst_add_piaa_stage - Add upstream/downstream masks and PIAA optics 0 and 1.
 * @elem:       Current element index cursor
 * @size:       2D wavefront array dimension [pix]
 * @beamradpix: Nominal beam radius [pix]
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
errno_t piaacmc_optsyst_add_piaa_stage(
    long  *elem,
    long   size,
    double beamradpix)
{
    FUNC_CHECK_RETURN(add_pre_piaa0_mask(elem));

    if (piaacmcopticaldesign.PIAAmode == 1)
    {
        long idx = *elem;

        snprintf(piaacmcopticalsystem.name[idx], STRINGMAXLEN_OPTSYST_ELEMNAME, "PIAA optics 0");
        piaacmcopticalsystem.elemtype[idx]       = 3;
        piaacmcopticalsystem.elemarrayindex[idx] = 1;
        piaacmcopticalsystem.elemZpos[idx]       = piaacmcopticaldesign.PIAA0pos;
        printf("============ (2) PIAA0pos = %f ==================\n",
               piaacmcopticalsystem.elemZpos[idx]);

        FUNC_CHECK_RETURN(bind_piaa_surface(idx, 1, "piaam0z", "piaar0zsag", 0));
        (*elem)++;
    }

    add_post_piaa0_mask(elem);

    if (piaacmcopticaldesign.PIAAmode == 1)
    {
        long idx = *elem;

        snprintf(piaacmcopticalsystem.name[idx], STRINGMAXLEN_OPTSYST_ELEMNAME, "PIAA optics 1");
        piaacmcopticalsystem.elemtype[idx]       = 3;
        piaacmcopticalsystem.elemarrayindex[idx] = 2;
        piaacmcopticalsystem.elemZpos[idx]       =
            piaacmcopticaldesign.PIAA0pos + piaacmcopticaldesign.PIAAsep;

        FUNC_CHECK_RETURN(bind_piaa_surface(idx, 2, "piaam1z", "piaar1zsag", 1));
        (*elem)++;

        idx = *elem;
        snprintf(
            piaacmcopticalsystem.name[idx],
            STRINGMAXLEN_OPTSYST_ELEMNAME,
            "opaque mask at PIAA elem 1");
        piaacmcopticalsystem.elemtype[idx] = 1;

        imageID ID = -1;
        load_fits("piaa1mask.fits", "piaa1mask", LOADFITS_ERRMODE_IGNORE, &ID);
        if (ID == -1)
        {
            ID = make_disk(
                "piaa1mask", size, size, 0.5 * size, 0.5 * size,
                piaacmcopticaldesign.r1lim * beamradpix);
        }
        piaacmcopticalsystem.elemarrayindex[idx] = ID;
        piaacmcopticalsystem.elemZpos[idx]       = piaacmcopticalsystem.elemZpos[idx - 1];
        (*elem)++;
    }

    return RETURN_SUCCESS;
}

/**
 * piaacmc_optsyst_add_fpm - Synthesize complex focal plane mask and register OPTSYST FPM element.
 * @elem: Current element index cursor
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
errno_t piaacmc_optsyst_add_fpm(long *elem)
{
    if (variable_ID("PIAACMC_NOFPM") != -1)
    {
        return RETURN_SUCCESS;
    }

    long idx = *elem;
    snprintf(
        piaacmcopticalsystem.name[idx],
        STRINGMAXLEN_OPTSYST_ELEMNAME,
        "post focal plane mask pupil");
    piaacmcopticalsystem.elemtype[idx]       = 5;
    piaacmcopticalsystem.elemarrayindex[idx] = 0;

    printf("=========== MAKE FOCAL PLANE MASK ===========\n");
    int        savefpm = 0;
    variableID IDv     = variable_ID("PIAACMC_SAVE_fpm");
    if (IDv != -1)
    {
        savefpm = (int) (coffee_dcvar_f(IDv) + 0.001);
    }

    FUNC_CHECK_RETURN(
        mkFocalPlaneMask(
            "fpmzmap",
            "piaacmcfpm",
            piaacmcparams.focmMode,
            savefpm,
            &(piaacmcopticalsystem.FOCMASKarray[0].fpmID)));

    mk_reim_from_complex("piaacmcfpm", "piaacmcfpm_re", "piaacmcfpm_im", 0);
    FUNC_CHECK_RETURN(save_fits("piaacmcfpm_re", "./testdir/test_piaacmcfpm_re.fits"));
    FUNC_CHECK_RETURN(save_fits("piaacmcfpm_im", "./testdir/test_piaacmcfpm_im.fits"));
    FUNC_CHECK_RETURN(delete_image_ID("piaacmcfpm_re", DELETE_IMAGE_ERRMODE_WARNING));
    FUNC_CHECK_RETURN(delete_image_ID("piaacmcfpm_im", DELETE_IMAGE_ERRMODE_WARNING));

    piaacmcopticalsystem.FOCMASKarray[0].zfactor = piaacmcopticaldesign.fpzfactor;
    piaacmcopticalsystem.elemZpos[idx]           = piaacmcopticalsystem.elemZpos[idx - 1];
    printf("=========== FOCAL PLANE MASK : DONE ===========\n");
    (*elem)++;

    return RETURN_SUCCESS;
}

/**
 * append_inv_piaa_pair - Append inverse PIAA optics (element 1 then element 0) to OPTSYST.
 * @elem: Current element index cursor
 */
static void append_inv_piaa_pair(long *elem)
{
    int  etype = (piaacmcopticaldesign.PIAAmaterial_code == 0) ? 3 : 4;
    long idx   = *elem;

    snprintf(piaacmcopticalsystem.name[idx], STRINGMAXLEN_OPTSYST_ELEMNAME, "invPIAA optics 1");
    piaacmcopticalsystem.elemtype[idx]       = etype;
    piaacmcopticalsystem.elemarrayindex[idx] = 2;
    piaacmcopticalsystem.elemZpos[idx]       = 0.0;
    (*elem)++;

    idx = *elem;
    snprintf(piaacmcopticalsystem.name[idx], STRINGMAXLEN_OPTSYST_ELEMNAME, "invPIAA optics 0");
    piaacmcopticalsystem.elemtype[idx]       = etype;
    piaacmcopticalsystem.elemarrayindex[idx] = 1;
    piaacmcopticalsystem.elemZpos[idx]       = piaacmcopticaldesign.PIAAsep;
    (*elem)++;
}

/**
 * piaacmc_optsyst_add_lyot_and_backend - Append inverse PIAA, Lyot stops, and back-end stop.
 * @elem: Current element index cursor
 * @size: 2D wavefront array dimension [pix]
 */
void piaacmc_optsyst_add_lyot_and_backend(
    long *elem,
    long  size)
{
    if ((piaacmcopticaldesign.PIAAmode == 1) && (piaacmcopticaldesign.invPIAAmode == 2))
    {
        append_inv_piaa_pair(elem);
    }

    for (long i = 0; i < piaacmcopticaldesign.NBLyotStop; i++)
    {
        long idx = *elem;

        snprintf(piaacmcopticalsystem.name[idx], STRINGMAXLEN_OPTSYST_ELEMNAME, "Lyot mask %ld", i);
        piaacmcopticalsystem.elemtype[idx]       = 1;
        piaacmcopticalsystem.elemarrayindex[idx] = piaacmcopticaldesign.IDLyotStop[i];
        printf("elem %ld  Lyot mask %ld : %ld\n", idx, i, piaacmcopticaldesign.IDLyotStop[i]);
        piaacmcopticalsystem.elemZpos[idx] = piaacmcopticaldesign.LyotStop_zpos[i];
        (*elem)++;
    }

    if ((piaacmcopticaldesign.PIAAmode == 1) && (piaacmcopticaldesign.invPIAAmode == 1))
    {
        append_inv_piaa_pair(elem);
    }

    if (piaacmcopticaldesign.PIAAmode == 1)
    {
        long idx = *elem;

        snprintf(
            piaacmcopticalsystem.name[idx],
            STRINGMAXLEN_OPTSYST_ELEMNAME,
            "back end pupil stop  (rad = %f)",
            piaacmcopticaldesign.pupoutmaskrad);
        piaacmcopticalsystem.elemtype[idx] = 1;
        piaacmcopticalsystem.elemarrayindex[idx] = make_disk(
            "pupoutmask",
            size,
            size,
            0.5 * size,
            0.5 * size,
            piaacmcopticaldesign.pupoutmaskrad * piaacmcopticaldesign.beamrad /
                piaacmcopticaldesign.pixscale);
        piaacmcopticalsystem.elemZpos[idx] = piaacmcopticalsystem.elemZpos[idx - 1];
        (*elem)++;
    }
}
