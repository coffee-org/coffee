/**
 * @file    exec_compute_image.c
 * @brief   On-axis PSF, shared-memory loop, or multi-point scene synthesis (Mode 0)
 */

#include <semaphore.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"
#include "PIAACMCsimul_computePSF.h"
#include "exec_compute_image.h"
#include "init_piaacmcopticaldesign.h"

#include "PIAAshape/makePIAAshapes.h"

/**
 * compute_image_init_design - Read Mode 0 CLI variables and initialize optical design and PIAA.
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t compute_image_init_design(void)
{
    double     fpmradld        = 0.95;
    double     centobs0        = 0.3;
    double     centobs1        = 0.2;
    int        PIAACMC_WFCmode = 0;
    variableID IDv;

    if ((IDv = variable_ID("PIAACMC_centobs0")) != -1)
    {
        centobs0 = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_centobs1")) != -1)
    {
        centobs1 = coffee_dcvar_f(IDv);
    }
    if ((IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = coffee_dcvar_f(IDv);
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }

    piaacmcparams.PIAACMC_fpmtype = 0;
    if ((IDv = variable_ID("PIAACMC_fpmtype")) != -1)
    {
        piaacmcparams.PIAACMC_fpmtype = (int) (coffee_dcvar_f(IDv) + 0.1);
    }
    printf("PIAACMC_fpmtype = %d\n", piaacmcparams.PIAACMC_fpmtype);

    if ((IDv = variable_ID("PIAACMC_WFCmode")) != -1)
    {
        PIAACMC_WFCmode = (int) (coffee_dcvar_f(IDv) + 0.1);
    }
    printf("PIAACMC_WFCmode = %d\n", PIAACMC_WFCmode);

    piaacmcparams.FORCE_CREATE_fpmza = 1;
    {
        uint64_t initflag = INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__READCONF;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF;

        if (piaacmcparams.PIAACMC_fpmtype == 1)
        {
            initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL;
        }
        if (PIAACMC_WFCmode == 1)
        {
            initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__WSCMODE;
        }

        FUNC_CHECK_RETURN(
            init_piaacmcopticaldesign(fpmradld, centobs0, centobs1, initflag, NULL));
    }

    FUNC_CHECK_RETURN(makePIAAshapes());
    piaacmcopticalsystem.FOCMASKarray[0].mode = 1;

    return RETURN_SUCCESS;
}

/**
 * copy_cropped_psf_cube - Copy central cropped region of full PSF cube into output stream.
 * @IDsrc:    Full-frame 3D PSF image ID
 * @IDdst:    Cropped 3D shared-memory PSF image ID
 * @sizecrop: Spatial crop dimension [pix]
 */
static void copy_cropped_psf_cube(
    imageID IDsrc,
    imageID IDdst,
    long    sizecrop)
{
    long fullsize = piaacmcopticaldesign.size;
    long offset   = (fullsize - sizecrop) / 2;

    for (long k = 0; k < piaacmcopticaldesign.nblambda; k++)
    {
        for (long ii1 = 0; ii1 < sizecrop; ii1++)
        {
            for (long jj1 = 0; jj1 < sizecrop; jj1++)
            {
                long ii  = ii1 + offset;
                long jj  = jj1 + offset;
                long dst = k * sizecrop * sizecrop + jj1 * sizecrop + ii1;
                long src = k * fullsize * fullsize + jj * fullsize + ii;

                dcimg[IDdst].array.F[dst] = dcimg[IDsrc].array.F[src];
            }
        }
    } // for k
}

/**
 * compute_image_loopmode - Stream cropped PSF cubes synchronized to shared-memory opderr updates.
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t compute_image_loopmode(void)
{
    printf("RUNNING PSF LOOP COMPUTATION\n");

    uint32_t *sizearray = (uint32_t *) malloc(sizeof(uint32_t) * 2);
    if (sizearray == NULL)
    {
        FUNC_RETURN_FAILURE("malloc error");
    }
    sizearray[0] = piaacmcopticaldesign.size;
    sizearray[1] = piaacmcopticaldesign.size;

    imageID IDopderrC;
    errno_t err = create_image_ID("opderr", 2, sizearray, _DATATYPE_FLOAT, 1, 0, 0, &IDopderrC);
    free(sizearray);
    FUNC_CHECK_RETURN(err);
    COREMOD_MEMORY_image_set_createsem("opderr", 10);

    long sizecrop = piaacmcopticaldesign.size / 16;
    sizearray     = (uint32_t *) malloc(sizeof(uint32_t) * 3);
    if (sizearray == NULL)
    {
        FUNC_RETURN_FAILURE("malloc error");
    }
    sizearray[0] = sizecrop;
    sizearray[1] = sizecrop;
    sizearray[2] = piaacmcopticaldesign.nblambda;

    imageID IDpsfi0;
    err = create_image_ID("psfiout0", 3, sizearray, _DATATYPE_FLOAT, 1, 0, 0, &IDpsfi0);
    free(sizearray);
    FUNC_CHECK_RETURN(err);

    long iter = 0;
    while (iter < 10)
    {
        double cval = 0.0;
        FUNC_CHECK_RETURN(
            PIAACMCsimul_computePSF(0.0, 0.0, 0, piaacmcopticalsystem.NBelem, 1, 0, 0, 1, &cval));

        imageID ID = image_ID("psfi0", dcimg, dcnimg);
        dcimg[IDpsfi0].md[0].write = 1;
        copy_cropped_psf_cube(ID, IDpsfi0, sizecrop);
        COREMOD_MEMORY_image_set_sempost_byID(IDpsfi0, -1);
        dcimg[IDpsfi0].md[0].cnt0++;
        dcimg[IDpsfi0].md[0].write = 0;

        COREMOD_MEMORY_image_set_semwait("opderr", 0);
        while (sem_trywait(dcimg[IDopderrC].semptr[0]) == 0)
        {
        }
    } // while (iter < 10)

    return RETURN_SUCCESS;
}

/**
 * compute_image_scene - Sum weighted off-axis PSFs listed in SCENE.txt and write scene.fits.
 * @fpscene: Open file stream positioned at start of SCENE.txt
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
static errno_t compute_image_scene(FILE *fpscene)
{
    int     initscene = 0;
    imageID IDscene   = -1;
    double  xpos;
    double  ypos;
    double  fval;

    while (fscanf(fpscene, "%lf %lf %lf\n", &xpos, &ypos, &fval) == 3)
    {
        printf("COMPUTING PSF AT POSITION %lf %lf, flux  = %g\n", xpos, ypos, fval);
        FUNC_CHECK_RETURN(
            PIAACMCsimul_computePSF(xpos, ypos, 0, piaacmcopticalsystem.NBelem, 1, 0, 0, 1, NULL));

        imageID  ID    = image_ID("psfi0", dcimg, dcnimg);
        uint32_t xsize = dcimg[ID].md[0].size[0];
        uint32_t ysize = dcimg[ID].md[0].size[1];
        uint32_t zsize = dcimg[ID].md[0].size[2];

        if (initscene == 0)
        {
            initscene = 1;
            FUNC_CHECK_RETURN(create_3Dimage_ID("scene", xsize, ysize, zsize, &IDscene));
        }

        ID = image_ID("psfi0", dcimg, dcnimg);
        uint64_t nelem = (uint64_t) xsize * ysize * zsize;
        for (uint64_t ii = 0; ii < nelem; ii++)
        {
            dcimg[IDscene].array.F[ii] += fval * dcimg[ID].array.F[ii];
        }
    } // while fscanf

    save_fits("scene", "scene.fits");
    return RETURN_SUCCESS;
}

/**
 * exec_compute_image - Compute on-axis PSF, closed-loop stream, or astronomical scene (Mode 0).
 *
 * Return: RETURN_SUCCESS on success, error code on failure.
 */
errno_t exec_compute_image(void)
{
    DEBUG_TRACE_FSTART();

    printf("=================================== mode 000 ===================================\n");

    FUNC_CHECK_RETURN(compute_image_init_design());

    FILE *fp = fopen("LOOPMODE.txt", "r");
    if (fp != NULL)
    {
        fclose(fp);
        FUNC_CHECK_RETURN(compute_image_loopmode());
        DEBUG_TRACE_FEXIT();
        return RETURN_SUCCESS;
    }

    FILE *fpscene = fopen("SCENE.txt", "r");
    if (fpscene != NULL)
    {
        errno_t err = compute_image_scene(fpscene);
        fclose(fpscene);
        FUNC_CHECK_RETURN(err);
        DEBUG_TRACE_FEXIT();
        return RETURN_SUCCESS;
    }

    double  valref = 0.0;
    errno_t fret   = PIAACMCsimul_computePSF(
        0.0, 0.0, 0, piaacmcopticalsystem.NBelem, 1, 0, 0, 1, &valref);
    if (fret != RETURN_SUCCESS)
    {
        FUNC_RETURN_FAILURE("Call to PIAACMCsimul_computePSF failed");
    }
    printf("valref = %g\n", valref);

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
