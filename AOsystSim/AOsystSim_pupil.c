/**
 * @file    AOsystSim_pupil.c
 * @brief   Parametric telescope pupil synthesis on DM grid and illumination geometry fitting
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "milk_runnable.h"
#include "coffee_compat.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "info/info.h"

#include "AOsystSim/AOsystSim.h"
#include "AOsystSim_internal.h"

/**
 * @brief Candidate or best-fit pupil geometry parameter set.
 */
struct pupfit_state
{
    double xc;
    double yc;
    double pupPA;
    double rout;
    double stretchx;
    double spiderPA;
    double spideroffset;
    double rms;
};

/**
 * @brief Search interval bounds for coarse-to-fine pupil geometry grid search.
 */
struct pupfit_bounds
{
    double xc_min;
    double xc_max;
    double yc_min;
    double yc_max;
    double pupPA_min;
    double pupPA_max;
    double rout_min;
    double rout_max;
    double stretchx_min;
    double stretchx_max;
    double spiderPA_min;
    double spiderPA_max;
    double spideroffset_min;
    double spideroffset_max;
};

/**
 * pup_classify_spider - Classify subpixel coordinates into one of four spider-bounded sectors
 * @x:            Stretched/rotated X coordinate in DM actuator units
 * @y:            Rotated Y coordinate in DM actuator units
 * @spiderPA:     Spider opening half-angle [rad]
 * @spideroffset: Spider lateral offset [actuators]
 * @spiderthick:  Spider vane width [actuators]
 *
 * Return: Sector index 1..4 inside illuminated quadrants, or 0 if obscured by spiders.
 */
static long pup_classify_spider(
    double x,
    double y,
    double spiderPA,
    double spideroffset,
    double spiderthick)
{
    long   index = 0;
    double x1;
    double PA;

    x1 = x - spideroffset - spiderthick / 2.0;
    PA = atan2(y, x1);
    if (fabs(PA) < spiderPA)
    {
        index = 1;
    }

    x1 = x + spideroffset + spiderthick / 2.0;
    PA = atan2(y, -x1);
    if (fabs(PA) < spiderPA)
    {
        index = 2;
    }

    x1 = x + spideroffset - spiderthick / 2.0;
    PA = atan2(x1, y);
    if ((fabs(PA) < M_PI / 2.0 - spiderPA) && (x < 0.0))
    {
        index = 3;
    }

    x1 = -x + spideroffset - spiderthick / 2.0;
    PA = atan2(x1, y);
    if ((fabs(PA) < M_PI / 2.0 - spiderPA) && (x > 0.0))
    {
        index = 3;
    }

    x1 = x + spideroffset - spiderthick / 2.0;
    PA = atan2(x1, -y);
    if ((fabs(PA) < M_PI / 2.0 - spiderPA) && (x < 0.0))
    {
        index = 4;
    }

    x1 = -x + spideroffset - spiderthick / 2.0;
    PA = atan2(x1, -y);
    if ((fabs(PA) < M_PI / 2.0 - spiderPA) && (x > 0.0))
    {
        index = 4;
    }

    return index;
}

/**
 * AOsystSim_mkTelPupDM - Synthesize parametric obstructed/spidered telescope pupil on DM grid
 * @ID_name:      Output 2D pupil image name
 * @msize:        DM grid linear dimension [actuators]
 * @xc:           Pupil center X coordinate [actuators]
 * @yc:           Pupil center Y coordinate [actuators]
 * @rin:          Central obstruction radius [actuators]
 * @rout:         Outer pupil radius [actuators]
 * @pupPA:        Pupil position angle [rad]
 * @spiderPA:     Spider vane angle [rad]
 * @spideroffset: Spider vane lateral offset [actuators]
 * @spiderthick:  Spider vane thickness [actuators]
 * @stretchx:     Anamorphic X stretch factor
 *
 * Return: Image ID of synthesized pupil map.
 */
long AOsystSim_mkTelPupDM(
    const char *ID_name,
    long        msize,
    double      xc,
    double      yc,
    double      rin,
    double      rout,
    double      pupPA,
    double      spiderPA,
    double      spideroffset,
    double      spiderthick,
    double      stretchx)
{
    long    binfact = 8;
    long    size    = msize * binfact;
    imageID ID;
    imageID IDz;
    imageID IDindex;
    imageID IDi;

    aosim_ensure_wdir();

    create_2Dimage_ID(ID_name, msize, msize, &ID);
    create_3Dimage_ID("TPind", msize, msize, 5, &IDi);
    create_2Dimage_ID("telpupDMz", size, size, &IDz);
    create_2Dimage_ID("telpupDMzindex", size, size, &IDindex);

    for (long ii = 0; ii < size; ii++)
    {
        for (long jj = 0; jj < size; jj++)
        {
            double val = 1.0;
            double x   = (1.0 * ii / binfact) - xc;
            double y   = (1.0 * jj / binfact) - yc;
            double r   = sqrt(x * x + y * y);
            double PA  = atan2(y, x) + pupPA;
            long   index;
            long   ii1 = ii / binfact;
            long   jj1 = jj / binfact;

            x = r * cos(PA) * stretchx;
            y = r * sin(PA);
            if ((r > rout) || (r < rin))
            {
                val = 0.0;
            }

            index = pup_classify_spider(x, y, spiderPA, spideroffset, spiderthick);
            if (index == 0)
            {
                val = 0.0;
            }

            dcimg[IDz].array.F[jj * size + ii]     = (float) val;
            dcimg[IDindex].array.F[jj * size + ii] = (float) (index * val);
            dcimg[ID].array.F[jj1 * msize + ii1] += (float) (val / binfact / binfact);

            if (val > 0.5)
            {
                dcimg[IDi].array.F[jj1 * msize + ii1]                       = 1.0f;
                dcimg[IDi].array.F[index * msize * msize + jj1 * msize + ii1] = 1.0f;
            }
        }
    }

    delete_image_ID("telpupDMz", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("telpupDMzindex", DELETE_IMAGE_ERRMODE_WARNING);
    save_fits("TPind", "AOsystSim_wdir/TPind.fits");

    return ID;
}

/**
 * pupfit_norm_percentiles - Rescale image intensities linearly between 10th and 90th percentiles
 * @ID_name: Image name for percentile query
 * @ID:      Image identifier
 * @verbose: Non-zero to print percentile values
 */
static void pupfit_norm_percentiles(
    const char *ID_name,
    imageID     ID,
    int         verbose)
{
    double   vp10 = img_percentile_float(ID_name, 0.1);
    double   vp90 = img_percentile_float(ID_name, 0.9);
    uint64_t npix = (uint64_t) dcimg[ID].md[0].size[0] * dcimg[ID].md[0].size[1];

    if (verbose)
    {
        printf("%f %f\n", vp10, vp90);
    }
    for (uint64_t ii = 0; ii < npix; ii++)
    {
        dcimg[ID].array.F[ii] = (float) ((dcimg[ID].array.F[ii] - vp10) / (vp90 - vp10));
    }
}

/**
 * pupfit_compensate_gradient - Estimate and flatten 2D linear illumination gradient across DM map
 * @ID_name: Measured DM response image name
 * @ID:      Measured DM response image ID
 * @size:    Grid linear dimension
 */
static void pupfit_compensate_gradient(
    const char *ID_name,
    imageID     ID,
    long        size)
{
    imageID ID1;
    double  val1    = 1.0e15;
    double  coeffx1 = 0.0;
    double  coeffy1 = 0.0;

    pupfit_norm_percentiles(ID_name, ID, 1);

    create_2Dimage_ID("tmpftpim", size, size, &ID1);
    for (double coeffx = -1.5; coeffx < 1.5; coeffx += 0.05)
    {
        for (double coeffy = -1.5; coeffy < 1.5; coeffy += 0.05)
        {
            double val;

            for (long ii = 0; ii < size; ii++)
            {
                for (long jj = 0; jj < size; jj++)
                {
                    double x = 1.0 * ii / size;
                    double y = 1.0 * jj / size;

                    dcimg[ID1].array.F[jj * size + ii] =
                        (float) (dcimg[ID].array.F[jj * size + ii] *
                                 (1.0 + coeffx * (x - 0.5)) * (1.0 + coeffy * (y - 0.5)));
                }
            }
            val = img_percentile_float("tmpftpim", 0.9) - img_percentile_float("tmpftpim", 0.7);
            if (val < val1)
            {
                val1    = val;
                coeffx1 = coeffx;
                coeffy1 = coeffy;
            }
            printf("     %f %f   %g       ( %f %f %g )\n",
                   coeffx, coeffy, val, coeffx1, coeffy1, val1);
        }
    }

    printf("COEFF : %f %f\n", coeffx1, coeffy1);
    for (long ii = 0; ii < size; ii++)
    {
        for (long jj = 0; jj < size; jj++)
        {
            double x = 1.0 * ii / size;
            double y = 1.0 * jj / size;

            dcimg[ID].array.F[jj * size + ii] =
                (float) (dcimg[ID].array.F[jj * size + ii] *
                         (1.0 + coeffx1 * (x - 0.5)) * (1.0 + coeffy1 * (y - 0.5)));
        }
    }

    pupfit_norm_percentiles(ID_name, ID, 1);

    for (long ii = 0; ii < size * size; ii++)
    {
        if (dcimg[ID].array.F[ii] > 0.05f)
        {
            dcimg[ID].array.F[ii] = sqrtf(dcimg[ID].array.F[ii]);
        }
        else
        {
            dcimg[ID].array.F[ii] = 0.0f;
        }
    }
    pupfit_norm_percentiles(ID_name, ID, 0);
}

/**
 * pupfit_eval_pose_grid - Sweep (xc, yc, pupPA) inner grid for fixed outer pupil/spider parameters
 * @ID:   Normalized measured DM response image ID
 * @size: Grid linear dimension
 * @b:    Current search bounds
 * @cand: Candidate state with outer parameters populated
 * @best: Best-fit state updated in place when RMS improves
 */
static void pupfit_eval_pose_grid(
    imageID                     ID,
    long                        size,
    const struct pupfit_bounds *b,
    const struct pupfit_state  *cand,
    struct pupfit_state        *best)
{
    long   nbstep = 2;
    double eps    = 1.0e-8;
    double dx     = (b->xc_max - b->xc_min) / nbstep;
    double dy     = (b->yc_max - b->yc_min) / nbstep;
    double dpa    = (b->pupPA_max - b->pupPA_min) / nbstep;

    for (double xc = b->xc_min; xc < b->xc_max + eps; xc += dx)
    {
        for (double yc = b->yc_min; yc < b->yc_max + eps; yc += dy)
        {
            for (double pupPA = b->pupPA_min; pupPA < b->pupPA_max + eps; pupPA += dpa)
            {
                double  rin          = 0.305 * cand->rout;
                double  spiderthick  = 0.06 * cand->rout;
                double  spideroffset = 0.165 * cand->rout;
                double  rms          = 0.0;
                imageID IDt;
                FILE   *fp;

                AOsystSim_mkTelPupDM("testpup", size, xc, yc, rin, cand->rout, pupPA,
                                     cand->spiderPA, spideroffset, spiderthick, cand->stretchx);
                IDt = image_ID("testpup", dcimg, dcnimg);

                for (long ii = 0; ii < size * size; ii++)
                {
                    double v1 = dcimg[ID].array.F[ii] - dcimg[IDt].array.F[ii];

                    rms += v1 * v1;
                }
                rms = sqrt(rms / size / size);

                if (rms < best->rms)
                {
                    best->rms          = rms;
                    best->rout         = cand->rout;
                    best->xc           = xc;
                    best->yc           = yc;
                    best->pupPA        = pupPA;
                    best->stretchx     = cand->stretchx;
                    best->spiderPA     = cand->spiderPA;
                    best->spideroffset = spideroffset;
                }
                delete_image_ID("testpup", DELETE_IMAGE_ERRMODE_WARNING);
                printf("%f %f %f  %f -> rms = %g\n", cand->rout, xc, yc, pupPA, rms);
                fp = fopen("pupfit.txt", "a");
                fprintf(fp, "%f %f %f %f %g\n", cand->rout, xc, yc, pupPA, rms);
                fclose(fp);
            }
        }
    }
}

/**
 * pupfit_grid_iter - Execute one coarse-to-fine 7D grid search pass and contract bounds by 1/3
 * @ID:   Normalized measured DM response image ID
 * @size: Grid linear dimension
 * @iter: Iteration index
 * @b:    Search bounds contracted in place around best solution
 * @best: Best-fit state updated in place
 */
static void pupfit_grid_iter(
    imageID               ID,
    long                  size,
    long                  iter,
    struct pupfit_bounds *b,
    struct pupfit_state  *best)
{
    long   nbstep = 2;
    double eps    = 1.0e-8;
    double d_spa  = (b->spiderPA_max - b->spiderPA_min) / nbstep;
    double d_soff = (b->spideroffset_max - b->spideroffset_min) / nbstep;
    double d_stx  = (b->stretchx_max - b->stretchx_min) / nbstep;
    double d_rout = (b->rout_max - b->rout_min) / nbstep;
    double range;
    double xc_range;

    for (double spa = b->spiderPA_min; spa < b->spiderPA_max + eps; spa += d_spa)
    {
        for (double soff = b->spideroffset_min; soff < b->spideroffset_max + eps; soff += d_soff)
        {
            for (double stx = b->stretchx_min; stx < b->stretchx_max + eps; stx += d_stx)
            {
                for (double rout = b->rout_min; rout < b->rout_max + eps; rout += d_rout)
                {
                    struct pupfit_state cand = { 0 };

                    cand.spiderPA     = spa;
                    cand.spideroffset = soff;
                    cand.stretchx     = stx;
                    cand.rout         = rout;
                    pupfit_eval_pose_grid(ID, size, b, &cand, best);
                }
            }
        }
    }

    printf("ITERATION %ld    %f %f  %f  %f  ->  %f\n",
           iter, best->rout, best->xc, best->yc, best->pupPA, best->rms);

    range               = (b->spideroffset_max - b->spideroffset_min) / 3.0;
    b->spideroffset_min = best->spideroffset - 0.5 * range;
    b->spideroffset_max = best->spideroffset + 0.5 * range;

    range           = (b->spiderPA_max - b->spiderPA_min) / 3.0;
    b->spiderPA_min = best->spiderPA - 0.5 * range;
    b->spiderPA_max = best->spiderPA + 0.5 * range;

    range           = (b->stretchx_max - b->stretchx_min) / 3.0;
    b->stretchx_min = best->stretchx - 0.5 * range;
    b->stretchx_max = best->stretchx + 0.5 * range;

    range       = (b->rout_max - b->rout_min) / 3.0;
    b->rout_min = best->rout - 0.5 * range;
    b->rout_max = best->rout + 0.5 * range;

    xc_range  = (b->xc_max - b->xc_min) / 3.0;
    b->xc_min = best->xc - 0.5 * xc_range;
    b->xc_max = best->xc + 0.5 * xc_range;

    range     = (b->yc_max - b->yc_min) / 3.0;
    b->yc_min = best->yc - 0.5 * range;
    b->yc_max = best->yc + 0.5 * range;

    range        = (b->xc_max - b->xc_min) / 3.0;
    b->pupPA_min = best->xc - 0.5 * xc_range;
    b->pupPA_max = best->xc + 0.5 * range;
}

/**
 * AOsystSim_fitTelPup - Fit telescope pupil registration and spider geometry to measured DM map
 * @ID_name:       Input measured DM response image name
 * @IDtelpup_name: Output fitted telescope pupil image name
 *
 * Return: Image ID of input DM response map.
 */
long AOsystSim_fitTelPup(
    const char *ID_name,
    const char *IDtelpup_name)
{
    struct pupfit_state  best   = { 0 };
    struct pupfit_bounds bounds = { 0 };
    imageID              ID     = image_ID(ID_name, dcimg, dcnimg);
    long                 size   = dcimg[ID].md[0].size[0];
    long                 NBiter = 5;
    FILE                *fp;

    best.rout         = 22.15;
    best.pupPA        = -0.105;
    best.spiderPA     = 0.85;
    best.spideroffset = 0.15 * best.rout;
    best.xc           = 24.629630;
    best.yc           = 23.518519;
    best.stretchx     = 1.048148;
    best.rms          = 1.0e21;

    pupfit_compensate_gradient(ID_name, ID, size);

    bounds.rout_min         = best.rout * 0.9;
    bounds.rout_max         = best.rout * 1.1;
    bounds.xc_min           = best.xc - 2.0;
    bounds.xc_max           = best.xc + 2.0;
    bounds.yc_min           = best.yc - 2.0;
    bounds.yc_max           = best.yc + 2.0;
    bounds.pupPA_min        = best.pupPA - 0.1;
    bounds.pupPA_max        = best.pupPA + 0.1;
    bounds.stretchx_min     = best.stretchx * 0.95;
    bounds.stretchx_max     = best.stretchx * 1.05;
    bounds.spideroffset_min = best.spideroffset * 0.9;
    bounds.spideroffset_max = best.spideroffset * 1.1;
    bounds.spiderPA_min     = best.spiderPA - 0.1;
    bounds.spiderPA_max     = best.spiderPA + 0.1;

    fp = fopen("pupfit.txt", "w");
    fclose(fp);

    for (long iter = 0; iter < NBiter; iter++)
    {
        pupfit_grid_iter(ID, size, iter, &bounds, &best);
    }

    printf("BEST SOLUTION: \n");
    printf("     xc            %f\n", best.xc);
    printf("     yc            %f\n", best.yc);
    printf("     rout          %f\n", best.rout);
    printf("     pupPA         %f\n", best.pupPA);
    printf("     spiderPA      %f\n", best.spiderPA);
    printf("     spideroffset  %f\n", best.spideroffset);
    printf("     stretchx      %f\n", best.stretchx);

    AOsystSim_mkTelPupDM(IDtelpup_name, size, best.xc, best.yc, 0.315 * best.rout, best.rout,
                         best.pupPA, best.spiderPA, best.spideroffset, 0.06 * best.rout,
                         best.stretchx);

    return ID;
}
