/**
 * @file    PIAAACMCsimul_mkFPM_zonemap.c
 * @brief   PIAA-type coronagraph design, initialize
 *
 */

#include <math.h>
#include <stdio.h>

// milk includes

#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "image_gen/image_gen.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"

/**
 * @param[out]  IDname  Name of output image
 */
errno_t mkFPM_zonemap(const char *__restrict__ IDname, imageID *outID)
{
    DEBUG_TRACE_FSTART();

    FILE     *fp;
    char      fname[STRINGMAXLEN_FULLFILENAME];
    uint32_t *sizearray;

    uint_fast32_t *nbsector;
    uint_fast32_t *nbsectorcumul;
    double         eps = 1.0e-6;

    uint_fast64_t cnt, cnt1;
    uint_fast32_t nbzonescc;

    double        hexstep = 1.0;
    double        hexgap  = -0.0001;
    double        hexsteppix;
    double        hex_x[10000];
    double        hex_y[10000];
    long          hex_ring[10000];
    uint_fast32_t hex_number[10000];
    uint_fast32_t hcnt;
    uint_fast32_t hindex;
    uint_fast32_t hindexMax;
    imageID       ID1;

    printf("function PIAACMCsimul_mkFPM_zonemap\n");
    fflush(stdout);

    sizearray = (uint32_t *) malloc(sizeof(uint32_t) * 2);
    if(sizearray == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }
    sizearray[0] = piaacmcopticaldesign.fpmarraysize;
    sizearray[1] = piaacmcopticaldesign.fpmarraysize;
    imageID ID;
    create_image_ID(IDname, 2, sizearray, _DATATYPE_UINT16, 0, 0, 0, &ID);
    free(sizearray);

    nbsector = (uint_fast32_t *) malloc(sizeof(uint_fast32_t) *
                                        piaacmcopticaldesign.NBrings);
    if(nbsector == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }

    nbsectorcumul = (uint_fast32_t *) malloc(sizeof(uint_fast32_t) *
                    piaacmcopticaldesign.NBrings);
    if(nbsectorcumul == NULL)
    {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }

    switch(piaacmcparams.PIAACMC_FPMsectors)
    {

        case 0: // rings
            nbsectorcumul[0] = 1;
            for(long ring = 1; ring < piaacmcopticaldesign.NBrings; ring++)
            {
                nbsectorcumul[ring] = nbsectorcumul[ring - 1] + 1;
            }
            break;

        case 1: // rings broken in sectors
            WRITE_FULLFILENAME(fname,
                               "%s/fpm_zonescoord_%d_%03ld.txt",
                               piaacmcparams.piaacmcconfdir,
                               piaacmcparams.PIAACMC_FPMsectors,
                               piaacmcopticaldesign.NBrings);
            fp = fopen(fname, "w");
            //NBzones = 0;
            cnt = 0;
            fprintf(fp, "0 0\n");
            nbsector[0]      = 1;
            nbsectorcumul[0] = 1;
            for(long ring = 1; ring < piaacmcopticaldesign.NBrings; ring++)
            {

                nbsector[ring]      = 2 * (ring + 1);
                nbsectorcumul[ring] = nbsectorcumul[ring - 1] + nbsector[ring];
                for(cnt1 = 0; cnt1 < nbsector[ring]; cnt1++)
                {
                    cnt++;
                    fprintf(fp, "%ld %ld\n", cnt, ring);
                }
            }

            fclose(fp);
            for(long ring = 0; ring < piaacmcopticaldesign.NBrings; ring++)
            {
                printf("ring %ld : %ld %ld\n",
                       ring,
                       nbsector[ring],
                       nbsectorcumul[ring]);
            }
            break;

        case 2: // rings of hexagons
            WRITE_FULLFILENAME(fname,
                               "%s/fpm_zonescoord_%d_%03ld.txt",
                               piaacmcparams.piaacmcconfdir,
                               piaacmcparams.PIAACMC_FPMsectors,
                               piaacmcopticaldesign.NBrings);
            fp = fopen(fname, "w");
            fprintf(fp, "# focal plane mask zones geometry\n");
            fprintf(fp, "# hexagonal tiling\n");
            fprintf(fp, "# col 1: hexagon index\n");
            fprintf(fp, "# col 2: ring index\n");
            fprintf(fp, "# col 3: hexagon center x coordinate\n");
            fprintf(fp, "# col 4: hexagon center y coordinate\n");
            fprintf(fp,
                    "# Note: unit = hexagon center to tip radius (circumradius)\n");
            fprintf(fp, "# \n");
            nbsector[0]      = 1;
            nbsector[0]      = 1;
            nbsectorcumul[0] = 1;
            for(long ring = 1; ring < piaacmcopticaldesign.NBrings; ring++)
            {
                nbsector[ring]      = 0;
                nbsectorcumul[ring] = 0;
            }
            hindex = 0;
            // hegagon side = s = ring unit
            long ii1max = (long)(piaacmcopticaldesign.NBrings / 3 + 2);
            long jj1max = (long)(piaacmcopticaldesign.NBrings / sqrt(3.0) + 2);
            for(long ii1 = -ii1max; ii1 < ii1max; ii1++)
                for(long jj1 = -jj1max; jj1 < jj1max; jj1++)
                {
                    double hx   = hexstep * ii1 * 3;
                    double hy   = hexstep * sqrt(3.0) * jj1;
                    long   ring = (long) sqrt(hx * hx + hy * hy);
                    if(ring < piaacmcopticaldesign.NBrings)
                    {
                        nbsector[ring]++;
                        hex_x[hindex]    = hx;
                        hex_y[hindex]    = hy;
                        hex_ring[hindex] = ring;
                        hindex++;
                    }

                    hx += hexstep * 1.5;
                    hy += hexstep * sqrt(3.0) / 2.0;
                    ring = (long) sqrt(hx * hx + hy * hy);
                    if(ring < piaacmcopticaldesign.NBrings)
                    {
                        nbsector[ring]++;
                        hex_x[hindex]    = hx;
                        hex_y[hindex]    = hy;
                        hex_ring[hindex] = ring;
                        hindex++;
                    }
                }
            hindexMax = hindex;

            fprintf(fp, "%5ld %5ld  %11.6f %11.6f\n", (long) 0, (long) 0, 0.0, 0.0);
            hcnt = 1;
            for(long ring = 1; ring < piaacmcopticaldesign.NBrings; ring++)
            {
                for(hindex = 0; hindex < hindexMax; hindex++)
                    if(hex_ring[hindex] == ring)
                    {
                        hex_number[hindex] = hcnt;
                        fprintf(fp,
                                "%5ld %5ld  %11.6f %11.6f\n",
                                hcnt,
                                ring,
                                hex_x[hindex],
                                hex_y[hindex]);
                        hcnt++;
                    }
                if(ring > 0)
                {
                    nbsectorcumul[ring] = nbsectorcumul[ring - 1] + nbsector[ring];
                }
            }
            fclose(fp);
            for(long ring = 0; ring < piaacmcopticaldesign.NBrings; ring++)
            {
                printf("ring %ld : %ld %ld\n",
                       ring,
                       nbsector[ring],
                       nbsectorcumul[ring]);
            }
            break;

        default:
            printf("ERROR: FPMsector mode (%d) not recognized\n",
                   piaacmcparams.PIAACMC_FPMsectors);
            exit(0);
            break;
    }

    if(piaacmcopticaldesign.NBringCentCone > 0)
    {
        nbzonescc = nbsectorcumul[piaacmcopticaldesign.NBringCentCone - 1];
    }
    else
    {
        nbzonescc = 0;
    }

    //  ID = create_2Dimage_ID(IDname, piaacmcopticaldesign.fpmarraysize, piaacmcopticaldesign.fpmarraysize);

    if((piaacmcparams.PIAACMC_FPMsectors == 0) ||
            (piaacmcparams.PIAACMC_FPMsectors == 1))
    {
        for(long ii = 0; ii < piaacmcopticaldesign.fpmarraysize; ii++)
            for(long jj = 0; jj < piaacmcopticaldesign.fpmarraysize; jj++)
            {
                double x =
                    (2.0 * ii - 1.0 * piaacmcopticaldesign.fpmarraysize) /
                    piaacmcopticaldesign.fpmarraysize /
                    piaacmcparams.FPMSCALEFACTOR;

                double y =
                    (2.0 * jj - 1.0 * piaacmcopticaldesign.fpmarraysize) /
                    piaacmcopticaldesign.fpmarraysize /
                    piaacmcparams.FPMSCALEFACTOR;

                double r   = sqrt(x * x + y * y);
                double PA  = atan2(y, x);
                double PAf = 0.5 * ((PA / M_PI) + 1.0);
                if(PAf < eps)
                {
                    PAf = eps;
                }
                if(PAf > 1.0 - eps)
                {
                    PAf = 1.0 - eps;
                }

                long zi = (long) ceil((1.0 - r) * piaacmcopticaldesign.NBrings);

                if(zi < 0.1)
                {
                    zi = 0;
                }
                if(zi > piaacmcopticaldesign.NBrings)
                {
                    zi = piaacmcopticaldesign.NBrings;
                }

                long ring = piaacmcopticaldesign.NBrings -
                            zi; // 0 for inner disk, increases outward
                if(zi == 0)
                {
                    ring = -1;
                }

                long zoneindex;
                if(piaacmcparams.PIAACMC_FPMsectors == 0)
                {
                    zoneindex = piaacmcopticaldesign.NBrings - zi + 1;
                    if(zi == 0)
                    {
                        zoneindex = 0;
                    }
                }
                else
                {
                    if(ring == -1)
                    {
                        zoneindex = 0;
                    }
                    else
                    {
                        if(ring == 0)  // inner disk
                        {
                            zoneindex = 1;
                        }
                        else
                        {
                            zoneindex = nbsectorcumul[ring - 1] + 1;
                            zoneindex += (long)(PAf * nbsector[ring]);
                        }

                        if(piaacmcopticaldesign.NBrings > 1)
                        {
                            if(ring < piaacmcopticaldesign.NBringCentCone)
                            {
                                zoneindex = 0;
                            }
                            else
                            {
                                zoneindex -= nbzonescc;
                            }
                        }
                    }
                }

                data.image[ID]
                .array.UI16[jj * piaacmcopticaldesign.fpmarraysize + ii] =
                    zoneindex;
            }

        if(piaacmcparams.PIAACMC_FPMsectors == 0)
        {
            if(piaacmcopticaldesign.NBrings > 1)
            {
                piaacmcopticaldesign.focmNBzone =
                    piaacmcopticaldesign.NBrings - nbzonescc;
            }
            else
            {
                piaacmcopticaldesign.focmNBzone = 1;
            }
        }
        else
        {
            if(piaacmcopticaldesign.NBrings > 1)
            {
                piaacmcopticaldesign.focmNBzone =
                    nbsectorcumul[piaacmcopticaldesign.NBrings - 1] - nbzonescc;
            }
            else
            {
                piaacmcopticaldesign.focmNBzone = 1;
            }
        }
    }

    if(piaacmcparams.PIAACMC_FPMsectors == 2)
    {
        hexsteppix = 0.5 * piaacmcopticaldesign.fpmarraysize /
                     piaacmcopticaldesign.NBrings *
                     piaacmcparams.FPMSCALEFACTOR;
        for(hindex = 0; hindex < hindexMax; hindex++)
        {
            ID1 = make_hexagon("_TMPhex",
                               piaacmcopticaldesign.fpmarraysize,
                               piaacmcopticaldesign.fpmarraysize,
                               0.5 * piaacmcopticaldesign.fpmarraysize +
                               hex_x[hindex] * hexsteppix,
                               0.5 * piaacmcopticaldesign.fpmarraysize +
                               hex_y[hindex] * hexsteppix,
                               hexsteppix * (1.0 - hexgap) * (sqrt(3.0) / 2.0));

            for(long ii = 0; ii < piaacmcopticaldesign.fpmarraysize *
                    piaacmcopticaldesign.fpmarraysize;
                    ii++)
            {
                if(data.image[ID1].array.F[ii] > 0.5)
                {
                    data.image[ID].array.UI16[ii] =
                        (unsigned int) hex_number[hindex] + 1;
                }
            }
            delete_image_ID("_TMPhex", DELETE_IMAGE_ERRMODE_WARNING);
        }

        if(piaacmcopticaldesign.NBrings > 1)
        {
            piaacmcopticaldesign.focmNBzone =
                nbsectorcumul[piaacmcopticaldesign.NBrings - 1];
        }
        else
        {
            piaacmcopticaldesign.focmNBzone = 1;
        }
    }

    printf(
        "[%d] piaacmcopticaldesign.focmNBzone  =  %ld %ld    %ld %ld   ->  %ld "
        "  (%ld)\n",
        piaacmcparams.PIAACMC_FPMsectors,
        piaacmcopticaldesign.NBrings,
        nbsectorcumul[piaacmcopticaldesign.NBrings - 1],
        piaacmcopticaldesign.NBringCentCone,
        nbzonescc,
        piaacmcopticaldesign.focmNBzone,
        piaacmcopticaldesign.NBrings);

    if(piaacmcparams.PIAACMC_FPMsectors != 0)
    {
        printf("Saving %s ....\n", IDname);
        save_fits(IDname, "__test_zonemap_00.fits"); //TEST
        // sleep(100000);
    }

    free(nbsector);
    free(nbsectorcumul);

    printf("NUMBER OF ZONES = %ld\n", piaacmcopticaldesign.focmNBzone); //TEST

    if(outID != NULL)
    {
        *outID = ID;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
