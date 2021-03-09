/**
 * @file    PIAAACMCsimul_mkFPM_zonemap.c
 * @brief   PIAA-type coronagraph design, initialize
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


/* =============================================================================================== */
/* =============================================================================================== */
/*                                        HEADER FILES                                             */
/* =============================================================================================== */
/* =============================================================================================== */

// System includes
#include <stdio.h>
#include <math.h>


// milk includes

#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "image_gen/image_gen.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"




/* =============================================================================================== */
/* =============================================================================================== */
/*                                  GLOBAL DATA DECLARATION                                        */
/* =============================================================================================== */
/* =============================================================================================== */

extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTSYST *optsyst;

extern OPTPIAACMCDESIGN *piaacmc;





/* =============================================================================================== */
/* =============================================================================================== */
/*                                    FUNCTIONS SOURCE CODE                                        */
/* =============================================================================================== */
/* =============================================================================================== */

/**
 * @param[out]  IDname  Name of output image
 */
long PIAACMCsimul_mkFPM_zonemap(
	const char *IDname
	)
{
    FILE *fp;
    char fname[500];
    long NBzones;
    long ID;
    double x, y, r, PA;
    uint_fast64_t ii, jj;
    uint_fast32_t zi;
    uint32_t *sizearray;

    uint_fast32_t ring;
    uint_fast32_t *nbsector;
    uint_fast32_t *nbsectorcumul;
    double PAf;
    double eps = 1.0e-6;

    uint_fast16_t zoneindex;
    uint_fast64_t cnt, cnt1;
    uint_fast32_t nbzonescc;

    double hexstep = 1.0;
    double hexgap = -0.0001;
    double hexsteppix;
    int_fast64_t ii1, jj1, ii1max, jj1max;
    double hx, hy;
    double hex_x[10000];
    double hex_y[10000];
    uint_fast32_t hex_ring[10000];
    uint_fast32_t hex_number[10000];
    uint_fast32_t hcnt;
    uint_fast32_t hindex;
    uint_fast32_t hindexMax;
    long ID1;
	
	
	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif
	
	printf("function PIAACMCsimul_mkFPM_zonemap\n");
	fflush(stdout);


    sizearray = (uint32_t*) malloc(sizeof(uint32_t)*2);
    sizearray[0] = piaacmc[0].fpmarraysize;
    sizearray[1] = piaacmc[0].fpmarraysize;
    ID = create_image_ID(IDname, 2, sizearray, _DATATYPE_UINT16, 0, 0);
    free(sizearray);

    nbsector = (uint_fast32_t*) malloc(sizeof(uint_fast32_t)*piaacmc[0].NBrings);
    nbsectorcumul = (uint_fast32_t*) malloc(sizeof(uint_fast32_t)*piaacmc[0].NBrings);


    switch ( piaacmcsimul_var.PIAACMC_FPMsectors ) {

    case 0: // rings
        nbsectorcumul[0] = 1;
        for(ring=1; ring<piaacmc[0].NBrings; ring++)
            nbsectorcumul[ring] = nbsectorcumul[ring-1] + 1;
        break;

    case 1: // rings broken in sectors
        sprintf(fname, "%s/fpm_zonescoord_%d_%03ld.txt", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.PIAACMC_FPMsectors, piaacmc[0].NBrings);
        fp = fopen(fname, "w");
        NBzones = 0;
        cnt = 0;
        fprintf(fp, "0 0\n");
        nbsector[0] = 1;
        nbsectorcumul[0] = 1;
        for(ring=1; ring<piaacmc[0].NBrings; ring++)
        {

            nbsector[ring] = 2*(ring+1);
            nbsectorcumul[ring] = nbsectorcumul[ring-1]+nbsector[ring];
            for(cnt1=0; cnt1<nbsector[ring]; cnt1++)
            {
                cnt++;
                fprintf(fp, "%ld %ld\n", cnt, ring);
            }
        }

        fclose(fp);
        for(ring=0; ring<piaacmc[0].NBrings; ring++)
            printf("ring %ld : %ld %ld\n", ring, nbsector[ring], nbsectorcumul[ring]);
        break;

    case 2: // rings of hexagons
         sprintf(fname, "%s/fpm_zonescoord_%d_%03ld.txt", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.PIAACMC_FPMsectors, piaacmc[0].NBrings);
        fp = fopen(fname, "w");
        fprintf(fp, "# focal plane mask zones geometry\n");
        fprintf(fp, "# hexagonal tiling\n");
        fprintf(fp, "# col 1: hexagon index\n");
        fprintf(fp, "# col 2: ring index\n");
        fprintf(fp, "# col 3: hexagon center x coordinate\n");
        fprintf(fp, "# col 4: hexagon center y coordinate\n");
        fprintf(fp, "# Note: unit = hexagon center to tip radius (circumradius)\n");
        fprintf(fp, "# \n");
        nbsector[0] = 1;
        nbsector[0] = 1;
        nbsectorcumul[0] = 1;
       for(ring=1; ring<piaacmc[0].NBrings; ring++)
        {
            nbsector[ring] = 0;
            nbsectorcumul[ring] = 0;
        }
        hindex = 0;
        // hegagon side = s = ring unit
        ii1max = (long) (piaacmc[0].NBrings/3+2);
        jj1max = (long) (piaacmc[0].NBrings/sqrt(3.0)+2);
        for(ii1=-ii1max; ii1<ii1max; ii1++)
            for(jj1=-jj1max; jj1<jj1max; jj1++)
            {
                hx = hexstep*ii1*3;
                hy = hexstep*sqrt(3.0)*jj1;
                ring = (long) sqrt(hx*hx+hy*hy);
                if(ring<piaacmc[0].NBrings)
                {
                    nbsector[ring] ++;
                    hex_x[hindex] = hx;
                    hex_y[hindex] = hy;
                    hex_ring[hindex] = ring;
                    hindex++;
                }


                hx += hexstep*1.5;
                hy += hexstep*sqrt(3.0)/2.0;
                ring = (long) sqrt(hx*hx+hy*hy);
                if(ring<piaacmc[0].NBrings)
                {
                    nbsector[ring] ++;
                    hex_x[hindex] = hx;
                    hex_y[hindex] = hy;
                    hex_ring[hindex] = ring;
                    hindex++;
                }
            }
        hindexMax = hindex;
        
        fprintf(fp, "%5ld %5ld  %11.6f %11.6f\n", (long) 0, (long) 0, 0.0, 0.0);
        hcnt = 1;
        for(ring=1; ring<piaacmc[0].NBrings; ring++)
        {
            for(hindex=0; hindex<hindexMax; hindex++)
                if(hex_ring[hindex]==ring)
                {
                    hex_number[hindex] = hcnt;
                    fprintf(fp, "%5ld %5ld  %11.6f %11.6f\n", hcnt, ring, hex_x[hindex], hex_y[hindex]);
                    hcnt++;
                }
            if(ring>0)
                nbsectorcumul[ring] = nbsectorcumul[ring-1]+nbsector[ring];
        }
        fclose(fp);
         for(ring=0; ring<piaacmc[0].NBrings; ring++)
            printf("ring %ld : %ld %ld\n", ring, nbsector[ring], nbsectorcumul[ring]);
        break;

    default:
        printf("ERROR: FPMsector mode (%d) not recognized\n", piaacmcsimul_var.PIAACMC_FPMsectors);
        exit(0);
        break;
    }



    if(piaacmc[0].NBringCentCone>0)
        nbzonescc = nbsectorcumul[piaacmc[0].NBringCentCone-1];
    else
        nbzonescc = 0;



    //  ID = create_2Dimage_ID(IDname, piaacmc[0].fpmarraysize, piaacmc[0].fpmarraysize);


    if( ( piaacmcsimul_var.PIAACMC_FPMsectors == 0 ) || ( piaacmcsimul_var.PIAACMC_FPMsectors == 1 ) )
    {
        for(ii=0; ii<piaacmc[0].fpmarraysize; ii++)
            for(jj=0; jj<piaacmc[0].fpmarraysize; jj++)
            {
                x = (2.0*ii-1.0*piaacmc[0].fpmarraysize)/piaacmc[0].fpmarraysize / piaacmcsimul_var.FPMSCALEFACTOR;
                y = (2.0*jj-1.0*piaacmc[0].fpmarraysize)/piaacmc[0].fpmarraysize / piaacmcsimul_var.FPMSCALEFACTOR;
                r = sqrt(x*x+y*y);
                PA = atan2(y,x);
                PAf = 0.5*((PA/M_PI)+1.0);
                if(PAf<eps)
                    PAf = eps;
                if(PAf>1.0-eps)
                    PAf = 1.0-eps;

                zi = (long) ceil((1.0-r)*piaacmc[0].NBrings);


                if(zi<0.1)
                    zi = 0;
                if(zi>piaacmc[0].NBrings)
                    zi = piaacmc[0].NBrings;

                ring = piaacmc[0].NBrings-zi; // 0 for inner disk, increases outward
                if(zi==0)
                    ring = -1;


                if( piaacmcsimul_var.PIAACMC_FPMsectors == 0 )
                {
                    zoneindex = (unsigned short int) (piaacmc[0].NBrings-zi+1);
                    if(zi==0)
                        zoneindex = 0;
                }
                else
                {
                    if(ring==-1)
                        zoneindex = 0;
                    else
                    {
                        if(ring==0) // inner disk
                            zoneindex = 1;
                        else
                        {
                            zoneindex = (unsigned short int) nbsectorcumul[ring-1]+1;
                            zoneindex += (unsigned short int) (PAf*nbsector[ring]);
                        }

                        if(piaacmc[0].NBrings>1)
                        {
                            if(ring<piaacmc[0].NBringCentCone)
                                zoneindex = 0;
                            else
                                zoneindex -= nbzonescc;
                        }
                    }
                }

                data.image[ID].array.UI16[jj*piaacmc[0].fpmarraysize+ii] = zoneindex;

            }

        if( piaacmcsimul_var.PIAACMC_FPMsectors == 0 )
        {
            if(piaacmc[0].NBrings>1)
                piaacmc[0].focmNBzone = piaacmc[0].NBrings - nbzonescc;
            else
                piaacmc[0].focmNBzone = 1;
        }
        else
        {
            if(piaacmc[0].NBrings>1)
                piaacmc[0].focmNBzone = nbsectorcumul[piaacmc[0].NBrings-1] - nbzonescc;
            else
                piaacmc[0].focmNBzone = 1;
        }
    }


    if( piaacmcsimul_var.PIAACMC_FPMsectors == 2 )
        {
            hexsteppix = 0.5*piaacmc[0].fpmarraysize/piaacmc[0].NBrings * piaacmcsimul_var.FPMSCALEFACTOR;
            for(hindex=0;hindex<hindexMax;hindex++)
            {
                ID1 = make_hexagon("_TMPhex", piaacmc[0].fpmarraysize, piaacmc[0].fpmarraysize, 0.5*piaacmc[0].fpmarraysize + hex_x[hindex]*hexsteppix, 0.5*piaacmc[0].fpmarraysize + hex_y[hindex]*hexsteppix, hexsteppix*(1.0-hexgap)*(sqrt(3.0)/2.0));
                
                for(ii=0; ii<piaacmc[0].fpmarraysize*piaacmc[0].fpmarraysize; ii++)                   
                        {
                            if(data.image[ID1].array.F[ii] > 0.5)
                                data.image[ID].array.UI16[ii] = (unsigned int) hex_number[hindex]+1;
                        }
                delete_image_ID("_TMPhex");
            }
            
            if(piaacmc[0].NBrings>1)
                piaacmc[0].focmNBzone = nbsectorcumul[piaacmc[0].NBrings-1];
            else
                piaacmc[0].focmNBzone = 1;
        }


    printf("[%d] piaacmc[0].focmNBzone  =  %ld %ld    %ld %ld   ->  %ld   (%ld)\n", piaacmcsimul_var.PIAACMC_FPMsectors, piaacmc[0].NBrings, nbsectorcumul[piaacmc[0].NBrings-1], piaacmc[0].NBringCentCone, nbzonescc, piaacmc[0].focmNBzone, piaacmc[0].NBrings);


   if( piaacmcsimul_var.PIAACMC_FPMsectors != 0 )
    {
      printf("Saving %s ....\n", IDname);
      save_fits(IDname, "!__test_zonemap_00.fits"); //TEST
       // sleep(100000);
    }

    free(nbsector);
    free(nbsectorcumul);

	printf("NUMBER OF ZONES = %ld\n", piaacmc[0].focmNBzone); //TEST


    return ID;
}
