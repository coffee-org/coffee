/**
 * @file    PIAACMCsimul_geomProp.c
 * @brief   Beam geometrical propagation based on local slopes
 * 
 * Uses local slopes of sag map to propagate (un-corrected intensity only) a beam
 *  
 * @author  O. Guyon
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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>



// milk includes
#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "image_filter/image_filter.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"


extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTSYST *optsyst;

extern OPTPIAACMCDESIGN *piaacmc;





/* =============================================================================================== */
/* =============================================================================================== */
/*                                    FUNCTION(S) SOURCE CODE                                      */
/* =============================================================================================== */
/* =============================================================================================== */

/**
 * @brief Lyot stops positions from zmin to zmax relative to current, working back (light goes from 0 to zmax)
 * 
 * @param[in]	IDin_name		image : Input 2D intensity
 * @param[in]	IDsag_name	    image : 2D sag
 * @param[out]  IDout_name		image : 2D propagated intensity
 * @param[out]  IDoutcnt_name	image : 2D rays counter
 * @param[in]   drindex         float : refractive index (2 for mirror)
 * @param[in]	pscale			float : pixel scale [m]
 * @param[in]   zprop          	float : propagation distance [m]
 * @param[in]	krad			float : kernel radius used to evaluate slope [pixel]
 * @param[in]	kstep			float : step size in input pupil [pixel]
 * @param[in]	rlim			float : clear aperture radius (don't compute outside this value) [pixel]
*/

long PIAACMCsimul_geomProp(
	const char *IDin_name,			
	const char *IDsag_name, 
	const char *IDout_name,
	const char *IDoutcnt_name, 
	float drindex,
	float pscale, 
	float zprop, 
	float krad,
	float kstep,
	float rlim 
	)
{


	long IDin = image_ID(IDin_name);
	long xsize = data.image[IDin].md[0].size[0];
	long ysize = data.image[IDin].md[0].size[1];
	
	long IDsag = image_ID(IDsag_name);

	long IDout = create_2Dimage_ID(IDout_name, xsize, ysize);
	long IDoutcnt = create_2Dimage_ID(IDoutcnt_name, xsize, ysize);
	
	printf("kstep = %f\n", kstep);
	
	float x, y;
	for ( x = 0.5*xsize-rlim; x < 0.5*xsize+rlim; x += kstep )
		for ( y = 0.5*ysize-rlim; y < 0.5*ysize+rlim; y += kstep )
		{
			long ii0 = (long) (x+0.5);
			long jj0 = (long) (y+0.5);
			
			double sumvx = 0.0;
			double sumvy = 0.0;
			double sumvx0 = 0.0;
			double sumvy0 = 0.0;
			double sumv = 0.0;
			double sumc = 0.0;
			
			long ii, jj;
			long iimin, iimax, jjmin, jjmax;
			iimin = ii0 - ((long) krad+1);
			iimax = ii0 + ((long) krad+2);
			jjmin = jj0 - ((long) krad+1);
			jjmax = jj0 + ((long) krad+2);
			
			if(iimin<0){
				printf("ERROR %s line %d : iimin = %ld < 0  ii0 = %ld\n", __FILE__, __LINE__, iimin, ii0);
				exit(0);
			}
			if(iimax>xsize-1){
				printf("ERROR %s line %d : iimax = %ld > %ld  ii0 = %ld\n", __FILE__, __LINE__, iimax, xsize-1, ii0);
				exit(0);
			}
			
			if(jjmin<0){
				printf("ERROR %s line %d : jjmin = %ld < 0  jj0 = %ld\n", __FILE__, __LINE__, jjmin, jj0);
				exit(0);
			}
			if(jjmax>ysize-1){
				printf("ERROR %s line %d : jjmax = %ld > %ld  jj0 = %ld\n", __FILE__, __LINE__, jjmax, ysize-1, jj0);
				exit(0);
			}

			for (ii = iimin; ii < iimax; ii++)
				for (jj = jjmin; jj < jjmax; jj++)
				{
					float dx = 1.0*ii - x;
					float dy = 1.0*jj - y;
					
					float dr = sqrt(dx*dx+dy*dy)/krad;
					if(dr>1.0)
						dr = 1.0;
					float coeff = 0.5 + 0.5*cos(dr*M_PI);
					
					sumv += coeff * data.image[IDsag].array.F[jj*xsize+ii];
					sumvx += coeff * dx * data.image[IDsag].array.F[jj*xsize+ii];
					sumvy += coeff * dy * data.image[IDsag].array.F[jj*xsize+ii];
					sumvx0 += coeff * dx;
					sumvy0 += coeff * dy;
					sumc += coeff;
				}
			
			// local slopes [unitless]
			float slx = ( (sumvx / sumc) - (sumvx0 * (sumv/sumc) / sumc) ) / pscale;
			float sly = ( (sumvy / sumc) - (sumvy0 * (sumv/sumc) / sumc) ) / pscale; 
			
			// displacements [pix]
			float dii = (slx * zprop * drindex) / pscale;
			float djj = (sly * zprop * drindex) / pscale;

			long ii1 = (long) ( x + dii + 0.5 );
			long jj1 = (long) ( y + djj + 0.5 );
			
			if((ii1>0)&&(ii1<xsize)&&(jj1>0)&&(jj1<ysize))
			{			
				data.image[IDout].array.F[jj1*xsize + ii1] += data.image[IDin].array.F[jj0*xsize + ii0];
				data.image[IDoutcnt].array.F[jj1*xsize + ii1] += 1.0;
			}
		}
	
	long ii1, jj1;
	for(ii1 = 0; ii1 < xsize; ii1++)
		for(jj1 = 0; jj1 < ysize; jj1++)
			{
				if (data.image[IDoutcnt].array.F[jj1*xsize + ii1] > 0.1)
					data.image[IDout].array.F[jj1*xsize + ii1] /= data.image[IDoutcnt].array.F[jj1*xsize + ii1];
			}
		

    return(IDout);
}



