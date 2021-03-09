/**
 * @file    PIAACMCsimul_mkPIAAMshapes_from_RadSag.c
 * @brief   PIAA-type coronagraph design, make PIAA shapes from radial sag
 * 
 * Can design both APLCMC and PIAACMC coronagraphs
 *  
 * 
 * ## Change log
 * - 20180323  Guyon   minor syntax cleanup  
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
#include <malloc.h>


#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "PIAACMCsimul/PIAACMCsimul.h"




/* =============================================================================================== */
/* =============================================================================================== */
/*                                  GLOBAL DATA DECLARATION                                        */
/* =============================================================================================== */
/* =============================================================================================== */

extern OPTPIAACMCDESIGN *piaacmc;







/* =============================================================================================== */
/* =============================================================================================== */
/*                                    FUNCTIONS SOURCE CODE                                        */
/* =============================================================================================== */
/* =============================================================================================== */


/**
 * @brief Make PIAA OPD screens from radial sag profile
 * 
 */
 
int PIAACMCsimul_mkPIAAMshapes_from_RadSag(
	const char *fname, 
	const char *ID_PIAAM0_name, 
	const char *ID_PIAAM1_name
	)
{
    FILE *fp;
    long size;
    long ii, jj;
    long ID_PIAAM0, ID_PIAAM1;

    long k;

    double x, y, r, r1;

    double *r0array;
    double *z0array;
    double *r1array;
    double *z1array;

    double alpha;
    double r00, r01;
    double val;

    double beamradpix;
    int ret;

	#ifdef PIAASIMUL_LOGFUNC0
		PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__, "");
	#endif


    size = piaacmc[0].size;
    beamradpix = piaacmc[0].beamrad/piaacmc[0].pixscale;
    printf("SIZE = %ld, beamrad = %f pix, sep = %f m\n", size, beamradpix, piaacmc[0].PIAAsep);
    fflush(stdout);



    r0array = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
    z0array = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
    r1array = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);
    z1array = (double*) malloc(sizeof(double)*piaacmc[0].NBradpts);

    fp = fopen(fname, "r");
    for(k=0; k<piaacmc[0].NBradpts; k++)
        ret = fscanf(fp,"%lf %lf %lf %lf\n", &r0array[k], &z0array[k], &r1array[k], &z1array[k]);
    fclose(fp);

    //  for(k=0;k<nbpt;k++)
    //  printf("%ld %.8lf %.8lf %.8lf %.8lf\n", k, r0array[k], z0array[k], r1array[k], z1array[k]);


    for(k=0; k<piaacmc[0].NBradpts; k++)
        z1array[k] -= piaacmc[0].PIAAsep;




    ID_PIAAM0 = create_2Dimage_ID(ID_PIAAM0_name, size, size);
    ID_PIAAM1 = create_2Dimage_ID(ID_PIAAM1_name, size, size);

    printf("\n\n");

# ifdef HAVE_LIBGOMP
    #pragma omp parallel default(shared) private(ii, jj, x, y, r, k, r00, r01, alpha, val)
    {
# endif


# ifdef HAVE_LIBGOMP
        #pragma omp for
# endif
        for(ii=0; ii<size; ii++)
        {
            //      printf("\r %ld / %ld     ", ii, size);
            //fflush(stdout);


            for(jj=0; jj<size; jj++)
            {
                x = (1.0*ii-0.5*size)/beamradpix;
                y = (1.0*jj-0.5*size)/beamradpix;
                r = sqrt(x*x+y*y)*piaacmc[0].beamrad;

                if(r<piaacmc[0].r0lim*piaacmc[0].beamrad)
                {
                    k = 1;
                    while((r0array[k]<r)&&(k<piaacmc[0].NBradpts-2))
                        k++;
                    r00 = r0array[k-1];
                    r01 = r0array[k];
                    alpha = (r-r00)/(r01-r00);
                    if(alpha>1.0)
                        alpha = 1.0;
                    val = (1.0-alpha)*z0array[k-1] + alpha*z0array[k];
                    data.image[ID_PIAAM0].array.F[jj*size+ii] = val;
                }
                else
                    data.image[ID_PIAAM0].array.F[jj*size+ii] = 0.0;

                if(r<piaacmc[0].r1lim*piaacmc[0].beamrad)
                {
                    k = 1;
                    while((r1array[k]<r)&&(k<piaacmc[0].NBradpts-2))
                        k++;
                    r00 = r1array[k-1];
                    r01 = r1array[k];
                    alpha = (r-r00)/(r01-r00);
                    if(alpha>1.0)
                        alpha = 1.0;
                    val = (1.0-alpha)*z1array[k-1] + alpha*z1array[k];
                    data.image[ID_PIAAM1].array.F[jj*size+ii] = -val;//-piaacmc[0].PIAAsep);
                }
                else
                    data.image[ID_PIAAM1].array.F[jj*size+ii] = 0.0;
            }
        }
# ifdef HAVE_LIBGOMP
    }
# endif



    printf("\n\n");

    free(r0array);
    free(z0array);
    free(r1array);
    free(z1array);

    return 0;
}


