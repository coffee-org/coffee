/**
 * @file mkPIAAMshapes_from_RadSag.c
 * @brief PIAA-type coronagraph design, make PIAA shapes from radial sag
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>

# ifdef _OPENMP
# include <omp.h>
# endif


#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"




/**
 * @brief Make PIAA OPD screens from radial sag profile
 *
 * @param[in] piaa1Dsagfname
 *
 * @param[in] PIAAsep
 *      Physical separation between PIAA surfaces [m]
 *
 * @param[in] beamrad
 *      Beam radius [m]
 *
 * @param[in] r0limfact
 *      Input beam extrapolation radius limit
 *
 * @param[in] r1limfact
 *      Output beam extrapolation radius limit
 *
 * @param size
 *      Output image linear size
 *
 * @param beamradpix
 *      Beam radius in pixel
 *
 * @param NBradpts
 *      Number of sample points
 *
 * @param ID_PIAAM0_name
 *
 * @param ID_PIAAM1_name
 *
 * @return errno_t
 */
errno_t mkPIAAMshapes_from_RadSag(
    const char *__restrict__ piaa1Dsagfname,
    double PIAAsep,
    double beamrad,
    double r0limfact,
    double r1limfact,
    uint32_t size,
    double beamradpix,
    long   NBradpts,
    const char *__restrict__ ID_PIAAM0_name,
    const char *__restrict__ ID_PIAAM1_name
)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %s %u %f %ld %s %s",
                     piaa1Dsagfname,
                     size,
                     beamradpix,
                     NBradpts,
                     ID_PIAAM0_name,
                     ID_PIAAM1_name
                    );

    imageID ID_PIAAM0, ID_PIAAM1;

    double *r0array;
    double *z0array;
    double *r1array;
    double *z1array;




    r0array = (double *) malloc(sizeof(double) * NBradpts);
    if(r0array == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }

    z0array = (double *) malloc(sizeof(double) * NBradpts);
    if(z0array == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }

    r1array = (double *) malloc(sizeof(double) * NBradpts);
    if(r1array == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }

    z1array = (double *) malloc(sizeof(double) * NBradpts);
    if(z1array == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort(); // or handle error in other ways
    }

    {
        // Read ASCII sag file
        DEBUG_TRACEPOINT_LOG("read sag radial profile %s", piaa1Dsagfname);

        FILE * fp = fopen(piaa1Dsagfname, "r");
        if(fp != NULL)
        {
            for(long k = 0; k < NBradpts; k++)
            {
                int ret = fscanf(fp, "%lf %lf %lf %lf\n",
                                 &r0array[k], &z0array[k], &r1array[k], &z1array[k]);
                (void) ret;
            }
            fclose(fp);
        }
        else
        {
            FUNC_RETURN_FAILURE("Cannot read file %s", piaa1Dsagfname);
        }

    }

    //  for(k=0;k<nbpt;k++)
    //  printf("%ld %.8lf %.8lf %.8lf %.8lf\n", k, r0array[k], z0array[k], r1array[k], z1array[k]);


    for(long k = 0; k < NBradpts; k++)
    {
        z1array[k] -= PIAAsep;
    }



    FUNC_CHECK_RETURN(
        create_2Dimage_ID(ID_PIAAM0_name, size, size, &ID_PIAAM0)
    );

    FUNC_CHECK_RETURN(
        create_2Dimage_ID(ID_PIAAM1_name, size, size, &ID_PIAAM1)
    );

    printf("\n\n");


# ifdef _OPENMP
    #pragma omp parallel
    {
# endif


# ifdef _OPENMP
        #pragma omp for
# endif
        for(long ii = 0; ii < size; ii++)
        {
            //      printf("\r %ld / %ld     ", ii, size);
            //fflush(stdout);


            for(long jj = 0; jj < size; jj++)
            {
                double x = (1.0 * ii - 0.5 * size) / beamradpix;
                double y = (1.0 * jj - 0.5 * size) / beamradpix;
                double r = sqrt(x * x + y * y) * beamrad;

                if(r < r0limfact * beamrad)
                {
                    long k = 1;
                    while((k < NBradpts - 2)
                            && (r0array[k] < r))
                    {
                        k++;
                    }
                    double r00 = r0array[k - 1];
                    double r01 = r0array[k];
                    double alpha = (r - r00) / (r01 - r00);
                    if(alpha > 1.0)
                    {
                        alpha = 1.0;
                    }
                    double val = (1.0 - alpha) * z0array[k - 1] + alpha * z0array[k];
                    data.image[ID_PIAAM0].array.F[jj * size + ii] = val;
                }
                else
                {
                    data.image[ID_PIAAM0].array.F[jj * size + ii] = 0.0;
                }

                if(r < r1limfact * beamrad)
                {
                    long k = 1;
                    while((k < NBradpts - 2)
                            && (r1array[k] < r))
                    {
                        k++;
                    }
                    double r00 = r1array[k - 1];
                    double r01 = r1array[k];
                    double alpha = (r - r00) / (r01 - r00);
                    if(alpha > 1.0)
                    {
                        alpha = 1.0;
                    }
                    double val = (1.0 - alpha) * z1array[k - 1] + alpha * z1array[k];
                    data.image[ID_PIAAM1].array.F[jj * size + ii] = -val; //-piaacmc[0].PIAAsep);
                }
                else
                {
                    data.image[ID_PIAAM1].array.F[jj * size + ii] = 0.0;
                }
            }
        }
# ifdef _OPENMP
    }
# endif


    {   //TEST
        char fnametest[STRINGMAXLEN_FILENAME];
        WRITE_FILENAME(fnametest, "test_piaam0z.fits");
        FUNC_CHECK_RETURN(save_fl_fits(ID_PIAAM0_name, fnametest));
    }

    printf("\n\n");

    free(r0array);
    free(z0array);
    free(r1array);
    free(z1array);

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}


