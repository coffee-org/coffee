/**
 * @file    PIAACMCsimul_makePIAAshapes.c
 * @brief   PIAA-type coronagraph design, make PIAA shapes
 *
 * Can design both APLCMC and PIAACMC coronagraphs
 *
 */


#include <stdlib.h>
#include <stdio.h>



#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"

#include "linopt_imtools/linopt_imtools.h"
#include "OpticsMaterials/OpticsMaterials.h"
#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"




extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTPIAACMCDESIGN *piaacmc;



///
/// mirror sag shapes:  piaam0z, piaam1z
///
/// piaa0Cmodescoeff -> piaa0Cz
/// piaa0Fmodescoeff -> piaa0Fz
/// piaa0Cz + piaa0Fz -> piaam0z
///

int PIAACMCsimul_makePIAAshapes(
    OPTPIAACMCDESIGN *design,
    long index
)
{
#ifdef PIAASIMUL_LOGFUNC0
    PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__,
                                 "");
#endif


    long size = piaacmc[0].size;

    // ============ construct PIAA shapes from fitting coefficients ==================



    if(piaacmc[0].PIAAmode == 1)
    {

        piaacmcsimul_var.MAKE_PIAA0shape = 0;
        if(piaacmcsimul_var.FORCE_MAKE_PIAA0shape == 0)
        {
            imageID ID = image_ID("piaam0z");
            if(ID == -1)
            {
                piaacmcsimul_var.MAKE_PIAA0shape = 1;
            }
        }
        else
        {
            piaacmcsimul_var.MAKE_PIAA0shape = 1;
        }


        piaacmcsimul_var.MAKE_PIAA1shape = 0;
        if(piaacmcsimul_var.FORCE_MAKE_PIAA1shape == 0)
        {
            imageID ID = image_ID("piaam1z");
            if(ID == -1)
            {
                piaacmcsimul_var.MAKE_PIAA1shape = 1;
            }
        }
        else
        {
            piaacmcsimul_var.MAKE_PIAA1shape = 1;
        }






        if(piaacmcsimul_var.MAKE_PIAA0shape == 1)
        {
            // assemble piaa0z and piaa1z images
            imageID ID0 = linopt_imtools_image_construct("Cmodes", "piaa0Cmodescoeff", "piaa0Cz");
            imageID ID1 = linopt_imtools_image_construct("Fmodes", "piaa0Fmodescoeff", "piaa0Fz");
            imageID ID = image_ID("piaam0z");
            if(ID == -1)
            {
                ID = create_2Dimage_ID("piaam0z", size, size);
            }

            uint32_t size0 = data.image[ID0].md[0].size[0];
            uint32_t size1 = data.image[ID1].md[0].size[0];
            for(long ii = 0; ii < size * size; ii++)
            {
                data.image[ID].array.F[ii] = 0.0;
            }

            printf("========================== STEP 01a  %ld %ld %ld\n", ID, ID0, ID1);


            for(uint32_t ii = 0; ii < size0; ii++)
                for(uint32_t jj = 0; jj < size0; jj++)
                {
                    data.image[ID].array.F[(jj + (size - size0) / 2)*size + (ii +
                                           (size - size0) / 2)] += data.image[ID0].array.F[jj * size0 + ii];
                }
            for(uint32_t ii = 0; ii < size1; ii++)
                for(uint32_t jj = 0; jj < size1; jj++)
                {
                    data.image[ID].array.F[(jj + (size - size1) / 2)*size + (ii +
                                           (size - size1) / 2)] += data.image[ID1].array.F[jj * size1 + ii];
                }

            if(piaacmcsimul_var.PIAACMC_save == 1)
            {
                char fname[STRINGMAXLEN_FULLFILENAME];

                WRITE_FULLFILENAME(fname, "%s/piaa0Cz.fits", piaacmcsimul_var.piaacmcconfdir);
                save_fits("piaa0Cz", fname);

                WRITE_FULLFILENAME(fname, "%s/piaa0Fz.fits", piaacmcsimul_var.piaacmcconfdir);
                save_fits("piaa0Fz", fname);

                WRITE_FULLFILENAME(fname, "%s/piaam0z.fits", piaacmcsimul_var.piaacmcconfdir);
                save_fits("piaam0z", fname);
            }
            delete_image_ID("piaa0Cz", DELETE_IMAGE_ERRMODE_WARNING);
            delete_image_ID("piaa0Fz", DELETE_IMAGE_ERRMODE_WARNING);

            imageID IDpiaam0z = ID;

            // make lense shapes if applicable
            if(design[index].PIAAmaterial_code != 0) // refractive PIAA
            {
                // if piaar0zsag does not exist or is wrong size, create it
                imageID IDpiaar0zsag = image_ID("piaar0zsag");
                int mkpiaar0zsag = 0;
                if(IDpiaar0zsag == -1)
                {
                    mkpiaar0zsag = 1;
                }
                else
                {
                    if((data.image[IDpiaar0zsag].md[0].size[0] != size)
                            || (data.image[IDpiaar0zsag].md[0].size[1] !=
                                size)) //||(data.image[IDpiaar0zsag].md[0].size[2] != design[index].nblambda))
                    {
                        delete_image_ID("piaar0zsag", DELETE_IMAGE_ERRMODE_WARNING);
                        mkpiaar0zsag = 1;
                    }
                }
                if(mkpiaar0zsag == 1)
                {
                    IDpiaar0zsag = create_2Dimage_ID("piaar0zsag", size,
                                                     size);    //, design[index].nblambda);
                }


                FILE *fpri;
                {
                    char fname[STRINGMAXLEN_FULLFILENAME];
                    WRITE_FULLFILENAME(fname, "%s/ri_array.txt", piaacmcsimul_var.piaacmcconfdir);
                    if((fpri = fopen(fname, "w")) == NULL)
                    {
                        printf("ERROR: cannot open file \"%s\"\n", fname);
                        exit(0);
                    }
                }
                double ri0 = OpticsMaterials_n(design[index].PIAAmaterial_code,
                                        design[index].lambda); // refractive index at central lambda
                double sag2opd_coeff0 = (ri0 - 1.0) / 2.0;
                for(long k = 0; k < design[index].nblambda; k++)
                {
                    // sag to OPD coeff
                    double ri = OpticsMaterials_n(design[index].PIAAmaterial_code,
                                           piaacmc[0].lambdaarray[k]); // refractive index
                    double sag2opd_coeff = (ri - 1.0) / 2.0;
                    fprintf(fpri, "%g %.16f %.16f %.16f %.16f\n", piaacmc[0].lambdaarray[k], ri,
                            ri0, sag2opd_coeff, sag2opd_coeff / sag2opd_coeff0);
                    //                for(ii=0; ii<size*size; ii++)
                    //                  data.image[IDpiaar0zsag].array.F[k*size*size+ii] = data.image[IDpiaam0z].array.F[ii] * sag2opd_coeff/sag2opd_coeff0; //sag2opd_coeff * data.image[IDpiaam0z].array.F[ii] / sag2opd_coeff0;
                }
                fclose(fpri);

                for(long ii = 0; ii < size * size; ii++)
                {
                    data.image[IDpiaar0zsag].array.F[ii] =
                        data.image[IDpiaam0z].array.F[ii] / sag2opd_coeff0;
                }

                {
                    char fname[STRINGMAXLEN_FULLFILENAME];
                    WRITE_FULLFILENAME(fname, "%s/piaar0zsag.fits", piaacmcsimul_var.piaacmcconfdir);
                    if(piaacmcsimul_var.PIAACMC_save == 1)
                    {
                        save_fl_fits("piaar0zsag", fname);
                    }
                    printf("Saved piaar0zsag to %s\n", fname);
                }
            }

        }



        if(piaacmcsimul_var.MAKE_PIAA1shape == 1)
        {
            imageID ID0 = linopt_imtools_image_construct("Cmodes", "piaa1Cmodescoeff", "piaa1Cz");
            imageID ID1 = linopt_imtools_image_construct("Fmodes", "piaa1Fmodescoeff", "piaa1Fz");
            imageID ID = image_ID("piaam1z");
            if(ID == -1)
            {
                ID = create_2Dimage_ID("piaam1z", size, size);
            }
            for(long ii = 0; ii < size * size; ii++)
            {
                data.image[ID].array.F[ii] = 0.0;
            }
            uint32_t size0 = data.image[ID0].md[0].size[0];
            uint32_t size1 = data.image[ID1].md[0].size[0];
            for(uint32_t ii = 0; ii < size0; ii++)
                for(uint32_t jj = 0; jj < size0; jj++)
                {
                    data.image[ID].array.F[(jj + (size - size0) / 2)*size + (ii +
                                           (size - size0) / 2)] += data.image[ID0].array.F[jj * size0 + ii];
                }
            for(uint32_t ii = 0; ii < size1; ii++)
                for(uint32_t jj = 0; jj < size1; jj++)
                {
                    data.image[ID].array.F[(jj + (size - size1) / 2)*size + (ii +
                                           (size - size1) / 2)] += data.image[ID1].array.F[jj * size1 + ii];
                }

            if(piaacmcsimul_var.PIAACMC_save == 1)
            {
                char fname[STRINGMAXLEN_FULLFILENAME];

                WRITE_FULLFILENAME(fname, "%s/piaa1Cz.fits", piaacmcsimul_var.piaacmcconfdir);
                save_fits("piaa1Cz", fname);

                WRITE_FULLFILENAME(fname, "%s/piaa1Fz.fits", piaacmcsimul_var.piaacmcconfdir);
                save_fits("piaa1Fz", fname);

                WRITE_FULLFILENAME(fname, "%s/piaam1z.fits", piaacmcsimul_var.piaacmcconfdir);
                save_fits("piaam1z", fname);
            }
            delete_image_ID("piaa1Cz", DELETE_IMAGE_ERRMODE_WARNING);
            delete_image_ID("piaa1Fz", DELETE_IMAGE_ERRMODE_WARNING);

            imageID IDpiaam1z = ID;


            // make lense shapes if applicable
            if(design[index].PIAAmaterial_code != 0) // refractive PIAA
            {

                // if piaar1zsag does not exist or is wrong size, create it
                imageID IDpiaar1zsag = image_ID("piaar1zsag");
                int mkpiaar1zsag = 0;
                if(IDpiaar1zsag == -1)
                {
                    mkpiaar1zsag = 1;
                }
                else
                {
                    if((data.image[IDpiaar1zsag].md[0].size[0] != size)
                            || (data.image[IDpiaar1zsag].md[0].size[1] !=
                                size)) //||(data.image[IDpiaar1zsag].md[0].size[2] != design[index].nblambda))
                    {
                        delete_image_ID("piaar1zsag", DELETE_IMAGE_ERRMODE_WARNING);
                        mkpiaar1zsag = 1;
                    }
                }
                if(mkpiaar1zsag == 1)
                {
                    IDpiaar1zsag = create_2Dimage_ID("piaar1zsag", size,
                                                     size);    //, design[index].nblambda);
                }

                FILE *fpri;
                {
                    char fname[STRINGMAXLEN_FULLFILENAME];
                    WRITE_FULLFILENAME(fname, "%s/ri_array.txt", piaacmcsimul_var.piaacmcconfdir);
                    if((fpri = fopen(fname, "w")) == NULL)
                    {
                        printf("ERROR: cannot open file \"%s\"\n", fname);
                        exit(0);
                    }
                }
                double ri0 = OpticsMaterials_n(design[index].PIAAmaterial_code,
                                        design[index].lambda); // refractive index at central lambda
                double sag2opd_coeff0 = (ri0 - 1.0) / 2.0;
                for(long k = 0; k < piaacmc[0].nblambda; k++)
                {
                    // sag to OPD coeff
                    double ri = OpticsMaterials_n(design[index].PIAAmaterial_code,
                                           piaacmc[0].lambdaarray[k]); // refractive index
                    double sag2opd_coeff = (ri - 1.0) / 2.0;
                    fprintf(fpri, "%g %.16f %.16f %.16f %.16f\n", piaacmc[0].lambdaarray[k], ri,
                            ri0, sag2opd_coeff, sag2opd_coeff / sag2opd_coeff0);
                    // for(ii=0; ii<size*size; ii++)
                    //    data.image[IDpiaar1zsag].array.F[k*size*size+ii] = sag2opd_coeff * data.image[IDpiaam1z].array.F[ii] / sag2opd_coeff0;
                }
                fclose(fpri);

                for(long ii = 0; ii < size * size; ii++)
                {
                    data.image[IDpiaar1zsag].array.F[ii] = data.image[IDpiaam1z].array.F[ii] /
                                                           sag2opd_coeff0;
                }

                {
                    char fname[STRINGMAXLEN_FULLFILENAME];
                    WRITE_FULLFILENAME(fname, "%s/piaar1zsag.fits", piaacmcsimul_var.piaacmcconfdir);
                    if(piaacmcsimul_var.PIAACMC_save == 1)
                    {
                        save_fl_fits("piaar1zsag", fname);
                    }
                    printf("Saved piaar1zsag to %s\n", fname);
                }
            }
        }
    }

    return RETURN_SUCCESS;
}

