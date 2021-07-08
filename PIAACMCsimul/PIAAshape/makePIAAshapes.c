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

#include "PIAACMCsimul.h"






//OPTPIAACMCDESIGN *design,

/** @brief Make PIAA shapes
 *
 * Creates mirror sag shapes:  piaam0z, piaam1z
 *  and lens shapes if applicable
 *
 * Shapes are constructed from mode coefficients
 *
 * piaa0Cmodescoeff -> piaa0Cz
 * piaa0Fmodescoeff -> piaa0Fz
 * piaa0Cz + piaa0Fz -> piaam0z
 *
 */
errno_t makePIAAshapes()
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("PIAA material code = %d", piaacmcopticaldesign.PIAAmaterial_code);



    // ============ construct PIAA shapes from fitting coefficients ==================

    DEBUG_TRACEPOINT("PIAAmode %d", piaacmcopticaldesign.PIAAmode);
    if(piaacmcopticaldesign.PIAAmode == 1)
    {   // if using PIAA optics

        // check if PIAA shapes should be made
        //
        piaacmcparams.MAKE_PIAA0shape = 0;

        // PIAA shapes must be computed if one of the following conditions are met :
        // - FORCE set to 1
        // - piaam0z missing
        // - if refractive optics, piaar0zsag missing
        //
        if( (piaacmcparams.FORCE_MAKE_PIAA0shape == 1)
                || (image_ID("piaam0z") == -1)
                || ( (piaacmcopticaldesign.PIAAmaterial_code != 0) && (image_ID("piaar0zsag") == -1)))
        {
            piaacmcparams.MAKE_PIAA0shape = 1;
        }




        piaacmcparams.MAKE_PIAA1shape = 0;

        if( (piaacmcparams.FORCE_MAKE_PIAA1shape == 1)
                || (image_ID("piaam1z") == -1)
                || ( (piaacmcopticaldesign.PIAAmaterial_code != 0) && (image_ID("piaar1zsag") == -1)))
        {
            piaacmcparams.MAKE_PIAA1shape = 1;
        }



        DEBUG_TRACEPOINT("piaacmcparams.MAKE_PIAA0shape %d", piaacmcparams.MAKE_PIAA0shape);
        if(piaacmcparams.MAKE_PIAA0shape == 1)
        {   // if PIAA shapes should be computed

            // assemble piaa0z and piaa1z images from their linear coefficients

            // cosine radial modes
            imageID ID0;
            linopt_imtools_image_construct("Cmodes", "piaa0Cmodescoeff", "piaa0Cz", &ID0);

            // Fourier 2D modes
            imageID ID1;
            linopt_imtools_image_construct("Fmodes", "piaa0Fmodescoeff", "piaa0Fz", &ID1);

            imageID ID = image_ID("piaam0z");
            if(ID == -1)
            {
                FUNC_CHECK_RETURN(
                    create_2Dimage_ID("piaam0z", piaacmcopticaldesign.size, piaacmcopticaldesign.size, &ID)
                );
            }

            uint32_t size0 = data.image[ID0].md[0].size[0];
            uint32_t size1 = data.image[ID1].md[0].size[0];
            for(long ii = 0; ii < piaacmcopticaldesign.size * piaacmcopticaldesign.size; ii++)
            {
                data.image[ID].array.F[ii] = 0.0;
            }


            for(uint32_t ii = 0; ii < size0; ii++)
                for(uint32_t jj = 0; jj < size0; jj++)
                {
                    data.image[ID].array.F[(jj + (piaacmcopticaldesign.size - size0) / 2)*piaacmcopticaldesign.size + (ii +
                                           (piaacmcopticaldesign.size - size0) / 2)] += data.image[ID0].array.F[jj * size0 + ii];
                }
            for(uint32_t ii = 0; ii < size1; ii++)
                for(uint32_t jj = 0; jj < size1; jj++)
                {
                    data.image[ID].array.F[(jj + (piaacmcopticaldesign.size - size1) / 2)*piaacmcopticaldesign.size + (ii +
                                           (piaacmcopticaldesign.size - size1) / 2)] += data.image[ID1].array.F[jj * size1 + ii];
                }

            DEBUG_TRACEPOINT("piaacmcparams.PIAACMC_save %d", piaacmcparams.PIAACMC_save);
            if(piaacmcparams.PIAACMC_save == 1)
            {
                char fname[STRINGMAXLEN_FULLFILENAME];

                WRITE_FULLFILENAME(fname, "%s/mkPIAAshapes_Cmodes.fits", piaacmcparams.piaacmcconfdir);
                FUNC_CHECK_RETURN(save_fits("Cmodes", fname));

                WRITE_FULLFILENAME(fname, "%s/mkPIAAshapes_piaa0Cmodescoeff.fits", piaacmcparams.piaacmcconfdir);
                FUNC_CHECK_RETURN(save_fits("piaa0Cmodescoeff", fname));

                WRITE_FULLFILENAME(fname, "%s/mkPIAAshapes_piaa0Cz.fits", piaacmcparams.piaacmcconfdir);
                FUNC_CHECK_RETURN(save_fits("piaa0Cz", fname));

                WRITE_FULLFILENAME(fname, "%s/mkPIAAshapes_piaa0Fz.fits", piaacmcparams.piaacmcconfdir);
                FUNC_CHECK_RETURN(save_fits("piaa0Fz", fname));

                WRITE_FULLFILENAME(fname, "%s/mkPIAAshapes_piaam0z.fits", piaacmcparams.piaacmcconfdir);
                FUNC_CHECK_RETURN(save_fits("piaam0z", fname));
            }
            FUNC_CHECK_RETURN(delete_image_ID("piaa0Cz", DELETE_IMAGE_ERRMODE_IGNORE));
            FUNC_CHECK_RETURN(delete_image_ID("piaa0Fz", DELETE_IMAGE_ERRMODE_IGNORE));

            imageID IDpiaam0z = ID;

            // make lense shapes if applicable
            DEBUG_TRACEPOINT("PIAA material code = %d", piaacmcopticaldesign.PIAAmaterial_code);
            if(piaacmcopticaldesign.PIAAmaterial_code != 0)
            {   // refractive PIAA
                // if piaar0zsag does not exist or is wrong size, create it
                imageID IDpiaar0zsag = image_ID("piaar0zsag");
                int mkpiaar0zsag = 0;
                if(IDpiaar0zsag == -1)
                {
                    mkpiaar0zsag = 1;
                }
                else
                {
                    if((data.image[IDpiaar0zsag].md[0].size[0] != piaacmcopticaldesign.size)
                            || (data.image[IDpiaar0zsag].md[0].size[1] !=
                                piaacmcopticaldesign.size)) //||(data.image[IDpiaar0zsag].md[0].size[2] != design[index].nblambda))
                    {
                        FUNC_CHECK_RETURN(delete_image_ID("piaar0zsag", DELETE_IMAGE_ERRMODE_IGNORE));
                        mkpiaar0zsag = 1;
                    }
                }
                if(mkpiaar0zsag == 1)
                {
                    FUNC_CHECK_RETURN(
                        create_2Dimage_ID(
                            "piaar0zsag",
                            piaacmcopticaldesign.size,
                            piaacmcopticaldesign.size,
                            &IDpiaar0zsag)
                    );
                }


                FILE *fpri;
                {
                    char fname[STRINGMAXLEN_FULLFILENAME];
                    WRITE_FULLFILENAME(fname, "%s/ri_array.txt", piaacmcparams.piaacmcconfdir);
                    if((fpri = fopen(fname, "w")) == NULL)
                    {
                        FUNC_RETURN_FAILURE("Cannot create file %s", fname);
                    }
                }
                double ri0 = OpticsMaterials_n(piaacmcopticaldesign.PIAAmaterial_code,
                                               piaacmcopticaldesign.lambda); // refractive index at central lambda
                double sag2opd_coeff0 = (ri0 - 1.0) / 2.0;
                for(long k = 0; k < piaacmcopticaldesign.nblambda; k++)
                {
                    // sag to OPD coeff
                    double ri = OpticsMaterials_n(piaacmcopticaldesign.PIAAmaterial_code,
                                                  piaacmcopticaldesign.lambdaarray[k]); // refractive index
                    double sag2opd_coeff = (ri - 1.0) / 2.0;
                    fprintf(fpri, "%g %.16f %.16f %.16f %.16f\n", piaacmcopticaldesign.lambdaarray[k], ri,
                            ri0, sag2opd_coeff, sag2opd_coeff / sag2opd_coeff0);
                    //                for(ii=0; ii<size*size; ii++)
                    //                  data.image[IDpiaar0zsag].array.F[k*size*size+ii] = data.image[IDpiaam0z].array.F[ii] * sag2opd_coeff/sag2opd_coeff0; //sag2opd_coeff * data.image[IDpiaam0z].array.F[ii] / sag2opd_coeff0;
                }
                fclose(fpri);

                for(long ii = 0; ii < piaacmcopticaldesign.size * piaacmcopticaldesign.size; ii++)
                {
                    data.image[IDpiaar0zsag].array.F[ii] =
                        data.image[IDpiaam0z].array.F[ii] / sag2opd_coeff0;
                }

                {
                    char fname[STRINGMAXLEN_FULLFILENAME];
                    WRITE_FULLFILENAME(fname, "%s/piaar0zsag.fits", piaacmcparams.piaacmcconfdir);
                    if(piaacmcparams.PIAACMC_save == 1)
                    {
                        FUNC_CHECK_RETURN(save_fl_fits("piaar0zsag", fname));
                        DEBUG_TRACEPOINT("Saved piaar0zsag to %s", fname);
                    }
                }
            }

        }



        if(piaacmcparams.MAKE_PIAA1shape == 1)
        {
            imageID ID0;
            linopt_imtools_image_construct("Cmodes", "piaa1Cmodescoeff", "piaa1Cz", &ID0);

            imageID ID1;
            linopt_imtools_image_construct("Fmodes", "piaa1Fmodescoeff", "piaa1Fz", &ID1);

            imageID ID = image_ID("piaam1z");
            if(ID == -1)
            {
                create_2Dimage_ID("piaam1z", piaacmcopticaldesign.size, piaacmcopticaldesign.size, &ID);
            }
            for(long ii = 0; ii < piaacmcopticaldesign.size * piaacmcopticaldesign.size; ii++)
            {
                data.image[ID].array.F[ii] = 0.0;
            }
            uint32_t size0 = data.image[ID0].md[0].size[0];
            uint32_t size1 = data.image[ID1].md[0].size[0];
            for(uint32_t ii = 0; ii < size0; ii++)
                for(uint32_t jj = 0; jj < size0; jj++)
                {
                    data.image[ID].array.F[(jj + (piaacmcopticaldesign.size - size0) / 2)*piaacmcopticaldesign.size + (ii +
                                           (piaacmcopticaldesign.size - size0) / 2)] += data.image[ID0].array.F[jj * size0 + ii];
                }
            for(uint32_t ii = 0; ii < size1; ii++)
                for(uint32_t jj = 0; jj < size1; jj++)
                {
                    data.image[ID].array.F[(jj + (piaacmcopticaldesign.size - size1) / 2)*piaacmcopticaldesign.size + (ii +
                                           (piaacmcopticaldesign.size - size1) / 2)] += data.image[ID1].array.F[jj * size1 + ii];
                }

            if(piaacmcparams.PIAACMC_save == 1)
            {
                char fname[STRINGMAXLEN_FULLFILENAME];

                WRITE_FULLFILENAME(fname, "%s/piaa1Cz.fits", piaacmcparams.piaacmcconfdir);
                FUNC_CHECK_RETURN(save_fits("piaa1Cz", fname));

                WRITE_FULLFILENAME(fname, "%s/piaa1Fz.fits", piaacmcparams.piaacmcconfdir);
                FUNC_CHECK_RETURN(save_fits("piaa1Fz", fname));

                WRITE_FULLFILENAME(fname, "%s/piaam1z.fits", piaacmcparams.piaacmcconfdir);
                FUNC_CHECK_RETURN(save_fits("piaam1z", fname));
            }
            FUNC_CHECK_RETURN(delete_image_ID("piaa1Cz", DELETE_IMAGE_ERRMODE_IGNORE));
            FUNC_CHECK_RETURN(delete_image_ID("piaa1Fz", DELETE_IMAGE_ERRMODE_IGNORE));

            imageID IDpiaam1z = ID;


            // make lense shapes if applicable
            if(piaacmcopticaldesign.PIAAmaterial_code != 0) // refractive PIAA
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
                    if((data.image[IDpiaar1zsag].md[0].size[0] != piaacmcopticaldesign.size)
                            || (data.image[IDpiaar1zsag].md[0].size[1] !=
                                piaacmcopticaldesign.size)) //||(data.image[IDpiaar1zsag].md[0].size[2] != design[index].nblambda))
                    {
                        FUNC_CHECK_RETURN(delete_image_ID("piaar1zsag", DELETE_IMAGE_ERRMODE_IGNORE));
                        mkpiaar1zsag = 1;
                    }
                }
                if(mkpiaar1zsag == 1)
                {
                    FUNC_CHECK_RETURN(
                        create_2Dimage_ID(
                            "piaar1zsag",
                            piaacmcopticaldesign.size,
                            piaacmcopticaldesign.size,
                            &IDpiaar1zsag)
                    );
                }

                FILE *fpri;
                {
                    char fname[STRINGMAXLEN_FULLFILENAME];
                    WRITE_FULLFILENAME(fname, "%s/ri_array.txt", piaacmcparams.piaacmcconfdir);
                    if((fpri = fopen(fname, "w")) == NULL)
                    {
                        PRINT_ERROR("cannot open file \"%s\"", fname);
                        abort();
                    }
                }
                double ri0 = OpticsMaterials_n(piaacmcopticaldesign.PIAAmaterial_code,
                                               piaacmcopticaldesign.lambda); // refractive index at central lambda
                double sag2opd_coeff0 = (ri0 - 1.0) / 2.0;
                for(long k = 0; k < piaacmcopticaldesign.nblambda; k++)
                {
                    // sag to OPD coeff
                    double ri = OpticsMaterials_n(piaacmcopticaldesign.PIAAmaterial_code,
                                                  piaacmcopticaldesign.lambdaarray[k]); // refractive index
                    double sag2opd_coeff = (ri - 1.0) / 2.0;
                    fprintf(fpri, "%g %.16f %.16f %.16f %.16f\n", piaacmcopticaldesign.lambdaarray[k], ri,
                            ri0, sag2opd_coeff, sag2opd_coeff / sag2opd_coeff0);
                    // for(ii=0; ii<size*size; ii++)
                    //    data.image[IDpiaar1zsag].array.F[k*size*size+ii] = sag2opd_coeff * data.image[IDpiaam1z].array.F[ii] / sag2opd_coeff0;
                }
                fclose(fpri);

                for(long ii = 0; ii < piaacmcopticaldesign.size * piaacmcopticaldesign.size; ii++)
                {
                    data.image[IDpiaar1zsag].array.F[ii] =
                        data.image[IDpiaam1z].array.F[ii] / sag2opd_coeff0;
                }

                {
                    char fname[STRINGMAXLEN_FULLFILENAME];
                    WRITE_FULLFILENAME(fname, "%s/piaar1zsag.fits", piaacmcparams.piaacmcconfdir);
                    if(piaacmcparams.PIAACMC_save == 1)
                    {
                        FUNC_CHECK_RETURN(save_fl_fits("piaar1zsag", fname));
                        DEBUG_TRACEPOINT("Saved piaar1zsag to %s", fname);
                    }
                }
            }
        }
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}

