/**
 * @file    PIAAACMCsimul_exec_optimize_PIAA_shapes_fpmtransm.c
 * @brief   Optimize PIAA shapes and focal plane mask transmission
 *
 */

// System includes
#include <stdio.h>
#include <stdlib.h>

// milk includes
#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "OptSystProp/OptSystProp.h"

#include "PIAACMCsimul.h"

#include "PIAACMCsimul_computePSF.h"
#include "init_piaacmcopticaldesign.h"

#include "PIAAshape/makePIAAshapes.h"

/**
 * ---
 *
 * ## Mode 40: Optimize PIAA optics shapes (and focal plane mask transmission for idealized PIAACMC)
 *
 *
 */

errno_t exec_optimize_PIAA_shapes_fpmtransm()
{
    DEBUG_TRACE_FSTART();

    imageID IDv;
    double  fpmradld = 0.95; // default
    double  centobs0 = 0.3;
    double  centobs1 = 0.2;

    imageID ID_CPAfreq;

    printf(
        "=================================== mode 040 "
        "===================================\n");
    //		piaacmcparams.FORCE_CREATE_fpmza = 1;

    if((IDv = variable_ID("PIAACMC_centobs0")) != -1)
    {
        centobs0 = data.variable[IDv].value.f;
    }
    if((IDv = variable_ID("PIAACMC_centobs1")) != -1)
    {
        centobs1 = data.variable[IDv].value.f;
    }
    if((IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }

    piaacmcparams.PIAACMC_fpmtype = 0; // idealized (default)
    if((IDv = variable_ID("PIAACMC_fpmtype")) != -1)
    {
        piaacmcparams.PIAACMC_fpmtype =
            (int)(data.variable[IDv].value.f + 0.1);
    }
    printf("PIAACMC_fpmtype = %d\n", piaacmcparams.PIAACMC_fpmtype);

    piaacmcparams.FORCE_CREATE_fpmza = 1;
    {
        uint64_t initflag = INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__READCONF;
        initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF;

        if(piaacmcparams.PIAACMC_fpmtype == 1)
        {
            initflag |= INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL;
        }
        FUNC_CHECK_RETURN(init_piaacmcopticaldesign(fpmradld,
                          centobs0,
                          centobs1,
                          initflag,
                          NULL));
    }
    //       printf("data.image[piaacmcopticaldesign.zoneaID].array.D[0] = %lf\n", data.image[piaacmcopticaldesign.zoneaID].array.D[0]);
    //        sleep(10);

    FUNC_CHECK_RETURN(makePIAAshapes());
    piaacmcopticalsystem.FOCMASKarray[0].mode = 1; // use 1-fpm

    if(0)  // TEST
    {
        double valref = 0.0;

        FUNC_CHECK_RETURN(PIAACMCsimul_computePSF(0.0,
                          0.0,
                          0,
                          piaacmcopticalsystem.NBelem,
                          1,
                          0,
                          0,
                          1,
                          &valref));
        printf("valref = %g\n", valref);
        printf("EXEC CASE 40 COMPLETED\n");
    }
    else
    {
        piaacmcparams.LINOPT = 1; // perform linear optimization

        if((IDv = variable_ID("PIAACMC_nbiter")) != -1)
        {
            piaacmcparams.linopt_NBiter =
                (long) data.variable[IDv].value.f + 0.01;
        }
        else
        {
            piaacmcparams.linopt_NBiter = 1000;
        }

        long kmaxC =
            data.image[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0];
        if((IDv = variable_ID("PIAACMC_maxoptCterm")) != -1)
        {
            kmaxC = (long) data.variable[IDv].value.f + 0.01;
        }
        if(kmaxC >
                data.image[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0])
        {
            kmaxC =
                data.image[piaacmcopticaldesign.piaa0CmodesID].md[0].size[0];
        }

        long kmaxF =
            data.image[piaacmcopticaldesign.piaa0FmodesID].md[0].size[0];
        if((IDv = variable_ID("PIAACMC_maxoptFterm")) != -1)
        {
            kmaxF = (long) data.variable[IDv].value.f + 0.01;
        }
        if(kmaxF >
                data.image[piaacmcopticaldesign.piaa0FmodesID].md[0].size[0])
        {
            kmaxF =
                data.image[piaacmcopticaldesign.piaa0FmodesID].md[0].size[0];
        }

        // PIAA shapes regularization

        piaacmcparams.linopt_REGPIAASHAPES = 0; // default
        if((IDv = variable_ID("REGPIAASHAPES")) != -1)
        {
            piaacmcparams.linopt_REGPIAASHAPES =
                (long) data.variable[IDv].value.f + 0.01;
        }

        piaacmcparams.linopt_piaa0C_regcoeff = 0.0e-7; // regularization coeff
        piaacmcparams.linopt_piaa1C_regcoeff = 0.0e-7; // regularization coeff
        if((IDv = variable_ID("REGPIAA_C_COEFF")) != -1)
        {
            piaacmcparams.linopt_piaa0C_regcoeff = data.variable[IDv].value.f;
            piaacmcparams.linopt_piaa1C_regcoeff = data.variable[IDv].value.f;
        }

        piaacmcparams.linopt_piaa0C_regcoeff_alpha =
            1.0; // regularization coeff power
        piaacmcparams.linopt_piaa1C_regcoeff_alpha =
            1.0; // regularization coeff power
        if((IDv = variable_ID("REGPIAA_C_ALPHA")) != -1)
        {
            piaacmcparams.linopt_piaa0C_regcoeff_alpha =
                data.variable[IDv].value.f;
            piaacmcparams.linopt_piaa1C_regcoeff_alpha =
                data.variable[IDv].value.f;
        }

        piaacmcparams.linopt_piaa0F_regcoeff = 0.0e-7; // regularization coeff
        piaacmcparams.linopt_piaa1F_regcoeff = 0.0e-7; // regularization coeff
        if((IDv = variable_ID("REGPIAA_F_COEFF")) != -1)
        {
            piaacmcparams.linopt_piaa0F_regcoeff = data.variable[IDv].value.f;
            piaacmcparams.linopt_piaa1F_regcoeff = data.variable[IDv].value.f;
        }

        piaacmcparams.linopt_piaa0F_regcoeff_alpha =
            1.0; // regularization coeff power
        piaacmcparams.linopt_piaa1F_regcoeff_alpha =
            1.0; // regularization coeff power
        if((IDv = variable_ID("REGPIAA_F_ALPHA")) != -1)
        {
            piaacmcparams.linopt_piaa0F_regcoeff_alpha =
                data.variable[IDv].value.f;
            piaacmcparams.linopt_piaa1F_regcoeff_alpha =
                data.variable[IDv].value.f;
        }

        if(piaacmcparams.linopt_REGPIAASHAPES == 1)
        {
            printf("loading CPA modes frequ\n");
            fflush(stdout);
            ID_CPAfreq = image_ID("cpamodesfreq");
            if(ID_CPAfreq == -1)
            {
                load_fits("cpamodesfreq.fits",
                          "cpamodesfreq",
                          LOADFITS_ERRMODE_EXIT,
                          &ID_CPAfreq);
            }
        }

        // FPM SAG regularization

        piaacmcparams.linopt_REGFPMSAG = 1; // default
        if((IDv = variable_ID("REGFPMSAG")) != -1)
        {
            piaacmcparams.linopt_REGFPMSAG =
                (long) data.variable[IDv].value.f + 0.01;
        }

        piaacmcparams.linopt_number_param = 0;

        if(piaacmcparams.PIAACMC_fpmtype == 0)  // ideal mask
        {
            piaacmcparams.linopt_paramtype[piaacmcparams.linopt_number_param] =
                _DATATYPE_DOUBLE;
            piaacmcparams.linopt_paramval[piaacmcparams.linopt_number_param] =
                &data.image[piaacmcopticaldesign.zoneaID].array.D[0];
            piaacmcparams.linopt_paramdelta[piaacmcparams.linopt_number_param] =
                1.0e-3;
            piaacmcparams
            .linopt_parammaxstep[piaacmcparams.linopt_number_param] =
                2.0e-1;
            piaacmcparams.linopt_parammin[piaacmcparams.linopt_number_param] =
                -1.0;
            piaacmcparams.linopt_parammax[piaacmcparams.linopt_number_param] =
                1.0;
            piaacmcparams.linopt_number_param++;
        }
        else // real physical mask
        {
            if(variable_ID("PIAACMC_mzOPT") != -1)  // optimize zones
            {
                for(uint32_t mz = 0;
                        mz <
                        data.image[piaacmcopticaldesign.zonezID].md[0].size[0];
                        mz++)
                {
                    piaacmcparams
                    .linopt_paramtype[piaacmcparams.linopt_number_param] =
                        _DATATYPE_DOUBLE;
                    piaacmcparams
                    .linopt_paramval[piaacmcparams.linopt_number_param] =
                        &data.image[piaacmcopticaldesign.zonezID].array.D[mz];
                    piaacmcparams
                    .linopt_paramdelta[piaacmcparams.linopt_number_param] =
                        1.0e-9;
                    piaacmcparams.linopt_parammaxstep
                    [piaacmcparams.linopt_number_param] = 5.0e-8;
                    piaacmcparams
                    .linopt_parammin[piaacmcparams.linopt_number_param] =
                        -2.0e-6;
                    piaacmcparams
                    .linopt_parammax[piaacmcparams.linopt_number_param] =
                        2.0e-6;
                    piaacmcparams.linopt_number_param++;
                }
            }
        }

        for(long k = 0; k < kmaxC; k++)
        {
            piaacmcparams.linopt_paramtype[piaacmcparams.linopt_number_param] =
                _DATATYPE_FLOAT;
            piaacmcparams.linopt_paramvalf[piaacmcparams.linopt_number_param] =
                &data.image[piaacmcopticaldesign.piaa0CmodesID].array.F[k];
            piaacmcparams.linopt_paramdelta[piaacmcparams.linopt_number_param] =
                1.0e-10;
            piaacmcparams
            .linopt_parammaxstep[piaacmcparams.linopt_number_param] =
                1.0e-7;
            piaacmcparams.linopt_parammin[piaacmcparams.linopt_number_param] =
                -1.0e-3;
            piaacmcparams.linopt_parammax[piaacmcparams.linopt_number_param] =
                1.0e-3;
            piaacmcparams.linopt_number_param++;
        }

        for(long k = 0; k < kmaxC; k++)
        {
            piaacmcparams.linopt_paramtype[piaacmcparams.linopt_number_param] =
                _DATATYPE_FLOAT;
            piaacmcparams.linopt_paramvalf[piaacmcparams.linopt_number_param] =
                &data.image[piaacmcopticaldesign.piaa1CmodesID].array.F[k];
            piaacmcparams.linopt_paramdelta[piaacmcparams.linopt_number_param] =
                1.0e-10;
            piaacmcparams
            .linopt_parammaxstep[piaacmcparams.linopt_number_param] =
                1.0e-7;
            piaacmcparams.linopt_parammin[piaacmcparams.linopt_number_param] =
                -1.0e-3;
            piaacmcparams.linopt_parammax[piaacmcparams.linopt_number_param] =
                1.0e-3;
            piaacmcparams.linopt_number_param++;
        }

        for(long k = 0; k < kmaxF; k++)
        {
            piaacmcparams.linopt_paramtype[piaacmcparams.linopt_number_param] =
                _DATATYPE_FLOAT;
            piaacmcparams.linopt_paramvalf[piaacmcparams.linopt_number_param] =
                &data.image[piaacmcopticaldesign.piaa0FmodesID].array.F[k];
            piaacmcparams.linopt_paramdelta[piaacmcparams.linopt_number_param] =
                1.0e-10;
            piaacmcparams
            .linopt_parammaxstep[piaacmcparams.linopt_number_param] =
                1.0e-7;
            piaacmcparams.linopt_parammin[piaacmcparams.linopt_number_param] =
                -1.0e-3;
            piaacmcparams.linopt_parammax[piaacmcparams.linopt_number_param] =
                1.0e-3;
            piaacmcparams.linopt_number_param++;
        }

        for(long k = 0; k < kmaxF; k++)
        {
            piaacmcparams.linopt_paramtype[piaacmcparams.linopt_number_param] =
                _DATATYPE_FLOAT;
            piaacmcparams.linopt_paramvalf[piaacmcparams.linopt_number_param] =
                &data.image[piaacmcopticaldesign.piaa1FmodesID].array.F[k];
            piaacmcparams.linopt_paramdelta[piaacmcparams.linopt_number_param] =
                1.0e-10;
            piaacmcparams
            .linopt_parammaxstep[piaacmcparams.linopt_number_param] =
                1.0e-7;
            piaacmcparams.linopt_parammin[piaacmcparams.linopt_number_param] =
                -1.0e-3;
            piaacmcparams.linopt_parammax[piaacmcparams.linopt_number_param] =
                1.0e-3;
            piaacmcparams.linopt_number_param++;
        }

        piaacmcparams.FORCE_MAKE_PIAA0shape = 1;
        piaacmcparams.FORCE_MAKE_PIAA1shape = 1;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}
