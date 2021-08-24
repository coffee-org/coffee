/**
 * @file    PIAACMCsimul_run.c
 * @brief   PIAA-type coronagraph design, run
 *
 * Can design both APLCMC and PIAACMC coronagraphs
 *
 *
 */



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>


// milk includes
#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "statistic/statistic.h"

#include "OpticsMaterials/OpticsMaterials.h"
#include "OptSystProp/OptSystProp.h"


#include "PIAACMCsimul.h"

#include "PIAACMCsimul_exec.h"
#include "init_piaacmcopticaldesign.h"
#include "PIAACMCsimul_loadsavepiaacmcconf.h"







// forward declaration
errno_t PIAACMCsimul_run(
    const char *confindex,
    long mode
);



// Local variables pointers
static char *confindex;
static uint32_t *mode;


static CLICMDARGDEF farg[] =
{
    {
        CLIARG_STR, ".confindex", "configuration index", "000",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &confindex
    },
    {
        CLIARG_LONG, ".mode", "mode", "0",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &mode
    }
};

static CLICMDDATA CLIcmddata =
{
    "run",
    "Simulate PIAACMC",
    CLICMD_FIELDS_DEFAULTS
};


// detailed help
static errno_t help_function()
{
    return RETURN_SUCCESS;
}





static errno_t compute_function()
{
    DEBUG_TRACE_FSTART();

    INSERT_STD_PROCINFO_COMPUTEFUNC_START

    FUNC_CHECK_RETURN(
        PIAACMCsimul_run(confindex,*mode)
    );

    INSERT_STD_PROCINFO_COMPUTEFUNC_END

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}




INSERT_STD_FPSCLIfunctions

// Register function in CLI
errno_t CLIADDCMD_PIAACMCsimul__run()
{
    INSERT_STD_CLIREGISTERFUNC
    return RETURN_SUCCESS;
}












static errno_t PIAACMCsimul_setparam_variables(
    const char *confindex,
    long mode
)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %s %ld", confindex, mode);

    variableID IDv = -1;

    // PIAACMC design mask radius in l/D
    if((IDv = variable_ID("PIAACMC_MASKRADLD")) != -1)
    {
        piaacmcparams.PIAACMC_MASKRADLD = data.variable[IDv].value.f;
    }

    // sectors
    if((IDv = variable_ID("PIAACMC_FPMsectors")) != -1)
    {
        piaacmcparams.PIAACMC_FPMsectors = (long) data.variable[IDv].value.f + 0.01;
    }
    printf("PIAACMC_FPMsectors = %d\n", piaacmcparams.PIAACMC_FPMsectors);


    if((IDv = variable_ID("SCORINGMASKTYPE")) != -1)
    {
        piaacmcparams.SCORINGMASKTYPE = (long) data.variable[IDv].value.f + 0.01;
    }
    printf("SCORINGMASKTYPE = %d\n", piaacmcparams.SCORINGMASKTYPE);


    if((IDv = variable_ID("PIAACMC_save")) != -1)
    {
        piaacmcparams.PIAACMC_save = (long) data.variable[IDv].value.f + 0.01;
    }
    printf("PIAACMC_save = %d\n", piaacmcparams.PIAACMC_save);



    if((IDv = variable_ID("PIAACMC_resolved")) != -1)
    {
        piaacmcparams.computePSF_ResolvedTarget = (long)(data.variable[IDv].value.f +
                0.01);
    }


    if((IDv = variable_ID("PIAACMC_extmode")) != -1)
    {
        piaacmcparams.computePSF_ResolvedTarget_mode = (long)(
                    data.variable[IDv].value.f + 0.01);
    }

    printf("mode = %ld\n", mode);
    fflush(stdout);



    double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;

    // set the result directories
    WRITE_DIRNAME(piaacmcparams.piaacmcconfdir, "%s", confindex);

    // load some more cli variables
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
        piaacmcparams.PIAACMC_fpmtype = (int)(data.variable[IDv].value.f + 0.1);
    }

    {
        uint64_t initflags = INIT_PIAACMCOPTICALDESIGN_MODE__DEFAULT;

        initflags |= INIT_PIAACMCOPTICALDESIGN_MODE__READCONF;
        initflags |= INIT_PIAACMCOPTICALDESIGN_MODE__LOADPIAACMCCONF;

        if(piaacmcparams.PIAACMC_fpmtype == 1)
        {
            initflags |= INIT_PIAACMCOPTICALDESIGN_MODE__FPMPHYSICAL;
        }

        FUNC_CHECK_RETURN(
            init_piaacmcopticaldesign(
                fpmradld,
                centobs0,
                centobs1,
                initflags,
                NULL
            )
        );
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}



static double read_searchtime()
{
    double searchtime = 3600.0 * 10.0; // [second] default 10 hours
    FILE *fp;
    char fname[STRINGMAXLEN_FILENAME];

    WRITE_FILENAME(fname,  "searchtime.txt");
    fp = fopen(fname, "r");
    if(fp != NULL)
    {
        int r = fscanf(fp, "%lf\n", &searchtime);
        (void) r;
        fclose(fp);
    }

    printf("searchtime = %f sec\n", searchtime);
    fflush(stdout);

    return searchtime;
}




/*
    entry point for PIAACMCsimul from the cli
*/
errno_t PIAACMCsimul_run(
    const char *confindex,		/// @param[in] confindex	configuration index (sets name of directory for results)
    long mode					/// @param[in] mode			operation to be executed
)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %s mode %ld", confindex, mode);

#ifdef PIAASIMUL_LOGFUNC0
    PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__,
                                 "");
#endif

    //piaacmc = NULL; // set the pointer to the piaacmc structure to null


    // read various cli variables, possibly setting globals
    FUNC_CHECK_RETURN(
        PIAACMCsimul_setparam_variables(confindex, mode)
    );



    // mode 13: optimize focal plane mask zones only, setting the sag values for each mask zone
    // This outer loop is to choose more different starting points for the exec loop
    if(mode == 13) // loop to keep looking for optimal solution
    {
        double bestval = 1.0;
        int loopin = 0;


        double searchtime = read_searchtime();

        int loopOK;
        if(searchtime < 0.1)
        {
            loopOK = 0;
        }
        else
        {
            loopOK = 1;
        }
        printf("loopOK = %d\n", loopOK);
        fflush(stdout);

        struct timeval start;
        gettimeofday(&start, NULL);


        // set the name of the stopfile
        char stopfile[STRINGMAXLEN_FILENAME];
        WRITE_FILENAME(stopfile,
                       "%s/stoploop13.txt",
                       piaacmcparams.piaacmcconfdir
                      );

        int loopiter = 0;
        imageID IDbestsol = -1; // data array index of current best solution

        double sag0; // sag to change OPD by 1 wave at central wavelength
        // compute the material thickness producing a lambda phase shift at the center of the spectral band
        sag0 = 1.0;

        // while not exceed searchtime or no stop file
        // EXECUTE_SYSTEM_COMMAND("touch start.loop.ttxt");
        while((loopOK == 1) && (loopiter < 1000000))
        {
            printf("LOOP start\n");
            fflush(stdout);

            //EXECUTE_SYSTEM_COMMAND("touch start.iter%05ld.ttxt", i);

            // read in the real searchtime nominally set by the bash script
            searchtime = read_searchtime();


            //EXECUTE_SYSTEM_COMMAND("touch step00.iter%05ld.ttxt", i);

            loopin = 1; // loop has been initialized
            if((loopiter < 1))
            {
                piaacmcparams.MODampl = 0.0;    // MODampl is a global
            }
            else
            {
                piaacmcparams.MODampl = 1.0e-6 * pow(ran1(),
                                                     8.0);    // pick amplitude for random optimization starting point
            }


            // compute the material thickness producing a lambda phase shift at the center of the spectral band
            //sag0 = 1.0e-6;



            // after the first iteration, half the time set zeroST = 0
            // 1/4th the time when a best solution exists set zeroST = 1
            // 1/4th the time set zeroST = 2
            // zeroST is an information-only flag that does not control activity: it reflects
            //  the settings in each conditional
            // zeroST = 0 => starting point uses previous solution
            // zeroST = 1 => starting point is a mask that has no sag: blank focal plane mask
            // zeroST = 2 => starting point is best solution found so far
            //
            long cnt00 = 0;
            long cnt0 = 0;
            int zeroST = 0; // default
            if((loopiter > 1) && (ran1() > 0.5))
            {
                if((ran1() > 0.5) && (IDbestsol != -1))
                {
                    zeroST = 2; // starting point = optimal solution
                    printf("[%d] Adopting best solution as starting point\n", __LINE__);
                    fflush(stdout);
                    // copy the best solution to the current zoneID of array of sags
                    for(uint32_t k = 0; k < data.image[piaacmcopticaldesign.zonezID].md[0].size[0]; k++)
                    {
                        data.image[piaacmcopticaldesign.zonezID].array.D[k] = data.image[IDbestsol].array.D[k];
                    }
                }
                else
                {
                    zeroST = 1; // starting point = 0
                    printf("[%d] Adopting zero as starting point\n", __LINE__);
                    printf("piaacmcopticaldesign.zonezID = %ld\n", piaacmcopticaldesign.zonezID);
                    fflush(stdout);
                    // zero out the current zoneID of array of sags
                    for(uint32_t k = 0; k < data.image[piaacmcopticaldesign.zonezID].md[0].size[0]; k++)
                    {
                        data.image[piaacmcopticaldesign.zonezID].array.D[k] = 0.0;
                    }
                }
            }

            //EXECUTE_SYSTEM_COMMAND("touch step01.iter%05ld.ttxt", i);

            // zeroST = 3 => starting point is best solution found so far.  Same as zeroST=2
            // this flags that it's this value 'cause it's third iteration
            if(loopiter == 3)
            {
                if(IDbestsol != -1)
                {
                    printf("[%d] Adopting best solution as starting point\n", __LINE__);
                    fflush(stdout);

                    zeroST = 3;
                    for(uint32_t k = 0; k < data.image[piaacmcopticaldesign.zonezID].md[0].size[0]; k++)
                    {
                        data.image[piaacmcopticaldesign.zonezID].array.D[k] = data.image[IDbestsol].array.D[k];
                    }
                    piaacmcparams.MODampl = 0.0;
                }
                else
                {
                    zeroST = 1; // starting point = 0
                    printf("[%d] Adopting zero as starting point\n", __LINE__);
                    printf("piaacmcopticaldesign.zonezID = %ld\n", piaacmcopticaldesign.zonezID);
                    fflush(stdout);
                    // zero out the current zoneID of array of sags
                    for(uint32_t k = 0; k < data.image[piaacmcopticaldesign.zonezID].md[0].size[0]; k++)
                    {
                        data.image[piaacmcopticaldesign.zonezID].array.D[k] = 0.0;
                    }
                }
            }

            if(loopiter > 0)
            {
                EXECUTE_SYSTEM_COMMAND("echo \"%g  %ld\" > sag0.txt", sag0,
                                       (long) data.image[piaacmcopticaldesign.zonezID].md[0].size[0]);

                // randomly select regions that are abs()>sag0/2 and push them back toward zero
                // probability that each zone is pushed back toward zero
                double prob1 = pow(ran1(), 8.0);

                for(uint32_t k = 0; k < data.image[piaacmcopticaldesign.zonezID].md[0].size[0]; k++)
                {
                    if(data.image[piaacmcopticaldesign.zonezID].array.D[k] > sag0 / 2.0)
                    {
                        cnt00++;
                        if(ran1() < prob1)
                        {
                            data.image[piaacmcopticaldesign.zonezID].array.D[k] -= sag0;
                            cnt0++;
                        }
                    }
                    if(data.image[piaacmcopticaldesign.zonezID].array.D[k] < -sag0 / 2.0)
                    {
                        cnt00++;
                        if(ran1() < prob1)
                        {
                            data.image[piaacmcopticaldesign.zonezID].array.D[k] += sag0;
                            cnt0++;
                        }
                    }
                }


                printf("Write sag values to file fpsagtest.txt\n");
                fflush(stdout);

                {
                    FILE *fp;

                    fp = fopen("fpsagtest.txt", "w");
                    fprintf(fp, "# %9.6f\n", sag0 * 1.0e6);
                    fprintf(fp, "#    %5ld    %5ld    %5ld\n", cnt0, cnt00,
                            (long) data.image[piaacmcopticaldesign.zonezID].md[0].size[0]);
                    for(uint32_t k = 0; k < data.image[piaacmcopticaldesign.zonezID].md[0].size[0]; k++)
                    {
                        fprintf(fp,
                                "%5ld %9.6f\n",
                                (long) k,
                                data.image[piaacmcopticaldesign.zonezID].array.D[k] * 1.0e6);
                    }
                    fclose(fp);
                }


            }

            //EXECUTE_SYSTEM_COMMAND("touch step02.iter%05ld.ttxt", i);

            // Perform the optmization
            printf("Execute optimization\n");
            {
                FUNC_CHECK_RETURN(PIAACMCsimul_exec(confindex, mode));
            }


            printf("%g m  -> %g rad\n",
                   sag0,
                   (double) OpticsMaterials_pha_lambda(piaacmcopticaldesign.fpmmaterial_code, sag0, piaacmcopticaldesign.lambda)
                  );
            fflush(stdout);

            sag0 = sag0 / (OpticsMaterials_pha_lambda(piaacmcopticaldesign.fpmmaterial_code, sag0,
                           piaacmcopticaldesign.lambda) / 2.0 / M_PI);

            printf("======================= sag0 = %g m  -> %g rad\n",
                   sag0,
                   (double) OpticsMaterials_pha_lambda(piaacmcopticaldesign.fpmmaterial_code, sag0, (double) piaacmcopticaldesign.lambda)
                  );





            // if there is no best _solution_, load the current solution
            if(IDbestsol == -1)
            {
                char fnamebestsol[STRINGMAXLEN_FILENAME];

                PIAACMCsimul_update_fnamedescr();

                WRITE_FILENAME(fnamebestsol,
                               "%s/fpm_zonez.%s.best.fits",
                               piaacmcparams.piaacmcconfdir,
                               piaacmcparams.fnamedescr
                              );

                printf("LOADING \"%s\"...\n", fnamebestsol);
                fflush(stdout);

                FUNC_CHECK_RETURN(
                    load_fits(fnamebestsol,
                              "fpmbestsol",
                              LOADFITS_ERRMODE_IGNORE,
                              &IDbestsol
                             )
                );
            }


            //EXECUTE_SYSTEM_COMMAND("touch step03.iter%05ld.ttxt", i);


            PIAACMCsimul_update_fnamedescr();

            char fnamebestval[STRINGMAXLEN_FILENAME];
            WRITE_FILENAME(fnamebestval,
                           "%s/mode13.%s.bestval.txt",
                           piaacmcparams.piaacmcconfdir,
                           piaacmcparams.fnamedescr
                          );


            // on first iteration load the best _value_ if it exists
            if(loopiter == 0)
            {
                FILE *fp;
                printf("READING FILE \"%s\"\n", fnamebestval);
                fflush(stdout);

                fp = fopen(fnamebestval, "r");
                if(fp != NULL)
                {
                    int r = fscanf(fp, "%lf",
                                   &bestval); // this loads only the first value on the line
                    (void) r;
                    fclose(fp);
                }
            }

            //EXECUTE_SYSTEM_COMMAND("touch step04.iter%05ld.ttxt", i);

            printf("\n\n\n\n======= val = %g [%g]\n",
                   piaacmcparams.PIAACMCSIMUL_VAL,
                   bestval);
            fflush(stdout);


            // Check if the best solution should be updated
            // piaacmcparams.PIAACMCSIMUL_VAL was set in PIAACMCsimul_exec()
            int bOK = 0; // initialize have better value flag for printing "best" in a nice place
            if(piaacmcparams.PIAACMCSIMUL_VAL < bestval)
            {
                //EXECUTE_SYSTEM_COMMAND("touch step05.iter%05ld.ttxt", i);

                // we have a better solution!
                bOK = 1;
                bestval = piaacmcparams.PIAACMCSIMUL_VAL; // record it
                printf("========================================= \
                       SAVING BEST MASK SOLUTION -> fpm_zonez.[name].best.fits\n");



                {
                    PIAACMCsimul_update_fnamedescr();

                    char fnamebestsol[STRINGMAXLEN_FILENAME];
                    WRITE_FILENAME(fnamebestsol,
                                   "%s/fpm_zonez.%s.best.fits",
                                   piaacmcparams.piaacmcconfdir,
                                   piaacmcparams.fnamedescr
                                  );

                    // if a previous best solution has not been identified with an index, set its index
                    // by loading the current best solution.  This probably never happens
                    if(IDbestsol == -1)
                    {
                        load_fits(fnamebestsol,
                                  "fpmbestsol",
                                  LOADFITS_ERRMODE_IGNORE,
                                  &IDbestsol
                                 );
                    }
                    else // otherwise load the temporary best solution.  This is probably what always happens
                    {
                        imageID IDbestsoltmp;
                        load_fits(fnamebestsol,
                                  "fpmbestsoltmp",
                                  LOADFITS_ERRMODE_IGNORE,
                                  &IDbestsoltmp
                                 );

                        for(uint32_t k = 0; k < data.image[IDbestsol].md[0].size[0]; k++)
                        {
                            data.image[IDbestsol].array.D[k] = data.image[IDbestsoltmp].array.D[k];
                        }

                        delete_image_ID("fpmbestsoltmp",
                                        DELETE_IMAGE_ERRMODE_WARNING);
                    }
                }




                {
                    // fname1 is the name of the current solution, which is now the best solution
                    PIAACMCsimul_update_fnamedescr();

                    char fname1[STRINGMAXLEN_FILENAME];
                    WRITE_FILENAME(fname1,
                                   "%s/fpm_zonez.%s.fits",
                                   piaacmcparams.piaacmcconfdir,
                                   piaacmcparams.fnamedescr
                                  );


                    // fnamebestsol is the name of the stored best solution, should always be the same
                    // as the name in line 8599 (if(IDbestsol==-1)...)
                    PIAACMCsimul_update_fnamedescr();

                    char fnamebestsol[STRINGMAXLEN_FILENAME];
                    WRITE_FILENAME(fnamebestsol,
                                   "%s/fpm_zonez.%s.best.fits",
                                   piaacmcparams.piaacmcconfdir,
                                   piaacmcparams.fnamedescr
                                  );

                    // copy the current solution to the best solution
                    EXECUTE_SYSTEM_COMMAND("cp %s %s", fname1, fnamebestsol);
                    EXECUTE_SYSTEM_COMMAND(
                        "echo \"cp %s %s\" > cmdlogtest.txt",
                        fname1,
                        fnamebestsol
                    );


                    WRITE_FILENAME(fname1,
                                   "%s/fpm_sagmapHR.%s.fits",
                                   piaacmcparams.piaacmcconfdir,
                                   piaacmcparams.fnamedescr
                                  );

                    WRITE_FILENAME(fnamebestsol,
                                   "%s/fpm_sagmapHR.%s.best.fits",
                                   piaacmcparams.piaacmcconfdir,
                                   piaacmcparams.fnamedescr
                                  );

                    EXECUTE_SYSTEM_COMMAND("cp %s %s", fname1, fnamebestsol);
                    EXECUTE_SYSTEM_COMMAND(
                        "echo \"cp %s %s\" > cmdlogtest1.txt",
                        fname1,
                        fnamebestsol
                    );
                }

                {
                    // write new best value in file
                    FILE *fp;

                    fp = fopen(fnamebestval, "w");
                    fprintf(fp,
                            "%30g %d %04ld %02ld %03ld %5.2f %02d %d %s %02d %07.3f\n",
                            bestval,
                            piaacmcparams.PIAACMC_FPMsectors,
                            (long)(1.0e9 * piaacmcopticaldesign.lambda + 0.1),
                            (long)(1.0 * piaacmcopticaldesign.lambdaB + 0.1),
                            piaacmcopticaldesign.NBrings,
                            piaacmcparams.PIAACMC_MASKRADLD,
                            piaacmcparams.computePSF_ResolvedTarget,
                            piaacmcparams.computePSF_ResolvedTarget_mode,
                            piaacmcopticaldesign.fpmmaterial_name,
                            piaacmcopticaldesign.nblambda,
                            piaacmcopticaldesign.fpmsagreg_coeff
                           );
                    fclose(fp);
                }

                // advertise the existence of new best solution via file signaling.  Currently no listeners?
                EXECUTE_SYSTEM_COMMAND("touch %s/newbestsol.txt",
                                       piaacmcparams.piaacmcconfdir);

                //EXECUTE_SYSTEM_COMMAND("touch step06.iter%05ld.ttxt", i);
            }


            // Add current solution (possibly not best) to the mode13...opt.txt file
            if (1)
            {
                FILE *fp;

                PIAACMCsimul_update_fnamedescr();

                char fname[STRINGMAXLEN_FILENAME];
                WRITE_FILENAME(fname,
                               "%s/mode13.%s.opt.txt",
                               piaacmcparams.piaacmcconfdir,
                               piaacmcparams.fnamedescr
                              );

                printf("Add current solution (possibly not best) to file: %s\n", fname);

                //EXECUTE_SYSTEM_COMMAND("touch step07.iter%05ld.ttxt", i);

                // open mode13...opt.txt for adding and write current value
                fp = fopen(fname, "a");
                fprintf(fp,
                        "%10d %20.5g   %16.5g -> %16.5g   (%16.5g) %d %5ld/%5ld [%12g %2d %12g %12g  %12g]",
                        loopiter,
                        piaacmcparams.MODampl,
                        piaacmcparams.PIAACMCSIMUL_VALREF,
                        piaacmcparams.PIAACMCSIMUL_VAL,
                        bestval,
                        zeroST,
                        cnt0,
                        cnt00,
                        piaacmcparams.CnormFactor,
                        piaacmcopticaldesign.nblambda,
                        piaacmcopticalsystem.flux[0],
                        piaacmcparams.SCORINGTOTAL,
                        piaacmcparams.PIAACMCSIMUL_VAL0
                       );
                if(bOK == 1) // mark it as best if it is
                {
                    fprintf(fp, " BEST\n");
                }
                else
                {
                    fprintf(fp, "\n");
                }
                fclose(fp);
            }


            loopiter++; // increment iteration counter

            // stop iterations if stopfile exists
            if(file_exists(stopfile) == 1)
            {
                printf("FILE \"%s\" found -> stopping\n", stopfile);
                loopOK = 0;
                EXECUTE_SYSTEM_COMMAND("touch stop.stopfile.txt");
                EXECUTE_SYSTEM_COMMAND("rm %s", stopfile);
            }
            else
            {
                printf("File \"%s\" not found -> continuing\n", stopfile);
            }


            // update timing : write used time in timeused.txt
            // and check to see if time has run out
            if (1)
            {
                //EXECUTE_SYSTEM_COMMAND("touch step08.iter%05ld.ttxt", i);
                long secs_used, micros_used;

                struct timeval end;
                gettimeofday(&end, NULL);

                printf("start: %ld secs, %ld usecs\n", (long) start.tv_sec,
                       (long) start.tv_usec);
                printf("end: %ld secs, %ld usecs\n", (long) end.tv_sec, (long) end.tv_usec);

                secs_used = (end.tv_sec - start.tv_sec); //avoid overflow by subtracting first
                micros_used = ((secs_used * 1000000) + end.tv_usec) - (start.tv_usec);

                FILE *fp;
                fp = fopen("timeused.txt", "w");
                fprintf(fp, "# Time used vs. time limit\n");
                fprintf(fp, "# Time limit can be written in file searchtime.txt\n");
                fprintf(fp, "\n");
                fprintf(fp, "%12.3f    %12.3f\n", 1.0e-6 * micros_used, searchtime);
                fclose(fp);


                //EXECUTE_SYSTEM_COMMAND("touch step09.iter%05ld.ttxt", i);

                // check to see if time has run out
                if(micros_used > 1000000.0 * searchtime) // searchtime is in seconds
                {
                    loopOK = 0; // stop loop flag
                    EXECUTE_SYSTEM_COMMAND("touch stop.time_elapsed.ttxt");
                }
                //EXECUTE_SYSTEM_COMMAND("touch step10.iter%05ld.ttxt", i);
            }

            printf("End of loop\n");
        }



        // initialize loop.  loopin is always set to 1 above.
        if(loopin == 1)
        {
            printf("piaacmcconfdir              : %s\n", piaacmcparams.piaacmcconfdir);
            printf("computePSF_ResolvedTarget   : %d\n",
                   piaacmcparams.computePSF_ResolvedTarget);
            printf("computePSF_ResolvedTarget_mode   : %d\n",
                   piaacmcparams.computePSF_ResolvedTarget_mode);
            printf("PIAACMC_FPMsectors   : %d\n", piaacmcparams.PIAACMC_FPMsectors);
            printf("(long) (10.0*PIAACMC_MASKRADLD+0.1)   : %ld\n",
                   (long)(10.0 * piaacmcparams.PIAACMC_MASKRADLD + 0.1));
            printf("piaacmcopticaldesign.NBrings   : %ld\n", piaacmcopticaldesign.NBrings);
            printf("piaacmcopticaldesign.nblambda   : %d\n", piaacmcopticaldesign.nblambda);
            fflush(stdout);


            {
                // copy current solution to best solution ************************** why?
                char fname[STRINGMAXLEN_FILENAME];
                char fnamebestsol[STRINGMAXLEN_FILENAME];

                PIAACMCsimul_update_fnamedescr();

                WRITE_FILENAME(fname,
                               "%s/fpm_zonez.%s.fits",
                               piaacmcparams.piaacmcconfdir,
                               piaacmcparams.fnamedescr
                              );

                PIAACMCsimul_update_fnamedescr();

                WRITE_FILENAME(fnamebestsol,
                               "%s/fpm_zonez.%s.best.fits",
                               piaacmcparams.piaacmcconfdir,
                               piaacmcparams.fnamedescr
                              );

                EXECUTE_SYSTEM_COMMAND("cp %s %s", fnamebestsol, fname);
            }
        }
    }
    else
    {
        {
            DEBUG_TRACEPOINT("running exec function");
            FUNC_CHECK_RETURN(
                PIAACMCsimul_exec(confindex, mode)
            );
        }
    }


    //free(piaacmcoptdesign);

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}


