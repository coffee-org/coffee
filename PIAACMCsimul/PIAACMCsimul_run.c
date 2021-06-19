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
#include "PIAACMCsimul/PIAACMCsimul.h"







extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTSYST *optsyst;

extern OPTPIAACMCDESIGN *piaacmc;


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
    errno_t rval = 0;

    INSERT_STD_PROCINFO_COMPUTEFUNC_START

    rval = PIAACMCsimul_run(confindex,*mode);

    INSERT_STD_PROCINFO_COMPUTEFUNC_END

    return rval;
}




INSERT_STD_FPSCLIfunctions

// Register function in CLI
errno_t CLIADDCMD_PIAACMCsimul__run()
{
    INSERT_STD_CLIREGISTERFUNC
    return RETURN_SUCCESS;
}
















/*
    entry point for PIAACMCsimul from the cli
*/
errno_t PIAACMCsimul_run(
    const char *confindex,		/// @param[in] confindex	configuration index (sets name of directory for results)
    long mode					/// @param[in] mode			operation to be executed
)
{
    long i;
    FILE *fp;
    char fname[1500];
    char fname1[1500];
    char fnamebestval[1500];
    double bestval = 1.0;
    long k;
    int fOK = 0;
    int bOK = 0;
    int zeroST = 0;
//    int r;
    int loopOK;
    char stopfile[500];
    long IDv;
    long IDbestsoltmp, IDbestsol;
    char fnamebestsol[1500];
    int loopin = 0;
    struct timeval start, end;
    long secs_used, micros_used;

    double prob1;
    double sag0; // sag to change OPD by 1 wave at central wavelength
    long cnt00, cnt0;

    double searchtime = 3600.0 * 10.0; // [second] default 10 hours

#ifdef PIAASIMUL_LOGFUNC0
    PIAACMCsimul_logFunctionCall("PIAACMCsimul.fcall.log", __FUNCTION__, __LINE__,
                                 "");
#endif

    piaacmc = NULL; // set the pointer to the piaacmc structure to null


    IDbestsol = -1; // data array index of current best solution


    // read various cli variables, possibly setting globals
    if((IDv = variable_ID("PIAACMC_MASKRADLD")) != -1)
    {
        piaacmcsimul_var.PIAACMC_MASKRADLD = data.variable[IDv].value.f;
    }


    if((IDv = variable_ID("PIAACMC_FPMsectors")) != -1)
    {
        piaacmcsimul_var.PIAACMC_FPMsectors = (long) data.variable[IDv].value.f + 0.01;
    }
    printf("PIAACMC_FPMsectors = %d\n", piaacmcsimul_var.PIAACMC_FPMsectors);

    if((IDv = variable_ID("SCORINGMASKTYPE")) != -1)
    {
        piaacmcsimul_var.SCORINGMASKTYPE = (long) data.variable[IDv].value.f + 0.01;
    }
    printf("SCORINGMASKTYPE = %d\n", piaacmcsimul_var.SCORINGMASKTYPE);

    if((IDv = variable_ID("PIAACMC_save")) != -1)
    {
        piaacmcsimul_var.PIAACMC_save = (long) data.variable[IDv].value.f + 0.01;
    }
    printf("PIAACMC_save = %d\n", piaacmcsimul_var.PIAACMC_save);



    if((IDv = variable_ID("PIAACMC_resolved")) != -1)
    {
        piaacmcsimul_var.computePSF_ResolvedTarget = (long)(data.variable[IDv].value.f +
                0.01);
    }
    if((IDv = variable_ID("PIAACMC_extmode")) != -1)
    {
        piaacmcsimul_var.computePSF_ResolvedTarget_mode = (long)(
                    data.variable[IDv].value.f + 0.01);
    }


    printf("mode = %ld\n", mode);
    fflush(stdout);



    double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;

    // set the result directories
    sprintf(piaacmcsimul_var.piaacmcconfdir, "%s", confindex);

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

    piaacmcsimul_var.PIAACMC_fpmtype = 0; // idealized (default)
    if((IDv = variable_ID("PIAACMC_fpmtype")) != -1)
    {
        piaacmcsimul_var.PIAACMC_fpmtype = (int)(data.variable[IDv].value.f + 0.1);
    }

    PIAACMCsimul_initpiaacmcconf(piaacmcsimul_var.PIAACMC_fpmtype, fpmradld,
                                 centobs0, centobs1, 0, 1);



    // compute the material thickness producing a lambda phase shift at the center of the spectral band
    sag0 = 1.0;


    // mode 13: optimize focal plane mask zones only, setting the sag values for each mask zone
    // This outer loop is to choose more different starting points for the exec loop
    if(mode == 13) // loop to keep looking for optimal solution
    {
        sprintf(fname, "searchtime.txt");
        fp = fopen(fname, "r");
        if(fp != NULL)
        {
            int r = fscanf(fp, "%lf\n", &searchtime);
            (void) r;
            fclose(fp);
        }


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


        gettimeofday(&start, NULL);
        i = 0;

        // while not exceed searchtime or no stop file
        // EXECUTE_SYSTEM_COMMAND("touch start.loop.ttxt");
        while((loopOK == 1) && (i < 1000000))
        {
            printf("LOOP start\n");
            fflush(stdout);

            //EXECUTE_SYSTEM_COMMAND("touch start.iter%05ld.ttxt", i);

            // read in the real searchtime nominally set by the bash script
            sprintf(fname, "searchtime.txt");
            fp = fopen(fname, "r");
            if(fp != NULL)
            {
                int r = fscanf(fp, "%lf\n", &searchtime);
                (void) r;
                fclose(fp);
            }

            printf("searchtime = %f sec\n", searchtime);
            fflush(stdout);

            //EXECUTE_SYSTEM_COMMAND("touch step00.iter%05ld.ttxt", i);

            loopin = 1; // loop has been initialized
            if((i < 1))
            {
                piaacmcsimul_var.MODampl = 0.0;    // MODampl is a global
            }
            else
            {
                piaacmcsimul_var.MODampl = 1.0e-6 * pow(ran1(),
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


            cnt00 = 0;
            cnt0 = 0;
            if((i > 1) && (ran1() > 0.5))
            {
                if((ran1() > 0.5) && (IDbestsol != -1))
                {
                    zeroST = 2; // starting point = optimal solution
                    printf("[%d] Adopting best solution as starting point\n", __LINE__);
                    fflush(stdout);
                    // copy the best solution to the current zoneID of array of sags
                    for(k = 0; k < data.image[piaacmc[0].zonezID].md[0].size[0]; k++)
                    {
                        data.image[piaacmc[0].zonezID].array.D[k] = data.image[IDbestsol].array.D[k];
                    }
                }
                else
                {
                    zeroST = 1; // starting point = 0
                    printf("[%d] Adopting zero as starting point\n", __LINE__);
                    printf("piaacmc[0].zonezID = %ld\n", piaacmc[0].zonezID);
                    list_image_ID();
                    fflush(stdout);
                    // zero out the current zoneID of array of sags
                    for(k = 0; k < data.image[piaacmc[0].zonezID].md[0].size[0]; k++)
                    {
                        data.image[piaacmc[0].zonezID].array.D[k] = 0.0;
                    }
                }
            }
            else
            {
                zeroST = 0;
            }

            //EXECUTE_SYSTEM_COMMAND("touch step01.iter%05ld.ttxt", i);

            // zeroST = 3 => starting point is best solution found so far.  Same as zeroST=2
            // this flags that it's this value 'cause it's third iteration
            if(i == 3)
            {
                if(IDbestsol != -1)
                {
                    printf("[%d] Adopting best solution as starting point\n", __LINE__);
                    fflush(stdout);

                    zeroST = 3;
                    for(k = 0; k < data.image[piaacmc[0].zonezID].md[0].size[0]; k++)
                    {
                        data.image[piaacmc[0].zonezID].array.D[k] = data.image[IDbestsol].array.D[k];
                    }
                    piaacmcsimul_var.MODampl = 0.0;
                }
                else
                {
                    zeroST = 1; // starting point = 0
                    printf("[%d] Adopting zero as starting point\n", __LINE__);
                    printf("piaacmc[0].zonezID = %ld\n", piaacmc[0].zonezID);
                    list_image_ID();
                    fflush(stdout);
                    // zero out the current zoneID of array of sags
                    for(k = 0; k < data.image[piaacmc[0].zonezID].md[0].size[0]; k++)
                    {
                        data.image[piaacmc[0].zonezID].array.D[k] = 0.0;
                    }
                }
            }

            if(i > 0)
            {
                EXECUTE_SYSTEM_COMMAND("echo \"%g  %ld\" > sag0.txt", sag0,
                                       (long) data.image[piaacmc[0].zonezID].md[0].size[0]);

                // randomly select regions that are abs()>sag0/2 and push them back toward zero
                prob1 = pow(ran1(),
                            8.0); // probability that each zone is pushed back toward zero

                for(k = 0; k < data.image[piaacmc[0].zonezID].md[0].size[0]; k++)
                {
                    if(data.image[piaacmc[0].zonezID].array.D[k] > sag0 / 2.0)
                    {
                        cnt00++;
                        if(ran1() < prob1)
                        {
                            data.image[piaacmc[0].zonezID].array.D[k] -= sag0;
                            cnt0++;
                        }
                    }
                    if(data.image[piaacmc[0].zonezID].array.D[k] < -sag0 / 2.0)
                    {
                        cnt00++;
                        if(ran1() < prob1)
                        {
                            data.image[piaacmc[0].zonezID].array.D[k] += sag0;
                            cnt0++;
                        }
                    }
                }


                printf("Write sag values to file fpsagtest.txt\n");
                fflush(stdout);

                fp = fopen("fpsagtest.txt", "w");
                fprintf(fp, "# %9.6f\n", sag0 * 1.0e6);
                fprintf(fp, "#    %5ld    %5ld    %5ld\n", cnt0, cnt00,
                        (long) data.image[piaacmc[0].zonezID].md[0].size[0]);
                for(k = 0; k < data.image[piaacmc[0].zonezID].md[0].size[0]; k++)
                {
                    fprintf(fp, "%5ld %9.6f\n", k,
                            data.image[piaacmc[0].zonezID].array.D[k] * 1.0e6);
                }
                fclose(fp);


            }

            //EXECUTE_SYSTEM_COMMAND("touch step02.iter%05ld.ttxt", i);

            // actually do the optmization
            printf("Execute optimization\n");
            fflush(stdout);
            {
                errno_t fret = PIAACMCsimul_exec(confindex, mode);
                if(fret != RETURN_SUCCESS)
                {
                    FUNC_RETURN_FAILURE("Call to function PIAACMCsimul_exec failed");
                }
            }

            bOK = 0; // initialize have better value flag for printing "best" in a nice place

            printf("%g m  -> %g rad\n", sag0,
                   (double) OpticsMaterials_pha_lambda(piaacmc[0].fpmmaterial_code, sag0,
                           piaacmc[0].lambda));
            fflush(stdout);

            sag0 = sag0 / (OpticsMaterials_pha_lambda(piaacmc[0].fpmmaterial_code, sag0,
                           piaacmc[0].lambda) / 2.0 / M_PI);
            printf("======================= sag0 = %g m  -> %g rad\n", sag0,
                   (double) OpticsMaterials_pha_lambda(piaacmc[0].fpmmaterial_code, sag0,
                           (double) piaacmc[0].lambda));





            // if there is no best _solution_, load the current solution
            if(IDbestsol == -1)
            {
                PIAACMCsimul_update_fnamedescr();
                sprintf(fnamebestsol, "%s/fpm_zonez.%s.best.fits",
                        piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.fnamedescr);
                printf("LOADING \"%s\"...\n", fnamebestsol);
                fflush(stdout);
                load_fits(fnamebestsol, "fpmbestsol", 0, &IDbestsol);
            }

            // set the name of the stopfile
            sprintf(stopfile, "%s/stoploop13.txt", piaacmcsimul_var.piaacmcconfdir);

            //EXECUTE_SYSTEM_COMMAND("touch step03.iter%05ld.ttxt", i);


            // on first iteration load the best _value_ if it exists
            if(i == 0)
            {
                PIAACMCsimul_update_fnamedescr();
                sprintf(fnamebestval, "%s/mode13.%s.bestval.txt",
                        piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.fnamedescr);

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

            printf("\n\n\n\n======= val = %g [%g]\n", piaacmcsimul_var.PIAACMCSIMUL_VAL,
                   bestval);
            fflush(stdout);

            if(piaacmcsimul_var.PIAACMCSIMUL_VAL <
                    bestval) // piaacmcsimul_var.PIAACMCSIMUL_VAL was set in PIAACMCsimul_exec()
            {
                //EXECUTE_SYSTEM_COMMAND("touch step05.iter%05ld.ttxt", i);

                // we have a better solution!
                bOK = 1;
                bestval = piaacmcsimul_var.PIAACMCSIMUL_VAL; // record it
                printf("============================================================   SAVING BEST MASK SOLUTION -> fpm_zonez.[name].best.fits\n");
                fflush(stdout);

                // if a previous best solution has not been identified with an index, set its index
                // by loading the current best solution.  This probably never happens
                if(IDbestsol == -1)
                {
                    PIAACMCsimul_update_fnamedescr();
                    sprintf(fnamebestsol, "%s/fpm_zonez.%s.best.fits",
                            piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.fnamedescr);

                    load_fits(fnamebestsol, "fpmbestsol", 0, &IDbestsol);
                }
                else // otherwise load the temporary best solution.  This is probably what always happens
                {
                    load_fits(fnamebestsol, "fpmbestsoltmp", 0, &IDbestsoltmp);
                    for(k = 0; k < data.image[IDbestsol].md[0].size[0]; k++)
                    {
                        data.image[IDbestsol].array.D[k] = data.image[IDbestsoltmp].array.D[k];
                    }
                    delete_image_ID("fpmbestsoltmp");
                }

                // fname1 is the name of the current solution, which is now the best solution
                PIAACMCsimul_update_fnamedescr();
                sprintf(fname1, "%s/fpm_zonez.%s.fits", piaacmcsimul_var.piaacmcconfdir,
                        piaacmcsimul_var.fnamedescr);

                // fnamebestsol is the name of the stored best solution, should always be the same
                // as the name in line 8599 (if(IDbestsol==-1)...)
                PIAACMCsimul_update_fnamedescr();
                sprintf(fnamebestsol, "%s/fpm_zonez.%s.best.fits",
                        piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.fnamedescr);

                // copy the current solution to the best solution
                {
                    char cmdstring[5000];
                    sprintf("cp %s %s", fname1, fnamebestsol);
                    EXECUTE_SYSTEM_COMMAND("%s", cmdstring);
                    EXECUTE_SYSTEM_COMMAND("echo \"%s\" > cmdlogtest.txt", cmdstring);
                }

                sprintf(fname1, "%s/fpm_sagmapHR.%s.fits", piaacmcsimul_var.piaacmcconfdir,
                        piaacmcsimul_var.fnamedescr);
                sprintf(fnamebestsol, "%s/fpm_sagmapHR.%s.best.fits",
                        piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.fnamedescr);
                {
                    char cmdstring[5000];
                    sprintf(cmdstring, "cp %s %s", fname1, fnamebestsol);
                    EXECUTE_SYSTEM_COMMAND("%s", cmdstring);
                    EXECUTE_SYSTEM_COMMAND("echo \"%s\" > cmdlogtest1.txt", cmdstring);
                }


                // write new best value in file
                fp = fopen(fnamebestval, "w");
                fprintf(fp, "%30g %d %04ld %02ld %03ld %5.2f %02d %d %s %02d %07.3f\n", bestval,
                        piaacmcsimul_var.PIAACMC_FPMsectors, (long)(1.0e9 * piaacmc[0].lambda + 0.1),
                        (long)(1.0 * piaacmc[0].lambdaB + 0.1), piaacmc[0].NBrings,
                        piaacmcsimul_var.PIAACMC_MASKRADLD, piaacmcsimul_var.computePSF_ResolvedTarget,
                        piaacmcsimul_var.computePSF_ResolvedTarget_mode, piaacmc[0].fpmmaterial_name,
                        piaacmc[0].nblambda, piaacmc[0].fpmsagreg_coeff);
                fclose(fp);

                // advertise the existence of new best solution via file signaling.  Currently no listeners?
                EXECUTE_SYSTEM_COMMAND("touch %s/newbestsol.txt",
                                       piaacmcsimul_var.piaacmcconfdir);

                //EXECUTE_SYSTEM_COMMAND("touch step06.iter%05ld.ttxt", i);
            }

            // Add current solution (possibly not best) to the mode13...opt.txt file
            PIAACMCsimul_update_fnamedescr();
            sprintf(fname, "%s/mode13.%s.opt.txt", piaacmcsimul_var.piaacmcconfdir,
                    piaacmcsimul_var.fnamedescr);
            printf("Add current solution (possibly not best) to file: %s\n", fname);
            fflush(stdout);

            // first time through, open mode13...opt.txt for additional writing
            // for additional writing.  possibly redundant on next line.
            if(fOK == 0)
            {
                fp = fopen(fname, "a");
                fclose(fp);
                fOK = 1;
            }

            //EXECUTE_SYSTEM_COMMAND("touch step07.iter%05ld.ttxt", i);

            // open mode13...opt.txt for adding and write current value
            fp = fopen(fname, "a");
            fprintf(fp,
                    "%10ld %20.5g   %16.5g -> %16.5g   (%16.5g) %d %5ld/%5ld [%12g %2d %12g %12g  %12g]",
                    i, piaacmcsimul_var.MODampl, piaacmcsimul_var.PIAACMCSIMUL_VALREF,
                    piaacmcsimul_var.PIAACMCSIMUL_VAL, bestval, zeroST, cnt0, cnt00,
                    piaacmcsimul_var.CnormFactor, piaacmc[0].nblambda, optsyst[0].flux[0],
                    piaacmcsimul_var.SCORINGTOTAL, piaacmcsimul_var.PIAACMCSIMUL_VAL0);
            if(bOK == 1) // mark it as best if it is
            {
                fprintf(fp, " BEST\n");
            }
            else
            {
                fprintf(fp, "\n");
            }
            fclose(fp);

            // if( piaacmcsimul_var.PIAACMCSIMUL_VAL > piaacmcsimul_var.PIAACMCSIMUL_VALREF )
            //	exit(0);

            i++; // increment iteration counter (!!)

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


            fflush(stdout);

            //EXECUTE_SYSTEM_COMMAND("touch step08.iter%05ld.ttxt", i);

            gettimeofday(&end, NULL);

            printf("start: %ld secs, %ld usecs\n", (long) start.tv_sec,
                   (long) start.tv_usec);
            printf("end: %ld secs, %ld usecs\n", (long) end.tv_sec, (long) end.tv_usec);

            secs_used = (end.tv_sec - start.tv_sec); //avoid overflow by subtracting first
            micros_used = ((secs_used * 1000000) + end.tv_usec) - (start.tv_usec);

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
            printf("End of loop\n");
            fflush(stdout);
        }



        // initialize loop.  loopin is always set to 1 above.
        if(loopin == 1)
        {
            printf("piaacmcconfdir              : %s\n", piaacmcsimul_var.piaacmcconfdir);
            fflush(stdout);

            printf("computePSF_ResolvedTarget   : %d\n",
                   piaacmcsimul_var.computePSF_ResolvedTarget);
            fflush(stdout);


            printf("computePSF_ResolvedTarget_mode   : %d\n",
                   piaacmcsimul_var.computePSF_ResolvedTarget_mode);
            fflush(stdout);

            printf("PIAACMC_FPMsectors   : %d\n", piaacmcsimul_var.PIAACMC_FPMsectors);
            fflush(stdout);

            printf("(long) (10.0*PIAACMC_MASKRADLD+0.1)   : %ld\n",
                   (long)(10.0 * piaacmcsimul_var.PIAACMC_MASKRADLD + 0.1));
            fflush(stdout);

            printf("piaacmc[0].NBrings   : %ld\n", piaacmc[0].NBrings);
            fflush(stdout);

            printf("piaacmc[0].nblambda   : %d\n", piaacmc[0].nblambda);
            fflush(stdout);

            // copy current solution to best solution ************************** why?
            PIAACMCsimul_update_fnamedescr();
            sprintf(fname1, "%s/fpm_zonez.%s.fits", piaacmcsimul_var.piaacmcconfdir,
                    piaacmcsimul_var.fnamedescr);

            PIAACMCsimul_update_fnamedescr();
            sprintf(fnamebestsol, "%s/fpm_zonez.%s.best.fits",
                    piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.fnamedescr);

            EXECUTE_SYSTEM_COMMAND("cp %s %s", fnamebestsol, fname1);
        }
    }
    else
    {
        {
            errno_t fret = PIAACMCsimul_exec(confindex, mode);
            if(fret != RETURN_SUCCESS)
            {
                FUNC_RETURN_FAILURE("Call to function PIAACMCsimul_exec failed");
            }
        }
    }


    free(piaacmc);

    return RETURN_SUCCESS;
}


