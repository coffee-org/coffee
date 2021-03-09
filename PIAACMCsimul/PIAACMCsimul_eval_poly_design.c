/**
 * @file    PIAAACMCsimul_eval_poly_design.c
 * @brief   Evaluate polychromatic design
 * 
 *  
 * @author  O. Guyon
 * @date    25 nov 2017
 *
 * 
 * @bug No known bugs.
 * 
 */



// System includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// milk includes
#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_tools/COREMOD_tools.h"

#include "OptSystProp/OptSystProp.h"
#include "PIAACMCsimul/PIAACMCsimul.h"



extern PIAACMCsimul_varType piaacmcsimul_var;

extern OPTSYST *optsyst;

extern OPTPIAACMCDESIGN *piaacmc;





int PIAACMCsimul_eval_poly_design()
{
    long IDopderrC;
    long nbOPDerr;
    long IDv;
    double fpmradld = 0.95;  // default
    double centobs0 = 0.3;
    double centobs1 = 0.2;
    double ldoffset;
    double valref;

    long ID;
    char fname[1500];
    long xsize, ysize, zsize;

    float *peakarray;
    float avpeak;

    long kk, ii, jj;
    double val;

    double focscale;
    double xc, yc, rc;

    double *eval_contrastCurve;
    double *eval_COHcontrastCurve;
    double *eval_INCcontrastCurve;
    double *eval_contrastCurve_cnt;
    double eval_sepstepld = 0.2; // in l/D
    double eval_sepmaxld = 20.0; // in l/D
    long eval_sepNBpt;

    long ri;
    FILE *fp;
    long IDps, IDps_re, IDps_im, IDps_COH, IDps_INC, IDre, IDim, IDopderr, IDrc;

    float *rcarray;
    float drc;
    long rcbin;
    long rcbinmax;
    float rcmin, rcmax;
    long rc_cnt;
    float p50, p20;

    long NBpt;

    double aveC;
    double aveC_COH;
    double aveC_INC;
    long aveCcnt;

    long OPDmode;
    long size;

    printf("=================================== mode 100 ===================================\n");

    if( (IDv = variable_ID("PIAACMC_centobs0")) != -1)
        centobs0 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_centobs1")) != -1)
        centobs1 = data.variable[IDv].value.f;
    if( (IDv = variable_ID("PIAACMC_fpmradld")) != -1)
    {
        fpmradld = data.variable[IDv].value.f;
        printf("MASK RADIUS = %lf lambda/D\n", fpmradld);
    }



    // measure sensitivity to errors
    //	printf("Loading (optional) OPDerr file\n");
    //	fflush(stdout);
    // load an error if it exists
    IDopderrC = image_ID("OPDerrC");
    if(IDopderrC == -1)
        IDopderrC = load_fits("OPDerrC.fits", "OPDerrC", 0);


    if(IDopderrC != -1)
    {
        nbOPDerr = data.image[IDopderrC].md[0].size[2];  // number of error arrays
        //printf("INCLUDING %ld OPD ERROR MODES\n", nbOPDerr);
        //fflush(stdout);
    }
    else
    {
        //printf("NO OPD ERROR MODES\n");
        //fflush(stdout);
        nbOPDerr = 0;
    }



    printf("Will add optional OPD error modes (%ld modes)\n", nbOPDerr);
    fflush(stdout);



    piaacmcsimul_var.PIAACMC_fpmtype = 0; // idealized (default)
    if((IDv=variable_ID("PIAACMC_fpmtype"))!=-1)
        piaacmcsimul_var.PIAACMC_fpmtype = (int) (data.variable[IDv].value.f + 0.1);




    piaacmcsimul_var.FORCE_CREATE_fpmza = 1;
    
    // debug from Justin
    //printf("Initializing PIAACMC configuration\n");
    //fflush(stdout);
    //sleep(10);
    PIAACMCsimul_initpiaacmcconf(piaacmcsimul_var.PIAACMC_fpmtype, fpmradld, centobs0, centobs1, 0, 1);
    
    


    PIAACMCsimul_makePIAAshapes(piaacmc, 0);
    optsyst[0].FOCMASKarray[0].mode = 1; // use 1-fpm

    //        valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 1, 0);
    //       printf("valref = %g\n", valref);

    //exit(0);




    ldoffset = 0.01; // default
    if((IDv=variable_ID("PIAACMC_ldoffset"))!=-1)
        ldoffset = data.variable[IDv].value.f;

    printf("ldoffset = %f\n", ldoffset);

    // compute off-axis POINT source
    // debug from Justin
    printf("Computing off-axis point spread function\n");
    fflush(stdout);
    //sleep(10);
    valref = PIAACMCsimul_computePSF(5.0, 0.0, 0, optsyst[0].NBelem, 1, 0, 0, 0);
    sprintf(fname,"!%s/psfi0_x50_y00.fits", piaacmcsimul_var.piaacmcconfdir);
    save_fits("psfi0", fname);
    //load_fits(fname, "psfi");




    ID = image_ID("psfi0");
    xsize = data.image[ID].md[0].size[0];
    ysize = data.image[ID].md[0].size[1];
    zsize = data.image[ID].md[0].size[2];
    peakarray = (float*) malloc(sizeof(float)*zsize);
    for(kk=0; kk<zsize; kk++)
    {
        peakarray[kk] = 0.0;
        for(ii=0; ii<xsize*ysize; ii++)
        {
            val = data.image[ID].array.F[kk*xsize*ysize+ii];
            if(val>peakarray[kk])
                peakarray[kk] = val;
        }
    }
    avpeak = 0.0;
    for(kk=0; kk<zsize; kk++)
    {
        printf("peak %02ld  %10lf\n", kk, peakarray[kk]);
        avpeak += peakarray[kk];
    }
    avpeak /= zsize;
    free(peakarray);
    delete_image_ID("psfi0");


    // compute on-axis POINT source
    // debug from Justin
    printf("Computing on-axis point spread function\n");
    fflush(stdout);
    //sleep(10);
    valref = PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 1, 0, 0, 1);
    sprintf(fname,"!%s/psfi0_x00_y00.fits", piaacmcsimul_var.piaacmcconfdir);
    save_fits("psfi0", fname);
    //load_fits(fname, "psfi");


    ID = image_ID("psfi0");




    /// compute contrast curve
    focscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size;
    printf("focscale = %f\n", focscale);
    eval_sepNBpt = (long) (eval_sepmaxld/eval_sepstepld);
    eval_contrastCurve = (double*) malloc(sizeof(double)*eval_sepNBpt);
    eval_contrastCurve_cnt = (double*) malloc(sizeof(double)*eval_sepNBpt);
    for(ri=0; ri<eval_sepNBpt; ri++)
    {
        eval_contrastCurve[ri] = 0.0;
        eval_contrastCurve_cnt[ri] = 0.0;
    }

    aveC = 0.0;
    aveCcnt = 0;

    for(kk=0; kk<zsize; kk++)
        for(ii=0; ii<xsize; ii++)
            for(jj=0; jj<ysize; jj++)
            {
                xc = 1.0*ii-0.5*xsize;
                yc = 1.0*jj-0.5*ysize;
                xc *= focscale;
                yc *= focscale;
                rc = sqrt(xc*xc+yc*yc);
                ri = (long) (rc/eval_sepstepld-0.5);
                if(ri<0)
                    ri = 0;
                if(ri<eval_sepNBpt)
                {
                    eval_contrastCurve[ri] += data.image[ID].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                    eval_contrastCurve_cnt[ri] += 1.0;
                }
                if((rc>2.0)&&(rc<6.0))
                {
                    aveC += data.image[ID].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                    aveCcnt++;
                }
            }


	PIAACMCsimul_update_fnamedescr();
    sprintf(fname, "%s/ContrastCurve_tt000.%s.txt", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.fnamedescr);
    fp = fopen(fname, "w");
    fprintf(fp, "# Contrast curve\n");
    fprintf(fp, "# col 1: Angular separation [l/D]\n");
    fprintf(fp, "# col 2: Contrast value\n");
    fprintf(fp, "# col 3: Number of points used for contrast measurement\n");
    fprintf(fp, "\n");    
    for(ri=0; ri<eval_sepNBpt; ri++)
    {
        eval_contrastCurve[ri] /= eval_contrastCurve_cnt[ri]+0.000001;
        fprintf(fp, "%10f %10g %10g\n", eval_sepstepld*ri, eval_contrastCurve[ri], eval_contrastCurve_cnt[ri]);
    }
    fclose(fp);
    
    free(eval_contrastCurve);
    free(eval_contrastCurve_cnt);
    
	PIAACMCsimul_update_fnamedescr();
    sprintf(fname, "%s/ContrastVal_tt000.%s.txt", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.fnamedescr);
    fp = fopen(fname, "w");
    fprintf(fp, "%10g %10g %4.2f %4.2f %4.2f %4.2f %04ld %02ld %d %03ld %03ld %02d %03ld %02d %d\n", valref, aveC/aveCcnt, piaacmc[0].fpmaskradld, piaacmcsimul_var.PIAACMC_MASKRADLD, piaacmc[0].centObs0, piaacmc[0].centObs1, (long) (piaacmc[0].lambda*1e9), (long) (piaacmc[0].lambdaB+0.1), piaacmcsimul_var.PIAACMC_FPMsectors, piaacmc[0].NBrings, piaacmc[0].focmNBzone, piaacmc[0].nblambda, (long) 0, piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode);
    fclose(fp);




    // measure pointing sensitivity
    IDps = create_3Dimage_ID("starim", piaacmc[0].size, piaacmc[0].size, zsize);

    IDps_re = create_3Dimage_ID("starim_re", piaacmc[0].size, piaacmc[0].size, zsize);
    IDps_im = create_3Dimage_ID("starim_im", piaacmc[0].size, piaacmc[0].size, zsize);

    IDps_COH = create_3Dimage_ID("starimCOH", piaacmc[0].size, piaacmc[0].size, zsize);
    IDps_INC = create_3Dimage_ID("starimINC", piaacmc[0].size, piaacmc[0].size, zsize);

    NBpt = 0;

    valref = 0.25*PIAACMCsimul_computePSF(ldoffset, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);
    sprintf(fname,"!%s/psfi0_p0.fits", piaacmcsimul_var.piaacmcconfdir);
    save_fits("psfi0", fname);
    ID = image_ID("psfi0");
    IDre = image_ID("psfre0");
    IDim = image_ID("psfim0");
    for(ii=0; ii<xsize*ysize*zsize; ii++) {
        data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];
        data.image[IDps_re].array.F[ii] += data.image[IDre].array.F[ii];
        data.image[IDps_im].array.F[ii] += data.image[IDim].array.F[ii];
    }
    NBpt++;
    delete_image_ID("psfre0");
    delete_image_ID("psfim0");

    valref += 0.25*PIAACMCsimul_computePSF(-ldoffset, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);
    sprintf(fname,"!%s/psfi0_m0.fits", piaacmcsimul_var.piaacmcconfdir);
    save_fits("psfi0", fname);
    ID = image_ID("psfi0");
    IDre = image_ID("psfre0");
    IDim = image_ID("psfim0");
    for(ii=0; ii<xsize*ysize*zsize; ii++) {
        data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];
        data.image[IDps_re].array.F[ii] += data.image[IDre].array.F[ii];
        data.image[IDps_im].array.F[ii] += data.image[IDim].array.F[ii];
    }
    NBpt++;
    delete_image_ID("psfre0");
    delete_image_ID("psfim0");

    valref += 0.25*PIAACMCsimul_computePSF(0.0, ldoffset, 0, optsyst[0].NBelem, 0, 0, 0, 0);
    sprintf(fname,"!%s/psfi0_0p.fits", piaacmcsimul_var.piaacmcconfdir);
    save_fits("psfi0", fname);
    ID = image_ID("psfi0");
    IDre = image_ID("psfre0");
    IDim = image_ID("psfim0");
    for(ii=0; ii<xsize*ysize*zsize; ii++) {
        data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];
        data.image[IDps_re].array.F[ii] += data.image[IDre].array.F[ii];
        data.image[IDps_im].array.F[ii] += data.image[IDim].array.F[ii];
    }
    NBpt++;
    delete_image_ID("psfre0");
    delete_image_ID("psfim0");

    valref += 0.25*PIAACMCsimul_computePSF(0.0, -ldoffset, 0, optsyst[0].NBelem, 0, 0, 0, 0);
    sprintf(fname,"!%s/psfi0_0m.fits", piaacmcsimul_var.piaacmcconfdir);
    save_fits("psfi0", fname);
    ID = image_ID("psfi0");
    IDre = image_ID("psfre0");
    IDim = image_ID("psfim0");
    for(ii=0; ii<xsize*ysize*zsize; ii++) {
        data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];
        data.image[IDps_re].array.F[ii] += data.image[IDre].array.F[ii];
        data.image[IDps_im].array.F[ii] += data.image[IDim].array.F[ii];
    }
    NBpt++;
    delete_image_ID("psfre0");
    delete_image_ID("psfim0");




    // measure sensitivity to errors

    //	printf("Adding optional OPD error modes (%ld modes)\n", nbOPDerr);
    //	fflush(stdout);
    // add error modes if any
    for(OPDmode=0; OPDmode < nbOPDerr; OPDmode++)
    {
        size = data.image[IDopderrC].md[0].size[0];
        IDopderr = create_2Dimage_ID("opderr", size, size);
        // "opderr" is a standard name read by PIAACMCsimul_init
        for(ii=0; ii<size*size; ii++)
            data.image[IDopderr].array.F[ii] = data.image[IDopderrC].array.F[size*size*OPDmode + ii];
        PIAACMCsimul_init(piaacmc, 0, 0.0, 0.0); // add error to the data
        PIAACMCsimul_computePSF(0.0, 0.0, 0, optsyst[0].NBelem, 0, 0, 0, 0);
        sprintf(fname, "!%s/psfi0_opderr%02ld.fits", piaacmcsimul_var.piaacmcconfdir, OPDmode);
        save_fits("psfi0", fname);
        delete_image_ID("opderr");
        ID = image_ID("psfi0");
        IDre = image_ID("psfre0");
        IDim = image_ID("psfim0");
        for(ii=0; ii<xsize*ysize*zsize; ii++) {
            data.image[IDps].array.F[ii] += data.image[ID].array.F[ii];
            data.image[IDps_re].array.F[ii] += data.image[IDre].array.F[ii];
            data.image[IDps_im].array.F[ii] += data.image[IDim].array.F[ii];
        }
        NBpt++;
        delete_image_ID("psfre0");
        delete_image_ID("psfim0");
    }


    for(ii=0; ii<xsize*ysize*zsize; ii++) {
        data.image[IDps].array.F[ii] /= NBpt;
        data.image[IDps_re].array.F[ii] /= NBpt; // average re
        data.image[IDps_im].array.F[ii] /= NBpt; // average im
    }


    for(ii=0; ii<xsize*ysize*zsize; ii++) {
        data.image[IDps_COH].array.F[ii] = data.image[IDps_re].array.F[ii]*data.image[IDps_re].array.F[ii] + data.image[IDps_im].array.F[ii]*data.image[IDps_im].array.F[ii];
        data.image[IDps_INC].array.F[ii] = data.image[IDps].array.F[ii] - data.image[IDps_COH].array.F[ii];
    }


	// same as psfi0_extsrc 
	//    sprintf(fname, "!%s/psfi0_starim.%s.fits", piaacmcsimul_var.piaacmcconfdir, piaacmcsimul_var.fnamedescr);
	//    save_fits("starim", fname);



	PIAACMCsimul_update_fnamedescr();
    sprintf(fname, "!%s/psfi0_extsrc%2ld_sm%d.%s.fits", piaacmcsimul_var.piaacmcconfdir, (long) (-log10(ldoffset)*10.0+0.1), piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.fnamedescr);
    save_fits("starim", fname);

	PIAACMCsimul_update_fnamedescr();
    sprintf(fname, "!%s/psfi0INC_extsrc%2ld_sm%d.%s.fits", piaacmcsimul_var.piaacmcconfdir, (long) (-log10(ldoffset)*10.0+0.1), piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.fnamedescr);
    save_fits("starimINC", fname);

	PIAACMCsimul_update_fnamedescr();
    sprintf(fname, "!%s/psfi0COH_extsrc%2ld_sm%d.%s.fits", piaacmcsimul_var.piaacmcconfdir, (long) (-log10(ldoffset)*10.0+0.1), piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.fnamedescr);
    save_fits("starimCOH", fname);




    /// compute contrast curve
    /// measure average contrast value, 2-6 lambda/D
    focscale = (2.0*piaacmc[0].beamrad/piaacmc[0].pixscale)/piaacmc[0].size;
    printf("focscale = %f\n", focscale);
    eval_sepNBpt = (long) (eval_sepmaxld/eval_sepstepld);
    eval_contrastCurve = (double*) malloc(sizeof(double)*eval_sepNBpt);
    eval_COHcontrastCurve = (double*) malloc(sizeof(double)*eval_sepNBpt);
    eval_INCcontrastCurve = (double*) malloc(sizeof(double)*eval_sepNBpt);
    eval_contrastCurve_cnt = (double*) malloc(sizeof(double)*eval_sepNBpt);
    for(ri=0; ri<eval_sepNBpt; ri++)
    {
        eval_contrastCurve[ri] = 0.0;
        eval_COHcontrastCurve[ri] = 0.0;
        eval_INCcontrastCurve[ri] = 0.0;
        eval_contrastCurve_cnt[ri] = 0.0;
    }

    aveC = 0.0;
    aveC_INC = 0.0;
    aveC_COH = 0.0;
    aveCcnt = 0;


    IDrc = create_3Dimage_ID("_tmp_rc", xsize, ysize, zsize);
    rcarray = (float*) malloc(sizeof(float)*xsize*ysize*zsize);


    for(kk=0; kk<zsize; kk++)
        for(ii=0; ii<xsize; ii++)
            for(jj=0; jj<ysize; jj++)
            {
                xc = 1.0*ii-0.5*xsize;
                yc = 1.0*jj-0.5*ysize;
                xc *= focscale;
                yc *= focscale;
                rc = sqrt(xc*xc+yc*yc);
                data.image[IDrc].array.F[kk*xsize*ysize+jj*xsize+ii] = rc;

                ri = (long) (rc/eval_sepstepld-0.5);
                if(ri<0)
                    ri = 0;
                if(ri<eval_sepNBpt)
                {
                    eval_contrastCurve[ri] += data.image[IDps].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                    eval_COHcontrastCurve[ri] += data.image[IDps_COH].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                    eval_INCcontrastCurve[ri] += data.image[IDps_INC].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                    eval_contrastCurve_cnt[ri] += 1.0;
                }
                if((rc>2.0)&&(rc<6.0))
                {
                    aveC += data.image[IDps].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                    aveC_COH += data.image[IDps_COH].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                    aveC_INC += data.image[IDps_INC].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                    aveCcnt++;
                }
            }

	PIAACMCsimul_update_fnamedescr();
    sprintf(fname, "%s/ContrastCurveP_extsrc%2ld_sm%d.%s.txt", piaacmcsimul_var.piaacmcconfdir, (long) (-log10(ldoffset)*10.0+0.1), piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.fnamedescr);

    fp = fopen(fname, "w");
    drc = 0.2;
    rcbinmax = 100;
    for(rcbin=0; rcbin<rcbinmax; rcbin++)
    {
        rcmin = drc * rcbin;
        rcmax = drc * (rcbin+1);

        rc_cnt = 0;

        for(kk=0; kk<zsize; kk++)
            for(ii=0; ii<xsize; ii++)
                for(jj=0; jj<ysize; jj++)
                {
                    rc = data.image[IDrc].array.F[kk*xsize*ysize+jj*xsize+ii];
                    if((rc>rcmin)&&(rc<rcmax))
                    {
                        rcarray[rc_cnt] = data.image[IDps].array.F[kk*xsize*ysize+jj*xsize+ii]/avpeak;
                        rc_cnt++;
                    }
                }
        quick_sort_float(rcarray, rc_cnt);

        p50 = rcarray[(long) (0.5*rc_cnt)];
        p20 = rcarray[(long) (0.2*rc_cnt)];
        fprintf(fp, "%5.2f %12g %12g\n", drc*(rcbin+0.5), p50, p20);
    }


    free(rcarray);
    delete_image_ID("_tmp_rc");
    fclose(fp);

    delete_image_ID("starimINC");
    delete_image_ID("starimCOH");


	PIAACMCsimul_update_fnamedescr();
    sprintf(fname, "%s/ContrastCurve_extsrc%2ld_sm%d.%s.txt", piaacmcsimul_var.piaacmcconfdir, (long) (-log10(ldoffset)*10.0+0.1), piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.fnamedescr);


    fp = fopen(fname, "w");
    fprintf(fp, "# Contrast curves\n");
    fprintf(fp, "# col 1 : separation [l/D]\n");
    fprintf(fp, "# col 2 : contrast, total intensity, radially averaged\n");
    fprintf(fp, "# col 3 : contrast, Coherent component, radially averaged\n");
    fprintf(fp, "# col 4 : contrast, incoherent component, radially averaged\n");
    fprintf(fp, "\n");
    for(ri=0; ri<eval_sepNBpt; ri++)
    {
        eval_contrastCurve[ri] /= eval_contrastCurve_cnt[ri]+0.000001;
        eval_COHcontrastCurve[ri] /= eval_contrastCurve_cnt[ri]+0.000001;
        eval_INCcontrastCurve[ri] /= eval_contrastCurve_cnt[ri]+0.000001;
        fprintf(fp, "%10f %10g %10g %10g %10g\n", eval_sepstepld*ri, eval_contrastCurve[ri], eval_contrastCurve_cnt[ri], eval_COHcontrastCurve[ri], eval_INCcontrastCurve[ri]);
    }
    fclose(fp);

    free(eval_contrastCurve);
    free(eval_COHcontrastCurve);
    free(eval_INCcontrastCurve);
    free(eval_contrastCurve_cnt);

	PIAACMCsimul_update_fnamedescr();
    sprintf(fname, "%s/ContrastVal_extsrc%2ld_sm%d.%s.txt", piaacmcsimul_var.piaacmcconfdir, (long) (-log10(ldoffset)*10.0+0.1), piaacmcsimul_var.SCORINGMASKTYPE, piaacmcsimul_var.fnamedescr);
    fp = fopen(fname, "w");
	fprintf(fp, "# Contrast value\n");
	fprintf(fp, "# valref, aveC/aveCcnt, piaacmc[0].fpmaskradld, piaacmcsimul_var.PIAACMC_MASKRADLD, piaacmc[0].centObs0, piaacmc[0].centObs1, (long) (piaacmc[0].lambda*1e9), (long) (piaacmc[0].lambdaB+0.1), piaacmcsimul_var.PIAACMC_FPMsectors, piaacmc[0].NBrings, piaacmc[0].focmNBzone, piaacmc[0].nblambda, (long) (1000.0*ldoffset), piaacmc[0].fpmsagreg_coeff, piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode\n");
	fprintf(fp, "#\n");
    fprintf(fp, "%10g %10g %4.2f %4.2f %4.2f %4.2f %04ld %02ld %d %03ld %03ld %02d %03ld %7.3f %02d %d\n", valref, aveC/aveCcnt, piaacmc[0].fpmaskradld, piaacmcsimul_var.PIAACMC_MASKRADLD, piaacmc[0].centObs0, piaacmc[0].centObs1, (long) (piaacmc[0].lambda*1e9), (long) (piaacmc[0].lambdaB+0.1), piaacmcsimul_var.PIAACMC_FPMsectors, piaacmc[0].NBrings, piaacmc[0].focmNBzone, piaacmc[0].nblambda, (long) (1000.0*ldoffset), piaacmc[0].fpmsagreg_coeff, piaacmcsimul_var.computePSF_ResolvedTarget, piaacmcsimul_var.computePSF_ResolvedTarget_mode);
    fclose(fp);

    return 0;
}











