
# Design steps: Monochromatic design ( #step < 100 ) {#page_PIAACMCsimul_DesignSteps_Mono}




## 1. Overview


The design proceeds in discrete steps. The example scripts show the individual steps, and each step is executed with a separate command line with :

	./runPIAACMCdesign -m <STEPNUMBER+1>

Note that the above command will execute all steps up to and including step STEPNUMBER. For example, to run steps 0 to 5 (included), you must run `./runPIAACMCdesign -m 6`.

Steps from 0 to 99 are executed sequentially to design a monochromatic PIAACMC. For these steps, the user may run multiple steps with a single command. For example, running step 18 will execute all steps from 0 to 17 included. If a step has already been completed, it will not be re-run.







The monochromatic PIAACMC design process is as follows:

* design an idealized monochromatic PIAACMC for a centrally obscured aperture (steps 1-4)

* modify the design for the pupil aperture (steps 5-)


---

---



## 2. STEP 000 (MODE=0): Create an idealized centrally obscured apodized PIAACMC monochromatic design

This is meant as a starting point for the PIAACMC, which will then be optimized further. This step takes a few minutes, and upon normal completion, displays: 

~~~
./runPIAACMC REACHED STATE EXIT POINT (1)
~~~


The prolate function for a centrally obscured pupili is created, along with the idealized focal plane mask. The diffraction propagation code then computes complex amplitude in each plane.

This step takes a few minutes - most of the time is spent on iterations to compute the 2D apodization prolate function, executing discrete Fourier transforms (DFTs).

After this step, the contrast will likely be around 1e-7 to 1e-8. The next step to improve this nominal design is to find the optimal locations for the Lyot stops.


Output file                             | Description
----------------------------------------|--------------------------------------------------------------------------
ref/Cmodes_1024.fits                    | Cosine modes applied for PIAA optics shapes optimization
ref/Fmodes_1024.fits                    | Fourier modes applied for PIAA optics shapes optimization
piaacmcconf_i000/APLCapo.<X>.<Y>.info   | Apodization function info
piaacmcconf_i000/WFamp0_xxx.fits        | Amplitude in plane xxx
piaacmcconf_i000/WFpha0_xxx.fits        | Phase in plane xxx
piaacmcconf_i000/conjugations.txt       | List of planes and conjugation distance to reference
piaacmcconf_i000/apo2Drad.fits          | Amplitude apodization (entirely allocated to PIAA optics)
piaacmcconf_i000/PIAA_Mshapes.txt       | Aspheric optics shapes (r0 z0 r1 z1), unit [m]
piaacmcconf_i000/piaam0z.fits           | First PIAA mirror sag [m]
piaacmcconf_i000/piaam1z.fits           | Second PIAA mirror sag [m]
piaacmcconf_i000/psfi0_step000.fits     | Coronagraphic PSF (flux normalized to total=1 without coronagraph)
log/varlog.txt                          | List of coronagraph parameters



The conjugations.txt file should contain number, location and description of each plane:

~~~
00  0.000000    input pupil
01  0.000000    TT mirror
02  -1.102609   pupil plane apodizer
03  1.199997    PIAA optics 0
04  -1.102609   PIAA optics 1
05  -1.102609   opaque mask at PIAA elem 1
06  -1.102609   post focal plane mask pupil
07  0.000000    Lyot mask 0
08  0.000000    back end pupil stop  (rad = 0.920000)
~~~

The complex amplitude for each of these planes (files WFamp0_xxx.fits) should show the pupil apodization and coronagraphic effect induced by the focal plane mask which moves light into the central obstruction.



[Amplitude in each of the 9 planes](../src/PIAACMCsimul/doc/figures/examplePIAA_step0_WFamp.png)


The script provides some approximate performance metrics upon normal exit:

~~~
TOTAL = 0.002426
    FLUX   0    114332.0000 1.000000
    FLUX   1    114332.0000 1.000000
    FLUX   2    114331.9872 1.000000
    FLUX   3    114331.9668 1.000000
    FLUX   4    114331.9357 0.999999
    FLUX   5    114331.8634 0.999999
    FLUX   6     23794.1440 0.208114
    FLUX   7       551.2827 0.004822
    FLUX   8       277.3360 0.002426
COMPUTING UNRESOLVED SOURCE PSF -*- [0.000000 x 0.000000]
SCORINGTOTAL = 2436.000000  2436.000000
Peak constrast (rough estimate)= 26231.4 -> 2.00672e-06
optsyst[0].flux[0]  = 114332
SCORINGMASKTYPE = 0
[0] Total light in scoring field = 8.58033e+06, peak PSF = -1, SCOTINGTOTAL = 2436   -> Average contrast = 2.69458e-07
~~~


Runtime 00:01:36 (size 1024)


---


## STEP 001: Propagate solution and compute PSF

The previous solution is propagated. This should yield the same result as step 0. 

Output file                             | Description
----------------------------------------|--------------------------------------------------------------------------
piaacmcconf_i000/WFamp0_xxx.fits        | Amplitude in plane xxx
piaacmcconf_i000/WFpha0_xxx.fits        | Phase in plane xxx
piaacmcconf_i000/psfi0_step001.fits     | On-axis PSF 


~~~
TOTAL = 0.002387
    FLUX   0    114332.0000 1.000000
    FLUX   1    114332.0000 1.000000
    FLUX   2    114331.9872 1.000000
    FLUX   3    114331.9668 1.000000
    FLUX   4    114331.9356 0.999999
    FLUX   5    114331.9105 0.999999
    FLUX   6     23785.0351 0.208035
    FLUX   7       544.8861 0.004766
    FLUX   8       272.9253 0.002387
COMPUTING UNRESOLVED SOURCE PSF -*- [0.000000 x 0.000000]
SCORINGTOTAL = 2436.000000  2436.000000
Peak constrast (rough estimate)= 25740.7 -> 1.96918e-06
optsyst[0].flux[0]  = 114332
SCORINGMASKTYPE = 0
[0] Total light in scoring field = 8.18737e+06, peak PSF = -1, SCOTINGTOTAL = 2436   -> Average contrast = 2.57118e-07
saving contrast value [2.57118e-07] -> piaacmcconf_i000/contrast_ptsrc_sm0_s0_l0565_sr10_nbr001_mr000_minsag-10000_maxsag010000_fpmreg001000_ssr00_ssm0_Mirror_wb01.txt
~~~


Runtime 00:00:07 (size 1024)

Runtime 00:00:-- (size 2048)

---




## STEP 002: Specify input pupil geometry

The pupil geometry is copied to file `piaacmcconf_i000/pupa0_1024.fits`

Runtime 00:00:00 (size 1024)

Runtime 00:00:-- (size 2048)

---





## STEP 003 (mode = 0): compute on-axis PSF for new pupil geometry

Same as step 1, but taking into account the pupil geometry. With spiders or gaps in the pupil, the contrast will not be as good.

Output file                             | Description
----------------------------------------|--------------------------------------------------------------------------
piaacmcconf_i000/WFamp0_xxx.fits        | Amplitude in plane xxx
piaacmcconf_i000/WFpha0_xxx.fits        | Phase in plane xxx
piaacmcconf_i000/psfi0_step003.fits     | On-axis PSF 


~~~
TOTAL = 0.140551
    FLUX   0    103004.0000 1.000000
    FLUX   1    103004.0000 1.000000
    FLUX   2    103003.9828 1.000000
    FLUX   3    103003.9639 1.000000
    FLUX   4    103003.9359 0.999999
    FLUX   5    103003.6175 0.999996
    FLUX   6     27580.8543 0.267765
    FLUX   7     15763.2885 0.153036
    FLUX   8     14477.3372 0.140551
COMPUTING UNRESOLVED SOURCE PSF -*- [0.000000 x 0.000000]
SCORINGTOTAL = 2436.000000  2436.000000
Peak constrast (rough estimate)= 4.65338e+07 -> 0.00438592
optsyst[0].flux[0]  = 103004
SCORINGMASKTYPE = 0
[0] Total light in scoring field = 3.77881e+09, peak PSF = -1, SCOTINGTOTAL = 2436   -> Average contrast = 0.000146207
~~~

Runtime 00:00:07 (size 1024)

Runtile 00:01:-- (size 2048)

---


## STEP 004 (mode = 5): Compute Lyot stops shapes and locations, 1st pass

For circular centrally obscured pupils, the default Lyot stops configuration consists of two stops: one that masks the outer part of the beam, and one that masks the central obstruction.

Fine optimization of the stops locations is done with a separate command. The optional lsoptrange value is the range (unit: m) for the mask position search.


When the optimization completes, the best solutions are listed:

~~~
BEST SOLUTION: 5.646381230377 / 12.274741805168    0.000000942068 / 0.000021030879  -> 0.700625751251  0.039432599420
BEST SOLUTION: 5.891876066481 / 12.274741805168    0.000000942068 / 0.000021030879  -> 0.700355617330  0.039385771373
BEST SOLUTION: 8.837814099721 / 12.274741805168    0.000001130481 / 0.000021030879  -> 0.721432709799  0.039346632150
BEST SOLUTION: 9.083308935825 / 12.274741805168    0.000001130481 / 0.000021030879  -> 0.716777237262  0.038827027960
BEST SOLUTION: 9.328803771928 / 12.274741805168    0.000001130481 / 0.000021030879  -> 0.712835211313  0.038398704677
BEST SOLUTION: 9.574298608031 / 12.274741805168    0.000001130481 / 0.000021030879  -> 0.708858067020  0.037977986529
BEST SOLUTION: 9.819793444135 / 12.274741805168    0.000001130481 / 0.000021030879  -> 0.705292934047  0.037609975658
BEST SOLUTION: 10.065288280238 / 12.274741805168    0.000001130481 / 0.000021030879  -> 0.701842146347  0.037262905586
BEST SOLUTION: 11.292762460755 / 12.274741805168    0.000001356577 / 0.000021030879  -> 0.703163640577  0.036896957877
BEST SOLUTION: 12.029246969065 / 12.274741805168    0.000001627893 / 0.000021030879  -> 0.701096376349  0.036559137596
~~~


The optimal Lyot stop(s) conjugation(s) is written in file `piaacmcconf_i000/piaacmcparams_step004.conf`

Runtime 00:53:49 (size 1024)

Runtime 04:24:-- (size 2048)

---


## STEP 005 (mode = 2): Optimize focal plane mask transmission, 1st pass

See source code in PIAACMCsimul_exec_optimize_fpmtransmission()


Optimizes the focal plane mask transmission. Results are written in file `piaacmcconf_i000/result_fpmt.log`

~~~
 -0.05940301  3.24056e-05         0          0.3          0.1
 +0.04059699  1.92755e-05         0          0.3          0.1
 +0.14059699  9.55978e-06         0          0.3          0.1
 +0.24059699  3.25834e-06         0          0.3          0.1
 +0.34059699  3.71207e-07         0          0.3          0.1
 +0.44059699  8.98387e-07         0          0.3          0.1
 +0.54059699  4.83988e-06         0          0.3          0.1

 +0.25059699  2.81598e-06         1         0.09         0.03
 +0.28059699  1.69377e-06         1         0.09         0.03
 +0.31059699  8.78845e-07         1         0.09         0.03
 +0.34059699  3.71207e-07         1         0.09         0.03
 +0.37059699  1.70859e-07         1         0.09         0.03
 +0.40059699  2.77799e-07         1         0.09         0.03
 +0.43059699  6.92026e-07         1         0.09         0.03

 +0.34359699  3.37345e-07         2        0.027        0.009
 +0.35259699  2.54193e-07         2        0.027        0.009
 +0.36159699  1.98699e-07         2        0.027        0.009
 +0.37059699  1.70859e-07         2        0.027        0.009
 +0.37959699  1.70676e-07         2        0.027        0.009
 +0.38859699  1.98148e-07         2        0.027        0.009
 +0.39759699  2.53276e-07         2        0.027        0.009

 +0.37149699  1.69596e-07         3       0.0081       0.0027
 +0.37419699  1.67467e-07         3       0.0081       0.0027
 +0.37689699  1.67827e-07         3       0.0081       0.0027
 +0.37959699  1.70676e-07         3       0.0081       0.0027
 +0.38229699  1.76014e-07         3       0.0081       0.0027
 +0.38499699  1.83841e-07         3       0.0081       0.0027
 +0.38769699  1.94156e-07         3       0.0081       0.0027

 +0.37176699  1.69271e-07         4      0.00243      0.00081
 +0.37257699  1.68446e-07         4      0.00243      0.00081
 +0.37338699  1.67844e-07         4      0.00243      0.00081
 +0.37419699  1.67467e-07         4      0.00243      0.00081
 +0.37500699  1.67314e-07         4      0.00243      0.00081
 +0.37581699  1.67384e-07         4      0.00243      0.00081
 +0.37662699  1.67679e-07         4      0.00243      0.00081

 +0.37427799  1.67442e-07         5     0.000729     0.000243
 +0.37452099  1.67379e-07         5     0.000729     0.000243
 +0.37476399  1.67336e-07         5     0.000729     0.000243
 +0.37500699  1.67314e-07         5     0.000729     0.000243
 +0.37524999  1.67311e-07         5     0.000729     0.000243
 +0.37549299  1.67329e-07         5     0.000729     0.000243
 +0.37573599  1.67367e-07         5     0.000729     0.000243
~~~

Runtime 00:03:56 (size 1024)

Runtime 00:09 (size 2048)

---


## STEP 006 (mode = 5): Compute Lyot stops shapes and locations, 2nd pass, 70% throughput

Source code: 

~~~
BEST SOLUTION: 30.006465758591 / 50.010776264318    0.000000319995 / 0.000021330725  -> 0.602147984170  0.001772430096
BEST SOLUTION: 31.006681283877 / 50.010776264318    0.000000319995 / 0.000021330725  -> 0.601618355971  0.001755040537
BEST SOLUTION: 32.006896809163 / 50.010776264318    0.000000319995 / 0.000021330725  -> 0.601141999997  0.001739932914
BEST SOLUTION: 33.007112334450 / 50.010776264318    0.000000319995 / 0.000021330725  -> 0.600760087273  0.001728177742
BEST SOLUTION: 34.007327859736 / 50.010776264318    0.000000319995 / 0.000021330725  -> 0.600397821167  0.001717338698
BEST SOLUTION: 44.009483112599 / 50.010776264318    0.000000383994 / 0.000021330725  -> 0.602141259615  0.001712595775
BEST SOLUTION: 45.009698637886 / 50.010776264318    0.000000383994 / 0.000021330725  -> 0.601238947002  0.001692309270
BEST SOLUTION: 46.009914163172 / 50.010776264318    0.000000383994 / 0.000021330725  -> 0.600243611414  0.001670430013
~~~



Runtime 00:02:50 (size 1024)

Runtime 00:19:-- (size 2048)


---


## STEP 007 (mode = 40): Tune PIAA shapes and focal plane mask transm, 10 cosine modes, 5 Fourier modes



Progress can be tracked by watching file :

~~~
tail -f linoptval.txt
~~~

The file shows contrast value improving. Lines starting with '##' are evaluations along lines of steepest gradient, and lines starting with `#` are computations of local derivatives.

Last lines of file `linoptval.txt`:

~~~
##  [     7 /    20 ]   1.000      0.046011   2.78809e-08   (reg =  5.68825e-09 [1]  1.99516e-14 [1]   contrast =          2.21926e-08)       [0] [31] bestval =  2.74617e-08
##  [     7 /    20 ]   1.000      0.056413   2.78643e-08   (reg =  5.69589e-09 [1]  1.99516e-14 [1]   contrast =          2.21684e-08)       [0] [31] bestval =  2.74617e-08
##  [     7 /    20 ]   1.000      0.068895   2.78449e-08   (reg =  5.70507e-09 [1]  1.99516e-14 [1]   contrast =          2.21398e-08)       [0] [31] bestval =  2.74617e-08
##  [     7 /    20 ]   1.000      0.083874   2.78232e-08   (reg =   5.7161e-09 [1]  1.99516e-14 [1]   contrast =          2.21071e-08)       [0] [31] bestval =  2.74617e-08
##  [     7 /    20 ]   1.000      0.101849   2.77964e-08   (reg =  5.72934e-09 [1]  1.99516e-14 [1]   contrast =          2.20671e-08)       [0] [31] bestval =  2.74617e-08
##  [     7 /    20 ]   1.000      0.123419   2.77655e-08   (reg =  5.74525e-09 [1]  1.99516e-14 [1]   contrast =          2.20203e-08)       [0] [31] bestval =  2.74617e-08
##  [     7 /    20 ]   1.000      0.149303   2.77303e-08   (reg =  5.76437e-09 [1]  1.99516e-14 [1]   contrast =          2.19659e-08)       [0] [31] bestval =  2.74617e-08
##  [     7 /    20 ]   1.000      0.180363   2.76929e-08   (reg =  5.78736e-09 [1]  1.99516e-14 [1]   contrast =          2.19055e-08)       [0] [31] bestval =  2.74617e-08
##  [     7 /    20 ]   1.000      0.217636   2.76508e-08   (reg =  5.81502e-09 [1]  1.99516e-14 [1]   contrast =          2.18358e-08)       [0] [31] bestval =  2.74617e-08
##  [     7 /    20 ]   1.000      0.262363   2.76037e-08   (reg =  5.84828e-09 [1]  1.99516e-14 [1]   contrast =          2.17553e-08)       [0] [31] bestval =  2.74617e-08
##  [     7 /    20 ]   1.000      0.316036   2.75583e-08   (reg =  5.88834e-09 [1]  1.99516e-14 [1]   contrast =          2.16699e-08)       [0] [31] bestval =  2.74617e-08
##  [     7 /    20 ]   1.000      0.380443   2.75127e-08   (reg =  5.93659e-09 [1]  1.99516e-14 [1]   contrast =          2.15761e-08)       [0] [31] bestval =  2.74617e-08
##  [     7 /    20 ]   1.000      0.457732   2.74739e-08   (reg =  5.99475e-09 [1]  1.99516e-14 [1]   contrast =          2.14791e-08)       [0] [31] bestval =  2.74617e-08
##  [     7 /    20 ]   1.000      0.550478   2.74509e-08   (reg =  6.06494e-09 [1]  1.99516e-14 [1]   contrast =          2.13859e-08)       [0] [31]  -> BEST VECTOR =======
##  [     7 /    20 ]   1.000      0.661774   2.74541e-08   (reg =  6.14971e-09 [1]  1.99516e-14 [1]   contrast =          2.13043e-08)       [0] [31] bestval =  2.74509e-08
->     7             2.74509e-08 <-          2.14234e-07 
~~~

Runtime 01:11:45 (size 1024)

Runtime 02:20:-- (size 2048)

----


## STEP 008 (mode = 40): Tune PIAA shapes and focal plane mask transm,  20 cosine modes, 20 Fourier modes


Progress can be tracked by watching file :

~~~
tail -f linoptval.txt
~~~

Last lines of file `linoptval.txt`:

~~~
##  [     3 /  1000 ]   0.600      0.002400   1.41799e-08   (reg =   7.6046e-09 [1]  1.99516e-14 [1]   contrast =          6.57528e-09)       [0] [81] bestval =  1.41792e-08
##  [     3 /  1000 ]   0.800      0.001000   1.41798e-08   (reg =  7.60463e-09 [1]  1.99516e-14 [1]   contrast =          6.57515e-09)       [0] [81] bestval =  1.41792e-08
##  [     3 /  1000 ]   0.800      0.002400   1.41796e-08   (reg =  7.60488e-09 [1]  1.99516e-14 [1]   contrast =          6.57471e-09)       [0] [81] bestval =  1.41792e-08
##  [     3 /  1000 ]   0.800      0.004080   1.41791e-08   (reg =  7.60517e-09 [1]  1.99516e-14 [1]   contrast =          6.57389e-09)       [0] [81]  -> BEST VECTOR =======
##  [     3 /  1000 ]   0.800      0.006096   1.41793e-08   (reg =  7.60552e-09 [1]  1.99516e-14 [1]   contrast =          6.57372e-09)       [0] [81] bestval =  1.41791e-08
##  [     3 /  1000 ]   1.000      0.001000     1.418e-08   (reg =  7.60475e-09 [1]  1.99516e-14 [1]   contrast =          6.57521e-09)       [0] [81] bestval =  1.41791e-08
##  [     3 /  1000 ]   1.000      0.002400   1.41794e-08   (reg =  7.60516e-09 [1]  1.99516e-14 [1]   contrast =          6.57427e-09)       [0] [81] bestval =  1.41791e-08
##  [     3 /  1000 ]   1.000      0.004080    1.4179e-08   (reg =  7.60565e-09 [1]  1.99516e-14 [1]   contrast =          6.57337e-09)       [0] [81]  -> BEST VECTOR =======
##  [     3 /  1000 ]   1.000      0.006096   1.41796e-08   (reg =  7.60623e-09 [1]  1.99516e-14 [1]   contrast =          6.57337e-09)       [0] [81] bestval =   1.4179e-08
->     3              1.4179e-08 <-          2.74509e-08 
~~~

Runtime 00:34:35 (size 1024)

Runtime 04:38:-- (size 2048)

---


## STEP 009 (mode = 5): Compute Lyot stops shapes and locations, 2nd pass, 70% throughput


~~~
BEST SOLUTION: 0.000000000000 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.654580018549  0.003217549471
BEST SOLUTION: 4.482217103616 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.654559468105  0.003212765094
BEST SOLUTION: 4.930438813977 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.654514836838  0.003203476539
BEST SOLUTION: 5.378660524339 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.654457119322  0.003192321029
BEST SOLUTION: 5.826882234700 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.654421840552  0.003186018405
BEST SOLUTION: 6.275103945062 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.654381291083  0.003179298807
BEST SOLUTION: 6.723325655424 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.654282076189  0.003164067279
BEST SOLUTION: 7.171547365785 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.654214540383  0.003154407998
BEST SOLUTION: 7.619769076147 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.654100992460  0.003139033972
BEST SOLUTION: 8.067990786508 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.654018192997  0.003128464132
BEST SOLUTION: 8.516212496870 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.653877197371  0.003111544800
BEST SOLUTION: 8.964434207231 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.653767183605  0.003098887074
BEST SOLUTION: 9.412655917593 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.653635689943  0.003084587083
BEST SOLUTION: 9.860877627955 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.653461347660  0.003066496777
BEST SOLUTION: 10.309099338316 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.653316127525  0.003052160309
BEST SOLUTION: 10.757321048678 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.653059187548  0.003027656847
BEST SOLUTION: 11.205542759039 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.652798175033  0.003003926502
BEST SOLUTION: 11.653764469401 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.652556795472  0.002982867822
BEST SOLUTION: 12.101986179762 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.652348492130  0.002965292945
BEST SOLUTION: 12.550207890124 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.652048657232  0.002940968576
BEST SOLUTION: 12.998429600486 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.651596715831  0.002905644761
BEST SOLUTION: 13.446651310847 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.651116702176  0.002869351600
BEST SOLUTION: 13.894873021209 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.650737856102  0.002841584172
BEST SOLUTION: 14.343094731570 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.650347660667  0.002813970087
BEST SOLUTION: 14.791316441932 / 22.411085518079    0.000000642797 / 0.000020663874  -> 0.650006287003  0.002790509936
BEST SOLUTION: 20.618198676632 / 22.411085518079    0.000000771357 / 0.000020663874  -> 0.653227445945  0.002786372082
BEST SOLUTION: 21.066420386994 / 22.411085518079    0.000000771357 / 0.000020663874  -> 0.652015770540  0.002728302947
BEST SOLUTION: 21.514642097355 / 22.411085518079    0.000000771357 / 0.000020663874  -> 0.651040104043  0.002682478046
BEST SOLUTION: 23.307528938802 / 22.411085518079    0.000000925628 / 0.000020663874  -> 0.650876386774  0.002653676594
~~~


Runtime 00:02:44 (size 1024)

Runtime 00:23:-- (size 2048)

---


## STEP 010 (mode = 1): Tune Lyot stops conjugations

To view result:

~~~
tail -f result_LMpos.log
~~~

Output:

~~~
 -3.059896 2.61744e-06
 -2.459896 1.93909e-06
 -1.859896 1.3068e-06
 -1.259896 7.31229e-07
 -0.659896 2.27333e-07
 -0.059896 7.75573e-09
 0.540104 2.56943e-07
 1.140104 7.75529e-07
 1.740104 1.42284e-06
 2.340104 2.10727e-06

 -0.959896 4.65961e-07
 -0.659896 2.27333e-07
 -0.359896 8.27041e-08
 -0.059896 7.75573e-09
 0.240104 9.46098e-08
 0.540104 2.56943e-07
 0.840104 4.87981e-07

 -0.329896 6.78089e-08
 -0.239896 3.03972e-08
 -0.149896 1.15237e-08
 -0.059896 7.75573e-09
 0.030104 1.16957e-08
 0.120104 2.55954e-08
 0.210104 7.44832e-08

 -0.140896 1.07657e-08
 -0.113896 9.08773e-09
 -0.086896 8.00901e-09
 -0.059896 7.75573e-09
 -0.032896 8.53905e-09
 -0.005896 9.88908e-09
~~~



~~~
TOTAL = 0.002857
    FLUX   0    103004.0000 1.000000
    FLUX   1    103004.0000 1.000000
    FLUX   2    103003.9828 1.000000
    FLUX   3    103003.9639 1.000000
    FLUX   4    103003.9360 0.999999
    FLUX   5    103003.6270 0.999996
    FLUX   6     34709.3163 0.336971
    FLUX   7       308.9348 0.002999
    FLUX   8       294.2806 0.002857
COMPUTING UNRESOLVED SOURCE PSF -*- [0.000000 x 0.000000]
SCORINGTOTAL = 2436.000000  2436.000000
Peak constrast (rough estimate)= 4750.5 -> 4.47746e-07
optsyst[0].flux[0]  = 103004
SCORINGMASKTYPE = 0
[0] Total light in scoring field = 255588, peak PSF = -1, SCOTINGTOTAL = 2436   -> Average contrast = 9.88908e-09
BEST SOLUTION :   -0.059896 7.75573e-09
~~~


Runtime 00:00:54 (size 1024)

Runtime 00:08:-- (size 2048)

---




## STEP 011 (mode = 40): Tune PIAA shapes and focal plane mask transm,  20 cosine modes, 20 Fourier modes

To view result:

~~~
tail -f result_LMpos.log
~~~

Output:

~~~
##  [     3 /    20 ]   0.000      0.001000   1.44695e-08   (reg =  7.71092e-09 [1]  1.99516e-14 [1]   contrast =          6.75855e-09)       [0] [81] ===== START POINT =====
##  [     3 /    20 ]   0.000      0.002400   1.44694e-08   (reg =  7.71092e-09 [1]  1.99516e-14 [1]   contrast =          6.75851e-09)       [0] [81]  -> BEST VECTOR =======
##  [     3 /    20 ]   0.000      0.004080   1.44695e-08   (reg =  7.71092e-09 [1]  1.99516e-14 [1]   contrast =          6.75859e-09)       [0] [81] bestval =  1.44694e-08
##  [     3 /    20 ]   0.200      0.001000   1.44696e-08   (reg =  7.71092e-09 [1]  1.99516e-14 [1]   contrast =          6.75864e-09)       [0] [81] bestval =  1.44694e-08
##  [     3 /    20 ]   0.200      0.002400   1.44696e-08   (reg =  7.71092e-09 [1]  1.99516e-14 [1]   contrast =          6.75864e-09)       [0] [81] bestval =  1.44694e-08
##  [     3 /    20 ]   0.400      0.001000   1.44696e-08   (reg =  7.71092e-09 [1]  1.99516e-14 [1]   contrast =          6.75866e-09)       [0] [81] bestval =  1.44694e-08
##  [     3 /    20 ]   0.400      0.002400   1.44695e-08   (reg =  7.71092e-09 [1]  1.99516e-14 [1]   contrast =          6.75852e-09)       [0] [81] bestval =  1.44694e-08
##  [     3 /    20 ]   0.400      0.004080   1.44695e-08   (reg =  7.71092e-09 [1]  1.99516e-14 [1]   contrast =          6.75855e-09)       [0] [81] bestval =  1.44694e-08
##  [     3 /    20 ]   0.600      0.001000   1.44696e-08   (reg =  7.71093e-09 [1]  1.99516e-14 [1]   contrast =          6.75864e-09)       [0] [81] bestval =  1.44694e-08
##  [     3 /    20 ]   0.600      0.002400   1.44696e-08   (reg =  7.71096e-09 [1]  1.99516e-14 [1]   contrast =          6.75863e-09)       [0] [81] bestval =  1.44694e-08
##  [     3 /    20 ]   0.800      0.001000   1.44695e-08   (reg =  7.71096e-09 [1]  1.99516e-14 [1]   contrast =          6.75856e-09)       [0] [81] bestval =  1.44694e-08
##  [     3 /    20 ]   0.800      0.002400   1.44701e-08   (reg =  7.71104e-09 [1]  1.99516e-14 [1]   contrast =          6.75903e-09)       [0] [81] bestval =  1.44694e-08
##  [     3 /    20 ]   1.000      0.001000   1.44698e-08   (reg =  7.71099e-09 [1]  1.99516e-14 [1]   contrast =          6.75881e-09)       [0] [81] bestval =  1.44694e-08
##  [     3 /    20 ]   1.000      0.002400   1.44703e-08   (reg =   7.7111e-09 [1]  1.99516e-14 [1]   contrast =          6.75918e-09)       [0] [81] bestval =  1.44694e-08
->     3             1.44694e-08 <-          1.53614e-08
~~~


~~~
TOTAL = 0.002730
    FLUX   0    103004.0000 1.000000
    FLUX   1    103004.0000 1.000000
    FLUX   2    103003.9828 1.000000
    FLUX   3    103003.9639 1.000000
    FLUX   4    103003.9359 0.999999
    FLUX   5    103003.6261 0.999996
    FLUX   6     34717.4800 0.337050
    FLUX   7       298.1019 0.002894
    FLUX   8       281.2208 0.002730
COMPUTING UNRESOLVED SOURCE PSF -*- [0.000000 x 0.000000]
SCORINGTOTAL = 2436.000000  2436.000000
Peak constrast (rough estimate)= 3113.01 -> 2.93408e-07
optsyst[0].flux[0]  = 103004
SCORINGMASKTYPE = 0
[0] Total light in scoring field = 174677, peak PSF = -1, SCOTINGTOTAL = 2436   -> Average contrast = 6.75851e-09

~~~


Runtime 00:37:29 (size 1024)

Runtime 03:19:00 (size 2048)

---






## STEP 012 (mode = 40): Tune PIAA shapes and focal plane mask transm,  40 cosine modes, 150 Fourier modes


To view result:

~~~
tail -f result_LMpos.log
~~~

Output:

~~~
##  [     3 /    10 ]   0.000      0.001000   1.20892e-08   (reg =  7.67807e-09 [1]  1.99516e-14 [1]   contrast =          4.41108e-09)       [0] [281] ===== START POINT =====
##  [     3 /    10 ]   0.000      0.002400   1.20892e-08   (reg =  7.67807e-09 [1]  1.99516e-14 [1]   contrast =          4.41106e-09)       [0] [281]  -> BEST VECTOR =======
##  [     3 /    10 ]   0.000      0.004080   1.20892e-08   (reg =  7.67807e-09 [1]  1.99516e-14 [1]   contrast =          4.41106e-09)       [0] [281] bestval =  1.20892e-08
##  [     3 /    10 ]   0.200      0.001000   1.20892e-08   (reg =  7.67807e-09 [1]  1.99516e-14 [1]   contrast =          4.41112e-09)       [0] [281] bestval =  1.20892e-08
##  [     3 /    10 ]   0.200      0.002400   1.20892e-08   (reg =  7.67807e-09 [1]  1.99516e-14 [1]   contrast =           4.4111e-09)       [0] [281] bestval =  1.20892e-08
##  [     3 /    10 ]   0.200      0.004080   1.20892e-08   (reg =  7.67807e-09 [1]  1.99516e-14 [1]   contrast =          4.41107e-09)       [0] [281] bestval =  1.20892e-08
##  [     3 /    10 ]   0.200      0.006096   1.20891e-08   (reg =  7.67807e-09 [1]  1.99516e-14 [1]   contrast =          4.41103e-09)       [0] [281]  -> BEST VECTOR =======
##  [     3 /    10 ]   0.200      0.008515   1.20891e-08   (reg =  7.67807e-09 [1]  1.99516e-14 [1]   contrast =          4.41105e-09)       [0] [281] bestval =  1.20891e-08
##  [     3 /    10 ]   0.400      0.001000   1.20892e-08   (reg =  7.67807e-09 [1]  1.99516e-14 [1]   contrast =          4.41109e-09)       [0] [281] bestval =  1.20891e-08
##  [     3 /    10 ]   0.400      0.002400   1.20892e-08   (reg =  7.67807e-09 [1]  1.99516e-14 [1]   contrast =          4.41107e-09)       [0] [281] bestval =  1.20891e-08
##  [     3 /    10 ]   0.400      0.004080   1.20892e-08   (reg =  7.67807e-09 [1]  1.99516e-14 [1]   contrast =          4.41107e-09)       [0] [281] bestval =  1.20891e-08
##  [     3 /    10 ]   0.600      0.001000   1.20894e-08   (reg =  7.67806e-09 [1]  1.99516e-14 [1]   contrast =          4.41128e-09)       [0] [281] bestval =  1.20891e-08
##  [     3 /    10 ]   0.600      0.002400   1.20891e-08   (reg =  7.67805e-09 [1]  1.99516e-14 [1]   contrast =          4.41104e-09)       [0] [281]  -> BEST VECTOR =======
##  [     3 /    10 ]   0.600      0.004080   1.20892e-08   (reg =  7.67803e-09 [1]  1.99516e-14 [1]   contrast =          4.41113e-09)       [0] [281] bestval =  1.20891e-08
##  [     3 /    10 ]   0.800      0.001000   1.20891e-08   (reg =  7.67805e-09 [1]  1.99516e-14 [1]   contrast =          4.41107e-09)       [0] [281] bestval =  1.20891e-08
##  [     3 /    10 ]   0.800      0.002400    1.2089e-08   (reg =  7.67801e-09 [1]  1.99516e-14 [1]   contrast =          4.41097e-09)       [0] [281]  -> BEST VECTOR =======
##  [     3 /    10 ]   0.800      0.004080   1.20893e-08   (reg =  7.67796e-09 [1]  1.99516e-14 [1]   contrast =          4.41131e-09)       [0] [281] bestval =   1.2089e-08
##  [     3 /    10 ]   1.000      0.001000   1.20891e-08   (reg =  7.67802e-09 [1]  1.99516e-14 [1]   contrast =           4.4111e-09)       [0] [281] bestval =   1.2089e-08
##  [     3 /    10 ]   1.000      0.002400   1.20893e-08   (reg =  7.67796e-09 [1]  1.99516e-14 [1]   contrast =          4.41135e-09)       [0] [281] bestval =   1.2089e-08
->     3              1.2089e-08 <-          1.44694e-08 
~~~


~~~
TOTAL = 0.002739
    FLUX   0    103004.0000 1.000000
    FLUX   1    103004.0000 1.000000
    FLUX   2    103003.9828 1.000000
    FLUX   3    103003.9639 1.000000
    FLUX   4    103003.9361 0.999999
    FLUX   5    103003.6283 0.999996
    FLUX   6     34705.9637 0.336938
    FLUX   7       299.4217 0.002907
    FLUX   8       282.1562 0.002739
COMPUTING UNRESOLVED SOURCE PSF -*- [0.000000 x 0.000000]
SCORINGTOTAL = 2436.000000  2436.000000
Peak constrast (rough estimate)= 806.251 -> 7.5991e-08
optsyst[0].flux[0]  = 103004
SCORINGMASKTYPE = 0
[0] Total light in scoring field = 114004, peak PSF = -1, SCOTINGTOTAL = 2436   -> Average contrast = 4.41097e-09

~~~


Runtime 01:31:26 (size 1024)

Runtime -- (size 2048)

---



## STEP 013 (mode = 5): Compute Lyot stops shapes and locations, 3nd pass, 70% throughput

~~~
BEST SOLUTION: 9.601575812043 / 15.002462206317    0.000001111027 / 0.000020668962  -> 0.704840156107  0.005970769686
BEST SOLUTION: 9.901625056169 / 15.002462206317    0.000001111027 / 0.000020668962  -> 0.704051224148  0.005889714915
BEST SOLUTION: 10.201674300296 / 15.002462206317    0.000001111027 / 0.000020668962  -> 0.703399562778  0.005824872470
BEST SOLUTION: 10.501723544422 / 15.002462206317    0.000001111027 / 0.000020668962  -> 0.702606826377  0.005748303532
BEST SOLUTION: 10.801772788548 / 15.002462206317    0.000001111027 / 0.000020668962  -> 0.701929916566  0.005684742522
BEST SOLUTION: 11.101822032675 / 15.002462206317    0.000001111027 / 0.000020668962  -> 0.701368943834  0.005633479310
BEST SOLUTION: 11.401871276801 / 15.002462206317    0.000001111027 / 0.000020668962  -> 0.700809895288  0.005583827596
BEST SOLUTION: 11.701920520927 / 15.002462206317    0.000001111027 / 0.000020668962  -> 0.700082557238  0.005520894724
BEST SOLUTION: 14.402363718064 / 15.002462206317    0.000001333232 / 0.000020668962  -> 0.700753966114  0.005407358350
~~~

Runtime 00:03:23 (size 1024)

---



## STEP 014 (mode = 1): Tune Lyot stops conjugations

~~~
TOTAL = 0.003362
    FLUX   0    103004.0000 1.000000
    FLUX   1    103004.0000 1.000000
    FLUX   2    103003.9828 1.000000
    FLUX   3    103003.9639 1.000000
    FLUX   4    103003.9361 0.999999
    FLUX   5    103003.6283 0.999996
    FLUX   6     34705.9637 0.336938
    FLUX   7       371.9539 0.003611
    FLUX   8       346.3072 0.003362
COMPUTING UNRESOLVED SOURCE PSF -*- [0.000000 x 0.000000]
SCORINGTOTAL = 2436.000000  2436.000000
Peak constrast (rough estimate)= 2900.95 -> 2.73421e-07
optsyst[0].flux[0]  = 103004
SCORINGMASKTYPE = 0
[0] Total light in scoring field = 582659, peak PSF = -1, SCOTINGTOTAL = 2436   -> Average contrast = 2.25439e-08
~~~

Runtime 00:01:16 (size 1024)

---




## STEP 015 (mode = 40): tune PIAA shapes and focal plane mask transm, 20 cosine modes, 20 Fourier modes


To view result:

~~~
tail -f result_LMpos.log
~~~

Output:

~~~
##  [     3 /    20 ]   0.000      0.001000   1.57523e-08   (reg =   7.8125e-09 [1]  1.99516e-14 [1]   contrast =          7.93981e-09)       [0] [81] ===== START POINT =====
##  [     3 /    20 ]   0.000      0.002400   1.57524e-08   (reg =   7.8125e-09 [1]  1.99516e-14 [1]   contrast =          7.93983e-09)       [0] [81] bestval =  1.57523e-08
##  [     3 /    20 ]   0.200      0.001000   1.57523e-08   (reg =   7.8125e-09 [1]  1.99516e-14 [1]   contrast =           7.9398e-09)       [0] [81]  -> BEST VECTOR =======
##  [     3 /    20 ]   0.200      0.002400   1.57523e-08   (reg =   7.8125e-09 [1]  1.99516e-14 [1]   contrast =          7.93982e-09)       [0] [81] bestval =  1.57523e-08
##  [     3 /    20 ]   0.400      0.001000   1.57523e-08   (reg =   7.8125e-09 [1]  1.99516e-14 [1]   contrast =          7.93978e-09)       [0] [81]  -> BEST VECTOR =======
##  [     3 /    20 ]   0.400      0.002400   1.57523e-08   (reg =   7.8125e-09 [1]  1.99516e-14 [1]   contrast =          7.93982e-09)       [0] [81] bestval =  1.57523e-08
##  [     3 /    20 ]   0.600      0.001000   1.57527e-08   (reg =   7.8125e-09 [1]  1.99516e-14 [1]   contrast =          7.94014e-09)       [0] [81] bestval =  1.57523e-08
##  [     3 /    20 ]   0.600      0.002400   1.57527e-08   (reg =  7.81252e-09 [1]  1.99516e-14 [1]   contrast =          7.94012e-09)       [0] [81] bestval =  1.57523e-08
##  [     3 /    20 ]   0.600      0.004080   1.57528e-08   (reg =  7.81254e-09 [1]  1.99516e-14 [1]   contrast =          7.94025e-09)       [0] [81] bestval =  1.57523e-08
##  [     3 /    20 ]   0.800      0.001000   1.57527e-08   (reg =  7.81252e-09 [1]  1.99516e-14 [1]   contrast =          7.94015e-09)       [0] [81] bestval =  1.57523e-08
##  [     3 /    20 ]   0.800      0.002400   1.57531e-08   (reg =  7.81256e-09 [1]  1.99516e-14 [1]   contrast =          7.94051e-09)       [0] [81] bestval =  1.57523e-08
##  [     3 /    20 ]   1.000      0.001000   1.57528e-08   (reg =  7.81254e-09 [1]  1.99516e-14 [1]   contrast =          7.94027e-09)       [0] [81] bestval =  1.57523e-08
##  [     3 /    20 ]   1.000      0.002400   1.57524e-08   (reg =  7.81261e-09 [1]  1.99516e-14 [1]   contrast =          7.93977e-09)       [0] [81] bestval =  1.57523e-08
##  [     3 /    20 ]   1.000      0.004080    1.5753e-08   (reg =  7.81268e-09 [1]  1.99516e-14 [1]   contrast =          7.94025e-09)       [0] [81] bestval =  1.57523e-08
->     3             1.57523e-08 <-          1.72109e-08 
~~~


~~~
TOTAL = 0.003110
    FLUX   0    103004.0000 1.000000
    FLUX   1    103004.0000 1.000000
    FLUX   2    103003.9828 1.000000
    FLUX   3    103003.9639 1.000000
    FLUX   4    103003.9359 0.999999
    FLUX   5    103003.6294 0.999996
    FLUX   6     34709.1463 0.336969
    FLUX   7       348.9913 0.003388
    FLUX   8       320.3876 0.003110
COMPUTING UNRESOLVED SOURCE PSF -*- [0.000000 x 0.000000]
SCORINGTOTAL = 2436.000000  2436.000000
Peak constrast (rough estimate)= 1003.22 -> 9.45553e-08
optsyst[0].flux[0]  = 103004
SCORINGMASKTYPE = 0
[0] Total light in scoring field = 205208, peak PSF = -1, SCOTINGTOTAL = 2436   -> Average contrast = 7.93978e-09
~~~

Runtime 00:32:31 (size 1024)

---






## STEP 016 (mode = 40): tune PIAA shapes and focal plane mask transm, 40 cosine modes, 150 Fourier modes
In this example step 17 is skipped by typing:

~~~
touch piaacmcconf_i000/step016.txt
~~~



## STEP 017 (mode = 40): tune PIAA shapes and focal plane mask transm, 40 cosine modes, 625 Fourier modes
In this example step 17 is skipped by typing:

~~~
touch piaacmcconf_i000/step017.txt
~~~
