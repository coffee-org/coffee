

# C code description {#page_PIAACMCsimul_Ccode}

Source code is in file PIAACMCsimul.c

## PIAACMC_designcodes


The main function in the source code is PIAACMCsimul_exec(), which takes two arguments: the configuration index (usually a 3 digit integer) and the mode (integer) which describes the operation to be performed to the PIAACMC design.

Mode		| Description
------------|---------------------------------------------------------------------------------------------------------------------
0           | Compute on-axis propagation for specified configuration. If configuration does not exist, create idealized monochromatic PIAACMC (intended to create a new index) and compute on-axis propagation
1           | Optimize Lyot stop(s) locations (scan) 
2           | Optimize focal plane mask transmission for idealized monochromatic PIAACMC (scan)
3           | Run without focal plane mask (for testing and calibration)
4           | Linear optimization around current design, free parameters = PIAA optics cosines shapes (track progress by looking at val.opt file)
5           | Optimize Lyot stops shapes and positions
10          | Setup polychromatic optimization
11          | Compute polychromatic response to zones, store result in FPMresp
12          | Search for best mask solution using FPMresp, random search
13          | Optimize PIAA optics shapes and focal plane mask transmission (idealized PIAACMC)
40          | Optimize PIAA optics shapes and focal plane mask transmission (idealized PIAACMC)
41          | Optimize PIAA optics shapes and focal plane mask zones (polychromatic)
100         | Evaluate current design: polychromatic contrast, pointing sensitivity
101         | Transmission as a function of angular separation
200         | Make focal plane mask OPD map




The following variables can be set :

Variable		    | Description
--------------------|-------------------------------------
PIAACMC_size        | Array size (1024, 2048, etc...)
PIAACMC_pixscale    | Pixel scale [m]
PIAACMC_dftgrid     | Sampling interval in DFTs
PIAACMC_centobs0    | Input central obstruction
PIAACMC_centobs1    | Output central obstruction
PIAACMC_nblambda    | Number of wavelength points
PIAACMC_resolved    | 1 if resolved source (3 points at r = 0.01 l/D, 120 deg apart)
PIAACMC_fpmtype     | 1 if physical mask, 0 if idealized mask
PIAACMC_FPMsectors  | Number of sectors in focal plane mask
PIAACMC_NBrings     | Number of rings in focal plane mask
PIAACMC_fpmradld    | Focal plane mask outer radius



## Initialization rules (function  PIAAsimul_initpiaacmc() )


* if configuration directory exists, use it and load configuration file ( function  PIAAsimul_loadpiaacmcconf ), otherwise, create it

* load/create Cmodes 

* load/create Fmodes

* load mode coefficients for piaa shapes if they exist. If not:
	
	* create radial apodization for centrally obscured idealized monochromatic PIAACMC
	
	* fit / extrapolate radial apodization profile with cosines
	
	* using above fit, create 2D radial sag for both PIAA optics ( -> PIAA_Mshapes.txt)
	
	* make 2D sag maps for both optics ( -> piaa0z.fits, piaa1z.fits)
	
	* fit 2D sag maps with Cmodes and Fmodes coefficients ( -> piaa0Cmodes, piaa0Fmodes, piaa1Cmodes, piaa1Fmodes )

* load/create focal plane mask zone map. This is the map that defines the geometry (which ring is where)

* load/create focal plane mask thickness array

* load/create focal plane mask transmission array

* load/create Lyot stops






### Code breakdown

- PIAACMCsimul_exec() :
	- PIAAsimul_initpiaacmcconf(): Load/Creates/initializes piaacmcconf structure and directory
		- perform default initialization
		- PIAAsimul_loadpiaacmcconf(): Loading PIAACMC configuration from "piaacmcconfxxx/piaacmcparams.conf" if it exists
		- Creating/loading Cmodes and Fmodes
		- IMPORT / CREATE PIAA SHAPES
			- create 2D prolate iteratively
			- PIAACMCsimul_load2DRadialApodization():fit PIAA shapes with Cosine modes
			- PIAACMCsimul_init_geomPIAA_rad(): compute radial PIAA sag from cosine apodization fit
			- PIAACMCsimul_mkPIAAMshapes_from_RadSag(): Make 2D sag shapes from radial PIAA sag
		- MAKE FOCAL PLANE MASK
		- MAKE LYOT STOPS
		- PIAAsimul_savepiaacmcconf(): save configuration
	- PIAACMCsimul_makePIAAshapes(): construct PIAA shapes from fitting coefficients
		- construct 2D PIAA mirror shapes from piaa0Cmodescoeff, piaa0Fmodescoeff, piaa1Cmodescoeff, piaa1Fmodescoeff
	- PIAACMCsimul_computePSF(): Compute PSF

### Output Files


APLCmaskCtransm.txt ??
fpm_ampl.fits
fpm_pha.fits
FPmask.tmp.fits 


#### Configuration : 

Output file	                       | Notes
-----------------------------------|-------------------------------------
./piaaconfxxx/piaacmcparams.conf   | Configuration parameters
./piaaconfxxx/conjugations.txt	   | Conjugations
./piaaconfxxx/lambdalist.txt       | list of wavelength values
./piaaconfxxx/pupa0_[size].fits	   | input pupil (created by default if does not exist)

#### Wavefront Modes :

Output file                      | Notes
---------------------------------|-------------------------------------
Cmodes.fits	                     | circular radial cosine modes (40 modes, hard coded)
Fmodes.fits	                     | Fourier modes (625 modes = 10 CPA, hard coded)
./piaaconfxxx/ModesExpr_CPA.txt  | modes definition
./piaaconfxxx/APOmodesCos.fits	 | Cosine modes for fitting 2D apodization profile


#### PIAA mirrors, apodization, fits:

Output file	                                | Notes
--------------------------------------------|-------------------------------------
./piaaconfxxx/APLCapo.1.400.0.300.info      | file written by prolate generation function coronagraph_make_2Dprolate in coronagraphs.c
./piaaconfxxx/apo2Drad.fits	                | idealized PIAACMC 2D apodization
./piaaconfxxx/piaam0z.fits                  | PIAA M0 shape (2D sag)
./piaaconfxxx/piaam1z.fits	                | PIAA M1 shape (2D sag)
./piaaconfxxx/PIAA_Mshapes.txt              | PIAA shapes (radial txt file, cols: r0, z0, r1, z1)
./piaaconfxxx/piaa0Fz.fits                  | PIAA M0 shape, Fourier components (2D file)	
./piaaconfxxx/piaa1Fz.fits                  | PIAA M1 shape, Fourier components (2D file)
./piaaconfxxx/piaa0Cmodes.fits              | idealized PIAACMC mirror 0 cosine modes (copied from ./piaaref/)
./piaaconfxxx/piaa0Fmodes.fits              | idealized PIAACMC mirror 0 Fourier modes (copied from ./piaaref/)
./piaaconfxxx/piaa1Cmodes.fits              | idealized PIAACMC mirror 1 cosine modes (copied from ./piaaref/)
./piaaconfxxx/piaa1Fmodes.fits              | idealized PIAACMC mirror 1 Fourier modes (copied from ./piaaref/)
./piaaconfxxx/piaa0Cres.fits                | idealized PIAA M0 cosine fit residual
./piaaconfxxx/piaa1Cres.fits                | idealized PIAA M1 cosine fit residual
./piaaconfxxx/piaa0Cz.fits                  | idealized PIAA M0 cosine fit sag
./piaaconfxxx/piaa1Cz.fits                  | idealized PIAA M1 cosine fit sag
./piaaconfxxx/piaa0Fz.fits                  | idealized PIAA M0 Fourier fit sag
./piaaconfxxx/piaa1Fz.fits                  | idealized PIAA M1 Fourier fit sag


#### Idealized PIAACMC reference point:

Output file	                                  | Notes
----------------------------------------------|-------------------------------------
./piaaconfxxx/piaaref/APLCmaskCtransm.txt     | idealized PIAACMC focal plane mask transmission
./piaaconfxxx/piaaref/apo2Drad.fits	          | idealized PIAACMC output apodization
./piaaconfxxx/piaaref/piaa0Cmodes.fits        | idealized PIAACMC mirror 0 cosine modes
./piaaconfxxx/piaaref/piaa0Fmodes.fits        | idealized PIAACMC mirror 0 Fourier modes
./piaaconfxxx/piaaref/piaa1Cmodes.fits        | idealized PIAACMC mirror 1 cosine modes
./piaaconfxxx/piaaref/piaa1Fmodes.fits        | idealized PIAACMC mirror 1 Fourier modes


#### Focal plane mask:

Focal plane mask design defined by :
- [s] Sectors flag (0: no sectors, 1: sectors), variable PIAACMC_FPMsectors
- [r] Resolved target flag (0: point source, 1: resolved source)
- [mr] Mask radius in units of 0.1 l/D
- [rrr] number of rings, variable piaacmc[0].NBrings
- [zzz] number of zones

Output file                              | Notes
-----------------------------------------|-------------------------------------
fpmzmap[s]_[rrr]_[zzz].fits              | Zones map
fpm_zonea[r][s]_[mr]_[rrr]_[zzz].fits    | amplitude for each zone
fpm_zonea[r][s]_[mr]_[rrr]_[zzz].fits    | thickness for each zone


IDEALIZED OR PHYSICAL MASK\n

Idealized mask is a single zone mask with thickness adjusted for lambda/2 phase shift and a (non-physical) partial transmission.
Physical mask consist of 1 or more zones with full transmission. Each zone can have a different thickness.

By default, a non-physical mask is first created with transmission piaacmc[0].fpmaskamptransm read from piaacmcparams.conf.
Computations indices using a physical mask:
- set transmission to 1.0:  piaacmc[0].fpmaskamptransm = 1.0.
- set focal plane mask radius to larger value: piaacmc[0].fpmRad = 0.5*(LAMBDASTART+LAMBDAEND)*piaacmc[0].Fratio * PIAACMC_MASKRADLD


#### Lyot stops:

Output file	                         | Notes
-------------------------------------|-------------------------------------
./piaaconfxxx/LyotStop0.fits         | Lyot Stop 0
./piaaconfxxx/LyotStop1.fits         | Lyot Stop 1


#### Amplitude & Phase at planes:

Files are /piaaconfxxx/WFamp_nnn.fits and WFpha_nnn.fits, where nnn is the plane index.\n
Complex amplitude is shown AFTER the element has been applied, in the plane of the element.\n

Plane index                          | description
-------------------------------------|-------------------------------------
000	                                 | Input pupil
001	                                 | Fold mirror used to induce pointing offsets
002	                                 | PIAA M0
003	                                 | PIAA M1
004	                                 | PIAAM1 edge opaque mask
005	                                 | post-focal plane mask pupil
006	                                 | Lyot Stop 0
007	                                 | Lyot Stop 1
008	                                 | invPIAA1
009	                                 | invPIAA0
010	                                 | back end mask


#### Performance Evaluation:

Plane index                          | description
-------------------------------------|-------------------------------------
./piaaconfxxx/scoringmask0.fits      | Evaluation points in focal plane, hardcoded in PIAACMCsimul_computePSF()
./piaaconfxxx/CnormFactor.txt        | PSF normalization factor used to compute contrast
./piaaconfxxx/flux.txt               | total intensity at each plane









## Low-level C code

### Key functions

Every time a PSF is computed, the following 3 functions are created in that order:

- `PIAAsimul_initpiaacmcconf`
	- initialize piaacmc
	- create modes for aspheric optical surfaces description
	- CREATE DMs (Cmodes and Fmodes)
	- IMPORT / CREATE PIAA SHAPES
	- if does not exist: Creating 2D apodization for idealized circular monochromatic PIAACMC
	- compute radial PIAA sag, make 2D sag maps
	- MAKE FOCAL PLANE MASK
	- MAKE LYOT STOPS (LOADING/CREATING LYOT MASK)
	
- `PIAACMCsimul_init(piaacmc, 0, xld, yld)` : initializes optical system to piaacmc design

- `PIAACMCsimul_makePIAAshapes(piaacmc, 0)` : create 2D PIAA shapes (`piaam0z` and `piaam1z`) from coefficient values and modes

- `OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcconfdir, 0)` : perform optical propagation



Note that the last 3 functions are wrapped togeter in the computePSF function, which allows multi-point (extended sources) :

- `PIAACMCsimul_computePSF(float xld, float yld, long startelem, long endelem, int savepsf, int sourcesize, int extmode, int outsave)` : compute PSF
	- `PIAACMCsimul_init(piaacmc, 0, xld+rad1, yld)`
	- `PIAACMCsimul_makePIAAshapes(piaacmc, 0)`
	- `OptSystProp_run(optsyst, 0, startelem, optsyst[0].NBelem, piaacmcconfdir, 0)`


### Full PIAACMC design in monochromatic light

The first step of the optimization process is to compute the PIAACMC optics and focal plane mask in the idealized case: monochromatic light, point source, and ideal focal plane mask (circular, phase-shifting and partially transmissive). Steps are shown in the figure below and correspond to steps in the `runPIAACMC` state machine.


@dot
digraph runPIAACMCsteps {
  graph [fontsize=12 labelloc="t" label="" splines=true overlap=false];
  size ="12,12";
  ratio = auto;
  
  "step000" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white"> STEP 000 </font></td></tr>
  <tr><td align="left"> piaacmcsimrun code 000</td></tr>
  <tr><td align="left"> make monochromatic PIAACMC idealized design</td></tr>
  </table>> ];
 
  "step001" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white"> STEP 001 </font></td></tr>
  <tr><td align="left"> piaacmcsimrun code 000</td></tr>
  <tr><td align="left"> repeat above step to compute on-axis PSF</td></tr>
  </table>> ];

  "step002" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white"> STEP 002 </font></td></tr>
  <tr><td align="left"> - </td></tr>
  <tr><td align="left"> Specify input pupil geometry </td></tr>
  </table>> ];

  "step003" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white"> STEP 003 </font></td></tr>
  <tr><td align="left"> piaacmcsimrun code 000</td></tr>
  <tr><td align="left"> compute on-axis PSF </td></tr>
  </table>> ];

  "step004" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white"> STEP 004 </font></td></tr>
  <tr><td align="left"> piaacmcsimrun code 005</td></tr>
  <tr><td align="left"> compute Lyot stops shapes and locations</td></tr>
  <tr><td align="left"> 2nd pass, LStransm2 throughput</td></tr>
  </table>> ];

  "step005" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white"> STEP 005 </font></td></tr>
  <tr><td align="left"> piaacmcsimrun code 002</td></tr>
  <tr><td align="left"> optimize focal plane mask transm</td></tr>
  </table>> ];

  "step006" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white"> STEP 006 </font></td></tr>
  <tr><td align="left"> piaacmcsimrun code 005</td></tr>
  <tr><td align="left"> compute Lyot stops shapes and locations</td></tr>
  <tr><td align="left"> 2nd pass, LStransm0 throughput</td></tr>
  </table>> ];

  "step007" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white"> STEP 007 </font></td></tr>
  <tr><td align="left"> piaacmcsimrun code 040</td></tr>
  <tr><td align="left"> tune PIAA shapes and focal plane mask transm</td></tr>
  <tr><td align="left"> 10 cosine modes, 5 Fourier modes</td></tr>
  </table>> ];

  "step008" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white"> STEP 008 </font></td></tr>
  <tr><td align="left"> piaacmcsimrun code 040</td></tr>
  <tr><td align="left"> tune PIAA shapes and focal plane mask transm</td></tr>
  <tr><td align="left"> 20 cosine modes, 10 Fourier modes</td></tr>
  </table>> ];

  "step009" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white"> STEP 009 </font></td></tr>
  <tr><td align="left"> piaacmcsimrun code 005 </td></tr>
  <tr><td align="left"> compute Lyot stops shapes and locations</td></tr>
  <tr><td align="left"> 3nd pass, LStransm1 throughput</td></tr>
  </table>> ];

  "step010" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white"> STEP 010 </font></td></tr>
  <tr><td align="left"> piaacmcsimrun code 001 </td></tr>
  <tr><td align="left"> Tune Lyot stops conjugations</td></tr>
  </table>> ];

  "step011" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white"> STEP 011 </font></td></tr>
  <tr><td align="left"> piaacmcsimrun code 040</td></tr>
  <tr><td align="left"> tune PIAA shapes and focal plane mask transm</td></tr>
    <tr><td align="left"> 20 cosine modes, 20 Fourier modes</td></tr>
  </table>> ];

  "step012" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white"> STEP 012 </font></td></tr>
  <tr><td align="left"> piaacmcsimrun code 040</td></tr>
  <tr><td align="left"> tune PIAA shapes and focal plane mask transm</td></tr>   
  <tr><td align="left"> 40 cosine modes, 150 Fourier modes</td></tr>
  </table>> ];

  "step013" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white"> STEP 013 </font></td></tr>
  <tr><td align="left"> piaacmcsimrun code 005 </td></tr>
  <tr><td align="left"> compute Lyot stops shapes and locations</td></tr>
  <tr><td align="left"> 4rth pass, LStransm2 throughput</td></tr>
  </table>> ];

  "step014" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white"> STEP 014 </font></td></tr>
  <tr><td align="left"> piaacmcsimrun code 001 </td></tr>
  <tr><td align="left"> Tune Lyot stops conjugations</td></tr>
  </table>> ];

  "step015" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white"> STEP 015 </font></td></tr>
  <tr><td align="left"> piaacmcsimrun code 040</td></tr>
  <tr><td align="left"> tune PIAA shapes and focal plane mask transm</td></tr>
  <tr><td align="left"> 20 cosine modes, 20 Fourier modes</td></tr>
  </table>> ];

  "step016" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white"> STEP 016 </font></td></tr>
  <tr><td align="left"> piaacmcsimrun code 040</td></tr>
  <tr><td align="left"> tune PIAA shapes and focal plane mask transm</td></tr>
  <tr><td align="left"> 40 cosine modes, 150 Fourier modes</td></tr>
  </table>> ];

  "step017" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white"> STEP 017 </font></td></tr>
  <tr><td align="left"> piaacmcsimrun code 040</td></tr>
  <tr><td align="left"> tune PIAA shapes and focal plane mask transm</td></tr>
  <tr><td align="left"> 40 cosine modes, 625 Fourier modes</td></tr>
  </table>> ];


  
  step000 -> step001 [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
  step001 -> step002 [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
  step002 -> step003 [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
  step003 -> step004 [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
  step004 -> step005 [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
  step005 -> step006 [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
  step006 -> step007 [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
  step007 -> step008 [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
  step008 -> step009 [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
  step009 -> step010 [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
  step010 -> step011 [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
  step011 -> step012 [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
  step012 -> step013 [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
  step013 -> step014 [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
  step014 -> step015 [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
  step015 -> step016 [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
  step016 -> step017 [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];

  {rank=same step000 step006 step012}
}


@enddot





### APLC design

The scripts produce an APLC design by setting `PIAAmode = 0`. In that case, only `step000` in the figure above is executed.


### Focal plane mask optimization

#### Overall flow

FPM optimization is done as `mode = 13` in the main routine. It is called from the `runPIAACMC` script.  
The `PIAACMCsimul_run(char *confindex, long mode)` function loop-calls the `PIAACMCsimul_exec(char *confindex, long mode)` function with mode 13 until search time is reached.

Each search proceeds as follows :

- a starting point is picked. A random integer `zeroST` encodes the starting poing:
	- 0: 
	- 1: start at zero
	- 2: start at best solution
	
- run `PIAACMCsimul_exec` in mode 13
	- randomly move next to starting point
	- compute value and store it to valref variable
	- enter gradient search loop
		- compute local derivatives
		- compute SVD of Jacobian matrix to identify vectors of steepest descent
		- for each alphareg = 0, 0.25, 0.5, 0.75, 1.0:
			- do a linescan descent along the vector identified by SVD, with the alphareg regularization coefficient
			- the best value is stored as `bestval` and the corresponding parameters into `data.image[IDoptvec].array.F[i]`
			- store best values in `PIAACMCSIMUL_VAL`, and reference value in `PIAACMCSIMUL_VALREF`
		

- `computePSF_FAST_FPMresp = 1` : uses pre-computed complex amplitude zones response to evaluate FPM solutions
- `PIAACMC_FPM_FASTDERIVATIVES = 1`


To stop the optimization and have the current best solution adopted, enter a small value in the `searchtime.txt` file, so that the optimization process completes.


#### Regularization

There are two regularization:

- Focal plane mask zone sags can be included in the linear optimization as extra components to the optimization vector
- PIAA shapes coefficients are simply added to the evaluation metric


#### Key files

- `fpm_zonez_.....best.fits` is the optimal solution
- `mode13_....bestval.txt` is the corresponding value

To restart the optimization from scratch, remove these two files.

Average contrast = 5.964e-05
0.001 -> Average contrast = 8.57575e-07, 1.9e-7
0.01  -> Average contrast = 7.99085e-07, 2.0e-7
7.8294e-07




