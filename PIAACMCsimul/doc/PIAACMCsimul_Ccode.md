# C Code Architecture {#page_PIAACMCsimul_Ccode}

`libcoffeePIAACMCsimul.so` implements the computational engine behind `coffeepiaacmcsim.*` in
`PIAACMCsimul/` and its three domain subdirectories (`FocalPlaneMask/`, `LyotStop/`, `PIAAshape/`).

---

## 1. Execution Modes (`PIAACMCsimul_run` / `PIAACMCsimul_exec`)

`coffeepiaacmcsim.run <confindex> <mode>` dispatches `PIAACMCsimul_run(confindex, mode)` ->
`PIAACMCsimul_exec(confindex, mode)`:

| Mode | Implementation | Operation |
|---|---|---|
| `0` | `exec_compute_image.c` | Initialize design (if new) and compute on-axis PSF / scene |
| `1` | `exec_optimize_lyot_stop_position.c` | Scan axial conjugations `z` `[m]` of Lyot stops |
| `2` | `exec_optimize_fpmtransmission.c` | Scan idealized monochromatic FPM transmission |
| `3` | `exec_computePSF_no_fpm.c` | Unocculted reference PSF (`CnormFactor.txt`) |
| `4` | `exec_optimize_PIAA_shapes.c` | Linear optimization of radial cosine PIAA modes |
| `5` | `exec_optimize_lyot_stops_shapes_positions.c` | Joint Lyot stop shape & Z optimization |
| `10` | `exec_multizone_fpm_calib.c` | Setup polychromatic physical FPM zone calibration |
| `11` | `exec_multizone_fpm_calib.c` | Compute polychromatic zone response slice (`FPMresp`) |
| `13` | `exec_optimize_fpm_zones.c` | Regularized SVD + analytical-derivative FPM sag search |
| `40` | `exec_optimize_PIAA_shapes_fpmtransm.c` | Monochromatic joint PIAA + FPM optimization |
| `41` | `exec_optimize_PIAA_shapes_fpmtransm.c` | Polychromatic joint PIAA + FPM optimization |
| `100` | `PIAACMCsimul_eval_poly_design.c` | Polychromatic contrast, jitter & stellar-size eval |
| `101` | `PIAACMCsimul_measure_transm_curve.c` | Off-axis planetary throughput vs. separation |
| `300..303` | `PIAACMCsimul_exec.c` | Export physical FPM OPD/sag FITS maps |

### Runtime Tuning Variables (`@calc`)

Optional overrides read from `mc_var` (`@calc "VAR = VAL"` in `milk-cli`):

| Variable | Description |
|---|---|
| `PIAACMC_size` | 2D wavefront grid size (`1024`, `2048`, ...) `[pix]` |
| `PIAACMC_pixscale` | Collimated pupil pixel scale `[m/pix]` |
| `PIAACMC_dftgrid` | Focal-plane DFT sub-sampling stride (`0` = full resolution) |
| `PIAACMC_centobs0` / `PIAACMC_centobs1` | Entrance / remapped central obscuration `[beamrad]` |
| `PIAACMC_nblambda` | Number of spectral channels |
| `PIAACMC_lambda` / `PIAACMC_lambdaB` | Central wavelength `[nm]` and bandwidth `[%]` |
| `PIAACMC_resolved` / `PIAACMC_extmode` | Extended stellar source flag and quadrature (`0`/`1`) |
| `PIAACMC_fpmtype` | `0`: idealized monochromatic mask, `1`: physical multi-zone mask |
| `PIAACMC_FPMsectors` | `0`: concentric rings, `1`: full hex, `2`: 60-deg sectors |
| `PIAACMC_NBrings` | Concentric hexagonal ring count in physical FPM |
| `PIAACMC_MASKRADLD` | Physical FPM outer radius `[lambda/D]` (`0` selects idealized FPM) |
| `PIAACMC_FPMresp_mp` / `PIAACMC_FPMresp_thread` | Worker count and worker index in Mode `11` |

---

## 2. Initialization & PSF Propagation Call Flow

Every simulation call executes four core C stages:

1. **`init_piaacmcopticaldesign(...)`** (`init_piaacmcopticaldesign.c`):
   - Creates or loads `piaacmcconf_i<idx>/piaacmcparams.conf` (`PIAACMCsimul_loadpiaacmcconf`).
   - Builds or loads 2D modal bases `Cmodes_<size>.fits` (radial cosines) and `Fmodes_<size>.fits`
     (2D Fourier modes).
   - When `PIAAmode == 1`, synthesizes the 2D prolate apodization (`apo2Drad.fits`), fits radial
     cosines (`load2DRadialApodization`), integrates 1D ray-mapping sags (`init_geomPIAA_rad` ->
     `PIAA_Mshapes.txt`), interpolates 2D sags (`mkPIAAMshapes_from_RadSag`), and projects onto
     `piaa[0-1][CF]modes.fits`.
   - Builds/loads FPM zone map (`mkFPM_zonemap`), zone sags/amplitudes (`fpm_zonez`, `fpm_zonea`),
     and Lyot stop masks (`LyotStop*.fits`).
2. **`init_piaacmcopticalsystem(xld, yld)`** (`init_piaacmcopticalsystem.c`):
   - Populates `piaacmcopticalsystem` (`OPTSYST`) with entrance pupil, pointing tip-tilt mirror,
     PIAA elements 0/1, opaque pupil masks, polychromatic focal plane mask (`mkFocalPlaneMask`),
     multi-plane Lyot stops, optional inverse PIAA optics, and back-end pupil stop.
3. **`makePIAAshapes()`** (`PIAAshape/makePIAAshapes.c`):
   - Reconstructs 2D PIAA surface sags (`piaam0z`, `piaam1z`) from modal coefficients.
4. **`OptSystProp_run(...)`** (`OptSystProp/OptSystProp_run.c`):
   - Propagates the multi-wavelength complex wavefront cube element-by-element via Fresnel
     diffraction and zoomed focal-plane DFTs.

`PIAACMCsimul_computePSF(...)` wraps stages 2–4 to evaluate point or extended stellar sources (3- or
6-point quadrature) and compute normalized dark-hole contrast.

---

## 3. Monochromatic & Polychromatic Optimization Flow

### Monochromatic Pipeline (`step000`–`step017`)
- **APLCMC (`PIAAmode = 0`)**: Runs `step000` only to initialize the pupil and optical prescription,
  then proceeds directly to polychromatic multi-zone FPM calibration.
- **PIAACMC (`PIAAmode = 1`)**:
  - `step000`–`001` (`Mode 0`): Idealized circular prolate + initial PIAA sags.
  - `step002`–`003`: Import true pupil (`pup_<size>.fits`) and evaluate initial PSF.
  - `step004`, `006`, `009`, `013` (`Mode 5`): Multi-plane Lyot stop shape & conjugation synthesis.
  - `step005` (`Mode 2`): Idealized FPM transmission optimization.
  - `step010`, `014` (`Mode 1`): Fine Lyot conjugation scan.
  - `step007`, `008`, `011`, `012`, `015`–`017` (`Mode 40`): Progressive regularized linear
    optimization of PIAA Cosine + Fourier modes (`10..40` Cosine, `5..625` Fourier modes).

### Polychromatic Multi-Zone FPM Optimization (`Mode 11` & `Mode 13`)
1. **Zone Response Calibration (`Mode 11`)**: Computes the complex electric field at dark-hole
   scoring pixels (`scoringmask0.fits`) for unit impulse on each FPM zone across all spectral
   slices, assembling `FPMresp1.<fnamedescr_conf>.fits`.
2. **Regularized SVD Descent (`Mode 13`)**:
   - `PIAACMCsimul_run` loops calls to `PIAACMCsimul_exec(confindex, 13)` until `searchtime.txt`
     expires, alternating clean starts and randomized perturbations around `*.best.fits`.
   - `exec_optimize_fpm_zones` evaluates analytical derivatives
     (`PIAACMCsimul_achromFPMsol_eval_zonezderivative`), computes the Jacobian SVD, and performs
     line searches across regularization weights `alphareg` in `[0, 1]`, saving improvements to
     `fpm_zonez.<fnamedescr>.best.fits` and `mode13.<fnamedescr>.bestval.txt`.

---

## 4. Key Workspace Files (`piaacmcconf_i<idx>/`)

| File Pattern | Contents |
|---|---|
| `piaacmcparams.conf` | Active optical prescription (`OPTPIAACMCDESIGN`) |
| `conjugations.txt` / `lambdalist.txt` | Element conjugations `[m]` and channel wavelengths `[m]` |
| `apo2Drad.fits` / `PIAA_Mshapes.txt` | 2D prolate apodizer and 1D radial PIAA sag table |
| `piaam0z.fits`, `piaam1z.fits` | 2D PIAA surface sag maps `[m]` |
| `piaa[0-1][CF]modes.fits` | PIAA Cosine and 2D Fourier modal coefficients |
| `LyotStop[0-9].fits` | Optimized binary Lyot stop masks |
| `fpmzmap*.fits` | High-resolution 2D integer FPM zone index map |
| `FPMresp1.*.fits` | Polychromatic FPM zone complex response matrix |
| `fpm_zonez.*.best.fits` | Optimal physical FPM zone sags `[m]` |
| `scoringmask0.fits` / `CnormFactor.txt` | Dark-hole scoring pixel mask and unocculted peak flux |
| `WFamp0_nnn.fits`, `WFpha0_nnn.fits` | Optional per-plane complex amplitude & phase cubes |
