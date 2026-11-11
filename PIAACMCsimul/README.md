# Module PIAACMCsimul {#page_module_PIAACMCsimul}

Diffraction-based simulation, multi-step optimization, and performance evaluation of
Phase-Induced Amplitude Apodization Complex Mask Coronagraphs (**PIAACMC**) and Apodized Pupil
Lyot Complex Mask Coronagraphs (**APLCMC**).

---

## 1. Design Pipeline Architecture

A complete coronagraph design proceeds in three stages inside a dedicated working directory:

1. **Monochromatic Idealized Design (Steps `000`–`017`)**:
   - **Step `000`–`001` (Mode `0`)**: Computes a 2D prolate spheroidal apodization for a centrally
     obscured circular aperture, solves radial and 2D PIAA mirror/lens sags (`piaam0z.fits`,
     `piaam1z.fits`), fits radial Cosine (`Cmodes`) and 2D Fourier (`Fmodes`) bases, and propagates
     monochromatic on-axis light through an idealized single-zone phase-shifting FPM.
     *(In APLCMC mode, `PIAAmode = 0`, steps `001`–`017` are automatically skipped.)*
   - **Step `002`–`003`**: Imports the true telescope pupil (`pup_<size>.fits` ->
     `pupa0_<size>.fits`, including spiders and segment gaps) and evaluates the initial on-axis PSF.
   - **Step `004`–`017` (Modes `5`, `2`, `1`, `40`)**: Alternates 3D off-axis incoherent Lyot stop
     geometry/conjugation synthesis (Mode `5` at progressive throughputs `LStransm0 -> LStransm1 ->
     LStransm2`), Lyot axial conjugation refinement (Mode `1`), FPM transmission scan (Mode `2`),
     and regularized SVD linear optimization of PIAA Cosine + Fourier modes (Mode `40`, up to 40
     radial cosine modes and 625 Fourier modes).
2. **Polychromatic Multi-Zone FPM Calibration & Optimization (`--poly` / `-m 200`)**:
   - **Calibration (`step699` / Mode `3`)**: Propagates without FPM to record peak flux
     normalization in `CnormFactor.txt`.
   - **Zone Response Matrix (`step100` / Mode `11`)**: Partitioned across `NBthreads` parallel
     background worker processes (logged in `log/step100_thread*.log`), computes the complex
     electric-field response of every hexagonal FPM zone across all `nblambda` spectral channels on
     the focal-plane scoring pixels, and atomically merges partial files into `FPMresp1.*.fits`.
   - **Multi-Zone Sag Optimization (`step101` / Mode `13`)**: Uses `FPMresp` with analytical
     zone-sag derivatives (`PIAACMCsimul_achromFPMsol_eval_zonezderivative`) and multi-alpha SVD
     line searches from randomized perturbations around the best solution (`*.best.fits`,
     `mode13.*.bestval.txt`) within the configured `searchtime` budget.
3. **Polychromatic Evaluation & Lithographic Quantization**:
   - **PSF & Contrast (`eval` / `-m 300`, `400`, `401`, `402`)**: Full Fresnel polychromatic
     propagation for point sources, resolved stellar disks (`ssize`), off-axis throughput curves
     (Mode `101`), and pointing jitter sensitivity (`eval0/`, `eval1/`, `eval2/`).
   - **Quantization (`quantize` / `-q <NBexp>`)**: Discretizes continuous FPM zone sags into
     $2^{\text{NBexp}}$ lithographic etch levels via `coffeepiaacmcsim.fpmprocess`.

---

## 2. Unified `coffee-piaacmc` CLI

`coffee-piaacmc` provides a single entry point backed by a human-readable `piaacmc.conf` file that
stays bidirectionally synchronized with `./conf/conf_<key>.txt`.

```bash
# Initialize from a built-in example (0..7) or SPEED preset
coffee-piaacmc init --example 1
coffee-piaacmc init --example 6          # Clear-aperture APLCMC example
coffee-piaacmc init --speed1024          # SPEED 1024 PIAACMC preset

# View, edit, or validate configuration
coffee-piaacmc config show
coffee-piaacmc config set lambda=600 mlambdaB=10 nblambda=8 NBthreads=4
coffee-piaacmc config check

# Run monochromatic (steps 000..017) or full polychromatic optimization
coffee-piaacmc run --to-step 5
coffee-piaacmc run --mono
coffee-piaacmc run --poly --jobs 4 --search-time 600

# Evaluate performance, quantize FPM sags, inspect status, or clean outputs
coffee-piaacmc eval --psf
coffee-piaacmc eval --level 0
coffee-piaacmc quantize 4
coffee-piaacmc status
coffee-piaacmc clean                     # Keep pupil + piaacmc.conf, clean outputs
coffee-piaacmc clean --all               # Full workspace reset
```

### Legacy Script Compatibility

All existing workflow scripts remain fully supported and interoperate with `coffee-piaacmc`:
- `coffee-runPIAACMCdesign -e <0..7> -l` — Initialize example and list parameters
- `coffee-runPIAACMCdesign -a` / `-p` — Switch between APLCMC and PIAACMC modes
- `coffee-runPIAACMCdesign -m <1..18|200|300|400|401|402>` — Run design or evaluation step
- `coffee-runPIAACMCdesign -q <NBexp>` — Quantize `.best.fits` FPM sag solution
- `coffee-runPIAACMCdesign -c` / `-C` — Clean outputs (`-c`) or full directory (`-C`)

Environment overrides recognized by both CLIs:
- `COFFEE_JOBS=<N>`: Number of parallel worker processes for `step100` (Mode `11`).
- `COFFEE_SEARCHTIME=<sec>`: Wall-clock search budget in seconds for Mode `13`.
- `COFFEE_NBPROPSTEP=<N>`: Number of Fresnel planes in Mode `5` Lyot conjugation search.
- `COFFEE_NBITER_SA=<N>` / `COFFEE_NBITER_DS=<N>`: Mode `13` iteration overrides.

---

## 3. `piaacmc.conf` Parameter Reference

| Key | Unit | Default | Description |
|---|---|---|---|
| `PIAAmode` | flag | `1` | `1`: PIAACMC (with PIAA optics), `0`: APLCMC (no PIAA) |
| `size` | pix | `1024` | 2D complex wavefront grid dimension (`pup_<size>.fits`) |
| `beamrad` | m | `0.01` | Physical pupil beam radius in simulation space |
| `pscale` | m/pix | `5e-5` | Pupil pixel scale (must satisfy `2 * beamrad / pscale < size`) |
| `Fratio` | — | `80.0` | Focal ratio at focal plane mask |
| `lambda` | nm | `550` | Central design wavelength |
| `mlambdaB` | % | `10` | Polychromatic spectral bandwidth ($\Delta\lambda / \lambda_0$) |
| `nblambda` | count | `10` | Number of discrete wavelengths across spectral band |
| `coin` | $R_{\text{beam}}$ | `0.30` | Input pupil central obscuration fraction |
| `coout` | $R_{\text{beam}}$ | `0.20` | Remapped output central obscuration fraction |
| `pupoutmaskrad` | $R_{\text{beam}}$ | `0.95` | Output pupil mask radius fraction |
| `PIAAmaterial` | name | `Mirror` | PIAA optics material (`Mirror`, `SiO2`, `CaF2`) |
| `invPIAAmode` | flag | `0` | `0`: no inverse PIAA, `1`: inverse PIAA after Lyot stop |
| `PIAAcirc` | flag | `0` | `0`: Cosine + 2D Fourier modes, `1`: circular Cosine modes |
| `PIAA0pos` | m | `0.0` | Axial conjugation of PIAA element 0 from pupil |
| `PIAAsep` | m | `1.00` | Axial separation between PIAA elements 0 and 1 |
| `PIAAr0lim` | $R_{\text{beam}}$ | `1.15` | Outer radius limit of PIAA element 0 |
| `PIAAr1lim` | $R_{\text{beam}}$ | `1.50` | Outer radius limit of PIAA element 1 |
| `piaaNBCmodesmax` | modes | `40` | Maximum radial cosine modes for PIAA sag fitting |
| `piaaCPAmax` | CPA | `10.0` | Maximum spatial frequency for 2D Fourier PIAA modes |
| `NBls` | count | `3` | Number of Lyot stop planes |
| `LStransm` | 0..1 | `0.85` | Target geometric Lyot stop throughput |
| `LyotZmin` / `LyotZmax` | m | `-1.5`/`1.5` | Axial conjugation search range for Lyot stops |
| `fpmmaterial` | name | `Mirror` | FPM material (`Mirror`, `SiO2`, `CaF2`, `CstPha`) |
| `fpmradld` | $\lambda/D$ | `1.00` | Idealized monochromatic prolate FPM radius |
| `maskradld` | $\lambda/D$ | `2.00` | Physical multi-zone FPM outer radius |
| `NBrings` | rings | `16` | Number of concentric hexagonal rings in physical FPM |
| `FPMsectors` | mode | `2` | Zone partitioning (`0`: rings, `1`: hex, `2`: 60-deg sectors) |
| `fpmminsag` / `fpmmaxsag` | m | $\pm 10^{-6}$ | Minimum and maximum physical FPM zone sag bounds |
| `fpmregsag_coeff` | weight | `1.0` | FPM zone sag regularization weight |
| `fpmregsag_alpha` | exp | `1.0` | FPM zone sag regularization exponent |
| `ssize` | 0.1 mas | `10` | Stellar angular diameter (`10` = `1.0 mas`) |
| `extmode` | mode | `0` | Extended stellar disk model (`0`: 3-point, `1`: 6-point) |
| `NBthreads` | procs | `4` | Parallel worker processes for Mode `11` FPM calibration |
| `searchtime` | s | `600` | Wall-clock time per Mode `13` optimization pass |
| `NBpropStep` | planes | `200` | Axial Fresnel planes in Mode `5` Lyot search |

---

## 4. C Engine CLI Commands & Execution Modes

Inside `coffee` (`milk-cli`), `libcoffeePIAACMCsimul.so` registers commands under
`coffeepiaacmcsim.`:

| CLI Command | Arguments | Description |
|---|---|---|
| `coffeepiaacmcsim.run` | `<confindex> <mode>` | Primary simulation & optimization driver |
| `coffeepiaacmcsim.fpmprocess` | `<in> <zmap> <nbexp> <out>` | Quantize FPM sags (`2^nbexp`) |
| `coffeepiaacmcsim.ring2sect` | `<zin> <zmap> <fout>` | Expand ring FPM design into sectors |
| `coffeepiaacmcsim.fpmresprs` | `<in> <out> <nbl>` | Resample FPM response along wavelength |
| `coffeepiaacmcsim.fpmresprm` | `<in> <out> <nbz>` | Trim outer rings from FPM response |
| `coffeepiaacmcsim.piaacmcgeomprop` | `<surf> <sag2opd> <z>` | Geometric ray propagation |

### `coffeepiaacmcsim.run <confdir> <mode>` Modes

| Mode | Implementation | Operation |
|---|---|---|
| `0` | `exec_compute_image.c` | Initialize design (if new) and compute PSF |
| `1` | `exec_optimize_lyot_stop_position.c` | Fine scan of Lyot stop axial conjugations |
| `2` | `exec_optimize_fpmtransmission.c` | Scan idealized monochromatic FPM transmission |
| `3` | `exec_computePSF_no_fpm.c` | Unocculted PSF without FPM (`CnormFactor.txt`) |
| `4` | `exec_optimize_PIAA_shapes.c` | Linear optimization of radial cosine PIAA modes |
| `5` | `exec_optimize_lyot_stops_shapes_positions.c` | Synthesize Lyot stop masks and positions |
| `11` | `exec_multizone_fpm_calib.c` | Compute polychromatic zone response (`FPMresp`) |
| `13` | `exec_optimize_fpm_zones.c` | SVD + analytical-derivative FPM zone sag search |
| `40` | `exec_optimize_PIAA_shapes_fpmtransm.c` | Joint optimization of PIAA modes + FPM transm |
| `100` | `PIAACMCsimul_eval_poly_design.c` | Polychromatic contrast & pointing evaluation |
| `101` | `PIAACMCsimul_measure_transm_curve.c` | Off-axis planet throughput vs. separation |
| `302` | `PIAACMCsimul_exec.c` | Export physical FPM OPD/sag maps (`fpmOPD.fits`) |

---

## 5. Key Output Files (`piaacmcconf_i000/`)

- `piaacmcparams.conf` / `piaacmcparams_idealmonomask.conf`: Active and reference optical state.
- `conjugations.txt` / `lambdalist.txt`: Optical plane conjugations `[m]` and channel wavelengths.
- `apo2Drad.fits`, `PIAA_Mshapes.txt`, `piaam0z.fits`, `piaam1z.fits`: Prolate apodization and 2D
  PIAA surface sags `[m]`.
- `piaa[0-1][CF]modes.fits`: Modal coefficients on `Cmodes_<size>.fits` and `Fmodes_<size>.fits`.
- `LyotStop0.fits .. LyotStop<N-1>.fits`: Optimized Lyot stop binary transmission masks.
- `fpmzmap*.fits`, `fpm_zonez.*.best.fits`, `fpm_zonea*.fits`: Multi-zone FPM geometry, optimal
  zone sags `[m]`, and zone amplitudes.
- `FPMresp1.*.fits`: Pre-computed polychromatic complex field response cube `(2*Neval x (Nzones+1)
  x Nlambda)`.
- `CnormFactor.txt`, `psfi0.fits`, `eval0/Contrast*.txt`: Unocculted normalization, polychromatic
  PSF cube, and contrast curves.

---

## 6. Source Files & Submodules

Top-level files in `PIAACMCsimul/` (see submodule tables in
[FocalPlaneMask/README.md](FocalPlaneMask/README.md), [LyotStop/README.md](LyotStop/README.md), and
[PIAAshape/README.md](PIAAshape/README.md)):

| File | Description |
|------|-------------|
| `PIAACMCsimul.h` / `.c` | Core `PIAACMCSIMUL_PARAMS`, `OPTPIAACMCDESIGN` structs & module table |
| `PIAACMCsimul_run.h` / `.c` | `coffeepiaacmcsim.run` CLI wrapper, `mc_var` setup & Mode 13 loop |
| `PIAACMCsimul_exec.h` / `.c` | Mode `0..303` dispatcher & regularized linear perturbation loop |
| `init_piaacmcopticaldesign.h` / `.c` | Load/create `piaacmcparams.conf`, modes, PIAA, FPM & Lyot |
| `init_piaacmcopticalsystem.h` / `.c` | Populate `OptSystProp` `OPTSYST` element sequence |
| `PIAACMCsimul_loadsavepiaacmcconf.h` / `.c` | Read/write `piaacmcparams.conf` & filename tags |
| `PIAACMCsimul_computePSF.h` / `.c` | Point-source & extended-star polychromatic PSF propagation |
| `exec_compute_image.h` / `.c` | On-axis PSF or multi-point scene synthesis (`Mode 0`) |
| `exec_computePSF_no_fpm.h` / `.c` | Unocculted reference PSF calibration (`Mode 3`) |
| `exec_optimize_PIAA_shapes_fpmtransm.h` / `.c` | Joint PIAA + FPM linear setup (`Modes 40/41`) |
| `PIAACMCsimul_CA2propCubeInt.h` / `.c` | Axial Fresnel intensity cube builder for Lyot design |
| `PIAACMCsimul_achromFPMsol_eval.h` / `.c` | Fast dark-hole contrast evaluation from `FPMresp` |
| `PIAACMCsimul_achromFPMsol_eval_zonezderivative.h` / `.c` | Analytical zone-sag field derivative |
| `PIAACMC_f_evalmask.h` / `.c` | GSL Nelder-Mead / simulated-annealing merit callback |
| `PIAACMCsimul_eval_poly_design.h` / `.c` | Polychromatic contrast & jitter eval (`Mode 100`) |
| `PIAACMCsimul_measure_transm_curve.h` / `.c` | Off-axis planet throughput curve (`Mode 101`) |
| `PIAACMCsimul_free.c` | Release dynamically allocated optimization working buffers |
| `scripts/coffee-piaacmc` | Unified single-command CLI backed by `piaacmc.conf` |
| `scripts/coffee-runPIAACMC*` | Underlying step-by-step monochromatic & polychromatic drivers |
