# Design Steps: Polychromatic Design {#page_PIAACMCsimul_DesignSteps_Poly}

---

## 1. Polychromatic FPM Zone Response Calibration (`step100` / Mode `11`)

```bash
coffee-piaacmc run --poly --jobs 4 --search-time 600
# or via legacy runner:
COFFEE_JOBS=4 COFFEE_SEARCHTIME=600 coffee-runPIAACMCdesign -m 200
```

In `step100` (`exec_multizone_fpm_calib.c`), the focal plane mask is partitioned into concentric
hexagonal rings (`NBrings`) and angular sectors (`FPMsectors`). The complex electric-field response
of every mask zone across all `nblambda` wavelengths on the focal-plane scoring pixels is computed:

- `NBthreads` (`--jobs` / `COFFEE_JOBS`) background worker processes are launched in parallel
  without requiring `tmux`.
- Each worker logs progress to `log/step100_thread<idx>.log` and writes partial zone responses to
  `FPMresp..._thread<idx>.fits.tmp`, atomically renaming to `.fits` upon completion.
- The coordinator process merges all partial zone responses into `FPMresp1.*.fits` and subtracts
  zone contributions from the unocculted field to isolate light passing outside the physical mask.
- If `FPMresp1.*.fits` already exists on disk, `step100` is automatically skipped on subsequent
  optimization passes.

---

## 2. Multi-Zone Sag Optimization (`step101` / Mode `13`)

Once `FPMresp1.*.fits` is available, `exec_optimize_fpm_zones.c` evaluates polychromatic focal-plane
electric fields and exact analytical derivatives with respect to zone sags in milliseconds:

1. **Outer Loop (Multi-Start Perturbation)**: Selects a starting zone-sag vector (zero sag, random
   perturbation, or neighborhood of the current `.best.fits` solution).
2. **Inner Loop (Analytical Jacobian + Multi-Alpha SVD Line Search)**:
   - Evaluates analytical derivatives around the current sag vector via
     `PIAACMCsimul_achromFPMsol_eval_zonezderivative`.
   - Computes an SVD pseudo-inverse search direction across multiple regularization weights
     ($\alpha_{\text{reg}} \in \{0, 0.25, 0.5, 0.75, 1.0\}$) including explicit FPM sag
     regularization (`fpmregsag_coeff`, `fpmregsag_alpha`).
   - Performs a line-scan along each search direction and records improved solutions into
     `fpm_zonez.*.best.fits` and `mode13.*.bestval.txt`.
   - Repeats until convergence (or 50 inner steps), then draws a new outer-loop starting point until
     `searchtime` (`--search-time` / `COFFEE_SEARCHTIME`) expires.

To stop Mode `13` cleanly after the active inner loop finishes:

```bash
touch piaacmcconf_i000/stoploop13.txt
```
