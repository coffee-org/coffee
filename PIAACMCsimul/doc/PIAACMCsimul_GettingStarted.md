# Getting Started {#page_PIAACMCsimul_GettingStarted}

*Source: `PIAACMCsimul_GettingStarted.md` in `PIAACMCsimul/doc`*

---

## 1. Environment Setup

Ensure `<milk-install-prefix>/bin` is on your `PATH`:

```bash
export PATH="<milk-install-prefix>/bin:$PATH"
coffee-piaacmc --help
```

No manual script copying (`syncscripts`) or `tmux` installation is required.

---

## 2. Configuration Setup

### 2.1. Unified `coffee-piaacmc` CLI (Recommended)

Create a clean working directory and initialize one of the built-in examples (`0..7`) or `SPEED`
presets:

```bash
mkdir work_piaacmc && cd work_piaacmc
coffee-piaacmc init --example 1
```

This creates the input pupil `pup_<size>.fits`, `./conf/conf_*.txt`, and a commented single-file
configuration `piaacmc.conf`. Inspect, modify, or validate parameters at any time:

```bash
coffee-piaacmc config show
coffee-piaacmc config set mlambdaB=10 nblambda=8 NBthreads=4
coffee-piaacmc config check
```

For an APLC/APLCMC design (no PIAA optics), pass `--aplc` during `init` or set `PIAAmode=0`:

```bash
coffee-piaacmc init --example 6 --aplc
```

### 2.2. Legacy `coffee-runPIAACMCdesign` Script

All legacy flags remain supported and stay synchronized with `piaacmc.conf`:

```bash
coffee-runPIAACMCdesign -e 1 -l
coffee-runPIAACMCdesign -a    # optional: switch to APLCMC mode
```

---

## 3. Running the Design

### Monochromatic Idealized Design (Steps `000`–`017`)

```bash
coffee-piaacmc run --mono
# or stop at an intermediate step (e.g., step 5):
coffee-piaacmc run --to-step 5
```

### Polychromatic Multi-Zone FPM Design (`step100` + `step101`)

```bash
coffee-piaacmc run --poly --jobs 4 --search-time 600
```

This automatically:
1. Completes any remaining monochromatic steps (`000`–`017`).
2. Calibrates unocculted flux normalization (`CnormFactor.txt`).
3. Computes the polychromatic FPM zone response matrix (`FPMresp1.*.fits`) in parallel across
   `--jobs` background worker processes (logs in `log/step100_thread*.log`).
4. Runs analytical-derivative + SVD multi-zone sag optimization (Mode `13`) for `--search-time`
   seconds.

---

## 4. Monitoring & Controlling Optimization

Track live optimization progress with:

```bash
coffee-piaacmc status
tail -f piaacmcconf_i000/linoptval.txt
tail -f piaacmcconf_i000/mode13*.opt.txt
```

To request a clean early exit from Mode `13` after the current inner line-search completes:

```bash
touch piaacmcconf_i000/stoploop13.txt
```

---

## 5. Evaluating Performance & Quantizing Masks

```bash
# Compute polychromatic on-axis PSF cube -> psfi0.fits
coffee-piaacmc eval --psf

# Full evaluation: level 0 (contrast + off-axis transmission), level 1 (stellar size),
# or level 2 (stellar size + pointing jitter)
coffee-piaacmc eval --level 0

# Quantize optimal continuous FPM sags (.best.fits) to 2^4 = 16 lithographic levels
coffee-piaacmc quantize 4
```
