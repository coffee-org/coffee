# Module coronagraphs {#page_module_coronagraphs}

Mathematical coronagraph design algorithms (`libcoffeecoronagraphs.so`), including iterative 2D
prolate spheroidal wave function synthesis for centrally obscured apertures, Apodized Pupil Lyot
Coronagraph (APLC) mask compilation, PIAA radial ray-mapping, and multi-architecture coronagraph
PSF / throughput simulation (`AIC`, `4QPM`, `BL4`, `BL8`, `OVC`, `RRPM`, `CPA`, `PIAA`, `PIAAC`,
`MULTISTEP_APLC`).

Also installs the top-level `coffee` launcher script, which starts `milk-cli` with
`coffeeOptSystProp`, `coffeecoronagraphs`, `coffeePIAACMCsimul`, and `coffeeAOsystSim` preloaded.

---

## CLI Commands (`coro.` Namespace)

Inside `coffee` (`milk-cli`), run `@help coro` or `<cmd> -h` to inspect arguments:

| CLI Command | Arguments | Description |
|---|---|---|
| `coro.cormk2Dprolate` | `<fpmpix> <beampix> <co> <out> <size>` | 2D prolate apodizer (pix) |
| `coro.cormk2Dprolateld` | `<fpmld> <beampix> <co> <out> <size>` | 2D prolate apodizer (`l/D`) |
| `coro.corup2Dprolate` | `<fpmld> <beampix> <co> <zfactor>` | Refine 2D prolate apodizer |
| `coro.corprolcomp` | *(none)* | Compile APLC prolate tables |
| `coro.corsimpsf` | `<xld> <yld> <psfname> <corotype>` | Simulate coronagraph PSF |
| `coro.coroscanpiaacmcperf` | `<obs0>` | Scan PIAACMC contrast vs. `obs0` |

Optional numeric tuning variables (set via `@calc "VAR = VAL"` before calling prolate routines):
- `PROLATE_NBiter`: Override DFT iteration count in 2D prolate generation.
- `CORONAGRAPHS_DFTZFACTOR`: Override discrete Fourier transform zoom factor.

---

## Source Files

| File | Description |
|------|-------------|
| `coronagraphs.h` | Public C API for 2D prolate generation, PIAA remapping, and PSF simulators |
| `coronagraphs_internal.h` | Internal state (`coro_state_t`), constants, and cross-file helpers |
| `coronagraphs.c` | Shared module state definitions and common pupil, FFT, and PSF utilities |
| `coronagraphs_prolate.c` | FFT-based and circularly symmetric 2D prolate apodizer synthesis |
| `coronagraphs_prolate_dft.c` | Zoomed-DFT 2D prolate spheroidal eigenmode iteration engine |
| `coronagraphs_prolate_update.c` | Iterative zoomed-DFT refinement of existing 2D prolate masks |
| `coronagraphs_apofit_eval.c` | Analytical multi-basis radial apodizer profile evaluation |
| `coronagraphs_apofit.c` | GSL Nelder-Mead analytical fitting of 1D prolate apodizer profiles |
| `coronagraphs_aplc_compile.c` | Precomputed APLC prolate FITS lookup-cube compiler |
| `coronagraphs_piaa_init.c` | Radial PIAA apodization profile loading and ray-mapping setup |
| `coronagraphs_piaa_remap.c` | PIAA mirror shape integration, forward remapping, and unwrapping |
| `coronagraphs_piaa_beam.c` | PIAA cleanup, Subaru pupil mask generator, and Zernike sensitivity |
| `coronagraphs_masks_init.c` | Focal-plane and pupil mask builders (`CPA`, `ODC`, `BL`, `OVC`) |
| `coronagraphs_sim_phase.c` | Phase-mask and interferometric PSF simulators (`AIC`, `4QPM`, etc.) |
| `coronagraphs_sim_bandlim.c` | Band-limited and optical differentiation PSF simulators |
| `coronagraphs_sim_vortex_rrpm.c` | Roddier & Roddier phase mask, optimizer, and vortex PSFs |
| `coronagraphs_sim_apodized.c` | Apodized (`CPA`, `PPA`, `STRIPC`, `SIMXY`) and unocculted PSFs |
| `coronagraphs_sim_piaac.c` | `PIAA`, `PIAAC`, and hybrid `AIC_PIAAC` coronagraph PSF simulators |
| `coronagraphs_sim_aplc.c` | Multi-step APLC/PIAACMC setup and external occulter simulators |
| `coronagraphs_sim_aplc_prop.c` | Multi-step APLC focal-plane cascade and Lyot-plane PSF readout |
| `coronagraphs_eval.c` | Top-level `coronagraph_simulPSF` dispatcher and user hook |
| `coronagraphs_transm.c` | Off-axis planetary throughput and IWA50 curve evaluator |
| `coronagraphs_limitcoeff.c` | Fundamental Taylor-series inner-working-angle limit calculator |
| `coronagraphs_scan_compute.c` | Grid PSF and radial profile synthesis for PIAACMC design scans |
| `coronagraphs_scan.c` | PIAACMC central-obscuration performance scanner and FITS summary |
| `coronagraphs_cli.c` | `MILK_RUNNABLE` CLI wrappers and `MILK_MODULE_RUNNABLES` registration |
| `scripts/coffee` | Top-level `milk-cli` launcher preloading all `coffee` optical modules |

---

## Public Header & Dependencies

- **Public Header**: `coronagraphs/coronagraphs.h`
- **Permitted Dependencies**: `CLIcore`, `milkCOREMODmemory`, `milkCOREMODiofits`,
  `milkCOREMODarith`, `milkfft`, `milkimagegen`, `milkinfo`, `milkstatistic`, `milkkdtree`,
  `GSL`, `FFTW3`, `OpenMP`.
