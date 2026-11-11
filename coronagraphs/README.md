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
| `coronagraphs.c` | Core prolate spheroidal math, APLC fitting, PIAA init, and coronagraph models |
| `coronagraphs_cli.c` | `MILK_RUNNABLE` CLI wrappers and `MILK_MODULE_RUNNABLES` registration |
| `scripts/coffee` | Top-level `milk-cli` launcher preloading all `coffee` optical modules |

---

## Public Header & Dependencies

- **Public Header**: `coronagraphs/coronagraphs.h`
- **Permitted Dependencies**: `CLIcore`, `milkCOREMODmemory`, `milkCOREMODiofits`,
  `milkCOREMODarith`, `milkfft`, `milkimagegen`, `milkinfo`, `milkstatistic`, `milkkdtree`,
  `GSL`, `FFTW3`, `OpenMP`.
