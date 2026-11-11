# Module coronagraphs {#page_module_coronagraphs}

Mathematical coronagraph design algorithms (`libcoffeecoronagraphs.so`), including iterative 2D
prolate spheroidal wave function synthesis for centrally obscured apertures, Apodized Pupil Lyot
Coronagraph (APLC) mask compilation, and classical/apodized coronagraph PSF simulation.

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

## Public Header & Dependencies

- **Public Header**: `coronagraphs/coronagraphs.h`
- **Dependencies**: `CLIcore`, `milkCOREMODmemory`, `milkCOREMODiofits`, `milkCOREMODarith`,
  `milkfft`, `milkimagegen`, `milkinfo`, `milkstatistic`, `milkkdtree`, `GSL`, `FFTW3`.
