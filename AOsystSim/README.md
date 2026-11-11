# Module AOsystSim {#page_module_AOsystSim}

Adaptive Optics system, Pyramid wavefront sensor, deformable mirror, Lyot low-order WFS (LOWFS),
and focal-plane wavefront sensing (FPWFS) simulator (`libcoffeeAOsystSim.so`).

---

## CLI Commands (`aosystsim.` Namespace)

Inside `coffee` (`milk-cli`), run `@help aosystsim` or `<cmd> -h` to inspect parameters:

| CLI Command | Arguments | Description |
|---|---|---|
| `aosystsim.AOsystsim` | `<syncmode> <dmindex> <delayus>` | Coupled DM-to-WFS loop simulator |
| `aosystsim.AOsimfilt` | `<inwf> <outwf>` | Open-loop Fourier spatial/temporal AO filter |
| `aosystsim.AOsystexaosim` | *(none)* | Analytical extreme-AO error-budget & contrast |
| `aosystsim.AOsystsfitpup` | `<dmresp> <telpup>` | Fit pupil geometry to measured DM response |
| `aosystsim.AOsimmkWF` | `<conffile>` | Multi-layer Fresnel atmospheric WF generator |
| `aosystsim.AOsimPyrWFS` | `<conffile>` | Modulated Pyramid WFS shared-memory daemon |
| `aosystsim.AOsimDM` | `<conffile>` | Deformable mirror surface OPD daemon |
| `aosystsim.AOsimcoroLOWFS` | `<conffile>` | Coronagraph & Lyot LOWFS shared-memory daemon |
| `aosystsim.AOsystmkABprobes` | `<prA> <prB> <nx> <ny> ...` | Generate FPWFS DM probes |
| `aosystsim.AOsystFPWFSan` | `<mapmode> <mode> <opt> <N>` | FPWFS coherence & sensitivity |

---

## Source Files

| File | Description |
|------|-------------|
| `AOsystSim.h` | Public `EXAOSIMCONF` struct and AO/WFS/DM/LOWFS/FPWFS C API declarations |
| `AOsystSim.c` | Core simulation routines (AO filter, pupil fit, WF, PyrWFS, DM, LOWFS, FPWFS) |
| `AOsystSim_cli.c` | `MILK_RUNNABLE` CLI wrappers and `MILK_MODULE_RUNNABLES` registration |

---

## Public Header & Dependencies

- **Public Header**: `AOsystSim/AOsystSim.h`
- **Permitted Dependencies**: `CLIcore`, `milkCOREMODmemory`, `milkCOREMODiofits`,
  `milkCOREMODarith`, `milkinfo`, `milkimagegen`, `milkimagefilter`, `milkstatistic`, `milkfft`,
  `milkOpticsMaterials`, `milkWFpropagate`, `coffeeOptSystProp`, `GSL`, `FFTW3`, `OpenMP`.
