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
| `AOsystSim_internal.h` | Module-internal state, context structs, and subsystem helper prototypes |
| `AOsystSim.c` | Shared module state and coupled real-time DM-to-WFS loop (`AOsystSim_run`) |
| `AOsystSim_filter.c` | Open-loop Fourier spatial/temporal AO filter (`AOsystSim_simpleAOfilter`) |
| `AOsystSim_extremeAO.c` | Analytical extreme-AO contrast error-budget simulator |
| `AOsystSim_pupil.c` | Telescope pupil generator and DM-to-pupil alignment fitter |
| `AOsystSim_wf_conf.c` | Atmospheric wavefront daemon configuration and pupil setup |
| `AOsystSim_wf.c` | Multi-layer Fresnel atmospheric wavefront streaming loop (`AOsystSim_mkWF`) |
| `AOsystSim_wfs_pyramid.c` | Modulated Pyramid WFS optical propagator and WFS dispatcher |
| `AOsystSim_pyrwfs.c` | Real-time shared-memory Pyramid WFS daemon (`AOsystSim_PyrWFS`) |
| `AOsystSim_dm_shape.c` | Sparse-cached DM influence convolution and DM template config writer |
| `AOsystSim_dm.c` | Real-time shared-memory deformable mirror OPD daemon (`AOsystSim_DM`) |
| `AOsystSim_lowfs_init.c` | Coronagraph & Lyot LOWFS configuration, stream attachment, and setup |
| `AOsystSim_lowfs.c` | Real-time coronagraph & LOWFS propagation loop (`AOsystSim_coroLOWFS`) |
| `AOsystSim_fpwfs_imsimul.c` | GSL objective (`f_eval`) and probed focal-plane PSF cube simulator |
| `AOsystSim_fpwfs_probes.c` | Quadrature DM speckle-nulling phase probe synthesis and test PSFs |
| `AOsystSim_fpwfs_solver.c` | Per-pixel Nelder-Mead simplex electric-field estimator |
| `AOsystSim_fpwfs_sens.c` | Monte Carlo focal-plane WFS coherence & sensitivity orchestrator |
| `AOsystSim_cli.c` | `MILK_RUNNABLE` CLI wrappers and `MILK_MODULE_RUNNABLES` registration |

---

## Public Header & Dependencies

- **Public Header**: `AOsystSim/AOsystSim.h`
- **Permitted Dependencies**: `CLIcore`, `milkCOREMODmemory`, `milkCOREMODiofits`,
  `milkCOREMODarith`, `milkinfo`, `milkimagegen`, `milkimagefilter`, `milkstatistic`, `milkfft`,
  `milkOpticsMaterials`, `milkWFpropagate`, `coffeeOptSystProp`, `GSL`, `FFTW3`, `OpenMP`.
