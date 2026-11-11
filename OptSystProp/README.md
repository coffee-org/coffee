# Module OptSystProp {#page_module_OptSystProp}

Sequential multi-element, multi-wavelength Fresnel beam propagation engine
(`libcoffeeOptSystProp.so`) used by `PIAACMCsimul` and `AOsystSim`.

---

## Architecture & Capabilities

`OptSystProp` models an optical train in equivalent collimated beam space (`OPTSYST` in
`OptSystProp/OptSystProp.h`). Each element carries an axial conjugation `elemZpos` `[m]` relative
to the entrance pupil (`z = 0`), and beam propagation between non-coincident elements is executed
via `WFpropagate` (`Fresnel_propagate_cube`).

Supported optical element types (`elemtype`):
- **`1` — Static Amplitude / Apodizing Mask**: Multiplies complex amplitude across all spectral
  slices (`elemarrayindex` = image ID).
- **`3` — Reflective Aspheric Surface (`ASPHSURFM`)**: Converts 2D mirror sag `surfID` `[m]` into
  per-wavelength phase delay `4 * pi * sag / lambda`.
- **`4` — Refractive Aspheric Surface (`ASPHSURFR`)**: Converts 2D lens/plate sag `surfID` `[m]`
  into dispersive OPD using material indices (`OpticsMaterials`).
- **`5` — Focal Plane Mask (`FOCMASK`)**: Applies a polychromatic complex focal plane mask via
  oversampled zoomed 2D discrete Fourier transforms (`fft_DFTinsertFPM`), supporting both direct
  (`mode = 0`) and Babinet `1 - FPM` (`mode = 1`) formulations.

---

## Public C API

Declared in `OptSystProp/OptSystProp.h`:
- `OptSystProp_propagateCube(...)`: Propagate a `(size x size x nblambda)` complex amplitude/phase
  wavefront cube across axial distance `zprop` `[m]`.
- `OptSystProp_run(...)`: Propagate sequentially from `elemstart` to `elemend`, optionally writing
  per-plane diagnostic FITS cubes (`WFamp0_xxx.fits`, `WFpha0_xxx.fits`) and computing the final
  focal-plane PSF (`psfi0`, `psfa0`, `psfp0`).

---

## Source Files

| File | Description |
|------|-------------|
| `OptSystProp.h` | Public `OPTSYST`, `DM_SIM`, `ASPHSURFM`, `ASPHSURFR`, `FOCMASK` structs & API |
| `OptSystProp_internal.h` | Internal optical element application context & helper declarations |
| `OptSystProp.c` | `MILK_MODULE_RUNNABLES` registration and module initialization |
| `OptSystProp_elements.c` | Amplitude mask, mirror/refractive OPD cubes, and zoomed FPM DFT |
| `OptSystProp_propagateCube.c` | Multi-wavelength complex cube Fresnel propagation wrapper |
| `OptSystProp_run.c` | Sequential multi-element Fresnel propagation and final PSF orchestrator |

---

## Dependencies

- **Public Header**: `OptSystProp/OptSystProp.h`
- **Permitted Dependencies**: `CLIcore`, `milkCOREMODmemory`, `milkCOREMODiofits`,
  `milkCOREMODarith`, `milkfft`, `milkOpticsMaterials`, `milkWFpropagate`, `OpenMP`, `FFTW3`.
