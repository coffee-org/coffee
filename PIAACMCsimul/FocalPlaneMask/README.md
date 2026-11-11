# PIAACMCsimul / FocalPlaneMask

Multi-zone focal plane mask (FPM) geometry generation, polychromatic complex-amplitude synthesis,
zone impulse-response calibration, analytical-derivative SVD sag optimization, and discrete
lithographic quantization for `PIAACMCsimul`.

---

## Source Files

| File | Description |
|------|-------------|
| `mkFPM_zonemap.h` / `.c` | Build 2D hexagonal ring/sector integer zone index map (`fpmzmap`) |
| `mkFocalPlaneMask.h` / `.c` | Synthesize 3D polychromatic complex-amplitude FPM cube |
| `exec_optimize_fpmtransmission.h` / `.c` | Optimize idealized monochromatic FPM transmission |
| `exec_multizone_fpm_calib.h` / `.c` | Calibrate polychromatic FPM zone response (`Mode 11`) |
| `exec_optimize_fpm_zones.h` / `.c` | Regularized SVD + line-search FPM sag optimizer (`Mode 13`) |
| `FPMresp_resample.h` / `.c` | Interpolate/decimate `FPMresp` cube (`coffeepiaacmcsim.fpmresprs`) |
| `FPMresp_rmzones.h` / `.c` | Collapse outer rings of `FPMresp` (`coffeepiaacmcsim.fpmresprm`) |
| `rings2sectors.h` / `.c` | Map radial ring sags onto sector zones (`coffeepiaacmcsim.ring2sect`) |
| `FPM_process.h` / `.c` | Quantize continuous sags into $2^N$ etch levels (`fpmprocess`) |

---

## Public Headers & Dependencies

- **Headers**: Installed under `include/PIAACMCsimul/` alongside `PIAACMCsimul.h`.
- **Permitted Dependencies**: `PIAACMCsimul` core context, `OptSystProp`, `milkOpticsMaterials`,
  `milklinoptimtools`, `milkCOREMODmemory`, `milkCOREMODiofits`, `OpenMP`.
