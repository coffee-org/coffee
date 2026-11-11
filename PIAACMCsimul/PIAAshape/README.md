# PIAACMCsimul / PIAAshape

Radial and 2D aspheric Phase-Induced Amplitude Apodization (PIAA) mirror/lens sag generation,
cosine apodization fitting, and modal surface optimization for `PIAACMCsimul`.

---

## Source Files

| File | Description |
|------|-------------|
| `load2DRadialApodization.h` / `.c` | Project 2D prolate apodizer onto radial cosine basis |
| `init_geomPIAA_rad.h` / `.c` | Integrate 1D geometric ray-mapping ODEs -> `PIAA_Mshapes.txt` |
| `mkPIAAMshapes_from_RadSag.h` / `.c` | Interpolate 1D radial sags onto 2D pupil grids |
| `makePIAAshapes.h` / `.c` | Reconstruct 2D PIAA surfaces from Cosine (`Cmodes`) + Fourier modes |
| `exec_optimize_PIAA_shapes.h` / `.c` | Linear optimizer for radial cosine PIAA modes (`Mode 4`) |

---

## Public Headers & Dependencies

- **Headers**: Installed under `include/PIAACMCsimul/` alongside `PIAACMCsimul.h`.
- **Permitted Dependencies**: `PIAACMCsimul` core context, `OptSystProp`, `milkOpticsMaterials`,
  `milklinoptimtools`, `milkimagegen`, `milkCOREMODmemory`, `milkCOREMODiofits`.
