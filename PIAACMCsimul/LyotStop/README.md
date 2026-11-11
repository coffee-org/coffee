# PIAACMCsimul / LyotStop

Multi-plane Lyot stop mask synthesis, slope-driven geometric pupil intensity propagation, and axial
conjugation optimization for `PIAACMCsimul`.

---

## Source Files

| File | Description |
|------|-------------|
| `mkSimpleLyotStop.h` / `.c` | Generate circular annular Lyot stop mask (`rin <= r <= rout`) |
| `mkLyotMask.h` / `.c` | Synthesize binary Lyot mask from off-axis/on-axis pupil ratio maps |
| `optimizeLyotStop.h` / `.c` | Search 3D Fresnel intensity cubes for optimal multi-Lyot Z planes |
| `optimizeLyotStop_internal.h` | Internal context & helpers for multi-plane Lyot search |
| `optimizeLyotStop_eval.c` | On-axis Fresnel smoothing and per-zone plane selection |
| `exec_optimize_lyot_stop_position.h` / `.c` | Axial conjugation scan for Lyot stops (`Mode 1`) |
| `exec_optimize_lyot_stops_shapes_positions.h` / `.c` | Joint Lyot shape & Z solver (`Mode 5`) |
| `geomProp.h` / `.c` | Local-slope geometric ray propagation (`coffeepiaacmcsim.piaacmcgeomprop`) |

---

## Public Headers & Dependencies

- **Headers**: Installed under `include/PIAACMCsimul/` alongside `PIAACMCsimul.h`.
- **Permitted Dependencies**: `PIAACMCsimul` core context, `OptSystProp`, `milkWFpropagate`,
  `milkCOREMODmemory`, `milkCOREMODiofits`, `milkCOREMODarith`, `OpenMP`.
