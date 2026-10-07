---
description: Naming conventions for files, variables, functions, and structures.
---

# Naming Conventions

Maintain consistency with the existing `coffee` codebase naming patterns.

## 1. Files & Folders
- Use lowercase `snake_case` or module-prefixed names for C source and header files:
  - `OptSystProp/` (optical propagation and Fresnel diffraction)
  - `AOsystSim/` (adaptive optics system simulations)
  - `coronagraphs/` (coronagraphic masks, Lyot stops, PSF calculation)
  - `PIAACMCsimul/` (PIAACMC optics and focal plane mask optimization)
- Keep file names descriptive and under 600 lines.

## 2. Functions
- Public API functions declared in headers: `<Module>_<verb>_<object>`:
  - Examples: `OptSystProp_FresnelProp`, `coronagraphs_init`, `PIAACMCsimul_run`.
- Static helper functions: lowercase `snake_case` descriptive names.

## 3. Variables
- Local variables: short, lowercase `snake_case`.
- Loop indices:
  - Inner grid/pixel loops: use doubled letters `ii`, `jj`, `kk`.
  - Outer or layer loops: use descriptive names like `layer_idx`, `step_idx`, `frame_idx`.
  - Avoid single-character indices (`i`, `j`, `k`) in non-trivial loops to remain searchable.
- Dimension variables:
  - `xsize`, `ysize`: 2D image dimensions
  - `zsize`: 3D cube depth / number of frames or wavelengths
  - `npupil`: pupil diameter / dimension in pixels
- Pointers: descriptive name, optionally suffixed with `_ptr` or typed arrays (`*pha`, `*amp`).

## 4. Structs & Types
- Struct names: use `UPPER_CASE` or `CamelCase` matching `milk` conventions:
  - `OPTSYSTPROP_CONFIG`, `PIAACMC_MASK_ZONE`, etc.
