---
trigger: always_on
---

# Architecture Principles

- **Minimize Cross-Module Dependencies**: Avoid introducing circular or unstructured dependencies
  between modules (`OptSystProp`, `coronagraphs`, `PIAACMCsimul`, `AOsystSim`). Refactor towards
  a layered architecture where modules interact only through well-defined public APIs declared
  in their top-level headers.
- **Pure Compute vs. CLI Wrappers**: Keep simulation engines (Fresnel propagation, prolate
  spheroidal wave functions, PIAACMC optimization, AO filtering) pure and decoupled from
  interactive CLI terminal commands and user prompts.
- **Adhere to Module Boundaries**: Each submodule directory contains its own code and documentation
  defining its role and public headers. Respect these boundaries when adding new features:
  - `OptSystProp`: General optical propagation engine and wavefront cubes.
  - `coronagraphs`: Coronagraph design mathematics (APLC, prolate masks).
  - `PIAACMCsimul`: PIAA mirror shapes, complex focal plane masks, Lyot stops.
  - `AOsystSim`: Adaptive optics system simulation and WFS modeling.
