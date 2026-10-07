# coffee Developer & Agent Guidelines (AGENTS.md)

Welcome to the `coffee` (`coffee-org/coffee`) codebase. This document outlines
core instructions, architecture boundaries, and tools for AI coding assistants and developers.

---

## 1. Project Overview & Architecture Principles

`coffee` (Coronagraph Optimization For Fast Exoplanets Exploration) is a high-performance
optical propagation and coronagraph design plugin suite for the `milk` framework
(`framework-dev` branch). It provides advanced coronagraph simulation (APLC, PIAACMC),
Fresnel wavefront propagation, deformable mirror and AO simulation, and real-time
shared memory image streaming.

### Modular Architecture:
- **`OptSystProp/`**: Optical system Fresnel propagation engine (`coffeeOptSystProp`),
  element-by-element beam propagation, multi-wavelength complex wavefront cubes.
- **`coronagraphs/`**: Coronagraph design algorithms (`coffeecoronagraphs`), prolate spheroidal
  functions, APLC optimization, focal plane mask and Lyot stop calculation.
- **`PIAACMCsimul/`**: Phase-Induced Amplitude Apodization Complex Mask Coronagraph suite
  (`coffeePIAACMCsimul`), multi-zone focal plane masks, Lyot stop geometry optimization,
  PIAA mirror shape generation, polychromatic PSF evaluation.
- **`AOsystSim/`**: Adaptive Optics system simulator (`coffeeAOsystSim`), filtering,
  pupil fitting, pyramid wavefront sensor simulation.
- **`coffee_compat.h`**: Backward-compatibility shims for legacy APIs and macros.
- **Standalone & In-Tree Builds**: Dual-mode CMake architecture supporting out-of-tree
  builds and direct integration into the `milk` framework plugin ecosystem.

### Key Constraints:
- **Zero Allocations in Simulation Loops**: Never call `malloc()`, `calloc()`, or `realloc()`
  inside hot wavefront propagation loops or PSF evaluation passes. Pre-allocate working arrays.
- **FFTW Plan Caching**: Never create and destroy FFTW plans inside per-frame loops. Plan once,
  reuse across iterations, and destroy during teardown.
- **Resource Management**: Always free allocated memory, release ImageStreamIO shared memory
  handles, close FITS file pointers, and follow the Linux kernel `goto cleanup` pattern.

---

## 2. Strict Code Style Rules

All C source code, header files, and markdown documentation must strictly comply with:

1. **Brace Style: Allman**:
   Opening braces must be on their own line for both functions and control flow
   (`if`, `for`, `while`, `switch`).
2. **Line Length Limit: $\le 100$ Characters**:
   Applies to all `.c`, `.h`, `.md`, and script files. Do not exceed 100 characters per line.
3. **Column-Aligned Function Prototypes**:
   Multi-line prototypes and definitions must use column-aligned parameter names per
   `.agents/rules/parameter-alignment.md`.
4. **Header Hygiene**:
   Every `.c` file must include exactly the headers it uses (no implicit includes).
5. **No Implicit Double Promotions**:
   Use single-precision float functions (`sqrtf`, `sinf`, `cosf`) and float literals (`0.5f`)
   when working with 32-bit floats.
6. **Code Size & Function Scope**:
   Keep files $\le 600$ lines (hard limit $1000$) and functions $\le 60$ lines (hard limit $150$).
   Enforced via `scripts/check_code_size.sh` and tracked in `scripts/code_size_baseline.txt`.
   See `.agents/rules/code-size-and-intent.md`.

---

## 3. Specialized Skills (`.agents/skills/`)

Activate these skills when working on specialized areas:
- `advanced-math-patterns`: FFTW3/FFTW3F, 2D phase grids, complex wavefronts, and SIMD math.
- `diagnose-build-failure`: CMake, OpenMP, GSL, FFTW, and milk linkage troubleshooting.
- `feature-planner`: Guidelines for scoping and planning new features or refactorings.
- `imagestream-internals`: ImageStreamIO shared memory layout and semaphore synchronization.
- `optimize-compute-function`: Checklists and guidelines for hot propagation and PSF loops.
- `pr-preparation`: Pull request validation checklist and disclosure note requirements.
- `refactor-c-source`: Safely split large simulation files into modular units (<600 lines).
- `simd-optimization`: Writing and benchmarking SIMD compute kernels with scalar fallbacks.

---

## 4. Build & Test Commands

```bash
# Build standalone coffee libraries in _build
make -C _build -j$(nproc)

# Build in-tree within milk-framework-dev
make -C ../milk-framework-dev/_build coffee-all -j$(nproc)

# Verify code size and function scope ratchet
./scripts/check_code_size.sh
```

---

## 5. Pull Request Standards & Pre-Merge Invariants

Before submitting or merging a Pull Request:
1. **Zero Warnings**: Code must compile cleanly with `-Wall -Wextra`.
2. **Ratchet Compliance**: `./scripts/check_code_size.sh` must pass with zero violations.
3. **Dual Build Support**: Both standalone build and in-tree milk build must succeed.
4. **Agentic Tool Disclosure**: If agentic tools were used, disclose in the PR description with:
   `Implemented by <model name>. Reviewed and signed off by O. Guyon.`
   Followed by a concise technical summary of what task the model performed.
