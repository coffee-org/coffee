---
name: diagnose-build-failure
description: Troubleshooting CMake and compilation failures in coffee.
---

# Diagnose Build Failures

A guide for troubleshooting build and link errors in `coffee`.

## 1. CMake Configuration Failures
- **Missing packages:** Look for pkg-config or CMake errors:
  - "fftw3 or fftw3f not found" -> Install `libfftw3-dev`.
  - "GSL not found" -> Install `libgsl-dev`.
  - "OpenMP not found" -> Install `libomp-dev`.
  - "cfitsio not found" -> Install `libcfitsio-dev`.
  - "CLIcore not found" -> Ensure `milk` is installed, or set `MILK_SOURCE_DIR` / `MILK_ROOT`.

## 2. Compiler Errors
- **Undefined references / symbols:** Check if headers are missing or if the source file was not
  added to the corresponding module's `CMakeLists.txt` (`OptSystProp`, `AOsystSim`, `coronagraphs`,
  `PIAACMCsimul`).
- **Implicit declaration warnings:** Occur when a function is called without a matching header
  include. Fix by explicitly adding the `#include` of the header defining the function.
- **Legacy OpticsMaterials / milk macros:** Ensure `coffee_compat.h` is included for legacy shims.
- **Size ratchet failures:** If compilation passes but `./scripts/check_code_size.sh` fails,
  decompose oversized functions (<60 lines) or files (<600 lines).

## 3. Linker Errors
- **Missing library linkage:** Ensure that the target library and executables link against
  `${FFTW_LIBRARIES}`, `${FFTWF_LIBRARIES}`, `GSL::gsl`, `OpenMP::OpenMP_C`, and `m`.
