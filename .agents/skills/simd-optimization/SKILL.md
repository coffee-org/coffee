---
name: simd-optimization
description: Methodology for writing, auditing, and benchmarking SIMD compute kernels in
  coffee.
---

# SIMD Optimization Guide

`coffee` performs optical propagation, coronagraph simulations, and AO reconstruction.
Hot compute loops must saturate CPU vector execution pipelines and cache bandwidth.

This skill documents the SIMD architecture in `coffee`, identifying where explicit intrinsics
are essential for maximum throughput and how kernels are organized across ISAs.

---

## 1. High-Value SIMD Opportunities in coffee

1. **Fresnel Propagation Kernel Products (`OptSystProp`):**
   * *Mechanism:* Complex multiplication of pupil fields by quadratic Fresnel phase factors.
   * *Files:* `OptSystProp/OptSystProp_FresnelProp.c`.
2. **Coronagraph Mask Application (`coronagraphs`):**
   * *Mechanism:* Vectorized amplitude masking and phase modulation across 2D optical grids.
   * *Files:* `coronagraphs/coronagraphs.c`.
3. **Multi-Zone Complex Transmission (`PIAACMCsimul`):**
   * *Mechanism:* Complex exponential calculations $e^{i \phi(\lambda)}$ across annular zones.
   * *Files:* `PIAACMCsimul/PIAACMCsimul_init.c`.
4. **AO Matrix-Vector Multiplication (`AOsystSim`):**
   * *Mechanism:* GEMV dot products for wavefront reconstruction from modal/zonal commands.
   * *Files:* `AOsystSim/AOsystSim.c`.

---

## 2. Architecture & File Layout

Explicit SIMD in `coffee` follows the ISA Strategy Separation pattern:
- **`<kernel>_simd.h`**: Public API declarations and function prototypes.
- **`<kernel>_simd_scalar.c`**: Portable C scalar fallback (always compiled).
- **`<kernel>_simd_avx2.c`**: AVX2 + FMA vectorized implementation (compiled with `-mavx2 -mfma`).
- **`<kernel>_simd_avx512.c`**: AVX-512 implementation (compiled with `-mavx512f -mavx512dq -mfma`).
- **`<kernel>_simd_dispatch.c`**: Runtime CPU feature detection via `__builtin_cpu_supports` and
  function pointer binding.

---

## 3. Kernel Guidelines

- **Latency Hiding:** Modern x86 FMA units require multiple independent accumulators to prevent
  pipeline stalls. Unroll reduction loops with 4 to 8 accumulators.
- **Unaligned Memory Loads:** Use `_mm256_loadu_ps` and `_mm512_loadu_ps` unless allocations
  are explicitly aligned to 32 or 64 bytes with `posix_memalign()`.
- **Mandatory Scalar Fallback:** Every handwritten SIMD function must have a matching scalar
  implementation to guarantee portability on non-x86 hardware.
- **Benchmark Verification:** Validate vectorized speedups against `benchmark_perf.c` in `tests/`.
