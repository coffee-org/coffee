---
description: Inspect compiler assembly outputs.
---

# Inspect Machine Code

To inspect the generated assembly for hot loop vectorization in `coffee`:
1. Compile the target source file to assembly:
   ```bash
   gcc -S -O3 -mavx2 -mfma -I. -IOptSystProp \
       OptSystProp/OptSystProp_FresnelProp.c -o OptSystProp_FresnelProp.s
   ```
2. Inspect the output `.s` file for vector instructions (e.g. `ymm` or `zmm` registers,
   `vfmadd213ps`, `vmulps`, `vaddps`).
3. Confirm that auto-vectorization did not fall back to scalar instructions in hot regions.
