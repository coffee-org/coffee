---
description: Perform compile and verification testing on coffee.
---

# Compile Test

To verify the build:
1. Compile standalone build:
   ```bash
   make -C _build -j$(nproc)
   ```
2. Or compile in-tree with milk:
   ```bash
   make -C /home/oguyon/src/milk-framework-dev/_build coffee-all -j$(nproc)
   ```
3. Verify there are no warnings or errors.
4. Check code size ratchet constraints:
   ```bash
   ./scripts/check_code_size.sh
   ```
5. Verify shared library exports:
   ```bash
   nm -D _build/libcoffeeOptSystProp.so | grep "T OptSystProp"
   nm -D _build/libcoffeeAOsystSim.so | grep "T AOsystSim"
   nm -D _build/libcoffeecoronagraphs.so | grep "T coronagraphs"
   nm -D _build/libcoffeePIAACMCsimul.so | grep "T PIAACMCsimul"
   ```
