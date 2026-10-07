---
description: Testing guidelines for verification of coffee modules.
---

# Testing Practices

Ensure all changes are validated for regressions and correctness before completing tasks.

## 1. Local Compile & Test
After making changes to source or CMake configuration:
1. Rebuild the project:
   ```bash
   make -C _build -j$(nproc)
   ```
2. Verify in-tree build with milk:
   ```bash
   make -C /home/oguyon/src/milk-framework-dev/_build coffee-all -j$(nproc)
   ```
3. Run code size ratchet check:
   ```bash
   ./scripts/check_code_size.sh
   ```
4. Verify shared libraries export valid symbols:
   ```bash
   nm -D _build/libcoffeeOptSystProp.so | grep "T OptSystProp"
   nm -D _build/libcoffeecoronagraphs.so | grep "T coronagraphs"
   ```

## 2. Regression Testing
When fixing a bug:
1. Create or identify a test case that triggers the bug.
2. Confirm the bug reproduces prior to code changes.
3. Confirm the bug is fully resolved after your changes, without introducing build warnings
   or breaking other tests.
