# Always compile after editing source code

After modifying any C/CMake source file in the `coffee` tree, you **must** run the
compile-test workflow to verify the build still succeeds before considering the task complete.

If the build fails, fix the errors and rebuild until it passes.

Verification steps:
```bash
# Standalone build
make -C _build -j$(nproc)

# Or in-tree build (from milk build directory)
make -C /home/oguyon/src/milk-framework-dev/_build coffee-all -j$(nproc)
```

Also verify the code size ratchet:
```bash
./scripts/check_code_size.sh
```
