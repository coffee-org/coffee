## Summary of Changes
<!-- Provide a clear, concise overview of what changed and why. -->

## Changes Checklist
- [ ] Code compiles cleanly with zero warnings (`-Wall -Wextra`)
- [ ] Code size ratchet passes (`./scripts/check_code_size.sh`)
- [ ] Adheres to C code style guide (Allman braces, <= 100 char lines, column-aligned params)
- [ ] No allocations (`malloc`) in hot simulation loops (Fresnel propagation, PSF compute)
- [ ] No resource/memory leaks (verified with ASan/UBSan or Valgrind)
- [ ] New source files reflected in `CMakeLists.txt` and corresponding module documentation

## Testing Performed
<!-- Describe how these changes were tested (commands executed, scripts run, output verified). -->
```bash
# Example:
make -C _build -j$(nproc)
./scripts/check_code_size.sh
```

---
<!--
If agentic AI tools were used, disclose below per project guidelines:
Implemented by <model name>. Reviewed and signed off by O. Guyon.
Followed by a concise technical summary of what task the model performed.
-->
Implemented by . Reviewed and signed off by O. Guyon.
