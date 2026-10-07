---
trigger: always_on
---

# Code Size, Function Scope & Architectural Intent

To keep the codebase maintainable, readable, and cleanly structured across architectural layers,
all C source files and functions must adhere to strict size and intent constraints.

## 1. Size Constraints

Excluding vendored third-party code and build artifacts, limits are enforced as follows:

| Metric | Soft Limit (refactor when editing) | Hard Limit (CI fails) |
|---|---|---|
| File length | 600 lines | 1000 lines |
| Function body | 60 lines | 150 lines |
| `main()` body | 40 lines | 80 lines |
| Nesting depth | 3 levels | 5 levels |
| Parameters | 6 arguments | 8 arguments (use a struct) |
| Cyclomatic complexity | 10 branches | 20 branches |

A ratchet baseline (`scripts/code_size_baseline.txt`) tracks legacy violations. Any new violation
or growth in an existing baseline entry fails the check via `scripts/check_code_size.sh`.

## 2. Function Intent & Readability

- **Single Purpose**: Every non-trivial function must begin with a Kernel-Doc comment whose
  first line defines a single clear responsibility. If the summary requires "and" to chain
  multiple tasks, split the function into cohesive helpers.
- **Naming Pattern**: Use `<module>_<verb>_<object>` (e.g., `OptSystProp_propagate_cube`). Avoid
  uninformative verbs (`do_`, `process_`, `handle_`) without a precise object.
- **Orchestrator Functions**: High-level workflow functions must read like a table of contents:
  a linear sequence of clear step calls, keeping mathematical loops in leaf compute modules.
- **Role of Drivers**: CLI and standalone runners must remain minimal: parse CLI options into
  a configuration struct, invoke the primary runner function, and map return codes to status.
- **Directory Documentation**: Every source directory should document its purpose, public
  headers, and permitted dependencies.

## 3. Exemptions

- Performance-critical SIMD kernels unrolled for vector throughput may exceed the function length
  limit if preceded on the line above by `/* size-exempt: <reason> */`.
- Prefer separating ISA specializations into dedicated files (`_scalar.c`, `_avx2.c`, `_avx512.c`)
  before applying size exemptions.

## 4. Remediation

When editing or creating code that nears these thresholds, consult the `refactor-c-source` skill
for step-by-step function decomposition and file reorganization patterns.
