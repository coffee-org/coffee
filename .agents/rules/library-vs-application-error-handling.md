---
description: Design principles for error boundaries between core simulation library and
  application files.
---

# Library vs Application Error Handling

Maintain a strict boundary between library compute functions (e.g., in `OptSystProp/`,
`coronagraphs/`, `PIAACMCsimul/`, `AOsystSim/`) and CLI / FPS application drivers.

## 1. Core Simulation Library Functions
- Must **never** call `exit()` or terminate the process.
- Must return error codes (`-1`, non-zero `errno_t`, or `NULL` pointers) to the caller.
- Clean up any temporary buffers before returning error codes.
- Should log warnings or errors to `stderr` only on severe, unrecoverable states.

## 2. CLI, FPS Drivers, and Application Layers
- Handle error codes returned by library functions.
- Present clean user-facing error messages to `stderr`.
- Decide whether to abort execution, log details, or exit cleanly using `exit(EXIT_FAILURE)`.
