---
description: Unified help message format for coffee standalone binaries.
---

# Help Message Standard

Every standalone FPS executable in `coffee` MUST support a help screen via `-h` or `--help`.

## Help Message Layout
The help output must follow this section order:
1. **NAME** — program name and brief description
2. **USAGE** — command-line syntax synopsis
3. **DESCRIPTION** — brief description of functionality and parameters
4. **OPTIONS** — detailed list of supported flags, defaults, and arguments
5. **EXAMPLES** — 1-3 concrete usage examples

## Options & Arguments Formatting
Align positional arguments, keywords, and descriptions cleanly:
```
  [0] .beam_diam  FLOAT64    Beam diameter [m] (default: 8.0)
  [1] .wavelength FLOAT64    Wavelength [m] (default: 1.6e-6)
  [2] .zprop      FLOAT64    Propagation distance [m] (default: 100.0)
```

## Color Formatting and Suppression
Apply standard ANSI colors where terminal output is supported:
* **Section Headings** (`NAME`, `USAGE`, etc.): **Bold** or **Bold Cyan**
* **Commands / Executables**: **Bold Green**
* **Option Flags / Keywords**: **Regular Green**
* **Variables / Placeholders**: **Regular Magenta**
* **Errors**: **Bold Red**

### Suppression Standard
Support the [no-color.org](https://no-color.org/) standard: if the `NO_COLOR` environment
variable is set, suppress all ANSI escape sequence outputs.

## Error Handling on Arguments
If required command-line arguments are missing or invalid:
1. Print a clear error message to `stderr`.
2. Print the usage synopsis.
3. Exit with a non-zero code.
