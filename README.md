# coffee — Coronagraph Optimization For Fast Exoplanets Exploration

[![License: GPL v3][gpl-badge]][gpl-url]

`coffee` is a high-performance optical propagation, coronagraph design, and extreme-AO simulation
plugin suite for the `milk` (`feat/bash-engine`) real-time framework.

[gpl-badge]: https://img.shields.io/badge/License-GPL%20v3-blue.svg
[gpl-url]: http://www.gnu.org/licenses/gpl-3.0

---

## Modular Architecture

| Module | Library | CLI Prefix | Purpose |
|---|---|---|---|
| `OptSystProp/` | `libcoffeeOptSystProp.so` | — (C API) | Fresnel beam propagation engine |
| `coronagraphs/` | `libcoffeecoronagraphs.so` | `coro.` | 2D prolate & APLC/PIAA scans |
| `PIAACMCsimul/` | `libcoffeePIAACMCsimul.so` | `coffeepiaacmcsim.` | PIAACMC/APLCMC design suite |
| `AOsystSim/` | `libcoffeeAOsystSim.so` | `aosystsim.` | Pyramid WFS, DM & AO simulation |

Upstream `milk` optical plugins required by `coffee`:
- `OpticsMaterials` (`libmilkOpticsMaterials.so`): Refractive index and dispersion lookup (`SiO2`,
  `CaF2`, `Mirror`, gases, resists).
- `WFpropagate` (`libmilkWFpropagate.so`): Angular-spectrum and Fresnel wavefront propagation.

---

## Building and Installing

### Prerequisites (Ubuntu / Debian)

```bash
sudo apt install build-essential cmake pkg-config libcfitsio-dev \
    libgsl-dev libfftw3-dev libreadline-dev libncurses-dev
```

### Out-of-Tree Build Against Installed `milk` (`feat/bash-engine`)

Ensure `milk` (`feat/bash-engine`), `OpticsMaterials`, and `WFpropagate` are built and installed
to `<milk-install-prefix>`, then build and install `coffee`:

```bash
cmake -S . -B _build \
    -DMILK_ROOT=<milk-install-prefix> \
    -DCMAKE_INSTALL_PREFIX=<milk-install-prefix>
make -C _build -j$(nproc) install
```

### In-Tree Build Inside `milk`

When symlinked or fetched under `<milk-src>/plugins/coffee-src/`:

```bash
make -C <milk-src>/_build -j$(nproc) install
```

---

## Interactive CLI Usage (`milk-cli` Bash Engine)

Add `<milk-install-prefix>/bin` to `PATH` and launch `coffee`:

```bash
export PATH="<milk-install-prefix>/bin:$PATH"
coffee
```

`coffee` preloads `coffeeOptSystProp`, `coffeecoronagraphs`, `coffeePIAACMCsimul`, and
`coffeeAOsystSim` into the `milk-cli` Bash engine:

```bash
# List loaded modules and commands inside milk-cli
@help
@help coffeepiaacmcsim
@help coro
@help aosystsim

# Inspect command arguments
coffeepiaacmcsim.run -h
coro.cormk2Dprolate -h

# Set numeric CLI variables via @calc and load/save FITS images
@calc "PIAACMC_size = 1024"
@calc "PIAACMC_lambda = 550.0"
iofits.loadfits "pup_1024.fits" pupin 2
exit
```

---

## PIAACMC & APLCMC Design Quickstart

Use the unified `coffee-piaacmc` CLI in a clean working directory:

```bash
mkdir my_coronagraph && cd my_coronagraph

# 1. Initialize workspace from built-in example (e.g., WFIRST/Roman 1024 mirror PIAACMC)
coffee-piaacmc init --example 1

# 2. Inspect or customize parameters in piaacmc.conf
coffee-piaacmc config show
coffee-piaacmc config set mlambdaB=10 nblambda=8 NBthreads=4

# 3. Run monochromatic design (steps 000..017)
coffee-piaacmc run --mono

# 4. Run polychromatic multi-zone FPM calibration + sag optimization
coffee-piaacmc run --poly --jobs 4 --search-time 600

# 5. Evaluate polychromatic contrast and off-axis transmission
coffee-piaacmc eval --level 0
coffee-piaacmc status
```

See [PIAACMCsimul/README.md](PIAACMCsimul/README.md) and
[coronagraphs/README.md](coronagraphs/README.md) for complete reference documentation.
