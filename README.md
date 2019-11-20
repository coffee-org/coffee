 [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/4bd280f0915f4174823fa89dc4758100)](https://www.codacy.com/app/oguyon/coffee?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=coffee-org/coffee&amp;utm_campaign=Badge_Grade)
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;dev branch: [![Build Status dev](https://travis-ci.org/coffee-org/coffee.svg?branch=dev)](https://travis-ci.org/coffee-org/coffee)
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;master branch: [![Build Status](https://travis-ci.org/coffee-org/coffee.svg?branch=master)](https://travis-ci.org/coffee-org/coffee)

---

IMPORTANT NOTE: coffee uses git submodules. Use `git clone --recursive` (see Downloading and Installing section)

---


# Coronagraph Optimization For Fast Exoplanet Exploration (coffee)

## Overview

Coronagraph optimization routines accessible through a command line interface (CLI). Holds images in RAM, with image stream support (shared memory with low-latency IPC support).

Written in C, optimized for performance.

Executable launches a command line interface (CLI). Type "help" in the CLI to get started.

coffee uses [milk](https://github.com/milk-org/milk)


---



### Install steps


```bash
git clone --recursive https://github.com/coffee-org/coffee coffee
cd coffee
mkdir _build
cd _build
cmake ..
# If you use NVIDIA GPUs, install cuda and magma libraries, and use "cmake .. -DUSE_MAGMA=ON"
make
sudo make install
```



### Post-installation 

You may need to add /usr/local/lib to LD_LIBRARY_PATH environment variable:
```bash
echo "/usr/local/lib" > usrlocal.conf
sudo mv usrlocal.conf /etc/ld.so.conf.d/
sudo ldconfig -v
```


Unless you have a separate install of milk on your system, a symbolic link to milk is required for milk scripts that are included in cacao:

```bash
sudo ln -s /usr/local/bin/coffee /usr/local/bin/milk
```

Create local bin directory and sym links to all cacao- scripts:
```bash
mkdir bin
cd bin
find /home/scexao/src/coffee/ -executable -type f -name "milk-*" -print0 | xargs -0 -I {} ln -s {} .
find /home/scexao/src/coffee/ -executable -type f -name "cacao-*" -print0 | xargs -0 -I {} ln -s {} .
```

Add milk and cacao executable scripts to PATH environment variable. Add this line to the .bashrc file (change source code location as needed):
```bash
export PATH=$PATH:/home/myname/src/coffee/bin
export PATH=$PATH:/home/myname/src/coffee/src/CommandLineInterface/scripts
```

OPTIONAL: Create tmpfs disk for high performance I/O:
```bash
echo "tmpfs /milk/shm tmpfs rw,nosuid,nodev" | sudo tee -a /etc/fstab
sudo mkdir -p /milk/shm
sudo mount /milk/shm
```



---



## Reporting bugs, issues

Report bugs and issues on [this page]( https://github.com/coffee-org/coffee/issues )


## Contributing to project


See [coding standards]( https://coffee-org.github.io/coffee/page_coding_standards.html ) 


---


## Documentation

[Online documentation]( http://coffee-org.github.io/coffee/index.html ) 


## Libraries

The following libraries are used:

- readline, for reading the command line input
- ncurses-dev
- flex, for parsing the command line input
- bison, to interpret the command line input
- fftw, for performing Fourier Transforms
- gsl, for math functions and tools
- fitsio, for reading and writing FITS image files
- CUDA, CuBLAS, MAGMA for GPU acceleration (optional)

If you use NVIDIA GPUs, install cuda and magma libraries, and add "--enable-cuda and --enable-magma" options to the configure command.



## Getting Started

All functions are accessible from the command line interface (CLI). Enter the CLI by typing "coffee" and type "help" for instructions.


## LICENCE


[GNU General Public License v3.0]( https://github.com/coffee-orga/coffee/blob/master/LICENCE.txt )

