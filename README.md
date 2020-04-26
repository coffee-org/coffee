 [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/4bd280f0915f4174823fa89dc4758100)](https://www.codacy.com/app/oguyon/coffee?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=coffee-org/coffee&amp;utm_campaign=Badge_Grade)
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;dev branch: [![Build Status dev](https://travis-ci.org/coffee-org/coffee.svg?branch=dev)](https://travis-ci.org/coffee-org/coffee)
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;master branch: [![Build Status](https://travis-ci.org/coffee-org/coffee.svg?branch=master)](https://travis-ci.org/coffee-org/coffee)

---

IMPORTANT NOTE: coffee uses git submodules. Use `git clone --recursive` (see Downloading and Installing section)

---


# Coronagraph Optimization For Fast Exoplanet Exploration (coffee)

---


## Installing cacao


coffee runs on Linux/x86 systems.

&#x26A0;
**coffee requires milk**: Install milk-package (ideally with GPU support) prior to installing coffee. See instructions on the [milk page](https://github.com/milk-org/milk-package).





### Download

	git clone --recursive https://github.com/coffee-org/coffee coffee


### Compile

	cd coffee
	mkdir _build; cd _build
	cmake ..
	make
	

### Install


	sudo make install

Will install milk in /usr/local/milk-&lt;version&gt;. Multiple versions of coffee can coexist in separate coffee-&lt;version&gt; directories. To select the version to be used:

	sudo ln -s /usr/local/coffee-<version> /usr/local/coffee

	
Add environment variables. Add to .bashrc file or similar :

	export COFFEE_ROOT=${HOME}/src/coffee  # point to source code directory. Edit as needed.
	export COFFEE_INSTALLDIR=/usr/local/coffee
	export PATH=${PATH}:${COFFEE_INSTALLDIR}/bin






---



## Reporting bugs, issues

Report bugs and issues on [this page]( https://github.com/coffee-org/coffee/issues )


## Contributing to project


See [coding standards]( https://coffee-org.github.io/coffee/page_coding_standards.html ) 


---


## Documentation

[Online documentation]( http://coffee-org.github.io/coffee/index.html ) 



## Getting Started

All functions are accessible from the command line interface (CLI). Enter the CLI by typing "coffee" and type "help" for instructions.


## LICENCE


[GNU General Public License v3.0]( https://github.com/coffee-orga/coffee/blob/master/LICENCE.txt )

