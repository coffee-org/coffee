# Getting Started {#page_PIAACMCsimul_GettingStarted}

*Source: PIAACMCsimul_GettingStarted.md in PIAACMCsimul/doc*


## 1. Copy scripts to local directory

Create symbolic link to CLI in system or user path:
	
	sudo ln -s /<srcdir>/bin/coffee /usr/local/bin/coffee

Create symbolic link to syncscripts into work directory:

	cd <workdir>
	ln -s /<srcdir>/src/PIAACMCsimul/scripts/syncscripts .

Import and setup scripts in work directory (run from <workdir>):

	./syncscripts

Examples can be loaded in `./examples/` with:

	./syncscripts -x


***



## 2. Configuration Setup


### 2.1. Setup from main script pre-defined configurations

Print help for top level script:

	./runPIAACMCdesign -h

The quickest way to get started is to setup and modify one of the example design scripts provided:

	./runPIAACMCdesign -e 0

This will configure all necessary files for a specific configuration. You can view parameters with :

	./runPIAACMCdesign -l
	
If optimizing in APLC mode (no PIAA optics, type (after the -e command above):

	./runPIAACMCdesign -a
	
	
### 2.2. Setup from example script

If examples are loaded, one of the `example/setup_XXXX` script can be copied in the working directory and executed to setup the corresponding configuration
	

***



## 3. Full Design

Design is done with the "-m" option. This command will execute all design steps from 0 to #step-1:

	./runPIAACMCdesign -m <#step>

To run the full polychromatic design process:

	./runPIAACMCdesign -m 200

Note that this may take a long time to run ... (days ?). The polychromatic design will:

- Construct response to focal plane mask zones if it does not already exist. This step runs in tmux sessions, and will write a FPMresp...fits file
- Execute a search for best solution.


***


## 4. Focal plane optimization

In the focal plane optimization mode, the search status appears in files:

	tail -f piaacmcconf_i000/linoptval.txt
	tail -f piaacmcconf_i000/mode13*.opt.txt
	
	

In the search mode, the code runs for a pre-determined amount of time. You can track progress in file `timeused.txt`. The first number is the search time elapsed, and the second number is the amount of total search time. To change / shorten the search time : 

	cat "3600" > searchtime.txt  # sets total search time to 1hr

You can also stop the search at anytime by:

	touch piaacmcconf_i000/stoploop13.txt


***


## 5. Inspecting results

Modes 300 and 400-402 are used to evaluate solution performance:

	./runPIAACMCdesign -m 300
             #         300  : compute polychromatic on-axis PSF -> psfi0.fits
             #         400  : evaluation, level 0 (on-axis PSF)
             #         401  : evaluation, level 1 (level0 + extended source with OPD errors)
             #         402  : evaluation, level 2 (level 1 + transmission curve)

