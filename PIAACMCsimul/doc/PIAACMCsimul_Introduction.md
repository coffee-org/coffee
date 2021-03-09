
# Introduction {#page_PIAACMCsimul_Introduction}


## Scope

Diffraction-based PIAACMC simulation / optimization

- Uses Fresnel propagation engine ( OptSystProp.c OptSystProp.h ) between optical elements
- Computes linear perturbations around current design for optimization
- Automatically computes multiple Lyot stop to follow variable conjugation in different parts of the beam
- Polychromatic propagations
- Fits aspheric PIAA shapes on basis of radial cosines and 2-D Fourier modes


## Usage

Scripts to run the software are located within the source code directory:

	./src/PIAACMCsimul/scripts/

The scripts can be linked to your working directory by executing the following command:

	ln -s $PWD/syncscripts /myworkdirectory/syncscripts

Then, execute in your work directory:

	./syncscripts
	
This will install all required scripts in workdirectory and install any packages required.

Code is composed of a several layers (from high to **low**) :


Script            |     Description
------------------|-----------------------------------------------------------
runPIAACMCdesign  | Top level script ("-h" for help)
run               |  Main script, calls **runopt** script
runopt            |  Optimize a PIAACMC design or run an existing design, calls **sim** script
sim               |  calls **runPIAACMC** script
runPIAACMC        |  lower-level script calling C-written executable



@dot
digraph scripts_flow {
	size ="4,4";
	ratio = auto;
  
  "runPIAACMCdesign" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white">runPIAACMCdesign</font></td></tr>
  <tr><td align="left">Top level script</td></tr>
  <tr><td align="left">Executes lower scripts</td></tr>
  </table>> ];

  "run" [ style = "filled, bold" penwidth = 5 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white">run</font></td></tr>
  <tr><td align="left">Main script</td></tr>
  <tr><td align="left">Reads input variables</td></tr>
  </table>> ];
  
  "runopt" [ style = "filled" penwidth = 1 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white">runopt</font></td></tr>
  <tr><td align="left">Optimize a PIAACMC design or</td></tr>
  <tr><td align="left"> run an existing design</td></tr>
  </table>> ];

  "sim" [ style = "filled" penwidth = 1 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white">sim</font></td></tr>
  <tr><td align="left"> </td></tr>
  </table>> ];

  "runPIAACMC" [ style = "filled" penwidth = 1 fillcolor = "white" fontname = "Courier New" shape = "Mrecord" 
  label =<<table border="0" cellborder="0" cellpadding="3" bgcolor="white"><tr><td bgcolor="black" align="center" colspan="2"><font color="white">runPIAACMC</font></td></tr>
  <tr><td align="left">lower-level script calling C-written executable</td></tr>
  </table>> ];


  
  runPIAACMCdesign-> run [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
  run-> runopt [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
  runopt -> sim [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
  sim -> runPIAACMC [ penwidth = 5 fontsize = 28 fontcolor = "black" label = "" ];
}
@enddot





