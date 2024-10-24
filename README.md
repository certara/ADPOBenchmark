# ADPOBenchmark 
## Overview:
ADPO is an abbreviation for Automatic differentiation parameter optimization. The method has been developed by Certara
and may present an opportunity for improved performance of parameter optimization in first order conditional 
method for mixed effects models. Briefly, rather than the finite difference method being used to calculate the
gradients WRT ETAs in the conditional step (the "inner loop") a closed form solution to the derivative is calculated using the chain rule. This approach has been adapted from a similar approach used in neural networks.  [(https://en.wikipedia.org/wiki/Automatic_differentiation].\

## Software requirements:

Running the suite of benchmarks requires:

NONMEM (https://www.iconplc.com/solutions/technologies/nonmem)

NONMEM is run using the nmfe??.bat command, where ?? is the version of NONMEM (e.g., 7.5 -> nmfe75.bat)

NLME engine (https://www.certara.com/software/phoenix-nlme/)

  A complimentary 30 day NLME engine license is avaiable from support@Certara.com.
  
## ADPOBenchmark
ADPOBenchmark runs performance benchmarks in NONMEM (Icon PLC - https://www.iconplc.com/solutions/technologies/nonmem) and NLME (Certara - https://www.certara.com/app/uploads/2020/06/BR_PhoenixNLME-v4.pdf) for 72 ODE models (ADVAN6 and DVERK ODE solvers). These models
are fit to simulated data. The simulate data is generate from the same 72 models, so the estimation model is fitting, to the extent possible, the "true" model in all 72 cases. Further, the initial estimates for each model are those used for the simulation.


It is likely that some customization of the script will be needed. At minimum, the path the nmfe??.bat (batch file for executing NONMEM) and the paths for NLME for Windows and/or Linux must be provided. These are specified in the first few lines of the
"RunNONMEM.R" and "RunNLME.R" file and default values are given below

For RunNONMEM.R:

  Windowsnmfe_path <- "C:/nm74g64/util/nmfe74.bat" 

For RunNLME.R: 

  Windowsnlme_dir <- "C:/Program Files/Certara/NLME_Engine"
  
  Linuxnlme_dir <- "/home/user/InstallDirNLME/" 
  
  gcc_dir <- "C:\\Program Files\\Certara\\mingw64"


Three R script can be run. First run RunNONMEM.R (after confirming that the nmfe??.bat path is correct). After this is done, the user may want to reboot the computer, to insure that there is no sequence effect. Certara testing suggests that for the present # of models there is no signficant unreleased memory, but the user should feel free to reboot.
Next, run RunNLME.R, again after confirming that the paths are correct.
When these are done, the MakePlots.R can be run to generate the plot and table.

Please feel free to contact support@Certara.com for a complimentary 30 day NLME license to the NLME engine, or any questions.
