# ADPOBenchmark 
## Overview:
ADPO is an abreviation for Automatic differentiation parameter optimization. The method has been developed by Certara
and may present an opportunity for improved performance of parameter optimization in first order conditional 
method for mixed effects models. Briefly, rather than the finite difference method being used to calculate the
gradients WRT ETAs in the conditional step (the "inner loop") a close form solution to the derivative is calculated. 
This approach has been successful in neural network [(https://en.wikipedia.org/wiki/Automatic_differentiation#:~:text=Automatic%20differentiation%20is%20particularly%20important,without%20a%20manually%2Dcomputed%20derivative.)].
Certara has adapted this method for mixed effects models.
## ADPOBenchmark
ADPOBenchmark runs performance benchmarks in NONMEM (Icon PLC - https://www.iconplc.com/solutions/technologies/nonmem) and NLME (Certara - https://www.certara.com/app/uploads/2020/06/BR_PhoenixNLME-v4.pdf) for 72 ODE models (ADVAN6 and DVERK ODE solvers). These models
are fit to simulated data. The simulate data is generate from the same 72 models, so the estimation model is fitting, to the extent possible, the "true" model in all 72 cases. Further, the initial estimates for each model are those used for the simulation.
Running the benchmark requires the following software:
NONMEM 
NLME engine.  A 30 day complimentary license is available from Certara (Support@certara.com)

Python and pyDarwin (https://certara.github.io/pyDarwin/html/index.html) are requred if the control files and simulated data are to be generated. In principle the NONMEM control file, NLME metamodel files and simulated data sets (sim*.csv) in this 
repository could be used, but the path to the data sets, in the controls is likely to be incorrect.

It is likely that some customization of the script will be needed. At minimum, the path the nmfe??.bat (batch file for executing NONMEM) and the paths for NLME for Windows and/or Linux must be provided. These are specified in the first few lines of the
"RunFirst_setup.r" file and default values are given below

WindowspyDarwinInterpreter <- "C:/Users/msale/AppData/Local/Programs/Python/Python310/python.exe"
Windowsgcc_dir <- "C:\\Program Files\\Certara\\mingw64"
Windowsnlme_dir <- "C:\\Program Files\\Certara\\NLME_Engine"
Windowsnmfe_path <- "C:/nm75g64/util/nmfe75.bat" 
Linuxnlme_dir <- "/home/user/InstallDirNLME/"
LinuxpyDarwinInterpreter <- "/home/user/venv/bin/python3"
Linuxnmfe_path <- "/opt/nm751/util/nmfe75"

