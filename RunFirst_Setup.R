rm(list=ls())
## path to python.exe can be found with where python.exe command from DOS command line
WindowspyDarwinInterpreter <- "C:/Users/msale/AppData/Local/Programs/Python/Python310/python.exe"
Windowsnmfe_path <- "C:/nm75g64/util/nmfe75.bat"
WindowsINSTALLDIRADPO <- "D:/NLME_Engine_ADPO"
WindowsINSTALLDIRnoHessian <- "D:/NLME_Engine_noHessian"
Windowsnlme_dir <- "C:\\Program Files\\Certara\\NLME_Engine"
Windowsgcc_dir <- "C:\\Program Files\\Certara\\mingw64"
Linuxnlme_dir <- "/home/user/InstallDirNLME/"
LinuxpyDarwinInterpreter <- "/home/user/venv/bin/python3"
Linuxnmfe_path <- "/opt/nm751/util/nmfe75"


home_dir <- getwd()

# used by NLME
if (Sys.info()['sysname'] == "Linux") {
  pyDarwinInterpreter <- LinuxpyDarwinInterpreter
  nmfe_path <- Linuxnmfe_path

  gcc_dir <- dirname(system("which gcc", intern = TRUE))
  #nlme_dir <- "/home/user/InstallDirNLME/"
  if (grepl("Ubuntu", Sys.info()["version"])) {
    Sys.setenv("PML_BIN_DIR" = "UBUNTU2204")
  }

  INSTALLDIRStandard <- Linuxnlme_dir
 # INSTALLDIRADPO <- "home/user/NLME_Engine_ADPO"
#  INSTALLDIRnoHessian <- "home/user/NLME_Engine_noHessian"

} else {
  pyDarwinInterpreter <- WindowspyDarwinInterpreter
  nmfe_path <- Windowsnmfe_path

  # for Windows NLME engine
  gcc_dir <- Windowsgcc_dir
  nlme_dir <- Windowsnlme_dir

  INSTALLDIRStandard <- nlme_dir
 # INSTALLDIRADPO <- WindowsINSTALLDIRADPO
#  INSTALLDIRnoHessian <- WindowsINSTALLDIRnoHessian
}

nlme_dirs <- c(standard = INSTALLDIRStandard) #,
 #              ADPO = INSTALLDIRADPO,
#               NoHessian = INSTALLDIRnoHessian)

library(Metrics)
#library(ggplot2)
library(readr)
library(stringr)
library(dplyr)
library(data.table)
library(lubridate)
library(readtext)
library(Certara.RsNLME)
library(Certara.RDarwin)
#library(xml2)

source(file.path(home_dir,"make_data.R"))
#source(file.path(home_dir,"check_data.R"))
#source(file.path(home_dir,"ReadNM_xml.R"))
#source(file.path(home_dir,"CompileResultsNONMEM.R"))
source(file.path(home_dir,"CleanUp.R"))
##source(file.path(home_dir,"CompileResultsNLME.R"))
#source(file.path(home_dir,"GetTrueParms.R"))

make_data(home_dir, pyDarwinInterpreter, nmfe_path, nlme_dir, gcc_dir)
