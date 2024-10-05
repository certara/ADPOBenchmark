rm(list=ls())
## path to python.exe can be found with where python.exe command from DOS command line
#WindowspyDarwinInterpreter <- "C:/Users/msale/AppData/Local/Programs/Python/Python310/python.exe"
Windowsnmfe_path <- "C:/nm75g64/util/nmfe75.bat"
Linuxnmfe_path <- "/opt/nm751/util/nmfe75"
home_dir <- getwd()

if (Sys.info()['sysname'] == "Linux") {
  nmfe_path <- Linuxnmfe_path
} else {
  nmfe_path <- Windowsnmfe_path
}

library(Metrics)
library(ggplot2)
library(readr)
library(stringr)
library(dplyr)
library(data.table)
library(lubridate)
library(readtext)
library(xml2)

source(file.path(home_dir,"check_data.R"))
source(file.path(home_dir,"ReadNM_xml.R"))
source(file.path(home_dir,"CompileResultsNONMEM.R"))
source(file.path(home_dir,"CleanUp.R"))
source(file.path(home_dir,"GetTrueParms.R"))
for(i in 1:72){
  check_data(home_dir,i)
}

run_NONMEM(home_dir, nmfe_path)
