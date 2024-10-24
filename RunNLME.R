rm(list=ls())
Linuxnlme_dir <- "/home/user/InstallDirNLME/"
Windowsnlme_dir <- "C:/Program Files/Certara/NLME_Engine"
home_dir <- getwd()

if (Sys.info()['sysname'] == "Linux") {
  gcc_dir <- dirname(system("which gcc", intern = TRUE))
  if (grepl("Ubuntu", Sys.info()["version"])) {
    Sys.setenv("PML_BIN_DIR" = "UBUNTU2204")
  }
  INSTALLDIRStandard <- Linuxnlme_dir

} else {
  # for Windows NLME engine
  gcc_dir <- "C:\\Program Files\\Certara\\mingw64"

  INSTALLDIRStandard <- Windowsnlme_dir
}

nlme_dir <-INSTALLDIRStandard
library(R.utils)
library(Metrics)
library(readr)
library(stringr)
library(dplyr)
library(data.table)
library(lubridate)
library(readtext)
library(Certara.RsNLME)

source(file.path(home_dir,"check_data.R"))
source(file.path(home_dir,"CleanUp.R"))
source(file.path(home_dir,"CompileResultsNLME.R"))

message("Veryfying data sets")
for(i in 1:72){
  check_data(home_dir,i)
}

run_nlme(home_dir, nlme_dir)
