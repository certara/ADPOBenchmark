# 1.0 Setup environment ----
# Note, the paths assigned below is the default NLME-Engine installation directory for Linux and Windows
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

# 2.0 Load required packages ----
library(R.utils)
library(Metrics)
library(readr)
library(stringr)
library(dplyr)
library(data.table)
library(lubridate)
library(readtext)
library(Certara.RsNLME)

# 3.0 Source functions -----
source(file.path(home_dir,"check_data.R"))
source(file.path(home_dir,"CleanUp.R"))
source(file.path(home_dir,"CompileResultsNLME.R"))

# 4.0 Verify data ----
message("Veryfying data sets")
for(i in 1:72){
  check_data(home_dir,i)
}

# 5.0 Authenticate NLME License ----
obtain_NLMELicense()

# 6.0 Run NLME benchmark analysis ----
run_nlme(home_dir, nlme_dir)
