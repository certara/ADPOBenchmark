rm(list=ls())
Windowsnmfe_path <- "C:/nm74g64/util/nmfe74.bat"
home_dir <- getwd()
if (Sys.info()['sysname'] == "Linux") {
  nmfe_path <- "/opt/nm751/util/nmfe75"
} else {
  nmfe_path <- Windowsnmfe_path
}

library(Metrics)
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

message("Veryfying data sets")
for(i in 1:72){
  check_data(home_dir,i)
}

run_NONMEM(home_dir, nmfe_path)
message("Done with NONMEM, consider rebooting before running RunNLME.R")
