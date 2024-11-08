# 1.0 Setup environment ----
# Note, update the NONMEM installation paths given your installation location and applicable version
rm(list=ls())
home_dir <- getwd()
if (Sys.info()['sysname'] == "Linux") {
  nmfe_path <- "/opt/nm751/util/nmfe75"
} else {
  nmfe_path <- "C:/nm74g64/util/nmfe74.bat"
}

# 2.0 Load required packages ----
library(Metrics)
library(readr)
library(stringr)
library(dplyr)
library(data.table)
library(lubridate)
library(readtext)
library(xml2)

# 3.0 Source functions -----
source(file.path(home_dir,"check_data.R"))
source(file.path(home_dir,"ReadNM_xml.R"))
source(file.path(home_dir,"CompileResultsNONMEM.R"))
source(file.path(home_dir,"CleanUp.R"))

# 4.0 Verify data ----
message("Veryfying data sets")
for(i in 1:72){
  check_data(home_dir,i)
}

# 6.0 Run NONMEM benchmark analysis ----
run_NONMEM(home_dir, nmfe_path)
message("Done with NONMEM, consider rebooting before running additional benchmark analysis.")
