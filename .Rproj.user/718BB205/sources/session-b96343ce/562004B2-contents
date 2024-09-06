rm(list=ls())
home_dir <- getwd()
source(file.path(home_dir,"ReadNM_xml.R"))
source(file.path(home_dir,"CompileResultsNONMEM.R"))
# still need to add ipred table to nlme!!
source(file.path(home_dir,"CompileResultsNLME.r"))
source(file.path(home_dir,"CalcRSME.r"))
run_NONMEM(home_dir)
run_nlme(home_dir,"standard")
