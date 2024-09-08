library(stringr)
library(readr)
CopyNLMEControls <- function(home_dir){
  this_model <- 1
  for(this_model in 1:72){
    file <- file.path(home_dir,"NLME","run","0",str_pad(this_model, 2, pad = "0"),paste0("NLME_0_",str_pad(this_model, 2, pad = "0"),".mmdl"))
    control <- readtext::readtext(file, verbosity  = 0)$text
    control <- str_replace(control,"numIterations = 10","numIterations = 9999")
    new_dir <- file.path(home_dir,"NLME",this_model)
    if(!file.exists(new_dir)) dir.create(new_dir)
    new_file <- file.path(new_dir,paste0("Run",this_model,".mmdl"))
    writeLines(control[[1]], new_file)
    }
}
