library(stringr)
library(readr)
CopyNONMEMControls <- function(home_dir){
  this_model <- 1
  for(this_model in 1:72){
    file <- file.path(home_dir,"NONMEM","run","0",str_pad(this_model, 2, pad = "0"),paste0("NM_0_",str_pad(this_model, 2, pad = "0"),".mod"))
    control <- readtext::readtext(file, verbosity  = 0)$text
    control <- str_replace(control,";;\\$EST","$EST")
    control <- str_replace(control,";;\\$COV","$COV")
    new_dir <- file.path(home_dir,"NONMEM",this_model)
    if(!file.exists(new_dir)) dir.create(new_dir)
    new_file <- file.path(new_dir,paste0("Run",this_model,".mod"))
    writeLines(control[[1]], new_file)
    }
}
