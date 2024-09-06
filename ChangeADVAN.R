library(readr)
library(stringr)
setwd("d:/git/dual_numbers/compileresults/nonmem")
i = 1
for(i in 1:72){
  control <- read_file(file.path(i,paste0("Run",i,".mod")))
  control <- str_replace(control,"ADVAN6","ADVAN13")
  write_file(control,file.path(i,paste0("Run",i,".mod")))
}
 
