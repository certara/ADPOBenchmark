library(readr)
library(stringr)
setwd("c:/git/adpoBenchmark/nonmem")
i = 1
for(i in 1:72){
  control <- read_file(file.path(i,paste0("Run",i,".mod")))
  control <- str_replace(control,"c:/git/adpo_speed/data/","c:/git/adpoBenchmark/data/")
  write_file(control,file.path(i,paste0("Run",i,".mod")))
}

setwd("c:/git/adpoBenchmark/nlme")
i = 1
for(i in 1:72){
  control <- read_file(file.path(i,paste0("run",i,".mmdl")))
  control <- str_replace(control,"c:/git/adpo_speed/data/","c:/git/adpoBenchmark/data/")
  write_file(control,file.path(i,paste0("run",i,".mmdl")))
}
