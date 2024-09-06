library(readr)
library(stringr)
setwd("c:/git/adpoBenchmark/nonmem")
i = 1
for(i in 1:72){
  control <- read_file(file.path(i,paste0("Run",i,".mod")))
  #control <- paste0(control,"\n$TABLE ID TIME IPRED DV NOAPPEND ONEHEADER NOPRINT FILE= Preds.dat")
  #control <- str_replace(control,"ERROR ","ERROR \n IPRED = F")
  #control <- str_replace(control,"Preds.dat","Run1Preds.dat")
  control <- str_replace_all(control,"(100,1200,3000)","100,4000")
  control <- str_replace_all(control,"(1,15,40)","1,50")
  write_file(control,file.path(i,paste0("Run",i,".mod")))

  }

setwd("c:/git/adpoBenchmark/nlme")
i = 1
for(i in 1:72){
  control <- read_file(file.path(i,paste0("run",i,".mmdl")))
  control <- str_replace(control,"c:/git/adpo_speed/data/","c:/git/adpoBenchmark/data/")
  write_file(control,file.path(i,paste0("run",i,".mmdl")))
}
