#run below to simulate data and control files
#
library(dplyr)
library(readr)
rm(list=ls())
LLOQ <- 0.1
home_dir <- "c:/git/adpobenchmark/"
set.seed(1)
nsubs <-  60
dose <-  100000 # mcg
WT_mean <- 70
WT_CV <- 0.2
mean_times <- c(0.5,2,6,12,24,48,96,120)
cv_times <- 0.1
data <- data.frame(matrix(-999, nrow=nsubs*8, ncol=7))
colnames(data) <- c("ID","TIME","AMT","DV","EVID","BQL","WT")
WTs <- exp(rnorm(nsubs,log(WT_mean),WT_CV))
for(this_id in 1:nsubs){
  times <- exp(rnorm(7,log(mean_times),cv_times))
  this_data <- matrix(-999,nrow=8,ncol=7)
  this_data[,1] <-  this_id
  this_data[,2] <- c(0,round(times,2))
  this_data[,3] <- c(dose,rep(".",7))
  this_data[,4] <- "."
  this_data[,5] <- c(1,rep(2,7))
  this_data[,6] <- 0 ## BQL
  this_data[,7] <- round(WTs[this_id],1) ## BQL
  data[(((this_id-1)*8+1):(this_id*8)),] <- this_data
}
write.csv(data,file=file.path(home_dir,"data","pydarwin_makedata","data_sim.csv"),quote = FALSE,row.names = FALSE)
## run pydarwin here, need to run before the NONMEM script
if(file.exists(file.path(home_dir,"data","pydarwin_makedata","run"))){
  unlink(file.path(home_dir,"data","pydarwin_makedata","run"),recursive = TRUE)
}
if(file.exists(file.path(home_dir,"data","OUT_0_0_0_0.dat"))){
  files <- list.files(file.path(home_dir,"data"),
                      full.names = TRUE,
                      recursive = FALSE)
  files <- files[grep("OUT_[0-9].*",files)]
  for(file in files){
    file.remove(file)
  }
}

setwd(file.path(home_dir,"data","pydarwin_makedata"))
shell("python -m darwin.run_search_in_folder .")
# ## set BQL values
i <- j <- l <- m <-  0

for(i in 0:1){ #
  for(j in 0:5){
    for(l in 0:1){
      for(m in 0:1){
        file <- file.path(home_dir,"data",paste0("OUT_",i,"_",j,"_",l,"_",m,".DAT"))
        sim_data <- read_table(file,skip=0,
                   col_names = FALSE)
        colnames(sim_data) <-  c("ID","TIME","AMT","IOBS","EVID","WT")
        sim_data <- sim_data %>%
          mutate(IOBS = if_else(IOBS>0,IOBS,LLOQ/2 )) %>%       # can't be < zero
          mutate(BQL = if_else(IOBS < LLOQ,1,0))  %>%           # set BQL
          mutate(IOBS = if_else(EVID==1,".",as.character(IOBS))) %>% # IOBS missing if dose
          mutate(BQL = if_else(EVID==1,".",as.character(BQL)))  %>% # BQL == "." if dose
          mutate(IOBSBQL = IOBS)
        sim_data <- sim_data %>%
          mutate(IOBSBQL = if_else(IOBSBQL < LLOQ,as.character(LLOQ/2),as.character(IOBSBQL))) %>% #
          mutate(IOBSBQL = if_else(EVID==1,".",as.character(IOBSBQL))) %>%
          select(ID, TIME, AMT, IOBS, EVID, BQL, IOBSBQL, WT )
        write_csv(sim_data,file=file.path(home_dir, "data", paste0("sim_",i,"_",j,"_",l,"_",m,".csv")))
    }
    }
  }
}

