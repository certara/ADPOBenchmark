
library(tidyvpc)
library(dplyr)
library(readtext)
library(stringr)  
GetNPDE <- function(curFolder, this_model, nsamp){
 tryCatch({
    
  CmaxRows <- rep(c(FALSE,TRUE,TRUE, FALSE,FALSE,FALSE,FALSE),60*10)
  CminRows <- rep(c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE),60*10) 
  model_name <- sprintf("%02d", this_model)  
  setwd(curFolder) 
  fileStem <-  paste0("NM_0_",model_name)  
  if (file.exists("OBS_DAT")){file.remove("OBS_DAT")} 
  if (file.exists("SIM_DAT")){file.remove("SIM_DAT")} 
  count <- 0
  if(file.exists("RunNPDE.mod")){file.remove("RunNPDE.mod")}
  while(!file.exists("RunNPDE.mod") & count<10){
    file.copy(paste0(fileStem,".mod"),"RunNPDE.mod")
    Sys.sleep(0.5)
    count <-  count + 1
  }
 
  count <- 0
  if(file.exists("RunNPDE.xml")){file.remove("RunNPDE.xml")}
  while(!file.exists("RunNPDE.xml") & count<10){
      OK <- file.copy(paste0(fileStem,".xml"),"RunNPDE.xml")
      Sys.sleep(0.5)
      count <-  count + 1
  }
 
  count <- 0
  if(file.exists("RunNPDE.lst")){file.remove("RunNPDE.lst")}
  while(!file.exists("RunNPDE.lst") & count<10){
    OK <- file.copy(paste0(fileStem,".lst"),"RunNPDE.lst")
    Sys.sleep(0.5)
    count <-  count + 1
  } 
  if(file.exists("RunNPDE.mod")){
    tryCatch({
      rval <- shell(paste("update_inits RunNPDE.mod"), intern = TRUE)
      if(file.exists("RunNPDE.mod")){
        file <- readtext("RunNPDE.mod", verbosity = 0)$text
        file <- str_replace(file,"\\$EST",";;$EST") 
        file <- str_replace(file,"\\$COV",";;$COV") 
        file <- str_replace(file,"NOSIGMABOUNDTEST",";; NOSIGMABOUNDTEST")  
        
        file <- str_replace(file,"\\$ERROR",paste("$ERROR \n  REP = IREP \n  IPRED = F")) 
        file <- c(file,
                  paste("$SIM ONLYSIM (1233) SUBPROB =",nsamp,"\n",
                        "$TABLE ID TIME DV MDV NOPRINT NOHEADER NOAPPEND FILE=OBS_DAT\n",
                        "$TABLE ID REP TIME DV IPRED PRED MDV NOPRINT NOHEADER NOAPPEND FILE=SIM_DAT"))
        file <- unlist(strsplit(file,"\n")) 
        fileConn <- file("RunNPDE.mod")
        writeLines(file, fileConn)
        close(fileConn)
        rval <- shell("nmfe75 RunNPDE.mod RunNPDE.lst",intern = TRUE) 
      }
    }
      )}
  if(file.exists(file.path(curFolder, "OBS_DAT")) & file.exists(file.path(curFolder, "SIM_DAT"))){
    obs_data <- read.table(file.path(curFolder, "OBS_DAT"))
    obs_data <- obs_data[1:480,]
    colnames(obs_data) <- c("ID","TIME","DV","MDV")
    obs_data <- obs_data %>% filter(MDV==0)
    sim_data <- read.table(file.path(curFolder, "SIM_DAT")) 
    colnames(sim_data) <- c("ID","REP","TIME","DV","IPRED","PRED","MDV")
    sim_data <- sim_data %>% filter(MDV==0)
    obs_data$PRED <- sim_data %>% filter(REP==1) %>% select(PRED) 
    npde <- observed(obs_data, x=NULL, y=DV) %>%
      simulated(sim_data, y=DV) %>%
      npde(id=ID)
    CmaxNPDEs <- npde$npdesim$npde[CmaxRows]
    CmaxMean <- mean(CmaxNPDEs)
    CmaxSD <- sd(CmaxNPDEs)
    CminNPDEs <- npde$npdesim$npde[CminRows]
    CminMean <- mean(CminNPDEs)
    CminSD <- sd(CminNPDEs)  
    xml_file <- file.path(curFolder,"RunNPDE.xml") 
    data <- read_xml(xml_file, encoding = "ASCII") 
    
    runtime_node <- xml_find_all(data, "//nm:total_cputime") 
    NPDETime <- as.numeric(xml_text(runtime_node))
    message("NPDE Cmax =", paste(round(CmaxMean,4),
                                  "(",round(CmaxSD,4),
                                  "), Cmin = ", round(CminMean,4),
                                  "(",round(CminSD,4),")"))
    
  }else{
    message("#### Failed NPDE, no obs_dat and or sim_dat for ", this_model, " ######") 
      model <-  this_model
      CmaxMean <-  -999
      CmaxSD <-  -999
      CminMean <-  -999
      CminSD <-  -999
      NPDETime <- -999
  }
  
  return(list(
    model = this_model,
    CmaxMean = CmaxMean,
    CmaxSD = CmaxSD,
    CminMean = CminMean,
    CminSD = CminSD,
    NPDETime = NPDETime
  ))
 },
 error = function(e){ 
   message("#### Failed NPDE ", this_model, ", message ", e, " ######") 
   return(list(
     model = this_model,
     CmaxMean = -999,
     CmaxSD = -999,
     CminMean = -999,
     CminSD = -999,
     NPDETime = -999
   ))
 }
)
}
