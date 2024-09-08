library("readtext")
library(stringr)
make_nm_control <-  function(directory, model){
      options(warn=-1)
      file <- file.path(directory,paste0("NM_0_",model,".mod"))
      control <- readtext(paste0(file))$text
      control_sav <- control
      control <- str_replace(control,"\\$DATA      ","$DATA test")
      control <- str_replace(control,"c:/dual_numbers/pydarwin/data.csv","")
      control <- str_replace(control,"\\^$SIM ","$EST METHOD=COND INTER NOHABORT MAX=9999\n#$COV UNCOND PRECOND = 2")
      control <- str_split_1(control,"\n")
      genotype <- str_replace(control[grep("^;; Genotype:",control)],";; Genotype: \\[","")
      genotype <- str_replace(genotype,"\\]","")
      genotype <- str_replace_all(genotype,", ","_")
      data_file_name <- file.path("c:/dual_numbers/data",paste0("sim_",genotype,".csv"))
      new_control_name <- paste0("NM_",genotype,".mod") # extra _0 for no bql
      # replace $data and $SIM

      control <- str_replace(control,"\\$DATA      ",paste0("$DATA ",data_file_name))
      control <- str_replace(control,"c:/dual_numbers/pydarwin/data.csv","")
      control <- str_replace(control,"IPRED =F",";;")

      control <- str_replace(control,"\\$INPUT       ID TIME AMT DV EVID DROP WT","$INPUT ID TIME AMT DV DROP DROP DROP WT")
      control <- str_replace(control,"\\$SIM ","$EST METHOD=COND NOHABORT INTER MAX=9999\n$COV UNCOND PRECOND = 2  ")
      control <- str_replace(control,"ONLYSIM ",";;")
      control <- str_replace(control,"\\$TABLE ",";;")
      control <- str_replace(control," NOPRINT ",";; NOPRINT")
      control <- str_replace(control,"DUM=\\(LLOQ-IPRED\\)/\\(W\\)",";;")
      control <- str_replace(control,"CUMD=PHI\\(DUM\\) ",";;")
      control <- str_replace(control,"IOBS = IPRED+",";;")
      control <- str_replace(control,"Y=IOBS","Y = IPRED+EPS(1)*W")
      control <- str_replace(control,"BQL=0",";;")
      control_NObql <- str_split_1(control,"\n")

      fileConn <- file(file.path(home_dir,"NONMEM",new_control_name))
      writeLines(control, fileConn)
      close(fileConn)
  }

