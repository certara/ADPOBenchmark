check_data <- function(home_dir, curr_model){
  # check if model an data file agree with model #
  control <- readtext::readtext(file.path(home_dir,"NONMEM",curr_model,paste0("Run",curr_model,".mod")),verbosity =0)$text
  control <- str_split(control,"\n")
  phenotype <- grep(";;; Model Identifier", unlist(control))
  phenotype <- control[[1]][phenotype[[1]]]
  phenotype <- str_trim(str_replace(phenotype,"  ;;; Model Identifier = ",""))
  phenotype <- str_split(phenotype,",")[[1]]
  datafile <-  control[[1]][grep(tolower("^\\$DATA"), tolower(unlist(control)))]
  datafile <- tolower(datafile)
  datafile <- str_replace(datafile,"..\\\\\\\\..\\\\\\\\data\\\\\\\\sim_" ,"")
  datafile <- str_trim(str_replace(datafile,"^\\$data" ,""))
  datafile <- str_trim(str_replace(datafile,".csv ignore=@",""))
  datafile <- str_split(datafile,"_")[[1]]
  #curr_data_type <- c(as.character(this_comp),as.character(this_eta),as.character(this_vwt),as.character(this_gamma))
  if(!(all(datafile==phenotype))){
    message("NONMEN Phenotype not matched in ", datafile, " current model = ", curr_model)
  }

  # same for NLME
  control <- readtext::readtext(file.path(home_dir,"NLME",curr_model,paste0("Run",curr_model,".mmdl")),verbosity =0)$text
  control <- str_split(control,"\n")
  phenotype <- grep("## Genotype: ", unlist(control))
  phenotype <- control[[1]][phenotype[[1]]]
  phenotype <- str_trim(str_replace(phenotype,"## Genotype: \\[",""))
  phenotype <- str_trim(str_replace(phenotype,"\\]",""))
  phenotype <- str_split(phenotype,", ")[[1]]
  datafile <-  control[[1]][grep("^##DATA ", unlist(control))]
  datafile <- str_trim(str_replace(datafile,"##DATA",""))
  datafile <- tolower(datafile)
  datafile <- str_trim(str_replace(datafile,tolower("..\\\\\\\\..\\\\\\\\data\\\\\\\\sim_"),""))
  datafile <- str_trim(str_replace(datafile,".csv",""))
  datafile <- str_split(datafile,"_")[[1]]
  if(!(all(datafile==phenotype))){
    message("NLME Phenotype not matched in ", datafile, " current model = ", curr_model)
  }

  # no way to check true parameters data, written by the pydarwin generated file above
  # but, preds, out and parms all have the same nameing e.g,
  # Preds_{COMP[6]}_{ETAs[7]}_{V~WT[3]}_{GAMMA[3]}.DAT
}

