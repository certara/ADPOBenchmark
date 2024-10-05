make_data <- function(home_dir, pyDarwinInterpreter, nmfe_path, nlme_dir, gcc_dir){
  ############## Create data_sim.csv ##############
  LLOQ <- 0.1
  set.seed(1)
  nsubs <-  60
  dose <-  100000 # mcg
  WT_mean <- 70
  WT_CV <- 0.2
  mean_times <- c(0.5,2,6,12,24,48,96,120)
  cv_times <- 0.1
  data <- data.frame(matrix(-999, nrow=nsubs*9, ncol=7))
  colnames(data) <- c("ID","TIME","AMT","DV","EVID","BQL","WT")
  WTs <- exp(rnorm(nsubs,log(WT_mean),WT_CV))
  for(this_id in 1:nsubs){
    times <- exp(rnorm(8,log(mean_times),cv_times))
    this_data <- matrix(-999,nrow=9,ncol=7)
    this_data[,1] <-  this_id
    this_data[,2] <- c(0,round(times,2))
    this_data[,3] <- c(dose,rep(".",8))
    this_data[,4] <- "."
    this_data[,5] <- c(1,rep(2,8))
    this_data[,6] <- 0 ## BQL
    this_data[,7] <- round(WTs[this_id],1) ## BQL
    data[(((this_id-1)*9+1):(this_id*9)),] <- this_data
  }

  data <- data %>%
    mutate(ID = as.integer(ID),TIME = as.numeric(TIME) ) %>%
    arrange(ID,TIME)

  pyDarwin_makedataDir <-
    file.path(home_dir,"data","pyDarwin_makedata")
  if (!dir.exists(pyDarwin_makedataDir)) {
    stop("tokens and template files are expected in ", pyDarwin_makedataDir)
  }

  pyDarwin_makedata_data_sim <-
    file.path(pyDarwin_makedataDir,"data_sim.csv")

  write.csv(data,
            file=pyDarwin_makedata_data_sim,
            quote = FALSE,
            row.names = FALSE)

  ############## Create .dat files ##############
  ## run pyDarwin here, need to run before the NONMEM script
  if(dir.exists(file.path(pyDarwin_makedataDir,"run"))){
    unlink(file.path(pyDarwin_makedataDir,"run"),recursive = TRUE)
  }

  if(file.exists(file.path(pyDarwin_makedataDir, "OUT_0_0_0_0.dat"))){
    files <- list.files(pyDarwin_makedataDir,
                        pattern = "OUT_[0-9].*",
                        full.names = TRUE,
                        recursive = FALSE)

    file.remove(files)
  }


  message("Generating 72 NONMEM control files for benchmark")
  message("NONMEM models are not run at this point, only constructing control files")
  message("example control file in nonmem/run/0/01/nm_0_01.mod")
  pyDarwinOptionsSet <-
    create_pyDarwinOptions(author = "Certara NONMEM",
                           algorithm = "EX",
                           working_dir = file.path(pyDarwin_makedataDir, "working"),
                           output_dir = file.path(pyDarwin_makedataDir, "output"),
                           temp_dir = file.path(pyDarwin_makedataDir, "run"),
                           data_dir = pyDarwin_makedataDir,
                           remove_run_dir = FALSE,
                           model_run_timeout = 3600,
                           engine_adapter = "nonmem",
                           nmfe_path = nmfe_path)
  write_pyDarwinOptions(pyDarwinOptionsSet,
                        file = file.path(pyDarwin_makedataDir, "options.json"))

  # tokens.json and template.txt should be already in pyDarwin_makedataDir!
  run_pyDarwin(
    InterpreterPath = pyDarwinInterpreter,
    Flags = c("-u", "-m"),
    DirectoryPath = pyDarwin_makedataDir,
    TemplatePath = "template.txt",
    TokensPath = "tokens.json",
    OptionsPath = "options.json"
  )

  ############## Preparing sim_.csv ##############

  message("Creating 72 simulation data sets, data sets will be in data/sim_x_x_x_x.csv")
  # sequence of token groups doesn't matter here, we do the same to all
  # but to confirm that the phenotype = data set = model number
  # ;;; Model Identifier =  {COMP[6]},{ETAs[7]},{V~WT[3]},{GAMMA[3]}
  # from control file, print out model #'s corresponding to token groups
  unlink(file.path(pyDarwin_makedataDir,"working"))
  unlink(file.path(pyDarwin_makedataDir,"output"))
  curr_model <- 0
  this_gamma <- this_vwt <- this_eta <- this_comp <- 0
  for (this_gamma in 0:1) {
    for (this_vwt in 0:1) {
      for (this_comp in 0:2) {
        for (this_eta in 0:5) {
          curr_model <- curr_model + 1
          file <- file.path(
            pyDarwin_makedataDir,
            paste0(
              "OUT_",
              this_comp,
              "_",
              this_eta,
              "_",
              this_vwt,
              "_",
              this_gamma,
              ".DAT"
            )
          )
          sim_data <- read_table(
            file,
            skip = 0,
            col_names = FALSE,
            show_col_types = FALSE
          )
          colnames(sim_data) <-  c("ID", "TIME", "AMT", "IOBS", "EVID", "WT")
          sim_data <- sim_data %>%
            mutate(IOBS = if_else(IOBS > 0, IOBS, LLOQ / 2)) %>%       # can't be < zero
            mutate(BQL = if_else(IOBS < LLOQ, 1, 0))  %>%           # set BQL
            mutate(IOBS = if_else(EVID == 1, ".", as.character(IOBS))) %>% # IOBS missing if dose
            mutate(BQL = if_else(EVID == 1, ".", as.character(BQL)))  %>% # BQL == "." if dose
            mutate(IOBSBQL = IOBS)
          sim_data <- sim_data %>%
            mutate(IOBSBQL = if_else(
              IOBSBQL < LLOQ,
              as.character(LLOQ / 2),
              as.character(IOBSBQL)
            )) %>% #
            mutate(IOBSBQL = if_else(EVID == 1, ".", as.character(IOBSBQL))) %>%
            select(ID, TIME, AMT, IOBS, EVID, BQL, IOBSBQL, WT)
          write_csv(sim_data, file = file.path(
            home_dir,
            "data",
            paste0(
              "sim_",
              this_comp,
              "_",
              this_eta,
              "_",
              this_vwt,
              "_",
              this_gamma,
              ".csv"
            )
          ))

        }
      }
    }
  }

  ############## NONMEM simulation with given parameters ##############
  pyDarwinOptionsSet <-
    create_pyDarwinOptions(author = "Certara NONMEM",
                           algorithm = "EX",
                           data_dir = dirname(pyDarwin_makedataDir),
                           working_dir = file.path(home_dir, "NONMEM", "working"),
                           output_dir = file.path(home_dir, "NONMEM", "output"),
                           temp_dir = file.path(home_dir, "NONMEM", "run"),
                           remove_run_dir = FALSE,
                           remove_temp_dir = FALSE,
                           model_run_timeout = 999,
                           engine_adapter = "nonmem",
                           nmfe_path = nmfe_path)
  write_pyDarwinOptions(pyDarwinOptionsSet,
                        file = file.path(home_dir, "Pydarwin_nonmem", "options.json"))

  # tokens.json and template.txt should be already!
  run_pyDarwin(
    InterpreterPath = pyDarwinInterpreter,
    Flags = c("-u", "-m"),
    DirectoryPath = file.path(home_dir, "Pydarwin_nonmem"),
    TemplatePath = "template.txt",
    TokensPath = "tokens.json",
    OptionsPath = "options.json",
  )

  message("Editting NONMEM control files for benchmarking")
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

  ############## NLME simulation with given parameters ##############
  message("Creating 72 NLME meta model files with pyDarwin")
  message("NLME models are not run at this point, only constructing control files")
  message("See example mmdl file in nlme/run/0/01/nlme_0_01.mmdl")
  pyDarwinOptionsSet <-
    create_pyDarwinOptions(author = "Certara NLME",
                           algorithm = "EX",
                           data_dir = dirname(pyDarwin_makedataDir),
                           working_dir = file.path(home_dir, "NLME", "working"),
                           output_dir = file.path(home_dir, "NLME", "output"),
                           temp_dir = file.path(home_dir, "NLME", "run"),
                           remove_run_dir = FALSE,
                           remove_temp_dir = FALSE,
                           model_run_timeout = 999,
                           engine_adapter = "nlme",
                           gcc_dir = gcc_dir,
                           nlme_dir = nlme_dir
                           )
  write_pyDarwinOptions(pyDarwinOptionsSet,
                        file = file.path(home_dir, "pydarwin_nlme", "options.json"))

  # tokens.json and template.txt should be already!
  run_pyDarwin(
    InterpreterPath = pyDarwinInterpreter,
    Flags = c("-u", "-m"),
    DirectoryPath = file.path(home_dir, "pydarwin_nlme"),
    TemplatePath = "template.txt",
    TokensPath = "tokens.json",
    OptionsPath = "options.json",
  )

  message("Editting NLME for benchmarking")
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
