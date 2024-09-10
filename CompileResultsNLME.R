
run_nlme <- function(home_dir, which_version = "standard"){
  nreps <- 200
  if(which_version == "ADPO"){
     Sys.setenv("NLME_HASH" = 1770978959)
     installationDirectory <- "c:/program files/Certara/nnlme_engine_dualmay21_old"
     syngrads <- TRUE
     outputfilename <-  "NLMEResults_ADPO.csv"
  }else{
    installationDirectory <- "c:/program files/Certara/nlme_engine"
    outputfilename <-  "NLMEResults_Standard.csv"
    syngrads <- FALSE
  }
  ETANOMEGA <- c(1,2,3,4,3,6) # no simple way to get n_omega from fit object
  Results <- data.frame(
    StartTime = as.character(),
    EndTime = as.character(),
    ModelNum = as.integer(),
    this_comp = as.integer(),
    this_eta = as.integer(),
    this_vwt = as.integer(),
    this_gamma = as.integer(),
    n_theta = as.integer(),
    n_omega = as.integer(),
    Est_time = as.numeric(),
    Cov_time = as.numeric(),
    Success = as.logical(),
    Covar = as.logical(),
    Iterations = as.integer(),
    Algorithm = as.character(),
    Good_inits = as.logical(),
    log_path = as.character(),
    control_file = as.character(),
    data_set = as.character(),
    nparms = as.integer(),
    messages = as.integer(),
    TVVmax = as.numeric(),
    TVKM = as.numeric(),
    TVV = as.numeric(),
    TVKA = as.numeric(),
    CVVmax = as.numeric(),
    RMSE = as.numeric(),
    MAE = as.numeric(),
    crash = as.logical()

  )
  this_gamma <- this_vwt <- this_eta <- this_comp <- 0
  Curr_model <- 0


  for(this_gamma in 0:1){
    for(this_vwt in 0:1){
      for(this_comp in 0:2){
        for(this_eta in 0:5){
          Curr_model <- Curr_model + 1
          setwd(file.path(home_dir,"nlme",Curr_model))
          model_meta <- create_model_from_metamodel(file.path(home_dir, "nlme",
                                                              Curr_model, paste0("Run",Curr_model,".mmdl")))

          model_orig <- model_meta$model
          message("############ Running ",model_orig@modelInfo@workingDir, " at ",
                  format(Sys.time(), format = "%F %R %Z") , " ############")
          modelNum <- str_sub(rev(setdiff(strsplit(model_orig@modelInfo@workingDir,"/|\\\\")[[1]], ""))[1],
                              start= 0)
          if(file.exists(file.path(model_orig@modelInfo@workingDir,"err2.txt"))){
            file.remove(file.path(model_orig@modelInfo@workingDir,"err2.txt"))
          }
          StartTime <- Sys.time()
          fit <- fitmodel(model_orig,
                          numIterations = 1000,
                          numCores = 1,
                          ODE = "DVERK",
                          sort = FALSE,
                          method = "FOCE-ELS",
                          workingDir = model_orig@modelInfo@workingDir,
                          stdErr = "Auto-Detect",
                          maxStepsODE = 100000,
                          allowSyntheticGradient = syngrads,
                          installationDirectory = installationDirectory,
                          runInBackground = TRUE)

          EndTime <- Sys.time()
          if(file.exists(file.path(model_orig@modelInfo@workingDir,"err2.txt"))){
            EndTime <- Sys.time()
            control_file <- file.path(home_dir,"nlme",Curr_model,paste0("Run",Curr_model,".mmdl"))
            file_conn <- file(model_orig@modelInfo@workingDirm,"err2.txt")
            messages <- readLines(file_conn)
            close(file_conn)
            This_Result <- data.frame(
              StartTime = strptime(StartTime,format = "%Y-%m-%d %H:%M"),
              EndTime = EndTime,
              Model_num = Curr_model,
              this_comp = this_comp,
              this_eta = this_eta,
              this_vwt = this_vwt,
              this_gamma = this_gamma,
              n_theta = -999,
              n_omega = -999,
              Est_time = -999,
              Cov_time = -999,
              Success = FALSE,
              Covar = FALSE,
              Iterations = -999,
              Algorithm = "NLME",
              Good_inits = TRUE,
              log_path = log_path,
              control_file = control_file,
              data_set = data_file,
              nparms = -999,
              messages = messages,
              TVVmax = -99,
              TVKM = -99,
              TVV = -99,
              TVKA = -99,
              CVVmax =-99,
              RMSE = -99,
              MAE = -99,
              crash = TRUE,

            )
            }else{
               if(fit$Overall$RetCode < 4){
                 Success <-  TRUE
               }else{
                 Success <- FALSE
               }
              n_theta <- dim(fit$theta)[1] - 1 # one theta for residual error
              n_omega <- ETANOMEGA[this_eta + 1] # dim(fit$omega)[1]
              logfile <- fit$nlme7engine.log
              Est_time <- logfile[grep(pattern = "engine runtime",x = logfile)]
              Est_time <-  as.numeric(str_trim(str_replace(Est_time," engine runtime \\(secs\\) =","")))
              Cov_time <- logfile[grep(pattern = "stderr runtime",x = logfile)]
              Cov_time <-  as.numeric(str_trim(str_replace(Cov_time," stderr runtime \\(secs\\) =","")))
              Overall <- colnames(read.csv(file.path(model_orig@modelInfo@workingDir,"overall.csv")))
              Covar <- "Condition" %in% Overall
              TVVmax = fit$theta$Estimate[2]
              TVKM = fit$theta$Estimate[3]
              TVV = fit$theta$Estimate[4]
              TVKA = fit$theta$Estimate[5]
              CVVmax = fit$omega$nVmax
              Iterations <- max(fit$ConvergenceData$Iter)
              # data set only in mmdl file
              control_file <- file.path(home_dir,"nlme",Curr_model,paste0("Run",Curr_model,".mmdl"))
              file_conn <- file(control_file)
              control <- readLines(file_conn)
              close(file_conn)
              data_file <- str_replace(control[1],"##DATA ","")
              log_path <- file.path(model_orig@modelInfo@workingDir, "nlme7engine.log")
              DV <- log(fit$residuals$DV)
              IPRED <- log(fit$residuals$IPRED)
              RMSE = rmse(DV, IPRED)
              MAE = mae(DV, IPRED)
              This_Result = data.frame(
                StartTime = as.character(StartTime),
                EndTime = as.character(EndTime),
                Model_num = Curr_model,
                this_comp = this_comp,
                this_eta = this_eta,
                this_vwt = this_vwt,
                this_gamma = this_gamma,
                n_theta = n_theta,
                n_omega = n_omega,
                Est_time = Est_time,
                Cov_time = Cov_time,
                Success = Success,
                Covar = Covar,
                Iterations = Iterations,
                Algorithm = "NLME",
                Good_inits = TRUE,
                log_path = log_path,
                control_file = file.path(home_dir,"nlme",Curr_model,paste0("Run",Curr_model,".mmdl")),
                data_set = data_file,
                nparms = n_theta + n_omega,
                messages = fit$Overall$RetCode,
                TVVmax = TVVmax,
                TVKM = TVKM,
                TVV = TVV,
                TVKA = TVKA,
                CVVmax = CVVmax,
                RMSE = RMSE,
                MAE = MAE,
                crash = FALSE
              )
            }
          Results <- rbind(Results, This_Result)
          write.csv(Results,file.path(home_dir,outputfilename),quote= FALSE, row.names = FALSE)
        }
      }
    }
  }

  message("Done at ", format(Sys.time(), format = "%F %R %Z"))
}
