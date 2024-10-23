
run_nlme <- function(home_dir, which_version = "standard", nlme_dirs) {
  Sys.setenv("INSTALLDIR" = nlme_dirs[which_version])
  timeout <-  3600*6
  message("Starting NLME benchmarking")
  # BAK files are just so we can look at intermediate results with locking the file
  if (which_version == "ADPO") {
    # ADPO and NoHessian are internal only, not available in public repo
    hash_file <- "d:/users/hash.txt"
    txt <- as.numeric(readtext(hash_file)$text)
    Sys.setenv("NLME_HASH" = as.numeric(txt))
    syngrads <- TRUE
    outputfilename <-  "NLMEResults_ADPO.csv"
    backupfilename <-  "NLMEResults_ADPOBAK.csv"

  } else if (which_version == "NoHessian") {
    hash_file <- "d:/users/hash.txt"
    txt <- as.numeric(readtext(hash_file)$text)
    Sys.setenv("NLME_HASH" = txt)
    syngrads <- FALSE
    outputfilename <-  "NLMEResults_NoHessian.csv"
    backupfilename <-  "NLMEResults_NoHessianBAK.csv"
  } else if (which_version == "standard") {
    syngrads <- FALSE
    outputfilename <-  "NLMEResults_Standard.csv"
    backupfilename <-  "NLMEResults_StandardBAK.csv"
  }  else if (which_version == "noCache") {
    hash_file <- "d:/users/hash.txt"
    txt <- as.numeric(readtext(hash_file)$text)
    Sys.setenv("NLME_HASH" = as.numeric(txt))
    syngrads <- FALSE
    outputfilename <-  "NLMEResults_NoCache.csv"
    backupfilename <-  "NLMEResults_NoCacheBAK.csv"
  }

  ETANOMEGA <- c(1, 2, 3, 4, 3, 6) # no simple way to get n_omega from fit object
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
    Status = as.integer(),
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
  Curr_model <- 0
  for (this_gamma in 0:1) {
    for (this_vwt in 0:1) {
      for (this_comp in 0:2) {
        for (this_eta in 0:5) {
          Curr_model <- Curr_model + 1
          StartTime <- Sys.time()
          runOk <- TRUE
          tryCatch({
            model_meta <- create_model_from_metamodel(file.path(
              home_dir,
              "NLME",
               Curr_model,
              paste0("Run", Curr_model, ".mmdl")
            ),
            directoryToRun = file.path(home_dir, "NLME", Curr_model, paste0("Run", Curr_model)))

            model_orig <- model_meta$model
            message(
              "############ Running ",
              model_orig@modelInfo@workingDir,
              " at ",
              format(Sys.time(), format = "%F %R %Z") ,
              " with ",
              which_version,
              " ############"
            )
            modelNum <- str_sub(rev(setdiff(
              strsplit(model_orig@modelInfo@workingDir, "/|\\\\")[[1]],
              ""
            ))[1], start = 0)
            if (file.exists(file.path(model_orig@modelInfo@workingDir, "err2.txt"))) {
              file.remove(file.path(model_orig@modelInfo@workingDir, "err2.txt"))
            }
            if(exists("fit")) rm(fit)
            if(file.exists(model_orig@modelInfo@workingDir)){
             unlink(model_orig@modelInfo@workingDir, recursive = TRUE)
          }
              Rval <- withTimeout(fit <- fitmodel(
                model_orig,
                numIterations = 9999,
                numCores = 1,
                ODE = "DVERK",
                sort = FALSE,
                method = "FOCE-ELS",
                workingDir = model_orig@modelInfo@workingDir,
                stdErr = "Auto-Detect",
                maxStepsODE = 100000,
                allowSyntheticGradient = syngrads,
                #  installationDirectory = installationDirectory,
                runInBackground = FALSE
              ),
            substitute=TRUE,
            envir=parent.frame(),
            timeout = timeout + 60,
            cpu=99999,
            elapsed = timeout + 60, # 60 seconds for compiling
            onTimeout="warning")
          }, error = function(cond) {
            runOk <- FALSE # this fails if no license
            message("Error!!\nNo License???")
          })
          EndTime <- Sys.time()

          if (!(runOk) |
              file.exists(file.path(model_orig@modelInfo@workingDir, "err2.txt")) |
            !exists("fit")) {
            EndTime <- Sys.time()
            control_file <- file.path(home_dir,
                                      "nlme",
                                      Curr_model,
                                      paste0("Run", Curr_model, ".mmdl"))
            if (file.exists(file.path(model_orig@modelInfo@workingDir, "err2.txt"))) {
              file_conn <- file(model_orig@modelInfo@workingDirm, "err2.txt")
              messages <- readLines(file_conn)
              close(file_conn)
            } else{
              if(exists("fit")){
              messages <- "Run Failed, unknown error"
              }else{
                messages <- "timed out"
              }
            }
            This_Result <- data.frame(

              StartTime = strptime(StartTime, format = "%Y-%m-%d %H:%M:%S"),
              EndTime = strptime(EndTime, format = "%Y-%m-%d %H:%M:%S"),
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
              Status = -999,
              Covar = FALSE,
              Iterations = -999,
              Algorithm = which_version, #"NLME",
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
              CVVmax = -99,
              RMSE = -99,
              MAE = -99,
              crash = TRUE
            )
          } else{ # all OK
            if (fit$Overall$RetCode < 4 & fit$Overall$RetCode > 0) {
              Success <-  TRUE
            } else{
              Success <- FALSE
            }
            status <- fit$Overall$RetCode
            n_theta <- dim(fit$theta)[1] - 1 # one theta for residual error
            n_omega <- ETANOMEGA[this_eta + 1] # dim(fit$omega)[1]
            logfile <- fit$nlme7engine.log
            Est_time <- logfile[grep(pattern = "engine runtime", x = logfile)]
            Est_time <-  as.numeric(str_trim(
              str_replace(Est_time, " engine runtime \\(secs\\) =", "")
            ))
            Cov_time <- logfile[grep(pattern = "stderr runtime", x = logfile)]
            Cov_time <-  as.numeric(str_trim(
              str_replace(Cov_time, " stderr runtime \\(secs\\) =", "")
            ))

            Covar <- !is.na(fit$Overall$Condition)
            TVVmax = fit$theta$Estimate[fit$theta$Parameter == 'tvVmax']
            TVKM = fit$theta$Estimate[fit$theta$Parameter == 'tvKm']
            TVV = fit$theta$Estimate[fit$theta$Parameter == 'tvV']
            TVKA = fit$theta$Estimate[fit$theta$Parameter == 'tvKa']
            CVVmax = fit$omega$nVmax[1]
            Iterations <-  -999
            try({
              if (!is.null(fit$ConvergenceData) &&
                  all(!is.na(fit$ConvergenceData$Iter))) {
                Iterations <- max(fit$ConvergenceData$Iter)
              } else{
                Iterations <-  -999
              }
            })

            control_file <- file.path(home_dir,
                                      "NLME",
                                      Curr_model,
                                      paste0("Run", Curr_model, ".mmdl"))
            file_conn <- file(control_file)
            control <- readLines(file_conn)
            close(file_conn)
            data_file <- str_replace(control[1], "##DATA ", "")
            log_path <- file.path(model_orig@modelInfo@workingDir,
                                  "nlme7engine.log")
            if (!is.null(fit$residuals)) {
              DV <- log(fit$residuals$DV)
              IPRED <- log(fit$residuals$IPRED)
              RMSE <-  rmse(DV, IPRED)
              MAE <-  mae(DV, IPRED)
            } else{
              DV <-  IPRED <-   RMSE <- MAE  <- -999
            }
            This_Result = data.frame(
              StartTime = strptime(StartTime, format = "%Y-%m-%d %H:%M:%S"),
              EndTime = strptime(EndTime, format = "%Y-%m-%d %H:%M:%S"),
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
              Status = status,
              Covar = Covar,
              Iterations = Iterations,
              Algorithm = which_version, #"NLME",
              Good_inits = TRUE,
              log_path = log_path,
              control_file = file.path(
                home_dir,
                "nlme",
                Curr_model,
                paste0("Run", Curr_model, ".mmdl")
              ),
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
          CleanUpNLME(model_orig@modelInfo@workingDir)
          Results <- rbind(Results, This_Result)
          write.csv(
            Results,
            file.path(home_dir, outputfilename),
            quote = FALSE,
            row.names = FALSE
          )
          try({
            write.csv(
              Results,
              file.path(home_dir, backupfilename),
              quote = FALSE,
              row.names = FALSE
            )
          })
        }
      }
    }
  }


  +message("Done at ", format(Sys.time(), format = "%F %R %Z"))
 }
