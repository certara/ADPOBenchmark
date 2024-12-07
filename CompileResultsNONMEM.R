timeout <- 3600 * 6
run_NONMEM <- function(home_dir, nmfe_path) {
  working_dir <- file.path(home_dir, "NONMEM")
  #setwd(working_dir)

  # run all models for time, reboot computer first
  # define values for the models/data sets to be used
  CompNTHETA <- c(0, 2, 4)
  VWTNTHETA <- c(0, 1)
  GAMMANTHETA <- c(0, 1)
  ETANOMEGA <- c(1, 2, 3, 4, 3, 6)
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
    crash = as.logical(),
    return_code = as.integer()
  )

  Curr_model <- 0
  Start_model <- 0
  nsamp <- 200
  # set all models/data sets to initial values
  this_gamma <- this_vwt <- this_eta <- this_comp <- 0

  for (this_gamma in 0:1) {
    for (this_vwt in 0:1) {
      for (this_comp in 0:2) {
        for (this_eta in 0:5) {
          Curr_model <- Curr_model + 1
         # Curr_model <- 72
          # skip stuck models
          # if (Curr_model %in% c(9, 12, 27, 45, 48, 63, 67)) next

          if (Curr_model >= Start_model) {
            tryCatch({
              message(
                "############ Starting model ",
                Curr_model ,
                " at ",
                strptime(Sys.time(), format = "%Y-%m-%d %H:%M") ,
                " ############"
              )
              Est_time <- Cov_time <- iterations <-  messages <- -99999
              Rval <- list(RMSE = -99, MAE = -99)
              Covar <- success <- FALSE
              wd <- file.path(working_dir, Curr_model)
              #setwd(wd)
              filenameStem <- paste0("Run", Curr_model)

              xml_file <- file.path(wd, paste0(filenameStem, ".xml"))
              FMSG_file <- file.path(wd, "FMSG")
              if (file.exists(xml_file))    file.remove(xml_file)
              if (file.exists(FMSG_file))    file.remove(FMSG_file)

              # below doesn't work in windows, changed to use rundir in nmfe command
              # command <- paste("cd",
              #                  wd,
              #                  "&&",
              #                  nmfe_path,
              #                  paste0(filenameStem, ".mod "),
              #                  paste0(filenameStem, ".lst"))
              command <- paste(nmfe_path,
                               paste0(filenameStem, ".mod "),
                               paste0(filenameStem, ".lst"),
                               paste0("-rundir=", wd))

              #setwd(home_dir)
              StartTime <- Sys.time() # as.ITime(Sys.time())
              #rval <- system(command, timeout = time)
              if(exists("rval")) rm(rval)
              if(file.exists(xml_file)) file.remove(xml_file)
              rval <- system(command, timeout = timeout) # 3600 * 6 for execution
              if (file.exists(FMSG_file))    file.remove(FMSG_file)
              # can't tell if timeout or crashed??
              if (file.exists(xml_file) & rval !=124) { # if timed out xml file is corrupted, rval is 124
                parms <- GetNMParms(xml_file)
                data_file <- parms$dataFile
                iterations <- parms$iterations
                messages <- parms$messages
                success <- parms$success
                status <- parms$status
                Covar <- parms$covar
                NTHETA <- length(parms$theta)
                if (rval != 124) {
                  # return code for timeout
                  Est_time <-  parms$EstTime
                  Cov_time <- parms$CovTime
                } else{
                  Est_time <-  timeout
                  Cov_time <- -999
                }

                TVVmax <- as.numeric(parms$theta[1])
                TVKM <- as.numeric(parms$theta[2])
                TVV <- as.numeric(parms$theta[3])
                TVKA <- as.numeric(parms$theta[4])
                CVVmax <- sqrt(parms$omega[1])
                data <- read.table(
                  file.path(wd, "Run1Preds.dat"),
                  skip = 1,
                  header = TRUE
                ) %>%
                  filter(TIME > 0) %>%
                  mutate(IPRED = log(IPRED), DV = log(DV))
                # log to get proportional
                RMSE = rmse(data$DV, data$IPRED)
                MAE = mae(data$DV, data$IPRED)
                message(
                  "############ Estimation time for ",
                  wd ,
                  " = ",
                  Est_time,
                  " seconds ############"
                )
              } else{ # timeout from system
                if(rval == 124){
                    Est_time <- Cov_time <-  timeout
                    message_char <- "timed out"
                    converge <- Covar <- success <- finished <- FALSE
                    status <- iterations <- RMSE <- MAE <-   -99999
                    }else{ # crash
                      Est_time <- Cov_time <-  timeout
                      message_char <- "crash"
                      converge <- Covar <- success <- finished <- FALSE
                      status <- iterations <- RMSE <- MAE <-   -99999
                }
                }
            }, error = function(e) {
              Run_time <- 9999999
              converge <- Covar <- success <- finished <- FALSE
              status <- iterations <- RMSE <- MAE <-   -99999
              if (rval == 124) {
                message_char <- "timed out"
              } else{ # probably failed to parsexml??
                message_char <- "crash"
              }
            })

            EndTime <- Sys.time()
            CleanUpNM(wd)
            NOMEGA <- ETANOMEGA[this_eta + 1]

            This_Result = data.frame(
              StartTime = strptime(StartTime, format = "%Y-%m-%d %H:%M:%S"),
              EndTime = strptime(EndTime, format = "%Y-%m-%d %H:%M:%S"),
              Model_num = Curr_model,
              this_comp = this_comp,
              this_eta = this_eta,
              this_vwt = this_vwt,
              this_gamma = this_gamma,
              n_theta = NTHETA,
              n_omega = NOMEGA,
              Est_time = Est_time,
              Cov_time = Cov_time,
              Success = success,
              Status = status,
              Covar = Covar,
              Iterations = iterations,
              Algorithm = "NONMEM",
              Good_inits = TRUE,
              log_path = xml_file,
              control_file = file.path(
                home_dir,
                "",
                "NONMEM",
                Curr_model,
                paste0("Run", Curr_model, ".mod ")
              ),
              data_set = data_file,
              nparms = NTHETA + NOMEGA ,
              messages = paste(messages, collapse = "-"),
              TVVmax = TVVmax,
              TVKM = TVKM,
              TVV = TVV,
              TVKA = TVKA,
              CVVmax = CVVmax,
              RMSE = RMSE,
              MAE = MAE,
              crash = FALSE
            )
            Results <- rbind(Results, This_Result)
            write.csv(
              Results,
              file.path(home_dir, "NONMEMResults.csv"),
              quote = FALSE,
              row.names = FALSE
            )
            try({
              file.copy(
                file.path(home_dir, "NONMEMResults.csv"),
                file.path(home_dir, "NONMEMResults_bak.csv")
               )
            })
          }
        }
      }
    }
  }

  +#setwd(home_dir)
  message("Done at ", strptime(Sys.time(), format = "%Y-%m-%d %H:%M"))
}
