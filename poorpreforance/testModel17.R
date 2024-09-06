library(stringr)
library(dplyr)
library(Certara.RsNLME)
library(data.table)  
library(tidyvpc)  
model_meta <- create_model_from_metamodel("c:/git/dual_numbers/compileresults/nlme/17/Run17.mmdl") 
model_orig <- model_meta$model 
fit <- fitmodel(model_orig,  
                numIterations = 1000, 
                numCores = 1,
                ODE = "DVERK",  
                sort = FALSE,
                method = "FOCE-ELS", 
                workingDir = model_orig@modelInfo@workingDir,
                stdErr = "Auto-Detect",
                maxStepsODE = 100000,
                allowSyntheticGradient = FALSE,  
                runInBackground = TRUE, 
                installationDirectory = "C:/Program Files/Certara/NLME_Engine")

logfile <- fit$nlme7engine.log
Est_time <- logfile[grep(pattern = "engine runtime",x = logfile)]  
StdEst_time <-  as.numeric(str_trim(str_replace(Est_time," engine runtime \\(secs\\) =","")))
Cov_time <- logfile[grep(pattern = "stderr runtime",x = logfile)]  
StdCov_time <-  as.numeric(str_trim(str_replace(Cov_time," stderr runtime \\(secs\\) =","")))
# ADPO

Sys.setenv("NLME_HASH" = 1788678031) # don't need for standard NLME
Sys.time()
fit <- fitmodel(model_orig,  
                numIterations = 1000, 
                numCores = 1,
                ODE = "DVERK",  
                sort = FALSE,
                method = "FOCE-ELS", 
                workingDir = model_orig@modelInfo@workingDir,
                stdErr = "Auto-Detect",
                maxStepsODE = 100000,
                allowSyntheticGradient = TRUE,  
                runInBackground = TRUE, 
                installationDirectory = "c:/program files/Certara/nlme_engine_3Jun")

logfile <- fit$nlme7engine.log
Est_time <- logfile[grep(pattern = "engine runtime",x = logfile)]  
ADPOEst_time <-  as.numeric(str_trim(str_replace(Est_time," engine runtime \\(secs\\) =","")))
Cov_time <- logfile[grep(pattern = "stderr runtime",x = logfile)]  
ADPOCov_time <-  as.numeric(str_trim(str_replace(Cov_time," stderr runtime \\(secs\\) =","")))
message("ADPO EST time/FD EST time  = ", ADPOEst_time/StdEst_time)
# check if dual