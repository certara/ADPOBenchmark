GetTrueParms <- function(home_dir){
  trueValues <- data.frame(
    Model = as.integer(),
    TVVmax = as.numeric(),
    TVKM = as.numeric(),
    TVV = as.numeric(),
    TVKA = as.numeric()
  )
  this_gamma <- this_vwt <- this_eta <- this_comp <- 0
  Curr_model <- 0
  for(this_gamma in 0:1){
    for(this_vwt in 0:1){
      for(this_comp in 0:2){
        for(this_eta in 0:5){
          Curr_model <- Curr_model + 1
          parms <- read.table(file.path(home_dir,"DATA", paste0("PARMS_", this_comp,"_",this_eta,"_",this_vwt,"_" ,this_gamma,".DAT")),
                            skip=1,header=TRUE) %>%
            dplyr::summarise(VMAX = exp(mean(log(VMAX))),
                             KM = exp(mean(log(KM))),
                             V2 = exp(mean(log(V2))),
                             KA = exp(mean(log(KA))))
          parms$Model <- Curr_model
          trueValues <- rbind(trueValues,parms)
        }
      }
    }
  }
  return(trueValues)
}
