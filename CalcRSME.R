
CalcRMSE <- function(DV,IPRED){
  RMSE <- rmse(DV, IPRED)
  MAE <- mae(DV, IPRED)
  return <- list(
    RMSE = RMSE,
    MAE = MAE
  )
}
