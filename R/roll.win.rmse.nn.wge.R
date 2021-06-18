roll.win.rmse.nn.wge = function(series, horizon = 1, s = 0, d = 0, fit_model)
{
  
  lags = fit_model$lags
  trainingSize = max((lags + s + d + 2),9) #Long story but forecast::ets does not like series smaller than 8 unit and nnfor calles "trendcheck" when difforder = 0 which returns error for lenth < 8.
  numwindows = length(series)-(trainingSize + horizon) + 1
  RMSEHolder = numeric(numwindows)
  MSEHolder = numeric(numwindows)
  
  print(paste("Please Hold For a Moment, TSWGE is processing the Rolling Window RMSE with", numwindows, "windows."))
  
  for( i in 1:numwindows)
  {
    fit = mlp(ts(series[i:(i+(trainingSize-1))]), model = fit_model)
    forecasts = predict(fit,h = horizon)$mean
    MSE = mean((series[(trainingSize+i):(trainingSize+ i + (horizon) - 1)] - forecasts)^2)
    MSEHolder[i] = MSE
  }
  
  RMSEHolder = sqrt(MSEHolder)
  hist(RMSEHolder, main = "RMSEs for Individual Windows")
  WindowedMSE = mean(MSEHolder)
  WindowedRMSE = mean(RMSEHolder)
  
  print("The Summary Statistics for the Rolling Window RMSE Are:")
  print(summary(RMSEHolder))
  print(paste("The Rolling Window RMSE is: ",round(WindowedRMSE,3)))
  invisible(list(rwRMSE = WindowedRMSE, rwMSE = WindowedMSE, numwindows = numwindows, horizon = horizon, s = s, d = d))
}
