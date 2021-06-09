roll.win.rmse.nn.wge = function(series, horizon = 1, s = 0, d = 0, fit_model)
{
  
  lags = fit_model$lags
  print(lags)
  trainingSize = max((lags + s + d + 2),9) #Long story but forecast::ets does not like series smaller than 8 unit and nnfor calles "trendcheck" when difforder = 0 which returns error for lenth < 8.
  numwindows = length(series)-(trainingSize + horizon) + 1
  RMSEHolder = numeric(numwindows)
  
  print(paste("Please Hold For a Moment, TSWGE is processing the Rolling Window RMSE with", numwindows, "windows."))
  
  for( i in 1:numwindows)
  {
    print("0")
    print((i+(trainingSize-1)))
    fit = mlp(ts(series[i:(i+(trainingSize-1))]), model = fit_model)
    print("1")  
    forecasts = predict(fit,h = horizon)$mean
    print("2")  
    print(forecasts)
    print(series[(trainingSize+i):(trainingSize+ i + (horizon) - 1)])
    #RMSE = sqrt(mean((series[(trainingSize+i):(trainingSize+ i + (horizon) - 1)] - forecasts)^2))
    RMSE = mean((series[(trainingSize+i):(trainingSize+ i + (horizon) - 1)] - forecasts)^2)
    print("3")
    RMSEHolder[i] = RMSE
    print("4")
  }
  
  RMSEHolder
  hist(RMSEHolder, main = "RMSEs for Individual Windows")
  WindowedRMSE = mean(RMSEHolder)
  
  print("The Summary Statistics for the Rolling Window RMSE Are:")
  print(summary(RMSEHolder))
  print(paste("The Rolling Window RMSE is: ",round(WindowedRMSE,3)))
  return(list(rwRMSE = WindowedRMSE, numwindows = numwindows, horizon = horizon, s = s, d = d))
}
